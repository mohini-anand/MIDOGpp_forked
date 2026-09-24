"""
    How the Otsu gate behaves when its window is a human-drawn box instead of the fixed 51 px (DECISIONS.md D12).

    A tight drawn box is simulated as a square of side (51 px component size + 2 * margin), centred on the component the
    51 px gray gate accepted, for every agreement-pool mitotic annotation on the 14 original ROIs. Also checks that
    otsu_window=50 and 51 give identical results (read_padded_patch forces odd windows). Simulated, square and centred,
    gray channel only; precision is not measured. Writes results/box_window_probe.csv, one row per (annotation, margin).
    Usage: box_window_probe.py
"""

import sys
from collections import Counter
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from skimage.measure import label, regionprops

REPO = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO))
from midog_utils import channels as ch  # noqa: E402
from midog_utils import dataset as ds  # noqa: E402
from midog_utils import seed_selection as ss  # noqa: E402
from midog_utils import template_match as tm  # noqa: E402

MARGINS = (0, 2, 5, 10, 15)
OUT_PATH = REPO / "results" / "box_window_probe.csv"


def reason(patch, min_area=50, max_area_frac=0.85, min_solidity=0.5):
    """
        Why tighten_box_otsu accepts or refuses a patch: a mirror of its steps that names the failing one.

        patch (np.ndarray): the gate window.
        min_area (int): tighten_box_otsu's min_area.
        max_area_frac (float): tighten_box_otsu's max_area_frac.
        min_solidity (float): tighten_box_otsu's min_solidity.

        Returns str: "accepted" or the refusal reason.
    """
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    _, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    labels = label(binary, connectivity=2)
    cy, cx = patch.shape[0] // 2, patch.shape[1] // 2
    lab = labels[cy, cx]
    if lab == 0:
        return "centre_not_fg"
    r = next(p for p in regionprops(labels) if p.label == lab)
    y0, x0, y1, x1 = r.bbox
    if not (y0 <= cy < y1 and x0 <= cx < x1):
        return "centre_outside_bbox"
    if r.area < min_area:
        return "too_small"
    if r.area > max_area_frac * patch.size:
        return "too_large"
    if r.solidity < min_solidity:
        return "low_solidity"
    return "accepted"


def main():
    """
        Run the probe over the 14 original ROIs, print the acceptance table and write the CSV.

        Returns None.
    """
    images, annotations = ds.load_annotations(str(REPO / "databases" / "MIDOG++.json"))
    rows, same_50_51, n_checked = [], 0, 0
    for path in sorted((REPO / "images" / "extra_valid").glob("*.tiff")):
        rgb = ds.load_roi(path)
        gray = ch.to_channel(rgb, "gray_inverted")
        gt = ds.image_annotations(annotations, path.name)
        pool, _ = ss.agreement_pool(gt[gt["category_id"] == ds.MITOTIC])
        pool = ss.border_filter(pool, tm.PATCH_SIZE // 2, rgb.shape)
        for _, a in pool.iterrows():
            cx, cy = float(a["cx"]), float(a["cy"])
            g51 = ss.tightened_template_box(gray, cx, cy, otsu_window=51)
            g50 = ss.tightened_template_box(gray, cx, cy, otsu_window=50)
            n_checked += 1
            same_50_51 += int(g51 == g50)
            if g51 is None:
                continue
            b, tx, ty = g51
            for m in MARGINS:
                side = min(b + 2 * m, tm.PATCH_SIZE)
                patch = tm.read_padded_patch(gray, tx, ty, side)
                got = ss.tightened_template_box(gray, tx, ty, otsu_window=side)
                rows.append(dict(roi=path.name, ann_id=int(a["ann_id"]), b51=b, margin=m, side=side, patch_side=patch.shape[0], accepted=got is not None, why=reason(patch), new_b=None if got is None else got[0]))
        print(path.name, "done", flush=True)

    df = pd.DataFrame(rows)
    OUT_PATH.parent.mkdir(exist_ok=True)
    df.to_csv(OUT_PATH, index=False)
    print(f"\notsu_window 50 vs 51 identical on {same_50_51}/{n_checked} annotations")
    print(f"annotations accepted at 51: {df[df.margin == 0].shape[0]} of {n_checked}; b51 median {df[df.margin == 0].b51.median()}, range {df.b51.min()}-{df.b51.max()}")
    print("\nacceptance by margin around the component (box side = b51 + 2*margin):")
    for m, g in df.groupby("margin"):
        c = Counter(g.why)
        print(f"  margin {m:>2}px: accepted {g.accepted.mean():6.1%}  ({g.accepted.sum()}/{len(g)})  reasons: {dict(c.most_common())}")
        acc = g[g.accepted]
        if len(acc):
            print(f"               new base_size vs b51: same {np.mean(acc.new_b == acc.b51):.0%}, smaller {np.mean(acc.new_b < acc.b51):.0%}, larger {np.mean(acc.new_b > acc.b51):.0%}")
    print(f"-> {OUT_PATH}")


if __name__ == "__main__":
    main()
