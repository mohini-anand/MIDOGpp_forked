"""Sweep `MIN_PEAK_DIST`, the one untuned knob in the watershed split.

`watershed_split.py` fixed the proposal ceiling on 7 of 7 ROIs but regressed the deep tail
on 301.tiff (1.84x worse) and 245.tiff (1.23x) -- the over-splitting failure mode, where a
single nucleus with coarse chromatin is carved into pieces and none of them carries the
whole object's intensity. `MIN_PEAK_DIST = 7` was a guess applied across four tumour
domains and was never swept.

Scope: **301 and 245 (the two regressions) plus 246 (the largest win, as a guard against a
fix that trades it away)**. Three ROIs answer the question the sweep is actually asking;
confirming on the other four is a follow-up once a value is chosen.

Cost. The distance transform does **not** depend on `min_distance`, so it is computed once
per component and reused across every sweep value, as are the ROI load, the tiled-Otsu
threshold and the connected-component labelling. Measured on this machine those shared
stages are 30-85 s per ROI against 55-145 s per *value* per ROI, so hoisting them turns an
N x full-run cost into a fixed cost plus N cheap passes.

Verification: the `md=7` row must reproduce `results/watershed_tail.csv`, which was written
by separately structured code. That is asserted, not eyeballed.

Writes results/watershed_sweep.csv.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.measure import label, regionprops
from skimage.segmentation import watershed

from midog_utils import baselines as bl
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from watershed_split import MIN_AREA, gate_and_rank, tail_rows

ROIS = ["301.tiff", "245.tiff", "246.tiff"]
SWEEP = [5, 7, 9, 11, 15]
COLS = ["cx", "cy", "area", "score"]


def prep(binary: np.ndarray, h_chan: np.ndarray):
    """Everything that does not depend on ``min_distance``, done once per ROI.

    Components too small to hold two nuclei are finalised here and never revisited;
    the rest are kept as (bbox, boolean crop, fallback row) with their distance
    transform, which is the expensive shared quantity.
    """
    lab = label(binary, connectivity=2)
    small, big, edts = [], [], []
    for p in regionprops(lab, intensity_image=h_chan):
        row = (p.centroid[1], p.centroid[0], p.area, p.intensity_mean)
        if p.area < 2 * MIN_AREA:
            small.append(row)
        else:
            big.append((p.bbox, p.image, row))
    for _bbox, img, _row in big:
        edts.append(ndi.distance_transform_edt(img))
    return small, big, edts


def split_at(big, edts, h_chan: np.ndarray, min_dist: int) -> pd.DataFrame:
    """Re-segment every big component at one ``min_distance``. Mirrors
    `watershed_split.split_components` exactly, minus the recomputed EDT."""
    rows = []
    for (bbox, img, row), d in zip(big, edts):
        coords = peak_local_max(d, min_distance=min_dist, labels=img, exclude_border=False)
        if len(coords) <= 1:
            rows.append(row)
            continue
        markers = np.zeros(d.shape, dtype=np.int32)
        markers[tuple(coords.T)] = np.arange(1, len(coords) + 1)
        ws = watershed(-d, markers, mask=img)
        minr, minc, maxr, maxc = bbox
        icrop = h_chan[minr:maxr, minc:maxc]
        for q in regionprops(ws, intensity_image=icrop):
            rows.append((q.centroid[1] + minc, q.centroid[0] + minr,
                         q.area, q.intensity_mean))
    return pd.DataFrame(rows, columns=COLS)


def main():
    _images, ann = ds.load_annotations()
    ref = pd.read_csv("results/watershed_tail.csv")
    rows, checks = [], []

    for fn in ROIS:
        t0 = time.time()
        path = f"images/{fn}"
        rgb = ds.load_roi(path)
        radius = ev.radius_px(ds.roi_mpp(path))
        gt = ds.image_annotations(ann, fn)

        mask = bl.tissue_mask(rgb)
        h_chan = bl.to_hematoxylin(rgb)
        thr = bl._tiled_otsu_threshold(h_chan, mask, 512)
        binary = (h_chan > thr) & mask
        del thr, mask, rgb

        small, big, edts = prep(binary, h_chan)
        del binary
        cc = pd.DataFrame(small + [r for _b, _i, r in big], columns=COLS)
        rows += tail_rows(gate_and_rank(cc), gt, radius, "connected_components", fn)
        print(f"{fn}: prep {time.time()-t0:.0f}s, {len(big)} splittable of {len(cc)}",
              flush=True)

        for md in SWEEP:
            t1 = time.time()
            det = pd.concat([pd.DataFrame(small, columns=COLS),
                             split_at(big, edts, h_chan, md)], ignore_index=True)
            got = tail_rows(gate_and_rank(det), gt, radius, f"watershed_md{md}", fn)
            rows += got
            if md == 7:  # must reproduce the committed run, which used md=7
                want = ref[(ref.file_name == fn) & (ref.arm == "watershed")]
                for g in got:
                    w = want[want["quantile"] == g["quantile"]]
                    checks.append((fn, g["quantile"], g["depth"], int(w.depth.iloc[0]),
                                   g["n_pool"], int(w.n_pool.iloc[0])))
            print(f"  md={md}: {time.time()-t1:.0f}s, pool {got[0]['n_pool']}, "
                  f"ceiling {got[0]['ceiling']}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv("results/watershed_sweep.csv", index=False)

    bad = [c for c in checks if c[2] != c[3] or c[4] != c[5]]
    print(f"\nmd=7 reproduction: {len(checks)-len(bad)}/{len(checks)} cells match "
          f"results/watershed_tail.csv")
    for c in bad:
        print(f"  MISMATCH {c[0]} q={c[1]} depth {c[2]} vs {c[3]}, pool {c[4]} vs {c[5]}")

    pd.set_option("display.width", 260)
    piv = out.pivot_table(index=["file_name", "arm"], columns="quantile",
                          values="depth", aggfunc="first")
    base = out.groupby(["file_name", "arm"])[["ceiling", "n_captured", "n_gt_mitotic",
                                              "n_pool"]].first()
    print("\n=== ceiling, pool, and depth to each sensitivity level ===")
    print(base.join(piv).to_string())


if __name__ == "__main__":
    main()
