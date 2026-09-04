"""Split merged nuclei with a watershed, and re-measure the ceiling and the tail.

`miss_attribution.py` establishes that **all eleven** problem annotations -- the six mitotic
figures the proposal stage misses outright, and the five matched only at 26.7-28.4 px
against a ~30.2 px match radius -- sit inside a connected component that exists and is in
the 71st-98th percentile of chromatin density. They are not faint. Three sit in components
over the 4,000 px gate; the rest sit in components that pass the gate but whose *centroid*
is displaced, because `label(connectivity=2)` fuses touching nuclei into one object.

That single defect explains the ceiling *and* the cliff: 246.tiff's two deep-tail
annotations (ranks 17,320 and 20,492) are fused components whose mean hematoxylin is diluted
by whatever they merged with. A re-ranker cannot fix either, because the object that should
have been ranked was never proposed.

This measures whether a distance-transform watershed fixes it, on three numbers: the
proposal **ceiling**, the **depth** to each sensitivity level, and the **candidate-pool
size**. The third is not optional -- an aggressive split doubles the pool, which improves
`frac_of_pool` while making the reader's job worse, so depth is compared in absolute
candidates.

Splitting is done **per component inside its own bounding box**, and only for components
large enough to hold two nuclei (>= 2 x `MIN_AREA`). A full-image distance transform,
`peak_local_max` and watershed over a 5412 x 7215 ROI costs 5-10 minutes each; restricting
the work to the foreground objects that could actually be merged costs seconds, and leaves
components that are already single nuclei bit-identical to the connected-components arm.

Writes results/watershed_tail.csv.
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
from midog_utils import chromatin as cm
from midog_utils import dataset as ds
from midog_utils import evaluate as ev

ROIS = ["301.tiff", "300.tiff", "246.tiff", "245.tiff", "201.tiff", "202.tiff", "405.tiff"]
QUANTS = [0.50, 0.90, 0.95, 0.98, 0.99, 1.00]
MIN_AREA, MAX_AREA = 80, 4000
# Half the diameter of the smallest component worth separating. An 80 px blob is ~10 px
# across, so two EDT maxima closer than this belong to one nucleus and must not seed two
# basins. This is the knob that trades ceiling against pool size.
MIN_PEAK_DIST = 7


def split_components(binary: np.ndarray, h_chan: np.ndarray) -> pd.DataFrame:
    """Connected components, with any component big enough to be two nuclei watershedded."""
    lab = label(binary, connectivity=2)
    rows = []
    for p in regionprops(lab, intensity_image=h_chan):
        minr, minc, maxr, maxc = p.bbox
        if p.area < 2 * MIN_AREA:
            rows.append((p.centroid[1], p.centroid[0], p.area, p.intensity_mean))
            continue
        sub = p.image
        d = ndi.distance_transform_edt(sub)
        coords = peak_local_max(d, min_distance=MIN_PEAK_DIST, labels=sub,
                                exclude_border=False)
        if len(coords) <= 1:
            rows.append((p.centroid[1], p.centroid[0], p.area, p.intensity_mean))
            continue
        markers = np.zeros(d.shape, dtype=np.int32)
        markers[tuple(coords.T)] = np.arange(1, len(coords) + 1)
        ws = watershed(-d, markers, mask=sub)
        icrop = h_chan[minr:maxr, minc:maxc]
        for q in regionprops(ws, intensity_image=icrop):
            rows.append((q.centroid[1] + minc, q.centroid[0] + minr,
                         q.area, q.intensity_mean))
    return pd.DataFrame(rows, columns=["cx", "cy", "area", "score"])


def gate_and_rank(df: pd.DataFrame, max_area: int = MAX_AREA) -> pd.DataFrame:
    out = df[(df["area"] >= MIN_AREA) & (df["area"] <= max_area)]
    return out.sort_values("score", ascending=False).reset_index(drop=True)


def tail_rows(det, gt, radius, tag, fn):
    _d, gt_out = ev.bucket_detections(det.reset_index(drop=True),
                                      gt.reset_index(drop=True), radius)
    mit = gt_out[gt_out["category_id"] == ds.MITOTIC]
    ranks = np.sort(mit.loc[mit["found"], "matched_rank"].to_numpy())
    n_gt, n_found, n_pool = len(mit), len(ranks), len(det)
    rows = []
    for q in QUANTS:
        if not n_found:
            continue
        i = int(np.ceil(q * n_found)) - 1
        rows.append(dict(file_name=fn, arm=tag, n_gt_mitotic=n_gt, n_captured=n_found,
                         ceiling=round(n_found / n_gt, 4), n_pool=n_pool, quantile=q,
                         depth=int(ranks[i]) + 1,
                         frac_of_pool=round((int(ranks[i]) + 1) / n_pool, 4)))
    return rows


def main():
    _images, ann = ds.load_annotations()
    rows = []
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

        lab = label(binary, connectivity=2)
        cc = pd.DataFrame(
            [(p.centroid[1], p.centroid[0], p.area, p.intensity_mean)
             for p in regionprops(lab, intensity_image=h_chan)],
            columns=["cx", "cy", "area", "score"])
        del lab
        rows += tail_rows(gate_and_rank(cc), gt, radius, "connected_components", fn)
        # isolates "the split worked" from "the upper area gate was the only problem"
        rows += tail_rows(gate_and_rank(cc, 10000), gt, radius, "cc_maxarea10k", fn)

        ws = split_components(binary, h_chan)
        rows += tail_rows(gate_and_rank(ws), gt, radius, "watershed", fn)
        rows += tail_rows(gate_and_rank(ws, 10000), gt, radius, "watershed_maxarea10k", fn)
        print(f"{fn}: {time.time()-t0:.0f}s", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv("results/watershed_tail.csv", index=False)
    pd.set_option("display.width", 260)
    piv = out.pivot_table(index=["file_name", "arm"], columns="quantile",
                          values="depth", aggfunc="first")
    base = out.groupby(["file_name", "arm"])[["ceiling", "n_captured", "n_gt_mitotic",
                                              "n_pool"]].first()
    print("\n=== ceiling, pool size, and depth to each sensitivity level ===")
    print(base.join(piv).to_string())


if __name__ == "__main__":
    main()
