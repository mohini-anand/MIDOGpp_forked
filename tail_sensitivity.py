"""What does 100% sensitivity actually cost? The tail statistic nobody has measured.

Every metric this project has used so far -- recall@budget, read-to-50% -- is a *median*.
The product requirement "we definitely want to hit all the mitotic cells" is a **tail**
statistic: the reading depth at which the *worst-ranked* true mitotic figure appears.
AUC 0.986 says nothing about it. A single mitosis at rank 14,000 of 20,000 makes 100%
sensitivity cost the entire slide while leaving AUC untouched.

This measures, per (ROI, ranker), on the seedless `nucleus_blobs` candidate set:

* **the proposal ceiling** -- what fraction of annotated mitoses are in the candidate set
  at all. No ranking can exceed this, and on 246.tiff it is already 96.6%.
* **the depth to reach 50 / 90 / 95 / 99 / 100% of the mitoses the proposals captured**,
  in candidates a reader must work through, and as a fraction of the pool.
* how many of those candidates are annotated look-alikes vs unannotated.

Rankers compared, all seedless and all already computed by the existing code:
`blob_score` (the Otsu component's mean hematoxylin), `chromatin_od` (mean of the darkest
10% of a 51 px window on unclipped hematoxylin OD), and `rankfuse` (mean of the two ranks).

Also diagnoses every *missed* mitosis -- the ones the proposal stage never emits -- against
the pre-area-filter component map, so the cause is attributable to the area filter, the
tissue mask, or a merge with a neighbour.

Writes results/tail_sensitivity.csv and results/tail_misses.csv.
"""

from __future__ import annotations

import time

import cv2
import numpy as np
import pandas as pd
from skimage.measure import label, regionprops

from midog_utils import baselines as bl
from midog_utils import chromatin as cm
from midog_utils import dataset as ds
from midog_utils import evaluate as ev

# Every ROI on disk with >= 14 annotated mitotic figures. The prior 7-ROI set omitted
# 300.tiff (181) and 245.tiff (90), the second and fourth densest available -- four of the
# seven it did use have n < 15, which the 2026-09-01 audit flagged as too thin to settle a
# question this close.
ROIS = ["301.tiff", "300.tiff", "246.tiff", "245.tiff", "201.tiff", "202.tiff", "405.tiff"]

QUANTS = [0.50, 0.90, 0.95, 0.98, 0.99, 1.00]


def all_components(rgb: np.ndarray, tile: int = 512):
    """`nucleus_blobs`' component map with **no** area filter, for miss attribution."""
    mask = bl.tissue_mask(rgb)
    h_chan = bl.to_hematoxylin(rgb)
    thr = bl._tiled_otsu_threshold(h_chan, mask, tile)
    binary = (h_chan > thr) & mask
    lab = label(binary, connectivity=2)
    rows = [{"cx": p.centroid[1], "cy": p.centroid[0], "area": int(p.area)}
            for p in regionprops(lab)]
    return pd.DataFrame(rows), mask


def depth_table(det: pd.DataFrame, gt: pd.DataFrame, radius: float, key: str):
    """Sort by ``key``, match, and report the depth to reach each sensitivity quantile.

    ``bucket_detections`` returns ``gt_out.matched_rank``: the position, in the order the
    detections were handed to it, of the detection that claimed each annotation. Sorting
    by ``key`` first therefore makes ``matched_rank`` literally "how far down this ranking
    the reader must go to see this mitosis".
    """
    order = det.sort_values(key, ascending=False, na_position="last", kind="mergesort")
    order = order.reset_index(drop=True)
    det_out, gt_out = ev.bucket_detections(order, gt.reset_index(drop=True), radius)

    mit = gt_out[gt_out["category_id"] == ds.MITOTIC]
    found = mit[mit["found"]]
    ranks = np.sort(found["matched_rank"].to_numpy())
    n_gt, n_found, n_pool = len(mit), len(ranks), len(order)

    buckets = det_out["bucket"].to_numpy()
    rows = []
    for q in QUANTS:
        if n_found == 0:
            continue
        # depth counted in candidates read, so rank 0 costs one read
        i = int(np.ceil(q * n_found)) - 1
        depth = int(ranks[i]) + 1
        top = buckets[:depth]
        rows.append(dict(
            quantile=q,
            sens_of_captured=round((i + 1) / n_found, 4),
            sens_of_all_gt=round((i + 1) / n_gt, 4),
            depth=depth,
            frac_of_pool=round(depth / n_pool, 4),
            fp_lookalike=int((top == ev.HUMAN_REJECTED_LABEL).sum()),
            fp_unannotated=int((top == ev.NON_HUMAN_FINDINGS).sum()),
            fp_per_tp=round((depth - (i + 1)) / (i + 1), 1),
        ))
    return rows, n_gt, n_found, n_pool, mit


def main():
    _images, ann = ds.load_annotations()
    tail_rows, miss_rows = [], []

    for fn in ROIS:
        t0 = time.time()
        path = f"images/{fn}"
        rgb = ds.load_roi(path)
        radius = ev.radius_px(ds.roi_mpp(path))
        gt = ds.image_annotations(ann, fn)
        hem_od = cm.hematoxylin_od(rgb)

        blobs = bl.nucleus_blobs(rgb)
        det = cm.score_detections(blobs, hem_od)
        # rank fusion: mean of the two 0-based ranks, negated so "higher is better" holds
        r_blob = det["score"].rank(ascending=False, method="average")
        r_od = det["od"].rank(ascending=False, method="average", na_option="bottom")
        det = det.assign(rankfuse=-(r_blob + r_od) / 2.0)

        for key in ("blob_score", "chromatin_od", "rankfuse"):
            col = {"blob_score": "score", "chromatin_od": "od", "rankfuse": "rankfuse"}[key]
            rows, n_gt, n_found, n_pool, mit = depth_table(det, gt, radius, col)
            for r in rows:
                tail_rows.append(dict(file_name=fn, ranker=key, n_gt_mitotic=n_gt,
                                      n_captured=n_found, n_pool=n_pool,
                                      ceiling=round(n_found / n_gt, 4), **r))
            if key == "blob_score":  # miss set is ranker-independent
                missed = mit[~mit["found"]]
                if len(missed):
                    comps, tmask = all_components(rgb)
                    cxy = comps[["cx", "cy"]].to_numpy()
                    for _, m in missed.iterrows():
                        d = np.hypot(cxy[:, 0] - m["cx"], cxy[:, 1] - m["cy"])
                        j = int(np.argmin(d))
                        inside = bool(tmask[int(round(m["cy"])), int(round(m["cx"]))])
                        area = int(comps.iloc[j]["area"])
                        cause = ("outside_tissue_mask" if not inside else
                                 "no_component_within_radius" if d[j] > radius else
                                 "area_below_min" if area < 80 else
                                 "area_above_max" if area > 4000 else "unknown")
                        miss_rows.append(dict(
                            file_name=fn, ann_id=int(m["ann_id"]),
                            cx=round(float(m["cx"]), 1), cy=round(float(m["cy"]), 1),
                            nearest_comp_dist=round(float(d[j]), 1),
                            nearest_comp_area=area, match_radius=round(radius, 1),
                            in_tissue_mask=inside, cause=cause))

        print(f"{fn}: {len(det)} candidates, {int((gt.category_id==ds.MITOTIC).sum())} "
              f"mitotic, {time.time()-t0:.0f}s", flush=True)

    tail = pd.DataFrame(tail_rows)
    tail.to_csv("results/tail_sensitivity.csv", index=False)
    miss = pd.DataFrame(miss_rows)
    miss.to_csv("results/tail_misses.csv", index=False)

    pd.set_option("display.width", 250)
    print("\n=== depth to reach each sensitivity level (candidates read) ===")
    piv = tail.pivot_table(index=["file_name", "ranker"], columns="quantile",
                           values="depth", aggfunc="first")
    ceil = tail.groupby(["file_name", "ranker"])[["ceiling", "n_pool", "n_gt_mitotic"]].first()
    print(ceil.join(piv).to_string())
    print("\n=== misses (proposal stage, ranker-independent) ===")
    print(miss.to_string(index=False) if len(miss) else "none")


if __name__ == "__main__":
    main()
