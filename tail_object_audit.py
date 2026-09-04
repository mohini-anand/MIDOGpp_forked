"""Are the proposal-stage misses and the deep-tail mitoses the same failure mode?

`tail_sensitivity.py` reports two separate-looking problems on 246.tiff: 4 mitoses the
proposal stage never emits (nearest component 1-17 px, i.e. tiled Otsu found almost no
foreground), and a cliff where the 111th and 112th of 112 captured mitoses sit at ranks
~17k and ~20k. If both populations are *faint* objects, that is **one** failure mode, and
the fix -- a proposal/ranking signal that is not chromatin intensity -- is a single
intervention rather than two.

For the 8 deepest captured mitoses per ROI this reports the rank, the percentile of the
matched detection's `blob_score` and area, the raw chromatin OD, and the distance from
annotation to matched detection. The distance column also catches the greedy-match
artifact: `bucket_detections` assigns in the order it is handed, so two annotations closer
than the match radius can swap claimants under a different sort key. With a headline that
turns on two objects, that is worth ruling out explicitly.

Writes results/tail_object_audit.csv.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from midog_utils import baselines as bl
from midog_utils import chromatin as cm
from midog_utils import dataset as ds
from midog_utils import evaluate as ev

ROIS = ["301.tiff", "300.tiff", "246.tiff", "245.tiff"]
N_DEEPEST = 8


def main():
    _images, ann = ds.load_annotations()
    rows = []
    for fn in ROIS:
        path = f"images/{fn}"
        rgb = ds.load_roi(path)
        radius = ev.radius_px(ds.roi_mpp(path))
        gt = ds.image_annotations(ann, fn).reset_index(drop=True)
        hem_od = cm.hematoxylin_od(rgb)

        blobs = bl.nucleus_blobs(rgb)          # already sorted by component mean, desc
        det = cm.score_detections(blobs, hem_od).reset_index(drop=True)
        det_out, gt_out = ev.bucket_detections(det, gt, radius)

        # percentile ranks within this ROI's pool, so values are comparable across ROIs
        n = len(det)
        pct_score = det["score"].rank(pct=True).to_numpy()
        pct_area = det["area"].rank(pct=True).to_numpy()
        pct_od = det["od"].rank(pct=True, na_option="bottom").to_numpy()

        mit = gt_out[gt_out["category_id"] == ds.MITOTIC]
        found = mit[mit["found"]].sort_values("matched_rank")
        deep = found.tail(N_DEEPEST)
        for _, m in deep.iterrows():
            k = int(m["matched_rank"])
            d = float(np.hypot(det.at[k, "cx"] - m["cx"], det.at[k, "cy"] - m["cy"]))
            rows.append(dict(
                file_name=fn, kind="deep_captured", ann_id=int(m["ann_id"]),
                rank=k, rank_pct=round(k / n, 4),
                blob_score=round(float(det.at[k, "score"]), 4),
                blob_score_pct=round(float(pct_score[k]), 4),
                area=int(det.at[k, "area"]), area_pct=round(float(pct_area[k]), 4),
                od51=round(float(det.at[k, "od"]), 4), od51_pct=round(float(pct_od[k]), 4),
                ann_to_det_px=round(d, 1), match_radius=round(radius, 1)))

        # the misses, scored the same way at the annotation's own coordinates
        for _, m in mit[~mit["found"]].iterrows():
            od = cm.chromatin_density(hem_od, float(m["cx"]), float(m["cy"]))
            rows.append(dict(
                file_name=fn, kind="missed", ann_id=int(m["ann_id"]),
                rank=-1, rank_pct=np.nan, blob_score=np.nan, blob_score_pct=np.nan,
                area=-1, area_pct=np.nan, od51=round(float(od), 4),
                od51_pct=round(float((det["od"] < od).mean()), 4),
                ann_to_det_px=np.nan, match_radius=round(radius, 1)))

        # reference: the median captured mitosis, for scale
        med = found.iloc[len(found) // 2]
        k = int(med["matched_rank"])
        rows.append(dict(
            file_name=fn, kind="median_captured", ann_id=int(med["ann_id"]),
            rank=k, rank_pct=round(k / n, 4),
            blob_score=round(float(det.at[k, "score"]), 4),
            blob_score_pct=round(float(pct_score[k]), 4),
            area=int(det.at[k, "area"]), area_pct=round(float(pct_area[k]), 4),
            od51=round(float(det.at[k, "od"]), 4), od51_pct=round(float(pct_od[k]), 4),
            ann_to_det_px=np.nan, match_radius=round(radius, 1)))
        print(f"{fn} done", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv("results/tail_object_audit.csv", index=False)
    pd.set_option("display.width", 260)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
