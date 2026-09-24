"""
    Extras outside V1-V10: whether load_annotations' `unanimous` column agrees with agreement_pool's n_mitotic_votes ==
    n_votes on the contested-tier ROI, and one positive multi-row run_roi (two HUMAN rows, seed chosen by seed_object_id).
"""

from __future__ import annotations

import json

import common as C
from common import ds, rr
import checks
from midog_utils import annotation_frame as af


def main():
    """
        Run both extras and write out/full/extra.json.

        Returns None.
    """
    rois, annotations = rr.build_rois()
    m = ds.image_annotations(annotations, "201.tiff", category_id=ds.MITOTIC)
    agree = int((m["unanimous"] == (m["n_mitotic_votes"] == m["n_votes"])).sum())
    out = dict(unanimous_vs_votes_201=dict(n_mitotic=len(m), n_agree=agree))

    roi = rois[rois["file_name"] == "300.tiff"].iloc[0]
    rgb, mpp = ds.load_roi(roi["path"]), ds.roi_mpp(roi["path"])
    gt = ds.image_annotations(annotations, "300.tiff")
    a, b = (gt[gt["ann_id"] == i].iloc[0] for i in (14499, 14578))
    single = af.human_seed_frame(a, "300.tiff")
    two = af.human_seed_frame(b, "300.tiff").pipe(lambda f: __import__("pandas").concat([f, single], ignore_index=True))
    f1, i1 = af.run_roi(rgb, single, mpp, refinement_channel="gray_inverted")
    f2, i2 = af.run_roi(rgb, two, mpp, seed_object_id=14499, refinement_channel="gray_inverted")
    p1, p2 = f1[f1["class"] == "PREDICTION"], f2[f2["class"] == "PREDICTION"]
    same_preds = C.compare_frames(p2, p1, checks.COMPARED_FRAME_COLUMNS, label="multi-row vs single")
    ref1, ref2 = f1.iloc[1].drop("current_time"), f2.iloc[2].drop("current_time")
    out["multi_row"] = dict(n_rows=len(f2), expected=2 + 1 + i2["n_detections"], input_rows_first=list(f2["object_id"].iloc[:2]), refined_row_equal=bool(ref1.equals(ref2)), seed_object_id=i2["seed_object_id"], click_center=i2["click_center"], preds_values=same_preds["n_values"], preds_diverged=same_preds["n_diverged"])
    (C.OUT / "full" / "extra.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
