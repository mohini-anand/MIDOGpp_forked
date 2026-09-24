"""
    Step 0 gate: B1 must reproduce R1 (per-run) and R2 (top 30), and B2's default_51 runs at seeds 0-2 must reproduce R4's
    tm_score top 30. R1 and R2 rows are rebuilt from the raw B1 pickles with the driver's own code (rr.score_top_k for buckets
    and tp_at_k, rerun_bbox3way_postD11.py:233-252 for the rest). Compares only the ROIs whose pickles exist, against the same
    ROIs' reference rows. Usage: step0_compare.py [--require-all]
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

import common as C
from common import ds, ev, rr


def main():
    """
        Rebuild R1/R2 from B1 and R4's tm_score top 30 from B2, compare exactly, and write the summary.

        Returns None.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-all", action="store_true", help="fail unless all 49 ROIs are captured")
    args = parser.parse_args()

    rois, _ = rr.build_rois()
    have = [fn for fn in rois["file_name"] if (C.STEP0 / f"{fn}.pkl").exists()]
    if args.require_all:
        assert len(have) == 49, f"only {len(have)}/49 ROIs captured"
    r1, r2 = C.read_ref(C.R1_PATH), C.read_ref(C.R2_PATH)
    r4_top = C.read_ref(C.R4_TOP30_PATH)

    rows, tops, b2_mine, b2_ref = [], [], [], []
    for _, roi in rois[rois["file_name"].isin(have)].iterrows():
        d = C.load_pickle(C.STEP0 / f"{roi['file_name']}.pkl")
        for seed_index, ann, specs, refusals, n_seed_pool, digests in d["draws"]:
            for condition in specs:
                detections, info = d["runs"][(seed_index, condition, "chromatin_od")]
                row, top = C.driver_row(roi, seed_index, ann, specs, refusals, n_seed_pool, digests, condition, detections, info, d["gt"], d["pool"], d["flagged"], d["n_agreement"], d["n_joint_valid"], d["mpp"])
                rows.append(row)
                tops.append(top)
            if roi["subset"] == "original_14":
                detections, _ = d["runs"][(seed_index, "default_51", "tm_score")]
                gt_eval = d["gt"][d["gt"]["ann_id"] != int(ann["ann_id"])].reset_index(drop=True)
                det_out, _ = rr.score_top_k(detections, gt_eval, ev.radius_px(d["mpp"]))
                b2_mine.append(det_out.head(C.TOP_N).assign(file_name=roi["file_name"], seed_index=seed_index))
                ref = r4_top[(r4_top["file_name"] == roi["file_name"]) & (r4_top["seed_index"] == seed_index) & (r4_top["condition"] == "default_51") & (r4_top["rank_key"] == "tm_score")]
                b2_ref.append(ref)

    per_run = C.finish_per_run(rows, rois)
    r1_sub = r1[r1["file_name"].isin(have)].reset_index(drop=True)
    r1_cols = [c for c in r1.columns if not c.startswith("t_")]
    assert list(per_run.columns) == r1_cols, f"rebuilt R1 columns differ:\n{list(per_run.columns)}\n{r1_cols}"
    key = ["file_name", "seed_index", "condition"]
    assert per_run[key].equals(r1_sub[key]), "R1 key order differs"
    res_r1 = C.compare_frames(per_run, r1_sub, r1_cols, label="B1 vs R1")

    top30 = pd.concat(tops, ignore_index=True)
    r2_sub = r2[r2["file_name"].isin(have)].reset_index(drop=True)
    assert list(top30.columns) == list(r2.columns), f"rebuilt R2 columns differ:\n{list(top30.columns)}\n{list(r2.columns)}"
    key2 = ["file_name", "seed_index", "condition", "rank"]
    assert top30[key2].reset_index(drop=True).equals(r2_sub[key2]), "R2 key order differs"
    res_r2 = C.compare_frames(top30, r2_sub, list(r2.columns), float32_cols=("score",), label="B1 vs R2")

    res_b2 = None
    if b2_mine:
        b2m = pd.concat(b2_mine, ignore_index=True)
        b2r = pd.concat(b2_ref, ignore_index=True)
        b2_cols = ["file_name", "seed_index", "rank", "cx", "cy", "score", "bucket", "matched_ann_id", "matched_category"]
        res_b2 = C.compare_frames(b2m, b2r, b2_cols, float32_cols=("score",), label="B2 default_51 vs R4 tm_score top30")
        res_b2["n_runs"] = len(b2_mine)

    summary = dict(n_rois=len(have), rois=have, r1=dict(n_runs=len(per_run), **res_r1), r2=res_r2, b2_vs_r4=res_b2)
    name = "step0_compare_all.json" if len(have) == 49 else f"step0_compare_{len(have)}rois.json"
    (C.OUT / name).write_text(json.dumps(summary, indent=1, default=str))
    for tag in ("r1", "r2", "b2_vs_r4"):
        r = summary[tag]
        if r is None:
            continue
        print(f"{tag}: rows={r['n_rows']} columns={r['n_columns']} values={r['n_values']} diverged={r['n_diverged']}")
        for line in r["detail"]:
            print("   ", line)
    print(f"-> {C.OUT / name}")


if __name__ == "__main__":
    main()
