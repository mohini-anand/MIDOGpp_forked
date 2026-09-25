"""
    Parity gate: the precomputed-channel search in common.py must reproduce production and the recorded post-D11 runs
    before any saliency number is trusted.

    Checks, on the ROIs given (default 300.tiff and 128.tiff):
      1. common.search + rank_by == production.run_production_pipeline, both axes, one click (every column).
      2. default_51 / gray_bbox / hem_bbox at chromatin_od: base_size, template centre, tpl_sha1, n_peaks, n_detections,
         n_blanked_px and TP@10/20/30 equal the bbox3way postD11 per_run rows.
      3. default_51 at tm_score: TP@10/20/30 equal the fillgeom postD11 per_run rows (6 seeds).
      4. common.gray_component's bbox equals tighten_box_otsu's on every border-filtered pool annotation.
      5. common.production_spec at window 51 equals the gray_bbox spec.
    Usage: parity.py [ROI ...]
"""

import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

import common as C


def main(fns):
    """
        Run every parity check on the named ROIs; raises AssertionError on the first mismatch.

        fns (list[str]): ROI file names.

        Returns None.
    """
    rois = C.roi_table().set_index("file_name")
    b3 = pd.read_csv(C.RESULTS / "precision_at_k_49roi_3seed_chromatin_bbox3way_postD11_per_run.csv")
    fg = pd.read_csv(C.RESULTS / "precision_at_k_49roi_6seed_fillgeom_postD11_per_run.csv")
    channel_of = {"default_51": None, "gray_bbox": "gray", "hem_bbox": "hem"}
    for fn in fns:
        roi = C.load_roi(rois.loc[fn, "path"])
        gt = roi["gt"]
        # 1. direct equality with production on the first recorded click
        first = b3[(b3.file_name == fn) & (b3.condition == "gray_bbox")].iloc[0]
        seed = SimpleNamespace(base_size=int(first.base_size), template_xy=(float(first.tpl_cx), float(first.tpl_cy)))
        det, _ = C.search(roi, seed.base_size, *seed.template_xy)
        for rank_key in ("tm_score", "chromatin_od"):
            ref, _ = C.prod.run_production_pipeline(roi["rgb"], seed, roi["mpp"], rank_key=rank_key, always_od=True)
            mine = C.rank_by(det, rank_key)
            pd.testing.assert_frame_equal(mine[ref.columns].reset_index(drop=True), ref.reset_index(drop=True), check_exact=True)
        print(f"{fn}: common.search == run_production_pipeline on both axes")

        # 2. bbox3way postD11 at chromatin_od
        n2 = 0
        for _, r in b3[b3.file_name == fn].iterrows():
            cx, cy = float(r.click_cx), float(r.click_cy)
            name = channel_of[r.condition]
            got = (51, cx, cy) if name is None else C.ss.tightened_template_box(roi[name], cx, cy, otsu_window=51)
            assert (int(got[0]), float(got[1]), float(got[2])) == (int(r.base_size), float(r.tpl_cx), float(r.tpl_cy)), (fn, r.condition, got)
            assert C.template_sha1(roi["hem"], *got) == r.tpl_sha1
            det, info = C.search(roi, *got)
            s = C.score(C.rank_by(det, "chromatin_od"), gt, int(r.seed_ann_id), (cx, cy), roi["mpp"])
            for k in C.BUDGETS:
                assert s[f"tp_at_{k}"] == int(r[f"tp_at_{k}"]), (fn, r.condition, k, s, r[f"tp_at_{k}"])
            for c in ("n_peaks", "n_detections", "n_blanked_px"):
                assert int(info[c]) == int(r[c]), (fn, r.condition, c)
            if r.condition == "gray_bbox":
                p = C.production_spec(roi["gray"], cx, cy, 51, roi["hem"].shape)
                assert (p["base_size"], p["tx"], p["ty"]) == (int(r.base_size), float(r.tpl_cx), float(r.tpl_cy)) and p["status"] == "refined"
            n2 += 1
        print(f"{fn}: {n2} bbox3way postD11 rows reproduced (spec, sha1, counts, TP@10/20/30)")

        # 3. fillgeom default_51 at tm_score
        n3 = 0
        for _, r in fg[(fg.file_name == fn) & (fg.condition == "default_51") & (fg.rank_key == "tm_score")].iterrows():
            det, _ = C.search(roi, 51, float(r.click_cx), float(r.click_cy))
            s = C.score(C.rank_by(det, "tm_score"), gt, int(r.seed_ann_id), (float(r.click_cx), float(r.click_cy)), roi["mpp"])
            for k in C.BUDGETS:
                assert s[f"tp_at_{k}"] == int(r[f"tp_at_{k}"]), (fn, "fillgeom", k)
            assert s["n_near_click"] == int(r.n_near_click)
            n3 += 1
        print(f"{fn}: {n3} fillgeom default_51 tm_score rows reproduced")

        # 4. gray_component == tighten_box_otsu on the pool
        pool, _ = C.ss.agreement_pool(gt[gt.category_id == C.ds.MITOTIC])
        pool = C.ss.border_filter(pool, C.HALF_PATCH, roi["hem"].shape)
        n4 = 0
        for _, a in pool.iterrows():
            patch = C.tm.read_padded_patch(roi["gray"], float(a.cx), float(a.cy), 51)
            mine = C.gray_component(patch)
            ref = C.ss.tighten_box_otsu(patch)
            assert (mine is None and ref is None) or (mine is not None and ref is not None and mine[1] == ref), (fn, a.ann_id)
            n4 += 1
        print(f"{fn}: gray_component == tighten_box_otsu on {n4} pool annotations")


if __name__ == "__main__":
    main(sys.argv[1:] or ["300.tiff", "128.tiff"])
    print("PARITY OK")
