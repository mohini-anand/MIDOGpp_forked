"""Two fixes, measured across all seven domain ROIs.

1. **Recall.** `FSConfig.score_threshold` is a constant in `TM_CCOEFF_NORMED` units, but
   that statistic's null distribution depends on template size, which varies 35-51 px
   between images. A fixed 0.5 therefore means a different depth of search on every ROI
   -- 246.tiff returned 5 detections against 115 mitotic figures while 301.tiff returned
   13,353. Replacing it with a floor in units of each map's own robust z (median and MAD,
   `template_match._robust_z`) makes one number mean the same thing everywhere.

2. **False positives.** Rank by chromatin density (`chromatin.rerank`) instead of by the
   correlation score, which is provably blind to it. See `midog_utils/chromatin.py`.

With a single augmentation `_robust_z` is a strictly monotone rescaling of the fused map,
so the peak set and their order are unchanged and only the *threshold* moves. This script
therefore matches once per ROI at a deep floor and derives every arm by filtering and
re-sorting that one candidate list -- which is exact, not an approximation, for the same
reason `Research Logs/design_choices.md` section 9 gives: greedy score-ordered NMS means a
peak above t can only ever be suppressed by a higher-scoring peak, which is also above t.

Writes results/od_experiment_{arms,workload,detections}.csv.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd

from midog_utils import channels as ch
from midog_utils import chromatin as cm
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from midog_utils import experiment as ex
from midog_utils import find_and_suppress as fs
from midog_utils import seed_selection as ss
from midog_utils import template_match as tm
from midog_utils.nms import nms_by_distance

EXTRACT_FLOOR = 0.25   # deep candidate floor, below every arm reported below
Z_LEVELS = (3.0, 2.5, 2.0, 1.5, 1.0)
CHANNEL = "rgb"        # the search channel; chromatin density is always read off hematoxylin


def robust_stats(fused: np.ndarray, stride: int = 8):
    """Median and MAD-scale of a response map -- the same statistics `_robust_z` uses."""
    s = fused[::stride, ::stride]
    s = s[s > -1.5]
    med = float(np.median(s))
    scale = float(1.4826 * np.median(np.abs(s - med)))
    return med, scale


def main():
    images, ann = ds.load_annotations()
    sel = ex.select_domain_images(images, ann)
    cfg = fs.FSConfig(channel=CHANNEL, score_threshold=EXTRACT_FLOOR)

    arms, workload, kept_det = [], [], []
    for _, row in sel.iterrows():
        fn = row["file_name"]
        t0 = time.time()
        path = f"images/{fn}"
        rgb = ds.load_roi(path)
        mpp = ds.roi_mpp(path)
        area = ds.check_roi_scale(path, rgb.shape)
        radius = ev.radius_px(mpp)

        gt = ds.image_annotations(ann, fn)
        gray_inv = ch.to_gray_inverted(rgb)
        # Per-image stream from (0, image_id) -- the same construction run_experiment
        # uses, so these seeds match the saved fs_simple_* runs. A fresh default_rng(0)
        # inside the loop is the bug the 2026-08-27 re-run fixed.
        rng = np.random.default_rng([0, int(row["image_id"])])
        seed, _ = ss.pick_seed(gt[gt["category_id"] == ds.MITOTIC], gray_inv,
                               rng, cfg.patch_size // 2, rgb.shape)
        gt_eval = gt[gt["ann_id"] != seed["ann_id"]].reset_index(drop=True)
        base = ss.tightened_base_size(gray_inv, seed["cx"], seed["cy"])
        run_cfg = fs.FSConfig(channel=CHANNEL, base_size=base, score_threshold=EXTRACT_FLOOR)

        # --- one match, one NMS, at the deep floor -------------------------------
        img = ch.to_channel(rgb, CHANNEL)
        patch = tm.read_padded_patch(img, seed["cx"], seed["cy"], run_cfg.patch_size)
        templates, _ = tm.build_augmentations(patch, base, run_cfg.scales,
                                              run_cfg.n_angles, run_cfg.flips)
        fused, _, valid = tm.fused_response(img, templates)
        med, scale = robust_stats(fused)
        centers, scores = tm.extract_peaks(fused, valid, run_cfg.peak_min_distance,
                                           EXTRACT_FLOOR, run_cfg.max_peaks)
        keep = nms_by_distance(centers, scores, radius)
        centers, scores = centers[keep], scores[keep]
        d_seed = np.hypot(centers[:, 0] - seed["cx"], centers[:, 1] - seed["cy"])
        centers, scores = centers[d_seed > run_cfg.self_hit_radius], scores[d_seed > run_cfg.self_hit_radius]

        pool = pd.DataFrame({"rank": np.arange(len(centers)), "cx": centers[:, 0],
                             "cy": centers[:, 1], "score": scores})
        hem = cm.hematoxylin_od(rgb)   # unclipped; to_hematoxylin saturates dense chromatin
        pool = cm.score_detections(pool, hem)
        n_mit = int((gt_eval["category_id"] == ds.MITOTIC).sum())
        print(f"{fn}: base={base}px  pool={len(pool)}  med={med:.3f} mad_scale={scale:.3f}  "
              f"z(0.5)={(0.5 - med) / scale:.2f}  z(floor {EXTRACT_FLOOR})="
              f"{(EXTRACT_FLOOR - med) / scale:.2f}  [{time.time() - t0:.0f}s]", flush=True)

        def record(arm, det, thr_note):
            det = det.reset_index(drop=True).assign(rank=np.arange(len(det)))
            det_out, _, m, _ = ev.evaluate_run(det, gt_eval, radius, area, roi_shape=rgb.shape)
            arms.append({"file_name": fn, "tumor_type": row["tumor_type"], "arm": arm,
                         "threshold": thr_note, "base_size": base,
                         "n_detections": len(det), "n_gt_mitotic": n_mit,
                         "mitotic_recall": m["mitotic_recall"], "recall_at_k": m["recall_at_k"],
                         "all_tp": m["all_tp"], "all_fp_lookalike": m["all_fp_lookalike"],
                         "all_fp_unannotated": m["all_fp_unannotated"],
                         "fp_per_tp": round(m["all_fp_unannotated"] / max(m["all_tp"], 1), 1),
                         "coverage_frac": m["coverage_frac"],
                         # True when the z cut sits below EXTRACT_FLOOR, i.e. the arm is
                         # limited by the extraction floor rather than by its own z level.
                         "floor_limited": floor_limited})
            return det_out

        # arm A: the current pipeline -- fixed 0.5 in correlation units
        floor_limited = False
        record("A_fixed_0.5", pool[pool.score >= 0.5].sort_values("score", ascending=False), "score>=0.5")

        for z in Z_LEVELS:
            cut = med + z * scale
            floor_limited = cut < EXTRACT_FLOOR
            sub = pool[pool.score >= cut]
            # arm B: per-image z floor, still ranked by correlation score
            record(f"B_z{z:g}", sub.sort_values("score", ascending=False), f"score>={cut:.3f}")
            # arm C: same candidates, ranked by chromatin density
            dc = record(f"C_z{z:g}_od", cm.rerank(sub), f"score>={cut:.3f}")

            # workload curve on arm C: how deep must a human read for X% of the mitoses?
            tp = (dc["bucket"] == ev.HUMAN_CORRECT_LABEL).to_numpy().cumsum()
            for frac in (0.5, 0.7, 0.8, 0.9, 1.0):
                want = int(np.ceil(frac * n_mit))
                i = int(np.searchsorted(tp, want))
                workload.append({"file_name": fn, "arm": f"C_z{z:g}_od", "target_recall": frac,
                                 "n_gt_mitotic": n_mit,
                                 "candidates_to_read": i + 1 if i < len(tp) else np.nan,
                                 "reachable": bool(i < len(tp))})
            if z == 2.5:
                kept_det.append(dc.assign(file_name=fn))

    pd.DataFrame(arms).to_csv("results/od_experiment_arms.csv", index=False)
    pd.DataFrame(workload).to_csv("results/od_experiment_workload.csv", index=False)
    pd.concat(kept_det).to_csv("results/od_experiment_detections.csv", index=False)
    print("\nwrote results/od_experiment_{arms,workload,detections}.csv")


if __name__ == "__main__":
    main()
