"""Does the reference's positional NMS ordering cost anything?

`bbox tuning code reference/bbox_tuning.py:483` (`nms_with_area`) resolves an
overlapping cluster by keeping whichever box sits **furthest right**: the frame is
sorted by ``['humanMade', 'nms applied', 'x top left', 'y top left']`` at
bbox_tuning.py:785 and the loop then takes ``last`` -- the largest ``x top left`` -- and
suppresses the earlier (smaller-x) members. The match score is never consulted, and is
dropped from the frame entirely at bbox_tuning.py:798.

`midog_utils/nms.py` instead suppresses in **descending score** order. This script
measures what that choice is worth, holding everything else fixed: same fused
correlation map, same extracted peaks, same suppression radius (each image's evaluation
match radius), same self-hit removal. Only the order the greedy loop walks in changes.

One arm is emitted, ``x_then_rank_by_score``: survivors chosen by position, then
**ranked by score** for evaluation. That is the isolating experiment -- it charges the
reference rule only for choosing the wrong cluster representative, and keeps every metric
well defined.

What the reference *actually* produces -- survivors chosen by position and emitted in x
order with no score at all -- is deliberately not scored, and no row for it exists in the
output. ``recall@K``, FROC and ``sens@Xfp`` all consume a ranked list, so applying them to
an x-ordered one would not yield a worse number, it would answer a different question. Its
detection list is the same set as ``x_then_rank_by_score``'s, by construction, differing
only in output order -- so there is nothing separate to report.

Run: PYTHONPATH=. python nms_ordering_probe.py
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from sklearn.neighbors import KDTree

from midog_utils import channels as ch
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from midog_utils import experiment as ex
from midog_utils import find_and_suppress as fs
from midog_utils import seed_selection as ss
from midog_utils import template_match as tm

CHANNEL = "rgb"
SCORE_THRESHOLD = 0.5      # the adopted floor (design_choices.md section 6/9)


def nms_by_x(centers: np.ndarray, radius: float) -> np.ndarray:
    """The reference's ordering, ported to points: greedy in descending x (then y).

    Mirrors `nms_with_area`'s loop exactly in the only respect that differs from
    `nms.nms_by_distance` -- the order the survivors are picked in. The suppression
    criterion is the same distance test, deliberately: changing the geometry too would
    confound ordering with criterion and answer neither question.
    """
    centers = np.asarray(centers, dtype=np.float64)
    if len(centers) == 0:
        return np.zeros(0, dtype=int)
    neighbours = KDTree(centers).query_radius(centers, r=radius)
    # ascending sort on (x, y), then take the last -- i.e. descending (x, y).
    order = np.lexsort((centers[:, 1], centers[:, 0]))[::-1]
    suppressed = np.zeros(len(centers), dtype=bool)
    keep = []
    for idx in order:
        if suppressed[idx]:
            continue
        keep.append(idx)
        suppressed[neighbours[idx]] = True
    return np.asarray(keep, dtype=int)


def survivor_deficit(centers, scores, keep_idx, radius):
    """How much score each survivor gave up against the best peak in its own radius.

    For every kept point, the highest-scoring *original* peak within ``radius`` of it is
    the representative a score-ordered rule would have preferred locally. Returns the
    per-survivor score gap and the distance to that better peak.
    """
    if len(keep_idx) == 0:
        return np.zeros(0), np.zeros(0)
    kept = centers[keep_idx]
    cands = KDTree(centers).query_radius(kept, r=radius)
    deficit, displacement = [], []
    for k, c in zip(keep_idx, cands):
        best = c[int(np.argmax(scores[c]))]
        deficit.append(float(scores[best] - scores[k]))
        displacement.append(float(np.hypot(*(centers[best] - centers[k]))))
    return np.asarray(deficit), np.asarray(displacement)


def main():
    images, annotations = ds.load_annotations("databases/MIDOG++.json")
    selection = ex.select_domain_images(images, annotations)
    cfg = fs.FSConfig(channel=CHANNEL, score_threshold=SCORE_THRESHOLD)

    rows, flip_rows = [], []
    for _, sel in selection.iterrows():
        fn = sel["file_name"]
        t0 = time.time()
        path = f"images/{fn}"
        rgb = ds.load_roi(path)
        mpp = ds.roi_mpp(path)
        area_mm2 = ds.check_roi_scale(path, rgb.shape)
        radius = ev.radius_px(mpp)

        gt = ds.image_annotations(annotations, fn)
        gray_inv = ch.to_gray_inverted(rgb)
        rng = np.random.default_rng([0, int(sel["image_id"])])
        seed, _ = ss.pick_seed(gt[gt["category_id"] == ds.MITOTIC], gray_inv, rng,
                               cfg.patch_size // 2, rgb.shape)
        gt_eval = gt[gt["ann_id"] != seed["ann_id"]].reset_index(drop=True)
        base = ss.tightened_base_size(gray_inv, seed["cx"], seed["cy"])

        img = ch.to_channel(rgb, cfg.channel)
        del gray_inv
        patch = tm.read_padded_patch(img, seed["cx"], seed["cy"], cfg.patch_size)
        templates, _ = tm.build_augmentations(patch, base, cfg.scales, cfg.n_angles, cfg.flips)
        fused, _, valid = tm.fused_response(img, templates, cfg.scale_normalize)
        del img
        centers, scores = tm.extract_peaks(fused, valid, cfg.peak_min_distance,
                                           cfg.score_threshold, cfg.max_peaks)
        del fused, valid
        n_peaks = len(centers)

        out = {}
        for name, keep in (("score", fs.nms_by_distance(centers, scores, radius)),
                           ("x_then_rank_by_score", nms_by_x(centers, radius))):
            c, s = centers[keep], scores[keep]
            deficit, displacement = survivor_deficit(centers, scores, keep, radius)
            # same self-hit removal as the pipeline
            self_hit = np.hypot(c[:, 0] - seed["cx"], c[:, 1] - seed["cy"]) <= cfg.self_hit_radius
            c, s = c[~self_hit], s[~self_hit]
            deficit, displacement = deficit[~self_hit], displacement[~self_hit]
            rank = np.argsort(s)[::-1]           # both variants ranked by score
            det = pd.DataFrame({"rank": np.arange(len(rank)), "cx": c[rank, 0],
                                "cy": c[rank, 1], "score": s[rank]})
            det_out, gt_out, metrics, _ = ev.evaluate_run(det, gt_eval, radius, area_mm2,
                                                          roi_shape=rgb.shape)
            out[name] = det_out
            rows.append({
                "file_name": fn, "tumor_type": sel["tumor_type"], "variant": name,
                "match_radius_px": round(radius, 1), "n_peaks_pre_nms": n_peaks,
                "n_after_nms": len(det_out), "K": metrics["k"],
                "recall_at_k": round(metrics["recall_at_k"], 4),
                "topk_tp": metrics["topk_tp"], "all_tp": metrics["all_tp"],
                "mitotic_recall": round(metrics["mitotic_recall"], 4),
                "sens@8fp_mm2": round(metrics["sens@8fp_mm2"], 4),
                "sens@64fp_mm2": round(metrics["sens@64fp_mm2"], 4),
                "coverage_frac": metrics["coverage_frac"],
                "median_survivor_score_deficit": round(float(np.median(deficit)), 4) if len(deficit) else np.nan,
                "frac_survivors_not_local_max": round(float((deficit > 1e-9).mean()), 4) if len(deficit) else np.nan,
                "median_displacement_px": round(float(np.median(displacement[deficit > 1e-9])), 2)
                    if (deficit > 1e-9).any() else 0.0,
                "p90_displacement_px": round(float(np.percentile(displacement[deficit > 1e-9], 90)), 2)
                    if (deficit > 1e-9).any() else 0.0,
                "mean_survivor_score": round(float(s.mean()), 4) if len(s) else np.nan,
            })

        # TP flips, both directions, at the top-K budget each variant is scored at.
        k = int(rows[-1]["K"])
        def topk_hits(d):
            t = d.head(k)
            return set(t.loc[t["bucket"] == ev.HUMAN_CORRECT_LABEL, "matched_ann_id"])
        a, b = topk_hits(out["score"]), topk_hits(out["x_then_rank_by_score"])
        flip_rows.append({"file_name": fn, "tumor_type": sel["tumor_type"], "K": k,
                          "tp_score_nms": len(a), "tp_x_nms": len(b),
                          "lost_by_x_nms": len(a - b), "gained_by_x_nms": len(b - a)})
        del rgb
        print(f"{fn}: peaks={n_peaks}  score-NMS {rows[-2]['n_after_nms']} dets "
              f"recall@K={rows[-2]['recall_at_k']:.3f} | x-NMS {rows[-1]['n_after_nms']} dets "
              f"recall@K={rows[-1]['recall_at_k']:.3f}  ({time.time()-t0:.0f}s)", flush=True)

    df = pd.DataFrame(rows).sort_values(["tumor_type", "variant"]).reset_index(drop=True)
    flips = pd.DataFrame(flip_rows)
    df.to_csv("results/fs_nms_ordering_metrics.csv", index=False)
    flips.to_csv("results/fs_nms_ordering_tp_flips.csv", index=False)
    pd.set_option("display.width", 240); pd.set_option("display.max_columns", 40)
    print("\n", df.to_string(index=False))
    print("\n", flips.to_string(index=False))
    print("\nmean recall@K by variant:")
    print(df.groupby("variant")[["recall_at_k", "mitotic_recall", "n_after_nms",
                                 "mean_survivor_score"]].mean().round(4).to_string())


if __name__ == "__main__":
    main()
