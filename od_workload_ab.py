"""The Fix-2 headline table: human workload under both ranking keys, same candidates.

Generates `results/od_workload_ab.csv`. Previously this table had no generator in the
repo -- an audit gap, since it carries the headline number.

For each of the 7 domain ROIs: match once, NMS, take the ``z >= 2.5`` candidate pool,
then rank it two ways -- by the correlation score and by chromatin density -- and ask how
far down a reader must go to reach 50/80/100% of the mitotic figures.

Also records what the ranking does to the *look-alike* bucket. Chromatin density promotes
pathologist-rejected look-alikes nearly as hard as it promotes mitoses
(AUC(od, look-alike vs unannotated) is close to AUC(od, TP vs unannotated)), so the
187:1 unannotated-to-look-alike ratio of the full list does **not** survive to the top of
a chromatin-ranked list. `lk_in_topk` is the number that shows it.
"""

from __future__ import annotations

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

EXTRACT_FLOOR, Z, CHANNEL = 0.25, 2.5, "rgb"


def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    if not len(pos) or not len(neg):
        return float("nan")
    return float((pos[:, None] > neg[None, :]).mean() + 0.5 * (pos[:, None] == neg[None, :]).mean())


def main():
    images, ann = ds.load_annotations()
    rows = []
    for _, row in ex.select_domain_images(images, ann).iterrows():
        fn = row["file_name"]
        path = f"images/{fn}"
        rgb = ds.load_roi(path)
        radius = ev.radius_px(ds.roi_mpp(path))
        area = ds.check_roi_scale(path, rgb.shape)
        gray = ch.to_gray_inverted(rgb)
        gt = ds.image_annotations(ann, fn)
        seed, _ = ss.pick_seed(gt[gt["category_id"] == ds.MITOTIC], gray,
                               np.random.default_rng([0, int(row["image_id"])]),
                               tm.PATCH_SIZE // 2, rgb.shape)
        gt_eval = gt[gt["ann_id"] != seed["ann_id"]].reset_index(drop=True)
        n_mit = int((gt_eval["category_id"] == ds.MITOTIC).sum())
        base = ss.tightened_base_size(gray, seed["cx"], seed["cy"])
        cfg = fs.FSConfig(channel=CHANNEL, base_size=base, score_threshold=EXTRACT_FLOOR)

        img = ch.to_channel(rgb, CHANNEL)
        templates, _ = tm.build_augmentations(
            tm.read_padded_patch(img, seed["cx"], seed["cy"], cfg.patch_size),
            base, cfg.scales, cfg.n_angles, cfg.flips)
        fused, _, valid = tm.fused_response(img, templates)
        samp = fused[::8, ::8]
        samp = samp[samp > -1.5]
        med = float(np.median(samp))
        mad = float(1.4826 * np.median(np.abs(samp - med)))
        centers, scores = tm.extract_peaks(fused, valid, cfg.peak_min_distance,
                                           EXTRACT_FLOOR, cfg.max_peaks)
        keep = nms_by_distance(centers, scores, radius)
        centers, scores = centers[keep], scores[keep]
        ok = np.hypot(centers[:, 0] - seed["cx"], centers[:, 1] - seed["cy"]) > cfg.self_hit_radius
        pool = pd.DataFrame({"cx": centers[ok, 0], "cy": centers[ok, 1], "score": scores[ok]})

        cut = med + Z * mad
        pool = cm.score_detections(pool[pool.score >= cut], cm.hematoxylin_od(rgb))
        for name, d in (("score", pool.sort_values("score", ascending=False, kind="mergesort")),
                        ("od", cm.rerank(pool))):
            det, _, _, _ = ev.evaluate_run(d.reset_index(drop=True), gt_eval, radius, area)
            b = det["bucket"].to_numpy()
            tp = (b == ev.HUMAN_CORRECT_LABEL).cumsum()
            r = {"image": fn, "n_mit": n_mit, "rank_by": name, "candidates": len(det),
                 "z2.5_cut": round(cut, 4), "floor_limited": bool(cut < EXTRACT_FLOOR)}
            for f in (0.5, 0.8, 1.0):
                i = int(np.searchsorted(tp, int(np.ceil(f * n_mit))))
                r[f"read_for_{int(f * 100)}pct"] = i + 1 if i < len(tp) else np.nan
            top = b[:n_mit]
            r["tp_in_topk"] = int((top == ev.HUMAN_CORRECT_LABEL).sum())
            r["lk_in_topk"] = int((top == ev.HUMAN_REJECTED_LABEL).sum())
            r["un_in_topk"] = int((top == ev.NON_HUMAN_FINDINGS).sum())
            r["n_lookalike_total"] = int((b == ev.HUMAN_REJECTED_LABEL).sum())
            r["n_unannotated_total"] = int((b == ev.NON_HUMAN_FINDINGS).sum())
            if name == "od":
                o = det["od"].to_numpy()
                r["auc_tp_vs_unann"] = round(auc(o[b == ev.HUMAN_CORRECT_LABEL],
                                                 o[b == ev.NON_HUMAN_FINDINGS]), 3)
                r["auc_lk_vs_unann"] = round(auc(o[b == ev.HUMAN_REJECTED_LABEL],
                                                 o[b == ev.NON_HUMAN_FINDINGS]), 3)
                r["auc_tp_vs_lk"] = round(auc(o[b == ev.HUMAN_CORRECT_LABEL],
                                              o[b == ev.HUMAN_REJECTED_LABEL]), 3)
            rows.append(r)
        print(f"{fn} done", flush=True)

    pd.DataFrame(rows).to_csv("results/od_workload_ab.csv", index=False)
    print("wrote results/od_workload_ab.csv")


if __name__ == "__main__":
    main()
