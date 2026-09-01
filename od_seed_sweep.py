"""Five seeds per ROI: does the chromatin-density re-rank survive seed variance?

Every number in `Research Logs/2026-08-31-chromatin-density-rerank.md` comes from one
seed per image. The 2026-08-27 corrected re-run showed seed choice alone moved 246.tiff's
recall@K from 0.191 to 0.043 -- larger than several effects that were being treated as
real -- so no single-seed result here should be quoted until this has run.

For each of the 7 domain ROIs and each of 5 seeds:

* pick the seed annotation from a ``(seed_index, image_id)`` stream, the same
  construction `experiment.run_experiment` uses;
* match once at a deep floor, NMS, drop the self-hit;
* record the fixed-0.5 arm and the per-image ``z >= 2.5`` arm;
* on the z arm, compute the human-workload curve under both ranking keys.

The image is loaded and colour-converted once per ROI and reused across its 5 seeds --
`rgb2hed` over a 39-megapixel ROI dominates everything else here.

Writes results/od_seed_sweep.csv.
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

EXTRACT_FLOOR = 0.25
Z = 2.5
N_SEEDS = 5
CHANNEL = "rgb"


def workload(det: pd.DataFrame, n_mit: int) -> dict:
    """Candidates a reader must work through to reach 50/80/100% of the mitoses."""
    tp = (det["bucket"] == ev.HUMAN_CORRECT_LABEL).to_numpy().cumsum()
    out = {}
    for frac in (0.5, 0.8, 1.0):
        i = int(np.searchsorted(tp, int(np.ceil(frac * n_mit))))
        out[f"read_{int(frac * 100)}"] = i + 1 if i < len(tp) else np.nan
    return out


def main():
    images, ann = ds.load_annotations()
    rows = []
    for _, row in ex.select_domain_images(images, ann).iterrows():
        fn = row["file_name"]
        path = f"images/{fn}"
        t0 = time.time()
        rgb = ds.load_roi(path)
        radius = ev.radius_px(ds.roi_mpp(path))
        area = ds.check_roi_scale(path, rgb.shape)
        gray = ch.to_gray_inverted(rgb)
        hem = cm.hematoxylin_od(rgb)          # unclipped OD; see chromatin.hematoxylin_od
        img = ch.to_channel(rgb, CHANNEL)     # the search channel, computed once
        gt = ds.image_annotations(ann, fn)
        gt_mit = gt[gt["category_id"] == ds.MITOTIC]
        print(f"{fn}: loaded [{time.time() - t0:.0f}s]", flush=True)

        for s in range(N_SEEDS):
            t1 = time.time()
            rng = np.random.default_rng([s, int(row["image_id"])])
            seed, _ = ss.pick_seed(gt_mit, gray, rng, tm.PATCH_SIZE // 2, rgb.shape)
            gt_eval = gt[gt["ann_id"] != seed["ann_id"]].reset_index(drop=True)
            n_mit = int((gt_eval["category_id"] == ds.MITOTIC).sum())
            base = ss.tightened_base_size(gray, seed["cx"], seed["cy"])
            cfg = fs.FSConfig(channel=CHANNEL, base_size=base, score_threshold=EXTRACT_FLOOR)

            patch = tm.read_padded_patch(img, seed["cx"], seed["cy"], cfg.patch_size)
            templates, _ = tm.build_augmentations(patch, base, cfg.scales,
                                                  cfg.n_angles, cfg.flips)
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
            rec = {"file_name": fn, "tumor_type": row["tumor_type"], "seed_index": s,
                   "seed_ann_id": int(seed["ann_id"]), "base_size": base,
                   "n_gt_mitotic": n_mit, "map_median": round(med, 4),
                   "mad_scale": round(mad, 4), "z_of_0.5": round((0.5 - med) / mad, 2),
                   "z2.5_cut": round(cut, 4),
                   # True when z=2.5 falls below the extraction floor, i.e. the arm is
                   # limited by the floor rather than by its own z level.
                   "floor_limited": bool(cut < EXTRACT_FLOOR)}

            # --- arm A: the current pipeline, fixed 0.5 in correlation units ---
            a = pool[pool.score >= 0.5].sort_values("score", ascending=False).reset_index(drop=True)
            det_a, _, m_a, _ = ev.evaluate_run(a, gt_eval, radius, area)
            rec.update(n_det_fixed05=len(a), recall_fixed05=m_a["mitotic_recall"],
                       recall_at_k_fixed05=m_a["recall_at_k"])
            rec.update({f"fixed05_{k}": v for k, v in workload(det_a, n_mit).items()})

            # --- arm B/C: per-image z floor, ranked by score then by chromatin ---
            sub = cm.score_detections(pool[pool.score >= cut], hem)
            rec["n_det_z2.5"] = len(sub)
            for name, d in (("z25_score", sub.sort_values("score", ascending=False)),
                            ("z25_od", cm.rerank(sub))):
                det, _, m, _ = ev.evaluate_run(d.reset_index(drop=True), gt_eval, radius, area)
                rec[f"recall_{name}"] = m["mitotic_recall"]
                rec[f"recall_at_k_{name}"] = m["recall_at_k"]
                rec.update({f"{name}_{k}": v for k, v in workload(det, n_mit).items()})

            rows.append(rec)
            print(f"  seed {s} ann={seed['ann_id']} base={base}px "
                  f"z(0.5)={rec['z_of_0.5']:.2f} n@0.5={rec['n_det_fixed05']} "
                  f"n@z2.5={rec['n_det_z2.5']} read50 score={rec['z25_score_read_50']} "
                  f"od={rec['z25_od_read_50']} [{time.time() - t1:.0f}s]", flush=True)

    pd.DataFrame(rows).to_csv("results/od_seed_sweep.csv", index=False)
    print("\nwrote results/od_seed_sweep.csv")


if __name__ == "__main__":
    main()
