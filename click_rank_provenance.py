"""Is it the click, or is it the representation?

`click_rank_stage3.py` found that ranking by cosine to the click in the MIDOG++-trained
FCOS backbone's feature space beats every seedless arm by 2-3x on reading depth, while the
same ranking in frozen self-supervised pathology encoders is 20-50x worse. That leaves one
question unanswered, and it is the one the product depends on:

**a mitosis-trained backbone has a mitosis-selective feature space, so cosine to *any*
mitotic-looking cell may rank mitoses high.** If a click on a pathologist-rejected
look-alike, or on a random nucleus, ranks nearly as well as a click on a genuine mitotic
figure, then the click is not carrying the information -- the encoder is -- and the same
list could be produced with no click at all by a classifier head.

This is the premise test's seed-provenance control, re-run in embedding space. The
evaluation ground truth is held identical across provenances (it always drops the *mitotic*
seed's annotation), the candidate set is identical, and only the query embedding changes.

Writes `results/click_rank_provenance.csv`.
"""

from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

from click_rank_stage3 import OUT, ROIS, SELF_HIT_RADIUS, auc_and_ci, raw_pixel_embeddings
from midog_utils import compare as cp
from midog_utils import evaluate as ev


def main():
    rows, aucs, checks = [], [], []
    for fn in ROIS:
        meta = json.load(open(f"{OUT}/meta_{fn}.json"))
        cands = pd.read_csv(f"{OUT}/cands_{fn}.csv")
        gt = pd.read_csv(f"{OUT}/gt_{fn}.csv")
        seeds = pd.read_csv(f"{OUT}/seeds_{fn}.csv")
        radius, mpp, image_id = meta["radius"], meta["mpp"], meta["image_id"]

        embs = {}
        for path in sorted(glob.glob(f"{OUT}/emb_{fn}_*.npy")):
            name = os.path.basename(path)[len(f"emb_{fn}_"):-4]
            embs[name] = np.load(path)
        crops = np.load(f"{OUT}/crops_{fn}.npy", mmap_mode="r")
        embs["rawpixel_ccoeff"] = raw_pixel_embeddings(np.asarray(crops))
        cl2n = {}
        for name, e in embs.items():
            x = e - e.mean(axis=0, keepdims=True)
            n = np.linalg.norm(x, axis=1, keepdims=True)
            cl2n[name] = (x / np.where(n > 0, n, 1.0)).astype(np.float32)

        # one canonical bucketing of the whole candidate set, to name each candidate's class
        canon, _ = ev.bucket_detections(cands, gt, radius)
        b = canon["bucket"].to_numpy()
        lk_pool = np.flatnonzero(b == ev.HUMAN_REJECTED_LABEL)
        un_pool = np.flatnonzero(b == ev.NON_HUMAN_FINDINGS)
        print(f"[{fn}] {len(cands)} candidates; look-alike pool {len(lk_pool)}, "
              f"unannotated pool {len(un_pool)}", flush=True)

        for _, srow in seeds.iterrows():
            s = int(srow["seed_index"])
            seed_ann = int(srow["ann_id"])
            gt_eval = gt[gt["ann_id"] != seed_ann].reset_index(drop=True)

            d = np.hypot(cands["cx"].to_numpy() - srow["cx"], cands["cy"].to_numpy() - srow["cy"])
            mit_idx = int(d.argmin())
            lk_idx = int(np.random.default_rng([s, image_id, 2]).choice(lk_pool))
            un_idx = int(np.random.default_rng([s, image_id, 3]).choice(un_pool))

            keep = d > SELF_HIT_RADIUS
            keep[[mit_idx, lk_idx, un_idx]] = False     # identical candidate set for all three
            sub = cands[keep].reset_index(drop=True)

            canon_s, _ = ev.bucket_detections(sub, gt_eval, radius)
            bs = canon_s["bucket"].to_numpy()
            tp_m, lk_m, un_m = (bs == ev.HUMAN_CORRECT_LABEL, bs == ev.HUMAN_REJECTED_LABEL,
                                bs == ev.NON_HUMAN_FINDINGS)

            arms, stats = [], {}
            for prov, qi in (("mitotic", mit_idx), ("lookalike", lk_idx), ("random", un_idx)):
                for name, e in cl2n.items():
                    cos = (e[keep] @ e[qi]).astype(np.float64)
                    key = f"{name}__{prov}"
                    frame = sub.assign(**{key: cos})
                    arms.append(cp.Arm(key, (lambda dd=frame: dd), rank_key=key, seeded=True,
                                       coverage_key=f"{fn}:{s}:cands",
                                       extra={"encoder": name, "provenance": prov}))
                    stats[key] = cos
            arms.append(cp.Arm("blob_native", (lambda dd=sub: dd), rank_key="score",
                               coverage_key=f"{fn}:{s}:cands",
                               extra={"encoder": "otsu_component_mean", "provenance": "none"}))
            arms.append(cp.Arm("chromatin", (lambda dd=sub: dd), rank_key="od",
                               coverage_key=f"{fn}:{s}:cands",
                               extra={"encoder": "chromatin_density_51px", "provenance": "none"}))

            ctx = {"file_name": fn, "image_id": image_id, "seed_index": s,
                   "seed_ann_id": seed_ann, "query_variant": "cl2n_snap",
                   "n_candidates": int(len(sub))}
            rows.append(cp.evaluate_arms(arms, gt_eval, radius, roi_shape=tuple(meta["roi_shape"]),
                                         mpp=mpp, context=ctx, checks=checks))
            for key, v in stats.items():
                for contrast, pos, neg in (("tp_vs_lookalike", v[tp_m], v[lk_m]),
                                           ("tp_vs_unannotated", v[tp_m], v[un_m])):
                    a, lo, hi, n1, n2 = auc_and_ci(pos, neg)
                    aucs.append(dict(file_name=fn, seed_index=s, statistic=key,
                                     contrast=contrast, auc=round(a, 4), ci_lo=round(lo, 4),
                                     ci_hi=round(hi, 4), n_pos=n1, n_neg=n2))
            print(f"  seed {s}: queries mitotic={mit_idx} lookalike={lk_idx} random={un_idx}",
                  flush=True)

    m = pd.concat(rows, ignore_index=True)
    m.to_csv("results/click_rank_provenance.csv", index=False)
    pd.DataFrame(aucs).to_csv("results/click_rank_provenance_auc.csv", index=False)
    bad = [c for c in checks if not c.get("passed", True)]
    print(f"\ninvariant checks: {len(checks)} run, {len(bad)} failed")

    d = m[m["budget"] == 500]
    piv = d.pivot_table(index=["file_name", "arm"], columns="seed_index", values="read_50")
    piv["median"] = piv.median(axis=1)
    piv["worst"] = piv.max(axis=1)
    pd.set_option("display.width", 220)
    print("\nread_50 by seed provenance (lower is better):")
    print(piv.round(0).to_string())


if __name__ == "__main__":
    sys.exit(main())
