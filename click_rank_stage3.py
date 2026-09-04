"""Step 3: score every ranker on one identical candidate set, per ROI and per seed.

The design rules this run has to obey, all of them from mistakes this project already made
once (`Research Logs/2026-08-31-premise-test-results.md`, and the §7 protocol of
`2026-09-01-one-click-retrieval-literature.md`):

* **One candidate set per (ROI, seed); arms differ only in the sort key.** So a difference
  between arms is a statement about the ranker and nothing else.
* **`gt_eval` drops the seed's own annotation**, and the same `gt_eval` goes to every arm --
  including the seedless ones, which do not know a seed exists.
* **The clicked cell is dropped from every arm** -- candidates within `self_hit_radius` of
  the click, plus the single nearest candidate to it whatever its distance -- so it cannot
  be handed back to the reader as a discovery, and cannot occupy rank 1 at cosine 1.0.
  `find_and_suppress` uses the same tight radius, and for the same reason: a match-radius
  exclusion would also delete a genuine neighbouring object (annotations get as close as
  26.2 px).
* **Two query constructions per encoder.** ``click``: embed a crop centred on the clicked
  pixel. ``snap``: use the embedding of the candidate the click lands on. A blob centroid
  sits a few pixels off the annotated centre, so a click-centred query is framed
  differently from every candidate it is compared against; snapping removes that
  systematic mismatch, and is what a real tool would do anyway (click -> nearest detected
  object).
* **The decision metric is `read_50`, not AUC.** With ~110 look-alikes against ~20,000
  unannotated candidates, reading depth is dominated by the unannotated contrast; an
  AUC-vs-look-alike gate would reject a ranker that wins on the metric a reader pays.
  The two AUCs are reported as diagnostics that explain *why* read_50 moved.
* **AUC labels come from one canonical bucketing** (the candidate set in blob order), not
  from each arm's own ranking, because `bucket_detections` matches greedily in row order
  and would otherwise label the same candidate differently for different arms.
* **Every encoder is scored twice**: plain cosine, and SimpleShot's CL2N -- subtract the
  ROI's mean embedding, re-normalise, then take the cosine. Centering is not a detail in
  the nearest-neighbour literature: without it every embedding shares a large mean
  direction and the cosine is dominated by it. Reporting only the uncentred number would
  understate what the click can do.
* **Report every seed, and gate on the worst one.** The premise test found read_50 varying
  2.2x with which cell was clicked.

Two passes are run. ``1click`` is the product interaction as specified. ``3click`` averages
three seeds' embeddings into a SimpleShot-style prototype and drops all three annotations
from ``gt_eval``, so "one click is simply too little signal" cannot explain a negative
result. Arms are comparable *within* a pass and not across them -- the denominators differ
by the number of seed annotations removed.

Writes `results/click_rank_metrics.csv` (one row per ROI x seed x arm x budget) and
`results/click_rank_auc.csv` (diagnostics).
"""

from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

from midog_utils import compare as cp
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from midog_utils import find_and_suppress as fs

ROIS = ("301.tiff", "246.tiff")
OUT = os.environ.get("CLICK_RANK_DIR", os.path.expanduser("~/.cache/annotatedx_click_rank"))
SELF_HIT_RADIUS = fs.FSConfig().self_hit_radius   # 5.0 px


def auc_and_ci(pos: np.ndarray, neg: np.ndarray):
    """Mann-Whitney AUC with a Hanley-McNeil standard error and normal 95% CI.

    Identical to `lookalike_auc_unbiased.auc_and_ci`; duplicated rather than imported so
    this script has no dependency on a top-level module that may be edited independently.
    """
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    n1, n2 = len(pos), len(neg)
    if not n1 or not n2:
        return float("nan"), float("nan"), float("nan"), n1, n2
    gt_ = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    a = (gt_ + 0.5 * eq) / (n1 * n2)
    q1 = a / (2 - a)
    q2 = 2 * a * a / (1 + a)
    se = np.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n2 - 1) * (q2 - a * a)) / (n1 * n2))
    return float(a), float(a - 1.96 * se), float(a + 1.96 * se), n1, n2


def raw_pixel_embeddings(crops: np.ndarray, side: int = 32, frac: float = 0.5) -> np.ndarray:
    """`TM_CCOEFF_NORMED` as an embedding: mean-centred, L2-normalised grey-inverted crop.

    Cosine between two of these *is* Pearson correlation, which is exactly what
    `cv2.matchTemplate(..., TM_CCOEFF_NORMED)` computes. This arm is therefore the
    incumbent representation -- the one the premise test measured the click through -- put
    on the same footing as the learned encoders. If a learned encoder beats chromatin and
    this does not, the gain is attributable to the representation rather than to the click.
    """
    import cv2
    p = crops.shape[1]
    k = int(round(p * frac))
    o = (p - k) // 2
    sub = crops[:, o:o + k, o:o + k, :]
    out = np.empty((len(sub), side * side), np.float32)
    for i, c in enumerate(sub):
        g = 255.0 - cv2.cvtColor(c, cv2.COLOR_RGB2GRAY).astype(np.float32)  # gray-inverted
        v = cv2.resize(g, (side, side), interpolation=cv2.INTER_AREA).ravel()
        v -= v.mean()
        n = np.linalg.norm(v)
        out[i] = v / n if n > 0 else v
    return out


def main():
    metrics, aucs, checks = [], [], []
    for fn in ROIS:
        meta = json.load(open(f"{OUT}/meta_{fn}.json"))
        cands = pd.read_csv(f"{OUT}/cands_{fn}.csv")
        gt = pd.read_csv(f"{OUT}/gt_{fn}.csv")
        seeds = pd.read_csv(f"{OUT}/seeds_{fn}.csv")
        radius, mpp = meta["radius"], meta["mpp"]
        roi_shape = tuple(meta["roi_shape"])

        # encoders present on disk, plus the raw-pixel control computed here
        embs = {}
        for path in sorted(glob.glob(f"{OUT}/emb_{fn}_*.npy")):
            name = os.path.basename(path)[len(f"emb_{fn}_"):-4]
            e = np.load(path)
            se = np.load(f"{OUT}/embseed_{fn}_{name}.npy")
            if len(e) != len(cands):
                raise ValueError(f"{name}: {len(e)} embeddings for {len(cands)} candidates")
            embs[name] = (e, se)
        crops = np.load(f"{OUT}/crops_{fn}.npy", mmap_mode="r")
        seedcrops = np.load(f"{OUT}/seedcrops_{fn}.npy")
        embs["rawpixel_ccoeff"] = (raw_pixel_embeddings(np.asarray(crops)),
                                   raw_pixel_embeddings(seedcrops))
        # SimpleShot CL2N: centre on the ROI's mean embedding, then re-normalise. The mean
        # is over all candidates in this ROI and so does not depend on which cell was
        # clicked -- it is ingest-time information, available before any click.
        cl2n = {}
        for name, (e, se) in embs.items():
            mu = e.mean(axis=0, keepdims=True)
            def _c(x):
                x = x - mu
                n = np.linalg.norm(x, axis=1, keepdims=True)
                return (x / np.where(n > 0, n, 1.0)).astype(np.float32)
            cl2n[name] = (_c(e), _c(se))
        print(f"[{fn}] {len(cands)} candidates, encoders: {sorted(embs)}", flush=True)

        n_seeds = len(seeds)
        passes = [("1click", [[s] for s in range(n_seeds)]),
                  ("3click", [[s, (s + 1) % n_seeds, (s + 2) % n_seeds]
                              for s in range(n_seeds)])]
        for pass_name, seed_sets in passes:
          for s, idxs in enumerate(seed_sets):
            srows = seeds.iloc[idxs]
            seed_ann = sorted({int(a) for a in srows["ann_id"]})
            gt_eval = gt[~gt["ann_id"].isin(seed_ann)].reset_index(drop=True)

            dists = np.stack([np.hypot(cands["cx"].to_numpy() - r["cx"],
                                       cands["cy"].to_numpy() - r["cy"])
                              for _, r in srows.iterrows()])
            d_seed = dists.min(axis=0)
            snap_idx = dists.argmin(axis=1)          # candidate each click lands on
            keep = d_seed > SELF_HIT_RADIUS
            keep[snap_idx] = False
            sub = cands[keep].reset_index(drop=True)

            # one canonical bucketing for the AUC diagnostics (see module docstring)
            canon, _ = ev.bucket_detections(sub, gt_eval, radius)
            b = canon["bucket"].to_numpy()
            tp_m = b == ev.HUMAN_CORRECT_LABEL
            lk_m = b == ev.HUMAN_REJECTED_LABEL
            un_m = b == ev.NON_HUMAN_FINDINGS

            arms = [
                cp.Arm("chromatin", (lambda d=sub: d), rank_key="od",
                       coverage_key=f"{fn}:{s}:cands",
                       extra={"encoder": "chromatin_density_51px"}),
                cp.Arm("blob_native", (lambda d=sub: d), rank_key="score",
                       coverage_key=f"{fn}:{s}:cands",
                       extra={"encoder": "otsu_component_mean"}),
            ]
            stats = {"chromatin": sub["od"].to_numpy(), "blob_native": sub["score"].to_numpy()}
            for name, (e, se) in embs.items():
                ec, sec = cl2n[name]
                for variant, ev_, qsrc in (("", e, se[idxs]),
                                           ("_cl2n", ec, sec[idxs]),
                                           ("_cl2n_snap", ec, ec[snap_idx])):
                    q = np.asarray(qsrc).mean(axis=0)
                    q = q / (np.linalg.norm(q) or 1.0)   # prototype; identity for k=1
                    cos = (ev_[keep] @ q).astype(np.float64)
                    key = f"cos_{name}{variant}"
                    frame = sub.assign(**{key: cos})
                    arms.append(cp.Arm(key, (lambda d=frame: d), rank_key=key, seeded=True,
                                       coverage_key=f"{fn}:{s}:cands",
                                       extra={"encoder": name,
                                              "variant": variant.strip("_") or "cosine",
                                              "query": "snap" if "snap" in variant else "click"}))
                    stats[key] = cos

            ctx = {"file_name": fn, "image_id": meta["image_id"], "mpp": round(mpp, 4),
                   "pass": pass_name, "n_clicks": len(idxs),
                   "seed_index": s, "seed_ann_id": "|".join(map(str, seed_ann)),
                   "n_candidates_after_selfhit": int(len(sub)),
                   "n_selfhit_dropped": int((~keep).sum()),
                   "dist_seed_to_nearest_cand": round(float(dists.min()), 2),
                   "crop_um": meta["crop_um"]}
            got = cp.evaluate_arms(arms, gt_eval, radius, roi_shape=roi_shape, mpp=mpp,
                                   context=ctx, checks=checks)
            metrics.append(got)

            for stat_name, v in stats.items():
                for contrast, pos, neg in (("tp_vs_lookalike", v[tp_m], v[lk_m]),
                                           ("tp_vs_unannotated", v[tp_m], v[un_m])):
                    a, lo, hi, n1, n2 = auc_and_ci(pos, neg)
                    aucs.append(dict(file_name=fn, seed_index=s, pass_name=pass_name,
                                     n_clicks=len(idxs), statistic=stat_name,
                                     contrast=contrast, auc=round(a, 4),
                                     ci_lo=round(lo, 4), ci_hi=round(hi, 4),
                                     n_pos=n1, n_neg=n2))
            print(f"  [{pass_name}] seed {s}: ann {seed_ann}, {len(sub)} candidates "
                  f"({(~keep).sum()} self-hit dropped), "
                  f"{tp_m.sum()} mitotic / {lk_m.sum()} look-alike captured", flush=True)

    m = pd.concat(metrics, ignore_index=True)
    a = pd.DataFrame(aucs)
    os.makedirs("results", exist_ok=True)
    m.to_csv("results/click_rank_metrics.csv", index=False)
    a.to_csv("results/click_rank_auc.csv", index=False)
    bad = [c for c in checks if not c.get("passed", True)]
    print(f"\ninvariant checks: {len(checks)} run, {len(bad)} failed")

    pd.set_option("display.width", 220)
    for pn in m["pass"].unique():
        d = m[(m["budget"] == 500) & (m["pass"] == pn)]
        piv = d.pivot_table(index=["file_name", "arm"], columns="seed_index",
                            values="read_50", dropna=False)
        piv["median"] = piv.median(axis=1)
        piv["worst"] = piv.max(axis=1)
        print(f"\nread_50 per seed, pass={pn} (lower is better; NaN = 50% never reached):")
        print(piv.round(0).to_string())


if __name__ == "__main__":
    sys.exit(main())
