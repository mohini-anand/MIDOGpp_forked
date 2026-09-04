"""Leave-one-*domain*-out for the cheap FP filter, with per-arm re-matching.

Two corrections to `fp_filter_probe.py`, both of which change the numbers:

**1. Whole domains are held out, not single ROIs.** The seven ROIs come in domain pairs --
301/300 mast cell, 246/245 lymphosarcoma, 201/202 lung, 405 soft tissue -- sharing tumour
type *and* scanner. Holding out 301 alone leaves 300 in training, so leave-one-ROI-out
measures "generalises to another ROI of a tumour I have already seen", not the deployment
condition. Both protocols are reported side by side so the size of the leak is visible.

**2. Each arm is re-matched in its own ranking order.** `fp_filter_probe.py` fixed the
true-positive labels from one greedy match taken in `blob_score` order and reused them for
every arm. That is not neutral: `evaluate.greedy_match` lets each annotation be claimed by
the best-ranked detection within the match radius, so an arm that reorders the list should
be re-matched in that order. Reusing blob_score's assignment penalises every other arm --
on 201.tiff it reports `chromatin_od51` needing 12,124 candidates to reach 95% where
re-matching gives 823. This re-matches per arm, matching `tail_sensitivity.py`, so the two
tables are directly comparable.

Reads results/fp_filter_features.csv. Writes results/fp_filter_domain.csv.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from fp_filter_probe import FEATURES, QUANTS


def fit_score(Xtr, ytr, Xte, kind: str) -> np.ndarray:
    if kind == "logreg":
        m = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=2000, class_weight="balanced"))
    else:
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.08,
                                           max_leaf_nodes=15, l2_regularization=1.0,
                                           class_weight="balanced", random_state=0)
    m.fit(Xtr, ytr)
    return m.predict_proba(Xte)[:, 1]


def tail_rematched(sub: pd.DataFrame, score: np.ndarray, gt: pd.DataFrame,
                   radius: float) -> list[dict]:
    """Sort by ``score``, re-run the greedy match, and read depths off `matched_rank`."""
    order = sub.assign(_s=score).sort_values("_s", ascending=False, kind="mergesort")
    order = order.reset_index(drop=True)
    _det, gt_out = ev.bucket_detections(order, gt.reset_index(drop=True), radius)
    mit = gt_out[gt_out["category_id"] == ds.MITOTIC]
    ranks = np.sort(mit.loc[mit["found"], "matched_rank"].to_numpy())
    n_gt, n_found, n_pool = len(mit), len(ranks), len(order)
    rows = []
    for q in QUANTS:
        if not n_found:
            continue
        i = int(np.ceil(q * n_found)) - 1
        rows.append(dict(quantile=q, depth=int(ranks[i]) + 1, n_found=n_found,
                         n_gt_mitotic=n_gt, ceiling=round(n_found / n_gt, 4),
                         frac_of_pool=round((int(ranks[i]) + 1) / n_pool, 4)))
    return rows


def main():
    df = pd.read_csv(os.environ.get("FP_FEATCACHE", ".cache_tail/fp_filter_features.csv"))
    X = np.nan_to_num(df[FEATURES].to_numpy(dtype=np.float64),
                      nan=0.0, posinf=0.0, neginf=0.0)
    y = df["y"].to_numpy()
    _images, ann = ds.load_annotations()

    rows = []
    for fn in df["file_name"].unique():
        te = (df["file_name"] == fn).to_numpy()
        sub = df[te].reset_index(drop=True)
        dom = sub["domain"].iloc[0]
        gt = ds.image_annotations(ann, fn)
        radius = ev.radius_px(ds.roi_mpp(f"images/{fn}"))

        arms = {"blob_score": sub["score"].to_numpy(),
                "chromatin_od51": np.nan_to_num(sub["od51"].to_numpy(), nan=-1e9)}
        for proto, tr in (("loro", ~te),
                          ("lodo", (~te) & (df["domain"] != dom).to_numpy())):
            for kind in ("logreg", "hgb"):
                arms[f"{kind}_{proto}"] = fit_score(X[tr], y[tr], X[te], kind)

        for name, s in arms.items():
            for r in tail_rematched(sub, np.asarray(s, dtype=np.float64), gt, radius):
                rows.append(dict(file_name=fn, domain=dom, ranker=name,
                                 n_pool=len(sub), **r))
        print(f"{fn} ({dom}) done", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv("results/fp_filter_domain.csv", index=False)
    pd.set_option("display.width", 260)
    piv = out.pivot_table(index=["file_name", "ranker"], columns="quantile",
                          values="depth", aggfunc="first")
    base = out.groupby(["file_name", "ranker"])[["n_pool", "n_found", "ceiling"]].first()
    print("\n=== candidates read to reach each sensitivity level (of captured) ===")
    print(base.join(piv).to_string())


if __name__ == "__main__":
    main()
