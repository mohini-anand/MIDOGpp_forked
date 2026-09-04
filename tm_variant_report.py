"""Turn a `tm_variant_sweep.py` grid into the tables section 5 of the log needs.

Kept separate from the sweep so re-reading a result never risks re-running it, and so every
table in the log is reproducible from the committed CSV by one command.

The reading order is fixed here rather than left to the reader, because three of these numbers
are only interpretable next to another one: `ceiling` next to `coverage_frac` and a matched
pool size, `read_100` next to `n_unreached`, and any median next to its `cv`.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

RESULTS = "results"


def load(stem):
    df = pd.read_csv(f"{RESULTS}/{stem}.csv")
    var = pd.read_csv(f"{RESULTS}/{stem}_seed_variance.csv")
    return df, var


def tail_table(var: pd.DataFrame, metric="read_100", decision_only=True) -> pd.DataFrame:
    """Worst-of-N depth per (ROI, arm), with the two columns that qualify it.

    `worst` is NaN whenever any seed failed to reach the quantile -- deliberately, since a
    worst-of-5 computed over the 3 seeds that happened to succeed is not a worst-of-5. Read
    `n_unreached` first; a high `cv` means the arm's advantage is not something a pathologist
    with one click can rely on.
    """
    v = var[var["metric"] == metric]
    if decision_only:
        v = v[v["decision_grade"] & (v["n_seeds"] >= 5)]
    # `pivot`, not `pivot_table`: the latter builds the full tumour_type x file_name
    # cross-product and fills it with NaN, which reads as "this arm failed here" for 60
    # ROI/domain pairs that do not exist.
    return v.pivot(index=["tumor_type", "file_name"], columns="arm", values="worst").round(0)


def spread_table(var: pd.DataFrame, metric="read_100") -> pd.DataFrame:
    """Per-arm summary over the decision-grade ROIs, reach reported before depth.

    ``cells_reached`` comes first because ``read_100`` is frequently unreachable: it is NaN
    whenever that arm's candidate pool never contained the last mitotic figure at all, which
    is a *ceiling* failure and not a deep read. Ranking arms by depth without showing how
    often the depth exists would flatter whichever arm gives up earliest.

    ``cv_median`` is the median across ROIs of the within-ROI coefficient of variation over
    seeds -- the number that says whether an arm's advantage is something a pathologist with
    one click can actually rely on.
    """
    v = var[(var["metric"] == metric) & var["decision_grade"] & (var["n_seeds"] >= 5)]
    g = v.groupby("arm")
    out = pd.DataFrame({
        "n_rois": g.size(),
        "cells": g["n_seeds"].sum(),
        "cells_reached": g["n_seeds"].sum() - g["n_unreached"].sum(),
        "median_of_roi_medians": g["median"].median().round(1),
        "cv_median": g["cv"].median().round(3),
        "cv_max": g["cv"].max().round(3),
        "best_roi": g["min"].min().round(0),
        "worst_roi": g["max"].max().round(0),
    })
    return out.sort_values("median_of_roi_medians")


def paired_summary(stem, a, b, decision_only=True) -> pd.DataFrame:
    """Sign counts and Wilcoxon over the per-seed deltas, which is where exp1's claim lives.

    Both arms saw the identical click on the identical ROI, so seed variance -- the dominant
    noise term in every prior experiment here -- cancels exactly. A negative delta means arm
    ``a`` reads *fewer* candidates, i.e. is better.
    """
    from scipy.stats import wilcoxon
    d = pd.read_csv(f"{RESULTS}/{stem}_paired_{a}_vs_{b}.csv")
    if decision_only:
        d = d[d["decision_grade"]]
    rows = []
    for metric in ["read_50_delta", "read_95_delta", "read_100_delta",
                   "full_list_recall_delta", "recall_at_250_delta"]:
        v = d[metric].dropna()
        better = "greater" if "recall" in metric else "less"
        wins = int((v > 0).sum() if "recall" in metric else (v < 0).sum())
        p = float(wilcoxon(v).pvalue) if len(v) > 1 and v.abs().sum() > 0 else np.nan
        rows.append({"metric": metric, "n_pairs": len(v), "n_missing": int(d[metric].isna().sum()),
                     f"{a}_better": wins, "ties": int((v == 0).sum()),
                     "median_delta": round(float(v.median()), 4),
                     "wilcoxon_p": round(p, 5) if np.isfinite(p) else np.nan,
                     "direction_of_better": better})
    return pd.DataFrame(rows)


def click_worth(df: pd.DataFrame, decision_only=True) -> pd.DataFrame:
    """How much of each method's performance the click actually causes.

    The ratio of the click-seeded arm's depth to the click-free disc template's, on the same
    ROI and the same method. A ratio near 1.0 means that variant is a blob detector wearing a
    template, and its wins do not belong to the pathologist -- which is the outcome section 6
    of the log pre-registers as a real possible answer.
    """
    b = df[df["budget"] == df["budget"].min()].copy()
    if decision_only:
        b = b[b["decision_grade"]]
    b = b[b["pool_scope"] == "full"]
    seeded = b[b["arm"].str.startswith("tm_") & ~b["arm"].str.endswith("_od")]
    disc = b[b["arm"].str.startswith("disc_") & ~b["arm"].str.endswith("_od")]
    s = seeded.groupby(["file_name", "method"])[["read_50", "read_95", "full_list_recall"]].median()
    d = disc.groupby(["file_name", "method"])[["read_50", "read_95", "full_list_recall"]].median()
    j = s.join(d, lsuffix="_click", rsuffix="_disc", how="inner")
    j["read_50_ratio"] = (j["read_50_disc"] / j["read_50_click"]).round(2)
    j["read_95_ratio"] = (j["read_95_disc"] / j["read_95_click"]).round(2)
    return j.reset_index().sort_values(["method", "file_name"])


def full_ceiling(df: pd.DataFrame, decision_only=True) -> pd.DataFrame:
    """Goal (1) on each arm's own untruncated pool, median over seeds.

    Read next to `ceiling_table`, never instead of it: this number is inflated by however
    much of the ROI the pool tiles, which `coverage_frac` in that table quantifies.
    """
    b = df[(df["budget"] == df["budget"].min()) & (df["pool_scope"] == "full")].copy()
    if decision_only:
        b = b[b["decision_grade"]]
    return (b.pivot_table(index="file_name", columns="arm", values="full_list_recall",
                          aggfunc="median").round(4))


def ceiling_table(df: pd.DataFrame, decision_only=True) -> pd.DataFrame:
    """Goal (1), read only where it is evidence: at matched pool size, with coverage beside it.

    An un-truncated pool tiles the ROI and reaches 1.0 by geometry -- the claim this repo
    retracted once already. `coverage_frac` is carried because equal `n_detections` is still
    not equal coverage: blob centroids cluster, so the same list length covers less ROI.
    """
    b = df[(df["budget"] == df["budget"].min()) & (df["pool_scope"] == "matched")].copy()
    if decision_only:
        b = b[b["decision_grade"]]
    b = b[~b["arm"].str.startswith(("disc_", "lookalike_"))]
    return (b.groupby(["tumor_type", "file_name", "arm"])
            .agg(ceiling=("full_list_recall", "median"),
                 n=("n_detections", "median"),
                 coverage=("coverage_frac", "median"),
                 geo_cap=("geometric_ceiling", "first"))
            .round(4).reset_index())


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("stem")
    ap.add_argument("--pair", nargs=2, metavar="ARM")
    ap.add_argument("--all-rois", action="store_true",
                    help="include ROIs below the decision-grade bar")
    args = ap.parse_args(argv)
    dec = not args.all_rois
    df, var = load(args.stem)

    for metric in ("read_50", "read_95", "read_100"):
        print(f"### Seed spread of {metric} (decision-grade ROIs, 5 delivered seeds)\n")
        print(spread_table(var, metric).to_string())
        print()
    print("\n### Worst-of-5 read_100 per ROI\n")
    print(tail_table(var, "read_100", dec).to_string())
    print("\n### Worst-of-5 read_50 per ROI\n")
    print(tail_table(var, "read_50", dec).to_string())
    print("\n### Ceiling on the FULL pool (why read_100 is often unreachable)\n")
    print(full_ceiling(df, dec).to_string())
    print("\n### Ceiling at matched pool size\n")
    print(ceiling_table(df, dec).to_string(index=False))
    cw = click_worth(df, dec)
    if len(cw):
        print("\n### What the click is worth (disc / click depth; 1.0 = click is worthless)\n")
        print(cw.to_string(index=False))
    if args.pair:
        print(f"\n### Paired per-seed delta: {args.pair[0]} minus {args.pair[1]}\n")
        print(paired_summary(args.stem, *args.pair, decision_only=dec).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
