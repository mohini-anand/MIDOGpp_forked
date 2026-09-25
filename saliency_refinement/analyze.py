"""
    Statistics for out/precision_per_run.csv: pooled precision per arm, paired contrasts with ROI as the unit (per-ROI
    click-summed TP deltas, exact sign-flip p, 95% bootstrap CI over ROIs, Holm within each family, win/loss/tie by ROI),
    template sizes, the seed's own figure leaking into the top 30, and the saliency-vs-production contrast split by
    whether both pick the same object. Writes out/contrasts.csv and out/arm_summary.csv and prints both.
"""

import numpy as np
import pandas as pd

import common as C

N_BOOT = 10000
FAMILIES = {
    "midog_box": [("prod_51", "default_51"), ("sal_51", "prod_51"), ("sal1_51", "prod_51"), ("sal_51", "default_51"), ("sal1_51", "default_51")],
    "drawn_box": [("sal_box", "prod_box"), ("sal1_box", "prod_box")],
}


def sign_flip_p(d: np.ndarray) -> tuple[float, int]:
    """
        Exact two-sided sign-flip p on integer per-ROI deltas, by convolution over the nonzero ones.

        d (np.ndarray): integer deltas, one per ROI.

        Returns tuple[float, int]: p and the number of nonzero deltas.
    """
    nz = np.abs(d[d != 0]).astype(np.int64)
    if len(nz) == 0:
        return 1.0, 0
    total = int(nz.sum())
    counts = np.zeros(2 * total + 1, dtype=np.float64)
    counts[total] = 1.0
    for v in nz:
        counts = np.roll(counts, v) + np.roll(counts, -v)  # support never wraps: |sum| <= total
    sums = np.arange(-total, total + 1)
    return float(counts[np.abs(sums) >= abs(int(d.sum()))].sum() / 2.0 ** len(nz)), len(nz)


def holm(p) -> np.ndarray:
    """
        Holm step-down adjusted p-values.

        p (array-like): raw p-values of one family.

        Returns np.ndarray: adjusted p-values in input order.
    """
    p = np.asarray(p, float)
    order = np.argsort(p, kind="stable")
    adj = np.maximum.accumulate((len(p) - np.arange(len(p))) * p[order]).clip(max=1.0)
    out = np.empty_like(adj)
    out[order] = adj
    return out


def contrast(runs: pd.DataFrame, a: str, b: str, k: int, rng) -> dict:
    """
        One paired contrast at one budget, clicks paired on (file_name, seed_ann_id, ref, pad).

        runs (pd.DataFrame): per_run rows of one rank key and click set, both arms present.
        a (str): arm.
        b (str): baseline arm.
        k (int): budget.
        rng (np.random.Generator): bootstrap stream.

        Returns dict: n_roi, n_clicks, pooled precision of each arm, pooled delta, CI, sign-flip p, W/L/T.
    """
    keys = ["file_name", "seed_ann_id", "ref", "pad"]
    ra = runs[runs.arm == a].set_index(keys)[f"tp_at_{k}"]
    rb = runs[runs.arm == b].set_index(keys)[f"tp_at_{k}"]
    if a.endswith("_51") and b.endswith("_51"):
        ra, rb = ra.droplevel(["ref", "pad"]), rb.droplevel(["ref", "pad"])
    both = pd.concat([ra.rename("a"), rb.rename("b")], axis=1, join="inner").reset_index()
    per_roi = both.groupby("file_name").agg(a=("a", "sum"), b=("b", "sum"), n=("a", "size"))
    d = (per_roi.a - per_roi.b).to_numpy()
    n = per_roi.n.to_numpy()
    p, nnz = sign_flip_p(d)
    idx = rng.integers(len(d), size=(N_BOOT, len(d)))
    boot = d[idx].sum(axis=1) / (n[idx].sum(axis=1) * k)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return dict(pair=f"{a} - {b}", K=k, n_roi=len(d), n_clicks=int(n.sum()), prec_a=per_roi.a.sum() / (n.sum() * k), prec_b=per_roi.b.sum() / (n.sum() * k), delta=d.sum() / (n.sum() * k), ci_lo=lo, ci_hi=hi, p=p, n_nonzero=nnz, wins=int((d > 0).sum()), losses=int((d < 0).sum()), ties=int((d == 0).sum()))


def main():
    """
        Print and write every table.

        Returns None.
    """
    runs = pd.read_csv(C.OUT_DIR / "precision_per_run.csv")
    runs["ref"] = runs["ref"].fillna("")
    runs["pad"] = runs["pad"].fillna(-1.0)
    geo = pd.read_csv(C.OUT_DIR / "geometry_midog_box.csv")
    rng = np.random.default_rng(20260924)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 400)

    # arm summary
    summ = runs.groupby(["click_set", "rank_key", "arm", "ref", "pad"]).agg(n=("tp_at_10", "size"), p10=("tp_at_10", "mean"), p20=("tp_at_20", "mean"), p30=("tp_at_30", "mean"), base_med=("base_size", "median"), base_q10=("base_size", lambda s: s.quantile(0.1)), base_q90=("base_size", lambda s: s.quantile(0.9)), fallback=("status", lambda s: float(s.isin(["refused", "unreadable", "no_foreground"]).mean())), near_click=("n_near_click", "sum")).reset_index()
    for k in C.BUDGETS:
        summ[f"p{k}"] = summ[f"p{k}"] / k
    summ.to_csv(C.OUT_DIR / "arm_summary.csv", index=False)
    print(summ.round(3).to_string(index=False))

    rows = []
    for click_set in ("joint", "refusal", "all"):
        cs = runs if click_set == "all" else runs[runs.click_set == click_set]
        for rank_key in ("tm_score", "chromatin_od"):
            rk = cs[cs.rank_key == rank_key]
            for fam, pairs in FAMILIES.items():
                if fam == "midog_box":
                    fam_rows = [dict(click_set=click_set, rank_key=rank_key, family=fam, ref="", pad=np.nan, **contrast(rk, a, b, k, rng)) for a, b in pairs for k in C.BUDGETS]
                else:
                    fam_rows = [dict(click_set=click_set, rank_key=rank_key, family=fam, ref=ref, pad=pad, **contrast(rk[(rk["ref"] == ref) & (rk["pad"] == pad)], a, b, k, rng)) for ref in ("sal", "gray") for pad in (0.15, 0.5, 1.0) for a, b in pairs for k in C.BUDGETS if len(rk[(rk["ref"] == ref) & (rk["pad"] == pad)])]
                adj = holm([r["p"] for r in fam_rows])
                for r, h in zip(fam_rows, adj):
                    r["holm_p"] = h
                rows += fam_rows
    con = pd.DataFrame(rows)
    con.to_csv(C.OUT_DIR / "contrasts.csv", index=False)
    show = ["click_set", "rank_key", "family", "ref", "pad", "pair", "K", "n_roi", "n_clicks", "prec_a", "prec_b", "delta", "ci_lo", "ci_hi", "p", "holm_p", "wins", "losses", "ties"]
    print(con[show].round(4).to_string(index=False))

    # saliency vs production at the MIDOG box, split by whether both pick the same object
    m = runs[runs.arm.isin(["sal_51", "prod_51"])].merge(geo[["file_name", "ann_id", "best_in_gray", "prod_status", "sal_n_kept"]], left_on=["file_name", "seed_ann_id"], right_on=["file_name", "ann_id"], how="left")
    m["object"] = np.where(m.prod_status != "refined", "gray refused", np.where(m.best_in_gray >= 0.5, "same object", "different object"))
    piv = m.pivot_table(index=["rank_key", "file_name", "seed_ann_id", "object", "sal_n_kept"], columns="arm", values=["tp_at_10", "tp_at_20", "tp_at_30"]).reset_index()
    print("\nsal_51 - prod_51 by object class (mean TP delta per click):")
    for (rk, obj), g in piv.groupby(["rank_key", "object"]):
        print(f"  {rk:13s} {obj:17s} n={len(g):4d}  " + "  ".join(f"K={k}: {(g[('tp_at_' + str(k), 'sal_51')] - g[('tp_at_' + str(k), 'prod_51')]).mean():+.3f}" for k in C.BUDGETS))

    # worst-ROI view at K=10: per-ROI mean precision per arm, MIDOG box
    j = runs[(runs.click_set == "joint") & runs.arm.isin(["default_51", "prod_51", "sal_51", "sal1_51"])]
    per = j.groupby(["rank_key", "arm", "file_name"]).tp_at_10.mean().div(10).unstack("arm")
    print("\nper-ROI precision@10 (joint clicks), worst 5 ROIs per arm and ROIs where an arm is >= 0.10 below prod_51:")
    for rk, g in per.groupby(level="rank_key"):
        g = g.droplevel("rank_key")
        print(f"  {rk}: min " + ", ".join(f"{a} {g[a].min():.2f}" for a in g.columns) + " | 5th-worst " + ", ".join(f"{a} {g[a].nsmallest(5).iloc[-1]:.2f}" for a in g.columns))
        for a in ("sal_51", "sal1_51"):
            print(f"    {a}: {int((g[a] - g['prod_51'] <= -0.10).sum())} ROIs >= 0.10 worse, {int((g[a] - g['prod_51'] >= 0.10).sum())} ROIs >= 0.10 better")


if __name__ == "__main__":
    main()
