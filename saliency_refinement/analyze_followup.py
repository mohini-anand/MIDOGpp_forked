"""
    Statistics for out/precision_followup_per_run.csv, with the joint clicks' default_51 / prod_51 rows taken from
    out/precision_per_run.csv. Same machinery as analyze.py (ROI as the unit, exact sign-flip, bootstrap CI, Holm within
    each click set x rank key). Writes out/contrasts_followup.csv and prints it.
"""

import numpy as np
import pandas as pd

import common as C
from analyze import contrast, holm

PAIRS = {
    "joint": [("hem_mask_fill_1um", "default_51"), ("sal_fill_51", "default_51"), ("sal_fill_51", "hem_mask_fill_1um"), ("sal_fill_51", "prod_51")],
    "gray_refused": [("sal_51", "default_51"), ("sal1_51", "default_51"), ("sal_fill_51", "default_51")],
}


def main():
    """
        Print and write the follow-up contrasts, arm means and the different-object clicks.

        Returns None.
    """
    fu = pd.read_csv(C.OUT_DIR / "precision_followup_per_run.csv")
    main_runs = pd.read_csv(C.OUT_DIR / "precision_per_run.csv")
    base = main_runs[(main_runs.click_set == "joint") & main_runs.arm.isin(["default_51", "prod_51"])]
    runs = pd.concat([fu, base[fu.columns.intersection(base.columns)]], ignore_index=True)
    runs["ref"], runs["pad"] = "", -1.0
    rng = np.random.default_rng(20260925)
    pd.set_option("display.width", 250)

    print(runs.groupby(["click_set", "rank_key", "arm"]).agg(n=("tp_at_10", "size"), p10=("tp_at_10", "mean"), p20=("tp_at_20", "mean"), p30=("tp_at_30", "mean"), base_med=("base_size", "median")).assign(p10=lambda d: d.p10 / 10, p20=lambda d: d.p20 / 20, p30=lambda d: d.p30 / 30).round(3).to_string())
    rows = []
    for click_set, pairs in PAIRS.items():
        for rank_key in ("tm_score", "chromatin_od"):
            rk = runs[(runs.click_set == click_set) & (runs.rank_key == rank_key)]
            fam = [dict(click_set=click_set, rank_key=rank_key, **contrast(rk.assign(arm=rk.arm.str.replace("hem_mask_fill_1um", "hmf_51")), a.replace("hem_mask_fill_1um", "hmf_51"), b.replace("hem_mask_fill_1um", "hmf_51"), k, rng)) for a, b in pairs for k in C.BUDGETS]
            for r, h in zip(fam, holm([r["p"] for r in fam])):
                r["holm_p"] = h
            rows += fam
    con = pd.DataFrame(rows)
    con.to_csv(C.OUT_DIR / "contrasts_followup.csv", index=False)
    print(con[["click_set", "rank_key", "pair", "K", "n_roi", "n_clicks", "prec_a", "prec_b", "delta", "ci_lo", "ci_hi", "p", "holm_p", "wins", "losses", "ties"]].round(4).to_string(index=False))

    d = runs[runs.click_set == "different_object"].pivot_table(index=["file_name", "seed_ann_id", "rank_key"], columns="arm", values="tp_at_10")
    print("\ndifferent-object clicks, TP@10:\n" + d.to_string())
    g = runs[runs.click_set == "gray_refused"]
    print(f"\ngray_refused: {g.groupby(['file_name', 'seed_ann_id']).ngroups} clicks on {g.file_name.nunique()} ROIs")


if __name__ == "__main__":
    main()
