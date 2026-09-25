"""
    D12's own reversal test: does a box-sized window (the drawn box as the gate window and fallback) do worse than the
    fixed 51 px window at the same clicks? prod_box - prod_51 and, for comparison, sal1_box - sal1_51, on the 294 joint
    clicks; Holm within each rank key x pair (2 references x 3 paddings x 3 K = 18 tests); D5 = CI excludes 0 and a
    majority of all 49 ROIs moves in the pooled direction. Writes out/box_vs_click.csv.
"""

import numpy as np
import pandas as pd

import common as C
from analyze import contrast, holm


def main():
    """
        Compute, print and write the box-vs-click contrasts.

        Returns None.
    """
    r = pd.read_csv(C.OUT_DIR / "precision_per_run.csv")
    r["ref"], r["pad"] = r["ref"].fillna(""), r["pad"].fillna(-1.0)
    r = r[r.click_set == "joint"]
    rng = np.random.default_rng(6)
    rows = []
    for rk in ("tm_score", "chromatin_od"):
        for box, click in (("prod_box", "prod_51"), ("sal1_box", "sal1_51")):
            fam = []
            for ref in ("sal", "gray"):
                for pad in (0.15, 0.5, 1.0):
                    b = r[(r.rank_key == rk) & (r.arm == box) & (r["ref"] == ref) & (r["pad"] == pad)].assign(arm="x_51", ref="", pad=-1.0)
                    c = r[(r.rank_key == rk) & (r.arm == click)].assign(arm="y_51")
                    for k in C.BUDGETS:
                        fam.append({**contrast(pd.concat([b, c]), "x_51", "y_51", k, rng), "rank_key": rk, "pair": f"{box} - {click}", "ref": ref, "pad": pad})
            for f, h in zip(fam, holm([f["p"] for f in fam])):
                f["holm_p"] = h
                same = f["wins"] if f["delta"] > 0 else f["losses"]
                f["d5_met"] = bool((f["ci_lo"] > 0 or f["ci_hi"] < 0) and same > f["n_roi"] / 2)
            rows += fam
    out = pd.DataFrame(rows)
    out.to_csv(C.OUT_DIR / "box_vs_click.csv", index=False)
    pd.set_option("display.width", 250)
    print(out[["rank_key", "pair", "ref", "pad", "K", "prec_a", "prec_b", "delta", "ci_lo", "ci_hi", "holm_p", "wins", "losses", "ties", "d5_met"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
