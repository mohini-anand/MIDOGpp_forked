"""F7 stage 2: apply the size prune across all arms and print every table.

Run from the repository root or from `f7_size_prune/` -- all shared artefact paths are
absolute. Needs `results/f7_repro_gate.csv` and `results/f7_candidate_sizes.npz`, which
`f7_size_features.py` writes.

    /Users/mohinianand/anaconda3/bin/python3.11 f7_size_prune/f7_analysis.py

Takes ~2 minutes on an idle machine. It is CPU-bound in the greedy re-match, so it crawls
if anything else large is running.
"""
import os, sys, time
import numpy as np, pandas as pd
from scipy.stats import binomtest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
import f7_prune_eval as f7

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 60)


def main():
    t0 = time.time()
    gate, sizes, pool, ledger, ann = f7.load_all()
    print(f"loaded in {time.time()-t0:.0f}s :: {len(gate)} cells, {gate.file_name.nunique()} ROIs, "
          f"{gate.tumor_type.nunique()} domains", flush=True)
    print(f"REPRODUCTION GATE: ranks exact {int(gate.ranks_exact.sum())}/{len(gate)}; "
          f"n_gt_mitotic agrees {int((gate.n_gt_mitotic==gate.n_gt_mitotic_csv).sum())}/{len(gate)}")
    assert gate.ranks_exact.all(), "gate failed -- every number below would be void"
    print("distinct clicks/ROI (seed draw is with replacement):",
          gate.groupby('file_name').seed_ann_id.nunique().to_dict())

    print("\n=== INSTRUMENT: on_nucleus, the cap on what the knob can act on ===")
    print(gate.groupby(["tumor_type", "file_name"]).agg(
        n_pool=("n_pool", "mean"), on_nucleus=("frac_on_nucleus", "mean"),
        n_comp=("n_components", "first")).round(3).to_string())

    print("\n=== SIZES (equivalent diameter, um) ===")
    rows = []
    for _, g in gate.iterrows():
        fn, si = g.file_name, int(g.seed_index)
        a = sizes[f"{fn}|{si}|eqd_um"]; a = a[~np.isnan(a)]
        t = f7.tp_sizes_for_cell(sizes, ledger, fn, si)
        rows.append(dict(file_name=fn, tumor_type=g.tumor_type, n_tp_measurable=len(t),
                         n_tp_total=int(g.n_gt_mitotic), tp_min=np.min(t),
                         tp_p05=np.percentile(t, 5), tp_med=np.median(t),
                         tp_p95=np.percentile(t, 95), tp_max=np.max(t),
                         pool_med=np.median(a), pool_p95=np.percentile(a, 95)))
    print(pd.DataFrame(rows).groupby(["tumor_type", "file_name"])
          .median(numeric_only=True).round(2).to_string())

    # the fast matcher must agree with ev.bucket_detections before it is used for anything
    for fn, si in [("013.tiff", 0), ("301.tiff", 2), ("459.tiff", 4)]:
        g = gate[(gate.file_name == fn) & (gate.seed_index == si)].iloc[0]
        cx, cy = pool[f"{fn}|{si}|cx"], pool[f"{fn}|{si}|cy"]
        for lo, hi in [(-np.inf, np.inf), (f7.APRIORI_LO, f7.APRIORI_HI)]:
            f7.verify_fast_path(cx, cy, sizes[f"{fn}|{si}|eqd_um"],
                                f7.gt_eval_for(ann, fn, int(g.seed_ann_id)),
                                float(g.match_radius_px), int(g.n_gt_mitotic), lo, hi)
    print("\nfast matcher verified against ev.bucket_detections on 3 cells x 2 arms", flush=True)

    tps = {(g.file_name, int(g.seed_index)):
           f7.tp_sizes_for_cell(sizes, ledger, g.file_name, int(g.seed_index))
           for _, g in gate.iterrows()}
    doms = gate.groupby("tumor_type")["file_name"].unique().to_dict()
    lodo_pool = {d: np.concatenate([v for (f2, s2), v in tps.items() if f2 not in rn])
                 for d, rn in doms.items()}

    print("\n=== RUNNING ALL ARMS ===", flush=True)
    t1, rows = time.time(), []
    for _, g in gate.iterrows():
        fn, si = g.file_name, int(g.seed_index)
        cx, cy = pool[f"{fn}|{si}|cx"], pool[f"{fn}|{si}|cy"]
        sz = sizes[f"{fn}|{si}|eqd_um"]
        ge = f7.gt_eval_for(ann, fn, int(g.seed_ann_id)); n_mit = int(g.n_gt_mitotic)
        neigh, gt_xy, det_xy = f7.cell_neighbours(cx, cy, ge, float(g.match_radius_px))
        gt_cls = ge["category_id"].to_numpy()
        base = dict(file_name=fn, tumor_type=g.tumor_type, seed_index=si, n_mit=n_mit,
                    n_pool=int(g.n_pool), frac_on_nucleus=float(g.frac_on_nucleus))
        arms = [("baseline", -np.inf, np.inf, np.nan),
                ("apriori", f7.APRIORI_LO, f7.APRIORI_HI, np.nan)]
        for r in f7.QUANTS:
            arms.append((f"insample@{r:.2f}", *f7.bounds_from_tp_sizes(tps[(fn, si)], r), r))
            arms.append((f"lodo@{r:.2f}", *f7.bounds_from_tp_sizes(lodo_pool[g.tumor_type], r), r))
        for name, lo, hi, r in arms:
            rows.append({**base, "arm": name, "lo_um": lo, "hi_um": hi, "retention_nominal": r,
                         **f7.score_arm_fast(neigh, gt_xy, det_xy, gt_cls,
                                             f7.keep_mask(sz, lo, hi), n_mit)})
        print(f"  {fn} s{si} ({time.time()-t1:.0f}s)", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(f7.RES, "f7_prune_results.csv"), index=False)
    print(f"all arms in {time.time()-t1:.0f}s")

    print("\n=== SUMMARY BY ARM (cell-level medians) ===")
    S = []
    for arm in (["apriori"] + [f"insample@{r:.2f}" for r in f7.QUANTS]
                + [f"lodo@{r:.2f}" for r in f7.QUANTS]):
        t = f7.paired_table(res, arm)
        S.append(dict(arm=arm, lo=t.lo_um.median(), hi=t.hi_um.median(),
                      frac_kept=t.frac_kept.median(), frac_kept_min=t.frac_kept.min(),
                      d_recall_med=t.d_recall.median(), d_recall_worst=t.d_recall.min(),
                      cells_losing_tp=int((t.d_recall < 0).sum()),
                      d_d90=t.d_depth_90.median(), d_d95=t.d_depth_95.median(),
                      d_d100=t.d_depth_100.median()))
    print(pd.DataFrame(S).round(4).to_string(index=False))

    for P in ["lodo@1.00", "apriori", "insample@1.00", "lodo@0.98"]:
        t = f7.paired_table(res, P)
        roi = t.groupby(["tumor_type", "file_name"]).median(numeric_only=True)
        for q in (90, 95, 100):
            roi[f"d{q}"] = roi[f"depth_{q}"] - roi[f"depth_{q}_base"]
        print(f"\n\n##### ARM {P}   bound [{t.lo_um.median():.2f}, {t.hi_um.median():.2f}] um #####")
        print(roi[["n_pool", "frac_kept", "recall_base", "recall", "d_recall",
                   "depth_90_base", "depth_90", "depth_95_base", "depth_95",
                   "depth_100_base", "depth_100"]].round(3).to_string())
        for q in (90, 95, 100):
            dd = roi[f"d{q}"].dropna(); b = int((dd < 0).sum()); n = len(dd)
            print(f"  depth_{q}: shorter on {b}/{n} ROIs, median {dd.median():+.0f} candidates, "
                  f"p={binomtest(b, n, 0.5).pvalue:.4f} [floor {2/2**n:.4f}]")
        print(f"  ROIs losing >=1 mitosis: {int((roi.d_recall < 0).sum())}/{len(roi)}; "
              f"worst {roi.d_recall.min():+.4f}")
        c = roi[["depth_100_base", "d100"]].dropna().copy()
        c["rel"] = -c.d100 / c.depth_100_base
        print(f"  Spearman(depth100_base, relative saving) = "
              f"{c.depth_100_base.corr(c.rel, method='spearman'):+.3f}  "
              f"(negative => pays least where the burden is worst)")
        roi.to_csv(os.path.join(f7.RES, f"f7_roi_{P.replace('@','_').replace('.','')}.csv"))

    print("\n\n=== LODO FOLD BOUNDS (the universality check) ===")
    F = []
    for dm, rn in doms.items():
        held = np.concatenate([v for (f, s), v in tps.items() if f in rn])
        lo, hi = f7.bounds_from_tp_sizes(lodo_pool[dm], 1.00)
        F.append(dict(domain=dm, n_tp_held=len(held), lo_um=lo, hi_um=hi,
                      held_min=held.min(), held_max=held.max(),
                      frac_held_deleted=float(np.mean((held < lo) | (held > hi)))))
    print(pd.DataFrame(F).round(3).to_string(index=False))

    print("\n=== RETENTION CURVE (ROI-level) ===")
    C = []
    for r in f7.QUANTS:
        for kind in ("insample", "lodo"):
            t2 = f7.paired_table(res, f"{kind}@{r:.2f}")
            r2 = t2.groupby("file_name").median(numeric_only=True)
            C.append(dict(regime=kind, retention=r, lo=t2.lo_um.median(), hi=t2.hi_um.median(),
                          frac_kept=r2.frac_kept.median(), worst_d_recall=r2.d_recall.min(),
                          rois_losing_tp=int((r2.d_recall < 0).sum()),
                          d_d90=r2.d_depth_90.median(), d_d95=r2.d_depth_95.median(),
                          d_d100=r2.d_depth_100.median(),
                          rois_shorter100=int((r2.d_depth_100 < 0).sum())))
    cv = pd.DataFrame(C)
    cv.to_csv(os.path.join(f7.RES, "f7_retention_curve.csv"), index=False)
    print(cv.round(4).to_string(index=False))
    print("\nDONE")


if __name__ == "__main__":
    main()
