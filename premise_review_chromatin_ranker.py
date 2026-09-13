#!/usr/bin/env python
"""Premise review: the chromatin-ranker adoption case (tm_score vs chromatin_od vs od_contrast).

Run with the anaconda interpreter from the repo root:

    /Users/mohinianand/anaconda3/bin/python3 premise_review_chromatin_ranker.py

Regenerates every table behind `Research Logs/2026-09-12-premise-review-chromatin-ranker`.
No claim here is supported by this repo's documentation; every number is recomputed from the
per-item CSVs in `results/` and `threshold_maxpeaks_ablation/`, from a synthetic microbenchmark,
or from statistical theory applied to those measurements.

The premises under test
-----------------------
P1  chromatin_od / od_contrast beat tm_score on precision@K by a measurable margin.
P2  max_peaks=50 is the right cap for a "minimal overhead" framing.
P3  od_contrast is dominated by chromatin_od at that cap.
P4  the ~40 ms re-ranking overhead is an intrinsic cost of the criterion.
P5  od51 is the right chromatin axis (five were measured, three were compared).
P6  the single-seed n=14 result is merely "not yet at D5's bar", i.e. under-powered but unbiased.
P7  D5's prior null and this positive result are about the same estimand.
"""

import os
import time

import cv2
import numpy as np
import pandas as pd
from scipy import stats

REPO = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(REPO, 'results')
ABL = os.path.join(REPO, 'threshold_maxpeaks_ablation')

CAP50_PREC = os.path.join(ABL, 'chromatin_od_ranker_max_peaks50_precision.csv')
CAP50_TIME = os.path.join(ABL, 'chromatin_od_ranker_max_peaks50_timing.csv')
CAP100_PREC = os.path.join(ABL, 'chromatin_od_ranker_precision.csv')
CAP100_TIME = os.path.join(ABL, 'chromatin_od_ranker_timing.csv')
UNBOUNDED = os.path.join(OUT, 'precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv')
UNBOUNDED_PREFIX = os.path.join(OUT, 'precision_at_k_14roi_prodseed_chromatin_raw.csv')
F5 = os.path.join(OUT, 'f5_nms_radius_ablation.csv')          # 14 ROI x 5 seed x 3 NMS radii
F1 = os.path.join(OUT, 'f1_seed_sweep.csv')                   # 7 ROI x 5 seed -- D5's dataset


def ci(d, alpha=0.05):
    """Two-sided t interval for the mean of a paired-difference vector."""
    d = np.asarray(pd.Series(d).dropna(), dtype=float)
    n = len(d)
    m = d.mean()
    se = d.std(ddof=1) / np.sqrt(n)
    tc = stats.t.ppf(1 - alpha / 2, n - 1)
    p = stats.ttest_1samp(d, 0).pvalue if n > 1 else np.nan
    return m, m - tc * se, m + tc * se, p, n


def paired(df, arm_a, arm_b, k, arm_col='arm', val='precision_at_budget', key='file_name'):
    a = df[(df[arm_col] == arm_a) & (df.budget == k)].set_index(key)[val]
    b = df[(df[arm_col] == arm_b) & (df.budget == k)].set_index(key)[val]
    idx = a.index.intersection(b.index)
    return (a.loc[idx] - b.loc[idx]).dropna()


# ---------------------------------------------------------------------------------
# Table 1 -- pooled precision@K under three candidate caps, same 14 ROIs / seed / z
# ---------------------------------------------------------------------------------
def table_pooled():
    cap50 = pd.read_csv(CAP50_PREC)
    cap100 = pd.read_csv(CAP100_PREC)
    unb = pd.read_csv(UNBOUNDED)

    # Pooling check: precision@K pooled (sum tp / sum delivered) equals the UNWEIGHTED mean of
    # per-ROI precision exactly, because every ROI delivers exactly K. Pooling is therefore not
    # doing any ROI-size weighting, and a paired per-ROI test is the matching inference.
    assert (cap50.budget_delivered == cap50.budget).all()
    assert (unb.budget_delivered == unb.budget).all()

    rows = []
    for tag, df, caps in [('cap50', cap50, 50), ('cap100', cap100, 100), ('unbounded', unb, np.nan)]:
        for k in sorted(set(df.budget) & {10, 20, 30, 50}):
            for arm in df.arm.unique():
                s = df[(df.arm == arm) & (df.budget == k)]
                if not len(s):
                    continue
                pooled = s.tp_at_budget.sum() / s.budget_delivered.sum()
                rows.append(dict(cap=tag, max_peaks=caps, arm=arm, K=k, n_roi=len(s),
                                 pooled_precision=pooled,
                                 unweighted_mean=s.precision_at_budget.mean(),
                                 pooled_equals_mean=bool(np.isclose(pooled, s.precision_at_budget.mean()))))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'premise_review_chromatin_ranker_pooled_precision.csv'), index=False)
    return t


# ---------------------------------------------------------------------------------
# Table 2 -- P1 at the single production seed: ROI-level and domain-level clustering
# ---------------------------------------------------------------------------------
def table_singleseed():
    cap50 = pd.read_csv(CAP50_PREC)
    unb = pd.read_csv(UNBOUNDED)
    dom = cap50.drop_duplicates('file_name').set_index('file_name').tumor_type

    rows = []
    for tag, df, ks in [('cap50', cap50, (10, 20, 30)), ('unbounded', unb, (10, 20, 30, 50))]:
        for a, b in [('chromatin_od', 'tm_score'), ('od_contrast', 'tm_score'),
                     ('od_contrast', 'chromatin_od')]:
            if a not in df.arm.unique():
                continue
            for k in ks:
                d = paired(df, a, b, k)
                m, lo, hi, p, n = ci(d)
                # domain-level: 7 tumour types x 2 ROIs -- the coarser exchangeable unit
                g = d.groupby(dom.reindex(d.index)).mean()
                gm, glo, ghi, gp, gn = ci(g)
                nz = d[d != 0]
                rows.append(dict(
                    dataset=tag, contrast=f'{a} - {b}', K=k,
                    roi_mean=m, roi_lo=lo, roi_hi=hi, roi_p=p, n_roi=n,
                    up=int((d > 0).sum()), down=int((d < 0).sum()), tie=int((d == 0).sum()),
                    wilcoxon_p=(stats.wilcoxon(nz).pvalue if len(nz) >= 3 else np.nan),
                    domain_mean=gm, domain_lo=glo, domain_hi=ghi, domain_p=gp, n_domain=gn))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'premise_review_chromatin_ranker_paired_singleseed.csv'), index=False)
    return t


# ---------------------------------------------------------------------------------
# Table 3 -- P6/P7: the multi-seed estimate, which is the estimand the product needs
# ---------------------------------------------------------------------------------
def table_multiseed():
    """f5 at the production NMS radius r7.5: 14 ROIs x 5 seeds, both criteria, paired.

    A user's click is a draw from the ROI's mitotic annotations, which is exactly what
    seed_index varies. The estimand for a product decision is therefore the expectation over
    seeds, not the value at one seed.
    """
    f5 = pd.read_csv(F5)
    f5 = f5[f5.arm.str.endswith('@r7.5')].copy()
    f5['crit'] = f5.arm.str.split('@').str[0]
    f5['prec'] = f5.tp_at_budget / f5.budget_delivered.replace(0, np.nan)

    rows = []
    for z in sorted(f5.z.unique()):
        for k in sorted(f5.budget.unique()):
            a = f5[(f5.crit == 'chromatin_od') & (f5.z == z) & (f5.budget == k)]
            b = f5[(f5.crit == 'tm_score') & (f5.z == z) & (f5.budget == k)]
            a = a.set_index(['file_name', 'seed_index']).prec
            b = b.set_index(['file_name', 'seed_index']).prec
            idx = a.index.intersection(b.index)
            d = (a.loc[idx] - b.loc[idx]).dropna()
            if len(d) < 5:
                continue
            per_roi = d.groupby(level=0).mean()     # average over seeds within ROI first
            m, lo, hi, p, n = ci(per_roi)
            per_seed = d.groupby(level=1).mean()
            rows.append(dict(
                z=z, K=k, n_pairs=len(d), n_roi=n, n_seed=per_seed.size,
                allseed_mean=m, allseed_lo=lo, allseed_hi=hi, allseed_p=p,
                seed0=per_seed.get(0, np.nan),
                seed_min=per_seed.min(), seed_max=per_seed.max(),
                seed_sd=per_seed.std(ddof=1),
                seed0_rank_of_n=int((per_seed >= per_seed.get(0, np.nan)).sum())))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'premise_review_chromatin_ranker_multiseed.csv'), index=False)
    return t


# ---------------------------------------------------------------------------------
# Table 4 -- the mechanism behind the seed effect: which criterion is seed-sensitive
# ---------------------------------------------------------------------------------
def table_seed_variance():
    f5 = pd.read_csv(F5)
    f5 = f5[f5.arm.str.endswith('@r7.5')].copy()
    f5['crit'] = f5.arm.str.split('@').str[0]
    f5['prec'] = f5.tp_at_budget / f5.budget_delivered.replace(0, np.nan)

    rows = []
    for k in [10, 25, 50, 100]:
        sub = f5[(f5.z == 1.0) & (f5.budget == k)]
        piv = {c: sub[sub.crit == c].pivot_table(index='file_name', columns='seed_index',
                                                 values='prec')
               for c in ['tm_score', 'chromatin_od']}
        if any(p.empty for p in piv.values()):
            continue
        # within-ROI seed residuals: removes the ROI main effect, leaving click-to-click spread
        res = {c: (piv[c].sub(piv[c].mean(axis=1), axis=0)).values.ravel() for c in piv}
        lev = stats.levene(res['tm_score'], res['chromatin_od'], center='median')
        worst = {c: piv[c].min(axis=1) for c in piv}
        dw = worst['chromatin_od'] - worst['tm_score']
        wm, wlo, whi, wp, _ = ci(dw)
        rows.append(dict(
            K=k,
            tm_withinroi_seed_sd=piv['tm_score'].std(axis=1, ddof=1).mean(),
            chrom_withinroi_seed_sd=piv['chromatin_od'].std(axis=1, ddof=1).mean(),
            tm_resid_var=np.var(res['tm_score'], ddof=1),
            chrom_resid_var=np.var(res['chromatin_od'], ddof=1),
            var_ratio=np.var(res['chromatin_od'], ddof=1) / np.var(res['tm_score'], ddof=1),
            levene_p=lev.pvalue,
            tm_worstseed=worst['tm_score'].mean(),
            chrom_worstseed=worst['chromatin_od'].mean(),
            worstseed_delta=wm, worstseed_lo=wlo, worstseed_hi=whi, worstseed_p=wp,
            worstseed_up=int((dw > 0).sum()), worstseed_down=int((dw < 0).sum())))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'premise_review_chromatin_ranker_seed_variance.csv'), index=False)
    return t


# ---------------------------------------------------------------------------------
# Table 5 -- P2: the cap 50 vs 100 trade, precision against measured latency
# ---------------------------------------------------------------------------------
def table_cap_tradeoff():
    cap50 = pd.read_csv(CAP50_PREC)
    cap100 = pd.read_csv(CAP100_PREC)
    t50 = pd.read_csv(CAP50_TIME)
    t100 = pd.read_csv(CAP100_TIME)

    rows = []
    for k in (10, 20, 30):
        d = paired(cap100, 'chromatin_od', 'chromatin_od', k)  # placeholder, replaced below
        a = cap100[(cap100.arm == 'chromatin_od') & (cap100.budget == k)].set_index('file_name').precision_at_budget
        b = cap50[(cap50.arm == 'chromatin_od') & (cap50.budget == k)].set_index('file_name').precision_at_budget
        d = (a - b).dropna()
        m, lo, hi, p, n = ci(d)
        rows.append(dict(K=k, contrast='chromatin_od cap100 - cap50',
                         mean=m, lo=lo, hi=hi, p=p, n=n,
                         up=int((d > 0).sum()), down=int((d < 0).sum()), tie=int((d == 0).sum()),
                         cap50_pooled=b.mean(), cap100_pooled=a.mean()))
    t = pd.DataFrame(rows)

    # latency side: the od51 loop is the only cap-dependent term; the pad is fixed cost
    lat = pd.DataFrame([
        dict(cap='cap50', n_detections_mean=t50.n_detections.mean(),
             pad_ms=t50.t6a_od_pad_s.mean() * 1e3,
             od51_loop_ms=t50.t6b_od51_loop_s.mean() * 1e3,
             overhead_ms=t50.t_chromatin_od_overhead_s.mean() * 1e3,
             shared_stages_ms=t50.t_shared_stages_s.mean() * 1e3),
        dict(cap='cap100', n_detections_mean=t100.n_detections.mean(),
             pad_ms=t100.t6a_od_pad_s.mean() * 1e3,
             od51_loop_ms=t100.t6b_od51_loop_s.mean() * 1e3,
             overhead_ms=t100.t_chromatin_od_overhead_s.mean() * 1e3,
             shared_stages_ms=t100.t_shared_stages_s.mean() * 1e3),
    ])
    lat['pad_share_of_overhead'] = lat.pad_ms / lat.overhead_ms
    t = t.merge(lat, how='cross')
    t.to_csv(os.path.join(OUT, 'premise_review_chromatin_ranker_cap_tradeoff.csv'), index=False)
    return t, lat


# ---------------------------------------------------------------------------------
# Table 6 -- P4: is the overhead intrinsic, or is it one avoidable whole-ROI copy?
# ---------------------------------------------------------------------------------
def table_pad_microbench(reps=7):
    """`chromatin_density` reads a `window`-sized box around each candidate. The notebooks pay
    for that by padding the ENTIRE ROI once (cv2.copyMakeBorder, BORDER_REPLICATE), then indexing
    into the padded copy. With only 50-100 candidates that is the wrong granularity: for every
    candidate at least `window//2` px from the border a direct slice is VALUE-IDENTICAL to the
    padded read, and for the minority nearer the border one replicate-pads the 51x51 patch
    (microseconds) instead of a 156 MB array. The clamped slice timed below is therefore a lower
    bound that is exact for interior candidates; it is not a drop-in for border candidates, which
    is why the correct fix is "pad the patch, not the ROI" rather than "drop the pad".

    Measured on a synthetic array of the production ROI's shape and dtype.
    """
    H, W = 5412, 7215
    rows = []
    for dtype in (np.float32, np.float64):
        a = np.random.rand(H, W).astype(dtype)
        pad_ts = []
        for _ in range(reps):
            t = time.perf_counter()
            p = cv2.copyMakeBorder(a, 25, 25, 25, 25, cv2.BORDER_REPLICATE)
            pad_ts.append(time.perf_counter() - t)
            del p
        rng = np.random.default_rng(0)
        for n in (50, 100):
            xs = rng.integers(0, W, n)
            ys = rng.integers(0, H, n)
            slice_ts = []
            for _ in range(reps):
                t = time.perf_counter()
                for x, y in zip(xs, ys):
                    x0, x1 = max(0, x - 25), min(W, x + 26)
                    y0, y1 = max(0, y - 25), min(H, y + 26)
                    w = a[y0:y1, x0:x1].ravel()
                    kk = max(1, int(round(0.10 * w.size)))
                    np.partition(w, -kk)[-kk:].mean()
                slice_ts.append(time.perf_counter() - t)
            rows.append(dict(dtype=np.dtype(dtype).name, array_mb=a.nbytes / 1e6, n_candidates=n,
                             copyMakeBorder_pad25_ms=np.median(pad_ts) * 1e3,
                             clamped_slice_ms=np.median(slice_ts) * 1e3,
                             pad_overhead_factor=np.median(pad_ts) / np.median(slice_ts)))
        del a
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'premise_review_chromatin_ranker_pad_microbench.csv'), index=False)
    return t


# ---------------------------------------------------------------------------------
# Table 7 -- P5: all five chromatin axes measured on the same 14 ROIs, not just three
# ---------------------------------------------------------------------------------
def table_axis_family():
    unb = pd.read_csv(UNBOUNDED)
    rows = []
    for arm in unb.arm.unique():
        if arm == 'tm_score':
            continue
        for k in (10, 20, 30, 50):
            d_tm = paired(unb, arm, 'tm_score', k)
            d_od51 = paired(unb, arm, 'chromatin_od', k)
            m1, lo1, hi1, p1, _ = ci(d_tm)
            m2, lo2, hi2, p2, _ = ci(d_od51)
            s = unb[(unb.arm == arm) & (unb.budget == k)]
            rows.append(dict(arm=arm, K=k,
                             pooled_precision=s.tp_at_budget.sum() / s.budget_delivered.sum(),
                             vs_tm_mean=m1, vs_tm_lo=lo1, vs_tm_hi=hi1, vs_tm_p=p1,
                             vs_od51_mean=m2, vs_od51_lo=lo2, vs_od51_hi=hi2, vs_od51_p=p2))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'premise_review_chromatin_ranker_axis_family.csv'), index=False)
    return t


# ---------------------------------------------------------------------------------
# Table 8 -- P7: reproduce D5's stated null from its own raw data, then change the metric
# ---------------------------------------------------------------------------------
def table_d5_reproduction():
    f1 = pd.read_csv(F1)
    f1['prec'] = f1.tp_at_budget / f1.budget_delivered.replace(0, np.nan)
    rows = []
    for k in sorted(f1.budget.unique()):
        for metric, col in [('recall', 'recall_at_budget'), ('precision', 'prec')]:
            a = f1[(f1.arm == 'chromatin_od') & (f1.z == 1.0) & (f1.budget == k)]
            b = f1[(f1.arm == 'tm_score') & (f1.z == 1.0) & (f1.budget == k)]
            a = a.set_index(['file_name', 'seed_index'])[col]
            b = b.set_index(['file_name', 'seed_index'])[col]
            idx = a.index.intersection(b.index)
            d = (a.loc[idx] - b.loc[idx]).dropna().groupby(level=0).mean()
            if len(d) < 3:
                continue
            m, lo, hi, p, n = ci(d)
            rows.append(dict(metric=metric, K=k, z=1.0, n_roi=n, n_seed=5,
                             mean=m, lo=lo, hi=hi, p=p, roi_positive=int((d > 0).sum())))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'premise_review_chromatin_ranker_d5_reproduction.csv'), index=False)
    return t


def table_seed_transport():
    """How much of the single-seed effect survives is NOT identified by the data available.

    f5 (5 seeds) and the halfpixfix/cap runs (1 seed) are not the same operating point: at the
    SAME seed and the same 14 ROIs they disagree by 2x on the same contrast (f5 seed0 +0.064 vs
    halfpixfix seed0 +0.121), and f5's tm_score baseline runs 0.079 higher. So f5's seed-0
    optimism can be transported to the production config additively or multiplicatively, and the
    two give answers ~4x apart. This table states that fork rather than hiding it.
    """
    f5 = pd.read_csv(F5)
    f5 = f5[f5.arm.str.endswith('@r7.5')].copy()
    f5['crit'] = f5.arm.str.split('@').str[0]
    f5['prec'] = f5.tp_at_budget / f5.budget_delivered.replace(0, np.nan)
    cap50 = pd.read_csv(CAP50_PREC)

    rows = []
    for k_f5, k_prod in [(10, 10), (25, 20), (50, 30)]:
        a = f5[(f5.crit == 'chromatin_od') & (f5.z == 1.0) & (f5.budget == k_f5)]
        b = f5[(f5.crit == 'tm_score') & (f5.z == 1.0) & (f5.budget == k_f5)]
        a = a.set_index(['file_name', 'seed_index']).prec
        b = b.set_index(['file_name', 'seed_index']).prec
        idx = a.index.intersection(b.index)
        d = (a.loc[idx] - b.loc[idx]).dropna()
        per_seed = d.groupby(level=1).mean()
        f5_seed0, f5_all = per_seed.get(0, np.nan), per_seed.mean()
        obs = paired(cap50, 'chromatin_od', 'tm_score', k_prod).mean()
        rows.append(dict(
            f5_K=k_f5, prod_K=k_prod,
            f5_seed0=f5_seed0, f5_allseed=f5_all,
            seed0_optimism_additive=f5_seed0 - f5_all,
            seed0_retention_multiplicative=(f5_all / f5_seed0 if f5_seed0 else np.nan),
            observed_cap50_seed0=obs,
            transported_additive=obs - (f5_seed0 - f5_all),
            transported_multiplicative=obs * (f5_all / f5_seed0 if f5_seed0 else np.nan)))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, 'premise_review_chromatin_ranker_seed_transport.csv'), index=False)
    return t


def main():
    pd.set_option('display.width', 200)
    pd.set_option('display.max_columns', 40)

    print('=' * 78)
    print('P1/P3  pooled precision@K under three caps (same 14 ROIs, seed 0, z=-1.5)')
    print('=' * 78)
    t = table_pooled()
    print(t.pivot_table(index=['cap', 'K'], columns='arm', values='pooled_precision').round(4).to_string())
    print('\npooled == unweighted mean on every row:', bool(t.pooled_equals_mean.all()),
          '  (every ROI delivers exactly K, so pooling applies no ROI weighting)')

    print('\n' + '=' * 78)
    print('P1  paired inference at the single production seed (ROI- and domain-clustered)')
    print('=' * 78)
    ss = table_singleseed()
    print(ss[['dataset', 'contrast', 'K', 'roi_mean', 'roi_lo', 'roi_hi', 'roi_p',
              'up', 'down', 'tie', 'domain_mean', 'domain_lo', 'domain_hi']].round(4).to_string(index=False))

    print('\n' + '=' * 78)
    print('P6  the multi-seed estimate (f5 @ r7.5, 14 ROIs x 5 seeds) -- the product estimand')
    print('=' * 78)
    ms = table_multiseed()
    show = ms[ms.K.isin([10, 25, 50]) & ms.z.isin([0.5, 1.0, 2.0])]
    print(show[['z', 'K', 'n_pairs', 'allseed_mean', 'allseed_lo', 'allseed_hi', 'allseed_p',
                'seed0', 'seed_min', 'seed_max', 'seed_sd', 'seed0_rank_of_n']].round(4).to_string(index=False))
    print('\nseed0_rank_of_n == 1 means the production seed was the MOST favourable of the five.')

    print('\n' + '=' * 78)
    print('P6 mechanism  which criterion is seed-sensitive?')
    print('=' * 78)
    sv = table_seed_variance()
    print(sv.round(5).to_string(index=False))

    print('\n' + '=' * 78)
    print('P2  cap 50 vs cap 100: precision cost against measured latency saving')
    print('=' * 78)
    ct, lat = table_cap_tradeoff()
    print(ct[['K', 'contrast', 'cap50_pooled', 'cap100_pooled', 'mean', 'lo', 'hi', 'p',
              'up', 'down', 'tie']].drop_duplicates().round(4).to_string(index=False))
    print('\nlatency, mean over 14 ROIs (ms):')
    print(lat.round(3).to_string(index=False))

    print('\n' + '=' * 78)
    print('P4  is the overhead intrinsic? synthetic pad vs clamped-slice microbenchmark')
    print('=' * 78)
    print(table_pad_microbench().round(3).to_string(index=False))

    print('\n' + '=' * 78)
    print('P5  all five chromatin axes on the same 14 ROIs (unbounded pool)')
    print('=' * 78)
    af = table_axis_family()
    print(af[af.K.isin([10, 20, 30])][
        ['arm', 'K', 'pooled_precision', 'vs_tm_mean', 'vs_tm_p', 'vs_od51_mean',
         'vs_od51_lo', 'vs_od51_hi']].round(4).to_string(index=False))

    print('\n' + '=' * 78)
    print("P7  reproduce D5's null from f1_seed_sweep, then vary the metric")
    print('=' * 78)
    d5 = table_d5_reproduction()
    print(d5[d5.K.isin([25, 50, 100, 250])].round(4).to_string(index=False))
    r250 = d5[(d5.metric == 'recall') & (d5.K == 250)].iloc[0]
    print(f"\nreproduced Delta recall@250 = {r250['mean']:+.3f} "
          f"95% CI [{r250['lo']:+.3f}, {r250['hi']:+.3f}] p={r250['p']:.2f} "
          f"positive on {int(r250['roi_positive'])} of {int(r250['n_roi'])} ROIs")
    print("D5 as quoted in midog_utils/chromatin.py: +0.032, CI [-0.047, +0.112], p=0.36, 3 of 7")

    print('\n' + '=' * 78)
    print('P6 transport  the single-seed effect under two defensible transports of f5 optimism')
    print('=' * 78)
    tr = table_seed_transport()
    print(tr.round(4).to_string(index=False))
    print('\nThe two transports disagree by ~4x. The LEVEL of the gain at the production config'
          ' is therefore not identified by any data in this repo; only the mechanism transports.')

    print('\nwrote results/premise_review_chromatin_ranker_*.csv')


if __name__ == '__main__':
    main()
