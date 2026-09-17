"""
    Compare the post-D11 re-run of the 49-ROI x 3-seed bbox-refinement experiment against the
    committed pre-D11 reference CSVs. Reads both, writes nothing over the reference. Run from
    the repo root with /Users/mohinianand/anaconda3/bin/python3.

    Answers, in order: (1) the n_detections delta distribution against the mechanism the
    stored n_self_hits column predicts; (2) mean precision@K per arm, old vs new; (3) every
    (ROI, seed, condition, K) cell whose TP count moved, with the rows that entered and left
    the top K; (4) whether the three arms' ranking and the sign of their differences survive.
"""

from __future__ import annotations

import itertools
import math
import sys
from fractions import Fraction

import numpy as np
import pandas as pd

sys.path.insert(0, '.')

pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 80)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_colwidth', 140)

REF_STEM = 'results/precision_at_k_49roi_3seed_chromatin_bbox3way'
NEW_STEM = f'{REF_STEM}_postD11'
CONDITIONS = ('default_51', 'gray_bbox', 'hem_bbox')
PAIRS = [('gray_bbox', 'default_51'), ('hem_bbox', 'default_51'), ('hem_bbox', 'gray_bbox')]
BUDGETS = (10, 20, 30)
SEED_INDICES = (0, 1, 2)
SUBSETS = ('all_49', 'original_14', 'testing_35')
HOLM_FAMILIES = ('all_49', 'testing_35')
N_BOOT = 10_000
BOOT_SEED = 20260916
KEY = ['file_name', 'seed_index', 'condition']


def sign_flip_null_counts(abs_deltas):
    """
        Exact null distribution of sum(+/-|d_i|) by convolution over integer offsets.

        abs_deltas (array-like): nonzero absolute integer deltas, one per ROI.

        Returns np.ndarray: int64 counts of each sum -T..T (index = sum + T, T = sum |d_i|).
    """
    d = [int(v) for v in abs_deltas]
    assert all(v > 0 for v in d) and len(d) <= 62, 'nonzero integer deltas only, and 2 ** n must fit int64'
    total = sum(d)
    counts = np.zeros(2 * total + 1, dtype=np.int64)
    counts[total] = 1
    for v in d:
        shifted = np.zeros_like(counts)
        shifted[v:] += counts[:len(counts) - v]
        shifted[:len(counts) - v] += counts[v:]
        counts = shifted
    return counts


def exact_sign_flip_p(deltas):
    """
        Two-sided exact sign-flip p on integer per-ROI deltas, by dynamic-programming convolution.

        deltas (array-like): one integer delta per ROI; zeros carry no sign.

        Returns tuple[float, int]: the p-value, and the number of nonzero deltas.
    """
    d = np.asarray(deltas)
    assert np.all(d == np.round(d)), 'the convolution needs integer deltas'
    d = d.astype(np.int64)
    nz = np.abs(d[d != 0])
    if len(nz) == 0:
        return 1.0, 0
    counts = sign_flip_null_counts(nz)
    total = int(nz.sum())
    observed = abs(int(d.sum()))
    sums = np.arange(-total, total + 1)
    return float(Fraction(int(counts[np.abs(sums) >= observed].sum()), 2 ** len(nz))), len(nz)


def brute_force_sign_flip_p(deltas):
    """
        Two-sided sign-flip p by enumerating every sign pattern, kept only to validate the convolution.

        deltas (array-like): one delta per ROI; zeros carry no sign and are dropped.

        Returns float: the p-value.
    """
    d = np.asarray(deltas, dtype=np.float64)
    nz = np.abs(d[np.abs(d) > 1e-12])
    if len(nz) == 0:
        return 1.0
    signs = np.array(list(itertools.product([1.0, -1.0], repeat=len(nz))))
    return float((np.abs(signs @ nz) >= abs(float(d.sum())) - 1e-9).mean())


def holm_adjust(pvalues):
    """
        Holm step-down adjusted p-values, valid under any dependence between the tests.

        pvalues (array-like): raw p-values of one family of tests.

        Returns np.ndarray: adjusted p-values, in the input order.
    """
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p, kind='stable')
    adjusted = np.maximum.accumulate((len(p) - np.arange(len(p))) * p[order]).clip(max=1.0)
    out = np.empty_like(adjusted)
    out[order] = adjusted
    return out


def subset_frame(frame, subset):
    """
        The rows of a frame that belong to one analysis subset.

        frame (pd.DataFrame): a frame with a subset column.
        subset (str): "all_49", "original_14" or "testing_35".

        Returns pd.DataFrame: every row for all_49, else that subset's rows.
    """
    return frame if subset == 'all_49' else frame[frame['subset'] == subset]


def delta_per_roi(results, roi_order, roi_meta):
    """
        The ROI-level paired statistic: TP@K summed over the three clicks, condition a minus condition b.

        results (pd.DataFrame): per-run rows.
        roi_order (list): file names in analysis order.
        roi_meta (pd.DataFrame): file_name-indexed subset and domain.

        Returns pd.DataFrame: one row per (ROI, pair, K).
    """
    tp_seed = results.set_index(KEY)[[f'tp_at_{k}' for k in BUDGETS]]
    rows = []
    for a, b in PAIRS:
        for k in BUDGETS:
            for fn in roi_order:
                d_seed = [int(tp_seed.loc[(fn, s, a), f'tp_at_{k}'] - tp_seed.loc[(fn, s, b), f'tp_at_{k}']) for s in SEED_INDICES]
                rows.append(dict(subset=roi_meta.loc[fn, 'subset'], domain=roi_meta.loc[fn, 'domain'], file_name=fn, pair=f'{a} - {b}', K=k, delta_tp_sum=sum(d_seed), **{f'delta_tp_s{s}': v for s, v in zip(SEED_INDICES, d_seed)}))
    return pd.DataFrame(rows)


def boot_ci_exact(int_values, idx, denominator):
    """
        Percentile 95% CI of a pooled ratio over resampled ROIs, from exact integer numerators.

        int_values (array-like): one integer numerator per ROI.
        idx (np.ndarray): (n_boot, n_roi) resampled ROI indices.
        denominator (int): the pooled denominator per resample.

        Returns tuple[float, float]: CI low and CI high.
    """
    sums = np.asarray(int_values, dtype=np.int64)[idx].sum(axis=1)
    lo, hi = np.percentile(sums / denominator, [2.5, 97.5])
    return float(lo), float(hi)


def delta_stats(results, roi_order, roi_meta, rois_by_subset):
    """
        The pooled delta, win/loss/tie counts, exact sign-flip p, Holm p and bootstrap CI for every pair and K.

        results (pd.DataFrame): per-run rows.
        roi_order (list): file names in analysis order.
        roi_meta (pd.DataFrame): file_name-indexed subset and domain.
        rois_by_subset (dict): subset -> its file names in analysis order.

        Returns tuple[pd.DataFrame, pd.DataFrame]: the per-(subset, pair, K) statistics, and the per-ROI deltas.
    """
    dpr = delta_per_roi(results, roi_order, roi_meta)
    boot_idx = {subset: np.random.default_rng([BOOT_SEED, i]).integers(len(rois_by_subset[subset]), size=(N_BOOT, len(rois_by_subset[subset]))) for i, subset in enumerate(SUBSETS)}
    delta_ix = dpr.set_index(['file_name', 'pair', 'K'])
    rows = []
    for subset in SUBSETS:
        order = rois_by_subset[subset]
        for a, b in PAIRS:
            for k in BUDGETS:
                d = delta_ix.loc[[(fn, f'{a} - {b}', k) for fn in order], 'delta_tp_sum'].to_numpy()
                p, n_nz = exact_sign_flip_p(d)
                denominator = len(order) * len(SEED_INDICES) * k
                lo, hi = boot_ci_exact(d, boot_idx[subset], denominator)
                pooled = d.sum() / denominator
                direction = np.sign(pooled)
                n_same = int((d > 0).sum()) if direction > 0 else (int((d < 0).sum()) if direction < 0 else 0)
                rows.append(dict(subset=subset, pair=f'{a} - {b}', K=k, n_roi=len(order), pooled_delta_precision=pooled, wins=int((d > 0).sum()), losses=int((d < 0).sum()), ties=int((d == 0).sum()), n_nonzero=n_nz, exact_p=p, ci_low=lo, ci_high=hi, n_same_direction=n_same, ci_excludes_0=bool(lo > 0 or hi < 0), majority_of_all_rois=bool(n_same > len(order) / 2)))
    stats = pd.DataFrame(rows)
    stats['d5_bar_met'] = stats['ci_excludes_0'] & stats['majority_of_all_rois']
    stats['holm_p'] = np.nan
    for family in HOLM_FAMILIES:
        in_family = stats['subset'] == family
        assert int(in_family.sum()) == len(PAIRS) * len(BUDGETS)
        stats.loc[in_family, 'holm_p'] = holm_adjust(stats.loc[in_family, 'exact_p'])
    return stats, dpr


def topk_membership(top, k):
    """
        The set of (cx, cy) positions each run puts in its top K, with each position's bucket and matched annotation.

        top (pd.DataFrame): top-30 rows for every run.
        k (int): the budget.

        Returns dict: (file_name, seed_index, condition) -> {(cx, cy): (rank, bucket, matched_ann_id)}.
    """
    out = {}
    for key, g in top[top['rank'] < k].groupby(KEY, sort=False):
        out[key] = {(float(r.cx), float(r.cy)): (int(r.rank), r.bucket, r.matched_ann_id) for r in g.itertuples()}
    return out


def col(frame, name):
    """
        One merged column, whether or not the join suffixed it.

        frame (pd.DataFrame): the merged old/new frame.
        name (str): the unsuffixed column name.

        Returns pd.Series: the column, preferring the pre-D11 side when both sides carry it.
    """
    return frame[name] if name in frame else frame[f'{name}_old']


def main():
    """
        Run the whole comparison and print it.

        Returns None.
    """
    old = pd.read_csv(f'{REF_STEM}_per_run.csv', float_precision='round_trip')
    new = pd.read_csv(f'{NEW_STEM}_per_run.csv', float_precision='round_trip')
    old_top = pd.read_csv(f'{REF_STEM}_top30.csv', float_precision='round_trip')
    new_top = pd.read_csv(f'{NEW_STEM}_top30.csv', float_precision='round_trip')
    assert len(old) == len(new) == 441, f'expected 441 runs each, got old {len(old)}, new {len(new)}'

    m = old.merge(new, on=KEY, suffixes=('_old', '_new'), validate='one_to_one')
    assert len(m) == 441, 'per-run rows did not join one-to-one'

    print('=' * 100)
    print('0. THE CLICKS AND TEMPLATES ARE THE SAME RUN')
    print('=' * 100)
    same = {c: int((m[f'{c}_old'].astype(str) == m[f'{c}_new'].astype(str)).sum()) for c in ('seed_ann_id', 'base_size', 'tpl_cx', 'tpl_cy', 'tpl_sha1', 'n_gt_mitotic', 'match_radius_px')}
    print(f'identical on all 441 runs: {same}')
    assert all(v == 441 for v in same.values()), 'the two runs are not the same experiment'
    print('-> every difference below is caused by the search-channel blanking alone, not by a different click or template.')

    print('\n' + '=' * 100)
    print('1. n_detections DELTA')
    print('=' * 100)
    m['d_det'] = m['n_detections_new'] - m['n_detections_old']
    print(f"delta distribution over all 441 runs: {m['d_det'].value_counts().sort_index().to_dict()}")
    print('\nper arm:')
    print(m.groupby('condition')['d_det'].value_counts().unstack(fill_value=0).to_string())
    print('\nstratified by the stored pre-D11 n_self_hits (the mechanism: removing the seed peak frees exactly one of the 100 pre-NMS slots):')
    print(pd.crosstab(col(m, 'n_self_hits'), m['d_det'], rownames=['n_self_hits (old)'], colnames=['n_detections delta']).to_string())
    out_of_family = m[~m['d_det'].isin([0, 1, 2])]
    print(f"\nruns outside the expected +0/+1/+2 family: {len(out_of_family)}")
    if len(out_of_family):
        print(out_of_family[KEY + [c for c in ('n_self_hits', 'n_self_hits_old') if c in out_of_family] + [ 'n_peaks_old', 'n_peaks_new', 'n_detections_old', 'n_detections_new', 'd_det']].to_string(index=False))
    print(f"\nn_peaks: old {m['n_peaks_old'].value_counts().sort_index().to_dict()}, new {m['n_peaks_new'].value_counts().sort_index().to_dict()}")
    print(f"max_peaks_binding new: {int(m['max_peaks_binding_new'].sum())}/441; n_detections new min {new['n_detections'].min()}, max {new['n_detections'].max()}")
    print(f"budget_delivered (n_detections >= 30) new: {int((new['n_detections'] >= 30).sum())}/441")

    print('\n' + '=' * 100)
    print('2. MEAN PRECISION@K PER ARM, OLD vs NEW')
    print('=' * 100)
    rows = []
    for subset in SUBSETS:
        mo, mn = subset_frame(old, subset), subset_frame(new, subset)
        for condition in CONDITIONS:
            co, cn = mo[mo['condition'] == condition], mn[mn['condition'] == condition]
            for k in BUDGETS:
                p_old = co[f'tp_at_{k}'].sum() / (k * len(co))
                p_new = cn[f'tp_at_{k}'].sum() / (k * len(cn))
                rows.append(dict(subset=subset, condition=condition, K=k, n_clicks=len(co), tp_old=int(co[f'tp_at_{k}'].sum()), tp_new=int(cn[f'tp_at_{k}'].sum()), precision_old=p_old, precision_new=p_new, delta_pts=100 * (p_new - p_old)))
    prec = pd.DataFrame(rows)
    for subset in SUBSETS:
        print(f'\n{subset}: pooled precision@K (= mean precision@K over that subset\'s clicks; every click delivers K)')
        print(prec[prec['subset'] == subset][['condition', 'K', 'n_clicks', 'tp_old', 'tp_new', 'precision_old', 'precision_new', 'delta_pts']].to_string(index=False, float_format=lambda v: f'{v:.4f}'))
    print(f"\nlargest |delta| in percentage points, any subset/arm/K: {prec['delta_pts'].abs().max():.4f}")
    print(f"signs of the 9 all_49 arm x K deltas: {prec[prec['subset'] == 'all_49']['delta_pts'].apply(np.sign).value_counts().to_dict()}  (a systematic one-directional move across all three arms would show as a single sign)")

    print('\n' + '=' * 100)
    print('3. EVERY (ROI, seed, condition, K) CELL WHOSE TP COUNT CHANGED')
    print('=' * 100)
    old_tp = old.set_index(KEY)[[f'tp_at_{k}' for k in BUDGETS]]
    new_tp = new.set_index(KEY)[[f'tp_at_{k}' for k in BUDGETS]]
    changed = []
    for k in BUDGETS:
        d = (new_tp[f'tp_at_{k}'] - old_tp[f'tp_at_{k}'])
        for key, v in d[d != 0].items():
            changed.append((key, k, int(v)))
    print(f'cells with a TP@K change: {len(changed)} of {441 * len(BUDGETS)} checked')
    if changed:
        by_delta = pd.Series([c[2] for c in changed]).value_counts().sort_index().to_dict()
        print(f'TP@K deltas: {by_delta}')
        print(f"by arm: {pd.Series([c[0][2] for c in changed]).value_counts().to_dict()}")
        print(f"by K:   {pd.Series([c[1] for c in changed]).value_counts().sort_index().to_dict()}")
    old_mem = {k: topk_membership(old_top, k) for k in BUDGETS}
    new_mem = {k: topk_membership(new_top, k) for k in BUDGETS}
    detail = []
    for key, k, dtp in changed:
        o, n = old_mem[k][key], new_mem[k][key]
        entered = sorted(set(n) - set(o), key=lambda p: n[p][0])
        left = sorted(set(o) - set(n), key=lambda p: o[p][0])
        detail.append(dict(file_name=key[0], seed_index=key[1], condition=key[2], K=k, tp_delta=dtp, n_entered=len(entered), n_left=len(left),
                           entered='; '.join(f'rank{n[p][0]}@({p[0]:.0f},{p[1]:.0f}) {n[p][1]} ann={n[p][2]}' for p in entered),
                           left='; '.join(f'rank{o[p][0]}@({p[0]:.0f},{p[1]:.0f}) {o[p][1]} ann={o[p][2]}' for p in left)))
    DETAIL = pd.DataFrame(detail)
    if len(DETAIL):
        DETAIL.to_csv(f'{NEW_STEM}_tp_changes.csv', index=False)
        for _, r in DETAIL.iterrows():
            print(f"\n{r['file_name']} s{r['seed_index']} {r['condition']} K={r['K']}  TP {r['tp_delta']:+d}")
            print(f"    entered top {r['K']}: {r['entered'] or '(none)'}")
            print(f"    left top {r['K']}:    {r['left'] or '(none)'}")
        print(f"\n-> {NEW_STEM}_tp_changes.csv")

    print('\n' + '=' * 100)
    print('4. DOES THE THREE-WAY RANKING SURVIVE?')
    print('=' * 100)
    roi_order = list(dict.fromkeys(new['file_name']))
    roi_meta = new.drop_duplicates('file_name').set_index('file_name')[['subset', 'domain']]
    rois_by_subset = {s: [fn for fn in roi_order if s == 'all_49' or roi_meta.loc[fn, 'subset'] == s] for s in SUBSETS}
    for n in range(1, 11):
        rng = np.random.default_rng(n)
        d = rng.integers(-6, 7, size=n)
        assert exact_sign_flip_p(d)[0] == brute_force_sign_flip_p(d), f'convolution != enumeration on {d}'
    print('sign-flip convolution validated against brute-force enumeration on random integer vectors.')

    for label, frame in (('OLD (stored reference)', old), ('NEW (post-D11)', new)):
        print(f'\n--- arm ranking, {label} ---')
        piv = []
        for subset in SUBSETS:
            f = subset_frame(frame, subset)
            for k in BUDGETS:
                r = {c: f[f['condition'] == c][f'tp_at_{k}'].sum() / (k * len(f[f['condition'] == c])) for c in CONDITIONS}
                piv.append(dict(subset=subset, K=k, **{c: round(r[c], 4) for c in CONDITIONS}, best=max(r, key=r.get), worst=min(r, key=r.get), order=' > '.join(sorted(r, key=r.get, reverse=True))))
        print(pd.DataFrame(piv).to_string(index=False))

    stats_old, _ = delta_stats(old, roi_order, roi_meta, rois_by_subset)
    stats_new, _ = delta_stats(new, roi_order, roi_meta, rois_by_subset)
    cmp = stats_old.merge(stats_new, on=['subset', 'pair', 'K'], suffixes=('_old', '_new'))
    cmp['delta_pts_old'] = (100 * cmp['pooled_delta_precision_old']).round(3)
    cmp['delta_pts_new'] = (100 * cmp['pooled_delta_precision_new']).round(3)
    cmp['sign_old'] = np.sign(cmp['pooled_delta_precision_old']).astype(int)
    cmp['sign_new'] = np.sign(cmp['pooled_delta_precision_new']).astype(int)
    cmp['sign_kept'] = cmp['sign_old'] == cmp['sign_new']
    cmp['wlt_old'] = cmp['wins_old'].astype(str) + '/' + cmp['losses_old'].astype(str) + '/' + cmp['ties_old'].astype(str)
    cmp['wlt_new'] = cmp['wins_new'].astype(str) + '/' + cmp['losses_new'].astype(str) + '/' + cmp['ties_new'].astype(str)
    cmp['ci_old'] = [f'[{100 * lo:+.2f}, {100 * hi:+.2f}]' for lo, hi in zip(cmp['ci_low_old'], cmp['ci_high_old'])]
    cmp['ci_new'] = [f'[{100 * lo:+.2f}, {100 * hi:+.2f}]' for lo, hi in zip(cmp['ci_low_new'], cmp['ci_high_new'])]
    cmp['d5_kept'] = cmp['d5_bar_met_old'] == cmp['d5_bar_met_new']
    cols = ['pair', 'K', 'delta_pts_old', 'delta_pts_new', 'sign_kept', 'ci_old', 'ci_new', 'wlt_old', 'wlt_new', 'exact_p_old', 'exact_p_new', 'holm_p_old', 'holm_p_new', 'd5_bar_met_old', 'd5_bar_met_new', 'd5_kept']
    for subset in SUBSETS:
        print(f'\n--- paired ROI-level deltas, {subset} ---')
        print(cmp[cmp['subset'] == subset][cols].to_string(index=False, float_format=lambda v: f'{v:.5g}'))
    cmp.to_csv(f'{NEW_STEM}_delta_stats_comparison.csv', index=False)
    print(f'\n-> {NEW_STEM}_delta_stats_comparison.csv')
    print(f"\nsign of the pooled delta preserved on {int(cmp['sign_kept'].sum())}/{len(cmp)} (subset, pair, K) cells")
    print(f"D5 CI-and-majority verdict preserved on {int(cmp['d5_kept'].sum())}/{len(cmp)}")
    flipped = cmp[~cmp['sign_kept'] | ~cmp['d5_kept']]
    if len(flipped):
        print('\ncells where a sign or a D5 verdict moved:')
        print(flipped[['subset'] + cols].to_string(index=False, float_format=lambda v: f'{v:.5g}'))

    print('\n' + '=' * 100)
    print('5. THE ARM-DEPENDENT BLANK FOOTPRINT (a confound D11 introduces, not present pre-D11)')
    print('=' * 100)
    print('Pre-D11 the removed region was a 5.0 px disc, the same in all three arms. Post-D11 it is')
    print('base_size x base_size, and base_size differs by arm, so default_51 blanks the most pixels.')
    print(new.groupby('condition')['n_blanked_px'].agg(['min', 'median', 'max']).to_string())
    print(f"\nsmallest click-to-other-annotation distance anywhere in the 441 runs: {old['nearest_other_ann_px'].min():.1f} px")
    print("default_51's blank square reaches 25 px along each axis and 35.4 px at its corners, so a")
    print('neighbouring annotation can in principle be blanked in default_51 but not in a tightened arm.')
    at_risk = m[col(m, 'nearest_other_ann_px') <= 40]
    print(f'\nruns with another annotation within 40 px of the click: {len(at_risk)}')
    if len(at_risk):
        print(at_risk[KEY + ['base_size_old', 'nearest_other_ann_px_old' if 'nearest_other_ann_px_old' in at_risk else 'nearest_other_ann_px', 'n_blanked_px', 'n_near_click_old', 'n_near_click_new', 'd_det'] + [f'tp_at_{k}_old' for k in BUDGETS] + [f'tp_at_{k}_new' for k in BUDGETS]].to_string(index=False))
    print(f"\ncandidates inside one match radius of the click ('seed annulus'): old {int((old['n_near_click'] > 0).sum())} runs, new {int((new['n_near_click'] > 0).sum())} runs")
    leak = new[new['n_near_click'] > 0]
    if len(leak):
        print(leak[KEY + ['base_size', 'nearest_other_ann_px', 'n_near_click', 'near_click_detail']].to_string(index=False))


if __name__ == '__main__':
    main()
