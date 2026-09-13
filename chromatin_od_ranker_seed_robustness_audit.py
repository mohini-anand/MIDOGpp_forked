"""Independent audit of threshold_maxpeaks_ablation/chromatin_od_ranker_seed_robustness.ipynb.

Run with the anaconda interpreter from the repo root:

    /Users/mohinianand/anaconda3/bin/python3 chromatin_od_ranker_seed_robustness_audit.py

Regenerates every table in Research Logs/2026-09-12-chromatin-od-ranker-seed-robustness-audit.md.
Writes results/chromatin_od_ranker_seed_robustness_audit_*.csv. Reads only; never modifies the
notebook, its CSVs, DECISIONS.md, any existing Research Log, or any midog_utils module.

Tiers, per the auditor protocol:
  Tier A  -- recomputation from the persisted artifacts (the notebook's four CSVs, the committed
             seed-0 oracle CSVs, and databases/MIDOG++.json). Written in plain numpy/pandas/scipy;
             midog_utils helpers are read for their math but are NOT used as the oracle for any
             statistic, with the single deliberate exception of the Tier B pixel recompute (which
             must exercise the production code path to be a check on it at all).
  Tier B  -- two pixel-level recomputes (301.tiff/seed1, and 013.tiff/seed0 as an oracle
             cross-check) plus a cv2.copyMakeBorder microbenchmark. ~40 s of compute.
"""

import itertools
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

REPO = os.path.dirname(os.path.abspath(__file__))
NBDIR = os.path.join(REPO, 'threshold_maxpeaks_ablation')
OUT = os.path.join(REPO, 'results')
PRE = 'chromatin_od_ranker_seed_robustness_audit'

NB = os.path.join(NBDIR, 'chromatin_od_ranker_seed_robustness.ipynb')
TIMING_CSV = os.path.join(NBDIR, 'chromatin_od_ranker_seed_robustness_timing.csv')
PRECISION_CSV = os.path.join(NBDIR, 'chromatin_od_ranker_seed_robustness_precision.csv')
POOLED_CSV = os.path.join(NBDIR, 'chromatin_od_ranker_seed_robustness_per_seed_pooled.csv')
STABILITY_CSV = os.path.join(NBDIR, 'chromatin_od_ranker_seed_robustness_stability.csv')
ORACLE_TIMING = os.path.join(NBDIR, 'chromatin_od_ranker_timing.csv')
ORACLE_PRECISION = os.path.join(NBDIR, 'chromatin_od_ranker_precision.csv')

BUDGETS = (10, 20, 30)
SEEDS = (0, 1, 2, 3, 4)


def banner(s):
    print('\n' + '=' * 78)
    print(s)
    print('=' * 78)


def save(df, name):
    p = os.path.join(OUT, f'{PRE}_{name}.csv')
    df.to_csv(p, index=False)
    print(f'  -> results/{PRE}_{name}.csv  ({df.shape[0]}x{df.shape[1]})')


def sh(cmd):
    return subprocess.run(cmd, shell=True, cwd=REPO, capture_output=True, text=True).stdout.strip()


# =====================================================================================
# GATE 1 -- execution coherence
# =====================================================================================
def gate_execution():
    banner('GATE 1 -- execution coherence')
    nb = json.load(open(NB))
    code = [c for c in nb['cells'] if c['cell_type'] == 'code']
    ecs = [c.get('execution_count') for c in code]
    contiguous = ecs == list(range(1, len(ecs) + 1))
    errs = [i for i, c in enumerate(nb['cells'])
            for o in c.get('outputs', []) if o.get('output_type') == 'error']
    unrun = [i for i, c in enumerate(nb['cells'])
             if c['cell_type'] == 'code' and c.get('execution_count') is None]
    pngs = sum(1 for c in nb['cells'] for o in c.get('outputs', [])
               if 'image/png' in o.get('data', {}))
    rows = [
        dict(check='execution_count contiguous from 1', value=str(ecs), passed=contiguous),
        dict(check='n code cells', value=len(code), passed=True),
        dict(check='error outputs', value=len(errs), passed=len(errs) == 0),
        dict(check='unrun code cells', value=len(unrun), passed=len(unrun) == 0),
        dict(check='last cell run', value=code[-1].get('execution_count'),
             passed=code[-1].get('execution_count') is not None),
        dict(check='embedded PNGs', value=pngs, passed=True),
    ]
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))
    save(df, 'gate_execution')
    return df


# =====================================================================================
# GATE 2 -- provenance
# =====================================================================================
def gate_provenance():
    banner('GATE 2 -- provenance')
    mods = ['midog_utils/' + m for m in
            ['channels.py', 'chromatin.py', 'compare.py', 'dataset.py', 'evaluate.py',
             'find_and_suppress.py', 'invariants.py', 'nms.py', 'seed_selection.py',
             'template_match.py']]
    arts = ['threshold_maxpeaks_ablation/chromatin_od_ranker_seed_robustness.ipynb',
            'threshold_maxpeaks_ablation/chromatin_od_ranker_seed_robustness_timing.csv',
            'threshold_maxpeaks_ablation/chromatin_od_ranker_seed_robustness_precision.csv',
            'threshold_maxpeaks_ablation/chromatin_od_ranker_seed_robustness_per_seed_pooled.csv',
            'threshold_maxpeaks_ablation/chromatin_od_ranker_seed_robustness_stability.csv',
            'threshold_maxpeaks_ablation/chromatin_od_ranker_timing.csv',
            'threshold_maxpeaks_ablation/chromatin_od_ranker_precision.csv',
            'databases/MIDOG++.json']
    rows = []
    for f in mods + arts:
        p = os.path.join(REPO, f)
        mt = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(p)))
        last = sh(f"git log -1 --format='%h %ad %an' --date=short -- '{f}'")
        porcelain = sh(f"git status --porcelain -- '{f}'")
        tracked = sh(f"git ls-files --error-unmatch '{f}' 2>/dev/null") != ''
        rows.append(dict(path=f, kind='module' if f.startswith('midog_utils') else 'artifact',
                         mtime=mt, tracked=tracked, last_commit=last or '(none)',
                         status=porcelain or '(clean)'))
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))
    newest_mod = max(os.path.getmtime(os.path.join(REPO, m)) for m in mods)
    art_mt = min(os.path.getmtime(os.path.join(REPO, a)) for a in arts[1:5])
    print(f'\nNewest module mtime: {time.ctime(newest_mod)}')
    print(f'Oldest notebook artifact mtime: {time.ctime(art_mt)}')
    print(f'All artifacts newer than every module: {art_mt > newest_mod}')
    dirty = [r for r in rows if r['kind'] == 'module' and r['status'] != '(clean)']
    print(f'Uncommitted module modifications: {len(dirty)}')
    save(df, 'provenance')
    return df


# =====================================================================================
# GATE 3 -- composition
# =====================================================================================
def gate_composition():
    banner('GATE 3 -- composition')
    T = pd.read_csv(TIMING_CSV)
    P = pd.read_csv(PRECISION_CSV)
    PS = pd.read_csv(POOLED_CSV)
    S = pd.read_csv(STABILITY_CSV)

    imgs = pd.DataFrame(json.load(open(os.path.join(REPO, 'databases/MIDOG++.json')))['images'])
    imgs['file_name'] = imgs['file_name'].astype(str)
    roi_files = sorted(f for f in os.listdir(os.path.join(REPO, 'images/extra_valid'))
                       if f.endswith('.tiff'))
    extra = imgs[imgs['file_name'].isin(roi_files)]
    per_type = extra.groupby('tumor_type').size().rename('n_rois_in_extra_valid')
    all_per_type = imgs.groupby('tumor_type').size().rename('n_rois_in_full_dataset')
    strata = pd.concat([per_type, all_per_type], axis=1).reset_index()
    print('ROIs per tumor_type (the stratum), measured from databases/MIDOG++.json:')
    print(strata.to_string(index=False))

    rows = []

    def chk(name, got, want, note=''):
        rows.append(dict(check=name, observed=got, expected=want, passed=(got == want), note=note))

    chk('timing rows', len(T), 5 * 14, '5 seeds x 14 ROIs')
    chk('precision rows', len(P), 5 * 14 * 2 * 3, '5 seeds x 14 ROIs x 2 arms x 3 budgets')
    chk('per_seed_pooled rows', len(PS), 5 * 3 * 2, '5 seeds x 3 budgets x 2 arms')
    chk('stability rows', len(S), 3, '3 budgets')
    chk('distinct ROIs in timing', T.file_name.nunique(), 14)
    chk('distinct ROIs in precision', P.file_name.nunique(), 14)
    chk('distinct seeds', T.seed_index.nunique(), 5)
    chk('distinct arms', P.arm.nunique(), 2)
    chk('distinct budgets', P.budget.nunique(), 3)
    chk('tumor_types in extra_valid', extra.tumor_type.nunique(), 7)
    chk('ROIs per tumor_type (min)', int(per_type.min()), 2)
    chk('ROIs per tumor_type (max)', int(per_type.max()), 2)
    chk('duplicate (seed,file) keys in timing', int(T.duplicated(['seed_index', 'file_name']).sum()), 0)
    chk('duplicate (seed,file,arm,budget) keys in precision',
        int(P.duplicated(['seed_index', 'file_name', 'arm', 'budget']).sum()), 0)
    chk('NaN cells in precision', int(P.isna().sum().sum()), 0)
    chk('NaN cells in timing', int(T.isna().sum().sum()), 0)
    chk('rows where budget_delivered != budget', int((P.budget_delivered != P.budget).sum()), 0)
    chk('rows where n_peaks != 100', int((T.n_peaks != 100).sum()), 0)
    chk('rows where chromatin_od_nan_rate != 0', int((T.chromatin_od_nan_rate != 0).sum()), 0)
    # every ROI present at every (seed, arm, budget)
    full = len(list(itertools.product(roi_files, SEEDS, ['tm_score', 'chromatin_od'], BUDGETS)))
    present = P.groupby(['file_name', 'seed_index', 'arm', 'budget']).size()
    chk('complete (roi,seed,arm,budget) grid', len(present), full)
    # n_gt_mitotic must agree with the annotation DB minus the one seed annotation
    anns = pd.DataFrame(json.load(open(os.path.join(REPO, 'databases/MIDOG++.json')))['annotations'])
    fmap = imgs.set_index('id')['file_name'] if 'id' in imgs.columns else None
    df = pd.DataFrame(rows)
    print()
    print(df.to_string(index=False))
    save(df, 'composition')
    save(strata.rename(columns={'index': 'tumor_type'}), 'strata')
    return T, P, PS, S, strata


# =====================================================================================
# CLAIM 1 -- seed_index=0 reproduces the committed seed-0 oracle, EVERY shared column
# =====================================================================================
def claim1(T, P):
    banner('CLAIM 1 -- seed_index=0 vs the committed seed-0 oracle (all shared columns)')
    ot = pd.read_csv(ORACLE_TIMING)
    op = pd.read_csv(ORACLE_PRECISION)
    # the oracle timing file renames two columns; map them onto this notebook's names
    ren = {'t6c_variant_rank_s': 't6c_chromatin_od_rank_s',
           't6c_variant_rank_ms': 't6c_chromatin_od_rank_ms',
           't_variant_pipeline_s': 't_chromatin_od_pipeline_s',
           't_variant_pipeline_ms': 't_chromatin_od_pipeline_ms'}
    ot = ot.rename(columns=ren)

    t0 = T[T.seed_index == 0].set_index('file_name').sort_index()
    o0 = ot.set_index('file_name').sort_index()

    # Timing file: which columns are deterministic (config/result) vs wall-clock?
    TIME_COLS = [c for c in o0.columns if c.startswith('t_') or c.startswith('t1_')
                 or c.startswith('t2_') or c.startswith('t3_') or c.startswith('t4_')
                 or c.startswith('t5_') or c.startswith('t6')]
    shared = [c for c in o0.columns if c in t0.columns]
    det_cols = [c for c in shared if c not in TIME_COLS]
    print(f'Oracle timing columns: {len(ot.columns)}; shared with this run: {len(shared)}; '
          f'deterministic (non-wall-clock): {len(det_cols)}')
    print(f'deterministic columns compared: {det_cols}')

    rows = []
    n_cmp = n_div = 0
    for c in det_cols:
        a, b = t0[c], o0[c]
        if a.dtype.kind in 'fc':
            bad = ~np.isclose(a.astype(float), b.astype(float), rtol=0, atol=1e-5)
        else:
            bad = a.astype(str).values != b.astype(str).values
        n_cmp += len(a)
        n_div += int(np.sum(bad))
        rows.append(dict(source='timing', column=c, n_compared=len(a),
                         n_divergent=int(np.sum(bad)),
                         divergent_rois=','.join(t0.index[np.asarray(bad)].tolist())))

    # Precision file: every shared column, joined on (file_name, arm, budget)
    p0 = P[P.seed_index == 0].copy()
    key = ['file_name', 'arm', 'budget']
    m = p0.merge(op, on=key, suffixes=('_this', '_oracle'), how='outer', indicator=True)
    print(f'\nprecision join: {len(p0)} this x {len(op)} oracle -> {len(m)} rows; '
          f'merge indicator = {dict(m._merge.value_counts())}')
    pcols = [c for c in p0.columns if c in op.columns and c not in key]
    print(f'precision columns compared: {pcols}')
    for c in pcols:
        a, b = m[f'{c}_this'], m[f'{c}_oracle']
        if a.dtype.kind in 'fc':
            bad = ~np.isclose(a.astype(float), b.astype(float), rtol=0, atol=1e-9)
        else:
            bad = a.astype(str).values != b.astype(str).values
        n_cmp += len(a)
        n_div += int(np.sum(bad))
        rows.append(dict(source='precision', column=c, n_compared=len(a),
                         n_divergent=int(np.sum(bad)),
                         divergent_rois=','.join(
                             (m.file_name + '/' + m.arm + '/K' + m.budget.astype(str))[
                                 np.asarray(bad)].tolist())))
    df = pd.DataFrame(rows)
    print()
    print(df.to_string(index=False))
    print(f'\nTOTAL (deterministic columns only): {n_cmp} values compared, {n_div} divergences')
    print(f"Notebook's own claim: 84 values compared, 0 divergences "
          f"(5 timing cols x 14 + 2 arms x 3 budgets x 14 = 70 + 84 = it counted only the tp block)")
    save(df, 'claim1_column_by_column')

    # what the notebook actually checked
    nb_cols = ['seed_ann_id', 'base_size', 'n_detections', 'map_median', 'mad_scale']
    print(f'\nColumns the notebook compared on timing: {nb_cols} ({len(nb_cols)*14} values)')
    print(f'Deterministic timing columns it did NOT compare: '
          f'{[c for c in det_cols if c not in nb_cols]}')
    print(f'Precision columns it compared: [tp_at_budget] ({2*3*14} values)')
    print(f'Precision columns it did NOT compare: {[c for c in pcols if c != "tp_at_budget"]}')
    return df, n_cmp, n_div


# =====================================================================================
# CLAIM 2 -- MAX_PEAKS=100 binds on all 70 draws
# =====================================================================================
def claim2(T):
    banner('CLAIM 2 -- MAX_PEAKS=100 binds on all 70 (seed, ROI) draws')
    print(f'rows with n_peaks == 100 : {(T.n_peaks == 100).sum()}/{len(T)}')
    print(f'n_peaks value counts     : {dict(T.n_peaks.value_counts())}')
    print(f'n_detections (post-NMS)  : min={T.n_detections.min()} max={T.n_detections.max()} '
          f'median={T.n_detections.median()}')
    print(f'rows where n_detections == 100 (would collide with the cap value): '
          f'{(T.n_detections == 100).sum()}')
    print('\nNOTE: n_peaks in this CSV is written by the same run that asserted it, so this is a '
          'consistency check, not an independent one. The independent leg is the Tier B pixel '
          'recompute below, which counts the peaks BEFORE the cap is applied.')
    df = pd.DataFrame([dict(check='n_peaks == MAX_PEAKS on every draw',
                            n_rows=len(T), n_passing=int((T.n_peaks == 100).sum()),
                            n_detections_min=int(T.n_detections.min()),
                            n_detections_max=int(T.n_detections.max()))])
    save(df, 'claim2_cap_binds')
    return df


# =====================================================================================
# CLAIM 3 -- precision delta, per seed and across seeds; plus the ROI-clustered version
# =====================================================================================
def pooled_precision(P, seed, k, arm):
    s = P[(P.seed_index == seed) & (P.budget == k) & (P.arm == arm)]
    return s.tp_at_budget.sum() / s.budget_delivered.sum()


def exact_signflip(x):
    """Exact two-sided sign-flip p-value for H0: mean(x)=0, enumerating all 2^G patterns."""
    x = np.asarray(x, dtype=float)
    G = len(x)
    obs = abs(x.mean())
    signs = np.array(list(itertools.product([-1, 1], repeat=G)), dtype=float)
    means = np.abs(signs @ x) / G
    return float((means >= obs - 1e-12).mean()), G


def claim3(P, PS):
    banner('CLAIM 3 -- precision change: per-seed delta, n=5 CI, and the ROI-clustered version')

    # --- 3a. per-seed pooled precision, recomputed from the long precision CSV
    rows = []
    for s in SEEDS:
        for k in BUDGETS:
            for arm in ['tm_score', 'chromatin_od']:
                mine = pooled_precision(P, s, k, arm)
                theirs = PS[(PS.seed_index == s) & (PS.budget == k)
                            & (PS.arm == arm)].pooled_precision.iloc[0]
                rows.append(dict(seed_index=s, budget=k, arm=arm,
                                 pooled_recomputed=mine, pooled_notebook=theirs,
                                 abs_diff=abs(mine - theirs)))
    pp = pd.DataFrame(rows)
    print(f'per-seed pooled precision: {len(pp)} values compared, '
          f'{(pp.abs_diff > 1e-9).sum()} divergences (max |diff| = {pp.abs_diff.max():.2e})')
    save(pp, 'pooled_precision_recompute')

    # pooled == plain ROI mean, because budget_delivered == K on every row
    piv = pp.pivot_table(index=['seed_index', 'budget'], columns='arm', values='pooled_recomputed')
    mean_of_roi = (P.groupby(['seed_index', 'budget', 'arm']).precision_at_budget.mean()
                   .unstack('arm'))
    print(f'pooled == unweighted ROI mean on all {len(piv)*2} cells: '
          f'{np.allclose(piv.values, mean_of_roi.loc[piv.index].values)}')

    # --- 3b. the notebook's estimator: 5 seeds as the unit
    seed_rows = []
    for k in BUDGETS:
        d = np.array([pooled_precision(P, s, k, 'chromatin_od')
                      - pooled_precision(P, s, k, 'tm_score') for s in SEEDS])
        m, sd = d.mean(), d.std(ddof=1)
        se = sd / np.sqrt(len(d))
        tcrit = stats.t.ppf(0.975, len(d) - 1)
        t_stat, t_p = stats.ttest_1samp(d, 0.0)
        p_sf, G = exact_signflip(d)
        seed_rows.append(dict(budget=k, unit='seed', G=len(d),
                              per_seed_deltas=np.round(d, 4).tolist(),
                              mean=m, ci_lo=m - tcrit * se, ci_hi=m + tcrit * se,
                              t_p=float(t_p), exact_signflip_p=p_sf,
                              signflip_floor=1 / 2 ** (len(d) - 1),
                              n_positive=int((d > 0).sum()), n_units=len(d)))
    sd_df = pd.DataFrame(seed_rows)
    print('\n--- 3b. seed as the unit of replication (the notebook\'s estimator, n=5) ---')
    print(sd_df[['budget', 'G', 'mean', 'ci_lo', 'ci_hi', 't_p', 'exact_signflip_p',
                 'n_positive']].to_string(index=False))
    print('notebook reported:  K=10 mean=+0.0957 CI=[+0.0128,+0.1786]; '
          'K=20 mean=+0.0857 CI=[+0.0291,+0.1423]; K=30 mean=+0.0729 CI=[+0.0287,+0.1170]')
    for _, r in sd_df.iterrows():
        print(f"  K={r.budget:<3} per-seed deltas = {r.per_seed_deltas}")

    # --- 3c. D5's own estimator: the ROI is the exchangeable unit (D5 / F5 s8)
    roi_rows = []
    for metric in ['precision_at_budget', 'recall_at_budget']:
        for k in BUDGETS:
            sub = P[P.budget == k].pivot_table(index=['file_name', 'seed_index'],
                                               columns='arm', values=metric)
            delta = (sub['chromatin_od'] - sub['tm_score']).rename('delta').reset_index()
            per_roi = delta.groupby('file_name').delta.mean()
            x = per_roi.to_numpy()
            m, sd = x.mean(), x.std(ddof=1)
            se = sd / np.sqrt(len(x))
            tcrit = stats.t.ppf(0.975, len(x) - 1)
            t_stat, t_p = stats.ttest_1samp(x, 0.0)
            p_sf, G = exact_signflip(x)
            roi_rows.append(dict(metric=metric, budget=k, unit='ROI', G=len(x),
                                 mean=m, ci_lo=m - tcrit * se, ci_hi=m + tcrit * se,
                                 t_p=float(t_p), exact_signflip_p=p_sf,
                                 signflip_floor=1 / 2 ** (len(x) - 1),
                                 n_positive=int((x > 0).sum()), n_units=len(x),
                                 majority_positive=bool((x > 0).sum() > len(x) / 2)))
    rd_df = pd.DataFrame(roi_rows)
    print('\n--- 3c. ROI as the exchangeable unit (D5 / F5 s8; G=14, exact sign-flip) ---')
    print(rd_df.to_string(index=False))

    both = pd.concat([sd_df.drop(columns=['per_seed_deltas']).assign(metric='precision_at_budget'),
                      rd_df], ignore_index=True, sort=False)
    save(both, 'effect_size_by_unit')
    save(sd_df.assign(per_seed_deltas=sd_df.per_seed_deltas.astype(str)), 'per_seed_delta')

    # --- 3d. per-ROI mean deltas, so the majority-positive criterion is visible
    tab = []
    for k in BUDGETS:
        sub = P[P.budget == k].pivot_table(index=['file_name', 'seed_index'],
                                           columns='arm', values='precision_at_budget')
        d = (sub['chromatin_od'] - sub['tm_score']).groupby('file_name').mean()
        tab.append(d.rename(f'precision_delta_K{k}'))
        sub = P[P.budget == k].pivot_table(index=['file_name', 'seed_index'],
                                           columns='arm', values='recall_at_budget')
        d = (sub['chromatin_od'] - sub['tm_score']).groupby('file_name').mean()
        tab.append(d.rename(f'recall_delta_K{k}'))
    per_roi_tab = pd.concat(tab, axis=1).reset_index()
    dom = P[['file_name', 'tumor_type']].drop_duplicates()
    per_roi_tab = per_roi_tab.merge(dom, on='file_name')
    print('\n--- 3d. per-ROI mean delta over the 5 seeds (chromatin_od - tm_score) ---')
    print(per_roi_tab.to_string(index=False))
    save(per_roi_tab, 'per_roi_delta')

    # --- 3e. worst-seed / worst-ROI behaviour (product metric, Step 4.8)
    ws = []
    for k in BUDGETS:
        d = np.array([pooled_precision(P, s, k, 'chromatin_od')
                      - pooled_precision(P, s, k, 'tm_score') for s in SEEDS])
        sub = P[P.budget == k].pivot_table(index=['file_name', 'seed_index'],
                                           columns='arm', values='precision_at_budget')
        cell = (sub['chromatin_od'] - sub['tm_score'])
        ws.append(dict(budget=k, worst_seed_delta=d.min(), worst_seed=int(SEEDS[int(d.argmin())]),
                       n_seeds_negative=int((d < 0).sum()),
                       worst_cell_delta=cell.min(), worst_cell=str(cell.idxmin()),
                       n_cells_negative=int((cell < 0).sum()), n_cells=len(cell),
                       n_cells_tied=int((cell == 0).sum())))
    ws_df = pd.DataFrame(ws)
    print('\n--- 3e. worst-case behaviour ---')
    print(ws_df.to_string(index=False))
    save(ws_df, 'worst_case')
    return sd_df, rd_df, per_roi_tab, ws_df


# =====================================================================================
# CLAIM 4 -- stability / variance ratio
# =====================================================================================
def claim4(P, S):
    banner('CLAIM 4 -- stability: Levene reproduction and a cluster-aware replacement')
    rows = []
    for k in BUDGETS:
        piv = P[P.budget == k].pivot_table(index=['file_name', 'seed_index'],
                                           columns='arm', values='precision_at_budget')
        resid = piv - piv.groupby('file_name').transform('mean')
        tm_r = resid['tm_score'].dropna().to_numpy()
        ch_r = resid['chromatin_od'].dropna().to_numpy()
        lev_stat, lev_p = stats.levene(tm_r, ch_r)
        v_tm, v_ch = tm_r.var(ddof=1), ch_r.var(ddof=1)
        nbr = S[S.budget == k].iloc[0]
        rows.append(dict(budget=k, var_tm_mine=v_tm, var_ch_mine=v_ch,
                         ratio_mine=v_ch / v_tm, levene_stat_mine=lev_stat, levene_p_mine=lev_p,
                         var_tm_nb=nbr.var_tm_score, var_ch_nb=nbr.var_chromatin_od,
                         ratio_nb=nbr.variance_ratio, levene_p_nb=nbr.levene_p,
                         ratio_absdiff=abs(v_ch / v_tm - nbr.variance_ratio),
                         p_absdiff=abs(lev_p - nbr.levene_p)))
    lv = pd.DataFrame(rows)
    print('--- 4a. Levene, exactly as the notebook computes it (reproduction) ---')
    print(lv[['budget', 'var_tm_mine', 'var_ch_mine', 'ratio_mine', 'levene_p_mine',
              'ratio_nb', 'levene_p_nb', 'ratio_absdiff', 'p_absdiff']].to_string(index=False))
    save(lv, 'claim4_levene_reproduction')

    # --- 4b. the cluster-aware, PAIRED version: per-ROI variance, 14 paired units
    rows = []
    per_roi_var_rows = []
    for k in BUDGETS:
        piv = P[P.budget == k].pivot_table(index=['file_name', 'seed_index'],
                                           columns='arm', values='precision_at_budget')
        v = piv.groupby('file_name').var(ddof=1)
        n_ch_lower = int((v['chromatin_od'] < v['tm_score']).sum())
        n_tie = int((v['chromatin_od'] == v['tm_score']).sum())
        n = len(v) - n_tie
        # exact two-sided sign test on the 14 paired per-ROI variances
        p_sign = float(stats.binomtest(n_ch_lower, n, 0.5).pvalue) if n else np.nan
        # exact sign-flip on the paired log-variance-ratio (permutation of arm labels within ROI)
        eps = 1e-12
        lr = np.log((v['chromatin_od'] + eps) / (v['tm_score'] + eps)).to_numpy()
        p_sf, G = exact_signflip(lr)
        w = stats.wilcoxon(v['chromatin_od'], v['tm_score'])
        rows.append(dict(budget=k, n_rois=len(v),
                         mean_var_tm=v['tm_score'].mean(), mean_var_ch=v['chromatin_od'].mean(),
                         median_var_tm=v['tm_score'].median(),
                         median_var_ch=v['chromatin_od'].median(),
                         ratio_of_mean_var=v['chromatin_od'].mean() / v['tm_score'].mean(),
                         n_rois_ch_lower=n_ch_lower,
                         sign_test_p=p_sign,
                         exact_signflip_p_logratio=p_sf,
                         wilcoxon_p=float(w.pvalue),
                         levene_p_notebook=float(S[S.budget == k].levene_p.iloc[0])))
        vv = v.reset_index().assign(budget=k)
        per_roi_var_rows.append(vv)
    cl = pd.DataFrame(rows)
    print('\n--- 4b. cluster-aware, arm-paired version (G=14 ROIs) ---')
    print(cl[['budget', 'n_rois', 'mean_var_tm', 'mean_var_ch', 'ratio_of_mean_var',
              'n_rois_ch_lower', 'sign_test_p', 'exact_signflip_p_logratio', 'wilcoxon_p',
              'levene_p_notebook']].to_string(index=False))
    save(cl, 'claim4_cluster_aware')
    save(pd.concat(per_roi_var_rows, ignore_index=True), 'per_roi_variance')

    # --- 4c. how many residuals are actually free?
    print('\n--- 4c. degrees of freedom in the Levene sample ---')
    print('Levene is fed 70 residuals per arm, but they are 14 ROI blocks of 5 that each sum to '
          'zero by construction, so the free count is 14 x 4 = 56 per arm, and the two arms are '
          'paired (same ROI, same seed, same candidate pool), which Levene assumes they are not.')
    for k in BUDGETS:
        piv = P[P.budget == k].pivot_table(index=['file_name', 'seed_index'],
                                           columns='arm', values='precision_at_budget')
        resid = piv - piv.groupby('file_name').transform('mean')
        r = stats.pearsonr(resid['tm_score'], resid['chromatin_od'])
        print(f'  K={k}: pearson r between the two arms\' paired residuals = {r[0]:+.4f} '
              f'(p={r[1]:.4f}) -- non-zero pairing means Levene\'s independence assumption fails')
    return lv, cl


# =====================================================================================
# CLAIM 5 -- timing
# =====================================================================================
def claim5(T):
    banner('CLAIM 5 -- timing: is the chromatin_od overhead stable across seeds?')
    T = T.copy()
    T['order'] = np.arange(len(T))
    ov = T.t_chromatin_od_overhead_s * 1000
    print(f'overhead over all {len(T)} draws: mean={ov.mean():.3f}ms median={ov.median():.3f}ms '
          f'std={ov.std():.3f}ms min={ov.min():.3f}ms max={ov.max():.3f}ms')
    print('notebook printed: mean=37.141 median=34.081 std=12.881 min=24.896 max=105.704')

    g = T.groupby('seed_index').agg(
        baseline_ms=('t_baseline_pipeline_s', lambda s: s.mean() * 1000),
        chrom_ms=('t_chromatin_od_pipeline_s', lambda s: s.mean() * 1000),
        overhead_ms=('t_chromatin_od_overhead_s', lambda s: s.mean() * 1000),
        t6a_pad_ms=('t6a_od_pad_s', lambda s: s.mean() * 1000),
        t6b_loop_ms=('t6b_od51_loop_s', lambda s: s.mean() * 1000),
        t6c_rank_ms=('t6c_chromatin_od_rank_s', lambda s: s.mean() * 1000),
        overhead_median_ms=('t_chromatin_od_overhead_s', lambda s: s.median() * 1000),
    )
    g['pct'] = 100 * g.overhead_ms / g.baseline_ms
    print('\n--- 5a. per-seed means, decomposed ---')
    print(g.round(3).to_string())
    print(f'\nspread of per-seed MEAN overhead: {g.overhead_ms.min():.2f} - {g.overhead_ms.max():.2f} ms '
          f'= {g.overhead_ms.max()/g.overhead_ms.min():.2f}x')
    print(f'spread of per-seed MEDIAN overhead: {g.overhead_median_ms.min():.2f} - '
          f'{g.overhead_median_ms.max():.2f} ms = '
          f'{g.overhead_median_ms.max()/g.overhead_median_ms.min():.2f}x')
    print(f'spread of per-seed BASELINE pipeline: {g.baseline_ms.min():.0f} - '
          f'{g.baseline_ms.max():.0f} ms = {g.baseline_ms.max()/g.baseline_ms.min():.2f}x '
          '<-- seed index is run order; this is machine drift, not a seed effect')
    save(g.reset_index(), 'timing_per_seed')

    # --- 5b. is seed 0's high mean carried by t6a pad spikes?
    print('\n--- 5b. pad-spike attribution ---')
    med_pad = T.t6a_od_pad_s.median() * 1000
    T['pad_ms'] = T.t6a_od_pad_s * 1000
    T['pad_over_median'] = T.pad_ms / med_pad
    spikes = T.nlargest(8, 'pad_ms')[['seed_index', 'file_name', 'pad_ms', 'pad_over_median',
                                      't6b_od51_loop_s', 't6c_chromatin_od_rank_s',
                                      't_chromatin_od_overhead_s', 'roi_pixels', 'n_detections']]
    spikes = spikes.assign(loop_ms=spikes.t6b_od51_loop_s * 1000,
                           rank_ms=spikes.t6c_chromatin_od_rank_s * 1000,
                           overhead_ms=spikes.t_chromatin_od_overhead_s * 1000).drop(
        columns=['t6b_od51_loop_s', 't6c_chromatin_od_rank_s', 't_chromatin_od_overhead_s'])
    print(f'median t6a_od_pad = {med_pad:.3f} ms; the 8 largest pad calls:')
    print(spikes.round(3).to_string(index=False))
    save(spikes.round(4), 'pad_spikes')

    # per-seed mean overhead with the two visible outliers removed
    thresh = med_pad * 2
    clean = T[T.pad_ms <= thresh]
    gc_ = clean.groupby('seed_index').t_chromatin_od_overhead_s.mean() * 1000
    print(f'\nper-seed mean overhead after dropping the {int((T.pad_ms > thresh).sum())} draws whose '
          f'pad call exceeded 2x the median pad time:')
    print(gc_.round(3).to_string())
    print(f'spread now: {gc_.min():.2f} - {gc_.max():.2f} ms = {gc_.max()/gc_.min():.2f}x '
          f'(was {g.overhead_ms.max()/g.overhead_ms.min():.2f}x)')

    # --- 5b-bis. cell 19 claims the seed-0 mean overhead "matches" the committed oracle.
    #     Verification A only ever compared DETERMINISTIC columns; timing is wall clock.
    print('\n--- 5b-bis. does the seed-0 timing actually match the committed oracle? ---')
    O = pd.read_csv(ORACLE_TIMING)
    t0 = T[T.seed_index == 0]
    a = O.set_index('file_name').t_chromatin_od_overhead_s.sort_index()
    b = t0.set_index('file_name').t_chromatin_od_overhead_s.sort_index()
    orows = [dict(quantity='mean t_chromatin_od_overhead_ms',
                  oracle=round(float(O.t_chromatin_od_overhead_s.mean() * 1000), 3),
                  this_run=round(float(t0.t_chromatin_od_overhead_s.mean() * 1000), 3)),
             dict(quantity='mean t6a_od_pad_ms',
                  oracle=round(float(O.t6a_od_pad_s.mean() * 1000), 3),
                  this_run=round(float(t0.t6a_od_pad_s.mean() * 1000), 3)),
             dict(quantity='mean t_baseline_pipeline_ms',
                  oracle=round(float(O.t_baseline_pipeline_s.mean() * 1000), 3),
                  this_run=round(float(t0.t_baseline_pipeline_s.mean() * 1000), 3))]
    od = pd.DataFrame(orows)
    od['abs_diff'] = (od.oracle - od.this_run).abs()
    od['pct_diff'] = (100 * od.abs_diff / od.oracle).round(2)
    print(od.to_string(index=False))
    print(f'  per-ROI overhead identical across the two runs: {np.allclose(a, b)}; '
          f'max per-ROI abs diff = {float((a - b).abs().max() * 1000):.3f} ms')
    print('  Cell 19 prints "seed_index=0 mean overhead here: 41.988ms -- matches '
          'chromatin_od_ranker_timing.csv". It does not: the oracle mean is 45.094 ms. And the '
          'parenthetical ("Verification A already confirmed the per-ROI values this is computed '
          'from are identical") cites a determinism check on seed_ann_id / base_size / '
          'n_detections / map_median / mad_scale -- none of which is a timing column.')
    save(od, 'timing_vs_oracle')

    # --- 5b-ter. t6a's share of the overhead
    share = (T.groupby('seed_index').t6a_od_pad_s.mean()
             / T.groupby('seed_index').t_chromatin_od_overhead_s.mean())
    print(f'\nt6a (full-ROI pad) share of the total overhead, per seed: '
          f'{[round(x, 3) for x in share]} -- {share.min():.1%} to {share.max():.1%}')

    # --- 5c. does the overhead track n_detections, as the readout claims?
    print('\n--- 5c. the readout claims the overhead "depends on pool composition (n_detections), '
          'not on which click was drawn" ---')
    corr = []
    for y, yl in [('t_chromatin_od_overhead_s', 'total overhead'),
                  ('t6a_od_pad_s', 't6a pad'), ('t6b_od51_loop_s', 't6b od51 loop'),
                  ('t6c_chromatin_od_rank_s', 't6c rank')]:
        for x in ['n_detections', 'roi_pixels', 'order', 'base_size']:
            r = stats.spearmanr(T[x], T[y])
            corr.append(dict(stage=yl, predictor=x, spearman_rho=r.statistic, p=r.pvalue))
    cdf = pd.DataFrame(corr)
    print(cdf.round(4).to_string(index=False))
    print(f'\nn_detections range across the 70 draws: {T.n_detections.min()}-{T.n_detections.max()} '
          f'(sd={T.n_detections.std():.2f}) -- a 13-count spread, so "depends on n_detections" is '
          'barely testable at this config')
    save(cdf, 'timing_correlations')

    # --- 5d. repeat counts
    print('\n--- 5d. measurement design ---')
    print('Each of the 70 draws is a SINGLE shot: run_roi has no repeat loop, so every per-draw '
          'stage time is n=1. Per-seed "spread" is across 14 different ROIs at 14 different '
          'clicks, not measurement precision.')
    print(f'n_retries across the 70 draws: {dict(T.n_retries.value_counts().sort_index())}')
    return g, cdf


# =====================================================================================
# CLAIM 6 -- seed_annulus_empty, root-cause and arm symmetry (Tier B, from pixels)
# =====================================================================================
def claim6(T, P):
    banner('CLAIM 6 -- seed_annulus_empty on 301.tiff/seed1: root cause and arm symmetry')
    sys.path.insert(0, REPO)
    import cv2
    from midog_utils import channels as ch
    from midog_utils import chromatin as cm
    from midog_utils import dataset as ds
    from midog_utils import evaluate as ev
    from midog_utils import find_and_suppress as fs
    from midog_utils import seed_selection as ss
    from midog_utils import template_match as tm
    from midog_utils.nms import nms_by_distance

    CHANNEL, METHOD = 'hematoxylin_od', cv2.TM_CCOEFF
    PEAK_MIN_DISTANCE, SELF_HIT_RADIUS, DEEP_FLOOR_Z, MAX_PEAKS = 7, 5.0, -1.5, 100
    OD51_WINDOW, OD_PAD_51 = 51, 25
    CFG = fs.FSConfig(channel=CHANNEL, base_size=tm.BASE_SIZE, scales=(1.0,), n_angles=1,
                      flips=(False,), peak_min_distance=PEAK_MIN_DISTANCE,
                      self_hit_radius=SELF_HIT_RADIUS)
    BORDER, OTSU_WINDOW = CFG.patch_size // 2, tm.BASE_SIZE

    images, annotations = ds.load_annotations(os.path.join(REPO, 'databases/MIDOG++.json'))
    meta = images.set_index('file_name')[['image_id', 'tumor_type']]

    def recompute(fn, seed_index, uncapped=False):
        t0 = time.time()
        path = os.path.join(REPO, 'images/extra_valid', fn)
        rgb = ds.load_roi(path)
        mpp = ds.roi_mpp(path)
        roi_shape = rgb.shape
        match_radius = ev.radius_px(mpp, ev.MIDOG_RADIUS_UM)
        nms_radius = ev.radius_px(mpp, ev.MIDOG_RADIUS_UM)
        hem = ch.to_channel(rgb, CHANNEL)
        gray_inv = ch.to_gray_inverted(rgb)
        H, W = hem.shape[:2]
        del rgb
        gt = ds.image_annotations(annotations, fn)
        pool_seeds = ss.agreement_pool(gt[gt['category_id'] == ds.MITOTIC])[0]
        pool_seeds = ss.border_filter(pool_seeds, BORDER, roi_shape)
        image_id = int(meta.loc[fn, 'image_id'])
        rng = np.random.default_rng([seed_index, image_id])
        working, retries = pool_seeds.copy(), 0
        seed = tpl = None
        while len(working) > 0:
            idx = int(rng.integers(len(working)))
            row = working.iloc[idx]
            r = ss.tightened_template_box(gray_inv, float(row['cx']), float(row['cy']),
                                          otsu_window=OTSU_WINDOW)
            if r is not None and tm.read_padded_patch(hem, r[1], r[2], CFG.patch_size) is not None:
                seed, tpl = row, r
                break
            working = working.drop(working.index[idx])
            retries += 1
        base_size, tpl_cx, tpl_cy = tpl
        tpl_xy = (float(tpl_cx), float(tpl_cy))
        seed_ann_id = int(seed['ann_id'])
        gt_eval = gt[gt['ann_id'] != seed_ann_id].reset_index(drop=True)

        patch = tm.read_padded_patch(hem, *tpl_xy, CFG.patch_size)
        templates, _ = tm.build_augmentations(patch, base_size, CFG.scales, CFG.n_angles, CFG.flips)
        PAD = max((t.shape[0] - 1) // 2 for t in templates)
        hem_p = cv2.copyMakeBorder(hem, PAD, PAD, PAD, PAD, borderType=cv2.BORDER_REPLICATE)
        fused_p, _, valid_p = tm.fused_response(hem_p, templates, CFG.scale_normalize, method=METHOD)
        fused = fused_p[PAD:PAD + H, PAD:PAD + W]
        valid = valid_p[PAD:PAD + H, PAD:PAD + W]
        med, mad = tm.robust_stats(fused, valid)
        cut = med + DEEP_FLOOR_Z * mad
        cap = 2_000_000 if uncapped else MAX_PEAKS
        centers, scores = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE, cut, cap)
        n_peaks_uncapped = len(centers)
        centers, scores = centers[:MAX_PEAKS], scores[:MAX_PEAKS]
        keep = nms_by_distance(centers, scores, nms_radius)
        c, s = centers[keep], scores[keep]
        ok = np.hypot(c[:, 0] - tpl_xy[0], c[:, 1] - tpl_xy[1]) > SELF_HIT_RADIUS
        c, s = c[ok], s[ok]
        pool = pd.DataFrame({'cx': c[:, 0], 'cy': c[:, 1], 'score': s})
        hem_pad51 = cv2.copyMakeBorder(hem, OD_PAD_51, OD_PAD_51, OD_PAD_51, OD_PAD_51,
                                       cv2.BORDER_REPLICATE)
        pool['od51'] = [cm.chromatin_density(hem_pad51, x + OD_PAD_51, y + OD_PAD_51,
                                             window=OD51_WINDOW)
                        for x, y in zip(pool.cx.to_numpy(), pool.cy.to_numpy())]
        return dict(pool=pool, gt=gt, gt_eval=gt_eval, tpl_xy=tpl_xy, base_size=base_size,
                    seed_ann_id=seed_ann_id, mpp=mpp, match_radius=match_radius,
                    nms_radius=nms_radius, n_peaks_uncapped=n_peaks_uncapped, med=med, mad=mad,
                    n_retries=retries, elapsed=time.time() - t0, roi_shape=(H, W))

    # --- independent greedy one-to-one matcher, written here, not imported ---
    def greedy_tp(ranked_xy, gt_xy, gt_cls, radius):
        """Return per-rank bucket: 1 = matched a mitosis, 2 = matched a look-alike, 0 = none."""
        taken = np.zeros(len(gt_xy), dtype=bool)
        out = np.zeros(len(ranked_xy), dtype=int)
        for i, (x, y) in enumerate(ranked_xy):
            if len(gt_xy) == 0:
                continue
            d = np.hypot(gt_xy[:, 0] - x, gt_xy[:, 1] - y)
            d[taken] = np.inf
            j = int(np.argmin(d))
            if d[j] <= radius:
                taken[j] = True
                out[i] = int(gt_cls[j])
        return out

    def prec_at(pool, key, gt_eval, radius, drop_idx=None):
        p = pool if drop_idx is None else pool.drop(index=drop_idx)
        r = p.sort_values(key, ascending=False, na_position='last', kind='mergesort')
        b = greedy_tp(r[['cx', 'cy']].to_numpy(), gt_eval[['cx', 'cy']].to_numpy(),
                      gt_eval['category_id'].to_numpy(), radius)
        cum = np.cumsum(b == 1)
        return {k: int(cum[min(k, len(cum)) - 1]) for k in BUDGETS}, r

    rows = []
    # ---- 6a. 301.tiff / seed1 : the failing draw
    R = recompute('301.tiff', 1)
    print(f'recomputed 301.tiff/seed1 from pixels in {R["elapsed"]:.1f}s')
    print(f"  base_size={R['base_size']} (csv {int(T[(T.file_name=='301.tiff')&(T.seed_index==1)].base_size.iloc[0])}) "
          f"seed_ann_id={R['seed_ann_id']} (csv {int(T[(T.file_name=='301.tiff')&(T.seed_index==1)].seed_ann_id.iloc[0])}) "
          f"n_detections={len(R['pool'])} (csv {int(T[(T.file_name=='301.tiff')&(T.seed_index==1)].n_detections.iloc[0])})")
    print(f"  map_median={R['med']:.5f} (csv {T[(T.file_name=='301.tiff')&(T.seed_index==1)].map_median.iloc[0]}) "
          f"mad_scale={R['mad']:.5f} (csv {T[(T.file_name=='301.tiff')&(T.seed_index==1)].mad_scale.iloc[0]})")
    print(f"  n_retries={R['n_retries']} (csv {int(T[(T.file_name=='301.tiff')&(T.seed_index==1)].n_retries.iloc[0])})")

    pool = R['pool']
    d_seed = np.hypot(pool.cx - R['tpl_xy'][0], pool.cy - R['tpl_xy'][1])
    near = pool[d_seed <= R['match_radius']]
    print(f"\n  match_radius = {R['match_radius']:.3f}px (mpp={R['mpp']:.4f}), self_hit = 5.0px")
    print(f"  candidates inside match_radius of the template centre: {len(near)}")
    print(near.assign(dist=d_seed[near.index]).round(3).to_string())

    echo = near.index[0]
    e = pool.loc[echo]
    # where does the echo rank in each arm?
    rank_tm = int((pool.score > e.score).sum())
    rank_od = int((pool.od51 > e.od51).sum())
    print(f"\n  echo candidate at ({e.cx:.0f},{e.cy:.0f}) score={e.score:.4f} od51={e.od51:.4f}")
    print(f"  rank in tm_score ordering    : {rank_tm} (0-based) -> list position {rank_tm+1}/{len(pool)}")
    print(f"  rank in chromatin_od ordering: {rank_od} (0-based) -> list position {rank_od+1}/{len(pool)}")

    # is it a forced FP? nearest GT of either category
    g = R['gt']
    dg = np.hypot(g.cx - e.cx, g.cy - e.cy)
    nearest = g.assign(dist=dg).nsmallest(3, 'dist')[['ann_id', 'category_id', 'cx', 'cy', 'dist']]
    print('\n  nearest GT (either category) to the echo candidate:')
    print(nearest.round(2).to_string(index=False))
    forced_fp = bool(nearest.iloc[0].ann_id == R['seed_ann_id']
                     and (dg <= R['match_radius']).sum() == 1)
    print(f"  nearest GT ann_id == seed_ann_id ({R['seed_ann_id']}): "
          f"{bool(nearest.iloc[0].ann_id == R['seed_ann_id'])}")
    print(f"  n GT of either category within match_radius: {int((dg <= R['match_radius']).sum())}")
    print(f"  -> the echo candidate is a FORCED false positive: {forced_fp}")
    print("     (its only reachable ground truth is the seed annotation, which gt_eval excludes)")

    # --- arm-asymmetry: precision with and without the echo candidate
    tp_with, _ = prec_at(pool, 'score', R['gt_eval'], R['match_radius'])
    tp_wo, _ = prec_at(pool, 'score', R['gt_eval'], R['match_radius'], drop_idx=echo)
    tp_with_od, _ = prec_at(pool, 'od51', R['gt_eval'], R['match_radius'])
    tp_wo_od, _ = prec_at(pool, 'od51', R['gt_eval'], R['match_radius'], drop_idx=echo)
    for k in BUDGETS:
        csv_tm = int(P[(P.file_name == '301.tiff') & (P.seed_index == 1) & (P.budget == k)
                       & (P.arm == 'tm_score')].tp_at_budget.iloc[0])
        csv_od = int(P[(P.file_name == '301.tiff') & (P.seed_index == 1) & (P.budget == k)
                       & (P.arm == 'chromatin_od')].tp_at_budget.iloc[0])
        rows.append(dict(
            roi='301.tiff', seed=1, budget=k,
            tm_tp_csv=csv_tm, tm_tp_recomputed=tp_with[k], tm_tp_without_echo=tp_wo[k],
            od_tp_csv=csv_od, od_tp_recomputed=tp_with_od[k], od_tp_without_echo=tp_wo_od[k],
            tm_prec=tp_with[k] / k, od_prec=tp_with_od[k] / k,
            delta_with_echo=(tp_with_od[k] - tp_with[k]) / k,
            delta_without_echo=(tp_wo_od[k] - tp_wo[k]) / k,
            echo_bias_on_this_cell=((tp_with_od[k] - tp_with[k]) - (tp_wo_od[k] - tp_wo[k])) / k,
            echo_bias_on_seed1_pooled=((tp_with_od[k] - tp_with[k])
                                       - (tp_wo_od[k] - tp_wo[k])) / (k * 14)))
    echo_df = pd.DataFrame(rows)
    print('\n--- 6b. does the echo candidate bias the two arms differently? ---')
    print(echo_df.round(5).to_string(index=False))
    save(echo_df, 'claim6_echo_bias')

    # ---- 6b-bis. Is the "echo" actually SECONDARY? Where is the click's own self-match?
    print("\n--- 6b-bis. the notebook calls the annulus candidate a SECONDARY local maximum "
          "(an \"echo\") of the click's own correlation surface. Is there a primary? ---")
    import cv2 as _cv2
    selfrows = []
    for fn, si in [('301.tiff', 1), ('013.tiff', 0)]:
        path = os.path.join(REPO, 'images/extra_valid', fn)
        rgb = ds.load_roi(path)
        mpp = ds.roi_mpp(path)
        roi_shape = rgb.shape
        mr = ev.radius_px(mpp, ev.MIDOG_RADIUS_UM)
        hem = ch.to_channel(rgb, CHANNEL)
        gray_inv = ch.to_gray_inverted(rgb)
        H, W = hem.shape[:2]
        del rgb
        gt = ds.image_annotations(annotations, fn)
        ps = ss.border_filter(ss.agreement_pool(gt[gt['category_id'] == ds.MITOTIC])[0],
                              BORDER, roi_shape)
        rng = np.random.default_rng([si, int(meta.loc[fn, 'image_id'])])
        w = ps.copy()
        seed = tpl = None
        while len(w):
            i = int(rng.integers(len(w)))
            row = w.iloc[i]
            r = ss.tightened_template_box(gray_inv, float(row['cx']), float(row['cy']),
                                          otsu_window=OTSU_WINDOW)
            if r is not None and tm.read_padded_patch(hem, r[1], r[2], CFG.patch_size) is not None:
                seed, tpl = row, r
                break
            w = w.drop(w.index[i])
        bs, tcx, tcy = tpl
        txy = (float(tcx), float(tcy))
        patch = tm.read_padded_patch(hem, *txy, CFG.patch_size)
        tps, _ = tm.build_augmentations(patch, bs, CFG.scales, CFG.n_angles, CFG.flips)
        PAD = max((t.shape[0] - 1) // 2 for t in tps)
        hp = _cv2.copyMakeBorder(hem, PAD, PAD, PAD, PAD, borderType=_cv2.BORDER_REPLICATE)
        fp_, _, vp_ = tm.fused_response(hp, tps, CFG.scale_normalize, method=METHOD)
        fused = fp_[PAD:PAD + H, PAD:PAD + W]
        valid = vp_[PAD:PAD + H, PAD:PAD + W]
        med, mad = tm.robust_stats(fused, valid)
        C, Sc = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE, med + DEEP_FLOOR_Z * mad, MAX_PEAKS)
        dd = np.hypot(C[:, 0] - txy[0], C[:, 1] - txy[1])
        gy, gx = np.unravel_index(int(np.argmax(fused)), fused.shape)
        # the fused score AT the click's own template-box centre
        self_score = float(fused[int(round(txy[1])), int(round(txy[0]))])
        keep = nms_by_distance(C, Sc, mr)
        c2, s2 = C[keep], Sc[keep]
        d2 = np.hypot(c2[:, 0] - txy[0], c2[:, 1] - txy[1])
        selfrows.append(dict(
            file_name=fn, seed_index=si, base_size=bs, tpl_xy=str(txy),
            match_radius_px=round(mr, 2),
            peak_rank100_score_cut=round(float(Sc.min()), 4),
            top_peak_score=round(float(Sc[0]), 4),
            top_peak_dist_from_click_px=round(float(dd[0]), 2),
            fused_value_at_click=round(self_score, 4),
            click_makes_top100=bool(self_score >= float(Sc.min())),
            n_peaks_within_selfhit_5px=int((dd <= 5).sum()),
            n_peaks_within_match_radius=int((dd <= mr).sum()),
            nms_kept=int(len(keep)),
            n_removed_by_selfhit=int((d2 <= 5.0).sum()),
            n_detections=int((d2 > 5.0).sum())))
        del hem, gray_inv, hp, fp_, vp_, fused, valid
    sdf = pd.DataFrame(selfrows)
    print(sdf.to_string(index=False))
    print("\nReading: on 013.tiff/seed0 the click's own location IS the global maximum "
          "(0.50px away, score 3.0023) and self-hit suppression removes it -- the normal case. "
          "On 301.tiff/seed1 there is NO peak within 5px of the click at all: the fused value at "
          "the click's own template-box centre falls below the rank-100 score cut, so the click's "
          "own location never entered the pool. The annulus candidate is therefore not a "
          "SECONDARY maximum -- there is no primary -- and n_detections=100 on that draw is "
          "explained by NMS removing nothing and self-hit having nothing to remove.")
    save(sdf, 'claim6_self_peak')

    # ---- 6c. the pre-cap peak count (independent evidence the cap genuinely bound)
    print('\n--- 6c. Tier B: does MAX_PEAKS=100 genuinely bind? (uncapped recompute) ---')
    unc = []
    for fn, si in [('301.tiff', 1), ('013.tiff', 0)]:
        Ru = recompute(fn, si, uncapped=True)
        csv = T[(T.file_name == fn) & (T.seed_index == si)]
        unc.append(dict(file_name=fn, seed_index=si,
                        n_peaks_uncapped=Ru['n_peaks_uncapped'],
                        n_peaks_in_csv=int(csv.n_peaks.iloc[0]),
                        cap_binds=Ru['n_peaks_uncapped'] > 100,
                        n_detections_recomputed=len(Ru['pool']),
                        n_detections_in_csv=int(csv.n_detections.iloc[0]),
                        seed_ann_id_recomputed=Ru['seed_ann_id'],
                        seed_ann_id_in_csv=int(csv.seed_ann_id.iloc[0]),
                        base_size_recomputed=Ru['base_size'],
                        base_size_in_csv=int(csv.base_size.iloc[0]),
                        map_median_recomputed=round(float(Ru['med']), 5),
                        map_median_in_csv=float(csv.map_median.iloc[0]),
                        elapsed_s=round(Ru['elapsed'], 1)))
        if fn == '013.tiff':
            # full independent precision cross-check against the committed oracle
            for k in BUDGETS:
                tpm, _ = prec_at(Ru['pool'], 'score', Ru['gt_eval'], Ru['match_radius'])
                tpo, _ = prec_at(Ru['pool'], 'od51', Ru['gt_eval'], Ru['match_radius'])
            print(f'  013.tiff/seed0 independent greedy match: tm tp@{BUDGETS} = '
                  f'{[tpm[k] for k in BUDGETS]}, chromatin_od tp@{BUDGETS} = '
                  f'{[tpo[k] for k in BUDGETS]}')
            for arm, tps in [('tm_score', tpm), ('chromatin_od', tpo)]:
                for k in BUDGETS:
                    c = int(P[(P.file_name == '013.tiff') & (P.seed_index == 0) & (P.budget == k)
                              & (P.arm == arm)].tp_at_budget.iloc[0])
                    print(f'    {arm:14s} K={k:<3} recomputed={tps[k]}  csv={c}  '
                          f'{"MATCH" if tps[k] == c else "DIVERGENCE"}')
    ud = pd.DataFrame(unc)
    print()
    print(ud.to_string(index=False))
    save(ud, 'claim6_tierb_from_pixels')
    return echo_df, ud


# =====================================================================================
# CHECK 7 -- effective units: duplicate clicks across seeds
# =====================================================================================
def check_effective_units(T, P):
    banner('CHECK 7 -- effective units (duplicate clicks shrink observed seed variance)')
    n_distinct = T.groupby('file_name').seed_ann_id.nunique()
    print('distinct clicks per ROI across the 5 seeds:')
    print(n_distinct.to_string())
    print(f'\ntotal distinct (ROI, click) pairs: {n_distinct.sum()} of {len(T)} draws '
          f'-- {len(T) - n_distinct.sum()} draws are exact repeats of another seed on the same ROI')
    dups = []
    for fn, grp in T.groupby('file_name'):
        seen = {}
        for _, r in grp.iterrows():
            seen.setdefault(int(r.seed_ann_id), []).append(int(r.seed_index))
        for ann, ss_ in seen.items():
            if len(ss_) > 1:
                dups.append(dict(file_name=fn, seed_ann_id=ann, seeds=str(ss_), n=len(ss_)))
    dd = pd.DataFrame(dups)
    print('\nduplicated clicks:')
    print(dd.to_string(index=False) if len(dd) else '(none)')
    # how much do duplicates shrink the observed within-ROI variance?
    rows = []
    for k in BUDGETS:
        piv = P[P.budget == k].pivot_table(index=['file_name', 'seed_index'],
                                           columns='arm', values='precision_at_budget')
        v_all = piv.groupby('file_name').var(ddof=1).mean()
        keep = T[['file_name', 'seed_index', 'seed_ann_id']].drop_duplicates(
            ['file_name', 'seed_ann_id'])
        idx = pd.MultiIndex.from_frame(keep[['file_name', 'seed_index']])
        v_uniq = piv.loc[piv.index.intersection(idx)].groupby('file_name').var(ddof=1).mean()
        rows.append(dict(budget=k,
                         mean_var_tm_all70=v_all['tm_score'], mean_var_tm_unique64=v_uniq['tm_score'],
                         mean_var_ch_all70=v_all['chromatin_od'],
                         mean_var_ch_unique64=v_uniq['chromatin_od']))
    ef = pd.DataFrame(rows)
    print('\nmean per-ROI variance, all 70 draws vs. the 64 unique (ROI, click) pairs only:')
    print(ef.round(5).to_string(index=False))
    save(ef, 'effective_units')
    if len(dd):
        save(dd, 'duplicate_clicks')
    return n_distinct, dd, ef


# =====================================================================================
# CHECK 8 -- the figure (step 1a: re-derive the plotted series)
# =====================================================================================
def check_figure(P, PS):
    banner('CHECK 8 -- figure (cell 25): re-derive the six plotted series')
    rows = []
    for k in BUDGETS:
        for arm in ['tm_score', 'chromatin_od']:
            plotted = PS[(PS.budget == k) & (PS.arm == arm)].sort_values(
                'seed_index').pooled_precision.to_numpy()
            mine = np.array([pooled_precision(P, s, k, arm) for s in SEEDS])
            rows.append(dict(budget=k, arm=arm,
                             plotted=np.round(plotted, 4).tolist(),
                             recomputed=np.round(mine, 4).tolist(),
                             max_abs_diff=float(np.max(np.abs(plotted - mine))),
                             ymin=float(mine.min()), ymax=float(mine.max())))
    fd = pd.DataFrame(rows)
    print(fd.to_string(index=False))
    allv = np.concatenate([np.array(r) for r in fd.recomputed])
    print(f'\nplotted value range over all six series: {allv.min():.4f} to {allv.max():.4f} '
          '-- a shared y-axis spanning this range clips nothing')
    print(f'total plotted points: {len(allv)} (= 3 budgets x 2 arms x 5 seeds); '
          'the unit is the seed, and 5 seeds are drawn, so no per-unit over-plotting')
    save(fd.assign(plotted=fd.plotted.astype(str), recomputed=fd.recomputed.astype(str)),
         'figure_series')
    return fd


# =====================================================================================
# CHECK 9 -- config drift vs FSConfig defaults / invariants / DECISIONS
# =====================================================================================
def check_config():
    banner('CHECK 9 -- config drift')
    sys.path.insert(0, REPO)
    import cv2
    from midog_utils import find_and_suppress as fs
    from midog_utils import evaluate as ev
    from midog_utils import template_match as tm
    from midog_utils import chromatin as cm
    d = fs.FSConfig()
    rows = [
        dict(param='channel', notebook='hematoxylin_od',
             fsconfig_default=f'{d.channel} (INERT: ch.to_channel is called directly)',
             decision='D3 (no rescale after deconvolution)', agrees=True),
        dict(param='tm_method',
             notebook=f'cv2.TM_CCOEFF={cv2.TM_CCOEFF}',
             fsconfig_default=f'{getattr(d, "tm_method", "n/a")} '
                              f'(=TM_CCOEFF_NORMED; INERT: method= is passed to fused_response)',
             decision='D1 (moved to DECISIONS_UNVERIFIED.md 2026-09-12) selects TM_CCOEFF',
             agrees=True),
        dict(param='peak_min_distance', notebook=7, fsconfig_default=d.peak_min_distance,
             decision='-', agrees=d.peak_min_distance == 7),
        dict(param='self_hit_radius', notebook=5.0, fsconfig_default=d.self_hit_radius,
             decision='-', agrees=d.self_hit_radius == 5.0),
        dict(param='base_size', notebook='tm.BASE_SIZE', fsconfig_default=d.base_size,
             decision='-', agrees=d.base_size == tm.BASE_SIZE),
        dict(param='max_peaks', notebook=100, fsconfig_default=getattr(d, 'max_peaks', 'n/a'),
             decision='D9 (max_peaks=100)', agrees=True),
        dict(param='NMS_RADIUS_UM', notebook='ev.MIDOG_RADIUS_UM', fsconfig_default='n/a',
             decision=f'D7 (7.5 um) -- module constant is {ev.MIDOG_RADIUS_UM}',
             agrees=ev.MIDOG_RADIUS_UM == 7.5),
        dict(param='MATCH_RADIUS_UM', notebook='ev.MIDOG_RADIUS_UM', fsconfig_default='n/a',
             decision='D7', agrees=ev.MIDOG_RADIUS_UM == 7.5),
        dict(param='od51 window', notebook=51, fsconfig_default=f'tm.BASE_SIZE={tm.BASE_SIZE}',
             decision='D5 (51px "is reading the neighbours")', agrees=tm.BASE_SIZE == 51),
        dict(param='od51 frac', notebook='cm.DEFAULT_FRAC', fsconfig_default=cm.DEFAULT_FRAC,
             decision='-', agrees=True),
        dict(param='DEEP_FLOOR_Z', notebook=-1.5, fsconfig_default='n/a',
             decision='D9 ("DEEP_FLOOR_Z is unchanged at -1.5")', agrees=True),
        dict(param='images dir', notebook='images/extra_valid (14 ROIs)', fsconfig_default='n/a',
             decision='D5 standing constraint: must use the 14 ROIs of images/extra_valid',
             agrees=True),
    ]
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))
    print(f'\nFSConfig defaults not overridden but unused by this notebook '
          f'(find_and_suppress() is not the caller): max_peaks={getattr(d,"max_peaks","?")}, '
          f'max_detections={getattr(d,"max_detections","?")}, '
          f'score_threshold={getattr(d,"score_threshold","?")}')
    save(df, 'config_drift')
    return df


# =====================================================================================
# CHECK 10 -- pad microbenchmark (Tier B): is the pad-noise attribution reasonable?
# =====================================================================================
def check_pad_bench():
    banner('CHECK 10 -- Tier B: cv2.copyMakeBorder repeat-timing microbenchmark')
    sys.path.insert(0, REPO)
    import cv2
    from midog_utils import channels as ch
    from midog_utils import dataset as ds
    rgb = ds.load_roi(os.path.join(REPO, 'images/extra_valid/013.tiff'))
    hem = ch.to_channel(rgb, 'hematoxylin_od')
    del rgb
    ts = []
    for _ in range(40):
        t = time.perf_counter()
        p = cv2.copyMakeBorder(hem, 25, 25, 25, 25, cv2.BORDER_REPLICATE)
        ts.append((time.perf_counter() - t) * 1000)
        del p
    ts = np.array(ts)
    df = pd.DataFrame([dict(roi='013.tiff', shape=str(hem.shape), n_repeats=len(ts),
                            mean_ms=ts.mean(), median_ms=np.median(ts), min_ms=ts.min(),
                            max_ms=ts.max(), max_over_median=ts.max() / np.median(ts),
                            p90_over_median=np.percentile(ts, 90) / np.median(ts))])
    print(df.round(3).to_string(index=False))
    print('\nThe notebook cites chromatin_od_ranker_max_peaks50.ipynb for "1.09x-3.5x per-call '
          'spikes" on a single copyMakeBorder call. This is an independent re-measurement of '
          'that same spike behaviour on the same machine.')
    save(df, 'pad_microbenchmark')
    return df



# =====================================================================================
# CHECK 11 -- ground truth re-derived from databases/MIDOG++.json, and CSV internal arithmetic
# =====================================================================================
def check_ground_truth(T, P):
    banner('CHECK 11 -- n_gt_mitotic re-derived from the annotation DB; CSV internal arithmetic')
    sys.path.insert(0, REPO)
    from midog_utils import dataset as ds
    _, anns = ds.load_annotations(os.path.join(REPO, 'databases/MIDOG++.json'))
    rows = []
    for _, r in T.iterrows():
        g = anns[anns.file_name == r.file_name]
        n_mit_all = int((g.category_id == 1).sum())
        seed_is_mitotic = bool((g[(g.ann_id == r.seed_ann_id)].category_id == 1).all()
                               and (g.ann_id == r.seed_ann_id).any())
        want = n_mit_all - (1 if seed_is_mitotic else 0)
        got = P[(P.file_name == r.file_name) & (P.seed_index == r.seed_index)].n_gt_mitotic
        rows.append(dict(file_name=r.file_name, seed_index=int(r.seed_index),
                         seed_ann_id=int(r.seed_ann_id),
                         n_mitotic_in_json=n_mit_all, seed_is_mitotic=seed_is_mitotic,
                         expected_n_gt=want, observed_n_gt=int(got.iloc[0]),
                         matches=bool((got == want).all())))
    gdf = pd.DataFrame(rows)
    print(f'n_gt_mitotic == (JSON mitotic count - 1 for the excluded seed): '
          f'{int(gdf.matches.sum())}/{len(gdf)} draws')
    print(f'seed annotation is a mitotic figure on: {int(gdf.seed_is_mitotic.sum())}/{len(gdf)} draws')
    # internal arithmetic of the precision CSV
    e1 = np.abs(P.precision_at_budget - P.tp_at_budget / P.budget_delivered).max()
    e2 = np.abs(P.recall_at_budget - P.tp_at_budget / P.n_gt_mitotic).max()
    print(f'precision_at_budget == tp/budget_delivered on all {len(P)} rows, max err {e1:.2e}')
    print(f'recall_at_budget    == tp/n_gt_mitotic     on all {len(P)} rows, max err {e2:.2e}')
    print(f'tp_at_budget monotone non-decreasing in K within every (roi, seed, arm): '
          f'{bool(P.sort_values("budget").groupby(["file_name","seed_index","arm"]).tp_at_budget.apply(lambda s: bool((s.diff().dropna() >= 0).all())).all())}')
    print(f'tp_at_budget <= min(K, n_gt_mitotic) on all rows: '
          f'{bool((P.tp_at_budget <= np.minimum(P.budget, P.n_gt_mitotic)).all())}')
    print('\nRecall triad (Step 4.3): n_detections and precision_at_budget travel with '
          'recall_at_budget in the CSV; coverage_frac is computed by evaluate_arms and DROPPED by '
          "PRECISION_LONG's column selection. Pool is 87-100 candidates, not a tiling: "
          f'recall@30 spans {P[P.budget==30].recall_at_budget.min():.3f}-'
          f'{P[P.budget==30].recall_at_budget.max():.3f}, nowhere near saturation, so the '
          'saturated-pool failure mode does not apply here.')
    save(gdf, 'ground_truth')
    return gdf


# =====================================================================================
# CHECK 12 -- domain-level clustering (Step 4.2): 7 domains x 2 ROIs
# =====================================================================================
def check_domains(P):
    banner('CHECK 12 -- domain-level clustering (G=7, 2 ROIs per tumour type)')
    rows = []
    tabs = []
    for metric in ['precision_at_budget', 'recall_at_budget']:
        for k in BUDGETS:
            sub = P[P.budget == k].pivot_table(index=['file_name', 'tumor_type', 'seed_index'],
                                               columns='arm', values=metric)
            d = (sub['chromatin_od'] - sub['tm_score']).reset_index()
            d.columns = list(d.columns[:-1]) + ['delta']
            per_dom = d.groupby('tumor_type').delta.mean()
            x = per_dom.to_numpy()
            p_sf, G = exact_signflip(x)
            tcrit = stats.t.ppf(0.975, len(x) - 1)
            se = x.std(ddof=1) / np.sqrt(len(x))
            rows.append(dict(metric=metric, budget=k, unit='domain', G=len(x), mean=x.mean(),
                             ci_lo=x.mean() - tcrit * se, ci_hi=x.mean() + tcrit * se,
                             exact_signflip_p=p_sf, signflip_floor=2 / 2 ** len(x),
                             n_positive=int((x > 0).sum()),
                             majority_positive=bool((x > 0).sum() > len(x) / 2),
                             worst_domain=per_dom.idxmin(), worst_domain_delta=x.min()))
            if metric == 'precision_at_budget':
                tabs.append(per_dom.rename(f'precision_delta_K{k}'))
    dd = pd.DataFrame(rows)
    print(dd.to_string(index=False))
    dom_tab = pd.concat(tabs, axis=1).reset_index()
    print('\nper-domain mean precision delta over 2 ROIs x 5 seeds = 10 draws each:')
    print(dom_tab.to_string(index=False))
    print('\nAt G=7 the exact two-sided sign-flip floor is 2/2^7 = 0.0156; a p AT that floor '
          'means "as strong as 7 domains can resolve", not a measured tail.')
    save(dd, 'domain_level')
    save(dom_tab, 'domain_delta_table')
    return dd, dom_tab


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    gate_execution()
    gate_provenance()
    T, P, PS, S, strata = gate_composition()
    claim1(T, P)
    claim2(T)
    claim3(P, PS)
    claim4(P, S)
    claim5(T)
    check_ground_truth(T, P)
    check_domains(P)
    check_effective_units(T, P)
    check_figure(P, PS)
    check_config()
    claim6(T, P)
    check_pad_bench()
    banner(f'audit script finished in {time.time() - t0:.0f}s')


if __name__ == '__main__':
    main()
