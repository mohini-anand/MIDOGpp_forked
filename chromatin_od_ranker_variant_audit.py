"""Independent audit of `threshold_maxpeaks_ablation/chromatin_od_ranker_variant.ipynb`.

Run with the anaconda interpreter from the repo root:

    /Users/mohinianand/anaconda3/bin/python3 chromatin_od_ranker_variant_audit.py

Regenerates every table in `Research Logs/2026-09-12-chromatin-od-ranker-variant-audit.md`
into `results/chromatin_od_ranker_variant_audit_*.csv`.

Nothing here reads the notebook's printed output. Tier A recomputes from the two CSVs the
notebook wrote, from `results/precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv`
(the already-committed unbounded-pool run at the same seeds and the same config),
`threshold_maxpeaks_ablation/max_peaks_100_*.csv`, `latency_profiling/
chromatin_od_latency_per_roi.csv`, `databases/MIDOG++.json` and the TIFF headers.
Tier B (section G) recomputes `tp_at_budget` for three ROIs straight from the pixels with
its own peak extraction, NMS, od51, ranking and greedy matcher -- `midog_utils.compare`,
`midog_utils.evaluate` and `midog_utils.nms` are read but never called.
"""

from __future__ import annotations

import itertools
import json
import math
import os
import sys
import time

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.abspath(__file__))
ABL = os.path.join(REPO, 'threshold_maxpeaks_ablation')
OUT = os.path.join(REPO, 'results')
STEM = 'chromatin_od_ranker_variant_audit'
BUDGETS = (10, 20, 30)
RUN_TIER_B = os.environ.get('AUDIT_SKIP_TIER_B', '0') != '1'


def out(name, df):
    path = os.path.join(OUT, f'{STEM}_{name}.csv')
    df.to_csv(path, index=False)
    print(f'  -> {os.path.relpath(path, REPO)}  ({len(df)} rows)')


def head(t):
    print('\n' + '=' * 78 + f'\n{t}\n' + '=' * 78)


# ---------------------------------------------------------------------------------
# A. Composition gate on the two artifacts the notebook wrote
# ---------------------------------------------------------------------------------
head('A. Composition gate')

PREC = pd.read_csv(os.path.join(ABL, 'chromatin_od_ranker_precision.csv'))
TIME = pd.read_csv(os.path.join(ABL, 'chromatin_od_ranker_timing.csv'))

comp = []
comp.append(('precision rows', len(PREC), 14 * 2 * 3))
comp.append(('precision unique ROIs', PREC['file_name'].nunique(), 14))
comp.append(('precision unique branches', PREC['branch'].nunique(), 2))
comp.append(('precision unique budgets', PREC['budget'].nunique(), 3))
comp.append(('precision duplicate (roi,branch,budget) keys',
             int(PREC.duplicated(['file_name', 'branch', 'budget']).sum()), 0))
comp.append(('precision fully duplicate rows', int(PREC.duplicated().sum()), 0))
comp.append(('precision NaN cells', int(PREC.isna().sum().sum()), 0))
comp.append(('timing rows', len(TIME), 14))
comp.append(('timing duplicate file_name', int(TIME.duplicated(['file_name']).sum()), 0))
comp.append(('tumour domains in precision', PREC['tumor_type'].nunique(), 7))
comp.append(('ROIs per domain (min)', int(PREC.groupby('tumor_type')['file_name'].nunique().min()), 2))
comp.append(('ROIs per domain (max)', int(PREC.groupby('tumor_type')['file_name'].nunique().max()), 2))
comp.append(('rows with budget_delivered != budget',
             int((PREC['budget_delivered'] != PREC['budget']).sum()), 0))
comp.append(('rows where tp > budget_delivered',
             int((PREC['tp_at_budget'] > PREC['budget_delivered']).sum()), 0))
comp.append(('timing n_peaks != 100 (cap must bind pre-NMS)',
             int((TIME['n_peaks'] != 100).sum()), 0))
comp.append(('timing n_detections outside 89-99',
             int((~TIME['n_detections'].between(89, 99)).sum()), 0))
comp.append(('od51 nan_rate != 0', int((TIME['chromatin_od_nan_rate'] != 0).sum()), 0))
comp.append(('od51 largest_tie_block != 1',
             int((TIME['chromatin_od_largest_tie_block'] != 1).sum()), 0))

# precision_at_budget must equal tp / budget_delivered, and recall tp / n_gt_mitotic
p_err = float((PREC['precision_at_budget'] - PREC['tp_at_budget'] / PREC['budget_delivered']).abs().max())
r_err = float((PREC['recall_at_budget'] - PREC['tp_at_budget'] / PREC['n_gt_mitotic']).abs().max())
comp.append(('max |precision - tp/delivered| < 1e-12', float(p_err < 1e-12), 1.0))
comp.append(('max |recall - tp/n_gt| < 1e-12', float(r_err < 1e-12), 1.0))
print(f'  (float residuals: precision {p_err:.2e}, recall {r_err:.2e})')

COMP = pd.DataFrame(comp, columns=['check', 'observed', 'expected'])
COMP['pass'] = COMP['observed'] == COMP['expected']
print(COMP.to_string(index=False))
out('composition', COMP)

# per-ROI n_gt_mitotic and the domain map, straight from the annotation database
anns = json.load(open(os.path.join(REPO, 'databases', 'MIDOG++.json')))
img_by_id = {im['id']: im for im in anns['images']}
name_by_id = {im['id']: im['file_name'] for im in anns['images']}
roi_files = sorted(f for f in os.listdir(os.path.join(REPO, 'images', 'extra_valid'))
                   if f.endswith('.tiff'))
id_by_name = {v: k for k, v in name_by_id.items()}
gt_rows = []
for a in anns['annotations']:
    fn = name_by_id.get(a['image_id'])
    if fn in roi_files:
        x0, y0, x1, y1 = a['bbox']
        gt_rows.append(dict(file_name=fn, ann_id=a['id'], category_id=a['category_id'],
                            cx=(x0 + x1) / 2.0, cy=(y0 + y1) / 2.0,
                            n_votes=len(a.get('labels', [])),
                            n_mitotic_votes=sum(1 for l in a.get('labels', []) if l == 1)))
GT = pd.DataFrame(gt_rows)
dom = {im['file_name']: im.get('tumortype', im.get('tumor_type'))
       for im in anns['images'] if im['file_name'] in roi_files}

# n_gt_mitotic in the CSV is *after* the seed annotation is removed
seed_ids = dict(zip(TIME['file_name'], TIME['seed_ann_id']))
ngt_chk = []
for fn in roi_files:
    g = GT[GT['file_name'] == fn]
    n_mit_all = int((g['category_id'] == 1).sum())
    n_mit_eval = int(((g['category_id'] == 1) & (g['ann_id'] != seed_ids[fn])).sum())
    csv_n = int(PREC[PREC['file_name'] == fn]['n_gt_mitotic'].iloc[0])
    ngt_chk.append(dict(file_name=fn, domain_json=dom[fn],
                        domain_csv=PREC[PREC['file_name'] == fn]['tumor_type'].iloc[0],
                        n_mitotic_json=n_mit_all, n_mitotic_minus_seed=n_mit_eval,
                        n_gt_mitotic_csv=csv_n, match=n_mit_eval == csv_n,
                        seed_is_unanimous_mitotic=bool(
                            ((g['ann_id'] == seed_ids[fn]) & (g['category_id'] == 1) &
                             (g['n_votes'] == g['n_mitotic_votes'])).any())))
NGT = pd.DataFrame(ngt_chk)
NGT['domain_match'] = NGT['domain_json'] == NGT['domain_csv']
print('\nn_gt_mitotic and domain, re-derived from databases/MIDOG++.json:')
print(NGT.to_string(index=False))
print(f"  all n_gt match: {NGT['match'].all()} | all domains match: {NGT['domain_match'].all()} | "
      f"all seeds unanimous mitotic: {NGT['seed_is_unanimous_mitotic'].all()}")
out('gt_composition', NGT)

# ---------------------------------------------------------------------------------
# B. Claim 2 -- pooled precision, win counts, and the ROI-level test the notebook omits
# ---------------------------------------------------------------------------------
head('B. Claim 2 -- pooled precision and win counts (Tier A)')

wide = PREC.pivot_table(index=['file_name', 'tumor_type', 'budget'], columns='branch',
                        values=['tp_at_budget', 'precision_at_budget', 'budget_delivered'])
wide.columns = [f'{a}_{b}' for a, b in wide.columns]
wide = wide.reset_index()
wide['d_precision'] = wide['precision_at_budget_variant'] - wide['precision_at_budget_baseline']
wide['d_tp'] = wide['tp_at_budget_variant'] - wide['tp_at_budget_baseline']

rows = []
for k in BUDGETS:
    s = wide[wide['budget'] == k]
    pb = s['tp_at_budget_baseline'].sum() / s['budget_delivered_baseline'].sum()
    pv = s['tp_at_budget_variant'].sum() / s['budget_delivered_variant'].sum()
    d = s['d_precision'].to_numpy()
    up, down, tie = int((d > 0).sum()), int((d < 0).sum()), int((d == 0).sum())
    g_eff = up + down
    # exact two-sided sign test over the non-tied ROIs
    k_up = up
    p_sign = 2.0 * sum(math.comb(g_eff, i) for i in range(min(k_up, g_eff - k_up) + 1)) / 2 ** g_eff
    p_sign = min(1.0, p_sign)
    # exact two-sided sign-flip (permutation) test on the delta magnitudes, all 14 ROIs
    dd = wide[wide['budget'] == k]['d_precision'].to_numpy()
    obs = dd.mean()
    n = len(dd)
    cnt = 0
    for signs in itertools.product([1, -1], repeat=n):
        if abs((dd * np.array(signs)).mean()) >= abs(obs) - 1e-12:
            cnt += 1
    p_flip = cnt / 2 ** n
    # mean per-ROI delta with a t interval on G-1 df
    sd = dd.std(ddof=1)
    from scipy import stats as sps
    tcrit = sps.t.ppf(0.975, n - 1)
    lo, hi = obs - tcrit * sd / np.sqrt(n), obs + tcrit * sd / np.sqrt(n)
    rows.append(dict(budget=k, pooled_baseline=round(pb, 4), pooled_variant=round(pv, 4),
                     pooled_delta=round(pv - pb, 4),
                     pooled_pct_change=round((pv - pb) / pb * 100, 1),
                     mean_per_roi_delta=round(obs, 4),
                     t_ci95_lo=round(lo, 4), t_ci95_hi=round(hi, 4),
                     up=up, down=down, tied=tie, G_effective=g_eff,
                     p_exact_sign=round(p_sign, 5), p_signflip_exact=round(p_flip, 5),
                     signflip_floor=round(2.0 / 2 ** n, 6)))
SUM = pd.DataFrame(rows)
print(SUM.to_string(index=False))
print('\ndomain-level (Step 4.2) -- 2 ROIs per domain, G=7:')
print('\npooled precision == unweighted mean of per-ROI precision (budget_delivered == K everywhere):')
for k in BUDGETS:
    s = wide[wide['budget'] == k]
    print(f'  K={k}: baseline pooled {s["tp_at_budget_baseline"].sum()/s["budget_delivered_baseline"].sum():.6f} '
          f'vs mean-of-ROIs {s["precision_at_budget_baseline"].mean():.6f} | '
          f'variant pooled {s["tp_at_budget_variant"].sum()/s["budget_delivered_variant"].sum():.6f} '
          f'vs mean-of-ROIs {s["precision_at_budget_variant"].mean():.6f}')
out('pooled_precision', SUM)
out('per_roi_delta', wide.sort_values(['budget', 'file_name']))

# Step 4.2 -- the same test clustered at the DOMAIN instead of the ROI. 2 ROIs per domain,
# so G falls from 14 to 7 and the two-sided sign-flip floor rises to 2/2**7 = 0.0156.
drows = []
for k in BUDGETS:
    s = wide[wide['budget'] == k]
    g = s.groupby('tumor_type')['d_precision'].mean()
    up, down = int((g > 0).sum()), int((g < 0).sum())
    n = len(g)
    obs = g.mean()
    cnt = sum(1 for sg in itertools.product([1, -1], repeat=n)
              if abs((g.to_numpy() * np.array(sg)).mean()) >= abs(obs) - 1e-12)
    for dom_name, val in g.items():
        drows.append(dict(budget=k, tumor_type=dom_name, domain_mean_delta=round(val, 4)))
    print(f'  K={k}: domain-level G={n}, {up} up / {down} down, mean delta {obs:+.4f}, '
          f'exact sign-flip p={cnt / 2 ** n:.4f} (floor {2 / 2 ** n:.4f}); '
          f'domains going down: {sorted(g[g < 0].index)}')
print('domain-level means (2 ROIs each):')
DOM = pd.DataFrame(drows)
print(DOM.pivot(index='tumor_type', columns='budget', values='domain_mean_delta').to_string())
out('domain_level_delta', DOM)

# ---------------------------------------------------------------------------------
# C. Claim 3 -- every (ROI, budget) pair where the variant lost
# ---------------------------------------------------------------------------------
head('C. Claim 3 -- decliners, enumerated over all 42 (ROI, budget) pairs')

DEC = wide[wide['d_precision'] < 0][['file_name', 'tumor_type', 'budget',
                                     'tp_at_budget_baseline', 'tp_at_budget_variant',
                                     'precision_at_budget_baseline',
                                     'precision_at_budget_variant', 'd_precision']]
DEC = DEC.sort_values(['budget', 'file_name'])
print(f'{len(wide)} (ROI, budget) pairs total; {len(DEC)} with variant < baseline:')
print(DEC.to_string(index=False))
print('\nunchanged pairs:')
print(wide[wide['d_precision'] == 0][['file_name', 'budget', 'tp_at_budget_baseline',
                                      'tp_at_budget_variant']].to_string(index=False))
lymph = sorted(PREC[PREC['tumor_type'] == 'canine lymphosarcoma']['file_name'].unique())
print(f"\ncanine lymphosarcoma ROIs in the set: {lymph}; "
      f"present among decliners: {sorted(set(DEC['file_name']) & set(lymph))}")
out('decliners', DEC)

# ---------------------------------------------------------------------------------
# D. Claim 1 -- timing decomposition and the outlier attribution
# ---------------------------------------------------------------------------------
head('D. Claim 1 -- timing (Tier A on the notebook timing CSV)')

T = TIME.set_index('file_name')
crit = T['t6b_od51_loop_s'] + T['t6c_variant_rank_s']
stage6 = T['t_chromatin_od_overhead_s']
trows = [
    dict(quantity='t6a_od_pad_ms (pad only)', mean=T['t6a_od_pad_ms'].mean(),
         median=T['t6a_od_pad_ms'].median(), sd=T['t6a_od_pad_ms'].std(),
         min=T['t6a_od_pad_ms'].min(), max=T['t6a_od_pad_ms'].max()),
    dict(quantity='t6b_od51_loop_ms (loop only)', mean=T['t6b_od51_loop_ms'].mean(),
         median=T['t6b_od51_loop_ms'].median(), sd=T['t6b_od51_loop_ms'].std(),
         min=T['t6b_od51_loop_ms'].min(), max=T['t6b_od51_loop_ms'].max()),
    dict(quantity='t6c_variant_rank_ms (resort only)', mean=T['t6c_variant_rank_ms'].mean(),
         median=T['t6c_variant_rank_ms'].median(), sd=T['t6c_variant_rank_ms'].std(),
         min=T['t6c_variant_rank_ms'].min(), max=T['t6c_variant_rank_ms'].max()),
    dict(quantity='t6b+t6c_ms (THE CRITERION)', mean=crit.mean() * 1e3,
         median=crit.median() * 1e3, sd=crit.std() * 1e3,
         min=crit.min() * 1e3, max=crit.max() * 1e3),
    dict(quantity='t6a+t6b+t6c_ms (stage-6 overhead)', mean=stage6.mean() * 1e3,
         median=stage6.median() * 1e3, sd=stage6.std() * 1e3,
         min=stage6.min() * 1e3, max=stage6.max() * 1e3),
    dict(quantity='t6_baseline_ms (sort by tm_score)', mean=T['t6_baseline_ms'].mean(),
         median=T['t6_baseline_ms'].median(), sd=T['t6_baseline_ms'].std(),
         min=T['t6_baseline_ms'].min(), max=T['t6_baseline_ms'].max()),
]
TD = pd.DataFrame(trows).round(3)
print(TD.to_string(index=False))
d_pipe = T['t_variant_pipeline_ms'].mean() - T['t_baseline_pipeline_ms'].mean()
print(f"\nfull pipeline mean: baseline {T['t_baseline_pipeline_ms'].mean():.2f}ms -> "
      f"variant {T['t_variant_pipeline_ms'].mean():.2f}ms  "
      f"delta {d_pipe:+.2f}ms = {d_pipe / T['t_baseline_pipeline_ms'].mean() * 100:+.3f}%")
print(f"stage-6 overhead / baseline pipeline = "
      f"{stage6.mean() / T['t_baseline_pipeline_s'].mean() * 100:.3f}%")
print(f"share of stage-6 overhead that is the pad: {T['t6a_od_pad_s'].mean() / stage6.mean() * 100:.1f}%; "
      f"criterion: {crit.mean() / stage6.mean() * 100:.1f}%")
out('timing_decomposition', TD)

# ROI pixel dimensions, from the TIFF headers
import tifffile
dims = []
for fn in roi_files:
    with tifffile.TiffFile(os.path.join(REPO, 'images', 'extra_valid', fn)) as tf:
        s = tf.series[0].shape
    dims.append(dict(file_name=fn, height=int(s[0]), width=int(s[1]),
                     megapixels=round(s[0] * s[1] / 1e6, 2)))
DIM = pd.DataFrame(dims).set_index('file_name')

LAT = pd.read_csv(os.path.join(REPO, 'latency_profiling',
                               'chromatin_od_latency_per_roi.csv')).set_index('file_name')
PAD = DIM.join(T[['t6a_od_pad_ms', 't6b_od51_loop_ms', 't_setup_ms', 't3_template_matching_ms',
                  'n_detections']])
PAD['pad_ms_latency_run'] = LAT['t7a_od_pad_s'] * 1e3
PAD['od51_loop_ms_latency_run'] = LAT['t7b_od51_loop_s'] * 1e3
PAD['n_det_latency_run'] = LAT['n_detections']
PAD['shape'] = PAD['height'].astype(str) + 'x' + PAD['width'].astype(str)
PAD['rank_pad_this_run'] = PAD['t6a_od_pad_ms'].rank(ascending=False).astype(int)
PAD['rank_pad_latency_run'] = PAD['pad_ms_latency_run'].rank(ascending=False).astype(int)
PAD = PAD.sort_values('t6a_od_pad_ms', ascending=False)
print('\nROI pixel dimensions against pad cost, in both independent runs:')
print(PAD[['shape', 'megapixels', 't6a_od_pad_ms', 'rank_pad_this_run',
           'pad_ms_latency_run', 'rank_pad_latency_run', 'n_detections',
           'n_det_latency_run']].to_string())

big = PAD[PAD['shape'] == '5412x7215']
small = PAD[PAD['shape'] != '5412x7215']
print(f"\n  5412x7215 ROIs (n={len(big)}): {sorted(big.index)}")
print(f"  this run   -- large pad ms min {big['t6a_od_pad_ms'].min():.2f} vs "
      f"small pad ms max {small['t6a_od_pad_ms'].max():.2f}  -> "
      f"{'PERFECT separation' if big['t6a_od_pad_ms'].min() > small['t6a_od_pad_ms'].max() else 'overlap'}")
print(f"  latency run -- large pad ms min {big['pad_ms_latency_run'].min():.2f} vs "
      f"small pad ms max {small['pad_ms_latency_run'].max():.2f}  -> "
      f"{'PERFECT separation' if big['pad_ms_latency_run'].min() > small['pad_ms_latency_run'].max() else 'overlap'}")
n_big, n_tot = len(big), len(PAD)
p_perm = 1.0 / math.comb(n_tot, n_big)
print(f"  exact permutation p for a perfect {n_big}-vs-{n_tot - n_big} split by dimension: "
      f"{p_perm:.5f} (one-sided), in each run independently")
print(f"  spike ratio within the 5412x7215 group, this run: "
      f"{big['t6a_od_pad_ms'].max() / big['t6a_od_pad_ms'].median():.2f}x; "
      f"latency run: {big['pad_ms_latency_run'].max() / big['pad_ms_latency_run'].median():.2f}x")
out('pad_vs_dimensions', PAD.reset_index())

# pool-size independence of the pad, and pool-size dependence of the loop
print(f"\npad, this run (pool {T['n_detections'].mean():.0f}) mean {T['t6a_od_pad_ms'].mean():.1f}ms "
      f"vs latency run (pool {LAT['n_detections'].mean():.0f}) mean {PAD['pad_ms_latency_run'].mean():.1f}ms")
print(f"od51 loop, this run mean {T['t6b_od51_loop_ms'].mean():.2f}ms "
      f"vs latency run mean {PAD['od51_loop_ms_latency_run'].mean():.1f}ms  "
      f"(pool ratio {LAT['n_detections'].mean() / T['n_detections'].mean():.0f}x, "
      f"time ratio {PAD['od51_loop_ms_latency_run'].mean() / T['t6b_od51_loop_ms'].mean():.0f}x)")
print(f"per-candidate od51 cost: this run "
      f"{T['t6b_od51_loop_ms'].sum() / T['n_detections'].sum() * 1000:.1f}us, latency run "
      f"{PAD['od51_loop_ms_latency_run'].sum() / LAT['n_detections'].sum() * 1000:.1f}us")

# ---------------------------------------------------------------------------------
# E. Claim 4 -- the cascade framing, tested against the unbounded-pool run
# ---------------------------------------------------------------------------------
head('E. Claim 4 -- capped cascade vs the already-committed unbounded pool')

UNB = pd.read_csv(os.path.join(OUT, 'precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv'))
UNB = UNB[UNB['arm'].isin(['tm_score', 'chromatin_od']) & UNB['budget'].isin(BUDGETS)]

# provenance: same seeds, same response map, only the cap differs
seed_chk = (UNB[UNB['arm'] == 'tm_score'][['file_name', 'seed_ann_id', 'base_size',
                                           'map_median', 'mad_scale', 'n_gt_mitotic']]
            .drop_duplicates().set_index('file_name'))
this_chk = T[['seed_ann_id', 'base_size', 'map_median', 'mad_scale']].copy()
this_chk['n_gt_mitotic'] = PREC.groupby('file_name')['n_gt_mitotic'].first()
J = this_chk.join(seed_chk, lsuffix='_capped', rsuffix='_unbounded')
same = pd.DataFrame({c: J[f'{c}_capped'] == J[f'{c}_unbounded']
                     for c in ['seed_ann_id', 'base_size', 'map_median', 'mad_scale', 'n_gt_mitotic']})
print('identical seed / template / response-map statistics between the two runs:')
print(same.all().to_string())
print(f'  -> every column identical on all {len(J)} ROIs: {bool(same.values.all())}')

casc = []
for k in BUDGETS:
    row = dict(budget=k)
    for arm, lab in [('tm_score', 'tm'), ('chromatin_od', 'od51')]:
        u = UNB[(UNB['arm'] == arm) & (UNB['budget'] == k)]
        row[f'unbounded_{lab}'] = round(u['tp_at_budget'].sum() / u['budget_delivered'].sum(), 4)
        b = 'baseline' if arm == 'tm_score' else 'variant'
        c = PREC[(PREC['branch'] == b) & (PREC['budget'] == k)]
        row[f'capped_{lab}'] = round(c['tp_at_budget'].sum() / c['budget_delivered'].sum(), 4)
    row['od51_gain_unbounded'] = round(row['unbounded_od51'] - row['unbounded_tm'], 4)
    row['od51_gain_capped'] = round(row['capped_od51'] - row['capped_tm'], 4)
    row['gain_attributable_to_cap'] = round(row['od51_gain_capped'] - row['od51_gain_unbounded'], 4)
    row['share_of_capped_gain_present_uncapped_pct'] = round(
        row['od51_gain_unbounded'] / row['od51_gain_capped'] * 100, 1)
    casc.append(row)
CASC = pd.DataFrame(casc)
print('\npooled precision@K, same 14 ROIs, same seeds, same config -- only max_peaks differs:')
print(CASC.to_string(index=False))
out('cascade_vs_unbounded', CASC)

# per-ROI: does the capped pool's od51 list differ from the unbounded pool's?
perroi = []
for k in BUDGETS:
    u = UNB[(UNB['arm'] == 'chromatin_od') & (UNB['budget'] == k)].set_index('file_name')
    c = PREC[(PREC['branch'] == 'variant') & (PREC['budget'] == k)].set_index('file_name')
    for fn in roi_files:
        perroi.append(dict(file_name=fn, budget=k, tp_unbounded_od51=int(u.loc[fn, 'tp_at_budget']),
                           tp_capped_od51=int(c.loc[fn, 'tp_at_budget']),
                           delta=int(c.loc[fn, 'tp_at_budget']) - int(u.loc[fn, 'tp_at_budget'])))
PR = pd.DataFrame(perroi)
print('\nper-ROI tp@K, unbounded od51 vs capped od51:')
for k in BUDGETS:
    s = PR[PR['budget'] == k]
    print(f"  K={k}: cap helps {int((s['delta'] > 0).sum())}, hurts {int((s['delta'] < 0).sum())}, "
          f"ties {int((s['delta'] == 0).sum())} of 14   (sum delta {int(s['delta'].sum())})")
out('cascade_per_roi', PR)

# largest_tie_block, capped vs unbounded -- claim 5
TIE = pd.DataFrame({
    'capped_largest_tie_block': T['chromatin_od_largest_tie_block'],
    'capped_n_detections': T['n_detections'],
    'unbounded_largest_tie_block': UNB[(UNB['arm'] == 'chromatin_od') & (UNB['budget'] == 10)]
        .set_index('file_name')['largest_tie_block'],
    'unbounded_n_detections': UNB[(UNB['arm'] == 'chromatin_od') & (UNB['budget'] == 10)]
        .set_index('file_name')['n_detections'],
    'unbounded_coverage_frac': UNB[(UNB['arm'] == 'chromatin_od') & (UNB['budget'] == 10)]
        .set_index('file_name')['coverage_frac'],
})
print('\nClaim 5 -- od51 largest_tie_block, capped pool vs the unbounded oracle:')
print(TIE.to_string())
print(f"  capped: {TIE['capped_largest_tie_block'].min()}-{TIE['capped_largest_tie_block'].max()}; "
      f"unbounded: {TIE['unbounded_largest_tie_block'].min()}-{TIE['unbounded_largest_tie_block'].max()}")
out('tie_blocks', TIE.reset_index())

# ---------------------------------------------------------------------------------
# F. Verification 1, re-derived without the notebook
# ---------------------------------------------------------------------------------
head('F. Verification 1, re-derived independently')

OT = pd.read_csv(os.path.join(ABL, 'max_peaks_100_timing.csv')).set_index('file_name')
OP = pd.read_csv(os.path.join(ABL, 'max_peaks_100_precision.csv'))
OPV = OP[OP['branch'] == 'variant'].set_index(['file_name', 'budget'])

vrows, ncmp, ndiv = [], 0, 0
for fn in roi_files:
    for col, ocol in [('seed_ann_id', 'seed_ann_id'), ('base_size', 'base_size'),
                      ('n_detections', 'n_detections_variant'), ('map_median', 'map_median'),
                      ('mad_scale', 'mad_scale')]:
        a, b = T.loc[fn, col], OT.loc[fn, ocol]
        ncmp += 1
        ok = np.isclose(a, b, rtol=0, atol=1e-9)
        ndiv += (not ok)
        vrows.append(dict(file_name=fn, column=col, this=a, oracle=b, equal=bool(ok)))
    for k in BUDGETS:
        a = int(PREC[(PREC['file_name'] == fn) & (PREC['branch'] == 'baseline') &
                     (PREC['budget'] == k)]['tp_at_budget'].iloc[0])
        b = int(OPV.loc[(fn, k), 'tp_at_budget'])
        ncmp += 1
        ndiv += (a != b)
        vrows.append(dict(file_name=fn, column=f'tp_at_{k}', this=a, oracle=b, equal=a == b))
VER = pd.DataFrame(vrows)
print(f'{ncmp} values compared against max_peaks_100_*.csv, {ndiv} divergences')
out('verification1', VER)

# the tm_score baseline also has a second, unbounded oracle
ncmp2 = ndiv2 = 0
v2 = []
for fn in roi_files:
    for k in BUDGETS:
        a = int(PREC[(PREC['file_name'] == fn) & (PREC['branch'] == 'baseline') &
                     (PREC['budget'] == k)]['tp_at_budget'].iloc[0])
        b = int(UNB[(UNB['arm'] == 'tm_score') & (UNB['budget'] == k) &
                    (UNB['file_name'] == fn)]['tp_at_budget'].iloc[0])
        ncmp2 += 1
        ndiv2 += (a != b)
        v2.append(dict(file_name=fn, budget=k, capped_tm_tp=a, unbounded_tm_tp=b, equal=a == b))
V2 = pd.DataFrame(v2)
print(f'{ncmp2} tm_score tp@K values compared against the UNBOUNDED run, {ndiv2} divergences '
      f'(D9 claims the cap is free through K=30)')
out('verification2_unbounded_tm', V2)

# ---------------------------------------------------------------------------------
# G. Tier B -- tp_at_budget straight from the pixels, and a re-timed pad
# ---------------------------------------------------------------------------------
if RUN_TIER_B:
    head('G. Tier B -- recomputed from pixels with an independent matcher')
    import cv2
    sys.path.insert(0, REPO)
    from midog_utils import channels as ch
    from midog_utils import dataset as ds
    from midog_utils import seed_selection as ss
    from midog_utils import template_match as tm

    MPP = {}
    IMAGES = os.path.join(REPO, 'images', 'extra_valid')

    def my_od51(hem_pad, xs, ys, window=51, frac=0.10):
        half = window // 2
        k = max(1, int(frac * window * window))
        vals = []
        for x, y in zip(xs, ys):
            ix, iy = int(round(x)), int(round(y))
            p = hem_pad[iy - half:iy + half + 1, ix - half:ix + half + 1]
            if p.shape != (window, window):
                vals.append(float('nan'))
                continue
            f = p.ravel()
            vals.append(float(np.partition(f, -k)[-k:].mean()))
        return np.array(vals)

    def my_nms(cx, cy, sc, radius):
        order = np.argsort(-sc, kind='stable')
        keep, alive = [], np.ones(len(sc), dtype=bool)
        for i in order:
            if not alive[i]:
                continue
            keep.append(i)
            alive &= np.hypot(cx - cx[i], cy - cy[i]) > radius
            alive[i] = False
        return np.array(keep, dtype=int)

    def my_greedy_tp(order_idx, cx, cy, gx, gy, gcat, radius):
        """Greedy one-to-one match in rank order; returns cumulative mitotic TP."""
        taken = np.zeros(len(gx), dtype=bool)
        tp = []
        run = 0
        for i in order_idx:
            d = np.hypot(gx - cx[i], gy - cy[i])
            cand = np.where((d <= radius) & (~taken))[0]
            if len(cand):
                g = cand[int(np.argmin(d[cand]))]
                taken[g] = True
                if gcat[g] == 1:
                    run += 1
            tp.append(run)
        return np.array(tp)

    tierb, padrows = [], []
    TIER_B_ROIS = ['245.tiff', '246.tiff', '301.tiff']
    PAD_ROIS = ['013.tiff', '403.tiff', '548.tiff', '301.tiff', '094.tiff']
    t_tierb0 = time.time()
    for fn in sorted(set(TIER_B_ROIS) | set(PAD_ROIS)):
        path = os.path.join(IMAGES, fn)
        rgb = ds.load_roi(path)
        mpp = ds.roi_mpp(path)
        hem = ch.to_channel(rgb, 'hematoxylin_od')
        H, W = hem.shape[:2]
        # --- re-timed pad, 5 repeats
        reps = []
        for _ in range(5):
            t0 = time.perf_counter()
            hp = cv2.copyMakeBorder(hem, 25, 25, 25, 25, cv2.BORDER_REPLICATE)
            reps.append((time.perf_counter() - t0) * 1e3)
            del hp
        padrows.append(dict(file_name=fn, shape=f'{H}x{W}', megapixels=round(H * W / 1e6, 2),
                            pad_ms_min=round(min(reps), 2), pad_ms_median=round(float(np.median(reps)), 2),
                            pad_ms_max=round(max(reps), 2),
                            notebook_t6a_ms=float(T.loc[fn, 't6a_od_pad_ms']),
                            latency_run_pad_ms=round(float(LAT.loc[fn, 't7a_od_pad_s'] * 1e3), 2)))
        if fn in TIER_B_ROIS:
            gray_inv = ch.to_gray_inverted(rgb)
            g = GT[GT['file_name'] == fn].reset_index(drop=True)
            sid = int(T.loc[fn, 'seed_ann_id'])
            srow = g[g['ann_id'] == sid].iloc[0]
            r = ss.tightened_template_box(gray_inv, float(srow['cx']), float(srow['cy']),
                                          otsu_window=51)
            base_size, tcx, tcy = r
            assert base_size == int(T.loc[fn, 'base_size']), fn
            patch = tm.read_padded_patch(hem, float(tcx), float(tcy), 73)
            c0 = patch.shape[0] // 2
            hb = base_size // 2
            tmpl = np.ascontiguousarray(patch[c0 - hb:c0 + hb + 1, c0 - hb:c0 + hb + 1])
            P = (base_size - 1) // 2
            hem_p = cv2.copyMakeBorder(hem, P, P, P, P, cv2.BORDER_REPLICATE)
            fused = np.asarray(cv2.matchTemplate(hem_p, tmpl, cv2.TM_CCOEFF), dtype=np.float32)
            assert fused.shape == (H, W), (fn, fused.shape, (H, W))
            samp = fused[::8, ::8].ravel()
            samp = samp[np.isfinite(samp)]
            med = float(np.median(samp))
            mad = float(1.4826 * np.median(np.abs(samp - med)))
            cut = med - 1.5 * mad
            kk = 2 * 7 + 1
            dil = cv2.dilate(fused, np.ones((kk, kk), np.uint8))
            mask = (fused >= dil) & (fused >= cut)
            ys, xs = np.nonzero(mask)
            sc = fused[ys, xs]
            order = np.lexsort((ys, xs, -sc))[:100]
            pcx, pcy, psc = xs[order].astype(float), ys[order].astype(float), sc[order]
            nms_r = 7.5 / mpp
            keep = my_nms(pcx, pcy, psc, nms_r)
            kx, ky, ks = pcx[keep], pcy[keep], psc[keep]
            ok = np.hypot(kx - tcx, ky - tcy) > 5.0
            kx, ky, ks = kx[ok], ky[ok], ks[ok]
            hem_pad = cv2.copyMakeBorder(hem, 25, 25, 25, 25, cv2.BORDER_REPLICATE)
            od = my_od51(hem_pad, kx + 25, ky + 25)
            ge = g[g['ann_id'] != sid].reset_index(drop=True)
            gx, gy, gc = ge['cx'].to_numpy(), ge['cy'].to_numpy(), ge['category_id'].to_numpy()
            radius = 7.5 / mpp
            res = dict(file_name=fn, n_detections=len(kx),
                       n_detections_csv=int(T.loc[fn, 'n_detections']),
                       med_this=round(med, 5), med_csv=float(T.loc[fn, 'map_median']),
                       mad_this=round(mad, 5), mad_csv=float(T.loc[fn, 'mad_scale']),
                       od51_nan=int(np.isnan(od).sum()),
                       od51_largest_tie=int(pd.Series(od).value_counts().iloc[0]))
            for key, vals, branch in [('tm', ks, 'baseline'), ('od51', od, 'variant')]:
                idx = np.argsort(-vals, kind='mergesort')
                tpc = my_greedy_tp(idx, kx, ky, gx, gy, gc, radius)
                for k in BUDGETS:
                    dlv = min(k, len(tpc))
                    mine = int(tpc[dlv - 1])
                    theirs = int(PREC[(PREC['file_name'] == fn) & (PREC['branch'] == branch) &
                                      (PREC['budget'] == k)]['tp_at_budget'].iloc[0])
                    res[f'tp{k}_{key}_audit'] = mine
                    res[f'tp{k}_{key}_csv'] = theirs
            tierb.append(res)
        del rgb, hem
    TB = pd.DataFrame(tierb)
    PADB = pd.DataFrame(padrows)
    print(f'Tier B wall clock: {time.time() - t_tierb0:.0f}s')
    print('\nre-timed cv2.copyMakeBorder(pad=25) on hematoxylin_od, 5 repeats each:')
    print(PADB.sort_values('megapixels', ascending=False).to_string(index=False))
    print(f"  repeat-stabilised pad over these {len(PADB)} ROIs: "
          f"mean-of-medians {PADB['pad_ms_median'].mean():.2f}ms, mean-of-mins "
          f"{PADB['pad_ms_min'].mean():.2f}ms, against the notebook's single-shot mean "
          f"{T['t6a_od_pad_ms'].mean():.2f}ms / median {T['t6a_od_pad_ms'].median():.2f}ms")
    print(f"  within-ROI max/min ratio across 5 repeats: "
          f"{(PADB['pad_ms_max'] / PADB['pad_ms_min']).min():.2f}x-"
          f"{(PADB['pad_ms_max'] / PADB['pad_ms_min']).max():.2f}x; ROIs that hit >90ms on at "
          f"least one repeat: {sorted(PADB[PADB['pad_ms_max'] > 90]['file_name'])}")
    print(f"  ms per megapixel (median repeat): "
          f"{(PADB['pad_ms_median'] / PADB['megapixels']).min():.3f}-"
          f"{(PADB['pad_ms_median'] / PADB['megapixels']).max():.3f}")
    # what the notebook's headline becomes on repeat-stabilised pad timings
    pad_rob = PADB['pad_ms_median'].mean()
    stage6_rob = pad_rob + crit.mean() * 1e3
    print(f"  -> stage-6 overhead on repeat-stabilised pad: {stage6_rob:.1f}ms = "
          f"{stage6_rob / T['t_baseline_pipeline_ms'].mean() * 100:.2f}% of the baseline pipeline "
          f"(notebook reports {stage6.mean()*1e3:.1f}ms = "
          f"{stage6.mean() / T['t_baseline_pipeline_s'].mean() * 100:.2f}%)")

    # Step 4 mode 10 -- selection breadth: chromatin_od's rank among the six axes already
    # scored on these same 14 ROIs at these same budgets, on the unbounded pool
    six = pd.read_csv(os.path.join(OUT, 'precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv'))
    six = six[six['budget'].isin(BUDGETS)]
    pool6 = (six.groupby(['budget', 'arm'])
             .apply(lambda g: g['tp_at_budget'].sum() / g['budget_delivered'].sum(),
                    include_groups=False).unstack().round(4))
    print('\nselection breadth -- pooled precision@K of all six axes already scored on these '
          '14 ROIs (unbounded pool):')
    print(pool6.T.to_string())
    print('  rank of chromatin_od (1 = best of six): '
          + ', '.join(f'K={k}: {int(pool6.T.rank(ascending=False).loc["chromatin_od", k])}'
                      for k in BUDGETS))
    out('selection_breadth_six_axes', pool6.T.reset_index())
    print('\ntp_at_budget recomputed from pixels with an independent matcher:')
    cols = ['file_name', 'n_detections', 'n_detections_csv', 'od51_nan', 'od51_largest_tie'] + \
           [f'tp{k}_{a}_{s}' for k in BUDGETS for a in ('tm', 'od51') for s in ('audit', 'csv')]
    print(TB[cols].to_string(index=False))
    ndiff = sum(int(TB[f'tp{k}_{a}_audit'].ne(TB[f'tp{k}_{a}_csv']).sum())
                for k in BUDGETS for a in ('tm', 'od51'))
    nvals = len(TB) * len(BUDGETS) * 2
    print(f'  {nvals} tp_at_budget values recomputed from pixels, {ndiff} divergences')
    print(f'  n_detections: {int(TB["n_detections"].ne(TB["n_detections_csv"]).sum())} divergences')
    out('tier_b_from_pixels', TB)
    out('tier_b_pad_retimed', PADB)
else:
    print('\n[Tier B skipped -- AUDIT_SKIP_TIER_B=1]')

head('done')
