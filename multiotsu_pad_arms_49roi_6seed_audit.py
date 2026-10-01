"""
    Independent audit of production_hematoxylin_only/multiotsu_pad_arms/multiotsu_pad_arms_49roi_6seed.ipynb.

    Tier A re-derives every persisted number from the notebook's CSVs, the window sweep's reference CSVs, the
    fill-geometry CSVs and databases/MIDOG++.json, with matching, statistics and flags written here -- no
    midog_utils helper is imported anywhere in this file.
    Tier B goes back to the pixels: it re-derives the multi-Otsu box on all 294 clicks, rebuilds every template
    (4 arms + 3 references) and checks its SHA-1, re-runs a re-implemented search on the 17 stress clicks, runs the
    production (D11) blank on a fixed 7-ROI subset, and re-times one search.

    Writes results/multiotsu_pad_arms_49roi_6seed_audit_*.csv.
    Run: /Users/mohinianand/anaconda3/bin/python3 multiotsu_pad_arms_49roi_6seed_audit.py
"""

import hashlib
import json
import os
import subprocess
import time
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from scipy import stats as sps
from scipy.spatial import cKDTree
from skimage.color import rgb2hed
from skimage.filters import threshold_multiotsu

T0 = time.time()
REPO = Path(__file__).resolve().parent
os.chdir(REPO)
NB = 'production_hematoxylin_only/multiotsu_pad_arms/multiotsu_pad_arms_49roi_6seed.ipynb'
MO = 'results/precision_at_k_49roi_6seed_mopad_postD11'
SW = 'results/precision_at_k_49roi_6seed_windowsweep_postD11'
FG = 'results/precision_at_k_49roi_6seed_fillgeom_postD11'
OUT = 'results/multiotsu_pad_arms_49roi_6seed_audit'
SUBSET_DIRS = {'original_14': 'images/extra_valid', 'testing_35': 'images/extra_valid/testing_set'}

ARMS = ('mo_box', 'mo_pad51', 'mo_pad61', 'mo_ring5')
REFS = ('default_51', 'fill_51', 'fill_61')
ALL_ARMS = ARMS + REFS
AXES = ('chromatin_od', 'tm_score')
KEY = ['file_name', 'seed_index', 'condition', 'rank_key']
TESTED = {'mo_pad51 - mo_box': 'pad', 'mo_pad61 - mo_box': 'pad', 'mo_ring5 - mo_box': 'pad', 'mo_pad51 - fill_51': 'vs_incumbent', 'mo_pad61 - fill_61': 'vs_incumbent'}
DESCRIPTIVE = ['mo_box - default_51', 'mo_pad51 - default_51', 'mo_pad61 - default_51', 'mo_ring5 - default_51', 'fill_51 - default_51', 'fill_61 - default_51', 'mo_pad61 - mo_pad51', 'mo_ring5 - mo_pad51', 'mo_ring5 - mo_pad61', 'fill_61 - fill_51']
METRIC_K = {'tp_at_10': 10, 'tp_at_10_excl_border': 10, 'tp_at_20': 20, 'tp_at_30': 30}
N_BOOT, BOOT_SEED = 10_000, 20260917
HALF_BLANK, EDGE_BAND, SELF_PX, NEAR_PX = 25, 35, 2, 60
LEDGER = []


def ledger(section, claim, n_compared, n_diverged, tier, note=''):
    """
        Record one reproduction check.

        section (str): where in the audit the check lives.
        claim (str): what was checked.
        n_compared (int): values compared.
        n_diverged (int): values that disagree.
        tier (str): A, B or consistency.
        note (str): anything a reader needs.

        Returns None.
    """
    LEDGER.append(dict(section=section, claim=claim, n_compared=int(n_compared), n_diverged=int(n_diverged), tier=tier, note=note))
    print(f'  [{tier}] {claim}: {n_compared} compared, {n_diverged} diverged {note}')


def save(df, name):
    """
        Write one audit table.

        df (pd.DataFrame): the table.
        name (str): the suffix after the audit stem.

        Returns None.
    """
    df.to_csv(f'{OUT}_{name}.csv', index=False)


def contrast_arms(name):
    """
        Split 'a - b' into its two arms.

        name (str): the contrast.

        Returns tuple[str, str]: minuend and subtrahend.
    """
    a, b = name.split(' - ')
    return a, b


# ------------------------------------------------------------------------------------------------------------------
print('=== 0. execution coherence and provenance ===')
nb = json.load(open(NB))
code = [c for c in nb['cells'] if c['cell_type'] == 'code']
counts = [c.get('execution_count') for c in code]
errors = sum(1 for c in code for o in c.get('outputs', []) if o.get('output_type') == 'error')
exec_gate = pd.DataFrame([dict(n_code_cells=len(code), execution_counts=' '.join(map(str, counts)), contiguous_from_1=counts == list(range(1, len(code) + 1)), n_error_outputs=errors, n_unrun=sum(v is None for v in counts), last_cell_run=counts[-1] is not None, n_png=sum(1 for c in code for o in c.get('outputs', []) if 'image/png' in o.get('data', {})))])
print(exec_gate.T.to_string())
save(exec_gate, 'execution_gate')


def git(*args):
    """
        Run a git command in the repo.

        args (str): git arguments.

        Returns str: stdout, stripped.
    """
    return subprocess.run(['git', *args], capture_output=True, text=True, cwd=REPO).stdout.strip()


prov_paths = ['midog_utils/template_match.py', 'midog_utils/evaluate.py', 'midog_utils/nms.py', 'midog_utils/chromatin.py', 'midog_utils/channels.py', 'midog_utils/dataset.py', 'midog_utils/production.py', 'midog_utils/invariants.py', 'midog_utils/seed_selection.py', 'production_hematoxylin_only/multiotsu_pad_arms/multiotsu_template_refinement_14roi.ipynb', 'production_hematoxylin_only/multiotsu_pad_arms/multiotsurefine14_summary.csv', f'{FG}_clicks.csv', f'{FG}_per_run.csv', f'{FG}_delta_stats.csv', 'production_hematoxylin_only/window_sweep_51_61_71_49roi_6seed.ipynb', f'{SW}_per_run.csv', f'{SW}_top30.csv', f'{SW}_delta_stats.csv', 'MULTIOTSU_PAD_ARMS_PROMPT.md', NB] + [f'{MO}_{n}.csv' for n in ('templates', 'per_run', 'top30', 'per_roi', 'delta_per_roi', 'delta_stats', 'selfmatch_sensitivity')]
prov = pd.DataFrame([dict(path=p, sha256_12=hashlib.sha256(open(p, 'rb').read()).hexdigest()[:12], mtime=time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(p))), last_commit=git('log', '-1', '--format=%h %ad', '--date=format:%Y-%m-%d %H:%M', '--', p), status=git('status', '--porcelain', '--', p)[:2].strip() or 'clean') for p in prov_paths])
print(prov.to_string(index=False))
save(prov, 'provenance')

# ------------------------------------------------------------------------------------------------------------------
print('\n=== 1. inputs, parsed without midog_utils ===')
raw = json.load(open('databases/MIDOG++.json'))
IMG = pd.DataFrame([dict(image_id=im['id'], file_name=im['file_name'], width=im['width'], height=im['height'], tumor_type=im['tumor_type']) for im in raw['images']])
id2fn = dict(zip(IMG['image_id'], IMG['file_name']))
GT = pd.DataFrame([dict(ann_id=a['id'], file_name=id2fn[a['image_id']], cx=(a['bbox'][0] + a['bbox'][2]) / 2.0, cy=(a['bbox'][1] + a['bbox'][3]) / 2.0, category_id=a['category_id']) for a in raw['annotations']])
GT_BY_FILE = {fn: g.reset_index(drop=True) for fn, g in GT.groupby('file_name')}
WH = IMG.set_index('file_name')[['width', 'height']]

CLICKS = pd.read_csv(f'{FG}_clicks.csv', float_precision='round_trip')
FG_RUN = pd.read_csv(f'{FG}_per_run.csv', float_precision='round_trip')
FG_STATS = pd.read_csv(f'{FG}_delta_stats.csv', float_precision='round_trip')
SW_RUN = pd.read_csv(f'{SW}_per_run.csv', float_precision='round_trip')
SW_TOP = pd.read_csv(f'{SW}_top30.csv', float_precision='round_trip')
SW_STATS = pd.read_csv(f'{SW}_delta_stats.csv', float_precision='round_trip')
RUN = pd.read_csv(f'{MO}_per_run.csv', float_precision='round_trip')
TOP = pd.read_csv(f'{MO}_top30.csv', float_precision='round_trip')
TPL = pd.read_csv(f'{MO}_templates.csv', float_precision='round_trip')
PROI = pd.read_csv(f'{MO}_per_roi.csv', float_precision='round_trip')
DPR = pd.read_csv(f'{MO}_delta_per_roi.csv', float_precision='round_trip')
DST = pd.read_csv(f'{MO}_delta_stats.csv', float_precision='round_trip')
SMS = pd.read_csv(f'{MO}_selfmatch_sensitivity.csv', float_precision='round_trip')
ADOPT = pd.read_csv('production_hematoxylin_only/multiotsu_pad_arms/multiotsurefine14_summary.csv', float_precision='round_trip')
ROI_ORDER = list(pd.unique(FG_RUN['file_name']))
DOMAIN = CLICKS.drop_duplicates('file_name').set_index('file_name')['domain']
print(f'{len(IMG)} images, {len(GT)} annotations; ROI order from the fill-geometry per_run ({len(ROI_ORDER)} ROIs)')

# ------------------------------------------------------------------------------------------------------------------
print('\n=== 2. composition gate ===')
REF_RUN = SW_RUN[SW_RUN['condition'].isin(REFS)].copy()
REF_TOP = SW_TOP[SW_TOP['condition'].isin(REFS)].copy()
tt = IMG.set_index('file_name')['tumor_type'].replace({'canine lymphoma': 'canine lymphosarcoma'})
comp = []


def comp_row(what, expected, observed):
    """
        One composition check.

        what (str): the property.
        expected: the value the design implies.
        observed: the value in the artifact.

        Returns None.
    """
    comp.append(dict(check=what, expected=str(expected), observed=str(observed), ok=str(expected) == str(observed)))


comp_row('clicks rows (49 ROIs x 6 seeds)', 294, len(CLICKS))
comp_row('ROIs in clicks', 49, CLICKS['file_name'].nunique())
comp_row('ROIs per domain (7 domains x 7)', {d: 7 for d in sorted(CLICKS['domain'].unique())}, CLICKS.groupby('domain')['file_name'].nunique().sort_index().to_dict())
comp_row('domain label == MIDOG++.json tumor_type on every ROI', 49, int((CLICKS.drop_duplicates('file_name').set_index('file_name')['domain'] == tt.reindex(CLICKS.drop_duplicates('file_name')['file_name'])).sum()))
comp_row('tumour types in MIDOG++.json', 7, tt.nunique())
comp_row('per_run rows (294 clicks x 4 arms x 2 axes)', 2352, len(RUN))
comp_row('per_run duplicate keys', 0, int(RUN.duplicated(KEY).sum()))
comp_row('per_run ROIs', 49, RUN['file_name'].nunique())
comp_row('per_run rows per (arm, axis)', [294], sorted(RUN.groupby(['condition', 'rank_key']).size().unique().tolist()))
comp_row('top30 rows (2352 lists x 30)', 70560, len(TOP))
comp_row('top30 lists with ranks exactly 0..29', 2352, int(TOP.groupby(KEY)['rank'].apply(lambda r: sorted(r.tolist()) == list(range(30))).sum()))
comp_row('templates rows (294 clicks x 4 arms)', 1176, len(TPL))
comp_row('templates accepted', 1176, int(TPL['accepted'].sum()))
comp_row('per_roi rows (49 x 7 arms x 2 axes)', 686, len(PROI))
comp_row('delta_per_roi rows (60 specs x 49)', 2940, len(DPR))
comp_row('delta_stats rows (60 recorded + 34 reruns)', 94, len(DST))
comp_row('selfmatch_sensitivity rows (10 tests x 2 versions)', 20, len(SMS))
comp_row('reference per_run rows (294 x 3 x 2)', 1764, len(REF_RUN))
comp_row('reference per_run duplicate keys', 0, int(REF_RUN.duplicated(KEY).sum()))
comp_row('reference top30 rows', 52920, len(REF_TOP))
comp_row('clicks equal between clicks CSV and per_run (seed_ann_id, click)', 2352, int(RUN.merge(CLICKS[['file_name', 'seed_index', 'seed_ann_id', 'click_cx', 'click_cy']], on=['file_name', 'seed_index'], suffixes=('', '_c')).eval('seed_ann_id == seed_ann_id_c and click_cx == click_cx_c and click_cy == click_cy_c').sum()))
comp_row('per_run (file, seed) set == clicks CSV set', True, set(zip(RUN['file_name'], RUN['seed_index'])) == set(zip(CLICKS['file_name'], CLICKS['seed_index'])))
COMP = pd.DataFrame(comp)
print(COMP.to_string(index=False))
save(COMP, 'composition')
ledger('composition', 'row counts, keys, strata', len(COMP), int((~COMP['ok']).sum()), 'A')

# ------------------------------------------------------------------------------------------------------------------
print('\n=== 3. matching re-derived from the top 30 against MIDOG++.json ===')


def my_greedy(det_xy, gt_xy, r):
    """
        Greedy one-to-one matching, detections in rank order, each claiming its nearest free GT within r.

        det_xy (np.ndarray): (n, 2) ranked detections.
        gt_xy (np.ndarray): (m, 2) ground truth.
        r (float): match radius in px (inclusive).

        Returns np.ndarray: GT index per detection, -1 when none.
    """
    taken = np.zeros(len(gt_xy), bool)
    out = np.full(len(det_xy), -1, int)
    for i, (x, y) in enumerate(det_xy):
        d = np.hypot(gt_xy[:, 0] - x, gt_xy[:, 1] - y)
        cand = np.where((d <= r) & ~taken)[0]
        if len(cand):
            g = cand[np.lexsort((cand, d[cand]))[0]]
            out[i] = g
            taken[g] = True
    return out


MPP = CLICKS.drop_duplicates('file_name').set_index('file_name')['mpp']
SEED_ANN = CLICKS.set_index(['file_name', 'seed_index'])['seed_ann_id']
ALLTOP = pd.concat([TOP.assign(source='this run'), REF_TOP.assign(source='window sweep')], ignore_index=True).sort_values(KEY + ['rank'], kind='stable').reset_index(drop=True)
buckets, anns = np.empty(len(ALLTOP), object), np.full(len(ALLTOP), -1, int)
for (fn, s), idx in ALLTOP.groupby(['file_name', 'seed_index']).groups.items():
    g = GT_BY_FILE[fn]
    g = g[g['ann_id'] != SEED_ANN[(fn, s)]].reset_index(drop=True)
    gxy, gcat, gid = g[['cx', 'cy']].to_numpy(float), g['category_id'].to_numpy(), g['ann_id'].to_numpy()
    r = 7.5 / MPP[fn]
    sub = ALLTOP.loc[idx]
    for _, lst in sub.groupby(['condition', 'rank_key']):
        lst = lst.sort_values('rank')
        m = my_greedy(lst[['cx', 'cy']].to_numpy(float), gxy, r)
        buckets[lst.index] = np.where(m < 0, 'non_human_findings', np.where(gcat[np.maximum(m, 0)] == 1, 'human_correct_label', 'human_rejected_label'))
        anns[lst.index] = np.where(m < 0, -1, gid[np.maximum(m, 0)])
ALLTOP['my_bucket'], ALLTOP['my_ann'] = buckets, anns
n_b = int((ALLTOP['my_bucket'] != ALLTOP['bucket']).sum())
n_a = int((ALLTOP['my_ann'] != ALLTOP['matched_ann_id']).sum())
ledger('matching', 'top-30 bucket per detection (4 arms + 3 references)', len(ALLTOP), n_b, 'A')
ledger('matching', 'top-30 matched_ann_id per detection', len(ALLTOP), n_a, 'A')
ALLTOP['hit'] = ALLTOP['my_bucket'] == 'human_correct_label'
my_tp = pd.concat({f'tp_at_{k}': ALLTOP[ALLTOP['rank'] < k].groupby(KEY)['hit'].sum() for k in (10, 20, 30)}, axis=1)
both_run = pd.concat([RUN.assign(source='this run'), REF_RUN.assign(source='window sweep')], ignore_index=True).set_index(KEY)
cmp_tp = both_run[['tp_at_10', 'tp_at_20', 'tp_at_30']].join(my_tp, rsuffix='_mine')
n_div = sum(int((cmp_tp[c] != cmp_tp[f'{c}_mine']).sum()) for c in ('tp_at_10', 'tp_at_20', 'tp_at_30'))
ledger('matching', 'TP@10/20/30 per list from own matching (2352 + 1764 lists x 3)', 3 * len(cmp_tp), n_div, 'A')
save(pd.DataFrame([dict(rows=len(ALLTOP), bucket_div=n_b, ann_div=n_a, tp_values=3 * len(cmp_tp), tp_div=n_div, lookalike_top10_slots_this_run=int(((ALLTOP['my_bucket'] == 'human_rejected_label') & (ALLTOP['rank'] < 10) & (ALLTOP['source'] == 'this run')).sum()))]), 'matching')

# ------------------------------------------------------------------------------------------------------------------
print('\n=== 4. reference provenance: the window-sweep references against the committed fill-geometry run ===')
fg = FG_RUN.set_index(KEY)
rp = []
for sw_arm, fg_arm in (('default_51', 'default_51'), ('fill_51', 'hem_rect_fill')):
    a = REF_RUN[REF_RUN['condition'] == sw_arm].set_index(['file_name', 'seed_index', 'rank_key'])
    b = FG_RUN[FG_RUN['condition'] == fg_arm].set_index(['file_name', 'seed_index', 'rank_key'])
    j = a.join(b, rsuffix='_fg', how='inner')
    for c in ('tp_at_10', 'tp_at_20', 'tp_at_30', 'n_detections', 'n_peaks', 'tpl_sha1', 'n_blanked_px'):
        rp.append(dict(window_sweep_arm=sw_arm, fillgeom_arm=fg_arm, column=c, n=len(j), n_diverged=int((j[c] != j[f'{c}_fg']).sum())))
sha = CLICKS.set_index(['file_name', 'seed_index'])
for arm, col in (('default_51', 'sha_default_51'), ('fill_51', 'sha_hem_rect_fill')):
    a = REF_RUN[(REF_RUN['condition'] == arm)]
    rp.append(dict(window_sweep_arm=arm, fillgeom_arm=f'clicks CSV {col}', column='tpl_sha1', n=len(a), n_diverged=int((a['tpl_sha1'].to_numpy() != sha.loc[list(zip(a['file_name'], a['seed_index'])), col].to_numpy()).sum())))
RP = pd.DataFrame(rp)
print(RP.to_string(index=False))
save(RP, 'reference_provenance')
ledger('provenance', 'window-sweep default_51 and fill_51 == committed fill-geometry default_51 / hem_rect_fill', int(RP['n'].sum()), int(RP['n_diverged'].sum()), 'A')

# ------------------------------------------------------------------------------------------------------------------
print('\n=== 5. per-ROI sums, pooled precision ===')
ALLRUN = both_run.reset_index()
ALLRUN['tp_at_10_excl_border'] = ALLRUN['tp_at_10_excl_border'].astype(int)
roi = ALLRUN.groupby(['rank_key', 'condition', 'file_name'])[list(METRIC_K)].sum()
nclk = ALLRUN.groupby(['rank_key', 'condition', 'file_name'])['seed_index'].nunique()
pr = PROI.set_index(['rank_key', 'condition', 'file_name'])
n_d = sum(int((pr[m] != roi.reindex(pr.index)[m]).sum()) for m in METRIC_K) + int((pr['n_clicks'] != nclk.reindex(pr.index)).sum())
ledger('per_roi', 'per_roi TP sums and n_clicks (686 rows x 5)', 5 * len(pr), n_d, 'A')
pooled = ALLRUN.groupby(['rank_key', 'condition'])[list(METRIC_K)].sum()
for m, k in METRIC_K.items():
    pooled[m] = pooled[m] / (294 * k)
pooled = pooled.rename(columns=lambda m: m.replace('tp_at', 'precision_at')).reset_index()
print(pooled.to_string(index=False, float_format=lambda v: f'{v:.4f}'))
save(pooled, 'pooled_precision')
QUOTED_POOLED = {('tm_score', 'mo_box'): 0.4602, ('tm_score', 'mo_ring5'): 0.5102, ('tm_score', 'mo_pad51'): 0.5248, ('tm_score', 'mo_pad61'): 0.5418, ('tm_score', 'fill_51'): 0.5265, ('tm_score', 'fill_61'): 0.5425}
pp = pooled.set_index(['rank_key', 'condition'])['precision_at_10']
ledger('pooled', 'prose pooled precision@10 values (cell 28)', len(QUOTED_POOLED), sum(round(pp[k], 4) != v for k, v in QUOTED_POOLED.items()), 'A')
# the printed cell-19 table: 7 arms x 2 axes x 4 metrics, at 4 decimals
PRINTED_19 = {'chromatin_od': {'mo_box': (0.5429, 0.5578, 0.4662, 0.4019), 'mo_pad51': (0.5704, 0.5915, 0.5099, 0.4437), 'mo_pad61': (0.5626, 0.5847, 0.5090, 0.4430), 'mo_ring5': (0.5755, 0.5908, 0.5058, 0.4397), 'default_51': (0.5667, 0.5844, 0.5029, 0.4351), 'fill_51': (0.5680, 0.5861, 0.5073, 0.4437), 'fill_61': (0.5622, 0.5827, 0.5066, 0.4426)}, 'tm_score': {'mo_box': (0.4602, 0.4612, 0.3855, 0.3418), 'mo_pad51': (0.5248, 0.5235, 0.4571, 0.4086), 'mo_pad61': (0.5418, 0.5429, 0.4685, 0.4162), 'mo_ring5': (0.5102, 0.5088, 0.4366, 0.3908), 'default_51': (0.5085, 0.5065, 0.4361, 0.3906), 'fill_51': (0.5265, 0.5276, 0.4595, 0.4096), 'fill_61': (0.5425, 0.5439, 0.4724, 0.4188)}}
pv = pooled.set_index(['rank_key', 'condition'])
n19 = 0
for ax, d in PRINTED_19.items():
    for arm, vals in d.items():
        for col, v in zip(('precision_at_10', 'precision_at_10_excl_border', 'precision_at_20', 'precision_at_30'), vals):
            n19 += round(pv.loc[(ax, arm), col], 4) != v
ledger('pooled', 'cell-19 pooled table (7 arms x 2 axes x 4 metrics)', 56, n19, 'A')

# ------------------------------------------------------------------------------------------------------------------
print('\n=== 6. statistics, re-implemented ===')


def exact_signflip_p(d):
    """
        Two-sided exact sign-flip p over the nonzero integer deltas, by counting sign patterns with a DP over sums.

        d (np.ndarray): integer deltas, one per unit.

        Returns tuple[float, int]: p and the number of nonzero units.
    """
    nz = [abs(int(v)) for v in d if int(v) != 0]
    if not nz:
        return 1.0, 0
    ways = {0: 1}
    for v in nz:
        nxt = {}
        for s, c in ways.items():
            nxt[s + v] = nxt.get(s + v, 0) + c
            nxt[s - v] = nxt.get(s - v, 0) + c
        ways = nxt
    obs = abs(int(np.sum(d)))
    return sum(c for s, c in ways.items() if abs(s) >= obs) / 2 ** len(nz), len(nz)


def holm(p):
    """
        Holm step-down adjustment.

        p (list[float]): raw p-values.

        Returns np.ndarray: adjusted p-values in input order.
    """
    p = np.asarray(p, float)
    o = np.argsort(p, kind='stable')
    adj = np.minimum(1.0, np.maximum.accumulate([(len(p) - i) * p[j] for i, j in enumerate(o)]))
    out = np.empty(len(p))
    out[o] = adj
    return out


def my_stats(d, n, k, stream=0):
    """
        Every statistic of one contrast from per-ROI integer deltas and per-ROI click counts.

        d (np.ndarray): integer delta per ROI, in ROI order.
        n (np.ndarray): included clicks per ROI.
        k (int): the budget.
        stream (int): bootstrap stream (0 = the reported one).

        Returns dict: pooled delta, W/L/T, exact p, bootstrap CI (same RNG, a consistency check), t CI, D5.
    """
    d, n = np.asarray(d, np.int64), np.asarray(n, np.int64)
    G = len(d)
    pooled = d.sum() / (k * n.sum())
    w, l, t = int((d > 0).sum()), int((d < 0).sum()), int((d == 0).sum())
    p, nnz = exact_signflip_p(d)
    idx = np.random.default_rng([BOOT_SEED, stream]).integers(G, size=(N_BOOT, G))
    lo, hi = np.percentile(d[idx].sum(1) / (k * n[idx].sum(1)), [2.5, 97.5])
    x = d / (k * n)
    half = sps.t.ppf(0.975, G - 1) * x.std(ddof=1) / np.sqrt(G)
    same = w if pooled > 0 else (l if pooled < 0 else 0)
    ci_ex = bool(lo > 0 or hi < 0)
    return dict(n_roi=G, n_clicks=int(n.sum()), pooled_delta_precision=pooled, wins=w, losses=l, ties=t, n_nonzero=nnz, exact_p=p, min_attainable_p=2.0 ** (1 - nnz) if nnz else 1.0, ci_low=lo, ci_high=hi, t_ci_low=x.mean() - half, t_ci_high=x.mean() + half, ci_excludes_0=ci_ex, t_ci_excludes_0=bool(x.mean() - half > 0 or x.mean() + half < 0), n_same_direction=same, d5_bar_met=bool(ci_ex and same > G / 2))


CLICK = ALLRUN.set_index(['rank_key', 'condition', 'file_name', 'seed_index']).sort_index()


def roi_deltas(axis, metric, name, exclude=frozenset()):
    """
        Per-ROI deltas of one contrast, summed over each ROI's clicks not excluded.

        axis (str): ranking axis.
        metric (str): click-level TP column.
        name (str): 'a - b'.
        exclude (frozenset): (file_name, seed_index) clicks left out.

        Returns tuple[np.ndarray, np.ndarray, list]: deltas, click counts, ROI names (ROIs with no click left drop).
    """
    a, b = contrast_arms(name)
    ds_, ns_, fs_ = [], [], []
    for fn in ROI_ORDER:
        seeds = [s for s in range(6) if (fn, s) not in exclude]
        if not seeds:
            continue
        va = CLICK.loc[(axis, a, fn), metric].reindex(seeds)
        vb = CLICK.loc[(axis, b, fn), metric].reindex(seeds)
        ds_.append(int(va.sum() - vb.sum()))
        ns_.append(len(seeds))
        fs_.append(fn)
    return np.array(ds_), np.array(ns_), fs_


# 6a. the fill-geometry validation rows (cell 8), from the committed per_run
fgr = FG_RUN.groupby(['rank_key', 'condition', 'file_name'])['tp_at_10'].sum()
val = []
for r in FG_STATS.itertuples():
    a, b = contrast_arms(r.pair)
    d = np.array([fgr[(r.rank_key, a, f)] - fgr[(r.rank_key, b, f)] for f in ROI_ORDER])
    st = my_stats(d, np.full(49, 6), int(r.K))
    val.append(dict(rank_key=r.rank_key, pair=r.pair, **{c: st[c] for c in ('pooled_delta_precision', 'exact_p', 'ci_low', 'ci_high', 't_ci_low', 't_ci_high', 'wins', 'losses', 'ties', 'd5_bar_met')}, **{f'{c}_rec': getattr(r, c) for c in ('pooled_delta_precision', 'exact_p', 'ci_low', 'ci_high', 't_ci_low', 't_ci_high', 'wins', 'losses', 'ties', 'd5_bar_met')}))
VAL = pd.DataFrame(val)
VAL['holm_p'] = VAL.groupby(FG_STATS['holm_family'].to_numpy())['exact_p'].transform(lambda p: holm(p.to_numpy()))
VAL['holm_p_rec'] = FG_STATS['holm_p'].to_numpy()
vcols = ('pooled_delta_precision', 'exact_p', 'ci_low', 'ci_high', 't_ci_low', 't_ci_high', 'holm_p')
nv = sum(int((~np.isclose(VAL[c], VAL[f'{c}_rec'], rtol=1e-9, atol=1e-12)).sum()) for c in vcols) + sum(int((VAL[c] != VAL[f'{c}_rec']).sum()) for c in ('wins', 'losses', 'ties', 'd5_bar_met'))
ledger('statistics', "cell-8 claim: the 10 fill-geometry delta_stats rows reproduce (own stats code)", 11 * len(VAL), nv, 'A', '(bootstrap columns: consistency, same RNG)')
save(VAL, 'fillgeom_validation')

# 6b. all 60 recorded rows and 34 flag-excluded reruns
TFLAG = TPL.copy()
TPL_BOX = TPL[TPL['condition'] == 'mo_box'].set_index(['file_name', 'seed_index'])
TFLAG['click_ix'] = [int(round(v)) for v in TFLAG['click_cx']]
TFLAG['click_iy'] = [int(round(v)) for v in TFLAG['click_cy']]
TFLAG['off'] = np.maximum((TFLAG['tcx'] - TFLAG['click_ix']).abs(), (TFLAG['tcy'] - TFLAG['click_iy']).abs())
TFLAG['my_box_outside_blank'] = TFLAG['off'] + (TFLAG['bs'] - 1) // 2 > HALF_BLANK
TFLAG['my_template_outside_window'] = TFLAG['off'] + (TFLAG['S'] - 1) // 2 > HALF_BLANK
TFLAG['my_S_expected'] = [bs if arm == 'mo_box' else (bs + 10 if arm == 'mo_ring5' else (bs if bs >= int(arm[-2:]) else int(arm[-2:]))) for arm, bs in zip(TFLAG['condition'], TFLAG['bs'])]
TFLAG['my_box_ge_S'] = TFLAG['condition'].isin(['mo_pad51', 'mo_pad61']) & (TFLAG['bs'] >= TFLAG['S_nominal'])
TFLAG['my_mass_weight_nonpositive'] = (TFLAG['n_filled_px'] > 0) & (TFLAG['mass_weight'] <= 0)
TFLAG['my_n_filled'] = TFLAG['S'] ** 2 - np.where(TFLAG['S'] == TFLAG['bs'], TFLAG['S'] ** 2, TFLAG['bs'] ** 2)
nf = sum(int((TFLAG[f'my_{f}'] != TFLAG[f]).sum()) for f in ('box_outside_blank', 'template_outside_window', 'box_ge_S', 'mass_weight_nonpositive')) + int((TFLAG['my_S_expected'] != TFLAG['S']).sum()) + int((TFLAG['my_n_filled'] != TFLAG['n_filled_px']).sum()) + int((TFLAG['off'] != TFLAG['offset_cheb_px']).sum())
ledger('flags', 'flags, S, n_filled and Chebyshev offset re-derived from box geometry (1176 x 7)', 7 * len(TFLAG), nf, 'A')
FLAG_COUNTS = TFLAG.groupby('condition')[['my_box_outside_blank', 'my_template_outside_window', 'my_box_ge_S', 'my_mass_weight_nonpositive']].sum().reindex(ARMS)
print(FLAG_COUNTS.to_string())
QUOTED_FLAGS = {('mo_box', 'my_box_outside_blank'): 14, ('mo_pad51', 'my_box_outside_blank'): 14, ('mo_pad61', 'my_box_outside_blank'): 14, ('mo_ring5', 'my_box_outside_blank'): 14, ('mo_box', 'my_template_outside_window'): 14, ('mo_pad51', 'my_template_outside_window'): 276, ('mo_pad61', 'my_template_outside_window'): 294, ('mo_ring5', 'my_template_outside_window'): 85, ('mo_pad51', 'my_box_ge_S'): 6, ('mo_pad61', 'my_box_ge_S'): 0, ('mo_pad61', 'my_mass_weight_nonpositive'): 3, ('mo_pad51', 'my_mass_weight_nonpositive'): 0, ('mo_ring5', 'my_mass_weight_nonpositive'): 0}
ledger('flags', 'prose flag table (cell 28)', len(QUOTED_FLAGS), sum(int(FLAG_COUNTS.loc[a, f]) != v for (a, f), v in QUOTED_FLAGS.items()), 'A')
save(TFLAG.groupby('condition')[['my_box_outside_blank', 'my_template_outside_window', 'my_box_ge_S', 'my_mass_weight_nonpositive']].sum().reset_index(), 'flag_counts')
FLAG_OF = {f: set(zip(TFLAG.loc[TFLAG[f'my_{f}'], 'file_name'], TFLAG.loc[TFLAG[f'my_{f}'], 'seed_index'], TFLAG.loc[TFLAG[f'my_{f}'], 'condition'])) for f in ('box_outside_blank', 'template_outside_window', 'box_ge_S', 'mass_weight_nonpositive')}
FLAG_OF['roi_edge_read'] = set()
FLAG_SETS = {'box_outside_blank': ('box_outside_blank',), 'template_outside_window': ('template_outside_window',), 'roi_edge_read': ('roi_edge_read',), 'box_ge_S': ('box_ge_S',), 'mass_weight_nonpositive': ('mass_weight_nonpositive',), 'any_non_informational': ('box_outside_blank', 'roi_edge_read', 'box_ge_S', 'mass_weight_nonpositive')}

rows = []
for axis in AXES:
    for name in list(TESTED) + DESCRIPTIVE:
        d, n, _ = roi_deltas(axis, 'tp_at_10', name)
        rows.append(dict(rank_key=axis, metric='tp_at_10', contrast=name, version='recorded', **my_stats(d, n, 10)))
    for m in ('tp_at_10_excl_border', 'tp_at_20', 'tp_at_30'):
        for name in TESTED:
            d, n, _ = roi_deltas(axis, m, name)
            rows.append(dict(rank_key=axis, metric=m, contrast=name, version='recorded', **my_stats(d, n, METRIC_K[m])))
for label_, fset in FLAG_SETS.items():
    for axis in AXES:
        for name in TESTED:
            arms = [a for a in contrast_arms(name) if a in ARMS]
            excl = frozenset((fn, int(s)) for f in fset for (fn, s, c) in FLAG_OF[f] if c in arms)
            if not excl:
                continue
            d, n, fs_ = roi_deltas(axis, 'tp_at_10', name, excl)
            if not len(d):
                continue
            rows.append(dict(rank_key=axis, metric='tp_at_10', contrast=name, version=f'excl {label_}', n_clicks_excluded=len(excl), **my_stats(d, n, 10)))
MINE = pd.DataFrame(rows)
tested_mask = (MINE['metric'] == 'tp_at_10') & MINE['contrast'].isin(TESTED) & (MINE['version'] == 'recorded')
MINE['holm_family'] = np.where(tested_mask, MINE['rank_key'] + ' | ' + MINE['contrast'].map(TESTED).fillna(''), '')
MINE['holm_p'] = np.nan
for fam in MINE.loc[tested_mask, 'holm_family'].unique():
    sel = MINE['holm_family'] == fam
    MINE.loc[sel, 'holm_p'] = holm(MINE.loc[sel, 'exact_p'])
j = MINE.merge(DST, on=['rank_key', 'metric', 'contrast', 'version'], suffixes=('', '_rec'), how='outer', indicator=True)
assert (j['_merge'] == 'both').all(), j.loc[j['_merge'] != 'both', ['rank_key', 'metric', 'contrast', 'version', '_merge']]
float_cols = ['pooled_delta_precision', 'exact_p', 'min_attainable_p', 't_ci_low', 't_ci_high']
boot_cols = ['ci_low', 'ci_high']
int_cols = ['n_roi', 'n_clicks', 'wins', 'losses', 'ties', 'n_nonzero', 'n_same_direction', 'ci_excludes_0', 't_ci_excludes_0', 'd5_bar_met']
div_ind = sum(int((~np.isclose(j[c], j[f'{c}_rec'], rtol=1e-9, atol=1e-12)).sum()) for c in float_cols) + sum(int((j[c] != j[f'{c}_rec']).sum()) for c in int_cols)
div_holm = int((~np.isclose(j.loc[j['holm_family'] != '', 'holm_p'], j.loc[j['holm_family'] != '', 'holm_p_rec'], rtol=1e-9)).sum())
div_boot = sum(int((~np.isclose(j[c], j[f'{c}_rec'], rtol=1e-9, atol=1e-12)).sum()) for c in boot_cols)
ledger('statistics', 'delta_stats: pooled delta, exact p, floor, t CI, W/L/T, n_nonzero, D5 (94 rows)', len(j) * (len(float_cols) + len(int_cols)), div_ind, 'A')
ledger('statistics', 'delta_stats: Holm p on the 10 tested rows', 10, div_holm, 'A')
ledger('statistics', 'delta_stats: bootstrap CI bounds (same BOOT_SEED)', len(j) * 2, div_boot, 'consistency')
save(j.drop(columns='_merge'), 'delta_stats_recomputed')

# 6c. delta_per_roi
dpr_n = 0
for (axis, metric, name), g in DPR.groupby(['rank_key', 'metric', 'contrast']):
    d, n, fs_ = roi_deltas(axis, metric, name)
    gg = g.set_index('file_name').reindex(fs_)
    dpr_n += int((gg['delta_tp_sum'].to_numpy() != d).sum())
ledger('statistics', 'delta_per_roi delta_tp_sum (2940 rows)', len(DPR), dpr_n, 'A')

# 6d. prose table of the 10 tests (cell 28), at the printed precision
PROSE = {('tm_score', 'mo_pad51 - mo_box'): (6.46, 4.35, 8.67, 4.20, 8.72, '38/7/4'), ('tm_score', 'mo_pad61 - mo_box'): (8.16, 5.71, 10.65, 5.57, 10.76, '36/8/5'), ('tm_score', 'mo_ring5 - mo_box'): (5.00, 3.37, 6.67, 3.26, 6.74, '37/10/2'), ('tm_score', 'mo_pad51 - fill_51'): (-0.17, -0.95, 0.65, -0.99, 0.65, '16/19/14'), ('tm_score', 'mo_pad61 - fill_61'): (-0.07, -0.82, 0.68, -0.84, 0.70, '15/17/17'), ('chromatin_od', 'mo_pad51 - mo_box'): (2.76, 0.95, 4.56, 0.86, 4.65, '26/8/15'), ('chromatin_od', 'mo_pad61 - mo_box'): (1.97, -0.31, 3.98, -0.28, 4.23, '27/11/11'), ('chromatin_od', 'mo_ring5 - mo_box'): (3.27, 1.90, 4.66, 1.82, 4.71, '27/8/14'), ('chromatin_od', 'mo_pad51 - fill_51'): (0.24, -0.48, 1.02, -0.54, 1.01, '13/9/27'), ('chromatin_od', 'mo_pad61 - fill_61'): (0.03, -0.61, 0.82, -0.71, 0.78, '10/11/28')}
mt = MINE[tested_mask].set_index(['rank_key', 'contrast'])
npd = 0
for key, (dl, lo, hi, tlo, thi, wlt) in PROSE.items():
    r = mt.loc[key]
    npd += sum(round(100 * r[c], 2) != v for c, v in zip(('pooled_delta_precision', 'ci_low', 'ci_high', 't_ci_low', 't_ci_high'), (dl, lo, hi, tlo, thi)))
    npd += f"{r['wins']}/{r['losses']}/{r['ties']}" != wlt
ledger('prose', 'cell-28 table of the 10 tests (delta, boot CI, t CI, W/L/T)', 60, npd, 'A')

# 6e. bootstrap re-seed stability (cell 21): same streams 1..20 -- consistency
stab = []
for r in MINE[tested_mask & MINE['ci_excludes_0']].itertuples():
    d, n, _ = roi_deltas(r.rank_key, 'tp_at_10', r.contrast)
    bnds = []
    for s in range(20):
        idx = np.random.default_rng([BOOT_SEED, 1 + s]).integers(len(d), size=(N_BOOT, len(d)))
        bnds.append(np.percentile(d[idx].sum(1) / (10 * n[idx].sum(1)), [2.5, 97.5]))
    stab.append(dict(rank_key=r.rank_key, contrast=r.contrast, all_exclude_0=all(lo > 0 or hi < 0 for lo, hi in bnds), closest_bound_pts=round(100 * min(abs(lo) if r.pooled_delta_precision > 0 else abs(hi) for lo, hi in bnds), 3)))
STAB = pd.DataFrame(stab)
print(STAB.to_string(index=False))
save(STAB, 'reseed_stability')
ledger('statistics', 'cell-21 re-seed stability: 5 CIs exclude 0 on all 20 streams; closest bound 0.850', 6, int((~STAB['all_exclude_0']).sum()) + int(STAB['closest_bound_pts'].min() != 0.85), 'consistency')

# 6f. the fill_61 - fill_51 join check (cell 21) against the window sweep's saved row
jc = 0
for axis in AXES:
    mine = MINE[(MINE['rank_key'] == axis) & (MINE['metric'] == 'tp_at_10') & (MINE['contrast'] == 'fill_61 - fill_51') & (MINE['version'] == 'recorded')].iloc[0]
    saved = SW_STATS[(SW_STATS['rank_key'] == axis) & (SW_STATS['metric'] == 'tp_at_10') & (SW_STATS['contrast'] == 'fill_61 - fill_51')].iloc[0]
    jc += sum(not np.isclose(mine[c], saved[c], rtol=1e-12) for c in ('pooled_delta_precision', 'exact_p', 't_ci_low', 't_ci_high')) + sum(mine[c] != saved[c] for c in ('wins', 'losses', 'ties'))
ledger('statistics', 'fill_61 - fill_51 equals the window sweep saved row (both axes)', 14, jc, 'A')

# 6g. the prose's specific rerun, floor and mass-weight numbers (cell 28)
ms = MINE.set_index(['rank_key', 'metric', 'contrast', 'version'])
QUOTED_RERUN = {('tm_score', 'mo_pad61 - mo_box', 'excl any_non_informational'): (8.41, 5.82, 11.03), ('tm_score', 'mo_pad51 - fill_51', 'excl any_non_informational'): (-0.11, -0.93, 0.72), ('chromatin_od', 'mo_pad61 - mo_box', 'excl any_non_informational'): (2.09, -0.33, 4.23), ('tm_score', 'mo_ring5 - mo_box', 'excl template_outside_window'): (4.88, 2.86, 7.00), ('chromatin_od', 'mo_ring5 - mo_box', 'excl template_outside_window'): (4.21, 2.38, 6.19), ('chromatin_od', 'mo_pad61 - fill_61', 'excl mass_weight_nonpositive'): (-0.10, None, None), ('chromatin_od', 'mo_pad61 - fill_61', 'excl any_non_informational'): (-0.07, None, None), ('chromatin_od', 'mo_pad61 - mo_box', 'recorded_excl_border'): (2.69, 0.85, 4.52)}
nq = 0
for (axis, name, ver), vals in QUOTED_RERUN.items():
    r = ms.loc[(axis, 'tp_at_10_excl_border', name, 'recorded')] if ver == 'recorded_excl_border' else ms.loc[(axis, 'tp_at_10', name, ver)]
    nq += sum(v is not None and round(100 * r[c], 2) != v for c, v in zip(('pooled_delta_precision', 'ci_low', 'ci_high'), vals))
tow = ms.loc[('tm_score', 'tp_at_10', 'mo_pad51 - mo_box', 'excl template_outside_window')]
nq += int(tow['n_clicks'] != 18) + int(tow['n_roi'] != 15)
ledger('prose', 'cell-28 rerun numbers (union, template_outside_window, mass-weight, edge band) and 18 clicks / 15 ROIs', 22, nq, 'A')
rec_sign = MINE[tested_mask].set_index(['rank_key', 'contrast'])['pooled_delta_precision']
ni = MINE[MINE['version'].isin([f'excl {f}' for f in ('box_outside_blank', 'roi_edge_read', 'box_ge_S', 'mass_weight_nonpositive', 'any_non_informational')])]
sign_changes = ni[np.sign(ni['pooled_delta_precision'].to_numpy()) != np.sign(rec_sign.reindex(list(zip(ni['rank_key'], ni['contrast']))).to_numpy())][['rank_key', 'contrast', 'version', 'pooled_delta_precision']]
d5_changes = ni[ni['d5_bar_met'].to_numpy() != MINE[tested_mask].set_index(['rank_key', 'contrast'])['d5_bar_met'].reindex(list(zip(ni['rank_key'], ni['contrast']))).to_numpy()]
print('non-informational reruns with a sign change:\n', sign_changes.to_string(index=False))
ledger('prose', "'excluding any non-informational flag leaves every D5 verdict unchanged; only one sign changes (chromatin_od mo_pad61 - fill_61)'", len(ni), len(d5_changes) + int(set(zip(sign_changes['rank_key'], sign_changes['contrast'])) != {('chromatin_od', 'mo_pad61 - fill_61')}), 'A')
floor = (MINE[tested_mask]['min_attainable_p'] * MINE[tested_mask]['contrast'].map(TESTED).map({'pad': 3, 'vs_incumbent': 2})).max()
ledger('prose', 'largest Holm-adjusted floor 1.91e-06 (no tested row floor-blocked)', 1, int(f'{floor:.3g}' != '1.91e-06'), 'A')
mwq = TPL[TPL['n_filled_px'] > 0].groupby('condition')['mass_weight'].min()
nm = int(round(mwq['mo_pad51'], 5) != 0.00053) + int(round(mwq['mo_pad61'], 5) != -0.00503) + int(round(mwq['mo_ring5'], 5) != 0.00445)
neg = TPL[(TPL['condition'] == 'mo_pad61') & (TPL['mass_weight'] <= 0)].set_index(['file_name', 'seed_index'])
nm += int(sorted(neg.index.tolist()) != sorted([('301.tiff', 2), ('245.tiff', 2), ('284.tiff', 3)]))
nm += sum(round(neg.loc[k, 'mass_weight'], 4) != v for k, v in ((('301.tiff', 2), -0.0050), (('245.tiff', 2), -0.0004), (('284.tiff', 3), -0.0012)))
nm += int(not neg['bs'].between(25, 33).all()) + int(neg['box_outside_blank'].any())
ledger('prose', 'minimum mass weights 0.00053 / -0.00503 / 0.00445; the three non-positive clicks, values, bs 25-33, none box_outside_blank', 9, nm, 'A')
nbr = 0
for r in CLICKS.itertuples():
    g = GT_BY_FILE[r.file_name]
    g = g[(g['ann_id'] != r.seed_ann_id) & (g['category_id'] == 1)]
    rad = 7.5 / r.mpp
    tc = TPL_BOX.loc[(r.file_name, r.seed_index)]
    dc = np.hypot(g['cx'] - r.click_cx, g['cy'] - r.click_cy)
    dt = np.hypot(g['cx'] - tc['tcx'], g['cy'] - tc['tcy'])
    nbr += int(((dc <= 2 * rad + 3) | (dt <= 2 * rad + 3)).any())
ledger('prose', 'cell-25: 12 clicks have another mitosis within 2r+3 of the click or box centre', 1, int(nbr != 12), 'A')

# ------------------------------------------------------------------------------------------------------------------
print('\n=== 7. diagnostics re-derived from the top 30 ===')
T10 = ALLTOP[ALLTOP['rank'] < 10].copy()
cen = TPL[['file_name', 'seed_index', 'condition', 'tcx', 'tcy', 'bs', 'S', 'click_cx', 'click_cy']]
T30 = ALLTOP.merge(CLICKS[['file_name', 'seed_index', 'click_cx', 'click_cy']], on=['file_name', 'seed_index']).merge(cen.drop(columns=['click_cx', 'click_cy']), on=['file_name', 'seed_index', 'condition'], how='left')
T30['tcx'] = T30['tcx'].fillna(T30['click_cx'])
T30['tcy'] = T30['tcy'].fillna(T30['click_cy'])
T30['S'] = T30['S'].fillna(T30['condition'].map({'default_51': 51, 'fill_51': 51, 'fill_61': 61}))
nt = ~T30['hit']
cheb_click = np.maximum((T30['cx'] - T30['click_cx']).abs(), (T30['cy'] - T30['click_cy']).abs())
cheb_cent = np.maximum((T30['cx'] - T30['tcx']).abs(), (T30['cy'] - T30['tcy']).abs())
T30['self_match'] = nt & ((cheb_click <= SELF_PX) | (cheb_cent <= SELF_PX))
ix = np.array([int(round(v)) for v in T30['click_cx']])
iy = np.array([int(round(v)) for v in T30['click_cy']])
off = np.maximum((T30['tcx'] - ix).abs(), (T30['tcy'] - iy).abs())
T30['reach'] = HALF_BLANK + (T30['S'].astype(int) - 1) // 2 + off.astype(int)
T30['near_reach'] = np.maximum(np.abs(T30['cx'] - ix), np.abs(T30['cy'] - iy)) <= T30['reach']
T30['near60c'] = cheb_click <= NEAR_PX
T30['near60e'] = np.hypot(T30['cx'] - T30['click_cx'], T30['cy'] - T30['click_cy']) <= NEAR_PX
w_ = T30['file_name'].map(WH['width'])
h_ = T30['file_name'].map(WH['height'])
T30['band'] = np.minimum.reduce([T30['cx'], T30['cy'], w_ - 1 - T30['cx'], h_ - 1 - T30['cy']]) < EDGE_BAND
top10 = T30['rank'] < 10
diag = T30[top10].groupby(KEY).agg(sm10=('self_match', 'sum'), nr=('near_reach', 'sum'), c60=('near60c', 'sum'), e60=('near60e', 'sum'), band=('band', 'sum'))
diag['sm30'] = T30.groupby(KEY)['self_match'].sum()
diag['nr_fp'] = T30[top10 & nt].groupby(KEY)['near_reach'].sum().reindex(diag.index, fill_value=0)
diag['c60_fp'] = T30[top10 & nt].groupby(KEY)['near60c'].sum().reindex(diag.index, fill_value=0)
diag['e60_fp'] = T30[top10 & nt].groupby(KEY)['near60e'].sum().reindex(diag.index, fill_value=0)
rr = RUN.set_index(KEY)
dj = rr.join(diag)
pairs_ = [('n_selfmatch_top10', 'sm10'), ('n_selfmatch_top30', 'sm30'), ('n_near_seed_top10_reach', 'nr'), ('n_near_seed_fp_top10_reach', 'nr_fp'), ('n_near_seed_top10_cheb', 'c60'), ('n_near_seed_fp_top10_cheb', 'c60_fp'), ('n_near_seed_top10_euclid', 'e60'), ('n_near_seed_fp_top10_euclid', 'e60_fp'), ('n_border_band_top10', 'band'), ('near_seed_reach_px', None)]
ndg = 0
for a, b in pairs_:
    if b is None:
        continue
    ndg += int((dj[a] != dj[b]).sum())
ledger('diagnostics', 'self-match, near-seed (reach, 60 px Cheb/Euclid) and edge-band counts per list (2352 x 9)', 9 * len(dj), ndg, 'A')
# reach itself
rj = rr.join(T30[T30['source'] == 'this run'].groupby(KEY)['reach'].first())
ledger('diagnostics', 'near_seed_reach_px = 25 + (S-1)//2 + off', len(rj), int((rj['near_seed_reach_px'] != rj['reach']).sum()), 'A')
# references
ref_diag = T30[T30['source'] == 'window sweep']
print('reference self-match top30:', int(ref_diag['self_match'].sum()), '; reference near-seed 60px Cheb non-TP top10:', int((ref_diag['near60c'] & ~ref_diag['hit'] & (ref_diag['rank'] < 10)).sum()))
# the whole-list self-match count only exists as a column; check it is 0 and is consistent with the top 30
print('n_selfmatch_all total (column only; the full list is not persisted):', int(RUN['n_selfmatch_all'].sum()))
# MAX_PEAKS and NMS
first = RUN[RUN['rank_key'] == 'chromatin_od']
nms = first.groupby('condition').agg(n_peaks_min=('n_peaks', 'min'), n_peaks_max=('n_peaks', 'max'), supp_mean=('n_nms_suppressed', 'mean'), n_det_min=('n_detections', 'min'))
print(nms.to_string())
ledger('diagnostics', 'MAX_PEAKS binds (n_peaks == 100) on every row; NMS suppression mean 0.918-1.72', len(RUN) + 2, int((RUN['n_peaks'] != 100).sum()) + int(round(nms['supp_mean'].min(), 3) != 0.918) + int(round(nms['supp_mean'].max(), 2) != 1.72), 'A')
ledger('diagnostics', 'n_nms_suppressed == n_peaks - n_detections', len(RUN), int((RUN['n_nms_suppressed'] != RUN['n_peaks'] - RUN['n_detections']).sum()), 'A')
edge = T30[top10].groupby(['rank_key', 'condition'])['band'].mean() * 100
print('edge-band share of top-10 slots, %:\n', edge.round(2).to_string())
EDGE_Q = {('chromatin_od', 'mo_box'): 5.41, ('chromatin_od', 'mo_pad61'): 6.46, ('tm_score', 'mo_ring5'): 2.65, ('tm_score', 'mo_pad61'): 3.06, ('chromatin_od', 'default_51'): 6.02, ('chromatin_od', 'fill_61'): 6.56, ('tm_score', 'default_51'): 3.10, ('tm_score', 'fill_61'): 3.30}
ledger('diagnostics', 'prose edge-band shares (endpoints of the quoted ranges)', len(EDGE_Q), sum(round(edge[k], 2) != v for k, v in EDGE_Q.items()), 'A')
# tp_at_10_excl_border from the top 30, where >= 10 non-band detections are in it
xb = []
for key, lst in T30.groupby(KEY):
    lst = lst.sort_values('rank')
    kept = lst[~lst['band']]
    if len(kept) < 10:
        xb.append((key, np.nan))
        continue
    fn, s = key[0], key[1]
    g = GT_BY_FILE[fn]
    g = g[g['ann_id'] != SEED_ANN[(fn, s)]].reset_index(drop=True)
    m = my_greedy(kept[['cx', 'cy']].to_numpy(float)[:10], g[['cx', 'cy']].to_numpy(float), 7.5 / MPP[fn])
    xb.append((key, int(sum(1 for q in m if q >= 0 and g['category_id'].iat[q] == 1))))
XB = pd.Series(dict(xb))
XB.index.names = KEY
xj = both_run[['tp_at_10_excl_border']].join(XB.rename('mine'))
ok = xj['mine'].notna()
ledger('diagnostics', 'tp_at_10_excl_border re-matched from the top 30 (lists with >= 10 non-band detections in it)', int(ok.sum()), int((xj.loc[ok, 'tp_at_10_excl_border'] != xj.loc[ok, 'mine']).sum()), 'A', f'({int((~ok).sum())} lists not recomputable from the top 30)')

# box statistics, tripwire, identity checks
box = TPL[TPL['condition'] == 'mo_box']
bstats = dict(bs_min=int(box['bs'].min()), bs_median=float(box['bs'].median()), bs_max=int(box['bs'].max()), bs_ge_51=int((box['bs'] >= 51).sum()), base=int((box['mode'] == 'base').sum()), fallback=int((box['mode'] == 'fallback').sum()), off_median=round(float(np.hypot(box['center_x'] - box['click_cx'], box['center_y'] - box['click_cy']).median()), 2), off_max=round(float(np.hypot(box['center_x'] - box['click_cx'], box['center_y'] - box['click_cy']).max()), 2), pad_px_all_2=bool((box['pad_px'] == 2).all()))
print(bstats)
ledger('box', 'bs 23/33/55, bs>=51 on 6, base 259 / fallback 35, offset median 2.50 max 9.66, pad_px 2', 9, sum(bstats[k] != v for k, v in dict(bs_min=23, bs_median=33.0, bs_max=55, bs_ge_51=6, base=259, fallback=35, off_median=2.5, off_max=9.66, pad_px_all_2=True).items()), 'A')
trip = box.merge(ADOPT, left_on=['file_name', 'seed_ann_id'], right_on=['file_name', 'ann_id'], suffixes=('', '_a'))
ledger('box', 'tripwire: 19 adoption clicks reproduce accepted / base_size / centre', 19 * 4 + 1, int(len(trip) != 19) + int((trip['bs'] != trip['base_size_px']).sum() + (trip['center_x'] != trip['center_x_a']).sum() + (trip['center_y'] != trip['center_y_a']).sum() + (trip['accepted'] != trip['accepted_a']).sum()), 'A')
geS = TFLAG[TFLAG['my_box_ge_S']]
nge = 0
for r in geS.itertuples():
    for axis in AXES:
        a = TOP[(TOP['file_name'] == r.file_name) & (TOP['seed_index'] == r.seed_index) & (TOP['condition'] == r.condition) & (TOP['rank_key'] == axis)].sort_values('rank')[['cx', 'cy', 'score', 'od', 'bucket']].reset_index(drop=True)
        b = TOP[(TOP['file_name'] == r.file_name) & (TOP['seed_index'] == r.seed_index) & (TOP['condition'] == 'mo_box') & (TOP['rank_key'] == axis)].sort_values('rank')[['cx', 'cy', 'score', 'od', 'bucket']].reset_index(drop=True)
        nge += int(not a.equals(b))
ledger('box', 'the 6 box_ge_S top-30 lists equal mo_box on both axes', 2 * len(geS), nge, 'A')
sha_g = TPL[TPL['condition'].isin(['mo_pad51'])].merge(box[['file_name', 'seed_index', 'tpl_sha1']], on=['file_name', 'seed_index'], suffixes=('', '_box'))
ledger('box', 'box_ge_S templates carry mo_box SHA; every other pad template does not', len(sha_g), int(((sha_g['tpl_sha1'] == sha_g['tpl_sha1_box']) != sha_g['box_ge_S']).sum()), 'A')
# smoke vs full
SMK = pd.read_csv(f'{MO}_smoke_per_run.csv', float_precision='round_trip').set_index(KEY)
full_sm = rr.loc[rr.index.intersection(SMK.index)]
sm_cols = ['tp_at_10', 'tp_at_20', 'tp_at_30', 'tp_at_10_excl_border', 'n_detections', 'tpl_sha1', 'fill_value', 'n_selfmatch_all']
nsm = sum(int((~((full_sm[c] == SMK.loc[full_sm.index, c]) | (full_sm[c].isna() & SMK.loc[full_sm.index, c].isna()))).sum()) for c in sm_cols)
ledger('determinism', 'smoke run (337, 509) rows == full-run rows', len(full_sm) * len(sm_cols), nsm, 'A', f'({len(SMK)} smoke rows)')
SMT = pd.read_csv(f'{MO}_smoke_top30.csv', float_precision='round_trip')
mt30 = TOP.merge(SMT, on=KEY + ['rank'], suffixes=('', '_s'))
ledger('determinism', 'smoke top-30 rows == full-run top-30 rows (cx, cy, score, od, bucket)', 5 * len(mt30), int(sum((mt30[c] != mt30[f'{c}_s']).sum() for c in ('cx', 'cy', 'score', 'od', 'bucket'))), 'A')

# ------------------------------------------------------------------------------------------------------------------
print('\n=== 8. what the notebook did not run: LOO-ROI, domain-level null, cross-family Holm, product metrics, LOO selection ===')
loo = []
for axis in AXES:
    for name in TESTED:
        d, n, fs_ = roi_deltas(axis, 'tp_at_10', name)
        full = my_stats(d, n, 10)
        pooled_loo, d5_loo, lo_loo = [], [], []
        for i in range(len(d)):
            keep = np.arange(len(d)) != i
            st = my_stats(d[keep], n[keep], 10)
            pooled_loo.append(st['pooled_delta_precision'])
            d5_loo.append(st['d5_bar_met'])
            lo_loo.append(st['ci_low'] if full['pooled_delta_precision'] > 0 else -st['ci_high'])
        i_max = int(np.argmax(np.abs(np.array(pooled_loo) - full['pooled_delta_precision'])))
        keep = np.arange(len(d)) != i_max
        refit = my_stats(d[keep], n[keep], 10)
        loo.append(dict(rank_key=axis, contrast=name, pooled_pts=round(100 * full['pooled_delta_precision'], 2), d5_full=full['d5_bar_met'], loo_pooled_min_pts=round(100 * min(pooled_loo), 2), loo_pooled_max_pts=round(100 * max(pooled_loo), 2), most_influential_roi=fs_[i_max], its_domain=DOMAIN[fs_[i_max]], its_delta_tp=int(d[i_max]), without_it_pooled_pts=round(100 * refit['pooled_delta_precision'], 2), without_it_ci_pts=f"[{100 * refit['ci_low']:+.2f}, {100 * refit['ci_high']:+.2f}]", without_it_wlt=f"{refit['wins']}/{refit['losses']}/{refit['ties']}", without_it_d5=refit['d5_bar_met'], loo_d5_flips=int(sum(v != full['d5_bar_met'] for v in d5_loo)), flip_rois=' '.join(fs_[i] for i, v in enumerate(d5_loo) if v != full['d5_bar_met']), loo_min_bound_toward_0_pts=round(100 * min(lo_loo), 2), n_same_direction=full['n_same_direction'], majority_needed=25, t_ci_excl0_equals_boot_excl0=full['t_ci_excludes_0'] == full['ci_excludes_0'], max_abs_t_minus_boot_bound_pts=round(100 * max(abs(full['t_ci_low'] - full['ci_low']), abs(full['t_ci_high'] - full['ci_high'])), 2)))
LOO = pd.DataFrame(loo)
print(LOO.to_string(index=False))
save(LOO, 'loo_roi')

dom = []
for axis in AXES:
    fam_rows = []
    for name in TESTED:
        d, n, fs_ = roi_deltas(axis, 'tp_at_10', name)
        dd = pd.Series(d, index=fs_).groupby(DOMAIN).sum()
        nn = pd.Series(n, index=fs_).groupby(DOMAIN).sum()
        p, nnz = exact_signflip_p(dd.to_numpy())
        x = dd / (10 * nn)
        half = sps.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
        fam_rows.append(dict(rank_key=axis, contrast=name, family=TESTED[name], n_domains=len(dd), domains_positive=int((dd > 0).sum()), domains_negative=int((dd < 0).sum()), domain_deltas_tp=' '.join(str(int(v)) for v in dd.to_numpy()), exact_p_domain=p, floor_domain=2.0 ** (1 - nnz) if nnz else 1.0, t6_ci_low_pts=100 * (x.mean() - half), t6_ci_high_pts=100 * (x.mean() + half)))
    fr = pd.DataFrame(fam_rows)
    for fam in fr['family'].unique():
        sel = fr['family'] == fam
        fr.loc[sel, 'holm_p_domain'] = holm(fr.loc[sel, 'exact_p_domain'])
    dom.append(fr)
DOM = pd.concat(dom, ignore_index=True)
print(DOM.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
save(DOM, 'domain_level')

allp = MINE[tested_mask][['rank_key', 'contrast', 'exact_p', 'holm_p']].copy()
allp['holm_p_all10'] = holm(allp['exact_p'])
print(allp.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
save(allp, 'cross_family_holm')

# product metrics: per-click deltas, worst click, worst seed per ROI, zero lists
prod = []
for axis in AXES:
    for name in list(TESTED) + ['mo_box - default_51', 'mo_pad51 - default_51', 'mo_pad61 - default_51', 'fill_51 - default_51', 'fill_61 - default_51']:
        a, b = contrast_arms(name)
        va = CLICK.loc[(axis, a), 'tp_at_10']
        vb = CLICK.loc[(axis, b), 'tp_at_10']
        dc = (va - vb.reindex(va.index))
        wa = va.groupby(level='file_name').min()
        wb = vb.groupby(level='file_name').min()
        dw = (wa - wb.reindex(wa.index)).reindex(ROI_ORDER)
        pw, _ = exact_signflip_p(dw.to_numpy())
        prod.append(dict(rank_key=axis, contrast=name, clicks_better=int((dc > 0).sum()), clicks_worse=int((dc < 0).sum()), clicks_tied=int((dc == 0).sum()), worst_click_delta_tp=int(dc.min()), best_click_delta_tp=int(dc.max()), clicks_losing_3plus=int((dc <= -3).sum()), clicks_gaining_3plus=int((dc >= 3).sum()), worst_seed_precision_a=round(wa.mean() / 10, 4), worst_seed_precision_b=round(wb.mean() / 10, 4), worst_seed_delta_pts=round(100 * dw.mean() / 10, 2), worst_seed_W_L_T=f'{int((dw > 0).sum())}/{int((dw < 0).sum())}/{int((dw == 0).sum())}', worst_seed_exact_p=pw, zero_lists_a=int((va == 0).sum()), zero_lists_b=int((vb == 0).sum())))
PROD = pd.DataFrame(prod)
print(PROD.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
save(PROD, 'product_metrics')

# LOO selection of the pad (mode 10): pick the best pad on 48 ROIs, score it on the held-out one
sel = []
for axis in AXES:
    for pool_name, pool in (('3 pads', ('mo_pad51', 'mo_pad61', 'mo_ring5')), ('all 7 arms', ALL_ARMS)):
        tot = {arm: CLICK.loc[(axis, arm), 'tp_at_10'].groupby(level='file_name').sum() for arm in pool}
        insample_best = max(pool, key=lambda arm: tot[arm].sum())
        held, picks = 0, []
        for fn in ROI_ORDER:
            best = max(pool, key=lambda arm: (tot[arm].sum() - tot[arm][fn], -pool.index(arm)))
            picks.append(best)
            held += tot[best][fn]
        sel.append(dict(rank_key=axis, pool=pool_name, insample_best=insample_best, insample_precision=round(tot[insample_best].sum() / 2940, 4), loo_selected_precision=round(held / 2940, 4), picks=pd.Series(picks).value_counts().to_dict()))
SEL = pd.DataFrame(sel)
print(SEL.to_string(index=False))
save(SEL, 'loo_selection')

# context: length-matched random-list precision null at K = 10 (mode 4), per click from the annotation database
nullr = []
for fn, cl in CLICKS.groupby('file_name'):
    g = GT_BY_FILE[fn]
    g = g[g['category_id'] == 1].reset_index(drop=True)
    W, H = WH.loc[fn, 'width'], WH.loc[fn, 'height']
    ys, xs = np.mgrid[0:H:8, 0:W:8]
    dist, nn = cKDTree(g[['cx', 'cy']].to_numpy(float)).query(np.c_[xs.ravel(), ys.ravel()], k=2)
    for r in cl.itertuples():
        own = int(np.flatnonzero(g['ann_id'].to_numpy() == r.seed_ann_id)[0])
        d = np.where(nn[:, 0] == own, dist[:, 1], dist[:, 0])
        hit = float((d <= 7.5 / r.mpp).mean())
        n_gt = len(g) - 1
        nullr.append(dict(file_name=fn, seed_index=r.seed_index, n_gt=n_gt, hit_frac=hit, null_p10=(n_gt * (1 - (1 - hit / n_gt) ** 10) / 10) if n_gt else 0.0))
NULL = pd.DataFrame(nullr)
print(f"length-matched null precision@10, pooled over clicks: {NULL['null_p10'].mean():.4f}; clicks with no other mitosis in the ROI: {int((NULL['n_gt'] == 0).sum())}")
save(NULL, 'precision_null')

# ------------------------------------------------------------------------------------------------------------------
print('\n=== 9. figures: the plotted series and the gallery captions ===')
order = list(TESTED) + DESCRIPTIVE
fig1 = MINE[(MINE['metric'] == 'tp_at_10') & (MINE['version'] == 'recorded')].set_index(['rank_key', 'contrast']).loc[[(a, c) for a in AXES for c in order], ['pooled_delta_precision', 'ci_low', 'ci_high', 'wins', 'losses', 'ties', 'holm_p']].reset_index()
fig1['annotation'] = [f"{100 * r.pooled_delta_precision:+.2f}  W/L/T {r.wins}/{r.losses}/{r.ties}" + (f'  Holm {r.holm_p:.2g}' if r.contrast in TESTED else '') for r in fig1.itertuples()]
save(fig1, 'figure1_series')
print(f"figure 1 data span: {100 * fig1['ci_low'].min():+.2f} to {100 * fig1['ci_high'].max():+.2f}")
ledger('figures', 'figure-1 printed span -7.35 to +10.65', 2, int(round(100 * fig1['ci_low'].min(), 2) != -7.35) + int(round(100 * fig1['ci_high'].max(), 2) != 10.65), 'A')
# gallery picks, following cell 27's rule on the templates CSV
TA = TPL.copy()
clk = TA[TA['condition'] == 'mo_box'].reset_index(drop=True)
picks = []


def pick(i, why):
    """
        Add one click to the gallery unless present, up to 8.

        i (int): row of clk.
        why (str): the reason.

        Returns None.
    """
    k = (clk.loc[i, 'file_name'], int(clk.loc[i, 'seed_index']))
    if k not in [p[0] for p in picks] and len(picks) < 8:
        picks.append((k, why))


pick(int(clk['bs'].idxmin()), 'smallest bs')
pick(int(clk['bs'].idxmax()), 'largest bs')
pick(int(clk['offset_px'].idxmax()), 'largest centre offset')
for flag in ('box_outside_blank', 'template_outside_window', 'roi_edge_read', 'box_ge_S', 'mass_weight_nonpositive'):
    hit = TA[TA[flag]]
    fresh = [k for k in zip(hit['file_name'], hit['seed_index']) if k not in [p[0] for p in picks]]
    if fresh:
        pick(int(clk.index[(clk['file_name'] == fresh[0][0]) & (clk['seed_index'] == fresh[0][1])][0]), flag)
for q in (0.25, 0.5, 0.75, 0.1, 0.9):
    target = clk['bs'].quantile(q)
    for i in (clk['bs'] - target).abs().sort_values(kind='stable').index:
        before = len(picks)
        pick(int(i), f'bs near q{int(100 * q)}')
        if len(picks) > before or len(picks) == 8:
            break
gal = []
for (fn, s), why in picks:
    for arm in ARMS:
        r = TA[(TA['file_name'] == fn) & (TA['seed_index'] == s) & (TA['condition'] == arm)].iloc[0]
        gal.append(dict(file_name=fn, seed_index=s, why=why, arm=arm, S=int(r['S']), bs=int(r['bs']), offset_cheb=int(r['offset_cheb_px']), fill=None if np.isnan(r['fill_value']) else round(r['fill_value'], 3), mass_wt=None if np.isnan(r['mass_weight']) else round(r['mass_weight'], 3), flags='; '.join(f for f in ('box_outside_blank', 'template_outside_window', 'roi_edge_read', 'box_ge_S', 'mass_weight_nonpositive') if r[f])))
GAL = pd.DataFrame(gal)
print(GAL.to_string(index=False))
save(GAL, 'figure2_captions')
PRINTED_PICKS = ['337.tiff s0', '425.tiff s1', '128.tiff s4', '007.tiff s3', '300.tiff s0', '128.tiff s0', '301.tiff s2', '322.tiff s0']
ledger('figures', 'figure-2 gallery picks (8 clicks, in order)', 8, sum(a != f'{fn} s{s}' for a, ((fn, s), _) in zip(PRINTED_PICKS, picks)), 'A')

# ------------------------------------------------------------------------------------------------------------------
print('\n=== 10. Tier B: back to the pixels ===')


def load_rgb(path):
    """
        Full-resolution RGB for one ROI TIFF, alpha dropped.

        path (str): TIFF path.

        Returns np.ndarray: H x W x 3 uint8.
    """
    with tifffile.TiffFile(path) as tf:
        s = tf.series[0]
        a = s.levels[0].asarray() if getattr(s, 'levels', None) else s.asarray()
    return np.ascontiguousarray(a[:, :, :3])


def tiff_mpp(path):
    """
        Microns per pixel from the TIFF resolution tags.

        path (str): TIFF path.

        Returns float: um / px.
    """
    with tifffile.TiffFile(path) as tf:
        p = tf.pages[0]
        num, den = p.tags['XResolution'].value
        unit = int(p.tags['ResolutionUnit'].value)
    return {2: 25400.0, 3: 10000.0}[unit] / (num / den)


def hem_od(rgb):
    """
        Hematoxylin optical density, the D3 channel.

        rgb (np.ndarray): uint8 RGB.

        Returns np.ndarray: float32 H x W.
    """
    return rgb2hed(rgb.astype(np.float32) / 255.0)[:, :, 0].astype(np.float32)


def window(img, cx, cy, side):
    """
        A side x side window centred on the rounded point, None when it leaves the image.

        img (np.ndarray): the channel.
        cx (float): x.
        cy (float): y.
        side (int): odd side.

        Returns np.ndarray or None: a float32 copy.
    """
    h = side // 2
    x, y = int(round(cx)), int(round(cy))
    if x - h < 0 or y - h < 0 or x + h >= img.shape[1] or y + h >= img.shape[0]:
        return None
    return np.array(img[y - h:y + h + 1, x - h:x + h + 1], dtype=np.float32, copy=True)


def my_box(hem, cx, cy, mpp):
    """
        The multi-Otsu box, re-implemented from the method's description with scipy.ndimage for the morphology.

        hem (np.ndarray): hematoxylin_od channel.
        cx (float): click x.
        cy (float): click y.
        mpp (float): um / px.

        Returns dict: accepted, mode, bs, center_x, center_y.
    """
    p = window(hem, cx, cy, 51).astype(np.float64)
    thr = threshold_multiotsu(p, classes=3)
    top = ndi.binary_dilation(p > thr[-1], structure=ndi.generate_binary_structure(2, 1), iterations=1)
    lab, nlab = ndi.label(top)
    pieces = {}
    for L in range(1, nlab + 1):
        ys, xs = np.nonzero(lab == L)
        pieces[L] = dict(area=len(ys) * mpp ** 2, cen=np.hypot(ys.mean() - 25, xs.mean() - 25) * mpp, ext=max(ys.max() - ys.min() + 1, xs.max() - xs.min() + 1) * mpp)
    c = int(lab[25, 25])
    if c and pieces[c]['area'] >= 1.0 and pieces[c]['ext'] >= 5.0:
        mode, chosen = 'base', [c]
    else:
        mode, chosen = 'fallback', [L for L, q in pieces.items() if q['area'] >= 1.0 and q['cen'] <= 3.75]
    if not chosen:
        return dict(accepted=False, mode=mode)
    mask = np.isin(lab, chosen)
    if mask.sum() * mpp ** 2 < 3.0 or mask.sum() > 0.85 * mask.size:
        return dict(accepted=False, mode=mode)
    ys, xs = np.nonzero(mask)
    y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
    side = max(y1 - y0, x1 - x0)
    side = side + 1 if side % 2 == 0 else side
    side = max(5, side)
    pad = int(round(0.5 / mpp))
    return dict(accepted=True, mode=mode, bs=int(side + 2 * pad), center_x=int(round(cx)) - 25 + (x0 + x1 - 1) / 2.0, center_y=int(round(cy)) - 25 + (y0 + y1 - 1) / 2.0)


def keep_fill(raw, keep):
    """
        Keep a region and set every other pixel to the float32 mean (in float64) of those other pixels.

        raw (np.ndarray): float32 window.
        keep (np.ndarray): boolean keep mask.

        Returns tuple[np.ndarray, float]: the template and the fill (nan when nothing is filled).
    """
    T = np.array(raw, dtype=np.float32, copy=True)
    if keep.all():
        return T, float('nan')
    f = np.float32(raw[~keep].astype(np.float64).mean())
    T[~keep] = f
    return T, float(f)


def sha1(T):
    """
        SHA-1 of a template's float32 bytes.

        T (np.ndarray): template.

        Returns str: hex digest.
    """
    return hashlib.sha1(np.ascontiguousarray(T, dtype=np.float32).tobytes()).hexdigest()


def my_templates(hem, tcx, tcy, bs):
    """
        The four arms' templates at the box centre.

        hem (np.ndarray): channel.
        tcx (int): box centre x.
        tcy (int): box centre y.
        bs (int): box side.

        Returns dict: arm -> (T, fill, m_keep).
    """
    out = {}
    for arm in ARMS:
        S = {'mo_box': bs, 'mo_pad51': max(51, bs), 'mo_pad61': max(61, bs), 'mo_ring5': bs + 10}[arm]
        raw = window(hem, tcx, tcy, S)
        o = (S - bs) // 2
        keep = np.zeros((S, S), bool)
        keep[o:o + bs, o:o + bs] = True
        T, f = keep_fill(raw, keep)
        out[arm] = (T, f, float(raw[keep].astype(np.float64).mean()))
    return out


def my_ref_templates(hem, cx, cy, rect):
    """
        The three references, click-centred, from the recorded 51-frame rectangle.

        hem (np.ndarray): channel.
        cx (float): click x.
        cy (float): click y.
        rect (tuple): (y0, y1, x0, x1) half-open in the 51 frame.

        Returns dict: arm -> T.
    """
    out = {}
    w51 = window(hem, cx, cy, 51)
    out['default_51'] = w51
    y0, y1, x0, x1 = rect
    k = np.zeros((51, 51), bool)
    k[y0:y1, x0:x1] = True
    out['fill_51'] = keep_fill(w51, k)[0]
    w61 = window(hem, cx, cy, 61)
    k = np.zeros((61, 61), bool)
    k[y0 + 5:y1 + 5, x0 + 5:x1 + 5] = True
    out['fill_61'] = keep_fill(w61, k)[0]
    return out


def my_search(hem, T, blank_xy, blank_side, mr, od_pad_img, near_xy=None):
    """
        One search, re-implemented: blank a square to the channel minimum, replicate-pad, TM_CCOEFF, robust floor, 15 px local maxima, top 100, greedy distance NMS, od.

        hem (np.ndarray): the unblanked channel.
        T (np.ndarray): the S x S template.
        blank_xy (tuple): the blank's centre (x, y), rounded here.
        blank_side (int): the blank's odd side.
        mr (float): NMS radius = match radius, px.
        od_pad_img (np.ndarray): the channel replicate-padded by 25 for od.
        near_xy (tuple or None): a point (the click) around which to report the largest response within 80 px (Chebyshev).

        Returns pd.DataFrame or tuple[pd.DataFrame, dict]: detections in score order with cx, cy, score, od; with near_xy, also the near-click maximum and the 100th peak's score.
    """
    H, W = hem.shape
    ch = hem.copy()
    hb = blank_side // 2
    bx, by = int(round(blank_xy[0])), int(round(blank_xy[1]))
    ch[max(0, by - hb):min(H - 1, by + hb) + 1, max(0, bx - hb):min(W - 1, bx + hb) + 1] = float(hem.min())
    pad = (T.shape[0] - 1) // 2
    R = cv2.matchTemplate(cv2.copyMakeBorder(ch, pad, pad, pad, pad, cv2.BORDER_REPLICATE), T, cv2.TM_CCOEFF)
    assert R.shape == (H, W)
    samp = R[::8, ::8].ravel()
    med = float(np.median(samp))
    thr = med - 1.5 * 1.4826 * float(np.median(np.abs(samp - med)))
    loc = (R >= ndi.maximum_filter(R, size=15, mode='nearest')) & (R >= thr)
    ys, xs = np.nonzero(loc)
    sc = R[ys, xs]
    o = np.lexsort((ys, xs, -sc))[:100]
    xs, ys, sc = xs[o].astype(float), ys[o].astype(float), sc[o]
    keep, dead = [], np.zeros(len(xs), bool)
    for i in np.argsort(-sc.astype(np.float64), kind='stable'):
        if dead[i]:
            continue
        keep.append(i)
        dead |= np.hypot(xs - xs[i], ys - ys[i]) <= mr
    keep = np.array(keep)
    det = pd.DataFrame(dict(cx=xs[keep], cy=ys[keep], score=sc[keep]))
    od = []
    for x, y in zip(det['cx'], det['cy']):
        p = od_pad_img[int(y):int(y) + 51, int(x):int(x) + 51].ravel()
        od.append(float(np.partition(p, -260)[-260:].mean()))
    det = det.assign(od=od)
    if near_xy is None:
        return det
    nx, ny = int(round(near_xy[0])), int(round(near_xy[1]))
    return det, dict(near_click_max_response=float(R[max(0, ny - 80):ny + 81, max(0, nx - 80):nx + 81].max()), peak_cut_score=float(sc.min()))


def rank_match(det, axis, gxy, gcat, gid, mr):
    """
        Rank a detection list on one axis and match it.

        det (pd.DataFrame): unranked detections.
        axis (str): chromatin_od or tm_score.
        gxy (np.ndarray): GT centres.
        gcat (np.ndarray): GT categories.
        gid (np.ndarray): GT ids.
        mr (float): match radius.

        Returns pd.DataFrame: ranked detections with hit and ann columns.
    """
    r = det.sort_values('od' if axis == 'chromatin_od' else 'score', ascending=False, kind='mergesort', na_position='last').reset_index(drop=True)
    m = my_greedy(r[['cx', 'cy']].to_numpy(float), gxy, mr)
    return r.assign(rank=np.arange(len(r)), hit=[q >= 0 and gcat[q] == 1 for q in m], ann=[int(gid[q]) if q >= 0 else -1 for q in m])


STRESS = sorted(set(zip(TFLAG.loc[TFLAG['my_box_outside_blank'] | TFLAG['my_mass_weight_nonpositive'], 'file_name'], TFLAG.loc[TFLAG['my_box_outside_blank'] | TFLAG['my_mass_weight_nonpositive'], 'seed_index'])))
PB_ROIS = [next(f for f in ROI_ORDER if DOMAIN[f] == d) for d in sorted(DOMAIN.unique())]
PB_ARMS = ('mo_box', 'mo_pad51', 'mo_pad61', 'mo_ring5', 'fill_61')
print(f'stress clicks ({len(STRESS)}): {STRESS}')
print(f'production-blank ROIs (first ROI of each domain in ROI order): {PB_ROIS}')
TPL_IX = TPL.set_index(['file_name', 'seed_index', 'condition'])
box_rows, tpl_rows, search_rows, pb_rows, top_div_rows = [], [], [], [], []
tb0 = time.time()
n_searches = 0
for fn in ROI_ORDER:
    cl = CLICKS[CLICKS['file_name'] == fn].sort_values('seed_index')
    path = f"{SUBSET_DIRS[cl['subset'].iat[0]]}/{fn}"
    mpp = tiff_mpp(path)
    hem = hem_od(load_rgb(path))
    need_search = fn in PB_ROIS or any(f == fn for f, _ in STRESS)
    od_img = cv2.copyMakeBorder(hem, 25, 25, 25, 25, cv2.BORDER_REPLICATE) if need_search else None
    g = GT_BY_FILE[fn]
    mr = 7.5 / mpp
    for r in cl.itertuples():
        s = int(r.seed_index)
        b = my_box(hem, r.click_cx, r.click_cy, mpp)
        rec = TPL_IX.loc[(fn, s, 'mo_box')]
        box_rows.append(dict(file_name=fn, seed_index=s, mpp_tiff=mpp, mpp_recorded=r.mpp, accepted=b['accepted'], mode=b['mode'], mode_rec=rec['mode'], bs=b.get('bs'), bs_rec=int(rec['bs']), center_x=b.get('center_x'), center_x_rec=rec['center_x'], center_y=b.get('center_y'), center_y_rec=rec['center_y']))
        tcx, tcy, bs = int(rec['tcx']), int(rec['tcy']), int(rec['bs'])
        mt_ = my_templates(hem, tcx, tcy, bs)
        for arm, (T, f, mk) in mt_.items():
            q = TPL_IX.loc[(fn, s, arm)]
            tpl_rows.append(dict(file_name=fn, seed_index=s, arm=arm, sha_match=sha1(T) == q['tpl_sha1'], fill_match=(np.isnan(f) and np.isnan(q['fill_value'])) or np.float32(f) == np.float32(q['fill_value']), mass_weight_mine=(mk - f) if not np.isnan(f) else np.nan, mass_weight_rec=q['mass_weight'], m_keep_match=np.isclose(mk, q['m_keep'], rtol=1e-12)))
        rt = my_ref_templates(hem, r.click_cx, r.click_cy, (int(r.rect_y0), int(r.rect_y1), int(r.rect_x0), int(r.rect_x1)))
        for arm, T in rt.items():
            q = REF_RUN[(REF_RUN['file_name'] == fn) & (REF_RUN['seed_index'] == s) & (REF_RUN['condition'] == arm)]['tpl_sha1'].iat[0]
            tpl_rows.append(dict(file_name=fn, seed_index=s, arm=arm, sha_match=sha1(T) == q))
        if not need_search:
            continue
        ge = g[g['ann_id'] != r.seed_ann_id].reset_index(drop=True)
        gxy, gcat, gid = ge[['cx', 'cy']].to_numpy(float), ge['category_id'].to_numpy(), ge['ann_id'].to_numpy()
        if (fn, s) in STRESS:
            for arm in ARMS + ('fill_51', 'fill_61'):
                det = my_search(hem, mt_[arm][0] if arm in ARMS else rt[arm], (r.click_cx, r.click_cy), 51, mr, od_img)
                n_searches += 1
                tm_rk = rank_match(det, 'tm_score', gxy, gcat, gid, mr)
                cxr, cyr = (tcx, tcy) if arm in ARMS else (r.click_cx, r.click_cy)
                selfm = [int(((~tm_rk['hit']) & ((np.maximum(np.abs(tm_rk['cx'] - r.click_cx), np.abs(tm_rk['cy'] - r.click_cy)) <= SELF_PX) | (np.maximum(np.abs(tm_rk['cx'] - cxr), np.abs(tm_rk['cy'] - cyr)) <= SELF_PX))).sum())]
                full_near = int((np.maximum(np.abs(det['cx'] - int(round(r.click_cx))), np.abs(det['cy'] - int(round(r.click_cy)))) <= HALF_BLANK).sum())
                for axis in AXES:
                    rk = rank_match(det, axis, gxy, gcat, gid, mr)
                    src_top = TOP if arm in ARMS else REF_TOP
                    sv = src_top[(src_top['file_name'] == fn) & (src_top['seed_index'] == s) & (src_top['condition'] == arm) & (src_top['rank_key'] == axis)].sort_values('rank').reset_index(drop=True)
                    h30 = rk.head(30)
                    div = dict(cx=int((h30['cx'].to_numpy() != sv['cx'].to_numpy()).sum()), cy=int((h30['cy'].to_numpy() != sv['cy'].to_numpy()).sum()), score=int((h30['score'].to_numpy().astype(np.float32) != sv['score'].to_numpy().astype(np.float32)).sum()), od=int((~np.isclose(h30['od'].to_numpy(), sv['od'].to_numpy(), rtol=1e-9, atol=0)).sum()), ann=int((h30['ann'].to_numpy() != sv['matched_ann_id'].to_numpy()).sum()))
                    src_run = RUN if arm in ARMS else REF_RUN
                    saved = src_run[(src_run['file_name'] == fn) & (src_run['seed_index'] == s) & (src_run['condition'] == arm) & (src_run['rank_key'] == axis)].iloc[0]
                    search_rows.append(dict(file_name=fn, seed_index=s, arm=arm, rank_key=axis, n_detections=len(det), n_detections_rec=int(saved['n_detections']), tp10=int(rk['hit'].head(10).sum()), tp10_rec=int(saved['tp_at_10']), tp30=int(rk['hit'].head(30).sum()), tp30_rec=int(saved['tp_at_30']), selfmatch_full_list=int(selfm[0]), dets_inside_blank=full_near, **{f'top30_div_{k}': v for k, v in div.items()}))
        if fn in PB_ROIS:
            for arm in PB_ARMS:
                if arm == 'fill_61':
                    T, bxy, side = rt['fill_61'], (r.click_cx, r.click_cy), 61
                else:
                    T = mt_[arm][0]
                    bxy, side = (tcx, tcy), T.shape[0]
                det, dg = my_search(hem, T, bxy, side, mr, od_img, near_xy=(r.click_cx, r.click_cy))
                n_searches += 1
                for axis in AXES:
                    rk = rank_match(det, axis, gxy, gcat, gid, mr)
                    base_ = (RUN if arm in ARMS else REF_RUN)
                    saved = base_[(base_['file_name'] == fn) & (base_['seed_index'] == s) & (base_['condition'] == arm) & (base_['rank_key'] == axis)].iloc[0]
                    near = rk.head(10)
                    nfp = int(((~near['hit']) & (np.maximum(np.abs(near['cx'] - r.click_cx), np.abs(near['cy'] - r.click_cy)) <= NEAR_PX)).sum())
                    src_top = TOP if arm in ARMS else REF_TOP
                    sv = src_top[(src_top['file_name'] == fn) & (src_top['seed_index'] == s) & (src_top['condition'] == arm) & (src_top['rank_key'] == axis)].sort_values('rank')
                    same30 = bool(np.array_equal(rk.head(30)[['cx', 'cy']].to_numpy(), sv[['cx', 'cy']].to_numpy()))
                    pb_rows.append(dict(file_name=fn, domain=DOMAIN[fn], seed_index=s, arm=arm, rank_key=axis, blank_side=side, blank_at='box centre' if arm in ARMS else 'click', tp10_prodblank=int(rk['hit'].head(10).sum()), tp10_notebook=int(saved['tp_at_10']), n_det_prodblank=len(det), n_det_notebook=int(saved['n_detections']), top30_identical_to_notebook=same30, near_seed_fp_top10_60=nfp, **dg))
    print(f'  {fn}: {time.time() - tb0:.0f} s, {n_searches} searches', flush=True)
BOXB = pd.DataFrame(box_rows)
TPLB = pd.DataFrame(tpl_rows)
SRB = pd.DataFrame(search_rows)
PBB = pd.DataFrame(pb_rows)
save(BOXB, 'tierb_boxes')
save(TPLB, 'tierb_templates')
save(SRB, 'tierb_searches')
save(PBB, 'tierb_prodblank')
nbx = int((BOXB['mode'] != BOXB['mode_rec']).sum() + (BOXB['bs'] != BOXB['bs_rec']).sum() + (BOXB['center_x'] != BOXB['center_x_rec']).sum() + (BOXB['center_y'] != BOXB['center_y_rec']).sum() + (~BOXB['accepted']).sum())
ledger('tier B', 'multi-Otsu box re-implemented (scipy morphology) on all 294 clicks: accepted, mode, bs, centre', 5 * len(BOXB), nbx, 'B')
ledger('tier B', 'mpp from the TIFF tags == recorded mpp', len(BOXB), int((BOXB['mpp_tiff'] != BOXB['mpp_recorded']).sum()), 'B')
ledger('tier B', 'template SHA-1 rebuilt from pixels: 4 arms x 294 + 3 references x 294', len(TPLB), int((~TPLB['sha_match']).sum()), 'B')
tf_ = TPLB[TPLB['arm'].isin(ARMS)]
ledger('tier B', 'fill value and m_keep rebuilt from pixels (1176 x 2)', 2 * len(tf_), int((~tf_['fill_match'].astype(bool)).sum() + (~tf_['m_keep_match'].astype(bool)).sum()), 'B')
mwn = tf_[tf_['mass_weight_mine'].notna()]
ledger('tier B', 'mass weight sign rebuilt from pixels', len(mwn), int(((mwn['mass_weight_mine'] <= 0) != (mwn['mass_weight_rec'] <= 0)).sum()), 'B', f"(min {mwn.groupby('arm')['mass_weight_mine'].min().round(5).to_dict()})")
cols = [c for c in SRB.columns if c.startswith('top30_div_')]
ledger('tier B', 'own search on the 17 stress clicks x (4 arms + fill_51 + fill_61) x 2 axes: top-30 cx, cy, float32 score, od, matched ann', len(SRB) * 30 * len(cols), int(SRB[cols].to_numpy().sum()), 'B')
ledger('tier B', 'own search: n_detections, TP@10, TP@30 on the stress clicks', 3 * len(SRB), int((SRB['n_detections'] != SRB['n_detections_rec']).sum() + (SRB['tp10'] != SRB['tp10_rec']).sum() + (SRB['tp30'] != SRB['tp30_rec']).sum()), 'B')
print(f"full-list self-matches on the stress clicks (tm_score list; the list is the same set on both axes): {int(SRB[SRB['rank_key'] == 'tm_score']['selfmatch_full_list'].sum())}; detections inside the 51 px blank: {int(SRB['dets_inside_blank'].sum())}")

# production blank: effect on TP@10 and on the tested contrasts, on the 7-ROI subset
pbc = PBB.assign(changed=PBB['tp10_prodblank'] != PBB['tp10_notebook'], d=PBB['tp10_prodblank'] - PBB['tp10_notebook'])
pb_arm = pbc.groupby(['rank_key', 'arm']).agg(clicks=('changed', 'size'), clicks_changed=('changed', 'sum'), tp_gained=('d', lambda v: int(v[v > 0].sum())), tp_lost=('d', lambda v: int(-v[v < 0].sum())), net_pts=('d', lambda v: round(100 * v.sum() / (10 * len(v)), 2)), near_fp60=('near_seed_fp_top10_60', 'sum')).reset_index()
print(pb_arm.to_string(index=False))
pb_mech = pbc[pbc['rank_key'] == 'tm_score']
print(f"production-blank searches: {len(pb_mech)}; top-30 identical to the notebook's on {int(pbc['top30_identical_to_notebook'].sum())} of {len(pbc)} (search, axis) lists; n_detections equal on {int((pb_mech['n_det_prodblank'] == pb_mech['n_det_notebook']).sum())}; max response within 80 px of the click below the 100th peak's score on {int((pb_mech['near_click_max_response'] < pb_mech['peak_cut_score']).sum())} of {len(pb_mech)} searches (median ratio {float((pb_mech['near_click_max_response'] / pb_mech['peak_cut_score']).median()):.2f}, max {float((pb_mech['near_click_max_response'] / pb_mech['peak_cut_score']).max()):.2f})")
ledger('tier B', 'production (D11) blank on 7 ROIs x 6 clicks x 5 arms x 2 axes: TP@10 unchanged from the notebook blank', len(pbc), int(pbc['changed'].sum()), 'B', '(scope test: the notebook blank is not D11)')
save(pb_arm, 'tierb_prodblank_by_arm')
sub = ALLRUN[ALLRUN['file_name'].isin(PB_ROIS)]
pbw = pbc.pivot_table(index=['rank_key', 'file_name', 'seed_index'], columns='arm', values='tp10_prodblank')
nbw = sub.pivot_table(index=['rank_key', 'file_name', 'seed_index'], columns='condition', values='tp_at_10')
pbw['fill_51'] = nbw['fill_51']
pbw['default_51'] = nbw['default_51']
pbc_rows = []
for axis in AXES:
    for name in TESTED:
        a, b = contrast_arms(name)
        for label_, tab in (('notebook blank (51 at click)', nbw), ('production blank (S at template centre)', pbw)):
            t_ = tab.loc[axis]
            dd = (t_[a] - t_[b]).groupby(level='file_name').sum().reindex(PB_ROIS)
            p, _ = exact_signflip_p(dd.to_numpy())
            pbc_rows.append(dict(rank_key=axis, contrast=name, blank=label_, n_roi=len(dd), pooled_pts=round(100 * dd.sum() / (10 * 6 * len(dd)), 2), W_L_T=f'{int((dd > 0).sum())}/{int((dd < 0).sum())}/{int((dd == 0).sum())}', exact_p_7roi=p))
PBC = pd.DataFrame(pbc_rows)
print(PBC.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
save(PBC, 'tierb_prodblank_contrasts')
print(f'Tier B pixel work: {(time.time() - tb0) / 60:.1f} min, {n_searches} searches')

# re-timing: mo_box vs mo_pad61, 3 repeats each, on one click
rt_rows = []
r = CLICKS[(CLICKS['file_name'] == '300.tiff') & (CLICKS['seed_index'] == 0)].iloc[0]
path = f"{SUBSET_DIRS[r.subset]}/300.tiff"
mpp = tiff_mpp(path)
hem = hem_od(load_rgb(path))
od_img = cv2.copyMakeBorder(hem, 25, 25, 25, 25, cv2.BORDER_REPLICATE)
rec = TPL_IX.loc[('300.tiff', 0, 'mo_box')]
mt_ = my_templates(hem, int(rec['tcx']), int(rec['tcy']), int(rec['bs']))
load_avg = os.getloadavg()[0]
for rep in range(3):
    for arm in ('mo_box', 'mo_pad61'):
        t0 = time.time()
        my_search(hem, mt_[arm][0], (r.click_cx, r.click_cy), 51, 7.5 / mpp, od_img)
        rt_rows.append(dict(arm=arm, rep=rep, S=mt_[arm][0].shape[0], seconds=round(time.time() - t0, 3), load_avg_1min=round(load_avg, 2)))
RT = pd.DataFrame(rt_rows)
print(RT.groupby('arm')['seconds'].describe().to_string())
save(RT, 'retiming')
nb_med = RUN[RUN['rank_key'] == 'chromatin_od'].groupby('condition')['t_search_s'].agg(['median', 'min', 'max', 'std'])
print('notebook t_search_s per arm (one run, across 294 clicks):\n', nb_med.to_string())
save(nb_med.reset_index(), 'notebook_timing')
ledger('timing', 'per-search medians 2.48 / 2.90 / 3.29 / 2.63 s from the per_run column', 4, sum(round(nb_med.loc[a, 'median'], 2) != v for a, v in (('mo_box', 2.48), ('mo_pad51', 2.90), ('mo_pad61', 3.29), ('mo_ring5', 2.63))), 'A')

# ------------------------------------------------------------------------------------------------------------------
print('\n=== 11. ties at the K = 10 boundary, and the agent-file appendix facts, recounted ===')
tt_ = TOP.assign(key=np.where(TOP['rank_key'] == 'tm_score', TOP['score'], TOP['od'])).sort_values(KEY + ['rank'])
b9 = tt_[tt_['rank'] == 9].set_index(KEY)['key']
b10 = tt_[tt_['rank'] == 10].set_index(KEY)['key']
ties = dict(boundary_ties_rank9_rank10=int((b9 == b10).sum()), lists=len(b9), duplicate_keys_within_top30=int(tt_.groupby(KEY)['key'].apply(lambda v: v.duplicated().sum()).sum()))
print(ties)
save(pd.DataFrame([ties]), 'ties')
import glob
import re
nbs = glob.glob('**/*.ipynb', recursive=True)
fs_src = open('midog_utils/find_and_suppress.py').read()
ds_src = open('midog_utils/dataset.py').read()
ev_src = open('midog_utils/evaluate.py').read()
csvs = glob.glob('results/*.csv')
untracked = [l for l in git('status', '--porcelain', '--untracked-files=all', '--', 'results/').splitlines() if l.startswith('??') and l.endswith('.csv') and not l[3:].startswith(OUT)]
csvs = [c for c in csvs if not c.startswith(OUT)]
app = pd.DataFrame([dict(fact='notebooks found by the recount glob', appendix='59', now=str(len(nbs))), dict(fact='embedded PNGs across notebooks', appendix='217', now=str(sum(1 for q in nbs for c in json.load(open(q))['cells'] if c['cell_type'] == 'code' for o in c.get('outputs', []) if 'image/png' in o.get('data', {})))), dict(fact='results/*.csv total / untracked', appendix='233 / 0 (2026-09-12)', now=f'{len(csvs)} / {len(untracked)}'), dict(fact='FSConfig.tm_method default', appendix='TM_CCOEFF_NORMED', now=re.search(r'tm_method: int = (\S+)', fs_src).group(1)), dict(fact='FSConfig.nms_radius default', appendix="None = this image's match radius", now=re.search(r'nms_radius: float = None\s*#\s*(.*)', fs_src).group(1).strip()), dict(fact='image_annotations category_id default', appendix='None', now=re.search(r'def image_annotations\(annotations: pd.DataFrame, file_name: str, category_id=(\w+)\)', ds_src).group(1)), dict(fact='MIDOG_RADIUS_UM', appendix='7.5', now=re.search(r'MIDOG_RADIUS_UM = ([0-9.]+)', ev_src).group(1))])
print(app.to_string(index=False))
save(app, 'appendix_drift')

# ------------------------------------------------------------------------------------------------------------------
LED = pd.DataFrame(LEDGER)
save(LED, 'ledger')
print('\n=== ledger ===')
print(LED.to_string(index=False))
print(f"\nTOTAL: {LED['n_compared'].sum()} values compared, {LED['n_diverged'].sum()} diverged; by tier {LED.groupby('tier')[['n_compared', 'n_diverged']].sum().to_dict('index')}")
print(f'audit runtime {(time.time() - T0) / 60:.1f} min')
