"""
    Independent audit of production_hematoxylin_only/multiotsu_pad_arms/tm_border_floor_49roi_6seed.ipynb.

    Everything is re-derived from the notebook's persisted artifacts (results/tm_border_floor_49roi_6seed_*.csv),
    from databases/MIDOG++.json, from the ROI TIFFs, and -- for the pool-level claims the notebook never persisted --
    from its search pickle in the system temp directory. Nothing in the notebook's helpers is imported: matching,
    sign-flip, Holm, the search and the synthetic demonstrations are re-implemented here. The only repo code used is
    midog_utils.dataset.load_roi / roi_mpp (TIFF reading) and midog_utils.channels.to_channel (colour deconvolution).

    Run: /Users/mohinianand/anaconda3/bin/python3 tm_border_floor_49roi_6seed_audit.py
    Env: AUDIT_TIERB_ROIS=comma list overrides the Tier B search ROIs (default: ten damage-weighted ROIs).
"""

import hashlib
import json
import os
import pickle
import subprocess
import sys
import tempfile
import time
from fractions import Fraction

import cv2
import numpy as np
import pandas as pd
from scipy import stats as sps
from skimage.filters import threshold_multiotsu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from midog_utils import channels as ch   # colour deconvolution only (a stated dependency)
from midog_utils import dataset as ds    # TIFF reading only (a stated dependency)

T_START = time.time()
ROOT = os.path.dirname(os.path.abspath(__file__))
os.chdir(ROOT)
NB = 'production_hematoxylin_only/multiotsu_pad_arms/tm_border_floor_49roi_6seed.ipynb'
R = 'results/tm_border_floor_49roi_6seed'
OUT = 'results/tm_border_floor_49roi_6seed_audit'
PICKLE = os.path.join(tempfile.gettempdir(), 'tm_border_floor_49roi_6seed_search.pkl')
POOL_CACHE = f'{OUT}_pool_from_pickle.csv'
SUB = {'original_14': 'images/extra_valid', 'testing_35': 'images/extra_valid/testing_set'}
TIERB_DEFAULT = ['291.tiff', '401.tiff', '339.tiff', '300.tiff', '201.tiff', '546.tiff', '400.tiff', '499.tiff', '100.tiff', '462.tiff']
TIERB_ROIS = [s for s in os.environ.get('AUDIT_TIERB_ROIS', '').split(',') if s] or TIERB_DEFAULT
TP, LA, NH = 'human_correct_label', 'human_rejected_label', 'non_human_findings'
WHITE_GRAY, WHITE_SHARE, PAD, BAND = 215, 0.20, 30, 30
N_BOOT, BOOT_SEED = 10_000, 20260917
CON = {'A3 - A0': ('A3', 'A0'), 'A1 - A0': ('A1', 'A0'), 'A2 - A0': ('A2', 'A0'), 'A3 - A1': ('A3', 'A1'), 'A3 - A2': ('A3', 'A2'), 'A2m - A2': ('A2m', 'A2'), 'A1 - A1v': ('A1', 'A1v'), 'A3 - A3v': ('A3', 'A3v'), 'A2m - A0': ('A2m', 'A0')}
CHECKS = []   # (section, what, n_compared, n_divergent, tier, note)


def check(section, what, n, bad, tier, note=''):
    """
        Record one comparison block for the Part 0 table.

        section (str): audit section.
        what (str): what was compared.
        n (int): values compared.
        bad (int): divergences.
        tier (str): A, B or coherence.
        note (str): free text.

        Returns None.
    """
    CHECKS.append(dict(section=section, what=what, n_compared=int(n), n_divergent=int(bad), tier=tier, note=note))
    print(f'[{section}] {what}: {n} compared, {bad} divergent {note}')


def say(*a):
    """
        Print with flush.

        a (tuple): things to print.

        Returns None.
    """
    print(*a, flush=True)


# ---------------------------------------------------------------- 1. execution coherence and provenance
say('\n== 1. coherence and provenance')
nb = json.load(open(NB))
code = [c for c in nb['cells'] if c['cell_type'] == 'code']
ecs = [c.get('execution_count') for c in code]
n_err = sum(1 for c in code for o in c.get('outputs', []) if o['output_type'] == 'error')
contig = ecs == list(range(1, len(code) + 1))
say(f'code cells {len(code)}; execution counts contiguous from 1: {contig}; error outputs {n_err}; unrun {sum(e is None for e in ecs)}; last cell type {nb["cells"][-1]["cell_type"]}')
check('coherence', 'execution_count contiguous 1..N, no error, no unrun code cell', len(code), int(not contig) + n_err + sum(e is None for e in ecs), 'coherence')
out33 = ''.join(''.join(o.get('text', '')) for o in nb['cells'][33].get('outputs', []) if o['output_type'] == 'stream')
src19 = ''.join(nb['cells'][19]['source'])
adv = pd.read_csv(f'{R}_adversarial_log.csv')
adv_rows = []
for r in adv.itertuples():
    # the printed log truncates each column at 180 characters; compare a 150-character prefix of each field
    same = {fld: getattr(r, fld)[:150] in out33 for fld in ('concern', 'test', 'result', 'action')}
    adv_rows.append(dict(item=r.item, **{f'{k}_matches_printed': v for k, v in same.items()}))
ADVCMP = pd.DataFrame(adv_rows)
n_adv_bad = int((~ADVCMP[[c for c in ADVCMP if c.endswith('_printed')]]).sum().sum())
check('coherence', 'adversarial_log.csv fields vs the executed cell-33 printout', ADVCMP.shape[0] * 4, n_adv_bad, 'coherence', '(post-run hand patch, admitted in cell 34)')
say(ADVCMP[~ADVCMP.drop(columns='item').all(axis=1)].to_string(index=False))
say('cell 19 source now says "its 273 is edge_px <= 2":', "its 273 is edge_px <= 2" in src19, '; executed cell-33 output still says "matches neither exactly":', 'matches neither exactly' in out33)
mt = {p: os.path.getmtime(p) for p in [NB, f'{R}_per_run.csv', f'{R}_adversarial_log.csv', f'{R}_delta_stats.csv']}
say('mtimes:', {k: time.strftime('%H:%M:%S', time.localtime(v)) for k, v in mt.items()})
git = lambda *a: subprocess.run(['git', *a], capture_output=True, text=True).stdout.strip()
say('git status of midog_utils (empty = clean):', repr(git('status', '--porcelain', 'midog_utils')))
say('git status of the notebook and its outputs:', git('status', '--porcelain', NB, *[f'{R}_{n}.csv' for n in ('per_run', 'top30', 'delta_stats')]).replace('\n', ' | '))
say('last commit touching midog_utils:', git('log', '-1', '--format=%h %ad %s', '--', 'midog_utils'))
inputs = ['results/precision_at_k_49roi_6seed_fillgeom_postD11_clicks.csv', 'results/precision_at_k_49roi_6seed_fillgeom_postD11_per_run.csv', 'results/precision_at_k_49roi_6seed_fillgeom_postD11_delta_stats.csv', 'results/deep_pool_mo_pad61_49roi_6seed_postD11_per_run.csv', 'results/deep_pool_mo_pad61_49roi_6seed_postD11_top1000.csv', 'results/precision_at_k_49roi_6seed_mopad_postD11_templates.csv', 'results/precision_at_k_49roi_6seed_mopad_postD11_per_run.csv', 'results/precision_at_k_49roi_6seed_mopad_postD11_top30.csv', 'databases/MIDOG++.json']
out2 = ''.join(''.join(o.get('text', '')) for o in nb['cells'][2].get('outputs', []) if o['output_type'] == 'stream')
sha_ok = [hashlib.sha256(open(p, 'rb').read()).hexdigest()[:12] in out2 for p in inputs]
check('provenance', 'input SHA-256 prefixes vs those printed in cell 2', len(inputs), len(inputs) - sum(sha_ok), 'coherence')

# ---------------------------------------------------------------- 2. composition gate
say('\n== 2. composition')
pr = pd.read_csv(f'{R}_per_run.csv')
top = pd.read_csv(f'{R}_top30.csv')
raw = json.load(open('databases/MIDOG++.json'))
IMGS = {im['file_name']: im for im in raw['images']}
BY_IMG = {}
for a in raw['annotations']:
    BY_IMG.setdefault(a['image_id'], []).append(a)
ORDER = list(pd.unique(pd.read_csv('results/precision_at_k_49roi_6seed_fillgeom_postD11_per_run.csv', usecols=['file_name'])['file_name']))
DOM = pr.drop_duplicates('file_name').set_index('file_name')['domain']
json_dom = {f: IMGS[f]['tumor_type'].replace('canine lymphoma', 'canine lymphosarcoma') for f in ORDER}
comp = dict(rows=len(pr) == 49 * 6 * 7 * 2, dup_keys=not pr.duplicated(['file_name', 'seed_index', 'arm', 'rank_key']).any(), rois=pr.file_name.nunique() == 49, clicks6=(pr.groupby('file_name').seed_index.nunique() == 6).all(), domains7x7=(pr.groupby('domain').file_name.nunique() == 7).all() and pr.domain.nunique() == 7, domain_matches_json=all(DOM[f] == json_dom[f] for f in ORDER), top30_rows=len(top) == 4116 * 30, roi_set_equals_order=set(ORDER) == set(pr.file_name))
check('composition', 'per_run 4116 rows, unique keys, 49 ROIs x 6 clicks, 7 domains x 7 ROIs (domain = JSON tumor_type), top30 = 4116 x 30', len(comp), sum(not v for v in comp.values()), 'A', str(comp))
roi_per_type = pd.Series(json_dom).value_counts()
say('ROIs per tumour type in the 49 (stratum vs exchangeable unit):', roi_per_type.to_dict())

# ---------------------------------------------------------------- 3. independent re-match of every top-30 list
say('\n== 3. independent re-match of all 4116 top-30 lists against MIDOG++.json')


def gt_of(f, seed_ann):
    """
        One ROI's ground truth from the raw JSON, both categories, the seed annotation dropped.

        f (str): file name.
        seed_ann (int): the click's own annotation id.

        Returns tuple[np.ndarray, np.ndarray, np.ndarray]: ids, bbox centres (TLBR midpoints), categories.
    """
    anns = [a for a in BY_IMG[IMGS[f]['id']] if a['id'] != seed_ann]
    return np.array([a['id'] for a in anns]), np.array([((a['bbox'][0] + a['bbox'][2]) / 2.0, (a['bbox'][1] + a['bbox'][3]) / 2.0) for a in anns]), np.array([a['category_id'] for a in anns])


def greedy(xy, gid, gxy, gcat, radius):
    """
        Rank-order greedy one-to-one match: each detection takes the nearest free GT within radius (inclusive).

        xy (np.ndarray): (n, 2) detections best-first.
        gid (np.ndarray): GT ids.
        gxy (np.ndarray): GT centres.
        gcat (np.ndarray): GT categories.
        radius (float): match radius in px.

        Returns tuple[np.ndarray, np.ndarray]: matched category (0 = none) and matched id (-1 = none) per detection.
    """
    taken = np.zeros(len(gid), bool)
    cat, mid = np.zeros(len(xy), int), np.full(len(xy), -1)
    for i, (x, y) in enumerate(xy):
        d = np.hypot(gxy[:, 0] - x, gxy[:, 1] - y)
        ok = np.where((d <= radius) & ~taken)[0]
        if len(ok):
            g = ok[np.argmin(d[ok])]
            taken[g] = True
            cat[i], mid[i] = gcat[g], gid[g]
    return cat, mid


KEY = ['file_name', 'seed_index', 'arm', 'rank_key']
PRK = pr.set_index(KEY)
rows, n_bucket_bad, n_mid_bad, n_order_bad, n_edge_bad = [], 0, 0, 0, 0
for k, g in top.groupby(KEY, sort=False):
    g = g.sort_values('rank')
    r = PRK.loc[k]
    W, H = IMGS[k[0]]['width'], IMGS[k[0]]['height']
    gid, gxy, gcat = gt_of(k[0], r.seed_ann_id)
    rad = 7.5 / r.mpp
    cat, mid = greedy(g[['cx', 'cy']].to_numpy(), gid, gxy, gcat, rad)
    n_bucket_bad += int((np.where(cat == 1, TP, np.where(cat == 2, LA, NH)) != g.bucket.to_numpy()).sum())
    n_mid_bad += int((mid != g.matched_ann_id.to_numpy()).sum())
    n_order_bad += int((np.diff(g['score' if k[3] == 'tm_score' else 'od'].to_numpy()) > 0).sum())
    x, y = g.cx.to_numpy(), g.cy.to_numpy()
    e = np.minimum.reduce([x, y, W - 1 - x, H - 1 - y])
    n_edge_bad += int((e != g.edge_px.to_numpy()).sum())
    tp, wh = cat == 1, g.white_frac_61.to_numpy() > WHITE_SHARE
    dist = np.hypot(x - r.click_cx, y - r.click_cy)
    near = ~tp & (dist >= rad) & (dist <= 2 * rad)
    rows.append(dict(zip(KEY, k), tp_at_10=int(tp[:10].sum()), tp_at_20=int(tp[:20].sum()), tp_at_30=int(tp[:30].sum()), n_fp_top10=int((~tp[:10]).sum()), n_white_fp_top10=int((~tp & wh)[:10].sum()), n_tp_top30=int(tp[:30].sum()), n_white_tp_top30=int((tp & wh)[:30].sum()), n_fp_top30=int((~tp[:30]).sum()), n_white_fp_top30=int((~tp & wh)[:30].sum()), n_edge30_top10=int((e[:10] < BAND).sum()), n_nearseed_fp_top10=int(near[:10].sum()), n_nearseed_fp_top30=int(near[:30].sum())))
RE = pd.DataFrame(rows).set_index(KEY)
check('rematch', 'top-30 buckets re-derived from JSON (own greedy matcher)', len(top), n_bucket_bad, 'A')
check('rematch', 'top-30 matched_ann_id re-derived', len(top), n_mid_bad, 'A')
check('rematch', 'top-30 rank order non-increasing in its own axis key', len(top) - 4116, n_order_bad, 'A')
check('rematch', 'top-30 edge_px re-derived from JSON width/height', len(top), n_edge_bad, 'A')
cmpc = list(RE.columns)
bad = (PRK.loc[RE.index, cmpc].to_numpy() != RE[cmpc].to_numpy())
check('rematch', 'per_run TP@10/20/30 and top-10/30 white/edge/near-seed counts from the re-match', bad.size, int(bad.sum()), 'A')
RE.reset_index().to_csv(f'{OUT}_rematch_per_run.csv', index=False)

# ---------------------------------------------------------------- 4. the pool (pickle-only evidence)
say('\n== 4. full post-NMS pools from the search pickle')
if os.path.exists(PICKLE):
    with open(PICKLE, 'rb') as fh:
        PK = pickle.load(fh)
    num = [c for c in pr.columns if pr[c].dtype.kind in 'if']
    same = np.isclose(PK['PER_RUN'][num].to_numpy(float), pr[num].to_numpy(float), equal_nan=True, rtol=1e-12, atol=0)
    check('pickle', 'pickle PER_RUN numeric columns vs per_run.csv', same.size, int((~same).sum()), 'A', '(the pickle is the run the CSV came from)')
    POOL = PK['POOL']
    clk = pr[(pr.arm == 'A0') & (pr.rank_key == 'tm_score')].set_index(['file_name', 'seed_index'])
    parts, pb = [], 0
    for (f, s, arm), g in POOL.groupby(['file_name', 'seed_index', 'arm'], sort=False):
        g = g.sort_values('rank')
        c = clk.loc[(f, s)]
        gid, gxy, gcat = gt_of(f, c.seed_ann_id)
        cat, mid = greedy(g[['cx', 'cy']].to_numpy(), gid, gxy, gcat, 7.5 / c.mpp)
        pb += int((np.where(cat == 1, TP, np.where(cat == 2, LA, NH)) != g.bucket.to_numpy()).sum())
        parts.append(g.assign(cat=cat, mid=mid))
    PM = pd.concat(parts, ignore_index=True)
    check('pickle', 'pool buckets (tm_score order) re-derived from JSON', len(PM), pb, 'A')
    PM.to_csv(POOL_CACHE, index=False)
    PK['LEAK'].to_csv(f'{OUT}_leak_from_pickle.csv', index=False)
    POOL_SOURCE = f'{PICKLE} (re-matched, cached to {POOL_CACHE})'
else:
    PM = pd.read_csv(POOL_CACHE)
    POOL_SOURCE = f'{POOL_CACHE} (the temp pickle no longer exists; this is the audit-time copy)'
say('pool source:', POOL_SOURCE)
clk = pr[(pr.arm == 'A0') & (pr.rank_key == 'tm_score')].set_index(['file_name', 'seed_index'])
pg = PM.groupby(['file_name', 'seed_index', 'arm']).agg(n_detections=('cx', 'size'), n_mit_pool=('cat', lambda c: int((c == 1).sum())), n_lookalike_pool=('cat', lambda c: int((c == 2).sum())), n_edge0=('edge_px', lambda e: int((e == 0).sum())), n_edge2=('edge_px', lambda e: int((e <= 2).sum())), n_edge30=('edge_px', lambda e: int((e < BAND).sum())))
near_pool = {}
for (f, s, arm), g in PM.groupby(['file_name', 'seed_index', 'arm']):
    c = clk.loc[(f, s)]
    rad = 7.5 / c.mpp
    dd = np.hypot(g.cx - c.click_cx, g.cy - c.click_cy)
    near_pool[(f, s, arm)] = int(((g.cat != 1) & (dd >= rad) & (dd <= 2 * rad)).sum())
pg['n_nearseed_fp_pool'] = pd.Series(near_pool)
tmk = pr[pr.rank_key == 'tm_score'].set_index(['file_name', 'seed_index', 'arm'])[list(pg.columns)]
bad = tmk.to_numpy() != pg.reindex(tmk.index).to_numpy()
check('pickle', 'per_run pool columns (n_det, pool mitoses/look-alikes, edge 0/<=2/<30, near-seed) from the pool', bad.size, int(bad.sum()), 'A')
a0p = PM[PM.arm == 'A0']
EDGE_LT2 = int((a0p.edge_px < 2).sum())
say(f'A0 pools: edge_px == 0 {int((a0p.edge_px == 0).sum())}, < 2 {EDGE_LT2}, <= 2 {int((a0p.edge_px <= 2).sum())}')

# ---------------------------------------------------------------- 5. statistics, re-implemented
say('\n== 5. delta_stats re-derived (own exact sign-flip by dictionary DP, own Holm, t on 48 df; bootstrap replicated with the notebook seed = consistency)')


def signflip_p(d):
    """
        Exact two-sided sign-flip p over the nonzero integer deltas, by a dictionary of attainable sums.

        d (array-like): one integer delta per unit.

        Returns tuple[float, int]: p, number of nonzero units.
    """
    v = [abs(int(x)) for x in d if int(x) != 0]
    if not v:
        return 1.0, 0
    dist = {0: 1}
    for a in v:
        nd = {}
        for s, c in dist.items():
            nd[s + a] = nd.get(s + a, 0) + c
            nd[s - a] = nd.get(s - a, 0) + c
        dist = nd
    obs = abs(int(np.sum(d)))
    return float(Fraction(sum(c for s, c in dist.items() if abs(s) >= obs), 2 ** len(v))), len(v)


def holm(p):
    """
        Holm step-down adjusted p-values.

        p (array-like): raw p-values.

        Returns np.ndarray: adjusted, input order.
    """
    p = np.asarray(p, float)
    o = np.argsort(p, kind='stable')
    out, run = np.empty(len(p)), 0.0
    for i, j in enumerate(o):
        run = max(run, min(1.0, (len(p) - i) * p[j]))
        out[j] = run
    return out


IDX = np.random.default_rng([BOOT_SEED, 0]).integers(49, size=(N_BOOT, 49))
ROI = pr.groupby(['rank_key', 'arm', 'file_name'])[['tp_at_10', 'tp_at_20', 'tp_at_30']].sum()


def roi_vec(axis, arm, col):
    """
        One per-ROI click-summed column in the bootstrap ROI order.

        axis (str): rank_key.
        arm (str): arm.
        col (str): column.

        Returns np.ndarray: 49 integers.
    """
    return ROI.loc[(axis, arm)].reindex(ORDER)[col].to_numpy(np.int64)


dstats = pd.read_csv(f'{R}_delta_stats.csv')
rows = []
for r in dstats.itertuples():
    a, b = CON[r.pair]
    d = roi_vec(r.rank_key, a, f'tp_at_{r.K}') - roi_vec(r.rank_key, b, f'tp_at_{r.K}')
    den = 6 * 49 * r.K
    p, n = signflip_p(d)
    lo, hi = np.percentile(d[IDX].sum(1) / den, [2.5, 97.5])
    x = d / (6 * r.K)
    half = sps.t.ppf(0.975, 48) * x.std(ddof=1) / np.sqrt(49)
    w, l, t = int((d > 0).sum()), int((d < 0).sum()), int((d == 0).sum())
    pooled = d.sum() / den
    maj = (w if pooled > 0 else l if pooled < 0 else 0) > 24.5
    dd = pd.Series(d, index=ORDER).groupby(DOM.reindex(ORDER).to_numpy()).sum()
    pdom, ndom = signflip_p(dd.to_numpy())
    rows.append(dict(holm_family=r.holm_family, rank_key=r.rank_key, pair=r.pair, K=r.K, pooled=pooled, wins=w, losses=l, ties=t, exact_p=p, n_nonzero=n, ci_low=lo, ci_high=hi, t_low=x.mean() - half, t_high=x.mean() + half, d5=bool((lo > 0 or hi < 0) and maj), domains_pos=int((dd > 0).sum()), domains_neg=int((dd < 0).sum()), domains_zero=int((dd == 0).sum()), domain_signflip_p=pdom, domain_floor=2.0 ** (1 - ndom) if ndom else 1.0))
ST = pd.DataFrame(rows)
tested = ST.holm_family != 'descriptive'
ST.loc[tested, 'holm_p'] = ST[tested].groupby('holm_family')['exact_p'].transform(lambda s: holm(s.to_numpy()))
mp = {'pooled': 'pooled_delta_precision', 'wins': 'wins', 'losses': 'losses', 'ties': 'ties', 'exact_p': 'exact_p', 'n_nonzero': 'n_nonzero', 'ci_low': 'ci_low', 'ci_high': 'ci_high', 't_low': 't_ci_low', 't_high': 't_ci_high', 'd5': 'd5_bar_met', 'holm_p': 'holm_p'}
nb_bad = 0
for m, t_ in mp.items():
    a_, b_ = ST[m].to_numpy(), dstats[t_].to_numpy()
    if m == 'd5':
        nb_bad += int((a_.astype(bool) != b_.astype(bool)).sum())
    else:
        nb_bad += int((~(np.isclose(a_.astype(float), b_.astype(float), rtol=1e-9, atol=1e-15) | (np.isnan(a_.astype(float)) & np.isnan(b_.astype(float))))).sum())
check('stats', 'delta_stats: pooled, W/L/T, exact p, n_nonzero, bootstrap CI, t CI, D5 bar, Holm (54 rows x 12)', len(dstats) * len(mp), nb_bad, 'A', '(bootstrap columns are consistency: same seed and resample)')
# primary family under a domain-level null, Holm over its 3 K
prim = ST[(ST.holm_family == 'primary')]
say('primary A3 - A0 tm: domains +/-/0 per K', prim[['K', 'domains_pos', 'domains_neg', 'domains_zero', 'domain_signflip_p']].to_dict('records'), '; Holm over 3 K at domain level:', holm(prim.domain_signflip_p.to_numpy()).round(4).tolist())
ST.to_csv(f'{OUT}_delta_stats_rederived.csv', index=False)

# ---------------------------------------------------------------- 6. M1-M4, pooled P, generator, spread, domains, adversarial numbers
say('\n== 6. M1-M4 and every printed table')
TM = pr[pr.rank_key == 'tm_score']


def rsum(axis, arm, col):
    """
        A per_run column summed per ROI, bootstrap order.

        axis (str): rank_key.
        arm (str): arm.
        col (str): column.

        Returns np.ndarray: 49 values.
    """
    return pr[(pr.rank_key == axis) & (pr.arm == arm)].groupby('file_name')[col].sum().reindex(ORDER).to_numpy()


white = pd.read_csv(f'{R}_white.csv')
wb = 0
wn = 0
for axis in ('tm_score', 'chromatin_od'):
    for a in ('A3', 'A2', 'A1', 'A2m', 'A3v'):
        na, da, nb_, db = rsum(axis, a, 'n_white_fp_top10'), rsum(axis, a, 'n_fp_top10'), rsum(axis, 'A0', 'n_white_fp_top10'), rsum(axis, 'A0', 'n_fp_top10')
        diff = na[IDX].sum(1) / da[IDX].sum(1) - nb_[IDX].sum(1) / db[IDX].sum(1)
        mine = [na.sum() / da.sum(), nb_.sum() / db.sum(), *np.percentile(diff, [2.5, 97.5])]
        w = white[(white.check == 'M2') & (white.rank_key == axis) & (white.contrast == f'{a} - A0')].iloc[0]
        wb += int((~np.isclose(mine, [w.share_a, w.share_b, w.ci_low, w.ci_high], rtol=1e-12)).sum())
        wn += 4
        dm3 = rsum(axis, a, 'n_white_tp_top30') - rsum(axis, 'A0', 'n_white_tp_top30')
        p3, _ = signflip_p(dm3)
        lo3, hi3 = np.percentile(dm3[IDX].sum(1), [2.5, 97.5])
        w3 = white[(white.check == 'M3') & (white.rank_key == axis) & (white.contrast == f'{a} - A0')].iloc[0]
        mine3 = [rsum(axis, a, 'n_white_tp_top30').sum(), rsum(axis, 'A0', 'n_white_tp_top30').sum(), dm3.sum(), (dm3 > 0).sum(), (dm3 < 0).sum(), (dm3 == 0).sum(), p3, lo3, hi3]
        wb += int((~np.isclose(mine3, [w3.count_a, w3.count_b, w3.total, w3.wins, w3.losses, w3.ties, w3.exact_p, w3.ci_low, w3.ci_high], rtol=1e-12)).sum())
        wn += 9
m3h = white[(white.check == 'M3') & (white.rank_key == 'tm_score') & white.contrast.isin(['A2 - A0', 'A3 - A0'])]
wb += int((~np.isclose(holm(m3h.exact_p.to_numpy()), m3h.holm_p.to_numpy())).sum())
wn += 2
check('M2/M3', 'white.csv: M2 shares + bootstrap CIs, M3 counts/total/WLT/exact p/CI, M3 Holm', wn, wb, 'A')
nsd = pd.read_csv(f'{R}_nearseed.csv')
nbad, nn = 0, 0
for r in nsd.itertuples():
    col = f'n_nearseed_fp_{r.scope}'
    if ' - ' not in r.contrast:
        nbad += int(rsum('tm_score', r.contrast, col).sum() != r.total)
        nn += 1
    else:
        a, b = r.contrast.split(' - ')
        d = rsum('tm_score', a, col) - rsum('tm_score', b, col)
        p, _ = signflip_p(d)
        nbad += int(d.sum() != r.total) + int(not np.isclose(p, r.exact_p))
        nn += 2
check('M4', 'nearseed.csv totals and deltas (pool counts from the re-matched pool, top 10/30 from the re-match)', nn, nbad, 'A')
edge = pd.read_csv(f'{R}_edge.csv').set_index('arm')
eg = TM.groupby('arm')[['n_edge0', 'n_edge2', 'n_edge30', 'n_edge30_top10', 'n_detections']].sum()
eg['od10'] = pr[pr.rank_key == 'chromatin_od'].groupby('arm')['n_edge30_top10'].sum()
eb = (eg[['n_edge0', 'n_edge2', 'n_edge30', 'n_edge30_top10', 'od10', 'n_detections']].to_numpy() != edge.loc[eg.index, ['n_edge0', 'n_edge2', 'n_edge30', 'n_edge30_top10_tm', 'n_edge30_top10_od', 'n_detections']].to_numpy())
check('M1', 'edge.csv (edge 0, <=2, <30, top-10 band slots both axes, detections) per arm', eb.size, int(eb.sum()), 'A')
POOLED = {(axis, a): [rsum(axis, a, f'tp_at_{k}').sum() / (294 * k) for k in (10, 20, 30)] for axis in ('tm_score', 'chromatin_od') for a in ('A0', 'A1', 'A1v', 'A2', 'A3', 'A3v', 'A2m')}
out23 = ''.join(''.join(o.get('text', '')) for o in nb['cells'][23].get('outputs', []) if o['output_type'] == 'stream')
pooled_printed = {}
for line in out23.split('pooled precision per arm (all clicks):')[1].split('decision rule')[0].strip().split('\n')[3:]:
    t_ = line.split()
    if t_ and t_[0] in ('A0', 'A1', 'A1v', 'A2', 'A3', 'A3v', 'A2m'):
        v = list(map(float, t_[1:]))
        pooled_printed[t_[0]] = {'chromatin_od': v[0::2], 'tm_score': v[1::2]}
pb_ = sum(int(abs(round(POOLED[(ax, a)][i], 4) - pooled_printed[a][ax][i]) > 6e-5) for a in pooled_printed for ax in ('tm_score', 'chromatin_od') for i in range(3))
check('pooled', 'pooled P@10/20/30 per arm and axis vs cell-23 printout (4 dp)', 42, pb_, 'A')
# per_roi.csv
proi = pd.read_csv(f'{R}_per_roi.csv').set_index(['rank_key', 'arm', 'file_name'])
g = pr.groupby(['rank_key', 'arm', 'file_name']).agg(n_clicks=('tp_at_10', 'size'), tp_at_10=('tp_at_10', 'sum'), tp_at_20=('tp_at_20', 'sum'), tp_at_30=('tp_at_30', 'sum'), tp10_click_sd=('tp_at_10', 'std'), tp10_worst_click=('tp_at_10', 'min'), n_mit_pool=('n_mit_pool', 'sum'), n_white_fp_top10=('n_white_fp_top10', 'sum'), n_fp_top10=('n_fp_top10', 'sum'), n_white_tp_top30=('n_white_tp_top30', 'sum'))
bad = ~np.isclose(g.to_numpy(float), proi.loc[g.index, g.columns].to_numpy(float), rtol=1e-12)
check('per_roi', 'per_roi.csv click sums, SD, worst click, pool mitoses, white counts', bad.size, int(bad.sum()), 'A')

# ---------------------------------------------------------------- 7. claims that exist only in cell 34 prose (and in cell outputs), recomputed
say('\n== 7. numbers stated in the closing cell, recomputed')
P = []


def claim(cid, text, nbv, mine, tol=0.0):
    """
        Record one prose number against its recomputation.

        cid (str): claim id.
        text (str): what the claim says.
        nbv (float or str): the notebook's value.
        mine (float or str): the recomputed value.
        tol (float): absolute tolerance for numeric values (rounding in prose).

        Returns None.
    """
    ok = (nbv == mine) if isinstance(nbv, str) else bool(abs(float(nbv) - float(mine)) <= tol + 1e-12)
    P.append(dict(claim=cid, text=text, notebook=nbv, recomputed=mine, match=ok))


def con(axis, a, b, k):
    """
        Per-ROI integer TP@K delta for one contrast.

        axis (str): rank_key.
        a (str): arm.
        b (str): arm.
        k (int): budget.

        Returns np.ndarray: 49 deltas.
    """
    return rsum(axis, a, f'tp_at_{k}') - rsum(axis, b, f'tp_at_{k}')


d10 = con('tm_score', 'A3', 'A0', 10)
claim('C2', 'A3-A0 tm K=10 pts', 3.71, 100 * d10.sum() / 2940, 0.005)
claim('C2', 'P@10 A0', 0.542, POOLED[('tm_score', 'A0')][0], 0.0005)
claim('C2', 'P@10 A3', 0.579, POOLED[('tm_score', 'A3')][0], 0.0005)
claim('C3', 'od A3-A0 K=10/20/30 pts', '+1.80/+0.97/+0.78', '/'.join(f'{100 * con("chromatin_od", "A3", "A0", k).sum() / (294 * k):+.2f}' for k in (10, 20, 30)))
claim('C6', 'A2-A0 tm K=10 pts', 3.20, 100 * con('tm_score', 'A2', 'A0', 10).sum() / 2940, 0.005)
claim('C6', 'A1-A0 tm K=10 pts', 0.41, 100 * con('tm_score', 'A1', 'A0', 10).sum() / 2940, 0.005)
od_band = top[(top.arm == 'A0') & (top.rank_key == 'chromatin_od') & (top['rank'] < 10) & (top.edge_px < BAND)]
claim('C7', 'A0 od top-10 band entries', 174, len(od_band))
claim('C7', 'share of them unlabelled', 0.76, (od_band.bucket == NH).mean(), 0.005)
tm_band = top[(top.arm == 'A0') & (top.rank_key == 'tm_score') & (top['rank'] < 10) & (top.edge_px < BAND)]
claim('C22', 'A0 band mitoses in tm top 10s', 17, int((tm_band.bucket == TP).sum()))
claim('C22', 'A0 band slots in tm top 10s', 75, len(tm_band))
claim('C22', 'A0 tm top-30 band TPs (adv log)', 44, int(((top.arm == 'A0') & (top.rank_key == 'tm_score') & (top.edge_px < BAND) & (top.bucket == TP)).sum()))
claim('C9', 'A0 top-10 FPs (tm)', 1347, int(rsum('tm_score', 'A0', 'n_fp_top10').sum()))
claim('C9', 'A3 top-10 FPs (tm)', 1238, int(rsum('tm_score', 'A3', 'n_fp_top10').sum()))
claim('C9', 'A0 white-adjacent top-10 FPs', 322, int(rsum('tm_score', 'A0', 'n_white_fp_top10').sum()))
claim('C9', 'A3 white-adjacent top-10 FPs', 240, int(rsum('tm_score', 'A3', 'n_white_fp_top10').sum()))
mast = TM[(TM.arm == 'A0') & (TM.domain == 'canine cutaneous mast cell tumor')]
claim('C20', 'mast cell A0 white share of top-10 non-mitoses', 0.82, mast.n_white_fp_top10.sum() / mast.n_fp_top10.sum(), 0.005)
dd = pd.Series(d10, index=ORDER)
claim('C19', 'largest single-ROI gain (291.tiff), TP over 6 clicks', 13, int(dd['291.tiff']))
claim('C19', 'is it the largest', '291.tiff', dd.idxmax())
claim('C19', 'pooled K=10 gain without 291 (pts)', 3.33, 100 * (dd.sum() - dd['291.tiff']) / (48 * 60), 0.005)
roi_sd = pr.groupby(['rank_key', 'arm', 'file_name'])['tp_at_10'].agg(['std', 'min'])
claim('C19', 'mean within-ROI SD TP@10 A0 (tm)', 1.03, roi_sd.loc[('tm_score', 'A0')]['std'].mean(), 0.005)
claim('C19', 'mean within-ROI SD TP@10 A3 (tm)', 0.99, roi_sd.loc[('tm_score', 'A3')]['std'].mean(), 0.005)
claim('C19', 'mean worst click TP@10 A0', 4.06, roi_sd.loc[('tm_score', 'A0')]['min'].mean(), 0.005)
claim('C19', 'mean worst click TP@10 A3', 4.55, roi_sd.loc[('tm_score', 'A3')]['min'].mean(), 0.005)
z0 = sorted(roi_sd.loc[('tm_score', 'A0')].index[roi_sd.loc[('tm_score', 'A0')]['min'] == 0])
z3 = sorted(roi_sd.loc[('tm_score', 'A3')].index[roi_sd.loc[('tm_score', 'A3')]['min'] == 0])
claim('C19', 'count of ROIs whose worst click finds 0, A0 -> A3', '4 -> 4', f'{len(z0)} -> {len(z3)}')
claim('C19', '... and are they the same four ROIs', 'implied same', 'same' if z0 == z3 else f'A0 {z0} vs A3 {z3}')
pm = {a: rsum('tm_score', a, 'n_mit_pool') for a in ('A0', 'A1', 'A2', 'A3')}
claim('C18', 'pool mitoses per click A0', 24.67, pm['A0'].sum() / 294, 0.005)
claim('C18', 'pool mitoses per click A3', 26.21, pm['A3'].sum() / 294, 0.005)
claim('C18', 'pool mitoses per click A2', 26.15, pm['A2'].sum() / 294, 0.005)
claim('C18', 'pool mitoses per click A1', 24.68, pm['A1'].sum() / 294, 0.005)
dp = pm['A3'] - pm['A0']
claim('C18', 'A3-A0 pool mitoses W/L/T', '44/4/1', f'{(dp > 0).sum()}/{(dp < 0).sum()}/{(dp == 0).sum()}')
claim('C18', 'A3-A0 pool mitoses exact p', 4e-10, signflip_p(dp)[0], 0.5e-10)
a0c = TM[TM.arm == 'A0'].set_index(['file_name', 'seed_index'])
a2c = TM[TM.arm == 'A2'].set_index(['file_name', 'seed_index'])
a1c = TM[TM.arm == 'A1'].set_index(['file_name', 'seed_index'])
bf = a2c.floor / a0c.fill_value
claim('C10', 'clicks with b_otsu >= own template fill', 156, int((bf >= 1).sum()))
claim('C10', 'median b_otsu / fill', 1.03, bf.median(), 0.005)
claim('C10', 'range b_otsu / fill', '0.47-3.33', f'{bf.min():.2f}-{bf.max():.2f}')
claim('C10', 'clipped / recorded fill, median', 1.30, (a2c.fill_value / a0c.fill_value).median(), 0.005)
claim('C10', 'clipped / recorded fill, max', 3.40, (a2c.fill_value / a0c.fill_value).max(), 0.005)
claim('C10', 'clipped / recorded mass weight, median', 0.89, (a2c.mass_weight / a0c.mass_weight).median(), 0.005)
claim('C22', 'mass_weight <= 0: clipped templates / A0 templates', '2/3', f"{int(a2c.mass_weight_nonpositive.sum())}/{int(a0c.mass_weight_nonpositive.sum())}")
claim('C22', 'threshold mean A0', -0.956, a0c.threshold.mean(), 0.0005)
claim('C22', 'threshold mean A2', -0.595, a2c.threshold.mean(), 0.0005)
claim('C22', 'MAX_PEAKS binds on every click of every arm', 'yes', 'yes' if TM.max_peaks_binding.all() and (TM.n_peaks == 100).all() else 'no')
mr, Mr = a1c['median'] / a0c['median'], a1c['mad'] / a0c['mad']
claim('C22', 'A1/A0 median ratio median (range)', '1.004 (0.94-1.03)', f'{mr.median():.3f} ({mr.min():.2f}-{mr.max():.2f})')
claim('C22', 'A1/A0 MAD ratio median (range)', '1.001 (1.00-1.14)', f'{Mr.median():.3f} ({Mr.min():.2f}-{Mr.max():.2f})')
claim('C15', 'additivity gap K=10 (pts)', 0.10, 100 * (con('tm_score', 'A3', 'A0', 10) - con('tm_score', 'A1', 'A0', 10) - con('tm_score', 'A2', 'A0', 10)).sum() / 2940, 0.005)
claim('C15', 'A3-A1 tm K=10 (pts)', 3.30, 100 * con('tm_score', 'A3', 'A1', 10).sum() / 2940, 0.005)
claim('C16', 'A2-A0 od K=10 (pts)', -0.34, 100 * con('chromatin_od', 'A2', 'A0', 10).sum() / 2940, 0.005)
desc = ST[ST.holm_family == 'descriptive']
claim('C23', 'max |A1-A1v|, |A3-A3v| over K and axis (pts) <= 0.1', 'yes', 'yes' if (desc[desc.pair.isin(['A1 - A1v', 'A3 - A3v'])].pooled.abs() * 100).max() <= 0.1 else 'no')
mf = ST[ST.holm_family == 'median_floor']
claim('C10', 'A2m vs A2: "all Holm p >= 0.21"', 'min Holm p 0.21', f'min Holm p {mf.holm_p.min():.2f}; min raw p {mf.exact_p.min():.2f}')
fl = pd.read_csv(f'{R}_floors.csv').set_index('file_name')
claim('C10', 'ROIs with b_otsu clip share > 60 %', 20, int((fl.clip_share_otsu > 0.6).sum()))
claim('C10', 'white share 529.tiff', 0.01, fl.loc['529.tiff', 'white_share_roi'], 0.005)
claim('C10', 'white share 401.tiff', 0.55, fl.loc['401.tiff', 'white_share_roi'], 0.005)
cv = fl.groupby('domain').b_otsu.agg(lambda v: v.std() / v.mean())
claim('C22', 'b_otsu CV within domain, min (mast) .. max (neuro)', '0.14..0.44', f'{cv.min():.2f}..{cv.max():.2f}')
claim('C22', 'every hem.min() is 0', 'yes', 'yes' if (fl.hem_min == 0).all() else 'no')
mit = [(a, IMGS[f]) for f in ORDER for a in BY_IMG[IMGS[f]['id']] if a['category_id'] == 1]
me = [min((a['bbox'][0] + a['bbox'][2]) / 2, (a['bbox'][1] + a['bbox'][3]) / 2, im['width'] - 1 - (a['bbox'][0] + a['bbox'][2]) / 2, im['height'] - 1 - (a['bbox'][1] + a['bbox'][3]) / 2) for a, im in mit]
claim('C22', 'annotated mitoses within 30 px of the edge / all, 49 ROIs', '28/2782', f'{int(np.sum(np.array(me) < 30))}/{len(me)}')
claim('C22', 'edge_px < 2 in A0 pools (cell-1 concern text)', 264, EDGE_LT2)
sp_w = sps.spearmanr(fl.loc[ORDER, 'white_share_roi'], con('tm_score', 'A2', 'A0', 10))
sp_c = sps.spearmanr(fl.loc[ORDER, 'clip_share_otsu'], con('tm_score', 'A2', 'A0', 10))
claim('C9', 'Spearman white share vs A2-A0 (rho, p)', '+0.08 (0.60)', f'{sp_w.statistic:+.2f} ({sp_w.pvalue:.2f})')
claim('C9', 'Spearman clip share vs A2-A0 (rho, p)', '-0.26 (0.08)', f'{sp_c.statistic:+.2f} ({sp_c.pvalue:.2f})')
# per-domain table (tm_score, P@10)
DT = []
for dname in pd.unique(DOM.reindex(ORDER)):
    fs = [f for f in ORDER if DOM[f] == dname]
    row = dict(domain=dname, clip_share=fl.loc[fs, 'clip_share_otsu'].mean(), white_share=fl.loc[fs, 'white_share_roi'].mean())
    for a in ('A0', 'A1', 'A2', 'A3', 'A2m'):
        row[a] = ROI.loc[('tm_score', a)].loc[fs, 'tp_at_10'].sum() / (60 * len(fs))
    dlt = ROI.loc[('tm_score', 'A3')].loc[fs, 'tp_at_10'] - ROI.loc[('tm_score', 'A0')].loc[fs, 'tp_at_10']
    row['A3-A0 pts'] = 100 * (row['A3'] - row['A0'])
    row['A2-A0 pts'] = 100 * (row['A2'] - row['A0'])
    row['WLT'] = f'{(dlt > 0).sum()}/{(dlt < 0).sum()}/{(dlt == 0).sum()}'
    row['od A3-A0 pts'] = 100 * (ROI.loc[('chromatin_od', 'A3')].loc[fs, 'tp_at_10'].sum() - ROI.loc[('chromatin_od', 'A0')].loc[fs, 'tp_at_10'].sum()) / (60 * len(fs))
    row['od A1-A0 pts'] = 100 * (ROI.loc[('chromatin_od', 'A1')].loc[fs, 'tp_at_10'].sum() - ROI.loc[('chromatin_od', 'A0')].loc[fs, 'tp_at_10'].sum()) / (60 * len(fs))
    DT.append(row)
DT = pd.DataFrame(DT)
say(DT.to_string(index=False, float_format=lambda v: f'{v:.3f}'))
DT.to_csv(f'{OUT}_domain_table.csv', index=False)
nb_dom = {'canine cutaneous mast cell tumor': (0.64, 0.215, 0.648, 0.652, 0.664, 0.669, 0.660, 2.1, '5/1/1'), 'canine lung cancer': (0.50, 0.073, 0.376, 0.374, 0.419, 0.417, 0.421, 4.0, '4/1/2'), 'canine lymphosarcoma': (0.39, 0.039, 0.605, 0.605, 0.664, 0.674, 0.674, 6.9, '6/1/0'), 'canine soft tissue sarcoma': (0.55, 0.040, 0.576, 0.579, 0.607, 0.610, 0.605, 3.3, '5/1/1'), 'human breast cancer': (0.65, 0.096, 0.631, 0.636, 0.650, 0.652, 0.650, 2.1, '4/2/1'), 'human melanoma': (0.56, 0.016, 0.393, 0.395, 0.402, 0.410, 0.410, 1.7, '4/1/2'), 'human neuroendocrine tumor': (0.59, 0.154, 0.564, 0.581, 0.610, 0.621, 0.600, 5.7, '6/0/1')}
dbad = 0
for r in DT.to_dict('records'):
    v = nb_dom[r['domain']]
    mine = (round(r['clip_share'], 2), round(r['white_share'], 3), round(r['A0'], 3), round(r['A1'], 3), round(r['A2'], 3), round(r['A3'], 3), round(r['A2m'], 3), round(r['A3-A0 pts'], 1))
    dbad += sum(int(abs(a - b) > 1e-9) for a, b in zip(mine, v[:8])) + int(r['WLT'] != v[8])
check('domain', 'cell-34 per-domain table (clip, white, P@10 x5 arms, A3-A0, W/L/T)', 7 * 9, dbad, 'A')
claim('C20', 'od A3-A0 by domain min (melanoma) .. max (neuro), pts', '-0.5..+6.0', f"{DT['od A3-A0 pts'].min():+.1f}..{DT['od A3-A0 pts'].max():+.1f}")
LEAK = pd.read_csv(f'{OUT}_leak_from_pickle.csv')
claim('C22', 'band leak: clicks where A1 (A3) pre-NMS top 100 differs from the spec version', '68 (67)', f"{int((LEAK.A1_peaks_only_floored > 0).sum())} ({int((LEAK.A3_peaks_only_floored > 0).sum())})")
claim('C22', 'band leak: clicks where the A1 (A3) post-NMS pool differs', '68 (67)', f"{int((~LEAK.A1_pool_equal.astype(bool)).sum())} ({int((~LEAK.A3_pool_equal.astype(bool)).sum())})")
claim('C25', 'median seconds per click (all 7 arms, 3 correlations)', 6.9, float(TM.t_click_s.median()), 0.05)
rng_det = TM.groupby('arm').n_detections.agg(['min', 'max'])
claim('C18', 'pool size "83-100 detections per click" (all arms)', '83-100', f"{int(rng_det['min'].min())}-{int(rng_det['max'].max())} (A0 {int(rng_det.loc['A0', 'min'])}-{int(rng_det.loc['A0', 'max'])}, A3 {int(rng_det.loc['A3', 'min'])}-{int(rng_det.loc['A3', 'max'])})")
claim('C20', 'A3-A0 positive in every domain (tm)', 'yes', 'yes' if (DT['A3-A0 pts'] > 0).all() else 'no')
say(pd.DataFrame(P).to_string(index=False, max_colwidth=70))

# ---------------------------------------------------------------- 8. what "changes" means: churn by coordinate vs by object, and where the FP reduction sits
say('\n== 8. churn re-analysis')


def churn10(a, b, axis, tol):
    """
        Top-10 composition change from arm a to arm b, matching slots by Chebyshev distance <= tol px, and mitoses also by annotation id.

        a (str): from-arm.
        b (str): to-arm.
        axis (str): rank_key.
        tol (int): coordinate tolerance in px (0 = exact equality, the notebook's apparent definition).

        Returns dict: counts.
    """
    t10 = top[(top['rank'] < 10) & (top.rank_key == axis)]
    A, B = t10[t10.arm == a], t10[t10.arm == b]
    gb = dict(list(B.groupby(['file_name', 'seed_index'])))
    c = dict(arm_from=a, arm_to=b, axis=axis, tol_px=tol, slots=0, slots_changed=0, mit_leave=0, mit_enter=0, fp_leave=0, fp_enter=0, white_fp_leave=0, white_fp_enter=0, mit_leave_by_ann_id=0, mit_enter_by_ann_id=0)
    for k, ga in A.groupby(['file_name', 'seed_index']):
        gk = gb[k]
        D = np.max(np.abs(ga[['cx', 'cy']].to_numpy()[:, None, :] - gk[['cx', 'cy']].to_numpy()[None, :, :]), axis=2)
        ain, bin_ = (D <= tol).any(1), (D <= tol).any(0)
        ta, tb = (ga.bucket == TP).to_numpy(), (gk.bucket == TP).to_numpy()
        wa, wb_ = (ga.white_frac_61 > WHITE_SHARE).to_numpy(), (gk.white_frac_61 > WHITE_SHARE).to_numpy()
        c['slots'] += len(ga)
        c['slots_changed'] += int((~ain).sum())
        c['mit_leave'] += int((~ain & ta).sum())
        c['mit_enter'] += int((~bin_ & tb).sum())
        c['fp_leave'] += int((~ain & ~ta).sum())
        c['fp_enter'] += int((~bin_ & ~tb).sum())
        c['white_fp_leave'] += int((~ain & ~ta & wa).sum())
        c['white_fp_enter'] += int((~bin_ & ~tb & wb_).sum())
        ia, ib = set(ga.matched_ann_id[ta]), set(gk.matched_ann_id[tb])
        c['mit_leave_by_ann_id'] += len(ia - ib)
        c['mit_enter_by_ann_id'] += len(ib - ia)
    return c


CH10 = pd.DataFrame([churn10('A0', b, 'tm_score', tol) for b in ('A3', 'A2') for tol in (0, 1, 2, 3, 5)])
say(CH10.to_string(index=False))
CH10.to_csv(f'{OUT}_churn_top10.csv', index=False)
c0 = CH10[(CH10.arm_to == 'A3') & (CH10.tol_px == 0)].iloc[0]
claim('C9', 'top-10 slots that change A0->A3 (exact coordinates)', '1193 of 2940', f"{c0.slots_changed} of {c0.slots}")
claim('C9', 'mitoses leave / enter (exact coordinates)', '560/669', f'{c0.mit_leave}/{c0.mit_enter}')


def churn_pool(a, b, tol):
    """
        Whole-pool composition change, by coordinate tolerance and (mitoses) by annotation id.

        a (str): from-arm.
        b (str): to-arm.
        tol (int): Chebyshev tolerance px.

        Returns dict: counts.
    """
    A, B = PM[PM.arm == a], PM[PM.arm == b]
    gb = dict(list(B.groupby(['file_name', 'seed_index'])))
    c = dict(arm_from=a, arm_to=b, tol_px=tol, mit_leave=0, mit_enter=0, nonmit_leave=0, nonmit_enter=0, mit_leave_by_ann_id=0, mit_enter_by_ann_id=0)
    for k, ga in A.groupby(['file_name', 'seed_index']):
        gk = gb[k]
        D = np.max(np.abs(ga[['cx', 'cy']].to_numpy()[:, None, :] - gk[['cx', 'cy']].to_numpy()[None, :, :]), axis=2)
        ain, bin_ = (D <= tol).any(1), (D <= tol).any(0)
        ma, mb = ga.cat.to_numpy() == 1, gk.cat.to_numpy() == 1
        c['mit_leave'] += int((~ain & ma).sum())
        c['mit_enter'] += int((~bin_ & mb).sum())
        c['nonmit_leave'] += int((~ain & ~ma).sum())
        c['nonmit_enter'] += int((~bin_ & ~mb).sum())
        ia, ib = set(ga.mid[ma]), set(gk.mid[mb])
        c['mit_leave_by_ann_id'] += len(ia - ib)
        c['mit_enter_by_ann_id'] += len(ib - ia)
    return c


CHP = pd.DataFrame([churn_pool('A0', 'A2', tol) for tol in (0, 1, 3)] + [churn_pool('A0', 'A3', 1)])
say(CHP.to_string(index=False))
CHP.to_csv(f'{OUT}_churn_pool.csv', index=False)
p0 = CHP[(CHP.arm_to == 'A2') & (CHP.tol_px == 0)].iloc[0]
for nm, v in (('mit_leave', 2520), ('nonmit_leave', 10310), ('mit_enter', 2960), ('nonmit_enter', 9790)):
    claim('C17', f'pool A0->A2 {nm} ("about", exact-coordinate definition; prose rounds to 10)', v, int(p0[nm]), 5)
od10 = top[(top.rank_key == 'chromatin_od') & (top['rank'] < 10)]
fl_rows = []
for tol in (0, 1, 3):
    r = dict(tol_px=tol, enter_mit=0, enter_nonmit=0, leave_mit=0, leave_nonmit=0, enter_mit_new_ann_id=0, leave_mit_gone_ann_id=0)
    PA, PB = dict(list(PM[PM.arm == 'A0'].groupby(['file_name', 'seed_index']))), dict(list(PM[PM.arm == 'A2'].groupby(['file_name', 'seed_index'])))
    TA, TB = dict(list(od10[od10.arm == 'A0'].groupby(['file_name', 'seed_index']))), dict(list(od10[od10.arm == 'A2'].groupby(['file_name', 'seed_index'])))
    for k in PA:
        pa, pb_, ta, tb = PA[k], PB[k], TA[k], TB[k]
        new = ~(np.max(np.abs(tb[['cx', 'cy']].to_numpy()[:, None, :] - pa[['cx', 'cy']].to_numpy()[None, :, :]), axis=2) <= tol).any(1)
        gone = ~(np.max(np.abs(ta[['cx', 'cy']].to_numpy()[:, None, :] - pb_[['cx', 'cy']].to_numpy()[None, :, :]), axis=2) <= tol).any(1)
        mb, ma = (tb.bucket == TP).to_numpy(), (ta.bucket == TP).to_numpy()
        r['enter_mit'] += int((new & mb).sum())
        r['enter_nonmit'] += int((new & ~mb).sum())
        r['leave_mit'] += int((gone & ma).sum())
        r['leave_nonmit'] += int((gone & ~ma).sum())
        ida, idb = set(pa.mid[pa.cat == 1]), set(pb_.mid[pb_.cat == 1])
        r['enter_mit_new_ann_id'] += int(sum(1 for i in tb.matched_ann_id[mb] if i not in ida))
        r['leave_mit_gone_ann_id'] += int(sum(1 for i in ta.matched_ann_id[ma] if i not in idb))
    fl_rows.append(r)
ODF = pd.DataFrame(fl_rows)
say('od top-10 entries that are new to the pool (A0 -> A2):')
say(ODF.to_string(index=False))
ODF.to_csv(f'{OUT}_od_top10_flows.csv', index=False)
o0 = ODF[ODF.tol_px == 0].iloc[0]
for nm, v in (('enter_mit', 631), ('enter_nonmit', 571), ('leave_mit', 604), ('leave_nonmit', 507)):
    claim('C17', f'od top 10, A0->A2 {nm} (exact-coordinate definition; +/-1 tolerated)', v, int(o0[nm]), 1)
# FP reduction by white adjacency, per domain, and covariates at ROI level
rows = []
for dname in list(pd.unique(DOM.reindex(ORDER))) + ['ALL']:
    fs = ORDER if dname == 'ALL' else [f for f in ORDER if DOM[f] == dname]
    s = lambda a, c: int(TM[(TM.arm == a) & TM.file_name.isin(fs)][c].sum())
    fp0, w0, fp3, w3 = s('A0', 'n_fp_top10'), s('A0', 'n_white_fp_top10'), s('A3', 'n_fp_top10'), s('A3', 'n_white_fp_top10')
    rows.append(dict(domain=dname, A0_fp=fp0, A0_white_fp=w0, A0_white_share=w0 / fp0, A3_fp=fp3, A3_white_fp=w3, net_fp=fp3 - fp0, net_white_fp=w3 - w0, net_nonwhite_fp=(fp3 - w3) - (fp0 - w0), white_fp_removed_frac=(w0 - w3) / w0 if w0 else np.nan, nonwhite_fp_removed_frac=((fp0 - w0) - (fp3 - w3)) / (fp0 - w0)))
WD = pd.DataFrame(rows)
say(WD.to_string(index=False, float_format=lambda v: f'{v:.3f}'))
WD.to_csv(f'{OUT}_white_fp_decomposition.csv', index=False)
cov = []
w0r, f0r = rsum('tm_score', 'A0', 'n_white_fp_top10'), rsum('tm_score', 'A0', 'n_fp_top10')
for cname, cv_ in (('ROI pixel white share', fl.loc[ORDER, 'white_share_roi'].to_numpy()), ('ROI clip share', fl.loc[ORDER, 'clip_share_otsu'].to_numpy()), ('A0 white-adjacent top-10 FP count', w0r), ('A0 white-adjacent share of top-10 FPs', w0r / np.where(f0r == 0, np.nan, f0r)), ('A0 top-10 FP count (headroom)', f0r)):
    for cn in ('A2 - A0', 'A3 - A0'):
        a, b = CON[cn]
        d = con('tm_score', a, b, 10)
        ok = ~np.isnan(cv_.astype(float))
        sr = sps.spearmanr(cv_[ok], d[ok])
        cov.append(dict(covariate=cname, contrast=cn, rho=sr.statistic, p=sr.pvalue, n=int(ok.sum())))
COV = pd.DataFrame(cov)
say(COV.to_string(index=False, float_format=lambda v: f'{v:+.3f}'))
COV.to_csv(f'{OUT}_white_covariates.csv', index=False)
# per click, the product view
pc = TM.set_index(['file_name', 'seed_index', 'arm'])['tp_at_10'].unstack('arm')
dpc = pc['A3'] - pc['A0']
say(f'per click A3 - A0 TP@10: {int((dpc > 0).sum())} better / {int((dpc < 0).sum())} worse / {int((dpc == 0).sum())} same; worst loss {int(dpc.min())}; 245.tiff s1 A0 {int(pc.loc[("245.tiff", 1), "A0"])} -> A3 {int(pc.loc[("245.tiff", 1), "A3"])}')
pd.DataFrame(dict(A0=pc['A0'], A3=pc['A3'], delta=dpc)).reset_index().to_csv(f'{OUT}_per_click_tp10.csv', index=False)

# ---------------------------------------------------------------- 9. Tier B: pixels
say('\n== 9. Tier B: floors on all 49 ROIs, an independent search on the Tier B ROIs, white fraction, M5')
FLOORX = np.float32(-3.0e38)
clicks = pd.read_csv('results/precision_at_k_49roi_6seed_fillgeom_postD11_clicks.csv')
tpl = pd.read_csv('results/precision_at_k_49roi_6seed_mopad_postD11_templates.csv').query("condition == 'mo_pad61'").set_index(['file_name', 'seed_index'])


def win(c_, cx, cy, S):
    """
        S x S window centred on integer (cx, cy), from a 36 px replicate-padded copy.

        c_ (np.ndarray): channel.
        cx (int): centre x.
        cy (int): centre y.
        S (int): odd side.

        Returns np.ndarray: float32 copy.
    """
    h = S // 2
    p = cv2.copyMakeBorder(c_, 36, 36, 36, 36, cv2.BORDER_REPLICATE)
    return p[cy + 36 - h:cy + 36 - h + S, cx + 36 - h:cx + 36 - h + S].astype(np.float32).copy()


def template61(c_, tcx, tcy, bs):
    """
        The mo_pad61 template: a 61 px window with everything outside the centred bs box set to the outside mean.

        c_ (np.ndarray): channel.
        tcx (int): centre x.
        tcy (int): centre y.
        bs (int): box side.

        Returns tuple[np.ndarray, float, float]: template, fill value, mean of the kept box.
    """
    Wn = win(c_, tcx, tcy, 61)
    keep = np.zeros((61, 61), bool)
    o = (61 - bs) // 2
    keep[o:o + bs, o:o + bs] = True
    f = np.float32(Wn[~keep].astype(np.float64).mean())
    T = Wn.copy()
    T[~keep] = f
    return np.ascontiguousarray(T), float(f), float(Wn[keep].astype(np.float64).mean())


def corr(c_, T, cx, cy, side=51):
    """
        TM_CCOEFF response of T over the blanked, 30 px replicate-padded channel, in ROI coordinates.

        c_ (np.ndarray): channel.
        T (np.ndarray): 61 x 61 template.
        cx (float): blank centre x (the click in the notebook; the template centre under D11's text).
        cy (float): blank centre y.
        side (int): odd blank side (51 in the notebook; 61 = the mo_pad61 footprint under D11's text).

        Returns np.ndarray: (h, w) float32.
    """
    s = c_.copy()
    ix, iy, hb = int(round(cx)), int(round(cy)), side // 2
    s[iy - hb:iy + hb + 1, ix - hb:ix + hb + 1] = c_.min()
    res = cv2.matchTemplate(cv2.copyMakeBorder(s, PAD, PAD, PAD, PAD, cv2.BORDER_REPLICATE), T, cv2.TM_CCOEFF).astype(np.float32)
    assert res.shape == c_.shape
    return res


def peaks_nms(res, band, mpp):
    """
        Robust threshold, 15 x 15 local maxima, top 100 by (score desc, x, y), greedy distance NMS at 7.5 um.

        res (np.ndarray): response map.
        band (np.ndarray or None): border band mask (floored and invalid) or None.
        mpp (float): microns per pixel.

        Returns tuple[pd.DataFrame, float]: detections ranked by score, and the threshold.
    """
    valid = np.ones(res.shape, bool)
    if band is not None:
        valid = ~band
        res = np.where(band, FLOORX, res)
    smp = res[::8, ::8][valid[::8, ::8]]
    med = np.median(smp)
    thr = med - 1.5 * 1.4826 * np.median(np.abs(smp - med))
    mx = cv2.dilate(res, np.ones((15, 15), np.uint8))
    ys, xs = np.nonzero((res >= mx) & valid & (res >= thr))
    sc = res[ys, xs]
    o = np.lexsort((ys, xs, -sc))[:100]
    xs, ys, sc = xs[o].astype(float), ys[o].astype(float), sc[o]
    rad = 7.5 / mpp
    keep, sup = [], np.zeros(len(xs), bool)
    for i in np.argsort(-sc.astype(np.float64), kind='stable'):
        if not sup[i]:
            keep.append(i)
            sup |= np.hypot(xs - xs[i], ys - ys[i]) <= rad
    k = np.array(keep)
    det = pd.DataFrame(dict(cx=xs[k], cy=ys[k], score=sc[k]))
    return det.sort_values('score', ascending=False, kind='mergesort').reset_index(drop=True), float(thr)


def od_top10pct(hem_p, cx, cy):
    """
        Mean of the darkest 10 % of a 51 px window of the unclipped channel (25 px replicate pad).

        hem_p (np.ndarray): channel padded by 25.
        cx (np.ndarray): x.
        cy (np.ndarray): y.

        Returns np.ndarray: od per detection.
    """
    out = []
    for x, y in zip(cx.astype(int), cy.astype(int)):
        flat = hem_p[y:y + 51, x:x + 51].ravel()
        out.append(float(np.partition(flat, -260)[-260:].mean()))
    return np.array(out)


floor_rows, tb_rows, wf_rows, timing_rows, d11_rows = [], [], [], [], []
D11_ROIS = ('300.tiff', '201.tiff', '291.tiff')   # the D11-footprint blank (61 px at the template centre) instead of the notebook's 51 px at the click
t_tb = time.time()
FLOOR_ROIS = TIERB_ROIS if os.environ.get('AUDIT_QUICK') else ORDER   # AUDIT_QUICK only for debugging; the reported run is the full one
for f in FLOOR_ROIS:
    c0_ = clicks[clicks.file_name == f].iloc[0]
    path = f"{SUB[c0_.subset]}/{f}"
    rgb = ds.load_roi(path)
    mpp = ds.roi_mpp(path)
    hem = ch.to_channel(rgb, 'hematoxylin_od')
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    b_o = np.float32(threshold_multiotsu(hem, classes=3)[0])
    b_m = np.float32(np.median(hem))
    floor_rows.append(dict(file_name=f, hem_min=float(hem.min()), zero_share=float((hem == 0).mean()), b_otsu=float(b_o), b_med=float(b_m), clip_share_otsu=float((hem < b_o).mean()), clip_share_med=float((hem < b_m).mean()), white_share_roi=float((gray > WHITE_GRAY).mean())))
    if f in TIERB_ROIS:
        h, w = hem.shape
        yy, xx = np.arange(h), np.arange(w)
        band = np.minimum(np.minimum(yy, h - 1 - yy)[:, None], np.minimum(xx, w - 1 - xx)[None, :]) < BAND
        hem_p = cv2.copyMakeBorder(hem, 25, 25, 25, 25, cv2.BORDER_REPLICATE)
        gpad = cv2.copyMakeBorder(gray, 30, 30, 30, 30, cv2.BORDER_REPLICATE)
        chans = {None: hem, 'b_otsu': np.maximum(hem, b_o), 'b_med': np.maximum(hem, b_m)}
        for r in clicks[clicks.file_name == f].sort_values('seed_index').itertuples():
            t_ = tpl.loc[(f, r.seed_index)]
            gid, gxy, gcat = gt_of(f, r.seed_ann_id)
            for fk, arms in ((None, (('A0', False), ('A1', True))), ('b_otsu', (('A2', False), ('A3', True))), ('b_med', (('A2m', False),))):
                T, fill, mk = template61(chans[fk], int(t_.tcx), int(t_.tcy), int(t_.bs))
                sha = hashlib.sha1(np.ascontiguousarray(T, dtype=np.float32).tobytes()).hexdigest()
                tc0 = time.perf_counter()
                res = corr(chans[fk], T, float(r.click_cx), float(r.click_cy))
                t_corr = time.perf_counter() - tc0
                for arm, use_band in arms:
                    tp0 = time.perf_counter()
                    det, thr = peaks_nms(res, band if use_band else None, mpp)
                    timing_rows.append(dict(file_name=f, seed_index=r.seed_index, arm=arm, t_corr_s=t_corr, t_peaks_nms_s=time.perf_counter() - tp0))
                    if f in D11_ROIS and arm in ('A0', 'A3'):
                        det61, _ = peaks_nms(corr(chans[fk], T, int(t_.tcx), int(t_.tcy), side=61), band if use_band else None, mpp)
                        cat61, _ = greedy(det61[['cx', 'cy']].to_numpy(), gid, gxy, gcat, 7.5 / mpp)
                        ref = PRK.loc[(f, r.seed_index, arm, 'tm_score')]
                        tt = top[(top.file_name == f) & (top.seed_index == r.seed_index) & (top.arm == arm) & (top.rank_key == 'tm_score')].sort_values('rank')
                        d11_rows.append(dict(file_name=f, seed_index=r.seed_index, arm=arm, top30_xy_equal=bool(np.array_equal(det61.head(30)[['cx', 'cy']].to_numpy(), tt[['cx', 'cy']].to_numpy())), **{f'tp_at_{k}_d11': int((cat61[:k] == 1).sum()) for k in (10, 20, 30)}, **{f'tp_at_{k}_nb': int(ref[f'tp_at_{k}']) for k in (10, 20, 30)}))
                    det['od'] = od_top10pct(hem_p, det.cx.to_numpy(), det.cy.to_numpy())
                    for axis, keycol in (('tm_score', 'score'), ('chromatin_od', 'od')):
                        rk = det.sort_values(keycol, ascending=False, kind='mergesort').reset_index(drop=True)
                        cat, _ = greedy(rk[['cx', 'cy']].to_numpy(), gid, gxy, gcat, 7.5 / mpp)
                        ref = PRK.loc[(f, r.seed_index, arm, axis)]
                        tt = top[(top.file_name == f) & (top.seed_index == r.seed_index) & (top.arm == arm) & (top.rank_key == axis)].sort_values('rank')
                        tb_rows.append(dict(file_name=f, seed_index=r.seed_index, arm=arm, rank_key=axis, tpl_sha1_equal=sha == ref.tpl_sha1, fill_equal=bool(np.isclose(fill, ref.fill_value, rtol=1e-6)), threshold_equal=bool(np.isclose(thr, ref.threshold, rtol=1e-9)), n_det_equal=len(det) == ref.n_detections, top30_xy_equal=bool(np.array_equal(rk.head(30)[['cx', 'cy']].to_numpy(), tt[['cx', 'cy']].to_numpy())), top30_score_equal=bool(np.array_equal(rk.head(30).score.to_numpy(np.float32), tt.score.to_numpy(np.float32))), top30_od_close=bool(np.allclose(rk.head(30).od.to_numpy(), tt.od.to_numpy(), rtol=1e-6)), tp_equal=all(int((cat[:k] == 1).sum()) == int(ref[f'tp_at_{k}']) for k in (10, 20, 30))))
                        if axis == 'tm_score' and arm in ('A0', 'A3'):
                            xs_, ys_ = tt.cx.to_numpy().astype(int), tt.cy.to_numpy().astype(int)
                            wf = np.array([(gpad[y:y + 61, x:x + 61] > WHITE_GRAY).mean() for x, y in zip(xs_, ys_)])
                            wf_rows.append(dict(file_name=f, seed_index=r.seed_index, arm=arm, n=len(wf), n_equal=int(np.isclose(wf, tt.white_frac_61.to_numpy(), atol=1e-12).sum())))
        say(f'  Tier B search {f}: {time.time() - t_tb:.0f} s')
    del rgb, hem, gray
FLB = pd.DataFrame(floor_rows).set_index('file_name')
bad = 0
for c in ('hem_min', 'b_otsu', 'b_med', 'clip_share_otsu', 'clip_share_med', 'white_share_roi'):
    bad += int((~np.isclose(FLB[c].to_numpy(), fl.loc[FLB.index, c].to_numpy(), rtol=1e-6, atol=1e-9)).sum())
check('tierB', f'floors.csv (hem min, b_otsu, b_med, clip shares, white share) re-derived from pixels, {len(FLB)} ROIs', len(FLB) * 6, bad, 'B', '(to_channel is a dependency)')
say(f'share of hematoxylin_od pixels that are exactly 0 (the channel already has a bottom clip at 0): median {FLB.zero_share.median():.3f}, range {FLB.zero_share.min():.3f}-{FLB.zero_share.max():.3f}')
FLB.reset_index().to_csv(f'{OUT}_tierb_floors.csv', index=False)
TB = pd.DataFrame(tb_rows)
TB.to_csv(f'{OUT}_tierb_search.csv', index=False)
flags = ['tpl_sha1_equal', 'fill_equal', 'threshold_equal', 'n_det_equal', 'top30_xy_equal', 'top30_score_equal', 'top30_od_close', 'tp_equal']
check('tierB', f'independent search re-implementation, {TB.file_name.nunique()} ROIs x 6 clicks x 5 arms x 2 axes: template SHA, fill, threshold, n_det, top-30 xy/score/od, TP@10/20/30', TB[flags].size, int((~TB[flags]).sum().sum()), 'B')
say((~TB[flags]).groupby(TB.arm).sum().to_string())
TIM = pd.DataFrame(timing_rows)
TIM.to_csv(f'{OUT}_tierb_timing.csv', index=False)
say(f"re-timing (this machine, load as found): correlation median {TIM.drop_duplicates(['file_name', 'seed_index', 't_corr_s']).t_corr_s.median():.2f} s, peaks+NMS median {TIM.t_peaks_nms_s.median():.2f} s per arm")
D11B = pd.DataFrame(d11_rows)
if len(D11B):
    D11B.to_csv(f'{OUT}_tierb_d11_footprint_blank.csv', index=False)
    tpd = sum(int((D11B[f'tp_at_{k}_d11'] != D11B[f'tp_at_{k}_nb']).sum()) for k in (10, 20, 30))
    check('tierB', f'D11-footprint blank (61 px at the template centre) vs the notebook 51 px at the click, A0 and A3, {D11B.file_name.nunique()} ROIs: TP@10/20/30 changes', len(D11B) * 3, tpd, 'B', f'(top-30 identical on {int(D11B.top30_xy_equal.sum())}/{len(D11B)} lists)')
WF = pd.DataFrame(wf_rows)
check('tierB', 'white_frac_61 recomputed from RGB (top 30, A0 and A3, tm_score)', int(WF.n.sum()), int((WF.n - WF.n_equal).sum()), 'B')
# M5, re-derived with this file's own template and correlation
f5 = '300.tiff'
c5 = clicks[(clicks.file_name == f5) & (clicks.seed_index == 0)].iloc[0]
rgb5 = ds.load_roi(f"{SUB[c5.subset]}/{f5}")
hem5 = ch.to_channel(rgb5, 'hematoxylin_od')
gray5 = cv2.cvtColor(rgb5, cv2.COLOR_RGB2GRAY)
b5 = np.float32(threshold_multiotsu(hem5, classes=3)[0])
t5 = tpl.loc[(f5, 0)]
T0, fill0, mk0 = template61(hem5, int(t5.tcx), int(t5.tcy), int(t5.bs))
Tc, fillc, _ = template61(np.maximum(hem5, b5), int(t5.tcx), int(t5.tcy), int(t5.bs))
t_tis = float(np.median(hem5[gray5 <= WHITE_GRAY]))
glass = float(np.median(hem5[gray5 > WHITE_GRAY]))
t_ab = float(np.median(hem5[(gray5 <= WHITE_GRAY) & (hem5 >= b5)]))
dr = (int(t5.bs) - 4) / 2.0
del rgb5, hem5, gray5
SIDE, C = 301, 150


def disc(t, xy):
    """
        Flat synthetic channel at level t with one flat disc at the template's own box mean.

        t (float): background level.
        xy (tuple): disc centre.

        Returns np.ndarray: float32 image.
    """
    im = np.full((SIDE, SIDE), t, np.float32)
    yy, xx = np.mgrid[:SIDE, :SIDE]
    im[(xx - xy[0]) ** 2 + (yy - xy[1]) ** 2 <= dr ** 2] = mk0
    return im


def resp(im, T, const=None):
    """
        TM_CCOEFF response over a 30 px padded synthetic image (replicate, or constant at a level).

        im (np.ndarray): image.
        T (np.ndarray): template.
        const (float or None): constant pad level; None = replicate.

        Returns np.ndarray: (SIDE, SIDE) response.
    """
    p = cv2.copyMakeBorder(im, PAD, PAD, PAD, PAD, cv2.BORDER_REPLICATE) if const is None else cv2.copyMakeBorder(im, PAD, PAD, PAD, PAD, cv2.BORDER_CONSTANT, value=float(const))
    return cv2.matchTemplate(p, T, cv2.TM_CCOEFF)


def pk(res, x, y, r=15):
    """
        Max response within r px (Chebyshev) of (x, y), clipped to the image.

        res (np.ndarray): response.
        x (int): x.
        y (int): y.
        r (int): half-width.

        Returns tuple[float, int]: value and its y.
    """
    y0, x0 = max(0, y - r), max(0, x - r)
    sub = res[y0:y + r + 1, x0:x + r + 1]
    iy, ix = np.unravel_index(int(np.argmax(sub)), sub.shape)
    return float(sub[iy, ix]), int(y0 + iy)


m5 = {}
for lvl, t in (('t_tissue', t_tis), ('t_above_b', t_ab)):
    for gap in (30, 20):
        base = disc(t, (C, C))
        wim = base.copy()
        wim[:, C + gap:] = glass
        pu0, pw0 = pk(resp(base, T0), C, C)[0], pk(resp(wim, T0), C, C)[0]
        puc, pwc = pk(resp(np.maximum(base, b5), Tc), C, C)[0], pk(resp(np.maximum(wim, b5), Tc), C, C)[0]
        m5[f'a {lvl} gap{gap} A0 ratio'] = pw0 / pu0
        m5[f'a {lvl} gap{gap} clip ratio'] = pwc / puc
        m5[f'a {lvl} gap{gap} excess frac'] = (pwc - puc) / (pw0 - pu0)
m5['a predicted frac t_above_b'] = (t_ab - float(b5)) / (t_ab - glass)
pe, ye = pk(resp(disc(t_tis, (C, 0)), T0), C, 0)
pi = pk(resp(disc(t_tis, (C, 60)), T0), C, 60)[0]
m5['b edge/interior'] = pe / pi
m5['b edge peak y'] = ye
for yc in (0, -6):
    im = disc(t_tis, (C, yc))
    m5[f'b replicate/constant y{yc}'] = pk(resp(im, T0), C, 0)[0] / pk(resp(im, T0, const=t_tis), C, 0)[0]
dx = C + 25 + 20
for nm, fv in (('none', None), ('min', 0.0), ('b', float(b5))):
    im = disc(t_tis, (dx, C))
    if fv is not None:
        im[C - 25:C + 26, C - 25:C + 26] = fv
    m5[f'c square={nm}'] = pk(resp(im, T0), dx, C)[0]
nb_m5 = {'a t_tissue gap30 A0 ratio': 1.0384, 'a t_tissue gap30 clip ratio': 1.0, 'a t_tissue gap20 A0 ratio': 1.1493, 'a t_above_b gap30 excess frac': 0.27026, 'a t_above_b gap20 excess frac': 0.36062, 'a predicted frac t_above_b': 0.48026, 'b edge/interior': 0.87473, 'b edge peak y': 0, 'b replicate/constant y0': 1.4131, 'b replicate/constant y-6': 2.0414, 'c square=none': 1.1811, 'c square=min': 1.328, 'c square=b': 1.1726}
M5 = pd.DataFrame([dict(item=k, notebook=v, recomputed=m5[k], match=bool(abs(m5[k] - v) <= 5e-4 * max(1, abs(v)))) for k, v in nb_m5.items()])
say(M5.to_string(index=False))
M5.to_csv(f'{OUT}_m5.csv', index=False)
check('tierB', 'M5 synthetic demonstrations (a, b, b2/b3, c) re-derived with an independent template and correlation', len(M5), int((~M5.match).sum()), 'B')
say(f'M5 levels: t_tissue {t_tis:.4f}, b_otsu {float(b5):.4f}, glass {glass:.4f}, t_above_b {t_ab:.4f}, A0 fill {fill0:.4f}, clipped fill {fillc:.4f}')

# ---------------------------------------------------------------- summary
PROSE = pd.DataFrame(P)
PROSE.to_csv(f'{OUT}_prose_claims.csv', index=False)
check('prose', 'numbers stated in cell 34 / cell outputs (incl. the churn and flow counts), recomputed', len(PROSE), int((~PROSE.match).sum()), 'A', '(see prose_claims.csv for which)')
CHK = pd.DataFrame(CHECKS)
CHK.to_csv(f'{OUT}_checks.csv', index=False)
say('\n== summary of checks')
say(CHK.to_string(index=False, max_colwidth=110))
say('\nprose claims that do not match:')
say(PROSE[~PROSE.match].to_string(index=False, max_colwidth=90))
say(f'\ntotal values compared {int(CHK.n_compared.sum())}, divergent {int(CHK.n_divergent.sum())}; runtime {(time.time() - T_START) / 60:.1f} min')
