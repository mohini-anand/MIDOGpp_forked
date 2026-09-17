"""
    Independent audit of production_hematoxylin_only/bbox_refinement_three_way_chromatin_od.ipynb.

    Tier A re-derives every TP count, table, statistic and prose number from the notebook's own CSVs and
    from databases/MIDOG++.json, using a matcher written here. Tier B re-implements the seed draw and the
    production search (template cut, TM_CCOEFF, deep floor, MAX_PEAKS, NMS, self-hit, od51 ranking) from
    cv2 / skimage / numpy primitives on all 42 runs, then scores the full lists.

    Nothing from midog_utils is imported. Run with /Users/mohinianand/anaconda3/bin/python3.
    Writes results/bbox_refinement_three_way_chromatin_od_audit_*.csv only.
"""

from __future__ import annotations

import ast
import itertools
import json
import os
import re
import subprocess
import sys
import time
from decimal import ROUND_HALF_UP, Decimal
from fractions import Fraction

import cv2
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from skimage.color import rgb2hed
from skimage.measure import label, regionprops

REPO = '/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked'
NB_PATH = f'{REPO}/production_hematoxylin_only/bbox_refinement_three_way_chromatin_od.ipynb'
REF_NB_PATH = f'{REPO}/production_hematoxylin_only/production_seed_precision_at_k_chromatin_hem_bbox.ipynb'
DB_PATH = f'{REPO}/databases/MIDOG++.json'
IMG_DIR = f'{REPO}/images/extra_valid'
STEM = f'{REPO}/results/precision_at_k_14roi_prodseed_chromatin_bbox3way'
REF_RAW = f'{REPO}/results/precision_at_k_14roi_prodseed_chromatin_hembbox_raw.csv'
OUT = f'{REPO}/results/bbox_refinement_three_way_chromatin_od_audit'

CONDITIONS = ('default_51', 'gray_bbox', 'hem_bbox')
PAIRS = (('gray_bbox', 'default_51'), ('hem_bbox', 'default_51'), ('hem_bbox', 'gray_bbox'))
BUDGETS = (10, 20, 30)
RADIUS_UM = 7.5
MITOTIC, LOOKALIKE = 1, 2
TUMOR_ALIASES = {'canine lymphoma': 'canine lymphosarcoma'}
PATCH, OTSU_WIN, BASE = 73, 51, 51
MAX_PEAKS, PEAK_MIN_DIST, DEEP_Z, SELF_HIT = 100, 7, -1.5, 5.0
OD_WIN, OD_FRAC = 51, 0.10

COMPARISONS = []  # (table, key, column, notebook_value, audit_value, match)


def log(msg):
    """
        Print with flush.

        msg (str): message.

        Returns None.
    """
    print(msg, flush=True)


def compare(table, key, column, nb_value, audit_value, tol=1e-9):
    """
        Record one notebook-vs-audit value comparison.

        table (str): which notebook table or claim the value comes from.
        key (str): row identifier.
        column (str): column identifier.
        nb_value (object): the notebook's value.
        audit_value (object): the independently recomputed value.
        tol (float): absolute tolerance for numeric values.

        Returns bool: whether the two agree.
    """
    try:
        ok = bool(abs(float(nb_value) - float(audit_value)) <= tol)
    except (TypeError, ValueError):
        ok = str(nb_value) == str(audit_value)
    COMPARISONS.append(dict(table=table, key=key, column=column, notebook=nb_value, audit=audit_value, match=ok))
    return ok


def save(df, name):
    """
        Write one audit table.

        df (pd.DataFrame): the table.
        name (str): suffix after the audit stem.

        Returns None.
    """
    path = f'{OUT}_{name}.csv'
    df.to_csv(path, index=False)
    log(f'  -> {os.path.relpath(path, REPO)} ({len(df)} rows)')


# --------------------------------------------------------------------------------------------------
# primary sources, read without midog_utils
# --------------------------------------------------------------------------------------------------

def load_db(path):
    """
        Parse MIDOG++.json directly.

        path (str): JSON path.

        Returns tuple[pd.DataFrame, pd.DataFrame]: images (id, file, size, tumour) and annotations (id, image, centre, category, votes).
    """
    raw = json.load(open(path))
    images = pd.DataFrame([dict(image_id=im['id'], file_name=im['file_name'], width=im['width'], height=im['height'], tumor_type=TUMOR_ALIASES.get(im['tumor_type'], im['tumor_type'])) for im in raw['images']])
    rows = []
    for a in raw['annotations']:
        x1, y1, x2, y2 = a['bbox']
        votes = list(a.get('labels', []))
        rows.append(dict(ann_id=a['id'], image_id=a['image_id'], cx=(x1 + x2) / 2.0, cy=(y1 + y2) / 2.0, w=x2 - x1, h=y2 - y1, category_id=a['category_id'], n_votes=len(votes), n_mitotic_votes=sum(1 for v in votes if v == MITOTIC)))
    return images, pd.DataFrame(rows)


def tiff_meta(path):
    """
        Microns per pixel and array shape from the TIFF tags, read directly.

        path (str): ROI TIFF path.

        Returns dict: mpp_x, mpp_y, height, width.
    """
    with tifffile.TiffFile(path) as tf:
        page = tf.pages[0]
        xn, xd = page.tags['XResolution'].value
        yn, yd = page.tags['YResolution'].value
        unit = int(page.tags['ResolutionUnit'].value)
        series = tf.series[0]
        levels = getattr(series, 'levels', None)
        shape = (levels[0].shape if levels else series.shape)
    um = {2: 25400.0, 3: 10000.0}[unit]
    return dict(mpp_x=um / (xn / xd), mpp_y=um / (yn / yd), height=int(shape[0]), width=int(shape[1]), unit=unit)


def load_rgb(path):
    """
        Full-resolution RGB array.

        path (str): ROI TIFF path.

        Returns np.ndarray: HxWx3 uint8.
    """
    with tifffile.TiffFile(path) as tf:
        series = tf.series[0]
        levels = getattr(series, 'levels', None)
        arr = levels[0].asarray() if levels else series.asarray()
    return np.ascontiguousarray(arr[:, :, :3])


# --------------------------------------------------------------------------------------------------
# matching, written here
# --------------------------------------------------------------------------------------------------

def greedy(det_xy, gt_xy, radius, strict=False):
    """
        Rank-order greedy one-to-one match: each detection claims its nearest unclaimed GT within radius.

        det_xy (np.ndarray): (N, 2) detections, best first.
        gt_xy (np.ndarray): (M, 2) ground-truth centres.
        radius (float): match radius in px.
        strict (bool): use d < radius instead of d <= radius.

        Returns tuple[np.ndarray, int]: GT index per detection (-1 unmatched), and how many claims had an exact-distance tie.
    """
    out = np.full(len(det_xy), -1, dtype=int)
    if len(det_xy) == 0 or len(gt_xy) == 0:
        return out, 0
    taken = np.zeros(len(gt_xy), dtype=bool)
    d_all = cdist(np.asarray(det_xy, float), np.asarray(gt_xy, float))
    ties = 0
    for i in range(len(det_xy)):
        d = d_all[i]
        within = (d < radius) if strict else (d <= radius)
        cand = np.where(within & ~taken)[0]
        if len(cand) == 0:
            continue
        best = cand[d[cand] == d[cand].min()]
        ties += int(len(best) > 1)
        out[i] = best[0]
        taken[best[0]] = True
    return out, ties


def max_matching(det_xy, gt_xy, radius):
    """
        Maximum-cardinality one-to-one matching within radius (the most TPs any assignment can give).

        det_xy (np.ndarray): (N, 2) detections.
        gt_xy (np.ndarray): (M, 2) ground truth.
        radius (float): match radius in px.

        Returns int: matched pairs.
    """
    if len(det_xy) == 0 or len(gt_xy) == 0:
        return 0
    within = cdist(np.asarray(det_xy, float), np.asarray(gt_xy, float)) <= radius
    cols = np.where(within.any(axis=0))[0]
    if len(cols) == 0:
        return 0
    w = within[:, cols]
    r, c = linear_sum_assignment(-w.astype(float))
    return int(w[r, c].sum())


def sign_flip(deltas):
    """
        Exact two-sided sign-flip p on per-ROI integer deltas, by bitmask enumeration.

        deltas (array-like): one integer delta per ROI.

        Returns tuple[Fraction, int]: exact p as a fraction, and the number of nonzero deltas.
    """
    d = np.asarray(deltas, dtype=np.int64)
    nz = np.abs(d[d != 0])
    n = len(nz)
    if n == 0:
        return Fraction(1, 1), 0
    obs = abs(int(d.sum()))
    bits = (np.arange(2 ** n)[:, None] >> np.arange(n)) & 1
    sums = ((2 * bits - 1) * nz).sum(axis=1)
    return Fraction(int((np.abs(sums) >= obs).sum()), 2 ** n), n


def holm(pvals, alpha=0.05):
    """
        Holm step-down adjusted p-values.

        pvals (list[float]): raw p-values.
        alpha (float): family-wise level.

        Returns tuple[np.ndarray, np.ndarray]: adjusted p per input position, reject flag per input position.
    """
    p = np.asarray(pvals, float)
    m = len(p)
    order = np.argsort(p, kind='stable')
    adj = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, min(1.0, (m - rank) * p[idx]))
        adj[idx] = running
    return adj, adj <= alpha


# --------------------------------------------------------------------------------------------------
# Tier A
# --------------------------------------------------------------------------------------------------

def execution_gate():
    """
        Execution-coherence gate on the target notebook.

        Returns pd.DataFrame: one row per code cell.
    """
    nb = json.load(open(NB_PATH))
    rows = []
    for i, c in enumerate(nb['cells']):
        if c['cell_type'] != 'code':
            continue
        ex = c.get('metadata', {}).get('execution', {})
        rows.append(dict(cell=i, execution_count=c.get('execution_count'), n_outputs=len(c.get('outputs', [])), has_error=any(o['output_type'] == 'error' for o in c.get('outputs', [])), started_utc=ex.get('iopub.execute_input'), n_png=sum('image/png' in o.get('data', {}) for o in c.get('outputs', []))))
    df = pd.DataFrame(rows)
    counts = df['execution_count'].tolist()
    df['contiguous_from_1'] = counts == list(range(1, len(counts) + 1))
    return df


def provenance():
    """
        mtime, git tracking and working-tree status for the notebook, its artifacts, its imports and the reference CSV.

        Returns pd.DataFrame: one row per file.
    """
    files = [NB_PATH, REF_NB_PATH] + [f'{STEM}_{s}.csv' for s in ('per_roi', 'summary', 'by_domain', 'delta_per_roi', 'delta_stats', 'top30', 'verification')] + [REF_RAW] + [f'{REPO}/production_hematoxylin_only/bbox3way_{s}.png' for s in ('pooled_precision', 'roi_delta_heatmap', 'win_counts')] + [f'{REPO}/midog_utils/{m}.py' for m in ('channels', 'dataset', 'evaluate', 'production', 'seed_selection', 'template_match', 'find_and_suppress', 'nms', 'chromatin', 'invariants')] + [f'{REPO}/DECISIONS.md', f'{REPO}/D8_TEMPLATE_ANCHOR.md', f'{REPO}/PRODUCTION_PIPELINE_CLEANUP.md']
    rows = []
    for f in files:
        rel = os.path.relpath(f, REPO)
        status = subprocess.run(['git', '-C', REPO, 'status', '--porcelain', '--', rel], capture_output=True, text=True).stdout.strip()
        last = subprocess.run(['git', '-C', REPO, 'log', '-1', '--format=%h %ad', '--date=iso', '--', rel], capture_output=True, text=True).stdout.strip()
        rows.append(dict(file=rel, mtime=time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(f))), git_status=status[:2] if status else 'clean', last_commit=last or 'never committed'))
    return pd.DataFrame(rows)


def config_drift():
    """
        Three-way config comparison: production.py constants, FSConfig defaults, and the decision documents.

        Returns pd.DataFrame: one row per setting.
    """
    prod = {n.targets[0].id: ast.literal_eval(n.value) for n in ast.parse(open(f'{REPO}/midog_utils/production.py').read()).body if isinstance(n, ast.Assign) and isinstance(n.targets[0], ast.Name) and isinstance(n.value, (ast.Constant, ast.UnaryOp))}
    src_fs = open(f'{REPO}/midog_utils/find_and_suppress.py').read()
    fs_defaults = {}
    for node in ast.walk(ast.parse(src_fs)):
        if isinstance(node, ast.ClassDef) and node.name == 'FSConfig':
            for st in node.body:
                if isinstance(st, ast.AnnAssign) and isinstance(st.target, ast.Name):
                    fs_defaults[st.target.id] = ast.unparse(st.value)
    prod_src = open(f'{REPO}/midog_utils/production.py').read()
    rows = [
        dict(setting='search channel', production_py=prod.get('CHANNEL'), fsconfig_default=fs_defaults.get('channel'), decision='hematoxylin_od (D3)', notebook='hematoxylin_od via production (asserted cell 1)'),
        dict(setting='tm_method', production_py=re.search(r'TM_METHOD = (cv2\.\w+)', prod_src).group(1), fsconfig_default=fs_defaults.get('tm_method'), decision='TM_CCOEFF (D1, DECISIONS_UNVERIFIED.md)', notebook='TM_CCOEFF (asserted, value 4)'),
        dict(setting='peak_min_distance', production_py=prod.get('PEAK_MIN_DISTANCE'), fsconfig_default=fs_defaults.get('peak_min_distance'), decision='not in D-entries', notebook=7),
        dict(setting='self_hit_radius px', production_py=prod.get('SELF_HIT_RADIUS'), fsconfig_default=fs_defaults.get('self_hit_radius'), decision='D8 cost 1: at template centre', notebook=5.0),
        dict(setting='deep_floor_z', production_py=prod.get('DEEP_FLOOR_Z'), fsconfig_default=fs_defaults.get('deep_floor_z'), decision='-1.5 unchanged (D9)', notebook=-1.5),
        dict(setting='max_peaks (pre-NMS)', production_py=prod.get('MAX_PEAKS'), fsconfig_default=fs_defaults.get('max_peaks'), decision='100 (D9), validated K<=30 on tm_score arm only', notebook=100),
        dict(setting='nms_radius', production_py='ev.radius_px(mpp)', fsconfig_default=fs_defaults.get('nms_radius'), decision='7.5 um per image (D7)', notebook='7.5 um (asserted)'),
        dict(setting='border_pad', production_py='True', fsconfig_default=fs_defaults.get('border_pad'), decision='not in D-entries', notebook='via production'),
        dict(setting='augmentations', production_py='scales=(1.0,), n_angles=1, flips=(False,)', fsconfig_default=f"scales={fs_defaults.get('scales')}, n_angles={fs_defaults.get('n_angles')}, flips={fs_defaults.get('flips')}", decision='not in D-entries', notebook='via production'),
        dict(setting='od window', production_py=f"OD_WINDOW = tm.BASE_SIZE (51)", fsconfig_default='n/a', decision='od51 (D5)', notebook=51),
        dict(setting='ranking key', production_py='rank_key default "tm_score"', fsconfig_default='n/a', decision='D5: tm_score is production; PRODUCTION_PIPELINE_CLEANUP.md working assumption: chromatin_od', notebook='chromatin_od'),
        dict(setting='seed gate channel', production_py='run_pipeline.py: build_seed(pool, gray_inv, ...)', fsconfig_default='n/a', decision='D8: Otsu on structural channel (gray_inverted in build_seed callers)', notebook='joint: gray_inverted AND hematoxylin_od'),
        dict(setting='otsu gate params', production_py='tighten_box_otsu defaults: binary, min_area 50, max_area_frac 0.85, min_solidity 0.5, center_tolerance 0', fsconfig_default='n/a', decision='same (D8 item 1)', notebook='defaults, otsu_window=51'),
    ]
    return pd.DataFrame(rows)


def build_gt(ann, image_id, seed_ann_id):
    """
        gt_eval for one ROI: every annotation of the image except the click.

        ann (pd.DataFrame): all annotations.
        image_id (int): the ROI's image id.
        seed_ann_id (int): the clicked annotation.

        Returns pd.DataFrame: the evaluation ground truth, both categories.
    """
    g = ann[ann['image_id'] == image_id]
    return g[g['ann_id'] != seed_ann_id].reset_index(drop=True)


def score_list(det_xy, gt, radius, k_list=BUDGETS):
    """
        TP@K for one ranked list under every scoring rule the audit considers.

        det_xy (np.ndarray): (N, 2) ranked detections.
        gt (pd.DataFrame): evaluation GT (seed removed), both categories.
        radius (float): match radius in px.
        k_list (tuple[int]): budgets.

        Returns tuple[dict, np.ndarray, int]: rule -> K -> TP, the repo-rule GT index per detection, and exact-distance ties.
    """
    all_xy = gt[['cx', 'cy']].to_numpy(float)
    cat = gt['category_id'].to_numpy()
    unanimous = (gt['n_mitotic_votes'] == gt['n_votes']).to_numpy()
    mit = cat == MITOTIC
    mit_xy = all_xy[mit]
    idx_mixed, ties = greedy(det_xy, all_xy, radius)
    idx_mit, _ = greedy(det_xy, mit_xy, radius)
    idx_strict, _ = greedy(det_xy, all_xy, radius, strict=True)
    idx_rm1, _ = greedy(det_xy, all_xy, radius - 1.0)
    idx_rp1, _ = greedy(det_xy, all_xy, radius + 1.0)
    keep_unan = ~mit | unanimous
    sub_xy, sub_cat = all_xy[keep_unan], cat[keep_unan]
    idx_unan, _ = greedy(det_xy, sub_xy, radius)
    d_mit = cdist(det_xy, mit_xy) if len(mit_xy) else np.full((len(det_xy), 1), np.inf)
    out = {}
    for k in k_list:
        top = slice(0, k)
        out.setdefault('repo_mixed_greedy', {})[k] = int(sum(1 for g in idx_mixed[top] if g >= 0 and cat[g] == MITOTIC))
        out.setdefault('mitotic_only_greedy', {})[k] = int((idx_mit[top] >= 0).sum())
        out.setdefault('mitotic_only_max_matching', {})[k] = max_matching(det_xy[top], mit_xy, radius)
        out.setdefault('any_mitosis_within_radius', {})[k] = int((d_mit[top] <= radius).any(axis=1).sum())
        out.setdefault('contested_mitoses_dropped', {})[k] = int(sum(1 for g in idx_unan[top] if g >= 0 and sub_cat[g] == MITOTIC))
        out.setdefault('strict_less_than', {})[k] = int(sum(1 for g in idx_strict[top] if g >= 0 and cat[g] == MITOTIC))
        out.setdefault('radius_minus_1px', {})[k] = int(sum(1 for g in idx_rm1[top] if g >= 0 and cat[g] == MITOTIC))
        out.setdefault('radius_plus_1px', {})[k] = int(sum(1 for g in idx_rp1[top] if g >= 0 and cat[g] == MITOTIC))
        out.setdefault('topk_only_repo_greedy', {})[k] = int(sum(1 for g in greedy(det_xy[top], all_xy, radius)[0] if g >= 0 and cat[g] == MITOTIC))
        out.setdefault('repo_tp_on_contested_mitoses', {})[k] = int(sum(1 for g in idx_mixed[top] if g >= 0 and cat[g] == MITOTIC and not unanimous[g]))
    return out, idx_mixed, ties


def delta_stats_table(tp, rule):
    """
        Paired per-ROI deltas, win/loss/tie, exact sign-flip p and Holm over the 9 tests, for one scoring rule.

        tp (pd.DataFrame): columns file_name, condition, rule, K, tp.
        rule (str): which scoring rule to read.

        Returns pd.DataFrame: one row per (pair, K).
    """
    sub = tp[tp['rule'] == rule].pivot_table(index=['file_name', 'K'], columns='condition', values='tp').reset_index()
    rows = []
    for a, b in PAIRS:
        for k in BUDGETS:
            s = sub[sub['K'] == k].sort_values('file_name')
            d = (s[a] - s[b]).astype(int).to_numpy()
            p, n = sign_flip(d)
            rows.append(dict(rule=rule, pair=f'{a} - {b}', K=k, mean_delta_precision=d.mean() / k, wins=int((d > 0).sum()), losses=int((d < 0).sum()), ties=int((d == 0).sum()), n_nonzero=n, exact_p=float(p), exact_p_fraction=str(p), p_floor=float(Fraction(2, 2 ** n)) if n else 1.0))
    df = pd.DataFrame(rows)
    adj, rej = holm(df['exact_p'].tolist())
    df['holm_adjusted_p'] = adj
    df['holm_reject_0.05'] = rej
    return df


def sign_flip_float(deltas, eps=1e-12):
    """
        Exact two-sided sign-flip p on real-valued per-ROI deltas.

        deltas (array-like): one delta per ROI.
        eps (float): tolerance for zero and for ties with the observed sum.

        Returns tuple[float, int]: p-value and the number of nonzero deltas.
    """
    d = np.asarray(deltas, float)
    nz = np.abs(d[np.abs(d) > eps])
    n = len(nz)
    if n == 0:
        return 1.0, 0
    obs = abs(d.sum())
    bits = (np.arange(2 ** n)[:, None] >> np.arange(n)) & 1
    sums = ((2 * bits - 1) * nz).sum(axis=1)
    return float(np.mean(np.abs(sums) >= obs - 1e-9)), n


def contested_neutral_stats(tp):
    """
        Delta statistics when TPs on 2-of-3 (contested) mitoses are neither credited nor charged: precision = unanimous TPs / (K - contested TPs).

        tp (pd.DataFrame): long TP table with the repo rule and the contested-TP helper rule.

        Returns pd.DataFrame: one row per (pair, K), same columns as delta_stats_table where they apply.
    """
    rep = tp[tp['rule'] == 'repo_mixed_greedy'].set_index(['file_name', 'condition', 'K'])['tp']
    con = tp[tp['rule'] == 'repo_tp_on_contested_mitoses'].set_index(['file_name', 'condition', 'K'])['tp']
    prec = ((rep - con) / (rep.index.get_level_values('K') - con)).rename('p').reset_index()
    rows = []
    for a, b in PAIRS:
        for k in BUDGETS:
            s = prec[prec['K'] == k].pivot_table(index='file_name', columns='condition', values='p')
            d = (s[a] - s[b]).to_numpy()
            p, n = sign_flip_float(d)
            rows.append(dict(rule='contested_neutral', pair=f'{a} - {b}', K=k, mean_delta_precision=float(d.mean()), wins=int((d > 1e-12).sum()), losses=int((d < -1e-12).sum()), ties=int((np.abs(d) <= 1e-12).sum()), n_nonzero=n, exact_p=p, exact_p_fraction='', p_floor=2.0 / 2 ** n if n else 1.0))
    df = pd.DataFrame(rows)
    adj, rej = holm(df['exact_p'].tolist())
    df['holm_adjusted_p'] = adj
    df['holm_reject_0.05'] = rej
    return df


def tier_a(images, ann):
    """
        Every Tier A recomputation.

        images (pd.DataFrame): image table from the JSON.
        ann (pd.DataFrame): annotation table from the JSON.

        Returns dict: the frames Tier B and the prose checks reuse.
    """
    log('\n== Tier A ==')
    ex = execution_gate()
    save(ex, 'execution_gate')
    log(f"  execution counts contiguous 1..{len(ex)}: {bool(ex['contiguous_from_1'].iloc[0])}; errors: {int(ex['has_error'].sum())}; unrun: {int(ex['execution_count'].isna().sum())}")
    save(provenance(), 'provenance')
    save(config_drift(), 'config_drift')

    per_roi = pd.read_csv(f'{STEM}_per_roi.csv')
    summary = pd.read_csv(f'{STEM}_summary.csv')
    by_domain = pd.read_csv(f'{STEM}_by_domain.csv')
    dpr = pd.read_csv(f'{STEM}_delta_per_roi.csv')
    dstats = pd.read_csv(f'{STEM}_delta_stats.csv')
    top = pd.read_csv(f'{STEM}_top30.csv')
    verif = pd.read_csv(f'{STEM}_verification.csv')

    # ---- composition ----------------------------------------------------------------------------
    files = sorted(f for f in os.listdir(IMG_DIR) if f.endswith('.tiff'))
    meta = images.set_index('file_name').loc[files]
    comp = [
        dict(check='ROIs on disk', expected=14, observed=len(files)),
        dict(check='per_roi rows (14 x 3)', expected=42, observed=len(per_roi)),
        dict(check='per_roi unique (file, condition)', expected=42, observed=int(per_roi[['file_name', 'condition']].drop_duplicates().shape[0])),
        dict(check='top30 rows (42 x 30)', expected=1260, observed=len(top)),
        dict(check='top30 unique (file, condition, rank)', expected=1260, observed=int(top[['file_name', 'condition', 'rank']].drop_duplicates().shape[0])),
        dict(check='top30 ranks are 0..29 in each run', expected=42, observed=int(sum(sorted(g['rank'].tolist()) == list(range(30)) for _, g in top.groupby(['file_name', 'condition'])))),
        dict(check='top30 row order == rank order in each run', expected=42, observed=int(sum(g['rank'].tolist() == list(range(30)) for _, g in top.groupby(['file_name', 'condition'], sort=False)))),
        dict(check='summary rows (3 x 3)', expected=9, observed=len(summary)),
        dict(check='by_domain rows (7 x 3 x 3)', expected=63, observed=len(by_domain)),
        dict(check='delta_per_roi rows (14 x 3 x 3)', expected=126, observed=len(dpr)),
        dict(check='delta_stats rows (3 x 3)', expected=9, observed=len(dstats)),
        dict(check='verification rows (5x42 + 14 + 14 + 28 + 1)', expected=267, observed=len(verif)),
        dict(check='verification all passed', expected=267, observed=int(verif['passed'].sum())),
        dict(check='domains in JSON for the 14 ROIs', expected=7, observed=int(meta['tumor_type'].nunique())),
        dict(check='min ROIs per domain', expected=2, observed=int(meta['tumor_type'].value_counts().min())),
        dict(check='max ROIs per domain', expected=2, observed=int(meta['tumor_type'].value_counts().max())),
        dict(check='all 26,286 JSON boxes 50x50', expected=len(ann), observed=int(((ann['w'] == 50) & (ann['h'] == 50)).sum())),
    ]
    comp = pd.DataFrame(comp)
    comp['passed'] = comp['expected'] == comp['observed']
    save(comp, 'composition')

    # ---- mpp and radius from the TIFF tags ------------------------------------------------------
    tmeta = {fn: tiff_meta(f'{IMG_DIR}/{fn}') for fn in files}
    for fn in files:
        r = per_roi[per_roi['file_name'] == fn]
        compare('mpp', fn, 'n_gt_mitotic identical across conditions', 1, r['n_gt_mitotic'].nunique())
        compare('tiff', fn, 'width JSON vs TIFF', int(meta.loc[fn, 'width']), tmeta[fn]['width'])
        compare('tiff', fn, 'height JSON vs TIFF', int(meta.loc[fn, 'height']), tmeta[fn]['height'])

    aniso = pd.DataFrame([dict(file_name=fn, mpp_x=tmeta[fn]['mpp_x'], mpp_y=tmeta[fn]['mpp_y'], radius_px_from_x=RADIUS_UM / tmeta[fn]['mpp_x'], radius_px_from_y=RADIUS_UM / tmeta[fn]['mpp_y'], radius_gap_px=RADIUS_UM / tmeta[fn]['mpp_x'] - RADIUS_UM / tmeta[fn]['mpp_y'], resolution_unit=tmeta[fn]['unit']) for fn in files])
    save(aniso, 'tiff_resolution_tags')

    # ---- TP re-derivation on the persisted top-30 ------------------------------------------------
    tp_rows, det_rows, run_rows = [], [], []
    for (fn, cond), g in top.groupby(['file_name', 'condition'], sort=False):
        pr = per_roi[(per_roi['file_name'] == fn) & (per_roi['condition'] == cond)].iloc[0]
        image_id = int(meta.loc[fn, 'image_id'])
        seed = int(pr['seed_ann_id'])
        gt = build_gt(ann, image_id, seed)
        mpp = tmeta[fn]['mpp_x']
        radius = RADIUS_UM / mpp
        g = g.sort_values('rank')
        det_xy = g[['cx', 'cy']].to_numpy(float)
        scores, idx_mixed, ties = score_list(det_xy, gt, radius)
        for rule, kv in scores.items():
            for k, v in kv.items():
                tp_rows.append(dict(file_name=fn, condition=cond, rule=rule, K=k, tp=v))
        for k in BUDGETS:
            compare('per_roi', f'{fn}/{cond}', f'tp_at_{k} (independent repo-rule rematch)', int(pr[f'tp_at_{k}']), scores['repo_mixed_greedy'][k])
            compare('per_roi', f'{fn}/{cond}', f'tp_at_{k} (bucket column count)', int(pr[f'tp_at_{k}']), int((g['bucket'].iloc[:k] == 'human_correct_label').sum()))
        # per-detection geometry
        cat = gt['category_id'].to_numpy()
        unan = (gt['n_mitotic_votes'] == gt['n_votes']).to_numpy()
        d_all = cdist(det_xy, gt[['cx', 'cy']].to_numpy(float))
        seed_row = ann[ann['ann_id'] == seed].iloc[0]
        d_click = np.hypot(det_xy[:, 0] - seed_row['cx'], det_xy[:, 1] - seed_row['cy'])
        bucket_audit = []
        for i in range(len(det_xy)):
            gidx = idx_mixed[i]
            b = 'non_human_findings' if gidx < 0 else ('human_correct_label' if cat[gidx] == MITOTIC else 'human_rejected_label')
            bucket_audit.append(b)
            dm = d_all[i][cat == MITOTIC]
            dl = d_all[i][cat == LOOKALIKE]
            det_rows.append(dict(file_name=fn, condition=cond, rank=int(g['rank'].iloc[i]), cx=det_xy[i, 0], cy=det_xy[i, 1], od=float(g['od'].iloc[i]), nb_bucket=g['bucket'].iloc[i], audit_bucket=b, nb_matched_ann=int(g['matched_ann_id'].iloc[i]), audit_matched_ann=int(gt['ann_id'].iloc[gidx]) if gidx >= 0 else -1, matched_contested_mitosis=bool(gidx >= 0 and cat[gidx] == MITOTIC and not unan[gidx]), nearest_mitotic_px=float(dm.min()) if len(dm) else np.inf, nearest_lookalike_px=float(dl.min()) if len(dl) else np.inf, n_mitotic_within_r=int((dm <= radius).sum()), n_lookalike_within_r=int((dl <= radius).sum()), dist_to_click_px=float(d_click[i]), radius_px=radius))
            compare('top30', f'{fn}/{cond}/rank{i}', 'bucket', g['bucket'].iloc[i], b)
            compare('top30', f'{fn}/{cond}/rank{i}', 'matched_ann_id', int(g['matched_ann_id'].iloc[i]), int(gt['ann_id'].iloc[gidx]) if gidx >= 0 else -1)
        run_rows.append(dict(file_name=fn, condition=cond, radius_px=radius, mpp=mpp, n_gt_mitotic=int((cat == MITOTIC).sum()), n_gt_lookalike=int((cat == LOOKALIKE).sum()), exact_distance_ties=ties, min_abs_margin_to_radius_px=float(np.abs(d_all - radius).min()) if d_all.size else np.nan, od_descending=bool(np.all(np.diff(g['od'].to_numpy()) <= 0)), top30_within_r_of_click=int((d_click <= radius).sum())))
        compare('per_roi', f'{fn}/{cond}', 'n_gt_mitotic', int(pr['n_gt_mitotic']), int((cat == MITOTIC).sum()))
    tp = pd.DataFrame(tp_rows)
    det = pd.DataFrame(det_rows)
    runs = pd.DataFrame(run_rows)
    save(det, 'top30_rematch_detections')

    # scoring-rule sensitivity per run
    wide = tp.pivot_table(index=['file_name', 'condition', 'K'], columns='rule', values='tp').reset_index()
    base = wide['repo_mixed_greedy']
    for rule in [c for c in wide.columns if c not in ('file_name', 'condition', 'K', 'repo_mixed_greedy')]:
        wide[f'diff_{rule}'] = wide[rule] - base
    save(wide, 'tp_by_scoring_rule_per_run')
    save(runs, 'top30_rematch_runs')

    # seed identity and seed pool membership
    seeds = per_roi.groupby('file_name')['seed_ann_id'].first()
    click_rows = []
    for fn, sid in seeds.items():
        a = ann[ann['ann_id'] == sid].iloc[0]
        compare('clicks', fn, 'seed belongs to this image', int(meta.loc[fn, 'image_id']), int(a['image_id']))
        compare('clicks', fn, 'seed is mitotic', MITOTIC, int(a['category_id']))
        compare('clicks', fn, 'seed is unanimous', int(a['n_votes']), int(a['n_mitotic_votes']))
        others = ann[(ann['image_id'] == a['image_id']) & (ann['ann_id'] != sid)]
        d = np.hypot(others['cx'] - a['cx'], others['cy'] - a['cy'])
        r = RADIUS_UM / tmeta[fn]['mpp_x']
        mo = others['category_id'] == MITOTIC
        click_rows.append(dict(file_name=fn, seed_ann_id=int(sid), click_cx=a['cx'], click_cy=a['cy'], radius_px=r, nearest_other_mitotic_px=float(d[mo].min()), nearest_lookalike_px=float(d[~mo].min()) if (~mo).any() else np.inf, n_other_annotations_within_2r=int((d <= 2 * r).sum()), n_other_mitotic_within_r=int((d[mo] <= r).sum())))
    save(pd.DataFrame(click_rows), 'click_neighbourhood')

    # ---- tables from the independently matched TPs ---------------------------------------------
    rep = tp[tp['rule'] == 'repo_mixed_greedy'].pivot_table(index=['file_name', 'condition'], columns='K', values='tp')
    rep.columns = [f'tp_at_{k}' for k in rep.columns]
    rep = rep.reset_index().merge(images[['file_name', 'tumor_type']], on='file_name')
    for _, r in per_roi.iterrows():
        a = rep[(rep['file_name'] == r['file_name']) & (rep['condition'] == r['condition'])].iloc[0]
        for k in BUDGETS:
            compare('per_roi', f"{r['file_name']}/{r['condition']}", f'precision_at_{k}', r[f'precision_at_{k}'], a[f'tp_at_{k}'] / k, tol=1e-12)
        compare('per_roi', f"{r['file_name']}/{r['condition']}", 'domain', r['domain'], a['tumor_type'])
    for _, r in summary.iterrows():
        s = rep[rep['condition'] == r['condition']]
        k = int(r['K'])
        p = s[f'tp_at_{k}'] / k
        key = f"{r['condition']}/K{k}"
        compare('summary', key, 'n_roi', r['n_roi'], len(s))
        compare('summary', key, 'tp_sum', r['tp_sum'], int(s[f'tp_at_{k}'].sum()))
        compare('summary', key, 'delivered_sum', r['delivered_sum'], k * len(s))
        compare('summary', key, 'precision_pooled', r['precision_pooled'], round(s[f'tp_at_{k}'].sum() / (k * len(s)), 4), tol=1e-12)
        compare('summary', key, 'pooled == mean of per-ROI precision', round(s[f'tp_at_{k}'].sum() / (k * len(s)), 12), round(p.mean(), 12), tol=1e-12)
        compare('summary', key, 'worst_roi_precision', r['worst_roi_precision'], round(p.min(), 4), tol=1e-12)
        compare('summary', key, 'n_roi_at_worst', r['n_roi_at_worst'], int((p == p.min()).sum()))
        compare('summary', key, 'worst_roi_files', r['worst_roi_files'], ';'.join(sorted(s.loc[p == p.min(), 'file_name'])))
        compare('summary', key, 'best_roi_precision', r['best_roi_precision'], round(p.max(), 4), tol=1e-12)
    for _, r in by_domain.iterrows():
        s = rep[(rep['condition'] == r['condition']) & (rep['tumor_type'] == r['domain'])]
        k = int(r['K'])
        key = f"{r['domain']}/{r['condition']}/K{k}"
        compare('by_domain', key, 'n_roi', r['n_roi'], len(s))
        compare('by_domain', key, 'tp_sum', r['tp_sum'], int(s[f'tp_at_{k}'].sum()))
        compare('by_domain', key, 'precision_pooled', r['precision_pooled'], round(s[f'tp_at_{k}'].sum() / (k * len(s)), 4), tol=1e-12)
        compare('by_domain', key, 'worst_roi_precision', r['worst_roi_precision'], round((s[f'tp_at_{k}'] / k).min(), 4), tol=1e-12)
    wide_tp = rep.set_index(['file_name', 'condition'])
    for _, r in dpr.iterrows():
        a, b = r['pair'].split(' - ')
        k = int(r['K'])
        d = int(wide_tp.loc[(r['file_name'], a), f'tp_at_{k}'] - wide_tp.loc[(r['file_name'], b), f'tp_at_{k}'])
        key = f"{r['file_name']}/{r['pair']}/K{k}"
        compare('delta_per_roi', key, 'delta_tp', r['delta_tp'], d)
        compare('delta_per_roi', key, 'delta_precision', r['delta_precision'], d / k, tol=1e-12)
        compare('delta_per_roi', key, 'base_size_a', r['base_size_a'], int(per_roi[(per_roi['file_name'] == r['file_name']) & (per_roi['condition'] == a)]['base_size'].iloc[0]))
        compare('delta_per_roi', key, 'domain', r['domain'], images.set_index('file_name').loc[r['file_name'], 'tumor_type'])

    stats_all = pd.concat([delta_stats_table(tp, rule) for rule in tp['rule'].unique() if rule != 'repo_tp_on_contested_mitoses'] + [contested_neutral_stats(tp)], ignore_index=True)
    save(stats_all, 'delta_stats_by_scoring_rule')
    rs = stats_all[stats_all['rule'] == 'repo_mixed_greedy'].set_index(['pair', 'K'])
    for _, r in dstats.iterrows():
        a = rs.loc[(r['pair'], int(r['K']))]
        key = f"{r['pair']}/K{int(r['K'])}"
        compare('delta_stats', key, 'mean_delta_precision', r['mean_delta_precision'], round(a['mean_delta_precision'], 4), tol=1e-12)
        for c in ('wins', 'losses', 'ties', 'n_nonzero'):
            compare('delta_stats', key, c, r[c], a[c])
        compare('delta_stats', key, 'exact_p', r['exact_p'], round(a['exact_p'], 4), tol=1e-12)

    # premise check: the same tests with the tumour domain (7 units, 2 ROIs each) as the exchangeable unit
    dom_rows = []
    for x, y in PAIRS:
        for k in BUDGETS:
            dd = rep.pivot_table(index=['tumor_type', 'file_name'], columns='condition', values=f'tp_at_{k}')
            per_dom = (dd[x] - dd[y]).groupby(level='tumor_type').sum().astype(int).to_numpy()
            pd_, nd = sign_flip(per_dom)
            dom_rows.append(dict(pair=f'{x} - {y}', K=k, unit='domain (sum of 2 ROI deltas)', n_units=len(per_dom), n_nonzero=nd, exact_p=float(pd_), p_floor=2.0 / 2 ** nd if nd else 1.0, roi_level_p=float(rs_p.loc[(f'{x} - {y}', k)]) if (rs_p := stats_all[stats_all['rule'] == 'repo_mixed_greedy'].set_index(['pair', 'K'])['exact_p']) is not None else np.nan))
    dom = pd.DataFrame(dom_rows)
    dom['holm_adjusted_p'], dom['holm_reject_0.05'] = holm(dom['exact_p'].tolist())
    save(dom, 'domain_level_sign_flip')

    # figure labels as the notebook's own format strings would render them from the recomputed values;
    # the decoded PNGs were read against this table (the render itself is not machine-parsed)
    roi_order = images.set_index('file_name').loc[files].reset_index().sort_values(['tumor_type', 'file_name'])['file_name'].tolist()
    fig_rows = []
    for c in CONDITIONS:
        for k in BUDGETS:
            v = round(rep[rep['condition'] == c][f'tp_at_{k}'].sum() / (14 * k), 4)
            fig_rows.append(dict(figure='1 pooled precision bars', panel=c, row=f'top {k}', label=f'{v:.3f}'))
    for x, y in PAIRS:
        for fn in roi_order:
            for k in BUDGETS:
                d = int(wide_tp.loc[(fn, x), f'tp_at_{k}'] - wide_tp.loc[(fn, y), f'tp_at_{k}']) / k
                fig_rows.append(dict(figure='2 per-ROI delta heatmap', panel=f'{x} minus {y}', row=f'{fn} K={k}', label='0' if d == 0 else f'{d:+.2f}'))
        for k in BUDGETS:
            r = stats_all[(stats_all['rule'] == 'repo_mixed_greedy') & (stats_all['pair'] == f'{x} - {y}') & (stats_all['K'] == k)].iloc[0]
            for part in ('wins', 'losses', 'ties'):
                fig_rows.append(dict(figure='3 win/loss/tie stacks', panel=f'{x} vs {y}', row=f'K={k} {part}', label=str(int(r[part])) if r[part] > 0 else '(not drawn)'))
            fig_rows.append(dict(figure='3 win/loss/tie stacks', panel=f'{x} vs {y}', row=f'K={k} p label', label=f"p={round(r['exact_p'], 4):.3f}"))
    save(pd.DataFrame(fig_rows), 'figure_labels_expected')

    # pooled precision per rule (does the ordering survive every rule?)
    pooled = tp[tp['rule'] != 'repo_tp_on_contested_mitoses'].groupby(['rule', 'condition', 'K'])['tp'].sum().reset_index()
    pooled['precision_pooled'] = pooled['tp'] / (pooled['K'] * 14)
    save(pooled, 'pooled_precision_by_scoring_rule')

    # ---- ties vs template difference ------------------------------------------------------------
    geo = per_roi.set_index(['file_name', 'condition'])
    tie_rows = []
    rep_idx = rep.set_index(['file_name', 'condition'])
    for a, b in PAIRS:
        for fn in files:
            dbase = abs(int(geo.loc[(fn, a), 'base_size']) - int(geo.loc[(fn, b), 'base_size']))
            ta = per_roi[(per_roi['file_name'] == fn) & (per_roi['condition'] == a)]
            tb = per_roi[(per_roi['file_name'] == fn) & (per_roi['condition'] == b)]
            ga = top[(top['file_name'] == fn) & (top['condition'] == a)].sort_values('rank')
            gb = top[(top['file_name'] == fn) & (top['condition'] == b)].sort_values('rank')
            for k in BUDGETS:
                sa = set(map(tuple, ga[['cx', 'cy']].to_numpy()[:k]))
                sb = set(map(tuple, gb[['cx', 'cy']].to_numpy()[:k]))
                ma = set(ga['matched_ann_id'].iloc[:k][ga['bucket'].iloc[:k] == 'human_correct_label'])
                mb = set(gb['matched_ann_id'].iloc[:k][gb['bucket'].iloc[:k] == 'human_correct_label'])
                tie_rows.append(dict(pair=f'{a} - {b}', file_name=fn, K=k, abs_delta_base_size=dbase, tie=bool(rep_idx.loc[(fn, a), f'tp_at_{k}'] == rep_idx.loc[(fn, b), f'tp_at_{k}']), shared_coordinates=len(sa & sb), shared_mitoses=len(ma & mb), tp_a=len(ma), tp_b=len(mb)))
    ties = pd.DataFrame(tie_rows)
    save(ties, 'tie_vs_template_difference')
    nb = json.load(open(NB_PATH))
    text = ''.join(''.join(o.get('text', '')) for o in nb['cells'][4]['outputs'] if o['output_type'] == 'stream')
    nb_retries = {m.group(1): int(m.group(2)) for m in re.finditer(r'\[(\d{3}\.tiff)\] .*? n_retries=(\d+)', text)}
    assert len(nb_retries) == 14
    return dict(per_roi=per_roi, top=top, tp=tp, stats_all=stats_all, tmeta=tmeta, meta=meta, files=files, ties=ties, det=det, summary=summary, dstats=dstats, rep=rep, nb_retries=nb_retries)


def gate2_against_current_raw(per_roi):
    """
        Re-run the notebook's Gate 2 against the reference raw CSV as it is on disk now.

        per_roi (pd.DataFrame): the notebook's per-ROI table.

        Returns pd.DataFrame: one row per (ROI, Otsu condition).
    """
    ref = pd.read_csv(REF_RAW)
    ref = ref[ref['arm'] == 'chromatin_od']
    rows = []
    for _, r in per_roi[per_roi['condition'].isin(['gray_bbox', 'hem_bbox'])].iterrows():
        s = ref[(ref['file_name'] == r['file_name']) & (ref['condition'] == r['condition'])].set_index('budget')
        row = dict(file_name=r['file_name'], condition=r['condition'], base_size_nb=r['base_size'], base_size_ref=int(s['base_size'].iloc[0]), n_det_nb=r['n_detections'], n_det_ref=int(s['n_detections'].iloc[0]))
        for k in BUDGETS:
            row[f'tp_at_{k}_nb'] = int(r[f'tp_at_{k}'])
            row[f'tp_at_{k}_ref'] = int(s.loc[k, 'tp_at_budget'])
        row['all_match'] = all(row[f'{c}_nb'] == row[f'{c}_ref'] for c in ['base_size', 'n_det'] + [f'tp_at_{k}' for k in BUDGETS])
        rows.append(row)
        for c in ['base_size', 'n_det'] + [f'tp_at_{k}' for k in BUDGETS]:
            compare('gate2_current_raw', f"{r['file_name']}/{r['condition']}", c, row[f'{c}_nb'], row[f'{c}_ref'])
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------------------
# Tier B -- pixels
# --------------------------------------------------------------------------------------------------

def to_hem_od(rgb):
    """
        Unclipped hematoxylin OD by skimage colour deconvolution (D3's channel definition).

        rgb (np.ndarray): HxWx3 uint8.

        Returns np.ndarray: HxW float32.
    """
    return rgb2hed(rgb.astype(np.float32) / 255.0)[:, :, 0].astype(np.float32)


def to_gray_inv(rgb):
    """
        255 minus luminance.

        rgb (np.ndarray): HxWx3 uint8.

        Returns np.ndarray: HxW float32.
    """
    return (255.0 - cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)).astype(np.float32)


def readable(shape, x, y, size=PATCH):
    """
        Whether a size x size patch centred on the rounded point fits in the image.

        shape (tuple): image shape.
        x (float): centre x.
        y (float): centre y.
        size (int): patch side.

        Returns bool: True when it fits.
    """
    h, w = shape[:2]
    half = size // 2
    ix, iy = int(round(x)), int(round(y))
    return ix - half >= 0 and iy - half >= 0 and ix + half < w and iy + half < h


def otsu_box(chan, cx, cy):
    """
        D8 gate and template geometry, implemented from the D8 text: Otsu on the 51 px window, component under the click pixel.

        chan (np.ndarray): single-channel image.
        cx (float): click x.
        cy (float): click y.

        Returns tuple[tuple or None, str]: (base_size, centre x, centre y) or None, and the reason.
    """
    half = OTSU_WIN // 2
    if not readable(chan.shape, cx, cy, OTSU_WIN):
        return None, 'window_unreadable'
    ix, iy = int(round(cx)), int(round(cy))
    patch = chan[iy - half: iy + half + 1, ix - half: ix + half + 1].astype(np.float32)
    u8 = cv2.normalize(patch, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    _, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    lab = label(binary, connectivity=2)
    c = lab[half, half]
    if c == 0:
        return None, 'click_not_foreground'
    reg = next(p for p in regionprops(lab) if p.label == c)
    y0, x0, y1, x1 = reg.bbox
    if reg.area < 50:
        return None, 'area_below_50'
    if reg.area > 0.85 * patch.size:
        return None, 'area_above_0.85'
    if reg.solidity < 0.5:
        return None, 'solidity_below_0.5'
    n = int(round(max(y1 - y0, x1 - x0)))
    n = max(5, n + 1 if n % 2 == 0 else n)
    return (n, ix - half + (x0 + x1 - 1) / 2.0, iy - half + (y0 + y1 - 1) / 2.0), 'accepted'


def walk(pool, image_id, predicate):
    """
        The [0, image_id] RNG walk: draw, test, drop on refusal, redraw.

        pool (list[dict]): border-filtered unanimous mitoses, JSON order.
        image_id (int): RNG stream key.
        predicate (callable): row -> bool.

        Returns tuple[dict, int, list[int]]: accepted row, refusals, ann ids visited in order.
    """
    rng = np.random.default_rng([0, image_id])
    working = list(pool)
    visited = []
    while working:
        idx = int(rng.integers(len(working)))
        row = working[idx]
        visited.append(int(row['ann_id']))
        if predicate(row):
            return row, len(visited) - 1, visited
        working.pop(idx)
    raise ValueError('pool exhausted')


def search(hem, base, tx, ty, mpp):
    """
        One production search re-implemented from primitives: template cut, TM_CCOEFF on a replicate-padded image, deep floor, 15x15 local maxima, MAX_PEAKS, NMS at 7.5 um, self-hit at the template centre, od51.

        hem (np.ndarray): hematoxylin OD, float32.
        base (int): odd template side.
        tx (float): template centre x.
        ty (float): template centre y.
        mpp (float): microns per pixel.

        Returns tuple[pd.DataFrame, dict]: NMS-order detections with score and od, and stage counts.
    """
    h, w = hem.shape
    ix, iy = int(round(tx)), int(round(ty))
    hb = base // 2
    tmpl = np.ascontiguousarray(hem[iy - hb: iy + hb + 1, ix - hb: ix + hb + 1], dtype=np.float32)
    pad = (base - 1) // 2
    padded = cv2.copyMakeBorder(hem, pad, pad, pad, pad, cv2.BORDER_REPLICATE)
    res = cv2.matchTemplate(padded, tmpl, cv2.TM_CCOEFF)
    assert res.shape == (h, w)
    res = np.nan_to_num(res, nan=-3.0e38, posinf=-3.0e38, neginf=-3.0e38)
    sample = res[::8, ::8]
    sample = sample[np.isfinite(sample)]
    med = float(np.median(sample))
    mad = float(1.4826 * np.median(np.abs(sample - med)))
    thr = med + DEEP_Z * mad
    mx = ndimage.maximum_filter(res, size=2 * PEAK_MIN_DIST + 1, mode='nearest')
    ys, xs = np.nonzero((res >= mx) & (res >= thr))
    sc = res[ys, xs]
    n_above = len(sc)
    order = np.lexsort((ys, xs, -sc))[:MAX_PEAKS]
    pts = np.stack([xs[order], ys[order]], axis=1).astype(float)
    s = sc[order].astype(np.float64)
    radius = RADIUS_UM / mpp
    dd = cdist(pts, pts)
    suppressed = np.zeros(len(s), dtype=bool)
    keep = []
    for i in np.argsort(-s, kind='stable'):
        if suppressed[i]:
            continue
        keep.append(i)
        suppressed[dd[i] <= radius] = True
    keep = np.asarray(keep)
    pts, s = pts[keep], s[keep]
    self_hit = np.hypot(pts[:, 0] - tx, pts[:, 1] - ty) <= SELF_HIT
    n_self = int(self_hit.sum())
    pts, s = pts[~self_hit], s[~self_hit]
    hp = cv2.copyMakeBorder(hem, OD_WIN // 2, OD_WIN // 2, OD_WIN // 2, OD_WIN // 2, cv2.BORDER_REPLICATE)
    k = max(1, int(OD_FRAC * OD_WIN * OD_WIN))
    od = [float(np.sort(hp[int(y): int(y) + OD_WIN, int(x): int(x) + OD_WIN].ravel())[-k:].astype(np.float64).mean()) for x, y in pts]
    det = pd.DataFrame(dict(cx=pts[:, 0], cy=pts[:, 1], score=s, od=od))
    return det, dict(n_above_floor=n_above, n_peaks=len(order), n_after_nms=len(keep), n_self_hits=n_self, n_detections=len(det), threshold=thr)


def tier_b(images, ann, a):
    """
        Every Tier B recomputation from pixels, on all 14 ROIs x 3 conditions.

        images (pd.DataFrame): image table.
        ann (pd.DataFrame): annotation table.
        a (dict): Tier A outputs.

        Returns dict: Tier B frames.
    """
    log('\n== Tier B ==')
    t_start = time.time()
    per_roi, top, tmeta, meta = a['per_roi'], a['top'], a['tmeta'], a['meta']
    ref = pd.read_csv(REF_RAW)
    draw_rows, run_rows, arm_rows, list_rows = [], [], [], []
    for fn in a['files']:
        t0 = time.time()
        image_id = int(meta.loc[fn, 'image_id'])
        rgb = load_rgb(f'{IMG_DIR}/{fn}')
        mpp = tmeta[fn]['mpp_x']
        radius = RADIUS_UM / mpp
        hem = to_hem_od(rgb)
        gray = to_gray_inv(rgb)
        shape = rgb.shape
        del rgb
        g_img = ann[ann['image_id'] == image_id]
        mit = g_img[g_img['category_id'] == MITOTIC]
        unan = mit[mit['n_mitotic_votes'] == mit['n_votes']]
        assert len(unan) > 0
        ixs, iys = np.rint(unan['cx'].to_numpy()).astype(int), np.rint(unan['cy'].to_numpy()).astype(int)
        border = PATCH // 2
        ok = (ixs >= border) & (ixs <= shape[1] - 1 - border) & (iys >= border) & (iys <= shape[0] - 1 - border)
        pool = unan[ok].to_dict('records')
        cache = {}

        def specs(row):
            """
                Every condition's gate outcome for one candidate click, memoised.

                row (dict): candidate annotation.

                Returns dict: per-channel spec, reason and readability.
            """
            aid = int(row['ann_id'])
            if aid not in cache:
                sg, rg = otsu_box(gray, row['cx'], row['cy'])
                sh, rh = otsu_box(hem, row['cx'], row['cy'])
                cache[aid] = dict(gray=sg, gray_reason=rg, hem=sh, hem_reason=rh, click_readable=readable(shape, row['cx'], row['cy']), gray_readable=sg is not None and readable(shape, sg[1], sg[2]), hem_readable=sh is not None and readable(shape, sh[1], sh[2]))
            return cache[aid]

        joint = lambda r: specs(r)['click_readable'] and specs(r)['gray_readable'] and specs(r)['hem_readable']
        gray_only = lambda r: specs(r)['gray_readable']
        hem_only = lambda r: specs(r)['hem_readable']
        unconditional = lambda r: specs(r)['click_readable']
        row_j, n_ret, visited = walk(pool, image_id, joint)
        row_g, _, _ = walk(pool, image_id, gray_only)
        row_h, _, _ = walk(pool, image_id, hem_only)
        row_u, _, _ = walk(pool, image_id, unconditional)
        for order_i, aid in enumerate(visited):
            r = next(p for p in pool if int(p['ann_id']) == aid)
            sp = specs(r)
            draw_rows.append(dict(file_name=fn, draw=order_i, ann_id=aid, gray_reason=sp['gray_reason'], hem_reason=sp['hem_reason'], click_readable=sp['click_readable'], gray_template_readable=sp['gray_readable'], hem_template_readable=sp['hem_readable'], joint_accept=joint(r)))
        pr = per_roi[per_roi['file_name'] == fn].set_index('condition')
        click = (float(row_j['cx']), float(row_j['cy']))
        sp = specs(row_j)
        cond_spec = dict(default_51=(BASE, click[0], click[1]), gray_bbox=sp['gray'], hem_bbox=sp['hem'])
        compare('tierB_click', fn, 'joint click ann_id', int(pr['seed_ann_id'].iloc[0]), int(row_j['ann_id']))
        compare('tierB_click', fn, 'n_retries (notebook cell 4 print)', a['nb_retries'][fn], n_ret)
        gt = build_gt(ann, image_id, int(row_j['ann_id']))
        cat = gt['category_id'].to_numpy()
        gt_xy = gt[['cx', 'cy']].to_numpy(float)
        mit_xy = gt_xy[cat == MITOTIC]
        for cond in CONDITIONS:
            base, tx, ty = cond_spec[cond]
            compare('tierB_geometry', f'{fn}/{cond}', 'base_size', int(pr.loc[cond, 'base_size']), int(base))
            compare('tierB_geometry', f'{fn}/{cond}', 'tpl_offset_px', float(pr.loc[cond, 'tpl_offset_px']), round(float(np.hypot(tx - click[0], ty - click[1])), 3), tol=1e-9)
            t1 = time.time()
            det, info = search(hem, int(base), float(tx), float(ty), mpp)
            t_search = time.time() - t1
            ranked_od = det.sort_values('od', ascending=False, kind='mergesort').reset_index(drop=True)
            ranked_score = det.reset_index(drop=True)  # NMS order is score order
            nb_top = top[(top['file_name'] == fn) & (top['condition'] == cond)].sort_values('rank')
            same_xy = np.array_equal(ranked_od[['cx', 'cy']].to_numpy()[:30], nb_top[['cx', 'cy']].to_numpy())
            n_shared_xy = len(set(map(tuple, ranked_od[['cx', 'cy']].to_numpy()[:30])) & set(map(tuple, nb_top[['cx', 'cy']].to_numpy())))
            od_maxdiff = float(np.max(np.abs(ranked_od['od'].to_numpy()[:30] - nb_top['od'].to_numpy()))) if same_xy else np.nan
            compare('tierB_pipeline', f'{fn}/{cond}', 'top30 coordinates identical', True, bool(same_xy))
            compare('tierB_pipeline', f'{fn}/{cond}', 'n_detections', int(pr.loc[cond, 'n_detections']), info['n_detections'])
            full_od = ranked_od[['cx', 'cy']].to_numpy(float)
            idx_full, _ = greedy(full_od, gt_xy, radius)
            hit_full = np.array([g >= 0 and cat[g] == MITOTIC for g in idx_full])
            idx_mit_full, _ = greedy(full_od, mit_xy, radius)
            row = dict(file_name=fn, condition=cond, largest_od_tie_block=int(ranked_od['od'].value_counts().max()), base_size=int(base), tpl_cx=tx, tpl_cy=ty, round_tpl_equals_round_click=(int(round(tx)), int(round(ty))) == (int(round(click[0])), int(round(click[1]))), n_above_floor=info['n_above_floor'], n_peaks=info['n_peaks'], n_after_nms=info['n_after_nms'], n_self_hits=info['n_self_hits'], n_detections=info['n_detections'], top30_identical_to_notebook=bool(same_xy), top30_shared_coordinates=n_shared_xy, od_max_abs_diff=od_maxdiff, full_list_within_r_of_click=int((np.hypot(full_od[:, 0] - click[0], full_od[:, 1] - click[1]) <= radius).sum()), t_search_s=round(t_search, 2))
            for k in BUDGETS:
                tp_full = int(hit_full[:k].sum())
                idx_k, _ = greedy(full_od[:k], gt_xy, radius)
                tp_k = int(sum(1 for g in idx_k if g >= 0 and cat[g] == MITOTIC))
                tp_mit_full = int((idx_mit_full[:k] >= 0).sum())
                tp_mit_k = int((greedy(full_od[:k], mit_xy, radius)[0] >= 0).sum())
                row[f'tp_at_{k}_full_list_match'] = tp_full
                row[f'tp_at_{k}_topk_only_match'] = tp_k
                row[f'tp_at_{k}_mitotic_only_full'] = tp_mit_full
                row[f'tp_at_{k}_mitotic_only_topk'] = tp_mit_k
                compare('tierB_pipeline', f'{fn}/{cond}', f'tp_at_{k} full-list match vs notebook', int(pr.loc[cond, f'tp_at_{k}']), tp_full)
                compare('tierB_prefix', f'{fn}/{cond}', f'tp_at_{k} full-list vs top-K-only (repo rule)', tp_full, tp_k)
                compare('tierB_prefix', f'{fn}/{cond}', f'tp_at_{k} full-list vs top-K-only (mitotic-only rule)', tp_mit_full, tp_mit_k)
            run_rows.append(row)
            # tm_score arm and K=50, for scope and as a second-artifact check against the reference raw CSV
            full_sc = ranked_score[['cx', 'cy']].to_numpy(float)
            idx_sc, _ = greedy(full_sc, gt_xy, radius)
            hit_sc = np.array([g >= 0 and cat[g] == MITOTIC for g in idx_sc])
            for k in (10, 20, 30, 50):
                arm_rows.append(dict(file_name=fn, condition=cond, arm='tm_score', K=k, tp=int(hit_sc[:k].sum())))
                arm_rows.append(dict(file_name=fn, condition=cond, arm='chromatin_od', K=k, tp=int(hit_full[:k].sum())))
                if cond != 'default_51':
                    for arm, hits in (('tm_score', hit_sc), ('chromatin_od', hit_full)):
                        r = ref[(ref['file_name'] == fn) & (ref['condition'] == cond) & (ref['arm'] == arm) & (ref['budget'] == k)]
                        compare('tierB_vs_reference_raw', f'{fn}/{cond}/{arm}', f'tp_at_{k}', int(r['tp_at_budget'].iloc[0]), int(hits[:k].sum()))
            list_rows.append(ranked_od.assign(file_name=fn, condition=cond, rank=np.arange(len(ranked_od))))
            log(f"  [{fn}/{cond:10s}] base={base:2d} n_det={info['n_detections']} top30_same={same_xy} tp@10/20/30={[row[f'tp_at_{k}_full_list_match'] for k in BUDGETS]} ({t_search:.1f}s)")
        draw_rows.append(dict(file_name=fn, draw='summary', ann_id=int(row_j['ann_id']), gray_reason=f"solo gray draw={int(row_g['ann_id'])}", hem_reason=f"solo hem draw={int(row_h['ann_id'])}", click_readable=f"unconditional draw={int(row_u['ann_id'])}", gray_template_readable=None, hem_template_readable=None, joint_accept=f'n_retries={n_ret}'))
        del hem, gray
        log(f'  [{fn}] {time.time() - t0:.1f}s')
    t_total = time.time() - t_start
    log(f'  Tier B wall clock {t_total:.0f}s')
    draws = pd.DataFrame(draw_rows)
    runs = pd.DataFrame(run_rows)
    arms = pd.DataFrame(arm_rows)
    save(draws, 'tier_b_seed_draws')
    save(runs, 'tier_b_pipeline_runs')
    save(arms, 'tier_b_tp_by_arm')
    full_lists = pd.concat(list_rows, ignore_index=True)
    save(full_lists, 'tier_b_full_lists_chromatin_od')
    # tm_score arm, three-way, same statistics as the notebook's chromatin_od table
    arm_stats = []
    for arm in ('tm_score', 'chromatin_od'):
        sub = arms[arms['arm'] == arm]
        for k in (10, 20, 30, 50):
            s = sub[sub['K'] == k].pivot_table(index='file_name', columns='condition', values='tp')
            for c in CONDITIONS:
                arm_stats.append(dict(arm=arm, K=k, row='pooled_precision', condition=c, value=s[c].sum() / (14 * k)))
            for x, y in PAIRS:
                d = (s[x] - s[y]).astype(int).to_numpy()
                p, n = sign_flip(d)
                arm_stats.append(dict(arm=arm, K=k, row=f'{x} - {y}', condition=f'W/L/T {(d > 0).sum()}/{(d < 0).sum()}/{(d == 0).sum()}', value=float(p)))
    save(pd.DataFrame(arm_stats), 'tier_b_arm_scope_stats')
    return dict(draws=draws, runs=runs, arms=arms, t_total=t_total)


def null_share_below(deltas, level=0.10):
    """
        Under the sign-flip null conditional on |deltas|, the probability the exact p lands strictly below a level.

        deltas (array-like): one integer delta per ROI.
        level (float): threshold.

        Returns float: probability over the 2^n equally likely sign patterns.
    """
    d = np.asarray(deltas, dtype=np.int64)
    nz = np.abs(d[d != 0])
    n = len(nz)
    if n == 0:
        return 0.0
    bits = (np.arange(2 ** n)[:, None] >> np.arange(n)) & 1
    sums = np.abs(((2 * bits - 1) * nz).sum(axis=1))
    srt = np.sort(sums)
    p_of = 1.0 - np.searchsorted(srt, sums, side='left') / len(srt)
    return float(np.mean(p_of < level))


def prose_checks(a, g2, b):
    """
        Every number and ordering the notebook's markdown and printed summaries state, recomputed.

        a (dict): Tier A outputs.
        g2 (pd.DataFrame): Gate 2 against the current reference CSV.
        b (dict or None): Tier B outputs, None when skipped.

        Returns pd.DataFrame: one row per claim.
    """
    per_roi, rep, st = a['per_roi'], a['rep'], a['stats_all']
    rs = st[st['rule'] == 'repo_mixed_greedy'].set_index(['pair', 'K'])
    rows = []

    def add(cell, claim, stated, recomputed, ok):
        rows.append(dict(cell=cell, claim=claim, stated=str(stated), recomputed=str(recomputed), holds=bool(ok)))

    pooled = {(c, k): rep[rep['condition'] == c][f'tp_at_{k}'].sum() / (14 * k) for c in CONDITIONS for k in BUDGETS}
    stated_pooled = {('default_51', 10): 0.671, ('gray_bbox', 10): 0.657, ('hem_bbox', 10): 0.636, ('default_51', 20): 0.632, ('gray_bbox', 20): 0.604, ('hem_bbox', 20): 0.564, ('default_51', 30): 0.567, ('gray_bbox', 30): 0.550, ('hem_bbox', 30): 0.524}
    for key, v in stated_pooled.items():
        add(24, f'pooled precision {key[0]} K={key[1]}', v, round(pooled[key], 4), abs(round(pooled[key], 3) - v) < 1e-9)
    add(24, 'default_51 > gray_bbox > hem_bbox at every K (pooled)', 'strict order', [(round(pooled[('default_51', k)], 4), round(pooled[('gray_bbox', k)], 4), round(pooled[('hem_bbox', k)], 4)) for k in BUDGETS], all(pooled[('default_51', k)] > pooled[('gray_bbox', k)] > pooled[('hem_bbox', k)] for k in BUDGETS))
    for pair, pts, wins, losses, ties, ps in [('hem_bbox - default_51', (3.6, 6.8, 4.3), (0, 1, 1), (5, 11, 10), None, (0.063, 0.003, 0.007)), ('gray_bbox - default_51', (1.4, 2.9, 1.7), (1, 0, 3), (3, 4, 5), (10, 10, 6), (0.63, 0.13, 0.30)), ('hem_bbox - gray_bbox', (2.1, 3.9, 2.6), (0, 3, 3), (3, 9, 7), (11, 2, 4), (0.25, 0.090, 0.098))]:
        for i, k in enumerate(BUDGETS):
            r = rs.loc[(pair, k)]
            add(24, f'{pair} K={k} points of precision lost', pts[i], round(-100 * r['mean_delta_precision'], 3), abs(round(-100 * r['mean_delta_precision'], 1) - pts[i]) < 1e-9)
            add(24, f'{pair} K={k} wins/losses', f'{wins[i]}/{losses[i]}', f"{r['wins']}/{r['losses']}", (r['wins'], r['losses']) == (wins[i], losses[i]))
            if ties is not None:
                add(24, f'{pair} K={k} ties', ties[i], r['ties'], r['ties'] == ties[i])
            dec = len(str(ps[i]).split('.')[1])
            half_up = float(Decimal(str(r['exact_p'])).quantize(Decimal(1).scaleb(-dec), rounding=ROUND_HALF_UP))
            add(24, f'{pair} K={k} exact p', ps[i], r['exact_p_fraction'] + f" = {r['exact_p']:.5f}", abs(half_up - ps[i]) < 1e-9)
    hd = a['tp'][(a['tp']['rule'] == 'repo_mixed_greedy')].pivot_table(index=['file_name', 'K'], columns='condition', values='tp').reset_index()
    win_files = sorted(set(hd[(hd['hem_bbox'] > hd['default_51'])]['file_name']))
    add(24, 'hem_bbox - default_51: the one win is 013.tiff', '013.tiff', win_files, win_files == ['013.tiff'])
    adj = rs['holm_adjusted_p']
    rej = rs[rs['holm_reject_0.05']].index.tolist()
    add(24, 'Holm over 9 rejects only hem_bbox - default_51 at K=20', "[('hem_bbox - default_51', 20)]", rej, rej == [('hem_bbox - default_51', 20)])
    add(24, 'Holm step 1: 0.0029 x 9 = 0.026', 0.026, round(9 * rs.loc[('hem_bbox - default_51', 20), 'exact_p'], 5), abs(round(9 * rs.loc[('hem_bbox - default_51', 20), 'exact_p'], 3) - 0.026) < 1e-9)
    add(24, 'Holm step 2: 0.0068 x 8 = 0.054', 0.054, round(8 * rs.loc[('hem_bbox - default_51', 30), 'exact_p'], 5), abs(round(8 * rs.loc[('hem_bbox - default_51', 30), 'exact_p'], 3) - 0.054) < 1e-9)
    gd = hd.copy()
    loss_files = gd[gd['gray_bbox'] < gd['default_51']].groupby('file_name')['K'].apply(list).to_dict()
    add(24, 'gray_bbox - default_51 losses cluster on 402, 529, 548, 246', '402, 529, 548, 246', loss_files, set(loss_files) >= {'402.tiff', '529.tiff', '548.tiff', '246.tiff'})
    ceil = sorted(per_roi[(per_roi['condition'] == 'gray_bbox') & (per_roi['base_size'] == 51)]['file_name'])
    add(24, 'gray_bbox fills the 51 px window on 3 ROIs', '201, 245, 403', ceil, ceil == ['201.tiff', '245.tiff', '403.tiff'])
    ties = a['ties']
    t403 = ties[(ties['file_name'] == '403.tiff') & (ties['pair'] == 'gray_bbox - default_51')]
    add(24, '403.tiff: top-30 list changes completely, TP counts at 10/20/30 do not', '0 shared, equal TP', f"shared coords {t403['shared_coordinates'].tolist()}, ties {t403['tie'].tolist()}, shared mitoses {t403['shared_mitoses'].tolist()}", bool((t403['shared_coordinates'] == 0).all() and t403['tie'].all()))
    for fn in ('201.tiff', '245.tiff'):
        t = ties[(ties['file_name'] == fn) & (ties['pair'] == 'gray_bbox - default_51') & (ties['K'] == 30)]
        add(24, f'{fn}: template identical to default_51 and so is the list', '30/30 shared', int(t['shared_coordinates'].iloc[0]), int(t['shared_coordinates'].iloc[0]) == 30)
    worst = all(rep[rep['condition'] == c].set_index('file_name')[f'tp_at_{k}'].idxmin() == '245.tiff' or rep[rep['condition'] == c].set_index('file_name')[f'tp_at_{k}'].min() == rep[(rep['condition'] == c) & (rep['file_name'] == '245.tiff')][f'tp_at_{k}'].iloc[0] for c in CONDITIONS for k in BUDGETS)
    add(24, '245.tiff worst ROI under every condition at every K', True, worst, worst)
    g30 = rep[rep['condition'] == 'gray_bbox'].set_index('file_name')['tp_at_30']
    add(24, 'tied with 529.tiff for gray_bbox at K=30', '245;529', sorted(g30[g30 == g30.min()].index), sorted(g30[g30 == g30.min()].index) == ['245.tiff', '529.tiff'])
    r245 = rep[rep['file_name'] == '245.tiff'].set_index('condition')
    add(24, '245.tiff hem_bbox 0.20 -> 0.10 at K=20, 0.23 -> 0.17 at K=30', '0.20->0.10, 0.23->0.17', f"{r245.loc['default_51', 'tp_at_20'] / 20:.2f}->{r245.loc['hem_bbox', 'tp_at_20'] / 20:.2f}, {r245.loc['default_51', 'tp_at_30'] / 30:.2f}->{r245.loc['hem_bbox', 'tp_at_30'] / 30:.2f}", (round(r245.loc['default_51', 'tp_at_20'] / 20, 2), round(r245.loc['hem_bbox', 'tp_at_20'] / 20, 2), round(r245.loc['default_51', 'tp_at_30'] / 30, 2), round(r245.loc['hem_bbox', 'tp_at_30'] / 30, 2)) == (0.20, 0.10, 0.23, 0.17))
    geo = per_roi.pivot_table(index='file_name', columns='condition', values=['base_size', 'tpl_offset_px'])
    med = (geo[('base_size', 'default_51')].median(), geo[('base_size', 'gray_bbox')].median(), geo[('base_size', 'hem_bbox')].median())
    add(24, 'median template size 51 -> 40 -> 28 px', '51/40/28', med, tuple(int(round(m)) for m in med) == (51, 40, 28))
    add(9, 'gray_bbox median base_size printed as 40 (f-string :.0f)', 40, geo[('base_size', 'gray_bbox')].median(), True)
    add(9, 'hem_bbox median base_size printed as 28 (f-string :.0f)', 28, geo[('base_size', 'hem_bbox')].median(), True)
    moff = (geo[('tpl_offset_px', 'gray_bbox')].median(), geo[('tpl_offset_px', 'hem_bbox')].median())
    add(24, 'median offset 2.5 px gray, 3.5 px hem', '2.5/3.5', moff, abs(moff[0] - 2.5) < 0.05 and abs(moff[1] - 3.5) < 0.05)
    add(9, 'gray smaller than 51 on 11/14; hem 14/14; hem smaller than gray 13/14, equal 1/14', '11,14,13,1', (int((geo[('base_size', 'gray_bbox')] < 51).sum()), int((geo[('base_size', 'hem_bbox')] < 51).sum()), int((geo[('base_size', 'hem_bbox')] < geo[('base_size', 'gray_bbox')]).sum()), int((geo[('base_size', 'hem_bbox')] == geo[('base_size', 'gray_bbox')]).sum())), (int((geo[('base_size', 'gray_bbox')] < 51).sum()), int((geo[('base_size', 'hem_bbox')] < 51).sum()), int((geo[('base_size', 'hem_bbox')] < geo[('base_size', 'gray_bbox')]).sum()), int((geo[('base_size', 'hem_bbox')] == geo[('base_size', 'gray_bbox')]).sum())) == (11, 14, 13, 1))
    add(24, 'both Otsu conditions reproduce the reference chromatin_od numbers on 28/28 rows (re-run against the CSV now on disk)', '28/28', f"{int(g2['all_match'].sum())}/{len(g2)}", bool(g2['all_match'].all()))
    ref_nb = json.load(open(REF_NB_PATH))
    ref_text = '\n'.join(''.join(c['source']) for c in ref_nb['cells'])
    add(24, "reference notebook calls its K=10 tie count 'a structural property of the production cap'", 'quoted', f"phrase present in current reference notebook: {'structural property' in ref_text}; current text attributes K=10 ties to object convergence: {'recover the same' in ref_text}", 'structural property' in ref_text)
    ref = pd.read_csv(REF_RAW)
    add(24, 'reference largest_tie_block of 1', 1, int(ref[ref['arm'] == 'chromatin_od']['largest_tie_block'].max()), int(ref[ref['arm'] == 'chromatin_od']['largest_tie_block'].max()) == 1)
    tie_sum = ties.groupby(['pair', 'K']).agg(ties=('tie', 'sum'), mean_abs_dbase=('abs_delta_base_size', 'mean')).reset_index()
    k10 = tie_sum[tie_sum['K'] == 10].set_index('pair')
    add(24, 'Delta ties follow how different the two templates are, not the cap alone', 'ties track template difference', '; '.join(f"{p} K={k}: ties {t}, mean |d base| {m:.1f}" for p, k, t, m in tie_sum[['pair', 'K', 'ties', 'mean_abs_dbase']].itertuples(index=False)), bool(k10.loc['hem_bbox - default_51', 'ties'] < k10.loc['gray_bbox - default_51', 'ties'] and k10.loc['hem_bbox - gray_bbox', 'ties'] < k10.loc['gray_bbox - default_51', 'ties']))
    nulls = []
    for (pair, k), _ in rs.iterrows():
        x, y = pair.split(' - ')
        d = (hd[hd['K'] == k].set_index('file_name')[x] - hd[hd['K'] == k].set_index('file_name')[y]).astype(int).to_numpy()
        nulls.append(null_share_below(d))
    add(16, 'even with no real effect, about one test in ten would land below 0.10', 'about 0.9 of 9', f'expected count under the conditional sign-flip null {sum(nulls):.3f} (per test {[round(v, 4) for v in nulls]})', abs(sum(nulls) - 0.9) < 0.3)
    retr = a['nb_retries']
    add(24, 'first draw refused on 013, 245, 300, 403 (n_retries=1)', '013,245,300,403', sorted(f for f, v in retr.items() if v == 1), sorted(f for f, v in retr.items() if v == 1) == ['013.tiff', '245.tiff', '300.tiff', '403.tiff'])
    if b is not None:
        runs, draws = b['runs'], b['draws']
        add(24, 'D9 cap bound on 42/42 runs (independent search)', '42/42', int((runs['n_peaks'] == 100).sum()), int((runs['n_peaks'] == 100).sum()) == 42)
        add(24, 'seed annulus empty on 42/42 (independent full lists)', '42/42', int((runs['full_list_within_r_of_click'] == 0).sum()), int((runs['full_list_within_r_of_click'] == 0).sum()) == 42)
        add(6, 'self-hits removed per run {1: 42}', '{1: 42}', runs['n_self_hits'].value_counts().to_dict(), runs['n_self_hits'].value_counts().to_dict() == {1: 42})
        r201 = runs[(runs['file_name'] == '201.tiff') & (runs['condition'] == 'gray_bbox')].iloc[0]
        add(24, '201.tiff 0.5 px offset rounds to the click pixel in read_padded_patch', True, bool(r201['round_tpl_equals_round_click']), bool(r201['round_tpl_equals_round_click']))
        summ = draws[draws['draw'] == 'summary'].set_index('file_name')
        solo_gray = {f: int(str(v).split('=')[1]) for f, v in summ['gray_reason'].items()}
        joint = summ['ann_id'].to_dict()
        diff_gray = sorted(f for f in solo_gray if solo_gray[f] != joint[f])
        add(24, "caveat: only default_51 is scored on a click it would not have drawn alone (the production gray-only draw is not named)", 'default_51 only', f'production gray-only build_seed draw differs from the joint click on {diff_gray}', len(diff_gray) == 0)
        add(0, 'default_51 adds no constraint beyond border_filter (joint predicate == reference predicate)', True, f"click readable on every visited draw: {bool(draws[draws['draw'] != 'summary']['click_readable'].astype(bool).all())}", bool(draws[draws['draw'] != 'summary']['click_readable'].astype(bool).all()))
    return pd.DataFrame(rows)


def main():
    """
        Run Tier A, Gate 2 against the current reference CSV, and Tier B, then write the comparison ledger.

        Returns None.
    """
    t0 = time.time()
    images, ann = load_db(DB_PATH)
    a = tier_a(images, ann)
    g2 = gate2_against_current_raw(a['per_roi'])
    save(g2, 'gate2_current_reference_raw')
    log(f"  Gate 2 against the reference raw CSV as it stands now: {int(g2['all_match'].sum())}/{len(g2)}")
    b = None if '--tier-a-only' in sys.argv else tier_b(images, ann, a)
    prose = prose_checks(a, g2, b)
    save(prose, 'prose_claims')
    log(f"  prose claims: {int(prose['holds'].sum())}/{len(prose)} hold")
    ledger = pd.DataFrame(COMPARISONS)
    ledger['notebook'] = ledger['notebook'].astype(str)
    ledger['audit'] = ledger['audit'].astype(str)
    save(ledger, 'comparisons')
    summ = ledger.groupby('table')['match'].agg(compared='size', divergences=lambda m: int((~m).sum())).reset_index()
    save(summ, 'comparison_counts')
    log(summ.to_string(index=False))
    log(f'total {len(ledger)} values compared, {int((~ledger["match"]).sum())} divergences; {time.time() - t0:.0f}s')


if __name__ == '__main__':
    main()
