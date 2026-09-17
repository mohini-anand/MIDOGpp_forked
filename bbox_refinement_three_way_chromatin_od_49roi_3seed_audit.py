"""
    Independent audit of production_hematoxylin_only/bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb.

    Every number in `Research Logs/2026-09-17-bbox-refinement-three-way-chromatin-od-49roi-3seed-audit.md`
    is produced by this script. Tier A recomputes from the notebook's persisted per-run and top-30
    tables plus `databases/MIDOG++.json` and the TIFF resolution tags, in plain numpy/pandas/scipy;
    nothing in Tier A imports `midog_utils`. Tier B re-derives the click draws and 15 pipeline runs
    from pixels with a re-implementation, and runs a production-draw sensitivity through the
    production entry point (labelled as such). Run with the anaconda interpreter from the repo root:

        /Users/mohinianand/anaconda3/bin/python3 bbox_refinement_three_way_chromatin_od_49roi_3seed_audit.py

    Writes results/bbox_refinement_three_way_chromatin_od_49roi_3seed_audit_*.csv. Never writes anywhere else.
"""

from __future__ import annotations

import itertools
import json
import math
import os
import subprocess
import sys
import time
from fractions import Fraction

import cv2
import numpy as np
import pandas as pd
import tifffile
from scipy import stats
from skimage.color import rgb2hed
from skimage.measure import label, regionprops

REPO = os.path.dirname(os.path.abspath(__file__))
os.chdir(REPO)
NB = 'production_hematoxylin_only/bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb'
STEM = 'results/precision_at_k_49roi_3seed_chromatin_bbox3way'
REF14 = 'results/precision_at_k_14roi_prodseed_chromatin_bbox3way'
DB = 'databases/MIDOG++.json'
OUT = 'results/bbox_refinement_three_way_chromatin_od_49roi_3seed_audit'
SUBSET_DIRS = {'original_14': 'images/extra_valid', 'testing_35': 'images/extra_valid/testing_set'}
CONDS = ('default_51', 'gray_bbox', 'hem_bbox')
PAIRS = (('gray_bbox', 'default_51'), ('hem_bbox', 'default_51'), ('hem_bbox', 'gray_bbox'))
KS = (10, 20, 30)
SEEDS = (0, 1, 2)
SUBSETS = ('all_49', 'original_14', 'testing_35')
RADIUS_UM = 7.5
T0 = time.time()

COMPARISONS = []   # (section, quantity, n_compared, n_divergent, max_abs_diff, independence)
PROSE = []         # (cell, claim, claimed, recomputed, match)


def save(frame, name):
    """
        Write one audit table.

        frame (pd.DataFrame): the table.
        name (str): suffix after the audit stem.

        Returns str: the path written.
    """
    path = f'{OUT}_{name}.csv'
    frame.to_csv(path, index=False)
    print(f'  -> {path} ({len(frame)} rows)')
    return path


def compare(section, quantity, mine, theirs, independence, tol=0.0):
    """
        Record a value-by-value comparison of two aligned arrays.

        section (str): audit section label.
        quantity (str): what is being compared.
        mine (array-like): recomputed values.
        theirs (array-like): the notebook's persisted values.
        independence (str): "independent" or "consistency".
        tol (float): absolute tolerance for numeric values.

        Returns int: number of divergent values.
    """
    a, b = np.asarray(mine), np.asarray(theirs)
    assert a.shape == b.shape, f'{section}/{quantity}: shape {a.shape} vs {b.shape}'
    if a.dtype.kind in 'fiub' and b.dtype.kind in 'fiub':
        a, b = a.astype(float), b.astype(float)
        both_nan = np.isnan(a) & np.isnan(b)
        diff = np.where(both_nan, 0.0, np.abs(a - b))
        diff = np.where(np.isnan(diff), np.inf, diff)
        bad = int((diff > tol).sum())
        mx = float(np.max(diff)) if diff.size else 0.0
    else:
        bad = int((a.astype(str) != b.astype(str)).sum())
        mx = float('nan')
    COMPARISONS.append(dict(section=section, quantity=quantity, n_compared=int(a.size), n_divergent=bad, max_abs_diff=mx, tolerance=tol, independence=independence))
    return bad


def prose(cell, claim, claimed, recomputed, tol=0.0):
    """
        Record one prose or printed claim against its recomputation.

        cell (int): notebook cell index holding the claim.
        claim (str): what the claim says.
        claimed (object): the value the notebook states.
        recomputed (object): the value this script derives.
        tol (float): tolerance for numeric claims.

        Returns bool: whether the claim matches.
    """
    if claimed == 'n/a':
        ok = True
    elif isinstance(claimed, (int, float, np.integer, np.floating)) and isinstance(recomputed, (int, float, np.integer, np.floating)):
        ok = abs(float(claimed) - float(recomputed)) <= tol
    else:
        ok = str(claimed) == str(recomputed)
    PROSE.append(dict(cell=cell, claim=claim, claimed=str(claimed), recomputed=str(recomputed), match=bool(ok)))
    return ok


# ----------------------------------------------------------------------------------------------
# 1. Execution coherence
# ----------------------------------------------------------------------------------------------
def execution_gate():
    """
        Execution-count contiguity, errors, unrun cells and kernel-timestamp monotonicity.

        Returns pd.DataFrame: one row per code cell.
    """
    nb = json.load(open(NB))
    rows = []
    for i, c in enumerate(nb['cells']):
        if c['cell_type'] != 'code':
            continue
        meta = c.get('metadata', {}).get('execution', {})
        rows.append(dict(cell=i, execution_count=c.get('execution_count'), n_outputs=len(c.get('outputs', [])), has_error=any(o['output_type'] == 'error' for o in c.get('outputs', [])), n_png=sum('image/png' in o.get('data', {}) for o in c.get('outputs', [])), started=meta.get('iopub.execute_input'), finished=meta.get('shell.execute_reply')))
    g = pd.DataFrame(rows)
    counts = g['execution_count'].tolist()
    contiguous = counts == list(range(1, len(counts) + 1))
    ts = pd.to_datetime(g['started'])
    monotonic = bool(ts.is_monotonic_increasing)
    print(f'[exec] {len(g)} code cells, execution_count contiguous 1..{len(g)}: {contiguous}; errors {int(g.has_error.sum())}; unrun {int(g.execution_count.isna().sum())}; last cell run: {g.execution_count.iloc[-1] is not None}; kernel timestamps monotonic: {monotonic}; {ts.iloc[0]} -> {pd.to_datetime(g.finished.iloc[-1])}')
    g['contiguous_all'] = contiguous
    g['timestamps_monotonic_all'] = monotonic
    return g


# ----------------------------------------------------------------------------------------------
# 2. Provenance
# ----------------------------------------------------------------------------------------------
def git(*args):
    """
        Run a git command in the repo and return stdout.

        args (str): git arguments.

        Returns str: stripped stdout.
    """
    return subprocess.run(['git', *args], capture_output=True, text=True, cwd=REPO).stdout.strip()


def provenance():
    """
        mtime, last commit and working-tree status of the notebook, its imported modules and its artifacts.

        Returns pd.DataFrame: one row per file.
    """
    files = [NB, 'BBOX3WAY_49ROI_3SEED_PROMPT.md'] + [f'midog_utils/{m}.py' for m in ('channels', 'chromatin', 'dataset', 'evaluate', 'find_and_suppress', 'invariants', 'nms', 'production', 'seed_selection', 'template_match')] + [f'{STEM}_{n}.csv' for n in ('per_run', 'per_roi', 'summary', 'by_domain', 'by_seed', 'delta_per_roi', 'delta_stats', 'bootstrap_ci', 'contested_sensitivity', 'top30', 'verification')] + [f'{REF14}_{n}.csv' for n in ('per_roi', 'top30', 'delta_per_roi', 'delta_stats')] + [DB]
    rows = []
    for f in files:
        st = os.stat(f)
        rows.append(dict(path=f, mtime=pd.Timestamp(st.st_mtime, unit='s', tz='UTC').tz_convert('America/New_York').strftime('%Y-%m-%d %H:%M:%S'), last_commit=git('log', '-1', '--format=%h %ad %an', '--date=iso', '--', f) or 'UNTRACKED/never committed', status=git('status', '--porcelain', '--', f) or 'clean'))
    p = pd.DataFrame(rows)
    mods = p[p.path.str.startswith('midog_utils/')]
    arts = p[p.path.str.startswith(STEM)]
    print(f'[prov] newest module mtime {mods.mtime.max()} (status: {sorted(set(mods.status))}); artifacts {arts.mtime.min()} .. {arts.mtime.max()} (status: {sorted(set(arts.status))}); notebook {p.loc[p.path == NB, "mtime"].iat[0]} ({p.loc[p.path == NB, "status"].iat[0]})')
    print(f"[prov] notebook author signal: {git('log', '-1', '--format=%an %ar', '--', NB) or 'no commit (untracked)'}; DECISIONS.md status: {git('status', '--porcelain', '--', 'DECISIONS.md')}")
    return p


# ----------------------------------------------------------------------------------------------
# 3. Data loading, independent of midog_utils
# ----------------------------------------------------------------------------------------------
def load_db():
    """
        Parse MIDOG++.json directly.

        Returns tuple[pd.DataFrame, pd.DataFrame]: images (id, file_name, width, height, tumor_type) and annotations (ann_id, file_name, cx, cy, category, votes, unanimous).
    """
    raw = json.load(open(DB))
    alias = {'canine lymphoma': 'canine lymphosarcoma'}
    imgs = pd.DataFrame([dict(image_id=im['id'], file_name=im['file_name'], width=im['width'], height=im['height'], tumor_type=alias.get(im['tumor_type'], im['tumor_type'])) for im in raw['images']])
    name = dict(zip(imgs.image_id, imgs.file_name))
    anns = pd.DataFrame([dict(ann_id=a['id'], image_id=a['image_id'], file_name=name[a['image_id']], cx=(a['bbox'][0] + a['bbox'][2]) / 2.0, cy=(a['bbox'][1] + a['bbox'][3]) / 2.0, category=a['category_id'], n_votes=len(a.get('labels', [])), n_mitotic_votes=sum(v == 1 for v in a.get('labels', [])), unanimous=len(set(a.get('labels', []))) == 1 if a.get('labels') else False) for a in raw['annotations']])
    return imgs, anns


def tiff_mpp(path):
    """
        Microns per pixel from TIFF XResolution/ResolutionUnit, read directly.

        path (str): TIFF path.

        Returns float: microns per pixel.
    """
    with tifffile.TiffFile(path) as tf:
        page = tf.pages[0]
        num, den = page.tags['XResolution'].value
        unit = int(page.tags['ResolutionUnit'].value)
    return {2: 25400.0, 3: 10000.0}[unit] / (num / den)


def greedy_match(det_xy, gt_xy, radius):
    """
        Rank-order one-to-one matching: each detection claims the nearest unclaimed GT within radius.

        det_xy (np.ndarray): (N, 2) detections, best first.
        gt_xy (np.ndarray): (M, 2) ground truth.
        radius (float): match radius, px.

        Returns np.ndarray: GT index claimed per detection, -1 for none.
    """
    claimed = np.zeros(len(gt_xy), dtype=bool)
    out = np.full(len(det_xy), -1, dtype=int)
    for i, (x, y) in enumerate(det_xy):
        d = np.hypot(gt_xy[:, 0] - x, gt_xy[:, 1] - y)
        ok = np.flatnonzero((d <= radius) & ~claimed)
        if len(ok):
            g = ok[np.argmin(d[ok])]
            out[i] = g
            claimed[g] = True
    return out


# ----------------------------------------------------------------------------------------------
# 4. Composition gate
# ----------------------------------------------------------------------------------------------
def composition(T, imgs):
    """
        Row counts, key uniqueness and strata against the design.

        T (dict): the notebook's artifact tables.
        imgs (pd.DataFrame): DB images.

        Returns pd.DataFrame: one row per check.
    """
    rows = []
    exp = dict(per_run=49 * 3 * 3, per_roi=49 * 3, summary=3 * 3 * 3, by_domain=7 * 3 * 3, by_seed=3 * 3 * 3 * 3, delta_per_roi=49 * 3 * 3, delta_stats=(3 + 7) * 9, bootstrap_ci=3 * (3 + 3) * 3, contested_sensitivity=3 * 3 * 3 * 3, top30=441 * 30, verification=441 * 5 + 147 + 49 * 3 + 42 + 42 + 17)
    for n, e in exp.items():
        rows.append(dict(check=f'rows:{n}', expected=e, observed=len(T[n]), passed=len(T[n]) == e))
    keys = dict(per_run=['file_name', 'seed_index', 'condition'], per_roi=['file_name', 'condition'], delta_per_roi=['file_name', 'pair', 'K'], top30=['file_name', 'seed_index', 'condition', 'rank'], summary=['subset', 'condition', 'K'], delta_stats=['group', 'pair', 'K'])
    for n, k in keys.items():
        d = int(T[n].duplicated(k).sum())
        rows.append(dict(check=f'duplicate_keys:{n}', expected=0, observed=d, passed=d == 0))
    pr = T['per_run']
    rows.append(dict(check='rois', expected=49, observed=pr.file_name.nunique(), passed=pr.file_name.nunique() == 49))
    rows.append(dict(check='rois_original_14', expected=14, observed=pr[pr.subset == 'original_14'].file_name.nunique(), passed=pr[pr.subset == 'original_14'].file_name.nunique() == 14))
    rows.append(dict(check='rois_testing_35', expected=35, observed=pr[pr.subset == 'testing_35'].file_name.nunique(), passed=pr[pr.subset == 'testing_35'].file_name.nunique() == 35))
    # the directories, listed independently
    on_disk = {s: sorted(f for f in os.listdir(d) if f.endswith('.tiff')) for s, d in SUBSET_DIRS.items()}
    for s, fs in on_disk.items():
        rows.append(dict(check=f'files_on_disk_match:{s}', expected=len(fs), observed=len(set(fs) & set(pr[pr.subset == s].file_name)), passed=set(fs) == set(pr[pr.subset == s].file_name)))
    firsts = pr.groupby('file_name').first()
    dom_db = pd.Series(firsts.index, index=firsts.index).map(imgs.set_index('file_name').tumor_type)
    rows.append(dict(check='domain_column_matches_db', expected=49, observed=int((dom_db == firsts.domain).sum()), passed=bool((dom_db == firsts.domain).all())))
    per_dom = dom_db.value_counts()
    rows.append(dict(check='rois_per_domain_all_7', expected='7x7', observed=';'.join(f'{d}={n}' for d, n in per_dom.items()), passed=len(per_dom) == 7 and bool((per_dom == 7).all())))
    annotated = imgs[~imgs.image_id.isin(range(151, 201))]
    db_dom = annotated.tumor_type.value_counts()
    rows.append(dict(check='db_images_per_tumor_type (context)', expected='n/a', observed=';'.join(f'{d}={n}' for d, n in db_dom.items()), passed=True))
    cells = pr.groupby(['file_name']).size()
    rows.append(dict(check='runs_per_roi_is_9', expected=9, observed=f'{cells.min()}..{cells.max()}', passed=bool((cells == 9).all())))
    t = T['top30'].groupby(['file_name', 'seed_index', 'condition'])['rank']
    rows.append(dict(check='top30_ranks_0_to_29_every_run', expected=441, observed=int((t.apply(lambda r: sorted(r) == list(range(30)))).sum()), passed=bool(t.apply(lambda r: sorted(r) == list(range(30))).all())))
    c = pd.DataFrame(rows)
    print(f'[comp] {int(c.passed.sum())}/{len(c)} composition checks pass; failing: {c.loc[~c.passed, "check"].tolist()}')
    return c


# ----------------------------------------------------------------------------------------------
# 5. Independent re-scoring from top30 + DB + TIFF tags
# ----------------------------------------------------------------------------------------------
def rescore(T, imgs, anns):
    """
        Re-derive buckets, matched ids and TP@K for all 441 runs; check seeds, geometry and mpp.

        T (dict): artifact tables.
        imgs (pd.DataFrame): DB images.
        anns (pd.DataFrame): DB annotations.

        Returns tuple[pd.DataFrame, pd.DataFrame]: per-run recomputation (with my tp_at_k) and per-detection recomputed top-30 rows.
    """
    pr, top = T['per_run'], T['top30'].sort_values(['file_name', 'seed_index', 'condition', 'rank']).reset_index(drop=True)
    paths = {fn: f'{SUBSET_DIRS[s]}/{fn}' for s, fn in pr[['subset', 'file_name']].drop_duplicates().itertuples(index=False)}
    mpp = {fn: tiff_mpp(p) for fn, p in paths.items()}
    wh = imgs.set_index('file_name')[['width', 'height']]
    by_img = {fn: g.reset_index(drop=True) for fn, g in anns.groupby('file_name')}
    run_rows, det_rows = [], []
    for (fn, s, cond), g in top.groupby(['file_name', 'seed_index', 'condition'], sort=False):
        run = pr[(pr.file_name == fn) & (pr.seed_index == s) & (pr.condition == cond)].iloc[0]
        gt = by_img[fn]
        gt_eval = gt[gt.ann_id != run.seed_ann_id].reset_index(drop=True)
        r = RADIUS_UM / mpp[fn]
        m = greedy_match(g[['cx', 'cy']].to_numpy(float), gt_eval[['cx', 'cy']].to_numpy(float), r)
        cat = np.where(m >= 0, gt_eval.category.to_numpy()[np.maximum(m, 0)], 0)
        bucket = np.where(m < 0, 'non_human_findings', np.where(cat == 1, 'human_correct_label', 'human_rejected_label'))
        mid = np.where(m >= 0, gt_eval.ann_id.to_numpy()[np.maximum(m, 0)], -1)
        hit = bucket == 'human_correct_label'
        contested = np.array([(i >= 0) and (not bool(gt_eval.unanimous.iat[i])) for i in m]) & hit
        mito = gt_eval[gt_eval.category == 1]
        seed_row = gt[gt.ann_id == run.seed_ann_id].iloc[0]
        # nearest mitotic GT distance per detection (for the 211.tiff question)
        dmin = np.array([np.hypot(mito.cx - x, mito.cy - y).min() for x, y in g[['cx', 'cy']].to_numpy(float)])
        det_rows.append(pd.DataFrame(dict(file_name=fn, seed_index=s, condition=cond, rank=g['rank'].to_numpy(), my_bucket=bucket, my_matched_ann_id=mid, nb_bucket=g.bucket.to_numpy(), nb_matched_ann_id=g.matched_ann_id.to_numpy(), my_contested_hit=contested, d_nearest_mitotic_px=dmin, od=g.od.to_numpy())))
        W, H = wh.loc[fn, 'width'], wh.loc[fn, 'height']
        ix, iy = int(round(float(run.tpl_cx))), int(round(float(run.tpl_cy)))
        cx_r, cy_r = int(np.rint(seed_row.cx)), int(np.rint(seed_row.cy))
        run_rows.append(dict(file_name=fn, seed_index=s, condition=cond, my_mpp=mpp[fn], nb_mpp=run.mpp, my_radius=r, nb_radius=run.match_radius_px, my_tp_10=int(hit[:10].sum()), my_tp_20=int(hit[:20].sum()), my_tp_30=int(hit[:30].sum()), nb_tp_10=run.tp_at_10, nb_tp_20=run.tp_at_20, nb_tp_30=run.tp_at_30, my_contested_10=int(contested[:10].sum()), my_contested_20=int(contested[:20].sum()), my_contested_30=int(contested[:30].sum()), my_n_gt_mitotic=int((gt_eval.category == 1).sum()), nb_n_gt_mitotic=run.n_gt_mitotic, my_nearest_other=round(float(np.hypot(gt_eval.cx - seed_row.cx, gt_eval.cy - seed_row.cy).min()), 2), nb_nearest_other=run.nearest_other_ann_px, seed_is_mitotic=seed_row.category == 1, seed_is_unanimous=bool(seed_row.n_votes > 0 and seed_row.n_mitotic_votes == seed_row.n_votes), seed_click_matches=bool(np.isclose(seed_row.cx, run.click_cx) and np.isclose(seed_row.cy, run.click_cy)), seed_border_ok=bool(cx_r >= 36 and cx_r <= W - 1 - 36 and cy_r >= 36 and cy_r <= H - 1 - 36), tpl_patch_readable=bool(ix - 36 >= 0 and iy - 36 >= 0 and ix + 36 < W and iy + 36 < H), my_offset=round(float(np.hypot(run.tpl_cx - run.click_cx, run.tpl_cy - run.click_cy)), 3), nb_offset=run.tpl_offset_px, od_descending=bool(np.all(np.diff(g.od.to_numpy()) <= 0)), od_nan=int(np.isnan(g.od.to_numpy()).sum()), base_odd=int(run.base_size) % 2 == 1, base_le_51=int(run.base_size) <= 51))
    R = pd.DataFrame(run_rows)
    D = pd.concat(det_rows, ignore_index=True)
    n = 0
    n += compare('rescore', 'mpp from TIFF tags', R.my_mpp, R.nb_mpp, 'independent', 1e-9)
    n += compare('rescore', 'match radius px', R.my_radius, R.nb_radius, 'independent', 1e-9)
    for k in KS:
        n += compare('rescore', f'tp_at_{k} (441 runs)', R[f'my_tp_{k}'], R[f'nb_tp_{k}'], 'independent')
    n += compare('rescore', 'bucket (13,230 top-30 rows)', D.my_bucket, D.nb_bucket, 'independent')
    n += compare('rescore', 'matched_ann_id (13,230 top-30 rows)', D.my_matched_ann_id, D.nb_matched_ann_id, 'independent')
    n += compare('rescore', 'n_gt_mitotic', R.my_n_gt_mitotic, R.nb_n_gt_mitotic, 'independent')
    n += compare('rescore', 'nearest_other_ann_px', R.my_nearest_other, R.nb_nearest_other, 'independent', 0.011)
    n += compare('rescore', 'tpl_offset_px', R.my_offset, R.nb_offset, 'independent', 0.0011)
    flags = ['seed_is_mitotic', 'seed_is_unanimous', 'seed_click_matches', 'seed_border_ok', 'tpl_patch_readable', 'od_descending', 'base_odd', 'base_le_51']
    for f in flags:
        n += compare('rescore', f'{f} (all True)', R[f].astype(bool), np.ones(len(R), dtype=bool), 'independent')
    n += compare('rescore', 'od NaN count (all 0)', R.od_nan, np.zeros(len(R)), 'independent')
    d51 = pr[pr.condition == 'default_51']
    n += compare('rescore', 'default_51 base 51 and template == click', ((d51.base_size == 51) & (d51.tpl_cx == d51.click_cx) & (d51.tpl_cy == d51.click_cy)).to_numpy(), np.ones(len(d51), dtype=bool), 'independent')
    per_click = pr.groupby(['file_name', 'seed_index']).agg(a=('seed_ann_id', 'nunique'), g=('n_gt_mitotic', 'nunique'))
    n += compare('rescore', 'Gate C: one ann and one n_gt per (ROI, click)', ((per_click.a == 1) & (per_click.g == 1)).to_numpy(), np.ones(len(per_click), dtype=bool), 'independent')
    distinct = pr.groupby('file_name').seed_ann_id.nunique()
    n += compare('rescore', 'Gate C: 3 distinct clicks per ROI', distinct.to_numpy(), np.full(49, 3), 'independent')
    print(f'[rescore] {n} divergences over the rescoring comparisons')
    return R, D


# ----------------------------------------------------------------------------------------------
# 6. Tables B, C, per-domain, per-seed, delta_per_roi
# ----------------------------------------------------------------------------------------------
def subset_rows(frame, subset):
    """
        Rows of a frame in one analysis subset.

        frame (pd.DataFrame): frame with a subset column.
        subset (str): all_49, original_14 or testing_35.

        Returns pd.DataFrame: the rows.
    """
    return frame if subset == 'all_49' else frame[frame.subset == subset]


def tables(T, R):
    """
        Recompute per_roi, summary, by_domain, by_seed and delta_per_roi from my own TP counts.

        T (dict): artifact tables.
        R (pd.DataFrame): my per-run rescoring.

        Returns tuple[pd.DataFrame, pd.DataFrame]: my per-run frame (with subset/domain) and my delta_per_roi.
    """
    pr = T['per_run'][['subset', 'domain', 'file_name', 'seed_index', 'condition']].merge(R[['file_name', 'seed_index', 'condition', 'my_tp_10', 'my_tp_20', 'my_tp_30', 'my_contested_10', 'my_contested_20', 'my_contested_30', 'my_n_gt_mitotic']], on=['file_name', 'seed_index', 'condition'], how='inner', validate='1:1')
    assert len(pr) == 441
    pr = pr.rename(columns={f'my_tp_{k}': f'tp{k}' for k in KS})
    # per_roi
    roi = pr.groupby(['file_name', 'condition']).agg(**{f'tp{k}': (f'tp{k}', 'sum') for k in KS}, n=('seed_index', 'nunique')).reset_index()
    nb = T['per_roi'].merge(roi, on=['file_name', 'condition'], validate='1:1')
    for k in KS:
        compare('tables', f'per_roi tp_at_{k}', nb[f'tp{k}'], nb[f'tp_at_{k}'], 'independent')
        compare('tables', f'per_roi precision_at_{k}', nb[f'tp{k}'] / (3 * k), nb[f'precision_at_{k}'], 'independent', 1e-12)
    # summary
    srows = []
    for sub in SUBSETS:
        runs = subset_rows(pr, sub)
        rois = runs.groupby(['file_name', 'condition']).agg(**{f'tp{k}': (f'tp{k}', 'sum') for k in KS}).reset_index()
        for c in CONDS:
            rc, rr = runs[runs.condition == c], rois[rois.condition == c]
            for k in KS:
                p_roi = rr[f'tp{k}'] / (3 * k)
                p_click = rc[f'tp{k}'] / k
                srows.append(dict(subset=sub, condition=c, K=k, my_tp_sum=int(rc[f'tp{k}'].sum()), my_pooled=round(rc[f'tp{k}'].sum() / (k * len(rc)), 4), my_worst_roi=round(float(p_roi.min()), 4), my_n_roi_at_worst=int((p_roi == p_roi.min()).sum()), my_worst_roi_files=';'.join(sorted(rr.loc[p_roi == p_roi.min(), 'file_name'])), my_best_roi=round(float(p_roi.max()), 4), my_worst_click=round(float(p_click.min()), 4), my_n_clicks_at_worst=int((p_click == p_click.min()).sum())))
    S = pd.DataFrame(srows).merge(T['summary'], on=['subset', 'condition', 'K'], validate='1:1')
    compare('tables', 'summary tp_sum', S.my_tp_sum, S.tp_sum, 'independent')
    compare('tables', 'summary precision_pooled', S.my_pooled, S.precision_pooled, 'independent', 1e-9)
    compare('tables', 'summary worst_roi_precision', S.my_worst_roi, S.worst_roi_precision, 'independent', 1e-9)
    compare('tables', 'summary n_roi_at_worst', S.my_n_roi_at_worst, S.n_roi_at_worst, 'independent')
    compare('tables', 'summary worst_roi_files', S.my_worst_roi_files, S.worst_roi_files.map(lambda v: ';'.join(sorted(str(v).split(';')))), 'independent')
    compare('tables', 'summary best_roi_precision', S.my_best_roi, S.best_roi_precision, 'independent', 1e-9)
    compare('tables', 'summary worst_click_precision', S.my_worst_click, S.worst_click_precision, 'independent', 1e-9)
    compare('tables', 'summary n_clicks_at_worst', S.my_n_clicks_at_worst, S.n_clicks_at_worst, 'independent')
    # by_domain
    drows = []
    for (d, c), g in pr.groupby(['domain', 'condition']):
        rr = g.groupby('file_name').agg(**{f'tp{k}': (f'tp{k}', 'sum') for k in KS})
        for k in KS:
            drows.append(dict(domain=d, condition=c, K=k, my_tp_sum=int(g[f'tp{k}'].sum()), my_pooled=round(g[f'tp{k}'].sum() / (k * len(g)), 4), my_worst_roi=round(float((rr[f'tp{k}'] / (3 * k)).min()), 4), my_worst_click=round(float((g[f'tp{k}'] / k).min()), 4)))
    BD = pd.DataFrame(drows).merge(T['by_domain'], on=['domain', 'condition', 'K'], validate='1:1')
    compare('tables', 'by_domain tp_sum', BD.my_tp_sum, BD.tp_sum, 'independent')
    compare('tables', 'by_domain precision_pooled', BD.my_pooled, BD.precision_pooled, 'independent', 1e-9)
    compare('tables', 'by_domain worst_roi_precision', BD.my_worst_roi, BD.worst_roi_precision, 'independent', 1e-9)
    compare('tables', 'by_domain worst_click_precision', BD.my_worst_click, BD.worst_click_precision, 'independent', 1e-9)
    # by_seed
    brows = []
    for sub in SUBSETS:
        for (s, c), g in subset_rows(pr, sub).groupby(['seed_index', 'condition']):
            for k in KS:
                brows.append(dict(subset=sub, seed_index=s, condition=c, K=k, my_tp_sum=int(g[f'tp{k}'].sum()), my_pooled=round(g[f'tp{k}'].sum() / (k * len(g)), 4)))
    BS = pd.DataFrame(brows).merge(T['by_seed'], on=['subset', 'seed_index', 'condition', 'K'], validate='1:1')
    compare('tables', 'by_seed tp_sum', BS.my_tp_sum, BS.tp_sum, 'independent')
    compare('tables', 'by_seed precision_pooled', BS.my_pooled, BS.precision_pooled, 'independent', 1e-9)
    # ordering claim: default > gray > hem within all_49, per click index, every K (27 cells = 9 x 3 conditions)
    order_ok = []
    for (s, k), g in BS[BS.subset == 'all_49'].groupby(['seed_index', 'K']):
        v = g.set_index('condition').my_tp_sum
        order_ok.append(bool(v['default_51'] > v['gray_bbox'] > v['hem_bbox']))
    prose(27, 'order default_51 > gray_bbox > hem_bbox holds in all_49 for each click index at every K (9 seed x K cells)', '9/9', f'{sum(order_ok)}/9')
    pooled_order = [bool(S[(S.subset == 'all_49') & (S.K == k)].set_index('condition').my_tp_sum.pipe(lambda v: v['default_51'] > v['gray_bbox'] > v['hem_bbox'])) for k in KS]
    prose(27, 'order default_51 > gray_bbox > hem_bbox holds at every K (all_49 pooled)', '3/3', f'{sum(pooled_order)}/3')
    # delta_per_roi
    wide = pr.pivot_table(index=['file_name', 'seed_index'], columns='condition', values=[f'tp{k}' for k in KS])
    drows = []
    meta = pr.groupby('file_name')[['subset', 'domain']].first()
    for a, b in PAIRS:
        for k in KS:
            d = (wide[(f'tp{k}', a)] - wide[(f'tp{k}', b)]).unstack('seed_index')
            for fn, row in d.iterrows():
                drows.append(dict(subset=meta.loc[fn, 'subset'], domain=meta.loc[fn, 'domain'], file_name=fn, pair=f'{a} - {b}', K=k, s0=int(row[0]), s1=int(row[1]), s2=int(row[2]), dsum=int(row.sum())))
    DP = pd.DataFrame(drows)
    m = DP.merge(T['delta_per_roi'], on=['file_name', 'pair', 'K'], validate='1:1')
    for s in SEEDS:
        compare('tables', f'delta_per_roi delta_tp_s{s}', m[f's{s}'], m[f'delta_tp_s{s}'], 'independent')
    compare('tables', 'delta_per_roi delta_tp_sum', m.dsum, m.delta_tp_sum, 'independent')
    compare('tables', 'delta_per_roi delta_precision_sum', m.dsum / m.K, m.delta_precision_sum, 'independent', 1e-12)
    return pr, DP


# ----------------------------------------------------------------------------------------------
# 7. Statistics: exact sign-flip, Holm, bootstrap, t intervals, majority
# ----------------------------------------------------------------------------------------------
def signflip_exact(deltas):
    """
        Exact two-sided sign-flip p on integer deltas, by polynomial convolution of (x^-|d| + x^|d|).

        deltas (array-like): integer ROI-level deltas.

        Returns tuple[float, int]: p and the number of nonzero deltas.
    """
    d = np.abs(np.asarray([int(v) for v in deltas if int(v) != 0], dtype=np.int64))
    if len(d) == 0:
        return 1.0, 0
    dist = np.array([1], dtype=np.int64)
    for v in d:
        kern = np.zeros(2 * v + 1, dtype=np.int64)
        kern[0] = kern[-1] = 1
        dist = np.convolve(dist, kern)
    total = int(d.sum())
    sums = np.arange(-total, total + 1)
    obs = abs(int(np.asarray(deltas, dtype=np.int64).sum()))
    return float(Fraction(int(dist[np.abs(sums) >= obs].sum()), 2 ** len(d))), len(d)


def signflip_mc(values, n_draws, seed):
    """
        Monte Carlo two-sided sign-flip p on real-valued ROI deltas.

        values (array-like): ROI-level deltas.
        n_draws (int): number of sign patterns.
        seed (int): RNG seed.

        Returns float: (1 + extreme) / (1 + n_draws).
    """
    v = np.asarray(values, dtype=float)
    v = v[v != 0]
    if len(v) == 0:
        return 1.0
    rng = np.random.default_rng(seed)
    obs = abs(v.sum())
    ext = 0
    for start in range(0, n_draws, 50_000):
        n = min(50_000, n_draws - start)
        sg = rng.choice([-1.0, 1.0], size=(n, len(v)))
        ext += int((np.abs(sg @ v) >= obs - 1e-9 * max(1.0, obs)).sum())
    return (1 + ext) / (1 + n_draws)


def holm(p):
    """
        Holm step-down adjusted p-values.

        p (array-like): raw p-values.

        Returns np.ndarray: adjusted p, input order.
    """
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p, kind='mergesort')
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * p[i]))
        adj[i] = running
    return adj


def inference(T, DP, zero_rois=()):
    """
        Recompute delta_stats and D5-bar columns independently; add t intervals, bootstrap-seed sensitivity and tie-excluded majorities.

        T (dict): artifact tables.
        DP (pd.DataFrame): my delta_per_roi.
        zero_rois (iterable): ROIs with 0 TP on every run, which tie by construction.

        Returns pd.DataFrame: one row per (group, pair, K).
    """
    ds_nb = T['delta_stats']
    groups = [('subset', s) for s in SUBSETS] + [('domain', d) for d in sorted(DP.domain.unique())]
    rows = []
    for scope, grp in groups:
        f = subset_rows(DP, grp) if scope == 'subset' else DP[DP.domain == grp]
        for a, b in PAIRS:
            for k in KS:
                g = f[(f.pair == f'{a} - {b}') & (f.K == k)].sort_values('file_name')
                d = g.dsum.to_numpy()
                p, nz = signflip_exact(d)
                G = len(d)
                pooled = d.sum() / (G * 3 * k)
                mean, sd = d.mean(), d.std(ddof=1)
                tcrit = stats.t.ppf(0.975, G - 1)
                tlo, thi = (mean - tcrit * sd / math.sqrt(G)) / (3 * k), (mean + tcrit * sd / math.sqrt(G)) / (3 * k)
                wins, losses, ties = int((d > 0).sum()), int((d < 0).sum()), int((d == 0).sum())
                same = wins if pooled > 0 else (losses if pooled < 0 else 0)
                n_zero = int(g.file_name.isin(list(zero_rois)).sum())
                rows.append(dict(n_structural_zero_rois=n_zero, majority_excl_structural_zero=bool(same > (G - n_zero) / 2), scope=scope, group=grp, pair=f'{a} - {b}', K=k, G=G, pooled_delta=pooled, wins=wins, losses=losses, ties=ties, n_nonzero=nz, exact_p=p, min_attainable_p=2.0 ** (1 - nz) if nz else 1.0, t_ci_low=tlo, t_ci_high=thi, t_excludes_0=bool(tlo > 0 or thi < 0), n_same=same, majority_all=bool(same > G / 2), majority_nonzero=bool(same > nz / 2)))
    I = pd.DataFrame(rows)
    I['holm_p'] = np.nan
    for fam in ('all_49', 'testing_35'):
        msk = (I.scope == 'subset') & (I.group == fam)
        I.loc[msk, 'holm_p'] = holm(I.loc[msk, 'exact_p'])
    m = I.merge(ds_nb, on=['scope', 'group', 'pair', 'K'], validate='1:1', suffixes=('', '_nb'))
    compare('inference', 'pooled_delta_precision (90 rows)', m.pooled_delta, m.pooled_delta_precision, 'independent', 1e-12)
    compare('inference', 'wins/losses/ties (270 values)', m[['wins', 'losses', 'ties']].to_numpy(), m[['wins_nb', 'losses_nb', 'ties_nb']].to_numpy(), 'independent')
    compare('inference', 'n_nonzero', m.n_nonzero, m.n_nonzero_nb, 'independent')
    compare('inference', 'exact_p (90 rows)', m.exact_p, m.exact_p_nb, 'independent', 1e-12)
    compare('inference', 'min_attainable_p', m.min_attainable_p, m.min_attainable_p_nb, 'independent', 1e-15)
    sub = m[m.scope == 'subset']
    fam = sub[sub.group.isin(['all_49', 'testing_35'])]
    compare('inference', 'holm_p (18 family rows)', fam.holm_p, fam.holm_p_nb, 'independent', 1e-12)
    compare('inference', 'n_same_direction', m.n_same, m.n_same_direction, 'independent')
    compare('inference', 'majority_of_all_rois (27 subset rows)', sub.majority_all.astype(bool), sub.majority_of_all_rois.astype(bool), 'independent')
    compare('inference', 'majority_of_nonzero_rois (27 subset rows)', sub.majority_nonzero.astype(bool), sub.majority_of_nonzero_rois.astype(bool), 'independent')
    compare('inference', 'mc_p within 5 SE + 2/N of exact (90 rows)', (np.abs(ds_nb.mc_p - ds_nb.exact_p) <= 5 * np.sqrt(ds_nb.exact_p * (1 - ds_nb.exact_p) / 200_000) + 2 / 200_000 + 1 / 200_001).to_numpy(), np.ones(len(ds_nb), dtype=bool), 'consistency')
    out = m[['scope', 'group', 'pair', 'K', 'G', 'pooled_delta', 'wins', 'losses', 'ties', 'n_nonzero', 'exact_p', 'min_attainable_p', 'holm_p', 'holm_p_nb', 'ci_low', 'ci_high', 'boot_share_ge_0', 'boot_share_le_0', 'ci_excludes_0', 'majority_of_all_rois', 'd5_bar_met', 't_ci_low', 't_ci_high', 't_excludes_0', 'n_same', 'majority_all', 'majority_nonzero', 'n_structural_zero_rois', 'majority_excl_structural_zero']].rename(columns={'holm_p_nb': 'nb_holm_p', 'ci_low': 'nb_boot_ci_low', 'ci_high': 'nb_boot_ci_high', 'ci_excludes_0': 'nb_boot_excludes_0', 'd5_bar_met': 'nb_d5_bar_met'})
    out['d5_bar_under_t'] = (out.t_excludes_0 & out.majority_all).where(out.scope == 'subset')
    print(f"[inference] max |mc_p - exact_p| in delta_stats: {float((ds_nb.mc_p - ds_nb.exact_p).abs().max()):.3e}")
    return out


def bootstrap(DP, I):
    """
        Reproduce the notebook's cluster bootstrap with its own RNG spec (determinism), then re-run it under 200 other seeds.

        DP (pd.DataFrame): my delta_per_roi.
        I (pd.DataFrame): inference table (gets bootstrap columns merged in).

        Returns tuple[pd.DataFrame, pd.DataFrame]: the inference table with reproduction and seed-sensitivity columns, and the per-seed long table.
    """
    B, SEED = 10_000, 20260916
    order = {s: sorted(subset_rows(DP[(DP.pair == PAIRS[0][0] + ' - ' + PAIRS[0][1]) & (DP.K == 10)], s).file_name) for s in SUBSETS}
    return_rows, seed_rows = [], []
    # the notebook's ROI order within a subset is domain, subset, file_name; reproduce it for the determinism check
    meta = DP.groupby('file_name')[['subset', 'domain']].first().reset_index().sort_values(['domain', 'subset', 'file_name'])
    nb_order = {s: subset_rows(meta, s).file_name.tolist() for s in SUBSETS}
    for i, s in enumerate(SUBSETS):
        G = len(nb_order[s])
        idx_nb = np.random.default_rng([SEED, i]).integers(G, size=(B, G))
        for a, b in PAIRS:
            for k in KS:
                v = DP[(DP.pair == f'{a} - {b}') & (DP.K == k)].set_index('file_name').loc[nb_order[s], 'dsum'].to_numpy(np.int64)
                den = G * 3 * k
                sums = v[idx_nb].sum(axis=1)
                lo, hi = np.percentile(sums / den, [2.5, 97.5])
                return_rows.append(dict(group=s, pair=f'{a} - {b}', K=k, repro_ci_low=lo, repro_ci_high=hi, repro_share_ge_0=float((sums >= 0).mean()), repro_share_le_0=float((sums <= 0).mean())))
                excl = []
                for rs in range(200):
                    idx = np.random.default_rng([777, rs, i]).integers(G, size=(B, G))
                    ss = v[idx].sum(axis=1) / den
                    l2, h2 = np.percentile(ss, [2.5, 97.5])
                    excl.append(bool(l2 > 0 or h2 < 0))
                    seed_rows.append(dict(group=s, pair=f'{a} - {b}', K=k, boot_seed=rs, ci_low=l2, ci_high=h2, excludes_0=excl[-1]))
                return_rows[-1]['share_of_200_seeds_ci_excludes_0'] = float(np.mean(excl))
    Rb = pd.DataFrame(return_rows)
    I = I.merge(Rb, on=['group', 'pair', 'K'], how='left')
    sub = I[I.scope == 'subset']
    compare('bootstrap', 'CI low/high with the notebook RNG spec (54 values)', sub[['repro_ci_low', 'repro_ci_high']].to_numpy(), sub[['nb_boot_ci_low', 'nb_boot_ci_high']].to_numpy(), 'consistency', 1e-12)
    compare('bootstrap', 'share >=0 / <=0 with notebook RNG spec (54 values)', sub[['repro_share_ge_0', 'repro_share_le_0']].to_numpy(), sub[['boot_share_ge_0', 'boot_share_le_0']].to_numpy(), 'consistency', 1e-12)
    I['d5_bar_seed_fragile'] = ((I.share_of_200_seeds_ci_excludes_0 > 0) & (I.share_of_200_seeds_ci_excludes_0 < 1)).where(I.scope == 'subset')
    return I, pd.DataFrame(seed_rows)


# ----------------------------------------------------------------------------------------------
# 8. Contested-mitosis sensitivity
# ----------------------------------------------------------------------------------------------
def contested(T, pr, DP):
    """
        Recompute the three contested-mitosis rules from my own buckets and the DB's vote labels.

        T (dict): artifact tables.
        pr (pd.DataFrame): my per-run frame with tp and contested counts.
        DP (pd.DataFrame): my delta_per_roi.

        Returns pd.DataFrame: one row per (subset, rule, pair, K).
    """
    rows = []
    ix = pr.set_index(['file_name', 'seed_index', 'condition'])
    meta = pr.groupby('file_name')[['subset', 'domain']].first()
    B, SEED = 10_000, 20260916
    nb_order_meta = meta.reset_index().sort_values(['domain', 'subset', 'file_name'])
    for si, sub in enumerate(SUBSETS):
        files = subset_rows(nb_order_meta, sub).file_name.tolist()
        G = len(files)
        idx = np.random.default_rng([SEED, si]).integers(G, size=(B, G))
        for rule in ('all_mitoses', 'contested_excluded', 'contested_as_fp'):
            for a, b in PAIRS:
                for k in KS:
                    fr = []
                    for fn in files:
                        tot = Fraction(0)
                        for s in SEEDS:
                            va, vb = ix.loc[(fn, s, a)], ix.loc[(fn, s, b)]
                            def val(r):
                                tp, ch = int(r[f'tp{k}']), int(r[f'my_contested_{k}'])
                                return Fraction(tp, k) if rule == 'all_mitoses' else (Fraction(tp - ch, k) if rule == 'contested_as_fp' else Fraction(tp - ch, k - ch))
                            tot += val(va) - val(vb)
                        fr.append(tot)
                    scale = math.lcm(*(x.denominator for x in fr))
                    ints = np.array([int(x * scale) for x in fr], dtype=np.int64)
                    if rule == 'contested_excluded':
                        p = signflip_mc(ints.astype(float), 1_000_000, 4242 + si * 100 + k + PAIRS.index((a, b)) * 1000)
                        test = 'mc_1e6_own_seed'
                    else:
                        p, _ = signflip_exact(ints)
                        test = 'exact'
                    sums = ints[idx].sum(axis=1)
                    lo, hi = np.percentile(sums / (scale * G * 3), [2.5, 97.5])
                    pooled = float(sum(fr) / (G * 3))
                    w, l = int((ints > 0).sum()), int((ints < 0).sum())
                    same = w if pooled > 0 else (l if pooled < 0 else 0)
                    rows.append(dict(subset=sub, rule=rule, pair=f'{a} - {b}', K=k, pooled_delta=pooled, wins=w, losses=l, ties=G - w - l, p=p, test=test, ci_low=lo, ci_high=hi, majority=bool(same > G / 2), d5=bool((lo > 0 or hi < 0) and same > G / 2)))
    C = pd.DataFrame(rows)
    C['holm_p'] = np.nan
    for fam in ('all_49', 'testing_35'):
        for rule in C.rule.unique():
            msk = (C.subset == fam) & (C.rule == rule)
            C.loc[msk, 'holm_p'] = holm(C.loc[msk, 'p'])
    m = C.merge(T['contested_sensitivity'], on=['subset', 'rule', 'pair', 'K'], validate='1:1', suffixes=('', '_nb'))
    compare('contested', 'pooled_delta_precision (81 rows)', m.pooled_delta, m.pooled_delta_precision, 'independent', 1e-12)
    compare('contested', 'wins/losses/ties (243 values)', m[['wins', 'losses', 'ties']].to_numpy(), m[['wins_nb', 'losses_nb', 'ties_nb']].to_numpy(), 'independent')
    ex = m[m.rule != 'contested_excluded']
    compare('contested', 'exact p, all_mitoses + contested_as_fp (54 rows)', ex.p, ex.p_nb, 'independent', 1e-12)
    exf = ex[ex.subset.isin(['all_49', 'testing_35'])]
    compare('contested', 'Holm p, exact rules (36 family rows)', exf.holm_p, exf.holm_p_nb, 'independent', 1e-12)
    mc = m[m.rule == 'contested_excluded']
    se = np.sqrt(np.maximum(mc.p, 1e-6) * (1 - mc.p) / 200_000)
    compare('contested', 'contested_excluded MC p (own 1e6 draws) within 5 SE(200k) + 2/200k of notebook MC p (27 rows)', (np.abs(mc.p - mc.p_nb) <= 5 * se + 2e-5).to_numpy(), np.ones(len(mc), dtype=bool), 'independent')
    compare('contested', 'bootstrap CI with notebook RNG spec (162 values)', m[['ci_low', 'ci_high']].to_numpy(), m[['ci_low_nb', 'ci_high_nb']].to_numpy(), 'consistency', 1e-12)
    compare('contested', 'd5_bar_met (81 rows)', m.d5.astype(bool), m.d5_bar_met.astype(bool), 'consistency')
    shares = pr.groupby('condition')[[f'my_contested_{k}' for k in KS] + [f'tp{k}' for k in KS]].sum()
    share_vals = [shares.loc[c, f'my_contested_{k}'] / shares.loc[c, f'tp{k}'] for c in CONDS for k in KS]
    prose(27, 'contested share of top-K TPs spans 23.1-25.8% across conditions and K', '0.231-0.258', f'{min(share_vals):.3f}-{max(share_vals):.3f}')
    a49 = m[(m.subset == 'all_49') & (m.K.isin([20, 30]))]
    prose(27, 'every all_49 K=20/30 contrast keeps Holm p <= 0.05 and D5 bar under all three rules', 'True', str(bool((a49.holm_p <= 0.05).all() and a49.d5.all())))
    t35 = m[(m.subset == 'testing_35') & (m.rule == 'contested_as_fp') & (m.pair == 'hem_bbox - gray_bbox') & (m.K == 20)]
    prose(27, 'testing_35 hem_bbox - gray_bbox K=20 under contested_as_fp Holm p', 0.065, round(float(t35.holm_p.iat[0]), 3), 0.0005)
    x = m[(m.subset == 'all_49') & (m.rule == 'contested_as_fp') & (m.pair == 'hem_bbox - default_51') & (m.K == 20)]
    prose(27, 'hem_bbox - default_51 K=20 under contested_as_fp (points)', -4.52, round(100 * float(x.pooled_delta.iat[0]), 2), 0.005)
    return m


# ----------------------------------------------------------------------------------------------
# 9. Table A geometry, refusals, Gate B, verification table
# ----------------------------------------------------------------------------------------------
def geometry_and_gates(T, D):
    """
        Recompute Table A's summary, the refusal counts, Gate B and the verification tallies.

        T (dict): artifact tables.
        D (pd.DataFrame): my top-30 rescoring (for Gate B list comparison).

        Returns pd.DataFrame: geometry summary rows (claimed vs recomputed are logged as prose).
    """
    pr = T['per_run']
    geo = pr.pivot_table(index=['subset', 'file_name', 'seed_index'], columns='condition', values=['base_size', 'tpl_offset_px']).reset_index()
    geo.columns = ['_'.join(c).strip('_') for c in geo.columns]
    clicks = pr[pr.condition == 'default_51'][['subset', 'file_name', 'seed_index', 'n_retries', 'refused_draws']].copy()
    clicks['refused_draws'] = clicks.refused_draws.fillna('')
    clicks['default_alone_differs'] = clicks.n_retries > 0
    clicks['gray_alone_differs'] = clicks.refused_draws.map(lambda r: bool(r) and any('gray_bbox' not in e.split(':')[1] for e in r.split(';')))
    claims = {'all_49': dict(n_clicks=147, gray_median=39.0, gray_range='23-51', gray_below=116, gray_off=2.55, hem_median=29.0, hem_range='19-51', hem_below=144, hem_off=2.5, hem_lt=135, hem_eq=12, hem_gt=0, ceiling=31, refusal=26, dflt=26, gray=10), 'original_14': dict(n_clicks=42, gray_median=40.0, gray_range='23-51', gray_below=32, gray_off=2.368, hem_median=29.0, hem_range='19-41', hem_below=42, hem_off=2.693, hem_lt=41, hem_eq=1, hem_gt=0, ceiling=10, refusal=12, dflt=12, gray=6), 'testing_35': dict(n_clicks=105, gray_median=39.0, gray_range='23-51', gray_below=84, gray_off=2.915, hem_median=29.0, hem_range='19-51', hem_below=102, hem_off=2.236, hem_lt=94, hem_eq=11, hem_gt=0, ceiling=21, refusal=14, dflt=14, gray=4)}
    rows = []
    for s in SUBSETS:
        g, c = subset_rows(geo, s), subset_rows(clicks, s)
        mine = dict(n_clicks=len(g), gray_median=float(g.base_size_gray_bbox.median()), gray_range=f'{int(g.base_size_gray_bbox.min())}-{int(g.base_size_gray_bbox.max())}', gray_below=int((g.base_size_gray_bbox < 51).sum()), gray_off=round(float(g.tpl_offset_px_gray_bbox.median()), 3), hem_median=float(g.base_size_hem_bbox.median()), hem_range=f'{int(g.base_size_hem_bbox.min())}-{int(g.base_size_hem_bbox.max())}', hem_below=int((g.base_size_hem_bbox < 51).sum()), hem_off=round(float(g.tpl_offset_px_hem_bbox.median()), 3), hem_lt=int((g.base_size_hem_bbox < g.base_size_gray_bbox).sum()), hem_eq=int((g.base_size_hem_bbox == g.base_size_gray_bbox).sum()), hem_gt=int((g.base_size_hem_bbox > g.base_size_gray_bbox).sum()), ceiling=int((g.base_size_gray_bbox == 51).sum()), refusal=int((c.n_retries > 0).sum()), dflt=int(c.default_alone_differs.sum()), gray=int(c.gray_alone_differs.sum()))
        for key, v in mine.items():
            ok = prose(10, f'Table A {s} {key}', claims[s][key], v, 1e-9 if isinstance(v, float) else 0)
            rows.append(dict(subset=s, quantity=key, claimed=claims[s][key], recomputed=v, match=ok))
    ref = pd.Series([e.split(':')[1] for r in clicks.refused_draws if r for e in r.split(';')]).value_counts().to_dict()
    prose(10, 'refused draws by refusing condition(s)', "{'gray_bbox+hem_bbox': 16, 'hem_bbox': 10, 'gray_bbox': 2}", str({k: ref[k] for k in ('gray_bbox+hem_bbox', 'hem_bbox', 'gray_bbox')}))
    # refusals are internally consistent with n_retries
    compare('geometry', 'n_retries == number of refused_draws entries (147 clicks)', clicks.n_retries.to_numpy(), clicks.refused_draws.map(lambda r: len(r.split(';')) if r else 0).to_numpy(), 'independent')
    # Gate B from my own top-30 rows
    top = T['top30']
    lists = {k: g.sort_values('rank')[['cx', 'cy', 'od']].to_numpy() for k, g in top.groupby(['file_name', 'seed_index', 'condition'])}
    gb = []
    for (fn, s), g in pr.groupby(['file_name', 'seed_index']):
        r = g.set_index('condition')
        for a, b in PAIRS:
            ca = (int(round(float(r.loc[a, 'tpl_cx']))), int(round(float(r.loc[a, 'tpl_cy']))))
            cb = (int(round(float(r.loc[b, 'tpl_cx']))), int(round(float(r.loc[b, 'tpl_cy']))))
            same_geo = int(r.loc[a, 'base_size']) == int(r.loc[b, 'base_size']) and ca == cb
            same_px = r.loc[a, 'tpl_sha1'] == r.loc[b, 'tpl_sha1']
            if same_geo or same_px:
                gb.append(dict(file_name=fn, seed_index=s, pair=f'{a} - {b}', same_geo=same_geo, same_px=same_px, offset=round(float(np.hypot(r.loc[a, 'tpl_cx'] - r.loc[b, 'tpl_cx'], r.loc[a, 'tpl_cy'] - r.loc[b, 'tpl_cy'])), 3), same_list=bool(np.array_equal(lists[(fn, s, a)], lists[(fn, s, b)])), same_tp=all(r.loc[a, f'tp_at_{k}'] == r.loc[b, f'tp_at_{k}'] for k in KS)))
    GB = pd.DataFrame(gb)
    prose(8, 'Gate B cases (identical template)', 17, len(GB))
    prose(8, 'Gate B gray_bbox = default_51 cases', 11, int((GB.pair == 'gray_bbox - default_51').sum()))
    prose(8, 'Gate B hem_bbox = gray_bbox cases', 6, int((GB.pair == 'hem_bbox - gray_bbox').sum()))
    prose(8, 'Gate B cases with identical top-30 and TP', 17, int((GB.same_list & GB.same_tp).sum()))
    prose(8, 'Gate B geometry and SHA tests agree', 17, int((GB.same_geo == GB.same_px).sum()))
    prose(27, 'Gate B sub-pixel centre offsets', '0.5,0.5,0.707', ','.join(str(x) for x in sorted(GB.offset[GB.offset > 0])))
    # verification table
    V = T['verification']
    halt_fail = int(((V.severity == 'halt') & ~V.passed).sum())
    prose(8, 'verification records / halt failures', '2600/0', f'{len(V)}/{halt_fail}')
    prose(6, 'self-hits removed per run {0: 11, 1: 430}', "{0: 11, 1: 430}", str(pr.n_self_hits.value_counts().sort_index().to_dict()))
    prose(6, 'n_detections min/max', '81/100', f'{pr.n_detections.min()}/{pr.n_detections.max()}')
    prose(6, 'n_peaks == 100 on all runs', 441, int((pr.n_peaks == 100).sum()))
    prose(6, 'runs with a candidate within one match radius of the click', 0, int((pr.n_near_click > 0).sum()))
    prose(6, 'min joint-valid pool (201.tiff;364.tiff)', '6 (201.tiff;364.tiff)', f"{pr.n_joint_valid.min()} ({';'.join(sorted(pr.loc[pr.n_joint_valid == pr.n_joint_valid.min(), 'file_name'].unique()))})")
    prose(6, 'contested seed tier used on 0 clicks', 0, int(pr.contested_seed_tier.sum()))
    # annulus check from my own top-30 geometry: nothing within one match radius of the click in the top 30
    near = []
    for (fn, s, c), g in top.groupby(['file_name', 'seed_index', 'condition']):
        run = pr[(pr.file_name == fn) & (pr.seed_index == s) & (pr.condition == c)].iloc[0]
        near.append(int((np.hypot(g.cx - run.click_cx, g.cy - run.click_cy) <= run.match_radius_px).sum()))
    compare('geometry', 'top-30 detections within one match radius of the click (all 0)', np.array(near), np.zeros(len(near)), 'independent')
    # Gate A against the committed 14-ROI artifacts (consistency with a prior artifact)
    ref = pd.read_csv(f'{REF14}_top30.csv', float_precision='round_trip')
    mine0 = top[(top.subset == 'original_14') & (top.seed_index == 0)]
    ga = mine0.merge(ref, on=['file_name', 'condition', 'rank'], suffixes=('', '_ref'), validate='1:1')
    compare('gateA', 'seed-0 original_14 top-30 cx,cy,od,bucket,matched_ann_id vs committed 14-ROI top30 (1260 rows x 5)', ga[['cx', 'cy', 'od', 'matched_ann_id']].to_numpy(float), ga[['cx_ref', 'cy_ref', 'od_ref', 'matched_ann_id_ref']].to_numpy(float), 'consistency', 0.0)
    compare('gateA', 'bucket vs committed 14-ROI top30', ga.bucket, ga.bucket_ref, 'consistency')
    return pd.DataFrame(rows), GB


# ----------------------------------------------------------------------------------------------
# 10. Prose claims of the closing summary and printed tables
# ----------------------------------------------------------------------------------------------
def closing_claims(T, I, pr, DP):
    """
        Check the quantified sentences of cell 27 and the headline tables against my recomputation.

        T (dict): artifact tables.
        I (pd.DataFrame): my inference table.
        pr (pd.DataFrame): my per-run frame.
        DP (pd.DataFrame): my delta_per_roi.

        Returns None: claims are appended to PROSE.
    """
    S = T['summary']
    pooled = {(c, k): round(pr[pr.condition == c][f'tp{k}'].sum() / (147 * k), 4) for c in CONDS for k in KS}
    claimed = {('default_51', 10): 0.5694, ('gray_bbox', 10): 0.5612, ('hem_bbox', 10): 0.5395, ('default_51', 20): 0.4969, ('gray_bbox', 20): 0.4735, ('hem_bbox', 20): 0.4405, ('default_51', 30): 0.4331, ('gray_bbox', 30): 0.4136, ('hem_bbox', 30): 0.3798}
    for key, v in claimed.items():
        prose(27, f'all_49 pooled precision {key}', v, pooled[key], 5e-5)
    head = [('gray_bbox - default_51', 10, -0.82, -2.52, 0.95, '15/18/16', 0.40, 0.40, False), ('gray_bbox - default_51', 20, -2.35, -3.78, -0.92, '12/32/5', 0.0025, 0.0100, True), ('gray_bbox - default_51', 30, -1.95, -3.11, -0.79, '9/29/11', 0.0014, 0.0068, True), ('hem_bbox - default_51', 10, -2.99, -5.78, 0.00, '11/25/13', 0.052, 0.10, False), ('hem_bbox - default_51', 20, -5.65, -8.06, -3.37, '8/37/4', 5.5e-6, 3.9e-5, True), ('hem_bbox - default_51', 30, -5.33, -7.30, -3.47, '7/36/6', 5.0e-7, 4.5e-6, True), ('hem_bbox - gray_bbox', 10, -2.18, -4.15, -0.34, '7/23/19', 0.034, 0.10, False), ('hem_bbox - gray_bbox', 20, -3.30, -5.20, -1.73, '8/31/10', 5.1e-5, 3.1e-4, True), ('hem_bbox - gray_bbox', 30, -3.38, -5.01, -1.97, '5/33/11', 2.2e-6, 1.7e-5, True)]
    a = I[(I.scope == 'subset') & (I.group == 'all_49')].set_index(['pair', 'K'])
    for pair, k, dpts, lo, hi, wlt, p, hp, d5 in head:
        r = a.loc[(pair, k)]
        prose(27, f'{pair} K={k} pooled delta pts', dpts, round(100 * r.pooled_delta, 2), 0.005)
        prose(27, f'{pair} K={k} CI pts', f'[{lo:+.2f}, {hi:+.2f}]', f'[{100 * r.repro_ci_low:+.2f}, {100 * r.repro_ci_high:+.2f}]')
        prose(27, f'{pair} K={k} W/L/T', wlt, f'{r.wins}/{r.losses}/{r.ties}')
        prose(27, f'{pair} K={k} exact p (2 sf)', p, float(f'{r.exact_p:.2g}'), abs(p) * 0.06)
        prose(27, f'{pair} K={k} Holm p (2 sf)', hp, float(f'{r.holm_p:.2g}'), abs(hp) * 0.06)
        prose(27, f'{pair} K={k} D5 bar', d5, bool(r.nb_d5_bar_met))
    prose(27, 'min_attainable_p at most 2.3e-10 in every all_49 row', '<=2.3e-10', f"<={a.min_attainable_p.max():.2g}")
    prose(27, 'hem_bbox - default_51 K=10: 2.5% of resamples at or above 0', 0.025, round(float(a.loc[('hem_bbox - default_51', 10)].repro_share_ge_0), 3), 0.0005)
    t = I[(I.scope == 'subset') & (I.group == 'testing_35')].set_index(['pair', 'K'])
    o = I[(I.scope == 'subset') & (I.group == 'original_14')].set_index(['pair', 'K'])
    same_sign = sum(np.sign(t.loc[(f'{x} - {y}', k)].pooled_delta) == np.sign(o.loc[(f'{x} - {y}', k)].pooled_delta) != 0 for x, y in PAIRS for k in KS)
    overlap = sum(max(t.loc[(f'{x} - {y}', k)].repro_ci_low, o.loc[(f'{x} - {y}', k)].repro_ci_low) <= min(t.loc[(f'{x} - {y}', k)].repro_ci_high, o.loc[(f'{x} - {y}', k)].repro_ci_high) for x, y in PAIRS for k in KS)
    prose(27, 'replication: same sign 9/9', 9, int(same_sign))
    prose(27, 'replication: CIs overlap 9/9', 9, int(overlap))
    for pair, k, hp in [('gray_bbox - default_51', 20, 0.0325), ('gray_bbox - default_51', 30, 0.0014), ('hem_bbox - default_51', 20, 0.0014), ('hem_bbox - default_51', 30, 2.8e-5), ('hem_bbox - gray_bbox', 20, 0.0037), ('hem_bbox - gray_bbox', 30, 5.4e-4)]:
        prose(27, f'testing_35 {pair} K={k} Holm p', hp, float(f'{t.loc[(pair, k)].holm_p:.2g}'), hp * 0.06)
        prose(27, f'testing_35 {pair} K={k} D5 bar met', True, bool(t.loc[(pair, k)].nb_d5_bar_met))
    prose(27, 'testing_35: no K=10 Holm p below 0.17', True, bool((t.xs(10, level='K').holm_p >= 0.17).all()))
    r = t.loc[('hem_bbox - gray_bbox', 10)]
    prose(27, 'testing_35 hem_bbox - gray_bbox K=10 CI and ROIs', '[-4.38, -0.10], 18/35', f'[{100 * r.repro_ci_low:+.2f}, {100 * r.repro_ci_high:+.2f}], {r.losses}/35')
    for pair, k, p in [('hem_bbox - default_51', 20, 0.013), ('hem_bbox - default_51', 30, 0.042), ('hem_bbox - gray_bbox', 20, 0.046), ('hem_bbox - gray_bbox', 30, 0.018)]:
        prose(27, f'original_14 {pair} K={k} exact p', p, round(float(o.loc[(pair, k)].exact_p), 3), 0.0005)
        prose(27, f'original_14 {pair} K={k} D5 bar met (bootstrap)', True, bool(o.loc[(pair, k)].nb_d5_bar_met))
    prose(27, 'original_14 gray_bbox - default_51 CI includes 0 at every K', True, bool(not o.xs('gray_bbox - default_51', level='pair').nb_boot_excludes_0.astype(bool).any()))
    dom = I[I.scope == 'domain']
    gd = dom[dom.pair == 'gray_bbox - default_51'].set_index(['group', 'K'])
    human = ['human breast cancer', 'human melanoma', 'human neuroendocrine tumor']
    prose(27, 'gray_bbox - default_51 negative at every K in all three human domains', True, bool(all(gd.loc[(h, k)].pooled_delta < 0 for h in human for k in KS)))
    vals = {h: '/'.join(f'{100 * gd.loc[(h, k)].pooled_delta:.2f}' for k in KS) for h in human + ['canine cutaneous mast cell tumor', 'canine lung cancer']}
    prose(27, 'human breast / melanoma / NET gray-default pts', '-1.43/-5.48/-5.24; -3.81/-4.76/-3.65; -5.71/-3.81/-3.33', '; '.join(vals[h] for h in human))
    prose(27, 'mast cell gray-default +0.48..+1.43 at every K; lung +1.43/0/0', '1.43/0.95/0.48; 1.43/0.00/0.00', f"{vals['canine cutaneous mast cell tumor']}; {vals['canine lung cancer']}")
    hd = dom[dom.pair == 'hem_bbox - default_51']
    prose(27, 'hem_bbox - default_51 negative in 6/7 domains at K=20 and 7/7 at K=30', '6/7;7/7', f"{int((hd[hd.K == 20].pooled_delta < 0).sum())}/7;{int((hd[hd.K == 30].pooled_delta < 0).sum())}/7")
    maxabs = dom.groupby('group').pooled_delta.apply(lambda v: v.abs().max())
    prose(27, 'canine lung is the one domain where no refinement moves precision by more than 1.9 pts at any K', 'canine lung cancer only, 1.9', f"{';'.join(maxabs[maxabs <= 0.0191].index)} only, {100 * maxabs['canine lung cancer']:.1f}")
    z = pr[pr.file_name == '211.tiff']
    prose(27, '211.tiff has 0 TPs in the top 30 on all 9 runs', 9, int((z.tp30 == 0).sum()))
    o14 = pr[pr.subset == 'original_14'].groupby(['file_name', 'condition']).tp20.sum().unstack() / 60
    prose(27, '245.tiff K=20 original_14 worst ROI: 0.0833/0.1167/0.1000', '0.0833/0.1167/0.1000', '/'.join(f'{o14.loc["245.tiff", c]:.4f}' for c in CONDS))
    worst = S[S.subset == 'original_14'].worst_roi_files
    prose(27, 'original_14 worst ROI is 245.tiff under every condition (every K)', True, bool((worst == '245.tiff').all()))
    geo = T['per_run'].pivot_table(index=['file_name', 'seed_index'], columns='condition', values='base_size')
    prose(27, 'median template size 51 / 39 / 29', '51/39/29', '/'.join(str(int(geo[c].median())) for c in CONDS))


# ----------------------------------------------------------------------------------------------
# 10b. Figures: the series each plotting cell passes in, and numbers read off the rendered PNGs
# ----------------------------------------------------------------------------------------------
def figure_checks(T, I, pr, DP):
    """
        Step 1a: re-derive each figure's plotted series; step 2: numbers read by eye off the decoded PNGs, checked against the recomputation.

        T (dict): artifact tables.
        I (pd.DataFrame): inference table with reproduced bootstrap columns.
        pr (pd.DataFrame): my per-run frame.
        DP (pd.DataFrame): my delta_per_roi.

        Returns None: comparisons and claims are recorded.
    """
    boot = T['bootstrap_ci']
    pp = boot[boot.quantity == 'pooled_precision'].copy()
    mine = []
    for r in pp.itertuples():
        runs = subset_rows(pr, r.subset)
        mine.append(runs[runs.condition == r.condition_or_pair][f'tp{r.K}'].sum() / (len(runs[runs.condition == r.condition_or_pair]) * r.K))
    compare('figures', 'Figure 1 (cell 23) bar heights: pooled precision, 27 bars', np.array(mine), pp.estimate.to_numpy(), 'independent', 1e-12)
    read_fig1 = {('all_49', 'default_51', 10): 0.569, ('all_49', 'gray_bbox', 10): 0.561, ('all_49', 'hem_bbox', 10): 0.539, ('all_49', 'default_51', 20): 0.497, ('all_49', 'gray_bbox', 20): 0.473, ('all_49', 'hem_bbox', 20): 0.440, ('all_49', 'default_51', 30): 0.433, ('all_49', 'gray_bbox', 30): 0.414, ('all_49', 'hem_bbox', 30): 0.380, ('original_14', 'default_51', 10): 0.676, ('original_14', 'gray_bbox', 10): 0.671, ('original_14', 'hem_bbox', 10): 0.652, ('original_14', 'default_51', 20): 0.625, ('original_14', 'gray_bbox', 20): 0.604, ('original_14', 'hem_bbox', 20): 0.561, ('original_14', 'default_51', 30): 0.564, ('original_14', 'gray_bbox', 30): 0.557, ('original_14', 'hem_bbox', 30): 0.514, ('testing_35', 'default_51', 10): 0.527, ('testing_35', 'gray_bbox', 10): 0.517, ('testing_35', 'hem_bbox', 10): 0.494, ('testing_35', 'default_51', 20): 0.446, ('testing_35', 'gray_bbox', 20): 0.421, ('testing_35', 'hem_bbox', 20): 0.392, ('testing_35', 'default_51', 30): 0.381, ('testing_35', 'gray_bbox', 30): 0.356, ('testing_35', 'hem_bbox', 30): 0.326}
    got = {(r.subset, r.condition_or_pair, r.K): v for r, v in zip(pp.itertuples(), mine)}
    prose(23, 'Figure 1 rendered bar labels (27, read off the PNG)', 27, sum(abs(round(got[k], 3) - v) < 1e-9 for k, v in read_fig1.items()))
    dm = DP.set_index(['file_name', 'pair', 'K']).dsum / DP.set_index(['file_name', 'pair', 'K']).index.get_level_values('K')
    nbv = T['delta_per_roi'].set_index(['file_name', 'pair', 'K']).delta_precision_sum
    compare('figures', 'Figure 2 (cell 24) heatmap cells: seed-summed precision delta, 441', dm.loc[nbv.index].to_numpy(), nbv.to_numpy(), 'independent', 1e-12)
    read_fig2 = {('293.tiff', 'hem_bbox - default_51', 10): 1.20, ('460.tiff', 'hem_bbox - default_51', 20): -1.20, ('344.tiff', 'hem_bbox - gray_bbox', 30): -0.63, ('072.tiff', 'gray_bbox - default_51', 30): -0.53, ('301.tiff', 'gray_bbox - default_51', 30): 0.33, ('460.tiff', 'hem_bbox - gray_bbox', 20): -1.05, ('343.tiff', 'hem_bbox - default_51', 10): -0.60, ('220.tiff', 'hem_bbox - default_51', 10): 0.40}
    prose(24, 'Figure 2 rendered cells spot-checked (8, read off the PNG)', 8, sum(abs(round(float(dm.loc[k]), 2) - v) < 1e-9 for k, v in read_fig2.items()))
    lim = float(np.abs(dm).max())
    prose(24, 'Figure 2 shared colour limit = max |cell| (colourbar reads to about +-1.2)', 1.2, round(lim, 2), 0.001)
    dom = I[I.scope == 'domain']
    compare('figures', 'Figure 3 (cell 25) cells: per-domain pooled delta, 63', (100 * dom.pooled_delta).round(1).to_numpy(), (100 * T['delta_stats'].set_index(['scope', 'group', 'pair', 'K']).loc[list(zip(dom.scope, dom.group, dom.pair, dom.K)), 'pooled_delta_precision']).round(1).to_numpy(), 'independent', 1e-9)
    read_fig3 = {('canine cutaneous mast cell tumor', 'hem_bbox - gray_bbox', 30): (-8.6, '0/7/0'), ('canine soft tissue sarcoma', 'hem_bbox - default_51', 20): (-9.8, '0/6/1'), ('canine lymphosarcoma', 'gray_bbox - default_51', 10): (1.0, '3/3/1'), ('human neuroendocrine tumor', 'gray_bbox - default_51', 20): (-3.8, '0/7/0'), ('canine lung cancer', 'hem_bbox - default_51', 10): (1.9, '2/1/4'), ('human breast cancer', 'hem_bbox - gray_bbox', 10): (0.0, '2/2/3')}
    di = dom.set_index(['group', 'pair', 'K'])
    prose(25, 'Figure 3 rendered cells spot-checked (6 deltas + W/L/T, read off the PNG)', 6, sum(abs(round(100 * di.loc[k].pooled_delta, 1) - v[0]) < 1e-9 and f'{di.loc[k].wins}/{di.loc[k].losses}/{di.loc[k].ties}' == v[1] for k, v in read_fig3.items()))
    a = I[(I.scope == 'subset') & (I.group == 'all_49')].set_index(['pair', 'K'])
    read_fig4 = {('gray_bbox - default_51', 10): ('15/18/16', '0.404'), ('gray_bbox - default_51', 20): ('12/32/5', '0.00998'), ('gray_bbox - default_51', 30): ('9/29/11', '0.00684'), ('hem_bbox - default_51', 10): ('11/25/13', '0.104'), ('hem_bbox - default_51', 20): ('8/37/4', '3.85e-05'), ('hem_bbox - default_51', 30): ('7/36/6', '4.47e-06'), ('hem_bbox - gray_bbox', 10): ('7/23/19', '0.102'), ('hem_bbox - gray_bbox', 20): ('8/31/10', '0.000308'), ('hem_bbox - gray_bbox', 30): ('5/33/11', '1.74e-05')}
    prose(26, 'Figure 4 rendered stacks and Holm labels (9 bars, read off the PNG)', 9, sum(f'{a.loc[k].wins}/{a.loc[k].losses}/{a.loc[k].ties}' == v[0] and f'{a.loc[k].holm_p:.3g}' == v[1] for k, v in read_fig4.items()))


# ----------------------------------------------------------------------------------------------
# 11. Extensions for Step 4
# ----------------------------------------------------------------------------------------------
def extensions(T, pr, DP, I):
    """
        Ceiling ties, worst-click product metrics, recall@K sensitivity, fragility, domain structure, 211.tiff.

        T (dict): artifact tables.
        pr (pd.DataFrame): my per-run frame.
        DP (pd.DataFrame): my delta_per_roi.
        I (pd.DataFrame): my inference table.

        Returns dict[str, pd.DataFrame]: extension tables.
    """
    out = {}
    # (a) precision ceiling: min(n_gt_eval_mitotic, K) / K
    rows = []
    for k in KS:
        cap = np.minimum(pr.my_n_gt_mitotic, k)
        for c in CONDS:
            m = pr.condition == c
            rows.append(dict(K=k, condition=c, clicks_with_ceiling_below_K=int((cap[m] < k).sum()), clicks_at_ceiling=int((pr.loc[m, f'tp{k}'] == cap[m]).sum()), clicks_within_1_of_ceiling=int((pr.loc[m, f'tp{k}'] >= cap[m] - 1).sum())))
    ceil = pd.DataFrame(rows)
    wide = pr.pivot_table(index=['file_name', 'seed_index'], columns='condition', values=[f'tp{k}' for k in KS] + ['my_n_gt_mitotic'])
    tie_rows = []
    for a, b in PAIRS:
        for k in KS:
            f = DP[(DP.pair == f'{a} - {b}') & (DP.K == k)].set_index('file_name')
            tied = f[f.dsum == 0].index
            both0 = [fn for fn in tied if (wide.loc[fn][(f'tp{k}', a)] == 0).all() and (wide.loc[fn][(f'tp{k}', b)] == 0).all()]
            capped = [fn for fn in tied if all(wide.loc[(fn, s)][(f'tp{k}', a)] == min(k, wide.loc[(fn, s)][('my_n_gt_mitotic', a)]) and wide.loc[(fn, s)][(f'tp{k}', b)] == min(k, wide.loc[(fn, s)][('my_n_gt_mitotic', a)]) for s in SEEDS)]
            per_click_identical = [fn for fn in tied if all(wide.loc[(fn, s)][(f'tp{k}', a)] == wide.loc[(fn, s)][(f'tp{k}', b)] for s in SEEDS)]
            tie_rows.append(dict(pair=f'{a} - {b}', K=k, roi_ties=len(tied), ties_all_zero_tp=len(both0), ties_both_at_ceiling_every_click=len(capped), ties_identical_on_every_click=len(per_click_identical), ties_cancelling_across_clicks=len(tied) - len(per_click_identical)))
    out['ceiling'] = ceil
    out['ties'] = pd.DataFrame(tie_rows)
    # (b) product metrics: per-click deltas, worst click per ROI, worst click per subset
    prod_rows = []
    for a, b in PAIRS:
        for k in KS:
            d = (wide[(f'tp{k}', a)] - wide[(f'tp{k}', b)])
            worst_a = wide[(f'tp{k}', a)].groupby('file_name').min()
            worst_b = wide[(f'tp{k}', b)].groupby('file_name').min()
            wd = (worst_a - worst_b).to_numpy()
            p_w, nz_w = signflip_exact(wd)
            prod_rows.append(dict(pair=f'{a} - {b}', K=k, clicks=len(d), clicks_worse=int((d < 0).sum()), clicks_better=int((d > 0).sum()), clicks_equal=int((d == 0).sum()), mean_click_delta_tp=float(d.mean()), worst_click_delta_tp=int(d.min()), best_click_delta_tp=int(d.max()), clicks_losing_ge_5_tp=int((d <= -5).sum()), clicks_gaining_ge_5_tp=int((d >= 5).sum()), worst_of_3_clicks_mean_precision_a=float(worst_a.mean() / k), worst_of_3_clicks_mean_precision_b=float(worst_b.mean() / k), worst_of_3_delta_W_L_T=f'{int((wd > 0).sum())}/{int((wd < 0).sum())}/{int((wd == 0).sum())}', worst_of_3_exact_p=p_w, worst_of_3_min_p=2.0 ** (1 - nz_w) if nz_w else 1.0))
    out['product'] = pd.DataFrame(prod_rows)
    sd_rows = []
    for c in CONDS:
        for k in KS:
            v = pr[pr.condition == c][f'tp{k}'] / k
            sd_rows.append(dict(condition=c, K=k, click_precision_sd=float(v.std(ddof=1)), click_precision_p10=float(v.quantile(0.10)), clicks_at_zero=int((v == 0).sum())))
    out['click_spread'] = pd.DataFrame(sd_rows)
    S = T['summary']
    wc = S[S.subset == 'original_14'].pivot(index='K', columns='condition', values='worst_click_precision')
    out['worst_click_original_14'] = wc.reset_index()
    # (c) recall@K sensitivity (D4's unit): per-click delta TP / n_gt_mitotic, summed per ROI; MC sign-flip
    rec_rows = []
    ngt = pr.groupby('file_name').my_n_gt_mitotic.first()
    for fam in ('all_49', 'testing_35'):
        ps = []
        for a, b in PAIRS:
            for k in KS:
                f = subset_rows(DP[(DP.pair == f'{a} - {b}') & (DP.K == k)], fam).set_index('file_name')
                v = (f.dsum / ngt.loc[f.index]).to_numpy(float)
                p = signflip_mc(v, 1_000_000, 9000 + KS.index(k) * 10 + PAIRS.index((a, b)) + (0 if fam == 'all_49' else 100))
                ps.append(p)
                rec_rows.append(dict(family=fam, pair=f'{a} - {b}', K=k, mean_roi_recall_delta_per_click=float((v / 3).mean()), recall_mc_p_1e6=p, wlt=f'{int((v > 0).sum())}/{int((v < 0).sum())}/{int((v == 0).sum())}'))
        hp = holm(ps)
        for i, r in enumerate(rec_rows[-9:]):
            r['recall_holm_p'] = hp[i]
    R = pd.DataFrame(rec_rows)
    prec_holm = I[(I.scope == 'subset') & I.group.isin(['all_49', 'testing_35'])][['group', 'pair', 'K', 'holm_p']].rename(columns={'group': 'family', 'holm_p': 'precision_holm_p'})
    R = R.merge(prec_holm, on=['family', 'pair', 'K'])
    R['verdict_same_at_0.05'] = (R.recall_holm_p <= 0.05) == (R.precision_holm_p <= 0.05)
    out['recall_sensitivity'] = R
    # (d) fragility: unit TP changes toward zero needed to lift Holm p above 0.05
    fr_rows = []
    for fam in ('all_49', 'testing_35'):
        famI = I[(I.scope == 'subset') & (I.group == fam)].reset_index(drop=True)
        for j, row in famI.iterrows():
            if row.holm_p > 0.05:
                continue
            f = subset_rows(DP[(DP.pair == row.pair) & (DP.K == row.K)], fam).sort_values('file_name')
            d = f.dsum.to_numpy(np.int64).copy()
            direction = -np.sign(d.sum())
            others = famI.exact_p.to_numpy().copy()
            steps = 0
            while True:
                p_now, _ = signflip_exact(d)
                others[j] = p_now
                if holm(others)[j] > 0.05 or steps > 400:
                    break
                best, best_p = None, -1
                for i in range(len(d)):
                    d[i] += direction
                    pi, _ = signflip_exact(d)
                    d[i] -= direction
                    if pi > best_p:
                        best, best_p = i, pi
                d[best] += direction
                steps += 1
            fr_rows.append(dict(family=fam, pair=row.pair, K=row.K, holm_p=row.holm_p, total_abs_delta_tp=int(np.abs(f.dsum).sum()), net_delta_tp=int(f.dsum.sum()), unit_changes_to_lose_holm_005=steps, tp_cells_in_family=len(f) * 3 * 2))
    out['fragility'] = pd.DataFrame(fr_rows)
    # (e) domain structure: species split of the descriptive human/canine contrast (7 domain units)
    dom = I[I.scope == 'domain']
    sp_rows = []
    human = {'human breast cancer', 'human melanoma', 'human neuroendocrine tumor'}
    for (pair, k), g in dom.groupby(['pair', 'K']):
        v = g.set_index('group').pooled_delta
        hs = [x for x in v.index if x in human]
        obs = v[hs].mean() - v[[x for x in v.index if x not in human]].mean()
        null = []
        for combo in itertools.combinations(v.index, 3):
            null.append(v[list(combo)].mean() - v[[x for x in v.index if x not in combo]].mean())
        null = np.array(null)
        sp_rows.append(dict(pair=pair, K=k, human_minus_canine_pts=100 * obs, one_sided_p_human_more_negative=float((null <= obs + 1e-12).mean()), two_sided_p=float((np.abs(null) >= abs(obs) - 1e-12).mean()), floor=1 / 35))
    out['species_split'] = pd.DataFrame(sp_rows)
    geo = T['per_run'].pivot_table(index=['file_name', 'seed_index', 'domain'], columns='condition', values='base_size').reset_index()
    kw_rows = []
    for c in ('gray_bbox', 'hem_bbox'):
        by_dom = [g[c].to_numpy() for _, g in geo.groupby('domain')]
        H, p = stats.kruskal(*by_dom)
        roi_mean = geo.groupby(['domain', 'file_name'])[c].mean()
        grand = geo[c].mean()
        between = sum(len(g) * (g[c].mean() - grand) ** 2 for _, g in geo.groupby('domain'))
        total = ((geo[c] - grand) ** 2).sum()
        kw_rows.append(dict(condition=c, kruskal_H_click_level=H, kruskal_p_click_level=p, between_domain_share_of_variance=between / total, domain_medians=';'.join(f'{d}={int(np.median(g[c]))}' for d, g in geo.groupby('domain'))))
    out['size_by_domain'] = pd.DataFrame(kw_rows)
    # (e2) domain as the exchangeable unit (premise check): 7 domain-summed deltas, exact sign-flip, Holm over 9
    du = []
    for fam in ('all_49', 'testing_35'):
        ps = []
        for a, b in PAIRS:
            for k in KS:
                f = subset_rows(DP[(DP.pair == f'{a} - {b}') & (DP.K == k)], fam)
                dsum = f.groupby('domain').dsum.sum().to_numpy()
                p, nz = signflip_exact(dsum)
                ps.append(p)
                du.append(dict(family=fam, pair=f'{a} - {b}', K=k, domain_deltas_tp=';'.join(str(int(v)) for v in dsum), n_nonzero=nz, exact_p=p, floor=2.0 ** (1 - nz) if nz else 1.0))
        hp = holm(ps)
        for i, r in enumerate(du[-9:]):
            r['holm_p'] = hp[i]
    out['domain_as_unit'] = pd.DataFrame(du)
    # (f) which ROIs carry the K=10 hem contrasts, under precision weighting vs recall weighting
    drv = []
    for pair in ('hem_bbox - default_51', 'hem_bbox - gray_bbox'):
        f = DP[(DP.pair == pair) & (DP.K == 10)].set_index('file_name')
        for fn, r in f.iterrows():
            drv.append(dict(pair=pair, K=10, file_name=fn, subset=r.subset, domain=r.domain, n_gt_mitotic=int(ngt.loc[fn]), delta_tp_sum=int(r.dsum), recall_weight_delta=float(r.dsum / ngt.loc[fn])))
    drv = pd.DataFrame(drv)
    drv['share_of_positive_tp'] = drv.groupby('pair').delta_tp_sum.transform(lambda v: v.clip(lower=0) / v.clip(lower=0).sum())
    out['k10_drivers'] = drv.sort_values(['pair', 'delta_tp_sum'])
    # (g) D10 exposure: runs whose seed peak was not removed (outside the 100-peak cap or >5 px), by condition
    ex = T['per_run'][T['per_run'].n_self_hits == 0][['file_name', 'seed_index', 'condition', 'base_size', 'n_detections', 'tp_at_10', 'tp_at_20', 'tp_at_30']]
    out['d10_exposure_no_self_hit_runs'] = ex
    return out


def precision_null(T, imgs, anns):
    """
        Length-matched random-list precision@K per ROI: n_gt * (1 - (1 - hit_frac / n_gt) ** K) / K, hit_frac = share of the ROI within one match radius of a mitosis.

        T (dict): artifact tables.
        imgs (pd.DataFrame): DB images.
        anns (pd.DataFrame): DB annotations.

        Returns pd.DataFrame: one row per (ROI, K) with the null and the observed pooled precision per condition.
    """
    from scipy.spatial import cKDTree
    pr = T['per_run']
    rows = []
    for fn, g in pr.groupby('file_name'):
        w, h = imgs.set_index('file_name').loc[fn, ['width', 'height']]
        mpp = float(g.mpp.iat[0])
        r = RADIUS_UM / mpp
        mito = anns[(anns.file_name == fn) & (anns.category == 1)]
        ys, xs = np.mgrid[0:h:16, 0:w:16]
        d, _ = cKDTree(mito[['cx', 'cy']].to_numpy(float)).query(np.stack([xs.ravel(), ys.ravel()], 1).astype(float), k=1)
        hit_frac = float((d <= r).mean())
        n_gt = len(mito) - 1
        for k in KS:
            row = dict(file_name=fn, K=k, n_gt_eval=n_gt, hit_frac=hit_frac, null_precision=n_gt * (1 - (1 - hit_frac / n_gt) ** k) / k)
            for c in CONDS:
                row[f'observed_{c}'] = float(g[g.condition == c][f'tp_at_{k}'].sum() / (3 * k))
            rows.append(row)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------------
# 12. Tier B -- pixels
# ----------------------------------------------------------------------------------------------
def load_rgb(path):
    """
        Level-0 RGB of a ROI TIFF, read directly with tifffile.

        path (str): TIFF path.

        Returns np.ndarray: HxWx3 uint8.
    """
    with tifffile.TiffFile(path) as tf:
        s = tf.series[0]
        arr = s.levels[0].asarray() if getattr(s, 'levels', None) else s.asarray()
    return np.ascontiguousarray(arr[:, :, :3])


def my_patch(img, cx, cy, size):
    """
        Square patch at Python-rounded centre, None if it leaves the image.

        img (np.ndarray): 2D image.
        cx (float): x.
        cy (float): y.
        size (int): odd side.

        Returns np.ndarray or None: the patch.
    """
    h = size // 2
    ix, iy = int(round(cx)), int(round(cy))
    H, W = img.shape[:2]
    if ix - h < 0 or iy - h < 0 or ix + h >= W or iy + h >= H:
        return None
    return img[iy - h:iy + h + 1, ix - h:ix + h + 1]


def my_otsu_box(channel, cx, cy):
    """
        Re-implementation of the D8 gate and anchor: Otsu on a 51 px window, component under the centre pixel, area/solidity gates, (x0+x1-1)/2 centre.

        channel (np.ndarray): single-channel image, object-high.
        cx (float): click x.
        cy (float): click y.

        Returns tuple[int, float, float] or None: (base_size, centre x, centre y).
    """
    p = my_patch(channel, cx, cy, 51)
    if p is None:
        return None
    u8 = cv2.normalize(p.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    _, bw = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    lab = label(bw, connectivity=2)
    lc = lab[25, 25]
    if lc == 0:
        return None
    reg = [r for r in regionprops(lab) if r.label == lc][0]
    y0, x0, y1, x1 = reg.bbox
    if reg.area < 50 or reg.area > 0.85 * p.size or reg.solidity < 0.5:
        return None
    n = max(y1 - y0, x1 - x0)
    n = n + 1 if n % 2 == 0 else n
    n = max(5, n)
    return n, int(round(cx)) - 25 + (x0 + x1 - 1) / 2.0, int(round(cy)) - 25 + (y0 + y1 - 1) / 2.0


def my_pipeline(hem, tx, ty, base, mpp, click, gt_eval):
    """
        Re-implementation of one production run: TM_CCOEFF on a replicate-padded channel, robust floor, 15 px local maxima, top-100 cap, greedy distance NMS, 5 px self-hit removal, od51 darkest-10% ranking, and scoring.

        hem (np.ndarray): hematoxylin OD channel, float32.
        tx (float): template centre x.
        ty (float): template centre y.
        base (int): template side.
        mpp (float): microns per pixel.
        click (tuple): click (x, y).
        gt_eval (pd.DataFrame): ground truth minus the click.

        Returns tuple[pd.DataFrame, dict]: ranked detections with buckets, and counts.
    """
    patch = my_patch(hem, tx, ty, 73)
    c, hb = 36, base // 2
    tmpl = np.ascontiguousarray(patch[c - hb:c + hb + 1, c - hb:c + hb + 1], dtype=np.float32)
    pad = (base - 1) // 2
    H, W = hem.shape
    padded = cv2.copyMakeBorder(hem, pad, pad, pad, pad, cv2.BORDER_REPLICATE)
    res = cv2.matchTemplate(padded, tmpl, cv2.TM_CCOEFF).astype(np.float32)
    fused = np.full((H + 2 * pad, W + 2 * pad), np.float32(-3e38))
    fused[pad:pad + res.shape[0], pad:pad + res.shape[1]] = res
    fused = fused[pad:pad + H, pad:pad + W]
    sample = fused[::8, ::8]
    med = float(np.median(sample))
    mad = 1.4826 * float(np.median(np.abs(sample - med)))
    thr = med - 1.5 * mad
    dil = cv2.dilate(fused, np.ones((15, 15), np.uint8))
    ys, xs = np.nonzero((fused >= dil) & (fused >= thr))
    sc = fused[ys, xs]
    order = np.lexsort((ys, xs, -sc))[:100]
    xy, sc = np.stack([xs[order], ys[order]], 1).astype(float), sc[order]
    r = RADIUS_UM / mpp
    keep, sup = [], np.zeros(len(xy), bool)
    for i in np.argsort(-sc, kind='stable'):
        if sup[i]:
            continue
        keep.append(i)
        sup |= np.hypot(xy[:, 0] - xy[i, 0], xy[:, 1] - xy[i, 1]) <= r
    xy, sc = xy[keep], sc[keep]
    self_hit = np.hypot(xy[:, 0] - tx, xy[:, 1] - ty) <= 5.0
    n_self = int(self_hit.sum())
    xy, sc = xy[~self_hit], sc[~self_hit]
    hp = cv2.copyMakeBorder(hem, 25, 25, 25, 25, cv2.BORDER_REPLICATE)
    od = []
    for x, y in xy:
        p = my_patch(hp, x + 25, y + 25, 51).ravel()
        kk = max(1, int(0.10 * p.size))
        od.append(float(np.partition(p, -kk)[-kk:].mean()))
    det = pd.DataFrame(dict(cx=xy[:, 0], cy=xy[:, 1], score=sc, od=od)).sort_values('od', ascending=False, kind='mergesort').reset_index(drop=True)
    m = greedy_match(det[['cx', 'cy']].to_numpy(), gt_eval[['cx', 'cy']].to_numpy(float), r)
    cat = np.where(m >= 0, gt_eval.category.to_numpy()[np.maximum(m, 0)], 0)
    det['hit'] = (m >= 0) & (cat == 1)
    det['score_rank'] = det.score.rank(ascending=False, method='first').astype(int) - 1
    return det, dict(n_peaks=len(order), n_detections=len(det), n_self_hits=n_self)


def tier_b(T, imgs, anns):
    """
        Pixel-level spot checks: every click draw on all 49 ROIs, 15 re-implemented pipeline runs, a production-draw sensitivity and a re-timing.

        T (dict): artifact tables.
        imgs (pd.DataFrame): DB images.
        anns (pd.DataFrame): DB annotations.

        Returns dict[str, pd.DataFrame]: Tier B tables.
    """
    sys.path.insert(0, REPO)
    from midog_utils import production as prod          # used only for the labelled production-draw sensitivity and re-timing
    from midog_utils import seed_selection as ss
    pr = T['per_run']
    top = T['top30']
    out = {}
    t_start = time.time()
    draw_rows, prod_draw_rows = [], []
    rerun_targets = {('211.tiff', 0), ('245.tiff', 1), ('425.tiff', 1), ('460.tiff', 1), ('293.tiff', 0)}
    pipe_rows, sens_rows = [], []
    for (sub, fn), g in pr.groupby(['subset', 'file_name'], sort=False):
        path = f'{SUBSET_DIRS[sub]}/{fn}'
        image_id = int(g.image_id.iat[0])
        rgb = load_rgb(path)
        gray = (255.0 - cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)).astype(np.float32)
        hem = rgb2hed(rgb.astype(np.float32) / 255.0)[:, :, 0].astype(np.float32)
        gt = anns[anns.file_name == fn].reset_index(drop=True)
        mito = gt[gt.category == 1]
        pool = mito[(mito.n_votes > 0) & (mito.n_mitotic_votes == mito.n_votes)]
        H, W = rgb.shape[:2]
        ixr, iyr = np.rint(pool.cx).astype(int), np.rint(pool.cy).astype(int)
        pool = pool[(ixr >= 36) & (ixr <= W - 37) & (iyr >= 36) & (iyr <= H - 37)]
        specs = {}
        for _, a in pool.iterrows():
            sp = {'default_51': (51, float(a.cx), float(a.cy)), 'gray_bbox': my_otsu_box(gray, a.cx, a.cy), 'hem_bbox': my_otsu_box(hem, a.cx, a.cy)}
            specs[int(a.ann_id)] = {c: (None if v is None or my_patch(hem, v[1], v[2], 73) is None else v) for c, v in sp.items()}
        n_joint = sum(all(v is not None for v in s.values()) for s in specs.values())
        accepted, own_draws = [], []
        for s in SEEDS:
            work = pool[~pool.ann_id.isin(accepted)].copy()
            rng = np.random.default_rng([s, image_id])
            refusals, first = [], None
            gray_alone, hem_alone = None, None
            while len(work):
                j = int(rng.integers(len(work)))
                a = work.iloc[j]
                sp = specs[int(a.ann_id)]
                if first is None:
                    first = int(a.ann_id)
                if gray_alone is None and sp['gray_bbox'] is not None:
                    gray_alone = int(a.ann_id)
                if hem_alone is None and sp['hem_bbox'] is not None:
                    hem_alone = int(a.ann_id)
                if all(v is not None for v in sp.values()):
                    break
                refusals.append(f"{int(a.ann_id)}:{'+'.join(c for c, v in sp.items() if v is None)}")
                work = work.drop(work.index[j])
            nb = g[g.seed_index == s].set_index('condition')
            for c in CONDS:
                draw_rows.append(dict(file_name=fn, seed_index=s, condition=c, my_ann=int(a.ann_id), nb_ann=int(nb.loc[c, 'seed_ann_id']), my_base=int(sp[c][0]), nb_base=int(nb.loc[c, 'base_size']), my_tx=sp[c][1], nb_tx=float(nb.loc[c, 'tpl_cx']), my_ty=sp[c][2], nb_ty=float(nb.loc[c, 'tpl_cy']), my_refused=';'.join(refusals), nb_refused='' if pd.isna(nb.loc[c, 'refused_draws']) else nb.loc[c, 'refused_draws'], my_joint_valid=n_joint, nb_joint_valid=int(nb.loc[c, 'n_joint_valid']), my_default_alone_differs=first != int(a.ann_id), my_gray_alone_differs=gray_alone != int(a.ann_id), my_hem_alone_differs=hem_alone != int(a.ann_id)))
            for c, own_ann in (('default_51', first), ('gray_bbox', gray_alone), ('hem_bbox', hem_alone)):
                if own_ann != int(a.ann_id):
                    own_draws.append((s, c, own_ann, specs[own_ann][c], 'own_draw'))
            if gray_alone != int(a.ann_id):
                own_draws.append((s, 'default_51', gray_alone, specs[gray_alone]['default_51'], 'gray_gate'))
                own_draws.append((s, 'gray_bbox', gray_alone, specs[gray_alone]['gray_bbox'], 'gray_gate'))
            # production build_seed on the gray channel, pool minus earlier accepted clicks, same stream
            gtm = mito[~mito.ann_id.isin(accepted)].rename(columns={'category': 'category_id'})
            try:
                seed = ss.build_seed(gtm, gray, np.random.default_rng([s, image_id]), rgb.shape)
                prod_draw_rows.append(dict(file_name=fn, seed_index=s, prod_gray_ann=seed.ann_id, my_gray_alone_ann=gray_alone, joint_ann=int(a.ann_id), prod_base=seed.base_size, prod_tx=seed.template_xy[0], prod_ty=seed.template_xy[1]))
            except ValueError as e:
                prod_draw_rows.append(dict(file_name=fn, seed_index=s, prod_gray_ann=-1, my_gray_alone_ann=gray_alone, joint_ann=int(a.ann_id), error=str(e)))
            accepted.append(int(a.ann_id))
        mpp = tiff_mpp(path)
        # own-draw policy arm: each condition on the click its own gate alone would have drawn (earlier clicks held at the joint choices), through the production entry point
        done = {}
        for s, c, own_ann, spec, arm in own_draws:
            if (s, c, own_ann) in done:
                sens_rows.append(dict(done[(s, c, own_ann)], arm=arm))
                continue
            sd = ss.Seed(ann_id=own_ann, click_xy=(0.0, 0.0), template_xy=(spec[1], spec[2]), base_size=int(spec[0]), recentred=c != 'default_51', offset_px=0.0, n_retries=0, agreement_flagged=False, n_agreement_pool=0, n_after_border=0)
            det, info = prod.run_production_pipeline(rgb, sd, mpp, rank_key='chromatin_od')
            gt_own = gt[gt.ann_id != own_ann].reset_index(drop=True)
            mm = greedy_match(det[['cx', 'cy']].to_numpy(float), gt_own[['cx', 'cy']].to_numpy(float), RADIUS_UM / mpp)
            cat = np.where(mm >= 0, gt_own.category.to_numpy()[np.maximum(mm, 0)], 0)
            hit = (mm >= 0) & (cat == 1)
            row = dict(subset=sub, file_name=fn, seed_index=s, condition=c, own_ann=own_ann, joint_ann=int(g[g.seed_index == s].seed_ann_id.iat[0]), base_size=int(spec[0]), n_detections=info['n_detections'], **{f'tp{k}': int(hit[:k].sum()) for k in KS})
            done[(s, c, own_ann)] = row
            sens_rows.append(dict(row, arm=arm))
        # re-implemented pipeline on the chosen (ROI, click)s
        for s in SEEDS:
            if (fn, s) not in rerun_targets:
                continue
            nb = g[g.seed_index == s].set_index('condition')
            ann = int(nb.seed_ann_id.iat[0])
            click = (float(nb.click_cx.iat[0]), float(nb.click_cy.iat[0]))
            gt_eval = gt[gt.ann_id != ann].reset_index(drop=True)
            for c in CONDS:
                t0 = time.time()
                det, info = my_pipeline(hem, float(nb.loc[c, 'tpl_cx']), float(nb.loc[c, 'tpl_cy']), int(nb.loc[c, 'base_size']), mpp, click, gt_eval)
                ref = top[(top.file_name == fn) & (top.seed_index == s) & (top.condition == c)].sort_values('rank')
                h30 = det.head(30)
                same_xy = bool(np.array_equal(h30[['cx', 'cy']].to_numpy(), ref[['cx', 'cy']].to_numpy()))
                od_rel = float(np.max(np.abs(h30.od.to_numpy() - ref.od.to_numpy()) / np.abs(ref.od.to_numpy()))) if same_xy else float('nan')
                pipe_rows.append(dict(file_name=fn, seed_index=s, condition=c, my_n_peaks=info['n_peaks'], nb_n_peaks=int(nb.loc[c, 'n_peaks']), my_n_detections=info['n_detections'], nb_n_detections=int(nb.loc[c, 'n_detections']), my_n_self_hits=info['n_self_hits'], nb_n_self_hits=int(nb.loc[c, 'n_self_hits']), top30_xy_identical=same_xy, top30_od_max_rel_diff=od_rel, **{f'my_tp_{k}': int(det.hit.head(k).sum()) for k in KS}, **{f'nb_tp_{k}': int(nb.loc[c, f'tp_at_{k}']) for k in KS}, tp_full_list=int(det.hit.sum()), tp_in_top30_by_tm_score=int(det[det.score_rank < 30].hit.sum()), n_mitotic_eval=int((gt_eval.category == 1).sum()), t_reimpl_s=round(time.time() - t0, 2)))
        del rgb, gray, hem
        print(f'  [tierB] {fn} done ({time.time() - t_start:.0f}s)', flush=True)
    DR = pd.DataFrame(draw_rows)
    compare('tierB', 'click draw: seed_ann_id (147 clicks x 3 conditions)', DR.my_ann, DR.nb_ann, 'independent')
    compare('tierB', 'template base_size (441)', DR.my_base, DR.nb_base, 'independent')
    compare('tierB', 'template centre x,y (882)', DR[['my_tx', 'my_ty']].to_numpy(), DR[['nb_tx', 'nb_ty']].to_numpy(), 'independent', 1e-9)
    compare('tierB', 'refused_draws strings (441)', DR.my_refused, DR.nb_refused, 'independent')
    compare('tierB', 'n_joint_valid (441)', DR.my_joint_valid, DR.nb_joint_valid, 'independent')
    d0 = DR[DR.condition == 'default_51']
    prose(10, 'default_51 alone would have clicked differently (re-derived from pixels)', 26, int(d0.my_default_alone_differs.sum()))
    prose(10, 'production gray_bbox alone would have clicked differently (re-derived from pixels)', 10, int(d0.my_gray_alone_differs.sum()))
    PD = pd.DataFrame(prod_draw_rows)
    compare('tierB', 'production build_seed(gray) ann vs my gray-alone walk (147 clicks, earlier clicks fixed)', PD.prod_gray_ann, PD.my_gray_alone_ann, 'consistency')
    out['draws'] = DR
    out['prod_draws'] = PD
    P = pd.DataFrame(pipe_rows)
    for k in KS:
        compare('tierB', f're-implemented pipeline tp_at_{k} (15 runs)', P[f'my_tp_{k}'], P[f'nb_tp_{k}'], 'independent')
    compare('tierB', 're-implemented pipeline n_peaks / n_detections / n_self_hits (45 values)', P[['my_n_peaks', 'my_n_detections', 'my_n_self_hits']].to_numpy(), P[['nb_n_peaks', 'nb_n_detections', 'nb_n_self_hits']].to_numpy(), 'independent')
    compare('tierB', 're-implemented pipeline top-30 coordinates identical (15 runs)', P.top30_xy_identical.to_numpy(), np.ones(len(P), dtype=bool), 'independent')
    out['pipeline_reimpl'] = P
    # own-draw policy arm: substitute each condition's own-gate click where it differs from the joint click
    SR = pd.DataFrame(sens_rows)
    out['own_draw_runs'] = SR
    jt = pr[['file_name', 'seed_index', 'condition', 'tp_at_10', 'tp_at_20', 'tp_at_30']]
    vj = SR[SR.arm == 'own_draw'].merge(jt, on=['file_name', 'seed_index', 'condition'], validate='1:1')
    for k in KS:
        vj[f'own_minus_joint_tp{k}'] = vj[f'tp{k}'] - vj[f'tp_at_{k}']
    out['own_draw_vs_joint_by_condition'] = vj.groupby('condition').agg(n_clicks=('file_name', 'size'), **{f'sum_own_minus_joint_tp{k}': (f'own_minus_joint_tp{k}', 'sum') for k in KS}, **{f'mean_own_minus_joint_tp{k}': (f'own_minus_joint_tp{k}', 'mean') for k in KS}).reset_index()
    prose(10, 'hem_bbox alone would have clicked differently (re-derived from pixels; not claimed by the notebook)', 'n/a', int(d0.my_hem_alone_differs.sum()))
    base = T['per_run'].pivot_table(index=['subset', 'file_name', 'seed_index'], columns='condition', values=[f'tp_at_{k}' for k in KS])
    base.columns = [f'{v}|{c}' for v, c in base.columns]
    base = base.reset_index()
    s_rows = []
    for arm in ('own_draw', 'gray_gate'):
        sub_runs = SR[SR.arm == arm]
        own = base.copy()
        for r in sub_runs.itertuples():
            msk = (own.file_name == r.file_name) & (own.seed_index == r.seed_index)
            for k in KS:
                own.loc[msk, f'tp_at_{k}|{r.condition}'] = getattr(r, f'tp{k}')
        pairs = PAIRS if arm == 'own_draw' else (('gray_bbox', 'default_51'),)
        for fam in ('all_49', 'testing_35'):
            famp, start = [], len(s_rows)
            for a, b in PAIRS:
                for k in KS:
                    oo = subset_rows(own if (a, b) in pairs else base, fam)  # the gray-gated arm leaves the hem pairs at their joint values
                    d = (oo[f'tp_at_{k}|{a}'] - oo[f'tp_at_{k}|{b}']).to_numpy()
                    roi = pd.Series(d, index=oo.file_name.to_numpy()).groupby(level=0).sum()
                    p, nz = signflip_exact(roi.to_numpy())
                    famp.append(p)
                    s_rows.append(dict(arm=arm, family=fam, pair=f'{a} - {b}', K=k, n_clicks_substituted=int(sub_runs[sub_runs.file_name.isin(oo.file_name)].groupby(['file_name', 'seed_index']).ngroups), pooled_delta_pts=100 * roi.sum() / (len(roi) * 3 * k), wlt=f'{int((roi > 0).sum())}/{int((roi < 0).sum())}/{int((roi == 0).sum())}', exact_p=p, pair_meaningful_for_arm=(a, b) in pairs))
            hp = holm(famp)  # Holm over the same 9 tests as the notebook, with the substituted values
            for i, r in enumerate(s_rows[start:]):
                r['holm_p'] = hp[i]
    out['own_draw_sensitivity'] = pd.DataFrame(s_rows)
    # re-timing of one production call, three repeats, same ROI
    rgb = load_rgb('images/extra_valid/300.tiff')
    row = pr[(pr.file_name == '300.tiff') & (pr.seed_index == 0) & (pr.condition == 'gray_bbox')].iloc[0]
    sd = ss.Seed(ann_id=int(row.seed_ann_id), click_xy=(row.click_cx, row.click_cy), template_xy=(row.tpl_cx, row.tpl_cy), base_size=int(row.base_size), recentred=True, offset_px=0.0, n_retries=0, agreement_flagged=False, n_agreement_pool=0, n_after_border=0)
    times = []
    for _ in range(3):
        t0 = time.time()
        prod.run_production_pipeline(rgb, sd, float(row.mpp), rank_key='chromatin_od')
        times.append(time.time() - t0)
    t0 = time.time()
    _ = rgb2hed(rgb.astype(np.float32) / 255.0)
    t_hem = time.time() - t0
    out['timing'] = pd.DataFrame([dict(run='300.tiff s0 gray_bbox', notebook_t_pipeline_s=row.t_pipeline_s, audit_repeat_1_s=times[0], audit_repeat_2_s=times[1], audit_repeat_3_s=times[2], audit_rgb2hed_s=t_hem, notebook_median_all_runs_s=float(pr.t_pipeline_s.median()), notebook_sum_t_pipeline_min=float(pr.t_pipeline_s.sum() / 60))])
    print(f'[tierB] finished in {(time.time() - t_start) / 60:.1f} min')
    return out


# ----------------------------------------------------------------------------------------------
def main():
    """
        Run every audit section and write every table.

        Returns None.
    """
    T = {n: pd.read_csv(f'{STEM}_{n}.csv', float_precision='round_trip') for n in ('per_run', 'per_roi', 'summary', 'by_domain', 'by_seed', 'delta_per_roi', 'delta_stats', 'bootstrap_ci', 'contested_sensitivity', 'top30', 'verification')}
    save(execution_gate(), 'execution_gate')
    save(provenance(), 'provenance')
    imgs, anns = load_db()
    save(composition(T, imgs), 'composition')
    R, D = rescore(T, imgs, anns)
    save(R, 'rescore_per_run')
    pr, DP = tables(T, R)
    zero_rois = sorted(fn for fn, g in pr.groupby('file_name') if (g.tp30 == 0).all())
    print(f'[inference] ROIs with 0 TP@30 on all 9 runs: {zero_rois}')
    I = inference(T, DP, zero_rois)
    I, BS = bootstrap(DP, I)
    figure_checks(T, I, pr, DP)
    save(I, 'inference')
    save(BS.groupby(['group', 'pair', 'K']).agg(share_excludes_0=('excludes_0', 'mean'), ci_high_min=('ci_high', 'min'), ci_high_max=('ci_high', 'max'), ci_low_min=('ci_low', 'min'), ci_low_max=('ci_low', 'max')).reset_index(), 'bootstrap_seed_sensitivity')
    save(contested(T, pr, DP), 'contested')
    geo, GB = geometry_and_gates(T, D)
    save(GB, 'gate_b')
    closing_claims(T, I, pr, DP)
    ext = extensions(T, pr, DP, I)
    for n, f in ext.items():
        save(f, f'ext_{n}')
    save(precision_null(T, imgs, anns), 'ext_precision_null')
    z = D[D.file_name == '211.tiff']
    save(z.groupby(['seed_index', 'condition']).d_nearest_mitotic_px.describe().reset_index(), 'ext_211_nearest_mitosis')
    tb = tier_b(T, imgs, anns)
    for n, f in tb.items():
        save(f, f'tierb_{n}')
    C = pd.DataFrame(COMPARISONS)
    save(C, 'tier_a_b_comparisons')
    P = pd.DataFrame(PROSE)
    save(P, 'prose_claims')
    print(f"\nTOTAL: {int(C.n_compared.sum())} values compared, {int(C.n_divergent.sum())} divergent; by independence: {C.groupby('independence')[['n_compared', 'n_divergent']].sum().to_dict('index')}")
    print(f'prose/printed claims: {int(P.match.sum())}/{len(P)} match; mismatches:')
    print(P[~P.match].to_string(index=False))
    print(C[C.n_divergent > 0].to_string(index=False))
    print(f'audit ran in {(time.time() - T0) / 60:.1f} min')


if __name__ == '__main__':
    main()
