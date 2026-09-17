"""
    Re-run the 49-ROI x 3-seed bbox-refinement experiment against post-D11 production code.

    Replicates the measurement logic of
    production_hematoxylin_only/bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb
    (cells 1-4 and 12) by calling production.run_production_pipeline directly, and writes to
    *_postD11.csv so the committed pre-D11 reference CSVs are never touched. Run from the
    repo root with /Users/mohinianand/anaconda3/bin/python3.

    Modes:
      --select-only   run the joint click gate only, no pipeline calls (pre-flight)
      --rois A,B,C    restrict to these file names
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import os
import sys
import time

import cv2
import numpy as np
import pandas as pd

sys.path.insert(0, '.')
from midog_utils import channels as ch
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from midog_utils import invariants as inv
from midog_utils import production as prod
from midog_utils import seed_selection as ss
from midog_utils import template_match as tm

pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 80)

SUBSET_DIRS = {'original_14': 'images/extra_valid', 'testing_35': 'images/extra_valid/testing_set'}
SUBSET_SIZES = {'original_14': 14, 'testing_35': 35}
N_DOMAINS = 7
N_PER_DOMAIN = 7
DB_PATH = 'databases/MIDOG++.json'
SEED_INDICES = (0, 1, 2)

CONDITIONS = {'default_51': None, 'gray_bbox': 'gray_inverted', 'hem_bbox': 'hematoxylin_od'}
PAIRS = [('gray_bbox', 'default_51'), ('hem_bbox', 'default_51'), ('hem_bbox', 'gray_bbox')]

RANK_KEY = 'chromatin_od'
BUDGETS = (10, 20, 30)

DEFAULT_BASE_SIZE = tm.BASE_SIZE    # 51
OTSU_WINDOW = tm.BASE_SIZE          # 51
PATCH_SIZE = tm.PATCH_SIZE          # 73
BORDER = PATCH_SIZE // 2            # 36

# The notebook's REFERENCE_CONFIG minus SELF_HIT_RADIUS, which D11 removed. Everything else must still hold.
REFERENCE_CONFIG = dict(CHANNEL='hematoxylin_od', TM_METHOD=cv2.TM_CCOEFF, PEAK_MIN_DISTANCE=7, DEEP_FLOOR_Z=-1.5, MAX_PEAKS=100, OD_WINDOW=51)
MAX_PEAKS = prod.MAX_PEAKS

REF_STEM = 'results/precision_at_k_49roi_3seed_chromatin_bbox3way'
OUT_STEM = f'{REF_STEM}_postD11'
OUT = {name: f'{OUT_STEM}_{name}.csv' for name in ('per_run', 'per_roi', 'top30', 'selection_gate')}
assert not any(os.path.abspath(p) == os.path.abspath(f'{REF_STEM}_{n}.csv') for n, p in OUT.items()), 'an output path would overwrite the committed pre-D11 run'


def roi_table(subset_dirs=SUBSET_DIRS):
    """
        Every ROI TIFF in each subset directory, listed non-recursively.

        subset_dirs (dict): subset name -> directory holding that subset's TIFFs.

        Returns pd.DataFrame: one row per file, with subset, file_name and path.
    """
    return pd.DataFrame([dict(subset=subset, file_name=fn, path=f'{d}/{fn}') for subset, d in subset_dirs.items() for fn in sorted(f for f in os.listdir(d) if f.endswith('.tiff'))])


def template_spec(condition, channels, cx, cy):
    """
        The template size and centre one condition cuts for one click.

        condition (str): a key of CONDITIONS.
        channels (dict): channel name -> this ROI's converted single-channel image.
        cx (float): click x, full-image pixels.
        cy (float): click y, full-image pixels.

        Returns tuple[int, float, float] or None: (base_size, template x, template y), or None if the gate refuses.
    """
    tighten_channel = CONDITIONS[condition]
    spec = (DEFAULT_BASE_SIZE, cx, cy) if tighten_channel is None else ss.tightened_template_box(channels[tighten_channel], cx, cy, otsu_window=OTSU_WINDOW)
    if spec is None or tm.read_padded_patch(channels['hematoxylin_od'], spec[1], spec[2], PATCH_SIZE) is None:
        return None
    return int(spec[0]), float(spec[1]), float(spec[2])


def template_digest(hem, base_size, tx, ty):
    """
        A digest of the exact template pixels find_and_suppress cuts for one condition.

        hem (np.ndarray): this ROI's hematoxylin_od channel, the search channel.
        base_size (int): the template side length.
        tx (float): template centre x.
        ty (float): template centre y.

        Returns str: SHA-1 of the float32 template bytes.
    """
    templates, _ = tm.build_augmentations(tm.read_padded_patch(hem, tx, ty, PATCH_SIZE), base_size, (1.0,), 1, (False,))
    return hashlib.sha1(np.ascontiguousarray(templates[0], dtype=np.float32).tobytes()).hexdigest()


def pool_gate_specs(pool, channels):
    """
        Every condition's template spec for every annotation in the border-filtered pool.

        pool (pd.DataFrame): agreement-tier, border-filtered mitotic annotations.
        channels (dict): channel name -> this ROI's converted single-channel image.

        Returns dict: ann_id -> {condition: (base_size, x, y) or None}.
    """
    return {int(row['ann_id']): {condition: template_spec(condition, channels, float(row['cx']), float(row['cy'])) for condition in CONDITIONS} for _, row in pool.iterrows()}


def draw_joint_click(pool, rng, gate_specs):
    """
        Draw one click every condition can build a template from, redrawing on the same RNG stream after a refusal.

        pool (pd.DataFrame): the border-filtered pool, minus the clicks earlier seeds accepted.
        rng (np.random.Generator): this seed's [seed_index, image_id] stream.
        gate_specs (dict): ann_id -> condition -> spec, from pool_gate_specs.

        Returns tuple[pd.Series, dict, list[str]]: the accepted annotation, condition -> spec, and one refusal entry per refused draw.
    """
    working, refusals = pool.copy(), []
    while len(working) > 0:
        idx = int(rng.integers(len(working)))
        row = working.iloc[idx]
        specs = gate_specs[int(row['ann_id'])]
        if all(spec is not None for spec in specs.values()):
            return row, specs, refusals
        refusals.append(f"{int(row['ann_id'])}:{'+'.join(c for c, spec in specs.items() if spec is None)}")
        working = working.drop(working.index[idx])
    raise ValueError('seed pool exhausted -- no annotation passed every condition')


def score_top_k(detections, gt_eval, match_radius):
    """
        Match a ranked detection list to ground truth and count true positives in each top-K.

        detections (pd.DataFrame): ranked detections, best first.
        gt_eval (pd.DataFrame): ground truth with the click's own annotation removed.
        match_radius (float): match radius in pixels.

        Returns tuple[pd.DataFrame, dict]: the detections with match buckets added, and K -> true positives in the top K.
    """
    det_out, _ = ev.bucket_detections(detections, gt_eval, match_radius)
    hit = (det_out['bucket'] == ev.HUMAN_CORRECT_LABEL).to_numpy()
    return det_out, {k: int(hit[:k].sum()) for k in BUDGETS}


def select_clicks(roi, annotations):
    """
        The joint-gate click selection for one ROI -- the half of run_roi that runs before any pipeline call.

        roi (pd.Series): one ROIS row (subset, file_name, path, image_id, domain).
        annotations (pd.DataFrame): the full MIDOG++ annotation table.

        Returns tuple: (rgb, mpp, gt, pool, flagged, n_agreement, n_joint_valid, draws) where draws is one (seed_index, ann, specs, refusals, n_seed_pool, digests) per seed.
    """
    fn, image_id = roi['file_name'], int(roi['image_id'])
    rgb = ds.load_roi(roi['path'])
    mpp = ds.roi_mpp(roi['path'])
    channels = {name: ch.to_channel(rgb, name) for name in ('gray_inverted', 'hematoxylin_od')}

    gt = ds.image_annotations(annotations, fn)
    pool, flagged = ss.agreement_pool(gt[gt['category_id'] == ds.MITOTIC])
    n_agreement = len(pool)
    pool = ss.border_filter(pool, BORDER, rgb.shape)
    gate_specs = pool_gate_specs(pool, channels)
    n_joint_valid = sum(all(spec is not None for spec in specs.values()) for specs in gate_specs.values())
    if n_joint_valid < len(SEED_INDICES):
        raise ValueError(f'{fn}: only {n_joint_valid} pool annotations pass the joint gate, fewer than the {len(SEED_INDICES)} distinct clicks needed')

    draws, accepted = [], []
    for seed_index in SEED_INDICES:
        seed_pool = pool[~pool['ann_id'].isin(accepted)] if accepted else pool
        rng = np.random.default_rng([seed_index, image_id])
        ann, specs, refusals = draw_joint_click(seed_pool, rng, gate_specs)
        digests = {condition: template_digest(channels['hematoxylin_od'], *spec) for condition, spec in specs.items()}
        draws.append((seed_index, ann, specs, refusals, len(seed_pool), digests))
        accepted.append(int(ann['ann_id']))
    del channels
    gc.collect()
    return rgb, mpp, gt, pool, flagged, n_agreement, n_joint_valid, draws


def selection_rows(roi, draws, pool, flagged, n_agreement, n_joint_valid):
    """
        One row per (seed, condition) describing the click and template that selection chose, with no pipeline output.

        roi (pd.Series): one ROIS row.
        draws (list): the draws select_clicks returned.
        pool (pd.DataFrame): the border-filtered pool.
        flagged (bool): whether the agreement pool fell back to the contested tier.
        n_agreement (int): agreement-tier pool size.
        n_joint_valid (int): pool annotations passing the joint gate.

        Returns list[dict]: the selection rows.
    """
    out = []
    for seed_index, ann, specs, refusals, n_seed_pool, digests in draws:
        for condition, (base_size, tx, ty) in specs.items():
            offset = float(np.hypot(tx - float(ann['cx']), ty - float(ann['cy'])))
            out.append(dict(subset=roi['subset'], file_name=roi['file_name'], domain=roi['domain'], image_id=int(roi['image_id']), seed_index=seed_index, condition=condition, seed_ann_id=int(ann['ann_id']), click_cx=float(ann['cx']), click_cy=float(ann['cy']), base_size=int(base_size), tpl_cx=tx, tpl_cy=ty, tpl_offset_px=round(offset, 3), tpl_sha1=digests[condition], n_agreement_pool=n_agreement, n_after_border=len(pool), n_joint_valid=n_joint_valid, n_seed_pool=n_seed_pool, n_retries=len(refusals), refused_draws=';'.join(refusals), contested_seed_tier=bool(flagged)))
    return out


def run_roi(roi, annotations):
    """
        Draw three distinct joint clicks on one ROI, then run the production pipeline once per (click, condition).

        roi (pd.Series): one ROIS row (subset, file_name, path, image_id, domain).
        annotations (pd.DataFrame): the full MIDOG++ annotation table.

        Returns tuple[list[dict], list[pd.DataFrame]]: one result row per (seed, condition), and each run's top-max(BUDGETS) scored detections.
    """
    fn = roi['file_name']
    rgb, mpp, gt, pool, flagged, n_agreement, n_joint_valid, draws = select_clicks(roi, annotations)

    match_radius = ev.radius_px(mpp)
    rows, tops = [], []
    for seed_index, ann, specs, refusals, n_seed_pool, digests in draws:
        click_xy = (float(ann['cx']), float(ann['cy']))
        seed_ann_id = int(ann['ann_id'])
        gt_eval = gt[gt['ann_id'] != seed_ann_id].reset_index(drop=True)
        n_gt = int((gt_eval['category_id'] == ds.MITOTIC).sum())
        nearest_other = float(np.hypot(gt_eval['cx'] - click_xy[0], gt_eval['cy'] - click_xy[1]).min())
        line = []
        for condition, (base_size, tx, ty) in specs.items():
            offset = float(np.hypot(tx - click_xy[0], ty - click_xy[1]))
            seed = ss.Seed(ann_id=seed_ann_id, click_xy=click_xy, template_xy=(tx, ty), base_size=base_size, recentred=CONDITIONS[condition] is not None, offset_px=offset, n_retries=len(refusals), agreement_flagged=bool(flagged), n_agreement_pool=n_agreement, n_after_border=len(pool))
            t0 = time.time()
            detections, info = prod.run_production_pipeline(rgb, seed, mpp, rank_key=RANK_KEY)
            det_out, tp = score_top_k(detections, gt_eval, match_radius)

            d_click = np.hypot(det_out['cx'] - click_xy[0], det_out['cy'] - click_xy[1]).to_numpy()
            d_tpl = np.hypot(det_out['cx'] - tx, det_out['cy'] - ty).to_numpy()
            near = np.flatnonzero(d_click <= match_radius)
            near_detail = ';'.join(f"rank{int(det_out['rank'].iat[i])}:{det_out['bucket'].iat[i]}:ann{int(det_out['matched_ann_id'].iat[i])}:d_click={d_click[i]:.1f}:d_tpl={d_tpl[i]:.1f}" for i in near)
            od = det_out['od'].to_numpy()
            rows.append(dict(subset=roi['subset'], file_name=fn, domain=roi['domain'], image_id=int(roi['image_id']), seed_index=seed_index, condition=condition, seed_ann_id=seed_ann_id, click_cx=click_xy[0], click_cy=click_xy[1], base_size=int(base_size), tpl_cx=tx, tpl_cy=ty, tpl_offset_px=round(offset, 3), tpl_sha1=digests[condition], n_agreement_pool=n_agreement, n_after_border=len(pool), n_joint_valid=n_joint_valid, n_seed_pool=n_seed_pool, n_retries=len(refusals), refused_draws=';'.join(refusals), contested_seed_tier=bool(flagged), n_gt_mitotic=n_gt, nearest_other_ann_px=round(nearest_other, 2), n_peaks=int(info['n_peaks']), max_peaks_binding=bool(info['max_peaks_binding']), n_detections=int(info['n_detections']), n_blanked_px=int(info['n_blanked_px']), n_near_click=len(near), near_click_detail=near_detail, od_nan=int(np.isnan(od).sum()), od_descending=bool(np.all(np.diff(od) <= 0)), match_radius_px=match_radius, mpp=mpp, t_pipeline_s=round(time.time() - t0, 1), **{f'tp_at_{k}': tp[k] for k in BUDGETS}))
            tops.append(det_out.head(max(BUDGETS)).assign(subset=roi['subset'], file_name=fn, seed_index=seed_index, condition=condition))
            line.append(f"{condition} base={base_size:2d} off={offset:4.1f} n_det={info['n_detections']:3d} blank={info['n_blanked_px']:5d} tp={tuple(tp.values())}")
        print(f"  [{fn} s{seed_index}] ann={seed_ann_id} retries={len(refusals)} | " + ' | '.join(line), flush=True)

    del rgb
    gc.collect()
    return rows, tops


SELECTION_COLS = ['seed_ann_id', 'base_size', 'tpl_cx', 'tpl_cy', 'tpl_offset_px', 'tpl_sha1', 'n_agreement_pool', 'n_after_border', 'n_joint_valid', 'n_seed_pool', 'n_retries', 'refused_draws']


def compare_selection(mine, ref):
    """
        Check that the joint click gate picked the same clicks and cut the same templates as the stored pre-D11 run.

        mine (pd.DataFrame): selection or per-run rows from this run.
        ref (pd.DataFrame): the stored pre-D11 per_run.csv.

        Returns pd.DataFrame: one row per (file_name, seed_index, condition) with a per-column match flag and a mismatch detail string.
    """
    key = ['file_name', 'seed_index', 'condition']
    ref = ref[ref['file_name'].isin(set(mine['file_name']))]  # a partial run compares only against its own ROIs
    merged = mine[key + SELECTION_COLS].merge(ref[key + SELECTION_COLS], on=key, how='outer', suffixes=('', '_ref'), indicator=True)
    assert (merged['_merge'] == 'both').all(), f"selection: {int((merged['_merge'] != 'both').sum())} rows unmatched against the reference"
    merged[['refused_draws', 'refused_draws_ref']] = merged[['refused_draws', 'refused_draws_ref']].fillna('')  # an empty refusal list reads back from CSV as NaN
    detail = []
    for _, r in merged.iterrows():
        bad = [c for c in SELECTION_COLS if not (pd.isna(r[c]) and pd.isna(r[f'{c}_ref'])) and str(r[c]) != str(r[f'{c}_ref'])]
        detail.append(';'.join(f'{c}={r[c]!r}!={r[f"{c}_ref"]!r}' for c in bad))
    merged['detail'] = detail
    merged['passed'] = [d == '' for d in detail]
    return merged.drop(columns='_merge')


def build_rois():
    """
        The 49 ROIs, ordered and checked the way the notebook orders and checks them.

        Returns tuple[pd.DataFrame, pd.DataFrame]: the ROI table with image_id and domain, and the full annotation table.
    """
    images, annotations = ds.load_annotations(DB_PATH)
    assert ((annotations['w'] == ds.BOX_SIZE) & (annotations['h'] == ds.BOX_SIZE)).all(), f'an annotation box is not {ds.BOX_SIZE}x{ds.BOX_SIZE}'
    assert set(annotations['category_id'].unique()) <= {ds.MITOTIC, ds.LOOKALIKE}, 'unexpected annotation category id'
    assert not images['file_name'].duplicated().any(), 'a file name appears twice in the annotation DB'
    meta_ix = images.set_index('file_name')[['image_id', 'tumor_type']]

    rois = roi_table()
    missing = sorted(set(rois['file_name']) - set(meta_ix.index))
    assert not missing, f'.tiff on disk absent from the annotation DB: {missing}'
    assert not rois['file_name'].duplicated().any(), 'file name in both subsets'
    assert rois['subset'].value_counts().to_dict() == SUBSET_SIZES, f"expected {SUBSET_SIZES}, found {rois['subset'].value_counts().to_dict()}"
    rois['image_id'] = rois['file_name'].map(meta_ix['image_id']).astype(int)
    rois['domain'] = rois['file_name'].map(meta_ix['tumor_type'])
    per_domain = rois.groupby('domain').size()
    assert len(rois) == sum(SUBSET_SIZES.values()) == N_DOMAINS * N_PER_DOMAIN and len(per_domain) == N_DOMAINS and (per_domain == N_PER_DOMAIN).all(), f'expected {N_PER_DOMAIN} ROIs in each of {N_DOMAINS} domains:\n{per_domain}'
    return rois.sort_values(['domain', 'subset', 'file_name']).reset_index(drop=True), annotations


def main():
    """
        Run the experiment (or just its click selection) and write the post-D11 CSVs.

        Returns None.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--select-only', action='store_true', help='run the joint click gate only, no pipeline calls')
    parser.add_argument('--rois', default='', help='comma-separated file names to restrict to')
    args = parser.parse_args()

    drift = pd.DataFrame([dict(setting=k, production_py=getattr(prod, k), reference_notebook=v, match=getattr(prod, k) == v) for k, v in REFERENCE_CONFIG.items()])
    print(drift.to_string(index=False))
    assert drift['match'].all(), 'midog_utils.production no longer matches the configuration the notebook ran at'
    assert ev.MIDOG_RADIUS_UM == 7.5, 'D7: NMS radius = match radius = 7.5 um'
    assert not hasattr(prod, 'SELF_HIT_RADIUS'), 'production.SELF_HIT_RADIUS should be gone post-D11'
    print(f'\npost-D11 re-run: MAX_PEAKS = {MAX_PEAKS}, rank_key = {RANK_KEY}, K = {BUDGETS}\n')

    rois, annotations = build_rois()
    if args.rois:
        wanted = [s.strip() for s in args.rois.split(',') if s.strip()]
        rois = rois[rois['file_name'].isin(wanted)].reset_index(drop=True)
        assert len(rois) == len(wanted), f'not all requested ROIs found: {sorted(set(wanted) - set(rois["file_name"]))}'
    ref = pd.read_csv(f'{REF_STEM}_per_run.csv', float_precision='round_trip')

    t_run = time.time()
    if args.select_only:
        sel_rows = []
        for i, roi in rois.iterrows():
            rgb, mpp, gt, pool, flagged, n_agreement, n_joint_valid, draws = select_clicks(roi, annotations)
            del rgb
            gc.collect()
            sel_rows.extend(selection_rows(roi, draws, pool, flagged, n_agreement, n_joint_valid))
            print(f"[{i + 1}/{len(rois)} {roi['file_name']}] joint-valid pool {n_joint_valid}; clicks {[int(d[1]['ann_id']) for d in draws]} [{(time.time() - t_run) / 60:.1f} min]", flush=True)
        sel = pd.DataFrame(sel_rows)
        cmp = compare_selection(sel, ref)
        cmp.to_csv(OUT['selection_gate'], index=False)
        print(f"\nselection gate: {int(cmp['passed'].sum())}/{len(cmp)} (ROI, seed, condition) rows reproduce the stored pre-D11 click and template exactly -> {OUT['selection_gate']}")
        if not cmp['passed'].all():
            print(cmp.loc[~cmp['passed'], ['file_name', 'seed_index', 'condition', 'detail']].to_string(index=False))
        return

    result_rows, top_frames = [], []
    for i, roi in rois.iterrows():
        t_roi = time.time()
        rows, tops = run_roi(roi, annotations)
        result_rows.extend(rows)
        top_frames.extend(tops)
        frame = pd.DataFrame(result_rows)
        frame.to_csv(OUT['per_run'], index=False)
        pd.concat(top_frames, ignore_index=True).to_csv(OUT['top30'], index=False)
        # fail fast: this ROI's clicks and templates must match the stored pre-D11 run
        cmp_roi = compare_selection(frame[frame['file_name'] == roi['file_name']], ref[ref['file_name'] == roi['file_name']])
        assert cmp_roi['passed'].all(), f"selection diverged on {roi['file_name']}:\n{cmp_roi.loc[~cmp_roi['passed'], ['seed_index', 'condition', 'detail']].to_string(index=False)}"
        print(f"[{i + 1}/{len(rois)} {roi['file_name']}] {roi['subset']}, {roi['domain']}: joint-valid pool {rows[0]['n_joint_valid']}, selection matches reference [{time.time() - t_roi:.0f}s; {(time.time() - t_run) / 60:.1f} min so far]", flush=True)

    results = pd.DataFrame(result_rows)
    for k in BUDGETS:
        results[f'precision_at_{k}'] = results[f'tp_at_{k}'] / k
    roi_order = rois['file_name'].tolist()
    results['_roi_o'] = results['file_name'].map({fn: i for i, fn in enumerate(roi_order)})
    results['_cond_o'] = results['condition'].map({c: i for i, c in enumerate(CONDITIONS)})
    results = results.sort_values(['_roi_o', 'seed_index', '_cond_o']).drop(columns=['_roi_o', '_cond_o']).reset_index(drop=True)
    results.to_csv(OUT['per_run'], index=False)

    per_roi = results.groupby(['subset', 'domain', 'file_name', 'condition'], sort=False).agg(n_seeds=('seed_index', 'nunique'), seed_ann_ids=('seed_ann_id', lambda v: ';'.join(map(str, v))), base_sizes=('base_size', lambda v: ';'.join(map(str, v))), tpl_offsets_px=('tpl_offset_px', lambda v: ';'.join(f'{x:g}' for x in v)), min_n_detections=('n_detections', 'min'), **{f'tp_at_{k}': (f'tp_at_{k}', 'sum') for k in BUDGETS}).reset_index()
    for k in BUDGETS:
        per_roi[f'precision_at_{k}'] = per_roi[f'tp_at_{k}'] / (per_roi['n_seeds'] * k)
    per_roi.to_csv(OUT['per_roi'], index=False)

    cmp = compare_selection(results, ref)
    cmp.to_csv(OUT['selection_gate'], index=False)
    print(f"\n{len(results)} pipeline runs in {(time.time() - t_run) / 60:.1f} min -> {OUT['per_run']}, {OUT['per_roi']}, {OUT['top30']}")
    print(f"selection gate: {int(cmp['passed'].sum())}/{len(cmp)} rows reproduce the stored pre-D11 click and template exactly")

    halt = dict(budget_delivered=int((results['n_detections'] >= max(BUDGETS)).sum()), od_no_nan=int((results['od_nan'] == 0).sum()), od_descending=int(results['od_descending'].sum()))
    print(f"halt checks (of {len(results)}): {halt}")
    print(f"n_blanked_px per run: {results['n_blanked_px'].value_counts().sort_index().to_dict()}")
    print(f"n_peaks per run: {results['n_peaks'].value_counts().sort_index().to_dict()}; n_detections min {results['n_detections'].min()}, max {results['n_detections'].max()}")


if __name__ == '__main__':
    main()
