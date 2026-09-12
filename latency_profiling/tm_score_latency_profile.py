#!/usr/bin/env python3
"""Wall-clock latency profile of the tm_score-only production pipeline, per stage,
across the 14 extra_valid ROIs. Disposable profiling script -- not a durable artifact.

Reuses, verbatim, the config and seed-draw logic of the accepted production notebook:
    production_seed_precision_at_k/production_seed_precision_at_k_chromatin_half_pix_fix.ipynb
(cells 1, 3, 5, 6; D8_TEMPLATE_ANCHOR.md, current as of 2026-09-10). Pure timing -- no
precision, no chromatin/shape features, no other ranking arms.

SETUP vs PIPELINE-stage-1 note: the notebook's own prose calls "drawing the GT row (the
click proxy)" test-harness setup, which could be read as pulling the whole seed draw out
of the timed region. But `draw_seed_with_retry` performs the RNG draw and the
`tightened_template_box` gate check as one inseparable retry loop -- splitting them would
break the RNG stream this script must match to stay comparable to the notebook. So SETUP
here stops at building `seed_pool` (loading the ROI, the channels, and the candidate
pool), and PIPELINE stage 1 is the entire `draw_seed_with_retry` call, retries included --
exactly as the "stage boundaries" section (the part called out to get right) specifies.
"""
import gc
import os
import sys
import time

import cv2
import numpy as np
import pandas as pd

sys.path.insert(0, '..')
from midog_utils import channels as ch
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from midog_utils import find_and_suppress as fs
from midog_utils import seed_selection as ss
from midog_utils import template_match as tm
from midog_utils.nms import nms_by_distance

pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 60)

# ---------------------------------------------------------------------------------------
# Config -- copied verbatim from cell 1 of the accepted notebook. Do not change.
# ---------------------------------------------------------------------------------------
IMAGES_DIR = '../images/extra_valid'
SEED_INDEX = 0

CHANNEL = 'hematoxylin_od'
METHOD = cv2.TM_CCOEFF
PEAK_MIN_DISTANCE = 7
SELF_HIT_RADIUS = 5.0
DEEP_FLOOR_Z = -1.5
MAX_PEAKS = 2_000_000

NMS_RADIUS_UM = ev.MIDOG_RADIUS_UM
MATCH_RADIUS_UM = ev.MIDOG_RADIUS_UM

CFG = fs.FSConfig(channel=CHANNEL, base_size=tm.BASE_SIZE, scales=(1.0,),
                  n_angles=1, flips=(False,), peak_min_distance=PEAK_MIN_DISTANCE,
                  self_hit_radius=SELF_HIT_RADIUS)
BORDER = CFG.patch_size // 2
OTSU_WINDOW = tm.BASE_SIZE

BUDGET = 100  # this task's top-K; the notebook's own BUDGETS=(10,20,30,50) don't apply

ORACLE_RAW_CSV = '../results/precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv'
OUT_PER_ROI = 'tm_score_latency_per_roi.csv'
OUT_SUMMARY = 'tm_score_latency_summary.csv'
OUT_BY_DOMAIN = 'tm_score_latency_by_domain.csv'


def draw_seed_with_retry(pool, rng, check_fn):
    """Draw a row via `rng.integers`; on failure drop it and redraw on the same stream.

    Copied verbatim from the notebook (cell 2).
    """
    working, retries = pool.copy(), 0
    while len(working) > 0:
        idx = int(rng.integers(len(working)))
        row = working.iloc[idx]
        result = check_fn(row)
        if result is not None:
            return row, result, retries
        working = working.drop(working.index[idx])
        retries += 1
    raise ValueError('seed pool exhausted -- no candidate passed check_fn')


def suppress(centers, scores, radius, ref_xy):
    """NMS at `radius`, then drop the template's own self-correlation.

    Reimplemented from the notebook's `suppress()` (cell 2) -- it is inline in the
    notebook, not part of any midog_utils module.
    """
    keep = nms_by_distance(centers, scores, radius)
    c, s = centers[keep], scores[keep]
    if len(c):
        ok = np.hypot(c[:, 0] - ref_xy[0], c[:, 1] - ref_xy[1]) > SELF_HIT_RADIUS
        c, s = c[ok], s[ok]
    return c, s


def roi_files(images_dir=IMAGES_DIR):
    return sorted(f for f in os.listdir(images_dir) if f.endswith('.tiff'))


def time_roi(fn, image_id, domain, anns):
    gc.disable()
    stages = {}

    # =================== SETUP (once per ROI; not part of click latency) ==============
    # Loads the ROI, builds the click-candidate pool, then -- test-harness only --
    # searches that pool (with retries) for a point tightened_template_box will accept.
    # A real click skips this search: the user clicks one point, already assumed valid.
    t0 = time.perf_counter()
    path = f'{IMAGES_DIR}/{fn}'
    rgb = ds.load_roi(path)
    mpp = ds.roi_mpp(path)
    roi_shape = rgb.shape
    nms_radius = ev.radius_px(mpp, NMS_RADIUS_UM)
    hem = ch.to_channel(rgb, CHANNEL)
    gray_inv = ch.to_gray_inverted(rgb)
    H, W = hem.shape[:2]
    del rgb

    gt = ds.image_annotations(anns, fn)
    seed_pool, flagged = ss.agreement_pool(gt[gt['category_id'] == ds.MITOTIC])
    seed_pool = ss.border_filter(seed_pool, BORDER, roi_shape)

    rng = np.random.default_rng([SEED_INDEX, image_id])

    def _check(row):
        r = ss.tightened_template_box(gray_inv, float(row['cx']), float(row['cy']),
                                      otsu_window=OTSU_WINDOW)
        if r is None:
            return None
        if tm.read_padded_patch(hem, r[1], r[2], CFG.patch_size) is None:
            return None
        return r

    seed, seed_tpl_probe, n_retries = draw_seed_with_retry(seed_pool, rng, _check)
    seed_ann_id = int(seed['ann_id'])
    click_cx, click_cy = float(seed['cx']), float(seed['cy'])
    t_setup = time.perf_counter() - t0

    # ======================= PIPELINE (click-to-top-100 latency) =======================
    # The box has been drawn: `seed` is already a *valid* click, established above.
    # Timing starts at the step that refines that valid seed's bounding box.
    t_outer0 = time.perf_counter()

    # ---- Stage 1: refine the bounding box of the valid seed (single call, no retry) ---
    t1 = time.perf_counter()
    r = ss.tightened_template_box(gray_inv, click_cx, click_cy, otsu_window=OTSU_WINDOW)
    assert r is not None, f'{fn}: seed was pre-validated in SETUP but refused here'
    border_probe = tm.read_padded_patch(hem, r[1], r[2], CFG.patch_size)
    assert border_probe is not None, f'{fn}: seed was pre-validated in SETUP but unreadable here'
    assert r == seed_tpl_probe, f'{fn}: single-shot refinement diverged from the SETUP search'
    base_size, tpl_cx, tpl_cy = r
    tpl_xy = (float(tpl_cx), float(tpl_cy))
    stages['t1_refine_seed_box_s'] = time.perf_counter() - t1

    # ---- Stage 2: patch + template build ----------------------------------------------
    t2 = time.perf_counter()
    patch = tm.read_padded_patch(hem, *tpl_xy, CFG.patch_size)
    templates, _ = tm.build_augmentations(patch, base_size, CFG.scales, CFG.n_angles, CFG.flips)
    stages['t2_patch_template_build_s'] = time.perf_counter() - t2

    # ---- Stage 3: template matching (suspected dominant cost) -------------------------
    t3 = time.perf_counter()
    PAD = max((t.shape[0] - 1) // 2 for t in templates)
    hem_p = cv2.copyMakeBorder(hem, PAD, PAD, PAD, PAD, borderType=cv2.BORDER_REPLICATE)
    fused_p, _, valid_p = tm.fused_response(hem_p, templates, CFG.scale_normalize, method=METHOD)
    stages['t3_template_matching_s'] = time.perf_counter() - t3

    # ---- Stage 4: threshold + peak extraction ------------------------------------------
    t4 = time.perf_counter()
    fused = fused_p[PAD:PAD + H, PAD:PAD + W]
    valid = valid_p[PAD:PAD + H, PAD:PAD + W]
    assert bool(valid.all()), f'{fn}: padding left part of the ROI unreachable (PAD={PAD})'
    med, mad = tm.robust_stats(fused, valid)
    cut = med + DEEP_FLOOR_Z * mad
    centers, scores = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE, cut, MAX_PEAKS)
    assert len(centers) < MAX_PEAKS, f'{fn}: MAX_PEAKS is binding, raise it'
    stages['t4_threshold_peak_extraction_s'] = time.perf_counter() - t4

    # ---- Stage 5: NMS + self-hit suppression -------------------------------------------
    t5 = time.perf_counter()
    c, s = suppress(centers, scores, nms_radius, tpl_xy)
    stages['t5_nms_selfhit_s'] = time.perf_counter() - t5
    assert bool(np.all(np.diff(s) <= 0)), f'{fn}: post-NMS pool is not score-descending'

    # ---- Stage 6: rank + top-100 (matches midog_utils/compare.py's _rank convention) ---
    t6 = time.perf_counter()
    pool = pd.DataFrame({'cx': c[:, 0], 'cy': c[:, 1], 'score': s})
    top = pool.sort_values('score', ascending=False, na_position='last',
                           kind='mergesort').head(BUDGET)
    stages['t6_rank_top100_s'] = time.perf_counter() - t6

    t_outer_total = time.perf_counter() - t_outer0
    gc.enable()

    n_detections = len(pool)
    n_top = len(top)
    stage_sum = sum(stages.values())

    del hem, gray_inv, hem_p, fused_p, valid_p, fused, valid, templates, patch, centers, scores, c, s
    gc.collect()

    meta = dict(
        file_name=fn, tumor_type=domain,
        n_detections=n_detections, base_size=base_size, n_retries=n_retries, n_top=n_top,
        seed_ann_id=seed_ann_id, map_median=round(float(med), 5), mad_scale=round(float(mad), 5),
        t_setup_s=round(t_setup, 5),
        **{k: round(v, 5) for k, v in stages.items()},
        t_stage_sum_s=round(stage_sum, 5),
        t_outer_total_s=round(t_outer_total, 5),
    )
    print(f"[{fn}] {domain:32s} base={base_size:2d} n_detections={n_detections:6d} "
          f"n_retries={n_retries} n_top={n_top:3d} "
          f"sum={stage_sum * 1000:7.1f}ms outer={t_outer_total * 1000:7.1f}ms "
          f"(setup={t_setup * 1000:.0f}ms)", flush=True)
    return meta


def cross_check_against_oracle(df):
    """These 14 ROIs should match `ORACLE_RAW_CSV` exactly: same RNG seed and config ->
    same template -> same fused map -> same post-NMS pool. Any mismatch means this
    profiling harness diverged from the accepted notebook's own pipeline -- a bug to
    report, not proceed past.
    """
    oracle = pd.read_csv(ORACLE_RAW_CSV)
    oracle_roi = (oracle.drop_duplicates('file_name')
                  .set_index('file_name')[['seed_ann_id', 'base_size', 'n_detections',
                                            'map_median', 'mad_scale']])
    cmp = df[['seed_ann_id', 'base_size', 'n_detections', 'map_median', 'mad_scale']].join(
        oracle_roi, lsuffix='_this', rsuffix='_oracle')

    mismatches = []
    for col in ['seed_ann_id', 'base_size', 'n_detections']:
        bad = cmp[f'{col}_this'].astype(int) != cmp[f'{col}_oracle'].astype(int)
        if bad.any():
            mismatches.append((col, cmp.index[bad].tolist()))
    for col in ['map_median', 'mad_scale']:
        bad = ~np.isclose(cmp[f'{col}_this'], cmp[f'{col}_oracle'], rtol=0, atol=1e-5)
        if bad.any():
            mismatches.append((col, cmp.index[bad].tolist()))

    if mismatches:
        print("\n!! MISMATCH vs oracle -- profiling harness diverged from production:")
        for col, rois in mismatches:
            print(f"   {col}: {rois}")
        print(cmp)
        raise AssertionError('profiling harness does not reproduce the accepted pipeline')
    print(f"\nall {len(cmp)} ROIs match {ORACLE_RAW_CSV} exactly on "
          f"seed_ann_id/base_size/n_detections/map_median/mad_scale -- harness reproduces "
          f"production.")


def main():
    images, annotations = ds.load_annotations('../databases/MIDOG++.json')
    ds.check_invariants(annotations)
    meta_ix = images.set_index('file_name')[['image_id', 'tumor_type']]

    files = roi_files()
    assert len(files) == 14, f'expected 14 ROIs in {IMAGES_DIR}/, found {len(files)}'

    rows = []
    for fn in files:
        image_id = int(meta_ix.loc[fn, 'image_id'])
        domain = meta_ix.loc[fn, 'tumor_type']
        rows.append(time_roi(fn, image_id, domain, annotations))
        gc.collect()

    df = pd.DataFrame(rows).set_index('file_name')
    df.to_csv(OUT_PER_ROI)
    print(f"\nwrote {OUT_PER_ROI}")
    print(df)

    cross_check_against_oracle(df)

    stage_cols = [c for c in df.columns if c.startswith('t') and c.endswith('_s')
                 and c not in ('t_setup_s',)]
    time_cols = ['t_setup_s'] + stage_cols

    summary_all = df[time_cols].agg(['mean', 'std', 'min', 'max']).T
    print("\n=== summary across all 14 ROIs (seconds) ===")
    print(summary_all)
    summary_all.to_csv(OUT_SUMMARY)
    print(f"wrote {OUT_SUMMARY}")

    by_domain = df.groupby('tumor_type')[time_cols].agg(['mean', 'std', 'min', 'max'])
    print("\n=== summary by tumor_type (seconds) ===")
    print(by_domain)
    by_domain.to_csv(OUT_BY_DOMAIN)
    print(f"wrote {OUT_BY_DOMAIN}")

    # cold-start check: OpenCV's parallel backend and file cache are cold on ROI #1
    print("\n=== stage 3 (template matching): ROI #1 (cold) vs remaining 13 ===")
    print(f"  ROI #1 ({df.index[0]}): {df['t3_template_matching_s'].iloc[0] * 1000:.1f} ms")
    rest = df['t3_template_matching_s'].iloc[1:]
    print(f"  remaining 13: mean={rest.mean() * 1000:.1f}ms std={rest.std() * 1000:.1f}ms "
          f"min={rest.min() * 1000:.1f}ms max={rest.max() * 1000:.1f}ms")

    starved = df[df['n_top'] < BUDGET]
    if len(starved):
        print(f"\n!! {len(starved)} ROI(s) delivered fewer than top-{BUDGET}: "
              f"{starved['n_top'].to_dict()}")
    else:
        print(f"\nall 14 ROIs delivered the full top-{BUDGET}.")


if __name__ == '__main__':
    main()
