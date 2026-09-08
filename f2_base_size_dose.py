"""F2: where is the template-size threshold, when the dose is *assigned* rather than observed?

Pre-registration: `Research Logs/2026-09-08-f2-preregistration.md` (committed `b571811`,
before this file existed). Every constant below is fixed by that document; changing one
invalidates the pre-registration.

DESIGN -- one factor, seven levels, every level in every cell.

  For each of 14 ROIs (`images/extra_valid/`, the set `DECISIONS.md` D5 requires) and each of
  5 seeds, draw ONE click and run the identical pipeline seven times, changing only the
  template size:

      b51  b43  b37  b33  b27  b19        assigned doses; b51 is the reference level
      lcc                                 observational: `_odd(max side)` of the largest Otsu
                                          component in the 51 px crop -- the shipped rule

  F1 could not do this: its dose was whatever `largest_cc_box` happened to pick, so the
  contrast was between cells and the analysis had to correlate a delta against an observed
  covariate. Here the contrast is *within* the cell, and the analysis is a mean per dose.

  Held fixed across all seven arms: the seed click, the evaluation GT, the channel, METHOD,
  peak_min_distance, self_hit_radius, the decoupled 5.0 um NMS radius, OD_PAD, the deep
  floor, Z_LEVELS and BUDGETS. NOT fixed, and it cannot be: `PAD = (base_size-1)//2` runs 9
  at 19 px against 25 at 51 px, so within ~25 px of the ROI edge the arms score against
  BORDER_REPLICATE fabrication of different widths (pre-registration section 3a).

TWO THRESHOLD RULES, both computed here rather than offline, because neither costs a second
template match:

  matched_z  -- `score >= med + z*mad`, six z levels. F1's rule; the shipped behaviour.
  matched_n  -- top-N by score, N = THIS CELL's b51 count at z = 1.0. Emitted once per
                (cell, dose) with z = NaN and `z_dependent=False`, which is what `compare.Arm`
                documents for arms that do not depend on z.

  IMPLEMENTATION TRAP, named in the pre-registration so it cannot happen: the dose loop must
  run b51 FIRST, because N is defined from it. `DOSES[0] == 51` is asserted at import.

  matched_n exists to decide audit finding F1.5. On `tm_score` it is a formality -- recall@K
  is exactly invariant to the threshold wherever the budget is delivered, measured 0 of 560
  in F1 and re-asserted by GATE 3 here -- so the rule bites only on `chromatin_od`, where the
  membership gate (`score`) and the ranking key (`od`) are different quantities.

SEEDS come from `tm_variant_sweep.draw_seeds`: without replacement, one RNG stream per seed
index. F1 drew with replacement and produced two byte-identical cells that its sign-flip null
then treated as independent strata. Verified before the pre-registration was written that
`draw_seeds`' seed 0 returns F1's exact `ann_id` on all 7 shared ROIs, which is what makes
GATE 1 possible.

Writes results/f2_base_size_dose.csv, _cells.csv and _verification.csv (incrementally, per
ROI, under a `.partial` suffix until the run completes), plus the raw pools to
.cache_tail/f2_pools/ so any further threshold rule can be evaluated without re-running.
"""
from __future__ import annotations
import argparse, gc, os, pathlib, time
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from skimage.measure import label, regionprops

from midog_utils import chromatin as cm, channels as ch, compare as cp, dataset as ds
from midog_utils import evaluate as ev, find_and_suppress as fs, invariants as inv
from midog_utils import template_match as tm
from midog_utils.nms import nms_by_distance
from tm_variant_sweep import draw_seeds

# --- pre-registered constants (sections 3a, 3d, 6.1, 8) ---------------------------------
IMAGES_DIR = "images/extra_valid"          # section 3c; the directory IS the 14-ROI selection
CHANNEL, METHOD = 'hematoxylin_od', cv2.TM_CCOEFF          # D3, and D1 at the call site
PEAK_MIN_DISTANCE, SELF_HIT_RADIUS = 7, 5.0
N_SEEDS = 5
DOSES = (51, 43, 37, 33, 27, 19)           # 51 FIRST -- matched_n's N is defined from it
REFERENCE_DOSE = 51
Z_LEVELS = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
DEEP_FLOOR_Z = -1.5
MAX_PEAKS = 2_000_000
NMS_RADIUS_UM, MATCH_RADIUS_UM = 5.0, ev.MIDOG_RADIUS_UM
OTSU_WINDOW, OD_PAD = tm.BASE_SIZE, tm.BASE_SIZE // 2
LARGEST_CC_KW = dict(min_area=50, max_area_frac=0.85, min_solidity=0.5)
AXES = ('tm_score', 'chromatin_od')
AXIS_RANK_KEY = {'tm_score': 'score', 'chromatin_od': 'od'}

# The protocol constants GATE 10 asserts against `f2_analysis.py` and the pre-registration.
PRIMARY_AXIS, Z_PRIMARY, K_PRIMARY = 'tm_score', 1.0, 250
HOLM_FAMILY = 5                            # the five dose tests; NOT ten (section 8)

POOL_DIR = '.cache_tail/f2_pools'
OUT_MAIN = 'results/f2_base_size_dose.csv'
OUT_CELLS = 'results/f2_base_size_dose_cells.csv'
OUT_VERIF = 'results/f2_base_size_dose_verification.csv'

# GATE 1: F1's saved run, and the arm-name map between the two schemas.
F1_CSV = 'results/f1_seed_sweep.csv'
F1_ARM_NAME = {'b51': 'base51', 'lcc': 'largest_cc'}
GATE_EXACT = ['n_detections', 'tp_at_budget', 'lookalike_at_budget']
GATE_CLOSE = ['full_list_recall', 'read_50', 'read_95', 'recall_at_budget',
              'coverage_frac', 'nan_rate', 'largest_tie_block', 'n_lookalike_in_list']

CFG = fs.FSConfig(channel=CHANNEL, base_size=tm.BASE_SIZE, scales=(1.0,), n_angles=1,
                  flips=(False,), peak_min_distance=PEAK_MIN_DISTANCE,
                  self_hit_radius=SELF_HIT_RADIUS)
BORDER = CFG.patch_size // 2               # 36 -- the seed's own rotation-safe margin

assert DOSES[0] == REFERENCE_DOSE, "b51 must run first: matched_n's N is defined from it"
assert Z_PRIMARY in Z_LEVELS and K_PRIMARY in cp.BUDGETS
assert HOLM_FAMILY == len(DOSES) - 1, "the Holm family is the non-reference doses"


def roi_files(images_dir: str = IMAGES_DIR):
    """The 14-ROI selection is the directory, not a predicate -- `MANIFEST.md` says so.

    Same two lines as `f5_nms_radius_ablation.roi_files`. There is no filter to call: the
    "2 per tumour type, n_mitotic >= 15, seed pool >= 5" guarantee was applied by hand when
    the folder was assembled, and lives in how it was built rather than in anything that
    re-checks it on read.
    """
    return sorted(f for f in os.listdir(images_dir) if f.endswith('.tiff'))


def _odd_local(n: int, minimum: int = 5) -> int:
    """Verbatim from `f1_seed_sweep.py`, itself verbatim from the largest-CC notebook."""
    n = int(round(n))
    if n % 2 == 0:
        n += 1
    return max(minimum, n)


def largest_cc_box(patch, min_area=50, max_area_frac=0.85, min_solidity=0.5):
    """Verbatim from `f1_seed_sweep.py` -- byte-identical so GATE 1 compares like with like.

    The `cv2.normalize` here is a *local* rescale inside a thresholding routine, which D3
    explicitly distinguishes from a global rescale of the search channel.
    """
    if patch.ndim != 2:
        raise ValueError(f"largest_cc_box needs a single-channel patch, got shape {patch.shape}")
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    _, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    regions = regionprops(label(binary, connectivity=2))
    if not regions:
        return None
    largest = max(regions, key=lambda r: r.area)
    if largest.area < min_area:
        return None
    if largest.area > max_area_frac * patch.size:
        return None
    if largest.solidity < min_solidity:
        return None
    y0, x0, y1, x1 = largest.bbox
    return int(y0), int(y1), int(x0), int(x1)


def build_pool(hem, H, W, seed_xy, base_size, nms_radius, ctx, checks, check_shortcut=False):
    """One search at one base_size. Line-for-line `f1_seed_sweep.run_arm`'s first half.

    Returns ``(pool, meta)``. The pool carries cx, cy, score and od and is the only thing
    downstream needs, which is why it is also what gets persisted.
    """
    patch73 = tm.read_padded_patch(hem, *seed_xy, CFG.patch_size)
    templates, _ = tm.build_augmentations(patch73, base_size, CFG.scales, CFG.n_angles, CFG.flips)

    PAD = max((t.shape[0] - 1) // 2 for t in templates)
    hem_p = cv2.copyMakeBorder(hem, PAD, PAD, PAD, PAD, borderType=cv2.BORDER_REPLICATE)
    fused_p, _, valid_p = tm.fused_response(hem_p, templates, CFG.scale_normalize, method=METHOD)
    fused, valid = fused_p[PAD:PAD + H, PAD:PAD + W], valid_p[PAD:PAD + H, PAD:PAD + W]
    assert valid.all(), f"{ctx['label']}: padding left the ROI incomplete (PAD={PAD})"   # GATE 8
    del fused_p, valid_p, hem_p

    med, mad = tm.robust_stats(fused, valid)
    deep_floor = med + DEEP_FLOOR_Z * mad
    centers, scores = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE, deep_floor, MAX_PEAKS)
    # GATE 7: the PRE-NMS length, which is the one that could actually equal the cap.
    checks.append(inv.check_no_cap(len(centers), (MAX_PEAKS,), label=f"{ctx['label']}/pre_nms"))

    keep = nms_by_distance(centers, scores, nms_radius)
    centers, scores = centers[keep], scores[keep]
    ok = np.hypot(centers[:, 0] - seed_xy[0], centers[:, 1] - seed_xy[1]) > CFG.self_hit_radius
    centers, scores = centers[ok], scores[ok]
    pool = pd.DataFrame({'cx': centers[:, 0], 'cy': centers[:, 1], 'score': scores})

    shortcut_ok = None
    if check_shortcut:                     # GATE 9 -- the one-match-many-z economy
        cut = med + Z_PRIMARY * mad
        c2, s2 = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE, cut, MAX_PEAKS)
        c2 = c2[nms_by_distance(c2, s2, nms_radius)]
        c2 = c2[np.hypot(c2[:, 0] - seed_xy[0], c2[:, 1] - seed_xy[1]) > CFG.self_hit_radius]
        f = pool.loc[pool['score'] >= cut, ['cx', 'cy']].to_numpy()
        a = c2[np.lexsort((c2[:, 1], c2[:, 0]))] if len(c2) else c2
        b = f[np.lexsort((f[:, 1], f[:, 0]))] if len(f) else f
        shortcut_ok = (a.shape == b.shape) and bool(np.allclose(a, b))
        assert shortcut_ok, f"{ctx['label']}: one-match-many-z shortcut failed"
        checks.append({'check': 'shortcut_prefix', 'label': ctx['label'], 'passed': True,
                       'n_fresh': int(len(a)), 'n_prefix': int(len(b))})
    del fused, valid

    # OD_PAD is deliberately independent of PAD: the chromatin window is fixed at 51 px, so
    # reusing a PAD that shrinks with base_size would silently NaN border candidates.
    hem_od_p = cv2.copyMakeBorder(hem, OD_PAD, OD_PAD, OD_PAD, OD_PAD, cv2.BORDER_REPLICATE)
    shifted = pool[['cx', 'cy']].assign(cx=pool['cx'] + OD_PAD, cy=pool['cy'] + OD_PAD)
    pool['od'] = cm.score_detections(shifted, hem_od_p)['od'].to_numpy()
    assert int(pool['od'].isna().sum()) == 0, f"{ctx['label']}: NaN od survived OD_PAD"  # GATE 8
    del hem_od_p

    # FULL PRECISION, deliberately. `arms_for_dose` computes `med + z*mad` from these, and
    # rounding first shifts every cut: at 094.tiff seed 0, med/mad rounded to 5 dp moves the
    # z=2.0 cut from 0.226041718850 to 0.22605 and drops one candidate. F1 rounded only for
    # its context dict and thresholded on the raw floats; GATE 1 caught the difference.
    meta = dict(med=med, mad=mad, deep_floor=deep_floor, pad_px=PAD, n_pool=len(pool),
                n_peaks_deep=int(len(keep)), shortcut_ok=shortcut_ok)
    return pool, meta


def meta_columns(meta):
    """The rounded, human-readable form of `meta` for the results CSV. Display only."""
    return dict(map_median=round(meta['med'], 5), mad_scale=round(meta['mad'], 5),
                deep_floor=round(meta['deep_floor'], 5), pad_px=meta['pad_px'],
                n_pool=meta['n_pool'], n_peaks_deep=meta['n_peaks_deep'],
                shortcut_ok=meta['shortcut_ok'])


def arms_for_dose(pool, meta, dose_tag, base_size, n_ref):
    """The 12 matched_z arms and the 2 matched_n arms for one (cell, dose).

    `matched_n` takes the top `n_ref` by score. It is emitted with ``z=nan`` and
    ``z_dependent=False`` -- `compare.Arm`'s own convention for an arm that does not depend
    on z, which the harness then emits once instead of duplicating at every z.
    """
    med, mad = meta['med'], meta['mad']            # raw, never the rounded CSV values
    ranked_by_score = pool.sort_values('score', ascending=False, kind='mergesort')
    matched_n = ranked_by_score.head(n_ref)
    arms = []
    for z in Z_LEVELS:
        cut = med + z * mad
        sub = pool[pool['score'] >= cut]
        for axis in AXES:
            arms.append(cp.Arm(f'{axis}@{dose_tag}', (lambda d=sub: d),
                               rank_key=AXIS_RANK_KEY[axis], seeded=True, z=z,
                               z_dependent=True, floor_limited=bool(cut < meta['deep_floor']),
                               nms_radius=None, caps=(MAX_PEAKS,),
                               extra={'axis': axis, 'dose_tag': dose_tag,
                                      'base_size': base_size, 'rule': 'matched_z',
                                      'n_target': -1, 'z_cut': round(cut, 5)}))
    for axis in AXES:
        arms.append(cp.Arm(f'{axis}@{dose_tag}', (lambda d=matched_n: d),
                           rank_key=AXIS_RANK_KEY[axis], seeded=True, z=float('nan'),
                           z_dependent=False, floor_limited=False,
                           nms_radius=None, caps=(MAX_PEAKS,),
                           extra={'axis': axis, 'dose_tag': dose_tag,
                                  'base_size': base_size, 'rule': 'matched_n',
                                  'n_target': int(n_ref), 'z_cut': float('nan')}))
    return arms


def gate_seed0_vs_f1(fn, dose_tag, got, checks):
    """GATE 1 -- `b51` and `lcc` at seed 0 must reproduce F1's saved run on the shared ROIs.

    Keyed on (axis, z, budget), not (axis, z): `read_50`/`read_95` pin only two points on the
    TP-cumulative curve, so a ranking divergence that preserved them would pass a
    budget-collapsed gate and silently corrupt `recall_at_budget`, the metric this run leads
    with. Runs inside the ROI loop, so a divergence aborts on the first ROI.

    Only `matched_z` rows are compared: F1 has no `matched_n` arm to compare against.
    """
    ref = pd.read_csv(F1_CSV)
    # `seed_index == 0` is not optional: F1 wrote all five seeds to one CSV, so without it
    # (arm, z, budget) selects five rows and `.loc` returns a Series instead of a scalar.
    ref = ref[(ref['file_name'] == fn) & (ref['arm_name'] == F1_ARM_NAME[dose_tag])
              & (ref['seed_index'] == 0)]
    if not len(ref):
        return None                        # one of the 7 ROIs F1 never ran -- nothing to gate
    ref = ref.set_index(['arm', 'z', 'budget']).sort_index()
    mine = got[got['rule'] == 'matched_z'].copy()
    mine = mine.set_index(['axis', 'z', 'budget']).sort_index()
    assert ref.index.is_unique and mine.index.is_unique, \
        f"{fn}/{dose_tag}: gate index is not unique -- the comparison would be ambiguous"
    expected = len(Z_LEVELS) * len(AXES) * len(cp.BUDGETS)
    assert len(mine) == expected, f"{fn}/{dose_tag}: {len(mine)} matched_z rows, want {expected}"
    shared = [k for k in mine.index if k in ref.index]
    assert len(shared) == expected, (
        f"{fn}/{dose_tag}: only {len(shared)} of {expected} cells exist in {F1_CSV} -- a "
        "Z_LEVEL or budget outside F1's grid would be silently ungated")
    bad = []
    for k in shared:
        for c in GATE_EXACT + GATE_CLOSE:
            a, b = ref.loc[k, c], mine.loc[k, c]
            if pd.isna(a) or pd.isna(b):
                same = bool(pd.isna(a) and pd.isna(b))
            elif c in GATE_EXACT:
                same = float(a) == float(b)
            else:
                same = bool(np.isclose(float(a), float(b), rtol=0, atol=1e-9))
            if not same:
                bad.append(f"{k} {c}: f1={a} here={b}")
    if bad:
        raise AssertionError(f"GATE 1 FAILED {fn}/{dose_tag} vs {F1_CSV} "
                             f"({len(bad)} divergences):\n  " + "\n  ".join(bad[:12]))
    n = len(shared) * (len(GATE_EXACT) + len(GATE_CLOSE))
    checks.append({'check': 'f1_seed0_reproduction', 'label': f'{fn}/{dose_tag}',
                   'passed': True, 'n_values': n})
    return n


def gate_prefix_property(df, checks):
    """GATE 3 -- the property section 2a's whole argument rests on, in the form that can fail.

    Wherever the budget is actually delivered, `recall_at_budget` on `tm_score` must not
    move with z. Asserting it for K <= 2000 instead would pass trivially: the 25 of 560 F1
    cells that vary at K = 5000 all have `budget_delivered` below the budget.
    """
    s = df[(df['axis'] == PRIMARY_AXIS) & (df['rule'] == 'matched_z')
           & (df['budget_delivered'] == df['budget'])]
    g = s.groupby(['file_name', 'seed_index', 'dose_tag', 'budget'])['recall_at_budget'].nunique()
    bad = g[g > 1]
    if len(bad):
        raise AssertionError(
            f"GATE 3 FAILED: recall_at_budget on {PRIMARY_AXIS} moves with z in {len(bad)} "
            f"fully-delivered (cell, dose, budget) groups -- the prefix property is broken:\n"
            f"{bad.head(10).to_string()}")
    checks.append({'check': 'prefix_property', 'label': 'all', 'passed': True,
                   'n_groups': int(len(g))})


def gate_pool_roundtrip(fn, si, dose_tag, pool_path, meta, base_size, n_ref, gt_eval,
                        match_radius, roi_shape, mpp, ctx, got, checks):
    """GATE 2 -- a persisted pool must re-evaluate to the rows the run wrote.

    Without this every future offline re-analysis of `.cache_tail/f2_pools/` is unverified,
    and the pools exist precisely so a further threshold rule need not re-run the sweep.
    """
    back = pd.read_parquet(pool_path).astype({'cx': 'float64', 'cy': 'float64',
                                              'score': 'float32', 'od': 'float64'})
    arms = arms_for_dose(back, meta, dose_tag, base_size, n_ref)
    out = cp.evaluate_arms(arms, gt_eval, match_radius, roi_shape=roi_shape, mpp=mpp,
                           budgets=cp.BUDGETS, context=dict(ctx))
    cols = ['n_detections', 'recall_at_budget', 'read_95', 'full_list_recall', 'coverage_frac']
    key = ['axis', 'rule', 'z', 'budget']
    a = got.sort_values(key).reset_index(drop=True)
    b = out.sort_values(key).reset_index(drop=True)
    assert len(a) == len(b), f"GATE 2 {fn} s{si} {dose_tag}: {len(a)} rows vs {len(b)}"
    for c in cols:
        if not np.allclose(a[c].astype(float), b[c].astype(float), rtol=0, atol=1e-9,
                           equal_nan=True):
            raise AssertionError(f"GATE 2 FAILED {fn} s{si} {dose_tag}: column {c} differs "
                                 "between the run and a re-evaluation of its persisted pool")
    checks.append({'check': 'pool_roundtrip', 'label': f'{fn}/s{si}/{dose_tag}',
                   'passed': True, 'n_values': int(len(a) * len(cols))})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--images-dir', default=IMAGES_DIR)
    ap.add_argument('--seeds', type=int, default=N_SEEDS)
    ap.add_argument('--smoke', nargs='?', const='201.tiff', default=None, metavar='ROI',
                    help='one ROI, one seed -- exercises GATES 1, 2, 7, 8, 9. Defaults to '
                         '201.tiff; name another to reproduce a gate failure on it.')
    args = ap.parse_args()

    t0 = time.time()
    images, anns = ds.load_annotations()
    ds.check_invariants(anns)
    files = roi_files(args.images_dir)
    n_seeds = args.seeds
    if args.smoke:
        files, n_seeds = [args.smoke], 1     # default is a ROI F1 ran, so GATE 1 fires
    Path('results').mkdir(exist_ok=True)
    pathlib.Path(POOL_DIR).mkdir(parents=True, exist_ok=True)
    pathlib.Path(POOL_DIR).parent.joinpath('.gitignore').write_text('*\n')

    n_expected = len(files) * n_seeds
    all_rows, cells, checks, completed = [], [], [], False
    suffix = '.smoke' if args.smoke else ''

    def _dump(final: bool):
        rows = pd.concat(all_rows, ignore_index=True)
        cf = pd.DataFrame(cells).assign(run_complete=final, n_cells_expected=n_expected)
        for path, frame in ((OUT_MAIN, rows), (OUT_CELLS, cf),
                            (OUT_VERIF, pd.DataFrame(checks))):
            frame.to_csv(path + suffix + '.partial', index=False)
            if final:                        # only a completed run earns the plain name
                os.replace(path + suffix + '.partial', path + suffix)

    try:
        for fn in files:
            path = f'{args.images_dir}/{fn}'
            rgb = ds.load_roi(path)
            roi_shape, mpp = rgb.shape, ds.roi_mpp(path)
            match_radius = ev.radius_px(mpp, MATCH_RADIUS_UM)
            nms_radius = ev.radius_px(mpp, NMS_RADIUS_UM)
            hem, gray_inv = ch.to_channel(rgb, CHANNEL), ch.to_gray_inverted(rgb)
            H, W = hem.shape[:2]
            del rgb
            gc.collect()

            gt = ds.image_annotations(anns, fn)
            gt_mit = gt[gt['category_id'] == ds.MITOTIC]
            iid, domain = int(gt_mit['image_id'].iloc[0]), images.loc[
                images['file_name'] == fn, 'tumor_type'].iloc[0]
            seeds, sinfo, seed_records = draw_seeds(gt_mit, BORDER, roi_shape, n_seeds)
            assert len(seeds) == n_seeds, f"{fn}: {len(seeds)} seeds, wanted {n_seeds}"
            checks.append(inv.check_distinct_seeds(seed_records, pool_size=sinfo.n_after_border,
                                                   label=fn))                        # GATE 6

            for si, seed in seeds:
                ts = time.time()
                seed_xy = (float(seed['cx']), float(seed['cy']))
                seed_ann = int(seed['ann_id'])
                gt_eval = gt[gt['ann_id'] != seed_ann].reset_index(drop=True)

                p = tm.read_padded_patch(gray_inv, *seed_xy, OTSU_WINDOW)
                box = largest_cc_box(p, **LARGEST_CC_KW) if p is not None else None
                assert box is not None, (
                    f"{fn} s{si}: largest_cc_box rejected the drawn seed. F1 recorded 0 "
                    "retries in all 35 cells; a rejection here means the draw diverged.")
                lcc_size = _odd_local(max(box[1] - box[0], box[3] - box[2]))

                ctx0 = dict(file_name=fn, image_id=iid, tumor_type=domain, mpp=round(mpp, 4),
                            seed_index=si, seed_ann_id=seed_ann, lcc_size=lcc_size,
                            nms_radius_px=round(nms_radius, 3),
                            match_radius_px=round(match_radius, 3))
                n_ref, timings, metas = None, {}, {}
                # `lcc` frequently lands on an assigned dose -- 20% of candidates tighten to
                # exactly 51 (pre-registration section 2b), and the grid covers six of the 17
                # odd sizes in range. Re-running an identical search would produce identical
                # rows at full cost, so the search is memoised by `base_size` within the cell.
                # This is an optimisation, not a protocol change: a cache hit IS the same
                # computation, which is why F1's null cells were byte-identical between arms.
                by_size = {}

                for dose_tag, bs in [(f'b{d}', d) for d in DOSES] + [('lcc', lcc_size)]:
                    td = time.time()
                    ctx = dict(ctx0, label=f'{fn}/s{si}/{dose_tag}')
                    reused = bs in by_size
                    if reused:
                        pool, meta = by_size[bs]
                    else:
                        pool, meta = build_pool(hem, H, W, seed_xy, bs, nms_radius, ctx,
                                                checks, check_shortcut=(si == 0))
                        by_size[bs] = (pool, meta)
                    if dose_tag == f'b{REFERENCE_DOSE}':
                        n_ref = int((pool['score'] >= meta['med']
                                     + Z_PRIMARY * meta['mad']).sum())
                    assert n_ref is not None, "b51 must run first"
                    assert len(pool) >= n_ref, (                                     # GATE 4
                        f"{ctx['label']}: deep pool {len(pool)} < N={n_ref}, so matched_n "
                        "would silently return a shorter list than the reference arm")

                    pool_path = f'{POOL_DIR}/{fn[:-5]}_s{si}_{dose_tag}.parquet'
                    pool.astype({'cx': 'int32', 'cy': 'int32',
                                 'score': 'float32', 'od': 'float64'}
                                ).to_parquet(pool_path, index=False)

                    full = dict(ctx0, dose_tag=dose_tag, base_size=bs, n_target_ref=n_ref,
                                search_reused=reused, **meta_columns(meta))
                    out = cp.evaluate_arms(arms_for_dose(pool, meta, dose_tag, bs, n_ref),
                                           gt_eval, match_radius, roi_shape=roi_shape,
                                           mpp=mpp, budgets=cp.BUDGETS, context=full,
                                           checks=checks)
                    got_n = out.loc[out['rule'] == 'matched_n', 'n_detections'].unique()
                    assert list(got_n) == [n_ref], (                                 # GATE 5
                        f"{ctx['label']}: matched_n delivered {got_n}, not N={n_ref}")

                    if si == 0 and dose_tag in F1_ARM_NAME:
                        gate_seed0_vs_f1(fn, dose_tag, out, checks)                   # GATE 1
                    if si == 0 and dose_tag == f'b{REFERENCE_DOSE}':
                        gate_pool_roundtrip(fn, si, dose_tag, pool_path, meta, bs, n_ref,
                                            gt_eval, match_radius, roi_shape, mpp, full,
                                            out, checks)                              # GATE 2
                    all_rows.append(out)
                    metas[dose_tag], timings[f't_{dose_tag}_s'] = meta, round(time.time() - td, 1)

                cells.append(dict(ctx0, n_ref=n_ref, n_gt_mitotic=int(
                    (gt_eval['category_id'] == ds.MITOTIC).sum()), seed_pool=sinfo.n_after_border,
                    **{f'n_pool_{k}': v['n_pool'] for k, v in metas.items()},
                    lcc_search_reused=bool(lcc_size in {d for d in DOSES}),
                    **timings, t_s=round(time.time() - ts, 1)))
                del by_size
                gc.collect()
                print(f"[{fn} s{si}] {domain:32s} lcc={lcc_size:2d} N={n_ref:6d} "
                      f"pools {[metas[f'b{d}']['n_pool'] for d in DOSES]} "
                      f"[{time.time() - ts:.0f}s]", flush=True)

            _dump(final=False)                                                       # GATE 12
            del hem, gray_inv
            gc.collect()
        completed = True
    finally:
        if all_rows:
            df = pd.concat(all_rows, ignore_index=True)
            if completed:
                gate_prefix_property(df, checks)                                     # GATE 3
                cp.assert_floor_not_limiting(df[df['rule'] == 'matched_z'])          # GATE 11
            _dump(final=completed)
            if completed:
                print(f"\nRUN COMPLETE -- {len(cells)}/{n_expected} cells, {len(df)} rows")
                print(f"  {OUT_MAIN}{suffix}\n  {OUT_CELLS}{suffix}\n  {OUT_VERIF}{suffix}")
            else:
                print(f"\n*** RUN INCOMPLETE *** {len(cells)}/{n_expected} cells. Output left "
                      f"at *{suffix}.partial; the plain names were NOT written, so nothing "
                      "downstream can mistake this for a finished run.")

    n_gate1 = sum(c.get('n_values', 0) for c in checks
                  if c['check'] == 'f1_seed0_reproduction')
    print(f"GATE 1: {n_gate1} values reproduced {F1_CSV} exactly")
    print(f"gates recorded: {len(checks)}; all passed = "
          f"{all(c.get('passed', True) for c in checks)}")
    print(f"total {time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
