"""F1: is largest-CC template tightening's effect resolvable once seed noise is paired out?

`Research Logs/2026-09-03-tm-axis-sweep-audit.md` Finding 5: every `read_95_delta` in
`tm_threshold_axis_sweep_largest_cc.ipynb`'s deliverable (d) is smaller than that ROI's own
seed-to-seed SD, and seed 0 sits on a distribution *endpoint* for 3 of 7 ROIs. So the
per-domain signs reported there are not measurements. This is the follow-up that log named F1.

DESIGN -- a *paired* ablation on `base_size`.

  For each (domain, seed_index in 0..4): draw ONE seed, then run the identical pipeline twice,
  changing only the template size:
      "base51"      -- fixed `tm.BASE_SIZE` (51 px), what v2 does
      "largest_cc"  -- `_odd(max side)` of the largest Otsu component in the 51 px click crop
  The seed annotation, the decoupled 5.0 um NMS radius, the chromatin `od` padding, the deep
  floor, `peak_min_distance` and `self_hit_radius` are identical between arms.

  NOT identical between arms, and it matters: the match padding `PAD = (template-1)//2` is 25
  for the 51 px arm but 12 or 9 for a 25/19 px template. Both arms reach `valid.all()`, so the
  *property* is held fixed, not the amount -- within ~25 px of the ROI edge the two arms score
  against `BORDER_REPLICATE` fabrication of different widths and the pairing is not clean there.

WHAT THE PRIMARY ANALYSIS MUST BE, decided before the run (measured from a pre-flight probe):

  The treatment is NOT constant within a domain. `301.tiff` draws tightened sizes 41, 51, 19,
  19, 27 across its five streams; `246.tiff` draws 41, 51, 47, 51, 37. A 19 px template is 14%
  of the 51 px area -- a qualitatively different search, and smaller than anything the largest-CC
  notebook ever ran (its seed-0 range was 25-51). Averaging 19 px and 41 px cells into one
  "mast cell tumor delta" produces a number with no mechanism behind it.

  Worse for the obvious framing: 10 of 35 cells are NULL cells (`tightened == 51`, the two arms
  are the same computation, delta exactly 0 by construction), and they are not evenly spread --

      201.tiff 1 informative pair | 246,402 3 | 301,459 4 | 094,548 5   (25 of 35 total)

  so a per-domain 5-seed sign test is attainable in only 2 of 7 domains, and is impossible for
  `201.tiff`. **The primary analysis is therefore a pooled dose-response**: paired delta against
  `tightened_size` (or area ratio `(tightened/51)^2`) over the 25 informative cells, with domain
  as a blocking factor. The heterogeneous treatment is an asset under that framing and a
  confound under any other.

  Lead with `recall_at_budget`, which `compare.py` designates primary and which is
  length-invariant while `n_detections >> K`. Do NOT lead with `read_95`: both arms are cut at
  the same per-map `z`, and `TM_CCOEFF` scales with template pixel count, so the arms are
  compared at the same z but NOT at a matched list length -- a `read_95` delta re-reports the
  candidate-volume difference under a new name, which is audit Finding 2 all over again. Report
  `read_95` second, beside the paired `n_detections` and `coverage_frac` deltas.

  Report the WITHIN-arm seed SD alongside every paired delta. The audit compared lcc's deltas
  against `tm_ccoeff_headtohead_seed_variance.csv`, whose SDs come from a different pipeline
  configuration and carry that caveat explicitly. This run measures the same-configuration SD
  for both arms and retires it. If the paired effect is real but small against that SD, the
  conclusion is "tightening is not a decision-relevant lever" -- and only this run can say so.

WHAT IS *NOT* A LIMITATION: seeds are drawn uniformly WITH REPLACEMENT, and two streams collide
(`201.tiff` draws ann 4457 at seeds 0 and 1; `094.tiff` draws ann 2494 at seeds 2 and 4). Those
are still iid draws -- `invariants.py`'s `check_distinct_seeds` says a collision "is expected,
not a defect" -- so `Var(mean) = sigma^2/5` already prices it and BOTH draws are kept. Dropping
them would condition on the realized sample and bias the estimate. The real limit is the
population: `201.tiff`'s entire border-filtered pool is 7 candidates.

TWO GATES, both executed rather than described:
  1. SEED-0 REPRODUCTION, checked inside the domain loop at `si == 0` so a divergence aborts on
     the first domain (~40 s) rather than after the full run. `base51` must reproduce
     `results/tm_ccoeff_threshold_axis_sweep_v2.csv` and `largest_cc` must reproduce
     `results/tm_ccoeff_threshold_axis_sweep_largest_cc.csv` on n_detections / full_list_recall
     / read_50 / read_95, both axes, every shared z. (A probe confirmed the seed draw and every
     tightened size match both saved CSVs with 0 retries in all 35 cells -- so `base51` is a
     genuine v2 replication at every seed. What this gate adds is the matching/evaluation path.)
  2. PER-CELL NULL CONTROL: in a null cell every column must be exactly equal between arms.
  Plus the one-match-many-z shortcut, re-verified at `si == 0` per domain.

Writes results/f1_seed_sweep.csv (row per domain x seed x arm x z x axis x budget, written
incrementally per domain) and results/f1_seed_sweep_cells.csv (row per domain x seed).
"""
from __future__ import annotations
import gc, os, time
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from skimage.measure import label, regionprops
from sklearn.neighbors import KDTree

from midog_utils import chromatin as cm, channels as ch, compare as cp, dataset as ds
from midog_utils import evaluate as ev, experiment as ex, find_and_suppress as fs
from midog_utils import seed_selection as ss, template_match as tm
from midog_utils.nms import nms_by_distance

CHANNEL, METHOD = 'hematoxylin_od', cv2.TM_CCOEFF
PEAK_MIN_DISTANCE, SELF_HIT_RADIUS = 7, 5.0
SEEDS = range(5)
# Includes 2.0-3.0: that is the only range in which the largest-CC notebook made a positive
# claim ("recall strictly better beyond the operating point"), so it is the range that most
# needs an error bar. Costs only extra evaluate_arms iterations on an already-extracted pool.
Z_LEVELS = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
CURRENT_Z, DEEP_FLOOR_Z = 1.0, -1.5
MAX_PEAKS = 2_000_000
NMS_RADIUS_UM, MATCH_RADIUS_UM = 5.0, ev.MIDOG_RADIUS_UM
OTSU_WINDOW, OD_PAD = tm.BASE_SIZE, tm.BASE_SIZE // 2
LARGEST_CC_KW = dict(min_area=50, max_area_frac=0.85, min_solidity=0.5)
AXES = ('tm_score', 'chromatin_od')
AXIS_RANK_KEY = {'tm_score': 'score', 'chromatin_od': 'od'}
GATE_EXACT = ['n_detections', 'tp_at_budget', 'lookalike_at_budget']
GATE_CLOSE = ['full_list_recall', 'read_50', 'read_95', 'recall_at_budget',
              'coverage_frac', 'nan_rate', 'largest_tie_block', 'n_lookalike_in_list']
REFERENCE = {'base51': 'results/tm_ccoeff_threshold_axis_sweep_v2.csv',
             'largest_cc': 'results/tm_ccoeff_threshold_axis_sweep_largest_cc.csv'}

CFG = fs.FSConfig(channel=CHANNEL, base_size=tm.BASE_SIZE, scales=(1.0,), n_angles=1,
                  flips=(False,), peak_min_distance=PEAK_MIN_DISTANCE,
                  self_hit_radius=SELF_HIT_RADIUS)
BORDER = CFG.patch_size // 2


def _odd_local(n: int, minimum: int = 5) -> int:
    n = int(round(n))
    if n % 2 == 0:
        n += 1
    return max(minimum, n)


def largest_cc_box(patch, min_area=50, max_area_frac=0.85, min_solidity=0.5):
    """Verbatim from tm_threshold_axis_sweep_largest_cc.ipynb cell 5."""
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


def draw_seed_with_retry(pool: pd.DataFrame, rng, check_fn):
    """Verbatim from tm_threshold_axis_sweep_largest_cc.ipynb cell 7."""
    working = pool.copy()
    retries = 0
    while len(working) > 0:
        idx = int(rng.integers(len(working)))
        row = working.iloc[idx]
        result = check_fn(row)
        if result is not None:
            return row, result, retries
        working = working.drop(working.index[idx])
        retries += 1
    raise ValueError("seed pool exhausted -- no candidate passed check_fn after all retries")


def run_arm(hem, H, W, seed_xy, base_size, gt_eval, match_radius, nms_radius, ctx,
            check_shortcut=False):
    """One search at one base_size -- line-for-line the v2 / largest-CC loop bodies."""
    patch73 = tm.read_padded_patch(hem, *seed_xy, CFG.patch_size)
    templates, _ = tm.build_augmentations(patch73, base_size, CFG.scales, CFG.n_angles, CFG.flips)

    PAD = max((t.shape[0] - 1) // 2 for t in templates)
    hem_p = cv2.copyMakeBorder(hem, PAD, PAD, PAD, PAD, borderType=cv2.BORDER_REPLICATE)
    fused_p, _, valid_p = tm.fused_response(hem_p, templates, CFG.scale_normalize, method=METHOD)
    fused, valid = fused_p[PAD:PAD + H, PAD:PAD + W], valid_p[PAD:PAD + H, PAD:PAD + W]
    assert valid.all(), f"{ctx['file_name']}: padding left the ROI incomplete (PAD={PAD})"
    del fused_p, valid_p, hem_p

    med, mad = tm.robust_stats(fused, valid)
    deep_floor = med + DEEP_FLOOR_Z * mad
    centers, scores = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE, deep_floor, MAX_PEAKS)
    keep = nms_by_distance(centers, scores, nms_radius)
    centers, scores = centers[keep], scores[keep]
    ok = np.hypot(centers[:, 0] - seed_xy[0], centers[:, 1] - seed_xy[1]) > CFG.self_hit_radius
    centers, scores = centers[ok], scores[ok]
    pool = pd.DataFrame({'cx': centers[:, 0], 'cy': centers[:, 1], 'score': scores})

    shortcut_ok = None
    if check_shortcut:  # the one-match-many-z economy this whole design rests on
        cut = med + CURRENT_Z * mad
        c2, s2 = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE, cut, MAX_PEAKS)
        c2 = c2[nms_by_distance(c2, s2, nms_radius)]
        c2 = c2[np.hypot(c2[:, 0] - seed_xy[0], c2[:, 1] - seed_xy[1]) > CFG.self_hit_radius]
        f = pool.loc[pool['score'] >= cut, ['cx', 'cy']].to_numpy()
        a = c2[np.lexsort((c2[:, 1], c2[:, 0]))] if len(c2) else c2
        b = f[np.lexsort((f[:, 1], f[:, 0]))] if len(f) else f
        shortcut_ok = (a.shape == b.shape) and bool(np.allclose(a, b))
        assert shortcut_ok, f"{ctx['file_name']}: one-match-many-z shortcut failed"
    del fused, valid

    hem_od_p = cv2.copyMakeBorder(hem, OD_PAD, OD_PAD, OD_PAD, OD_PAD, cv2.BORDER_REPLICATE)
    shifted = pool[['cx', 'cy']].assign(cx=pool['cx'] + OD_PAD, cy=pool['cy'] + OD_PAD)
    pool['od'] = cm.score_detections(shifted, hem_od_p)['od'].to_numpy()
    assert int(pool['od'].isna().sum()) == 0, f"{ctx['file_name']}: NaN od survived OD_PAD"
    del hem_od_p

    ps = pool.sort_values('score', ascending=False).reset_index(drop=True)
    det, _ = ev.bucket_detections(ps, gt_eval, match_radius)
    d_gt, _ = KDTree(gt_eval[['cx', 'cy']].to_numpy()).query(det[['cx', 'cy']].to_numpy(), k=1)
    inside = d_gt[:, 0] <= match_radius
    n_dup_fp = int((inside & (det['bucket'] == ev.NON_HUMAN_FINDINGS).to_numpy()).sum())

    arms = []
    for z in Z_LEVELS:
        cut = med + z * mad
        sub = pool[pool['score'] >= cut]
        for axis in AXES:
            arms.append(cp.Arm(axis, (lambda d=sub: d), rank_key=AXIS_RANK_KEY[axis], seeded=True,
                               z=z, z_dependent=True, floor_limited=bool(cut < deep_floor),
                               nms_radius=None, caps=(MAX_PEAKS,),
                               extra={'axis': axis, 'z_cut': round(cut, 5)}))
    full_ctx = dict(ctx, base_size=base_size, map_median=round(med, 5), mad_scale=round(mad, 5),
                    deep_floor=round(deep_floor, 5), n_pool=len(pool), pad_px=PAD,
                    n_dup_fp=n_dup_fp, shortcut_ok=shortcut_ok,
                    nms_radius_px=round(nms_radius, 3), match_radius_px=round(match_radius, 3))
    out = cp.evaluate_arms(arms, gt_eval, match_radius, roi_shape=(H, W), mpp=ctx['mpp'],
                           budgets=cp.BUDGETS, context=full_ctx)
    return out, len(pool)


def gate_seed0(fn, arm_name, got):
    """GATE 1 -- must reproduce the saved notebook CSV. Raises on the first divergence.

    Keyed on (arm, z, BUDGET), not (arm, z): `read_50`/`read_95` pin only two points on the
    TP-cumulative curve, so a ranking divergence that preserves them would pass a
    budget-collapsed gate and silently corrupt `recall_at_budget` -- the metric this run
    leads with. Also gates the four diagnostics `compare.py` calls mandatory rather than
    optional (`coverage_frac`, `nan_rate`, `largest_tie_block`, `n_lookalike_in_list`).
    """
    ref = pd.read_csv(REFERENCE[arm_name])
    ref = ref[ref['file_name'] == fn].set_index(['arm', 'z', 'budget'])
    mine = got.set_index(['arm', 'z', 'budget'])
    expected = len(Z_LEVELS) * len(AXES) * len(cp.BUDGETS)
    shared = [k for k in mine.index if k in ref.index]
    assert len(mine) == expected, f"{fn}/{arm_name}: produced {len(mine)} rows, expected {expected}"
    assert len(shared) == expected, (
        f"{fn}/{arm_name}: only {len(shared)} of {expected} cells exist in "
        f"{REFERENCE[arm_name]} -- a Z_LEVEL or budget outside the saved grid would be "
        "silently ungated")
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
                bad.append(f"{k} {c}: saved={a} here={b}")
    if bad:
        raise AssertionError(f"GATE 1 FAILED {fn}/{arm_name} vs {REFERENCE[arm_name]} "
                             f"({len(bad)} divergences):\n  " + "\n  ".join(bad[:12]))
    return len(shared) * (len(GATE_EXACT) + len(GATE_CLOSE))


def main():
    t0 = time.time()
    images, anns = ds.load_annotations()
    ds.check_invariants(anns)
    dom = ex.select_domain_images(images, anns, images_dir='images').sort_values('tumor_type')
    Path('results').mkdir(exist_ok=True)
    out_path, cells_path = 'results/f1_seed_sweep.csv', 'results/f1_seed_sweep_cells.csv'
    part_out, part_cells = out_path + '.partial', cells_path + '.partial'
    n_expected = len(dom) * len(SEEDS)
    completed = False

    def _dump(final: bool):
        rows = pd.concat(all_rows, ignore_index=True)
        cf = pd.DataFrame(cells).assign(run_complete=final, n_cells_expected=n_expected)
        rows.to_csv(part_out, index=False)
        cf.to_csv(part_cells, index=False)
        if final:                       # only a completed run earns the un-suffixed name
            os.replace(part_out, out_path)
            os.replace(part_cells, cells_path)
        return rows

    all_rows, cells, n_gated = [], [], 0
    try:
        for _, r in dom.iterrows():
            fn, iid, domain = r['file_name'], int(r['image_id']), r['tumor_type']
            rgb = ds.load_roi(f'images/{fn}')
            roi_shape, mpp = rgb.shape, ds.roi_mpp(f'images/{fn}')
            match_radius = ev.radius_px(mpp, MATCH_RADIUS_UM)
            nms_radius = ev.radius_px(mpp, NMS_RADIUS_UM)
            hem, gray_inv = ch.to_channel(rgb, CHANNEL), ch.to_gray_inverted(rgb)
            H, W = hem.shape[:2]
            del rgb
            gc.collect()
            gt = ds.image_annotations(anns, fn)
            pool_seed, _ = ss.agreement_pool(gt[gt['category_id'] == ds.MITOTIC])
            pool_seed = ss.border_filter(pool_seed, BORDER, roi_shape)

            for si in SEEDS:
                ts = time.time()
                rng = np.random.default_rng([si, iid])

                def _check(c, _g=gray_inv):
                    p = tm.read_padded_patch(_g, float(c['cx']), float(c['cy']), OTSU_WINDOW)
                    return None if p is None else largest_cc_box(p, **LARGEST_CC_KW)

                seed, box, n_retries = draw_seed_with_retry(pool_seed, rng, _check)
                y0, y1, x0, x1 = box
                tightened = _odd_local(max(y1 - y0, x1 - x0))
                seed_xy = (float(seed['cx']), float(seed['cy']))
                seed_ann = int(seed['ann_id'])
                gt_eval = gt[gt['ann_id'] != seed_ann].reset_index(drop=True)
                is_null = bool(tightened == tm.BASE_SIZE)

                ctx = dict(file_name=fn, image_id=iid, tumor_type=domain, mpp=round(mpp, 4),
                           seed_index=si, seed_ann_id=seed_ann, n_retries=n_retries,
                           tightened_size=tightened, null_cell=is_null,
                           area_ratio=round((tightened / tm.BASE_SIZE) ** 2, 4))
                got, npools = {}, {}
                for arm_name, bs in (('base51', tm.BASE_SIZE), ('largest_cc', tightened)):
                    out, npool = run_arm(hem, H, W, seed_xy, bs, gt_eval, match_radius,
                                         nms_radius, dict(ctx, arm_name=arm_name),
                                         check_shortcut=(si == 0))
                    got[arm_name], npools[arm_name] = out, npool
                    if si == 0:            # GATE 1 runs BEFORE the row is kept, so a failure
                        n_gated += gate_seed0(fn, arm_name, out)   # leaves no trace in the CSV
                    all_rows.append(out)

                if is_null:                                        # GATE 2
                    a = got['base51'].drop(columns=['arm_name', 'base_size']).reset_index(drop=True)
                    b = got['largest_cc'].drop(columns=['arm_name', 'base_size']).reset_index(drop=True)
                    assert a.equals(b), f"{fn} seed {si}: null cell but arms differ"

                def _n_unreach(d):     # right-censored read_95: the arm cannot reach 95%
                    return int(d.drop_duplicates(['arm', 'z'])['read_95'].isna().sum())
                cells.append(dict(**ctx, n_pool_base51=npools['base51'],
                                  n_pool_largest_cc=npools['largest_cc'],
                                  read95_unreachable_base51=_n_unreach(got['base51']),
                                  read95_unreachable_largest_cc=_n_unreach(got['largest_cc']),
                                  n_gt_mitotic=int((gt_eval['category_id'] == ds.MITOTIC).sum()),
                                  t_s=round(time.time() - ts, 1)))
                print(f"[{fn} s{si}] {domain:32s} tight={tightened:2d} "
                      f"{'NULL' if is_null else '    '} pool {npools['base51']:6d}->"
                      f"{npools['largest_cc']:6d} [{time.time() - ts:.0f}s]", flush=True)
                gc.collect()

            _dump(final=False)                                        # incremental, .partial
            del hem, gray_inv
            gc.collect()
        completed = True
    finally:
        if all_rows:
            _dump(final=completed)
            if completed:
                print(f"\nRUN COMPLETE -- wrote {out_path} and {cells_path} "
                      f"({len(cells)}/{n_expected} cells)")
            else:
                print(f"\n*** RUN INCOMPLETE *** {len(cells)}/{n_expected} cells. Output left "
                      f"at {part_out} / {part_cells}; the un-suffixed names were NOT written, "
                      "so nothing downstream can mistake this for a finished run.")

    print(f"GATE 1: {n_gated} (arm, z, metric) values reproduced the saved CSVs exactly")
    print(f"total {time.time() - t0:.0f}s")
    c = pd.DataFrame(cells)
    print(f"null cells {c.null_cell.sum()}/{len(c)}, informative {(~c.null_cell).sum()}/{len(c)}")
    print(c.groupby('tumor_type').agg(informative=('null_cell', lambda s: int((~s).sum())),
                                      distinct_seeds=('seed_ann_id', 'nunique'),
                                      tightened=('tightened_size', list)).to_string())


if __name__ == '__main__':
    main()
