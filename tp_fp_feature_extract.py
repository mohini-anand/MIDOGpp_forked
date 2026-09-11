"""Per-candidate features for the TP-vs-FP separability question.

Runs the committed search configuration -- `hematoxylin_od` channel, `cv2.TM_CCOEFF`,
largest-CC template tightening, distance NMS -- once per domain ROI, then measures a
superset of candidate features on every surviving detection and labels each one against
the ground truth. `tp_fp_separability.ipynb` consumes the cache this writes; nothing here
decides which features the analysis keeps.

The search is line-for-line `f1_seed_sweep.run_arm`'s `largest_cc` arm at `seed_index=0`,
including the `[0, image_id]` seed stream, so the candidate pool is the same pool that
run's gate checked against `results/tm_ccoeff_threshold_axis_sweep_largest_cc.csv`.

Two deliberate departures, both about *measurement* rather than search:

* **Peaks are extracted at the deep floor** (`z = -1.5`), not at the operating point.
  `f1_seed_sweep` verifies that thresholding the deep pool at `med + z*mad` reproduces a
  fresh extraction at that z exactly (its "one-match-many-z shortcut"), so one extraction
  serves every z. `map_median` and `mad_scale` are written out; the notebook subsets.
* **Shape features are gate-free.** `seed_selection.tighten_box_otsu` and
  `f1_seed_sweep.largest_cc_box` both reject a component on `min_area` /
  `max_area_frac` / `min_solidity` and return None. Those gates are correct for choosing
  a *template* and wrong for *measuring a population*: they truncate the very features
  under comparison, and they drop 9-12% of detections non-randomly (see
  `chromatin.py`'s module docstring, and `morphometric_diagnostic_audit.ipynb` sec 11).
  Here the largest component under the click is measured whatever its shape, and
  `shape_ok` records only the cases where Otsu found no component at all.

Selection and outputs
---------------------
``--images-dir`` / ``--select`` choose the ROI set. ``--select per-domain`` (the default, over
``images/``) is the original rule -- the densest ROI per tumour domain, 7 of them -- and
regenerates the committed artefacts. ``--select all`` takes every ``.tiff`` in the directory,
which over ``images/extra_valid`` is the balanced 2-per-domain draw of 14 (that folder's
MANIFEST.md).

**Output paths are derived from the selection, not fixed.** The default run writes
``.cache_tail/tp_fp_candidate_features.csv`` (a pure cache -- delete to regenerate) and
``results/tp_fp_extract_summary.csv``; any other selection gets a ``_<n>roi`` suffix, so a
14-ROI run cannot silently overwrite the 7-ROI inputs that `tp_fp_separability.ipynb`,
`chromatin_threshold_selection.ipynb` and `verify_chromatin_ranker.py` all read from the fixed
paths without asserting a row count. ``TPFP_CACHE`` / ``TPFP_SUMMARY`` override both.
"""
from __future__ import annotations

import gc
import os
import pathlib
import time

import cv2
import numpy as np
import pandas as pd
from skimage.measure import label, regionprops

from midog_utils import channels as ch
from midog_utils import chromatin as cm
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from midog_utils import experiment as ex
from midog_utils import seed_selection as ss
from midog_utils import template_match as tm
from midog_utils.nms import nms_by_distance

from f1_seed_sweep import CFG, BORDER, largest_cc_box, draw_seed_with_retry, _odd_local

CHANNEL, METHOD = 'hematoxylin_od', cv2.TM_CCOEFF
PEAK_MIN_DISTANCE, SELF_HIT_RADIUS = 7, 5.0
DEEP_FLOOR_Z, MAX_PEAKS = -1.5, 2_000_000
NMS_RADIUS_UM = 5.0
LARGEST_CC_KW = dict(min_area=50, max_area_frac=0.85, min_solidity=0.5)
OD_WINDOWS = (31, 51, 81)
CTX_WINDOW, CTX_FRAC = 121, 0.50
SHAPE_WINDOW = tm.BASE_SIZE
SEED_INDEX = 0
SHAPE_COLS = ('area', 'solidity', 'extent', 'eccentricity', 'perimeter',
              'major_axis_length', 'minor_axis_length', 'circularity',
              'tightened_size', 'mask_od_mean', 'mask_od_max', 'area_frac_of_window')

DEFAULT_CACHE = '.cache_tail/tp_fp_candidate_features.csv'
DEFAULT_SUMMARY = 'results/tp_fp_extract_summary.csv'


def output_paths(n_rois: int, is_default_selection: bool):
    """Fixed paths for the committed 7-ROI run, suffixed paths for anything else.

    Four consumers read `DEFAULT_CACHE` / `DEFAULT_SUMMARY` by hard-coded string and none of
    them asserts an ROI count, so a differently-scoped run writing there would change their
    numbers without raising anything. The suffix is the cheapest thing that makes that
    impossible. `TPFP_CACHE` / `TPFP_SUMMARY` still win, for a caller that means it.
    """
    if is_default_selection:
        cache, summary = DEFAULT_CACHE, DEFAULT_SUMMARY
    else:
        cache = DEFAULT_CACHE.replace('.csv', f'_{n_rois}roi.csv')
        summary = DEFAULT_SUMMARY.replace('.csv', f'_{n_rois}roi.csv')
    return os.environ.get('TPFP_CACHE', cache), os.environ.get('TPFP_SUMMARY', summary)


def roi_files(images_dir: str):
    """Every ROI in the directory, sorted -- `f5_nms_radius_ablation.roi_files`'s convention.

    Used with `--select all` for `images/extra_valid`, which is already the balanced
    2-per-domain draw (see that folder's MANIFEST.md). `--select per-domain` over `images/`
    -- and only that pair -- is the original one-densest-ROI-per-domain rule, so it is the
    combination that regenerates the committed 7-ROI artefacts.
    """
    return sorted(f for f in os.listdir(images_dir) if f.endswith('.tiff'))


def shape_features(chan: np.ndarray, cx: float, cy: float, window: int = SHAPE_WINDOW):
    """Largest Otsu component under a point, measured with no accept/reject gate.

    Same thresholding as `f1_seed_sweep.largest_cc_box` (per-patch min-max to uint8, then
    binary Otsu, then the largest connectivity-2 component) with every gate removed, so
    the returned distributions are not truncated at the gate boundaries. Returns None
    only when the window is unreadable or Otsu produces no foreground at all.
    """
    patch = tm.read_padded_patch(chan, cx, cy, window)
    if patch is None:
        return None
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    _, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    regions = regionprops(label(binary, connectivity=2), intensity_image=patch)
    if not regions:
        return None
    r = max(regions, key=lambda x: x.area)
    y0, x0, y1, x1 = r.bbox
    minor = r.axis_minor_length
    return {
        'area': float(r.area),
        'solidity': float(r.solidity),
        'extent': float(r.extent),
        'eccentricity': float(r.eccentricity),
        'perimeter': float(r.perimeter),
        'major_axis_length': float(r.axis_major_length),
        'minor_axis_length': float(minor),
        # `circularity` and `tightened_size` are scale-free / native-size restatements of
        # the same segmentation; kept so the notebook can choose, dropped there if redundant.
        'circularity': float(4 * np.pi * r.area / max(r.perimeter, 1e-6) ** 2),
        'tightened_size': float(max(y1 - y0, x1 - x0)),
        'mask_od_mean': float(r.intensity_mean),
        'mask_od_max': float(r.intensity_max),
        # how much of the window the object claims -- the censoring diagnostic
        'area_frac_of_window': float(r.area) / float(patch.size),
    }


def run_roi(fn: str, image_id: int, domain: str, anns: pd.DataFrame, images_dir: str):
    t0 = time.time()
    path = f'{images_dir}/{fn}'
    rgb = ds.load_roi(path)
    roi_shape, mpp = rgb.shape, ds.roi_mpp(path)
    match_radius = ev.radius_px(mpp)
    nms_radius = ev.radius_px(mpp, NMS_RADIUS_UM)
    hem = ch.to_channel(rgb, CHANNEL)
    gray_inv = ch.to_gray_inverted(rgb)
    H, W = hem.shape[:2]
    del rgb
    gc.collect()

    gt = ds.image_annotations(anns, fn)
    pool_seed, _ = ss.agreement_pool(gt[gt['category_id'] == ds.MITOTIC])
    pool_seed = ss.border_filter(pool_seed, BORDER, roi_shape)
    rng = np.random.default_rng([SEED_INDEX, image_id])

    def _check(c):
        p = tm.read_padded_patch(gray_inv, float(c['cx']), float(c['cy']), tm.BASE_SIZE)
        return None if p is None else largest_cc_box(p, **LARGEST_CC_KW)

    seed, box, n_retries = draw_seed_with_retry(pool_seed, rng, _check)
    y0, y1, x0, x1 = box
    tightened = _odd_local(max(y1 - y0, x1 - x0))
    seed_xy = (float(seed['cx']), float(seed['cy']))
    seed_ann = int(seed['ann_id'])
    gt_eval = gt[gt['ann_id'] != seed_ann].reset_index(drop=True)

    # --- search: f1_seed_sweep.run_arm's largest_cc arm ------------------------------
    patch = tm.read_padded_patch(hem, *seed_xy, CFG.patch_size)
    templates, _ = tm.build_augmentations(patch, tightened, CFG.scales, CFG.n_angles, CFG.flips)
    PAD = max((t.shape[0] - 1) // 2 for t in templates)
    hem_p = cv2.copyMakeBorder(hem, PAD, PAD, PAD, PAD, borderType=cv2.BORDER_REPLICATE)
    fused_p, _, valid_p = tm.fused_response(hem_p, templates, CFG.scale_normalize, method=METHOD)
    fused, valid = fused_p[PAD:PAD + H, PAD:PAD + W], valid_p[PAD:PAD + H, PAD:PAD + W]
    assert valid.all(), f'{fn}: padding left the ROI incomplete (PAD={PAD})'
    del fused_p, valid_p, hem_p
    gc.collect()

    med, mad = tm.robust_stats(fused, valid)
    deep_floor = med + DEEP_FLOOR_Z * mad
    centers, scores = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE, deep_floor, MAX_PEAKS)
    keep = nms_by_distance(centers, scores, nms_radius)
    centers, scores = centers[keep], scores[keep]
    ok = np.hypot(centers[:, 0] - seed_xy[0], centers[:, 1] - seed_xy[1]) > SELF_HIT_RADIUS
    centers, scores = centers[ok], scores[ok]
    del fused, valid
    gc.collect()
    t_search = time.time() - t0

    pool = pd.DataFrame({'cx': centers[:, 0], 'cy': centers[:, 1], 'score': scores})
    pool = pool.sort_values('score', ascending=False).reset_index(drop=True)

    # --- chromatin density at three windows plus a local-context contrast -------------
    t1 = time.time()
    pad = CTX_WINDOW // 2
    hem_pad = cv2.copyMakeBorder(hem, pad, pad, pad, pad, cv2.BORDER_REPLICATE)
    px, py = pool['cx'].to_numpy() + pad, pool['cy'].to_numpy() + pad
    for w in OD_WINDOWS:
        pool[f'od{w}'] = [cm.chromatin_density(hem_pad, x, y, window=w)
                          for x, y in zip(px, py)]
    pool['od_ctx'] = [cm.chromatin_density(hem_pad, x, y, window=CTX_WINDOW, frac=CTX_FRAC)
                      for x, y in zip(px, py)]
    del hem_pad
    gc.collect()
    # "darker than its own neighbourhood" and "concentrated rather than diffuse" -- both
    # within-image contrasts, so they survive the fact that OD is uncalibrated between ROIs.
    pool['od_contrast'] = pool['od51'] - pool['od_ctx']
    pool['od_falloff'] = pool['od31'] - pool['od81']
    n_od_nan = int(pool[[f'od{w}' for w in OD_WINDOWS] + ['od_ctx']].isna().sum().sum())
    assert n_od_nan == 0, f'{fn}: {n_od_nan} NaN od survived the {pad} px pad'
    t_od = time.time() - t1

    # --- gate-free shape -------------------------------------------------------------
    t2 = time.time()
    spad = SHAPE_WINDOW // 2
    hem_spad = cv2.copyMakeBorder(hem, spad, spad, spad, spad, cv2.BORDER_REPLICATE)
    rows = [shape_features(hem_spad, x + spad, y + spad)
            for x, y in zip(pool['cx'].to_numpy(), pool['cy'].to_numpy())]
    del hem_spad, hem
    gc.collect()
    shape_ok = np.array([r is not None for r in rows])
    for c in SHAPE_COLS:
        pool[c] = [r[c] if r is not None else np.nan for r in rows]
    pool['shape_ok'] = shape_ok
    t_shape = time.time() - t2

    # --- labels ----------------------------------------------------------------------
    det, _ = ev.bucket_detections(pool, gt_eval, match_radius)
    det['file_name'] = fn
    det['tumor_type'] = domain
    det['image_id'] = image_id
    det['map_median'] = med
    det['mad_scale'] = mad
    det['deep_floor'] = deep_floor
    det['z'] = (det['score'] - med) / mad
    det['tightened_size_template'] = tightened
    det['seed_ann_id'] = seed_ann
    det['mpp'] = mpp

    n_gt = int((gt_eval['category_id'] == ds.MITOTIC).sum())
    vc = det['bucket'].value_counts()
    summary = dict(
        file_name=fn, tumor_type=domain, image_id=image_id, roi_h=H, roi_w=W,
        mpp=round(mpp, 4), match_radius_px=round(match_radius, 2),
        nms_radius_px=round(nms_radius, 2), seed_ann_id=seed_ann, n_retries=n_retries,
        template_size=tightened, pad_px=PAD, map_median=round(float(med), 6),
        mad_scale=round(float(mad), 6), n_candidates=len(det), n_gt_mitotic=n_gt,
        n_tp=int(vc.get(ev.HUMAN_CORRECT_LABEL, 0)),
        n_lookalike=int(vc.get(ev.HUMAN_REJECTED_LABEL, 0)),
        n_unannotated=int(vc.get(ev.NON_HUMAN_FINDINGS, 0)),
        full_list_recall=round(int(vc.get(ev.HUMAN_CORRECT_LABEL, 0)) / max(n_gt, 1), 4),
        shape_fail_rate=round(float(1 - shape_ok.mean()), 5),
        t_search_s=round(t_search, 1), t_od_s=round(t_od, 1), t_shape_s=round(t_shape, 1),
    )
    print(f"[{fn}] {domain:32s} tmpl={tightened:2d} cand={len(det):6d} "
          f"TP={summary['n_tp']:4d}/{n_gt} look={summary['n_lookalike']:4d} "
          f"shape_fail={summary['shape_fail_rate']:.3%} "
          f"[{time.time() - t0:.0f}s]", flush=True)
    return det, summary


def main(images_dir: str = 'images', select: str = 'per-domain'):
    images, anns = ds.load_annotations()
    ds.check_invariants(anns)
    if select == 'per-domain':
        dom = ex.select_domain_images(images, anns, images_dir=images_dir)
    else:
        on_disk = roi_files(images_dir)
        dom = images[images['file_name'].isin(on_disk)]
        missing = sorted(set(on_disk) - set(dom['file_name']))
        if missing:      # an intersection drops these silently otherwise
            print(f'WARNING: {len(missing)} .tiff on disk absent from the annotation DB, '
                  f'skipped: {missing}', flush=True)
        # `select_domain_images` filters n_mitotic > 0; the `all` branch has to do it too, or
        # a zero-mitosis ROI (001.tiff) dies in `draw_seed_with_retry` partway through the run.
        n_mit = (anns[anns['category_id'] == ds.MITOTIC].groupby('file_name').size())
        seedable = dom['file_name'].map(n_mit).fillna(0) > 0
        if not seedable.all():
            print(f'WARNING: dropping {int((~seedable).sum())} ROI(s) with no mitotic '
                  f'annotation: {sorted(dom.loc[~seedable, "file_name"])}', flush=True)
        dom = dom[seedable]
    dom = dom.sort_values(['tumor_type', 'file_name']).reset_index(drop=True)
    if len(dom) == 0:
        raise SystemExit(f'no seedable ROI found in {images_dir}/ with --select {select}')
    cache, summary_path = output_paths(len(dom), images_dir == 'images' and select == 'per-domain')
    frames, summaries = [], []
    t0 = time.time()
    print(f'{len(dom)} ROIs from {images_dir}/ ({select})\n  cache   -> {cache}\n'
          f'  summary -> {summary_path}', flush=True)
    for _, r in dom.iterrows():
        det, s = run_roi(r['file_name'], int(r['image_id']), r['tumor_type'], anns,
                         images_dir=images_dir)
        frames.append(det)
        summaries.append(s)
        gc.collect()
    out = pd.concat(frames, ignore_index=True)
    d = os.path.dirname(cache)
    if d:
        os.makedirs(d, exist_ok=True)
        pathlib.Path(d, '.gitignore').write_text('*\n')
    out.to_csv(cache, index=False)
    ds_ = os.path.dirname(summary_path)
    if ds_:                                  # symmetric with the cache: TPFP_SUMMARY may point
        os.makedirs(ds_, exist_ok=True)      # anywhere, and failing here loses the whole run
    summary_df = pd.DataFrame(summaries)
    summary_df.insert(0, 'images_dir', images_dir)   # provenance: which set produced this file
    summary_df.insert(1, 'select', select)
    summary_df.to_csv(summary_path, index=False)
    print(f'\n{len(out)} candidates across {len(summaries)} ROIs in {time.time() - t0:.0f}s')
    print(f'  cache   -> {cache}')
    print(f'  summary -> {summary_path}')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--images-dir', default='images')
    ap.add_argument('--select', default='per-domain', choices=('per-domain', 'all'))
    a = ap.parse_args()
    main(images_dir=a.images_dir, select=a.select)
