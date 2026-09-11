"""Build the recall-workload ledger: every mitosis's rank in its own deep candidate pool.

Why a ledger of *ranks* and not another z grid
----------------------------------------------
`tm_threshold_axis_sweep_largest_cc_high_z.ipynb` (below: *high_z*) establishes the prefix
identity: the deep pool is score-descending, a `score >= cut` filter is exactly a prefix of
it, and `evaluate.greedy_match` on a prefix is the first N steps of the full pass. One
consequence that notebook states but does not exploit -- `recall_at_budget` on the `tm_score`
axis is *exactly* invariant to the cutoff for every z at or below `z_max`, measured at
280/280 (ROI, seed, budget) cells. So recall@K cannot choose a threshold on this axis: below
`z_max` it is flat by construction, and above `z_max` high_z does not look at all.

The quantity that is not flat is the pair **(candidates retained, recall achieved)** as the
cutoff rises *through and past* `z_max`. That curve is a step function which changes value
only at a true positive, so it is completely determined by, for each mitosis, its **rank in
the score-descending deep pool**. Storing those m ranks per (ROI, seed) is a lossless
representation of the entire curve at every achievable operating point -- and it makes both
readings of the curve exact and free:

* fix the recall tolerance, read off the cutoff and the list length  (this script's Table 1)
* fix the list length K, read off the recall                          (recall@K, D4's metric)

Both are the *same* step function; `compare._read_depths` already computes six points on it
(`read_50` ... `read_100`) and `evaluate.recall_at_k` the other reading. Neither is stored
over seeds anywhere in this repo, and no artefact stores the curve itself.

Scope
-----
`images/extra_valid` -- 2 ROIs per tumour type over all 7 domains, the set `DECISIONS.md` D5
names and `f5_nms_radius_ablation.py` already uses. high_z ran one ROI per domain, which left
the *ROI* half of D4's "worst ROI, worst click, per domain" rule untested while its stated
rationale is precisely that ROIs cluster by tumour type. All 7 of high_z's ROIs are inside
this 14, so `z_max` on those must reproduce bit-for-bit -- that is the reproduction gate the
notebook checks first.

Ranking is on `tm_score` only. `chromatin_od` is not computed: D5 keeps it off the arbiter,
high_z showed its one measurable advantage evaporates at every deployable cutoff, and
skipping the per-detection OD window is most of this script's runtime.

Everything else is copied from high_z rather than re-derived -- channel, method,
largest-CC tightening, the seed-draw rule, the border/NMS fixes, `DEEP_FLOOR_Z`. The seed
draw in particular is high_z's `draw_seed_with_retry` (with replacement across seed indices),
*not* `tm_variant_sweep.draw_seeds` (without replacement), because keeping it is what makes
the 7 shared ROIs reproduce exactly. The cost is that a repeated annotation makes a
worst-of-5 secretly a worst-of-4; `n_distinct_clicks` is recorded per ROI so the notebook can
say how often that happened instead of assuming it did not.

Outputs
-------
``results/tm_recall_workload_cells.csv``   one row per (ROI, seed): config, med/mad, pool size,
                                           n at z=1.0, deep recall, z_max.
``results/tm_recall_workload_tp_ledger.csv`` one row per (ROI, seed, mitosis): its rank in the
                                           deep pool, z, score, and its annotation's votes.
``results/tm_recall_workload_pool_z.npz``  each cell's deep-pool ``z``, ``cx``, ``cy``
                                           (score-descending), so the notebook can evaluate
                                           the retained count and `coverage_frac` at ANY
                                           cutoff exactly, with no grid interpolation.
"""

from __future__ import annotations

import gc
import os
import time

import cv2
import numpy as np
import pandas as pd
from skimage.measure import label, regionprops

from midog_utils import channels as ch
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from midog_utils import find_and_suppress as fs
from midog_utils import seed_selection as ss
from midog_utils import template_match as tm
from midog_utils.nms import nms_by_distance

# --------------------------------------------------------------------------------------
# Configuration -- byte-identical to high_z except IMAGES_DIR and the dropped od axis
# --------------------------------------------------------------------------------------
IMAGES_DIR = "images/extra_valid"
SEEDS = tuple(range(5))

CHANNEL = "hematoxylin_od"
METHOD = cv2.TM_CCOEFF
PEAK_MIN_DISTANCE = 7
SELF_HIT_RADIUS = 5.0
CURRENT_Z = 1.0                    # the operating point every saving is measured against
DEEP_FLOOR_Z = -1.5                # high_z's min(Z_LEVELS) - 0.5, pinned as a constant
MAX_PEAKS = 2_000_000

NMS_RADIUS_UM = 5.0
MATCH_RADIUS_UM = ev.MIDOG_RADIUS_UM

CFG = fs.FSConfig(channel=CHANNEL, base_size=tm.BASE_SIZE, scales=(1.0,),
                  n_angles=1, flips=(False,), peak_min_distance=PEAK_MIN_DISTANCE,
                  self_hit_radius=SELF_HIT_RADIUS)
BORDER = CFG.patch_size // 2       # 36
OTSU_WINDOW = tm.BASE_SIZE         # 51
LARGEST_CC_KW = dict(min_area=50, max_area_frac=0.85, min_solidity=0.5)

OUT_CELLS = "results/tm_recall_workload_cells.csv"
OUT_LEDGER = "results/tm_recall_workload_tp_ledger.csv"
OUT_POOLZ = "results/tm_recall_workload_pool_z.npz"


# --------------------------------------------------------------------------------------
# Seed drawing -- copied verbatim from high_z (which copied it from lcc)
# --------------------------------------------------------------------------------------
def _odd_local(n: int, minimum: int = 5) -> int:
    n = int(round(n))
    if n % 2 == 0:
        n += 1
    return max(minimum, n)


def largest_cc_box(patch, min_area=50, max_area_frac=0.85, min_solidity=0.5):
    """Otsu-threshold `patch`, return the bbox of its LARGEST connected component by area.

    (y0, y1, x0, x1), patch-local, half-open, or None. Mirrors
    `seed_selection.tighten_box_otsu`'s "binary" path with only component selection changed.
    Kept local to the caller rather than promoted into `midog_utils/seed_selection.py` for
    high_z's reason: this rule was probed and explicitly *not* adopted as production logic
    (`Research Logs/2026-08-25-bbox-tuning-domain-mismatch.md`, Experiment 5).
    """
    if patch.ndim != 2:
        raise ValueError(f"largest_cc_box needs a single-channel patch, got {patch.shape}")
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    _, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    regions = regionprops(label(binary, connectivity=2))
    if not regions:
        return None
    largest = max(regions, key=lambda r: r.area)
    if largest.area < min_area or largest.area > max_area_frac * patch.size:
        return None
    if largest.solidity < min_solidity:
        return None
    y0, x0, y1, x1 = largest.bbox
    return int(y0), int(y1), int(x0), int(x1)


def draw_seed_with_retry(pool: pd.DataFrame, rng, check_fn):
    """Draw a row via `rng.integers`; on failure drop it and redraw on the same stream."""
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
    raise ValueError("seed pool exhausted -- no candidate passed check_fn")


def roi_files(images_dir: str = IMAGES_DIR):
    return sorted(f for f in os.listdir(images_dir) if f.endswith(".tiff"))


# --------------------------------------------------------------------------------------
# Main sweep
# --------------------------------------------------------------------------------------
def main():
    t_all = time.time()
    images, annotations = ds.load_annotations()
    ds.check_invariants(annotations)
    meta = images.set_index("file_name")[["image_id", "tumor_type"]]

    files = roi_files()
    print(f"{len(files)} ROIs x {len(SEEDS)} seeds from {IMAGES_DIR}/", flush=True)

    cell_rows, tp_rows, pool_z = [], [], {}

    for fn in files:
        image_id = int(meta.loc[fn, "image_id"])
        domain = meta.loc[fn, "tumor_type"]
        t_dom = time.time()

        rgb = ds.load_roi(f"{IMAGES_DIR}/{fn}")
        mpp = ds.roi_mpp(f"{IMAGES_DIR}/{fn}")
        match_radius = ev.radius_px(mpp, MATCH_RADIUS_UM)
        nms_radius = ev.radius_px(mpp, NMS_RADIUS_UM)
        hem = ch.to_channel(rgb, CHANNEL)
        gray_inv = ch.to_gray_inverted(rgb)
        H, W = hem.shape[:2]
        roi_shape = rgb.shape
        gt = ds.image_annotations(annotations, fn)
        gt_mit = gt[gt["category_id"] == ds.MITOTIC]
        seed_pool, flagged = ss.agreement_pool(gt_mit)
        seed_pool = ss.border_filter(seed_pool, BORDER, roi_shape)
        del rgb
        gc.collect()

        def _check(cand_row, _gi=gray_inv):
            p = tm.read_padded_patch(_gi, float(cand_row["cx"]), float(cand_row["cy"]),
                                     OTSU_WINDOW)
            return None if p is None else largest_cc_box(p, **LARGEST_CC_KW)

        for si in SEEDS:
            t0 = time.time()
            rng = np.random.default_rng([si, image_id])
            seed, seed_box, n_retries = draw_seed_with_retry(seed_pool, rng, _check)
            seed_xy = (float(seed["cx"]), float(seed["cy"]))
            seed_ann_id = int(seed["ann_id"])
            y0, y1, x0, x1 = seed_box
            base_size = _odd_local(max(y1 - y0, x1 - x0), minimum=5)
            gt_eval = gt[gt["ann_id"] != seed_ann_id].reset_index(drop=True)

            patch73 = tm.read_padded_patch(hem, *seed_xy, CFG.patch_size)
            templates, _ = tm.build_augmentations(patch73, base_size, CFG.scales,
                                                  CFG.n_angles, CFG.flips)
            PAD = max((t.shape[0] - 1) // 2 for t in templates)
            hem_p = cv2.copyMakeBorder(hem, PAD, PAD, PAD, PAD, borderType=cv2.BORDER_REPLICATE)
            fused_p, best_p, valid_p = tm.fused_response(hem_p, templates,
                                                         CFG.scale_normalize, method=METHOD)
            fused = fused_p[PAD:PAD + H, PAD:PAD + W]
            valid = valid_p[PAD:PAD + H, PAD:PAD + W]
            assert bool(valid.all()), f"{fn} s{si}: padding did not make the whole ROI valid"
            del hem_p, fused_p, best_p, valid_p
            gc.collect()

            med, mad = tm.robust_stats(fused, valid)
            centers, scores = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE,
                                               med + DEEP_FLOOR_Z * mad, MAX_PEAKS)
            keep = nms_by_distance(centers, scores, nms_radius)
            centers, scores = centers[keep], scores[keep]
            ok = np.hypot(centers[:, 0] - seed_xy[0],
                          centers[:, 1] - seed_xy[1]) > CFG.self_hit_radius
            centers, scores = centers[ok], scores[ok]

            pool = pd.DataFrame({"cx": centers[:, 0], "cy": centers[:, 1], "score": scores})
            pool["z"] = (pool["score"] - med) / mad

            # The prefix identity's step 1, asserted rather than assumed: the pool must
            # already be score-descending or "rank in the deep pool" is not well defined and
            # a score cut is not a prefix.
            assert bool(np.all(np.diff(pool["score"].to_numpy()) <= 0)), \
                f"{fn} s{si}: pool is not score-descending -- the rank ledger would be void"
            # mergesort, matching `compare._rank` and `extract_peaks`' documented tie
            # convention, so this frame IS the ranking every downstream arm would produce.
            pool = pool.sort_values("score", ascending=False,
                                    kind="mergesort").reset_index(drop=True)

            det_out, _ = ev.bucket_detections(pool, gt_eval, match_radius)
            n_mit = int((gt_eval["category_id"] == ds.MITOTIC).sum())
            tp = det_out[det_out["bucket"] == ev.HUMAN_CORRECT_LABEL]
            deep_recall = len(tp) / n_mit
            n_at_current = int((pool["z"].to_numpy() >= CURRENT_Z).sum())

            # THE LEDGER. `det_out` preserves `pool`'s row order, so a TP's positional index
            # is its 0-based rank in the score-descending deep pool; +1 makes it the list
            # LENGTH a reader would have to be given to receive it, which is the cost axis.
            assert len(tp) > 0, f"{fn} s{si}: deep pool found no mitosis at all"
            # `pool` carries a RangeIndex and `bucket_detections` copies it through, so a TP's
            # index IS its 0-based position in the score-descending pool. Asserted, because
            # the whole ledger is wrong (not merely noisy) if that ever stops holding.
            ranks0 = tp.index.to_numpy()
            assert np.array_equal(ranks0, np.flatnonzero(
                (det_out["bucket"] == ev.HUMAN_CORRECT_LABEL).to_numpy())), \
                f"{fn} s{si}: TP index is not its position in the pool"
            votes = gt_eval.set_index("ann_id")
            for rank0, (_, r_tp) in zip(ranks0, tp.iterrows()):
                ann = votes.loc[int(r_tp["matched_ann_id"])]
                tp_rows.append({
                    "file_name": fn, "tumor_type": domain, "seed_index": si,
                    "matched_ann_id": int(r_tp["matched_ann_id"]),
                    "rank": int(rank0) + 1,
                    "z": float(r_tp["z"]), "score": float(r_tp["score"]),
                    "unanimous": bool(ann["unanimous"]),
                    "n_votes": int(ann["n_votes"]),
                    "n_mitotic_votes": int(ann["n_mitotic_votes"])})

            z_max = float(tp["z"].min())
            cell_rows.append({
                "file_name": fn, "tumor_type": domain, "image_id": image_id,
                "seed_index": si, "seed_ann_id": seed_ann_id, "n_retries": n_retries,
                "contested_seed_tier": bool(flagged), "base_size": base_size,
                "mpp": mpp, "roi_h": int(roi_shape[0]), "roi_w": int(roi_shape[1]),
                "match_radius_px": match_radius, "nms_radius_px": nms_radius, "pad_px": PAD,
                # med / mad / z_max carry FULL precision on purpose. z_max is a minimum, and
                # round-to-nearest moves a minimum up in about half of all cells, which then
                # drops the very mitosis it was computed to protect.
                "map_median": med, "mad_scale": mad,
                "n_pool": len(pool), "n_gt_mitotic": n_mit, "n_tp_found": len(tp),
                "deep_recall": deep_recall, "n_at_current_z": n_at_current,
                "z_max": z_max, "tp_score_min": float(tp["score"].min())})

            # The pool's z (descending) makes the retained count exact at ANY cutoff via
            # searchsorted; cx/cy make `coverage_frac` re-computable at any cutoff too.
            # float64, NOT float32. Every cutoff this ledger is read at is some cell's own
            # `z_max`, so for that cell the cutoff sits EXACTLY on a stored value. Rounding
            # the store to float32 moves it below the cutoff about half the time, and
            # `searchsorted` then drops the binding mitosis and reports recall (m-1)/m --
            # a systematic off-by-one in precisely the cell that defines the threshold.
            pool_z[f"{fn}|{si}|z"] = pool["z"].to_numpy(dtype=np.float64)
            pool_z[f"{fn}|{si}|cx"] = pool["cx"].to_numpy(dtype=np.int32)
            pool_z[f"{fn}|{si}|cy"] = pool["cy"].to_numpy(dtype=np.int32)

            print(f"[{fn}] s{si} base={base_size:2d} n_pool={len(pool):6d} "
                  f"TP={len(tp)}/{n_mit} z_max={z_max:6.3f} "
                  f"n(z=1)={n_at_current:6d} n(z_max)={int((pool['z'] >= z_max).sum()):6d} "
                  f"[{time.time() - t0:.0f}s]", flush=True)

            del fused, valid, patch73, templates, pool, det_out
            gc.collect()

        del hem, gray_inv, _check
        gc.collect()
        print(f"  -- {fn} ({domain}) done in {time.time() - t_dom:.0f}s", flush=True)

    cells = pd.DataFrame(cell_rows)
    ledger = pd.DataFrame(tp_rows)
    os.makedirs("results", exist_ok=True)
    cells.to_csv(OUT_CELLS, index=False)
    ledger.to_csv(OUT_LEDGER, index=False)
    np.savez_compressed(OUT_POOLZ, **pool_z)

    print(f"\ntotal {time.time() - t_all:.0f}s")
    print(f"  {len(cells)} (ROI, seed) cells, {len(ledger)} ledger rows")
    print(f"  deep_recall == 1.0 in {int((cells['deep_recall'] == 1.0).sum())}/{len(cells)}")
    print(f"  distinct clicks: {cells.groupby('file_name')['seed_ann_id'].nunique().sum()}"
          f" of {len(cells)} draws")
    print(f"  wrote {OUT_CELLS}, {OUT_LEDGER}, {OUT_POOLZ}")


if __name__ == "__main__":
    main()
