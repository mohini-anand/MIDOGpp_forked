"""F5: what does shrinking the NMS radius do, with border padding held fixed?

Pre-registered in `Research Logs/2026-09-04-f5-preregistration.md` (revision 3, converged
after two rounds of adversarial review). Read that document before this file; every design
choice here has a numbered justification there and this module only implements it.

The one-line summary of why this exists
---------------------------------------
`tm_threshold_axis_sweep_v2.ipynb` changed **two** things at once relative to v1 -- border
padding and a smaller NMS radius -- and says so in its own cost table: *"v1 differs from this
run by padding and radius, so this delta cannot be attributed to the radius alone."* This run
attributes it. Padding is held ON in every arm; the NMS radius is the only thing that moves.

The three arms (pre-registration section 3a)
--------------------------------------------
``r7.5``  7.5 um = `evaluate.MIDOG_RADIUS_UM`, 29.61-33.14 px -- the repo default, the only
          value satisfying `invariants.check_nms_radius`, and **the arm that has never
          actually been run**: v1 had this radius but no padding.
``r5.9``  5.9 um, 23.29-26.07 px -- the largest um radius strictly below the closest real
          annotation pair anywhere in the 14-ROI set (5.955 um on 403.tiff).
``r5.0``  5.0 um, 19.74-22.09 px -- v2's value, so every number here is comparable to v2's.

These are a **dose-response on one mechanism**, not a mechanism decomposition. Measured over
all 14 ROIs, *zero* annotation pairs lie between 5.0 and 5.9 um, so the two shrink levels are
identical with respect to un-merging neighbouring ground truth; only one mitotic pair in the
whole set (26.57 px, 245.tiff) is even close enough for the *default* radius to merge. Three
levels is the minimum that carries the registered falsifier: if the effect is the same-object
secondary peak v2 traced, it must be monotone across the dose.

What makes this "very controlled" (section 3b)
----------------------------------------------
Per (ROI, seed) the padded fused response and the deep peak pool are computed **once**, and
all three radii are applied to the *identical* ``(centers, scores)`` arrays. Nothing upstream
of `nms_by_distance` is radius-dependent -- checked in code, not docstrings: the seed draw,
`template_match.robust_stats`, the deep floor and `template_match.extract_peaks` all take no
radius, and `nms.nms_by_distance` never mutates its inputs. So the arms differ by exactly one
argument and every comparison is a within-(ROI, seed, z, axis) paired difference.

Two consequences that are easy to get wrong and are therefore load-bearing here:

* **The shared arrays are never rebound.** v2's single-arm loop does
  ``centers, scores = centers[keep], scores[keep]``; doing that here would make arm 2 suppress
  arm 1's survivors. Each arm reads the shared pool into new local names.
* **`od` is hoisted out of the arm loop**, because chromatin density depends only on (cx, cy)
  and the image, so it is identical for any point shared between arms. Scoring it once on the
  deep pool and subsetting turns the extra `r5.9` arm into a net *saving*. The equivalence is
  asserted rather than assumed (`--gates`, gate 7), following F4 section 9 gate 7.

Outputs
-------
``results/f5_nms_radius_ablation.csv``               the tidy long frame from `compare.evaluate_arms`
``results/f5_nms_radius_ablation_verification.csv``  every invariant and gate record
``results/f5_nms_radius_ablation_ledger.csv``        per-(ROI, seed, radius) win/loss + gap + churn
``results/f5_nms_radius_ablation_spacing.csv``       the annotation-spacing geometry, section 2a

Run
---
    python f5_nms_radius_ablation.py --spacing   # section 2a's table alone, ~10 s
    python f5_nms_radius_ablation.py --gates     # pre-run gates, one ROI, ~40 s
    python f5_nms_radius_ablation.py             # the full run, ~9-12 min
    python f5_nms_radius_ablation.py --smoke     # one ROI, one seed, ~20 s

Analysis lives in `f5_nms_radius_ablation.ipynb`, which reads these CSVs and does no matching
-- the same split as `tm_variant_sweep.py` / `tm_variant_report.py`, for the same reason.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from sklearn.neighbors import KDTree

from midog_utils import channels as ch
from midog_utils import chromatin as cm
from midog_utils import compare as cp
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from midog_utils import find_and_suppress as fs
from midog_utils import invariants as inv
from midog_utils import template_match as tm
from midog_utils.nms import nms_by_distance
from tm_variant_sweep import draw_seeds

# ---------------------------------------------------------------------------------------
# Configuration -- every value traceable to a numbered section of the pre-registration
# ---------------------------------------------------------------------------------------

IMAGES_DIR = "images/extra_valid"     # section 3c: 2 ROIs per tumour type, all 7 domains
N_SEEDS = 5                           # section 3c; `tm_variant_sweep.N_SEEDS`

CHANNEL = "hematoxylin_od"
METHOD = cv2.TM_CCOEFF
PEAK_MIN_DISTANCE = 7
SELF_HIT_RADIUS = 5.0

# section 3a -- the one factor. Physically scaled in um so each is one physical distance on
# scanners whose mpp runs 0.2263-0.2533, the same way `evaluate.radius_px` is already defined.
CONTROL_UM = ev.MIDOG_RADIUS_UM       # 7.5
RADII_UM = (7.5, 5.9, 5.0)
RADIUS_TAG = {7.5: "r7.5", 5.9: "r5.9", 5.0: "r5.0"}

MATCH_RADIUS_UM = ev.MIDOG_RADIUS_UM  # scoring radius, pinned in every arm -- never varies

# section 3d. z is a nuisance axis. -1.0/-0.5/0.0 are dropped as coverage-saturated; 2.5 and
# 3.0 are kept because they are the ONLY region where coverage_frac is not saturated
# (0.26-0.63 against 0.88-0.99 at z <= 1.0) and so the only place full-list recall can carry
# evidence at all.
Z_LEVELS = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
HEADLINE_Z = 1.0                      # the operating point set in tp_fp_score_distribution.ipynb

# section 3d: pinned as a CONSTANT, deliberately decoupled from min(Z_LEVELS). v2 defines it
# as min(Z_LEVELS) - 0.5; keeping that coupling while trimming the z grid would silently move
# the deep pool and break the seed-0 reproduction gate (section 7.9). The decoupling is free
# for the same reason the one-match-many-z shortcut works: `nms_by_distance` walks descending
# score and writes `suppressed[]` only from already-visited, strictly higher-scoring indices,
# so a peak at -1.4 can never suppress a peak at +0.6 -- every peak below min(Z_LEVELS) is
# inert with respect to the kept set at every reported z.
DEEP_FLOOR_Z = -1.5
MAX_PEAKS = 2_000_000

AXES = ("tm_score", "chromatin_od")
AXIS_RANK_KEY = {"tm_score": "score", "chromatin_od": "od"}
PRIMARY_AXIS = "chromatin_od"         # section 3d: the shipped ranker (commit 7c3af93)

# section 6.1: a strict superset of `compare.BUDGETS`, so rows stay comparable with every
# earlier CSV. Dense because it costs nothing (pure indexing into tp_cum) and because the
# prior 2-ROI sweep predicts a crossover a sparse grid would step over.
BUDGETS_F5 = (10, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500,
              750, 1000, 1500, 2000, 3000, 5000, 10000)
HEADLINE_K = 250

OD_PAD = tm.BASE_SIZE // 2            # 25 -- the chromatin window's OWN margin (v2 cell 14)

CFG = fs.FSConfig(channel=CHANNEL, base_size=tm.BASE_SIZE, scales=(1.0,), n_angles=1,
                  flips=(False,), peak_min_distance=PEAK_MIN_DISTANCE,
                  self_hit_radius=SELF_HIT_RADIUS)
BORDER = CFG.patch_size // 2          # 36 -- governs which annotation may be the SEED

OUT = Path("results")
OUT_MAIN = OUT / "f5_nms_radius_ablation.csv"
OUT_VERIF = OUT / "f5_nms_radius_ablation_verification.csv"
OUT_LEDGER = OUT / "f5_nms_radius_ablation_ledger.csv"
OUT_SPACING = OUT / "f5_nms_radius_ablation_spacing.csv"

V2_CSV = "results/tm_ccoeff_threshold_axis_sweep_v2.csv"


# ---------------------------------------------------------------------------------------
# Section 2a -- the annotation-spacing geometry
# ---------------------------------------------------------------------------------------

def roi_files(images_dir: str = IMAGES_DIR):
    return sorted(f for f in os.listdir(images_dir) if f.endswith(".tiff"))


def _min_pair_distance(df: pd.DataFrame) -> float:
    """Closest pair within one annotation frame, in pixels (inf when fewer than 2)."""
    if len(df) < 2:
        return float("inf")
    xy = df[["cx", "cy"]].to_numpy()
    d, _ = KDTree(xy).query(xy, k=2)
    return float(d[:, 1].min())


def _pairs_within(df: pd.DataFrame, radius: float) -> int:
    if len(df) < 2:
        return 0
    xy = df[["cx", "cy"]].to_numpy()
    n = KDTree(xy).query_radius(xy, r=radius, count_only=True)
    return int((n.sum() - len(xy)) // 2)


def spacing_table(images, annotations, images_dir: str = IMAGES_DIR) -> pd.DataFrame:
    """Section 2a. Closest annotation pairs, and how many fall inside each candidate radius.

    Reported in **microns** as well as pixels. Microns is the only unit in which a um-scaled
    radius can be compared across scanners: comparing one ROI's px spacing against another
    ROI's px radius is a units error, and it is the one revision 1 of the pre-registration
    made (it read 5.945 um where the binding constraint is 5.955).
    """
    rows = []
    for fn in roi_files(images_dir):
        im = images[images["file_name"] == fn].iloc[0]
        mpp = ds.roi_mpp(f"{images_dir}/{fn}")
        gt = ds.image_annotations(annotations, fn)
        mit = gt[gt["category_id"] == ds.MITOTIC]
        row = {"file_name": fn, "tumor_type": im["tumor_type"], "mpp": round(mpp, 5),
               "n_mitotic": int(len(mit)), "n_annotations": int(len(gt))}
        d_mit, d_all = _min_pair_distance(mit), _min_pair_distance(gt)
        row["min_mit_mit_px"] = round(d_mit, 2)
        row["min_any_any_px"] = round(d_all, 2)
        row["min_mit_mit_um"] = round(d_mit * mpp, 4)
        row["min_any_any_um"] = round(d_all * mpp, 4)
        for um in RADII_UM:
            tag, r = RADIUS_TAG[um], ev.radius_px(mpp, um)
            row[f"{tag}_px"] = round(r, 2)
            row[f"{tag}_pairs_mit"] = _pairs_within(mit, r)
            row[f"{tag}_pairs_any"] = _pairs_within(gt, r)
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------
# The shared per-(ROI, seed) work -- everything radius-independent, computed exactly once
# ---------------------------------------------------------------------------------------

def build_pool(hem, seed_xy):
    """The padded response, the deep peak pool, and `od` -- all radius-independent.

    Returns ``(centers, scores, od, meta)``. `od` is scored **once here** on the deep pool
    and subset per arm downstream: chromatin density depends only on (cx, cy) and the image,
    so it is identical for any point shared between arms, and scoring it per radius would
    recompute identical numbers three times.
    """
    H, W = hem.shape[:2]
    patch = tm.read_padded_patch(hem, seed_xy[0], seed_xy[1], CFG.patch_size)
    templates, _ = tm.build_augmentations(patch, CFG.base_size, CFG.scales,
                                          CFG.n_angles, CFG.flips)

    # Padding, held ON in every arm. PAD is derived from the bank, not hardcoded, so it stays
    # correct if BASE_SIZE or the augmentation scales ever change. With an odd template side
    # t and PAD = (t-1)//2, matchTemplate on a PAD-padded array returns a map of shape
    # (H + 2*PAD - t + 1, W + ...) = (H, W), which `fused_response` offsets to (PAD, PAD) --
    # so the crop below is all-valid by construction, asserted per cell rather than trusted.
    PAD = max((t.shape[0] - 1) // 2 for t in templates)
    hem_p = cv2.copyMakeBorder(hem, PAD, PAD, PAD, PAD, borderType=cv2.BORDER_REPLICATE)
    fused_p, _, valid_p = tm.fused_response(hem_p, templates, CFG.scale_normalize, method=METHOD)
    fused = fused_p[PAD:PAD + H, PAD:PAD + W]
    valid = valid_p[PAD:PAD + H, PAD:PAD + W]
    del hem_p, fused_p, valid_p

    med, mad = tm.robust_stats(fused, valid)
    deep_floor = med + DEEP_FLOOR_Z * mad
    centers, scores = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE, deep_floor, MAX_PEAKS)

    # `od` needs its OWN padding: `chromatin.score_detections` reads a fixed 51 px window, so
    # every candidate the border padding newly admits would score NaN and `compare._rank`
    # would sort it last. This is v2's own bug, which overstated its cost table by ~116x.
    hem_od = cv2.copyMakeBorder(hem, OD_PAD, OD_PAD, OD_PAD, OD_PAD,
                                borderType=cv2.BORDER_REPLICATE)
    shifted = pd.DataFrame({"cx": centers[:, 0] + OD_PAD, "cy": centers[:, 1] + OD_PAD})
    od = cm.score_detections(shifted, hem_od)["od"].to_numpy()
    od_nan_naive = int(pd.isna(cm.score_detections(
        pd.DataFrame({"cx": centers[:, 0], "cy": centers[:, 1]}), hem)["od"]).sum())
    del hem_od

    meta = {"PAD": PAD, "valid_coverage": float(valid.mean()), "valid_all": bool(valid.all()),
            "map_median": float(med), "mad_scale": float(mad), "deep_floor": float(deep_floor),
            "n_peaks_deep": int(len(centers)), "od_nan_naive": od_nan_naive,
            "od_nan_fixed": int(np.isnan(od).sum())}
    return centers, scores, od, meta, (fused, valid)


def apply_radius(centers, scores, od, radius, seed_xy):
    """One arm's candidate pool: NMS at ``radius``, then the self-hit filter.

    Order matters and is pinned to v2's (NMS first, self-hit second) so the comparison to v2
    holds. That order leaves a detection-free hole of radius ~= ``radius`` around the seed,
    which shrinks with the radius and could in principle mimic a recall gain -- but measured
    over all 14 ROIs x 5 seeds, **zero** evaluation annotations fall inside that hole even at
    the largest radius, so the confound cannot fire here. Asserted per cell in `run_roi`.

    Never rebinds the caller's arrays: the shared pool must survive for the next arm.
    """
    keep = nms_by_distance(centers, scores, radius)
    c_keep, s_keep, od_keep = centers[keep], scores[keep], od[keep]
    ok = np.hypot(c_keep[:, 0] - seed_xy[0], c_keep[:, 1] - seed_xy[1]) > CFG.self_hit_radius
    pool = pd.DataFrame({"cx": c_keep[ok, 0], "cy": c_keep[ok, 1],
                         "score": s_keep[ok], "od": od_keep[ok]})
    return pool, c_keep


# ---------------------------------------------------------------------------------------
# Section 6.5 / 6.6 -- the ledger: win/loss, the same-object gap, and churn
# ---------------------------------------------------------------------------------------

def suppressor_index(centers, scores, radius):
    """For each peak, the index of the kept peak that suppressed it (itself, if kept).

    Re-traces `nms_by_distance`'s own greedy order rather than a proxy -- the same walk v2
    cell 10 does, generalised to any radius.
    """
    tree = KDTree(centers)
    neighbours = tree.query_radius(centers, r=radius)
    order = np.argsort(-scores, kind="stable")
    suppressed = np.zeros(len(centers), dtype=bool)
    sup_of = np.full(len(centers), -1, dtype=int)
    for idx in order:
        if suppressed[idx]:
            continue
        for nb in neighbours[idx]:
            if nb != idx and not suppressed[nb]:
                sup_of[nb] = idx
        suppressed[neighbours[idx]] = True
        sup_of[idx] = idx
    return sup_of


def same_object_gaps(centers, scores, radius, gt_mit_eval):
    """Section 6.5. Distance from each mitotic annotation's nearest deep-pool local maximum
    to the peak that suppressed it, at this radius.

    This is the quantity the whole dose grid is positioned against, and until now it was
    measured on **one ROI and one seed** (v2 cells 9-11, the 23-29 px band). Computed once
    per (ROI, seed, radius): which peak suppressed which is entirely upstream of ranking, so
    it is axis-free and running it per axis would duplicate identical numbers.
    """
    if len(centers) == 0 or len(gt_mit_eval) == 0:
        return np.zeros(0)
    tree = KDTree(centers)
    sup_of = suppressor_index(centers, scores, radius)
    _, i = tree.query(gt_mit_eval[["cx", "cy"]].to_numpy(), k=1)
    own = i[:, 0]
    sup = sup_of[own]
    gaps = np.hypot(centers[sup, 0] - centers[own, 0], centers[sup, 1] - centers[own, 1])
    return gaps[sup != own]          # 0 gap means the peak was itself kept -- not a gap


def ledger_row(pool, gt_eval, match_radius, axis, ctrl_top=None):
    """Found-annotation identities, duplicate accounting, and top-K churn for one arm.

    ``ctrl_top`` is the control arm's top-K (cx, cy) set at the same cell; when given, the
    Jaccard overlap against it is returned as ``topk_churn``. Without a churn measure a null
    at K is uninterpretable -- an inert intervention and an active-but-neutral one look the
    same (F4 section 6's rationale, here as two columns off a table already being built).
    """
    ranked = pool.sort_values(AXIS_RANK_KEY[axis], ascending=False, na_position="last",
                              kind="mergesort").reset_index(drop=True)
    det, gt_out = ev.bucket_detections(ranked, gt_eval, match_radius)

    found = set(gt_out[(gt_out["category_id"] == ds.MITOTIC) & gt_out["found"]]["ann_id"])
    gt_mit_out = gt_out[gt_out["category_id"] == ds.MITOTIC]
    claim = {int(a): int(r) for a, r in zip(gt_mit_out["ann_id"].to_numpy(),
                                            gt_mit_out["matched_rank"].to_numpy()) if r >= 0}

    # Duplicate FPs, reported WITH their denominator: a second detection inside an already
    # claimed object's match circle is bucketed as an unannotated false positive, and the raw
    # count rises mechanically with list length, so it is meaningless without n_within.
    if len(det) and len(gt_eval):
        d_gt, _ = KDTree(gt_eval[["cx", "cy"]].to_numpy()).query(
            det[["cx", "cy"]].to_numpy(), k=1)
        inside = d_gt[:, 0] <= match_radius
        dup = inside & (det["bucket"] == ev.NON_HUMAN_FINDINGS).to_numpy()
        n_within, n_dup = int(inside.sum()), int(dup.sum())
    else:
        n_within = n_dup = 0

    top = ranked.head(HEADLINE_K)[["cx", "cy"]]
    top_set = set(map(tuple, np.round(top.to_numpy(), 3)))
    churn = float("nan")
    if ctrl_top is not None:
        union = len(top_set | ctrl_top)
        churn = float(1.0 - len(top_set & ctrl_top) / union) if union else float("nan")

    return {"found": found, "claim": claim, "n_within_match_radius": n_within,
            "n_dup_fp": n_dup, "n_detections": int(len(ranked)),
            "topk_churn": churn, "top_set": top_set}


# ---------------------------------------------------------------------------------------
# The main loop
# ---------------------------------------------------------------------------------------

def run_roi(fn, images, annotations, images_dir, n_seeds, checks, ledger_rows,
            verify_shortcut, quiet=False):
    """One ROI: load once, then n_seeds x 3 radii off one response map per seed."""
    im = images[images["file_name"] == fn].iloc[0]
    image_id, domain = int(im["image_id"]), im["tumor_type"]
    path = f"{images_dir}/{fn}"

    t_roi = time.time()
    rgb = ds.load_roi(path)
    mpp = ds.roi_mpp(path)
    hem = ch.to_channel(rgb, CHANNEL)          # hoisted: 3 s per ROI, not per seed
    roi_shape = rgb.shape
    match_radius = ev.radius_px(mpp, MATCH_RADIUS_UM)
    radii = {RADIUS_TAG[um]: ev.radius_px(mpp, um) for um in RADII_UM}

    gt = ds.image_annotations(annotations, fn)
    gt_mit = gt[gt["category_id"] == ds.MITOTIC]
    seeds, _, seed_records = draw_seeds(gt_mit, BORDER, roi_shape, n_seeds)
    checks.append(inv.check_distinct_seeds(seed_records, label=fn))

    frames = []
    for si, seed in seeds:
        t0 = time.time()
        seed_xy = (float(seed["cx"]), float(seed["cy"]))
        gt_eval = gt[gt["ann_id"] != int(seed["ann_id"])].reset_index(drop=True)
        gt_mit_eval = gt_eval[gt_eval["category_id"] == ds.MITOTIC]

        centers, scores, od, meta, (fused, valid) = build_pool(hem, seed_xy)

        # -- section 7.1/7.2/7.6: the running checks on the shared, radius-independent work
        assert meta["valid_all"], f"{fn}/si{si}: padding did not make the whole ROI valid"
        assert meta["od_nan_fixed"] == 0, (
            f"{fn}/si{si}: {meta['od_nan_fixed']} candidates still NaN after the od padding")
        # Applied to the PRE-NMS pool: `compare.evaluate_arms` checks the post-NMS,
        # post-z-filter length, which against a 2,000,000 cap cannot fail, so on its own it
        # records a pass it could not have withheld.
        checks.append(inv.check_no_cap(meta["n_peaks_deep"], (MAX_PEAKS,),
                                       label=f"{fn}/si{si}/deep_pool"))

        # -- section 7.10: the self-hit hole, pre-measured at 0/70 and asserted anyway so a
        # future change to the seed rule cannot silently re-open the confound.
        d_seed = np.hypot(gt_eval["cx"] - seed_xy[0], gt_eval["cy"] - seed_xy[1])
        n_hole = int((d_seed <= max(radii.values())).sum())
        assert n_hole == 0, (
            f"{fn}/si{si}: {n_hole} evaluation annotations sit inside the self-hit hole; the "
            "confound the pre-registration measured at zero has re-opened")

        ctx_base = {"file_name": fn, "image_id": image_id, "tumor_type": domain,
                    "seed_index": si, "seed_ann_id": int(seed["ann_id"]),
                    "mpp": round(mpp, 5), "map_median": round(meta["map_median"], 5),
                    "mad_scale": round(meta["mad_scale"], 5),
                    "deep_floor": round(meta["deep_floor"], 5),
                    "n_peaks_deep": meta["n_peaks_deep"], "pad_px": meta["PAD"],
                    "match_radius_px": round(match_radius, 3),
                    "n_gt_within_seed_hole": n_hole}

        ctrl_tops, ctrl_found, ctrl_claim = {}, None, None
        pool_sizes = []
        for um in RADII_UM:                       # control FIRST, so churn has its reference
            tag, radius = RADIUS_TAG[um], radii[RADIUS_TAG[um]]
            pool, kept_centers = apply_radius(centers, scores, od, radius, seed_xy)
            pool_sizes.append(len(pool))

            # section 7.4: direct proof the suppression claimed is the suppression that ran
            checks.append(inv.check_min_separation(kept_centers, radius,
                                                   label=f"{fn}/si{si}/{tag}"))

            # ---- side diagnostics at the headline z (sections 6.4-6.6) ----
            cut_h = meta["map_median"] + HEADLINE_Z * meta["mad_scale"]
            pool_h = pool[pool["score"] >= cut_h]
            gaps = same_object_gaps(centers, scores, radius, gt_mit_eval)
            for axis in AXES:
                led = ledger_row(pool_h, gt_eval, match_radius, axis,
                                 ctrl_tops.get(axis) if um != CONTROL_UM else None)
                if um == CONTROL_UM:
                    ctrl_tops[axis] = led["top_set"]
                    if axis == PRIMARY_AXIS:
                        ctrl_found, ctrl_claim = led["found"], led["claim"]
                gained = lost = claim_change = -1
                if (um != CONTROL_UM and axis == PRIMARY_AXIS
                        and ctrl_found is not None and ctrl_claim is not None):
                    gained = len(led["found"] - ctrl_found)
                    lost = len(ctrl_found - led["found"])
                    claim_change = sum(1 for a, r in led["claim"].items()
                                       if ctrl_claim.get(a, -1) != r)
                ledger_rows.append({
                    **{k: ctx_base[k] for k in ("file_name", "tumor_type", "seed_index")},
                    "nms_radius_um": um, "nms_radius_px": round(radius, 3), "axis": axis,
                    "z": HEADLINE_Z, "n_detections": led["n_detections"],
                    "n_gt_mitotic": int(len(gt_mit_eval)),
                    "n_found": len(led["found"]), "gained_vs_control": gained,
                    "lost_vs_control": lost, "gt_claim_change": claim_change,
                    "n_within_match_radius": led["n_within_match_radius"],
                    "n_dup_fp": led["n_dup_fp"],
                    "dup_fp_frac": round(led["n_dup_fp"] / max(led["n_detections"], 1), 5),
                    "topk_churn": led["topk_churn"],
                    "gap_n": int(len(gaps)),
                    "gap_median_px": round(float(np.median(gaps)), 2) if len(gaps) else np.nan,
                    "gap_p05_px": round(float(np.percentile(gaps, 5)), 2) if len(gaps) else np.nan,
                    "gap_p95_px": round(float(np.percentile(gaps, 95)), 2) if len(gaps) else np.nan,
                })

            # ---- section 7.5: the one-match-many-z shortcut, si == 0 only ----
            if verify_shortcut and si == 0:
                for z in (HEADLINE_Z, max(Z_LEVELS)):
                    cut = meta["map_median"] + z * meta["mad_scale"]
                    c2, s2 = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE, cut, MAX_PEAKS)
                    k2 = nms_by_distance(c2, s2, radius)
                    c2 = c2[k2]
                    ok2 = np.hypot(c2[:, 0] - seed_xy[0],
                                   c2[:, 1] - seed_xy[1]) > CFG.self_hit_radius
                    direct = c2[ok2]
                    filt = pool.loc[pool["score"] >= cut, ["cx", "cy"]].to_numpy()
                    a = direct[np.lexsort((direct[:, 1], direct[:, 0]))] if len(direct) else direct
                    b = filt[np.lexsort((filt[:, 1], filt[:, 0]))] if len(filt) else filt
                    passed = (a.shape == b.shape) and bool(np.allclose(a, b))
                    checks.append({"check": "shortcut", "label": f"{fn}/si{si}/{tag}/z{z}",
                                   "n_direct": len(direct), "n_filtered": len(filt),
                                   "passed": passed})
                    assert passed, f"{fn}/si{si}/{tag}: shortcut failed at z={z}"

            # ---- the arms ----
            ctx = dict(ctx_base, nms_radius_um=um, nms_radius_tag=tag,
                       nms_radius_px=round(radius, 3), n_pool=len(pool))
            arms = []
            for z in Z_LEVELS:
                cut = meta["map_median"] + z * meta["mad_scale"]
                sub = pool[pool["score"] >= cut]
                limited = bool(cut < meta["deep_floor"])
                for axis in AXES:
                    arms.append(cp.Arm(
                        f"{axis}@{tag}", (lambda d=sub: d), rank_key=AXIS_RANK_KEY[axis],
                        seeded=True, z=z, z_dependent=True, floor_limited=limited,
                        # The control passes the invariant as a POSITIVE control that it
                        # really is the repo default; the shrunk arms pass None because the
                        # check is deliberately inapplicable, and record the radius in extra.
                        nms_radius=(radius if um == CONTROL_UM else None),
                        caps=(MAX_PEAKS,),
                        coverage_key=f"{fn}|{si}|{tag}|{z}",   # set, not sort order
                        extra={"axis": axis, "z_cut": round(cut, 5),
                               "nms_radius_um": um, "nms_radius_tag": tag}))

            frames.append(cp.evaluate_arms(arms, gt_eval, match_radius, roi_shape=roi_shape,
                                           mpp=mpp, budgets=BUDGETS_F5, context=ctx,
                                           coverage_cache={}, checks=checks))

        if not quiet:
            print(f"  [{fn} si{si}] deep={meta['n_peaks_deep']:6d} "
                  f"pools={'/'.join(str(n) for n in pool_sizes)} "
                  f"od_nan {meta['od_nan_naive']}->0 [{time.time() - t0:.0f}s]", flush=True)
        del centers, scores, od, fused, valid

    del rgb, hem
    if not quiet:
        print(f"[{fn}] {domain:34s} done [{time.time() - t_roi:.0f}s]", flush=True)
    return frames


def reproduction_gate(results: pd.DataFrame) -> pd.DataFrame:
    """Section 7.9. `r5.0` at seed 0 must reproduce v2's CSV on the 7 overlapping ROIs.

    `draw_seeds(s=0)` is byte-identical to v2's inline draw (agreement_pool -> border_filter(36)
    -> default_rng([0, image_id]) -> .iloc[integers(...)]), and DEEP_FLOOR_Z is pinned at v2's
    -1.5, so this must hold exactly at every shared (z, budget). It is the only check that can
    catch a silent divergence in the whole pipeline rather than in one function.

    One mapping is needed: v2's ``arm`` column is the axis, F5's is ``{axis}@{radius}``.
    """
    if not Path(V2_CSV).exists():
        return pd.DataFrame([{"check": "v2_reproduction", "label": V2_CSV,
                              "passed": False, "note": "v2 CSV absent"}])
    v2 = pd.read_csv(V2_CSV)
    rows = []
    mine = results[(results["seed_index"] == 0) & (results["nms_radius_tag"] == "r5.0")]
    for axis in AXES:
        a = (v2[v2["arm"] == axis][["file_name", "z", "budget",
                                    "recall_at_budget", "n_detections", "full_list_recall"]]
             .drop_duplicates(["file_name", "z", "budget"]))
        b = (mine[mine["arm"] == f"{axis}@r5.0"][["file_name", "z", "budget",
                                                  "recall_at_budget", "n_detections",
                                                  "full_list_recall"]]
             .drop_duplicates(["file_name", "z", "budget"]))
        m = a.merge(b, on=["file_name", "z", "budget"], suffixes=("_v2", "_f5"))
        if not len(m):
            rows.append({"check": "v2_reproduction", "label": axis, "n_cells": 0,
                         "passed": False, "note": "no overlapping cells"})
            continue
        same_n = bool((m["n_detections_v2"] == m["n_detections_f5"]).all())
        same_r = bool(np.allclose(m["recall_at_budget_v2"], m["recall_at_budget_f5"],
                                  atol=1e-9, equal_nan=True))
        rows.append({"check": "v2_reproduction", "label": axis, "n_cells": int(len(m)),
                     "n_rois": int(m["file_name"].nunique()),
                     "max_abs_recall_diff": float(
                         np.nanmax(np.abs(m["recall_at_budget_v2"] - m["recall_at_budget_f5"]))),
                     "max_abs_ndet_diff": int(
                         (m["n_detections_v2"] - m["n_detections_f5"]).abs().max()),
                     "passed": bool(same_n and same_r)})
    return pd.DataFrame(rows)


def od_hoist_gate(images, annotations, images_dir=IMAGES_DIR) -> dict:
    """Gate 7 (F4 section 9 gate 7, borrowed with the optimisation it licenses).

    `od` scored once on the deep pool and subset must equal `od` scored after NMS. If this
    fails the hoist in `build_pool` is invalid and every chromatin-ranked number is wrong.
    """
    fn = roi_files(images_dir)[0]
    rgb = ds.load_roi(f"{images_dir}/{fn}")
    mpp = ds.roi_mpp(f"{images_dir}/{fn}")
    hem = ch.to_channel(rgb, CHANNEL)
    gt = ds.image_annotations(annotations, fn)
    gt_mit = gt[gt["category_id"] == ds.MITOTIC]
    seeds, _, _ = draw_seeds(gt_mit, BORDER, rgb.shape, 1)
    _, seed = seeds[0]
    seed_xy = (float(seed["cx"]), float(seed["cy"]))
    centers, scores, od, _, _ = build_pool(hem, seed_xy)

    radius = ev.radius_px(mpp, 5.0)
    keep = nms_by_distance(centers, scores, radius)
    hoisted = od[keep]

    hem_od = cv2.copyMakeBorder(hem, OD_PAD, OD_PAD, OD_PAD, OD_PAD,
                                borderType=cv2.BORDER_REPLICATE)
    after = cm.score_detections(
        pd.DataFrame({"cx": centers[keep, 0] + OD_PAD, "cy": centers[keep, 1] + OD_PAD}),
        hem_od)["od"].to_numpy()
    passed = bool(np.allclose(hoisted, after, equal_nan=True))
    return {"check": "od_hoist_equivalence", "label": fn, "n": int(len(keep)),
            "max_abs_diff": float(np.nanmax(np.abs(hoisted - after))) if len(keep) else 0.0,
            "passed": passed}


# ---------------------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spacing", action="store_true", help="section 2a's table only")
    ap.add_argument("--gates", action="store_true", help="pre-run gates only")
    ap.add_argument("--smoke", action="store_true", help="one ROI, one seed")
    ap.add_argument("--images-dir", default=IMAGES_DIR)
    ap.add_argument("--seeds", type=int, default=N_SEEDS)
    args = ap.parse_args(argv)

    t00 = time.time()
    OUT.mkdir(exist_ok=True)
    images, annotations = ds.load_annotations()
    ds.check_invariants(annotations)

    sp = spacing_table(images, annotations, args.images_dir)
    sp.to_csv(OUT_SPACING, index=False)
    print(f"wrote {OUT_SPACING} {sp.shape}")
    print(f"  closest any-any pair : {sp['min_any_any_um'].min():.4f} um "
          f"({sp.loc[sp['min_any_any_um'].idxmin(), 'file_name']})")
    print(f"  closest mit-mit pair : {sp['min_mit_mit_um'].min():.4f} um "
          f"({sp.loc[sp['min_mit_mit_um'].idxmin(), 'file_name']})")
    for um in RADII_UM:
        tag = RADIUS_TAG[um]
        print(f"  {tag}: {sp[f'{tag}_px'].min():.2f}-{sp[f'{tag}_px'].max():.2f} px, "
              f"GT pairs inside: {sp[f'{tag}_pairs_mit'].sum()} mit / "
              f"{sp[f'{tag}_pairs_any'].sum()} any")
    if args.spacing:
        return 0

    if args.gates:
        g = od_hoist_gate(images, annotations, args.images_dir)
        print(f"gate od_hoist_equivalence: passed={g['passed']} "
              f"n={g['n']} max_abs_diff={g['max_abs_diff']:.3e}")
        return 0 if g["passed"] else 1

    files = roi_files(args.images_dir)
    n_seeds = args.seeds
    if args.smoke:
        files, n_seeds = files[:1], 1

    checks, ledger_rows, frames = [], [], []
    checks.append(od_hoist_gate(images, annotations, args.images_dir))
    assert checks[-1]["passed"], "od hoist equivalence gate failed -- the optimisation is invalid"

    for fn in files:
        frames += run_roi(fn, images, annotations, args.images_dir, n_seeds,
                          checks, ledger_rows, verify_shortcut=True)

    results = pd.concat(frames, ignore_index=True)
    cp.assert_floor_not_limiting(results)
    results.to_csv(OUT_MAIN, index=False)
    pd.DataFrame(ledger_rows).to_csv(OUT_LEDGER, index=False)

    verif = pd.concat([pd.DataFrame(checks), reproduction_gate(results)], ignore_index=True)
    verif.to_csv(OUT_VERIF, index=False)

    print(f"\nwrote {OUT_MAIN} {results.shape}")
    print(f"wrote {OUT_LEDGER} ({len(ledger_rows)} rows)")
    print(f"wrote {OUT_VERIF} ({len(verif)} rows)")
    failed = verif[~verif["passed"].astype(bool)] if "passed" in verif else verif.iloc[:0]
    print(f"checks: {len(verif) - len(failed)}/{len(verif)} passed")
    if len(failed):
        print(failed.to_string())
    print(f"total {time.time() - t00:.0f}s")
    return 1 if len(failed) else 0


if __name__ == "__main__":
    sys.exit(main())
