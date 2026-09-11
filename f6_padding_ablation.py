"""F6: what does border padding actually do, and is the effect bigger than the click?

Pre-registered in `Research Logs/2026-09-08-f6-preregistration.md` (revision 4, CONVERGED
after three rounds of adversarial review). Read that document before this file; every
constant and every column here has a numbered justification there and this module only
implements it. Where the plan is silent this file says so in a comment rather than
choosing quietly.

The one-line summary of why this exists
---------------------------------------
`tm_threshold_axis_sweep_v2.ipynb` turned border padding ON and shrank the NMS radius in
the same commit. F5 attributed the radius half. This attributes the padding half, and --
because both of the mechanisms padding was *justified* by are pre-measured at zero
(section 2b: reachability ceiling 0 of 7 annotations, re-normalisation <= 2.27e-3 MAD) --
it also has to name the mechanism that actually moves the numbers.

The three arms (section 3)
--------------------------
``pad_off``            unpadded response; `valid` is the interior rectangle only. This is
                       v1's configuration, and gate A re-derives v1's committed CSV from it.
``pad_on_dropband``    the padded response, **the same** extraction, NMS and `od` scoring as
                       `pad_on`, with border-band candidates then dropped from the ranked
                       list. Border peaks competed during grey dilation and during NMS but
                       none of them is in the list.
``pad_on``             the padded response, cropped back, `valid.all()`. The shippable
                       configuration, and gate B asserts it is F5's `r7.5` arm exactly.

Read across, the contrasts isolate the two live mechanisms (section 3):

* ``pad_on_dropband - pad_off``   = suppression of interior peaks by border peaks, plus the
                                    inert re-normalisation channel (`pad_on_dropband`
                                    inherits `pad_on`'s med/mad/deep_floor/z_cut).
* ``pad_on - pad_on_dropband``    = the band candidates' own contribution: the TPs they claim
                                    minus the reading slots they consume.
* ``pad_on - pad_off``            = the total, and the **primary contrast** -- the only one
                                    between two configurations that could ship.

Three implementation facts that are load-bearing and easy to get wrong
---------------------------------------------------------------------
1. **`pad_on_dropband` is one boolean mask on `pad_on`'s candidate frame.** No second
   response map, no second extraction, no second NMS, no second `od` scoring. Section 3's
   whole point is that the arm costs none of those; draft 1's third arm cost all four and
   measured a channel pinned at zero.
2. **There really are two response maps per seed** -- `fused_response` on the padded image
   *and* on the unpadded image. The unpadded one is **not** emulated by masking the padded
   one. The binding reason is gate A, not the grey dilation (section 3): the padded map's
   interior differs from the unpadded map's by up to 3.8e-6 because OpenCV's DFT block
   tiling depends on image size, and gate A compares `n_detections` on exact integer
   equality -- a `_FLOOR`-written emulation still loses 1 peak of 20,855 on 459.tiff. The
   dilation hazard is separately real (a naive `valid`-mask emulation loses 228-267 peaks)
   but is defeatable with `np.where(valid, fused, _FLOOR)`, so it is not what forces the
   second `matchTemplate`. Recorded here so a future reader does not delete it on the
   strength of an argument that does not hold.
3. **The shared arrays are never rebound.** `pad_on` and `pad_on_dropband` read the same
   post-NMS pool; `pool_on[mask]` produces a new frame and nothing mutates in place. F5's
   docstring has the failure mode this rule exists to prevent.

Held constant in every arm (section 3, last paragraph): NMS radius = scoring radius =
7.5 um, so `invariants.check_nms_radius` passes on all three; `hematoxylin_od`;
`TM_CCOEFF`; a single 51 px template; `peak_min_distance = 7`; `self_hit_radius = 5.0`
after NMS; `DEEP_FLOOR_Z = -1.5` pinned as a constant; the seed draw; the ground truth;
and `chromatin.score_detections` on the `OD_PAD = 25` padded channel in **all three** arms.

Outputs (section 10)
--------------------
``results/f6_padding_ablation.csv``               the tidy long frame from `compare.evaluate_arms`
``results/f6_padding_ablation_verification.csv``  every invariant and all 12 section 8 gates
``results/f6_padding_ablation_ledger.csv``        section 6.1's mechanism instrumentation
``results/f6_padding_ablation_seedvar.csv``       section 6.2's Q2 spread, no CV
``results/f6_padding_ablation_ceiling.csv``       section 2b/6.2's border-annotation geometry

This module writes **only** those five paths. `results/f5_nms_radius_ablation*.csv` and
`results/tm_ccoeff_threshold_axis_sweep*.csv` are read-only inputs to the reproduction
gates and must never be written; `main` asserts the output paths mechanically.

Run
---
    python f6_padding_ablation.py --spacing   # section 2b's ceiling geometry alone, ~15 s
    python f6_padding_ablation.py --gates     # pre-run gates, writes nothing, ~60 s
    python f6_padding_ablation.py --smoke     # one ROI, one seed, ~1 min
    python f6_padding_ablation.py             # the full run, ~45-60 min

Analysis lives in `f6_padding_ablation.ipynb`, which reads these CSVs and does no matching
-- the same split as `f5_nms_radius_ablation.py` / `.ipynb`, for the same reason.
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

IMAGES_DIR = "images/extra_valid"     # section 4: 2 ROIs per tumour type, all 7 domains
N_SEEDS = 5                           # section 4; 70 paired cells in 14 ROI clusters

CHANNEL = "hematoxylin_od"
METHOD = cv2.TM_CCOEFF
PEAK_MIN_DISTANCE = 7
SELF_HIT_RADIUS = 5.0

# section 3: the NMS radius and the scoring radius are the SAME 7.5 um in every arm. F5's
# decoupled radius is deliberately not used here -- one factor moves, and it is padding.
NMS_RADIUS_UM = ev.MIDOG_RADIUS_UM    # 7.5
MATCH_RADIUS_UM = ev.MIDOG_RADIUS_UM  # 7.5

# section 3. `pad_off` first so the ledger's reference arm is built before the arms that
# difference against it; `pad_on` before `pad_on_dropband` because the latter is a mask on
# the former's frame. The emitted order follows section 3's table.
ARM_TAGS = ("pad_off", "pad_on_dropband", "pad_on")
CONTROL_TAG = "pad_off"

# section 5. z is a nuisance axis; the headline is z = 1.0, and the grid is F5's so cells
# line up row for row with `results/f5_nms_radius_ablation.csv`.
Z_LEVELS = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
HEADLINE_Z = 1.0

# section 5: pinned as a CONSTANT, deliberately decoupled from min(Z_LEVELS). This is what
# makes section 8.9's two reproduction gates possible -- v1 and F5 both extracted at
# med - 1.5*mad, and deriving the floor from the z grid instead would silently move the
# deep pool and break both gates for a reason that has nothing to do with padding.
DEEP_FLOOR_Z = -1.5
MAX_PEAKS = 2_000_000

AXES = ("tm_score", "chromatin_od")   # section 5: co-primary, neither secondary
AXIS_RANK_KEY = {"tm_score": "score", "chromatin_od": "od"}

# section 6.1: a strict superset of `compare.BUDGETS`, so rows line up with v1 and F5.
BUDGETS_F6 = (10, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500,
              750, 1000, 1500, 2000, 3000, 5000, 10000)
HEADLINE_K = 250

# section 6.1: band detections in the top-K are counted at these three budgets, and P4's
# `margin_tp_density` is the TP count in ranks K-25 ... K-1 of the *`pad_off`* ranking at
# the same three. MARGIN_WINDOW is 25 because that is the window section 7's pre-measured
# table used (ranks 225-249 at K = 250).
BAND_K = (50, 250, 1000)
MARGIN_WINDOW = 25

# section 6.2: the headline spread is at K = 250; the rest are "the same spread at
# K in {50, 100, 500, 1000}", carried for shape.
SEEDVAR_K = (50, 100, 250, 500, 1000)

OD_PAD = tm.BASE_SIZE // 2            # 25 -- the chromatin window's OWN margin (v2 cell 14)

CFG = fs.FSConfig(channel=CHANNEL, base_size=tm.BASE_SIZE, scales=(1.0,), n_angles=1,
                  flips=(False,), peak_min_distance=PEAK_MIN_DISTANCE,
                  self_hit_radius=SELF_HIT_RADIUS)
BORDER = CFG.patch_size // 2          # 36 -- governs which annotation may be the SEED

# The band half-width. Derived per cell from the template bank as `(t-1)//2` and asserted
# equal to this, which is the constant the standalone ceiling table (which builds no
# templates) has to use. Section 8.1's `valid.sum()` identity is provable from PAD = (t-1)//2
# only because the bank is single-size, so `build_templates` asserts that too.
PAD_CONST = tm.BASE_SIZE // 2         # 25

OUT = Path("results")
OUT_MAIN = OUT / "f6_padding_ablation.csv"
OUT_VERIF = OUT / "f6_padding_ablation_verification.csv"
OUT_LEDGER = OUT / "f6_padding_ablation_ledger.csv"
OUT_SEEDVAR = OUT / "f6_padding_ablation_seedvar.csv"
OUT_CEILING = OUT / "f6_padding_ablation_ceiling.csv"
OUT_PATHS = (OUT_MAIN, OUT_VERIF, OUT_LEDGER, OUT_SEEDVAR, OUT_CEILING)

# Read-only gate inputs. Nothing in this module may write to either.
V1_CSV = "results/tm_ccoeff_threshold_axis_sweep.csv"          # section 8.9 gate A
F5_CSV = "results/f5_nms_radius_ablation.csv"                  # section 8.9 gate B

INTERIOR_IDENTITY_ATOL = 1e-4          # section 8.7; measured max difference 3.8e-6


def roi_files(images_dir: str = IMAGES_DIR):
    return sorted(f for f in os.listdir(images_dir) if f.endswith(".tiff"))


# ---------------------------------------------------------------------------------------
# Geometry: the border band, and distances to the ROI edge / the interior rectangle
# ---------------------------------------------------------------------------------------

def band_mask(cx, cy, H, W, pad):
    """Section 3's border band: exactly the region the unpadded search cannot place a peak in.

    ``cx < PAD or cx >= W-PAD or cy < PAD or cy >= H-PAD``. `fused_response` writes a
    ``t``-sized template's map at offset ``((t-1)//2, (t-1)//2)`` and it is ``H-t+1`` rows
    tall, so on the unpadded image the reachable rows are ``PAD ... H-PAD-1`` inclusive --
    the complement of this mask, exactly.
    """
    cx = np.asarray(cx, dtype=np.float64)
    cy = np.asarray(cy, dtype=np.float64)
    return (cx < pad) | (cx >= W - pad) | (cy < pad) | (cy >= H - pad)


def dist_to_interior(cx, cy, H, W, pad):
    """Euclidean distance from a point to the nearest pixel of the interior rectangle.

    The interior rectangle is ``[pad, W-1-pad] x [pad, H-1-pad]`` -- the pixels the
    *unpadded* search can place a detection on. Section 2b's ceiling is this distance
    against the ROI's match radius: a detection on the first valid row can still claim an
    annotation up to one match radius outside it, so a band annotation is only genuinely
    unreachable when this exceeds the radius. Reproduces section 2b's seven values exactly.
    """
    cx = np.asarray(cx, dtype=np.float64)
    cy = np.asarray(cy, dtype=np.float64)
    dx = np.maximum(0.0, np.maximum(pad - cx, cx - (W - 1 - pad)))
    dy = np.maximum(0.0, np.maximum(pad - cy, cy - (H - 1 - pad)))
    return np.hypot(dx, dy)


def dist_to_edge(cx, cy, H, W):
    """Distance from a point to the nearest ROI edge -- section 6.1's attribution column.

    A gained or lost annotation close to zero here is a band story; one far from zero is an
    interior story, and the ledger has to be able to tell them apart.
    """
    cx = np.asarray(cx, dtype=np.float64)
    cy = np.asarray(cy, dtype=np.float64)
    return np.minimum(np.minimum(cx, cy), np.minimum(W - 1 - cx, H - 1 - cy))


# ---------------------------------------------------------------------------------------
# Section 2b / 6.2 -- the reachability ceiling table
# ---------------------------------------------------------------------------------------

def ceiling_table(images, annotations, images_dir: str = IMAGES_DIR,
                  n_seeds: int = N_SEEDS) -> pd.DataFrame:
    """Section 6.2 / 2b. Every evaluation mitosis in the border band, per (ROI, seed).

    One row per band mitosis, carrying its distance to the interior rectangle and whether
    that exceeds the ROI's match radius -- the two numbers section 2b's ceiling of **zero**
    is made of. Cells with no band mitosis emit one placeholder row (``ann_id = -1``) so
    that every (ROI, seed) appears and ``n_band_mitoses`` is recoverable for all 70 cells
    rather than only for the ones that happen to have a border annotation.

    Needs no response map, so it runs standalone in ~15 s under ``--spacing``. The ROI
    shape comes from the images table; `run_roi` asserts it equals the loaded array's.
    """
    rows = []
    for fn in roi_files(images_dir):
        im = images[images["file_name"] == fn].iloc[0]
        H, W = int(im["height"]), int(im["width"])
        mpp = ds.roi_mpp(f"{images_dir}/{fn}")
        match_radius = ev.radius_px(mpp, MATCH_RADIUS_UM)
        gt = ds.image_annotations(annotations, fn)
        gt_mit = gt[gt["category_id"] == ds.MITOTIC]
        seeds, _, _ = draw_seeds(gt_mit, BORDER, (H, W), n_seeds)
        for si, seed in seeds:
            gt_eval = gt[gt["ann_id"] != int(seed["ann_id"])]
            mit = gt_eval[gt_eval["category_id"] == ds.MITOTIC]
            cx, cy = mit["cx"].to_numpy(), mit["cy"].to_numpy()
            inb = band_mask(cx, cy, H, W, PAD_CONST)
            d_int = dist_to_interior(cx, cy, H, W, PAD_CONST)
            d_edge = dist_to_edge(cx, cy, H, W)
            base = {"file_name": fn, "image_id": int(im["image_id"]),
                    "tumor_type": im["tumor_type"], "seed_index": si,
                    "seed_ann_id": int(seed["ann_id"]), "roi_h": H, "roi_w": W,
                    "pad_px": PAD_CONST, "match_radius_px": round(match_radius, 3),
                    "n_gt_mitotic_eval": int(len(mit)),
                    "n_band_mitoses": int(inb.sum()),
                    "granularity_1_over_n": round(1.0 / len(mit), 6) if len(mit) else np.nan}
            if not inb.any():
                rows.append({**base, "ann_id": -1, "cx": np.nan, "cy": np.nan,
                             "dist_to_interior_px": np.nan, "dist_to_edge_px": np.nan,
                             "exceeds_match_radius": False})
                continue
            for i in np.nonzero(inb)[0]:
                rows.append({**base, "ann_id": int(mit["ann_id"].iloc[i]),
                             "cx": float(cx[i]), "cy": float(cy[i]),
                             "dist_to_interior_px": round(float(d_int[i]), 3),
                             "dist_to_edge_px": round(float(d_edge[i]), 3),
                             "exceeds_match_radius": bool(d_int[i] > match_radius)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------
# The per-(ROI, seed) work: two response maps, two candidate frames, one mask
# ---------------------------------------------------------------------------------------

def build_templates(hem, seed_xy):
    """The template bank and the padding half-width, built **once** per (ROI, seed).

    The template is a property of the seed, not of the arm (section 3), so building it
    inside the arm loop would be both wasteful and a place for the arms to diverge.

    ``PAD`` is derived from the bank rather than hardcoded so it stays correct if
    `BASE_SIZE` or the augmentation scales ever change -- but section 8.1's `valid.sum()`
    identity is only *provable* from ``PAD = (t-1)//2`` when every template shares a size,
    so that is asserted rather than assumed.
    """
    patch = tm.read_padded_patch(hem, seed_xy[0], seed_xy[1], CFG.patch_size)
    templates, _ = tm.build_augmentations(patch, CFG.base_size, CFG.scales,
                                          CFG.n_angles, CFG.flips)
    shapes = {t.shape[:2] for t in templates}
    assert len(shapes) == 1, f"multi-size template bank {shapes}: PAD is no longer provable"
    (th, tw), = shapes
    assert th == tw and th % 2 == 1, f"template {th}x{tw} is not odd-square"
    return templates, (th - 1) // 2


def padded_response(hem, templates, pad):
    """`fused_response` on the REPLICATE-padded image, cropped back to the ROI.

    With an odd template side ``t`` and ``pad = (t-1)//2``, `matchTemplate` on a
    ``pad``-padded array returns a map of shape ``(H + 2*pad - t + 1, ...) = (H, W)``, which
    `fused_response` offsets to ``(pad, pad)`` -- so the crop below is all-valid by
    construction. Section 8.1 asserts it per cell rather than trusting the arithmetic.
    """
    H, W = hem.shape[:2]
    hem_p = cv2.copyMakeBorder(hem, pad, pad, pad, pad, borderType=cv2.BORDER_REPLICATE)
    fused_p, _, valid_p = tm.fused_response(hem_p, templates, CFG.scale_normalize,
                                            method=METHOD)
    fused = np.ascontiguousarray(fused_p[pad:pad + H, pad:pad + W])
    valid = np.ascontiguousarray(valid_p[pad:pad + H, pad:pad + W])
    del hem_p, fused_p, valid_p
    return fused, valid


def build_frame(hem_od, fused, valid, seed_xy, radius):
    """One response map's deep pool, `od`, and its post-NMS candidate frame.

    `od` is scored **once** on the deep pool and subset by the NMS keep mask: chromatin
    density depends only on (cx, cy) and the image, so it is identical for any point that
    survives, and scoring it after NMS would recompute identical numbers. The equivalence
    is asserted, not assumed (`od_hoist_gate`, run before every full run).

    ``hem_od`` is the `OD_PAD`-padded hematoxylin channel, hoisted to the ROI level and
    used in **all three arms** (section 3, last paragraph). Without it every candidate the
    border padding newly admits scores NaN and `compare._rank` sorts it last -- v2's own
    bug, which overstated its cost table by ~116x.

    Order is NMS first, self-hit second, pinned to v1/v2's so both reproduction gates hold.

    Returns ``(centers, scores, od, pool, kept_centers, meta)``. Nothing returned aliases a
    mutable array the caller then rebinds; `pool` is a fresh frame.
    """
    med, mad = tm.robust_stats(fused, valid)
    deep_floor = med + DEEP_FLOOR_Z * mad
    centers, scores = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE, deep_floor, MAX_PEAKS)

    shifted = pd.DataFrame({"cx": centers[:, 0] + OD_PAD, "cy": centers[:, 1] + OD_PAD})
    od = cm.score_detections(shifted, hem_od)["od"].to_numpy()

    keep = nms_by_distance(centers, scores, radius)
    c_keep, s_keep, od_keep = centers[keep], scores[keep], od[keep]
    ok = np.hypot(c_keep[:, 0] - seed_xy[0], c_keep[:, 1] - seed_xy[1]) > CFG.self_hit_radius
    pool = pd.DataFrame({"cx": c_keep[ok, 0], "cy": c_keep[ok, 1],
                         "score": s_keep[ok], "od": od_keep[ok]}).reset_index(drop=True)

    meta = {"valid_coverage": float(valid.mean()), "valid_sum": int(valid.sum()),
            "valid_all": bool(valid.all()), "map_median": float(med), "mad_scale": float(mad),
            "deep_floor": float(deep_floor), "n_peaks_deep": int(len(centers)),
            "od_nan": int(np.isnan(od).sum()), "n_pool": int(len(pool))}
    return centers, scores, od, pool, c_keep, meta


def coord_set(cx, cy):
    """A hashable identity for a peak. Peaks sit on integer pixels, so this is exact."""
    return set(zip(np.asarray(cx).tolist(), np.asarray(cy).tolist()))


# ---------------------------------------------------------------------------------------
# Section 6.1 -- the mechanism ledger
# ---------------------------------------------------------------------------------------

def ledger_row(pool_h, gt_eval, match_radius, axis, H, W, pad):
    """One arm's found-annotation set, band occupancy of the top-K, and duplicate accounting.

    ``pool_h`` is that arm's candidate frame already filtered to the headline z. The frame
    is sorted into its final ranking here with exactly `compare._rank`'s settings
    (``ascending=False, na_position="last", kind="mergesort"``) so the ledger describes the
    same ordering the main CSV scores.

    Three quantities come out of it:

    * ``found`` -- the mitotic ``ann_id`` set, which the caller differences against
      `pad_off`'s to get section 6.1's gained/lost.
    * ``band_top[K]`` -- **band detections** in the top K. Section 6.1 keeps band
      *detections* and band *annotations* strictly apart: the band is 25 px and the match
      radius ~30 px, so a band detection at ``cx = 20`` can claim an *interior* annotation
      at ``cx = 45``, and a detection at ``cx = 26`` -- not in the band -- can claim a band
      annotation at ``cx = 2``. These are detections.
    * ``margin[K]`` -- P4's TP density in ranks ``K-25 ... K-1`` (0-based ``[K-25 : K]``).
      Only the `pad_off` value is meaningful, because P4's product is "band detections
      *entering* the top-K x the TP density at the margin **of the list they displace**";
      the caller writes `pad_off`'s value onto every arm's row so the pivot is a straight
      join. NaN when the list is shorter than K -- section 2b puts the minimum
      `n_detections` at z = 3.0 at 2,301, so this should never fire, and if it does it must
      be visible rather than silently a partial count.
    """
    ranked = pool_h.sort_values(AXIS_RANK_KEY[axis], ascending=False, na_position="last",
                                kind="mergesort").reset_index(drop=True)
    det, gt_out = ev.bucket_detections(ranked, gt_eval, match_radius)

    found = set(int(a) for a in
                gt_out[(gt_out["category_id"] == ds.MITOTIC) & gt_out["found"]]["ann_id"])

    inb = band_mask(ranked["cx"].to_numpy(), ranked["cy"].to_numpy(), H, W, pad)
    band_top = {k: int(inb[:k].sum()) for k in BAND_K}

    tp = (det["bucket"].to_numpy() == ev.HUMAN_CORRECT_LABEL)
    margin = {k: (int(tp[k - MARGIN_WINDOW:k].sum()) if len(tp) >= k else float("nan"))
              for k in BAND_K}

    # Duplicate FPs, reported WITH their denominator: a second detection inside an already
    # claimed object's match circle buckets as an unannotated false positive, and the raw
    # count rises mechanically with list length, so it is meaningless without n_within.
    if len(det) and len(gt_eval):
        d_gt, _ = KDTree(gt_eval[["cx", "cy"]].to_numpy()).query(
            det[["cx", "cy"]].to_numpy(), k=1)
        inside = d_gt[:, 0] <= match_radius
        dup = inside & (det["bucket"] == ev.NON_HUMAN_FINDINGS).to_numpy()
        n_within, n_dup = int(inside.sum()), int(dup.sum())
    else:
        n_within = n_dup = 0

    return {"found": found, "band_top": band_top, "margin": margin,
            "n_within_match_radius": n_within, "n_dup_fp": n_dup,
            "n_detections": int(len(ranked)),
            "n_band_candidates_h": int(inb.sum())}


def interior_peak_split(centers_off, centers_on, pool_off_h, pool_on_h, H, W, pad):
    """Section 6.1's "interior peaks lost, split by mechanism".

    **The plan under-specifies this** -- section 2c quotes 99-122 interior peaks lost on
    single cells, section 3 quotes 228-267 for a *different* object (an emulation), and
    section 6.1 quotes 267 pre-NMS / 119 post-NMS on 402.tiff/s0. Rather than pick one, all
    four set differences below are emitted and the exact algebra is written down here, so
    the analysis is not locked into this file's attribution. Let

        ``Ioff``  = interior coords of `pad_off`'s DEEP pool (post-extraction, pre-NMS)
        ``Ion``   = interior coords of `pad_on`'s  DEEP pool
        ``Foff``  = interior coords of `pad_off`'s pool at z = HEADLINE_Z (post-NMS, post-self-hit)
        ``Fon``   = interior coords of `pad_on`'s  pool at z = HEADLINE_Z

    then

        ``lost_prenms``      = |Ioff \\ Ion|                 -- the grey-dilation channel, raw
        ``gained_prenms``    = |Ion \\ Ioff|                 -- its counterpart, for context
        ``lost_total_h``     = |Foff \\ Fon|                 -- what actually reaches the ranking
        ``lost_dilation_h``  = |(Foff \\ Fon) \\ Ion|        -- of those, never extracted padded
        ``lost_postnms_z``   = |(Foff \\ Fon) & Ion|         -- extracted, then lost to NMS or the z cut

    and ``lost_dilation_h + lost_postnms_z == lost_total_h`` exactly, which is the split
    section 6.1 asks for ("identity of the suppressing peak is only defined for the NMS
    case"). Everything here is axis-free -- which peak suppressed which is entirely
    upstream of ranking -- so it is computed once per cell and replicated onto both axes'
    ledger rows.

    Caveat, stated because it bounds the interpretation: this is an **exact-coordinate**
    set difference. The two maps agree only to ~3.8e-6 (section 8.7), so a local maximum
    that shifts by one pixel between them counts as one loss and one gain rather than as a
    move. `gained_prenms` is emitted so that channel stays visible.
    """
    def interior(cx, cy):
        keep = ~band_mask(cx, cy, H, W, pad)
        return coord_set(np.asarray(cx)[keep], np.asarray(cy)[keep])

    i_off = interior(centers_off[:, 0], centers_off[:, 1])
    i_on = interior(centers_on[:, 0], centers_on[:, 1])
    f_off = interior(pool_off_h["cx"].to_numpy(), pool_off_h["cy"].to_numpy())
    f_on = interior(pool_on_h["cx"].to_numpy(), pool_on_h["cy"].to_numpy())

    lost_h = f_off - f_on
    return {"interior_lost_prenms": len(i_off - i_on),
            "interior_gained_prenms": len(i_on - i_off),
            "interior_lost_total_h": len(lost_h),
            "interior_lost_dilation_h": len(lost_h - i_on),
            "interior_lost_postnms_z": len(lost_h & i_on)}


# ---------------------------------------------------------------------------------------
# The main loop
# ---------------------------------------------------------------------------------------

def run_roi(fn, images, annotations, images_dir, n_seeds, checks, ledger_rows,
            pool_sizes, verify_shortcut=True, quiet=False):
    """One ROI: load once, then n_seeds x (2 response maps -> 3 arms) off one template bank."""
    im = images[images["file_name"] == fn].iloc[0]
    image_id, domain = int(im["image_id"]), im["tumor_type"]
    path = f"{images_dir}/{fn}"

    t_roi = time.time()
    rgb = ds.load_roi(path)
    mpp = ds.roi_mpp(path)
    hem = ch.to_channel(rgb, CHANNEL)          # hoisted: 7.5 s per ROI, not per seed
    roi_shape = rgb.shape
    H, W = int(roi_shape[0]), int(roi_shape[1])
    assert (H, W) == (int(im["height"]), int(im["width"])), (
        f"{fn}: images table says {(int(im['height']), int(im['width']))} but the array is "
        f"{(H, W)} -- the standalone ceiling table and the run disagree on the ROI shape")
    del rgb                                    # 117 MB back before the maps are allocated

    match_radius = ev.radius_px(mpp, MATCH_RADIUS_UM)
    nms_radius = ev.radius_px(mpp, NMS_RADIUS_UM)

    # The `od` channel, hoisted to the ROI: it depends on the image alone, and it is the
    # SAME padded channel in all three arms (section 3, last paragraph).
    hem_od = cv2.copyMakeBorder(hem, OD_PAD, OD_PAD, OD_PAD, OD_PAD,
                                borderType=cv2.BORDER_REPLICATE)

    gt = ds.image_annotations(annotations, fn)
    gt_mit = gt[gt["category_id"] == ds.MITOTIC]
    seeds, info, seed_records = draw_seeds(gt_mit, BORDER, roi_shape, n_seeds)
    pool_sizes[fn] = int(info.n_after_border)          # section 4: reported beside every Q2 number
    checks.append(inv.check_distinct_seeds(seed_records, pool_size=info.n_after_border,
                                           label=fn))  # section 8.6

    frames = []
    for si, seed in seeds:
        t0 = time.time()
        seed_xy = (float(seed["cx"]), float(seed["cy"]))
        gt_eval = gt[gt["ann_id"] != int(seed["ann_id"])].reset_index(drop=True)
        gt_mit_eval = gt_eval[gt_eval["category_id"] == ds.MITOTIC]

        # -- section 8.11: no evaluation annotation may sit in the self-hit hole
        d_seed = np.hypot(gt_eval["cx"] - seed_xy[0], gt_eval["cy"] - seed_xy[1])
        n_hole = int((d_seed <= nms_radius).sum())
        checks.append({"check": "gt_within_seed_hole", "label": f"{fn}/si{si}",
                       "n_items": n_hole, "passed": bool(n_hole == 0)})
        assert n_hole == 0, (
            f"{fn}/si{si}: {n_hole} evaluation annotations sit inside the self-hit hole")

        templates, pad = build_templates(hem, seed_xy)
        assert pad == PAD_CONST, f"{fn}/si{si}: derived PAD {pad} != PAD_CONST {PAD_CONST}"

        keep_maps = bool(si == 0)     # section 8.7 and 8.10 run at si == 0 only

        # ---------------- arm 1: pad_off -- the unpadded response ----------------
        # NOT emulated from the padded map. Section 3: gate A compares `n_detections` on
        # exact integer equality and a `_FLOOR`-written emulation still loses 1 peak of
        # 20,855 on 459.tiff, so an emulated pad_off would not be v1.
        fused_off, _, valid_off = tm.fused_response(hem, templates, CFG.scale_normalize,
                                                    method=METHOD)
        centers_off, scores_off, od_off, pool_off, kept_off, meta_off = build_frame(
            hem_od, fused_off, valid_off, seed_xy, nms_radius)

        # section 8.1: provable from PAD = (t-1)//2 and a single-size bank, not empirical
        want_valid = (H - 2 * pad) * (W - 2 * pad)
        checks.append({"check": "valid_geometry", "label": f"{fn}/si{si}/pad_off",
                       "n_items": meta_off["valid_sum"], "n": want_valid,
                       "passed": bool(meta_off["valid_sum"] == want_valid)})
        assert meta_off["valid_sum"] == want_valid, (
            f"{fn}/si{si}: pad_off valid.sum()={meta_off['valid_sum']} != {want_valid}")

        if not keep_maps:
            del fused_off, valid_off          # section 11: one 39 MP map at a time

        # ---------------- arms 2/3: pad_on and its dropband mask ----------------
        fused_on, valid_on = padded_response(hem, templates, pad)
        centers_on, scores_on, od_on, pool_on, kept_on, meta_on = build_frame(
            hem_od, fused_on, valid_on, seed_xy, nms_radius)

        checks.append({"check": "valid_geometry", "label": f"{fn}/si{si}/pad_on",
                       "n_items": meta_on["valid_sum"], "n": H * W,
                       "passed": bool(meta_on["valid_all"])})
        assert meta_on["valid_all"], f"{fn}/si{si}: padding did not make the whole ROI valid"

        # `pad_on_dropband` is ONE BOOLEAN MASK on pad_on's candidate frame -- no second
        # extraction, no second NMS, no second `od` scoring (section 3). `pool_on` is never
        # rebound; the mask produces a new frame.
        band_on = band_mask(pool_on["cx"].to_numpy(), pool_on["cy"].to_numpy(), H, W, pad)
        pool_drop = pool_on.loc[~band_on].reset_index(drop=True)
        kept_drop = kept_on[~band_mask(kept_on[:, 0], kept_on[:, 1], H, W, pad)]

        # section 8.12: strict subset, asserted on the CANDIDATE SETS -- ranks shift by
        # construction when rows are removed, so a rank-level assertion would be false.
        s_on = coord_set(pool_on["cx"].to_numpy(), pool_on["cy"].to_numpy())
        s_drop = coord_set(pool_drop["cx"].to_numpy(), pool_drop["cy"].to_numpy())
        subset_ok = s_drop < s_on if len(pool_drop) < len(pool_on) else s_drop == s_on
        checks.append({"check": "dropband_subset", "label": f"{fn}/si{si}",
                       "n": len(s_on), "n_items": len(s_drop),
                       "n_dropped": int(band_on.sum()),
                       "passed": bool(s_drop <= s_on and subset_ok)})
        assert s_drop <= s_on, f"{fn}/si{si}: dropband is not a subset of pad_on"

        # section 8.1 corollary, asserted because the whole decomposition rests on it:
        # pad_off cannot place a candidate in the band, by construction.
        n_band_off = int(band_mask(pool_off["cx"].to_numpy(), pool_off["cy"].to_numpy(),
                                   H, W, pad).sum())
        assert n_band_off == 0, f"{fn}/si{si}: pad_off has {n_band_off} band candidates"

        # ---------------- section 8.7: interior identity, once per ROI ----------------
        if keep_maps:
            a = fused_on[pad:H - pad, pad:W - pad]
            b = fused_off[pad:H - pad, pad:W - pad]
            d = float(np.max(np.abs(a.astype(np.float64) - b.astype(np.float64))))
            ok = bool(np.allclose(a, b, atol=INTERIOR_IDENTITY_ATOL, rtol=0.0))
            checks.append({"check": "interior_identity", "label": f"{fn}/si{si}",
                           "n": int(a.size), "max_abs_diff": d, "atol": INTERIOR_IDENTITY_ATOL,
                           "passed": ok})
            assert ok, f"{fn}/si{si}: interior identity failed, max |diff| = {d:.3e}"

        # ---------------- section 8.2/8.4/8.5: per-arm invariants ----------------
        pools = {"pad_off": pool_off, "pad_on_dropband": pool_drop, "pad_on": pool_on}
        metas = {"pad_off": meta_off, "pad_on_dropband": meta_on, "pad_on": meta_on}
        kepts = {"pad_off": kept_off, "pad_on_dropband": kept_drop, "pad_on": kept_on}
        # `pad_on_dropband` shares pad_on's map, so its od_nan and deep-pool cap records are
        # pad_on's by construction; both are still emitted per arm so section 8's "all three
        # arms" is literally satisfied in the CSV.
        for tag in ARM_TAGS:
            m = metas[tag]
            checks.append({"check": "od_nan", "label": f"{fn}/si{si}/{tag}",
                           "n_items": m["od_nan"], "passed": bool(m["od_nan"] == 0)})
            assert m["od_nan"] == 0, f"{fn}/si{si}/{tag}: {m['od_nan']} candidates still NaN"
            # section 8.5: on the PRE-NMS pool. `compare.evaluate_arms` checks the post-NMS,
            # post-z length, which against a 2,000,000 cap cannot fail, so on its own it
            # records a pass it could not have withheld.
            checks.append(inv.check_no_cap(m["n_peaks_deep"], (MAX_PEAKS,),
                                           label=f"{fn}/si{si}/{tag}/deep_pool"))
            # section 8.4
            checks.append(inv.check_min_separation(kepts[tag], nms_radius,
                                                   label=f"{fn}/si{si}/{tag}"))

        # ---------------- section 6.1: the ledger, at the headline z ----------------
        cuts = {tag: metas[tag]["map_median"] + HEADLINE_Z * metas[tag]["mad_scale"]
                for tag in ARM_TAGS}
        pools_h = {tag: pools[tag][pools[tag]["score"] >= cuts[tag]] for tag in ARM_TAGS}

        split = interior_peak_split(centers_off, centers_on,
                                    pools_h["pad_off"], pools_h["pad_on"], H, W, pad)

        ann_xy = {int(r.ann_id): (float(r.cx), float(r.cy)) for r in gt_eval.itertuples()}

        for axis in AXES:
            leds = {tag: ledger_row(pools_h[tag], gt_eval, match_radius, axis, H, W, pad)
                    for tag in ARM_TAGS}
            margin = leds[CONTROL_TAG]["margin"]      # P4: the displaced list is pad_off's
            ctrl_found = leds[CONTROL_TAG]["found"]
            for tag in ARM_TAGS:
                led = leds[tag]
                gained = sorted(led["found"] - ctrl_found)
                lost = sorted(ctrl_found - led["found"])

                def _edges(ids):
                    return [round(float(dist_to_edge(*ann_xy[a], H, W)), 2) for a in ids]

                def _in_band(ids):
                    return int(sum(bool(band_mask(*ann_xy[a], H, W, pad)) for a in ids))

                row = {"file_name": fn, "image_id": image_id, "tumor_type": domain,
                       "seed_index": si, "seed_ann_id": int(seed["ann_id"]),
                       "arm_tag": tag, "arm": f"{axis}@{tag}", "axis": axis, "z": HEADLINE_Z,
                       "roi_h": H, "roi_w": W, "pad_px": pad,
                       "match_radius_px": round(match_radius, 3),
                       "map_median": round(metas[tag]["map_median"], 6),
                       "mad_scale": round(metas[tag]["mad_scale"], 6),
                       "deep_floor": round(metas[tag]["deep_floor"], 6),
                       "z_cut": round(cuts[tag], 6),
                       "n_peaks_deep": metas[tag]["n_peaks_deep"],
                       "n_pool": len(pools[tag]),
                       "n_detections": led["n_detections"],
                       "n_gt_mitotic": int(len(gt_mit_eval)),
                       "n_found": len(led["found"]),
                       "n_gained_vs_pad_off": len(gained), "n_lost_vs_pad_off": len(lost),
                       "n_gained_in_band": _in_band(gained),
                       "n_lost_in_band": _in_band(lost),
                       "gained_ann_ids": ";".join(str(a) for a in gained),
                       "gained_edge_dists_px": ";".join(str(d) for d in _edges(gained)),
                       "lost_ann_ids": ";".join(str(a) for a in lost),
                       "lost_edge_dists_px": ";".join(str(d) for d in _edges(lost)),
                       # section 6.1: band candidates in the FULL pool, and band DETECTIONS
                       # in the top-K -- detections, never annotations.
                       "n_band_candidates_pool": int(
                           band_mask(pools[tag]["cx"].to_numpy(), pools[tag]["cy"].to_numpy(),
                                     H, W, pad).sum()),
                       "n_band_candidates_h": led["n_band_candidates_h"],
                       "n_within_match_radius": led["n_within_match_radius"],
                       "n_dup_fp": led["n_dup_fp"],
                       "dup_fp_frac": round(led["n_dup_fp"] / max(led["n_detections"], 1), 6)}
                for k in BAND_K:
                    row[f"band_in_top{k}"] = led["band_top"][k]
                    # pad_off's margin density, written onto every arm's row so P4's product
                    # (band_in_topK x margin_tp_density_k) is a straight per-row multiply.
                    row[f"margin_tp_density_k{k}"] = margin[k]
                row.update(split)                     # axis-free, replicated onto both axes
                ledger_rows.append(row)

        # ---------------- section 8.10: the one-match-many-z shortcut ----------------
        if verify_shortcut and keep_maps:
            for tag, fu, va, pl, mt in (("pad_off", fused_off, valid_off, pool_off, meta_off),
                                        ("pad_on", fused_on, valid_on, pool_on, meta_on)):
                for z in (HEADLINE_Z, max(Z_LEVELS)):
                    cut = mt["map_median"] + z * mt["mad_scale"]
                    c2, s2 = tm.extract_peaks(fu, va, PEAK_MIN_DISTANCE, cut, MAX_PEAKS)
                    k2 = nms_by_distance(c2, s2, nms_radius)
                    c2k = c2[k2]
                    ok2 = np.hypot(c2k[:, 0] - seed_xy[0],
                                   c2k[:, 1] - seed_xy[1]) > CFG.self_hit_radius
                    direct = c2k[ok2]
                    filt = pl.loc[pl["score"] >= cut, ["cx", "cy"]].to_numpy()
                    a2 = direct[np.lexsort((direct[:, 1], direct[:, 0]))] if len(direct) else direct
                    b2 = filt[np.lexsort((filt[:, 1], filt[:, 0]))] if len(filt) else filt
                    passed = (a2.shape == b2.shape) and bool(np.allclose(a2, b2))
                    checks.append({"check": "shortcut", "label": f"{fn}/si{si}/{tag}/z{z}",
                                   "n_direct": len(direct), "n_filtered": len(filt),
                                   "passed": passed})
                    assert passed, f"{fn}/si{si}/{tag}: shortcut failed at z={z}"
        if keep_maps:
            del fused_off, valid_off        # freed early at si > 0, held here for 8.7/8.10
        del fused_on, valid_on

        # ---------------- the arms ----------------
        ctx_base = {"file_name": fn, "image_id": image_id, "tumor_type": domain,
                    "seed_index": si, "seed_ann_id": int(seed["ann_id"]),
                    "mpp": round(mpp, 5), "pad_px": pad,
                    "match_radius_px": round(match_radius, 3),
                    "nms_radius_px": round(nms_radius, 3),
                    "n_gt_within_seed_hole": n_hole,
                    "seed_pool_size": int(info.n_after_border)}
        cov_cache = {}                    # section 5: keyed (file, seed, arm_tag, z)
        for tag in ARM_TAGS:
            m, pool = metas[tag], pools[tag]
            ctx = dict(ctx_base, arm_tag=tag,
                       map_median=round(m["map_median"], 6),
                       mad_scale=round(m["mad_scale"], 6),
                       deep_floor=round(m["deep_floor"], 6),
                       n_peaks_deep=m["n_peaks_deep"], n_pool=len(pool),
                       n_band_candidates_pool=int(band_mask(
                           pool["cx"].to_numpy(), pool["cy"].to_numpy(), H, W, pad).sum()))
            arms = []
            for z in Z_LEVELS:
                cut = m["map_median"] + z * m["mad_scale"]
                sub = pool[pool["score"] >= cut]
                limited = bool(cut < m["deep_floor"])
                for axis in AXES:
                    arms.append(cp.Arm(
                        f"{axis}@{tag}", (lambda d=sub: d), rank_key=AXIS_RANK_KEY[axis],
                        seeded=True, z=z, z_dependent=True, floor_limited=limited,
                        # section 8.3: the radius is 7.5 um in EVERY arm, so every arm is a
                        # positive control for `check_nms_radius` -- F6 does not decouple it.
                        nms_radius=nms_radius, caps=(MAX_PEAKS,),
                        coverage_key=f"{fn}|{si}|{tag}|{z}",   # a set, not a sort order
                        extra={"axis": axis, "z_cut": round(cut, 6), "arm_tag": tag}))
            frames.append(cp.evaluate_arms(arms, gt_eval, match_radius, roi_shape=roi_shape,
                                           mpp=mpp, budgets=BUDGETS_F6, context=ctx,
                                           coverage_cache=cov_cache, checks=checks))

        if not quiet:
            print(f"  [{fn} si{si}] deep off/on={meta_off['n_peaks_deep']:6d}/"
                  f"{meta_on['n_peaks_deep']:6d} pools="
                  f"{'/'.join(str(len(pools[t])) for t in ARM_TAGS)} "
                  f"band={int(band_on.sum())} [{time.time() - t0:.0f}s]", flush=True)
        del centers_off, scores_off, od_off, centers_on, scores_on, od_on
        del pool_off, pool_on, pool_drop, pools, pools_h

    del hem, hem_od
    if not quiet:
        print(f"[{fn}] {domain:34s} done [{time.time() - t_roi:.0f}s]", flush=True)
    return frames


# ---------------------------------------------------------------------------------------
# Section 8.8 / 8.9 -- the post-hoc gates
# ---------------------------------------------------------------------------------------

def join_integrity_gate(results: pd.DataFrame) -> pd.DataFrame:
    """Section 8.8. Join integrity, **not** "additivity".

    ``(A-C) == (A-B) + (B-C)`` holds for any three reals and cannot fail except on NaN or a
    broken join; draft 1 presented it as evidence that the decomposition holds together. It
    is kept, demoted to what it is: every ``(file, seed, axis, z, budget)`` must carry all
    three arms exactly once, with no NaN in ``recall_at_budget``, and the dedup key of
    section 5 must be unique. The residual is reported so the reader can see it is 0.0 by
    arithmetic rather than by measurement.
    """
    rows = []
    key = ["file_name", "seed_index", "axis", "z", "budget"]
    n_dup = int(results.duplicated(key + ["arm_tag"]).sum())
    rows.append({"check": "dedup_key", "label": "file|seed|arm|z|budget",
                 "n_items": n_dup, "passed": bool(n_dup == 0)})

    cnt = results.groupby(key, dropna=False)["arm_tag"].nunique()
    n_bad = int((cnt != len(ARM_TAGS)).sum())
    rows.append({"check": "join_integrity", "label": "all three arms present",
                 "n": int(len(cnt)), "n_items": n_bad, "passed": bool(n_bad == 0)})

    n_nan = int(results["recall_at_budget"].isna().sum())
    rows.append({"check": "join_integrity", "label": "recall_at_budget not NaN",
                 "n_items": n_nan, "passed": bool(n_nan == 0)})

    wide = results.pivot_table(index=key, columns="arm_tag", values="recall_at_budget")
    if set(ARM_TAGS).issubset(wide.columns) and len(wide):
        resid = ((wide["pad_on"] - wide["pad_off"])
                 - ((wide["pad_on"] - wide["pad_on_dropband"])
                    + (wide["pad_on_dropband"] - wide["pad_off"]))).abs()
        rows.append({"check": "join_integrity", "label": "(A-C) == (A-B)+(B-C)",
                     "n": int(len(resid)), "max_abs_diff": float(np.nanmax(resid)),
                     "n_items": int(resid.isna().sum()),
                     "passed": bool(np.nanmax(resid) == 0.0 and not resid.isna().any())})
    else:
        rows.append({"check": "join_integrity", "label": "(A-C) == (A-B)+(B-C)",
                     "passed": False, "note": "not all three arms in the pivot"})
    return pd.DataFrame(rows)


def _repro(results, ref_csv, ref_arm_of, my_tag, check_name, ref_filter=None):
    """Shared body of section 8.9's two reproduction gates.

    **Scope, carried over from F5's implementation and stated so the pass is not read as
    wider than it is:** the gate compares ``recall_at_budget`` (allclose) and
    ``n_detections`` (exact integer equality). ``full_list_recall`` is merged but never
    compared. That is still a whole-pipeline check -- a silent divergence anywhere from the
    seed draw to the bucketing would move one of the two.

    Mapping: v1's ``arm`` column is the bare axis, F5's and F6's are ``{axis}@{tag}``.
    """
    cols = ["file_name", "z", "budget", "recall_at_budget", "n_detections", "full_list_recall"]
    if not Path(ref_csv).exists():
        return pd.DataFrame([{"check": check_name, "label": ref_csv, "passed": False,
                              "note": f"{ref_csv} absent -- gate could not be evaluated"}])
    ref = pd.read_csv(ref_csv)
    if ref_filter is not None:
        ref = ref_filter(ref)
    mine = results[(results["seed_index"] == 0) & (results["arm_tag"] == my_tag)]
    rows = []
    for axis in AXES:
        a = ref[ref["arm"] == ref_arm_of(axis)][cols].drop_duplicates(
            ["file_name", "z", "budget"])
        b = mine[mine["arm"] == f"{axis}@{my_tag}"][cols].drop_duplicates(
            ["file_name", "z", "budget"])
        m = a.merge(b, on=["file_name", "z", "budget"], suffixes=("_ref", "_f6"))
        if not len(m):
            rows.append({"check": check_name, "label": axis, "n_cells": 0, "passed": False,
                         "note": "no overlapping cells"})
            continue
        same_n = bool((m["n_detections_ref"] == m["n_detections_f6"]).all())
        same_r = bool(np.allclose(m["recall_at_budget_ref"], m["recall_at_budget_f6"],
                                  atol=1e-9, equal_nan=True))
        rows.append({"check": check_name, "label": axis, "n_cells": int(len(m)),
                     "n_rois": int(m["file_name"].nunique()),
                     "max_abs_recall_diff": float(np.nanmax(
                         np.abs(m["recall_at_budget_ref"] - m["recall_at_budget_f6"]))),
                     "max_abs_ndet_diff": int(
                         (m["n_detections_ref"] - m["n_detections_f6"]).abs().max()),
                     "passed": bool(same_n and same_r)})
    return pd.DataFrame(rows)


def reproduction_gate_a(results: pd.DataFrame) -> pd.DataFrame:
    """Section 8.9 A. `pad_off` at si == 0 reproduces v1 on the 7 overlapping ROIs.

    v1 (`tm_ccoeff_threshold_axis_sweep.csv`) is the unpadded, 7.5 um-NMS, SEED_INDEX = 0
    run. `draw_seeds(s=0)` was verified against its committed ``seed_ann_id`` on all 7 ROIs
    before this run, and `DEEP_FLOOR_Z` is pinned at v1's -1.5, so this must hold at every
    shared (z, budget) -- v1's z grid is a superset of F6's and its budget grid a subset.
    """
    return _repro(results, V1_CSV, lambda ax: ax, "pad_off", "v1_reproduction_A")


def reproduction_gate_b(results: pd.DataFrame) -> pd.DataFrame:
    """Section 8.9 B. `pad_on` at si == 0 reproduces F5's `r7.5` arm on all 14 ROIs.

    F5's `r7.5` is padding ON at the 7.5 um NMS radius with `od` on the `OD_PAD`-padded
    channel -- F6's `pad_on` by construction. The two runs share the seed draw, the z grid,
    the budget grid and the deep floor, so this is a 14-ROI whole-pipeline equality.
    """
    return _repro(results, F5_CSV, lambda ax: f"{ax}@r7.5", "pad_on", "f5_r75_reproduction_B",
                  ref_filter=lambda d: d[d["seed_index"] == 0])


# ---------------------------------------------------------------------------------------
# Section 6.2 -- Q2, the click-noise table
# ---------------------------------------------------------------------------------------

def seedvar_table(results: pd.DataFrame, pool_sizes: dict) -> pd.DataFrame:
    """Section 6.2. Per (ROI, arm, axis, K): the seed values' mean, SD, min, max, range.

    **No CV.** It is meaningless for a bounded metric near its ceiling -- 013.tiff has mean
    1.0 and SD 0 at K = 250, giving CV 0, which does not mean "less variable than" an ROI at
    0.30. SD is reported with its mean instead.

    ``seed_pool_size`` (7-184 across these ROIs) rides on every row because 201.tiff's 5
    seeds are 71 % of its eligible population and 301.tiff's are 3 %; those estimates are
    not on the same sampling footing and the table must not pretend otherwise.

    ``granularity_1_over_n`` is one TP expressed as recall -- 0.0042 (238 mitoses) to 0.059
    (17 mitoses) across the 14 ROIs. An SD at that scale is quantisation, not click
    sensitivity, and a 14x difference in the smallest expressible SD is exactly the kind of
    artefact a bare SD column invites.

    ``sd`` is the sample SD (``ddof=1``), the convention `tm_variant_sweep.seed_variance`
    already uses. Zero-variance ROIs report 0.0 and are never used as a denominator here;
    section 6.3's practical ruler is where the ``n/a`` rule applies.
    """
    sub = results[(results["z"] == HEADLINE_Z) & (results["budget"].isin(SEEDVAR_K))]
    rows = []
    keys = ["file_name", "tumor_type", "arm_tag", "axis", "budget"]
    for (fn, dom, tag, axis, k), g in sub.groupby(keys, dropna=False):
        v = pd.to_numeric(g["recall_at_budget"], errors="coerce").dropna()
        n_eval = int(g["n_gt_mitotic"].iloc[0])
        rows.append({
            "file_name": fn, "tumor_type": dom, "arm_tag": tag,
            "arm": f"{axis}@{tag}", "axis": axis, "budget": int(k),
            "n_seeds": int(len(v)),
            "mean": float(v.mean()) if len(v) else np.nan,
            "sd": float(v.std(ddof=1)) if len(v) > 1 else np.nan,
            "min": float(v.min()) if len(v) else np.nan,
            "max": float(v.max()) if len(v) else np.nan,
            "range": float(v.max() - v.min()) if len(v) else np.nan,
            "seed_pool_size": int(pool_sizes.get(fn, -1)),
            "n_gt_mitotic": n_eval,
            "granularity_1_over_n": (1.0 / n_eval) if n_eval else np.nan,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------
# Pre-run gates (--gates): nothing here writes a file
# ---------------------------------------------------------------------------------------

def od_hoist_gate(images, annotations, images_dir=IMAGES_DIR) -> dict:
    """`od` scored once on the deep pool and subset must equal `od` scored after NMS.

    F4 section 9 gate 7, borrowed with the optimisation it licenses. If this fails the
    hoist in `build_frame` is invalid and every `chromatin_od`-ranked number in the run is
    wrong. Run on the padded map, which is the one with band candidates in it.
    """
    fn = roi_files(images_dir)[0]
    rgb = ds.load_roi(f"{images_dir}/{fn}")
    mpp = ds.roi_mpp(f"{images_dir}/{fn}")
    hem = ch.to_channel(rgb, CHANNEL)
    del rgb
    gt = ds.image_annotations(annotations, fn)
    gt_mit = gt[gt["category_id"] == ds.MITOTIC]
    seeds, _, _ = draw_seeds(gt_mit, BORDER, hem.shape, 1)
    _, seed = seeds[0]
    seed_xy = (float(seed["cx"]), float(seed["cy"]))
    templates, pad = build_templates(hem, seed_xy)
    hem_od = cv2.copyMakeBorder(hem, OD_PAD, OD_PAD, OD_PAD, OD_PAD,
                                borderType=cv2.BORDER_REPLICATE)
    fused, valid = padded_response(hem, templates, pad)
    radius = ev.radius_px(mpp, NMS_RADIUS_UM)
    centers, scores, od, _, _, _ = build_frame(hem_od, fused, valid, seed_xy, radius)

    keep = nms_by_distance(centers, scores, radius)
    hoisted = od[keep]
    after = cm.score_detections(
        pd.DataFrame({"cx": centers[keep, 0] + OD_PAD, "cy": centers[keep, 1] + OD_PAD}),
        hem_od)["od"].to_numpy()
    passed = bool(np.allclose(hoisted, after, equal_nan=True))
    return {"check": "od_hoist_equivalence", "label": fn, "n": int(len(keep)),
            "max_abs_diff": float(np.nanmax(np.abs(hoisted - after))) if len(keep) else 0.0,
            "passed": passed}


def interior_identity_gate(images, annotations, images_dir=IMAGES_DIR) -> dict:
    """Section 8.7, run standalone on one ROI before the full run commits 45 minutes.

    The padded map cropped to the interior must equal the unpadded map on the interior to
    ``atol = 1e-4``. Measured max difference 3.8e-6 (OpenCV's DFT block tiling differs with
    image size), so the two arms are identical to ~1e-6 -- **not** "byte-identical", which
    draft 1 claimed and which would have made gate A fail.
    """
    fn = roi_files(images_dir)[0]
    rgb = ds.load_roi(f"{images_dir}/{fn}")
    hem = ch.to_channel(rgb, CHANNEL)
    H, W = hem.shape[:2]
    del rgb
    gt = ds.image_annotations(annotations, fn)
    gt_mit = gt[gt["category_id"] == ds.MITOTIC]
    seeds, _, _ = draw_seeds(gt_mit, BORDER, hem.shape, 1)
    seed_xy = (float(seeds[0][1]["cx"]), float(seeds[0][1]["cy"]))
    templates, pad = build_templates(hem, seed_xy)
    fused_on, valid_on = padded_response(hem, templates, pad)
    fused_off, _, valid_off = tm.fused_response(hem, templates, CFG.scale_normalize,
                                                method=METHOD)
    a = fused_on[pad:H - pad, pad:W - pad].astype(np.float64)
    b = fused_off[pad:H - pad, pad:W - pad].astype(np.float64)
    d = float(np.max(np.abs(a - b)))
    want_valid = (H - 2 * pad) * (W - 2 * pad)
    return {"check": "interior_identity", "label": fn, "n": int(a.size),
            "max_abs_diff": d, "atol": INTERIOR_IDENTITY_ATOL,
            "n_items": int(valid_off.sum()),
            "passed": bool(d <= INTERIOR_IDENTITY_ATOL and valid_on.all()
                           and int(valid_off.sum()) == want_valid)}


def seed_identity_gate(images, annotations, images_dir=IMAGES_DIR) -> pd.DataFrame:
    """Pre-flight, not a section 8 gate: does `draw_seeds(si=0)` match both gate inputs?

    Both reproduction gates carry ``seed_ann_id``. If the draw has diverged, gates A and B
    fail for a reason that has nothing to do with padding -- and without this check that
    would only surface 45 minutes later. Costs nothing: the ROI shape comes from the images
    table, so no pixels are read.
    """
    rows = []
    refs = {}
    if Path(V1_CSV).exists():
        refs["v1"] = pd.read_csv(V1_CSV).groupby("file_name")["seed_ann_id"].first()
    if Path(F5_CSV).exists():
        f5 = pd.read_csv(F5_CSV)
        refs["f5"] = f5[f5["seed_index"] == 0].groupby("file_name")["seed_ann_id"].first()
    for name, ref in refs.items():
        bad, n = [], 0
        for fn in roi_files(images_dir):
            if fn not in ref.index:
                continue
            im = images[images["file_name"] == fn].iloc[0]
            gt_mit = ds.image_annotations(annotations, fn)
            gt_mit = gt_mit[gt_mit["category_id"] == ds.MITOTIC]
            seeds, _, _ = draw_seeds(gt_mit, BORDER,
                                     (int(im["height"]), int(im["width"])), 1)
            n += 1
            if int(seeds[0][1]["ann_id"]) != int(ref[fn]):
                bad.append(fn)
        rows.append({"check": "preflight_seed_identity", "label": name, "n": n,
                     "n_items": len(bad), "passed": bool(n > 0 and not bad),
                     "note": ";".join(bad)})
    if not rows:
        rows.append({"check": "preflight_seed_identity", "label": "-", "passed": False,
                     "note": "neither reproduction-gate input is present"})
    return pd.DataFrame(rows)


def ceiling_formula_gate(ceiling: pd.DataFrame) -> dict:
    """Pre-flight: the band predicate and `dist_to_interior` must reproduce section 2b.

    Section 2b lists seven band mitoses at si == 0 across v1's 7 ROIs, with their distances,
    and states there are no others. Both the predicate and the distance formula are used by
    the ceiling CSV and by every band count in the ledger, so validating them against a
    committed table is a five-second check on the geometry the whole experiment rests on.
    """
    want = {(201, 4494): 12.0, (201, 4502): 8.0, (246, 6612): 9.0, (246, 6618): 10.0,
            (301, 14779): 16.0, (301, 15022): 3.0, (402, 20275): 8.0}
    v1_rois = {"094.tiff", "201.tiff", "246.tiff", "301.tiff", "402.tiff", "459.tiff",
               "548.tiff"}
    got = {}
    sub = ceiling[(ceiling["seed_index"] == 0) & ceiling["file_name"].isin(v1_rois)
                  & (ceiling["ann_id"] >= 0)]
    for r in sub.itertuples():
        got[(int(r.file_name.split(".")[0]), int(r.ann_id))] = round(float(r.dist_to_interior_px), 1)
    return {"check": "preflight_ceiling_formula", "label": "section 2b's seven distances",
            "n": len(want), "n_items": len(got), "passed": bool(got == want),
            "note": "" if got == want else f"got {sorted(got.items())}"}


# ---------------------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spacing", action="store_true",
                    help="section 2b's ceiling geometry only, no response maps")
    ap.add_argument("--gates", action="store_true", help="pre-run gates only; writes nothing")
    ap.add_argument("--smoke", action="store_true", help="one ROI, one seed")
    ap.add_argument("--images-dir", default=IMAGES_DIR)
    ap.add_argument("--seeds", type=int, default=N_SEEDS)
    args = ap.parse_args(argv)

    # Mechanical guard: `results/f5_nms_radius_ablation*.csv` and
    # `results/tm_ccoeff_threshold_axis_sweep*.csv` are committed INPUTS to the reproduction
    # gates. This module must never write anything but its own five files.
    assert all(p.name.startswith("f6_padding_ablation") and p.parent == OUT
               for p in OUT_PATHS), f"an output path escaped the f6 namespace: {OUT_PATHS}"

    t00 = time.time()
    images, annotations = ds.load_annotations()
    ds.check_invariants(annotations)

    if args.spacing:
        OUT.mkdir(exist_ok=True)
        ceil = ceiling_table(images, annotations, args.images_dir, args.seeds)
        ceil.to_csv(OUT_CEILING, index=False)
        band = ceil[ceil["ann_id"] >= 0]
        print(f"wrote {OUT_CEILING} {ceil.shape}")
        print(f"  cells                       : {len(ceil.groupby(['file_name', 'seed_index']))}")
        print(f"  band mitoses (rows)         : {len(band)}")
        print(f"  unreachable (dist > radius) : {int(band['exceeds_match_radius'].sum())}")
        if len(band):
            print(f"  dist to interior, px        : "
                  f"{band['dist_to_interior_px'].min():.1f}-{band['dist_to_interior_px'].max():.1f}")
        print(f"  granularity 1/n_eval        : {ceil['granularity_1_over_n'].min():.4f}-"
              f"{ceil['granularity_1_over_n'].max():.4f}")
        g = ceiling_formula_gate(ceil)
        print(f"gate {g['check']}: passed={g['passed']} {g['note']}")
        return 0 if g["passed"] else 1

    if args.gates:
        # Writes NOTHING. The ceiling frame is built in memory only so its formula gate can
        # run; --spacing is what persists it.
        ceil = ceiling_table(images, annotations, args.images_dir, args.seeds)
        gs = [ceiling_formula_gate(ceil), od_hoist_gate(images, annotations, args.images_dir),
              interior_identity_gate(images, annotations, args.images_dir)]
        gs += seed_identity_gate(images, annotations, args.images_dir).to_dict("records")
        for g in gs:
            extra = ""
            if g.get("max_abs_diff") is not None and not pd.isna(g.get("max_abs_diff", np.nan)):
                extra = f" max_abs_diff={g['max_abs_diff']:.3e}"
            print(f"gate {g['check']:28s} [{g.get('label', '')}] passed={g['passed']}{extra}"
                  f" {g.get('note', '') or ''}".rstrip())
        print(f"total {time.time() - t00:.0f}s")
        return 0 if all(g["passed"] for g in gs) else 1

    OUT.mkdir(exist_ok=True)
    files = roi_files(args.images_dir)
    n_seeds = args.seeds
    if args.smoke:
        files, n_seeds = files[:1], 1

    checks, ledger_rows, frames, pool_sizes = [], [], [], {}

    ceil = ceiling_table(images, annotations, args.images_dir, n_seeds)
    ceil.to_csv(OUT_CEILING, index=False)
    checks.append(ceiling_formula_gate(ceil) if not args.smoke else
                  {"check": "preflight_ceiling_formula", "label": "skipped under --smoke",
                   "passed": True, "note": "smoke covers one ROI; the gate needs all 7"})
    checks += seed_identity_gate(images, annotations, args.images_dir).to_dict("records")
    checks.append(od_hoist_gate(images, annotations, args.images_dir))
    assert checks[-1]["passed"], "od hoist equivalence gate failed -- the optimisation is invalid"

    for fn in files:
        frames += run_roi(fn, images, annotations, args.images_dir, n_seeds,
                          checks, ledger_rows, pool_sizes, verify_shortcut=True)

    results = pd.concat(frames, ignore_index=True)
    cp.assert_floor_not_limiting(results)
    results.to_csv(OUT_MAIN, index=False)

    ledger = pd.DataFrame(ledger_rows)
    ledger.to_csv(OUT_LEDGER, index=False)

    sv = seedvar_table(results, pool_sizes)
    sv.to_csv(OUT_SEEDVAR, index=False)

    verif = pd.concat([pd.DataFrame(checks), join_integrity_gate(results),
                       reproduction_gate_a(results), reproduction_gate_b(results)],
                      ignore_index=True)
    verif.to_csv(OUT_VERIF, index=False)

    print(f"\nwrote {OUT_MAIN} {results.shape}")
    print(f"wrote {OUT_LEDGER} {ledger.shape}")
    print(f"wrote {OUT_SEEDVAR} {sv.shape}")
    print(f"wrote {OUT_CEILING} {ceil.shape}")
    print(f"wrote {OUT_VERIF} {verif.shape}")
    failed = verif[~verif["passed"].astype(bool)] if "passed" in verif else verif.iloc[:0]
    print(f"checks: {len(verif) - len(failed)}/{len(verif)} passed")
    print(verif[verif["check"].isin(["v1_reproduction_A", "f5_r75_reproduction_B",
                                     "interior_identity", "dropband_subset",
                                     "join_integrity", "dedup_key",
                                     "preflight_seed_identity", "preflight_ceiling_formula",
                                     "od_hoist_equivalence"])]
          .groupby("check")["passed"].agg(["size", "sum"]).to_string())
    if len(failed):
        print(failed.to_string())
    print(f"total {time.time() - t00:.0f}s")
    return 1 if len(failed) else 0


if __name__ == "__main__":
    sys.exit(main())
