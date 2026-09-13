"""
    Production pipeline: one gated seed in, ranked detections on the ROI out.

    Wires DECISIONS.md D1 (TM_CCOEFF), D3 (unclipped hematoxylin_od), D7 (NMS radius =
    evaluate.radius_px(mpp)), D8 (seed_selection.build_seed's gated/recentred template) and
    D9 (max_peaks=100) into one callable, replacing the by-hand reimplementation every
    production_seed_precision_at_k* notebook carried. `find_and_suppress`, `build_seed`,
    `chromatin_density` and `radius_px` are unchanged; this module only names their
    production configuration as code.

    `chromatin_od` is an opt-in ranking axis (`rank_key='chromatin_od'`), not the default.
    DECISIONS.md D5 bars a chromatin statistic from becoming the *default* ranker without
    further validation on images/extra_valid; it does not bar it from existing as a
    selectable axis, and this module does not change that ledger.
"""

from __future__ import annotations

import cv2
import numpy as np
import pandas as pd

from . import channels as ch
from . import chromatin as cm
from . import evaluate as ev
from . import invariants as inv
from . import template_match as tm
from .find_and_suppress import FSConfig, find_and_suppress

CHANNEL = "hematoxylin_od"        # D3
TM_METHOD = cv2.TM_CCOEFF         # D1
PEAK_MIN_DISTANCE = 7
SELF_HIT_RADIUS = 5.0
DEEP_FLOOR_Z = -1.5
MAX_PEAKS = 100                   # D9, 2026-09-12 -- replaces the 2,000,000 left un-wired
                                   # in production_seed_precision_at_k_chromatin_half_pix_fix.ipynb
OD_WINDOW = tm.BASE_SIZE          # chromatin.chromatin_density's own default window (od51)

# rank_key -> the ranked-detections column that axis sorts on, descending. D5: 'tm_score'
# is the production default; 'chromatin_od' is opt-in only.
AXES = {"tm_score": "score", "chromatin_od": "od"}


def run_production_pipeline(rgb: np.ndarray, seed, mpp: float, *,
                            rank_key: str = "tm_score",
                            max_peaks: int = MAX_PEAKS) -> tuple[pd.DataFrame, dict]:
    """
        One click's gated seed -> the ranked candidate list for its ROI.

        ``seed`` is duck-typed on ``seed.template_xy``/``seed.base_size`` only -- typically a
        `seed_selection.build_seed` result. ``template_xy``/``base_size`` are where and how
        large the search template is cut (D8's gated/recentred point); this function never
        reads ``seed.click_xy`` or any other `Seed` field -- ground truth exclusion and
        match-radius computation against the click stay the caller's job, same as everywhere
        else in this repo. ``mpp`` is that ROI's own microns-per-pixel (`dataset.roi_mpp`);
        the NMS radius is derived from it (D7), never hardcoded.

        Returns ``(detections, info)``. ``detections`` is ranked best-first by
        ``AXES[rank_key]``, columns ``rank, cx, cy, score, angle, flip, scale``, plus ``od``
        when ``rank_key='chromatin_od'`` (``info['od_computed']`` records whether it was).
        Caller does ``.head(K)`` for a budget -- no separate top-K stage here, matching
        `compare.evaluate_arms`'s truncate-at-the-end convention.
    """
    if rank_key not in AXES:
        raise KeyError(f"rank_key must be one of {sorted(AXES)}, got {rank_key!r}")

    hem = ch.to_channel(rgb, CHANNEL) # run skimage's color deconvolution and keep only the hematoxylin channel
    nms_radius = ev.radius_px(mpp) # determine the NMS radius in pixels based on a selected ROI's microns-per-pixel value

    # This is the main location to tune any hyperparameters for the production pipeline
    cfg = FSConfig(channel=CHANNEL, base_size=seed.base_size, scales=(1.0,), n_angles=1,
                  flips=(False,), peak_min_distance=PEAK_MIN_DISTANCE,
                  max_peaks=max_peaks, nms_radius=nms_radius, self_hit_radius=SELF_HIT_RADIUS,
                  deep_floor_z=DEEP_FLOOR_Z, border_pad=True, tm_method=TM_METHOD)

    detections, info = find_and_suppress(hem, seed.template_xy, cfg, nms_radius=nms_radius)

    # D9: this cap is *meant* to bind on nearly every ROI (~17k pre-NMS peaks -> 100) and
    # is validated at zero cost through K=30 -- record it, don't invariant-check it.
    info["max_peaks_binding"] = bool(info["n_peaks"] == max_peaks)
    inv.check_nms_radius(nms_radius, mpp, label="production_pipeline")

    # A dormant tripwire: max_detections (10**9 by default) is never meant to bind, so
    # this can never fire at shipped defaults. It is unrelated to compare.py's separate,
    # already-documented Arm.caps vacuity -- this module doesn't use compare.Arm at all.
    inv.check_no_cap(int(info["n_detections"]), caps=(cfg.max_detections,), label="production_pipeline")

    # 'od' costs an extra padded-window pass, so only pay for it when ranking needs it.
    if rank_key == "chromatin_od":
        pad = OD_WINDOW // 2  # replicate-pad so a border detection's window is still fully readable
        hem_padded = cv2.copyMakeBorder(hem, pad, pad, pad, pad, borderType=cv2.BORDER_REPLICATE)
        shifted = detections.assign(cx=detections["cx"] + pad, cy=detections["cy"] + pad)
        scored = cm.score_detections(shifted, hem_padded, window=OD_WINDOW)
        detections = detections.assign(od=scored["od"].to_numpy())
    info["od_computed"] = bool(rank_key == "chromatin_od")

    key = AXES[rank_key]
    # Stable sort, NaN last -- same convention as compare._rank / chromatin.rerank,
    # duplicated locally (rather than importing chromatin.rerank) so that module's own
    # "nothing in the pipeline calls this" docstring line (D5) stays accurate.
    detections = (detections.sort_values(key, ascending=False, na_position="last", kind="mergesort")
                            .reset_index(drop=True))
    detections = detections.assign(rank=np.arange(len(detections)))
    info["rank_key"] = rank_key
    return detections, info
