"""
    One gated seed in, ranked detections on the ROI out -- the production click-to-verify
    pipeline. Wraps `find_and_suppress` with the production channel, matcher, and radii.
    `rank_key="chromatin_od"` is an opt-in second ranking axis; `"tm_score"` is the default.
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
MAX_PEAKS = 100                   # D9
OD_WINDOW = tm.BASE_SIZE

AXES = {"tm_score": "score", "chromatin_od": "od"}


def run_production_pipeline(rgb: np.ndarray, seed, mpp: float, *, rank_key: str = "tm_score", max_peaks: int = MAX_PEAKS) -> tuple[pd.DataFrame, dict]:
    """
        One click's gated seed -> the ranked candidate list for its ROI.

        rgb (np.ndarray): the ROI's full-resolution RGB image.
        seed: a `seed_selection.build_seed` result; reads only template_xy and base_size.
        mpp (float): the ROI's microns-per-pixel (`dataset.roi_mpp`).
        rank_key (str): "tm_score" (default) or "chromatin_od" -- which column ranks the list.
        max_peaks (int): pre-NMS peak cap passed to `find_and_suppress`.

        Returns tuple[pd.DataFrame, dict]: ranked detections (rank, cx, cy, score, angle,
        flip, scale, plus od when rank_key is "chromatin_od") and a per-stage info dict.
    """
    if rank_key not in AXES:
        raise KeyError(f"rank_key must be one of {sorted(AXES)}, got {rank_key!r}")

    hem = ch.to_channel(rgb, CHANNEL)
    nms_radius = ev.radius_px(mpp)

    cfg = FSConfig(base_size=seed.base_size, scales=(1.0,), n_angles=1, flips=(False,), peak_min_distance=PEAK_MIN_DISTANCE, max_peaks=max_peaks, nms_radius=nms_radius, self_hit_radius=SELF_HIT_RADIUS, deep_floor_z=DEEP_FLOOR_Z, border_pad=True, tm_method=TM_METHOD)

    detections, info = find_and_suppress(hem, seed.template_xy, cfg, nms_radius=nms_radius)

    info["max_peaks_binding"] = bool(info["n_peaks"] == max_peaks)
    inv.check_nms_radius(nms_radius, mpp, label="production_pipeline")

    if rank_key == "chromatin_od":
        pad = OD_WINDOW // 2
        hem_padded = cv2.copyMakeBorder(hem, pad, pad, pad, pad, borderType=cv2.BORDER_REPLICATE)
        shifted = detections.assign(cx=detections["cx"] + pad, cy=detections["cy"] + pad)
        scored = cm.score_detections(shifted, hem_padded, window=OD_WINDOW)
        detections = detections.assign(od=scored["od"].to_numpy())
    info["od_computed"] = bool(rank_key == "chromatin_od")

    key = AXES[rank_key]
    detections = detections.sort_values(key, ascending=False, na_position="last", kind="mergesort").reset_index(drop=True)
    detections = detections.assign(rank=np.arange(len(detections)))
    info["rank_key"] = rank_key
    return detections, info
