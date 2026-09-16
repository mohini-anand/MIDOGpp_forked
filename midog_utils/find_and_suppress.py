"""
    Template-match search from one seed point. `FSConfig` holds the search settings;
    `find_and_suppress` cuts a template, correlates it against the ROI, extracts peaks,
    suppresses overlaps, drops the seed's own self-hit, and returns detections ranked
    best-first. See `midog_utils/FIND_AND_SUPPRESS_REFERENCE_DIFFS.md` for how this differs
    from the reference implementation it replaces.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Optional

import cv2
import numpy as np
import pandas as pd

from . import template_match as tm
from .nms import nms_by_distance


@dataclass
class FSConfig:
    channel: str = "gray_inverted" # production overrides this to "hematoxylin_od"; this default is only live in experiment.py's bare FSConfig() calls
    base_size: int = tm.BASE_SIZE # longer side of refined template, actual template becomes square 
    patch_size: int = tm.PATCH_SIZE # size of the patch to read from the image, ceil(BASE_SIZE * sqrt(2))
    scales: tuple = (1.0,)
    n_angles: int = 1
    flips: tuple = (False,)
    peak_min_distance: int = 7  # local-maxima window: k = 2 * peak_min_distance + 1
    score_threshold: float = 0.5
    deep_floor_z: Optional[float] = None  # robust-z floor (median + z*MAD); None = use score_threshold
    max_peaks: int = 250000
    nms_radius: float = None  # None = caller must supply, normally evaluate.radius_px(mpp)
    self_hit_radius: float = 5.0
    border_pad: bool = False  # replicate-pad by half the template size so `valid` reaches the ROI edge
    max_detections: int = 10 ** 9
    scale_normalize: bool = False  # see template_match.fused_response
    tm_method: int = cv2.TM_CCOEFF_NORMED

    @property
    def n_augmentations(self) -> int:
        return len(self.scales) * self.n_angles * len(self.flips)


def find_and_suppress(img_channel: np.ndarray, seed_xy, cfg: FSConfig = None, nms_radius=None):
    """
        Run one single-pass search from one seed point.

        img_channel (np.ndarray): the converted search channel for the whole ROI, 2D float32.
        seed_xy (tuple[float, float]): (x, y) pixel to centre the template on.
        cfg (FSConfig or None): search settings; defaults to FSConfig().
        nms_radius (float or None): overrides cfg.nms_radius; one of the two must resolve.

        Returns tuple[pd.DataFrame, dict]: detections ranked best-first (columns rank, cx,
        cy, score, angle, flip, scale) and a per-stage info dict of counts and timings.
    """
    cfg = cfg or FSConfig()
    seed_x, seed_y = float(seed_xy[0]), float(seed_xy[1])
    nms_radius = cfg.nms_radius if nms_radius is None else nms_radius
    if nms_radius is None:
        raise ValueError(
            "nms_radius is unset. Pass the image's evaluation match radius "
            "(evaluate.radius_px(mpp)), or set FSConfig(nms_radius=...) explicitly."
        )
    info = {"nms_radius": round(float(nms_radius), 2)}

    patch = tm.read_padded_patch(img_channel, seed_x, seed_y, cfg.patch_size)
    if patch is None:
        raise ValueError(
            f"seed ({seed_x:.0f}, {seed_y:.0f}) is within {cfg.patch_size // 2} px of the "
            "ROI border, so no rotation-safe template can be read"
        )

    t0 = time.time()
    templates, metas = tm.build_augmentations(patch, cfg.base_size, cfg.scales, cfg.n_angles, cfg.flips)
    info["n_augmentations"] = len(templates)

    if cfg.border_pad:
        pad = max((t.shape[0] - 1) // 2 for t in templates)
        h, w = img_channel.shape[:2]
        padded = cv2.copyMakeBorder(img_channel, pad, pad, pad, pad, borderType=cv2.BORDER_REPLICATE)
        fused_p, best_p, valid_p = tm.fused_response(padded, templates, cfg.scale_normalize, method=cfg.tm_method)
        fused, best, valid = (fused_p[pad:pad + h, pad:pad + w],
                              best_p[pad:pad + h, pad:pad + w],
                              valid_p[pad:pad + h, pad:pad + w])
        assert valid.all(), "border_pad left part of the ROI unreachable -- pad derivation is wrong"
        info["pad_px"] = pad
    else:
        fused, best, valid = tm.fused_response(img_channel, templates, cfg.scale_normalize, method=cfg.tm_method)
    info["t_match_s"] = round(time.time() - t0, 2)

    threshold = cfg.score_threshold
    if cfg.deep_floor_z is not None:
        med, mad = tm.robust_stats(fused, valid)
        threshold = med + cfg.deep_floor_z * mad
        info["deep_floor_median"], info["deep_floor_mad"] = round(med, 5), round(mad, 5)
    info["score_threshold_used"] = round(float(threshold), 5)

    t0 = time.time()
    centers, scores = tm.extract_peaks(fused, valid, cfg.peak_min_distance, threshold, cfg.max_peaks)
    info["n_peaks"] = len(centers)
    info["max_peak_score"] = float(scores[0]) if len(scores) else float("nan")

    keep = nms_by_distance(centers, scores, nms_radius)
    centers, scores = centers[keep], scores[keep]
    info["n_after_nms"] = len(centers)

    if len(centers):
        d_seed = np.hypot(centers[:, 0] - seed_x, centers[:, 1] - seed_y)
        self_hit = d_seed <= cfg.self_hit_radius
        info["n_self_hits"] = int(self_hit.sum())
        info["seed_self_score"] = float(scores[self_hit].max()) if self_hit.any() else float("nan")
        centers, scores = centers[~self_hit], scores[~self_hit]
    else:
        info["n_self_hits"] = 0
        info["seed_self_score"] = float("nan")

    centers, scores = centers[: cfg.max_detections], scores[: cfg.max_detections]
    info["n_detections"] = len(centers)
    info["max_detection_score"] = float(scores[0]) if len(scores) else float("nan")
    info["t_postprocess_s"] = round(time.time() - t0, 2)

    aug = [metas[int(best[int(y), int(x)])] for x, y in centers]
    detections = pd.DataFrame({
        "rank": np.arange(len(centers)),
        "cx": centers[:, 0],
        "cy": centers[:, 1],
        "score": scores,
        "angle": [a.angle for a in aug],
        "flip": [a.flip for a in aug],
        "scale": [a.scale for a in aug],
    })
    return detections, info
