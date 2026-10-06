"""
    Template-match search from one seed point. `FSConfig` holds the search settings;
    `find_and_suppress` cuts a template, blanks its own footprint out of a copy of the
    search channel so it can't match itself (and, on request, further squares: the boxes
    earlier session iterations stored), correlates against that copy, extracts peaks,
    suppresses overlaps, and returns detections ranked best-first. See
    `midog_utils/FIND_AND_SUPPRESS_REFERENCE_DIFFS.md` for how this differs from the
    reference implementation it replaces.
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
    base_size: int = tm.BASE_SIZE # longer side of refined template, actual template becomes square 
    patch_size: int = tm.PATCH_SIZE # size of the patch to read from the image, ceil(BASE_SIZE * sqrt(2))
    scales: tuple = (1.0,)
    n_angles: int = 1
    flips: tuple = (False,)
    peak_min_distance: int = 7  # local-maxima window: k = 2 * peak_min_distance + 1
    deep_floor_z: Optional[float] = None  # robust-z floor (median + z*MAD); required, find_and_suppress raises when None
    max_peaks: int = 250000
    nms_radius: float = None  # None = caller must supply, normally evaluate.radius_px(mpp)
    border_pad: bool = False  # replicate-pad by half the template size so `valid` reaches the ROI edge
    scale_normalize: bool = False  # see template_match.fused_response
    tm_method: int = cv2.TM_CCOEFF

    @property
    def n_augmentations(self) -> int:
        return len(self.scales) * self.n_angles * len(self.flips)


def blank_squares(search_channel: np.ndarray, squares, fill_value: float) -> tuple[int, int]:
    """
        Blank extra odd squares into a search channel, in place, with `template_match.blank_seed_square`'s footprint: ``side``
        pixels centred on ``(round(cx), round(cy))`` (round half to even), clipped to the image. For the boxes earlier
        iterations of a session stored (`annotation_frame.run_iteration`), written into the copy blank_seed_square returned.

        search_channel (np.ndarray): the search copy; modified in place.
        squares (iterable[tuple[float, float, int]]): one (cx, cy, side) per square; centres unrounded, side a positive odd integer.
        fill_value (float): the value written: blank_seed_square's fill, the unblanked channel's minimum.

        Returns tuple[int, int]: the number of squares, and the number of distinct pixels they cover.
        Raises ValueError on a non-finite centre or a side that is not a positive odd integer.
    """
    h, w = search_channel.shape[:2]
    mask = np.zeros((h, w), dtype=bool)
    n_squares = 0
    for cx, cy, side in squares:
        cx, cy = float(cx), float(cy)
        if not (np.isfinite(cx) and np.isfinite(cy)):
            raise ValueError(f"a blanked square's centre must be finite, got ({cx}, {cy})")
        if float(side) != int(side) or int(side) < 1 or int(side) % 2 != 1:
            raise ValueError(f"a blanked square's side must be a positive odd integer, got {side}")
        half = int(side) // 2
        ix, iy = int(round(cx)), int(round(cy))
        y0, y1 = max(0, iy - half), min(h - 1, iy + half)
        x0, x1 = max(0, ix - half), min(w - 1, ix + half)
        if y0 <= y1 and x0 <= x1:
            mask[y0:y1 + 1, x0:x1 + 1] = True
        n_squares += 1
    search_channel[mask] = fill_value
    return n_squares, int(mask.sum())


def find_and_suppress(img_channel: np.ndarray, seed_xy, cfg: FSConfig = None, nms_radius=None, *, extra_blanks=None):
    """
        Run one single-pass search from one seed point.

        img_channel (np.ndarray): the converted search channel for the whole ROI, 2D float32.
        seed_xy (tuple[float, float]): (x, y) pixel to centre the template on.
        cfg (FSConfig or None): search settings, deep_floor_z required; defaults to FSConfig().
        nms_radius (float or None): overrides cfg.nms_radius; one of the two must resolve.
        extra_blanks (list[tuple[float, float, int]] or None): further (cx, cy, side) squares to blank in the search copy, after the template is cut and beside the seed's own blank (`blank_squares`); None leaves the search and info exactly as before.

        Returns tuple[pd.DataFrame, dict]: detections ranked best-first (columns rank, cx,
        cy, score, angle, flip, scale) and a per-stage info dict of counts and timings;
        n_extra_blanks and n_extra_blanked_px are added only when extra_blanks is given.
    """
    cfg = cfg or FSConfig()
    seed_x, seed_y = float(seed_xy[0]), float(seed_xy[1])
    nms_radius = cfg.nms_radius if nms_radius is None else nms_radius
    if nms_radius is None:
        raise ValueError(
            "nms_radius is unset. Pass the image's evaluation match radius "
            "(evaluate.radius_px(mpp)), or set FSConfig(nms_radius=...) explicitly."
        )
    if cfg.deep_floor_z is None:
        raise ValueError("deep_floor_z is unset. The peak threshold is median + deep_floor_z * MAD of the response map; set FSConfig(deep_floor_z=...) explicitly (production uses -1.5).")
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

    # The template above is cut from the ORIGINAL, unblanked img_channel. Only now, after
    # that cut, do we blank a COPY for the search -- img_channel itself (and anything else
    # the caller derives from it, e.g. production.py's chromatin_od ranking) is untouched.
    search_channel, n_blanked = tm.blank_seed_square(img_channel, seed_x, seed_y, cfg.base_size)
    info["n_blanked_px"] = n_blanked
    if extra_blanks is not None:
        # Same copy, same fill: the boxes earlier session iterations stored, blanked after the cut so the template stays clean.
        info["n_extra_blanks"], info["n_extra_blanked_px"] = blank_squares(search_channel, extra_blanks, float(img_channel.min()))

    if cfg.border_pad:
        pad = max((t.shape[0] - 1) // 2 for t in templates)
        h, w = search_channel.shape[:2]
        padded = cv2.copyMakeBorder(search_channel, pad, pad, pad, pad, borderType=cv2.BORDER_REPLICATE)
        fused_p, best_p, valid_p = tm.fused_response(padded, templates, cfg.scale_normalize, method=cfg.tm_method)
        fused, best, valid = (fused_p[pad:pad + h, pad:pad + w],
                              best_p[pad:pad + h, pad:pad + w],
                              valid_p[pad:pad + h, pad:pad + w])
        assert valid.all(), "border_pad left part of the ROI unreachable -- pad derivation is wrong"
        info["pad_px"] = pad
    else:
        fused, best, valid = tm.fused_response(search_channel, templates, cfg.scale_normalize, method=cfg.tm_method)
    info["t_match_s"] = round(time.time() - t0, 2)

    med, mad = tm.robust_stats(fused, valid)
    threshold = med + cfg.deep_floor_z * mad
    info["deep_floor_median"], info["deep_floor_mad"] = round(med, 5), round(mad, 5)
    info["score_threshold_used"] = round(float(threshold), 5)

    t0 = time.time()
    centers, scores = tm.extract_peaks(fused, valid, cfg.peak_min_distance, threshold, cfg.max_peaks)
    info["n_peaks"] = len(centers)

    keep = nms_by_distance(centers, scores, nms_radius)
    centers, scores = centers[keep], scores[keep]
    info["n_after_nms"] = len(centers)

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
