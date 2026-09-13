"""
    Template-match search from one seed point.

    Two things live here. `FSConfig` is the dataclass of everything one search run needs --
    which channel to search in, how the template is cut and augmented (scale / rotation /
    flip), the peak-extraction floor (a fixed `score_threshold` or an adaptive
    `deep_floor_z`), the NMS and self-hit radii, optional border padding, and which
    `cv2.TM_*` similarity to match with.

    `find_and_suppress` is the one function that runs a
    search: cut a template around `seed_xy`, build its augmentation bank, correlate against
    the whole image, extract peaks, suppress overlapping detections, drop the seed's own
    self-hit, and return the survivors ranked best-first.

    See `midog_utils/FIND_AND_SUPPRESS_REFERENCE_DIFFS.md` for how this deliberately differs
    from the reference `bbox tuning code reference/bbox_tuning.py:757` implementation it
    replaces.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional

import cv2
import numpy as np
import pandas as pd

from . import template_match as tm
from .nms import nms_by_distance


@dataclass
class FSConfig:
    channel: str = "gray_inverted"
    base_size: int = tm.BASE_SIZE
    patch_size: int = tm.PATCH_SIZE

    # Single-scale by default: the template is the annotation box itself, 51 px odd-sized
    # so the click is the centre pixel. Multi-scale fusion was measured and rejected --
    # taking a raw max across 31/41/51 px templates gave AUC 0.906 (mitosis vs. ordinary
    # nucleus, 002.tiff) against 0.940 for the 51 px template alone, because
    # TM_CCOEFF_NORMED has a wider null for smaller templates and so the small one wins the
    # max by chance. See results/fs_fusion_variants.csv.
    scales: tuple = (1.0,)

    # No rotation/flip augmentation by default: the simplest configuration for the next
    # run, deliberately -- not a claim that it scores better than the previous 12-angle
    # x 2-flip default. `experiment.AUGMENTATION_VARIANTS` still compares this against
    # 12 angles x 2 flips and a 4-angle/2-flip middle ground; on the one seed measured
    # so far (002.tiff ann 17) rot90_4angles_2flips actually had the best discrimination
    # of the three (0.068), the old 12x2 default second (0.037), and this no-augmentation
    # default worst (-0.034). See `Research Logs/design_choices.md`, section 6.
    n_angles: int = 1
    flips: tuple = (False,)
    peak_min_distance: int = 7 # window size for cheap grey-dilation to find local maxima; k = 2 * peak_min_distance + 1, before NMS

    # Raised from 0.25: the lower floor let the detection list tile most of the ROI
    # (coverage_frac 0.88-0.97 in the first single-pass run), at which point full-list
    # recall stops being evidence about the detector -- see evaluate.py's module
    # docstring. `evaluate.threshold_sweep`'s default range was trimmed to match.
    # `Research Logs/design_choices.md`, section 6.
    score_threshold: float = 0.5

    deep_floor_z: Optional[float] = None  # robust-z floor (median + z*MAD) in place of score_threshold; None = unchanged

    max_peaks: int = 250000     # above the ~173k theoretical max, so the cap never binds


    # None = "use this image's evaluation match radius", which the caller must supply --
    # see find_and_suppress(). A fixed radius is wrong here: the match radius is derived
    # per image from microns-per-pixel and runs 29.6-33.1 px across these scanners, so a
    # constant 25.0 let two detections 26 px apart both survive while both sat inside one
    # ground-truth object's radius. One was then credited as a true positive and the other
    # as a false positive on the very same object, inflating the FROC's FP axis with
    # duplicates.
    nms_radius: float = None

    self_hit_radius: float = 5.0

    border_pad: bool = False  # replicate-pad by half the template size so `valid` reaches the ROI edge (F6 pad_on); False = unchanged

    max_detections: int = 10 ** 9  # no truncation: report every detection the search found

    scale_normalize: bool = False  # see template_match.fused_response; changes score semantics

    # Which cv2.TM_* similarity the search uses. The default is the only one this repo had
    # ever run before `Research Logs/2026-09-01-tm-variant-sweep.md`. Changing it changes the
    # *units* of `score_threshold` -- TM_CCORR on hematoxylin OD spans [0.34, 10.7] where 0.5
    # is a permissive floor, and negated TM_SQDIFF spans [-14.7, -1.9e-6] where 0.5 admits
    # nothing at all. Any non-default method should set the floor in per-map z units via
    # `template_match.robust_stats`, which is what `tm_variant_sweep.py` does.
    tm_method: int = cv2.TM_CCOEFF_NORMED

    @property
    def n_augmentations(self) -> int:
        return len(self.scales) * self.n_angles * len(self.flips)


def find_and_suppress(img_channel: np.ndarray, seed_xy, cfg: FSConfig = None, nms_radius=None):
    """
        Run one single-pass search from one seed point.

        ``img_channel`` -- ``np.ndarray`` (2D, float32): the already-converted search
        channel for the whole ROI, e.g. ``channels.to_channel(rgb, cfg.channel)``.
        ``cfg.channel`` itself is never read here -- conversion is the caller's job.

        ``seed_xy`` -- ``(float, float)``: the ``(x, y)`` pixel to centre the template
        on. Usually ``seed_selection.build_seed(...).template_xy``.

        ``cfg`` -- ``FSConfig`` or ``None`` (defaults to ``FSConfig()``): every other
        search knob -- template size, augmentation, thresholds, radii, matcher.

        ``nms_radius`` -- ``float`` or ``None``: overrides ``cfg.nms_radius``. One of
        the two must resolve to a value, normally ``evaluate.radius_px(mpp)``.

        Returns ``(detections, info)``:

        ``detections`` -- ``pandas.DataFrame``, one row per surviving detection, ranked
        best-first (row 0 = best). Columns:

        * ``rank``  -- ``int``, 0-indexed position in the ranked list
        * ``cx``, ``cy`` -- ``float``, pixel coordinates in ``img_channel``
        * ``score`` -- ``float``, raw ``cfg.tm_method`` similarity at that peak
        * ``angle`` -- ``float`` degrees, which rotation augmentation won there
        * ``flip``  -- ``bool``, whether the winning augmentation was flipped
        * ``scale`` -- ``float``, which scale augmentation won there

        ``info`` -- ``dict``: per-stage counts and timings (peaks found, NMS survivors,
        threshold used, etc.), so a run can be audited without re-running it.

        ``cfg.deep_floor_z`` replaces ``score_threshold`` with an adaptive robust-z floor
        when set; 
        ``cfg.border_pad`` replicate-pads before matching so edge pixels are
        reachable (F6 ``pad_on``) when set. Both default to current behaviour.
    """
    cfg = cfg or FSConfig()
    seed_x, seed_y = float(seed_xy[0]), float(seed_xy[1]) # this is the CENTER, and NOT: bottom left, top right
    nms_radius = cfg.nms_radius if nms_radius is None else nms_radius
    if nms_radius is None:
        raise ValueError(
            "nms_radius is unset. Pass the image's evaluation match radius "
            "(evaluate.radius_px(mpp)), or set FSConfig(nms_radius=...) explicitly."
        )
    info = {"nms_radius": round(float(nms_radius), 2)} # gathering information about this pipeline starting with the NMS radius

    patch = tm.read_padded_patch(img_channel, seed_x, seed_y, cfg.patch_size) # usually the default or the tightened bbox dims just for safety before a rotation based augmentation 
    if patch is None:
        raise ValueError(
            f"seed ({seed_x:.0f}, {seed_y:.0f}) is within {cfg.patch_size // 2} px of the "
            "ROI border, so no rotation-safe template can be read"
        )

    t0 = time.time()
    templates, metas = tm.build_augmentations(
        patch, cfg.base_size, cfg.scales, cfg.n_angles, cfg.flips
    )
    info["n_augmentations"] = len(templates)

    if cfg.border_pad:  # F6 pad_on: reach the edge instead of leaving a valid=False border band; now we apply template matching 
        pad = max((t.shape[0] - 1) // 2 for t in templates)
        h, w = img_channel.shape[:2]
        padded = cv2.copyMakeBorder(img_channel, pad, pad, pad, pad, borderType=cv2.BORDER_REPLICATE) # ensure border pixels aren't ignored
        fused_p, best_p, valid_p = tm.fused_response(padded, templates, cfg.scale_normalize,
                                                     method=cfg.tm_method) # run template matching on the 'padded' image
        fused, best, valid = (fused_p[pad:pad + h, pad:pad + w],
                              best_p[pad:pad + h, pad:pad + w],
                              valid_p[pad:pad + h, pad:pad + w])
        assert valid.all(), "border_pad left part of the ROI unreachable -- pad derivation is wrong"
        info["pad_px"] = pad # record the amount of padding applied to the image for later information
    else:
        fused, best, valid = tm.fused_response(img_channel, templates, cfg.scale_normalize,
                                               method=cfg.tm_method)
    info["t_match_s"] = round(time.time() - t0, 2)

    threshold = cfg.score_threshold # not really used anymore, previously a fixed value but the next block now dynamically adjusts the threshold based on the ROI
    if cfg.deep_floor_z is not None:  # adaptive robust-z floor, computed on the cropped map
        med, mad = tm.robust_stats(fused, valid)
        threshold = med + cfg.deep_floor_z * mad
        info["deep_floor_median"], info["deep_floor_mad"] = round(med, 5), round(mad, 5)
    info["score_threshold_used"] = round(float(threshold), 5)

    t0 = time.time()
    centers, scores = tm.extract_peaks(
        fused, valid, cfg.peak_min_distance, threshold, cfg.max_peaks
    ) # trim down to our pre NMS candidate list, currently set to 100 for demo purposes
    info["n_peaks"] = len(centers)
    # Pre-removal, so this is always the seed's own self-correlation at ~1.0. Kept only
    # as an audit trail; `max_detection_score` below is the informative number.
    info["max_peak_score"] = float(scores[0]) if len(scores) else float("nan")

    keep = nms_by_distance(centers, scores, nms_radius)
    centers, scores = centers[keep], scores[keep]
    info["n_after_nms"] = len(centers)

    # Drop the seed's own detection. It is guaranteed to be there at score ~1.0 and it
    # is not a discovery.
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
    # The best *genuine* match -- i.e. what the search actually discovered, as opposed to
    # the seed rediscovering itself.
    info["max_detection_score"] = float(scores[0]) if len(scores) else float("nan")
    info["t_postprocess_s"] = round(time.time() - t0, 2)

    aug = [metas[int(best[int(y), int(x)])] for x, y in centers]
    detections = pd.DataFrame(
        {
            "rank": np.arange(len(centers)),
            "cx": centers[:, 0],
            "cy": centers[:, 1],
            "score": scores,
            "angle": [a.angle for a in aug],
            "flip": [a.flip for a in aug],
            "scale": [a.scale for a in aug],
        }
    )
    return detections, info
