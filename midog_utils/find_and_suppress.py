"""The MIDOG++ port of `find_and_supress`.

One human click becomes a template; the template is matched across the whole ROI with
rotation/scale augmentation; overlapping peaks are suppressed; what survives is a
ranked detection list.

Differences from `bbox tuning code reference/bbox_tuning.py:757` that change behaviour
rather than style:

* **No image masking.** The reference blacks out every existing annotation before
  correlating (`apply_mask_to_image(..., mode='black')`). On H&E that creates
  zero-variance windows -- where `TM_CCOEFF_NORMED` is undefined -- and hard edges that
  rotated templates correlate against. The seed's own detection is dropped afterwards
  instead, which has the same effect with no artefact. (Note the reference's `'none'`
  mode does not disable masking; it sets masked pixels to 1.)
* **The self-hit is dropped with a tight radius**, not the match radius. The minimum
  spacing between two MIDOG++ annotations anywhere in the dataset is 26.2 px (403.tiff;
  245.tiff is 26.6), below the ~30 px match radius, so dropping everything within a match
  radius of the seed could delete a legitimate detection of a neighbouring ground-truth
  object.
* **Scores are kept and used to rank.** The reference computes `match score` and then
  drops the column at bbox_tuning.py:798-801, before an NMS that is ordered by
  `x top left` rather than by score.
* **Boundary refinement, the stability filter, and the aspect-ratio / PSNR filters are
  not run.** Evaluation here is point-based, so box geometry moves no metric; and the
  stability check passes 0/9 objects on H&E (see the bbox-tuning research log), so
  keeping it would delete essentially every detection.
* **Nothing is random.** The reference picks its template with an unseeded
  `sample(frac=1)` and, when it has too many matches, keeps `max_templates` of them
  with an unseeded `random.sample` -- discarding the best matches at random. Here the
  seed is chosen by the caller and truncation is by score.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

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
    n_angles: int = 12
    flips: tuple = (False, True)
    peak_min_distance: int = 7
    score_threshold: float = 0.25
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
    max_detections: int = 10 ** 9  # no truncation: report every detection the search found
    scale_normalize: bool = False  # see template_match.fused_response; changes score semantics

    @property
    def n_augmentations(self) -> int:
        return len(self.scales) * self.n_angles * len(self.flips)


def find_and_suppress(img_channel: np.ndarray, seed_xy, cfg: FSConfig = None, nms_radius=None):
    """Run one single-pass search from one seed point.

    ``nms_radius`` overrides ``cfg.nms_radius``; one of the two must be set. There is
    deliberately no fallback default -- the correct value is the image's evaluation match
    radius, which depends on that ROI's microns-per-pixel and is therefore not knowable
    from the config alone. `experiment.run_one_image` passes it.

    Returns ``(detections, info)``. ``detections`` is ranked best-first with columns
    ``rank, cx, cy, score, angle, flip, scale``; ``info`` carries stage counts and
    timings so a run can be audited without re-running it.
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
    templates, metas = tm.build_augmentations(
        patch, cfg.base_size, cfg.scales, cfg.n_angles, cfg.flips
    )
    info["n_augmentations"] = len(templates)

    fused, best, valid = tm.fused_response(img_channel, templates, cfg.scale_normalize)
    info["t_match_s"] = round(time.time() - t0, 2)

    t0 = time.time()
    centers, scores = tm.extract_peaks(
        fused, valid, cfg.peak_min_distance, cfg.score_threshold, cfg.max_peaks
    )
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
