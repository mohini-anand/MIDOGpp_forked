"""Chromatin density: the signal `TM_CCOEFF_NORMED` ranks far more weakly than it could.

`cv2.matchTemplate(..., TM_CCOEFF_NORMED)` mean-centres and L2-normalises *both* the
template and the window, so it is invariant to ``I -> aI + b``. A pale, low-contrast
structure with the same spatial pattern as the seed scores identically to a dark, dense
one -- verified directly: a 30%-contrast copy of a patch scores 0.999999 against the
original's 1.0.

*Calibrated*: "blind" overstates this and was retracted on 2026-08-31. The invariance is
exact only inside the 0-255 range -- clipping breaks it, zero-variance white windows score
0 (a crude brightness gate), and empirically the correlation score still reaches AUC
0.61-0.81 on TP-vs-unannotated. It is a much weaker ranker than chromatin density, not a
blind one.

That matters because `results/morph_diag_bhattacharyya.csv` measures ``mean_intensity``
as the strongest mitotic-vs-ordinary-nucleus feature in **all seven** domains
(0.93-3.70). ``solidity`` is runner-up in only 2 of the 7 domains, so "next-best
solidity at 0.15-1.05" -- an earlier version of this line -- was wrong and is retracted.
The search score therefore under-weights the best available discriminator, which is why every
operating-point experiment to date (score threshold, template size, channel,
augmentation count, NMS ordering) came back neutral: none of them added information to
a ranking function that cannot see the signal.

This module supplies that signal as a cheap post-hoc statistic over the detections the
search already produced.

Why "mean of the darkest `frac` of the window" and not the Otsu component
------------------------------------------------------------------------
The obvious choice is to reuse `seed_selection.tighten_box_otsu`, segment the object
under the detection, and average it. Measured against a gate-free window statistic on
the saved detections (`od_rerank_probe.py`), that is the *worse* option:

    AUC, true positive vs unannotated FP, every detection given a value

    image      comp (Otsu, top-10% fallback)   top-10% of window
    301.tiff             0.895                       0.914
    405.tiff             0.864                       0.959
    002.tiff             0.938                       0.987
    506.tiff             0.918                       0.972

`tighten_box_otsu` returns ``None`` on 9-12% of detections (the click's pixel is not
foreground, or the component fails the size/shape sanity check), and those are not a
random subset -- scoring only the ones it accepts flatters it (0.943 on 301.tiff), while
giving the rejects a fallback value costs it the same advantage again. The window
statistic is defined for every detection, needs no gate, and so costs no recall.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from . import template_match as tm

DEFAULT_FRAC = 0.10


def hematoxylin_od(rgb: np.ndarray) -> np.ndarray:
    """Unclipped hematoxylin optical density from colour deconvolution.

    `channels.to_hematoxylin` rescales to 0-255 between the ROI's 0.5/99.5 percentiles
    and clips, which saturates exactly the dense-chromatin pixels this module needs to
    rank. This returns `rgb2hed`'s hematoxylin channel as-is, so no value is clipped.
    """
    from skimage.color import rgb2hed
    return rgb2hed(rgb.astype(np.float32) / 255.0)[:, :, 0].astype(np.float32)


def chromatin_density(structural: np.ndarray, cx: float, cy: float,
                      window: int = tm.BASE_SIZE, frac: float = DEFAULT_FRAC) -> float:
    """Mean of the darkest ``frac`` of pixels in a ``window``-sized box around a point.

    ``structural`` is a single-channel map in "more object -> higher value" convention,
    so the densest chromatin is the *highest*-valued tenth. Use `hematoxylin_od`, **not**
    `channels.to_hematoxylin`: the latter min-max rescales by the ROI's own 0.5/99.5
    percentiles, which puts exactly 0.5% of pixels at the 255 ceiling. A window over
    dense chromatin then lies entirely above that percentile and this statistic returns
    255.000 exactly -- 66 of 002.tiff's detections tie at the ceiling, and its whole
    top-100 saturates, which makes the ranking there arbitrary rather than informative.
    `hematoxylin_od` is unclipped and cuts the largest tie block from 66 to 10.

    Note this is a *within-image* quantity. Optical density is not calibrated across
    scanners, so values are comparable within one ROI's ranking and not between ROIs.

    Returns ``nan`` when the point sits too close to the ROI border to read a full
    window, which the caller should rank last rather than drop.
    """
    patch = tm.read_padded_patch(structural, cx, cy, window)
    if patch is None:
        return float("nan")
    flat = patch.ravel()
    k = max(1, int(frac * flat.size))
    return float(np.partition(flat, -k)[-k:].mean())


def score_detections(detections: pd.DataFrame, structural: np.ndarray,
                     window: int = tm.BASE_SIZE, frac: float = DEFAULT_FRAC) -> pd.DataFrame:
    """Add an ``od`` column to a detection frame. Does not reorder it."""
    od = [chromatin_density(structural, float(cx), float(cy), window, frac)
          for cx, cy in zip(detections["cx"].to_numpy(), detections["cy"].to_numpy())]
    return detections.assign(od=od)


def rerank(detections: pd.DataFrame, structural: Optional[np.ndarray] = None,
           window: int = tm.BASE_SIZE, frac: float = DEFAULT_FRAC) -> pd.DataFrame:
    """Re-sort a detection list by chromatin density, best-first, and renumber ``rank``.

    Pass ``structural`` to compute the ``od`` column, or omit it if `score_detections`
    already added one. ``nan`` (border) sorts last.

    This *replaces* the correlation score as the ranking key rather than blending with
    it. Measured on the saved detections, a rank-sum of the two never beats ``od`` alone
    and is worse on 405.tiff (0.857 vs 0.864) and 002.tiff (0.921 vs 0.957), despite the
    two being largely independent (Spearman +0.18 to +0.37) -- the correlation score's
    independent contribution is noise with respect to the mitotic/non-mitotic
    distinction. This demotes the search as a *ranker*.

    It does **not** establish that the search earns its place as the *candidate
    generator*, which an earlier version of this line claimed. `premise_test.py`
    measured that question directly and the answer went the other way: at matched budget
    on 301.tiff and 246.tiff the seeded pipeline **wins on 0 of 5 seeds** against the best
    seedless generator (`nucleus_blobs`, or a chromatin-ranked lattice) at every budget --
    i.e. it loses on all five. That comparison picks the seedless arm post hoc per
    (ROI, budget), so it is an upper bound on the seedless side; the sign is not in doubt
    but the size is. See `Research Logs/2026-09-01-premise-test-audit.md` (M1) and
    `results/premise_report.md` section 3.
    """
    if structural is not None:
        detections = score_detections(detections, structural, window, frac)
    if "od" not in detections.columns:
        raise ValueError("no 'od' column -- pass structural, or call score_detections first")
    # Stable sort: ties keep the caller's incoming order (normally descending score)
    # rather than quicksort's arbitrary one, so a run is reproducible.
    out = detections.sort_values("od", ascending=False, na_position="last",
                                 kind="mergesort").reset_index(drop=True)
    return out.assign(rank=np.arange(len(out)))
