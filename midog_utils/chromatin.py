"""Chromatin density: absolute darkness under a detection, as a second ranking axis.

.. warning::

   **This is NOT the production ranker. See `DECISIONS.md` D5 (2026-09-04).** Detections are
   ordered by the `TM_CCOEFF` match score; chromatin density is reported beside it, never
   instead of it, and no experiment may declare it a *primary* axis without measuring on its
   own data that it beats `tm_score`. Nothing in the pipeline calls `rerank` -- its only
   callers are the three superseded 2026-08-31 probes (`od_experiment.py`,
   `od_workload_ab.py`, `od_seed_sweep.py`).

   **The premise the rest of this docstring is written on has expired.** Everything below the
   next paragraph argues against `TM_CCOEFF_NORMED`, which `DECISIONS.md` D1 retired on
   2026-09-04, four days after this module was committed (`7c3af93`). `TM_CCOEFF` is *not*
   invariant to contrast: at contrast 0.7 / 0.3 / 0.1 it returns 0.700 / 0.300 / 0.101 of the
   full-contrast score where `TM_CCOEFF_NORMED` returns 0.999 / 0.992 / 0.930
   (`verify_chromatin_ranker.py` section 1). D1's distinction still holds -- the matcher reads
   *contrast*, this module reads *absolute darkness*, so they are correlated rather than
   duplicates -- but the justification below names a matcher no longer in use, and the measured
   marginal advantage over the current one is not distinguishable from zero: paired
   Delta(recall@250) = +0.032, 95% CI [-0.047, +0.112], p = 0.36, positive on 3 of 7 ROIs.

   The text below is kept verbatim as the record of why this was built. Read it as history.

Historical rationale (2026-08-31, superseded by D1 and D5)
----------------------------------------------------------
`cv2.matchTemplate(..., TM_CCOEFF_NORMED)` mean-centres and L2-normalises *both* the
template and the window, so it is invariant to ``I -> aI + b``. A pale, low-contrast
structure with the same spatial pattern as the seed scores identically to a dark, dense
one -- verified directly: a 30%-contrast copy of a patch scores 0.999999 against the
original's 1.0.

That matters because `results/morph_diag_bhattacharyya.csv` measures ``mean_intensity``
as the strongest mitotic-vs-ordinary-nucleus feature in **all seven** domains
(0.93-3.70, against next-best ``solidity`` at 0.15-1.05). The search score therefore
discards the best available discriminator by construction, which is why every
operating-point experiment to date (score threshold, template size, channel,
augmentation count, NMS ordering) came back neutral: none of them added information to
a ranking function that cannot see the signal.

    D5 correction: that file has **two** columns and this paragraph quotes one. Against
    *look-alikes* -- the class that survives to the operating point -- ``mean_intensity``
    reads 0.192 (mast cell, n=305), 0.161 (lymphosarcoma, n=169) and 0.332 (lung, n=30).
    The rationale above was established on the easy contrast.

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

    D5 correction: this argument is against the *gated* component only, and the gate is not
    load-bearing. `tp_fp_feature_extract.shape_features` measures the largest Otsu component
    with every gate removed and fails on **0.0%** of candidates in all seven ROIs. Gate-free,
    ``mask_od_mean`` beats ``od51`` on 2-class AUC (6 wins, 1 tie), on the look-alike contrast
    (5/7) and on ``read_95`` (6/7). The window statistic's advantage was over a gate, not over
    the component.

A third thing the 2026-08-31 work never swept: the window size
--------------------------------------------------------------
``window`` defaults to `tm.BASE_SIZE` = 51 because that is the annotation box, not because 51
was measured against anything. It is a poor choice on both counts available. A 31 px window
beats it on 2-class AUC in 6 of 7 ROIs and cuts median ``read_95`` from 4,488 to 1,336; and
between candidate pairs <= 25 px apart -- whose 51 px windows share about half their pixels --
``od51`` correlates 0.72-0.84 against ``od31``'s 0.44-0.63 and the match score's 0.45-0.64, so
at 51 px the statistic is substantially reading the neighbour rather than the object.

Every number in these three corrections comes from the 7-ROI densest-per-domain draw at one
seed. **D5 requires the 14 ROIs of `images/extra_valid/` for any re-measurement**, including
any run that would restore a chromatin statistic to primary. See `verify_chromatin_ranker.py`.
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

    **Not the production ranker** (`DECISIONS.md` D5). Nothing in the pipeline calls this; the
    live experiments rank by `score` and report `od` as a second axis. It is kept because the
    three 2026-08-31 probes call it and because the comparison is worth being able to re-run.

    Pass ``structural`` to compute the ``od`` column, or omit it if `score_detections`
    already added one. ``nan`` (border) sorts last.

    This *replaces* the correlation score as the ranking key rather than blending with
    it. **D5 note: every number in the rest of this docstring was measured against
    `TM_CCOEFF_NORMED`, which D1 retired.** Under `TM_CCOEFF` the two axes are
    indistinguishable at recall@250 (p = 0.36 clustered by ROI), so "the correlation score's
    independent contribution is noise" does not carry over and has not been re-measured.
    Measured on the saved detections, a rank-sum of the two never beats ``od`` alone
    and is worse on 405.tiff (0.857 vs 0.864) and 002.tiff (0.921 vs 0.957), despite the
    two being largely independent (Spearman +0.18 to +0.37) -- the correlation score's
    independent contribution is noise with respect to the mitotic/non-mitotic
    distinction. The search still earns its place as the *candidate generator*; it is
    only demoted as a ranker.
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
