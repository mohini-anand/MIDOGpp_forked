"""
    Chromatin density: absolute darkness under a detection, as a second ranking axis.
    Not the production ranker -- see `DECISIONS.md` D5. Nothing in the production pipeline
    calls `rerank`; `chromatin_density`/`score_detections` back the opt-in `chromatin_od`
    axis in `production.run_production_pipeline`.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from . import template_match as tm

DEFAULT_FRAC = 0.10


def hematoxylin_od(rgb: np.ndarray) -> np.ndarray:
    """Unclipped hematoxylin optical density from colour deconvolution."""
    from skimage.color import rgb2hed
    return rgb2hed(rgb.astype(np.float32) / 255.0)[:, :, 0].astype(np.float32)


def chromatin_density(structural: np.ndarray, cx: float, cy: float, window: int = tm.BASE_SIZE, frac: float = DEFAULT_FRAC) -> float:
    """
        Mean of the darkest ``frac`` of pixels in a ``window``-sized box around a point.

        structural (np.ndarray): single-channel map, "more object -> higher value" (use
            `hematoxylin_od`, not `channels.to_hematoxylin` -- the latter clips).
        cx, cy (float): centre point.
        window (int): box side length.
        frac (float): fraction of darkest pixels averaged.

        Returns float: the statistic, or nan if the window can't be read at the border.
    """
    patch = tm.read_padded_patch(structural, cx, cy, window)
    if patch is None:
        return float("nan")
    flat = patch.ravel()
    k = max(1, int(frac * flat.size))
    return float(np.partition(flat, -k)[-k:].mean())


def score_detections(detections: pd.DataFrame, structural: np.ndarray, window: int = tm.BASE_SIZE, frac: float = DEFAULT_FRAC) -> pd.DataFrame:
    """Add an ``od`` column to a detection frame. Does not reorder it."""
    od = [chromatin_density(structural, float(cx), float(cy), window, frac)
          for cx, cy in zip(detections["cx"].to_numpy(), detections["cy"].to_numpy())]
    return detections.assign(od=od)


def rerank(detections: pd.DataFrame, structural: Optional[np.ndarray] = None, window: int = tm.BASE_SIZE, frac: float = DEFAULT_FRAC) -> pd.DataFrame:
    """
        Re-sort a detection list by chromatin density, best-first, and renumber ``rank``.
        Not the production ranker (`DECISIONS.md` D5).

        detections (pd.DataFrame): a ranked detection list.
        structural (np.ndarray or None): pass to compute ``od``, or omit if already present.
        window, frac: see `chromatin_density`.

        Returns pd.DataFrame: re-sorted by ``od`` descending, ``nan`` last.
    """
    if structural is not None:
        detections = score_detections(detections, structural, window, frac)
    if "od" not in detections.columns:
        raise ValueError("no 'od' column -- pass structural, or call score_detections first")
    out = detections.sort_values("od", ascending=False, na_position="last", kind="mergesort").reset_index(drop=True)
    return out.assign(rank=np.arange(len(out)))
