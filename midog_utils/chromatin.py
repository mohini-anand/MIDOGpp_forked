"""
    Chromatin density: absolute darkness under a detection, as a second ranking axis.
    Not the production ranker -- see `DECISIONS.md` D5. `chromatin_density`/`score_detections`
    back the opt-in `chromatin_od` axis in `production.run_production_pipeline`.
"""

from __future__ import annotations

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
            the unclipped `hematoxylin_od`).
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
