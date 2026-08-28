"""Turning an RGB H&E crop into the single-channel image the matcher searches.

The reference implementation assumed bright particles on a dark background and its
whole edge/threshold core is written for that polarity. H&E is the opposite: dark
chromatin on light tissue. Both channels here return "more object -> higher value",
which is the convention `explore_dataset.ipynb` adopted for the same reason.

Everything returns float32 because `cv2.matchTemplate` with `TM_CCOEFF_NORMED`
requires it.
"""

from __future__ import annotations

import cv2
import numpy as np
from skimage.color import rgb2hed


def to_gray_inverted(rgb: np.ndarray) -> np.ndarray:
    """255 - grayscale. Cheap, stain-agnostic, and the notebook's convention."""
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    return (255.0 - gray).astype(np.float32)


def to_hematoxylin(rgb: np.ndarray) -> np.ndarray:
    """Hematoxylin channel via colour deconvolution, rescaled to 0-255.

    Chromatin-selective and largely insensitive to eosin, so it should travel better
    across the four scanners in this dataset than raw grayscale does.
    """
    hed = rgb2hed(rgb.astype(np.float32) / 255.0)
    h = hed[:, :, 0].astype(np.float32)
    lo, hi = np.percentile(h, [0.5, 99.5])
    if hi <= lo:
        return np.zeros(h.shape, dtype=np.float32)
    return (np.clip((h - lo) / (hi - lo), 0.0, 1.0) * 255.0).astype(np.float32)


def to_rgb(rgb: np.ndarray) -> np.ndarray:
    """Raw RGB, unchanged apart from dtype -- the third search-channel variant.

    No stain separation and no inversion, so scanner-specific colour variation enters
    the match directly (`Research Logs/design_choices.md`, section 3). Requires
    `template_match.fused_response`'s multi-channel `cv2.matchTemplate` path, which
    collapses to the single-channel case unchanged when a caller passes a 2-D image.
    """
    return rgb.astype(np.float32)


CHANNELS = {
    "gray_inverted": to_gray_inverted,
    "hematoxylin": to_hematoxylin,
    "rgb": to_rgb,
}


def to_channel(rgb: np.ndarray, name: str) -> np.ndarray:
    if name not in CHANNELS:
        raise KeyError(f"unknown channel {name!r}; options are {sorted(CHANNELS)}")
    return CHANNELS[name](rgb)
