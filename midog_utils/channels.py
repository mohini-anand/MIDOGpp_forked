"""Turning an RGB H&E crop into the single-channel image the matcher searches. Every
channel here returns "more object -> higher value", in float32."""

from __future__ import annotations

import cv2
import numpy as np
from skimage.color import rgb2hed

from . import chromatin as cm


def to_gray_inverted(rgb: np.ndarray) -> np.ndarray:
    """255 - grayscale."""
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    return (255.0 - gray).astype(np.float32)


def to_hematoxylin(rgb: np.ndarray) -> np.ndarray:
    """Hematoxylin channel via colour deconvolution, rescaled to 0-255 and clipped."""
    hed = rgb2hed(rgb.astype(np.float32) / 255.0)
    h = hed[:, :, 0].astype(np.float32)
    lo, hi = np.percentile(h, [0.5, 99.5])
    if hi <= lo:
        return np.zeros(h.shape, dtype=np.float32)
    return (np.clip((h - lo) / (hi - lo), 0.0, 1.0) * 255.0).astype(np.float32)


def to_hematoxylin_od(rgb: np.ndarray) -> np.ndarray:
    """Unclipped hematoxylin optical density -- `chromatin.hematoxylin_od`, as a channel."""
    return cm.hematoxylin_od(rgb)


def to_rgb(rgb: np.ndarray) -> np.ndarray:
    """Raw RGB, unchanged apart from dtype."""
    return rgb.astype(np.float32)


CHANNELS = {
    "gray_inverted": to_gray_inverted,
    "hematoxylin": to_hematoxylin,
    "hematoxylin_od": to_hematoxylin_od,
    "rgb": to_rgb,
}


def to_channel(rgb: np.ndarray, name: str) -> np.ndarray:
    if name not in CHANNELS:
        raise KeyError(f"unknown channel {name!r}; options are {sorted(CHANNELS)}")
    return CHANNELS[name](rgb)
