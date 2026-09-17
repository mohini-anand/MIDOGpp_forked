"""Turning an RGB H&E crop into the single-channel image the matcher searches. Every
channel here returns "more object -> higher value", in float32."""

from __future__ import annotations

import cv2
import numpy as np

from . import chromatin as cm


def to_gray_inverted(rgb: np.ndarray) -> np.ndarray:
    """255 - grayscale."""
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    return (255.0 - gray).astype(np.float32)


def to_hematoxylin_od(rgb: np.ndarray) -> np.ndarray:
    """Unclipped hematoxylin optical density -- `chromatin.hematoxylin_od`, as a channel."""
    return cm.hematoxylin_od(rgb)


CHANNELS = {
    "gray_inverted": to_gray_inverted,
    "hematoxylin_od": to_hematoxylin_od,
}


def to_channel(rgb: np.ndarray, name: str) -> np.ndarray:
    if name not in CHANNELS:
        raise KeyError(f"unknown channel {name!r}; options are {sorted(CHANNELS)}")
    return CHANNELS[name](rgb)
