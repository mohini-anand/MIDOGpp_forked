"""Overlays for eyeballing a run before trusting its numbers."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Rectangle

from .dataset import BOX_SIZE, LOOKALIKE, MITOTIC
from .evaluate import FP_LOOKALIKE, FP_UNANNOTATED, TP

BUCKET_COLORS = {TP: "#00c853", FP_LOOKALIKE: "#ff9100", FP_UNANNOTATED: "#2979ff"}


def overlay(rgb, gt, detections, seed_xy=None, ax=None, downsample=4, title="", top_n=None,
            radius=None):
    """Ground truth (red = mitotic, yellow = look-alike) with detections coloured by bucket.

    ``radius`` is the evaluation match radius in full-resolution pixels; each detection is
    drawn as a circle of exactly that radius, so what the picture shows is what the scorer
    counts as a hit. The previous fixed radius was written ``8 * s * downsample / 2``, which
    algebraically collapses to the constant 4 in display units regardless of ``downsample``
    -- roughly half the match radius, so detections looked further from the ground truth
    than the scorer considered them.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(14, 10))

    small = rgb[::downsample, ::downsample]
    ax.imshow(small)
    s = 1.0 / downsample
    half = BOX_SIZE / 2

    for _, a in gt.iterrows():
        color = "#d50000" if a["category_id"] == MITOTIC else "#ffea00"
        ax.add_patch(
            Rectangle(
                ((a["cx"] - half) * s, (a["cy"] - half) * s),
                BOX_SIZE * s, BOX_SIZE * s,
                fill=False, edgecolor=color, linewidth=0.8,
            )
        )

    det = detections if top_n is None else detections.head(top_n)
    r_disp = (BOX_SIZE / 2 if radius is None else float(radius)) * s
    for _, d in det.iterrows():
        color = BUCKET_COLORS.get(d.get("bucket", FP_UNANNOTATED), "#2979ff")
        ax.add_patch(Circle((d["cx"] * s, d["cy"] * s), radius=r_disp,
                            fill=False, edgecolor=color, linewidth=0.7))

    if seed_xy is not None:
        ax.plot(seed_xy[0] * s, seed_xy[1] * s, marker="*", markersize=18,
                color="white", markeredgecolor="black", markeredgewidth=0.8)

    ax.set_title(title, fontsize=10)
    ax.axis("off")
    return ax


def legend_handles():
    from matplotlib.lines import Line2D

    return [
        Line2D([], [], color="#d50000", lw=2, label="GT mitotic figure"),
        Line2D([], [], color="#ffea00", lw=2, label="GT look-alike (cat 2)"),
        Line2D([], [], color=BUCKET_COLORS[TP], lw=2, label="detection: TP"),
        Line2D([], [], color=BUCKET_COLORS[FP_LOOKALIKE], lw=2, label="detection: FP look-alike"),
        Line2D([], [], color=BUCKET_COLORS[FP_UNANNOTATED], lw=2, label="detection: FP unannotated"),
        Line2D([], [], color="white", marker="*", markeredgecolor="black", lw=0, label="seed"),
    ]


def froc_plot(curves, ax=None, title="FROC"):
    """curves: list of (label, fp_per_mm2, sensitivity)."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4.5))
    for label, fp, sens in curves:
        ax.plot(fp, sens, label=label, lw=1.5)
    ax.set_xscale("symlog", linthresh=1)
    ax.set_xlabel("false positives per mm$^2$")
    ax.set_ylabel("sensitivity (mitotic GT recovered)")
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title(title, fontsize=10)
    return ax
