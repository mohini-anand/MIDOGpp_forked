"""Overlays for eyeballing a run before trusting its numbers."""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

from .dataset import BOX_SIZE, MITOTIC
from .evaluate import HUMAN_CORRECT_LABEL, HUMAN_REJECTED_LABEL, NON_HUMAN_FINDINGS

BUCKET_COLORS = {
    HUMAN_CORRECT_LABEL: "#00c853",
    HUMAN_REJECTED_LABEL: "#ff9100",
    NON_HUMAN_FINDINGS: "#2979ff",
}


def overlay(rgb, gt, detections, seed_xy=None, tpl_xy=None, ax=None, downsample=4, title="", top_n=None, radius=None):
    """
        Ground truth (red = mitotic, yellow = look-alike) with detections coloured by bucket.

        radius (float): evaluation match radius in full-resolution pixels; each detection
        is drawn as a circle of that radius, so the picture matches what the scorer counts
        as a hit.
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
        color = BUCKET_COLORS.get(d.get("bucket", NON_HUMAN_FINDINGS), "#2979ff")
        ax.add_patch(Circle((d["cx"] * s, d["cy"] * s), radius=r_disp,
                            fill=False, edgecolor=color, linewidth=0.7))

    if seed_xy is not None:
        ax.plot(seed_xy[0] * s, seed_xy[1] * s, marker="*", markersize=18,
                color="white", markeredgecolor="black", markeredgewidth=0.8)

    if tpl_xy is not None:  # the recentred template centre, drawn distinct from the click
        ax.plot(tpl_xy[0] * s, tpl_xy[1] * s, marker="P", markersize=11,
                color="#00e5ff", markeredgecolor="black", markeredgewidth=0.7)

    ax.set_title(title, fontsize=10)
    ax.axis("off")
    return ax


def draw_box(ax, box, color, lw=1.6, label=None, scale=1.0):
    """
        Draw a half-open ``(y0, y1, x0, x1)`` bbox, e.g. from `seed_selection.tighten_box_otsu`.

        scale (float): maps the box's own pixel frame to display coordinates -- 1.0 for a
        box already in the axes' frame, or 1/downsample to draw a full-resolution box on
        an `overlay()`-style thumbnail.
    """
    if box is None:
        return
    y0, y1, x0, x1 = box
    ax.add_patch(Rectangle(((x0 - 0.5) * scale, (y0 - 0.5) * scale),
                           (x1 - x0) * scale, (y1 - y0) * scale,
                           fill=False, edgecolor=color, linewidth=lw, label=label))
