"""Comparators that make the "matched nothing" bucket interpretable.

Without a baseline, a precision number from one-shot template matching says nothing:
the real question is whether one human click beats *any dark blob*, and only a nucleus
detector answers it.

Note there was nothing to reuse from `explore_dataset.ipynb` here. `tighten_box_otsu`
takes a crop *and a known centre* and returns a refined box -- it is a refinement
function, not a detector. And the obvious extension does not work either: Otsu over a
whole 2 mm^2 H&E ROI separates dark tissue from pale tissue, not nucleus from cytoplasm,
and returns one enormous connected component. Nucleus segmentation therefore uses a
*tiled* Otsu, below.
"""

from __future__ import annotations

import cv2
import numpy as np
import pandas as pd
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops

from .channels import to_gray_inverted, to_hematoxylin

# Slide glass under a brightfield scanner saturates near white; stained tissue, however
# pale, does not. 220 was validated across all 14 downloaded ROIs: it excludes 0 of the
# 690 mitotic annotations in them and yields tissue fractions of 0.67-1.00.
TISSUE_GRAY_MAX = 220


def tissue_mask(rgb: np.ndarray, gray_max: int = TISSUE_GRAY_MAX, close_px: int = 9) -> np.ndarray:
    """Tissue vs. slide glass, by a fixed brightness cut on grayscale.

    **Not Otsu.** Otsu assumes the histogram is bimodal with a real background mode, and
    these ROIs do not have one -- they are near-solid tumour, with only 0.7-3.7% of
    pixels above gray 220 (001-506 measured). Global Otsu on inverted grayscale
    consequently lands at gray 119-169 and splits *tissue into dark and pale halves*,
    keeping 38-73% of the ROI. On 301.tiff that put 31 of 218 mitotic figures outside the
    "tissue" mask, where `nucleus_blobs` cannot see them, and turned `random_in_tissue`
    from a uniform floor into a chromatin-biased one.

    Only a morphological *closing* is applied, to fill small bright holes inside tissue.
    The reference-style opening is deliberately absent: an opening with a 9 px element
    erases isolated dark objects smaller than the element -- exactly what a nucleus
    detector is looking for.
    """
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    mask = (gray < gray_max).astype(np.uint8)
    if close_px:
        k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (close_px, close_px))
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, k)
    return mask.astype(bool)


def _tiled_otsu_threshold(chan: np.ndarray, mask: np.ndarray, tile: int = 512) -> np.ndarray:
    """Per-tile Otsu, bilinearly upsampled into a full-resolution threshold map."""
    h, w = chan.shape
    ny, nx = int(np.ceil(h / tile)), int(np.ceil(w / tile))
    grid = np.full((ny, nx), np.nan, dtype=np.float32)
    for j in range(ny):
        for i in range(nx):
            sub = chan[j * tile:(j + 1) * tile, i * tile:(i + 1) * tile]
            sub_mask = mask[j * tile:(j + 1) * tile, i * tile:(i + 1) * tile]
            vals = sub[sub_mask]
            if vals.size < 1000 or vals.max() <= vals.min():
                continue
            grid[j, i] = threshold_otsu(vals)
    if np.isnan(grid).all():
        grid[:] = threshold_otsu(chan[mask]) if mask.any() else 0.0
    else:
        grid[np.isnan(grid)] = np.nanmedian(grid)
    return cv2.resize(grid, (w, h), interpolation=cv2.INTER_LINEAR)


def nucleus_blobs(
    rgb: np.ndarray,
    mask: np.ndarray = None,
    tile: int = 512,
    min_area: int = 80,
    max_area: int = 4000,
    max_detections: int = None,
) -> pd.DataFrame:
    """Ranked "any dark nucleus" detector: hematoxylin -> tiled Otsu -> components.

    Components are ranked by mean hematoxylin intensity (darkest chromatin first) so the
    output is a ranked list and `recall@K` is defined for it. A connected-component set
    is otherwise unordered, which would make it incomparable with a scored detector.

    ``max_detections`` is ``None`` (no truncation) by default. It used to default to
    20000, which bound silently on the two densest ROIs: it made the "blob detector
    truncated to the matcher's budget" comparison impossible to satisfy (the matcher
    emits up to 23143), and it capped ``experiment.score_probe``'s nucleus population at
    exactly 20000 -- a truncation that keeps only the *darkest*, hardest competitors and
    so biases every rank and base-rate statistic derived from it.
    """
    mask = tissue_mask(rgb) if mask is None else mask
    h_chan = to_hematoxylin(rgb)
    thr = _tiled_otsu_threshold(h_chan, mask, tile)
    binary = (h_chan > thr) & mask

    lab = label(binary, connectivity=2)
    rows = []
    for p in regionprops(lab, intensity_image=h_chan):
        if not (min_area <= p.area <= max_area):
            continue
        cy, cx = p.centroid
        rows.append({"cx": cx, "cy": cy, "score": float(p.mean_intensity), "area": int(p.area)})

    df = pd.DataFrame(rows, columns=["cx", "cy", "score", "area"])
    if len(df):
        df = df.sort_values("score", ascending=False)
        if max_detections is not None:
            df = df.head(max_detections)
        df = df.reset_index(drop=True)
    df.insert(0, "rank", np.arange(len(df)))
    return df


def random_in_tissue(mask: np.ndarray, n: int, rng=None) -> pd.DataFrame:
    """Uniform random points restricted to tissue -- the floor.

    Restricting to tissue matters because uniform-over-ROI points can land on slide
    glass, which a real detector never does. On these particular ROIs that is a small
    correction -- they are 92-100% tissue under `tissue_mask` -- but it is the correct
    comparison, and it is the only thing that makes "lift over random" mean anything.

    Sampling is over ``np.flatnonzero(mask)`` rather than ``np.nonzero(mask)``, which
    returns *two* 35-million-element index arrays where one will do. (``rng.choice(N, k,
    replace=False)`` itself is cheap at this ratio -- numpy's Generator uses Floyd's
    algorithm, measured at 0.00 s and no extra allocation for k=20000 out of N=35e6 -- so
    the saving is in the index arrays, not in the draw.)
    """
    rng = np.random.default_rng(0) if rng is None else rng
    flat = np.flatnonzero(mask.ravel())
    if flat.size == 0:
        return pd.DataFrame(columns=["rank", "cx", "cy", "score"])
    n = min(n, flat.size)
    pick = flat[rng.choice(flat.size, size=n, replace=False)] if n < flat.size else flat
    ys, xs = np.divmod(pick, mask.shape[1])
    return pd.DataFrame(
        {
            "rank": np.arange(n),
            "cx": xs.astype(float),
            "cy": ys.astype(float),
            "score": np.linspace(1.0, 0.0, n),  # arbitrary but stable ordering
        }
    )


def grid_lattice(roi_shape, step: int) -> pd.DataFrame:
    """A regular lattice of candidate points over the ROI -- no seed, no detector.

    The point of this comparator is that it contains no image evidence whatsoever: it is
    the geometry of the evaluation alone. Any generator that does not beat it at matched
    budget is not contributing detection, and a lattice ranked by a per-point statistic
    isolates how much of a pipeline's performance is the *ranker* rather than the search.

    The lattice is offset by ``step // 2`` so it does not start on the ROI edge, and no
    tissue restriction is applied -- these ROIs are 92-100% tissue under `tissue_mask` and
    a lattice point on glass carries a low chromatin value, so it sinks in any ranking
    rather than needing to be excluded.

    ``step`` is a real hyperparameter and must be swept, not fixed: measured on 301.tiff a
    *finer* lattice is worse (step 20: recall@budget 0.825, read-50 177; step 10: 0.765,
    356), because a denser lattice finds a better-optimised maximum of the chromatin field
    inside each match-radius disc and that field's true maxima are dense stromal and
    nuclear clumps rather than mitoses. A single step reported as "the grid baseline" is a
    tuned hyperparameter wearing a baseline's clothes.

    Returns ``cx, cy`` only. The caller attaches the ranking statistic and applies
    `nms.nms_by_distance` at that image's match radius, ordered by the same key the arm is
    ranked by -- suppression ordered by anything else silently reranks the arm.
    """
    h, w = int(roi_shape[0]), int(roi_shape[1])
    off = int(step) // 2
    xs = np.arange(off, w, int(step), dtype=np.float64)
    ys = np.arange(off, h, int(step), dtype=np.float64)
    gx, gy = np.meshgrid(xs, ys)
    return pd.DataFrame({"cx": gx.ravel(), "cy": gy.ravel()})
