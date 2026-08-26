"""Comparators that make the "matched nothing" bucket interpretable.

Without a baseline, a precision number from one-shot template matching says nothing:
the real question is whether one human click beats *any dark blob*, and only a nucleus
detector answers it.

Note there was nothing to reuse from `explore_dataset.ipynb` here. `tighten_box_otsu`
takes a crop *and a known centre* and returns a refined box -- it is a refinement
function, not a detector. And the obvious extension does not work either: global Otsu
over a whole 2 mm^2 H&E ROI separates tissue from white space, not nucleus from
cytoplasm, and returns one enormous connected component. Global Otsu is the right tool
for the tissue mask, which is exactly where it is used below.
"""

from __future__ import annotations

import cv2
import numpy as np
import pandas as pd
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops

from .channels import to_gray_inverted, to_hematoxylin


def tissue_mask(rgb: np.ndarray, close_px: int = 9) -> np.ndarray:
    """Tissue vs. slide background, by global Otsu on inverted grayscale."""
    inv = to_gray_inverted(rgb)
    thr = threshold_otsu(inv)
    mask = (inv > thr).astype(np.uint8)
    if close_px:
        k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (close_px, close_px))
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, k)
        mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, k)
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
    max_detections: int = 20000,
) -> pd.DataFrame:
    """Ranked "any dark nucleus" detector: hematoxylin -> tiled Otsu -> components.

    Components are ranked by mean hematoxylin intensity (darkest chromatin first) so the
    output is a ranked list and `recall@K` is defined for it. A connected-component set
    is otherwise unordered, which would make it incomparable with a scored detector.
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
        df = df.sort_values("score", ascending=False).head(max_detections).reset_index(drop=True)
    df.insert(0, "rank", np.arange(len(df)))
    return df


def random_in_tissue(mask: np.ndarray, n: int, rng=None) -> pd.DataFrame:
    """Uniform random points restricted to tissue -- the floor.

    Restricting to tissue matters: uniform-over-ROI points land in white space, which a
    real detector never does, so lift over that would look impressive and mean nothing.
    """
    rng = np.random.default_rng(0) if rng is None else rng
    ys, xs = np.nonzero(mask)
    if len(ys) == 0:
        return pd.DataFrame(columns=["rank", "cx", "cy", "score"])
    pick = rng.choice(len(ys), size=min(n, len(ys)), replace=False)
    return pd.DataFrame(
        {
            "rank": np.arange(len(pick)),
            "cx": xs[pick].astype(float),
            "cy": ys[pick].astype(float),
            "score": np.linspace(1.0, 0.0, len(pick)),  # arbitrary but stable ordering
        }
    )
