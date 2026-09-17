"""
    Seed selection under pathologist agreement and bbox-tightening filters. `build_seed`
    is the production entry point: draws a mitotic annotation, gates it against its own
    Otsu-thresholded connected component, and returns the template's size and centre.
"""

from __future__ import annotations

from dataclasses import dataclass

import cv2
import numpy as np
import pandas as pd
from skimage.measure import label, regionprops

from . import template_match as tm


def _odd(n: int, minimum: int = 5) -> int:
    n = int(round(n))
    if n % 2 == 0:
        n += 1
    return max(minimum, n)


def agreement_pool(gt_mitotic: pd.DataFrame):
    """
        Split an image's mitotic annotations into the agreement tier to draw from.

        gt_mitotic (pd.DataFrame): the image's mitotic ground truth.

        pool (pd.DataFrame): the unanimous tier when non-empty, else the 2-of-3
            contested tier.
        flagged (bool): True when ``pool`` is the contested tier.
    """
    unanimous = gt_mitotic[gt_mitotic["n_mitotic_votes"] == gt_mitotic["n_votes"]]
    if len(unanimous):
        return unanimous, False
    contested = gt_mitotic[
        (gt_mitotic["n_votes"] > 0)
        & (gt_mitotic["n_mitotic_votes"] / gt_mitotic["n_votes"] >= 2.0 / 3.0)
        & (gt_mitotic["n_mitotic_votes"] < gt_mitotic["n_votes"])
    ]
    return contested, True


def border_filter(df: pd.DataFrame, border: int, roi_shape) -> pd.DataFrame:
    """
        Annotations far enough from the ROI edge to read a full padded patch, tested on
        the rounded centre.

        df (pd.DataFrame): candidate annotations, with cx/cy columns.
        border (int): margin in pixels from every ROI edge.
        roi_shape (tuple): the ROI's array shape.

        Returns pd.DataFrame: rows whose rounded centre clears the margin on all sides.
    """
    h, w = roi_shape[:2]
    ix = np.rint(df["cx"].to_numpy()).astype(int)
    iy = np.rint(df["cy"].to_numpy()).astype(int)
    ok = (ix >= border) & (ix <= w - 1 - border) & (iy >= border) & (iy <= h - 1 - border)
    return df[ok]


def tighten_box_otsu(patch: np.ndarray, min_area: int = 50, max_area_frac: float = 0.85, min_solidity: float = 0.5):
    """
        Two-class Otsu threshold ``patch`` and return the connected component under its centre.

        patch (np.ndarray): single-channel patch, "more object -> higher value".
        min_area (int): reject a component smaller than this.
        max_area_frac (float): reject a component larger than this fraction of the patch.
        min_solidity (float): reject a component less convex than this (filled/hull area).

        Returns tuple[int, int, int, int] or None: (y0, y1, x0, x1) half-open, patch-local
        bbox, or None when the centre isn't foreground or the component fails the gate.
    """
    if patch.ndim != 2:
        raise ValueError(
            f"tighten_box_otsu needs a single-channel patch, got shape {patch.shape} -- "
            "an RGB array was passed instead of a single structural channel (e.g. gray_inverted)"
        )
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    _, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    labels = label(binary, connectivity=2)
    cy, cx = patch.shape[0] // 2, patch.shape[1] // 2
    center_label = labels[cy, cx]
    if center_label == 0:
        return None

    region = next(p for p in regionprops(labels) if p.label == center_label)
    min_row, min_col, max_row, max_col = region.bbox
    if not (min_row <= cy < max_row and min_col <= cx < max_col):
        return None
    if region.area < min_area:
        return None
    if region.area > max_area_frac * patch.size:
        return None
    if region.solidity < min_solidity:
        return None

    y0, x0, y1, x1 = region.bbox  # skimage bbox order is (min_row, min_col, max_row, max_col)
    return int(y0), int(y1), int(x0), int(x1)


def tightened_template_box(structural_channel: np.ndarray, cx: float, cy: float, otsu_window: int = tm.BASE_SIZE, minimum: int = 5):
    """
        The template size and centre for a seed: both taken from the accepted component.
        D8's production seed/template constructor (`D8_TEMPLATE_ANCHOR.md`).

        structural_channel (np.ndarray): single-channel image for the Otsu gate.
        cx (float): click centre, x.
        cy (float): click centre, y.
        otsu_window (int): window size for the Otsu gate.
        minimum (int): smallest odd size to return.

        base_size (int): the odd template size.
        center_x (float): accepted component's bbox pixel centre, x -- ``(x0 + x1 - 1) / 2``,
            not ``(x0 + x1) / 2``, since a half-open bbox covers pixels x0..x1-1.
        center_y (float): accepted component's bbox pixel centre, y (same convention).
        Returns None when the patch can't be read or the gate refuses.

        Caller owns: self-hit/seed-annulus removal must reference the returned centre, not
        (cx, cy); border readability at the returned centre must be re-checked; ground
        truth stays keyed to (cx, cy), never the returned centre.
    """
    patch = tm.read_padded_patch(structural_channel, cx, cy, otsu_window)
    if patch is None:
        return None
    bbox = tighten_box_otsu(patch)
    if bbox is None:
        return None
    y0, y1, x0, x1 = bbox
    base_size = _odd(max(y1 - y0, x1 - x0), minimum=minimum)
    half = otsu_window // 2
    ix, iy = int(round(cx)), int(round(cy))  # matches read_padded_patch's own rounding
    # convert the patch's origin to the recentered templates center, in image coordinates
    center_x = ix - half + (x0 + x1 - 1) / 2.0
    center_y = iy - half + (y0 + y1 - 1) / 2.0
    return base_size, center_x, center_y


@dataclass(frozen=True)
class Seed:
    """
        One drawn seed: where its template is cut, and where its ground truth stays.
        ``click_xy`` and ``template_xy`` can differ, since the template is recentred on its
        accepted component; nothing downstream may substitute one for the other.
    """

    ann_id: int # id of the seleced annotation
    click_xy: tuple # center of the click
    template_xy: tuple # cener of the tightened template
    base_size: int # size of the tightened template
    recentred: bool # always True: the template is recentred on the accepted component
    offset_px: float # offset between the click and template 
    n_retries: int # number of retries
    agreement_flagged: bool # whether the seed is flagged; True for non-unanimous annotations
    n_agreement_pool: int # number of seeds in the agreement pool
    n_after_border: int # number of seeds after the border filter


def _patch_readable(roi_shape, cx: float, cy: float, patch_size: int) -> bool:
    """
        `template_match.read_padded_patch`'s own bounds predicate, on the rounded centre.

        roi_shape (tuple): the ROI's array shape.
        cx (float): centre pixel, x.
        cy (float): centre pixel, y.
        patch_size (int): the patch's side length.

        Returns bool: True when the patch fits inside the ROI.
    """
    half = patch_size // 2
    ix, iy = int(round(cx)), int(round(cy))
    h, w = roi_shape[:2]
    return not (ix - half < 0 or iy - half < 0 or ix + half >= w or iy + half >= h)


def build_seed(gt_mitotic: pd.DataFrame, structural_channel: np.ndarray, rng, roi_shape, patch_size: int = tm.PATCH_SIZE, otsu_window: int = tm.BASE_SIZE, border: int = None) -> Seed:
    """
        Draw a seed and everything needed to cut its template, centred on the accepted
        component's bbox centre. Production entry point (`D8_TEMPLATE_ANCHOR.md`).

        gt_mitotic (pd.DataFrame): the image's mitotic ground truth, the draw pool.
        structural_channel (np.ndarray): single-channel image for the Otsu gate.
        rng (np.random.Generator): drawn without replacement, retrying on a refused candidate.
        roi_shape (tuple): the ROI's array shape.
        patch_size (int): full rotation-safe patch size for the border check.
        otsu_window (int): window size for the Otsu gate.
        border (int): edge margin; defaults to patch_size // 2.

        Returns Seed: the drawn seed and everything needed to cut its template.
        Raises ValueError naming the stage that emptied the pool.

        Caller still owns: exclude seed.ann_id from the evaluation set and compute every
        match radius against seed.click_xy, never seed.template_xy.
    """
    border = patch_size // 2 if border is None else border
    pool, flagged = agreement_pool(gt_mitotic)
    n_pool = len(pool)
    pool = border_filter(pool, border, roi_shape)
    n_border = len(pool)

    working, retries = pool.copy(), 0
    while len(working) > 0:
        idx = int(rng.integers(len(working)))
        row = working.iloc[idx]
        cx, cy = float(row["cx"]), float(row["cy"])
        got = tightened_template_box(structural_channel, cx, cy, otsu_window=otsu_window)
        spec = None if got is None else (got[0], got[1], got[2])
        if spec is not None and _patch_readable(roi_shape, spec[1], spec[2], patch_size):
            base_size, tx, ty = spec
            return Seed(
                ann_id=int(row["ann_id"]), 
                click_xy=(cx, cy),
                template_xy=(tx, ty),
                base_size=int(base_size), 
                recentred=True,
                offset_px=float(np.hypot(tx - cx, ty - cy)),
                n_retries=retries,
                agreement_flagged=bool(flagged),
                n_agreement_pool=n_pool,
                n_after_border=n_border
                )
        working = working.drop(working.index[idx])
        retries += 1
    # if the entire pool is exhausted, record why that happened and raise an error
    stage = "agreement" if n_pool == 0 else ("border" if n_border == 0 else "gate")
    raise ValueError(
        f"no seed candidates left (emptied at the {stage} stage); agreement_flagged={flagged}, "
        f"pool sizes: agreement={n_pool}, border={n_border}, refused={retries}"
    )
