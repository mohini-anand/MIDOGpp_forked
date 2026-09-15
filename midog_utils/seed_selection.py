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
from skimage.filters import threshold_multiotsu
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

        Returns tuple[pd.DataFrame, bool]: (pool, flagged) -- pool is the unanimous tier
        when non-empty, else the 2-of-3 contested tier with flagged=True.
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
    """Annotations far enough from the ROI edge to read a full padded patch, tested on
    the rounded centre."""
    h, w = roi_shape[:2]
    ix = np.rint(df["cx"].to_numpy()).astype(int)
    iy = np.rint(df["cy"].to_numpy()).astype(int)
    ok = (ix >= border) & (ix <= w - 1 - border) & (iy >= border) & (iy <= h - 1 - border)
    return df[ok]


def _nearest_label_within(labels: np.ndarray, cy: int, cx: int, tolerance: int) -> int:
    """The foreground label closest to (cy, cx) within an L-inf tolerance, or 0."""
    h, w = labels.shape
    y0, y1 = max(0, cy - tolerance), min(h, cy + tolerance + 1)
    x0, x1 = max(0, cx - tolerance), min(w, cx + tolerance + 1)
    window = labels[y0:y1, x0:x1]
    ys, xs = np.nonzero(window)
    if len(ys) == 0:
        return 0
    dy, dx = (ys + y0 - cy), (xs + x0 - cx)
    nearest = np.argmin(dy * dy + dx * dx)
    return int(window[ys[nearest], xs[nearest]])


def tighten_box_otsu(patch: np.ndarray, min_area: int = 50, max_area_frac: float = 0.85, min_solidity: float = 0.5, method: str = "binary", center_tolerance: int = 0, headroom_frac: float = None):
    """
        Threshold ``patch`` and return the connected component under (or near) its centre.

        patch (np.ndarray): single-channel patch, "more object -> higher value".
        min_area (int): reject a component smaller than this.
        max_area_frac (float): reject a component larger than this fraction of the patch.
        min_solidity (float): reject a component less convex than this (filled/hull area).
        method (str): "binary" (two-class Otsu), "multiotsu" (3-class, brightest kept), or
            "headroom" (Otsu threshold raised by headroom_frac of the headroom to max).
        center_tolerance (int): widen the centre check to the nearest foreground pixel
            within this L-inf half-width, instead of requiring the exact click pixel.
        headroom_frac (float): required in [0, 1] when method="headroom".

        Returns tuple[int, int, int, int] or None: (y0, y1, x0, x1) half-open, patch-local
        bbox, or None when the centre isn't foreground or the component fails the gate.
    """
    if patch.ndim != 2:
        raise ValueError(
            f"tighten_box_otsu needs a single-channel patch, got shape {patch.shape} -- "
            "pass the structural (gray_inverted/hematoxylin) channel, not channels.to_rgb"
        )
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    if method == "binary":
        _, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    elif method == "multiotsu":
        if len(np.unique(u8)) < 3:
            return None
        try:
            thresholds = threshold_multiotsu(u8, classes=3)
        except ValueError:
            return None
        binary = np.where(u8 > thresholds[-1], np.uint8(255), np.uint8(0))
    elif method == "headroom":
        if headroom_frac is None or not (0 <= headroom_frac <= 1):
            raise ValueError(
                f"headroom_frac must be a float in [0, 1] for method='headroom', "
                f"got {headroom_frac!r}"
            )
        otsu_thresh, _ = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        threshold = otsu_thresh + headroom_frac * (255 - otsu_thresh)
        binary = np.where(u8 > threshold, np.uint8(255), np.uint8(0))
    else:
        raise ValueError(f"method must be 'binary', 'multiotsu', or 'headroom', got {method!r}")
    labels = label(binary, connectivity=2)
    cy, cx = patch.shape[0] // 2, patch.shape[1] // 2
    center_label = labels[cy, cx]
    if center_label == 0 and center_tolerance > 0:
        center_label = _nearest_label_within(labels, cy, cx, center_tolerance)
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


def foreground_filter(df: pd.DataFrame, structural_channel: np.ndarray, otsu_window: int = tm.BASE_SIZE, method: str = "binary", center_tolerance: int = 0, headroom_frac: float = None) -> pd.DataFrame:
    """
        Keep only annotations whose click lands inside (or near) its own Otsu component.

        structural_channel (np.ndarray): single-channel image `tighten_box_otsu` runs on.
        otsu_window (int): window size around each click.
        method, center_tolerance, headroom_frac: passed through to `tighten_box_otsu`.

        Returns pd.DataFrame: the surviving rows of ``df``.
    """
    keep = np.zeros(len(df), dtype=bool)
    for i, (_, row) in enumerate(df.iterrows()):
        patch = tm.read_padded_patch(structural_channel, row["cx"], row["cy"], otsu_window)
        keep[i] = patch is not None and tighten_box_otsu(
            patch, method=method, center_tolerance=center_tolerance, headroom_frac=headroom_frac
        ) is not None
    return df[keep]


@dataclass
class SeedInfo:
    agreement_flagged: bool
    n_agreement_pool: int
    n_after_border: int
    n_after_foreground: int


def pick_seed(gt_mitotic: pd.DataFrame, structural_channel: np.ndarray, rng, border: int, roi_shape, otsu_window: int = tm.BASE_SIZE, method: str = "binary", center_tolerance: int = 0, headroom_frac: float = None, tighten_bbox: bool = True):
    """
        Pick a mitotic seed under the pathologist-agreement and bbox-tightening filters.
        `build_seed` is the entry point for new work; this is kept for callers that only
        need the drawn row.

        gt_mitotic (pd.DataFrame): the image's mitotic ground truth.
        structural_channel (np.ndarray): single-channel image for the foreground filter.
        rng: numpy Generator, drawn from without replacement.
        border (int): edge margin, usually FSConfig.patch_size // 2.
        roi_shape: the ROI's array shape.
        otsu_window, method, center_tolerance, headroom_frac: passed through to
            `tighten_box_otsu`/`foreground_filter`.
        tighten_bbox (bool): when False, skip the foreground filter entirely.

        Returns tuple[pd.Series, SeedInfo]: the drawn row and per-stage pool sizes; raises
        ValueError naming the stage that emptied the pool.
    """
    pool, flagged = agreement_pool(gt_mitotic)
    n_pool = len(pool)
    pool = border_filter(pool, border, roi_shape)
    n_border = len(pool)
    if tighten_bbox:
        pool = foreground_filter(pool, structural_channel, otsu_window, method=method,
                                 center_tolerance=center_tolerance, headroom_frac=headroom_frac)
    n_fg = len(pool)

    info = SeedInfo(flagged, n_pool, n_border, n_fg)
    if n_fg == 0:
        stage = "agreement" if n_pool == 0 else ("border" if n_border == 0 else "foreground")
        raise ValueError(
            f"no seed candidates left (emptied at the {stage} filter); "
            f"agreement_flagged={flagged}, pool sizes: agreement={n_pool}, "
            f"border={n_border}, foreground={n_fg}"
        )
    seed = pool.iloc[int(rng.integers(len(pool)))]
    return seed, info


def tightened_base_size(structural_channel: np.ndarray, cx: float, cy: float, otsu_window: int = tm.BASE_SIZE, minimum: int = 5, method: str = "binary", center_tolerance: int = 0, headroom_frac: float = None):
    """
        The native template size for a seed after bbox tightening, click-centred (size
        only, no recentring). See `tightened_template_box` for the recentred variant
        `build_seed` uses by default.

        Returns int or None: the odd template size, or None if ungated.
    """
    patch = tm.read_padded_patch(structural_channel, cx, cy, otsu_window)
    if patch is None:
        return None
    bbox = tighten_box_otsu(patch, method=method, center_tolerance=center_tolerance,
                            headroom_frac=headroom_frac)
    if bbox is None:
        return None
    y0, y1, x0, x1 = bbox
    size = max(y1 - y0, x1 - x0)
    return _odd(size, minimum=minimum)


def tightened_template_box(structural_channel: np.ndarray, cx: float, cy: float, otsu_window: int = tm.BASE_SIZE, minimum: int = 5, method: str = "binary", center_tolerance: int = 0, headroom_frac: float = None):
    """
        The template size and centre for a seed: both taken from the accepted component.
        D8's production seed/template constructor (`D8_TEMPLATE_ANCHOR.md`).

        The centre is the accepted component's bounding-box pixel centre --
        ``(x0 + x1 - 1) / 2``, not ``(x0 + x1) / 2``, since a half-open bbox covers pixels
        x0..x1-1.

        Returns tuple[int, float, float] or None: (base_size, center_x, center_y), or None
        when the patch can't be read or the gate refuses.

        Caller owns: self-hit/seed-annulus removal must reference the returned centre, not
        (cx, cy); border readability at the returned centre must be re-checked; ground
        truth stays keyed to (cx, cy), never the returned centre.
    """
    patch = tm.read_padded_patch(structural_channel, cx, cy, otsu_window)
    if patch is None:
        return None
    bbox = tighten_box_otsu(patch, method=method, center_tolerance=center_tolerance,
                            headroom_frac=headroom_frac)
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
    """One drawn seed: where its template is cut, and where its ground truth stays.
    ``click_xy`` and ``template_xy`` differ under ``recentred=True``; nothing downstream
    may substitute one for the other."""

    ann_id: int # id of the seleced annotation
    click_xy: tuple # center of the click
    template_xy: tuple # cener of the tightened template
    base_size: int # size of the tightened template
    recentred: bool # whether the template is recentered based on the component
    offset_px: float # offset between the click and template 
    n_retries: int # number of retries
    agreement_flagged: bool # whether the seed is flagged; True for non-unanimous annotations
    n_agreement_pool: int # number of seeds in the agreement pool
    n_after_border: int # number of seeds after the border filter


def _patch_readable(roi_shape, cx: float, cy: float, patch_size: int) -> bool:
    """`template_match.read_padded_patch`'s own bounds predicate, on the rounded centre."""
    half = patch_size // 2
    ix, iy = int(round(cx)), int(round(cy))
    h, w = roi_shape[:2]
    return not (ix - half < 0 or iy - half < 0 or ix + half >= w or iy + half >= h)


def build_seed(gt_mitotic: pd.DataFrame, structural_channel: np.ndarray, rng, roi_shape, patch_size: int = tm.PATCH_SIZE, otsu_window: int = tm.BASE_SIZE, recentre: bool = True, border: int = None, method: str = "binary", center_tolerance: int = 0, headroom_frac: float = None) -> Seed:
    """
        Draw a seed and everything needed to cut its template. Production entry point
        (`D8_TEMPLATE_ANCHOR.md`).

        gt_mitotic (pd.DataFrame): the image's mitotic ground truth, the draw pool.
        structural_channel (np.ndarray): single-channel image for the Otsu gate.
        rng: numpy Generator; drawn without replacement, retrying on a refused candidate.
        roi_shape: the ROI's array shape.
        patch_size (int): full rotation-safe patch size for the border check.
        otsu_window (int): window size for the Otsu gate.
        recentre (bool): True (default, production) centres the template on the accepted
            component's bbox centre; False keeps it on the click and tightens size only.
        border (int): edge margin; defaults to patch_size // 2.
        method, center_tolerance, headroom_frac: passed through to `tighten_box_otsu`.

        Returns Seed. Raises ValueError naming the stage that emptied the pool.

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
        kw = dict(otsu_window=otsu_window, method=method,
                  center_tolerance=center_tolerance, headroom_frac=headroom_frac)
        if recentre:
            got = tightened_template_box(structural_channel, cx, cy, **kw)
            spec = None if got is None else (got[0], got[1], got[2])
        else:
            got = tightened_base_size(structural_channel, cx, cy, **kw)
            spec = None if got is None else (got, cx, cy)
        if spec is not None and _patch_readable(roi_shape, spec[1], spec[2], patch_size):
            base_size, tx, ty = spec
            return Seed(
                ann_id=int(row["ann_id"]), 
                click_xy=(cx, cy),
                template_xy=(tx, ty),
                base_size=int(base_size), 
                recentred=bool(recentre),
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
