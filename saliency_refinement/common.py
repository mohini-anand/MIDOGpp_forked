"""
    Shared pieces for testing the design notes' saliency box refinement ("Mitosis annotation co-pilot pipeline improvement
    notes.pdf", Stage 1b, "Recommended implementation") against this repo's production Otsu gate (D8, D12).

    The saliency rule is implemented as the notes specify it: 3-class Otsu on hematoxylin OD inside the box, score each
    component by darkness * sqrt(area) * gaussian(distance to box centre, sigma = 0.30 * box diagonal), keep the winner and
    every peer at >= 0.45 of its score, optionally dilate the union by 1 um. Every length is set in micrometres.

    The template the repo then cuts is production's: a square of the kept region's bounding-box size, centred on its
    bounding-box pixel centre, cut raw from hematoxylin OD. So saliency changes only which object and what size square;
    the search, blanking, NMS and ranking are production's, unchanged.
"""

from __future__ import annotations

import hashlib
import math
import os
import sys
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from skimage.filters import threshold_multiotsu
from skimage.measure import label, regionprops

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from midog_utils import annotation_frame as af  # noqa: E402
from midog_utils import channels as ch  # noqa: E402
from midog_utils import chromatin as cm  # noqa: E402
from midog_utils import dataset as ds  # noqa: E402
from midog_utils import evaluate as ev  # noqa: E402
from midog_utils import production as prod  # noqa: E402
from midog_utils import seed_selection as ss  # noqa: E402
from midog_utils import template_match as tm  # noqa: E402
from midog_utils.find_and_suppress import FSConfig, find_and_suppress  # noqa: E402

SUBSET_DIRS = {"original_14": ROOT / "images" / "extra_valid", "testing_35": ROOT / "images" / "extra_valid" / "testing_set"}
DB_PATH = ROOT / "databases" / "MIDOG++.json"
RESULTS = ROOT / "results"
OUT_DIR = Path(__file__).resolve().parent / "out"
BUDGETS = (10, 20, 30)
PATCH = tm.PATCH_SIZE          # 73
HALF_PATCH = PATCH // 2        # 36
OD_PAD = prod.OD_WINDOW // 2   # 25

# The notes' parameters (p.8 "Recommended implementation", p.9 parameter sweep, p.2 size floor)
SIGMA_FRAC = 0.30              # gaussian sigma as a fraction of the box diagonal
PEER_FRAC = 0.45               # keep every component scoring at least this fraction of the winner
MIN_AREA_UM2 = 1.0             # component size floor; the notes sweep 0.4-1.6 um^2 and use ~15 px at 0.25 um/px
DILATE_UM = 1.0                # "union(keep) dilated by 1 um"
SEARCH_MARGIN_UM = 1.0         # "let the mask search extend slightly beyond the drawn box"
GUARD_MIN_FRAC = 0.05          # "flag when the kept mask occupies less than about 5% of the box area"


def roi_table() -> pd.DataFrame:
    """
        Every ROI TIFF in the two validation subsets with its domain.

        Returns pd.DataFrame: subset, file_name, path, image_id, domain.
    """
    images, _ = ds.load_annotations(str(DB_PATH))
    rows = [dict(subset=s, file_name=p.name, path=p) for s, d in SUBSET_DIRS.items() for p in sorted(d.glob("*.tiff"))]
    out = pd.DataFrame(rows).merge(images[["file_name", "image_id", "tumor_type"]], on="file_name").rename(columns={"tumor_type": "domain"})
    return out


def load_roi(path):
    """
        One ROI's image, channels, resolution and annotations.

        path (Path): the ROI TIFF.

        Returns dict: rgb, hem, gray, hem_od_padded, mpp, gt (all categories), fn.
    """
    _, annotations = ds.load_annotations(str(DB_PATH))
    rgb = ds.load_roi(path)
    hem = ch.to_channel(rgb, "hematoxylin_od")
    gray = ch.to_channel(rgb, "gray_inverted")
    od_padded = cv2.copyMakeBorder(hem, OD_PAD, OD_PAD, OD_PAD, OD_PAD, borderType=cv2.BORDER_REPLICATE)
    return dict(fn=path.name, rgb=rgb, hem=hem, gray=gray, od_padded=od_padded, mpp=ds.roi_mpp(path), gt=ds.image_annotations(annotations, path.name))


def readable(shape, cx: float, cy: float, size: int = PATCH) -> bool:
    """
        Whether a size-px patch fits at the rounded centre, read_padded_patch's own test.

        shape (tuple): the ROI's shape.
        cx (float): centre x.
        cy (float): centre y.
        size (int): the patch side.

        Returns bool: True when the patch is readable.
    """
    return ss._patch_readable(shape, cx, cy, size)


def gray_component(patch: np.ndarray):
    """
        tighten_box_otsu's accepted component as a mask, re-derived with its own internals and gates.

        patch (np.ndarray): the gate window, single channel.

        Returns tuple[np.ndarray, tuple] or None: the component mask and its (y0, y1, x0, x1) half-open bbox, or None when
            tighten_box_otsu would refuse.
    """
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    _, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    labels = label(binary, connectivity=2)
    c = patch.shape[0] // 2
    lab = labels[c, c]
    if lab == 0:
        return None
    r = next(p for p in regionprops(labels) if p.label == lab)
    y0, x0, y1, x1 = r.bbox
    if not (y0 <= c < y1 and x0 <= c < x1) or r.area < 50 or r.area > 0.85 * patch.size or r.solidity < 0.5:
        return None
    return labels == lab, (int(y0), int(y1), int(x0), int(x1))


def production_spec(gray: np.ndarray, cx: float, cy: float, window: int, shape):
    """
        Production's template for a seed (annotation_frame._template_choice with a precomputed gray channel): the gray
        Otsu gate's component in a window-px window, or the window itself centred on the seed when the gate refuses.

        gray (np.ndarray): the ROI's gray_inverted channel.
        cx (float): seed centre x (the box centre).
        cy (float): seed centre y.
        window (int): the gate window and fallback size, gate_window_side of the seed's box.
        shape (tuple): the ROI's shape.

        Returns dict: base_size, tx, ty, status ("refined", "refused" or "unreadable").
    """
    got = ss.tightened_template_box(gray, cx, cy, otsu_window=window)
    if got is None:
        return dict(base_size=window, tx=cx, ty=cy, status="refused")
    base_size, tx, ty = int(got[0]), float(got[1]), float(got[2])
    if not readable(shape, tx, ty):
        return dict(base_size=window, tx=cx, ty=cy, status="unreadable")
    return dict(base_size=base_size, tx=tx, ty=ty, status="refined")


def box_in_patch(x_tl: float, y_tl: float, width: float, height: float, ix: int, iy: int):
    """
        A frame box's pixels in the coordinates of the 73 px patch centred on (ix, iy).

        x_tl (float): the box's first pixel centre, x (frame convention: pixels x_tl .. x_tl + width).
        y_tl (float): the box's first pixel centre, y.
        width (float): the box's extent, pixel count - 1.
        height (float): the box's extent on y.
        ix (int): the patch centre pixel, x.
        iy (int): the patch centre pixel, y.

        Returns tuple[int, int, int, int]: (y0, y1, x0, x1) half-open, patch-local, clipped to the patch.
    """
    ox, oy = ix - HALF_PATCH, iy - HALF_PATCH
    x0, y0 = int(round(x_tl)) - ox, int(round(y_tl)) - oy
    x1, y1 = x0 + int(round(width)) + 1, y0 + int(round(height)) + 1
    return max(0, y0), min(PATCH, y1), max(0, x0), min(PATCH, x1)


def saliency_mask(hem_patch: np.ndarray, box, mpp: float):
    """
        The notes' saliency rule inside a drawn box, on a 73 px hematoxylin OD patch.

        hem_patch (np.ndarray): 73x73 hematoxylin OD patch centred on the box centre pixel.
        box (tuple): (y0, y1, x0, x1) half-open, patch-local drawn box.
        mpp (float): microns per pixel.

        Returns dict: mask (undilated kept union, 73x73 bool, or None when the box holds no foreground), mask_1um (dilated),
            best (the winning component's mask), n_comp, n_above_floor, n_kept, below_floor (winner chosen with no component
            over the floor), guard_small, guard_edge.
    """
    by0, by1, bx0, bx1 = box
    m = int(round(SEARCH_MARGIN_UM / mpp))
    sy0, sy1, sx0, sx1 = max(0, by0 - m), min(PATCH, by1 + m), max(0, bx0 - m), min(PATCH, bx1 + m)
    region = hem_patch[sy0:sy1, sx0:sx1].astype(np.float64)
    empty = dict(mask=None, mask_1um=None, best=None, n_comp=0, n_above_floor=0, n_kept=0, below_floor=False, guard_small=True, guard_edge=False)
    try:
        t = threshold_multiotsu(region, classes=3)[-1]
    except ValueError:
        return empty
    fg = region > t
    lab = label(fg, connectivity=2)
    props = regionprops(lab, intensity_image=region)
    if not props:
        return empty
    floor_px = MIN_AREA_UM2 / (mpp * mpp)
    above = [p for p in props if p.area >= floor_px]
    pool = above if above else props
    cy_box, cx_box = (by0 + by1 - 1) / 2.0 - sy0, (bx0 + bx1 - 1) / 2.0 - sx0
    sigma = SIGMA_FRAC * math.hypot(by1 - by0, bx1 - bx0)
    scores = np.array([p.intensity_mean * math.sqrt(p.area) * math.exp(-((p.centroid[0] - cy_box) ** 2 + (p.centroid[1] - cx_box) ** 2) / (2 * sigma * sigma)) for p in pool])
    best = int(np.argmax(scores))
    kept = [p for p, s in zip(pool, scores) if s >= PEER_FRAC * scores[best]]
    local = np.isin(lab, [p.label for p in kept])
    mask = np.zeros((PATCH, PATCH), bool)
    mask[sy0:sy1, sx0:sx1] = local
    best_mask = np.zeros((PATCH, PATCH), bool)
    best_mask[sy0:sy1, sx0:sx1] = lab == pool[best].label
    r = int(round(DILATE_UM / mpp))
    se = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1))
    mask_1um = cv2.dilate(mask.astype(np.uint8), se).astype(bool)
    inside = np.zeros((PATCH, PATCH), bool)
    inside[by0 + 1:by1 - 1, bx0 + 1:bx1 - 1] = True  # the box minus its own border pixels
    box_area = (by1 - by0) * (bx1 - bx0)
    return dict(mask=mask, mask_1um=mask_1um, best=best_mask, n_comp=len(props), n_above_floor=len(above), n_kept=len(kept), below_floor=not above, guard_small=bool(mask.sum() < GUARD_MIN_FRAC * box_area), guard_edge=bool((mask & ~inside).any()))


def mask_bbox(mask: np.ndarray):
    """
        The half-open bounding rectangle of a non-empty boolean mask.

        mask (np.ndarray): boolean mask.

        Returns tuple[int, int, int, int]: (y0, y1, x0, x1).
    """
    rows, cols = np.flatnonzero(mask.any(axis=1)), np.flatnonzero(mask.any(axis=0))
    return int(rows[0]), int(rows[-1]) + 1, int(cols[0]), int(cols[-1]) + 1


def spec_from_bbox(bbox, ix: int, iy: int):
    """
        Production's template spec for a patch-local bbox (seed_selection.tightened_template_box's own arithmetic), side
        capped at the 73 px patch.

        bbox (tuple): (y0, y1, x0, x1) half-open, patch-local.
        ix (int): the patch centre pixel, x.
        iy (int): the patch centre pixel, y.

        Returns tuple[int, float, float]: (base_size, tx, ty) in image coordinates.
    """
    y0, y1, x0, x1 = bbox
    base_size = min(ss._odd(max(y1 - y0, x1 - x0), minimum=5), PATCH)
    return base_size, ix - HALF_PATCH + (x0 + x1 - 1) / 2.0, iy - HALF_PATCH + (y0 + y1 - 1) / 2.0


def saliency_spec(hem: np.ndarray, x_tl: float, y_tl: float, width: float, height: float, mpp: float, window: int, *, dilate: bool):
    """
        The saliency rule's template for a frame box, with production's fallback (the box's window centred on its centre)
        when the box holds no foreground or the template is unreadable.

        hem (np.ndarray): the ROI's hematoxylin OD channel.
        x_tl (float): frame box first pixel x.
        y_tl (float): frame box first pixel y.
        width (float): frame box extent x.
        height (float): frame box extent y.
        mpp (float): microns per pixel.
        window (int): gate_window_side of the box, the fallback template size.
        dilate (bool): use the 1 um dilated union (the notes' literal recipe) rather than the kept union itself.

        Returns tuple[dict, dict]: the spec (base_size, tx, ty, status) and the saliency_mask output with the patch origin.
    """
    cx, cy = af.box_center(x_tl, width), af.box_center(y_tl, height)
    ix, iy = int(round(cx)), int(round(cy))
    patch = tm.read_padded_patch(hem, cx, cy, PATCH)
    assert patch is not None, f"({cx}, {cy}): 73 px patch not readable"
    box = box_in_patch(x_tl, y_tl, width, height, ix, iy)
    sal = saliency_mask(patch, box, mpp)
    sal.update(ix=ix, iy=iy, box=box)
    key = "mask_1um" if dilate else "mask"
    if sal[key] is None:
        return dict(base_size=window, tx=cx, ty=cy, status="no_foreground"), sal
    base_size, tx, ty = spec_from_bbox(mask_bbox(sal[key]), ix, iy)
    if not readable(hem.shape, tx, ty):
        return dict(base_size=window, tx=cx, ty=cy, status="unreadable"), sal
    return dict(base_size=base_size, tx=tx, ty=ty, status="below_floor" if sal["below_floor"] else "refined"), sal


def template_sha1(hem: np.ndarray, base_size: int, tx: float, ty: float) -> str:
    """
        SHA-1 of the exact template find_and_suppress cuts, the recorded runs' tpl_sha1.

        hem (np.ndarray): the ROI's hematoxylin OD channel.
        base_size (int): template side.
        tx (float): template centre x.
        ty (float): template centre y.

        Returns str: the hex digest.
    """
    templates, _ = tm.build_augmentations(tm.read_padded_patch(hem, tx, ty, PATCH), base_size, (1.0,), 1, (False,))
    return hashlib.sha1(np.ascontiguousarray(templates[0], dtype=np.float32).tobytes()).hexdigest()


def search(roi: dict, base_size: int, tx: float, ty: float) -> tuple[pd.DataFrame, dict]:
    """
        production.run_production_pipeline's search on a precomputed hematoxylin channel, with od always computed so both
        ranking axes rank the same detections.

        roi (dict): a load_roi output.
        base_size (int): template side.
        tx (float): template centre x.
        ty (float): template centre y.

        Returns tuple[pd.DataFrame, dict]: detections in find_and_suppress order with an od column, and its info.
    """
    nms_radius = ev.radius_px(roi["mpp"])
    cfg = FSConfig(base_size=base_size, scales=(1.0,), n_angles=1, flips=(False,), peak_min_distance=prod.PEAK_MIN_DISTANCE, max_peaks=prod.MAX_PEAKS, nms_radius=nms_radius, deep_floor_z=prod.DEEP_FLOOR_Z, border_pad=True, tm_method=prod.TM_METHOD)
    det, info = find_and_suppress(roi["hem"], (tx, ty), cfg, nms_radius=nms_radius)
    shifted = det.assign(cx=det["cx"] + OD_PAD, cy=det["cy"] + OD_PAD)
    det = det.assign(od=cm.score_detections(shifted, roi["od_padded"], window=prod.OD_WINDOW)["od"].to_numpy())
    return det, info


def rank_by(det: pd.DataFrame, rank_key: str) -> pd.DataFrame:
    """
        Production's sort for one axis.

        det (pd.DataFrame): search output.
        rank_key (str): "tm_score" or "chromatin_od".

        Returns pd.DataFrame: ranked detections, rank reassigned.
    """
    out = det.sort_values(prod.AXES[rank_key], ascending=False, na_position="last", kind="mergesort").reset_index(drop=True)
    return out.assign(rank=np.arange(len(out)))


def score(det_ranked: pd.DataFrame, gt: pd.DataFrame, seed_ann_id: int, click_xy, mpp: float) -> dict:
    """
        TP at each budget against all annotations but the seed's own, plus the seed's own figure leaking into the top 30.

        det_ranked (pd.DataFrame): ranked detections.
        gt (pd.DataFrame): the ROI's annotations, all categories.
        seed_ann_id (int): the seed annotation, excluded from ground truth.
        click_xy (tuple): the seed's annotated centre, for the leak count.
        mpp (float): microns per pixel.

        Returns dict: tp_at_K for each budget and n_near_click (top-30 detections within one match radius of the click).
    """
    radius = ev.radius_px(mpp)
    gt_eval = gt[gt["ann_id"] != seed_ann_id].reset_index(drop=True)
    out, _ = ev.bucket_detections(det_ranked, gt_eval, radius)
    hit = (out["bucket"] == ev.HUMAN_CORRECT_LABEL).to_numpy()
    top = det_ranked.head(30)
    near = np.hypot(top["cx"].to_numpy() - click_xy[0], top["cy"].to_numpy() - click_xy[1]) <= radius
    return {**{f"tp_at_{k}": int(hit[:k].sum()) for k in BUDGETS}, "n_near_click": int(near.sum())}


def midog_box(ann: pd.Series):
    """
        A MIDOG annotation's frame box (x_tl, y_tl, width, height), as human_seed_frame stores it.

        ann (pd.Series): one load_annotations row.

        Returns tuple[float, float, float, float]: the box.
    """
    return float(ann["cx"]) - float(ann["w"]) / 2, float(ann["cy"]) - float(ann["h"]) / 2, float(ann["w"]), float(ann["h"])
