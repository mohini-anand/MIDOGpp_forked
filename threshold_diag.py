"""Core diagnostic function used by every arm of the threshold-sweep analysis."""
import numpy as np
import cv2
from skimage.filters import threshold_multiotsu
from skimage.measure import label, regionprops
from midog_utils import template_match as tm

OTSU_WINDOW = tm.BASE_SIZE  # 51px, pinned explicitly everywhere (subagent fix #1)


def _otsu_threshold(u8):
    """Same value cv2.threshold(..., THRESH_OTSU) computes -- exposed standalone
    so every arm (binary / frac sweep / multiotsu) can share one code path."""
    t, _ = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    return float(t)


def diagnose(patch, threshold_value, min_area=50, max_area_frac=0.85, min_solidity=0.5,
             center_tolerance=0):
    """Full diagnostic version of seed_selection.tighten_box_otsu's body: same
    normalize -> threshold -> label -> centre-pixel -> sanity-check pipeline, but
    returns *why* a candidate was rejected, not just None.

    ``threshold_value`` is applied with a strict '>' (matching cv2.THRESH_BINARY's
    own semantics exactly, verified against cv2.threshold below) -- everything else
    (8-connectivity labeling, centre-pixel lookup, bbox-containment guard, the three
    sanity-check bounds) is identical to the production function.

    center_tolerance is accepted for signature parity but pinned to 0 by every
    caller in this analysis (subagent fix), so _nearest_label_within is not needed.
    """
    if patch.ndim != 2:
        raise ValueError("expected a single-channel patch")
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    binary = (u8 > threshold_value).astype(np.uint8) * 255
    labels_arr = label(binary, connectivity=2)

    cy, cx = patch.shape[0] // 2, patch.shape[1] // 2
    center_label = labels_arr[cy, cx]

    out = dict(center_foreground=bool(center_label != 0), area=None, solidity=None,
               area_frac=None, longer_side=None, bbox=None, box=None,
               reject_reason=None, accepted=False)

    if center_label == 0:
        out["reject_reason"] = "no_foreground_at_click"
        return out

    region = next(p for p in regionprops(labels_arr) if p.label == center_label)
    min_row, min_col, max_row, max_col = region.bbox
    area = int(region.area)
    solidity = float(region.solidity)
    area_frac = area / patch.size
    longer_side = max(max_row - min_row, max_col - min_col)

    out.update(area=area, solidity=solidity, area_frac=area_frac, longer_side=int(longer_side),
               bbox=(int(min_row), int(max_row), int(min_col), int(max_col)))

    # centre_tolerance=0 -> click's own pixel is foreground -> trivially inside its
    # own component's bbox, so the bbox-containment guard never fires here.
    if area < min_area:
        out["reject_reason"] = "min_area"
    elif area > max_area_frac * patch.size:
        out["reject_reason"] = "max_area_frac"
    elif solidity < min_solidity:
        out["reject_reason"] = "min_solidity"
    else:
        out["accepted"] = True
        out["box"] = out["bbox"]

    return out


def diagnose_binary(patch, **kw):
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    t = _otsu_threshold(u8)
    d = diagnose(patch, t, **kw)
    d["threshold_used"] = t
    return d


def diagnose_headroom(patch, frac, **kw):
    """threshold = otsu + frac * (255 - otsu) -- a fraction of the headroom between
    Otsu's own split and the crop's max (always 255 post-normalize)."""
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    otsu_t = _otsu_threshold(u8)
    t = otsu_t + frac * (255.0 - otsu_t)
    d = diagnose(patch, t, **kw)
    d["threshold_used"] = t
    return d


def diagnose_multiotsu(patch, **kw):
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    if len(np.unique(u8)) < 3:
        d = dict(center_foreground=False, area=None, solidity=None, area_frac=None,
                  longer_side=None, bbox=None, box=None, reject_reason="near_flat", accepted=False,
                  threshold_used=None)
        return d
    try:
        thresholds = threshold_multiotsu(u8, classes=3)
    except ValueError:
        return dict(center_foreground=False, area=None, solidity=None, area_frac=None,
                    longer_side=None, bbox=None, box=None, reject_reason="near_flat", accepted=False,
                    threshold_used=None)
    t = float(thresholds[-1])
    d = diagnose(patch, t, **kw)
    d["threshold_used"] = t
    return d
