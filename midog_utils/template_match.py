"""
    Rotation- and scale-augmented template matching over a whole ROI: build an
    augmentation bank from a seed patch, correlate it against the image, and extract
    ranked peaks.
"""

from __future__ import annotations

from dataclasses import dataclass

import cv2
import numpy as np

BASE_SIZE = 51   # annotation box size, forced odd so the click is the centre pixel
PATCH_SIZE = 73  # ceil(BASE_SIZE * sqrt(2)), also odd

_FLOOR = np.float32(-3.0e38)
_SENTINEL_CUT = -1.5e38


def robust_stats(fused: np.ndarray, valid: np.ndarray, stride: int = 8):
    """
        Median and MAD-scale of a response map, over reachable pixels only.

        fused (np.ndarray): the response map.
        valid (np.ndarray): boolean mask of reachable pixels.
        stride (int): subsampling stride.

        median (float): the sample median.
        mad_scale (float): 1.4826 * the sample MAD, a robust std-dev estimate.
    """
    sample = fused[::stride, ::stride][valid[::stride, ::stride]]
    sample = sample[np.isfinite(sample)]
    if sample.size == 0:
        return float("nan"), float("nan")
    med = float(np.median(sample))
    return med, float(1.4826 * np.median(np.abs(sample - med)))


@dataclass(frozen=True)
class Augmentation:
    index: int
    angle: float
    flip: bool
    scale: float
    size: int


def _odd(n: int, minimum: int = 5) -> int:
    n = int(round(n))
    if n % 2 == 0:
        n += 1
    return max(minimum, n)


def read_padded_patch(img: np.ndarray, cx: float, cy: float, patch_size: int = PATCH_SIZE):
    """
        Square patch centred on a point, large enough to survive arbitrary rotation.

        img (np.ndarray): the source image.
        cx (float): centre pixel, x.
        cy (float): centre pixel, y.
        patch_size (int): the patch's side length.

        Returns np.ndarray or None: the patch, or None if too close to the border.
    """
    half = patch_size // 2
    ix, iy = int(round(cx)), int(round(cy))
    h, w = img.shape[:2]
    if ix - half < 0 or iy - half < 0 or ix + half >= w or iy + half >= h:
        return None
    return np.ascontiguousarray(img[iy - half: iy + half + 1, ix - half: ix + half + 1])


def build_augmentations(patch: np.ndarray, base_size: int = BASE_SIZE, scales=(0.6, 0.8, 1.0), n_angles: int = 12, flips=(False, True)):
    """
        Expand one padded patch into the full (flip x angle x scale) template bank.

        patch (np.ndarray): padded source patch, centred on the click.
        base_size (int): template side length before scaling.
        scales (tuple[float]): scale factors applied to each rotated/flipped base.
        n_angles (int): number of rotation angles, evenly spaced over 360 degrees.
        flips (tuple[bool]): which flip states to generate.

        templates (list[np.ndarray]): the augmentation bank.
        metas (list[Augmentation]): metadata paired with ``templates``.
    """
    patch = patch.astype(np.float32, copy=False)
    c = patch.shape[0] // 2
    half_base = base_size // 2
    lo, hi = c - half_base, c + half_base + 1

    templates, metas = [], []
    for flip in flips:
        flipped = np.ascontiguousarray(patch[:, ::-1] if flip else patch)
        for k in range(n_angles):
            angle = 360.0 * k / n_angles
            if angle == 0.0:
                rotated = flipped
            else:
                m = cv2.getRotationMatrix2D((float(c), float(c)), angle, 1.0)
                rotated = cv2.warpAffine(
                    flipped, m, (patch.shape[1], patch.shape[0]),
                    flags=cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE,
                )
            base = np.ascontiguousarray(rotated[lo:hi, lo:hi])
            for scale in scales:
                size = _odd(base_size * scale)
                if size == base_size:
                    tmpl = base
                else:
                    interp = cv2.INTER_AREA if size < base_size else cv2.INTER_CUBIC
                    tmpl = cv2.resize(base, (size, size), interpolation=interp)
                templates.append(np.ascontiguousarray(tmpl, dtype=np.float32))
                metas.append(Augmentation(len(templates) - 1, angle, flip, scale, size))
    return templates, metas


def _robust_z(res: np.ndarray, stride: int = 8) -> np.ndarray:
    """
        Standardise a response map by its own median and MAD, on a strided subsample.

        res (np.ndarray): a single template's response map.
        stride (int): subsampling stride used to estimate median and MAD.

        Returns np.ndarray: the standardised map, or ``res`` unchanged if the sample
            is too small or too uniform to scale.
    """
    sample = res[::stride, ::stride]
    sample = sample[sample > _SENTINEL_CUT]  # drop the sentinel written for zero-variance windows
    if sample.size < 1000:
        return res
    med = np.median(sample)
    mad = np.median(np.abs(sample - med))
    scale = 1.4826 * mad
    if scale <= 1e-6:
        return res
    return ((res - med) / scale).astype(np.float32)


def fused_response(img: np.ndarray, templates, scale_normalize: bool = False, method: int = cv2.TM_CCOEFF):
    """
        Element-wise max of every augmentation's response, in image-centre coordinates.

        img (np.ndarray): the search image.
        templates (list[np.ndarray]): the augmentation bank.
        scale_normalize (bool): put each map on a robust z-scale before fusing.
        method (int): the cv2.TM_* method passed to cv2.matchTemplate.

        fused (np.ndarray): fused score map.
        best (np.ndarray): winning augmentation index per pixel.
        valid (np.ndarray): mask of pixels any template could reach.
    """
    h, w = img.shape[:2]
    img = np.ascontiguousarray(img, dtype=np.float32)

    fused = np.full((h, w), _FLOOR, dtype=np.float32)
    best = np.full((h, w), -1, dtype=np.int16)
    valid = np.zeros((h, w), dtype=bool)

    for i, tmpl in enumerate(templates):
        th, tw = tmpl.shape[:2]
        res = np.asarray(cv2.matchTemplate(img, tmpl, method), dtype=np.float32)
        np.nan_to_num(res, copy=False, nan=float(_FLOOR), posinf=float(_FLOOR), neginf=float(_FLOOR))
        if scale_normalize:
            res = _robust_z(res)

        oy, ox = (th - 1) // 2, (tw - 1) // 2  # map is indexed by template top-left; shift to centre
        rows, cols = res.shape
        view = fused[oy: oy + rows, ox: ox + cols]
        newer = res > view
        np.copyto(view, res, where=newer)
        np.copyto(best[oy: oy + rows, ox: ox + cols], np.int16(i), where=newer)
        valid[oy: oy + rows, ox: ox + cols] = True

    return fused, best, valid


def extract_peaks(fused, valid, min_distance=7, score_threshold=0.5, max_peaks=20000):
    """
        Local maxima of the fused map, returned best-first.

        fused (np.ndarray): fused score map.
        valid (np.ndarray): mask of pixels any template could reach.
        min_distance (int): minimum pixel separation enforced between peaks.
        score_threshold (float): minimum fused score to keep a peak.
        max_peaks (int): cap on the number of peaks returned.

        centers (np.ndarray): (x, y) pixel coordinates, score-descending.
        scores (np.ndarray): fused score at each centre, score-descending.
    """
    k = 2 * int(min_distance) + 1
    dilated = cv2.dilate(fused, np.ones((k, k), np.uint8))
    mask = (fused >= dilated) & valid & (fused >= score_threshold)

    ys, xs = np.nonzero(mask)
    if len(ys) == 0:
        return np.zeros((0, 2)), np.zeros(0, dtype=np.float32)

    scores = fused[ys, xs]
    order = np.lexsort((ys, xs, -scores))[:max_peaks]  # global tie-break so ties are stable across subset sizes
    centers = np.stack([xs[order], ys[order]], axis=1).astype(np.float64)
    return centers, scores[order]


def blank_seed_square(channel: np.ndarray, cx: float, cy: float, base_size: int):
    """
        Blank the seed's own refined-template footprint out of a COPY of a search channel,
        so its own match can never form as a correlation peak.

        channel (np.ndarray): the search channel for the whole ROI; never mutated.
        cx (float): template centre, x, unrounded.
        cy (float): template centre, y, unrounded.
        base_size (int): the odd template side length; the blanked square is exactly
            base_size x base_size, centred on ``(round(cx), round(cy))`` -- the same
            rounding convention ``read_padded_patch`` uses, so the blanked square lines up
            exactly with the footprint the template itself was cut from.

        Returns tuple[np.ndarray, int]: a blanked copy of ``channel``, and the exact pixel
        count blanked (clipped to the image; ``base_size**2`` unless the seed sits within
        half a template of the ROI edge, which does not occur for a readable seed since
        ``read_padded_patch``'s own border check already requires a half-``patch_size``
        margin, larger than half a template).

        The fill value is ``channel.min()`` over the WHOLE, unblanked ROI -- computed before
        this function writes anything, so a caller that reuses ``channel`` (e.g. for
        ``chromatin_od`` ranking) sees only the original, untouched array.
    """
    fill_value = float(channel.min())
    half = base_size // 2
    ix, iy = int(round(cx)), int(round(cy))
    h, w = channel.shape[:2]
    y0, y1 = max(0, iy - half), min(h - 1, iy + half)
    x0, x1 = max(0, ix - half), min(w - 1, ix + half)
    out = channel.copy()
    out[y0:y1 + 1, x0:x1 + 1] = fill_value
    n_blanked = (y1 - y0 + 1) * (x1 - x0 + 1)
    return out, n_blanked


def plant_and_recover(templates, metas, canvas=257, noise_frac=0.25, tolerance=1.0, rng=None, scale_normalize=False, method: int = cv2.TM_CCOEFF):
    """
        Coordinate round-trip gate: plant each augmentation into a noise canvas and check
        the fused map's peak lands within tolerance of the known centre.

        templates (list[np.ndarray]): the augmentation bank.
        metas (list[Augmentation]): metadata paired with ``templates``.
        canvas (int): side length of the synthetic canvas.
        noise_frac (float): background noise std, as a fraction of the template's own std.
        tolerance (float): allowed peak offset from the planted centre, in pixels.
        rng (np.random.Generator or None): defaults to a fixed seed when None.
        scale_normalize (bool): put each map on a robust z-scale before fusing.
        method (int): the cv2.TM_* method passed to cv2.matchTemplate.

        Returns list[tuple]: (Augmentation, dx, dy, score) per augmentation.
        Raises AssertionError when any augmentation's peak lands outside tolerance.
    """
    rng = np.random.default_rng(0) if rng is None else rng
    results = []
    cx = cy = canvas // 2

    for meta, tmpl in zip(metas, templates):
        loc = float(np.mean(tmpl))
        sd = float(np.std(tmpl)) * float(noise_frac)
        img = rng.normal(loc, sd if sd > 0 else 1e-6, size=(canvas, canvas)).astype(np.float32)
        half = tmpl.shape[0] // 2
        img[cy - half: cy + half + 1, cx - half: cx + half + 1] = tmpl

        fused, _, valid = fused_response(img, templates, scale_normalize=scale_normalize, method=method)
        fused = np.where(valid, fused, _FLOOR)
        py, px = np.unravel_index(int(np.argmax(fused)), fused.shape)
        dx, dy = px - cx, py - cy
        if abs(dx) > tolerance or abs(dy) > tolerance:
            raise AssertionError(f"augmentation {meta} recovered at offset ({dx}, {dy}) px, tolerance {tolerance}")
        results.append((meta, dx, dy, float(fused[py, px])))
    return results
