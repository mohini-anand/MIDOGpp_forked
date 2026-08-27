"""Rotation- and scale-augmented normalised cross-correlation over a whole ROI.

Three deliberate departures from `bbox tuning code reference/bbox_tuning.py`:

1. **The rotations are real.** `shortlist_augmentations` (bbox_tuning.py:323-324)
   overwrites its resized+rotated template with the un-rotated one on the very next
   line, so all four rotations become byte-identical, the cosine-similarity filter
   collapses them to one, and the rotation search silently does nothing. That is fatal
   here -- mitotic figures have no canonical orientation.

2. **Pad, then rotate, then centre-crop.** A 73x73 patch is read from the ROI, rotated,
   then cropped back to 51x51. No replicated border pixels enter the correlation and no
   box content is thrown away. The bound that makes that true: the far corner of the
   inner 51x51 sits 25*sqrt(2) = 35.36 px from the centre, and the patch reaches 36 px,
   so every source pixel an arbitrary rotation needs is real. (73 is chosen as
   ceil(51*sqrt(2)), which is a looser argument for the same conclusion -- it bounds the
   patch's own diagonal rather than the inner square's.)

3. **One fused response map.** The per-augmentation maps are combined by element-wise
   maximum in a common *centre* coordinate frame and peaks are picked once, instead of
   concatenating per-augmentation detections and suppressing afterwards. Cheaper, and
   it removes the duplicate-peak class entirely.

All template sizes are forced odd so a template's centre is an exact pixel and the
map-to-image offset is `(size - 1) // 2` with no half-pixel rounding.
"""

from __future__ import annotations

from dataclasses import dataclass

import cv2
import numpy as np

BASE_SIZE = 51   # the 50x50 annotation box, made odd so the click is the centre pixel
PATCH_SIZE = 73  # ceil(BASE_SIZE * sqrt(2)), also odd

_FLOOR = np.float32(-2.0)  # below TM_CCOEFF_NORMED's [-1, 1] range


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
    """Square patch centred on a point, large enough to survive arbitrary rotation.

    Returns ``None`` when the point sits too close to the ROI border for a full patch.
    Callers should log and skip those rather than pad, since padding would feed
    fabricated pixels into the template.
    """
    half = patch_size // 2
    ix, iy = int(round(cx)), int(round(cy))
    h, w = img.shape[:2]
    if ix - half < 0 or iy - half < 0 or ix + half >= w or iy + half >= h:
        return None
    return np.ascontiguousarray(img[iy - half: iy + half + 1, ix - half: ix + half + 1])


def build_augmentations(
    patch: np.ndarray,
    base_size: int = BASE_SIZE,
    scales=(0.6, 0.8, 1.0),
    n_angles: int = 12,
    flips=(False, True),
):
    """Expand one padded patch into the full (flip x angle x scale) template bank."""
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
    """Standardise a response map by its own median and MAD.

    `TM_CCOEFF_NORMED` normalises by the template's pixel count, so a smaller template
    has a higher-variance null distribution and therefore a systematically higher
    *maximum* by chance alone. A raw element-wise max across template sizes is
    consequently biased toward the smallest scale. Putting every map on a common z-scale
    first removes that bias.

    Statistics are taken on a strided subsample -- 39 megapixels is far more than needed
    for a median and a MAD.
    """
    sample = res[::stride, ::stride]
    sample = sample[sample > -1.5]  # drop the NaN sentinel from zero-variance windows
    if sample.size < 1000:
        return res
    med = np.median(sample)
    mad = np.median(np.abs(sample - med))
    scale = 1.4826 * mad
    if scale <= 1e-6:
        return res
    return ((res - med) / scale).astype(np.float32)


def fused_response(img: np.ndarray, templates, scale_normalize: bool = False):
    """Element-wise max of every augmentation's response, in image-centre coordinates.

    ``cv2.matchTemplate`` returns a map indexed by template *top-left*; entry ``(r, c)``
    of a ``t``-sized template's map describes a window centred at
    ``(c + (t-1)//2, r + (t-1)//2)``. Each map is therefore shifted by half of *its own*
    size before the maximum is taken -- getting this wrong shifts every detection by
    10-25 px, which at a 30 px match radius looks exactly like "template matching does
    not work on H&E".

    ``scale_normalize`` puts each map on a robust z-scale before fusing, which corrects
    the small-template bias described in ``_robust_z``. It is **off by default** so that
    the fused value stays an interpretable correlation coefficient; with it on, scores
    are z-scores and any threshold must be reinterpreted accordingly.

    Returns ``(fused, best_aug, valid)``: the fused score map, the index of the winning
    augmentation per pixel, and a mask of pixels any template could reach.
    """
    h, w = img.shape[:2]
    img = np.ascontiguousarray(img, dtype=np.float32)

    fused = np.full((h, w), _FLOOR, dtype=np.float32)
    best = np.full((h, w), -1, dtype=np.int16)
    valid = np.zeros((h, w), dtype=bool)

    for i, tmpl in enumerate(templates):
        th, tw = tmpl.shape
        res = cv2.matchTemplate(img, tmpl, cv2.TM_CCOEFF_NORMED)
        # Uniform windows (saturated white background) give zero variance and NaN here.
        np.nan_to_num(res, copy=False, nan=-2.0, posinf=-2.0, neginf=-2.0)
        if scale_normalize:
            res = _robust_z(res)

        oy, ox = (th - 1) // 2, (tw - 1) // 2
        rows, cols = res.shape
        view = fused[oy: oy + rows, ox: ox + cols]
        newer = res > view
        np.copyto(view, res, where=newer)
        np.copyto(best[oy: oy + rows, ox: ox + cols], np.int16(i), where=newer)
        valid[oy: oy + rows, ox: ox + cols] = True

    return fused, best, valid


def extract_peaks(fused, valid, min_distance=7, score_threshold=0.25, max_peaks=20000):
    """Local maxima of the fused map, returned best-first.

    A grey dilation is used rather than `skimage.feature.peak_local_max` for speed on a
    39-megapixel map. This only thins the candidate list; the real de-duplication is the
    score-ordered distance NMS that follows.
    """
    k = 2 * int(min_distance) + 1
    dilated = cv2.dilate(fused, np.ones((k, k), np.uint8))
    mask = (fused >= dilated) & valid & (fused >= score_threshold)

    ys, xs = np.nonzero(mask)
    if len(ys) == 0:
        return np.zeros((0, 2)), np.zeros(0, dtype=np.float32)

    scores = fused[ys, xs]
    order = np.argsort(scores)[::-1][:max_peaks]
    centers = np.stack([xs[order], ys[order]], axis=1).astype(np.float64)
    return centers, scores[order]


def plant_and_recover(templates, metas, canvas=257, noise=8.0, tolerance=1.0, rng=None,
                      scale_normalize=False):
    """Coordinate round-trip gate -- run this before anything touches a real image.

    Each augmentation is planted into a noise canvas at a known centre, the *whole*
    bank is then matched against that canvas, and the peak of the **fused** map must
    land within ``tolerance`` px of where it was planted.

    Asserting on the fused map rather than per-augmentation maps is the point: a 31 px
    template yields a (H-30, W-30) map and a 51 px template a (H-50, W-50) map, and each
    must be offset by half of its own size before they are combined. Per-augmentation
    recovery can pass while the fusion is silently misaligned.

    Returns a list of ``(Augmentation, dx, dy, score)``; raises on the first failure.
    """
    rng = np.random.default_rng(0) if rng is None else rng
    results = []
    cx = cy = canvas // 2

    for meta, tmpl in zip(metas, templates):
        img = rng.normal(120.0, noise, size=(canvas, canvas)).astype(np.float32)
        half = tmpl.shape[0] // 2
        img[cy - half: cy + half + 1, cx - half: cx + half + 1] = tmpl

        fused, _, valid = fused_response(img, templates, scale_normalize=scale_normalize)
        fused = np.where(valid, fused, _FLOOR)
        py, px = np.unravel_index(int(np.argmax(fused)), fused.shape)
        dx, dy = px - cx, py - cy
        if abs(dx) > tolerance or abs(dy) > tolerance:
            raise AssertionError(
                f"augmentation {meta} recovered at offset ({dx}, {dy}) px, "
                f"tolerance {tolerance}"
            )
        results.append((meta, dx, dy, float(fused[py, px])))
    return results
