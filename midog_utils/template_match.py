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

# Below every score any `cv2.TM_*` method can produce on any channel this repo uses, and far
# enough above float32's -3.4e38 minimum that nothing overflows. It replaces an earlier -2.0,
# which was only below range for TM_CCOEFF_NORMED: TM_CCORR on hematoxylin OD spans
# [0.34, 10.7] and sign-flipped TM_SQDIFF [-14.7, -1.9e-6], where -2.0 is a *valid* score, and
# on a 0-255 channel both reach ~1.7e8. Pixels that no template could reach are identified by
# the `valid` mask, never by comparing against this value.
_FLOOR = np.float32(-3.0e38)
_SENTINEL_CUT = -1.5e38  # anything below this is the sentinel, not a score

# `higher is better` for every method, so peak extraction, NMS and ranking are unchanged.
# TM_SQDIFF and TM_SQDIFF_NORMED are distances; everything else is a similarity.
_SIGN = {cv2.TM_SQDIFF: -1.0, cv2.TM_SQDIFF_NORMED: -1.0}

METHODS = {
    "sqdiff": cv2.TM_SQDIFF,
    "sqdiff_normed": cv2.TM_SQDIFF_NORMED,
    "ccorr": cv2.TM_CCORR,
    "ccorr_normed": cv2.TM_CCORR_NORMED,
    "ccoeff": cv2.TM_CCOEFF,
    "ccoeff_normed": cv2.TM_CCOEFF_NORMED,
}


def method_sign(method: int) -> float:
    """+1.0 for similarity methods, -1.0 for the two distance methods."""
    return _SIGN.get(int(method), 1.0)


def robust_stats(fused: np.ndarray, valid: np.ndarray, stride: int = 8):
    """Median and MAD-scale of a response map, over reachable pixels only.

    The `valid` mask replaces the older ``sample > -1.5`` test, which was a
    TM_CCOEFF_NORMED-specific way of dropping the unreachable border and the NaN sentinel.
    That test silently deletes real detections under any other method -- every score
    sign-flipped TM_SQDIFF produces on hematoxylin OD is below -1.5.

    Returns ``(median, 1.4826 * MAD)``. The pair defines the per-map z scale that makes one
    extraction floor mean the same search depth for all six methods; no raw-score threshold
    can, since the methods' ranges differ by eight orders of magnitude.
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
    # Drop the sentinel written for zero-variance (NaN) windows. Compared against half the
    # sentinel rather than a fixed -1.5 so this stays correct for methods whose genuine
    # scores are large and negative, e.g. sign-flipped TM_SQDIFF.
    sample = sample[sample > _SENTINEL_CUT]
    if sample.size < 1000:
        return res
    med = np.median(sample)
    mad = np.median(np.abs(sample - med))
    scale = 1.4826 * mad
    if scale <= 1e-6:
        return res
    return ((res - med) / scale).astype(np.float32)


def fused_response(img: np.ndarray, templates, scale_normalize: bool = False,
                   method: int = cv2.TM_CCOEFF_NORMED):
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

    ``method`` is any ``cv2.TM_*``. The two distance methods (``TM_SQDIFF``,
    ``TM_SQDIFF_NORMED``) are negated on the way in, so "higher is better" holds for every
    method and peak extraction, NMS, ranking and evaluation are all unchanged. Two things
    that matter when leaving the default:

    * **Only rotations and flips may be fused across for the unnormalised methods.**
      ``TM_CCORR``, ``TM_CCOEFF`` and ``TM_SQDIFF`` scale with the template's pixel count, so
      an element-wise max over a bank of *different sizes* is decided by size rather than by
      fit. Rotations and flips at one scale all share a size, so the max stays honest; a
      multi-scale bank needs ``scale_normalize=True`` or it is meaningless. (The normalised
      three have a milder version of the same problem -- see ``_robust_z``.)
    * **Scores are not comparable between methods**, only within one map: ``TM_CCORR`` on
      hematoxylin OD spans [0.34, 10.7] and negated ``TM_SQDIFF`` [-14.7, -1.9e-6]. Use
      `robust_stats` to put a threshold in per-map z units, which is monotone and so changes
      no ranking.

    Returns ``(fused, best_aug, valid)``: the fused score map, the index of the winning
    augmentation per pixel, and a mask of pixels any template could reach.
    """
    h, w = img.shape[:2]
    img = np.ascontiguousarray(img, dtype=np.float32)

    fused = np.full((h, w), _FLOOR, dtype=np.float32)
    best = np.full((h, w), -1, dtype=np.int16)
    valid = np.zeros((h, w), dtype=bool)

    sign = np.float32(method_sign(method))
    for i, tmpl in enumerate(templates):
        th, tw = tmpl.shape[:2]  # [:2]: also correct for a 3-channel (RGB) template
        res = np.asarray(cv2.matchTemplate(img, tmpl, method), dtype=np.float32)
        if sign < 0:
            res *= sign  # a distance becomes a similarity; nothing downstream changes
        # Uniform windows (saturated white background) give zero variance, which is NaN for
        # the three _NORMED methods. The unnormalised three cannot produce it.
        np.nan_to_num(res, copy=False, nan=float(_FLOOR),
                      posinf=float(_FLOOR), neginf=float(_FLOOR))
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


def extract_peaks(fused, valid, min_distance=7, score_threshold=0.5, max_peaks=20000):
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
    # Ordered by a *global* key -- score descending, then x, then y -- rather than by
    # `np.argsort(scores)[::-1]`, whose default quicksort is unstable and therefore breaks
    # exact ties differently depending on how many peaks were passed in. That is not
    # cosmetic here: `premise_test`-style economy re-uses one deep pool for every z level,
    # and the argument that filtering the pool at t equals re-extracting at t holds only if
    # two peaks with byte-identical scores are always ordered the same way. Measured on
    # 350.tiff it did not: two adjacent peaks tied at 0.5803979 swapped, so 2 of 8747
    # coordinates disagreed between the two routes. A global key restricts to any subset
    # unchanged, which makes the shortcut exact rather than exact-up-to-ties.
    order = np.lexsort((ys, xs, -scores))[:max_peaks]
    centers = np.stack([xs[order], ys[order]], axis=1).astype(np.float64)
    return centers, scores[order]


def plant_and_recover(templates, metas, canvas=257, noise_frac=0.25, tolerance=1.0, rng=None,
                      scale_normalize=False, method: int = cv2.TM_CCOEFF_NORMED):
    """Coordinate round-trip gate -- run this before anything touches a real image.

    Each augmentation is planted into a noise canvas at a known centre, the *whole*
    bank is then matched against that canvas, and the peak of the **fused** map must
    land within ``tolerance`` px of where it was planted.

    Asserting on the fused map rather than per-augmentation maps is the point: a 31 px
    template yields a (H-30, W-30) map and a 51 px template a (H-50, W-50) map, and each
    must be offset by half of its own size before they are combined. Per-augmentation
    recovery can pass while the fusion is silently misaligned.

    Run it once per ``method``. Besides the offset bug it also catches a sign error: with
    ``TM_SQDIFF`` un-negated the planted location is the map's *minimum*, so the argmax lands
    on arbitrary noise and this raises immediately rather than 65 ROI-runs later.

    ``noise_frac`` is the canvas noise as a fraction of the template's **own** standard
    deviation (it replaces an absolute ``noise=8.0``, which silently assumed a 0-255 channel;
    see the comment on the canvas below).

    Returns a list of ``(Augmentation, dx, dy, score)``; raises on the first failure.
    """
    rng = np.random.default_rng(0) if rng is None else rng
    results = []
    cx = cy = canvas // 2

    for meta, tmpl in zip(metas, templates):
        # The canvas is drawn from the template's own distribution, not a fixed
        # N(120, 8). An absolute canvas only works for the methods that normalise it
        # away: under TM_CCORR, `sum(T*I)` is maximised wherever the image is brightest,
        # so planting a hematoxylin-OD template (values ~0.02) into a mean-120 canvas
        # puts the argmax on arbitrary background and the gate fails for the right
        # arithmetic reason. Matching the mean leaves the planted patch ahead by
        # `n * Var(T)`, which `noise_frac` keeps well outside the noise floor.
        loc = float(np.mean(tmpl))
        sd = float(np.std(tmpl)) * float(noise_frac)
        img = rng.normal(loc, sd if sd > 0 else 1e-6,
                         size=(canvas, canvas)).astype(np.float32)
        half = tmpl.shape[0] // 2
        img[cy - half: cy + half + 1, cx - half: cx + half + 1] = tmpl

        fused, _, valid = fused_response(img, templates, scale_normalize=scale_normalize,
                                         method=method)
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
