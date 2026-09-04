"""Seed selection under pathologist agreement and bbox-tightening filters.

Two independent filters narrow the candidate pool below what `experiment.pick_seed`
draws from uniformly at random, applied in this order:

1. **Pathologist agreement** (`agreement_pool`): prefer a mitotic annotation every rater
   who saw it called mitotic (`n_mitotic_votes == n_votes`). If none exist in the image,
   fall back to the majority-mitotic contested set (>= 2/3 of votes mitotic) and flag the
   image -- a seed the experts themselves argued about is plausibly a morphologically
   atypical example, and using it as the *only* template risks building a poor template
   from the start.
2. **Bbox tightening** (`tighten_box_otsu`, `foreground_filter`): Otsu-threshold the
   padded patch and require the click to land inside the connected component that
   produces. An annotation whose click sits outside its own object's foreground would
   need the largest-connected-component fallback the reference implementation
   (`bbox tuning code reference/`) uses to recover a box at all -- refused here, not
   applied; such annotations are dropped from the candidate pool instead of being built
   into a template from the wrong object.

`tightened_base_size` turns the same Otsu+CC bbox into a `FSConfig(base_size=...)`
value for the actual template-matching pass -- the tightened box's *native size*, kept
centred on the click rather than recentred to the component's own centroid, since only
the scale is being corrected here.

See `Research Logs/design_choices.md`, sections 1-2.
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
    """Split an image's mitotic annotations into the agreement tier to draw from.

    Returns ``(pool, flagged)``. ``pool`` is every annotation every rater who saw it
    called mitotic (``n_mitotic_votes == n_votes``), when that set is non-empty;
    otherwise it is the majority-mitotic contested set (>= 2/3 of votes mitotic -- which,
    given only 2- and 3-rater annotations exist, means the 2-of-3 case) and ``flagged``
    is True, meaning the caller should record that this image had no fully-agreed
    annotation to seed from.
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
    """Annotations far enough from the ROI edge to read a full padded patch.

    Same predicate as `experiment.pick_seed` -- duplicated rather than imported from
    there so this module has no dependency on `experiment`, and tested on the *rounded*
    centre for the same reason documented there: the reader's own predicate is on the
    rounded pixel, and testing the unrounded float could pass an annotation the reader
    then rejects.
    """
    h, w = roi_shape[:2]
    ix = np.rint(df["cx"].to_numpy()).astype(int)
    iy = np.rint(df["cy"].to_numpy()).astype(int)
    ok = (ix >= border) & (ix <= w - 1 - border) & (iy >= border) & (iy <= h - 1 - border)
    return df[ok]


def _nearest_label_within(labels: np.ndarray, cy: int, cx: int, tolerance: int) -> int:
    """The foreground label closest to ``(cy, cx)`` within an L-inf ``tolerance``, or 0.

    Ties (equal squared-distance pixels with different labels) resolve to whichever
    `np.nonzero` visits first -- row-major order -- which is an arbitrary but stable
    choice; ties are rare enough at these tolerances (1-3 px) not to warrant more.
    """
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


def tighten_box_otsu(
    patch: np.ndarray,
    min_area: int = 50,
    max_area_frac: float = 0.85,
    min_solidity: float = 0.5,
    method: str = "binary",
    center_tolerance: int = 0,
    headroom_frac: float = None,
):
    """Threshold ``patch`` and return the connected component under (or near) its centre.

    ``patch`` is a single-channel array in "more object -> higher value" convention (a
    `channels.to_gray_inverted` or `channels.to_hematoxylin` crop -- not `channels.to_rgb`,
    which this raises on). Returns ``(y0, y1, x0, x1)`` -- half-open, patch-local -- or
    ``None`` when the centre pixel itself is not foreground, or when the component it
    lands in fails the size/shape sanity check below. Not-foreground is the case the
    reference implementation resolves by falling back to the largest component in the
    crop; that fallback is refused here, so a ``None`` means the caller should treat the
    annotation as unusable as a seed.

    ``method`` selects the foreground threshold:

    * ``"binary"`` (default): the original two-class `cv2.THRESH_OTSU` split -- the
      patch's whole non-background range becomes foreground.
    * ``"multiotsu"``: `skimage.filters.threshold_multiotsu` with 3 classes, keeping only
      the brightest class as foreground. Binary Otsu can under-separate a component that
      bridges into a neighbouring structure through a lighter, swept-in fringe -- a
      visibly darker core distinct from that fringe reads as a threshold-placement
      problem the middle multi-Otsu class absorbs instead of merging into the core. See
      `Research Logs/design_choices.md`, section 7. Raises no further than returning
      ``None`` when the patch has fewer unique intensities than classes (a near-flat
      window multi-Otsu can't threshold) -- refused rather than silently falling back to
      binary Otsu for that one patch, consistent with the rest of this function's
      "exclude, don't patch over" handling.
    * ``"headroom"``: binary Otsu's own split, raised by a fraction of the headroom
      between that split and the crop's max (always 255, since the crop is normalized
      to [0, 255] before thresholding either way) -- ``threshold = otsu_thresh +
      headroom_frac * (255 - otsu_thresh)``, foreground is ``u8 > threshold`` (strict,
      matching ``cv2.THRESH_BINARY``'s own semantics). A dialable middle ground between
      ``"binary"`` and ``"multiotsu"``: on the window-saturating merges ``"multiotsu"``
      was built for (a component spanning several adjacent nuclei rather than one, e.g.
      `245.tiff` ann 6245's entire 51x51 window at area_frac=0.58), raising the
      threshold un-merges most of them without multiotsu's separate 3-class fit -- but
      the same intervention also shrinks or drops components that were never merged in
      the first place, and that collateral cost lands in the same dense domains as the
      benefit. ``frac=0`` reproduces ``"binary"`` exactly (verified byte-for-byte).
      ``frac=0.15`` resolves roughly half the window-saturating population at about a
      third of multiotsu's collateral exclusion rate in the dense domains that matter;
      ``frac=0.30`` lands within ~3 points of multiotsu on every metric measured. See
      `Research Logs/2026-09-03-bbox-threshold-sweep.md` for the full sweep (frac in
      {0.05, 0.10, 0.15, 0.20, 0.30} vs. multiotsu, 933 candidates across 18 ROIs) and
      its visual audit -- including cases where the flagged-large component was not
      actually a merge, just a legitimately large chromatin mass the size-based flag
      can't tell apart from one. ``headroom_frac`` must be a float in ``[0, 1]``;
      raises ``ValueError`` otherwise (mirroring the bad-``method`` check below), since
      an unset or out-of-range fraction almost certainly means the caller forgot to
      pass it, not that some sensible default was intended. ``0.0`` is accepted (not
      just documented as equivalent to it): Otsu's threshold is always integer-valued
      on a ``uint8`` image, so ``threshold = otsu_thresh`` under this branch's strict
      ``>`` gives byte-identical output to ``method="binary"`` -- deliberately, so a
      sweep over ``headroom_frac`` (e.g. `Research Logs/2026-09-03-bbox-threshold-sweep.md`'s
      own {0.05, ..., 0.30} grid) can include 0 as its own left endpoint rather than
      needing a workaround value to stand in for it.

    ``center_tolerance`` (default 0, exact-pixel-only, unchanged) widens the centre
    check to the nearest foreground pixel within an L-inf square of this half-width
    around the click, rather than requiring the exact rounded click pixel itself to be
    foreground. A pathologist's click is a point annotation at the object's nominal
    centre, not necessarily its brightest pixel; under `method="multiotsu"`'s tighter
    top-class boundary, calibration found roughly half of the clicks it dropped had a
    foreground pixel 1 px away and the large majority had one within 3 px -- i.e. the
    exact-pixel test itself, not the absence of the object, was responsible for a
    meaningful share of those losses. See `Research Logs/design_choices.md`, section 7a.
    A real cost of widening this, also measured there: the accepted component's mask
    never contains the click by construction (if it did, tolerance wouldn't have been
    needed) -- calibration found its *bounding box* doesn't even reach the click for
    about 1 in 6 candidates tolerance alone recovers, i.e. that component more plausibly
    belongs to a neighbouring structure than the annotated one. Guarded against below:
    the accepted component is rejected outright when the click falls outside its own
    bounding box (see the bbox check just after ``center_label`` is resolved). This closes
    the bbox-containment failure completely, but not the softer version of the same risk
    -- section 7b found roughly 1 in 5 guard-recovered candidates still have the click
    farther from the component's centroid than the component's own approximate radius
    (an elongated neighbour whose bbox reaches the click while its mass sits elsewhere
    would still pass). `Research Logs/design_choices.md`, section 7b has the numbers.

    The sanity check mirrors `baselines.nucleus_blobs`'s ``min_area``/``max_area`` gate,
    since an oversized or malformed component can pass the centre-pixel check above and
    still not be a single nucleus:

    * ``min_area`` rejects a degenerate sliver -- Otsu noise the click's rounded centre
      happened to land on, not real chromatin.
    * ``max_area_frac`` (of the window's own pixel count) rejects a component that has
      grown implausibly large for the window.
    * ``min_solidity`` (filled area / convex-hull area) rejects a markedly non-convex
      shape -- the signature of two touching nuclei bridged into one component by Otsu.
      Eccentricity is deliberately *not* used for this: a real mitotic figure can be
      legitimately elongated or irregular (anaphase/telophase chromatin, for instance),
      which reads as high eccentricity on a single, genuine object; solidity does not
      flag a convex ellipse like that; it only flags the concave, dumbbell-shaped
      bridge a merge produces.

    Defaults are calibrated on ~450 unanimous mitotic candidates across 10 downloaded
    ROIs (dense and sparse domains both) -- see `Research Logs/design_choices.md`.
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
    # A tolerance-recovered label is, by construction, one the click's own pixel is NOT
    # part of -- it can belong to a neighbouring structure rather than the annotated one
    # (`Research Logs/design_choices.md`, section 7a). Requiring the click to at least
    # fall inside the accepted component's own bounding box is a cheap, partial identity
    # check: it does nothing when the exact pixel was already foreground (that click is
    # trivially inside its own bbox), and only excludes tolerance-recovered candidates
    # whose nearest component doesn't actually reach the click.
    min_row, min_col, max_row, max_col = region.bbox
    if not (min_row <= cy < max_row and min_col <= cx < max_col):
        return None
    if region.area < min_area:
        return None
    if region.area > max_area_frac * patch.size:
        return None
    if region.solidity < min_solidity:
        return None

    y0, x0, y1, x1 = region.bbox  # skimage bbox is already half-open: (min_row, min_col, max_row, max_col)
    return int(y0), int(y1), int(x0), int(x1)


def foreground_filter(
    df: pd.DataFrame,
    structural_channel: np.ndarray,
    otsu_window: int = tm.BASE_SIZE,
    method: str = "binary",
    center_tolerance: int = 0,
    headroom_frac: float = None,
) -> pd.DataFrame:
    """Keep only annotations whose click lands inside (or near) its own Otsu component.

    ``structural_channel`` is the single-channel image `tighten_box_otsu` runs on --
    independent of whatever channel `FSConfig.channel` later searches in, since this is a
    structural test of the click's location, not a matching score.

    ``method``, ``center_tolerance``, and ``headroom_frac`` are passed straight through
    to `tighten_box_otsu` -- ``method`` is ``"binary"`` (default), ``"multiotsu"``, or
    ``"headroom"``; ``center_tolerance`` (default 0) is the exact-pixel-vs-nearby-pixel
    centre check; ``headroom_frac`` (default ``None``) is `"headroom"`'s required
    ``[0, 1]`` fraction, unused by the other two methods.

    ``otsu_window`` defaults to the 51 px annotation box itself, *not*
    `template_match.PATCH_SIZE` (73 px) -- despite the design doc's "Otsu-threshold the
    padded crop" wording, which reads that way. Measured on 246.tiff (lymphosarcoma,
    dense: 97 unanimous candidates), the 73 px window saturates 60/97 components at the
    patch edge (mean tightened size 61 px, i.e. barely smaller than the padded read
    itself -- not a tightening) because dense tissue gives Otsu enough neighbouring
    nuclei to bridge into within 73 px. The 51 px window keeps sizes inside the box
    (mean 43 px) at the cost of 8/97 clicks whose component doesn't cross threshold at
    all in the smaller window (excluded here, not fallen back on). See the notebook's
    side-by-side comparison.
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


def pick_seed(
    gt_mitotic: pd.DataFrame,
    structural_channel: np.ndarray,
    rng,
    border: int,
    roi_shape,
    otsu_window: int = tm.BASE_SIZE,
    method: str = "binary",
    center_tolerance: int = 0,
    headroom_frac: float = None,
    tighten_bbox: bool = True,
):
    """Pick a mitotic seed under the pathologist-agreement and bbox-tightening filters.

    Order: agreement tiering first (on the image's full mitotic set), then the border
    filter, then the foreground filter. Any of the three can empty the pool; unlike
    `experiment.pick_seed`'s single border check, this raises naming the stage that
    emptied it rather than silently returning from whatever tier happened to survive --
    see `Research Logs/design_choices.md`'s caveat about the pool shrinking to zero on
    sparse domains.

    ``border`` is independent of ``otsu_window`` -- it should still be the caller's full
    rotation-safe half-patch (`FSConfig.patch_size // 2`), since `find_and_suppress` reads
    that much regardless of how large the tightened template turns out to be.

    ``method``, ``center_tolerance``, and ``headroom_frac`` are passed straight through
    to `tighten_box_otsu`/`foreground_filter` -- ``method`` is ``"binary"`` (default),
    ``"multiotsu"``, or ``"headroom"``; ``center_tolerance`` (default 0) is the
    exact-pixel-vs-nearby-pixel centre check; ``headroom_frac`` (default ``None``) is
    `"headroom"`'s required ``[0, 1]`` fraction, unused by the other two methods.

    ``tighten_bbox`` (default True, unchanged) toggles the foreground filter entirely --
    when False, seed selection is pathologist agreement + the border filter only, no
    Otsu/CC step at all (`method`/``center_tolerance``/``headroom_frac`` are then
    unused). ``n_fg`` is left
    equal to ``n_border`` in that case rather than a separate, always-identical number, so
    `SeedInfo` still means "pool size after this stage" consistently whether or not the
    stage actually ran.

    Returns ``(seed, info)`` where ``info`` is a `SeedInfo` recording the pool size after
    each stage, for logging.
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


def tightened_base_size(structural_channel: np.ndarray, cx: float, cy: float,
                         otsu_window: int = tm.BASE_SIZE, minimum: int = 5,
                         method: str = "binary", center_tolerance: int = 0,
                         headroom_frac: float = None):
    """The native template size for a seed after bbox tightening.

    Otsu-thresholds the ``otsu_window``-sized box around the click, takes the connected
    component under it, and returns the odd size spanning its longer side -- for
    `FSConfig(base_size=...)`. Kept centred on the click rather than recentred to the
    component's own centroid, since only the *scale* is being corrected
    (`Research Logs/design_choices.md`: "scale is dropped -- matching uses each
    tightened box's native size").

    ``otsu_window`` defaults to `template_match.BASE_SIZE` (51 px, the annotation box),
    not `template_match.PATCH_SIZE` (73 px) -- see `foreground_filter`'s docstring for
    why the wider window was rejected.

    ``method``, ``center_tolerance``, and ``headroom_frac`` must match whatever
    `pick_seed` used to accept this seed -- passing different ones here can retighten a
    box under settings that never actually validated this click's foreground membership.
    ``method`` is ``"binary"`` (default), ``"multiotsu"``, or ``"headroom"``;
    ``headroom_frac`` (default ``None``) is `"headroom"`'s required ``[0, 1]`` fraction,
    unused by the other two methods.

    Returns ``None`` when the patch can't be read or the click isn't foreground;
    `pick_seed`'s foreground filter means a seed it returned will not hit the second case
    when the same ``method``/``center_tolerance``/``headroom_frac`` are passed to both.
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
