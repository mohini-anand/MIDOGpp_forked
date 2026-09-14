"""Assert-style checks that catch the bug classes a regression test cannot.

`Research Logs/2026-08-31-next-steps-plan.md`, Step 0 rule 11: a test that re-runs seed 0
and asserts previously committed numbers cannot catch the 2026-08-25 bug class *by
construction* -- those bugs were present when the baseline was committed, so such a test
would lock them in. The four checks below assert structural properties instead, each one
tied to a specific defect that actually occurred in this repo:

1. `check_no_cap` -- a list length exactly equal to a configured maximum is almost never a
   coincidence. `baselines.nucleus_blobs` used to default to ``max_detections=20000`` and
   silently bound on the two densest ROIs, truncating to the *darkest* competitors and so
   biasing every rank statistic derived from it.
2. `check_tissue_mask_covers_gt` -- global Otsu as a tissue mask put 31 of 301.tiff's 218
   mitotic figures outside "tissue", where `baselines.nucleus_blobs` cannot see them.
   The fixed brightness cut is supposed to exclude zero ground-truth annotations; this
   asserts it on every ROI it is used on rather than trusting the one-off validation.
   Running it found that the *full-coverage* form of that claim is false -- one category-2
   look-alike on 506.tiff is excluded -- so the check is fatal on mitotic annotations and
   records the rest. See its docstring; the exception is real and is reported, not hidden.
3. `check_distinct_seeds` -- an RNG rebuilt inside a loop (``np.random.default_rng(0)``
   per iteration) returns the same draw every time, so a "5-seed sweep" silently becomes
   one seed measured five times. Note this is asserted as "distinct ``seed_index`` values
   use distinct RNG *streams*", not the tempting "distinct ``seed_index`` gives distinct
   ``seed_ann_id``" -- see that function's docstring for why the latter is not an
   invariant of a with-replacement draw and fails on correct runs.
4. `check_nms_radius` -- a constant 25.0 px NMS radius let two detections 26 px apart both
   survive while both sat inside one ground-truth object's match radius (29.6-33.1 px
   across these scanners), inflating the FROC's false-positive axis with duplicates. The
   suppression radius must be the image's own evaluation match radius.

Every function raises `InvariantError` (an `AssertionError`) and returns a small dict of
what it measured, so a caller can log the evidence as well as the pass.
"""

from __future__ import annotations

from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd

from .evaluate import radius_px


class InvariantError(AssertionError):
    """A structural property this experiment relies on has stopped holding."""


def check_no_cap(n_items: int, caps: Sequence, label: str = "") -> dict:
    """Invariant 1: no arm's list length may equal a configured cap.

    ``caps`` is every truncation limit that could have bound on this arm
    (``FSConfig.max_peaks``, ``FSConfig.max_detections``, ``nucleus_blobs``'s
    ``max_detections``, ...). ``None`` entries -- meaning "no limit configured" -- are
    skipped. An exact equality is the signal; a length one below a cap is fine.
    """
    hits = [int(c) for c in caps if c is not None and int(c) == int(n_items)]
    if hits:
        raise InvariantError(
            f"{label}: candidate list length {n_items} exactly equals configured cap(s) "
            f"{hits} -- the list was almost certainly truncated, so every rank statistic "
            "derived from it is conditioned on the truncation"
        )
    return {"check": "no_cap", "label": label, "n_items": int(n_items),
            "caps": tuple(int(c) for c in caps if c is not None), "passed": True}


def check_tissue_mask_covers_gt(mask: np.ndarray, gt: pd.DataFrame, label: str = "",
                                strict_categories=(1,)) -> dict:
    """Invariant 2: `baselines.tissue_mask` must exclude 0 ground-truth annotations.

    Tested on the *rounded* click pixel, which is the pixel a detector reading the mask
    would actually consult. Points outside the array bounds count as excluded.

    ``strict_categories`` is which categories make an exclusion fatal, and defaults to
    category 1 (mitotic) alone. **Measured, not assumed:** the full-coverage form of this
    invariant is false on this dataset. 506.tiff's annotation 24249 -- a *category-2*
    look-alike at x = 13, i.e. 13 px from the ROI's left edge, where the click pixel
    reads gray 232 against the 220 cut and only 30% of its 51 px window is inside the
    mask -- is excluded. It is the only exclusion across all seven ROIs (1 of 726
    annotations pooled; 0 of 391 mitotic). `baselines.tissue_mask`'s own docstring claims
    coverage of *mitotic* annotations specifically, and that claim holds exactly.

    The distinction is not cosmetic. Every recall number in this experiment is computed
    over category-1 ground truth, so a mitotic figure outside the mask silently caps what
    a mask-restricted arm (`nucleus_blobs`, `random_in_tissue`) can score -- that must
    abort. A look-alike outside the mask affects only the ``n_lookalike_in_list``
    diagnostic and caps it by one. Both counts are returned either way, so the softer
    case is recorded and reported rather than swallowed; pass ``strict_categories=None``
    to require full coverage of every category.
    """
    empty = {"check": "tissue_mask_covers_gt", "label": label, "n_gt": 0,
             "n_excluded": 0, "n_excluded_strict": 0, "excluded_ann_ids": "",
             "passed": True}
    if len(gt) == 0:
        return empty
    h, w = mask.shape[:2]
    ix = np.rint(gt["cx"].to_numpy()).astype(int)
    iy = np.rint(gt["cy"].to_numpy()).astype(int)
    inside = (ix >= 0) & (ix < w) & (iy >= 0) & (iy < h)
    covered = np.zeros(len(gt), dtype=bool)
    covered[inside] = mask[iy[inside], ix[inside]]

    out = gt[~covered]
    if strict_categories is None:
        fatal = out
    else:
        fatal = out[out["category_id"].isin(list(strict_categories))]
    if len(fatal):
        raise InvariantError(
            f"{label}: tissue_mask excludes {len(fatal)} of {len(gt)} ground-truth "
            f"annotations in categories {strict_categories} (ann_ids "
            f"{sorted(fatal['ann_id'].tolist())}) -- every seedless arm restricted to "
            "this mask is being scored against ground truth it was never allowed to reach"
        )
    return {"check": "tissue_mask_covers_gt", "label": label, "n_gt": int(len(gt)),
            "n_excluded": int(len(out)), "n_excluded_strict": 0,
            "excluded_ann_ids": ";".join(
                f"{int(r.ann_id)}(cat{int(r.category_id)})" for r in out.itertuples()),
            "passed": True}


def check_distinct_seeds(seed_records: Iterable, pool_size: Optional[int] = None,
                         label: str = "") -> dict:
    """Invariant 3: distinct ``seed_index`` values must draw from distinct RNG streams.

    ``seed_records`` is an iterable of ``(seed_index, seed_ann_id, rng_spec)`` for one
    image, where ``rng_spec`` is the seed sequence handed to `np.random.default_rng`.

    **This is deliberately not the literal check "distinct seed_index gives distinct
    seed_ann_id".** That statement is not an invariant of this experiment and asserting it
    would fail on correct runs: seed selection draws *uniformly with replacement* from a
    pool of accepted annotations, so a collision is expected, not a defect. The pools here
    run 3-162 candidates and 5 draws collide with probability 1 - prod(1 - i/p), which is
    6% at p=162 (301.tiff) and certain at p<5 (350.tiff has 3 mitotic annotations after
    the seed is excluded). `results/od_seed_sweep.csv`, produced by already-committed
    code with the exact RNG construction this experiment uses, shows 4 of 7 ROIs drawing
    only 4 distinct annotations across 5 seed indices -- correct behaviour that the literal
    check would report as a bug.

    What is asserted instead is the property the literal check was a proxy for, in two
    parts, both of which an RNG rebuilt inside the loop violates with certainty:

    1. **The streams themselves differ.** Each ``rng_spec`` is replayed and its first 8
       draws compared; two seed indices sharing a stream is the ``default_rng(0)``-inside-
       the-loop bug exactly, and is deterministic rather than probabilistic.
    2. **Not every draw is the same annotation**, whenever the pool offers a choice.

    ``n_distinct_ann`` is returned so legitimate collisions stay visible in the audit
    trail instead of being silently tolerated.
    """
    recs = [(int(i), int(a), tuple(np.atleast_1d(spec).tolist()))
            for i, a, spec in seed_records]

    streams = {}
    for idx, _, spec in recs:
        draws = tuple(np.random.default_rng(list(spec)).integers(0, 2 ** 62, 8).tolist())
        if draws in streams and streams[draws] != idx:
            raise InvariantError(
                f"{label}: seed_index {streams[draws]} and {idx} produce an identical RNG "
                f"stream from specs {spec!r} -- the generator is being rebuilt per "
                "iteration, so this is one seed measured several times, not a sweep"
            )
        streams[draws] = idx

    n_distinct = len({a for _, a, _ in recs})
    if len(recs) > 1 and n_distinct == 1 and (pool_size is None or pool_size > 1):
        raise InvariantError(
            f"{label}: all {len(recs)} seed indices drew the same seed_ann_id from a pool "
            f"of {pool_size} candidates -- the draw is not varying with seed_index"
        )
    return {"check": "distinct_seeds", "label": label, "n_seeds": len(recs),
            "n_distinct_ann": n_distinct, "pool_size": pool_size,
            "n_distinct_streams": len(streams), "passed": True}


def check_nms_radius(used_radius: float, mpp: float, label: str = "",
                     tol: float = 1e-6) -> dict:
    """Invariant 4: the NMS radius used must equal `evaluate.radius_px(mpp)`."""
    want = radius_px(mpp)
    if not np.isfinite(used_radius) or abs(float(used_radius) - want) > tol:
        raise InvariantError(
            f"{label}: NMS radius {used_radius} != evaluate.radius_px(mpp={mpp}) = "
            f"{want:.4f} -- suppression and scoring are using different neighbourhoods"
        )
    return {"check": "nms_radius", "label": label, "used_radius": float(used_radius),
            "expected_radius": float(want), "mpp": float(mpp), "passed": True}


def check_min_separation(centers: np.ndarray, radius: float, label: str = "") -> dict:
    """Post-NMS geometry: no two kept points may lie within ``radius`` of each other.

    Not one of the four required invariants -- this is the direct check that an arm
    declared "NMS'd at the match radius" actually was, used on the lattice arms where the
    suppression is applied by the caller rather than inside `find_and_suppress`.
    """
    centers = np.asarray(centers, dtype=np.float64).reshape(-1, 2)
    if len(centers) < 2:
        return {"check": "min_separation", "label": label, "n": int(len(centers)),
                "min_distance": float("inf"), "radius": float(radius), "passed": True}
    from sklearn.neighbors import KDTree
    dist, _ = KDTree(centers).query(centers, k=2)
    d_min = float(dist[:, 1].min())
    if d_min <= radius:
        raise InvariantError(
            f"{label}: two kept points are {d_min:.3f} px apart, within the "
            f"{radius:.3f} px NMS radius -- suppression was not applied at the match radius"
        )
    return {"check": "min_separation", "label": label, "n": int(len(centers)),
            "min_distance": d_min, "radius": float(radius), "passed": True}
