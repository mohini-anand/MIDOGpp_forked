"""
    Assert-style structural checks that catch bug classes a regression test cannot. See
    `INVARIANTS_HISTORY.md` for the defect each check was written against. Every function
    raises `InvariantError` (an `AssertionError`) and returns a dict of what it measured.
"""

from __future__ import annotations

from typing import Iterable, Optional

import numpy as np

from .evaluate import radius_px


class InvariantError(AssertionError):
    """A structural property this experiment relies on has stopped holding."""


def check_distinct_seeds(seed_records: Iterable, pool_size: Optional[int] = None, label: str = "") -> dict:
    """
        Invariant: distinct ``seed_index`` values must draw from distinct RNG streams.
        Deliberately not "distinct seed_index gives distinct seed_ann_id" -- that is not
        an invariant of a with-replacement draw and fails on correct runs.

        seed_records (Iterable): (seed_index, seed_ann_id, rng_spec) tuples for one image.
        pool_size (int or None): candidate pool size, for the all-identical-draw check.
        label (str): identifies the caller in the error message.

        Returns dict: stream/annotation counts, or raises InvariantError on a stream collision.
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


def check_nms_radius(used_radius: float, mpp: float, label: str = "", tol: float = 1e-6) -> dict:
    """Invariant: the NMS radius used must equal `evaluate.radius_px(mpp)`."""
    want = radius_px(mpp)
    if not np.isfinite(used_radius) or abs(float(used_radius) - want) > tol:
        raise InvariantError(
            f"{label}: NMS radius {used_radius} != evaluate.radius_px(mpp={mpp}) = "
            f"{want:.4f} -- suppression and scoring are using different neighbourhoods"
        )
    return {"check": "nms_radius", "label": label, "used_radius": float(used_radius),
            "expected_radius": float(want), "mpp": float(mpp), "passed": True}
