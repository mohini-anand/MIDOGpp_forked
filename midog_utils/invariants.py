"""
    Assert-style structural checks that catch bug classes a regression test cannot. See
    `INVARIANTS_HISTORY.md` for the defect each check was written against. Every function
    raises `InvariantError` (an `AssertionError`) and returns a dict of what it measured.
"""

from __future__ import annotations

from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd

from .evaluate import radius_px


class InvariantError(AssertionError):
    """A structural property this experiment relies on has stopped holding."""


def check_no_cap(n_items: int, caps: Sequence, label: str = "") -> dict:
    """
        Invariant: no arm's list length may equal a configured cap.

        n_items (int): the list length to check.
        caps (Sequence): every truncation limit that could have bound on this arm; None
            entries (no limit configured) are skipped.
        label (str): identifies the caller in the error message.

        Returns dict: what was checked, or raises InvariantError on an exact-cap hit.
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


def check_tissue_mask_covers_gt(mask: np.ndarray, gt: pd.DataFrame, label: str = "", strict_categories=(1,)) -> dict:
    """
        Invariant: `baselines.tissue_mask` must exclude 0 ground-truth annotations in
        ``strict_categories`` (default: category 1, mitotic).

        mask (np.ndarray): the tissue mask.
        gt (pd.DataFrame): ground-truth annotations to check, tested on the rounded click pixel.
        label (str): identifies the caller in the error message.
        strict_categories (tuple or None): which categories make an exclusion fatal; None
            requires full coverage of every category.

        Returns dict: exclusion counts, or raises InvariantError on a fatal exclusion.
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


def check_min_separation(centers: np.ndarray, radius: float, label: str = "") -> dict:
    """Post-NMS geometry: no two kept points may lie within ``radius`` of each other."""
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
