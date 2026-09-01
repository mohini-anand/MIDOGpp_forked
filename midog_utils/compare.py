"""One harness for comparing candidate generators and rankers on identical ground truth.

Every operating-point experiment in this repo before now built its own evaluation loop
inline, which is how it reached 2026-08-31 without ever running a seedless arm
(`Research Logs/2026-08-31-next-steps-plan.md`, Step 0). This module exists so that
adding a comparator is a two-line data declaration rather than a new script.

The contract
------------
* An `Arm` is **data**: a name, a callable returning a candidate frame with ``cx``, ``cy``
  and a ranking-key column, plus the metadata needed to interpret it.
* `evaluate_arms` takes **one** ``gt_eval`` and applies it to every arm, so arms differ
  only in what is under test. The seed annotation must already be removed from
  ``gt_eval`` by the caller -- and the *same* ``gt_eval`` must be passed for all arms in
  one comparison, including arms whose candidates do not depend on the seed at all.
* Each arm is sorted into its final ranking **before** bucketing.
  `evaluate.bucket_detections` matches greedily in row order, so bucketing a frame that
  is not yet in its final order silently scores a different ranking than the one
  reported.
* Sorting is always ``kind="mergesort"`` (stable) with ``na_position="last"``, so a tie
  block keeps the caller's incoming order and a run is reproducible.

Diagnostics that are mandatory rather than optional (Step 0 rule 8): ``floor_limited``,
``coverage_frac``, ``largest_tie_block`` and ``nan_rate``. The last two are computed on
the column that arm is actually ranked by -- which differs between arms -- not on one
column uniformly, since their whole purpose is to say how much of the reported ordering
is arbitrary.

Metric choice (Step 0 rule 4): ``recall_at_budget`` is primary. ``read_50`` and friends
read from the top of the list, so truncating a list barely moves them, which makes them a
poor control when the thing under comparison is the candidate *generator*. They are
reported because reading depth is the quantity a human pays, not because they arbitrate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np
import pandas as pd

from . import invariants as inv
from .evaluate import (
    HUMAN_CORRECT_LABEL,
    HUMAN_REJECTED_LABEL,
    bucket_detections,
    coverage_fraction,
)
from .dataset import LOOKALIKE, MITOTIC

BUDGETS = (500, 1000, 2000, 5000)
READ_FRACTIONS = (0.5, 0.8, 1.0)


@dataclass
class Arm:
    """One comparator, declared as data.

    ``candidates`` is a callable returning a frame with ``cx``, ``cy`` and ``rank_key``.
    A callable rather than a frame so a seedless arm can be built once per ROI and handed
    to every seed as a zero-cost closure, while a seeded arm can be built lazily.

    ``caps`` is every truncation limit that could have bound on this arm; the harness
    asserts the returned length equals none of them (`invariants.check_no_cap`).

    ``nms_radius`` is the suppression radius this arm's candidates were produced with, if
    any. When set, the harness asserts it equals ``evaluate.radius_px(mpp)``.

    ``z`` is ``nan`` for arms that do not depend on the extraction z level at all; those
    also carry ``z_dependent=False``. They are emitted **once** per (ROI, seed) rather
    than duplicated at every z, because duplicating them manufactures rows that look like
    independent measurements and are not.
    """

    name: str
    candidates: Callable[[], pd.DataFrame]
    rank_key: str = "od"
    seeded: bool = False
    z: float = float("nan")
    z_dependent: bool = False
    floor_limited: bool = False
    caps: Sequence = ()
    nms_radius: Optional[float] = None
    coverage_key: Optional[str] = None  # cache identity; None disables caching
    extra: dict = field(default_factory=dict)


def _rank(df: pd.DataFrame, key: str) -> pd.DataFrame:
    """Sort into the final ranking, best-first, stably, NaN last."""
    out = df.sort_values(key, ascending=False, na_position="last",
                         kind="mergesort").reset_index(drop=True)
    return out.assign(rank=np.arange(len(out)))


def _tie_and_nan(values: pd.Series) -> tuple:
    """Largest block of exactly equal ranking-key values, and the NaN rate.

    The two are reported separately and deliberately do not overlap: NaNs are excluded
    from the tie count, because ``nan_rate * n_detections`` already gives the size of the
    NaN block and folding it in would mask a genuine tie among real values behind a
    border-window artefact. Both mean the same thing for the reader -- that many list
    positions in a row are ordered arbitrarily -- but only one of them is a property of
    the ranking key itself.
    """
    if len(values) == 0:
        return 0, float("nan")
    nan_rate = float(values.isna().mean())
    counts = values.value_counts(dropna=True)
    return (int(counts.iloc[0]) if len(counts) else 0), nan_rate


def _read_depths(tp_cumulative: np.ndarray, n_mitotic: int) -> dict:
    """Candidates a reader must work through to reach each target sensitivity."""
    out = {}
    for frac in READ_FRACTIONS:
        label = f"read_{int(frac * 100)}"
        if n_mitotic == 0:
            out[label] = float("nan")
            continue
        want = int(np.ceil(frac * n_mitotic))
        i = int(np.searchsorted(tp_cumulative, want))
        out[label] = float(i + 1) if i < len(tp_cumulative) else float("nan")
    return out


def evaluate_arms(
    arms: Sequence[Arm],
    gt_eval: pd.DataFrame,
    radius: float,
    roi_shape=None,
    mpp: Optional[float] = None,
    budgets: Sequence[int] = BUDGETS,
    context: Optional[dict] = None,
    coverage_cache: Optional[dict] = None,
    checks: Optional[list] = None,
) -> pd.DataFrame:
    """Score every arm against one evaluation ground truth. Returns a tidy long frame.

    One row per (context..., arm, z, budget). Arm-level quantities (``read_50``,
    ``coverage_frac``, ``n_detections``, ...) repeat across that arm's budget rows;
    ``recall_at_budget`` is the only column that varies within them.

    ``coverage_cache`` is a caller-owned dict keyed by ``Arm.coverage_key``. Use it only
    for arms whose candidate set is genuinely identical between calls -- i.e. the seedless
    arms within one ROI. Keying it wrong reports one seed's coverage for all of them.

    ``checks`` collects the invariant records so a caller can log the evidence.
    """
    context = dict(context or {})
    checks = checks if checks is not None else []
    coverage_cache = coverage_cache if coverage_cache is not None else {}

    n_mit = int((gt_eval["category_id"] == MITOTIC).sum()) if len(gt_eval) else 0
    n_look = int((gt_eval["category_id"] == LOOKALIKE).sum()) if len(gt_eval) else 0

    rows = []
    label_stem = str(context.get("file_name", "?"))
    for arm in arms:
        det = arm.candidates()
        if arm.rank_key not in det.columns:
            raise KeyError(
                f"arm {arm.name!r} is ranked by {arm.rank_key!r}, which its candidate "
                f"frame does not have (columns: {list(det.columns)})"
            )
        checks.append(inv.check_no_cap(len(det), arm.caps,
                                       label=f"{label_stem}/{arm.name}"))
        if arm.nms_radius is not None and mpp is not None:
            checks.append(inv.check_nms_radius(arm.nms_radius, mpp,
                                               label=f"{label_stem}/{arm.name}"))

        ranked = _rank(det, arm.rank_key)
        largest_tie, nan_rate = _tie_and_nan(ranked[arm.rank_key])

        det_out, _ = bucket_detections(ranked, gt_eval, radius)
        b = det_out["bucket"].to_numpy()
        tp_cum = np.cumsum(b == HUMAN_CORRECT_LABEL)
        lk_cum = np.cumsum(b == HUMAN_REJECTED_LABEL)

        if roi_shape is None:
            coverage = float("nan")
        elif arm.coverage_key is not None and arm.coverage_key in coverage_cache:
            coverage = coverage_cache[arm.coverage_key]
        else:
            coverage = coverage_fraction(ranked[["cx", "cy"]].to_numpy(), roi_shape, radius)
            if arm.coverage_key is not None:
                coverage_cache[arm.coverage_key] = coverage

        base = dict(context)
        base.update({
            "arm": arm.name,
            "rank_key": arm.rank_key,
            "seeded": arm.seeded,
            "z": arm.z,
            "z_dependent": arm.z_dependent,
            "floor_limited": bool(arm.floor_limited),
            "n_gt_mitotic": n_mit,
            "n_gt_lookalike": n_look,
            "n_detections": int(len(ranked)),
            "coverage_frac": round(float(coverage), 4) if np.isfinite(coverage) else coverage,
            "n_lookalike_in_list": int(lk_cum[-1]) if len(lk_cum) else 0,
            "largest_tie_block": largest_tie,
            "nan_rate": round(nan_rate, 6) if np.isfinite(nan_rate) else nan_rate,
            "full_list_recall": float(tp_cum[-1] / n_mit) if (len(tp_cum) and n_mit) else (
                0.0 if n_mit else float("nan")),
            "match_radius_px": round(float(radius), 3),
        })
        base.update(_read_depths(tp_cum, n_mit))
        base.update(arm.extra)

        for k in budgets:
            delivered = min(int(k), len(tp_cum))
            tp_at = int(tp_cum[delivered - 1]) if delivered else 0
            lk_at = int(lk_cum[delivered - 1]) if delivered else 0
            row = dict(base)
            row.update({
                "budget": int(k),
                "budget_delivered": delivered,
                "tp_at_budget": tp_at,
                "lookalike_at_budget": lk_at,
                "recall_at_budget": float(tp_at / n_mit) if n_mit else float("nan"),
            })
            rows.append(row)

    return pd.DataFrame(rows)


def assert_floor_not_limiting(df: pd.DataFrame) -> None:
    """`floor_limited` must be uniformly False when the floor is in per-image z units.

    With an extraction floor at ``med + 0.5*mad`` and every reported z >= 1.0, no cell can
    be limited by the floor. This is a verification, not a formality: a True here means
    the floor and the threshold are no longer in the same units, which is the exact bug
    class this experiment exists to characterise.
    """
    bad = df[df["floor_limited"].astype(bool)]
    if len(bad):
        raise inv.InvariantError(
            f"{len(bad)} rows are floor_limited, which cannot happen when the extraction "
            f"floor sits below every reported z: {bad[['arm', 'z']].drop_duplicates()}"
        )


def summarise(df: pd.DataFrame, by=("file_name", "arm", "z", "budget"),
              value="recall_at_budget") -> pd.DataFrame:
    """Median / IQR / n over seeds. Never report a single-seed number (Step 0 rule 5)."""
    g = df.groupby(list(by), dropna=False)[value]
    out = g.agg(median="median", q1=lambda s: s.quantile(0.25),
                q3=lambda s: s.quantile(0.75), n="size",
                n_unreachable=lambda s: int(s.isna().sum()),
                min="min", max="max").reset_index()
    out["iqr"] = out["q3"] - out["q1"]
    return out
