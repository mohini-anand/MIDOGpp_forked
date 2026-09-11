"""F7 stage 2: apply the size prune and score it. Imported by `f7_size_prune.ipynb`.

Every function here is a pure deletion on the cached, score-descending pool followed by a
re-match. Nothing is re-ranked: the pool is already TM_CCOEFF-descending, so the survivors
keep their order and "prune then re-rank by TM_CCOEFF" is the identity on them.

Re-matching after the prune is *required*, not optional: `evaluate.greedy_match` lets each
annotation be claimed by the best-ranked detection within the match radius, so deleting a
candidate that was a claimant can hand its annotation to a deeper one -- or lose it. Reusing
the committed ledger's ranks after a deletion would silently assume that never happens.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from midog_utils import dataset as ds
from midog_utils import evaluate as ev

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # repo root -- every shared artefact path is absolute,
RES = os.path.join(HERE, "results")   # so the notebook can run from its own directory

APRIORI_LO, APRIORI_HI = 4.0, 18.0      # pre-registered in PREREGISTRATION.md section 6
QUANTS = (1.00, 0.99, 0.98, 0.95)       # TP-retention levels the bound is set at
LEVELS = (0.90, 0.95, 1.00)             # recall levels the depth is read at


def load_all():
    """Every artefact F7 needs, plus `gt_eval` per cell. No images are read."""
    gate = pd.read_csv(os.path.join(RES, "f7_repro_gate.csv"))
    sizes = np.load(os.path.join(RES, "f7_candidate_sizes.npz"))
    pool = np.load(os.path.join(ROOT, "results", "tm_recall_workload_pool_z.npz"))
    ledger = pd.read_csv(os.path.join(ROOT, "results", "tm_recall_workload_tp_ledger.csv"))
    _, ann = ds.load_annotations(os.path.join(ROOT, "databases", "MIDOG++.json"))
    return gate, sizes, pool, ledger, ann


def gt_eval_for(ann: pd.DataFrame, fn: str, seed_ann_id: int) -> pd.DataFrame:
    """Rebuilt exactly as `recall_workload_ledger.py` builds it."""
    gt = ds.image_annotations(ann, fn)
    return gt[gt["ann_id"] != seed_ann_id].reset_index(drop=True)


def keep_mask(size_um: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """`lo <= size <= hi`, with **undefined sizes kept** (PREREGISTRATION.md section 5).

    Pruning a candidate whose size could not be measured would be a border/nucleus rule
    smuggled in as a second knob.
    """
    return np.isnan(size_um) | ((size_um >= lo) & (size_um <= hi))


def depths_from_tp_ranks(tp_depths: np.ndarray, n_mit: int, levels=LEVELS) -> dict:
    """Depth (1-based list length) at which each recall level is reached; NaN if unreachable.

    `tp_depths` must be the sorted 1-based positions of the true positives in the ranked
    list. The denominator is `n_mit`, fixed across arms, so a lost mitosis shows as a rising
    depth or a NaN rather than being hidden by a shrinking denominator.
    """
    out = {}
    d = np.sort(np.asarray(tp_depths))
    for q in levels:
        need = int(np.ceil(q * n_mit))
        out[f"depth_{int(q*100)}"] = float(d[need - 1]) if len(d) >= need >= 1 else np.nan
    return out


def score_arm(cx, cy, size_um, gt_eval, match_radius, n_mit, lo, hi) -> dict:
    """Prune, re-match, and read the depth curve. `lo=-inf, hi=inf` gives the baseline."""
    keep = keep_mask(size_um, lo, hi)
    det, _ = ev.bucket_detections(pd.DataFrame({"cx": cx[keep], "cy": cy[keep]}),
                                  gt_eval, match_radius)
    tp = det[det["bucket"] == ev.HUMAN_CORRECT_LABEL]
    # `det` preserves the incoming row order and carries a RangeIndex, so a TP's index is its
    # 0-based rank among the survivors; +1 is the list length a reader would be handed.
    tp_depths = tp.index.to_numpy() + 1
    row = dict(n_kept=int(keep.sum()), n_tp=len(tp_depths),
               recall=len(tp_depths) / n_mit if n_mit else np.nan)
    row.update(depths_from_tp_ranks(tp_depths, n_mit))
    return row


def tp_sizes_for_cell(sizes, ledger, fn, si) -> np.ndarray:
    """Sizes of the candidates that claimed a mitosis in the UNPRUNED pool, NaNs dropped.

    Uses the committed ledger's 1-based ranks, which the reproduction gate has verified.
    """
    r = ledger[(ledger.file_name == fn) & (ledger.seed_index == si)]["rank"].to_numpy() - 1
    v = sizes[f"{fn}|{si}|eqd_um"][r]
    return v[~np.isnan(v)]


def bounds_from_tp_sizes(tp_sizes: np.ndarray, retention: float) -> tuple[float, float]:
    """Symmetric quantile bound retaining `retention` of the measurable TPs.

    At retention 1.00 this is [min, max]; note that is a min/max order statistic and one
    atypical mitosis pins it -- which is exactly why the whole curve is reported.
    """
    if len(tp_sizes) == 0:
        return -np.inf, np.inf
    a = (1.0 - retention) / 2.0
    return float(np.quantile(tp_sizes, a)), float(np.quantile(tp_sizes, 1.0 - a))


def run_regimes(gate, sizes, pool, ledger, ann) -> pd.DataFrame:
    """One row per (cell, arm). Arms: baseline, a-priori, in-sample, LODO -- the last two
    at every retention level in QUANTS.

    LODO holds out a whole **domain** (both its ROIs), per `fp_filter_domain.py`: the 14 ROIs
    come in domain pairs sharing tumour type *and* scanner, so leave-one-ROI-out would measure
    "generalises to another ROI of a tumour I have already seen", not deployment.
    """
    cellinfo, tps = {}, {}
    for _, g in gate.iterrows():
        fn, si = g.file_name, int(g.seed_index)
        cellinfo[(fn, si)] = g
        tps[(fn, si)] = tp_sizes_for_cell(sizes, ledger, fn, si)

    domains = gate.groupby("tumor_type")["file_name"].unique().to_dict()

    rows = []
    for (fn, si), g in cellinfo.items():
        cx, cy = pool[f"{fn}|{si}|cx"], pool[f"{fn}|{si}|cy"]
        size_um = sizes[f"{fn}|{si}|eqd_um"]
        gt_eval = gt_eval_for(ann, fn, int(g.seed_ann_id))
        n_mit = int(g.n_gt_mitotic)
        mr = float(g.match_radius_px)
        base = dict(file_name=fn, tumor_type=g.tumor_type, seed_index=si, n_mit=n_mit,
                    n_pool=int(g.n_pool), frac_on_nucleus=float(g.frac_on_nucleus))

        arms = [("baseline", -np.inf, np.inf, np.nan),
                ("apriori", APRIORI_LO, APRIORI_HI, np.nan)]
        for r in QUANTS:
            lo, hi = bounds_from_tp_sizes(tps[(fn, si)], r)
            arms.append((f"insample@{r:.2f}", lo, hi, r))
            # LODO: bound set on every cell whose ROI is NOT in this cell's domain
            other = np.concatenate([v for (f2, s2), v in tps.items()
                                    if f2 not in domains[g.tumor_type]] or [np.zeros(0)])
            lo2, hi2 = bounds_from_tp_sizes(other, r)
            arms.append((f"lodo@{r:.2f}", lo2, hi2, r))

        for name, lo, hi, r in arms:
            rows.append({**base, "arm": name, "lo_um": lo, "hi_um": hi, "retention_nominal": r,
                         **score_arm(cx, cy, size_um, gt_eval, mr, n_mit, lo, hi)})
    return pd.DataFrame(rows)


def paired_table(df: pd.DataFrame, arm: str) -> pd.DataFrame:
    """Per-cell paired deltas against the baseline, for one arm."""
    b = df[df.arm == "baseline"].set_index(["file_name", "seed_index"])
    a = df[df.arm == arm].set_index(["file_name", "seed_index"])
    out = a[["tumor_type", "n_mit", "lo_um", "hi_um", "n_kept", "recall"]].copy()
    out["n_pool"] = b["n_kept"]
    out["frac_kept"] = a["n_kept"] / b["n_kept"]
    out["recall_base"] = b["recall"]
    out["d_recall"] = a["recall"] - b["recall"]
    for q in LEVELS:
        c = f"depth_{int(q*100)}"
        out[f"{c}_base"] = b[c]
        out[c] = a[c]
        out[f"d_{c}"] = a[c] - b[c]
    return out.reset_index()


# --------------------------------------------------------------------------------------
# Fast path: cache the radius query per cell instead of rebuilding it for every arm.
#
# `evaluate.bucket_detections` builds a KDTree over the ground truth and calls
# `query_radius` over every candidate on EVERY call. F7 scores ten arms per cell, and the
# candidate coordinates never change -- only which of them survive. Querying once per cell
# and replaying the same greedy rule over the survivors is the identical computation with
# the query hoisted out of the loop. `verify_fast_path` asserts that equality against
# `ev.bucket_detections` rather than assuming it.
# --------------------------------------------------------------------------------------
from sklearn.neighbors import KDTree                                    # noqa: E402


def cell_neighbours(cx, cy, gt_eval, radius):
    """For every candidate, the indices of ground-truth points within `radius`."""
    gt_xy = gt_eval[["cx", "cy"]].to_numpy(dtype=np.float64)
    det_xy = np.stack([np.asarray(cx, float), np.asarray(cy, float)], axis=1)
    return KDTree(gt_xy).query_radius(det_xy, r=radius), gt_xy, det_xy


def score_arm_fast(neigh, gt_xy, det_xy, gt_cls, keep, n_mit) -> dict:
    """`score_arm` with the radius query supplied. Same greedy rule as `evaluate.greedy_match`:
    walk detections best-first, each claims the nearest ground truth still free."""
    gt_to_det = np.full(len(gt_xy), -1, dtype=int)
    tp_depths, pos = [], 0
    for d in np.flatnonzero(keep):
        pos += 1                                    # 1-based rank among the survivors
        cands = neigh[d]
        if len(cands) == 0:
            continue
        free = [g for g in cands if gt_to_det[g] == -1]
        if not free:
            continue
        dist = np.hypot(gt_xy[free, 0] - det_xy[d, 0], gt_xy[free, 1] - det_xy[d, 1])
        g = free[int(np.argmin(dist))]
        gt_to_det[g] = pos
        if gt_cls[g] == ds.MITOTIC:
            tp_depths.append(pos)
    row = dict(n_kept=int(keep.sum()), n_tp=len(tp_depths),
               recall=len(tp_depths) / n_mit if n_mit else np.nan)
    row.update(depths_from_tp_ranks(np.array(tp_depths, dtype=float), n_mit))
    return row


def verify_fast_path(cx, cy, size_um, gt_eval, radius, n_mit, lo, hi) -> bool:
    """The fast path must agree with `ev.bucket_detections` exactly, or it is not usable."""
    slow = score_arm(cx, cy, size_um, gt_eval, radius, n_mit, lo, hi)
    neigh, gt_xy, det_xy = cell_neighbours(cx, cy, gt_eval, radius)
    gt_cls = gt_eval["category_id"].to_numpy()
    fast = score_arm_fast(neigh, gt_xy, det_xy, gt_cls, keep_mask(size_um, lo, hi), n_mit)
    same = all((np.isnan(slow[k]) and np.isnan(fast[k])) if isinstance(slow[k], float)
               and np.isnan(slow[k]) else slow[k] == fast[k] for k in slow)
    if not same:
        raise AssertionError(f"fast path disagrees:\n  slow={slow}\n  fast={fast}")
    return True
