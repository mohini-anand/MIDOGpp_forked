"""
    Scoring detections against MIDOG++ ground truth by centre distance, greedily matched
    best-detection-first so a top-K prefix's matches never depend on lower-ranked
    detections. Three buckets: a detection matches a mitotic figure (human_correct_label),
    a look-alike (human_rejected_label), or nothing (non_human_findings).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.neighbors import KDTree

from .dataset import LOOKALIKE, MITOTIC

HUMAN_CORRECT_LABEL = "human_correct_label"
HUMAN_REJECTED_LABEL = "human_rejected_label"
NON_HUMAN_FINDINGS = "non_human_findings"

MIDOG_RADIUS_UM = 7.5  # MIDOG's own operating point


def radius_px(mpp: float, radius_um: float = MIDOG_RADIUS_UM) -> float:
    return radius_um / mpp


def greedy_match(det_xy, gt_xy, radius):
    """
        Match detections (already ranked best-first) to ground truth, one-to-one.

        Returns tuple[np.ndarray, np.ndarray]: (det_to_gt, gt_to_det) -- for each detection
        the GT index it claimed (-1 if none), and for each GT the rank that claimed it.
    """
    det_xy = np.asarray(det_xy, dtype=np.float64).reshape(-1, 2)
    gt_xy = np.asarray(gt_xy, dtype=np.float64).reshape(-1, 2)

    det_to_gt = np.full(len(det_xy), -1, dtype=int)
    gt_to_det = np.full(len(gt_xy), -1, dtype=int)
    if len(det_xy) == 0 or len(gt_xy) == 0:
        return det_to_gt, gt_to_det

    tree = KDTree(gt_xy)
    neighbours = tree.query_radius(det_xy, r=radius)
    for d_idx, cands in enumerate(neighbours):
        free = [g for g in cands if gt_to_det[g] == -1]
        if not free:
            continue
        d = np.hypot(gt_xy[free, 0] - det_xy[d_idx, 0], gt_xy[free, 1] - det_xy[d_idx, 1])
        g = free[int(np.argmin(d))]
        det_to_gt[d_idx] = g
        gt_to_det[g] = d_idx
    return det_to_gt, gt_to_det


def coverage_fraction(det_xy, roi_shape, radius: float, stride: int = 16) -> float:
    """Fraction of the ROI lying within ``radius`` of some detection -- the saturation
    check for any un-budgeted recall number."""
    det_xy = np.asarray(det_xy, dtype=np.float64).reshape(-1, 2)
    if len(det_xy) == 0:
        return 0.0
    h, w = int(roi_shape[0]), int(roi_shape[1])
    ys, xs = np.mgrid[0:h:stride, 0:w:stride]
    probes = np.stack([xs.ravel(), ys.ravel()], axis=1).astype(np.float64)
    dist, _ = KDTree(det_xy).query(probes, k=1)
    return float((dist[:, 0] <= radius).mean())


def bucket_detections(detections: pd.DataFrame, gt: pd.DataFrame, radius: float):
    """
        Annotate a ranked detection frame with its match and bucket.

        detections (pd.DataFrame): ranked detections with cx, cy columns.
        gt (pd.DataFrame): evaluation ground truth (seed annotation already excluded).
        radius (float): match radius in pixels.

        Returns tuple[pd.DataFrame, pd.DataFrame]: (detections, gt) each with match columns added.
    """
    det_xy = detections[["cx", "cy"]].to_numpy() if len(detections) else np.zeros((0, 2))
    gt_xy = gt[["cx", "cy"]].to_numpy() if len(gt) else np.zeros((0, 2))
    det_to_gt, gt_to_det = greedy_match(det_xy, gt_xy, radius)

    gt_cls = gt["category_id"].to_numpy() if len(gt) else np.zeros(0, dtype=int)
    gt_ann = gt["ann_id"].to_numpy() if len(gt) else np.zeros(0, dtype=int)

    buckets, matched_ann, matched_cls = [], [], []
    for g in det_to_gt:
        if g < 0:
            buckets.append(NON_HUMAN_FINDINGS)
            matched_ann.append(-1)
            matched_cls.append(0)
        else:
            buckets.append(HUMAN_CORRECT_LABEL if gt_cls[g] == MITOTIC else HUMAN_REJECTED_LABEL)
            matched_ann.append(int(gt_ann[g]))
            matched_cls.append(int(gt_cls[g]))

    out = detections.copy()
    out["bucket"] = buckets
    out["matched_ann_id"] = matched_ann
    out["matched_category"] = matched_cls

    gt_out = gt.copy()
    gt_out["matched_rank"] = gt_to_det
    gt_out["found"] = gt_to_det >= 0
    return out, gt_out


def recall_at_k(det_buckets, n_gt_mitotic: int, k: int = None) -> float:
    """Fraction of mitotic GT recovered by the top-K detections, K = number of mitotic GT."""
    if n_gt_mitotic == 0:
        return float("nan")
    k = n_gt_mitotic if k is None else k
    top = np.asarray(det_buckets)[:k]
    return float(np.sum(top == HUMAN_CORRECT_LABEL) / n_gt_mitotic)


def _found_within(df: pd.DataFrame, k: int = None):
    """Boolean 'this GT was claimed', optionally restricted to the top-K detections."""
    if k is None:
        return df["found"]
    return (df["matched_rank"] >= 0) & (df["matched_rank"] < k)


def lookalike_attraction_rate(gt_out: pd.DataFrame, k: int = None) -> dict:
    """How often the search lands on a pathologist-rejected look-alike, at the same
    budget as mitotic recall."""
    look = gt_out[gt_out["category_id"] == LOOKALIKE]
    mito = gt_out[gt_out["category_id"] == MITOTIC]

    suffix = "" if k is None else "_at_k"
    return {
        "n_lookalike_gt": len(look),
        f"lookalike_attraction{suffix}": float(_found_within(look, k).mean()) if len(look) else float("nan"),
        "n_mitotic_gt": len(mito),
        f"mitotic_recall{suffix}": float(_found_within(mito, k).mean()) if len(mito) else float("nan"),
    }


def topk_composition(det_out: pd.DataFrame, k: int) -> dict:
    """What the top-K detections actually are -- the bucket split at the K budget."""
    top = det_out.head(k)
    return {
        "k": k,
        "topk_tp": int((top["bucket"] == HUMAN_CORRECT_LABEL).sum()),
        "topk_fp_lookalike": int((top["bucket"] == HUMAN_REJECTED_LABEL).sum()),
        "topk_fp_unannotated": int((top["bucket"] == NON_HUMAN_FINDINGS).sum()),
    }


def recall_by_agreement(gt_out: pd.DataFrame, k: int = None) -> dict:
    """Mitotic recall split by whether the experts agreed (``unanimous`` votes)."""
    mito = gt_out[gt_out["category_id"] == MITOTIC]
    unan, cont = mito[mito["unanimous"]], mito[~mito["unanimous"]]
    found_unan, found_cont = _found_within(unan, k), _found_within(cont, k)
    suffix = "" if k is None else "_at_k"
    return {
        "n_unanimous": len(unan),
        f"recall_unanimous{suffix}": float(found_unan.mean()) if len(unan) else float("nan"),
        "n_contested": len(cont),
        f"recall_contested{suffix}": float(found_cont.mean()) if len(cont) else float("nan"),
    }


def evaluate_run(detections, gt_eval, radius, area_mm2, roi_shape=None):
    """Everything above, bundled. ``gt_eval`` must already exclude the seed. Pass
    ``roi_shape`` (h, w) to also get ``coverage_frac``."""
    det_out, gt_out = bucket_detections(detections, gt_eval, radius)
    n_mit = int((gt_eval["category_id"] == MITOTIC).sum())

    metrics = {
        "match_radius_px": round(float(radius), 1),
        "n_gt_mitotic_eval": n_mit,
        "recall_at_k": recall_at_k(det_out["bucket"].to_numpy(), n_mit),
        "roi_area_mm2": round(float(area_mm2), 2),
    }
    metrics.update(topk_composition(det_out, n_mit))
    metrics.update(lookalike_attraction_rate(gt_out, k=n_mit))
    metrics.update(lookalike_attraction_rate(gt_out))
    metrics["n_detections_total"] = len(det_out)
    metrics.update(recall_by_agreement(gt_out, k=n_mit))
    metrics.update(recall_by_agreement(gt_out))
    metrics["coverage_frac"] = (
        round(coverage_fraction(det_out[["cx", "cy"]].to_numpy(), roi_shape, radius), 4)
        if roi_shape is not None else float("nan")
    )
    return det_out, gt_out, metrics
