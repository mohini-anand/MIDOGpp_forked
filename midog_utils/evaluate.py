"""Scoring detections against MIDOG++ ground truth.

Ground truth is a point click, so matching is by centre distance, not IoU. Because
every annotation box is a synthetic 50x50 around that point, IoU would just be a
monotone function of centre distance anyway.

Matching is **greedy in descending detection score**, once, over the full ranked list.
That is deliberate: optimal assignment recomputed at each score threshold can
re-assign earlier detections and produce a non-monotone FROC, where lowering the
threshold *decreases* sensitivity. With one greedy pass every threshold is a truncation
of the ranked list, so the curve is monotone by construction. This is also what
`evalutils.scorers.score_detection` -- the function `evaluation.py` uses -- does.
`optimal_assignment_disagreement` is provided as a one-off cross-check at a fixed
threshold.

The three detection buckets:

``TP``               matched a category-1 GT -- a real mitotic figure found from one click
``FP_lookalike``     matched a category-2 GT -- a structure a pathologist examined and rejected
``FP_unannotated``   matched nothing

Both FP buckets are false positives for the mitosis-detection task. The split is a
diagnostic of *what kind* of mistake the algorithm makes, not a separate metric family.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.neighbors import KDTree

from .dataset import LOOKALIKE, MITOTIC

TP = "TP"
FP_LOOKALIKE = "FP_lookalike"
FP_UNANNOTATED = "FP_unannotated"

# MIDOG's own operating point: evaluation.py calls score_detection(radius=7.5E-3), in
# millimetres. Converted per image via that ROI's microns-per-pixel.
MIDOG_RADIUS_UM = 7.5


def radius_px(mpp: float, radius_um: float = MIDOG_RADIUS_UM) -> float:
    return radius_um / mpp


def greedy_match(det_xy, gt_xy, radius):
    """Match detections (already ranked best-first) to ground truth, one-to-one.

    Returns ``(det_to_gt, gt_to_det)``: for each detection the index of the GT it
    claimed (-1 if none), and for each GT the rank of the detection that claimed it
    (-1 if none).
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


def bucket_detections(detections: pd.DataFrame, gt: pd.DataFrame, radius: float):
    """Annotate a ranked detection frame with its match and bucket.

    ``gt`` must be the *evaluation* ground truth -- i.e. with the seed annotation
    already removed -- and must pool both categories, so one detection cannot claim a
    mitotic figure and a look-alike at the same time.
    """
    det_xy = detections[["cx", "cy"]].to_numpy() if len(detections) else np.zeros((0, 2))
    gt_xy = gt[["cx", "cy"]].to_numpy() if len(gt) else np.zeros((0, 2))
    det_to_gt, gt_to_det = greedy_match(det_xy, gt_xy, radius)

    gt_cls = gt["category_id"].to_numpy() if len(gt) else np.zeros(0, dtype=int)
    gt_ann = gt["ann_id"].to_numpy() if len(gt) else np.zeros(0, dtype=int)

    buckets, matched_ann, matched_cls = [], [], []
    for g in det_to_gt:
        if g < 0:
            buckets.append(FP_UNANNOTATED)
            matched_ann.append(-1)
            matched_cls.append(0)
        else:
            buckets.append(TP if gt_cls[g] == MITOTIC else FP_LOOKALIKE)
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
    """Fraction of mitotic GT recovered by the top-K detections, K = number of mitotic GT.

    This is the headline metric because it needs no threshold. Peak counts above any
    fixed TM_CCOEFF_NORMED value differ by orders of magnitude between a dense domain
    like lymphosarcoma and a sparse one like melanoma, which makes precision at a fixed
    threshold incomparable across domains.
    """
    if n_gt_mitotic == 0:
        return float("nan")
    k = n_gt_mitotic if k is None else k
    top = np.asarray(det_buckets)[:k]
    return float(np.sum(top == TP) / n_gt_mitotic)


def froc(det_buckets, n_gt_mitotic: int, area_mm2: float):
    """Cumulative (false positives per mm^2, sensitivity) over the ranked list."""
    b = np.asarray(det_buckets)
    if n_gt_mitotic == 0 or len(b) == 0:
        return np.zeros(0), np.zeros(0)
    tp = np.cumsum(b == TP)
    fp = np.cumsum(b != TP)  # both FP buckets count against the mitosis task
    return fp / area_mm2, tp / n_gt_mitotic


def sensitivity_at_fp(fp_per_mm2, sensitivity, levels=(1, 2, 4, 8, 16, 32, 64)):
    out = {}
    for lv in levels:
        ok = fp_per_mm2 <= lv
        out[lv] = float(sensitivity[ok].max()) if ok.any() else 0.0
    return out


def threshold_metrics(det: pd.DataFrame, gt: pd.DataFrame, score_threshold: float) -> dict:
    """Precision / recall / F1 and FP composition at one score threshold."""
    sel = det[det["score"] >= score_threshold]
    n_gt_mitotic = int((gt["category_id"] == MITOTIC).sum())
    tp = int((sel["bucket"] == TP).sum())
    fp_look = int((sel["bucket"] == FP_LOOKALIKE).sum())
    fp_un = int((sel["bucket"] == FP_UNANNOTATED).sum())
    fp = fp_look + fp_un
    eps = 1e-9
    return {
        "score_threshold": score_threshold,
        "n_detections": len(sel),
        "tp": tp,
        "fp_lookalike": fp_look,
        "fp_unannotated": fp_un,
        "fn": n_gt_mitotic - tp,
        "precision": tp / (tp + fp + eps),
        "recall": tp / (n_gt_mitotic + eps),
        "f1": 2 * tp / (2 * tp + fp + (n_gt_mitotic - tp) + eps),
        "fp_composition_lookalike": fp_look / (fp + eps),
    }


def lookalike_attraction_rate(gt_out: pd.DataFrame, k: int = None) -> dict:
    """How often the search lands on a pathologist-rejected look-alike.

    Reported *at a detection budget*, and at the same budget as mitotic recall, because
    neither number means anything on its own -- a long enough detection list eventually
    covers every annotation. With ``k`` set, only detections ranked above ``k`` count.

    If attraction and recall are similar at the same budget, the matcher is not
    discriminating mitosis from mimic at all. That comparison is the single most
    informative number in this experiment.
    """
    look = gt_out[gt_out["category_id"] == LOOKALIKE]
    mito = gt_out[gt_out["category_id"] == MITOTIC]

    def _found(df):
        if k is None:
            return df["found"]
        return (df["matched_rank"] >= 0) & (df["matched_rank"] < k)

    suffix = "" if k is None else "_at_k"
    return {
        "n_lookalike_gt": len(look),
        f"lookalike_attraction{suffix}": float(_found(look).mean()) if len(look) else float("nan"),
        "n_mitotic_gt": len(mito),
        f"mitotic_recall{suffix}": float(_found(mito).mean()) if len(mito) else float("nan"),
    }


def topk_composition(det_out: pd.DataFrame, k: int) -> dict:
    """What the top-K detections actually are -- the bucket split at the K budget."""
    top = det_out.head(k)
    return {
        "k": k,
        "topk_tp": int((top["bucket"] == TP).sum()),
        "topk_fp_lookalike": int((top["bucket"] == FP_LOOKALIKE).sum()),
        "topk_fp_unannotated": int((top["bucket"] == FP_UNANNOTATED).sum()),
    }


def recall_by_agreement(gt_out: pd.DataFrame) -> dict:
    """Mitotic recall split by whether the experts agreed, using the `labels` votes."""
    mito = gt_out[gt_out["category_id"] == MITOTIC]
    unan, cont = mito[mito["unanimous"]], mito[~mito["unanimous"]]
    return {
        "n_unanimous": len(unan),
        "recall_unanimous": float(unan["found"].mean()) if len(unan) else float("nan"),
        "n_contested": len(cont),
        "recall_contested": float(cont["found"].mean()) if len(cont) else float("nan"),
    }


def optimal_assignment_disagreement(det: pd.DataFrame, gt: pd.DataFrame, radius: float) -> dict:
    """Cross-check greedy matching against optimal assignment at one operating point.

    Only meaningful at a fixed threshold -- see the module docstring for why optimal
    assignment must not drive the ranked metrics. `linear_sum_assignment` rejects
    `np.inf`, so forbidden pairs get a large finite cost and are filtered afterwards.
    """
    if len(det) == 0 or len(gt) == 0:
        return {"greedy_matches": 0, "optimal_matches": 0, "disagreements": 0}

    d = det[["cx", "cy"]].to_numpy()
    g = gt[["cx", "cy"]].to_numpy()
    cost = np.hypot(d[:, None, 0] - g[None, :, 0], d[:, None, 1] - g[None, :, 1])
    big = radius * 1000.0
    cost_masked = np.where(cost <= radius, cost, big)

    rows, cols = linear_sum_assignment(cost_masked)
    ok = cost[rows, cols] <= radius
    optimal = {(int(r), int(c)) for r, c in zip(rows[ok], cols[ok])}

    det_to_gt, _ = greedy_match(d, g, radius)
    greedy = {(i, int(gi)) for i, gi in enumerate(det_to_gt) if gi >= 0}

    return {
        "greedy_matches": len(greedy),
        "optimal_matches": len(optimal),
        "disagreements": len(greedy ^ optimal) // 2 if greedy != optimal else 0,
    }


def evaluate_run(detections, gt_eval, radius, area_mm2, score_threshold=None):
    """Everything above, bundled. ``gt_eval`` must already exclude the seed."""
    det_out, gt_out = bucket_detections(detections, gt_eval, radius)
    n_mit = int((gt_eval["category_id"] == MITOTIC).sum())

    fp_mm2, sens = froc(det_out["bucket"].to_numpy(), n_mit, area_mm2)
    metrics = {
        "match_radius_px": round(float(radius), 1),
        "n_gt_mitotic_eval": n_mit,
        "recall_at_k": recall_at_k(det_out["bucket"].to_numpy(), n_mit),
        "roi_area_mm2": round(float(area_mm2), 2),
    }
    metrics.update({f"sens@{k}fp_mm2": v for k, v in sensitivity_at_fp(fp_mm2, sens).items()})
    metrics.update(topk_composition(det_out, n_mit))
    # Both at the same K budget, so they are directly comparable, plus the full-list
    # versions for reference.
    metrics.update(lookalike_attraction_rate(gt_out, k=n_mit))
    metrics.update(lookalike_attraction_rate(gt_out))
    metrics["n_detections_total"] = len(det_out)
    metrics.update(recall_by_agreement(gt_out))
    metrics.update(full_list_breakdown(det_out, gt_out))
    if score_threshold is not None:
        metrics.update(threshold_metrics(det_out, gt_eval, score_threshold))
    return det_out, gt_out, metrics, (fp_mm2, sens)


def full_list_breakdown(det_out: pd.DataFrame, gt_out: pd.DataFrame) -> dict:
    """Bucket counts over the *entire* detection list, with no top-K truncation.

    `recall@K` answers "given a budget, how good are the best guesses". This answers the
    complementary question: of everything the pipeline flagged, how much landed on an
    annotation at all?

    **Caveat, measured.** At a low score floor the detection list saturates the ROI --
    at 0.25 roughly 91-97% of each image lies within a match radius of some detection --
    so ``all_mitotic_gt_found`` becomes a coverage statistic rather than a detection
    result (uniform random points reach full-list recall 1.000 on two of the seven ROIs).
    ``all_precision_mitotic`` is likewise a property of the score floor, which swings the
    list length by three orders of magnitude between images. Use this alongside
    ``threshold_sweep`` and against the random control, not on its own.
    """
    n = len(det_out)
    tp = int((det_out["bucket"] == TP).sum())
    look = int((det_out["bucket"] == FP_LOOKALIKE).sum())
    un = int((det_out["bucket"] == FP_UNANNOTATED).sum())
    mit_gt = gt_out[gt_out["category_id"] == MITOTIC]
    look_gt = gt_out[gt_out["category_id"] == LOOKALIKE]
    eps = 1e-9
    return {
        "all_n_detections": n,
        "all_tp": tp,
        "all_fp_lookalike": look,
        "all_fp_unannotated": un,
        "all_matched_any_gt": tp + look,
        "all_matched_frac": (tp + look) / (n + eps),
        "all_precision_mitotic": tp / (n + eps),
        "all_mitotic_gt_found": int(mit_gt["found"].sum()),
        "all_mitotic_gt_missed": int((~mit_gt["found"]).sum()),
        "all_lookalike_gt_found": int(look_gt["found"].sum()),
        "all_lookalike_gt_missed": int((~look_gt["found"]).sum()),
    }


def threshold_sweep(det_out: pd.DataFrame, gt_eval: pd.DataFrame,
                    thresholds=(0.25, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) -> pd.DataFrame:
    """Bucket counts and P/R/F1 as the score cutoff moves.

    Included because "all the results that were found" is only well defined relative to the
    score floor the search ran at. This makes that dependence explicit instead of hiding it
    in a config value.
    """
    n_mit = int((gt_eval["category_id"] == MITOTIC).sum())
    rows = []
    for t in thresholds:
        sel = det_out[det_out["score"] >= t]
        tp = int((sel["bucket"] == TP).sum())
        look = int((sel["bucket"] == FP_LOOKALIKE).sum())
        un = int((sel["bucket"] == FP_UNANNOTATED).sum())
        eps = 1e-9
        rows.append({
            "score_threshold": t, "n_detections": len(sel), "tp": tp,
            "fp_lookalike": look, "fp_unannotated": un,
            "precision": round(tp / (len(sel) + eps), 4),
            "recall": round(tp / (n_mit + eps), 4),
            "f1": round(2 * tp / (2 * tp + look + un + (n_mit - tp) + eps), 4),
        })
    return pd.DataFrame(rows)
