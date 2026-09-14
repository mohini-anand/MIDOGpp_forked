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

``human_correct_label``   matched a category-1 GT -- a real mitotic figure found from one click
``human_rejected_label``  matched a category-2 GT -- a structure a pathologist examined and rejected
``non_human_findings``    matched nothing

Both of the latter two are false positives for the mitosis-detection task. The split is
a diagnostic of *what kind* of mistake the algorithm makes, not a separate metric
family. (Renamed from the original `TP`/`FP_lookalike`/`FP_unannotated` -- see
`Research Logs/design_choices.md`, section 5; the underlying TP/FP counting logic below
is unchanged, only these three label strings and their names.)

**Read ``coverage_frac`` before reading any full-list number.** A detection list long
enough to tile the ROI answers "is this annotation within the match radius of some
detection?" by geometry rather than by evidence, and every metric here that is not
budgeted at K is exactly that question. ``coverage_frac`` measures the saturation
directly so it cannot be mistaken for a result.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.neighbors import KDTree

from .dataset import LOOKALIKE, MITOTIC

HUMAN_CORRECT_LABEL = "human_correct_label"
HUMAN_REJECTED_LABEL = "human_rejected_label"
NON_HUMAN_FINDINGS = "non_human_findings"

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


def coverage_fraction(det_xy, roi_shape, radius: float, stride: int = 16) -> float:
    """Fraction of the ROI lying within ``radius`` of *some* detection.

    This is the sanity check on every un-budgeted metric in this module. "Was this
    annotation found?" is scored as "is there a detection within ``radius`` of it?", so
    if a large fraction of *arbitrary* ROI locations already satisfy that, a full-list
    recall near 1.0 says nothing about the detector.

    At the 0.25 score floor this run uses, the answer is 0.91-0.97 -- the detection list
    tiles the ROI more finely than the metric can resolve. Sampled on a strided grid;
    stride 16 gives ~150k probe points on a 39-megapixel ROI, far more than enough for
    two decimal places.
    """
    det_xy = np.asarray(det_xy, dtype=np.float64).reshape(-1, 2)
    if len(det_xy) == 0:
        return 0.0
    h, w = int(roi_shape[0]), int(roi_shape[1])
    ys, xs = np.mgrid[0:h:stride, 0:w:stride]
    probes = np.stack([xs.ravel(), ys.ravel()], axis=1).astype(np.float64)
    dist, _ = KDTree(det_xy).query(probes, k=1)
    return float((dist[:, 0] <= radius).mean())


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
    return float(np.sum(top == HUMAN_CORRECT_LABEL) / n_gt_mitotic)


def froc(det_buckets, n_gt_mitotic: int, area_mm2: float):
    """Cumulative (false positives per mm^2, sensitivity) over the ranked list."""
    b = np.asarray(det_buckets)
    if n_gt_mitotic == 0 or len(b) == 0:
        return np.zeros(0), np.zeros(0)
    tp = np.cumsum(b == HUMAN_CORRECT_LABEL)
    fp = np.cumsum(b != HUMAN_CORRECT_LABEL)  # both FP buckets count against the mitosis task
    return fp / area_mm2, tp / n_gt_mitotic


def sensitivity_at_fp(fp_per_mm2, sensitivity, levels=(1, 2, 4, 8, 16, 32, 64)):
    out = {}
    for lv in levels:
        ok = fp_per_mm2 <= lv
        out[lv] = float(sensitivity[ok].max()) if ok.any() else 0.0
    return out


def _found_within(df: pd.DataFrame, k: int = None):
    """Boolean 'this GT was claimed', optionally restricted to the top-K detections."""
    if k is None:
        return df["found"]
    return (df["matched_rank"] >= 0) & (df["matched_rank"] < k)


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
    """Mitotic recall split by whether the experts agreed, using the `labels` votes.

    Takes the same ``k`` budget as `lookalike_attraction_rate`, and for the same reason.
    Without it this reads the full detection list, where coverage saturation forces both
    halves of the split to 1.0 and the comparison -- "did we miss an obvious mitosis or a
    borderline one?" -- carries no information at all.
    """
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


def optimal_assignment_disagreement(det: pd.DataFrame, gt: pd.DataFrame, radius: float) -> dict:
    """Cross-check greedy matching against optimal assignment at one operating point.

    Only meaningful at a fixed threshold -- see the module docstring for why optimal
    assignment must not drive the ranked metrics. `linear_sum_assignment` rejects
    `np.inf`, so forbidden pairs get a large finite cost and are filtered afterwards.

    ``greedy_only`` / ``optimal_only`` are the (detection, GT) pairs each scheme makes
    and the other does not, counted directly. An earlier version reported
    ``len(greedy ^ optimal) // 2``, which is only correct when every disagreement is a
    swap between two pairs: a detection matched by one scheme and unmatched by the other
    contributes 1 to the symmetric difference and was floored away to 0.
    """
    empty = {"greedy_matches": 0, "optimal_matches": 0, "greedy_only": 0,
             "optimal_only": 0, "n_detections_differing": 0}
    if len(det) == 0 or len(gt) == 0:
        return empty

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
        "greedy_only": len(greedy - optimal),
        "optimal_only": len(optimal - greedy),
        "n_detections_differing": len({i for i, _ in greedy ^ optimal}),
    }


def full_list_breakdown(det_out: pd.DataFrame, gt_out: pd.DataFrame) -> dict:
    """Bucket counts over the *entire* detection list, with no top-K truncation.

    `recall@K` answers "given a budget, how good are the best guesses". This answers the
    complementary question: of everything the pipeline flagged, how much landed on an
    annotation at all?

    These numbers are **not** "the cost of the method at its natural operating point".
    There is no natural operating point -- they are the cost of whatever
    `FSConfig.score_threshold` happens to be, and what that floor means varies by orders
    of magnitude between images (at 0.70: 1 detection on 301.tiff, 2601 on 405.tiff).
    Read them next to `coverage_frac`, and next to the `random_in_tissue` row at the same
    list length, which is the floor they have to beat.
    """
    n = len(det_out)
    tp = int((det_out["bucket"] == HUMAN_CORRECT_LABEL).sum())
    look = int((det_out["bucket"] == HUMAN_REJECTED_LABEL).sum())
    un = int((det_out["bucket"] == NON_HUMAN_FINDINGS).sum())
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


def evaluate_run(detections, gt_eval, radius, area_mm2, roi_shape=None):
    """Everything above, bundled. ``gt_eval`` must already exclude the seed.

    ``roi_shape`` is ``(h, w)``; supply it to get ``coverage_frac``, without which the
    un-budgeted metrics below cannot be interpreted.
    """
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
    metrics.update(recall_by_agreement(gt_out, k=n_mit))
    metrics.update(recall_by_agreement(gt_out))
    metrics.update(full_list_breakdown(det_out, gt_out))
    metrics["coverage_frac"] = (
        round(coverage_fraction(det_out[["cx", "cy"]].to_numpy(), roi_shape, radius), 4)
        if roi_shape is not None else float("nan")
    )
    return det_out, gt_out, metrics, (fp_mm2, sens)


def threshold_sweep(det_out: pd.DataFrame, gt_eval: pd.DataFrame,
                    thresholds=(0.5, 0.6, 0.7, 0.8, 0.9)) -> pd.DataFrame:
    """Bucket counts and P/R/F1 as the score cutoff moves.

    Included because "all the results that were found" is only well defined relative to the
    score floor the search ran at. This makes that dependence explicit instead of hiding it
    in a config value.

    ``thresholds`` starts at 0.5, matching `FSConfig.score_threshold`'s new floor -- values
    below it are below the score the search itself now stops reporting at, so sweeping them
    here would describe detections `find_and_suppress` no longer returns. See
    `Research Logs/design_choices.md`, section 6.
    """
    n_mit = int((gt_eval["category_id"] == MITOTIC).sum())
    rows = []
    for t in thresholds:
        sel = det_out[det_out["score"] >= t]
        tp = int((sel["bucket"] == HUMAN_CORRECT_LABEL).sum())
        look = int((sel["bucket"] == HUMAN_REJECTED_LABEL).sum())
        un = int((sel["bucket"] == NON_HUMAN_FINDINGS).sum())
        eps = 1e-9
        rows.append({
            "score_threshold": t, "n_detections": len(sel), "tp": tp,
            "fp_lookalike": look, "fp_unannotated": un,
            "precision": round(tp / (len(sel) + eps), 4),
            "recall": round(tp / (n_mit + eps), 4),
            "f1": round(2 * tp / (2 * tp + look + un + (n_mit - tp) + eps), 4),
        })
    return pd.DataFrame(rows)
