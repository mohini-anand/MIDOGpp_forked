"""AUC(mitosis vs look-alike) measured without the detection selection bias.

The AUC table in `Research Logs/2026-08-31-chromatin-density-rerank.md` is computed on the
*pipeline's own* candidate set, which contains only a minority of the annotated look-alikes
-- and which set of look-alikes it contains is decided by the same chromatin statistic being
scored. That is a selection correlated with the feature under test, and it inflates the
estimate.

This re-measures on the `nucleus_blobs` candidate set, which captures very nearly all of
both classes, so the comparison is over the full annotated population rather than the
subset the pipeline happened to find. Two statistics are compared on that identical set:

* `chromatin.chromatin_density` -- mean of the darkest 10% of a 51 px window;
* the mean of the Otsu component `nucleus_blobs` itself selected on.

Writes results/lookalike_auc_unbiased.csv.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from midog_utils import baselines as bl
from midog_utils import chromatin as cm
from midog_utils import dataset as ds
from midog_utils import evaluate as ev


ROIS = ["301.tiff", "246.tiff"]  # the only two ROIs with enough of both classes


def auc_and_ci(pos: np.ndarray, neg: np.ndarray):
    """Mann-Whitney AUC with a Hanley-McNeil standard error and normal 95% CI."""
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    n1, n2 = len(pos), len(neg)
    if not n1 or not n2:
        return float("nan"), float("nan"), float("nan"), n1, n2
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    a = (gt + 0.5 * eq) / (n1 * n2)
    q1 = a / (2 - a)
    q2 = 2 * a * a / (1 + a)
    se = np.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n2 - 1) * (q2 - a * a)) / (n1 * n2))
    return float(a), float(a - 1.96 * se), float(a + 1.96 * se), n1, n2


def main():
    _images, ann = ds.load_annotations()
    rows = []
    for fn in ROIS:
        path = f"images/{fn}"
        rgb = ds.load_roi(path)
        radius = ev.radius_px(ds.roi_mpp(path))
        gt = ds.image_annotations(ann, fn)
        hem_od = cm.hematoxylin_od(rgb)          # unclipped, the ranking channel

        # nucleus_blobs takes RGB and does its own hematoxylin conversion internally
        blobs = bl.nucleus_blobs(rgb)
        det, _gt_out = ev.bucket_detections(
            blobs.reset_index(drop=True), gt.reset_index(drop=True), radius)

        # chromatin density at each blob centre, and the blob's own component mean
        det = cm.score_detections(det, hem_od)
        native = det["score"].to_numpy() if "score" in det.columns else np.full(len(det), np.nan)

        b = det["bucket"].to_numpy()
        tp_mask = b == ev.HUMAN_CORRECT_LABEL
        lk_mask = b == ev.HUMAN_REJECTED_LABEL
        un_mask = b == ev.NON_HUMAN_FINDINGS

        for stat_name, v in (("chromatin_density", det["od"].to_numpy()),
                             ("blob_component_score", native)):
            for contrast, pos, neg in (
                ("tp_vs_lookalike", v[tp_mask], v[lk_mask]),
                ("tp_vs_unannotated", v[tp_mask], v[un_mask]),
                ("lookalike_vs_unannotated", v[lk_mask], v[un_mask]),
            ):
                a, lo, hi, n1, n2 = auc_and_ci(pos, neg)
                rows.append(dict(file_name=fn, statistic=stat_name, contrast=contrast,
                                 auc=round(a, 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                                 n_pos=n1, n_neg=n2))

        # capture rates -- the whole point of using this candidate set
        rows.append(dict(file_name=fn, statistic="_capture", contrast="mitotic_captured",
                         auc=round(tp_mask.sum() / max((gt.category_id == ds.MITOTIC).sum(), 1), 3),
                         ci_lo=np.nan, ci_hi=np.nan,
                         n_pos=int(tp_mask.sum()), n_neg=int((gt.category_id == ds.MITOTIC).sum())))
        rows.append(dict(file_name=fn, statistic="_capture", contrast="lookalike_captured",
                         auc=round(lk_mask.sum() / max((gt.category_id == ds.LOOKALIKE).sum(), 1), 3),
                         ci_lo=np.nan, ci_hi=np.nan,
                         n_pos=int(lk_mask.sum()), n_neg=int((gt.category_id == ds.LOOKALIKE).sum())))
        print(f"{fn}: {len(blobs)} blobs, captured {tp_mask.sum()} mitotic / "
              f"{lk_mask.sum()} look-alike", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv("results/lookalike_auc_unbiased.csv", index=False)
    pd.set_option("display.width", 200)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
