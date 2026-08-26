"""Driver for the single-pass, one-seed-per-image experiment.

Scope of this first run, deliberately narrow: **one image per tumour domain, one
randomly chosen mitotic-figure seed each, a single forward pass** -- no iterative
bootstrapping, no multi-seed variance estimate, no hematoxylin variant. Those are
wired up but not run here.

Numbers from a single seed per image are illustrative, not a domain ranking. With one
50x50 template the run-to-run spread from seed choice is expected to be large; treat
the per-domain ordering as provisional until the 5-seed sweep is run.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import baselines as bl
from . import channels as ch
from . import dataset as ds
from . import evaluate as ev
from . import find_and_suppress as fs
from . import template_match as tm


def select_domain_images(images: pd.DataFrame, annotations: pd.DataFrame, images_dir="images"):
    """One ROI per tumour domain: the downloaded image with the most mitotic figures.

    Maximising category-1 count per image matters at this sample size -- with one seed
    and one pass, an image with 4 mitotic figures gives an evaluation set of 3. It also
    rules out 001.tiff automatically, which has zero category-1 annotations and so
    cannot be seeded with a mitotic figure at all. Ties break to the lower image id.
    """
    on_disk = {p.name for p in Path(images_dir).glob("*.tiff")}
    counts = (
        annotations[annotations["category_id"] == ds.MITOTIC]
        .groupby("file_name").size().rename("n_mitotic")
    )
    df = images[images["file_name"].isin(on_disk)].merge(
        counts, left_on="file_name", right_index=True, how="left"
    )
    df["n_mitotic"] = df["n_mitotic"].fillna(0).astype(int)
    df = df[df["n_mitotic"] > 0]
    df = df.sort_values(["tumor_type", "n_mitotic", "image_id"], ascending=[True, False, True])
    return df.groupby("tumor_type", as_index=False).first().reset_index(drop=True)


def pick_seed(gt_mitotic: pd.DataFrame, rng, border: int, roi_shape) -> pd.Series:
    """A random mitotic annotation, far enough from the border to build a full template.

    Annotations within half a padded-patch of the edge cannot yield a rotation-safe
    template; they are excluded here rather than padded, since padding would feed
    fabricated pixels into the template.
    """
    h, w = roi_shape[:2]
    ok = gt_mitotic[
        (gt_mitotic["cx"] >= border) & (gt_mitotic["cx"] < w - border)
        & (gt_mitotic["cy"] >= border) & (gt_mitotic["cy"] < h - border)
    ]
    if len(ok) == 0:
        raise ValueError("no mitotic annotation far enough from the ROI border to seed with")
    return ok.iloc[int(rng.integers(len(ok)))]


def run_one_image(
    file_name: str,
    annotations: pd.DataFrame,
    images_dir="images",
    cfg: fs.FSConfig = None,
    rng=None,
    run_baselines=True,
) -> dict:
    """Full pipeline + baselines for one ROI. Returns everything needed to plot or re-score."""
    cfg = cfg or fs.FSConfig()
    rng = np.random.default_rng(0) if rng is None else rng
    path = f"{images_dir}/{file_name}"

    rgb = ds.load_roi(path)
    mpp = ds.roi_mpp(path)
    area_mm2 = ds.check_roi_scale(path, rgb.shape)
    radius = ev.radius_px(mpp)

    gt = ds.image_annotations(annotations, file_name)
    seed = pick_seed(gt[gt["category_id"] == ds.MITOTIC], rng, cfg.patch_size // 2, rgb.shape)
    gt_eval = gt[gt["ann_id"] != seed["ann_id"]].reset_index(drop=True)

    img = ch.to_channel(rgb, cfg.channel)
    t0 = time.time()
    det, info = fs.find_and_suppress(img, (seed["cx"], seed["cy"]), cfg)
    info["t_total_s"] = round(time.time() - t0, 1)

    det_out, gt_out, metrics, curve = ev.evaluate_run(
        det, gt_eval, radius, area_mm2, score_threshold=None
    )
    metrics = {"file_name": file_name, "method": "find_and_suppress",
               "seed_ann_id": int(seed["ann_id"]), "channel": cfg.channel,
               "mpp": round(mpp, 4), **metrics, **info}

    out = {
        "file_name": file_name,
        "rgb": rgb,
        "seed": seed,
        "gt": gt,
        "gt_eval": gt_eval,
        "radius": radius,
        "area_mm2": area_mm2,
        "detections": det_out,
        "gt_scored": gt_out,
        "curve": curve,
        "metrics": [metrics],
    }

    if run_baselines:
        mask = bl.tissue_mask(rgb)
        out["tissue_fraction"] = float(mask.mean())

        blobs = bl.nucleus_blobs(rgb, mask)
        b_det, b_gt, b_metrics, b_curve = ev.evaluate_run(blobs, gt_eval, radius, area_mm2)

        # `recall@K` already compares at an equal budget, but the full-list recall does
        # not -- the blob detector emits thousands more candidates than the matcher. Score
        # it again truncated to the matcher's list length so the two are comparable at a
        # second, much looser operating point as well.
        capped = blobs.head(len(det_out)).reset_index(drop=True)
        capped["rank"] = np.arange(len(capped))
        _, _, capped_metrics, _ = ev.evaluate_run(capped, gt_eval, radius, area_mm2)
        b_metrics["mitotic_recall_at_matcher_budget"] = capped_metrics["mitotic_recall"]
        b_metrics["matcher_budget"] = len(det_out)

        out["baseline_blobs"] = b_det
        out["baseline_blobs_curve"] = b_curve
        out["metrics"].append(
            {"file_name": file_name, "method": "nucleus_blobs", "seed_ann_id": -1,
             "channel": "hematoxylin", "mpp": round(mpp, 4), **b_metrics}
        )

        rnd = bl.random_in_tissue(mask, n=len(det_out), rng=rng)
        r_det, r_gt, r_metrics, r_curve = ev.evaluate_run(rnd, gt_eval, radius, area_mm2)
        out["baseline_random"] = r_det
        out["baseline_random_curve"] = r_curve
        out["metrics"].append(
            {"file_name": file_name, "method": "random_in_tissue", "seed_ann_id": -1,
             "channel": "-", "mpp": round(mpp, 4), **r_metrics}
        )

    return out


def run_experiment(
    selection: pd.DataFrame,
    annotations: pd.DataFrame,
    images_dir="images",
    cfg: fs.FSConfig = None,
    seed: int = 0,
    run_baselines=True,
    keep_images=False,
    verbose=True,
):
    """Loop `run_one_image` over the selected ROIs.

    ``keep_images=False`` drops the decoded RGB arrays from the returned results, since
    seven 39-megapixel ROIs will not comfortably sit in a notebook's memory at once.
    """
    cfg = cfg or fs.FSConfig()
    results, metric_rows, detection_rows = {}, [], []

    for _, row in selection.iterrows():
        fn = row["file_name"]
        if verbose:
            print(f"=== {fn}  ({row['tumor_type']}, {row['n_mitotic']} mitotic) ", flush=True)
        t0 = time.time()
        res = run_one_image(
            fn, annotations, images_dir, cfg,
            rng=np.random.default_rng(seed), run_baselines=run_baselines,
        )
        for m in res["metrics"]:
            m["tumor_type"] = row["tumor_type"]
            metric_rows.append(m)

        d = res["detections"].copy()
        d.insert(0, "file_name", fn)
        d.insert(1, "method", "find_and_suppress")
        detection_rows.append(d)

        if not keep_images:
            res.pop("rgb", None)
        results[fn] = res
        if verbose:
            m = res["metrics"][0]
            print(
                f"    recall@K={m['recall_at_k']:.3f}  "
                f"attraction@K={m['lookalike_attraction_at_k']:.3f}  "
                f"({time.time() - t0:.0f}s)",
                flush=True,
            )

    metrics = pd.DataFrame(metric_rows)
    detections = pd.concat(detection_rows, ignore_index=True) if detection_rows else pd.DataFrame()
    return results, metrics, detections


SUMMARY_COLUMNS = [
    "tumor_type", "file_name", "method", "n_gt_mitotic_eval", "n_lookalike_gt",
    "recall_at_k", "lookalike_attraction_at_k", "topk_tp", "topk_fp_lookalike",
    "topk_fp_unannotated", "sens@8fp_mm2", "sens@64fp_mm2",
    "mitotic_recall", "lookalike_attraction", "n_detections_total",
    "recall_unanimous", "recall_contested",
]

FULL_LIST_COLUMNS = [
    "tumor_type", "file_name", "all_n_detections", "all_tp", "all_fp_lookalike",
    "all_fp_unannotated", "all_matched_any_gt", "all_matched_frac", "all_precision_mitotic",
    "all_mitotic_gt_found", "all_mitotic_gt_missed",
    "all_lookalike_gt_found", "all_lookalike_gt_missed",
]


def full_list_table(metrics: pd.DataFrame, method="find_and_suppress") -> pd.DataFrame:
    """Every detection the search returned, bucketed -- no top-K budget applied."""
    cols = [c for c in FULL_LIST_COLUMNS if c in metrics.columns]
    sel = metrics if method is None else metrics[metrics["method"] == method]
    return sel[cols].sort_values("tumor_type").reset_index(drop=True)


def summary_table(metrics: pd.DataFrame, method="find_and_suppress") -> pd.DataFrame:
    cols = [c for c in SUMMARY_COLUMNS if c in metrics.columns]
    sel = metrics if method is None else metrics[metrics["method"] == method]
    return sel[cols].sort_values("tumor_type").reset_index(drop=True)


def scale_usage(metrics: pd.DataFrame, detections: pd.DataFrame, scales=(0.6, 0.8, 1.0)) -> pd.DataFrame:
    """Which template scale actually won each detection, versus uniform expectation.

    This is a check on `fused_response`, not on the biology. `TM_CCOEFF_NORMED`
    normalises by the template's pixel count, so a smaller template has a
    higher-variance null and therefore a systematically higher *maximum* by chance. A
    raw element-wise max across sizes is biased toward the smallest scale, and if the
    smallest scale dominates then the pipeline is effectively matching a small chromatin
    patch rather than the annotated object -- which would confound any claim about *why*
    template matching underperforms.

    Compare each column against 1/len(scales). `template_match.fused_response(...,
    scale_normalize=True)` is the corrective.
    """
    k_by_file = {
        r["file_name"]: int(r["n_gt_mitotic_eval"])
        for _, r in metrics[metrics["method"] == "find_and_suppress"].iterrows()
    }
    rows = []
    for fn, k in k_by_file.items():
        d = detections[detections["file_name"] == fn].sort_values("rank")
        top = d.head(k)["scale"].value_counts(normalize=True)
        allv = d["scale"].value_counts(normalize=True)
        row = {"file_name": fn, "K": k}
        row.update({f"topK_scale_{s}": round(float(top.get(s, 0.0)), 3) for s in scales})
        row.update({f"all_scale_{s}": round(float(allv.get(s, 0.0)), 3) for s in scales})
        rows.append(row)
    return pd.DataFrame(rows).sort_values("K", ascending=False).reset_index(drop=True)


def score_probe(file_name, annotations, seed_ann_id, images_dir="images", cfg=None,
                n_random=3000, pad=5, rng=None):
    """Does the response map score mitotic figures above ordinary nuclei at all?

    Detection metrics conflate two questions -- whether the map carries signal, and whether
    the ranking can surface it against the base rate. This separates them by reading the
    fused map directly at annotated locations, at detected nucleus centroids, and at random
    tissue points.

    ``auc_mitosis_vs_nucleus`` is P(a random mitotic figure outscores a random ordinary
    nucleus); 0.5 means no signal. ``median_mitosis_rank`` = ``(1 - auc) * n_nuclei`` is
    roughly where the typical mitotic figure lands in a ranking of every nucleus in the ROI
    -- compare it against ``K`` to see how far short the ranking falls.

    A local max over a +-``pad`` window is taken rather than the exact centre pixel, since
    the annotated click need not sit on the correlation peak.
    """
    cfg = cfg or fs.FSConfig()
    rng = np.random.default_rng(1) if rng is None else rng
    path = f"{images_dir}/{file_name}"

    rgb = ds.load_roi(path)
    img = ch.to_channel(rgb, cfg.channel)
    gt = ds.image_annotations(annotations, file_name)
    seed = gt[gt["ann_id"] == seed_ann_id].iloc[0]
    gt_eval = gt[gt["ann_id"] != seed_ann_id]
    mit = gt_eval[gt_eval["category_id"] == ds.MITOTIC]
    look = gt_eval[gt_eval["category_id"] == ds.LOOKALIKE]

    patch = tm.read_padded_patch(img, seed["cx"], seed["cy"], cfg.patch_size)
    templates, _ = tm.build_augmentations(patch, cfg.base_size, cfg.scales, cfg.n_angles, cfg.flips)
    fused, _, valid = tm.fused_response(img, templates, cfg.scale_normalize)
    fused = np.where(valid, fused, -1e9)

    mask = bl.tissue_mask(rgb)
    blobs = bl.nucleus_blobs(rgb, mask)
    rand = bl.random_in_tissue(mask, n_random, rng)

    def peaks(xs, ys):
        return np.array([
            fused[max(0, int(round(y)) - pad): int(round(y)) + pad + 1,
                  max(0, int(round(x)) - pad): int(round(x)) + pad + 1].max()
            for x, y in zip(xs, ys)
        ])

    s_mit, s_look = peaks(mit["cx"], mit["cy"]), peaks(look["cx"], look["cy"])
    s_nuc, s_rnd = peaks(blobs["cx"], blobs["cy"]), peaks(rand["cx"], rand["cy"])
    s_seed = float(peaks([seed["cx"]], [seed["cy"]])[0])
    auc = lambda s: float(np.mean([(s_nuc < v).mean() for v in s])) if len(s) else float("nan")
    auc_mit = auc(s_mit)

    del rgb, img, fused, valid  # ~150 MB each; everything needed is already sampled
    return {
        "file_name": file_name,
        "n_mitotic": len(mit), "n_nuclei": len(blobs),
        "mitotic_base_rate_pct": round(100 * len(mit) / max(len(blobs), 1), 3),
        "seed_self_score": round(s_seed, 4),
        "median_score_mitotic": round(float(np.median(s_mit)), 4),
        "median_score_lookalike": round(float(np.median(s_look)), 4),
        "median_score_nucleus": round(float(np.median(s_nuc)), 4),
        "median_score_random": round(float(np.median(s_rnd)), 4),
        "p90_mitotic": round(float(np.percentile(s_mit, 90)), 4),
        "p99.9_nucleus": round(float(np.percentile(s_nuc, 99.9)), 4),
        "auc_mitosis_vs_nucleus": round(auc_mit, 3),
        "auc_lookalike_vs_nucleus": round(auc(s_look), 3),
        "median_mitosis_rank": int((1 - auc_mit) * len(blobs)),
        "K": len(mit),
    }
