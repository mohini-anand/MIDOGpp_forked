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

from dataclasses import replace

from . import baselines as bl
from . import channels as ch
from . import dataset as ds
from . import evaluate as ev
from . import find_and_suppress as fs
from . import seed_selection as ss
from . import template_match as tm


def select_domain_images(images: pd.DataFrame, annotations: pd.DataFrame, images_dir="images"):
    """One ROI per tumour domain: the downloaded image with the most mitotic figures.

    Maximising category-1 count per image matters at this sample size -- with one seed
    and one pass, an image with 4 mitotic figures gives an evaluation set of 3. It also
    rules out 001.tiff automatically, which has zero category-1 annotations and so
    cannot be seeded with a mitotic figure at all. Ties break to the lower image id.

    **This rule is optimistic for the method being measured, and the results should be
    read that way.** `recall@K` uses K = the ROI's own mitotic count, so as mitosis
    density rises the base rate climbs against a roughly constant nucleus population and
    the metric gets easier. Picking the densest ROI per domain therefore yields
    upper-ish bounds for each domain, not representative draws from it.
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
    # `.head(1)`, not `.first()`: GroupBy.first() takes the first *non-null* value in each
    # column independently and can therefore assemble a row that does not exist in the
    # data. It happens to be safe here only because nothing is NaN after the fillna above.
    return df.groupby("tumor_type", as_index=False).head(1).reset_index(drop=True)


def pick_seed(gt_mitotic: pd.DataFrame, rng, border: int, roi_shape) -> pd.Series:
    """A random mitotic annotation, far enough from the border to build a full template.

    Annotations within half a padded-patch of the edge cannot yield a rotation-safe
    template; they are excluded here rather than padded, since padding would feed
    fabricated pixels into the template.

    The test is applied to the *rounded* centre, which is the predicate
    `template_match.read_padded_patch` actually uses. Testing the unrounded float was
    marginally looser and could pass an annotation the reader then rejected -- a seed at
    ``cx = w - border - 0.4`` satisfied ``cx < w - border`` but rounds up to
    ``w - border``, for which ``ix + half >= w``.
    """
    h, w = roi_shape[:2]
    ix = np.rint(gt_mitotic["cx"].to_numpy()).astype(int)
    iy = np.rint(gt_mitotic["cy"].to_numpy()).astype(int)
    ok = gt_mitotic[
        (ix >= border) & (ix <= w - 1 - border) & (iy >= border) & (iy <= h - 1 - border)
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
    bbox_method: str = "binary",
    bbox_center_tolerance: int = 0,
    bbox_headroom_frac: float = None,
    tighten_bbox: bool = True,
) -> dict:
    """Full pipeline + baselines for one ROI. Returns everything needed to plot or re-score.

    ``bbox_method`` selects `seed_selection.tighten_box_otsu`'s foreground threshold --
    ``"binary"`` (default, unchanged pipeline behaviour), ``"multiotsu"`` (`Research
    Logs/design_choices.md`, section 7), or ``"headroom"`` (`Research
    Logs/2026-09-03-bbox-threshold-sweep.md`). ``bbox_center_tolerance`` (default 0,
    unchanged) widens the centre-pixel check to a small neighbourhood -- see
    `tighten_box_otsu`'s ``center_tolerance``, section 7's tolerance follow-up.
    ``bbox_headroom_frac`` (default ``None``, unchanged) is `"headroom"`'s required
    ``(0, 1]`` fraction of the headroom between binary Otsu's split and the crop's max;
    unused by the other two methods. All three are passed identically to `pick_seed` and
    `tightened_base_size` so the seed a given setting accepts is retightened under that
    same setting, not a different one.

    ``tighten_bbox`` (default True, unchanged) toggles Otsu/CC bbox tightening off
    entirely -- when False, seed selection is pathologist agreement + the border filter
    only (`bbox_method`/``bbox_center_tolerance``/``bbox_headroom_frac`` are then unused,
    since there's no Otsu/CC step left to configure), and the template keeps
    ``cfg.base_size``'s native size (`FSConfig`'s own default, `template_match.BASE_SIZE`
    = 51px) instead of being resized to a tightened box. See `Research
    Logs/design_choices.md`, section 8.
    """
    cfg = cfg or fs.FSConfig()
    rng = np.random.default_rng(0) if rng is None else rng
    path = f"{images_dir}/{file_name}"

    rgb = ds.load_roi(path)
    mpp = ds.roi_mpp(path)
    area_mm2 = ds.check_roi_scale(path, rgb.shape)
    radius = ev.radius_px(mpp)

    gt = ds.image_annotations(annotations, file_name)
    # The structural channel for seed selection (pathologist agreement + Otsu/CC bbox
    # tightening) is always gray_inverted, independent of `cfg.channel` -- it is a test
    # of where the click sits, not a matching score. See seed_selection.py.
    gray_inv = ch.to_gray_inverted(rgb)
    seed, seed_info = ss.pick_seed(
        gt[gt["category_id"] == ds.MITOTIC], gray_inv, rng, cfg.patch_size // 2, rgb.shape,
        method=bbox_method, center_tolerance=bbox_center_tolerance,
        headroom_frac=bbox_headroom_frac, tighten_bbox=tighten_bbox,
    )
    gt_eval = gt[gt["ann_id"] != seed["ann_id"]].reset_index(drop=True)

    if tighten_bbox:
        # The template's native size after bbox tightening, replacing cfg.base_size's
        # fixed 51 px default. `pick_seed`'s foreground filter already guarantees this
        # seed's component passes the size/shape sanity check, so this cannot return
        # None here.
        tightened_size = ss.tightened_base_size(gray_inv, seed["cx"], seed["cy"], method=bbox_method,
                                                center_tolerance=bbox_center_tolerance,
                                                headroom_frac=bbox_headroom_frac)
        run_cfg = replace(cfg, base_size=tightened_size)
    else:
        # No tightening requested: keep cfg.base_size as-is (51px by default) and skip
        # the Otsu/CC step entirely -- `tightened_base_size` is never called.
        tightened_size = cfg.base_size
        run_cfg = cfg

    img = ch.to_channel(rgb, run_cfg.channel)
    t0 = time.time()
    # Suppress at the same radius the scoring uses, so one true object cannot be counted
    # as a hit and a duplicate false positive at the same time.
    det, info = fs.find_and_suppress(img, (seed["cx"], seed["cy"]), run_cfg, nms_radius=radius)
    info["t_total_s"] = round(time.time() - t0, 1)

    det_out, gt_out, metrics, curve = ev.evaluate_run(
        det, gt_eval, radius, area_mm2, roi_shape=rgb.shape
    )
    metrics = {"file_name": file_name, "method": "find_and_suppress",
               "seed_ann_id": int(seed["ann_id"]), "channel": run_cfg.channel,
               "mpp": round(mpp, 4), "tightened_base_size": tightened_size,
               "bbox_method": bbox_method, "bbox_center_tolerance": bbox_center_tolerance,
               "bbox_headroom_frac": bbox_headroom_frac,
               "tighten_bbox": tighten_bbox,
               "seed_agreement_flagged": seed_info.agreement_flagged,
               "seed_n_agreement_pool": seed_info.n_agreement_pool,
               "seed_n_after_border": seed_info.n_after_border,
               "seed_n_after_foreground": seed_info.n_after_foreground,
               **metrics, **info}

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
        tissue_fraction = round(float(mask.mean()), 4)
        out["tissue_fraction"] = tissue_fraction
        metrics["tissue_fraction"] = tissue_fraction

        blobs = bl.nucleus_blobs(rgb, mask)
        b_det, b_gt, b_metrics, b_curve = ev.evaluate_run(
            blobs, gt_eval, radius, area_mm2, roi_shape=rgb.shape
        )

        # `recall@K` already compares at an equal budget, but the full-list recall does
        # not -- the blob detector emits a different number of candidates than the
        # matcher. Score it again truncated to the matcher's list length so the two are
        # comparable at a second, much looser operating point as well. Record what was
        # actually delivered: if the blob detector has fewer components than the matcher
        # has detections, this second comparison is not equal-budget either.
        capped = blobs.head(len(det_out)).reset_index(drop=True)
        capped["rank"] = np.arange(len(capped))
        _, _, capped_metrics, _ = ev.evaluate_run(
            capped, gt_eval, radius, area_mm2, roi_shape=rgb.shape
        )
        b_metrics["mitotic_recall_at_matcher_budget"] = capped_metrics["mitotic_recall"]
        b_metrics["matcher_budget"] = len(det_out)
        b_metrics["budget_delivered"] = len(capped)

        out["baseline_blobs"] = b_det
        out["baseline_blobs_curve"] = b_curve
        out["metrics"].append(
            {"file_name": file_name, "method": "nucleus_blobs", "seed_ann_id": -1,
             "channel": "hematoxylin", "mpp": round(mpp, 4),
             "tissue_fraction": tissue_fraction, **b_metrics}
        )

        rnd = bl.random_in_tissue(mask, n=len(det_out), rng=rng)
        r_det, r_gt, r_metrics, r_curve = ev.evaluate_run(
            rnd, gt_eval, radius, area_mm2, roi_shape=rgb.shape
        )
        out["baseline_random"] = r_det
        out["baseline_random_curve"] = r_curve
        out["metrics"].append(
            {"file_name": file_name, "method": "random_in_tissue", "seed_ann_id": -1,
             "channel": "-", "mpp": round(mpp, 4),
             "tissue_fraction": tissue_fraction, **r_metrics}
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
    bbox_method: str = "binary",
    bbox_center_tolerance: int = 0,
    bbox_headroom_frac: float = None,
    tighten_bbox: bool = True,
):
    """Loop `run_one_image` over the selected ROIs.

    ``keep_images=False`` drops the decoded RGB arrays from the returned results, since
    seven 39-megapixel ROIs will not comfortably sit in a notebook's memory at once.

    Each image gets its **own** generator, spawned from ``(seed, image_id)``. Building
    ``default_rng(seed)`` fresh inside the loop -- as this used to -- makes every image
    consume the same first draw, so `pick_seed` lands on the same relative position in
    every eligible list (it picked the ~80th percentile of all seven). That is one fixed
    quantile repeated, not seven independent draws, and a multi-seed sweep built on it
    would repeat the same quantile sequence for every seed.

    ``bbox_method``, ``bbox_center_tolerance``, ``bbox_headroom_frac``, and
    ``tighten_bbox`` are passed straight through to every `run_one_image` call -- see
    that function's docstring.
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
            rng=np.random.default_rng([seed, int(row["image_id"])]),
            run_baselines=run_baselines,
            bbox_method=bbox_method, bbox_center_tolerance=bbox_center_tolerance,
            bbox_headroom_frac=bbox_headroom_frac,
            tighten_bbox=tighten_bbox,
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
                f"coverage={m['coverage_frac']:.2f}  "
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
    "n_unanimous", "recall_unanimous_at_k", "n_contested", "recall_contested_at_k",
    "max_detection_score", "n_detections_total",
]

FULL_LIST_COLUMNS = [
    "tumor_type", "file_name", "all_n_detections", "coverage_frac", "all_tp",
    "all_fp_lookalike", "all_fp_unannotated", "all_matched_any_gt", "all_matched_frac",
    "all_precision_mitotic", "all_mitotic_gt_found", "all_mitotic_gt_missed",
    "all_lookalike_gt_found", "all_lookalike_gt_missed",
]


def full_list_table(metrics: pd.DataFrame, method="find_and_suppress") -> pd.DataFrame:
    """Every detection the search returned, bucketed -- no top-K budget applied.

    ``coverage_frac`` is placed second on purpose: it is the number that says how much of
    the rest of the row is explained by the detection list tiling the ROI.
    """
    cols = [c for c in FULL_LIST_COLUMNS if c in metrics.columns]
    sel = metrics if method is None else metrics[metrics["method"] == method]
    return sel[cols].sort_values("tumor_type").reset_index(drop=True)


def summary_table(metrics: pd.DataFrame, method="find_and_suppress") -> pd.DataFrame:
    cols = [c for c in SUMMARY_COLUMNS if c in metrics.columns]
    sel = metrics if method is None else metrics[metrics["method"] == method]
    return sel[cols].sort_values("tumor_type").reset_index(drop=True)


def scale_usage(metrics: pd.DataFrame, detections: pd.DataFrame, scales=(0.6, 0.8, 1.0)) -> pd.DataFrame:
    """Which template scale actually won each detection, versus uniform expectation.

    Only meaningful with ``len(FSConfig.scales) > 1``; at the single-scale default every
    column is 0 or 1 by construction.

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


def _rank_stats(s_target: np.ndarray, s_ref: np.ndarray):
    """(AUC, mean rank, median rank) of ``s_target`` against the reference population.

    AUC is ``P(ref < t) + 0.5 * P(ref == t)`` averaged over ``t`` -- the Mann-Whitney
    convention, which splits ties rather than scoring them as losses.

    "Rank" is how many reference scores beat a target score, i.e. where that target lands
    in a descending ranking of the reference population (0 = top). The mean rank is
    exactly ``(1 - AUC) * n_ref``; the median is reported separately because it is what
    the phrase "where the typical mitotic figure lands" actually means and is far less
    sensitive to the tail.
    """
    if len(s_target) == 0 or len(s_ref) == 0:
        return float("nan"), float("nan"), float("nan")
    srt = np.sort(np.asarray(s_ref, dtype=np.float64))
    n = len(srt)
    lo = np.searchsorted(srt, s_target, side="left")
    hi = np.searchsorted(srt, s_target, side="right")
    ranks = (n - hi) + 0.5 * (hi - lo)   # strictly-greater, ties split
    auc = float(1.0 - ranks.mean() / n)
    return auc, float(ranks.mean()), float(np.median(ranks))


def _probe_variants(file_name, annotations, seed_ann_id, variants, images_dir="images",
                    n_random=3000, pad=5, rng=None):
    """Score one or more matcher configurations against the same sampled populations.

    The expensive shared work -- decoding the ROI, building the tissue mask, segmenting
    every nucleus, drawing the random points -- is done **once** and reused across
    variants, so a fusion comparison costs one nucleus segmentation rather than one per
    variant, and every variant is scored against a byte-identical reference population.

    ``variants`` is a list of ``(name, FSConfig)``. Returns one row per variant.
    """
    from sklearn.neighbors import KDTree

    rng = np.random.default_rng(1) if rng is None else rng
    path = f"{images_dir}/{file_name}"

    rgb = ds.load_roi(path)
    mpp = ds.roi_mpp(path)
    radius = ev.radius_px(mpp)
    gt = ds.image_annotations(annotations, file_name)
    seed = gt[gt["ann_id"] == seed_ann_id].iloc[0]
    gt_eval = gt[gt["ann_id"] != seed_ann_id]
    mit = gt_eval[gt_eval["category_id"] == ds.MITOTIC]
    look = gt_eval[gt_eval["category_id"] == ds.LOOKALIKE]

    mask = bl.tissue_mask(rgb)
    blobs = bl.nucleus_blobs(rgb, mask)
    n_blobs_all = len(blobs)
    if n_blobs_all:
        # Drop every blob sitting on an annotation, so "nucleus" means "unannotated".
        near = KDTree(ds.points(gt)).query_radius(blobs[["cx", "cy"]].to_numpy(), r=radius)
        blobs = blobs[[len(c) == 0 for c in near]].reset_index(drop=True)
    rand = bl.random_in_tissue(mask, n_random, rng)
    del mask

    rows = []
    for name, cfg in variants:
        img = ch.to_channel(rgb, cfg.channel)
        patch = tm.read_padded_patch(img, seed["cx"], seed["cy"], cfg.patch_size)
        if patch is None:
            raise ValueError(f"seed {seed_ann_id} is too close to the border of {file_name}")
        templates, _ = tm.build_augmentations(
            patch, cfg.base_size, cfg.scales, cfg.n_angles, cfg.flips
        )
        fused, _, valid = tm.fused_response(img, templates, cfg.scale_normalize)
        fused = np.where(valid, fused, -1e9)

        def peaks(xs, ys):
            return np.array([
                fused[max(0, int(round(y)) - pad): int(round(y)) + pad + 1,
                      max(0, int(round(x)) - pad): int(round(x)) + pad + 1].max()
                for x, y in zip(xs, ys)
            ])

        s_mit, s_look = peaks(mit["cx"], mit["cy"]), peaks(look["cx"], look["cy"])
        s_nuc, s_rnd = peaks(blobs["cx"], blobs["cy"]), peaks(rand["cx"], rand["cy"])
        s_seed = float(peaks([seed["cx"]], [seed["cy"]])[0])
        del img, fused, valid  # ~150 MB each; everything needed is already sampled

        auc_mit, mean_rank, median_rank = _rank_stats(s_mit, s_nuc)
        auc_look, _, _ = _rank_stats(s_look, s_nuc)
        rows.append({
            "file_name": file_name,
            "variant": name,
            "n_templates": len(templates),
            "n_mitotic": len(mit),
            "n_nuclei_all_blobs": n_blobs_all,
            "n_nuclei": len(blobs),           # ordinary = unannotated
            "mitotic_base_rate_pct": round(100 * len(mit) / max(n_blobs_all, 1), 3),
            "seed_self_score": round(s_seed, 4),
            "median_score_mitotic": round(float(np.median(s_mit)), 4),
            "median_score_lookalike": round(float(np.median(s_look)), 4),
            "median_score_nucleus": round(float(np.median(s_nuc)), 4),
            "median_score_random": round(float(np.median(s_rnd)), 4),
            "p90_mitotic": round(float(np.percentile(s_mit, 90)), 4),
            "p99.9_nucleus": round(float(np.percentile(s_nuc, 99.9)), 4),
            "auc_mitosis_vs_nucleus": round(auc_mit, 3),
            "auc_lookalike_vs_nucleus": round(auc_look, 3),
            "discrimination": round(auc_mit - auc_look, 3),
            "mean_mitosis_rank": int(round(mean_rank)),
            "median_mitosis_rank": int(round(median_rank)),
            "K": len(mit),
        })
    del rgb
    return rows


def score_probe(file_name, annotations, seed_ann_id, images_dir="images", cfg=None,
                n_random=3000, pad=5, rng=None):
    """Does the response map score mitotic figures above ordinary nuclei at all?

    Detection metrics conflate two questions -- whether the map carries signal, and whether
    the ranking can surface it against the base rate. This separates them by reading the
    fused map directly at annotated locations, at detected nucleus centroids, and at random
    tissue points.

    ``auc_mitosis_vs_nucleus`` is P(a random mitotic figure outscores a random *ordinary*
    nucleus); 0.5 means no signal. "Ordinary" is enforced: blob centroids within the match
    radius of any annotation -- mitotic, look-alike, or the seed itself -- are removed from
    the reference population, so the comparison is against unannotated nuclei rather than
    against a set that silently contains the very objects being scored. Ties are split
    (Mann-Whitney) rather than scored as losses.

    ``mean_mitosis_rank`` = ``(1 - AUC) * n_nuclei`` and ``median_mitosis_rank`` say where
    a mitotic figure lands in a ranking of every ordinary nucleus in the ROI -- compare
    them against ``K`` to see how far short the ranking falls. The two used to be
    conflated: ``(1 - AUC) * n`` is an expectation, and was reported under the name
    "median".

    A local max over a +-``pad`` window is taken rather than the exact centre pixel, since
    the annotated click need not sit on the correlation peak.
    """
    cfg = cfg or fs.FSConfig()
    row = _probe_variants(file_name, annotations, seed_ann_id, [("default", cfg)],
                          images_dir, n_random, pad, rng)[0]
    row.pop("variant")
    row.pop("n_templates")
    return row


# The fusion comparison behind FSConfig.scales defaulting to a single 51 px template.
FUSION_VARIANTS = [
    ("fused_max_multiscale", dict(scales=(0.6, 0.8, 1.0), scale_normalize=False)),
    ("fused_z_normalised", dict(scales=(0.6, 0.8, 1.0), scale_normalize=True)),
    ("scale_1.0_only", dict(scales=(1.0,), scale_normalize=False)),
    ("scale_0.6_only", dict(scales=(0.6,), scale_normalize=False)),
]


def fusion_variants(file_names, annotations, seed_ann_ids, images_dir="images", cfg=None,
                    variants=None, n_random=3000, pad=5, rng=None) -> pd.DataFrame:
    """Compare response-map fusion strategies on the same images and reference nuclei.

    This is the evidence for `FSConfig.scales` defaulting to a single 51 px template, and
    it belongs in the notebook rather than in a comment: the conclusion is not the same on
    every image, and a table showing only one of them reads as more settled than it is.

    ``discrimination`` (mitosis-vs-nucleus AUC minus look-alike-vs-nucleus AUC) is
    reported alongside, because it -- not the raw AUC the default was chosen on -- is the
    quantity the experiment's central claim turns on.
    """
    base = cfg or fs.FSConfig()
    variants = FUSION_VARIANTS if variants is None else variants
    named = [(name, replace(base, **kw)) for name, kw in variants]
    rows = []
    for fn, ann_id in zip(file_names, seed_ann_ids):
        rows.extend(_probe_variants(fn, annotations, int(ann_id), named,
                                    images_dir, n_random, pad, rng))
    return pd.DataFrame(rows)


# The augmentation-footprint comparison from `Research Logs/design_choices.md`, section
# 4 -- deliberately not the 72-template (12 angles x 2 flips x 3 scales) grid `plant_and
# _recover`'s coordinate gate uses; that grid was already measured and rejected on the
# scale axis alone (see FUSION_VARIANTS above). `n_angles=4` gives exactly the angles
# {0, 90, 180, 270} (`build_augmentations` computes `angle = 360 * k / n_angles`), so no
# change to `template_match.py` was needed for either variant below.
#
# `no_augmentation` is also `FSConfig`'s own default as of section 6 -- listed here
# explicitly anyway (rather than left implicit in `base`) so this comparison keeps
# working regardless of what the pipeline default happens to be later.
AUGMENTATION_VARIANTS = [
    ("no_augmentation", dict(n_angles=1, flips=(False,))),
    ("angles12_flips2", dict(n_angles=12, flips=(False, True))),
    ("rot90_4angles_2flips", dict(n_angles=4, flips=(False, True))),
]


def augmentation_variants(file_names, annotations, seed_ann_ids, images_dir="images", cfg=None,
                          variants=None, n_random=3000, pad=5, rng=None) -> pd.DataFrame:
    """Compare augmentation footprints on the same images and reference nuclei.

    Same machinery as `fusion_variants` (same seed and sampled reference population per
    image across every variant, so only augmentation count differs) applied to the
    angle/flip axis instead of scale. Tests the doc's hypothesis: fusing more augmented
    maps via element-wise max inflates spurious high scores under noise -- the same
    mechanism the multi-scale fusion comparison rejected -- so `median_score_nucleus` /
    `p99.9_nucleus` should fall as augmentation count drops, at some cost to
    `auc_mitosis_vs_nucleus` on off-angle mitoses.
    """
    base = cfg or fs.FSConfig()
    variants = AUGMENTATION_VARIANTS if variants is None else variants
    named = [(name, replace(base, **kw)) for name, kw in variants]
    rows = []
    for fn, ann_id in zip(file_names, seed_ann_ids):
        rows.extend(_probe_variants(fn, annotations, int(ann_id), named,
                                    images_dir, n_random, pad, rng))
    return pd.DataFrame(rows)
