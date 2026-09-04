"""Which `cv2.TM_*` variant retrieves mitotic figures from one click, on hematoxylin OD.

`Research Logs/2026-09-01-tm-variant-sweep.md`. Before this script the entire project ran one
similarity: `grep -roh "TM_[A-Z_]*"` over every .py/.md/.ipynb returned 61 hits, all
`TM_CCOEFF_NORMED`. That is the one method invariant to `I -> aI + b`, i.e. blind to mean
intensity -- which `results/morph_diag_bhattacharyya.csv` measures as the strongest
mitotic-vs-ordinary-nucleus feature in all seven domains (0.925-3.703). The six OpenCV
variants are a 2x3 factorial over exactly the two invariances that matter, so "normed vs
unnormed" *is* "can the matcher see chromatin magnitude".

Three experiments, in order, each readable on its own:

  --exp1  TM_CCOEFF vs TM_CCOEFF_NORMED, paired on the identical click. The only
          single-variable contrast in the grid: same function, gain normalisation the only
          difference. Carries `pipeline_as_shipped` as a reference row so the comparison is
          answered both ways -- against a matched-everything twin, and against the
          configuration od_experiment.py/premise_test.py actually run today.
  --exp2  all six methods, plus the two controls that decide whether the click is doing the
          work: a click-free disc template, and a click on a pathologist-rejected look-alike.
  --exp3  rotation augmentation on the two best methods from exp2.

What makes the methods comparable at all is that no raw score threshold is ever used. Each
response map is reduced to its own median and MAD (`template_match.robust_stats`) and the
extraction floor is `med + FLOOR_Z*mad`. A robust z is monotone, so it changes no ranking; it
only makes one floor mean the same search depth for methods whose ranges differ by eight
orders of magnitude (TM_CCORR on OD spans [0.34, 10.7], negated TM_SQDIFF [-14.7, -1.9e-6]).

Two economies, both asserted rather than assumed:

* **One match per (ROI, seed, method).** `extract_peaks`'s local maxima are
  threshold-independent and greedy score-ordered NMS means a peak above t can only be
  suppressed by a higher-scoring peak also above t, so every z level is an exact *filter* of
  one deep pool. `verify_shortcut` checks it per method by coordinate-set equality.
* **One `nucleus_blobs` and one disc template per ROI**, reused across seeds -- neither
  depends on the click.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import cv2
import numpy as np
import pandas as pd

from midog_utils import baselines as bl
from midog_utils import channels as ch
from midog_utils import chromatin as cm
from midog_utils import compare as cp
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from midog_utils import find_and_suppress as fs
from midog_utils import invariants as inv
from midog_utils import seed_selection as ss
from midog_utils import template_match as tm
from midog_utils.nms import nms_by_distance

# --- configuration ----------------------------------------------------------------------
CHANNEL = "hematoxylin_od"   # unclipped OD; `to_hematoxylin` saturates 0.5% of pixels at 255
BASE_SIZE = tm.BASE_SIZE     # 51 px, fixed -- no Otsu tightening, see the module note below
FLOOR_Z = 0.5                # extraction floor in units of that map's own MAD
N_SEEDS = 5
DECISION_MIN_MITOTIC = 15    # `premise_test.DECISION_GRADE`'s bar, applied by count not name
SELF_HIT_RADIUS = 5.0
EXP1_METHODS = ("ccoeff", "ccoeff_normed")
ALL_METHODS = tuple(tm.METHODS)
SMOKE_ROI = "246.tiff"
RESULTS = Path("results")

# Template size is held at 51 px rather than Otsu-tightened. Tightening is not a *method*-axis
# confound (all arms share one template per run), but it varies 25-51 px with which cell is
# clicked and is a documented driver of huge list-length variance -- 301.tiff swings 1985 to
# 15484 detections from the seed alone. With 5 seeds, dropping the nuisance variable is worth
# more than matching the shipped default, and exp1 keeps the shipped configuration as its own
# reference arm anyway.


def roi_table(images_dir="images") -> pd.DataFrame:
    """Every downloaded ROI that can be seeded, with its domain and decision grade.

    001.tiff is excluded automatically: it has zero category-1 annotations, so it can neither
    be seeded with a mitotic figure nor scored for mitotic recall. This deliberately does NOT
    use `experiment.select_domain_images`, which keeps only the densest ROI per domain -- two
    ROIs per domain is what makes leave-one-domain-out possible, and leave-one-ROI-out leaks
    (`fp_filter_domain.py`).
    """
    images, ann = ds.load_annotations()
    on_disk = {p.name for p in Path(images_dir).glob("*.tiff")}
    counts = (ann[ann["category_id"] == ds.MITOTIC].groupby("file_name").size()
              .rename("n_mitotic"))
    df = images[images["file_name"].isin(on_disk)].merge(
        counts, left_on="file_name", right_index=True, how="left")
    df["n_mitotic"] = df["n_mitotic"].fillna(0).astype(int)
    df = df[df["n_mitotic"] > 0].copy()
    df["decision_grade"] = df["n_mitotic"] >= DECISION_MIN_MITOTIC
    return df.sort_values(["tumor_type", "file_name"]).reset_index(drop=True)


def geometric_ceiling(gt_mitotic: pd.DataFrame, radius: float) -> dict:
    """Annotation crowding, reported as a diagnostic. **`geometric_ceiling` is NOT a bound.**

    It was originally written as one, on the argument that two mitotic figures closer together
    than the NMS radius cannot both survive suppression. That argument is false, and an
    independent audit caught it: the detections do not sit on the ground-truth centres, so two
    surviving peaks can be more than a radius apart from *each other* while each still lies
    within a match radius of a *different* annotation. Empirically, 245.tiff -- the only ROI
    with such a pair (26.6 px against a 30.2 px radius) -- reaches `full_list_recall` 1.0000 in
    22 cells, with both members of the pair matched and neither ever the seed.

    It is also on a different denominator from everything it would be compared against: this
    counts all mitotic annotations, while every metric downstream is measured against
    `gt_eval`, which excludes the seed.

    Both returned values are **reported columns only** -- never a filter, a cap, or an input to
    any metric -- so they corrupt nothing. `n_mitotic_crowded` is the honest half: it says how
    many annotations sit inside a suppression radius of another, which is worth knowing.
    """
    from sklearn.neighbors import KDTree
    p = gt_mitotic[["cx", "cy"]].to_numpy()
    if len(p) < 2:
        return {"n_mitotic_crowded": 0, "geometric_ceiling": 1.0}
    d = KDTree(p).query(p, k=2)[0][:, 1]
    crowded = int((d < radius).sum())
    # A crowded pair can still yield one of its two members, so the loss is at most half.
    return {"n_mitotic_crowded": crowded,
            "geometric_ceiling": float(1.0 - np.ceil(crowded / 2) / len(p))}


def disc_template(hem: np.ndarray, size: int = BASE_SIZE) -> np.ndarray:
    """A generic dark round blob carrying no click information, scaled to this ROI's OD.

    The null for "does the click matter". Inside/outside values are this ROI's 99.5th and 50th
    OD percentiles, so the disc has a nucleus-like magnitude but none of the clicked cell's
    texture. Magnitude-matching is not cosmetic: TM_SQDIFF is `sum((T-I)^2)` and so depends on
    the template's absolute scale, unlike the other five where a constant factor on T leaves
    the ranking unchanged.

    Built once per ROI -- it does not depend on the seed.
    """
    lo, hi = np.percentile(hem[::8, ::8], [50.0, 99.5])
    yy, xx = np.mgrid[-(size // 2):size // 2 + 1, -(size // 2):size // 2 + 1]
    inside = (xx * xx + yy * yy) <= (size // 2) ** 2
    return np.where(inside, hi, lo).astype(np.float32)


class DegenerateMapError(RuntimeError):
    """A response map from which no honest candidate pool can be built. See `match_pool`."""


def match_pool(img, hem, template, radius, seed_xy, method, cfg, keep_map=False):
    """One match, one NMS, at the deep floor -- the pool every z level is derived from.

    ``seed_xy`` may be None (the disc control has no click, so nothing to self-suppress).
    Returns ``(pool, diag, maps_or_None)``; ``pool`` carries ``cx, cy, score, od``.
    """
    fused, _, valid = tm.fused_response(img, template, cfg.scale_normalize,
                                        method=tm.METHODS[method])
    med, mad = tm.robust_stats(fused, valid)
    floor = med + FLOOR_Z * mad

    # A degenerate response map, caught loudly instead of silently producing a pool.
    # `TM_SQDIFF_NORMED` is clamped to [0,1] by OpenCV, and on hematoxylin OD 83% of the map
    # saturates there -- so median = -1.0 (post sign-flip), MAD = 0.0, and `med + z*mad` lands
    # exactly ON the clamp for every z. The floor then admits ~28M tied local maxima, the
    # `max_peaks` cap keeps the first 250k, and since `extract_peaks` breaks ties by x the
    # survivors are a left-edge strip. Every downstream number for that arm is an artefact of
    # the tie-break, not a measurement. Worse, `n_matched` is a min over arms, so one such arm
    # silently truncates every other arm in the same cell.
    if not np.isfinite(mad) or mad <= 0:
        raise DegenerateMapError(
            f"method {method!r}: response map has MAD={mad} (median={med}); the extraction "
            f"floor is degenerate and any pool built from it would be a tie-break artefact. "
            f"{float(np.isclose(fused[valid], med).mean()):.1%} of reachable pixels sit on the "
            "median. This arm cannot be measured on this image."
        )

    centers, scores = tm.extract_peaks(fused, valid, cfg.peak_min_distance, floor,
                                       cfg.max_peaks)
    n_peaks = len(centers)
    if n_peaks >= cfg.max_peaks:
        # `invariants.check_no_cap` is applied downstream to the post-NMS length, which is not
        # the quantity that binds -- the cap bites here, inside extraction. That gap is why 228
        # capped rows passed every gate in the first exp2 run.
        raise DegenerateMapError(
            f"method {method!r}: extract_peaks returned exactly max_peaks={cfg.max_peaks}, so "
            "the candidate list is a truncation of the score ordering rather than everything "
            "above the floor. Raise max_peaks or the floor; do not report this arm."
        )
    keep = nms_by_distance(centers, scores, radius)
    centers, scores = centers[keep], scores[keep]
    n_self = 0
    if seed_xy is not None and len(centers):
        ok = np.hypot(centers[:, 0] - seed_xy[0], centers[:, 1] - seed_xy[1]) > SELF_HIT_RADIUS
        n_self = int((~ok).sum())
        centers, scores = centers[ok], scores[ok]

    pool = pd.DataFrame({"cx": centers[:, 0], "cy": centers[:, 1], "score": scores})
    pool = cm.score_detections(pool, hem)
    diag = {"map_median": med, "mad_scale": mad, "extract_floor": floor,
            "n_peaks": n_peaks, "n_pool": len(pool), "n_self_hits": n_self}
    return pool, diag, ((fused, valid) if keep_map else None)


def verify_shortcut(fused, valid, radius, seed_xy, cfg, pool, cut) -> dict:
    """Assert that filtering the deep pool at ``cut`` equals extracting at ``cut``.

    Strict: the two coordinate sets must be identical after sorting, not merely the same size.
    The argument holds for any monotone score, negated TM_SQDIFF included -- which is exactly
    why it is re-asserted per method instead of taken on faith from `premise_test`.
    """
    centers, scores = tm.extract_peaks(fused, valid, cfg.peak_min_distance, cut, cfg.max_peaks)
    centers = centers[nms_by_distance(centers, scores, radius)]
    if seed_xy is not None and len(centers):
        centers = centers[np.hypot(centers[:, 0] - seed_xy[0],
                                   centers[:, 1] - seed_xy[1]) > SELF_HIT_RADIUS]
    filt = pool.loc[pool["score"] >= cut, ["cx", "cy"]].to_numpy()
    a = centers[np.lexsort((centers[:, 1], centers[:, 0]))] if len(centers) else centers
    b = filt[np.lexsort((filt[:, 1], filt[:, 0]))] if len(filt) else filt
    return {"check": "one_match_many_z", "n_direct": int(len(a)), "n_filtered": int(len(b)),
            "passed": bool(a.shape == b.shape and np.array_equal(a, b))}


def top(df: pd.DataFrame, key: str, n: int | None = None) -> pd.DataFrame:
    """Rank best-first by ``key`` (stable, NaN last), optionally truncated to ``n``."""
    out = df.sort_values(key, ascending=False, na_position="last",
                         kind="mergesort").reset_index(drop=True)
    if n is not None:
        out = out.head(n).reset_index(drop=True)
    return out.assign(rank=np.arange(len(out)))


def draw_seeds(gt_mitotic, border, roi_shape, n_seeds):
    """Up to ``n_seeds`` *distinct* mitotic annotations, one RNG stream per seed index.

    Agreement tiering then the border filter, matching `seed_selection.pick_seed` -- but with
    the Otsu foreground filter skipped, because this experiment does not tighten the template
    and running a gate whose only job is to size a box we are not going to size would silently
    change which cells are eligible.

    Drawn **without replacement**, unlike `pick_seed`. `check_distinct_seeds` documents that
    with-replacement collisions are expected and correct for that function, but here the
    headline is worst-of-N and a repeated annotation would make a worst-of-5 secretly a
    worst-of-4. Each seed index still gets its own `default_rng([s, image_id])` stream, so the
    invariant check still sees five independent streams; only the pool shrinks between draws.

    Returns ``(seeds, info)`` -- fewer than ``n_seeds`` when the pool is smaller, which is
    expected on the sparse ROIs (505.tiff has 2 mitotic figures in total).
    """
    pool, flagged = ss.agreement_pool(gt_mitotic)
    n_pool = len(pool)
    pool = ss.border_filter(pool, border, roi_shape)
    remaining, seeds, records = pool, [], []
    for s in range(n_seeds):
        if len(remaining) == 0:
            break
        rng = np.random.default_rng([s, int(gt_mitotic["image_id"].iloc[0])])
        i = int(rng.integers(len(remaining)))
        row = remaining.iloc[i]
        seeds.append((s, row))
        records.append((s, int(row["ann_id"]), [s, int(gt_mitotic["image_id"].iloc[0])]))
        remaining = remaining.drop(remaining.index[i])
    info = ss.SeedInfo(flagged, n_pool, len(pool), len(pool))
    return seeds, info, records


def _method_arms(pools, blobs, n_matched, tag, fn, s, cfg, caps):
    """Full and pool-matched arms for one seed's set of method pools, plus the blob bar.

    Every arm is emitted twice, ranked by its own similarity and by chromatin density, because
    both come free off the same pool and `results/od_workload_ab.csv` puts a median 19.4x
    between them under the incumbent method. `coverage_key` is shared between an arm's two
    rank keys -- coverage is a property of the candidate *set*, not of the sort order, and
    recomputing it per key would double the most expensive diagnostic in the run.

    Matched arms exist because `ceiling` is otherwise not evidence: a deep pool tiles the ROI
    and reaches 1.0 by geometry, which is the claim this repo already retracted (coverage_frac
    0.91-0.97; `random_in_tissue` hitting full-list recall 1.000). Truncating every arm to the
    smallest pool in the comparison is what makes "TM's ceiling beats the blob detector's"
    testable rather than vacuous.
    """
    arms = []
    for m, pool in pools.items():
        ck = f"{fn}:{s}:{tag}:{m}"
        for key, suffix in (("score", ""), ("od", "_od")):
            arms.append(cp.Arm(f"{tag}_{m}{suffix}", (lambda d=pool, k=key: top(d, k)),
                               rank_key=key, seeded=(tag != "disc"), caps=caps,
                               coverage_key=f"{ck}:full",
                               extra={"method": m, "pool_scope": "full",
                                      "generator": f"tm_{m}"}))
        arms.append(cp.Arm(f"{tag}_{m}@matched",
                           (lambda d=pool: top(d, "score", n_matched)),
                           rank_key="score", seeded=(tag != "disc"),
                           coverage_key=f"{ck}:matched",
                           extra={"method": m, "pool_scope": "matched",
                                  "generator": f"tm_{m}"}))
    if blobs is not None:
        for key, name in (("score", "blob_native"), ("od", "blob_od")):
            arms.append(cp.Arm(name, (lambda d=blobs, k=key: top(d, k)), rank_key=key,
                               coverage_key=f"{fn}:blob:full",
                               extra={"method": "-", "pool_scope": "full",
                                      "generator": "nucleus_blobs"}))
        # The seed index MUST be in this key. `n_matched` is the smallest pool in that
        # (ROI, seed)'s comparison and varies between seeds -- 6727-8278 within 201.tiff -- so
        # a key of f"{fn}:blob:matched" cached one seed's coverage and served it for all five.
        # The full-pool blob arms above are genuinely seed-independent and correctly share a
        # key. This is precisely the mis-keying `compare.evaluate_arms`' docstring warns about,
        # and it landed on the column used to adjudicate the ceiling comparison.
        arms.append(cp.Arm("blob_native@matched",
                           (lambda d=blobs: top(d, "score", n_matched)), rank_key="score",
                           coverage_key=f"{fn}:{s}:blob:matched",
                           extra={"method": "-", "pool_scope": "matched",
                                  "generator": "nucleus_blobs"}))
    return arms


def run_roi(row, methods, *, n_seeds=N_SEEDS, n_angles=1, flips=(False,),
            controls=False, shipped=False, tag="tm"):
    """Everything for one ROI. Returns ``(rows, checks, verification)``.

    The ROI is loaded, deconvolved and blob-segmented **once**; the per-seed loop then costs
    only the matches themselves. `nucleus_blobs` is ~30 s and does not depend on the click, so
    computing it per seed would triple the run for no information.
    """
    fn, image_id = row["file_name"], int(row["image_id"])
    path = f"images/{fn}"
    t_roi = time.time()

    rgb = ds.load_roi(path)
    mpp = ds.roi_mpp(path)
    radius = ev.radius_px(mpp)
    hem = np.ascontiguousarray(cm.hematoxylin_od(rgb))
    img = ch.to_channel(rgb, CHANNEL)
    gray = ch.to_gray_inverted(rgb)
    gt = ds.image_annotations(ds.load_annotations()[1], fn)
    gt_mit = gt[gt["category_id"] == ds.MITOTIC]
    gt_look = gt[gt["category_id"] == ds.LOOKALIKE]

    cfg = fs.FSConfig(channel=CHANNEL, base_size=BASE_SIZE, scales=(1.0,),
                      n_angles=n_angles, flips=flips)
    caps = (cfg.max_peaks, cfg.max_detections)
    border = cfg.patch_size // 2

    checks, verif, frames, unmeasurable = [], [], [], []
    mask = bl.tissue_mask(rgb)
    checks.append(inv.check_tissue_mask_covers_gt(mask, gt, label=fn))
    blobs = cm.score_detections(bl.nucleus_blobs(rgb, mask), hem)

    seeds, seed_info, records = draw_seeds(gt_mit, border, rgb.shape, n_seeds)
    geo = geometric_ceiling(gt_mit, radius)
    base_ctx = {"file_name": fn, "image_id": image_id, "tumor_type": row["tumor_type"],
                "decision_grade": bool(row["decision_grade"]), "mpp": round(mpp, 4),
                "channel": CHANNEL, "base_size": BASE_SIZE, "n_angles": n_angles,
                "n_flips": len(flips), "n_seeds_delivered": len(seeds),
                "seed_agreement_flagged": seed_info.agreement_flagged,
                "n_gt_mitotic_total": int(len(gt_mit)), **geo}
    print(f"[{fn}] {row['tumor_type']} n_mit={len(gt_mit)} blobs={len(blobs)} "
          f"seeds={len(seeds)}/{n_seeds} geo_ceiling={geo['geometric_ceiling']:.3f} "
          f"[{time.time() - t_roi:.0f}s]", flush=True)

    # --- the click-free disc control: one template, one match per method, per ROI ----------
    disc_pools = {}
    if controls:
        disc = [disc_template(hem)]
        for m in methods:
            try:
                disc_pools[m], d, _ = match_pool(img, hem, disc, radius, None, m, cfg)
            except DegenerateMapError as exc:
                unmeasurable.append({"file_name": fn, "seed_index": -1, "method": f"disc_{m}",
                                     "reason": str(exc)})
                print(f"  [{fn}] disc/{m:14s} UNMEASURABLE", flush=True)
                continue
            print(f"  [{fn}] disc/{m:14s} pool={d['n_pool']}", flush=True)

    cov_cache: dict = {}
    for s, seed in seeds:
        t_seed = time.time()
        seed_xy = (float(seed["cx"]), float(seed["cy"]))
        gt_eval = gt[gt["ann_id"] != int(seed["ann_id"])].reset_index(drop=True)
        patch = tm.read_padded_patch(img, *seed_xy, cfg.patch_size)
        templates, _ = tm.build_augmentations(patch, BASE_SIZE, cfg.scales, n_angles, flips)

        pools, diags = {}, {}
        for m in methods:
            keep = (s == 0)
            try:
                pools[m], diags[m], maps = match_pool(img, hem, templates, radius, seed_xy, m,
                                                      cfg, keep_map=keep)
            except DegenerateMapError as exc:
                # Record and carry on. A method whose response map cannot yield an honest pool
                # on this image should drop out of the grid, not abort the other five -- and it
                # must be *visible* as unmeasurable rather than absent, or a reader will take
                # the gap for a missing run. It is also excluded from `n_matched` below, so a
                # degenerate arm can no longer set the truncation depth for everyone else.
                unmeasurable.append({"file_name": fn, "seed_index": s, "method": m,
                                     "reason": str(exc)})
                print(f"  [{fn}] seed {s} {m}: UNMEASURABLE -- {str(exc)[:70]}", flush=True)
                continue
            if keep:
                cut = diags[m]["map_median"] + 2.0 * diags[m]["mad_scale"]
                v = verify_shortcut(maps[0], maps[1], radius, seed_xy, cfg, pools[m], cut)
                v.update({"file_name": fn, "seed_index": s, "method": m, "z": 2.0})
                verif.append(v)
                del maps

        if not pools:
            print(f"  [{fn}] seed {s}: every method unmeasurable, skipping", flush=True)
            continue
        n_matched = min([len(p) for p in pools.values()] + [len(blobs)])
        arms = _method_arms(pools, blobs, n_matched, tag, fn, s, cfg, caps)
        if controls and disc_pools:
            arms += _method_arms(disc_pools, None, n_matched, "disc", fn, s, cfg, ())

        if shipped:
            arms.append(_shipped_arm(rgb, gray, hem, seed_xy, radius, cfg, fn, s))

        ctx = dict(base_ctx, seed_index=s, seed_ann_id=int(seed["ann_id"]),
                   seed_cx=round(seed_xy[0], 1), seed_cy=round(seed_xy[1], 1),
                   n_matched_pool=n_matched,
                   map_median=round(diags[methods[0]]["map_median"], 6),
                   mad_scale=round(diags[methods[0]]["mad_scale"], 6))
        got = cp.evaluate_arms(arms, gt_eval, radius, roi_shape=rgb.shape, mpp=mpp,
                               context=ctx, coverage_cache=cov_cache, checks=checks)
        for m in methods:
            sel = got["method"] == m
            got.loc[sel, "extract_floor"] = round(diags[m]["extract_floor"], 6)
            got.loc[sel, "n_peaks"] = diags[m]["n_peaks"]
        frames.append(got)
        print(f"  [{fn}] seed {s} ann={int(seed['ann_id'])} matched={n_matched} "
              f"pools={ {m: diags[m]['n_pool'] for m in methods} } "
              f"[{time.time() - t_seed:.0f}s]", flush=True)

    if records:
        checks.append(inv.check_distinct_seeds(records, pool_size=seed_info.n_after_border,
                                               label=fn))

    # --- the look-alike-click control: one seed per ROI, all methods ----------------------
    if controls and len(gt_look):
        lk = ss.border_filter(gt_look, border, rgb.shape)
        if len(lk):
            rng = np.random.default_rng([image_id, 2])
            sd = lk.iloc[int(rng.integers(len(lk)))]
            # Drop BOTH the mitotic seed and the look-alike seed, so the two arms are scored
            # against identical ground truth and neither can rediscover its own template.
            mit0 = int(seeds[0][1]["ann_id"]) if seeds else -1
            gt_prov = gt[~gt["ann_id"].isin([mit0, int(sd["ann_id"])])].reset_index(drop=True)
            p = tm.read_padded_patch(img, float(sd["cx"]), float(sd["cy"]), cfg.patch_size)
            if p is not None:
                tpl, _ = tm.build_augmentations(p, BASE_SIZE, cfg.scales, n_angles, flips)
                lk_pools = {}
                for m in methods:
                    lk_pools[m], _, _ = match_pool(img, hem, tpl, radius,
                                                   (float(sd["cx"]), float(sd["cy"])), m, cfg)
                n_m = min([len(x) for x in lk_pools.values()] + [len(blobs)])
                ctx = dict(base_ctx, seed_index=0, seed_ann_id=int(sd["ann_id"]),
                           seed_cx=round(float(sd["cx"]), 1),
                           seed_cy=round(float(sd["cy"]), 1), n_matched_pool=n_m)
                frames.append(cp.evaluate_arms(
                    _method_arms(lk_pools, None, n_m, "lookalike", fn, "lk", cfg, caps),
                    gt_prov, radius, roi_shape=rgb.shape, mpp=mpp, context=ctx,
                    coverage_cache=cov_cache, checks=checks))

    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if unmeasurable:
        checks.extend({"check": "degenerate_map", "passed": False, **u} for u in unmeasurable)
    print(f"[{fn}] done, {len(out)} rows, {len(unmeasurable)} unmeasurable "
          f"[{time.time() - t_roi:.0f}s]", flush=True)
    return out, checks, verif


def _shipped_arm(rgb, gray, hem, seed_xy, radius, cfg, fn, s):
    """`TM_CCOEFF_NORMED` in the configuration this repo actually runs today.

    `od_experiment.py`, `od_seed_sweep.py`, `od_workload_ab.py` and `premise_test.py` all set
    ``CHANNEL = "rgb"`` and take the template's size from `seed_selection.tightened_base_size`.
    So exp1's paired arms answer "does removing the gain normalisation help, holding everything
    else fixed", and this row answers the different question "how does that compare with the
    thing on disk". It moves two variables at once and is therefore a reference, never part of
    the paired test.
    """
    img_rgb = ch.to_channel(rgb, "rgb")
    base = ss.tightened_base_size(gray, seed_xy[0], seed_xy[1]) or BASE_SIZE
    patch = tm.read_padded_patch(img_rgb, *seed_xy, cfg.patch_size)
    tpl, _ = tm.build_augmentations(patch, base, (1.0,), 1, (False,))
    pool, _, _ = match_pool(img_rgb, hem, tpl, radius, seed_xy, "ccoeff_normed", cfg)
    return cp.Arm("pipeline_as_shipped", (lambda d=pool: top(d, "score")), rank_key="score",
                  seeded=True, coverage_key=f"{fn}:{s}:shipped",
                  extra={"method": "ccoeff_normed", "pool_scope": "full",
                         "generator": "shipped_rgb_tightened", "shipped_base_size": int(base)})


def seed_variance(df: pd.DataFrame) -> pd.DataFrame:
    """Spread of each arm's tail across the seeds, on one row per (ROI, arm, metric).

    Reported because seed choice is the dominant noise term and medians have repeatedly hidden
    it: `read_50` on 301.tiff at the previously committed z=2.5 reads `247, -, 734, 313, -`
    across five seeds -- 2 of 5 never reach 50% -- and the detection list on that ROI swings
    1985-15484 purely from which annotation was clicked.

    ``n_unreached`` counts seeds where the quantile was never reached. Those are **not**
    imputed and **not** silently dropped: the spread is computed over reached seeds only and
    carried next to the count, because a mean over 3 of 5 seeds is not comparable to a mean
    over 5. That is precisely the defect the premise-test audit found, where a reported "28"
    was the median of {28, 28, 200}.

    ``cv`` is the cross-arm comparable number. A method 1.3x better at the median but 3x more
    seed-sensitive is not an improvement for a tool where the pathologist gets one click.
    """
    metrics = ["read_100", "read_95", "read_50", "n_detections", "full_list_recall"]
    at250 = (df[df["budget"] == 250].set_index(["file_name", "arm", "seed_index"])
             ["recall_at_budget"].rename("recall_at_250"))
    base = (df[df["budget"] == df["budget"].min()]
            .set_index(["file_name", "arm", "seed_index"]).sort_index())
    wide = base[metrics].join(at250).sort_index()
    rows = []
    for (fn, arm), g in wide.groupby(level=[0, 1]):
        meta = base.loc[(fn, arm)].iloc[0]
        for metric in metrics + ["recall_at_250"]:
            v = pd.to_numeric(g[metric], errors="coerce")
            ok = v.dropna()
            q1, q3 = (np.percentile(ok, [25, 75]) if len(ok) else (np.nan, np.nan))
            rows.append({
                "file_name": fn, "arm": arm, "metric": metric,
                "tumor_type": meta["tumor_type"], "decision_grade": meta["decision_grade"],
                "method": meta.get("method", "-"), "pool_scope": meta.get("pool_scope", "-"),
                "n_seeds": int(len(v)), "n_unreached": int(v.isna().sum()),
                "mean": float(ok.mean()) if len(ok) else np.nan,
                "sd": float(ok.std(ddof=1)) if len(ok) > 1 else np.nan,
                "cv": float(ok.std(ddof=1) / ok.mean()) if len(ok) > 1 and ok.mean() else np.nan,
                "min": float(ok.min()) if len(ok) else np.nan,
                "median": float(ok.median()) if len(ok) else np.nan,
                "max": float(ok.max()) if len(ok) else np.nan,
                "iqr": float(q3 - q1) if len(ok) else np.nan,
                # The headline. For a depth, worst = the deepest read; for a recall, the
                # lowest. Undefined when any seed never reached the quantile at all.
                "worst": (np.nan if v.isna().any() else
                          (float(ok.max()) if metric.startswith(("read", "n_")) else float(ok.min()))),
            })
    return pd.DataFrame(rows)


def paired_delta(df: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    """Per-(ROI, seed) delta between two arms that saw the identical click.

    This is the statistic exp1's conclusion rests on. Because both arms share the seed, the
    seed variance quantified by `seed_variance` cancels exactly, which no unpaired comparison
    on 13 ROIs could achieve.
    """
    metrics = ["read_100", "read_95", "read_50", "full_list_recall", "n_detections"]
    base = df[df["budget"] == df["budget"].min()]
    idx = ["file_name", "tumor_type", "decision_grade", "seed_index"]
    A = base[base["arm"] == a].set_index(idx)[metrics]
    B = base[base["arm"] == b].set_index(idx)[metrics]
    at250 = df[df["budget"] == 250].set_index(idx + ["arm"])["recall_at_budget"].unstack("arm")
    out = (A - B).add_suffix("_delta").join(
        (at250[a] - at250[b]).rename("recall_at_250_delta"))
    return out.assign(arm_a=a, arm_b=b).reset_index()


# --- gates ---------------------------------------------------------------------------------
def run_gates(gate_roi="002.tiff") -> pd.DataFrame:
    """Gates 1, 3, 4 and 7 -- everything checkable without running the grid.

    Gates 2 (`verify_shortcut`), 5 (the SQDIFF_NORMED clamp) and 6 (tissue mask, NMS radius,
    distinct seeds, no-cap) are per-ROI and are asserted inside `run_roi` instead, where the
    data they need already exists.
    """
    rows = []
    rgb = ds.load_roi(f"images/{gate_roi}")
    hem = np.ascontiguousarray(cm.hematoxylin_od(rgb))
    ann = ds.load_annotations()[1]
    seed = ds.image_annotations(ann, gate_roi).query("category_id == 1").iloc[0]
    patch = tm.read_padded_patch(hem, seed["cx"], seed["cy"], tm.PATCH_SIZE)

    for label, (na, fl) in {"single": (1, (False,)), "rot4x2": (4, (False, True))}.items():
        templates, metas = tm.build_augmentations(patch, BASE_SIZE, (1.0,), na, fl)
        for name, m in tm.METHODS.items():
            # Gate 1: the fused map's peak must land where the template was planted. Catches
            # the per-map half-size offset bug, and any sign error -- un-negated TM_SQDIFF
            # puts the planted patch at the map's minimum and this fails immediately.
            try:
                res = tm.plant_and_recover(templates, metas, method=m,
                                           rng=np.random.default_rng(0))
                off = max(max(abs(dx), abs(dy)) for _, dx, dy, _ in res)
                rows.append({"gate": "plant_and_recover", "bank": label, "method": name,
                             "detail": f"max|offset|={off}px over {len(res)} augmentations",
                             "passed": True})
            except AssertionError as e:
                rows.append({"gate": "plant_and_recover", "bank": label, "method": name,
                             "detail": str(e)[:160], "passed": False})

    T = np.ascontiguousarray(templates[0])
    rng = np.random.default_rng(0)
    H, W = hem.shape[0] - T.shape[0], hem.shape[1] - T.shape[1]
    rs, cs = rng.integers(0, H, 200), rng.integers(0, W, 200)
    for name, m in tm.METHODS.items():
        res = np.asarray(cv2.matchTemplate(hem, T, m), dtype=np.float32)
        errs = []
        for r, c in zip(rs, cs):
            Wd = hem[r:r + T.shape[0], c:c + T.shape[1]].astype(np.float64)
            Td = T.astype(np.float64)
            if name == "sqdiff":
                ref = ((Td - Wd) ** 2).sum()
            elif name == "sqdiff_normed":
                ref = min(((Td - Wd) ** 2).sum()
                          / np.sqrt((Td ** 2).sum() * (Wd ** 2).sum()), 1.0)
            elif name == "ccorr":
                ref = (Td * Wd).sum()
            elif name == "ccorr_normed":
                ref = (Td * Wd).sum() / np.sqrt((Td ** 2).sum() * (Wd ** 2).sum())
            elif name == "ccoeff":
                ref = ((Td - Td.mean()) * (Wd - Wd.mean())).sum()
            else:
                a, b = Td - Td.mean(), Wd - Wd.mean()
                d = np.sqrt((a * a).sum() * (b * b).sum())
                ref = (a * b).sum() / d if d > 0 else 0.0
            errs.append(abs(float(res[r, c]) - ref) / (abs(ref) + 1e-12))
        worst = float(np.max(errs))
        # Gate 3: OpenCV uses a DFT path for large templates, so the unnormalised methods
        # -- whose values are sums, not ratios -- are the ones worth checking for float32
        # cancellation. Measured they are in fact the *more* precise: ~1e-7 against ~1e-3
        # for TM_CCOEFF_NORMED near a vanishing denominator.
        rows.append({"gate": "matchtemplate_precision", "bank": "-", "method": name,
                     "detail": f"median rel err {np.median(errs):.2e}, worst {worst:.2e}",
                     "passed": bool(worst < 1e-2)})

        # Gate 4: the sign flip and the robust-z are monotone, so they must not reorder
        # anything. If either did, every read-depth in the run would be measuring a
        # different ranking than the one named in the CSV.
        from scipy.stats import spearmanr
        sub = np.asarray(res[::37, ::37], dtype=np.float64).ravel()
        sub = sub[np.isfinite(sub)]
        flipped = tm.method_sign(m) * sub
        med = np.median(flipped)
        mad = 1.4826 * np.median(np.abs(flipped - med))
        z = (flipped - med) / mad if mad > 0 else flipped
        rho = spearmanr(flipped, z).statistic if len(np.unique(flipped)) > 1 else 1.0
        rows.append({"gate": "z_is_monotone", "bank": "-", "method": name,
                     "detail": f"spearman(sign*raw, z) = {rho:.6f} over {len(sub)} points",
                     "passed": bool(abs(rho - 1.0) < 1e-9)})

    rows.append(_gate_regression_tie())
    return pd.DataFrame(rows)


def _gate_regression_tie(fn="002.tiff") -> dict:
    """Gate 7: the refactor must not have moved `TM_CCOEFF_NORMED` by one pixel.

    Reconstructs the pre-refactor `fused_response` inline and requires identical scores on
    reachable pixels and an identical post-NMS detection list. Without this, exp1's whole
    claim -- that removing the gain normalisation changed something -- could be an artefact of
    the sentinel and signature changes made to accommodate the other five methods.
    """
    rgb = ds.load_roi(f"images/{fn}")
    img = ch.to_channel(rgb, "gray_inverted")
    ann = ds.load_annotations()[1]
    seed = ds.image_annotations(ann, fn).query("category_id == 1").iloc[0]
    patch = tm.read_padded_patch(img, seed["cx"], seed["cy"], tm.PATCH_SIZE)
    templates, _ = tm.build_augmentations(patch, tm.BASE_SIZE, (1.0,), 1, (False,))

    h, w = img.shape[:2]
    old = np.full((h, w), np.float32(-2.0), np.float32)
    valid_o = np.zeros((h, w), bool)
    for tpl in templates:
        th, tw = tpl.shape[:2]
        r = cv2.matchTemplate(np.ascontiguousarray(img, np.float32), tpl, cv2.TM_CCOEFF_NORMED)
        np.nan_to_num(r, copy=False, nan=-2.0, posinf=-2.0, neginf=-2.0)
        oy, ox = (th - 1) // 2, (tw - 1) // 2
        view = old[oy:oy + r.shape[0], ox:ox + r.shape[1]]
        np.copyto(view, r, where=r > view)
        valid_o[oy:oy + r.shape[0], ox:ox + r.shape[1]] = True

    new, _, valid_n = tm.fused_response(img, templates)
    radius = ev.radius_px(ds.roi_mpp(f"images/{fn}"))
    co, so = tm.extract_peaks(old, valid_o, 7, 0.5, 250000)
    cn, sn = tm.extract_peaks(new, valid_n, 7, 0.5, 250000)
    ko = co[nms_by_distance(co, so, radius)]
    kn = cn[nms_by_distance(cn, sn, radius)]
    ok = (np.array_equal(valid_o, valid_n) and np.array_equal(old[valid_o], new[valid_n])
          and np.array_equal(ko, kn))
    return {"gate": "regression_tie_ccoeff_normed", "bank": "single", "method": "ccoeff_normed",
            "detail": f"{fn}: {len(ko)} detections, scores and centres identical to "
                      "the pre-refactor implementation", "passed": bool(ok)}


# --- clamp gate (gate 5), applied to the grid output ----------------------------------------
CLAMP_MAX_FRAC = 0.05   # pre-committed before the grid existed; see the research log


def clamp_gate(df: pd.DataFrame) -> pd.DataFrame:
    """Gate 5: mark any arm whose ranking is mostly one tie block as arbitrary, not measured.

    `TM_SQDIFF_NORMED` is clamped to [0, 1] by OpenCV -- verified against the hand formula,
    which it matches wherever that is below 1 and saturates wherever it exceeds it. The clamp
    lands at the *worst* end of the ranking, which is exactly where `read_100` reads. Measured
    on candidate centroids it is not benign: 19.9% of 301.tiff's 19919 competing candidates sit
    at exactly 1.0, against 0.26% on 246.tiff.

    The 5% threshold is pre-committed rather than chosen once the numbers landed. It applies to
    every arm, not just the SQDIFF one -- an unexpected tie block anywhere means the same
    thing.
    """
    frac = df["largest_tie_block"] / df["n_detections"].replace(0, np.nan)
    return df.assign(tie_frac=frac.round(4), tail_arbitrary=(frac > CLAMP_MAX_FRAC).fillna(False))


def report(df: pd.DataFrame, checks, verif, stem: str, pairs=()):
    """Write the grid, its seed spread, any paired deltas, and the verification evidence."""
    RESULTS.mkdir(exist_ok=True)
    df = clamp_gate(df)
    df.to_csv(RESULTS / f"{stem}.csv", index=False)
    seed_variance(df).to_csv(RESULTS / f"{stem}_seed_variance.csv", index=False)
    for a, b in pairs:
        paired_delta(df, a, b).to_csv(RESULTS / f"{stem}_paired_{a}_vs_{b}.csv", index=False)

    ev_rows = [dict(c) for c in checks] + [dict(v) for v in verif]
    pd.DataFrame(ev_rows).to_csv(RESULTS / f"{stem}_verification.csv", index=False)
    failed = [r for r in ev_rows if not r.get("passed", True)]
    arb = df.loc[df["tail_arbitrary"], "arm"].value_counts()
    print(f"\n{stem}: {len(df)} rows, {len(ev_rows)} checks, {len(failed)} failed")
    for f in failed:
        print("  FAIL", f)
    if len(arb):
        print(f"  tail_arbitrary (tie block > {CLAMP_MAX_FRAC:.0%} of the pool):")
        for arm, n in arb.items():
            print(f"    {arm}: {n} rows")
    return df


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gates", action="store_true", help="gates only, no grid")
    ap.add_argument("--smoke", action="store_true", help="one ROI, one seed -- wiring only")
    ap.add_argument("--exp1", action="store_true", help="TM_CCOEFF vs TM_CCOEFF_NORMED")
    ap.add_argument("--exp2", action="store_true", help="all six methods plus controls")
    ap.add_argument("--exp3", nargs=2, metavar="METHOD",
                    help="two methods to re-run with rot4 x flip2")
    ap.add_argument("--rois", nargs="*", help="restrict to these file names")
    args = ap.parse_args(argv)

    if args.gates:
        g = run_gates()
        RESULTS.mkdir(exist_ok=True)
        g.to_csv(RESULTS / "tm_variant_gates.csv", index=False)
        print(g.to_string(index=False))
        bad = int((~g["passed"]).sum())
        print(f"\n{len(g)} gates, {bad} failed")
        return 1 if bad else 0

    rois = roi_table()
    if args.rois:
        rois = rois[rois["file_name"].isin(args.rois)]

    if args.smoke:
        rois = rois[rois["file_name"] == SMOKE_ROI]
        cfgs = [("tm_variant_smoke", ALL_METHODS, 1, 1, (False,), True, True)]
    elif args.exp1:
        cfgs = [("tm_ccoeff_headtohead", EXP1_METHODS, N_SEEDS, 1, (False,), True, True)]
    elif args.exp2:
        cfgs = [("tm_variant_sweep", ALL_METHODS, N_SEEDS, 1, (False,), True, False)]
    elif args.exp3:
        cfgs = [("tm_variant_stage_b", tuple(args.exp3), N_SEEDS, 4, (False, True), False, False)]
    else:
        ap.error("pick one of --gates / --smoke / --exp1 / --exp2 / --exp3")

    for stem, methods, n_seeds, na, fl, controls, shipped in cfgs:
        t0 = time.time()
        frames, checks, verif = [], [], []
        for _, row in rois.iterrows():
            f, c, v = run_roi(row, list(methods), n_seeds=n_seeds, n_angles=na, flips=fl,
                              controls=controls, shipped=shipped)
            frames.append(f)
            checks += c
            verif += v
        df = pd.concat(frames, ignore_index=True)
        pairs = [(f"tm_{methods[0]}", f"tm_{methods[1]}")] if len(methods) == 2 else []
        report(df, checks, verif, stem, pairs=pairs)
        print(f"{stem}: {time.time() - t0:.0f}s over {len(rois)} ROIs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
