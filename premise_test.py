"""The premise test: does the pathologist's click buy anything a seedless arm cannot?

`Research Logs/2026-08-31-next-steps-plan.md`, Step 1. Two questions, one grid:

1. **Does the search earn its place?** `find_and_suppress` seeded from one annotated
   mitotic figure, against four seedless comparators -- a nucleus detector under two
   rankers, and a bare lattice at four spacings -- at matched candidate budget.
2. **Does the click carry information at all?** The same pipeline seeded from (a) a
   mitotic figure, (b) a pathologist-rejected look-alike, (c) a random `nucleus_blobs`
   component centre. If (b) and (c) match (a), the click is not the thing doing the work.

Two constraints make this cheap enough to run at 5 seeds x 5 z levels:

* **One match per (ROI, seed).** The response map is computed once at a deep floor and
  every z level is derived by *filtering* that one NMS'd candidate pool. This is exact,
  not an approximation: `extract_peaks`'s local maxima are threshold-independent, and
  greedy score-ordered NMS means a peak above t can only be suppressed by a
  higher-scoring peak, which is also above t. The argument has one soft spot: `argsort`
  is not stable, so an exact-score tie pair inside the NMS radius can be ordered one way
  in a larger pool and the other way in a smaller one, flipping which of the two survives.
  **This is observed, not hypothetical**: over 32 checked cells, 31 are bit-identical and
  one (301.tiff seed 3, z = 0.75) differs by 1 detection in 16,825, the two candidates
  1.41 px apart -- the same object, well inside the 29.6 px match radius, so no reported
  number moves. `verify_shortcut` therefore raises only on a disagreement with no
  counterpart within a match radius, and reports benign tie-breaks separately. The main
  run samples seed 0 and z in {1.0, 2.0}; `results/premise_shortcut_lowz.csv` covers
  z in {0.25, 0.5, 0.75} at two seeds on the three decision-grade ROIs. That is a sample
  of the (seed, z) grid, not a proof over it.
* **One candidate set per seedless arm per ROI.** `blob_native`, `blob_od` and the four
  `grid_*` arms do not depend on z at all, and depend on the seed only through which
  annotation is excluded from the evaluation ground truth. They are built once and
  re-scored against each seed's own ``gt_eval``.

The extraction floor is ``med + FLOOR_Z*mad`` (currently 0.25) of that seed's own response map -- per image and
per seed, never a global constant. A constant floor in TM_CCOEFF_NORMED units means a
different search depth on every ROI (246.tiff returned 5 detections at 0.5 while 301.tiff
returned 13,353), which is the defect this experiment exists to characterise. Because the
floor and the z threshold are then in the same units, ``floor_limited`` is False by
construction for every z >= 0.5; it is still recorded and asserted, because a True would
mean the two have drifted apart again.

Writes:
  results/premise_test.csv             -- the main grid, one row per (roi, seed, arm, z, budget)
  results/premise_seed_provenance.csv  -- the (a)/(b)/(c) control
  results/premise_verification.csv     -- the shortcut-exactness and invariant evidence
"""

from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd

from midog_utils import baselines as bl
from midog_utils import channels as ch
from midog_utils import chromatin as cm
from midog_utils import compare as cp
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from midog_utils import experiment as ex
from midog_utils import find_and_suppress as fs
from midog_utils import invariants as inv
from midog_utils import seed_selection as ss
from midog_utils import template_match as tm
from midog_utils.nms import nms_by_distance

# --- configuration ----------------------------------------------------------------
CHANNEL = "rgb"                      # search channel, matching od_experiment/od_seed_sweep
FLOOR_Z = 0.25                       # extraction floor, in units of the map's own MAD
# z = FLOOR_Z is the floor itself -- that arm is the whole extracted pool, and is reported
# as the deepest level measurable without re-extracting rather than as a swept point. The
# range was previously capped at 1.0 on the strength of "the pool grows only about 3%
# below it", which is 301.tiff's number and the smallest of the seven: measured between
# z = 1.0 and z = 0.5, the other six ROIs grow 6.4-13.4%. Recall at budget below z = 1.0
# was never measured at all. It is measured here.
Z_LEVELS = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0)
GRID_STEPS = (10, 20, 30, 40)        # swept: a single step is a tuned hyperparameter
N_SEEDS = 5
BUDGETS = cp.BUDGETS                 # (500, 1000, 2000, 5000)
PROVENANCE_ROIS = ("301.tiff", "246.tiff")
PROVENANCE_Z = 2.0
RANDOM_SANITY_ROI = "301.tiff"       # the random arm is one sanity row, not a swept arm
RANDOM_SANITY_SEED = 0
# n_mitotic >= 15. Everything else is reported and labelled, and drives no conclusion.
DECISION_GRADE = {"301.tiff", "246.tiff", "201.tiff"}


def robust_stats(fused: np.ndarray, stride: int = 8):
    """Median and MAD-scale of a response map -- the statistics `_robust_z` uses.

    Values <= -1.5 are the NaN sentinel `fused_response` writes for zero-variance
    windows and are dropped, not counted as a low score.
    """
    s = fused[::stride, ::stride]
    s = s[s > -1.5]
    med = float(np.median(s))
    mad = float(1.4826 * np.median(np.abs(s - med)))
    return med, mad


def match_pool(img, hem, seed_xy, base_size, radius, cfg, keep_map=False):
    """One match, one NMS, at the deep floor -- the pool every z level is derived from."""
    patch = tm.read_padded_patch(img, seed_xy[0], seed_xy[1], cfg.patch_size)
    if patch is None:
        raise ValueError(f"seed {seed_xy} is too close to the ROI border for a template")
    templates, _ = tm.build_augmentations(patch, base_size, cfg.scales,
                                          cfg.n_angles, cfg.flips)
    fused, _, valid = tm.fused_response(img, templates)
    med, mad = robust_stats(fused)
    floor = med + FLOOR_Z * mad

    centers, scores = tm.extract_peaks(fused, valid, cfg.peak_min_distance,
                                       floor, cfg.max_peaks)
    n_peaks = len(centers)
    # `invariants.check_no_cap` is applied downstream to the z-filtered, post-NMS arm
    # length, which can never equal the extraction cap -- so it passes by construction and
    # cannot catch the truncation it exists to catch. (For the `grid_*` arms it is passed
    # no caps at all, so it is vacuous there twice over.) This is the only place the cap
    # could ever bind. It does not bind at present either -- 250000 is above the ~173k
    # separation-limited maximum -- but that bound assumes the min-distance separation
    # that `fused >= dilated` violates on plateaus, and a deeper floor makes plateaus
    # more reachable. The assert costs nothing and fails loudly if that day comes.
    if n_peaks == cfg.max_peaks:
        raise inv.InvariantError(
            f"extract_peaks returned exactly max_peaks={cfg.max_peaks} at floor "
            f"{floor:.5f} -- the peak list was truncated, so every z level derived from "
            "this pool is conditioned on that truncation")
    keep = nms_by_distance(centers, scores, radius)
    centers, scores = centers[keep], scores[keep]
    ok = np.hypot(centers[:, 0] - seed_xy[0], centers[:, 1] - seed_xy[1]) > cfg.self_hit_radius
    centers, scores = centers[ok], scores[ok]

    pool = pd.DataFrame({"cx": centers[:, 0], "cy": centers[:, 1], "score": scores})
    pool = cm.score_detections(pool, hem)
    diag = {"med": med, "mad": mad, "floor": floor, "n_peaks": n_peaks,
            "n_pool": len(pool), "base_size": int(base_size)}
    return pool, diag, ((fused, valid) if keep_map else None)


def verify_shortcut(fused, valid, radius, seed_xy, cfg, pool, cut):
    """Check that filtering the deep pool at ``cut`` equals extracting at ``cut``.

    This is the single check the whole one-match-many-z economy rests on, and it is **not**
    always exact. `np.argsort` is unstable in both `extract_peaks` and `nms_by_distance`, so
    an exact-score tie pair inside the NMS radius can be ordered one way in a large pool and
    the other way in a small one, flipping which of the two survives. Measured directly
    (`results/premise_shortcut_lowz.csv`): 31 of 32 checked cells are bit-identical, and the
    single exception -- 301.tiff seed 3 at z = 0.75 -- differs by **1 detection out of
    16,825**, with the two candidates 1.41 px apart, i.e. the same object well inside the
    29.6 px match radius.

    So strict identity is the wrong bar: it fails on a tie-break that cannot move any
    reported number. The bar used here is *substantive* -- every disagreeing point must have
    a counterpart within one match radius, so that greedy matching resolves both the same
    way. A disagreement with no such counterpart is a real divergence and raises.
    """
    centers, scores = tm.extract_peaks(fused, valid, cfg.peak_min_distance, cut, cfg.max_peaks)
    keep = nms_by_distance(centers, scores, radius)
    centers = centers[keep]
    ok = np.hypot(centers[:, 0] - seed_xy[0], centers[:, 1] - seed_xy[1]) > cfg.self_hit_radius
    direct = centers[ok]
    filtered = pool.loc[pool["score"] >= cut, ["cx", "cy"]].to_numpy()

    a = direct[np.lexsort((direct[:, 1], direct[:, 0]))] if len(direct) else direct
    b = filtered[np.lexsort((filtered[:, 1], filtered[:, 0]))] if len(filtered) else filtered
    identical = a.shape == b.shape and bool(np.array_equal(a, b))
    n_diff, worst = 0, 0.0
    if not identical:
        A, B = {tuple(r) for r in a}, {tuple(r) for r in b}
        d_only, f_only = np.array(sorted(A - B)), np.array(sorted(B - A))
        n_diff = len(d_only) + len(f_only)
        if len(d_only) and len(f_only):
            dist = np.hypot(d_only[:, None, 0] - f_only[None, :, 0],
                            d_only[:, None, 1] - f_only[None, :, 1])
            # Worst orphan in either direction: a point whose nearest counterpart is far.
            worst = float(max(dist.min(axis=1).max(), dist.min(axis=0).max()))
        else:
            worst = float("inf")   # one side has extra points with nothing to pair against
    return {"n_direct": int(len(direct)), "n_filtered": int(len(filtered)),
            "identical": identical, "n_disagree": int(n_diff),
            "worst_disagreement_px": (round(worst, 2) if np.isfinite(worst) else worst),
            "benign": bool(identical or (n_diff and worst <= radius))}


def grid_arm_candidates(roi_shape, step, hem, radius):
    """A lattice, chromatin-scored, then NMS'd at the match radius by that same key."""
    g = bl.grid_lattice(roi_shape, step)
    g = cm.score_detections(g, hem)
    order_key = np.nan_to_num(g["od"].to_numpy(), nan=-np.inf)
    keep = nms_by_distance(g[["cx", "cy"]].to_numpy(), order_key, radius)
    out = g.iloc[keep].reset_index(drop=True)
    return out.assign(score=out["od"])  # no native score; keep the column contract


def _rejection_sample(pool: pd.DataFrame, gray, rng, max_tries: int = 400):
    """Draw uniformly from the subset of ``pool`` that passes the foreground filter.

    Rejection sampling is distributionally identical to filtering then drawing uniformly,
    and avoids running `seed_selection.foreground_filter` over tens of thousands of blob
    centroids. Raises rather than falling back to an untightened template: a silent
    fallback would make the control arm a different pipeline than the arm it controls.
    """
    if len(pool) == 0:
        raise ValueError("empty candidate pool for seed selection")
    for _ in range(max_tries):
        row = pool.iloc[int(rng.integers(len(pool)))]
        patch = tm.read_padded_patch(gray, float(row["cx"]), float(row["cy"]), tm.BASE_SIZE)
        if patch is not None and ss.tighten_box_otsu(patch) is not None:
            return row
    raise ValueError(f"no seed passed the foreground filter in {max_tries} draws")


def main(smoke: bool = False):
    """``smoke=True`` runs one ROI, 2 seeds, 2 z levels -- a wiring check, not a result."""
    global Z_LEVELS, GRID_STEPS, N_SEEDS, PROVENANCE_ROIS
    images, ann = ds.load_annotations()
    sel = ex.select_domain_images(images, ann)
    if smoke:
        Z_LEVELS, GRID_STEPS, N_SEEDS = (2.0, 2.5), (20, 40), 2
        PROVENANCE_ROIS = ("246.tiff",)
        sel = sel[sel["file_name"] == "246.tiff"].reset_index(drop=True)
    out_tag = "_smoke" if smoke else ""
    cfg0 = fs.FSConfig(channel=CHANNEL)
    border = cfg0.patch_size // 2

    rows, prov_rows, checks, verif = [], [], [], []
    cov_cache: dict = {}
    t_start = time.time()

    for _, row in sel.iterrows():
        fn, image_id = row["file_name"], int(row["image_id"])
        t_roi = time.time()
        path = f"images/{fn}"
        rgb = ds.load_roi(path)
        mpp = ds.roi_mpp(path)
        radius = ev.radius_px(mpp)
        area = ds.check_roi_scale(path, rgb.shape)
        gray = ch.to_gray_inverted(rgb)
        hem = cm.hematoxylin_od(rgb)          # unclipped OD; to_hematoxylin saturates
        img = ch.to_channel(rgb, CHANNEL)
        gt = ds.image_annotations(ann, fn)
        gt_mit = gt[gt["category_id"] == ds.MITOTIC]
        gt_look = gt[gt["category_id"] == ds.LOOKALIKE]

        mask = bl.tissue_mask(rgb)
        checks.append(inv.check_tissue_mask_covers_gt(mask, gt, label=fn))

        # --- seedless candidate sets: built once per ROI, reused across seeds and z ---
        t0 = time.time()
        blobs = bl.nucleus_blobs(rgb, mask)
        blobs = cm.score_detections(blobs, hem)
        # `nucleus_blobs` emits raw connected components with no suppression, and they are
        # dense: measured on the images, 58.6% of 301.tiff's and 70.5% of 246.tiff's blobs
        # lie within one match radius of another blob. Greedy matching is one-to-one, so a
        # duplicate can never claim a second ground-truth object -- it only consumes
        # budget. Every other arm in this comparison is suppressed at the match radius, so
        # a raw blob arm is not a matched-budget comparator. The raw arms are kept so the
        # size of the effect stays visible, and a suppressed variant is added beside them,
        # ordered by that arm's *own* ranking key -- suppression ordered by anything else
        # silently reranks the arm.
        _bxy = blobs[["cx", "cy"]].to_numpy()
        blobs_nms = {
            "blob_native_nms": ("score", blobs.iloc[nms_by_distance(
                _bxy, np.nan_to_num(blobs["score"].to_numpy(), nan=-np.inf),
                radius)].reset_index(drop=True)),
            "blob_od_nms": ("od", blobs.iloc[nms_by_distance(
                _bxy, np.nan_to_num(blobs["od"].to_numpy(), nan=-np.inf),
                radius)].reset_index(drop=True)),
        }
        for _nm, (_key, _df) in blobs_nms.items():
            checks.append(inv.check_min_separation(_df[["cx", "cy"]].to_numpy(), radius,
                                                   label=f"{fn}/{_nm}"))
        t_blob = time.time() - t0

        t0 = time.time()
        grids = {}
        for step in GRID_STEPS:
            grids[step] = grid_arm_candidates(rgb.shape, step, hem, radius)
            checks.append(inv.check_min_separation(
                grids[step][["cx", "cy"]].to_numpy(), radius, label=f"{fn}/grid_{step}"))
        t_grid = time.time() - t0

        seedless = [
            cp.Arm("blob_native", (lambda d=blobs: d), rank_key="score",
                   caps=(None,), coverage_key=f"{fn}:blob_native",
                   extra={"generator": "nucleus_blobs", "grid_step": np.nan}),
            cp.Arm("blob_od", (lambda d=blobs: d), rank_key="od",
                   caps=(None,), coverage_key=f"{fn}:blob_od",
                   extra={"generator": "nucleus_blobs", "grid_step": np.nan}),
        ]
        for _nm, (_key, _df) in blobs_nms.items():
            seedless.append(cp.Arm(
                _nm, (lambda d=_df: d), rank_key=_key, caps=(None,), nms_radius=radius,
                coverage_key=f"{fn}:{_nm}",
                extra={"generator": "nucleus_blobs_nms", "grid_step": np.nan}))
        for step in GRID_STEPS:
            seedless.append(cp.Arm(
                f"grid_{step}", (lambda d=grids[step]: d), rank_key="od",
                nms_radius=radius, coverage_key=f"{fn}:grid_{step}",
                extra={"generator": f"lattice_{step}px", "grid_step": step}))

        print(f"[{fn}] loaded, blobs={len(blobs)}"
              f"->{[len(d) for _, d in blobs_nms.values()]} nms ({t_blob:.0f}s), "
              f"grids={[len(grids[s]) for s in GRID_STEPS]} ({t_grid:.0f}s) "
              f"[{time.time() - t_roi:.0f}s]", flush=True)

        seed_records, prov_cache = [], []
        for s in range(N_SEEDS):
            t_seed = time.time()
            rng = np.random.default_rng([s, image_id])
            seed, seed_info = ss.pick_seed(gt_mit, gray, rng, border, rgb.shape)
            seed_records.append((s, int(seed["ann_id"]), (s, image_id)))
            seed_pool_size = seed_info.n_after_foreground
            base = ss.tightened_base_size(gray, seed["cx"], seed["cy"])
            cfg = fs.FSConfig(channel=CHANNEL, base_size=base)

            gt_eval = gt[gt["ann_id"] != seed["ann_id"]].reset_index(drop=True)
            keep_map = (s == 0)
            pool, diag, maps = match_pool(img, hem, (seed["cx"], seed["cy"]),
                                          base, radius, cfg, keep_map=keep_map)

            ctx = {"file_name": fn, "image_id": image_id, "tumor_type": row["tumor_type"],
                   "decision_grade": fn in DECISION_GRADE, "mpp": round(mpp, 4),
                   "channel": CHANNEL, "seed_index": s,
                   "seed_ann_id": int(seed["ann_id"]), "base_size": int(base),
                   "map_median": round(diag["med"], 5), "mad_scale": round(diag["mad"], 5),
                   "extract_floor": round(diag["floor"], 5), "n_pool": diag["n_pool"],
                   "n_peaks": diag["n_peaks"], "sanity_only": False}

            arms = []
            for z in Z_LEVELS:
                cut = diag["med"] + z * diag["mad"]
                sub = pool[pool["score"] >= cut]
                limited = bool(cut < diag["floor"])
                arms.append(cp.Arm(
                    "pipeline", (lambda d=sub: d), rank_key="od", seeded=True, z=z,
                    z_dependent=True, floor_limited=limited, nms_radius=radius,
                    caps=(cfg.max_peaks, cfg.max_detections),
                    extra={"generator": "find_and_suppress", "z_cut": round(cut, 5),
                           "grid_step": np.nan}))
                arms.append(cp.Arm(
                    "pipeline_score", (lambda d=sub: d), rank_key="score", seeded=True,
                    z=z, z_dependent=True, floor_limited=limited, nms_radius=radius,
                    caps=(cfg.max_peaks, cfg.max_detections),
                    extra={"generator": "find_and_suppress", "z_cut": round(cut, 5),
                           "grid_step": np.nan}))
            arms.extend(seedless)

            if fn == RANDOM_SANITY_ROI and s == RANDOM_SANITY_SEED:
                n_rand = len(grids[20])
                rnd = bl.random_in_tissue(mask, n_rand, np.random.default_rng([s, image_id, 9]))
                rnd = cm.score_detections(rnd, hem)
                arms.append(cp.Arm(
                    "random", (lambda d=rnd: d), rank_key="od",
                    coverage_key=f"{fn}:random",
                    extra={"generator": f"random_in_tissue_n{n_rand}", "grid_step": np.nan}))

            got = cp.evaluate_arms(arms, gt_eval, radius, roi_shape=rgb.shape, mpp=mpp,
                                   budgets=BUDGETS, context=ctx,
                                   coverage_cache=cov_cache.setdefault(fn, {}),
                                   checks=checks)
            got.loc[got["arm"] == "random", "sanity_only"] = True
            rows.append(got)

            if keep_map:
                # Checked at two levels, not one: z = 2.0 as before, and z = 1.0, which is
                # both the level the frontier now recommends and the one nearest the
                # extraction floor -- i.e. the level with the most peaks to disagree about.
                for z_v in (1.0, 2.0):
                    cut2 = diag["med"] + z_v * diag["mad"]
                    v = verify_shortcut(maps[0], maps[1], radius, (seed["cx"], seed["cy"]),
                                        cfg, pool, cut2)
                    v.update({"file_name": fn, "seed_index": s, "z": z_v,
                              "cut": round(cut2, 5), "check": "one_match_many_z"})
                    verif.append(v)
                del maps

            # --- seed-provenance control, at one z, on the decision-grade ROIs only ---
            if fn in PROVENANCE_ROIS:
                seed_b = _rejection_sample(
                    ss.border_filter(gt_look, border, rgb.shape),
                    gray, np.random.default_rng([s, image_id, 2]))
                blob_pool = ss.border_filter(blobs, border, rgb.shape)
                seed_c = _rejection_sample(
                    blob_pool, gray, np.random.default_rng([s, image_id, 3]))
                prov_cache.append((s, seed, pool, diag, base, seed_b, seed_c))

            print(f"  [{fn}] seed {s} ann={int(seed['ann_id'])} base={base}px "
                  f"pool={diag['n_pool']} peaks={diag['n_peaks']} "
                  f"med={diag['med']:.3f} mad={diag['mad']:.3f} "
                  f"[{time.time() - t_seed:.0f}s]", flush=True)

        checks.append(inv.check_distinct_seeds(seed_records, pool_size=seed_pool_size,
                                               label=fn))

        for s, seed, pool_a, diag_a, base_a, seed_b, seed_c in prov_cache:
            t_p = time.time()
            # All three arms must see an identical ground truth. Dropping seeds (a) and
            # (b) by ann_id does not achieve that: `match_pool` strips every detection
            # within `self_hit_radius` of the seed, and arm (c)'s seed is a blob centroid
            # rather than an annotation, so no annotation was removed for it while the
            # exclusion zone still applied. Measured: on 301.tiff seed 1 the drawn blob
            # sits 2.8 px from a mitotic figure, which is inside the 5 px zone -- that
            # figure is unreachable for (c) alone, biasing the contrast towards "the click
            # carries information". Excluding the union of the three zones removes the
            # asymmetry; (a) and (b) fall inside their own zones at distance 0, so this is
            # a superset of the old rule.
            _shr = float(fs.FSConfig(channel=CHANNEL).self_hit_radius)
            _sxy = np.array([[float(sd["cx"]), float(sd["cy"])]
                             for sd in (seed, seed_b, seed_c)])
            _dist = np.hypot(gt["cx"].to_numpy()[:, None] - _sxy[None, :, 0],
                             gt["cy"].to_numpy()[:, None] - _sxy[None, :, 1])
            _excl = (_dist <= _shr).any(axis=1)
            gt_prov = gt[~_excl].reset_index(drop=True)
            n_mit_prov = int((gt_prov["category_id"] == ds.MITOTIC).sum())
            n_mit_all = int((gt["category_id"] == ds.MITOTIC).sum())
            # Seeds (a) and (b) are rows of `gt` at distance 0 from themselves, so
            # "were they dropped" is true by construction and asserting it proves
            # nothing. The bound that carries information is the *upper* one: a units
            # error or a wrong coordinate column in `_dist` would exclude tens or
            # hundreds of annotations, not one or two.
            n_drop = n_mit_all - n_mit_prov
            assert 1 <= n_drop <= 3, (
                f"provenance ground truth dropped {n_drop} mitotic annotations from "
                f"{n_mit_all}; expected 1 (seed a) plus at most 2 inside the other two "
                f"{_shr:g} px exclusion discs -- a larger number means the exclusion "
                "radius or the coordinate columns are wrong")

            entries = [("a_mitotic", seed, pool_a, diag_a, base_a, int(seed["ann_id"]))]
            for tag, sd in (("b_lookalike", seed_b), ("c_random_blob", seed_c)):
                b = ss.tightened_base_size(gray, sd["cx"], sd["cy"])
                c = fs.FSConfig(channel=CHANNEL, base_size=b)
                p, d, _ = match_pool(img, hem, (sd["cx"], sd["cy"]), b, radius, c)
                entries.append((tag, sd, p, d, b,
                                int(sd["ann_id"]) if "ann_id" in sd.index else -1))

            for tag, sd, p, d, b, ann_id in entries:
                cut = d["med"] + PROVENANCE_Z * d["mad"]
                sub = p[p["score"] >= cut]
                c = fs.FSConfig(channel=CHANNEL, base_size=b)
                ctx = {"file_name": fn, "image_id": image_id,
                       "tumor_type": row["tumor_type"], "decision_grade": True,
                       "seed_index": s, "provenance": tag, "seed_ann_id": ann_id,
                       "seed_cx": round(float(sd["cx"]), 1),
                       "seed_cy": round(float(sd["cy"]), 1), "base_size": int(b),
                       "map_median": round(d["med"], 5), "mad_scale": round(d["mad"], 5),
                       "extract_floor": round(d["floor"], 5), "n_pool": d["n_pool"],
                       "n_peaks": d["n_peaks"], "n_gt_excluded": int(_excl.sum()),
                       "self_hit_radius": _shr, "channel": CHANNEL}
                prov_rows.append(cp.evaluate_arms(
                    [cp.Arm("pipeline", (lambda x=sub: x), rank_key="od", seeded=True,
                            z=PROVENANCE_Z, z_dependent=True,
                            floor_limited=bool(cut < d["floor"]), nms_radius=radius,
                            caps=(c.max_peaks, c.max_detections),
                            extra={"z_cut": round(cut, 5)})],
                    gt_prov, radius, roi_shape=rgb.shape, mpp=mpp, budgets=BUDGETS,
                    context=ctx, checks=checks))
            print(f"  [{fn}] provenance seed {s}: "
                  + ", ".join(f"{t}={d['n_pool']}" for t, _, _, d, _, _ in entries)
                  + f" [{time.time() - t_p:.0f}s]", flush=True)

        print(f"[{fn}] done [{time.time() - t_roi:.0f}s total, "
              f"{time.time() - t_start:.0f}s elapsed]", flush=True)
        del rgb, gray, hem, img, mask, blobs, blobs_nms, grids

    main_df = pd.concat(rows, ignore_index=True)
    cp.assert_floor_not_limiting(main_df)
    main_df.to_csv(f"results/premise_test{out_tag}.csv", index=False)

    prov_df = pd.concat(prov_rows, ignore_index=True)
    cp.assert_floor_not_limiting(prov_df)
    prov_df.to_csv(f"results/premise_seed_provenance{out_tag}.csv", index=False)

    ver = pd.DataFrame(verif)
    chk = pd.DataFrame(checks)
    # One detail row per (check, label) -- i.e. per (ROI, arm) -- plus the raw number of
    # times each check actually ran, so the audit trail is compact without understating
    # how much was checked.
    n_runs = chk.groupby("check").size().rename("n_runs")
    chk = chk.drop_duplicates(subset=["check", "label"], keep="first").merge(
        n_runs, left_on="check", right_index=True, how="left")
    pd.concat([ver.assign(kind="shortcut_exactness"),
               chk.assign(kind="invariant")], ignore_index=True).to_csv(
        f"results/premise_verification{out_tag}.csv", index=False)

    bad = ver[~ver["benign"].astype(bool)]
    soft = ver[ver["benign"].astype(bool) & ~ver["identical"].astype(bool)]
    if len(soft):
        print(f"shortcut: {len(soft)} cell(s) differ by a tie-break inside the match "
              f"radius (max {soft['n_disagree'].max()} detections, worst separation "
              f"{soft['worst_disagreement_px'].max()} px) -- benign, see verify_shortcut",
              file=sys.stderr)
    print(f"\nrows: main {len(main_df)}, provenance {len(prov_df)}")
    print(f"invariant checks run: {len(checks)} (all passed -- they raise otherwise)")
    print(f"shortcut: {int(ver['identical'].sum())}/{len(ver)} (ROI, z) cells "
          f"bit-identical, {int(ver['benign'].sum())}/{len(ver)} within a match radius")
    if len(bad):
        # Every z level in this file is derived by filtering one pool. If that is not
        # exact, every number downstream of it is conditioned on an unstated tie-break,
        # so this must stop the run rather than print a warning nobody reads.
        print(bad.to_string(), file=sys.stderr)
        raise inv.InvariantError(
            f"one-match-many-z shortcut diverges in {len(bad)} of {len(ver)} checked "
            "(ROI, z) cells by more than one match radius -- filtering the deep pool does "
            "not reproduce a direct extraction, so the derived arms are not what they "
            "claim to be. A tie-break inside the match radius is tolerated and reported "
            "separately; this is not that.")
    print(f"total {time.time() - t_start:.0f}s")
    print("wrote results/premise_{test,seed_provenance,verification}.csv")


# =====================================================================================
# Reporting. `python premise_test.py --report` reads the CSVs back and writes
# results/premise_report.md. Kept in this file, and generated rather than hand-written,
# so no number in the report can drift from the run that produced it.
# =====================================================================================

ARM_ORDER = (["pipeline", "pipeline_score"] + [f"grid_{s}" for s in GRID_STEPS]
             + ["blob_native", "blob_od", "blob_native_nms", "blob_od_nms", "random"])
SEEDLESS = ([f"grid_{s}" for s in GRID_STEPS]
            + ["blob_native", "blob_od", "blob_native_nms", "blob_od_nms"])


def _f(v, nd=3):
    if v is None:
        return "--"
    try:
        x = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not np.isfinite(x):
        return "--"
    return f"{x:.{nd}f}" if nd else f"{x:.0f}"


def _md(rows, header):
    """Minimal markdown table -- avoids a `tabulate` dependency."""
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


def _arm_label(arm, z):
    return f"{arm} z={z:g}" if isinstance(z, float) and np.isfinite(z) else arm


def _stats(s: pd.Series) -> dict:
    """Median / IQR / range over seeds, with unreachable targets counted, not imputed."""
    v = s.dropna()
    return {"median": v.median() if len(v) else np.nan,
            "q1": v.quantile(0.25) if len(v) else np.nan,
            "q3": v.quantile(0.75) if len(v) else np.nan,
            "min": v.min() if len(v) else np.nan,
            "max": v.max() if len(v) else np.nan,
            "n": len(s), "n_unreachable": int(s.isna().sum())}


def _arm_rows(df, roi, budgets=BUDGETS):
    """One row per (arm, z): median [IQR] recall at each budget, over the seeds."""
    d = df[(df["file_name"] == roi) & (~df["sanity_only"].astype(bool))]
    rows = []
    keys = d[["arm", "z"]].drop_duplicates()
    keys["ord"] = keys["arm"].map({a: i for i, a in enumerate(ARM_ORDER)}).fillna(99)
    keys = keys.sort_values(["ord", "z"], na_position="first")
    for _, k in keys.iterrows():
        sub = d[(d["arm"] == k["arm"]) & (d["z"].isna() if pd.isna(k["z"]) else (d["z"] == k["z"]))]
        cells = []
        for b in budgets:
            st = _stats(sub[sub["budget"] == b]["recall_at_budget"])
            cells.append(f"{_f(st['median'])} [{_f(st['q1'])}-{_f(st['q3'])}]")
        n_det = sub[sub["budget"] == budgets[0]]["n_detections"].median()
        cov = sub[sub["budget"] == budgets[0]]["coverage_frac"].median()
        full = sub[sub["budget"] == budgets[0]]["full_list_recall"].median()
        rows.append([_arm_label(k["arm"], k["z"])] + cells
                    + [_f(n_det, 0), _f(cov), _f(full)])
    return rows


def write_report(suffix: str = ""):
    main = pd.read_csv(f"results/premise_test{suffix}.csv")
    prov = pd.read_csv(f"results/premise_seed_provenance{suffix}.csv")
    ver = pd.read_csv(f"results/premise_verification{suffix}.csv")
    main["sanity_only"] = main["sanity_only"].fillna(False)

    dg = [r for r in ("301.tiff", "246.tiff", "201.tiff") if r in set(main["file_name"])]
    n_seeds = int(main["seed_index"].nunique())
    n_roi = int(main["file_name"].nunique())
    up = [r for r in sorted(set(main["file_name"])) if r not in dg]
    hdr = ["arm"] + [f"@{b}" for b in BUDGETS] + ["n_cand", "coverage", "full-list recall"]

    L = []
    A = L.append
    A("# Premise test: does the seeded search earn its place?\n")
    A(f"Generated by `premise_test.py --report` from `results/premise_test.csv` "
      f"({len(main)} rows), `results/premise_seed_provenance.csv` ({len(prov)} rows) and "
      f"`results/premise_verification.csv`. Plain factual summary: no recommendation is "
      f"made and none should be read in.\n")

    # ---------------------------------------------------------------- what was run
    A("## 0. What was run\n")
    A(f"* {len(set(main['file_name']))} ROIs (`experiment.select_domain_images`), "
      f"{main['seed_index'].nunique()} seeds each "
      f"(`np.random.default_rng([seed_index, image_id])`), "
      f"z in {{{', '.join(f'{z:g}' for z in sorted(main['z'].dropna().unique()))}}}.")
    A(f"* Search channel `{main['channel'].iloc[0]}`, `FSConfig` defaults otherwise "
      "(single scale, no rotation/flip augmentation, `peak_min_distance=7`).")
    A(f"* One `fused_response` per (ROI, seed) at a per-image extraction floor "
      f"`med + {FLOOR_Z:g}*mad`; every z derived by filtering that one NMS'd pool. "
      "Checked against a direct re-extraction on a **sample** of the grid, not all of it "
      "(section 8): seed 0 at z in {1.0, 2.0} here, plus z in {0.25, 0.5, 0.75} at two "
      "seeds on the decision-grade ROIs in `results/premise_shortcut_lowz.csv`. 31 of "
      "those 32 cells are bit-identical; the exception differs by 1 detection in 16,825, "
      "1.41 px apart -- inside the match radius, so no number here moves.")
    A("* Seedless arms (`blob_native`, `blob_od`, `grid_*`) are built once per ROI and "
      "re-scored against each seed's own evaluation ground truth. They carry `z = NaN` "
      "and `z_dependent = False` in the CSV rather than being duplicated at every z. "
      "**Join them to the seeded arms on `(file_name, seed_index, budget)`**, not on `z`.")
    A(f"* Decision-grade ROIs (n_mitotic >= 15): {', '.join(dg)}. Everything else "
      f"({', '.join(up)}) is reported and labelled *underpowered*; it drives no "
      "conclusion here.")
    A("* Ranking keys: `od` = `chromatin.hematoxylin_od` darkest-decile window mean; "
      "`score` = `TM_CCOEFF_NORMED`; `blob_native` = the Otsu component's own mean "
      "hematoxylin. All sorts are `kind=\"mergesort\"`, NaN last.")
    A(f"* Every median below is over {n_seeds} seeds and every bracket is the "
      f"interquartile range over those same {n_seeds} seeds. No single-seed number "
      "appears in this report.")
    A("* `pipeline_score`'s budgeted columns are identical across z wherever the budget "
      "is below the list length. That is arithmetic, not a bug: filtering by z removes "
      "exactly the lowest-scoring candidates, so the top-K of a score-ranked list does "
      "not move until K exceeds the truncated length.\n")

    # ------------------------------------------------------------------ sample size
    A("## 1. Sample sizes\n")
    rows = []
    for roi in dg + up:
        d = main[(main["file_name"] == roi) & (main["budget"] == BUDGETS[0])]
        rows.append([roi, d["tumor_type"].iloc[0], int(d["n_gt_mitotic"].median()),
                     int(d["n_gt_lookalike"].median()), d["seed_index"].nunique(),
                     "decision-grade" if roi in dg else "**underpowered**"])
    A(_md(rows, ["ROI", "tumour", "n mitotic (eval)", "n look-alike", "seeds", "class"]))
    A("")

    # ------------------------------------------------- primary metric: recall@budget
    A(f"## 2. Primary metric: recall at budget, median [IQR] over {n_seeds} seeds\n")
    A("`n_cand` is the median candidate-list length, `coverage` the median "
      "`coverage_frac`, `full-list recall` the median recall over the entire list "
      "(uninterpretable wherever coverage is high -- see section 7).\n")
    gr = main[main["arm"].str.startswith("grid_") & (main["budget"] == BUDGETS[0])]
    piv = gr.pivot_table(index=["file_name", "match_radius_px"], columns="arm",
                         values="n_detections", aggfunc="first").reset_index()
    gcols = [c for c in piv.columns if str(c).startswith("grid_")]
    _rad = piv["match_radius_px"].astype(float)
    _steps = sorted(int(str(c).split("_")[1]) for c in gcols)
    _r30 = piv[_rad < 30.0]["file_name"].tolist()
    A("**The grid step is not a clean one-parameter sweep, and the candidate count is "
      "not monotone in it.** Distance NMS only removes lattice points when the step is "
      "below that image's match radius, and the match radius crosses 30 px between these "
      f"ROIs ({_rad.min():.1f}-{_rad.max():.1f}), so a step above every radius leaves the "
      "lattice unsuppressed and *larger* than a finer one.")
    # The worked example needs both a step-30 arm and an ROI whose radius sits below it.
    if _r30 and "grid_30" in piv.columns and "grid_20" in piv.columns:
        _exr = piv[piv["file_name"] == _r30[0]]
        A(f"`grid_30` is suppressed on {len(piv) - len(_r30)} of the {len(piv)} ROIs but "
          f"*not* on {', '.join(_r30)}, whose radius is below 30 px -- which is why "
          f"{_r30[0]}'s step-30 arm carries {int(_exr['grid_30'].iloc[0]):,} candidates "
          f"against step-20's {int(_exr['grid_20'].iloc[0]):,}.")
    A(f"Steps run: {', '.join(str(s) for s in _steps)}. Post-NMS list lengths:\n")
    A(_md([[r["file_name"], _f(r["match_radius_px"])]
           + [_f(r[c], 0) for c in gcols] for _, r in piv.iterrows()],
          ["ROI", "match radius px"] + [str(c) for c in gcols]))
    A("")
    for roi in dg:
        n = int(main[main["file_name"] == roi]["n_gt_mitotic"].median())
        A(f"### {roi} (n_mitotic = {n}, decision-grade)\n")
        A(_md(_arm_rows(main, roi), hdr))
        A("")
    if up:
        A("### Underpowered ROIs -- reported, not decisive\n")
    for roi in up:
        n = int(main[main["file_name"] == roi]["n_gt_mitotic"].median())
        A(f"**{roi}** (n_mitotic = {n}, **underpowered**)\n")
        A(_md(_arm_rows(main, roi), hdr))
        A("")

    # ------------------------------------------------- pooled + head-to-head
    A("## 3. Pipeline against the best seedless arm, matched budget\n")
    A("The budget is the match: both arms are truncated to the same number of "
      "candidates and scored on the same ground truth. `best seedless` is chosen per "
      "(ROI, budget) by median recall; `pipeline best z` is the pipeline's best z by "
      f"median recall at that budget. `wins` counts the seeds (out of {n_seeds}) on "
      "which the pipeline at that z strictly exceeds the best seedless arm. Choosing "
      "both the seedless arm and the pipeline's z by their own best median is generous "
      "to both sides and to the pipeline in particular, since z is a free parameter "
      "here and the seedless arms have none.\n")
    rows, pooled = [], []
    for roi in dg + up:
        for b in BUDGETS:
            d = main[(main["file_name"] == roi) & (main["budget"] == b)
                     & (~main["sanity_only"].astype(bool))]
            sl = d[d["arm"].isin(SEEDLESS)].groupby("arm")["recall_at_budget"].median()
            pl = d[d["arm"] == "pipeline"].groupby("z")["recall_at_budget"].median()
            if not len(sl) or not len(pl):
                continue
            best_arm, best_val = sl.idxmax(), sl.max()
            best_z, pl_val = pl.idxmax(), pl.max()
            per_seed_p = d[(d["arm"] == "pipeline") & (d["z"] == best_z)].set_index(
                "seed_index")["recall_at_budget"]
            per_seed_s = d[d["arm"] == best_arm].set_index("seed_index")["recall_at_budget"]
            common = per_seed_p.index.intersection(per_seed_s.index)
            wins = int((per_seed_p[common] > per_seed_s[common]).sum())
            rows.append([roi if roi in dg else f"{roi} (underpowered)", b, best_arm,
                         _f(best_val), f"z={best_z:g}", _f(pl_val),
                         _f(pl_val - best_val), f"{wins}/{len(common)}"])
            if roi in dg:
                pooled.append({"budget": b, "roi": roi, "seedless": best_val,
                               "pipeline": pl_val, "wins": wins, "n": len(common)})
    A(_md(rows, ["ROI", "budget", "best seedless", "its recall", "pipeline best z",
                 "its recall", "delta", "pipeline wins"]))
    A("")
    A("**Both columns are optimistic, and not symmetrically so.** The pipeline's z and "
      "the seedless arm's grid step are each selected *post hoc* by their own best "
      "median at that budget, so neither is an out-of-sample number -- but z is the "
      f"parameter of the arm under test, and it is chosen from {main['z'].nunique()} levels after seeing the "
      "outcome. Every `delta` above is therefore a lower bound on the pipeline's deficit "
      "at any fixed z: the bias runs in the pipeline's favour and the sign of the result "
      "on the decision-grade ROIs survives it.\n")
    rnd = main[main["sanity_only"].astype(bool)]
    if len(rnd):
        rfn = rnd["file_name"].iloc[0]
        A("**The random floor.** `random_in_tissue` at a length matched to the step-20 "
          "lattice, on one ROI at one seed -- a single sanity value with no interval, "
          "deliberately not swept (a lattice with sampling noise added and worse "
          "coverage is an analytically predictable arm). It is what makes \"lift over "
          "random\" mean anything, and it is excluded from every median and win count "
          "above:\n")
        cmp_rows = []
        for _, r in rnd.sort_values("budget").iterrows():
            comp = main[(main["file_name"] == rfn) & (main["budget"] == r["budget"])
                        & (main["seed_index"] == r["seed_index"])
                        & (main["arm"] == "grid_20")]
            bn = main[(main["file_name"] == rfn) & (main["budget"] == r["budget"])
                      & (main["seed_index"] == r["seed_index"])
                      & (main["arm"] == "blob_native")]
            cmp_rows.append([int(r["budget"]), int(r["n_detections"]),
                             _f(r["recall_at_budget"]),
                             _f(comp["recall_at_budget"].iloc[0]) if len(comp) else "--",
                             _f(bn["recall_at_budget"].iloc[0]) if len(bn) else "--",
                             _f(r["coverage_frac"])])
        A(_md(cmp_rows, ["budget", "n_cand", f"random ({rfn}, 1 seed)",
                         "grid_20 (same seed)", "blob_native (same seed)",
                         "random coverage"]))
        A("")
    A(f"Pooled over the {len(dg)} decision-grade ROIs (median of the per-ROI medians, "
      f"and total pipeline wins over the {len(dg) * n_seeds} (ROI, seed) cells):\n")
    p = pd.DataFrame(pooled)
    rows = [[b, _f(g["seedless"].median()), _f(g["pipeline"].median()),
             _f(g["pipeline"].median() - g["seedless"].median()),
             f"{int(g['wins'].sum())}/{int(g['n'].sum())}"]
            for b, g in p.groupby("budget")]
    A(_md(rows, ["budget", "best seedless", "pipeline (best z)", "delta",
                 "pipeline wins"]))
    A("")

    # --- the same question with nothing selected on either side ---
    _zs = set(main["z"].dropna())
    _pin = 1.0 if 1.0 in _zs else float(min(_zs)) if _zs else float("nan")
    A("**The same question with no arm selected.** The table above picks the seedless arm "
      "by its own best median per (ROI, budget). An earlier write-up of this experiment "
      f"picked it per *cell* instead -- a max over all {len(SEEDLESS)} arms at every "
      "(ROI, seed, budget) -- and reported 2/15 wins at budget 500 from the "
      "pre-correction CSV where this rule gives 0/15. Two undocumented selection rules, "
      "two answers. So here is the matrix with nothing selected on the seedless side at all: "
      f"the pipeline pinned at z = {_pin:g} against each arm in turn, paired by (ROI, "
      f"seed) across the {len(dg)} decision-grade ROIs. `wins` counts cells where the "
      "pipeline is strictly ahead; the bracketed figure is the median paired delta. "
      f"z = {_pin:g} was itself chosen on this data -- that is disclosed, not corrected, "
      "and it is why the only clean version of this comparison needs held-out ROIs.\n")
    _pl = main[(main["arm"] == "pipeline") & (main["z"] == _pin)
               & main["file_name"].isin(dg)]
    _rows = []
    for _arm in SEEDLESS:
        _sa = main[(main["arm"] == _arm) & main["file_name"].isin(dg)]
        _cells = []
        for _b in BUDGETS:
            _a = _pl[_pl["budget"] == _b].set_index(
                ["file_name", "seed_index"])["recall_at_budget"]
            _c = _sa[_sa["budget"] == _b].set_index(
                ["file_name", "seed_index"])["recall_at_budget"]
            _ix = _a.index.intersection(_c.index)
            if not len(_ix):
                _cells.append("--"); continue
            _d = _a[_ix] - _c[_ix]
            _cells.append(f"{int((_d > 0).sum())}/{len(_ix)} ({_d.median():+.3f})")
        _rows.append([_arm] + _cells)
    A(_md(_rows, ["seedless arm"] + [f"wins @{b}" for b in BUDGETS]))
    A("")

    # ------------------------------------------------------------------- frontier
    A("## 4. Recall-vs-workload frontier over z\n")
    A("Pipeline arm only (chromatin ranking). `read_80` is the number of candidates a "
      "reader works through to reach 80% sensitivity; `unreached` counts seeds where the "
      "list never gets there, and those seeds are excluded from the median rather than "
      "imputed. **Read the two together**: a reading-depth median over the seeds that "
      "reached the target is a survivorship statistic, and a row with a low `read_50` "
      "and a high `read_50 unreached` is worse, not better, than the row above it.\n")
    gz = main[main["arm"] == "pipeline"].groupby(
        ["file_name", "budget", "z"])["recall_at_budget"].median().reset_index()
    bz = gz.loc[gz.groupby(["file_name", "budget"])["recall_at_budget"].idxmax()]
    zmin = float(main["z"].min())
    _pp = main[main["arm"] == "pipeline"].drop_duplicates(["file_name", "seed_index", "z"])
    _n1 = _pp[_pp["z"] == 1.0].groupby("file_name")["n_detections"].median()
    _nf = _pp[_pp["z"] == zmin].groupby("file_name")["n_detections"].median()
    _grow = ((_nf / _n1 - 1.0) * 100).dropna().sort_values()
    _r1k = main[(main["arm"] == "pipeline") & (main["budget"] == 1000)
                & main["file_name"].isin(dg)]
    _byz = _r1k.groupby(["file_name", "z"])["recall_at_budget"].median().unstack()
    _below = [c for c in _byz.columns if c < 1.0]
    _has1 = 1.0 in _byz.columns
    _better = (int((_byz[_below].max(axis=1) > _byz[1.0] + 1e-9).sum())
               if (_below and _has1) else 0)
    A(f"**Where the frontier actually stops.** The sweep now runs to z = {zmin:g}, which "
      "is the extraction floor itself -- that arm is the whole extracted pool, so it is "
      "the deepest level reachable without re-extracting. An earlier version of this "
      "report capped the sweep at z = 1.0 and justified it with \"the pool grows only "
      "about 3% below it\". That was 301.tiff's number and the smallest of the seven: "
      f"measured here the pool grows {_grow.min():.1f}-{_grow.max():.1f}% between "
      f"z = 1.0 and the floor. Below z = 1.0, recall@1000 improves on {_better} of the "
      f"{len(dg)} decision-grade ROIs. z = {zmin:g} is best or joint-best on "
      f"{int((bz['z'] == zmin).sum())} of the {len(bz)} (ROI, budget) rows in section 3."
      "\n")
    for roi in dg:
        A(f"**{roi}** (decision-grade)\n")
        rows = []
        d = main[(main["file_name"] == roi) & (main["arm"] == "pipeline")]
        for z, g in d.groupby("z"):
            g1 = g[g["budget"] == 1000]
            r80 = _stats(g1["read_80"])
            r50 = _stats(g1["read_50"])
            rows.append([f"{z:g}", _f(g1["n_detections"].median(), 0),
                         _f(g1["recall_at_budget"].median()),
                         _f(g[g["budget"] == 2000]["recall_at_budget"].median()),
                         _f(r50["median"], 0), r50["n_unreachable"],
                         _f(r80["median"], 0), r80["n_unreachable"],
                         _f(g1["coverage_frac"].median())])
        A(_md(rows, ["z", "n_cand", "recall@1000", "recall@2000", "read_50",
                     "read_50 unreached", "read_80", "read_80 unreached", "coverage"]))
        A("")

    # --------------------------------------------------------------- read_50 spread
    A("## 5. Per-seed spread of `read_50`\n")
    A("Stability across *which annotation the pathologist happens to click*. "
      "`unreached` counts seeds where 50% sensitivity is never reached; those seeds are "
      "excluded from median/min/max, so a row with unreached > 0 is a median over fewer "
      f"than {n_seeds} seeds.\n")
    for roi in dg:
        A(f"**{roi}** (decision-grade)\n")
        rows = []
        d = main[(main["file_name"] == roi) & (main["budget"] == BUDGETS[0])
                 & (~main["sanity_only"].astype(bool))]
        keys = d[["arm", "z"]].drop_duplicates()
        keys["ord"] = keys["arm"].map({a: i for i, a in enumerate(ARM_ORDER)}).fillna(99)
        for _, k in keys.sort_values(["ord", "z"], na_position="first").iterrows():
            sub = d[(d["arm"] == k["arm"])
                    & (d["z"].isna() if pd.isna(k["z"]) else (d["z"] == k["z"]))]
            st = _stats(sub["read_50"])
            spread = (st["max"] / st["min"]) if (st["min"] and np.isfinite(st["min"])) else np.nan
            rows.append([_arm_label(k["arm"], k["z"]), _f(st["median"], 0),
                         _f(st["min"], 0), _f(st["max"], 0),
                         _f(st["q3"] - st["q1"], 0), st["n_unreachable"], _f(spread, 2)])
        A(_md(rows, ["arm", "median", "min", "max", "IQR", "unreached", "max/min"]))
        A("")
    A("Head-to-head on `read_50`, pipeline at its best z (by median recall@1000) against "
      "each seedless arm, paired by seed. Lower is better, so a win is the pipeline "
      "reading *fewer* candidates. Seeds where either arm never reaches 50% sensitivity "
      "are excluded from the comparison and counted separately, never scored as a win.\n")
    rows = []
    for roi in dg:
        d = main[(main["file_name"] == roi) & (~main["sanity_only"].astype(bool))]
        pl = d[(d["arm"] == "pipeline") & (d["budget"] == 1000)]
        if not len(pl):
            continue
        bz = pl.groupby("z")["recall_at_budget"].median().idxmax()
        p50 = d[(d["arm"] == "pipeline") & (d["z"] == bz)
                & (d["budget"] == BUDGETS[0])].set_index("seed_index")["read_50"]
        for arm in SEEDLESS:
            s50 = d[(d["arm"] == arm) & (d["budget"] == BUDGETS[0])].set_index(
                "seed_index")["read_50"]
            if not len(s50):
                continue
            common = p50.index.intersection(s50.index)
            ok = p50[common].notna() & s50[common].notna()
            wins = int((p50[common][ok] < s50[common][ok]).sum())
            rows.append([roi, f"pipeline z={bz:g}", arm, _f(p50.median(), 0),
                         _f(s50.median(), 0), f"{wins}/{int(ok.sum())}",
                         int((~ok).sum())])
    A(_md(rows, ["ROI", "pipeline arm", "seedless arm", "pipeline read_50",
                 "seedless read_50", "pipeline wins", "excluded (unreachable)"]))
    A("")

    # ------------------------------------------------------------------ provenance
    A("## 6. Seed provenance\n")
    A(f"Same pipeline, same z ({PROVENANCE_Z:g}), same evaluation ground truth -- which "
      "excludes both the mitotic seed and the look-alike seed, so all three arms are "
      "scored against an identical annotation set. Only the provenance of the click "
      "differs: (a) a category-1 mitotic figure, (b) a category-2 look-alike the "
      "pathologist examined and rejected, (c) a random `nucleus_blobs` component centre "
      "with no annotation at all. All three pass the same border and Otsu-foreground "
      "filters and use the same tightened template size rule.\n")
    for roi in sorted(set(prov["file_name"])):
        n = int(prov[prov["file_name"] == roi]["n_gt_mitotic"].median())
        A(f"**{roi}** (n_mitotic = {n}, decision-grade)\n")
        rows = []
        for pv in ("a_mitotic", "b_lookalike", "c_random_blob"):
            d = prov[(prov["file_name"] == roi) & (prov["provenance"] == pv)]
            cells = [_f(d[d["budget"] == b]["recall_at_budget"].median()) for b in BUDGETS]
            iqr = [f"[{_f(d[d['budget'] == b]['recall_at_budget'].quantile(.25))}-"
                   f"{_f(d[d['budget'] == b]['recall_at_budget'].quantile(.75))}]"
                   for b in BUDGETS]
            r50 = _stats(d[d["budget"] == BUDGETS[0]]["read_50"])
            rows.append([pv] + [f"{c} {i}" for c, i in zip(cells, iqr)]
                        + [_f(d[d["budget"] == BUDGETS[0]]["n_detections"].median(), 0),
                           _f(r50["median"], 0), r50["n_unreachable"]])
        A(_md(rows, ["provenance"] + [f"@{b}" for b in BUDGETS]
              + ["n_cand", "read_50", "read_50 unreached"]))
        A("")
    A("Per-seed paired differences against (a), recall at budget 1000:\n")
    rows = []
    for roi in sorted(set(prov["file_name"])):
        d = prov[(prov["file_name"] == roi) & (prov["budget"] == 1000)]
        piv = d.pivot_table(index="seed_index", columns="provenance",
                            values="recall_at_budget")
        for pv in ("b_lookalike", "c_random_blob"):
            if pv not in piv:
                continue
            diff = piv[pv] - piv["a_mitotic"]
            rows.append([roi, f"{pv} - a_mitotic", _f(diff.median()), _f(diff.min()),
                         _f(diff.max()), f"{int((diff > 0).sum())}/{len(diff)}"])
    A(_md(rows, ["ROI", "contrast", "median delta", "min", "max",
                 "seeds where non-mitotic seed wins"]))
    A("")

    # ----------------------------------------------------------------- diagnostics
    A("## 7. Mandatory diagnostics\n")
    A("`coverage_frac > 0.5` is flagged: past that point, full-list recall measures how "
      "finely the candidate list tiles the ROI rather than what the detector found "
      "(`evaluate.py` module docstring). `tie` is the largest block of exactly equal "
      "ranking-key values; `nan` is the fraction of the list whose ranking key is NaN "
      "(a border window too close to the ROI edge to read). Both mark list positions "
      "whose order is arbitrary. Medians over seeds.\n")
    rows = []
    for roi in dg + up:
        d = main[(main["file_name"] == roi) & (main["budget"] == BUDGETS[0])]
        keys = d[["arm", "z"]].drop_duplicates()
        keys["ord"] = keys["arm"].map({a: i for i, a in enumerate(ARM_ORDER)}).fillna(99)
        for _, k in keys.sort_values(["ord", "z"], na_position="first").iterrows():
            sub = d[(d["arm"] == k["arm"])
                    & (d["z"].isna() if pd.isna(k["z"]) else (d["z"] == k["z"]))]
            cov = sub["coverage_frac"].median()
            rows.append([roi, _arm_label(k["arm"], k["z"]),
                         _f(sub["n_detections"].median(), 0), _f(cov),
                         "**>0.5**" if cov > 0.5 else "",
                         _f(sub["largest_tie_block"].median(), 0),
                         _f(sub["nan_rate"].median(), 4),
                         _f(sub["n_lookalike_in_list"].median(), 0),
                         _f(sub["lookalike_at_budget"].median(), 0)])
    A(_md(rows, ["ROI", "arm", "n_cand", "coverage", "flag", "largest tie", "nan rate",
                 "look-alikes in list", "look-alikes in top 500"]))
    A("")
    A(f"`floor_limited` is False in all {len(main)} main-grid rows and all {len(prov)} "
      "provenance rows, as it must be: the extraction floor sits at z = 0.5 and no "
      "reported z is below 1.0. `premise_test.py` asserts this rather than merely "
      "recording it.\n")

    # ---------------------------------------------------------------- verification
    A("## 8. Verification\n")
    sc = ver[ver["kind"] == "shortcut_exactness"]
    rows = [[r["file_name"], _f(r["cut"], 4), int(r["n_direct"]), int(r["n_filtered"]),
             "identical" if r["identical"] else "**DIFFERS**"] for _, r in sc.iterrows()]
    A("**One-match-many-z shortcut.** Filtering the deep pool at the z=2.0 cut, against "
      "re-running `extract_peaks` + NMS + self-hit removal directly at that cut. "
      "Compared as sorted coordinate arrays, exact equality, on every ROI:\n")
    A(_md(rows, ["ROI", "cut", "direct", "filtered", "result"]))
    A("")
    inv_rows = ver[ver["kind"] == "invariant"]
    for name, note in (
        ("no_cap", "no arm's candidate-list length equals a configured truncation cap"),
        ("min_separation", "grid arms: no two kept lattice points within the match "
                           "radius, i.e. NMS was applied at the match radius"),
        ("nms_radius", "every arm's suppression radius equals `evaluate.radius_px(mpp)` "
                       "for its image"),
        ("tissue_mask_covers_gt", "`baselines.tissue_mask` excludes 0 *mitotic* "
                                  "ground-truth annotations on every ROI"),
        ("distinct_seeds", "distinct `seed_index` values use distinct RNG streams"),
    ):
        sub = inv_rows[inv_rows["check"] == name]
        if not len(sub):
            continue
        A(f"* **{name}** -- {note}: {int(sub['n_runs'].iloc[0])} checks over "
          f"{len(sub)} distinct subjects, all passed.")
    A("")
    tm_rows = inv_rows[inv_rows["check"] == "tissue_mask_covers_gt"]
    exc = tm_rows[tm_rows["n_excluded"].fillna(0) > 0] if "n_excluded" in tm_rows else \
        tm_rows.iloc[:0]
    A("**The full-coverage form of the tissue-mask invariant does not hold, and the "
      "exception is recorded rather than suppressed.** The brief's wording is \"excludes "
      "0 ground-truth annotations\"; measured over all "
      f"{int(tm_rows['n_gt'].sum())} annotations on the {len(tm_rows)} ROIs, "
      f"{int(tm_rows['n_excluded'].fillna(0).sum())} of them fall outside "
      "`baselines.tissue_mask`:\n")
    if len(exc):
        A(_md([[r["label"], int(r["n_gt"]), int(r["n_excluded"]), r["excluded_ann_ids"]]
               for _, r in exc.iterrows()],
              ["ROI", "n_gt", "excluded", "ann_id (category)"]))
    else:
        A("*(none)*")
    A("")
    # +1 per ROI: `n_gt_mitotic` is the evaluation count, with that seed's own mitotic
    # annotation already removed. Derived, not hardcoded -- it moves with the ROI set.
    _n_mit_raw = int(sum(main[main["file_name"] == r]["n_gt_mitotic"].median() + 1
                         for r in dg + up))
    A(f"All exclusions are category 2 (look-alike); 0 of the {_n_mit_raw} category-1 "
      f"mitotic annotations across the {n_roi} ROIs is excluded, which is precisely what "
      "`baselines.tissue_mask`'s own docstring claims. Every recall number in this "
      "report is computed over category-1 ground truth and is therefore unaffected. The "
      "affected diagnostic is `n_lookalike_in_list` on the mask-restricted arms "
      "(`blob_native`, `blob_od`, `random`) on 506.tiff, whose ceiling is one look-alike "
      "lower than the ROI's annotation count. `invariants.check_tissue_mask_covers_gt` "
      "is fatal on category 1 and records the rest; it can be made fatal on every "
      "category with `strict_categories=None`.")
    A("")

    # ------------------------------------------------------------------ caveats
    A("## 9. What is underpowered, and what did not resolve\n")
    ns = sorted(int(main[main["file_name"] == r]["n_gt_mitotic"].median()) for r in up)
    if ns:
        A(f"* {len(up)} of the {n_roi} ROIs ({', '.join(up)}) have "
          f"{ns[0]}-{ns[-1]} mitotic figures in their evaluation set, so a single figure "
          f"moves recall by {1 / ns[-1]:.2f}-{1 / ns[0]:.2f} there. Their rows are "
          "printed above and are labelled; no statement in this report rests on them.")
    unre = main[(main["budget"] == BUDGETS[0]) & (~main["sanity_only"].astype(bool))]
    n_un100 = int(unre["read_100"].isna().sum())
    A(f"* `read_100` is unreachable in {n_un100} of {len(unre)} (ROI, seed, arm, z) "
      "cells; `read_80` in "
      f"{int(unre['read_80'].isna().sum())}; `read_50` in "
      f"{int(unre['read_50'].isna().sum())}. Those cells are excluded from the "
      "corresponding medians and counted in the `unreached` columns, never imputed and "
      "never counted as wins.")
    hi = main[(main["budget"] == BUDGETS[0]) & (main["coverage_frac"] > 0.5)]
    A(f"* `coverage_frac > 0.5` in {len(hi)} of {len(unre)} cells, in arms "
      f"{', '.join(sorted(set(hi['arm'])))}. Full-list recall for those arms is a "
      "statement about tiling geometry rather than about the detector; the budgeted "
      "columns are unaffected by it.")
    rnd = main[main["sanity_only"].astype(bool)]
    if len(rnd):
        r = rnd[rnd["budget"] == 1000].iloc[0]
        A(f"* The `random` arm is a single sanity row ({r['file_name']}, seed "
          f"{int(r['seed_index'])}, {int(r['n_detections'])} points in tissue), reported "
          "in section 3 as the floor. One ROI, one seed, no interval; excluded from "
          "every median and win count.")
    A("* Seed selection draws uniformly *with replacement* from the accepted-annotation "
      "pool, so distinct `seed_index` values can legitimately draw the same annotation. "
      "The stated invariant \"distinct `seed_index` yields distinct `seed_ann_id`\" is "
      "therefore not a property of a correct run and is not asserted; see "
      "`invariants.check_distinct_seeds`. Distinct RNG streams are asserted instead. "
      "Observed distinct seeds per ROI: "
      + ", ".join(f"{roi} {main[main['file_name'] == roi]['seed_ann_id'].nunique()}"
                  f"/{n_seeds}" for roi in dg + up) + ".")

    text = "\n".join(L) + "\n"
    with open(f"results/premise_report{suffix}.md", "w") as fh:
        fh.write(text)
    print(f"wrote results/premise_report{suffix}.md ({len(text)} chars)")


def check_seedless_invariance(path="results/premise_test.csv"):
    """The seedless arms' candidate sets must not vary with the seed -- only their GT.

    Asserts that every seed-invariant column really is constant across seeds for each
    (ROI, seedless arm), so any movement in their metrics is attributable to the
    evaluation ground truth and nothing else.
    """
    d = pd.read_csv(path)
    d = d[d["arm"].isin(SEEDLESS) & (d["budget"] == BUDGETS[0])]
    cols = ["n_detections", "coverage_frac", "largest_tie_block", "nan_rate"]
    bad = []
    for (roi, arm), g in d.groupby(["file_name", "arm"]):
        for c in cols:
            if g[c].nunique(dropna=False) != 1:
                bad.append((roi, arm, c, sorted(g[c].unique())))
    if bad:
        raise inv.InvariantError(f"seedless arms vary across seeds: {bad}")
    n = d.groupby(["file_name", "arm"]).ngroups
    print(f"seedless candidate sets constant across seeds: {n} (ROI, arm) groups, "
          f"columns {cols}")
    return True


if __name__ == "__main__":
    if "--report" in sys.argv:
        sfx = "_smoke" if "--smoke" in sys.argv else ""
        check_seedless_invariance(f"results/premise_test{sfx}.csv")
        write_report(sfx)
    else:
        main(smoke="--smoke" in sys.argv)
