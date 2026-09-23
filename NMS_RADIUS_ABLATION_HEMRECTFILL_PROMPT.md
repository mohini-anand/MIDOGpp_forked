# Task: does widening NMS beyond the match radius fix same-object duplicate candidates in `hem_rect_fill`'s top-10/20?

## Goal

Compare 4 arms on the same 147 clicks (49 ROIs x 3), template = **`hem_rect_fill` only** (not
`default_51`), ranked by `chromatin_od`:

| arm | NMS radius |
|---|---|
| `nms_1.00x` (baseline) | `1.00 * ev.radius_px(mpp)` -- current production |
| `nms_1.25x` | `1.25 * ev.radius_px(mpp)` |
| `nms_1.50x` | `1.50 * ev.radius_px(mpp)` |
| `nms_2.00x` | `2.00 * ev.radius_px(mpp)` |

The **match radius (scoring)** stays `ev.radius_px(mpp)` -- unchanged, in every arm. Only the NMS
suppression radius moves. Question: does precision@10/20 improve, hold, or get worse as NMS
suppresses more aggressively than the radius it's currently pinned to -- and at what cost to
genuinely distinct nearby mitoses?

Repo: `/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked`, branch `find-and-suppress-midog`.

## Why this experiment, in plain words

A prior session found that two ranked candidates on `353.tiff` (a `human_neuroendocrine_tumor`
ROI, outside this notebook's own 49-ROI set but the same mechanism applies everywhere) were 39.05
px apart -- both within the ~32.98 px match radius of the *same* ground-truth mitotic figure
(ann_id 18303: 16.76 px from one candidate, 26.08 px from the other), but current NMS (radius =
match radius) doesn't suppress either, because they're farther apart from *each other* than the
NMS radius. The higher-ranked one (rank 3, by `chromatin_od`) claims the annotation; the
lower-ranked one (rank 8) reads as `non_human_findings` even though it's sitting on the same real
object, not a genuine miss.

`DECISIONS.md`'s **D7** (2026-09-09) set NMS radius = match radius deliberately, and names this
exact failure mode as the reason: *"If NMS keeps two survivors closer together than the match
radius, both can sit inside one ground-truth object's radius -- one is credited a true positive
and the other a false positive on the same object."* But D7's own ablation (`f5_nms_radius_ablation.py`)
only ever compared two *shared* radius values (5.0 µm vs 7.5 µm, i.e. NMS radius and match radius
moved together, both <= 7.5 µm) -- **NMS radius has never been set *larger* than match radius in
this codebase.** `invariants.check_nms_radius` currently enforces strict equality
(`abs(used - want) > tol` raises), so this decoupled configuration is new. This notebook is that
measurement -- restricted to `hem_rect_fill`, not `default_51`, per the user's explicit choice, and
restricted to what the top-10/20 looks like, not deep-pool recall (D7's own concern; out of scope
here on purpose).

## Read first

- `midog_utils/find_and_suppress.py` (current, post-D11 pipeline), `production.py`,
  `template_match.py`, `seed_selection.py`, `chromatin.py`, `evaluate.py`, `invariants.py`.
- `DECISIONS.md`'s **D7** entry in full, and its cross-referenced **D11** application commit
  (`git show 9e22eb9`) -- `production.SELF_HIT_RADIUS` no longer exists; D11 replaced the
  post-NMS self-hit-disc filter with `template_match.blank_seed_square`, which blanks the seed's
  own template footprint out of a *copy* of the search channel before correlation, so the seed's
  own peak never forms. There is no separate self-hit step to configure in this notebook.
- Code cells of `production_hematoxylin_only/hem_rect_fill_vs_default_chromatin_od_49roi_3seed.ipynb`
  ("the hem_rect_fill notebook") -- `cut_templates` is copied from it unchanged (given verbatim below).
- **This notebook's baseline (`nms_1.00x`) has no pre-D11 recorded reference to reproduce.**
  Every notebook in the `hem_rect_fill_vs_default_*` / `*_by_domain` / `human_neuroendocrine_followup`
  family predates D11 (their recorded `per_run.csv` files have a populated `n_self_hits` column, a
  D11-deleted info key) -- do not try to bit-match against them. The correctness anchor here is
  `production.run_production_pipeline` itself (see "Parity gate" below), not historical CSVs.

If a step can't run as written, or is wrong, stop and report. Don't work around a halt.

## Inputs (read-only)

- `HRF = results/precision_at_k_49roi_3seed_chromatin_hemrectfill`. Read
  `{HRF}_per_run.csv` with `pd.read_csv(..., float_precision='round_trip')`.
- `databases/MIDOG++.json`, `images/extra_valid/*.tiff`, `images/extra_valid/testing_set/*.tiff`.

## Before running

- `git status --short midog_utils/` prints nothing. Don't modify `midog_utils`.
- `{HRF}_per_run.csv` passes: 294 rows; `condition` is exactly `default_51`, `hem_rect_fill`, 147
  rows each; 49 `file_name`s, each with `seed_index` 0, 1, 2 under both conditions.

## Notebook: `production_hematoxylin_only/hem_rect_fill_nms_radius_ablation_chromatin_od_49roi_3seed.ipynb`

nbconvert runs it from its own directory, so every path starts with `../`.

### Cell 0 (markdown), written before the run

State the goal, the 4 arms, the "why" section above (condensed), and what's identical across arms:
window position/size, search image, peak settings, ranking, scoring. Only NMS radius moves; match
radius does not. State these reading rules, fixed now:

- **Bootstrap CI excludes 0, delta > 0 (vs `nms_1.00x`):** widening NMS to this multiple raises precision at this K.
- **Bootstrap CI excludes 0, delta < 0:** it lowers precision.
- **CI includes 0:** no difference detected. Quote the CI as the bound.

Caveats to state:
- `chromatin_od` ranking only; `hem_rect_fill` template only (not `default_51`); the 353.tiff
  pattern that motivated this occurred in **1 case out of 132 (click, arm) pairs** in an unrelated,
  smaller follow-up sample -- this notebook's 147 clicks give more power, but the underlying event
  may still be rare enough that a null result here is not strong evidence of "no effect," just "no
  effect detected at this n."
- Deep-pool recall and the deployable score cutoff (`z_max`) are **out of scope** here by the
  user's explicit choice -- D7's own concern, not re-litigated. K=30 is reported descriptively only
  (not part of the Holm family) since the user only cares about K=10/20.
- Clicks are inherited from the recorded `hemrectfill` run's `default_51` rows (same as every
  notebook in this family) -- not redrawn.

### Cell 1: config

- **Imports:** `gc, hashlib, itertools, math, sys, time`; `from fractions import Fraction`; `cv2`,
  `numpy as np`, `pandas as pd`, `from scipy import stats as sps`, `matplotlib.pyplot as plt`;
  `sys.path.insert(0, '..')`; from `midog_utils`: `channels as ch`, `chromatin as cm`,
  `dataset as ds`, `evaluate as ev`, `invariants as inv`, `production as prod`,
  `seed_selection as ss`, `template_match as tm`; `from midog_utils.nms import nms_by_distance`;
  `from midog_utils.find_and_suppress import FSConfig, find_and_suppress`.
- **Config checks:** `REFERENCE_CONFIG = dict(CHANNEL='hematoxylin_od', TM_METHOD=cv2.TM_CCOEFF, PEAK_MIN_DISTANCE=7, DEEP_FLOOR_Z=-1.5, MAX_PEAKS=100, OD_WINDOW=51)`
  -- **no `SELF_HIT_RADIUS` key**, D11 removed it from `production.py`. Build `CONFIG_DRIFT` the
  same way the rest of this notebook family does and assert it all matches. Also assert
  `not hasattr(prod, 'SELF_HIT_RADIUS')` (documents that this notebook knowingly targets the
  post-D11 pipeline; if this ever fires, D11 was reverted and the whole approach below needs
  re-checking). Assert `ev.MIDOG_RADIUS_UM == 7.5`, `tm.BASE_SIZE == 51`, `tm.PATCH_SIZE == 73`.
- **Constants:**
  - `SMOKE_N_ROIS = None`, `N_ROIS = 49 if SMOKE_N_ROIS is None else SMOKE_N_ROIS`.
  - `RADIUS_MULTIPLIERS = {'nms_1.00x': 1.00, 'nms_1.25x': 1.25, 'nms_1.50x': 1.50, 'nms_2.00x': 2.00}`.
    `ARMS = tuple(RADIUS_MULTIPLIERS)`. `BASELINE_ARM = 'nms_1.00x'`. Every delta is `arm - BASELINE_ARM`.
  - `BUDGETS = (10, 20, 30)`; `PRIMARY_BUDGETS = (10, 20)` -- K=30 computed but excluded from the
    Holm family and from the D5-style bar, reported descriptively only.
  - `SEED_INDICES = (0, 1, 2)`, `N_BOOT = 10_000`, `BOOT_SEED = <pick today's date as YYYYMMDD, document it>`.
  - `HRF = '../results/precision_at_k_49roi_3seed_chromatin_hemrectfill'`.
  - `SUBSET_DIRS = {'original_14': '../images/extra_valid', 'testing_35': '../images/extra_valid/testing_set'}`.
- **Outputs:** `OUT_STEM = '../results/precision_at_k_49roi_3seed_chromatin_hemrectfill_nms_radius_ablation' + ('' if SMOKE_N_ROIS is None else '_smoke')`;
  `OUT = {name: f'{OUT_STEM}_{name}.csv' for name in ('clicks', 'per_run', 'top30', 'per_roi', 'delta_per_roi', 'delta_stats', 'verification', 'duplicates', 'collateral')}`.
  `FIG_PATHS = {...}` -- your choice of what to visualize (at minimum: a forest-style plot of the
  per-K, per-arm delta with its CI, following the dataviz conventions the rest of this repo's
  notebooks use -- load the `dataviz` skill if it's available to you before writing chart code).
  Assert no `OUT` value collides with `{HRF}_*` or the plain `hemrectfill_by_domain_*` files.

### Cell 2: clicks

- `REC = pd.read_csv(f'{HRF}_per_run.csv', float_precision='round_trip')`.
- `CLICKS` = the `REC` rows with `condition == 'default_51'`, in CSV row order (these are the
  click coordinates -- both arms of the *original* notebook share them; you are not using its
  `hem_rect_fill` per-run rows, since this notebook recomputes everything at the D11 pipeline).
- `ROI_ORDER` = the first `N_ROIS` distinct `file_name`s in that order.
- `images, annotations = ds.load_annotations('../databases/MIDOG++.json')`. Verify every
  `seed_ann_id` against the JSON exactly as every prior notebook in this family does (category
  MITOTIC, cx/cy match).
- Assert `len(CLICKS) == 3 * N_ROIS`.

### Cell 3: functions

Style: one-line `def`s and calls; docstrings per "Writing rules" below.

1. **`cut_templates(hem, cx, cy)`** -- copied verbatim from the hem_rect_fill notebook:
   ```python
   def cut_templates(hem, cx, cy):
       patch = tm.read_padded_patch(hem, cx, cy, tm.PATCH_SIZE)
       assert patch is not None, f'({cx}, {cy}): the 73 px patch is not readable'
       templates, _ = tm.build_augmentations(patch, tm.BASE_SIZE, (1.0,), 1, (False,))
       T = templates[0]
       W = tm.read_padded_patch(hem, cx, cy, tm.BASE_SIZE)
       assert np.array_equal(W, T)
       rect = ss.tighten_box_otsu(W)
       assert rect is not None, f'({cx}, {cy}): tighten_box_otsu refused the window'
       y0, y1, x0, x1 = rect
       outside = np.ones((51, 51), bool)
       outside[y0:y1, x0:x1] = False
       fill = np.float32(T[outside].astype(np.float64).mean())
       F = T.copy()
       F[outside] = fill
       return dict(T=T, F=F, rect=rect, fill=fill)
   ```

2. **`run_hem_rect_fill(hem, shape, cx, cy, mpp, nms_radius)`** -- reimplements
   `find_and_suppress`'s *current* (post-D11) internals for a single 51 px template, substituting
   the `hem_rect_fill` template for the one `find_and_suppress` would cut itself, and taking
   `nms_radius` as an independent parameter (not `ev.radius_px(mpp)`). Match `find_and_suppress`'s
   own order of operations exactly -- re-read it before writing this, don't work from memory of it:
   1. `cut = cut_templates(hem, cx, cy)`; `template = cut['F']`.
   2. **Blank a COPY of `hem` at the seed's own footprint**, matching D11:
      `search, n_blanked = tm.blank_seed_square(hem, cx, cy, tm.BASE_SIZE)`. `hem` itself is
      untouched (needed unblanked for `chromatin_od` scoring below -- this mirrors
      `production.py`'s own documented behavior: *"production.py still hands find_and_suppress
      the unblanked hem, which it reuses for chromatin_od ranking"*).
   3. `pad = tm.BASE_SIZE // 2` (= 25, matching `find_and_suppress`'s
      `max((t.shape[0]-1)//2 for t in templates)` for one unscaled 51 px template).
      `h, w = shape`. `padded = cv2.copyMakeBorder(search, pad, pad, pad, pad, borderType=cv2.BORDER_REPLICATE)`.
   4. `fused_p, best_p, valid_p = tm.fused_response(padded, [template], False, method=prod.TM_METHOD)`.
      `fused, valid = fused_p[pad:pad+h, pad:pad+w], valid_p[pad:pad+h, pad:pad+w]`. Assert `valid.all()`.
   5. `med, mad = tm.robust_stats(fused, valid)`. `threshold = med + prod.DEEP_FLOOR_Z * mad`.
   6. `centers, scores = tm.extract_peaks(fused, valid, prod.PEAK_MIN_DISTANCE, threshold, prod.MAX_PEAKS)`.
      `n_peaks = len(centers)`.
   7. **NMS at the swept radius**, not `ev.radius_px(mpp)`:
      `keep = nms_by_distance(centers, scores, nms_radius)`; `centers, scores = centers[keep], scores[keep]`.
   8. No self-hit step -- D11 already prevented the seed's own peak from forming in step 2.
   9. `det = pd.DataFrame({'rank': np.arange(len(centers)), 'cx': centers[:, 0], 'cy': centers[:, 1], 'score': scores})`.
   10. **`chromatin_od` off the unblanked, padded `hem`** (not `search`/`padded` from step 3):
       `hem_padded = cv2.copyMakeBorder(hem, pad, pad, pad, pad, borderType=cv2.BORDER_REPLICATE)`;
       `shifted = det.assign(cx=det['cx'] + pad, cy=det['cy'] + pad)`;
       `det = det.assign(od=cm.score_detections(shifted, hem_padded, window=prod.OD_WINDOW)['od'].to_numpy())`.
   11. `det = det.sort_values('od', ascending=False, na_position='last', kind='mergesort').reset_index(drop=True)`;
       `det = det.assign(rank=np.arange(len(det)))`.
   12. Returns `(det, {'n_peaks': n_peaks, 'n_detections': len(det), 'n_blanked_px': n_blanked, 'nms_radius': nms_radius})`.

3. **`score_top_k(detections, gt_eval, match_radius)`** -- copied from the hem_rect_fill notebook
   (unchanged; scoring logic is untouched by D11 or by this experiment).

### Cell 4: parity gate (smoke this first, on 2 ROIs, before the full run)

**Before trusting `run_hem_rect_fill` on the `hem_rect_fill` template, validate it against real
production on `default_51`** -- the one case where a direct comparison is possible, since
`production.run_production_pipeline` can't take a custom template.

1. Build a `seed_selection.Seed` (read its dataclass fields in `seed_selection.py` first; don't
   guess them) with `template_xy=(cx, cy)`, `base_size=51`, for one click.
2. Run `prod.run_production_pipeline(rgb, seed, mpp, rank_key='chromatin_od', max_peaks=prod.MAX_PEAKS)`.
3. Run your own `run_hem_rect_fill`-style pipeline but with `template = cut['T']` (the plain
   `default_51` crop, not `F`) and `nms_radius = ev.radius_px(mpp)` (i.e. the `nms_1.00x` setting).
4. Assert these match **exactly**: `n_peaks`, `n_detections`, the full ranked `cx`/`cy`/`score`/`od`
   sequences (score at float32, od at `rtol=1e-12`). If they don't, the reimplementation has a bug
   -- stop and report the first divergence; do not proceed to the full run.
5. Only once this passes on 2 ROIs x 3 seeds does the `hem_rect_fill` arm (swapping `T` for `F`)
   get to be trusted as "the same mechanism, different template pixels."

Record this gate's result in `OUT['verification']` and in Cell 10's markdown.

### Cell 5: run loop

For each `file_name` in `ROI_ORDER`, each of its 3 clicks, each of the 4 `ARMS`:

1. Load `rgb`, `mpp`, `hem = ch.to_channel(rgb, 'hematoxylin_od')` once per ROI (not per arm).
2. `gt = ds.image_annotations(annotations, file_name)`; `match_radius = ev.radius_px(mpp)`
   (**fixed**, independent of the arm); `gt_eval = gt[gt['ann_id'] != seed_ann_id]`.
3. For each arm: `nms_radius = RADIUS_MULTIPLIERS[arm] * match_radius`; run `run_hem_rect_fill(...)`;
   score with `score_top_k(det, gt_eval, match_radius)` (match radius the same in every arm).
4. Write `per_run` rows (subset, domain, file_name, image_id, seed_index, arm, seed_ann_id,
   click_cx/cy, mpp, match_radius_px, nms_radius_px, n_gt_mitotic, `n_peaks`, `n_detections`,
   `n_blanked_px`, `od_nan`, `od_descending`, `tp_at_10/20/30`) and `top30` rows (rank, cx, cy,
   score, od, bucket, matched_ann_id, matched_category, file_name, seed_index, arm) -- same
   pattern as every prior notebook.
5. Persist incrementally, free memory, print one progress line per ROI (elapsed minutes, TP@K per
   arm for that ROI's clicks).

**Runtime:** ~2x the original hem_rect_fill notebook's 294 runs (this is 4 arms x 147 clicks = 588
template searches). That notebook measured 1.1-1.7 s/run; budget 30-60 min. Smoke-test with
`SMOKE_N_ROIS = 2` first; run the full pass in the background.

### Cell 6: verification

- **Halt:** `n_detections >= 30` every run; `od_nan == 0`; `od_descending`; every (file_name,
  seed_index) has exactly 4 rows, one per arm, equal `seed_ann_id`; `n_detections` is
  non-increasing as the radius multiplier increases, within each click (wider NMS can only remove
  survivors, never add them -- if this fails, something is wrong with the reimplementation).
- **Report:** `n_peaks == prod.MAX_PEAKS` rate per arm; `n_detections` distribution per arm
  (confirms wider NMS is actually thinning the pool); the Cell 4 parity gate result.

### Cell 7: duplicate-resolution diagnostic

For **each arm independently** (not just baseline), detect the "hub-and-spoke" pattern that
motivated this notebook: a top-30 candidate `C` with `bucket == 'non_human_findings'` such that
some mitotic annotation `m` (excluding the click's own `seed_ann_id`) is within `match_radius` of
`C`, **and** `m.ann_id` is already claimed (`matched_ann_id == m.ann_id`) by a *better-ranked*
candidate in that same click's list. Report, per arm: the count of such instances (total, and
separately within rank < 10 and rank < 20), and the full list (file_name, seed_index, C's rank,
claimer's rank, ann_id, distance C-to-m). This count should trend toward 0 as the radius multiplier
increases -- confirm whether it reaches 0 at `nms_2.00x` (the theoretical guarantee point: two
points each within `match_radius` of a shared third point are at most `2 * match_radius` apart, so
NMS at that radius must suppress one of them). Write `OUT['duplicates']`.

### Cell 8: collateral-damage diagnostic

For each click and each of `nms_1.25x`/`nms_1.50x`/`nms_2.00x`, compare the set of `matched_ann_id`
values with `bucket == 'human_correct_label'` in that click's **top-10** (and separately top-20)
against the same set at `nms_1.00x`. Flag any `ann_id` matched in the baseline's top-K but absent
from the wider-radius arm's top-K. For each flagged case, check whether that `ann_id` still appears
matched *anywhere* in the wider arm's full top-30 list (still found, just pushed down in rank -- a
ranking effect) versus not found at all (its peak was actually suppressed by the wider NMS --
true collateral damage). Report counts of each kind, per arm, per K. Write `OUT['collateral']`.

### Cell 9: precision and statistics

Same machinery as every notebook in this family -- copy these functions verbatim (they're already
validated in this session; re-derive nothing):

```python
def sign_flip_null_counts(abs_deltas):
    d = [int(v) for v in abs_deltas]
    assert all(v > 0 for v in d) and len(d) <= 62
    total = sum(d)
    counts = np.zeros(2 * total + 1, dtype=np.int64)
    counts[total] = 1
    for v in d:
        shifted = np.zeros_like(counts)
        shifted[v:] += counts[:len(counts) - v]
        shifted[:len(counts) - v] += counts[v:]
        counts = shifted
    return counts

def exact_sign_flip_p(deltas):
    d = np.asarray(deltas).astype(np.int64)
    nz = np.abs(d[d != 0])
    if len(nz) == 0:
        return 1.0, 0
    counts = sign_flip_null_counts(nz)
    total = int(nz.sum())
    observed = abs(int(d.sum()))
    sums = np.arange(-total, total + 1)
    return float(Fraction(int(counts[np.abs(sums) >= observed].sum()), 2 ** len(nz))), len(nz)

def holm_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p, kind='stable')
    adjusted = np.maximum.accumulate((len(p) - np.arange(len(p))) * p[order]).clip(max=1.0)
    out = np.empty_like(adjusted)
    out[order] = adjusted
    return out

def boot_ci_exact(int_values, idx, denominator):
    sums = np.asarray(int_values, dtype=np.int64)[idx].sum(axis=1)
    lo, hi = np.percentile(sums / denominator, [2.5, 97.5])
    return float(lo), float(hi)

def pair_stats(d, n_roi, k, boot_idx):
    pooled_delta_precision = d.sum() / (3 * n_roi * k)
    wins, losses, ties = int((d > 0).sum()), int((d < 0).sum()), int((d == 0).sum())
    exact_p, n_nonzero = exact_sign_flip_p(d)
    ci_low, ci_high = boot_ci_exact(d, boot_idx, 3 * n_roi * k)
    x = d / (3 * k)
    if n_roi > 1 and x.std(ddof=1) > 0:
        half = sps.t.ppf(0.975, n_roi - 1) * x.std(ddof=1) / math.sqrt(n_roi)
        t_ci_low, t_ci_high = float(x.mean() - half), float(x.mean() + half)
    else:
        t_ci_low = t_ci_high = pooled_delta_precision
    ci_excludes_0 = bool(ci_low > 0 or ci_high < 0)
    n_same_direction = wins if pooled_delta_precision > 0 else (losses if pooled_delta_precision < 0 else 0)
    majority = bool(n_same_direction > n_roi / 2)
    return dict(n_roi=n_roi, delta_tp_sum=int(d.sum()), pooled_delta_precision=pooled_delta_precision, wins=wins, losses=losses, ties=ties, n_nonzero=n_nonzero, exact_p=exact_p, ci_low=ci_low, ci_high=ci_high, t_ci_low=t_ci_low, t_ci_high=t_ci_high, ci_excludes_0=ci_excludes_0, majority=majority, bar_met=bool(ci_excludes_0 and majority))
```

Validate `exact_sign_flip_p` against brute-force enumeration on random small vectors first (same
pattern as every prior notebook), before trusting it on real data.

- `per_roi`: `tp_at_K` summed over 3 clicks, per (file_name, arm). `precision_at_K = tp_at_K / (3*K)`.
- Pooled precision per arm and K: `sum(tp_at_K) / (3 * N_ROIS * K)`, with a cluster-bootstrap CI
  (unpaired) using one `BOOT_IDX = np.random.default_rng([BOOT_SEED, 0]).integers(N_ROIS, size=(N_BOOT, N_ROIS))`
  shared across every CI in this notebook.
- `delta_per_roi`: for each of `nms_1.25x`/`nms_1.50x`/`nms_2.00x`, `delta_tp_sum = arm_tp - baseline_tp`
  per ROI, per K (K=10, 20, 30 all computed; 30 flagged `primary=False` in the output).
- `delta_stats`: `pair_stats` per (arm, K). **Holm family = the 3 arms x 2 primary K (10, 20) = 6
  tests, corrected together.** K=30 rows get their own `exact_p` but no Holm column populated (or
  populated over its own, separate 3-test family across arms -- your choice, just don't merge it
  into the primary family). A domain-style `bar_met` (CI excludes 0 and majority of the 49 ROIs
  agree) is fine to report per row alongside Holm, exactly like `hem_rect_fill_vs_default_chromatin_od_by_domain.ipynb`.

### Cell 10: figure(s)

At minimum, one chart showing the per-K (10, 20; 30 visually distinguished as secondary), per-arm
pooled delta with its bootstrap CI vs the `nms_1.00x` baseline -- readers should be able to see at
a glance whether any arm's CI clears 0, and in which direction. If you have access to the
`dataviz` skill, load it before writing chart code and use its validated palette conventions
(diverging blue/red for CI-exclusion direction, neutral gray for "CI includes 0") -- this repo's
other notebooks in this family already follow that convention; match it rather than inventing a new one.

### Cell 11 (markdown), written after the full run

Every number must be printed by a cell above. Report, in this order:
1. The Cell 4 parity gate result (this validates everything after it).
2. Report checks (n_peaks binding rate, n_detections distribution per arm).
3. Pooled precision per arm and K, with CIs.
4. The delta_stats table (K=10/20 primary with Holm; K=30 descriptive).
5. The duplicate-resolution counts per arm (Cell 7) -- does `nms_2.00x` actually reach 0?
6. The collateral-damage findings per arm (Cell 8) -- any true suppression losses, not just rank shuffling?
7. A plain-language verdict: at which multiplier (if any) does precision@10/20 clear its CI in the
   positive direction, survive Holm, *and* show no material collateral damage? If none, say so
   plainly -- don't round a null result up to "promising."

Read the result with Cell 0's rules only. Make no mechanism claims beyond what's measured here.

## Running

- **Interpreter:** `/Users/mohinianand/anaconda3/bin/python3` and
  `/Users/mohinianand/anaconda3/bin/jupyter` only -- the bare `python3` on this machine can't
  `import cv2`.
- **Command:** `/Users/mohinianand/anaconda3/bin/jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=10800 production_hematoxylin_only/hem_rect_fill_nms_radius_ablation_chromatin_od_49roi_3seed.ipynb`.
- **Smoke run first:** `SMOKE_N_ROIS = 2`. Every cell, including the Cell 4 parity gate, must pass
  before the full run.
- **Full run:** delete the smoke CSVs, set `SMOKE_N_ROIS = None`, run in the background (30-60+ min expected).
- **Finish:** write Cell 11. Final state: `SMOKE_N_ROIS = None`, execution counts contiguous from
  1, no error outputs.

## Writing rules

- Function definitions and calls go on one line, never wrapped.
- Docstrings: opening `"""` alone on its line, content indented one level deeper, one-line summary,
  one `name (type): description` line per parameter, one `Returns type: description` line.

## Don't

- Commit.
- Modify `midog_utils` or any `{HRF}_*` / `*_by_domain_*` file.
- Use `default_51` as the tested template -- `hem_rect_fill` only, per the user's explicit choice.
- Change the match radius, peak settings, or ranking key.
- Redraw clicks.
- Call `inv.check_nms_radius` on any arm except `nms_1.00x` -- it will (correctly) raise on the
  other three, since they deliberately violate the invariant it checks. Don't work around this by
  loosening the invariant itself; just don't call it on those arms, and say why in Cell 0.
- Treat a null result (no arm clears CI + Holm + no collateral damage) as a reason to keep
  searching for a multiplier that "works" -- report what's measured.

## Report back

- The Cell 4 parity gate: passed or not, and the first divergence if not.
- Report checks (n_peaks binding, n_detections distribution per arm).
- Pooled precision per arm and K.
- The delta_stats table (K=10/20 primary + Holm; K=30 descriptive).
- Duplicate-resolution counts per arm -- does `nms_2.00x` reach 0?
- Collateral-damage findings per arm -- true suppression losses vs rank shuffling.
- The plain-language verdict from Cell 11.
- Runtime, files written.
- Any deviation from this prompt, and why.
