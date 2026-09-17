# AnnotateDx click-to-verify pipeline: stage-by-stage breakdown

Written 2026-09-12, ahead of packaging a shippable pipeline for GitHub. Scope: the
**find-and-suppress / click-to-verify pipeline** in `midog_utils/`, not the upstream
MIDOG++ RetinaNet training/inference code (`main.py`, `training.py`, `inference.py`,
`evaluation.py`, `configs/all.yaml`) — that is a separate, unrelated project sharing this
repo.

**Read this before the rest of the doc.** There is no single `pipeline.py`. "Production"
is a convention — a set of decisions in `DECISIONS.md` and `D8_TEMPLATE_ANCHOR.md` —
implemented with varying fidelity across notebooks. The most current, most complete
embodiment is `production_seed_precision_at_k/production_seed_precision_at_k_chromatin_half_pix_fix.ipynb`
(cell 1 config, `run_roi` for the per-ROI loop), and that notebook is what this doc
walks through. But **its literal constants and the decision record disagree in places**
— those disagreements are called out inline, at the stage they affect, not saved for the
end. If you package this for GitHub without resolving them, you are shipping the
disagreement, not the decision.

**Status (2026-09-16):** superseded as a description of production by `midog_utils/production.py`
(see `production_pipeline/PRODUCTION_PIPELINE_WALKTHROUGH.md`). All six punch-list items are
resolved. Code removed in the `c066829` cleanup is marked below. Line refs date from 2026-09-12.

Every stage below is given as three columns of truth:
- **Decided** — what `DECISIONS.md` / `D8_TEMPLATE_ANCHOR.md` says should happen
- **Library default** (`midog_utils/find_and_suppress.py`'s `FSConfig`, or the relevant
  module) — what you get if you construct the config with no overrides
- **Notebook value** — what the production notebook's cell 1 actually sets

---

## Stage 0 — Channel and template-match method

| | Decided | Library default (`FSConfig`) | Notebook value |
|---|---|---|---|
| Channel | `hematoxylin_od` (D3) | `"gray_inverted"` (`find_and_suppress.py:49`; field removed in `c066829`) | `CHANNEL = 'hematoxylin_od'` |
| TM method | `cv2.TM_CCOEFF` (D1) | `cv2.TM_CCOEFF_NORMED` (`find_and_suppress.py:93`; `cv2.TM_CCOEFF` since `c066829`) | `METHOD = cv2.TM_CCOEFF` |

`hematoxylin_od` (`midog_utils/chromatin.py:102-110`) = `rgb2hed(rgb/255.0)[:,:,0]` —
**unclipped** optical density from skimage's colour deconvolution, hematoxylin channel
only. D3's whole point: never rescale/clip this to 0–255 (that operation was
`to_hematoxylin`, removed in `c066829`), because clipping saturates exactly
the densest-chromatin pixels at 255 and destroys the contrast signal `TM_CCOEFF` needs.

**Gap (resolved in `c066829`: the default is now `cv2.TM_CCOEFF`):** D1 was applied only where a call site explicitly passed `METHOD = cv2.TM_CCOEFF`.
`FSConfig.tm_method`'s dataclass default was `TM_CCOEFF_NORMED`
(`DECISIONS_UNVERIFIED.md:105-109` — D1 itself was moved out of `DECISIONS.md` on
2026-09-12 as *never independently re-checked*, not reversed). **Anyone who constructs
`FSConfig()` with defaults and calls `find_and_suppress()` directly gets the
non-D1 method.** If this ships as a library call rather than notebook constants, either
change the dataclass default or make the caller pass `tm_method=cv2.TM_CCOEFF` explicitly
and say so in the API.

---

## Stage 1 — Seed draw and gate (`midog_utils/seed_selection.py`)

1. **`agreement_pool`** (`:48-66`) — prefer unanimous-mitotic annotations; fall back to
   the ≥2/3-mitotic contested tier (flagged).
2. **`border_filter`** (`:69-82`) — keep annotations at least
   `BORDER = CFG.patch_size // 2 = 36` px from the ROI edge.
3. **RNG draw** (notebook `run_roi`): `rng = np.random.default_rng([SEED_INDEX, image_id])`,
   drawing without replacement via `draw_seed_with_retry`.
4. **`tighten_box_otsu`** — the gate (`:104-263`). Otsu-threshold a 51 px window
   (`otsu_window = tm.BASE_SIZE = 51`) around the click; take the connected component
   under the click's own rounded pixel (no tolerance). Reject (redraw) unless:
   `min_area=50`, `max_area_frac=0.85`, `min_solidity=0.5` (all `seed_selection.py:107-109`
   defaults), and the click falls inside that component's own bbox.
5. **`tightened_template_box`** (`:431-475`) — the D8-current seed constructor:
   ```python
   base_size = _odd(max(y1 - y0, x1 - x0), minimum=minimum)
   half = otsu_window // 2
   ix, iy = int(round(cx)), int(round(cy))
   center_x = ix - half + (x0 + x1 - 1) / 2.0
   center_y = iy - half + (y0 + y1 - 1) / 2.0
   ```
   The `- 1` is the D8 half-pixel correction: a skimage bbox is half-open
   (`x0..x1-1`), and `read_padded_patch` rounds its centre argument to a pixel index, so
   the *pixel*-centre formula is `(x0+x1-1)/2`, not the *area*-centre `(x0+x1)/2` used
   before 2026-09-10. Verified on 896/993 accepted seeds (14 ROIs): the corrected formula
   fully contains the Otsu component 100.0% of the time (0/896 crops) vs. 60.5% for the
   old formula and 28.2% for anchoring on the raw click. D8_TEMPLATE_ANCHOR.md is explicit
   that this is a **mechanism** claim only — it does not claim a precision/recall
   improvement, and the one outcome-level read that exists is single-seed and mixed.
6. Ground truth (`gt_eval`) stays keyed to the **raw click** `(cx, cy)` throughout; only
   the *template* is built from the gated/recentred box. This split (click for scoring,
   gated-and-recentred box for the template) is D8's point.

**Gap — D8 is "partly in code."** `tightened_template_box` computes the corrected
centre, but nothing routes self-hit removal to that returned centre except `build_seed`
— and nothing in the production notebook calls `build_seed`; the notebook's own
`run_roi`/`draw_seed_with_retry` inlines the equivalent logic. Confirm before shipping
that the inlined path in the notebook actually matches `tightened_template_box`'s formula
line for line (it does, per direct comparison against `seed_selection.py:431-475`, but
this is exactly the kind of copy-paste drift the repo's own audits flag elsewhere).

**Gap — a stale contradiction exists in the repo.**
`pipeline_debug_visuals/seed_refinement_variants.ipynb` (commit `24d2429`, dated
2026-09-11, one day *after* `D8_TEMPLATE_ANCHOR.md` went current) states in its own
markdown cell that the click-centred, pre-D8 formula is "the production path today."
No `DECISIONS.md` entry re-reverses back to click-centred after 2026-09-10 — this reads
as a stale notebook cell, not a real re-reversal, but **flag it to your postdoc and fix
or delete that cell before shipping**, since a reader who opens that notebook first will
walk away with the wrong answer.

**The "2×2 gate × centre" grid**, for context on why gating and recentring are treated
as separate, independently-justified axes (`pipeline_debug_visuals/seed_refinement_variants.ipynb`):

| | centred on **the click** | centred on **the component bbox** |
|---|---|---|
| **ungated** | `ungated_click` | `ungated_recentred` |
| **gated** | `gated_click` (pre-D8) | `gated_recentred` — **D8, current** |

Two ROIs isolate the two axes: on 245.tiff the gate is inert (all four variants pick the
same annotation; only centring differs, a 4.3px shift); on 403.tiff the gate changes
*which* annotation is drawn at all (ungated accepts a click sitting on background and
borrows an unrelated component's size).

---

## Stage 2 — Template augmentation bank

| | Decided/library default | Notebook value |
|---|---|---|
| `scales` | `(1.0,)` (`FSConfig` default) | `(1.0,)` |
| `n_angles` | `1` (`FSConfig` default) | `1` |
| `flips` | `(False,)` (`FSConfig` default) | `(False,)` |

**Single-scale, single-angle, no flip — this is the current default, not a simplification
made for one notebook.** `FSConfig`'s own comment (`find_and_suppress.py:59-65`)
documents this as a deliberate "simplest configuration for the next run," explicitly
*not* a claim that it outperforms the previous 12-angle × 2-flip default — a 3-way
augmentation-footprint comparison (`n_angles=1` vs. `4×2` vs. `12×2`,
`midog_utils_full/experiment.py:565-569`) exists and, on the one seed measured so far, actually ranked the
4-angle/2-flip config best and the current no-augmentation default worst. That comparison
was never promoted to a decision. `base_size = tm.BASE_SIZE = 51`,
`patch_size = tm.PATCH_SIZE = 73`.

**Verified for this doc:** the D9 `max_peaks` ablation (below) was measured under the
*identical* `FSConfig(...)` construction line as production — same single-aug config —
so D9's result transfers cleanly on this axis; it is not confounded by an augmentation
mismatch between what was measured and what ships.

---

## Stage 3 — Template matching → fused response map (`template_match.fused_response`, `:183-246`)

The ROI is padded by `PAD = max((t.shape[0]-1)//2 for t in templates)` with
`cv2.BORDER_REPLICATE` before matching (`find_and_suppress.py:135-136`), so that no
part of the actual ROI is unreachable after cropping back. For each augmented template,
`cv2.matchTemplate(img, tmpl, cv2.TM_CCOEFF)` is run, and the per-pixel **element-wise
max across the whole augmentation bank** is kept, in a shared centre-coordinate frame
(each template's raw correlation output is shifted by its own `((t.shape[0]-1)//2,
(t.shape[1]-1)//2)` offset before the max — the module's own docstring flags this offset
as easy to get wrong, worth a 10–25 px shift if mishandled). `scale_normalize=False`
(`FSConfig` default) — raw, unbounded `TM_CCOEFF` correlation values are kept, not
z-scored; this is what makes them only comparable *within* one image/run, never across
runs or images.

With a single-augmentation bank (Stage 2), the "max across augmentations" step is a
no-op in the production notebook today — it is architecture for the general case, not
something currently exercised.

**The `valid` mask** (constructed inside `fused_response`, not by any caller): a pixel is
`True` wherever *some* augmentation's `matchTemplate` output actually reached it once
shifted into the shared centre frame. This is a **border-reachability artifact of
correlation geometry**, not a tissue/foreground test — there is no brightness or tissue
check anywhere in this function (that idea was explicitly rejected, D2). Because the ROI
is pre-padded exactly enough to make every real ROI pixel reachable, `valid` ends up
trivially all-`True` over the ROI in the production/ablation notebooks; the notebooks
assert this (`assert bool(valid.all())`).

---

## Stage 4 — Deep-floor peak extraction

```python
med, mad = tm.robust_stats(fused, valid)              # template_match.py:66-83
cut = med + DEEP_FLOOR_Z * mad                          # DEEP_FLOOR_Z = -1.5
centers, scores = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE, cut, MAX_PEAKS)
```

| Constant | Value | Notes |
|---|---|---|
| `robust_stats` | median, `1.4826·MAD` of `fused` on a stride-8 sample of `valid` pixels | replaces an older `sample > -1.5` heuristic |
| `PEAK_MIN_DISTANCE` | 7 | grey-dilation pre-thinning: `k = 2·7+1 = 15`; local maxima kept where `fused >= dilated` |
| `DEEP_FLOOR_Z` | −1.5 | score cutoff = `median − 1.5·MAD` — deliberately permissive; the label "deep floor" means almost everything above noise survives to this pool |
| `MAX_PEAKS` (pre-NMS cap) | **`2,000,000`** in the production notebook cell 1 (never binds); **`100`** only in `threshold_maxpeaks_ablation/*` (D9) | **see gap below** |

Peak ordering is a global lexicographic sort, `order = np.lexsort((ys, xs, -scores))[:max_peaks]`
(`template_match.py:264-274`) — score descending, then x, then y — chosen specifically so
"filter the pool at a threshold" and "re-extract at that threshold" agree exactly,
including tie order.

**Gap — D9's decided value is not wired into the notebook this doc walks through.**
*(Resolved: `production.MAX_PEAKS = 100`.)*
`DECISIONS.md` (D9, 2026-09-12) adopts `max_peaks = 100` before NMS as production,
citing zero measured precision/recall cost through K=30 on 14 ROIs and a ~99.7% NMS
runtime reduction. That value lives only in
`threshold_maxpeaks_ablation/max_peaks_100_variant.ipynb` and its sibling ablation
notebooks. The notebook this doc otherwise treats as "current production"
(`production_seed_precision_at_k_chromatin_half_pix_fix.ipynb`) still has
`MAX_PEAKS = 2_000_000` in cell 1. **Since you've decided to ship with `max_peaks=100`,
this is the one constant you must actually change at the point you assemble the
shippable pipeline — it is not already there.** D9 itself flags its own result as
untested past K=30 and single-seed (short of D4/D5's multi-seed bar) — worth restating
in whatever doc/README ships with the code, not just carrying forward silently as settled.

---

## Stage 5 — NMS (distance-based) + self-hit removal

**NMS** (`midog_utils/nms.py:22-49`, `nms_by_distance`): greedy, score-descending,
KD-tree radius query.
```python
order = np.argsort(-scores, kind="stable")
for idx in order:
    if suppressed[idx]: continue
    keep.append(idx)
    suppressed[neighbours[idx]] = True   # KDTree radius query
```
Radius = `evaluate.radius_px(mpp, 7.5)` (D7: `NMS_RADIUS_UM = evaluate.MIDOG_RADIUS_UM =
7.5` µm) — **computed per-image from that image's microns-per-pixel**, 29.6–33.1 px
across the scanners in this dataset, never a hardcoded pixel constant. D7 is enforced,
not just documented: `invariants.check_nms_radius` (`invariants.py:177-187`) raises
`InvariantError` unless the radius used exactly equals `evaluate.radius_px(mpp)`. A
stray hardcoded `5.0` µm radius that propagated into one earlier sweep notebook is
explicitly named in `DECISIONS.md` as "a bug, not a configuration choice" — at 7.5 µm vs.
5.0 µm the deep pool loses 6 mitoses (of ~5,000+) across 4 ROIs, and the deployable
z cutoff moves from 1.308 to 0.438. If you fork any config for experimentation, this
invariant is your tripwire — do not disable it.

**Self-hit removal** (notebook `suppress()`; replaced in production on 2026-09-17 by blanking the seed's own template footprint before correlation, see `SELF_HIT_MASKING_PLAN.md` and `DECISIONS.md` D11):
```python
keep = nms_by_distance(centers, scores, radius)
c, s = centers[keep], scores[keep]
ok = np.hypot(c[:,0]-ref_xy[0], c[:,1]-ref_xy[1]) > SELF_HIT_RADIUS
c, s = c[ok], s[ok]
```
`SELF_HIT_RADIUS = 5.0` px — deliberately **tighter** than the ~30 px match radius,
because the two closest real MIDOG++ annotations anywhere in the dataset are 26.2 px
apart (403.tiff); using the full match radius to drop the seed's self-hit could delete a
legitimate neighbouring detection. `ref_xy` is the D8-recentred template centre
(`tpl_xy`), not the raw click.

---

## Stage 6 — Chromatin-family features on the NMS survivor pool

Computed on `pool` (post-NMS, post-self-hit), all on the same `hematoxylin_od` channel:

| Feature | Formula | Constants |
|---|---|---|
| `od{31,51,81}` | `chromatin.chromatin_density(hem, x, y, window=w)` = mean of the darkest `frac` fraction of pixels in a `w`×`w` window (`chromatin.py:113-137`) | `frac=0.10` (`DEFAULT_FRAC`, `chromatin.py:99`) for all three windows |
| `od_ctx` | same function, wider/softer window | `window=121, frac=0.50` |
| `od_contrast` | `od51 − od_ctx` | candidate's own darkest-10%-of-51px density minus the darkest half of a 121px neighbourhood; a **within-image difference**, constructed so the OD scale's lack of cross-run calibration cancels out (D5, 2026-09-08 amendment) |
| `od_falloff` | `od31 − od81` | |
| `mask_od_mean` | mean intensity of the **largest** Otsu component under the candidate | gate-free — no accept/reject, unlike Stage 1's gate |

`chromatin.py`'s own module docstring states in bold that this module **"is NOT the
production ranker"** — `chromatin.rerank`, which would swap `score` for an OD axis as the
sort key, was removed in `c066829`; production sorts by `od` inline when `rank_key="chromatin_od"`.

---

## Stage 7 — Ranking (`midog_utils/compare.py`, `Arm` + `_rank`)

*`compare.py` was removed in `c066829`; production ranks with the same stable sort inside
`production.run_production_pipeline`.*

Each candidate axis is declared as an `Arm`:
```python
Arm(name, candidates=lambda: pool, rank_key=key, caps=(MAX_PEAKS,),
    nms_radius=nms_radius, coverage_key=fn)
```
and sorted by:
```python
out = df.sort_values(key, ascending=False, na_position="last", kind="mergesort")
```
`kind="mergesort"` is **stable** — this is the tie-break mechanism. Since `pool` arrives
already `tm_score`-descending (the post-NMS invariant), any tie on a chromatin axis
breaks in `tm_score` order. **`tm_score` is the implicit secondary sort key for every
arm, including `chromatin_od`.**

This module already contains almost exactly the "knob" you're asking for — the axis
registry:
```python
AXES = {
    'tm_score':     'score',
    'chromatin_od': 'od51',
    'od31':         'od31',
    'od_falloff':   'od_falloff',
    'mask_od_mean': 'mask_od_mean',
    'od_contrast':  'od_contrast',
}
```
— but **redeclared inline, per notebook**, not a shared library-level registry. There's
no config file or CLI flag anywhere that selects `rank_key`; every notebook that wants to
compare axes copy-pastes this dict and the feature-computation code that produces its
columns. *(Resolved: `production.AXES` is the shared registry; `run_pipeline.py --rank-key` selects from it.)*

**Gap — `check_no_cap` is structurally vacuous at `max_peaks=100`.**
*(Resolved: `check_no_cap` was removed in `c066829`; `info["max_peaks_binding"]` records the cap instead.)*
`invariants.check_no_cap(len(det), arm.caps)` compares `len(det)` — `det = arm.candidates()`,
the **post-NMS** pool the arm actually ranks — against `caps=(MAX_PEAKS,)`, which is the
**pre-NMS** extraction cap. At `MAX_PEAKS=100`, the cap binds pre-NMS on every ROI (14/14
in the ablation), but NMS + self-hit removal always strip at least the seed's own
detection, so the post-NMS count lands at 89–99, never exactly 100 — `check_no_cap` can
never fire at this configuration. This was independently found and documented in
`Research Logs/2026-09-12-chromatin-od-ranker-variant-audit.md` ("Finding 6"): the
check "passes by construction," which the notebook does acknowledge, but a printed
"70/70 checks passed" reads stronger than what it means. This is a real property of how
the check is wired (comparing the wrong stage's length against a pre-NMS cap), not a
crash anyone has hit — but it means **the safety net that's supposed to catch "the pool
silently got truncated" cannot catch that at your chosen `max_peaks` value.** If you want
this invariant to mean something once shipped, it needs a cap on the post-NMS length
too, or an explicit acknowledgment in the shipped docs that it's inert here by
construction.

---

## Stage 8 — The top-K budget step (answers "what happens at the top-10 budget")

`compare.evaluate_arms` (`compare.py:133-231`; removed in `c066829`, production's equivalent is
`run_pipeline.summarize`) did the truncation:
```python
for k in budgets:                      # BUDGETS = (10, 20, 30, 50) in the notebook
    delivered = min(int(k), len(tp_cum))
    tp_at = int(tp_cum[delivered - 1]) if delivered else 0
    ...
    recall_at_budget = tp_at / n_mit
```
**There is no separate top-K selection stage.** The final K-length list *is*
`ranked.head(K)` of the arm's full stable re-sort from Stage 7 — no re-scoring, no
dedup/diversity pass, nothing else. "Top-10" = `rank_key`-descending (tm_score-tie-broken
for every arm), first 10 rows of the post-NMS, post-self-hit pool.

**How "true positive" is decided at that budget** — `evaluate.greedy_match`
(`evaluate.py:57-82`), the function every precision/recall number in this repo ultimately
rests on:
```python
def greedy_match(det_xy, gt_xy, radius):
    tree = KDTree(gt_xy)
    neighbours = tree.query_radius(det_xy, r=radius)
    for d_idx, cands in enumerate(neighbours):        # detections, best-first, in rank order
        free = [g for g in cands if gt_to_det[g] == -1]
        if not free: continue
        g = free[int(np.argmin(distance to each free candidate))]
        det_to_gt[d_idx] = g; gt_to_det[g] = d_idx
    return det_to_gt, gt_to_det
```
- **Radius**: the same 7.5 µm / `evaluate.radius_px(mpp)` value as NMS (D7) — two
  separate `radius_px(mpp, ...)` calls (`MATCH_RADIUS_UM`, `NMS_RADIUS_UM`) held equal by
  decision and by `invariants.check_nms_radius`, not the same Python variable.
- **One-to-one**: yes — each GT can be claimed by at most one detection
  (`gt_to_det[g] == -1` gate), each detection claims at most one GT.
- **Order matters and is deliberate**: detections are processed **best-first, in rank
  order** — the exact order Stage 7 produced. A higher-ranked detection claims a
  contested ground-truth object before a lower-ranked one ever sees it. The module
  docstring states this is deliberate: a single greedy pass (vs. recomputing an optimal
  assignment at every threshold) means a top-K prefix's matches never depend on
  lower-ranked detections — recomputing per-threshold could reassign an earlier
  detection and produce a non-monotone curve. `evaluate.optimal_assignment_disagreement`,
  a purely diagnostic Hungarian-algorithm cross-check, was removed in `c066829`; it never
  drove any reported precision@K.

Precision at budget = `tp_at_budget / budget_delivered`. At `max_peaks=100`, note D9's
own caveat: the capped pool can under-deliver even a K=100 ask (89–99 candidates
survive NMS), so `budget_delivered` can be less than the nominal budget for larger K —
irrelevant at K=10 (no ROI in the ablation dropped below 89 survivors) but worth keeping
in mind if the shipped budgets list is ever extended upward.

---

## On the tm_score / chromatin_od knob

What you're asking for — `tm_score` default, `chromatin_od` selectable — is
**consistent with D5** as it stands: D5 says `chromatin_od` may not become the *default*
without the outstanding 5-seed, `images/extra_valid`-measured, ROI-clustered paired
comparison beating `tm_score` (not yet done — at the last checked configuration the two
axes were statistically indistinguishable: paired Δrecall@250 = +0.032, 95% CI
[−0.047, +0.112]). Offering it as an **opt-in second axis**, not the default, doesn't
trigger that bar.

Two things worth documenting alongside the knob, not burying:
1. **The production seed (`SEED_INDEX=0`) is not a neutral draw for this comparison.**
   A 5-seed sweep (`chromatin_od_ranker_seed_robustness.ipynb`, 2026-09-12) found seed 0
   ranks best-of-5 for `tm_score`-vs-`chromatin_od` at every budget checked, and the
   all-5-seed mean effect is 2–6× smaller than the seed-0-only value. If your postdoc's
   users mostly experience seed 0 (the hardcoded default everywhere), whatever apparent
   edge `chromatin_od` shows in your own testing is inflated relative to what a random
   click will show.
2. **`chromatin_od`'s tie-break is still `tm_score`** (Stage 7) — switching the knob
   changes the primary sort key, not the fallback, which is worth one line in whatever
   docstring/README documents the flag so nobody's surprised that two candidates with
   identical `od51` don't tie-break arbitrarily.

Implementation-wise, per Stage 7: the natural place for this knob is to promote the
inline `AXES` dict into a shared `midog_utils` module (it doesn't exist as shared code
today — every notebook redeclares it) and expose `rank_key` as the parameter the knob
sets, rather than inventing a new mechanism — `Arm.rank_key` already does exactly this.
*(Done: `production.AXES` and `run_production_pipeline(rank_key=...)`.)*

---

## Punch list before this ships

*All resolved; see "Punch-list resolution" in `production_pipeline/PRODUCTION_PIPELINE_WALKTHROUGH.md`.*

1. Set `MAX_PEAKS = 100` at the point you assemble the shippable entrypoint — it is
   `2,000,000` in the notebook this doc calls "production" today (Stage 4).
2. Decide whether `FSConfig.tm_method`'s default changes to `cv2.TM_CCOEFF`, or the
   shippable entrypoint passes it explicitly and documents why (Stage 0).
3. Fix or delete the stale "click-centred is production" cell in
   `pipeline_debug_visuals/seed_refinement_variants.ipynb`, or add a dated correction —
   it contradicts `D8_TEMPLATE_ANCHOR.md` and will mislead the next reader (Stage 1).
4. Promote the per-notebook `AXES` dict into one shared registry in `midog_utils`, and
   wire the tm_score/chromatin_od knob through it rather than a new mechanism (Stage 7,
   knob section).
5. Decide what to do about `check_no_cap`'s vacuity at `max_peaks=100` — leave it
   documented as inert-by-construction here, or add a genuine post-NMS cap check
   (Stage 7).
6. Carry D9's own caveats (single-seed, untested past K=30) and D5's seed-0-bias caveat
   into whatever user-facing docs ship with the knob — both are currently only in
   `DECISIONS.md`/`Research Logs/`, not anywhere a downstream user of the shipped code
   would see them.
