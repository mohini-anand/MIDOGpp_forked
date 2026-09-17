# Audit of `bbox_refinement_three_way_chromatin_od.ipynb`: every TP count, table and p-value reproduces from pixels with 0 divergences, and the TP/FP scoring survives every alternative convention tried; the tie paragraph's replacement mechanism is refuted by the third pair it computed and did not cite, and the click caveat misses that production `gray_bbox` is also off its own click on 3 of 14 ROIs

**Scope.** `production_hematoxylin_only/bbox_refinement_three_way_chromatin_od.ipynb` (26 cells, 15 code),
its seven CSVs `results/precision_at_k_14roi_prodseed_chromatin_bbox3way_{per_roi,summary,by_domain,delta_per_roi,delta_stats,top30,verification}.csv`,
its three figures `production_hematoxylin_only/bbox3way_{pooled_precision,roi_delta_heatmap,win_counts}.png`, the
reference CSV it gates against (`results/precision_at_k_14roi_prodseed_chromatin_hembbox_raw.csv`), `DECISIONS.md` D1–D9,
`D8_TEMPLATE_ANCHOR.md`, `PRODUCTION_PIPELINE_CLEANUP.md`. Audit script: `bbox_refinement_three_way_chromatin_od_audit.py`
(repo root), tables `results/bbox_refinement_three_way_chromatin_od_audit_*.csv` (23 files).
**Everything below is re-derived from `databases/MIDOG++.json`, the TIFFs in `images/extra_valid/`, and the notebook's
persisted top-30 lists, with a matcher, seed draw and search pipeline written in the audit script. `midog_utils` is
never imported.**

## Conflict of interest

Signals I can see: the notebook and all its artifacts are **untracked** (`??`, never committed). The notebook was
written today (mtime 2026-09-16 17:40:31, execution metadata 21:35:15–21:40:29 UTC). The invoker says a different
Claude session wrote it. Every imported module is committed and clean (`git status --porcelain midog_utils/` is empty).
No git author can be attached to an untracked file. The work and this audit may come from the same project hand; only
my fresh context separates them. What limits the conflict: no number here is read from a notebook output table and
used as its own check. TP counts come from my matcher against the JSON. Clicks come from my own RNG walk and my own
implementation of the Otsu gate. Candidate lists come from my own `matchTemplate` → peaks → NMS → self-hit → od51
implementation. What this cannot cover is the design choices themselves: the ranker, the joint gate, K ≤ 30, and one
click per ROI. Those go to `premise-reviewer` and to a reader who did not make them.

---

## Part 0 — what reproduces

**The engineering is sound. The statistics the notebook publishes are the statistics the code in the working tree
computes, and an independent re-implementation from pixels reproduces all 42 top-30 lists coordinate-for-coordinate.**
Across the whole audit, 4,895 values were compared with 0 divergences (`..._comparison_counts.csv`).

| check | tier | values | divergences |
|---|---|---|---|
| Execution coherence: `execution_count` 1…15 contiguous from 1; 0 error outputs; 0 unrun cells; last cell run; 3 embedded PNGs. Not vendored (imports `midog_utils`, reads `../images/extra_valid`) | A | 15 cells | — |
| Composition: 42 = 14×3 per-ROI rows, 1,260 = 42×30 top-30 rows with ranks 0..29 in row order, 9/63/126/9/267 rows in summary/by_domain/delta_per_roi/delta_stats/verification. 7 domains × exactly 2 ROIs. 0 duplicate keys. 267/267 verification records passed. All 26,286 JSON boxes 50×50 | A | 17 checks | 0 |
| **TP re-match of every top-30 detection**: bucket and matched `ann_id` for 1,260 detections against `gt_eval` built from the JSON | A | 2,520 | **0** |
| Per-ROI table (TP@K, precision@K, `n_gt_mitotic`, domain) | A | 462 | 0 |
| Summary (tp_sum, delivered, pooled, pooled = mean-of-ROI identity, worst, n_at_worst, worst files, best) | A | 81 | 0 |
| By-domain table | A | 252 | 0 |
| Delta per ROI (integer and precision delta, base sizes, domain) | A | 504 | 0 |
| Delta stats (mean, W/L/T, n_nonzero, exact p by my own bitmask enumeration) | A | 54 | 0 |
| Gate 2 re-run against the reference CSV **as it is on disk now** (see T3-6) | A | 140 | 0 |
| Seed identity (right image, mitotic, unanimous); JSON vs TIFF dimensions | A | 70 | 0 |
| **Seed draw from pixels**: my RNG walk + my Otsu gate reproduce all 14 `seed_ann_id` and all 14 `n_retries` | B | 28 | 0 |
| Template geometry from pixels (`base_size`, `tpl_offset_px`, 3 conditions) | B | 84 | 0 |
| **Search from pixels**: top-30 coordinates identical on 42/42 runs (max od diff 4.3e-8), `n_detections` 42/42, TP@10/20/30 from full-list matching 126/126 | B | 210 | 0 |
| Prefix property: TP@K from matching the full list = matching only the top K (both matching conventions) | B | 252 | 0 |
| My `tm_score` and `chromatin_od` TP@10/20/30/50 for `gray_bbox`/`hem_bbox` vs the reference raw CSV | B | 224 | 0 |
| Figures, step 1a: every bar label (9), every heatmap annotation (126), every stack count and p label (36): the script emits the labels the notebook's format strings produce from recomputed values (`..._figure_labels_expected.csv`), and I read the decoded PNGs against that table by eye (the render is not machine-parsed); step 2: all three renders read | A + render (manual) | 171 | 0 |
| Prose: 72 numeric or ordering claims recomputed (`..._prose_claims.csv`) | A/B | 72 | 5 do not hold (Part 1) |

**Gates.** Provenance applied: modules clean, artifacts untracked, one oddity (T3-6). Composition applied. The
execution gate applied and passed. The vendored exemption does not apply. No claim was left without an artifact
underneath it: every printed number traces to a persisted CSV, the JSON, or the TIFFs. So nothing needed a
`cannot check`, and the pixel-only claims (annulus, self-hits, rounding on 201.tiff, the `n_retries` pattern) went
to Tier B.

### The TP/FP calculation, sub-question by sub-question (the invoker's priority)

| question | ruling | evidence |
|---|---|---|
| **Match radius and mpp conversion** | **Correct.** 7.5 µm / mpp, mpp from the TIFF `XResolution` tag with `ResolutionUnit` honoured (inch on 12 ROIs, cm on 529/548). Radii 29.61–33.14 px. | My own tag read reproduces every radius. Moving the radius ±1 px changes TP on **0/126** (run × K) cells. The nearest detection-to-GT distance to the boundary is 0.36 px, and it changes nothing. *Observed, immaterial:* `YResolution` ≠ `XResolution` on 7 ROIs (6 canine at 0.2482 vs 0.2491 µm/px, a 0.12 px radius gap; 094 at 0.002 px). `..._tiff_resolution_tags.csv` |
| **Greedy one-to-one order** | **Correct, and not costing anything here.** Rank order (the `chromatin_od` list), each detection claiming its nearest unclaimed GT, `d ≤ r`. | Maximum-cardinality matching (Hungarian) gives the same TP on 126/126 cells. Strict `d < r` gives the same on 126/126. No exact-distance tie between two free GTs on any of 1,260 detections. |
| **What `gt_eval` contains** | **Correct.** All annotations of the image of both categories, minus the click's `ann_id`. Look-alikes (15–123 per ROI) sit in the pool but are never credited. | Rebuilt from the JSON: `n_gt_mitotic` matches on 42/42. No other annotation lies within one radius of any click; the nearest other annotation is 62.0 px away (094.tiff). 159 of 1,260 top-30 detections are bucketed as look-alike hits, 412 as unannotated, 689 as TP. |
| **A detection matching a look-alike, or nearer a look-alike than a mitosis** | **Moot here.** **0 of 1,260** top-30 detections have both a mitosis and a look-alike within `r`. So mixed-pool greedy, mitotic-only greedy, mitotic-only maximum matching and "any mitosis within r" all give identical TP on **126/126** cells. | `..._top30_rematch_detections.csv`, `..._tp_by_scoring_rule_per_run.csv` (every `diff_*` column is 0 except `diff_contested_mitoses_dropped`, which is nonzero by construction). As a convention, the mixed pool can only turn a TP into an FP when a look-alike annotation is nearer. That is defensible and conservative, but mitotic-only matching followed by look-alike labelling of the unmatched detections is the cleaner definition of precision against mitotic GT. Worth adopting before a deeper list makes the two diverge. |
| **TP@K from the full list = TP@K from the top K alone** | **Correct** (cell 2's argument holds). | 252/252 from independently generated full lists (94–99 detections), under both matching conventions. Also 126/126 on the persisted top-30. |
| **Contested (2-of-3) mitoses** | **Sensible.** They are counted as mitotic GT, which is the dataset's consensus `category_id = 1`. The click pool is unanimous-only. They make up **165 of 689 (24 %)** top-30 TPs, evenly spread (56/238, 56/231, 53/220). | Sensitivity in T2-3: the headline survives, and strengthens, if contested hits are neither credited nor charged. It fails Holm only if they are charged as false positives. |
| **precision = TP/K** | **Correct.** Every list delivers ≥ 94 > 30. | `n_detections` min 97 / 96 / 94 per condition, reproduced from pixels. |
| **Off-by-one in top-K slicing** | **None.** `rank` is 0-based and contiguous; `hit[:k]` takes exactly K rows. | Composition gate: ranks 0..29 in row order on 42/42 runs. |
| **Self-hit / annulus** | **Correct.** Self-hit removal is at `seed.template_xy` (5 px). The annulus check is measured from the click. | From pixels: 1 self-hit on 42/42 runs, 0 full-list detections within `r` of the click on 42/42. |

### Items 2, 3 and 5 of the brief

- **The joint click gate, and "`default_51` adds nothing beyond `border_filter`": correct.** `border_filter` tests
  `rint(c) ∈ [36, dim−1−36]`, and `read_padded_patch(…, 73)` refuses when `round(c)±36` leaves the image. Both round
  half-to-even, so they are the same predicate. All 18 visited draws are readable at the click. The joint predicate
  therefore equals the reference notebook's `_check`, and my walk reproduces its clicks 14/14.
- **The conditions differ only in template size and centre: correct.** Reading `production.py`, only `seed.base_size`
  and `seed.template_xy` are read. The channel is `hematoxylin_od` and `OD_WINDOW = 51` for all three conditions;
  self-hit is at `template_xy`; ground truth is keyed to the click. Empirically, my pipeline varies only
  `(base, tx, ty)` and reproduces all 42 lists. What varies downstream (template pixels, which top-100 peaks survive)
  follows from those two inputs. The deep floor never binds (20,595–69,995 peaks above it, capped to 100 on 42/42).
- **Config: matches production.py, D1, D3, D7, D8 and D9** (`..._config_drift.csv`). The known three-way lag is
  still there: `FSConfig` still defaults to `TM_CCOEFF_NORMED`, `max_peaks = 250000`, `border_pad = False`,
  `deep_floor_z = None`, `nms_radius = None`, and production overrides all of them. **Declared deviations:** ranking
  by `chromatin_od` (D5 names `tm_score`; the notebook cites `PRODUCTION_PIPELINE_CLEANUP.md`), and a joint seed gate
  where `run_pipeline.py` gates on `gray_inverted` alone.

---

## Part 1 — findings

### Tier 1

#### T1-1 — The tie paragraph's replacement mechanism, "delta ties follow how different the two templates are", is refuted by the third pair, which the notebook computed and did not cite

Cell 24: *"That notebook called its `chromatin_od` K = 10 tie count (11/14) 'a structural property of the production
cap.' This run argues against that as a general explanation. At the same cap, `hem_bbox - default_51` has only 2 ties
at K = 20, and `gray_bbox - default_51` has 10. **Delta ties follow how different the two templates are, not the cap
alone.**"*

The paragraph has two halves, and they do not share a verdict.

**The negative half reproduces.** At the same cap, `hem − default` ties 2/14 and `gray − default` ties 10/14 at
K = 20. The cap alone cannot explain tie counts that differ that much between pairs.

**The positive half does not.** Here is "how different the two templates are", measured as mean |Δbase_size| over the
14 ROIs, next to all three pairs:

| pair | mean \|Δbase\| px | ties K=10 | ties K=20 | ties K=30 | top-10 shared coordinates (mean) | top-10 shared mitoses (mean) |
|---|---|---|---|---|---|---|
| `gray_bbox − default_51` | 11.9 | 10 | 10 | 6 | 1.57 | 6.14 |
| `hem_bbox − default_51` | **21.7** (most different) | 9 | 2 | 3 | 0.07 | 5.29 |
| `hem_bbox − gray_bbox` | **9.9** (most similar) | 11 | **2** | 4 | 0.14 | 5.21 |

The notebook compares two pairs, 21.7 px → 2 ties and 11.9 px → 10 ties, and on those two its claim holds. The third
pair decides it. `hem − gray` has the **smallest** size gap (9.9 px) and ties only **2** times at K = 20 and 4 at
K = 30. If ties followed template similarity, it should tie most. What the three pairs actually show is that any pair
involving `hem_bbox` rarely ties at K = 20/30, which is the headline effect (hem loses), not a similarity gradient.
At K = 10, the budget the quoted reference claim was about, tie counts are 10 / 9 / 11 and do not move with template
difference at all.

What does fit is the object-convergence account. The 1.57 shared coordinates for `gray − default` come entirely from
the two identical-template ROIs (201, 245: 10/10 each). On the other 12 ROIs the three pairs share 2, 1 and 2 top-10
coordinates *in total*, yet at K = 10 their lists claim 5.2–6.1 of the same mitoses on average. This is the round-3
audit's F1 mechanism (the chromatin axis converges on the same objects whatever template produced the candidates), and
it holds for all three pairs here. The paragraph also answers a K = 10 claim with K = 20 counts.

*Caveats on my numbers:* |Δbase| is one proxy for template difference. I did not test the distance between the two
templates' centres as a proxy. Per-pair means are over 14 ROIs, and the tie counts are 14-ROI integers.

**Fix:** keep the negative half; replace the positive sentence with the measured one: *"At K = 20/30, ties are rare
for any pair involving `hem_bbox` (2–4/14), including `hem − gray`, whose templates differ least in size, so ties do
not track template size difference. At K = 10 all three pairs tie on 9–11/14 while sharing ~0 candidate coordinates
and 5.2–6.1 of the same mitoses, consistent with the chromatin axis converging on the same objects."*
`..._tie_vs_template_difference.csv`.

### Tier 2

#### T2-1 — The click caveat names only `default_51` as displaced; the production condition `gray_bbox` is displaced on 3 of the 4 retried ROIs, including the worst ROI

Cell 24, caveats: *"Here the first draw was refused on `013.tiff`, `245.tiff`, `300.tiff` and `403.tiff`… **On
those four, `default_51` is scored on a click it would not have drawn alone.**"* Cell 0 describes `gray_bbox` as
*"what `build_seed` does in production today (D8)."*

From my own walk over the pixels (`..._tier_b_seed_draws.csv`):

| ROI | first draw | `gray_inverted` gate | `hematoxylin_od` gate | production (`build_seed` on gray) would click | joint click |
|---|---|---|---|---|---|
| 013.tiff | 254 | **accepted** | refused (click not foreground) | **254** | 249 |
| 245.tiff | 6274 | **accepted** | refused (click not foreground) | **6274** | 6317 |
| 300.tiff | 14581 | **accepted** | refused (click not foreground) | **14581** | 14499 |
| 403.tiff | 20334 | refused | refused | 20356 | 20356 |

On 013, 245 and 300, the hematoxylin gate alone refused the first draw. `gray_bbox` there is D8's template
construction on a click production would not have drawn. 245.tiff is the ROI the "Worst ROI" paragraph reports
(0.20 → 0.10, 0.23 → 0.17). The reference notebook states this itself (joint = gray-only draw on 11/14). This
notebook's caveat does not. It moves no number and no paired delta, since all three conditions share the click, but
it narrows what "`gray_bbox` = production" means on 3/14 ROIs.

**Fix:** *"On 013, 245 and 300 only the `hematoxylin_od` gate refused the first draw, so `gray_bbox` (production)
is also scored on a click production would not have drawn; on 403 both gates refused it and only `default_51` is
displaced."*

#### T2-2 — Context the reader needs beside the headline (extensions, not defects)

**Robustness of "`hem_bbox` vs `default_51` resolves"** (K = 20, Holm-adjusted p = 9 × 3/1024 = 0.0264):

| variation | K=20 exact p | K=30 exact p | Holm survivors (of 9) |
|---|---|---|---|
| notebook (consensus labels, mixed pool) | 3/1024 = 0.0029 | 7/1024 = 0.0068 | K=20 only |
| mitotic-only matching / maximum matching / r ± 1 px | identical | identical | K=20 only |
| contested hits neither credited nor charged: P = unanimous TP / (K − contested TP) | 0.0034 | 0.0032 | **K=20 and K=30** |
| contested mitoses removed from GT (their hits become FPs) | 33/512 = 0.064 | 0.037 | **none** |
| domain as the unit (7 domain sums of 2 ROI deltas) | 0.031 (floor 0.031) | 0.0625 | **none** |

The consensus-label convention the notebook uses is the right one. The headline survives the neutral treatment of
contested figures and fails only if 2-of-3 mitoses are charged as errors. It depends on treating ROI, not domain, as
the exchangeable unit: at 7 domain units the floor is 1/32. See Part 4.

**`tm_score` arm (my Tier B, no notebook artifact).** Validated against the reference raw CSV on 112 `tm_score`
values, 0 divergences. The `default_51` column is audit-derived only:

| K | pooled P `default_51` | `gray_bbox` | `hem_bbox` | `gray − default` W/L/T, p | `hem − default` W/L/T, p | worst ROI P (d / g / h) |
|---|---|---|---|---|---|---|
| 10 | 0.629 | 0.536 | 0.521 | 3/6/5, 0.17 | 4/6/4, 0.17 | 0.10 / 0.10 / **0.20** |
| 20 | 0.568 | 0.471 | 0.443 | 2/8/4, 0.078 | 2/11/1, 0.061 | 0.10 / 0.10 / **0.20** |
| 30 | 0.512 | 0.443 | 0.407 | 2/8/4, 0.0625 | 2/11/1, 0.054 | 0.10 / 0.10 / **0.17** |

The caveat *"Nothing here carries over to `tm_score`"* is a fair hedge. Measured, the pooled ordering does carry over,
with gaps 2–7× larger (the D8 production refinement costs 7–10 points against the plain 51 px box). On the worst ROI,
the metric D4 cares about, the order reverses: `hem_bbox` is best, driven by 245.tiff (0.30 / 0.35 / 0.30 vs 0.10).
One click per ROI, no test survives multiplicity, and this is not a claim the notebook makes. It is the number a
reader weighing D8 or D5 would want next to this one.

### Tier 3

- **T3-1 — `default_51` is labelled "the pre-D8 default" (cell 0 table). It was not.** Before D8, production cut an
  Otsu-gated, tightened-size, click-centred template: `experiment.run_one_image(tighten_bbox=True)` →
  `tightened_base_size`, and D8's superseded text says *"Only the size was corrected before."* An untightened 51 px
  box is the pre-*tightening* default (`tighten_bbox=False`, `FSConfig.base_size`). **Fix:** "untightened 51 px box
  (pre-tightening default)", and note the click-centred tightened-size arm is absent.
- **T3-2 — "`check_nms_radius` passed 42 times" (cells 5 and 24) is vacuous evidence.** `production.py:48` sets
  `nms_radius = ev.radius_px(mpp)` and `:55` checks that same value against `ev.radius_px(mpp)`, so it cannot fail as
  wired. D7 compliance holds by construction (verified by reading, and by Tier B using an independent 7.5 µm NMS).
  **Fix:** cite the construction, not the check.
- **T3-3 — Holm step 2 is "0.0068 × 8 = 0.054"; exact is 7/1024 × 8 = 0.0547 → 0.055.** The conclusion is unchanged.
  Relatedly, Figure 3 prints `p=0.062` for 1/16, where the prose says 0.063.
- **T3-4 — "about one test in ten would land below 0.10" (cell 16) overstates chance.** Under the sign-flip null
  conditional on each test's |deltas|, the expected count below 0.10 is **0.42 of 9**. Three tests cannot go below
  0.10 at all: p-floors 1/8, 1/8 and 1/4 for `gray − default` K = 10, 20 and `hem − gray` K = 10. The error is in the
  conservative direction. **Fix:** report `n_nonzero` and the floor 2/2^n beside each p.
- **T3-5 — "the budgets D9's `MAX_PEAKS = 100` was validated through" (cell 0) holds for K, not for this ranker.**
  D9 was measured on the `tm_score` arm. Under `chromatin_od` the cap decides which 100 score-ranked peaks get
  re-ranked, and that was never validated.
- **T3-6 — Provenance.** (a) Gate 2 compared against a version of `..._hembbox_raw.csv` that was overwritten at
  18:12:22, 32 minutes after this notebook ran. Re-run against the current file it still passes 28/28 on all seven
  fields, plus my 224 values. (b) `bbox3way_win_counts.png` on disk has mtime 17:33:37, which predates cell 23's
  recorded execution (21:40:28 UTC = 17:40:28). Its content matches the embedded figure number for number. (c) Every
  artifact is untracked. **Fix:** commit the notebook with its CSVs and figures, and record the reference CSV's hash
  in Gate 2.
- **T3-7 — "The losses cluster on 402, 529, 548 and 246"** omits 459.tiff (a loss at K = 30). "Cluster" still
  describes it fairly.
- **T3-8 — The reference text the tie paragraph quotes has since been replaced.** When this notebook ran
  (21:35–21:40 UTC), the reference notebook's summary did call the K = 10 ties a property of the production cap; the
  round-3 audit (17:39 local) quotes it. The reference was re-executed at 22:09–22:18 UTC and now attributes them to
  object convergence, and "structural property" no longer appears in it. The quote was accurate when written. A reader
  who opens the reference today will not find it. **Fix:** cite the reference's execution time, or drop the quote.

---

## Part 2 — verdict per conclusion

| # | cell | conclusion | verdict |
|---|---|---|---|
| 1 | 24 | Pooled precision table (9 values) | reproduces |
| 2 | 24 | `default_51 > gray_bbox > hem_bbox` at every K; neither refinement beats 51 px at any K | reproduces (pooled; per ROI `gray` wins 1/0/3) |
| 3 | 24 | `hem − default`: −3.6/−6.8/−4.3 pts, L 5/11/10, W 0/1/1 (013), p 0.063/0.003/0.007 | reproduces |
| 4 | 24 | Holm over 9 rejects only `hem − default` K=20; K=30 fails at the next step | reproduces (step-2 figure 0.054 → 0.055, T3-3; sensitivity in T2-2) |
| 5 | 24 | `gray − default` unresolved, mostly ties; −1.4/−2.9/−1.7, T 10/10/6, L 3/4/5, W 1/0/3, p 0.63/0.13/0.30 | reproduces |
| 6 | 24 | Ties partly by construction: 3 ceiling ROIs; 245 (0 px) and 201 (0.5 px rounds to same pixel) identical lists | reproduces (201: `round(tpl) = round(click)`, 30/30 shared) |
| 7 | 24 | 403 the only real same-size centring shift; top-30 changes completely, TP does not | reproduces (0/30 shared coordinates, yet the **same** 5/10/14 mitoses) |
| 8 | 24 | `hem − gray` reproduces the reference exactly: −2.1/−3.9/−2.6, W/L/T 0/3/11, 3/9/2, 3/7/4, p 0.25/0.090/0.098 | reproduces (against the current reference CSV too) |
| 9 | 24 | Worst ROI 245 under every condition at every K (tied 529 gray K=30); hem 0.20→0.10, 0.23→0.17 | reproduces; on 245 `gray_bbox` is not the production click (T2-1) |
| 10 | 24 | Precision falls in the order of median template size 51→40→28; median offset 2.5/3.5; not a dose-response | reproduces (hedged correctly) |
| 11a | 24 | "At the same cap, `hem − default` has 2 ties at K=20 and `gray − default` 10", so the cap alone does not explain ties | reproduces |
| 11b | 24 | "Delta ties follow how different the two templates are" | **does not reproduce**: refuted by `hem − gray` (smallest size gap, 2 ties at K=20) (T1-1); the quoted reference text has since been replaced (T3-8) |
| 12 | 24 | Caveat: `default_51` displaced on the four retried ROIs | reproduces, overstated by omission: `gray_bbox` is displaced on 013/245/300 (T2-1) |
| 13 | 24 | Caveats: one click, D5 bar unmet; `chromatin_od` only; K ≤ 30, domain steps 1/(2K) | reproduces |
| 14 | 5, 24 | D9 cap bound 42/42; K delivered; annulus empty 42/42; `check_nms_radius` passed 42 | reproduces from pixels; the `check_nms_radius` item is vacuous (T3-2) |
| 15 | 7–8 | Gates 1–3 pass (14/14, 28/28, 1/1) | reproduces (independently: 14/14 clicks from my walk; 28/28 against the current CSV; 245 and 201 identical) |
| 16 | 0 | Only template size and centre change; same channel and od window; self-hit at template centre; GT keyed to click | reproduces (code reading plus 42/42 lists from a re-implementation varying only `(base, tx, ty)`) |
| 17 | 0 | `default_51` has no gate beyond readability, which `border_filter` guarantees; same predicate as the reference | reproduces |
| 18 | 0 | `default_51` is "the pre-D8 default" | does not reproduce (T3-1) |
| 19 | 0 | K ≤ 30 = budgets D9 validated through | reproduces, overstated: D9 validated `tm_score` only (T3-5) |
| 20 | 2 | Matching the full list = matching the top K alone | reproduces (252/252 from full lists) |
| 21 | 9 | Table A prints (gray median 40, 25–51, <51 on 11/14, offset 2.50; hem median 28, 23–41, 14/14, 3.47; hem < gray 13/14, = 1/14) | reproduces |
| 22 | 16 | ~1 test in 10 below 0.10 by chance | reproduces, overstated (0.42 of 9; T3-4) |
| 23 | 18 | Sub-pixel offsets can round to the same patch | reproduces (201) |
| 24 | 21–23 | Figures 1–3 | reproduces (step 1a: 171 plotted values; step 2: renders read, legend/colour/stack order/scale correct) |

`cannot check` was not needed for any claim.

---

## Part 3 — what was re-run versus read

- **Tier A** (JSON, TIFF tags, the notebook's seven CSVs, the reference raw CSV): every table cell, every statistic,
  every prose number, plus alternative scoring conventions (`..._delta_stats_by_scoring_rule.csv`,
  `..._pooled_precision_by_scoring_rule.csv`), the tie analysis, and domain-level sign-flip. Sign-flip p-values use
  exact enumeration over 2^n_nonzero patterns (n ≤ 12, G = 14), not a bootstrap.
- **Tier B, from pixels, all 14 ROIs × 3 conditions** (~3.5–4.5 min wall clock per full run; three full runs): the
  RNG walk and Otsu gate from the D8 text (cv2 Otsu, skimage label/regionprops); `rgb2hed` hematoxylin OD;
  `cv2.matchTemplate(TM_CCOEFF)` on a replicate-padded image; robust-z floor; `scipy.ndimage.maximum_filter` 15×15
  local maxima; the lexsort tie rule; top 100; greedy NMS at 7.5 µm; self-hit ≤ 5 px at the template centre; od51 as
  the mean of the top 10 % in a 51 px replicate-padded window; stable sort. Also `tm_score`-ranked lists and K = 50.
  I chose these checks because the load-bearing questions (right click, right list, right TP) all live here.
- **Tier C:** not run. No divergence needed explaining, and the Tier B re-implementation already covers the whole
  measurement path.
- **Read, not re-run:** `production.py`, `find_and_suppress.py`, `template_match.py`, `seed_selection.py`,
  `evaluate.py`, `dataset.py`, `channels.py`, `chromatin.py`, `nms.py`, `invariants.py`,
  `production_pipeline/run_pipeline.py`, `experiment.py` (the pre-D8 path), and the reference notebook in full.
- **Audit script:** `bbox_refinement_three_way_chromatin_od_audit.py` ran start to finish, clean (exit 0, 234 s on
  the final run) and regenerates all 23 `results/bbox_refinement_three_way_chromatin_od_audit_*.csv`.
  `--tier-a-only` skips pixels.
- **Beyond the named Step 4 modes, I looked for:** a look-alike stealing a TP (none within `r` of both classes); X/Y
  resolution anisotropy (present, immaterial); annotations near the click that the self-hit NMS disc could erase
  (nearest 62 px); tied od values (largest tie block 1 on 42/42 full lists); whether Gate 2's file still exists as
  read (it does not; still passes); a stale figure on disk (T3-6); which gate refused each retried click (T2-1);
  whether the headline depends on the ranker (T2-2); and whether "resolves" depends on contested labels or the unit
  of exchange (T2-2).
- **Step 4 modes.** (1) Unit: ROI-level exact sign-flip, correct; floors in T3-4, domain-level in T2-2. (2) Domain
  confound: N/A, because treatment is within-ROI. (3) Recall triad: precision@K only. `n_detections` is reported;
  recall@K is algebraically redundant with precision@K here (identical `n_gt` across conditions); ~97 points cover under 1 % of the ROI at one match radius, so saturation is not live. (4) Length-matched null: N/A, no comparison against chance is
  made. (5) Scope: T1-1 (a K = 10 claim answered with K = 20 counts), T3-5, T2-2. (6)/(7) Pre-registration / decision table: N/A, none exists
  (`Research Logs/*preregistration*` covers F1, F2, F4, F5, F6 only). (8) Product: worst ROI reported; worst-seed
  unmeasurable at one click; T2-2 `tm_score` worst-ROI reversal. (9) Machinery: none superfluous. (10) In-sample
  selection: no operating point is chosen, since the three conditions are fixed a priori. The ranker is 1 of the 3
  arms the reference computed. (11) Timing: `t_pipeline_s` and "ran in 307 s" are single unreplicated runs and are not
  interpreted; my search stage took 1.8–10 s per run on this machine. (12) Above.
- **Family.** Established by grepping notebook sources for artifact stems (`bbox3way`, `hembbox`); git history is
  unavailable because the whole folder is untracked (`git log --all -- production_hematoxylin_only/*` is empty).
  Sibling arm: `production_seed_precision_at_k_chromatin_hem_bbox.ipynb` (hem vs gray × 3 rankers × K ≤ 50; this
  notebook reuses its gate and its two conditions and adds `default_51`). Excluded as not-an-arm:
  `bbox_refinement_three_way_detail_14roi.ipynb` (visual companion reading this notebook's CSVs),
  `tightening_process_hem_vs_gray.ipynb` and `production_seed_tightening_detail_14roi.ipynb` (mechanism figures,
  no precision output).
- **Prior audits, read last.** None covers this notebook. `2026-09-16-hembbox-precision-at-k-audit-round3.md`
  covers the reference notebook's shared material. T1-1 **extends** its F1 (object convergence) to three pairs; T3-8 records that the reference text
  this notebook quoted was replaced, after this notebook ran, in response to that F1. T2-1 **confirms** its Part 0 draw
  table (013→254, 245→6274, 300→14581) and **extends** it to this notebook's caveat. The `tm_score` worst-ROI
  reversal **confirms** its F4 and adds the `default_51` column. Its §I says passing `category_id=MITOTIC` "would in
  fact be wrong"; I **disagree in part**: both conventions are defensible, and here they are identical on 126/126
  cells. Its F7 (K = 50 outside D9) does not recur, since this notebook stops at 30. Rounds 1–2 and the
  tightening-process audit predate this notebook, and I did not check them against it.
- **Appendix facts, checked:** `FSConfig.tm_method` still defaults to `TM_CCOEFF_NORMED` (holds);
  `image_annotations(category_id=None)` (holds); radius 29.6–33.1 px (holds). The per-pass Tier B cost is in the
  appendix's range. Not re-verified: repo-wide PNG, notebook and CSV counts. New house fact: `timeout` is not
  installed on this machine.

---

## Part 4 — premises this audit inherited

ROI-per-stratum count in this sample: **exactly 2 ROIs in each of 7 tumour domains** (JSON). The exchangeable unit is
**not** collinear with the stratum, but within-stratum replication is only 2. Effective units per test
(`n_nonzero`): 4, 4, 8, 5, 12, 11, 3, 12, 10.

| premise | source | what would falsify it |
|---|---|---|
| ROI is the exchangeable unit | D5; F5 §8 | ROIs within a domain or scanner behaving as one unit. At the domain level, `hem − default` K=20 is p = 0.031 at its 1/32 floor and survives no Holm step (T2-2). Worth a premise review: 14 ROIs from 7 domains × 5 scanners. |
| Match radius = NMS radius = 7.5 µm | D7; `ev.MIDOG_RADIUS_UM` | TP counts moving with radius. Here ±1 px moves 0/126. |
| Consensus `category_id` is truth; contested mitoses are mitoses | MIDOG++ labels; `evaluate.py` | Conclusions flipping when 2-of-3 figures are charged as FPs. They do (T2-2); under the neutral reading they strengthen. |
| `chromatin_od` is the ranker to evaluate | `PRODUCTION_PIPELINE_CLEANUP.md` (contradicts D5) | The ordering reversing under `tm_score`. Pooled it does not; at the worst ROI it does (T2-2). |
| precision@K, K ≤ 30, worst ROI are the product metrics | D4 (asks recall@K and worst click per domain) | Recall@K diverging from precision@K (it cannot here: same `n_gt`), or worst-click behaviour differing from worst-ROI (unmeasurable at one click). |
| `MAX_PEAKS = 100` leaves precision unchanged through K = 30 | D9 (`tm_score` arm only) | An uncapped `chromatin_od` run giving different top-30 lists. Not run. |
| One click per ROI (`seed_index = 0`) is representative | this notebook family | A 5-seed sweep (D5's bar) moving the paired deltas. |
| Prior audits' findings | `Research Logs/2026-09-16-hembbox-*` | Same premises and the same day, so not independent corroboration. |

**Which verdicts would change if a premise failed:** verdicts 3–4, "`hem_bbox` vs `default_51` resolves" under Holm,
depend on ROI exchangeability and on not charging contested mitoses as errors. Verdicts 1–2 and 9, the ordering and
the worst ROI, depend on `chromatin_od` being the ranker, and the worst-ROI ordering reverses under `tm_score`. Every
arithmetic verdict in Part 0 is independent of these premises.
