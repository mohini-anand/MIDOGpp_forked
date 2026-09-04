# Audit: the three `tm_threshold_axis_sweep` notebooks

Date: 2026-09-03
Scope: `tm_threshold_axis_sweep.ipynb` (v1), `tm_threshold_axis_sweep_v2.ipynb` (v2),
`tm_threshold_axis_sweep_largest_cc.ipynb` (lcc), and their three saved CSVs. Everything
below is re-derived from `results/tm_ccoeff_threshold_axis_sweep{,_v2,_largest_cc}.csv`
plus two targeted re-runs on `301.tiff`; no sweep was re-executed.

These notebooks are unusually self-critical and most of their machinery is correct (see
"Checked and sound" at the end). The findings are graded so the ones that change a stated
conclusion are not flattened in with the ones that don't.

---

## Tier 1 -- changes a stated conclusion

### 1. The headline metric is saturated by geometry, and the harness's own guard was never read

`evaluate.py`'s module docstring: *"Read `coverage_frac` before reading any full-list
number. A detection list long enough to tile the ROI answers 'is this annotation within
the match radius of some detection?' by geometry rather than by evidence."* All three
CSVs carry the column. None of the three notebooks prints or mentions it.

`coverage_frac` (fraction of arbitrary ROI locations already within the match radius of
some detection), `tm_score` axis:

| domain | v1 @ z=-1 | v2 @ z=-1 | lcc @ z=-1 | lcc @ z=1.0 (operating point) |
|---|---|---|---|---|
| mast cell tumor | 0.9496 | 0.9936 | 0.9944 | 0.9739 |
| lymphosarcoma | 0.9505 | 0.9982 | 0.9987 | 0.9921 |
| breast cancer | 0.9442 | 0.9964 | 0.9994 | 0.9417 |
| neuroendocrine | 0.9371 | 0.9918 | 0.9990 | 0.9824 |
| soft tissue sarcoma | 0.9265 | 0.9675 | 0.9954 | 0.9649 |
| melanoma | 0.9266 | 0.9737 | 0.9967 | 0.9390 |
| lung cancer | 0.9293 | 0.9805 | 0.9805 | 0.8957 |

v2 exists to answer "does every domain now reach full recall", and lcc's first verdict
bullet is "recall unchanged at the operating point, strictly better beyond it". Both rest
entirely on `full_list_recall`, measured on lists that cover 90-99.9% of the ROI. At that
saturation `full_list_recall == 1.0` is close to a tautology, and the correct framing of
v2's result is "the two fixes raised ROI coverage from ~0.95 to ~0.99, which is why the
last few annotations stopped being missed" -- not "the pipeline now finds everything".

### 2. lcc reports one effect twice, as two independent findings

lcc's summary calls the candidate-volume increase "the one clean, confound-free effect in
this run", and separately presents "recall >= v2 in every (domain, z) cell" as a benefit
of tightening. They are the same effect. Every positive recall delta tracks the coverage
and volume delta:

| domain | z | n_det ratio (lcc/v2) | coverage delta | recall delta |
|---|---|---|---|---|
| soft tissue sarcoma | 3.0 | 2.06 | +0.229 | +0.046 |
| mast cell tumor | 3.0 | 1.69 | +0.170 | +0.152 |
| neuroendocrine | 3.0 | 1.91 | +0.183 | +0.039 |
| lung cancer (null control) | 3.0 | 1.00 | 0.000 | 0.000 |

Where the template did not change (`201.tiff`, `untightened_51px=True`) all three columns
are exactly zero at every z. "Tightening never costs recall and recovers some at aggressive
cutoffs" is therefore not independent evidence -- it is the volume increase, scored through
a metric the volume increase saturates.

### 3. v2's cost table is dominated by a bug v2 itself introduced, and was never corrected

`nan_rate` for `arm='chromatin_od'`, across all rows:

* v1: **0.000000** (min = median = max)
* v2: **0.0167 - 0.0301** (61-746 candidates per (domain, z) cell; 566-746 per domain at the deep floor)
* lcc: 0.000000

The cause is v2's own Fix 1. Padding made the ROI border reachable for the *match*, but
`chromatin.score_detections` still read the un-padded channel, so every newly-reachable
border candidate got `od = NaN`; `compare._rank` sorts NaN last (behaviour already
documented in `chromatin.rerank`: *"nan (border) sorts last"*). v2's Fix 1 removed the
border protection that had been keeping `od` NaN-free by accident.

v2 cell 22/23 attributes the whole read-depth blow-up to the smaller NMS radius
("shrinking the radius is not a free lunch"). Its largest single entry is canine lung
cancer, `chromatin_od`, `read_95`: 653 -> 19,373 (+18,720). But `201.tiff` in lcc uses the
same seed and a byte-identical 51 px template with `od` fixed:

| 201.tiff, chromatin_od, z=0.5 | read_90 | read_95 |
|---|---|---|
| v1 (no padding, no NaN) | -- | 653 |
| v2 (padded, od NaN) | 16,142 | 19,373 |
| lcc (padded, od fixed) | 164 | 814 |

The real cost of the radius change on that cell is about **+161**, not +18,720 -- a ~116x
overstatement. `tp_fp_chromatin_and_raw_score_distribution.ipynb` did spot the NaN gap
("the NaN gap the border fix reopens for this second signal") and lcc notes v2's
`chromatin_od` read depths are contaminated -- but nobody went back and retracted v2's
stated conclusion, which is the one a reader of v2 walks away with.

### 4. v1's cross-domain summary silently changes its domain set as z rises

`compare._read_depths` returns NaN when a recall fraction is unreachable. v1 cell 21
hand-rolls `.median()` over those rows instead of using `compare.summarise`, which carries
`n_unreachable` for exactly this case.

Domains contributing to `median_read_95` (of 7): z <= 2.0 -> 7; **z = 2.5 -> 5; z = 3.0 -> 4**.

So the apparent improvement in reading depth at z=2.5 (`tm_score` 3158 -> 2591;
`chromatin_od` 2613 -> 1768) and the non-monotone bounce at z=3.0 (1768 -> 2040.5) are
survivorship: the domains needing the deepest reading drop out of the median precisely
because they can no longer reach 95% recall at all. The notebook calls this "the single
table to defend a `z` choice from".

### 5. Every read-depth delta interpreted in v2 and lcc sits inside single-seed noise

All three notebooks use `SEED_INDEX = 0` only, one ROI per domain.
`compare.summarise`'s docstring: *"Never report a single-seed number (Step 0 rule 5)."*

`results/tm_ccoeff_headtohead_seed_variance.csv` -- same 7 ROIs, same `TM_CCOEFF`,
5 seeds, `decision_grade=True`, `pool_scope='full'` -- gives `read_95` seed variance:

| ROI | mean | sd | min | max | v1's seed-0 value |
|---|---|---|---|---|---|
| 301.tiff | 6606 | 2208 | 4753 | 10389 | **10389 (the max)** |
| 246.tiff | 4379 | 1635 | 2591 | 6546 | **2591 (the min)** |
| 548.tiff | 4282 | 1325 | 3158 | 6130 | **3158 (the min)** |
| 459.tiff | 7620 | 1199 | 6280 | 9021 | 7047 |
| 402.tiff | 2327 | 593 | 1707 | 3276 | 2446 |
| 094.tiff | 1330 | 260 | 943 | 1658 | 1440 |
| 201.tiff | 4395 | 2744 | 2181 | 7465 | 3539 |

**The load-bearing part is the three endpoint hits.** Seed 0's `read_95` lands exactly on
the *maximum* of the 5-seed distribution for `301.tiff` and exactly on the *minimum* for
both `246.tiff` and `548.tiff`. Three exact coincidences with distribution endpoints
cannot be an artefact of configuration differences -- they identify seed 0 as an
unrepresentative draw, in opposite directions depending on the domain, whatever else
changed between the runs.

The SD comparison supports that: lcc's `read_95_delta` values (+1416, +527, -93, -717,
-1199, -1547, 0) are every one of them smaller than their own ROI's seed SD. lcc's "mixed,
no consistent direction ... no domain property visible in this run that predicts the sign"
is describing sign noise, and the same applies to v2's cell 23. (The variance CSV comes
from the earlier `tm_variant_sweep` configuration -- no padding, NMS at the match radius,
a different pool scope -- so treat the SDs as a scale reference, not an exact error bar.
The endpoint argument above does not depend on them.)

---

## Tier 2 -- real slips, conclusion survives

### 6. v2's leak accounting is off by one and conflates "unreachable" with "missed"

Re-ran `301.tiff` at the deep floor. v1 misses exactly **6** mitoses
(0.97235 x 217 = 211 found), and v2's own recall delta agrees (0.0276 x 217 = 5.99).
v2's prose says 2 border-stranded + 5 NMS-suppressed = 7.

Measured: 2 evaluation mitoses have `valid=False` at their click pixel (14779, 15022) --
but only **14779** was actually missed; 15022 was still claimed by a detection inside the
29.6 px match radius. The correct split is **1 border + 5 NMS = 6**. "Stranded" was used
to mean "unreachable at the click pixel", which does not imply "not found". The 5 NMS
cases are exactly v2's `_missed_ids`, and the fix itself is sound.

### 7. lcc's coincidence probability uses the wrong denominator

`coincide_rate = n_coincide_match / n_tighten_ok`, but the seed is drawn uniformly from
the whole border-filtered pool and `draw_seed_with_retry` only rejects on `largest_cc_box`
failure -- so P(a drawn seed is a coincide case) is `n_coincide / n_pool`:

| domain | n_pool | n_coincide | P |
|---|---|---|---|
| 301 | 162 | 130 | 0.8025 |
| 201 | 7 | 6 | 0.8571 |
| 246 | 97 | 82 | 0.8454 |
| 459 | 121 | 115 | 0.9504 |
| 094 | 49 | 43 | 0.8776 |
| 548 | 184 | 165 | 0.8967 |
| 402 | 65 | 60 | 0.9231 |

Product = **0.4015**, not the reported 0.8699. 7/7 is still unremarkable at 40%, so the
"not a coincidence needing separate explanation" reading holds; the number is 2.2x
optimistic.

### 8. v2 re-violates `check_nms_radius` by design and reports only half its cost

`invariants.py`'s docstring names the exact defect: *"a constant 25.0 px NMS radius let
two detections 26 px apart both survive while both sat inside one ground-truth object's
match radius (29.6-33.1 px across these scanners), inflating the FROC's false-positive
axis with duplicates."* v2 passes `nms_radius=None` to skip the check. The decoupling
argument is sound, and v2 does warn about more candidates and deeper reading -- but not
that the surviving duplicates are bucketed as **unannotated false positives**.

Measured on `301.tiff`'s deep pool, detections inside the match radius of some GT that
bucket as `non_human_findings`:

* v1 config (NMS = match radius): 111
* v2 config (NMS = 5 um): **281**

That is 1.1% of the FP pool -- it does not overturn any number -- but it feeds
`precision_full` in v1's composition figure and `frac_fp_in_tp_range` in lcc's
deliverable (f), neither of which adjusts for it.

---

## Tier 3 -- worth noting

### 9. `assert_floor_not_limiting` cannot fire in any of the three

`DEEP_FLOOR_Z = min(Z_LEVELS) - 0.5` makes `cut < deep_floor` impossible for every swept
z, so `floor_limited` is False by construction. v2's summary lists it among the things
"verified, not assumed".

### 10. v1's strongest, most decision-relevant result is never stated

`recall_at_budget` is the metric `compare.py`'s own docstring designates as **primary**.
On the `tm_score` axis at K=100 it is *exactly identical* at all nine z levels in every
domain (mast cell 0.2627 x9, lung 0.5882 x9, lymphosarcoma 0.4174 x9, ...). That is a
mathematical identity -- a score-descending list's top-K cannot change when you delete
only lower-scoring entries, *provided the surviving list is still longer than K*. At
K=100 it always is: the smallest pool anywhere in v1's grid is 3,070 (canine lung cancer
at z=3.0), 31x the budget, so K=100 is never near the truncation boundary at any swept z
in any domain. v1's cell-19 figure therefore asks an empirical question that, at the
budget it plots, has no other possible answer.

The identity is budget-scoped, not universal: at K=5000 the pool does truncate at z=3.0
(3,070-3,424 < 5,000) and `recall_at_budget` moves in 5 of 7 domains, by up to 0.0922
(mast cell tumor 0.8387 -> 0.7465). So the correct statement is "constant across z
wherever `n_detections >= budget`", which covers every budget in `cp.BUDGETS` except
5000. On the `chromatin_od` axis, where
the question is genuinely live, the entire 4-unit z sweep moves recall@100 by at most
0.0096 (under 1 percentage point) in 4 of 7 domains, and by exactly 0 in the other 3.

The actionable sentence -- "across the whole swept range the extraction floor is
irrelevant to the product metric" -- is nowhere in the notebook, which instead defers the
call to "whoever is setting the product's reading-depth budget".

### 11. lcc downgrades a measured check to an assertion (it does still hold)

v1 *measured* axis-independence of full-list counts (0/63 cells) precisely because greedy
rank-order matching can flip contested GT. lcc asserts it "since the same candidates just
get re-sorted", without re-checking at the new radius and pool. Re-checked here: **0/63 in
v2 and 0/63 in lcc** -- the assertion is true, just unverified where it was made.

### 12. lcc's volume mechanism is only partly borne out

"A smaller template has a broader, less-selective correlation response (more of the map
clears the same per-map z cut)". At the deep floor `301.tiff`'s pools are essentially
identical (v2 25,449 vs lcc 25,447) despite the template shrinking 51 -> 41 px; the
divergence there appears only at higher z. `459.tiff` does behave as described
(18,829 -> 26,152 at the deep floor). At low z the count is bounded by peak-spacing and
NMS packing geometry, not by the threshold, so the stated mechanism operates in some
domains and not others.

---

## Checked and sound

* **v1 -> v2 z units are genuinely comparable.** `map_median`/`mad_scale` change by
  <0.5% (mad ratios 0.996-1.000), so "the same z" really is the same cut across v1 and v2;
  padding does not perturb `robust_stats`. (lcc vs v2 is different: mad ratios 0.41-1.00,
  which lcc does flag for the separability metric.)
* **The NMS suppression diagnostic holds in all 5 cases.** No exact ties -- the closest
  (ann 14736) is 0.814998 vs 0.814613, a rounding artefact in the printed table, not a
  tie broken by sort order. Every suppressor genuinely outscores the own-peak, and
  `same_annotation=True` throughout.
* **The 36.8 px minimum annotation spacing claim verifies exactly** (548.tiff). Dataset-wide
  minimum mitotic-to-mitotic spacing is 26.6 px, consistent with `find_and_suppress.py`'s
  documented 26.2 px over all annotation pairs. Every `nms_radius_px` (19.7-22.1) sits
  below both.
* **The one-match-many-z shortcut is exact by construction**, given
  `extract_peaks`' global `np.lexsort((ys, xs, -scores))` key and `nms_by_distance`'s
  `kind="stable"`. Verified 7/7 in each notebook (at CURRENT_Z only, but the argument
  covers every z).
* **lcc's transcription of the 7 v2 reference separability rows** from
  `tp_fp_chromatin_and_raw_score_distribution.ipynb`: all 14 values match that notebook's
  executed output exactly.
* **The `201.tiff` null control is exact.** All 13 shared `tm_score` columns
  (`n_detections`, `read_50/80/90/95/99/100`, `coverage_frac`, `n_lookalike_in_list`,
  `map_median`, `mad_scale`, `n_pool`, `full_list_recall`) are identical between v2 and
  lcc at z=1.0.
* **`OD_PAD = tm.BASE_SIZE // 2 = 25` is exactly right** -- `chromatin.score_detections`
  defaults `window=tm.BASE_SIZE`, so 25 px is the minimum sufficient margin, and giving it
  its own constant rather than reusing the template-dependent `PAD` was the correct call.
* **Fix 1's `valid.all()`-by-construction argument is correct** and measured
  (0.981983 -> 1.000000 on 301.tiff).
* **`largest_cc_box` faithfully mirrors `tighten_box_otsu`'s binary path** with only the
  component-selection step changed, and the synthetic retry test covers both branches.

## Reproducing

CSV-only checks (1-5, 7, 9-12) run directly against
`results/tm_ccoeff_threshold_axis_sweep{,_v2,_largest_cc}.csv` and
`results/tm_ccoeff_headtohead_seed_variance.csv`. Checks 6, 8 and 13 need one re-run of
`301.tiff` under both configurations (~10 s of matching each, plus ROI load).
