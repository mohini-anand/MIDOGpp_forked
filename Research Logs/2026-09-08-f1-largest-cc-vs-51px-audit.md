# Audit: the largest-CC vs fixed-51 px template comparison — is it apples to apples?

Date: 2026-09-08.
Scope: the two runs that compare a largest-CC-tightened template against the fixed 51 px box
before template matching, and specifically whether their comparison is like-for-like.

* `tm_threshold_axis_sweep_largest_cc.ipynb` (lcc) vs `tm_threshold_axis_sweep_v2.ipynb` —
  single seed per domain, 7 ROIs. Already audited twice
  (`Research Logs/2026-09-03-tm-axis-sweep-audit.md`,
  `Research Logs/2026-09-04-largest-cc-z-headroom-audit.md`); those findings were spot-checked
  here and hold. Not re-derived.
* **`f1_seed_sweep.py` + `f1_analysis.py`** — the paired 5-seed ablation, and the only run that
  is *designed* as a controlled comparison. This audit is about that one.

Everything below is re-derived from `results/f1_seed_sweep.csv`,
`results/f1_seed_sweep_cells.csv` and the committed source. `f1_analysis.py` was re-executed;
no sweep was re-run.

---

## The answer to "is it apples to apples"

Yes on every axis except one, and the exception is the one that moves.

**Held fixed and verified:** the seed click (same `ann_id` in both arms — the same
`draw_seed_with_retry` draw feeds both), the evaluation GT with that seed removed, the
channel (`hematoxylin_od`), the method (`TM_CCOEFF`), `peak_min_distance`, `self_hit_radius`,
the decoupled 5.0 µm NMS radius, the chromatin `OD_PAD = 25` window, the deep floor
(`DEEP_FLOOR_Z = -1.5`), the budget grid and the matching radius. Only `base_size` differs.
GATE 2 confirms it: in all 10 null cells (`tightened == 51`) the two arms are equal on every
column.

**Not held fixed:**

1. `PAD = (base_size-1)//2` — 25 px for the 51 px arm, 9-12 px for a 19/25 px template. Both
   arms reach `valid.all()`, so the *property* is fixed and the *amount* is not; within ~25 px
   of the ROI edge the arms score against `BORDER_REPLICATE` fabrication of different widths.
   **Disclosed in both the docstring and the pre-registration.** Correctly handled.
2. **List length.** Both arms are cut at `median + z·MAD` of *their own* response map. Equal `z`
   is not equal selectivity: over the 25 informative cells the `largest_cc` arm's candidate list
   is longer in every single one, mean **+3,851** candidates (min +356, max +7,460). The
   comparison is matched on **z**, not on **n**.

So the honest one-line answer: **it is a clean paired ablation matched on the z threshold, not
on candidate volume.** A genuinely like-for-like contrast would match list length or coverage.
That correction was never run — §2 below runs an approximation of it, and the conclusion
survives.

---

## Tier 1 — changes reported evidence

### F1.1 The Robustness table in `2026-09-04-f1-results.md` reports a test the pre-registration did not specify

Every one of its five p-values is scipy's **asymptotic** Spearman p, not the pre-registered
**sign-flip permutation** p. Five-for-five to four decimals:

| subset | n | rho | log's p | scipy asymptotic | **pre-registered sign-flip** |
|---|---:|---:|---:|---:|---:|
| all cells (primary) | 35 | +0.5986 | 0.0001 | 0.0001 | **0.0012** |
| informative only | 25 | +0.5591 | 0.0037 | 0.0037 | **0.0019** |
| all zero deltas dropped | 23 | +0.5157 | 0.0118 | 0.0118 | **0.0043** |
| nulls + seed 0 dropped | 19 | +0.4537 | 0.0510 | 0.0510 | **0.0261** |
| all zero deltas + seed 0 dropped | 17 | +0.3619 | 0.1535 | 0.1535 | **0.0824** |

Two consequences:

* The row **"nulls + seed 0 dropped" flips sides**: 0.0510 (fails 0.05) → **0.0261** (clears).
* The log's headline caveat — *"the strictest subset does not clear 0.05 at n=17"* — remains
  true but is overstated roughly 2×: p = **0.0824**, not 0.1535.

This also explains the log's internal contradiction, flagged during this audit: the headline
says the primary is p = 0.0012 and the Robustness table's first row restates the *same test* at
p = 0.0001. With `N_PERM = 10,000` the permutation p floor is 1/10001 ≈ 0.0001, so 0.0012 (≈11
exceedances) and 0.0001 (zero exceedances) cannot be the same computation. They are not: one is
the permutation test, the other is scipy's asymptotic approximation.

**The verdict is unaffected.** The two pre-registered decision-rule inputs — primary p = 0.0012
and seed-0-excluded p = 0.0119 — both reproduce exactly from `f1_analysis.py` and both clear
0.05, so the run still lands in row 1 of the decision table. The defect is in the *reported
evidence*, not the conclusion.

### F1.2 The Robustness table has no committed provenance

`f1_analysis.py` calls `report_primary` exactly twice — all-35 and seed-0-excluded. It never
computes the informative-only, zero-deltas-dropped, or nulls-and-seed-0-dropped subsets. The
whole Robustness section, and the Wilcoxon in the effect-size section, were produced ad hoc
outside the committed script. That is how a non-pre-registered test got into a
pre-registered write-up without anyone noticing.

The Wilcoxon does not reproduce either: the log says p = 0.0101; the data give **0.0121**
(zeros dropped, n = 23) or **0.0145** (`zero_method="zsplit"`). Same provenance bucket — small,
but it is a published number no committed code produces.

**Fix:** move all five robustness subsets into `f1_analysis.py` using `signflip_spearman`, and
re-issue the table.

### F1.3 The pre-registration's justification for `recall_at_budget` is false on the axis it made primary

The pre-registration rejects `read_95` and picks `recall_at_budget` because the latter is
*"length-invariant while `n_detections >> K`"*, so that a delta cannot *"re-report the
candidate-volume difference under a new name"*.

That is true on `tm_score`, where the membership gate and the ranking key are the same quantity
(the audit of 2026-09-03 measured `recall_at_budget` as *exactly* constant across all nine z at
K=100 and K=250). It is **not** true on `chromatin_od`, the pre-registered primary axis: there
`score >= cut` decides membership and `od` decides order, so a longer list injects more
competitors into the top-K. Measured at z=1.0, K=250:

```
rho(area_ratio,      Δrecall@250) = +0.5986   <- the pre-registered primary
rho(Δn_detections,   Δrecall@250) = -0.7732   <- stronger
```

The metric chosen *because* it was immune to candidate volume is more strongly associated with
candidate volume than with the dose it was chosen to measure.

**It also mis-cites `compare.py`.** `midog_utils/compare.py:29-32` designates
`recall_at_budget` primary for the *opposite* reason: `read_*` metrics *"read from the top of
the list, so truncating a list barely moves them, which makes them a poor control when the
thing under comparison is the candidate generator."* `compare.py` picks it for its
**generator-sensitivity**. The pre-registration cites it for length-*invariance*.

### F1.4 …but the dose-response is not a volume artefact — checked three ways, and it survives

The obvious inference from F1.3 — that the whole result is dilution — is **wrong**, and it is
worth recording that it was tested rather than assumed.

**(a) Dilution slope, estimated where the treatment is not varying.** Fitting
`d(recall@250)/d(n_detections)` *within* each arm across the z grid gives a median base51 slope
of **−9.3e-07 recall per candidate**. Applied to each cell's own Δn it predicts a mean delta of
**−0.0036** against an observed **−0.0180**: volume explains ~**20%** of the effect. Subtracting
the dilution prediction cell by cell leaves `rho(area_ratio, residual) = +0.559` on the 25
informative cells — indistinguishable from the raw +0.5591.

**(b) Length-matching, both directions.** Using the z grid to pick, per cell, the cut that
equalises the two arms' list lengths (median residual length error 5-7%):

| matching | mean Δ | rho(area_ratio, Δ) |
|---|---:|---:|
| none (as published) | −0.0180 | +0.559 |
| raise z on `largest_cc` | −0.0151 | +0.487 |
| lower z on `base51` | −0.0179 | +0.539 |

**(c) What a partial correlation says, and why it is the wrong tool here.** Residualising both
variables on Δn flips the sign: `partial rho(area_ratio, Δrecall | Δn) = −0.089` on the
informative cells. This is **over-adjustment**, not a refutation: Δn is a *post-treatment
mediator* — shrinking the template is what makes the list longer — and it is 0.88-collinear
with `area_ratio`. Conditioning on it removes most of the treatment itself. Recorded here so the
next reader who computes it does not mistake it for a finding. (a) and (b) are the valid tests.

**Net:** F1.3 stands as a defect in the pre-registration's stated *rationale*; the conclusion
it was used to reach is independently sound. The "don't tighten below ~36 px" threshold is not
an artefact of candidate volume.

### F1.5 The "axis specificity" reading has a simpler competing explanation that was never considered

The results log reads the dose-response appearing on `chromatin_od` (rho +0.599) and not on
`tm_score` (+0.164, p = 0.399) as *"only the monotone structure is specific to `chromatin_od`,
the ranker that can see chromatin density."*

The competing explanation is that `chromatin_od` is the only axis whose recall@K is
length-sensitive at all. The evidence does not cleanly separate them, and this audit's own
measurements sit on both sides:

* For the volume account: `rho(Δn, Δrecall)` is −0.773 on `chromatin_od` but only −0.342 on
  `tm_score`, while `rho(area_ratio, Δn)` is *identical* (−0.882) on both — the arms share a
  pool, so the axes differ only in length-sensitivity.
* Against it: the dilution slope of §F1.4(a) accounts for only 20% of `chromatin_od`'s effect.

**The chromatin-density reading is therefore not established, and the log states it without a
hedge.** Separating the two needs a length-matched run, not a re-analysis.

---

## Tier 2 — real slips, conclusions survive

### F1.6 Two of the 35 cells are byte-identical draws, and the sign-flip null treats them as independent

`094.tiff` draws ann 2494 at seeds 2 **and** 4, and `201.tiff` draws ann 4457 at seeds 0 **and**
1. These are not merely correlated cells — every column is identical, and so is the delta
(−0.049383 twice; 0.0 twice). Under the actual randomization they are **one** observation and
must flip as a block; `signflip_spearman` flips them independently, which narrows the null and
is anti-conservative. n = 35 overstates the independent information (33), and there are 9
distinct null cells, not 10.

The pre-registration's defence — *"`Var(mean) = σ²/5` already prices it"* — is about estimating
a **mean** from draws with replacement. It says nothing about stratum exchangeability in a
permutation null. Different argument, and it does not cover this.

**Impact, measured, with duplicates collapsed to one cell each:**

| subset | published | duplicates collapsed |
|---|---|---|
| all cells | n=35, rho +0.5986, p 0.0012 | n=33, rho +0.5921, **p 0.0025** |
| informative only | n=25, rho +0.5591, p 0.0019 | n=24, rho +0.5732, **p 0.0016** |
| strictest (zero-deltas + seed 0 dropped) | n=17, rho +0.3619, p 0.0824 | n=16, rho +0.3661, **p 0.1020** |

Real violation, no conclusion change. It should be a stated limitation rather than a claimed
non-limitation.

### F1.7 Population mislabelling — one error, three places

Each is harmless alone; together they make the log's `n` column untrustworthy.

1. **Crossover table.** The `< 0.50` row (n=14) is the full pre-registered population; the
   `>= 0.50` (n=11) and `>= 0.58` (n=7) rows silently **drop the 10 null cells**, which all sit
   at `area_ratio = 1.0` with delta 0. On the pre-registered 35-cell population the rows are
   n=21, mean **+0.0025**, 4 neg / 5 pos / **12** zero; and n=17, mean **+0.0040**, 1 neg /
   4 pos / 12 zero. The conclusion (indistinguishable above the crossover) is unchanged.
2. **Cost side.** *"largest_cc's list is longer in every informative cell … mean +2,751"* — the
   mean over the 25 informative cells is **+3,851**. +2,751 is the mean over all 35 including
   the ten zero-delta nulls (3851.4 × 25/35 = 2751.0 exactly). Population and label disagree.
3. **The likely seed of both.** `f1_analysis.py:83` prints
   `f"n={len(w)} cells ({zero} null)"` where `zero = (w.delta == 0).sum()`. It labels
   **zero-delta** cells as **null** cells — 12 against the actual 10. Rename to `zero-delta`.

---

## Tier 3 — worth noting

### F1.8 The unstable sort is present but cannot reach a published number

`f1_seed_sweep.py:199` — `ps = pool.sort_values('score', ascending=False)` — is pandas' default
quicksort, exactly what the 2026-09-04 audit's Finding D flagged in the notebook, with
`kind='mergesort'` as the one-word fix. `compare._rank` uses mergesort and documents why.

**But it is inert here.** `ps` feeds only `ev.bucket_detections` for `n_dup_fp`, a context
diagnostic. Every arm frame is built from `pool`, and `compare.evaluate_arms` re-ranks it with
`_rank`/mergesort before scoring. `largest_tie_block` is 2-3, so at most one adjacent pair could
swap in `n_dup_fp` anyway. Hygiene fix, not a correction.

Worth stating separately: **GATE 1 structurally cannot catch a defect of this class**, because
it validates byte-exactness against reference CSVs produced by code carrying the same defect.

### F1.9 What "largest CC" actually varies — size only

`largest_cc_box` returns a bbox, but only `_odd_local(max(y1-y0, x1-x0))` is used; the template
is then a **click-centred** crop of the 73 px patch via `build_augmentations`. An off-centre
nucleus yields a template centred on the click, not on the component. The notebook is explicit
about this (cell 18: *"centred on the click, never recentred to the component's own
centroid"*), so it is a documented shared design choice, not an F1 divergence — but it bounds
what the result means: **F1 answers "does shrinking the template help?", not "does fitting the
template to the nucleus help?"** The "~36 px threshold" headline should not be read as the
latter.

Relatedly, the connected component is found on `gray_inv` while the matching runs on
`hematoxylin_od` — the size is derived from a different image than the one searched. Documented
in `largest_cc_box`'s docstring, consistent between the notebook and F1, and defensible (size is
a morphological property), but it is an assumption nothing has tested.

### F1.10 The primary axis is outside D5's standing constraint

`DECISIONS.md` D5 requires any run ranking by a chromatin statistic — *"as primary axis, as a
compared axis, or as a re-measurement of this decision"* — to use the 14 ROIs in
`images/extra_valid/`. F1 ran on the 7-ROI draw and made `chromatin_od` primary.

`f1_analysis.py:14-23` addresses this head-on and the reasoning is right: F1 was complete and
published before D5 existed, and re-keying a finished pre-registered result to an axis chosen
after seeing the data is exactly what pre-registration prevents. **No change wanted.** Noted
only so that a forward use of F1's `chromatin_od` numbers carries the caveat — and F1.5 above
is the concrete reason it matters.

---

## Re-checked and confirmed correct

Re-derived from the CSVs during this audit; no need to re-verify:

* **The pre-registered decision-rule inputs.** Primary rho = +0.5986, sign-flip p = **0.0012**;
  seed-0-excluded rho = +0.5569, p = **0.0119**. Both reproduce from `f1_analysis.py` exactly as
  the headline reports them. Verdict row 1 of the decision table is correct.
* **GATE arithmetic.** 11 gated columns × 96 (arm, z, budget) cells × 7 domains × 2 arms =
  **14,784**, matching the log. `len(Z_LEVELS)=6 × len(AXES)=2 × len(BUDGETS)=8 = 96`. The gate
  is keyed on (arm, z, **budget**), which is the right key — collapsing budgets would let a
  ranking divergence that preserves `read_50`/`read_95` pass silently.
* **The paired location shift.** mean −0.0129, SD 0.0275, **t = −2.772 (p = 0.0090)**,
  sign-flip **p = 0.0083**. All three as published (the log's "t = −2.77" is a rounding of
  −2.772). Only the Wilcoxon differs — see F1.2.
* **The effect-size table.** All seven per-domain within-arm `base51` seed SDs and mean deltas
  reproduce, including the two ratios above 1: human melanoma 0.0066 / −0.0261 = 3.92×, human
  breast cancer 0.0110 / −0.0346 = 3.13×. The claim "the effect is smaller than seed noise is
  false in 2 of 7 domains" is correct.
* **The crossover direction.** `< 0.50`: n=14, mean −0.0360, 13 neg / 1 pos. Exact.
* **The `tm_score` secondary.** mean delta −0.0153, 14 neg / 7 pos; cross-axis per-cell delta
  correlation 0.352 (Pearson). Exact.
* **`read_95` censoring by z.** Zero censoring at z ≤ 2.0 in both arms; 6 vs 2 at z = 2.5;
  12 vs 8 at z = 3.0; totals 36 vs 20. The log's reading — that the asymmetry is list length,
  not ranking quality — is right.
* **Cost side.** `read_95` deeper (worse) for `largest_cc` in **20 of 25** informative cells at
  z = 1.0, median **+355**. Exact.
* **Null-cell bookkeeping.** 10 null cells; per domain 201 → 1 informative, 246/402 → 3,
  301/459 → 4, 094/548 → 5, total 25. Exact, and it is a genuine reason the per-domain sign
  test is unavailable.
* **Closure discipline.** `cp.Arm(axis, (lambda d=sub: d), ...)` binds `sub` as a default
  argument, so the late-binding bug that would give every z the last cut's frame is avoided.
* **Deep floor.** F1's `DEEP_FLOOR_Z = -1.5` equals the notebook's `min(Z_LEVELS) - 0.5`, so the
  pools are built to the same depth in both runs.
* **`OD_PAD = 25` is independent of `PAD`** and correct: `nan_rate` is 0.0 on every row and the
  `pool['od'].isna()` assertion never fires, so no candidate's 51 px chromatin window was
  fabricated or dropped when the match padding shrank below 25.

---

## Recommended actions, in order

1. **Re-issue the Robustness table** in `Research Logs/2026-09-04-f1-results.md` with sign-flip
   p-values (F1.1), and move its computation into `f1_analysis.py` (F1.2). Note that
   "nulls + seed 0 dropped" clears 0.05 under the pre-registered test.
2. **Correct the pre-registration's `recall_at_budget` rationale** (F1.3) and add F1.4's three
   checks as the evidence that the conclusion survives it. This is the single most useful edit:
   right conclusion, wrong stated reason.
3. **Hedge the axis-specificity claim** (F1.5) or drop it pending a length-matched run.
4. **Move the duplicate-draw argument** from "not a limitation" to "limitations", with the
   measured impact (F1.6).
5. **Fix the three population labels** (F1.7), including the one-word `f1_analysis.py` print.
6. `kind='mergesort'` in `f1_seed_sweep.py:199` (F1.8) — hygiene only.
7. **If one more run is affordable**, the run worth having is a **length-matched** arm
   comparison — equal `n_detections` or equal `coverage_frac` rather than equal `z` — on the 14
   ROIs of `images/extra_valid` per D5. It would close F1.3, F1.5 and the "is it apples to
   apples" question in one pass.
