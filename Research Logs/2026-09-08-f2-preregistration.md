# Pre-registration F2: where is the template-size threshold, when the dose is assigned rather than observed?

Date: 2026-09-08. **Written and committed before `f2_base_size_dose.py` exists.**
Follow-up to `Research Logs/2026-09-08-f1-largest-cc-vs-51px-audit.md`, which is committed at
`18dc53b` together with `DECISIONS.md` D5 and `verify_chromatin_ranker.py`, so every citation
below points at something a reader of this repository can find. F1–F5 are taken; this is F2 by
the audit's own numbering of the follow-up it recommends.

---

## 1. The question, in the terms it was asked

`Research Logs/2026-09-04-f1-results.md` closed with:

> **The useful result is a threshold, not a winner: do not tighten below ~36 px.**

and immediately conceded how it was located:

> The threshold (~36 px) is estimated from 25 informative cells with a crossover located by
> eye at the midpoint. Locating it properly wants either more seeds or a deliberate sweep of
> `base_size` at fixed seed -- the latter is cheaper and removes the seed draw entirely as a
> source of dose variation.

F2 runs that sweep. The question is **where the threshold is**, on the 14-ROI set `DECISIONS.md`
D5 requires, with the dose **assigned** to every cell rather than read off whatever the
largest-CC rule happened to pick.

## 2. What is already known, so this run is not re-answering it

| established | where | consequence for this design |
|---|---|---|
| Tightening below ~36 px costs recall; at or above it the arms are indistinguishable | F1 results | The *direction* is not under test. The **location** is. |
| The effect is real and not a candidate-volume artefact | Audit F1.4 — within-arm dilution slope explains ~20%; length-matching in both directions leaves it intact | Volume is a covariate to report, not the finding to re-litigate. |
| The largest-CC rule changes template **size** only — the crop stays click-centred | `tm_threshold_axis_sweep_largest_cc.ipynb` cell 18; audit F1.9 | A `base_size` sweep is the *same intervention*, cleanly dosed. This is not a proxy. |
| `chromatin_od` is not the production ranker | D5; `verify_chromatin_ranker.py` | `tm_score` arbitrates (§3d). `chromatin_od` is reported, which is what triggers D5's 14-ROI requirement. |
| `recall_at_budget` is primary; `read_*` are diagnostics | D4 | §6.1 / §6.3. |
| Thresholds must be per-image, in robust-z units | D1, D3 | `med + z·MAD` per arm, as F1. |
| `TM_CCOEFF` is set at the call site, not in `FSConfig` | D1 amendment | `METHOD = cv2.TM_CCOEFF` explicitly. |
| Pool size must be reported beside every depth figure | D2 | §6.3, §6.5. |

### 2a. The measurement that reshaped this design, made before writing it

The audit recommended a **length-matched** run, because F1's two arms were matched on `z` and
not on candidate count (informative-cell mean Δn = +3,851). Re-derived from
`results/f1_seed_sweep.csv` before designing F2:

* On `tm_score`, `recall_at_budget` is **exactly invariant** to the threshold — **0 of 560**
  (cell, budget) groups vary across the six z levels at K ≤ 2000, and `read_50`/`read_95`
  never differ in **70 of 70** cells. Top-K-by-score of `{score ≥ cut}` *is* top-K of the pool,
  so thresholding cannot move it.
* On `chromatin_od`, **416 of 560** groups vary, because `score ≥ cut` gates membership while
  `od` decides order.

**So length-matching is vacuous on the ranker this repository ships**, and F1's `tm_score`
comparison was already matched on volume by construction. Matching bites only on
`chromatin_od`, where it is the discriminator for audit finding F1.5. It is therefore carried
as a secondary (§3e, §6.4) rather than as the design's purpose, and the purpose becomes the
dose sweep.

The K ≤ 2000 bound in that sentence is an artefact worth naming: invariance fails in 25 of 560
groups at K = 5000, and in every one of them `budget_delivered` has fallen below the budget
(minimum 2,379). The property is "constant across z **wherever the budget is delivered**", and
§7.3 gates it in that form rather than in the form that passes trivially.

### 2b. The dose the shipped rule actually produces, measured on all 14 ROIs

`largest_cc_box` + `_odd_local(max side)` over each ROI's border-filtered agreement pool
(14 ROIs, 993 candidates, run 2026-09-08 before the grid was chosen):

| statistic | value |
|---|---|
| range across all pools | 19 – 51 px |
| median of the 14 per-ROI medians | **37 px** |
| per-ROI medians | 27, 29, 33, 33, 33, 35, 35, 39, 41, 43, 47, 48, 51, 51 |
| fraction landing exactly at 51 (no tightening) | 0.20 mean per ROI, 0.0 – 0.6 |

The grid in §3a is chosen against this, not against taste.

---

## 3. Design

### 3a. One factor, seven levels

| arm | `base_size` | area ratio vs 51 | `PAD` | why this level |
|---|---:|---:|---:|---|
| `b19` | 19 | 0.139 | 9 | the low tail F1 measured as costly; smaller than anything the largest-CC notebook ran |
| `b27` | 27 | 0.280 | 13 | mid-low |
| `b33` | 33 | 0.419 | 16 | **below** F1's estimated crossover |
| `b37` | 37 | 0.526 | 18 | **above** it, and the median dose the shipped rule picks (§2b) |
| `b43` | 43 | 0.717 | 21 | mild tightening |
| `b51` | 51 | 1.000 | 25 | **reference level** — `tm.BASE_SIZE`, no tightening |
| `lcc` | rule's own | varies | varies | observational: what the shipped rule delivers |

Six assigned levels bracket the threshold with 33 and 37 straddling it. `lcc` is not a dose
level and never enters the primary test; it is the link back to F1 and to the shipped pipeline,
read as a residual against the assigned curve (§6.6).

**What the ladder can show:** the sign and size of Δ(recall@250) at each dose, and therefore an
interval containing the threshold. **What it cannot:** resolve the threshold finer than one grid
step (33 → 37). A point estimate is not claimed anywhere in this document.

**Held constant in every arm:** channel `hematoxylin_od`; `METHOD = cv2.TM_CCOEFF` passed to
`template_match.fused_response` at the call site (D1 amendment); `CFG.patch_size = 73`;
`peak_min_distance = 7`; `self_hit_radius = 5.0`; NMS radius 5.0 µm decoupled from the match
radius; `MATCH_RADIUS_UM = ev.MIDOG_RADIUS_UM`; `OD_PAD = 25`; `DEEP_FLOOR_Z = -1.5`;
`MAX_PEAKS = 2_000_000`; `Z_LEVELS = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)`; `cp.BUDGETS`; and the seed
click, which is drawn once per cell and handed to all seven arms.

**Not constant, and it cannot be:** `PAD = (base_size - 1) // 2`, which is 9 at 19 px against 25
at 51 px. Every arm reaches `valid.all()`, so the *property* is held fixed and the *amount* is
not; within ~25 px of the ROI edge the arms score against `BORDER_REPLICATE` fabrication of
different widths. Same disclosed limitation as F1, and it is unavoidable without changing what
the intervention is.

### 3b. The pairing — what makes this a controlled contrast

Every cell runs all seven arms against the **same drawn click** and the **same evaluation GT**
(that click's annotation removed). The dose contrast is therefore entirely within-cell, and the
seed draw is not a source of dose variation. This is the whole point of the design and it is
what F1 could not do.

**Implementation trap, named so it cannot happen.** The dose loop must run `b51` **first** and
record its `n_detections` at z = 1.0, because §3e's `matched_n` rule is defined from it. A loop
that computed `matched_n` before the reference existed would silently fall back to the arm's own
count and turn the secondary into an identity.

### 3c. Sample

14 ROIs (`images/extra_valid/`) × 5 seeds = **70 cells**, 490 arm-runs.
**The ROI is the design cluster, not the cell** — see §8.

Seeds come from `tm_variant_sweep.draw_seeds`, the house function F5 uses:
`np.random.default_rng([seed_index, image_id])` drawing **without replacement** from
`agreement_pool` → `border_filter(pool, 36, roi_shape)`. F1 drew *with* replacement and produced
two byte-identical cells that its sign-flip null then treated as independent strata (audit F1.6).
That failure mode is removed by construction, and it would have been worse here: 5 of the 14 ROIs
have seed pools of 7–13.

**Verified 2026-09-08, before this document was written, because §7.1 depends on it:**
`draw_seeds`' seed 0 returns the same `ann_id` as F1's `draw_seed_with_retry` on all 7 shared
ROIs — 094→2512, 201→4457, 246→6548, 301→14969, 402→20254, 459→22441, 548→25804 — and all 14
ROIs yield 5 distinct seeds. The code paths differ (F1 wraps the draw in a `largest_cc_box`
retry loop) but F1 recorded 0 retries in all 35 cells, so the realized seed is identical. Seeds
1–4 legitimately differ and are not gated.

### 3d. The primary axis, declared once

**`PRIMARY_AXIS = 'tm_score'`.** Per D5: chromatin density is not the production ranker, and no
experiment may declare `chromatin_od` primary without measuring on its own data that it beats
`tm_score`.

`chromatin_od` is computed and reported on every arm. That reporting is itself what triggers
D5's standing constraint — *"as primary axis, **as a compared axis**, or as a re-measurement"* —
which is why F2 runs on `images/extra_valid/` and not on the 7-ROI densest-per-domain draw.

**The head-to-head is pre-registered as an unconditional secondary (§6.4), with its consequence
fixed now.** F5 results §2.1 measured `chromatin_od − tm_score` = +0.0869, CI [+0.0419,
+0.1332], p = 0.0020, 12/14 ROIs on these same 14 ROIs — but ran that test after seeing both
axes' answers, which is why F5 §2.6 declined to re-declare the axis. F2 fixes the rule before
the data exists: **if S3 replicates that result, both axes are reported with neither labelled
secondary, and the Holm family stays at 5 (the five dose tests on `tm_score`).** The arbiter
does not move. Writing this down is what stops it becoming axis-shopping.

Per F5's Correction 2 §2.8, the primary axis appears in five coupled places and revision 4
changed two of them. In this document it is §3d, §6.1, §8 and §12, plus `PRIMARY_AXIS` in
`f2_base_size_dose.py` and `f2_analysis.py`. §7.10 gates that they agree.

### 3e. Two threshold rules, both computed in the run

| rule | definition | `z` | what it answers |
|---|---|---|---|
| `matched_z` | `score ≥ med + z·MAD`, six z levels | as recorded | continuity with F1; the shipped behaviour |
| `matched_n` | top-N by score, **N = this cell's `b51` count at z = 1.0** | `NaN` | §6.4 / audit F1.5 |

`matched_n` costs no template matching — it is a second slice of a pool already built — so it is
a first-class output rather than an offline re-analysis. Pools are *also* persisted (§5) so any
further rule can be evaluated later without re-running, and §7.2 gates that round trip.

**Matched-coverage is deliberately not a rule.** `coverage_frac` at z = 1.0 runs 0.90–0.999, so
matching on a near-saturated quantity is ill-conditioned and the matched value would be
dominated by rounding. Coverage is *reported* beside every depth figure instead, which is what
D2 requires anyway.

---

## 4. Predictions, registered before the run

Stated so they can fail.

* **P1 (the primary).** Δ(recall@250) on `tm_score` is ≤ 0 at every dose and **monotone
  increasing** in dose, crossing 0 between 33 and 43 px. *Confidence: moderate on the sign and
  monotonicity, low on the location.* The live question is the location; F1 estimated it on 7
  ROIs from an observed dose, and this run's grid could put it below 33 or above 43.
* **P2 (the F1.5 discriminator).** The `chromatin_od` dose curve is **steeper** than
  `tm_score`'s under `matched_z`, and **converges toward it** under `matched_n`. *Confidence:
  moderate.* If it does **not** converge, F1's chromatin-specific reading is supported and the
  audit's competing explanation is wrong — that is a real outcome, not a failure.
* **P3 (volume).** `n_detections` rises monotonically as dose falls, in every cell.
  *Confidence: high* — it held in 25 of 25 informative F1 cells.
* **P4 (censoring).** No `read_95` censoring at z = 1.0 in any arm. *Confidence: high* — F1 had
  0 of 70. If it appears, §6.3's depth comparison is not clean and will be reported as censored,
  never dropped.

**What would make F2 report "there is no usable threshold"** — bound to §8's arbiter, not to a
descriptive: if every dose's Holm-adjusted test fails to reject and the per-dose CIs all contain
0, including `b19`, then template size does not move recall@250 on `tm_score` at this sample
size, F1's threshold does not replicate on the 14-ROI set, and the recommendation reverts to
"template size is not a decision-relevant lever". That is a clean negative and will be reported
as one.

---

## 5. Deliverables

* `f2_base_size_dose.py` — the runner. Committed before it is executed.
* `f2_analysis.py` — committed **before results exist**. Every number in the results log must
  come from it; F1's Robustness table did not come from `f1_analysis.py` and that is audit
  finding F1.2.
* `results/f2_base_size_dose.csv` — one row per (ROI, seed, arm, rule, z, budget).
* `results/f2_base_size_dose_cells.csv` — one row per (ROI, seed), with `run_complete`,
  `n_cells_expected`, per-dose timings and the `lcc` realized size.
* `results/f2_base_size_dose_verification.csv` — every gate record. **Read it before any table.**
* `.cache_tail/f2_pools/*.parquet` — the persisted pools, gitignored.
* `Research Logs/2026-09-08-f2-results.md` — after the run.

---

## 6. Metrics, and which one arbitrates

### 6.1 Primary — recall@K (D4)

Budgets `cp.BUDGETS = (25, 50, 100, 250, 500, 1000, 2000, 5000)`. **K = 250 arbitrates**, the
same K F1 used and for the same measured reason: at z = 1.0 the reference arm's recall@250 is
unsaturated in every domain, while K ≤ 100 gives discretised near-zero deltas and K = 5000
begins to run out of list.

**Reporting frame (D4), descriptive and without a p-value:** recall@250 per tumour domain, worst
ROI and worst seed, per dose. D4's "worst ROI, worst click, per domain" is a *reporting* rule;
the quantity that arbitrates is §8's and is a different statistic. F5 revision 2 conflated the
two and had to correct it.

**Two guards, mandatory before any cell is read:** `coverage_frac` printed adjacent to every
recall figure (§6.5), and `n_pool` / `n_detections` printed adjacent to every depth figure (D2).

### 6.2 The threshold, reported as an interval

**Pre-committed definition:** the threshold is reported as the **interval between the largest
dose whose Holm-adjusted test rejects 0 and the smallest dose whose test does not**, together
with the full per-dose CI table. No point estimate, and no equivalence test — failure to reject
is not equivalence, and inventing a margin here would be machinery the question does not need.
If the rejecting and non-rejecting doses are not adjacent, that non-monotonicity is reported as
the finding it is.

### 6.3 Reading depth — diagnostic, and a correction to why F1 demoted it

`read_50 … read_100` are reported per dose with `n_pool`, `n_detections` and `coverage_frac`
adjacent, and with censored cells counted, never dropped. They **do not arbitrate** — D4 retires
them, having found that "the entire `read_50` result reproduces with no click at all".

That is the correct reason. F1's stated reason was different and was wrong: it demoted `read_95`
because a delta "would re-report the candidate-volume difference under a new name". On
`tm_score` that cannot happen — `_read_depths` is `np.searchsorted` into a prefix of the same
`tp_cum` array, so a threshold change leaves every finite `read_*` identical or turns it NaN
(audit F1.3, and §2a's 70-of-70 measurement). Between *arms* the pools genuinely differ, but a
longer list does not mechanically inflate a finite `read_95`; it only removes censoring, and F1
had none at z = 1.0. So depth is a fair between-arm comparison here, reported as a diagnostic
under D4 rather than suppressed under a mistaken confound.

### 6.4 The two secondaries that carry the audit's open questions

* **S1 / S2 — the `chromatin_od` dose curve under `matched_z` and under `matched_n`.** Their
  difference is the F1.5 discriminator: if the curve flattens toward `tm_score`'s under
  `matched_n`, `chromatin_od`'s steeper dose-response was list length; if it does not, it was
  chromatin sensitivity.
* **S3 — the axis head-to-head.** Δ(recall@250) `chromatin_od − tm_score` on the `b51` arm at
  z = 1.0, 14 ROI clusters, same estimator as §8. Reported whatever it shows, with the
  consequence fixed in §3d.

### 6.5 Carried diagnostics

`n_pool`, `n_detections`, `coverage_frac`, `largest_tie_block`, `nan_rate`, `full_list_recall`,
`n_lookalike_in_list`, `floor_limited`, `z_cut`, `pad_px`, `map_median`, `mad_scale`,
`n_gt_mitotic`, `n_retries`, `tightened_size` (the `lcc` arm), per-arm wall-clock.

`full_list_recall` is **never a headline**: at z = 1.0 coverage is 0.90–0.999, so
`full_list_recall == 1.0` is close to geometry (the 2026-09-03 audit's finding 1). It is always
printed with `coverage_frac`.

### 6.6 The `lcc` residual — a reported diagnostic, not a test

Each `lcc` cell's Δ against the assigned-dose curve linearly interpolated at that cell's own
realized size. Pooled over 70 cells it should centre on 0; a systematic offset would mean the
rule does something beyond setting a size, contradicting the notebook's own account (audit
F1.9). Its x differs per cell, so there is no single curve point to test it against and none is
claimed.

---

## 7. Verifications and gates — asserted in the run, not claimed in prose

1. **Seed-0 reproduction.** `b51` and `lcc` at `seed_index == 0` must reproduce
   `results/f1_seed_sweep.csv` on the 7 overlapping ROIs, across 11 columns × 96 (arm, z,
   budget) cells, under `matched_z`. Checked **inside the ROI loop** so a divergence aborts on
   the first ROI rather than after hours. It is the only check that can catch a silent
   divergence in the whole pipeline rather than in one function.
2. **Pool round-trip.** For one cell per ROI, re-evaluating the persisted parquet under
   `matched_z` reproduces that cell's in-run rows exactly. Without it every future offline
   re-analysis of those pools is unverified.
3. **The prefix property, in the form that can fail.** For every (cell, dose, budget) where
   `budget_delivered == budget`, `recall_at_budget` on `tm_score` is constant across z. Not
   "K ≤ 2000", which passes trivially. §2a's whole argument rests on this.
4. **`matched_n` feasibility.** Every dose's deep pool ≥ N, asserted per cell; otherwise the rule
   silently returns a shorter list and the "matched" comparison is not matched.
5. **`matched_n` is actually matched.** Every `matched_n` arm's `n_detections` equals N exactly.
6. `invariants.check_distinct_seeds` per ROI — guaranteed by `draw_seeds`, asserted anyway.
7. `invariants.check_no_cap` on the **pre-NMS** pool length, not the post-NMS length
   `compare.py:172` checks: against a 2,000,000 cap the latter records a pass it could not have
   withheld (F5 §7.6).
8. `valid.all()` per arm, and `pool['od'].isna().sum() == 0` after `OD_PAD` — the NaN class the
   v2 padding bug introduced.
9. One-match-many-z shortcut re-verified at `seed_index == 0` per (ROI, dose): thresholding the
   deep pool at `med + z·MAD` equals a fresh `extract_peaks` at that floor.
10. **Protocol coherence.** `PRIMARY_AXIS`, `Z_PRIMARY`, `K_PRIMARY` and the Holm family size
    are asserted equal between `f2_base_size_dose.py` and `f2_analysis.py`, and against the
    values written in §3d / §6.1 / §8. F5's revision 4 changed the axis in 2 of 5 places and
    produced a protocol with no coherent reading.
11. `compare.assert_floor_not_limiting` over the finished frame.
12. **Failure protocol.** Output carries a `.partial` suffix and `run_complete=False` until the
    run finishes; only a completed run earns the un-suffixed name. Dumped once per ROI, so a
    crash at ROI 13 of 14 leaves twelve ROIs of usable, clearly-labelled partial data. F5 wrote
    everything at the end and would have lost 37 minutes to a late assertion.

Every record lands in `results/f2_base_size_dose_verification.csv`.

---

## 8. Statistics, fixed before the data exists

> **Primary estimand.** For each dose *d* ∈ {19, 27, 33, 37, 43}: form the paired difference
> `recall@250(d) − recall@250(b51)` within each of the 70 cells, at **z = 1.0**, axis
> **`tm_score`**, rule **`matched_z`**; average the 5 seeds within each ROI; treat the resulting
> **14 ROI means as the clustered units**.

* **Test:** exact sign-flip over all 2¹⁴ = 16,384 sign assignments of the 14 ROI means,
  two-sided. Exact, so no permutation-count arbitrariness and no Monte-Carlo p floor — the
  defect that let F1's results log publish 0.0001 and 0.0012 for the same test (audit F1.1).
* **Interval:** cluster bootstrap over ROIs, 10,000 resamples, `seed=20260908`, percentile CI.
* **Multiplicity: Holm across the 5 dose tests.** The family is 5, not 10, because §3d declares
  `tm_score` primary. "Holm across 10" and "the arbiter is `tm_score`" cannot both stand.
* **Effect size** is reported in recall points and against **each ROI's own within-arm seed SD**
  under this exact configuration — not against `tm_ccoeff_headtohead_seed_variance.csv`, whose
  SDs come from a different pipeline.

**Why the ROI and not the cell.** 70 cells are clustered in 14 ROIs and 7 domains; cells within
an ROI share the tissue, the scanner and the annotation pool. Treating them as 70 independent
units would be anti-conservative. Averaging within the ROI first is the same estimator F5's
arbiter uses.

**Why not F1's Spearman-over-dose.** F1 had to correlate the delta against an *observed* dose
because it had no control over which dose each cell received. That forced sign-flip permutation
over cells, an argument about null-cell anchoring, and a duplicate-draw correction — and it
still could not locate the threshold. With the dose assigned, the direct measurement is a mean
per dose with an interval. No correlation, no permutation over an observed covariate, no
machinery to defend.

**Secondaries** (S1, S2, S3 of §6.4, plus P3's monotonicity and §6.6's residual) use the same
estimator and are reported **unconditionally, whatever they show**, without multiplicity
adjustment and labelled as unadjusted. They are pre-specified because they answer different
questions off the same paired matrix — not so the better-looking one can be chosen afterwards.

**Dedup key: `(file_name, seed_index, arm, rule, z)`.** D4's data trap — the CSV emits one row
per budget, so every arm appears 8 times with identical `read_*` — multiplies any row-level
count or p-value by 8. F5 upgraded D4's three-tuple to include `z`; F2 adds `rule`. Valid
because arm names embed axis and dose (`tm_score@b37`).

---

## 9. Circularity disclosure

P1's 33–43 px bracket comes from F1's **completed** data, which I have read: its crossover sits
at `area_ratio ≈ 0.50`, i.e. ~36 px. Two further contaminations, stated rather than denied:

1. **7 of the 14 ROIs are F1's ROIs**, and seed 0 is the identical draw (§3c). So 7 of the 70
   cells have counterparts whose F1 numbers I have seen.
2. The audit that motivated F2 re-derived F1's numbers in detail, including the per-dose split
   table.

**How it is handled.** The primary is a **per-dose test with a pre-committed reporting rule
(§6.2)**, not a test of the bracket — every dose's result is reported whatever it is, so a
bracket in the wrong place is visible rather than absorbable. Two decontaminated robustness
readings are pre-committed and will be reported beside the primary: **(a)** the primary
restricted to the **7 ROIs F1 never ran** (013, 233, 245, 300, 403, 460, 529), and **(b)** the
primary with `seed_index == 0` excluded. Neither is decisive — at 7 clusters (a) is
underpowered by construction — but a bracket that survives in neither is a bracket generated by
the cells that suggested it, and that will be said.

Five of the seven dose levels (19, 27, 33, 43, and the `b51` reference at four of the seven
shared ROIs) were never run by F1 on any ROI, so most of the grid is new even on the shared
cells.

---

## 10. Limitations

### What is *not* a limitation

**Reusing F1's ROIs and seed 0.** It looks like non-independence and is not, for the primary
question: the dose contrast is **within-cell**, and F1 ran only two of the seven levels. What the
reuse buys is §7.1 — the only gate that can catch a silent divergence across the whole pipeline
— and giving it up to obtain a nominally fresh sample would trade a real check for a cosmetic
one. §9 handles the part that genuinely is contamination.

**Seed pools of 7–13 on 5 of the 14 ROIs.** The 5 draws are a large fraction of the finite pool
there, which *narrows* the estimand rather than biasing it: on those ROIs the run nearly
enumerates the population of clicks a pathologist could make, which is the population the
product cares about. `MANIFEST.md`'s own bar is seed pool ≥ 5, and 201.tiff's pool of 7 is the
binding case — sampling without replacement is what makes 5 draws well-defined there at all.

### Real limitations

* 14 ROIs, 7 domains, 2 per domain. Not exchangeable; every pooled statement is about this sample.
* `PAD` varies with dose (9 at 19 px, 25 at 51 px), so the arms are not comparable within ~25 px
  of the ROI edge. Unavoidable without changing the intervention (§3a).
* The grid cannot resolve the threshold finer than 33 → 37.
* `n_gt_mitotic` is 18, 18 and 20 on 013/233/529, so recall@250 there moves in steps of ~0.056
  and those ROI means are coarse.
* `lcc`'s realized size is not randomised — it correlates with whatever makes a nucleus large in
  that ROI, so §6.6's residual is descriptive only.
* D5's own evidence point 3 (`chromatin_od` and `tm_score` indistinguishable, +0.032, p = 0.36)
  was measured on 7 ROIs at one seed, and F5 §2.1 measured the opposite on these 14. §3d holds
  the arbiter on `tm_score` by D5's letter while pre-registering S3 to measure it here; a reader
  should know the constraint and the evidence behind it are in tension.

---

## 11. Cost

**~1.9 h estimated**, from F1's own per-template-size wall-clock (19 px ≈ 28 s, 27 ≈ 15, 33 ≈ 6,
37 ≈ 10, 43/51 ≈ 13 → ~98 s per cell × 70 cells). **Budget 3–4 h.** F5's cost section predicted
9–12 min for a run that took 37, and this estimate is built from timings whose spread within one
template size is 3× — treat 1.9 h as a lower bound, not a forecast.

Machine: CPU-only Intel i7-8750H, 6 threads, no GPU path. Interpreter `~/anaconda3/bin/python`
(numpy 1.26.4, cv2 4.8.1, pyarrow 16.1.0); the system `python3` carries numpy 2.4.6 and cannot
`import cv2`, so `midog_utils` will not import there.

Storage: ~12.7 M pool rows, ~254 MB as parquet (int32 `cx`/`cy`, float32 `score`, float64 `od` -- float32 would be lossy and break section 7.2) against
~635 MB as CSV. Parquet is a deviation — nothing in the active codebase writes it — taken for
the 3× size and recorded here rather than in a commit message.

---

## 12. Decision rule, fixed now

Read on the §8 arbiter: Holm-adjusted, exact sign-flip, `tm_score`, recall@250, z = 1.0,
`matched_z`, 14 ROI clusters.

| outcome | conclusion |
|---|---|
| ≥ 1 dose rejects, and the rejecting doses are the **smaller** ones, contiguously | **P1 supported.** The threshold lies in the interval of §6.2 and F1's finding replicates on the D5 set with the dose assigned. This is the deliverable. |
| ≥ 1 dose rejects, but not contiguously in dose order | **Non-monotone.** Reported as such; the threshold framing does not survive and the per-dose table is the result. F1's monotone reading would be sample-specific. |
| no dose rejects after Holm, and every CI contains 0 | **Clean negative.** Template size is not a decision-relevant lever at this sample size; F1's threshold does not replicate on 14 ROIs. Reported as a result, not as a failure. |
| `b19` rejects but the CI at every larger dose contains 0 | **Threshold below 27 px.** The practical recommendation weakens to "do not tighten below ~27 px" and F1's ~36 px was too conservative. |

Secondary S1/S2 decide audit finding F1.5 independently of the above, and are reported whichever
way the primary falls.

**A null is a real result and closes this question honestly. I am not entitled to go looking for
a different (axis, z, K, rule) afterwards.** Any post-hoc reading is labelled exploratory in the
results log, as §8 requires for the secondaries.
