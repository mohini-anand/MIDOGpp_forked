# Independent audit: F1, and whether bbox tightening beats the fixed 51 px template

Date: 2026-09-08. Second, independent pass over the same material as
`Research Logs/2026-09-08-f1-largest-cc-vs-51px-audit.md`, written without editing it.

Scope: `f1_seed_sweep.py`, `f1_analysis.py`, `Research Logs/2026-09-04-f1-preregistration.md`,
`Research Logs/2026-09-04-f1-results.md`, and the two result CSVs. Everything below is
re-derived from `results/f1_seed_sweep.csv`, `results/f1_seed_sweep_cells.csv`,
`results/tm_ccoeff_threshold_axis_sweep_v2.csv`,
`results/tm_ccoeff_threshold_axis_sweep_largest_cc.csv` and the seven ROI images. No sweep was
re-run. The one new measurement that touches pixels is §N4, which re-derives all 35 seed boxes.

---

## Part 0 — what reproduces

Re-run independently, not read off the log:

| check | result |
|---|---|
| **GATE 1** re-executed against both reference CSVs, 11 columns × 96 (arm, z, budget) × 7 ROIs × 2 arms | **14,784 values compared, 0 divergences** |
| primary sign-flip test, all 35 cells | rho **+0.5986**, p **0.0012** — exact |
| seed-0-excluded, 28 cells | rho **+0.5569**, p **0.0119** — exact |
| informative-only, 25 cells | rho **+0.5591**, p **0.0019** — exact (matches the prior audit's F1.1 correction) |
| `recall_at_budget == tp_at_budget / n_gt_mitotic` | true on every row |
| budget actually delivered at K=250 | 70/70 rows; `n_detections` 12,807–31,198, so `n >> K` holds |
| `nan_rate`, `largest_tie_block` | 0.0 everywhere; tie blocks ≤ 3 |
| F1's `tightened_size` vs a fresh `largest_cc_box` on all 35 seeds | **35/35 identical** |
| null cells | 10, delta exactly 0 in all of them (GATE 2 holds) |

**The engineering is sound and the published statistics are the statistics the committed code
computes.** Every defect below is in the *inference* and in *what the numbers are taken to
mean*, not in the run.

---

## Part 1 — new findings

### N1 — [Tier 1] The pre-registered null assumes 35 independent observations; the design has 7

`signflip_spearman` flips each cell's delta independently. The pre-registration defends this as
*"each cell is its own stratum, so domain clustering is handled by construction rather than by a
modelling assumption."* That is the wrong direction: stratifying by cell does not absorb
clustering, it **assumes it away**. The 35 cells are 5 seeds drawn on each of 7 ROIs; cells
sharing an ROI share the image, the GT, the response-map statistics and the candidate pool.
Independent flips across correlated units narrow the null and the test is anti-conservative.

Flipping at the **ROI** — the exchangeable unit `DECISIONS.md` D5 and F5 §8 both declare — gives:

| test | published (cell-level flip) | ROI-level flip |
|---|---:|---:|
| primary, all 35 cells | p **0.0012** | p **0.0669** |
| seed-0-excluded, 28 cells (**decisive** per the pre-registration) | p **0.0119** | p **0.1009** |
| informative-only, 25 cells | p 0.0019 | p **0.0333** |

The pre-registered decision rule requires the primary **and** the seed-0-excluded version to
clear 0.05. Under a cluster-robust null **neither does**, so the run lands in **row 3** of its own
decision table — *"no resolvable dose-response at n=35"* — not row 1.

And the same correction applies to the location shift the results log reports as
*"resolvable: mean −0.0129, SD 0.0275, t = −2.77, sign-flip p = 0.0083, Wilcoxon p = 0.0101"*:

| axis | population | mean Δ recall@250 | ROI-clustered 95 % CI | ROI sign-flip p | ROIs negative |
|---|---|---:|---|---:|---|
| `chromatin_od` | all 35 | **−0.0129** | [−0.0281, **+0.0016**] | 0.221 | 4/7 |
| `chromatin_od` | informative 25 | **−0.0180** | [−0.0386, **+0.0022**] | 0.219 | 4/7 |
| `tm_score` | all 35 | **−0.0153** | [−0.0450, **+0.0179**] | 0.498 | 4/7 |
| `tm_score` | informative 25 | **−0.0214** | [−0.0568, **+0.0275**] | 0.496 | 4/7 |

**Two honest caveats on this finding.** (a) It is conservative and coarse. Only 6 of the 7 ROIs
contribute any signal — `201.tiff`'s delta is exactly 0.0 in all five of its cells — so the
attainable p-floor is 1/2⁶ = 0.0156, and p = 0.0669 means roughly 8 of 64 sign patterns beat the
observed one. The correct reading is **"not resolvable with 7 units"**, not "refuted". (b) All
three cluster tests rest on the *same* 6 effective units, so the informative-25 subset clearing
0.05 while the other two do not is a statement about the null cells, not about power.

Which is itself worth recording: the pre-registration justified keeping the null cells partly
because *"they contribute no permutation entropy either way."* Under the cell-level flip they
help (p 0.0019 → 0.0012); under the ROI-level flip they **hurt** (0.0333 → 0.0669). They are not
inert under the correct null.

**The internal inconsistency is the real finding.** `DECISIONS.md` D5, written on 2026-09-04 from
*this same CSV*, clusters at the ROI when comparing the two axes and reports Δ = +0.032,
CI [−0.047, +0.112], **p = 0.36**. F1, on the same 70 cells, does not cluster when comparing the
two arms. One project, one dataset, two inference standards.

### N2 — [Tier 1] The dose is confounded with the domain, and the promised blocking never happened

`f1_seed_sweep.py`'s docstring and the pre-registration both specify a pooled dose-response
*"with domain as a blocking factor"*. `f1_analysis.py` never blocks: `report_primary` computes one
pooled Spearman over all cells.

It matters, because the treatment is close to a domain label. Sorted by mean dose:

| domain | tightened sizes | mean `area_ratio` | mean Δ |
|---|---|---:|---:|
| human melanoma | 29, 37, 33, 33, 33 | 0.421 | −0.0261 |
| canine cutaneous mast cell tumor | 41, 51, 19, 19, 27 | 0.441 | +0.0055 |
| human breast cancer | 25, 37, 35, 37, 35 | 0.447 | −0.0346 |
| human neuroendocrine tumor | 29, 51, 51, 31, 27 | 0.595 | −0.0442 |
| canine soft tissue sarcoma | 33, 47, 39, 51, 45 | 0.726 | −0.0031 |
| canine lymphosarcoma | 41, 51, 47, 51, 37 | 0.804 | +0.0122 |
| canine lung cancer | 51, 51, 51, 47, 51 | 0.970 | 0.0000 |

Kruskal–Wallis of `area_ratio` on domain: H = 14.68, **p = 0.023**. **45.7 %** of the dose variance
is between domains (**58.8 %** among the 25 informative cells). The three domains that never
tighten much are the three canine domains whose Δ is ≈ 0 or positive; the two most-tightened are
the two human domains with the largest negative Δ. A between-domain dose-response and "tightening
hurts human tumours and not canine ones" are the same fit.

**But the within-domain evidence is real, and points the same way.** Spearman rho computed
*inside* each domain:

| domain | doses | rho (all cells) | rho (informative) |
|---|---:|---:|---:|
| human neuroendocrine tumor | 4 | +0.895 | +0.500 |
| human melanoma | 3 | +0.894 | +0.894 |
| canine soft tissue sarcoma | 5 | +0.667 | +0.800 |
| human breast cancer | 3 | +0.649 | +0.649 |
| canine cutaneous mast cell tumor | 4 | +0.359 | +0.632 |
| canine lymphosarcoma | 4 | **−1.000** | **−1.000** |
| canine lung cancer | 2 | undefined (1 informative cell) | — |

**5 of 6 positive, mean +0.411.** This points the same way as the pooled result and is *not* a
domain effect — but it must carry the same deflation as everything else in this document: a
domain-level sign test on 5 of 6 is two-sided **p ≈ 0.22**, on the same 6 effective units. Domain-centring both variables and pooling gives rho = +0.465 on all 35
and **+0.370 on the 25 informative cells** — reported here as *descriptive*, since centring within
domain and then flipping at cell level is a hybrid null, not a clean test.

**Net:** the pooled rho of +0.599 is roughly half between-domain and cannot carry the weight the
log puts on it; the within-domain version is weaker but real. The defect is that only the pooled
number was ever computed or published.

### N3 — [Tier 1] "~36 px" is not an estimate, and the cell that most locates it is counted twice

The results log's own follow-up section says the threshold is *"estimated from 25 informative
cells with a crossover located by eye at the midpoint."* That is literal: the crossover table
splits `area_ratio` at **0.50**, and 51·√0.50 = **36.06 px**. 36 is where the midpoint of the
area-ratio scale falls, not where a fitted changepoint sits. There is no changepoint fit, no
confidence interval, and no pre-registered threshold estimator.

Scanning every attainable split of `tightened_size` over the 25 informative cells:

| cut (px) | n below | mean Δ below | n at/above | mean Δ at/above | separation |
|---:|---:|---:|---:|---:|---:|
| 33 | 8 | −0.0387 | 17 | −0.0083 | 0.0304 |
| 35 | 12 | −0.0338 | 13 | −0.0035 | 0.0302 |
| **37** | 14 | −0.0360 | 11 | +0.0048 | **0.0408** |
| 39 | 18 | −0.0288 | 7 | +0.0097 | 0.0385 |
| 41 | 19 | −0.0277 | 6 | +0.0125 | 0.0402 |
| 45 | 21 | −0.0223 | 4 | +0.0041 | 0.0264 |

Separation does peak at a cut of 37 px (i.e. "≤ 35 px is harmed"), so the published number is not
arbitrary — but 39 and 41 are within 0.0023 of it, and 33/35/45 within 0.015. On 6 effective
clusters **the whole 33–45 px band fits equally well.** "Do not tighten below ~36 px" should read
"below somewhere in the 33–45 px range, if the effect is a step at all."

**New, and it links to the prior audit's F1.6:** the boundary that most looks like the crossover
is 35 → 37 px, and it is carried almost entirely by `094.tiff` (35, 35, 37, 37 px → −0.0494,
−0.0494, −0.0247, −0.0123). Its two 35 px cells are **the same annotation drawn twice** — ann 2494
at seeds 2 and 4, byte-identical rows, delta −0.049383 in both. The duplicate draw that F1.6
flagged as an independence violation is also the pair that most sharply positions the threshold.

### N4 — [Tier 1] "Indistinguishable above 36 px" is a failure to reject presented as equivalence

| population | n | mean Δ | bootstrap 95 % CI |
|---|---:|---:|---|
| informative, < 36 px | 14 | −0.0360 | [−0.0493, −0.0226] |
| informative, ≥ 36 px | 11 | **+0.0048** | **[−0.0062, +0.0165]** |
| all cells ≥ 36 px (incl. the 10 nulls) | 21 | +0.0025 | [−0.0031, +0.0088] |

The ≥ 36 px interval cannot exclude a 1.5-point harm or a 1.5-point benefit. No equivalence
margin was pre-specified, and none of this is clustered — clustered it is wider still. The
log's *"the two arms are indistinguishable"* is true only in the sense that nothing was
distinguished.

### N5 — [Tier 1] The low end of the dose axis is built from seeds the shipped rule would refuse

F1's arm is `largest_cc_box`: `max(regions, key=area)` over the Otsu components of the 51 px crop.
The method the repo ships is `seed_selection.tighten_box_otsu`, which takes the component **under
the click** and returns `None` otherwise — explicitly refusing the largest-component fallback that
`BBOX_TUNING.md` calls *"a weak heuristic ... in a crowded crop it can grab an unrelated, larger
structure."*

Re-derived on all 35 F1 seeds:

* the click lands in a component that is **not** the largest: **0 of 35**. So wherever the shipped
  rule answers at all, the two rules are **byte-identical — 30 of 30**. F1 is *not* ablating a
  rejected heuristic; that framing would be wrong.
* the click pixel is **background**, so `tighten_box_otsu` returns `None` and the production
  pipeline would discard the seed and re-draw: **5 of 35**.

Those five are the problem:

| ROI | seed | tightened | Δ recall@250 |
|---|---:|---:|---:|
| 301.tiff | 2 | **19 px** | +0.0138 |
| 301.tiff | 3 | **19 px** | −0.0230 |
| 402.tiff | 4 | 27 px | **−0.0769** (the single most negative cell) |
| 201.tiff | 3 | 47 px | 0.0000 |
| 301.tiff | 1 | 51 px | 0.0000 (null) |

**Both 19 px cells — the two lowest doses in the entire experiment, the ones the pre-registration
singles out as *"smaller than anything the largest-CC notebook ever ran"* — are seeds the shipped
rule refuses.** Their "template size" is the extent of a component the pathologist did not click
on. `draw_seed_with_retry` would have retried and replaced all five had `check_fn` been the shipped
rule; with `largest_cc_box` as `check_fn`, retries were 0 in all 35 cells.

**This does not flip the conclusion** — dropping the five *raises* rho, to +0.7176 on the remaining
30 and +0.7103 on the 21 informative. It bounds what the dose axis means at its low end.

It also disposes of an apparent non-monotonicity rather than establishing one. Mean Δ by size
(informative cells) is −0.0046 at 19 px, −0.0370 at 25, −0.0408 at 27, −0.0574 at 29, −0.0673 at
31 — the smallest templates look *least* harmed, which contradicts H1 as stated. But n = 2 at
19 px, both from one domain, both refused seeds. The honest statement: **whether the curve is
monotone below 25 px is unmeasurable in this run, not measured-and-non-monotone.** Relatedly,
dropping `301.tiff` — the only domain supplying a 19 px dose — raises rho from +0.559 to **+0.827**,
so the monotone reading is strongest exactly when the lowest doses are removed.

### N6 — [Tier 1] There is no dose-response on the ranker the repo actually ships

`DECISIONS.md` D5 makes `tm_score` the production ranker and adds a standing constraint that no
run may make a chromatin statistic primary on this 7-ROI draw. F1's primary axis is
`chromatin_od`. `f1_analysis.py:14-23` gives the right reason for not retro-fitting it, and this
audit agrees with that call. But the forward consequence has to be stated:

All eight budgets in `compare.BUDGETS`, 25 informative cells, cell-level flips (the
anti-conservative null of N1):

| K | `chromatin_od` mean Δ | rho | p | `tm_score` mean Δ | rho | p | `tm_score` signs (−/+/0) |
|---:|---:|---:|---:|---:|---:|---:|---|
| 25 | −0.0063 | +0.189 | 0.264 | −0.0167 | +0.199 | 0.364 | 14/4/7 |
| 50 | −0.0095 | −0.119 | 0.536 | −0.0247 | +0.416 | 0.029 | 16/6/3 |
| 100 | −0.0088 | +0.562 | 0.010 | −0.0320 | +0.097 | 0.613 | 18/3/4 |
| **250** | **−0.0180** | **+0.559** | **0.002** | **−0.0214** | **+0.016** | **0.938** | 14/7/4 |
| 500 | −0.0181 | +0.472 | 0.011 | −0.0199 | −0.200 | 0.359 | 15/9/1 |
| 1000 | −0.0129 | +0.454 | 0.005 | **+0.0044** | −0.191 | 0.407 | 10/13/2 |
| 2000 | −0.0022 | +0.491 | 0.018 | **+0.0187** | **−0.458** | 0.038 | **3/16/6** |
| 5000 | −0.0099 | +0.643 | 0.001 | **+0.0078** | +0.014 | 0.944 | 5/15/5 |

Two things fall out.

**(a) On `tm_score` the dose-response is absent at every budget but one of eight — and at K=2000
its sign is the opposite of `chromatin_od`'s** (rho −0.458 against +0.491 on the same cells).
A monotone structure that inverts between the two rankers at the same operating point is not a
property of template size. **The "~36 px threshold" is a property of `chromatin_od` only** — the axis D5 demoted four days later, on evidence drawn
from this very CSV. The direction (tightening does not help) replicates on both axes; the
*threshold* does not replicate on the shipped one.

**(b) On `tm_score`, tightening flips from harmful to helpful past K ≈ 1000** — +0.0187 at K=2000,
where it wins in **16 of 25** informative cells against 3 losses. ROI-clustered this is not
significant (p = 0.217), and it is the same mechanism as the `read_95` censoring the results log
already attributes to list length: a smaller template yields a longer list, which reaches deeper.
It is worth stating because it is the one budget range where the arm looks good — and it sits
outside the product's operating point by an explicit decision, D4 having moved the reporting
metric onto short lists (`compare.py`: *"the interesting budgets are tens, not thousands"*). Real,
and off-target.

The results log reads this as *"only the monotone structure is specific to `chromatin_od`, the
ranker that can see chromatin density."* The prior audit's F1.5 already noted the competing
length-sensitivity explanation. Adding to it: `rho(Δn_detections, Δrecall@250)` on the 25
informative cells is **−0.866**, against `rho(area_ratio, Δrecall@250)` of +0.559 — on the
population the crossover table actually uses, the volume association is materially stronger than
the prior audit's all-35 figure of −0.773, and stronger than the dose association it is being
credited to.

---

## Part 2 — confirms the 2026-09-08 audit

Re-derived and agreed, not restated in detail: **F1.3** (the `recall_at_budget` rationale is false
on `chromatin_od`; my informative-cell figure of −0.866 is stronger than its −0.773), **F1.4** (the
effect is not purely a volume artefact), **F1.5**, **F1.6** (both duplicate pairs confirmed
byte-identical), **F1.7**, **F1.8**, **F1.9**, **F1.10**. Its Part-0 reproduction list also holds
against my own re-derivation.

## Part 3 — what stands, and needs no clustering

The cost side is not a statistical claim and does not depend on any of the above:

| quantity, 25 informative cells at z=1.0 | largest_cc vs base51 |
|---|---|
| `n_detections` | **longer in 25 of 25**, mean **+3,851** |
| `read_95` | **deeper (worse) in 20 of 25**, mean **+524** |
| `coverage_frac` | never worse (mean +0.03) |
| sign of Δ recall@250, `chromatin_od` | 17 neg / 6 pos / 2 zero |
| sign of Δ recall@250, `tm_score` | 14 neg / 7 pos / 4 zero |

Tightening has never been observed to *help* in this run on either ranker; every pooled and
per-axis point estimate at every budget from 25 to 500 is negative or zero; and it unambiguously
costs more candidates to read.

## Part 4 — declared, not chased

The `BORDER_REPLICATE` padding asymmetry (`PAD` = 25 px at 51 px against 9–12 px at 19–25 px) is
disclosed in the docstring, the pre-registration and the prior audit, and its magnitude has never
been quantified by anyone. Doing so needs per-candidate positions the CSVs do not store. Recorded
as an open, unquantified limitation; not worth re-running the sweep for.

**One latent code bug, provably inert here.** `f1_seed_sweep.draw_seed_with_retry` retries with
`working = working.drop(working.index[idx])` — a drop by *label*, not by position. The pool comes
from `ss.agreement_pool` → `ss.border_filter`, which carry the annotation frame's index through;
if that index were ever non-unique, one retry would drop every row sharing the label and the pool
would shrink faster than `retries` records. It cannot bite in this run — `n_retries` is 0 in all
35 cells, so the loop never iterates — and F2 does not inherit it (it draws via
`tm_variant_sweep.draw_seeds`). Fix is `working.iloc[[i for i in range(len(working)) if i != idx]]`
or a positional mask, and it is worth doing before the function is reused.

---

## Verdict

**Does bbox tightening beat the fixed 51 px template? No — and the run cannot show the reverse
either.**

* At the ROI level, on either ranker, the recall difference is **not distinguishable from zero**:
  `tm_score` −0.0153 [−0.045, +0.018]; `chromatin_od` −0.0129 [−0.028, +0.002].
* **At every budget the product actually reads — 25 to 500 — the point estimate on both axes is
  negative or zero.** Past K ≈ 1000 on `tm_score` it turns positive (+0.019 at K=2000, winning in
  16 of 25 cells), which is the longer list reaching deeper, ROI-clustered non-significant
  (p = 0.22), and outside the operating point D4 fixed. Direction is otherwise consistently
  against tightening; magnitude is not resolvable with 7 ROIs.
* The **workload cost is measured and unambiguous**: +3,851 candidates in 25 of 25 cells, `read_95`
  deeper in 20 of 25.
* The **"~36 px threshold" is a hypothesis, not a measurement**: located by eye at the midpoint of
  the area-ratio scale, flat across 33–45 px, positioned largely by a duplicated draw, absent on
  the shipped ranker, and resting on a dose axis whose low end comes from seeds production would
  refuse.

So the operational call — keep the fixed 51 px template — is right, but it is a **decision under
uncertainty backed by a measured cost**, not a measurement that tightening is worse.

### Scope of that verdict — it is narrower than "bbox tightening"

Three fences, all of which a reader could easily jump.

1. **`tighten_bbox` does two things and F1 varies one of them.** `experiment.run_find_and_suppress`
   couples (a) narrowing the seed pool through `seed_selection.foreground_filter` inside
   `pick_seed` with (b) resizing the template through `tightened_base_size`
   (`midog_utils/experiment.py:142-154`). F1 holds (a) fixed — `largest_cc_box` is its `check_fn`
   in **both** arms — and varies only (b). **The verdict above is about (b) alone.** No controlled
   run in this repo has isolated (a); `find_and_suppress_midog_raw_seed.ipynb` runs
   `tighten_bbox=False` throughout and names the coupling as the reason its own
   True-vs-False contrast would be confounded.
2. **F1 never tested *fitting*, only *sizing*.** The template stays click-centred at every size
   (prior audit F1.9). "Does moving the box onto the nucleus help?" is untested by anything.
3. **`51 px` is not the current harness default.** `experiment.py:97` is
   `tighten_bbox: bool = True`. The recent TM work (`tm_threshold_axis_sweep_v2`, `_largest_cc`,
   F1, F5) bypasses that entry point and configures `fs.FSConfig(base_size=tm.BASE_SIZE)` directly.
   So "keep 51 px" is a **change** to the harness default, and setting `tighten_bbox=False` to get
   it would also switch off the seed-pool filter — the half F1 did not measure.

### And the same segmentation has live evidence *for* it, in a different role

`DECISIONS.md` D5: `mask_od_mean` — the mean OD **inside the largest Otsu component**, i.e. a
tightening-derived region — beats the fixed 51 px window statistic `od51` on 2-class AUC (6 wins,
1 tie), on the look-alike contrast (5/7) and on `read_95` (6/7), computed **gate-free** with a
0.0 % failure rate on all seven ROIs (`tp_fp_feature_extract.shape_features`). D5 records that no
chromatin-family axis has had a seed sweep and that any re-measurement owes the 14-ROI
`images/extra_valid` draw.

So the Otsu component is not a dead end. It is unhelpful as a **template size**, and prima facie
useful as a **ranking region**. Reading F1 as "stop segmenting the nucleus" would be reading it
several fences past where it stops.

**`f2_base_size_dose.py`** (pre-registered `b571811`, currently untracked) is the design that
settles it, and it independently fixes most of the above: dose **assigned** rather than observed,
seven levels **within** each cell so the contrast is within-ROI, 14 ROIs per D5, seeds drawn
**without** replacement, and a `matched_n` threshold rule alongside `matched_z` that closes the
volume confound of F1.3/N6. Its primary is *"a mean per dose with a cluster interval and nothing
else"* (commit `b571811`), which is N1 already answered. The one thing this audit would add: if
`tighten_box_otsu` is the rule that ships, F2's observational `lcc` arm should record how often
the shipped rule's **refusal** would have changed which seed entered the pool (N5) — on F1's
7-ROI draw that was 5 of 35, and it selected the two extreme low doses.
