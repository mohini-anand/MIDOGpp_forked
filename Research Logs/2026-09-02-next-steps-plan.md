# Next steps: what to do Monday, and why

**Status: plan only. No experiment was run to write this.** Every number below is either
re-derived by me from the CSVs already in `results/` (marked *[verified]*), taken from a prior
log and *not* re-checked (marked *[unverified]*), or a design estimate (marked *[estimate]*).

**The metric this plan adopts was validated before it was used.** `recall_at_budget` /
`tp_at_budget` had never been audited — no prior log reads them and neither audit covered them — so
adopting them as the standard (M6) would otherwise have been one more unchecked adoption from a
script with 15 known defects. Four checks, all three `tm_*` CSVs *[verified]*:
`recall_at_budget == tp_at_budget / n_gt_mitotic` on every row (0 violations); monotone
non-decreasing in `budget` within every cell (0); and — the one that actually discriminates, since
`read_50` and the budget curve are independently computed views of the same ranking — **recall
crosses 0.5 exactly between the grid points bracketing `read_50` in 4,276 of 4,276 cells (0
violations)**. The two views agree.

The one real caveat: `budget_delivered < budget` where the pool is shorter than the budget, so
recall there is "recall at end of list", not "recall at N". This affects **0 rows in
`tm_ccoeff_headtohead.csv`** (every table in sections 1 and 2), **0 decision-grade rows in
`tm_variant_stage_b.csv`** (the 5 affected rows are all 405.tiff), and 69 decision-grade rows in
`tm_variant_sweep.csv`, from which this plan quotes no table.

This plan was written against a briefing produced by the agent that did the 2026-09-01/02
template-matching work, which documents its own framing bias in
`Research Logs/2026-09-01-tm-variant-sweep.md` §9. I treated every claim in that briefing as an
unverified hypothesis and re-derived what I relied on. Two independent audits ran alongside this
plan — one re-deriving the sweep's numbers from the CSVs without reading the log, one auditing the
four never-logged CSVs and their scripts. Their findings are attributed where used.

**Eight claims do not survive in the form they were stated** (U1–U8 in section 1), including the
session's headline recommendation, the file that framed the session, and the file that supplied its
load-bearing prior. **One came back stronger than claimed** (U0) — the direction of the corrections
is not uniformly deflationary, and where the evidence favours the previous session's conclusion I
say so.

**The correction runs in both directions.** The briefing's own over-correction (§5e, "the
click-seeded pipeline beats the click-free blob detector") is also unsupported at this sample size —
see U2. Neither "blobs win" nor "TM wins" is a conclusion these 10 ROIs can carry.

The single most decision-relevant finding is not in any existing log:

> **In the configuration the previous session recommends shipping — TM peaks as the candidate
> generator, chromatin density as the ranker — the pathologist's click buys nothing measurable.**
> A synthetic disc template with no click produces the same median recall at every reading budget on
> 10 decision-grade ROIs (3–6 wins of 10, median difference +0.002 to +0.009 recall, sign-test p
> never below 0.125); so does a click on a pathologist-*rejected* look-alike. **The click-free disc
> is no worse than the real click at any budget under any of three seed aggregations**, and on the
> worst-of-5-seed gate — the one this project's own product notes say to use — the direction is
> consistently *against* the click at every budget from 100 up. Individual p-values there reach
> 0.016–0.039 but **do not survive correction for the 24 comparisons run**, so read this as *no
> click advantage established*, not as a demonstrated click penalty. The null is what matters and
> the null is solid. *[verified, section 2]*

That is not "the product premise is dead". It is "the product premise is not supported by *this
representation*", and the repo already contains direct evidence that the same click is highly
informative in a different one. Section 4 is about which representation to try next; section 3
puts the experiment that would settle it first.

---

## 0. Terms

Same as `2026-09-01-tm-variant-sweep.md` §0, with two additions this plan leans on.

**`recall_at_budget`** — the fraction of an ROI's mitotic figures found in the top *N* candidates.
Already computed and sitting unused in every `tm_*` CSV at N = 25, 50, 100, 250, 500, 1000, 2000,
5000. It is a better product metric than `read_50` for two reasons: it is defined for every cell
(no `NaN`-when-unreached bookkeeping), and *N* is the quantity the product actually budgets. I use
it throughout.

**Decision-grade** — `n_mitotic >= 15`. With the four ROIs downloaded on 2026-09-02 this is now
**10 ROIs, not 6**: 094, 201, 202, 245, 246, 300, 301, 402, 459, 548 — and **three are human**
(094 breast/Hamamatsu S360, 402 neuroendocrine/Hamamatsu XR, 548 melanoma/Hamamatsu XR), plus one
new canine (459 soft-tissue sarcoma/3D Histech). *[verified]* The tm-variant log's repeated caveat
that "no human ROI is decision grade" and "roughly three effective domains arbitrate" is **stale**:
it describes the 13-ROI state, and `tm_ccoeff_headtohead.csv` and `tm_variant_stage_b.csv` both
already contain all 17 ROIs. The log was not updated after the re-run. Nothing in it is wrong
because of this; it is simply more conservative than the data now warrants.

> **Correction (2026-09-04).** Five more ROIs were downloaded on 2026-09-04 — 013 (human breast /
> Hamamatsu XR), 233 (canine lung / 3D Histech), 403 (human neuroendocrine / Hamamatsu XR), 460
> (canine soft tissue sarcoma / 3D Histech) and 529 (human melanoma / Hamamatsu XR) — bringing
> `images/` to **23 ROIs**. **Decision-grade is now 15, not 10**: the ten named above plus those
> five. Six are human (013, 094, 402, 403, 529, 548) and all 7 domains are represented.
>
> A second bar governed that draw, which this section does not define. **Valid** = `n_mitotic >= 15`
> **and** a seed pool — `agreement_pool` → `border_filter(36, roi_shape)`, the same pool
> `tm_variant_sweep.draw_seeds` uses — of **>= 5**. On the 23 downloaded ROIs that gives **14, not
> 15**: 202 is the sole ROI clearing the mitosis bar but not the seed bar (16 mitotic figures, of
> which only 4 are unanimous and survive the border filter). So `tm_variant_report.py`'s
> `decision_grade & (n_seeds >= 5)` is load-bearing rather than belt-and-braces — it drops exactly
> one ROI. **Every domain now has exactly 2 valid ROIs**, which is what the 2026-09-04 draw was for.
>
> Everything computed below is unaffected: those numbers come from `tm_ccoeff_headtohead.csv` and
> `tm_variant_stage_b.csv` as they stood on 2026-09-02, and have not been recomputed. Where this
> document says "10 decision-grade ROIs" of a *result*, that remains the correct description of
> what was measured.

---

## 1. What is established, what is suggestive, what is unsupported

Re-derived from the CSVs, at ROI level (10 decision-grade ROIs), because 5 seeds within an ROI are
not independent replicates.

### Established

**E1. Ranking a dense candidate set by absolute chromatin density is the single largest effect in
the project, and it is click-free.** *[verified]* Every arm that ranks by chromatin lands in a tight
band, and every arm that ranks by a TM score is worse or much worse. Pooled median
`recall_at_budget` over the 10 decision-grade ROIs, single-template run
(`tm_ccoeff_headtohead.csv`):

| budget | `tm_ccoeff_od` | `disc_ccoeff_od` | `lookalike_ccoeff_od` | `tm_ccoeff_normed_od` | `tm_ccoeff` | `disc_ccoeff` | `blob_native` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 25 | 0.153 | 0.132 | 0.144 | 0.153 | 0.161 | 0.027 | 0.096 |
| 100 | 0.439 | 0.437 | 0.432 | 0.446 | 0.432 | 0.087 | 0.323 |
| 250 | 0.657 | 0.644 | 0.662 | 0.666 | 0.613 | 0.273 | 0.570 |
| 500 | 0.795 | 0.782 | 0.790 | 0.790 | 0.727 | 0.415 | 0.777 |
| 1000 | 0.864 | 0.872 | 0.872 | 0.864 | 0.798 | 0.534 | 0.858 |
| 5000 | 0.961 | 0.960 | 0.968 | 0.946 | 0.945 | 0.853 | 0.943 |

The first four columns are indistinguishable. They differ in whether there was a click, whether the
click was on a mitotic figure, and which `cv2.TM_*` constant generated the peaks.

**And this survives rotation augmentation, which closes the one loophole that could have rescued the
click.** The tm log §5e states that under augmentation "ranking augmented `TM_CCOEFF` by its own
match score is as good as chromatin ranking (125 vs 126)" — which would matter enormously, because
score ranking *is* strongly click-dependent (E4) while chromatin ranking is not. That statement rests
on median `read_50` over 6 ROIs. On 10 decision-grade ROIs and `recall_at_budget` it does not hold
*[verified]*: chromatin ranking is better on 5–6 of 10 at the median, and on the worst-of-5-seed gate
it is better on **8 of 9 non-tied ROIs at budgets 500 and 1000, p = 0.039**. So there is no
configuration measured here in which a click-dependent ranker is the right choice.

**E2. Goal (2) is not met, and is not close.** *[verified]* At 5,000 candidates the best
configuration has **96.5%** of the mitotic figures. The remaining ~3.5% costs the other ~10,000
candidates in the pool: augmented `tm_ccoeff` median `read_100` runs 2,602–14,888 against pools of
13,756–17,249, i.e. **19%–86% of the entire candidate list**, and many seeds never reach 100% at
all. "Reach full recall in as few candidates as possible" currently reads as "reach full recall by
reading most of the list."

I lead with recall@5000 rather than `read_100` deliberately: `read_100` is hostage to individual
unreachable annotations (300.tiff ann 14512 sits 3 px from the image edge and pins that ROI's
ceiling regardless of method — *[unverified, log §5f]*), whereas the recall@5000 statement does not
depend on any single object.

**E3. Both candidate generators tile most of the ROI, and there is now a usable rule for when a
ceiling means anything.** *[verified]* `coverage_frac` at full pool: TM arms 0.882–0.953, disc arms
0.826–0.935, `blob_native` 0.643–0.992. An independent audit of `results/blob_coverage.csv` — the
only one of the four unlogged CSVs that comes out sound — gives the criterion this project has been
missing: **a full-list ceiling of 1.0 carries essentially no information once `cov^n_gt ≳ 0.05`**,
i.e. once a random point pattern of the same density would have hit every annotation anyway. Pool
size is *not* a usable proxy (34,910 candidates gives coverage 0.992 on 245 while 40,829 gives 0.951
on 201); the test must be on measured coverage.

Applied to the 10 decision-grade ROIs *[verified]*, this kills exactly the ceilings you would expect
and one you would not:

| ROI | `n_gt` | TM coverage | `cov^n_gt` | TM ceiling informative? | blob `cov^n_gt` | blob informative? |
|---|---:|---:|---:|:--|---:|:--|
| 201 | 17 | 0.893 | **0.146** | **no** | 0.0089 | yes |
| 202 | 15 | 0.909 | **0.241** | **no** | 0.0144 | yes |
| 245 | 89 | 0.953 | 0.013 | yes | **0.507** | **no** |
| 246 | 115 | 0.948 | 0.0023 | yes | 4.5e-5 | yes |
| 094 / 300 / 301 / 402 / 459 / 548 | 81–238 | 0.882–0.950 | ≤3.7e-5 | yes | ≤7e-7 | yes |

Adopt this as the standing rule and report `cov^n_gt` beside every ceiling.

**U4 (consequence, and it contradicts a claim on record).** The tm log's P4 survives on exactly one
ROI — "only **201.tiff** survives the coverage objection: +0.12 ceiling for +0.07 coverage" (§5).
201.tiff is the ROI where TM's ceiling is **least** informative of the ten: `n_gt` = 17 and
`cov^n_gt` = 0.146, so the coverage null cannot be rejected there at all. *[verified]* **The single
ROI carrying the P4 claim is the single ROI where the claim cannot be made.** P4 should be recorded
as unresolved, not as "confirmed on 201.tiff". (Symmetrically, `blob_native`'s ceiling of 1.0 on
245.tiff is also uninformative, `cov^n_gt` = 0.507 — the tm log half-noticed this when it found blob
dropping from 5/5 to 0/5 there at matched pool.)

**E4. The click is worth a great deal when the TM score is the ranking key.** *[verified]*
`tm_ccoeff` vs `disc_ccoeff`, chromatin *not* used: 10/10 ROIs at budgets 25–500,
sign-test p = 0.0020, median recall difference **+0.40 at budget 250**. This is real and large. It
is also irrelevant to the recommended pipeline, which does not rank by the TM score — see U1.

### Suggestive but not established

**S1. TM peaks are a better candidate generator than blob centroids at small budgets.**
*[verified]* Augmented `tm_ccoeff_od` vs `blob_native`, ROI-level paired, 10 ROIs:

| budget | TM better | worse | tied | sign p | median Δrecall |
|---:|---:|---:|---:|---:|---:|
| 25 | 7 | 1 | 2 | 0.070 | +0.061 |
| 100 | 8 | 2 | 0 | 0.109 | +0.084 |
| 250 | 8 | 2 | 0 | 0.109 | +0.082 |
| 500 | 7 | 2 | 1 | 0.180 | +0.043 |
| 1000 | 6 | 3 | 1 | 0.508 | +0.028 |
| **2000** | **5** | **5** | 0 | 1.000 | +0.015 |
| 5000 | 6 | 3 | 1 | 0.508 | +0.025 |

Never significant at any budget. The advantage is real in direction and confined to the head; by
budget 2,000 it is a dead heat and `blob_native`'s pooled median (0.930) is *above* TM's (0.906).
**Click-seeded wins the head, seedless wins the middle tail** — which, together with the
disjoint-failure pattern in §5c of the tm log, is the actual case for a union generator.

**S2. TM's advantage may be a coverage confound.** *[verified, inconclusive]* Across the 10 ROIs,
Spearman(Δcoverage, Δrecall@250) between TM and blob = **0.382, p = 0.276**; Pearson 0.464,
p = 0.176. TM's biggest wins (094 +0.30 recall, 201 +0.35) are also its biggest coverage
advantages (+0.24, +0.14). Ten points cannot resolve this. It is the most important unresolved
threat to S1 and it needs an equal-coverage comparison, not an equal-*n* one.

**S3. The blob detector and TM fail on disjoint ROIs, and the under-segmentation half of that story
now has direct support.** *[log §5c unverified; the watershed half verified by independent audit]*
The mechanism claimed is `label(connectivity=2)` fusing touching nuclei. An audit of
`results/watershed_sweep.csv` confirms the fix is **real detection, not carpeting**: on 246.tiff the
four annotations the connected-component arm missed and the watershed arm found match at
**2.4–6.5 px** at ranks 89–624, and four of the connected-component arm's carpeting-signature matches
(26.7–28.2 px, two at ranks beyond 17,000) become 3.9–10.9 px at ranks 49–3,081 under splitting. No
annotations are lost. Splitting genuinely recovers fused nuclei.

But the same audit says the sweep **cannot choose the constant**: no `MIN_PEAK_DIST` value
dominates, every value still regresses depth-to-100% on 301 and 245 versus not splitting at all
(best case 1.21x and 1.17x), and 246's depth-to-100% swings 44,663 → 3,082 between `md`=5 and `md`=7
— i.e. it is set by the rank of a single annotation. Also, splitting *is* the knob for the ceiling,
not `md`: every value from 5 to 15 fixes 246's ceiling to 1.0.

So the union recommendation is better motivated than the tm log establishes, and still imports the
blob detector's four untuned constants plus a fifth (`MIN_PEAK_DIST`) that this sweep is not
entitled to set. Gate it behind P2.

### Unsupported as stated

(U4 is stated above, next to the coverage criterion it follows from.)

**U0 (for balance — the claim that came back *stronger* than stated).** An independent re-derivation
confirms claim 1 and finds it understated: on `read_50`, `TM_CCOEFF` beats `TM_CCOEFF_NORMED` on
**49/49 decision-grade cells and 10/10 ROIs, p = 0.00195, median ratio 4.32×** (148 vs 710), rising
to 6.23× over all 17 ROIs. `read_50` has **zero NaN in either arm**, so NaN handling cannot flip it;
and where NaN does appear (`read_100`) the normed arm has *more* of it, so dropping NaN pairs removes
`TM_CCOEFF` wins — the bias runs against the claim, not for it. The coverage objection also fails
here in the claim's favour: `TM_CCOEFF_NORMED` has *higher* coverage (0.925 vs 0.885) and still reads
4.3× deeper, which is the opposite of what a tiling artefact predicts. **`read_50` is solid. The same
comparison's `full_list_recall` is not** — median ratio exactly 1.000, 24 of 49 cells tied, both arms
at ceiling; that half is geometry-dominated and should be dropped. U1 is about *relevance*, not
validity.

**U1. "Switch `TM_CCOEFF_NORMED` to `TM_CCOEFF`; there is no argument for keeping the
normalisation" (§5d item 1) is contradicted by the same log's §5b and by the CSVs.** *[verified]*
Under chromatin ranking — the ranker the log itself recommends — the two methods are a dead heat:
`tm_ccoeff_od` 0.153 / 0.439 / 0.657 / 0.795 against `tm_ccoeff_normed_od` 0.153 / 0.446 / 0.666 /
0.790 at budgets 25/100/250/500. No consistent direction at any budget. The log's §5b already says
this ("median delta +1 candidate, p = 0.20") and §6 says it again; §5d item 1 then states the
opposite as a recommendation. **The 8.6x is a fact about the TM score as a ranking key, not about
the method.** The change is still worth making — it is one constant and it is strictly better in
the score-ranked configuration you might fall back to — but it must not be described as an 8.6x
product win.

**U2. "Augmented TM_CCOEFF beats a click-free blob detector on depth in 8/10 ROIs" understates how
weak that is, and the log's framing of it as the correction to its own bias overshoots.**
*[verified]* 8/10 gives p = 0.109. With 10 ROIs a sign test cannot go below p = 0.109 without
reaching 9/10. The honest statement is that TM and blob are **not distinguishable** at this sample
size at any budget, in either direction — which is neither the original "blobs win" nor the §5e
correction "the click-seeded pipeline beats the click-free blob detector."

Three further corrections to that claim, from independent re-derivation *[verified by audit]*:
* **The "8/10 on depth" belongs to the chromatin-reranked arm, not to template matching.**
  Score-ranked augmented `tm_ccoeff` vs `blob_native` on `read_50` is **6 better / 3 worse / 1 tied,
  p = 0.508**. Only `tm_ccoeff_od` gives 8/2.
* **The ceiling half of the claim is vacuous.** `full_list_recall`, `coverage_frac` and
  `n_detections` are **bit-identical** between `tm_ccoeff` and `tm_ccoeff_od` — same pool, only the
  sort key differs. So "5 better / 2 worse / 3 tied on ceiling" is *the same number for both arms*
  and carries no information about which one was used. It is claim 6 (chromatin ranking) restated,
  not evidence about template matching.
* **The briefing's "8/10 ROIs" is an extension the log does not make.** The log says "5 of 6
  decision-grade ROIs" and its §5e table is explicitly labelled *chromatin-ranked*; that table
  reproduces exactly. On its own 6 ROIs the result is 5/1, p = 0.219 — also not significant.

**U3. The recommended configuration has never been compared against a click-free version of
itself.** *[verified]* `tm_variant_stage_b.csv` contains exactly nine arms — `blob_native`,
`blob_native@matched`, `blob_od`, `tm_ccoeff`, `tm_ccoeff@matched`, `tm_ccoeff_od`, `tm_ccorr`,
`tm_ccorr@matched`, `tm_ccorr_od`. **There is no `disc_*` and no `lookalike_*` arm in the augmented
experiment.** Every statement about what the click is worth — the 6.9x, the look-alike null — comes
from the *single-template* run. The configuration actually recommended for shipping has no
click-value measurement at all. This is the highest-value gap in the project and P1 below closes it.

**U5. "`TM_SQDIFF`'s click performs worse than a click-free disc" reports the depth half and omits
the ceiling half, where it reverses.** *[verified by audit]* On `read_50` it holds at **seed 0 only**
(2/10, p = 0.039); over all 5 seeds it is 4/9, **p = 0.267 — not significant**. And on
`full_list_recall` over 5-seed medians the **click wins 9/0/4, p = 0.0039** (0.835 vs 0.767). The
mechanism is the coverage confound again: `disc_sqdiff` covers 0.508 of the ROI against the click's
0.818, so the disc pool reads faster because it is smaller and finds less. The correct statement is
that under `TM_SQDIFF` the click trades head depth for tail coverage — not that it is worthless.
(This does not rehabilitate `TM_SQDIFF`, which is last on every other measure.)

**U6. The `TM_SQDIFF_NORMED` withdrawal is right, but two of its three supporting numbers cannot be
checked from the CSVs and one is confusable with an unrelated figure.** *[verified by audit]* The
tie-block finding is **supported and method-specific** — `tie_frac` p90 0.148 / max 0.495 for
`sqdiff_normed` against max 0.0005 for all five other methods, and `tail_arbitrary=True` occurs on
that method and no other (129 rows). The `n_peaks ≥ 250000` cap on 228 rows across 10 ROIs, and
`sqdiff_normed` setting the minimum `n_detections` in 20 of 58 cells, both re-derive exactly. But
"83% of the map saturates" and "MAD collapses to exactly 0" come from a pilot probe on the raw
OpenCV response map and are **not derivable from any CSV column** — and there is a coincidental
**82.76%** in the CSVs (the share of `disc_sqdiff_normed` cells with `tail_arbitrary=True`, 48/58)
that is easy to mistake for it. Keep the withdrawal; re-source the 83% or drop it.

### The two unaudited files that framed this work — now audited

Both were flagged in the briefing as un-logged and never reviewed. I had an independent audit run on
them and the scripts that produced them. *[verified by that audit; I have not personally re-run its
arithmetic, and say so.]*

**U7. `tm_vs_blob_depth.csv` — the "4.3–216x" bar that framed the whole session. Direction survives;
three of its four quoted properties do not.**
* The ratios reproduce exactly (4.26–216.07 at read-50, 1.90–12.94 at 95%, 1.01–5.68 at 100%), and
  the fairness objections I expected do **not** apply: both arms score the *identical* cached pool,
  so `n_pool`, `ceiling` and `n_found` are bit-identical across arms on every ROI. This is not a
  short-list-vs-long-list artefact. No ground truth enters either ranking function.
* **The 100% column must be withdrawn.** `tm_vs_blob_depth.py:126-127` assigns a sentinel score of
  −2.0 to any candidate whose 51 px patch crosses the ROI edge, while `blob_score` scores those same
  candidates normally. On 246, 300 and 301 the reported depth-to-100% lands *inside* that sentinel
  block — so those cells measure the sentinel, not template matching. That includes the headline
  "1.01x".
* **"At every quantile" over-counts.** On 201, 202 and 405, `n_found` is small enough that the 95%
  and 100% depths are literally the same measurement.
* **"7/7 ROIs" is 4 domains**, two of which sit at or below the decision-grade bar (405 has 13
  mitotic figures and fails it). ~21 apparent confirmations are ~8 independent ones.
* **TM was not run as it ships**, and the gap is large: 24 templates against the shipped 72, no
  scale bank, no robust-z fusion, and a raw 51 px crop rotated with `BORDER_REPLICATE` (rotated
  corners are smeared border, not tissue) instead of the padded-then-cropped patch
  `template_match.py` exists to provide.
* **The deepest problem is one nobody had named: TM was scored *at blob centroids*.** Holding
  candidate *generation* constant does not hold candidate *representation* constant. A region mean
  is translation-insensitive; a 51 px correlation is not, and fused connected components have
  displaced centroids — the exact defect `miss_attribution.py` documents. The comparison was
  structurally tilted toward the blob detector by its evaluation geometry, independently of the
  framing bias §9 confesses to.

  **Net:** a 216x gap is not closed by these defects, so the read-50 and 95% *direction* stands. But
  the file was never entitled to frame an experiment, and its 100% column should be struck from the
  record.

**U8. `click_darkness.csv` is broken as evidence, and this reopens a question everyone treated as
closed.** This file is the "load-bearing prior" for prediction P1 in the tm log §3 and the source of
the claim that click-conditioned darkness loses to plain darkness by 30–345x.
* The quoted range does not reproduce under any aggregation tried; the script's own worst-of-5 gives
  **28.0–366.6**, and it collapses to **1.81–15.14 at 95%** and 1.00–15.14 at 100%. The headline was
  a head-only effect reported as a general one.
* **The steel-man arm is a mathematical no-op.** `od51_onesided` — documented as "keeps monotonicity
  while still using the click to set the scale" — is provably identical to the click-free `od51_abs`
  in **210 of 210 rows, maximum difference 0**. Its formula is strictly increasing in candidate
  darkness on both branches, so it is `od51_abs` re-sorted.
* Therefore the only two arms that actually differ from the incumbent (`od51_sim`, `feat_sim`) are
  **the same defect class**: two-sided distance to a single draw. The experiment contains **zero**
  arms testing a monotonicity-preserving use of the click. "The click adds nothing beyond a dark
  nucleus-shaped object" rests on one defect tested twice, and **is not established**.

  This matters for what to do next: it means the obvious cheap use of the click as a *ranker* — one
  that cannot rank a candidate below the clicked cell for being darker — has never been tried. P1b
  below tries it.

---

## 2. The finding that should change Monday

Under chromatin ranking, on 10 decision-grade ROIs including 3 human ones, a click on a real
mitotic figure is indistinguishable from **no click at all**.

**What is and is not in scope here.** The tm log §4 fixed the template at 51 px with no Otsu
tightening, deliberately, so that only the method varied. Both arms below therefore use a 51 px
template, and **what is measured is the click's *texture*, with the click's *size* affordance held
out of the comparison by design.** Template size is one of the two things a click supplies that a
seedless method structurally cannot get (§9 item 5 of the tm log makes exactly this point), and the
shipped pipeline does tighten. So this section does not close the question "is a click worth
anything"; it closes "is the clicked cell's appearance worth anything, at fixed template size". P1b
is where the size affordance gets its test.

**`tm_ccoeff_od` (real mitotic click) vs `disc_ccoeff_od` (synthetic disc, no click)**, ROI-level
paired, single-template: *[verified]*

| budget | click better | worse | tied | sign p | median Δrecall |
|---:|---:|---:|---:|---:|---:|
| 25 | 5 | 3 | 2 | 0.727 | +0.0021 |
| 100 | 3 | 5 | 2 | 0.727 | −0.0028 |
| 250 | 6 | 1 | 3 | 0.125 | +0.0087 |
| 500 | 6 | 2 | 2 | 0.289 | +0.0069 |
| 1000 | 2 | 4 | 4 | 0.688 | +0.0000 |
| 5000 | 3 | 3 | 4 | 1.000 | +0.0000 |

**`tm_ccoeff_od` vs `lookalike_ccoeff_od`** (click on a pathologist-rejected non-mitotic object):
4/3, 6/4, 1/6, 3/3, 4/2, 2/4 across the same budgets — no direction, no significance, and at budget
250 the *look-alike* click is ahead on 6 ROIs of 10. Ceiling: click vs disc is 3 better, 1 worse, 6
tied, median difference exactly 0.

An independent re-derivation on `read_50` pushes this further than "no difference": for
**`TM_CCOEFF_NORMED` the look-alike click is significantly *better* than the real one** — 2/11 over
13 ROIs (p = 0.0225) and **0/6 over the decision-grade ROIs (p = 0.031)**, winning on the two largest
(300: 1762 vs 3496; 301: 1690 vs 4110). *[verified by audit]* The honest reading is that this is 1 of
6 uncorrected tests (Bonferroni-adjusted 0.135) at one seed per ROI, so it is not a finding on its
own; but it is one more control pointing the same way, and none point the other way.

**Robustness of the section-2 result to the aggregation choice** *[verified]* — because "median over
seeds" was a free parameter I could have chosen after seeing the answer. Click vs disc, chromatin-
ranked, per-ROI aggregation swapped:

| aggregation | budget 250 | budget 1000 | budget 2000 |
|---|---|---|---|
| median | 6/1, p = 0.125 | 2/4, p = 0.688 | 1/6, p = 0.125 |
| mean | 7/3, p = 0.344 | 2/7, p = 0.180 | 1/8, **p = 0.039** |
| min (worst seed) | 3/5, p = 0.727 | 0/7, **p = 0.016** | 1/8, **p = 0.039** |

No aggregation produces a significant result in the click's favour at any budget. Two of three
produce significant results *against* it at budget 2000.

### And on the product's own gate, the click is worse than no click

The project memory records the correct product gate as **worst-of-5 seeds, not the median** — "a
ranker that wins on median and fails on one click in five is exactly the product defect already
identified". On that gate the comparison reverses and becomes significant: *[verified]*

**`tm_ccoeff_od` vs `disc_ccoeff_od`, worst of 5 seeds, ROI-level paired, 10 ROIs:**

| budget | click better | worse | tied | sign p | median Δrecall |
|---:|---:|---:|---:|---:|---:|
| 25 | 4 | 4 | 2 | 1.000 | +0.000 |
| 100 | 2 | 6 | 2 | 0.289 | −0.032 |
| 250 | 3 | 5 | 2 | 0.727 | −0.008 |
| 500 | 2 | 6 | 2 | 0.289 | −0.015 |
| **1000** | **0** | **7** | 3 | **0.016** | −0.023 |
| **2000** | **1** | **8** | 1 | **0.039** | −0.017 |
| 5000 | 1 | 7 | 2 | 0.070 | −0.014 |

Negative at every budget from 100 up. **But apply the same multiplicity standard I applied to the
look-alike result above, or this is a double standard running the other way.** These are 3 nominal
hits out of **24 tests** (8 budgets × 3 aggregations); ~1.2 false positives are expected at α = 0.05,
Bonferroni across 24 gives α = 0.002, and even correcting within the min-aggregation family alone
(8 budgets, α = 0.00625) neither 0.016 nor 0.039 clears. **No individual p-value here survives
correction.**

What does survive is the *direction*, which is negative at every budget from 100 up under all three
aggregations — consistent direction across correlated tests is real evidence even when the
individual p-values are not. So the defensible claim is **"the click-free arm is no worse, and
nothing suggests the click is better"**, not "the click is worse". The distinction matters because
the plan's recommendations rest on a null, and a null does not need the stronger claim.

The mechanism behind the direction is nonetheless visible in the data: `disc_ccoeff_od` produces
**one** candidate pool per ROI regardless of seed (its recall spread across seeds is 0.000–0.067,
and that residual is only the evaluation ground truth dropping a different annotation each time),
whereas `tm_ccoeff_od` produces five different pools with a recall spread up to **0.125** on
201.tiff and 402.tiff. *[verified]*

**State the structural caveat with it.** Worst-of-N necessarily penalises the arm that varies, and
the click-free arm is near-deterministic by construction — the tm log's §9 item 3 makes exactly this
objection in the opposite direction. So this is not evidence that a click *harms* retrieval. It is
evidence that, at this effect size, **a one-click product is worse than a no-click product for the
unlucky pathologist**, which is the risk the product must actually price. It is also the strongest
available argument for multi-click (section 4c): averaging over 3 clicks is the standard fix for
exactly this variance.

**What I deliberately did not compute:** the same worst-seed contrast against the look-alike arm.
The look-alike control has **one seed per ROI** *[verified]*, so its "worst of 5" is a single draw,
and comparing it to a genuine worst-of-5 would be biased against the mitotic click. The median-level
look-alike comparison above is the only valid form of that control at present, and section 6 lists
fixing the look-alike seed count as a cheap change.

Three things make this stronger than the equivalent statement in §6 of the tm log (which was
"disc ties or beats the real click on 3 of 5 ROIs" on worst-of-5 `read_50`):

1. **10 ROIs, not 5**, and three of them human — the domain caveat that qualified everything in
   that log does not apply here.
2. **The evaluation denominators are equal.** I checked: `n_gt_mitotic` is `total − 1` for *every*
   arm on every ROI, click-free arms included, so the seed annotation is dropped symmetrically and
   the comparison is not tilted either way. *[verified]*
3. **It holds at every budget**, not at one aggregate.
4. **The disc control is genuinely click-free.** I read `disc_template` in `tm_variant_sweep.py`:
   it is a fixed 51 px binary disc whose inside/outside values are *this ROI's* 50th and 99.5th OD
   percentiles, built once per ROI and explicitly independent of the seed; `match_pool` takes
   `seed_xy=None` for it. *[verified]* Note what that implies — the disc is **stain-adaptive without
   a click**. Per-slide stain calibration, one of the two things the click was supposed to buy that
   a seedless method structurally cannot have, is available from the image's own OD percentiles.
   That leaves per-case *morphology* adaptation as the only remaining click-specific claim, and the
   look-alike control is the direct test of it — which it fails.

**Say it as a pool statement, not as "the click is useless".** The click *does* change the candidate
set: TM pools are 14,234–17,503 against disc's 12,625–15,956, and TM coverage is +0.02 to +0.07
higher. The click enlarges and reshapes the candidate set by 10–20% — and none of that buys
measurable recall at any reading budget, because the ranker cannot express what the click knows.

**Why this is mechanistically unsurprising, and why it is not a death sentence.** Chromatin density
is a function of the candidate patch alone. It is invariant to everything the click carries except
darkness, so a chromatin ranker *cannot* express click information by construction. The result is
therefore **representation-bound, not interaction-bound**. The counterweight is already in the repo:
`2026-09-01-click-ranking-experiment.md` reports that in a mitosis-task-trained feature space one
click cut reading depth to 60–80 candidates on 246.tiff against the best seedless arm's 193, and
that a look-alike click cost 1,083–10,640 on the same ROI. *[unverified, and that log states plainly
that its encoder was trained on MIDOG++ including 246 and 301, so the numbers are an upper bound
reproducible by pure memorisation.]* But it is direct evidence that the click carries rankable
information that *some* representation can extract, and that neither chromatin nor any frozen
off-the-shelf encoder tested so far is that representation.

---

## 3. Experiments, prioritised

Each has: the decision it settles, cost, a pre-committed success/failure criterion, and what result
would falsify the current direction. Timings marked *[estimate]* unless taken from the tm log's own
measured runtimes.

### P0. Commit everything. 15 minutes. Not an experiment.

Nothing from 2026-09-01/02 is in git: `tm_variant_sweep.py`, `tm_variant_report.py`,
`new_domain_roi_check.py`, every `results/tm_*` CSV, and the research log are untracked, and five
`midog_utils/` modules are modified-uncommitted. No number in the tm log ties to a code version, and
the log itself flags this as defect M10 from the premise-test audit repeating. **Do this before
running anything else**, or P1's results will not be attributable either. Commit the working tree as
it stands first (so the existing CSVs tie to the code that made them), then start P1.

### P1. Add `disc_*` and `lookalike_*` arms to the augmented run. ~25 min. **Do this first.**

**Decision:** does the pathologist's click contribute anything to the configuration we would ship?

The augmented run (`--exp3`) has no click-free control (U3). The code for both controls exists and
runs in `--exp2`; this is an arm-list change plus a re-run over the 10 decision-grade ROIs. The tm
log measures `--exp3` at 20 min for 17 ROIs with 3 arms; adding 2 arms over 10 ROIs is the same
order. *[estimate]*

**Registered prediction, before the run:** rotation-and-flip augmentation takes an element-wise max
over 8 orientations of the same template, which moves the fused response toward rotational symmetry
— i.e. toward what the disc already is. So augmentation should **reduce** click specificity, not
increase it. I predict augmented disc matches augmented click at least as closely as the
single-template pair does.

**Pre-committed criterion — two-dimensional, deliberately.** Compare `tm_ccoeff_od` against
`disc_ccoeff_od` at budgets 100, 250, 500, ROI-level paired over the 10 decision-grade ROIs, on
**both** the median-seed and the worst-of-5-seed curve. "Wins" below means ≥9 of 10 ROIs at ≥2 of
the 3 budgets.

**Gating on the worst-seed curve alone would be rigged, and I nearly did it.** Section 2 establishes
that `disc_ccoeff_od` produces one pool per ROI while `tm_ccoeff_od` produces five — so worst-of-N
structurally penalises the arm under test, exactly the objection the tm log's §9 item 3 raises in the
opposite direction. A gate known to be biased against the hypothesis is not a test of it. So:

| median curve | worst curve | verdict | what to do |
|---|---|---|---|
| click wins | click does not lose | **click justified** | ship the click-seeded generator |
| click wins | click loses | **click has signal, unusable variance** | build the 3-click prototype (§4c) — this is a *positive* result for multi-click, not a negative for the click |
| click does not win | either | **click not justified in this representation** | the premise restates; go to P3/P4 |

Section 2 predicts row 3. Row 2 is the outcome a one-dimensional gate would have misreported as a
failure, and it points at a different Monday.

**If row 3 — the product restates.** Not "one click, then read a ranked list", but "a seedless ranked
list, and the click means confirm/reject". That is a different product with a different value
proposition, and it should be said out loud rather than discovered later. It also makes P3 (a real
detector) the main line rather than a side quest.

**Also run 5 look-alike seeds per ROI, not 1**, in the same pass. It is the same script change and it
converts the look-alike control from "can refute a large effect" into a measurement.

**What would falsify the current direction:** this experiment coming back "click not justified" is
itself the falsification of the click-as-search premise *in this representation*. It does not
falsify the premise generally — P4 is the test that does.

### P1b. A monotonicity-preserving click-conditioned ranker. ~2 h. Run alongside P1.

**Decision:** can the click improve the *ranker* at all, cheaply, without a learned representation?

U8 shows this has never actually been tested. The one arm designed to test it was algebraically
identical to the click-free incumbent in all 210 rows, and the two arms that did differ were both
two-sided distances to a single draw — which fail for a reason that is now well understood: darkness
is monotone with respect to the label, so a symmetric similarity penalises every candidate *darker*
than whichever mitosis was clicked.

The untested family is a ranker that uses the click for **scale and shape but never inverts the
monotone direction**. Concrete candidates, all cheap and all computable on the pools P1 already
builds:
* **One-sided in darkness, two-sided in everything else.** Rank by chromatin density as now, but
  break ties (or apply a bounded multiplicative adjustment) using similarity in click-derived
  features that are *not* monotone with the label — eccentricity, area, texture energy, the
  Otsu-tightened size the click gives for free.
* **A per-slide threshold set by the click.** The click supplies "at least this dark counts". Use it
  to set the operating point rather than the ordering — a use the disc control cannot replicate,
  because the disc's percentile-based scaling is a property of the ROI, not of a cell a pathologist
  endorsed.
* **Negative-click hard negatives**, once any candidate is rejected (section 4c).

**Add the size affordance as its own arm.** Section 2's null is measured at a fixed 51 px template,
so it says nothing about the Otsu-tightened size the click supplies for free — one of only two
things a seedless method structurally cannot obtain. The cheapest test is a `disc_ccoeff_od` variant
built at the *click's tightened size* rather than at 51 px: if that beats the fixed-size disc, the
click's value is size, not appearance, and the product claim changes accordingly.

**Pre-committed criterion:** the same two-dimensional gate as P1 (median curve and worst curve
reported separately, ≥9 of 10 ROIs at ≥2 of 3 budgets), so that "signal present but too variable"
is distinguishable from "no signal" here too.

**Why run it even though I expect it to fail:** it is two hours, it closes a question the record
currently answers with a broken experiment, and a negative here is what makes the "the click belongs
in verification, not search" conclusion safe to act on rather than merely likely.

### P2. Download held-out ROIs and re-run the decision table on them. ~2 h wall clock, mostly download.

**Decision:** does anything measured here generalise beyond the 18 training ROIs?

`datasets_xvalidation.csv` in this repo defines the official split: **111 test slides across all 7
tumour types** (canine mast cell 11, canine lung 10, canine lymphoma 12, canine STS 22, human breast
33, human melanoma 11, human neuroendocrine 12). **All 18 ROIs currently in `images/` are `train`.**
*[verified]* There is no held-out evaluation anywhere in this project and the split has been sitting
in the repo the whole time. This is the cheapest large improvement available.

> **Correction (2026-09-04).** `images/` now holds **23** ROIs, not 18 — and the claim survives the
> change: **all 23 are still `train`**. The five added on 2026-09-04 were drawn uniformly at random
> (`np.random.default_rng([20260904, i])`, `i` indexing the alphabetically sorted domains that
> lacked a second valid ROI) over that domain's candidates — the ROIs meeting the valid bar and not
> yet in `images/`, sorted by file name, taken at `rng.integers(len(candidates))` — with no
> restriction on split; the candidate pools
> are train-heavy, so all five landed in `train`. That satisfies rule (a) below but not P2 itself.
> Rule (b)'s target — "2 per tumour type, 14 ROIs" — is now numerically met by the *valid* set
> (2 per domain, 14 total), but by downloaded training ROIs, not the held-out draw this section
> asks for. **P2 remains open.** Held-out candidates were available at the time of the draw: 5, 5,
> 8, 4 and 1 `test`-split qualifying ROIs in canine lung, canine STS, human breast, human melanoma
> and human neuroendocrine respectively — note human neuroendocrine has only **one**, so a
> test-only draw is nearly forced there.

**Two selection rules to pre-commit, because the current set violates both.** The four ROIs added on
2026-09-02 were chosen as the *densest* qualifying image in each domain, which biases every depth
metric optimistically. For the held-out draw: **(a)** sample uniformly at random among slides
meeting the inclusion bar, with the RNG seed recorded in the commit message; **(b)** fix the target
count *before* looking at any result — I suggest 2 per tumour type, 14 ROIs, which roughly doubles
the decision-grade set. Also re-state the inclusion bar as a rule rather than an accident: it
currently excludes 405.tiff by a single annotation, which is arbitrary. Either lower it to
`n_mitotic >= 10` or keep 15 and say why, but decide before you draw.

**Pre-committed criterion:** report the P1 and S1 contrasts on the held-out ROIs alone, and on the
pooled set, separately. If the held-out direction differs in sign from the training direction on
either contrast, every conclusion from the 2026-09-01/02 sessions is provisional and must be
re-litigated.

### P3. Run the MIT-licensed MIDOG++-trained FCOS detector as the candidate generator. ~1 day.

**Decision:** is this whole line of work optimising the wrong stage?

`jonas-amme/FCOS_Inference_CLI` is MIT-licensed with public weights, trained on MIDOG++, and reports
F1 0.737–0.753 on the MIDOG 2022 test set; the same checkpoint is the MIDOG 2025 Track 1 reference.
*[unverified — from `2026-09-01-one-click-retrieval-literature.md` §A1, which cites the repos; I did
not download or run it.]* Measured detector timing on this laptop: `fcos_resnet50_fpn` at 1024 px
tiles ≈ 87 s for a full ROI *[unverified, same log §3]* — an ingest-time cost, not a click-time one.

The framing to interrogate: **every experiment in this project compares two ways of producing a
~15,000-candidate list.** A detector at F1 ≈ 0.75 produces a few hundred. That is two orders of
magnitude off the axis everything has been measured on, and it changes the click's job from
"re-rank 20,000 things" to "re-rank 300" — a different and much easier problem, where a weak
click signal might actually be usable.

**Sequencing is mandatory here.** These weights were trained on MIDOG++, which contains all 18
current ROIs. Any recall number computed with them on 094–548 is contaminated and must not be
reported. **P3 is only measurable after P2**, and even then only on held-out *slides* — note the
split is by slide, and a leave-slide-out split does not guarantee the encoder never saw a
neighbouring ROI from the same case. Say so when reporting.

> **Correction (2026-09-04).** "all 18 current ROIs" is now **23**. The argument is unaffected and
> if anything stronger: MIDOG++ contains all 23, so the contamination warning covers the five ROIs
> added on 2026-09-04 as well.

**Pre-committed criterion:** on held-out ROIs, recall@250 of the FCOS candidate list vs recall@250
of augmented `tm_ccoeff_od`. If FCOS at 250 candidates matches or beats TM at 250, the candidate-
generation question is settled and the blob-vs-TM comparison should be retired, not extended.

### P4. The representation experiment: does the click work in a feature space that has never seen these ROIs? ~2–3 days.

**Decision:** is the click-as-search premise viable at all, or only viable in a contaminated
representation?

This is `2026-09-01-one-click-retrieval-literature.md` §7 re-run without the contamination that
made `2026-09-01-click-ranking-experiment.md`'s positive result uninterpretable. The construction is
already written (`click_rank_embed.py`, `lookalike_auc_unbiased.py`); what is missing is an encoder
that is mitosis-aware and has not seen these images. Two routes:

* **Cheap:** LoRA- or linear-probe-adapt a permissively licensed backbone on MIDOG++ minus the
  domains under test. ResNet-18 @ 64 px is 61 s end-to-end for 20k candidates on this laptop
  *[unverified, lit log §3]*, so this is CPU-feasible.
* **Cleaner:** retrain the FCOS backbone leave-these-domains-out using the MIT training code.

**Pre-committed criterion, inherited from the lit log §7 and kept:** cosine-to-click must beat
chromatin's reading depth on the **worst of 5 seeds**, on every decision-grade held-out ROI. Median
wins do not count — the premise test already found reading depth varying 2.2x with which cell was
clicked, and a ranker that fails on one click in five is a product defect. Drop the seed annotation
from its own positive set (cosine scores it 1.0 by construction).

**This is what falsifies the whole direction.** If a mitosis-trained, uncontaminated representation
*also* shows the click adding nothing over a click-free ranker, then one-click retrieval is not a
viable interaction for this task, and the honest product is a seedless detector with a
confirm/reject UI. That conclusion would be worth reaching in three days rather than three months.

### P5. Equal-coverage comparison instead of equal-*n*. ~1 h, bundled into P1's run.

**Decision:** is S1 (TM beats blob at small budgets) a detection result or a tiling artefact?

Every "matched pool" comparison in this project truncates arms to equal list *length*. At equal
length TM covers 0.88–0.95 of ROI area and blob covers 0.64–0.99 — not the same experiment. Truncate
each arm to equal `coverage_frac` instead (binary-search the depth per arm per ROI).

**Cost correction on my own first estimate:** `tm_variant_sweep.py` writes only aggregate rows —
the candidate pools themselves are never persisted *[verified]* — so this cannot be done as a
re-analysis of the CSVs on disk and needs the matching re-run. Bundle it into P1's pass, where the
pools exist in memory anyway, rather than paying for a separate run.

**Pre-committed criterion:** if TM's recall@250 advantage over blob drops below +0.02 median once
coverage is equalised, S1 is a tiling artefact and the union-generator plan (§5d item 4 of the tm
log) loses its motivation.

### P6. Re-run `--exp2` on the current 17 ROIs. ~70 min. Low priority; do only if the six-method grid is needed.

`tm_variant_sweep.csv` predates the four new ROIs *[verified: it has 13 ROIs, the other two files
have 17]* and predates the `DegenerateMapError` fix, so `TM_SQDIFF_NORMED` still sets `n_matched` in
part of it and truncates every other arm's `@matched` figures. **Do not re-run it to "complete the
grid".** Given E1 and U1 — the method choice does not matter under the ranker we would use — the
six-method comparison has no live decision attached to it. Re-run it only if something downstream
actually needs `TM_CCORR`.

> **Correction (2026-09-04).** This heading's "the current 17 ROIs" is stale, and it matters here
> because P6 is an instruction rather than a result: `images/` now holds **23** ROIs, of which
> **22** are seedable (001.tiff has no mitotic annotations, so it has no seed pool). A re-run today
> would therefore cover 22 ROIs, not 17, at proportionally more than the ~70 min quoted. The
> paragraph below is unchanged and still correct: `tm_variant_sweep.csv` does hold 13 ROIs, and the
> other two files 17.

### Explicitly not worth doing

* **A full 10-ROI NMS-radius sweep.** The existing sweep covers only the two ROIs where TM was
  losing, which is a biased sample, and it buys the tail by paying at the head. Extending it to 10
  ROIs would still be underpowered for the +0.03-recall effect it is chasing (see section 6). Set
  the radius from a stated principle, or leave it.
* **Further operating-point tuning of the correlation search.** Six such moves — score threshold,
  template size, search channel, NMS ordering, bbox tightening, augmentation count — came back
  neutral, for the reason the chromatin log identified: they are monotone moves on a fixed ranking
  function. (Rotation augmentation is the one that later stopped being neutral, and only for the
  ceiling, and only for the unnormalised method. *[unverified, tm log §5c]*)
* **Per-domain tuning of the blob constants** (`min_area`, `max_area`, `tile`, `gray_max`). There is
  no ground truth at inference time to tune them against, so any gain is unrealisable in the product.

---

## 4. Approaches nobody here has pursued

The framing in every log so far is "template matching vs a blob detector as candidate generators".
That framing has three unexamined assumptions.

**(a) It assumes the candidate list should be ~15,000 long.** Both arms produce pools covering
80–95% of the ROI. At that density the ranker is doing all the work and the "generator" is close to
an exhaustive tiling — which is why the disc control matches the click (section 2) and why the
method choice vanishes under chromatin ranking (E1). P3 attacks this directly. A second angle
nobody has tried: **a cascade** — cheap chromatin threshold to a few thousand, then an expensive
per-candidate score on those only. The lit log's latency analysis makes the expensive stage
affordable precisely because the candidate count drops.

**(b) It assumes the click should influence the candidates.** Four possibilities have never been
separated: click influences candidates only; ranking only; both; neither. The evidence now says
**candidates-only is where the click currently acts, and it buys nothing** (section 2), while
ranking-only is where the one positive result in the project lives (the task-trained encoder), and
that result is contaminated. This decomposition should be the axis of the next experiment, not
`cv2.TM_*` constants. `midog_utils` should be split into `propose()` and `rank(click, candidates)`
so the two are separately measurable — the lit log recommended this and it was not done.

**(c) It assumes one click, then a static ranked list.** Untested alternatives, roughly in order of
expected value per unit of effort:

* **Multi-click / iterative relevance feedback.** The click-ranking log found three clicks better on
  every ROI, every encoder and every metric it measured, and that it fixes the one-bad-click-in-five
  failure mode *[unverified, that log's §"What follows" item 3]*. The pathologist is already reading
  the list and marking things; feeding accepted and rejected candidates back as positives and
  negatives is nearly free at the interaction level and turns a one-shot similarity into an online
  classifier. Given that per-seed variance is the documented product defect, this may be worth more
  than any generator change. **Nobody has run it as an experiment on the TM/chromatin path at all.**
* **Negative clicks specifically.** Every control in this project shows look-alikes rank like
  mitoses. A rejected look-alike is a *labelled hard negative on this slide* — exactly the
  information no seedless method can have and the one thing the click demonstrably carries that
  chromatin cannot express. This is the most under-exploited signal in the repo.
* **Field of view.** The template is fixed at 51 px throughout. The click-ranking log found FOV worth
  2–3x at one click and up to 9x at three, and says it was not in any plan *[unverified]*. It is a
  one-parameter sweep and it has never been run on the TM path.
* **Abandoning the ranked list for a coverage-ordered tour.** Goal (1) is "every mitotic figure
  findable". A list ordered by *similarity* revisits the same morphology repeatedly; a list ordered
  to maximise marginal coverage of the ROI's feature space (facility-location / max-marginal-
  relevance selection over the candidate embeddings) reaches the *unusual* mitoses sooner. Given
  that the last 3.5% of recall costs ~10,000 candidates (E2), and that the misses are by definition
  atypical, diversity-aware ordering targets exactly the failure mode. Cheap to test on cached
  scores.
* **Calibrated stopping instead of full recall.** The product goal as written ("every mitotic figure
  findable") may be unachievable at any acceptable budget — E2 says it currently is not. A
  defensible alternative is a calibrated estimate: "you have read 300 candidates; the estimated
  remaining count is 4 ± 2". Mitotic *count* is what the grading actually needs, and a
  capture-recapture or score-distribution estimator can deliver a count without exhaustive reading.
  This is a product-design question that should be put to the pathologist before more engineering.

**(d) Learned components, CPU-only.** The hardware note is a real constraint but a softer one than it
looks: neither proposal generation nor candidate embedding depends on the click, so both belong at
ingest. What must be click-time is one patch embedding plus a cosine against a cached matrix —
milliseconds. So "CPU-only" rules out ViT-g foundation models on this laptop (10–25 h/ROI
*[unverified]*) but does not rule out a small trained CNN, which is the P4 route. The frozen-
pathology-FM path is already refuted by the click-ranking log's own measurements and should not be
retried.

---

## 5. Methodology debt

Five things should change about how work here is validated. The first is not optional.

**M1. Commit, and tie every CSV to a commit.** Nothing from two sessions is in git. This has now
been flagged as a defect in two consecutive audits (`2026-09-01-premise-test-audit.md` M10, then the
tm log's own header) and repeated anyway. Every result script should stamp `git rev-parse HEAD` and
a dirty-tree flag into its output CSV as a column. A result that cannot be tied to code is not a
result.

**M2. Self-written invariants do not substitute for adversarial review, and should stop being
reported as if they do.** 4,257 checks passed; three independent audits then found 15 defects,
including one that invalidated a method outright and one where the invariant meant to catch it was
guarding the wrong quantity. The pattern is that gates catch what their author anticipated. Two
concrete changes: **(a)** stop quoting check counts as evidence of correctness in log headers — they
measure effort, not validity; **(b)** make one adversarial audit by a reader who has not seen the log
a required step before any result is quoted, and budget for it in the experiment's cost.

**M3. Every CSV needs a held-out set and a selection rule.** See P2. Add to this: the four ROIs added
on 2026-09-02 were selected as the densest qualifying image per domain, and that selection rule was
not pre-registered. Any future ROI addition states its rule first.

**M4. A specific trap in the current CSVs: every cell is duplicated 8×.** `tm_variant_stage_b.csv`
and the others emit one row per `budget` level, so each (ROI, seed, arm) appears 8 times with
identical `read_*`, `ceiling`, and `coverage_frac` values. *[verified]* Any row-level test on these
files silently octuples *n*. Medians survive it; counts, *p*-values, and correlations do not. This
is precisely the class of defect the invariant suite was meant to catch and did not. Either emit
budget columns wide, or add an assertion that any statistical routine deduplicates on
(file_name, seed_index, arm) first.

Four more traps in the same CSVs, all found by independent re-derivation and none by the invariant
suite *[verified by audit]*:
* **A live, uncorrected data defect.** `blob_native@matched`'s `coverage_frac` is frozen at the
  seed-0 value in the shipped CSVs — wrong in **23 of 29 cells** in both `tm_variant_sweep.csv` and
  `tm_variant_stage_b.csv` (16 of 29 in `tm_ccoeff_headtohead.csv`), maximum error **0.182**
  (301.tiff seed 3: stored 0.639 against a correct 0.821). The tm log records this as fixed *in the
  code* and recomputed in a side file, but **the main CSVs still carry the wrong column**. Anyone
  reading blob coverage from them gets the bug. Regenerate or annotate the files.
* **`mad_scale` and `map_median` are bit-identical across all six methods within every cell** — they
  are per-(ROI, seed) image statistics broadcast to every row, not per-method map statistics. No
  claim phrased as "method X's MAD" can be sourced from these columns, and the `TM_SQDIFF_NORMED`
  MAD-collapse claim in particular cannot.
* **7 of the sweep's 13 ROIs carry almost no information.** Effective target counts are 1 (505,
  which also delivers only 2 seeds), 3 (350), 3 (351), 5 (406), 8 (002) and 11 (506). All are
  non-decision-grade, so no headline rests on them — but they are two-thirds of the file, and 62
  cells in it have `full_list_recall == 0` outright.
* **The look-alike arms have one seed**, so their `cv` is NaN throughout the seed-variance files:
  there is no seed-stability estimate for any click control. P1 fixes this.

**M5. There is ground-truth leakage in the shared baseline code, and it has never been logged.**
*[verified by independent audit of the four unlogged scripts]* No *ranking function* consumes ground
truth — that part is clean. The leaks are in pool construction and seed placement, and they fall in
three classes that need different responses:

* **(a) A pool constant selected by counting GT mitoses on the evaluation images.**
  `midog_utils/baselines.py:26-28` sets `TISSUE_GRAY_MAX = 220`, justified in its own comment as
  *"validated across all 14 downloaded ROIs: it excludes 0 of the 690 mitotic annotations in them"*,
  and line 38 records global Otsu being **rejected** because it put 31 of 218 mitotic figures on
  301.tiff outside the tissue mask. This is the clearest leak in the repo. It affects every arm
  equally, so it does not bias arm-vs-arm comparisons — but it inflates every **absolute** ceiling
  and recall number, and those are the numbers a product claim would quote. Re-derive the constant
  on held-out ROIs (P2) or from a label-free criterion, and state it in the log.
* **(b) Seed placement conditioned on GT — in two scripts but, importantly, *not* in the main TM
  sweep.** `tm_vs_blob_depth.py:105-107` and `click_darkness_probe.py:59-63` draw the simulated click
  from `y == 1`, i.e. from *blob detections that matched a mitosis*. Three consequences: the click is
  snapped to a blob centroid rather than a pathologist's pixel; mitoses the proposal stage **missed**
  can never be clicked, so the seeded protocol is conditioned on the seedless arm's own success; and
  the click always lands inside the area gate. **I checked `tm_variant_sweep.py` separately and it is
  clean** *[verified]* — `draw_seeds` draws from `gt_mitotic` annotations under agreement tiering and
  a border filter, never from detections. So every number in section 2 is free of this defect, and
  the two probe files are not. Record the distinction rather than tarring all of it.
* **(c) Selection on a prior GT-scored outcome.** `watershed_sweep.py:41` fixes
  `ROIS = ["301","245","246"]`, chosen as "the two regressions plus the largest win" from an earlier
  scored run — so any constant chosen from that sweep is validated on ROIs picked because of how
  that constant behaved. Same defect shape as the NMS-radius sweep covering only the two ROIs where
  TM was losing. **Standing rule: the ROI set is fixed before the sweep, or the sweep does not
  select a constant.**

**M6. Report the product metric, not the proxy.** `read_50` has driven every decision in this
project, and section 2 shows it is measured exactly in the budget range where the click's apparent
advantage lives and the tail's cost is invisible. `recall_at_budget` is already computed in every
CSV and unused in every log. Make the standard reporting unit a recall-vs-budget curve at
25/100/250/500/1000/5000, at ROI level, with the worst-of-5-seeds curve beside the median.

---

## 6. Statistical power — what is and is not worth attempting

With 10 decision-grade ROIs and a two-sided sign test at ROI level (the correct unit, since 5 seeds
within an ROI are not independent): *[verified]*

| wins / 10 | p |
|---:|---:|
| 6 | 0.754 |
| 7 | 0.344 |
| 8 | 0.109 |
| 9 | 0.021 |
| 10 | 0.0020 |

**Only near-unanimity is detectable.** Concretely:

* An effect must be **consistent in direction on 9 of 10 ROIs** to clear p < 0.05. Given that per-ROI
  recall differences between the arms under comparison are +0.02 to +0.08 while between-ROI variation
  in absolute recall spans 0.05 to 0.93, only a mechanism that helps *every* ROI will do this.
* **Anything expected to move recall by less than ~0.05 consistently is not worth running as a
  comparison at this sample size.** That rules out: the NMS-radius sweep at scale, further template-
  size tuning, `TM_CCOEFF` vs `TM_CCORR`, and any second-order operating-point move.
* **Per-cell Wilcoxon tests over 29–50 cells are not a substitute.** They treat 5 seeds per ROI as 5
  replicates and overstate evidence by roughly the design effect. The tm log's own audit found this
  (§8 item 4) and the correction should be standing policy: **ROI-level sign test as the headline,
  per-cell tests as description only.**
* **ROIs are themselves clustered by domain, so the effective *n* is below 10.** In the S1 contrast
  TM wins on 8 ROIs and loses on exactly two — 300 and 301 — which are the *same tumour type on the
  same scanner* (canine cutaneous mast cell, Aperio CS2). *[verified]* So "8 of 10 ROIs" is really
  "every domain except one", and the two losses are one observation, not two. A sign test over 10
  ROIs drawn from 7 domains overstates its own evidence for the same reason a Wilcoxon over 50 cells
  drawn from 10 ROIs does. When P2 lands, report a domain-level as well as an ROI-level test.
  (Worth noting the 300/301 loss is *not* a coverage artefact in the S2 sense: on both, TM has
  **higher** coverage than blob (+0.042, +0.014) and still lower recall, so the blob detector is
  genuinely better on that domain.) *[verified]*
* **The only purchase on power is more ROIs** — which is P2, and which is also the only way to get
  domain generalisation. Doubling to ~24 decision-grade ROIs makes an 8/10-equivalent (19/24) reach
  p ≈ 0.007 and brings +0.05 effects into range. Prioritise *breadth of domain* over count: a second
  ROI from an already-represented tumour type buys much less than a first from a new one.
* One genuinely underpowered contrast worth naming: the **look-alike control is one seed per ROI**.
  It can refute a large specificity effect and cannot measure a small one. If look-alike
  discrimination matters to the product, it needs 5 look-alike seeds per ROI, which is a cheap change
  to the same script.

---

## 7. Summary — the Monday sequence

1. **P0** — commit the working tree (15 min). Precondition for everything.
2. **P1 + P1b + P5, one pass** (~2–3 h total). The augmented run with `disc_*` and 5-seed
   `lookalike_*` arms, a monotonicity-preserving click-conditioned ranker, and equal-coverage
   truncation — all three need the same pools, so build them once. Settles whether the click is
   worth anything in the configuration we would ship, whether it can be made worth something
   cheaply, and whether TM's head advantage is detection or tiling. My prediction is P1 row 3
   ("click not justified in this representation"), P1b negative, and "partly tiling" — but P1's gate
   is built so that row 2, "signal present but too variable", is reported as the multi-click result
   rather than misread as another negative.
3. **P2** — draw held-out ROIs from `datasets_xvalidation.csv` under a pre-registered random rule
   (~2 h). Fixes the deepest methodology hole and is the only route to statistical power.
4. **Housekeeping while downloads run** (~1 h): strike the `tm_vs_blob_depth.csv` 100% column from
   the record (U7), re-derive `TISSUE_GRAY_MAX` without GT (M5 item a), and add `cov^n_gt` beside every
   ceiling in the reporting code (E3).
5. **P3** — FCOS as the candidate generator, measured on held-out ROIs only (~1 day). Interrogates
   the framing that has governed every experiment so far.
6. **P4** — the uncontaminated representation experiment (~2–3 days). The one that decides whether
   one-click retrieval is a viable product at all.

If P1 and P4 both come back negative, the honest conclusion is that the click belongs in the
verification loop rather than the search, and the product is a seedless detector with a
confirm/reject interface plus negative-click hard-negative mining. That is still a product, and it
is better to know in a week than after another month of generator comparisons.
