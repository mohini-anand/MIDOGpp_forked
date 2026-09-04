# The premise test: what the pathologist's click is worth

Date: 2026-08-31
Plan: `2026-08-31-next-steps-plan.md`
Code: `midog_utils/compare.py`, `midog_utils/invariants.py`, `premise_test.py`
Results: `results/premise_{test,seed_provenance,verification}.csv`, `results/premise_report.md`
Run: 7 ROIs x 5 seeds x z in {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0} x 11 arms, 91 min.
Invariant checks pass (197 distinct checks logged in `premise_verification.csv`; the
run's own counter reports 1,727 including repeats). The one-match-many-z shortcut is
**not** exact — see `shortcut_probe.py`.
Corrected and re-run 2026-09-01, then corrected again 2026-09-02 after a second audit found
errors in the first round's corrections. Both ledgers are at the foot of this file.
Decision-grade ROIs (n_mitotic >= 15): 301 (217), 246 (115), 201 (17).

## The answer, in one paragraph

**This experiment cannot settle the question it was built to ask, and the reason is
structural rather than statistical.** It varies how candidates are *generated* while holding
the *ranker* fixed at chromatin density for every arm — and on identical candidate sets the
ranker is worth far more than the generator. Swapping chromatin density for the correlation
score costs +0.706, +0.461 and +0.488 recall@1000 on 201, 246 and 301; the entire spread
across the eight generator arms is 0.353, 0.104 and 0.226. So "does the click earn its
place?" is being adjudicated on the smaller term, inside a design where the arms cannot
separate much, at a sample size that cannot resolve small differences: 3 usable images, 5
clicks each drawn **with replacement** (4 of 7 ROIs reuse an annotation), where the best
attainable one-sided p-value from a 5-seed sign test is 0.031.

**What the run does establish, and it is worth having:**

* **Candidate generation barely matters once a good ranker is applied.** Eight structurally
  different generators — a click-seeded template match, a dark-nucleus detector, and lattices
  at four spacings — land within roughly ±0.08 of each other. That is the real finding.
* **The pipeline does not beat the best seedless arm on the densest image.** On 301.tiff it
  wins 0 of 5 seeds at all four budgets. On 246.tiff it wins 1 of 5 at budget 5000 and
  otherwise loses; on 201.tiff it wins 3 of 5 at budget 1000 and ties 3 of 5 at the other
  three. "Loses everywhere" would be wrong.
* **The click carries information** — but the contrast is confounded (see below).
* **z = 2.5 is a bad operating point**, and the frontier is flat below z = 1.0.

**What it does not establish, and the write-up previously implied:** that the search fails to
earn its place. Against a *fixed* comparator the pipeline is near parity; against the best of
eight, chosen after the fact, it is behind. Neither number is out of sample.

**And the number nobody measured.** At budget 1000 the pipeline reads **57.8 false positives
per true positive** on 201.tiff, 9.2 on 246.tiff, 5.2 on 301.tiff. Reading 1,000 candidates
to find 17 mitotic figures is worse than reading the slide. `recall_at_budget` cannot see
this, and precision is the axis the product actually lives or dies on.

## Does the click carry information? Yes — a little

Identical evaluation ground truth, z = 2.0, 5 seeds, same pipeline, only the seed's
provenance changes (`results/premise_seed_provenance.csv`):

| ROI | seeded from a mitotic figure | from a look-alike | from a random nucleus |
|---|---:|---:|---:|
| 246.tiff | **0.809** | 0.757 | 0.748 |
| 301.tiff | **0.710** | 0.659 | 0.668 |

Paired per-seed at budget 1000:

| contrast | 246.tiff | 301.tiff |
|---|---|---|
| mitotic > look-alike | 4/5 seeds, +0.052 | 3/5 seeds, +0.051 |
| mitotic > random nucleus | 5/5 seeds, +0.070 | 3/5 seeds, +0.014 |

**This contrast is confounded, and the confound surfaced only on the second audit.** Arm
(a) is drawn by `seed_selection.pick_seed`, whose first stage is `agreement_pool` — it
prefers annotations that **every rater who saw them called mitotic**. Arms (b) and (c) are
drawn by `premise_test._rejection_sample`, which applies a foreground filter and no agreement
filter. So this is not "mitosis vs look-alike vs nucleus"; it is "*unanimously agreed*
mitosis vs arbitrary look-alike vs arbitrary blob". Seed provenance and seed quality move
together and this design cannot separate them. The tiering is arguably realistic — a
pathologist does click things they are sure of — but "only the seed's provenance changes",
as this log said, is false.

Consistent in sign everywhere, decisive only on 246.tiff. On 301.tiff, 3 of 5 is not
distinguishable from a coin flip: with 5 seeds the best attainable one-sided p from a sign
test is 0.031, so 3 of 5 is p ≈ 0.5 — no evidence either way. **The premise is not dead — but the effect is at the edge
of what 5 seeds on 2 ROIs can resolve, and it is smaller than the pipeline's deficit against
methods that use no click at all.**

## Does the search earn its place? Not against the best seedless arm; near parity against a fixed one

**Read the floor first.** `random_in_tissue` — uniform points in tissue, ranked by the same
chromatin statistic every arm uses — scores 0.502 / 0.567 / 0.641 / 0.710 at the four budgets
on 301.tiff. The pipeline scores 0.747 at budget 1000 and the best seedless arm 0.935. So the
whole argued range sits between 0.57 and 0.94, and the pipeline-vs-seedless gap being debated
is roughly 0.04 of that. The floor was measured on **one ROI at one seed** and is excluded
from every median and win count below. That exclusion is disclosed in the report, but it
means no arm here is credited against the only comparator that has no image evidence at all.

Pipeline against the best seedless arm, chosen per (ROI, budget) by its own best median;
the pipeline's z chosen the same way. **Both sides are selected post hoc; neither is out of
sample.**

| budget | best seedless | pipeline (best z) | delta | pipeline wins |
|---:|---:|---:|---:|---:|
| 500 | 0.811 | 0.774 | −0.037 | 0/15 |
| 1000 | 0.935 | 0.852 | −0.083 | 3/15 |
| 2000 | 0.968 | 0.930 | −0.037 | 0/15 |
| 5000 | 0.991 | 0.965 | −0.026 | 1/15 |

The win column hides ties, and the ties matter. Per ROI, counting both:

| ROI | @500 | @1000 | @2000 | @5000 |
|---|---|---|---|---|
| 301.tiff | 0 win, 0 tie | 0, 0 | 0, 0 | 0, 0 |
| 246.tiff | 0, 0 | 0, 0 | 0, 0 | **1 win**, 0 |
| 201.tiff | 0 win, **3 tie** | **3 win**, 1 tie | 0, **3 tie** | 0, **3 tie** |

So: a clean loss on 301.tiff at every budget, effectively a draw on 201.tiff, and 246.tiff
in between. "Wins 0 of 5 at every budget on both powered ROIs" — which an earlier version of
this log and of `chromatin.py`'s docstring both asserted — is false at 246.tiff @5000.

**Against a single fixed comparator it is near parity.** Pinning the pipeline at z = 1.0 and
pairing against each seedless arm by (ROI, seed) across the decision-grade ROIs:

| seedless arm | wins @500 | wins @1000 | wins @2000 | wins @5000 |
|---|---|---|---|---|
| grid_10 | 10/15 (+0.046) | 10/15 (+0.023) | 4/15 (+0.000) | 1/15 (−0.026) |
| grid_20 | 2/15 (−0.009) | 3/15 (−0.017) | 5/15 (−0.023) | 0/15 (−0.035) |
| grid_30 | 8/15 (+0.009) | 8/15 (+0.014) | 7/15 (+0.000) | 4/15 (−0.009) |
| grid_40 | 11/15 (+0.017) | 3/15 (+0.000) | 4/15 (−0.014) | 4/15 (−0.009) |
| blob_native | 7/15 (+0.000) | 5/15 (−0.009) | 5/15 (−0.026) | 7/15 (+0.000) |
| blob_od | 12/15 (+0.078) | 8/15 (+0.043) | 9/15 (+0.017) | 10/15 (+0.035) |
| blob_native_nms | 5/15 (−0.026) | 5/15 (−0.017) | 5/15 (−0.026) | 7/15 (+0.000) |
| blob_od_nms | 8/15 (+0.026) | 7/15 (+0.000) | 6/15 (+0.000) | 9/15 (+0.017) |

The pipeline loses to `grid_20` and to `blob_native_nms`, beats `blob_od` and `grid_10`, and
is roughly level with the rest. Median paired deltas stay inside ±0.08 throughout. Two
caveats on this table: z = 1.0 was chosen on this same data (and under the extended sweep the
pipeline's own argmax is now z = 0.25 in most cells, so the pin is not even current), and
"which arm is the one to beat" changes with budget — on 301.tiff it is `blob_native_nms` at
500 and 1000 and `blob_native` beyond; on 246.tiff it moves between a lattice and a blob arm
across the four budgets.

**How the earlier numbers were produced.** The first version of this log printed 2/15, 3/15,
0/15, 0/15 at −0.017/−0.026/−0.052/−0.041. Those came from a per-cell max over the six
seedless arms that existed then, taken at every (ROI, seed, budget), with the pipeline
likewise maxed over z. `results/premise_report.md` §3 applies a milder rule to the same data
and prints different numbers. Neither rule was stated. The log's remark that the pipeline's
post-hoc z is "a handicap in its own favour" was *incomplete* rather than backwards — the
report says the same thing — but it omitted that the seedless side is maxed too.

## The reading-depth win does not survive either correction

Candidates read to reach 50% sensitivity, median over seeds, pipeline pinned at z = 1.0,
every cell fully reached:

| ROI | pipeline z=1.0 | best seedless (0/5 unreached) | |
|---|---:|---:|---|
| 201.tiff | 38 | 38 (`blob_od_nms`) | **tie** |
| 246.tiff | **111** | 138 (`grid_20`) | pipeline |
| 301.tiff | 188 | **164** (`blob_native_nms`) | seedless |

**1 win, 1 tie, 1 loss — not 2 of 3.** Two separate errors produced the earlier "2 of 3".
First, the original table read `201.tiff: 28 (z=2.5) vs 40` and picked a different z for each
row; that 28 is a median over the **3 of 5 seeds that reached 50% at all**, against an arm
that reached it on 5 of 5. Second, the corrected table still named the *raw* blob arms as
comparators — `blob_od` at 40 on 201.tiff and `blob_native` at 168 on 301.tiff — after this
same log had added suppressed versions precisely because the raw arms carried duplicates that
waste reading budget. Scored against the arms the correction created, 201.tiff is a tie and
301.tiff a clear loss.

**And at the sensitivity this log itself calls defensible, the pipeline loses 2 of 3.** The
paragraph below has always said 50% "is not a defensible clinical operating point". `read_80`
was in the CSV the whole time, fully reached on both sides in every cell:

| ROI | pipeline z=1.0 | best seedless | |
|---|---:|---:|---|
| 201.tiff | **123** | 128 (`grid_20`) | pipeline, by 5 candidates |
| 246.tiff | 652 | **529** (`blob_native_nms`) | seedless |
| 301.tiff | 1499 | **458** (`blob_native_nms`) | seedless, by 3.3× |

The two metrics still disagree, and the disagreement is structural rather than noise: the
pipeline's list is short and front-loaded, so it reaches the first half of the mitoses quickly
and then runs out; the seedless lists are long and even, so they keep going. That is exactly
why the gap widens as the sensitivity target rises — and why picking 50% flattered the
pipeline. At `read_100` the pipeline is unreached in the large majority of cells.

## z = 2.5 is a bad operating point, and the frontier below z = 1.0 is flat

`read_50` on 301.tiff, per seed:

| arm | s0 | s1 | s2 | s3 | s4 | unreached |
|---|---:|---:|---:|---:|---:|---:|
| pipeline z=2.5 | 247 | — | 734 | 313 | — | **2/5** |
| pipeline z=2.0 | 156 | 203 | 176 | 172 | 349 | 0/5 |
| pipeline z=1.0 | 174 | 195 | 172 | 188 | 203 | 0/5 |
| pipeline z=0.5 | 176 | 195 | 172 | 190 | 211 | 0/5 |
| pipeline z=0.25 (floor) | 177 | 195 | 172 | 190 | 210 | 0/5 |
| blob_native (no seed) | 168 | 167 | 168 | 167 | 168 | 0/5 |
| blob_native_nms (no seed) | 164 | 163 | 164 | 163 | 164 | 0/5 |
| grid_20 (no seed) | 184 | 181 | 184 | 181 | 184 | 0/5 |

At z = 2.5 the pipeline **fails to reach 50% sensitivity on 2 of 5 seeds** and spans 247-734
on the rest. The seedless arms are deterministic to within 1%. z = 1.0 is the best stable
setting (172-203).

**The frontier below z = 1.0 has now been measured, and it has resolved.** An earlier version
of this log said "z = 1.0 is best on all 12 cells — i.e. the frontier is still improving at
the edge of the swept range", and that the sweep was capped at 1.0 because an earlier probe
showed the pool grows only 3% between z = 1.0 and z = 0.5. Both statements needed fixing. The
3% was 301.tiff's number and the **smallest of the seven** — the other six grow 6.4-13.4%.
And the level was free to measure: at an extraction floor of `med + 0.5·mad`, a z = 0.5 arm
*is* the whole extracted pool. This run drops the floor to `med + 0.25·mad` and sweeps
z ∈ {0.25, 0.5, 0.75, 1.0, …}. Median recall@1000:

| ROI | z=0.25 | z=0.5 | z=0.75 | **z=1.0** | z=1.5 | z=2.0 | z=2.5 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 201.tiff | 1.000 | 1.000 | 1.000 | **1.000** | 0.941 | 0.824 | 0.529 |
| 246.tiff | 0.843 | 0.843 | 0.843 | **0.852** | 0.852 | 0.809 | 0.730 |
| 301.tiff | 0.742 | 0.742 | 0.742 | **0.747** | 0.742 | 0.710 | 0.525 |

At budget 1000 recall does not improve below z = 1.0 on any decision-grade ROI, and
`read_50` improves in 0 of 60 cells (301.tiff: 188 at z = 1.0, 190 at the floor). **But it
is a plateau, not an interior optimum.** At larger budgets the argmax slides below 1.0:
pooled decision-grade medians are 0.9043 at z ≤ 0.75 against 0.8957 at z = 1.0 for budget
2000, and 0.9565 against 0.9478 at budget 5000; 12 of 60 (ROI, seed, budget) cells improve
below z = 1.0. The report's own §4 finds z = 0.25 best or joint-best on 22 of 28 rows.

What the extension actually bought is narrower than "the frontier resolved", and worth
stating exactly: every difference across z ∈ [0.25, 1.0] is ≤ 0.01 recall, so **the choice
inside that range does not matter**, and the sweep now terminates *at the extraction floor*
rather than at an arbitrary cut — an edge argmax at the floor is a different and much weaker
objection than an edge argmax at a number someone picked. z = 1.0 remains a defensible
operating point because the plateau is flat, not because it is optimal. The "one confirming
level below 1.0 is owed" item is discharged; "where the frontier stops improving" is
answered only in the sense that it stops mattering.

**On "the committed default should move to z = 1.0".** That sentence, in the first version of
this log, presumed a default that does not exist. There is no z-based threshold anywhere in
`midog_utils/`: `FSConfig.score_threshold = 0.5` is in raw `TM_CCOEFF_NORMED` units,
`extract_peaks` carries its own `0.5`, and `compare.Arm.z` is a provenance label with no
default. z exists only as a constant inside the experiment scripts
(`od_experiment.py`, `od_seed_sweep.py`, `od_workload_ab.py`, `premise_test.py`). Adopting
z = 1.0 is therefore either a change of *units* for `FSConfig.score_threshold` — machinery
nobody has written — or an edit to four script constants. **This is not a done deal and the
recommendation has to say which.**

A related trap, not yet fixed: `od_experiment.py` sets its floor at 0.25 in *correlation*
units while `premise_test.py` sets 0.25 in *MAD* units. In `results/od_experiment_arms.csv`
the z = 1.0 arms are `floor_limited` on **6 of 7 ROIs**, so those rows measure the
correlation floor, not z = 1.0. The z = 1.0 evidence rests on this file alone, and the two
must never be pooled.

None of this touches the chromatin-vs-correlation result, which was a same-candidate-set
comparison — but that result was measured at a badly chosen operating point.

## Two specification errors the implementation caught

Both were in the brief, both are real, and one found a genuine data defect:

1. **"Distinct seed_index gives distinct seed_ann_id" is not an invariant.** Seed selection
   draws uniformly *with replacement*, so collisions are expected on correct runs. Verified
   against committed code predating this work: `results/od_seed_sweep.csv` shows 4 of 7 ROIs
   with only 4 distinct annotations across 5 seed indices, and 350.tiff has 3 candidates in
   its pool so collisions are certain. Replaced with the property it was a proxy for —
   distinct seed indices use distinct RNG *streams* — plus a reported `n_distinct_ann`.
2. **`tissue_mask` excludes 1 of 726 annotations**, not 0: 506.tiff ann 24249, **category 2**,
   at x = 13 px from the ROI edge, click pixel grey 232 against the 220 cut. Confirmed
   independently. **0 of 391 mitotic annotations are excluded**, which is exactly what
   `baselines.tissue_mask`'s docstring claims, so the check is now fatal on category 1 and
   records category-2 exclusions rather than aborting.

## A confound in the grid sweep, reported not hidden

NMS only removes lattice points when the step is below the match radius, and that radius
crosses 30 px between ROIs (29.6-33.1). So `grid_40` is always the raw lattice, and
`grid_30` is suppressed on six ROIs but *not* on 301.tiff — 34,615 candidates there against
step-20's 15,736. The step sweep therefore confounds spacing with whether suppression bites
at all. Per-ROI candidate counts are reported so this is visible.

Also: `coverage_frac > 0.5` in 413 of 561 cells, so full-list recall is largely tiling
geometry across this run. The budgeted columns — which carry every conclusion above — are
unaffected.

## The metric ignores the thing the product is sold on

Every comparison above is `recall_at_budget`: of the mitotic figures present, what fraction
appears in the top K. Nothing here measures **precision**, and precision is what a reader
experiences. At budget 1000 with the pipeline at z = 1.0:

| ROI | mitotic figures | found in top 1000 | false positives per true positive |
|---|---:|---:|---:|
| 201.tiff | 17 | 17 | **57.8** |
| 246.tiff | 115 | 98 | 9.2 |
| 301.tiff | 217 | 162 | 5.2 |

On 201.tiff the tool asks a pathologist to reject 983 non-mitoses to confirm 17 — plainly
worse than reading the slide. Recall at budget is blind to this by construction, and every
arm in this experiment shares the blindness, so the comparison between arms is unaffected
while the comparison between *the tool and doing nothing* never happens.

This also exposes an assumption the design never states: **that every candidate costs the
reader the same.** A lattice point landing on stroma is dismissed at a glance; a
high-scoring template match on a dark, mitosis-like nucleus is not. Matched candidate budget
is not matched reader time, and reader time is the axis the product competes on. A
cost-weighted budget would likely move the arms relative to each other, in an unknown
direction.

For a mitotic-count product the natural primary metric is not recall at an invented budget
but **count agreement** — predicted versus true count over a fixed area, and how often the
error crosses a grading boundary. That metric was never computed.

## What this design cannot settle, and why

Stating the limits together, because individually each looks survivable and collectively
they are not:

* **The wrong term is varied.** Every arm shares the chromatin ranker; only generation
  differs. The ranker is worth 2–4× the entire spread across generators.
* **n = 3, not 15.** The five seeds per ROI are five clicks on the *same* image scored
  against the *same* ground truth. They share all image-level variance, and image is where
  the variance lives. Pooling them into "wins out of 15" is pseudoreplication; the effective
  sample size is the number of images.
* **The clicks are not even five.** `pick_seed` draws with replacement, and 4 of 7 ROIs
  reuse an annotation. On 201.tiff seeds 0 and 1 are the same click — and that is the ROI
  carrying the pipeline's only decision-grade win.
* **Nothing is out of sample.** Eight z levels, four lattice spacings and the ranking
  statistic itself were all chosen on these same seven images.
* **Two of three decision-grade ROIs are near ceiling** at the budgets under discussion, so
  the dynamic range is thin exactly where the argument is made.

The cheapest repair costs no compute: choose z and the lattice step on the four
non-decision-grade ROIs and report the three decision-grade ones as held out. Every row
already exists in `results/premise_test.csv`. The second cheapest is to promote
`random_in_tissue` from a single sanity value to a full arm across all ROIs, seeds and
budgets — it is seconds of work and it recalibrates every number in this log.

## What follows

* **Move the default to z = 1.0** and re-state the chromatin-vs-correlation result at it.
* **The seedless baselines are now the thing to beat**, and one of them (`blob_native`) needs
  no seed, no template, no correlation, and 7-11 s. Any further work on the correlation
  search has to clear that bar first.
* **The click's +0.05 is real and small.** If the one-click premise is to be rescued, the
  route is not a better candidate generator — it is using the click for something the
  seedless arms cannot do, which at present nothing in the pipeline does.
* Texture (GLCM/LBP) remains open and is now better motivated: chromatin density and the blob
  component-mean are stuck at AUC 0.69-0.77 against look-alikes and are statistically
  indistinguishable from each other, with 105-123 look-alikes available to measure against.

## What the 2026-09-01 audit changed

An independent audit re-derived every number from the raw data and attacked each claim, then
re-tested its own findings before acting on them. Five of the audit's own ten findings were
wrong or overstated and were retracted before any code changed; that ledger is in
`2026-09-01-premise-test-audit.md` Part 3a. Corrections forced here, all applied above:

1. **The headline win-count table was not reproducible from its own report generator.** The
   printed 2/15, 3/15, 0/15, 0/15 (−0.017/−0.026/−0.052/−0.041) came from a per-cell max over
   the six seedless arms that existed then, with the pipeline likewise maxed over z.
   `results/premise_report.md` §3, from the same CSV, applies a milder rule and prints
   different numbers. Neither rule was stated anywhere. (A best-achievable-vs-best-achievable
   comparison is defensible on its own terms — the objection is the silent divergence between
   two rules, not that the rule is indefensible.) Replaced with the report's rule *and* a full
   per-arm matrix, so no selection is needed to read the result.
2. **The selection bias was described incompletely.** "The pipeline is allowed to pick its
   best z post hoc — a handicap in its own favour" omitted that the seedless side is maxed
   too. *(A first version of this ledger claimed the log and the generated report contradicted
   each other on the direction of the bias. They do not — the report's "the bias runs in the
   pipeline's favour" says the same thing. That claim was wrong and is withdrawn.)*
3. **The reading-depth table compared a survivorship median against a complete one.** The
   `201.tiff: 28 (z=2.5)` cell is a median over the 3 of 5 seeds that reached 50% at all, and
   each row used a different z. Rebuilt at one z with an explicit `unreached` column. **The
   2-of-3 headline does not survive**, for two further reasons found on the second audit: the
   rebuilt table still scored against the *raw* blob arms this same ledger had just
   disqualified (item 5), which turns 201.tiff into a tie and 301.tiff into a loss; and at
   `read_80`, the sensitivity this log itself calls defensible, the pipeline wins 1 of 3.
4. **"A blob detector does better" hid a reversal.** True on 2 of 3 decision-grade ROIs, and
   false by 0.41 on the third. Split per ROI.
5. **The blob baselines were the only arms never suppressed.** `nucleus_blobs` emits raw
   connected components, 16.4-70.5% of which sit within one match radius of another, so a
   duplicate consumed budget without being able to claim a second object. (An earlier version
   of this line said "8-36%" — that is the *suppression* rate, a different quantity.)
   Suppressed variants added
   (`blob_native_nms`, `blob_od_nms`). On *paired* per-cell deltas — not the difference of
   two independently-computed medians, which overstates the native arm threefold — the od
   arm gains +0.051 at budget 500 and about +0.009 at every larger budget, and the native
   arm +0.009 at budget 500, +0.005 at 1000, and nothing beyond — **on the three
   decision-grade ROIs only**. Pooled over all seven the median paired delta is exactly
   0.0000 at every budget for both arms, since the four small ROIs sit at ceiling. Small,
   and restricted, but it moves the comparison **against** the pipeline. The raw arms are
   kept so the size is visible. Their duplication is larger than first reported:
   **16.4–70.5%** of blobs have another blob inside one match radius. The "8–36%" quoted
   earlier is the *suppression* rate — a different quantity, and the two were conflated.
6. **`coverage_frac` is a full-list number and was being used to defend budgeted claims.**
   `coverage_at_budget` now exists, and it does **not** say what a first draft of this
   ledger claimed. That draft quoted a 0.080-0.085 spread at budget 1000 on 301.tiff — which
   was a five-arm subset that happened to exclude the four arms breaking it. The full range
   is 0.0638-0.0852, a **1.34x** spread, against a full-list spread of the same order. So
   coverage at budget is *not* better matched across arms than full-list coverage is. What it
   does show is the **level**: a reader at budget 1000 sees about 6-9% of the **ROI** — the
   statistic is a fraction of the ROI, not of tissue, which a first draft also got wrong —
   where full-list coverage runs 0.715-1.000 for the seedless arms and 0.816-0.951 for the
   pipeline. The correct
   conclusion is therefore weaker than "the budgeted comparison is fair": it is click-matched
   but **not** coverage-matched, and the ~34% coverage spread between arms is a live caveat
   on every comparison in this log. Recorded rather than buried because the cherry-pick was
   mine.
7. **The z sweep stopped at 1.0 for a reason that was 301.tiff's number and the smallest of
   seven.** Floor deepened, sweep extended to z = 0.25. The frontier below z = 1.0 is a flat
   plateau — differences within 0.01 at median level — so the choice inside z ∈ [0.25, 1.0]
   does not matter. It is **not** an interior optimum, which an earlier draft of this ledger
   and the section heading above both claimed while the section body said the opposite: at
   budgets 2000 and 5000 the argmax slides to the floor under one pooling rule, though not
   under the report's §3 rule. What the extension bought is that the sweep now terminates at
   the extraction floor rather than at a number someone chose.
8. **"The committed default should move to z = 1.0" presumed a default that does not exist.**
   Rewritten to say what would actually have to change.
9. **The seed-provenance control was asymmetric.** Arms (a) and (b) had their seed annotation
   removed from the ground truth; arm (c)'s blob-centroid seed did not, while its 5 px
   exclusion zone still applied. Now all three score against the union-excluded ground truth.
   **This changed nothing measurable** — the paired contrasts are identical to three decimals
   before and after (4/5 +0.052 and 5/5 +0.070 on 246.tiff; 3/5 +0.051 and 3/5 +0.014 on
   301.tiff). Reported because the asymmetry was real, not because it mattered.
10. **Four retracted claims were still live in committed code** (`chromatin.py` ×3,
    `od_experiment.py` ×1), plus a fifth found later: `grid_lattice`'s docstring quoted recall
    and read-depth figures measured on a *different* lattice. All corrected.
11. **`verify_shortcut` never asserted anything, and the property it checks is not exact.**
    The module docstring claimed it "asserts it directly"; the code printed to stderr and
    continued. Making it raise exposed the real finding. A widened probe
    (`shortcut_probe.py`, 40 cells against the main run's 14) finds **38 of 40 bit-identical
    and two genuine violations** — 301.tiff seed 3 at z = 0.75 *and at z = 1.0*, each
    differing by 2 detections out of ~16,800, the disputed candidates 1.41 px apart.
    `np.argsort` is unstable, so a tie pair inside the NMS radius can flip between a large
    pool and a small one. The z = 1.0 case matters most: that is the operating point this
    log recommends, and the main run checked it at seed 0 only, which is clean. No number
    here moves — 1.41 px against a 29.6 px match radius is the same object — but "this is
    exact, not an approximation" was false. The criterion now requires the two lists to be
    the same **length** and every disputed point to have a counterpart within
    `peak_min_distance` (7 px), not merely within the match radius: two points 25 px apart
    are inside a 29.6 px radius yet can match different objects.

12. **"5 seeds" is 5 RNG streams, not 5 distinct clicks.** `pick_seed` draws with
    replacement, and 4 of the 7 ROIs reuse an annotation: 201 (seeds 0 and 1 are the same
    annotation), 350, 405 and 506. This lands inside the only positive decision-grade result
    in the report — 201.tiff at budget 1000, "3 of 5 wins" — where seeds 0 and 1 are the same
    click and both are wins. That result is 2 distinct clicks out of 4, not 3 of 5. The log
    already noted that seeds are drawn with replacement; it did not notice that the
    duplication concentrates in the one cell carrying a positive finding.

Claims that survived unchanged: the click carries information, in sign and size, on both
powered ROIs (identical to four decimals before and after the ground-truth fix); the pipeline
loses to the best seedless arm at matched budget on both, winning 0 of 5 seeds at every
budget; the two metrics genuinely disagree; z = 2.5 is unstable and reaches 50% on only 3 of
5 seeds on 301.tiff.

**Not editable in place:** commit `3117402`'s message reproduces the pre-correction figures
("even when allowed to pick its best z post hoc"; "201: 28 vs 40") and attributes the z = 2.5
default to `7574b09`, which added a channel-comparison slide deck and has nothing to do with
z. History is not rewritten; this section is the correction of record.

## What the 2026-09-02 audit changed

Five independent agents re-checked the 2026-09-01 corrections above against the raw data.
The arithmetic held — the per-arm matrix reproduces in all 32 cells, the pooled table, the
z-frontier table, the provenance contrasts and the 1,960-row identity check all verify
exactly. What did not hold was a set of claims *about* that arithmetic. Corrections applied,
beyond the amendments already folded into the ledger above:

1. **"Wins 0 of 5 seeds at every budget on both powered ROIs" is false.** 246.tiff at budget
   5000 is 1 of 5. The same sentence had been copied into `midog_utils/chromatin.py`'s
   docstring; both are fixed. Ties were also suppressed: on 201.tiff the pipeline ties the
   best seedless arm on 3 of 5 seeds at three of four budgets, which reads as a draw, not a
   defeat.
2. **The reading-depth conclusion is gone**, twice over — see ledger item 3 above.
3. **The blob-suppression gain is decision-grade only.** Over all seven ROIs the median
   paired delta is exactly 0.0000 at every budget for both arms. The restriction was never
   stated.
4. **The seed-provenance contrast is confounded by rater agreement** (arm (a) is drawn from
   unanimous annotations, arms (b) and (c) are not). This lands inside the one claim the
   first audit certified as surviving unchanged.
5. **Precision was never measured**, and it is the product's actual constraint: 57.8 false
   positives per true positive on 201.tiff at budget 1000. New section added.
6. **The experiment varies candidate generation while holding the ranker fixed**, and the
   ranker is worth 2–4× the entire spread across generators. New section added.
7. **The shortcut is inexact at z = 1.0, not only at z = 0.75.** A widened probe (40 cells
   against the main run's 14) found a second violation, at the operating point this log
   recommends. The benign criterion was also too permissive — it never compared list lengths
   and allowed displacement up to the full match radius; it now requires equal length and a
   counterpart within `peak_min_distance`.
8. **`results/premise_verification.csv` predates the code that would produce it.**
   `verify_shortcut`'s return schema changed after the 91-minute run, so that file lacks the
   `benign` / `n_disagree` / `worst_disagreement_px` columns. Its 14 `one_match_many_z` rows
   are all `identical`, so the criterion change does not affect them, but the file and the
   function no longer correspond. `results/premise_shortcut_probe.csv` is the current
   evidence. This is the same code-outlives-its-artifact defect the audit was written to
   catch, committed by the audit.
9. **Assorted stale prose**, all from fixes applied in one place and not another: the section
   heading and a ledger item still said "interior optimum" while the body said the opposite;
   `coverage_at_budget` was added, declared "the one to read", and never printed in the
   report; the invariant-count and coverage figures in
   `2026-09-01-premise-test-audit.md` Part 2 still described the pre-change run; and the
   five-arm coverage cherry-pick that item 6 retracts was still present in
   `2026-08-31-next-steps-plan.md`, added by the same commit that retracted it.

**The pattern is the finding.** Nearly every correction in the first round was directionally
right and quantitatively overstated, and three of them reintroduced in a new place the exact
defect they were fixing. The response has been to weaken the conclusions to what the evidence
supports rather than to patch the numbers one at a time — hence the two new sections above on
what this design cannot settle.
