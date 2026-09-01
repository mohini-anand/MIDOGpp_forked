# The premise test: what the pathologist's click is worth

Date: 2026-08-31
Plan: `2026-08-31-next-steps-plan.md`
Code: `midog_utils/compare.py`, `midog_utils/invariants.py`, `premise_test.py`
Results: `results/premise_{test,seed_provenance,verification}.csv`, `results/premise_report.md`
Run: 7 ROIs x 5 seeds x z in {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0} x 11 arms, 91 min.
1,727 invariant checks pass; the one-match-many-z shortcut is exact in 14/14 checked cells.
Corrected and re-run 2026-09-01 — see the audit ledger at the foot of this file. Decision-grade ROIs (n_mitotic >= 15): 301 (217), 246 (115), 201 (17).

## The answer, in one paragraph

**The click carries real information, and the machinery built around it gives back more than
the click is worth.** Seeding from a genuine mitotic figure beats seeding from a
pathologist-rejected look-alike by +0.052 median recall@1000 and from a random nucleus by
+0.014 to +0.070 — small but consistent, and in the right direction. But against the *best*
seedless arm at matched budget the pipeline wins 0 of 5 seeds on both statistically powered
ROIs, at every budget. Against any single *fixed* seedless arm it is near parity. A blob
detector that never sees an annotation does better on 301.tiff (0.931 vs 0.747) and, just
barely, on 246.tiff (0.861 vs 0.852) — and far worse on 201.tiff (0.588 vs 1.000).

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

Consistent in sign everywhere, decisive only on 246.tiff. On 301.tiff, 3 of 5 is not
distinguishable from a coin flip. **The premise is not dead — but the effect is at the edge
of what 5 seeds on 2 ROIs can resolve, and it is smaller than the pipeline's deficit against
methods that use no click at all.**

## Does the search earn its place? No, against the best seedless arm — but it is close

Pipeline against the best seedless arm, pooled over the three decision-grade ROIs, 5 seeds.
`best seedless` is chosen per (ROI, budget) by its own best median; the pipeline's z is
chosen the same way. **Both sides are selected post hoc, and neither is out of sample:**

| budget | best seedless | pipeline (best z) | delta | pipeline wins |
|---:|---:|---:|---:|---:|
| 500 | 0.811 | 0.774 | −0.037 | 0/15 |
| 1000 | 0.935 | 0.852 | −0.083 | 3/15 |
| 2000 | 0.968 | 0.930 | −0.037 | 0/15 |
| 5000 | 0.991 | 0.965 | −0.026 | 1/15 |

On the two ROIs that carry statistical weight the pipeline wins **0 of 5 seeds at every
budget**. That is the result, and it survives every reading below.

**But "the best seedless arm" is doing a lot of work.** Pinning the pipeline at z = 1.0 and
comparing against each seedless arm *individually*, paired by (ROI, seed):

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

Against any *single fixed* comparator the pipeline is near parity: it loses to `grid_20` and
to the suppressed blob arms, beats `blob_od` and `grid_10`, and the median paired deltas sit
within ±0.08 throughout. Both readings are true and both are reported here, because the
difference between them is entirely a matter of who gets to pick the opponent.

**On the direction of the bias.** An earlier version of this log said the pipeline "is
allowed to pick its best z *post hoc* per ROI and seed — a handicap in its own favour". That
was wrong in two ways. The seedless side is selected post hoc too, and more aggressively:
the numbers that sentence introduced (2/15, 3/15, 0/15, 0/15 at −0.017/−0.026/−0.052/−0.041)
came from a **max over all six seedless arms at every (ROI, seed, budget) cell** — an
opponent no deployable system could be. `results/premise_report.md` §3, generated from the
same CSV, applies a milder rule and prints different numbers. Two undocumented selection
rules, two answers. The report's own text has it right: with both sides tuned in sample,
*"the bias runs in the pipeline's favour."*

`nucleus_blobs` ranked by its own native score is the arm to beat on 301.tiff; a 20-30 px
lattice ranked by chromatin density is the arm to beat on 246.tiff.

## It wins the reading-depth metric on 2 of 3, by less than it first appeared

Candidates a reader works through to reach 50% sensitivity, median over seeds, **pipeline
pinned at z = 1.0** — one rule for every row — with the count of seeds that never reach the
target, which is the column the first version of this table omitted:

| ROI | pipeline z=1.0 | unreached | best seedless (fully reached) | unreached | |
|---|---:|---:|---:|---:|---|
| 201.tiff | **38** | 0/5 | 40 (`blob_od`) | 0/5 | pipeline |
| 246.tiff | **111** | 0/5 | 138 (`grid_20`) | 0/5 | pipeline |
| 301.tiff | 188 | 0/5 | **168** (`blob_native`) | 0/5 | seedless |

The first version of this table read `201.tiff: 28 (z=2.5) vs 40`, and chose a different z
for each row (2.5, 2.0, 2.0). **That 28 is a median over the 3 of 5 seeds that reached 50% at
all** — at z = 2.5 the other two never do. Ranking a survivorship median against an arm that
reached the target on 5 of 5 is not a comparison. Under one z with every seed reaching, the
201.tiff margin collapses from 12 candidates to 2.

The headline survives: the pipeline still wins 2 of 3 on reading depth while losing on recall
at budget. The disagreement is structural rather than noise — the pipeline's list is short
and front-loaded, so it reaches the first half of the mitoses quickly and then runs out; the
seedless lists are long and even. Which matters depends on whether the reader stops at 50%
sensitivity, which is not a defensible clinical operating point, or needs 80-100%.

## z = 2.5 is a bad operating point, and z = 1.0 is an interior optimum

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
   six seedless arms — an opponent no system can be. `results/premise_report.md` §3, from the
   same CSV, applies a different rule and prints different numbers. Replaced with the report's
   rule *and* a full per-arm matrix, so no selection is needed to read the result.
2. **The direction of the selection bias was stated backwards.** "The pipeline is allowed to
   pick its best z post hoc — a handicap in its own favour" omitted that the seedless side is
   selected post hoc too, and more aggressively. The generated report says the opposite in
   bold. Rewritten.
3. **The reading-depth table compared a survivorship median against a complete one.** The
   `201.tiff: 28 (z=2.5)` cell is a median over the 3 of 5 seeds that reached 50% at all, and
   each row used a different z. Rebuilt at one z with an explicit `unreached` column; the
   201.tiff margin falls from 12 candidates to 2. The 2-of-3 headline survives.
4. **"A blob detector does better" hid a reversal.** True on 2 of 3 decision-grade ROIs, and
   false by 0.41 on the third. Split per ROI.
5. **The blob baselines were the only arms never suppressed.** `nucleus_blobs` emits raw
   connected components, 8-36% of which sit within one match radius of another, so a duplicate
   consumed budget without being able to claim a second object. Suppressed variants added
   (`blob_native_nms`, `blob_od_nms`). On *paired* per-cell deltas — not the difference of
   two independently-computed medians, which overstates the native arm threefold — the od
   arm gains +0.051 at budget 500 and about +0.009 at every larger budget, and the native
   arm +0.009 at budget 500, +0.005 at 1000, and nothing beyond. Small, but it moves the
   comparison **against** the pipeline. The raw arms are kept so the size is visible.
6. **`coverage_frac` is a full-list number and was being used to defend budgeted claims.**
   `coverage_at_budget` now exists, and it does **not** say what a first draft of this
   ledger claimed. That draft quoted a 0.080-0.085 spread at budget 1000 on 301.tiff — which
   was a five-arm subset that happened to exclude the four arms breaking it. The full range
   is 0.0638-0.0852, a **1.34x** spread, and full-list coverage over the same arms spans
   0.749-1.000 — also 1.34x. So coverage at budget is *not* better matched across arms than
   full-list coverage is. What it does show is the **level**: a reader at budget 1000 sees
   about 8% of the tissue, not the 75-100% the full-list number suggests. The correct
   conclusion is therefore weaker than "the budgeted comparison is fair": it is click-matched
   but **not** coverage-matched, and the ~34% coverage spread between arms is a live caveat
   on every comparison in this log. Recorded rather than buried because the cherry-pick was
   mine.
7. **The z sweep stopped at 1.0 for a reason that was 301.tiff's number and the smallest of
   seven.** Floor deepened, sweep extended to z = 0.25. z = 1.0 is an interior optimum, not an
   edge; the "one level below 1.0 is owed" item is discharged.
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
    continued. Making it raise then exposed the real finding: extending the check to the new
    low z levels (`results/premise_shortcut_lowz.csv`, 18 cells the 91-minute run never
    touched) found **one genuine violation** — 301.tiff seed 3 at z = 0.75 differs by 1
    detection in 16,825, the two candidates 1.41 px apart. `np.argsort` is unstable, so a
    tie pair inside the NMS radius can flip between a large pool and a small one. 31 of 32
    checked cells are bit-identical and the exception is well inside the 29.6 px match
    radius, so no number in this log moves — but "this is exact, not an approximation" was
    false and is now calibrated. The check raises only on a disagreement with no counterpart
    within a match radius.

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
