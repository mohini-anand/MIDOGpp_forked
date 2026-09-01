# The premise test: what the pathologist's click is worth

Date: 2026-08-31
Plan: `2026-08-31-next-steps-plan.md`
Code: `midog_utils/compare.py`, `midog_utils/invariants.py`, `premise_test.py`
Results: `results/premise_{test,seed_provenance,verification}.csv`, `results/premise_report.md`
Run: 7 ROIs x 5 seeds x z in {1.0, 1.5, 2.0, 2.5, 3.0} x 9 arms, 9.6 min. 1,153 invariant
checks pass. Decision-grade ROIs (n_mitotic >= 15): 301 (217), 246 (115), 201 (17).

## The answer, in one paragraph

**The click carries real information, and the machinery built around it gives back more than
the click is worth.** Seeding from a genuine mitotic figure beats seeding from a
pathologist-rejected look-alike by +0.052 median recall@1000 and from a random nucleus by
+0.014 to +0.070 — small but consistent, and in the right direction. But the whole pipeline
loses to *seedless* baselines by 0.026-0.052 on the same metric, a deficit larger than the
advantage the click buys. A blob detector that never sees an annotation does better.

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

## Does the search earn its place? No, on recall at budget

Pipeline against the best seedless arm, pooled over the three decision-grade ROIs, 5 seeds.
The pipeline is allowed to pick its best z *post hoc* per ROI and seed — a handicap in its
own favour — and still:

| budget | pipeline wins | median delta |
|---:|---:|---:|
| 500 | 2/15 | −0.017 |
| 1000 | 3/15 | −0.026 |
| 2000 | 0/15 | −0.052 |
| 5000 | 0/15 | −0.041 |

`nucleus_blobs` ranked by its own native score is the arm to beat on 301.tiff
(0.802/0.931/0.968/0.986 at the four budgets); a 20-30 px lattice ranked by chromatin density
is the arm to beat on 246.tiff.

## But it wins the reading-depth metric on 2 of 3

Candidates a reader works through to reach 50% sensitivity, median over seeds:

| ROI | pipeline (best z) | best seedless | |
|---|---:|---:|---|
| 201.tiff | **28** (z=2.5) | 40 (`blob_od`) | pipeline |
| 246.tiff | **106** (z=2.0) | 138 (`grid_20`) | pipeline |
| 301.tiff | 176 (z=2.0) | **168** (`blob_native`) | seedless |

The two metrics genuinely disagree, and the disagreement is structural rather than noise:
the pipeline's candidate list is short and front-loaded, so it reaches the first half of the
mitoses quickly and then runs out. The seedless arms are long and even, so they keep going.
Which matters depends on whether the reader stops at 50% sensitivity — which is not a
defensible clinical operating point — or needs 80-100%.

## The committed operating point is bad, and unstable

`read_50` on 301.tiff, per seed:

| arm | s0 | s1 | s2 | s3 | s4 | unreached |
|---|---:|---:|---:|---:|---:|---:|
| **pipeline z=2.5** (committed default) | 247 | — | 734 | 313 | — | **2/5** |
| pipeline z=2.0 | 156 | 203 | 176 | 172 | 349 | 0/5 |
| pipeline z=1.0 | 174 | 195 | 172 | 188 | 203 | 0/5 |
| blob_native (no seed) | 168 | 167 | 168 | 167 | 168 | 0/5 |
| grid_20 (no seed) | 184 | 181 | 184 | 181 | 184 | 0/5 |

At the committed z = 2.5 the pipeline **fails to reach 50% sensitivity on 2 of 5 seeds** and
spans 247-734 on the rest. The seedless arms are deterministic to within 1%. z = 1.0 is both
the best and by far the most stable pipeline setting (172-203). **The committed default
should move to z = 1.0.** This does not touch the 19.4x chromatin-vs-correlation result,
which was a same-candidate-set comparison, but it does mean that result was measured at a
badly chosen operating point.

**z = 1.0 is best on all 12 (decision-grade ROI x budget) cells** — i.e. the frontier is
still improving at the edge of the swept range. The range was capped at 1.0 on the basis of
an earlier probe showing the candidate pool grows only 3% between z=1.0 and z=0.5 with no
metric movement, so this is probably saturation rather than truncation, but it is not
established by this run and one confirming level below 1.0 is owed.

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
