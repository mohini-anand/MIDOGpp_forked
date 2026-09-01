# Next steps: the premise test

Date: 2026-08-31
Follows: `2026-08-31-chromatin-density-rerank.md` (commit 7c3af93)
Status: **converged plan after adversarial review. Nothing in Step 1+ has been run.**

This replaces a first draft whose Phase 1 would have returned a maximally negative verdict
for the wrong reason. The review that caught it ran its own cut-down experiment; its
working files are at `scratchpad/plancheck/` (`fiveseed.csv` is the load-bearing one) and
every number below was re-derived from them independently.

---

## What the review found, and what it means

A **regular 20 px lattice with no seed and no detector**, NMS'd at the match radius and
ranked by chromatin density, is competitive with the whole pipeline. Medians over 5 seeds,
2 decision-grade ROIs, matched candidate budget:

**301.tiff (n=217)** — recall at budget / read-to-50%

| arm | @500 | @1000 | @2000 | @5000 | read-50 |
|---|---:|---:|---:|---:|---:|
| `nucleus_blobs`, native score — **no seed** | **0.802** | **0.931** | **0.968** | **0.986** | **167** |
| grid @20px + chromatin — **no seed, no detector** | 0.700 | 0.816 | 0.899 | 0.972 | 177 |
| `find_and_suppress` z>=2.0 + chromatin | 0.641 | 0.710 | 0.793 | 0.853 | 176 |
| `find_and_suppress` **z>=2.5** (current default) | 0.470 | 0.525 | 0.544 | 0.544 | 313 |

**246.tiff (n=115)**

| arm | @500 | @1000 | @2000 | @5000 | read-50 |
|---|---:|---:|---:|---:|---:|
| grid @20px + chromatin — no seed, no detector | **0.783** | **0.887** | 0.930 | **0.983** | 144 |
| `find_and_suppress` z>=1.0 + chromatin | 0.774 | 0.852 | 0.922 | 0.965 | 111 |
| `nucleus_blobs` native — no seed | 0.765 | 0.861 | **0.948** | 0.948 | 193 |
| `find_and_suppress` z>=2.0 + chromatin | 0.757 | 0.809 | 0.878 | 0.913 | **106** |
| `find_and_suppress` **z>=2.5** (current default) | 0.704 | 0.730 | 0.774 | 0.783 | 115 |

Coverage is matched across arms (0.09-0.11 on 301 at budget 1258), so this is not a tiling
artefact. Three consequences:

1. **z = 2.5 — the operating point we just committed — is the worst configuration in the
   table on both ROIs.** The 19.4x ranker improvement stands (it is a same-candidate-set
   comparison and nothing here touches it), but the *absolute* operating point it was
   measured at is badly chosen. z = 2.0 lifts 301's recall@1000 from 0.525 to 0.710 and its
   read-50 from 313 to 176.
2. **The seeded pipeline does not beat seedless baselines on recall at budget** on either
   ROI, at any z tested.
3. **The pipeline converges to the grid as z falls**, mechanically: at low z the candidate
   set is "every local correlation maximum", which is everywhere.

### Where the review overreached — conceded on re-measurement

It concluded "a grid matches or beats the pipeline on both ROIs." On **read-to-50%** that is
false on 246.tiff: the pipeline wins at every z (106-115) against the grid's 144 and the
blob detector's 193 — a 26-36% advantage, stable across 5 seeds. On 301 the tuned pipeline
(176) ties the grid (177) and slightly loses to blob-native (167).

Per-seed, the pipeline beats the grid on 246's read-50 in **5 of 5 seeds** at both z=2.0
(104/121/121/88/106 against 144/153/144/153/144) and z=1.0. A clean sweep, not noise. The
defensible statement is narrower and still damning: **on recall at budget the click buys
nothing; on reading depth to 50% sensitivity it buys a real advantage on the dense
lymphosarcoma ROI and nothing on the dense mast-cell ROI.** Not "the pipeline loses."

### The finding neither of us was looking for: variance

301.tiff `read_50` per seed:

| arm | s0 | s1 | s2 | s3 | s4 |
|---|---:|---:|---:|---:|---:|
| blob_native (no seed) | 167 | 167 | 167 | 167 | 167 |
| grid step 20 (no seed) | 177 | 172 | 177 | 172 | 177 |
| pipeline z>=2.0 | 156 | 203 | 176 | 172 | **349** |
| pipeline z>=2.5 (committed) | 247 | **unreached** | **734** | 313 | **unreached** |

The seedless arms are deterministic to within 3%. The pipeline varies 2.2x at its best z and
**fails to reach 50% sensitivity on 2 of 5 seeds at the committed z**. No metric in either
draft captured this. A one-click tool whose reading burden depends 2x on which cell the
pathologist happens to click is a product problem independent of its median, so per-seed
spread of read-50 is now a first-class reported column.

### The grid is a tunable knob, not a floor

Refining the lattice makes it **worse**: on 301.tiff step 20 gives recall@budget 0.825 and
read-50 177, step 10 gives 0.765 and 356. A denser lattice finds a better-optimised maximum
of the chromatin field inside each match-radius disc, and that field's true maxima are dense
stromal and nuclear clumps rather than mitoses -- optimising the statistic harder moves away
from the target. The step must therefore be swept (10/20/30/40), not fixed, or the "baseline"
is a tuned hyperparameter. The 301 conclusion survives -- step 20 beats every seeded arm --
but it is one point on an uncharacterised curve.

### One finding neither of us anticipated

`nucleus_blobs` ranked by its **own native score** beats `nucleus_blobs` + chromatin on both
ROIs at every budget (301 @1000: 0.931 vs 0.760). The chromatin re-rank helps the
correlation pipeline and *hurts* the blob detector. Plausible reason: on the blob candidate
set the Otsu component *is* the candidate, so the component-mean statistic has none of the
9-12% segmentation-failure problem that made the window statistic the better choice for
correlation peaks. This means **normalising all arms to one ranker biases the comparison in
favour of the pipeline — and the pipeline still loses on recall at budget.** State it that way.

### Corrections this forces to the committed log

Both are in `2026-08-31-chromatin-density-rerank.md` and are being fixed:

* *"Look-alikes are now the majority problem"* — they are **23.4%** of top-K false positives
  (46 of 197), up from 0.58% on the full list: a 40x enrichment, not a majority. Perfect
  look-alike suppression moves FPs 197 -> 151.
* *"not enough to separate mitosis from mimic"* — too harsh, but the first correction to it
  was too generous. The AUCs in that log are measured on the pipeline's own candidate set,
  which captures only 29-48% of annotated look-alikes — a selection correlated with the
  feature being tested. Re-measured on the `nucleus_blobs` set (98-100% of look-alikes,
  97-100% of mitoses), the bias is worth 0.06-0.10: **246.tiff 0.767 [0.706, 0.829]** on
  109 v 123, **301.tiff 0.692 [0.634, 0.750]** on 215 v 105, against the biased 0.831 and
  0.797. Both exclude 0.5 decisively. The defensible sentence is *0.69-0.77 against 0.94-0.99
  for TP-vs-unannotated — a real but much weaker signal*. Four of the log's seven rows (405,
  002, 506, 350) rest on 2-11 true positives with CIs spanning 0.5 and are noise.
* A hypothesis this kills: the darkest-decile window statistic is **not** doing something
  cleverer than a component mean on the hard class (0.692 vs 0.685 on 301, indistinguishable).
  Its advantage is confined to the easy class — which is exactly why `nucleus_blobs` ranked by
  its own native score beats the same detector re-ranked by chromatin density.

---

## Step 0 — the comparison protocol

`midog_utils/compare.py`: one entry point, arms declared as data, returns a tidy frame.

1. **Never cite a stored CSV as a comparison baseline.** Every arm is re-run from current
   code in the same process.
2. **Hold the seed fixed across arms**, and hold the excluded seed annotation identical so
   the evaluation set is the same.
3. **Rankers are compared on an identical candidate set**; **generators at matched budget.**
4. **recall@budget is the primary metric.** read-to-50/80/100% is secondary: it reads from
   the top, so truncation barely moves it (every arm's read-50 is 106-193 on lists of
   1,258-21,058), which makes it an inconsistent control for generator comparisons.
5. **>= 5 seeds per ROI, always.** Median and IQR and win count; never a point estimate.
   Runs where no arm reaches the target are reported separately, never counted as wins.
6. **Every comparison carries a seedless arm** — grid and `nucleus_blobs`-native. Their
   absence is why this project reached 2026-08-31 without knowing the above.
7. **Every constant that crosses ROIs must be per-image or in physical units.** This one
   rule catches the 0.5 score threshold, the old fixed 25 px NMS radius, a global extraction
   floor, and GLCM offsets specified in pixels.
8. **Mandatory per-arm diagnostics:** `floor_limited`, `coverage_frac`, largest tie block in
   the ranking key, nan rate.
9. **Permutation check on any new ranking key** — re-rank tied blocks 500x, report the range
   of the headline metric. This is what caught the hematoxylin ceiling.
10. **"Can this add information?" gate.** Every rejected experiment in this repo was a
    monotone operating-point move on a fixed ranking function. If a proposed change is
    monotone in the existing key, predict it analytically and skip the run.
11. **Invariant assertions, not a regression test.** A test that re-runs seed 0 and asserts
    committed numbers cannot catch the 2026-08-25 bug class *by construction* — those bugs
    were present when the baseline was committed, so it would lock them in. Assert instead:
    no arm's list length equals a configured cap; `tissue_mask` excludes 0 of the GT
    annotations; distinct `seed_index` gives distinct `seed_ann_id`; `nms_radius ==
    evaluate.radius_px(mpp)` for that image.

---

## Step 1 — the premise test (~15 min, one script)

Merges the draft's Phases 1 and 2, which had a dependency error: Phase 1 ran the pipeline at
z=2.5 while Phase 2 existed to choose z, so the arm-A verdict was a function of an unchosen
parameter.

**One match per (ROI, seed)** at a per-image floor `med + 0.5*mad`; every z derived by
filtering that one pool — the "one match, many arms" shortcut `od_experiment.py` already
proves exact. A global constant floor is forbidden by Rule 7 and is also wrong on the
numbers: 246.tiff needs <= 0.0968 for a clean z=1.0 on its worst seed, against the draft's
proposed 0.12.

**Arms, at every z in {1.0, 1.5, 2.0, 2.5, 3.0}:**

| | seeded? | candidates | ranker |
|---|---|---|---|
| pipeline | yes | `find_and_suppress` at z | chromatin |
| pipeline-score | yes | same | correlation score |
| blob-native | **no** | `nucleus_blobs` | its own score |
| blob-chromatin | **no** | `nucleus_blobs` | chromatin |
| grid | **no** | 20 px lattice, NMS'd at match radius | chromatin |
| random | **no** | `random_in_tissue` | chromatin |

**Plus the seed-provenance control**, which is the experiment that can *rescue* the premise
rather than merely confirm it is beaten. At the best z, matched budget, 5 seeds: seed from
(a) a category-1 mitotic figure, (b) a category-2 pathologist-rejected look-alike, (c) a
random `nucleus_blobs` component. If (b) and (c) match (a), the click carries no information
and no amount of generator tuning will help.

5 seeds, all 7 ROIs. Report recall@{500, 1000, 2000, 5000}, read-50/80/100, `coverage_frac`,
`floor_limited`, and look-alikes-in-list per cell.

**Decision rules, all decidable:**

* *Does the click carry information?* Compare (a) against (b)/(c) on recall@budget, median
  over 5 seeds x 2 decision-grade ROIs. If the CIs overlap, the premise is dead — say so.
* *Does the search earn its place?* Pipeline vs the best seedless arm at matched budget. It
  must win on the **median over ROIs with n >= 15** (301, 246, 201 only), not on a
  small-n ROI.
* *What z?* Report the recall-vs-workload frontier. Bind on a stated sensitivity target —
  **80%**, not 50%, since 50% is already met at 1.7-3.3 candidates per mitosis on the three
  powered ROIs and does not discriminate. Add a **coverage ceiling of 0.5**: past that,
  `evaluate.py`'s own docstring says full-list recall is geometry, not evidence, and the
  grid arm reaches 0.95-0.99.
* *Is it stable enough to ship?* Per-seed IQR and unreachable-count on read-50, against the
  seedless arms' ~3%. A median that hides a 2x spread is not a result.

**Trims, all justified by measurement:** `random_in_tissue` drops to a single sanity row (it
is the grid plus sampling noise; Rule 10 forbids spending compute on a predictable outcome);
the seedless arms are computed once per ROI and evaluated against each seed's `gt_eval`
rather than recomputed per z (an 8x saving on the two slowest arms); and z stops at 1.0,
below which the pool grows 3% and nothing moves. `floor_limited` is false by construction
once the floor is per-image in the same units as the threshold — keep the column as a bug
detector, not as a filter.

---

## Step 2 — widen the evidence base (promoted from "parallel")

Only 301.tiff (n=217) and 246.tiff (n=115) carry statistical weight; 201 (n=17) is
marginal; the other four are anecdotes (n=13, 11, 8, 3). Every headline in this project
rests on 7 ROIs chosen as the densest per domain — a rule `select_domain_images` documents
as optimistic. Add non-densest ROIs per domain and report per-domain medians with n.

---

## Step 3 — texture, contingent on Step 1

Only if Step 1 leaves a pipeline worth improving. Two corrections to the draft:

* **It targets 23% of the FP mass, not the majority.** Perfect look-alike suppression moves
  FPs 197 -> 151, while z=2.5 -> 2.0 already moved 301's recall@1000 by 0.185. Texture is
  worth less than choosing z properly. Rank it accordingly.
* **Run it on the `nucleus_blobs` candidate set specifically** — not merely "a set that
  carries look-alikes." It is the only one measured that captures ~100% of *both* classes,
  which is what makes the AUC unbiased and a DeLong test meaningful. The grid at step 20
  misses an unmeasured look-alike fraction.
* **But it is better motivated than the first draft of this plan allowed.** Unbiased on the
  blob set, chromatin density reaches AUC 0.692 (301) / 0.767 (246) against look-alikes, and
  the blob component-mean reaches 0.685 / 0.742 — statistically indistinguishable. Both
  available statistics are stuck in the same place, with 105-123 look-alikes to measure
  against. That is a real, well-powered gap for a different feature family to close.

GLCM offsets and window in **µm, converted per ROI** (measured mpp is 0.2269-0.2533, an 11%
spread). Sweep 1-8 px equivalent rather than asserting 1-3: a condensed chromosome is
~0.5-1.5 µm, and GLCM contrast peaks near half the period. Restrict GLCM to the darkest-decile
mask — over the full 51 px window it measures the neighbourhood, not the object, which is
exactly what `chromatin_density` avoids via its implicit segmentation.

State a **paired DeLong test** and a minimum detectable effect: at SE 0.033-0.038 the powered
cells need a ~0.09-0.11 AUC gain to detect at 80% power. Two confounds to state, not fix:
the AUC uses only *detected* look-alikes (32 of 109 on 301), a selection correlated with the
feature under test; and category 2 is pathologist-curated and non-exhaustive, so some of the
151 "unannotated" top-K FPs are unannotated look-alikes.

---

## Step 4 — contingent

* **`design_choices.md` §6 augmentation revisit** *before* any multi-seed template bank.
  `n_angles=1` was measured **worst of three** variants (discrimination -0.034 against +0.068
  for `rot90_4angles_2flips`) on a single seed and flagged "revisit"; never done. Now that
  the search is only a candidate generator, augmentation moves candidate-set recall directly
  — the quantity Step 1 shows is losing. ~4 s/run.
* **`peak_min_distance=7` vs `nms_radius`.** It sets maximum candidate density, which is what
  binds in the low-z regime; the grid saturates at ~17k after NMS and so does the pipeline at
  z<=1.0. That is not a coincidence and should be understood before the frontier is read.
* **Multi-seed template bank** fused by mean/median, not element-wise max.
* **Look-alike margin ranking** on the chromatin+texture vector, not on NCC.

---

## Cost, measured not guessed

Load + `rgb2hed` ~3 s; `matchTemplate` over 39 MP ~4 s; `extract_peaks` at a deep floor 0.3 s
(27-38k peaks, `max_peaks=250000` never binds); NMS 0.1-0.2 s; `chromatin.score_detections`
30-33 µs/detection; `nucleus_blobs` 7-11 s; grid 2-3 s. **The whole 7-ROI x 5-seed x 5-z
frontier is ~15 minutes.** Compute is not a constraint on this project and the plan should
stop rationing it — the draft's too-shallow extraction floor was justified on compute grounds
that do not exist.
