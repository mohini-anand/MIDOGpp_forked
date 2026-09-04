# Audit of the 2026-08-31 work (chromatin re-rank + premise test)

Date: 2026-09-01
Audits: commit `7c3af93` and the uncommitted premise-test work
(`midog_utils/compare.py`, `midog_utils/invariants.py`, `premise_test.py`,
`results/premise_*.csv`, `Research Logs/2026-08-31-{next-steps-plan,premise-test-results}.md`)
Method: every headline number re-derived from `results/premise_test.csv`,
`results/premise_seed_provenance.csv`, `results/premise_verification.csv`,
`results/od_seed_sweep.csv`, `results/od_workload_ab.csv` and
`results/morph_diag_bhattacharyya.csv`. Three claims re-measured from the raw images.

---

## Part 1 — what yesterday's work actually did, in plain words

Two pieces of work happened on 2026-08-31: a committed change to how detections are ranked,
and an uncommitted experiment that questioned the whole approach.

### The setup

A pathologist clicks **one** mitotic figure in a tissue image. The tool cuts a small square
around that click and slides it over the whole image looking for places that look similar
(template matching). Wherever the match is good enough, it plants a candidate. Overlapping
candidates are thinned out (non-maximum suppression). The pathologist then reads down the
resulting list, and the question that matters is: **how many candidates must they look at
before they have found half the real mitotic figures?** Fewer is better.

The problem has been false positives — the list is long and mostly wrong.

### Piece 1 (committed, `7c3af93`): rank by darkness, not by similarity

The insight is a mathematical one. The similarity score being used, `TM_CCOEFF_NORMED`,
subtracts the mean and divides by the spread of *both* the template and the patch it is
compared against. That makes it deliberately indifferent to brightness and contrast — a
pale washed-out copy of the seed scores essentially the same as a dark dense one.

But a separate measurement taken back in August (`morph_diag_bhattacharyya.csv`) had already
established that **mean intensity — how dark the nucleus is — is the single strongest
feature separating mitotic figures from ordinary nuclei**, in all seven tumour types.

So the search was throwing away the best available signal, by design. That explains a run of
six earlier experiments (score threshold, multi-scale, augmentation, colour channel, box
tightening, NMS ordering) that all came back neutral: each was just a different place to cut
the *same* ranking, and none of them added new information.

Two changes followed:

* **A new ranking statistic: "chromatin density."** For each candidate, take the darkest 10%
  of pixels in a 51-pixel window around it and average them, read off an *unclipped*
  hematoxylin channel. Keep exactly the same candidates the search found; only re-sort them.
* **A score floor in relative rather than absolute units.** The old cut-off `score >= 0.5`
  meant wildly different things on different images — 1.7 standard deviations on one image,
  4.5 on another, which is why one image returned 3-22 candidates against 115 real mitotic
  figures. Replacing it with "2.5 standard deviations above this image's own noise" made the
  search depth comparable across images.

**Result:** the number of candidates a reader must work through to reach 50% sensitivity fell
by a median factor of **19.4x**, winning 29 of the 29 runs where both rankings got there
(35 runs = 7 images x 5 seeds; 6 runs never reached 50% under either and were excluded).

An adversarial audit of that work caught a genuine defect on the way: the original
hematoxylin channel was clipped at 255, so dozens of the darkest candidates all scored
*exactly* 255 and their order was decided by an unstable sort. Fixing it improved the
result. I re-derived all of these numbers and they are correct (Part 2).

### Piece 2 (uncommitted): the premise test — is the click worth anything?

The second piece asks the question nobody had asked: **the machinery is built around one
human click, but does the click actually buy anything?**

Two comparisons, run over 7 images x 5 seeds x 5 threshold levels (~10 minutes of compute):

**(a) Does the click carry information?** Run the identical pipeline seeded from three
different things: a real mitotic figure, a look-alike the pathologist examined and
*rejected*, and a random nucleus. Answer: **yes, but only a little.** Seeding from a real
mitotic figure beats a rejected look-alike by about +0.05 recall and a random nucleus by
+0.01 to +0.07. Consistent in direction, but small.

**(b) Does the search earn its place?** Compare the whole pipeline against comparators that
never see any annotation at all: a nucleus detector (tiled Otsu → connected components), and
a **plain regular grid of points** ranked by the same chromatin-density statistic. Answer:
**the pipeline does not clearly win.** A bare lattice with no seed, no template and no
correlation is competitive with it. That is a much bigger finding than (a).

Two other things came out of it:

* **The committed threshold (z = 2.5) is a bad operating point.** On the largest image it
  fails to reach 50% sensitivity on 2 of 5 seeds, and where it does it swings between 247 and
  734 candidates. At z = 1.0 it reaches the target on 5 of 5 and lands within 172-203 every
  time. The default should move to z = 1.0.
* **Which cell the pathologist happens to click changes the reading burden by ~2x.** The
  seedless comparators are deterministic to within 1%. That is a product-level problem the
  median hides.

The experiment also caught two errors in its own brief: "distinct seeds must give distinct
annotations" is not true when the draw is with replacement, and the tissue mask excludes one
annotation, not zero (a category-2 look-alike 13 px from the image edge).

The engineering here is genuinely good — a reusable comparison harness (`compare.py`), four
structural invariants (`invariants.py`) run 1,153 times, and a report generated from the CSVs
so no number can drift. The mistakes below are almost all in how the results are **compared
and described**, not in how they were measured.

---

## Part 2 — what survived the audit unchanged

Re-derived independently and correct to the last digit.

**From the premise test:**

* The seed-provenance table (246: 0.809 / 0.757 / 0.748; 301: 0.710 / 0.659 / 0.668) and
  every paired win count and delta (4/5 +0.052, 3/5 +0.051, 5/5 +0.070, 3/5 +0.014).
* The 301.tiff `read_50` per-seed table (z=2.5: 247, —, 734, 313, —, **2 of 5 unreached**;
  z=2.0: 156/203/176/172/349; z=1.0: 174/195/172/188/203) and the seedless arms'
  determinism (blob_native 167-168, grid_20 181-184).
* "z = 1.0 is best on all 12 decision-grade cells" — z=1.0 is the argmax in 12/12.
* `coverage_frac > 0.5` in **413 of 561** arm-cells. Row counts 2244 / 120.
* **1,153** invariant checks (no_cap 591 + nms_radius 520 + min_separation 28 +
  tissue_mask 7 + distinct_seeds 7).
* `tissue_mask` excludes **1 of 726** annotations — 506.tiff ann 24249, category 2 — and
  **0 of 391** mitotic. Both totals check out.
* Shortcut exactness: 7/7 ROIs, `n_direct == n_filtered` and coordinate-identical.
* The grid-suppression confound is real and was disclosed: `grid_30` on 301.tiff is the
  **raw** lattice (34,615 points, min separation 30.0 px against a 29.61 px radius) while it
  is suppressed on the other six ROIs.
* Diagnostics confirm the saturation fix worked: largest tie block is 2-3 everywhere,
  `nan_rate` < 2.8%, `floor_limited` False in all 561 cells.

**From commit `7c3af93`:**

* The whole Fix-2 workload table reproduces from `od_workload_ab.csv`, and the median of the
  seven ratios is **19.36 → 19.4x**.
* Over the 35-run sweep: median **19.36x**, IQR **10.72-39.56**, **29/29** wins,
  **6** runs unreached under both rankings and excluded, **22.62x / 27/27** excluding the two
  floor-limited runs, median full-list recall **0.455 → 0.727**, list-length spread
  **5,161x → 6.4x**, `recall@K` ≥ in **35/35** and strictly better in **29**. Every one
  matches. **The 19.4x claim does *not* have the survivorship defect described in M2 below —
  the exclusions are counted and disclosed.**
* Top-K composition: TP 187 / look-alike 46 / unannotated 151 under chromatin ranking, and
  87 / 19 / 278 under score. The full-list ratio 28,284 : 166 = **170:1**. All exact.
* `mean_intensity` is the top Bhattacharyya feature in all 7 domains (0.925-3.703); the
  runner-up is `tightened_size` in 4, `solidity` in 2, `area` in 1, max 1.051. The log's
  correction of its own earlier "next-best solidity" claim is itself correct.
* The z-of-0.5 table (301.tiff = 1.72) reproduces from `od_seed_sweep.csv`.

**Cross-validation between the two experiments.** `od_seed_sweep.csv` and
`premise_test.csv` were written by separately authored scripts on separate days. On 301.tiff
at z=2.5 they agree exactly on all five seeds — same seed annotation ids (14976, 14754,
14817, 14787, 14816), same `read_50` (247, unreached, 734, 313, unreached). That is strong
evidence both pipelines are doing what they claim.

---

## Part 3 — mistakes

### M1 — "the pipeline loses to seedless" is measured against a per-cell oracle (material)

`2026-08-31-premise-test-results.md` reports 2/15, 3/15, 0/15, 0/15 wins and median deltas
−0.017 / −0.026 / −0.052 / −0.041, and frames the protocol as generous to the pipeline:
*"The pipeline is allowed to pick its best z post hoc per ROI and seed — a handicap in its
own favour."*

Those numbers are arithmetically correct, but the seedless side is **also** maximised post
hoc, over **six** arms, per (ROI, seed, budget) cell — a wider selection than the pipeline's
five z levels. The winning seedless arm changes by ROI and by budget (301: `blob_native` in
20/20 cells; 201: `grid_10` in 19/20; 246: `grid_30`, then `blob_native`, then `grid_20` as
the budget rises). No deployable system can make that choice.

The full matrix — pipeline **pinned at z = 1.0**, the level the log itself now recommends, so
neither side gets a post-hoc choice. Cell is (wins out of 15, median paired delta):

| seedless arm | @500 | @1000 | @2000 | @5000 |
|---|---|---|---|---|
| `grid_10` | **10, +0.046** | **10, +0.023** | 4, 0.000 | 1, −0.026 |
| `grid_20` | 2, −0.009 | 3, −0.017 | 5, −0.023 | 0, −0.035 |
| `grid_30` | **8, +0.009** | **8, +0.014** | 7, 0.000 | 4, −0.009 |
| `grid_40` | **11, +0.017** | 3, 0.000 | 4, −0.014 | 4, −0.009 |
| `blob_native` | 7, 0.000 | 5, −0.009 | 5, −0.026 | 7, 0.000 |
| `blob_od` | **12, +0.078** | **8, +0.043** | **9, +0.017** | **10, +0.035** |
| *oracle, max of 6* | *0, −0.026* | *3, −0.026* | *0, −0.052* | *0, −0.041* |

The pipeline beats four of six seedless arms at budgets ≤ 1000 and loses decisively only to
`grid_20`. It loses to the oracle at every budget — but the oracle is worth 0.019-0.053 of
median recall over the best single arm (@1000: 0.931 oracle vs 0.878 for `grid_30`), which is
**larger than the deficit being reported**.

Per ROI against `grid_30`, pinned at z = 1.0: 301 wins at 500/1000/2000, 246 loses at all
four, 201 wins at all four. So the defensible statement is *"against one pre-committed
seedless arm the pipeline wins on 2 of 3 decision-grade ROIs at budgets ≤ 2000"* — not
"the pipeline loses."

*Fix.* Publish the whole matrix above rather than one comparator, keep the oracle row
clearly labelled as an upper bound, and delete the "handicap in its own favour" sentence — it
describes the opposite of what the procedure does.

### M2 — the reading-depth win on 201.tiff is a 3-of-5-seed survivorship median (material)

The log's table:

| ROI | pipeline (best z) | best seedless |
|---|---|---|
| 201.tiff | **28** (z=2.5) | 40 (`blob_od`) |

At z=2.5 on 201.tiff, **2 of 5 seeds never reach 50% sensitivity**. `28` is the median of
{28, 28, 200}. This is the exact failure mode the same document condemns two sections later
on 301.tiff, and which the plan's Step 0 rule 5 forbids ("runs where no arm reaches the
target are reported separately, never counted as wins"). The z per ROI is also picked by
whichever minimises `read_50` (2.5 on 201, 2.0 on 246 and 301) rather than by one rule.

The auto-generated `results/premise_report.md` §5 gets this right — one rule (best z by median
recall@1000 → z=1 everywhere), an explicit `unreached` column, and a warning that a low
`read_50` beside a high `unreached` is *worse* than the row above it. The hand-written log
silently disagrees with the report generated from the same CSV.

Under the report's rule the verdict is the same but the margins are not:

| ROI | log | correct (z=1, 0 unreached) | seedless | paired wins |
|---|---:|---:|---:|---|
| 201.tiff | 28 | **38** | 40 (`blob_od`) | 4/5 |
| 246.tiff | 106 | **111** | 138 (`grid_20`) | 5/5 |
| 301.tiff | 176 | **188** | 168 (`blob_native`) | 0/5 |

201 goes from a 30% advantage to a 5% one; 301 from 5% behind to 12% behind.

### M3 — the blob arms are the only arms never suppressed (measured)

`grid_*` arms get `nms_radius=radius` and a `check_min_separation`; the pipeline is NMS'd
inside `match_pool`. `blob_native` / `blob_od` get neither, and `Arm.nms_radius` is left
`None`, so `check_nms_radius` silently skips them.

Measured on the raw images:

| ROI | n blobs | min separation | share within match radius of another blob | after NMS |
|---|---:|---:|---:|---:|
| 301.tiff | 20,370 | 3.50 px | **58.6%** | 14,518 (71.3%) |
| 246.tiff | 21,059 | 1.44 px | **70.5%** | 13,413 (63.7%) |

Duplicates waste budget, so this **handicaps the blob arms** — the ones that win. The headline
is therefore conservative, but "matched budget" is not "matched suppression". Re-measured
with NMS at the match radius (seed 0):

| arm | 301 @500 | 301 @1000 | 301 read_50 | 246 @500 | 246 @1000 | 246 read_50 |
|---|---:|---:|---:|---:|---:|---:|
| `blob_native` raw | 0.802 | 0.931 | 168 | 0.765 | 0.861 | 193 |
| `blob_native` NMS'd | **0.811** | **0.935** | **164** | **0.791** | **0.870** | **183** |
| `blob_od` raw | 0.631 | 0.760 | 193 | 0.678 | 0.774 | 211 |
| `blob_od` NMS'd | **0.682** | **0.770** | **164** | **0.730** | **0.817** | **148** |

Note this changes M2's table: NMS'd `blob_od` reaches 164 on 301, better than `blob_native`'s
168, so the 301 seedless cell becomes 164. **Fix M3 before republishing the M2 table.**

### M4 — `check_no_cap` is vacuous on the arms it matters for

`premise_test.py` bypasses `find_and_suppress` and calls `tm.extract_peaks` /
`nms_by_distance` directly, so:

* `cfg.max_detections` (1e9) is **never applied** on this code path.
* `cfg.max_peaks` (250,000) **is** applied, inside `extract_peaks` — but to the *pre-NMS*
  peak count `diag["n_peaks"]`, while `check_no_cap` is handed `len(sub)`, the z-filtered
  post-NMS length. A filtered length can never equal the extraction cap, so the check passes
  by construction.

`n_peaks` is printed to stdout and written to no CSV, so it cannot be audited after the fact
either. This is invariant #1 — written precisely to catch `nucleus_blobs`' old silent
`max_detections=20000` — defeated in the code meant to enforce it. No numbers are wrong:
250,000 is above the ~173k theoretical maximum, so nothing truncated.

### M5 — the provenance control's arm (c) carries an exclusion the others don't

`gt_prov` drops seed (a)'s and seed (b)'s `ann_id`, so those two arms' `self_hit_radius`
exclusion zones sit on annotations that are gone from the ground truth. Arm (c)'s seed is a
`nucleus_blobs` centroid — not an annotation — so **nothing is dropped for it**, while
`match_pool` still strips every detection within 5 px of it.

Measured from the stored `seed_cx`/`seed_cy`: on **301.tiff seed 1 the drawn blob is 2.8 px
from a mitotic figure**, inside the exclusion zone. One of ten provenance draws.

Direction: biases (c) down, i.e. **inflates** "the click carries information". Materiality
here is nil — one mitosis is 0.0046 recall on 301, and that seed is already a loss for (a)
(a = 0.682, c = 0.747), so the reported +0.0138 median and 3/5 win count do not move. Latent
defect, immaterial in this run, will bite in a larger one.

### M6 — `coverage_frac` is a full-list number used to defend budgeted ones

`compare.evaluate_arms` computes coverage on the **whole** ranked list. Both documents then
make a claim about coverage *at budget*: the plan says "coverage is matched across arms
(0.09-0.11 on 301 at budget 1258)", the log says "the budgeted columns — which carry every
conclusion above — are unaffected". Neither is supported by anything in this run: all 210
seedless arm-cells have full-list coverage 0.75-1.00, no top-K coverage is stored, and the
detections are not saved, so it cannot be recovered post hoc.

### M7 — four retracted claims still in committed code

The audit corrected the prose and left the source asserting what it retracted:

| location | text | status |
|---|---|---|
| `midog_utils/chromatin.py:1` | "the signal `TM_CCOEFF_NORMED` is mathematically blind to" | retracted, log audit item 8 ("much weaker ranker", AUC 0.61-0.81) |
| `midog_utils/chromatin.py:11` | "against next-best `solidity` at 0.15-1.05" | retracted, item 7 — runner-up is `tightened_size` in 4 of 7 domains |
| `midog_utils/chromatin.py:111` | "The search still earns its place as the *candidate generator*" | now contradicted by the premise test itself |
| `od_experiment.py:11` | "provably blind to it" | same as row 1 |

### M8 — "a blob detector that never sees an annotation does better" is true on one ROI

Median recall@1000, pipeline at z=1.0 vs `blob_native`:

| ROI | pipeline z=1 | `blob_native` |
|---|---:|---:|
| 201.tiff | **1.000** | 0.588 |
| 246.tiff | 0.852 | 0.861 |
| 301.tiff | 0.747 | **0.931** |

The pipeline is not merely competitive on 201.tiff, it is 0.41 ahead. The median-over-ROIs
framing hides a sign flip between the two ROIs that carry the most weight.

### M9 — the plan's grid numbers do not reconcile with the run

Every `pipeline` / `blob_*` number in `2026-08-31-next-steps-plan.md` reproduces exactly in
`premise_test.csv`. Every **grid** number differs:

| cell | plan | run |
|---|---:|---:|
| 301 grid_20 @500/1000/2000/5000 | 0.700 / 0.816 / 0.899 / 0.972 | 0.696 / 0.802 / 0.876 / 0.968 |
| 246 grid_20 @500/1000/2000/5000 | 0.783 / 0.887 / 0.930 / 0.983 | 0.765 / 0.861 / 0.922 / 0.991 |
| 246 grid_20 read-50 | 144 | 138 |

Almost certainly the `step // 2` offset `baselines.grid_lattice` added after the plancheck
script. The plan is uncorrected, and still cites `scratchpad/plancheck/fiveseed.csv` as "the
load-bearing one" — that directory lives under `/private/tmp/claude-501/...`, is ephemeral,
and is not in the repo. Same reproducibility gap the previous audit's item 6 flagged.

### M10 — none of the premise-test work is committed

Untracked: `midog_utils/compare.py`, `midog_utils/invariants.py`, `premise_test.py`,
`results/premise_{test,seed_provenance,verification}.csv`, `results/premise_report.md`, both
new research logs. Modified: `midog_utils/baselines.py`, the chromatin log.

### Minor

* `2026-08-31-chromatin-density-rerank.md` and the commit message say "**30** distinct seed
  annotations"; `od_seed_sweep.csv` has **31** (5+4+5+5+4+4+4). Off by one.
* `seed_pool_size` leaks the last seed's value out of the loop into `check_distinct_seeds`.
  Harmless — `pick_seed`'s three filters are seed-independent — but it reads as a bug.
* `Arm.nms_radius=radius` is declared on `grid_30` / `grid_40` where suppression removes
  nothing (lattice spacing already exceeds the radius). The field then means "the radius I
  would have used".
* `premise_test.py`'s module docstring says `verify_shortcut` "asserts it directly against a
  re-extraction on every ROI" — true, but only at z=2.0 and seed 0. The newly recommended
  z=1.0, the level nearest the extraction floor, is the unverified one.

---

## Part 4 — next steps

### P0 — cheap code fixes first, then one ~10 min re-run

The reporting fixes depend on these, so they come first (M3 changes M2's table).

1. NMS the blob arms at the match radius, ordered by each arm's own key, and set
   `nms_radius` so the check fires (M3).
2. Assert `diag["n_peaks"] != cfg.max_peaks` in `match_pool`; add `n_peaks` to `ctx` (M4).
3. Build one common exclusion set across the three provenance arms — drop every annotation
   within `self_hit_radius` of *any* of the three seeds (M5).
4. Emit `coverage_at_budget` per budget row; save detections for the three decision-grade
   ROIs (M6).
5. Verify the shortcut at z=1.0 as well as z=2.0.
6. Re-run `premise_test.py` (~10 min) and regenerate the report.

### P1 — make the write-up honest (~1 hour, no compute)

7. Publish the full 6-arm matrix pinned at one z; drop the "handicap in its own favour"
   framing (M1).
8. Rewrite the reading-depth table under one z rule with an `unreached` column (M2).
9. Split the "blob detector does better" sentence per ROI (M8).
10. Fix the four retracted docstrings (M7) and the 30→31 count.
11. Reconcile or retire the plan's grid table; copy `plancheck/` into the repo or delete the
    citation (M9).
12. Commit everything (M10).

### P2 — the science, in order

13. **Settle the comparator protocol before anything else.** The load-bearing question of the
    whole project — does the click earn its place — currently flips on whether the seedless
    comparator is an oracle or a fixed arm, and flips again per ROI. Pre-commit one baseline,
    state it in `design_choices.md`, and re-derive. Nothing downstream is worth doing until
    this is decided.
14. **One z level below 1.0.** The frontier is still improving at the edge of the swept range
    on all 12 decision-grade cells; the plan already concedes this is owed.
15. **Widen the ROI base** (plan Step 2). Four of seven ROIs have n < 15 mitotic figures, and
    the two ROIs that carry weight disagree in *sign* on the headline (M8). Seven
    densest-per-domain ROIs cannot settle a question this close.
16. **Texture last.** It targets 23% of the top-K false-positive mass, and the log already
    ranks it below choosing z properly. Not worth starting until 13-15 are done.
