# Pre-registration F5: what does shrinking the NMS radius do, with border padding held fixed?

Date: 2026-09-04. **Written before `f5_nms_radius_ablation.py` is run.** Status: PLAN — revision 3,
after two rounds of adversarial review. Not yet executed.

F1–F4 are taken (`Research Logs/2026-09-03-tm-axis-sweep-edit-plan.md` §Follow-ups;
`Research Logs/2026-09-04-f4-preregistration.md`). This is the next free number.

> **Revision 2 (same day).** Ten changes forced by review, each recorded in §11 with what it
> would have cost. The largest: the two shrink levels were originally justified as a *mechanism
> decomposition* (5.9 µm = "stop merging distinct annotations", 5.0 µm = "additionally stop the
> same-object secondary peak"). Measurement kills that: **zero** annotation pairs on the 14-ROI
> set sit between the two radii, so they are identical with respect to GT un-merging. The ladder
> is now registered as a one-mechanism **dose–response**, and the prediction that depended on the
> decomposition is withdrawn rather than reinterpreted after the fact.
>
> **Revision 3 (same day).** Four further changes, §12. One is a regression revision 2 introduced
> while fixing a revision-1 defect: the primary estimand became `worst-of-5(shrink) −
> worst-of-5(control)`, which takes each arm's minimum **independently** and so throws away the
> pairing §3b exists to create. One is a factual error carried since revision 1:
> `results/tm_ccoeff_threshold_axis_sweep_v2.csv` is the **r5.0** arm, so the
> `full_list_recall = 1.0` at z ≤ 1.0 that P1 leaned on is the *shrunk* arm's number, never the
> control's. The control has never been run at all, which is the point of F5 and makes P1 more
> live than revision 2 conceded.

---

## 1. The question, in the terms it was asked

`tm_threshold_axis_sweep_v2.ipynb` changed **two** things relative to
`tm_threshold_axis_sweep.ipynb` (v1) and reported their combined effect:

1. **border padding**, so `cv2.matchTemplate` has a full window everywhere and the whole ROI
   becomes reachable (v1's `valid` mask was `False` in a ~25 px margin), and
2. **a decoupled, smaller NMS radius** — 5.0 µm (19.7–22.1 px) instead of the evaluation match
   radius of 7.5 µm (29.6–33.1 px).

v2 says so itself, in its own cost table:

> Note "the cost of *both* fixes": v1 differs from this run by padding *and* radius, so this
> delta cannot be attributed to the radius alone.

**F5 attributes it.** Padding is held **ON in every arm**; the NMS radius is the only thing that
moves. The primary question, in the user's words: does shrinking the radius change *how many true
positives we surface in the fewest candidates*, and does the candidate list carry fewer false
positives?

---

## 2. What is already known, so this run is not re-answering it

| established | where | consequence for this design |
|---|---|---|
| Padding fixes a real leak: on `301.tiff` ann 14779 sits 9 px from the ROI edge and produces no candidate at any threshold without it. `valid.all()` after padding is provable from `PAD = (t-1)//2`, not empirical | v2 cells 6–8 | Padding is **not** under test. It is held on, and `valid.all()` is re-asserted per cell as a running check, not as a result. |
| The NMS mechanism is real and traced: on `301.tiff`, 5 mitoses had their on-centre peak (4.5–10 px from the click) suppressed by an off-centre peak 23–29 px away **on the same nucleus** | v2 cells 9–11 | This is the **only** live mechanism (see §2a). It is measured on 1 ROI and 1 seed, which is why §6.5 re-measures the gap band on all 14 ROIs. |
| The earlier radius sweep, read **paired within seed** on `301.tiff`: `read_50` rises monotonically as the radius shrinks on **all 5 seeds** (seed 0: 198 → 250 → 265 → 277 → 293 at 29.6 → 24 → 20 → 16 → 12 px), while every seed that had a miss reaches 0 misses at **every** shrink level | `results/tm_nms_radius_sweep.csv` | **The expected direction is a trade, not a win**: the head cost is dose-dependent and universal across seeds; the tail gain is real and saturates immediately. So there is a crossover K, and a design that samples only K = 250 can land on one side of it and miss it. |
| That sweep covered **2 ROIs, both ones where TM was losing** — a deliberately biased sample — and its own log says "do not read this as a recommendation to set the radius to 24" | `Research Logs/2026-09-01-tm-variant-sweep.md` §5f | F5 exists to replace that biased 2-ROI sample with the balanced 14-ROI one. |
| `full_list_recall` is coverage-saturated at these list lengths (88.1–99.8 % of *arbitrary* ROI locations already sit within a match radius of some detection) | v2 cell 18; `evaluate.py` module docstring | A smaller radius keeps more peaks → higher coverage → mechanically higher full-list recall. Full recall **cannot** be the headline; and the same guard must be extended to `read_99`/`read_100` (§6.4). |
| The reporting metric is **recall@K**, worst-ROI and worst-click, per tumour domain | `DECISIONS.md` D4 | Fixes the metric family before any data is seen. |
| Results CSVs emit one row per budget, so any count/test not deduplicated multiplies n | `DECISIONS.md` D4, "one data trap" | Dedup key fixed in §8; note F5's arms are named differently from v2's, which changes the key. |
| `read_50` is the depth at which methods matter least — it reproduces with no click at all | `DECISIONS.md` D4 §5 | `read_50` is a diagnostic here, never an arbiter. The prior sweep's headline was `read_50`; F5's is not. |
| `read_95` has a seed-to-seed SD of 260–2744 candidates | `results/tm_ccoeff_headtohead_seed_variance.csv` | Single-seed depth deltas are **smaller than seed noise**. 5 seeds is not optional. |

### 2a. The annotation-spacing geometry, measured on all 14 ROIs — and what it rules out

v2 cell 0 and the F4 pre-registration state the minimum true annotation spacing as 26.2 px
documented / 36.8 px measured; the 36.8 figure is over the **7-ROI** set. Recomputed over the
14-ROI `images/extra_valid` set, in microns (the only unit in which a µm-scaled radius can be
compared across scanners whose mpp runs 0.2263–0.2533):

| quantity | value | where |
|---|---|---|
| closest annotation→annotation pair | **5.955 µm** (26.25 px) | `403.tiff` |
| closest mitotic→mitotic pair | **6.594 µm** (26.57 px) | `245.tiff` |
| evaluation match radius | 7.5 µm (29.61–33.14 px) | all 14 |

Pair counts inside each candidate radius, over all 14 ROIs:

| radius | mitotic↔mitotic pairs | any↔any pairs |
|---|---:|---:|
| 7.5 µm (control) | **1** (245.tiff) | **7** (1 on 245, 6 on 403) |
| 5.9 µm | 0 | 0 |
| 5.0 µm | 0 | 0 |

**Three things this settles, two of them against the first draft of this document.**

1. **The GT-un-merging story is negligible, not pervasive.** 12 of 14 ROIs have *zero* annotation
   pairs inside even the default radius, and across ~1,325 mitotic annotations exactly **one**
   mitotic pair is close enough for the default suppression to merge. Recall is computed over
   category 1, so the entire un-merging mechanism has a ceiling of **one annotation**. Revision 1
   of this document said "the default NMS radius exceeds the closest real annotation pair on
   every single ROI" — true as a cross-ROI pixel statement, and misleading as a mechanism claim.
2. **F4's claim survives the re-measurement intact.** F4 states that shrinking cannot recover a
   merged pair of ground-truth mitoses; it says this about the **5.0 µm** radius, and
   6.594 µm > 5.0 µm holds on the 14-ROI set. Revision 1 presented this as a correction to F4. It
   is not one, and the overstatement is withdrawn.
3. **No radius level can separate the two mechanisms.** Zero pairs lie between 5.0 and 5.9 µm, so
   the two shrink levels are *identical* with respect to GT un-merging; and the one mitotic pair
   that exists (26.57 px) falls **inside** the 23–29 px same-object suppression band, so the two
   mechanisms are not separable by any choice of radius on this data. The ladder is therefore a
   dose–response on one mechanism, and is registered as such in §3a.

---

## 3. Design

### 3a. The one factor, three levels — a dose–response, not a decomposition

Every radius is physically scaled in µm — the same way `evaluate.radius_px` and v2's NMS radius
are already defined.

| arm tag | radius | px over the 14 ROIs | *predicted* position vs. the 23–29 px same-object gap band — **to be confirmed by §6.5, not assumed** | why this value |
|---|---|---|---|---|
| `r7.5` **(control — "variant 1")** | 7.5 µm = `ev.MIDOG_RADIUS_UM` | 29.61 – 33.14 | **above** the band — suppresses the whole of it | The repo default, and the only value satisfying `invariants.check_nms_radius`. This is v1's suppression with v2's padding: **the arm that has never actually been run.** |
| `r5.9` (mid dose — "variant 2a") | 5.9 µm | 23.29 – 26.07 | **inside** the band — suppresses part of it | The largest µm radius strictly below the closest real annotation pair anywhere in the set (5.955 µm on `403.tiff`), rounded down. So it is the mildest shrink that cannot merge two annotations — and, given §2a, that property buys nothing on its own, which is exactly why it is registered as a dose rather than a mechanism arm. |
| `r5.0` (full dose — "variant 2b") | 5.0 µm | 19.74 – 22.09 | **below** the band — suppresses none of it | v2's value, so every F5 number is directly comparable to v2's. Well above `peak_min_distance = 7`, so genuinely redundant peaks on one object still merge. |

**The band column is a prediction, not a property.** The 23–29 px band comes from **one ROI and
one seed** (v2 cells 9–11). Writing "above / inside / below the band" as though it were a fact
about the ladder would let §6.5 only ever confirm the placement that was used to justify it. The
column is therefore registered as *predicted*, and §6.5 re-measures the band on all 14 ROIs; if it
runs 12–18 px elsewhere, `r5.0` is above it too and the ladder is mis-centred — which is a
reportable finding, not a failure to be quietly absorbed.

**What the ladder can and cannot show.** It is a monotone dose on one mechanism — the same-object
secondary peak. If the effect is mechanical, `r5.9` should sit between `r7.5` and `r5.0` on every
metric; a non-monotone result is evidence the mechanism is not what v2 traced. **Three levels is
the minimum that carries that falsifier at all** — with two, there is no monotonicity to test, and
§4's replacement for the withdrawn P4 would have no registered way to fail. It **cannot** attribute
anything to GT un-merging, and no design on this dataset can.

Held constant in every arm: border padding ON (`PAD = max((t-1)//2)`, `BORDER_REPLICATE`);
`OD_PAD = 25` chromatin padding ON; `hematoxylin_od`; `TM_CCOEFF`; single 51 px template;
`peak_min_distance = 7`; `self_hit_radius = 5.0` **applied after NMS** (v2's order — see §7.10);
and **the scoring match radius pinned at 7.5 µm in all three arms** — suppression varies, scoring
never does.

### 3b. The pairing — what makes this "very controlled"

Per `(ROI, seed)`: compute the padded fused response **once**, extract the deep peak pool
**once**, then apply all three NMS radii to the *identical* `(centers, scores)` arrays. Nothing
upstream of `nms_by_distance` is radius-dependent — verified in code, not docstrings: the seed
draw (`tm_variant_sweep.py:252-265`), `robust_stats` (`template_match.py:66-83`), the deep floor,
and `extract_peaks` (`template_match.py:249-276`) all take no radius, and `nms_by_distance`
(`nms.py:32-57`) never mutates its inputs.

So the three arms differ by exactly one argument, and every reported comparison is a
within-`(ROI, seed, z, axis)` paired difference with the response map, template, seed draw, floor,
peak set and ground truth bit-identical between arms.

**Implementation trap, named so it cannot happen:** the pool arrays must never be rebound
(`centers, scores = centers[keep], scores[keep]`, as v2's single-arm loop does) or arm 2 would
suppress arm 1's survivors. Each arm reads the shared pool and writes a new local name.

**`chromatin_od` is hoisted out of the arm loop for the same reason the peak pool is.** `od`
depends only on `(cx, cy)` and the padded image, so it is identical for any point shared between
arms. It is scored **once on the deep pool per (ROI, seed)** and subset per arm, which turns
`r5.9` from a 3.5-minute addition into a net *saving* of ~126 s against scoring it three times.
The equivalence is asserted, not assumed — F4 §9 gate 7 already asserts exactly this
("`od`-once-per-template equals `od`-after-NMS"), and F5 borrows the gate with the optimisation.

This is also why F5 is **one run with three arms and not two notebooks**. Two separate runs would
re-do the matching (3× the compute) and, worse, would make the comparison unpaired for no gain.
The user's "two variants" framing is preserved in the *reporting*: variant 1 is the `r7.5` control,
variant 2 is the `r5.9`/`r5.0` pair.

### 3c. Sample

* **14 ROIs** — `images/extra_valid`, 2 per tumour type across all 7 domains. Fixes v2's own
  stated limitation ("one ROI per domain").
* **5 seeds per ROI** — `tm_variant_sweep.draw_seeds` (agreement pool → border filter at 36 px →
  drawn **without replacement**, one `default_rng([s, image_id])` stream per seed index). Every
  one of the 14 has a seed pool ≥ 7, so all 5 draws exist everywhere. Fixes v2's other stated
  limitation (`SEED_INDEX = 0` only).
* **70 cells**, clustered in **14 ROIs**. The ROI is the design cluster, not the cell — see §8.

### 3d. Nuisance axes, carried but not under test

* **`z ∈ (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)`, headline at z = 1.0** — the operating point set in
  `tp_fp_score_distribution.ipynb`, *not* v2's `CURRENT_Z = 0.5`, which v2 itself flags as
  retained only because its before/after snapshot was taken there. v2's z = −1.0, −0.5 and 0.0
  levels are dropped: in the v2 run `full_list_recall` is 1.0 on every ROI at every z ≤ 1.0 while
  `coverage_frac` is 0.96–0.998, so those cells are pure saturation and cannot discriminate.
  (That v2 run is the **r5.0** arm — see P1 — so it bounds what the *shrunk* arm does there, not
  the control.) 0.5 is kept for v2 comparability. **z = 2.5 and 3.0 are kept for a stronger reason
  than "the control fails there": they are the only region in the entire design where
  `coverage_frac` is not saturated** — measured on the v2 CSV it runs 0.880–0.988 at z = 1.0,
  0.621–0.774 at z = 2.0, but **0.434–0.625 at z = 2.5 and 0.262–0.469 at z = 3.0**. That is the
  only place §6.3's guard stops binding and full-list recall can carry evidence at all, which makes
  those two levels the sole testable region for P1 rather than a robustness afterthought.
  (`budget_delivered < budget` bites there only at K ≥ 5000 — 201.tiff falls to 3,202 detections at
  z = 3.0 — never at the K = 250 headline, and §6.1's guard catches it.)
* **`DEEP_FLOOR_Z = −1.5`, pinned as a constant and deliberately decoupled from `min(Z_LEVELS)`.**
  v2 defines it as `min(Z_LEVELS) − 0.5`; keeping that coupling while trimming the z grid would
  silently change the deep pool and break §7.9's reproduction gate. The decoupling is **free, and
  for the same reason the one-match-many-z shortcut works**: in `nms.py:52-57` the loop walks
  descending score and `suppressed[]` is written only by already-visited, strictly higher-scoring
  indices, so a peak at −1.4 can never suppress a peak at +0.6. Every peak below `min(Z_LEVELS)`
  is inert with respect to the kept set at every reported z. `min(Z_LEVELS) − 0.5` was a
  convenience, never a semantic coupling. `assert_floor_not_limiting` still holds trivially.
* **Ranking axis ∈ {`tm_score`, `chromatin_od`}**, both reported; `chromatin_od` is primary
  (it is the shipped ranker, commit `7c3af93`, and F4's primary axis).
* **Arm naming: `{axis}@{radius_tag}`** (e.g. `chromatin_od@r5.0`), so `arm` is unique per
  (axis, radius) and the D4 dedup key still works. `coverage_key = (file, seed, radius, z)` —
  coverage depends on the candidate *set*, so on radius and z, but not on the sort order.

---

## 4. Predictions, registered before the run

Stated so they can fail.

* **P1 (tail).** Full-list recall is non-decreasing across the ladder.
  **Revision 3 corrects the evidence this was resting on.** Revision 2 said full recall "is
  already 1.0 for the control at z ≤ 1.0" and cited v2's CSV. That CSV records
  `nms_radius_px` ∈ [19.74, 22.09] — it **is the r5.0 arm**. The 1.0 belongs to the *shrunk* arm.
  The control (padded, r7.5) has never been run, and v1 (unpadded, r7.5) was *below* 1.0. So P1 is
  a genuine prediction at every z, not a near-vacuous one — the error was attributing a measured
  value to an unmeasured arm.
  It is still true that a full-recall gain at z ≤ 2.0 is largely geometric: `coverage_frac` there
  is 0.62–0.99, so §6.3's guard governs how it may be read. **z ≥ 2.5 (coverage 0.26–0.63) is
  where P1 can carry evidence**, and it is the only such region.
* **P2 (head).** recall@K at small K (≤ 250) is **flat or worse** under the shrink, dose-dependently,
  because every extra surviving duplicate on an already-found object consumes a rank slot without
  adding a TP. `n_detections` rises monotonically across the ladder. *Confidence: high — the prior
  sweep's `read_50` rises monotonically with the dose on all 5 seeds of 301.tiff, and on all 5 of
  300.tiff (paired Δ at 29.6 → 20.0 px: +3, +20, +4, +6, +6).* The sign at
  K = 250 is the live question.
* **P3 (crossover).** There exists a K\* at which the shrunk arms overtake the control. If
  K\* > 5000 the shrink is irrelevant to the product, since no pathologist reads 5000 candidates.
  **The decision-relevant output of F5 is K\*, per domain**, not a single recall@250 number.
* **P4 — withdrawn.** Revision 1 registered a mechanism-decomposition inference ("if `r5.9` ≈
  `r7.5` and only `r5.0` moves, the mechanism is the same-object peak, not GT merging"). §2a shows
  the data cannot support that inference at any radius. Withdrawn before the run rather than
  reinterpreted after it. What replaces it is weaker and honest: **monotonicity across the dose**
  is evidence the traced mechanism is what is operating; non-monotonicity is evidence it is not.

**What would make F5 report "shrinking is not worth it"** — bound to §8's arbiter, not to a
descriptive: the mean-of-5 paired Δ(recall@250) at z = 1.0 on `chromatin_od` has a cluster
bootstrap CI containing 0 at both doses, K\* is absent below 5000, and `n_detections` is up ≥ 20 %.
That combination is a clean negative and will be reported as one. The worst-of-10 per-domain table
is reported beside it and can qualify the reading, but does not overturn it.

---

## 5. Deliverables

1. `f5_nms_radius_ablation.py` — the run. Writes:
   * `results/f5_nms_radius_ablation.csv` — the tidy long frame from `compare.evaluate_arms`,
     one row per `(file_name, seed_index, arm, z, budget)`.
   * `results/f5_nms_radius_ablation_verification.csv` — every invariant and gate record (§7).
   * `results/f5_nms_radius_ablation_ledger.csv` — the per-annotation win/loss ledger and the
     same-object gap measurement (§6.5).
   * `results/f5_nms_radius_ablation_spacing.csv` — §2a's table.
2. `f5_nms_radius_ablation.ipynb` — reads those CSVs; no matching. Analysis, plots, prose.

Split this way because the run is ~15 min of CPU and the analysis will be re-run many times —
the same reason `tm_variant_sweep.py` and `tm_variant_report.py` are already split.

---

## 6. Metrics, and which one arbitrates

### 6.1 Primary — recall@K (D4)
`recall_at_budget` over a **dense** budget grid, a strict superset of `compare.BUDGETS` so rows
stay comparable with every earlier CSV:

```
BUDGETS_F5 = (10, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500,
              750, 1000, 1500, 2000, 3000, 5000, 10000)
```

Dense because it costs nothing (pure indexing into `tp_cum`) and because §2 predicts a crossover
that a sparse grid would step over.

**Reporting frame (D4): recall@250 per domain, worst ROI and worst seed (worst-of-10), per arm,
on the `chromatin_od` axis — reported without a p-value.** This is a *descriptive* frame, and D4's
"per tumour domain, worst ROI, worst click" is a reporting rule, not an inference rule. The
quantity that arbitrates is declared in §8 and is a different statistic; revision 2 conflated the
two by asserting they were the same, and they are not.

Two guards, both mandatory before any cell is read:

* **`budget_delivered == budget`.** The control has the *shortest* list, so it exhausts first; a K
  where its list has run out makes it look flat for a reason unrelated to suppression. (At z = 1.0
  lists run ~15k–26k, so this bites only at z ≥ 2.5 and K ≥ 5000 — but it is checked, not assumed.)
* **Ceiling saturation at K = 250.** `201`, `233`, `013` and `529` carry 17–19 evaluation mitoses,
  so recall@250 can be pinned at or near 1.0 for *every* arm, contributing a structural Δ = 0 that
  dilutes any pooled test. Flagged per ROI and reported alongside, as F4 §7 does.

### 6.2 The crossover K\*
Smallest K in the grid at which a shrunk arm's recall@K exceeds the control's, computed per
`(ROI, seed, axis)` at z = 1.0 and reported as a distribution, with "absent below 10 000" counted
explicitly. This is P3, and it is what "the most TPs in the fewest candidates" actually asks.

### 6.3 Full-list recall — reported, never as a headline
Always printed adjacent to `coverage_frac` for the same cell, with the sentence that a shrink
raises pool size → raises coverage → mechanically raises full-list recall. A full-recall gain that
arrives with a coverage rise is not evidence about suppression.

### 6.4 False positives — three genuinely different questions, kept apart
The user asked whether "our candidate list has fewer FPs overall". That is three questions:

1. **At fixed K:** `FP@K = budget_delivered − tp_at_budget`, split into `lookalike_at_budget`
   (a structure a pathologist examined and rejected) and unannotated. When the list is not
   exhausted — which is every headline cell — `precision@K = tp_at_budget / budget_delivered =
   recall@K · n_mit / K`, an affine restatement of recall@K, i.e. **the same question**, and it
   will be labelled as such rather than presented as independent corroboration. When
   `budget_delivered < K` the denominator differs per arm and the equivalence fails; those cells
   are flagged by §6.1's guard and precision is reported over `budget_delivered`, never over `K`.
2. **Over the full list:** `n_detections − TP − look-alikes`. A shrink **necessarily** increases
   this; the honest number is the `n_detections` ratio to control, per domain. If the question is
   read as "fewer FPs in the whole list", the answer is known in advance to be no, and F5 will say
   so in one line rather than dress it up.
3. **FPs paid to reach a recall target** — the dual of recall@K and the literal form of "most TPs
   per candidate read": `FP_to_r = read_r − ⌈r · n_mit⌉` for r ∈ {0.8, 0.9, 0.95, 0.99, 1.0},
   reported with `n_unreached`. **`FP_to_0.99` and `FP_to_1.0` carry the same coverage guard as
   §6.3** — `read_100` is finite only because the list tiles the ROI, so an unguarded improvement
   there is the same geometry §6.3 rules out, arriving under a different name. `FP_to_0.8` and
   `FP_to_0.9` are the ones that can carry a claim.

Plus **`n_dup_fp` reported with its denominator `n_within_match_radius`** (both as v2 cell 15
records them): a second detection on an already-claimed object, bucketed `non_human_findings`.
This is the defect `invariants.check_nms_radius` exists to prevent, and it rises mechanically with
list length — so the raw count is a list-length artefact without the denominator. F4 §8.5 demotes
this metric for a related reason (it misses duplication on ordinary nuclei, which is most of it);
F5 reports it as a *ratio* and does not let it arbitrate.

### 6.5 The ledger — does a shrink *gain* mitoses, or just shuffle them? And is the dose grid centred?
Greedy NMS is **not monotone in the radius**: a point kept at radius R can be lost at r < R,
because a smaller radius lets an intermediate competitor survive and suppress it. Concretely, and
verified by running `nms_by_distance` itself — A(10), B(9), C(8) with |AB| = 25, |BC| = 15,
|AC| = 30 gives `keep(29.6) = {A, C}` but `keep(20) = {A, B}`. So "shrinking can only help recall"
is **false as stated**, and F5 measures rather than assumes it.

At z = 1.0, per `(ROI, seed)`, two things off one walk of the annotations — the first on both
axes, the second once per radius (see below):

* **Win/loss:** the set of mitotic `ann_id`s found by each arm, differenced against the control →
  **gained**, **lost**, net. A net gain of 5 made of 8 gains and 3 losses is a materially
  different finding from one made of 5 clean gains.
* **The same-object gap:** for each mitotic annotation, the distance from its nearest deep-pool
  local maximum to the peak that suppressed it. The 23–29 px band is measured on **one ROI and one
  seed** (v2 cells 9–11) and now carries the entire justification for where the dose grid sits. If
  the band runs 12–18 px on other ROIs, `r5.0` is above it too and the ladder is mis-centred — a
  fact that would otherwise only be discoverable post hoc. **Computed once per
  `(ROI, seed, radius)`, not per axis**: which peak suppressed which is entirely upstream of
  ranking, so this quantity is axis-free and running it on both axes would duplicate identical
  numbers. The *win/loss* half above stays on both axes, because greedy rank-order matching
  genuinely can hand a contested GT to a different detection (v2 cell 18 measures rather than
  assumes this).

### 6.6 Was the intervention measurable at all?
`topk_churn` — the fraction of the top-250 list that changes between control and each shrink arm —
and `gt_claim_change`, the number of mitotic annotations whose claiming detection changes rank
bucket, both at z = 1.0. Without these, a null at K = 250 is uninterpretable: an inert
intervention and an active-but-neutral one look identical. F4 §6 built a whole pre-flight for this
reason; here it is two columns off tables already being computed.

### 6.7 Carried diagnostics
`coverage_frac`, `largest_tie_block`, `nan_rate`, `floor_limited`, `n_pool`, `n_detections`,
`nms_radius_px`, `match_radius_px`, `read_{50,80,90,95,99,100}`, `n_unreached`,
`n_gt_within_seed_hole` (§7.10).

---

## 7. Verifications and gates — asserted in the run, not claimed in prose

1. `valid.all()` after padding, per `(ROI, seed)`. Padding's provable consequence; a running check.
2. `od_nan == 0` after `OD_PAD` padding, per `(ROI, seed, radius)`, with the naive count recorded
   beside it. This is v2's own bug, which corrupted its cost table by ~116×.
3. **Control passes `invariants.check_nms_radius`** — a *positive* control that the control arm
   really is the repo default. The shrunk arms pass `nms_radius=None` (the check is deliberately
   inapplicable) and record their radius in `extra`, as v2 does.
4. `invariants.check_min_separation(kept_centers, radius)` on **every** arm at its own radius —
   direct proof the suppression that was claimed is the suppression that ran.
5. The one-match-many-`z` shortcut re-verified at **each of the three radii**, at **z = 1.0 and
   z = 3.0**, on **`si == 0` only — 14 ROIs × 3 radii × 2 z = 84 checks, ~34 s**. The shortcut is a
   property of the code, not of the data, so the other 56 cells teach nothing; running all 70 would
   cost ~170 s, ~19 % of the run. This is exactly F4 §9 gate 3's scope. z = 1.0 alone is the weakest possible check — it drops only ~15 % of
   the pool — and the one historical failure of this shortcut was found at z = 2 on 350.tiff
   (`2026-09-01-tm-variant-sweep.md` §7). The argument for the shortcut does hold at any radius
   (`nms.py:49`: `suppressed[idx]` is written only by strictly higher-scoring indices, so a
   point's fate never depends on a lower-scoring one) and survives both `extract_peaks`
   parameters, but it is checked where it could actually break.
6. `check_no_cap` applied to the **pre-NMS pool length** against `MAX_PEAKS`, not only to the
   post-NMS post-filter length the harness checks (`compare.py:172`) — against a 2,000,000 cap the
   latter cannot fail, so on its own it is a check that records a pass it could not have withheld.
7. `check_distinct_seeds` on the 5 RNG streams per ROI.
8. **Nesting is NOT asserted** — see §6.5. A run asserting `keep(r) ⊇ keep(R)` would assert
   something false.
9. **Seed-0 reproduction gate.** `draw_seeds(s=0)` is byte-identical to v2's inline draw (same
   `agreement_pool` → `border_filter(36)` → `default_rng([0, image_id])` → `.iloc[integers(...)]`),
   and `DEEP_FLOOR_Z` is pinned at v2's −1.5, so the `r5.0` arm at seed 0 **must** reproduce
   `results/tm_ccoeff_threshold_axis_sweep_v2.csv` exactly on the 7 overlapping ROIs at every
   shared `(z, budget)`. Free, and it is the only check that can catch a silent divergence in the
   whole pipeline rather than in one function. One mapping is needed: v2's `arm` column is the axis
   (`tm_score` / `chromatin_od`) while F5's is `{axis}@{radius}` (§3d), so the gate compares v2's
   `arm == a` against F5's `arm == f'{a}@r5.0'`.
10. **`self_hit_radius` ordering pinned; the confound it creates is measured, and it is zero.**
    v2 applies the self-hit filter *after* NMS, so the seed's own peak suppresses everything within
    the NMS radius and is then deleted, leaving a detection-free hole of radius ≈ `nms_radius`
    around the seed — 29.6–33.1 px in the control, 19.7–22.1 px at `r5.0`. Part of a measured recall
    gain could then be the hole shrinking rather than suppression relaxing. **Pre-measured before
    the run**, replicating `draw_seeds` exactly over all 14 ROIs × 5 seeds and counting `gt_eval`
    annotations inside the hole at the *largest* radius (r7.5): **0 in every one of the 70 cells.**
    Every seed this design will draw is more than 33.14 px from every other annotation, so the
    confound cannot fire here. F5 keeps v2's order (so the comparison to v2 holds) and asserts
    `n_gt_within_seed_hole == 0` at runtime — which now also guards a future change to the seed
    rule. Revision 2's conditional re-analysis path (recompute the headline on `gt_eval` restricted
    to > 33.14 px) is **deleted**: it was a branch that could not be taken.

Every record lands in `results/f5_nms_radius_ablation_verification.csv`. Read it before any table.

---

## 8. Statistics, fixed before the data exists

**Primary estimand, declared once and bound to the decision rule.**

> **The arbiter is the mean-of-5 paired Δ(recall@250) per ROI**, at z = 1.0 on the `chromatin_od`
> axis: for each of the 70 cells form the *paired* difference (shrink − control, same ROI, same
> seed, same peak pool), average the 5 within an ROI, and treat the resulting **14 values as the
> clustered units**. One test per (dose × axis).

**Why not the worst-of-5, which revision 2 declared.** `worst-of-5(shrink) − worst-of-5(control)`
takes each arm's minimum *independently*, so the argmin seed need not be the same seed in both
arms — which discards exactly the pairing §3b exists to create and re-admits the seed variance the
design was built to cancel. This is not hypothetical; on `results/tm_nms_radius_sweep.csv`,
`read_50`, 29.6 → 20.0 px:

| ROI | per-seed paired Δ | worst **paired** Δ | difference of the two **worsts** |
|---|---|---:|---:|
| `300.tiff` | +3, +20, +4, +6, +6 | **+20** (seed 1) | **+17** (worst A is seed 2, worst B is seed 1) |
| `301.tiff` | +67, +83, +140, +104, +15 | **+140** (seed 2) | **+104** (both seed 3) |

Different argmin seeds, and the two statistics disagree by 3 and by 36 on a five-seed sample.
The bias also has a direction: the minimum of 5 draws is a decreasing function of variance, so at
*identical mean recall* the arm with more seed-to-seed variance gets the lower worst-of-5 — and
the shrunk arm has the longer list and hence the more seed-dependent top-250. So
`min(A) − min(B)` has non-zero expectation under "the radius changes nothing on average". The
sign-flip test would stay valid against the *sharp* null (every cell identical, Δ ≡ 0) while being
read against the *weak* one. The mean-of-5 paired Δ has neither problem.

* **Test: exact cluster sign-flip permutation over the 14 ROIs** (2¹⁴ = 16,384 sign assignments,
  enumerated, not sampled), plus a **14-ROI cluster bootstrap** percentile CI on the mean paired Δ.
  This follows F4 §7 — same repo, same day, same reason: five seeds share an ROI (one tissue, one
  ground truth, one response surface) and two ROIs share a domain, so the 70 cells are not 70
  independent units. The treatment acts on the response surface, which is an ROI property, so the
  intraclass correlation is high by construction and a cell-level test is anti-conservative by an
  unknown factor.
* **Secondary, clustered, and explicitly not the arbiter — the worst-case *paired* reading.** Two
  forms, both of which keep the pairing: `minᵢ Δᵢ` within each ROI ("how bad is the treatment where
  it does least good"), and Δ evaluated **at the control's argmin seed** ("the control is weakest
  on this click; what does the shrink do *there*"). The second is what D4's "worst click" rule is
  actually reaching for. Neither arbitrates, because the minimum of 5 carries no useful uncertainty
  at n = 14.
* **Descriptive reporting frame, no p-value:** §6.1's worst-of-10 recall@250 per domain, per arm,
  with the per-seed paired Δ distribution printed beside it. D4's "worst ROI, worst click, per
  domain" is a *reporting* rule; it is discharged as one.

* **Reported alongside, explicitly labelled as the anti-conservative bound:** the cell-level
  Wilcoxon signed-rank on all 70 paired Δ. It is what revision 1 proposed as primary; it is kept
  only so the gap between the two is visible.
* **Multiplicity: Holm across the 2 arbiter tests** — `r5.9` and `r5.0`, both on `chromatin_od`.
  The family is 2, not 4, because §3d declares `chromatin_od` primary (it is the shipped ranker,
  commit `7c3af93`, and F4's primary axis) and the estimand box above scopes the arbiter to it.
  `tm_score` is reported as a **secondary axis**, uncorrected and labelled as such; it does not
  enter the arbiter's family. Writing this down closes a real forking path: any raw p in
  (0.0125, 0.025] is significant under a family of 2 and not under a family of 4, so "Holm across
  4" and "the arbiter is chromatin_od" cannot both stand.
* **Both secondary statistics are reported unconditionally**, whatever they show, and neither
  carries a corrected p-value into the arbiter's family. They are pre-specified because they answer
  different questions off the same paired Δ matrix — not so that the better-looking one can be
  chosen afterwards.
* **Dedup key: `(file_name, seed_index, arm, z)`** — valid only because §3d embeds the axis in the
  arm name. With v2's naming (`arm` == axis) this key would collapse the two axes and silently
  halve n.
* Effect sizes in **recall points and candidates**, never as ratios of ratios.
* A per-domain depth delta smaller than that ROI's known seed-to-seed SD (260–2744 candidates for
  `read_95`) is reported as **not resolvable at this sample size**, not as a result.

---

## 9. Threats to validity, named now

| threat | handling |
|---|---|
| Coverage saturation makes full recall uninformative | §6.3 — never a headline; always printed with `coverage_frac`. |
| The *same* saturation arriving under another name in `read_99`/`read_100`/`FP_to_1.0` | §6.4.3 — same guard extended to the reach-target family. |
| Crossover missed by a sparse budget grid | §6.1 — 18 budget levels; §6.2 makes K\* itself an output. |
| Control's shorter list exhausts before the largest K | §6.1 — `budget_delivered` guard. |
| recall@250 ceiling-pinned on the four 17–19-mitosis ROIs | §6.1 — flagged per ROI and reported. |
| A null that is really an inert intervention | §6.6 — `topk_churn` / `gt_claim_change`. |
| Seed noise swamps depth deltas | §3c 5 seeds; §8 the "not resolvable" rule. |
| Pseudo-replication from 5 seeds per ROI | §8 — ROI-clustered permutation, n = 14. |
| Domain imbalance from unequal ROI counts | `images/extra_valid` is 2 per domain by construction. |
| The `od` NaN bug that corrupted v2's cost table | §7.2 asserted per cell, naive count recorded. |
| "Shrinking can only help" taken on faith | §6.5 — counter-example verified against the code; measured as a ledger. |
| The dose grid centred on a band measured from one ROI and one seed | §6.5 — the band is re-measured on all 14. |
| The self-hit hole shrinking with the radius, mimicking a recall gain | §7.10 — pre-measured at 0/70 cells; cannot fire; asserted at runtime anyway. |
| An unpaired order statistic re-admitting the seed variance the pairing cancels | §8 — the arbiter is the mean-of-5 **paired** Δ; worst-of-N is descriptive only. |
| A measured value attributed to an arm that was never run | §4 P1 — v2's CSV is the r5.0 arm (`nms_radius_px` 19.74–22.09); the control is unrun. |
| Reading recall@250 as if it settled the FP question | §6.4 — labelled as the same question restated. |
| Silent pipeline divergence from v2 | §7.9 — seed-0 exact reproduction gate. |

## 10. Cost

~9–12 min single-threaded, after the four economies below (the ~15 min figure in revision 1 was
for 3,780 arms with none of them). Measured, not guessed, on `246.tiff` (4933 × 6577):
`load_roi` 0.4 s + `to_channel` 2.6 s **per ROI** (hoisted out of the seed loop);
`fused_response` 1.4 s + `extract_peaks` 0.2 s **per seed**; per radius, NMS 0.2 s +
`score_detections` 0.9 s + `coverage_fraction` 0.16 s + `bucket_detections` 0.06 s.
14 ROIs × 5 seeds × 3 radii × 6 z × 2 axes = 2,520 arms, with `coverage_frac` computed once per
(ROI, seed, radius, z) rather than per arm, `od` once per (ROI, seed) rather than per radius
(§3b), the shortcut check on `si == 0` only (§7.5), and §6.5's gap once per radius rather than per
axis. Those four economies together are worth roughly 6 minutes of a 15-minute run, and none of
them changes a reported number.

## 11. What review changed, and what each change would have cost

| # | change | what the unrevised version would have produced |
|---|---|---|
| 1 | The ladder is a dose–response, not a mechanism decomposition; **P4 withdrawn** | A mechanism claim the geometry cannot support — 0 annotation pairs lie between the two shrink radii. |
| 2 | §2a rewritten; "corrects F4" withdrawn; 5.945 → **5.955 µm** | A false correction of an adjacent pre-registration, and a units error (comparing one ROI's px against another's). |
| 3 | Axis embedded in the arm name; dedup key and `coverage_key` fixed | Half the rows silently dropped from every test, and one radius's coverage reported for all of them. |
| 4 | ROI-clustered exact sign-flip + cluster bootstrap replaces the 70-cell Wilcoxon as primary | An anti-conservative p — plausibly p < 0.05 at n = 70 where the clustered test gives p > 0.1. |
| 5 | One primary estimand declared and bound to the headline | A tested quantity (pooled median) different from the reported one (worst-of-10), with no uncertainty on the headline. |
| 6 | Seed-0 reproduction gate against v2's CSV added; `DEEP_FLOOR_Z` pinned to −1.5 | No end-to-end check at all; and trimming the z grid would have silently moved the deep pool. |
| 7 | Coverage guard extended to `read_99` / `read_100` / `FP_to_1.0` | The saturated regime nominated (§6.4.3) as "where a shrink can genuinely win". |
| 8 | `n_within_match_radius` restored beside `n_dup_fp`; the ratio, not the count, is reported | A list-length artefact read as a duplication result. |
| 9 | Same-object gap re-measured on 14 ROIs; ledger runs on both axes | The entire dose placement resting on one ROI and one seed; the ledger run on the axis that is not shipped. |
| 10 | `self_hit_radius` order pinned, seed-hole confound measured; shortcut checked at z = 3.0; `check_no_cap` moved to the pre-NMS pool; z grid trimmed to 6; §2's "198 → 266" re-cited as the paired seed-0 **198 → 250** | A radius-dependent confound unnamed; a shortcut check at the z least able to fail it; a vacuous invariant; three saturated z levels; and a motivating number that mixed two seeds. |

## 12. What the second review round changed

| # | change | what revision 2 would have produced |
|---|---|---|
| 1 | **§8: the arbiter is the mean-of-5 *paired* Δ per ROI**, not `worst-of-5(shrink) − worst-of-5(control)`. Worst-of-N becomes descriptive; `minᵢ Δᵢ` and "Δ at the control's argmin seed" are secondary, and both stay paired | An **unpaired** primary statistic — a regression introduced while fixing revision 1's estimand defect. It takes each arm's minimum at possibly different seeds (verified: `300.tiff` worst paired Δ = +20 at seed 1 vs difference-of-worsts +17 across seeds 1 and 2; `301.tiff` +140 vs +104) and is biased under the weak null, because the minimum of 5 falls with variance and the shrunk arm has the more seed-dependent top-250. |
| 2 | §6.1 and §8 reconciled: the reporting frame (worst-of-10 per domain, no p-value) and the arbiter (mean-of-5 paired Δ, n = 14 clusters) are now named as **different statistics**, and §4's negative-result rule is re-bound to the arbiter | Two different quantities asserted to be one, with the decision rule pointing at the untested one. |
| 3 | **P1's evidence corrected.** `results/tm_ccoeff_threshold_axis_sweep_v2.csv` records `nms_radius_px` ∈ [19.74, 22.09] — it **is** the r5.0 arm, so its `full_list_recall = 1.0` at z ≤ 1.0 is the shrunk arm's number. The control has never been run; v1 (r7.5, unpadded) was below 1.0 | A prediction registered as "near-vacuous" on the strength of a number belonging to a different arm — understating how live P1 actually is, in a document whose whole purpose is that this arm is unmeasured. Also surfaced the stronger reason to keep z = 2.5/3.0: coverage there is **0.26–0.63**, the only unsaturated region in the design and so the only place P1 can carry evidence. |
| 4 | §3a's band-position column marked **predicted, to be confirmed by §6.5** | §6.5 could only ever confirm the dose placement that the band was used to justify — a circle, from a band measured on one ROI and one seed. |
| 5 | Economies adopted with the above: `od` hoisted to once per (ROI, seed) with F4 gate 7's equivalence assert; §7.5's shortcut check cut to `si == 0`; §7.10's conditional re-analysis deleted (pre-measured 0/70); §6.5's gap computed once per radius rather than per axis; §7.9 arm-name mapping added | ~6 min of a 15-min run spent re-deriving identical numbers, plus one conditional analysis path that could not be taken and one gate that would have failed on a column-name mismatch. |
