# Pre-registration F4: does shrinking the NMS radius with the template rescue largest-CC tightening?

Date: 2026-09-04. **Written before `f4_nms_shrink_sweep.py` is run.**
Status: PLAN, revision 2 — revised after an adversarial review whose findings are folded in
throughout and summarised in §14. Not yet executed.

Follow-up to `Research Logs/2026-09-04-f1-results.md`. F1–F3 are the follow-ups named in
`Research Logs/2026-09-03-tm-axis-sweep-edit-plan.md`; this is the next free number.

> ### Amendment, 2026-09-04 — the primary ranking axis is `tm_score`
>
> **`f4_nms_shrink_sweep.py` does not exist and F4 has not run** — no file of that name is in
> the repo, and no `results/f4_*` output exists. This amendment is prompted by re-analysis of
> F1's *completed* data, not by any F4 data, which is the line that separates it from the
> post-hoc axis switch a pre-registration exists to prevent.
>
> §7.1's filter, §7.2's estimand and §11.1's decision tree now arbitrate on **`tm_score`**;
> `chromatin_od` stays a mandatory companion axis, reported in full and uncorrected. The reason
> is the one **§7.2, §10.9 and §14 already argue** and then decline to act on: `chromatin_od`'s
> justification predates D1's matcher switch, and its advantage under the current configuration
> is +0.032 recall. Two things settle it that this document did not have:
>
> * **`chromatin_od` was never "the ranker the repo ships"** (§7.2's words). `chromatin.rerank`'s
>   only callers are three superseded 2026-08-31 probes; nothing in the pipeline calls it. The
>   inheritance was a documentation error tracing to commit `7c3af93`'s message.
> * **The +0.032 is not an effect.** Clustered at the ROI — F4 §7.3's own unit — it is
>   95 % CI [−0.047, +0.112], p = 0.36, positive on 3 of 7 ROIs. §7.2 read it cell-level
>   (40 of 70), which is the anti-conservative reading §7.3 exists to reject.
>
> **Neither axis won**, so the tie-break goes to the shipped default rather than to a non-default
> post-hoc statistic. §10.9's point — that the axis is *not* neutral with respect to the
> treatment, because §2.2's redundancy mechanism runs through the `od` window — is strengthened
> by this change, not weakened: the arbiter now sits off the mechanism under test.
>
> See `DECISIONS.md` D5 and `verify_chromatin_ranker.py`. D5 also binds any future chromatin
> re-measurement to the 14 ROIs of `images/extra_valid/`, which §5 already uses.

---

## 1. The question, in the terms it was asked

* **Method A — "simple v2".** `tm_threshold_axis_sweep_v2.ipynb`: template is the raw, fixed
  `tm.BASE_SIZE` = 51 px box centred on the click; border padding so the whole ROI is reachable;
  NMS de-duplication at a decoupled **5.0 µm** radius `R₅` (19.7–22.1 px here), separate from the
  7.5 µm evaluation match radius.
* **Method B — "largest-CC + shrunken NMS radius".** The largest-CC size-refined template (Otsu
  the 51 px click crop, take the largest connected component regardless of the click, use its
  bounding box's odd longer side as `base_size = b`) — **and**, new here, the NMS radius shrunk
  in the same proportion as the template.

A vs B is the **primary** contrast. Everything else exists to make that number interpretable.

**Note what is *not* new.** The 7.5 µm → 5.0 µm decoupling is v2's Fix 2, and the largest-CC
notebook inherits it verbatim. Method B is therefore a **second** shrink on top of an existing
one, not the introduction of decoupling.

---

## 2. Prior state, and an honest directional prior

| established | where |
|---|---|
| The 5.0 µm radius fixes two real leaks: unreachable border pixels, and a GT's own on-centre peak eaten by an off-centre **same-object** peak 23–29 px away | v2 cells 8–11 |
| Largest-CC tightening **at the fixed `R₅`** hurts, dose-dependently, crossing zero at `b` ≈ 36 px: below, −2 to −8 recall points; at or above, indistinguishable. ρ = +0.5986, p = 0.0012, n = 35 | `2026-09-04-f1-results.md` |
| Largest-CC's list is longer in **every** informative cell (mean +2,751 detections); `read_95` deeper in 20 of 25 cells (median +355) | same |
| `full_list_recall` is coverage-saturated (0.88–0.999) and discriminates nothing | audit Findings 1–2 |
| Reporting metric is recall@K, worst-ROI/worst-click per domain | `DECISIONS.md` D4 |
| **Minimum GT-to-GT spacing on these 14 ROIs is 26.6 px** (`245.tiff`); max `R₅` is 22.09 px. No GT pair is merged at `R₅` on any ROI | measured during review |
| **A pure NMS-radius sweep at a fixed 51 px template already exists**, on `300.tiff` and `301.tiff` — **both in `extra_valid`** — 5 seeds × radii {29.6, 24, 20, 16, 12} px | `results/tm_nms_radius_sweep.csv` (2026-09-02) |
| **A 150-cell template-size dose-response against recall@250 already exists** on 10 ROIs (**9 in `extra_valid`**) × 5 seeds × 3 bbox methods, sizes 13–51: pooled ρ = **−0.130, p = 0.114** — no monotone relationship, nominally opposite in sign to F1's +0.599 | `results/bbox_headroom_end_to_end.csv` (2026-09-03) |

### 2.0 What the saved radius sweep already settles

`results/tm_nms_radius_sweep.csv` measures **this experiment's own variable** — NMS radius at a
fixed 51 px template — on two main-set ROIs. Re-derived here rather than quoted:

* **Volume: `pool` rises as `r` shrinks in 10 of 10 streams, 40 of 40 adjacent radius steps.**
  `301.tiff` seed 1 goes 16,929 → 45,107 detections from 29.6 → 12 px. H_dup's volume mechanism
  is **established, not open.**
* **Head cost: `read_50` rises (worse) monotonically in 9 of 10 streams** (`300.tiff` seed 0 dips
  135 → 133 at the first step, then rises). `301.tiff` seed 3: 304 → 500.
* **Tail gain saturates at the *first* shrink step, above `R₅`.** Every stream whose `missed`
  moves at all does so between 29.6 and 24 px and is then **flat at 20, 16 and 12** — `301.tiff`
  seeds 0/2/4 go 1/1/2 → 0 at 24 px and stay there; `300.tiff` seed 2 goes 2 → 1. **Nothing below
  24 px buys further ceiling**, and `R₅` = 20 px already sits below that.

The third bullet is §2.1's conclusion — *nothing left to rescue* — as a **measurement on two of
this run's own ROIs at five radii**, not as a geometric argument plus one held-out click. It is
the strongest single piece of prior evidence in this document, and F4 must not re-report it as
new. §8.1's `Δ_radius|51 = D − A` is therefore **partially pre-measured** on `300`/`301`, which
makes it a free external consistency check on arm D — otherwise the least-gated arm in the design.

### 2.1 The rescue hypothesis, and why it is the weaker of the two

**H_rescue.** v2's Leak 2 was an off-centre same-object peak outranking and suppressing the
on-centre one. That gap is a property of the correlation surface, which is a property of the
template; a 25–41 px template should produce proportionally closer secondary peaks, plausibly
inside the 20 px radius calibrated for a 51 px template. Shrinking `r` in proportion would
restore the intended behaviour.

**This is retained as falsifiable, but it is not 50/50, and the plan will not pretend it is.**
The Leak-2 diagnostic was run at the **old ~30 px radius**, and the suppressor was nearest the
*same* annotation in every case. A suppression only costs a ground-truth point when the surviving
peak lands **outside the ~30 px match radius** of the click. Take both branches:

* *If the gap scales with `b`* (the premise): at `b = 19` the gap is ~9–11 px, so the suppressor
  sits ~14–21 px from the click — well inside the match radius. The GT is still claimed. Nothing
  to rescue.
* *If the gap does not scale with `b`* (it is set by nuclear structure, which is what v2 cell 11
  actually measured): it is 23–29 px, already **above** `R₅` = 20 px, so `R₅` never suppresses
  across it. Nothing to rescue.

**The one residual channel H_rescue retains is rank, not claim — and it is signed *positive*,
which makes it the only route by which Method B can win.** At `R₅` the survivor is the *off-centre*
peak, because it suppressed by having the higher `score`. Ranking is by `od`, and an `od` window
centred 23–29 px off the nucleus captures less dense chromatin than an on-centre one, so the
off-centre survivor ranks **deeper**. Admitting the on-centre peak at a smaller `r` therefore puts
a **higher-`od`** candidate into the list, and since `greedy_match` claims in rank order
(`evaluate.py:57-81`), the ground truth is claimed **earlier**. From this channel alone recall@K
can only improve or stay equal. What is unsigned is the **net**, once it is offset against §2.2's
slot consumption. §8.6 measures it directly rather than leaving it as an argument.

(The numbers above are v2 cell 10's stored output, columns `own_peak_dist_to_click` = 4.5, 5.0,
7.1, 8.9, 10.0 px and `gap_own_to_suppressor_px` = 23.3, 25.2, 26.9, 27.3, 29.2 px, with
`same_annotation = True` in all five, at an old radius of 29.6 px.)

### 2.2 The duplication hypothesis, which is favoured

**H_dup.** A smaller radius suppresses less, so the list lengthens, and each extra detection on an
already-claimed object consumes a rank slot at fixed K without adding a true positive. Method B
then loses on both counts — smaller template *and* longer list.

**A mechanism that makes this worse, and that revision 1 missed.** `chromatin_density` averages
the darkest 10% of a **fixed 51 px window** (`midog_utils/chromatin.py:64-88`), regardless of
template size. Two detections 10–18 px apart have heavily overlapping windows and therefore
near-identical `od`. A detection that a smaller radius re-admits beside an already-kept one
enters the `od` ranking at essentially the same rank; `evaluate.greedy_match` is one-to-one
(`evaluate.py:57-81`), so it is bucketed an unannotated FP occupying a slot next to the true
positive it duplicates. **Shrinking the NMS radius below the `od` window cannot remove
`od`-redundancy — it can only manufacture it.** `largest_tie_block` will not catch this: these
are near-tied floats, not exact ties.

**The trade is two-sided, and only one side is visible at fixed K.** §2.0 shows `read_50` rising
(head worse) and `missed` falling (tail better) as `r` shrinks. §8.7 predicts recall@250 falling
monotonically with `eff`, which is the **head** half. That prediction is defensible only because
§2.0 also shows the tail half is **already banked above `R₅`** and, beyond that, would surface
only in `full_list_recall`, which the audit established is coverage-saturated and cannot
discriminate. **§8.8's partial AUC over `cp.BUDGETS ≤ 500` is the K-robustness check that keeps
this from being right for the wrong reason** — the sign of a radius effect is K-dependent, and a
design sampling only K = 250 can land on one side of a crossover and miss it.

**And this mechanism is axis-specific, which makes the choice of ranking axis part of the
treatment rather than a neutral reporting choice.** Under `tm_score`, a re-admitted duplicate is
a *different peak on the correlation surface* and carries a genuinely different score, so it does
not necessarily land adjacent to the detection it duplicates. Under `chromatin_od` it is read
through an overlapping 51 px window and lands adjacent almost by construction. So H_dup predicts
a **larger** penalty on `chromatin_od` than on `tm_score` — a differential prediction the run can
check, and a reason §7.2 makes `tm_score` a mandatory companion rather than a footnote.

### 2.3 Disclosure — a held-out probe run during review, before the main run

An adversarial reviewer ran a **forced-`b`** probe on `405.tiff` (held out; not in
`extra_valid`), one click, template size forced rather than drawn, NMS at `R₅` vs `r(b)` on the
*same* peak set, `z = 1.0`, `od`-ranked, K = 250, `n_gt` = 13:

| `b` | `r(b)` px | pool @`R₅` | pool @`r(b)` | ratio | top-250 **set** diff | tp@250 | redundant@250 |
|---:|---:|---:|---:|---:|---:|---|---|
| 19 | 7.51 | 21674 | 29763 | 1.373 | 48 | 10 → **9** | 9 → **15** |
| 29 | 11.46 | 14093 | 15573 | 1.105 | 21 | 12 → **11** | 6 → **8** |
| 39 | 15.41 | 11526 | 11859 | 1.029 | 7 | 12 → 12 | 3 → 4 |
| 45 | 17.78 | 10693 | 10816 | 1.012 | 6 | 12 → 12 | 2 → 3 |
| 51 | 20.15 | 9781 | 9781 | 1.000 | 0 | 12 → 12 | 1 → 1 |

One ROI, one click, one domain, forced doses — **not evidence about the methods**, and cited only
for its direction and its mechanism, never its magnitude. It establishes three things: the
intervention **is** measurable (revision 1's motivating fear was wrong); it is measurable **only
at small `b`**; and where measurable it is **harmful**, via exactly the redundancy of §2.2.

**Consequence for this document.** The prior favours H_dup. The run still earns its keep — the
user asked for a *quantified* difference, the magnitude is unknown, arm C on the 7 new ROIs is a
genuine replication, and §8.7's radius grid turns a single rule into a dose curve for the price
of a few evaluations.

---

## 3. The shrink rule, fixed now

```
r(b) = radius_px(mpp, 5.0) * (b / 51.0)          # NOT R5 * b / 51.0 -- see §9 gate 2
r(51) == R5                                       # asserted at construction
```

Equivalently `5.0 × b/51` µm, reported in both units.

**Justification, corrected.** Revision 1 claimed the rule "holds *radius / template side*
constant at 20/51 ≈ 0.39, which is what Method A already uses." **That provenance was false.**
v2 set `NMS_RADIUS_UM = 5.0` as a *physical* radius chosen to sit below a **measured 23–29 px
same-object gap**; that it also equals 0.39 × 51 px is a coincidence of two independently-set
quantities. The defensible argument is the one that was missing: the cross-correlation main lobe
of a `b`-wide template is ~`b` wide, so under linear lobe scaling `gap(b) ∝ b`, and
`r(b) = R₅·b/51` preserves `r/gap` — which *is* what the 5.0 µm calibration fixed.

**Secondary properties, useful but not the argument:** it degenerates to exactly `R₅` at `b = 51`
(so a null cell is a four-way identity — a gate convenience), and `b ≤ 51` always because the
Otsu crop is 51 px wide.

**Alternatives considered and not taken:**

* `r = b/2` or `r = b` px — the standard "no two detections closer than one template width". At
  `b = 51` this gives 25.5 px ≠ `R₅`, so it does *not* reduce to Method A; the degeneracy
  property is a gate convenience, not a reason. Named as the road not taken.
* **Measure `gap(b)` per template size on the held-out ROIs** and set `r` from it. Revision 1
  rejected this as circular; that was wrong — it is circular only if measured on the main ROIs,
  and §6 already establishes a held-out probe set. It is the only rule with a directly measured
  mechanism, and §6.3 now **measures `gap(b)` as a by-product** so the `b/51` scaling can be
  checked rather than assumed. The rule itself is not changed on the strength of it; that would
  be a different experiment.
* Floor `r` at the peak-separation limit — rejected as it would put a discontinuity exactly where
  the interesting small templates live. The degenerate region is flagged and dose-graded instead.
* Area scaling `r ∝ (b/51)²` — rejected: NMS is a distance operation. Area is retained only as
  the dose covariate for response *magnitude*, per F1.

### 3.1 The degeneracy threshold is 8 px, not 7, and it bites inside prior experience

`template_match.extract_peaks` masks `fused >= dilate(fused, 15×15)`, an L∞ ≥ 8 constraint on
distinct-valued peaks. **Measured** on `405.tiff`: peaks with a nearest neighbour under 8 px are
0 of 21,745 (`b`=45) and 2 of 18,847 (`b`=51) — 0.01%. (The "zero by construction" form is
*false*: minimum separation came back 1.414 px on tie plateaus. This is stated as measured.)

So NMS below ~8 px is a no-op to within 0.01% of the peak set. `r(b) < 8` at **`b < 18.5`**
(human, mpp 0.226–0.230) to **`b < 20.7`** (canine, mpp 0.2533). **F1 drew `b = 19` twice, both on
`301.tiff` (mpp 0.2533 → `r` = 7.35 px), an ROI in `extra_valid`.** The degenerate region is
therefore *inside* prior experience, not outside it as revision 1 claimed, and probe row `b`=19
confirms arm B there runs with NMS effectively disabled.

### 3.2 Effective suppression fraction — the real dose

The annulus in which suppression can actually act is `r² − 64` against `R₅² − 64`:

```
eff(b) = max(0, r(b)^2 - 64) / (R5^2 - 64)
```

Over F1's realised sizes: `b`=19 → **0.00**, 25 → 0.12, 29 → 0.22, 31 → 0.27, 33 → 0.31,
37 → 0.44, 41 → 0.58, 47 → 0.82, 51 → 1.00. **7 of F1's 35 cells received under 25% of the
nominal treatment.** `eff(b)` is recorded per cell and is the dose covariate for the radius
factor.

**It does *not* break the collinearity, and revision 2 was wrong to say so.** `radius_ratio = b/51`
and `area_ratio = (b/51)²` are both deterministic functions of the single variable `b`, related by
a bijection on (0, 1]; non-affinity is not identifying variation. **Within arm B the two doses are
perfectly collinear**, exactly as revision 1 said. The only genuinely new variation `eff(b)`
introduces is that it depends on `R₅` and hence on **mpp**, so at fixed `b` it differs by scanner
(`b` = 31: canine 0.251, human 0.275) — real, tiny, and **fully confounded with scanner and
therefore with domain**, so unusable for identification. The identifying variation comes from arms
C and D and from §8.7's grid. That is the honest argument for all three.

**Implied radii** (`r(b) < 8` shaded by the flag, not by a floor):

| `b` | canine (mpp 0.248–0.253) | human (mpp 0.226–0.230) |
|---:|---:|---:|
| 51 | 19.7–20.2 px | 22.0–22.1 px |
| 41 | 15.9–16.2 | 17.7–17.8 |
| 31 | 12.0–12.3 | 13.4 |
| 19 | 7.4–7.5 **(degenerate)** | 8.2 |

---

## 4. Arms — a 2×2 factorial

Two template matches per cell (`b = 51` shared by A and D; `b = lcc` shared by B and C); only NMS,
`od` subsetting and evaluation are repeated.

| arm | `base_size` | NMS radius | role |
|---|---|---|---|
| **A** `base51` | 51 | `R₅` | **Method A.** Reproduces v2 / F1's `base51`. |
| **B** `lcc_shrunk` | `b` | `r(b)` | **Method B.** The comparison asked for. |
| **C** `lcc` | `b` | `R₅` | Decomposition. Reproduces the largest-CC notebook / F1's `largest_cc`. |
| **D** `base51_shrunk` | 51 | `r(b)` | Decomposition. Counterfactual only. |

**C and D are not candidate methods.** D in particular is not deployable — its radius is set by a
tightening step it does not perform. They exist to split `Δ(A→B) = Δ_template + Δ_radius +
interaction` along both paths. The radius factor is **dose-matched within cell** between `(B−C)`
and `(D−A)` (same `b`, same `r(b)`), so the 2×2 is identified and `Δ(A→B) = (C−A) + (B−C)` is
exact with no residual.

---

## 5. Data and design

* **All 14 ROIs in `images/extra_valid/`** — `n_mitotic ≥ 15` and border-filtered unanimous seed
  pool ≥ 5; byte-identical APFS clones of `images/`, so `load_roi` / `roi_mpp` are identical and
  results are directly comparable with everything saved.
* **5 seeds per ROI**, `si ∈ 0..4`, exactly F1's rule: `rng = np.random.default_rng([si,
  image_id])` over `agreement_pool → border_filter(36)`, then `draw_seed_with_retry` with
  `largest_cc_box` as the check. One seed per cell, **shared by all four arms**.
* **70 cells** (14 × 5) × 4 arm-runs.
* Seeds drawn **with replacement across indices** (F1's rule and rationale: collisions are iid
  draws, not defects — `invariants.check_distinct_seeds`). `201.tiff`'s pool is 7 and
  `013.tiff`'s is 11, so collisions are expected on the sparse ROIs.

**Overlap with F1, stated up front.** `extra_valid` contains **all 7 of F1's ROIs**
(`301, 201, 246, 459, 094, 548, 402`), and the seed draw is deterministic in `image_id`, so on
those 7 this run **redraws F1's exact 35 cells**. That is what makes GATE 1 work, and it means:

* the primary A-vs-B contrast is genuinely new (arm B is new everywhere);
* **`Δ_template = C − A` is 50% recycled from F1** and must never be reported as independent
  confirmation of it. See §8.2.

**Overlap with F5, which is running concurrently on the same lever.**
`Research Logs/2026-09-04-f5-preregistration.md` pre-registers *"what does shrinking the NMS
radius do, with border padding held fixed?"* — arms at 7.5 / 5.9 / 5.0 µm on a fixed 51 px
template. That is F4's A↔D factor and §8.7's grid, run **upward** from `R₅` where F4 runs
**downward**, on an overlapping ROI set. Both use `default_rng([si, image_id])` over the same
`agreement_pool → border_filter(36)`, so the cells coincide. **Verified, not inferred:** F5's
smoke run (`results/f5_nms_radius_ablation.csv`) drew `013.tiff` `si` = 0 → `seed_ann_id` **254**;
F4's rule draws **ann 254, `b` = 31, 0 retries** on the same ROI and index.

**How much actually overlaps, stated precisely rather than as a blanket.** F5 varies the radius at
a **fixed 51 px** template. So:

| F4 component | template | overlaps F5? |
|---|---|---|
| arm D (`base51` + `r(b)`) | 51 px | **Yes** — same factor, same cells, different radius levels |
| §8.7 grid at **null** cells (`b` = 51) | 51 px | **Yes**, on ~29% of cells |
| §8.7 grid at **informative** cells | largest-CC `b` | **No** — F5 never runs a tightened template |
| arms A, B, C | — | **No** — A is F5's baseline, B and C are tightened |

> **Pre-committed:** where the table says Yes, F4 and F5 are **never** reported as independent
> evidence about the NMS radius — they are the same cells under two parameterisations of one
> factor, and whichever document is written second cites the first and reports the union, not the
> sum. Where it says No, F4 is measuring something F5 does not: the radius on a *tightened*
> template, which is the whole point of §8.7 and cannot be recovered from F5 at any radius.
>
> Arm D is therefore retained (it is 4 evaluations, ~10 min, and it is the D→A path the 2×2 needs)
> but is reported as **corroborating F5 and `tm_nms_radius_sweep.csv`, not as new evidence**.

**On worst-of-10.** 2 ROIs per domain makes D4's "worst ROI *and* worst click" a worst-of-10
within this run. It is **not** an improvement over F1's worst-of-5 in the comparable sense —
worst-of-K is biased downward in K, so the two are not on the same scale. It is a fair
**within-run, paired** A-vs-B contrast (identical seeds), and that is all it is claimed to be.

### Held fixed across every arm and cell

Channel `hematoxylin_od` (D3), `TM_CCOEFF` (D1), no `tissue_mask` (D2), `scales=(1.0,)`,
`n_angles=1`, `flips=(False,)`, `peak_min_distance=7`, `self_hit_radius=5.0`, `patch_size=73`,
`BORDER=36`, `DEEP_FLOOR_Z=-1.5`, `MAX_PEAKS=2_000_000`, `Z_LEVELS=(0.5,1.0,1.5,2.0,2.5,3.0)`,
`OTSU_WINDOW=51`, `LARGEST_CC_KW=dict(min_area=50, max_area_frac=0.85, min_solidity=0.5)`,
`OD_PAD=25`, match radius `radius_px(mpp)`, `cp.BUDGETS`, both axes.

`Z_LEVELS` is F1's, a subset of both saved notebook grids — that is what lets GATE 1 work, and
why it is not re-chosen.

### Varying, and only this

`base_size ∈ {51, b}` × `nms_radius ∈ {R₅, r(b)}`, plus §8.7's radius grid at `z = 1.0` only.

---

## 6. Step 0 — a forced-dose pre-flight on held-out ROIs

**Why revision 1's version would have failed — corrected, because the first version of this
paragraph was itself wrong.** Revision 2 asserted that a drawn probe "would very likely" median to
zero because `405.tiff` draws four nulls in five. **Measured on all four probe ROIs, that
generalisation is false:**

| ROI | pool | `b` by `si` 0–4 |
|---|---:|---|
| `202.tiff` | 4 | 37, 35, 27, 27, 35 |
| `405.tiff` | 10 | **51**, 45, 51, 51, 51 |
| `506.tiff` | 8 | 43, 39, 31, 43, 43 |
| `002.tiff` | 7 | 35, 29, 27, 37, 29 |

Revision 1's 2-seeds × 4-ROIs probe would have drawn `b` = {29, 35, 35, 37, 39, 43, 45, 51} —
**1 of 8 null, not "most"**. The four-null claim is true of `405.tiff` alone, which is the only ROI
either party had measured when the claim was made.

**The real defect is dose coverage, and it is stronger.** Those eight draws span `b ∈ [29, 45]`,
i.e. `eff(b) ≈ 0.22–0.74`. **Not one of them reaches the regime the experiment is about** — `b ≤ 27`
where F1's penalty lives, `b ≈ 19` where §3.1's degeneracy bites, or the doses at which §2.3's
probe actually found an effect (set difference 48 at `b`=19, 21 at `b`=29, 7 at `b`=39). A drawn
probe does not merely risk nulls; it **systematically fails to sample the dose range under test**,
and adding seeds cannot fix it because the click distribution is what it is. That is what forces
§6.1's design.

**6.1 Design.** Forced `b ∈ {19, 25, 31, 37, 43, 51}` at one drawn click per ROI. One
`fused_response` + `extract_peaks` per `b`, then NMS at `R₅` and at `r(b)` on the *same* peak set.
This removes the draw as a source of dose variation entirely. ~15 s per (ROI, `b`).

**6.2 ROIs.** `202.tiff` (16 mitoses), `405.tiff` (14), `506.tiff` (12), `002.tiff` (9) — all
**outside** `extra_valid` and outside the main run. They are below the `n_mitotic ≥ 15`
decision-grade bar and are **not evidence about the methods**. They span only **4 of 7 domains**
(canine lung, canine soft tissue sarcoma, human melanoma, human breast); mast cell,
lymphosarcoma and neuroendocrine — including the two densest ROIs in the main set — are unprobed.

**6.3 Measured.**

1. `topk_setdiff` — **set** difference of the top-250 candidate positions between the two radii.
   Defined as a set, not positionally: one insertion at rank 10 shifts 240 positions, which would
   make a positional threshold vacuous. The set difference also bounds the effect directly:
   `|Δrecall@250| ≤ topk_setdiff / n_gt`.
2. `redundant@250` — top-250 detections inside a GT match circle already claimed by an
   earlier detection (§2.2's mechanism, measured).
3. `tp@250` and `n_detections` ratio.
4. **`gap(b)`** — over GTs whose nearest peak is *actually suppressed* at `R₅`, the distance from
   that peak to its suppressor, per `b`. Reported as a **distribution with the count of such
   GTs**, plus the fraction of GTs having any same-object secondary peak at all. Defined this way
   because v2 cell 10 sets `_suppressor_of[idx] = idx` for kept peaks, so a median over *all* GTs
   is ≈0 by construction; the v2 diagnostic computed it only over the 5 mitoses actually missed.
   This tests §3's linear-lobe scaling assumption on held-out data. **If `gap(b) < 8` at small
   `b`, no radius can act on it** — which is itself the answer to §3's question.

`gt_claim_change` is **dropped** as a gate criterion: over ~8 cells its median is 0 under almost
any effect ≤ 0.5 GT/cell, so it cannot discriminate.

**6.4 The branch — one statistic, exhaustive rows, in §11's units.**

> **Gate statistic**, defined per dose: `S(b) = median over the 4 probe ROIs of
> (topk_setdiff / n_gt)`, evaluated at each informative `b ∈ {19, 25, 31, 37, 43}`.

It is expressed as a **fraction of `n_gt`, not a raw count**, for two reasons: §6.3's bound
`|Δrecall@250| ≤ topk_setdiff / n_gt` is already in those units, and the probe ROIs have
`n_gt` = 8–15 against the main set's 17–238, so a fixed count of 5 would mean `|Δ| ≤ 0.33` on a
probe ROI and `|Δ| ≤ 0.02` on `548.tiff` — three orders of severity apart. The threshold is
§11's own resolution threshold, 0.02.

| outcome | action |
|---|---|
| `S(b) ≥ 0.02` at **any** informative dose | Run §4 as specified. Primary is A vs B. |
| `S(b) < 0.02` at **every** informative dose | **The pre-flight is the answer.** Do not run a 65-minute A-vs-B that §11 would read as "resolution-bounded". Report that the `b/51` rule cannot produce a difference above §11's resolution threshold at any achievable dose, and promote §8.7's radius grid to primary in its §8.7.1 form. |

The two rows are exhaustive and mutually exclusive on one statistic — revision 2's were not (a
per-dose profile of 0.08, 0.02, 0.015, 0.01, 0.008 satisfied neither its pooled-median row nor its
every-dose row, leaving the analyst to choose, which is the forking path §6 exists to close).

**The threshold is not tunable, and saying so is the honest derivation.** Probe `n_gt` after seed
removal is 15 (`202`), 13 (`405`), 11 (`506`), 8 (`002`), so the smallest attainable non-zero
`S(b)` is 1/15 = **0.067 — more than 3× the 0.02 threshold**. Row 2 therefore fires **iff the
median `topk_setdiff` is exactly 0 at every informative dose.** That severity is right, because
row 2's consequence is drastic (abandon the primary contrast) — but the `n_gt` normalisation is
then doing no work while appearing to, so it is stated in §11's units only to make the transfer
visible, not because it is a tuned quantity. **Raw `topk_setdiff` is reported alongside**, with
the bound it implies on the main set's sparsest (÷17) and densest (÷238) ROI, which is the number
that actually transfers.

**Necessary, not sufficient.** The bound is a *ceiling*, so a large `S(b)` proves only that an
effect is not excluded; and the probe ROIs' `n_gt` is 1–2 orders below the dense main ROIs. Row 1
licenses running the experiment; it licenses no expectation about the result.

**6.5 Disclosure.** Revision 1 claimed "no cell in the main analysis is inspected before the main
run." That was **false**: the 21,694-peak figure and the whole of §13's compute budget come from
`460.tiff`, which **is** in `extra_valid`. It is a list-mechanics and timing measurement, not an
outcome measurement, but the claim is corrected rather than repeated.

**6.6 Corrections to this document's own review.** The reviewer's round-1 finding that the drawn
pre-flight would median to zero was over-generalised from `405.tiff`, the only ROI measured at the
time, and revision 2 adopted it verbatim. Measuring the other three refuted it (§6 opening). The
remedy — forced `b` — survives on the dose-coverage argument, which is the stronger one. Recorded
here rather than silently softened, per this repo's standing practice.

---

## 7. Primary analysis

### 7.1 The frame, and the assert that protects it

`DECISIONS.md:318-323` documents the trap verbatim: the results CSV emits one row per budget, so
each (ROI, click, method) appears 8×, and *"any row-level count, p-value or correlation computed
without deduplicating on `(file_name, seed_index, arm)` silently multiplies the sample size by
8."* F1 hit a related `pivot_table` cartesian-product bug mid-run. F4 adds `arm_name` with **four**
levels on top of that, giving 4 × 6 z × 8 budgets × 2 axes × 70 cells = **26,880 rows** in the
factorial file.

**§8.7's radius grid must not land in this frame.** It is evaluated at `z = 1.0` on both axes, and
`evaluate_arms` emits one row per budget regardless of what varies, so it would add
70 × 5 levels × 2 axes × 8 budgets = **5,600 rows** whose `arm_name` is neither of the four — and
§7.1's pivot would silently return **9 arm columns instead of 4** (4 factorial + 5 grid levels). It therefore goes to its
**own file**, and the filter names the four arms explicitly rather than relying on the file split
alone.

> **Construction, fixed now.** Read `results/f4_nms_shrink_sweep.csv` only. Filter
> `z == 1.0 & budget == 250 & arm == 'tm_score' &
> arm_name.isin(['base51','lcc_shrunk','lcc','base51_shrunk'])`, pivot on `arm_name`.
> **Assert** the frame has exactly **70 rows**, over 14 × 5 unique `(file_name, seed_index)`, with
> exactly **4** non-null arm columns — *before* any statistic is computed. The same assert runs for
> `arm == 'chromatin_od'` (§7.2's companion axis). *Amended: the two axis names are swapped
> relative to revision 3; see the amendment at the head of this document.*

### 7.2 The estimand

> **Mean paired `Δ = recall@250(B) − recall@250(A)`, at `z = 1.0`, ranked by `tm_score`, over
> all 70 cells.**

Every coordinate is inherited, not chosen on this data: `recall_at_budget` (`compare.py` primary,
D4), **`K = 250`** and `z = 1.0` (verbatim from F1's pre-registration). **The axis is the one
coordinate that is *not* inherited** — F1 pre-registered `chromatin_od`, and the amendment at the
head of this document replaces it with `tm_score` for reasons the next two paragraphs already
give. `chromatin_od` remains a mandatory companion axis (§10.9).

**On inheriting `K`.** It is the right anti-forking discipline, but honesty requires naming what
does *not* transfer: F1 justified `K = 250` on 7 dense ROIs ("0 of 7 domains saturated"). This
dataset adds four ROIs with **17–19 GT after seed removal** (`201`, `233`, `013`, `529`), where
the **minimum non-zero |Δ| is 1/17 = 0.059 — 4.5× F1's pooled effect of −0.013**. Those four
(2 of 7 domain clusters) can contribute only exact zeros or large outliers. `K` is **not**
re-picked; instead §7.4 adds a granularity-free companion.

**On inheriting the *axis*, which is a weaker inheritance than it looks.** `chromatin_od` is the
ranker the repo ships (commit `7c3af93`) and F1's pre-registered primary, so inheriting it is the
same anti-forking discipline as inheriting `K`. But its evidence base was measured under a
**different matcher**: `chromatin.py`'s own docstring justifies it as *"the signal
`TM_CCOEFF_NORMED` is mathematically blind to"*, and `DECISIONS.md` **D1 replaced that matcher
with `TM_CCOEFF`, which reads contrast** — D1 says so itself (*"re-ranked afterwards by chromatin
density, the two methods are a dead heat"*). Re-measured on F1's saved run under the current
configuration (recall@250, `z` = 1.0, 70 arm-cells): `chromatin_od` 0.6403 vs `tm_score` 0.6081,
**+0.032 mean, better in only 40 of 70 cells and 3 of 7 domains**. And on **`read_50` — the metric the original *"median 19.4×, 7/7"* was measured on
(`Research Logs/2026-08-31-chromatin-density-rerank.md:121`) — the two available measurements
under the current matcher *disagree*, so the advantage is not established either way:

| source | structure | `read_50(od)/read_50(tm)` | od better in |
|---|---|---:|---|
| `results/f1_seed_sweep.csv`, z = 1.0 | 7 ROIs × 5 seeds, median | **1.021** (IQR 0.841–1.202) | 33 of 70 cells, **3 of 7** domains |
| `results/tp_fp_reading_depth.csv`, z ≥ 1.0 | 7 ROIs, single pool | **0.635** | **5 of 7** ROIs |

They agree in sign on 5 of 7 ROIs but flip on `094` (0.635 vs 1.194) and `548` (0.883 vs 1.023),
and diverge in magnitude on `459` (0.269 vs 0.892). Single-seed pool against 5-seed median is the
likely cause — which is itself the point.

So, stated with the orientation explicit and without picking a side: **on recall@250 the advantage
survives but is marginal and domain-inconsistent; on `read_50` the two artifacts disagree on
direction and magnitude, and the 19.4× does not reproduce in either.**

> **Concurrency caveat.** `results/tp_fp_reading_depth.csv` was written at 17:28 on 2026-09-04 by
> a session running in parallel with this one and is mid-flight. Re-read before the run; the
> conclusion drawn from it (carry both axes) does not depend on its exact values.

> **Therefore `tm_score` at the identical coordinates is a mandatory companion axis**, on the
> same "hold on both or report as split" rule as §7.4.1. Both axes come free off the same pool —
> the run already emits them — so this costs nothing but the discipline of saying in advance that
> the conclusion must survive the axis choice. F1 supports the expectation that it will:
> *"tightening hurts on both rankers; only the monotone structure is specific to
> `chromatin_od`."* §2.2 gives the reason it might not.

### 7.3 Tests

* **Primary p:** exact cluster sign-flip permutation, two-sided, statistic = mean Δ. Cells within
  an ROI share tissue and ground truth, so signs flip **per cluster**. Reported at **both**
  clustering levels, with the conclusion required to hold at the conservative one:
  * ROI-clustered, 2¹⁴ = 16,384 assignments, **enumerated exactly**;
  * **domain-clustered, 2⁷ = 128 assignments** — the conservative version, since `extra_valid`
    is 7 deliberate domain *pairs* and F1 measured per-domain means spanning −0.044 to +0.012.

> **Attainable-minimum-p check, run and reported before any p-value is interpreted.** A cluster
> whose Δ sums to exactly 0 contributes nothing under either sign, so the exact test's support is
> `2^(n_nonzero clusters)`, **not** `2^(n clusters)`, and the minimum attainable two-sided p is
> `2 / 2^n_nonzero`. **Measured on F1's own data** (Δ = C−A, z = 1.0, K = 250, `chromatin_od`),
> `201.tiff`'s cluster sum is **exactly 0.00000 with 0 of 5 cells non-zero**, giving 6 non-zero
> clusters and `p_min` = 2/2⁶ = **0.031** — no margin. One more zero cluster gives
> `p_min` = **0.0625 > α**, at which point **the domain-clustered test cannot reject at 0.05 no
> matter how large the effect is.**
>
> The risk is identifiable in advance: canine lung pairs `201.tiff` and `233.tiff`, both with 17
> GT after seed removal — the **only** domain in `extra_valid` whose members are both sparse,
> where every other domain has one dense ROI, and `201.tiff` was identically zero throughout F1.
>
> **Rule.** Report `n_nonzero_clusters` and `p_min` at both levels first. If `p_min > α` at the
> domain level, that level is **declared uninformative by arithmetic**, the ROI-clustered test
> becomes the conservative-available one, and the fact is stated in the results. A non-rejection
> that was forced by the support size is never reported as evidence of no effect. The
> domain-clustered test is additionally run on §7.4.2's `Δ_count`, which is far less likely to be
> exactly zero on a sparse cluster.
* **Primary CI:** obtained by **inverting the exact sign-flip test** — no bootstrap assumption,
  and the natural companion to the p-value.
* **Secondary CI:** BCa cluster bootstrap (not percentile: 14 clusters with a zero-inflated,
  granularity-heterogeneous Δ gives materially sub-nominal percentile coverage). Report the
  **number of clusters with a non-zero Δ**; clusters that are identically zero contribute no
  resampling variability, so the effective cluster count is below 14 (F1's `201.tiff` was zero in
  all 5 cells at K = 250).
* **Cross-check:** Wilcoxon signed-rank.

### 7.4 Companion estimands, both mandatory

1. **Informative-cell mean** (`b ≠ 51`), same machinery. Null cells are ~29% structural zeros
   (10/35 in F1) that mechanically shrink both the point estimate and the CI; retaining them is
   right for the deployment estimand (§10.7) and wrong for an equivalence claim. **The §11
   conclusion must hold on both the pooled and the informative-cell estimand, or be reported as
   split.**

   *Why this needs no multiplicity correction.* Requiring the conclusion on **both** is an
   **intersection–union test**: the operative statistic is the *maximum* of the two p-values, and
   its size is bounded by α with no correction. (The design that would need one is the opposite —
   where *either* estimand suffices — which is union–intersection.) The two estimands are
   moreover **nested**, informative cells being a subset of all cells, so their p-values are
   strongly positively dependent and the rule is only mildly conservative; the "split" row should
   fire rarely.
2. **Per-mitosis, GT-weighted count statistic**, granularity-free:
   `Δ_count = (Σ_cells tp@250(B) − Σ_cells tp@250(A)) / Σ_cells n_gt_mitotic`.
   Computed from rows the run already emits; immune to the 1/17 quantisation of §7.2. **Named
   "GT-weighted" because it is:** `548.tiff` (238 GT) carries 14× the weight of `201.tiff` (17),
   which cures the granularity problem but partly undoes D4's per-domain concern by construction.
   The **per-domain version is reported beside the pooled one** for exactly that reason.

Reported beside every headline: `n_null`, `n_degenerate` (`r(b) < 8`), mean `eff(b)`,
`Δn_detections`, `Δcoverage_frac`, `Δread_95`, `Δredundant@250`.

### 7.5 Estimation primary vs. the D4 decision frame — they are different questions

D4 (`DECISIONS.md:296-311`) is explicit that a median/mean **hides the failure it cares about**
and mandates worst-case per domain. This plan's primary is a pooled mean — the statistic D4
rejects as a *decision tool*. Both are correct for their own purpose and the plan states so
rather than eliding it:

* **Estimation primary** (§7.2): the user asked to *quantify the difference*. That is an
  estimation target, and a mean with an interval is its answer.
* **D4 decision frame** (§11.2): per-domain **worst-of-10** recall@250 for A and for B, with its
  **own** pre-registered decision rule. No §11 conclusion may be phrased as "on the metric D4
  designates" unless it is derived from this frame.

---

## 8. Secondary analyses

**Multiplicity, stated.** α = 0.05 two-sided applies to the §7 primary alone. Everything in §8 is
**estimation, not testing**: intervals and effect sizes, with p-values reported as descriptive and
uncorrected, and no §8 result may on its own change the §11 verdict.

1. **Decomposition.** Δ_template = C − A; Δ_radius|lcc = B − C; Δ_radius|51 = D − A;
   interaction = (B − C) − (D − A).
2. **Named deliverable — replication of F1, on the 7 *new* ROIs only.** *Promoted out of the
   estimation-only banner:* this is the **only external validation of F1's headline that exists
   anywhere**. F1's ρ = +0.5986 is a 7-ROI, one-region-per-domain result whose strictest subset
   did not clear 0.05 (p = 0.1535 at n = 17), and `Δ_template = C − A` on 7 untouched ROIs is a
   genuine out-of-sample test of it. It gets its own section in the results document and its own
   pre-registered reading: **replicated** if the sign and the sub-/supra-36 px split both hold on
   the new 35 cells; **downgraded to a 7-ROI finding** if they do not.

   **"Replicated" means same channel, same ranker, same radius regime, same tightening rule — and
   that qualifier is load-bearing.** `results/bbox_headroom_end_to_end.csv` is a 150-cell
   size-vs-recall@250 sweep on nine of these fourteen ROIs and returns ρ = **−0.130, p = 0.114**,
   nominally opposite to F1. It does **not** contradict F1, because it differs on all four:
   `gray_inverted` (D3: not comparable to `hematoxylin_od`), NMS at the **match** radius
   (29.6–33.1 px, pre-v2), click-aligned `tighten_box_otsu` variants (`binary` / `headroom` /
   `multiotsu`) rather than `largest_cc_box`, and match-score ranking — the `tm_score` analogue,
   where **F1 itself found ρ = +0.164, p = 0.399**. A null there is what F1 predicts, so this is a
   third independent line of evidence that the dose-response is **specific to `chromatin_od`**,
   which is exactly what §10.9 flags. Without the qualifier, §8.2's verdict is over-readable in
   both directions. Details — `300, 233, 245, 460, 013, 529, 403`, 35
   cells. Revision 1 called the 70-cell version "an external replication on 7 new regions"; that
   was **wrong**, because 35 of those 70 cells *are* F1's cells. Reported separately from the
   pooled version, which is explicitly labelled "includes F1's cells".
3. **Does F1's ~36 px crossover hold?** Split at `area_ratio` 0.50 — F1's split, not a new search
   — on the **7 new ROIs**, for the same reason as §8.2.
4. **Dose–response.** Spearman ρ of paired Δ against `area_ratio` (template dose) and separately
   against `eff(b)` (radius dose), cluster sign-flip.
5. **Border sensitivity.** `PAD = (b−1)/2` differs by arm (25 px at `b`=51, 9–20 px tightened), so
   near the ROI edge arms score against `BORDER_REPLICATE` fabrication of different widths —
   named in `f1_seed_sweep.py`'s docstring and never measured. Recompute recall@250 against
   `gt_eval` restricted to annotations ≥ 25 px from the edge. If the conclusion survives, the
   caveat is retired.
6. **Mechanism, led by `Δn_detections` (C→B)** — the direct net readout of the radius change,
   valid precisely because C and B share `med`/`mad` (§10.1). **Not** led by `n_dup_fp`, which
   counts duplicates only inside *annotated* match circles and so misses duplication on ordinary
   nuclei, which is most of it; it is reported as the annotated slice it is.

   **Plus the rank channel, measured** (§2.1): **median `matched_rank` per claimed GT**, taken from
   `bucket_detections`' `gt_out` (`evaluate.py:114-137`), paired A-vs-B and C-vs-B. A negative
   shift in median `matched_rank` with no recall gain is the positive rank channel firing and
   being cancelled by slot consumption — a reportable finding rather than a hand-wave, and the
   only direct evidence either way on H_rescue's surviving mechanism.
7. **Radius dose grid**, spaced in `eff`, not in `R₅`. On the largest-CC template, at
   **`z = 1.0` only**, both axes, written to its own file (§7.1):

   ```
   r(eff) = sqrt(64 + eff * (R5^2 - 64))     for eff in {0.25, 0.50, 0.75, 1.00}
   plus  eff = 0  ->  NO NMS AT ALL (nms_radius_px = 0.0, no_nms = True)
   ```

   The `eff = 0` anchor is **"no NMS at all", not `r = sqrt(64) = 8 px`.** The two differ — NMS at
   exactly 8 px still collapses tie plateaus (2 of 18,847 peaks on `405.tiff`) — and only "no NMS"
   is the interpretable reference for §2.2 (*what the list looks like with no de-duplication*). It
   is recorded as `nms_radius_px = 0.0` with a `no_nms` flag so the value is defined for `extra`
   and for gate 10's key tuple. It is free of *NMS* cost, not of evaluation cost: it still needs a
   rank, a bucket and a coverage, all already inside the 4–5 s/cell below. `check_no_cap` cannot
   fire on it — `MAX_PEAKS` = 2,000,000 tests exact equality and the largest peak set measured
   anywhere is 48,858.

   giving `r/R₅` ≈ {0.60, 0.75, 0.89, 1.00} — canine 12.2, 15.3, 17.9, 20.2 px; human 13.0, 16.6,
   19.6, 22.1 px.

   *Why `eff` units.* An `R₅`-spaced grid is measurably wrong here: **`0.25 × R₅` is below the 8 px
   floor on every scanner in the dataset** (4.93–5.52 px), so it is a complete no-op that
   reproduces the un-suppressed pool — one wasted level — while leaving `eff` between 0.49 and
   1.00, half the informative range, unsampled. `eff` spacing also puts the grid and arm B on the
   **same dose axis**, which is what lets the grid say where on this curve the `b/51` rule sits.

   *Why `z = 1.0` only is safe despite §10.1.* Every grid level sits on the **same** largest-CC
   template, so all of them share one `med`/`mad`, one deep pool and one `z` cut. The grid varies
   `r` and **nothing else** — it is the only place in the design where the radius is manipulated
   independently of the template, which makes it cleaner than B−C, where `r` is a function of `b`.
   **Rule that follows: never compare a grid level to arm A**, which has a different cut.

   *The null cells are the grid's best cells.* Where `b` = 51 the "largest-CC template" **is**
   Method A's template, so the grid there is a pure radius manipulation on Method A's own search,
   with the template factor removed entirely. The ~29% of cells that are structural zeros for the
   primary are the **cleanest cells in the whole design** for the radius question.

   *Pre-registered shape, so the curve can disconfirm rather than merely describe:* §2.2 predicts
   recall@250 **monotonically decreasing as `eff` decreases**, and more steeply on `chromatin_od`
   than on `tm_score`.

   **8.7.1 Promoted form** (used only if §6.4 row 2 fires). Estimand: mean paired Δrecall@250
   between adjacent `eff` levels, and Spearman ρ of recall@250 against `eff`, over the same 70
   cells. Test and CI: the identical cluster sign-flip machinery of §7.3, including the
   attainable-minimum-p check. Decision rows: §11.1's tree, read against `eff` rather than against
   arm B. Written out here so the branch cannot become a forking path.

   Cost ≈ 4–5 s/cell → **~10 min**, the search being already cached.
8. **Product-side reading metrics**, pre-specified so they are not a post-hoc pick from §8.9:
   `250 − redundant@250` (distinct objects verified per 250 candidates read — the throughput
   quantity AnnotateDx actually pays), and the normalised partial area under recall-vs-K over
   `cp.BUDGETS ≤ 500`.
9. **Full grids** — all budgets × z × axes, robustness only. ~384 cells per contrast; picking
   from it after the fact is the forking path this document exists to close.

---

## 9. Gates

| # | gate | scope |
|---|---|---|
| 1 | **Seed-0 reproduction.** Arm **A** reproduces `results/tm_ccoeff_threshold_axis_sweep_v2.csv`, arm **C** reproduces `results/tm_ccoeff_threshold_axis_sweep_largest_cc.csv`, on `n_detections`/`tp_at_budget`/`lookalike_at_budget` (exact) and 8 further columns (atol 1e-9), keyed on (arm, z, **budget**). Reuses `f1_seed_sweep.gate_seed0`. Verified during review to be passable: `select_domain_images` still returns exactly these 7, and both CSVs carry the fixed `od` padding and cover all 6 z-levels × 8 budgets. | **7 of 14 ROIs**, `si == 0`. Runs before the row is kept. |
| 2 | **Null-cell four-way identity.** `b == 51` ⇒ arms A, B, C, D exactly equal in every column. Requires `r(b)` computed as `R₅ * (b/51.0)`: **measured during review, `R₅ * 51 / 51.0 != R₅` for mpp 0.2268 (`529.tiff`) and 0.2298 (`094.tiff`)**, both in `extra_valid`. `assert r(51) == R₅` at construction. | every null cell |
| 3 | **From-scratch arm-B recompute** — the gate revision 1 lacked. On **one informative cell per ROI**, recompute arm B with a fresh `fused_response`, `od` scored on the post-NMS pool directly, and `coverage_key=None`; assert exact equality with the cached path on all 11 gate columns. ~3 s × 14. **Fallback:** an ROI with no informative cell (possible — `405.tiff` draws 4 nulls in 5, `201.tiff` gave F1 one informative cell in 5) falls back to a **null** cell, which still exercises the arm-B plumbing though more weakly. The **count of ROIs gated on an informative vs a null cell is reported**, so the load-bearing gate's coverage is visible in the output rather than assumed. | 14 cells |
| 4 | `valid.all()` after crop-back | every arm-run |
| 5 | `od` NaN = 0 after `OD_PAD`; `check_no_cap` vs `MAX_PEAKS` | every arm-run |
| 6 | `check_distinct_seeds` | every ROI |
| 7 | One-match-many-`z` shortcut, per (template, radius) | `si == 0`, all arms |
| 8 | `od`-once-per-template equals `od`-after-NMS | `si == 0`, both templates |
| 9 | **§8.7's `eff = 1.00` grid arm equals arm C** at `z = 1.0` on all 11 gate columns. On the largest-CC template `r(eff=1) = R₅` **is** arm C, so recomputing it is either waste or a check; it is made a check. Run with **`coverage_key=None`** so `coverage_frac` is verified from scratch rather than served from a shared cache — otherwise that column is decorative. **70 free assertions**, wider coverage than gate 3. Cost is already inside §8.7's 10 min (the arm is a grid level; the gate is a comparison). | every cell |
| 10 | **Coverage-cache cardinality.** Record the `coverage_key` tuple in `extra` and assert distinct-key counts. **§8.7 gets its own `coverage_cache` dict**, so the two blocks are counted separately with no cross-block arithmetic: factorial cache **24** in an informative cell (4 arms × 6 z) and **6** in a null one (all four arms share one candidate set per z); grid cache **5** (5 `eff` levels at `z` = 1.0 only). Catches the realistic bug — keying on `(template, z)` and dropping the radius — on the A/D half gate 3 does not recompute. | every cell |

**Why gate 3 is the load-bearing one.** GATE 1 covers **A and C only**, on **7 of 14** ROIs. GATE 2
covers B and D **only where they are identical to A by construction**. So without gate 3, arms B
and D — where all three of F4's new optimisations live (the shared search cache, `od`-pre-NMS
subsetting, and the `coverage_key` cache §13's budget depends on) — would be **gated by nothing in
any informative cell**. A miskeyed `coverage_cache` (`compare.py:150-153, 188-193`) survives both
GATE 1 and GATE 2, because at a null cell all four arms are genuinely identical.

**What gates 7 and 8 do *not* cover, restated honestly.** Both are true *independently of the
plumbing they appear to protect*, so neither may be cited as evidence the new code is correct:

* Gate 7's proof is radius-independent — `extract_peaks`' mask is monotone in the threshold, and
  `nms_by_distance` is greedy in descending score with a stable sort (`nms.py:22-49`), so a peak
  is only ever suppressed by a **higher**-scored kept peak. That is exactly why it carries no
  information about the radius change.
* Gate 8 validates the *mathematics* (`chromatin.score_detections` is a pure per-row map,
  `chromatin.py:91-95`), not the *indexing*. **Concrete hazard it would miss:** `pool.iloc[keep]`
  retains the original index while `assign(od=...)` from a list is positional — subset **before**
  assigning, or `reset_index(drop=True)` first.

Both are kept because they are cheap and this repo checks rather than assumes.

**Free external cross-check (directional, not a gate).** §8.7's grid on `300.tiff` and `301.tiff`
can be compared against `results/tm_nms_radius_sweep.csv`'s `pool` and `read_50` at the nearest
radii (20, 16, 12 px) for monotonicity and rough magnitude. The template-selection path differs,
so it is not an exact reproduction — but agreeing in direction with an independently produced CSV
is worth reporting, and disagreeing would be a real signal.

`check_nms_radius` is **deliberately not run** — every arm violates it by design, as v2 does. Arms
carry `nms_radius=None`; `nms_radius_px`, its µm equivalent and `eff(b)` are recorded in `extra`.

### Flags on every cell, retained in the primary, broken out separately

`null_cell` (`b == 51`); `nms_degenerate` (`r(b) < 8`, §3.1 — **not** 7);
`size_below_f1_range` (`b < 19`); `eff(b)`; `n_retries`; `seed_ann_id`; `tightened_size`;
`area_ratio`; `radius_ratio`.

---

## 10. Confounds and limits

1. **`robust_stats` is computed per template** (`template_match.py:66-81`), so `z = 1.0` is a
   **different absolute score cut** in arms {A, D} than in {B, C}. The A-vs-B primary therefore
   changes template **and** cut **and** radius. This is why §8.6 leads the mechanism story with
   `Δn_detections` on **C→B**, which shares `med`/`mad`, and not on A→B.
2. **`recall_at_budget` is not "length-invariant"** — revision 1 asserted this and it is wrong.
   Adding candidates changes *which* candidates occupy the top K, which is the entire effect under
   test. The true statement is F1's: recall@K is not *mechanically monotone* in list length the
   way `read_95` is.
3. **The two dose variables are *perfectly* collinear within arm B** (§3.2) — both are
   deterministic functions of `b`, and `eff(b)`'s only extra variation is scanner-driven and so
   confounded with domain. A vs B alone therefore **cannot** attribute the effect, at all. That is
   why arms C and D and §8.7's grid are pre-registered rather than optional: they are the sole
   source of identifying variation in the design.
4. **`b` is measured on `to_gray_inverted`, not on the matching channel.** §2's mechanism reasons
   about a `hematoxylin_od` template while `b` is an extent measured on a different channel.
   Inherited, and required for GATE 1, but named.
5. **The seed population is conditioned on largest-CC success.** `draw_seed_with_retry` redraws on
   failure, so arm A's marginal seed distribution is not v2's. Harmless for a paired delta.
   Per-cell retry counts recorded.
6. **Padding fabrication differs by arm** near the edge — §8.5 measures it.
7. **`b` is not a controlled dose**; it is whatever the click draws. F1's lesson stands: *"the
   pooled mean is a property of this seed draw, not of the method."* The pooled mean is still the
   deployment estimand — for AnnotateDx the pathologist's click *is* what sets the dose — but is
   reported as an average over the click distribution of these 14 ROIs, never as a context-free
   "B beats A by x".
8. **`n_detections` is not strictly monotone in `r`.** Greedy-NMS keep-sets are not nested:
   with `z > y > x`, `d(z,y)=5`, `d(y,x)=3`, `d(z,x)=8`, `r=7` keeps `{z,x}` while `r=4` keeps
   `{z,y}` — the smaller radius *drops* `x`. Realizable at ≥8 px separations, so §8.6 says "net".
9. **The primary ranking axis is not neutral with respect to the treatment.** `chromatin_od`'s
   justification predates D1's switch to `TM_CCOEFF` and its measured advantage under the current
   configuration is +0.032 recall / 40-of-70 cells, and on `read_50` — the very metric it was
   adopted on — two repo artifacts **disagree on the direction** (§7.2). Two measurements
   disagreeing about which ranker wins is the strongest available argument for carrying both. Worse, §2.2's redundancy mechanism runs *through* the `od` window, so the
   axis most implicated in the treatment is the one pre-registered as primary. Handled by making
   `tm_score` a mandatory companion, not by switching — switching after seeing this would be the
   forking path.
10. **Population.** 14 ROIs, 2 per domain, all `train` split, unanimous clicks only. Not held out.
11. **One method family.** `TM_CCOEFF` on `hematoxylin_od`, single scale, no rotation.
12. **`select_domain_images`' tie is one download from flipping**: `201.tiff` and `233.tiff` both
    have 18 mitoses and `201` wins on the `image_id` tie-break (`experiment.py:32-59`). A future
    canine-lung ROI with ≥19 mitoses silently changes GATE 1's reference set.

---

## 11. Decision rules

### 11.1 Estimation primary (§7.2) — a tree, so the branches cannot overlap

Revision 2 stated these as a flat table whose rows were **not mutually exclusive**: a CI of
[0.005, 0.018] satisfied both "excludes 0, Δ > 0" and "within ±0.02". A decision table with
overlapping rows is not a decision rule. Evaluated in order, on Δ = B − A, recall@250, `z` = 1.0,
`tm_score` (amended; §7.2), using the **ROI-clustered inverted-sign-flip interval** as the
headline (§7.3):

```
0.  Do the two estimands (pooled, informative-cell) agree on whether the CI excludes 0?
    NO  -> SPLIT. Report both, with n_null and mean eff(b). No single verdict.
    YES -> continue.

1.  Does the CI exclude 0?

    YES -> 1a. |mean Delta| >= 0.02 and Delta > 0
               -> RESCUED. Shrinking the radius with the template recovers the
                  largest-CC deficit. Report the gain and the reading cost beside it.
           1b. |mean Delta| >= 0.02 and Delta < 0
               -> WORSE. With F1, largest-CC tightening is rejected in both radius
                  regimes and the matter is closed.
           1c. |mean Delta| < 0.02
               -> REAL BUT SUB-RESOLUTION. A difference that is statistically
                  resolvable and smaller than the threshold below which this design
                  declines to call a winner. Report the sign, the interval, and that
                  it is not a decision-relevant lever at this magnitude.

    NO  -> 2a. CI contained in +/- 0.02
               -> RESOLUTION-BOUNDED NULL. The data rule out a difference larger
                  than 0.02 recall in either direction.
           2b. CI wider than +/- 0.02
               -> UNDERPOWERED. Report the interval, name no winner, state what n
                  would resolve it.
```

**Which estimand steps 1–2 use, fixed now: the pooled one** (§10.7's deployment estimand), with
the informative-cell values reported beside it. Step 0 has already required the two to agree on
the excludes-0 question, which is what the intersection–union rule is for; leaving the 0.02 test
unassigned would be a forking path, because `|pooled| ≤ |informative|` **always** and the two can
land on opposite sides of 0.02.

**Sign disagreement is structurally impossible, so the tree needs no branch for it.** Null cells
contribute exactly 0, so `mean_pooled = (n_informative / 70) × mean_informative` — the two point
estimates are proportional with a positive constant, always share a sign, and satisfy
`|pooled| ≤ |informative|`. Step 0 can therefore only fire on an excludes-0 disagreement, never
on a sign flip.

Branch 1c is a real and likely outcome that revision 2 had no row for.

**On 0.02.** It is a **resolution threshold, not a clinical one.** D4 leaves the clinically
meaningful difference explicitly open pending clinical input, so no equivalence claim is made
beyond "smaller than half the effect F1 measured below 36 px (−0.036)". Every terminal node is
reported next to `n_null`, `n_degenerate`, mean `eff(b)`, `n_nonzero_clusters` and the attainable
`p_min` (§7.3), so no "equivalent" reading can be assembled from cells that were never treated or
from a test that could not have rejected.

### 11.2 D4 decision frame, ruled on separately

Per-domain **worst-of-10** recall@250. **Rule:** Method B is preferred under D4 only if its
worst-of-10 is ≥ Method A's in **every** domain and strictly greater in at least one. Any domain
where B's worst-of-10 falls below A's is a D4 failure regardless of §11.1.

**Reported beside it, because the unpaired minima throw away this design's main strength.**
`min_i recall_B(i)` and `min_i recall_A(i)` are generally attained at *different* clicks, so
comparing them compares two order statistics from different cells, and a min over 10 has no
sampling interval. So also report **`min_i (recall_B(i) − recall_A(i))`** — the worst *paired*
degradation across the 10 cells — and name the cell that attains it. That is the number a
deployment decision actually needs: how bad can this change make any single click.

### 11.3 Decomposition — intervals, no verdicts

Revision 1 pre-committed that "Δ_radius|lcc and Δ_radius|51 disagreeing ⇒ genuinely an
interaction." That was wrong: `(B−C) − (D−A)` is a difference of differences with ~4× the
primary's variance, on ~10 effectively-informative clusters, ~29% of them exact zeros —
**disagreement is the default outcome of noise at this n**. The interaction is therefore reported
as an interval with **no verdict attached**, and §8.1's components likewise.

---

## 12. Outputs

| file | content | when |
|---|---|---|
| `Research Logs/2026-09-04-f4-preregistration.md` | this document | before the run |
| `f4_preflight.py` | §6, forced-`b` on held-out ROIs | before the main run |
| `f4_nms_shrink_sweep.py` | the run; extends `f1_seed_sweep.py` (`run_arm` already takes `nms_radius`) with a cached search stage | before the run |
| `f4_analysis.py` | every test in §7–§8, committed **before results exist** (F1's convention), including §7.1's assert | before the run |
| `results/f4_preflight.csv` | §6.3's four measurements incl. `gap(b)` | Step 0 |
| `results/f4_nms_shrink_sweep.csv` | row per ROI × seed × arm × z × axis × budget, `.partial` until complete | during |
| `results/f4_radius_grid.csv` | §8.7's grid, kept **out** of the factorial file so §7.1's assert stays meaningful | during |
| `results/f4_nms_shrink_cells.csv` | row per ROI × seed: seed id, retries, `b`, radii, `eff(b)`, `coverage_key` cardinality, gate-3 informative/null flag, flags, per-arm pool sizes | during |
| `Research Logs/2026-09-04-f4-results.md` | findings against §11 | after |

## 13. Compute budget

Interpreter `/Users/mohinianand/anaconda3/bin/python` (3.11.5, numpy 1.26.4, cv2 4.8.1) — the
numpy < 2 pin `requirements-midog-utils.txt` documents; the repo default `python` is **Python 2**
and system `python3` cannot `import cv2`.

Anchored on F1's **measured 26 s/cell for 2 arm-runs** (912 s / 35 cells), not on the per-stage
probe alone — revision 1 projected 16 s/cell for *four* arm-runs, which was inconsistent with F1's
own measurement. F4 adds 2 further evaluations per cell (~10 s each), less the `coverage_key`
saving of ~1.7 s per arm-run from sharing coverage across the two axes at each
(template, radius, z):

* main run ≈ **40–45 s/cell** → 70 cells ≈ **45–55 min**
* §8.7 radius grid — 5 levels × 2 axes, each needing `_rank` over ~30k rows plus a
  `bucket_detections` KDTree, and 5 coverage computations ≈ **4–5 s/cell** → **~10 min**
  (revision 2's "2 s/cell → 7 min" undercounted the ranking and bucketing)
* §6 pre-flight ≈ 6 doses × 4 ROIs × 15 s ≈ **6 min**
* gates 3, 9, 10 ≈ 1–2 min

**Total ≈ 65–75 min.** Memory: `fused` is ~130 MB per template at 4933×6577 — free it after
`robust_stats`/`extract_peaks` or peak RSS doubles under the sharing scheme.

Compute is not the constraint, which is why arms C and D and §8.7 are affordable and why cutting
them would be a false economy.

---

## 14. Changes from revisions 1, 2 and 3

Revision 3 was reviewed a third time. Three further blocking defects, all found by checking the
plan against the **repository** rather than against the code and F1:

1. **`results/tm_nms_radius_sweep.csv` was uncited** — a pure NMS-radius sweep at a fixed 51 px
   template on `300.tiff` and `301.tiff`, **both in `extra_valid`**, 5 seeds × 5 radii, dated
   2026-09-02. §2 exists to stop this run re-answering measured questions, and it was about to
   re-answer this one. → new §2.0, with the three facts re-derived: pool rises in 10/10 streams
   and 40/40 steps; `read_50` worsens in 9/10; and **the tail gain saturates at 24 px, above
   `R₅`** — which is §2.1's conclusion turned from an argument into a measurement.
2. **F5 is a sibling pre-registration on the same lever, running concurrently, and the seed draws
   coincide** — verified: `013.tiff` `si` = 0 gives ann 254 under both rules. Neither document
   mentioned the other. → the "Overlap with F5" pre-commitment in §5, and the open question of
   which document owns the fixed-template radius factor.
3. **§7.2's `read_50` ratio was stated in the wrong orientation**, making "the advantage has not
   reversed" false: `read_50(od)/read_50(tm)` is **median 1.021**, so `od` reads *deeper*, and is
   better in only 33 of 70 cells. The 19.4× reduction has become a 1.02× increase. → §7.2, §10.9
   and this section corrected. The design consequence (`tm_score` as mandatory companion) was
   right and is now properly argued.

Also in this pass: the head-vs-tail trade named in §2.2 with §8.8's partial AUC leaning on it;
§11.1's tree now says steps 1–2 use the **pooled** estimand and notes that sign disagreement is
structurally impossible (`mean_pooled = (n_inf/70) × mean_inf`); §8.7 gets its **own**
`coverage_cache` and gate 9 runs with `coverage_key=None` so its coverage column is a real check;
the `eff = 0` anchor defined unambiguously as **no NMS at all** (`nms_radius_px = 0.0`), not
`r = 8 px`; §6.4's 0.02 threshold derived honestly as equivalent to `setdiff = 0` given probe
`n_gt` ≤ 15, with raw `setdiff` reported alongside; "13 arm columns" corrected to 9; and the
directional cross-check of §8.7's grid against `tm_nms_radius_sweep.csv` added to §9.

## Changes from revisions 1 and 2

Revision 2 was reviewed a second time. Four further blocking defects, and their fixes:

1. **§6's motivating claim was false as measured.** It said a drawn pre-flight would "very
   likely" median to zero; measuring all four probe ROIs shows **1 of 8 cells null, not most**.
   The claim was over-generalised from `405.tiff`, the only ROI either party had measured. The
   real and stronger defect is **dose coverage** — the drawn doses span `b ∈ [29, 45]` and never
   reach the `b ≤ 27` regime the experiment is about. → §6 opening and §6.6.
2. **§7.1's assert collided with §8.7's rows**, which went to the same file: the pivot would have
   returned 13 arm columns, not 4. → separate `results/f4_radius_grid.csv`, and the filter names
   the four arms explicitly.
3. **§6.4's two branches did not partition the outcome space** — different statistics in the two
   rows meant both could be false, leaving the analyst to choose. → one per-dose statistic in
   §11's own units, "any" vs "all".
4. **The domain-clustered test can be arithmetically unable to reject.** Support is
   `2^(n_nonzero clusters)`; F1 had **6 of 7** (`201.tiff` identically zero), giving `p_min` =
   0.031 with no margin, and 5 non-zero gives 0.0625 > α — while §11.1 conditioned every row on
   that level. The risk is identifiable: canine lung pairs the only two sparse ROIs in the set.
   → the attainable-minimum-p rule in §7.3.

Also in this pass: the `eff`-spaced radius grid, replacing an `R₅`-spaced one whose `0.25 × R₅`
level was **below the 8 px floor on every scanner** (§8.7); §8.7.1, so the promoted form is not a
forking path; the collinearity statement **restored** to revision 1's, which was right and which
revision 2 replaced with a wrong one (§3.2); the rank channel restated as **signed positive** and
made measurable via `matched_rank` (§2.1, §8.6); §11.1 turned into a tree with a branch for the
real-but-sub-resolution outcome its flat table had no row for; §11.2's paired worst-case;
§8.2 promoted to a named deliverable; gate 3's fallback and coverage reporting; new gates 9 and
10; the intersection–union justification replacing a wrong multiplicity claim (§7.4.1); `Δ_count`
labelled GT-weighted with a per-domain companion; `gap(b)` redefined over suppressed peaks only;
and the compute budget corrected to 65–75 min.

### Changes from revision 1 (retained for the record)

Revision 1 was reviewed adversarially and found **not ready to run**. The four blocking defects
and what fixed them:

1. **The primary frame's dedup key was never stated** and F4 adds a fourth key level — D4's own
   documented 8× trap. → §7.1, with an assert.
2. **Arms B and D were gated by nothing in any informative cell**, and that is where all three new
   optimisations live. → gate 3.
3. **§6's pre-flight median would have been contaminated by null cells** (`405.tiff` draws `b`=51
   in 4 of 5 seeds), `topk_churn` was undefined, and the failing branch had no consequence. → §6
   rebuilt as a forced-dose probe with a real branch.
4. **"External replication" was not external** — `extra_valid` contains all 7 of F1's ROIs and F4
   redraws F1's exact 35 cells on them. → §5, §8.2, §8.3 restricted to the 7 new ROIs.

Added in the same pass, from a question raised after the review: the **ranking axis** is now
treated as part of the treatment rather than a neutral reporting choice (§2.2, §7.2, §10.9) —
`chromatin_od`'s justification predates D1's matcher switch; under the current configuration its
advantage is +0.032 recall on recall@250 and is **contested on `read_50`** — the metric the 19.4×
was measured on — by two repo artifacts that disagree in direction; and §2.2's redundancy
mechanism runs through the `od` window itself. `tm_score` is now a mandatory companion axis.

Also corrected: the shrink rule's false provenance (§3); the degeneracy threshold, 8 px not 7,
and it bites inside prior experience (§3.1); `eff(b)` as the real dose (§3.2); the `od`-window
redundancy mechanism (§2.2); the directional prior (§2.1–2.3); null-cell dilution of an
equivalence claim (§7.4); recall@250's 1/17 granularity on four ROIs (§7.2, §7.4.2); the D4
conflict (§7.5, §11.2); float-exactness of `r(51)` (gate 2); `robust_stats` per template (§10.1);
"length-invariant" (§10.2); non-nested NMS keep-sets (§10.8); and the compute budget (§13).
