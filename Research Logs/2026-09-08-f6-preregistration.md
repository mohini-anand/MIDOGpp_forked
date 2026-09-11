# Pre-registration F6: what does border padding actually do, at 14 ROIs × 5 seeds — and is it bigger than the click?

Date: 2026-09-08. **Written before `f6_padding_ablation.py` is run.** Status: PLAN — revision 4,
CONVERGED after three rounds of adversarial review. Not yet executed.

F1–F5 are taken. This is the next free number.

> **Revision 2, same day.** Draft 1 decomposed padding into *reachability* and *re-normalisation*.
> Both channels are now measured at approximately **zero** before the run — the reachability
> ceiling is 0 of 7 annotations and the re-normalisation shift is ≤ 2.27e-3 MAD, which moves
> nothing at K = 250 — while the mechanism that demonstrably produces the observed losses was
> neither named nor isolated. The arms, the predictions and §2's mechanism paragraph are all
> replaced. Draft 1's error is recorded in §13 rather than quietly overwritten.
>
> **Revision 3, same day.** Review then *built* revision 2's arms on the two cells where the total
> is non-zero, and falsified two of its four predictions. Suppression contributes **exactly 0 %**
> at K = 250 and the band candidates' own occupancy **100 %**, with the sign opposite to the one
> registered. The mechanism names were right; their magnitudes were inverted. Re-signed in §7,
> with what revision 2 would have claimed recorded in §13.
>
> **Revision 4, same day — converged.** Two last items. The tie-aware floor rule is applied
> *prospectively* to F6's own `tm_score` test, where it shows that test is likely unable to reject
> at all (§9). And P4's mechanism is replaced by one that is quantitatively predictive rather than
> merely directional (§7).

---

## 1. The two questions

**Q1 — what does border padding do?** With the NMS radius held at its default (7.5 µm = the
evaluation match radius), what does turning padding **on** do to retrieval, and *through which
mechanism*?

**Q2 — the noise floor it must clear.** Using recall@K at a fixed budget, how much does *which cell
the pathologist clicks* move performance, per ROI and per tumour domain?

Q2 is the ruler for Q1. `results/tm_ccoeff_headtohead_seed_variance.csv` records a seed-to-seed SD
of 260–2744 candidates on `read_95`, larger than every per-domain depth delta F5 reported. An
effect smaller than the click is not a result.

---

## 2. What is already measured, and what it rules out

The padding question was answered once on 2026-09-08 against
`results/tm_ccoeff_threshold_axis_sweep.csv` (v1) — **7 ROIs, `SEED_INDEX = 0`**. Everything below
has been reproduced independently at least twice.

### 2a. The headline numbers

| | |
|---|---|
| Evaluation mitoses inside the 25 px margin | **7 of 902** (0.78 %) |
| Δ recall@250, `chromatin_od` | **−0.0071** — 0 of 7 ROIs better, 4 worse; p = 0.125, which **is that test's attainable floor** (4 non-tied clusters ⇒ 2/2⁴) |
| Δ recall@250, `tm_score` | **+0.0187** — 90 % of it is 201.tiff's 2 annotations on a 17-mitosis ROI; **excluding 201 it is +0.0022** (mean over the remaining 6) |
| Pooled Δ TP@250, both axes | **−3** |
| Δ `n_detections` | **+248 (+1.75 %), positive on 7 of 7, p = 0.0156** — floor-limited, and the *smallest* p in the comparison |
| Δ `read_50`, `chromatin_od` | **+3.7**, p = 0.031 — also floor-limited (6 non-tied ⇒ 2/2⁶) |
| Δ `coverage_frac` | +0.0082 |

Both p-values sit exactly on their own floors, so each means *"every moving cluster agreed"*, not
*"significant"*. §9's tie-aware rule is applied to this table, not only to F6's own output.

### 2b. Two mechanisms are dead before the run

**Reachability has a ceiling of zero.** For each of the 7 border mitoses, the distance from its
click to the nearest pixel of the unpadded valid rectangle, against that ROI's match radius:

```
201/4494 12.0px   201/4502  8.0px   246/6612  9.0px   246/6618 10.0px
301/14779 16.0px  301/15022 3.0px   402/20275 8.0px      radii 29.61 - 33.14 px
```

**0 of 7 are further than one match radius from the valid region.** The unpadded search starts
25 px in but the match radius is ~30 px, so a detection on the first valid row can still claim an
annotation up to 30 px outside it. Nothing on this data is *unreachable*; the 2 annotations padding
recovered (301, 246) were already geometrically claimable and were recovered by a better-placed
peak — claim churn, not new reach. **v2's stated justification for padding has no headroom here.**

**Re-normalisation is inert at the operating point.** Draft 1 read "`map_median` moves up to 1.77 %
relative" as a global threshold shift. That is a near-zero-denominator artefact: the absolute move
is 0.0001 against a MAD of 0.194. Expressed in the unit that matters — MAD, the unit `z` is
measured in — the z = 1.0 cut moves by at most **2.27e-3 MAD** across the 7 ROIs, on a grid whose
step is 0.5. And the top-250 sits far above any cut in the grid: the minimum `n_detections` at
z = 3.0 across all 140 F5 `r7.5` cells is **2,301**, so the 250th-ranked candidate is above z = 3.0.
A 0.002-MAD shift relocates nothing there. **Re-normalisation is therefore folded into the noise
and given no arm.**

### 2c. The mechanism that is left, and that the data actually supports

Replicate-padded border pixels produce **real-valued peaks**. Those peaks do two things that have
nothing to do with reaching a border annotation:

* **they suppress interior local maxima** — via `extract_peaks`' grey dilation within
  `peak_min_distance = 7` px, and via NMS within ~30 px. Measured on single cells: **99–122
  interior peaks present unpadded and absent padded**, ~98 % of them within 30 px of the band.
* **they consume top-K reading slots** — measured: **5–12 of the top 250** are border-band
  detections, i.e. 2–5 % of the reading budget spent on replicated pixels.

The decisive evidence that this, and not the two dead channels, is what moves the numbers:
**459.tiff has zero border mitoses, no measurable re-normalisation, and still loses a TP at
K = 250** on `chromatin_od`. Draft 1's syllogism was "not reachability ⟹ re-normalisation"; a third
mechanism sits in that gap and it is the one F6 must isolate.

**And it is axis-specific, which is itself a clue.** Per-ROI ΔTP@250: 402 is −3 on `chromatin_od`
and **0** on `tm_score`; 459 is −1 and **0**. An independent ranking key promotes border artefacts
into the top 250; a score-ranked list, sorted by the same quantity NMS sorts by, does not.

---

## 3. Design

Per `(ROI, seed)`, build the template **once** from the seed patch — a property of the seed, not of
the arm — then run the search two ways and split one of them:

**Two response maps, and the binding reason is gate A — not the dilation.** The padded map's
interior differs from the unpadded map's by up to **3.8e-6** (OpenCV's DFT block tiling depends on
image size), and gate A compares `n_detections` on *exact integer equality*. Measured: a
`_FLOOR`-written emulation still loses 1 peak of 20,855 on 459.tiff and 6 of 31,104 on 402.tiff —
enough to flip that comparison. An emulated `pad_off` would not be v1, and the gate that makes
`pad_off` trustworthy would be testing the emulation instead of the pipeline.

The grey-dilation hazard is separately real but is **not** the binding reason, and revision 2 was
wrong to rest on it: `extract_peaks` computes `dilated = cv2.dilate(fused, ...)` before applying
`& valid`, so a naive `valid`-mask emulation loses **228–267 peaks of 21k–31k** — but writing
`np.where(valid, fused, _FLOOR)` first, the pattern `template_match.plant_and_recover` already
uses, defeats it entirely. Recorded so a future reader does not delete the second `matchTemplate`
on the strength of an argument that does not hold.

### The three arms

| arm | how it is produced | reachability | suppression by border peaks | border candidates occupy top-K |
|---|---|---|---|---|
| **`pad_off`** | unpadded response; `valid` = interior only | — | — | — |
| **`pad_on_dropband`** | padded response, full extraction + NMS, **then border-band candidates dropped from the ranked list** | — | **yes** | — |
| **`pad_on`** | padded response, cropped back, `valid.all()` | yes | yes | yes |

The "border band" is `cx < PAD or cx ≥ W−PAD or cy < PAD or cy ≥ H−PAD` with `PAD = (t−1)//2 = 25`
— exactly the region the unpadded search could not place a detection in.

Read across, the two contrasts isolate the two live mechanisms:

* **`pad_on_dropband` − `pad_off` = suppression, plus the inert re-normalisation channel.** Border
  peaks competed during dilation and NMS, but none is in the list. Calling this "suppression alone"
  would repeat draft 1's error one arm to the left: `pad_on_dropband` inherits `pad_on`'s
  `med`/`mad`/`deep_floor`/`z_cut` while `pad_off` uses the unpadded map's. §2b gives that channel
  no arm because it is bounded at 2.27e-3 MAD and is zero at K = 250; §5's per-row statistics make
  that checkable on F6's own data.
* **`pad_on` − `pad_on_dropband` = the border candidates' own contribution** — the TPs they claim
  minus the reading slots they consume. **Pre-measured on two cells, this is the whole effect**
  (§7 P2/P3).
* **`pad_on` − `pad_off` = the total**, and the only contrast between two configurations that could
  ship. **This is the primary contrast.**

`pad_on_dropband` needs **no extra response map, extraction, NMS or `od` scoring** — it is one
boolean mask on `pad_on`'s ranked frame. Draft 1's third arm cost all four and measured a channel
pinned at zero.

**Held constant in every arm:** NMS radius = scoring radius = 7.5 µm (so
`invariants.check_nms_radius` passes on all three — F5's decoupled radius is *not* used);
`hematoxylin_od`; `TM_CCOEFF`; single 51 px template; `peak_min_distance = 7`;
`self_hit_radius = 5.0` after NMS; `DEEP_FLOOR_Z = −1.5` pinned; the seed draw; the ground truth;
and — a change from draft 1 — **`chromatin.score_detections` on the `OD_PAD = 25` padded channel in
all three arms.** Draft 1 scored `pad_off` on the unpadded channel for fidelity to v1; that was
verified to be a no-op (every `pad_off` candidate is ≥ 25 px from the edge, so its 51 px window
covers identical pixels either way, and the padded-scored arm still reproduced v1's
`chromatin_od` recall@250 to six decimals). Using one channel everywhere removes an asymmetry that
costs nothing to remove.

---

## 4. Sample

* **14 ROIs** — `images/extra_valid`, 2 per tumour type across all 7 domains.
* **5 seeds per ROI** — `tm_variant_sweep.draw_seeds`: agreement pool → `border_filter(36)` → drawn
  **without replacement**, one `default_rng([s, image_id])` stream per index.
* **70 paired cells in 14 ROI clusters.** The ROI is the design unit.
* **The seed pool size (7–184) is reported beside every Q2 number.** 201.tiff's 5 seeds are 71 % of
  its eligible population; 301.tiff's are 3 %. Those estimates are not on the same sampling footing
  and the table must not pretend otherwise.

**Power, stated before the run rather than after.** From the 7-ROI per-ROI deltas that already
exist, the between-ROI SD of Δrecall@250 on `chromatin_od` is ≈ 0.011, so at 14 clusters
SE ≈ 0.0029 and **MDE₈₀ ≈ 0.008** against an expected effect of −0.0071. **F6 is marginally
powered for the total contrast and under-powered for the decomposition.** This is registered so
that a null is read as "under-powered", not as "no effect".

---

## 5. Nuisance axes

* **`z ∈ (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)`, headline z = 1.0** — F5's grid, so cells line up.
* **`DEEP_FLOOR_Z = −1.5`** pinned as a constant, not derived from `min(Z_LEVELS)` — this is what
  makes §8's reproduction gates possible.
* **Both ranking axes co-primary, neither secondary** (F5 Correction 2 §2.6) — plus §9's registered
  within-F6 head-to-head, so this does not rest on an imported post-hoc result.
* **Arm naming `{axis}@{arm_tag}`**; dedup key `(file_name, seed_index, arm, z)`.
* **`coverage_key = (file, seed, arm_tag, z)`.**
* `map_median`, `mad_scale`, `deep_floor` and `z_cut` recorded on every row, so §2b's claim that
  the shift is inert is checkable on F6's own data rather than inherited.

---

## 6. Metrics

### 6.1 Q1 — the padding effect and its mechanism

**Primary: recall@250** (D4), z = 1.0, both axes, over

```
BUDGETS_F6 = (10, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500,
              750, 1000, 1500, 2000, 3000, 5000, 10000)
```

a strict superset of `compare.BUDGETS`, so rows line up with v1 and F5.

**Mechanism instrumentation — the part draft 1 omitted, and the only part that can explain a signed
Δ now that both of its named channels are dead:**

* **found-annotation ledger** per (ROI, seed, arm, axis): the set of mitotic `ann_id`s found,
  differenced against `pad_off` → **gained / lost**, and for each gained or lost annotation its
  distance to the ROI edge, so a change can be attributed to the band or to the interior.
* **border-band detections in the top-K**, per arm and axis, at K ∈ {50, 250, 1000} — the direct
  product cost of padding, and the quantity §2c measures at 5–12 of 250.
* **interior peaks lost, split by mechanism** — two different causes with two different fixes, and
  lumping them hides both: **pre-NMS** losses to `extract_peaks`' grey dilation within
  `peak_min_distance = 7` px (measured 267 on 402.tiff/s0) and **post-NMS + z** losses (119 on the
  same cell). "Identity of the suppressing peak" is only defined for the NMS case.
* **Band *detections* and band *annotations* counted separately, never conflated.** The band is
  25 px and the match radius ~30 px, so a detection at `cx = 26` can claim an annotation at
  `cx = 2`, and a band detection at `cx = 20` can claim an interior annotation at `cx = 45`.
  §6.1's top-K counts are about detections; P5 is about annotations.
* **border-band candidate count** in the full list, per arm.

Also carried: `full_list_recall` **printed only beside `coverage_frac`**; `n_detections`; FP@K; the
full-list FP count; `FP_to_r` for r ∈ {0.8, 0.9} (0.95/0.99/1.0 reported but coverage-guarded);
`n_dup_fp` with its denominator `n_within_match_radius`; `read_{50,80,90,95}`; `largest_tie_block`;
`nan_rate`.

### 6.2 Q2 — seed variance at a fixed budget

**Computed from `results/f5_nms_radius_ablation.csv` BEFORE the run**, for the `r7.5` arm, which
§8's gate B asserts is byte-identical to F6's `pad_on`. Q2 for the shippable configuration needs no
new compute; F6 adds the same table for `pad_off` and `pad_on_dropband`, which is new.

The statistic, chosen deliberately:

* **Click noise = the mean of the two within-ROI seed SDs of recall@250 in a domain.** This is what
  "how much does the click move performance" means.
* **The pooled-10 SD (2 ROIs × 5 seeds) is reported *separately and labelled as ROI
  heterogeneity*, never as click noise.** Draft 1 used it as the noise floor; on the existing data
  it inflates the within-ROI figure by **10× on canine lymphosarcoma and 41× on human breast
  cancer**, because the 14-ROI draw is unbalanced within domains (melanoma's two ROIs differ 11.95×
  in mitosis count). It measures between-ROI signal, not the click.
* **Zero-variance ROIs are reported as zero, never as a denominator.** At least 013.tiff and
  233.tiff have within-ROI SD exactly 0 at K = 250 on `chromatin_od` (recall pinned at 1.0 on all
  5 seeds). Any ratio with such a denominator is reported as `n/a`, not as ∞ or a large number.
* **A granularity baseline accompanies every SD**: one TP is `1/n_eval` of recall, which runs
  **0.0042 (238 mitoses) to 0.059 (17 mitoses)** across the 14 ROIs — a 14× difference in the
  smallest SD that is even expressible. An SD at that scale is quantisation, not click sensitivity.
* **No CV.** It is meaningless for a bounded metric near its ceiling (013 has mean 1.0, SD 0 → CV 0,
  which does not mean "less variable than" an ROI at 0.30). SD is reported with its mean.

Secondary, for shape: the same spread at K ∈ {50, 100, 500, 1000}.

### 6.3 Putting Q1 next to Q2 — two rulers, kept apart

Draft 1 divided a *paired* effect by an *unpaired* SD. Those are not on the same footing: the
numerator's variance is tiny because both arms share the tissue, the ground truth and (to ~1e-6)
the response map, while the denominator is the spread of the *levels*. On the existing numbers that
rule fires "smaller than the click" wherever the click matters and "bigger than the click" wherever
recall is pinned at 1.0 — backwards. So:

* **Statistical ruler (same footing):** the **seed-to-seed SD of the paired difference** within each
  ROI. This says whether the padding effect is stable across clicks, and it is the quantity §9's
  cluster test already implies.
* **Practical ruler (labelled as such):** |mean paired Δ| against the within-ROI SD of the *level*,
  reported per domain with `n/a` where the denominator is 0. It is a **relevance heuristic and
  cannot override or be overridden by** §9's arbiter; §10's decision rule is bound to the arbiter
  alone.

---

## 7. Registered predictions

Each carries a number, so each can fail.

* **P1 — the total.** `pad_on − pad_off` at z = 1.0, K = 250: negative on `chromatin_od` with a
  cluster-bootstrap CI whose upper bound is below +0.005, and within ±0.005 of zero on `tm_score`.
* **P2 — suppression is inert at the operating point.** `pad_on_dropband − pad_off` is within
  ±0.002 of zero on both axes and accounts for **less than 20 %** of the total on `chromatin_od`,
  *despite* the interior peaks it demonstrably destroys — those losses land in the tail, on
  border-adjacent low-scoring peaks that were never in the top 250. **Revision 2 registered the
  opposite** (≥ 50 % of the total, negative on both axes); building the arms measured **exactly
  0.00000, 0 % of the total**, on 2 of 2 non-zero cells and both axes.
* **P3 — the band candidates are the whole cost, and it is negative.** `pad_on − pad_on_dropband`
  is **negative** on `chromatin_od`, magnitude up to 0.03. Revision 2 registered it as *positive
  and below 0.01*; 402.tiff/s0 measured **−0.02885** — three times that bound, and the other way.
  Band candidates sit on replicated pixels, so they are near-certain false positives: they consume
  reading slots and claim almost nothing (2 annotations recovered across the entire prior 7-ROI
  comparison).
* **P4 — the effect is `chromatin_od`-only, and the mechanism is a product of two measurable
  factors.** Total, suppression and band contributions are all within ±0.002 of zero on `tm_score`
  **even though band detections do sit in its top 250** — 5 of 250 on both probe cells, with
  total = 0. So the cost is *not* "band detections enter the reading budget", and an earlier
  formulation saying a score-ranked list keeps them out is simply false: they enter.

  Adding `N` band detections into the top-K pushes the **last `N` entries out**, whatever the key.
  The TP cost is therefore

  > **(band detections entering the top-K) × (TP density at the list margin)**

  and `chromatin_od` loses on *both* factors. Measured in `pad_off` over ranks 225–249:

  | cell | axis | band dets entering | TPs in ranks 225–249 | predicted loss | observed ΔTP |
  |---|---|---:|---:|---:|---:|
  | 402/s0 | `chromatin_od` | 12 | 5 of 25 | 2.4 | **−3** |
  | 402/s0 | `tm_score` | 5 | 1 of 25 | 0.2 | **0** |
  | 459/s0 | `chromatin_od` | 8 | 2 of 25 | 0.6 | **−1** |
  | 459/s0 | `tm_score` | 5 | 1 of 25 | 0.2 | **0** |

  The product predicts the observed loss on 4 of 4. **Registered as a quantitative prediction:**
  across the 70 cells, the correlation between `(band_in_topK × margin_TP_density)` and the
  observed ΔTP@250 is positive with a cluster-bootstrap CI excluding 0. F5 Correction 2 §2.7 gives
  the base rate for the axis asymmetry: 9 of 14 ROIs exactly zero on `tm_score`.

  Instrumented by two ledger columns F6 already nearly has: **band detections in the top-K** (§6.1)
  and **`margin_tp_density`** — the TP count in ranks `K−25 … K−1` of the `pad_off` ranking, per
  (cell, axis, K).
* **P5 — reachability contributes essentially nothing.** Of the annotations gained by `pad_on` over
  `pad_off`, **fewer than half** lie in the band. **Reported untestable if fewer than 5 annotations
  are gained in total** — §2b's ceiling of 0 unreachable annotations plus the prior recovery of
  exactly 2 makes a small count the expected case, and "fewer than half of 1" is a coin flip.

**Withdrawn from draft 1:** its P2 (re-normalisation dominates) — the channel is measured at zero,
so the claim was arithmetically unreachable; and its P4 (seed SD varies ≥3× and the largest-SD
domain is a low-count one) — checked against the existing CSV, the largest-SD domain is canine
**lymphosarcoma**, the third-densest, on 3 of 4 readings, and a ≥3× spread arises from pure
sampling noise 42 % of the time under the within-ROI statistic. Draft 1 also mis-stated 013/529 as
a domain; they are not (013 is breast, paired with 094; 529 is melanoma, paired with 548).

**What would make F6 report "padding is a real gain"** — bound to §9's arbiter, not to §6.3's
heuristic: the total contrast positive on both axes with cluster-bootstrap CIs excluding 0, and
`pad_on − pad_on_dropband` accounting for the majority of it. Anything less is reported as
"padding closes a leak with zero headroom, at a measurable reading cost."

**Registered expectation: `pad_on_dropband`'s recall@250 column will be zeros.** On 4 of 4
(cell × axis) combinations pre-measured it reproduced `pad_off` exactly. That is the *expected*
result and is what proves the attribution — not a bug and not a broken join (gate 12 covers that).
The arm's informative outputs are the ledger, the full-list level, and the finding that **dropping
the band from the ranked list recovers 100 % of the loss**, which if it holds at 14 × 5 is a
one-line shippable fix rather than an argument about whether to pad.

---

## 8. Verifications and gates

1. `valid.all()` on the padded arms; on `pad_off`, `valid.sum() == (H−2·PAD)·(W−2·PAD)` exactly —
   provable from `PAD = (t−1)//2`, not empirical.
2. `od_nan == 0` in all three arms.
3. `invariants.check_nms_radius` passes on **all three** arms (F6 does not decouple the radius).
4. `invariants.check_min_separation` on every arm's kept centres at 7.5 µm.
5. `invariants.check_no_cap` on the **pre-NMS** pool against `MAX_PEAKS`.
6. `invariants.check_distinct_seeds` on the 5 RNG streams per ROI.
7. **Interior identity**, checked once per ROI at `si == 0`: the padded map cropped to the interior
   equals the unpadded map on the interior to `atol = 1e-4`. Measured max difference **3.8e-6**
   (OpenCV's DFT block tiling differs with image size), so the two arms are *identical to ~1e-6* —
   **not** "byte-identical", which draft 1 claimed and which would have made gate A fail.
8. **Join integrity, not "additivity".** `(A−C) == (A−B) + (B−C)` holds for any three reals; it
   cannot fail except on NaN or a broken join, and draft 1 presented it as evidence that the
   decomposition holds together. Kept, demoted to what it is.
9. **Two reproduction gates**, both free, both whole-pipeline:
   * **A.** `pad_off` at `si == 0` reproduces `results/tm_ccoeff_threshold_axis_sweep.csv` on the
     **7 overlapping ROIs** at shared `(z, budget)`, both axes.
   * **B.** `pad_on` at `si == 0` reproduces `results/f5_nms_radius_ablation.csv`'s `r7.5` arm on
     **all 14 ROIs** at shared `(z, budget)`, both axes.
   Both were pre-verified on 2 ROIs during review and reproduce exactly. **Scope, carried over from
   F5's implementation:** the gate compares `recall_at_budget` (allclose) and `n_detections` (exact
   integer); `full_list_recall` is merged but not compared. Mapping: v1's `arm` is the bare axis,
   F5's and F6's are `{axis}@{tag}`.
10. The one-match-many-`z` shortcut is re-verified at z = 1.0 **and z = 3.0**, on `si == 0` only,
    for both extracted arms.
11. `n_gt_within_seed_hole == 0` per cell.
12. **`pad_on_dropband`'s candidate set is a strict subset of `pad_on`'s** — asserted on the
    candidate sets, not the ranked frames, since ranks shift by construction when rows are removed.

---

## 9. Statistics

> **Arbiter: the mean-of-5 paired Δ(recall@250) per ROI** at z = 1.0 for `pad_on − pad_off`,
> computed **separately on each axis** — 14 clustered units per axis.

* **Test:** exact cluster sign-flip permutation over the 14 ROIs (2¹⁴ = 16,384, enumerated), plus a
  14-ROI cluster bootstrap percentile CI (`seed = 20260908`, 10,000 draws).
* **The tie-aware attainable floor is printed next to every p.** With `k` non-zero clusters the
  smallest reachable two-sided p is **2/2ᵏ**, not 2/2¹⁴. A p equal to its floor means "every moving
  cluster agreed", not "significant"; a floor above 0.05 means the test could not have rejected at
  all. §2a applies this rule to the prior numbers too.
* **Multiplicity: Holm across the 2 axes** — one formal contrast per axis.
* **The `tm_score` arm of the arbiter is expected to be floor-limited, and this is registered
  before the run rather than discovered after it.** The floor rule above is not only retrospective:
  applied *prospectively* to F6's own primary test, it says that with `k` non-zero clusters the
  test can only reject at Holm's α/2 = 0.025 if **2/2ᵏ ≤ 0.025, i.e. k ≥ 7**. F5's base rate on
  `tm_score` is 9 of 14 ROIs exactly zero — `k = 5`, floor `2/2⁵ = 0.0625`, **above 0.05, so the
  test could not reject at any α it would be judged at**. Padding's zero rate looks *higher* still:
  exactly zero on 3 of 3 ROIs built at seed 0. If it lands at 12 of 14, `k = 2` and the floor is
  0.5.
  **So: if `k < 7` on `tm_score`, its formal test is reported as UNINFORMATIVE, not as a null.**
  The descriptive cluster-bootstrap CI carries that axis instead, and the Holm-adjusted
  `chromatin_od` p is printed beside its unadjusted value so the α spent on an unrejectable test is
  visible. Without this, P1's "within ±0.005 of zero on `tm_score`" would be confirmed by a test
  that could not have failed — the same error §2a exists to correct, one document later and on
  F6's own primary test.
* **The two decomposition contrasts get cluster-bootstrap CIs but no p and no Holm family.** They
  are outside the multiplicity family, and giving them intervals is what makes P2 falsifiable —
  draft 1 called one of them "the live claim" while refusing it an interval.
* **D5 reconciliation, registered rather than imported:** the `chromatin_od − tm_score` head-to-head
  on `pad_on` at z = 1.0, K = 250, is computed **within F6** as a separately named descriptive
  contrast with a cluster-bootstrap CI. `DECISIONS.md` D5 requires that any experiment treating
  `chromatin_od` as primary measure it beating `tm_score` **on that run's own data**; F5's
  head-to-head was run after seeing both axes' answers and is not carried here as authority.
  **Contingency, registered now because the entire padding effect is `chromatin_od`-only and this
  must not be decided once the answers are visible: if that CI includes 0, `chromatin_od` is
  reported as secondary and the arbiter is `tm_score` alone** — under which, on the pre-measured
  cells, the total padding effect is zero.
* **Dedup key `(file_name, seed_index, arm, z)`** before any count, correlation or test.
* **Reporting frame (D4):** worst-of-10 recall@250 per domain, per arm, no p-value — a reporting
  rule, a different statistic from the arbiter.
* **Q2 is reported before Q1's test**, so the noise floor is on the page before any effect is
  claimed against it.

---

## 10. Deliverables

1. `f6_padding_ablation.py` → `results/f6_padding_ablation.csv`,
   `..._verification.csv`, `..._ledger.csv` (§6.1's mechanism instrumentation),
   `..._seedvar.csv` (§6.2), `..._ceiling.csv` (per-(ROI, seed) border-annotation counts and
   distance-to-interior).
2. `f6_padding_ablation.ipynb` — reads those CSVs, no matching. Verification, then Q2's noise floor,
   then Q1 against it, then the mechanism decomposition.

## 11. Cost

Measured on this hardware, not extrapolated: `load_roi` + `to_channel` **7.5 s per ROI** (hoisted
out of the seed loop); `fused_response` **2.7 s**, `extract_peaks` **0.4 s**, NMS 0.2 s,
`score_detections` 0.9 s. Per seed: **2** response maps, **2** extractions, **2** NMS, **2** `od`
scorings — `pad_on_dropband` is a mask on `pad_on`'s frame and adds none of these. Arms:
3 × 6 z × 2 axes of `bucket_detections`, with `coverage_fraction` cached per (arm, z). From F5's
real 31.7 s/seed: **≈ 38 s/seed → ~44 min**, plus ~2 min of ROI loading. Budget an hour. `del` the
unpadded map after the `si == 0` identity check — two 39 MP float32 maps plus `rgb` and `hem` peaks
around 1–1.5 GB on the largest ROI.

## 12. Threats

| threat | handling |
|---|---|
| Coverage saturation makes `full_list_recall` uninformative | §6.1 — printed only with `coverage_frac`; guard extended to `read_99`/`read_100`/`FP_to_1.0` |
| A signed Δ with no attributable mechanism | §6.1's ledger and band counts — the whole point of revision 2 |
| Seed noise swamps the effect | §6.2 measured first, §6.3's two rulers, §4's MDE |
| Pseudo-replication from 5 seeds per ROI | §9 — ROI-clustered, n = 14 |
| Between-ROI heterogeneity read as click noise | §6.2 — within-ROI SD is the statistic; pooled-10 is labelled heterogeneity |
| Zero-variance denominators | §6.2 — reported `n/a`, never as a ratio |
| An axis declaration that decides the headline | §9's within-F6 head-to-head, per D5 |
| Under-power read as a null | §4's MDE, registered before the run |
| A tie-limited p read as significance | §9's floor rule, applied to §2a as well |

## 13. What review changed, and what draft 1 would have produced

| # | change | what draft 1 would have produced |
|---|---|---|
| 1 | `pad_on_offstats` → `pad_on_dropband`; draft 1's P2 withdrawn | An arm measuring a channel bounded at 2.27e-3 MAD — 3 candidates of 13,331 on the worst cell, **exactly 0 at K = 250** — presented as "the live claim", at the cost of a third response-map pipeline per seed. |
| 2 | §2c names suppression and top-K displacement | The syllogism "not reachability ⟹ re-normalisation", with a third mechanism in the gap; F6 would have attributed an artefact-competition effect to reachability. |
| 3 | §2b's ceiling corrected from ~4 to **0** | A ceiling four times too high on the mechanism v2 justified padding by. |
| 4 | Q2 computed from the existing CSV first; within-ROI SD replaces pooled-10 | A 45-minute re-derivation of a committed CSV, reporting a number inflated 10–41× by between-ROI heterogeneity as "click noise". |
| 5 | draft 1's P4 withdrawn | A prediction already false on 3 of 4 readings of existing data, whose "≥3×" clause pure noise satisfies 42 % of the time, and which mis-stated 013/529 as a domain. |
| 6 | §6.3's ruler split into statistical and practical | A rule that fires "smaller than the click" wherever the click matters and "bigger" wherever recall is pinned at 1.0. |
| 7 | §9's D5 head-to-head registered within F6 | An axis declaration resting on a head-to-head run *after* seeing both axes' answers — the forking path F5 §2.6 exists to close. |
| 8 | §6.1's ledger and band counts added | A signed Δ with no attribution, in a design whose two named mechanisms are both pre-measured at zero. |
| 9 | §4's MDE stated before the run | A marginally-powered null (MDE₈₀ ≈ 0.008 vs an expected −0.0071) read as "no effect". |
| 12 | **Revision 4:** the tie-aware floor applied *prospectively* to F6's own `tm_score` test (k ≥ 7 needed to reject at Holm's α/2; F5's base rate implies k = 5, floor 0.0625); P4's mechanism replaced by the quantitative product (band detections entering × TP density at the margin), which predicts the observed ΔTP on 4 of 4 pre-measured combinations. | P1's `tm_score` clause "confirmed" by a test that could not have failed; and a mechanism sentence that was **false as written** — it claimed a score-ranked list keeps band detections out of the top 250, and 5 of 250 are measurably in it. |
| 11 | **Revision 3:** P2 and P3 re-signed after review built the arms — suppression measured at 0 % of the total (registered ≥50 %), band occupancy at 100 % and negative (registered positive, <0.01; measured −0.02885). §3's no-emulation argument re-based on gate A's exact `n_detections` comparison rather than the dilation, which is defeatable with `np.where(valid, fused, _FLOOR)`. `dropband − pad_off` relabelled to include the inert re-normalisation channel. Interior-peak losses split pre/post-NMS; band detections separated from band annotations; P5 floored at 5 gained; gate 12 scoped to candidate sets; D5 contingency registered. | Two predictions with the wrong sign and magnitude, an irreducibility argument that does not hold, a mislabelled contrast repeating draft 1's error one arm to the left, and an axis declaration with no fallback in an experiment whose whole effect is on one axis. |
| 10 | Factual fixes: ex-201 **+0.0022** not +0.0019 (divided by 7, not 6); Δ`n_detections` p = **0.0156** is the smallest p, so "only nominally significant result anywhere" was false; both prior p-values relabelled floor-limited; the 402/459/548 TP losses scoped to `chromatin_od`; "byte-identical" → "identical to ~1e-6"; additivity demoted to a join check; F5's gate-scoping note restored; `od` padded in all arms; an un-edited draft artifact deleted from §9 | Five wrong or overstated numbers in a document meant to be binding. |
