# Audit of `chromatin_od_ranker_variant.ipynb`: every number in it reproduces, including from raw pixels — but the "two-stage cascade" that is supposed to reconcile it with D5 is not what produces the precision gain, and the pad-timing outlier is misattributed

**Scope.** Target: `threshold_maxpeaks_ablation/chromatin_od_ranker_variant.ipynb` (untracked,
executed in place 2026-09-12 15:17). Artifacts: `threshold_maxpeaks_ablation/chromatin_od_ranker_timing.csv`,
`threshold_maxpeaks_ablation/chromatin_od_ranker_precision.csv` (both untracked),
`threshold_maxpeaks_ablation/max_peaks_100_{timing,precision}.csv`,
`results/precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv`,
`latency_profiling/chromatin_od_latency_per_roi.csv`, `databases/MIDOG++.json`, the
`images/extra_valid/*.tiff` headers and pixels. Decisions read: `DECISIONS.md` D1 (moved to
`DECISIONS_UNVERIFIED.md` today), D4, D5 + its 2026-09-08 amendment, D6, D7, D9 (dated today);
`D8_TEMPLATE_ANCHOR.md` read through its two superseded sections. Audit script:
`chromatin_od_ranker_variant_audit.py`; tables `results/chromatin_od_ranker_variant_audit_*.csv`.
**Everything below is re-derived from those artifacts and, for three ROIs, from the ROI pixels
themselves — never from the notebook's printed output.**

---

## Conflict of interest

`git log -1 --format='%an %ar' -- threshold_maxpeaks_ablation/chromatin_od_ranker_variant.ipynb`
returns nothing: the notebook has never been committed. `git status --porcelain` lists it and both
its CSVs as `??`, and the last five commits are all `mohini-anand` within the last hour. So the
work under audit and this audit are the same hand, minutes apart, and only fresh context separates
them.

What limits it: every table here is re-derived by `chromatin_od_ranker_variant_audit.py` from the
raw artifacts, never from the notebook's output cells; Section G re-implements peak extraction,
NMS, self-hit removal, `od51`, ranking and the greedy one-to-one matcher in plain
numpy/cv2 and never calls `midog_utils.compare`, `midog_utils.evaluate` or `midog_utils.nms`; the
cascade test in Section E uses an artifact (`precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv`)
committed on 2026-09-10, before this notebook existed.

What it cannot cover: the design premises listed in Part 4 — ROI as the exchangeable unit, 7.5 µm
as both NMS and match radius, precision@K on one click per ROI as the read-out. Those go to
`premise-reviewer` and to a reader who did not make them.

---

## Part 0 — what reproduces

The engineering is sound. **The published statistics are the statistics the code in the working
tree computes, and — for the three ROIs I re-ran from pixels — the statistics the pixels support.**
Every finding in Part 1 is about inference or attribution, not about the run.

| gate / check | result |
|---|---|
| Execution coherence | 14 code cells, `execution_count` **1…14 contiguous**, 0 error outputs, 0 unrun cells, no unrun last cell. Not vendored (imports `midog_utils`, reads `../images/extra_valid`, `../databases/MIDOG++.json`). |
| Composition gate | 20 checks, **20 pass**. 84 precision rows = 14 ROIs × 2 branches × 3 budgets; 14 timing rows; 0 duplicate keys; 0 duplicate rows; 0 NaN; 7 domains × exactly 2 ROIs each; `budget_delivered == budget` on all 84; `n_peaks == 100` on 14/14 (the pre-NMS cap binds, as intended); `n_detections` 89–99 on 14/14. |
| Ground truth, re-derived from `databases/MIDOG++.json` | 14/14 ROIs: `n_gt_mitotic` in the CSV equals the JSON's mitotic count **minus the seed annotation**; 14/14 `tumor_type` strings match; all 14 seeds are unanimous mitotic annotations. |
| Internal arithmetic of the CSV | `precision_at_budget == tp/budget_delivered` to 6.9e-17; `recall_at_budget == tp/n_gt_mitotic` to 9.7e-17, over all 84 rows. |
| **Verification 1, re-derived** | **112 values compared** against `max_peaks_100_{timing,precision}.csv` (`seed_ann_id`, `base_size`, `n_detections`, `map_median`, `mad_scale`, `tp_at_{10,20,30}` × 14 ROIs), **0 divergences.** The notebook's own Verification-1 cell printed its success branch at `ec=5`; this is that check re-run without it. |
| **Second, unclaimed oracle** | **42 values**: the `baseline` (tm_score) `tp_at_budget` also matches the **unbounded-pool** run `precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv` exactly on all 14 ROIs × 3 budgets, **0 divergences** — independent confirmation of D9's "the cap is free through K=30". |
| Pooled precision + win counts (claim 2) | reproduces exactly; see Part 2. |
| Decliner enumeration (claim 3) | reproduces exactly, all 4 pairs, and they are the only 4 of 42; see Part 2. |
| Timing decomposition (claim 1) | every mean/median/sd/min/max in cells 14, 24 and 26 reproduces from the timing CSV. |
| Sanity-check cell | `chromatin_od_nan_rate == 0.0` on 14/14; `chromatin_od_largest_tie_block == 1` on 14/14. Assert at `ec=6` executed and passed; `CHECKS` 70/70 at `ec=4`. |
| **Tier B — from pixels** | **18 `tp_at_budget` values** recomputed on `245.tiff`, `246.tiff`, `301.tiff` (all three decliner ROIs) with my own extraction/NMS/od51/matcher: **0 divergences**, both arms, all three budgets. `n_detections` 89/96/98 exact; `od51` NaN count 0; `od51` largest tie block 1. |
| Config three-way check | notebook cell 1 vs `FSConfig` defaults vs `DECISIONS.md`: **no stale override.** `NMS_RADIUS_UM = MATCH_RADIUS_UM = ev.MIDOG_RADIUS_UM = 7.5` (D7 — **not** the 5.0 µm that propagated elsewhere); `CHANNEL='hematoxylin_od'` → `channels.to_hematoxylin_od` → `chromatin.hematoxylin_od`, the *unclipped* OD the `chromatin_density` docstring requires (D3); `METHOD=cv2.TM_CCOEFF` passed **explicitly** to `tm.fused_response`, so `FSConfig.tm_method`'s default is never consulted; `PEAK_MIN_DISTANCE=7`, `SELF_HIT_RADIUS=5.0`, `scales=(1.0,)`, `n_angles=1`, `flips=(False,)`, `OTSU_WINDOW=51`, `DEEP_FLOOR_Z=-1.5` all match `D8_TEMPLATE_ANCHOR.md` and the sibling notebooks. Seed refinement is `ss.tightened_template_box` — the *current*, half-pixel-corrected anchor — and `suppress()` references the **returned centre** `tpl_xy`, which is D8's cost-1 requirement. |
| Gates not applicable | **Figure protocol: 0 embedded PNGs in the notebook** — there are no figures, so neither figure step applies. **Tier C: not run**; no divergence required it. **Leave-one-ROI-out (mode 10): no operating point is chosen on this data inside the notebook** — the two rankers are fixed a priori — so there is nothing to hold out. See finding 5 for the selection that happened *upstream*. |
| Provenance | `midog_utils/` is **clean** in `git status`, last touched `6d5cc01` (2026-09-10 13:28), two days before the notebook ran (15:17) — no module changed after the artifacts. Mtimes are ordered correctly: modules → notebook (15:17:03) → CSVs (15:17:02, same write). **But the notebook and both CSVs are untracked and have never been committed**, so they carry no commit provenance. Tier A against them is therefore *consistency with an uncommitted artifact* — which Section G upgrades to *independent reproduction* for the three ROIs re-run from pixels. No `.partial` file is read anywhere. |
| Claims with nothing persisted underneath | none. Every printed number traces to one of the two CSVs, and every CSV column traces to `databases/MIDOG++.json` or the ROI pixels. |

**`pool` is never reordered — the specific concern raised is unfounded.** `pool['od51'] = [...]`
adds a column; `pool.sort_values(...)` at `t6` and `t6c` return new frames (`top_tm`, `top_od51`)
used only for timing and never evaluated; `compare._rank` returns a fresh stable-mergesort copy.
Both arms receive the same `pool` object via a default-arg lambda, and `pool` is in NMS
score-descending order — asserted in-cell (`np.all(np.diff(s) <= 0)`). The stable tie-break is
intact, and with `largest_tie_block == 1` there are no ties to break.

---

## Part 1 — findings

### Tier 1

**Finding 1 — the "two-stage cascade" does not explain the precision gain. 74–122 % of it is
already present on the *uncapped* pool, in an artifact this repo committed two days ago.**

The notebook's closing claim (cells 22 and 26) is that the result "measures a two-stage cascade —
high tm_score AND high od51 — not od51 alone" and therefore "does not contradict D5's finding that
the two axes are statistically indistinguishable on the *unbounded* pool."

`results/precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv` is the *same pipeline with
the cap removed*: I verified `seed_ann_id`, `base_size`, `map_median`, `mad_scale` and
`n_gt_mitotic` are **identical on all 14 ROIs**, and its config cell pins `CHANNEL='hematoxylin_od'`,
`METHOD=cv2.TM_CCOEFF`, `PEAK_MIN_DISTANCE=7`, `SELF_HIT_RADIUS=5.0`, `DEEP_FLOOR_Z=-1.5`,
`NMS_RADIUS_UM=ev.MIDOG_RADIUS_UM`, `scales=(1.0,)`, `n_angles=1`, `flips=(False,)`,
`SEED_INDEX=0` — every value the same, with `MAX_PEAKS = 2_000_000` the only difference.

| K | unbounded tm | capped tm | unbounded od51 | capped od51 | od51 gain **uncapped** | od51 gain **capped** | gain the cap adds | share of the capped gain already present uncapped |
|---|---|---|---|---|---|---|---|---|
| 10 | 0.5000 | 0.5000 | **0.6214** | 0.6643 | **+0.1214** | +0.1643 | +0.0429 | **73.9 %** |
| 20 | 0.4536 | 0.4536 | **0.6107** | 0.5964 | **+0.1571** | +0.1428 | **−0.0143** | **110.0 %** |
| 30 | 0.4286 | 0.4286 | **0.5595** | 0.5357 | **+0.1309** | +0.1071 | **−0.0238** | **122.2 %** |

At K = 20 and K = 30 the cap makes the od51 ranking **worse**, not better. Per-ROI, cap-vs-no-cap
on the od51 arm: K=10 helps 4 / hurts 0 / ties 10; K=20 helps 4 / **hurts 6** / ties 4; K=30
helps 3 / **hurts 7** / ties 4.

Two consequences:

* **The mechanism half of the claim is true and the explanatory half is not.** `extract_peaks`
  does order by `np.lexsort((ys, xs, -scores))[:max_peaks]`, so the pool genuinely is the top-100
  by `tm_score`. But the cap is not what makes `od51` beat `tm_score` — `od51` beats it by roughly
  the same margin with no cap at all.
* **"Does not contradict D5" rests on a scope mismatch (Step 4 mode 5).** D5 point 3's null is
  **Δrecall@250**, on `results/f1_seed_sweep.csv`, **7 ROIs × 5 seeds**, at **z = 1.0**: Δ = +0.032,
  95 % CI [−0.047, +0.112], p = 0.36, positive on 3 of 7. D5 further says the advantage is "≈0 at
  budgets 25–100". None of that is a statement about **precision@10–30 on the 14 ROIs of
  `extra_valid` at the production seed**, which is this notebook's entire operating range — and on
  that measurement the repo's own committed unbounded-pool artifact already puts `od51` ahead by
  +0.121 / +0.157 / +0.131. There is a genuine tension to resolve, and the resolving variables are
  the ROI set, the seed count, the budget depth and z — not the cap.

**Correction to write.** Replace "it does not contradict D5's unbounded-pool null because the cap
pre-filters by the other axis" with something like: *"`od51` also beats `tm_score` at K = 10–30 on
the uncapped pool at these same seeds (`precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv`:
0.6214/0.6107/0.5595 against 0.5000/0.4536/0.4286), so the cap contributes at most +0.043 of the
K=10 gain and is net negative at K=20 and K=30. D5's null is Δrecall@250 on 7 ROIs × 5 seeds at
z = 1.0 and says nothing about precision@10–30; this result is outside its scope rather than
reconciled with it."*

**Caveat on my own finding.** The unbounded-pool comparison is single-seed (`SEED_INDEX=0`) on
14 ROIs, exactly like the notebook. It establishes that the cap is not the cause of the gain; it
does not establish that `od51` beats `tm_score` in general.

### Tier 2

**Finding 2 — the pad-timing outlier is attributed wrongly in two directions at once. ROI
dimensions explain the pad cost almost exactly, and the ~100 ms spike is not a property of the two
ROIs named.**

Cell 24 prints: *"Two ROIs run well above the rest: 013.tiff (113.1ms), 403.tiff (111.3ms) … Image
dimensions do not explain it: three other ROIs share the same exact pixel size as these two
(5412x7215) without the same spike, so this is reported unattributed rather than guessed at."*

The "three other ROIs" arithmetic is right — exactly 5 of 14 ROIs are 5412×7215 (013, 402, 403,
529, 548). The inference from it is wrong on both halves:

*Dimensions do explain the ordering.* All five 5412×7215 ROIs occupy the top five pad times in
**both** persisted runs:

| ROI | shape | Mpx | `t6a_od_pad_ms` (this run) | rank | pad ms (latency-profile run, pool 17 331) | rank |
|---|---|---|---|---|---|---|
| 013 | 5412×7215 | 39.05 | 107.68 | 1 | 95.38 | 1 |
| 403 | 5412×7215 | 39.05 | 106.63 | 2 | 90.42 | 2 |
| 402 | 5412×7215 | 39.05 | 36.87 | 3 | 76.66 | 4 |
| 529 | 5412×7215 | 39.05 | 33.60 | 4 | 77.40 | 3 |
| 548 | 5412×7215 | 39.05 | 31.23 | 5 | 76.55 | 5 |
| 094 | 5364×7103 | 38.10 | 29.69 | 6 | 29.54 | 6 |
| …9 smaller ROIs | 4835–4933 × 6447–6577 | 31.2–32.4 | 26.03–29.66 | 7–14 | 22.71–24.77 | 7–14 |

A perfect 5-vs-9 split by dimension has exact permutation probability **1/C(14,5) = 0.0005**, and
it occurred in **each run independently**. My own re-timing (5 repeats per ROI, 7 ROIs) puts the
median pad at **0.776–0.845 ms per megapixel** — essentially linear in area. So the three other
5412×7215 ROIs *did* run slower than every smaller ROI; they simply ran 31–37 ms rather than 107 ms.

*And the spike is not theirs.* Re-timed 5× per ROI, twice:

| ROI | Mpx | pad ms min | median | max | notebook single shot |
|---|---|---|---|---|---|
| 013 | 39.05 | 29.25 | 32.36 | **99.66** | 107.68 |
| 403 | 39.05 | 29.76 | 31.09 | **104.32** | 106.63 |
| **246** | **32.44** | 24.41 | 27.41 | **80.79** (96.06 in my first run) | 29.66 |
| 548 | 39.05 | 29.16 | 31.23 | 33.36 | 31.23 |
| 094 | 38.10 | 29.43 | 30.54 | 32.17 | 29.69 |
| 245 | 32.44 | 24.08 | 25.18 | 27.32 | 26.03 |
| 301 | 31.17 | 23.58 | 24.27 | 25.85 | 26.69 |

`246.tiff` — one of the nine *small* ROIs the notebook cites as spike-free — hits 81–96 ms on one
repeat in five, in both of my runs. The ~3.5× spike is an allocator/first-touch cost on a single
timed `copyMakeBorder` call, not a per-ROI property; the notebook's single-shot run simply did not
catch it on 246.

**Correction to write.** Replace the "unattributed" paragraph with: *"Pad cost tracks ROI area at
0.78–0.85 ms/Mpx: all five 5412×7215 ROIs pad slower than all nine smaller ones, in this run and
in `latency_profiling/chromatin_od_latency_per_roi.csv` independently (exact permutation
p = 0.0005 each). On top of that, a single `copyMakeBorder` call intermittently costs ~3.5× its
median — re-timed at 5 repeats, 013, 403 **and 246.tiff** all hit 81–104 ms on one repeat in five,
so the two ROIs singled out here are the ones the spike happened to land on in this single-shot
run, not ROIs with a property the others lack."*

**Finding 3 — every millisecond figure is n = 1 per ROI per stage, and the headline pad mean is
outlier-inflated by ~40 %.** (Step 4 mode 11.)

No stage is timed more than once. Across 5 repeats my within-ROI max/min ratio runs
**1.09×–3.51×**. Consequences for the reported figures:

| quantity | notebook (single shot) | repeat-stabilised (mean of per-ROI medians, 5 reps, 7 ROIs) |
|---|---|---|
| pad, `t6a` | **mean 40.39 ms** / median 29.44 ms | **28.87 ms** (mean of mins 27.10 ms) |
| stage-6 overhead `t6a+t6b+t6c` | **mean 45.47 ms** / median 34.53 ms | **≈ 33.9 ms** |
| stage-6 as % of baseline pipeline | **+3.28 %** | **+2.45 %** |

The notebook's own **median** (34.53 ms, 29.44 ms) is close to the repeat-stabilised figure; its
**mean** is not. Any summary that quotes a single figure should quote the median, or state that
the mean is two-ROI-driven. The criterion figure (`t6b+t6c`) is unaffected — sd 0.48 ms over 14
ROIs, so it is stable.

**Finding 4 — "the split held across two independent runs" and "bit-identical across both runs"
cannot be checked: only one run is persisted.** (Cell 26.) `threshold_maxpeaks_ablation/` holds
exactly three notebooks and their six CSVs, no `.partial`, no nbconvert log, no second timing file.
The determinism half is independently *supported* — Verification 1's 112/112, the unbounded run's
42/42, and my 18/18 from pixels all agree — but the second timing run itself is unverifiable from
the tree. Either persist it or drop the sentence.

**Finding 5 — selection breadth 6, upstream of this notebook.** (Step 4 mode 10.) `chromatin_od`
was not drawn from thin air: `precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv` already
scored **six** ranking axes on **these same 14 ROIs at these same budgets**:

| arm | P@10 | P@20 | P@30 |
|---|---|---|---|
| `od_contrast` | **0.6357** | **0.6179** | 0.5500 |
| `chromatin_od` (od51) | 0.6214 | 0.6107 | **0.5595** |
| `od31` | 0.5643 | 0.5500 | 0.5238 |
| `tm_score` | 0.5000 | 0.4536 | 0.4286 |
| `mask_od_mean` | 0.4071 | 0.4000 | 0.4119 |
| `od_falloff` | 0.3429 | 0.3464 | 0.3333 |

`chromatin_od` ranks 2nd of six at K=10 and K=20 and 1st at K=30. So the axis carried into this
notebook was chosen *in sample* from a six-way comparison on the same 14 ROIs the headline is then
reported on, and it is **not** the best of the six at two of the three budgets — D5's own
2026-09-08 amendment names `od_contrast` as the axis that satisfies the shortlist rule most
cheaply. That does not invalidate anything the notebook says (it compares exactly two arms and
says so), but a reader deciding "should we switch the ranker?" is entitled to know the comparison
set is six, not two, and that `od51` is not the winner of it at K=10/20.

**And the breadth is two-dimensional, not one.** The `chromatin_od_ranker_nms5um` sibling that
appeared mid-audit (see Part 3) measures the same quantity at a 5.0 µm suppression radius:
variant pooled precision **0.6536 / 0.5982 / 0.5333** against the target's
**0.6643 / 0.5964 / 0.5357**. So the headline 0.6643 at K = 10 is the **higher of two measured
values** for the same quantity, with no persisted statement anywhere of the comparison. The gap is
+0.011 and the direction reverses at K = 20 (5.0 µm is the higher there), so this is sizing rather
than a distortion — but the honest breadth is *six axes upstream × two suppression radii
laterally*, both scored on the same 14 ROIs, and at K = 10 the reported number is the maximum over
both.

### Tier 3

**Finding 6 — the `no_cap` invariant is vacuous as configured.** `caps=(MAX_PEAKS,) = (100,)` makes
`invariants.check_no_cap` assert that the **post-NMS** length is not 100. The cap is applied
**pre-NMS**, where it binds on 14/14 ROIs (`n_peaks == 100` in every row). The check therefore
passes by construction and contributes 28 of the 70 "checks passed". The notebook does state the
cap is meant to bind, so this is presentation, not a defect — but "70/70" reads stronger than it is.

**Finding 7 — `coverage_frac` is computed and then discarded.** `evaluate_arms` returns it;
`PRECISION_LONG`'s column selection drops it before `to_csv`, and no cell prints it. The recall
triad (D4 / D5) is otherwise complete — `n_detections`, `precision_at_budget` and
`recall_at_budget` all travel with recall. The omission is harmless here and favourable to the
notebook: 97 detections × π(31 px)² ≈ 0.8 % of a 32–39 Mpx ROI, against 0.930–0.975 on the
unbounded pool. Worth one line so a reader can see the pool is not tiling the ROI.

**Finding 8 — the notebook calls the private `cp._tie_and_nan`.** Works, matches
`evaluate_arms`'s own internal use, and my Tier B reimplementation agrees (tie block 1 on all
three ROIs). Cosmetic.

**Finding 9 — "+3.4 %" is not a number this notebook produces.** The notebook prints **+3.3 %**
twice (45.25/1384.86 = 3.267 %; 45.474/1384.86 = 3.284 %). Any summary quoting 3.4 % should say 3.3 %.

**Finding 10 — "~4.8 ms" is `t6b` alone, not the criterion.** `t6b_od51_loop_ms` mean = 4.680 ms;
`t6c_variant_rank_ms` mean = 0.398 ms; **`t6b + t6c` mean = 5.078 ms**, median 4.920 ms. The
notebook itself prints 5.078 and "~5.1ms/ROI" twice and is correct; a summary that labels 4.8 ms
as "computing od51 + re-sorting" has taken the loop-only figure.

**Finding 11 — the tie-block contrast is 1 vs 2–3, not 1 vs ~10.** On the unbounded pool at these
same 14 ROIs, `chromatin_od`'s `largest_tie_block` is **2 on 13 ROIs and 3 on 548.tiff**. The "~10"
figure comes from `chromatin.py`'s docstring, where it describes **002.tiff** (not in
`extra_valid`) and contrasts `channels.to_hematoxylin` (clipped, tie block 66) against
`chromatin.hematoxylin_od` (unclipped, tie block 10) — a different ROI and a different comparison
entirely. Separately: a `largest_tie_block` of 1 on a float statistic over 97 points means "no
exact ties", which is close to arithmetically inevitable and therefore carries little weight.

---

## Part 2 — verdict per conclusion

| # | Conclusion (cell) | Verdict | Evidence |
|---|---|---|---|
| N1 | The two branches rank one shared `max_peaks=100`, post-NMS pool; stages 1–5 run once (cell 0, 5) | **reproduces** | Code read; `n_peaks == 100`, `n_detections` 89–99 on 14/14; both arms take the same `pool` object. |
| N2 | Everything else held at D8/production defaults (cell 0) | **reproduces** | Three-way config check vs `FSConfig` and `DECISIONS.md`; no stale 5.0 µm radius, `TM_CCOEFF` explicit, `tightened_template_box` + self-hit on the returned centre. |
| N3 | The `baseline` branch reproduces `max_peaks_100_variant.ipynb` exactly (cell 9) | **reproduces** | 112 values re-compared, 0 divergences. Plus 42 more against the unbounded run, 0 divergences. |
| N4 | `od51` finite for every candidate; `largest_tie_block` 1–1 (cell 11) | **reproduces** | 14/14 in the CSV; confirmed from pixels on 3 ROIs. |
| N5 | Stage-6 means; +45.25 ms, +3.3 % full pipeline (cells 14, 22, 24, 26) | **reproduces, overstated** | Arithmetic exact. n = 1 per ROI per stage; my 5-repeat re-timing gives 28.87 ms for the pad against the reported 40.39 ms mean, i.e. +2.45 % rather than +3.28 %. Findings 2, 3. |
| N6 | Pooled precision 0.5000→0.6643 / 0.4536→0.5964 / 0.4286→0.5357; 10-1-3, 11-1-2, 12-2-0 (cells 20, 22, 26) | **reproduces** | Recomputed from the CSV; all three decliner ROIs re-derived from pixels. Strengthened by the ROI-level test below — and qualified by the domain-level one, where the two decliner ROIs turn out to be one entire domain. |
| N7 | "Image dimensions do not explain it … reported unattributed" (cell 24) | **does not reproduce** | Dimensions split pad time 5-vs-9 perfectly in two independent runs (p = 0.0005 each), at 0.78–0.85 ms/Mpx; the spike reproduces on 246.tiff too. Finding 2. |
| N8 | The pad is pool-size-independent; cross-check 45.5 ms vs 597.5 ms (cells 24, 26) | **reproduces** | 40.4 ms at pool 97 vs 45.5 ms at pool 17 331 (means); 29.4 vs 27.2 (medians). od51 loop 4.68 ms vs 597.5 ms, a 128× time ratio on a 179× pool ratio. |
| N9 | "the ~89 %/11 % split held across two independent runs"; "bit-identical across both runs" (cell 26) | **cannot check** | Only one run is persisted anywhere in the tree. The 89 %/11 % split itself reproduces (88.8 %/11.2 %) for the one run that exists. Finding 4. |
| N10a | `max_peaks=100` selects the top-100 *by `tm_score`* pre-NMS (cells 22, 26) | **reproduces** | `template_match.extract_peaks`: `np.lexsort((ys, xs, -scores))[:max_peaks]`. |
| N10b | "What this measures is a two-stage cascade … not od51 alone" (cell 26) | **reproduces, overstated** | True as a description of the pool. False as an explanation of the effect: 73.9 %/110.0 %/122.2 % of the gain is present with no cap. Finding 1. |
| N10c | "it does not contradict D5's finding that the two axes are statistically indistinguishable on the unbounded pool" (cells 22, 26) | **does not reproduce** | D5's null is Δrecall@250, 7 ROIs × 5 seeds, z = 1.0. On the unbounded pool at K=10/20/30 on these 14 ROIs, `od51` leads `tm_score` by +0.121/+0.157/+0.131. Scope mismatch, not reconciliation. Finding 1. |
| N11 | Single-seed, n = 14; does not meet D5's bar; not a decision on switching rankers (cells 22, 26) | **reproduces** | Correctly and repeatedly hedged. |
| N12 | Neither branch under-delivers K = 30 on any ROI (cell 16) | **reproduces** | `budget_delivered == budget` on all 84 rows. |
| N13 | 70/70 invariant checks passed (cell 7) | **reproduces** | 28 `no_cap` + 28 `nms_radius` + 14 `seed_annulus_empty`. The `no_cap` leg is vacuous as configured — Finding 6. |

### The test the notebook does not run, which makes N6 stronger

The notebook reports win counts and never converts them to a test. Because `budget_delivered == K`
on every ROI, pooled precision **is** the unweighted mean of per-ROI precision (verified to 6 dp),
so the point estimate is already ROI-exchangeable; only the uncertainty is missing. At G = 14, both
exact procedures enumerate instantly:

| K | mean per-ROI Δ | t CI95 (13 df) | up/down/tied | G effective | exact sign test | exact ROI sign-flip | floor |
|---|---|---|---|---|---|---|---|
| 10 | +0.1643 | [+0.0748, +0.2538] | 10 / 1 / 3 | 11 | p = 0.0117 | **p = 0.0049** | 1.2e-4 |
| 20 | +0.1429 | [+0.0653, +0.2204] | 11 / 1 / 2 | 12 | p = 0.0064 | **p = 0.0044** | 1.2e-4 |
| 30 | +0.1071 | [+0.0537, +0.1606] | 12 / 2 / 0 | 14 | p = 0.0129 | **p = 0.0028** | 1.2e-4 |

Exact two-sided sign-flip over the 14 ROI-level deltas (2¹⁴ = 16 384 enumerated, repo precedent per
the F1 independent audit), with the sign test on the non-tied units beside it. The effect is
resolvable at n = 14 and the CI excludes 0 at all three budgets. **This does not meet D5's bar** —
that bar is 5 seeds × 14 ROIs with the paired delta clustered at the ROI, and this is one click per
ROI, so the ROI-level variance here contains no click-to-click component. Given
`results/tm_ccoeff_headtohead_seed_variance.csv` records a click-to-click SD of `read_95` of 2 153
candidates on 245.tiff alone, that omitted component is not small.

**Clustered at the domain instead (Step 4.2), the picture is thinner and one domain reverses.**
With 2 ROIs per `tumor_type`, G falls from 14 to 7 and the two-sided sign-flip floor rises to
2/2⁷ = 0.0156:

| domain (2 ROIs each) | ΔP@10 | ΔP@20 | ΔP@30 |
|---|---|---|---|
| canine cutaneous mast cell tumor | +0.05 | +0.10 | +0.150 |
| canine lung cancer | +0.25 | +0.15 | +0.083 |
| **canine lymphosarcoma** | +0.05 | **−0.10** | **−0.067** |
| canine soft tissue sarcoma | +0.15 | +0.30 | +0.183 |
| human breast cancer | +0.20 | +0.25 | +0.150 |
| human melanoma | +0.30 | +0.15 | +0.133 |
| human neuroendocrine tumor | +0.15 | +0.15 | +0.117 |
| **exact sign-flip, G = 7** | **p = 0.0156** (at the floor, 7/7 up) | **p = 0.0469** (6/7) | **p = 0.0312** (6/7) |

The effect still clears 0.05 at all three budgets, but at K = 10 it sits **exactly on the
resolution floor** — 7 domains cannot produce a smaller two-sided p than 0.0156 — so "p = 0.0156"
there means *"as strong as G = 7 can resolve"*, not a measured tail.

**And the two decliner ROIs are not independent failures: they are one whole domain.** 245.tiff and
246.tiff are the repo's only two canine lymphosarcoma ROIs, and canine lymphosarcoma is the only
domain whose mean delta is negative at K = 20 and K = 30. The notebook reports "1 down of 14" and
"2 down of 14" as scattered exceptions; clustered at the exchangeable-above-ROI level they are
100 % of one of the seven domains. That is worth naming in the notebook's own prose, because a
domain-shaped failure is the shape that matters for a tool meant to travel across scanners and
tumour types.

---

## Part 3 — what was re-run versus read

**Tier A (recomputed from persisted artifacts, in `chromatin_od_ranker_variant_audit.py`):** the
composition gate (20 checks); ground truth re-derived from `databases/MIDOG++.json` (14 ROIs ×
5 quantities); the CSV's internal arithmetic (84 × 2); Verification 1 re-derived against
`max_peaks_100_*.csv` (112 values); a second, unclaimed oracle against
`precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv` (42 values); pooled precision, win
counts, exact sign test and exact 2¹⁴ sign-flip; the full decliner enumeration over all 42 (ROI,
budget) pairs; the whole timing decomposition; the capped-vs-unbounded cascade comparison (pooled
and per-ROI); the tie-block comparison; the six-axis selection-breadth table; ROI pixel dimensions
from TIFF headers.

**Tier B (~2 minutes of compute, well inside the 20-minute budget):** (i) `tp_at_budget` for both
arms at all three budgets on `245.tiff`, `246.tiff` and `301.tiff` — deliberately the three ROIs
where the variant *lost*, because a wrong answer there would most damage claim 3 — recomputed from
the pixels with my own `matchTemplate`→dilate peak extraction→greedy NMS→self-hit→`od51`→stable
rank→greedy one-to-one matcher, calling none of `compare.py`, `evaluate.py` or `nms.py`: **18/18
exact**, `n_detections` 3/3 exact, `map_median`/`mad_scale` reproduced. (ii) `cv2.copyMakeBorder`
re-timed 5× on each of 7 ROIs spanning 31.2–39.05 Mpx, twice.

**Tier C: not run.** No divergence required it, and the invoker did not ask. The measurements that
touch pixels are Section G of the audit script.

**Read, not re-run:** `midog_utils/{compare,evaluate,nms,channels,chromatin,template_match,find_and_suppress,seed_selection,dataset}.py`;
`DECISIONS.md` D4/D5/D5-amendment/D6/D7/D9 and `D8_TEMPLATE_ANCHOR.md` including both superseded
sections; `production_seed_precision_at_k_chromatin_half_pix_fix.ipynb`'s config cell.

**Beyond Step 4's named modes, I looked for and did not find:** in-place mutation of `pool` between
the two arms (none — Part 0); a category filter left at its permissive default
(`ds.image_annotations(anns, fn)` is called with `category_id=None` *deliberately* here, and the
full both-category frame is what `bucket_detections` requires so one detection cannot claim a
mitosis and a look-alike at once — this is the correct call, not the `tp_fp_separability` bug);
the seed annotation leaking into `gt_eval` (excluded by `ann_id`, and `n_gt_mitotic` = JSON count
minus 1 on 14/14); a self-hit counted as a TP (`suppress` removes within 5.0 px of the **returned
template centre**, and `seed_annulus_empty` passes 14/14); NaN sorting to the top (`na_position='last'`,
and nan_rate is 0 anyway); `budget_delivered` silently below K (0/84); rows lost to a merge (no
merge in the analysis path); the two arms receiving different `gt_eval` (they share one frame);
whether the clipped `to_hematoxylin` was used by mistake (it was not — `hematoxylin_od` resolves to
the unclipped `chromatin.hematoxylin_od`, which is why `largest_tie_block` is 1 rather than
saturated); and whether the notebook is one of several arms (it is not — see below).

**Family.** `threshold_maxpeaks_ablation/` was created whole in `5bf02c8` (2026-09-12 14:48,
*"Ablate max_peaks and the deep-floor z, verify their equivalence, adopt max_peaks=100"*), which
also carries D9. The target postdates that commit, is untracked, and is the only file in the repo
that reads or writes `chromatin_od_ranker_*` (`grep -rl "chromatin_od_ranker"` returns the notebook
and this audit script and nothing else). `git log --diff-filter=A -- threshold_maxpeaks_ablation/`
and the commits either side show no competing arm. **This is a follow-on measurement taken after
the cap was adopted, not the best of N arms** — the headline is a result, not a selection *within
this notebook*. The upstream six-axis selection is Finding 5.

**Family, amended mid-audit.** At 15:49 — after the family check above and while this audit was
running — a second arm appeared: `threshold_maxpeaks_ablation/chromatin_od_ranker_nms5um.ipynb`
plus `chromatin_od_ranker_nms5um_{timing,precision}.csv`. It is the **same measurement with the
suppression radius as the free variable**: `NMS_RADIUS_UM_CONTROL = ev.MIDOG_RADIUS_UM` (7.5) and
`NMS_RADIUS_UM_ABLATION = 5.0`, match radius explicitly untouched at 7.5, `MAX_PEAKS = 100`,
channel, method, `PEAK_MIN_DISTANCE`, `SELF_HIT_RADIUS` and `DEEP_FLOOR_Z` all identical to the
target. That is a **labelled ablation, not a D7 violation** — D7 forbids `NMS_RADIUS_UM = 5.0`
inherited *silently* into a new script, which is not what this is. Its control arm reproduces the
target's baseline bit-identically (0.5000 / 0.4536 / 0.4286). **State, as observed:** the CSVs are
timestamped 15:47:53 and the notebook 15:49:27 with all 16 `execution_count` null — the signature
of a run whose outputs were stripped on save. So the saved notebook cannot be scored for execution
coherence, and its CSVs carry no in-notebook provenance. That is an observation about the sibling,
not a defect of the target, and nothing in Parts 0–2 changes. Its consequence for the target is
folded into Finding 5.

**Prior audits.** Deferred to the end per protocol, then checked: `grep -l` over
`Research Logs/*audit*.md` for `chromatin_od_ranker`, `threshold_maxpeaks`, `max_peaks_100` and
`chromatin_od_latency` returns **nothing**. No prior audit covers this notebook or its artifacts,
so there is nothing to confirm, extend or overturn.

**The audit script runs start to finish, clean, in 31 s** under
`/Users/mohinianand/anaconda3/bin/python3`, and regenerates all 13 tables cited above. Tier B can
be skipped with `AUDIT_SKIP_TIER_B=1` (Tier A alone completes in ~3 s).

**Appendix facts that no longer hold.**

* **D1 is no longer in `DECISIONS.md`.** `f835885` (today, 14:54 — 23 minutes before the notebook
  ran) moved it to `DECISIONS_UNVERIFIED.md` "pending independent verification". The agent file's
  "`DECISIONS.md` D1 selects `TM_CCOEFF`" is now a *pending* decision, which matters because this
  notebook's config comment cites "D1 — unnormalized, contrast-sensitive" as settled.
* **`DECISIONS.md` now carries a D9** (`max_peaks = 100` before NMS, dated today) that the agent
  file's D1–D8 framing predates. D9 is the direct parent of this notebook's fixed cap.
* **`FSConfig.tm_method` still defaults to `TM_CCOEFF_NORMED`** while D1 selects `TM_CCOEFF` — the
  drift the agent file names is still live. It does not bite here: `METHOD` is passed explicitly to
  `tm.fused_response` and the dataclass default is never consulted.
* **`results/` is no longer fully tracked.** My own new tables are untracked, and so are the
  target notebook and both of its CSVs — the first untracked `threshold_maxpeaks_ablation/`
  artifacts since that directory was created this afternoon.
* **`images/extra_valid` ROI geometry, for the record:** 5 ROIs at 5412×7215 (013, 402, 403, 529,
  548), 1 at 5364×7103 (094), 6 at 4933×6577 (201, 233, 245, 246, 459, 460), 2 at 4835×6447 (300,
  301). This is not in the appendix and Finding 2 depends on it.
* Confirmed still true: `MIDOG_RADIUS_UM = 7.5`; `dataset.image_annotations(anns, fn, category_id=None)`;
  `template_match.BASE_SIZE = 51`, `PATCH_SIZE = 73`; the anaconda interpreter is the only one that
  imports `cv2`.

---

## Part 4 — premises this audit inherited

Each line is a **project decision I applied, not a fact I verified**.

| premise | source | what would falsify it |
|---|---|---|
| **ROI is the exchangeable unit** — used for every p-value and CI in Part 2 | `DECISIONS.md` D5; `Research Logs/2026-09-04-f5-preregistration.md` §8 | Evidence that per-ROI deltas are not exchangeable — e.g. a strong domain effect making ROIs within a domain more alike than across, which at 2 ROIs per domain would need a domain-level null instead. |
| **NMS radius = match radius = 7.5 µm** — makes `n_detections` 89–99 and every TP count what it is | `DECISIONS.md` D7 | D7's own exposure: the 6-mitosis proposal-stage loss at 7.5 µm reaching the top of the list. Gate 4 says it does not at K = 20, on one click per ROI. |
| **precision@K with `n_detections` and `recall_at_budget` beside it is the right read-out** | `DECISIONS.md` D4, D5 | If `coverage_frac` were large the recall numbers would be uninformative — here it is ≈0.8 %, so this premise is safe in this instance. |
| **`tightened_template_box` (half-pixel-corrected) is the production seed anchor** — fixes which template every number descends from | `D8_TEMPLATE_ANCHOR.md` (2026-09-10), read through both superseded sections | D8's own open question: a 5 seeds × 14 ROIs paired comparison showing recentring neutral or harmful at precision@K. D8 states plainly it adopts on mechanism, not outcome. |
| **`max_peaks = 100` is free through K = 30** — the premise that makes this a clean ranking ablation | `DECISIONS.md` D9 | Falsified for K > 30 by D9's own admission (the pool delivers 89–99, never 100). **I verified the K ≤ 30 half independently**: the capped `tm_score` `tp_at_budget` matches the unbounded run on 42/42 values. |
| **Reading burden and worst-seed behaviour are the product metrics** | AnnotateDx framing; `DECISIONS.md` D4 | Not exercised: this notebook reports median/pooled precision, not `read_95` or worst-seed. On one seed per ROI, worst-*seed* behaviour is unmeasurable by construction — worth stating as the notebook's main product-relevance gap. |
| **Prior audits in `Research Logs/`** | — | Not leaned on: none cover this material. |

**Measured ROI-per-stratum count:** `databases/MIDOG++.json` restricted to `images/extra_valid`
gives **exactly 2 ROIs per `tumor_type`, 7 domains, 14 ROIs** — verified against the JSON, not
assumed. So ROI is **not** collinear with the stratum here, and Step 4.1 (cluster by ROI) and
Step 4.2 (cluster by domain) are genuinely two different tests. But within-stratum replication is
n = 2, i.e. 1 degree of freedom per domain, so a domain-level null would have at most G = 7 and a
sign-flip floor of 1/2⁶ = 0.0156 — coarser than the ROI-level test I ran. All 14 ROI deltas are
non-zero at K = 30 and 11–12 of 14 at K = 10/20, so the effective unit count is close to the
nominal one; nothing is carried by ROIs that contribute zero.

**Which verdicts would change if a premise were wrong.** The verdicts that survive every premise
above are the reproduction ones — N3, N4, N6, N12 and Tier B's 18/18 are arithmetic and pixel
facts. **Finding 1 is the most robust of the findings**, because it compares two runs that share
every premise: whatever D7, D8 or D9 turn out to be, both arms of that comparison inherit it
identically, so "the cap is not what produces the gain" holds under any of them. Findings 2, 3
and 4 are timing facts and depend on no analytic premise at all. What *would* move is the
significance table under N6: if the ROI were not the exchangeable unit — if the right unit were the
domain — G falls from 14 to 7 and the two-sided sign-flip floor rises to 0.0156. I computed that
case rather than assuming it: the effect still clears 0.05 (p = 0.0156 / 0.0469 / 0.0312) but
K = 10 lands exactly on the floor, and **canine lymphosarcoma reverses sign at K = 20 and K = 30**
— so the direction does *not* survive uniformly, and the "1–2 ROIs down of 14" framing becomes
"1 domain down of 7". The conclusion "od51 beats tm_score at K = 10–30 on this sample" survives the
premise change; the strength attached to it, and the reading of the decliners as isolated, do not.
That premise is worth a `premise-reviewer` pass on its own, since this repo's entire inferential
apparatus rests on it.
