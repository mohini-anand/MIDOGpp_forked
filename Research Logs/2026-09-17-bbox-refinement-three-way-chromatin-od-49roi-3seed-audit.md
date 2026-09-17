# Audit of `bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb`: every published number reproduces independently (50,008 values, 0 divergences) and the `hem_bbox` cost at K = 20/30 is robust; the `gray_bbox` (D8) cost against the 51 px box holds as a template effect at production's own gated click (1.90–2.24 points, Holm ≤ 0.014) but does not extend to dropping the Otsu gate, and "not measurably at K = 10" depends on how ROIs are weighted

**Scope.** Target: `production_hematoxylin_only/bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb` (29 cells,
18 code, 4 figures). Its 11 CSVs are `results/precision_at_k_49roi_3seed_chromatin_bbox3way_{per_run,per_roi,summary,by_domain,by_seed,delta_per_roi,delta_stats,bootstrap_ci,contested_sensitivity,top30,verification}.csv`
and its 4 saved figures are `production_hematoxylin_only/bbox3way_49roi3seed_*.png`. It also reads the committed
14-ROI CSVs `results/precision_at_k_14roi_prodseed_chromatin_bbox3way_{per_roi,top30,delta_per_roi,delta_stats}.csv`.
Documents read: its spec `BBOX3WAY_49ROI_3SEED_PROMPT.md`, `DECISIONS.md` D4, D5 (with its 2026-09-08 amendment), D7,
D9 and D10 (the uncommitted working-tree version), `D8_TEMPLATE_ANCHOR.md`, and both image manifests. Audit script:
`bbox_refinement_three_way_chromatin_od_49roi_3seed_audit.py` (repo root). Tables:
`results/bbox_refinement_three_way_chromatin_od_49roi_3seed_audit_*.csv`.
**Everything below is re-derived from `databases/MIDOG++.json`, the TIFF pixels and resolution tags in
`images/extra_valid/` and `images/extra_valid/testing_set/`, and the notebook's persisted `per_run` and `top30`
tables.** The matcher, the statistics, the click walk, the Otsu gate and a 15-run search pipeline are written in the
audit script. `midog_utils` is imported only inside Tier B, and only where labelled: production `build_seed` as a
consistency check on the gray walk, the own-draw and gray-gated sensitivity arms, and the re-timing.

## Conflict of interest

- **Signals I can see.** The notebook, its 11 CSVs, its 4 figures and its prompt are **untracked** (`??`, never
  committed), so no git author or age can be attached. mtimes: prompt 2026-09-16 21:05:56, kernel 22:25:06–22:54:11,
  notebook file 22:54:13. `DECISIONS.md` carries an uncommitted D10 dated today that names this notebook. Every
  imported module is committed (`c066829`) and clean. The work and this audit may come from the same project hand;
  only my fresh context separates them.
- **What limits the conflict.** No table is checked against itself. TP counts come from my matcher run against the
  JSON, never from the notebook's `bucket` column. p-values come from my own convolution. Clicks and template
  geometry come from my own RNG walk and my own Otsu gate on the pixels. 15 candidate lists come from my own
  `matchTemplate` → peaks → cap → NMS → self-hit → od51 implementation.
- **What it cannot cover.** The design choices themselves: the `chromatin_od` ranker, the joint gate, precision@K
  at K ≤ 30, three clicks, and the ROI as the unit. Those belong to `premise-reviewer` and to a reader who did not
  make them.

---

## Part 0 — what reproduces

**The engineering is sound. The published statistics are the statistics the code in the working tree computes.** An
independent re-implementation from pixels reproduces all 147 clicks and every template geometry, and reproduces the 15
candidate lists it re-ran coordinate for coordinate, with `od` bit-exact. Across the audit, **50,008 values were
compared with 0 divergences**: 43,120 independent and 6,888 consistency (`..._tier_a_b_comparisons.csv`). Of 174
printed or prose claims, **173 hold and 1 does not** (`..._prose_claims.csv`, F4; the file's 175th row is an
informational count the notebook does not claim).

| check | tier | values | divergences |
|---|---|---|---|
| Execution: 18 code cells, `execution_count` 1…18 contiguous; 0 errors; 0 unrun; last cell run; kernel timestamps monotonic 02:25:06 → 02:54:11 UTC (one session). Not vendored: it imports `midog_utils` and reads `images/`. | A | 18 cells | — |
| Provenance: newest module mtime 20:03:03, committed in `c066829` at 20:30:49, clean. Artifacts written 22:53:54–22:54:04 and figures 22:54:05–22:54:11, both after the code. **Tier A here is therefore an independent reproduction, not consistency with a stale artifact.** | A | 28 files | — |
| Composition: row counts 441/147/27/63/81/441/90/54/81/13,230/2,600 all equal the design (49 × 3 × 3 etc.); 0 duplicate keys; 14 + 35 ROIs equal the two directory listings; 7 domains × 7 ROIs, with domains confirmed from the JSON; ranks 0..29 in every run | A | 27 checks | 0 |
| **Re-scoring every top-30 detection** against `gt_eval` rebuilt from the JSON (all annotations minus the click), radius 7.5 µm / mpp from the TIFF tags (29.61–33.14 px): bucket and matched `ann_id`. My matcher saw only the persisted top 30 of each run while the notebook matched the full list, so this agreement is also the cell 2 prefix property, on 441/441 runs | A, independent | 26,460 | **0** |
| TP@10/20/30 for all 441 runs; `n_gt_mitotic`; nearest other annotation; `tpl_offset_px`; mpp; radius | A, independent | 3,528 | 0 |
| Seed and list properties: seed is a unanimous mitosis, click equals its box centre, border ≥ 36 px, 73 px patch readable at the template centre, odd base ≤ 51, `od` has no NaN and descends, `default_51` is 51 px at the click, Gate C (one click and one `n_gt` per run; 3 distinct clicks per ROI) | A, independent | 4,312 | 0 |
| `per_roi`, `summary` (pooled, worst/best ROI, worst click with ties), `by_domain`, `by_seed`, `delta_per_roi` | A, independent | 3,717 | 0 |
| `delta_stats`: pooled delta, W/L/T, `n_nonzero`, **exact p (my own polynomial convolution)**, floor, Holm (18 family rows), majority flags | A, independent | 792 | 0 |
| Monte Carlo p within 5 SE of exact (90); bootstrap CIs and shares re-drawn with the notebook's RNG spec (108) | A, consistency | 198 | 0 |
| Contested-mitosis sensitivity: my contested flag from JSON votes; pooled delta, W/L/T, exact p and Holm for two rules; my own 10⁶-draw Monte Carlo for `contested_excluded` | A, independent | 441 | 0 |
| Contested bootstrap CIs and D5 flags with the notebook's RNG spec | A, consistency | 243 | 0 |
| Refusal counts, `n_retries`, Gate B from persisted lists (17 cases), top-30 annulus from coordinates | A, independent | 588 | 0 |
| Gate A: seed-0 `original_14` top-30 against the committed 14-ROI top-30 | A, consistency | 6,300 | 0 |
| Figures, step 1a: bar heights (27), heatmap cells (441), domain cells (63), recomputed | A, independent | 531 | 0 |
| **Clicks from pixels**: my RNG walk + my Otsu gate + my `(x0+x1−1)/2` anchor reproduce `seed_ann_id`, `base_size`, centre x/y, `refused_draws` and `n_joint_valid` for all 147 clicks × 3 conditions | B, independent | 2,646 | 0 |
| **Search from pixels**, 15 runs (211 s0, 245 s1, 293 s0, 425 s1, 460 s1 × 3): `n_peaks`, `n_detections`, `n_self_hits`, identical top-30 coordinates, TP@K; `od` max relative difference 0.0 | B, independent | 105 | 0 |
| Production `build_seed(gray)` vs my gray-alone walk (147 clicks, earlier clicks fixed) | B, consistency | 147 | 0 |

**Gates.**
- *Execution:* applied, passes.
- *Provenance:* applied, clean for code. Everything the notebook wrote is untracked (Part 3).
- *Composition:* applied, passes.
- *Vendored exemption:* does not apply.

**Claims with no persisted artifact underneath.**
- *Template centres, refusals and joint-pool sizes.* These are persisted in `per_run`, and were re-derived from
  pixels anyway (Tier B).
- *Full-list properties* (`n_peaks = 100`, the annulus beyond rank 30, `n_self_hits`). These live only as per-run
  counts in `per_run`. Tier B checked them on 15 runs.
- *Figures.* Checked by the two-step protocol.
- *Runtime lines.* Re-timed.

No claim needed `cannot check`.

**Checked and clear, not findings.**
- *Ceilings do not drive the ties.* No click sits at its one-to-one precision ceiling `min(n_gt, K)/K` at K = 20. One
  `hem_bbox` click does at K = 30. Only `211.tiff` ties by construction, with 0 TP on all 9 runs
  (`..._ext_ceiling.csv`, `..._ext_ties.csv`).
- *Step 4.4.* The length-matched random-list precision is a median 0.0032 per ROI (max 0.021), against observed
  pooled 0.33–0.60 (`..._ext_precision_null.csv`). No claim is made against chance.
- *211.tiff's zero is real, not misalignment.* Its re-implemented full lists contain 2, 2 and 0 of 19 mitoses, and 0
  inside the top 30 even under `tm_score` ordering (`..._tierb_pipeline_reimpl.csv`).

---

## Part 1 — findings

No Tier 1 finding. The arithmetic and the measurement are right. What follows is about what the numbers support.

### Tier 2

#### F1 — The `gray_bbox` (D8) cost against `default_51` holds as a template effect at production's gated click; it does not extend to an ungated 51 px draw, where it shrinks by about 30% and loses Holm significance in `all_49`. Every `hem_bbox` contrast holds

**The claim.** Cell 27: *"`gray_bbox`, the production template (D8), costs 1.95-2.35 points against the untightened 51
px box"*. All three K = 20/30 contrasts *"resolve: Holm p <= 0.01 and D5's bar is met."*

**The caveat the notebook already gives.** It says clicks come from the joint gate, and counts how many would
differ: 26/147 for `default_51` alone, 10/147 for production `gray_bbox`. Both counts re-derive from pixels, and
production's own `build_seed` reproduces the gray walk on 147/147. What the notebook does not report is which way
that displacement moves the result, or by how much.

**Measured.** I re-ran the production entry point on each condition's own-gate click wherever it differs from the
joint click: 26 `default_51` runs, 10 `gray_bbox` and 2 `hem_bbox`. For seeds 1–2, earlier clicks are held at the
joint choices, as in the notebook's own caveat. I then recomputed the notebook's ROI-level exact sign-flip, with Holm
over the same 9 tests (`..._tierb_own_draw_sensitivity.csv`, `..._tierb_own_draw_runs.csv`):

| family | pair | K | joint gate (notebook): pts, W/L/T, p, Holm | own draw: pts, W/L/T, p, Holm |
|---|---|---|---|---|
| all_49 | `gray − default` | 20 | −2.35, 12/32/5, 0.0025, **0.0100** | −1.67, 14/28/7, 0.038, **0.18** |
| all_49 | `gray − default` | 30 | −1.95, 9/29/11, 0.0014, **0.0068** | −1.41, 13/30/6, 0.036, **0.18** |
| testing_35 | `gray − default` | 20 | −2.43, 8/23/4, 0.0081, 0.032 | −2.00, 10/21/4, 0.038, **0.15** |
| testing_35 | `gray − default` | 30 | −2.44, 4/21/10, 0.0002, 0.0014 | −2.03, 7/22/6, 0.0074, 0.037 |
| all_49 | `hem − default` | 20 / 30 | −5.65 / −5.33, Holm 3.9e-5 / 4.5e-6 | −4.73 / −4.51, Holm 7.2e-4 / 3.0e-4 |
| all_49 | `hem − gray` | 20 / 30 | −3.30 / −3.38, Holm 3.1e-4 / 1.7e-5 | −3.06 / −3.11, Holm 6.4e-4 / 1.7e-4 |

**Why it moves.** On the 26 clicks the gates refused but an ungated 51 px box would have taken, `default_51` does
worse than on the joint replacement: ΔTP@20 = −23 and ΔTP@30 = −33, a mean of −0.88 and −1.27 per click. So the
joint gate removes clicks that are poor seeds for the untightened box too, which flatters `default_51` when it is
read as a gate-free policy. `gray_bbox` on its own 10 clicks moves −3 and −9
(`..._tierb_own_draw_vs_joint_by_condition.csv`).

**The gated half survives.** A third arm puts both `gray_bbox` and `default_51` on production's own gray-gated click,
wherever that differs from the joint click: 10 clicks, 1 new pipeline run. That is the template contrast under
production's gate, and it stays resolved:

| family | K | pts | W/L/T | p | Holm |
|---|---|---|---|---|---|
| all_49 | 20 | −2.24 | 12/30/7 | 0.0034 | 0.014 |
| all_49 | 30 | −1.90 | 10/29/10 | 0.0015 | 0.0074 |
| testing_35 | 20 | −2.48 | 9/22/4 | 0.0092 | 0.037 |
| testing_35 | 30 | −2.44 | 5/21/9 | 0.0004 | 0.0026 |

So what fails is the comparison against a 51 px box drawn without any Otsu gate, not the template comparison
itself.

**Two corroborating fragilities for the same contrast.**
- *Fragility index.* The greedy count of single-click ±1 TP changes that lifts Holm above 0.05 is **6** (K = 20) and
  **9** (K = 30) in `all_49`, and **2** at K = 20 in `testing_35`. The `hem_bbox` contrasts need 16–54 in `all_49` and
  8–36 in `testing_35` (`..._ext_fragility.csv`). This is an upper bound on the true minimum. D10 (decided today, not applied) predicts
  exactly this class of perturbation under `chromatin_od`: ±1 TP on 9 of 770 stored rows.
- *Where the cost lives.* It sits in the three human domains: −1.43/−5.48/−5.24, −3.81/−4.76/−3.65 and
  −5.71/−3.81/−3.33 points. The four canine domains are near zero. Human minus canine is −4.96 / −4.09 / −3.72 points
  at K = 10/20/30. The one-sided domain-permutation p is 1/35 = 0.029 at every K, which is the floor for a 3-vs-4
  split (`..._ext_species_split.csv`). Species is largely confounded with scanner here (canine: 3D Histech and Aperio CS2;
  human: mostly Hamamatsu). The pooled −1.95 to −2.35 describes
  neither group.

**Caveats on my numbers.**
- *The own-draw arm unpairs up to 26 clicks.* That adds click-to-click noise, so part of the p increase is lost power,
  not bias. The direction evidence is the point estimate shrinking 29% (K = 20) and 28% (K = 30).
- *It approximates each policy.* Earlier clicks are held at the joint choices, and an "ungated `default_51`" is a
  hypothetical production path.
- *It is still ROI-level and exact.* The test is the same one, on 49 and 35 units.

**Supported statement.**
- *On gated clicks the cost holds.* When both templates are cut at a click that passes the Otsu gate (joint, or
  production's gray gate), `gray_bbox` costs 1.9–2.5 points against a click-centred 51 px box at K = 20/30, Holm
  ≤ 0.037.
- *Against an ungated 51 px pipeline it does not resolve.* With each drawing its own click, the cost is 1.41–1.67
  points in `all_49` at Holm 0.18. Only `testing_35` K = 30 stays below 0.05 (Holm 0.037).
- *So the ~2-point gain requires keeping the gate.* It belongs to "51 px with the Otsu refusal", not to "drop the
  gate and cut a 51 px box".
- *Scope.* The contrast is concentrated in human domains, and it is the fragile part of the headline.
- *The `hem_bbox` conclusions are not affected.*

#### F2 — "Not measurably at K = 10" is a statement about precision-weighted pooling. Under D4's recall@K, and under a t interval, `hem_bbox`'s K = 10 cost resolves

**The claim.** Cell 27: *"the template refinement matters at K = 20 and K = 30, and not measurably at K = 10 … No
contrast survives Holm, and none meets D5's bar."*

**Measured** (`..._ext_recall_sensitivity.csv`, `..._inference.csv`, `..._bootstrap_seed_sensitivity.csv`):

| K = 10 contrast | precision Holm (notebook) | recall@K Holm (Σ ΔTP / n_gt per ROI, 10⁶-draw sign-flip) | bootstrap CI, pts (notebook seed) | t CI, pts, G − 1 df | share of 200 other bootstrap seeds excluding 0 |
|---|---|---|---|---|---|
| all_49 `hem − default` | 0.10 | **0.0062** | [−5.78, **+0.00**] | [−5.97, **−0.02**] | **28%** |
| all_49 `hem − gray` | 0.10 | **0.0099** | [−4.15, −0.34] | [−4.15, −0.21] | 100% (fails majority, 23/49) |
| testing_35 `hem − default` | 0.19 | **0.012** | [−6.57, +0.48] | [−6.97, +0.49] | 0% |
| testing_35 `hem − gray` | 0.17 | **0.012** | [−4.38, −0.10] | [−4.56, −0.02] | 74% |
| all_49 `gray − default` | 0.40 | 0.074 | [−2.52, +0.95] | [−2.59, +0.96] | 0% |

- **The ROI signs do not change; the weights do.** W/L/T is identical under recall and precision.
- **Three dense canine ROIs carry the null.** At K = 10 the positive side of `hem − default` totals 32 TP, and 20 of
  those (62.5%) come from `293.tiff` (+12, n_gt 104), `220.tiff` (+4, 76) and `291.tiff` (+4, 86)
  (`..._ext_k10_drivers.csv`). The losses are spread over 25 ROIs. Recall weighting down-weights dense ROIs.
- **D5's bar at K = 10 in `all_49` turns on one bootstrap seed.** Under a t interval, or in 28% of other seeds,
  `hem − default` meets D5's bar (majority 25/49 holds).
- **The majority verdicts survive dropping 211.tiff.** It is the one structurally tied ROI, and removing it changes
  no majority verdict (`majority_excl_structural_zero`).

**Supported statement.** At K = 10, `gray_bbox` vs `default_51` is unresolved under every weighting. `hem_bbox` is
worse in direction (−2.2 to −3.2 points). Whether that cost is "measurable" depends on the metric: precision-pooled
Holm 0.10, recall@K Holm 0.006–0.012 (the reporting unit D4 names). Say "not resolved under precision@K pooling"
rather than "not measurably".

#### F3 — Two of the four `original_14` "meets D5's bar" verdicts rest on a percentile bootstrap at G = 14, and fail under a t interval

Cell 27 says: *"`hem_bbox - default_51` and `hem_bbox - gray_bbox` meet D5's bar at K = 20 and 30."* At G = 14 the
cluster bootstrap is below asymptotic validity. The notebook itself calls those CIs "coarse … descriptive", and then
uses them for the bar anyway.

| original_14 | exact p | bootstrap CI (pts) | t CI, 13 df (pts) | D5 bar, bootstrap → t |
|---|---|---|---|---|
| `hem − default` K=20 | 0.013 | [−12.62, −1.90] | [−12.66, −0.19] | met → met |
| `hem − default` K=30 | 0.042 | [−9.92, −0.95] | [−10.11, **+0.11**] | met → **not met** |
| `hem − gray` K=20 | 0.046 | [−9.88, −0.60] | [−9.74, **+1.17**] | met → **not met** |
| `hem − gray` K=30 | 0.018 | [−8.57, −1.19] | [−8.53, −0.05] | met → met (by 0.05 pts) |

The flips are estimator-driven, not seed-driven: all 200 alternative seeds exclude 0 for these four rows. The
`testing_35` `hem − gray` K = 10 "meets D5's bar" is the mirror case. Its bootstrap upper bound is −0.095 points
(= −1/1050, one TP), and 26% of other seeds include 0. Across the 27 subset rows, the t interval changes the D5
verdict on 3.

**Supported statement.** In `original_14`, `hem − default` K = 20 and `hem − gray` K = 30 meet the bar under both
estimators. The other two are borderline (exact p 0.042 and 0.046).

#### F4 — "`min_attainable_p` is at most 2.3e-10 in every row" is wrong by 8×

The largest `all_49` floor is **1.86e-9**, at `hem − gray` K = 10 with 30 nonzero ROIs; 2.3e-10 is the
`gray − default` K = 10 floor. The sentence's point, that ties do not limit these tests, still holds.

### Tier 3

- **T3-1 — `all_49` is not out-of-sample; `testing_35` is the confirmatory read.**
  - *Why:* seed 0 × `original_14` is the hypothesis-generating draw (Gate A bit-identity), and `all_49` includes it
    while being labelled "primary".
  - *What confirms:* `testing_35` is the only in-repo consumer of `images/extra_valid/testing_set/` (grep), and it
    confirms `hem_bbox` at K = 20/30 and `gray − default` at K = 20/30 (subject to F1).
  - *Fix:* call `testing_35` the confirmatory family.
- **T3-2 — "D5's bar" is D5's CI-and-majority clause only.** D5's full bar is 5 seeds × the 14 `extra_valid` ROIs,
  read at recall@K per D4. `D8_TEMPLATE_ANCHOR.md` restates it the same way. This run is 3 seeds × 49 ROIs at
  precision@K. Say "D5's CI-and-majority criterion".
- **T3-3 — Figure 1 invites the wrong reading.** Its unpaired per-condition CIs overlap heavily: `all_49` K = 20 is
  `default_51` [0.423, 0.569] against `hem_bbox` [0.370, 0.513]. Yet the paired `hem − default` CI is
  [−8.06, −3.37] points. A caption line pointing to Figure 4 and the paired table would prevent that.
- **T3-4 — Worst click.**
  - *The inversion is one cell.* In `original_14`, `hem_bbox`'s worst click is 0.10 against 0.00 at K = 10/20. That
    is `245.tiff` s1 alone.
  - *The D4-style metric confirms the headline.* The worst of each ROI's 3 clicks, averaged: `hem − default` at
    K = 30 is 0.331 vs 0.393, W/L/T 6/33/10, p 2.3e-5. `gray − default` is 0.420 vs 0.451 at K = 20 (p 0.0086) and
    0.370 vs 0.393 at K = 30 (p 0.015).
  - *The tail is concentrated.* 14 clicks lose ≥ 5 TP@30 under `hem_bbox` against 1 that gains. Per-click precision
    SD is unchanged across conditions (0.25–0.28) (`..._ext_product.csv`, `..._ext_click_spread.csv`).
  - *Worth adding.* This table is what D4's worst-click framing asks for, and the notebook does not report it.
- **T3-5 — Ties partly cancel across clicks.** At K = 30, 6 of 11 `gray − default` ties and 4–5 of the `hem` ties are
  ROIs whose clicks disagree and sum to 0. The ROI-level test treats these as no effect, which is correct for that
  unit but hides per-click disagreement.
- **T3-6 — Runtime lines.** "Pipeline runs took 28.7 min" covers the whole loop, including loads, channel conversion,
  the gate and scoring. `t_pipeline_s` sums to 26.2 min. It is one unreplicated wall-clock run. My three repeats of
  one call (300.tiff s0 `gray_bbox`) took 3.65 / 3.84 / 3.67 s in the final script run (3.51 / 3.08 / 3.04 s in an earlier run), against the notebook's 3.2 s, so run-to-run spread on this machine is about ±10%.
  Nothing rests on these numbers.
- **T3-7 — Shelf life.** D10 (uncommitted, 2026-09-17, decided and not applied) removes `SELF_HIT_RADIUS` and
  `n_self_hits`, which this notebook reads, and states that `chromatin_od` comparisons must not straddle the change.
  11 runs here removed no self-hit: `default_51` 2, `gray_bbox` 4, `hem_bbox` 5 (`..._ext_d10_exposure_no_self_hit_runs.csv`).
  Given F1's fragility index, re-run the `gray − default` contrast after D10 lands.

---

## Part 2 — verdict per conclusion

| # | cell | conclusion | verdict |
|---|---|---|---|
| 1 | 0, 1 | Only template size and centre change. Channel `hematoxylin_od`, od51, `run_production_pipeline` with no overrides (TM_CCOEFF, z = −1.5, MAX_PEAKS 100, NMS = match radius 7.5 µm, 5 px self-hit at the template centre); config-drift table all match | reproduces (code read; 15 lists re-implemented bit-exact from `(base, tx, ty)` alone) |
| 2 | 0 | `gray_bbox` is "what `build_seed` does in production today (D8)" | reproduces for the template: production `build_seed` agrees 147/147, and the anchor is `(x0+x1−1)/2`. Note that production's default ranker is `tm_score` (D5, `production.py`), and the joint gate moves gray's click on 10/147 |
| 3 | 0, 8 | Seed 0 on `original_14` is the 14-ROI draw; Gate A 42/42 rows and 1,260/1,260 top-30 rows, `od` bit-exact | reproduces (consistency with the committed artifact; the clicks also re-derive from pixels) |
| 4 | 0 | `testing_35` is the out-of-sample check | reproduces; `all_49` is not out-of-sample (T3-1) |
| 5 | 1, 4, 12 | 49 ROIs, 7 per domain; 441 runs | reproduces |
| 6 | 4, 28 | 441 runs in 28.7 min, median 3.6 s, max 5.5 s | reproduces as arithmetic; one run, and the loop total is not the pipeline sum (T3-6) |
| 7 | 6 | Halt checks: K delivered (`n_detections` 81–100), no NaN, `od`-descending, Gate C, distinct streams, joint pool ≥ 3 (min 6 on 201/364), contested tier on 0 clicks | reproduces (independent) |
| 8 | 6 | Report checks: cap binds 441/441; no candidate within one match radius 441/441; self-hits {0: 11, 1: 430} | reproduces (top-30 annulus on 441; full lists on 15) |
| 9 | 8 | Gate B: 17 identical-template cases (11 gray = default, 6 hem = gray), offsets {0: 14, 0.5: 2, 0.707: 1}, 0 divergences, geometry and SHA agree 17/17 | reproduces |
| 10 | 10 | Table A: medians 39/29; ranges 23–51 and 19–51; below-51 116/144; hem < gray 135, = 12, > 0; ceiling 31; per-subset values | reproduces (Tier A and pixels) |
| 11 | 10, 27 | `default_51` alone differs on 26/147 (12/14 by subset); production gray on 10/147 (6/4); refusals 16/10/2 | reproduces from pixels and production `build_seed` |
| 12 | 12 | Table B per-ROI precision | reproduces |
| 13 | 14 | Table C pooled precision, worst/best ROI, worst click with ties | reproduces |
| 14 | 15 | Per-domain and per-click-index tables | reproduces |
| 15 | 17 | Convolution equals brute force and the committed 14-ROI p (9/9) | reproduces (my convolution matches all 90 `delta_stats` p) |
| 16 | 18 | Monte Carlo agrees 90/90, max \|Δ\| 1.85e-3 | reproduces |
| 17 | 19, 27 | `all_49` table: pooled deltas, CIs, W/L/T, exact and Holm p, D5 flags | numbers reproduce; the K = 10 `hem − default` "D5 bar not met" is seed- and estimator-dependent (F2) |
| 18 | 19, 27 | `original_14`: `hem` pairs meet D5's bar at K = 20 and 30 (p 0.013/0.042; 0.046/0.018); `gray − default` CIs include 0 at every K | reproduces, overstated: 2 of the 4 fail under t, 13 df (F3) |
| 19 | 19, 27 | `testing_35`: all K = 20/30 contrasts resolve (Holm 0.0325/0.0014, 0.0014/2.8e-5, 0.0037/5.4e-4) with D5's bar met; no K = 10 Holm p below 0.17; `hem − gray` K = 10 meets D5's bar (CI [−4.38, −0.10], 18/35) | reproduces; the K = 10 bar is seed-fragile (26%, F3); `gray − default` K = 20 does not survive an ungated 51 px draw (Holm 0.15), though it holds on production-gray clicks (F1) |
| 20 | 19, 27 | Replication: same sign 9/9, CIs overlap 9/9 | reproduces |
| 21 | 21, 27 | Contested share 23.1–25.8%; all `all_49` K = 20/30 contrasts keep Holm ≤ 0.05 and D5's bar under all three rules; `testing_35` `hem − gray` K = 20 at Holm 0.065 under `contested_as_fp`; −5.65 → −4.52 | reproduces (exact p and Holm independent; MC within tolerance) |
| 22 | 23 | Figure 1: pooled precision with cluster-bootstrap CIs | reproduces (27 bars, 27 labels read off the render); presentation note T3-3 |
| 23 | 24 | Figure 2: per-ROI seed-summed delta heatmap | reproduces (441 cells; 8 rendered labels read; colour limit ±1.20) |
| 24 | 25 | Figure 3: per-domain delta with W/L/T | reproduces (63 cells; 6 rendered cells read) |
| 25 | 26 | Figure 4: `all_49` W/L/T with Holm labels | reproduces (9 stacks and 9 labels read off the render) |
| 26 | 27 | Order `default_51 > gray_bbox > hem_bbox` at every K, and per click index | reproduces (3/3 pooled, 9/9 per seed × K) |
| 27 | 27 | "Ties don't limit these tests: `min_attainable_p` ≤ 2.3e-10 in every row" | **does not reproduce** as a number: max 1.86e-9 (F4); the substance holds |
| 28 | 27 | **Headline, first half:** refinement matters at K = 20 and 30; all three contrasts resolve | reproduces for `hem − default` and `hem − gray`; **reproduces, overstated** for `gray − default`, which resolves on gated clicks (joint or production gray gate) but not against an ungated 51 px draw, is fragile, and is concentrated in human domains (F1) |
| 29 | 27 | **Headline, second half:** not measurably at K = 10 | **reproduces, overstated:** true of precision-pooled Holm; `hem_bbox`'s K = 10 cost resolves under recall@K (Holm 0.006–0.012) and under a t interval in `all_49` (F2) |
| 30 | 27 | `hem_bbox` costs 5.3–5.7 points vs `default_51` and 3.3–3.4 vs `gray_bbox` | reproduces (robust to own draws, contested rules and the fragility check) |
| 31 | 27 | `gray_bbox` (D8) costs 1.95–2.35 points vs the untightened 51 px box | reproduces, overstated in scope, not magnitude: as a one-factor template contrast on production's gray-gated click it is 1.90–2.24 points at Holm ≤ 0.014 (`all_49`), essentially the published number; it does not extend to an ungated 51 px draw, 1.41–1.67 points at Holm 0.18 (F1) |
| 32 | 27 | K = 10 details: `hem − default` CI upper end at 0, 2.5% of resamples ≥ 0; `hem − gray` excludes 0 but 23/49 | reproduces (0.0253); seed-specific (F2) |
| 33 | 27 | Per domain: human domains negative at every K (values), mast cell +0.48…+1.43, lung +1.43/0/0; `hem − default` negative 6/7 and 7/7; lung the only domain within 1.9 points | reproduces; the human/canine split sits at its 1/35 floor (F1) |
| 34 | 27 | 211.tiff 0 TP on all 9 runs, no sign information; 245.tiff worst in `original_14` (0.0833/0.1167/0.1000 at K = 20) | reproduces (211 confirmed from pixels) |
| 35 | 27 | Precision falls with median size 51/39/29; hem < gray on 135/147, never larger; offsets 2.55/2.50; not a dose-response | reproduces, correctly hedged. Also, `hem_bbox` size is strongly domain-structured (Kruskal–Wallis p = 1.3e-10, 37% of variance between domains) and `gray_bbox` size is not (p = 0.28, 5%) |
| 36 | 27 | Caveats: joint gate (26/147, 10/147), clicks without replacement, `chromatin_od` only, D9 validated on `tm_score`, K ≤ 30, bootstrap at 14 ROIs descriptive | reproduces; F1 measures the joint-gate caveat's size and direction; F3 shows the "descriptive" CIs were still used for verdicts |
| 37 | 2 | "Matching goes best-rank-first, so which annotations the top K claim never depends on anything ranked below K" | reproduces (independent): greedy matching of the top 30 alone gives the same bucket and matched `ann_id` as the notebook's full-list matching on all 13,230 rows of 441 runs (Part 0) |

`cannot check`: not used.

---

## Part 3 — what was re-run versus read

**Tiers.**
- **Tier A** (JSON, TIFF tags, the notebook's `per_run` and `top30`, the committed 14-ROI CSVs):
  - every table cell, statistic, flag and printed or prose number;
  - exact sign-flip by my own polynomial convolution;
  - Holm by my own step-down;
  - the bootstrap re-drawn under the notebook's RNG spec (labelled consistency) and under 200 other seeds;
  - t intervals on G − 1 df;
  - recall@K weighting (10⁶-draw sign-flip);
  - fragility, domain-as-unit, species split, ceilings and ties, the length-matched null, and the contested rules
    from JSON votes.
- **Tier B** (pixels, 4.9 min of compute in the final run):
  - All 49 ROIs: my `rgb2hed` and gray conversion, my Otsu gate and D8 anchor, and my RNG walk, giving all 147 clicks
    × 3 conditions. I chose this because every paired delta depends on the right click and the right template.
  - 15 re-implemented search runs, chosen for what a wrong answer would most damage: the structural zero (211 s0),
    the largest recentring (245 s1, 14.3 px), a near-ceiling ROI with all three at 51 px (425 s1), the largest
    `hem_bbox` loss (460 s1) and the largest `hem_bbox` gain (293 s0).
  - The own-draw and gray-gated arms through `production.run_production_pipeline`: 58 tagged rows in `..._tierb_own_draw_runs.csv` (38 own-draw, 20 gray-gated), which are 39 distinct pipeline calls after de-duplicating on (ROI, seed, condition, click). These are
    sensitivity analyses of the object under study, not reproductions.
  - Three repeats of one production call.
- **Tier C:** not run. No divergence needed explaining, and Tier B re-implements the whole measurement path on the
  runs that matter most.

**Script.** `bbox_refinement_three_way_chromatin_od_49roi_3seed_audit.py` ran start to finish on the final version
(exit 0, 5.8 min, of which Tier B 4.9 min; `TOTAL: 50008 values compared, 0 divergent`; claims 174/175, the one mismatch being F4) and regenerates every `results/bbox_refinement_three_way_chromatin_od_49roi_3seed_audit_*.csv`
cited here. Tier A does not import `midog_utils`.

### Step 3 — implementation review

**Config: no drift. The one deviation (rank key) is declared.**

| setting | notebook | `production.py` | `FSConfig` default | decision |
|---|---|---|---|---|
| search channel | `hematoxylin_od` | same | — | D3 |
| matcher | `TM_CCOEFF` | same | `TM_CCOEFF` | D1 (now in `DECISIONS_UNVERIFIED.md`) |
| peak min distance | 7 | 7 | 7 | — |
| deep floor z | −1.5 | −1.5 | `None` (required) | D9 "unchanged" |
| max peaks | 100 | 100 | 250,000 | D9 |
| NMS radius | `radius_px(mpp)`, asserts 7.5 µm | same | `None`, caller must supply | D7 |
| self-hit radius | 5.0 | 5.0 | 5.0 | D10 removes it (uncommitted, not applied) |
| border pad | — | `True` | `False` | — |
| od window | 51 | 51 | — | D5 notes it was never swept |
| template anchor | `tightened_template_box`, `(x0+x1−1)/2` | reads `template_xy` | — | D8 (`D8_TEMPLATE_ANCHOR.md`, current) |
| rank key | `chromatin_od` | default `tm_score` | — | D5: `tm_score` is the production ranker; declared in the title and caveats |

**Other implementation checks.**
- **Caps and floors.** `MAX_PEAKS` binds on 441/441. `n_detections` 81–100 ≥ 30, so K is always delivered.
- **Leakage.**
  - *The click's annotation:* excluded from `gt_eval`.
  - *Self-hit removal:* at the template centre.
  - *Annulus:* empty on the top-30 of all 441 runs (independent) and on full lists (15 runs).
- **Matching.** Rank-order greedy one-to-one with `d ≤ r`, per-image radius from `ResolutionUnit`-aware tags.
  Reproduced on 26,460 values.
- **RNG.**
  - *Streams:* `default_rng([s, image_id])`, with refused rows dropped and redrawn on the same stream, and earlier
    clicks removed.
  - *Pairing:* the same clicks serve all three arms.
  - *Reproduction:* the walk reproduces from pixels.
- **NaN and ties.** No NaN `od`. Stable `mergesort`. No `idxmax` on tied values anywhere.
- **Definitional drift, the high-yield check.** Nothing is mislabelled.
  - *`image_annotations(annotations, fn)`:* leaves `category_id=None`, so `gt_eval` holds both categories, as the
    prompt specifies. `n_gt_mitotic` filters to category 1, and only a category-1 match is `human_correct_label`.
  - *`unanimous`:* is `len(set(labels)) == 1`, so "contested" means a mitosis with a dissenting vote.
    `agreement_pool`'s `n_mitotic_votes == n_votes` is the same set for category 1.
  - *`n_near_click`:* is measured from the click, not the template.
  - *`gray_alone_differs`:* matches production `build_seed` on 147/147.
- **The notebook's own verification.** 2,600 records, all halt checks passed, report checks all true. They do not
  cover scoring correctness (Gate A is determinism against a run of the same code) or the click walk. Tier A and
  Tier B above cover both.

### Step 4 modes

1. **Unit.** ROI-level exact sign-flip.
   - `all_49`: G = 49, nonzero 30–45, floor ≤ 1.9e-9.
   - `testing_35`: G = 35, nonzero 22–31, floor ≤ 4.8e-7.
   - `original_14`: G = 14, nonzero 8–14, floor ≤ 7.8e-3.
   - Domain level: floor 1/64 at best (7 nonzero domains), 1/32 where a domain ties.
2. **Domain confound.** N/A for the treatment, which is assigned within the click. The size observation is
   domain-structured for `hem_bbox` (row 35).
3. **Recall triad.** No recall claim. `n_detections` is reported and there is no saturation (~98 candidates). I
   added recall@K as a sensitivity (F2).
4. **Length-matched null.** N/A, since nothing is compared with chance. Context in Part 0.
5. **Scope.** F1 (joint gate vs production draw), T3-1, T3-2. `chromatin_od` is not the default ranker, and the
   notebook declares that.
6. **Pre-registration.** The prompt (written 21:05, before the notebook existed at 21:26) serves as the spec.
   - *Adhered to:* conditions, configuration, K, gates A/B/C with halt vs report, DP convolution validated at 6 dp,
     MC ≥ 200k, B = 10,000, Holm over 9 on `all_49` and `testing_35`, the summed `contested_excluded`, per-domain and
     per-seed tables, the 4 figures, and the caveats.
   - *Deviations, both disclosed:* Gate B widened to rounded centres plus SHA (an improvement), and "majority" pinned
     to all ROIs.
7. **Decision table.** None pre-registered. D8's falsifier needs identical size, which this design does not have.
   The notebook applies D5's CI-and-majority clause (T3-2).
8. **Product.** The worst-of-3-click metric confirms the ordering. `hem_bbox` has a catastrophic-click tail
   (T3-4). The `gray_bbox` cost is concentrated in human domains (F1).
9. **Machinery.** Appropriate. The percentile bootstrap at G = 14 was used for verdicts (F3).
10. **In-sample selection.** No operating point is chosen; conditions and K were fixed in the prompt. The hypothesis
    came from `original_14` seed 0, and `testing_35` is clean (T3-1).
11. **Timing.** Operational printouts only (T3-6).
12. **Beyond the named modes.** I looked for and found: the direction of the joint-gate bias (F1); weighting
    dependence (F2); estimator and seed dependence of D5 verdicts (F3); fragility against D10's predicted
    perturbation; a structural-zero ROI (211, real); ceiling-driven ties (none); cancelling ties (T3-5);
    human/canine concentration; misalignment on 211 (none); contested MC agreement (yes).

**Family.**
- *How established.* `git log` (0982979 the 14-ROI three-way; a7279ba the `hembbox` rebuild; 2183a99 the testing
  set, 18:44), `ls results | grep prodseed`, and `grep -rl testing_set` / `bbox3way` over `.ipynb`, `.py` and `.md`.
- *Lineage.* This notebook is the only 49-ROI arm: a scale-up of the 14-ROI `bbox3way` notebook, itself preceded by
  `…_chromatin`, `…_halfpixfix` and `…_hembbox` at 14 ROIs × seed 0. The headline is not a best-of-N selection;
  conditions and K were fixed before the run.
- *Excluded as not-an-arm.* `bbox_refinement_three_way_detail_14roi.ipynb` is a visual companion.
  `tightening_process_hem_vs_gray.ipynb` and `production_seed_tightening_detail_14roi.ipynb` are mechanism figures
  with no precision output.

**Prior audits, read last.** None covers this notebook. Rounds 1–3 of the 14-ROI parent cover the committed CSVs it
gates against.
- *Confirms and extends round 1 T2-1* (production gray is off its own click): re-derived at 10/147, and F1 measures
  what it does to the result.
- *Extends round 1 T2-2*:
  - the 14-ROI `contested_as_fp` failure does not recur at 49 × 3; every `all_49` K = 20/30 contrast survives all
    three rules;
  - its domain-as-unit concern holds more strongly: at 7 domain units nothing can survive Holm (minimum 0.14,
    `..._ext_domain_as_unit.csv`).
- *Disagree in part with round 1 Part 3's* "recall@K is algebraically redundant with precision@K here". Per-ROI signs
  are identical, but across ROIs the 1/n_gt weights change the K = 10 verdicts (F2).
- *Confirms round 1 T3-5* (D9 validated on `tm_score` only); this notebook carries it as a caveat.
- *Confirms round 3 N1's fix:* Gate B now tests rounded centres and SHA.
- *Did not look:* round 1's `tm_score` arm (pooled ordering carries over; worst ROI reverses) at 49 ROIs, round 2's
  prose findings, and the `hembbox` and tightening-process audits (sibling notebooks).

**Appendix facts that no longer hold.**
- **`FSConfig.tm_method`** now defaults to `TM_CCOEFF`. It changed in `c066829` (20:30:49 on 2026-09-16), after
  round 1 recorded `TM_CCOEFF_NORMED`, which was true of its tree.
- **`FSConfig(nms_radius=None)`** now means "caller must supply", and `find_and_suppress` raises. It no longer means
  "this image's evaluation match radius".
- **The ~8 s single-seed pass** is now ~3 s per production call. The notebook's median is 3.6 s over 441 calls; my
  repeats took 3.0–3.8 s across two script runs, of which `rgb2hed` is 1.7–1.9 s. A re-implemented search after conversion takes 0.9–1.5 s.
- **"No held-out split anywhere"**: `images/extra_valid/testing_set/` (35 ROIs, added 2026-09-16) is held out from
  the 14.
- **Embedded PNGs**: 274 across 70 notebooks (the appendix said 217).
- **Results CSVs**: 496 `results/*.csv` on disk at the end of this audit (31 of them this audit's tables) and 367 tracked. Untracked ones include all 11 of this notebook's.
- **`DECISIONS.md` layout**: D1 now lives in `DECISIONS_UNVERIFIED.md` and D8 in `D8_TEMPLATE_ANCHOR.md`.
- *Still holds:* `image_annotations(category_id=None)`, and radii 29.6–33.1 px.

---

## Part 4 — premises this audit inherited

**ROI per stratum.** Exactly **7 ROIs in each of 7 tumour domains**: 2 in `original_14` and 5 in `testing_35`. So ROI
and domain are not collinear here, and Steps 4.1 and 4.2 are separate checks. The JSON's annotated images per tumour
type are: human breast 150, canine STS 100, canine lymphosarcoma 55, human NET 55, canine mast cell 50, human melanoma
49, canine lung 44. **Effective units:** `n_nonzero` 30–45 of 49, with one ROI (211.tiff) identically zero in every
test.

| premise | source | what would falsify it |
|---|---|---|
| ROI is the exchangeable unit | D5; F5 §8 | ROIs within a domain or scanner behaving as one unit. With domain as the unit, no contrast survives Holm (min 0.14 at 7 units), so every "resolves" verdict here depends on this premise. **Worth a premise review:** species and scanner are largely confounded, and F1's cost splits by species at the 1/35 floor. |
| precision@K pooled over clicks is the metric | the prompt; D9 | recall@K (D4's named unit) giving different verdicts. It does at K = 10 for `hem_bbox` (F2). |
| NMS radius = match radius = 7.5 µm | D7 | TP moving with radius. Not re-tested here; round 1 found 0/126 cells move at ±1 px. |
| Consensus `category_id` is truth; contested mitoses count | MIDOG++ labels | Verdicts flipping when contested hits are charged as FPs. At 49 × 3 they do not for `all_49` K = 20/30. |
| `chromatin_od` is the ranker worth evaluating | the prompt (D5 names `tm_score`) | The ordering reversing under `tm_score`. Not measured at 49 ROIs. |
| The joint gate is a fair common click set | the notebook's design | Own-draw results differing. They do for `gray − default` against an ungated 51 px draw, and do not on production-gray-gated clicks (F1). |
| `MAX_PEAKS = 100` does not change precision through K = 30 | D9 (`tm_score` arm only) | An uncapped `chromatin_od` run with different top-30 lists. Not run. |
| The pre-D10 self-hit filter is the pipeline | `production.py` today; D10 pending | D10 applied and `chromatin_od` cells moving. D10 predicts ±1 on about 1% of rows; F1's contrast has a fragility index of 2–9. |
| Prior audits' findings | `Research Logs/2026-09-16-*` | Same premises and the same day, so not independent corroboration. |

**Which verdicts would change if a premise failed.**
- *Domain as the unit:* every "resolves" verdict (rows 28, 30, 31) becomes unresolvable.
- *recall@K as the metric:* row 29's K = 10 null fails for `hem_bbox`.
- *An ungated 51 px pipeline as the comparator:* row 31 loses its Holm resolution in `all_49`, while on gated clicks it holds (F1).
- *D10 applied:* the `gray − default` rows are the ones at risk.
- *Unaffected:* every arithmetic verdict in Part 0.
