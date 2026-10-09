# Floor-clip stress test, Stage B: does the clip's gain survive stain and scanner change? (plan)

**Status:** plan v3, final, 2026-10-09. The reviewer judged it CONVERGED in round 3 (see §5). Not built, not run.

**DEFERRED (user decision, 2026-10-09) to a pre-ship robustness check.** It is not a precision lever, for four reasons:
- **The likely verdict is "keep 1.0×".** Stage A's predicted losses (≤ 0.41 pt P@10, about 12 % of the clip's 3.3-pt
  gain) sit well inside this plan's 1.1-pt tolerance.
- **The best case is small.** Switching floors could recover at most about 0.4 pt under H06.
- **The candidates are underpowered.** L150's minimum detectable edge (about 0.97 pt) is more than twice that.
  W128's predicted edge (about 0.09 pt) is below its own minimum (0.31 pt).
- **Unperturbed, the candidates already tie 1.0×:** 1.5× at +0.27 pt and W128 at +0.24, neither significant.

Run it before shipping if robustness to other labs must be certified. Real cross-scanner slides (Stage C) would answer
that better. The review scripts behind §5 were scratch files and were deleted on 2026-10-09; their numbers are
recorded in §5.

**Sources**
- `FUTURE_IDEAS.md` §4.1, Stage B: the spec.
- `TM_FLOOR_STRESS_STAGE_A_PLAN.md` and `tm_floor_stress_stageA_49roi.ipynb` (commit 3f4a0ef): the screen that triggered
  this stage.
- `tm_floor_definition_49roi_6seed.ipynb` and `tm_border_floor_49roi_6seed.ipynb`: the recorded search path, its arms
  and its statistics.

---

## 1. Purpose

**What is known.**
- The floor clip `hem_c = max(hem, b)` at `b = b_otsu` is the adopted change (D3 amendment).
- Unperturbed, on the 49 ROIs × 6 recorded clicks, it gains over no clip (`A3 − A1`, Intel):
  - P@10 +3.30 pts [2.11, 4.63];
  - P@20 +2.62;
  - +1.52 pool mitoses per click.

**What Stage A showed.** Stage A checked only where the floor lands; it ran no searches. Under realistic synthetic
perturbations it predicted ≥ 0.25 pt of P@10 loss in three cells:
- H06 (hematoxylin ×0.6) × `b_otsu`: 0.41;
- H06 × W128: 0.33;
- E15 (eosin ×1.5) × `b_otsu`: 0.28.

**What Stage B answers, by running the searches on the perturbed images:**
1. **Carry-over.** Under each perturbation, does the adopted clip still beat no clip, with its gain not shrunk by more
   than a third? (§4.1's decision.) There are two tiers of evidence:
   - **"Broken"** needs one of:
     - a shrink *significantly beyond* a third (C2-broken);
     - a gain no longer significantly above no clip (C1);
     - significantly grown per-click harm (C3).
     Only "broken" can conclude "rework the floor".
   - **"Weakened"** is a significant shrink whose point estimate lies beyond a third. It only opens the comparison with
     the alternatives, and an alternative must then beat the adopted floor significantly.
2. **Candidates.** Under each perturbation, does an alternative floor do significantly better than the adopted one? The
   alternatives:
   - W128, a 128 µm click-centred window, product-feasible;
   - L150, 1.5 × `b_otsu`, for headroom against undershoot.
   This decides whether to keep 1.0× `b_otsu`, or to recommend W128 or 1.5×.
3. **Stage A's assumption (descriptive).** On like-for-like boxes, is the measured loss close to Stage A's PL?

**What it cannot answer:**
- real cross-scanner data (Stage C);
- held-out ROIs;
- later loop rounds;
- chromatin od.

**The production path is not the one tested.** This tests the **research** search path. Production's
`annotation_frame` (on `feature/annotation-frame`) has no clip, no `mo_pad61` and no multi-Otsu refinement: its gate is
a 2-class Otsu on `gray_inverted`. Porting the clip will need its own gate.

## 2. Design

### 2.1 Units

- **ROIs:** the 49 ROIs of `tm_floor_definition`.
- **Clicks:** each ROI keeps its 6 recorded clicks (the `A0` rows of `results/tm_border_floor_49roi_6seed_per_run.csv`),
  plus up to 6 new clicks drawn by this fixed rule, so up to 12 per ROI:
  - **Pool:** the ROI's unanimous-agreement mitoses, excluding the 6 recorded seeds, sorted by `ann_id`. These are the
    recorded draw's own tier (`seed_selection.agreement_pool` returns the unanimous tier on every ROI).
  - **Eligibility:** each click is at the box centre `(cx, cy)` and must pass both gates:
    - **the recorded draw's joint gate,** as in `fill_geometry_precision_at_10_49roi_6seed.ipynb`:
      - `seed_selection.border_filter` at `PATCH_SIZE // 2`;
      - `seed_selection.tightened_template_box(·, cx, cy, otsu_window=51)` must succeed on both `gray_inverted` and
        `hematoxylin_od`, with each tightened box's 73 px patch readable;
      - the click's own 73 px patch must be readable;
    - **acceptance by `refine_template_box`** on the unperturbed `hem`.
  - **Order:** permuted with `np.random.default_rng([20261009, roi_index])`, where `roi_index` is the position in the
    sorted 49-file list. Walk the order and take eligible clicks until 6 are accepted or the pool runs out.
  - **Gate 10** checks that this joint gate reproduces the recorded `n_joint_valid` on 49/49 ROIs, and that all 294
    recorded seeds are eligible.
- **Expected yield:** 545 clicks (1.85×). 201 and 364 keep 6 clicks, and 7 ROIs have 6–8.
  - The pooled estimand is click-weighted. Review measured n_eff = 47.8 of 49 ROIs, domain weights of 12.7–15.4 %
    against 14.3 %, and a shift in G₀ of 0.07–0.16 pt. The margins stay absolute.
- **The roles differ by click set:**
  - the recorded 6 carry the Intel parity (gate 3) and the Stage A coupling (gate 4);
  - all clicks carry the decision contrasts.
- **Ground truth:** MIDOG++ annotations minus the click's own; match radius 7.5 µm (D7).
- **Search path:** `tm_floor_definition` cell 11 `search_arm`, unchanged:
  - the `mo_pad61` template from the clipped channel;
  - the D11 51 px blank filled with `channel.min()` (so it equals `b` when clipped);
  - `cv2.TM_CCOEFF`;
  - the 30 px edge band, invalid and floored;
  - median − 1.5·MAD threshold;
  - `MAX_PEAKS` 100;
  - NMS at `ev.radius_px(mpp)`;
  - ranking by `tm_score`.

### 2.2 Conditions (images)

Every perturbed image comes from the committed Stage A worker (`make_rgb_variant`, same per-ROI Macenko seeds), so it is
identical to Stage A's (gate 4).

| condition | role | definition (Stage A §3.2) |
|---|---|---|
| `UNPERT` | reference | the original ROI, re-searched on this Mac |
| `JPG95` | phase-1 calibration (near-null) | JPEG quality 95, 4:2:0 |
| `NULL2` | phase-1 harness gate | `hem' = 2·hem`, channel level; click 0 of every ROI |
| `H06` | decision; fired in Stage A | hematoxylin ×0.6 |
| `E15` | decision; fired in Stage A | eosin ×1.5 |
| `GAIN11` | decision; overshoot guard | exposure ×1.1, the highest `b_otsu` landing in Stage A (median 1.107, max 1.47). The ladder's harm rises above 1.0×, and L150's own failure mode is overshoot. |
| `GRAD` | decision; rank-altering (by rule) | thickness gradient 0.8 → 1.2 |
| `BLUR` | decision; rank-altering (by rule) | Gaussian σ = 1 px |
| `JPG` | decision; rank-altering (by rule) | JPEG quality 75, 4:2:0 |
| `H06_RECBOX`, `E15_RECBOX` | descriptive; recorded clicks only | H06 and E15 searched with each click's `UNPERT` box, for the like-for-like Stage A link (§2.6) and the box effect |

### 2.3 Arms (floors), per condition

`hem_p` is the condition's production `hematoxylin_od`, computed on 1024-row blocks.

| arm | floor on `hem_p` | role |
|---|---|---|
| `A1` | none | no clip |
| `A3` | `b_otsu(hem_p)` | the adopted floor |
| `L150` | `1.5 × b_otsu(hem_p)` | candidate |
| `W128` | `window_floor(hem_p, click, round(128/mpp) \| 1)` | candidate |

### 2.4 The template box

**Re-refined per condition.** `refine_template_box` (verbatim) runs on the condition's unclipped `hem_p` at the click,
as the product would, and all four arms of a click share that box. So within-condition contrasts are not confounded by
box changes. The perturbed-vs-original contrasts include them; `*_RECBOX` measures how much.

**On refusal:**
- The click uses a 51 px box at the rounded click, flagged. This approximates production's raw 51 px fallback.
- Stage A saw 1 refusal under H06 and 3 under JPG, out of 294.
- **Sensitivity:** flagged clicks are excluded *paired*, i.e. from that condition and from `UNPERT`, using
  `pair_stats_clicks`.

### 2.5 Metrics

Per click and arm:
- **The co-primaries:** TP@10 and TP@20 (as P@K), and the pool's mitoses (≤ 100 post-NMS detections).
- **Per-click harm vs no clip:** TP@10(arm) < TP@10(A1), same click and condition.
- **Per-click severe loss vs A3:** TP@10(f) ≤ TP@10(A3) − 2, and the mirror-image severe gain.
- **Descriptive:**
  - TP@30;
  - the worst click per ROI;
  - zero-TP@10 clicks;
  - top-10 FPs, and the white-adjacent ones: 61 px white share > 0.20, measured on the **original** RGB;
  - template SD, mass weight and share at the floor;
  - self-matches.

### 2.6 Statistics and decision rule (fixed before the smoke run)

**Notation.**
- `X_c(arm)` is the per-ROI click sum of a metric under condition `c`.
- `G_c(f) = X_c(f) − X_c(A1)`.
- `G_0 = G_UNPERT(A3)` is the validated gain, on this Mac and on all clicks.
- `n_r` is ROI r's click count.

**Estimates and two-sided CIs.** Every pooled estimate and two-sided 95 % CI (C1, PASS, the three-number report)
comes from per-ROI integer deltas through the verbatim `pair_stats_clicks` (tm_border_floor cell 11). That gives
click-weighted pooling and a ROI-cluster bootstrap (N_BOOT 10000, seed 20260917). Pool counts use k = 1.

**Every one-sided test** goes through the verbatim `sign_flip_null_counts` directly. That covers C2-weak, C2-broken,
C3, "better", "change", "worse" and the severe-loss check. The p-value is the null mass at or beyond the observed sum in
the stated direction. `pair_stats_clicks` is not used for these: its p is two-sided, and its pooled value and CI mean
nothing for a shifted vector.

The decision conditions are `p ∈ {H06, E15, GAIN11, GRAD, BLUR, JPG}`. The floors are `f ∈ {A3, L150, W128}`.

**C1, the clip still beats no clip.** It fails if the 95 % CI lower bound of `G_p(f)` is ≤ 0 on P@10, P@20 or pool.

**C2, the gain shrinks by at most a third: a level test against the validated gain.**
- `S_p(f) = G_p(f) − G_0`.
- For candidates this includes their unperturbed head start (L150 +0.27, W128 +0.24 pt at P@10).
- The margins m are a third of the validated gain: 1.1 pt (P@10), 0.9 pt (P@20) and 0.5 mitoses per click (pool).
- Every test is an exact one-sided sign-flip on the per-ROI integer deltas. **Each tier has its own Holm family per
  floor**, 6 conditions × 3 metrics = 18 tests, at 0.05. Two tests run per cell:
  - **C2-weak** (it opens the candidate comparison): the point estimate is below −m, and `S < 0` is Holm-significant.
    - At 12 clicks the significance part barely binds at P@10.
    - Its per-cell false-trigger rate is about 3 % at a true shrink of 0.41 pt, and about 20 % at 0.8 pt (within the
      margin).
    - That is acceptable only because a weak result cannot switch the floor by itself (BEATS is required).
  - **C2-broken** (the only C2 result that can conclude "rework"): `S < −m` is Holm-significant.
    - The per-ROI statistic is `S_r` (in TP counts) plus the margin over that ROI's clicks, scaled to an integer:
      - P@10: `100·S_r + 11·n_r`, since 1.1 pt × 10 × 100 = 11 per click;
      - P@20: `100·S_r + 18·n_r`, since 0.9 pt × 20 × 100 = 18 per click;
      - pool: `2·S_r + n_r`.
    - Scaling keeps the sign symmetry the sign-flip test needs. The convolution grows to about 10⁵ bins, which is still
      fast.
    - It fires rarely: at 12 clicks, at a true shrink of 1.5 pt, P@10 detects it in about 4 % of cells and pool in about
      29 %; at 2.0 pt, 35 % and 98 %.
    - The false rate is at most α at shrinks up to m.
- **PASS** if the two-sided 95 % CI lower bound is above −m. Otherwise the cell is **not shown robust**.

**C3, per-click harm does not grow (a paired test).**
- Per ROI, `D_r = 2·harm_p,r(f) − 3·harm_0,r(f)`. That is twice `harm_p − 1.5·harm_0`, kept integer so the verbatim
  integer statistics do not truncate it.
- It fails if the exact one-sided sign-flip test of `Σ D_r > 0` is Holm-significant at 0.05 within the floor's 6
  conditions.
- The 1.5× factor is §4.1's example; the user may change it. Each floor's own baseline is used (Intel: A3 20,
  L150 38, W128 21).

**Status of floor f under p:**
- **BROKEN** if C1 fails, C2-broken fires on any metric, or C3 fails.
- **WEAKENED** if it is not broken but C2-weak fires on any metric.
- **HOLDS** otherwise.

**Candidate contrasts (paired, within condition, `f − A3`).** Every test is an exact one-sided sign-flip on per-ROI
integers.
- **"Better" family:** `(f − A3)_p > 0` at P@10. Holm over 6 p × 2 candidates = 12 tests.
- **"Change" family:** the difference-in-differences `(f − A3)_p − (f − A3)_UNPERT > 0` at P@10. Holm over 12 tests.
  The per-ROI correlation ρ of `(f − A3)` between p and `UNPERT` is reported, because when ρ is low the
  difference-in-differences SE, SE·√(2(1−ρ)), can exceed the plain SE.
- **"Worse" family:** `(f − A3)_p < 0` on P@10, P@20 and pool. Holm over 36 tests.
- **f BEATS A3 under p** if it is Holm-significant in **both** the "better" and the "change" family. An unperturbed edge
  alone cannot make a candidate BEAT A3.
- A Holm-significant `UNPERT` contrast is **report-only**. At 12 clicks, W128's +0.24 pt unperturbed edge could reach
  significance, and the in-sample-switch rule still forbids changing the floor on unperturbed data.
- **Severe-loss check vs A3, per p.** Per ROI, `D_r` = (clicks where TP@10(f) ≤ TP@10(A3) − 2) − (clicks where
  TP@10(f) ≥ TP@10(A3) + 2). The exact one-sided sign-flip of `Σ D_r > 0` must not be Holm-significant over the floor's
  6 conditions.
  - It is defined against A3, not against no clip, because a switch is felt relative to the adopted floor.
  - Against no clip, "worse" grows with any distance from it even unperturbed: L150 has 38 such clicks against A3's
    20, yet against A3 directly it shows 4 severe losses and 3 severe gains unperturbed.
  - **Printed beside any L150 recommendation,** so §4.1's per-click criterion stays visible:
    - harm against no clip, for L150 and A3, in every condition;
    - the ≥ 1-TP churn against A3, worse and better clicks (unperturbed: L150 42 worse, 51 better; W128 9 worse,
      16 better).

**A3's status under p** is one of:
- **BROKEN**;
- **WEAKENED**;
- **BEATEN**: neither of the above, but some candidate BEATS A3;
- **HOLDS**.

**Recommendation rule** (to the user, not an automatic adoption):
1. If A3 HOLDS under every p: **keep 1.0× `b_otsu`.** List any cells not shown robust.
2. Otherwise, recommend each candidate f that meets all of:
   - f is BROKEN or WEAKENED under no p;
   - f BEATS A3 under every p where A3 is BROKEN, WEAKENED or BEATEN;
   - f is in no Holm-significant "worse" cell;
   - f passes the severe-loss check under every p.
   If both candidates qualify, report both; the user picks. The unmeasured multi-round concern about 1.5× goes beside
   any L150 recommendation.
3. **If A3 is BROKEN somewhere and no candidate qualifies:** report that floor estimation needs rework, e.g. stain
   normalisation before `b_otsu`.
4. **Otherwise** (A3 WEAKENED or BEATEN only, and no candidate qualifies): **keep 1.0×**, and report the cells where it
   is weakened or beaten.

**Stated power** (review; recorded data; no perturbed searches):
- **C2.** With up to 12 clicks, the half-width at P@10 is about 0.74 (1.01 at 6 clicks), so P(PASS | s = 0) ≈ 0.83.
- **Minimum detectable "better" effect** (one-sided, Holm rank 1 of 12): about 0.31 pt for W128 and 0.97 pt for L150.
  - Both figures assume the effects vary click by click. If they are ROI-specific (W128 − A3 is nonzero on only 16/49
    ROIs unperturbed), the 12-click gain is much smaller.
  - **At the edges Stage A predicts, L150 BEATS A3 is effectively unreachable.** This is a stated limitation, not a
    reason to loosen the rule.
- **The v2 rule simulated** (6 clicks; review c15): with no true shrink, A3 is kept in 97 % of runs. The v3 rule is stricter about switching, so this is a lower bound for v3.
- **C2-broken's level at the H0 boundary** (review c19; 1,000 draws): 0.047–0.054 at 0.05 and 0–0.007 at 0.05/18, for 6 and 12 clicks.

**Three-number report per p (§4.1), on P@10, P@20 and pool, with CIs:**
- clip vs no clip under p: `G_p(A3)`;
- clipped, perturbed vs original: `X_p(A3) − X_0(A3)`;
- no clip, perturbed vs original: `X_p(A1) − X_0(A1)`.
The same three are reported for L150 and W128.

**Stage A link (descriptive), on the recorded 6 clicks only:**
- On `H06_RECBOX` and `E15_RECBOX`, compare the pooled measured `−S` with Stage A's PL, for A3 and W128.
- Report the per-ROI correlation alongside its ceiling. The ceiling is about 0.08, because per-ROI PL has an SD of
  0.20–0.29 while per-ROI `S` has an SD of about 3.6 pts.

### 2.7 Harness gates and the action on each failure (all fixed before the smoke run)

| # | gate | on failure |
|---|---|---|
| 1 | Production constants unchanged; skimage 0.24.0; multi-Otsu nbins 256; every worker's thread pools at 1 | halt; fix the code; re-run from scratch |
| 2 | Verbatim copies asserted against their source cells. From `tm_border_floor` cells 6 and 8: the search helpers (`fill_outside` … `white_frac`). From `tm_floor_definition` cells 9 and 11: `window_floor`, `search_arm`, `tp_rows`. From `tm_border_floor` cells 10, 11 and 21: the statistics (incl. `pair_stats_clicks`, `sign_flip_null_counts`). And `refine_template_box` (SHA-1 `21c43a20…`). The cell-11 `parity()` is **not** reused, because it asserts template-SHA and threshold equality, and those differ across machines. | halt; fix; re-run |
| 3 | **Parity on outcomes**, recorded 6 clicks, `UNPERT`, against the Intel lists (A1/A3 from `tm_border_floor_49roi_6seed_search.pkl`; L150/W128 from `tm_floor_definition_49roi_6seed_pool.csv.gz`). Each list's TP@10/20/30 must be identical and its pool count within ±1. A mismatch is allowed only if a non-verbatim helper `trace_list` explains it as a near-tie; `trace_list` re-runs `fused_map` + `detect` and keeps the fused map and the pre-NMS peaks. Three kinds count: (i) a swap at a K boundary, where this Mac's scores of the two swapped peaks are within 2e-5 relative; (ii) a 1-px local-maximum shift, where an Intel-only coordinate lies within 1 px of a Mac coordinate with scores within 2e-5 relative (observed on 343 s0 L150, rank 93); (iii) a tail or NMS near-tie, where this Mac's fused map at the Intel-only coordinate is within 2e-5 relative of the competing peak at the 100-peak cap or in NMS. Reported: list-identity rate per arm ((cx, cy, bucket) and |Δscore| ≤ 2e-5 relative). Review: 48/48 lists identical on 529 and 343, max relative drift 3.8e-6. | halt; investigate; no rule change |
| 4 | Stage A coupling, recorded 6 clicks: every condition's floors (`b_otsu`, W128) and re-refined boxes equal Stage A's CSVs exactly (read `round_trip`) | halt; fix |
| 5 | `NULL2`, click 0 of each ROI. The boxes are re-refined on `2·hem` and must be unchanged. Floors and fill exactly ×2. Lists identical in (cx, cy, bucket, order). Scores, threshold, median and MAD exactly ×4. od ×2, with `od_padded` built from `2·hem`. The band sentinel is unchanged. | halt; fix |
| 6 | Clip identity: `clip_floor`'s own assertions; blank fill equals `b` | halt; fix |
| 7 | The statistics reproduce the recorded Intel `A3 − A1` rows at K = 10/20/30, and `pair_stats_clicks` equals `pair_stats` at 6 clicks | halt; fix |
| 8 | Smoke determinism: N = 3 equals N = 1 after a canonical sort | halt; fix |
| 9 | **`JPG95` calibration** (phase 1). Families are defined within `JPG95` alone. No floor may be BROKEN or WEAKENED: C2 Holm over its 3 metrics per floor, C3 one test per floor. No candidate may show a Holm-significant change-from-baseline: the difference-in-differences test, 2 tests. No "worse" cell either: 6 tests. | **stop before phase 2** and report the numbers to the user. The rule is not changed unilaterally. |
| 10 | Extra-click draw: the joint gate reproduces the recorded `n_joint_valid` on 49/49 ROIs, and all 294 recorded seeds are eligible | halt; fix |

The `JPG95` C2 half-width at P@10 is reported as the measured noise floor beside every C2 result.

## 3. Implementation

### 3.1 Files

None of these are committed until the user says so.

- `production_hematoxylin_only/multiotsu_pad_arms/tm_floor_stress_stageB_49roi.ipynb`
  - Cell 0 is the pre-registration (§1–§2); the reading goes in the last cell.
  - Figures are embedded only.
- `production_hematoxylin_only/multiotsu_pad_arms/tm_floor_stress_stageB_worker.py`
  - It holds the verbatim helpers, gated by the notebook, and the click-draw function.
  - It imports the Stage A worker for synthesis.
- **Outputs**, in `results/tm_floor_stress_stageB_49roi_*`:
  - `clicks.csv`: the click set and the extra draw's record;
  - `per_run.csv`: one row per condition × arm × click;
  - `top30.csv.gz`;
  - `unpert_pool.csv.gz`, for parity;
  - `contrasts.csv`: every contrast, including the three-number report;
  - `verdicts.csv`, `harm.csv`, `parity.csv`, `stageA_link.csv`, `checks.csv`, `decision.csv`, `adversarial_log.csv`;
  - `_smoke` copies.

### 3.2 Parallelism and memory

- **Pool:** a spawn pool with one task per (ROI, condition), largest first by pixels × searches. The draw of the extra
  clicks (one `hem` and refinement per ROI) is its own pool step before phase 1.
- **Per task:** load the image, synthesise, compute `hem_p`, the original-RGB white integral and the boxes. Then, for
  each distinct floor: build the clipped channel, run that floor's searches, and **free it** before the next. W128 has
  one floor per click.
- **Memory (review, 529.tiff, build-use-free):** 2.99 GB per UNPERT task and 3.16 GB with `NULL2`. Holding every clipped
  channel would need 6.6–8.6 GB.
- **Workers:** `N = min(16, floor(36 GB / smoke peak RSS))`, about 10–11.
- **Threads:** pinned as in Stage A (env vars before numpy; `threadpool_limits(1)`; `cv2.setNumThreads(1)`), and
  asserted.
- **GPU: declined.** Gates 3 and 5 need `cv2`'s CPU `matchTemplate`; torch and mlx are not installed.

### 3.3 Phases and runtime

The decision rule, the margins and the halt actions are frozen in this plan **before the smoke run**.

| step | contents | searches (≈) | wall (≈, 12 searches/s) |
|---|---|---|---|
| extra-click draw | per ROI: `hem` plus refinement over the pool | 0 | 1–2 min |
| smoke | phase-1 conditions only (`UNPERT`, `JPG95`, `NULL2`) on 529, 128 and 401, at N = 3 and N = 1; sets N; gate 8. **No decision condition is searched in the smoke run.** | ~400 ×2 | 5–8 min |
| phase 1 | `UNPERT` + `JPG95` (545 clicks × 4 arms each) + `NULL2` (49 × 4); gates 3, 5, 7 and 9 | ~4,550 | ~7 min |
| phase 2 | H06, E15, GAIN11, GRAD, BLUR, JPG (545 × 4 each) + `H06_RECBOX` and `E15_RECBOX` (294 × 4 each); gate 4 runs inside each task, before its searches | ~15,430 | ~22 min |

**Total:** about 20,000 searches, roughly 35–40 min of wall time.

**Stage A caveat carried forward.** H06 cannot lighten copied zero-channel pixels: 31 ROIs have some, at most 0.92 % of
pixels. The per-ROI count is reported beside the H06 results.

## 4. Adversarial concerns (registered in cell 1)

| concern | test |
|---|---|
| cross-machine drift | gate 3 checks outcomes with near-tie tracing; every contrast uses this Mac's own `UNPERT` |
| a FAIL or a recommendation fires on noise | "rework" needs C2-broken (significant beyond the margin) or C1/C3 failures; C2-weak only opens the comparison; BEATS needs Holm significance in both the "better" and "change" families; C3 is a paired integer test with Holm |
| C2 is underpowered | up to 12 clicks; a cell without PASS is reported as "not shown robust", and WEAKENED alone keeps the status quo; `JPG95` noise floor; power stated |
| L150's BEATS is unreachable | stated limitation (MDE about 0.97 pt at 12 clicks, if effects vary click by click) |
| L150 overshoot is never tested | `GAIN11` is a decision condition |
| box changes confound perturbed vs original | within-condition contrasts share the box; `H06_RECBOX` and `E15_RECBOX` |
| refusals | paired exclusion sensitivity with `pair_stats_clicks` |
| the rule is tuned on decision data | the rule is frozen before the smoke run; neither the smoke nor phase 1 searches a decision condition; every halt has a pre-registered action |
| the extra clicks differ from the recorded ones | the recorded draw's own tier and joint gate (gate 10), plus refinement acceptance; parity and Stage A coupling stay on the recorded 6 |
| research path vs production | stated in §1 |

## 5. Review log

### Round 1 (fresh reviewer agent, 2026-10-09; scripts in the session scratchpad `review_stageB/`)

No [BLOCKS] findings. The harness checks passed:
- **UNPERT parity:** 48/48 lists identical on 529 and 343.
- **`NULL2`:** exact.
- **Memory:** 3.0–3.2 GB per worker with build-use-free.
- **Throughput:** 12 searches/s.

| # | finding | verdict | action |
|---|---|---|---|
| 1 | FAIL fires on a point estimate, unioned over 15 cells; the recommendation needs only "≥ A3 at the point estimate" | **accepted** | FAIL needs the margin plus Holm significance; a candidate must BEAT A3 (Holm-significant) wherever A3 fails; INCONCLUSIVE is reported as "not shown robust" |
| 2 | C3 is noisy, and own-baseline scaling lets L150 carry about 2× A3's harm | **accepted, with a different candidate check** | C3 becomes a paired bootstrap test with Holm. The per-click check for candidates is severe losses **vs A3**, not harm vs no clip, which grows with any distance from no clip (L150 against A3 directly: 4 severe losses, 3 severe gains unperturbed). Open for round 2. |
| 3 | L150 contrast unpowered; "statistics hard-code 6 clicks" is wrong (`pair_stats_clicks` exists); 12 clicks feasible | **accepted** | up to 12 clicks, same unanimous tier, via `pair_stats_clicks`; L150's unreachable BEATS stated. The DiD candidate contrast is used only as the BEATS side condition, not as the test. |
| 4 | Overshoot (L150's failure mode) is never searched | **accepted** | `GAIN11` added as a decision condition |
| 5 | The smoke run searches decision conditions | **accepted** | the smoke runs phase-1 conditions only; the rule is frozen before the smoke |
| 6 | No pre-registered action for halts; `JPG95` candidate family unspecified | **accepted** | gate table with actions; `JPG95`'s own Holm families (2 / 6) |
| 7 | The refusal sensitivity uses the wrong denominator | **accepted** | paired exclusion with `pair_stats_clicks` |
| 8 | Gate 3: do not reuse `parity()`; use a relative tolerance; untraceable tail and NMS cases | **accepted** | outcome parity, 2e-5 relative, Mac-side fused-map tracing; `tp_rows` added to gate 2 |
| 9 | `NULL2` expectations must be stated | **accepted** | stated in gate 5; boxes re-refined on `2·hem` |
| 10 | The Stage A link cannot be answered as designed | **accepted** | `H06_RECBOX` and `E15_RECBOX` replace `JPG_RECBOX`; pooled comparison; per-ROI correlation ceiling reported |
| 11 | C1 is tested on P@10 only | **accepted** | C1 on all three co-primaries |
| 12 | One Holm family of 30 dilutes the "better" tests; the inherited unperturbed edge can read as BEATEN | **accepted** | separate "better" (12) and "worse" (36) families; BEATS also needs a positive change from baseline; the `UNPERT` contrast is reported |
| 13 | §4.1's perturbed-vs-original numbers are missing | **accepted** | the three-number report |
| 14 | Memory: free each clipped channel after use | **accepted** | build-use-free |
| 15 | The smoke run takes about 8–12 min at N = 1 | **accepted** | the smoke is now phase-1 conditions only; re-estimated |
| 16 | This tests the research path, not production | **accepted** | stated in §1 |
| 17 | The shared G₀ gives candidates a head start | **accepted** | C2 is labelled a level test |

### Round 2 (same reviewer, plan v2)

One [BLOCKS] and three [CHANGES-DECISION] findings. All four are fixed in v3, and the reviewer judged that the plan
converges after them.

| # | finding | verdict | action |
|---|---|---|---|
| 1 | [BLOCKS] C3's `D_r = harm_p − 1.5·harm_0` is half-integer, and `pair_stats_clicks` casts it to int64, silently turning A3's 1.5× allowance into about 1.15× | **accepted** | `D_r = 2·harm_p − 3·harm_0`; exact one-sided sign-flip |
| 2 | The BEATS side condition (positive change from baseline in point estimate) filters nothing: "better" equals BEATS within 0.001 in simulation | **accepted** | BEATS needs Holm significance in both the "better" and the difference-in-differences "change" family; gate 9 uses the change test; ρ is reported; a significant unperturbed W128 edge is report-only |
| 3 | FAIL's text ("significantly beyond a third") and its test (significant against 0, plus the point estimate) disagree; at 12 clicks a within-margin 0.8 pt shrink reaches "rework" about 35 % of the time | **accepted, option (b)** | two tiers. C2-weak (the old test) only opens the comparison. C2-broken, a significant `S < −m` with an integer-shifted statistic, is the only C2 route to "rework". §1 text updated. |
| 4 | Gate 9's C2 and C3 families are undefined for `JPG95` | **accepted** | defined within `JPG95` (3 metrics per floor; 1 C3 test; 2 change tests; 6 worse tests) |
| 5 | The severe-loss bootstrap degenerates at zero baselines (W128); §4.1's harm criterion should stay visible | **accepted** | exact one-sided sign-flip; harm and ≥ 1-TP churn printed beside any L150 recommendation |
| 6 | The extra-click gate differs from the recorded joint gate on 18 % of eligible clicks | **accepted** | joint gate plus refinement acceptance; gate 10; 545 clicks; estimand stated as click-weighted |
| 7 | Gate 3 wording: K-boundary swaps have no Intel-only coordinate; dilation ties are missing; tracing needs a non-verbatim helper | **accepted** | three tracing kinds; `trace_list` named |
| 8 | Stale MDEs and click count | **accepted** | 0.31 / 0.97 pt with the click-varying caveat; 545 |

### Round 3 (same reviewer, plan v3): CONVERGED

No [BLOCKS] or [CHANGES-DECISION] findings.

**Verified:**
- **C2-broken's integer statistics.**
  - Each sums to zero at exactly m; the sign symmetry holds with 6–12 clicks.
  - Measured level at the boundary: 0.047–0.054 at 0.05.
  - The P@10 and P@20 deltas are never zero.
  - The null arrays stay under 38,000 bins, about 0.5 ms per test.
- **Rules 1–4.** None switches on noise. Every switch needs BEATS, so a false switch is bounded by the change family's
  error rate (5 %). "Rework" on noise is bounded by α per floor.
- **BEATS and gate 9** are correct as written.
- **The extra-click rule** matches the reproduction: 49/49 ROIs, all 294 seeds eligible, 545 clicks.

Six [MINOR] wording items, all fixed in v3 final:
- the stale 1.17 pt in §4;
- §1's "only C2-broken concludes rework" (C1 and C3 also do);
- per-tier Holm families stated;
- one-sided tests carved out of `pair_stats_clicks`;
- the orphaned bootstrap sentence and the leftover "INCONCLUSIVE" removed;
- the 97 % figure labelled as a v2 simulation.
