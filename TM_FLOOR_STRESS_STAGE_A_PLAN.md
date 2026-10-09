# Floor-clip stress test, Stage A: where the floor lands under stain and scanner change (plan)

**Status:** plan v3 (final), 2026-10-08. The reviewer judged it converged after three rounds (see §6).
**Run 2026-10-08:** `production_hematoxylin_only/multiotsu_pad_arms/tm_floor_stress_stageA_49roi.ipynb`, with outputs in
`results/tm_floor_stress_stageA_49roi_*.csv`. All gates held. Three cells fired: H06 × `b_otsu` (PL 0.41 pt), H06 × W128
(0.33) and E15 × `b_otsu` (0.28). By the rule, Stage B would run on H06, E15, GRAD, BLUR, JPG and the null control.
**Stage B was planned** (`TM_FLOOR_STRESS_STAGE_B_PLAN.md`) **and deferred to a pre-ship check on 2026-10-09.**
**Source spec:** `FUTURE_IDEAS.md` §4.1, Stage A. This plan makes the spec executable and departs from it where §6 says
why.

---

## 1. Purpose

The floor clip `hem_c = max(hem, b)` with `b = b_otsu` (the lower 3-class multi-Otsu threshold of the ROI's
`hematoxylin_od`) is the adopted change (D3 amendment, 0401ae3).

The floor-level ladder (`tm_floor_definition_49roi_6seed.ipynb`, c × `b_otsu`) gives P@10 54.6 / 56.5 / 57.9 / 58.2 /
57.9 at c = 0 / 0.5 / 1.0 / 1.5 / 2.0.
- 0.5× is a significant loss against 1.0×.
- 1.5× is non-inferior on P@10, P@20 and pool.
- 2.0× is not non-inferior on P@10.
- W128, the same threshold over a 128 µm click-centred window, landed at 0.89–1.23× (5th–95th percentile) and was
  non-inferior.

The floor is estimated per image. On another lab's slides it lands wherever that slide's stain and scanner put it.

**Stage A asks:** under realistic stain and scanner changes, do the two decision-set estimators still land where the
searched ladder predicts no meaningful loss? It runs **no template searches**.

**It decides one thing: whether Stage B runs, and on which perturbations.** Stage B means about 8,200 searches (about
15–20 min of compute) plus a new parallel search runner.

**What Stage A cannot answer:**
- **precision or pool changes** (Stage B);
- **whether the knee moves on a perturbed image** (A1);
- **the effect of rank-altering perturbations (blur, JPEG, gradient) on the clip's gain.** Their level is checked here.
  Their reshuffle of pixel ranks is tested only if Stage B runs.
- **real cross-scanner data** (Stage C);
- **held-out ROIs;**
- **later loop rounds;**
- **chromatin od.**

## 2. Assumptions

- **A1. Equal tissue share clipped means the same position relative to the knee.**
  - A floor `b'` on the perturbed channel clips a share `p'` of the original tissue pixels. Its equivalent level
    `c_eq` is the multiple of the original `b_otsu` that clips the same share on the original channel.
  - This makes the *clipped set* comparable exactly for transforms that are monotone in `hem`. It does not make the
    TM response equal.
  - For rank-altering transforms it is a proxy: blur, JPEG, eosin leak, and the RGB tone curves (which are not
    monotone in `hem`).
  - Stage A cannot test A1. Stage B can.
  - A tissue-quantile floor would be invariant under `c_eq` by construction. That is the construct's blind spot, stated
    here rather than measured.
- **A2. The synthesised perturbations stand in for real variation at "realistic" strength.**
  - The strengths are a judgment call, inherited from §4.1 where it gives them.
  - They are fixed in this plan and not changed after any result.
- **A3. The per-ROI Macenko basis is a usable model of the slide's stains.**
  - The residual is kept.
  - Negative concentrations are model error, so they are not scaled.
  - Pixels with any channel at 0 lie outside the synthesis model's range: it caps OD at 2.41, while `rgb2hed` maps
    0 to 6. Those pixels are left unchanged by synthesis, and their count is reported.
  - Caveat: Macenko's E vector is the pinkest 1 % extreme, so E×1.5 likely overstates eosin leak into H.
- **A4. The tissue mask comes from the original image** (gray ≤ 215) as the `c_eq` denominator. Estimators that use a
  tissue mask (TOTSU, TW128) compute it on the perturbed image, because a shipped estimator would see that image.

## 3. Design

### 3.1 Units

- **ROIs:** the 49 ROIs of `tm_floor_definition` (14 `original_14` + 35 `testing_35`).
- **Clicks:** the recorded 294 clicks, from the `A0` rows of `results/tm_border_floor_49roi_6seed_per_run.csv`
  (`rank_key == 'tm_score'`): `click_cx`, `click_cy`, `tcx`, `tcy`, `bs`.
- **`roi_index`** is the ROI's position in the sorted 49-file list, including in smoke runs.

### 3.2 Image variants

`hem'` is always computed with the **production** `midog_utils.channels.to_channel(rgb', 'hematoxylin_od')`
(skimage 0.24.0 `rgb2hed`, clamp at 0).

**Synthesis model** (stain perturbations only):
- `OD = −log10((I + 1)/256)` in float64, and `I = clip(rint(256·10^(−OD) − 1), 0, 255)`.
- **Basis:** Macenko 2009, Algorithm 1, as in `floor_robustness/macenko_vectors.py`. That is β = 0.15 (log10 OD),
  α = 1 %, 500,000 sampled pixels, SVD through the origin, sign fixes, and H = the vector with the larger red component.
  The RNG is `np.random.default_rng([20261008, roi_index])`.
- **Concentrations:** `M = [h e]` (3×2) and `C = OD @ pinv(M).T`.
- **Perturbation:** `OD' = OD + (kH − 1)·max(cH, 0)·h + (kE − 1)·max(cE, 0)·e`. Equivalent to rescaling the positive
  part of the concentrations and keeping everything else.
- **Pixels with any channel = 0** are copied unchanged.
- Synthesis runs in 256-row blocks.
- The fully reconstructed form `(C·k) @ M.T + R` is used only for the `ID` gate.

| id | class | definition |
|---|---|---|
| `ID` | gate | the synthesis with k = (1, 1), reconstructed through `C @ M.T + R` → RGB must be bit-identical (encoding path only) |
| `NULL2` | gate | `hem' = 2·hem` at channel level; TOTSU and TW128 use the original gray |
| `GAMMA` | near-null calibration | equal-channel RGB gamma, `I' = rint(255·(I/255)^1.1)`. Under `rgb2hed` this scales OD by 1.1 exactly, apart from 8-bit rounding and the I = 0 floor. |
| `JPG95` | near-null calibration | JPEG quality 95, 4:2:0 (cv2, BGR) |
| `H06`, `H14` | realistic | kH = 0.6, 1.4 (kE = 1) |
| `E15` | realistic | kE = 1.5 (kH = 1) |
| `CDN` | realistic | RGB LUT `rint(12 + I·231/255)` (haze or flare) |
| `CUP` | realistic | RGB LUT `clip(rint((I − 12)·255/231), 0, 255)` (crush) |
| `GAIN09`, `GAIN11` | realistic | all channels `clip(rint(g·I), 0, 255)`, g = 0.9, 1.1 (exposure) |
| `WBR095`, `WBR105` | realistic | red channel only, gain 0.95, 1.05 (white balance) |
| `WBB095`, `WBB105` | realistic | blue channel only, gain 0.95, 1.05 |
| `BLUR` | realistic | `cv2.GaussianBlur(rgb, (0, 0), 1.0)`, `BORDER_REFLECT_101` |
| `JPG` | realistic | JPEG quality 75, 4:2:0 |
| `GRAD` | realistic | `OD' = g(x)·OD`, `g` linear 0.8 → 1.2 along the long axis (pale at index 0); zero-channel pixels unchanged |
| `SHIFT` | realistic, W128 only | original image; window moved per §3.3 |

Why GAIN and WB are included: the review measured that `rgb2hed`'s clamp breaks the "predicted null" for them (128.tiff,
22 % of tissue clamped).
- GAIN0.9 gives `b_otsu` a `c_eq` of 0.82.
- WB R×1.05 gives 1.36.

### 3.3 Floor estimators

Every estimator is the lower threshold of `threshold_multiotsu(·, classes=3)` (nbins 256), as float32.

| id | definition | role |
|---|---|---|
| `b_otsu` | over the whole `hem'` (`roi_floors`, verbatim) | **decision set:** the adopted estimator. The product keeps the default ROI size (user, 2026-10-08). |
| `W128` | over the click-centred `round(128/mpp) \| 1` square, truncated at the ROI edge (`window_floor`, verbatim) | **decision set:** product-feasible and search-validated |
| `TOTSU` | over `hem'[gray' ≤ 215]` | descriptive: a fallback candidate, never searched |
| `TW128` | the W128 window, tissue pixels only. NaN if fewer than 5 % of the window's pixels are tissue; NaNs are excluded from medians and counted | descriptive: a fallback candidate, never searched |

**`SHIFT`** (W128 only, original image):
- The offset is `d = round(64/mpp)` px, in eight directions: E, W, S, N, SE, SW, NE, NW, i.e. (±d, 0), (0, ±d),
  (±d, ±d).
- Each offset window is the W128 square centred at the shifted point, with the centre clipped into the ROI and the
  window truncated at the edge.
- Its glass share is `mean(gray0 > 215)` over that truncated window.
- The offset with the highest glass share is used. Ties go to the first in that order.

### 3.4 Metric

There is one primary metric. In the notation below:
- `T0` is the original tissue mask;
- `N_t = |T0|`;
- `b0` is the original ROI `b_otsu`;
- `v0 = sort(hem0[T0])`.

The metric is the equivalent level:

```
p'   = count(hem'[T0] < float32(b')) / N_t
c_eq = v0[min(round(p'·N_t), N_t − 1)] / b0
```

**Per-ROI value:**
- For `b_otsu` and `TOTSU`, the single value.
- For `W128` and `TW128`, the median over the ROI's 6 clicks.

**Also recorded,** because they matter for Stage B or for sanity:
- **Refinement box.** `refine_template_box` (verbatim, SHA-1 `21c43a20…`) on `hem'` at each click: is it accepted, and
  do `bs`, `tcx` and `tcy` equal the recorded box?
- **Synthesis checks.** The count of zero-channel pixels, `‖R‖/‖OD‖` on tissue, and the angles of `h` and `e` against
  `floor_robustness/macenko.csv`.
- **Bin width.** The multi-Otsu bin width over `b'` for `b_otsu`.
- **Spread.** The 5th–95th percentile of per-ROI `c_eq`, alongside the median.

### 3.5 Decision rule, fixed before the run

The measured ladder is used as a dose–response curve, and each landing is turned into the P@10 it is predicted to
cost. This uses the whole landing distribution. A median or a band edge cannot, which is why v2's rule was replaced
(§6, round 2).

**Ladder** (`results/tm_floor_definition_49roi_6seed_ladder.csv`, P@10 at `tm_score`):

| c | 0 | 0.25 | 0.5 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|---|---|
| P@10 | 0.545918 | 0.554422 | 0.564966 | 0.578912 | 0.581633 | 0.578912 |

**Definitions:**
- `L(c)` is piecewise-linear interpolation through the measured points from 0 to 1.0.
- The **loss of one landing** is `loss(c) = max(0, L(1.0) − L(c))` for c ≤ 1.0, and 0 for c > 1.0. Up to 2.0 the
  ladder is flat or higher in mean; above 2.0, OVER handles the ROI level.
- Landings above 1.0 are never credited as gains, because the 1.0 → 1.5 step of +0.27 pt is not significant.
- The linear 0.5 → 1.0 segment assumes a knee shape that was never measured.
  - The adjacent steps are 3.4, 4.2 and 2.8 pt per unit c, so the curve is close to linear.
  - A curve that saturates before 1.0 would make this estimate conservative, i.e. it would predict too much loss.

**Rules.** For each realistic perturbation × decision-set estimator (`b_otsu` or `W128`):
- **PL, the predicted loss.** `PL` is the mean of `loss(c_eq)` over the 294 click-landings. A `b_otsu` value counts
  once for each of its ROI's 6 clicks; `W128` uses each click's own floor.
  - **PL fires at ≥ 0.25 pt P@10.**
  - That is a quarter of the ladder's 1.0-pt non-inferiority margin, and the loss of a uniform landing at about 0.91.
- **OVER, the extrapolation guard.** It fires if the median per-ROI `c_eq` is > 1.5, or if any ROI is > 2.0.
  - Above 1.5 only L200 was searched, and it was not non-inferior at P@10.
  - Above 2.0 nothing was searched.

**Outcome:**
- **If any realistic perturbation fires,** Stage B runs on those perturbations, plus `GRAD`, `BLUR`, `JPG` and the
  null control (§4.1).
- **If none fires,** stop.
  - Stage B is not run.
  - The §1 limitation about rank-altering perturbations is reported as untested.
- There is no grey zone.
- Adding one perturbation to Stage B costs about 1,176 searches, about 84 s at the benchmarked 14 searches/s, so a
  sensitive screen costs little. Stage B's per-worker memory has not been measured yet.

**Calibration gate.**
- `NULL2`'s PL equals the unperturbed PL exactly. That is 0 for `b_otsu`, and 0.056 pt for W128 on all 49 ROIs, from `floors.csv`.
- `GAMMA` and `JPG95` must fire neither rule. If either does, the rule is miscalibrated, and the run halts with a
  report.
- Unperturbed W128's PL is also reported. The W128 search measured +0.24 pt [−0.10, +0.58] against A3, so its PL
  should be small (≲ 0.1).
  - This is the only available check of the screen against a search, and it is descriptive.

**Reported beside the decision, never used for it:**
- the median and the 5th–95th percentile of per-ROI `c_eq`;
- the share of clicks below 0.9;
- the ROIs that contribute most to PL, with their `b_otsu` bin width in `c_eq` units. A ROI whose landing is within
  one `b_otsu` bin of 1.0 is flagged "quantisation-limited", but it still counts;
- TOTSU and TW128;
- per-domain splits;
- box changes.

**Product condition on `b_otsu`: resolved.** The user decided on 2026-10-08 that the product keeps the default ROI size.
So `b_otsu` over the ROI, which is what Stage A measures, is the shippable estimator, and no frame-size window is
needed.

### 3.6 Pilot observations (seen during review; NOT predictions)

During review, the reviewer ran perturbations on hand-picked extreme ROIs: 5 in round 1, 7 in rounds 2–3. These are not 49-ROI estimates. The
rule in §3.5 is anchored to the measured ladder and the NI margin, not to these numbers.

- **`b_otsu` median `c_eq`:** H06 0.864, E15 0.877. On 401, H06 gives 0.671 for `b_otsu` while W128 holds at 1.048.
- **`W128` median `c_eq`:** H06 0.890, E15 0.921.
- **E15 lowers `c_eq`** even though `b'/b0` > 1, because rising stroma reduces the tissue share clipped.
- **SHIFT raised W128** on 401 and 128 rather than lowering it.
- **Unperturbed TW128** sits 0.1–0.3 above W128.
- **Implication:** with this rule, Stage B is likely to run for `b_otsu` on H06 and E15.

**Round-2 pilot** (v2 synthesis; 7 ROIs: 201, 343, 529, 128, 401, 211, 403)
- **H06, `b_otsu` per ROI:** 0.685, 0.829, 0.950, 0.875, 0.768, 0.783, 0.964.
  - Under §3.5 that is a PL of about 0.46 pt on these 7 ROIs, so it fires.
- **E15, `b_otsu` per ROI:** 0.790, 0.890, 0.964, 0.939, 0.850, 0.751, 1.008. That is a PL of about 0.32 pt, so it
  fires.
- **W128** stays within 0.82–1.47 per ROI.
- **Stage B is therefore near-certain.** What Stage A still decides is which marginal variants join it, and by how
  much each estimator misses.

## 4. Implementation

### 4.1 Files

None of these are committed until the user says so.

- **Notebook:** `production_hematoxylin_only/multiotsu_pad_arms/tm_floor_stress_stageA_49roi.ipynb`.
  - Cell 0 is the pre-registration (§1–§3.5).
  - The reading goes in the last cell.
  - Figures are embedded, with no PNG files.
- **Worker module:** `production_hematoxylin_only/multiotsu_pad_arms/tm_floor_stress_stageA_worker.py`.
  - It is needed because spawn workers cannot unpickle functions defined in a notebook cell.
  - Its sha256 is recorded as an input.
- **Outputs:** `results/tm_floor_stress_stageA_49roi_{clicks,roi,summary,boxes,checks,decision,adversarial_log}.csv`,
  plus `_smoke` copies. Every frame is sorted canonically before it is written.

### 4.2 Parallelism

- **Pool.** `multiprocessing.get_context('spawn').Pool(N, initializer=init_worker)`. One task per ROI, largest ROI
  first, through `imap_unordered(chunksize=1)`.
- **Threads.**
  - Set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and `VECLIB_MAXIMUM_THREADS` to 1 at the top
    of the worker module, before `import numpy`. Set them in `os.environ` before the pool is created too.
  - `init_worker` calls `threadpoolctl.threadpool_limits(1)` (threadpoolctl 3.7.0 is installed) and
    `cv2.setNumThreads(1)`. It then asserts that every pool in `threadpool_info()` reports `num_threads == 1`.
  - Do not gate on `cv2.getNumThreads()`: it reports 18 under the GCD backend even after the call.
  - Unpinned threads made one synthesis take ~950 s in review, against 1.44 s pinned.
- **Memory and worker count.**
  - Every `hem` is computed with the production function on 1024-row blocks.
  - Review round 2 found this bit-identical to the whole-image call on 529. It held for the original image at 1024-,
    777- and 1-row blocks, and for H06, JPG and BLUR. It is also slightly faster.
  - It cuts peak RSS per worker from 5.1 GB to 2.46 GB on the 39 MP 529.tiff.
  - A bit-identity assert stays on the smoke ROIs, for the original and one perturbed image.
  - The smoke set includes 529.tiff.
  - `N = min(16, floor(36 GB / smoke peak RSS))`, which gives ≈ 14.
- **Determinism.** Smoke with N = 1 and with N = 3 must give identical outputs after the canonical sort.
- **GPU: declined.**
  - The floor must come from the production CPU `rgb2hed` channel.
  - Synthesis is 1.5–1.95 s pinned, used by 4 variants plus `ID`. That is about 30 % of a worker's ~22–27 s per ROI,
    so a GPU would save under 1 min.
  - Neither torch nor mlx is installed, and either would need a separate env (numpy is pinned at 1.26.4).
  - A GPU could matter in Stage B's correlation.

### 4.3 Per-ROI worker (sketch)

```python
def run_roi(task):
    rgb0 = ds.load_roi(task.path)
    hem0 = ch.to_channel(rgb0, 'hematoxylin_od')
    gray0 = cv2.cvtColor(rgb0, cv2.COLOR_RGB2GRAY)
    T0 = gray0 <= 215
    b0 = roi_floors(hem0)['b_otsu']
    v0 = np.sort(hem0[T0])
    M = macenko_basis(rgb0, np.random.default_rng([20261008, task.roi_index]))
    rows, boxes, checks = [], [], [self_consistency(hem0, T0, v0, b0), synthesis_checks(rgb0, M)]
    for var in VARIANTS:
        hem1, gray1 = make_variant(var, rgb0, hem0, gray0, M)
        ht = hem1[T0]
        for est, seed_index, b in floors(var, hem1, gray1, gray0, task.clicks, task.mpp):
            rows.append(dict(var=var, est=est, seed_index=seed_index, b=b, c_eq=c_eq(ht, v0, b0, b)))
        boxes += [box_row(var, refine_template_box(hem1, r.click_cx, r.click_cy, task.mpp), r) for r in task.clicks]
    return task.file_name, rows, boxes, checks, peak_rss()


def c_eq(ht, v0, b0, b):
    p = np.count_nonzero(ht < np.float32(b)) / ht.size
    return float(v0[min(int(round(p * ht.size)), ht.size - 1)]) / b0
```

### 4.4 Gates (halt on failure)

1. **Configuration.**
   - `production` constants are unchanged; skimage is 0.24.0; `threshold_multiotsu` nbins is 256.
   - `hem0.min() ≥ 0`, so the clamp is present.
   - Every worker reports `threadpool_info` at 1 thread.
2. **Verbatim copies.**
   - `roi_floors`, from `tm_border_floor` cell 8.
   - `window_floor`, from `tm_floor_definition` cell 9.
   - `refine_template_box`, from `tm_border_floor` cell 4, matched by SHA-1.
3. **Reproduction on this Mac** (CSV read with `float_precision='round_trip'`; compared as float32).
   - `b0` and unperturbed W128 match the recorded values within one multi-Otsu bin. Report the float32-equal count and
     the max relative difference.
   - The review saw `b0` equal on 7/8 ROIs, and 343 off by about 4 ULP.
   - All 294 recorded refinement boxes reproduce from `hem0`.
4. **`ID`.** The RGB is bit-identical.
5. **Direction.** This catches k applied to the wrong column. On tissue pixels with every channel in [20, 250] before
   and after, and with `hem0` and `hem'` > 0.01:
   - Regress `Δhem` through the origin on `(k−1)·max(c,0)·(s @ hed_from_rgb[:, 0]) / 6`. Here s is the scaled stain
     vector, and the factor 1/6 converts log10 OD to `rgb2hed`'s `ln(·)/ln(1e-6)` units.
   - For H06, H14 and E15, the slope must lie in [0.9, 1.1].
   - Review measured 1.008, 1.011 and 1.005 on 201, 0.168 without the 1/6, and 0.04 or 12.9 with the columns swapped.
   - **Label assert, per ROI:** `M[:,0] @ hed_from_rgb[:,0] > M[:,1] @ hed_from_rgb[:,0]`, so H responds more in
     `rgb2hed`'s H channel than E does. The slope gate alone cannot catch a basis returned as [e, h].
6. **`NULL2`.** Every estimator's `b'` equals exactly 2× its unperturbed value, and every `c_eq` is unchanged.
7. **`c_eq` self-consistency.** On the original image, `|c_eq(c·b0) − c| < 0.01` for c ∈ {0.25, 0.5, 0.75, 0.9, 1.0,
   1.5, 2.0}.
8. **Calibration.** `GAMMA` and `JPG95` fire nothing; report their per-ROI jitter.
9. **Determinism.** N = 1 equals N = 3 on the smoke ROIs, after the canonical sort.
10. **Macenko sanity.** Report the angles against `macenko.csv` (review: < 0.1° on 529) and `‖R‖/‖OD‖` (review: 0.049).
    This is reported, not halted on.

### 4.5 Runtime

- **Measured in review:** 22.4 s per ROI on 529 (the largest) for 16 variants plus setup, single-threaded with blocked
  `hem`. Add about 3–5 s for boxes, TOTSU/TW128, SHIFT, NULL2 and ID.
- **Full run:** about 22 CPU-min, so **~2 min wall at N ≈ 14**. The earlier ~5 min estimate was conservative.
- **Smoke run** (529, 128, 401, at N = 3 and N = 1): about 3 min.
- **Building the module and gates** is the main cost.
- Keep the lid open and stay on AC power during the run.

**Stage B caveat carried forward.** H06 cannot lighten copied zero-channel pixels: 0.14 % of 211's pixels, 0 on the
other pilot ROIs. They stay as dark specks in lightened nuclei. The per-ROI count travels with any Stage B H06 search.

## 5. Adversarial concerns (registered in the notebook before the run)

| concern | test |
|---|---|
| A1 hides rank changes | stated as a limit; rank-altering perturbations go to Stage B if it runs |
| quantisation or noise fires the rule | PL averages over 294 landings (one ROI is 1/49); the calibration gate (`GAMMA`, `JPG95`, `NULL2`); bin width reported beside top PL contributors |
| the synthesis is wrong but the gates pass | the direction gate (5); `ID` covers only encoding |
| the rule was revised after seeing pilot data | PL is anchored to the measured ladder and the NI margin, not to pilot values; the pilot is reported in §3.6 |
| cross-machine float drift | gate 3 within one bin; exact count reported |
| thread oversubscription or memory exhaustion | `threadpool_info` assert; RSS measured on the largest ROI |
| worker scheduling changes results | per-ROI seeds; canonical sort; determinism gate |

## 6. Review log

### Round 1 (reviewer agent, 2026-10-08; scripts and outputs in `floor_robustness/stageA_review/`)

| # | finding | verdict | action |
|---|---|---|---|
| 1 | The `Rres` stop clause is unreachable for JPG/BLUR: the perturbation's own smear exceeds the 0.9× reference even with the floor placed perfectly | **accepted** | `Rres` dropped entirely; it cannot separate estimator error from perturbation smear |
| 2 | The grey zone is the likely outcome, and "user decides" is not a rule | **accepted** | no grey zone; one band [0.9, 1.5] anchored to search evidence |
| 3 | Estimators that cannot ship trigger Stage B; TW128 was never searched | **partly accepted** | TOTSU/TW128 are descriptive. **`b_otsu` stays in the decision set**: it is the adopted estimator (D3), the product computes it over its frame, and Stage B tests it as the 1.0× arm |
| 4 | No overshoot side; the ranking is undefined | **accepted** | the band has an upper end of 1.5, R2 covers > 2.0, and the ranking claim is dropped |
| 5 | GAIN and WB are not null under the `rgb2hed` clamp | **accepted** | added GAIN09/11 and WB R/B ±5 %; GAMMA added as the true RGB near-null |
| 6 | `M.T @ R` crashes, and the orthogonality gate cannot fail | **accepted** | replaced by the direction gate (5); I added the 1/6 unit factor the fix omitted |
| 7 | The dry run cannot fail | **accepted** | replaced by the calibration gate on `GAMMA`/`JPG95`/`NULL2` |
| 8 | T2's 0.25-drop clause never binds | **accepted** | removed; R2 is a plain band edge |
| 9 | The smoke run underestimates memory (5.08 GB on 529) | **accepted** | 529 is in the smoke set; N comes from its RSS; blocked `hem` is optional and gated |
| 10 | Thread limits need a hard check | **accepted** | `threadpoolctl` assert; no gate on `cv2.getNumThreads` |
| 11 | The determinism gate fails spuriously, and the smoke `roi_index` would differ | **accepted** | canonical sort; `roi_index` taken from the full list |
| 12 | The exact-match count misleads (parser, ULP drift) | **accepted** | `round_trip`, a float32 compare, a one-bin tolerance, and the count reported |
| 13 | The predictions were wrong (E15, SHIFT, TW128) | **accepted** | §3.6 now reports pilot observations, labelled as seen, not predicted |
| 14 | Synthesis details: zero channels, negative cH, SHIFT underspecified, NaN handling, A1 wording | **accepted** | zero-channel pixels unchanged, positive-part scaling, SHIFT fully specified, NaN rule, and A1 reworded |
| 15 | Scope: much of the plan cannot change the decision | **accepted** | cut `c_eq_local`, overlap/J_rank, the naive ratio, drift, ABS, TQ, bootstrap CIs, and scanner splits; kept box acceptance and the synthesis checks |

### Round 2 (same reviewer, plan v2)

No [BLOCKS] findings.

| # | finding | verdict | action |
|---|---|---|---|
| 1 | The 0.9 edge is a tail anchor (W128's per-click 5th percentile) applied to a median; a median at 0.9 puts about 46 % of clicks below anything searched. The upper edge has the same mismatch. | **accepted; I chose neither proposed fix** | The band was replaced by **PL**: the ladder as a dose–response curve, averaged over all 294 landings, firing at ≥ 0.25 pt, a quarter of the NI margin. Unlike (a), it is anchored to measured numbers. Unlike (b), it is not a bare quantile, and a single pigment-ROI bin snap cannot fire it. The over side keeps a plain extrapolation guard. |
| 2 | `b_otsu` transfers only if the product's search region is about ROI-sized | **accepted as a condition; resolved** | the user decided on 2026-10-08 to keep the default ROI size; W512 not needed |
| 3 | Gate 5 cannot catch a basis returned as [e, h] | **accepted** | label assert added; regression through the origin specified |
| 4 | `b_otsu` quantisation on pigment ROIs (one bin ≈ 0.2–0.3 in `c_eq`) | **accepted, reduced** | under PL, one ROI contributes 1/49 of the mean; the bin width is reported and flagged beside the top PL contributors |
| 5 | Zero-channel copy is asymmetric under H06 | **accepted** | carried as a Stage B caveat (§4.5) |
| 6 | The synthesis share is about 30 %, not 15 % | **accepted** | text fixed; the GPU decline holds |
| 7 | Stage B is near-certain | **accepted** | stated in §3.6 |
| Q6 | Blocked `hem` is bit-identical, and RSS falls to 2.46 GB | **adopted unconditionally** | the assert is kept; N ≈ 14 |

### Round 3 (same reviewer, plan v3): converged

No [CHANGES-DECISION] findings. PL is sound, and its approximations lean toward firing:
- **The pooled curve.** The per-ROI response does not correlate with landing covariates; the one near-signal, absolute
  `b_otsu`, points the conservative way. Per-ROI curves would be unbiased but have an SE of 0.09–0.26 pt.
- **The linear chord.** The measured curve is concave above 0.5, so the chord overstates the loss.
- **No overshoot credit.**
- **The check against a search.** PL predicts −0.056 pt for the searched W128 arm; the observed value was +0.24 pt
  [−0.10, +0.58].

| # | finding | verdict | action |
|---|---|---|---|
| 1 | [BLOCKS, text] "NULL2 gives PL = 0" is false for W128 (0.056) | **accepted** | the gate now reads "equals the unperturbed PL exactly" |
| 2 | `loss(c)` is undefined for c > 2 | **accepted** | loss = 0 for every c > 1; OVER handles the ROI level |
| 3 | stale text from the in-place edits (§3.6 band, §4.5 runtime and list order, §5 "raised the firing rate", "its cause", Stage B N) | **accepted** | all fixed |
| 4 | subtracting W128's unperturbed PL is optional | **declined** | the absolute rule is kept. It is more conservative, and on the pilot it does not change the Stage B set. |

**Pilot PL on 7 extreme ROIs** (not a 49-ROI estimate):
- **Fires:**
  - H06 × `b_otsu`: 0.456, still 0.334 with the two largest contributors dropped;
  - H06 × W128: 0.262, borderline;
  - E15 × `b_otsu`: 0.325, borderline.
- **Below the threshold:** CDN, GAIN09, WBR095 and JPG on `b_otsu` (0.12–0.22, driven by 128 and 211), and everything
  else at ≤ 0.10.
- **Calibration:** GAMMA and JPG95 give ≤ 0.083 and fire nothing.
- **OVER** fires nowhere.
