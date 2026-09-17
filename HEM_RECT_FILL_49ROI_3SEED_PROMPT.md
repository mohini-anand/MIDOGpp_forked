# Task: default 51 px template vs. the same template with its surround flattened outside the hematoxylin rectangle

## Goal
Compare two templates on the same 147 clicks (49 ROIs x 3), ranked by `chromatin_od`:
- `default_51`: the 51x51 `hematoxylin_od` window centred on the click, cut the way production cuts it.
- `hem_rect_fill`: the same 51x51 array. Pixels inside the hematoxylin Otsu rectangle keep their values. Every pixel outside the rectangle is set to the mean of those outside pixels.

Only the template's pixel values differ. These are identical in both arms:
- window position and template size;
- search image;
- peak, NMS and self-hit settings, with the self-hit centre at the click;
- the `chromatin_od` window;
- scoring.

Question: does precision@K change?

Repo: `/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked`, branch `find-and-suppress-midog`.

## Read first
- `midog_utils/production.py`, `find_and_suppress.py`, `template_match.py`, `seed_selection.py`, `chromatin.py`, `evaluate.py`.
- Code cells of `production_hematoxylin_only/bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb` ("the 49-ROI notebook"). Functions named below are copied from it unchanged.

If a step can't run as written, or is wrong, stop and report. Don't work around a halt.

## Inputs (read-only)
- `B3 = results/precision_at_k_49roi_3seed_chromatin_bbox3way`. Files: `{B3}_per_run.csv`, `{B3}_top30.csv`, `{B3}_delta_per_roi.csv`, `{B3}_delta_stats.csv`. Read each with `pd.read_csv(..., float_precision='round_trip')`.
- `databases/MIDOG++.json`, `images/extra_valid/*.tiff`, `images/extra_valid/testing_set/*.tiff`.

## Before running
- `git status --short midog_utils/` prints nothing. Don't modify `midog_utils`.
- `{B3}_per_run.csv` passes all of these:
  - 441 rows;
  - `condition` values are exactly `default_51`, `gray_bbox`, `hem_bbox`, with 147 rows each;
  - 49 `file_name`s, each with `seed_index` 0, 1 and 2 under every condition;
  - 3 distinct `seed_ann_id` per `file_name`.
- Already verified: a dry run of Cells 3-4 on `300.tiff` and `509.tiff` reproduced the recorded `default_51` results on 6/6 clicks. That covered TP@K, `n_peaks`, `n_detections`, `n_self_hits`, the top-30 `rank`/`cx`/`cy`/`bucket`/`matched_ann_id`, float32 `score`, and bit-exact `od`.

## Notebook: `production_hematoxylin_only/hem_rect_fill_vs_default_chromatin_od_49roi_3seed.ipynb`
nbconvert runs it from its own directory, so every path starts with `../`.

### Cell 0 (markdown), written before the run
State the goal, the two arms, and what is identical. Also state these reading rules, fixed now:
- **Bootstrap CI excludes 0, delta < 0:** flattening the surround's pattern lowers precision.
- **Bootstrap CI excludes 0, delta > 0:** flattening it raises precision.
- **CI includes 0:** no difference detected. Quote the CI as the bound. Don't claim equivalence.

Caveats:
- `chromatin_od` ranking only.
- The clicks are inherited from the 49-ROI joint gate.
- One fill rule: the mean of the outside pixels.
- The rectangle comes from `tighten_box_otsu` at its default gate settings.
- Each arm's peak threshold comes from its own response map, so the two arms' candidate lists can differ in length (by up to 3 in the dry run). Precision is compared at fixed K.
- The filled template keeps `default_51`'s overall mean, up to float rounding.

### Cell 1: config
- **Imports:**
  - `gc, hashlib, itertools, math, os, sys, time`;
  - `from fractions import Fraction`;
  - `cv2`, `numpy as np`, `pandas as pd`, `from scipy import stats as sps`, `matplotlib.pyplot as plt`, `from matplotlib.patches import Rectangle`;
  - `sys.path.insert(0, '..')`;
  - from `midog_utils`: `channels as ch`, `chromatin as cm`, `dataset as ds`, `evaluate as ev`, `invariants as inv`, `production as prod`, `seed_selection as ss`, `template_match as tm`;
  - `from midog_utils.nms import nms_by_distance`.
- **Config checks:**
  - Copy the 49-ROI notebook's `REFERENCE_CONFIG` / `CONFIG_DRIFT` block and its assert.
  - Assert `ev.MIDOG_RADIUS_UM == 7.5`, `tm.BASE_SIZE == 51`, `tm.PATCH_SIZE == 73` and `prod.OD_WINDOW // 2 == 25`.
- **Constants:**
  - `SMOKE_N_ROIS = None`, and `N_ROIS = 49 if SMOKE_N_ROIS is None else SMOKE_N_ROIS`.
  - `ARMS = ('default_51', 'hem_rect_fill')`. Every delta is `hem_rect_fill - default_51`.
  - `BUDGETS = (10, 20, 30)`, `SEED_INDICES = (0, 1, 2)`, `N_BOOT = 10_000`, `BOOT_SEED = 20260917`.
  - `B3 = '../results/precision_at_k_49roi_3seed_chromatin_bbox3way'`.
  - `SUBSET_DIRS = {'original_14': '../images/extra_valid', 'testing_35': '../images/extra_valid/testing_set'}`.
- **Outputs:**
  - `OUT_STEM = '../results/precision_at_k_49roi_3seed_chromatin_hemrectfill' + ('' if SMOKE_N_ROIS is None else '_smoke')`.
  - `OUT = {name: f'{OUT_STEM}_{name}.csv' for name in ('clicks', 'per_run', 'top30', 'per_roi', 'delta_per_roi', 'delta_stats', 'verification')}`.
  - `FIG_PATH = 'hemrectfill_49roi3seed_templates.png'`.
  - Assert no `OUT` value starts with `B3 + '_'` or `'../results/precision_at_k_14roi_'`.

### Cell 2: clicks and recorded reference
- **Recorded tables:**
  - `REC_RUN = pd.read_csv(f'{B3}_per_run.csv', float_precision='round_trip')`.
  - `REC_TOP = pd.read_csv(f'{B3}_top30.csv', float_precision='round_trip')`, keeping rows with `condition == 'default_51'`.
- **Clicks:**
  - `CLICKS` = the `REC_RUN` rows with `condition == 'default_51'`, in CSV row order.
  - `ROI_ORDER` = the first `N_ROIS` distinct `file_name`s in that order. Keep only `CLICKS` rows for those ROIs.
  - Merge in the `hem_bbox` rows' `base_size`, `tpl_cx`, `tpl_cy` on (`file_name`, `seed_index`), as `rec_hem_base`, `rec_hem_cx`, `rec_hem_cy`.
- **Annotations:** `images, annotations = ds.load_annotations('../databases/MIDOG++.json')`.
- **Asserts:**
  - Per click, exactly one `annotations` row has `ann_id == seed_ann_id` and `file_name == file_name`, and it has `category_id == ds.MITOTIC`, `cx == click_cx` and `cy == click_cy`.
  - `len(CLICKS) == 3 * N_ROIS`.

### Cell 3: functions
Style: one-line `def`s and calls, and the docstring format under Writing rules.
1. `cut_templates(hem, cx, cy)`:
   - `patch = tm.read_padded_patch(hem, cx, cy, 73)`. Assert not `None`.
   - `templates, _ = tm.build_augmentations(patch, 51, (1.0,), 1, (False,))`. `T = templates[0]`. Assert `T.shape == (51, 51)` and `T.dtype == np.float32`.
   - `W = tm.read_padded_patch(hem, cx, cy, 51)`. Assert `np.array_equal(W, T)`: the Otsu window and the template are the same pixels.
   - `bbox = ss.tighten_box_otsu(W)`, with default arguments. Assert not `None`. `y0, y1, x0, x1 = bbox`, half-open and template-local.
   - `outside = np.ones((51, 51), bool)`, then `outside[y0:y1, x0:x1] = False`. Assert `outside.any()`.
   - `fill = np.float32(T[outside].astype(np.float64).mean())`. `F = T.copy()`, then `F[outside] = fill`. Assert `F.dtype == np.float32` and `F.flags['C_CONTIGUOUS']`.
   - Returns `(T, F, (y0, y1, x0, x1), fill)`.
2. `run_template(hem_padded, shape, template, cx, cy, mpp)`. It performs the same operations, in the same order, as `find_and_suppress` inside `run_production_pipeline(rank_key='chromatin_od')`, for one 51 px template:
   1. `nms_radius = ev.radius_px(mpp)`. Call `inv.check_nms_radius(nms_radius, mpp, label='hem_rect_fill')`.
   2. `pad = 25`. `h, w = shape`.
   3. `fused_p, best_p, valid_p = tm.fused_response(hem_padded, [template], False, method=prod.TM_METHOD)`. `fused = fused_p[pad:pad + h, pad:pad + w]`, `valid = valid_p[pad:pad + h, pad:pad + w]`. Assert `valid.all()`.
   4. `med, mad = tm.robust_stats(fused, valid)`. `threshold = med + prod.DEEP_FLOOR_Z * mad`.
   5. `centers, scores = tm.extract_peaks(fused, valid, prod.PEAK_MIN_DISTANCE, threshold, prod.MAX_PEAKS)`. `n_peaks = len(centers)`.
   6. `keep = nms_by_distance(centers, scores, nms_radius)`. `centers, scores = centers[keep], scores[keep]`.
   7. `self_hit = np.hypot(centers[:, 0] - float(cx), centers[:, 1] - float(cy)) <= prod.SELF_HIT_RADIUS`. `n_self_hits = int(self_hit.sum())`. `centers, scores = centers[~self_hit], scores[~self_hit]`.
   8. `det = pd.DataFrame({'rank': np.arange(len(centers)), 'cx': centers[:, 0], 'cy': centers[:, 1], 'score': scores})`.
   9. `shifted = det.assign(cx=det['cx'] + 25, cy=det['cy'] + 25)`. `det = det.assign(od=cm.score_detections(shifted, hem_padded, window=prod.OD_WINDOW)['od'].to_numpy())`.
   10. `det = det.sort_values('od', ascending=False, na_position='last', kind='mergesort').reset_index(drop=True)`. `det = det.assign(rank=np.arange(len(det)))`.
   11. Returns `(det, {'n_peaks': n_peaks, 'n_detections': len(det), 'n_self_hits': n_self_hits})`.
3. `score_top_k(detections, gt_eval, match_radius)`: copied from the 49-ROI notebook.

### Cell 4: run loop
For each `file_name` in `ROI_ORDER`, with `rows` = its 3 `CLICKS` rows sorted by `seed_index`:
1. **Load:**
   - `path = f"{SUBSET_DIRS[rows['subset'].iat[0]]}/{file_name}"`. Assert `os.path.exists(path)`.
   - `rgb = ds.load_roi(path)`, `mpp = ds.roi_mpp(path)`. Assert `mpp == rows['mpp'].iat[0]`.
2. **Prepare:**
   - `hem = ch.to_channel(rgb, 'hematoxylin_od')`.
   - `hem_padded = cv2.copyMakeBorder(hem, 25, 25, 25, 25, borderType=cv2.BORDER_REPLICATE)`.
   - `gt = ds.image_annotations(annotations, file_name)`.
   - `match_radius = ev.radius_px(mpp)`. Assert `match_radius == rows['match_radius_px'].iat[0]`.
3. `for _, r in rows.iterrows():` with `cx = float(r.click_cx)` and `cy = float(r.click_cy)`:
   1. **Templates:**
      - `T, F, (y0, y1, x0, x1), fill = cut_templates(hem, cx, cy)`.
      - Assert `ss.tightened_template_box(hem, cx, cy) == (r.rec_hem_base, r.rec_hem_cx, r.rec_hem_cy)`, so the rectangle is the one `hem_bbox` used.
      - `sha_T = hashlib.sha1(T.tobytes()).hexdigest()`, `sha_F = hashlib.sha1(F.tobytes()).hexdigest()`. Assert `sha_T != sha_F`.
   2. **Ground truth:** `gt_eval = gt[gt['ann_id'] != r.seed_ann_id].reset_index(drop=True)`. Assert `int((gt_eval['category_id'] == ds.MITOTIC).sum()) == r.n_gt_mitotic`.
   3. **`default_51`:** time this step as `t_d`. `det_d, info_d = run_template(hem_padded, hem.shape, T, cx, cy, mpp)`, then `out_d, tp_d = score_top_k(det_d, gt_eval, match_radius)`.
   4. **Parity gate (assert; runs before the filled arm).** With `ref = REC_TOP[(REC_TOP['file_name'] == file_name) & (REC_TOP['seed_index'] == r.seed_index)].sort_values('rank')` and `top = out_d.head(30)`:
      - `(info_d['n_peaks'], info_d['n_detections'], info_d['n_self_hits'])` equals `(r.n_peaks, r.n_detections, r.n_self_hits)`;
      - `tp_d[k] == r[f'tp_at_{k}']` for every K;
      - `len(ref) == 30`;
      - `np.array_equal(top[c].to_numpy(), ref[c].to_numpy())` for `c` in `rank, cx, cy, bucket, matched_ann_id`;
      - `np.array_equal(top['score'].to_numpy(), ref['score'].to_numpy().astype(np.float32))`. The CSV stores float32's shortest decimal form, so compare at float32;
      - `np.isclose(top['od'].to_numpy(), ref['od'].to_numpy(), rtol=1e-12, atol=0).all()`. Record `n_od_bit_exact = int((top['od'].to_numpy() == ref['od'].to_numpy()).sum())`.
   5. **`hem_rect_fill`:** time this step as `t_f`. `det_f, info_f = run_template(hem_padded, hem.shape, F, cx, cy, mpp)`, then `out_f, tp_f = score_top_k(det_f, gt_eval, match_radius)`.
   6. **`per_run` rows**, one per arm, from that arm's `out`, `info`, `tp`, sha and time:
      - keys and click fields: `subset, domain, file_name, image_id, seed_index, condition, seed_ann_id, click_cx, click_cy, mpp, match_radius_px, n_gt_mitotic`;
      - `tpl_sha1`, `n_peaks`, `n_detections`, `n_self_hits`;
      - `n_near_click` = count of rows with `np.hypot(cx - click_cx, cy - click_cy) <= match_radius`;
      - `near_click_detail` = `';'.join(f"rank{rank}:{bucket}:ann{matched_ann_id}:d_click={d:.1f}")` over those rows;
      - `od_nan = int(out['od'].isna().sum())`;
      - `od_descending = bool(np.all(np.diff(out['od'].to_numpy()) <= 0))`;
      - `t_run_s`, `tp_at_10`, `tp_at_20`, `tp_at_30`;
      - `n_od_bit_exact`, only on the `default_51` row (NaN on the filled row).
   7. **`clicks` row:**
      - `file_name, seed_index, seed_ann_id`;
      - `rect_y0, rect_y1, rect_x0, rect_x1`, `rect_h = y1 - y0`, `rect_w = x1 - x0`;
      - `rect_touches_window_edge = y0 == 0 or x0 == 0 or y1 == 51 or x1 == 51`;
      - `n_filled_px = 2601 - rect_h * rect_w`, `filled_share = n_filled_px / 2601`, `fill_value`;
      - `tpl_mean_default = T.astype(np.float64).mean()`, `tpl_mean_filled = F.astype(np.float64).mean()`.
   8. **`top30` rows:** `out.head(30)` for each arm, columns `rank, cx, cy, score, od, bucket, matched_ann_id, matched_category`, plus `file_name, seed_index, condition`.
   9. **For the figure:** `TEMPLATES[(file_name, r.seed_index)] = (T, F, (y0, y1, x0, x1))`.
4. **Persist and free:**
   - Write `OUT['clicks']`, `OUT['per_run']` and `OUT['top30']` from all rows so far.
   - `del rgb, hem, hem_padded`, then `gc.collect()`.
   - Print one progress line per ROI: file, elapsed minutes, and `tp_d`/`tp_f` per click.

### Cell 5: verification
- **Halt** (assert after writing the CSV):
  - every run has `n_detections >= 30`, `od_nan == 0` and `od_descending`;
  - every (`file_name`, `seed_index`) has exactly one row per arm, with equal `seed_ann_id` and `n_gt_mitotic`;
  - every click has `np.isclose(tpl_mean_default, tpl_mean_filled, rtol=0, atol=1e-6)`.
- **Report** (recorded, no assert):
  - `n_peaks == prod.MAX_PEAKS`;
  - `n_near_click == 0`, with `near_click_detail`;
  - the `n_self_hits` value counts per arm;
  - the sum of `n_od_bit_exact` out of `30 * len(CLICKS)`.
- Write `OUT['verification']` with columns `check, severity, file_name, seed_index, condition, passed, detail`. Use `'ALL'` where a check has no such key.

### Cell 6: fill geometry (descriptive)
From the `clicks` table, print:
- median, min and max of `rect_h`, `rect_w`, `filled_share` and `fill_value`;
- the count with `rect_touches_window_edge`;
- the count with `max(rect_h, rect_w) >= 50`.

### Cell 7: precision
- Add `precision_at_{K} = tp_at_{K} / K` to `per_run`. Rewrite `OUT['per_run']`.
- `per_roi`: per (`file_name`, `condition`), `tp_at_{K}` summed over the 3 clicks, and `precision_at_{K} = tp_at_{K} / (3 * K)`. Write `OUT['per_roi']`.
- Print pooled precision per arm and K: `sum(tp_at_K) / (3 * N_ROIS * K)`.

### Cell 8: statistics
ROI is the unit. One pair, K = 10, 20, 30.
- **Copy unchanged from the 49-ROI notebook:** `sign_flip_null_counts`, `exact_sign_flip_p`, `brute_force_sign_flip_p`, `holm_adjust`, `boot_ci_exact`, `d5_robustness`, `BOOT_SHARE_BAND = 2 * math.sqrt(0.025 * 0.975 / N_BOOT)`.
- **Validation (assert):**
  - With `check_rng = np.random.default_rng(0)`, for n in 1..12 and 20 vectors each, `d = check_rng.integers(-6, 7, size=n)` gives `exact_sign_flip_p(d)[0] == brute_force_sign_flip_p(d)`.
  - For each K: take `d` = `delta_tp_sum` from `{B3}_delta_per_roi.csv` where `pair == 'hem_bbox - default_51'` and `K == K` (49 rows). `exact_sign_flip_p(d)[0]` equals, under `np.isclose(rtol=1e-12, atol=0)`, the `exact_p` of `{B3}_delta_stats.csv` where `scope == 'subset'`, `group == 'all_49'`, `pair == 'hem_bbox - default_51'` and `K == K`.
- **`delta_per_roi`:** one row per (`file_name`, K), in `ROI_ORDER`. Columns: `delta_tp_s0`, `delta_tp_s1`, `delta_tp_s2` (`hem_rect_fill - default_51` TP@K per click), and `delta_tp_sum`. Write `OUT['delta_per_roi']`.
- **Bootstrap indices:** `BOOT_IDX = np.random.default_rng([BOOT_SEED, 0]).integers(N_ROIS, size=(N_BOOT, N_ROIS))`, shared by every CI below.
- **Per K**, with `d` = the `N_ROIS` values of `delta_tp_sum`:
  - `pooled_delta_precision = d.sum() / (3 * N_ROIS * K)`;
  - `wins = (d > 0).sum()`, `losses = (d < 0).sum()`, `ties = (d == 0).sum()`;
  - `exact_p, n_nonzero = exact_sign_flip_p(d)`, and `min_attainable_p = 2.0 ** (1 - n_nonzero)` (1.0 if `n_nonzero == 0`);
  - `holm_p = holm_adjust` over the 3 `exact_p`;
  - `ci_low, ci_high, boot_share_ge_0, boot_share_le_0 = boot_ci_exact(d, BOOT_IDX, 3 * N_ROIS * K)`;
  - `boot_share_0_or_opposite = boot_share_ge_0 if pooled_delta_precision < 0 else boot_share_le_0`;
  - with `x = d / (3 * K)`, `half = sps.t.ppf(0.975, N_ROIS - 1) * x.std(ddof=1) / math.sqrt(N_ROIS)`, giving `t_ci_low = x.mean() - half` and `t_ci_high = x.mean() + half`;
  - `ci_excludes_0 = ci_low > 0 or ci_high < 0`, and `t_ci_excludes_0` the same for the t interval;
  - `n_same_direction` = `wins` if the delta > 0, `losses` if < 0, else 0;
  - `majority_of_all_rois = n_same_direction > N_ROIS / 2`;
  - `d5_bar_met = ci_excludes_0 and majority_of_all_rois`;
  - column `d5_robustness` = the return value of the function `d5_robustness(row)`. `row` is a `pd.Series` with keys `ci_excludes_0`, `t_ci_excludes_0`, `majority_of_all_rois` and `boot_share_0_or_opposite`. Don't assign the result to a variable named `d5_robustness`, or the function is overwritten.
  - Write `OUT['delta_stats']`.
- **Pooled precision CIs per arm and K:** `boot_ci_exact(per_roi tp_at_K in ROI_ORDER, BOOT_IDX, 3 * N_ROIS * K)`.
- **Descriptive (print only, no p-values):**
  - per domain and K: pooled delta in points, and ROI-level W/L/T;
  - per K: click-level W/L/T over `len(CLICKS)` clicks.

### Cell 9: figure
- **Which clicks:** sort the `clicks` rows by (`filled_share`, position in `ROI_ORDER`, `seed_index`). With `n = len(CLICKS)`, take positions `0, 1, n//2 - 1, n//2, n - 2, n - 1`.
- **Layout:** `plt.subplots(6, 2, figsize=(5, 15))`, one row per click.
  - Left panel: `T`. Right panel: `F`.
  - Both panels use `imshow(..., cmap='gray_r', vmin=T.min(), vmax=T.max())`.
  - Both panels get the rectangle: `Rectangle((x0 - 0.5, y0 - 0.5), x1 - x0, y1 - y0, fill=False, edgecolor='red')`.
- **Title:** left panel `f"{file_name} s{seed_index} rect {rect_h}x{rect_w} filled {filled_share:.0%}"`.
- **Save:** `fig.savefig(FIG_PATH, dpi=150, bbox_inches='tight')`. The full run overwrites the smoke run's file.

### Cell 10 (markdown), written after the full run
Every number must be printed by a cell above. Report:
- the parity gate (runs checked, bit-exact `od` count);
- the report checks;
- the fill geometry;
- pooled precision per arm and K, with CIs;
- the per-K table: pooled delta in points, bootstrap CI, t CI, W/L/T, exact p, Holm p, `min_attainable_p`, D5 bar, `d5_robustness`.

Then read the result with Cell 0's rules only. Make no mechanism claims beyond what is measured.

## Running
- **Interpreter:** `/Users/mohinianand/anaconda3/bin/python3` and `/Users/mohinianand/anaconda3/bin/jupyter` only.
- **Command:** `/Users/mohinianand/anaconda3/bin/jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=10800 production_hematoxylin_only/hem_rect_fill_vs_default_chromatin_od_49roi_3seed.ipynb`.
- **Smoke run:** set `SMOKE_N_ROIS = 2` and run the command. Every cell must finish without error.
- **Full run:**
  - delete `results/precision_at_k_49roi_3seed_chromatin_hemrectfill_smoke_*.csv`;
  - set `SMOKE_N_ROIS = None`;
  - run the command in the background.
- **Time:** 294 template runs at 1.1-1.7 s each (measured), plus 49 ROI loads and conversions. Budget 15-30 min.
- **Finish:** write Cell 10. Final state: `SMOKE_N_ROIS = None`, execution counts contiguous from 1, no error outputs.

## Writing rules
- Function definitions and calls go on one line, never wrapped.
- Docstrings:
  - the opening `"""` alone on its line;
  - content indented one level deeper;
  - a one-line summary;
  - one `name (type): description` line per parameter;
  - one `Returns type: description` line.

## Don't
- Commit.
- Modify `midog_utils` or any `{B3}_*` file.
- Change a pipeline setting, add rankers, or use K > 30.
- Redraw clicks.
- Read any notebook's printed output as data.

## Report back
- Parity gate: runs checked, mismatches, bit-exact `od` count.
- Report checks.
- Fill geometry.
- Pooled precision per arm and K.
- The per-K statistics table.
- Runtime.
- Files written.
- Any deviation from this prompt, and why.
