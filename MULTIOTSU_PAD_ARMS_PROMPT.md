# Multi-Otsu template, four pad arms: does padding the tightened box help at K = 10?

Build and run `production_hematoxylin_only/multiotsu_pad_arms_49roi_6seed.ipynb`.

Reason independently. This prompt fixes the **design**, not the conclusion. If a check fails, if a value here
contradicts the repo, or if a result contradicts what you expected, stop and report it. Do not adapt the design
to make anything pass. A null is a result.

## What changes, and what does not

The adopted template-refinement method (`refine_template_box`, committed in `154863e`,
`production_hematoxylin_only/multiotsu_template_refinement_14roi.ipynb`) returns, for one click, a square box:
3-class multi-Otsu on the 51 px `hematoxylin_od` window around the click, top class kept, 1 px dilation, piece
selection (base case or 3.75 um fallback), and a square of side `bs = odd(longer bbox side) + 2 * pad_px`
(0.5 um margin a side; `pad_px` = 2 on every one of the 49 ROIs), centred on the **selection's bounding box**.

This run turns that box into four templates and measures precision@10 on the same 294 clicks as the window sweep
(`window_sweep_51_61_71_49roi_6seed.ipynb`). **Two things change** against the sweep:

1. the kept region is the new box (the whole `bs x bs` square, raw pixels), not the recorded Otsu rectangle;
2. every template is **centred on the box centre, not the click**.

**Everything else stays exactly as in the window sweep:** the clicks, the channel, the blank (51 x 51, centred on
the click), `TM_CCOEFF`, peak extraction, `MAX_PEAKS`, NMS, `od`, matching, both ranking axes, K = 10, and the
statistics code.

## The four arms

Let `bs` be the box side and `(tcx, tcy) = (int(round(center_xy[0])), int(round(center_xy[1])))` the box centre,
rounded **once**, with `read_padded_patch`'s own convention. Every template is an `S x S` window of the
`hematoxylin_od` channel cut at `(tcx, tcy)`. Its **kept region is always the central `bs x bs` square**,
`o = (S - bs) // 2`, `mask[o:o+bs, o:o+bs] = True`; every other pixel is set to the mean of those other pixels,
computed on that arm's own window.

| arm | `S` | template |
|---|---|---|
| `mo_box` | `bs` | the box itself, raw pixels, no fill |
| `mo_pad51` | 51 | the box, padded evenly to 51 x 51; the pad set to the mean of the pad's own pixels |
| `mo_pad61` | 61 | the same, padded to 61 x 61 |
| `mo_ring5` | `bs + 10` | the box plus a 5 px ring on every side, whatever `bs` is; the ring set to its own mean |

**"51 px pad" and "61 px pad" mean a total template side of 51 and 61**, the window sweep's `W`. A 51 px pad on
each side would make templates over 127 px, far past the 73 px limit below. If anything in the repo suggests the
other reading, stop and report it.

4 arms x 294 clicks = 1176 searches. Each search's one detection list is ranked on both axes.

Fill convention: copy the window sweep's `build_templates` (cell 6) for the dtype, the fill value
(`np.float32(T[~mask].astype(np.float64).mean())`) and the SHA-1. `bs` and `bs + 10` are odd, so every `S` is odd.

## Reference arms, from the saved CSVs; no extra searches

`default_51`, `fill_51` and `fill_61` come from the window sweep's saved outputs:
`../results/precision_at_k_49roi_6seed_windowsweep_postD11_{per_run,top30}.csv`, keyed by
`(file_name, seed_index, condition, rank_key)`. **Include `rank_key` in every join**, or the two axes duplicate.
They were searched with the identical pipeline on the identical clicks, and the pipeline is bit-deterministic.

**Parity gate, before the full run.** Copy the sweep's search code verbatim. Rebuild `default_51` and `fill_61`
exactly as the sweep did (click-centred, the recorded rectangle), re-run them on `337.tiff` and `509.tiff` (12
clicks), and require a row-for-row match to the saved sweep rows on both axes: `tp_at_10/20/30`, `n_peaks`,
`n_detections`, `tpl_sha1`, and top-30 `rank, cx, cy, bucket, matched_ann_id`; `score` after an `np.float32` cast,
`od` at `rtol = 1e-12`. `fill_61` exercises the `S != 51` padding path. Halt on the first mismatch.

## Inputs; no click selection is run

- `../results/precision_at_k_49roi_6seed_fillgeom_postD11_clicks.csv`: 294 rows. Use `file_name`, `seed_index`,
  `seed_ann_id`, `click_cx`, `click_cy`, `mpp`, `subset`, `domain`, and (for the parity gate only) the rectangle
  `rect_y0, rect_y1, rect_x0, rect_x1`. The click equals its annotation's centre on all 294 rows (checked).
- `refine_template_box`: copy it **verbatim** from cell 3 of `multiotsu_template_refinement_14roi.ipynb`, with its
  constants from cell 1. Do not edit it.
- **Tripwire.** 19 of the 294 clicks are among the method's 24 adoption clicks (join `file_name` + `seed_ann_id`
  to `ann_id` in `production_hematoxylin_only/multiotsurefine14_summary.csv`). On those 19, `accepted`,
  `base_size_px`, `center_x` and `center_y` must match exactly. Halt on a mismatch.
- Ground truth: `images, annotations = ds.load_annotations('../databases/MIDOG++.json')`, then
  `ds.image_annotations(annotations, file_name)`; drop the click's own annotation by `seed_ann_id`.
- ROI images: `{'original_14': '../images/extra_valid', 'testing_35': '../images/extra_valid/testing_set'}[subset] + '/' + file_name`.
- Bootstrap ROI order: the order `file_name` first appears in `..._fillgeom_postD11_per_run.csv`.

## Refusals

`refine_template_box` can refuse a click. Run it on all 294 clicks **before any search** and print the refusals
with their reasons.

- **More than 15 refusals (5%): stop before searching and report.**
- Otherwise, a refused click is dropped from **every** arm and **every** reference. The design stays paired.
  Reference contrasts are recomputed on the accepted clicks only. A ROI's value is its TP@10 summed over its
  accepted clicks, and the pooled denominator is 10 x the number of accepted clicks. A ROI with no accepted click
  leaves the analysis, and D5's majority is counted over the ROIs that remain.

Scope caveat, to state in the conclusion: the 294 clicks were drawn under the old Otsu joint gate.

## Pipeline, per (click, arm)

1. `hem = ch.to_channel(rgb, 'hematoxylin_od')`
2. the template, as above, from `tm.read_padded_patch(hem, tcx, tcy, S)`
3. `search, n_blanked = tm.blank_seed_square(hem, click_cx, click_cy, 51)`: **the click, 51, every arm**
4. pad `search` by `PAD = (S - 1) // 2`, `cv2.BORDER_REPLICATE`
5. `fused, _, valid = tm.fused_response(search_padded, [template], False, method=4)`; crop back; assert `valid.all()`
6. `med, mad = tm.robust_stats(fused, valid)`; `threshold = med + (-1.5) * mad`
7. `centers, scores = tm.extract_peaks(fused, valid, 7, threshold, 100)`
8. `keep = nms_by_distance(centers, scores, ev.radius_px(mpp))`
9. `od` via `cm.score_detections(...)` with `window = 51` on the **unblanked** channel. **Pad the `od` channel by
   25 px on every arm**, never by `PAD`. For `mo_box`, `PAD` is 12-22 px, and a 51 px `od` window near an ROI
   edge would run off a narrower pad. The sweep showed `od` is bit-identical for any pad of at least 25.
10. rank once per axis from this one list: `tm_score -> score`, `chromatin_od -> od`, descending,
    `kind='mergesort'`, `na_position='last'`
11. `ev.bucket_detections(...)` against `gt_eval`; TP@K = count of `ev.HUMAN_CORRECT_LABEL` in the top K

The detection coordinates are template centres, so the search and matching need no change for the new centre.

## Hyperparameters

| name | value |
|---|---|
| channel | `hematoxylin_od` |
| `TM_METHOD` | `4` (`cv2.TM_CCOEFF`, unnormalised) |
| `PEAK_MIN_DISTANCE` | `7` |
| `DEEP_FLOOR_Z` | `-1.5` (inert; `MAX_PEAKS` binds) |
| `MAX_PEAKS` | `100` |
| `OD_WINDOW` | `51`, unblanked channel, pad 25 |
| blanked square | **51 x 51, centred on the click, every arm** |
| NMS radius = match radius | `ev.radius_px(mpp)` (7.5 um, D7) |
| `TEST_K` | `10`; 20 and 30 recorded, quoted in no verdict |
| `N_BOOT` / `BOOT_SEED` | `10000` / `20260917` |
| max template side | 73 (`tm.PATCH_SIZE`) |

Assert `midog_utils.production` still matches this configuration, as the sweep's cell 1 does.

## Inline asserts (halt)

1. `n_blanked_px == 2601` on every row.
2. For every fill arm: `n_filled_px == S*S - bs*bs`, and the fill value equals
   `np.float32(T[~mask].astype(np.float64).mean())`, recomputed from the raw window inside the assert.
3. The template's central `bs x bs` square equals the raw `bs x bs` box read at `(tcx, tcy)`, pixel for pixel, on
   every arm. `mo_box` equals that box exactly.
4. `n_detections >= 10`.

## Flags: report, never halt, run anyway

Count each flag per arm, list the flagged clicks, and **pre-register a descriptive rerun of every tested contrast
with the flagged clicks excluded.** This matters most for `mo_ring5`, whose 5 px ring ignores the box size.

- **`box_outside_blank`**: the box reaches past the click-centred 51 x 51 blank
  (`abs(tcx - round(click_cx)) + (bs - 1) // 2 > 25`, or the same in y). Its raw kept pixels are then unblanked
  image, and the template can match itself.
- **`template_outside_window`**: the `S x S` footprint reaches past the click-centred 51 px window. The pad is
  constant, so this is informational, but report it (`mo_ring5` and `mo_pad61` will carry it).
- **`roi_edge_read`**: `tm.read_padded_patch` returns `None` because the box-centred `S` window crosses the ROI edge
  (`BORDER` = 36 guarantees only a 73 px window centred on the click). Build that template from a
  `cv2.BORDER_REPLICATE`-padded read instead.
- **`box_ge_S`** (`mo_pad51`, `mo_pad61`): `bs >= S`. There is nothing to pad; use the `bs x bs` box unfilled (so
  that click equals `mo_box`).
- **`mass_weight_nonpositive`**: `m_k - f <= 0`, where `m_k` is the template mean over the kept square and `f` the
  fill. A narrow ring (`mo_ring5`) can push `f` up toward `m_k`. At the template's own position the response is
  `k (m_k - f)(min - mu_S)`, so a negative mass weight makes it **positive** there. Print the minimum mass weight
  per arm.

## Contamination rule, fixed before the run

- **Self-match**: a non-TP detection within 2 px (Chebyshev) of the template centre `(tcx, tcy)` **or** of the
  click. Count it per (arm, axis) in the top 10 and in the top 30.
- **The rule**: if an arm has any self-match in its top 10 on an axis, every tested contrast containing that arm on
  that axis is labelled **contaminated**. Its pre-registered reading is then:
  - at `tm_score`: **resolved** only if both the self-match-dropped version and the worst-case bound give the same
    D5 verdict as the recorded version; otherwise **unresolved**;
  - at `chromatin_od`: **unresolved**. No bound exists there, because ranking is by `od`, and the peak that a
    freed pool slot admits can enter the top 10 anywhere.
- **Self-match-dropped version**: drop every self-match detection from each saved top-30 list and take the top 10
  again. Assert that each dropped detection is `non_human_findings` with `matched_ann_id == -1`, so no other bucket
  moves.
- **Worst-case bound, `tm_score` only**: count every self-match top-10 slot as a displaced TP, plus one TP for
  every other mitosis within `2 r + 3 px` of the click. The self-match peak lies within 2 sqrt(2) px of it; NMS
  radius = match radius. The post-run cell of the window sweep (`precision_at_k_49roi_6seed_windowsweep_postD11_selfmatch_sensitivity.csv`,
  the cell after "Post-run: how much the self-match moves the contrasts it enters") implements both. Reuse its logic.
- **Never lift or override a flag yourself.** Report it; the user decides.

**Near-seed count (a diagnostic, not a rule)**: non-TP top-10 detections within each click's actual reach, the
Chebyshev distance from the click within which the template overlaps the blank. That reach is
`25 + (S - 1) // 2 + max(abs(tcx - round(click_cx)), abs(tcy - round(click_cy)))`. Report also the sweep's fixed
60 px count, for comparison with the references (0 on `default_51` and `fill_51`).

**ROI edge**: `tp_at_10_excl_border` with the sweep's band, **35 px**, so it compares directly with the references'
saved column. Drop every detection with `min(cx, cy, w-1-cx, h-1-cy) < 35` identically in every arm, re-match the
rest against the same ground truth, and count TP in the top 10.

## Statistics

Copy the window sweep's statistics cell (cell 4) verbatim, **including its validation**: the fill-geometry run's
10 recorded `delta_stats` rows must reproduce before any search.

The unit is the ROI. Per contrast: the pooled delta in points, a 95% cluster-bootstrap CI over ROIs, a t interval,
ROI W/L/T, `n_nonzero`, an exact sign-flip p by convolution, `min_attainable_p`, Holm within family, and D5 (the CI
excludes 0 **and** a majority of ROIs move in the delta's direction). The two axes are **two outcomes**. They are
never pooled, and a result on one is not a replication of the other.

**Four Holm families, 10 tests, all at K = 10:**

| family | tests |
|---|---|
| `tm_score \| pad` | `mo_pad51 - mo_box`, `mo_pad61 - mo_box`, `mo_ring5 - mo_box` |
| `tm_score \| vs_incumbent` | `mo_pad51 - fill_51`, `mo_pad61 - fill_61` |
| `chromatin_od \| pad` | the same three |
| `chromatin_od \| vs_incumbent` | the same two |

`vs_incumbent` compares the new box with the recorded Otsu-rectangle fill at the same window side. It moves three
things at once: the kept region, the centre (box, not click) and the fill value. Say so, and do not attribute a
difference to any one of them.

**Descriptive, in no family, supporting no verdict:**
- every arm minus `default_51`;
- `mo_pad61 - mo_pad51`, `mo_ring5 - mo_pad51`, `mo_ring5 - mo_pad61`;
- `fill_61 - fill_51` on the accepted clicks;
- `tp_at_10_excl_border`;
- K = 20 and K = 30;
- the flagged-clicks-excluded rerun;
- the self-match-dropped version.

## Figures

1. Every K = 10 contrast with its CI, both axes, tested rows solid, descriptive rows hollow, with W/L/T and Holm.
   **Share the x-axis across the two panels**, or print on the figure that the scales differ.
2. A template gallery: the four templates side by side for 8 clicks, chosen to include the smallest and the
   largest `bs`, the largest centre offset, and one click from each flag that fires, with the click and the blank
   outline marked.

## Do not

- Do not analyse by domain. n = 7 per domain puts the Holm floor across 7 domains above 0.05.
- Do not redraw clicks or edit `refine_template_box`.
- Do not move, resize or re-centre the blank for any arm.
- Do not re-search the reference arms beyond the parity gate.
- Do not quote K = 20 or K = 30, or any descriptive row, in a verdict.
- Do not lift a flag.

## Outputs

`../results/precision_at_k_49roi_6seed_mopad_postD11_{templates,per_run,per_roi,delta_per_roi,delta_stats,top30,selfmatch_sensitivity}.csv`
(note `../`, since `nbconvert` sets the working directory to the notebook's own).

- `templates`: one row per (click, arm) with `accepted`, `refusal_reason`, `bs`, `S`, `tcx`, `tcy`,
  `offset_px`, `fill_value`, `mass_weight`, `n_kept_px`, `n_filled_px`, `tpl_sha1`, and every flag.
- `per_run`: the sweep's columns, plus the self-match and near-seed counts.

Figures: `production_hematoxylin_only/mopad_{contrasts,templates}.png` [files since removed; the figures are embedded in `multiotsu_pad_arms/multiotsu_pad_arms_49roi_6seed.ipynb`]. Rewrite `per_run` and `top30` after every
ROI, so a crash leaves partial results.

Cell 0 is the pre-registration: this design, the four arms, the families, the flags and the contamination rule.
Write it before the full run. Every number in the conclusion must be printed by a cell above it.

## Environment

- Interpreter: `/Users/mohinianand/anaconda3/bin/python3` (3.11.5). No other interpreter here can import `cv2`.
- Execute: `/Users/mohinianand/anaconda3/bin/jupyter nbconvert --to notebook --execute --inplace
  --ExecutePreprocessor.timeout=14400 <path>`. The per-cell limit must cover the whole search loop.
- Smoke first: support `MOPAD_SMOKE_ROIS=337.tiff,509.tiff` in the environment, which runs the identical code on
  those ROIs and writes `_smoke`-suffixed outputs. The recorded run has it unset.
- Runtime: the sweep's 2058 searches took about 92 min, so 1176 here should take about 50-60 min. The machine is
  CPU-only and shared with another long-running job, so expect overrun.
- Never use a pattern-matched `pkill` or `killall`. To stop your own run, kill only the exact PID you started.
- Do not commit. Do not modify any existing file. Create only the notebook, its outputs and its figures.

## Code conventions

Docstrings: the opening `"""` on its own line; then a brief summary, one line per parameter (type and a short
description) and one line for the return, all indented one level deeper than the quotes. Function definitions
and calls go on a **single line**, never wrapped. Match the window sweep notebook's style.
