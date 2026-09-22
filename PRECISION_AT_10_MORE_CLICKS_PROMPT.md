# Task: resolve precision@10 for the surround-fill arms, with more clicks per ROI

## Why this exists

Three earlier runs measured whether flattening the seed template's surround helps, all at **3 clicks per ROI**.
At K = 20 and K = 30 the answer resolved; at **K = 10 it never did**, and not because the effect is absent — because
the measurement is too coarse there. A ROI's TP@10 summed over 3 clicks is out of 30, so one detection moves that
ROI's precision by 3.3 points and ROIs tie constantly. In the most recent run, `hem_rect_fill - default_51` at
K = 10 was **+2.45 points, bootstrap CI [+0.48, +4.69]** at `tm_score` — the largest point estimate of any K — yet
**17 of 49 ROIs tied**, only 21 of 49 moved in the delta's direction, and Holm came out at 0.124. At
`chromatin_od`, 31 of 49 tied.

K = 10 is the part of the list a pathologist actually reads, so "not resolvable" is not an acceptable resting
place. **Your job is to re-run the comparison with more clicks per ROI so the K = 10 granularity is fine enough to
resolve it.** Nothing else about the method changes.

## What to build

One new notebook in `production_hematoxylin_only/`, executed with
`/Users/mohinianand/anaconda3/bin/jupyter nbconvert --to notebook --execute --inplace
--ExecutePreprocessor.timeout=14400 <path>`. Run it from the repo root; nbconvert sets the kernel's cwd to the
notebook's own directory, which is why the existing notebooks use `sys.path.insert(0, '..')` and `../images/...`.
Use `/Users/mohinianand/anaconda3/bin/python3` for everything — the other interpreters on PATH fail at `import cv2`.

### The four arms

All four cut a **51x51 `hematoxylin_od` window centred on the click**, exactly as production does
(`tm.read_padded_patch(hem, cx, cy, 73)` then `tm.build_augmentations(patch, 51, (1.0,), 1, (False,))[0][0]`).
Only the pixel values differ. In every filled arm the fill value is **the mean of the pixels that arm fills** —
recomputed per arm, since each keeps a different region.

| arm | kept region | template |
|---|---|---|
| `default_51` | everything | the window, untouched |
| `hem_rect_fill` | `ss.tighten_box_otsu(W)` | pixels outside it set to the mean of those outside pixels |
| `hem_rect_fill_1um` | that rectangle grown by `int(round(1.0 / mpp))` px on every edge, **each edge clipped to the 51 px window independently** | pixels outside the grown rectangle set to the mean of those outside pixels |
| `hem_mask_fill_1um` | **the accepted Otsu component itself, dilated by `int(round(1.0 / mpp))` px** — a non-rectangular mask | pixels outside the dilated mask set to the mean of those outside pixels |

**The last two arms differ only in the corners.** For any structuring element,
`bbox(component ⊕ B) == bbox(component) ⊕ bbox(B)`, so the bounding rectangle of the dilated component is exactly
`hem_rect_fill_1um`'s rectangle — verified identical on 12/12 sample clicks including two that touch the window
edge. The new arm is therefore **not** "pad the component then re-box it", which would be a byte-identical no-op.
It keeps the **padded mask itself**, which covers a median **72%** of that rectangle; the remaining 28% is corner
that was never nucleus. Dropping those corners is the whole point of the arm — if you replace the mask with its
bounding box anywhere in the implementation, you have silently rebuilt `hem_rect_fill_1um` and wasted the run.

To build the mask you need the component, which `ss.tighten_box_otsu` does not return — re-derive it with that
function's own internals (normalise the window to uint8 with `cv2.normalize(..., NORM_MINMAX)`, Otsu-threshold,
`skimage.measure.label(binary, connectivity=2)`, take the label at the centre pixel `[25, 25]`) and **assert that
its `regionprops` bbox equals `ss.tighten_box_otsu(W)`**, which proves you selected the same component production
does. Dilate with `cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1))`, `r = int(round(1.0 / mpp))`
(4 px on all 49 ROIs).

**The rectangle arms' fill values are recorded** as `fill_value` and `fill_value_1um`; gate against them. The mask
arm has no recorded reference of any kind — it is entirely new, so its correctness rests on the internal gates
below.

Do **not** change the fill rule. "Keep the default area filling method" means `mean(outside)`, the value the
recorded runs used. A previous investigation established that `mean(kept)` (which would make the ring weightless
under `TM_CCOEFF`) is 6 points *worse*, so the current fill is the right one and is not in question here.

### Clicks: extend to 6 per ROI, do not redraw

The recorded runs draw clicks per ROI with `rng = np.random.default_rng([seed_index, image_id])` over the
agreement-tier, border-filtered mitotic pool, removing each accepted annotation before the next draw, and
accepting only candidates that pass the joint gate (both `gray_inverted` and `hematoxylin_od` Otsu boxes succeed
and the resulting centre has a readable 73 px patch). `rerun_bbox3way_postD11.py` holds that procedure —
`pool_gate_specs`, `draw_joint_click`, `select_clicks`. **Reuse it.**

**Extend `seed_index` from `(0, 1, 2)` to `(0, 1, 2, 3, 4, 5)`.** Because each seed has its own RNG stream and the
pool shrinks by the accepted annotations, seeds 0-2 reproduce the recorded clicks exactly and seeds 3-5 are new.
Assert that: your seeds 0-2 must match the recorded `seed_ann_id`, `click_cx`, `click_cy` on all 49 ROIs. If they
do not, stop and fix that before anything else.

**Six is the largest balanced choice.** The joint-valid pool per ROI is min 6, median 18, max 162; all 49 ROIs
have at least 6, but only 42 have 9 and 37 have 12. Going past 6 would either drop ROIs (changing the population
and forfeiting comparability with every recorded run) or make the design unbalanced (ROIs would carry different
weights, and the ROI-level exact sign-flip assumes a common per-ROI scale). Stay at 6 and keep all 49 ROIs.

If 6 clicks still leaves K = 10 unresolved, **report that with the tie counts** — do not escalate to more clicks,
a different gate, or extra budgets on your own. The next option (a hematoxylin-only gate, which enlarges the pool
but forfeits the parity gates) is a decision for the user, not for you.

### Both ranking axes, from one search

Measure at **`rank_key = "chromatin_od"` and `rank_key = "tm_score"`**. These are two different outcomes, not two
looks at one, so give each its own Holm family.

Run the template search **once** per (click, arm) and let both axes sort that same detection list — the axes differ
only in the final sort key (`production.AXES` maps `tm_score -> score`, `chromatin_od -> od`). This makes the axis
comparison exact rather than approximate and costs one extra sort instead of a second search. Assert it: the two
rankings must be permutations of the same `(cx, cy, score, od)` rows on every run.

`shape_vs_mass_template_terms_49roi_3seed.ipynb` already does this (`search_once` / `rank_by`). Copy that harness.

### Budget

**Only K = 10 is tested.** Every hypothesis test, Holm family and D5 verdict is at K = 10 alone. Record TP@20 and
TP@30 in the per-run CSV as free descriptive columns, but do not put them in any family or quote them as results.

## Hold these identical across arms

- the 294 clicks (49 ROIs x 6), in a fixed order;
- window position (centred on the click) and template size (51 px), so `tm.blank_seed_square` blanks the **same**
  51x51 square in every arm and no base-size confound enters;
- search image (`hematoxylin_od`, D11-blanked then replicate-padded by 25 px);
- `PEAK_MIN_DISTANCE` 7, `DEEP_FLOOR_Z` -1.5, `MAX_PEAKS` 100, NMS radius = match radius = `ev.radius_px(mpp)`;
- the `chromatin_od` window (51 px), scored against the **unblanked** channel;
- greedy centre-distance matching at 7.5 um with the click's own annotation removed from the ground truth;
- a fixed K, since each arm's peak threshold comes from its own response map and list lengths can differ.

## D11 — read this before writing any pipeline code

D11 (commit `9e22eb9`) replaced the 5 px post-NMS self-hit filter with `tm.blank_seed_square`, which blanks the
seed's own template footprint from a **copy** of the search channel before correlation. **`production.SELF_HIT_RADIUS`
no longer exists.** Consequences:

- The recorded `hemrectfill` and `hemrectfill1um` runs are **pre-D11** and cannot be reproduced by current code.
- The recorded `bbox3way_postD11` and `fillalpha_postD11` runs **are** post-D11 and can be.
- D11 changed outcomes but **not** templates: the same click's `default_51` template SHA-1 is identical pre and
  post. So template-level gates still work against the pre-D11 runs even though detection-level ones do not.

Your run is post-D11. Assert `hasattr(tm, 'blank_seed_square')` and `not hasattr(prod, 'SELF_HIT_RADIUS')`.

## Gates — all of them, before any statistic is read

1. **Click reproduction.** Seeds 0-2 must match the recorded `seed_ann_id`/`click_cx`/`click_cy` on all 49 ROIs
   (`results/precision_at_k_49roi_3seed_chromatin_bbox3way_postD11_per_run.csv`).
2. **Harness fidelity, pre-D11, subset.** On 2 ROIs x 3 clicks, run your inlined pipeline in *pre-D11 mode* — no
   blanking, plus a 5.0 px post-NMS self-hit filter at the click — and require it to reproduce, bit for bit:
   `default_51` from `..._chromatin_bbox3way_{per_run,top30}.csv` and `hem_rect_fill` from
   `..._chromatin_hemrectfill_{per_run,top30}.csv`. This proves your reimplementation is faithful before the D11
   switch is flipped.
3. **Production-path equivalence, post-D11, subset.** Your post-D11 `default_51` path must equal
   `prod.run_production_pipeline(rgb, seed_stub, mpp, rank_key=axis)` row for row, at **both** axes
   (`seed_stub = types.SimpleNamespace(base_size=51, template_xy=(cx, cy))`).
4. **Cross-notebook parity, seeds 0-2, full run, `chromatin_od`.** `default_51` must reproduce
   `..._bbox3way_postD11_*`, and `hem_rect_fill` must reproduce the `fill_alpha_1` arm of
   `..._chromatin_fillalpha_postD11_*` — counts, TP@10/20/30, template SHA-1, and the top-30
   rank/cx/cy/bucket/matched_ann_id/score/od.
5. **Template parity for the 1 um arm, seeds 0-2.** No post-D11 detection reference exists, so gate the geometry
   instead: your grown rectangle and fill value must equal `rectm_{y0,y1,x0,x1}` and `fill_value_1um` in
   `results/precision_at_k_49roi_3seed_chromatin_hemrectfill1um_clicks.csv`, and the tight rectangle and fill must
   equal `rect_*`/`fill_value` in `..._chromatin_hemrectfill_clicks.csv`.
6. **Margin geometry, all clicks.** Each edge grows by `min(margin_px, room to that border)` and no more; the
   grown rectangle never leaves the 51 px window; the filled pixel count equals `2601 - h*w` of the grown
   rectangle; the three arms' templates are pairwise distinct on every click.
7. **One detection set per (click, arm).** The two axes' rankings must be permutations of the same rows.
8. **The mask arm, which has no external reference.** On every click assert: the re-derived component's bbox
   equals `ss.tighten_box_otsu(W)`; the dilated mask's bounding box equals `hem_rect_fill_1um`'s grown rectangle
   **exactly** (this is the identity above, and it proves the two 1 um arms differ only in corners); the mask is a
   strict subset of that rectangle with strictly fewer kept pixels; the fill equals the mean of the pixels outside
   the mask; and all four arms' template SHA-1s are distinct. Record `n_kept_px` and the mask's share of its
   bounding rectangle per click, and report the median — it should land near 72%.

**Comparison of `score` against a recorded CSV must cast the reference to `float32` first** — the recorded files
store `score` at float32 repr, so `np.array_equal(mine, ref.astype(np.float32))` passes while a float64 comparison
at `rtol=1e-12` fails on rounding alone. `od` is stored at full precision and compares at `rtol=1e-12`.

## Statistics

Copy the audited machinery from `shape_vs_mass_template_terms_49roi_3seed.ipynb` verbatim, including its two
self-validations (convolution against brute-force enumeration; convolution reproducing the recorded bbox3way
`hem_bbox - default_51` exact p). Do not reimplement it.

- **Unit: the ROI (49).** A ROI's delta at K = 10 is its TP@10 **summed over its 6 clicks**.
- Pooled precision delta, cluster bootstrap 95% CI over ROIs (`BOOT_SEED = 20260917`, `N_BOOT = 10_000`; assert
  your ROI order equals the recorded runs' so the resampled ROI sets are the same ones), a t interval over the
  ROI-level deltas, ROI-level W/L/T, `n_nonzero`, `min_attainable_p`, exact sign-flip p by convolution, Holm
  within family, and the D5 bar (bootstrap CI excludes 0 **and** a majority of all 49 ROIs move in the delta's
  direction).
- **Pairs. Each (axis, family) is its own Holm family — four families, 10 rows in total.**
  - *primary, 3 tests per axis* — the question as posed, our method against the incumbent:
    - `hem_rect_fill - default_51`
    - `hem_rect_fill_1um - default_51`
    - `hem_mask_fill_1um - default_51`
  - *secondary, 2 tests per axis* — which ingredient, if any, is doing the work:
    - `hem_rect_fill_1um - hem_rect_fill` — does the 1 um margin add anything? (This is the contrast that was
      K = 10-only at 3 clicks: +0.75 pts [+0.27, +1.36], but 39 of 49 ROIs tied and D5 was not met.)
    - `hem_mask_fill_1um - hem_rect_fill_1um` — does dropping the rectangle's corners add anything? Both arms
      carry the same 1 um padding and the same bounding geometry, so this isolates keep-region **shape** with
      nothing else moving.

### Reading rules, fix them before the run

- CI excludes 0, delta > 0: `a` has higher precision than `b`. CI excludes 0, delta < 0: lower.
- CI includes 0: no difference detected; quote the CI as the bound. **No equivalence is claimed.**
- State every verdict with its t CI, exact p, Holm p, ROI W/L/T, `n_nonzero` and D5 bar alongside — never a CI
  on its own.

### The tie report is a required output, not a nicety

The entire point of this run is granularity. For every pair, print **ties / 49, `n_nonzero`, and
`min_attainable_p`**, and state them **against the 3-click baseline**: at K = 10 the 3-click run had 17/49 ties at
`tm_score` (32 nonzero) and 31/49 at `chromatin_od` (18 nonzero). Flag any row whose Holm-adjusted floor
(`min_attainable_p * family size`) exceeds 0.05 — such a row could not have reached significance at any effect
size, so it bounds nothing and is not a null.

## Outputs

Write to `results/precision_at_k_49roi_6seed_fillgeom_postD11_{clicks,per_run,top30,per_roi,delta_per_roi,delta_stats,verification}.csv`
and figures alongside the notebook. **Assert your output paths collide with no recorded stem** — never overwrite
`..._bbox3way*`, `..._hemrectfill*`, `..._fillalpha*` or `precision_at_k_14roi_*`.

## How to work

1. Read `production_hematoxylin_only/shape_vs_mass_template_terms_49roi_3seed.ipynb` first — it is the closest
   template: dual-axis ranking from one search, every gate above in working form, the statistics machinery, the
   figure style. `fill_value_alpha_sweep_chromatin_od_49roi_3seed.ipynb` is the second reference, and
   `hem_rect_fill_1um_margin_vs_default_chromatin_od_49roi_3seed.ipynb` holds the exact `margin_rect` and
   `fill_outside` implementations you must reproduce.
2. **Smoke-test on 3 ROIs first** via a `SMOKE_ROIS` constant and a throwaway copy of the notebook. The gates fire
   on real data there and will catch most mistakes in ~5 minutes rather than ~60. Delete the smoke copy and its
   outputs afterwards.
3. Then the full run: 49 ROIs x 6 clicks x 4 arms = 1176 searches, roughly 60-75 minutes on this machine (CPU
   only, no CUDA). The notebook should persist its CSVs after every ROI so a crash leaves partial results.
4. Write the closing "Reading the result" section from the run's own numbers, against the reading rules you fixed
   beforehand. Name which outcome happened; do not re-choose the criteria afterwards.

## Repo conventions

- Function docstrings: brief summary, one line per parameter (name, type, short description), one line for the
  return. Google-style.
- `def` and call signatures on a single line, never wrapped, however many parameters.
- Match the surrounding notebooks' comment density and idiom.

## What would make this run worthless

- Redrawing the clicks instead of extending them (seeds 0-2 must reproduce, or the parity gates are meaningless).
- Changing the fill rule, the window size, the channel, or the peak/NMS settings.
- Replacing `hem_mask_fill_1um`'s mask with its bounding rectangle at any point. That is provably identical to
  `hem_rect_fill_1um`, so the arm would cost a full run and measure nothing.
- Reporting a K = 10 null without its tie count and `min_attainable_p`.
- Reading the two ranking axes as confirming each other. They rank the same detections; one is not a replication
  of the other.
- Quoting the `chromatin_od` result as the headline. The mass term in these templates is nearly the same statistic
  as the `chromatin_od` ranker, so that axis flatters any arm carrying one. `tm_score` is the cleaner axis here;
  report both, and say so.
