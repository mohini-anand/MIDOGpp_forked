# Window sweep: does a wider mean-centring window help at K = 10?

Build and run `production_hematoxylin_only/window_sweep_51_61_71_49roi_6seed.ipynb`.

Reason independently. This prompt fixes the **design**, not the conclusion. If a check fails, if a value here
contradicts the repo, or if the result contradicts the hypothesis, stop and report it — do not adapt the
design to make anything pass. A null is a result.

## What this sweep actually manipulates — read before designing anything

For a fill template that equals `t_i` on a keep set of `k` pixels and a constant `f_W` elsewhere,
`cv2.matchTemplate(..., TM_CCOEFF)` reduces exactly to

```
R_W(p) = Σ_keep (t_i − f_W)(I_i − μ_W(p))
       = Σ_keep (t_i − m_k)·I_i  +  (m_k − f_W)·k·(μ_keep(p) − μ_W(p))
```

`μ_W(p)` is the plain box mean of the image over the W×W window; `m_k` is the template mean over the keep
region. **The first term does not depend on W.** So widening the window does not add "flattened surround":
the filled pixels enter only through the scalar `f_W` and through the *width of the background-mean box*.
The mass weight `(m_k − f_W)` moves only ~4% from W=51 to W=71 while the filled share moves ~28%.

**Do not call this a "ring width" or "dose" manipulation.** It manipulates the mean-centring window.

**It is not inert.** Pilot on 2 ROIs × 6 clicks (`337.tiff` small nuclei, `509.tiff` large), fill arm, keep
region frozen: TP@10 changed on **9 of 12 clicks**, mean top-10 overlap between W=51 and W=71 was **5.9 of
10**, net **+4 TP of 120 slots** — and the direction split by nucleus size, `337.tiff` net −2, `509.tiff` net
+6. Two ROIs is a hint, not a result; it is quoted so you know the run is worth doing and what size of effect
to expect.

## Questions

Q1 and Q2 are asked on **both** ranking axes, `tm_score` and `chromatin_od`, as separate outcomes.

**Q1.** Does widening the mean-centring window raise precision@10 for a fill template?
(`fill_61 − fill_51`, `fill_71 − fill_51`)

**Q2.** Is any gain specific to the fill template, or does an untouched window gain the same?
Difference-in-differences: `(fill_W − fill_51) − (default_W − default_51)`.

**Q3.** Does the benefit depend on nucleus size? **Descriptive only.**

**Q4.** Is the change driven by `μ_W` or by the `f_W` drift? (`fill_71_pinned − fill_71`) **Descriptive only.**

## The seven arms

A 3 × 2 factorial plus one mechanism probe. For `W ∈ {51, 61, 71}`, cut a `W × W` `hematoxylin_od` window
centred on the click:

| arm | template |
|---|---|
| `default_W` | the window, untouched |
| `fill_W` | keep the pixels inside the Otsu rectangle; set every other pixel to the mean of those other pixels, **recomputed on this W** |
| `fill_71_pinned` | as `fill_71` but with the fill value pinned to `fill_51`'s `f_51`, holding the mass weight exactly |

7 arms × 294 clicks = 2058 searches; both ranking axes come from one search each.

## Inputs — no click selection is run

Everything needed is already in the audited run's CSVs. Do **not** import `rerun_bbox3way_postD11`, rebuild
the ROI table, or redraw clicks.

- `../results/precision_at_k_49roi_6seed_fillgeom_postD11_clicks.csv` — 294 rows. Use `file_name`,
  `seed_index`, `seed_ann_id`, `click_cx`, `click_cy`, `mpp`, `subset`, `domain`, and the Otsu rectangle
  **`rect_y0, rect_y1, rect_x0, rect_x1`** (half-open, 51-window-local). The rectangle is read, not
  recomputed — `ss.tighten_box_otsu` is never called.
- Ground truth: `images, annotations = ds.load_annotations('../databases/MIDOG++.json')`, then
  `ds.image_annotations(annotations, file_name)`.
- ROI images: `{'original_14': '../images/extra_valid', 'testing_35': '../images/extra_valid/testing_set'}[subset] + '/' + file_name`.
- ROI order for the bootstrap: the order `file_name` first appears in
  `..._fillgeom_postD11_per_run.csv`. The seed only reproduces the recorded CIs in that order.

In a W window the keep mask is `mask[rect_y0+off : rect_y1+off, rect_x0+off : rect_x1+off] = True` with
`off = (W - 51) // 2`.

## Hyperparameters

| name | value |
|---|---|
| channel | `hematoxylin_od` |
| `TM_METHOD` | `4` (`cv2.TM_CCOEFF`, unnormalised) |
| `PEAK_MIN_DISTANCE` | `7` |
| `DEEP_FLOOR_Z` | `-1.5` — **inert**, the threshold never binds |
| `MAX_PEAKS` | `100` — **binds on 100% of rows** |
| `OD_WINDOW` | `51`, against the **unblanked** channel |
| blanked square | **51 for every arm, no exceptions** |
| NMS radius = match radius | `ev.radius_px(mpp)` (7.5 µm) |
| `TEST_K` | `10` — the only budget tested; 20 and 30 recorded, quoted in no verdict |
| `N_BOOT` / `BOOT_SEED` | `10000` / `20260917` |
| max window | **W ≤ 73** (`BORDER = tm.PATCH_SIZE // 2 = 36`). Do not exceed |

`MAX_PEAKS` binds everywhere, so the pool is a hard top-100 truncation and `DEEP_FLOOR_Z` never fires — do
not attribute anything to the threshold. A wider template gives a smoother map and may concentrate peaks, so
**report NMS suppression (`n_peaks − n_detections`) per arm.**

## Pipeline, per (click, arm)

1. `hem = ch.to_channel(rgb, 'hematoxylin_od')`
2. template window: `tm.read_padded_patch(hem, cx, cy, W)`
3. `search, n_blanked = tm.blank_seed_square(hem, cx, cy, 51)`
4. pad `search` and the **unblanked** `hem` by `PAD = (W-1)//2`, `cv2.BORDER_REPLICATE`
5. `fused, _, valid = tm.fused_response(search_padded, [template], False, method=4)`
6. crop `fused`/`valid` to `[PAD:PAD+h, PAD:PAD+w]`; assert `valid.all()`
7. `med, mad = tm.robust_stats(fused, valid)`; `threshold = med + (-1.5) * mad`
8. `centers, scores = tm.extract_peaks(fused, valid, 7, threshold, 100)`
9. `keep = nms_by_distance(centers, scores, ev.radius_px(mpp))` (`midog_utils/nms.py`)
10. `od` via `cm.score_detections(detections shifted by PAD, od_padded, window=51)`
11. rank **once per axis from this one detection list**: `tm_score → score`, `chromatin_od → od`, descending,
    `kind='mergesort'`, `na_position='last'`
12. `gt_eval = gt[gt['ann_id'] != seed_ann_id]`; `ev.bucket_detections(...)`; TP@K = count of
    `ev.HUMAN_CORRECT_LABEL` in the top K

`od` is bit-identical for PAD 25/30/35 (verified), so the wider pad is safe.

## Correctness — one gate, four inline asserts, two diagnostics

**The gate: W=51 parity.** `default_51` and `fill_51` must reproduce the recorded run row for row on **both**
axes. This single check proves the clicks, the rectangle, the fill value, the blanking and the whole pipeline
are right — which is why nothing else in this spec re-derives them.

Reference: `../results/precision_at_k_49roi_6seed_fillgeom_postD11_per_run.csv` (588 rows per condition,
conditions `default_51` and `hem_rect_fill`) and `_top30.csv` (17 640 rows per condition). `_top30.csv` is
keyed by `(file_name, seed_index, condition, rank_key)` — **include `rank_key` in the join** or the two axes
duplicate. Compare `tp_at_10/20/30`, `n_detections`, `tpl_sha1`, and top-30
`rank, cx, cy, bucket, matched_ann_id`; `score` after an `np.float32` cast, `od` at `rtol=1e-12`.

**Four one-line asserts in the loop**, each catching a bug parity cannot see because parity only covers W=51:

1. `n_blanked_px == 2601` on every row of every arm.
2. `n_filled_px == W*W - n_kept_px` and `fill_value == np.float32(T_W[~mask].astype(np.float64).mean())`.
   Without this, a bug where `fill_61`/`fill_71` reuse `f_51` passes everything else and silently turns the
   run into `fill_W_pinned` — indistinguishable at the 4% level.
3. the W-local keep mask, cropped back to its central 51×51, equals the 51-local mask (catches a wrong offset).
4. `n_detections >= 10`.

**Two reported diagnostics** — report, do not halt:

- **Near-seed.** The blank is filled with `channel.min()`, which lowers `μ_W` and therefore *inflates* `R` for
  tissue in the window. The keep region clears the 51 px blank at ~39 px, but the window overlaps it out to
  50 px at W=51, 55 at W=61, 60 at W=71 — so the bias runs in the direction the hypothesis predicts, and a
  match-radius check cannot see it. With the blank fixed at 51 by decision, **this is the only instrument
  that can detect it.** Report per arm: the count of top-10 detections within 60 px of the click that are
  **not** true positives. On the recorded W=51 run this is 0 — every top-10 detection within 60 px is a
  genuine mitosis. If it rises on the wide arms, say so prominently and treat every `fill_61`/`fill_71` gain
  as unresolved. Do not net it out yourself.
- **ROI edge.** `PAD` grows 25 → 30 → 35, so the band whose template overlaps replicate-padded fabricated
  content grows with W. In the recorded run 3.68% of top-10 slots lie within 25 px of an edge and 4.70%
  within 35 px, and that band's TP rate is **0.246 against 0.560 outside it** — a band growing by 1.5 points
  of slots moves precision ≈ **−0.47 points**, a quarter of the effect, against the wide arms. Add a column
  `tp_at_10_excl_border` computed from the same detection lists with detections within 35 px of any ROI edge
  dropped identically in every arm. Zero extra search cost; it does not touch the primary numbers.

## Statistics

Unit of analysis is the **ROI (49)**. A ROI's value is its **TP@10 summed over its 6 clicks** (out of 60);
pooled precision uses denominator `6 × 49 × 10 = 2940`.

Per pair: pooled delta in points, 95% cluster-bootstrap CI over ROIs (`N_BOOT=10000`, `BOOT_SEED=20260917`,
the ROI order above), a t interval over ROI deltas, ROI W/L/T, `n_nonzero`, exact sign-flip p by convolution,
`min_attainable_p = 2 ** (1 - n_nonzero)`, Holm within family, and the D5 bar (**CI excludes 0 AND more than
24 of 49 ROIs move in the delta's direction**).

**Four Holm families, 8 tests.** Each (axis, family) is its own family, matching the parent run's convention:

| family | tests |
|---|---|
| `tm_score \| window` | `fill_61 − fill_51`, `fill_71 − fill_51` |
| `tm_score \| fill_vs_raw` | `(fill_61−fill_51) − (default_61−default_51)`, same at 71 |
| `chromatin_od \| window` | `fill_61 − fill_51`, `fill_71 − fill_51` |
| `chromatin_od \| fill_vs_raw` | the two difference-in-differences |

The difference-in-differences is what answers Q2 — compute it per ROI from the four per-ROI numbers. **Do
not** answer Q2 by eyeballing two separately-tested effects.

Both are production ranking axes (`production.py`: `tm_score` default, `chromatin_od` opt-in). **They are two
outcomes, not two looks at one** — one search per (click, arm), differing only in the sort key. A result on
one axis is **not** a replication of the other. Never pool them.

**Descriptive, in no family, supporting no verdict:** `default_W − default_51`; `fill_71_pinned − fill_71`;
`fill_W − default_W`; `tp_at_10_excl_border`; K = 20 and K = 30; the axis comparison (there is no
prespecified threshold for "larger at `tm_score`", so do not declare a mechanism confirmed from it).

Prior for `default_W − default_51`: on `feature/annotation-frame`, "precision tracks footprint size" —
`prod_51` is 3.8 points **below** `default_51` at `tm_score` (Holm 0.001, D5 met). That predicts a *positive*
`default_W − default_51`, which is exactly why Q2 needs the difference-in-differences.

## Q3 — descriptive, and bounded

Because the keep region is frozen, 51 → 71 adds exactly **+10 px uniformly to every click**. A flat
association with nucleus size is therefore equally consistent with "one window fits all" and with "size
matters but is locally linear here". **This design cannot separate those, and cannot identify a saturation
point from 3 window values. Do not claim it does, and do not turn a non-significant association into a design
decision.**

Report two things and stop:

1. The per-ROI `fill_71 − fill_51` against ROI median nucleus size (`max(rect_h, rect_w)`), as a scatter plus
   a Spearman with a CI.
2. The same association computed **within** ROIs (ROI-demeaned click-level), which is free of the species
   confound.

State these two facts beside them: Spearman(ROI median nucleus px, canine indicator) = **−0.69, p = 3.7e-8**,
so a between-ROI size analysis is largely a species contrast; and the parent run's own ROI-level
Spearman(nucleus size, `hem_rect_fill − default_51`) is **−0.322 (p = 0.024)** at `tm_score` and **+0.158
(p = 0.28)** at `chromatin_od` — that is the baseline, not the +0.76 seven-domain-mean figure.

Worth watching: the parent's fill gain runs **negative** in nucleus size (small nuclei gain more), while the
pilot has the wider window helping the **large**-nucleus ROI. Different manipulations, so not a
contradiction — but if the full run reproduces that split it is the most informative thing here.

## Do not

- **Do not analyse by domain.** n = 7 ROIs/domain; the Holm floor across 7 domains is
  `2**(1-7) × 7 = 0.109 > 0.05`, so no domain can reach corrected significance at any effect size.
- Do not call this a "ring width" or "dose" manipulation.
- Do not redraw clicks or recompute the Otsu rectangle.
- Do not let the blanked square follow `W`, anywhere, for any arm.
- Do not exceed W = 73.
- Do not quote K = 20 or K = 30 in any verdict.

## Outputs

`../results/precision_at_k_49roi_6seed_windowsweep_postD11_{per_run,per_roi,delta_per_roi,delta_stats,top30}.csv`
— note `../`, since `nbconvert` sets the working directory to the notebook's own. `per_run` carries
`fill_value`, `n_filled_px`, `n_kept_px`, `n_peaks`, `n_detections`, `n_blanked_px`, `tp_at_10/20/30`,
`tp_at_10_excl_border`, and the near-seed count.

Every number in the conclusion must be printed by a cell above it.

## Environment

- Interpreter: `/Users/mohinianand/anaconda3/bin/python3` (3.11.5, cv2 4.8.1, numpy 1.26.4, pandas 2.2.2).
  No other interpreter here can import `cv2`.
- Execute: `/Users/mohinianand/anaconda3/bin/jupyter nbconvert --to notebook --execute --inplace
  --ExecutePreprocessor.timeout=7200 <path>`.
- Runtime: `matchTemplate` uses DFT and gets *faster* with a larger template (measured 3.02 / 2.56 / 1.98 s at
  W = 51 / 61 / 71). The parent's 1176 searches took 35.9 min; 2058 searches here cost the same or less each,
  so budget **~50–60 min**. CPU only, no CUDA.

## Code conventions

Google-style docstrings: brief summary, one line per parameter (type + short description), one line for the
return. Function definitions and calls on a **single line**, never wrapped.
