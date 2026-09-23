# Annotation DataFrame — the pipeline's input and output format

Status: design, agreed 2026-09-23. Not yet implemented. Supersedes the output format implied by
`CORNER_COORDINATE_REFACTOR_PROMPT.md`; see
`Research Logs/2026-09-23-corner-coordinate-refactor-prompt-audit.md` for why that prompt's
centre-based design was rejected.

One frame carries every box in a click-to-verify session — what the human drew, what the algorithm
refined it to, and what the algorithm proposes. The pipeline takes this frame in and returns a new
one with rows appended. Evaluation is a separate function over that output.

---

## 1. Two coordinate layers, and the rounding boundary between them

The half-pixel bugs this project has hit all come from one cause: mixing a continuous box with a
pixel-index box. This spec keeps them apart.

**Layer 1 — the DataFrame (continuous).** Every box is `x_tl, y_tl, width, height` as float64.
The centre is always

    cx = x_tl + width / 2      cy = y_tl + height / 2

There is **no rounding anywhere in this layer.** `width` is an extent, `x2 - x1`, matching
`dataset.py:68` and MIDOG's own arithmetic ("Scheme A").

**Layer 2 — pixel footprints (integer indices).** `read_padded_patch`, `blank_seed_square` and
`build_augmentations` address pixels. A footprint of odd side `base_size` centred on a continuous
centre covers, inclusively,

    ix = round(cx)   half = base_size // 2   pixels ix-half .. ix+half        (exactly base_size of them)

**Rounding happens exactly once, at the Layer 1 → Layer 2 boundary, through a shared helper.**
Nothing else in the codebase may call `round` on a coordinate.

### The `+1` that must never be forgotten

A footprint of `base_size = 51` pixels is written in Layer 1 as `width = 50.0`:

    width = base_size - 1        base_size = width + 1

This is the price of the seed row round-tripping MIDOG byte-exactly, and it is deliberate. State it
in the docstring of every function that converts between the two.

### Why this is safe (measured, not argued)

Under `cx = x_tl + width / 2`, all three row kinds recover their centre exactly:

| row | `x_tl` | `width` | recovered centre |
| --- | --- | --- | --- |
| human seed | MIDOG `x1`, verbatim (`.0` or `.5`) | 50.0 | `x1 + 25` — exact |
| refined seed | `round(tx) - base_size//2` (integer) | `base_size - 1` | `round(tx)` — exact |
| prediction | `cx - 25.0` (cx integer from `extract_peaks`) | 50.0 | `cx` — exact |

Facts behind this, re-derived from `databases/MIDOG++.json` (26,286 annotations) and
`results/precision_at_k_49roi_3seed_chromatin_bbox3way_postD11_per_run.csv` (441 runs):

- 2,328 MIDOG `x1` values are half-integers → `x_tl` **must** be float64, never int.
- 169 boxes have a negative top-left (min −30.0) → **never clamp** coordinates to the image.
- `w = x2 - x1 = 50` and `h = 50` on all 26,286 rows.
- 134 of 441 runs have a half-integer template centre `tx` → the Layer 2 rounding is live, not
  theoretical.

### The one thing the refined row does not preserve

The footprint is centred on `round(tx)`, so the unrounded analytic centre `tx` is **not recoverable
from the frame**. Click-to-template `offset_px` must be read from `info`, never recomputed from the
DataFrame. The half-open Otsu component bbox never appears in the frame at all — it stays inside
`tightened_template_box`, which confines the half-open/closed distinction to one function.

---

## 2. Schema

Column order is the frame's column order.

| column | dtype | human seed | refined seed | prediction |
| --- | --- | --- | --- | --- |
| `file_name` | str | ROI file name | same | same |
| `x_tl`, `y_tl` | float64 | MIDOG `x1, y1` verbatim | footprint TL (integer-valued) | `cx - 25.0`, `cy - 25.0` |
| `width`, `height` | float64 | 50.0 | `base_size - 1` | 50.0 |
| `class` | str | `MITOTIC_GT` | `MITOTIC_GT` | `PREDICTION` |
| `object_id` | int64 | MIDOG `ann_id` | **same** `ann_id` | 0, 1, 2, … by rank |
| `made_by` | str | `HUMAN` | `ALGORITHM` | `ALGORITHM` |
| `refined` | bool | `False` | `True`, or `False` on gate fallback (§4) | `False` |
| `current_time` | str, ISO-8601 UTC | click time | refinement finished | prediction stage finished |
| `score` | float64 | `NaN` | `NaN` | fused template-match score |
| `od` | float64 | `NaN` | `NaN` | chromatin optical density |

**Conventions.**

- `class` stays `MITOTIC_GT` on the refined row: it is the same object. `made_by` and `refined`
  carry the provenance.
- `object_id` is an **entity** id, not a row key. Seed and refined seed deliberately share it — that
  is the link between them. The row key is `(file_name, object_id, made_by, refined)`. Prediction
  ids are unique only within one ROI.
- `made_by` is uppercase in both values.
- `refined` is a real bool, never the strings `"TRUE"` / `"FALSE"`.
- `current_time` is one value per stage, repeated across that stage's rows. All three stage values
  differ within a run.
- `score` and `od` are **both always computed**, whatever ranks the list. `rank_key` chooses the
  sort order only. `info["od_computed"]` is therefore always `True`.
- Row order: human seed, refined seed, then predictions by rank ascending. Deterministic, so diffs
  are readable.
- Coordinates are full-resolution ROI pixels, the same frame as MIDOG ground truth. The internal
  `border_pad` and `OD_WINDOW // 2` padding must not leak into emitted coordinates.
- Read persisted frames back with `pd.read_csv(..., float_precision='round_trip')`, as
  `rerun_bbox3way_postD11.py:334` already does, or the half-integer coordinates drift.

Per-run metadata — `mpp`, `rank_key`, `refinement_channel`, `base_size`, counts, timings — lives in
`info`, not as a column repeated on every row. When concatenating runs across conditions, add the
discriminating column at concat time, the pattern `rerun_bbox3way_postD11.py:252` uses.

---

## 3. The production function

    run_roi(rgb, annotations_df, mpp, *, seed_object_id=None, rank_key="tm_score",
            refinement_channel="gray_inverted") -> tuple[pd.DataFrame, dict]

- **Returns a new frame.** Input rows are copied through unchanged; new rows are appended. The
  caller's frame is never mutated.
- **Never draws a seed.** In production the human already clicked. The validation path uses
  `draw_seed_row` (§6) to produce the seed row, then calls this.
- `seed_object_id` selects which input row is the seed. `None` means infer: exactly one row with
  `made_by == 'HUMAN'`; raise naming the count if zero or more than one.
- Emits **all** post-NMS survivors, not a truncated top-N. Evaluation takes `head(K)`.

**Row-count contract.** The returned frame has exactly `2 + n_detections` rows: the human seed, the
refined seed, and one row per post-NMS survivor. Measured across the 441 reference runs,
`n_detections` was 82–100 (median 99), so a frame is typically ~101 rows.

Post-NMS survivors are the only detections that reach the frame. The pre-NMS peak count `n_peaks` is
capped at `max_peaks = 100` (D9) inside `extract_peaks`, before NMS runs, and never becomes rows.

### What the pipeline consumes from the human row

Only its centre, `x_tl + width/2`. `width` and `height` are accepted and ignored — the template size
comes from the Otsu gate or the 51 px default. This is a deliberate, documented contract, not an
oversight: the human's own box is preserved in the frame whether or not the matcher uses it, so
switching from clicks to hand-drawn boxes later loses no recorded information.

`build_seed` and `tightened_template_box` take `(x_tl, y_tl, width, height)` and derive the centre
once through the shared helper. `read_padded_patch` keeps taking a literal point — its window is
inherently symmetric about a pixel, which is irreducible.

---

## 4. The Otsu gate fallback

In production the human clicks once and there is no pool to redraw from, so the current
`build_seed` behaviour (redraw on refusal, raise when the pool empties) cannot apply.

When the click's own Otsu component fails the gate — centre not foreground, or the component fails
`min_area` / `max_area_frac` / `min_solidity` — fall back to the default `BASE_SIZE = 51` template
centred on the click, and emit the refined-seed row with **`refined = False`**. The row is still
`made_by = 'ALGORITHM'` (the algorithm produced that box), and `info` records the refusal reason.

This is what the `refined` column is for; it must not be hardcoded `True`.

---

## 5. `info`

Returned alongside the frame, always. Carries everything the frame deliberately does not.

- **Identity / config**: `file_name`, `mpp`, `match_radius_px`, `nms_radius`, `rank_key`,
  `refinement_channel`, `channel`, `tm_method`, `base_size`, `patch_size`, `peak_min_distance`,
  `deep_floor_z`, `max_peaks`.
- **Seed**: `seed_object_id`, click centre, `refined`, gate refusal reason when applicable,
  `template_center` (the **unrounded** `tx, ty`), `offset_px`, `n_blanked_px`.
- **Search**: `n_augmentations`, `deep_floor_median`, `deep_floor_mad`, `score_threshold_used`,
  `n_peaks`, `max_peaks_binding`, `n_after_nms`, `n_detections`.
- **Timings**: `t_match_s`, `t_postprocess_s`.

`max_peaks_binding` means the peak pool hit the D9 cap of 100 **before** NMS ran. Measured on the
reference CSV: `n_peaks == 100` and `max_peaks_binding` is `True` on **all 441 runs**. At production
config it is therefore a constant, not an anomaly flag — a `True` is not a fault. What it records is
that every run's candidate pool is already a capped subset, so the post-NMS survivors in the frame
are never the complete detection set for the ROI. It must still survive into `info`, since it stops
being constant the moment `max_peaks` changes.

`n_blanked_px` must equal `base_size ** 2` exactly — verified true on all 441 runs of the reference
CSV, across `base_size` 19…51. It is the sharpest available check that the Layer 1 → Layer 2
conversion is right, because a wrong rounding changes it even when every centre still round-trips.

---

## 6. Evaluation

    evaluate(result_df, gt_df, mpp, roi_shape) -> tuple[pd.DataFrame, pd.DataFrame, dict]

Separate from `run_roi`. There is no `production_only` flag: the caller evaluates or does not.

It needs three things the production frame deliberately does not carry:

- **`gt_df`** — the ROI's full ground truth, mitotic **and** look-alike. Bucketing is three-way
  (`human_correct_label` / `human_rejected_label` / `non_human_findings`), so look-alikes are
  required for `lookalike_attraction`. A human in production never annotates these, so they do not
  belong in the production frame.
- **`mpp`** — the match radius is `7.5 µm / mpp` (D7), not derivable from the frame.
- **`roi_shape`** — for `coverage_fraction`.

Rules it must follow:

- Score only `class == 'PREDICTION'` rows. The two seed rows are never candidates.
- Exclude the seed's own annotation from ground truth before matching, by `object_id == ann_id`.
  Forgetting this inflates recall.
- Match against the **human** seed's geometry, never the refined one — D8: ground truth stays keyed
  to the click.
- Matching is by centre distance on `cx = x_tl + width/2`, unchanged from `evaluate.greedy_match`.
  Box overlap is not used anywhere; the match radius is a physical circle (D7).

---

## 7. Validation

Reference: `results/precision_at_k_49roi_3seed_chromatin_bbox3way_postD11_per_run.csv` —
`tp_at_10`, `tp_at_20`, `tp_at_30` over 49 ROIs × 3 seeds × 3 conditions, `rank_key='chromatin_od'`.

Two requirements, both easy to get wrong:

1. **The seed draw must match the reference's draw exactly.** That CSV is produced by a *joint*
   gate: a click is accepted only if all three conditions' Otsu gates pass it
   (`draw_joint_click`). A plain `build_seed` draw selects different clicks and the numbers will not
   match, for reasons unrelated to this refactor. `draw_seed_row` must reproduce the joint gate for
   validation runs.
2. **Exclude `current_time` from every regression comparison, by name.** It differs on every run by
   design. Nothing else in the frame may differ.

---

## 8. Deliberately deferred

- Multi-loop sessions. The schema does not block them: a second loop appends more `HUMAN` rows. No
  accept/reject verdict column is defined yet — when it is needed, a `HUMAN` row whose `object_id`
  matches a `PREDICTION` row is the natural encoding.
- Using the human's drawn `width`/`height` to size the template.
- Emitting the tight Otsu component bbox as its own row.
