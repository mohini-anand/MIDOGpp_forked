# Annotation DataFrame — the pipeline's input and output format

Status: implemented. The design was agreed 2026-09-23 and amended the same day after review — every
amendment is listed, and tagged as a correction or a design decision, in §9. The implementation
(`2c0d875`) passes every §7 check: 65 checks, 0 divergences, outputs in
`validation/annotation_frame/out/`. `click_frame` and `production_pipeline/run_frame.py` were added
on 2026-09-24 (`41ffa2e`, §3). Later the same day `run_frame.py` became the pipeline's single entry
point (`--source click|box|midog`), and a drawn box became a seed source whose size sets the Otsu
gate's window and the fallback template (`box_frame`, `gate_window_side`; §3, §4; `DECISIONS.md`
D12). Clicks and MIDOG rows, all 50 wide, keep the 51 px window and fallback. The full §7
harness was re-run after that change: 0 divergences on all 49 ROIs (960,825 values) and on the
named checks, every regenerated per-ROI result identical to the committed one apart from timings;
the only changed named output is the edge-error wording. Box-specific checks are in
`validation/annotation_frame/box_checks.py` (`out/box/box.json`). The
hand-off prompt is `ANNOTATION_FRAME_IMPLEMENTATION_PROMPT.md`; where the two disagree, this spec
governs.

Supersedes `CORNER_COORDINATE_REFACTOR_PROMPT.md`. That prompt changed the internal point
representation of the existing `midog_utils` functions. This spec leaves the existing functions'
behaviour unchanged — one of them gains one keyword-only argument whose default is today's
behaviour — and adds the frame as a new layer on top of them (§0).
`Research Logs/2026-09-23-corner-coordinate-refactor-prompt-audit.md` audits the old prompt.

One frame carries every box in a click-to-verify session — what the human drew, what the algorithm
refined it to, and what the algorithm proposes. The pipeline takes this frame in and returns a new
one with rows appended. Evaluation is a separate function over that output.

---

## 0. Layering: a new layer on top of frozen code

The experiments run in this repository since 2026-09-12 reach below the production entry point
rather than through it. Checked at the call sites:

- `fill_geometry_precision_at_10_49roi_6seed.ipynb` and `shape_vs_mass_template_terms_49roi_3seed.ipynb`
  copy `find_and_suppress`'s body with a custom template, calling `tm.blank_seed_square(hem, cx, cy,
  base_size)`, `tm.fused_response`, `tm.extract_peaks`, `nms_by_distance`, `cm.score_detections` and
  `ev.bucket_detections` on `cx`/`cy` frames. Both check that copy against
  `prod.run_production_pipeline(rgb, SimpleNamespace(base_size=..., template_xy=(cx, cy)), mpp, ...)`.
  `NMS_RADIUS_ABLATION_HEMRECTFILL_PROMPT.md` prescribes the same pattern.
- The hybrid gate prototype (`hybrid_gray_hem_tightening_14roi.ipynb`) calls
  `ss.tightened_template_box` (unpacking its 3-tuple), `ss.tighten_box_otsu` and `ss.border_filter`.
- `fill_geometry_precision_at_10_49roi_6seed.ipynb` imports `rerun_bbox3way_postD11` as a module for
  its click draw.

So the frame is a new module on top, and the layers underneath keep working exactly as today.

**Frozen** — signatures, defaults, return types, column names and numerics all unchanged:

- every file in `midog_utils/` other than `production.py`, including `__init__.py`,
  `FSConfig`, `Seed`, `build_seed`, and `tightened_template_box`'s `(base_size, tx, ty)` return;
- `production.py` apart from the one addition below;
- `rerun_bbox3way_postD11.py`, `production_pipeline/run_pipeline.py` and `midog_utils_full/`.
  (`run_pipeline.py` carries a deprecation note in its docstring since 2026-09-24; its code is
  unchanged, which `ast_equal.py` proves.)
  `rerun_bbox3way_postD11.main()` writes tracked CSVs that are this spec's validation reference
  (§7), so validation imports that script's functions and never calls `main()`.

**The one change to existing code.** `production.run_production_pipeline` gains a keyword-only
argument `always_od: bool = False`. When `True`, the `od` column is computed by production's own
`chromatin_od` block (`production.py:56-61`: same channel, pad and window) whatever `rank_key` is,
and the sort still uses `rank_key`'s column. `info["od_computed"]` becomes
`rank_key == "chromatin_od" or always_od`. With the default, every existing caller gets
bit-identical output. That is proven in §7, not assumed.

**New.** One module, `midog_utils/annotation_frame.py`, holds everything this spec defines: the
frame constants and converters, `click_frame`, `box_frame`, `gate_window_side` and `run_roi` (§3,
§4), `evaluate` (§6) and
`draw_seed_rows` (§7). It calls the frozen `midog_utils` functions and never re-implements them.
The one thing it ports is the driver's joint-gate click draw, because a package module cannot
import a repo-root script; §7 defines the port and requires it to reproduce the driver exactly. The
module is not added to `midog_utils/__init__.py`; import it as
`from midog_utils import annotation_frame`. `production_pipeline/run_frame.py` is the pipeline's
single entry point (§3). `run_pipeline.py` is deprecated, code unchanged, until the collaborator
repo has the new command line.

Research that needs a hook `run_roi` deliberately lacks — a template override, a custom gate, an
`nms_radius` or `max_peaks` override, an augmented bank, per-detection angle/flip/scale — keeps
calling the frozen layers directly, as it does today.

---

## 1. Two coordinate layers, and the rounding boundary between them

The half-pixel bugs this project has hit all come from one cause: mixing a continuous box with a
pixel-index box. This spec keeps them apart.

**Layer 1 — the DataFrame (continuous).** Every box is `x_tl, y_tl, width, height` as float64.
The centre is always

    cx = x_tl + width / 2      cy = y_tl + height / 2

There is **no rounding anywhere in this layer.** `width` is an extent, `x2 - x1`, matching
`dataset.py:67-68` and MIDOG's own arithmetic ("Scheme A").

**Where the origin is.** A coordinate is an array index: the value `k` is the *centre* of pixel
column `k` (row `k` for y), so pixel `k` spans `[k - 0.5, k + 0.5]`. This is the convention every
function in `midog_utils` already assumes (`read_padded_patch` slices the array at `round(cx)`).
A box's edges `x_tl` and `x_tl + width` therefore fall on the centres of its outermost pixels.
Drawn over an image, it outlines those pixels only at `x_tl - 0.5 .. x_tl + width + 0.5`, which is
why `viz.draw_box` subtracts 0.5. A client whose coordinates are pixel *edges* (pixel `k` spans
`[k, k + 1)`, as in COCO and most canvas APIs) must subtract 0.5 before building a row;
`click_frame(..., edge_origin=True)` does this (§3). The frame never stores edge coordinates.

**Layer 2 — pixel footprints (integer indices).** `read_padded_patch`, `blank_seed_square` and
`build_augmentations` address pixels. A footprint of odd side `base_size` centred on a continuous
centre covers, inclusively,

    ix = round(cx)   half = base_size // 2   pixels ix-half .. ix+half        (exactly base_size of them)

**The Layer 1 → Layer 2 conversion stays where it already is**, inside the frozen pixel-level
functions (§0):

- `border_filter`, via `np.rint`;
- `_patch_readable`, `read_padded_patch`, `tightened_template_box` and `blank_seed_square`, via
  `int(round(·))`.

Both are **round-half-to-even**, and that matters. D8's centre formula gives a half-integer
template centre on at least one axis in 207 of the 441 reference runs. On 103 of those,
round-half-up would pick a different pixel. `D8_TEMPLATE_ANCHOR.md` describes this tie behaviour.

New code passes continuous centres into the frozen functions and never rounds a coordinate itself,
with two exceptions:

- the refined-seed row's top-left (table below). It must be computed as
  `int(round(tx)) - base_size // 2`, the same expression `read_padded_patch` uses. Never use
  `math.floor(x + 0.5)` or any other round-half-up form for a centre-to-pixel conversion.
- `box_frame`'s snap of a drawn box (§3), which happens before the row exists: each side becomes
  `round(w)` pixels raised to odd, centred on the pixel `floor(x0 + w / 2)` that contains the drawn
  centre (`x0` a pixel edge). In the frame's centre coordinates that is `floor(c + 0.5)`, the
  round-half-up form banned above, and the carve-out is deliberate: the ban protects conversions
  that must agree with `read_padded_patch`; this one picks the pixel containing an edge-coordinate
  point, has no tie, does not depend on where the box sits, and yields an integer centre, so the
  frozen functions' half-to-even rounding never meets a tie on a box row. `round()` there would make
  the stored window depend on the box's position.

`gate_window_side` also rounds, but a size, not a coordinate: `round(max(width, height)) + 1`
pixels, raised to odd (§4).

### The `+1` that must never be forgotten

A footprint of `base_size = 51` pixels is written in Layer 1 as `width = 50.0`:

    width = base_size - 1        base_size = width + 1

This is the price of the seed row round-tripping MIDOG byte-exactly, and it is deliberate. State it
in the docstring of every function that converts between the two.

The identity holds only for **odd** `base_size`. Every size that reaches it is odd: `_odd()` in
`seed_selection` and `template_match` forces it, and the fallback is `gate_window_side`'s window,
odd by construction (51 for clicks and MIDOG rows). Converters must still
assert it, because for an even size, `round(tx) ± base_size // 2` covers `base_size + 1` pixels.

### Why this is safe (measured, not argued)

Under `cx = x_tl + width / 2`, all three row kinds recover their centre exactly:

| row | `x_tl` | `width` | recovered centre |
| --- | --- | --- | --- |
| human seed | MIDOG `x1`, computed as `cx - w / 2` from `load_annotations` (`.0` or `.5`) | `w` = 50.0 | `x1 + 25` — exact |
| refined seed | `int(round(tx)) - base_size // 2` (integer) | `base_size - 1` | `round(tx)` — exact |
| prediction | `cx - BOX_SIZE / 2`, i.e. `cx - 25.0` (cx integer from `extract_peaks`) | `BOX_SIZE` = 50.0 | `cx` — exact |

Facts behind this, re-derived from `databases/MIDOG++.json` (26,286 annotations) and
`results/precision_at_k_49roi_3seed_chromatin_bbox3way_postD11_per_run.csv` (441 runs):

- Every bbox value has fractional part `.0` or `.5`. For all 26,286 rows, both
  `cx - w / 2 == x1` and `(cx - w / 2) + w / 2 == cx` hold bit-exactly. The human row can therefore
  be built from `load_annotations`, which stays frozen and has no `x1` column.
- 2,328 MIDOG `x1` values are half-integers → `x_tl` **must** be float64, never int.
- The JSON mixes types: 95,842 bbox values are Python ints and 9,302 are floats. 87 of the 503
  annotated images have only int boxes.
  - `load_annotations`' `cx` and `cy` are float regardless, because of the `/ 2.0`.
  - Its `w` and `h` are float64 only because one column holds every image.
  - A frame built from the raw JSON of one of those 87 images comes out int64, so every converter
    casts the coordinate columns to float64 explicitly.
- 169 boxes have a negative top-left (min −30.0), and 105 have `x2 > W` or `y2 > H`. Two
  ground-truth centres lie outside the image (305.tiff ann 15363 at `cx = −2`; 309.tiff ann 15728
  at `cy = −5`) → **never clamp** coordinates to the image.
- `w = x2 - x1 = 50` and `h = 50` on all 26,286 rows.
- 134 of 441 runs have a half-integer `tx`, and 207 have one in `tx` or `ty` → the Layer 2
  rounding is live, not theoretical.

**Exactness comes from dyadic widths, not from the formula.** MIDOG's half-pixel grid makes every
value above exact.

- A future client sending arbitrary sub-pixel clicks with `width = 50` also round-trips exactly in
  the range `run_roi` accepts. Of 100,000 uniform clicks in `[25, 7215)`, 0 fail
  `(c - 25) + 25 == c`. The failures all sit at `c < 25`, below the 36 px border step 2
  enforces.
- The risk is a non-dyadic width, such as a hand-built row. With `width = 47.3`, about 230 in
  100,000 clicks fail `(c - w/2) + w/2 == c`; with `33.1`, about 60. `box_frame` snaps a drawn box
  to whole pixels, so its row has an integer width and centre and round-trips exactly.
- The error is one ulp, so it can only change a round-half-to-even result for a click within one
  ulp of a `.5` tie.
- Tests that assert bit-exact round trips should use clicks on the half-pixel grid and
  `width = 50`.

### The one thing the refined row does not preserve

The footprint is centred on `round(tx)`, so the unrounded analytic centre `tx` is **not recoverable
from the frame**. `info["template_center"]` carries the unrounded `(tx, ty)` and
`info["offset_px"]` carries the click-to-template offset (§5). Anything measured from the template
centre — `offset_px`, the driver's `d_tpl`, `tpl_cx`/`tpl_cy` — must be read from `info`, never
recomputed from the DataFrame. The half-open Otsu component bbox never appears in the frame at
all. It stays inside the frozen `tightened_template_box`, which confines the half-open/closed
distinction to that one function.

---

## 2. Schema

Column order is the frame's column order.

| column | dtype | human seed | refined seed | prediction |
| --- | --- | --- | --- | --- |
| `file_name` | str | ROI file name | same | same |
| `x_tl`, `y_tl` | float64 | MIDOG `x1, y1` verbatim, as `cx - w / 2`, `cy - h / 2`; a click: `cx - 25.0`; a drawn box: snapped by `box_frame` (§3), integer-valued | footprint TL `int(round(tx)) - base_size // 2` (integer-valued) | `cx - 25.0`, `cy - 25.0` |
| `width`, `height` | float64 | `w`, `h` (50.0; a click: 50.0; a drawn box: its odd pixel count − 1) | `base_size - 1` | `BOX_SIZE` (50.0) |
| `class` | str | `MITOTIC_GT` | `MITOTIC_GT` | `PREDICTION` |
| `object_id` | int64 | MIDOG `ann_id` in validation; in production, any integer unique among the ROI's `HUMAN` `MITOTIC_GT` rows | **same** as the human row | 0, 1, 2, … = production's `rank` |
| `made_by` | str | `HUMAN` | `ALGORITHM` | `ALGORITHM` |
| `refined` | bool | `False` | `True` when the gate accepted; `False` on fallback or when no refinement was requested (§4) | `False` |
| `current_time` | str, ISO-8601 UTC, microseconds | click time (validation: when `draw_seed_rows` built the row) | refinement finished | prediction stage finished |
| `score` | float64 | `NaN` | `NaN` | fused template-match score, upcast from production's float32 |
| `od` | float64 | `NaN` | `NaN` | chromatin optical density |

**Conventions.**

- `class` stays `MITOTIC_GT` on the refined row: it is the same object. `made_by` and `refined`
  carry the provenance.
- `object_id` is an **entity** id, not a row key, and it is scoped by `class`. An entity is
  `(file_name, class, object_id)`: ground truth and predictions are separate id spaces that share
  one column.
  - Seed and refined seed deliberately share their `(class, object_id)` — that is the link between
    them.
  - A prediction whose rank happens to equal the seed's id is a different entity, because its
    class differs. Without `class` in the key, the two would collide. On a gate fallback, for
    instance, the refined row would have the same `(object_id, made_by, refined)` as that
    prediction. In validation this is latent (the smallest seed `ann_id` among the reference
    clicks is 130; ids ≤ 100 exist only in images 1–6). In production it is near-certain for any
    scheme that numbers clicks from 0 or 1.
  - The row key is `(file_name, class, object_id, made_by, refined)`, unique within one run's
    frame.
  - Prediction ids are unique only within one run. The same physical detection gets different ids
    under the two rank keys, so join frames from different runs on coordinates, not on
    `object_id`. Coordinates are unique among one run's predictions, because NMS never keeps two
    detections closer than its radius.
- `made_by` is uppercase in both values.
- `refined` is a real bool, never the strings `"TRUE"` / `"FALSE"`.
- `class` is a Python keyword, so read it as `df["class"]`, never `df.class`.
- `current_time` has one value per stage, repeated across that stage's rows. Write it as
  `datetime.now(timezone.utc).isoformat(timespec="microseconds")`. Of the two stamps `run_roi`
  writes itself, the refinement stamp is never later than the prediction stamp. Input rows keep
  whatever time the caller gave them. Stamps are not guaranteed to be distinct, and every
  regression comparison excludes the column by name (§7).
- `score` and `od` are **both always computed**, whatever ranks the list. `rank_key` chooses the
  sort order only. `info["od_computed"]` is therefore always `True`.
- Row order: the input rows in input order, then the refined seed, then the predictions in exactly
  the order `run_production_pipeline` returned them; `object_id` is that position, which equals
  production's `rank`. `run_roi` never re-sorts. Production breaks ties with a stable mergesort
  over the NMS order, which is score descending, then x ascending, then y ascending (the lexsort in
  `extract_peaks`). The reference top-30s already contain 3 exact `od` ties, and any re-sort could
  reorder them.
- Coordinates are full-resolution ROI pixels, the same frame as MIDOG ground truth. The internal
  `border_pad` and `OD_WINDOW // 2` padding must not leak into emitted coordinates.
- Read persisted frames back with `pd.read_csv(..., float_precision='round_trip')`, as
  `rerun_bbox3way_postD11.py:334` already does. Half-integer coordinates are safe either way. The
  columns at risk are `od` and the other long decimals (`mpp`, radii): with the default parser,
  9,052 of the 13,230 `od` values in the committed top30 CSV read back different in the last bits.
- `score` is float32 inside production and float64 in the frame. The upcast is lossless, but a
  CSV written from the float32 column holds float32 text, e.g. `4.7094283`; the committed top30 CSV
  does. Compare scores against such a file only after casting both sides to float32. Compared as
  float64, 0 of 13,230 rows match.

Per-run metadata — `mpp`, `rank_key`, `refinement_channel`, `base_size`, counts, timings — lives in
`info`, not as a column repeated on every row. When concatenating runs across conditions, add the
discriminating column at concat time, the pattern `rerun_bbox3way_postD11.py:252` uses.

---

## 3. The production function

    run_roi(rgb, annotations_df, mpp, *, seed_object_id=None, rank_key="tm_score",
            refinement_channel="gray_inverted") -> tuple[pd.DataFrame, dict]

- `rgb` is the ROI's full-resolution RGB array (`dataset.load_roi`) and `mpp` its microns per pixel
  (`dataset.roi_mpp`).
- `rank_key` must be a key of `production.AXES` (`"tm_score"` or `"chromatin_od"`).
  `refinement_channel` must be `None`, `"gray_inverted"` or `"hematoxylin_od"`. Anything else
  raises before any work is done.
- **Returns a new frame.** Input rows are copied through with their values unchanged; only the
  coordinate columns are cast to float64, which is lossless. New rows are appended. The caller's
  frame is never mutated.
- **Never draws a seed.** In production the human already clicked, and `click_frame` builds that
  click's row (below). The validation path uses `draw_seed_rows` (§7) to produce the seed row, then
  calls this.
- **Seed selection.** If `seed_object_id` is given, the seed is the one input row with
  `class == "MITOTIC_GT"`, `made_by == "HUMAN"` and that `object_id`. If it is `None`, the seed is
  the one input row with `class == "MITOTIC_GT"` and `made_by == "HUMAN"`. Zero or several
  matches raise, naming the count.
- **Input rows.** Every input row must share the seed's `file_name`, and must be
  `made_by == "HUMAN"` with `class == "MITOTIC_GT"`. Their `object_id`s must be unique. Every
  row's box must give a valid gate window (`gate_window_side`, §4: finite, 19–73 px), checked after
  the column and class checks. Anything else raises. Feeding `run_roi` its own output back would
  duplicate the row key of every prediction row, and multi-loop input is deferred (§8).
- Emits **all** post-NMS survivors, not a truncated top-N. Evaluation takes `head(K)`.

**Algorithm.** Every numerical step is delegated to frozen code (§0).

1. Compute the gate window `gate_window_side(width, height)` (§4), and
   `cx = x_tl + width / 2` and `cy = y_tl + height / 2`, from the seed row, in float64.
2. **Border.** If `seed_selection._patch_readable(rgb.shape, cx, cy, tm.PATCH_SIZE)` is `False`,
   raise `ValueError`.
   - This is the same predicate as `border_filter`'s 36 px margin. Every existing draw applies
     `border_filter` to the click before gating it (`build_seed`, the joint gate), and this keeps
     `run_roi` consistent with that.
   - It also guarantees that the fallback template and the gate window, both centred on the seed
     and at most 73 px (36 px margin), are readable.
   - It does **not** mean nothing could be read. A refined template can move inward far enough to
     be readable. On 300.tiff, the click `(25.0, 290.0)` gates on gray to `(33, 37.0, 281.0)`, and
     `(25.0, 715.0)` gates on hematoxylin to `(31, 37.0, 718.5)`; both templates are readable.
     Such clicks are refused anyway, by design (§4, §8).
3. **Template.** Choose the template spec `(base_size, tx, ty)` and set `refined` as §4 specifies.
4. **Search.** Call
   `prod.run_production_pipeline(rgb, SimpleNamespace(base_size=base_size, template_xy=(tx, ty)), mpp, rank_key=rank_key, always_od=(rank_key != "chromatin_od"))`.
   - Production reads only those two seed fields. The fill-geometry and shape-vs-mass parity
     checks already pass the same duck-typed object.
   - `always_od` is set only where it changes anything. `chromatin_od` computes `od` anyway, so
     passing the default there keeps every `chromatin_od` run on the unmodified default path, the
     one existing callers use. §7's V3 exercises that path on all 441 runs.
5. **Rows.** Append the refined-seed row (§1 table) and one prediction row per returned
   detection, in returned order. Each prediction row gets `x_tl = cx - BOX_SIZE / 2` and
   `y_tl = cy - BOX_SIZE / 2`, from production's integer `cx` and `cy`, and `object_id = rank`.
   Stamp both stages.
6. **Info.** Assemble `info` (§5).

For `refinement_channel="hematoxylin_od"`, the hematoxylin channel is computed twice: once here for
the gate, and once inside production, whose frozen signature takes `rgb`. That costs seconds per
call. It is accepted to keep production frozen.

**Row-count contract.** The returned frame has exactly `len(annotations_df) + 1 + n_detections`
rows: the input rows, the refined seed, and one row per post-NMS survivor. In the single-click case
that is `2 + n_detections`; the single-click case is the validation path and today's production.
Across the 441 reference runs, `n_detections` was 82–100 (median 99), so a frame typically has
about 101 rows.

Post-NMS survivors are the only detections that reach the frame. The pre-NMS peak count `n_peaks` is
capped at `max_peaks = 100` (D9) inside `extract_peaks`, before NMS runs, and never becomes rows.

### What the pipeline consumes from the human row

Its centre, `x_tl + width/2`, and its size: `gate_window_side(width, height)` is the Otsu gate's
window and the fallback template's side (§4). Every click and MIDOG row is 50 wide, so window 51,
exactly the behaviour before 2026-09-24. A drawn box's row is stored snapped (`box_frame`, below):
it is exactly the gate window (for a square box) and the fallback footprint, but the raw drawn box
is not recoverable from the frame — neither its sub-pixel coordinates nor the difference between an
even drawn size and the next odd one (decision, user 2026-09-24). A client that needs the raw box
(e.g. for a later radius-based refinement) must keep it itself.

`run_roi` derives the click centre once, in step 1, and hands it to `_patch_readable` and
`tightened_template_box`. The template centre `(tx, ty)` from step 3 goes to production. There
`read_padded_patch` and `blank_seed_square` receive it; `(tx, ty)` is the click itself on fallback
and when no refinement is requested. Each of these frozen functions does its own
round-half-to-even (§1). `build_seed` is not used: it draws from a pool, and `run_roi` never
draws. It stays in the package for `run_pipeline.py` and for research.

### Building the input row, and the command line

`run_roi` loads nothing itself. Its caller loads `rgb` and `mpp` with `dataset.load_roi` and
`dataset.roi_mpp`, as `run_pipeline.py` does, and builds the one-row input frame with one of four
functions:

| source | builder | `object_id` | scored by `run_frame.py` |
| --- | --- | --- | --- |
| a raw click (production) | `click_frame(file_name, cx, cy, object_id, *, edge_origin=False)` | the caller's id | no |
| a drawn box (production) | `box_frame(file_name, x0, y0, w, h, object_id)` | the caller's id | no |
| a MIDOG annotation | `human_seed_frame(ann, file_name)` | its `ann_id` | yes |
| the validation draw | `draw_seed_rows(rgb, gt, image_id, *, seed_indices)` (§7.1) | the drawn `ann_id` | yes |

- `click_frame` writes a `CLICK_BOX = 50.0` box centred on the click: `x_tl = cx - 25.0`,
  `width = 50.0`, `class = MITOTIC_GT`, `made_by = HUMAN`, `refined = False`, `score` and `od` NaN.
  With `edge_origin=True` it subtracts 0.5 from both coordinates first (§1).
- At a MIDOG annotation's centre, with `object_id = ann_id`, it returns exactly
  `human_seed_frame`'s row. Checked on all 26,286 annotations, and with `edge_origin=True` at the
  centre plus 0.5.
- `evaluate` accepts a `click_frame` row only when its `object_id` is the mitotic `ann_id` at that
  exact centre (§6). A production click's id is arbitrary, so `run_frame.py` does not score it.
- `box_frame` takes a drawn box in pixel-edge coordinates (left, top, width, height, as in COCO
  and canvas APIs) and snaps each side (§1): `n = round(w)` pixels, raised to odd; centre pixel
  `c = floor(x0 + w / 2)`; `x_tl = c - n // 2`, `width = n - 1`. The stored size depends only on
  `w`, the position only on the drawn centre. It raises unless the snapped box gives a 19–73 px
  window (a drawn box of 18–73 px on its longer side). Not for MIDOG boxes, which are already in
  pixel-centre coordinates.
- `box_frame(x0, y0, 51, 51)` equals `click_frame(x0 + 25, y0 + 25)` (all columns but
  `current_time`) for **integer** origins. Elsewhere the box snaps to the pixel containing its
  centre while a click keeps its sub-pixel centre, e.g. `x0 = 11.5`: box centre 37, click 36.5
  (which the frozen functions round to 36). No equivalence is claimed there.

`production_pipeline/run_frame.py` is the pipeline's single entry point (since 2026-09-24;
`run_pipeline.py` is deprecated):

- `--source {click,box,midog}` is required and chooses the builder. `--source click` takes
  `--click X Y` (pixel-centre coordinates unless `--edge-origin`) and `--object-id`; `--source box`
  takes `--box X Y W H` (pixel-edge coordinates) and `--object-id`; `--source midog` takes
  `--ann-id`, or `--seed-index k` (default 0) for the validation draw. Another source's option, or a
  missing required one, is an argparse error before anything loads, and so is an invalid `--box`.
  `--click-id` was renamed `--object-id`.
- Click and box seeds are not scored; MIDOG seeds are. `--save-fig` and `--plot` draw the overview:
  ground truth, the top-N predictions, the seed and template centre, the gate window (white, when a
  gate ran) and the template footprint (cyan), with a full-resolution zoom on the seed.
- For `--seed-index k` it draws `range(k + 1)` and keeps the last. Draws are sequential (§7.1), so
  drawing `(k,)` alone gives a different click for any `k > 0`.
- The joint-gate draw is not `build_seed`'s, so the same `--seed-index` can pick a different click
  in `run_pipeline.py`. On 300.tiff, seed 0 is 14581 there and 14499 here: 14581 passes the gray
  gate but not the hematoxylin one.
- It finds an ROI in `--images-dir` or in its `testing_set/` subfolder, and `--save-csv` writes
  the whole frame.

Checked 2026-09-24:

- `run_frame_on_roi` reproduces R1's `seed_ann_id`, `n_detections` and `tp_at_{10,20,30}` on 18
  runs: 300.tiff and 322.tiff, seeds 0–2, all three conditions, `rank_key="chromatin_od"`.
- On 300.tiff click 14499, `--ann-id`, `--click` and `--edge-origin` give identical frames.
- The full §7 harness was not re-run for that addition; it was re-run after the box seed and
  `--source` changes (status line above).
- Only `rank_key="chromatin_od"` (R1, R2) and `default_51` under `tm_score` (R4) have a committed
  reference. `run_frame.py`'s defaults, `tm_score` with the gray gate, are covered by B2 and V4
  only.

---

## 4. Template choice and the gate fallback

In production the human clicks once and there is no pool to redraw from, so the current
`build_seed` behaviour (redraw on refusal, raise when the pool empties) cannot apply.

By this point `run_roi`'s step 2 has already confirmed the click is at least 36 px
(`tm.PATCH_SIZE // 2`) from every edge; otherwise it raised. Below,
`channel = ch.to_channel(rgb, refinement_channel)` and `W = gate_window_side(width, height)` of the
seed row: the pixel count of the box's longer side, `round(max(width, height)) + 1`, raised to odd.
`W` is 51 for every click and MIDOG row (width 50), and must lie in 19–73 (`MIN_WINDOW_PX`,
`MAX_WINDOW_PX = tm.PATCH_SIZE`), or the input raises. Every remaining case is one of these:

| case | template `(base_size, tx, ty)` | `refined` | `info["refusal_reason"]` |
| --- | --- | --- | --- |
| `refinement_channel is None` | `(W, cx, cy)` | `False` | `None` |
| the gate accepts: `tightened_template_box(channel, cx, cy, otsu_window=W)` returns `(b, tx, ty)` and `_patch_readable(rgb.shape, tx, ty, 73)` holds | `(b, tx, ty)` | `True` | `None` |
| `tightened_template_box` returns `None` | `(W, cx, cy)` | `False` | `"gate_refused"` |
| it returns `(b, tx, ty)`, but the template centre is within 36 px of an edge | `(W, cx, cy)` | `False` | `"template_unreadable"` |

- For a non-square box the window and the fallback are the enclosing square on the longer side,
  centred on the box centre. The gate's window is always square.
- `W` bounds the window and the fallback, not the accepted template: the gate can return templates
  below 19 px (its own floor is 9 px in effect: `min_area = 50` needs a bounding box of at least
  8 px, made odd), as it could at 51.
- The window is not in `info`: it is recoverable from the frame's human row with
  `gate_window_side`, and the fallback size is `info["base_size"]`.

- The first case is the untightened `default_51` arm. It is the control in nearly every recent
  experiment and one third of the validation reference (147 of 441 runs).
  `info["refinement_requested"]` tells it apart from a fallback.
- After step 2, any window up to 73 px (which needs at most a 36 px margin) is readable, and
  `_template_choice` checks it, raising `RuntimeError` otherwise, so this holds by construction.
  So `tightened_template_box` returning `None` always means the gate refused: the click's own pixel
  is not foreground, or its component fails `min_area` / `max_area_frac` / `min_solidity`.
  `tighten_box_otsu` is frozen and does not say which. Finer reasons are deferred (§8).
- The fallback template is the human's box: `W` px, centred on the unrounded seed centre, with
  its blanked square there. For clicks and MIDOG rows that is the `default_51` template exactly. So
  on a refused seed, `run_roi` with a refinement channel must return exactly the same predictions
  as `run_roi(refinement_channel=None)`. §7 tests this for clicks, `box_checks.py` for boxes.
- These cases reproduce the driver's `template_spec` (`rerun_bbox3way_postD11.py:79-94`). Its
  `None` for a tightened condition covers both refusal rows. Its hematoxylin readability check for
  `default_51` is the step-2 predicate.
- A click within 36 px of the edge raises, as the border filter in `build_seed` and in the joint
  gate already does. A client must catch the `ValueError` and tell the user.

The refined-seed row is always emitted, with `made_by = 'ALGORITHM'`, because the algorithm chose
that box. `refined` says whether the box came from an accepted component. This is what the
`refined` column is for; it must never be hardcoded `True`.

---

## 5. `info`

Returned alongside the frame, always. Carries everything the frame deliberately does not.

**Everything production returns.** `info` starts as a copy of `run_production_pipeline`'s info
dict, with every key kept unchanged:

- `nms_radius`, which `find_and_suppress` rounds to 2 decimals. Use `match_radius_px` for any
  computation.
- `n_augmentations`, `n_blanked_px`, `pad_px`, `t_match_s`, `deep_floor_median`,
  `deep_floor_mad`, `score_threshold_used`, `n_peaks`, `n_after_nms`, `n_detections`,
  `max_detection_score`, `t_postprocess_s`, `max_peaks_binding`, `od_computed`, `rank_key`.

**Added by `run_roi`.** None of these may overwrite a production key.

- **Identity / config.**
  - `file_name`, `mpp`, and `match_radius_px = ev.radius_px(mpp)` (unrounded).
  - `refinement_channel`, and `refinement_requested = refinement_channel is not None`.
  - Production's own constants: `channel = prod.CHANNEL`, `tm_method = prod.TM_METHOD`,
    `peak_min_distance = prod.PEAK_MIN_DISTANCE`, `deep_floor_z = prod.DEEP_FLOOR_Z`,
    `max_peaks = prod.MAX_PEAKS`, `patch_size = tm.PATCH_SIZE`.
  - `base_size`, the template side actually used.
- **Seed.**
  - `seed_object_id`.
  - `click_center`, the `(cx, cy)` from step 1.
  - `refined`.
  - `refusal_reason`: one of the values in §4's table, which is the only place they are defined.
  - `template_center`, the **unrounded** `(tx, ty)` that was passed to production.
  - `offset_px = hypot(tx - cx, ty - cy)`, unrounded. The driver's `tpl_offset_px` is this
    rounded to 3 decimals.
- **Timings.** `t_refine_s` for step 3 and `t_total_s` for the whole call, alongside production's
  `t_match_s` and `t_postprocess_s`. Every timing key starts with `t_`, and every comparison
  excludes timing keys.

`max_peaks_binding` means the peak pool hit the D9 cap of 100 **before** NMS ran. Measured on the
reference CSV: `n_peaks == 100` and `max_peaks_binding` is `True` on **all 441 runs**. At production
config it is therefore a constant, not an anomaly flag — a `True` is not a fault. What it records is
that every run's candidate pool is already a capped subset, so the post-NMS survivors in the frame
are never the complete detection set for the ROI. It must still survive into `info`, since it stops
being constant the moment `max_peaks` changes.

`n_blanked_px` must equal `base_size ** 2` exactly — verified true on all 441 runs of the reference
CSV, across `base_size` 19…51. It checks the blanked footprint's **size**, not its position.

- A centre shifted by a pixel still blanks `base_size ** 2` pixels. `blank_seed_square` at
  `cx = 100.0, 100.5, 101.0, 101.5` blanks 2,601 pixels each time, in columns 75–125, 75–125,
  76–126 and 77–127.
- So it catches size errors, like the old `size / 2.0` bug that blanked 53 × 53, and nothing else.
- Position errors are caught instead by the template fingerprint `tpl_sha1` (§7 V1) and by
  prediction-level equality (V3, V6).

---

## 6. Evaluation

    evaluate(result_df, gt_df, mpp, roi_shape, *, seed_object_id=None, budgets=(10, 20, 30))
        -> tuple[pd.DataFrame, pd.DataFrame, dict]

Separate from `run_roi`. There is no `production_only` flag: the caller evaluates or does not.

It needs three things the production frame deliberately does not carry:

- **`gt_df`** — the ROI's full ground truth, mitotic **and** look-alike, in the frozen
  `dataset.load_annotations` format: the ROI's rows of its annotations frame, e.g.
  `ds.image_annotations(annotations, file_name)`.
  - It is not in frame format. Bucketing is three-way (`human_correct_label` /
    `human_rejected_label` / `non_human_findings`), so it needs look-alikes, and the frame defines
    no look-alike class. `recall_by_agreement` needs the vote columns.
  - A human in production never annotates look-alikes, so they do not belong in the production
    frame.
- **`mpp`** — the match radius is `7.5 µm / mpp` (D7), not derivable from the frame.
- **`roi_shape`** — `(h, w)` or `(h, w, c)`, for `coverage_fraction` and the area.

Steps:

1. **Seed.** Find the seed row with `run_roi`'s selection rule (§3), using the `seed_object_id`
   argument only to choose among rows. Set `seed_id = int(seed_row["object_id"])`.
   - `evaluate` is for MIDOG ground truth only, and it requires `seed_id` to be the seed's MIDOG
     `ann_id`. `draw_seed_rows` guarantees that. A production click's id is arbitrary (§2), and
     there is no ground truth to evaluate against anyway.
   - Assert it: `gt_df` must hold exactly one row with `ann_id == seed_id`, with
     `category_id == ds.MITOTIC`, and with `(cx, cy)` exactly equal to the seed row's centre.
     Otherwise raise.
   - Without that check, a frame whose clicks are numbered 1, 2, … would, on images 1–6, delete
     an unrelated real annotation and keep the clicked one.
2. **Exclude the seed's own annotation.**
   `gt_eval = gt_df[gt_df["ann_id"] != seed_id].reset_index(drop=True)`, exactly as the driver
   builds it (`rerun_bbox3way_postD11.py:235`). Forgetting this inflates recall.
   - Never filter on the `seed_object_id` argument. It is `None` in the default single-click
     call, and `ann_id != None` is `True` on every row, so nothing would be excluded.
   - Exclude by the seed row's id only — never by every `object_id` in the frame. Prediction ids
     run 0…~99, so doing that would silently delete real ground truth whose `ann_id` is small.
3. **Candidates.** Take the rows with `class == "PREDICTION"`, stable-sorted by `object_id`
   ascending. Greedy matching depends on rank order, so row order is never trusted. The seed rows
   are never candidates.
4. **Convert.** Build the centre-based frame the frozen evaluator reads: `rank = object_id`,
   `cx = x_tl + width / 2`, `cy = y_tl + height / 2`, `score`, `od`.
5. **Score.** `ev.evaluate_run(detections, gt_eval, ev.radius_px(mpp), area_mm2, roi_shape)`, with
   `area_mm2 = h * w * mpp * mpp / 1e6`. That is the exact expression of
   `midog_utils_full.dataset.roi_area_mm2` (`:197`). Don't write `mpp ** 2`: it differs in the
   last bit on 29 of the 49 ROIs, although `evaluate_run` rounds the reported value to 2 decimals.
6. **Budgets.** Add `tp_at_{k}` to the metrics for each `k` in `budgets`: the number of
   `human_correct_label` buckets among the first `k` candidates. That is exactly the driver's
   `score_top_k`. Also add `seed_object_id = seed_id`.

Returns:

- the candidate rows in frame format, sorted, with `bucket`, `matched_ann_id` and
  `matched_category` added;
- `evaluate_run`'s `gt_out`. Its `matched_rank` is a position among the sorted candidates. That
  equals the prediction's `object_id` only when the candidates are the complete, contiguous
  `0…n−1` that `run_roi` emitted; `evaluate` asserts this and raises otherwise;
- the metrics dict. Its `match_radius_px` is `evaluate_run`'s, rounded to 1 decimal; the
  unrounded radius is in `run_roi`'s `info`.

Rules:

- Ground truth stays keyed to the click (D8). The click enters evaluation only through step 2's
  exclusion, and `greedy_match` never sees seed geometry. Diagnostics measured from the click —
  the driver's `d_click`, `n_near_click` and `nearest_other_ann_px` — use the human row's centre
  (`info["click_center"]`), never the refined row. Diagnostics measured from the template
  (`d_tpl`) use `info["template_center"]`.
- Matching is by centre distance on `cx = x_tl + width/2`, unchanged from `evaluate.greedy_match`.
  Box overlap is not used anywhere; the match radius is a physical circle (D7).

---

## 7. Validation

### 7.1 The validation draw: `draw_seed_rows`

    draw_seed_rows(rgb, gt, image_id, *, seed_indices=(0, 1, 2)) -> list[tuple[pd.DataFrame, dict]]

The reference was drawn through a **joint** gate: a click is accepted only if every condition can
build a template from it. A plain `build_seed` draw picks different clicks, and the numbers then
fail to match for reasons that have nothing to do with this code.

`draw_seed_rows` is a port of the draw in `rerun_bbox3way_postD11.select_clicks` (`:161-194`) and
its helpers `template_spec`, `pool_gate_specs` and `draw_joint_click`. It is plural because draws
are sequential: each seed excludes the clicks earlier seeds accepted. One call therefore produces
all of an ROI's seeds.

The draw picks rows **by position**, so pool order is load-bearing. The port must do exactly this,
in this order, with the frozen functions:

1. `gt` is the ROI's annotations as `ds.image_annotations(annotations, file_name)` returns them,
   all categories, in that row order. Assert that `gt` holds one `file_name` and one `image_id`,
   and that the latter equals the `image_id` argument. That argument seeds the RNG, so a mismatch
   would silently draw a different stream.
2. Convert `gray = ch.to_channel(rgb, "gray_inverted")` and `hem = ch.to_channel(rgb, "hematoxylin_od")`.
3. Build the pool:
   - `pool, flagged = ss.agreement_pool(gt[gt["category_id"] == ds.MITOTIC])`;
   - `n_agreement_pool = len(pool)`;
   - `pool = ss.border_filter(pool, tm.PATCH_SIZE // 2, rgb.shape)`.
4. For each pool row, in pool order, compute one spec per condition in
   `JOINT_CONDITIONS = {"default_51": None, "gray_bbox": "gray_inverted", "hem_bbox": "hematoxylin_od"}`.
   These are the driver's names, in the driver's order; they appear verbatim in `refused_draws`.
   Each spec is computed exactly as `template_spec` does:
   - start from `(51, cx, cy)` when the condition's channel is `None`, and from
     `tightened_template_box(channel, cx, cy, otsu_window=51)` otherwise;
   - the spec is `None` if that is `None` or if `read_padded_patch(hem, tx, ty, 73)` is `None`;
   - otherwise it is `(int(b), float(tx), float(ty))`.
5. `n_joint_valid` is the number of rows whose specs are all non-`None`. Raise `ValueError` if it is
   below `len(seed_indices)`.
6. For each `seed_index` in `seed_indices`, in order:
   - `seed_pool` is the pool minus the rows whose `ann_id` an earlier seed accepted;
   - `rng = np.random.default_rng([seed_index, image_id])`, and `working = seed_pool.copy()`;
   - loop:
     - `idx = int(rng.integers(len(working)))` and `row = working.iloc[idx]`;
     - accept `row` if all its specs are non-`None`;
     - otherwise append `f"{int(row['ann_id'])}:{'+'.join(names of the None conditions, in JOINT_CONDITIONS order)}"`
       to the refusals and drop `working.index[idx]`. The `int()` matters, as in the driver
       (`:141`): a row taken from an all-numeric frame would otherwise print `22673.0`.
   - raise if `working` empties.
7. Return, per seed:
   - the one-row human-seed frame, as specified in §2, with `file_name` taken from `gt`;
   - `draw_info` with `seed_index`, `seed_ann_id`, `n_agreement_pool`, `n_after_border` (the
     border-filtered pool size), `n_joint_valid`, `n_seed_pool`, `n_retries` (the number of
     refusals), `refused_draws` (the refusals joined by `";"`), `contested_seed_tier` (`flagged`)
     and `specs` (condition → spec).

### 7.2 What must be shown

**References.** These are tracked files; nothing may ever write to them.

- **R1** — `results/precision_at_k_49roi_3seed_chromatin_bbox3way_postD11_per_run.csv`: 441 runs
  = 49 ROIs × 3 seeds × 3 conditions, `rank_key="chromatin_od"`.
- **R2** — `results/precision_at_k_49roi_3seed_chromatin_bbox3way_postD11_top30.csv`: each run's
  top 30 detections (13,230 rows).
- **R3** — `results/precision_at_k_49roi_3seed_chromatin_bbox3way_per_run.csv`: the pre-D11
  selection that `rerun_bbox3way_postD11.compare_selection` gates against.
- **R4** — from the fill-geometry run: the `condition == "default_51"` rows of
  `results/precision_at_k_49roi_6seed_fillgeom_postD11_per_run.csv` and `..._top30.csv`, plus all
  294 rows of `..._clicks.csv`, which has no `condition` column (one row per ROI × seed). This
  covers 49 ROIs × 6 seeds under both rank keys, with `od` in the `tm_score` top 30 too.
  - They were committed in a9293d0 on 2026-09-22, after the last `midog_utils` commit (9e22eb9,
    D11).
  - That notebook drew its clicks with `rr.select_clicks` at seeds 0–5, so seeds 0–2 are R1's 147
    clicks.
  - It searched with its own copy of `find_and_suppress`, which it checked against
    `run_production_pipeline`.
  - It covers the untightened template only, not gray or hematoxylin.

R4 is the committed post-D11 reference for `tm_score` that this spec uses. The shape-vs-mass run
(`results/precision_at_k_49roi_3seed_shapemass_postD11_{per_run,top30}.csv`, 0dfda5d) holds an
identical `default_51` × `tm_score` arm on the same 147 clicks: 4,410 of 4,410 top-30 rows match
on `cx`, `cy` and float32 `score`. It adds no independent coverage. The older seed-robustness CSV
predates D11, and its notebook no longer imports against the current package: it needs
`midog_utils.compare`, which exists only in `midog_utils_full/`.

`tm_score` is therefore covered by R4, by a baseline captured from the unedited code (B2 below),
and by rank-key independence (V4).

**Baselines,** captured from the unedited code before any edit, with every run's seed built
exactly as `rerun_bbox3way_postD11.run_roi` builds it:

- **B1** — every survivor and the full `info` of `run_production_pipeline(rank_key="chromatin_od")`,
  for all 441 runs.
- **B2** — the same for `rank_key="tm_score"`, on the 126 runs of the 14 `original_14` ROIs.

Before anything is edited:

- B1 must reproduce R1 and R2;
- B2's 42 `default_51` runs (14 ROIs × 3 seeds) must reproduce R4's `tm_score` top 30 on `rank`,
  `cx`, `cy` and float32 `score`, and on `rr.score_top_k`'s buckets compared with R4's `bucket`,
  `matched_ann_id` and `matched_category`. R4 has no `angle`, `flip` or `scale`, and B2 has no
  `od`.

Together these confirm the environment for both rank keys.

**Rules for every check.**

- Every comparison is exact. `score` is compared after casting both sides to float32 (§2).
- Two missing values compare equal. A float NaN equals a float NaN: 18 runs have a NaN
  `recall_contested`, because 6 reference clicks leave no contested mitotic annotation in
  `gt_eval`. String columns are compared after `.fillna("")`, as `compare_selection` already does
  (`rerun_bbox3way_postD11.py:277`). A CSV reads an empty string back as NaN: `near_click_detail`
  is NaN on 437 of R1's 441 rows and `refused_draws` on 363.
- `current_time`, and every `t_*` column or key (including R1's `t_pipeline_s`), are excluded by
  name.
- CSVs are read with `float_precision="round_trip"`.
- A divergence stops the run and is reported. It is never absorbed into a tolerance.

**Checks.**

- **V1 — the draw.**
  - `draw_seed_rows` must equal `rerun_bbox3way_postD11.select_clicks` field by field on all 49
    ROIs: accepted `ann_id`, specs, refusals and `n_seed_pool` per seed, plus `n_agreement_pool`,
    `n_after_border`, `n_joint_valid` and `flagged`.
  - With `seed_indices=range(6)`, it must reproduce R4's clicks CSV on all 294 draws:
    `seed_ann_id`, `n_agreement_pool`, `n_after_border`, `n_joint_valid`, `n_seed_pool`,
    `n_retries`, `refused_draws` and `contested_seed_tier`. This tests the sequential exclusion
    past seed 2.
  - **The contested-tier path.** `contested_seed_tier` is `False` on all 441 reference runs, so
    that path is otherwise untested. Pick an ROI whose contested tier has at least 3 joint-valid
    annotations, and remove its unanimous mitotic rows from a copy of the annotation table.
    `draw_seed_rows` on that ROI must then equal `rr.select_clicks(roi, modified_table)`, with
    `flagged` `True`. `select_clicks` takes the annotation table as an argument, so no module
    global is touched.
  - `compare_selection(rows, R3)` must pass 441/441. The rows are built from `draw_info` plus
    `run_roi`'s `info`: `base_size`, `tpl_cx`/`tpl_cy` from `template_center`, and
    `tpl_offset_px = round(offset_px, 3)`. `tpl_sha1` is
    `template_digest(hem, info["base_size"], *info["template_center"])`, which fingerprints the
    template `run_roi` actually searched with.
- **V2 — the frozen layer is unchanged.**
  - `git diff <base> -- midog_utils/` touches only `production.py`, which gains exactly the
    `always_od` argument. `git status --short midog_utils/` lists only that modification and the
    new, untracked `annotation_frame.py`; `git diff` alone does not show untracked files. Every
    other tracked file is byte-identical to the base commit.
  - After the edit, `run_production_pipeline(rank_key="tm_score")` with the default `always_od`
    must equal B2 on all 126 runs: every column, and `info` minus `t_*`.
  - With `always_od=True`, on the same 126 runs, it must return the same rows, order, `score` and
    `info` (except `od_computed`), plus an `od` column.
  - The default `chromatin_od` path needs no separate re-run. `run_roi` passes
    `always_od=False` for `chromatin_od` (§3 step 4), so V3's 441 runs go through exactly the path
    existing callers use.
- **V3 — the adapter, `chromatin_od`, all 441 runs.**
  - `run_roi`'s prediction rows must equal B1 converted by §2, on every row, not just the top 30.
  - The production part of `run_roi`'s `info` (§5) must equal B1's `info` key for key and value for
    value, minus `t_*`.
  - The frame must have `2 + n_detections` rows.
  - The refined row must equal
    `(int(round(tpl_cx)) - b // 2, int(round(tpl_cy)) - b // 2, b - 1, b - 1)`, with `b = base_size`,
    against R1's `tpl_cx`, `tpl_cy` and `base_size` on all 441 runs. This is the one new rounding
    site (§1), and 103 of these runs tell round-half-even from round-half-up.
  - Every R1 column must be recomputed and must match R1:
    - `subset`, `domain` and `image_id` come from `rr.build_rois()`'s ROI table;
    - `precision_at_{k}` is `tp_at_{k} / k`, as in the driver's `main()`;
    - everything else comes from the frame, `info`, `draw_info` and `evaluate`. That includes each
      diagnostic, computed with the driver's own formulas (`:235-251`): `n_near_click`,
      `near_click_detail`, `nearest_other_ann_px`, `od_nan` and `od_descending`.
  - R2 must be reproduced from each frame's first 30 predictions plus `evaluate`'s buckets. At
    production config R2's `angle`, `flip` and `scale` are the constants `0.0`, `False` and `1.0`;
    confirm that, then skip those columns.
- **V4 — rank-key independence, all 441 runs.**
  - `run_roi(rank_key="tm_score")` must keep the same set of prediction coordinates as
    `run_roi(rank_key="chromatin_od")`, with identical `score` and `od` at each coordinate.
  - On the 126 runs it must also equal B2 converted by §2, order included, on every frame column
    except `od`. `angle`, `flip` and `scale` are not frame columns. B2 has no `od`, because the
    unedited `tm_score` path doesn't compute it. `od` in `tm_score` order is checked by the bullet
    above and by R4, below.
  - On the 147 `default_51` runs (seeds 0–2), `run_roi(refinement_channel=None, rank_key="tm_score")`
    must reproduce R4's `tm_score` top 30: `rank`, `cx`, `cy`, `score` (as float32), `od`,
    `bucket`, `matched_ann_id` and `matched_category`. It must also reproduce R4's per-run
    `tp_at_{10,20,30}`, `n_peaks`, `n_detections` and `n_blanked_px`. These runs are already part
    of this check, so the comparison costs nothing extra.
  - This was measured on one ROI before the spec was amended (300.tiff, seed 0, `default_51` and
    `hem_bbox`): identical 100-survivor sets, and the committed top 30 reproduced exactly. That is
    a feasibility check, not a result.
- **V5 — fallback and no refinement.**
  - Take clicks one gate refused in the reference draws. These are candidates; re-verify each
    before use:
    - 462.tiff ann 22673 (gray refused);
    - 529.tiff ann 24934 (gray refused);
    - 300.tiff ann 14581 (hematoxylin refused);
    - 245.tiff ann 6274 (hematoxylin refused).
  - `run_roi` with the refused channel must return `refined=False` and
    `refusal_reason == "gate_refused"` exactly, and predictions identical to
    `run_roi(refinement_channel=None)` on the same click.
  - With the other channel it must return `refined=True`.
  - **The `template_unreadable` branch.** Use a synthetic click on 300.tiff at `(6410.0, 250.0)`
    (width 6447, so `x = w - 37`), which step 2 accepts.
    - Both gates accept it: gray gives `(51, 6414.0, 250.0)` and hematoxylin gives
      `(27, 6414.5, 257.0)`. Both template centres are too close to the right edge to read.
    - For each channel, `run_roi` must return `refined=False`,
      `refusal_reason == "template_unreadable"` exactly, and predictions identical to
      `run_roi(refinement_channel=None)` on the same click.
    - An implementation that forgets the readability check crashes inside `find_and_suppress`
      here, and one that mislabels the case fails the string check.
    - Re-verify the click before use. If the frozen code ever stops producing this case there,
      scan 300.tiff's clicks 36–40 px from its right edge for another.
- **V6 — half-integer clicks.** No reference click is half-integer (0 of 147), so pin two. These
  are candidates; re-verify before use:
  - 072.tiff ann 1835 at (3590.5, 3139.5);
  - 007.tiff ann 144 at (6922.5, 2004.0).

  On both, round-half-even and round-half-up pick different x pixels for the click (3590 vs 3591;
  6922 vs 6923). On 007.tiff they also differ on the hematoxylin template's `ty = 2004.5`. For
  each:
  - the human row's `x_tl` and `y_tl` must equal the JSON `x1` and `y1` exactly;
  - for all three refinement settings, the refined row's `x_tl` and `y_tl` must equal
    `int(round(tx)) - base_size // 2` and `int(round(ty)) - base_size // 2`;
  - `n_blanked_px` must equal `base_size ** 2`. That is a size check only (§5);
  - for all three refinement settings and both rank keys, the prediction rows must equal the §2
    conversion of `run_production_pipeline(rgb, seed, mpp, rank_key=k, always_od=(k != "chromatin_od"))`,
    called with a hand-built seed from `template_spec` — the call `run_roi` makes.
- **V7 — edges and inputs.**
  - On a real ROI, a click 35 px from an edge raises and one 36 px from it does not.
  - Zero or two matching `HUMAN` `MITOTIC_GT` rows raise.
  - An input row that is not `HUMAN` `MITOTIC_GT` raises (§3 "Input rows"), for instance
    `run_roi`'s own output fed back in.
  - An unknown `rank_key` or `refinement_channel` raises.
  - Input rows with duplicate `object_id`s raise.
  - `evaluate` raises on a seed row whose `object_id` is not a mitotic `ann_id` at that centre
    (§6 step 1). Test it by renumbering a valid seed row's `object_id`.
- **V8 — schema.**
  - Column order and dtypes are as specified.
  - The row key is unique.
  - Prediction `object_id`s are exactly `0…n−1`.
  - The refinement stamp is no later than the prediction stamp.
  - The caller's frame is unmutated.
  - An input frame with int64 coordinate columns comes back as float64 with the same values. Build
    it from the raw JSON of one of the 87 all-int images.
- **V9 — augmented config.** The repo requires every behaviour-preservation check to cover an
  augmented bank.
  - Run `find_and_suppress` with `scales=(0.8, 1.2)`, `n_angles=2` and `flips=(False, True)` on
    1024 px crops of 5 ROIs, before and after the edit. The outputs must be identical.
  - Be clear about what this can show. `find_and_suppress.py` is frozen, so if V2's byte-identity
    holds, V9 cannot fail, and production never reaches an augmented bank anyway
    (`production.py:49` hard-codes one template).
  - It is kept as the repo's standing rule and as a cheap tripwire in case frozen code is edited
    anyway. It is not independent evidence.
- **V10 — evaluation.** On all 441 runs, calling `evaluate(frame, gt, mpp, roi_shape)` with
  `seed_object_id` left at its default:
  - every key `ev.evaluate_run` returns must have the same value in `evaluate`'s metrics as in
    `ev.evaluate_run` called directly on B1's detections, which are already centre-based, with the
    same `gt_eval`, radius, area and shape;
  - `evaluate`'s `tp_at_{10,20,30}` must equal R1's.

---

## 8. Deliberately deferred

- **Multi-loop sessions.** The schema does not block them, but `run_roi` does not accept a
  previous run's output as input (§3 "Input rows"). A second loop passes only the session's
  `HUMAN` rows, and `seed_object_id`, since inference needs exactly one `HUMAN` `MITOTIC_GT` row.
  The caller concatenates each run's `ALGORITHM` rows into the session frame, adding a run
  discriminator column at concat time (§2): predictions from two runs share `object_id`s.
  No accept/reject verdict column is defined yet. When one is added, a verdict must name its
  prediction by run and by `(class == "PREDICTION", object_id)`; `object_id` alone is ambiguous
  across id spaces (§2).
- **A look-alike or rejected class** in the frame. Today ground truth for evaluation stays in the
  `load_annotations` format (§6).
- **Finer refusal reasons.** The frozen `tighten_box_otsu` returns `None` without saying whether the
  click missed the foreground or its component failed `min_area`, `max_area_frac` or
  `min_solidity`.
- **Clicks within 36 px of an edge.** Today they raise (§3 step 2, §4), matching `border_filter` in
  every existing draw. Two cases differ here:
  - a click 25–35 px from the edge whose gate accepts it and whose refined template moves inward
    could be searched without touching frozen code. `run_roi` would lower step 2's margin to the
    `W // 2` px the Otsu window needs (25 for clicks; see §3 step 2's examples). Its fallback would
    still be unreadable.
  - a click closer than 25 px, or any fallback template, would need a padded read, which is a
    change to frozen code.

  **Decided 2026-09-23 (user): neither is supported now.** Every seed draw used so far applies
  the same 36 px margin (`build_seed`, the joint gate, the neuro follow-up's draw), so no
  experiment has used such a click. Revisit only when the product needs edge clicks.
- **Avoiding the second hematoxylin conversion** when `refinement_channel="hematoxylin_od"` (§3).
  That needs production to accept a precomputed channel.
- Using the human's drawn `width`/`height` to size the template: **partly done 2026-09-24** — the
  box sets the gate window and the fallback (§4, D12); an accepted template still comes from the
  Otsu component. Still deferred: rectangular gate windows (the gate uses the enclosing square),
  windows above 73 px (a padded read or a larger `PATCH_SIZE`) or below 19 px, and a radius-based
  refinement around the box centre. Windows 53–73 and templates that large are unvalidated
  (D12).
- Emitting the tight Otsu component bbox as its own row.

---

## 9. Amendments (2026-09-23)

The review behind these amendments found three structural problems:

- the design never said the existing functions stay as they are, yet every recent experiment calls
  them directly;
- it could not run the `default_51` arm;
- its validation checked `tp_at_k` only, which cannot see most of the bugs this document warns
  about.

Each amendment is tagged. **Correction**: the original was wrong or incomplete, and there was only
one right fix. **Decision**: a design choice where a reasonable alternative exists. Decisions are
the ones to review.

| § | amendment | tag |
| --- | --- | --- |
| 0 | Existing functions are frozen; the frame is a new module on top. Rationale: the recent call sites. | decision |
| 0 | `run_production_pipeline` gains `always_od=False`, so that `od` comes from production's own code under either rank key. Rejected alternative: run `chromatin_od` once and re-sort for `tm_score`, which reconstructs tie order by argument rather than by measurement. | decision |
| 0, 7 | `draw_seed_row` becomes `draw_seed_rows`, a port of the driver's joint gate, and it must reproduce the driver exactly. Alternative: have validation call the driver directly and keep the draw out of the package. | decision, confirmed by user 2026-09-23 |
| 1 | Pixel-centre origin stated; clients with edge coordinates subtract 0.5. | correction |
| 1 | "Rounding happens exactly once, through a shared helper" replaced. Rounding stays in the frozen functions, round-half-to-even is pinned, and the one new rounding site is named. | correction |
| 1 | Odd `base_size` is a precondition of `width = base_size − 1`. | correction |
| 1 | Facts added: how the human row is derived from `load_annotations`; mixed int/float JSON and the float64 cast; out-of-image centres; 207 as well as 134; sub-pixel exactness caveat. | correction |
| 2 | Entity identity is class-scoped, and `class` joins the row key. Alternative: a disjoint numeric id range for predictions. | decision |
| 2 | Prediction order is production's returned order, never re-sorted; the tie-break is documented. | correction |
| 2 | Timestamp format pinned. "All three stage values differ" dropped, since it is not guaranteed. | correction |
| 2 | The reason for `round_trip` corrected; the float32 `score` comparison rule added. | correction |
| 3 | Full algorithm; input validation; seed selection requires `class == "MITOTIC_GT"`; row count `len(input) + 1 + n`. | correction |
| 3 | `build_seed` no longer said to take a box. It is unused by `run_roi`, since `run_roi` never draws. | correction |
| 4 | `refinement_channel=None` added (the `default_51` arm), still emitting a refined row with `refined=False`. | decision, confirmed by user 2026-09-23 |
| 4 | "template unreadable" added as a fallback case beside "gate refused"; the refusal strings enumerated. | correction |
| 4 | A click within 36 px of an edge raises, matching `border_filter`. Alternatives (§8): accept clicks 25–35 px from the edge whose refined template is readable, with no change to frozen code; or add a padded read for anything closer and for the fallback. | decision, confirmed by user 2026-09-23 |
| 5 | `info` carries production's whole dict plus named additions; `refinement_requested`, `click_center`, `refusal_reason` and the extra timings defined. | correction |
| 6 | `evaluate` inputs, steps and outputs specified: seed exclusion by the seed row only, candidates sorted by `object_id`, and the area formula. "Match against the human seed's geometry" replaced by what the click actually does. | correction |
| 7 | Validation rebuilt: `draw_seed_rows` defined; references R1–R3; baselines B1–B2; checks V1–V10. | correction |
| 8 | Verdict encoding updated for class-scoped ids; four new deferrals. | correction |

**Second and third rounds, from adversarial verification:**

| § | amendment | tag |
| --- | --- | --- |
| 3, 7 | `run_roi` passes `always_od=(rank_key != "chromatin_od")`. As first written it passed `True`, so nothing exercised production's default `chromatin_od` path after the edit, and an edit that broke that path would have passed every check. | correction |
| 3 | Input rows other than `HUMAN` `MITOTIC_GT` raise. Alternative: copy them through and define row-key handling for repeated runs. | decision |
| 3, 8 | The border step's rationale corrected: it matches `border_filter` in every existing draw. It is not true that no template could be read there. | correction |
| 5 | `n_blanked_px` described as a size check only; position errors are caught by `tpl_sha1` and prediction equality. | correction |
| 6 | `evaluate` excludes by the seed row's own id. Filtering on a `None` argument excludes nothing. | correction |
| 6 | Area uses `mpp * mpp`, bit-identical to the historical formula. | correction |
| 7 | NaN-equals-NaN and `fillna("")` comparison rules. Without them, Step 0 would halt on the first ROI. | correction |
| 7 | V5 exercises `template_unreadable` with a named click, and requires the exact refusal strings. | correction |
| 7 | R4 (the fill-geometry post-D11 run) added; "no post-D11 `tm_score` reference exists" was false. V1 extended to seeds 3–5 and to the contested-tier path. | correction |
| 7 | V3 checks the refined row on all 441 runs; V6 checks both axes; V9's evidential weight stated honestly. | correction |
| 1 | Sub-pixel caveat, the 105 count, the dtype sentence and a line citation corrected. | correction |
| 3 | Input `object_id`s must be unique; duplicates raise. | correction |
| 6 | `evaluate` asserts that the seed's `object_id` is a mitotic `ann_id` at the seed's centre. It is MIDOG-only by construction. | correction |
| 7.1 | `draw_seed_rows` asserts that its `image_id` argument matches `gt`, because it seeds the RNG. The refusal string formats the id with `int()`, as the driver does. | correction |
| 7.2 | V2 adds `git status` for the untracked module; V3 names the sources of R1's columns; V4 notes that B2 has no `od`; V6 pins `rank_key` and `always_od`; V7 adds the input-row, duplicate-id and `evaluate`-seed checks; V10 calls `evaluate` with its default argument. | correction |
| 7.2 | R4's definition fixed, since the clicks CSV has no `condition` column. The shape-vs-mass run is noted as an identical copy of R4's `tm_score` arm. B2-vs-R4 added to the pre-edit environment gate. | correction |
| 8 | Multi-loop text reconciled with the input-row rule: pass only `HUMAN` rows, and add a run discriminator when concatenating. | correction |

**Additions after implementation (2026-09-24):**

| § | addition | tag |
| --- | --- | --- |
| 0, 1, 3 | `click_frame` builds the input row for a raw click; `edge_origin=True` does §1's 0.5 subtraction. Alternative: leave row construction, and the 0.5 subtraction, to each client. | decision, confirmed by user 2026-09-24 |
| 0, 3 | `production_pipeline/run_frame.py` added as the frame path's command line; `run_pipeline.py` stays frozen. | decision, confirmed by user 2026-09-24 |
| — | Status updated: implemented and validated. | correction |
| 0, 3 | `run_frame.py` becomes the single entry point with a required `--source {click,box,midog}`; `--click-id` renamed `--object-id`; `run_pipeline.py` deprecated by a docstring note, code unchanged, deleted once the collaborator repo has the new CLI. Alternative: keep two entry points. | decision, confirmed by user 2026-09-24 |
| 1, 3 | `box_frame` builds the input row for a drawn box, snapped to whole pixels (odd size, centre = pixel containing the drawn centre); the snap's `floor` is a stated carve-out from §1's round-half-up ban. The raw drawn box is not kept, reversing §3's earlier promise that the human's box is preserved. Alternative: store the raw box and snap only when running. | decision, confirmed by user 2026-09-24 |
| 3, 4 | The seed row's box sets the Otsu window and the fallback template (`gate_window_side`: longer side, odd, 19–73 px, raising outside); every input row is checked. 51 for clicks and MIDOG rows, so their runs are unchanged. Alternative: keep the fixed 51 px window and use the box only as a click. | decision, confirmed by user 2026-09-24 (`DECISIONS.md` D12) |
| 3 | `run_frame.py`: a pinned MIDOG annotation that the gate refuses falls back (no raise), and `--source midog` without `--ann-id` uses the joint-gate draw. | decision, confirmed by user 2026-09-24 |
| 4 | `_template_choice` raises `RuntimeError` if the gate window is unreadable, so a refusal always means the gate refused. | correction |
