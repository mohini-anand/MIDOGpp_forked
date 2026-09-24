# Plan: one pipeline entry point, and a drawn box as a seed source

Status: **implemented and verified 2026-09-24; not yet committed (awaiting the user's go).**
Results: `Research Logs/2026-09-24-single-entry-point-and-box-seed.md` (V-A 19/19 identical; full
tier 0 divergences, 49/49 pickles identical; `box_checks.py` 0 divergences; Step 0 CLI frames
identical). Written 2026-09-24. Reviewed adversarially in two rounds (round 1:
2 MAJOR + 17 MINOR/NIT, all accepted; round 2: AGREE on v2, 4 MINOR + 5 NIT, applied as written in v3).
v4 records the user's decisions (§1); it was not re-reviewed. Changelogs at the end.

Worktree: `/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked_annotation_frame`, branch
`feature/annotation-frame` at `d224c80`. Python: `/Users/mohinianand/anaconda3/bin/python3`.

## 0. What the user asked for, and what is already decided

1. **One entry point** that takes a flag for the seed source and builds the seed row accordingly.
   Today: `production_pipeline/run_pipeline.py` (legacy: `ss.build_seed`, gray-only gate with redraw, the
   only path with a plot) and `production_pipeline/run_frame.py` (annotation frame: `click_frame` /
   `human_seed_frame` / `draw_seed_rows` -> `af.run_roi` -> `af.evaluate`; no plot). On 300.tiff ann 14499
   both give the same template and bit-identical detections and buckets (measured twice, independently).
2. **A new seed source: a box the human draws.** Decided with the user:
   - the box's centre is the proxy click (the gate's centre-pixel rule is unchanged; a radius-based
     refinement comes later);
   - the box is the Otsu gate's **window**: `tightened_template_box(..., otsu_window=<window side>)`;
   - the window is **square**, from the box's **longer** side;
   - a box larger than **73 px** (`tm.PATCH_SIZE`) **raises**;
   - on gate refusal (or an unreadable refined template) the fallback template is the human's box — for a
     non-square box, the enclosing square on its longer side, centred on the box centre.
3. The rule "Otsu window and fallback template = the human row's box" applies to **every** human row.
   Every click row (`CLICK_BOX = 50.0`) and every MIDOG row (all 26,286 raw boxes are 50x50;
   `load_annotations` gives w = h = 50.0) has frame width 50 -> window 51 and fallback 51, exactly today's
   `OTSU_WINDOW` / `DEFAULT_BASE_SIZE`. Every validation-harness input row also has width 50
   (`synthetic_seed`, V8's raw JSON, `human_seed_frame`). Click and MIDOG runs must stay bit-identical.

Probe backing for (2): `box_window_probe.py` (to be added at the repo root with its output, §3): 993
mitotic annotations, 14 original ROIs, gray channel, square windows of (51-px component size + 2*margin)
centred on that component: acceptance 95.5-96.7% at margins 0-15 px, refusals almost all "centre pixel not
foreground"; template size tracks the margin (margin 0: 45% of templates smaller than at 51, none
larger; margin 15: 59% larger). The
reviewer's probes (not in the repo; 300.tiff, 40 clicks x 2 channels): acceptance 28/80 at window 9, 44/80
at 11, 66/80 at 51; a 9 px template at 300.tiff/14499 re-detected the seed's own figure (rank 62, 9.5 px
away) past its 9x9 blank.

## 1. Decisions (user, 2026-09-24)

- **D-a. `run_pipeline.py`: deprecate now, delete later — confirmed.** A docstring-only deprecation note
  pointing at `run_frame.py`, at the docstring's 4-space indent. Code-identical, proven with
  `validation/annotation_frame/ast_equal.py`; only its `--help` description changes, because argparse builds
  it from `textwrap.dedent(__doc__)` (`run_pipeline.py:170`). Spec §0's freeze covers signatures, defaults,
  return types, columns and numerics, not help text. Delete once the collaborator repo
  (`~/Desktop/RA Work/AnnotateDx`, which ships `run_pipeline.py` at its root) has the new CLI. Kept today
  because: the user froze it on 2026-09-24 (spec §0/§9); it is the only CLI for `build_seed`'s gray-only
  draw (300.tiff seed 0: 14581, vs the joint draw's 14499); `SELF_HIT_MASKING_PLAN.md:673,714` and
  `bbox_refinement_three_way_chromatin_od_audit_round3.py:155-157,439` call it.
- **D-b. MIDOG with no `--ann-id` uses the joint-gate random draw `draw_seed_rows` — confirmed** (the draw
  R1/R2 were made with).
- **D-c. A pinned MIDOG annotation the gate refuses falls back, no raise — confirmed** (current `run_roi`).
- **D-d. Window side >= 19 px, i.e. drawn boxes of 18 px or more — confirmed.** 19 is the smallest
  `base_size` in the 441 validated reference runs (spec §5). Below it the gate mostly refuses and the
  fallback template is tiny, and a small template means a small D11 blank, so the seed's own figure can leak
  back (the 9 px example above). The mechanical minimum would be 9. The floor bounds the window and the
  fallback, not the accepted template: the gate can still return templates below 19 (hematoxylin accepted
  17 px templates at windows 19 and 21 on 300.tiff), as it already can at 51. Windows 53-73 are equally
  unvalidated; D12 says so.
- **D-e. Box geometry snapped to whole pixels — confirmed** (after a clarifying question). Size `n = round(w)` raised to odd, centre = the pixel containing the drawn centre (§2.1). The
  stored human row is then exactly the gate window and the fallback footprint. Nothing is cut off: an
  even-pixel box grows by one pixel, and the drawn centre moves by at most 0.5 px to a pixel centre. What
  is lost is information: the frame cannot tell a drawn 48 px box from a drawn 49 px one, nor recover
  sub-pixel coordinates. If the planned radius refinement needs the raw box, the client must keep it.
  This reverses spec §3's promise that the human's own box is preserved.
- **D-f. `--source {click,box,midog}` is required — confirmed.** Plain `run_frame.py 300.tiff` no longer
  runs the draw; `--click-id` becomes `--object-id`.
- **D-g. The box-size check runs on every HUMAN input row — confirmed** (after first choosing seed-only,
  the user switched once it was clear the check concerns other human rows in the same input, e.g. a
  multi-click session, not predictions). A bad non-seed row fails at input, like other row errors.

## 2. Code changes

Two files change: `midog_utils/annotation_frame.py` (new module, not frozen) and
`production_pipeline/run_frame.py`; plus the docstring-only note in `run_pipeline.py` (D-a). **No frozen
file changes.** `draw_seed_rows` / `_joint_spec` keep their hardcoded 51 (a port of the driver).

### 2.1 `midog_utils/annotation_frame.py`

Constants:

    MAX_WINDOW_PX = tm.PATCH_SIZE   # 73
    MIN_WINDOW_PX = 19              # D-d

Helper (public; all arguments and returns in the units named):

    def gate_window_side(width, height) -> int:
        # width, height: frame extents (x2 - x1), so a side's pixel count is extent + 1 (spec §1's +1).
        # ValueError unless both are finite and >= 0.
        n_px = int(round(max(float(width), float(height)))) + 1      # Python int; round is half-to-even
        side = n_px if n_px % 2 else n_px + 1
        # ValueError if side > MAX_WINDOW_PX or side < MIN_WINDOW_PX (message: box pixel count, window, limit)
        assert side % 2 == 1 and side <= tm.PATCH_SIZE
        return side

  Vectors, frame-width units: 50.0 -> 51; 51.0 -> 53 (a no-+1 version gives 51); 71.0 -> 73 (72 px);
  72.0 -> 73 (73 px); 73.0 raises (74 px -> 75); 17.0 -> 19 (18 px); 18.0 -> 19; 16.0 raises (17 px);
  non-square uses the longer; NaN / negative raise.
  - Odd is required: `footprint_box` raises on even sizes; `build_augmentations`, `blank_seed_square` and
    `read_padded_patch` use `2 * (n // 2) + 1` pixels, so an even value would silently act as n + 1.
  - side <= 73 guarantees the refined `base_size` (<= window) fits the 73 px patch `build_augmentations`
    crops from, and the window is readable after run_roi's 36 px border check (equality: 73 // 2 = 36).

Builder:

    def box_frame(file_name, x0, y0, w, h, object_id, *, current_time=None) -> pd.DataFrame:
        # x0, y0, w, h: left, top, width, height in pixel-edge coordinates (canvas / COCO xywh; pixel k spans
        # [k, k + 1)). Per axis:
        n = int(round(w)); ValueError if n < 1
        n_odd = n if n % 2 else n + 1
        c = math.floor(x0 + w / 2)            # the pixel containing the drawn centre; no rounding tie
        x_tl, width = float(c - n_odd // 2), float(n_odd - 1)
        # same for y; then gate_window_side(width, height) validates the size (raises)
        return _block(1, x_tl, y_tl, width, height, class=MITOTIC_GT, object_id, made_by=HUMAN, refined=False, ...)

  - Size depends only on w; position only on the drawn centre. The stored centre `x_tl + width / 2 = c`
    is an integer pixel, so the frozen functions' half-to-even rounding never meets a tie on box rows.
  - The stored row IS the gate window (square case) and the fallback footprint; on fallback the refined row
    equals the human row's geometry.
  - Not preserved (D-e): sub-pixel drawn coordinates, and the difference between an even size and the next
    odd one. Documented.
  - `box_frame(x0, y0, 51, 51) == click_frame(x0 + 25, y0 + 25)` (all columns but `current_time`) holds
    for **integer** origins only. For other origins the box snaps to the pixel containing its centre, while
    a click keeps its sub-pixel centre, which the frozen functions round half-to-even — usually the same
    pixel, but not when `x0 + w/2` is an integer (e.g. W = 51 at x0 = 11.5: box -> 37, click 36.5 -> 36).
    No equivalence is claimed there.
  - `floor(x0 + w/2)` is, in the frame's centre coordinates, `floor(c + 0.5)` — the round-half-up form spec
    §1 bans. The carve-out is deliberate and is written into §1 and the docstring: the ban protects
    centre-to-pixel conversions that must agree with `read_padded_patch`; this one picks the pixel
    containing an edge-coordinate point, has no tie, is position-independent, and yields integer centres
    so the frozen rounding never meets a tie. Changing it to `round()` would bring back the origin
    dependence of v1.
  - Not for MIDOG boxes (already pixel-centre; they go through `human_seed_frame`). Docstring says so.

`_template_choice(rgb, cx, cy, refinement_channel, window)` — `DEFAULT_BASE_SIZE`/`OTSU_WINDOW` become
`window`:

    None channel        -> (window, cx, cy, False, None)
    (new guard)         -> RuntimeError unless ss._patch_readable(rgb.shape, cx, cy, window)
    gate refuses        -> (window, cx, cy, False, "gate_refused")
    refined unreadable  -> (window, cx, cy, False, "template_unreadable")
    accepted            -> (b, tx, ty, True, None) with tightened_template_box(..., otsu_window=window)

  The guard keeps "None means refused" true by construction, not by 36 == 36.

`_input_rows`: also check `gate_window_side(width, height)` on **every** input row (D-g). Placed **after**
the column check and the HUMAN MITOTIC_GT check, so a missing column still raises the columns ValueError
and a fed-back output (whose refined row can be 9-17 px) still raises "not HUMAN MITOTIC_GT". No harness
outcome changes: every V7 input row is width 50, and `expect_raise` checks only that a call raised.
`run_roi` still computes the seed's window itself (below), so the value used is never a stale one.

`run_roi`: `window = gate_window_side(seed["width"], seed["height"])` right after `_select_seed` and the
file_name check (before any image work); step 2 unchanged; pass `window` to `_template_choice`; the edge
error message says "seed centre", not "click". **No new info keys** (`checks.py:386` treats unknown keys as
production keys and `compare_info` would flag them on all 441 runs); the window is recoverable from the
frame's human row via `gate_window_side`, and the fallback size is `info["base_size"]`.

Docstrings/comments: module docstring; `click_frame` (CLICK_BOX now also fixes the 51 px window and
fallback); `DEFAULT_BASE_SIZE`'s comment ("and the fallback" -> clicks/MIDOG only); `footprint_box` ("the
one rounding site" -> no longer the only one: `box_frame`, `gate_window_side`); `_template_choice`;
`run_roi` (new raises; `refinement_channel=None` means "the human's box, 51 px for clicks").

### 2.2 `production_pipeline/run_frame.py` — the single entry point

CLI:

    run_frame.py ROI --source {click,box,midog} [source options] [common options]
      click:  --click X Y [--object-id N] [--edge-origin]     X Y in pixel-centre coords unless --edge-origin
      box:    --box X Y W H [--object-id N]                   always pixel-edge coords (canvas / COCO xywh)
      midog:  [--ann-id N | --seed-index K]                   no --ann-id: joint-gate random draw, seed K (D-b)
      common: --rank-key, --refinement-channel {gray_inverted,hematoxylin_od,none}, --images-dir, --db-path,
              --show N, --save-csv PATH, --top-n N, --save-fig PATH, --plot

  - `--source` required (D-f). Source-specific options default to `None` in argparse (and in
    `build_seed_frame`), so "given" is detectable; `parser.error` on another source's option or a missing
    required one; 0 / False applied after validation.
  - A `--box` is validated with `af.box_frame` before any loading; its ValueError becomes `parser.error`.
  - `--click-id` renamed `--object-id` (click and box). Added in 41ffa2e; only docs mention it.
  - `--refinement-channel` help: "none = no tightening: the template is the human's box (51 px for clicks
    and MIDOG)".

Functions:

  - `build_seed_frame(fn, rgb, gt, image_id, source, *, click=None, box=None, object_id=None, edge_origin=None, ann_id=None, seed_index=None) -> (frame, scorable)`:
    ValueError on an unknown source or another source's argument (mirrors the CLI).
      click -> af.click_frame(fn, x, y, object_id or 0, edge_origin=bool(edge_origin)), False
      box   -> af.box_frame(fn, x, y, w, h, object_id or 0), False
      midog -> ann_id: af.human_seed_frame(the mitotic row), True; else draw_seed_rows(range(k+1))[-1][0], True
  - `run_frame_on_roi(fn, source, *, ..., keep_rgb=False) -> dict`: as today plus `source`; `rgb` only when
    `keep_rgb` (93 MB per result on 300.tiff).
  - `visualize(result, ax=None, downsample=4, top_n=10) -> Axes` (lazy matplotlib / viz imports):
      - seed row: the HUMAN MITOTIC_GT row with `object_id == info["seed_object_id"]`; refined row: the
        ALGORITHM MITOTIC_GT row with that id.
      - detections: the scored rows when scored, else PREDICTION rows; stable-sorted by object_id;
        `cx, cy = af.box_center(...)`; `bucket` only when scored.
      - ground truth: `gt_out` when scored; an empty frame with `cx, cy, category_id` otherwise.
      - `viz.overlay(result["rgb"], gt, dets, seed_xy=info["click_center"], tpl_xy=info["template_center"], ax=ax, downsample=downsample, top_n=top_n, radius=info["match_radius_px"], title=...)`
      - `viz.draw_box(..., scale=1.0 / downsample)` for (a) the gate window in white:
        `iy, ix = int(round(cy)), int(round(cx)); h = window // 2; (iy - h, iy + h + 1, ix - h, ix + h + 1)`
        with `window = af.gate_window_side(seed width, seed height)` — exact for every row type, including
        MIDOG rows with half-integer x_tl; and (b) the refined footprint in cyan from the refined row:
        `(y_tl, y_tl + height + 1, x_tl, x_tl + width + 1)`.
  - `summarize` adds `source`. `main`: `--save-fig` saves, `--plot` shows, neither -> no plot; sets
    `keep_rgb` when either is given. Update the example-invocation comment block (lines 14-20).

### 2.3 `production_pipeline/run_pipeline.py` (D-a)

Docstring-only deprecation note (module docstring) pointing at `run_frame.py --source midog`. Verify with
`ast_equal.py`. Deletion is a later, separate step after the collaborator port.

## 3. Documentation

- `ANNOTATION_DATAFRAME_SPEC.md`: status header (65 checks; "full §7 harness was not re-run" -> the
  re-run); §1 ("the one rounding site" -> add `box_frame`'s snap and `gate_window_side`, plus the explicit
  floor carve-out from the round-half-up ban at lines 110-113); §2 table (human row may be a snapped drawn
  box); §3 builder table (+ box), "What the pipeline consumes from the human row" (now centre and size),
  the sentence at lines 324-327 promising the human's own box is preserved (now: snapped, D-e), and CLI
  bullets (incl. line 359's `--click-id`); §4 table and bullets (fallback = human box, 51 for
  clicks/MIDOG; "the 51 px Otsu window (which needs a 25 px margin)"; the new raises and guard); §8
  (edge-click deferral: the margin would be `window // 2`; rectangular windows, >73 px and <19 px windows
  deferred; the "use the human's width/height" item now partly done); §9 new Additions rows (decision,
  confirmed by user 2026-09-24).
- `DECISIONS.md`: D12 "A human's box sets the Otsu window and the fallback template": rule; 51 px
  equivalence for clicks/MIDOG; 19..73 window limits (bounding the window and fallback; accepted templates
  can still be smaller, the gate's own floor being 9); probe numbers and caveats (simulated, square,
  centred, gray only, precision not measured); windows < 19 and 53-73 unvalidated; templates 53-73 are not
  rotation-safe inside the 73 px patch (no effect at production's single unrotated template).
- `D8_TEMPLATE_ANCHOR.md:14`: one-line pointer to D12.
- `Research Logs/2026-09-24-single-entry-point-and-box-seed.md`: probes, validation re-runs, the D-a
  kept-for-now list, the collaborator-port follow-up.
- `box_window_probe.py` at the repo root, path from `Path(__file__).resolve().parent`, writing its result
  table to `results/box_window_probe.csv`; both committed. `out/box/box.json` committed.
- This plan: status line updated to "implemented" with the commit hashes.

## 4. Verification (nothing committed without the user's go)

Total machine time is roughly 1.5-2 hours; most of it is the full tier.

Step 0, **before any edit**: work continues on `feature/annotation-frame` (user, 2026-09-24: no new
branch; `d224c80` is the restore point and the baseline every comparison below names); save
`git show d224c80:midog_utils/annotation_frame.py` to a scratch directory outside the repo; save HEAD's
CLI outputs with `--save-csv` for `--ann-id 14499`, the draw (seed 0), `--click 1411 3551`, and
`--ann-id 14499 --refinement-channel hematoxylin_od --rank-key chromatin_od` (300.tiff), and
`run_pipeline.py 300.tiff --ann-id 14499`'s printed output, to the same place.

V-A. **Old-vs-new differential.** Import the saved HEAD module as `midog_utils._af_head` (importlib
     `spec_from_file_location`; shares `ss`/`tm` with the live package). ROIs loaded via
     `rr.build_rois()`'s `path` column (462/072/007 are under `testing_set/`). Per case, run both `run_roi`s;
     compare frames on every column but `current_time`, and info with `C.same_value` on every non-`t_` key.
       - 300.tiff ann 14499: 3 refinement channels x 2 rank keys;
       - V5 refusals: 462.tiff/22673 gray, 529.tiff/24934 gray, 300.tiff/14581 hem, 245.tiff/6274 hem;
       - V5 unreadable: 300.tiff (6410.0, 250.0);
       - V6 half-integer clicks: 072.tiff/1835, 007.tiff/144;
       - a `click_frame` row, 300.tiff (1411, 3551), and a two-row input with `seed_object_id`.
     0 divergences required.
V-B. **Existing harness, unchanged code.** smoke tier; `named_checks.py --tier full`; `extra_checks.py`;
     full tier as **3 concurrent shards** (~47 min; ~2.3 h if serial). `run_tier.py` skips existing
     pickles, and `out/` is tracked: delete the tier's `rois/*.pkl` (git holds them), run, then before
     restoring compare every key of every regenerated pickle except `t_check_s` (tally, r1_rows,
     r2_frames, sel_rows, r4_tops, n_half_even_vs_half_up_differ; `t_*` excluded inside them) against
     `git show HEAD:<path>`; aggregate without `--git` (it cannot pass for any base now); then
     `git restore validation/annotation_frame/out`. Replacement for V2's repository half:
     `git diff --name-only d224c80 -- midog_utils/ production_pipeline/ rerun_bbox3way_postD11.py midog_utils_full/`
     lists only `midog_utils/annotation_frame.py`, `production_pipeline/run_frame.py` and (D-a)
     `production_pipeline/run_pipeline.py`, whose change `ast_equal.py` proves docstring-only.
V-C. **New box checks**, `validation/annotation_frame/box_checks.py` -> `out/box/box.json`. Expected values
     computed by hand from pixel counts, never with `af.gate_window_side`:
       1. `gate_window_side` vectors of §2.1 (frame-width units), incl. 51.0 -> 53; Python-int return type.
       2. `box_frame`: (0,0,51,51) -> x_tl 0, width 50; (0,0,52,52) -> width 52 (window 53);
          (0,0,40,40) -> 41; (0,0,73,73) -> 73; (0,0,74,74) raises; (0,0,17,17) raises; (0,0,18,18) -> 19;
          empty raises; shift behaviour, for widths 18-73 in 0.5 steps and x0 in {k, k+.25, k+.5, k+.75}
          with k even and odd: the side is the same across fractional offsets, and
          `box_frame(x0 + 1).x_tl == box_frame(x0).x_tl + 1`; integer-origin equivalence with
          `click_frame` row-for-row; non-square stores both axes, window = longer.
       3. The gate uses the box window, on 4+ accepted boxes with side != 51 and
          `tightened_template_box(side) != tightened_template_box(51)` (verified first): `info["base_size"]`
          and `template_center` equal an independent `ss.tightened_template_box(channel, cx, cy, otsu_window=side)`;
          `base_size <= side`; `n_blanked_px == base_size ** 2`; refined row == `footprint_box(...)`. Cases
          include an even-pixel box, a non-square box and a half-integer drawn origin.
       4. Fallback: a verified-refused box with side != 51 -> `(side, cx, cy)`, `refined False`,
          `gate_refused`; refined row geometry == human row for a square box, == the enclosing square
          (longer side, same centre) for a non-square one; predictions identical to
          `refinement_channel=None` on the same box. A `template_unreadable` box if a candidate can be
          found (search near V5's unreadable click; record the outcome either way).
       5. Edge: a 73 px box centred exactly 36 px from an edge runs; the guard never fires on valid input.
       6. `run_roi` raises on a hand-built seed row of width 80, and on a two-row input whose non-seed
          row has width 80, width 16 (window 17) or NaN, whichever row is the seed (D-g).
V-D. **CLI.** Each source on 300.tiff (`click 1411 3551`, `box` around ann 14499, `midog --ann-id 14499`,
     `midog --seed-index 0`) with `--save-fig`; PNGs inspected by eye (window, footprint, star, buckets).
     The Step 0 CSVs re-run under the new CLI (`--source midog ...`, `--source click ...`) are identical
     (read with `float_precision="round_trip"`, `current_time` excluded). Wrong-source options, missing
     options, an explicit `--seed-index 0` with `--source click`, and a 74 px `--box` exit with argparse
     errors before any loading. `run_pipeline.py --help` shows the deprecation note; `run_pipeline.py
     300.tiff --ann-id 14499` output unchanged from Step 0.

## 5. Commits (only when the user says so)

1. `annotation_frame`: box window/fallback, `box_frame`, `gate_window_side`, the guard, input-row check.
2. `run_frame.py`: `--source`, box branch, `visualize`; `run_pipeline.py` deprecation docstring.
3. Docs, D12, research log, probe + output, `box_checks.py` + `out/box/box.json`, this plan's status.

Follow-ups outside this plan: port the new CLI to the collaborator repo; then delete `run_pipeline.py`;
measure precision for box-sized templates (windows 19-73) before relying on them.

## Changelog v1 -> v2 (reviewer round 1)

1 vectors restated in frame-width units, discriminating vectors added. 2 position-independent snapping
(`round(w)` -> odd, centre = floor), row = window. 3 V-C gate/fallback cases forced off 51. 4 Step 0
baselines saved before editing. 5 ROI paths via build_rois; `same_value`; Python int. 6 pickle comparison
before restore, delete-not-move, concurrent shards, replacement git check. 7 readability guard. 8 visualize:
`scale=1/downsample`, window drawn from rounded centre, seed row by id. 9 stale docs listed. 10 D-d floor 19.
11 D-a deprecate-then-delete. 12 None defaults. 13 box validated before loading. 14 coordinate conventions
in help. 15 every input row checked. 16 `keep_rgb`. 17 rotation caveat in D12. 18 `gate_window_side` name,
odd/<=73 assert, "enclosing square" wording. 19 probe path and output committed; `box.json` committed.

## Changelog v2 -> v3 (reviewer round 2, verdict AGREE)

1 non-integer-origin equivalence claim corrected (x0 + w/2 integer case). 2 explicit carve-out for
`floor(x0 + w/2)` from spec §1's round-half-up ban, in spec and docstring. 3 spec §3's "human's own box is
preserved" sentence rewritten; raw box not recoverable (D-e consequence). 4 `_input_rows` size check placed
after the column and class checks. 5 19 px floor bounds the window/fallback, not accepted templates.
6 V-C.4 geometry check square-only / enclosing square. 7 V-B pickle comparison covers every key but
`t_check_s`. 8 D-a worded as code-identical with a changed `--help`; `--help` added to V-D. 9 V-C.2 shift
test made non-circular.

## Changelog v3 -> v4 (user decisions, 2026-09-24; not re-reviewed)

D-a, D-b, D-c, D-d, D-f, D-g confirmed as proposed (D-g after a brief detour to seed-only, reverted by the
user; the v3 text for `_input_rows` and V-C.6 stands). D-e confirmed after a clarifying question; its
wording now says what is lost (information about the drawn size and sub-pixel position), not "size".
Step 0 saves to a scratch directory outside the repo. Every decision is now settled.

## Implementation notes (2026-09-24): where the code departs from v4

- `run_frame.check_seed_args` holds the source/argument rules; `main`, `run_frame_on_roi` and
  `build_seed_frame` all call it, so a function caller also gets the errors before anything loads.
  It also rejects `--ann-id` together with `--seed-index`, and a negative `--seed-index`.
- `visualize` draws the white gate window only when a gate ran (`refinement_requested`); with
  `refinement_channel=None` there is no window. With `ax=None` it adds a full-resolution zoom panel
  on the seed (window and footprint are about 13 display px on the 4x-downsampled overview).
- `box_frame` re-raises `gate_window_side`'s error with the drawn and snapped sizes, so a 74 px box
  says "drawn box 74.0 x 74.0 px snaps to 75 x 75 px: ...".
- V-B: the regenerated pickles and summaries were restored from git as planned, but the two
  `named.json` files are kept regenerated: they differ from the committed ones only in timings and
  the four edge-error messages, whose wording changed ("click" -> "seed centre").
- V-C.3's cases were first drawn from one annotation (six of seven on 300.tiff/14459, differing from
  51 px only by a half-pixel centre); `box_checks.py` now takes seven different annotations, each
  verified to change the template **size** against the 51 px window.
- V-D: `run_pipeline.py` without `--save-fig` calls `plt.show()`, which blocks under the macOS
  backend; its unchanged-output check ran with `MPLBACKEND=Agg` (environment only).
