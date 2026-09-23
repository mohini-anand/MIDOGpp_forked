# Task: convert midog_utils's internal point representation to box corners (TL/BR)

## Goal

Convert the internal representation of seeds, templates, and detections in `midog_utils/`
from center points (`cx, cy`) to box corners (`x1, y1, x2, y2` — top-left / bottom-right,
matching `dataset.py`'s existing "TLBR, not xywh" convention), while leaving
`nms_by_distance` and `evaluate.greedy_match` **semantically and numerically untouched** —
they keep operating on centers, derived on the fly at the narrowest possible set of
boundaries. This is a pure internal-representation refactor: no algorithm changes, no
IoU-based NMS, no change to what gets suppressed or matched.

**Why**: a prior conversation considered replacing the custom point-radius NMS with a
standard box+IoU NMS and rejected it — axis-aligned box overlap is not equivalent to the
physically-motivated circular match radius that `evaluate.radius_px(mpp)` /
`invariants.check_nms_radius` (DECISIONS.md D7) enforce. This refactor is **not** a step
toward that; it is purely representational cleanup, motivated by wanting the *detections*
DataFrame (the CSV `find_and_suppress`/`production.py` emit) to carry corner coordinates as
its source of truth instead of a center. Do not let this scope creep into an NMS redesign.

**GT annotations (`dataset.load_annotations`) are explicitly NOT switched to corners-only —
see Step 2.** They stay backward-compatible. The corners-not-centers goal applies to the
*detections* representation; `dataset.py` gets corners added, not centers removed.

## Critical technical landmine — read before writing any code

At least **four** places in this codebase construct a box from (or reduce a box to) a center,
and they do not all use the same convention. Get this wrong and centers silently shift by up
to 0.5px in a way that shows up as small, plausible-looking recall/precision drift, not a
crash.

1. **GT boxes from `MIDOG++.json`** (`dataset.py:58,65`): `cx = (x1+x2)/2.0`, a naive average,
   no correction — **closed/symmetric convention**. Correct because `BOX_SIZE = 50`
   (`dataset.py:12`) is a symmetric "click ± 25" box, not a pixel-counting bbox — every
   annotation in the real 26,286-row `MIDOG++.json` has `w = x2-x1 = h = y2-y1 = 50` exactly
   (cross-checked directly against the file; also asserted at `rerun_bbox3way_postD11.py:294`).
2. **Otsu-derived template boxes** (`seed_selection.py:101,117-119,134-139`): skimage's
   half-open pixel-grid bbox (`region.bbox`, max exclusive) — **half-open convention**.
   `tightened_template_box` already documents this explicitly: `(x0 + x1 - 1) / 2`, not
   `(x0 + x1) / 2`, because a half-open bbox covers pixels `x0..x1-1`. `DECISIONS.md`
   D8 (~L583-585) states the same formula and explains the half-pixel difference in words —
   read it for the reasoning, but verify against the code, not the prose.
3. **The synthetic `click_box` you will construct in Step 3** around a raw annotation's
   `(cx, cy)` (e.g. for the `Seed.click_xy` replacement). Must use the **closed** convention
   (same as #1) to round-trip exactly, since it is built from the same closed-convention `cx`.
4. **The synthetic box you will construct in Step 5** around each surviving NMS center
   (`box_from_center`, built from `extract_peaks`' integer-valued float64 output,
   `template_match.py:202`, plus `cfg.base_size`). Must use the **closed** convention to
   round-trip exactly back to the original integer center — a half-open ("-1") box here would
   reintroduce the exact 0.5px bug this section warns about, at a brand-new site.

Resolve this **before** touching any caller of a geometry helper, by building two distinct,
clearly-named helper functions (`center_of_closed_box` / `box_from_center_closed` and
`center_of_half_open_box`), used only where each convention actually applies — do not try to
normalize everything to one convention, since #2's Otsu box is a real pixel-grid region with
a shape you don't control (see Step 3's note on `tightened_template_box`'s return value).

Document, in your final report, which sites use which helper and why. This must be backed by
the whole-dataset verification check in Verification item 1a below — a plausible-sounding
argument is not what makes this landmine defused, a passing check across all 26,286
annotations is.

## Non-goals / scope boundary

- Do **not** change `nms_by_distance`'s algorithm, `greedy_match`'s algorithm, or any
  evaluation semantics (match radius, bucket definitions, recall/precision formulas).
- Do **not** introduce IoU-based NMS or any box-overlap suppression logic anywhere.
- Do **not** touch the `find-and-suppress-midog` branch or anything in the original working
  directory (`/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked`) once your worktree
  exists — this includes all untracked research files (CSVs, notebooks, `Research Logs/`,
  `*_PROMPT.md` files, `production_hematoxylin_only/`, `human_neuroendocrine_followup/`,
  `debug_pics/`, everything currently untracked per `git status`).
- Do **not** attempt to migrate the ~130 downstream research scripts/notebooks that consume
  `cx`/`cy` columns from these APIs. This refactor proves functionality preservation for the
  core `midog_utils/` pipeline only, in isolation. Downstream migration is a separate future
  decision, not part of this task. **The one named exception** is
  `rerun_bbox3way_postD11.py` — see "Regression driver" below; adapting it is in-scope
  verification infrastructure, not downstream migration, and it is the *only* file outside
  the 10-file `midog_utils/` set (plus the new `geometry.py`) you may edit.
- Do **not** rewrite `DECISIONS.md` D7-D11 — their substance (radii, semantics) doesn't
  change. If you think a documentation note is warranted, propose it in your report; don't
  edit the file.
- Do **not** commit anything unless asked, and do not delete the worktree unless asked.

## Setup

1. From the current repo, on branch `find-and-suppress-midog`, create a git worktree in a
   sibling directory on a new branch:
   `git worktree add ../MIDOGpp_forked_corner_refactor -b refactor/corner-coordinates`
2. Confirm immediately afterward, in the **original** directory: `git status --short` is
   unchanged from before you ran that command (same untracked files, still on
   `find-and-suppress-midog`). If anything looks different, stop and report — do not proceed.
3. **`images/` is gitignored** (`.gitignore:137-138`, ~6+ GB of ROI TIFFs fetched separately
   via `Setup.ipynb`) and will **not** be populated by `git worktree add` — only tracked files
   are. Every ROI read (`dataset.load_roi`) will fail with `FileNotFoundError` until you fix
   this. Symlink it from the original directory before doing anything else:
   `ln -s /Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked/images ./images` (run from
   inside the new worktree directory). `databases/MIDOG++.json` is tracked, so it needs no
   such fix — confirm this assumption still holds (`git ls-files databases/` from inside the
   worktree) rather than trusting this document.
4. All work happens inside the new worktree directory from here on.
5. Use `/Users/mohinianand/anaconda3/bin/python3` for anything touching `cv2`/`skimage` — the
   bare `python3` on this machine cannot import `cv2` (confirmed: NumPy ABI mismatch).

## Regression driver — read this before "Read first" below

There are two candidate "regression reference" scripts in this repo. **Only one is real.**

- `chromatin_od_ranker_seed_robustness_audit.py` is the literal source of the "84/84,
  0 divergences" phrase referenced elsewhere. **It is already broken against the current,
  stock `midog_utils/`, for reasons unrelated to this task**: its `check_config()` reads
  `FSConfig().channel` and `FSConfig().self_hit_radius`, neither of which exists on the
  current `FSConfig` (`find_and_suppress.py:25-37`) — `self_hit_radius` was D10's mechanism,
  removed by D11. Confirm this yourself (`AttributeError` on both) and then **do not use it
  and do not fix it** — it's a pre-existing, unrelated staleness, out of scope, not your bug.
- `rerun_bbox3way_postD11.py` is the actually-live driver — it calls
  `production.run_production_pipeline(rgb, seed, mpp, rank_key=...)` matching today's real
  signature. **Use this as your Step 0 baseline driver.**
- `rerun_bbox3way_postD11.py` will itself break under your Steps 2/3 in more than one or two
  spots — don't assume "minimal" means "one line." Confirmed breakage sites: `:91`/`:98`
  (unpacks `tightened_template_box(...)` as a 3-tuple `(base_size, cx, cy)` — Step 3 changes
  this shape), `:241` (`ss.Seed(..., click_xy=click_xy, template_xy=(tx, ty), ...)` —
  Step 3 renames these fields to `click_box`/`template_box`), and `:294` (asserts
  `annotations['w'] == BOX_SIZE`, which Step 2's backward-compatible column addition should
  keep working — verify this, don't assume it). The same `(base_size, tx, ty)` spec shape
  and the `Seed` fields also thread through `template_spec`, `template_digest`,
  `pool_gate_specs`, `draw_joint_click`, `selection_rows`, and `run_roi` — expect to touch
  most of this file's spec-threading plumbing. **You are authorized to make mechanical edits
  to this one script** — adapting its read sites, tuple unpacking, and field/column names to
  your new APIs — so it keeps working as your verification driver. Do not change its ROI
  selection or scoring *logic* (what gets selected, what gets scored, how) — only how it
  reads coordinates out of the objects your refactor changed the shape of.

## Read first, fully, before editing anything

- All of `midog_utils/`: `nms.py`, `evaluate.py`, `dataset.py`, `seed_selection.py`,
  `template_match.py`, `chromatin.py`, `viz.py`, `find_and_suppress.py`, `production.py`,
  `invariants.py`. ~1300 LOC total, all in scope. Read every one in full before editing any —
  several functions below are only safe to change once you've seen every caller.
- `DECISIONS.md` D7 (NMS radius = eval radius, mpp-scaled), D8 / `D8_TEMPLATE_ANCHOR.md`
  (template anchoring, the half-open-bbox formula), D9 (`max_peaks`), D10/D11 (self-hit
  handling via `blank_seed_square`) — so you understand why the current center-based design
  looks the way it does before changing it.
- `rerun_bbox3way_postD11.py` in full — it's both your regression driver and a file you're
  authorized to edit; know what it actually does before touching it.

## Execution plan (Step 2 runs first in isolation; Steps 3-9 are grouped by file, not a strict
dependency chain — nothing in 3-9 calls `dataset.py`'s functions for geometry, so their
internal order doesn't matter, but every one of them must land before Verification runs)

**Step 0 — capture the baseline.** Before any edit, run `rerun_bbox3way_postD11.py` (see
"Regression driver" above) on its existing fixed seed-set/ROI-set against **stock** code, and
save: full detection DataFrames (not just summary metrics), the full GT annotations frame
(`load_annotations`'s output), and the **full per-run `info` dict** `find_and_suppress`
returns (not just the driver's own summary columns) — `info["n_blanked_px"]` specifically is
needed for Verification item 1c below. Save all of this to a scratch location. Everything
below is diffed against this.

**Step 1 — shared geometry helpers.** Add one small module, e.g. `midog_utils/geometry.py`,
with exactly the conversion primitives named in the landmine section: `center_of_closed_box`,
`box_from_center_closed` (its inverse), `center_of_half_open_box`. Every center derivation in
every later step must call through these — do not scatter ad hoc `(x1+x2)/2` math across 8
files. Unit-test each helper standalone, before touching any caller, against real numbers: a
`BOX_SIZE` GT box recovering the exact original click coordinate (closed), an Otsu half-open
bbox recovering the exact center `tightened_template_box` computes today (half-open), and a
round trip of `box_from_center_closed` → `center_of_closed_box` on a handful of known integer
centers.

**`box_from_center_closed`'s exact formula is pinned here, not left to judgment**:
`half = size // 2` (integer floor division, matching the existing `template_match.py:230`
convention exactly, not `size / 2.0`), `x1 = round(cx) - half`, `x2 = round(cx) + half`.
`cfg.base_size` is odd in every real production/template config (`_odd()` in both
`seed_selection.py` and `template_match.py` forces this) — `size / 2.0` true division
produces half-integer box corners for odd sizes, which then round to a **53×53 blanked
footprint instead of 51×51** once `blank_seed_square` converts the box back to array
indices (verify this yourself: `round(cx - 25.5)` and `round(cx + 25.5)` do not bracket a
51px span under Python's banker's rounding). This is a real, silent corruption of the
self-hit-blanking footprint, not just a half-pixel center wobble — and it will not be caught
by an unbalanced unit test built only around `BOX_SIZE=50` (even), which is exactly why the
formula is specified explicitly here rather than left for you to derive. Test the helper
against at least one odd size (e.g. 51) specifically, confirming the round trip AND that the
resulting box has width/height exactly `size` (not `size ± 1`).

**Step 2 — `dataset.py`.** **Additive, not a replacement**: keep `cx,cy,w,h` exactly as they
are today, and add `x1,y1,x2,y2` alongside them (`dataset.py:58-68`). This is deliberately
backward-compatible — `load_annotations` is a widely-used shared loader (the ~130 downstream
consumers this task doesn't migrate, plus `rerun_bbox3way_postD11.py:294`'s own assertion,
all read `cx`/`cy`/`w` today) and breaking it is unnecessary cost for zero benefit to this
task's actual goal, which is the *detections* representation (Step 5), not GT loading.

**Step 3 — `seed_selection.py`.**
- `tightened_template_box` (L105-140): return shape becomes `(base_size, box)` — `base_size`
  computed exactly as today (`_odd(max(y1-y0, x1-x0), minimum=minimum)`, L134, unchanged),
  `box` the **translated Otsu component bbox** in global coordinates — i.e. take the exact
  `(y0,y1,x0,x1)` shape `tighten_box_otsu` already returns (half-open, generally non-square,
  patch-local) and translate its origin to global coordinates. Do **not** build a synthetic
  `base_size`-square box around a rounded center — that would introduce a third, undocumented
  convention and reintroduce the landmine's 0.5px bug via rounding. This change **removes**
  the existing `-1` center-reconstruction math (L138-139); it does not add a step. The actual
  matched footprint used elsewhere stays derived from `template_xy` (now: `box`'s own
  half-open center) + `base_size` at each call site, exactly as today.
- `Seed` dataclass (L143-160): replace `click_xy`/`template_xy` with `click_box`/
  `template_box` (4-tuples). `click_box` is a **closed-convention** synthetic box (landmine
  item 3, built via `box_from_center_closed`); `template_box` is the **half-open** Otsu bbox
  from the point above (used as-is, already a box). Keep `base_size` as an explicit field
  (recommended over deriving it from `template_box`, since it's read directly by
  `production.py` and throughout `find_and_suppress.py` — removing it would touch more
  callers than keeping it). State in your report which you chose and why.
- `border_filter` (L47-62): derive `ix,iy` from the candidate pool's box columns via
  `center_of_closed_box`, once, vectorized over the whole pool.
- `_patch_readable` (L163-177): derive the scalar center via `center_of_closed_box` wherever
  `cx,cy` is currently read directly (L209 particularly) — this is a **closed**-box read (the
  raw candidate annotation).
- `build_seed` (L180-233): two *different* conventions apply here, do not use
  `center_of_closed_box` for both. The candidate draw (`cx, cy = row["cx"], row["cy"]`, L209)
  is closed-box, as above. But the two uses of `tightened_template_box`'s result — the
  border-readability check (`_patch_readable(roi_shape, spec[1], spec[2], patch_size)`,
  ~L212) and `offset_px = float(np.hypot(tx - cx, ty - cy))` (~L220) — read the **template**
  center, which is now half-open-box-derived: use `center_of_half_open_box` for `tx, ty`
  there, mirroring how Step 6 already treats an analogous half-open-box-to-center read. Get
  this wrong and `offset_px` silently corrupts by ~0.5px — this is exactly the sentinel
  Verification item 4 checks, so a mistake here should surface there, but don't rely on that;
  get it right the first time.

**Step 4 — `template_match.py`.**
- `blank_seed_square` (L206-238): change signature to accept a box directly instead of
  `(cx, cy, base_size)`. Lines 230-231 (`half = base_size // 2; ix, iy = round(cx), round(cy)`)
  are replaced by unpacking the incoming box. Lines 232-234 (reading `channel.shape` and
  clipping to array bounds) must be **kept** — the function still needs to clip whatever box
  it's given against the real array bounds; don't delete this part.
- `read_padded_patch` (L56-72): signature unchanged — it genuinely needs a literal center
  (its window is inherently symmetric about a point). Callers now derive that center via
  `center_of_closed_box`/`center_of_half_open_box` (whichever matches their box's own
  convention) before calling it.
- `extract_peaks`, `fused_response`, `build_augmentations`, `plant_and_recover`,
  `robust_stats`: unchanged. None of these touch box/center representation. Confirm this
  remains true after your edits.

**Step 5 — `find_and_suppress.py`** (the one hot-path file; two changes, not one). Verify
explicitly — with a real check in your own test code, not just by inspection — that
`centers`/`scores` flowing from `extract_peaks` (L105) through `nms_by_distance` (L108)
remain **completely unmodified, raw center arrays**.
1. L82: `tm.blank_seed_square(img_channel, seed_x, seed_y, cfg.base_size)` must change to
   build a box first (`box_from_center_closed(seed_x, seed_y, cfg.base_size)` or similar) and
   pass that box, matching Step 4's new signature. Easy to miss if you only look at the final
   DataFrame assembly — this call happens earlier in the same function and will `TypeError`
   immediately if left as-is.
2. L117-125: the DataFrame assembly — emit `x1,y1,x2,y2` (derived from each surviving center
   via `box_from_center_closed(cx, cy, cfg.base_size)`) instead of `cx,cy`.
`seed_x, seed_y` (L57) stays a literal center argument to this function; the caller derives
it once (Step 6).

**Step 6 — `production.py`.** L51: derive `seed_xy` (a literal center, for
`find_and_suppress`'s `seed_xy` argument) from `seed.template_box` via
`center_of_half_open_box`, once. L59 (`chromatin_od` branch): shift the detections' box
columns by `pad`, not `cx,cy`.

**Step 7 — `chromatin.py`.** `score_detections` (L43-47): derive `cx,cy` from the detections'
box columns via `center_of_closed_box` (vectorized), before calling `chromatin_density` per
row (which keeps taking a literal center — irreducible, since it calls `read_padded_patch`).

**Step 8 — `evaluate.py`.** `bucket_detections` (L68-104): derive `det_xy`/`gt_xy` from box
columns via `center_of_closed_box` **once** (L78-79), and stash the derived `cx,cy` onto the
returned `out`/`gt_out` frames so `evaluate_run`'s separate `coverage_fraction` call
(L181-184) reuses them instead of re-deriving from boxes a second time. `greedy_match`,
`coverage_fraction` themselves: signatures unchanged, still take raw point arrays —
irreducible.

**Step 9 — `viz.py`.** `overlay` (L18-61): draw GT boxes directly from `x1,y1,x2,y2` (removes
the `half`/`BOX_SIZE` reconstruction at L32,38). Detection circle markers still need a
derived center (a circle is inherently point+radius) — derive once via
`center_of_closed_box`, vectorized, before the plotting loop. `draw_box` (L64-77): unchanged,
already box-native.

**Step 10 — stale docstrings.** `Seed`'s class docstring (`seed_selection.py:147`,
"`click_xy` and `template_xy` can differ"), `build_seed`'s ("compute every match radius
against `seed.click_xy`, never `seed.template_xy`," `:197`), `tightened_template_box`'s
("ground truth stays keyed to `(cx, cy)`, never the returned centre," `:125`), and
`production.run_production_pipeline`'s ("reads only `template_xy` and `base_size`,"
`production.py:35`; its return-value doc listing `cx, cy`, `:40`) all name pre-refactor
fields/columns. Update these to the new field/column names — they're load-bearing comments a
future reader will trust.

**Step 11 — `invariants.py`, `nms.py`.** Confirm unchanged — don't assume, check.

## Verification — this is the actual bar for "done," not "it runs"

1. **Exact equivalence against Step 0's baseline**, in two parts:
   - **1a — GT, exhaustive.** For all 26,286 real annotations, assert
     `center_of_closed_box(x1,y1,x2,y2) == (original cx, original cy)` exactly. This is cheap
     (pure arithmetic, no image I/O) and is the single most direct test of the landmine — run
     it on the real dataset, not a handful of examples.
   - **1b — detections, per baseline run.** For every detection in Step 0's saved output,
     `center_of_closed_box(x1,y1,x2,y2)` must equal the baseline's original `(cx,cy)` exactly,
     or within float epsilon — state which, and why, if not bit-exact.
   - **1c — blanked footprint size.** For every run, `info["n_blanked_px"]` must equal
     `base_size**2` exactly (per `blank_seed_square`'s own docstring guarantee, unless the
     seed sits within half a template of the ROI edge, which shouldn't occur for a readable
     seed) and must match Step 0's saved baseline value exactly. This is the direct catch for
     an odd-size `box_from_center_closed` bug (see Step 1's pinned-formula note) — a wrong
     rounding choice changes this number even when 1a/1b still pass, since GT boxes are
     always even-sized and wouldn't otherwise exercise the bug.
2. **Identical NMS survivor set** — same detections by rank and identity, not just by count.
   This is the sharpest test that `nms_by_distance` saw byte-identical center/score arrays.
3. **Identical driver-level metrics.** `rerun_bbox3way_postD11.py` doesn't call
   `evaluate_run` directly — it reports its own metrics (`tp_at_k` and friends; read the
   script to get the exact names before writing this check, don't assume). Every one of
   those, for every ROI/seed the driver covers, must match baseline exactly. If you also want
   full `evaluate_run` parity (`recall_at_k`, `lookalike_attraction_rate`,
   `topk_composition`, `coverage_frac`), that requires adding a new call the driver doesn't
   currently make — treat that as optional extra evidence, not a requirement, since it goes
   beyond "adapt read sites."
4. **`Seed`-level sentinel.** For each seed drawn during the baseline run, diff the
   refactored pipeline's `click_box`/`template_box`-derived centers against the baseline's
   `click_xy`/`template_xy`, and independently recompute `offset_px` from both and compare.
   `offset_px = hypot(template_center - click_center)` is the one existing quantity in this
   codebase that already mixes the closed (`click_xy`) and half-open (`template_xy`)
   conventions — the cheapest available end-to-end sentinel for a convention mix-up, and not
   otherwise exercised by items 1-3.
5. **Run under an augmented config too**, not just production's `scales=(1.0,), n_angles=1,
   flips=(False,)` default — this repo's convention requires cleanups to test an augmented
   config. Specifically confirm `find_and_suppress.py`'s `aug` metadata lookup (L116, indexes
   `centers` directly, runs **before** your Step 5 box-conversion point) still executes on
   raw centers, unaffected by the box conversion that happens after it.
6. **Smoke check, not primary evidence:** `check_nms_radius` still passes, unmodified, at
   production config. (This is a scalar comparison unrelated to coordinate representation —
   it cannot meaningfully fail due to this refactor, so its passing is not evidence of
   correctness; only items 1-4 are.)
7. **The original working directory is untouched** — `git status --short` there, checked at
   the end, must match what it was before Setup step 1.

## Report back

- The landmine's four sites and which helper each uses, with Step 1's helper unit tests' and
  Verification item 1a's actual results (numbers, not "passed").
- Baseline-vs-refactored diff results for Verification items 1b, 1c, 2, 3, and 4, with real
  numbers.
- Confirmation the augmented-config test (item 5) passed, with what config you used.
- Confirmation of item 7 (original directory untouched).
- Which fields `Seed` ended up with, and why (the `base_size` decision from Step 3).
- What you changed in `rerun_bbox3way_postD11.py` and why, plus confirmation you touched no
  other file outside the 10-file `midog_utils/` set and `geometry.py`.
- One line each, at most, for anything else you noticed but didn't touch or fix.

## Don't

- Don't touch `find-and-suppress-midog` or anything in the original working directory.
- Don't change NMS or evaluation algorithm/semantics.
- Don't migrate downstream research scripts/notebooks (the one named exception is
  `rerun_bbox3way_postD11.py`, and only its read sites, per "Regression driver" above).
- Don't touch or attempt to fix `chromatin_od_ranker_seed_robustness_audit.py` — it's already
  broken, unrelated to this task, not your problem.
- Don't rewrite `DECISIONS.md` D7-D11.
- Don't add IoU-based NMS or box-overlap logic anywhere.
- Don't skip Step 0's baseline capture, or Verification item 1a's exhaustive GT check —
  without them, "verification" is just "it ran without crashing," which is not the bar here.
- Don't commit anything unless asked.
