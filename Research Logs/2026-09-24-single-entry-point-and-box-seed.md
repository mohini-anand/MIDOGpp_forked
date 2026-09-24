# One entry point for every seed source, and a drawn box whose size becomes the Otsu window: clicks and MIDOG runs unchanged bit for bit, 0 divergences on the full annotation-frame harness

**Scope.** Implements `SINGLE_ENTRY_POINT_PLAN.md` (reviewed adversarially in two rounds, then
settled with the user on 2026-09-24). Branch `feature/annotation-frame`; every comparison below is
against `d224c80`. Code: `midog_utils/annotation_frame.py` (`gate_window_side`, `box_frame`,
`_template_choice(..., window)`, the input-row size check), `production_pipeline/run_frame.py`
(`--source`, `visualize`), and a docstring-only deprecation note in `run_pipeline.py`. Decision
record: `DECISIONS.md` D12. Spec: `ANNOTATION_DATAFRAME_SPEC.md` §0–§4, §8, §9.

---

## What changed

- **One entry point.** `run_frame.py ROI --source {click,box,midog}`. `click` takes `--click X Y`
  (pixel-centre unless `--edge-origin`), `box` takes `--box X Y W H` (pixel-edge, as a canvas or
  COCO gives it), `midog` takes `--ann-id` or `--seed-index` (joint-gate draw, default seed 0).
  `--click-id` is now `--object-id`. A wrong-source option, a missing required one, or an invalid
  box is an argparse error before anything loads. `--save-fig` / `--plot` draw the overview plus a
  full-resolution zoom with the gate window (white) and the template footprint (cyan).
- **A drawn box is a seed.** Its centre is the click; its longer side, as a pixel count raised to
  odd, is the Otsu gate's window and the fallback template (`gate_window_side`). The stored row is
  snapped to whole pixels (`box_frame`), so it is the window and the fallback footprint. Windows
  must be 19–73 px, checked on every human input row.
- **Clicks and MIDOG rows are 50 wide, so their window is 51** — the old constant. That is why
  every existing result should be unchanged, and the checks below test exactly that.
- **`run_pipeline.py` is deprecated, not deleted.** It is still the only CLI for `build_seed`'s
  gray-only draw (300.tiff seed 0: 14581 there, 14499 in `run_frame.py`), and the collaborator repo
  ships it. Delete it once that repo has the new CLI.

## Checks, and what they showed

**V-A. Before/after differential.** The `d224c80` `annotation_frame.py`, loaded beside the new one,
and both `run_roi`s run on 19 cases: 300.tiff/14499 under all three refinement channels × both rank
keys; the four V5 gate refusals (462/22673 gray, 529/24934 gray, 300/14581 hem, 245/6274 hem); the
V5 `template_unreadable` click (300.tiff (6410, 250), both channels); the V6 half-integer clicks
(072/1835 gray and hem, 007/144 gray and untightened); a `click_frame` row; a two-row input seeded on
either row. **19/19 identical** on every frame column but `current_time` and every non-timing info
key, with the input frame unmutated.

**V-B. The existing harness, unchanged.** Smoke tier (300, 007) and the named smoke checks: 0
divergences; both regenerated per-ROI pickles identical to the committed ones on every key but
timings. The only change in `named.json` is the edge-error wording ("click (...)" → "seed centre
(...)"), on the four V7 edge cases; every other V7 raise has the same message as before, so each still
fails for its original reason and not for the new size check.

The full tier (49 ROIs, 441 runs, three concurrent shards, with `annotation_frame.py` checksummed
at launch and unchanged at the end): **0 divergences over 960,825 compared values**; all 49
regenerated pickles identical to the committed ones on every key but timings; the aggregate summary
identical apart from `v2_git`, which was not run (`--git` cannot pass for any base once
`annotation_frame.py` changes). The named full checks (33 checks) and `extra_checks.py`: 0
divergences, output identical to the committed JSON apart from the same four edge messages. In
place of V2's repository half: `git diff --name-only d224c80` over `midog_utils/`,
`production_pipeline/`, the driver and `midog_utils_full/` lists only the three intended files, and
`ast_equal.py` proves `run_pipeline.py`'s change docstring-only. The regenerated pickles and
summaries were restored from git (identical but for timings); the two `named.json` files are kept
as regenerated, so the tracked outputs carry the new edge wording.

**V-C. `validation/annotation_frame/box_checks.py`** (`out/box/box.json`), 0 divergences:

- `gate_window_side` on 23 frame-width vectors (e.g. 50 → 51, 51 → 53, 72 → 73, 17 → 19, a 50 × 16
  box → 51) and 8 raising ones (73, 16, 15 × 16, 80, NaN, inf, negative).
- `box_frame`: fixed boxes, 7 raising boxes, 999 shift cases (widths 18–73 in 0.5 steps × 8
  origins): the stored size depends only on the drawn size and the stored centre is the pixel
  containing the drawn centre; a 51 px box equals the click at its centre pixel for integer origins.
- The gate runs on the box's window: 7 boxes on 7 different annotations (even-pixel, non-square,
  half-integer and fractional origins, one hematoxylin, windows 31–73), each verified beforehand to
  give a **different template size** than the 51 px window would (e.g. a 31 px box on 14462: 29 px
  template where the 51 px window gave 51). `base_size`, template centre, refined row and
  `n_blanked_px` all match an independent `tightened_template_box` call.
- Fallback: a refused 41 px square box and a refused 25 × 57 box fall back to the human's box (the
  enclosing 57 px square for the non-square one), with predictions identical to
  `refinement_channel=None` on the same box; a 41 px box at V5's edge click gives
  `template_unreadable` and falls back the same way.
- A 73 px box centred exactly 36 px from the left or bottom edge runs; the new window guard raises
  `RuntimeError` when called directly on an unreadable window, and never fired on valid input.
- Width 80, width 16 and NaN rows raise, as the seed or as the non-seed row of a two-row input.

**V-D. CLI.** The four Step 0 baselines re-run under `--source` (`midog --ann-id 14499`, `midog`
(draw seed 0), `click 1411 3551`, and `midog --ann-id 14499` with hematoxylin and `chromatin_od`)
give frames identical to `d224c80`'s (read back with `round_trip`, `current_time` excluded) and the
same summary apart from the new `source` column. Thirteen bad invocations exit with argparse errors
before loading. `run_pipeline.py 300.tiff --ann-id 14499` prints what it printed before, and its
`--help` shows the deprecation note. Figures inspected by eye: MIDOG 14499 (51 px window, 49 px
template), a 30 × 44 box (45 px window, 27 px template), and a refused background box (the 41 px
box is the template).

**Probe.** `box_window_probe.py` → `results/box_window_probe.csv`; numbers in D12. At margin 0,
45% of accepted templates come out smaller than at 51 px and none larger; at margin 15, 59% come
out larger. Acceptance stays at 95.5–96.7%.

## What this does not show

- **Precision with box-sized windows was not measured.** Every precision number in the repo is
  at a 51 px window. Windows 53–73 px and templates that large are unvalidated.
- The probe's boxes are simulated (square, centred on the component), gray only.
- `run_pipeline.py` still calls `plt.show()` without `--save-fig`; under the macOS backend that
  opens a window and blocks, so its headless check ran with `MPLBACKEND=Agg` (no code change).

## Follow-ups

1. Port `run_frame.py`'s CLI to the collaborator repo (`~/Desktop/RA Work/AnnotateDx`), then delete
   `run_pipeline.py` and update `SELF_HIT_MASKING_PLAN.md:673,714` and
   `bbox_refinement_three_way_chromatin_od_audit_round3.py`, which call it.
2. Measure precision@K for drawn-box windows (19–73 px) against the 51 px window at the same clicks
   before relying on box-sized templates.
3. The planned radius-based refinement around the box centre; the client must keep the raw drawn
   box if it needs it, since the frame stores the snapped one.
