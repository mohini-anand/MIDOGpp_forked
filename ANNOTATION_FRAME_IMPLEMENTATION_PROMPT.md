# Task: implement the annotation frame specified in `ANNOTATION_DATAFRAME_SPEC.md`

## What governs

`ANNOTATION_DATAFRAME_SPEC.md` (the spec) is the specification. This prompt adds only four things:
the edit scope, the setup, the order of work, and the report format. It deliberately does not
restate the spec.

- Where this prompt and the spec disagree, the spec wins.
- Where the spec is ambiguous, or contradicts the code or data, **stop and report**. Don't choose
  an interpretation yourself.
- `CORNER_COORDINATE_REFACTOR_PROMPT.md` is superseded. Don't execute it or take instructions
  from it.

In one paragraph: add one module, `midog_utils/annotation_frame.py`, holding `run_roi` (spec §3–§5),
`evaluate` (§6) and `draw_seed_rows` (§7.1). Add one keyword-only argument, `always_od=False`, to
`production.run_production_pipeline` (§0). Change nothing else in the package. Then prove every
check in §7.2 (V1–V10) against the committed references and a baseline you capture from the
unedited code.

## Edit scope

**You may create or edit only:**

- `midog_utils/production.py`, and only to add `always_od` exactly as §0 describes. Four places
  change:
  - the signature;
  - the `if` condition of the `chromatin_od` block;
  - the `info["od_computed"]` line;
  - the docstring, gaining a line for the new argument and an updated Returns line (today it says
    `od` appears only for `chromatin_od`).

  Nothing else in that file changes.
- `midog_utils/annotation_frame.py`, which you create.
- Validation harness scripts and their outputs, under `validation/annotation_frame/` in the
  worktree. This is a new, untracked directory. Outputs go in `validation/annotation_frame/out/`,
  never under `results/`.

**Read-only — everything else**, in particular every file §0 freezes.

- `rerun_bbox3way_postD11.py` may be **imported** for its `select_clicks`, `template_spec`,
  `template_digest`, `compare_selection`, `build_rois` and `score_top_k`.
  - Run with the worktree root as the working directory, and do
    `sys.path.insert(0, ".")` before `import rerun_bbox3way_postD11 as rr`.
  - Importing it has only harmless side effects: a `sys.path` insert, two pandas display options
    and a module-level assert. It does no file I/O.
- Never edit it, never reassign its module globals, and never call `rr.main()`.
  - `main()` writes four tracked files: R1, R2, `…_postD11_per_roi.csv` and
    `…_postD11_selection_gate.csv`. It reads R3.
  - `fill_geometry_precision_at_10_49roi_6seed.ipynb` also imports this file.

**Style.** The new module follows the package's own conventions; `midog_utils/seed_selection.py`
and `midog_utils/dataset.py` are the models:

- a docstring with a brief summary, one line per parameter giving its type and a plain
  description, and one line for the return value;
- `def` lines and calls on a single line, never wrapped, however many arguments there are.

Per spec §1, every function that converts between `width` and `base_size` states the `+1` in its
docstring.

## Setup

1. Work from the original repo, on branch `find-and-suppress-midog`. **Precondition:** the spec
   and this prompt must already be committed there, because the worktree is created from `HEAD`.
   Check `git log -1 -- ANNOTATION_DATAFRAME_SPEC.md ANNOTATION_FRAME_IMPLEMENTATION_PROMPT.md`,
   and check that `git status --short` lists neither file. If either is modified or untracked,
   stop and ask.

   Record the output of `git status --short`. Then run
   `git worktree add ../MIDOGpp_forked_annotation_frame -b feature/annotation-frame`, and record
   `git rev-parse HEAD`: this is `<base>` for V2.
2. In the **original** directory, confirm `git status --short` is unchanged. If anything differs,
   stop and report.
3. **ROI images.** `images/` is not empty in a fresh worktree: `images/.gitkeep` and the two
   `MANIFEST.md` files are tracked. So `ln -s <orig>/images ./images` would create
   `images/images` and leave every ROI missing. Don't do that. Symlink the TIFFs file by file into
   the tracked subset directories, from the worktree root:

   ```
   ORIG=/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked
   for d in images/extra_valid images/extra_valid/testing_set; do for f in "$ORIG/$d"/*.tiff; do ln -s "$f" "$d/"; done; done
   ```

   Then confirm all of the following:
   - `ls images/extra_valid/*.tiff | wc -l` is 14;
   - `ls images/extra_valid/testing_set/*.tiff | wc -l` is 35;
   - `git status --short` in the worktree does not list them, because `.gitignore` covers
     `images/**/*.tiff`;
   - `rr.build_rois()` returns 49 ROIs without raising.
4. `databases/MIDOG++.json` and the reference CSVs R1–R4 (spec §7.2) are tracked. Confirm with
   `git ls-files` in the worktree; don't assume.
5. Use `/Users/mohinianand/anaconda3/bin/python3` for everything. The bare `python3` cannot import
   `cv2`.
6. Every run starts with the worktree root as the working directory, because `rr`'s paths are
   relative.

## Read first, fully, before editing anything

- The whole spec, including §9, which lists which parts are corrections and which are decisions.
- `midog_utils/production.py`, `find_and_suppress.py`, `seed_selection.py`, `template_match.py`,
  `chromatin.py`, `evaluate.py`, `dataset.py`, `channels.py`.
- `rerun_bbox3way_postD11.py`: all of `template_spec`, `template_digest`, `pool_gate_specs`,
  `draw_joint_click`, `score_top_k`, `select_clicks`, `run_roi`, `compare_selection` and
  `build_rois`. `draw_seed_rows` ports the draw, and your harness reuses the rest.
- `DECISIONS.md` D7, D9 and D11, and `D8_TEMPLATE_ANCHOR.md`, for why the frozen code looks the way
  it does.

## Order of work

**Step 0 — baseline, before any edit.** Capture B1 and B2 (§7.2) from the unedited code.
- Each run's seed is built exactly as `rr.run_roi` builds it (`:232-243`): `rr.select_clicks` for
  the draw, then `ss.Seed(...)` from each condition's spec.
- Save every survivor and the full `info`, as pickles, under `validation/annotation_frame/out/`.
- In the same pass, capture the V9 augmented-bank outputs.
- Then show that B1 reproduces R1 and R2. B1 is raw production output, so first rebuild each run's
  R1 row and R2 rows from it with the driver's own code: `rr.score_top_k` for the buckets and
  `tp_at_k`, and the formulas at `rerun_bbox3way_postD11.py:233-252` for everything else. Use
  §7.2's rules: exact comparison, float32 `score`, NaN equals NaN, strings `fillna("")`, `t_*`
  excluded, `round_trip` reads.
- Also show that B2's 42 `default_51` runs (14 ROIs × 3 seeds) reproduce R4's `tm_score` top 30
  on the columns the spec's §7.2 baselines paragraph lists.

If B1 does not reproduce R1 and R2, or B2 does not reproduce R4, **stop and report**. The
environment is wrong, and nothing after this step would mean anything.

**Step 1 — `always_od`.** Make the one edit to `production.py`. Show the diff.

**Step 2 — the module.** Implement `annotation_frame.py` to the spec. Before any long run, test the
pieces on one ROI:
- `draw_seed_rows` against `rr.select_clicks`;
- the §2 conversions;
- `run_roi` for each refinement setting.

**Step 3 — smoke tier.** Each check has its own inputs, so the smoke tier names them. Fix anything
it finds, then proceed.

- V1–V4 and V10 on `300.tiff` (original_14) and `007.tiff` (testing_35).
- V5 on its named clicks: 462, 529, 300 and 245, plus the `(6410.0, 250.0)` click on 300.
- V6 on 072 and 007.
- V7 on 300.
- V8 on one all-int ROI (211, 364 and 529 are among the 49).
- V9 on one of its 5 crops.
- V1's contested-tier test on its chosen ROI.

Every one of these is small except V1–V4 and V10. Those cost about 36 pipeline calls on
`300.tiff` (18 `run_roi` calls, plus V2's 18, since 300 is an original_14 ROI) and about 18 on
`007.tiff`.

**Step 4 — full tier.** Run V1–V10 on all 49 ROIs. To avoid re-reading each ROI, organise the
post-edit pass per ROI:
- load the ROI once;
- run the draw check (V1);
- run `run_roi` under both rank keys for every seed and condition (V3, V4);
- evaluate (V10);
- write per-ROI results to disk as you go, so an interruption loses at most one ROI.

Run long passes in the background with a log file.

**Runtime.** This machine is a CPU-only Intel Mac. R1's recorded pipeline time is 43 min for 441
runs, about 5.9 s each. That includes production's own hematoxylin conversion (1.6–3.3 s) and
excludes TIFF loads and the draw's channel conversions. `run_roi` adds a second hematoxylin
conversion on `hem_bbox` runs, about 10 min over the full tier.

| pass | calls | estimate |
| --- | --- | --- |
| Step 0: B1 and the 49 draws | 441 | ~1 h |
| Step 0: B2 | 126 | ~15 min |
| Step 4: V3 and V4 (`run_roi`, both rank keys) | 882 | ~1.5–2 h |
| Step 4: V2 (the 126 `tm_score` runs twice) | 252 | ~25 min |
| everything else | — | minutes |

Budget about 4 h of compute. Don't cut the run to fit: the full tier is the bar.

## Specifics the spec leaves to you, decided here

- **V5 and V6 clicks.** The spec lists candidate clicks. Re-verify each before relying on it:
  - V5: the named gate refuses it and the other accepts it;
  - V6: it is half-integer and both gates accept it.

  If one fails, pick a replacement with the same property from the same source (R1's
  `refused_draws` for V5; half-integer mitotic annotations in the 49 ROIs for V6), and report the
  substitution.
- **V7 edge clicks.** Use a synthetic human row on `300.tiff`: a click at `x = 35` raises; at
  `x = 36` it does not. Keep `y` well inside the image. Repeat at the right edge with
  `x = w - 36` (raises) and `x = w - 37` (does not).
- **V8 all-int image.** Take any of the 87 images whose JSON boxes are all Python ints. Build a
  one-row human frame from its raw JSON with int64 columns, and check the cast.
- **V9 crops.** Take 5 ROIs from `original_14`. Crop 1024 × 1024 around each ROI's seed-0 click,
  clamped inside the image. Call `find_and_suppress(crop, click_in_crop, FSConfig(scales=(0.8, 1.2), n_angles=2, flips=(False, True), deep_floor_z=-1.5, nms_radius=ev.radius_px(mpp), max_peaks=100))`
  on the hematoxylin channel. Every column and every non-`t_*` `info` value must be identical
  before and after the edit.

## Report back

1. The setup confirmations: `<base>`, the TIFF counts, the clean `git status` in both directories,
   and `rr.build_rois()` returning 49.
2. The `production.py` diff, verbatim.
3. Step 0: whether B1 reproduced R1 (441 rows) and R2 (13,230 rows), and whether B2's 42
   `default_51` runs reproduced R4's `tm_score` top 30 (1,260 rows). Give the number of values
   compared and the number that diverged for each.
4. For each of V1–V10: pass or fail, what was compared, and counts — runs, rows, values and
   divergences. For example, "V3: 441/441 runs, N prediction rows, 0 divergences". Numbers, not
   "passed".
5. Every file you created, with its path.
6. Any spec ambiguity you hit, and where you stopped because of it.
7. At most one line each for anything else you noticed but did not touch.

## Don't

- Don't edit any file outside the edit scope. That includes the spec, `DECISIONS.md`, the
  reference CSVs, `rerun_bbox3way_postD11.py` and every frozen `midog_utils` file. Propose
  changes in the report instead.
- Don't re-sort predictions, round a coordinate outside the one site §1 names, or recompute
  anything template-relative from the frame (§1, §2).
- Don't loosen any comparison into a tolerance. A divergence is a finding.
- Don't touch the original working directory or the `find-and-suppress-midog` branch once the
  worktree exists.
- Don't commit anything or delete the worktree unless asked.
