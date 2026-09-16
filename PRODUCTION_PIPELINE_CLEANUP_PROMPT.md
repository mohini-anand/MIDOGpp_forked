You are removing dead code from the MIDOG++ click-to-verify production pipeline at
`/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked` (a git repo). Work in the working tree.
Behaviour of the production pipeline and of the template-augmentation machinery must not change by
a single bit. Being slow and stopping early is fine; improvising is not.

## Decisions already made (do not re-decide these)

1. **Harness location:** `HARNESS_DIR = /Users/mohinianand/Desktop/AnnotateDx/cleanup_harness/`
   (outside the repo). Everything you create goes there, and nothing there gets deleted at the end.
2. **Stop after Step 2** (harness built, baseline captured, trust checks run). Report the results
   and wait for the go-ahead before editing any code.
3. **`FSConfig.deep_floor_z`:** keep it as `Optional[float] = None`, so `FSConfig()` still
   constructs. `find_and_suppress()` raises `ValueError` when it is `None`, since the
   `score_threshold` fallback is being removed.
4. **`seed_selection` parameters:** remove `method`, `center_tolerance` and `headroom_frac` from
   `tighten_box_otsu`, `tightened_template_box` and `build_seed`, and remove `recentre` from
   `build_seed`. Keep the `Seed.recentred` field; it is now always `True`. Nothing on the
   production path passes these parameters, and the `production_hematoxylin_only/` notebooks only
   mention them in prose. Older research code does use them (e.g. `bbox_headroom_end_to_end.py`,
   `bbox_threshold_sweep.py`, `morphometric_diagnostic.ipynb`) and will break, which is accepted.
5. **Matching method:** only `cv2.TM_CCOEFF` is ever used. Remove the sign handling and add nothing
   in its place: no SQDIFF guard, check or special case anywhere. Every matching-method default
   becomes `cv2.TM_CCOEFF`: `FSConfig.tm_method` (B3), and the `method` parameter of `fused_response`
   and `plant_and_recover` (B4). Today all three default to `cv2.TM_CCOEFF_NORMED`.
6. **The augmentation machinery stays fully intact and working.** That covers:
   - `template_match.py`: `Augmentation`, `_odd`, `BASE_SIZE`, `PATCH_SIZE`, `_FLOOR`,
     `_SENTINEL_CUT`, `read_padded_patch`, `build_augmentations`, `_robust_z`, `robust_stats`,
     `fused_response` (multi-template max-fusion, the `best` index map, the `valid` mask, the
     `scale_normalize` branch), `extract_peaks` (including its own `score_threshold` parameter;
     only the `FSConfig` field goes) and **`plant_and_recover`**. `plant_and_recover` is kept even
     though the cleanup doc lists it for removal; this prompt overrides the doc on that.
   - `find_and_suppress.py`: `FSConfig.base_size`, `patch_size`, `scales`, `n_angles`, `flips`,
     `scale_normalize`, `border_pad`, the `n_augmentations` property; inside `find_and_suppress()`,
     the `build_augmentations` call, the `border_pad` pad derived from the largest template, and
     the `best` → `metas` lookup that fills the `angle`/`flip`/`scale` columns.

   Production runs a single template (`scales=(1.0,)`, `n_angles=1`, `flips=(False,)`), so the
   production harness alone can't prove any of this. The augmentation checks in Step 1 do.
7. **No commits.** Don't commit or push unless explicitly told to.

## Ground rules

- **Interpreter:** use `/Users/mohinianand/anaconda3/bin/python3` and `/Users/mohinianand/anaconda3/bin/jupyter`.
  The bare `python3` fails at `import cv2`. Run all commands from the repo root.
- **Files you may edit:** only these, under `midog_utils/`: `__init__.py`, `dataset.py`,
  `seed_selection.py`, `find_and_suppress.py`, `production.py`, `template_match.py`, `chromatin.py`,
  `invariants.py`, `evaluate.py`, `viz.py`, `channels.py`.
- **Files you may delete:** only `midog_utils/{baselines,experiment,compare}.py`.
- **Off limits:** no other file, no `.md` file, no notebook, nothing in `midog_utils_full/`.
- **Git:** read-only commands only (`status`, `diff`, `log`). Never `checkout`, `restore`, `stash`,
  `reset`, `clean`, `add` or `commit`. The repo has many untracked research files that some of
  those commands would destroy. Undo changes with the file checkpoints described below.
- **Deferred items:** don't touch `production.py`'s `rank_key`/`AXES`, `run_pipeline.py`'s
  `--rank-key` flag, or its `tumor_type` selection.
- **When code and instructions disagree, stop.** If the code contradicts this prompt or
  `PRODUCTION_PIPELINE_CLEANUP.md` (a symbol has a live caller neither mentions, a line isn't where
  it's described, an edit doesn't fit), stop and report. Don't work around it.
- **Locate code by symbol name.** Line numbers in the cleanup doc are from before any edit and will
  drift.
- **Style:**
  - Every `def` header and every call stays on one line, never wrapped.
  - Docstrings match the file's existing style.
  - When you delete something, update in the same batch any in-code docstring, comment or error
    message in a kept file that names it.
  - Remove imports that a batch leaves unused.
- **Expected breakage outside scope:** about 65 research scripts and notebooks outside the
  production path import `experiment`, `baselines`, `compare` or functions removed here. That
  includes `verify_fixes.py`, `nms_ordering_probe.py`, `tm_variant_sweep.py` and the notebooks in
  `production_hematoxylin_only/`. They will break, and that is accepted. Don't run them (except
  `verify_fixes.py` as described below), don't fix them, and list them in the final report.

## Read first

1. `PRODUCTION_PIPELINE_CLEANUP.md` describes what to remove. Act only on "Confirmed for removal".
   This prompt overrides it in these places, all checked against the code:
   - The batch order.
   - `plant_and_recover` is kept (decision 6).
   - The `check_no_cap` call in `production.py` goes in B3, not with `invariants.py`. It reads
     `cfg.max_detections`, so it has to go no later than that field.
   - `FSConfig.score_threshold` is read unconditionally at `find_and_suppress.py` line
     `threshold = cfg.score_threshold`, so removing the field needs an edit (B3).
   - Removing `to_hematoxylin`/`to_rgb` also means removing their `CHANNELS` entries (B10).
   - `tm_variant_sweep.py` will crash after B4 (it uses `tm.METHODS` and `tm.method_sign` directly);
     it won't silently produce wrong scores as the doc suggests.
2. `production_pipeline/run_pipeline.py`, `midog_utils/production.py`, and
   `production_pipeline/PRODUCTION_PIPELINE_WALKTHROUGH.md` (current behaviour). Parts of the
   walkthrough, such as its `check_no_cap` sections, will go stale. Don't change code to match it.
3. Only if you need the reasoning behind a fixed constant: `DECISIONS.md` (D3, D7, D8, D9),
   `D8_TEMPLATE_ANCHOR.md`.

## Step 0: pre-flight

- Run `git status --short midog_utils production_pipeline`. It must show no modified tracked files.
  If it does, stop and report.
- Create `HARNESS_DIR` with subfolders `runs/` and `checkpoints/`, and a `LOG.md` you append to
  after every step.

## Step 1: build the harness (`HARNESS_DIR/harness.py`)

The harness must:

- Insert the repo root into `sys.path` (`production_pipeline/` has no `__init__.py`; importing it as
  a namespace package works). Call `matplotlib.use("Agg")` before any pyplot import.
- Import only:
  - `from production_pipeline.run_pipeline import run_pipeline_on_roi, summarize, visualize`
  - `from midog_utils import evaluate as ev, channels as ch, template_match as tm`
  - `from midog_utils.find_and_suppress import FSConfig, find_and_suppress`
  - `cv2`, `numpy`, `pandas` and the standard library.

  Use these as they are; don't reimplement anything.
- **`capture <label>`** writes `HARNESS_DIR/runs/<label>/` with everything below.

  **Production part.** For each of the 14 sorted `.tiff` files in `images/extra_valid/`, and each
  `rank_key` in (`"tm_score"`, `"chromatin_od"`), call `run_pipeline_on_roi(fn, seed_index=0, rank_key=rank_key)`
  once. Write:
  - `full.pkl`, keyed by `(fn, rank_key)`, holding:
    - `detections` (DataFrame)
    - `gt_eval` (DataFrame)
    - `dataclasses.asdict(seed)`
    - `info` with every key starting with `t_` dropped (those are timings)
    - `mpp` and `match_radius`
  - `summary.csv`: one row per (fn, rank_key, top_n) for top_n in (10, 20, 30, 50), holding
    `summarize(result, top_n)` plus the seed's `ann_id`, `base_size`, `template_xy`, `offset_px`,
    `n_retries`.
  - `evaluate_run.pkl`: for every result, call
    `ev.evaluate_run(result["detections"], result["gt_eval"], result["match_radius"], area_mm2, roi_shape=result["rgb"].shape[:2])`
    with `area_mm2 = h * w * mpp**2 / 1e6` computed inline (`ds.roi_area_mm2` is removed in B1).
    Store the returned tuple's length, `det_out`, `gt_out` and `metrics`.
  - `overlay.png`: `visualize()` on one result, saved, then `plt.close("all")`. This is the only
    coverage of `viz.overlay`/`draw_box`.

  **Augmentation part** (`aug.pkl`). For each of `013.tiff`, `245.tiff`, `403.tiff`, reuse that
  ROI's `tm_score` result (`seed`, `mpp`, `rgb`):
  - `hem = ch.to_channel(result["rgb"], "hematoxylin_od")`. Cut a 1024x1024 crop around
    `seed.template_xy`, clamped inside the ROI:
    `x0 = min(max(round(tx) - 512, 0), w - 1024)`, same for `y0`, and `seed_xy = (tx - x0, ty - y0)`.
    Then `patch = tm.read_padded_patch(crop, seed_xy[0], seed_xy[1], tm.PATCH_SIZE)`.
  - **One 8-template bank** (2 scales x 2 angles x 2 flips; covers downscaling, upscaling, a
    rotation through `warpAffine`, and flips):
    `templates, metas = tm.build_augmentations(patch, seed.base_size, (0.8, 1.2), 2, (False, True))`.
    Store each template's shape and the sha256 of its bytes, and `metas` as plain tuples. Use this
    same bank, or the same settings, for everything below. Don't enlarge it.
  - For `sn` in (False, True): `fused, best, valid = tm.fused_response(crop, templates, sn, method=cv2.TM_CCOEFF)`.
    Store the shape, dtype and sha256 of each array.
  - For `sn` in (False, True): `find_and_suppress(crop, seed_xy, cfg, nms_radius=ev.radius_px(mpp))` with
    `cfg = FSConfig(base_size=seed.base_size, scales=(0.8, 1.2), n_angles=2, flips=(False, True), peak_min_distance=7, max_peaks=100, deep_floor_z=-1.5, border_pad=True, scale_normalize=sn, tm_method=cv2.TM_CCOEFF)`.
    Never pass `channel`, `score_threshold` or `max_detections`; they are removed in B3. Store the
    detections DataFrame and `info` without `t_*` keys.
  - For `sn` in (False, True): store the list returned by
    `tm.plant_and_recover(templates, metas, method=cv2.TM_CCOEFF, scale_normalize=sn)` as plain
    tuples. If it raises, that is a failure.
  - This bank doesn't cover arbitrary-angle rotation (its only non-zero angle is 180 degrees) or
    the unscaled 1.0 template (the production part covers that). That's an accepted trade-off for
    speed.

  **API part** (`api.json`):
  - `str(inspect.signature(f))` for `tm.read_padded_patch`, `tm.build_augmentations`, `tm._robust_z`,
    `tm.robust_stats`, `tm.fused_response`, `tm.extract_peaks`, `tm.plant_and_recover`,
    `find_and_suppress`.
  - `[(f.name, repr(f.default)) for f in dataclasses.fields(FSConfig)]`.
  - `FSConfig().n_augmentations`.
  - Whether `tm.Augmentation`, `tm.BASE_SIZE` and `tm.PATCH_SIZE` exist.
- **`compare <label_a> <label_b>`:** exact comparison. Exit nonzero on any difference not covered by
  the two exceptions below, and print every difference found.
  - DataFrames: `pd.testing.assert_frame_equal(..., check_exact=True)`.
  - Dicts, lists, tuples and hashes: exactly equal, with NaN treated as equal to NaN.
  - `evaluate_run.pkl` exception: `det_out`/`gt_out` exact; `metrics` values equal on shared keys.
    Separately print (a) keys only in `a`, (b) keys only in `b`, (c) both tuple lengths, without
    failing on (a) or (c).
  - `api.json` exception: print every difference without failing on it. The batch instructions
    list the only acceptable ones; anything else is a failure you must act on.

A capture takes about 2-2.5 minutes. Measured: 3-10 s per production call; the whole
augmentation part took about 13 s for all three ROIs, including ROI loading (under 0.5 s per
`fused_response`, `find_and_suppress` or `plant_and_recover` call on a crop).

## Step 2: baseline and trust checks (then stop and report)

1. `capture baseline`, then `capture baseline_repeat`, then `compare baseline baseline_repeat`.
   Expect zero differences. The production pipeline was measured bit-deterministic on 28/28 calls
   on 2026-09-16. Any difference means exact-match checking is unusable: stop.
2. **Augmentation sanity.** In `runs/baseline/aug.pkl`, confirm:
   - for every ROI and both `scale_normalize` settings, `plant_and_recover` returned 8 entries, all
     with `dx == dy == 0`;
   - for every ROI and both settings, the `find_and_suppress` detections use both angles (0 and 180)
     and both `flip` values;
   - across all six runs together, both scales (0.8 and 1.2) appear. A single run often uses only
     one scale; that is normal.

   Measured on 2026-09-16 for 013/245/403: every item above held. `pad_px` was 18/28/30, so the
   border padding is derived from the largest template, and it varies with the seed's `base_size`.
3. **External reference.** Compare `runs/baseline/summary.csv` to
   `threshold_maxpeaks_ablation/chromatin_od_ranker_seed_robustness_precision.csv` filtered to
   `seed_index == 0`. That CSV comes from an independent notebook at production config
   (hematoxylin_od, TM_CCOEFF, deep floor -1.5, MAX_PEAKS=100). Match on `file_name`/`roi`,
   `arm`/`rank_key` and `budget`/`top_n` for K = 10/20/30.
   - Expected: all 84 rows match, with `n_detections` exact and `precision_at_budget`/`recall_at_budget`
     within 6e-4 (`summarize` rounds to 3 dp).
   - K=50 has no external reference.
   - **Do not** use `results/precision_at_k_14roi_prodseed_chromatin_halfpixfix_per_roi.csv`. It ran
     with MAX_PEAKS=2,000,000, so chromatin_od re-ranked about 18k candidates instead of about 100
     and only 21/56 of its chromatin_od cells match. That mismatch is expected and says nothing
     about the harness.
4. Run `python verify_fixes.py > HARNESS_DIR/runs/baseline/verify_fixes.txt 2>&1` from the repo
   root. It must end with `FAILURES: none` (it did on 2026-09-16, in about 16 s).
5. Run `python -m pyflakes midog_utils/*.py production_pipeline/run_pipeline.py > HARNESS_DIR/pyflakes_baseline.txt`.
   The expected output is only pre-existing lines: ten `midog_utils/__init__.py ... imported but unused`
   lines, `baselines.py:23 '.channels.to_gray_inverted' imported but unused`, and
   `experiment.py:453 undefined name 'fused'`.
6. **CLI smoke test:** `python production_pipeline/run_pipeline.py 403.tiff --rank-key chromatin_od --save-fig HARNESS_DIR/runs/baseline/cli_403.png`
   must finish and print a summary row.
7. **Demo notebook:** `jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=3600 --output-dir HARNESS_DIR/runs/baseline production_pipeline/production_pipeline_demo.ipynb`
   must finish without error. Never use `--inplace`; that would rewrite a tracked file.
8. Report the results of items 1-7 and wait.

## Step 3: per-batch protocol (apply to every batch in Step 4)

1. **Checkpoint.** Copy every file the batch will edit or delete to
   `HARNESS_DIR/checkpoints/<batch>/<repo-relative path>` (`cp -p`).
2. **Apply** the batch's edits, and nothing else.
3. **Static gate** (fast; do this before the harness):
   - pyflakes over the same files. Allowed output: the baseline lines (ten `__init__.py` lines
     through B8, eight from B9 on), minus lines for files that no longer exist. Any new line must be
     fixed within the batch (usually an import the batch made unused) or treated as a failure.
   - For each removed symbol, grep for it as code (e.g. `def points`, `\.points\(`) in
     `midog_utils/*.py`, `production_pipeline/run_pipeline.py` and
     `production_pipeline/production_pipeline_demo.ipynb`. The only acceptable hits are inside
     `baselines.py`, `experiment.py` and `compare.py` before B9 (latent references that die with
     those files). Note those hits in `LOG.md`.
4. **Harness:** `capture <batch>`, then `compare baseline <batch>`. There must be zero differences in
   the production part and the augmentation part. The only acceptable exceptions are the
   `evaluate_run` key/length changes (from B7) and the `api.json` changes (from B3 and B4) that the
   batches below list.
5. **`verify_fixes.py`:** do what the batch says.
6. **If anything fails:**
   - Copy every checkpointed file back (this also restores deleted files).
   - `capture <batch>_restored` and `compare baseline <batch>_restored` to confirm you are green again.
   - If the cause isn't obvious, re-apply the batch's items one at a time (checkpoint, verify each)
     to isolate it.
   - Then stop and report exactly what diverged: keys, rows, values. Don't try fixes this prompt
     doesn't specify.
7. **If everything passes:** save `git diff > HARNESS_DIR/checkpoints/<batch>/after.patch` and append
   the batch, the gate results and `git diff --stat` to `LOG.md`.

## Step 4: the batches, in this order

**B1 `dataset.py`.** Delete `load_slide_metadata`, `canonical_scanner`, `SCANNER_ALIASES`,
`check_invariants`, `points`, `roi_area_mm2`, `check_roi_scale`, `CATEGORY_NAMES`.
Keep everything the cleanup doc lists as kept. Expect `verify_fixes`: `FAILURES: none`.

**B2 `seed_selection.py`.**
- Delete `foreground_filter`, `SeedInfo`, `pick_seed`.
- In `build_seed`: remove the `recentre` parameter and the `recentre=False` branch; always call
  `tightened_template_box`; construct `Seed(..., recentred=True, ...)`. Then delete
  `tightened_base_size`, whose only in-package caller was that branch.
- In `tighten_box_otsu`: keep only the two-class Otsu path. Remove the `method` dispatch and its
  `"multiotsu"`/`"headroom"` branches and the `ValueError` for unknown methods. Remove the
  `center_tolerance > 0` branch, then delete `_nearest_label_within`.
- Remove `method`, `center_tolerance` and `headroom_frac` from the signatures of `tighten_box_otsu`,
  `tightened_template_box` and `build_seed`, along with `build_seed`'s `kw` pass-through (pass
  `otsu_window=otsu_window` directly, on one line).
- Remove the unused `threshold_multiotsu` import.
- Update docstrings. Fix the `tighten_box_otsu` error message that names `channels.to_rgb` (it is
  removed in B10): say an RGB array was passed instead of a single structural channel.
- The seed fields in `full.pkl` must be identical. 403.tiff exercises the retry path
  (`n_retries=1`).
- Expect `verify_fixes`: `FAILURES: none`.

**B3 `find_and_suppress.py` + `production.py`.**
- `FSConfig.channel`: delete the field and its comment, and delete `channel=CHANNEL` from the
  `FSConfig(...)` call in `production.py`. Keep `production.CHANNEL`; it is still used for `to_channel`.
- `FSConfig.score_threshold`: delete the field. In `find_and_suppress()`:
  - Right after the existing `nms_radius is None` check, raise `ValueError` if `cfg.deep_floor_z is None`.
  - Replace `threshold = cfg.score_threshold` plus the `if cfg.deep_floor_z is not None:` block with
    the unconditional robust-stats computation, still writing the info keys `deep_floor_median`,
    `deep_floor_mad` and `score_threshold_used`.
  - Update the `deep_floor_z` field comment (it currently says `None = use score_threshold`).
- `FSConfig.max_detections`: delete the field and the line
  `centers, scores = centers[: cfg.max_detections], scores[: cfg.max_detections]`. Also delete the
  `inv.check_no_cap(...)` line in `production.py`. The `invariants` import stays for
  `check_nms_radius`.
- `FSConfig.tm_method`: change its default to `cv2.TM_CCOEFF`. Production already passes
  `tm_method=TM_METHOD` explicitly; leave that call as it is.
- Leave every augmentation field and code path from decision 6 untouched.
- Acceptable `api.json` differences, and no others:
  - the `FSConfig` field list loses `channel`, `score_threshold` and `max_detections`;
  - the `tm_method` default repr changes from `5` to `4`.
- Expect `verify_fixes`: `FAILURES: none` (E6 calls `FSConfig().patch_size`, so `FSConfig()` must
  still construct).

**B4 `template_match.py`.**
- Delete `METHODS`, `method_sign` and `_SIGN`. Do **not** delete `plant_and_recover`.
- In `fused_response`, delete `sign = np.float32(method_sign(method))` and the
  `if sign < 0: res *= sign` block, and add nothing in their place.
- Change the `method` default of `fused_response` and `plant_and_recover` to `cv2.TM_CCOEFF`.
- In both docstrings, replace "any cv2.TM_* method; the two distance methods are negated on the way
  in" with a plain description of the method parameter. Drop any other SQDIFF or sign wording in
  kept code.
- Leave every other augmentation function from decision 6 untouched.
- Acceptable `api.json` differences, and no others: the `method` default repr in the
  `fused_response` and `plant_and_recover` signatures changes from `5` to `4`. `aug.pkl` must be
  identical (every harness call passes `method`/`tm_method` explicitly).
- Expect `verify_fixes`: `FAILURES: none`.

**B5 `chromatin.py`.** Delete `rerank` and the now-unused `Optional` import. Update the module
docstring's sentence about `rerank`. Expect `verify_fixes`: `FAILURES: none`.

**B6 `invariants.py`.** Delete `check_no_cap`, `check_tissue_mask_covers_gt`,
`check_min_separation`, and the imports this leaves unused (pyflakes will name them: expect
`pandas` and `Sequence`). Keep `InvariantError`, `check_distinct_seeds`, `check_nms_radius`.
Expect `verify_fixes`: `FAILURES: none`.

**B7 `evaluate.py` + `viz.py`**, everything except `optimal_assignment_disagreement`.
- Delete `froc`, `sensitivity_at_fp`, `full_list_breakdown`, `threshold_sweep`.
- In `evaluate_run`: delete the `fp_mm2, sens = froc(...)` line, the `sens@...` `metrics.update`
  line and the `full_list_breakdown` `metrics.update` line. Change the return to
  `return det_out, gt_out, metrics`.
- Update the module docstring's FROC wording.
- In `viz.py`, delete `froc_plot` and `legend_handles`.
- Expected in `compare`:
  - `evaluate_run` tuple length goes 4 to 3.
  - Keys only in baseline are exactly these 18 (check the printed set, not a prefix rule):
    `sens@1fp_mm2`, `sens@2fp_mm2`, `sens@4fp_mm2`, `sens@8fp_mm2`, `sens@16fp_mm2`,
    `sens@32fp_mm2`, `sens@64fp_mm2`, `all_n_detections`, `all_tp`, `all_fp_lookalike`,
    `all_fp_unannotated`, `all_matched_any_gt`, `all_matched_frac`, `all_precision_mitotic`,
    `all_mitotic_gt_found`, `all_mitotic_gt_missed`, `all_lookalike_gt_found`,
    `all_lookalike_gt_missed`.
  - No keys only in the new run.
  - Every shared key identical.
  - Everything else identical.
- Expect `verify_fixes`: `FAILURES: none`.

**B8 `evaluate.py`.** Delete `optimal_assignment_disagreement` and the now-unused
`linear_sum_assignment` import. Expect `verify_fixes` to print PASS for every E1-E6 check, print the
E7 header, then crash with an `AttributeError` naming `optimal_assignment_disagreement` (don't
match the exact message text; Python may append a suggestion). Anything else (a FAIL before E7, a
different error, a crash earlier) is a real regression. `verify_fixes.py` is not used after this
batch.

**B9 whole modules.** Delete `midog_utils/baselines.py`, `experiment.py` and `compare.py`, and, in
the same batch, remove `baselines` and `experiment` from `midog_utils/__init__.py`'s import list.
`experiment.py`'s augmentation experiment drivers (`augmentation_variants`, `fusion_variants`,
`scale_usage`) go with it. They are research drivers built on the removed `FSConfig.channel`, not
part of the augmentation machinery in decision 6. Then:
- `python -c "import midog_utils, midog_utils.production, midog_utils.invariants, midog_utils.chromatin, midog_utils.nms"`
  must succeed.
- `grep -nwE "baselines|experiment|compare" midog_utils/*.py` may only hit generic prose (today:
  `__init__.py`'s "histology experiment" and `invariants.py`'s "this experiment relies on").

**B10 `channels.py`.** Delete `to_hematoxylin`, `to_rgb`, their `"hematoxylin"` and `"rgb"` entries
in `CHANNELS`, and the now-unused `rgb2hed` import. Keep `to_gray_inverted` (run_pipeline.py gates
the seed on it), `to_hematoxylin_od` and `to_channel`. Update `chromatin.chromatin_density`'s
docstring, which names `channels.to_hematoxylin`.

## Step 5: final checks

- `compare baseline <B10 label>`: identical, apart from the B7 `evaluate_run` changes and the
  B3/B4 `api.json` changes listed above.
- Every kept augmentation symbol from decision 6 still exists (confirm from `api.json`), and the
  Step 2 augmentation sanity checks still hold.
- pyflakes: only `midog_utils/__init__.py ... imported but unused` lines (eight after B9).
- The B9 import one-liner succeeds.
- The CLI smoke test and the demo-notebook run from Step 2 both pass again, with output under
  `HARNESS_DIR/runs/final/`.

## Final report

1. A table of batches (B1-B10): applied or not; static gate, production harness, augmentation
   harness, `api.json` differences, and `verify_fixes` results. For any batch not applied: what
   diverged and why you stopped.
2. Final `git diff --stat`.
3. The Step 2 trust-check numbers, including the augmentation sanity checks.
4. Out-of-scope files that now reference removed code, from grep over `*.py` and `*.ipynb`
   excluding `midog_utils_full/`. Don't run or fix them. Report two groups separately, since they
   break for different reasons: (a) files importing the deleted `experiment`/`baselines`/`compare`
   modules; (b) files using a removed function, `FSConfig` field, `CHANNELS` entry or
   `seed_selection` parameter in code (not just in prose or comments). Group each by what they use.
5. Stale references in `.md` files (e.g. the walkthrough's `check_no_cap`/`max_detections`
   passages, and `PRODUCTION_PIPELINE_CLEANUP.md`, which still lists `plant_and_recover` for
   removal), listed and not edited.
6. Anything that surprised you.

Leave `HARNESS_DIR` in place.
