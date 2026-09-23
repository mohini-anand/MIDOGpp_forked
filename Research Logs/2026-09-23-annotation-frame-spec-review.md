# Review of the corner-refactor prompt, its audit, and the annotation-frame spec: the audit's main blockers hold but five of its claims are wrong and it missed nine problems; every number in the spec checks out, but as agreed the spec could not run its own control arm, its validation could not see the bugs it warned about, and it left unsaid whether the functions every recent experiment calls would change — all resolved in b52493b

**Scope.** Three documents, reviewed as they stood at `003b910`:

- `CORNER_COORDINATE_REFACTOR_PROMPT.md`, the prompt;
- its audit, `Research Logs/2026-09-23-corner-coordinate-refactor-prompt-audit.md`;
- `ANNOTATION_DATAFRAME_SPEC.md`, the spec, as agreed.

The questions: is the prompt correct; did the audit catch everything, and is the audit itself
right; does the spec preserve the production pipeline's behaviour exactly; and could the design
carry the kind of experiments this repository has run since 2026-09-12.

**Evidence.**

- `databases/MIDOG++.json`, read directly;
- the committed reference CSVs;
- all of `midog_utils/`, `rerun_bbox3way_postD11.py` and `production_pipeline/run_pipeline.py`;
- the call sites of the recent experiment notebooks;
- one real pipeline run: 300.tiff, seed 0, conditions `default_51` and `hem_bbox`, both rank keys.

**Every number below was re-derived from those sources, not taken from any document's prose.**

**Outcome.** The spec was amended and a new hand-off prompt,
`ANNOTATION_FRAME_IMPLEMENTATION_PROMPT.md`, was written. Both then went through three rounds of
adversarial verification. The old prompt is marked superseded. All of this is committed as
`b52493b`, and the spec's §9 lists every amendment. **Two claims in the first version of this
review were wrong**; they are corrected where they occur, and listed at the end.

---

## 1. Where the audit is wrong

1. **"None of the bbox values are Python ints" is false.** The same error is in the prose
   correction, "they are floats".
   - 95,842 of the 105,144 values are ints and 9,302 are floats.
   - 23,172 rows are all-int and 1,577 are mixed.
   - By image, 87 of 503 have only int boxes, and 416 mix ints and floats.
   - `load_annotations` returns float64 anyway: `cx` and `cy` pass through `/ 2.0`, and `w` and
     `h` share one column across every image. But a per-image frame built from raw JSON comes out
     int64, so any converter must cast.
2. **G8's "low severity" reasoning is wrong.** `compare_selection` fingerprints the template the
   driver computes itself (`template_digest(*spec)`, `rerun_bbox3way_postD11.py:189`). It never
   sees what production cuts from `Seed`.
   - Under the prompt's design, suppose `default_51`'s `template_box` were built with the closed
     convention and read back with `center_of_half_open_box`. Then 113 of the 147 `default_51`
     runs would move the template by one pixel.
   - Neither `compare_selection` nor `n_blanked_px` would see it.
3. **B3's proposed API couldn't reproduce the reference.** Its `run_roi` draws the seed through
   `build_seed`, but the reference was drawn through the joint gate (`draw_joint_click`). The
   agreed spec's §7 had caught this.
4. **B4 asks for a check the reference draw cannot supply**: a drawn click with a half-integer
   centre.
   - None of the 147 reference clicks has one.
   - Only 7 of the 2,075 annotations in the 49 ROIs' agreement pools do.
   - 295 of the 343 half-integer mitotic annotations are in the contested (3, 2) tier, which is
     drawn only when an ROI has no unanimous annotation.
   - The test has to pin an `ann_id`.
5. **G2's "breaks silently" is wrong.** Renaming `Seed.click_xy` and `Seed.template_xy` raises
   AttributeError, and unpacking a box as a centre raises ValueError. Both are loud.

**What the audit got right:** B1, the lossy `click_box`; B4, `build_seed` never being exercised by
the driver; G1, the Step 0 artifacts the edit scope forbids; and G4, the `cx`/`cy` stash. B2 is
right for odd sizes only (see §2, item 2).

## 2. What the audit missed in the prompt

1. **The `_patch_readable` instruction is mislabelled.** The prompt calls its input a closed-box
   read of the raw candidate. Its only caller (`seed_selection.py:212`) passes the half-open
   template centre, which contradicts the prompt's own `build_seed` instruction.
2. **The helper's width assertion depends on parity, and so does the audit's fix.**
   `box_from_center_closed` gives the same box for size 50 and size 51: `x2 − x1 = 50`, 51 pixels
   wide. So the prompt's "width == size" holds only for even sizes, and the audit's fix,
   `x2 − x1 == size − 1`, only for odd ones. The helper should reject even sizes.
3. **The field order is never pinned.** `template_box` would have been `(y0, y1, x0, x1)` and
   `click_box` `(x1, y1, x2, y2)`, on the same dataclass.
4. **Verification 5, the augmented-bank check, cannot be run as written.** It has no baseline in
   Step 0, and no driver: `run_production_pipeline` hard-codes one template. Nor does it name a
   config. The repo's rule names the 8-template bank on 1024 px crops.
5. **Verification 3 requires every metric to match exactly** but doesn't exclude the timing column
   `t_pipeline_s`, so it could never pass.
6. **The driver's output paths are tracked files.** Baseline and refactored runs would overwrite
   each other, and overwrite the reference.
7. **Checks 1a and 1b can only fail if the helper itself is wrong.** Only 1c, 2, 3 and 4 can tell
   a correct refactor from a broken one, and 2 and 3 need the full baselines that G1 notes are
   missing.
8. **Steps 4, 7 and 8 break research that runs today.** They change `blank_seed_square`,
   `score_detections` and `bucket_detections`. Two current callers use the old signatures: the
   NMS-radius ablation prompt committed the same day, and the post-D11
   `fill_geometry_precision_at_10_49roi_6seed.ipynb`.
9. **The setup step for images is broken.** *(Found after the first version of this review.)*
   - `images/.gitkeep` and the two `MANIFEST.md` files are tracked, so `git worktree add` creates
     `images/` as a real directory.
   - `ln -s <orig>/images ./images` then creates `images/images`, and every ROI read fails. This
     was verified in a scratch directory.
   - The audit marked this step correct. The new prompt symlinks the TIFFs one by one.

A smaller point: the detection boxes are sized from `cfg.base_size` and ignore each detection's
scale.

## 3. The spec as agreed (`003b910`)

### What checks out

| claim | re-derived |
| --- | --- |
| half-integer `x1` values | 2,328 (2,323 for `y1`); fractional parts are only `.0` and `.5` |
| boxes with a negative top-left | 169, min −30.0. Also: 105 have `x2 > W` or `y2 > H`, and two centres lie outside the image (305.tiff ann 15363 at `cx = −2`; 309.tiff ann 15728 at `cy = −5`) |
| `w = h = 50` | all 26,286 rows |
| `cx − w/2 == x1`, and `x_tl + w/2 == cx` | bit-exact on all 26,286 rows |
| runs with a half-integer template centre | 134 in `tx`; 207 in `tx` or `ty`. On 103, round-half-even and round-half-up pick different pixels |
| `n_blanked_px == base_size²` | 441/441, `base_size` 19–51 |
| `n_peaks == 100` and `max_peaks_binding` | 441/441 |
| `n_detections` | 82–100, median 99 |
| the current code reproduces the committed top 30 | 300.tiff, seed 0, two conditions: `rank`, `cx`, `cy` and `od` exact; `score` only when compared as float32 |
| survivors don't depend on the rank key | the same 100 survivors under both keys, on the same runs |

The design's central move is sound: the half-open Otsu bbox never reaches the frame. It stays
inside `tightened_template_box`, which removes most of the prompt's rounding traps.

### What was wrong or missing, and how b52493b resolved it

1. **Layering was unstated.**
   - The experiments since 09-12 all call the lower functions directly with `cx`/`cy`: local
     copies of `find_and_suppress` built from the `tm` primitives, `score_detections` and
     `bucket_detections`; `run_production_pipeline` with a duck-typed seed as a parity check;
     and the hybrid gate calling `tightened_template_box` and `border_filter`.
   - `fill_geometry_precision_at_10_49roi_6seed.ipynb` also imports `rerun_bbox3way_postD11` as a
     module.
   - **Resolved (§0):** every existing function is frozen. The frame lives in a new
     `midog_utils/annotation_frame.py` that converts at its boundary. The only change to existing
     code is a keyword-only `always_od=False` on `run_production_pipeline`.
2. **The `default_51` arm could not be run.** `refinement_channel` had no "none" value, and
   `default_51` is 147 of the reference's 441 runs and the control arm of nearly every recent
   experiment. **Resolved:** `refinement_channel=None` (§4).
3. **Validation checked only `tp_at_k`**, which cannot see most of the bugs the spec warns about.
   **Resolved (§7):** checks V1–V10 against references R1–R4 and a baseline captured before any
   edit. Comparisons are exact, `score` is compared as float32 (§3, item 8), NaN counts as equal
   to NaN, and a string compares equal after `fillna("")`.
4. **The rounding mode was unpinned.**
   - Round-half-up would change the template on 103 of 441 runs, and D8 relies on ties going to
     even.
   - "Rounding happens exactly once" was also false: six frozen sites round a coordinate.
   - **Resolved (§1):** the rounding stays in the frozen functions, round-half-to-even is pinned,
     and the one new rounding site is named.
5. **`object_id` mixed two id spaces**, prediction ranks and MIDOG `ann_id`s. On a gate fallback,
   the refined row's key would collide with the prediction whose rank equals the seed's id. This
   is latent in the reference (the smallest seed `ann_id` is 130), but near-certain in production.
   **Resolved (§2):** `class` joins the entity and the row key.
6. **The edge and readability paths were undefined.** The fallback covered only a gate refusal.
   **Resolved (§3–§4):**
   - a click within 36 px of an edge raises, matching every existing draw;
   - a refined template whose centre is too close to the edge falls back, with
     `refusal_reason = "template_unreadable"`;
   - the user decided on 2026-09-23 that edge clicks are not supported for now.
7. **Tie-breaking was unpinned.** The reference top 30s contain 3 exact `od` ties. **Resolved
   (§2):** predictions keep production's returned order and are never re-sorted.
8. **The float read-back advice gave the wrong reason.** `cx` reads back identically under the
   default parser. It is `od` that doesn't: 9,052 of 13,230 values change in the last bits.
   `score` is a separate problem: it was written as float32 text, so compared as float64, 0 of
   13,230 rows match. **Resolved (§2).**
9. **The inputs and outputs of `evaluate` were undefined**, and "match against the human seed's
   geometry" was vacuous: `greedy_match` never sees the seed. **Resolved (§6).**
10. **Smaller items.**
    - "All three stage values differ" is not guaranteed.
    - "Exactly 2 + n_detections rows" conflicted with copying input rows through.
    - "`build_seed` takes a box" contradicted "`run_roi` never draws".
    - The header called the prompt "centre-based" and "rejected".

    **All resolved.**

## 4. What adversarial verification found in the first amendment

There were three rounds, with evidence required for every finding. Four blocking gaps were found
in the first amended draft, and all are fixed in `b52493b`:

- **The default path went untested.** `run_roi` passed `always_od=True` for both rank keys, so
  nothing exercised production's default `chromatin_od` path after the edit. It now passes
  `always_od=(rank_key != "chromatin_od")`.
- **`evaluate` would have excluded nothing.** It filtered ground truth on an argument that
  defaults to `None`, and `ann_id != None` is true on every row. It now excludes by the seed row's
  own id, and asserts that id is a mitotic `ann_id` at the seed's centre.
- **Step 0 would have halted immediately.** The exact-comparison rules had no NaN provision:
  `near_click_detail` is NaN on 437 of R1's rows and `recall_contested` is NaN on 18 runs.
- **One fallback branch was never exercised.** Nothing reached the `template_unreadable` branch.
  A named click now does: 300.tiff at (6410.0, 250.0), where both gates accept but both template
  centres are unreadable.

Rounds 2 and 3 found no blocking issues, and the verifier's final verdict is "ready for hand-off".

## 5. Corrections to the first version of this review

1. **It said no post-D11 `tm_score` reference exists. That was wrong.**
   - The `default_51` arm of `results/precision_at_k_49roi_6seed_fillgeom_postD11_{per_run,top30,clicks}.csv`
     is one: `a9293d0`, committed after the last `midog_utils` change, `9e22eb9`.
   - It covers 49 ROIs × 6 seeds under both rank keys, with `od` in the `tm_score` top 30. Its
     seeds 0–2 are exactly R1's 147 clicks.
   - `results/precision_at_k_49roi_3seed_shapemass_postD11_*` holds an identical copy: 4,410 of
     4,410 top-30 rows match.
   - It is now the spec's reference R4.
2. **It did not flag the images symlink step.** See §2, item 9.
