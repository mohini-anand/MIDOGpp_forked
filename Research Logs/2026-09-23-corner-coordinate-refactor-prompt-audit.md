# Audit of `CORNER_COORDINATE_REFACTOR_PROMPT.md`: the line citations, the 53x53 warning and the `n_blanked_px` check all hold, but the one landmine site that *constructs* a click box rests on a premise 11.85% of the annotation file violates and no check in the prompt can catch it, the Step 1 width assertion is false for odd sizes and steers an agent back into the bug it warns about, the high-level `(ROI, df) -> frame` API the task calls for is absent, and `build_seed` — the function the refactor makes mix both conventions — is never called by the only authorised regression driver

**Scope.** Target: `CORNER_COORDINATE_REFACTOR_PROMPT.md` (untracked, 358 lines), read as a
hand-off document for another agent — is it correct, complete, and free of assumptions that would
silently corrupt the refactor it specifies. No file in the repository was edited; this audit is
read-only by instruction. Primary sources: `databases/MIDOG++.json` (26,286 annotations, read
directly, not via `dataset.load_annotations`), all 12 files of `midog_utils/`,
`rerun_bbox3way_postD11.py`, `production_pipeline/run_pipeline.py`,
`chromatin_od_ranker_seed_robustness_audit.py`, `.gitignore`, `DECISIONS.md` D8, and
`results/precision_at_k_49roi_3seed_chromatin_bbox3way_postD11_per_run.csv` (441 real runs) for the
empirical checks. **Every number below is re-derived from the JSON, the code, or that CSV — never
from the prompt's own prose.**

---

## Verdict

**Not ready to hand off.** Four blockers, all surgical. The prompt is unusually strong otherwise:
of ~30 line citations sampled, every one was accurate; the argument for holding `nms_by_distance`
and `evaluate.greedy_match` on raw centres is sound; and Verification item 1c is a genuinely sharp
check, not a formality. The defects below are fixes in place, not a rewrite.

---

## B1 — Landmine site #3 rests on a false premise, and 1a cannot catch it

The prompt (landmine item 3, and Step 3) requires `Seed.click_box` to be built via
`box_from_center_closed` from the annotation's `cx`, asserting it "round-trips exactly, since it is
built from the same closed-convention `cx`."

Re-derived from `databases/MIDOG++.json`:

| quantity | value |
| --- | --- |
| annotations | 26,286 |
| `bbox` entries that are Python ints | **none** — every value is a float (range −30.0 … 7233.0) |
| `w = x2−x1`, `h = y2−y1` | 50.0 on **all** 26,286 rows — the prompt is right here |
| `frac(cx) == 0.5` | **2,328** rows |
| `frac(cy) == 0.5` | 2,323 rows |
| non-integer centre in either axis | **3,114 rows = 11.85%** |
| …restricted to mitotic (the seed draw pool) | **343 of 11,937 = 2.9%** |
| `max abs((round(cx)−25 + round(cx)+25)/2 − cx)` | **0.5** |

The pinned formula applies `round(cx)` first, so `center_of_closed_box(box_from_center_closed(cx))`
returns `round(cx)`, losing exactly 0.5 px on 2,328 rows. `box_from_center_closed` is an inverse of
`center_of_closed_box` **only on integer centres**; the prompt calls it "its inverse" unqualified.

The verification gap is the sharper half. Of the four landmine sites, two *read* a box and two
*construct* one. Verification 1a reads the real JSON corners — it exercises the reader and passes
regardless of what the writer does. Writer site #4 (`blank_seed_square`) is covered by 1c. Writer
site #3 has a false premise **and no check**, yet the prompt calls 1a "the single most direct test
of the landmine."

**Fix.** Build `click_box` from the annotation's actual `x1,y1,x2,y2` — the columns Step 2 already
adds — never synthesised from `cx`. Lossless; resolves the unstated ambiguity about what feeds
`offset_px` after the refactor; makes Verification item 4 a real sentinel instead of a check on a
lossy round trip. Consequence to state in the prompt: `build_seed` gains a hard dependency on the
new columns (fine for any frame from `load_annotations`, breaks hand-built frames).

## B2 — Step 1's width assertion is false for odd sizes and steers back into the bug

Step 1 pins the formula correctly (`half = size // 2`, floor division) and then asks the agent to
confirm "the resulting box has width/height exactly `size` (not `size ± 1`)". Measured:

```
floor-div half (the pinned formula): size=51, cx=100 -> (75,125)   x2-x1=50, pixel span=51
true-div  half (the warned bug):     size=51, cx=100 -> (74,126)   x2-x1=52, pixel span=53
```

The 53×53 warning reproduces exactly. But under this repo's own definition of width —
`dataset.py:68` sets `w = x2 - x1`, and `rerun_bbox3way_postD11.py:294` asserts it equals 50 — the
**correct** box has width `size − 1`. An agent following Step 1 literally sees 50 ≠ 51, concludes
the pinned formula is wrong, and "fixes" it toward half-open: precisely the bug the section exists
to prevent.

**Fix.** Assert `x2 − x1 == size − 1` **and** `x2 − x1 + 1 == size`. Name the root cause: "closed
convention" is doing double duty for a continuous interval (GT, `w = 50`) and an inclusive pixel
range (template, 51 px).

## B3 — The high-level `(ROI, df) -> DataFrame` API is absent

The prompt never mentions it, and defines `width`/`height` for emitted detections nowhere — which
collides head-on with B2. `production_pipeline/run_pipeline.py:80` already carries ~90% of it as
`run_pipeline_on_roi(fn, ...)`; it returns a dict rather than a flat frame and emits no corners.

**Fix.** Specify: home in `midog_utils/production.py`; signature
`run_roi(roi_path, annotations, *, rank_key="tm_score", seed_index=0, ann_id=None)` deriving
`rgb`/`mpp` internally via `load_roi`/`roi_mpp` and filtering via `image_annotations`; draws its
seed internally through `build_seed` (which also closes B4); returns `(df, info)` per the existing
convention; columns `rank, x1, y1, x2, y2, width, height, score, angle, flip, scale`, plus `od`
when applicable and the match columns when GT scoring is requested; and
**`width = x2 − x1`, `height = y2 − y1`, pixel span `width + 1`** — forced by `dataset.py:68` and
the driver's assertion, so it must be written down, not left to judgment.

## B4 — `build_seed` has zero regression coverage under the named driver

Step 3 makes `build_seed` mix both conventions in one function — closed for the candidate draw at
`:209`, half-open for the template centre at `:212`/`:220` — and warns that getting it wrong
corrupts `offset_px` by ~0.5 px. But `rerun_bbox3way_postD11.py` **never calls `build_seed`**: it
reimplements the draw via `template_spec`/`pool_gate_specs`/`draw_joint_click` and constructs
`Seed(...)` by hand at `:241`. Verification item 4 therefore cannot be run against the only
authorised driver.

**Fix.** Add a purpose-built `build_seed` harness — cheap, no pipeline run needed: sweep the 49
ROIs × 3 seed indices through `build_seed` before and after, comparing `ann_id`, both centres,
`offset_px` and `n_retries`. It must include at least one ROI whose drawn annotation has a
half-integer centre, or B1 slips through a second time.

---

## Scope contradictions and coverage gaps

- **G1 — Step 0 demands artifacts the edit-scope rule forbids producing.** The driver writes only
  top-30 per run and a fixed column subset; capturing "full detection DataFrames" and "the full
  per-run `info` dict" means *adding output*, not adapting read sites. Verification item 2
  ("identical NMS survivor set … by rank and identity") needs full frames too. Authorise the
  additions or drop them. The stated justification is also wrong: `n_blanked_px` is already a
  driver CSV column (`:251`), so 1c needs no info-dict dump.
- **G2 — `production_pipeline/run_pipeline.py` is excluded but is the natural vehicle.** It is the
  one tracked script exercising `build_seed` + `viz.overlay` + `bucket_detections` end to end, and
  it is the API precedent for B3. It breaks silently at `seed.click_xy`/`seed.template_xy` and at
  the `box = (ty - half, …)` reconstruction in `visualize`.
- **G3 — Step 10's docstring list is incomplete.** Missing: `find_and_suppress.py:53-54` (documents
  the returned columns as `cx, cy`), `blank_seed_square`'s whole param block (its signature
  changes), `border_filter:52`, `seed_selection.py:4`, and `load_annotations`' column list at
  `dataset.py:36`.
- **G4 — Step 8's `cx`/`cy` stash works against the goal.** It saves one re-derive for
  `coverage_fraction` — negligible — while re-injecting `cx`/`cy` onto the returned frame, against
  "corners as source of truth". It also creates an unstated dependency: driver lines `246-247` keep
  working *only because* of the stash. Drop it; adapt those lines under the read-site
  authorisation.
- **G5 — `viz.overlay` gets no smoke test at all** under the current plan.
- **G6 — `midog_utils/` has 12 files, not 10.** The prompt omits `__init__.py` and `channels.py`,
  so as written an agent may not add `geometry` to `__init__.py`'s import list.
- **G7 — no runtime estimate.** `t_pipeline_s` sums to 43 min across the 441 runs in the postD11
  CSV, excluding TIFF loading and channel conversion — realistically ~1–1.5 h per pass, and the
  prompt requires two. It should name the driver's `--rois` and `--select-only` flags.
- **G8 — `template_spec` (`:91`) becomes a fifth landmine site**, since `default_51`'s spec is
  synthetic `(51, cx, cy)` while the other two come from the Otsu gate, so the box convention would
  vary by condition inside the one file the agent may edit. **Low severity**: `compare_selection`
  gates `tpl_cx`/`tpl_cy`/`tpl_sha1` against the reference CSV and fails fast per ROI at `:364`, so
  it costs a wasted run, not silent corruption. Minimal adaptation: keep `template_spec` returning
  `(base_size, tx, ty)`, deriving the centre inside with the right helper per condition.
- **G9 — Step 3's `border_filter` mandate** contradicts Step 2's own backward-compatibility
  rationale for zero numeric benefit (1a proves the values identical). Cosmetic, but inconsistent.
- **G10** — `midog_utils_full/` is a parallel copy of the package and will silently diverge.

## Two prose corrections

- Landmine #1 treats the JSON's `bbox` values as integers; they are floats.
- `chromatin_od_ranker_seed_robustness_audit.py` fails at `d.self_hit_radius` / `d.channel` inside
  `check_config()` (`:1001`), not at `FSConfig()` construction — it is called with no kwargs.

---

## What checked out — do not over-correct these

- **Line citations.** Every one sampled is accurate: `dataset.py:12/58/65`,
  `seed_selection.py:105-140/143-160/47-62/163-177/180-233/209/212/220/134/138-139`,
  `template_match.py:56-72/202/206-238/230-231/232-234`,
  `find_and_suppress.py:25-37/57/82/105/108/116/117-125`, `production.py:35/40/51/59`,
  `chromatin.py:43-47`, `evaluate.py:68-104/78-79/181-184`, `viz.py:18-61/32/38/64-77`,
  `.gitignore:137-138`, and D8's `(x0 + x1 - 1) / 2` at `DECISIONS.md:580-586`.
- **The 53×53 warning is real** and reproduces exactly (table in B2).
- **Verification 1c is a strong check.** All 441 runs in the postD11 CSV carry
  `n_blanked_px ∈ {361, 441, 529, …, 2601}` — every value a perfect square of an odd `base_size`,
  `base_size` itself ranging over 19…51. A `size / 2.0` bug moves this number even when 1a and 1b
  pass, exactly as the prompt claims.
- **The half-open handling is correctly specified and load-bearing**: 134 of 441 runs have a
  half-integer `tpl_cx`, so `(x0 + x1 - 1) / 2` is not a theoretical concern.
- **The regression-driver call is right.** `chromatin_od_ranker_seed_robustness_audit.py` is
  genuinely broken as described, and the "84 values compared, 0 divergences" phrase is literally at
  its `:262-263`.
- **`extract_peaks` emits integer-valued float64** (`template_match.py:202`), so the Step 5
  detection round trip is exact.
- **Setup holds.** `results/precision_at_k_49roi_3seed_chromatin_bbox3way_per_run.csv` — the
  reference the driver reads at `:334` and asserts against at `:364` — **is tracked**, so it
  survives `git worktree add`. The prompt does not say so, but it is true. `images/` needs the
  symlink as described (`.gitignore:137-138` ignores `images/*.tiff`, not the directory); both
  subset directories are populated (14 + 35 TIFFs).
