# 2026-10-06 — Multi-iteration sessions: store the page, blank it in the next search

Decision: `DECISIONS.md` D14. Spec: `ANNOTATION_DATAFRAME_SPEC.md` §0, §2 "The session frame", §3 "Sessions", §5, §8, §9.
Status: applied; code in `6257e49`.

## What was asked

The user wanted the pipeline to take an image and the previous iterations' frame, and return that frame with the current
iteration's rows appended, numbered by an `iteration` column starting at 0. Their design was an input wrapper and an
output wrapper around the existing pipeline.

## What changed in the design, and why

**Blanking.** The user asked for the stored seeds and predictions to be blanked before template matching. That raised
two problems, and both were resolved before any code was written.

- **`run_roi` returns every post-NMS survivor (82–100), not the K shown.** Blanking every stored prediction would erase
  mitoses nobody saw. With validation seed 0 at K = 10, the mitoses ranked 10–99 were:

  | ROI | mitoses at ranks 10–99 |
  |---|---|
  | 300.tiff | 64 |
  | 245.tiff | 8 |
  | 013.tiff | 3 |

  So the output wrapper takes `k` (default 10) and stores only the page.
- **Blanking the RGB image before `run_roi` reaches four places besides the search:**
  - the new seed's Otsu gate;
  - its template, which would be all fill if the human boxed an accepted candidate;
  - the refined-prediction gate;
  - `od`.

  So the blanks go into production's search copy right after the template cut, exactly where D11 blanks the current
  seed. This needs a keyword-only `extra_blanks=None` on `find_and_suppress` and `run_production_pipeline`, both frozen
  files; it is the same pattern as `always_od`.
- **The fill matches either way.** White RGB and D11's fill give the same hematoxylin value (0.0) on 300.tiff and
  245.tiff, so blanking either way would have used the same fill. The choice was about which steps see the blanks.

**What blanking cannot remove.** When a blanked object's peak is gone, it no longer suppresses peaks within the NMS
radius, so one of those can survive beside the blank. The uncovered half of a split figure can also survive. So
predictions within the match radius of a blank centre are flagged, and (user) skipped from the page.

**The guard (user).** The guard raises only when the new box *contains* a stored centre, bounds inclusive. A radius rule
was declined because the human is looking at the image.

**Deferred (user).** The blank's fill value (the channel minimum, which equals glass) and the template refinement are
meant to be replaceable parts.

## Measured

`validation/annotation_frame/session_checks.py`: S1–S6, 0 divergences (spec §3 "Sessions" lists them).

Repeats, logged and not judged: how many of the same seed's raw first 10 lie within the match radius of a stored box.

| | 300.tiff it. 1 | 300.tiff it. 2 | 245.tiff it. 1 | 245.tiff it. 2 |
|---|---|---|---|---|
| without the blanks | 0 | 3 | 1 | 2 |
| with the blanks | 0 | 0 | 0 | 1 (skipped; object 8) |

Each iteration blanked 11 squares (the refined seed and 10 refined predictions): 17,571 and 31,214 px on 300.tiff,
23,692 and 44,007 px on 245.tiff. No detection fell inside a blanked square in S3's production run.

**The §7 harness re-run: 0 divergences.** `find_and_suppress.py` and `production.py` are frozen (spec §0), so the
default path had to be shown unchanged.
- **Setup.**
  - It ran on code byte-identical to `6257e49`.
  - Outputs went to a scratch directory, so the tracked ones were untouched. `common.OUT` was patched after import, and
    the Step 0 inputs are bound at import.
  - The code checksum was the same at launch and at the end.
- **Results:**
  - `refine_checks.py`: 14 checks. `box_checks.py`: 17 checks.
  - Smoke tier: 300.tiff (36,571 values) and 007.tiff (13,005), and 33 named smoke checks.
  - Full tier: 49 ROIs in three shards, **960,825 compared values**, the same count as on 2026-09-24.
  - 33 named full checks, and `extra_checks.py`.
- **Against the committed outputs**, ignoring timings, every regenerated output is identical:
  - the 51 per-ROI result pickles;
  - `refine.json`, `box.json`, `extra.json` and the smoke summary.

  The rest differ only in ways that are not results:
  - The two `named.json` files differ only in V7's `seconds` field (31 values), a timing without the `t_` prefix.
  - The full summary lacks `v2_git`, which was not run. As on 2026-09-24, `--git` cannot pass once
    `annotation_frame.py` changes.
- **In V2's place:** `git diff --name-only 43f1aec` over `midog_utils/`, `production_pipeline/`, the driver and
  `midog_utils_full/` lists only the four intended files.
- **Timing:** the Mac slept for about 25 minutes mid-run (lid closed). That stretched one ROI's wall time and changed
  nothing else.

## Not done

- **Multi-round scoring.** `evaluate` excludes only the current seed and reads the raw list, so `run_frame.py` scores
  iteration 0 only.
- **The precision question.** Does blanking beat list exclusion page for page, and does the glass-like fill promote the
  neighbours of blanks? This is `FUTURE_IDEAS.md` H7; an A/B is about 882 re-searches.
- **`production_pipeline/run_frame_call_tree.html` is stale.** Its line numbers have moved, and it lacks the new calls.
- **The collaborator-repo port.**
