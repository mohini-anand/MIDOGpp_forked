# Task: re-run the three-way seed-template comparison on 49 ROIs x 3 seeds

## Goal
`production_hematoxylin_only/bbox_refinement_three_way_chromatin_od.ipynb` (committed in `0982979`)
compares three ways of cutting the seed template from one click, ranked by `chromatin_od`, on 14 ROIs
with one click each. Re-run exactly the same analysis on **49 ROIs x 3 seeds**. The question is how
broadly the choice of template refinement changes precision@K. Nothing else about the analysis
changes.

Repo: `/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked`, branch `find-and-suppress-midog`.

## Read first, fully
1. `production_hematoxylin_only/bbox_refinement_three_way_chromatin_od.ipynb`, every cell and its
   saved outputs. This is your template. Reuse its code and structure; don't rewrite from scratch.
2. `Research Logs/2026-09-16-bbox-refinement-three-way-chromatin-od-audit-round3.md`, and skim
   rounds 1 and 2. That notebook went through three independent audit rounds. The logs list the
   mistakes that were caught; don't reintroduce them. Round 3's finding N1 (the `int(round())`
   gate covering only 201.tiff's rounded-centre case) is directly relevant to Gate B below.
3. `DECISIONS.md` D1-D9, `D8_TEMPLATE_ANCHOR.md`, `midog_utils/production.py`,
   `midog_utils/seed_selection.py`, `midog_utils/find_and_suppress.py`.
4. `images/extra_valid/MANIFEST.md` and `images/extra_valid/testing_set/MANIFEST.md`.

## Think it through, don't just execute
This prompt is long and prescriptive, but it was written ahead of the run, not derived from it.
Before implementing each major piece below (the seed-removal rule, the DP-convolution sign-flip
test, the cluster bootstrap, the hard-vs-soft gate split, the `contested_excluded` handling), work
out *why* it's specified that way from the code and data you just read, not just *what* to type.
If your own reasoning turns up an inconsistency, a case the prompt didn't anticipate, or a mistake
in this prompt itself, say so, adjust the plan, and document the deviation and your reasoning in
the final report — the same way the three audit rounds on the committed notebook caught real bugs
by re-deriving things instead of trusting the docstring. This license does not extend to the things
marked non-negotiable elsewhere in this document: the three conditions, the pipeline configuration,
Gates A and C, the "don't tune anything" and "don't commit" rules. Inside those bounds, think for
yourself and correct course as you go rather than mechanically executing a step you've noticed is
wrong.

## What stays exactly the same
- **The three conditions:**
  - `default_51`: 51 px box centred on the click, no Otsu step.
  - `gray_bbox`: `seed_selection.tightened_template_box` on `gray_inverted`.
  - `hem_bbox`: the same on `hematoxylin_od`.
- **Search and ranking:**
  - Every condition runs through `midog_utils.production.run_production_pipeline(rgb, seed, mpp, rank_key="chromatin_od")`, unchanged: `TM_CCOEFF`, the `hematoxylin_od` search channel, deep floor z = -1.5, `MAX_PEAKS = 100`, NMS radius = match radius = 7.5 um, 5 px self-hit removal at the template centre.
  - Don't change any configuration. Don't pass `max_peaks` or other overrides.
- **Scoring:**
  - `evaluate.bucket_detections(detections, gt_eval, evaluate.radius_px(mpp))`, where `gt_eval` is all annotations (both categories) minus that seed's clicked annotation.
  - TP@K = count of `human_correct_label` in the top K. Precision@K = TP@K / K, with **K = 10, 20, 30 only**.
  - Keep the per-run checks, but they are not all equally hard invariants — see "Hard gates vs. reported properties" below for which of these may legitimately vary by ROI: `n_peaks == 100`, at least 30 detections, no NaN `od`, `od`-descending, and no candidate within one match radius of the click.
- **The joint click gate:**
  - Draw from `agreement_pool` then `border_filter(pool, 36, shape)`, retrying on the same RNG stream.
  - Accept an annotation only if `tightened_template_box` succeeds on both `gray_inverted` and `hematoxylin_od`, and all three templates' 73 px patches are readable at their own centres.
  - Record every refused draw with the condition(s) that refused it (`refused_draws`), as the committed notebook does. There is no committed CSV of `refused_draws`/`n_retries` for the 14-ROI run — it only exists in that notebook's `CLICKS` frame and printed output, which you are not allowed to read as data (see Don't). Recompute it fresh for all 147 (ROI, seed) draws in this run, seed 0 included; don't try to reconcile it against anything from the 14-ROI run.
- **Every analysis cell the committed notebook has, except "What the ties are" and "Same size, different centre"** (drop those two):
  - template geometry;
  - per-ROI precision;
  - pooled and per-domain precision with worst ROI;
  - paired deltas for the three pairs `gray_bbox - default_51`, `hem_bbox - default_51`, `hem_bbox - gray_bbox`;
  - wins/losses/ties;
  - the contested-mitosis sensitivity (`all_mitoses` / `contested_excluded` / `contested_as_fp`);
  - the figures.

## What changes
1. **ROIs.** All 14 `images/extra_valid/*.tiff` plus all 35 `images/extra_valid/testing_set/*.tiff`,
   49 in total, 7 per tumour domain.
   - Assert 49 files, 7 per domain, all present in `databases/MIDOG++.json`, no duplicate file names.
   - Keep a `subset` column: `original_14` or `testing_35`.
2. **Seeds: three distinct clicks per ROI.** `SEED_INDICES = (0, 1, 2)`, and each seed walks its own RNG stream `np.random.default_rng([seed_index, image_id])` under the same joint gate.
   - **The three seeds on a ROI must be different annotations.** Draw seeds in order 0, 1, 2. Before seed s walks its stream, remove the annotations already accepted by seeds < s from the pool.
     - Seed 0 has nothing removed, so it is exactly the committed notebook's draw (Gate A depends on this).
     - Removed annotations are not refusals and don't appear in `refused_draws`.
   - Hard-assert that the 3 `seed_ann_id`s per ROI are distinct.
   - Also assert distinct RNG streams with `invariants.check_distinct_seeds`.
   - Compute and report, per ROI, how many pool annotations pass the joint gate. If any ROI has fewer than 3, stop and report. Don't drop ROIs or reuse clicks.
   - **Performance:** load each ROI's pixels and convert its `gray_inverted`/`hematoxylin_od` channels once, then run all 3 seeds x 3 conditions against those same in-memory channels, mirroring how the committed notebook's `run_roi` already computes `channels` once per ROI before drawing. Don't let a per-seed loop trigger a per-seed channel conversion — that's a ~3x, easily-missed cost on top of the already-large 441-call runtime.
3. **Statistics: the ROI is the unit (D5).**
   - *Per ROI, pair and K:* the integer delta summed over the 3 seeds, `sum_s (tp_a - tp_b)`. ROI-level wins/losses/ties use that sum.
   - *Test:* an exact two-sided sign-flip on the 49 ROI-level integer deltas. Enumerating 2^49 is infeasible, so compute the exact null distribution of `sum(+/-|d_i|)` by dynamic-programming convolution over integer offsets.
     - Validate the implementation: run it on the committed notebook's 14-ROI seed-0 deltas (read from `results/precision_at_k_14roi_prodseed_chromatin_bbox3way_delta_per_roi.csv`, not retyped) and assert it reproduces that notebook's `exact_p` values **after rounding both to 6 decimal places** — the committed `delta_stats.csv` stores `exact_p` already rounded to 6dp, so raw float equality against your DP output will fail spuriously even when the implementation is correct.
     - Also report `min_attainable_p` (2^(1 - n_nonzero)).
     - Cross-check against a seeded Monte Carlo sign-flip (at least 200,000 draws).
   - *Multiplicity:* Holm-adjust across the 9 (pair, K) tests, computed on unrounded p.
   - *Effect-size CIs:* a cluster bootstrap resampling whole ROIs, keeping all 3 of a ROI's seeds together (B = 10,000, fixed RNG), giving a percentile 95% CI for the pooled-precision delta of each pair at each K.
   - *D5 bar:* state explicitly, per pair and K, whether both hold: the CI excludes 0, and a majority of ROIs move in the same direction.
   - *`contested_excluded`:* its per-(ROI, seed) values aren't integers, so the DP convolution doesn't apply — use the seeded Monte Carlo sign-flip for that rule only, and say so. Aggregate to the ROI level the same way as the main rule: **sum the per-seed delta over the 3 seeds before running the sign-flip**, not an average. (The sign-flip p-value and CI-excludes-zero conclusion are identical either way since summing is a positive rescaling per ROI; only the printed effect-size magnitude would differ by 3x if you picked the other convention, so pin it to "summed" for consistency with the main rule and with the delta heatmap below.)
4. **Replication split.** The effect was first seen on the original 14 ROIs at seed 0. Report pooled precision, pooled deltas, bootstrap CIs and ROI-level W/L/T separately for `original_14`, `testing_35`, and `all_49`. `testing_35` is the out-of-sample check. Run the Holm-adjusted tests on `all_49` (primary) and on `testing_35`.
5. **Per domain.** 7 ROIs per domain: pooled precision per condition and K, mean delta, ROI-level W/L/T, exact sign-flip p with its `min_attainable_p` (at least 2/128). These are descriptive; no Holm claims per domain.
6. **Per seed.** Pooled precision per condition and K for each seed separately, to show consistency across seeds. Descriptive only, never a test.

## Reproduction gates
Not all of these are true invariants — some are correctness properties that must hold no matter
which ROI you run on, and some are properties the 14-ROI notebook happened to measure, not
properties this run is guaranteed to reproduce. Conflating the two means either missing a real bug
or halting a 60-75 minute run over something that isn't one. Run all of them before any table;
only halt on the ones marked **(halt)**.

- **Gate A (halt):** seed 0 on the original 14 ROIs must reproduce the committed outputs exactly.
  - Against `results/precision_at_k_14roi_prodseed_chromatin_bbox3way_per_roi.csv`: `seed_ann_id`, `base_size`, `tpl_offset_px`, `n_detections`, `tp_at_10/20/30`.
  - Against `..._top30.csv`: `cx`, `cy`, `od`, `bucket`, `matched_ann_id` for all 42 runs.
- **Gate B (report, don't halt):** for every (ROI, seed) where the `gray_bbox` template covers the same pixels as `default_51`'s (base size 51 and the same `int(round())` centre), check whether both runs return an identical top-30 list and identical TP counts. They usually will, but this is not guaranteed: `find_and_suppress` (`midog_utils/find_and_suppress.py`) computes the self-hit distance from the **unrounded float** template centre (`seed_x, seed_y = float(seed_xy[0]), float(seed_xy[1])`), not the rounded pixel. Two conditions whose centres round to the same pixel and therefore produce byte-identical templates can still have their 5.0 px self-hit annulus centred a fraction of a pixel apart, which can flip a peak near that boundary in or out. Round 3 of the audit hand-checked one such case (201.tiff, 0.5 px offset) and found no divergence, but that is one instance, not a proof it can't happen across 147 (ROI, seed) pairs. For every divergence: record the float offset between the two centres and whether the differing detection(s) sit near the 5.0 px self-hit boundary. Only escalate to a real problem if a divergence shows up **without** matching centres rounding to the same pixel, or with the template pixels themselves differing — that would mean the template cut, not the self-hit float precision, is the cause.
- **Gate C (halt):** within each (ROI, seed), all three conditions share `seed_ann_id` and `n_gt_mitotic`. Across seeds within a ROI, the three `seed_ann_id`s are distinct.
- **`max_peaks_binds` / `n_peaks == 100` (report, don't halt):** D9 measured this cap landing on exactly 100 on all 14 original ROIs at `DEEP_FLOOR_Z = -1.5`; it was never characterized on the 35 new ROIs. A run that comes in under 100 peaks is a fact about that ROI's response-map statistics, not necessarily a pipeline defect. Record it; only treat it as a problem if it also breaks the `>= 30 detections` floor.
- **"No candidate within one match radius of the click" (report, don't halt):** `self_hit_radius` (5 px, fixed) is smaller than `match_radius` (`evaluate.radius_px(mpp)`, ~7.5 um -> ~30 px, mpp-scaled). A genuinely distinct nearby mitotic figure that falls between those two radii is a real, legitimate detection near the click, not a self-hit that leaked through — the committed 14-ROI notebook never happened to draw a click with that geometry, but 147 draws over 49 ROIs might. If this fires, report which annotation the nearby detection actually matched (via `bucket`/`matched_ann_id`) before assuming it's a bug.
- **At least 30 detections, no NaN `od`, `od`-descending (halt):** these are correctness properties of the pipeline's OD computation and ranking, not ROI-content-dependent measurements, and precision@30 is not computable without them. A failure here means something is actually broken.

If a gate marked **(halt)** fails, stop and report. Don't change configuration to make it pass.
For the "report, don't halt" checks, include the diagnostic detail above in the verification CSV
and in your final report, but continue the run.

## Current code state
`midog_utils` was leaned out in commit `c066829`:
- `dataset.check_invariants`, `compare.py`, `experiment.py` and `baselines.py` are gone.
- `tighten_box_otsu` has no `method` or `center_tolerance` parameters.
- The committed notebook already runs on this code.

Before running, check `git status --short midog_utils/`. If it shows modifications, another session may be editing the package: wait until it is committed or unchanged for 10 minutes, and say so in your report. Don't modify `midog_utils`.

## Outputs (never overwrite the committed 14-ROI files)
- Notebook: `production_hematoxylin_only/bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb`
- CSVs: `results/precision_at_k_49roi_3seed_chromatin_bbox3way_{per_run,per_roi,summary,by_domain,by_seed,delta_per_roi,delta_stats,bootstrap_ci,contested_sensitivity,top30,verification}.csv`
- Figures: `production_hematoxylin_only/bbox3way_49roi3seed_*.png`, using the committed notebook's colours and plotting style:
  1. pooled precision per condition per K with bootstrap 95% CIs, for `all_49` / `original_14` / `testing_35`;
  2. per-ROI seed-summed delta heatmap (49 rows, 3 pairs x 3 K, shared colour scale) — sum over the 3 seeds, matching the ROI-level statistic used everywhere else, not an average;
  3. per-domain pooled delta for each pair;
  4. ROI-level wins/losses/ties per pair and K, labelled with Holm p (`.3g` format).

## Running
- Use `/Users/mohinianand/anaconda3/bin/python3`, the only interpreter here that imports cv2. Execute with `/Users/mohinianand/anaconda3/bin/jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=10800 <notebook>`, run in the background.
- Expect roughly 60-75 minutes on this CPU-only machine: 441 pipeline calls (49 x 3 x 3), each recomputing the `hematoxylin_od` channel. Convert each ROI's channels once for the gate, not once per seed.
- Smoke-test the analysis cells on synthetic or partial data before the full run, so a bug in a late cell doesn't cost a second hour.

## Writing rules (the committed notebook's audits enforced these)
- Every number or claim in markdown must be printed or tabled by a cell above it.
- No mechanism claims (why a condition wins) that the notebook doesn't measure. Present template size as a confounded observation, not a cause.
- Quote exact p and Holm p, and give `min_attainable_p` wherever ties limit a test.
- State the contested-mitosis sensitivity next to any "resolved" claim.
- Caveats to include:
  - clicks come from the joint gate, not each condition's own production draw; report per subset how many (ROI, seed) pairs `default_51` alone and production `gray_bbox` alone would have clicked differently, from `refused_draws` (recomputed in this run — see "The joint click gate" above, there is no committed 14-ROI CSV of this to check it against);
  - `chromatin_od` ranking only;
  - D9's cap was validated with `tm_score` ranking, not `chromatin_od`.
- Code style:
  - function definitions and calls on one line, never wrapped;
  - docstrings with the opening `"""` alone on its line and content indented one level deeper: a one-line summary, one line per parameter (`name (type): description`), one `Returns type: description` line.

## Don't
- Don't commit.
- Don't add `tm_score`, other rankers, or K > 30.
- Don't tune anything.
- Don't read another notebook's printed output as data; read the committed CSVs.

## Report back
- Gate A/B/C results (for Gate B and the two "report, don't halt" checks: how many divergences, and whether the diagnostic points to the known float/self-hit-radius mechanisms above or to something unexplained), runtime, and per-ROI joint-valid pool sizes (all must be at least 3).
- The `all_49` pooled precision table.
- For each pair and K: mean delta, bootstrap CI, ROI-level W/L/T, exact p, Holm p, whether D5's bar is met, and whether `testing_35` agrees with `original_14`.
- The contested-mitosis sensitivity for any resolved contrast.
- Any point where you deviated from this prompt because your own reasoning (or something you found in the data) said it was wrong, and why.
- Anything surprising, and the list of files written.
