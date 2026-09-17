# Round-3 audit of `bbox_refinement_three_way_chromatin_od.ipynb`: all four round-2 fixes are correct, all 41 closing-summary claims hold, and the notebook runs against the cleaned `midog_utils` (now committed as `c066829`) and reproduces its own outputs exactly. The notebook is accurate and correct

**Scope.** `production_hematoxylin_only/bbox_refinement_three_way_chromatin_od.ipynb`, re-executed 20:13:11–20:16:53
(30 cells, 17 code). Its eight CSVs, `results/precision_at_k_14roi_prodseed_chromatin_bbox3way_*`, are dated
20:16:46–51, and its three figures 20:16:51–53. Prior logs:
`Research Logs/2026-09-16-bbox-refinement-three-way-chromatin-od-audit.md` and `…-audit-round2.md`. Audit script:
`bbox_refinement_three_way_chromatin_od_audit_round3.py`. Tables:
`results/bbox_refinement_three_way_chromatin_od_audit_round3_*.csv` (24).

**Everything below is re-derived from `databases/MIDOG++.json`, the TIFF resolution tags, the notebook's 20:16 CSVs
and printed outputs, round 1's persisted pixel re-implementation tables, and a fresh execution of the notebook.**
Sign-flip p-values and Holm adjustments are recomputed in exact rational arithmetic. `midog_utils` is imported only as
the object under test.

## Conflict of interest

The notebook, its CSVs and its figures are untracked. The coordinator who requested this round wrote the fixes, and
I wrote rounds 1 and 2, whose definitions are the yardstick for several checks. To limit that:
- every table value is recomputed from the JSON or pixels, not compared to my own earlier numbers;
- the round-2 comparison (§ identity) is a regression check only, not evidence of correctness.

During this audit the dead-code cleanup was committed: `c066829` at 20:30:49 (parent `af8bd0d`). The working-tree
`midog_utils` and `production_pipeline` now equal that commit. What this audit cannot cover is the design itself
(one click per ROI, `chromatin_od` only, K ≤ 30); those are disclosed in the notebook's caveats and belong to a
premise review.

---

## Part 0 — what reproduces

**937 values compared, 0 divergences** (`..._round3_comparison_counts.csv`). **41 of 41 closing-summary claims
hold**, and each claim's wording was found verbatim in cell 28 (`..._round3_closing_summary_claims.csv`).

| check | values | divergences |
|---|---|---|
| Execution gate: counts 1…17 contiguous, 0 errors, last cell run | 17 cells | — |
| Kernel start 20:13:11 is after the last module edit (20:03:03). The working tree equals `c066829`, so the run used exactly the committed cleanup | provenance | — |
| **Tier C:** a sibling copy executed with only `OUT_STEM`/`FIG_STEM` redirected to a temp dir. It ran against `c066829` in 337 s with counts 1…17 and 0 errors. The copy was then deleted | 1 run | — |
| Tier C: all 8 CSVs, compared cell by cell as text against the notebook's persisted CSVs | 19,881 cells | 0 |
| Tier C: printed output of all 17 code cells (24,570 chars, wall-clock timings masked), 3 embedded figures and 3 saved PNGs, pixel for pixel | 23 | 0 |
| F1: pre-cleanup `check_invariants` (`af8bd0d`) vs the inlined asserts, on the real table and 7 mutations (bad w, bad h, NaN w, category 3, category 0, look-alikes only, empty) | 8 cases | 0 |
| Top-30 lists vs the 17:40 run and vs round 1's pixel re-implementation | 8 checks over 1,260 rows | 0 |
| TP re-matched from the JSON; `per_roi`, `summary`, `by_domain`, `delta_per_roi` | 465 | 0 |
| `delta_stats` against exact fractions and exact Holm | 72 | 0 |
| `contested_sensitivity`: independent contested definition, exact sign flip, exact Holm within each rule | 135 | 0 |
| TIES: all 45 printed values, the renamed header, and the new denominator (`n_differ`, never 0: 12/14/14) | 47 | 0 |
| Refusal attribution: 14 `refused_draws` strings and both ROI lists | 16 | 0 |
| **Tier B:** production's own seed draw (`run_pipeline.select_annotation` → `build_seed`) on all 14 ROIs, on `gray_inverted` and `hematoxylin_od`, vs round 1's pixel walk; plus the set of ROIs where production gray clicks differently | 29 | 0 |
| Regression 19:23 → 20:16: every value round 2 compared that round 3 compares again | 858 | 0 changed |
| Stale references: 27 `midog_utils` names used in code (including the 7 read through `getattr(prod, k)`), plus all 50 backticked identifiers in prose, against the cleaned tree | 77 | 0 missing |

**Gates.**
- **Execution:** passes.
- **Provenance:** improved since round 2. The modules are now committed, and the run provably used them. The notebook
  and its outputs are still untracked (note N2).
- **Composition:** unchanged from rounds 1–2, and the Tier C run reproduces it: 42 runs = 14 ROIs × 3 conditions;
  1,260 top-30 rows. The 267 verification records are 5 × 42 + 14 + 14 + 28 + 1.
- **Not in the CSVs:** template centres and refused draws. They rest on round 1's persisted pixel tables and on this
  round's production seed draw. Two claims cite documents rather than data: C36 (the reference notebook's own
  reading) and C41 (`DECISIONS.md` D5/D9). Both were checked against the cited text.

**The published statistics are the statistics the code in the working tree computes,** and a fresh execution
produces them byte for byte.

## Part 1 — findings

**None, at any tier.** Two notes, neither of which is a defect:

- **N1 — Gate 3 is narrower than the sentence that reports it (optional).**
  - *The code:* cell 8's Gate 3 tests only ROIs where `gray_bbox`'s offset is exactly 0, which is `245`.
  - *The summary:* cell 28 says `default_51` returns the identical list "wherever its template is the same pixels as
    `gray_bbox`'s". That also covers `201` (0.5 px offset, same rounded pixel).
  - *Why it isn't a defect:* I checked `201` directly. Identical top-30 coordinates and identical TP at K = 10/20/30
    (C06, C28), so the sentence is true as written. The notebook's own gate does not establish it; this audit does.
  - *Optional fix:* widen the gate's condition to equal `round()`ed centres so it tests both ROIs.
- **N2 — commit the notebook with its outputs.**
  - `midog_utils` is committed as `c066829`, and the 20:13 run used exactly that code.
  - The notebook, its 8 CSVs and 3 figures are untracked, so nothing yet pins them to that commit.

## Part 2 — verdicts

### The four round-2 fixes

| round-2 finding | status | evidence |
|---|---|---|
| F1 — cell 4 calls a function the cleanup deleted | **Fixed, correct, complete.** The inlined asserts behave exactly like `af8bd0d`'s `check_invariants` in all 8 cases (logically equivalent: all(w = 50 and h = 50) is the negation of any(w ≠ 50 or h ≠ 50); the category test is identical). Only error-message detail is lost. `ds.BOX_SIZE`/`MITOTIC`/`LOOKALIKE` still exist, unchanged (50/1/2). The comment "removed by the production cleanup" matches `c066829`'s message. Cells 1 and 4 now run through the asserts, and the whole notebook runs (Tier C). | `..._f1_invariant_equivalence.csv`, `..._f1_constants.csv`, `..._runnability_now.csv`, `..._tier_c_run.csv` |
| F1 — cell 0 wording | **Fixed, accurate.** `tighten_box_otsu` takes the component under the patch centre, which is the rounded click pixel (`read_padded_patch` rounds). It refuses when that pixel is background or the component fails the area or solidity limits. | read |
| F2 — "5.2-6.1 of the same mitoses" mis-scoped | **Fixed, correct, complete.** The column now averages over ROIs where templates differ. All nine values reproduce, including the K = 20/30 `gray − default` rows that also changed (11.07 → 11.92, 14.86 → 16.00). The `hem` pairs are unchanged because their `n_differ` was already 14. "5.2-6.8" matches 5.21–6.75, and the TIES intro matches the code. | `..._ties_recomputed.csv`; C35 |
| F3 — heading contradicted its bullet | **Fixed, correct.** Holm ≤ 0.05 for `hem − default` at K = 20 under `all_mitoses`, at K = 20 and 30 under `contested_excluded`, and at no K under `contested_as_fp`. No other pair resolves. "Survives excluding … but not counting them as false positives" is right, slightly understated (it gains a K), which is fine for a heading. | C14, C15 |
| F4 — "the tightened conditions don't exist" | **Fixed, correct.** All 4 refused draws were refused by at least one tightened condition (`hem_bbox` ×3, `gray_bbox+hem_bbox` ×1), and `default_51` never refuses (click readable on all). | C39 |

### The closing summary, claim by claim (all reproduce)

| claims | verdict |
|---|---|
| C01–C06 "Every check passed": cap 42/42, K delivered (min `n_detections` 94), 0 candidates within r of the click (4,111 full-list candidates), 14 clicks, Gate 2 28/28, identical list wherever same pixels | reproduces (C06, see N1) |
| C07 the audit citation: round 1 imports no `midog_utils`, 42/42 pixel lists, 0 divergences; the 20:16 top-30 equals those pixel lists 1,260/1,260 | reproduces; the citation now covers this run through that identity |
| C08–C11 TP convention: radius 29.6–33.1 px; 126/126 cells identical under Hungarian (mitoses), mitoses-only greedy, and r ± 1 px; 0 detections within r of both classes; nearest other annotation ≥ radius + 29.4 px at every click | reproduces |
| C12–C13 pooled table and strict order `default_51 > gray_bbox > hem_bbox` at every K | reproduces |
| C14–C24 `hem − default`: resolves (only pair), survives `contested_excluded`, fails `contested_as_fp`; 3.6/6.8/4.3; losses 5/11/10, wins 0/1/1 (013); p 1/16, 3/1024, 7/1024; Holm 7/16, 27/1024, 7/128; min p 1/16 at K = 10; contested share 22.5–25.5%; Holm 117/4096 ≈ 0.029 at K = 20 and 30; p 0.064/0.037, Holm 0.52/0.33; negative 9/9 | reproduces |
| C25–C29 `gray − default`: −1.4/−2.9/−1.7; T/L/W 10/10/6, 3/4/5, 1/0/3; p 0.625/0.125/0.297; min p 1/8; loses on 246/402/459/529/548, wins on 094/233/301; 245 and 201 identical; 403 the only same-size centring shift (2.00 px), 0 of 30 coordinates shared, TP unchanged | reproduces |
| C30 `hem − gray` numbers, and equal to the reference notebook's printed `chromatin_od` rows (mean, W/L/T, exact p) | reproduces |
| C31–C33 worst ROI 245 everywhere (tied with 529 only for `gray_bbox` at K = 30); 0.20 → 0.10, 0.23 → 0.17; medians 51/40/28 px, offsets 2.50/3.47 px | reproduces |
| C34–C37 TIES prose: 9–11 ties; gaps 9.9–21.7 px; 1–2 shared coordinates; 5.2–6.8 mitoses; the reference notebook's own K = 10 reading ("recover the same mitotic objects", its 5.214 = this notebook's 5.21); K = 20 ties 10/2/2 | reproduces. C36 verifies the citation, not the inference it cites: the reference notebook's reading was itself rewritten after an earlier audit (its cell 29) |
| C38–C39 click caveats: `default_51` alone differs on 013/245/300/403; production gray differs on 013/245/300, confirmed by **running production's own seed draw**; "at least one tightened condition doesn't exist" | reproduces |
| C40–C41 domain steps of 1/(2K) with 2 ROIs per domain; D5 amendment (5 seeds × 14 ROIs, clustered at the ROI), D9 (`tm_score` arm), `SEED_INDEX = 0` | reproduces |

## Part 3 — what was re-run versus read

**Tiers.**
- **Tier A:** every CSV value and every summary claim.
- **Tier B:** production's seed draw on all 14 ROIs and both channels, about 1 min.
- **Tier C:** full execution of a sibling copy (`production_hematoxylin_only/.audit_tmp_bbox3way_round3.ipynb`) with
  `--ExecutePreprocessor.timeout=21600`. It ran 20:34:43 onward for 337 s against `c066829`.
  - Two strings in cell 1 were changed (`OUT_STEM`, `FIG_STEM`) so no existing result was overwritten. What makes the
    redirect safe is the script's assert that every `to_csv`/`savefig` in the notebook routes through those two names.
  - `REF_RAW` (the 18:12 reference CSV that Gate 2 reads) is only read, never written, so the copy's Gate 2 used the
    same input as the notebook's.
  - The copy and the temp dir were deleted.
  - Round 2 did not run Tier C, because the output paths are hard-coded; the redirect removes that obstacle.

**Script.** `bbox_refinement_three_way_chromatin_od_audit_round3.py` runs start to finish (exit 0, 421 s) and
regenerates all 24 `..._round3_*.csv`. `--skip-tier-b` and `--skip-tier-c` exist for development only.
- It imports helpers from rounds 1 and 2; round 1 imports no `midog_utils`.
- The `midog_utils` mtimes were identical at the script's start (20:33:19) and end (20:40:20).

**Round-2 checks not repeated: 57 values in 5 groups** (`..._round2_checks_not_repeated.csv`, one row per group):
- the 9 values of the old all-14-ROI TIES column, renamed on purpose;
- the 36 rule-definition cross-checks, subsumed because every `contested_sensitivity` value is unchanged 19:23 → 20:16;
- the 12-value production spot-check, superseded by Tier C.

**Beyond the fixes I looked for:**
- **HEAD moving mid-audit.** It did, at 20:30:49. The F1 comparison is therefore pinned to `af8bd0d`, not `HEAD`.
- **Stale prose references to removed names.** None: `dataset.check_invariants` appears only in the historical
  comment, `od51` is `DECISIONS.md`'s name for the statistic, and `precision_at_K` stands for the per-K columns.
- **Production's seed draw.** Whether production's own seed draw, not a re-implementation, matches the click caveat.
  It does.
- **The Gate 3 sentence.** Whether it covers a ROI the gate doesn't test (N1). True anyway.
- **Silent value changes.** Whether the fixes changed any number the summary quotes outside the TIES table. None did.
- **Error-message detail.** Whether the inlined asserts lose the pre-cleanup error-message detail. They do, which
  is harmless.

**Prior audits.**
- **Round 2:** all four findings are verified fixed. Its provenance caveat is resolved for the modules (committed,
  and the run used them) and remains for the notebook and its outputs (N2).
- **Round 1:** its nine findings were verified fixed in round 2, and nothing in this run reopens them. The 19:23 →
  20:16 identity and Tier C carry round 1's pixel-level verification forward to this run.

**Appendix facts.**
- **No longer holds:** the appendix says `FSConfig.tm_method` defaults to `TM_CCOEFF_NORMED`. `c066829` commits
  `TM_CCOEFF`.
- **Also invalidated by `c066829`:** `FSConfig` no longer has `channel`, `score_threshold` or `max_detections`.
  `dataset.py` no longer defines `check_invariants`, `points`, `roi_area_mm2`, `check_roi_scale`, `CATEGORY_NAMES` or
  `SCANNER_ALIASES`. `baselines.py`, `compare.py` and `experiment.py` are gone; the appendix's Tier A independence
  rule names `compare.py` as a module to read.
- **Still hold:** `image_annotations`' `category_id=None` default (`dataset.py` line 84) and `MIDOG_RADIUS_UM = 7.5`.

## Part 4 — premises inherited

Unchanged from rounds 1 and 2:
- **The ROI is the exchangeable unit.** Measured: 2 ROIs in each of 7 domains. Falsified if ROIs within a domain
  behave as one unit; at the domain level nothing survives Holm (round 1, T2-2).
- **7.5 µm match radius (D7).** Falsified if TP counts moved materially with the radius. They move by 0 cells at
  ± 1 px (C09).
- **Consensus labels as ground truth.** Falsified if contested mitoses behaved unlike unanimous ones. The notebook's
  own sensitivity cell tests this.
- **`chromatin_od` as the ranker, and one click per ROI.** Both are disclosed caveats, not conclusions.

**Which verdicts would change if a premise failed:** "resolves" and its Holm values (C14, C15, C19, C22) rest on ROI
exchangeability. C21–C23 rest on the label convention. The fix verifications, the runnability result and every
arithmetic reproduction are premise-free.
