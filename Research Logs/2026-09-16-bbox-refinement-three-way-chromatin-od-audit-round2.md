# Round-2 audit of `bbox_refinement_three_way_chromatin_od.ipynb`: all nine fixes landed and the two new cells compute what they say (915 values, 0 divergences); two summary sentences misstate their own tables, and the notebook no longer runs because a concurrent cleanup deleted a function it calls

**Scope.** `production_hematoxylin_only/bbox_refinement_three_way_chromatin_od.ipynb`, re-executed 19:18:23–19:23:58
(30 cells, 17 code). Its eight CSVs are `results/precision_at_k_14roi_prodseed_chromatin_bbox3way_*` (per_roi,
summary, by_domain, delta_per_roi, delta_stats, top30, verification, and the new contested_sensitivity). Also in
scope: its three figures, and the working-tree `midog_utils` modules it imports. Round-1 log:
`Research Logs/2026-09-16-bbox-refinement-three-way-chromatin-od-audit.md`. Audit script:
`bbox_refinement_three_way_chromatin_od_audit_round2.py`; tables: `results/bbox_refinement_three_way_chromatin_od_audit_round2_*.csv`.

**Everything below is re-derived from `databases/MIDOG++.json` and the notebook's top-30 lists.** Template centres and
refusals come from round 1's pixel re-implementation (its persisted `..._tier_b_*` tables). Sign-flip p-values are
recomputed in exact rational arithmetic. `midog_utils` is imported only as the object under test (§ runnability).

## Conflict of interest

The notebook and its artifacts are untracked. The coordinator who requested this round wrote the fixes, and I wrote
round 1, so round 1's definitions are the yardstick for several checks below. Where that matters (the contested
rules, the tie table), I recomputed from primary data rather than comparing to my own round-1 numbers, and I say
which is which. A separate process is editing `midog_utils` right now: dataset.py 19:19:06, seed_selection.py
19:25:17, production.py and find_and_suppress.py 19:29:49, all uncommitted. That matches commit `af8bd0d`, "Add
agent prompt for executing the production pipeline dead-code cleanup". The runnability finding is a snapshot as of
~19:35.

---

## Part 0 — what reproduces

**The published statistics are the statistics the data support.** 915 values compared, 0 divergences
(`..._round2_comparison_counts.csv`).

| check | values | divergences |
|---|---|---|
| Execution gate: `execution_count` 1…17 contiguous, 0 errors, last cell run; figure mtimes (19:23:51/54/57) now match the cells that write them | 17 cells | — |
| Value identity with the 17:40 run: top-30 `cx`, `cy`, `od`, `bucket` and `matched_ann_id` all equal, and top-30 coordinates equal round 1's pixel re-implementation | 8 checks over 1,260 rows each | 0 |
| TP re-matched from the JSON: bucket agreement (42 × 30), TP@K and precision@K per run | 294 | 0 |
| summary, by_domain, delta_per_roi | 297 | 0 |
| `delta_stats` incl. the new `min_attainable_p` and `holm_p`, against exact fractions and exact Holm | 72 | 0 |
| `contested_sensitivity`: my own definition of contested (category 1 and `n_mitotic_votes < n_votes`), my own matcher, exact rational sign-flip, exact Holm within each rule | 135 | 0 |
| Contested rules vs round 1's rules. `contested_as_fp` vs round-1 `contested_mitoses_dropped` is a genuinely different construction (round 1 re-matched with contested GT removed) that coincides because no detection has two annotations within r. `contested_excluded` vs round-1 `contested_neutral` is the same formula recomputed by the same hand, so it confirms arithmetic, not the definition | 36 | 0 |
| TIES table (every printed value), with `same_template` decided from pixel-derived centres | 45 | 0 |
| Refusal attribution: every `refused_draws` string and both printed ROI lists, against round 1's independent RNG walk and Otsu gate | 16 | 0 |
| Working tree: HEAD vs current `seed_selection.tightened_template_box` on all 993 border-filtered unanimous candidates × 2 channels | 1,986 calls | 0 |
| Working tree: current `production.run_production_pipeline` on 245 and 403, all three conditions: top-30 lists and `n_detections` | 12 | 0 |
| Figures: all three renders read; Figure 3's nine `.3g` p labels match `..._round2_figure3_p_labels_expected.csv`; Figures 1–2 unchanged from round 1 | 9 labels | 0 |
| Closing-summary claims (`..._round2_closing_summary_claims.csv`) | 21 | **2 do not hold** |

### The nine round-1 findings

| # | round-1 finding | status | evidence |
|---|---|---|---|
| 1 | Tie paragraph | **Fixed.** The "correcting the reference" framing is gone, and the mechanism claim is replaced by measured statements from the new TIES cell. One number in the rewrite is mis-scoped (F2). | TIES 45/45 |
| 2 | Click caveat omits `gray_bbox` | **Fixed, and the logic is correct and complete.** | 14/14 refusal strings; lists `[013, 245, 300, 403]` and `[013, 245, 300]` reproduce |
| 3 | "pre-D8 default" | **Fixed.** `FSConfig.base_size` = `tm.BASE_SIZE` = 51 | read |
| 4 | `check_nms_radius` vacuous | **Fixed**, and the new statement is still true of the current `production.py` | read |
| 5 | Holm 0.054 | **Fixed.** 0.4375 / 0.026367 / 0.054688, printed 0.44 / 0.026 / 0.055; Figure 3 labels `.3g` | exact |
| 6 | "one test in ten" | **Fixed.** `min_attainable_p` = 2^(1−n) is correct on 9/9 | exact |
| 7 | D9 scope | **Fixed** (cell 0 and caveats) | read |
| 8 | Provenance | **Fixed** for figures and Gate 2. A new provenance problem appeared (F1) | mtimes |
| 9 | Losses list | **Fixed.** Loses on 246/402/459/529/548; wins on 094/233/301 | recomputed |

### The new logic, item by item

- **Refusal attribution.** `draw_joint_click` evaluates all three conditions per draw with no short-circuit and
  records those returning `None`. `template_spec` returns `None` for an Otsu refusal *or* an unreadable 73 px patch
  at the template centre, and production `build_seed` applies the same pair of tests to gray alone. The test
  "production gray would click differently iff some refused draw was not refused by `gray_bbox`" is exactly right:
  the gray-only walk follows the joint walk (same RNG stream, same removals) until the first draw gray accepts.
  `'gray_bbox' not in …` is a safe substring test, since no condition name contains another. `default_51` can never
  appear in a refusal, because `border_filter` guarantees click readability.
- **TIES `same_template`.** "Equal `base_size` and equal `round()`ed centre" is the right identity test: same
  channel, odd size, and `read_padded_patch` rounds the same way. It classifies 201 and 245 as identical for
  `gray − default`, and every other ROI × pair as different (12/14/14), matching the pixels.
  `shared_topk_coordinates_where_templates_differ` is a total over differing ROIs.
  `mean_mitoses_found_by_both` is a mean over **all 14** ROIs. That difference is what F2 is about.
- **Contested definitions.** `contested_excluded` = unanimous TP / (K − contested hits) and `contested_as_fp` =
  unanimous TP / K, with matching unchanged. These are the definitions I would use. `contested_as_fp` gives values
  identical to round 1's stricter version that re-matched without contested GT, as expected when no detection has
  two annotations within the radius. `annotations['unanimous']` equals `n_mitotic_votes == n_votes` for category 1,
  whose only label multisets are (1,1) and (1,2,1)/(2,1,1). The smallest denominator reached is 7, so no division
  hazard.
- **Float-tolerant sign flip.** It equals the exact rational p on **27/27** sensitivity tests and 9/9 main tests.
  Zeros are dropped at |d| ≤ 1e-12 (equal rationals give identical doubles here). The 1e-9 absolute tolerance on the
  observed sum is safe at these denominators in practice. In principle, two distinct achievable sums can differ by
  1/lcm(denominators), which can be about 1e-12 when K − contested hits ranges up to 30, so a tolerance below that
  gap, or `Fraction` arithmetic, would make it provably exact. No effect on any number here.
- **`holm_adjust`.** `maximum.accumulate` of (m − rank)·p over ascending p, clipped at 1, is the standard step-down
  adjustment (clip and running max commute). It matches exact Holm on all 36 tests.

---

## Part 1 — findings

### F1 — bug (environmental, not introduced by the fixes): the notebook no longer runs; cell 4 calls a function deleted during its own execution

Cell 4, line 2: `ds.check_invariants(annotations)`. The working-tree `midog_utils/dataset.py` (mtime 19:19:06,
uncommitted) no longer defines it; `PRODUCTION_PIPELINE_CLEANUP.md` lists it for removal. The kernel imported
`dataset` at 19:18:23 and called it at 19:18:29, 37 s before the file changed, so the recorded run succeeded. Run
now, cell 4 raises `AttributeError: module 'midog_utils.dataset' has no attribute 'check_invariants'`
(`..._round2_runnability.csv`, cells 1 and 4 executed verbatim in a subprocess). The static name scan (`..._round2_runnability.csv`) finds no other missing name. The names cell 1 reads
through `getattr(prod, k)` are invisible to that scan, but the verbatim subprocess run executed cell 1, including its
`CONFIG_DRIFT` assert, and reached cell 4 before failing.

Three more imported modules changed after the run (seed_selection 19:25, production and find_and_suppress 19:29). I
checked that none of those changes moves a number:

- `tightened_template_box` is identical between HEAD and the working tree on 1,986 candidate × channel calls;
- the current `run_production_pipeline` reproduces the top-30 list and `n_detections` for all three conditions on
  245 and 403;
- the edits I read remove unreachable branches for this call path (`max_detections` = 1e9, `score_threshold` unused
  when `deep_floor_z` is set, `method`/`center_tolerance`/`headroom_frac` at their defaults).

**The published numbers stand; the notebook as written cannot reproduce them.** The recorded run is valid. The
notebook's citation of the round-1 audit is unaffected, because that audit re-derived from pixels and the JSON rather
than by executing the notebook. `PRODUCTION_PIPELINE_CLEANUP.md` removes `check_invariants` deliberately, so the
durable fix is on the notebook side.

**Fix:** delete the `ds.check_invariants(annotations)` line, or inline its two assertions (all boxes 50×50;
categories ⊆ {1, 2}). Re-run after the cleanup commit lands, and commit the notebook with the module state it ran
against. Cell 0's "`center_tolerance=0`" names a parameter that no longer exists; the behaviour is unchanged, so this
is wording only.

### F2 — wrong number (minor): cell 28, "What the ties are", K = 10

> *"On ROIs where the templates differ, the top-10 lists share only 1-2 exact candidate coordinates in total, **yet
> find 5.2-6.1 of the same mitoses per ROI on average**."*

The sentence scopes both numbers to ROIs where the templates differ. 5.2–6.1 is `TIES.mean_mitoses_found_by_both`,
a mean over all 14 ROIs, which for `gray − default` includes the two identical-template ROIs (201, 245):

| K = 10 pair | all 14 ROIs (the notebook's column) | ROIs where templates differ |
|---|---|---|
| `gray − default` | 6.14 | **6.75** (12 ROIs) |
| `hem − default` | 5.29 | 5.29 (14) |
| `hem − gray` | 5.21 | 5.21 (14) |

The point the sentence makes survives, and is slightly stronger. **Fix:** compute the column over differing ROIs
only and quote 5.2–6.8, or drop "on ROIs where the templates differ" from the second clause.

### F3 — misleading claim: the cell 28 heading contradicts a bullet further down the same section

Heading: *"`hem_bbox` vs. `default_51` is the one contrast that resolves, **but only under this TP convention**."*
Bullet: *"Excluding them from both numerator and denominator **keeps the result** (Holm p 0.029 at K = 20 and at
K = 30)."* The bullet is right; a reader who skims keeps the heading. Holm ≤ 0.05 at:

| rule | K where Holm ≤ 0.05 |
|---|---|
| `all_mitoses` (as run) | 20 |
| `contested_excluded` | **20 and 30** (0.0286 each) |
| `contested_as_fp` | none |

It resolves under two of the three conventions and fails only when contested mitoses are charged as errors.
**Fix:** "…the one contrast that resolves; it survives excluding contested mitoses, but not counting them as false
positives."

### F4 — misleading claim (minor): cell 28 caveats, "On clicks either gate refuses, the tightened conditions don't exist"

On 013, 245 and 300 only the hematoxylin gate refused the first draw, so `gray_bbox` did exist there. That is the
point of the bullet above it. **Fix:** "…at least one tightened condition doesn't exist, so the three-way comparison
isn't defined there."

---

## Part 2 — verdict per conclusion (rewritten closing summary)

| claim | verdict |
|---|---|
| Every check passed (cap 42/42, K delivered, annulus 42/42, gates 14/14, 28/28, identity) | reproduces |
| Audit citation: 42 top-30 lists and every TP count re-derived from pixels and JSON, 0 divergences | reproduces (values identical to the 17:40 run) |
| TP/FP convention description; identical TP in 126 cells under Hungarian, mitoses-only, r ± 1 px; no detection with both classes within r; nothing within r of a click | reproduces. Round 1 ran maximum-cardinality (Hungarian) matching over the mitotic subset only. It transfers to the mixed pool here because 0 of 1,260 detections are near both classes, so the bipartite graph splits into disjoint mitotic and look-alike parts. It would not transfer at a depth where a detection sits near both |
| Pooled table; ordering; neither refinement beats 51 px | reproduces |
| `hem − default`: 3.6/6.8/4.3; L 5/11/10, W 0/1/1 (013); p 0.0625/0.0029/0.0068; Holm 0.44/0.026/0.055; K=10 floor 0.0625 | reproduces |
| Heading "resolves, but only under this TP convention" | **does not reproduce, in part** (F3): "resolves" holds; "only under this TP convention" is false and understates robustness, since it also resolves with contested mitoses excluded |
| Contested share 22.5–25.5 %, even across conditions; excluded → Holm 0.029 at K=20 and K=30; as FP → p 0.064/0.037, Holm 0.52/0.33; direction negative under all three rules | reproduces |
| `gray − default`: −1.4/−2.9/−1.7; T 10/10/6, L 3/4/5, W 1/0/3; p 0.625/0.125/0.297; floor 0.125 at K=10/20; loses 246/402/459/529/548, wins 094/233/301; 245/201 identical, 403 the only same-size shift | reproduces |
| `hem − gray`: −2.1/−3.9/−2.6; W/L/T; p 0.25/0.090/0.098; Holm 0.75/0.54/0.54 | reproduces |
| Worst ROI 245; hem 0.20→0.10 and 0.23→0.17 | reproduces |
| Size-order observation, not dose-response | reproduces |
| K=10: ties 9–11 whatever the gap (9.9–21.7 px); 1–2 shared coordinates where templates differ | reproduces |
| K=10: "yet find 5.2–6.1 of the same mitoses per ROI" (scoped to differing templates) | **does not reproduce as scoped**: 5.2–6.75 (F2) |
| K=20: 10/2/2 ties at gaps 11.9/21.7/9.9; driver not established | reproduces |
| Click caveats: `default_51` differs on 013/245/300/403; production `gray_bbox` on 013/245/300 | reproduces |
| "On clicks either gate refuses, the tightened conditions don't exist" | reproduces, overstated (F4) |
| `chromatin_od` only; D9 validated on `tm_score`; K ≤ 30 | reproduces |
| The notebook is reproducible from the working tree | **does not reproduce** (F1) |

---

## Part 3 — what was re-run versus read

- **Tier A:** every CSV value and every summary claim, against the JSON, round 1's persisted pixel-derived tables,
  and exact fractions. **Tier B:** the HEAD-vs-working-tree `tightened_template_box` sweep (993 candidates × 2
  channels from pixels), and the current `run_production_pipeline` on 6 runs, about 2.5 min in total.
  **Tier C:** not run, because a full copy would overwrite the notebook's own results CSVs (its output paths are
  hard-coded). Runnability was tested by executing cells 1 and 4's first two statements verbatim, which write nothing.
- **Script:** `bbox_refinement_three_way_chromatin_od_audit_round2.py` runs start to finish (exit 0, 151 s) and
  regenerates all `results/bbox_refinement_three_way_chromatin_od_audit_round2_*.csv`. It imports helpers from the
  round-1 script (which imports no `midog_utils`) and writes nothing under round 1's names.
- **Correction to the request's description:** `delta_stats` did not only gain columns; `exact_p` also went from 4 to
  6 decimals (0.0029 → 0.002930). Harmless.
- **Beyond the fixes, I looked for:** whether the re-run changed any value (no); whether the new cells silently
  reuse round-1 audit numbers instead of computing them (no, all computed in-notebook); whether the contested share
  and rules use the same TP set as the main table (yes, the notebook asserts it and I confirm); whether figure files
  match their execution (yes); and whether the imported modules still support the notebook (no, F1).
- **Prior audits.** Round 1 (mine): all nine findings verified fixed. Of round 1's context items, the contested
  sensitivity was adopted and reproduces. The `tm_score` scope arm and the domain-as-unit check were not adopted;
  they were not requested fixes, and the summary's "nothing here carries over to `tm_score`" stays a hedge rather
  than a claim.
- **Appendix facts:** `FSConfig.tm_method` now defaults to `TM_CCOEFF` in the working tree (uncommitted). The
  appendix statement that it defaults to `TM_CCOEFF_NORMED` no longer holds as of 19:29.

## Part 4 — premises inherited

Unchanged from round 1, Part 4: ROI as the exchangeable unit, 2 ROIs in each of 7 domains; 7.5 µm match radius;
consensus labels; `chromatin_od` as ranker; single click. **Which verdicts would change:** F3 and the "resolves"
verdicts depend on ROI exchangeability; at the domain level nothing survives Holm (round 1, T2-2). The contested-rule
rows are the premise check for consensus labels, and the notebook now reports them itself.
