# Round-3 audit of the rebuilt `production_seed_precision_at_k_chromatin_hem_bbox.ipynb`: all five claimed fixes are genuinely applied and every number re-derives from pixels with zero divergences, but the closing summary's account of the K = 10 ties rests on a comparison the run itself declares invalid and misses the mechanism its own pools can measure

**Scope.** Target: `production_hematoxylin_only/production_seed_precision_at_k_chromatin_hem_bbox.ipynb`
(untracked, rebuilt and re-executed in place 2026-09-16 16:10:58, **after** both prior audit logs —
round 1 at 12:44 and round 2 at 14:28 — so this is a fresh, independent audit of the current
content, not a re-check of their findings). Artifacts:
`results/precision_at_k_14roi_prodseed_chromatin_hembbox_{raw,per_roi,by_domain,delta_per_roi,delta_by_domain,stat_context,verification}.csv`
(all seven written 16:09:50–16:10:49, all untracked); `databases/MIDOG++.json`; the 14
`images/extra_valid/*.tiff` pixels. Figures: all three embedded PNGs decoded and read. Decisions
read through their amendments (`grep -n '^## D\|^### Amendment' DECISIONS.md`): D1, D3, D5
(**+ its 2026-09-08 `od_contrast` amendment — the only amendment in the file**), D7, D8 (which
redirects to `D8_TEMPLATE_ANCHOR.md`), D9 (2026-09-12, no amendment, still operative). Prior audits
`Research Logs/2026-09-16-hembbox-precision-at-k-audit.md` and `-round2.md` read **only after**
Steps 1–4 were complete, and reconciled in Part 1 §R. Audit script:
`hembbox_precision_at_k_audit_round3.py`; tables `results/hembbox_precision_at_k_audit_round3_*.csv`
(21 of them). **Everything below is re-derived from those artifacts and, for all 14 ROIs and both
tightening conditions, from the ROI pixels themselves — never from the notebook's printed output,
never from its output CSVs where a raw one exists, and never from either prior round's
conclusions.**

---

## Conflict of interest

`git log -1 --format='%an %ar' -- production_hematoxylin_only/production_seed_precision_at_k_chromatin_hem_bbox.ipynb`
returns **nothing**: the notebook has never been committed. `git status --porcelain` lists the
notebook, its three PNGs, its seven CSVs, both prior audit logs, both prior audit scripts and their
~40 tables as `??`. The last commits before them (`ad268f0`, `74fbea0`, `bf532d1`) are all
`mohini-anand`'s from the same morning, and `74fbea0` (10:14) touched
`midog_utils/find_and_suppress.py`, which the notebook imports. **Signal seen: uncommitted target,
plus same-author commits hours earlier.** So the rebuild, the two prior audits and this audit are
all plausibly the same hand within one day, and only fresh context separates them.

What limits it: every table here is produced by `hembbox_precision_at_k_audit_round3.py` from
`..._hembbox_raw.csv`, `..._verification.csv`, `databases/MIDOG++.json` and the pixels — never from
the notebook's output cells, and never from its four derived CSVs except *as the thing being
checked*. The Tier B section re-implements, in plain numpy/pandas/scipy with only
cv2/skimage/tifffile/sklearn primitives: `rgb2hed` and `255 − gray` → the min-max/Otsu binary
threshold → 8-connected labelling → the area/`max_area_frac`/solidity gate → D8's
`(x0 + x1 − 1) / 2` anchor → `_odd(max(h, w))` → the RNG draw-and-retry walk → the 73 px patch and
`base_size` cut → replicate padding and `TM_CCOEFF` → the strided median/MAD deep floor →
**`scipy.ndimage.maximum_filter`** peak extraction (a different code path from the repo's
`cv2.dilate`) with the `lexsort` tie-break → greedy distance NMS → the 5 px self-hit filter →
`od51`/`od_ctx`/`od_contrast` → the stable descending rank → the greedy one-to-one matcher →
precision@K → the exact sign-flip enumeration. It calls **no** `midog_utils` function anywhere:
not `compare.evaluate_arms`, not `evaluate.bucket_detections`, not `nms.nms_by_distance`, not
`seed_selection.tightened_template_box`, not `template_match.extract_peaks`, not
`dataset.load_annotations`, not `chromatin.chromatin_density`. Those modules were **read** (all ten,
in full) for the implementation review in Part 1 §I, and read-only was respected everywhere: both
prior logs, both prior scripts, their tables, `DECISIONS.md`, `midog_utils/` and all seven of the
notebook's own outputs are untouched.

What it cannot cover: the project premises in Part 4 — ROI as the exchangeable unit, 7.5 µm as both
NMS and match radius, precision@K on one click per ROI as the read-out, reading burden and
worst-seed behaviour as the product metrics. Those go to `premise-reviewer` and to a reader who did
not make them, not to me.

---

## Part 0 — what reproduces

**The engineering is sound, and more than sound: this is the cleanest run this notebook family has
produced.** 4,982 Tier A values and 812 Tier B values re-derived, **0 divergences in both**. The
published statistics are the statistics the code in the working tree computes, and — on all 14 ROIs
and both conditions re-run from pixels with an independent implementation — the statistics the
pixels support. Every finding in Part 1 is about inference, attribution and prose, not about the
run.

| gate / check | result |
|---|---|
| **Execution coherence** | 29 cells, 17 code cells, `execution_count` **1…17 contiguous from 1**, 0 `output_type == 'error'`, 0 unrun code cells, last code cell run (ec 17), 3 embedded PNGs. Genuinely one kernel session. **Not vendored** — imports `midog_utils`, reads `../databases/MIDOG++.json` and `../images/extra_valid`, names the ROI files; the exemption does not apply. |
| **Provenance — modules** | `git status --porcelain midog_utils/` is **empty**. All ten imported modules committed and clean; newest is `find_and_suppress.py` (`74fbea0`, mtime 2026-09-16 10:14:46). Every artifact (16:09:50–16:10:49) is newer than every module: `artifacts_newer_than_every_module = True`. No `.partial` anywhere. **No module drift.** |
| **Provenance — notebook and artifacts** | Notebook 16:10:58; its seven CSVs 16:09:50–16:10:49 — consistent with one run, then `stat_context` written last, then the notebook saved. All eight files **untracked**, so Tier A alone would be *consistency with an uncommitted artifact*. **Tier B upgrades it to independent reproduction on 14/14 ROIs × 2/2 conditions.** |
| **Composition gate** | 24 checks, **24 pass**. 336 raw rows = 14 ROIs × 2 conditions × 3 arms × 4 budgets; 84 / 168 / 42 / 84 / 12 rows in Table A / Table B / DELTA_A / DELTA_B / STAT_CONTEXT, each exactly what the design implies; 14 distinct ROIs matching the disk and the annotation DB; 7 domains × **exactly 2 ROIs each**; 0 duplicate `(file, condition, arm, budget)` keys; 0 duplicate rows; 0 NaN in any analysis column; both inner merges lose nothing (42 = 14 × 3, 84 = 7 × 3 × 4); 242 verification records, all `passed`. |
| **Tier A — Table A** | 1,261 values — 84 rows × 14 numeric columns, each row's domain label, and the row order — re-derived from `..._raw.csv`: **0 divergences.** |
| **Tier A — Table B** | 1,344 values over 168 rows re-derived, including the markdown's algebraic claim that `precision_pooled` equals the simple mean of the domain's two ROIs (it does, on 168/168, to 1e-9) and the `worst_roi_tied` flags (12/168, and the arbitrary `worst_roi_file` pick is inside the tied set on all 12): **0 divergences.** |
| **Tier A — DELTA_A / DELTA_B** | 673 and 252 values (42 × 16 plus row order; 84 × 3): **0 divergences.** |
| **Tier A — STAT_CONTEXT** | All 12 × 9 values — mean delta, W/L/T, `n_nonzero`, the exact sign-flip p, ρ, its asymptotic p and its 2,000-permutation p — rebuilt with an independent numpy enumeration of all 2^k sign patterns and an independently seeded permutation loop: **0 divergences.** |
| **Tier A — RAW internal consistency** | 1,344 values — `precision_at_budget == tp/delivered`, `recall_at_budget == tp/n_gt`, `budget_delivered == min(K, n_detections)`, `tp ≤ delivered` on all 336 rows: **0 divergences.** |
| **Tier A — figures (step 1a)** | **639** plotted quantities re-derived: Figure 1's 168 heatmap cells (rendered 2 dp label, sign, and in-frame test against the *exact integer* delta), Figure 2's 84 bar heights and its shared y-limit, Figure 3's 12 stacks (each summing to 14): **0 mismatches.** |
| **Tier A — figures (step 2)** | All three PNGs decoded and read. Panel order, row order (domain then file), legend→colour mapping, stack order and `ylim(0, 15)` all correct. Two presentation notes in Part 1 (F3, F4). |
| **Tier B — all 14 ROIs × 2 conditions from pixels** | **812 values** (28 rows × `base_size`, `tpl_cx`, `tpl_cy`, `n_detections`, `n_gt_mitotic`; plus 336 `tp_at_budget` **and 336 `lookalike_at_budget`**): **0 divergences.** Every template size, every half-pixel anchor coordinate, every post-NMS pool size, every true-positive count and every look-alike count at every budget. Pre-NMS peak counts land on exactly 100 on 28/28. |
| **Tier B — the seed draw** | My own gate and my own RNG walk reproduce all **14** `seed_ann_id` values and the same `n_retries` pattern (1,0,0,0,1,0,1,0,0,1,0,0,0,0 — 18 candidates evaluated in total). Gate 1's headline numbers reproduce independently: joint == hem-only solo draw on **14/14**, == gray-only solo draw on **11/14**, with the three divergent ROIs and their gray-solo ann_ids exactly as printed (013→254, 245→6274, 300→14581). |
| **Tier B — the 403.tiff geometry** | One uncapped extraction per condition. `hem_bbox`: 53,746 pre-NMS peaks, exactly **one** annulus survivor, 33.106 px from the nearest kept point — **0.046 px outside** the 33.0599 px NMS radius — and 32.898 px from the D8 anchor. `gray_bbox`: 34,638 peaks, **zero** annulus survivors. The 0.046 px margin and the ~54k pool the notebook quotes reproduce exactly; the rank does not (F2). |
| **No category-filter leakage** | `ds.image_annotations(anns, fn)` is called **without** `category_id`, which is the repo's known trap — and here it is **correct**: `evaluate.bucket_detections` splits the matched GT by category and `tp_at_budget` counts only `HUMAN_CORRECT_LABEL`, so look-alikes occupy a GT slot under one-to-one matching but never score as true positives. My independent matcher, written to the same spec, reproduces all 336 counts. `n_gt_mitotic` is computed with an explicit `== ds.MITOTIC` filter, and the seed pool with another. |
| **All five claimed rebuild fixes** | Verified in the current source, not taken on faith. (1) joint gate: `_check` requires `tightened_template_box` under **both** channels *and* a readable 73 px patch at **both** resulting centres — and both patch reads are from `hem`, which is correct, since `hematoxylin_od` is the search channel for both conditions. (2) `MAX_PEAKS = 100`, passed positionally to `extract_peaks`, binding on 28/28. (3) the assert is inverted — `assert n_peaks == MAX_PEAKS`; no active `assert n_peaks < MAX_PEAKS` survives (the old form appears only inside a comment). (4) `caps=()`, with the reason in-comment; `check_no_cap` is therefore vacuous and the notebook says so. (5) `od_contrast` is a genuine third `Arm` with `rank_key='od_contrast'`, scored on the same pool. |
| **Config drift** | Ten settings diffed against `FSConfig` defaults and `DECISIONS.md`. Channel = `hematoxylin_od` (D3), method = `TM_CCOEFF` (D1), NMS radius = `ev.MIDOG_RADIUS_UM` = match radius (D7), `max_peaks = 100` (D9), `self_hit_radius = 5.0`, `peak_min_distance = 7`, `DEEP_FLOOR_Z = −1.5` (D9). **No stale override.** The 5.0 µm NMS radius that propagated through five earlier notebooks is absent. One scope gap only: K = 50 sits outside D9's tested range (Part 1 F7). |
| **Saturation (D4 / mode 3)** | `coverage_frac` = **0.0085** — 0.85 % of the ROI lies within a match radius of a detection. The "recall ≈ 1.0 because the pool tiles the ROI" failure mode is nowhere near live here. |

---

## Part 1 — findings

### Tier 1

#### F1 — The closing summary's account of the K = 10 ties has three separable clauses: one is uncheckable from this run, one does not reproduce, and one is unmeasured

Cell 27 says (clause labels mine):

> **`chromatin_od` at K=10 is close to fully tied (0 wins / 3 losses / 11 ties), and this is a
> property of the production cap, [a] *not evidence the two conditions rank alike*.** At
> `MAX_PEAKS = 100`, `chromatin_od` and `od_contrast` re-rank the top 100 raw-score peaks, [b] *not
> the ~18,000+ the uncapped notebook this is built on measured — far less room for either axis to
> move a detection in or out of a K=10 list*. … [c] *this is a real, structural consequence of
> running these axes as re-rankers under the D9 cap, worth stating as a finding*.

The tie counts are correct (Tier A, 0 divergences). The three clauses are not one claim and do not
share a verdict.

**[b] the comparative clause — "far less room … than the ~18,000+ candidate pool did" — is `cannot
check` from this run, and the reason is F5.** The re-ranked pool here is 94–99 candidates. Round 2's
Part 0 records the uncapped predecessor's post-NMS `n_detections` at **16,082–18,813**, so the pool
handed to the re-ranker did collapse by a factor of roughly 175 — the clause is prima facie strong.
But the only run that could settle it is the one cells 0 and 18 declare non-comparable, and I
declined the 28 uncapped passes it would take to rebuild it. What I can establish is that **ample
absolute room remains at the cap**, so the comparative premise, even if true, does not deliver the
conclusion: re-ranking the *same* capped 100-peak pool by `chromatin_od` instead of `tm_score`
replaces a **median 6 of the top 10**.

| within one condition, top-10 overlap with `tm_score`'s top-10 | mean | median | min | max | n |
|---|---|---|---|---|---|
| `chromatin_od` | 3.86 / 10 | 4 | 1 | 7 | 28 (ROI, condition) |
| `od_contrast` | 3.93 / 10 | 4 | 1 | 7 | 28 |

On one (ROI, condition) only **1** of the 10 positions survives the re-rank. Whatever the cap took
away, the axes still move six of ten list positions with what is left.

**[a] "not evidence the two conditions rank alike" does not reproduce — at the level that produces
the ties they do rank alike, and more alike than the template-dependent arm does.** This is the
Tier 1 core, and it is measurable from this run alone.

The two conditions' post-NMS pools share almost no candidate *points*: 66 of 1,364 coincide to the
pixel (mean Jaccard 0.030; exactly 0 on 10 of 14 ROIs), and the two top-10 lists under
`chromatin_od` share on average **0.14** of 10 points. Two entirely different ten-candidate lists
nonetheless return the same true-positive count — because they claim the same *objects*:

| K = 10, hem vs gray | mitoses hit by hem | by gray | **hit by both** | ROIs where `tp` ties |
|---|---|---|---|---|
| `chromatin_od` | 6.36 | 6.57 | **5.21** | **11 / 14** |
| `od_contrast` | 6.29 | 6.43 | **5.29** | 9 / 14 |
| `tm_score` | 5.21 | 5.36 | 3.71 | 3 / 14 |

`od51` and `od_ctx` are functions of the hematoxylin image **at a location** and carry no dependence
on the template; `tm_score` is the template's own correlation. So when the template changes, the two
pools land on different *pixels* but — because they cover largely the same objects (mean 69 % of
hem's points pair one-to-one with a gray point inside the match radius, range 19–93 %) — the
chromatin axes select the same *mitoses*, and `tm_score` does not. The tie ratio tracks that
exactly: 11/14 against 3/14.

**[c] "a real, structural consequence of running these axes as re-rankers under the D9 cap, worth
stating as a finding" reproduces, overstated.** The version stated is measured nowhere in the
notebook and, per [b], cannot be measured from this run. The version the notebook's own pools do
support is stronger and is absent from the summary: **at K = 10 the chromatin axes recover
substantially the same mitoses regardless of which of the two templates generated the candidates** —
a direct argument that the tightening channel barely matters for the chromatin arms, made from this
run and nothing else.

*Caveats on my own numbers.* The 69 % pairing figure uses greedy one-to-one matching at the ROI's
own match radius — the same radius the evaluation uses; at a tighter tolerance it falls. The
within-condition re-rank statistic is over 28 (ROI, condition) units, not 14 independent ones. And
[b] is genuinely open: if someone runs the uncapped sweep the comparative clause may well hold — it
would still not rescue [a].

#### F2 — Cell 8 states it "re-derives" the 403.tiff mechanism "rather than asserting it"; it derives one inequality, and the two quantitative claims it makes are inherited — one of which does not reproduce

Markdown cell 7 says the code cell below "re-derives why from the stored peak-extraction diagnostics
rather than asserting it." What cell 8 computes is `self_score >= cutoff_score` on two rows. The two
numbers it asserts are not computed anywhere in the notebook:

| claim (cell 8 comment) | this audit's uncapped re-run, this notebook's own configuration |
|---|---|
| "a survivor **0.046 px** outside the NMS radius" | **reproduces exactly** — 33.106 px from the nearest kept point vs a 33.0599 px NMS radius, on `hem_bbox` |
| "traced to the D8 anchor landing on a half-integer pixel" | **reproduces** — 403 `hem_bbox`'s anchor is (2225.0, **2974.5**); the survivor is 32.898 px from it, inside the match radius, while being 33.106 px from the kept self-peak |
| "that survivor's own rank was **~3282nd** of the **~17,600–54,000** pre-NMS candidates" | **does not reproduce.** The survivor is **4445th** (rank 4444, 0-indexed) of **53,746** peaks, against the quoted 3282nd. |
| the conclusion — the survivor is excluded by the cap before NMS or self-hit removal run | **reproduces, with margin**: rank 4444 ≫ 100 either way |

**The discrepancy localizes entirely to the rank's denominator.** Round 2's own Tier B trace reports
53,746 pre-NMS peaks and *"33.106 px from the self-peak (0.046 px outside the 33.060 px NMS radius),
32.898 px from the anchor"* — every geometric quantity identical to mine, so the two runs are
demonstrably reading the same response map. But it places that rank *"of 17,626"*, and 17,626 is not
a pre-NMS peak count of either condition of this run (`hem_bbox` 53,746; `gray_bbox` 34,638, and
`gray_bbox` has **zero** annulus survivors). So a rank measured in one list is quoted against a
candidate-count range spanning a different one, in a sentence that presents both as derived here. The conclusion the
paragraph exists to support is unaffected — 4444 is as far outside the top 100 as 3282 — so this is
a Tier 1 **attribution** finding rather than a Tier 1 *result* finding, but it is Tier 1 because
the prose's own claim about its evidentiary status ("checked directly here from the stored scores,
not assumed") is false for the only two numbers in it a reader would carry away.

The "still true on 8/14 ROIs" claim about half-pixel anchors **does** reproduce — for `hem_bbox`
(8/14). For `gray_bbox` it is 9/14, and for at-least-one-condition 12/14; the prose does not say
which condition it means (F6).

---

### Tier 2

#### F3 — "`gray_inverted`'s gate never refuses a candidate `hematoxylin_od` accepts" is true of the 18 candidates the RNG visited and false of the 993-candidate pool

Cell 27:

> `hematoxylin_od`'s own gate is the constraint that actually binds everywhere in this dataset;
> `gray_inverted`'s gate never refuses a candidate `hematoxylin_od` accepts, on these 14 ROIs and
> this RNG stream. (This is a measured property of this dataset and draw, not a general claim about
> the two channels.)

The trailing hedge is real and I credit it: the *scoped* sentence is true, and I reproduced it —
18 candidates were evaluated across the 14 walks (14 accepted + 4 retries) and gray refused none of
the hem-accepted ones. But "binds everywhere in this dataset" and "a measured property of this
dataset" invite the pool-level reading, and at pool level it is false. Evaluating both gates over
the **entire** border-filtered agreement pool on all 14 ROIs:

| | candidates | hem accepts | gray accepts | hem-accept, gray-refuse | gray-accept, hem-refuse |
|---|---|---|---|---|---|
| all 14 ROIs | **993** | 815 (82.1 %) | 896 (90.2 %) | **4** | **85** |

Four counterexamples exist, on three ROIs (245.tiff ×1, 246.tiff ×2, 529.tiff ×1). The *direction*
of the claim survives comfortably — hem's gate is the more restrictive one by 81 net candidates, and
it is the binding constraint on 85 of the 89 disagreements — so the supported statement is:
*"`hematoxylin_od`'s gate is the more restrictive of the two on 993 of 993 pool candidates measured
in aggregate (815 vs 896 accepted); of the 89 candidates the two gates disagree on, 85 are refused
by hem and 4 by gray."* This costs the notebook nothing and is a sharper number than the one it has.

#### F5 — The retired baseline is declared non-comparable in cells 0 and 18 and then used as the explanatory reference in cells 8, 17 and 27

Cell 0: *"Its own deltas are therefore not comparable to any previously-reported number from this
notebook family."* Cell 18 repeats it: *"no delta in this notebook is comparable to a
previously-reported one."* Then:

- cell 17 explains its tie block against *"the ~18,000+ candidate pool"* of the uncapped notebook;
- cell 27 repeats *"not the ~18,000+ the uncapped notebook this is built on measured"*;
- cell 8 quotes *"~3282nd of the ~17,600–54,000 pre-NMS candidates in the old uncapped run"*.

A baseline retired as non-comparable cannot be the yardstick that gives a number its meaning. The
cost is concrete: F1's mechanism claim and F2's rank claim are the two places this notebook is
wrong, and both are inherited from that baseline rather than measured here. My own uncapped
extraction on 403.tiff gives 53,746 and 34,638 pre-NMS peaks — so even the "~18,000+" figure is not
this notebook's own scale.

#### F6 — Cell 17's prose asserts an expectation that the number printed in the same f-string contradicts

The cell prints:

> `nan_rate max: 0.0 | largest_tie_block max: 1 (… The tie block is expected to run **larger** than
> the uncapped notebook this is built on — a 100-candidate pre-NMS pool leaves far less room for
> `chromatin_od`/`od_contrast` to break a tie than the ~18,000+ candidate pool did.)`

`largest_tie_block == 1` is the **minimum possible value**: it means the most frequent value of the
ranking key occurs exactly once, i.e. there are **no ties at all** in any of the 336 rows. The
sentence explaining why the number should be large sits next to a number that is as small as it can
be. (Round 2 measured `largest_tie_block` max = 2 on the previous version, so the prose is inherited
from a run where the expectation was at least directionally alive.) Both numbers are correct (Tier
A, 0 divergences); only the sentence is wrong. It also shares the "less room" framing that F1
refutes, applied to a different quantity — ties *within* a ranking key rather than ties in the
hem-minus-gray delta — which is the conflation that produced F1.

#### F7 — K = 50 is outside D9's tested range, and at this cap the entire reading-burden family is undefined

D9 (2026-09-12, no amendment) reads: *"The capped pool under-delivers a K=100 budget: post-NMS,
every ROI lands at 89–99 candidates, never the full 100. **Not tested past K=30.**"* The notebook
scores K ∈ {10, 20, 30, **50**}. Delivery is fine and the notebook checks it —
`budget_delivered == budget` on 336/336, `n_detections` 94–99 ≥ 50 on 28/28, both verified here —
so this is a scope note, not a defect.

The consequence that *is* worth stating is narrower: at `MAX_PEAKS = 100`, `compare.evaluate_arms`'s
reading-depth family is almost entirely **NaN** in the artifact this notebook ships. Of 14 ROIs,
`read_50` is finite on **2** and `read_80`, `read_90`, `read_95`, `read_99`, `read_100` on **0** —
a 97-candidate list cannot deliver 80 % recall of 17–238 mitoses. Those columns are in
`..._hembbox_raw.csv` "for provenance", and nothing should ever be mined from them at this
configuration. This is not a statement about product impact: in the precision@K framing the reading
burden *is* K, fixed at 10/20/30/50 by construction, and `read_95` belongs to the recall-workload
framing the notebook explicitly puts out of scope in cell 15.

---

### Tier 3

- **F4 — the summary reports the pooled mean only; on the minimum over the 14 ROIs the
  direction is not uniform.** Cell 27's claim is explicitly scoped — *"makes **pooled** precision@K
  worse … on all three arms at all four budgets, with no exceptions … (pooled across all 14 ROIs)"*
  — and it is exactly right: 12/12 mean deltas negative, 0 divergences. This note is about a
  different statistic, which the notebook asserts nothing about and which I introduce. The project's
  stated product metric is **worst-seed** behaviour, and this run is single-seed by construction, so
  worst-seed is not measurable here at all; the minimum over the 14 ROIs is my proxy for it, not the
  project's metric. Under that proxy:

  | arm | Δ min-over-ROIs precision (hem − gray), K = 10 / 20 / 30 / 50 |
  |---|---|
  | `tm_score` | **+0.100 / +0.100 / +0.067 / 0.000** |
  | `chromatin_od` | 0.000 / −0.100 / −0.067 / 0.000 |
  | `od_contrast` | 0.000 / −0.150 / −0.067 / +0.020 |

  So `tm_score` — the production ranker — has its worst ROI *improve* under `hem_bbox` at 3 of 4
  budgets, driven entirely by 245.tiff, the biggest template shrink in the set (51 → 23 px, −54.9 %),
  where `tm_score` gains +0.20 / +0.25 / +0.20 / +0.12. Heavy caveats: the minimum over ROIs is not
  a paired statistic (the argmin ROI differs between conditions), each cell is driven by one ROI,
  and there is no test attached. It overturns nothing; it is a number a reader weighing a production
  change would want beside the pooled mean, and neither the notebook nor this audit can supply the
  worst-*seed* version of it.

- **F8 — the markdown says "two panels"; both figures have three.** Cell 23: *"Figure 1's two panels
  share one colour scale and Figure 2's two panels share one y-axis"* — then *"the three arms'
  magnitudes are comparable at a glance"* in the next sentence. Stale from the two-arm version. The
  claims themselves are correct and I verified both: Figure 1's `lim` is computed over all three
  arms at once (0.2667) so the panels do share a scale, and Figure 2 sets the same
  `set_ylim(−0.1, +0.1)` on each panel.
- **F9 — Figure 2's y-limit is set to exactly `max|value|`, so 7 of its 84 bars are drawn flush with
  the frame** and read as clipped/off-scale (canine lung `tm_score` K=10 at +0.100; canine soft
  tissue `tm_score` and `od_contrast` K=10, human melanoma `tm_score` K=10, human neuroendocrine
  `tm_score` K=10 and K=20 and `chromatin_od` K=20, all at −0.100). Nothing is actually clipped; a
  1.05× margin would remove the ambiguity.
- **F10 — the verification CSV is written twice and the printed count describes the intermediate
  file.** Cell 8 prints *"228 verification records → …csv"* and *"passed: 144 / 144"*; cell 10
  appends 14 `tightened_box_reproduces` records and rewrites the same path. The file on disk has
  **242** records and **158** hard checks (242 − 84 vacuous `no_cap`). Both numbers are right when
  printed; a reader going to the CSV finds neither.
- **F11 — "8/14 ROIs" does not say which condition.** It is 8/14 for `hem_bbox`, 9/14 for
  `gray_bbox`, 12/14 for at-least-one. Both are re-derived from `..._raw.csv`'s `tpl_cx`/`tpl_cy`.
- **F14 — the internal control is measured on a 90 % subset of the clicks a gray-tightening
  production run would see.** The joint gate draws only from candidates *both* gates accept: 811 of
  993 pool candidates, which is **99.5 %** of hem's own accepted population (815) but only **90.5 %**
  of gray's (896). So `gray_bbox` here is gray tightening evaluated on the clicks hematoxylin would
  also have allowed, not on gray's full production population; `hem_bbox` is essentially unaffected.
  The notebook's "not template-matched" caveat covers the template, not the click population. It
  cuts the same way as the headline — restricting gray to a subset hem also accepts can only make
  gray look like its own easier sub-population — so it moves no verdict, but it is a scope fact a
  reader should have, from numbers already in Part 0.
- **F12 — the 122 s / 263 s timings are single unreplicated measurements** (mode 11). They are
  incidental to the argument and the notebook does not lean on them, but there is no repeat count
  and no statement of what else was running. My own independent pass over the same 14 ROIs and both
  conditions — which does strictly more work — took 200–249 s across four runs, so the order of
  magnitude is right.
- **F13 — "one in twelve is inside what chance alone would produce" treats 12 highly dependent tests
  as if they were independent.** The three arms re-rank the same pool and the four budgets are
  nested prefixes of the same list, so the effective number of independent tests is far below 12.
  The direction of the error is conservative — it makes the notebook *less* likely to claim a
  finding, and its conclusion is "not resolvable" — so this weakens nothing. Worth one line only.

---

### §I — implementation review: what I looked for and did not find

Read in full: `compare.py`, `evaluate.py`, `seed_selection.py`, `template_match.py`, `chromatin.py`,
`dataset.py`, `channels.py`, `nms.py`, `invariants.py`, `find_and_suppress.py`.

- **Definitional drift (the highest-yield check, and the repo's own largest prior finding).** Every
  helper the notebook calls was read with its call-site defaults. `ds.image_annotations` is called
  without `category_id` — the exact argument whose permissive default caused the
  `tp_fp_separability` defect — and here that is **correct**, because `bucket_detections` splits by
  category downstream and only `HUMAN_CORRECT_LABEL` counts as a TP; passing `category_id=MITOTIC`
  would in fact be *wrong*, because it would let a detection sitting on a pathologist-rejected
  look-alike claim a mitosis behind it under one-to-one matching. `cm.chromatin_density`'s `frac`
  default (0.10) is left at default for `od51` and explicitly overridden to 0.50 for `od_ctx`, both
  matching D5's amendment. `ss.tightened_template_box`'s `method`/`center_tolerance`/`headroom_frac`
  are all left at the D8 production defaults (`"binary"`, `0`, `None`). `ev.radius_px` is in µm→px
  and is used for both NMS and matching. No unit mix-up, no index misalignment after the two merges
  (both verified row-for-row), no silent `groupby` drop.
- **Caps and floors.** `MAX_PEAKS = 100` binds on 28/28 by design and is asserted. `max_detections`
  is never applied (the notebook does not call `find_and_suppress`). `budget_delivered == budget`
  on 336/336, so no budget column is labelled K without delivering K. `caps=()` makes
  `check_no_cap` vacuous, and the notebook says so rather than counting it as evidence.
- **Leakage and self-hits.** `gt_eval` drops the click's own `ann_id` before any scoring. Self-hit
  removal references `tpl_xy` — the recentred D8 anchor — not the click, per `tightened_template_box`'s
  own "caller owns" contract; ground truth stays keyed to the click. The `seed_annulus_empty` check
  then measures what the 5 px self-hit filter cannot: survivors inside the match radius of the
  anchor. 0/28.
- **Matching.** Greedy, one-to-one, best-detection-first, radius derived per image from the TIFF's
  own mpp (29.61–33.14 px across the seven scanners here). No GT double-crediting. My independent
  matcher agrees on all 336 counts. One residual: `evaluate.greedy_match` resolves an exact
  distance tie between two free GT by KD-tree neighbour order, mine by ascending index; no tie
  occurred.
- **RNG.** `np.random.default_rng([SEED_INDEX, image_id])` — one stream per ROI, distinct by
  construction. The retry loop drops the refused row and redraws on the *same* stream, so it does
  not restart the generator. It **does** bias the draw toward gate-passing candidates, which is
  inherent to the design (and is the mechanism behind F3); both conditions share one draw, so it
  cannot bias the comparison.
- **Dtype / NaN / ties.** `nan_rate` 0.0 and `largest_tie_block` 1 on all 336 rows — no ranking is
  order-dependent. The sort is `kind="mergesort", na_position="last"`, so stability is guaranteed
  and irrelevant here.
- **Invariants the notebook asserts.** Three hard ones, all passing, all re-derived. The notebook
  does not assert `check_min_separation` post-NMS; my pools satisfy it by construction of my own NMS.

### §M — mode-by-mode disposition (Step 4)

| mode | disposition |
|---|---|
| 1 unit of analysis | **Correct as done.** The test is an exact two-sided sign-flip over the 14 ROI-level paired deltas — ROI is the exchangeable unit (D5) — not a cell-level null. I re-derived all 12 p-values by independent enumeration: 0 divergences. Effective units and p-floors added in Part 4. |
| 2 treatment × domain confound | **N/A in the usual form** — the treatment is assigned within ROI (both conditions run on every ROI from one click), so it cannot be confounded with domain. Domain-pooled deltas exist in DELTA_B and the notebook flags them as more confounded, correctly. |
| 3 recall triad | **Checked, passes.** Recall is reported alongside `n_detections` and `precision` in RAW, and `coverage_frac` = 0.0085 rules out pool saturation. Recall@K is additionally *algebraically redundant* with precision@K in this paired design (`n_gt_mitotic` identical between conditions; max residual of `Δrecall − ΔP·K/n_gt` over all 168 pairs = 9.7e−17), so setting it aside costs nothing. The non-redundant members move the **same** way or are flat: `full_list_recall` 0.371 (hem) vs 0.409 (gray), Δ = −0.038, 5/9, p = 0.122; `n_detections` Δ = −0.21, p = 0.695; `n_lookalike_in_list` Δ = −0.21 (hem attracts *fewer*), p = 0.617. **Round 1's T2-1 — "the pool/recall context moves the other way" — does not recur.** |
| 4 length-matched precision null | **N/A.** The notebook reports no precision-vs-chance comparison; it compares two conditions to each other at equal K. |
| 5 claim/evidence scope | **F5** (retired baseline used as reference); **F7** (K=50 vs D9). |
| 6 pre-registration | **N/A** — `ls Research Logs/*preregistration*` shows F1, F2, F4, F5, F6; none covers this notebook or the hem-bbox question. There is no prereg to adhere to, which is itself worth knowing before any decision rests on this. |
| 7 decision-table placement | **N/A** — no pre-registered decision rule exists. |
| 8 product relevance | **F4.** Also F7: the reading-burden family is undefined at this cap. |
| 9 machinery that cannot resolve the question | **Nothing found.** The design is the direct paired measurement; no PCA, no clustering, no surrogate. |
| 10 in-sample selection | **Largely N/A** — nothing is swept and no operating point is chosen, so there is no selection to be optimistic about. The one place it bites is the *reporting* of "1 of 12 cells at p < 0.05": that is a minimum over 12 correlated tests, and the notebook itself declines to lean on it (F13). This notebook has no held-out split, consistent with the repo-wide pattern. |
| 11 measurement/timing | **F12.** |
| 12 anything else | I went looking for: (a) a category-filter leak at `image_annotations` — absent, and correctly so; (b) the search patch for the `gray_bbox` condition being cut from `gray_inverted` instead of `hematoxylin_od` — it is not, both `read_padded_patch` calls in `_check` and the template cut in `run_condition` read from `hem`; (c) the 61 px OD pad being too small for the 121 px context window — it is not (60 needed, 61 given), and `assert n_od_nan == 0` covers it; (d) `od_ctx` being computed on the *unpadded* image for one condition and padded for the other — it is the same padded array for both; (e) the joint gate silently making the gray condition invalid — it does not, the joint click passes gray's own gate by construction, which is why the 11/14 gray-solo divergences are *earlier* picks rather than refusals; (f) whether `precision_pooled`'s idxmax/idxmin tie handling repeats the repo's `precision_pooled ties bias win-counts` incident — the notebook flags ties explicitly with `worst_roi_tied` instead of silently picking, on all 12 affected rows; (g) whether the "Otsu on `hematoxylin_od` thresholds to the dense chromatin core" claim is measurable — **it is, and it holds**: the accepted component's area is smaller under hematoxylin on **14/14** ROIs, mean 453 px vs 877 px (−48 %), at comparable solidity (hem higher on 9/14). The notebook says "evidently" where it could have said 453 vs 877. |

### §R — reconciliation with the two prior rounds (read only after the above was fixed)

Round 1 raised T1-1…T1-4, T2-1…T2-5, T3-1…T3-7 against the pre-rebuild notebook; round 2 raised
V1…V7 against a correction pass that it found had not applied them. Against the **current** content:

| prior finding | status now |
|---|---|
| R1 T1-1 — the 403 mechanism is contradicted by the response map | **Fixed, and I confirm the replacement.** The notebook now gives the half-pixel-anchor account, and my uncapped trace reproduces its geometry exactly (0.046 px, 33.106 px, 32.898 px). What survives is F2 — the *rank* it quotes. |
| R1 T1-2 — "consistently worse" not resolvable at n = 11 | **Fixed.** The notebook now runs the ROI-level exact sign-flip itself and states "not resolvable at this design" in its own summary. I reproduce all 12 p-values. |
| R1 T1-3 — headline depends on excluding three ROIs | **Fixed.** The joint gate removes the subset split entirely; all 14 ROIs, one view. |
| R1 T1-4 — the template-shrink causal story is unsupported | **Fixed.** The notebook now reports it as unsupported for `tm_score` and declines to resolve it. I reproduce ρ and both p's on all 12 cells. |
| R1 T2-1 — precision reported alone; the pool/recall context moves the other way | **Does not recur** — see mode 3 above. In this run the set-aside context moves the same way or is flat. I **agree with round 1 that it needed checking** and **disagree that the finding still applies**. |
| R1 T2-2 — D9 config drift | **Fixed** (`MAX_PEAKS = 100`), with F7 as the residual scope note. |
| R1 T2-3 — `od_contrast` dropped | **Fixed** — genuine third arm. |
| R1 T2-4 — "`chromatin_od` essentially flat" is unresolvable, not shown | **Fixed in substance** and now the subject of F1: the notebook explains the flatness with the wrong mechanism rather than overclaiming it. |
| R1 T2-5 — Gate 1 cannot check the channel change | **Fixed.** Cell 9 now says so in as many words. And I extend it: Gate 1 as written calls `ss.tightened_template_box` to check what `ss.tightened_template_box` produced, which is determinism with no independent content. My re-implementation of the Otsu box from its published spec reproduces all 28 `(base_size, tpl_cx, tpl_cy)` triples — **that** is what makes Gate 1 evidence, and it is outside the notebook. |
| R1 T3-5/T3-6/T3-7 — win-count stacking, per-panel colour scaling, arbitrary worst-ROI label | **All fixed** — ties are now stacked and labelled, Figure 1's scale is shared across panels, `worst_roi_tied` flags all 12 arbitrary picks. |
| R2 V1 — the four pipeline changes never landed | **Overtaken.** All five land in this rebuild; verified in source and reproduced from pixels. |
| R2 V2 — `MAX_PEAKS = 100` as specified would crash the notebook | **Confirmed as a real hazard and correctly handled.** Both dependent edits are present (inverted assert, `caps=()`); had either been omitted the run would have raised. |
| R2 V4 — the float-vs-exact delta fix applied to only half the cells | **Fixed.** `delta_exact_at_K` now drives the Spearman *and* the win/loss/tie counts; `delta_precision_at_K` survives only in the Figure 1 labels, where I checked all 168 cells agree to 2 dp and in sign with the exact delta. |
| R2 V5 — cell 26 asserts numbers no cell computes | **Recurs in narrower form as F2** (the 403 rank) and **F1** (the ~18,000 pool). Everything else in the current cell 27 is computed above it; I checked all 25 of its numbers. |
| R2 V7 — Figure 2's shared scale | **Fixed** (`set_ylim(−delta_lim, delta_lim)` on every panel), with F9 as a residual. |

Neither prior round examined the candidate-pool overlap, the pool-wide gate matrix, the worst-ROI
delta or the reading-burden definedness, so F1, F3, F4 and F7 neither confirm nor overturn them;
F1's subject (the K=10 tie explanation) did not exist before this rebuild.

---

## Part 2 — verdict per conclusion

| # | cell | conclusion | verdict |
|---|---|---|---|
| C1 | 0 | The question: does moving the Otsu bbox-tightening gate from `gray_inverted` to `hematoxylin_od` change precision@K at `MAX_PEAKS = 100`, for `tm_score` / `chromatin_od` / `od_contrast`? | **reproduces** — Tier B confirms the swap is wired into the tightening gate and nowhere else; the search channel is `hematoxylin_od` for both conditions |
| C2 | 0 | "This version applies all four [five] together, in one rerun" | **reproduces** — all five verified in the current source and reproduced from pixels (Part 0) |
| C3 | 0 | Scope: 14 ROIs, 1 seed, D1/D2/D3/D7/D9, three arms, no z-sweep, "still not template-matched" | **reproduces** — with F7 (K=50 outside D9's tested range) |
| C4 | 2 | The channel swap is safe: both channels are more-object-higher, `tighten_box_otsu` never uses `_INV`, per-patch min-max normalisation removes the scale difference | **reproduces** — verified in source and by re-implementation; my own min-max→`THRESH_BINARY+THRESH_OTSU` path reproduces all 28 boxes |
| C5 | 6 | The 28-row per-(ROI, condition) table | **reproduces** — Tier B, 0 divergences on 28/28 |
| C6 | 6 | "336 rows in 122s" | **reproduces** (rows) / **F12** (timing, n = 1) |
| C7 | 7, 8 | Three pipeline invariants: pool ≥ max budget, `budget_delivered == budget`, `MAX_PEAKS` binds | **reproduces** — 336/336, 336/336, 28/28 |
| C8 | 7, 8 | `seed_annulus_empty` passes on all 28 | **reproduces** — 28/28, independently from pixels |
| C9 | 7, 8 | The 403.tiff mechanism, "re-derived … rather than asserted" | **reproduces, overstated** → **F2**. Conclusion holds (rank 4444 ≫ 100); the 0.046 px and ~54k reproduce exactly; "~3282nd" does not, and nothing in the notebook derives either |
| C10 | 8 | "228 verification records"; "passed: 144 / 144" | **reproduces** at the moment printed — **F10** (the saved file has 242 / 158) |
| C11 | 9, 10 | `tightened_template_box` reproduces inline on 14 × 2 | **reproduces** — and is upgraded from a determinism check to evidence by my independent Otsu implementation, 28/28 |
| C12 | 9, 10 | Joint click == hem-only draw 14/14, == gray-only draw 11/14 (013→254, 245→6274, 300→14581) | **reproduces** — my own gate and RNG walk give the same 14 clicks, the same 11/14, and the same three divergent ann_ids |
| C13 | 11, 12 | Table A — 84 rows of precision@K | **reproduces** — Tier A 0/1,261; Tier B from pixels 0/336 |
| C14 | 13, 14 | Table B — 168 rows; pooled == mean of the two ROIs; 12/168 worst-ROI ties flagged | **reproduces** — 0 divergences including the algebraic identity |
| C15 | 15 | The recall family "is not part of this analysis" | **reproduces** — and is defensible: recall@K is algebraically redundant here (residual 9.7e−17) and the non-redundant context moves the same way or is flat (mode 3) |
| C16 | 17 | `nan_rate` 0.0, `largest_tie_block` 1; "the tie block is expected to run larger" | **reproduces** (numbers) / **does not reproduce** (the prose) → **F6** |
| C17 | 18, 19 | `base_size`: 1/14 unchanged, median 40 → 28, median −13.3 %, 13 shrink / 0 grow | **reproduces** — all six numbers, and every per-ROI pair, from pixels |
| C18 | 18 | "Otsu on `hematoxylin_od` evidently thresholds to the dense chromatin core, not the fuller nuclear outline `gray_inverted` picks up" | **reproduces** — asserted without measurement in the notebook, but it holds: accepted-component area is smaller under hematoxylin on **14/14** ROIs, mean 453 vs 877 px, at comparable solidity |
| C19 | 18 | This run's `gray_bbox` is a fresh internal control, not the committed `half_pix_fix` baseline; no delta is comparable to a previously-reported one | **reproduces** — and is then contradicted in practice by cells 8/17/27 → **F5** |
| C20 | 20 | DELTA_A (42 rows) and DELTA_B (84 rows) | **reproduces** — 0 divergences including row order |
| C21 | 21, 22 | The exact ROI-level two-sided sign-flip is the right test; `delta_exact` is used for Spearman; the 12-row STAT_CONTEXT table | **reproduces** — all 108 values by independent enumeration; the unit is correct (D5) |
| C22 | 23, 24 | Figure 1 — per-ROI delta heatmap, shared colour scale | **reproduces** — 504 plotted quantities re-derived, 0 mismatches; scale genuinely shared |
| C23 | 23, 25 | Figure 2 — per-domain pooled-delta bars, shared y-axis | **reproduces** — 84 bar heights re-derived, 0 mismatches — **F9** (7 bars flush with the frame) |
| C24 | 26 | Figure 3 + printed win/loss/tie counts | **reproduces** — 12 stacks, all summing to 14, legend and colours correct |
| C25 | 23 | "Figure 1's **two** panels … Figure 2's **two** panels" | **does not reproduce** — both figures have three panels → **F8** |
| C26 | 27 | "Pooled precision@K worse on all three arms at all four budgets, with no exceptions": −1.4…−3.6 / −2.1…−3.9 / −1.4…−3.3 pp | **reproduces** — exact on all 12 (arm, budget) cells, and the claim is explicitly scoped to the pooled mean, which is the statistic verified. F4 records that a *different* statistic (the minimum over the 14 ROIs) is not uniform, as context rather than as a correction |
| C27 | 27 | 1 of 12 cells at p < 0.05 (`od_contrast` K=30, p = 0.047); others 0.07–0.83; "not resolvable at this design" | **reproduces** — all 12 p-values exact; **F13** is a conservative-direction caveat only |
| C28 | 27 | The K=10 tie structure "is a property of the production cap … far less room for either axis to move a detection in or out of a K=10 list" | **split, per F1** — the comparative clause ("than the ~18,000+ pool") is **cannot check**: only the run cells 0/18 declare non-comparable could settle it, and I did not rebuild it. "Not evidence the two conditions rank alike" **does not reproduce**: at object level they do (5.21 of ~6.36 shared mitoses vs `tm_score`'s 3.71 of ~5.21; ties 11/14 vs 3/14). "A real, structural consequence … worth stating as a finding" **reproduces, overstated**: the stated version is measured nowhere; the supportable version is template-invariance of a location-only ranking key |
| C29 | 27 | "`hematoxylin_od`'s own gate is the constraint that binds everywhere in this dataset; `gray_inverted`'s gate never refuses a candidate `hematoxylin_od` accepts" | **reproduces, overstated** → **F3**. True on the 18 candidates walked; false pool-wide (4 counterexamples in 993). Supported: *"hem accepts 815/993 vs gray's 896/993; of the 89 disagreements, 85 are hem refusing."* |
| C30 | 27 | The template-shrink causal story still has no support for `tm_score` (\|ρ\| ≤ 0.30, p ≥ 0.30); 2 of the other 8 cells nominally significant (ρ = 0.62, 0.61); "not a resolved mechanism either way" | **reproduces** — all 12 ρ, asymptotic p and permutation p exact; the characterisation is accurate and appropriately hedged |
| C31 | 27 | The 403 exception does not recur, "and the reason is checked, not asserted" | **reproduces, overstated** → **F2** (the reason is asserted; one inequality is checked) |
| C32 | 27 | "Free simplification": a real code simplification, not behaviourally free | **reproduces, overstated** — "not free" is right (`base_size` shrinks on 13/14, component area on 14/14); "pooled precision moves in the same (negative) direction on every arm and budget" is exact for the pooled mean (F4 notes the min-over-ROIs is not uniform) |
| C33 | 27 | Caveats: single-seed; not template-matched; small G and tie-heavy; fresh internal control; domain-pooled deltas confounded | **reproduces** — all five correct, correctly stated, and the notebook is candid about its own limits; the p-floor point in particular is right (see Part 4) |

---

## Part 3 — what was re-run versus read

**Re-run (Tier B, from pixels — ~200 s of the ~20-minute budget, four full passes over the day):**
all 14 ROIs × both tightening conditions, end to end, with no `midog_utils` call anywhere. Per ROI:
`tifffile` load, `rgb2hed` and `255 − gray`, both channels' Otsu gate over the **entire**
border-filtered agreement pool (993 candidates), the RNG draw-and-retry walk, D8's anchor, the 73 px
patch and `base_size` template cut, replicate-padded `TM_CCOEFF`, the strided median/MAD deep floor,
`scipy.ndimage.maximum_filter` peak extraction capped at 100, greedy NMS, the 5 px self-hit filter,
`od51`/`od_ctx`/`od_contrast`, three stable descending ranks, greedy one-to-one matching and
precision@K and look-alike counts at four budgets. Plus one **uncapped** extraction per condition
on 403.tiff (53,746 and 34,638 peaks) with KD-tree NMS, for the annulus trace. Compared: **812
values, 0 divergences** — 28 × (`base_size`, `tpl_cx`, `tpl_cy`, `n_detections`, `n_gt_mitotic`)
plus 336 `tp_at_budget` plus 336 `lookalike_at_budget`.

**Re-run (Tier A, from the raw per-item artifact):** Table A, Table B, DELTA_A, DELTA_B,
STAT_CONTEXT (including an independent 2^k enumeration of the sign-flip null and an independently
seeded 2,000-permutation Spearman), RAW's internal identities, the verification CSV's counts, all
three figures' plotted series, and 43 prose numbers. Compared: 4,982 + 639 values, **0 divergences**
on the tables and figures; 4 prose flags, of which 3 are real findings (F3, F6, F8) and 1 is
explained (F10).

**Read, not re-run:** all ten `midog_utils` modules (Part 1 §I); `DECISIONS.md` D1, D3, D5 + its
amendment, D7, D8, D9; `D8_TEMPLATE_ANCHOR.md`'s pointer; both prior audit logs (after Step 4, per
protocol).

**No notebook was executed.** Tier C was not reached and was not needed: Tier A and Tier B both
returned zero divergences, so there is no unexplained discrepancy for a re-execution to resolve.
**One measurement I declined and should name:** the 28 *uncapped* passes that would test F1's
comparative clause [b] against the predecessor's ~18,000-candidate pool. At roughly 30 s per
uncapped NMS that is ~15 min on top of the ~5 min the capped sweep costs, and it would rebuild a
baseline the notebook itself declares non-comparable — so I ran the two 403.tiff uncapped passes the
F2 trace needed and stopped. That is why [b] is `cannot check` rather than a verdict.
The measurements that touch pixels are §Tier B above and the 403 trace; nothing in the audit
depends on re-running the target, and the target was never copied, modified or executed.

**Audit script.** `hembbox_precision_at_k_audit_round3.py`, run start to finish under
`/Users/mohinianand/anaconda3/bin/python3` — `EXIT=0`, `TOTAL 249s` on its final clean run — and it
regenerates all 21 `results/hembbox_precision_at_k_audit_round3_*.csv` tables from which every
number above is taken. Round 1's and round 2's scripts and tables are untouched.

**Beyond the named modes**, I looked for the seven things listed in mode 12 of §M. Six came back
clean; the seventh (the "chromatin core" claim) came back *supported with a number the notebook did
not have*, which is why C18 reads "reproduces" rather than "reproduces, overstated".

**Appendix facts in the agent file that no longer hold** (2026-09-12 snapshot vs today):

- **59 notebooks → 67**; **217 embedded PNGs → 252** (~96 MB decoded then; more now).
- **"all 233 `results/*.csv` tracked"** → **381 CSVs, 287 tracked, 94 untracked**, including all
  seven of this notebook's. The appendix's own lesson holds: commit status swung back within four
  days, and keying Tier A to it would have been wrong again.
- **Non-contiguous execution "live in five notebooks"** — not recounted; this target passes cleanly
  (1…17), so the check found nothing here either way.
- **Still true and re-verified today:** `FSConfig.tm_method` still defaults to `TM_CCOEFF_NORMED`
  while D1 selects `TM_CCOEFF` (the notebook bypasses it by passing `method=METHOD` into
  `fused_response`, so this is a dataclass-vs-decision divergence, not a defect of this notebook);
  `FSConfig.max_peaks` still 250000 while D9 sets 100 at the call site;
  `dataset.image_annotations(annotations, file_name, category_id=None)` still has the permissive
  default; `evaluate.MIDOG_RADIUS_UM` still 7.5 and D7 carries no amendment; there is **no** `CLAUDE.md`.
- **Tier B cost.** The appendix's ~8 s single-seed pass is an over-estimate at this configuration:
  with `n_angles=1, scales=(1.0,)` and 23–51 px templates, one condition's full pass (match, extract,
  NMS, chromatin) ran 1.4–4.4 s in the notebook and comparably here; the dominant per-ROI cost is
  `rgb2hed` plus the TIFF read, not `matchTemplate`.

---

## Part 4 — premises this audit inherited

Each line is a **project decision I applied, not a fact I verified**. My verdicts lean on them.

| premise | source | what would falsify it |
|---|---|---|
| **ROI is the exchangeable unit** | D5; F5 §8 | Evidence that two ROIs from the same tumour domain (or the same slide) are more alike than two from different domains to a degree that swamps the within-ROI pairing — e.g. a variance decomposition putting most of the delta variance between domains. The notebook's sign-flip, my re-derivation of it, and C21/C26/C27/C29 all rest on this. |
| **ROI-vs-stratum collinearity** | measured here | **Not collinear: exactly 2 ROIs in each of the 7 tumour domains** (`databases/MIDOG++.json`, `images/extra_valid`). So within-stratum replication is 2, not 1 — ROI-level and domain-level clustering are *different* tests here, and the degenerate case the protocol warns about does not apply. It is still only 2, so a domain-clustered test would have 7 units, not 14. |
| **NMS radius = evaluation match radius = 7.5 µm** | D7 (no amendment) | A demonstration that suppressing at the match radius deletes legitimate neighbouring detections more often than it removes duplicates — e.g. two true mitoses closer than 7.5 µm being systematically merged. Everything about pool composition, and therefore every precision number here, moves if this moves. |
| **`max_peaks = 100` before NMS is the production configuration** | D9 (2026-09-12, no amendment) | Any product requirement for K > 30, or a multi-seed re-run showing precision moving at K ≤ 30. F7 notes that K = 50 is already outside D9's tested range and that the reading-burden family is undefined at this cap. |
| **Recall needs `n_detections`, `precision` and `coverage_frac` beside it** | project convention / D4 | Nothing here; I applied it and it passed (mode 3). Falsified if `coverage_frac` were shown not to detect pool saturation. |
| **Reading burden and worst-seed behaviour are the product metrics, not the median** | AnnotateDx framing | A product decision that mean-over-ROIs performance is what ships. Note what this audit could *not* apply it to: the run is single-seed, so worst-**seed** behaviour is unmeasurable here, and F4's min-over-14-ROIs is my proxy for it, not the project's metric. F4 is downstream of this premise and is recorded as Tier 3 context for that reason. |
| **The length-matched precision null formula** | project convention | Not used — no chance baseline is claimed. |
| **Prior audits in `Research Logs/` are not independent corroboration** | protocol | I treated both prior rounds as re-derivable artifacts, read them last, and overturned one of their findings (R1 T2-1 does not recur). Their 403 geometry I re-derived rather than inherited — which is how F2 surfaced. |

**Which verdicts would change if a premise here turned out to be wrong.** Two, and only two, are
load-bearing. If **ROI is not the exchangeable unit**, C21, C27 and C29's resolution claims all
move — a domain-clustered null would have 7 units and a p-floor of 2/2⁷ = 0.0156, and the one cell
currently at p = 0.047 would almost certainly stop clearing 0.05, strengthening rather than
weakening the notebook's "not resolvable" conclusion. If **the worst case rather than the
mean is the product metric**, F4 stops being context and becomes a caveat any production decision
has to carry — though even then it is a min-over-ROIs proxy for a worst-seed quantity this
single-seed design cannot measure. Everything else in Part 2 — C1–C20, C22–C26, C28, C30–C33 — is
arithmetic, geometry or prose, re-derived from the pixels and the
annotation database, and survives any reasonable change to the list above. In particular **F1 and F2,
the two Tier 1 findings, are premise-free**: they are measurements of the notebook's own candidate
pools and response map, and they would read the same under any unit of analysis and any product
metric.

**Effective units, stated honestly.** The sign-flip enumerates only the nonzero deltas. Per cell,
`n_nonzero` ranges **3 to 13** of 14, so the attainable p-floor (2/2ᵏ) ranges from 0.000244 to
**0.25**. `chromatin_od` at K = 10 has k = 3 and sits **exactly at its floor** (p = 0.2500 with all
three deltas the same sign): that cell could not have reached significance at any effect size, which
is precisely the point the notebook's own third caveat makes and makes correctly.
