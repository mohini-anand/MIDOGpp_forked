# Verification of the correction pass on `production_seed_precision_at_k_chromatin_hem_bbox.ipynb`: every number it now prints reproduces from pixels and 6 of the 16 blocks are applied correctly, but the one thing the corrections required — the combined rerun — did not happen, 5 blocks are unapplied, 4 more are wrong or incomplete, the declined correction was added anyway, and `MAX_PEAKS = 100` as specified would have crashed the notebook on all 14 ROIs

**Scope.** Target: `production_hematoxylin_only/production_seed_precision_at_k_chromatin_hem_bbox.ipynb`
(untracked, re-executed in place 2026-09-16 13:24 after the round-1 audit's last edit at 12:44),
verified against the 16 `CORRECTION TO BE MADE` blocks and the `Correction dependencies` section of
`Research Logs/2026-09-16-hembbox-precision-at-k-audit.md`. Artifacts:
`results/precision_at_k_14roi_prodseed_chromatin_hembbox_{raw,per_roi,by_domain,verification,stat_context}.csv`
and `..._vs_halfpixfix_{per_roi,by_domain}.csv` (all seven rewritten 13:24, all untracked; the
`_stat_context` file is new since round 1); the committed baseline
`results/precision_at_k_14roi_prodseed_chromatin_halfpixfix_{raw,per_roi,by_domain}.csv`
(`fa94b90`, 2026-09-10, clean in `git status`); `databases/MIDOG++.json`; the 14
`images/extra_valid/*.tiff` pixels. Round-1 tables re-derived rather than trusted:
`results/hembbox_precision_at_k_audit_{tier_b_per_roi,mechanism_spearman,inference,table_b_recomputed,tier_b_annulus,tier_b_403_selfpeak}.csv`.
Figures: all three embedded PNGs decoded and read. Decisions read through their amendments:
`DECISIONS.md` D1, D3, D5 (+2026-09-08 amendment), D7, D8, D9 (2026-09-12, **no amendments — still
operative**). Audit script: `hembbox_precision_at_k_audit_round2.py`; tables
`results/hembbox_precision_at_k_audit_round2_*.csv` (18 of them). **Everything below is re-derived
from those artifacts and, for all 14 ROIs, from the ROI pixels themselves — never from the
notebook's printed output and never from round 1's conclusions.**

---

## Conflict of interest

`git log -1 --format='%an %ar' -- production_hematoxylin_only/production_seed_precision_at_k_chromatin_hem_bbox.ipynb`
returns nothing: the notebook has never been committed. `git status --porcelain` lists it, its three
PNGs, its seven CSVs, the round-1 audit log, the round-1 audit script and its 23 tables all as `??`,
and the last commits before them are `mohini-anand`'s from the same morning. So the notebook, the
correction pass that edited it, the round-1 audit and this verification are all plausibly the same
hand within a few hours, and only fresh context separates them.

What limits it: every table here comes from `hembbox_precision_at_k_audit_round2.py`, built from the
raw CSVs, `databases/MIDOG++.json` and the pixels — never from the notebook's output cells and never
from round 1's saved tables except to *compare against* them. The Tier B section re-implements
`rgb2hed` → the Otsu containment gate → `tightened_template_box`'s `(x0+x1−1)/2` anchor → the 73 px
patch and `base_size` cut → `TM_CCOEFF` on a replicate-padded ROI → the median/MAD deep floor →
dilate-based peak extraction with the `lexsort` tie-break → KD-tree greedy NMS → the 5 px self-hit
filter → `od51` → the stable descending rank → the greedy one-to-one matcher, in plain
numpy/cv2/skimage/sklearn. It calls no `midog_utils.compare`, `evaluate`, `nms`, `seed_selection`,
`dataset` or `chromatin`. Read-only was respected on my side too: round 1's log (12:44), script
(11:58) and tables (11:58–12:00) are untouched and all predate the 13:24 notebook run; the baseline
halfpixfix CSVs and `midog_utils/` are clean in `git status`.

What it cannot cover: the project premises in Part 4, and whether the corrections round 1 *asked
for* were the right ones. I verify application, not desirability — except where a correction is
mechanically impossible as written, which is Finding V2.

---

## Part 0 — what reproduces

**The corrected run is arithmetically sound and its numbers are the numbers the pixels support.**
224 Tier B values recomputed end to end from all 14 ROIs, 0 divergences. Every finding below is
about which corrections were applied, and about two paragraphs of cell 26 whose provenance claim is
false — not about the run.

| gate / check | result |
|---|---|
| **Execution coherence** | 28 cells (was 26), 16 code cells (was 15), `execution_count` **1…16 contiguous from 1**, 0 error outputs, 0 unrun code cells, last cell run, 3 embedded PNGs. Not vendored. The notebook was genuinely re-executed as one session. |
| **Provenance — modules** | `git status --porcelain midog_utils/` clean; last touched `74fbea0` (2026-09-16 10:16), before both the round-1 audit and the correction pass. No module drift between the two runs. |
| **Provenance — artifacts** | Notebook 13:24:05; its seven CSVs 13:24:00–13:24:01; round-1 log 12:44, script 11:58, tables 11:58–12:00; baseline CSVs 2026-09-10, committed and clean. Ordering is consistent with "corrections applied, then one rerun." All seven outputs untracked, so Tier A against them is *consistency with an uncommitted artifact* — which Tier B upgrades to **independent reproduction on 14/14 ROIs**. No `.partial` anywhere. |
| **Composition** | 18 checks, **16 pass**. 112 raw rows = 14 ROIs × 2 arms × 4 budgets; 28 per-ROI; 56 by-domain; 28/56 delta rows (inner merges lose nothing); 14 distinct ROIs matching disk; 7 domains × exactly 2 ROIs; 0 duplicate `(file,arm,budget)` keys; 0 duplicate rows; 0 NaN `tp_at_budget`; 87 verification records; 8 `stat_context` rows. The **two failures are not run defects** — 28 `no_cap` rows still present (T3-3 unapplied) and 1 `seed_annulus_empty` violation still present (T2-2 unapplied, so the condition could not change). |
| **K = 50 delivery** | `budget_delivered == budget` on **112/112**; `n_detections` 16,082–18,813, so `n_detections ≥ 50` on 112/112. But see V1: this was measured at `MAX_PEAKS = 2_000_000`, so it is **not evidence the D9 cap is safe** — the risk was never incurred. |
| **Tier A — `stat_context` CSV** | All **64 values** of the notebook's new statistical-context table (8 rows × 8 columns: both mean deltas, both W/L/T strings, both exact sign-flip p's, ρ and its p) rebuilt independently from `..._hembbox_raw.csv`, `..._halfpixfix_raw.csv` and the two per-ROI tables: **0 divergences.** |
| **Tier A — figures** | **196 values** (168 Figure-1 heatmap cells + 14 star flags + 56 Figure-2 bars, minus overlap) re-derived from the delta tables: **0 divergences.** All three PNGs decoded and read (Figures 1, 2, 3). |
| **Tier A — prose** | 17 numbers cell 26 asserts that no cell prints. **16 agree exactly**, including all four numbers imported from the declined T2-1 (pool +398.6 / 10 of 14 / p = 0.0070; coverage +0.0066 / p = 0.0060; `read_95` −755.4 / p = 0.0273; 300.tiff and 301.tiff reaching full-list recall 1.000). The one that does not is Finding V5. |
| **Tier B — all 14 ROIs from pixels** | **224 values** (`base_size`, `seed_ann_id`, `n_gt_mitotic`, `n_detections`, `tp_at_{10,20,30,50}` × 28 rows): **0 divergences.** Every seed, every pool size, every true-positive count reproduces. Pre-NMS peak counts 24,975–69,995. |
| **Tier B — the 403 trace** | Every number cell 8's corrected comment quotes reproduces from the response map: self-peak score **8.401** at **0.5 px** from the D8 anchor and it **is** the global maximum of the map; it **survives NMS** and **is** inside the 5 px self-hit radius; survivor score **2.217** at **33.106 px** from the self-peak (0.046 px outside the 33.060 px NMS radius), **32.898 px** from the anchor, **33.287 px** from the click; ranks **3281 / 393** 0-indexed (**3282nd / 394th** ordinal) of **17,626**; **53,746** pre-NMS peaks, which is cell 8's "~54k". 8 of 8 quoted quantities correct. |
| **Round 1 re-derived** | Round 1's `tier_b_per_roi` agrees with the 13:24 CSVs on 28/28 rows to within its own 4-dp storage rounding (max |Δ| 3.3e-5 at K = 30), and with my independent pixel run. Round 1's `table_b_recomputed` tied-worst-ROI cells (3 of 56: canine-lung/`tm_score`/K=10, canine-soft-tissue/`tm_score`/K=20, human-breast/`chromatin_od`/K=10) reproduce exactly and match the notebook's new `worst_roi_tied` flags. Round 1's `inference` p-values reproduce to 4 dp. Round 1's `mechanism_spearman` reproduces exactly under its own delta definition — and is **superseded**, see Finding V4. |
| **Ties and NaNs** | `nan_rate` max 0.0, `largest_tie_block` max 2 on all 112 rows. Cell 20's claim that a same-`tp` tie always subtracts to exactly `0.0` **holds**: `(delta == 0).sum()` equals `(tp_new == tp_old).sum()` on 8/8 (arm, K) cells, max residual 0.0. The sign-flip's effective-unit counts are therefore uncontaminated. |
| **Config drift (3-way)** | 11 settings against `FSConfig` defaults and `DECISIONS.md` read through amendments. **9 agree** — including `NMS_RADIUS_UM = MATCH_RADIUS_UM = ev.MIDOG_RADIUS_UM` (D7, 7.5 µm, not the 5.0 µm that once propagated) and `METHOD = cv2.TM_CCOEFF` passed explicitly so `FSConfig.tm_method`'s still-stale `TM_CCOEFF_NORMED` (= 5) default is never consulted. **2 diverge, both unchanged since round 1**: `MAX_PEAKS = 2_000_000` against D9's 100, and a 2-arm `AXES` against D5's 2026-09-08 amendment. |
| **Tier C** | **Not run.** 224/224 from pixels, 0 divergences — no divergence required it, and the notebook's own last cell records 202 s. |

---

## Part 1 — findings

Numbering is `V` for this verification round, to keep it distinct from round 1's `T` findings.

### Tier 1

#### V1 — The `Correction dependencies` section's central requirement was not met: none of the four pipeline changes landed, and the "rerun" re-ran the identical pipeline

The dependency section required T1-3 (joint-validity seed draw), T2-2 (`MAX_PEAKS = 100`), T3-3
(drop `MAX_PEAKS` from `Arm(caps=…)`) and T2-3 (`od_contrast` as a third arm) to be applied
**together** and rerun **once** over all 14 ROIs, with Table A, Table B, the deltas, both heatmaps
and the win/loss/tie counts recomputed from that single output. Measured against the notebook now on
disk (`..._round2_correction_matrix.csv`, `..._round2_config_drift.csv`):

| required change | state in the corrected notebook |
|---|---|
| `MAX_PEAKS = 100` (T2-2) | still **`MAX_PEAKS = 2_000_000`**, cell 1 |
| joint-validity seed draw (T1-3) | **absent**. Cell 5's `_check` closure calls `ss.tightened_template_box(hem, …)` and `tm.read_padded_patch(hem, …)` only. No notebook cell mentions `joint`, converts `gray_inverted`, or references `_patch_readable`. The gate is still single-channel. |
| `od_contrast` third arm (T2-3) | **absent**. `AXES = {'tm_score': 'score', 'chromatin_od': 'od51'}`; the raw CSV carries exactly those two arms. `od_contrast` appears in 3 markdown cells and 0 code cells. |
| drop `Arm(caps=(MAX_PEAKS,))` (T3-3) | **still there**, cell 5; the verification CSV still carries **28 `no_cap` rows**. |

The notebook *was* re-executed — `execution_count` 1…16 contiguous, all seven CSVs rewritten at
13:24 — but with a byte-identical pipeline configuration, so it produced identical numbers. I
confirm that directly: the 13:24 raw CSV matches round 1's independent from-pixels Tier B on 28/28
rows to within round 1's own 4-dp storage rounding, and matches my own fresh 14-ROI pixel run on
**224/224** values. `seed_ann_id` still matches the baseline on **11/14** ROIs, not 14/14 —
013.tiff (249 vs 254), 245.tiff (6317 vs 6274) and 300.tiff (14499 vs 14581) still draw a different
click, exactly as before. `n_gt_mitotic` is still identical between the two runs on all 14 ROIs
(28/28 rows, `delta = 0`), so round 1's characterisation of how confounded those three ROIs are is
unchanged.

**What this costs.** Everything downstream that the dependency section said to recompute from the
combined rerun — T1-2's fresh sign-flip, T2-4's sentence, T3-1's win range, T1-1/T3-2's conditional
branch — was instead filled in from the un-rerun data. The closing summary is honest about the two
scope gaps it leaves open (it names D9 and `od_contrast` explicitly in a "Known scope gaps"
paragraph), so this is documented rather than silent. But the experiment is **not complete**: the
headline still rests on an 11-ROI subset with the seed-mismatch confound intact, and neither arm has
been measured at the production cap.

**Caveat on my finding.** "Identical numbers" is what I measure; I cannot prove the operator did not
*intend* a rerun. The evidence is that all four config values are unchanged and all 224 pixel-level
values are unchanged.

---

#### V2 — `MAX_PEAKS = 100` as round 1 specified it would raise an `AssertionError` on the first ROI, on all 14 — the correction spec named two coordinated edits where three are needed

This is a defect in the round-1 correction spec that nobody has hit yet because T2-2 was not
applied. Cell 5 contains:

```python
centers, scores = tm.extract_peaks(fused, valid, PEAK_MIN_DISTANCE, cut, MAX_PEAKS)
n_peaks = len(centers)
assert n_peaks < MAX_PEAKS, f'{fn}: MAX_PEAKS is binding, raise it'
```

`extract_peaks` truncates with `order = np.lexsort((ys, xs, -scores))[:max_peaks]`
(`template_match.py:226`), so `len(centers) == min(n_above_floor, max_peaks)`. My Tier B measured
`n_above_floor` on every ROI (`..._round2_cap_feasibility.csv`):

| quantity | measured |
|---|---|
| pre-NMS peaks clearing the deep floor, per ROI | **24,975 – 69,995** (min 460.tiff, max 094.tiff) |
| `n_peaks` at `MAX_PEAKS = 100` | **exactly 100 on 14/14 ROIs** |
| `assert n_peaks < MAX_PEAKS` | **fails on 14/14** — `100 < 100` is False |

So T2-2 requires **three** coordinated edits, not two: the constant, `Arm(caps=(MAX_PEAKS,))`
(T3-3), **and** this assert. The dependency section named only the first two. The assert's own
message — "MAX_PEAKS is binding, raise it" — is written for a cap that is not meant to bind; under
D9 the cap is *meant* to bind, so the assert must be deleted or inverted, not merely relaxed.

This also sharpens the dependency section's K = 50 warning in the notebook's favour: D9's own text
records that the capped pool lands at **89–99 post-NMS candidates** per ROI, so K = 50 would have
been fully delivered. The risk was real but would not have bitten — and correspondingly,
`check_no_cap` at `caps=(100,)` would compare 100 against a post-NMS length of 89–99 and **not**
raise, so T3-3's stated failure mode is the less dangerous of the two. The assert is the blocker.

---

#### V3 — T1-1's mechanism replacement substitutes the cause round 1 explicitly demoted, and drops the one that discriminates; the channel attribution also survives untouched in markdown cell 7

The old three-clause mechanism is gone from cell 8 — correct. What replaced it:

> *"The real cause is geometric, not a channel-change consequence: `self_hit_radius` is a fixed
> 5.0px while `match_radius` is mpp-scaled (~33px here) — this repo's own standing gap between the
> two … exposed by a coin-flip-thin 0.046px margin on this one ROI."*

and in cell 26: *"the survivor … escapes NMS by a 0.046px margin, landing in the gap between the
fixed 5px self-hit radius and the ~33px mpp-scaled match radius — this repo's own known gap."*

Round 1's T1-1 correction specified the opposite attribution: *"What opens the annulus is the 0.5 px
gap between the D8 half-pixel anchor `(x0+x1-1)/2` and the true integer-pixel peak — a fixed
anchor-geometry effect,"* and said of the self-hit gap that it *"is what leaves the event visible,
not what creates it."* The notebook promoted the demoted condition and deleted the primary one. The
words "half-pixel", "anchor" and `(x0+x1-1)/2` appear nowhere in cell 8 or cell 26.

**The discriminating number.** `SELF_HIT_RADIUS` is the constant `5.0` for every ROI, and
`match_radius` is 29.609–33.139 px across the 14. So the 5 px / ~33 px gap — 24.6 to 28.1 px wide —
is present on **14 of 14 ROIs**, while the event occurs on **1 of 14**. A condition present in every
unit cannot be the cause of an event in one of them. The anchor offset does discriminate
(`..._round2_anchor_geometry.csv`):

| ROI | D8 anchor `(tpl_cx, tpl_cy)` | half-pixel? | `n_near_seed_annulus` |
|---|---|---|---|
| 013 | (1711.0, **2860.5**) | yes | 0 |
| 094 | (**5557.5**, 1177.0) | yes | 0 |
| 201 | (4092.0, 1978.0) | no | 0 |
| 233 | (4419.0, 2138.0) | no | 0 |
| 245 | (4769.0, 2186.0) | no | 0 |
| 246 | (**3632.5**, 3619.0) | yes | 0 |
| 300 | (1412.0, 3555.0) | no | 0 |
| 301 | (5911.0, **3552.5**) | yes | 0 |
| 402 | (989.0, 4915.0) | no | 0 |
| **403** | (2225.0, **2974.5**) | **yes** | **1** |
| 459 | (1845.0, 4002.0) | no | 0 |
| 460 | (5149.0, **3996.5**) | yes | 0 |
| 529 | (**3414.5**, **2721.5**) | yes | 0 |
| 548 | (**5163.5**, 3868.0) | yes | 0 |

8 of 14 anchors are half-pixel; all 6 integer-anchor ROIs have `n_near_seed = 0`, which is what the
geometry forces. The proof round 1 gave holds and I re-verify it on 403's pixels: NMS suppresses
within `nms_radius` of the surviving **peak pixel** (2225, 2974), the diagnostic measures within
`match_radius` of the **anchor** (2225.0, 2974.5), and D7 makes those two radii equal — so at zero
offset the two discs coincide and `n_near_seed ≡ 0` for any post-NMS survivor. The 403 survivor sits
32.898 px from the anchor and 33.106 px from the peak; the 0.208 px difference is exactly the
crescent the 0.5 px offset opens, and it uses 0.208 of the 0.5 px available.

**A second symptom of the same substitution.** Cell 8 states, side by side and without reconciling
them, that the survivor is *"33.106px from that self-peak, 0.046px outside the 33.060px NMS radius"*
and *"only 32.898px from the template centre, inside the 33.060px annulus."* Both are correct — I
reproduce both — but under the notebook's own stated mechanism they are contradictory, because
nothing in it distinguishes "self-peak" from "template centre." The 0.5 px anchor offset is the only
thing that reconciles them, and it is the thing the corrected text removed.

**The channel attribution also survived.** Round 1's correction targeted cells 8 and 24 only.
Markdown cell 7 was not on that list and still reads: *"this is the first one to surface a genuine
exception, **which is itself informative about what the channel change does to template
distinctiveness**."* That is the same causal claim T1-1 was written to remove, one cell above the
comment that now disclaims it. The notebook therefore asserts the channel attribution and its
retraction in adjacent cells.

**Caveat on my finding.** n = 1 event. Both conditions are necessary and neither is sufficient: the
anchor offset is present on 8/14 without producing the event, and the 5 px filter is what keeps
`n_near_seed` at 0 in the normal case (the self-peak is itself inside `match_radius` of the anchor,
so without that filter every ROI would register ≥ 1). What I can show is that the condition cell 8
names as "the real cause" is shared by 14/14 ROIs and so explains nothing, while the one it deleted
narrows the field to 8/14 and is provably necessary.

---

#### V4 — I overturn round 1's Spearman table: the notebook's float diagnosis is correct and its fix is right — but it applied the fix to only 8 of the audit's 16 cells, and the dropped half contains the one result round 1 flagged

This is the place where the correction pass **improved on** the audit, and it deserves to be said
first. Cell 20 claims round 1's Spearman coefficients do not reproduce because
`delta_precision_at_K` is built by subtracting two independently-rounded floats, and that this
breaks otherwise-exact ties Spearman is sensitive to. I tested the claim at bit level
(`..._round2_rounding_artifact.csv`, `..._round2_spearman.csv`). It is right, and here is the
evidence, from `tm_score` at K = 10:

| ROI | `precision_at_10_new` | `precision_at_10_old` | their difference, exactly |
|---|---|---|---|
| 402.tiff | 0.3 | 0.4 | `-0.10000000000000003` |
| 529.tiff | 0.2 | 0.3 | `-0.09999999999999998` |
| 403.tiff | 0.4 | 0.5 | `-0.09999999999999998` |
| 548.tiff | 0.4 | 0.5 | `-0.09999999999999998` |

All four have `tp_new − tp_old = −1`, i.e. they are a four-way tie. IEEE subtraction splits 402.tiff
off at 5.6e-17, and `spearmanr` ranks it 3.0 against 5.0 for the other three. Across the 8 (arm, K)
cells, **IEEE noise splits a genuine tp-tie in 9 of them** (spread ~5.6e-17 to 1.1e-16), and
separately the deliberate `round(tp/K, 4)` in cell 12 splits ties at **1e-4 at K = 30 only** (3
cells) — three orders of magnitude larger, and a mechanism cell 20 does not mention. Both push the
same way, and the fix cell 21 adopted — computing `(tp_new − tp_old) / K` from the integer counts —
is the correct one.

The consequence: **round 1's `mechanism_spearman.csv` is superseded.** I reproduce it exactly under
its own (rounded-delta, 2,000-permutation) definition — all-14 `tm_score` ρ = −0.153 / −0.162 /
−0.306 / −0.078 against round 1's −0.1526 / −0.1620 / −0.3065 / −0.0780 — and it should be read as
noise-contaminated. Under the exact definition the all-14 `tm_score` coefficients are −0.104 /
−0.110 / −0.323 / −0.025. Round 1's *conclusion* (no relationship between template shrinkage and
precision change) is unaffected and if anything strengthens.

**Where the correction pass then fell short.** Round 1's table had **16 rows** — 2 subsets × 2 arms
× 4 budgets. Cell 21's `rows_stat` loop is arm × budget only: **8 rows, all-14 only**. The
same-seed-11 subset was dropped without comment. Recomputing it under the notebook's own exact-delta
method:

| subset | arm | K | ρ | p (asymptotic) | p (2,000-permutation) |
|---|---|---|---|---|---|
| all 14 (computed by the notebook) | worst cell | `chromatin_od`/50 | −0.398 | 0.158 | 0.198 |
| **same-seed 11 (dropped)** | | **`chromatin_od`/50** | **−0.595** | **0.054** | **0.032** |
| same-seed 11 (dropped) | | `chromatin_od`/10 | +0.410 | 0.210 | 0.192 |

The one nominally significant cell round 1 reported (ρ = −0.595, p = 0.040) **survives the float
fix** with an identical ρ and a permutation p of 0.032. Cell 26's sentence — *"every coefficient is
small and non-significant"* — is therefore false for the half of the grid the notebook stopped
computing, and loose even for the half it did compute (|ρ| up to 0.398). See V5.

**One undisclosed method change.** Cell 21 also switched the p-value from round 1's
2,000-permutation to `sps.spearmanr`'s asymptotic value. Cell 20 discloses only the delta-definition
change. It matters at exactly one cell — same-seed-11 `chromatin_od` at K = 50, where asymptotic
gives 0.054 and permutation 0.032, straddling 0.05.

**Caveat on my finding.** n = 14 and n = 11; one nominally significant cell in 16 uncorrected tests
is what you expect by chance, its sign is wrong for the notebook's own story, and it sits on the arm
the notebook says has no mechanical reason to move with `base_size`. The finding is about the
sentence's scope, not about a real effect.

---

### Tier 2

#### V5 — Two paragraphs of cell 26 carry a provenance claim that is false: the numbers are right, but no cell in the notebook computes them

Cell 26 opens: *"Every correction below was independently re-verified (not just copied from the
audit) before being written here — see the Statistical context section above and cell 8's
verification cell for the re-derivations."* Two paragraphs are not covered by either pointer.

**(a) The declined T2-1 numbers.** Round 1's T2-1 block reads, in full: *"**CORRECTION TO BE MADE.**
None. Out of scope by product decision: only precision@K is being optimized for right now;
pool/recall/reading-depth context is not being tracked."* Cell 26 nevertheless adds a paragraph
reporting pool size +399 (10/14 larger, p = 0.007), `coverage_frac` +0.007 (p = 0.006), `read_95`
−755 for `tm_score` (p = 0.027) and 300.tiff/301.tiff reaching full-list recall 1.000. Grepping all
16 code cells: `read_95`, `coverage_frac` and `n_lookalike_in_list` appear in **zero** of them;
`full_list_recall` appears only in markdown cell 15, which still says those columns *"are not part
of this analysis."* The numbers are transcribed from round 1's log. I recomputed all four from the
raw CSVs with an exact ROI-level sign-flip (`..._round2_recall_pool.csv`) and **all four are
correct** — +398.6 / 10 of 14 / p = 0.0070; +0.0066 / p = 0.0060; −755.4 / p = 0.0273; the two ROIs
reaching 1.000 are 300 and 301. So this is a scope violation plus a false provenance claim, not a
wrong number.

**(b) The 403 numbers.** Cell 8's comment says *"re-verified independently again here before
correcting it."* Cell 8's code computes exactly two things: `seed_annulus_passed` and
`n_annulus_violations`, both from the already-stored `ROI['n_near_seed_annulus']` column. The
quantities in its comment — 8.401, 2.217, 33.106, 32.898, 33.287, ranks 3281/393, "~54k-peak" — are
computed in no cell of the notebook. Again, all eight are correct: my 403 response-map trace
reproduces every one. The claim to have re-derived them "here" is what fails.

#### V6 — T1-3's reporting half is done; its root cause and its three named removals are not

The closing summary now reports both views — the 11-ROI subset and all 14 — which addresses the
substance of round 1's T1-3 complaint, and the all-14 numbers I recompute match it exactly (K = 10
mean +0.0214, 6/6/2, exact sign-flip p = 0.7515). But T1-3's correction was explicitly a root-cause
fix with four named consequences, and three remain:

- markdown cell 22 still carries the exclusion rationale T1-3.5 said to remove entirely: *"the more
  severe confound, since the click and ground truth both changed there, not just the template."*
- Figure 1 still stars the three different-seed ROIs (T1-3.3 said to drop the annotation). On the
  render, those three starred rows — 300, 245, 013 — are still the three most strongly red rows in
  the `tm_score` panel, directly above a summary that now, to its credit, says so in words.
- the win/loss/tie figure and counts are still restricted to `DELTA_A[DELTA_A['same_seed']]`, i.e.
  11 of 14.

These are all correct *given* that the joint-seed rerun did not happen — there is still a
different-seed subset to label. They become defects only once V1 is fixed.

#### V7 — T3-6 was answered with a markdown caveat instead of the specified code change, and Figure 2 has the same problem outside the correction's scope

T3-6 said: *"cells 19/21/22 — use one shared colour scale (same vmin/vmax) across both heatmap
panels instead of computing `lim = |mat|.max()` separately per panel."* Cell 23 (the renumbered
heatmap cell) still computes `lim = max(float(np.abs(mat).max()), 1e-6)` inside the per-arm loop.
What was added is a markdown caveat in cell 22 naming the audit finding and telling the reader to
read the printed values rather than the colours. On the render, the panels are ±0.30 and ±0.10 as
before. The caveat is honest and useful; it is not the correction that was specified.

Reading Figure 2 (`hembbox_domain_delta_heatmap.png`) turns up the same issue one step further out:
its two panels are on independent y-axes (±0.20 left, ±0.05 right), so a +0.05 domain bar is
one-quarter height on the left and full height on the right. T3-6's target list did not include it.
Legend-to-series mapping and the zero line are correct on both figures.

### Tier 3

- **V8 — T3-4 was not applied at all.** `nan_rate` and `largest_tie_block` appear in **zero** cells
  of the notebook, code or markdown. `compare.py`'s docstring still calls them "mandatory rather than
  optional", they are still written to `..._hembbox_raw.csv`, and no cell reads them. Measured here:
  `nan_rate` max **0.0** and `largest_tie_block` max **2** on all 112 rows, so there is nothing to
  report — which is exactly the one-line reassurance the correction asked for.
- **V9 — T3-2 is applied slightly better than asked.** Cell 8 now gives both forms: *"Its 0-indexed
  rank is 3281 under `tm_score` and 393 under `chromatin_od` (3282nd and 394th by ordinal count) of
  17,626 candidates."* All four numbers verified from pixels. (The dependency section's conditional
  resolved to the "still occurs" branch, since `MAX_PEAKS` never changed, so the fix was required.)
- **V10 — T3-7 is applied via the "tied marker" alternative.** Cell 14 adds a `worst_roi_tied`
  boolean and prints "arbitrary pick among tied ROIs on 3/56 rows"; `worst_roi_file` itself is still
  the first row in groupby order. The correction offered "(or a 'tied' marker)", so this satisfies
  it. The 3 flagged rows match round 1's 3 exactly.
- **V11 — T3-5 is applied correctly and the render confirms it.** Ties are a visible grey third
  segment, bar height is 11 at every K in both panels, `ylim` is `(0, 12)`, one subplot per arm so
  the colour semantics are constant, legend correct.
- **V12 — one sentence-shape number drifted.** Cell 26's *"1-3 wins vs. 6-8 losses"* (T3-1's
  fallback) is correct: I measure wins 3/1/2/2 and losses 6/7/8/8 across K = 10/20/30/50. The string
  "2-3 wins" also still appears in cell 26 — but only inside the retraction clause *"not '2-3 wins'
  as an earlier version said"*, which is correct usage.

---

## Part 2 — verdict per correction

Rolled up to one row per `CORRECTION TO BE MADE` block, in the four requested labels, with my finer
sub-item label in brackets. Sub-item detail is in
`results/hembbox_precision_at_k_audit_round2_correction_matrix.csv` (31 rows).

| # | block | verdict | note |
|---|---|---|---|
| 1 | **T1-1** (403 mechanism, 4 sub-items) | **applied-incorrectly** | 1-1.2 and 1-1.3 correct (both phrases deleted). 1-1.1 and 1-1.4 replace the false mechanism with the condition round 1 explicitly demoted, and omit the anchor-geometry one; markdown cell 7's channel attribution untouched (V3) |
| 2 | **T1-2** (rewrite cell 24 from a fresh test) | **applied-incorrectly** | test is computed (new cells 20/21) and its 64 values reproduce exactly — but on the un-jointly-gated 11/3 split, i.e. the pre-rerun p-values T1-2.1 said not to reuse. Precondition DEP-1 unmet. 1-2.3's caveats applied correctly |
| 3 | **T1-3** (joint-validity seed draw) | **not-applied** | no joint helper; `_check` still single-channel; 11/14 seed match unchanged. Reporting half addressed, three named removals outstanding (V6) |
| 4 | **T1-4** (delete the causal story) | **applied-correctly** | claim deleted, not softened; replaced with "has no support in this run's own data"; the supporting Spearman is computed, and computed better than round 1 did (V4) |
| 5 | **T2-1** (none — user declined) | **not-applicable, and violated** | the declined numbers were added to cell 26 anyway, computed in no cell, under a false re-verification claim (V5a) |
| 6 | **T2-2** (`MAX_PEAKS = 100`) | **not-applied** | still `2_000_000`. And as specified it would crash: the cell-5 assert fires on 14/14 (V2) |
| 7 | **T2-3** (`od_contrast` third arm) | **not-applied** | 2 arms in config and in all output tables; no measurement of the arm or its cost |
| 8 | **T2-4** (subsumed by T1-2) | **applied-incorrectly** | the sentence shape is filled — *"only 1 of 11 same-seed ROIs moves at all … 'flat' is unresolvable at this design, not confirmed"* — with the pre-rerun 11-ROI N that T2-4 said to replace. Defect inherited from DEP-1, not independent |
| 9 | **T2-5** (Gate 1 overclaim) | **applied-correctly** | markdown cell 9 retracts the claim explicitly and states what the gate does verify; no softened version |
| 10 | **T3-1** ("2-3" → "1-3 wins") | **applied-correctly** | fallback branch is the applicable one since the rerun was skipped; 1-3/6-8 verified |
| 11 | **T3-2** (rank off-by-one) | **applied-correctly** | both 0-indexed and ordinal given; 3281/393 and 3282/394 of 17,626 all verified from pixels |
| 12 | **T3-3** (drop `MAX_PEAKS` from `caps=`) | **not-applied** | `caps=(MAX_PEAKS,)` still in cell 5; 28 vacuous `no_cap` rows still in the verification CSV |
| 13 | **T3-4** (print the two diagnostics) | **not-applied** | neither string occurs anywhere in the notebook |
| 14 | **T3-5** (ties in the win-count bars) | **applied-correctly** | verified on the render; bar height 11 at every K |
| 15 | **T3-6** (shared heatmap colour scale) | **applied-incorrectly** | per-panel `lim` unchanged; a markdown caveat was substituted for the code fix (V7) |
| 16 | **T3-7** (tied worst-ROI label) | **applied-correctly** | `worst_roi_tied` flag added — the correction's own "or a tied marker" alternative; 3/56 rows flagged, matching round 1 |

**Dependency-section requirements** (not among the 16):

| requirement | verdict |
|---|---|
| **DEP-1** — T1-3 + T2-2 + T3-3 + T2-3 applied together, rerun **once** over 14 ROIs, Table A/B/deltas/figures/win-counts recomputed from that output | **not honoured** — none of the four applied; the rerun re-ran the identical pipeline (V1) |
| **DEP-2** — after the rerun, check `budget_delivered == 50` on every ROI and arm | **not-applicable** — risk never incurred. `budget_delivered == 50` on 112/112, but at the uncapped config, so it is no evidence about the cap. D9's own text (89–99 post-NMS) says K = 50 would have been delivered |
| **DEP-3** — T1-1/T3-2 conditional on whether the 403 event survives the cap | **condition unchanged** — the cap was never applied, so the event still occurs (403.tiff, `n_near_seed = 1`), and the "if it does" branch applies. T3-2 handled correctly under it (V9); T1-1 not (V3) |

---

## Ranked gap list — what is wrong and what still needs to happen

| rank | gap | what is wrong | what still needs to happen |
|---|---|---|---|
| **1** | **DEP-1 / T1-3 / T2-2 / T2-3 / T3-3 — the combined rerun** | None of the four pipeline changes was applied; the 13:24 rerun reproduced the identical pipeline and identical numbers (224/224 pixel-verified). The experiment's headline still rests on 11 of 14 ROIs with the seed-mismatch confound intact, and neither arm is measured at the D9 production cap | In one edit: (a) add a notebook-local joint-validity seed helper that accepts the first rng-ordered candidate passing `tightened_template_box` **and** `_patch_readable` under **both** `hematoxylin_od` and `gray_inverted`, and feed that `ann_id`/`click_xy` to both arms; (b) set `MAX_PEAKS = 100`; (c) delete `MAX_PEAKS` from `Arm(caps=…)`; (d) see gap 2; (e) add `od_contrast` to `AXES`. Then rerun **once** over 14 ROIs and rebuild Table A, Table B, both delta tables, both heatmaps, the win/loss/tie counts and the `stat_context` table from that single output. Verify the 11 already-agreeing ROIs draw an unchanged seed; report any ROI whose pool is exhausted under the joint gate, with that reason stated |
| **2** | **T2-2 is un-runnable as specified** | `assert n_peaks < MAX_PEAKS` in cell 5 fires on 14/14 ROIs at `MAX_PEAKS = 100`, because `extract_peaks` truncates to exactly 100 (measured pre-NMS counts 24,975–69,995) | Delete that assert, or replace it with a check that the cap bound *as intended* (`n_peaks == MAX_PEAKS`), and update its message, which currently says "raise it". Must land in the same edit as the constant |
| **3** | **T1-1 — the wrong mechanism was substituted** | Cell 8 and cell 26 name the fixed-5 px / mpp-scaled-match-radius gap as "the real cause". That gap is present on 14/14 ROIs (24.6–28.1 px wide) while the event occurs on 1/14, so it discriminates nothing. The half-pixel D8 anchor — which narrows the field to 8/14 and is provably necessary (`n_near_seed ≡ 0` at zero offset) — is absent from both cells. Cell 8 also leaves "33.106 px from the self-peak" and "32.898 px from the template centre" side by side with nothing to reconcile them | Rewrite both passages to name the 0.5 px gap between the D8 anchor `(x0+x1−1)/2` and the integer-pixel peak as what opens the annulus, keeping the 5 px self-hit radius as the *visibility* condition (a match-radius-referenced self-hit filter would have stripped the survivor). State the reconciliation explicitly: the two distances differ by 0.208 px because the anchor and the peak are 0.5 px apart |
| **4** | **Markdown cell 7's channel attribution** | Still reads *"which is itself informative about what the channel change does to template distinctiveness"* — the exact claim T1-1 exists to remove, one cell above the comment that disclaims it. Round 1's correction targeted only cells 8 and 24, so this was missed rather than refused | Delete that clause. The event is an anchor-geometry artifact reachable in any notebook in this family; it says nothing about the channel |
| **5** | **T2-1 was declined but added anyway** | Cell 26's fourth paragraph reports pool, `coverage_frac`, `read_95` and full-list-recall deltas — the correction the user explicitly declined as out of scope. No code cell computes any of them; cell 15 still says they *"are not part of this analysis"* | Either remove the paragraph (restoring the declined scope) **or** keep it and add the cell that computes it, so the claim and the code agree. Do not leave it transcribed. All four numbers are correct as stated |
| **6** | **Cell 26's re-verification claim is false for two paragraphs** | *"Every correction below was independently re-verified … see the Statistical context section above and cell 8's verification cell"* — cell 8 computes only `seed_annulus_passed` and `n_annulus_violations`; the 403 numbers (8.401, 2.217, 33.106, 32.898, 33.287, 3281/393, ~54k) are computed in no cell, as are the T2-1 numbers | Narrow the claim to what the notebook actually computes (the sign-flip, W/L/T and Spearman in cell 21), or add the cells that compute the rest. Every number involved is correct; only the provenance sentence is wrong |
| **7** | **Spearman scope and an undisclosed method switch** | Cell 21 computes 8 of the audit's 16 Spearman cells (all-14 only) and silently drops the same-seed-11 subset, which contains ρ = −0.595 at permutation p = 0.032 / asymptotic p = 0.054. Cell 26's *"every coefficient is small and non-significant"* is false for the dropped half and loose for the computed half (max \|ρ\| = 0.398). Cell 21 also switched from round 1's 2,000-permutation p to `spearmanr`'s asymptotic p without saying so | Restore the same-seed-11 rows to `rows_stat` (2 subsets × 2 arms × 4 budgets = 16), report both p methods or name the one used, and rewrite the sentence to *"no coefficient is significant on the all-14 set; on the same-seed-11 set one of eight cells (`chromatin_od`, K = 50, ρ = −0.595) is borderline, one hit in 16 uncorrected tests, wrong-signed for the proposed mechanism"* |
| **8** | **T1-3's three named removals** | Markdown cell 22 still carries the "more severe confound" exclusion rationale; Figure 1 still stars the three different-seed ROIs; the win/loss/tie figure is still same-seed-only | Correct *conditional on gap 1*. Once the joint-seed rerun exists there is no different-seed subset, so all three go together with it. Do not remove them before the rerun — they are accurate labels of the current data |
| **9** | **T3-6 answered with prose** | Cell 23 still computes `lim` per panel; a markdown caveat was substituted. Figure 2 has the same independent-y-axis issue and was outside T3-6's target list | Compute one `lim` over both arms' matrices before the loop and pass the same `vmin`/`vmax` to both `imshow` calls; do the same for Figure 2's `set_ylim`, or say in one line why per-arm scaling is wanted there |
| **10** | **T3-4 not applied** | `nan_rate` and `largest_tie_block` occur in zero cells; both are still computed by `evaluate_arms`, written to the raw CSV and discarded | Add one cell after the raw save printing `RAW['nan_rate'].max()` and `RAW['largest_tie_block'].max()`. Measured values are 0.0 and 2, so it is a one-line clean bill of health |

---

## Part 3 — what was re-run versus read

**Audit script.** `hembbox_precision_at_k_audit_round2.py`, at the repo root. It **runs start to
finish under `/Users/mohinianand/anaconda3/bin/python3` in 296 s, exit 0** (final clean run, after
the two fixes described below) and writes all 18 tables to
`results/hembbox_precision_at_k_audit_round2_*.csv`. Every number in this log comes from it. Round
1's `hembbox_precision_at_k_audit.py` and its 23 tables were read and diffed against, never
modified.

**Tier A — re-derived from persisted artifacts.** The notebook's whole `stat_context` table (64
values: both mean deltas, both W/L/T strings, both exact sign-flip p's, ρ and its p, over 8 rows),
the Figure-1 and Figure-2 plotted series (196 values), 17 prose numbers, the composition of all
seven output CSVs, the three-way config diff, and the declined-T2-1 recall/pool context. Exact
two-sided sign-flip over ROI deltas throughout, enumerating all 2^k patterns over the nonzero units
(k = 1–13 here, so exact enumeration is instant); I report the attainable floor as 2/2^k, matching
round 1's convention. **0 divergences except the one prose claim in V5/V4.**

**Tier B — recomputed from the ROI pixels (224 values, 0 divergences), ~95 s of the run.** All 14
ROIs, both arms, end to end, in plain numpy/cv2/skimage/sklearn/tifffile: `tifffile` level-0 read →
`skimage.color.rgb2hed`[:,:,0] → my own `cv2.normalize`/Otsu/`label`/`regionprops` gate with the
three sanity gates and the centre-containment gate → my own `(x0+x1−1)/2` anchor → 73 px patch,
`base_size` crop → `cv2.matchTemplate` with `TM_CCOEFF` on a replicate-padded ROI → median/MAD floor
at z = −1.5 → dilate-based peak extraction with the `lexsort` tie-break → KD-tree greedy NMS → 5 px
self-hit removal → `od51` at frac 0.10 over a 26 px replicate pad → stable descending rank → my own
greedy one-to-one matcher → precision@K. Plus the focused 403 response-map trace (global maximum,
nearest pre-NMS peak, NMS survival, the annulus survivor's four distances and its rank in both
arms), the 14-ROI anchor-geometry table, and the cap-feasibility calculation. No `midog_utils`
function appears anywhere in the Tier B path. That is well inside the twenty-minute budget; I spent
it on a complete re-run rather than spot-checks because the task turns on whether the 13:24 rerun
changed anything, which only a full comparison answers.

**Two bugs of my own, and what each one cost.** (1) My first `chromatin_density` used
`frac = 0.25`; `chromatin.DEFAULT_FRAC` is **0.10**. It changed the `chromatin_od` arm only, and
surfaced as the 403 annulus candidate ranking 774th instead of 393rd — i.e. as a disagreement with
both round 1 and the notebook, which is how I caught it. It was fixed before any Tier B number in
this log was produced. (2) More consequentially, and stated plainly because it is the same class of
imprecision V5 charges cell 26 with: **two verdicts in an earlier run of my own correction matrix
were wrong.** Two string tests were defeated by comment line-wrapping in cell 8 and cell 26, and
that run scored **T3-2 and T2-4 as unapplied**. Both are wrong — T3-2 is applied correctly and T2-4
is applied in fallback form. The matrix now strips comment markers and collapses whitespace before
testing, and the verdicts in Part 2 are from the corrected run. Anyone reading an intermediate
artifact of this audit rather than the shipped tables should discard those two rows.

**Tier C — not run.** No divergence required it: 224/224 from pixels and 64/64 on the notebook's own
new statistical table. The notebook's own last cell records 202 s. Instead I compared the notebook's
**printed cell outputs** against the CSVs on disk (cells 6, 8, 10, 12, 14, 16, 18, 19, 21, 25, 27)
and found them consistent, including the apparent mismatch between cell 8's printed "73 verification
records" and the 87-row CSV, which is correct: cell 10 appends 14 `build_seed` rows and rewrites the
file.

**Figures.** All three decoded from the notebook's own cell outputs and read as images, plus step-1a
re-derivation of their plotted series (196 values, 0 divergences). Figure 1 confirms the surviving
`*` annotation on 300/245/013 and the ±0.30 / ±0.10 per-panel scales; Figure 2 confirms the ±0.20 /
±0.05 independent y-axes; Figure 3 confirms ties are a visible third segment with bar height 11 at
every K.

**Read, not re-run.** `midog_utils/{channels,chromatin,compare,dataset,evaluate,find_and_suppress,
invariants,nms,seed_selection,template_match}.py`; `DECISIONS.md` D1/D3/D5(+amendment)/D7/D8/D9
checked for amendments with `grep -n '^## D\|^### Amendment'` — **D9 carries none, so its
`max_peaks = 100` is operative**; round 1's log, read in full because it is the specification under
verification; the sibling `Research Logs/2026-09-16-tightening-process-hem-vs-gray-audit.md` and
`production_hematoxylin_only/tightening_process_hem_vs_gray.ipynb`, read only far enough to establish
that the sibling's own correction note repeats the Spearman claim in the same all-14-only form
(*"every arm/budget non-significant, p>=0.16"*; my measured all-14 minimum is 0.158, so even the
quoted bound is rounded the wrong way) — gap 7 has propagated to a second notebook.

**What I looked for beyond the task's named checks.** (a) Whether the seed draw changed at all —
`seed_ann_id` and `n_retries` are identical on 14/14, so it did not. (b) Whether the notebook's own
`exact_sign_flip_p` is correct — it enumerates over the nonzero deltas and compares
`|Σ ± |d|| ≥ |Σd| − 1e-9`, which is the right two-sided statistic, and its 16 outputs match mine to
4 dp. (c) Whether cell 20's "ties always subtract to exactly 0.0" reassurance holds, since the
sign-flip's effective-unit count depends on it — it does, on 8/8 cells. (d) Whether cell 21's
`new_tp - old_tp` could misalign after `set_index('file_name')` — it cannot; both are unique-indexed
over the same 14 files. (e) Whether `STAT_CONTEXT`'s hardcoded `_11` / `_14` column names could
become wrong — they are accurate today but will silently lie after the joint-seed rerun, which is
worth one line in gap 1. (f) Whether the greedy matcher could double-credit a GT object or credit
the excluded seed — my independent matcher agrees on 112/112 `tp_at_budget`. (g) Whether any ROI
under-delivers K = 50 — none does, at this config. (h) Whether the round-1 tables I lean on are
themselves reproducible — `tier_b_per_roi`, `table_b_recomputed`, `inference` and
`mechanism_spearman` all re-derive, the last one only under its (superseded) delta definition, which
is V4.

**Appendix facts that no longer hold, or needed correcting.**
- **`FSConfig.tm_method` still defaults to `TM_CCOEFF_NORMED`** (integer 5) while D1 selects
  `TM_CCOEFF`. Still true at `74fbea0`. Inert here — `METHOD` is passed explicitly.
- **`FSConfig.deep_floor_z` defaults to `None`**, not to a number; the notebook sets −1.5 itself.
- **`chromatin.DEFAULT_FRAC = 0.10`**, not 0.25 — worth recording because it is the single most
  likely place an independent re-implementation of `od51` goes wrong, and it cost me one Tier B run.
- **`DECISIONS.md` D9 (2026-09-12) still carries no amendment**, so `max_peaks = 100` remains the
  operative production configuration and T2-2 is still live.
- **Tier B cost on this path: 5–9 s per ROI** end to end on this machine including the level-0
  `tifffile` read, `rgb2hed`, the match, ~18k `od51` windows and two greedy matches — consistent with
  round 1's 4–9 s and below the appendix's ~8 s for a 51 px template, because `base_size` here is
  23–41 px.
- **`results/` now holds 359 CSVs of which 72 are untracked** (measured after this round's 18 were
  written), all seven of this notebook's outputs among them. The appendix's 226/233 counts are long
  superseded, and its lesson — commit status is a provenance question, never the question of whether
  you recompute — held again.
- **The notebook family grew another member and another audit** since round 1:
  `production_hematoxylin_only/tightening_process_hem_vs_gray.ipynb` and
  `Research Logs/2026-09-16-tightening-process-hem-vs-gray-audit.md`, both corrected in the same
  13:24–13:41 window. The sibling notebook loads no `results/*.csv` in code — I checked its code
  cells — but its closing note **quotes this notebook's results in prose** ("13/14 ROIs shrink,
  median -13.3% per ROI over 14 production seeds") and repeats its corrected Spearman statement, so
  it is downstream in argument and inherits V1's scope gaps and gap 7's all-14-only framing.

---

## Part 4 — premises this verification inherited

Each line is a project decision I **applied rather than verified**.

| premise | source | what would falsify it | which verdicts lean on it |
|---|---|---|---|
| **ROI is the exchangeable unit** | D5; `Research Logs/2026-09-04-f5-preregistration.md` §8 | Two ROIs from the same domain moving together far more than two from different domains. **Measured here:** the analysis set is **2 ROIs per tumour domain, 7 domains, 14 total** (`..._round2_composition.csv`), so ROI and stratum are **not** collinear and there is within-stratum replication of 2. The full `MIDOG++.json` carries 553 images across those domains | every p-value and interval in V4/V5/V6, and my verification of the notebook's own `stat_context` |
| **Exact sign-flip over G units, not a cluster bootstrap, at G = 11–14** | the agent protocol's rule; F1's precedent at G = 7 | A demonstration that the bootstrap is calibrated at G = 11–14. Exact enumeration over 2^k is instant at k ≤ 13, so there is no cost argument for anything coarser | V4, V5, and the reproduction of all 16 notebook p-values |
| **NMS radius = evaluation match radius = 7.5 µm (D7)** | D7 | A domain packing annotations closer than 7.5 µm | **V3 primarily** — the anchor-crescent argument requires `nms_radius == match_radius`; at unequal radii the geometry changes and the whole 403 analysis would need redoing |
| **D8's anchor is the bbox pixel centre `(x0+x1−1)/2`, a half-pixel whenever the bbox side is even** | `D8_TEMPLATE_ANCHOR.md`, status current | A demonstration that half-pixel anchoring costs more than the containment guarantee it buys. **Measured here: 8 of 14 anchors are half-pixel** | **V3 primarily** |
| **`self_hit_radius = 5.0 fixed px` beside an mpp-scaled match radius** | `FSConfig`; `find_and_suppress.py`'s own comment on 26.2 px minimum annotation spacing | A measured case where 5 px deletes a legitimate neighbour, or where the 24.6–28.1 px annulus it leaves open costs a scored true positive. **Still worth a premise review** — but V3 shows it is a 14/14 condition, so it cannot explain a 1/14 event | V3 (secondary) |
| **D9's `max_peaks = 100` is the production configuration** | D9, 2026-09-12, **re-checked for amendments here: none** | A later amendment, or a multi-seed re-run showing the cap moves precision or recall at K ≤ 30 | V1, V2, and the T2-2/T3-3 verdicts |
| **`od_contrast` is a candidate ranking axis on D5's terms** | D5's 2026-09-08 amendment | A decision narrowing the axis set back to two | the T2-3 verdict |
| **Recall travels with `n_detections`, `precision`, `coverage_frac`; reading burden is a product metric** | D4; D5; the AnnotateDx framing | A product spec that fixes K and never asks for exhaustive review | V5a — but only as to whether the declined paragraph is *worth* having; my finding is about its provenance, which holds either way |
| **Round 1's audit is a specification to verify against, not settled fact** | the protocol's own rule that a prior audit is an artifact like any other | — | V4 overturns round 1's Spearman table; V2 finds a defect in its correction spec |

**Which verdicts would change if a premise here were wrong.** If D7's radius equality were dropped,
**V3's geometric argument would need redoing** — the crescent exists only because the NMS disc and
the diagnostic disc have the same radius — though the empirical discriminator (the 5 px gap present
on 14/14, the event on 1/14) survives any radius choice. If D8's anchor were rounded to the nearest
pixel, the 403 event would vanish and V3 would collapse from "the wrong cause is named" to "the
paragraph is moot." If D9 were amended away, **V1 and V2 would lose their T2-2 leg** — but DEP-1
would still be unmet on three of its four changes, and V1 would stand. If `od_contrast` were dropped
as a candidate axis, T2-3 would become not-applicable and gap 1's part (e) would fall away. **None
of Part 0 depends on any of these**: the 224 Tier B values, the 64 `stat_context` values and the 196
figure values are arithmetic and pixels, and they reproduce under any premise.
