# Audit of `production_seed_precision_at_k_chromatin_hem_bbox.ipynb`: every number reproduces from raw pixels, but the headline "consistently worse" survives neither an ROI-level test nor the three ROIs it excludes, and the 403 mechanism it offers is contradicted by the response map

**Scope.** Target: `production_seed_precision_at_k_chromatin_hem_bbox.ipynb`
(untracked, executed in place 2026-09-16 10:56). **It was invoked at
`production_seed_precision_at_k/` and moved to `production_hematoxylin_only/` at 11:39, while this
audit was being written** — unchanged (same 395,901 bytes, same mtime 10:56:42, same 26 cells,
same `execution_count` 1…15, same 3 PNGs), so nothing below is affected; the audit script now
locates it by filename anywhere under the repo rather than by folder. Artifacts:
`results/precision_at_k_14roi_prodseed_chromatin_hembbox_{raw,per_roi,by_domain,verification}.csv`,
`results/precision_at_k_14roi_prodseed_chromatin_hembbox_vs_halfpixfix_{per_roi,by_domain}.csv`
(all six untracked), the baseline
`results/precision_at_k_14roi_prodseed_chromatin_halfpixfix_{raw,per_roi,by_domain}.csv`
(committed `fa94b90`, 2026-09-11), the five sibling
`precision_at_k_14roi_prodseed*_verification.csv` files, `databases/MIDOG++.json`, and the
`images/extra_valid/*.tiff` headers and pixels. Figures: `hembbox_roi_delta_heatmap.png`,
`hembbox_domain_delta_heatmap.png`, `hembbox_win_counts.png` (decoded from the notebook's own
cell outputs, which are byte-identical claims). Decisions read: `DECISIONS.md` D1–D9 including
D5's 2026-09-08 amendment and D9 (2026-09-12); `D8_TEMPLATE_ANCHOR.md` read through both its
superseded sections. Prior audit reconciled: `Research Logs/2026-09-12-chromatin-od-ranker-variant-audit.md`.
Audit script: `hembbox_precision_at_k_audit.py`; tables `results/hembbox_precision_at_k_audit_*.csv`
(23 of them). **Everything below is re-derived from those artifacts and, for all 14 ROIs, from the
ROI pixels themselves — never from the notebook's printed output.**

---

## Conflict of interest

`git log -1 --format='%an %ar' -- production_seed_precision_at_k/production_seed_precision_at_k_chromatin_hem_bbox.ipynb`
returns nothing: the notebook has never been committed. `git status --porcelain` lists it, its
three PNGs and its six CSVs as `??`, and the three commits immediately before it are all
`mohini-anand` within the preceding 40 minutes (`ad268f0` 10:23, `74fbea0` 10:16 — the latter
touching `midog_utils/find_and_suppress.py`, which the notebook imports). So the work under audit
and this audit are almost certainly the same hand, minutes apart, and only fresh context separates
them. `git status --porcelain midog_utils/` is clean.

What limits it: every table here is produced by `hembbox_precision_at_k_audit.py` from the raw
CSVs and the JSON, never from the notebook's output cells; the Tier B section re-implements the
channel conversion, the Otsu bbox gate, `tightened_template_box`, the template cut, the
`matchTemplate` pass, peak extraction, NMS, self-hit removal, `od51`, the stable ranking and the
greedy one-to-one matcher in plain numpy/cv2/skimage, and never calls `midog_utils.compare`,
`midog_utils.evaluate`, `midog_utils.nms`, `midog_utils.seed_selection` or `midog_utils.dataset`.

What it cannot cover: the design premises listed in Part 4 — ROI as the exchangeable unit, 7.5 µm
as both NMS and match radius, precision@K on one click per ROI as the read-out, reading burden and
worst-seed behaviour as the product metrics. Those go to `premise-reviewer` and to a reader who
did not make them, not to me.

---

## Part 0 — what reproduces

**The engineering is sound. The published statistics are the statistics the code in the working
tree computes, and — on all 14 ROIs re-run from pixels — the statistics the pixels support.** Every
finding in Part 1 is about inference, attribution and scope, not about the run.

| gate / check | result |
|---|---|
| **Execution coherence** | 26 cells, 15 code cells, `execution_count` **1…15 contiguous from 1**, 0 error outputs, 0 unrun code cells, last cell run, 3 embedded PNGs. Not vendored (imports `midog_utils`, reads `../images/extra_valid` and `../databases/MIDOG++.json`). |
| **Composition gate** | 22 checks, **22 pass**. 112 raw rows = 14 ROIs × 2 arms × 4 budgets; 28 per-ROI; 56 by-domain; 28/56 delta rows (inner merges lose nothing); 14 distinct ROIs matching the disk exactly; 7 domains × exactly 2 ROIs; 0 duplicate `(file,arm,budget)` keys; 0 duplicate rows; 0 NaN `tp_at_budget`; `budget_delivered == budget` on all 112; `n_detections ≥ 50` on all 112. |
| **Provenance — modules** | All 10 imported modules are **clean** in `git status` and, stripped of docstrings/comments, **AST-identical** to their state at `e1b6176` (the commit the baseline artifacts were written under, 2026-09-10 15:51) — except `find_and_suppress.py`, whose only semantic changes are two **new** `FSConfig` fields (`deep_floor_z`, `border_pad`) and a `border_pad` branch inside `find_and_suppress()`. Neither notebook constructs those fields or calls that function. **The new-vs-old comparison is therefore not confounded by module drift.** |
| **Provenance — artifacts** | The notebook and all six of its CSVs are **untracked and never committed**, so they carry no commit provenance; mtimes are ordered correctly (modules → notebook 10:56:42 → CSVs 10:56:36, same write). Tier A against them is *consistency with an uncommitted artifact* — which Tier B upgrades to *independent reproduction* on 14/14 ROIs. The baseline CSVs are committed (`fa94b90`). No `.partial` file is read anywhere. |
| **Tier A — Table A** | **336 values** (`budget_delivered`, `tp_at_K`, `precision_at_K` × 28 rows × 4 budgets) rebuilt from `..._hembbox_raw.csv`: **0 divergences.** |
| **Tier A — Table B** | **336 values** (`n_roi`, `tp_sum`, `delivered_sum`, `precision_pooled`, `precision_worst_roi`, `worst_roi_file` × 56 cells) rebuilt from `..._hembbox_raw.csv`: **0 divergences.** |
| **Tier A — deltas** | **532 values** (`precision_at_K_new/_old/delta` × 28, `same_seed` × 28, the 56 pooled-delta triples) rebuilt from Table A/B and the committed baseline: **0 divergences.** |
| **Tier A total** | **1,204 values compared, 0 divergences.** |
| **Closing-summary means** | All eight reproduce exactly: `tm_score` −0.0455 / −0.0545 / −0.0576 / −0.0364; `chromatin_od` −0.0091 / 0.0000 / −0.0030 / +0.0073. |
| **Win/loss/tie counts** | Never persisted by the notebook; recomputed from the delta table and matching the printed figure exactly (3/6/2, 1/7/3, 2/8/1, 2/8/1 and 0/1/10, 2/1/8, 1/3/7, 2/0/9). |
| **Tier B — all 14 ROIs from pixels** | **308 values** (`seed_ann_id`, `base_size`, `n_gt_mitotic`, `n_detections`, `tp_at_K`, `precision_at_K`): **0 divergences.** My own `rgb2hed` → Otsu gate → `(x0+x1−1)/2` anchor → 73 px patch → `TM_CCOEFF` → dilate-peaks → KD-tree NMS → 5 px self-hit → `od51` → stable rank → greedy matcher reproduces every seed, every pool size, every true-positive count. **The hematoxylin channel is genuinely wired into the bbox-refinement gate** — feeding `gray_inverted` there instead does not reproduce these 14 `base_size` values. |
| **Figures** | Both heatmaps are step-1a: their plotted series are re-derived in `..._figure_series.csv` from the already-verified delta tables (168 + 56 cells) and the third figure's series from the recomputed win counts. All three rendered PNGs decoded and read. Step-2 observations are in Findings T3-5/T3-6. |
| **Prose numbers not printed by any cell** | 17 checked, **14 agree** (median `base_size` 39→28, median −13.3 %, 13 shrink / 1 grow, the three quoted per-ROI size pairs, the 403 distance 33.29 px, match radius 33.06 px, pool 17,626, and "does not reach K≤50"). The three that do not are Findings T3-1 and T3-2. |
| **Verification CSV** | 87 records: 28 `no_cap`, 28 `nms_radius`, 15 `seed_annulus_empty`, 14 `build_seed_matches_inline_draw`, 2 aggregate. All pass except the one informational `seed_annulus_empty` the notebook flags. The `no_cap` leg is vacuous — Finding T3-3. |
| **Family `seed_annulus_empty` history** | The notebook's claim *"Every prior notebook in this family passed it cleanly; this is the first one to surface a genuine exception"* **reproduces**: all five sibling verification CSVs carry 14 per-ROI rows each with 0 violations. |
| **Ties and NaNs** | `nan_rate == 0.0` on all 112 rows; `largest_tie_block == 2` on all 112. No NaN-sorting or tie-ordering artifact exists in either arm. Table B's `idxmin` worst-ROI label is arbitrary in **3 of 56 cells** (T3-7). |
| **Claims with nothing persisted underneath** | Two: the 403 mechanism narrative in cell 8's comment (traced by Tier B — Finding T1-1) and the win/loss/tie series (recomputed above). Everything else traces to a CSV or to `databases/MIDOG++.json`. |
| **Gates not applicable** | **Tier C: not run** — no divergence required it, and the notebook's own last cell records 155 s, so re-execution would add nothing Tier B did not. **Step 4.2 (treatment/domain confound): N/A** — the channel change is applied to every ROI, and the analysis set has **2 ROIs per tumour domain**, so ROI and stratum are not collinear. **Step 4.4 (length-matched precision null): N/A** — no precision null is computed; every comparison is paired at fixed K. **Step 4.6/4.7 (pre-registration, decision table): N/A** — no pre-registration for this notebook exists in `Research Logs/`. **Step 4.10 (in-sample operating-point selection): N/A** — one pre-specified change is tested, no threshold or radius is swept; the analogue that *did* happen is the same-seed subset, reported as T1-3. **Step 4.11 (timing): N/A** — no timing claim is made. |

---

## Part 1 — findings

### Correction dependencies — read before applying anything below

**One combined rerun, not four separate ones.** T1-3 (joint-validity seed draw), T2-2
(`MAX_PEAKS = 100`), T3-3 (drop `MAX_PEAKS` from `Arm(caps=...)` — must land with T2-2, not after),
and T2-3 (`od_contrast` third arm) all change the same pipeline. Apply all four together, rerun once
over all 14 ROIs, and recompute Table A, Table B, the deltas, both heatmap figures, and the
win/loss/tie counts from that single output.

**Downstream of that rerun:**
- T1-2 and T2-4 both rewrite cell 24 from the new data — T1-2 explicitly covers both `tm_score` and
  `chromatin_od`. Neither should reuse any p-value, win count, or ROI-movement fraction quoted
  anywhere in this document; all of it is pre-rerun.
- T3-1's literal "2–3 → 1–3" patch is superseded by T1-2's full-sentence rewrite. Treat it as a
  fallback only, for if the rerun is skipped and cell 24 is left on its current numbers.

**K = 50 is at risk under `MAX_PEAKS = 100`.** `extract_peaks` truncates to the top 100 peaks
*before* NMS and self-hit removal (`template_match.py:226`), so the post-NMS pool for every ROI and
every arm is now capped at ≤100 candidates — and D9 itself says the capped pool "under-delivers
K = 100" and was "not tested past K = 30" (T2-2). After the rerun, check `budget_delivered` at K = 50
for every ROI and arm (`chromatin_od` and `od_contrast` now re-rank whatever survives the cap, not
~18,000 candidates); if any ROI under-delivers, report it explicitly rather than leave K = 50 looking
fully delivered. This also puts two Part 0 rows at risk — C5 ("`budget_delivered == budget` at every
K", currently reproduces on 112/112) and the composition gate's "`n_detections ≥ 50` on all 112"
check — both verified against the uncapped config; re-verify rather than assume they still hold.

**T1-1 and T3-2 are conditional on the same cap.** The annulus candidate they correct ranks 3282nd of
17,626 by `tm_score` — nowhere near a top-100 pre-NMS cut. After the rerun, check whether 403.tiff's
`seed_annulus_empty` violation still occurs.
- If it does not: the whole consequence analysis in cell 8 is moot. Drop both the corrected mechanism
  paragraph (T1-1) and the rank fix (T3-2); replace with one line noting the event no longer
  reproduces under the corrected `max_peaks`.
- If it does: don't reuse the old numbers. Recompute the candidate's position, score, distance and
  rank from the new run — `n_detections` itself will be smaller, so T3-2's fix becomes "recompute the
  rank," not "add 1 to 3281."

Everything else — T1-4, T2-1 (none), T2-5, T3-4, T3-5, T3-6, T3-7 — is an independent text or
plotting edit and does not depend on the rerun.

---

### Tier 1

#### T1-1 — The stated mechanism for the `seed_annulus_empty` violation on 403.tiff is contradicted in every clause by the response map, and the closing summary uses it as evidence for the template-shrink story

Cell 8's comment asserts, and the closing summary repeats, a three-step chain:

> *"unnormalized TM_CCOEFF is maximized wherever mean-subtracted correlation is largest, which is
> not guaranteed to be the template's own source pixel — a nearby, higher-variance neighbour can
> outscore it. Since nms_radius == match_radius, NMS then keeps that neighbour and suppresses the
> true self-hit; self-hit removal (5px) finds nothing left there to strip."*

and in the closing summary: *"a 31px template, down from 51px, that **no longer correlates most
strongly with its own source pixel**."*

Re-run from 403.tiff's pixels (`results/hembbox_precision_at_k_audit_tier_b_403_selfpeak.csv`):

| quantity | measured |
|---|---|
| template centre (D8 anchor) | (2225.0, 2974.5) |
| nearest pre-NMS peak to it | (2225.0, 2974.0), **0.5 px away** |
| its score | **8.40** |
| global maximum of the whole response map | **8.40, at that same 0.5 px offset** |
| did the self-peak survive NMS? | **yes** — it is the top-scoring point, kept first |
| was it inside `self_hit_radius = 5.0 px`? | **yes** (0.5 ≤ 5.0), so the self-hit filter removed it, exactly as designed |
| the surviving annulus candidate | (2255.0, 2988.0), score **2.22** |
| does it outscore the self-peak? | **no** (2.22 < 8.40) |
| its distance to the self-peak | **33.106 px** |
| NMS radius | **33.060 px** |
| margin by which it escaped NMS | **0.046 px** |
| its distance to the template anchor, inside the annulus radius | **0.162 px** |

Every clause fails. The 31 px template **does** correlate most strongly with its own source pixel —
that peak is the global maximum of the entire 39-megapixel map. NMS did **not** suppress it. The
5 px self-hit filter did **not** find nothing to strip.

**What actually opened the annulus is D8's half-pixel anchor, and it is provable without data.**
NMS suppresses from the *peak pixel*; the `seed_annulus_empty` diagnostic measures from `tpl_xy`,
the D8 bbox-centre anchor. Here those two points are **0.5 px apart** — (2225.0, 2974.5) vs
(2225.0, 2974.0) — because `(y0 + y1 − 1)/2` lands on a half-pixel while a peak can only land on an
integer pixel. The set of points inside radius *R* of the anchor but outside *R* of the peak is a
crescent whose radial width is bounded by that offset: by the triangle inequality
|d(p, anchor) − d(p, peak)| ≤ 0.5. The survivor sits in exactly that crescent, using **0.208 px of
the 0.5 px available** (0.162 px inside the anchor's disc, 0.046 px outside the peak's). **At zero
anchor offset the crescent has zero width, so whenever the self-peak survives NMS, `n_near_seed` is
identically 0** — which is why every prior notebook in this family passed the check, and why it was
ever expected to pass at all.

Self-hit removal cannot be the cause, because it runs *after* NMS and so cannot un-suppress
anything: the suppression footprint was already fixed by the surviving self-peak. The fixed-pixel
`self_hit_radius` survives only in a weaker, downstream form — a self-hit filter referenced to
`tpl_xy` at the **match radius** rather than at 5.0 fixed px would have stripped this survivor
post-NMS and zeroed the diagnostic, so the 28.06 px gap between the two radii is what leaves the
event *visible*, not what creates it.

The notebook's conclusion *"this is a genuine, if rare, consequence of the channel change, not a
bug in this notebook"* is therefore unsupported twice over: the mechanism it names did not occur,
and the mechanism that did occur is a property of **D8's anchor geometry**, not of which channel
gates Otsu. The closing summary's use of it — *"consistent with … the precision drop for
`tm_score` (whose whole signal is that distinctiveness) and this run's one seed_annulus_empty
violation"* — is a second claim resting on the same false premise.

**What survives.** The notebook's *consequence* check is correct and I reproduce it: the candidate
is 33.287 px from the click against a 33.060 px match radius, so it cannot be credited against the
excluded seed annotation, and it ranks far outside every budget. It changes no number in Table A or
Table B.

**Caveat on my finding.** n = 1 event. I cannot exclude that a 51 px template would have produced a
peak geometry in which this candidate never appeared — the two runs are not counterfactuals of each
other. What I can show is that the specific causal chain the notebook asserts did not occur.

**CORRECTION TO BE MADE.** Target: `production_seed_precision_at_k_chromatin_hem_bbox.ipynb`
(currently under `production_hematoxylin_only/`), cell 8's code comment and the closing-summary
cell (cell 24).

1. Cell 8 — delete the three-step mechanism comment in full: *"unnormalized TM_CCOEFF is maximized
   wherever mean-subtracted correlation is largest, which is not guaranteed to be the template's own
   source pixel — a nearby, higher-variance neighbour can outscore it. Since nms_radius ==
   match_radius, NMS then keeps that neighbour and suppresses the true self-hit; self-hit removal
   (5px) finds nothing left there to strip."* Replace it with the measured mechanism: the self-peak
   is the global maximum of the response map (score 8.40, 0.5 px from the D8 anchor), it survives
   NMS, and the 5 px self-hit filter correctly strips it exactly as designed. What opens the annulus
   is the 0.5 px gap between the D8 half-pixel anchor `(x0+x1-1)/2` and the true integer-pixel peak
   — a fixed anchor-geometry effect, not a channel-driven scoring effect. The survivor (score 2.22)
   escapes NMS by only 0.046 px into that crescent. **Conditional on the rerun** (see "Correction
   dependencies" above): if the event still occurs post-rerun, keep the consequence analysis
   (33.29 px from the click vs. 33.06 px radius, zero effect on Table A/B) but recompute its numbers
   fresh; if it does not, drop this paragraph and the consequence analysis both.
2. Cell 24 — delete: *"a 31px template, down from 51px, that no longer correlates most strongly with
   its own source pixel."* This is false; the template does correlate most strongly with its own
   source pixel (it is the global response-map maximum). Do not replace with a softened version —
   remove the sentence.
3. Cell 24 — delete the clause crediting this event to the channel/template-shrink story: *"consistent
   with … this run's one seed_annulus_empty violation."* Remove the causal linkage entirely; do not
   append a caveat to it.
4. Cell 24 — replace *"this is a genuine, if rare, consequence of the channel change, not a bug in
   this notebook"* with: the event traces to D8's half-pixel template-anchor geometry, which can
   surface in any notebook in this family whenever the self-peak survives NMS at a nonzero anchor
   offset — it is not a consequence of the hematoxylin-channel switch, and it changes no reported
   number.

---

#### T1-2 — "`tm_score` is consistently worse … at every budget" is not resolvable at n = 11 on three of the four budgets

The notebook reports no test and no interval. The exchangeable unit is the ROI (D5, F5 §8), and the
11 same-seed ROIs give paired deltas, so an **exact sign-flip over the 11 ROI deltas** (enumerating
all 2^k sign patterns over the k non-zero units) is the right test and is cheap:

| arm | K | mean Δ | W/L/T | effective units | attainable p-floor | **exact sign-flip p** | exact sign-test p | 95 % t CI (10 df) |
|---|---|---|---|---|---|---|---|---|
| `tm_score` | 10 | −0.0455 | 3/6/2 | 9 | 0.0039 | **0.3125** | 0.5078 | [−0.121, +0.030] |
| `tm_score` | 20 | −0.0545 | 1/7/3 | 8 | 0.0078 | **0.0391** | 0.0703 | [−0.101, −0.008] |
| `tm_score` | 30 | −0.0576 | 2/8/1 | 10 | 0.0020 | **0.0957** | 0.1094 | [−0.122, +0.007] |
| `tm_score` | 50 | −0.0364 | 2/8/1 | 10 | 0.0020 | **0.1348** | 0.1094 | [−0.083, +0.011] |
| `chromatin_od` | 10 | −0.0091 | 0/1/10 | **1** | **1.0000** | 1.0000 | 1.0000 | [−0.029, +0.011] |
| `chromatin_od` | 20 | +0.0000 | 2/1/8 | 3 | 0.2500 | 1.0000 | 1.0000 | [−0.026, +0.026] |
| `chromatin_od` | 30 | −0.0030 | 1/3/7 | 4 | 0.1250 | 0.8750 | 0.6250 | [−0.022, +0.016] |
| `chromatin_od` | 50 | +0.0073 | 2/0/9 | 2 | 0.5000 | 0.5000 | 0.5000 | [−0.005, +0.020] |

The design *could* have resolved a tm_score effect — the p-floor is 0.002–0.008, well below 0.05 —
and it did not, on three budgets out of four. Only K = 20 crosses 0.05 (p = 0.039, uncorrected, 1
of 4 budgets × 2 arms reported), and it is the only budget whose t interval excludes zero.

**The supported claim** is: *"`tm_score` precision falls on 6–8 of the 11 same-seed ROIs at every
budget; the drop is resolvable at K = 20 (exact ROI-level sign-flip p = 0.039, uncorrected) and not
resolvable at K = 10, 30 or 50."* "Consistently worse … by roughly 4–6 percentage points" states a
point estimate as if it were an established effect.

Leave-one-ROI-out makes the fragility concrete: **246.tiff** alone carries the K = 30 headline. Drop
it and the mean moves from −0.0576 to −0.0367; the LOO range across the 11 ROIs is
[−0.0700, −0.0367]. The sign never flips within the same-seed subset — but see T1-3, where it does.

**Caveat on my finding.** k = 8–10 effective units, so the sign-flip is itself low-powered; a null
here is "not resolved", not "no effect".

**CORRECTION TO BE MADE.** Depends on T1-3's rerun below — do this second, after that rerun
produces the joint-seed 14-ROI dataset.

1. Rerun the exact sign-flip test (same method as the table above) on the post-rerun dataset, using
   all ROIs that received a jointly-valid seed (up to 14, fewer only if some ROI's pool is exhausted
   under the joint gate — see T1-3's correction, point 4). Do not reuse the p-values in the table
   above; they were computed on the pre-rerun, 11-ROI, non-jointly-gated seeds and no longer apply.
2. Rewrite cell 24's `tm_score`/`chromatin_od` claims to match whatever the new test shows, in the
   same form as this finding's supported-claim sentence: *"`tm_score` precision falls on X of N ROIs
   at every budget; the drop is resolvable at K = ? (exact ROI-level sign-flip p = ?) and not
   resolvable at K = ?."* Do not carry forward "consistently worse … by roughly 4–6 percentage
   points" regardless of outcome — that phrasing states a point estimate as an established effect,
   which a single-seed-draw design (even a jointly-gated one) still cannot support.
3. Note explicitly in the closing summary that this remains a single click per ROI: the joint-validity
   gate removes the seed-mismatch confound but does not replicate or denoise the estimate, so any
   resolved effect should be reported with its exact p-value and this caveat, not as a settled result.

---

#### T1-3 — The headline depends entirely on excluding three ROIs, all three of which favour the change, and at K = 10 the direction reverses over all 14

The notebook restricts the closing summary to the 11 same-seed ROIs and never reports the all-14
number. The effect of that restriction (`results/hembbox_precision_at_k_audit_subset_sensitivity.csv`):

| arm | K | mean Δ, same-seed (n = 11) — **what the notebook reports** | mean Δ, different-seed (n = 3) | mean Δ, all 14 | exact sign-flip p, all 14 |
|---|---|---|---|---|---|
| `tm_score` | 10 | **−0.0455** | **+0.2667** | **+0.0214** | 0.7515 |
| `tm_score` | 20 | **−0.0545** | +0.1500 | −0.0107 | 0.8145 |
| `tm_score` | 30 | **−0.0576** | +0.1111 | −0.0214 | 0.5510 |
| `tm_score` | 50 | **−0.0364** | +0.0800 | −0.0114 | 0.6702 |
| `chromatin_od` | 10–50 | −0.0091 … +0.0073 | +0.0333 … 0.0000 | 0.0000 … +0.0057 | 0.27–1.00 |

All three excluded ROIs (013, 245, 300) favour hem-based tightening on `tm_score` at **every**
budget — +0.30/+0.10/+0.03/+0.02, +0.20/+0.25/+0.20/+0.18, +0.30/+0.10/+0.10/+0.04. Over all 14
ROIs the `tm_score` penalty shrinks from 3.6–5.8 pp to 1.1–2.1 pp, **reverses sign at K = 10**, and
is nowhere near resolvable (p 0.55–0.81).

The notebook's *rationale* for the split — *"the more severe confound, since the click and ground
truth both changed there, not just the template"* — is weaker than it sounds, and its own raw CSV
says so. `n_gt_mitotic` is **identical between the two runs on all 14 ROIs** (28/28 rows,
`delta_n_gt_mitotic = 0`), so on a different-seed ROI `gt_eval` differs by exactly one annotation
swapped for another, not by size. As a fraction of ground truth that is **2 of 181 mitotic
annotations on 300.tiff (1.1 %)** and 2 of 90 on 245.tiff (2.2 %) — and 300.tiff is one of the three
strongly positive rows (+0.30/+0.10/+0.10/+0.04). On that ROI the "severe confound" is a ~1 % change
in the evaluation set; it is effectively a same-seed ROI with a different template, which is the
same situation as the eleven the notebook does report. Only 013.tiff (2 of 18, 11 %) is materially
confounded.

So the split is defensible in principle — when the click changes, `gt_eval` changes too, and the two
subsets estimate different things. What is not defensible is reporting only the subset that supports
the headline while never quantifying the other, and calling the excluded ROIs severely confounded
without measuring how much ground truth actually moved. The notebook's own Figure 1 shows the problem at a glance: in the `tm_score` panel, the
three starred rows are the three most strongly **red** (positive) rows on the plot, directly above
a closing summary that says the change is consistently worse.

The honest framing of the two estimands: *"holding the click fixed, the change costs 3.6–5.8 pp of
`tm_score` precision on 11 ROIs (resolvable only at K = 20); counting the three ROIs where the
change also moved the click — which is what adopting it in production would actually do — the
net effect over 14 ROIs is −2.1 to +2.1 pp and is not resolvable at all."*

**CORRECTION TO BE MADE.** Root-cause fix, not a reporting fix: eliminate the seed mismatch itself
by drawing one seed per ROI that both arms accept, instead of describing the split after the fact.

1. **Where the mismatch comes from.** `midog_utils/seed_selection.py:build_seed` (the production
   entry point both arms call) draws candidates from a shared-rng pool and accepts the first one for
   which `tightened_template_box` (via `tighten_box_otsu`) succeeds *under the channel it was
   called with*, retrying on refusal. Each arm calls it independently — same `[SEED_INDEX, image_id]`
   rng seed, but a fresh draw-and-retry sequence — so they only diverge when a candidate passes one
   channel's Otsu gate and fails the other's. That is the entire source of the 3-ROI split this
   finding documents.
2. **The fix.** Before either arm draws its seed, do one joint-gated draw per ROI: walk the same
   rng-ordered candidate pool and accept the first candidate for which `tightened_template_box`
   *and* `_patch_readable` both succeed under **both** `hematoxylin_od` and `gray_inverted` — not
   just one. Use that candidate's `ann_id`/`click_xy` as the shared seed identity fed to both arms.
   Implement this as a notebook-local helper (mirroring how the audit script keeps its own
   reimplementation separate from `midog_utils`) rather than modifying `build_seed` itself — this
   joint-validity requirement is specific to this A/B comparison, not a general production need.
   Downstream of the shared draw, change nothing: each arm still independently runs
   `tightened_template_box` on that click under its own channel to get its own `base_size` and
   `template_xy`, exactly as today — the fix only changes which candidate gets accepted as the seed.
3. **Rerun and recompute.** Rerun all 14 ROIs (the 11 that already agreed should draw the identical
   seed unchanged, since their existing seed already satisfies the stricter joint gate — treat any
   change on those 11 as a bug in the reimplementation, not an expected result). Recompute Table A,
   Table B, the per-ROI and per-domain deltas, the win/loss/tie counts, and both heatmap figures from
   the new run. Drop Figure 1's "starred" (different-seed) annotation entirely — there is no
   different-seed subset left.
4. **If some ROI's pool is exhausted under the joint gate** (no candidate passes both channels),
   report that ROI as excluded with the reason stated explicitly (pool exhausted under the joint
   gate, not "severe confound"), and note how many ROIs this affected. Do not silently drop it.
5. Rewrite cell 24's closing summary to report one number per arm per budget, over whatever N ROIs
   actually got a jointly-valid seed (up to 14) — not two estimands, not an 11-ROI subset. Remove the
   current exclusion rationale (*"the more severe confound, since the click and ground truth both
   changed there, not just the template"*) entirely; it no longer applies once there is no exclusion.

---

#### T1-4 — The causal story ("smaller template → less distinctive → worse `tm_score`") gets no support from the notebook's own per-ROI data

The closing summary states the mechanism as established: *"A smaller template is a less distinctive
one for unnormalized `TM_CCOEFF` matching — consistent with … the precision drop for `tm_score`."*
The notebook measured `delta_base_size` per ROI and never correlated it with `delta_precision`.
Doing so (Spearman, 2,000-permutation two-sided p, `..._mechanism_spearman.csv`):

| subset | arm | K = 10 | K = 20 | K = 30 | K = 50 |
|---|---|---|---|---|---|
| all 14 | `tm_score` | ρ = −0.153 (p 0.60) | −0.162 (0.57) | −0.307 (0.29) | −0.078 (0.80) |
| same-seed 11 | `tm_score` | ρ = −0.043 (p 0.91) | +0.067 (0.86) | −0.127 (0.70) | +0.266 (0.43) |

Nothing is significant in the eight `tm_score` cells, and the all-14 sign is the **wrong way round**
for the story: the ROIs whose templates shrank most tended to have the *less* negative precision
delta. For completeness — since showing only the supporting subset is the charge in T1-3 — the full
table has **16 rows** (2 subsets × 2 arms × 4 budgets) and exactly one is nominally significant:
same-seed `chromatin_od` at K = 50, **ρ = −0.595, p = 0.040**. That is one hit in 16 uncorrected
tests, on the arm the notebook itself says has no reason to move with `base_size`, and its sign is
also wrong for the notebook's story. The two extremes make it
concrete — **245.tiff** shrank the most (47 → 23 px, −51 %) and gained the second-most `tm_score`
precision (+0.20 at K = 10), while **246.tiff**, which alone carries the K = 30 headline
(−0.27), shrank by a mid-pack 34 % (41 → 27 px).

**Caveat on my finding.** n = 14 and n = 11; this is a weak null, not a refutation of the physics.
The finding is that the notebook asserts a mechanism as "consistent with" its data when its data
carries no signal for it, and that the check costs two lines on columns it already has.

**CORRECTION TO BE MADE.** Target: cell 24 (closing summary, "Why" explanation) of
`production_seed_precision_at_k_chromatin_hem_bbox.ipynb`.

1. Delete: *"A smaller template is a less distinctive one for unnormalized TM_CCOEFF matching —
   consistent with … the precision drop for tm_score."* This causal claim is untested and
   contradicted by the notebook's own per-ROI data: Spearman ρ between Δ`base_size` and
   Δprecision is null and wrong-signed at every budget, in both the all-14 and same-seed-11
   subsets.
2. Do not replace it with a softened version of the same causal claim. If a mechanism sentence is
   wanted, state only that no mechanism has been established from this data — the correlation
   between template shrinkage and precision change is not significant at any budget and, where
   nonzero, points the wrong way.

---

### Tier 2

#### T2-1 — Precision is reported alone, and the pool/recall context the notebook computed and set aside moves the *other* way

Cell 15 declares `recall_at_budget`, `full_list_recall` and `read_50..read_100` "an unavoidable
byproduct … not part of this analysis". At fixed K, recall@K is a monotone transform of
precision@K, so that one is genuinely redundant. The others are not, and they are informative
(`..._recall_pool_summary.csv`, exact ROI-level sign-flip, G = 14):

| column | mean, hem-bbox | mean, baseline | mean Δ | ROIs moving | exact sign-flip p |
|---|---|---|---|---|---|
| `n_detections` (pool size) | 17,729.6 | 17,331.1 | **+398.6** | 10/14 larger | **0.0070** |
| `coverage_frac` | 0.9630 | 0.9564 | **+0.0066** | 11/14 higher | **0.0060** |
| `full_list_recall` | 0.9992 | 0.9985 | **+0.0007** | 2 improve, **0 worsen** | 0.50 (k = 2) |
| `read_95`, `tm_score` | 3,186.6 | 3,941.9 | **−755.4** | 10/14 shallower | **0.0273** |
| `read_95`, `chromatin_od` | 3,591.6 | 3,663.7 | −72.1 | 4/14 shallower | 0.62 |
| `n_lookalike_in_list` | 69.07 | 68.93 | +0.14 | 2/14 more | 0.91 |

**The one resolvable, substantive move is reading depth.** For the production ranker, `read_95` —
the repo's own reading-burden metric, and a product metric under the AnnotateDx framing — falls by
**755 candidates on average, shallower on 10 of 14 ROIs, exact ROI-level sign-flip p = 0.0273**.
That is the leg this finding rests on.

The other legs are weaker, and I hold myself to the standard I apply to the notebook in T2-4. Pool
size and coverage move resolvably (p = 0.007, 0.006) but from 0.9564 to 0.9630 — both arms are near
saturation, so the gain is real and small, and does not by itself make un-budgeted recall
meaningful (the D4/D5 saturation caveat). Full-list recall has **k = 2 effective units** (p = 0.50):
two ROIs that the baseline could never fully recover — 300.tiff (0.9944) and 301.tiff (0.9954) —
reach 1.000 under hem-based tightening and none moves the other way, but at k = 2 that is
**suggestive, not resolvable**, exactly as I say of `chromatin_od` in T2-4.

The notebook's verdict "**No, this does not boost precision**" is true and narrow. Its adoption
recommendation — *"adopted knowingly as a template-sizing change with a measured `tm_score`
precision cost, not as a free channel-count reduction"* — is one-sided: there is also a measured,
ROI-resolvable **reading-depth gain** for the production ranker, sitting in the CSV the notebook
wrote and explicitly set aside.

**Caveat.** `read_95`/`read_100` are D4-secondary metrics and are NaN where an ROI never reaches the
target; I set those deltas to 0 for the sign-flip, which is conservative. `read_100` itself does not
resolve (p = 0.28).

**CORRECTION TO BE MADE.** None. Out of scope by product decision: only precision@K is being
optimized for right now; pool/recall/reading-depth context is not being tracked.

---

#### T2-2 — Config drift against D9, and the scope gap it opens for `chromatin_od`

Three-way check of the notebook's config cell against `FSConfig`'s dataclass defaults and
`DECISIONS.md` read through its amendments (`..._config_drift.csv`, 12 settings). Ten of twelve
agree, including the two the repo has a history with: `NMS_RADIUS_UM = MATCH_RADIUS_UM =
ev.MIDOG_RADIUS_UM = 7.5` (D7 — **not** the 5.0 µm that propagated through five notebooks), and
`METHOD = cv2.TM_CCOEFF` passed **explicitly** to `tm.fused_response`, so `FSConfig.tm_method`'s
still-stale `TM_CCOEFF_NORMED` default is never consulted. Seed refinement is
`ss.tightened_template_box` — the current, half-pixel-corrected anchor — and `suppress()` references
the returned template centre, which is D8's cost-1 requirement.

The one substantive divergence: **`MAX_PEAKS = 2_000_000`.** D9 (dated 2026-09-12, no amendments,
so operative on 2026-09-16) caps `extract_peaks`'s `max_peaks` at **100** pre-NMS, and names by
file this notebook's direct ancestor — *"the `MAX_PEAKS` constant in
`production_seed_precision_at_k_chromatin_half_pix_fix.ipynb` cell 1, currently `2_000_000`"* — as
the value to change. This notebook inherits the superseded value and runs at pools of 16,082–18,813.

This is **not** an internal-consistency defect: the baseline runs uncapped too, so the comparison is
apples-to-apples, and D9's own evidence says precision@10/20/30 is unchanged by the cap on 14/14
ROIs. It is a **scope** defect, and sharpened it bites harder than it first looks: D9's equivalence
evidence was measured on the **`tm_score` arm only**. Under a 100-peak pre-NMS cap, `chromatin_od`
stops being a re-rank of ~18,000 candidates and becomes a re-rank of the top-100 `tm_score` peaks —
a materially different arm. **The notebook's `chromatin_od` conclusion has no coverage at the D9
production configuration**, and its own framing ("the two arms that matter in production") does not
say so. D9 also states the capped pool under-delivers K = 100 and was "not tested past K = 30",
while this notebook reports K = 50.

**CORRECTION TO BE MADE.** Target: the config cell of
`production_seed_precision_at_k_chromatin_hem_bbox.ipynb`.

1. Set `MAX_PEAKS = 100` (matches D9, 2026-09-12) in place of `2_000_000`.
2. Rerun both arms, all 14 ROIs, under the corrected cap; recompute Table A, Table B, the deltas,
   the closing summary, and both figures. `tm_score`'s precision@10/20/30 is expected to hold
   (D9's own invariance evidence covers only these three budgets, not K = 50 — see "Correction
   dependencies" above for the K = 50 delivery risk). `chromatin_od` has never been run under this
   cap at all — report its actual numbers from this rerun rather than assuming they match the
   uncapped run.
3. This closes the scope gap by construction: once both arms run under the real production cap,
   drop the "no coverage at the D9 production configuration" caveat — it no longer applies.

---

#### T2-3 — "the two arms that matter in production" drops the arm with the highest precision@10 and @20

The notebook carries `tm_score` and `chromatin_od` and drops `od31`, `od_falloff`, `mask_od_mean`
and `od_contrast` as "out of scope". But `od_contrast` is a **candidate axis under D5's own
2026-09-08 amendment**, on the same terms as `chromatin_od`, and in the very baseline this notebook
measures against it is the best arm at two of four budgets (re-derived from
`..._halfpixfix_raw.csv`, pooled over all 14 ROIs, `..._dropped_arms.csv`):

| arm | P@10 | P@20 | P@30 | P@50 | carried into this notebook |
|---|---|---|---|---|---|
| `od_contrast` | **0.6357** | **0.6179** | 0.5500 | 0.4829 | no |
| `chromatin_od` | 0.6214 | 0.6107 | **0.5595** | **0.4886** | yes |
| `od31` | 0.5643 | 0.5500 | 0.5238 | 0.4657 | no |
| `tm_score` | 0.5000 | 0.4536 | 0.4286 | 0.3743 | yes |
| `mask_od_mean` | 0.4071 | 0.4000 | 0.4119 | 0.3843 | no |
| `od_falloff` | 0.3429 | 0.3464 | 0.3333 | 0.3143 | no |

The channel change's effect on `od_contrast` is untested, and it is the arm a reader of D5's
amendment would most want to see. This costs the notebook nothing in runtime (one extra
`chromatin_density` call per candidate at window 121). Confirms and extends
`Research Logs/2026-09-12-chromatin-od-ranker-variant-audit.md` Finding 5, which reported the same
six-arm selection breadth one level upstream.

**CORRECTION TO BE MADE.** Target: the arm set and Tables A/B/deltas/closing-summary/figures of
`production_seed_precision_at_k_chromatin_hem_bbox.ipynb`.

1. Add `od_contrast` as a third scored arm alongside `tm_score` and `chromatin_od`; rerun all 14
   ROIs, both channel settings, under this three-arm set (fold into the T2-2 rerun, since that one
   is happening anyway).
2. Recompute Table A, Table B, the deltas, the closing summary, and both figures with the third arm
   included.
3. Report `od_contrast`'s measured precision delta alongside its per-candidate compute cost (one
   extra `chromatin_density` call at window 121) — the point is to see whether the actual measured
   margin over `chromatin_od` is worth that cost, not to assume it is marginal in advance.

---

#### T2-4 — "`chromatin_od` is essentially flat" is unresolvable at this design, not confirmed

At K = 10 only **1 of 11** same-seed ROIs moves at all; at K = 50, **2 of 11**. The attainable
sign-flip p-floor is therefore 1.0 and 0.5 respectively — no observation at those budgets could have
produced a small p. The notebook's own explanation for the flatness (`OD_WINDOW` is fixed at 51 px
and unrelated to `base_size`) is mechanically sound and I do not dispute it — but the candidate
*pool* does change with the template, so flatness was not a foregone conclusion and is worth stating
as unresolved rather than shown.

**CORRECTION TO BE MADE.** Subsumed by T1-2's rewrite (see "Correction dependencies" above), which
covers `chromatin_od` too — apply T1-2, not this text verbatim. The "1–4 of 11" figure below is
pre-rerun; use it only as the sentence *shape* to fill in from the new data: *"`chromatin_od` moves
on X of N ROIs at every budget; the design cannot distinguish 'flat' from 'a small effect' there."*

---

#### T2-5 — Gate 1 is overclaimed: it cannot check the channel change

Markdown cell 9 says the `build_seed` reproduction *"doubles as a fidelity check on the channel
change itself: it isn't only confirming the inline reimplementation is correct in general, it's
confirming that reimplementation is still correct under the **new** channel argument."*

It cannot. Both the inline `_check` and `ss.build_seed` are handed the same `hem` array and both
call the same `ss.tightened_template_box`. The gate compares one path against the other, so it
would pass identically if `gray_inverted` had been passed to both, or if `to_channel` returned the
wrong array. It verifies the inline reimplementation — which is valuable — and nothing about the
channel.

The channel wiring **is** correct; that is established by Tier B here, not by Gate 1: my
independent re-run feeds `rgb2hed(...)[:,:,0]` to a locally written Otsu gate and reproduces all 14
`seed_ann_id` and all 14 `base_size` values exactly, none of which match the baseline's
`gray_inverted`-derived sizes.

**CORRECTION TO BE MADE.** Target: markdown cell 9 of
`production_seed_precision_at_k_chromatin_hem_bbox.ipynb`. Remove the claim that the `build_seed`
reproduction "doubles as a fidelity check on the channel change itself" (and the sentence
justifying it). Gate 1 verifies the inline reimplementation only; do not replace with a softened
version of the same claim.

---

### Tier 3

- **T3-1 — "2–3 wins vs. 6–8 losses at every K" is wrong at K = 20**, which has **1** win. The
  correct range is 1–3 wins vs. 6–8 losses. (Loss range is right.)
  **CORRECTION TO BE MADE:** superseded by T1-2's rewrite (see "Correction dependencies" above) —
  apply T1-2 instead. Fallback only, if the rerun is skipped: cell 24, change "2–3 wins vs. 6–8
  losses at every K" to "1–3 wins vs. 6–8 losses at every K".
- **T3-2 — off-by-one in the two quoted ranks.** Cell 8 says the 403 annulus candidate "ranks
  3281st (tm_score) / 393rd (chromatin_od)". Recomputed from pixels it is the **3282nd** and
  **394th** of 17,626. `compare._rank` assigns `rank = np.arange(n)`, i.e. 0-based; the prose
  reports those as ordinals. Harmless to the conclusion (both are ≫ 50) and worth fixing because
  0-based-rank-as-ordinal is a repeatable class of error here.
  **CORRECTION TO BE MADE:** conditional on the rerun (see "Correction dependencies" above) — this
  candidate is exactly what's at risk under `MAX_PEAKS = 100`. If the event no longer occurs
  post-rerun, this fix is moot. If it does: don't just add 1 to 3281/393 — `n_detections` will be
  smaller under the new cap, so recompute the rank fresh and report it as an ordinal (0-based rank
  + 1). Wherever `compare._rank`'s 0-based rank is reported as an ordinal elsewhere, apply the same
  +1.
- **T3-3 — the `no_cap` invariant is vacuous, and 28 of the 87 verification rows are it.**
  `Arm(caps=(MAX_PEAKS,))` makes `invariants.check_no_cap` assert that the **post-NMS** list length
  is not 2,000,000, while `MAX_PEAKS` is a **pre-NMS** cap that cannot bind (peak counts are ~10⁴).
  The check passes by construction. Confirms and extends prior-audit Finding 6, which found the same
  defect with the opposite binding behaviour (`caps=(100,)`, where the cap *did* bind pre-NMS).
  **CORRECTION TO BE MADE:** wherever the notebook builds `Arm(caps=(MAX_PEAKS,))` for the
  verification pass, drop `MAX_PEAKS` from `caps=` — `check_no_cap` tests post-NMS list length,
  which a pre-NMS cap can never bind, so the check passes by construction and should not be run
  this way. Apply together with T2-2's `MAX_PEAKS = 100` change, not after it: leaving
  `Arm(caps=(MAX_PEAKS,))` in place once `MAX_PEAKS` is 100 risks `check_no_cap` raising a false
  assertion error on any ROI whose post-NMS count happens to land exactly on 100 — the same failure
  mode prior-audit Finding 6 already documented at `caps=(100,)`.
- **T3-4 — the two mandatory diagnostics are written and never read.** `compare.py`'s docstring
  calls `largest_tie_block` and `nan_rate` "mandatory rather than optional"; both land in
  `..._hembbox_raw.csv` and no cell inspects them. Measured: `nan_rate = 0.0` on all 112 rows and
  `largest_tie_block = 2` on all 112 — so there is no silent tie-break or NaN-ordering artifact in
  either arm, which is a reassuring line the notebook could have written in one sentence.
  **CORRECTION TO BE MADE:** add one cell, after the raw table is built, printing
  `raw['nan_rate'].max()` and `raw['largest_tie_block'].max()` — a one-line confirmation that
  neither diagnostic flags a problem, in place of computing and discarding them.
- **T3-5 — Figure 3 (`hembbox_win_counts.png`) stacks losses on wins and omits ties**, so bar
  height reads as "cells that were decided", not "11 ROIs", despite the y-label *"ROI count (of 11
  same-seed)"*. At K = 50 the `chromatin_od` bar is 2 tall with 9 ties invisible; at K = 10 it is
  1 tall with 10 ties invisible. A reader cannot tell "few wins" from "almost nothing moved" —
  which is exactly the distinction T2-4 turns on. Legend-to-series mapping checked on the render
  and correct.
  **CORRECTION TO BE MADE:** cell 23 — add ties as a visible third segment in the stacked bar (or
  annotate each bar with its tie count) so bar height reflects all 11 ROIs, not only the decided
  ones.
- **T3-6 — the two heatmap panels use per-panel colour limits** (`lim = |mat|.max()`: ±0.30 left,
  ±0.10 right). Each panel carries its own colorbar and every cell is annotated with its value, so
  the figure is not misleading on close reading — but a +0.10 cell is pale pink on the left and
  near-saturated red on the right, so the two arms are not comparable at a glance.
  **CORRECTION TO BE MADE:** cells 19/21/22 — use one shared colour scale (same vmin/vmax) across
  both heatmap panels instead of computing `lim = |mat|.max()` separately per panel.
- **T3-7 — Table B's worst-ROI label is arbitrary in 3 of 56 cells.** `g.loc[g['precision'].idxmin()]`
  breaks ties to the first row in groupby order. Ties occur at canine-lung/`tm_score`/K=10 (both
  0.40), canine-soft-tissue/`tm_score`/K=20 (both 0.45) and human-breast/`chromatin_od`/K=10 (both
  0.50). The `precision_worst_roi` value is right in all three; only `worst_roi_file` is arbitrary.
  This is the repo's known `precision_pooled`-ties-bias-`idxmax` pattern in its `idxmin` form.
  **CORRECTION TO BE MADE:** cells 13/14 — when `g['precision'].idxmin()` ties, set `worst_roi_file`
  to the joined list of tied `file_name`s (or a "tied" marker) instead of silently keeping the
  first row in groupby order.

---

## Part 2 — verdict per conclusion

| # | cell | conclusion | verdict |
|---|---|---|---|
| C1 | 0 (title/question) | "Does moving the Otsu bbox-tightening step from `gray_inverted` onto `hematoxylin_od` change precision@K … the only call-site change is which channel is passed into `tightened_template_box`" | **reproduces** — Tier B confirms the swap is wired into the refinement gate and nowhere else; 14/14 seeds and base sizes reproduce from pixels with `hematoxylin_od` fed to the Otsu gate |
| C2 | 0, 2 | "This is a pure argument swap, confirmed safe … `tighten_box_otsu` per-patch min-max normalizes to uint8 before Otsu, so the raw scale difference doesn't matter, only the polarity does" | **reproduces** — `cv2.normalize(..., NORM_MINMAX)` then `THRESH_BINARY+THRESH_OTSU`, never `_INV`, verified in source and by re-implementation |
| C3 | 0, 17 | "Not built on `..._chromatin.ipynb` (stale, pre-D8-fix); built on `..._half_pix_fix.ipynb`" | **reproduces** — the baseline CSVs read are the halfpixfix ones (`fa94b90`), and all 10 modules are code-identical to their state when those were written |
| C4 | 6 | The 14-ROI per-ROI table (seed, base size, pool, radii, map stats) | **reproduces** — Tier B, 0 divergences on 14/14 |
| C5 | 8 | "the deep-floor pool clears the largest budget" and "`budget_delivered == budget` at every K" | **reproduces** — 112/112 |
| C6 | 7, 8 | "Every prior notebook in this family passed `seed_annulus_empty` cleanly; this is the first genuine exception" | **reproduces** — 5 sibling verification CSVs, 14 per-ROI rows each, 0 violations |
| C7 | 8 | The 403 violation's **mechanism**: self-hit outscored by a neighbour, NMS suppresses the self-hit, 5 px filter finds nothing | **does not reproduce** — self-peak = global max 8.40 at 0.5 px, survived NMS, was stripped by the 5 px filter; the survivor scores 2.22 and sits in the ≤ 0.5 px crescent opened by D8's half-pixel anchor, not in anything the channel did (T1-1) |
| C8 | 8 | The 403 violation's **consequence**: 33.29 px from the click vs 33.06 px radius, ranks far outside K ≤ 50, zero effect on Table A/B | **reproduces** (ranks off by one, T3-2) |
| C9 | 9, 10 | "Gate 1 … doubles as a fidelity check on the channel change itself" | **reproduces, overstated** — the gate is real and passes 14/14, but both paths share the channel argument, so it cannot test it (T2-5) |
| C10 | 11, 12 | Table A — precision@K per (ROI, arm) | **reproduces** — 336 values Tier A, 0 divergences; independently from pixels Tier B, 0 divergences |
| C11 | 13, 14 | Table B — pooled and worst-ROI precision per (domain, arm, K) | **reproduces** — 336 values, 0 divergences; worst-ROI *label* arbitrary in 3/56 (T3-7) |
| C12 | 17, 18 | "`base_size` differs on every one of the 14 ROIs, median ~13 % smaller (13 shrink, 1 grows)"; "no subset isolates a pure ranking-channel effect" | **reproduces** — median 39 → 28, median −13.3 %, 13/1/0; all quoted per-ROI pairs correct |
| C13 | 18 | "`seed_ann_id` matches on 11/14 ROIs" | **reproduces** |
| C14 | 19, 21 | Figure 1 — per-ROI delta heatmap, different-seed rows starred | **reproduces** — 168 plotted cells re-derived, 0 divergences; per-panel colour scaling noted (T3-6) |
| C15 | 19, 22 | Figure 2 — per-domain pooled-delta bars | **reproduces** — 56 plotted values re-derived, 0 divergences |
| C16 | 23 | Figure 3 + printed win/loss/tie counts, same-seed only | **reproduces** (counts) / **presentation defect** (ties omitted, losses stacked on wins, T3-5) |
| C17 | 24 | Closing-summary mean-delta table (8 values) | **reproduces** — all eight exact |
| C18 | 24 | "`tm_score` is consistently **worse** … by roughly 4–6 pp at every budget (2–3 wins vs 6–8 losses at every K)" | **reproduces, overstated** — arithmetic exact, but the ROI-level exact sign-flip resolves only K = 20 (p = 0.039); K = 10/30/50 give p = 0.31/0.10/0.13, and the win range is 1–3, not 2–3 (T1-2, T3-1). Supported: *"falls on 6–8 of 11 ROIs; resolvable at K = 20 only."* |
| C19 | 24 | "`chromatin_od` is essentially **flat**" | **reproduces, overstated** — 1–4 of 11 ROIs move; the p-floor at K = 10 is 1.0, so flatness is unresolvable, not shown (T2-4) |
| C20 | 24 | "**Why:** a smaller template is a less distinctive one for `TM_CCOEFF` — consistent with the precision drop and the 403 violation" | **does not reproduce** — Spearman(Δ`base_size`, Δprecision) is null and wrong-signed at every budget (T1-4), and the 403 leg is contradicted outright (T1-1) |
| C21 | 24 | "If production adopts this, adopt it knowingly as a template-sizing change with a measured `tm_score` precision cost, not a free channel-count reduction" | **reproduces, overstated** — the cost is real and one-sided reporting omits measured gains in pool coverage, full-list recall (2 ROIs to 1.000) and `tm_score` `read_95` (−755, p = 0.027) (T2-1); and the headline cost itself is subset-dependent (T1-3) |
| C22 | 24 | Caveats: single-seed; no template-matched subset; domain-pooled deltas more confounded | **reproduces** — all three are correct and correctly stated; the notebook is candid about its own limits |
| C23 | 0, 1 | "the two arms that matter in production" | **reproduces, overstated** — D5's 2026-09-08 amendment makes `od_contrast` a candidate axis, and it leads `chromatin_od` at P@10/P@20 in the very baseline used here (T2-3) |
| C24 | 0, 1 | "Everything else is held at D1/D2/D3/D7 per `DECISIONS.md`" | **reproduces** — and it is literally true, because D9 is not in the list; but D9 (2026-09-12) is operative and the notebook runs at the value D9 supersedes (T2-2) |

No conclusion earned **cannot check**.

---

## Part 3 — what was re-run versus read

**Audit script.** `hembbox_precision_at_k_audit.py`, at the repo root. It runs **start to finish
under `/Users/mohinianand/anaconda3/bin/python3` in 182 s** (final clean run, exit 0) and writes all
23 tables to `results/hembbox_precision_at_k_audit_*.csv`. Every number in this log comes from it.

**Tier A — re-derived from persisted artifacts (1,204 values, 0 divergences).** Table A from
`..._hembbox_raw.csv` (336); Table B pooled and worst-ROI (336); the per-ROI and per-domain deltas
against the committed baseline (532). Plus, on the same tier and not counted in that total: the
seed/`base_size` comparison, the win/loss/tie counts (which the notebook never persisted), the
closing-summary means, the recall/pool context table, the six-arm baseline pooled precision, the
family `seed_annulus_empty` history across five sibling verification CSVs, the 17 prose numbers, and
the two figure-series tables.

**Tier B — recomputed from the ROI pixels (308 values, 0 divergences), ~150 s of compute.** All 14
ROIs, both arms, end to end: `cv2.imread` → `skimage.color.rgb2hed` → my own `tighten_box_otsu`
(normalize/Otsu/label/regionprops + the three sanity gates + the containment gate) → my own
`tightened_template_box` with the `(x0+x1−1)/2` anchor → 73 px patch, `base_size` crop →
`cv2.matchTemplate` with `TM_CCOEFF` on a replicate-padded ROI → median/MAD floor at z = −1.5 →
dilate-based peak extraction with the `lexsort` tie-break → KD-tree greedy NMS → 5 px self-hit
removal → `od51` over a 26 px replicate pad → stable descending rank → my own greedy one-to-one
matcher → precision@K. No `midog_utils.compare`, `evaluate`, `nms`, `seed_selection` or `dataset`
call appears anywhere in the Tier B path. Plus the focused 403 response-map trace that produced
T1-1 (~6 s). That is well inside the twenty-minute budget; I spent it on a full re-run rather than
spot-checks because the notebook's own recorded runtime is 155 s, which made the complete
reproduction cheaper than choosing which parts to sample.

**Tier C — not run.** No divergence required it, and re-executing a 155 s notebook would add nothing
Tier B did not already establish independently. The measurements that touch pixels are the Tier B
section and the 403 trace.

**Read, not re-run.** `midog_utils/{channels,chromatin,compare,dataset,evaluate,find_and_suppress,
invariants,nms,seed_selection,template_match}.py` (source read for the implementation review, and
AST-diffed against `e1b6176`); `DECISIONS.md` D1–D9 including D5's amendment; `D8_TEMPLATE_ANCHOR.md`
through both superseded sections; the baseline notebook's config and `run_roi` cells, to establish
that the only difference is the `_check` channel argument and the arm set (the `OD_PAD` difference,
60 vs 26, is behaviourally inert: `BORDER_REPLICATE` replicates the edge pixel, so a 51 px window
reads identically at any pad ≥ 25); `Research Logs/2026-09-12-chromatin-od-ranker-variant-audit.md`,
read **after** Part 1 was formed.

**Reconciliation with the prior audit.** Its Finding 6 (`no_cap` vacuous) — **confirmed and
extended** here as T3-3, in the opposite binding regime. Its Finding 5 (selection breadth 6 upstream)
— **confirmed and extended** as T2-3, since this notebook inherits that selection and narrows it
further. Its Finding 7 (`coverage_frac` discarded) — **partially overturned in this notebook's
favour**: here `coverage_frac` *is* persisted in `..._hembbox_raw.csv`; it is simply never read,
which is T2-1's subject. Its Finding 11 (tie block 2–3, not ~10) — **independently confirmed**:
`largest_tie_block == 2` on all 112 rows here. Its Findings 1–4 and 8–10 concern
`chromatin_od_ranker_variant.ipynb` and its timing CSVs, which this notebook does not touch; **I did
not look**.

**What I looked for beyond Step 4's named modes.** (a) Whether the channel swap leaked into the
*search* channel or the `od51` ranking rather than the refinement gate — it did not; both notebooks
search on `hematoxylin_od` and score `od51` on it, and Tier B reproduces the seeds only when the
Otsu gate is fed `hematoxylin_od`. (b) Whether the differing `OD_PAD` (26 vs 60) could shift `od51`
— it cannot, per the replicate-padding argument above, and Tier B confirms it at pad 26. (c) Whether
module edits between 2026-09-10 and 2026-09-16 could have moved the baseline out from under the
comparison — AST-diff says no. (d) Whether the greedy matcher could double-credit one GT object or
credit the excluded seed — it cannot; `claimed` is one-to-one and `gt_eval` drops the seed's
`ann_id`, and my independent matcher agrees on 112/112 `tp_at_budget`. (e) Whether NaN or tie
ordering could be shaping the ranking — `nan_rate` 0.0 and tie block 2 on all 112 rows. (f) Whether
the response-map padding left any ROI pixel unreachable — the padded response is exactly H×W by
construction for a single unrotated template, so `valid.all()` is trivially true and the assert is
uninformative rather than wrong. (g) Whether the three different-seed ROIs differ systematically
from the eleven — they do, and strongly, which is T1-3. (h) Whether the seed swap on those three
ROIs changed the *size* of `gt_eval` rather than just its membership — it did not
(`delta_n_gt_mitotic = 0` on 28/28 rows), which is what lets T1-3 put a number on how confounded
they actually are. (i) Whether `n_near_seed > 0` is reachable at all when the self-peak survives
NMS — it is not, unless the template anchor is offset from the peak pixel, which is the proof
behind T1-1.

**Appendix facts that no longer hold, or needed correcting.**
- **`results/` is no longer partly untracked, but the *analysis* outputs are.** All six of this
  notebook's CSVs, the notebook itself and its three PNGs are `??`. The appendix's lesson — that
  commit status is a provenance question, never the question of whether you recompute — held
  exactly.
- **`FSConfig.tm_method` still defaults to `TM_CCOEFF_NORMED`** while D1 selects `TM_CCOEFF`. Still
  true as of `74fbea0` (2026-09-16). Inert here: the notebook passes `METHOD` explicitly.
- **`FSConfig` has gained two fields since the appendix was written** — `deep_floor_z` and
  `border_pad` (commit `f9b4a63`, 2026-09-14) — and `find_and_suppress()` a `border_pad` branch.
  Neither notebook uses them.
- **`image_annotations(anns, fn)`'s `category_id=None` default is unchanged**, and this notebook
  handles it correctly: it filters to `ds.MITOTIC` for the seed pool and deliberately keeps both
  categories in `gt_eval`, which is what `evaluate_arms` requires.
- **Tier B cost is lower than the appendix's ~8 s/ROI on this path**: 4–9 s per ROI end to end
  including `rgb2hed`, the match, ~18k `od51` windows and two greedy matches, because the single
  unrotated template at `base_size` 23–41 px is cheaper than the 51 px measurement in the appendix.
- **`DECISIONS.md` now carries D9** (2026-09-12), which the appendix predates and which is the
  live config-drift issue for this notebook family.
- **The appendix's own folder-reorganisation warning fired during this audit.** The target moved
  from `production_seed_precision_at_k/` to `production_hematoxylin_only/` at 11:39, 43 minutes
  after it was executed and while Part 1 was being written, together with its three PNGs. A
  hardcoded path in `hembbox_precision_at_k_audit.py` broke on the next clean run; it now resolves
  both notebooks with `REPO.rglob(<filename>)`. Nothing in this log derives a relationship from a
  folder name, and the file itself is byte-for-byte what was audited. A sibling notebook
  (`production_hematoxylin_only/tightening_process_hem_vs_gray.ipynb`, written 11:27) and its
  `tightening_process_summary.csv` appeared in the same window; they are **out of scope here** and
  were not read — but a reader should know the family grew a member after this audit began.

---

## Part 4 — premises this audit inherited

Each line is a project decision I **applied rather than verified**. Where a verdict leans on one, it
is named.

| premise | source | what would falsify it | which verdicts lean on it |
|---|---|---|---|
| **ROI is the exchangeable unit** | D5; `Research Logs/2026-09-04-f5-preregistration.md` §8 | Evidence that two ROIs from the same slide/scanner/domain are more alike than two from different ones to a degree that breaks exchangeability — e.g. the two ROIs within a domain moving together far more than across domains. **Measured here:** the analysis set has **2 ROIs per tumour domain, 7 domains, 14 total**; the full `MIDOG++.json` has 44–150 images per domain. So ROI and stratum are **not** collinear, there *is* within-stratum replication (n = 2), and Step 4.1 and Step 4.2 are genuinely two different checks. | T1-2, T1-3, T2-1 (every p and CI) |
| **Exact sign-flip / t on G−1 df, not a cluster bootstrap, at small G** | the agent protocol's own rule, and F1's precedent at G = 7 | A demonstration that the bootstrap is adequately calibrated at G = 11–14. I report both the exact sign-flip and the t interval and they agree on which budgets resolve. | T1-2, T2-1 |
| **NMS radius = evaluation match radius = 7.5 µm** | D7 | A domain whose annotations pack closer than 7.5 µm. Note this premise is load-bearing for T1-1 in a specific way: the 0.046 px escape margin is measured *against* the 7.5 µm radius, and a different radius would move it. | T1-1 (the margin), C8 |
| **`self_hit_radius = 5.0 fixed px` beside an mpp-scaled match radius** | `FSConfig` and `find_and_suppress.py`'s own comment (minimum real annotation spacing 26.2 px) | A measured case where a 5 px filter deletes a legitimate neighbouring detection, or where the 28 px annulus it leaves open costs a scored true positive. **Worth a premise review**, with T1-1's correction attached: the gap did not *create* the 403 event, it only left it visible after NMS. | T1-1 (secondary) |
| **D8's template anchor is the bbox pixel centre `(x0 + x1 − 1)/2`, which is a half-pixel whenever the bbox side is even** | `D8_TEMPLATE_ANCHOR.md`, status current | A demonstration that anchoring on a half-pixel costs more than the containment guarantee it buys. **This is the premise T1-1 actually leans on**: the 0.5 px offset between anchor and peak pixel is what makes `n_near_seed > 0` reachable at all. | T1-1 (primary) |
| **Recall must travel with `n_detections`, `precision` and `coverage_frac`** | D4; D5 | A design where recall@K at fixed K is the primary read-out and pool size is fixed by construction. | T2-1 |
| **Reading burden and worst-seed behaviour are the product metrics, not median performance** | the AnnotateDx framing | A product spec that fixes K and never asks for exhaustive review — in which case `read_95` is decorative and T2-1's second half weakens to a curiosity. | T2-1 |
| **D9's `max_peaks = 100` is the production configuration** | D9, 2026-09-12, no amendments | A later amendment or a multi-seed re-run showing the cap moves precision or recall at K ≤ 30. | T2-2 |
| **Prior audits in `Research Logs/` are not independent corroboration** | written by the same hand under these same premises | — | the reconciliation in Part 3 |

**Which verdicts would change if a premise here were wrong.** If ROI is *not* the exchangeable unit
— if, say, the domain were — then G falls from 11/14 to 7 and **T1-2 and T1-3 get weaker, not
stronger**: the p-floor rises to 1/2⁶ and fewer budgets resolve, so "consistently worse" stays
unsupported either way. If the D8 anchor were rounded to the nearest pixel — or if the self-hit
filter were referenced to `tpl_xy` at the match radius instead of 5 fixed px — the 403 event would
disappear entirely and **T1-1's verdict would change from "the mechanism is wrong" to "the event is
a geometry artifact that no longer occurs"** — the notebook's attribution to the channel change
would still be wrong, but the whole paragraph would become moot. If reading
burden is *not* a product metric, **T2-1 drops from Tier 2 to a footnote** and the notebook's
one-sided adoption recommendation becomes defensible. Nothing in Part 0 depends on any of these:
the 1,204 Tier A values and the 308 Tier B values are arithmetic and pixels, and they reproduce
under any premise.
