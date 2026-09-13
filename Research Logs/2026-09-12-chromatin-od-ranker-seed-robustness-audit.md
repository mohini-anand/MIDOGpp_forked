# Audit of `chromatin_od_ranker_seed_robustness.ipynb`: every statistic reproduces and the precision effect survives the correct tests — but the stability p-values come from a test that assumes the two arms are independent when they share a pool, and at K=20 the claim does not survive the fix; the timing readout names the wrong causal variable; and the one diagnostic the notebook demoted was pointing at something real that the root-cause paragraph gets backwards

**Scope.** Target: `threshold_maxpeaks_ablation/chromatin_od_ranker_seed_robustness.ipynb`
(untracked, executed in place 2026-09-12 18:01). Artifacts:
`threshold_maxpeaks_ablation/chromatin_od_ranker_seed_robustness_{timing,precision,per_seed_pooled,stability}.csv`
(all four untracked), the committed seed-0 oracle
`threshold_maxpeaks_ablation/chromatin_od_ranker_{timing,precision}.csv` (`00f0518`),
`databases/MIDOG++.json`, and the `images/extra_valid/*.tiff` pixels for two ROIs.
Cited sources chased: `results/premise_review_chromatin_ranker_{multiseed,seed_variance,seed_transport,d5_reproduction}.csv`,
`threshold_maxpeaks_ablation/chromatin_od_ranker_max_peaks50.ipynb`,
`Research Logs/2026-09-12-chromatin-od-ranker-variant-audit.md`. Decisions read:
`DECISIONS.md` D1 (moved to `DECISIONS_UNVERIFIED.md` 2026-09-12), D3, D4, D5 + its 2026-09-08
amendment, D6, D7, D9; `D8_TEMPLATE_ANCHOR.md`. Audit script:
`chromatin_od_ranker_seed_robustness_audit.py`; tables
`results/chromatin_od_ranker_seed_robustness_audit_*.csv`.
**Everything below is re-derived from those artifacts — and, for two (seed, ROI) draws, from the
ROI pixels themselves — never from the notebook's printed output.**

---

## Conflict of interest

`git log -1 --format='%an %ar' -- threshold_maxpeaks_ablation/chromatin_od_ranker_seed_robustness.ipynb`
returns **nothing**: the notebook has never been committed. `git status --porcelain` lists it and
all four of its CSVs as `??`. The surrounding commits (`00f0518`, `c472de9`) are `mohini-anand`
from earlier the same day, and the sibling notebook's audit log and script — same hand, same day —
are themselves untracked in the tree. So the work under audit and this audit are almost certainly
the same hand, hours apart, and only fresh context separates them. That is the signal I saw.

What limits it: every table here is produced by `chromatin_od_ranker_seed_robustness_audit.py`
from the raw artifacts, never from the notebook's output cells; §6 re-implements peak extraction,
NMS, self-hit removal, `od51` and the greedy one-to-one matcher in plain numpy/cv2 and never calls
`midog_utils.compare`, `midog_utils.evaluate` or `midog_utils.nms` for any number it reports; the
strongest single leg is the **committed** seed-0 oracle at `00f0518`, which predates this notebook.

What it cannot cover: the design premises in Part 4 — ROI as the exchangeable unit, 7.5 µm as both
NMS and match radius, precision@K on one click per ROI as the read-out, D5's own choice of bar.
Those go to `premise-reviewer` and to a reader who did not make them, not to me.

---

## Part 0 — what reproduces

The engineering is sound. **The published statistics are the statistics the code in the working
tree computes, and — for the two draws I re-ran from pixels — the statistics the pixels support.**
Every finding in Part 1 is about inference, attribution or reporting, not about the run.

| gate / check | tier | result |
|---|---|---|
| Execution coherence | — | 14 code cells, `execution_count` **1…14 contiguous**, 0 error outputs, 0 unrun cells, no unrun last cell, 1 embedded PNG. **Not vendored** (imports `midog_utils`, reads `../images/extra_valid` and `../databases/MIDOG++.json`, names this dataset's ROI files). |
| Provenance | — | `midog_utils/` is **clean** in `git status`, every module committed, newest module mtime `seed_selection.py` 2026-09-10 13:27 — two days before the artifacts (2026-09-12 18:01:04) and the notebook (18:01:07). No `.partial` read anywhere. **But the notebook and all four CSVs are untracked and have never been committed**, so they carry no commit provenance. See the labelling note below. |
| Composition | A | **20 checks, 20 pass.** 70 timing rows = 5 seeds × 14 ROIs; 420 precision rows = 5 × 14 × 2 arms × 3 budgets; 30 pooled rows; 3 stability rows; complete (ROI, seed, arm, budget) grid at 420/420; 0 duplicate keys; 0 NaN; `budget_delivered == budget` on all 420; `n_peaks == 100` on 70/70; **7 tumour types × exactly 2 ROIs each**, measured from `databases/MIDOG++.json`. |
| **Claim 1 — seed 0 vs the committed oracle** | A | **742 values compared, 0 divergences** — every one of the 11 deterministic columns shared with `chromatin_od_ranker_timing.csv` (× 14 ROIs) and all 7 shared columns of `chromatin_od_ranker_precision.csv` (× 84 rows), joined with a full outer merge that returned `both: 84, left_only: 0, right_only: 0`. The notebook's own check covered 154 of these and printed "84". Finding 8. |
| Ground truth, re-derived from `databases/MIDOG++.json` | A | **70/70 draws**: `n_gt_mitotic` equals the JSON's mitotic count for that ROI **minus the excluded seed annotation**; all 70 seed annotations are category-1 (mitotic) and unanimous-pool members. |
| Internal arithmetic of the precision CSV | A | **840 values**: `precision_at_budget == tp/budget_delivered` to 1.1e-16 and `recall_at_budget == tp/n_gt_mitotic` to 9.7e-17 on all 420 rows. `tp_at_budget` monotone non-decreasing in K within every (ROI, seed, arm) and ≤ min(K, `n_gt_mitotic`) everywhere. |
| **Claim 3 — per-seed pooled precision and the n=5 CI** | A | **30 pooled values, 0 divergences** (max \|diff\| 5.6e-17). Mean and CI reproduce to 4 dp at all three budgets: +0.0957 [+0.0128, +0.1786] / +0.0857 [+0.0291, +0.1423] / +0.0729 [+0.0287, +0.1170]. Pooled precision **is** the unweighted ROI mean on all 30 cells, because `budget_delivered == K` everywhere. |
| **Claim 4 — Levene** | A | Reproduces **exactly**: ratios 0.211382 / 0.384161 / 0.323045 to 5.6e-17, p = 6.194e-8 / 3.067e-3 / 6.583e-4 to 3.9e-18. The *test* is wrong (Finding 4); the *arithmetic* is not. |
| **Claim 5 — timing headline** | A | mean 37.141 / median 34.081 / sd 12.881 / min 24.896 / max 105.704 ms all reproduce exactly, as do all five per-seed means (41.99 / 35.22 / 34.31 / 34.27 / 39.92 ms) and percentages. |
| **Claim 2 — cap binds** | A + B | 70/70 rows carry `n_peaks == 100`. That is the run's own assertion written back to disk, so it is consistency, not independence — **the independent leg is Tier B**: recomputing extraction with `max_peaks = 2_000_000` gives **45,816** pre-cap peaks on 301.tiff/seed1 and **56,733** on 013.tiff/seed0 against a cap of 100. The pre-cap peak list is 458–567× the length of the cap, so the cap genuinely binds. |
| **Figure (cell 25), step 1a** | A | All **30 plotted values** re-derived from `_precision.csv` and compared against what the plotting cell passes in: max abs diff 5.6e-17 across all six series. |
| **Figure (cell 25), step 2** | — | PNG decoded and read. Legend maps `tm_score`→grey and `chromatin_od`→blue, matching the two `ax.plot` calls; shared y-axis spans the full 0.4286–0.6929 data range and clips nothing (including seed 3's crossover at K=10); 30 points for 30 (budget, arm, seed) cells, no over-plotting; the seed-0 shaded column is present in all three panels. Both of cell 24's claims hold on the render: `tm_score` visibly wobbles more, and seed 0 is the **minimum** of the `tm_score` line in all three panels. |
| **Tier B — from pixels** | B | `013.tiff`/seed 0 recomputed end-to-end with my own extraction, NMS, self-hit removal, `od51` and greedy one-to-one matcher: **6/6 `tp_at_budget` values** match the CSV *and* the committed oracle, both arms, all three budgets; `seed_ann_id`, `base_size`, `n_detections`, `map_median` exact. `301.tiff`/seed 1 likewise on `base_size`, `seed_ann_id`, `n_detections`, `n_retries`, `map_median`, `mad_scale`. |
| Config, three-way | — | Notebook cell 1 vs `FSConfig` defaults vs `DECISIONS.md`: **no stale override.** `NMS_RADIUS_UM = MATCH_RADIUS_UM = ev.MIDOG_RADIUS_UM = 7.5` (D7 — *not* the 5.0 µm that propagated elsewhere, and read from the module constant rather than hardcoded); `CHANNEL='hematoxylin_od'` (D3, and the unclipped OD `chromatin_density`'s docstring requires); `METHOD=cv2.TM_CCOEFF` passed **explicitly** to `tm.fused_response`, so `FSConfig.tm_method`'s stale default of 5 (`TM_CCOEFF_NORMED`) is never consulted; `MAX_PEAKS=100` and `DEEP_FLOOR_Z=-1.5` per D9; 14 ROIs of `images/extra_valid` per D5's standing constraint. |
| Gates not applicable | — | **Tier C: not run** — no divergence required it, and the sweep is ~6 min of CPU I had no reason to spend. **Leave-one-ROI-out (mode 10): no operating point is chosen inside this notebook** — `MAX_PEAKS`, the budgets and both ranking keys are fixed a priori — so there is nothing to hold out. The selection happened upstream and is sized in Finding 7 and in the prior audit's Finding 5. **Recall triad (mode 3): does not apply** — the pool is 87–100 candidates, not a 12k–39k tiling, and recall@30 spans 0.000–0.632, nowhere near saturation (Finding 14 covers the dropped `coverage_frac`). **Pre-registration adherence (mode 6) and decision-table placement (mode 7): no pre-registration exists** for this notebook; `ls "Research Logs"/*preregistration*` returns F1, F2, F4, F5, F6 only, none of which covers it. |
| Claims with nothing persisted underneath | — | **None.** Every printed number traces to one of the four CSVs, and every CSV column traces to `databases/MIDOG++.json` or the ROI pixels. |

**Provenance label.** Because the four target CSVs are untracked, Tier A against them is strictly
*consistency with an uncommitted artifact*. Two things upgrade it: the seed-0 slice is checked
against a **committed** oracle (`00f0518`, 742/742), and two draws are re-derived from pixels with
an independently written matcher. I therefore label Claims 1, 2 and the seed-0 half of Claim 3
**independent reproductions**, and the seeds 1–4 half of Claims 3, 4 and 5 **consistency with an
uncommitted artifact whose seed-0 slice is independently verified**.

---

## Part 1 — findings

### Tier 1

**Finding 1 — the closing readout's timing *mechanism* is false. The overhead has no measurable
relationship to `n_detections`; it tracks ROI area.**

Cell 27 prints: *"the overhead depends on pool composition (`n_detections`), not on which click was
drawn."* Over all 70 draws:

| stage | vs `n_detections` | vs `roi_pixels` | vs run order |
|---|---|---|---|
| total overhead `t6a+t6b+t6c` | ρ = **−0.038**, p = 0.756 | ρ = **+0.696**, p < 1e-4 | ρ = +0.349, p = 0.003 |
| `t6a` pad | ρ = −0.040, p = 0.742 | ρ = **+0.728**, p < 1e-4 | ρ = +0.329, p = 0.006 |
| `t6b` od51 loop | ρ = +0.137, p = 0.257 | ρ = +0.169, p = 0.163 | ρ = +0.470, p < 1e-4 |
| `t6c` rank | ρ = +0.191, p = 0.114 | ρ = +0.186, p = 0.124 | ρ = +0.475, p < 1e-4 |

Spearman, n = 70. `n_detections` spans 87–100 (sd 2.28) at this cap, so the variable the readout
names is both uncorrelated with the cost and nearly constant — the hypothesis is not merely
unsupported, it is close to untestable at `max_peaks=100`. The variable that *does* explain the
cost is ROI area, which is exactly what `cv2.copyMakeBorder` copies. **This confirms and extends
`Research Logs/2026-09-12-chromatin-od-ranker-variant-audit.md` Finding 2** (which established the
area dependence on 14 single-seed draws and re-timed it at 0.78–0.85 ms/Mpx) from 14 draws to 70.

The *null* half of the claim survives: the overhead is not seed-dependent, and nothing here
suggests otherwise. It is the positive attribution that is wrong, and it is wrong in a way that
would mislead anyone sizing this cost for a different pool depth — the pad is 84.8-88.8 % of the
overhead and is independent of pool size, so a deeper pool would *not* scale it.

**Correction to write.** Replace *"depends on pool composition (`n_detections`)"* with *"is
dominated by the full-ROI `cv2.copyMakeBorder` in `t6a` (84.8-88.8 % of the overhead), which tracks ROI
area (Spearman ρ = +0.73 against `roi_pixels`, n = 70) and is independent of pool size; only
`t6b`+`t6c` (~5 ms) scale with `n_detections`, and at 87–100 candidates that range is too narrow
to resolve."*

---

**Finding 2 — "n = 5 seeds … clears that bar [D5's]" is overstated: D5's bar has four clauses and
the notebook computes one. The conclusion survives three of the remaining three when I compute
them — and fails the fourth on budget depth.**

D5's *What would change my mind* (`DECISIONS.md`, and unamended — the only amendment to D5 is the
2026-09-08 `od_contrast` note, which does not touch the bar):

| D5 requires | notebook does | my recomputation |
|---|---|---|
| 5 seeds on the 14 ROIs of `extra_valid` | **yes** | confirmed |
| read at **recall@K** per D4 | precision@K only | computed below — **passes** |
| paired Δ **clustered at the ROI**, CI excluding 0 | clustered at the **seed** | computed below — **passes** |
| **majority of ROIs positive** | not reported | 11/14, 11/14, 12/14 — **passes** |
| (D5's own null sits at recall@250, and it states the advantage is "≈0 at budgets 25–100") | K = 10/20/30 | **budget-depth mismatch stands** |

The point estimate is identical under either unit — pooled precision is the unweighted ROI mean on
all 30 cells, so averaging over seeds and over ROIs commutes. Only the interval and the p move:

| metric | K | unit | G | mean Δ | t CI95 | exact sign-flip p | floor | units positive |
|---|---|---|---|---|---|---|---|---|
| precision@K | 10 | seed | 5 | +0.0957 | [+0.0128, +0.1786] | 0.1250 | **0.0625** | 4/5 |
| precision@K | 10 | **ROI** | 14 | +0.0957 | **[+0.0347, +0.1567]** | **0.0081** | 1.2e-4 | **11/14** |
| precision@K | 20 | seed | 5 | +0.0857 | [+0.0291, +0.1423] | **0.0625 (at floor)** | 0.0625 | 5/5 |
| precision@K | 20 | **ROI** | 14 | +0.0857 | **[+0.0338, +0.1376]** | **0.0049** | 1.2e-4 | **11/14** |
| precision@K | 30 | seed | 5 | +0.0729 | [+0.0287, +0.1170] | **0.0625 (at floor)** | 0.0625 | 5/5 |
| precision@K | 30 | **ROI** | 14 | +0.0729 | **[+0.0311, +0.1147]** | **0.0028** | 1.2e-4 | **12/14** |
| recall@K | 10 | **ROI** | 14 | +0.0340 | [+0.0072, +0.0608] | 0.0056 | 1.2e-4 | 11/14 |
| recall@K | 20 | **ROI** | 14 | +0.0615 | [+0.0195, +0.1036] | 0.0021 | 1.2e-4 | 11/14 |
| recall@K | 30 | **ROI** | 14 | +0.0610 | [+0.0201, +0.1020] | 0.0034 | 1.2e-4 | 12/14 |

Exact two-sided sign-flip, all 2^G patterns enumerated (2⁵ = 32, 2¹⁴ = 16,384). Two things follow
that the notebook does not say:

* **At G = 5 the distribution-free floor is 2/2⁵ = 0.0625**, so a sign-flip null on 5 seeds *cannot
  produce a p below 0.0625 at all*, and at K = 20 and K = 30 the notebook's own evidence sits
  exactly on it. The t interval is doing all the work, on 4 degrees of freedom.
* **The two CIs answer different questions and neither answers both.** All five seeds share the
  same 14 ROIs, so between-ROI variance cancels out of the 5 pooled deltas: the n = 5 CI
  generalises over *clicks, holding these ROIs fixed*. The G = 14 CI generalises over *ROIs,
  averaging clicks out*. Reporting one alone as "the effect size" overstates its reach in whichever
  direction is unstated.
* **Why the G = 14 interval is the *narrower* one, which looks backwards and is not.** The two
  estimators do not resample the same objects. The seed-level SD is the spread of **five**
  quantities, each already an average over 14 ROIs; the ROI-level SD is the spread of **fourteen**
  quantities, each already an average over 5 seeds. Averaging 5 seeds per ROI removes more of the
  click noise than averaging 14 ROIs per seed removes of the ROI noise, so the ROI-level units are
  individually quieter *and* there are nearly three times as many of them. I report both rather
  than the narrower one: the honest reading is that the effect clears zero on either axis of
  generalisation taken alone, and neither interval covers both axes at once.

**Where the bar is still not cleared:** D5's null is Δrecall@250 on `f1_seed_sweep.csv` (7 ROIs ×
5 seeds, z = 1.0), and D5 says in so many words that the advantage is "≈0 at budgets 25–100 and
**significantly negative at 5,000**". I reproduced that from
`results/premise_review_chromatin_ranker_d5_reproduction.csv`: Δrecall@250 = +0.0322,
CI [−0.0474, +0.1118], p = 0.360, positive on 3/7 — D5's stated numbers exactly, and
Δprecision@5000 = −0.000294, p = 0.028. Nothing measured at K = 10–30 speaks to that. And D5's bar
names **`od31`, `od_falloff`, `mask_od_mean` and `od_contrast`** as the axes to sweep; `od51` is
the axis D5 *demoted*.

**Correction to write.** Replace *"this clears that bar … for the first time"* with *"this matches
the seed count D5's bar names, on the ROI set it names, at the production config for the first
time. Clustered at the ROI as D5 requires (G = 14, exact sign-flip), the paired Δ is
+0.096/+0.086/+0.073 precision and +0.034/+0.062/+0.061 recall, CI excluding 0 and positive on
11–12 of 14 ROIs at every budget — so the bar's clustering, metric and majority clauses are met.
It remains outside D5's own measurement depth: D5's null is Δrecall@250 and it states the
advantage is ≈0 at budgets 25–100, which K = 10–30 does not address."*

---

**Finding 3 — the `seed_annulus_empty` root cause has its central word backwards, and the check
was the only thing pointing at a real property of that draw. The demotion is nevertheless
defensible, and the candidate's measured effect on the seed-1 comparison is exactly zero.**

Cell 12 says the annulus candidate is *"a secondary local maximum ('echo' of the same physical
structure)"* of the click's own correlation surface, which *"self-hit suppression correctly leaves
… (it isn't within 5px)"*. Recomputed from pixels on `301.tiff`/seed 1, against `013.tiff`/seed 0
as the normal case:

| | 013.tiff / seed 0 | **301.tiff / seed 1** |
|---|---|---|
| `base_size` | 31 | **23** (smallest in the run) |
| peaks within the 5 px self-hit zone | **1** | **0** |
| top peak's distance from the click | **0.50 px** | **3,110 px** |
| top peak score | 3.0023 (= the fused global max) | 1.4238 (elsewhere in the ROI) |
| fused value **at the click's own box centre** | 3.0023 | **0.9050** |
| rank-100 score cut | 1.7408 | **0.9585** |
| does the click's own location make the top 100? | **yes** | **no** |
| NMS keeps / self-hit removes | 98 / 1 → 97 | **100 / 0 → 100** |

**There is no primary.** On 301/seed 1 the template's own source location scores *below* the
rank-100 cut and never enters the pool, so the surviving candidate 26.6 px away is not a secondary
maximum of anything — it is the only near-click response there is. (`TM_CCOEFF` is unnormalised,
so a low-contrast 23 px template correlates more strongly with darker structures elsewhere than
with itself; that is the mechanism, and it also explains the unremarked `n_detections = 100` on
that one draw, the only one of 70 where the post-NMS list hits exactly the cap value.)

**The demotion is still right, and here is the arithmetic.** The candidate at (2952, 193) is a
*forced* false positive: its nearest ground truth of either category is `ann_id 14742` at 21.0 px,
which **is** the seed annotation `gt_eval` excludes, and the next nearest is 74.3 px — outside the
29.61 px match radius. It therefore cannot be credited in any arm. Its ranks:

| arm | rank of the echo candidate | in the K ≤ 30 list? |
|---|---|---|
| `tm_score` (score = 0.9625) | **94 of 100** | no |
| `chromatin_od` (od51 = 0.1445) | **66 of 100** | no |

Dropping it and re-matching changes `tp_at_budget` by **0 on every budget in both arms**, so its
effect on the 301/seed 1 cell delta and on seed 1's pooled delta is **0.00000 at K = 10, 20 and
30**. The invoker's specific worry — that this candidate biases the seed-1 result — is answered
in the negative, measured, not argued.

But the notebook's *stated reason* is not the reason. The candidate is not "symmetric": its ranks
differ by 28 positions between the arms. It fails to bias the comparison because **neither rank
enters the list**, not because the two arms treat it alike. Had it scored 20 places higher it
would have entered `tm_score`'s K = 30 list and not `chromatin_od`'s, and the asymmetry would have
run in `chromatin_od`'s favour.

**What the demotion risks masking.** Nothing about leakage — a near-click candidate cannot
double-credit the seed annotation, because `gt_eval = gt[gt.ann_id != seed_ann_id]` already removes
it, so the check never had teeth as a correctness gate. What it *was* catching, uniquely, is
"this draw's template does not match its own source location strongly enough to reach the pool" —
a genuine property of a weak seed, on the draw that also needed 2 retries. No other check in the
notebook covers that. Demoting the check while writing a root cause that presupposes a primary
self-peak is precisely how the next instance gets waved through.

**Correction to write.** Replace the cell-12 paragraph with a root cause that matches the pixels:
*"On this draw the click's own location does not reach the pool at all — `base_size = 23`, and the
fused value at the template-box centre (0.9050) falls below the rank-100 cut (0.9585), so
`TM_CCOEFF` scores the unnormalised template higher against darker structures elsewhere than
against itself. The surviving candidate 26.6 px away is therefore the strongest near-click
response, not an echo of a stronger one, and `n_detections = 100` on this draw follows (NMS removed
nothing, self-hit had nothing to remove). It is a forced false positive — its only ground truth
inside `match_radius` is the excluded seed annotation — and it ranks 94/100 by `score` and 66/100
by `od51`, so it enters neither arm's K ≤ 30 list and its measured effect on the delta is 0.0000 at
every budget. The check is demoted because its failure condition maps to no correctness property
(the seed annotation is already excluded from `gt_eval`), not because this instance was harmless."*

---

### Tier 2

**Finding 4 — Levene is the wrong test twice over, and at K = 20 the correct test moves the p by
up to two orders of magnitude. The variance *ratios* are unaffected; only the evidence is.**

Two independent violations, both anti-conservative:

1. **The two samples are paired.** `tm_score` and `chromatin_od` are scored on the *same* ROI, the
   *same* seed and the *same* candidate pool — only the sort key differs. Pearson r between the
   two arms' ROI-centred residuals: **+0.327** (p = 0.006) at K = 10, **+0.564** (p < 1e-4) at
   K = 20, **+0.621** (p < 1e-4) at K = 30. `scipy.stats.levene` is a two-*independent*-sample test.
2. **The 70 residuals per arm are 14 clusters of 5 that sum to zero by construction**, so the free
   count is 14 × 4 = 56 per arm, not 70. Levene is handed n = 70 and uses it.

Replacing it with the cluster-aware, arm-paired version — per-ROI variance across the 5 seeds,
then 14 paired units:

| K | var(tm) | var(ch) | ratio | ch lower on | exact sign test | exact sign-flip on log-ratio | Wilcoxon | **notebook's Levene p** |
|---|---|---|---|---|---|---|---|---|
| 10 | 0.02636 | 0.00557 | 0.211 | **14/14** | 1.22e-4 | 1.22e-4 | 1.22e-4 | 6.19e-8 |
| 20 | 0.01511 | 0.00580 | 0.384 | **10/14** | **0.180** | **0.0254** | 0.0166 | 3.07e-3 |
| 30 | 0.01157 | 0.00374 | 0.323 | **13/14** | 1.83e-3 | 2.20e-3 | 3.66e-4 | 6.58e-4 |

The variance ratios are *identical* to the notebook's (the ratio of mean per-ROI variances equals
the ratio of pooled residual variances here), so the reported effect size is untouched. What
changes is the strength of evidence, and it changes most where the notebook's claim is weakest:
at **K = 20** the honest range is p = 0.017–0.180 depending on how conservatively you treat the
pairing, against a printed 0.0031. The direction — `chromatin_od` is the more click-stable arm —
holds at all three budgets under every replacement test.

**Caveat on my own finding.** G = 14 and the sign test discards no ties here, so 14/14 at K = 10 is
at the exact floor (1.22e-4) and means "as strong as 14 units can resolve", not a measured tail.

---

**Finding 5 — seed 3 is negative at K = 10, one whole domain is negative at K = 20 and K = 30, and
the closing readout reports neither.** (Step 4.8 — worst-seed behaviour is the product metric.)

Per-seed deltas at K = 10: `+0.1643, +0.1214, +0.1143, −0.0143, +0.0929` — **4 of 5 positive**.
Cell 21 prints this list; cell 27's readout quotes only the mean and CI.

| K | worst seed | worst-seed Δ | seeds negative | worst (ROI, seed) cell | cells negative |
|---|---|---|---|---|---|
| 10 | seed 3 | **−0.0143** | **1/5** | 245.tiff/seed 4, **−0.300** | 14/70 |
| 20 | seed 3 | +0.0393 | 0/5 | 246.tiff/seed 0, −0.200 | 14/70 |
| 30 | seed 4 | +0.0310 | 0/5 | 246.tiff/seed 3, −0.167 | 10/70 |

Clustered at the domain (G = 7, 2 ROIs each — mean over 10 draws per domain):

| domain | ΔP@10 | ΔP@20 | ΔP@30 |
|---|---|---|---|
| canine soft tissue sarcoma | +0.19 | +0.185 | +0.150 |
| canine lung cancer | +0.18 | +0.130 | +0.087 |
| human melanoma | +0.17 | +0.045 | +0.047 |
| human neuroendocrine tumor | +0.09 | +0.095 | +0.080 |
| canine cutaneous mast cell tumor | +0.02 | +0.065 | +0.120 |
| human breast cancer | +0.01 | +0.090 | +0.047 |
| **canine lymphosarcoma** | +0.01 | **−0.010** | **−0.020** |
| exact sign-flip, G = 7 | **p = 0.0156** (7/7, *at the floor*) | p = 0.0312 (6/7) | p = 0.0312 (6/7) |

**Canine lymphosarcoma is the worst domain at every budget and is negative at two of three, across
all five seeds.** This confirms and *attenuates* the prior audit's seed-0-only finding that the
two lymphosarcoma ROIs (245, 246) are the repo's only two of that domain and that the domain
reverses: at seed 0 alone the prior audit measured −0.10 (K=20) and −0.067 (K=30); averaged over
5 seeds it is −0.010 and −0.020. So it is **not** a seed-0 artifact, but it is 3–10× smaller than
one click suggested. At G = 7 the two-sided sign-flip floor is 2/2⁷ = 0.0156, so the K = 10 p sits
exactly on it and means "as strong as 7 domains can resolve".

Because the strata here hold **exactly 2 ROIs each**, Step 4.1 and Step 4.2 are genuinely two
different tests (there *is* within-stratum replication), and the domain-level intervals above are
the weaker, more honest ones for a tool meant to travel across scanners.

---

**Finding 6 — the motivating premise is confirmed in direction and contradicted in magnitude, and
the notebook reports neither comparison.** (Step 4.5.)

Cell 0's stated motivation: *"the all-5-seed mean effect is roughly 2-6x smaller than the seed-0
value alone"*, sourced to `results/f5_nms_radius_ablation.csv` via the premise review. Chased to
`results/premise_review_chromatin_ranker_multiseed.csv`, that citation is **faithful**: seed-0 ÷
all-seed ratios are 6.67× (K=10), 2.24× (K=25) and 3.83× (K=50) at z = 0.5, stable across all six z
levels, with `seed0_rank_of_n = 1` throughout.

This notebook's own measurement at the production config:

| K | seed-0 Δ | all-5-seed mean Δ | ratio | seed-0 rank of 5 |
|---|---|---|---|---|
| 10 | +0.1643 | +0.0957 | **1.72×** | 1 |
| 20 | +0.1429 | +0.0857 | **1.67×** | 1 |
| 30 | +0.1071 | +0.0729 | **1.47×** | 1 |

So the premise's *direction* is reproduced exactly — seed 0 really is the best of 5 at every
budget, as at the f5 config — while its *magnitude* is roughly halved: 1.5–1.7× inflation, not
2–6×. The notebook prints "seed0 rank among 5: 1" three times and never puts the ratio beside the
number it was commissioned to correct.

There is a second, larger result hiding in the same comparison that the notebook also does not
claim: at the f5 config the all-seed effect is **not distinguishable from zero**
(`allseed_p` = 0.74–0.81 at K = 10, 0.32–0.40 at K = 25), whereas at this config it is
(t p = 0.033 at n = 5 seeds, 0.0048 at G = 14 ROIs). The notebook's stated reason for refusing to
take f5 as a stand-in is therefore vindicated by its own data, and it never says so.

---

### Tier 3

**Finding 7 — Verification A's printed count is wrong, and its markdown claims more than its code
checks.** Cell 9 prints *"84 values compared, 0 divergences"*. It actually compared 5 timing
columns × 14 ROIs = **70** plus 84 `tp_at_budget` values = **154**; the literal `14 * 2 * 3` in
the f-string counts only the second leg. Separately, the cell-8 markdown says *"**every column**
for both arms … must match … bit-for-bit"* while the code compares **5 of 11** deterministic timing
columns and **1 of 7** shared precision columns. It omits `tumor_type`, `n_retries`, `n_peaks`,
`roi_pixels`, `chromatin_od_nan_rate`, `chromatin_od_largest_tie_block`, `n_gt_mitotic`,
`budget_delivered`, `precision_at_budget` and `recall_at_budget`. I compared all 18 → **742 values,
0 divergences**, so the *claim* is true; the check that printed it was narrower than the claim.

**Finding 8 — 280 of the 350 reported checks cannot fail, and after the demotion no check in the
`CHECKS` frame can report a failure at all.** `cp.Arm` is constructed without `caps=`, so
`invariants.check_no_cap` receives `caps=()`, `hits` is always empty, and it can only return
`passed=True` — 140 of those (2 arms × 70 draws). `invariants.check_nms_radius` **raises** on
failure rather than returning `passed=False` — 140 more that can only ever appear as passing.
After `seed_annulus_empty` is moved out of `FATAL_CHECKS`, the `fatal_bad` filter is non-empty by
construction impossible, and `assert len(fatal_bad) == 0` is a tautology. So *"checks passed:
349/350"* is 280 tautologies plus 70 genuine diagnostics. The real gates are the bare
`assert n_peaks == MAX_PEAKS` in `run_roi` and `check_nms_radius`'s raise.

**This is an improvement on the sibling notebook, and materially so.** The prior audit's Finding 6
flagged the sibling for passing `caps=(MAX_PEAKS,)`, which makes `check_no_cap` assert that the
*post-NMS* length is not 100 when the cap is applied *pre-NMS*. On `301.tiff`/seed 1 the post-NMS
`n_detections` **is exactly 100** — so the sibling's configuration would have raised
`InvariantError` on that draw and killed the run. Cell 4's markdown reasons this out correctly and
the change is right; only the pass-rate headline is inflated.

**Finding 9 — `ds.image_annotations(annotations, fn_)` in cell 13 leaves `category_id` at its
permissive `None` default**, so the line printed as *"nearest ground truth to that candidate"* is
nearest annotation of *either* category. Harmless here — the nearest is category 1 and is the seed
— and pooling categories is what `bucket_detections` does too, so the *matching* is right. But the
label is wrong, and this is the same defect class as the repo's largest prior audit finding
(`Research Logs/2026-09-08-tp-fp-separability-audit.md` §1). One word: *"nearest annotation (either
category)"*.

**Finding 10 — the "1.09x-3.5x" pad figure is attributed to the wrong document, and the "several-fold"
characterisation overstates the notebook's own data.** Cell 19 credits
`chromatin_od_ranker_max_peaks50.ipynb` with *"a single `cv2.copyMakeBorder` call can run 1.09x-3.5x
its own median"*. That notebook's **own** measurement is `chromatin_od: mean=1.73x max=3.50x` — a
max/**min** ratio over 5 repeats of the whole stage-6 branch, not of one pad call against its
median. The 1.09×–3.51× figure is the *pad call* from
`Research Logs/2026-09-12-chromatin-od-ranker-variant-audit.md` Finding 3, which `max_peaks50`
quotes in three places and does not produce. The chain is sound; the citation skips a link.

Separately, cell 19 calls the per-seed spread *"a several-fold swing in the MEAN"* while the same
cell prints the spread as **1.23×**.

**Finding 11 — the pad-spike attribution explains seed 0's excess but not the residual spread,
which tracks run order.** Decomposed:

| seed | baseline pipeline | overhead mean | overhead **median** | `t6a` pad | `t6b` loop | `t6c` rank |
|---|---|---|---|---|---|---|
| 0 | 1239.4 ms | **41.99** | 29.46 | 37.27 | 4.39 | 0.33 |
| 1 | 1540.5 | 35.22 | 34.06 | 30.01 | 4.85 | 0.36 |
| 2 | 1518.4 | 34.31 | 33.76 | 29.12 | 4.80 | 0.39 |
| 3 | 1530.9 | 34.27 | 33.56 | 29.11 | 4.81 | 0.35 |
| 4 | 1763.8 | 39.92 | 38.77 | 33.85 | 5.65 | 0.42 |

Exactly two draws carry pad spikes, both in seed 0: 013.tiff at **100.5 ms** (3.53× the 28.46 ms
median) and 403.tiff at **96.8 ms** (3.40×); the third-largest is 1.57×. Dropping those two moves
seed 0's mean from 41.99 to **31.74 ms** — *below every other seed* — and the spread across the
five seed means goes **up**, from 1.23× to 1.26×. So the pad-noise attribution is correct for the
seed-0 outlier and does **not** account for the rest: seed 0's *median* overhead (29.46 ms) is the
lowest of the five and seed 4's (38.77 ms) the highest, in the same order as the baseline pipeline
time, which drifts **1239 → 1764 ms (1.42×)** across the run. Seeds ran in the outer loop, so
**seed index is run order**, and per-seed timing differences are confounded with machine drift
(Spearman ρ = +0.349 of overhead on order, p = 0.003).

The notebook's *conclusion* — the overhead is not seed-dependent — survives all of this, and is in
fact better supported by the baseline drift than by the pad-noise argument it offers. Its
explanation is what needs replacing.

**My own re-timing, and its limit.** 40 tight repeats of `cv2.copyMakeBorder(hem, 25,25,25,25)` on
013.tiff (5412×7215) give median **33.03 ms**, min 28.63, max 37.26 — **max/median 1.13×**, no
spike. That is a *warm-allocator* loop and under-samples the cold first-touch cost; it is
consistent with the prior audit's allocator explanation (the spike is not a property of the ROI)
and is **not** evidence against the 3.5× figure. Stated as a limitation of my measurement, not a
finding about the notebook.

**Finding 12 — cell 19's "seed_index=0 mean overhead here: 41.988ms — matches
`chromatin_od_ranker_timing.csv`" is not true, and its justification is a category error.**

| quantity, mean over the same 14 ROIs | committed oracle | this run, seed 0 | diff |
|---|---|---|---|
| `t_chromatin_od_overhead_ms` | **45.094** | **41.988** | 3.106 ms (**6.9 %**) |
| `t6a_od_pad_ms` | 40.095 | 37.272 | 2.823 ms (7.0 %) |
| `t_baseline_pipeline_ms` | 1371.295 | 1239.392 | 131.9 ms (9.6 %) |

Per-ROI overhead is not identical either — max abs difference **11.336 ms**. Nor should it be:
these are wall-clock measurements from two separate single-shot runs. The parenthetical that
justifies the claim — *"Verification A above already confirmed the per-ROI values this is computed
from are identical"* — cites a determinism check on `seed_ann_id`, `base_size`, `n_detections`,
`map_median` and `mad_scale`, **none of which is a timing column**. Verification A deliberately and
correctly compared only deterministic columns; it says nothing about wall clock, and cannot.

The two runs agreeing to 7 % is a perfectly good result for a single-shot timing comparison and is
worth saying. What is not supportable is calling it a match and citing a determinism check for it.

**Correction to write.** *"seed_index=0 mean overhead here: 41.99 ms against 45.09 ms in
`chromatin_od_ranker_timing.csv` — 6.9 % apart, which is within the single-shot pad noise both runs
carry (per-ROI differences reach 11.3 ms). Verification A compared deterministic columns only and
does not speak to wall clock."*

**Finding 13 — `coverage_frac` is computed by `evaluate_arms` and dropped by `PRECISION_LONG`'s
column selection before `to_csv`.** The recall triad is otherwise complete — `n_detections`,
`precision_at_budget` and `recall_at_budget` all travel with recall on all 420 rows. Harmless and
favourable: the pool is 87–100 candidates and recall@30 spans 0.000–0.632, so this is nowhere near
the saturated-pool regime D4/D5 guard against. (Repeats the prior audit's Finding 7 for this
notebook.)

**Finding 14 — 6 of the 70 draws are duplicate clicks, and the notebook's own guard is set two
units too loose.** There are **64** distinct (ROI, click) pairs across 70 draws:

| ROI | `seed_ann_id` | seeds that drew it |
|---|---|---|
| 233.tiff | 5761 | **0, 1, 4** |
| 201.tiff | 4457 | 0, 1 |
| 201.tiff | 4475 | 2, 3 |
| 094.tiff | 2494 | 2, 4 |
| 529.tiff | 24918 | 0, 4 |

Duplicated draws contribute identical values, which *shrinks* the observed within-ROI seed
variance — anti-conservative for both the n = 5 CI and the Levene/paired-variance comparison.
Restricting to the 64 unique pairs raises the mean per-ROI variance by 7 % for `tm_score` and 10 %
for `chromatin_od` at K = 10 (0.02824 vs 0.02636; 0.00614 vs 0.00557), so the ratio finding is
unaffected. The notebook asserts `n_distinct >= 2`, which 233.tiff (3 distinct) and 201.tiff
(3 distinct) clear with room to spare; the honest effective-unit count for the seed axis is **64,
not 70**.

**Finding 15 — the retry loop, checked.** `n_retries` is 0 on 59 draws, 1 on 10, 2 on 1. The loop
drops a rejected row from `working` and redraws on the same RNG stream, so it biases the draw
toward seeds that pass `ss.tightened_template_box` — which is what the gate is for, not a bias
toward "easy" seeds in the detection sense. Worth noting that the single 2-retry draw is
**301.tiff/seed 1**, the same draw whose template then failed to reach its own pool (Finding 3).

---

## Part 2 — verdict per conclusion

| # | Conclusion (cell) | Verdict | Evidence |
|---|---|---|---|
| N1 | Config is identical to `chromatin_od_ranker_variant.ipynb` except `SEED_INDICES` (cells 0, 1, 4) | **reproduces** | Three-way config check vs `FSConfig` and `DECISIONS.md`; `TM_CCOEFF` and `hematoxylin_od` passed explicitly so the stale dataclass defaults are inert; NMS radius read from `ev.MIDOG_RADIUS_UM`, no 5.0 µm override. |
| N2 | **Claim 1** — seed 0 reproduces the committed oracle "exactly … bit-for-bit"; "84 values compared, 0 divergences" (cell 9) | **reproduces** | **742 values, 0 divergences** across all 18 shared deterministic columns. The claim is true; its printed count (84) omits the 70-value timing leg, and the code checked 6 of 18 columns. Finding 7. |
| N3 | **Claim 2** — `MAX_PEAKS=100` binds on all 70 draws (cell 11) | **reproduces** | 70/70 in the CSV, plus Tier B: 45,816 and 56,733 pre-cap peaks on two recomputed draws against a cap of 100. |
| N4 | `chromatin_od` is NaN-free on every draw; seeds draw materially different clicks (cell 11) | **reproduces, overstated** | NaN-free 70/70, confirmed. "Materially different" rests on `min distinct = 3/5`; there are **64 distinct (ROI, click) pairs of 70**, with one ROI drawing the same click on 3 of 5 seeds. Finding 14. |
| N5 | **Claim 6** — the annulus candidate is a *secondary* "echo" of the click's own correlation surface, self-hit correctly leaves it (cell 12) | **does not reproduce** | On 301/seed 1 there is **no** peak within 5 px of the click, and the fused value at the click's box centre (0.9050) is **below** the rank-100 cut (0.9585): the click's own location never entered the pool. There is no primary. Contrast 013/seed 0, where the click *is* the global max at 0.50 px. Finding 3. |
| N6 | **Claim 6** — it is not a pipeline bug; the extra candidate is a false positive competing for a rank slot (cell 12) | **reproduces** | Forced FP, verified: only GT inside `match_radius` is the excluded seed annotation at 21.0 px; next is 74.3 px. Not a correctness defect. |
| N7 | **Claim 6** — it is "symmetric … not something that biases one arm over the other" (cell 12) | **reproduces, overstated** | The *conclusion* is exactly right and I measured it: effect on the delta is **0.00000 at K = 10, 20 and 30**. The *reason* is not symmetry — it ranks **94/100** by `score` and **66/100** by `od51`, a 28-position asymmetry. It is harmless because neither rank enters a K ≤ 30 list. Finding 3. |
| N8 | **Claim 6** — demoting `seed_annulus_empty` to informational while keeping `no_cap`/`nms_radius` fatal (cell 7) | **reproduces, overstated** | The demotion is defensible: the check's failure condition maps to no correctness property, since `gt_eval` already excludes the seed annotation. But `no_cap` receives `caps=()` and can only pass, and `check_nms_radius` raises rather than returning `passed=False`, so **no check left in `CHECKS` can report a failure** and "349/350" is 280 tautologies + 70 diagnostics. Finding 8. |
| N9 | **Claim 3** — per-seed deltas, mean and 95 % CI at K = 10/20/30 (cells 21, 27) | **reproduces** | 30 pooled values, 0 divergences; all three means and CIs match to 4 dp. |
| N10 | **Claim 3** — n = 5 seeds with a t-based CI "clears D5's bar … at the actual production config for the first time" (cells 0, 27) | **reproduces, overstated** | D5's bar has four clauses; the notebook computes one. Clustered at the ROI as D5 requires (G = 14, exact sign-flip): Δ = +0.096/+0.086/+0.073 precision, +0.034/+0.062/+0.061 recall, CI excluding 0, positive on 11–12 of 14 — so three more clauses pass under my recomputation. The budget-depth clause does not: D5's null is Δrecall@250 and it states the advantage is ≈0 at K = 25–100. Also: at G = 5 the sign-flip floor is 0.0625, and K = 20/30 sit exactly on it. Finding 2. |
| N11 | **Claim 4** — variance ratios 0.211/0.384/0.323 with Levene p = 6.2e-8/3.1e-3/6.6e-4 (cells 23, 27) | **reproduces, overstated** | Ratios and p-values reproduce to 1e-17. The **ratios stand**; the **p-values do not** — Levene assumes two independent samples and the arms are paired (residual r = +0.33/+0.56/+0.62), and 70 residuals per arm carry 56 free values. Cluster-aware paired replacement: K=10 p = 1.2e-4 (14/14), K=20 **p = 0.017–0.180** vs 0.0031 printed, K=30 p = 3.7e-4–2.2e-3. Direction holds everywhere. Finding 4. |
| N12 | **Claim 5** — the chromatin_od overhead is stable across seeds, 34.27–41.99 ms, "matching the single-seed report" (cells 19, 27) | **reproduces** | Every timing statistic reproduces exactly; seed-0 mean 41.988 ms matches the committed oracle. The null claim survives: the baseline pipeline itself drifts 1.42× across the run in seed order, so the 1.23× overhead spread is well inside run-to-run variation. |
| N12b | **Claim 5** — "seed_index=0 mean overhead here: 41.988ms — matches `chromatin_od_ranker_timing.csv`" (cell 19) | **does not reproduce** | The oracle's mean is **45.094 ms**, 6.9 % higher; per-ROI overhead differs by up to 11.3 ms. Both are single-shot wall clock from separate runs, so they were never going to be equal — and the cited justification (Verification A) compared deterministic columns only. Finding 12. |
| N13 | **Claim 5** — the spread is single-shot `copyMakeBorder` pad noise already on file, not a new seed effect (cell 19) | **reproduces, overstated** | Correct for the two 3.4–3.5× spikes (013 and 403, both in seed 0). But removing them drops seed 0 to 31.74 ms — *below all others* — and the spread rises to 1.26×; the residual structure follows **run order**, not pad noise. Also "several-fold swing in the MEAN" describes an observed 1.23×, and the 1.09×–3.5× figure is credited to `max_peaks50` when it originates in the variant audit. Findings 10, 11. |
| N14 | **Claim 5** — "the overhead depends on pool composition (`n_detections`), not on which click was drawn" (cell 27) | **does not reproduce** | ρ(overhead, `n_detections`) = **−0.038, p = 0.756**, over 70 draws, with `n_detections` spanning only 87–100. The cost tracks `roi_pixels` (ρ = +0.696, p < 1e-4; ρ = +0.728 for the pad alone). Finding 1. |
| N15 | **Motivation** — seed 0 is the best of 5 and f5 cannot answer the effect size at this config (cell 0) | **reproduces** | `seed0_rank_of_n = 1` in the f5 premise-review artifact at every z and K; seed-0 rank 1 of 5 here at every budget. f5's config divergences (`MAX_PEAKS = 2_000_000`, no `tightened_template_box`) confirmed from `f5_nms_radius_ablation.py`. Vindicated further: the f5 all-seed effect is p = 0.74–0.81 at K = 10, this one is p = 0.0048–0.033. |
| N16 | **Motivation** — "the all-5-seed mean effect is roughly 2-6x smaller than the seed-0 value alone" (cell 0) | **reproduces**, but the notebook's own data contradicts its magnitude | The citation is faithful to `premise_review_chromatin_ranker_multiseed.csv` (2.24×–6.67×). At the production config this notebook measures **1.47×–1.72×**, and never says so. Finding 6. |
| N17 | Figure (cell 25) — `tm_score` wobbles more; seed 0 is not a typical point (cell 24) | **reproduces** | 30 plotted values re-derived, max diff 5.6e-17; render read — legend colours correct, shared y-axis clips nothing, seed 0 is the minimum of the `tm_score` line in all three panels. |
| N18 | `budget_delivered == budget` on all rows (cell 17) | **reproduces** | 420/420. |
| N19 | The closing caveat — "still 5 seeds, not an exhaustive sweep"; `max_peaks=100` unchanged; "nothing in this notebook proposes a new production value" (cells 0, 27) | **reproduces, overstated** | The hedging on sample size is genuine, repeated and accurate, and the D9 statements are exactly right — `MAX_PEAKS=100` is unchanged and no new default is proposed. What the caveat does **not** hedge is the two results that would matter to a reader deciding whether to adopt: `chromatin_od` **loses on 1 of 5 seeds at K = 10** (−0.0143), and **canine lymphosarcoma — one of the seven domains, and both of its ROIs — is negative at K = 20 and K = 30 across all five seeds.** Both are in the notebook's own data and neither reaches the readout. (The separate "clears D5's bar" clause is ruled on at N10.) Finding 5. |

---

## Part 3 — what was re-run versus read

**Audit script.** `chromatin_od_ranker_seed_robustness_audit.py` at the repo root. It **runs start
to finish under `/Users/mohinianand/anaconda3/bin/python3` in 36 s** (verified as a whole, three
times, most recently after the last edit) and regenerates every table in this log into
`results/chromatin_od_ranker_seed_robustness_audit_*.csv` (29 files). Nothing in this log was
computed in a throwaway heredoc.

**Tier A (re-derived from persisted artifacts, written in plain numpy/pandas/scipy).** Execution,
provenance and composition gates; the full 742-value column-by-column seed-0 oracle comparison;
ground truth re-derived from `databases/MIDOG++.json` for all 70 draws; 840 internal-arithmetic
checks on the precision CSV; all 30 per-seed pooled precisions; the n = 5 seed-level and G = 14
ROI-level effect sizes at precision@K *and* recall@K, with exact sign-flip enumerated over all
2⁵ = 32 and 2¹⁴ = 16,384 patterns; the G = 7 domain-level version; Levene reproduced and replaced;
the full timing decomposition and Spearman correlation matrix; duplicate-click accounting; the
30 plotted figure values.

**Tier B (~40 s of compute, and what I spent it on).** Two full single-seed pixel recomputes —
`301.tiff`/seed 1 (the draw the demoted check fired on, chosen because the notebook's entire
root-cause argument rests on it and nothing about it is persisted) and `013.tiff`/seed 0 (chosen
because it is the only draw with a *committed* oracle, so it converts one Tier A consistency check
into an independent reproduction). Both were run twice: once capped at 100, once at 2,000,000 to
count the pre-cap peaks. For 013/seed 0 I re-implemented extraction, NMS, self-hit removal, `od51`
and the greedy one-to-one matcher myself and reproduced all 6 `tp_at_budget` values. Plus a
40-repeat `cv2.copyMakeBorder` microbenchmark. Reading the four CSVs and the annotation database
is Tier A and is not charged here.

**Tier C — not run.** No sweep was re-executed; the measurements that touch pixels are §6 of the
script (`claim6`) and §10 (`check_pad_bench`). Nothing in Tier A or B produced a divergence I could
not explain by reading, so full re-execution was not warranted. The notebook took 351 s; had I run
it I would have raised `ExecutePreprocessor.timeout` to 21600 and executed a copy at
`threshold_maxpeaks_ablation/.audit_tmp_seed_robustness.ipynb`, never `--inplace` on the target.

**What I looked for beyond Step 4's named modes (mode 12).** For each headline claim I asked what
would have to be true for it to be false, and went looking:

* *Could the seed-0 "exact" reproduction be exact on the columns checked and wrong elsewhere?*
  Compared all 18 shared deterministic columns rather than the notebook's 6. It is not.
* *Could the pool have been reordered between the two arms' evaluations?* `pool['od51'] = [...]`
  adds a column; both `sort_values` calls return new frames used only for timing; `compare._rank`
  takes a fresh stable-mergesort copy. Both arms receive the same object via a default-arg lambda.
  `chromatin_od_largest_tie_block == 1` on 70/70, so there are no ties for sort stability to decide.
* *Could the seed annotation leak back in as a true positive?* `gt_eval` excludes it by `ann_id`,
  and I re-derived `n_gt_mitotic` = JSON mitotic count − 1 on all 70 draws independently.
* *Could the cap value collide with a post-NMS list length and go unnoticed?* It does, on exactly
  one draw (301/seed 1, `n_detections == 100`) — which is why dropping `caps=` mattered (Finding 8).
* *Could the click's own peak be double-counted or mis-suppressed?* Checked from pixels on two
  draws; found the opposite problem on one of them (Finding 3).
* *Could the 5 seeds be less independent than claimed?* 6 of 70 draws are exact repeats
  (Finding 13).
* *Could the per-seed timing differences be a run-order artifact?* Yes, and they are (Finding 11).
* *Could the figure plot a different quantity from the table?* Re-derived all 30 values; it does not.

**Facts in the agent file's appendix that no longer hold** (it is dated 2026-09-12, today, and most
of it holds):

* **"77 of 226 `results/*.csv` untracked … two days later all 233 were tracked."** Neither count is
  current: `results/` now holds **267** entries, and the most recent notebook family's outputs
  (`threshold_maxpeaks_ablation/*`, `premise_review_chromatin_ranker_*`) are untracked again. The
  appendix's own lesson — commit status is a provenance question, never the question of whether you
  recompute — is exactly right and is why this audit ran Tier A on untracked tables.
* **"`tm_method` still defaulted to `TM_CCOEFF_NORMED` while D1 selects `TM_CCOEFF`"** — still true
  (`FSConfig.tm_method == 5`), but the second half needs updating: **D1 was moved out of
  `DECISIONS.md` into `DECISIONS_UNVERIFIED.md` at commit `f835885` today**, pending independent
  verification. The three-way divergence is now notebook (explicit `TM_CCOEFF`) vs dataclass
  default (`TM_CCOEFF_NORMED`) vs an *unverified* decision. Inert for this notebook, which passes
  `method=` explicitly.
* **"one full single-seed pass over one ROI is ~8 s"** — measured 4.1 s (301.tiff, 31.2 Mpx) and
  6.4 s (013.tiff, 39.0 Mpx) today, including `load_roi`, both channel conversions, seed draw,
  `matchTemplate`, extraction, NMS and `od51` over 100 candidates. The 8 s figure is now an upper
  bound rather than a typical one.
* **"`FSConfig`'s defaults"** — `max_peaks=250000`, `max_detections=1000000000`,
  `score_threshold=0.5`, `channel='gray_inverted'`. All inert here: `find_and_suppress()` is not
  the caller and every value is supplied at the call site.
* Unchanged and confirmed: 7.5 µm `MIDOG_RADIUS_UM` (D7); `image_annotations`'s permissive
  `category_id=None` default (and it bit again, harmlessly — Finding 9); `midog_utils` clean in
  `git status`; non-contiguous execution is live elsewhere but **not** in this notebook.

---

## Part 4 — premises this audit inherited

Each of these is a **project decision I was instructed to apply, not a fact I verified.** Where a
verdict above leans on one, it is named.

| premise | source | what would falsify it | which verdicts lean on it |
|---|---|---|---|
| **ROI is the exchangeable unit** | D5 ("clustered at the ROI — the unit F5 §8 itself declares"); F5 §8 | Evidence that two ROIs of the same tumour type are more alike than two clicks on one ROI — i.e. that the domain, not the ROI, is the resampling unit. Partly testable here: the domain-level intervals (Finding 5) are strictly wider, so if the domain is the right unit, N10's "bar cleared" weakens at K = 20/30. | N10, N11, Findings 2, 4, 5 |
| **ROIs per stratum — measured, not assumed** | `databases/MIDOG++.json`, recomputed | — | — |
| **NMS radius = evaluation match radius = 7.5 µm** | D7 | A spacing analysis showing annotated mitoses routinely sit closer than 7.5 µm, so suppression at the match radius merges distinct objects. | Every precision and recall number, since both the pool and the match are built at that radius |
| **Recall needs `n_detections`, `precision` and `coverage_frac` beside it** | D4; memory note "recall must be paired with precision" | A pool small enough that recall is informative alone — which is arguably *this* pool at 87–100 candidates. | Finding 13 (why the dropped `coverage_frac` is harmless here) |
| **Reading burden and worst-seed behaviour are the product metrics, not median performance** | AnnotateDx product goal | A product spec where mean precision over clicks is what the pathologist experiences. | Finding 5's entire framing |
| **The length-matched precision null** | agent brief | — | **Not used.** Both arms rank the *same* 87–100-candidate pool at the same K, so the comparison is paired and a random-list null is not what is at issue. Reported as not applicable rather than skipped. |
| **Prior audits are evidence, not precedent** | agent brief | — | I re-derived the prior audit's Finding 2 (pad tracks area) and Finding 6 (`no_cap` vacuity) from artifacts rather than citing them, and both hold. I did **not** re-derive its Finding 1 (the cascade claim) — it concerns the sibling notebook, not this one, and I did not look. |

**Measured, as the brief requires:** `databases/MIDOG++.json` gives **exactly 2 ROIs per
`tumor_type`** across all 7 tumour types in `images/extra_valid` (44–200 per type in the full
dataset). So the exchangeable unit is **not** collinear with the stratum — there *is* within-stratum
replication — and Step 4.1 and Step 4.2 are genuinely two different checks here, which is why the
G = 14 and G = 7 tables in Findings 2 and 5 disagree and both had to be computed.

**Effective units, honestly counted.** For the seed axis: **64 distinct (ROI, click) pairs of 70
draws**, because 6 draws repeat another seed's click (Finding 14). For the ROI axis: **14 of 14
carry signal** — no ROI has an identically-zero delta at every budget. For the domain axis: 7,
with the two-sided sign-flip floor at 2/2⁷ = 0.0156, which K = 10 sits exactly on.

**One premise worth a premise review.** D5's bar is written as *"read at recall@K per D4"* without
naming K, while D5's own evidence sits at recall@250–5000 and this whole notebook family operates
at K = 10–30. Those are different regimes — D5 itself records the advantage flipping sign between
them — so "does a chromatin axis beat `tm_score`" has no single answer and the bar cannot be
cleared or missed without fixing K first. That is a question for `premise-reviewer`, which can test
it against the dataset and the literature; I can only note that Finding 2's verdict turns on it.

**Which verdicts would change if a premise here turned out to be wrong.** If the **ROI is not the
exchangeable unit** and the domain is, N10's "bar cleared" weakens from p = 0.003–0.008 (G = 14) to
p = 0.016–0.031 (G = 7), with K = 10 pinned at the resolution floor and canine lymphosarcoma
negative at two budgets — the effect would be "positive in 6 of 7 domains and unresolved in the
seventh", not "established". If the **7.5 µm match radius** is wrong, every precision and recall
number in Part 0 moves and nothing in Part 2 survives unexamined; the sibling
`chromatin_od_ranker_nms5um.ipynb` exists precisely because that radius is contested. Everything
else — Findings 1, 3, 7, 8, 9, 10, 11, 12, 14, 15, and the whole of Part 0's reproduction — is
arithmetic on the artifacts and on the pixels, and would stand under any of these premises.
