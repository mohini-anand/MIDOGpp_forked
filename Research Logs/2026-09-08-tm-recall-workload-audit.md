# Audit: `tm_recall_workload_curve.ipynb` and `recall_workload_ledger.py`

**Date:** 2026-09-08 · **Auditor:** independent re-derivation from the stored artefacts
(`results/tm_recall_workload_{cells,tp_ledger}.csv`, `results/tm_recall_workload_pool_z.npz`,
`results/tm_ccoeff_high_z_*.csv`), not from the notebook's own output CSVs.

## Verdict

**The computation is correct. The prose is not.** Every table reproduces bit-for-bit from a
clean-room reimplementation; the defects are one missing measurement, one aggregation choice
that is wrong for this design, one mislabelled axis, and four false statements of fact in
markdown. Nothing here requires re-running `recall_workload_ledger.py`.

## Verified correct

| check | result |
|---|---|
| Tables 1 / 1b (all 5 tolerances, all 5 counts, global) | reproduce exactly, incl. `z_exact`, `z_ship`, `n_kept_*`, `cut_*`, `recall_worst` |
| LODO table (35 rows) | reproduces exactly, incl. the 3 failures and 32 successes |
| `n_at(k, 1.0)` vs `n_at_current_z` (CSV) | agree in 70/70 cells |
| `len(POOLZ[k])` vs `n_pool` (CSV) | agree in 70/70 cells |
| `n_at` / `tp_kept` / `cutoff_for` semantics | correct; `searchsorted(-pool, -z, 'right')` is exactly `#{z' >= z}` |
| prefix identity | genuinely valid — `greedy_match` is one ordered pass with no re-assignment, so a truncation is a prefix |
| `_read_depths` gate | `searchsorted(tp_cum, want) + 1` **is** the rank of the `want`-th TP; the re-derivation is a real cross-check |
| "280/280" citation of high_z | verified out of `tm_ccoeff_high_z_deployable_cutoffs.csv`: 280 cells, 0 moved |
| `roi_axis_dominates` on 3-dp-rounded vs exact values | 9 either way, no flips |
| "two ROIs an order of magnitude apart" (lymphosarcoma) | 14.5x at recall 0.5, 7.7x at 0.9, 3.5x at 1.0 — defensible |
| ledger construction (NMS order, lexsort tie key, float64 `z`, seed removal from `gt_eval`) | sound |

## F1 — `coverage_frac` is never computed, and every un-budgeted recall number needs it

`midog_utils/evaluate.py` opens with, in bold: *"Read `coverage_frac` before reading any
full-list number... every metric here that is not budgeted at K is exactly that question."*
`recall_workload_ledger.py` stores `cx`/`cy` in the npz *specifically* so the notebook can
recompute it ("so the notebook can evaluate the retained count and `coverage_frac` at ANY
cutoff exactly"). The notebook never reads `cx`/`cy` and never mentions coverage.

Tables 1, 1b, 2 and the LODO section are **all** un-budgeted full-list recall. Recomputed at
the shipped cutoffs (stride 16, each cell's own `match_radius_px`):

| cell | z=1.0 | z=1.3 (global t=1.00) | z=1.6 (t=0.98) | z=2.2 (t=0.90) | z=3.0 |
|---|---|---|---|---|---|
| 245.tiff s2 (lymphosarcoma) | 0.996 | **0.988** | 0.961 | 0.728 | 0.228 |
| 246.tiff s0 (lymphosarcoma) | 0.992 | **0.979** | 0.938 | 0.710 | 0.318 |
| 301.tiff s3 (mast cell)     | 0.987 | 0.976 | 0.962 | 0.890 | 0.617 |
| 300.tiff s1 (mast cell)     | 0.971 | 0.946 | 0.902 | 0.723 | 0.370 |
| 459.tiff s3 (STS, sets the global cutoff) | 0.901 | 0.865 | 0.812 | 0.632 | 0.294 |
| 013.tiff s0 (breast)        | 0.813 | 0.718 | 0.650 | 0.567 | 0.501 |

At the shipped global cutoff `z = 1.3`, 98–99% of arbitrary ROI locations on lymphosarcoma
are already within the match radius of some detection. `recall_worst = 1.0000` there is close
to a statement of geometry.

Second-order: the "cliff above `z_max`" is confounded with coverage decay. Lymphosarcoma's
coverage falls 0.996 -> 0.23 between z=1.0 and z=3.0 while breast cancer's falls only
0.81 -> 0.50 — the same ordering as the recall cliff the summary attributes to the domain's
difficulty. This does not refute *"a single global cutoff above 2 is a bet on the domain"*,
but the bet may be on candidate density rather than on discriminability.

## F2 — the domain and global medians land in the gap between the two ROIs

With 2 ROIs x 5 clicks and per-ROI clusters that do not overlap, the 10-cell median is the
midpoint of the 5th and 6th order statistics — i.e. the midpoint of the ROI gap.

Canine lymphosarcoma, `n_at(z = 1.3)` per cell:

```
246.tiff: 22119  23408  23999  24907  25174
245.tiff: 30322  32285  33067  33708  36251
reported n_kept_median = 27748 = (25174 + 30322) / 2
```

**No cell in the domain is within 5,000 of the reported median.** Canine lung cancer at
z = 2.9 is the same shape: 201.tiff spans 2,636–4,320, 233.tiff spans 6,237–11,962, reported
median 5,278. This affects every `*_median` column in Tables 1, 1b and 2, and LODO's
`cut_median_pct`. The global row is a median over 70 cells whose per-ROI medians span
13,118–33,067 at z=1.3, which is what the headline *"one annotation per ROI doubles the
median saving"* rests on.

## F3 — Table 3's "click spread" is mostly an ROI spread

`RK` is grouped by `(tumor_type, budget)` — 10 cells = **2 ROIs x 5 clicks** — and the derived
column is named `click_spread`, printed under *"spread across the 10 cells of each domain
(best - worst click)"* and the worst column under *"worst click (D4)"*.

Decomposed at the headline cell (lymphosarcoma, K = 1000):

```
245.tiff: 0.236 0.292 0.371 0.427 0.438     within-ROI spread 0.202
246.tiff: 0.852 0.844 0.861 0.800 0.887     within-ROI spread 0.087
reported "click_spread" 0.651 -> between-ROI 0.481, max within-ROI 0.202
```

The ranges are disjoint. The closing summary's *"the click moves it by up to two thirds of the
whole scale"* is wrong by ~3x; the largest genuine within-ROI click spread anywhere in the
table is 0.353 (canine lung cancer, K = 50/100). Cell 7 gets this right — it separates the
axes and names lymphosarcoma as the domain where the ROI axis dominates — and cell 15 then
pools them and calls the result "the click".

## F4 — four false statements in the markdown

1. **"24 mitoses missed against 12 budgeted, a 2x overrun."** The worst cell is
   `(245.tiff, seed 2)`, m = 89, kept 65. Its t = 0.90 budget is `89 - ceil(0.9 x 89) = 8`.
   **24 against 8 — a 3x overrun.** "12" is not this cell's budget (246.tiff's is 11). The
   error understates the notebook's own most decision-relevant finding.
2. **"`transfer_gap_z` is negative in 30 of 35 rows."** The printed and saved column shows
   **28** negative, 4 zero, 3 positive. 30/5 is the split on the *exact* gaps; flooring both
   sides to 1 dp before differencing collapses four signed gaps to 0.0.
3. **"The one domain where the outside cutoff sat above its own is precisely the one that
   failed."** On exact values, canine soft tissue sarcoma at t = 1.00 also has a positive gap
   (+0.0352: outside 1.3430 vs own 1.3079) and does not fail — flooring both to 1.3 rescues
   it. True at domain level, false at row level.
4. **"lymphosarcoma and neuroendocrine tumour throughout."** Each is 4 of 5, not 5 of 5
   (lymphosarcoma is False at t = 1.00, neuroendocrine at t = 0.90), and canine soft tissue
   sarcoma at t = 0.95 — the ninth True — goes unmentioned.

## F5 — the stated justification for integer `ceil` is fabricated

Cell 4: *"`int(np.ceil(0.99 * 100))` is `100`, not `99` -- `0.99*100` is `99.00000000000001`
in binary floating point."* Both halves are false: `0.99*100 == 99.0` exactly, and
`int(np.ceil(0.99*100)) == 99`. Swept over the five tolerances and m = 1..1000, the float and
integer forms **never** disagree. The `n_want` code is correct and harmless; only the reason
given for it is wrong. Self-test 7 already reports 0/350 and hedges to "belt-and-braces",
which contradicts the markdown two cells earlier.

## F6 — 7 duplicate draws are counted as independent cells in every median

`draw_seed_with_retry` is with-replacement: 233.tiff drew the same annotation 3 times
(seeds 0/1/4), and 094/201/300/403/529 twice each — 63 distinct of 70. The notebook discloses
this only as *"a worst-of-5 secretly a worst-of-4"*, which covers the min/max statistics.
It does not cover the medians, where a triplicated cell gets triple weight:

| canine lung cancer, t = 1.00 | over 10 cells | over 7 distinct draws |
|---|---|---|
| `n_kept_median` | 5,278 | 4,320 |
| `cut_vs_z1_median_pct` | 64.7% | 75.7% |

Smaller shifts in breast cancer (61.6% -> 57.4%) and neuroendocrine (24.5% -> 28.8%).

## F7 — latent code issues (do not bite on this data)

* **`tolerance_row` / `count_row` handle unreachable cells differently.** `tolerance_row`
  excludes the cell from the vote but still evaluates it and reports `unreachable_cells`;
  `count_row` returns `None` and silently drops the entire group row. Dormant only because
  `deep_recall == 1.0` in all 70 cells.
* **LODO `own` is unguarded.** `min(cutoff_for(k, n_want(k, t))[0] for k in inside)` will
  raise `TypeError` if any *inside* cell is unreachable; only the *outside* list is guarded.
* **The figure plots `rank`, the tables cost by `n_at(z)`.** `cutoff_for`'s docstring argues
  at length that these differ under score ties and that `rank` "would understate the
  workload". Identical here (`tie_rows == 0`), so the panel is numerically right, but it
  contradicts the stated rule.
* **560 of the 1680 recall@K comparisons are self-comparisons.** The baseline is
  `recall_at_k(k, K, CURRENT_Z)` and `CURRENT_Z` is one of the three cutoffs tested, so a
  third of the denominator is 0 by identity. 1120 comparisons are informative.
* **`recall_worst` is rounded to 4 dp before being asserted against `tolerance - 1e-9`.** Safe
  only because `z_ship <= z_exact` guarantees the inequality analytically.
* **`z_ship = floor(z_exact * 10) / 10` assumes the cutoff is shipped at exactly 1 dp.**
  `ship_dp` is a parameter in `tolerance_row`/`count_row` but hard-coded to `* 10 / 10` in the
  LODO cell.

## F8 — framing tensions (not errors)

* The intro rejects high_z's deliverable for being *"a minimum -- over mitoses, then over
  clicks, then over domains"*. The tolerance axis relaxes only the innermost minimum; the
  cross-cell aggregation is still a hard min, and 7 of the 10 global rows are set by a single
  cell (`245.tiff` s2 or `459.tiff` s3). Disclosed under "a minimum still carries no margin",
  but the opening implies more was fixed than was.
* The workload axis operates at 12k–36k candidates. "One annotation per ROI doubles the median
  saving" moves the median list from ~19,500 to ~15,750 — both far outside any regime a
  pathologist reads. The decision-relevant regime is Table 3's (K in the hundreds), where
  recall is 0.24–0.89 and the tolerance tables say nothing. The two halves of the notebook do
  not meet.
* LODO holds out a *domain*, but each domain is 2 ROIs, so domain-transfer and ROI-transfer
  are fully confounded. The notebook's own caveat ("a property of this particular 7-domain
  mixture") is right but does not name the confound.

## Recommended actions

1. Compute `coverage_frac` at every shipped cutoff and report it beside `recall_worst`
   (`cx`/`cy` are already in the npz; the code above runs in seconds).
2. Report per-ROI medians, or the min/max over ROIs, instead of a 10-cell median; or state
   explicitly that the median falls in the ROI gap.
3. Rename `click_spread` -> `cell_spread`, add a within-ROI decomposition, and fix the
   closing-summary sentence.
4. Fix the four prose errors (F4) and delete the false float example (F5).
5. Report medians over distinct draws, or switch to `tm_variant_sweep.draw_seeds` and drop the
   high_z gate to seed 0 only (where it already lives for `read_*`).
