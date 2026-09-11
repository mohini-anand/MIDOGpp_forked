# Round 3: the round-1 and round-2 findings applied, then re-audited

**Date:** 2026-09-08 · Rounds
[1](2026-09-08-tm-recall-workload-audit.md) and
[2](2026-09-08-tm-recall-workload-audit-round2.md).

**Conflict of interest, stated up front.** The same agent made these changes and audited them.
Three things limit the damage, and none of them removes it: every table was re-derived from the
raw artefacts (`*_cells.csv`, `*_tp_ledger.csv`, `*_pool_z.npz`, `*_coverage.npz`) and never
from the notebook's own output CSVs; `evaluate.coverage_fraction` was **re-implemented from
source** rather than imported, so the coverage check is a second opinion and not a tautology;
and every quantitative claim in the closing summary was matched against the artefacts
mechanically. What this cannot audit is the *design choices* made here -- the tie-break rule,
the decision to suppress two-value medians. Those need a reader who did not write them.

## What was changed

**`recall_workload_coverage.py`** -- gate widened from one point per curve (`z = 1.0`) to four
(`z = 1.0, 2.0, 3.0` and the full pool), plus a cdf assertion. Re-ran: 70 cells, 280 reference
calls, all pass, 78 s. Docstring cost figures corrected: `coverage_fraction` measures **0.27 s**
per call, not 5.4 s, and the notebook needs it at **17,194** (cell, cutoff) pairs, not 840 --
wrong in both directions, and the conclusion (precompute) survives only because the pair count
was understated by more than the per-call cost was overstated.

**Notebook** -- 16 of 23 cells edited, re-executed clean (11 code cells, no errors).

| round-2 finding | fix |
|---|---|
| N1 Table 4 exempted from the coverage gate by argument | coverage measured at every reading depth and printed; the prose now states the inverse-function argument and says which rows survive it |
| N1b Table 4's headline was median-ROI in a worst-cell notebook | both aggregations computed and quoted: 4.9x at the median ROI, **2.5x** at the worst cell |
| N2 Table 4's per-domain block medianed two ROIs | prints `roi_lo`/`roi_hi`/`worst` per target; `per_roi`'s docstring now names the trap |
| N2b figure's dashed `z = 1.0` line was a pooled median | one dashed line per ROI, in that ROI's colour |
| N3 Table 3's `roi_gap` double-counted repeated draws | routed through `per_roi`; `roi_gap` moves in 27 of 56 (domain, K) cells |
| N4 summary paired a K = 500 gap with a K = 1000 pooled range | both budgets named, and the cell prints the pooled range at the same K |
| N5 Table 2 lined up a min-recall cell with a max-coverage cell | `cov_at_worst_recall` (matched) and `cov_most_saturated` (max) are separate columns |
| N6 "coverage mostly above 0.8" | both columns counted: 9/32 and 24/32 |
| N7 `overrun_x` printed an unbounded overrun as 1.0x | reports `inf` |
| N8 one-point coverage gate; `rank`-vs-`n_at` note; `None` guards | four-point gate; note added in the gate cell; `read_depth()` returns `None` and LODO prints skipped folds instead of silently continuing |

## Four defects introduced by the fixes, and caught here

1. **The matched-cell tie-break was arbitrary.** `min(keys, key=(recall, k))` breaks a tie by
   file name -- and below `z_max` **all 70 cells tie at recall 1.000**, so the first version
   printed `013.tiff s0, coverage 0.718` as the global "worst cell" at `z = 1.3` when the
   quantity was a 70-way tie. Ties now go to the **most saturated** cell (`245.tiff s4`, 0.995),
   which is the conservative reading and makes the column agree with `cov_most_saturated`
   exactly where recall carries no information.
2. **The two-value median survived in four more places.** Fixing Table 4's per-domain block did
   not fix `cut_roi_median` in `tolerance.csv` / `tolerance_by_count.csv`, `cut_median_roi` and
   `cov_median_roi` in `z_sweep_summary.csv`, or `median` in `recall_at_k_spread.csv` -- all
   medians over a domain's two ROIs, i.e. midpoints, unprinted but saved. They are now `NaN`
   on any group with fewer than three ROIs, and Table 3's column is replaced by
   `roi_med_lo`/`roi_med_hi`. Verified: **0** non-null on any 2-ROI row.

3. **`recall_workload_coverage.py`'s corrected docstring was itself stale.** It was patched
   before the notebook could print `COV_CALLS`, so it carried an estimate -- "over 12,000 pairs
   … ~55 minutes" -- against a measured 17,194 / 77 minutes, while correction #16 in the
   notebook claimed the script now carried the true figures. Same defect class as the 5.4 s
   claim it replaced, in the same understating direction. Now 17,194 / ~77 min, sourced from
   the run.
4. **Table 4's summary header paired a matched column with an unmatched one.** `| median ROI |
   its coverage | worst cell | its coverage |` -- the worst-cell pair is matched (both are
   `245.tiff`'s, verified at all four targets), but the median-ROI pair is two *separate*
   medians over the 14 ROIs and need not come from the same ROI. That is N5 recurring one
   level up, in the section that made matched-cell reporting a principle. Relabelled, with the
   distinction spelled out.

Two stale numbers in the closing summary were also caught, both written before the tie-break
fix: breast cancer's matched coverage (0.813 -> 0.501 became **0.942 -> 0.665**) and the
coverage call count (15,934 -> **17,194**, i.e. 77 minutes, so "most of an hour" became "over
an hour").

## Re-audit result

Independent re-derivation from the raw artefacts, compared against every saved table:

| table | scope | mismatches |
|---|---|---|
| T1 tolerance | 40 rows x 11 columns | **0** |
| T1b count | 40 rows x 11 columns | **0** |
| T2 z-sweep summary | 72 rows x 7 columns | **0** |
| T3 recall@K spread | 56 rows x 5 columns | **0** |
| T4 reading depth | 490 rows x 2 columns | **0** |
| LODO | 35 rows x 10 columns | **0** |
| ROI table (dedup median, distinct-draw count) | 70 rows x 2 columns | **0** |
| coverage npz vs re-implemented `coverage_fraction` | 5 cells x 5 cutoffs | **0** (< 1e-12) |
| coverage curves: cdf, `cum[0] == 0`, `len == n_pool + 1` | 70 cells | **0** failures |
| closing-summary numeric claims | 30 claims | **30/30 verified** |

`read_100 == n_at(z_max)` on every cell. No cell ties for the deepest reading depth, so cell
17's key tie-break is never exercised. No `median` over fewer than three ROIs survives in any
saved artefact.

## What still stands, unfixed and correctly disclosed

* **The cutoffs are still sample minima with no margin**, and most global rows are set by one
  cell. Conformal risk control is named as the machinery that would change that; it is not done.
* **Two ROIs per domain does not bound the ROI axis**, and LODO's domain-transfer is confounded
  with ROI-transfer. Both are stated in the notebook.
* **Table 4's deep rows are geometry.** `t = 1.00` at the worst cell is coverage 0.986. The
  notebook now says so rather than claiming exemption, but the number itself does not improve.
* **63 distinct clicks, not 70.** Every median is over distinct draws; the draw rule is kept for
  the high_z gate.

## Residual cosmetic items, deliberately not changed

* `x_less_worst` prints "1.0x less" where the worst cell's depth does not move (t = 0.99). It is
  accurate; it reads awkwardly.
* `worst_recall_cell` calls `coverage` inside its sort key, so the printed `COV_CALLS` includes
  tie-break lookups. The count is honest about work actually done.
* The `_read_depths` gate compares a `rank` while Tables 1 and 4 cost by `n_at(z)`. They
  coincide only because `tie_rows == 0`; the gate cell now says so explicitly.
