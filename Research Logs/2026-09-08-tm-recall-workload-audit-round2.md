# Re-audit: `tm_recall_workload_curve.ipynb` after the round-1 corrections

**Date:** 2026-09-08 · Independent re-derivation from the stored artefacts, including a
from-scratch re-implementation of `evaluate.coverage_fraction` (the installed `midog_utils`
will not import under this NumPy, so the reference was rebuilt rather than called).

Round 1: [2026-09-08-tm-recall-workload-audit.md](2026-09-08-tm-recall-workload-audit.md).

## Verdict

**All seven round-1 findings are genuinely fixed, and the new coverage machinery is correct.**
The remaining defects are smaller and of one kind: **the fixes were applied to the tables that
round 1 named, and not to the two new sections built after it.** Table 4 is exempted from the
coverage gate by argument rather than measurement, and its per-domain row re-creates the
midpoint-of-the-ROI-gap defect that correction #2 was about.

## Verified correct

| check | result |
|---|---|
| `recall_workload_coverage.py`'s curve vs. an independent re-implementation of `evaluate.coverage_fraction` | **exact** (< 1e-12) at 16 (cell, cutoff) pairs the script itself never gated — z = 1.3, 2.2, 3.0, 6.0 on four cells |
| all 70 coverage curves | monotone, `cum[0] == 0`, `len == n_pool + 1` |
| Table 1 and 1b, per-domain and global, incl. every new ROI/coverage column | reproduce exactly |
| Table 2 global and per-domain | reproduce exactly |
| Table 3 decomposition (`max_click_spread` 0.353, `roi_gap` 0.498) | reproduce exactly |
| Table 4, all 7 targets | reproduce exactly; `read_100 == n_at(z_max)` verified on all 70 cells |
| LODO | reproduces exactly — 32/35, 3 failures all `245.tiff`, gap signs 28/4/3, 3.0x overrun, cell coverage 0.392 |
| high_z gate, 8 self-tests | pass |

### The O(1) coverage curve is a good piece of work

`coverage(n) = #{p : rank_min(p) <= n} / n_probes` is exactly right given the prefix identity,
and the reverse-order fill (`for i in range(len(neigh)-1, -1, -1)`) correctly makes the lowest
rank win. It turns a 76-minute recomputation into a lookup and it is bit-for-bit faithful.

### Round-1 findings, all closed

1. Coverage now reported in Tables 1, 1b, 2 and LODO. ✔
2. `per_roi` aggregates within an ROI over distinct draws, then across ROIs; `cut_roi_lo/hi`
   and `cov_roi_lo/hi` replace the pooled median. ✔ (but see N2)
3. `click_spread` decomposed into within-ROI and between-ROI. ✔ (but see N3)
4. All four prose errors corrected — 3.0x against 8 budgeted; 28 negative, stated explicitly as
   the shipped/floored column; ROI dominance enumerated per domain instead of "throughout".
   The added claim *"only a POSITIVE gap can fail"* is sound: on floored values a gap <= 0
   means the outside cutoff is at or below the domain's own, which meets `t` by construction. ✔
5. The fabricated `0.99*100` example is gone and the code is correctly described as defensive. ✔
6. Medians over distinct draws. ✔ (but see N3)
7. `ship()` unified; `group_rows` never blanks a row; LODO guards both cell lists; the figure
   costs by `n_at(z)`; the recall@K denominator reports 1,120 informative of 1,680. ✔

## N1 — Table 4 is exempted from the coverage gate by argument, not measurement

Stated three times, unconditionally: *"Table 4 is the section that escapes this, because a
reading depth is a top-K statement, not a full-list one"*; *"the one section the coverage caveat
does not touch"*; *"the geometry caveat that gates every number above does not apply"*.
Table 4 is the only table with no coverage column. Measured:

| target | median-ROI depth | coverage at it (14 ROI medians: median / lo / hi) | worst cell | its coverage |
|---|---:|---|---:|---:|
| 0.50 | 129 | 0.011 / 0.003 / 0.154 | 3,135 | 0.236 |
| 0.90 | 1,906 | 0.149 / 0.021 / **0.699** | 13,099 | 0.723 |
| 0.95 | 2,579 | 0.199 / 0.031 / **0.820** | 18,681 | 0.861 |
| 0.98 | 4,052 | 0.294 / 0.031 / **0.899** | 26,313 | 0.947 |
| 1.00 | 9,265 | 0.588 / 0.031 / **0.983** | 33,073 | **0.986** |

At the **median ROI** the headline holds: 9,265 -> 1,906 is coverage 0.588 -> 0.149. 0.588 is
well below saturation, though not clean -- `evaluate.py` calls 0.91-0.97 the regime where the
list tiles more finely than the metric resolves, and 0.588 still means 59% of arbitrary ROI
locations sit within a match radius of some detection. But the exemption is false for the deep
rows on the ROIs the narrative rests on: `245.tiff`'s full-recall reading depth sits at 0.98.

**The stated reason is a non-sequitur, and the bias has a direction.** `evaluate.py`'s caveat is
about list *saturation*, not about whether K is fixed or derived: a derived K of 33,073 tiles the
ROI exactly as thoroughly as the full list of 36,251. Steelmanned -- "a depth is not a recall, so
the concern does not transfer" -- it still fails, because `read_depth(t) = min{K : recall@K >= t}`
is the inverse of the very curve coverage gates. A denser list puts more detections near each
mitosis, the minimum rank among them drops, and the target is reached sooner. **Saturation makes
the reading depth look shorter than the evidence warrants**, so it flatters the headline rather
than merely muddying it.

**And the headline switches aggregation convention.** *"Giving up 10% of mitoses cuts reading
five-fold"* is a **median-ROI** number, in a notebook whose governing rule -- D4, quoted in cell
8 -- is worst ROI, worst click, and whose every safety number is worst-cell. At the worst cell
the same trade is **33,073 -> 13,099 = 2.5x**, at coverage **0.986 -> 0.723**. Both halves of the
Table 4 story hold only at the median ROI and weaken together at the worst cell: the saving
halves and the coverage exemption collapses. The `worst_cell` column is printed, so the numbers
are there -- the prose does not use them.

Reporting coverage beside Table 4 costs one `per_roi` call.

## N2 and N3 — two of the notebook's own seven "Corrections applied" are overstated

The closing section claims #2 *"Aggregation is now per ROI, then across ROIs"* and #6
*"Medians are now over distinct draws"*. Both are true of the tables round 1 named and false of
one section each.

## N2 — Table 4's per-domain row re-creates the round-1 midpoint defect

`med.groupby('tumor_type').median()` (cell 17) takes a median over exactly **2** ROI values,
which is their midpoint. `per_roi` is used to build those values and then discarded.

```
canine lymphosarcoma  t=1.00: ROIs [8,502 | 29,455]  -> printed 18,978
human breast cancer   t=1.00: ROIs [  366 |  4,448]  -> printed  2,407
human melanoma        t=1.00: ROIs [2,738 | 10,101]  -> printed  6,419
canine lymphosarcoma  t=0.90: ROIs [1,632 | 12,622]  -> printed  7,127
```

5 of 7 domains at `t = 1.00` and 4 of 7 at `t = 0.90` print a number no ROI is near — the exact
failure mode of round-1 F2. Tables 1 and 1b avoid it by printing `*_roi_lo` / `*_roi_hi`;
Table 4's per-domain block prints only the midpoint. (Its **global** row is a median over 14 ROI
medians, which is fine, and it does print `roi_lo` / `roi_hi`.)

The figure's dashed "today: z = 1.0" line has the same shape: a pooled 10-cell median, which for
lymphosarcoma is ~30,692 against per-ROI clusters at ~27,816 and ~33,567.

## N3 — Table 3 is the one place a median still double-counts repeated draws

Correction #6 says medians are now over distinct draws. Cell 15 uses a plain
`.agg(med='median')` over each ROI's 5 possibly-duplicated cells, and that median feeds
`roi_gap`. De-duplicating changes `roi_gap` in **27 of 56** (domain, K) cells:

```
canine lung cancer  K=25 : 0.118 -> 0.235     human melanoma K=100: 0.196 -> 0.249
canine lung cancer  K=500: 0.176 -> 0.059     breast cancer  K=500: 0.259 -> 0.216
```

Both headline numbers are safe: `max_click_spread` is a min/max and immune, and lymphosarcoma
has 5 distinct draws in both ROIs so its 0.498 is unaffected. The printed `roi_gap` table is
what is wrong.

## N4 — the closing summary contradicts the cell it summarises

> "the largest **between-ROI** gap is **0.498** (lymphosarcoma, K = 500). A pooled 10-cell range
> mixes the two and reads 0.651"

At K = 500 the pooled range is **0.605**; cell 15 prints exactly that (*"inside a pooled 10-cell
range of 0.605 at the same cell"*). 0.651 is the pooled range at K = 1000. Same class as round
1's "24 against 12": a number from one cell paired with a number from another.

## N5 — Table 2's confound table pairs mismatched cells

The recall block uses `aggfunc='min'`, the coverage block `aggfunc='max'`, printed as
*"worst-cell coverage_frac at the same cutoffs"* directly beneath *"worst-cell recall"*. They
are different cells:

| domain, z = 3.0 | min-recall cell | its recall | **its own** coverage | printed "worst-cell" coverage |
|---|---|---|---|---|
| canine lymphosarcoma | 245.tiff s2 | 0.584 | **0.228** | 0.386 (246.tiff s3) |
| human breast cancer | 013.tiff s0 | 1.000 | **0.501** | 0.665 (094.tiff s0) |

The confound argument is *stronger* with matched cells (0.228 vs 0.501, a wider separation than
0.386 vs 0.665), so this costs the notebook nothing but the row-wise presentation invites
reading two aggregates as one cell.

## N6 — "coverage mostly above 0.8" leans on the more flattering column

Of the 32 LODO successes, `cov_roi_hi > 0.8` in 24, but `cov_roi_lo > 0.8` in only **9**
(median `cov_roi_lo` 0.738). Both columns are printed in the table the sentence points at.

## N7 — `overrun_x` is undefined at a zero budget and prints as "1.0x"

`overrun_x = actually_lost / max(budgeted_lost, 1)`. The lymphosarcoma `t = 0.99` failure lost
1 mitosis against a budget of **0** and is priced at **1.0x** in a table headed *"priced in
mitoses"* — an unbounded overrun rendered as "on budget". The adjacent `budgeted_lost` /
`actually_lost` columns disclose it, so it is cosmetic, but `inf`/`n/a` would be honest.

## N8 — minor and latent

* `z_median_click` in `ROI_T` is still a pooled 5-cell median over duplicates (informational
  only; `roi_gap` correctly uses `z_worst_click`).
* Cell 17 indexes `cutoff_for(...)[1]` with no `None` guard; cell 19 `continue`s past an
  unreachable fold without recording that it did. Dormant while `deep_recall == 1.0`.
* Cell 8's note describes `cut_roi_lo/hi` as *"the two ROIs' own medians"*, but the GLOBAL row's
  version spans 14 ROIs.
* `recall_workload_coverage.py` gates each cell at `z = 1.0` only -- one point per curve --
  while its docstring says *"bit-for-bit identical to calling it -- verified per cell"*. The
  rest of the curve is exact (verified here independently), but the in-script gate is one point.
* Cell 3's gate compares `rank` against `_read_depths` while Table 4 reports `n_at(z)`. Equal
  only because `tie_rows == 0`.

## Recommended

1. Add a coverage column to Table 4 and replace the blanket exemption with the measured
   statement: it holds at the median ROI and at low targets, not at `t >= 0.95` on the deep ROIs.
2. Print `roi_lo` / `roi_hi` in Table 4's per-domain block instead of the 2-ROI midpoint; same
   for the figure's dashed line.
3. Route cell 15's `med` through `per_roi`.
4. Fix 0.651 -> 0.605 (or re-quote K = 1000 for both).
5. Label Table 2's coverage block as the max over cells, or report the min-recall cell's own
   coverage.
