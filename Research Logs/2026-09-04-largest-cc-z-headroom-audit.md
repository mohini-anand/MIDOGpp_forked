# Audit: `tm_threshold_axis_sweep_largest_cc.ipynb`, and the z-threshold headroom it leaves on the table

Date: 2026-09-04
Scope: a second, independent audit of `tm_threshold_axis_sweep_largest_cc.ipynb` (lcc),
prompted by the question "can `CURRENT_Z` be raised without losing mitoses?". Everything
below is re-derived from `results/tm_ccoeff_threshold_axis_sweep{_v2,_largest_cc}.csv`,
`results/tm_ccoeff_headtohead_seed_variance.csv` and the `midog_utils` source. No sweep was
re-executed for the audit findings.

`Research Logs/2026-09-03-tm-axis-sweep-audit.md` already covers this notebook thoroughly and
its findings 1-12 were re-checked and are correct. What follows is what that audit did **not**
catch, plus the answer to the z question.

---

## Tier 1 -- changes a stated conclusion

### A. The seed-noise claim that carries the reading-depth verdict is false for 2 of 7 ROIs

lcc's closing summary, third bullet:

> `read_95_delta` on `tm_score` ranges +1416 to -1547 across domains, and **every one of those
> is smaller than that ROI's own seed-to-seed SD** (260-2744 candidates,
> `results/tm_ccoeff_headtohead_seed_variance.csv`, same ROIs and method, 5 seeds).

The same sentence is in finding 5 of `2026-09-03-tm-axis-sweep-audit.md`, which is where lcc
took it from. Paired ROI by ROI against `arm='tm_ccoeff'`, `metric='read_95'`, it does not hold:

| ROI | domain | `read_95_delta` (lcc - v2) | that ROI's 5-seed SD | \|delta\| < SD |
|---|---|---:|---:|---|
| 301.tiff | mast cell tumor | +527 | 2208 | yes |
| 201.tiff | lung cancer | 0 | 2744 | yes (null control) |
| 246.tiff | lymphosarcoma | -717 | 1635 | yes |
| **459.tiff** | **soft tissue sarcoma** | **+1416** | **1199** | **no (1.18x)** |
| 094.tiff | breast cancer | -93 | 260 | yes |
| 548.tiff | melanoma | -1199 | 1325 | yes |
| **402.tiff** | **neuroendocrine tumor** | **-1547** | **593** | **no (2.61x)** |

The two failures are not marginal, and `402.tiff` is the **largest** reading-depth movement in
the whole table -- 2.6x its own ROI's seed SD.

The likely mechanical cause: the delta list and the SD list were compared as unordered sets
rather than joined on `file_name`. `548.tiff`'s delta (-1199) is numerically equal to
`459.tiff`'s SD (1199), which makes the mismatch easy to miss by eye.

**What survives.** The verdict "reading depth is not measurable at this sample size" is
probably still right, but not for the reason given. Its remaining support is (i) the three
exact endpoint hits (seed 0 lands on the 5-seed max for `301.tiff` and the min for `246.tiff`
and `548.tiff`), which is the load-bearing argument and is unaffected, and (ii) the fact --
stated in the audit log but dropped when the claim moved into the notebook -- that the variance
CSV comes from the earlier `tm_variant_sweep` configuration (no padding, NMS at the match
radius, different pool scope), so those SDs are a scale reference and not an error bar for
*this* run. **Both the notebook's closing summary and finding 5 of the 2026-09-03 audit log
need correcting**, or the next reader re-derives the wrong claim from the log.

### B. The notebook computes the answer to "how high can `z` go?" and discards it

Cell 38 computes `tp_z_min = float(tp['z'].min())` per domain, uses it only as the
`frac_fp_in_tp_range` cut point, and never prints it. That value **is** the exact answer, and
the derivation needs no grid at all:

1. `pool` is already score-descending when it is built (`extract_peaks` returns a globally
   lexsorted list, `nms_by_distance` preserves that with `kind="stable"`, and the boolean
   `ok` mask preserves row order).
2. Each `tm_score` arm's candidate frame is `pool[pool['score'] >= cut]`, re-sorted stably by
   score -- i.e. exactly a **prefix** of the deep pool.
3. `evaluate.greedy_match` walks detections in rank order, so matching a prefix is literally
   the first N steps of matching the whole list. `full_list_recall` at cut `z` is therefore
   `#{TP-bucketed detections with z >= cut} / n_gt_mitotic`, exactly.
4. The deep pool's own TP count equals `n_gt_mitotic` in all 7 domains (217, 17, 115, 130,
   81, 238, 104 -- see the main loop's own printout), so recall is 1.0 at the deep floor.

Hence **`z_max` = `tp['z'].min()`**, exactly, inclusive (the cut is `score >= cut`). The
`Z_LEVELS` grid brackets it (`full_list_recall == 1.0` at `z = 1.5` in all 7 domains) but never
resolves it, and for `201.tiff`, `094.tiff` and `402.tiff` the grid's own top (3.0) never
loses a mitosis at all, so the grid cannot even bracket those three.

The practical consequence: `CURRENT_Z = 1.0` is far below the binding constraint in every
domain, and the notebook's data already proves it.

---

## Tier 2 -- real slips, conclusion survives

### C. The `chromatin_od` arm is now outside `DECISIONS.md` D5's standing constraint

D5 (2026-09-04): *"Any run that ranks by a chromatin statistic -- as primary axis, **as a
compared axis**, or as a re-measurement of this decision -- must use the 14 ROIs in
`images/extra_valid/`."* lcc ran on the 7-ROI `select_domain_images(images_dir='images')`
draw at a single seed, which is the exact base D5 exists to stop being used.

This is not a defect in the run -- lcc executed at 12:48 and D5 was written at 17:42 the same
day -- but any *forward* use of lcc's `chromatin_od` numbers (deliverable (d)'s `read_95`
column, deliverable (f)'s `frac_fp_in_tp_range__chromatin_od`, the whole `sep_cmp` table) is
now non-decision-grade by the project's own rule, and nothing in the notebook says so.

### D. `pools[fn]` can order tied scores differently from the arm frames it is compared against

Main loop: `pool_sorted = pool.sort_values('score', ascending=False)` -- pandas' default
`kind='quicksort'`, which is **not** stable. `compare._rank`, which produces every row of the
saved CSV, uses `kind="mergesort"` and documents why: *"so a tie block keeps the caller's
incoming order and a run is reproducible."* `extract_peaks`' docstring makes the same argument
at length, having measured a real 2-of-8747 coordinate disagreement from exactly this.

`largest_tie_block` is 2 in every domain, so at most one adjacent pair can swap and the
practical risk is nil -- but `pools[fn]` is what deliverables (e) and (f) are computed from,
and it is the one frame in the notebook not built under the convention the rest of the repo
enforces. `kind='mergesort'` is the one-word fix.

---

## Tier 3 -- worth noting

### E. "All three deltas move together in every cell" has two counterexamples, hidden by rounding

Cell 30's prose. `301.tiff` at `z = -1.0` and `z = +0.5` has `n_detections` ratio 0.9999 and
0.9957 -- the list got **shorter** -- while `coverage_frac` rose by +0.0008 and +0.0006. The
printed grid rounds the ratios to `1.00` and the deltas to 4 dp, so neither is visible. The
recall delta is 0 in both cells, so the argument the sentence is making is untouched; the
sentence is just stronger than the data.

### F. `coincide_n`'s denominator does not match its own description

Cell 12's markdown: *"among candidates where **both** succeed, how often the two return the
byte-identical bbox (`coincide_match / coincide_n`)"* -- but `coincide_n = n_tighten_ok`.
A candidate where `tighten_box_otsu` succeeds and `largest_cc_box` fails (possible: the largest
component can breach `max_area_frac` or `min_solidity` when the click's own component does not)
sits in the denominator and can never reach the numerator. The corrected probability cell
already sidesteps this by using `n_pool` and `largest_cc_ok` as the two live denominators, so
nothing downstream depends on it.

---

## Re-checked and confirmed correct

Everything else examined held. Specifically re-verified in this pass, beyond the 2026-09-03
audit's own list:

* **The corrected coincidence probabilities.** 0.4014 (uniform over the border-filtered pool),
  0.4165 (uniform over the `largest_cc_ok` subset), 0.8699 (the superseded conditional), and
  the stated 2.2x overstatement -- all reproduce exactly.
* **Every arithmetic cell.** All 14 `sep_cmp` deltas, all seven `read_95_delta` and
  `n_detections_ratio` values, the 8-cell sub-1.0 recall count, and the deep-floor pool ratios
  quoted in the closing summary (25447 vs 25449 = 1.000; 246.tiff 1.017; 459.tiff 1.389;
  548.tiff 1.405) reproduce from the two CSVs.
* **The `201.tiff` null control is exact on *all* 39 shared CSV columns**, both arms, all nine
  z levels, all eight budgets -- 144 rows, zero differences. Stronger than the 13 columns the
  2026-09-03 audit checked.
* **The separability-vs-`read_95` sign-agreement claim** (3 of 6 agree: lymphosarcoma,
  melanoma, neuroendocrine; 3 disagree: mast cell, soft tissue sarcoma, breast) is exactly
  right.
* **`OD_PAD = 25` is the minimum sufficient margin and is correct.** `chromatin.score_detections`
  defaults `window=tm.BASE_SIZE=51`; `read_padded_patch` needs `ix-25 >= 0` and `ix+25 < w`,
  which a 25 px `BORDER_REPLICATE` margin supplies exactly at both edges. `nan_rate` is
  0.000000 on every row of the saved CSV.
* **`largest_cc_box` mirrors `tighten_box_otsu`'s `binary` path exactly** -- same
  `cv2.normalize` to uint8, same `THRESH_BINARY+THRESH_OTSU`, same `connectivity=2`, same three
  gates, same half-open `(y0, y1, x0, x1)` convention -- with only component selection changed.
  `_odd_local` is byte-equivalent to `seed_selection._odd`.
* **`PAD = (base_size-1)//2` makes `valid.all()` true by construction**, for any odd template
  size: `matchTemplate` on the padded map returns exactly `(H, W)` rows/cols and
  `fused_p[PAD:PAD+H, PAD:PAD+W]` is exactly the valid region.
* **`largest_tie_block == 2` and `nan_rate == 0.0` on every row** of the saved CSV.
* **Audit finding 10 reproduces on this notebook too.** `tm_score` `recall_at_budget` is
  *exactly* constant across all nine z levels at K=100 and K=250 in all 7 domains.
  `chromatin_od` is the only arm where z moves the product metric at all (largest movement:
  neuroendocrine recall@100, 0.3173 -> 0.3654 between z=1.0 and z=3.0).

---

## The z answer, and what raising z does and does not buy

Because `full_list_recall` on the `tm_score` axis is a prefix statistic (finding B), two
things follow immediately and neither is in the notebook:

* **`_read_depths` cannot improve.** It is `np.searchsorted` into a prefix of the *same*
  `tp_cum` array, so every `read_*` on the `tm_score` arm is either **identical** to its
  deep-pool value or **NaN** (target now unreachable). Raising z cannot shorten reading depth.
* **`recall_at_budget` cannot change either**, wherever `n_detections >= budget` -- audit
  finding 10, confirmed above for K=100 and K=250 across the full grid.

So raising z does **not** reduce what the pathologist reads. What it does buy:

1. **Full-list precision** -- the same TP count over a shorter list.
2. **Downstream compute** over the candidate set (feature extraction, FP filtering, any
   re-ranking pass), which is linear in pool size.
3. **The `chromatin_od` arm's top-K**, which is the one place z genuinely moves recall@K,
   because od-ranking promotes low-`tm_score` candidates that truncation removes.
4. **A much stronger recall claim.** `coverage_frac` at z=1.0 is 0.896-0.999, so the parent's
   `full_list_recall == 1.0` is close to the tautology the 2026-09-03 audit's finding 1 called
   out. At z=2.5 coverage is 0.434-0.771. "Full recall with coverage 0.55" is evidence;
   "full recall with coverage 0.99" is geometry.

`tm_threshold_axis_sweep_largest_cc_high_z.ipynb` measures the per-domain `z_max` exactly,
reports the order statistics around it (a threshold set at the minimum over TPs is by
construction hostage to one annotation -- D4's own argument), splits the binding annotations by
`unanimous`, and re-runs `tp_z_min` over 5 seeds so the number is not a single-seed number.

---

## Result (run 2026-09-04, `tm_threshold_axis_sweep_largest_cc_high_z.ipynb`)

The new notebook reproduces lcc's saved CSV **byte-identically** on all 1008 grid rows x 37
shared columns, so every number below is a clean addition to that run rather than a re-run of
it. `z_max` is verified as exactly binding: the arm cut at `score >= tp_score_min` reaches
`full_list_recall == 1.0` in all 7 domains, and the arm one notch higher
(`score > tp_score_min`) loses exactly one mitotic figure in all 7.

**Seed 0, per domain (the headroom, oracle):**

| domain | `z_max` | n at z=1.0 | n at `z_max` | cut | coverage 1.0 -> `z_max` |
|---|---:|---:|---:|---:|---|
| canine soft tissue sarcoma | 1.891 | 20,892 | 15,466 | 26.0% | 0.965 -> 0.867 |
| canine cutaneous mast cell tumor | 1.987 | 22,227 | 15,929 | 28.3% | 0.974 -> 0.871 |
| canine lymphosarcoma | 2.722 | 26,677 | 6,268 | 76.5% | 0.992 -> 0.439 |
| human melanoma | 2.915 | 20,810 | 7,949 | 61.8% | 0.939 -> 0.547 |
| human neuroendocrine tumor | 3.187 | 27,280 | 5,457 | 80.0% | 0.982 -> 0.381 |
| canine lung cancer | 4.076 | 16,529 | 665 | 96.0% | 0.896 -> 0.058 |
| human breast cancer | 4.768 | 26,147 | 4,052 | 84.5% | 0.942 -> 0.275 |

The coverage collapse is the point: at `z = 1.0`, `full_list_recall == 1.0` is largely
geometry (finding 1 of the 2026-09-03 audit); at these cutoffs it is not.

**Over 5 seeds (the deployable answer).** Worst of 35 (domain, seed) cells is
**`z = 1.3079`** (`459.tiff`, seed 3) -- seed 0's own per-domain minimum of 1.8908 overstates
the globally safe cutoff by 0.58 z. Per-domain worst-over-seeds (**floored**, never rounded to
nearest): soft tissue sarcoma 1.3078, mast cell 1.6656, neuroendocrine 1.8251, lymphosarcoma
1.8915, melanoma 2.2923, lung 2.9928, breast 3.4876.

Every cutoff below is **measured** -- each seed's pool is retained and re-evaluated through
`compare.evaluate_arms` at each cutoff -- not interpolated off a grid.

| cutoff rule | median pool cut vs `z = 1.0` | range | coverage at the cutoff |
|---|---:|---|---|
| global, worst over all 35 cells (`z = 1.308`, ship `1.3`) | 10% | 7-14% | 0.84-0.98 |
| per-domain, worst click within domain | 39% | 9-80% | 0.26-0.96 |
| seed-0 oracle (not available at inference) | 62% | 19-87% | 0.06-0.87 |

**Three reasons the deployable number is an upper bound, not a recommendation.** (i) It has no
margin and has not converged: the running minimum over seeds is 1.890 / 1.614 / 1.614 / 1.307 /
1.307, falling twice in five draws, most recently at k=4. (ii) The **ROI half of D4's rule is
untested** -- one ROI per domain here, and D4's stated reason for "worst ROI, worst click, per
domain" is precisely that the ROIs cluster by tumour type. (iii) `draw_seed_with_retry` samples
with replacement across seed indices (`f1_seed_sweep.py`'s rule, unlike
`tm_variant_sweep.draw_seeds`), so 2 of the 35 draws repeat a click: **33 distinct clicks, not
35**. Neither duplicate sets a domain minimum.

**A rounding trap worth recording.** `z_max` is a *minimum*, so any round-to-nearest applied to
it -- for storage or for display -- can move it up and silently drop the binding mitosis.
Storing it as `round(z_max, 4)` rounded **up** in 18 of 35 cells and made the cut lose a mitosis
in exactly the cells that set a minimum; an assertion in the deployable-cutoff cell caught it.
Every stored cutoff is now full precision and every displayed one is floored. `round(1.6657, 1)
= 1.7` would put 5 of the 7 per-domain recommendations above their measured safe value.

**Coverage, read correctly.** `full_list_recall == 1.0` at `z_max` is true **by construction**
-- `z_max` is defined as the lowest TP score -- so it carries no information at any coverage,
and the earlier reading ("the identical recall number carries more evidence at low coverage")
was wrong. The informative quantity is the *value* of `z_max`. Coverage does bite at the
**deployable** cutoff, and there it is unflattering: 0.84-0.98 at `z = 1.308`, squarely the
regime `evaluate.py` calls geometry rather than evidence.

**What it buys.** Measured over 10 cutoffs x 7 domains x 8 budgets: `recall_at_budget` on
`tm_score` moved in **zero** cells where the list was longer than the budget, and of 420
`read_*` cells **zero** changed and 15 became `NaN`. Raising the cutoff cannot reduce reading
burden. It buys full-list precision -- the only full-list quantity that moves, since recall
is pinned at 1.0 -- which goes from 0.10-1.14% to 0.84-2.99% (1.4x-24.8x, still under 3%
everywhere); downstream compute proportional to pool size; and very little else. The
`chromatin_od` top-K gain does **not** survive to a shippable cutoff: the +0.118 recall@100
headline is at the oracle `z_max`, while measured at the global cutoff the per-domain median
gain is +0.0000 in all seven domains at recall@100 (moving at all in 3 of 35 cells) and at most
+0.0046 at recall@250; at the per-domain cutoffs the largest median gain is +0.0174. This had
to be measured rather than bracketed -- the arm's `recall_at_budget` is not monotone in the
cutoff. Non-decision-grade under D5 regardless.

**The margin.** Giving up the single deepest mitosis per domain raises the cutoff by a further
0.11-1.25 z (median 0.31) -- **seed 0 only, and measured against the seed-0 oracle baseline**,
so it sizes the shape of the trade-off, not a shippable saving. The binding annotation is `unanimous` in 6 of 7 domains, so
"Still open, deliberately" #3's contested-tail hypothesis does not excuse crossing it here.

The `recall_at_budget` invariance check is applied **row-wise** (drop individually truncated
rows) rather than per cell: excluding a whole `(file, budget)` cell because one of its cutoffs
truncates would have skipped 5 of 56 cells, including every budget >= 1000 on `201.tiff`, the
domain whose pool shrinks most. Row-wise it covers 56/56 cells and 554/560 arm-rows, with 0
violations.

New artefacts: `results/tm_ccoeff_threshold_axis_sweep_largest_cc_high_z.csv`,
`results/tm_ccoeff_high_z_seed_table.csv`, `results/tm_ccoeff_high_z_tp_z_values.csv`,
`results/tm_ccoeff_high_z_volume_curve.csv` (candidate count vs cutoff at 0.05 z resolution
per (ROI, seed), so any proposed threshold is costable without re-running the sweep), and
`results/tm_ccoeff_high_z_deployable_cutoffs.csv` (all 1680 arm rows at the current, global and
per-domain cutoffs for every (ROI, seed), both axes).
