# Does raising the Otsu bbox-tightening threshold help end-to-end detection, or just geometry?

Date: 2026-09-03
Scope: `bbox_headroom_end_to_end.py`, 150 `find_and_suppress` runs (10 decision-grade ROIs x 5
seeds x 3 arms). Raw data: `results/bbox_headroom_end_to_end.csv`. Independently checked against
the raw CSV before this write-up; three corrections from that check are folded in below (marked
inline) rather than silently fixed.

## Why this run exists

`2026-09-03-bbox-threshold-sweep.md` found that raising `tighten_box_otsu`'s effective foreground
threshold un-merges a majority of the multi-nucleus connected components that plague dense-cellular
domains (lymphosarcoma, mast cell tumor) — `frac=0.15` (the `"headroom"` method) resolves about
half of the flagged `flagged_accepted` population at roughly a third of `multiotsu`'s collateral
cost in the domains that matter. But that sweep's own stated limit was explicit: it measured
seed-selection **geometry** only (box size/shape/solidity) — "downstream template-matching/
discrimination quality... remains completely untested." This run closes that gap with an actual
end-to-end `find_and_suppress` pass.

**Answer, stated up front**: raising the threshold — either `headroom@0.15` or `multiotsu` — does
not help end-to-end detection, and directionally hurts it, on this ROI set. No comparison reaches
p<0.05 at n=10. See "Conclusion" for the mechanism and why this tempers the geometry-only sweep's
optimism rather than contradicting it.

## Methodology

Matches `2026-09-02-next-steps-plan.md`'s established, already-audited practice (its "0. Terms"
and "6. Statistical power" sections) exactly, not a new design:

- **10 decision-grade ROIs** (`n_mitotic >= 15`, `tm_variant_sweep.DECISION_MIN_MITOTIC`, applied
  by count): computed fresh against the currently-downloaded `images/` via
  `dataset.load_annotations()` — **094, 201, 202, 245, 246, 300, 301, 402, 459, 548.tiff**.
  Matches the reference list exactly.
- **5 seeds per ROI**: `rng = np.random.default_rng(seed_index)`, a *fresh* generator per (ROI,
  seed_index, arm) triple, not one generator advanced across arms — so each arm's draw is the
  same nominal index against its own filtered pool, never a state-drifted one.
- **`fs.FSConfig()` defaults, unmodified** — confirmed directly from the dataclass in
  `find_and_suppress.py` (`score_threshold=0.5`, `n_angles=1`, `flips=(False,)`, single-scale),
  matching `design_choices.md` §6.
- **`experiment.run_one_image(..., run_baselines=False)`** called directly for all 150 (ROI,
  seed_index, arm) triples — no baseline/blob comparison, no reimplemented evaluation logic.
  `recall_at_budget` at budgets (100, 250, 500) derived via `evaluate.recall_at_k(det_out["bucket"],
  n_gt_mitotic_eval, k=budget)` on the already-bucketed `detections` frame `run_one_image` returns
  — algebraically identical to `compare.evaluate_arms`'s own `tp_at_budget / n_mit` with
  `delivered = min(budget, len(list))`.
- **Aggregation: ROI-level median across 5 seeds, per arm, per budget — not per-cell.** A
  two-sided sign test (`scipy.stats.binomtest(k, n=10, p=0.5, alternative='two-sided')`) across
  the 10 ROIs is the headline comparison, exactly as the Sep-2 plan's own audited standard
  requires; no per-cell Wilcoxon is used or reported as a substitute.

Wall clock: 1036.1s (17.3 min) for all 150 calls, ~6.9s/call. All 150 succeeded — no pool-emptying
exceptions, no errors.

## Seed-pool mismatch: tracked, and it turns out to matter a lot

`pick_seed`'s candidate pool depends on `bbox_method`, so the same `rng` draw at the same (ROI,
seed_index) can select a *different* annotation across arms. Per-ROI count of the 5 seed draws
where all three arms picked the same `seed_ann_id`:

| ROI | same (all 3 arms) | diverged |
|---|---:|---:|
| 094 | 1/5 | 4/5 |
| 201 | 5/5 | 0/5 |
| 202 | 4/5 | 1/5 |
| 245 | 0/5 | 5/5 |
| 246 | 0/5 | 5/5 |
| 300 | 1/5 | 4/5 |
| 301 | 3/5 | 2/5 |
| 402 | 4/5 | 1/5 |
| 459 | 2/5 | 3/5 |
| 548 | 1/5 | 4/5 |
| **total** | **21/50** | **29/50** |

Divergence is worst on the dense domains (245, 246: 0/5 same) — expected, since that is exactly
where the foreground-filter pool differs most between methods.

**This confound strengthens the "raising the threshold doesn't help" finding, it does not weaken
it.** Every ROI where `headroom` or `multiotsu` ever wins at any budget — {094, 245, 246, 300} —
has same-seed overlap of 0 or 1 out of 5 (minimal seed control: most of that "win" could be a
different, easier annotation, not a better template). Every ROI with real seed control — 202
(4/5), 301 (3/5), 402 (4/5) — favors `binary` at every budget; 201 (5/5 same) is an exact tie at
every budget for both arms **at the ROI median**, but the median hides a within-ROI split: at
budget=500 on 201's two `seed_ann_id=4477` seeds (2 of the 5), `binary` beats both arms 0.353 vs.
0.294; the other three seeds are genuine ties. So the single most seed-controlled ROI in the set
is, at the seed level, 2 wins for `binary` and 3 ties for the raised-threshold arms — not a null
result the way the median makes it look. None of the well-controlled ROIs favor the
raised-threshold arms.

A same-seed-only stress test makes this concrete: restricting to the 21 (ROI, seed) cells where
all three arms drew the *same* annotation, and comparing `recall_at_budget_500` seed-by-seed
(pseudoreplicated at the seed level, not a substitute for the ROI-level test — reported as an
informal check only):

| comparison | losses | ties | wins | n |
|---|---:|---:|---:|---:|
| headroom vs. binary | 13 | 6 | 2 | 21 |
| multiotsu vs. binary | 15 | 6 | 0 | 21 |

Same direction as the ROI-level result, if anything more one-sided — on matched templates,
`multiotsu` never beats `binary` once in 21 same-seed comparisons.

## ROI-level median recall_at_budget (binary / headroom@0.15 / multiotsu)

| ROI | @100 | @250 | @500 |
|---|---|---|---|
| 094 | 0.296 / 0.296 / 0.198 | 0.407 / 0.432 / 0.358 | 0.580 / 0.543 / 0.543 |
| 201 | 0.000 / 0.000 / 0.000 | 0.059 / 0.059 / 0.059 | 0.059 / 0.059 / 0.059 |
| 202 | 0.200 / 0.000 / 0.000 | 0.400 / 0.000 / 0.000 | 0.600 / 0.000 / 0.133 |
| 245 | 0.000 / 0.011 / 0.022 | 0.022 / 0.011 / 0.045 | 0.034 / 0.011 / 0.101 |
| 246\* | 0.148 / 0.296 / 0.226 | 0.148 / 0.348 / 0.330 | 0.148 / 0.478 / 0.391 |
| 300 | 0.022 / 0.083 / 0.078 | 0.078 / 0.150 / 0.144 | 0.122 / 0.244 / 0.189 |
| 301 | 0.147 / 0.065 / 0.051 | 0.221 / 0.115 / 0.115 | 0.309 / 0.175 / 0.175 |
| 402 | 0.212 / 0.135 / 0.106 | 0.327 / 0.269 / 0.144 | 0.404 / 0.375 / 0.212 |
| 459 | 0.123 / 0.108 / 0.108 | 0.231 / 0.177 / 0.200 | 0.277 / 0.223 / 0.262 |
| 548 | 0.118 / 0.101 / 0.088 | 0.193 / 0.151 / 0.172 | 0.277 / 0.223 / 0.244 |

\* 246.tiff is not an equal-budget comparison — see below. Kept in the table for completeness but
excluded from the headline win/loss framing in "Conclusion."

## Sign test (ROI-level, two-sided, headline per the Sep-2 plan's methodology)

| budget | comparison | wins / losses / ties | n=10 literal: p | tie-dropped n_eff: k/n_eff, p |
|---|---|---|---|---|
| 100 | headroom vs. binary | 3 / 5 / 2 | 0.3438 | 3/8, 0.7266 |
| 100 | multiotsu vs. binary | 3 / 6 / 1 | 0.3438 | 3/9, 0.5078 |
| 250 | headroom vs. binary | 3 / 6 / 1 | 0.3438 | 3/9, 0.5078 |
| 250 | multiotsu vs. binary | 3 / 6 / 1 | 0.3438 | 3/9, 0.5078 |
| 500 | headroom vs. binary | 2 / 7 / 1 | **0.1094** | 2/9, **0.1797** |
| 500 | multiotsu vs. binary | 3 / 6 / 1 | 0.3438 | 3/9, 0.5078 |

At every budget and for both arms, `binary` wins on more ROIs than it loses. **No comparison
reaches p<0.05.** Note on the budget=500/headroom row specifically: the literal `n=10` test gives
`binomtest(2, n=10)` = 0.1094; dropping the one tied ROI gives `binomtest(2, n=9)` = 0.1797 —
*weaker*, not stronger. Dropping ties reduces `n` and therefore reduces the best achievable power
(the n=9 floor, unanimous 9/0, is p=0.0039 vs. n=10's 10/0 floor of p=0.0020) — tie-dropping never
helps an argument for significance here, it only costs power.

Power ceiling at `n=10` (`scipy.stats.binomtest`, verified not assumed): the smallest two-sided p
achievable short of unanimous 10/0 (p=0.0020) is **9/10 → p=0.0215**. With the observed win counts
(2-3 out of 10), nothing here was ever going to clear that bar regardless of which way the effect
actually points.

### 246.tiff is not an equal-budget comparison, and drives most of the apparent "win"

All 5 `binary`-arm seeds on 246.tiff deliver far fewer than 500 total candidates (n_detections:
**2, 44, 54, 88, 99**) — `binary`'s occasionally-merged, multi-nucleus templates on this dense
domain are specific enough that very few windows clear `score_threshold=0.5` at all. `headroom`
delivers 362–2195 candidates there; `multiotsu` delivers 874–**8662**. So "`recall_at_budget_500`"
on 246.tiff is really "recall at the end of a 2-to-99-item list" for `binary` against a genuine
top-500-or-more for the other two arms — not the same budget. 246.tiff is also one of only two
ROIs (with 300) where `headroom` ever wins at budget=500, and one of only three (with 245, 300)
where `multiotsu` ever wins there.

**Excluding 246.tiff from the budget=500 headline changes the record substantially**: `headroom`
vs. `binary` becomes **1 win / 7 losses / 1 tie** (n=9 ROIs); `multiotsu` becomes 2 wins / 6 losses
/ 1 tie. An even more one-sided result than the already-negative 10-ROI number.

### 202.tiff: a concrete illustration of the collateral mechanism

Seed annotation 4534 (picked identically by all three arms on 3 of 202.tiff's 5 seed draws — a
matched-seed, apples-to-apples comparison, not a seed-pool artifact) shows the failure mode
directly:

| arm | tightened template size | n_detections_total | recall@100 / @250 / @500 |
|---|---:|---:|---|
| binary | 37px | 832 | 0.200 / 0.400 / 0.600 |
| headroom@0.15 | 15px | 12,812 | 0.000 / 0.000 / 0.000 |
| multiotsu | 13px | 15,184 | 0.000 / 0.000 / **0.133** |

Raising the threshold shrinks this click's template from 37px down to 15px (headroom) or 13px
(multiotsu) — geometrically "cleaner" (smaller, presumably closer to one nucleus), but a template
that small stops discriminating: the candidate list balloons by ~15-18x (832 → 12,812 / 15,184)
and the true mitotic figures get diluted into a much longer, noisier ranking. Recall collapses
from 0.6 to 0.0 (headroom) or 0.133 (multiotsu) at the same budget. This is exactly the
"collateral shrinkage" the `headroom` docstring in `seed_selection.py` already warns about for
seed-selection geometry — here it is visible in the actual detection output, not just box size.

### Domain composition

The 10 ROIs span 7 tumor-type/scanner domains (verified against `dataset.load_annotations()`'s
`tumor_type` column); three are represented twice (lung: 201/202, lymphosarcoma: 245/246, mast
cell: 300/301), the other four singly (breast/094, neuroendocrine/402, sarcoma/459,
melanoma/548). The threshold-raising arms' wins are concentrated in exactly two of those seven
domains (lymphosarcoma via 246, mast cell via 300) — and even within those pairs the result is
split, not consistent: 245 (same domain as 246) is a near-wash at best for `headroom` and a small
win for `multiotsu`, while 301 (same domain as 300) clearly favors `binary`. No domain is
unanimous for the raised-threshold arms; several (breast/094 partially, lung/201-202,
neuroendocrine/402, sarcoma/459, melanoma/548) favor `binary` outright or tie.

## Conclusion

This closes the gap `2026-09-03-bbox-threshold-sweep.md` explicitly left open ("downstream
template-matching/discrimination quality... remains completely untested"), and the answer tempers
that sweep's geometry-only optimism rather than confirming it: **raising the Otsu foreground
threshold — via `headroom@0.15` or `multiotsu` — does not improve end-to-end detection on this ROI
set, and directionally hurts it**, though n=10 keeps every comparison short of formal significance
(best case p=0.109 literal / 0.180 tie-dropped, both far from 0.05).

The mechanism is a real, identifiable tension, not noise: shrinking a bbox-tightened component can
be geometrically "correct" — smaller, more solid, visually closer to a single nucleus, exactly what
the threshold sweep measured and found improving — while simultaneously being template-matching
*worse*, because an over-tight template is less specific, not more. A smaller crop matches more
locations above a fixed correlation threshold, not fewer, which floods the ranked candidate list
(the 202.tiff case: 832 → 12,812–15,184 candidates from a 37px → 13-15px shrink) and dilutes true
positives deeper into a noisier ranking. The geometry sweep could not see this because it never
ran a search — it only checked whether the accepted box looked like one nucleus.

Two things make the case for "don't ship this" stronger than the raw 10-ROI sign test alone
suggests, both already detailed above and worth restating together: (1) the seed-pool mismatch
tracked throughout this run shows every apparent win for the raised-threshold arms sitting in the
*least* seed-controlled ROIs (0-1/5 same-seed draws), while every well-controlled ROI (3-5/5 same
draws) favors `binary` or ties it exactly — the confound inflates the raised-threshold arms'
apparent competitiveness, it does not muddy a real advantage; and (2) one of the two ROIs driving
`headroom`'s already-weak win record (246.tiff) is not even an equal-budget comparison, because
`binary`'s occasionally-merged templates there are so specific that very few candidates clear the
score floor at all — remove it and `headroom`'s record at budget=500 drops to 1 win / 7 losses out
of 9 ROIs.

**Not decided here**: this does not prove raising the threshold never helps end-to-end (n=10 is
underpowered for anything short of a near-unanimous effect, per the Sep-2 plan's own power
analysis), and it does not revisit whether `binary`'s occasional multi-nucleus merges are
themselves harmless — 246.tiff's own detection counts (2-99 candidates from a 115-mitotic-figure
ROI) suggest they might not be. What this run does establish is that the specific fix tested here
— raising the foreground threshold to un-merge components — trades a geometry problem for a
discrimination problem, and on the evidence collected so far the trade is not favorable.

Raw data: `results/bbox_headroom_end_to_end.csv` (150 rows). Reproducible via
`/Users/mohinianand/anaconda3/bin/python3 bbox_headroom_end_to_end.py` (the repo's own `python3`
has a broken numpy/cv2 ABI — do not use it). Re-run and diffed against 4 (ROI, seed, arm) rows
spanning all 3 arms and 4 different ROIs as part of this write-up's verification; all matched the
saved CSV to full floating-point precision.
