# Pre-registration: F1, does largest-CC template tightening have a resolvable effect?

Date: 2026-09-04. **Written and committed before `f1_seed_sweep.py` was run.**
Follow-up F1 of `Research Logs/2026-09-03-tm-axis-sweep-edit-plan.md`.

## Why pre-register this one

The output CSV will contain 2 axes x 6 z x 8 budgets of `recall_at_budget` plus 2 x 6 x 6
read-depth cells -- roughly **170 plausible "primary" tests**. Choosing after seeing them is a
far larger forking path than the seven per-domain sign tests this replaces. The audit this
follows exists partly because earlier work read per-domain signs off single-seed numbers; a
follow-up that picks its own test post hoc would repeat the error in a more sophisticated form.

## Hypothesis, and its honest provenance

**H1.** The paired effect of largest-CC tightening on ranking quality is *monotone in template
size*: the smaller the tightened template relative to 51 px, the worse the ranking at a fixed
reading budget.

**Mechanism** (stated independently of the data): the largest-CC notebook's own account is that
a smaller template produces a broader, less-selective correlation response, so more of the map
clears the same per-map `z` cut. A longer, less selective list ranks worse at fixed K. This is a
real prediction and not a tautology -- audit Finding 12 measured that the mechanism holds in
only some domains (`301.tiff`'s deep pool is unchanged at 51 -> 41 px, while `459.tiff` grows
1.39x), so it could fail.

**Circularity disclosure.** H1 was *suggested* by the seed-0 cells, which are 7 of the 35 this
run will produce, and which I have already seen: Spearman rho(`area_ratio`, paired delta in
`recall_at_budget` at z=1.0, K=250) = **+0.764 (p=0.046)** on `chromatin_od` and **+0.655
(p=0.111)** on `tm_score`, with the sign reversing from negative at `area_ratio <= 0.42` to
positive at 0.65. That is mild circularity and it is handled by pre-committing to the
seed-0-excluded version as decisive (below), not by pretending I had not looked.

## Why not a sign test

The obvious framing -- per-domain paired sign tests -- fails for two independent reasons,
both measured before the run:

1. **10 of 35 cells are null** (`tightened == 51`; the arms are the same computation and the
   delta is exactly 0 by construction), and they are not evenly spread: `201.tiff` yields **one**
   informative pair, `246`/`402` three, `301`/`459` four, `094`/`548` five. A five-seed sign test
   is attainable in **2 of 7** domains and impossible for `201.tiff`.
2. **The effect reverses sign with dose**, so a pooled sign or signed-rank test would see the
   negative and positive deltas cancel and return "no effect" while a strong monotone
   relationship is present. It is *mis-specified*, not merely underpowered.

## Primary test -- one, fixed now

> **Within-cell sign-flip permutation test (10,000 draws, two-sided) of Spearman rho between
> `area_ratio` and the paired delta (`largest_cc` - `base51`) in `recall_at_budget`, at
> `z = 1.0`, `K = 250`, axis = `chromatin_od`, over all 35 cells.**

* **Sign-flip, not label permutation**: independently flipping the sign of each paired delta is
  the exact randomization null for a paired design, and each cell is its own stratum, so domain
  clustering is handled by construction rather than by a modelling assumption.
* `area_ratio = (tightened_size / 51)^2` -- area, not side, because response magnitude scales
  with template pixel count (`template_match.fused_response`'s own docstring).
* **Null cells are included**, at `area_ratio = 1.0, delta = 0`. They are genuine observations of
  the treatment at dose 1, they anchor the boundary `delta -> 0` as `area_ratio -> 1`, and they
  contribute no permutation entropy either way.
* `z = 1.0` is `CURRENT_Z`, the operating point both TP/FP notebooks settled on -- not chosen here.
* `axis = chromatin_od` is the ranker the repo ships (commit `7c3af93`, "Rank detections by
  chromatin density instead of correlation score").
* `K = 250` is chosen from measured properties, not taste: at `z = 1.0` the `base51` arm's
  recall@250 spans 0.41-0.82 with **0 of 7 domains saturated**; K = 25/50 give discretised
  near-zero deltas that collapse power, and K >= 1000 begins to saturate. `n_gt_mitotic` after
  seed removal spans 17 (`201.tiff`) to 238 (`548.tiff`), so recall at large K has a
  domain-dependent ceiling.
* **`recall_at_budget`, not `read_95`.** `compare.py` designates it primary, and it is
  length-invariant while `n_detections >> K`. Both arms are cut at the same per-map `z` but
  **not** at a matched list length, so a `read_95` delta would re-report the candidate-volume
  difference under a new name -- audit Finding 2 all over again.

## Decision rule, fixed now

| outcome | conclusion |
|---|---|
| primary p < 0.05 **and** the seed-0-excluded version also p < 0.05 | H1 supported: tightening's effect is real and dose-dependent. lcc's per-domain signs were noise, but a monotone structure underlies them. |
| primary p < 0.05, seed-0-excluded p >= 0.05 | **Not supported.** Treat as generated by the cells that suggested it. The seed-0-excluded version is decisive. |
| primary p >= 0.05 | No resolvable dose-response at n=35. Combined with the effect being smaller than the within-arm seed SD, the conclusion is that **tightening is not a decision-relevant lever** and deliverable (d) should say so. |

A null is a real result here and closes deliverable (d) honestly. I am not entitled to go
looking for a different (axis, z, K) afterwards.

## Secondary and robustness, all reported regardless of the primary outcome

1. **Seed-0-excluded primary** (28 cells) -- decisive per the table above.
2. Same test on the `tm_score` axis.
3. **Within-arm seed SD under this exact configuration**, reported beside every paired delta.
   The audit could only compare lcc's deltas against `tm_ccoeff_headtohead_seed_variance.csv`,
   whose SDs come from a different pipeline and carry that caveat explicitly. This run retires
   it. If the paired effect is significant but small against this SD, the honest conclusion is
   still "not a decision-relevant lever".
4. **The full (axis, z, K) grid of dose-response slopes, as an explicitly labelled exploratory
   heatmap.** Cheaper and more honest than a post-hoc winner, since the CSV holds every cell.
5. `read_95` is **right-censored, not missing**: `compare._read_depths` returns NaN when 95%
   recall is unreachable, which means "deeper than the entire list". Measured on the two saved
   CSVs over the 84 paired cells in this z grid: 76 both-finite, **4 where `base51` is NaN and
   `largest_cc` is finite, 0 the reverse**. Dropping NaN pairs would discard exactly the cells
   where `base51` is strictly worse, biasing toward "tightening doesn't help" in the only z
   range where lcc claimed a benefit -- audit Finding 4 re-entering the follow-up to Finding 5.
   NaN is therefore ranked as greater than any finite value, never dropped, and
   `read95_unreachable_*` is recorded per cell.
6. Paired `n_detections` and `coverage_frac` deltas printed beside every read-depth number, so
   a volume change can never be read as a quality change.

## What is *not* a limitation

Seeds are drawn uniformly **with replacement**, and two streams collide (`201.tiff` draws ann
4457 at seeds 0 and 1; `094.tiff` draws ann 2494 at seeds 2 and 4). Those are still iid draws --
`invariants.py`'s `check_distinct_seeds` says a collision "is expected, not a defect" -- so
`Var(mean) = sigma^2 / 5` already prices it and **both draws are kept**. Dropping them would
condition on the realized sample.

## Real limitations

* 7 domains, one ROI each; not exchangeable. Every pooled statement is about this sample.
* `201.tiff`'s entire border-filtered agreement pool is 7 candidates.
* The match padding is *not* identical between arms (`PAD` = 25 for 51 px, 12 for 25 px, 9 for
  19 px), so within ~25 px of the ROI edge the arms score against `BORDER_REPLICATE`
  fabrication of different widths and the pairing is not clean there.
* `tightened_size` spans 19-47 across cells; 19 px is smaller than anything the largest-CC
  notebook ran (its seed-0 range was 25-51), so the low-dose end extrapolates beyond it.

## Gates the run must pass before any of this is read

1. At `seed_index = 0`, `base51` reproduces `results/tm_ccoeff_threshold_axis_sweep_v2.csv` and
   `largest_cc` reproduces `..._largest_cc.csv` on 11 columns x 96 (arm, z, budget) cells --
   1,056 values per domain-arm, checked *inside* the domain loop so a divergence aborts on the
   first domain. Self-tested against the real CSVs with seven injected faults, all caught.
2. In every null cell the two arms must be exactly equal on every column.
3. The one-match-many-`z` shortcut is re-verified per domain at `seed_index = 0`.

Failure of any gate invalidates the run; the output keeps a `.partial` suffix and a
`run_complete=False` marker so nothing downstream can mistake it for a finished result.
