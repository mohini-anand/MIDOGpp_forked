# D8 — The seed template is cut from the accepted component, not from the click

**Date:** 2026-09-10 · **Status:** current · **Supersedes:** D8 (2026-09-09) and its amendment of
the same date, both moved out of `DECISIONS.md` and reproduced verbatim at the bottom of this file.

---

## The decision

A seed's search template is fixed by three numbers: a side length, and a centre `(x, y)`. All
three now come from the Otsu component the click validated.

**1 · The click gates, and only gates.**
`seed_selection.tighten_box_otsu` thresholds the 51 px window around the click and takes the
connected component **the click's own rounded pixel falls inside** (`center_tolerance = 0`). That
component must clear `min_area = 50`, `max_area_frac = 0.85`, `min_solidity = 0.5`. If the click
is not inside such a component, the seed is **refused** and redrawn on the same RNG stream. This
is unchanged and is not up for revision here.

> **Not "the largest component".** The rule reads `labels[cy, cx]`; it never asks whether that
> component is the biggest thing in the window. On the 14 ROIs of `images/extra_valid` the two
> rules pick the same component on 878 of 896 accepted seeds (98.0%) — but where they differ,
> "largest" is the fallback this gate exists to refuse. `BBOX_TUNING.md` calls it *"a weak
> heuristic … in a crowded crop it can grab an unrelated, larger structure"*.

**2 · The size is the component's longer side.**
`base_size = _odd(max(y1 - y0, x1 - x0))` from the accepted component's bounding box. Unchanged.

**3 · The centre is the component's bounding box, in pixel coordinates.**

```
center_x = ix - half + (x0 + x1 - 1) / 2.0
center_y = iy - half + (y0 + y1 - 1) / 2.0
```

where `(ix, iy) = (round(cx), round(cy))` and `half = otsu_window // 2`. This replaces the click,
and it replaces the previous `(x0 + x1) / 2.0`.

---

## Why the click gates but does not centre

The click is a **point annotation**, not a segmentation. `BBOX_TUNING.md` establishes that every
MIDOG++ box is a synthetic `point ± 25` square and is "rarely centred on the actual chromatin
mass". The click's job is to say *which object is mitotic* — a label no image statistic can
supply. It is not evidence about where that object's centre is.

`tighten_box_otsu` already measures the object's real extent and passes it through three sanity
gates before anything downstream sees it. Using that measurement for size but not for position
throws away information the gate already paid for.

Measured cost of not using it: a click-centred template contains the whole validated component on
only **28.2%** of accepted seeds, and on `246.tiff` it misses about **a fifth** of the object the
gate had just certified.

---

## The half-pixel correction

`(x0 + x1) / 2` was the wrong formula for this consumer, and the error is systematic.

A skimage bbox is **half-open**: `x0 = 6, x1 = 51` covers columns 6…50. There are two honest
centres — the centre of the *area* it spans, `(6+51)/2 = 28.5`, which lands on the seam between
two pixels; and the centre of the *pixels* it covers, `(6+50)/2 = 28.0`, which lands on a pixel.

`template_match.read_padded_patch` computes `int(round(cx))` and then slices `± patch_size//2`.
It names a **pixel**; it has no way to centre a crop on a seam. So the anchor must be computed
the pixel way. The old formula was 0.5 px too far down and right, always — and because Python
rounds `.5` ties toward even, that half pixel becomes a **whole** pixel on about half of all
boxes, depending only on the parity of the box's true centre.

**What the correction buys is a guarantee, not a tendency.** Since
`base_size = _odd(max(h, w)) ≥ max(h, w)`, a square of that side centred on the box's own pixel
centre must contain the entire box. Every component pixel is inside the template. Under the old
formula that guarantee is void, and it fails hardest when `max(h, w)` is already odd — then
`_odd` rounds nothing up, the crop has zero slack, and one pixel of drift clips a row or column
off the object.

---

## The evidence

All figures on the 993 unanimous, border-filtered mitotic annotations of `images/extra_valid`;
`tighten_box_otsu` accepts **896 (90.2%)**, and recentring only ever applies to those.

**Geometry — how much of the validated component the template actually contains:**

| anchor | contains the whole component | worst case | crops the component |
|---|---|---|---|
| the click *(previous production)* | 28.2% | 0.556 | 643 / 896 (71.8%) |
| `(x0+x1)/2` *(previous D8)* | 60.5% | 0.932 | 354 / 896 (39.5%) |
| **`(x0+x1-1)/2` *(this decision)*** | **100.0%** | **1.000** | **0 / 896** |

Checked exhaustively over every box shape from 1×1 to 51×51 at four positions — 41,616
combinations: `(x0+x1)/2` clips on 25.6% of them, `(x0+x1-1)/2` on none.

**How far the template moves:** median **3.10 px** from the click, 95th percentile 8.14 px, max
13.87 px, against templates of roughly 19–51 px. Not a small correction.

**It reaches the seeds that carry every committed result.** On the 14 production seeds
(`seed_index = 0`), the old formula crops the component on **9 of 14**; the click-centred rule
crops it on **10 of 14**; this decision's anchor on **none**.

Derivation and every figure above: `pipeline_debug_visuals/template_anchor_halfpixel_fix.ipynb`,
which verifies each claim against the live `seed_selection` functions before using it and
reproduces the 896/993 accept rate as a gate.

---

## What it costs

1. **Self-hit removal moves with the centre.** The template's self-correlation peak lands at the
   returned centre, not at the click. Every caller doing self-hit or seed-annulus removal must
   reference the returned centre. `find_and_suppress.find_and_suppress` and the notebook-inline
   `suppress` helpers currently use the click and must be updated.
2. **Border safety moves with the centre.** `border_filter`'s margin (`patch_size // 2` = 36 px)
   is sized for a click-centred read; the returned centre can sit up to ~14 px further toward the
   ROI edge. Callers must re-check that the full rotation-safe patch is readable at the recentred
   point and refuse the seed if not.
3. **Ground truth does not move.** `(cx, cy)` stays the click for `gt_eval` exclusion and every
   MIDOG match-radius computation. Only where the template is *built from* changes; what counts
   as finding a mitotic figure does not.
4. **Committed numbers are tied to the old formula.** The four `production_seed_precision_at_k*`
   notebooks implement `(x0+x1)/2`; 12 of their 14 seeds move under this decision, so their
   outputs are a record of the old anchor and not comparable across the change. Every other
   committed `precision_at_k_*` result used click-centred templates and is likewise a record of
   the previous production path.
5. **Two docstrings become stale on adoption.** `seed_selection.tightened_base_size` says it is
   "now the production path" and `tightened_template_box` says it is "superseded" — both are
   assertions of the reversed 2026-09-09 amendment.

---

## What this decision does **not** claim

**That it improves precision@K or recall@K.** It does not. The evidence above is entirely
mechanism: what the template contains, and a formula that is provably right where the previous
one was provably wrong. No outcome-level comparison of this anchor against the click has been
run at the bar `DECISIONS.md` D5 sets — 5 seeds × 14 ROIs, paired delta clustered at the ROI.

The one outcome-level read that exists is single-seed and mixed. On `245.tiff` and `403.tiff`
(8-augmentation bank, everything else pinned, `SELF_HIT_RADIUS = 10.0` in all arms), this
decision's anchor lands between the click and the old formula: `tp@50` of 2 against the click's 3
on 245, and `tp@30` of 12 against the click's 11 on 403. That is ranking churn in a pool of
~19,000 candidates packed into a narrow score band, not a result.

**The honest summary: this entry corrects an implementation and adopts it on mechanism. Whether
aligning the template to the validated object helps the product is open.**

---

## What would change my mind

A paired sweep — **5 seeds × 14 ROIs, click-centred against this anchor, identical size, identical
`SELF_HIT_RADIUS`, paired delta clustered at the ROI** — showing this anchor neutral or harmful at
precision@K or recall@K. That is D5's bar for this class of decision and it has never been run for
centring. It is the natural next piece of work and the one thing this entry is exposed to.

Nothing would restore `(x0 + x1) / 2`. That formula is wrong for a consumer that rounds its input
to a pixel index, independently of whether recentring is adopted at all.

---

## Implementation status

**Not yet in code.** `seed_selection.tightened_template_box` still computes `(x0 + x1) / 2`. The
corrected anchor exists only as `tightened_template_box_fixed` inside
`pipeline_debug_visuals/template_anchor_halfpixel_fix.ipynb`. Adopting this decision requires:

- correcting the two centre lines in `seed_selection.tightened_template_box`;
- updating that function's and `tightened_base_size`'s docstrings, which currently assert the
  reversed amendment;
- routing `find_and_suppress` and the notebook-inline `suppress` helpers to the returned centre
  (cost 1) and adding the recentred border re-check (cost 2).

Recording this decision without the code change repeats the gap the superseded entry's own item 5
flagged. It is named here so it cannot be forgotten, not so it can be deferred.

---
---

# Superseded history

Kept verbatim because `DECISIONS.md`'s own policy is that the history of a decision is as useful
as the decision. Neither section below is current; both are reproduced exactly as they stood in
`DECISIONS.md` before being moved here.

## [SUPERSEDED] D8 — Bbox tightening recentres the template on the detected nucleus, not on the click

**Date:** 2026-09-09

### The decision

When `seed_selection.tighten_box_otsu` accepts a component for a seed click (the click's own
pixel, or a pixel within `center_tolerance`, lands inside it — the gate itself is unchanged by
this entry), the template built from that seed is now centred on **that component's own
bounding-box centre**, not on the raw click coordinate. Only the *size* was corrected before
(`tightened_base_size`); this decision corrects the *position* too.

In code: `seed_selection.tightened_template_box`, new, returns `(base_size, center_x,
center_y)` in full-image coordinates. `tightened_base_size` (size-only, click-centred) is kept,
unchanged, for callers that explicitly want it and for reproducing every already-committed
result in this repo — new work should call `tightened_template_box`.

**What this does not do.** It does not relax the containment gate. A component that the click
lands outside of is still refused as a seed (`None`), exactly as today — the template is never
built from a different, unannotated structure just because it happens to be the largest thing
in the window. That looser rule was considered and explicitly rejected in the course of
reaching this decision (see "Why, in plain words").

### Why, in plain words

MIDOG++'s ground-truth box is a synthetic `point ± 25` square (`BBOX_TUNING.md`), not a
segmentation — the click marks the object, it does not mark its centre. `tighten_box_otsu`
already does the work of finding that object's real extent and passing it through three sanity
gates (`min_area`, `max_area_frac`, `min_solidity`) before anything downstream sees it. Using
that measurement for size but not for position was discarding information the gate had already
paid for and validated.

This was visible directly on one example first (`precision_at_k_budgets_14roi_chromatin/
lcc_offset_demo.ipynb`): 094.tiff ann 2582's largest Otsu component sits 20.5 px from the click
and the click's own pixel is not even inside it. That example is illustrative, not evidentiary
for *this* decision, though — under `tighten_box_otsu`'s containment gate it would be refused
as a seed outright, not recentred onto a neighbouring object. It is what motivated checking the
gated population, not the measurement the decision rests on.

The measurement the decision rests on is the gated population: **it also loosened, then
tightened its own recentring rule during this session** — the first framing considered was "use
the single largest Otsu component in the window, click or no click," which would have recentred
onto exactly this kind of unannotated neighbour whenever the click missed its own object. That
was rejected in favour of keeping the containment gate exactly as it stands today and only
recentring within it.

### The evidence, and how strong it is

`lcc_offset_demo.ipynb`, same cell that produced the number above: 993 unanimous,
border-filtered mitotic annotations across the 14 ROIs of `images/extra_valid`.
`tighten_box_otsu` accepts **896/993 (90.2%)** under its existing containment gate — recentring
only ever applies to that 90.2%; the other 9.8% behave exactly as before (refused as a seed).

Among the 896 accepted, the offset between the raw click and the accepted component's own bbox
centre:

| stat | value |
|---|---|
| median | 3.10 px |
| 75th pct | 4.74 px |
| 90th pct | 6.76 px |
| 95th pct | 8.14 px |
| 99th pct | 10.63 px |
| max | 13.87 px |
| fraction > 5 px | 21.5% |

F1 (`Research Logs/2026-09-04-f1-results.md`, audited in
`Research Logs/2026-09-08-f1-largest-cc-vs-51px-audit.md`) found tightening the template down
to roughly 36 px helps and tightening further costs more than it gains — the same audit's
finding F1.9 states plainly that F1 "answers 'does shrinking the template help?', not 'does
fitting the template to the nucleus help?'" and that recentring was a documented, untested
scoping choice, not a measured-and-rejected one. A 3–14 px centring error is not small against
a ~19–45 px template (the range seen in the 30-annotation smoke test run for this entry).

**This is mechanism-level evidence, not outcome-level evidence, and that is a real gap.**
Nothing in this entry demonstrates a recall@K or precision@K improvement from recentring — only
that the correction moves the template by a measured, non-trivial amount on measured, real
data. D5's bar for promoting a chromatin ranker (5-seed sweep, paired delta clustered at the
ROI) has not been applied here, and this decision does not claim that bar is met.

### What it costs

1. **Self-hit removal moves with it.** The template's self-correlation peak now lands at the
   recentred point, not at the click. Any caller adopting `tightened_template_box` must use the
   returned centre — not `(cx, cy)` — as the reference for self-hit and seed-annulus removal
   (`find_and_suppress.find_and_suppress`, and every notebook's inline `suppress`, currently use
   the click for this and have not been updated).
2. **Border safety moves with it.** `border_filter`'s margin (`FSConfig.patch_size // 2` = 36 px)
   is sized for a click-centred read. The recentred point can sit up to ~14 px further toward
   the ROI edge (measured 99th percentile 10.6 px, max 13.9 px) than the click it was filtered
   on. No caller in this repo currently re-checks or widens for that; one adopting
   `tightened_template_box` near an ROI border must.
3. **Ground truth does not move, anywhere.** `(cx, cy)` stays the click for `gt_eval` exclusion
   and every MIDOG match-radius computation. Only where the search template is built from
   changes; what counts as finding a mitotic figure does not.
4. **No existing result is retroactively affected.** Every committed number in this repo —
   `f1_seed_sweep.csv`, every `precision_at_k_*` notebook, `tp_fp_candidate_features.csv` — used
   click-centred templates via `tightened_base_size` or its notebook-inline equivalents. None of
   those are re-run by this entry, and none should be read as evidence for or against it.
5. **Nothing in this repo calls `tightened_template_box` yet.** It exists in `seed_selection.py`,
   smoke-tested (sizes agree exactly with `tightened_base_size` on a 30-annotation sample; 8/9
   accepted seeds actually moved, 0.5–6.8 px). Wiring it into `find_and_suppress`, `pick_seed`,
   and the notebook-inline seed-draw helpers that currently duplicate this logic is follow-up
   work, not done by this entry.

### What would change my mind

A length-matched or seed-swept recall@K / precision@K comparison — recentred vs. click-centred,
identical size, identical everything else — showing recentring is neutral or harmful. That
comparison has not been run; it is the natural next step and the one this entry's "mechanism,
not outcome" gap is waiting on.

## [SUPERSEDED] Amendment, 2026-09-09 — recentring is reversed; keep the containment gate and the size correction, not the position correction

**The position correction is not adopted.** Production seed construction stays
**click-centred**: use the accepted Otsu component's size (as D8 always did) and the existing
containment gate (a click whose own pixel does not land inside its accepted component is refused
as a seed and redrawn — `tighten_box_otsu`'s gate, unchanged since before this entry existed),
but do **not** move the template's centre to that component's own bbox centre. In code: new work
calls `seed_selection.tightened_base_size` (size-only, click-centred, gate enforced) — not
`tightened_template_box` (size **and** position). `tightened_template_box` stays in
`seed_selection.py`, unused by the production path.

**Why.** This entry's own "What would change my mind" asked for exactly one thing: a
length-matched or seed-swept precision@K / recall@K comparison, recentred vs. click-centred,
identical everything else. A single-seed version of that comparison now exists:
`precision_at_k_budgets_14roi/precision_at_k_budgets_14roi_8aug.ipynb` (click-centred,
`largest_cc_box`, no containment gate) against
`production_seed_precision_at_k/production_seed_precision_at_k_8aug.ipynb` (recentred,
`tightened_template_box`, D8 as originally written) — same 14 ROIs, same `seed_index = 0`, same
8-augmentation bank, D1/D2/D3/D5/D7 all held fixed, so centring is the only free variable.

Comparing `precision_pooled` (`results/precision_at_k_14roi_8aug_by_domain.csv` against
`results/precision_at_k_14roi_prodseed_8aug_by_domain.csv`, both already committed) across the 7
domains x 4 budgets = 28 cells: recentring scores **lower in 17/28 cells, higher in 3/28, tied in
8/28** — mean pooled delta ≈ **-0.029** (about 3 points of precision). At the per-ROI level
(`precision_at_k_14roi_8aug_per_roi.csv` vs. `..._prodseed_8aug_per_roi.csv`), the worst
regressions are 301.tiff (mast cell tumour: -0.20 / -0.20 / -0.133 / -0.14 at K = 10/20/30/50) and
529.tiff (melanoma: -0.20 at K = 10). Excluding 403.tiff (below), the mean per-budget delta is
still negative throughout, -0.023 to -0.046.

**One ROI is not a clean read of centring alone.** 403.tiff's two runs disagree on *which*
annotation is the seed — the click-centred run's ungated `largest_cc_box` accepts a click sitting
on a background pixel and borrows an unrelated component's size for it, where the gated method
correctly refuses that draw and retries onto a real, click-containing component
(`production_seed_precision_at_k/bbox_refinement_comparison.ipynb`, Part 3). That is the
containment gate doing its job, not recentring — it is a reason to keep the gate, not evidence
about the position correction either way, and it is why 403.tiff is excluded from the per-ROI
figures above.

**How strong this is.** One seed per ROI — the same limitation both source notebooks flag in
their own scope sections — and it does not meet the bar D5 sets for this class of decision (a
5-seed sweep, paired delta clustered at the ROI). It is a first data point, not a resolved
question. But it is outcome-level evidence, which is exactly what the original entry above
admits it did not have, and the one measurement available now points against recentring, not for
it. Reverting now and re-opening if a fuller sweep disagrees costs less than keeping a correction
whose only support was mechanism-level.

**What is kept from the original decision.** The containment gate — refuse the seed if the click
does not land inside its accepted component, never fall back to "largest component regardless of
the click" — is not new to D8 and is not reversed here; D8 never weakened it and this amendment
does not either. Only the position half reverses: size still comes from the accepted component,
position reverts to the click.

**What this costs.** The four `production_seed_precision_at_k*` notebooks
(`production_seed_precision_at_k.ipynb`, `_8aug`, `_normed`, `_chromatin`) call
`tightened_template_box` and currently implement the position correction this amendment reverses.
Their committed numbers stand as the record of the recentred configuration — the same way item 4
above treated pre-D8 results — but they are not the production path going forward. Re-running
them against `tightened_base_size` is the natural follow-up; it is not required by this amendment.

**What would change my mind, again.** A multi-seed sweep — 5 seeds x 14 ROIs, paired delta
clustered at the ROI, D5's own bar for this class of decision — showing recentring is neutral or
positive at precision@K or recall@K. The single-seed result above is a first data point toward
that bar, not the bar itself.

### Why the 2026-09-09 amendment is superseded rather than simply reversed

Its evidence compared **ungated, click-centred** against **gated, recentred** — two free variables,
not one. It says so itself for 403.tiff and excludes that ROI, but the same confound is present on
the other thirteen, where `largest_cc_box` and the gate can pick different components. And the
recentred arm it measured carried the half-pixel error, which crops part of the component on 9 of
the 14 seeds and is 23% of the median recentring distance, pointing consistently one way.

The amendment's conclusion may still be correct. Its measurement was not clean enough to establish
it, and the anchor it measured is not the anchor above.
