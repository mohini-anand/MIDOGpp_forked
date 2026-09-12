# Decisions pending independent verification

`DECISIONS.md` records choices I have made and intend to keep. This file holds decisions
that were moved out of it for a narrower reason than being wrong or superseded: the entry's
own evidence has not yet been rigorously, independently re-verified — often because it was a
simple initial experiment, run once, without the re-derivation or adversarial check a
standing decision should have.

An entry here is not a weaker decision than one in `DECISIONS.md` — it is one whose
confidence label is honest. `DECISIONS.md` keeps a two-line stub at the entry's original
number pointing here, exactly the way it already points to `D8_TEMPLATE_ANCHOR.md` for D8.
Moving an entry back requires the same thing any `DECISIONS.md` entry requires: a dated
re-derivation, added to this file as an amendment, strong enough that the entry can carry its
original number back into `DECISIONS.md` unchanged.

Started 2026-09-12.

---

## D1 — Use `TM_CCOEFF` for template matching, not `TM_CCOEFF_NORMED`

**Moved from `DECISIONS.md` 2026-09-12.** Not reversed, not superseded — the evidence below
has not been rigorously, independently re-verified. It was a simple initial experiment and
should be read that way until it is.

**Date:** 2026-09-04

### The decision

The similarity measure used to compare the clicked cell's template against every position
in the image is `cv2.TM_CCOEFF`. Previously it was `cv2.TM_CCOEFF_NORMED`. In code this is
`FSConfig.tm_method` in `midog_utils/find_and_suppress.py`.

### Why, in plain words

Both of these compare a template to a patch by multiplying them together pixel by pixel and
summing. Both first subtract the average brightness from each, so **neither one cares how
bright the patch is overall** — that part is the same.

The difference is what happens next. `TM_CCOEFF_NORMED` then **divides by how much the
template and the patch vary internally**. That division is what makes the score land in
[−1, 1] and it makes the score *invariant to contrast*: a washed-out, low-contrast copy of
the clicked cell scores the same as a dark, high-contrast one. `TM_CCOEFF` does not divide,
so a patch with strong internal contrast — a dense dark chromatin mass sitting in paler
cytoplasm — produces a bigger number than a faint one of the same shape.

**Contrast is exactly what distinguishes a mitotic figure.** A cell in mitosis has its
chromosomes condensed into a dense dark mass. `results/morph_diag_bhattacharyya.csv` found
mean intensity to be the strongest single separator of mitotic from ordinary nuclei in all
seven tumour types. So the normalised version was deliberately discarding the best signal
available, by design, and that explains a run of six earlier experiments — score threshold,
multi-scale, augmentation, colour channel, box tightening, NMS ordering — that all came back
neutral. Each was a different way of cutting the same blinded ranking.

**Precision worth keeping straight, because it is easy to get wrong a year from now:**
`TM_CCOEFF` recovers **contrast**, not **absolute darkness**. Both measures subtract the
mean, so neither sees absolute brightness. Absolute darkness is a *separate* statistic —
that is what "chromatin density" measures, and it is computed on the candidate patch, not by
the matcher. The two are complementary, not duplicates.

### The evidence, and how strong it is

`TM_CCOEFF` beats `TM_CCOEFF_NORMED` on `read_50` (candidates read to find half the mitotic
figures) on **49 of 49 decision-grade cells and 10 of 10 ROIs, p = 0.00195, median ratio
4.32×** — 148 candidates against 710. Independently re-derived in
`Research Logs/2026-09-02-next-steps-plan.md` (U0), and the obvious objections were checked
and all fail in the claim's favour: there are no NaNs in either arm, and the normalised arm
actually has *higher* image coverage while still reading 4.3× deeper, which is the opposite
of what a tiling artefact would produce.

**Two limits on that number that I should not forget:**

1. **The 4.3× applies when the ranking key is the match score itself.** If candidates are
   re-ranked afterwards by chromatin density, the two methods are a dead heat (U1 in the same
   document). So this is a good default and a real improvement to the fallback configuration,
   but it is not a 4.3× product win.
2. **It is a `read_50` number, and D4 below moves the reporting metric away from `read_50`.**
   The audit in `Research Logs/2026-09-03-fp-reduction-framing.md` §3b showed `read_50` is the
   depth at which method choices matter least. **This comparison is owed a re-derivation at
   recall@K before it gets quoted again.** The decision still stands — it is one constant and
   nothing argues for the normalisation — but the size of the win is unestablished at the new
   metric.

### What it costs — the trap to watch for

**The score threshold no longer means what it used to mean.** `TM_CCOEFF_NORMED` always
returns a value in [−1, 1], so a fixed floor like `score_threshold = 0.5` was meaningful.
`TM_CCOEFF` is unnormalised, so its range depends on the template's contrast and the image's
contrast — it varies from image to image and from click to click. A fixed number like 0.5 is
meaningless and, depending on the channel, either admits everything or nothing.

**So the threshold must be set in robust-z units**, i.e. "this many robust standard deviations
above this image's own background noise level", using `template_match.robust_stats`. That is
what `tm_variant_sweep.py` already does. This is not optional; it is what makes the decision
safe. `FSConfig.tm_method`'s own docstring says the same thing.

### What would change my mind

A recall@K comparison, worst-ROI and worst-click, in which `TM_CCOEFF_NORMED` is equal or
better. Given U0's 10/10 result I do not expect this, but the metric change is real and the
re-derivation has not been done.

### Amendment, 2026-09-04 — where this decision actually lives in code, and one consequence

**It is not in `FSConfig.tm_method`.** That field still defaults to `cv2.TM_CCOEFF_NORMED`
(`midog_utils/find_and_suppress.py`). D1 is implemented at the *call sites*: `f1_seed_sweep.py`,
`f5_nms_radius_ablation.py` and `tp_fp_feature_extract.py` each set `METHOD = cv2.TM_CCOEFF` and
pass it to `template_match.fused_response` directly, bypassing the config. Every experiment since
2026-09-04 therefore runs `TM_CCOEFF`; `find_and_suppress()` itself does not.

**The default is deliberately left alone rather than flipped**, for the reason this entry's own
"What it costs" section gives: `FSConfig.score_threshold` is a constant in normalised units, and
the older callers that still use it (`od_experiment.py`, `nms_ordering_probe.py`,
`click_rank_stage3.py`) would silently change behaviour — a fixed `0.5` floor means nothing under
an unnormalised score. Flipping the default is safe only once those callers set their floor in
robust-z units. Until then the discrepancy is documented, not fixed.

**One consequence for D1's evidence, folded in from D5.** The "1." limit above — that the 4.3×
becomes a dead heat once candidates are re-ranked by chromatin density — is now weaker than it
reads, because the chromatin re-ranking it defers to is no longer the reference ranker (D5). The
comparison D1 is owed at recall@K should be run on the match score alone.

### Amendment, 2026-09-09 — the precision@K re-derivation is done; D1 holds

The "What would change my mind" line above — "the re-derivation has not been done" — is stale.
`precision_at_k_budgets_14roi_normed/precision_at_k_budgets_14roi_normed.ipynb` re-runs the same
14 ROIs and the same click per ROI as `precision_at_k_budgets_14roi.ipynb`, `TM_CCOEFF_NORMED` in
place of `TM_CCOEFF`, D2/D3/D5/D7 held fixed, and reads precision@K (K = 10/20/30/50) in place of
`read_50`. Its Table C compares directly against the `TM_CCOEFF` notebook's already-committed
output, with `seed_ann_id`/`base_size`/`n_gt_mitotic` asserted identical between the two runs
first, so the difference below is attributable to the method alone.

`TM_CCOEFF` wins 28/28 domain x budget cells and 54/56 ROI x budget cells — the 2 non-wins are
exact ties, zero reversals either way, the same completeness as the original 49/49. Magnitude is
comparable to the original 4.32×, not larger: pooled precision ratio runs 3.55× (K=10) to 4.42×
(K=50); the statistic actually comparable to the original median-of-ratios is 4.125× (median over
the 36/56 cells with a defined ratio — the other 20 have `precision_normed == 0` against a
nonzero `TM_CCOEFF` precision, checked to confirm each is a `TM_CCOEFF` win and not a 0/0 tie, so
4.125× understates the win rather than overstates it).

D1 stands, now confirmed at the metric family that replaced `read_50` (D4), on the match score
alone, per the scope D5 attached to this obligation.

### Why this is here despite two amendments and a 2026-09-09 re-derivation

Both amendments and the re-derivation are re-runs of the *same* comparison (`TM_CCOEFF` vs.
`TM_CCOEFF_NORMED` on this repo's own pipeline, by the same author, without an adversarial
check of the premise itself — e.g. whether `read_50`/precision@K is the right unit of
comparison, whether the six-earlier-experiments explanation is the actual cause or a
post-hoc story, whether `morph_diag_bhattacharyya.csv`'s "mean intensity is the strongest
separator" finding transfers to what the matcher is doing). None of that has had independent
scrutiny. That is the bar for staying in `DECISIONS.md`; this entry has not cleared it yet.
