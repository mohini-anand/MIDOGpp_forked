# Decisions

A standing record of choices I have made and intend to keep, with the reasoning written
out in plain words so that I can pick the thread back up after a long gap.

This is **not** a research log. Research logs record what an experiment found; this file
records what I decided to do about it. Each entry says what the decision is, why, what it
replaces, what it costs, and what would make me change my mind. If a future measurement
contradicts an entry here, add a dated amendment underneath it rather than silently editing
the entry — the history of a decision is as useful as the decision.

Started 2026-09-04. Entries are dated and appended, newest last.

---

## The problem, in one paragraph, for future me

A pathologist clicks one mitotic figure in a 2 mm² region of tumour. The tool returns a
list of candidate cells that might also be mitotic, ranked, and the pathologist works down
the list saying yes or no. Two things have to be true at once for this to be worth using:
the list has to **contain** the mitotic figures (so nothing is missed), and it has to be
**short** (so it saves time). A 2 mm² region contains roughly 20,000 nuclei and roughly
100–220 mitotic figures, so this is a needle-in-haystack problem at about 1% prevalence.

---

## D1 — Use `TM_CCOEFF` for template matching, not `TM_CCOEFF_NORMED`

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

---

## D2 — No `tissue_mask`

**Date:** 2026-09-04

### The decision

Stop applying `baselines.tissue_mask` as a pre-filter on the candidate pool. Every pixel of
the ROI is eligible; nothing is excluded for being "not tissue".

### Why, in plain words

`tissue_mask` marks a pixel as tissue if its greyscale value is below 220, the idea being
that empty slide glass under a brightfield scanner is nearly white while stained tissue,
however pale, is not. It exists to stop the detector proposing candidates in blank glass.

**The reason to drop it is that the constant 220 was chosen by looking at the answers.**
The comment in `midog_utils/baselines.py:26-28` says so outright: *"220 was validated across
all 14 downloaded ROIs: it excludes 0 of the 690 mitotic annotations in them."* That is
picking a hyperparameter by checking it against the ground-truth labels of the very images
the system is evaluated on. It is a small leak and it does not bias one method against
another — it affects every arm equally — but it inflates every **absolute** recall and
ceiling number, and those are exactly the numbers a product claim would quote. This is
recorded as M5(a) in `Research Logs/2026-09-02-next-steps-plan.md`.

It is also not quite true. The premise-test audit found it excludes **1 of 726**
annotations, not 0 — `506.tiff` annotation 24249, a category-2 look-alike sitting 13 px from
the image edge with a click pixel at grey 232 against the 220 cut.

### Why removing it is nearly free

**The mask is doing work the ranking already does.** Empty glass has no stain in it. On the
hematoxylin optical-density channel (see D3), no stain means an optical density near zero —
the lowest possible value. So a glass candidate is ranked *last* by construction and can
never appear near the top of the list where the pathologist is reading. The mask removes
candidates that the ranker was already going to bury.

The measured amount at stake is small: `tissue_mask`'s own docstring records that only
**0.7–3.7% of pixels** in these ROIs are above grey 220. These are 2 mm² regions selected as
tumour hotspots — they are near-solid tissue, with very little glass in them.

### What it costs

The candidate pool grows a little, by roughly the glass fraction of each ROI. Any absolute
pool-size numbers computed before this change are not comparable with numbers computed after
it, so **pool size must be reported alongside every depth figure** — a change that halves the
"fraction of pool read" while doubling the pool has made the pathologist's job worse, not
better.

Anything in the codebase that depended on the mask needs re-checking, in particular the
`random_in_tissue` baseline, which without a mask becomes uniform over the whole ROI rather
than over tissue.

### The thing this does *not* do

Removing the leak does not make the old numbers wrong in their comparisons between methods —
it makes the old **absolute** numbers optimistic. If a re-run gives slightly worse absolute
recall than the historical figures, that is the leak being removed, not a regression.

### What would change my mind

A whole-slide setting rather than pre-selected 2 mm² hotspots, where large genuinely empty
areas exist. Then some tissue-finding step earns its place again — but it should be derived
from the image without consulting any labels, or calibrated on images that are never
evaluated on.

---

## D3 — Do not rescale the image after colour deconvolution

**Date:** 2026-09-04

### The decision

After separating out the hematoxylin channel by colour deconvolution, use the raw optical
density values. Do **not** run `cv2.normalize`, a min–max rescale to 0–255, or any percentile
stretch on top of it.

In this codebase that means: use the `hematoxylin_od` channel (`channels.to_hematoxylin_od`,
which is `chromatin.hematoxylin_od`), **not** the `hematoxylin` channel
(`channels.to_hematoxylin`).

### Why, in plain words

**Colour deconvolution** un-mixes an H&E image into its two dyes, so you get a channel that
says "how much hematoxylin is at this pixel". Hematoxylin binds DNA, so that channel is
essentially a chromatin map. Its natural unit is **optical density** — how much light the
stain absorbed. The values are small and non-negative (measured on `002.tiff`: min 0.000,
median 0.0144, 99th percentile 0.0909, max 0.238).

The old `hematoxylin` channel then took those values and stretched them to fill 0–255, using
the image's own 0.5th and 99.5th percentiles as the ends, and **clipped everything outside**.
That last step is the problem. By construction it parks **0.5% of the pixels flat on the 255
ceiling** — and those are precisely the darkest, densest chromatin pixels, which is to say
precisely the mitotic figures. The information that distinguishes "dense" from "extremely
dense" is destroyed in the one place where it matters.

This is not hypothetical. An audit of the chromatin re-ranking work found dozens of the
darkest candidates all scoring **exactly 255**, so their relative order was decided by an
arbitrary sort rather than by the data. Fixing it improved the result
(`Research Logs/2026-09-01-premise-test-audit.md`, Part 1).

### Why this decision and D1 only make sense together

`TM_CCOEFF_NORMED` is invariant to `I → aI + b`, so it literally **cannot tell** whether the
channel was rescaled or not — under the old matcher this choice was invisible and harmless.
`TM_CCOEFF` (D1) reads magnitude, so it *can* tell, and rescaling would delete the signal it
was switched on to capture. `channels.to_hematoxylin_od`'s docstring makes exactly this point.

**So D1 and D3 are one decision in two parts.** Adopting D1 without D3 would switch on a
contrast-sensitive matcher and then feed it a channel with the contrast clipped out of the top
0.5%. That combination is worse than either half.

### What it costs, and the caveat to carry

Optical density is **not calibrated across scanners**. Two images of identical tissue on two
scanners will produce different absolute OD values. So:

* **Scores are comparable within one ROI's ranking, and not between ROIs.** That is fine —
  every metric in this project is a within-ROI rank — but it means no absolute score
  threshold can be shared across images. Thresholds must be per-image and in robust-z units,
  which is the same requirement D1 imposes.
* Numbers pulled from a `hematoxylin`-channel run and a `hematoxylin_od` run are not
  comparable. Label which channel produced any stored result.

There is one place where a *bounded* rescale is still legitimate: inside a small fixed crop
for Otsu thresholding, where the operation is local and the goal is a binary split rather than
a ranking. `bbox_threshold_sweep.py` relies on `cv2.normalize` running before thresholding so
that a crop's maximum is always 255. Do not remove it there without re-reading
`Research Logs/2026-09-03-bbox-threshold-sweep.md`. **The rule is: no global rescale of the
search channel; a local rescale inside a thresholding routine is a different thing.**

### What would change my mind

Nothing about saturation — that argument is structural. The only live question is whether
some *other* normalisation, one that is monotone and does not clip (a log or rank transform,
say), helps a downstream classifier. Monotone transforms cannot change a ranking, so this
would only matter for a learned model, not for the ranker.

---

## D4 — Report recall@K, not depth-to-target

**Date:** 2026-09-04

### The decision

The standard reporting unit becomes **recall@K**: fix a candidate budget K, and report the
fraction of mitotic figures found within the first K candidates.

Specifically: **the smallest K such that recall@K ≥ target on the worst ROI and the worst
click, reported per tumour domain.**

This replaces `read_50` (candidates read to find half the mitoses) and depth-to-100% as the
headline metrics. Both remain useful as diagnostics; neither drives a decision any more.

### Why, in plain words

The old metrics ask "how deep must the reader go to find X% of the mitoses?" The new one asks
"if the reader looks at K candidates, what fraction do they find?" These sound like the same
question read in two directions, and mathematically they are — but they behave completely
differently as decision tools.

**1. Depth is unbounded and hostage to single objects.** On `246.tiff`, reaching 98% of the
mitoses costs 1,966 candidates; reaching 99% costs **17,321**. One annotation, sitting deep in
the list, multiplies the workload by 8.8×. A metric that a single object can move by nine-fold
cannot be used to compare methods. recall@K is bounded in [0, 1] and moves smoothly.

**2. Depth is undefined whenever the target is unreachable.** On `246.tiff`, `read_100` is
`NaN` for every method and every click, because the proposal stage never finds 4 of the 116
mitoses — the ceiling is 0.965. You cannot compare methods on a metric that is missing for all
of them. recall@K is always defined.

**3. K is the decision the product actually makes.** The tool has to choose a list length
before the pathologist starts reading. "Depth to 100%" is not a setting anyone can dial.

**4. It makes both goals one number.** "Cover everything" and "keep the list short" stop being
two competing targets and become one curve, which is the only way to reason about the
trade-off honestly. The trade-off is severe: on `301.tiff`, the first half of the mitoses cost
about 1.5 candidate-reviews each, and the last four cost about **861 each**.

**5. `read_50` in particular is the depth where methods matter least.** The 2026-09-04 audit
found that in the one representation where clicking appeared to work, the entire `read_50`
result reproduces with **no click at all** — while the click's real contribution shows up only
at `read_80` (4–21% better) and `read_100` (32% better). Four months of work optimised the one
number that was least sensitive to the thing under test. That is the strongest single argument
for this change.

### Why "worst ROI, worst click, per domain" and not a median

Because the median hides the failure I actually care about. The requirement is that this works
on *any* domain and *any* region, and a method that averages well while collapsing on
lymphosarcoma is not acceptable. Two things this specifically guards against, both already
observed:

* **Reading burden varies about 2× depending on which cell the pathologist happens to click**,
  and in one earlier run the target was never reached at all on 2 of 5 clicks. A median over
  clicks hides that completely.
* **The ROIs cluster by tumour type**, so an average over ROIs over-weights whichever domain
  happens to have more images downloaded. Reporting per domain keeps that visible.

This is not a new invention — it is M6 in `Research Logs/2026-09-02-next-steps-plan.md`, and
the worst-of-5-clicks rule was already the pre-registered gate in
`Research Logs/2026-09-01-one-click-retrieval-literature.md` §7.

### What it costs

The number is already computed. `recall_at_budget` and `tp_at_budget` exist in every results
CSV in this project and have never been used in a write-up. Adopting this is a plotting and
reporting change, not a compute change.

**One data trap to remember:** the results CSVs emit one row per budget level, so each
(ROI, click, method) combination appears 8 times with identical `read_*` and `ceiling` values.
Any row-level count, p-value or correlation computed without deduplicating on
`(file_name, seed_index, arm)` silently multiplies the sample size by 8. Medians survive this;
statistics do not.

### The thing this decision does *not* settle

Choosing K is a clinical question, not a statistical one, and it should be put to a
pathologist rather than picked from a curve. The related open piece is **calibration**: given
a chosen K, being able to state "on an unseen slide, this list contains at least 90% of the
mitotic figures with 95% confidence." That is a solvable problem with standard machinery
(conformal risk control / risk-controlling prediction sets) and nothing in this repository
does it yet. Recorded here as the natural next decision, not as a decision.

### What would change my mind

Nothing about the metric. The open choices *within* it are the value of K and the target
recall, and both are product decisions waiting on clinical input.

---

## Cross-references

| decision | primary evidence |
|---|---|
| D1 `TM_CCOEFF` | `Research Logs/2026-09-02-next-steps-plan.md` U0, U1; `results/tm_ccoeff_headtohead.csv` |
| D2 no `tissue_mask` | `Research Logs/2026-09-02-next-steps-plan.md` M5(a); `Research Logs/2026-09-01-premise-test-audit.md`; `midog_utils/baselines.py:26-28` |
| D3 no rescale after deconvolution | `midog_utils/channels.py` (`to_hematoxylin` vs `to_hematoxylin_od`); `Research Logs/2026-09-01-premise-test-audit.md` Part 1 |
| D4 recall@K | `Research Logs/2026-09-02-next-steps-plan.md` M6; `Research Logs/2026-09-03-fp-reduction-framing.md` §3b, §7 |

## Still open, deliberately

These were raised alongside D1–D4 and are **not** decided. They are listed so that a year from
now it is clear they were considered and deferred, not overlooked. See
`Research Logs/2026-09-03-fp-reduction-framing.md` §§4, 5, 6b, 8.

1. **Is the false-positive mass redundant?** How many genuinely distinct visual types the
   ~7,400 ordinary-nucleus false positives fall into. If it is tens rather than thousands, the
   reading burden collapses through the interface rather than the model.
2. **Relevance feedback.** Using the pathologist's rejections as labelled negatives on this
   slide. Never tested on any path here.
3. **Is the deep tail made of contested annotations?** About a quarter of MIDOG++ mitotic
   labels are 2-of-3 majority votes; the deepest-ranked mitoses look enriched in them.
4. **Calibrated stopping.** Turning a chosen K into a coverage guarantee.
