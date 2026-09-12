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

## D5 — The production ranker is the `TM_CCOEFF` match score; chromatin density is a reported second axis, not the default

**Date:** 2026-09-04

### The decision

Detections are ordered by the search score. **Chromatin density (`chromatin.chromatin_density`,
`chromatin.rerank`) is not the production ranking key**, and no experiment may declare
`chromatin_od` its *primary* axis without measuring, on that run's own data, that it beats
`tm_score`. Both axes stay reported wherever both are cheap; the arbiter sits on the score.

This is a *confirmation* as much as a change. `chromatin.rerank` was never wired into the
pipeline — commit `7c3af93` said so itself ("neither altering existing pipeline defaults") and a
grep confirms its only callers are the three superseded 2026-08-31 probes `od_experiment.py`,
`od_workload_ab.py` and `od_seed_sweep.py`. What had drifted was the *documentation*: the module
docstring, the commit message, and then F1 → F4 → F5 in turn, each calling it "the shipped
ranker" on the authority of the commit hash rather than of a measurement.

### Why, in plain words

Four things, each re-derived from the raw data by `verify_chromatin_ranker.py` rather than read
off a published CSV.

**1. The reason the statistic was invented names a matcher we retired four days later.**
`chromatin.py`'s docstring justifies it as *"the signal `TM_CCOEFF_NORMED` is mathematically
blind to"*. D1 switched to `TM_CCOEFF`, which is **not** blind to it: scale a patch's contrast to
0.7 / 0.3 / 0.1 and `TM_CCOEFF` returns 0.700 / 0.300 / 0.101 of the full-contrast score, where
`TM_CCOEFF_NORMED` returns 0.999 / 0.992 / 0.930. D1's careful distinction still holds — the
matcher recovers *contrast*, not *absolute darkness*, so the two statistics are correlated, not
duplicates. But the stated justification is void, and no measurement replaced it.

**2. The rationale was established on the easy contrast.** The commit cites
`results/morph_diag_bhattacharyya.csv` making `mean_intensity` the strongest feature in all seven
domains. That file has **two** columns and the commit quoted one. Against *ordinary nuclei*:
0.93–3.70. Against *look-alikes* — the class that actually survives to the operating point —
0.192 (mast cell, n = 305), 0.161 (lymphosarcoma, n = 169), 0.332 (lung, n = 30). The larger
`vs_lookalike` values in that table sit on n = 7 and n = 8.

**3. At the configuration F4 and F5 declare primary, the two axes are indistinguishable.**
Recall@250 at z = 1.0 on `results/f1_seed_sweep.csv`: `chromatin_od` 0.6403 vs `tm_score` 0.6081,
reproduced exactly. But those 70 cells are 7 ROIs × 5 seeds × 2 arms. Clustered at the ROI — the
unit F5 §8 itself declares — the paired Δ is **+0.032, 95 % CI [−0.047, +0.112], p = 0.36,
positive on 3 of 7 ROIs**, and the aggregate rides on `201.tiff` (+0.200, n = 17 mitoses). Across
the full 48-cell (z × budget) grid the advantage is ≈0 at budgets 25–100 and **significantly
negative at 5,000** (p = 0.017–0.046, positive on 1 of 7).

**4. `od51` has the longest tail of any candidate ranker, which is the wrong shape for this
product.** Median `read_95` over the 7 ROIs: `od51` **4,488**, against `od31` 1,336,
`mask_od_mean` 1,748, `od_falloff` 1,759. It is beaten 6/7, 6/7 and 5/7. It wins at `read_50` and
loses everywhere deeper — and D4 already moved the reporting metric off `read_50` for exactly
this reason. A tool where a missed mitosis changes a grade cannot take that trade.

### Two findings that should not be lost, both currently in no committed file

* **51 px is reading the neighbours.** Between candidate pairs ≤ 25 px apart — whose 51 px
  windows share about half their pixels — `od51` correlates **0.72–0.84**, against `od31` at
  0.44–0.63 and the match score at 0.45–0.64. The window size comes from `tm.BASE_SIZE`, the
  annotation box, and was never swept. `od31` was first measured on 2026-09-04.
* **The rejected alternative wins, and the stated reason for rejecting it is void.**
  `chromatin.py` rejects the Otsu-component mean because it *"returns `None` on 9–12 % of
  detections"* and *"a gate can only discard"*. `tp_fp_feature_extract.shape_features` measures
  the largest component **gate-free**, and `shape_fail_rate = 0.0 on all 7 ROIs`. Gate-free,
  `mask_od_mean` beats `od51` on 2-class AUC (6 wins, 1 tie), on the look-alike contrast (5/7)
  and on `read_95` (6/7). The docstring's argument was against the *gated* version only.

### The standing constraint this decision carries

**Any run that ranks by a chromatin statistic — as primary axis, as a compared axis, or as a
re-measurement of this decision — must use the 14 ROIs in `images/extra_valid/`.** Every number
in this entry comes from `experiment.select_domain_images(images_dir='images')`, the
densest-per-domain draw of 7 that `select_domain_images` documents as optimistic for whatever is
being measured, at a single seed for everything except point 3. That is the same thin base
`chromatin_od` was promoted on, and it binds this entry as much as the one it corrects.
`images/extra_valid/MANIFEST.md` defines the set: 2 ROIs per tumour type, `n_mitotic >= 15`,
border-filtered unanimous seed pool >= 5.

### What it costs

**Nothing at runtime** — the pipeline never used it, so no result changes and no code path moves.
`chromatin.py` stays as it is, exercised by every experiment that reports both axes.

The real cost is that a genuine and possibly better signal is being held at "reported, not
default" while the evidence for it is thin in both directions. `od31`, `od_falloff`,
`mask_od_mean` and `od_contrast` (amendment below) all have prima facie evidence now and none of
them has a seed sweep. Deciding on today's data would repeat the mistake this entry exists to
correct.

### What would change my mind

`od31`, `od_falloff`, `mask_od_mean` and `od_contrast` added as axes in `f1_seed_sweep.py`'s
`AXES` / `AXIS_RANK_KEY` (two dicts, existing harness), swept over 5 seeds on the 14 ROIs of
`images/extra_valid`, and read at recall@K per D4 with the paired Δ clustered at the ROI. A
chromatin-family axis that beats `tm_score` there — with a CI excluding 0 and a majority of ROIs
positive — earns the default. Nothing less should move it, in either direction.

### Amendment, 2026-09-08 — `od_contrast` is a candidate axis, on the same terms as the other three

`tp_fp_separability_14roi.ipynb` measures a fourth chromatin-family axis this entry did not name:
**`od_contrast` = `od51 − od_ctx`** — the candidate's chromatin density minus the darkest half of a
121 px neighbourhood. Unlike `od51` it is a **within-image difference**, so the scanner's
uncalibrated optical-density scale appears in both terms and cancels. That is not a small
distinction here: on the 14 ROIs, `od51` is the *worst* of six axes at median `read_95` (3,709.5)
and `od_contrast` is the one that satisfies the shortlist rule most cheaply.

Under the operating rule set for the shortlist work — 95 % mean recall across ROIs with a **90 %
floor on every individual ROI** — `od_contrast` needs **K = 5,655** against `tm_score`'s **12,739**:
a 2.3x shorter list at the same guarantee, from a statistic that needs no training and no
calibration, only the ROI's own candidates. The greedy match was re-run in each ranker's own order,
per D6(b).

**Why this is a note and not a decision.**

* **One seed.** This entry's standing constraint is met on ROIs (the 14 of `images/extra_valid`)
  and *not* on seeds. `245.tiff` sets the 90 % floor for **all six** axes tested, so the ordering is
  a ranking of one ROI's depth curve at one click; that ROI's own click-to-click SD of `read_95` is
  **2,153** candidates (`results/tm_ccoeff_headtohead_seed_variance.csv`). The margin over
  `tm_score` is ~7,000 on that ROI and survives the band. The margin over a supervised
  leave-one-domain-out model is 423 — 0.2 SD — and does not; that comparison is not resolved.
* **Best-of-six on its own data.** `od_contrast` was selected from six axes scored on the same 14
  ROIs it is then reported on, with the selection turning on a single ROI.
* **It wins under a floor and loses on the median** (3,074 vs `od_falloff`'s 1,955). The two
  statistics disagree, and it is D4's "worst ROI, worst click" framing that picks `od_contrast`.
  Any sweep must read it at recall@K per D4, not at median depth.

So: add it to the sweep named above, do not promote it. And note the ceiling this shares with every
other axis in the family — against pathologist-marked look-alikes `od_contrast` reads AUC 0.58–0.89,
the same collapse `od51` and `od_falloff` show. It shortens the list; it does not separate mimickers.

---

## D6 — Two intrinsic candidate-prune criteria are closed; and no prune may be scored without re-matching

Date: 2026-09-08. Evidence: `f7_size_prune/RESULTS.md` (F7, 14 ROIs x 5 seeds, reproduction gate
70/70 exact); `Research Logs/2026-09-08-chromatin-vs-tm-and-residual-prune.md` (the residual probe,
7 ROIs, seed 0).

### The decision

**(a) Two specific criteria are closed.** Deleting candidates from the ranked list, post-NMS, by

* **hematoxylin purity** — the share of a candidate's unmixed optical density that is hematoxylin, and
* **physical size in µm** — the equivalent diameter of the nucleus its pixel falls in,

does not reduce reading burden and is not to be revisited without new evidence of the kind named
below. This closes two criteria on the `TM_CCOEFF` pool. It does **not** close intrinsic pruning as
a class: texture, solidity, and anything derived from the reader's own rejections are untested.

**(b) Any prune experiment must re-match after the deletion.** Reusing detection ranks recorded
before a prune is not an approximation, it is a different quantity. This is a standing requirement
on the method, independent of which criterion is being tested.

### Why, in plain words

**(a)** Both criteria fail for the same structural reason: **the search already encodes the
property, so the prune re-derives a filter that lives upstream of it.**

`cv2.matchTemplate` runs on `hematoxylin_od`, which is `rgb2hed(...)[:, :, 0]` — eosin and the
unmixing residual are projected out of the search image before the first correlation. Melanin, red
cells and ink are therefore *dim* in the channel being searched and cannot produce a high
correlation with a nucleus-shaped hematoxylin template. Measured: median hematoxylin fraction among
false positives **falls** from the best-ranked quartile to the deepest on 4 of 7 ROIs, is flat on 2
and rises on 1 — the head of the list is the *purest* chromatin. At zero true-positive loss the
criterion removes a median 2.4 % of the working region, and on 459.tiff its sign **inverts**.

Size fails on transfer rather than on mechanism. Fitted to a slide's own mitoses it works — shorter
on **14/14 ROIs** at 90 %, 95 % and 100 % recall (median −76 / −142 / −341, p = 0.0001). Held out by
domain it is inert: six of seven folds inherit the bound **[0.50, 99.85] µm**, a 200-fold interval
that excludes nothing, wrecked by two objects out of 5,175 measurable mitoses (a 0.50 µm
segmentation artefact and a 99.85 µm merged clump). The pre-committed a-priori bound of [4, 18] µm —
chosen from nuclear biology before any data was seen — **loses mitoses on 8 of 14 ROIs** and pushes
depth-to-95 % out by a median of 1,115 candidates. Mitotic size is learnable within a slide and does
not survive crossing a domain, which is exactly the property a "universal criterion" needed.

**(b)** `evaluate.greedy_match` gives each annotation to the **best-ranked** detection within the
match radius. Delete that detection and the annotation is re-claimed by a much deeper one — recall
is unchanged and the read gets far longer. On 233.tiff under the a-priori bound, recall stays at
1.000 while **depth-to-95 % goes 1,328 → 17,341**, a 13x increase, from deleting 8 % of the pool. A
pipeline that reuses stored ranks reports that as a free win.

### What it costs

Nothing at runtime; neither criterion was ever in the pipeline. The cost is that the "prune, then
re-rank by `TM_CCOEFF`" idea is now closed on both axes it was proposed on, and the two levers
`Research Logs/2026-09-03-fp-reduction-framing.md` §1 identified as genuinely different from ranking
— reader-derived negatives, and set-level selection — remain untested.

One instrument caveat, and its resolution. `on_nucleus` runs **0.22-0.71** per ROI, so only
**78.9 % (5,175 / 6,555)** of true positives have a measurable size and the rest were kept by rule —
which raised the question of whether the size null was really a statement about the candidate set.
It is not: `on_nucleus` is **0.990 in the top 100 candidates and 0.944 in the top 1,000**, against
0.463 over the whole pool, with the head-above-tail gap positive on **14 of 14** ROIs (median
+0.696, p = 0.0001). The list is nucleus-centred wherever it is actually read; the low pooled figure
is the unread tail (`f7_size_prune/RESULTS.md` §6b).

**That result generalises D6(a)'s mechanism, and is the reason to expect further criteria of this
shape to fail.** Stain purity, physical size and "is this even a nucleus" are each prunable exactly
where pruning is worthless — deep in a list nobody reaches. A criterion only pays if it separates
candidates *within the first few thousand*, and the first few thousand are already dark, round,
hematoxylin-positive, correctly-sized nuclei by construction of the search.

### What would change my mind

For **(a)**, on either axis: a criterion measured on a candidate set that the search channel does
*not* already encode, or a bound set by a method that is not a quantile of pooled true positives —
held out by whole domain per `fp_filter_domain.py`, scored on recall@K per D4, with the deletion
rate reported **restricted to candidates above the working depth** and the per-domain sign reported
so an inversion cannot hide in an average. Generalising D6(a) to intrinsic pruning as a *class*
would need the same treatment of at least a texture and a solidity axis; it has not been done.

Nothing changes **(b)**. It is arithmetic, not evidence.

---

## D7 — The NMS radius is the evaluation match radius, 7.5 µm

**Date:** 2026-09-09

### The decision

Distance NMS suppresses at **`evaluate.radius_px(mpp)` — 7.5 µm, converted per image**, which is
29.6–33.1 px across these scanners. This is what `FSConfig(nms_radius=None)` has always
documented; the decision is to stop overriding it. `NMS_RADIUS_UM = 5.0` in any new script is a
bug, not a configuration choice.

### Why, in plain words

**It is the only value that satisfies the repo's own invariant.** `invariants.check_nms_radius`
raises unless the suppression radius *equals* `radius_px(mpp)`, with the reason in its own error
string: *"suppression and scoring are using different neighbourhoods"*. If NMS keeps two
survivors closer together than the match radius, both can sit inside one ground-truth object's
radius — one is credited a true positive and the other a false positive **on the same object**,
which puts duplicates on the FROC's false-positive axis. `FSConfig.nms_radius`'s comment records
exactly that happening at a fixed 25.0 px.

**5.0 µm was never a considered choice — it propagated.** `tm_threshold_axis_sweep_v2.ipynb`
changed border padding *and* the radius in one step and said in its own cost table that the two
could not be separated. `..._largest_cc.ipynb` inherited it, `high_z` inherited it deliberately
(so its shared columns would reproduce byte-for-byte), and from there it was copied into
`recall_workload_ledger.py` and `tp_fp_feature_extract.py`.

**F5 separated the two factors and the radius was not the part that helped.**
`f5_nms_radius_ablation.py` held padding ON in every arm and moved only the radius: no benefit on
either ranking axis, at a cost of 23–37 % more candidates. Its §7 recommends the 7.5 µm default,
and Correction 2 leaves that unchanged while weakening the reasoning to *"no measured benefit at
real cost"*.

**And the head of the list barely notices.** NMS keeps the highest-scoring peaks at any radius, so
a score-ranked top-K is close to radius-invariant: `topk_churn` 0.0210 at K = 250 on `tm_score`
(1.1 % of membership), and at K = 20 over 14 ROIs the two radii give **the same 127 true positives
with 279 of 280 coordinates shared** (`find_and_suppress_high_threshold_precision.ipynb`, Gate 4).

### What it costs — and this is not free

Two costs, both measured on the 7-ROI × 5-seed paired re-run
(`tm_threshold_axis_sweep_largest_cc_high_z_r75.ipynb`), and neither priced by F5, which measured
`recall@K` rather than what the proposal stage finds at all:

* **The deep pool loses mitoses.** `deep_recall` is 1.0 in 35/35 cells at 5.0 µm and in **31/35**
  at 7.5 µm — **6 mitoses across 4 cells** (301 s4 and 246 s2 lose 2 each, 459 s4 and 548 s3 one
  each). This is *not* ground-truth merging: per `results/f5_nms_radius_ablation_spacing.csv`,
  only one mitotic pair in the whole 14-ROI set (245.tiff, 26.57 px) is close enough for the
  7.5 µm radius to merge, and none of the four affected ROIs has one — F5 §2a puts the ceiling of
  the ground-truth-merging mechanism at **one annotation** for the entire set, so it cannot account
  for six. **What does account for it has not been traced**, and the entry does not claim it has.
  The available hypothesis is peak shadowing: a survivor now occupies a disc the size of the match
  radius, so a strong detection can suppress the only peak that would have claimed a real mitosis
  just outside that mitosis's own match radius. F5 traced the adjacent — not identical — effect on
  301.tiff, where 5 mitoses had their on-centre peak suppressed by an off-centre peak 23–29 px away
  on the same nucleus; that geometry changes *which* detection claims the mitosis (and so moves
  `z_max`) but leaves it claimed, so it explains the second cost below and not the first. Tracing
  the recall loss on one affected cell is the cheapest open piece of work this entry names.
* **The deployable cutoff falls.** The global `z_max` — the highest cutoff keeping every mitosis
  the search finds — drops from **1.308 to 0.438**, because a mitosis that keeps its claim is
  sometimes claimed by a much weaker detection once the stronger one is suppressed
  (`results/tm_ccoeff_high_z_r75_vs_r50_tp_substitution.csv`).

The gain on the other side is volume: the 5.0 µm pool is **58 % larger** at the deep floor and
**44 % larger** at `z = 1.0` (means over the 35 cells) for those extra mitoses.

**Pool sizes are therefore not comparable across the change**, and any artefact written at 5.0 µm
must say so. The decision rests on the invariant and on duplicate-free scoring, *not* on the
radius being costless — it is not.

### What is already migrated, and what is not

`f6_padding_ablation.py`, `tm_threshold_axis_sweep_largest_cc_high_z_r75.ipynb` and
`find_and_suppress_high_threshold_precision.ipynb` use 7.5 µm. Still carrying the override:
`recall_workload_ledger.py:93`, `tp_fp_feature_extract.py:69`, and
`tm_precision_under_50_candidates.ipynb`. They are **not** wrong as published — each is internally
consistent and gate-locked to its predecessor — but anything built on them from now on should
re-run at 7.5 µm rather than inherit.

### What would change my mind

A domain whose annotations really are packed closer than 7.5 µm apart, where the merge the
invariant tolerates becomes common rather than one pair in fourteen ROIs. Or a demonstration that
the 6-mitosis proposal-stage loss above reaches the **top of the list** — Gate 4 says it does not
at K = 20, but that is one click per ROI on 14 ROIs, and the 4 affected cells are all at seed
indices this notebook did not run. That is the measurement this entry is most exposed to.

---

## D8 — The seed template is cut from the accepted component, not from the click

**Moved.** D8 lives in its own file: [`D8_TEMPLATE_ANCHOR.md`](D8_TEMPLATE_ANCHOR.md).

In one paragraph: the click **gates** the seed and nothing else — its own pixel must land inside an
Otsu component clearing `min_area`/`max_area_frac`/`min_solidity`, or the seed is refused and
redrawn. The template's size *and* centre then both come from that accepted component, with the
centre computed as `(x0 + x1 - 1) / 2` — the centre of the pixels the bounding box covers, not
`(x0 + x1) / 2`, the centre of the area it spans. The two differ by half a pixel; the consumer
(`read_padded_patch`) rounds its input to a pixel index, so only the former is correct, and only
the former guarantees the template contains the whole component.

The full entry carries the evidence, the costs, the open outcome question, and — kept verbatim per
this file's own history policy — the superseded 2026-09-09 decision and its amendment.

---

## D9 — `max_peaks = 100` before NMS

**Date:** 2026-09-12

### The decision

Cap `extract_peaks`'s `max_peaks` argument to **100** (pre-NMS, by raw score) in the production
call. Not `FSConfig.max_peaks` (still 250000, unused — `find_and_suppress()` is not the caller
in the current precision pipeline); it is the `MAX_PEAKS` constant in
`production_seed_precision_at_k_chromatin_half_pix_fix.ipynb` cell 1, currently `2_000_000`
(never binds). `DEEP_FLOOR_Z` is unchanged at −1.5.

### Why

Measured on the 14 ROIs of `images/extra_valid`, one seed each, `tm_score` arm: precision@10/20/30
is unchanged on 14/14 ROIs at every budget (pooled 0.5000 / 0.4536 / 0.4286, identical to the
uncapped pool), recall never lower. NMS runtime drops ~99.7 % (~100 candidates into NMS instead of
~17,300), for a 12.6–15.7 % reduction in full click-to-list latency across three independent
timing runs. Extraction cost is unaffected — it is set by the response-map size, not by
`max_peaks`. A per-ROI dynamically-solved z targeting the same 100-candidate count reproduces the
capped pool bit-for-bit, so 100 is not an arbitrary count — it is the tightest fixed-count cap the
data supports at zero measured cost through K=30.

### What it costs

The capped pool under-delivers a K=100 budget: post-NMS, every ROI lands at 89–99 candidates,
never the full 100. Not tested past K=30. Single seed per ROI — D4/D5's multi-seed bar for a
production default is not met.

### What would change my mind

Any product requirement for K > 30 — the cap would need headroom above 100 to avoid
under-delivery. A multi-seed re-run showing precision or recall moves at K≤30.

---

## Cross-references

| decision | primary evidence |
|---|---|
| D1 `TM_CCOEFF` | `Research Logs/2026-09-02-next-steps-plan.md` U0, U1; `results/tm_ccoeff_headtohead.csv`; `precision_at_k_budgets_14roi_normed/precision_at_k_budgets_14roi_normed.ipynb` Table C (the 2026-09-09 precision@K re-derivation) |
| D2 no `tissue_mask` | `Research Logs/2026-09-02-next-steps-plan.md` M5(a); `Research Logs/2026-09-01-premise-test-audit.md`; `midog_utils/baselines.py:26-28` |
| D3 no rescale after deconvolution | `midog_utils/channels.py` (`to_hematoxylin` vs `to_hematoxylin_od`); `Research Logs/2026-09-01-premise-test-audit.md` Part 1 |
| D4 recall@K | `Research Logs/2026-09-02-next-steps-plan.md` M6; `Research Logs/2026-09-03-fp-reduction-framing.md` §3b, §7 |
| D6 two prune criteria closed; re-match required | `f7_size_prune/RESULTS.md`; `f7_size_prune/PREREGISTRATION.md`; `Research Logs/2026-09-08-chromatin-vs-tm-and-residual-prune.md` |
| D5 `TM_CCOEFF` score is the ranker | `verify_chromatin_ranker.py` (re-derives all four); `results/f1_seed_sweep.csv`; `results/tp_fp_reading_depth.csv`; `tp_fp_separability_14roi.ipynb` and `results/tp_fp_14roi_budget_rule.csv` (the 2026-09-08 `od_contrast` amendment); `Research Logs/2026-09-08-tp-fp-separability-audit.md`; `results/morph_diag_bhattacharyya.csv`; `Research Logs/2026-08-31-chromatin-density-rerank.md` (the entry it corrects) |
| D7 NMS radius = 7.5 µm | `invariants.check_nms_radius`; `Research Logs/2026-09-04-f5-results.md` §7 and Correction 2; `tm_threshold_axis_sweep_largest_cc_high_z_r75.ipynb` and `results/tm_ccoeff_high_z_r75_vs_r50_seed_paired.csv` (the recall cost); `results/f5_nms_radius_ablation_spacing.csv` (annotation geometry); `find_and_suppress_high_threshold_precision.ipynb` Gate 4 (top-K invariance) |
| D8 seed template anchor | **[`D8_TEMPLATE_ANCHOR.md`](D8_TEMPLATE_ANCHOR.md)** — the entry itself, with the superseded 2026-09-09 decision and amendment kept verbatim; `pipeline_debug_visuals/template_anchor_halfpixel_fix.ipynb` (the half-pixel derivation and the anchor comparison, verified against the live `seed_selection` functions); `midog_utils/seed_selection.py` (`tighten_box_otsu` the gate, `tightened_base_size` and `tightened_template_box` the two superseded anchors) |
| D9 `max_peaks = 100` | `threshold_maxpeaks_ablation/max_peaks_100_variant.ipynb` (the measurement); `threshold_maxpeaks_ablation/z_floor_tightening_variant.ipynb` (`dynamic_z` branch — bit-exact equivalence proof) |

## Still open, deliberately

These were raised alongside D1–D4 and are **not** decided. They are listed so that a year from
now it is clear they were considered and deferred, not overlooked. See
`Research Logs/2026-09-03-fp-reduction-framing.md` §§4, 5, 6b, 8.

1. **Is the false-positive mass redundant?** How many genuinely distinct visual types the
   ~7,400 ordinary-nucleus false positives fall into. If it is tens rather than thousands, the
   reading burden collapses through the interface rather than the model.
2. **Relevance feedback.** Using the pathologist's rejections as labelled negatives on this
   slide. Never tested on any path here.
3. ~~**Is the deep tail made of contested annotations?**~~ **Closed 2026-09-08 — it is not.**
   Tested on `results/tm_recall_workload_tp_ledger.csv` (14 ROIs x 5 seeds, every captured
   mitosis's rank in its own deep pool), which the 2026-09-03 framing could not do at n = 32.
   Per ROI, the deepest 10 % of true positives are *less* unanimous on only **7 of 14** ROIs
   (sign test p = 1.0), and unanimity is flat across rank deciles (0.73-0.80, no trend).
   **Two caveats belong with this null.** The 6,555 rows are not independent — the same
   annotation recurs across up to 5 seeds — so the ROI-level sign test is the valid reading and
   the pooled decile counts are inflated. And per `dataset.py`, for category 1 the only label
   multisets present are `(1,1)` and `(1,1,2)`, so `unanimous` is exactly `n_votes == 2`: a
   proxy for "did this need a third reader", not a measure of label quality.
4. **Calibrated stopping.** Turning a chosen K into a coverage guarantee.
