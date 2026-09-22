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

**Moved.** D1 lives in its own file: [`DECISIONS_UNVERIFIED.md`](DECISIONS_UNVERIFIED.md#d1).

In one paragraph: `cv2.TM_CCOEFF` (unnormalised, contrast-sensitive) is used in place of
`cv2.TM_CCOEFF_NORMED`, on the strength of a `read_50` comparison (49/49 cells, 10/10 ROIs)
later re-derived at precision@K (28/28 domain x budget cells, 2026-09-09 amendment). Moved
out of this file not because it was reversed or found wanting, but because — unlike every
other entry here — it has not yet been rigorously, independently re-verified: it was a
simple initial experiment, re-run twice by the same method without an adversarial check of
its own premise.

The full entry carries the decision, the evidence, both amendments, and the open question
that keeps it out of this file.

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

## D10 — The self-hit filter is removed; the seed is excluded by masking one NMS radius before peak extraction

**Date:** 2026-09-17

**Status:** decided, not yet applied. The code change, its pre-registered numbers and its gates are
in [`SELF_HIT_MASKING_PLAN.md`](SELF_HIT_MASKING_PLAN.md). Record the application as a dated
amendment under this entry, not by editing it.

**Written for 7.5 µm, rewritten for 5 µm, then reverted to 7.5 µm, all on 2026-09-17 and
before it was applied.** The first version set the mask radius equal to the NMS radius
(7.5 µm). The user then set it to 5 µm, so that a neighbouring figure 5–7.5 µm from a masked
point could still be found, and this entry and the plan were rewritten for that. After
measuring how rarely such neighbours occur in MIDOG++ (below), the user reverted to 7.5 µm. The 5 µm variant stays fully measured: `../cleanup_harness/runs/selfmask_ref_post_5um`
and `../cleanup_harness/selfmask/superseded_5um/`.

### The decision

`find_and_suppress` no longer drops detections within `self_hit_radius = 5.0` px of the seed after
NMS. Instead, **after the deep-floor threshold and before `extract_peaks`**, every pixel within the
image's NMS radius of the template centre is marked unreachable in `valid`
(`template_match.mask_seed_disc`). The radius is D7's `evaluate.radius_px(mpp)`, 7.5 µm or
29.6–33.1 px, and the score map itself is never modified. There is no separate mask-radius
constant or `FSConfig` field.

Removed with the filter: `FSConfig.self_hit_radius`, `production.SELF_HIT_RADIUS`, and the
`max_peak_score`, `n_self_hits` and `seed_self_score` info keys. `n_seed_masked_px` is added.

### Why, in plain words

**The 5 px filter only ever removed the seed's own peak; NMS was what cleared its
neighbourhood.** The seed's own match is almost always the strongest peak around it. Greedy NMS
therefore deletes everything within one match radius of it before the filter runs, and the
filter then deletes the peak. On the 14 ROIs at seed 0:
- the seed's peak entered the 100-peak pool 14 of 14 times;
- it was the pool maximum 12 of 14 times.

So "nothing within one NMS radius of the seed" was already the behaviour, by accident of
ordering. The mask makes it explicit and drops a fixed-pixel constant that sat beside an
mpp-scaled radius.

**It covers the two cases the accident didn't.**
- *The seed's own peak falls below the `max_peaks` cutoff.* Nothing then suppresses its
  neighbourhood. On 301.tiff at seed 1, a secondary peak of the clicked cell 26.6 px away stayed
  in the list (rank 93 by `tm_score`, 65 by `chromatin_od`).
- *An augmented template bank moves the self-match off-centre.* With the cleanup harness's
  8-template bank, 013.tiff kept a self-match 12.66 px away at rank 7, and 233.tiff kept one
  5.39 px away at rank 0.

The mask left nothing within one NMS radius of the seed in every run measured:
- 28 production runs;
- 70 held-out testing-set runs;
- 112 runs at seeds 1–4;
- 84 recent-experiment runs;
- 28 augmented runs.

**Only `valid`, and the full radius, because the cheaper-looking versions add false positives.**
- *Blanking the score map in a ±7 px box* (the first draft of the plan) turns the slope of the
  seed's own match into a new peak one pixel outside the box. That peak landed 8.5 px from the
  seed at rank 39 on 403.tiff and at rank 11 on 460.tiff, for −1 TP at `chromatin_od`
  K = 10/20/30.
- *A ±7 px box on `valid` alone* still leaks secondary self-matches 9–30 px out, at ranks 0–10.

Leaving the score map intact lets `extract_peaks`'s own dilation keep the slope from becoming
peaks. The full radius covers the secondary matches.

**7.5 µm, not 5 µm: the neighbours a 5 µm mask would free almost never exist, and the rows it
lets back in are visible.** Both radii were run on the real code. They give bit-identical
detections on all 28 production-harness runs and all 84 recent-experiment condition runs.
- *What 5 µm would gain.* A second annotation 5–7.5 µm from the seed could be found from its
  own peak. In MIDOG++, 9 of 11,937 mitotic annotations have another annotation within 7.5 µm:
  four mitotic pairs 6.6–7.4 µm apart (245.tiff, 220.tiff, 255.tiff, 295.tiff) and one mitotic
  figure beside a look-alike (050.tiff). No two annotations are closer than 5.9 µm. Distances
  are between box centres. µm/px comes from the tumour type, because only 23 of 503 images
  are on disk, and the counts are the same at both ends of the human scanners' range
  (0.2263–0.2298). Three of the four mitotic pairs are within 0.3 µm of 7.5 µm, so a 7.5 µm
  mask may still leave those neighbours a peak just outside it.
- *What 5 µm would cost.* It keeps secondary matches of the clicked cell 5–7.5 µm from its
  centre, all unannotated:
  - 300.tiff seed 1 at `tm_score` rank 7, costing one TP at K = 50 under both rank keys;
  - 301.tiff seed 1 at rank 93 / 65 (the old filter's miss, which 5 µm doesn't fix);
  - 289.tiff (held-out) at `chromatin_od` rank 10;
  - with the 8-template bank, 245.tiff at ranks 4 and 8 and 548.tiff at ranks 16 and 23.

  In a click-to-verify tool, such a row is an extra box next to the cell just clicked.

**At `tm_score` it moves no precision number; at `chromatin_od` it moves a few by one.** Freeing
the seed's slot in D9's 100-peak cap admits one more candidate, almost always far from the seed.
- *Under `tm_score`* that candidate ranks last. No TP@K changed on any run measured: 14 ROIs × 5
  seeds, 35 held-out ROIs, and the bbox3way and hembbox runs.
- *Under `chromatin_od`* it can re-rank into the top K. Before the change, current code
  reproduces all 770 stored rows of the recent experiments exactly (bbox3way, hembbox,
  seed-robustness). After it, 9 of those rows change TP, 8 up and 1 down, each traced to that
  one extra row.

### What it costs

- **A real annotation within one match radius of the seed can't be found from its own peak.**
  This was already true whenever the seed's peak was in the pool.
  - It is rare: 9 of 11,937 mitotic annotations (0.08 %) have another annotation within 7.5 µm,
    the four mitotic pairs and one mitotic/look-alike pair above.
  - This retires the reason `midog_utils/FIND_AND_SUPPRESS_REFERENCE_DIFFS.md` gave for keeping
    the self-hit radius far below the match radius. NMS had been removing those neighbours all
    along.
- **Lists get one row longer.** `n_detections` rises by 1 on almost every run, so D9's post-NMS
  "89–99 candidates" becomes 90–100. Pool sizes from before and after the change are not
  comparable one for one. `tm_score` precision and recall are.
- **`chromatin_od` comparisons must not straddle the change.** A few cells move by ±1 TP in both
  directions, so re-run one side rather than mixing numbers from before and after.
- **One experiment breaks when re-run.**
  `production_hematoxylin_only/bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb` reads
  `prod.SELF_HIT_RADIUS` and `n_self_hits`. Its prompt, `BBOX3WAY_49ROI_3SEED_PROMPT.md`, pins
  "5 px self-hit removal".
- **Verification level, stated plainly.**
  - The design probes and every reference capture came from one session.
  - A separate review agent re-derived the production-harness numbers from the 7.5 µm captures
    and checked the code.
  - The held-out testing-set and 8-template-bank numbers, and the 5 µm comparison, come from
    running the real pre- and post-change code in that session. They have not been re-derived
    independently.
  - The close-pair counts come from the annotation database, with µm/px assigned by tumour
    type as described above.

### What would change my mind

- **Mitotic figures within 7.5 µm of each other turning out to matter.** For example, a product
  flow, dataset or multi-seed round where such pairs are common or clinically important;
  MIDOG++ has four. The 5 µm variant is measured and kept for that case. It would still need
  some other way to keep the seed's own secondary matches out of the top of the list.
- **A template bank whose self-match lands more than one NMS radius from the template centre**,
  for example larger rotations or scales than the 8-template bank tested. The disc would leak,
  and the mask would have to follow the fused map's own self-match instead of a fixed disc.
- **A multi-seed `chromatin_od` result where the extra candidate lowers precision
  systematically**, rather than moving single cells by ±1 in both directions.

### Amendment, 2026-09-17 — superseded

The user proposed a different mechanism the same day: blank the seed's own refined-template
footprint out of a copy of the search channel, before correlation, instead of clearing a
disc from `valid` after the deep-floor threshold. Recorded as **[D11](#d11)**, with its own
verification and gates; this entry's design was never applied to production. This body is
left otherwise unchanged, per this file's own rule that a decided entry is not rewritten.

---

## D11 — The self-hit filter is removed; the seed's own template footprint is blanked from the search channel before correlation, replacing D10's seed-disc mask

**Date:** 2026-09-17

**Status:** decided, not yet applied. The code change, pre-registered numbers and gates are in
[`SELF_HIT_MASKING_PLAN.md`](SELF_HIT_MASKING_PLAN.md). Record the application as a dated
amendment under this entry, not by editing it.

**Supersedes D10**, which is unchanged above except for its own amendment note. D10's disc-on-
`valid` design is fully built and gated (`../cleanup_harness/selfmask/`) but was never applied;
its assets remain untouched and still describe that design if it is ever revisited.

### The decision

`find_and_suppress` no longer drops detections within `self_hit_radius = 5.0` px of the seed
after NMS. Instead, after the template is cut from the search channel but before correlation
(`fused_response`) runs, a `base_size x base_size` square centred on the rounded template
centre — the same footprint the template itself was cut from — is blanked in a **copy** of the
search channel (`template_match.blank_seed_square`), filled with that channel's own whole-ROI
minimum computed before any blanking. The original, unblanked channel is never mutated, so
`chromatin_od` ranking (computed by the caller from the same channel reference) scores
candidates against real tissue, not the blanked patch. Unlike D10, `valid` is never touched —
the seed's own match simply has no source pixels left to correlate from.

Removed with the filter: `FSConfig.self_hit_radius`, `production.SELF_HIT_RADIUS`, and the
`max_peak_score`, `n_self_hits` and `seed_self_score` info keys — the same three D10 also
removes. `n_blanked_px` is added (D10's equivalent is `n_seed_masked_px`).

### Why, in plain words

D10's own reasoning for *why* the self-hit filter under-covers (the peak dominates its
neighbourhood via NMS, but only when it survives to the pool at all) applies unchanged here —
see D10 above. What differs is *how* the fix avoids the earlier rejected score-map-flooring
draft's failure (a false peak on the slope of the seed's own match, one pixel outside a masked
box): D10 keeps the score map intact and clears `valid` instead; this design blanks the image
itself, relying on a different argument — `TM_CCOEFF` (D1) against an exactly flat window
computes to (numerically) zero, not an extreme value, because the window's own local mean
equals the constant fill value, so the mean-subtracted term is exactly zero; and because
matching is a sliding window, the corrupted response tapers off with the window's shrinking
overlap rather than stepping sharply at a hard edge. Verified this session, both analytically
(a realistic-magnitude synthetic flat-window test: response ~3e-7, five orders below background
noise) and on real ROIs (403.tiff: the nearest pre-NMS peak to the blanked square sits 299 px
away; no boundary-artefact detection found anywhere measured).

**Why blank the image instead of reusing D10's disc.** This was the user's own proposal, tested
empirically rather than assumed to fail or succeed by analogy. It does not strictly dominate
D10 — see "What it costs" — but the blanked square is usually *smaller* than D10's full-NMS-
radius disc (`base_size` 23–51 px vs. `match_radius` 30–33 px), so it excludes less of the
image, at the cost of occasionally leaving a real secondary match reachable just outside the
square.

### What changes in the numbers (independently re-derived, not just re-read from the plan)

All of the following were recomputed this session by a fresh, from-scratch prototype
(`find_and_suppress_blanked`, modelled on `find_and_suppress.py`'s own body, not on the design
session's archived scripts) calling real, unmodified `midog_utils` primitives — not a
re-execution of the design session's own `blank_seed_variant_*.py` scripts, which were read
adversarially instead and used only as a cross-check.

- **49 ROIs, `chromatin_od`, pinned seeds:** pooled TP@K 270→269 / 460→461 / 589→590 at
  K=10/20/30 — exact agreement with the design session's own numbers, via an independent code
  path. The 3 (ROI, K) cells that move (128.tiff K10, 400.tiff K20/K30) and the 2 near-seed
  leaks (289.tiff 24.4 px rank 10, 548.tiff 28.0 px rank 56) reproduce exactly, including their
  precise scores and ranks. Zero bucket flips on shared rows across all 98 runs (49 ROIs × 2
  rank keys) — the TP movement is pure displacement, not ground-truth reassignment.
- **`tm_score`, same 49 ROIs (closing a gap the design session left open):** zero TP@K movement
  at every K. Production's actual default ranker (D5) is unaffected in every real,
  full-ROI, single-template run measured.
- **8-template augmented bank (`scales=(0.8,1.2)`, `n_angles=2`, `flips=(False,True)`),
  `harness.py`'s own crop methodology, plain `seed_index=0` draws (closing a second gap):**
  both of the two historically-documented self-hit-radius failures reproduce exactly —
  013.tiff (12.66 px, rank 7, score 20.196 — matching `probes/p3.log` to 5 decimal places) and
  233.tiff (5.39 px, rank 0, score 12.246) — and this design removes both. 245.tiff shows a new
  (smaller) near-seed leak at 24.2 px, the same "known cost" pattern as 289/548, not a
  self-hit-radius-style failure. No `tm_score` movement on any of 4 ROIs × 2 `scale_normalize`
  values with the correct (plain-draw) seeds.
- **Seeds 1–4, 14 canonical ROIs (closing a third gap):** 4 (ROI, seed, K) cells move —
  013.tiff seed 4 K20, 300.tiff seed 2 K50, 301.tiff seed 3 K30, 459.tiff seed 4 K50 — and all
  four match D10's own pre-registered seed-robustness numbers **exactly**, cell for cell,
  despite the two designs sharing no code. This is independent cross-validation that the shared
  "freed pool slot" mechanism behaves identically regardless of which masking approach frees
  the slot, in every case where both were measured.

### What it costs

- **Same class of cost as D10, narrower fix.** A real detection can survive within one match
  radius of the seed but outside the blanked square, whenever `base_size < match_radius`
  (true on the ROIs measured: base_size 23–51 px vs. match_radius ~30–33 px) — a structural
  property of this design, expected to recur, not a rare fluke. Measured: 289.tiff (24.4 px,
  548.tiff (28.0 px), 245.tiff under augmentation (24.2 px) — none reached top-10 or cost a TP
  in any run measured.
- **This design is strictly weaker than D10 on 301.tiff's seed 1.** D10's own §1 documents a
  secondary peak of the clicked cell 26.6 px from the seed surviving because it fell below the
  `max_peaks` cutoff, so nothing suppressed it. D10's full-match-radius disc (30-33 px) covers
  this; this design's `base_size`=23 px square (half=11) does not; re-verified this session at
  26.62 px, present unchanged in both the unpatched and patched code, in every run measured.
  It costs no TP@K here, but is the one documented case where D10's disc would do strictly
  better.
- **A non-production mechanism, found and fully explained, not a production risk.** Under
  `scale_normalize=True` (never used by production; `FSConfig`'s default and every production
  call use `False`) combined with `harness.py`'s 1024 px test crop (never used by production,
  which always searches the full ROI), blanking can, via `_robust_z`'s shared per-template
  normalisation statistics, perturb an *unrelated* peak's score by ~1e-3, occasionally flipping
  greedy NMS's survivor count among near-tied peaks in a "chain" suppression geometry (peak A
  suppresses B and C, but B does not suppress C) — moving `tm_score` TP@K on one non-canonical
  seed draw (K=20, 2→1). This pathway is structurally absent when `scale_normalize=False`
  (`_robust_z` is never called; raw `TM_CCOEFF` scores for any window not overlapping the
  blanked square are then provably bit-identical between arms), which is production's only
  configuration — consistent with zero `tm_score` movement measured across every real,
  full-ROI, `scale_normalize=False` run this session (49 ROIs + seeds 1-4 + the corrected
  augmented-bank rerun).
- **The `n_detections` unchanged-count cases are not what D10's text (written for its own
  design) would suggest.** Re-derived this session for the 3 of 49 ROIs where blanking doesn't
  change `n_detections`: on 401.tiff and 546.tiff the seed's own peak was never among the
  post-NMS survivors even in the unpatched run (`n_self_hits=0`), so freeing its slot changes
  nothing, by construction — not because `max_peaks_binding` was false (it was `True` on all
  three). 123.tiff is a different, coincidental case: the freed slot's replacement candidate
  happens to fall within NMS radius of another survivor, so post-NMS count drops by one
  (95→94) exactly offsetting the one self-hit the old filter used to remove.
- **Same downstream breakage as D10.** `production_hematoxylin_only/bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb`
  and its prompt read `prod.SELF_HIT_RADIUS`/`n_self_hits`; not ported, per D10's own
  disposition.

### Verification level, stated plainly

- The design's mechanism and pre-registered numbers (`SELF_HIT_MASKING_PLAN.md` sec 2-3) came
  from one prior session's real simulation.
- A separate verification session (this one) independently re-derived the mechanism claims
  (analytically and on real ROI data), the pooled §3a/§3b numbers (via a from-scratch
  prototype, not a re-execution of the prior session's scripts), and closed all three gaps the
  prior session left open (`tm_score` at scale, the augmented bank, seeds 1-4) — including
  finding and correcting one real bug in its own harness along the way (blanking before vs.
  after template extraction) and one real seed-selection bug in its own augmented-bank test.
- The real code patch, fresh bit-identical-by-design reference captures
  (`../cleanup_harness/runs/imageblank_ref_pre`, `imageblank_ref_post`), a design-specific
  `check`/`reference` script and negative tests were built in the same verification session,
  against a scratch copy — never applied to `midog_utils/`.

### What would change my mind

Same three conditions as D10 (see above): mitotic figures within 7.5 µm of each other turning
out to matter more broadly than MIDOG++'s 9 cases; a template bank whose self-match lands
further than this design's blanked square (which, unlike D10's disc, is *not* always
`nms_radius`-sized, so this risk is somewhat higher for augmentation configurations wider than
the 8-template bank tested); or a multi-seed `chromatin_od` result where the extra candidate
lowers precision systematically. Additionally: **evidence that 301.tiff seed 1's kind of gap
(a real neighbour outside a too-small blanked square) recurs often enough to matter** would
favour reverting to D10's full-radius disc instead.

### Amendment, 2026-09-17 — applied

Applied to `midog_utils/` via `SELF_HIT_MASKING_PLAN.md`. Harness and recent-experiment captures bit-identical to `imageblank_ref_post` (SHA-256 on all eight artefacts, the two `recent_*.pkl` included, not just the six `compare` diffs); `imageblank_check.py check`: PASS, 230 pairs, 9 TP@K changes, all attributed to the freed pool slot, plus 10 reported near-seed leaks accepted by design; `reference`: rows changed as pre-registered. `sanity`: PASS; `extref`: `n_detections` +1 uniformly on 84/84 rows with all 84 precision/recall rows within 0.0006, FAIL by design. Two documentation defects in the plan's own §5 text were found during application and both repaired: the call diagram still named the unblanked `hem` where the patched code correlates the blanked copy (§5.1 Edit 7 added the node but not its two downstream nodes), and §5.1 Edit 10's replacement text was affirmative inside a list headed "what `find_and_suppress` deliberately does NOT do". Edit 10 was first applied verbatim and the defect reported; on user authorisation it was repaired the same day and §5.1 Edit 10 corrected at source, so plan and walkthrough now agree. The Edit 7 repairs were not back-ported into §5.1 and live only in the execution log. Execution log: `../cleanup_harness/selfmask/image_blank_design/EXECUTION_LOG_2026-09-17.md`.

### Amendment, 2026-09-22 — near-seed leak checked against precision@10/20 on the tightened-template rerun, design kept

The 49-ROI × 3-seed × 3-arm (`default_51`/`gray_bbox`/`hem_bbox`) template-refinement experiment
was re-run post-D11 (commit `0af7ff1`, `rerun_bbox3way_postD11.py`), producing
`results/precision_at_k_49roi_3seed_chromatin_bbox3way_postD11_{per_run,tp_changes,top30}.csv`.
This is a different, later dataset from the one the 2026-09-17 amendment above verified (this one
uses the Otsu-tightened per-arm templates, `base_size` 23–51 px, not a single D8 template), so its
near-seed-leak count is measured separately here rather than folded into the "10" above.

**`n_near_click > 0` on 4 of 441 runs, all `non_human_findings`, none in `default_51`** (its 51 px
blank already exceeds `match_radius` ≈ 29.6–30.2 px): 289.tiff seed 0 `gray_bbox` rank 10
(d_click=24.6 px), 289.tiff seed 0 `hem_bbox` rank 7 (d_click=23.7 px), 548.tiff seed 0
`gray_bbox` rank 56 (d_click=30.5 px), 548.tiff seed 0 `hem_bbox` rank 55 (d_click=30.5 px). The
548.tiff pair sits far past any evaluated K.

**Checked directly against `..._tp_changes.csv` (the pre-D11-vs-post-D11 diff), not assumed:**
every K=10 and K=20 `tp_delta` row in all 441 runs has its entered/left detections thousands of
pixels from that run's click — the unrelated D9 freed-pool-slot mechanism, not this leak. Zero
K=10/K=20 precision changes are attributable to a near-seed leak anywhere in this rerun. The leak
changes a TP exactly once, at K=30: 289.tiff seed 0 `hem_bbox`, `tp_delta=-1`
(`entered` at d_click=23.7 px, matching the flagged leak exactly; precision@30 0.133→0.1).

**Decision: keep D11's blanking design as applied.** No revert to D10's full-`match_radius`
disc, no widening the blanked square. Basis: at the budgets this product cares about (K=10, 20)
the exposure this design accepts by construction (see "What it costs" above) is measured, not
theoretical, and it is invisible; it costs exactly one TP, at K=30, across 441 runs.

---

## Cross-references

| decision | primary evidence |
|---|---|
| D1 `TM_CCOEFF` | **[`DECISIONS_UNVERIFIED.md`](DECISIONS_UNVERIFIED.md#d1)** — moved 2026-09-12, pending independent verification; `Research Logs/2026-09-02-next-steps-plan.md` U0, U1; `results/tm_ccoeff_headtohead.csv`; `precision_at_k_budgets_14roi_normed/precision_at_k_budgets_14roi_normed.ipynb` Table C (the 2026-09-09 precision@K re-derivation) |
| D2 no `tissue_mask` | `Research Logs/2026-09-02-next-steps-plan.md` M5(a); `Research Logs/2026-09-01-premise-test-audit.md`; `midog_utils/baselines.py:26-28` |
| D3 no rescale after deconvolution | `midog_utils/channels.py` (`to_hematoxylin` vs `to_hematoxylin_od`); `Research Logs/2026-09-01-premise-test-audit.md` Part 1 |
| D4 recall@K | `Research Logs/2026-09-02-next-steps-plan.md` M6; `Research Logs/2026-09-03-fp-reduction-framing.md` §3b, §7 |
| D6 two prune criteria closed; re-match required | `f7_size_prune/RESULTS.md`; `f7_size_prune/PREREGISTRATION.md`; `Research Logs/2026-09-08-chromatin-vs-tm-and-residual-prune.md` |
| D5 `TM_CCOEFF` score is the ranker | `verify_chromatin_ranker.py` (re-derives all four); `results/f1_seed_sweep.csv`; `results/tp_fp_reading_depth.csv`; `tp_fp_separability_14roi.ipynb` and `results/tp_fp_14roi_budget_rule.csv` (the 2026-09-08 `od_contrast` amendment); `Research Logs/2026-09-08-tp-fp-separability-audit.md`; `results/morph_diag_bhattacharyya.csv`; `Research Logs/2026-08-31-chromatin-density-rerank.md` (the entry it corrects) |
| D7 NMS radius = 7.5 µm | `invariants.check_nms_radius`; `Research Logs/2026-09-04-f5-results.md` §7 and Correction 2; `tm_threshold_axis_sweep_largest_cc_high_z_r75.ipynb` and `results/tm_ccoeff_high_z_r75_vs_r50_seed_paired.csv` (the recall cost); `results/f5_nms_radius_ablation_spacing.csv` (annotation geometry); `find_and_suppress_high_threshold_precision.ipynb` Gate 4 (top-K invariance) |
| D8 seed template anchor | **[`D8_TEMPLATE_ANCHOR.md`](D8_TEMPLATE_ANCHOR.md)** — the entry itself, with the superseded 2026-09-09 decision and amendment kept verbatim; `pipeline_debug_visuals/template_anchor_halfpixel_fix.ipynb` (the half-pixel derivation and the anchor comparison, verified against the live `seed_selection` functions); `midog_utils/seed_selection.py` (`tighten_box_otsu` the gate, `tightened_base_size` and `tightened_template_box` the two superseded anchors) |
| D9 `max_peaks = 100` | `threshold_maxpeaks_ablation/max_peaks_100_variant.ipynb` (the measurement); `threshold_maxpeaks_ablation/z_floor_tightening_variant.ipynb` (`dynamic_z` branch — bit-exact equivalence proof) |
| D10 seed mask (one NMS radius) replaces the self-hit filter | **[`SELF_HIT_MASKING_PLAN.md`](SELF_HIT_MASKING_PLAN.md)** — the change, pre-registered numbers and gates; `../cleanup_harness/selfmask/` (the patch, `selfmask_check.py`, expected outputs, `probes/` including `analyze_realcode_output.txt`, `REFERENCE_PROVENANCE.md`); `../cleanup_harness/selfmask/superseded_5um/` (the measured 5 µm variant, with its patch, checker and `analyze_realcode_output_5um.txt`); `../cleanup_harness/runs/selfmask_ref_pre`, `selfmask_ref_post` and `selfmask_ref_post_5um` (captures); `../cleanup_harness/selfmask/SELF_HIT_MASKING_PLAN.midog_utils_full_draft.md` (the rejected ±7 px draft) |
| D11 image-blanking replaces D10's seed-disc mask | **[`SELF_HIT_MASKING_PLAN.md`](SELF_HIT_MASKING_PLAN.md)** (rewritten for this design) — the change, pre-registered and independently re-derived numbers, and gates; `../cleanup_harness/selfmask/image_blank_design/` (the patch, `imageblank_check.py`, expected outputs, the design session's own `blank_seed_variant_*.py`/logs/CSVs, this verification session's scripts); `../cleanup_harness/runs/imageblank_ref_pre`, `imageblank_ref_post` (captures); `../cleanup_harness/selfmask/SELF_HIT_MASKING_PLAN.valid_mask_disc_design.md` (D10's own design, archived unchanged) |

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

---

## Candidates raised for template refinement, not decided

Recorded for the record only. Not a decision, not evaluated against production's own
`gray_bbox`/`hem_bbox` arms or run through precision@K/NMS/matching, and not implemented in
`seed_selection.py`.

**2026-09-22 — a hybrid gray-gated, hem-anchored tightening prototype, with an erosion-first gray
gate.** `production_hematoxylin_only/hybrid_gray_hem_tightening_14roi.ipynb` prototypes a method
outside `seed_selection.py`: gate the hem (chromatin) component search to fall inside a
gray-channel component selected near the click. Its own closing notes flag two gaps: the
"completely inside" hem-containment rule is strict enough that eroding the gray gate instead of
dilating it collapses acceptance from 24/24 to 4/24 on a 24-click sample, and the fallback union
(no solidity/compactness gate) can span two separate nuclei in one box — observed on
`013.tiff`/243.
`production_hematoxylin_only/hybrid_erode50_vs_original_14roi.ipynb` tests a combined fix: erode
the gray channel's threshold *before* selecting a component (rather than selecting then dilating),
and loosen hem containment to `>=50%` (not 100%). Result on the same 24 clicks: 24/24 accepted,
23/24 land on the identical box as the original method, and the one that differs
(`245.tiff`/6243) is a visually confirmed fix of the fragmentation gap above. Whether this belongs
in `seed_selection.py`, and whether the hybrid approach beats `gray_bbox`/`hem_bbox` at all, is
unevaluated — see that notebook's own closing notes for the full numbers and caveats.
