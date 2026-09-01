# Ranking by chromatin density instead of correlation

Date: 2026-08-31
Code: `midog_utils/chromatin.py`, `od_experiment.py`, `od_workload_ab.py`, `od_seed_sweep.py`,
`od_rerank_probe.py`
Results: `results/od_experiment_{arms,workload,detections}.csv`, `results/od_workload_ab.csv`,
`results/od_seed_sweep.csv`

Every number below was independently re-derived by an adversarial audit (see "What the
audit changed"). Where the audit and the first draft disagreed, the audit won.

## Why every previous experiment came back neutral

Score threshold 0.25 -> 0.5 -> 0.75, multi-scale fusion, augmentation count, search
channel, bbox tightening, NMS ordering. Every one was measured; every one was rejected or
found neutral. That reads as a run of bad luck. It is not.

**Each is a monotone operating-point move on one fixed ranking function.** None could have
added information, so none did. The ranking function itself was never examined.

`cv2.matchTemplate(..., TM_CCOEFF_NORMED)` mean-centres and L2-normalises *both* template
and window, so it is invariant to `I -> aI + b`. Confirmed numerically from both sides: a
30%-contrast copy of a patch scores 0.999999 against the original's 1.0, and darkening only
the object on a uniform background leaves the score at 1.000000 at every darkness tested.

Meanwhile `results/morph_diag_bhattacharyya.csv` -- measured back in August, independent of
any template matching -- makes `mean_intensity` the **strongest** mitotic-vs-ordinary-nucleus
feature in all seven domains (0.925-3.703). The runner-up varies by domain
(`tightened_size` in 4, `solidity` in 2, `area` in 1) and never exceeds 1.051.

**The scoring function is invariant to the strongest measured discriminator in the dataset.**

*Calibrated caveat*: "blind" overstates it. The invariance is exact only inside the 0-255
range -- clipping breaks it, zero-variance white windows score 0 (a crude brightness gate),
and empirically the correlation score still reaches AUC 0.61-0.81 on TP-vs-unannotated.
It is a much weaker ranker than chromatin density, not a blind one.

## Fix 1: a score floor in units of each map's own noise

`FSConfig.score_threshold = 0.5` is a constant in `TM_CCOEFF_NORMED` units, but that
statistic's null width depends on template size, which varies 25-51 px with the seed.
Measured as robust z above each fused map's own median (`template_match._robust_z`), seed 0:

| image | template | MAD scale | **z of the 0.5 floor** |
|---|---:|---:|---:|
| 301.tiff | 35 px | 0.303 | **1.72** |
| 405.tiff | 51 px | 0.243 | 2.15 |
| 002.tiff | 35 px | 0.188 | 2.72 |
| 350.tiff | 41 px | 0.169 | 3.04 |
| 506.tiff | 43 px | 0.159 | 3.20 |
| 201.tiff | 51 px | 0.139 | 3.59 |
| 246.tiff | 51 px | 0.110 | **4.54** |

One constant searches 301.tiff at 1.7 sigma and 246.tiff at 4.5 sigma. That is the whole
explanation for 246.tiff returning **3-22 detections against 115 mitotic figures** on every
seed tried, while 301.tiff returned up to 15,484.

**The reparameterisation is what is defensible; the level still has to be chosen.** Swept
across the 7 ROIs (`results/od_experiment_arms.csv`):

| floor | mean full-list recall | list-length spread | floor-limited arms |
|---|---:|---:|---:|
| fixed `score >= 0.5` | 0.544 | 2 671x | 0 |
| `z >= 3` | **0.387** | **505x** | 0 |
| `z >= 2.5` | **0.763** | **4.8x** | 0 |
| `z >= 2` | 0.935 | 1.8x | 1 |
| `z >= 1.5` | 0.955 | 2.1x | 4 |
| `z >= 1` | 0.955 | 2.3x | 6 |

**At z = 3 the per-image floor is worse than the constant it replaces** (0.387 vs 0.544) and
the spread is 505x, not collapsed -- 301.tiff returns 6 detections. Robust-z equalises each
map's *null*, not its *tail*, so equal z is not equal list length. Both halves of "collapses
the spread and raises recall" are properties of **z = 2.5 specifically**, which is the lowest
swept level that stays above the 0.25 extraction floor on every ROI (246.tiff needs
z >= 2.27). The z >= 2/1.5/1 rows contain 1/4/6 floor-limited arms and are not clean
measurements. z = 2.5 is clean, and untuned beyond that constraint.

Also worth stating plainly: 301.tiff and 405.tiff *shorten* (10x, 2x) while 246.tiff
*lengthens* 5 -> 4,951 (990x) and 201.tiff 24.7x. Total candidates across the seven rise only
25,592 -> 28,694 (+12%), so the +0.219 mean recall is redistribution, not volume -- but the
per-image moves are large in both directions.

`coverage_frac`, which `evaluate.py`'s docstring says to read before any full-list number:
0.0005-0.836 in the fixed arm (301.tiff's 0.949 recall sits at 0.836), 0.107-0.450 at z=2.5.

## The statistic

`chromatin.chromatin_density`: mean of the darkest 10% of pixels in the 51 px window around
a detection, read off `chromatin.hematoxylin_od` -- the **unclipped** colour-deconvolved
hematoxylin channel.

Unclipped matters. `channels.to_hematoxylin` min-max rescales between the ROI's 0.5/99.5
percentiles, so exactly 0.50% of pixels sit at the 255 ceiling, and a window over dense
chromatin returns 255.000 exactly. The audit found 66-77 detections tied at the ceiling on
002.tiff -- its whole top-100 saturated -- making the ranking there arbitrary rather than
informative. Switching to the unclipped channel cuts the largest tie block from 66 to 10;
`chromatin.rerank` additionally uses a stable sort so residual ties keep the incoming score
order instead of quicksort's arbitrary one. **This was a real defect and it was suppressing
the result**: 002.tiff's read-to-50% went 50 -> 8 candidates after the fix.

The obvious alternative -- segment with `seed_selection.tighten_box_otsu` and average the
component -- is worse once every detection must get a value, because the segmenter returns
`None` on 9-12% of detections and those are not a random subset. It also costs recall, since
a gate can only discard. The window statistic needs no gate and is defined everywhere.

## Fix 2: rank by chromatin density, not by correlation

Same candidate set (z >= 2.5), only the sort key changes. Candidates a reader must work
through to reach 50% of the mitotic figures (`results/od_workload_ab.csv`, seed 0):

| image | n mitotic | candidates | score-ranked | **chromatin-ranked** | reduction |
|---|---:|---:|---:|---:|---:|
| 002.tiff | 8 | 5 206 | 315 | **8** | 39x |
| 201.tiff | 17 | 3 133 | 994 | **28** | 36x |
| 506.tiff | 11 | 6 101 | 3 122 | **138** | 23x |
| 405.tiff | 13 | 2 126 | 1 297 | **67** | 19x |
| 246.tiff | 115 | 4 951 | 1 794 | **95** | 19x |
| 301.tiff | 217 | 1 258 | 880 | **247** | 3.6x |
| 350.tiff | 3 | 5 919 | 2 531 | **2 362** | 1.1x |

Median **19.4x**, winning 7/7 -- though 350.tiff (n=3, target = 2 mitoses) is a tie dressed
as a win, not a seventh victory.

### It survives seed variance -- the main open risk, now closed

35 runs, 5 seeds per ROI, 30 distinct seed annotations (`results/od_seed_sweep.csv`).
Seed variance is as large as the 2026-08-27 log warned: on 301.tiff the fixed-0.5 list ranges
**1,985-15,484 detections** and the template 25-51 px, purely from which annotation is
clicked; on 201.tiff, 127-4,219.

| | seed 0 only | 5 seeds (35 runs) |
|---|---:|---:|
| median workload reduction at 50% | 19.4x | **19.4x** (IQR 10.7-39.6) |
| runs where chromatin ranking wins | 7/7 | **29/29** |
| excluding the 2 floor-limited runs | -- | 22.6x, **27/27** |
| median full-list recall, fixed -> z2.5 | -- | 0.455 -> **0.727** |
| list-length spread, fixed -> z2.5 | -- | **5 161x -> 6.4x** |

`recall@K` also improves once the saturation defect is fixed: chromatin ranking is >= score
in **35/35** runs and strictly better in 29/35. The 002/506 `recall@K` regressions reported
in the first draft were the tie lottery, not a property of the ranking.

Six of 35 runs reach 50% under neither ranking; those are excluded from the ratio, not
counted as wins.

## Where this leaves look-alikes -- the first draft got this backwards

The full-list false-positive composition is **170:1 unannotated to look-alike** (28,284 vs
166 pooled at z=2.5; 25,220 vs 135 in the earlier fixed-0.5 lists). That ratio is real, is
not driven by one ROI (excluding 301.tiff it *rises* to 348:1), but it **varies 15-fold
across images** -- 40:1 on 201.tiff to 601:1 on 002.tiff -- and the first draft quoted three
of seven and omitted both extremes.

More importantly, **the ratio does not survive to any operating point**, and Fix 2 is what
destroys it. Pooled top-K:

| where | ranker | TP | look-alike | unannotated | un:lk |
|---|---|---:|---:|---:|---:|
| full list | -- | -- | 166 | 28 284 | 170:1 |
| top-K | score | 87 | 19 | 278 | 14.6:1 |
| **top-K** | **chromatin** | **187** | **46** | **151** | **3.3:1** |

Chromatin density promotes look-alikes almost as hard as it promotes mitoses --
AUC(look-alike vs unannotated) 0.880-0.965 against AUC(TP vs unannotated) 0.799-0.997 -- so
the surviving false positives concentrate into exactly the pathologist-rejected class. In
absolute terms the trade is still strongly favourable (TP 87 -> 187, unannotated FP
278 -> 151, look-alike FP 19 -> 46), but *"you don't need to solve the look-alike problem"*
is wrong, and it is wrong specifically because this fix works.

### What chromatin density actually measures

`non_human_findings` means "matched no annotation" -- background and stroma included, not
"ordinary nuclei". Split the negatives and the statistic's character shows:

| image | n mitotic | AUC TP vs unannotated | AUC TP vs **look-alike** |
|---|---:|---:|---:|
| 201.tiff | 17 | 0.997 | 0.833 |
| 246.tiff | 115 | 0.988 | 0.831 |
| 301.tiff | 217 | 0.973 | 0.797 |
| 405.tiff | 13 | 0.969 | 0.700 |
| 002.tiff | 8 | 0.994 | 0.646 |
| 506.tiff | 11 | 0.925 | 0.527 |
| 350.tiff | 3 | 0.799 | **0.269** |

**It is much better at dense-object detection than at mitosis-vs-mimic.** But the table
above is measured on the *pipeline's own* candidate set, which captures only 29-48% of the
annotated look-alikes -- a selection correlated with the feature under test. Re-measured on
the `nucleus_blobs` candidate set, which captures 98-100% of look-alikes and 96.6-100% of
mitoses, the bias is worth 0.06-0.10 of AUC (`lookalike_auc_unbiased.py`,
`results/lookalike_auc_unbiased.csv`):

| ROI | AUC(TP vs look-alike), unbiased | 95% CI | n TP vs n look-alike | biased estimate above |
|---|---:|---|---|---:|
| 246.tiff | **0.768** | [0.707, 0.830] | 110 v 123 | 0.831 |
| 301.tiff | **0.693** | [0.635, 0.751] | 216 v 105 | 0.797 |

Both decisively exclude 0.5, so the honest reading is **0.69-0.77 against 0.94-0.99 for
TP-vs-unannotated -- a real but much weaker signal**, not "unable to separate mitosis from
mimic" (too harsh) and not "decent" (too generous). Four of the seven rows in the table above
(405, 002, 506, 350) rest on 2-11 true positives with CIs spanning 0.5; the ordering among
them is noise.

Two things this rules out. The window statistic is **not** doing something cleverer than a
component mean on the hard class -- chromatin density 0.693 vs the blob component-mean 0.687
on 301.tiff, statistically indistinguishable -- so its advantage is confined to the easy
class, where the component mean is markedly better (0.986 vs 0.935 on 301.tiff) -- which is
exactly why `nucleus_blobs` ranked by its own native score outperforms the same detector
re-ranked by chromatin density. And because both available statistics are stuck
at 0.69-0.77 with 105-123 look-alikes available to measure against, the gap for a texture
feature to close is genuine and well powered.

Median OD on 301.tiff: 0.137 (mitotic), 0.114 (look-alike), 0.090 (unannotated).

## Still open

- **`recall@K` is the wrong headline metric here and is reported only as a secondary check.**
  Its K is the ROI's own mitotic count -- the quantity being sought -- so it cannot define an
  operating budget in practice. The workload curve can. (The first draft used recall@K as
  evidence where it agreed and dismissed it where it did not; that contradiction is removed.)
- **z = 2.5 is the only clean level in the swept set, and it is not tuned.** A deeper
  extraction floor than 0.25 would let z = 1.5-2.0 be measured properly; the z >= 2 row
  (recall 0.935, spread 1.8x) is promising and currently contaminated by one floor-limited arm.
- **Look-alikes are enriched 40x at the operating point but are not the majority.**
  Top-K is TP 187 / look-alike 46 / unannotated 151, so look-alikes are **23.4%** of false
  positives, up from 0.58% on the full list. Perfect look-alike suppression moves FPs only
  197 -> 151. An earlier draft of this log called them "the majority problem"; that was
  wrong, and it matters because texture work aimed at this tier is worth less than choosing
  the score floor properly (z = 2.5 -> 2.0 moves 301.tiff's recall@1000 by 0.185).
- **Four of seven ROIs have n < 15 mitotic figures** (405 n=13, 506 n=11, 002 n=8, 350 n=3).
  The 002.tiff showpiece rests on n=8. Only 301.tiff (n=217) and 246.tiff (n=115) carry real
  statistical weight; both win (6.4x, 11.9x median over seeds), which is the reassuring part,
  and both win by less than the small-n ROIs, which is the honest part.
- **ROI selection is still the densest-per-domain rule**, documented in
  `select_domain_images` as optimistic for whatever is being measured.
- **Nothing here has been tried outside the 7 downloaded ROIs**, or with more than one seed
  used simultaneously.

## What the audit changed

An independent adversarial audit re-derived every number from the raw images and attacked
each claim. Corrections it forced, all applied above:

1. **A real defect in the statistic.** `channels.to_hematoxylin` saturates dense chromatin at
   its 255 ceiling, tying 66-77 detections on 002.tiff and 506.tiff and leaving `sort_values`'
   unstable quicksort to order the whole top-100 arbitrarily. The headline "002: 315 -> 50"
   was one draw from a 10-76 lottery. Fixed with `chromatin.hematoxylin_od` (unclipped) and a
   stable sort; every number above is post-fix, and the result improved (50 -> 8 on 002.tiff,
   median 12.6x -> 19.4x).
2. **Fix 1's recall claim reverses at z = 3** (0.387 vs the constant's 0.544). "The z-floor
   raises recall" was false as stated; it is true of z = 2.5 specifically. Rewritten.
3. **"You don't need to solve the look-alike problem" was wrong at every budget.** 170:1 is a
   full-list property; at top-K under chromatin ranking it is 3.3:1. Rewritten, with the
   AUC(look-alike vs unannotated) mechanism that explains it.
4. **`non_human_findings` is "matched no annotation", not "ordinary nuclei"**, and
   AUC(TP vs look-alike) is 0.269-0.833 against 0.799-0.997 for TP vs unannotated. The
   statistic is largely a dense-object detector. Relabelled and both columns now reported.
5. **`coverage_frac` was never quoted** despite `evaluate.py` instructing that it be read
   before any full-list number. Added.
6. **Reproducibility gaps.** `results/od_workload_ab.csv` -- the Fix-2 headline -- had no
   generator in the repo; `od_workload_ab.py` now produces it. `od_experiment.py` saved
   detections for z = 1.0 rather than the reported z = 2.5; fixed.
7. **Unsourced numbers.** A median-OD triple (165/144/99) matched no computable statistic and
   is replaced with measured values (0.137/0.114/0.090). "Next-best `solidity` 0.15-1.05" was
   wrong -- solidity is runner-up in only 2 of 7 domains.
8. **Overclaims trimmed.** "Blind to darkness" -> much weaker ranker (the score still reaches
   AUC 0.61-0.81). "Wins on all seven" -> 350.tiff is a 1.07x tie. Sample sizes now flagged
   for all four ROIs with n < 15, not two.
9. **An internal contradiction removed**: recall@K was used as Fix 2's headline and dismissed
   as invalid in the same document.

Claims that survived the attack unchanged: the z-of-0.5 table (re-derived cell by cell from
the images), the affine-invariance proof, `mean_intensity` as the top feature in all seven
domains, the 170:1 full-list ratio and its robustness to dropping 301.tiff, the
one-match-many-arms shortcut (exact -- reproduces the saved fixed-0.5 counts on all 7 ROIs),
the seed RNG construction, and -- checked directly -- that re-sorting before greedy matching
introduces **no bias**: zero detections lie within the match radius of two ground-truth
objects on any ROI, and an order-independent GT-centric metric reproduces the workload
numbers exactly on all 14 image x ranker cells.
