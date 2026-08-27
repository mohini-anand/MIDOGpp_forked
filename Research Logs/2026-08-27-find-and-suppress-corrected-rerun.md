# `find_and_suppress` on MIDOG++: corrected re-run

Date: 2026-08-27
Supersedes: `2026-08-25-find-and-suppress-midog-single-pass.md`
Audit that motivated it: `2026-08-25-find-and-suppress-code-review.md`
Code: `midog_utils/`, driven by `find_and_suppress_midog.ipynb`.
Environment: `requirements-midog-utils.txt` (numpy 1.26.4 / opencv 4.8.1, Python 3.11.5).
Verification: `verify_fixes.py` — 42 checks, all passing.

## What changed

Twenty-two fixes, of which six move numbers:

| | before | after |
|---|---|---|
| `recall_by_agreement` | read the full detection list → 1.0 in every cell | budgeted at K |
| `nucleus_blobs` cap | silently truncated at 20 000 | uncapped; delivers 7 577–21 059 |
| `tissue_mask` | global Otsu → 0.38–0.73 of the ROI, 31/218 mitoses excluded on 301.tiff | fixed cut at gray 220 → 0.78–1.00, **0 excluded on all 14 ROIs** |
| NMS radius | fixed 25 px against a 29.6–33.1 px match radius | the image's own match radius |
| seeds | `default_rng(seed)` rebuilt inside the loop → the same ~80th-percentile draw 7× | per-image stream from `(seed, image_id)` |
| `score_probe` | nucleus set capped and containing the annotated objects; AUC with ties as losses; the mean reported as "median" | annotations excluded, Mann-Whitney ties, mean and median both reported |

The rest — dead `threshold_metrics`, the `// 2` in the assignment cross-check, `max_score`
being the seed's own hit, `groupby().first()`, the `pick_seed` border off-by-one, the
overlay circle radius — are corrections that did not move a published number.

## Headline: `recall@K`, equal budget

K = that ROI's mitotic count minus the seed.

| domain | image | K | find_and_suppress | nucleus_blobs | random_in_tissue |
|---|---|---:|---:|---:|---:|
| canine cutaneous mast cell tumor | 301 | 217 | 0.097 | **0.599** | 0.018 |
| canine lymphosarcoma | 246 | 115 | 0.043 | **0.330** | 0.009 |
| canine lung cancer | 201 | 17 | 0.000 | **0.235** | 0.000 |
| canine soft tissue sarcoma | 405 | 13 | 0.000 | 0.000 | 0.000 |
| human breast cancer | 002 | 8 | **0.125** | 0.000 | 0.000 |
| human melanoma | 506 | 11 | 0.000 | 0.000 | 0.000 |
| human neuroendocrine tumor | 350 | 3 | 0.000 | 0.000 | 0.000 |

**The finding is unchanged and now better supported.** The "any dark blob" detector beats
one-shot template matching at equal budget on the three ROIs where either method scores at
all, ties on three, and loses on one. It does so *despite* the tissue-mask fix having
removed a handicap that was previously suppressing it — the blob detector was being denied
14% of the mitoses on 301.tiff and still won there.

Absolute numbers fell against the superseded run (246: 0.191 → 0.043; 201: 0.059 → 0.000).
That is the seed-choice variance the old code was concealing by drawing the same quantile
from every image, not a regression. It vindicates the original log's warning to treat the
per-domain ordering as provisional, and it argues that the multi-seed sweep is now the
next thing worth running rather than an optional refinement.

## Coverage saturation, now measured rather than inferred

`coverage_frac` — the fraction of the ROI within a match radius of some detection — is
**0.876–0.948**. At the 0.25 score floor the detection list tiles the ROI more finely than
the metric can resolve, so every un-budgeted recall is largely a statement about list
length. The floor is explicit in the same table: `random_in_tissue`, uniform noise, reaches
full-list mitotic recall of 0.46–1.00 at coverage 0.69–0.78.

The saturation is a property of the score floor, not the method. From the sweep, 301.tiff
goes 16 681 → 7 277 → **404** detections at 0.25 → 0.50 → 0.60. How arbitrary a fixed floor
is across images shows in the same row: at 0.60 there are 6 detections on 246.tiff and
6 479 on 405.tiff. Choosing a defensible operating point is the outstanding methodological
decision; the fixes make the problem visible but do not settle it.

## The agreement split, now that it says something

Previously 1.0 in every cell. Budgeted at K:

| image | n unanimous | recall | n contested | recall |
|---|---:|---:|---:|---:|
| 301 | 163 | 0.104 | 54 | 0.074 |
| 246 | 99 | 0.030 | 16 | **0.125** |
| 002 | 7 | 0.143 | 1 | 0.000 |

Too small to lean on, and it does not point one way — 301 recovers unanimous annotations
slightly more often, 246 recovers contested ones four times more often. Worth carrying into
the multi-seed sweep, where the counts become large enough to mean something.

## A conclusion the fixes reversed

The audit flagged that the notebook's "the single 51 px template wins outright" was
contradicted by its own saved CSV, where `fused_z_normalised` beat `scale_1.0_only` on
246.tiff (AUC 0.866 vs 0.853). **That is no longer true, because the CSV was wrong.** Its
AUCs were computed against a nucleus population capped at 20 000 and still containing the
annotated mitoses and look-alikes, with ties scored as losses. Recomputed correctly:

| image | variant | templates | AUC mitosis | AUC look-alike | discrimination | median rank |
|---|---|---:|---:|---:|---:|---:|
| 002 | fused_max_multiscale | 72 | 0.895 | 0.810 | 0.084 | 334 |
| 002 | fused_z_normalised | 72 | 0.942 | 0.841 | **0.101** | 87 |
| 002 | **scale_1.0_only** | 24 | **0.953** | 0.863 | 0.089 | **85** |
| 002 | scale_0.6_only | 24 | 0.734 | 0.667 | 0.067 | 2297 |
| 246 | fused_max_multiscale | 72 | 0.812 | 0.799 | 0.013 | 3472 |
| 246 | fused_z_normalised | 72 | 0.877 | 0.816 | 0.061 | 1656 |
| 246 | **scale_1.0_only** | 24 | **0.879** | 0.799 | **0.080** | **1103** |
| 246 | scale_0.6_only | 24 | 0.702 | 0.763 | −0.061 | 5661 |

The single 51 px template now wins on AUC and median rank on both images, so the default is
better justified than the audit concluded. Two caveats survive: the AUC margin on 246.tiff
is 0.002, so "wins" there rests on the rank; and `discrimination` — the criterion the
experiment's central claim actually turns on — prefers `fused_z_normalised` on 002.tiff.
The mechanism the default rests on is confirmed either way: `scale_0.6_only` is worst on
both images, the raw multi-scale max inherits most of that deficit despite having the 51 px
template available, and z-normalising before fusing recovers nearly all of it.

This is the clearest argument in the whole exercise for regenerating a comparison table
from fixed code rather than citing a stored one.

## A second claim that reversed

The superseded notebook read: *"the top 0.1% of ordinary nuclei outscore 90% of mitotic
figures"*, from `p99.9_nucleus > p90_mitotic`. Once the annotated objects are removed from
the nucleus population — they were the highest-scoring "nuclei" in it — the inequality holds
on 246.tiff (0.556 vs 0.517) and **reverses** on 002.tiff (0.651 vs 0.664).

The base-rate problem itself is untouched: mitotic figures are 0.09% and 0.55% of nuclei,
and the median mitotic figure ranks 85th of 8 796 ordinary nuclei on 002.tiff against a
budget of K=8, and 1 103rd of 20 609 on 246.tiff against K=115. Only that particular way of
phrasing it was an artefact.

Note also how far apart mean and median rank are — 418 vs 85 on 002.tiff. The superseded
run reported the mean under the name "median".

## Still open

- **The operating point.** 0.25 is a config value, not a chosen threshold. `coverage_frac`
  now makes the consequence visible; picking a better floor is a research decision.
- **One seed per image.** Seven draws, seven images. The spread between the two runs is the
  argument for the 5-seed sweep.
- **ROI selection is optimistic.** Densest ROI per domain makes `recall@K` easier; these are
  upper-ish bounds per domain, now documented in `select_domain_images`.
- **`nms_radius` vs `peak_min_distance`.** NMS now uses the match radius, but
  `peak_min_distance = 7` still governs the pre-NMS candidate thinning. It does not affect
  correctness — NMS is the real de-duplicator — but the two are worth reconciling.
