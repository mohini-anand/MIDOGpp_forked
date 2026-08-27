> **SUPERSEDED — 2026-08-27.** The numbers in this log come from a run whose evaluation
> layer carried five bugs, all of them catalogued in the "known issues" section at the
> foot of this file. Those are now fixed and the experiment re-run; see
> `2026-08-27-find-and-suppress-corrected-rerun.md` for the current numbers and
> `2026-08-25-find-and-suppress-code-review.md` for the full audit.
>
> Kept as-is because the reasoning, the port rationale, and the reference-code bug
> catalogue are all still correct — only the measured values changed. **Do not quote
> figures from this file.** One conclusion actually reverses under the corrected code
> (the fusion comparison; see §4a of the review log).

# `find_and_suppress` on MIDOG++: single-pass, one seed per domain

Date: 2026-08-25
Scope: exploratory. One ROI per tumour domain (7), one randomly-chosen mitotic-figure
seed each, a single forward pass. No iterative bootstrapping, no multi-seed variance
estimate, no hematoxylin variant.
Code: `midog_utils/`, driven by `find_and_suppress_midog.ipynb`.

> **Correction (2026-08-26).** A code review
> (`2026-08-25-find-and-suppress-code-review.md`) found that the full-list results below
> are invalidated by coverage saturation, and several narrative claims are wrong. The
> corrected readings are marked inline. The top-K results and the diagnosis are unaffected.

## The question

MIDOG++ annotates every mitotic figure inside a ~2 mm² ROI, which means the dataset can
answer a question the hydrogel reference code never could: hand the algorithm **one**
annotation, let it search the rest of the ROI, and score every detection against the
ground truth that was held back. How far does one human click get you?

## Setup

- **Seed**: one category-1 annotation per image, drawn with a fixed RNG, excluded from
  the evaluation set. Annotations within 36 px of the ROI border are ineligible (no
  rotation-safe template can be read there).
- **Images**: per domain, the downloaded ROI with the most category-1 annotations.
  This rule drops `001.tiff` on its own — it has zero category-1 annotations.
- **Matching**: 12 rotations × 2 flips × **one scale** (the 51 px annotation box itself,
  odd-sized so the click is the centre pixel) = 24 augmentations of `TM_CCOEFF_NORMED`,
  fused by element-wise max into one response map in image-centre coordinates, then
  peak-picked once. Multi-scale was measured and rejected — see below.
- **Detection list**: every local maximum scoring ≥ 0.25 that survives NMS, with no top-K
  truncation. 13 000–23 000 detections per ROI.
- **Suppression**: score-ordered greedy NMS at 25 px on detection centres.
- **Scoring**: greedy descending-score one-to-one matching at MIDOG's own radius of
  7.5 µm (≈30–33 px depending on scanner). Three detection buckets: `TP` (matched a
  category-1 GT), `FP_lookalike` (matched a category-2 GT — a structure a pathologist
  examined and rejected), `FP_unannotated` (matched nothing).

## Ported, not reused — and why

The reference `find_and_supress` was written for bright, near-isotropic particles on a
uniform dark background. Two of its bugs would have produced a wrong answer that reads
as a scientific finding rather than as a defect:

- **`shortlist_augmentations` (bbox_tuning.py:323-324) kills the rotation search.**
  `resized = cv2.resize(rotated, ...)` is overwritten on the very next line by
  `resized = flipped.astype(np.float32).ravel()`. All four rotations therefore become
  byte-identical, cosine similarity is 1.0, and the greedy de-duplication collapses them
  to one variant per flip. Rotation-invariant matching silently does nothing. Fatal on
  mitotic figures, which have no canonical orientation.
- **NMS is not score-ranked.** `nms_with_area` walks the DataFrame index and keeps the
  *last* element, so the survivor of a cluster is whichever box sits furthest right. The
  `match score` computed in `match_one_combination` is dropped at bbox_tuning.py:798-801
  and never ranks anything — which also means no FROC and no threshold sweep are
  possible from the reference's output.

Also fixed or dropped: overlap computed as `intersection / area[j]` with mixed ±1 pixel
conventions (can exceed 1.0); unseeded `sample(frac=1)` template choice and unseeded
`random.sample` truncation that discards the *best* matches; `df_temp` snapshotted once
in `iterative_find_and_supress` so its convergence test compares against the original
length forever; a `NameError` on undefined `extension` in its save path.

Deliberately not run: `refine_boundary`, `remove_unstable`, and the aspect-ratio / PSNR
filters. Evaluation here is point-based so box geometry moves no metric, and the
stability check passes 0/9 objects on H&E (see the bbox-tuning log of the same date) —
keeping it would delete essentially every detection.

**No image masking.** The reference blacks out existing annotations before correlating.
On H&E that creates zero-variance windows, where `TM_CCOEFF_NORMED` is undefined, and
hard edges that rotated templates correlate against. The seed's own detection is dropped
afterwards with a tight 5 px radius instead — not the match radius, because the minimum
spacing between two MIDOG++ annotations is 27 px and a match-radius drop would delete a
legitimate detection of a neighbouring ground-truth object.

## Dataset details confirmed along the way

- `ResolutionUnit` is **3 (centimetre)** for `505.tiff`/`506.tiff` and 2 (inch) for every
  other image. Assuming inches gives 0.578 µm/px for those two — implying a 13 mm² ROI —
  instead of the correct 0.227 µm/px and 2.02 mm². Asserting that every ROI computes to
  ~2 mm² is what catches this.
- `505/506` are also the only 7-level LZW pyramids; everything else is a flat
  uncompressed single page. Reading `series[0].levels[0]` handles both.
- `openslide` is not installed in the environment that has a working `cv2`
  (`~/anaconda3/bin/python`); the system `python3` has a `cv2` compiled against NumPy 1.x
  that fails to import under NumPy 2.3.4. `tifffile` avoids the problem entirely and is
  the better fit anyway — these are flat 2 mm² ROIs, not pyramidal WSIs.

## Verification gate

Multiscale + rotation means a response-map peak must be mapped back using **half of its
own template's size**: a 31 px template yields an `(H-30, W-30)` map, a 51 px template an
`(H-50, W-50)` map. An off-by-half-template error shifts every detection 10–25 px, which
at a ~33 px match radius destroys recall and presents as *"template matching doesn't work
on H&E"*.

Each of the 72 augmentations was planted into a noise canvas at a known centre and the
peak of the **fused** map required to come back within 1 px. Asserting on the fused map
rather than per-augmentation maps is the point — per-augmentation recovery can pass while
the fusion offsets are wrong. Result: **72/72 recovered at exactly 0 px, score 1.0000.**

## Results

Seven images, one seed each, 24 augmentations, full detection list.
Matching takes 37–59 s per ROI.

### The full list — every detection, no budget

| domain | image | detections | TP | FP look-alike | FP unannotated | precision | **mitotic GT missed** |
|---|---|---:|---:|---:|---:|---:|---:|
| mast cell tumor | 301 | 20 031 | 217 | 108 | 19 706 | 1.08% | **0** |
| lymphosarcoma | 246 | 23 143 | 115 | 123 | 22 905 | 0.50% | **0** |
| lung cancer | 201 | 19 492 | 17 | 49 | 19 426 | 0.09% | **0** |
| melanoma | 506 | 18 928 | 11 | 11 | 18 906 | 0.06% | **0** |
| soft tissue sarcoma | 405 | 13 370 | 13 | 19 | 13 338 | 0.10% | **0** |
| breast cancer | 002 | 16 605 | 8 | 7 | 16 590 | 0.05% | **0** |
| neuroendocrine | 350 | 20 313 | 3 | 16 | 20 294 | 0.01% | **0** |

> **RETRACTED.** This table showed `all_mitotic_gt_missed = 0` on every image and the
> original text read *"every mitotic figure in every ROI is found — this was never a
> detection failure"*. That conclusion is wrong. At the 0.25 score floor the detection list
> is dense enough that **91–97% of each ROI lies within a match radius of some detection**,
> so "was this annotation found" is answered by geometry, not by the response map.
>
> The control settles it: `random_in_tissue` — uniform noise at the same list length —
> reaches full-list mitotic recall of **1.000 on 002.tiff and 350.tiff**, and 0.571–0.923
> elsewhere. A recall of 1.000 here is a near-certainty available to any method, not a
> property of template matching.
>
> The precision column is likewise a property of the 0.25 floor rather than of the method:
> at threshold 0.70 the list is 1 detection on 301.tiff and 2 601 on 405.tiff. There is no
> "natural operating point" at which these numbers were measured.
>
> What survives: the **relative** counts at a common budget (the `recall@K` table below),
> and the threshold sweep, which shows how fast the list thins. What does not: any claim
> that all mitotic figures were detected, and `all_precision_mitotic` as a characterisation
> of the method. Re-measuring at a threshold where coverage is not saturated (~0.60 on
> 301.tiff, giving ~4% coverage) is the fix.

Pooled over all seven images (384 mitotic GT), as the score floor moves:

| score ≥ | detections | TP | FP look-alike | FP unannotated | precision | recall |
|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 131 882 | 384 | 333 | 131 165 | 0.29% | 1.000 |
| 0.40 | 88 911 | 368 | 307 | 88 236 | 0.41% | 0.958 |
| 0.50 | 46 817 | 301 | 223 | 46 293 | 0.64% | 0.784 |
| 0.60 | 13 538 | 90 | 71 | 13 377 | 0.67% | 0.234 |
| 0.70 | 3 177 | 24 | 20 | 3 133 | 0.76% | 0.062 |
| 0.80 | 419 | 6 | 7 | 406 | 1.43% | 0.016 |
| 0.90 | 1 | 0 | 0 | 1 | 0% | 0.000 |

**Precision never exceeds 1.5% anywhere in the operating range.** Recall holds at 0.958 down
to a threshold of 0.4 and only collapses past 0.6. There is no cutoff at which this method
is usable — the entire curve sits in a regime where 99+ detections in 100 are wrong.

### At a budget — recall@K

| domain | image | K | find_and_suppress | nucleus_blobs | random_in_tissue |
|---|---|---:|---:|---:|---:|
| lymphosarcoma | 246 | 115 | **0.191** | 0.357 | 0.000 |
| breast cancer | 002 | 8 | **0.125** | 0.000 | 0.000 |
| mast cell tumor | 301 | 217 | **0.097** | 0.599 | 0.005 |
| soft tissue sarcoma | 405 | 13 | **0.077** | 0.077 | 0.000 |
| lung cancer | 201 | 17 | **0.059** | 0.235 | 0.000 |
| melanoma | 506 | 11 | **0.000** | 0.000 | 0.000 |
| neuroendocrine | 350 | 3 | **0.000** | 0.000 | 0.000 |

Going single-scale improved `recall@K` on three images and left four unchanged; 246 nearly
doubled (0.096 → 0.191). It still loses to the blob baseline on the three images where that
baseline is non-zero, though the gap narrowed — 246 went from 3.7x worse to 1.9x worse.

### Mitosis versus mimic — the single-scale version does discriminate, a little

| image | recall@K | attraction@K | gap |
|---|---:|---:|---:|
| 246 | 0.191 | 0.065 | **+0.126** |
| 301 | 0.097 | 0.046 | +0.051 |
| 405 | 0.077 | 0.000 | +0.077 |
| 201 | 0.059 | 0.020 | +0.039 |
| 002 | 0.125 | 0.143 | −0.018 |
| 506 | 0.000 | 0.083 | −0.083 |

This is a genuine change from the multi-scale run, where the two were equal on every image.
Dropping the 31 px template — which on 246 *preferred* look-alikes — recovers a real
mitosis-over-mimic preference on four of six images with non-zero rates. It is small, and it
reverses on 002 and 506, but it is no longer flat.

### FROC

Sensitivity at 8 FP/mm² is 0.00–0.375, at 64 FP/mm² 0.06–0.625, against roughly 0.7 for a
trained MIDOG detector. Stated only to fix the scale.

The expert-agreement split is degenerate in this configuration: `recall_by_agreement` reads
the full-list `found` flag, which is saturated (above), so it reports 1.000 for both
unanimous and contested annotations everywhere. It needs the same `k` parameter
`lookalike_attraction_rate` already has. Note also that for category 1, `unanimous` is
*exactly* equivalent to `n_votes == 2` — the only two label multisets are `(1,1)` and
`(1,1,2)` — so it is a proxy for "needed a third reader", not independent signal.

## Multi-scale fusion measured and rejected

Found while checking why the top detections were almost all from the smallest template.

`TM_CCOEFF_NORMED` is a correlation over the template's own pixels, so a smaller template
has a wider null distribution and produces a higher *maximum* by chance. Taking a raw
element-wise max across 31 / 41 / 51 px therefore hands the win to the smallest one: the
51 px template — the annotation box itself — won only **4.7% pooled** of top-K detections
against a uniform expectation of 33%.

Probing the response map at annotated locations quantifies the cost
(`results/fs_fusion_variants.csv`). AUC = P(a random mitotic figure outscores a random
ordinary nucleus); *discrimination* = that AUC minus the same figure for look-alikes, i.e.
how much of the signal is mitosis-specific rather than generic chromatin:

| image | fusion | templates | AUC mitosis | AUC look-alike | discrimination | median mitosis rank |
|---|---|---:|---:|---:|---:|---:|
| 002 | max over 31/41/51 | 72 | 0.906 | 0.869 | +0.037 | 806 |
| 002 | z-normalised max | 72 | 0.928 | 0.886 | +0.042 | 619 |
| 002 | **51 px only** | 24 | **0.940** | 0.904 | +0.036 | **518** |
| 002 | 31 px only | 24 | 0.836 | 0.725 | +0.110 | 1413 |
| 246 | max over 31/41/51 | 72 | 0.841 | 0.838 | **+0.003** | 3178 |
| 246 | z-normalised max | 72 | **0.866** | 0.834 | +0.032 | **2672** |
| 246 | **51 px only** | 24 | 0.853 | 0.792 | **+0.061** | 2933 |
| 246 | 31 px only | 24 | 0.793 | 0.824 | **−0.031** | 4146 |

Readings:

- **The raw max is the worst option on the dense image.** On 246 its discrimination is
  +0.003 — to three decimal places the multi-scale pipeline is exactly as attracted to
  pathologist-rejected mimics as to real mitotic figures.
- **The 31 px template is the culprit.** On 246 its discrimination is *negative*: alone, it
  prefers look-alikes to mitotic figures. Letting it win the max by chance is what drags the
  fused result down.
- **The single 51 px template is the best practical choice.** Best AUC and best rank on 002,
  best discrimination on 246, and a third of the compute. Honest caveat: z-normalised
  multi-scale edges it on raw AUC and rank on 246 (0.866 / 2672 vs 0.853 / 2933), so it is
  not a clean sweep — but it costs 3x more for a worse mitosis-vs-mimic gap.

`FSConfig(scales=(1.0,))` is now the default. Multi-scale and `scale_normalize` remain
available and both pass the plant-and-recover gate 72/72.

## Diagnosis — why the true-positive count is so low

With the budget removed, every mitotic figure is found. So the question is not "why does it
miss them" but **"why can it not rank them"**. Probing the fused response map directly at
annotated locations (`results/fs_score_probe.csv`) answers it.

Median fused score on 002.tiff, by what is at that location:

| mitotic GT | look-alike GT | ordinary nucleus | random tissue |
|---:|---:|---:|---:|
| 0.756 | 0.710 | 0.491 | 0.308 |

A clean ordering: mitotic > look-alike > ordinary nucleus > background. As a ranking
statistic, AUC = P(a random mitotic figure outscores a random ordinary nucleus) = **0.940**,
against 0.904 for look-alikes.

### An AUC of 0.94 is still not enough

002.tiff holds **8 mitotic figures among 8 604 nuclei** — a base rate of **0.093%**. The
expected rank of the median mitotic figure among all nuclei is `(1 − AUC) x n_nuclei` =
**518**, against a top-K budget of **8**. Sixty-five times too deep.

The tail is where this bites. `recall@K` reads only the extreme upper end of the score
distribution, and there the two populations overlap almost exactly: the **99.9th percentile
of ordinary nuclei scores 0.789** against a **90th percentile of 0.788 for mitotic figures**.
With 8 604 nuclei the top 0.1% is about 9 objects — the entire budget. Those slots fill with
extreme-tail ordinary nuclei before a single mitosis appears.

### The ranking is good, just not good enough

Median rank of a true positive within each image's own detection list:

| image | K | detections | median TP rank | as a fraction |
|---|---:|---:|---:|---:|
| 002 | 8 | 16 605 | 57 | 0.003 |
| 405 | 13 | 13 370 | 547 | 0.041 |
| 350 | 3 | 20 313 | 915 | 0.045 |
| 246 | 115 | 23 143 | 1 336 | 0.058 |
| 201 | 17 | 19 492 | 1 190 | 0.061 |
| 506 | 11 | 18 928 | 2 034 | 0.107 |
| 301 | 217 | 20 031 | 2 916 | 0.146 |

Under a score carrying no information the fraction would be 0.5. Observed: **0.003–0.146**.
The ranking is genuinely informative — mitotic figures land in the top 5–15% of detections,
and in the top 0.3% on 002. But `recall@K` needs them in the top 0.05–1.1%, which is one to
two orders of magnitude tighter than what the score delivers.

(This is markedly better than the multi-scale configuration, where the same fractions sat at
0.3–0.5, i.e. near chance. Dropping the 31 px template is what recovered it.)

### Two causes, in order of size

1. **Base rate.** Mitotic figures are 0.09–1.1% of the nuclei in a 2 mm² ROI. Ranking a
   handful of objects above ~10 000 similar-looking distractors demands separation far out
   in the tail. An AUC of 0.94 leaves ~518 nuclei ahead of the typical mitosis. This is the
   dominant term and no tuning changes it.
2. **The score is only weakly mitosis-specific.** AUC 0.940 for mitotic figures against
   0.904 for look-alikes — a real gap, but a small one. `TM_CCOEFF_NORMED` against a
   chromatin template measures "dark, compact, similarly textured", which is precisely the
   property that makes a look-alike a look-alike.

## This is not a bug

Worth recording explicitly, because a low number from a ported pipeline invites the
assumption that something is mis-wired:

- The coordinate round-trip gate passed **72/72 augmentations at exactly 0 px offset**,
  asserted on the fused map.
- The seed's own detection comes back at score 1.0000 on every image (`seed_self_score`),
  which means the matcher and the coordinate mapping agree on real data too.
- The rotation search is verifiably alive here, unlike in the reference.

The pipeline does exactly what it is supposed to. The method is what falls short.

## Cost

100–172 s per image for the 72 correlations over a ~35-megapixel ROI, plus ~20 s for the
blob baseline. About 20 minutes for the seven-image sweep on 12 threads.

## Corrections carried from the code review

Verified against the saved results before recording:

- **Coverage saturation invalidates the full-list tables** (retraction above). Confirmed by
  the `random_in_tissue` control reaching 1.000 full-list recall on two ROIs.
- **`nucleus_blobs` silently truncates at 20 000** (`baselines.py`). The cap binds on
  301.tiff (20 000 delivered vs 20 031 requested) and 246.tiff (20 000 vs 23 143), so the
  "blob detector truncated to the matcher's budget" claim is false for exactly the two
  images where the blob baseline scores best. It also means `score_probe`'s `n_nuclei` reads
  20 000 on 246.tiff against a population truncated to the *darkest* blobs — i.e. biased to
  keep the hardest competitors. `mitotic_base_rate_pct`, `p99.9_nucleus` and
  `median_mitosis_rank` for 246.tiff are affected.
- **`tissue_mask` is not separating tissue from background.** These ROIs are near-solid
  tissue (pixels above gray 230 are 0.1–1.4%), so global Otsu splits dark tissue from pale
  tissue, not tissue from slide. `random_in_tissue` is therefore "uniform over the darker
  half" — a *stronger* floor than advertised, so the headline comparison is conservative,
  but the docstring's rationale is false as written. Separately, ~31 of 218 mitotic GT
  centres on 301.tiff fall outside the mask, and `nucleus_blobs` cannot place a candidate on
  them; whether that costs it anything at top-K is untested, because the blob detections are
  not retained.
- **The single-scale decision was drawn from n=1.** `results/fs_fusion_variants.csv` holds
  246.tiff too, and there `fused_z_normalised` beats `scale_1.0_only` on raw AUC (0.866 vs
  0.853) and rank (2672 vs 2933). The log's fusion table shows both; the notebook cell does
  not. The two selection criteria also disagree: on the `discrimination` column,
  `scale_0.6_only` wins on 002.tiff and `scale_1.0_only` on 246.tiff.
- **`max_score` in the metrics CSV is always the seed self-hit** — recorded before self-hit
  removal, so it reads 1.000 on every row instead of the interesting post-removal maximum
  (0.70 on 301.tiff to 0.91 on 405.tiff).
- **`recall@K` is easier on denser ROIs**, because the base rate climbs against a roughly
  constant nucleus count — and the ROI-selection rule picks the densest image per domain. The
  seven numbers are optimistic for their domains, not representative draws.

## Open items

- **Settled since the first pass:** multi-scale fusion was the wrong default and is gone;
  the low `recall@K` is a ranking failure against the base rate, not a detection failure.
- **Recompute the expert-agreement split at a budget.** With the full list it is degenerate
  (recall 1.000 for both unanimous and contested annotations everywhere).
- **Seed variance is unmeasured.** These are one seed per image. With a single 50x50
  template, run-to-run spread from seed choice could plausibly exceed every between-domain
  difference in the tables above. The 5-seed sweep is the next thing to run, and no domain
  ordering here should be quoted until it is.
- Hematoxylin-channel matching is implemented (`FSConfig(channel='hematoxylin')`) but not
  run. It is unlikely to change the conclusion — the blob baseline already uses that
  channel and wins — but it is cheap.
- Iterative bootstrapping is not implemented yet. Given Finding 1, seeding subsequent
  rounds with the matcher's own detections would mostly bootstrap on look-alikes and
  ordinary nuclei. Worth running to quantify the compounding, not as a candidate fix.
- The real conclusion points where the bbox-tightening log already pointed: a method with
  an actual notion of "nucleus" and "mitotic". The blob baseline outperforming template
  matching is evidence for the instance-segmentation tier (StarDist, CellPose, HoVer-Net)
  and the promptable-SAM tier (which would take the existing point annotation directly as
  a prompt), rather than for more correlation tuning.
