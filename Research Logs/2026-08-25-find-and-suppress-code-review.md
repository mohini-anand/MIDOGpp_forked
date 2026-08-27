> **RESOLVED — 2026-08-27.** Every finding below has been fixed, the experiment re-run,
> and the result recorded in `2026-08-27-find-and-suppress-corrected-rerun.md`. Fixes are
> verified by `verify_fixes.py` (42 checks). Per-finding status is at the foot of this file.
>
> **One finding in this review is wrong, and it is §4a.** The correction is inline there.

# Code review: `find_and_suppress_midog.ipynb` + `midog_utils/`

Date: 2026-08-25
Scope: correctness of the code, the assumptions behind it, and the fit between the
notebook's narrative and the numbers in `results/`.
Reviewed: all 36 notebook cells and all 9 modules in `midog_utils/`.

Verdict: the signal-processing core is correct — the coordinate handling, the peak
bound, the NMS ordering, the greedy-once matching and the resolution-tag handling all
check out, and the reference-code bugs the port set out to fix really are fixed. The
problems are concentrated in the **evaluation layer** and in **claims that the saved
results do not support**. One of them invalidates six of the notebook's tables.

---

## 1. Coverage saturation makes every full-list metric vacuous

**The finding.** The detection list is dense enough that the whole ROI is covered by
match radii, so "did the search find this annotation" is answered by geometry rather
than by the response map.

Measured, by sampling every 16th pixel of each ROI and asking the distance to the
nearest detection:

| ROI | detections | match radius | % of ROI within radius of *some* detection |
|---|---:|---:|---:|
| 002.tiff | 16 605 | 33.1 px | **91.1 %** |
| 506.tiff | 18 928 | 33.0 px | **92.6 %** |
| 350.tiff | 20 313 | 33.0 px | **95.9 %** |
| 246.tiff | 23 143 | 30.2 px | **97.3 %** |

This is a design consequence, not an accident, and the driver is the **`score_threshold
= 0.25` floor**, not the NMS radius. The cell-24 sweep shows it directly: on 301.tiff the
list is 20 031 detections at 0.25, 7 538 at 0.50 and **404 at 0.60**. At 404 detections
the covered area is ~1.1 Mpx against a 31 Mpx ROI — about 4 % — and every full-list
metric becomes meaningful again. `nms_radius = 25` being *smaller* than the 29.6–33.1 px
match radius contributes (surviving detections are packed closer than the distance that
counts as a hit) but is second-order: raising it to 33 thins the list by only ~(33/25)²
and still leaves ~37 Mpx of coverage.

**What it invalidates.** `mitotic_recall = 1.000` for all seven ROIs is a certainty, not
a result. The control proves it: `random_in_tissue` — uniform noise — reaches full-list
mitotic recall of **1.000 on 002.tiff and 350.tiff** and 0.57–0.92 elsewhere.

Affected: cells 18, 19, 21, 22, 32, 33, and the `mitotic_recall` / `lookalike_attraction`
columns of cell 14. Cell 22 tells the reader to "note `all_mitotic_gt_missed`
separately" — that column is 0 for every image by construction.

**1a. `recall_by_agreement` has no budget.** [`evaluate.py:207`] It reads full-list
`gt_out["found"]`, so cells 32–33 report `recall_unanimous = recall_contested = 1.0`
everywhere (and NaN where n=0). The expert-agreement split — a genuinely good idea, and
free from the `labels` votes — currently carries zero information. It needs the same
`k` parameter `lookalike_attraction_rate` already has (`matched_rank >= 0 & < k`).

**1b. The same gap inflates FROC false positives.** With NMS at 25 px and a 33 px match
radius, several detections can sit on one true object; only one can be the TP and the
rest are counted as FP. `nms_radius` is a fixed `25.0` in `FSConfig` while the match
radius is derived per image from `mpp` — so the fix is to derive the NMS radius from
`mpp` too (`ev.radius_px(mpp)`), not to edit the scalar. (I checked: at the top-K budget
it costs nothing in this run — every top-K detection within radius of a GT was matched —
but it is latent.)

**1c. `all_precision_mitotic` is a property of the floor, not of the method.** Cell 20
calls the full list "the raw yield of the method at its natural operating point" and
`full_list_breakdown`'s docstring calls its precision "the honest cost of the method at
its natural operating point". There is no natural operating point here — 0.25 is a
config value, and what it means swings by three orders of magnitude between images. From
the cell-24 sweep, at threshold 0.70: **1 detection on 301.tiff, 2 601 on 405.tiff**.
Cell 11 makes exactly this argument to justify preferring `recall@K` over a fixed
threshold; the `all_*` table then presents fixed-threshold numbers as a headline anyway.

**Fix.** Pick an operating point where coverage is not saturated (0.60 on 301.tiff, per
above) and report there; and where full-list numbers are kept, print them against the
coverage floor — what a matched-length random list achieves, which you already compute.

---

## 2. `nucleus_blobs` silently truncates at 20 000

[`baselines.py:64`] `max_detections=20000` binds on two ROIs and nothing warns.

- **Breaks the equal-budget claim.** Cell 17 says the blob detector is "truncated to the
  same number of candidates the matcher produced". It cannot be: 301.tiff needs 20 031
  and 246.tiff needs 23 143, and only 20 000 exist. Cell 18's `matcher_budget` column
  prints the requested number, not the delivered one.
- **Corrupts the score probe.** `score_probe` uses `len(blobs)` as `n_nuclei`, which
  reads exactly 20 000 for 246.tiff. So `mitotic_base_rate_pct`, `p99.9_nucleus` and
  `median_mitosis_rank` are computed against a nucleus population that is both truncated
  *and* biased to the darkest blobs (the list is sorted by mean hematoxylin before
  truncation) — i.e. the hardest competitors are kept and the easy ones discarded.

Raise the cap or record the delivered count.

---

## 3. `tissue_mask` is not separating tissue from background

[`baselines.py:26`] Global Otsu on inverted grayscale assumes the ROI contains slide
background. These ROIs are near-solid tissue: pixels with `gray > 230` are 0.1–1.4 % of
the image. The Otsu split lands at gray 119–169 and yields `tissue_fraction` of
0.38–0.73 pre-morphology — it is separating **dark tissue from pale tissue**.

The docstring's rationale is therefore false as written: "uniform-over-ROI points land
in white space, which a real detector never does". There is essentially no white space
to avoid. `random_in_tissue` is not a uniform floor; it is "uniform over the darker
half", i.e. biased toward chromatin — which makes it a stronger floor than advertised
(conservative, but mislabelled).

Reconstructing the mask (cv2 is not runnable in this env — see §7 — so this uses skimage
Otsu and morphology; the result is stable across square-9 and disk-4 structuring
elements, 31 vs 32) gives **31–32 of 218 mitotic GT centres on 301.tiff falling outside
the mask**, at a median gray of 150 against an image median of 118.7. Because
`nucleus_blobs` computes `binary = (h_chan > thr) & mask`, it cannot place a candidate
*on* those objects. Its full-list recall of 0.986 does not contradict this — an in-mask
blob within 30 px claims them anyway, which is §1 again. Whether it costs anything at
the top-K budget is **untested**: the blob detections are not retained, so there is no
way to check whether any of the 31 excluded GT sit among its top-K misses, and 301.tiff
is where the blob baseline scores its best `recall@K` (0.599). Retaining
`out["baseline_blobs"]` to CSV would settle it.

Direction matters for the headline: this **weakens** the baselines, so "the blob
detector beats the matcher on 3 of 7 ROIs" is if anything conservative.

Also note `MORPH_OPEN` with a 9 px element erases isolated dark objects smaller than the
element — the wrong operation to apply before a nucleus detector.

---

## 4. Narrative that the saved results contradict

**4a. "The single 51 px template wins outright"** (cell 7).

> **This finding did not survive the fixes — the CSV I checked against was itself wrong.**
> Its AUCs came from `score_probe`, which capped the nucleus population at 20 000, left the
> annotated mitoses and look-alikes inside that population, and scored ties as losses
> (§2, §5). Recomputed correctly, `scale_1.0_only` wins AUC *and* median rank on **both**
> images — 0.953/0.879 and ranks 85/1103 — so the notebook's claim was better supported than
> this review concluded. What survives is the narrower point: the cell showed one image and
> never displayed the CSV, the margin on 246.tiff is 0.002 in AUC, and `discrimination`
> still prefers `fused_z_normalised` on 002.tiff. The notebook now regenerates and displays
> the whole table. Full numbers in the corrected-run log.

The original (now-superseded) observation was that the CSV held a second image the cell
never shows, and that on it the rejected variant won:

| 246.tiff | AUC mitosis | median mitosis rank |
|---|---:|---:|
| `scale_1.0_only` | 0.853 | 2933 |
| `fused_z_normalised` | **0.866** | **2672** |

The conclusion is drawn from n=1 (002.tiff), and the CSV that would show the
disagreement is referenced in a source comment but never loaded into the notebook.

**4b. The selection criterion is not the one the notebook argues is binding.** The
config was chosen on `auc_mitosis_vs_nucleus`. But cell 27's own point 2 — the central
claim of the whole experiment — is that mitosis-vs-look-alike discrimination is what
limits the method. On the `discrimination` column already present in that CSV,
`scale_0.6_only` wins on 002.tiff (0.110 vs 0.036) while `scale_1.0_only` wins on
246.tiff (0.061 vs −0.031). The two criteria disagree, and so do the two images.

**4c. "~2 minutes per image: ~100 s for the 72 correlations"** (cell 9). The run used
**24** augmentations — cell 8 prints it, and `n_augmentations = 24` in the metrics CSV.
`t_match_s` is 30–105 s (mean 54) and the whole run was 481 s over 7 images (69 s each,
baselines included). "72" is a leftover from the rejected multi-scale config.

**4d. "The minimum spacing between two MIDOG++ annotations is 27 px (245.tiff)"**
(cell 34, and the `find_and_suppress.py` docstring). Measured: the global minimum is
**26.25 px on 403.tiff**; 245.tiff is 26.57 px. More important, **within the seven ROIs
actually used the minimum spacing is 40.7 px** (301.tiff), comfortably above both the
25 px NMS radius and the 29.6–33.1 px match radius. The ambiguity cell 34 says it is
checking cannot arise in this run — which is why cell 35 returns 0 disagreements on
every row. As run it is a vacuous check, not a validation.

**4e. "The top 0.1 % is about the size of the whole top-K budget"** (cell 27). True for
002.tiff (0.1 % of 8 604 = 8.6 ≈ K=8). Not for 246.tiff (0.1 % of 20 000 = 20 vs
K=115). The sentence generalises from one image.

---

## 5. Bugs that do not bite yet

- **`optimal_assignment_disagreement`** [`evaluate.py:245`] `len(greedy ^ optimal) // 2`
  assumes every disagreement is a swap. A detection matched in one scheme and unmatched
  in the other contributes 1 to the symmetric difference and is floored away. It
  reported 0 everywhere only because the two sets were exactly equal (see 4d).
- **`info["max_score"]` is always the seed self-hit.** [`find_and_suppress.py:104`] It is
  recorded before self-hit removal, so `fs_metrics.csv` has `max_score` identical to
  `seed_self_score` (1.000000, or 0.999999 on 002) on every row. It should be the max
  after removal — the interesting number, which the threshold sweep shows ranges from
  0.70 (301.tiff) to 0.91 (405.tiff).
- **`threshold_metrics` is dead code.** [`evaluate.py:145`] `evaluate_run` is called with
  `score_threshold=None` everywhere in `experiment.py`, so it never runs — confirmed by
  the absence of `precision` / `fp_composition_lookalike` from the metrics CSV.
  `experiment.scale_usage` is likewise never called (and is meaningless at one scale).
- **`pick_seed` border test is looser than the reader's.** [`experiment.py:59`] It accepts
  `cx < w - border`, while `read_padded_patch` requires `round(cx) + half <= w - 1`. A
  seed at `cx = w - 36.4` passes the first and raises in the second. Never triggered here.
- **`score_probe`'s `median_mitosis_rank` is a mean, not a median.** `(1 − AUC) · n` is
  E[number of nuclei outscoring a *random* mitotic figure]. The name and both docstrings
  say median. Two smaller biases in the same function: the AUC uses strict `<` so ties
  count as losses, and `s_nuc` includes the mitotic figures themselves (negligible at a
  0.1–0.6 % base rate, but it should be documented).
- **`select_domain_images` uses `groupby(...).first()`** [`experiment.py:48`], which takes
  the first non-null value *per column independently*. Safe here only because nothing is
  NaN after the `fillna`. `.head(1)` is the intended semantics.

---

## 6. Two design points that shape the results and are not stated

**6a. The ROI selection is optimistic for this method.** Picking the densest ROI per
domain is well justified on statistical-power grounds (cell 4), but the direction of the
bias goes unmentioned: `recall@K` with `K = n_mitotic` gets easier as mitosis density
rises, because the base rate climbs against a roughly constant nucleus count. The
results broadly show that ordering — 246.tiff (K=115) 0.191 and 301.tiff (K=217) 0.097
at the top, 506.tiff (K=11) and 350.tiff (K=3) at zero — though it is not monotone
(002.tiff, K=8, scores 0.125, above 201.tiff's 0.059 at K=17). These seven numbers are upper-ish
bounds for their domains, not representative draws.

**6b. The seven seeds are one fixed quantile, not seven random draws.**
`run_experiment` constructs `np.random.default_rng(seed)` **inside** the per-image loop,
so `pick_seed` consumes the same first value every time. The chosen index lands at the
same relative position in every image's eligible list:

    n=218 → idx 185 (q=0.85)   n=116 → 98 (0.85)   n=18 → 15 (0.83)   n=14 → 11 (0.79)
    n=12  → 10  (0.83)         n=9   → 7  (0.78)   n=4  → 3  (0.75)

The notebook says "one randomly-chosen mitotic-figure seed each". It is the ~80th
percentile of each list, seven times. Harmless for a single pass, but the per-domain
spread carries even less independent information than `experiment.py`'s caveat implies —
and the planned multi-seed sweep must advance the generator across images, or it will
repeat the same quantile sequence.

---

## 7. Reproducibility

`import midog_utils` fails under the default `python3` (3.11.5): the installed `cv2` was
built against numpy 1.x and numpy is 2.3.4 (`AttributeError: _ARRAY_API not found`).
`tifffile`'s LZW path — needed for 505/506 — fails the same way through `imagecodecs`.
The notebook's outputs exist, so it ran somewhere, but nothing records that environment:
`requirements.txt` is the fastai training pin, which `midog_utils/__init__.py` explicitly
says the package does not need. A small `requirements-midog-utils.txt` (numpy, opencv,
scikit-image, scikit-learn, scipy, tifffile, pandas) would close it.

---

## 8. What is correct

Worth stating, because it is most of the code and it is what makes the rest credible:

- **The coordinate gate is the right test and passes for the right reason.**
  `matchTemplate` returns a top-left-indexed map; shifting each by `(size − 1) // 2` of
  *its own* size before fusing is correct, and asserting on the **fused** map rather than
  per-augmentation maps is exactly the check that catches a misalignment that
  per-template recovery would hide. Forcing odd sizes to avoid half-pixel rounding is right.
- **The patch-size bound holds.** For `BASE_SIZE = 51`, an arbitrary rotation of the
  inner 51×51 draws on source pixels out to 25·√2 = 35.36 px, and `PATCH_SIZE // 2 = 36`.
  No replicated border pixels enter a template. (The docstring's `ceil(51·√2) = 73`
  argument is a different and looser one than the bound that actually holds, but the
  conclusion is right.)
- **`max_peaks = 250000` vs "~173k theoretical max" is correct reasoning.** Under
  `fused >= dilate(fused, 15×15)`, two non-tied maxima must be ≥ 15 px apart in Chebyshev
  distance, so density is bounded by 39 Mpx / 225 ≈ 173k. Observed peaks are 17k–38k;
  the cap never binds.
- **The small-template bias argument is sound.** `TM_CCOEFF_NORMED` normalises by the
  template's pixel count, so a smaller template has a wider null and a higher maximum by
  chance; a raw element-wise max across sizes does favour it. `_robust_z` is a reasonable
  correction — see 4a for the part that should change.
- **Score-ordered greedy NMS**, and the diagnosis of why the reference's index-ordered
  `nms_with_area` is wrong.
- **The dataset conventions.** `bbox` is TLBR, not COCO xywh; every box is a synthetic
  50×50 around a point click, so point matching rather than IoU is right;
  `ResolutionUnit` really is 3 (cm) on 505/506 and 2 (inch) elsewhere, and ignoring it
  would give 0.578 µm/px there. `check_roi_scale` catching that at 2 mm² is a good guard.
- **Greedy-once matching for a monotone FROC** is the correct choice and the stated
  reason (optimal assignment recomputed per threshold can make sensitivity fall as the
  threshold drops) is right; it also matches `evalutils.score_detection`.
- **`unanimous` is well-defined**, though worth knowing it is not independent
  information: category-1 annotations carry only two label multisets, `(1,1)` ×8 917 and
  `(1,1,2)` ×3 020, so for category 1 `unanimous` is *exactly* equivalent to
  `n_votes == 2` — a clean proxy for "needed a third reader", but not a separate signal.
  (Category 2 is not so clean: `(2,2)` ×8 810, `(2,2,2)` ×3 044, `(1,2,2)` ×2 495.)
- **The headline is honestly reported.** At equal budget the blob detector beats the
  matcher on 3 of 7 ROIs (0.599 vs 0.097 on 301.tiff, 0.357 vs 0.191 on 246.tiff, 0.235
  vs 0.059 on 201.tiff), ties on 3, and loses on 1. Cell 15 says in advance that this
  would be the finding, and cell 16 prints it without softening.

---

## Priority

1. §1a — give `recall_by_agreement` a `k`; it is a one-line fix that restores a table.
2. §1 — move the operating point off the 0.25 floor (coverage is ~4 % at 0.60 on
   301.tiff), then reframe or drop the full-list tables and add the coverage floor.
3. §2 — raise or record the 20 000 blob cap; it currently falsifies a stated claim.
4. §3 — replace `tissue_mask` with a fixed brightness cut (e.g. `gray < 220`) or Otsu
   constrained to a bimodal check; global Otsu on solid tissue does not mean what the
   docstring says.
5. §4 — correct the four narrative claims; show `fs_fusion_variants.csv` in the notebook.
6. §6b — advance the RNG across images before the multi-seed sweep.


---

## Resolution status (2026-08-27)

| § | finding | fix | evidence |
|---|---|---|---|
| 1a | `recall_by_agreement` un-budgeted | `k` parameter, mirrors `lookalike_attraction_rate` | table now varies (0.104/0.074, 0.030/0.125) instead of all 1.0 |
| 1 | coverage saturation invisible | `evaluate.coverage_fraction`, reported per method | 0.876–0.948, printed beside every full-list number |
| 1b | NMS radius below match radius | derived per image from µm/px | verified: no surviving pair closer than the radius |
| 1c | `all_*` framed as a "natural operating point" | docstrings and cells 20/22 rewritten | — |
| 2 | blob cap at 20 000 | uncapped; `budget_delivered` recorded | 20 370 and 21 059 now delivered — the cap was binding |
| 3 | `tissue_mask` splitting tissue, not glass | fixed cut at gray 220, opening dropped | 0.78–1.00 fraction; **0 of 690 mitotic GT excluded across all 14 ROIs** |
| 4a | "wins outright" | **finding was wrong** — see the inline correction | table regenerated and displayed in the notebook |
| 4b | criterion mismatch | `discrimination` column shown; disagreement stated | holds: criteria differ on 002.tiff |
| 4c | "72 correlations, ~2 min" | corrected to 24 augmentations, 26–47 s | measured |
| 4d | "27 px (245.tiff)" | corrected to 26.2 px (403.tiff); vacuity of the check stated | measured |
| 4e | "top 0.1% ≈ K" | stated per image; **and the underlying p90/p99.9 claim reversed** on 002.tiff once the nucleus set was cleaned | see corrected-run log |
| 5 | `// 2`, `max_score`, dead `threshold_metrics`, `pick_seed` border, mean-vs-median, `groupby().first()` | all fixed | unit-tested individually |
| 6a | ROI selection optimistic | documented in `select_domain_images` | — |
| 6b | seeds one fixed quantile | per-image RNG stream | quantiles now 0.21–0.89, previously 0.75–0.85 |
| 7 | environment unrecorded | `requirements-midog-utils.txt` | — |
| — | `viz.overlay` circle radius was the constant 4 | drawn at the match radius | not in the original plan; flagged when applied |

**Two conclusions reversed under the corrected code** — §4a (the fusion comparison) and the
`p99.9_nucleus > p90_mitotic` claim behind §4e. Both had the same root cause: a contaminated,
truncated nucleus population feeding the response-map statistics. That is the strongest
argument the exercise produced for regenerating stored comparison tables from fixed code
rather than citing them.
