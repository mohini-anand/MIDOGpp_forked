# Future ideas and directions to explore

**Started:** 2026-10-02. This is a living backlog, not a plan. Each item has:
- a status;
- the evidence that motivates it;
- the test that would decide it.

When an item is run, record the result here in one line and link to the notebook or log.

**Status labels:** `idea` · `spec'd` · `running` · `tested: helps` · `tested: no effect` · `tested: hurts` · `adopted` ·
`dropped`.

**Related documents:**
- `Mitosis annotation co-pilot pipeline improvement notes.pdf` (2026-09-22). On 245.tiff alone, it covers:
  - its "Stage 5: the feedback loop";
  - "Next experiments".
- `low_precision_rois/`, the 2026-10-02 analysis of why ROIs fail at P@10/P@20: theory, citations and montages.
- `CANDIDATE_REFINEMENT_RERANK_PLAN.md` and `CANDIDATE_FEATURE_EXTRACTION_PROMPT.md`.

---

## 0. The constraint everything below works against

Both ranking axes rank by darkness:
- `TM_CCOEFF` is an unnormalised zero-mean linear filter, so it scores centre-minus-surround hematoxylin OD.
- `chromatin_od` is darkness by construction.

The ROIs that fail have two properties (`low_precision_rois/`):
1. **A dark-object tail.** Dozens to hundreds of non-mitotic objects are darker than the ROI's few best mitoses. These
   are pigment, small dark round cells, apoptotic-looking nuclei, or crowded hyperchromatic nuclei.
2. **Few mitoses.** 14–27 in most of them.

Untouched MIDOG++ ROIs are sparser than the 49 used so far. Median mitoses per ROI:

| tumour type | median mitoses |
|---|---|
| mast cell | 10 |
| lung | 14 |
| lymphoma | 57 |
| sarcoma | 8 |
| breast | 6 |
| melanoma | 8.5 |
| neuroendocrine | 4 |

So in the product, ROIs like the failing nine will be the norm, not the exception.

---

## 1. Using human input across rounds: the priority

The loop: the pathologist draws a box, the tool shows K candidates, the pathologist accepts or rejects each one, then
draws another box. A learned detector cannot use that input; this tool can.

### 1.1 What feedback can and cannot change

**The ceiling cannot be raised inside a fixed region.** Page t has a ceiling P@K ≤ min(1, R_t / K), where R_t is the
number of mitoses not yet found in the searched region. Every accepted mitosis lowers R_t.

**What feedback *can* do:**
- **Bring precision closer to the ceiling.** This is the main lever.
- **Raise R_t, by widening the searched region.** The product shows the whole slide, so neighbouring tissue is the only
  way to find more mitoses once an ROI is depleted.
- **Tell the user when to stop.**

**Multi-round metrics to adopt:**
- per-page *ceiling-normalised* precision, P@K / min(1, R_t / K);
- cumulative recall against candidates reviewed;
- time to 50 / 80% recall under the 2026-09-30 cost model. A page pays when E = pK·T_scan / (r·(K·v + c_box + L)) > 1.
- Report the worst ROI, not only the mean.

**Prior evidence**

- **Max-aggregating confirmed templates** (the Sep-22 notes, Stage 5, on 245.tiff). Adding each confirmed figure as a
  template, with the maximum score as the aggregate, found more mitoses but did not raise precision: "taking the
  maximum keeps every template's own false positives".
  - Seed 6243 found 31 of 89 mitoses in 150 reviews, at 0–4 hits per round with no upward trend.
- **The 2026-09-30 loop simulation** (`deep_pool_mo_pad61`, list exclusion). Median per-page precision for ROIs with
  20 or fewer mitoses was .34 .18 .13 .09 .06 .04.
  - Per-ROI AUC re-weighting of the old features added +0.3 [−0.06, +0.72] mitoses over 6 pages. The features were too
    weak to re-weight.
- **FP families repeat across clicks** (2026-10-02, A3 top 10). Across 6 clicks, a failing ROI's top-10 false positives
  are only 15–32 distinct objects, each shown about 2.25 times.
  - Excluding already-shown objects removes the exact repeats. It does not reach the rest of the family: the next
    pigment clump, the next small dark round cell.
  - So feedback has to generalise from a rejected object to objects *like* it.

### 1.2 What each round gives you

| input | content | value |
|---|---|---|
| drawn box | a positive, plus its extent (size and scale) | a new template, or a better class mean |
| accepted candidates | positives found by the tool | the class mean; templates for other phases |
| rejected candidates | **hard negatives specific to this ROI** | under-used; failing ROIs' false positives come in a few families. Whether negatives are worth as much as positives here is open: in text retrieval "positive feedback also turns out to be much more valuable than negative feedback" [Manning 2008, ch. 9] |
| the rest of the ROI | 6k–35k nuclei, ≤ 2% mitoses (0.07–1.9% across the 49) | a free estimate of the background (clutter) distribution, with no labels needed |

### 1.3 Ideas

**H1. A closed-form discriminant template, updated every round.** `idea`
- **The formula.** Replace the raw template with w_t = (Σ + λI)⁻¹ (μ_pos − μ_neg):
  - Σ is the covariance of patches (or feature vectors) over the ROI's whole nucleus population, computed once at
    slide load;
  - μ_pos is the mean of the box plus the accepted positives, after alignment;
  - μ_neg is the background mean, shrunk toward the mean of the rejected candidates.
- **Theory.** This is the Hotelling observer, the optimal *linear* detector in correlated clutter [Li 2021, in
  `low_precision_rois/`]. The current template is the non-prewhitening special case, which is only optimal in white
  noise.
  - With one exemplar and a generic background covariance, this is "exemplar-LDA" [Hariharan 2012]:
    - "w = Σ−1(xmedoid −µ0)";
    - "we can learn a generic µ0 and Σ from unlabeled windows";
    - "Background statistics in the form of Σ and µ are computed just once, and training only involves computing the
      mean of the positive examples".
  - What Hariharan 2012 found:
    - Whitening is what makes it work. Their un-whitened "centered model", the closest analogue of plain template
      matching, "achieves an AP of less than 10%, indicating that whitening is crucial".
    - Exemplar-LDA showed "little or no drop in performance" against exemplar-SVMs (with a calibration step).
    - A **single** exemplar is the weak case: mean AP 4.84 for LDA on one exemplar, against 7.59 for LDA on a cluster.
      That fits the loop. Round 1 has one exemplar; every accepted mitosis moves μ_pos toward the stronger cluster
      case.
    - Σ needed regularising ("adding a small value (λ = .01) to its diagonal").
- **Why it should help the failing ROIs.** Clutter that is mostly dark blobs makes Σ largest along the mass direction.
  Σ⁻¹ down-weights that direction and moves weight onto outline and texture. That is the derivation in
  `low_precision_rois/` §3.4.
- **Cost:**
  - Σ is computed once per ROI at load;
  - each round is one solve of size d × d, with d ≈ 50–200 after PCA: milliseconds;
  - applying w is still a single linear filter (`matchTemplate` or an FFT), or a dot product over precomputed
    per-nucleus features.
- **Test.** Simulate rounds on the 294 recorded clicks, with ground truth as the reviewer. Compare ceiling-normalised
  per-page precision against the current template.
- **Risks:**
  - the choice of space (down-sampled H-OD pixels or a feature vector) and λ;
  - the alignment of positives, because rotation matters;
  - the estimation domain. Hariharan estimated Σ from generic windows at every offset, so his filter applies
    everywhere. If Σ comes from nucleus-centred patches, sliding w over every pixel applies it where Σ was never
    estimated. Either estimate Σ from random tissue windows and slide the filter, or keep the nucleus-centred Σ and score
    only at nucleus centres (the G7 path).

**H2. Rocchio-style relevance feedback in feature space.** `idea`
- Each round, q_{t+1} = α·q_t + β·mean(accepted) − γ·mean(rejected), applied over per-nucleus feature vectors
  precomputed at load [Manning 2008, ch. 9, see the source list below].
- It is H1 without the covariance: simpler, and weaker when the clutter is correlated.
- It is the right baseline for H1.
- Manning 2008, ch. 9:
  - typical weights are "α = 1, β = 0.75, and γ = 0.15";
  - relevance feedback "has been shown to be most useful for increasing recall";
  - it needs an initial query "at least somewhere close to the documents they desire";
  - it "does not work as well if the relevant documents are a multimodal class".
- Mitoses are multimodal: plate, ring, anaphase pair, prophase clump. So keep one prototype per phase family, not one
  mean (see H3).

**H3. Aggregate the template set with mean or top-k mean, not max.** `idea` (Sep-22 notes, items 1–2)
- Max acts like OR and favours recall; a mean acts like AND and favours precision.
- Weight each template by its own confirmation rate.
- Possibly use max on the first page after a new box and mean on the pages after.
- Average *within* a phase family, then take the max *across* families. Averaging a plate with an anaphase pair blurs
  both.

**H4. Suppress negative families.** `idea`
- Penalise each candidate by its similarity to the rejected candidates:
  score − γ·max_j sim(candidate, rejected_j).
- Similarity can be measured in feature space or by template correlation.
- In `211`, one rejected pigment clump should push every pigment clump down.
- The heavier version is one linear model per exemplar: "a separate linear SVM classifier for every exemplar … defined
  by a single positive instance and millions of negatives" [Malisiewicz 2011]. Here it would be applied to each rejected
  family.
- Test it against Manning's warning that negative feedback is usually worth less than positive feedback.
- **Test:** the same simulation as H1, plus the share of each FP family remaining by page 3.

**H5. Re-weight features per ROI from the first page's labels.** `tested: no effect` (+0.3, 2026-09-30) with weak
features.
- Revisit once the outline and texture features in §3 exist and are shown to be valid.

**H6. Spend some slots on learning, not only on yield.** `idea`
- Reserve 1–2 of the K slots for the candidates whose labels would change the model most. Examples: the top member of
  an FP family nobody has reviewed yet, or the most uncertain candidate under H1.
- This gives up precision now for precision later. Decide the split with the cost model.

**H7. Blank shown candidates in the channel, instead of excluding them from the list.** `idea`, untested (2026-09-30)
- An A/B costs about 882 re-searches, roughly 75 min.
- Watch for leak peaks at the blank edges.

**H8. Widen the search when an ROI is depleted.** `idea`
- This is the only lever on R_t.
- It needs tissue beyond the 2 mm² MIDOG++ ROIs to test: a whole slide, or neighbouring regions.

**H9. A stopping rule.** `idea`
- The candidates:
  - two consecutive pages with at most 1 hit (Sep-22 notes);
  - the model's predicted page yield falling below break-even (E < 1).
- Report the recall reached at the stop.
- **Adaptive page size.** The ceiling on page t is R_t / K. Shrinking K as the predicted yield falls is the one way to
  keep per-page precision up while the ROI depletes. The cost model decides the trade.

**H10. Speed architecture.** `idea` (consistent with `Research Logs/2026-09-01-one-click-retrieval-literature.md` §2)
- **At slide load, everything that doesn't depend on the click:**
  - the channels;
  - the nucleus list (`nucleus_blobs`);
  - the per-nucleus feature vectors;
  - Σ and its PCA basis.
- **Per round:**
  - update w (H1) or q (H2);
  - re-score every nucleus with one matrix-vector product (about 15k × d: milliseconds);
  - drop what was already shown;
  - show the top K.
- No re-search is needed per round.

**Evaluation harness for §1.**
- A round simulator on the 49 ROIs, with GT as the reviewer:
  - the 294 recorded clicks;
  - several click orders.
- Report the multi-round metrics in §1.1.
- Confirm on held-out sparse ROIs before adopting anything.
  - There are 13–44 unused valid ROIs per domain, and none in neuroendocrine (counted 2026-09-30).
  - Pull them using the pattern in `human_neuroendocrine_followup/pull_rois.py`. `download_testing_set.py` no longer
    runs.

---

## 2. More mitoses into the candidate pool (the generator)

**Target.** In the five worst ROIs, the top 30 / 100 / 300 by score hold 2.2 / 4.2 / 7.1 mitoses, against 13.7 /
27.1 / 39.3 elsewhere.
- Judge every item here by mitoses in the top 30 / 100 / 300 on the nine failing ROIs.
- Require no regression on the other 40.

**G1. Whiten the click's template.** `idea`
- This is H1 at round 1. Σ comes from the ROI's own patches; mitoses are at most 2% of its nuclei, so no labels are
  needed. Mind the estimation-domain point under H1's risks.

**G2. Partial whitening: a high-pass or Laplacian prefilter on H-OD before correlation.** `idea`
- Lewis (1995): "With typical (ρ ≈ 0.95) image correlation the best pre-filtering is approximately Laplacian rather than
  a pure whitening".
- **Caution:** fully deleting the mass term cost 12 pts pooled (`template-is-mostly-a-mass-detector`), because darkness
  works where mitoses *are* the darkest objects.
- G1 sets the weight per ROI; G2 sets it once. Test both.

**G3. Match on several channels.** `idea`
- Correlate the H-OD template together with templates on:
  - a high-frequency channel: white top-hat, small-scale LoG (~0.5–1 µm), or gradient magnitude;
  - the E channel (cytoplasm colour).
- Sum the per-channel scores, or whiten jointly (G1).
- The hair-like projections and the absence of a membrane live in the high-frequency channel.

**G4. A relative-darkness channel for crowded, hyperchromatic fields** (`524`, `527`, `245`). `idea`
- Divide H-OD by a local nuclear-OD baseline, for example the median OD of nuclei within about 50 µm. A mitosis that
  is 20% darker than its neighbours then scores like one in a pale field.
- First check it against D3. D3's rationale is clipping by global min-max rescaling; a local *ratio* is a different
  operation, but should be argued explicitly.

**G5. A pigment gate in the search channel** (§4). `idea`

**G6. A template bank over rotation, flip and scale.** `tested: mixed`.
- Rotations gave +4–6 pts at P@10 on 14 ROIs pre-D11 by `tm_score`, but were never adopted.
- The gain disappeared under chromatin ranking (p 0.33).
- Re-test with a mean or top-k aggregate (H3), scored by pool content.

**G7. Score every nucleus, not just TM peaks.** `idea`
- `nucleus_blobs` puts a nucleus within 7.5 µm of 77–100% of GT mitoses (90–100% in the failing nine).
- With features precomputed at load, the "generator" becomes a ranking over all nuclei. That removes peak-finding and
  NMS as places where mitoses are lost.

**G8. Frozen pretrained embeddings with LDA (H1), not cosine.** `idea`, last resort before training anything.
- The 2026-09-01 test used cosine-to-click (`lunit_dino_vits16`, `kaiko_vits16`). That is the embedding-space analogue
  of the non-prewhitening template, and it read 3–42× worse than the seedless blob detector.
- A whitened version against the ROI's own background was never tried.
- The commercial-licence constraint applies.

---

## 3. Features orthogonal to darkness: hair-like projections, no nuclear membrane

**Where these are measured.** At the candidate stage (top N, or every nucleus at load), not inside `matchTemplate`.
- A linear filter cannot compute a mask-based boundary statistic.
- It can carry part of the signal through a high-frequency channel (G3) or through whitening (G1).

### Already specified in `CANDIDATE_FEATURE_EXTRACTION_PROMPT.md`

- **Outline**, on the undilated `M_chr`, with a synthetic star / spiked-ellipse suite:
  - `chr_roughness`;
  - `chr_concavities_per_um`;
  - `chr_radial_hf`.
- **Membrane and clear zone**, with a blur test:
  - `rim_core_nuc`;
  - `ring_peak_prominence`;
  - `ring1_frac_nuclear`.
- **Texture**, with a shuffle test:
  - `h_chr_grad`;
  - `h_chr_lap`;
  - GLCM contrast and homogeneity.
- **Colour:** brown-in-surround features.

### Gaps worth adding before that spec runs

`spec'd` here, not yet in the prompt.

1. **Depth.** Its candidate set is the top 100. The five worst ROIs have only 4.2 mitoses there per click, so validity
   on exactly the ROIs that need these features rests on about 4 positives.
   - Extend to the top 300 (it is in `deep_pool_..._top1000.csv`), or add every GT mitosis as an extra row.
2. **Black pigment sits *inside* the chromatin mask.** The spec assumes pigment lives outside `M_chr`. That is true
   for brown melanin, which projects mostly to the D channel, but not for neutral black pigment.
   - Black pigment projects about 70% of its OD onto H, so it *is* the chromatin mask in `211`.
   - Add an OD-direction feature inside `M_chr`: the angle from the ROI's hematoxylin vector, or the blue/red OD ratio,
     read on non-saturated pixels.
   - Add `211`, `410`, `291`, `289` and `339` to the pigment checks. The listed melanoma ROIs (`524`, `527`, `546`)
     showed no pigment among their top-10 false positives.
3. **Montages of the failing ROIs** `211`, `242`, `546`, `524` and `245`, in addition to `300`, `201`, `013`, `529` and
   `425`.
4. **AUC conditional on darkness.** Compute each feature's within-click AUC inside darkness strata (deciles of
   `h_chr_mean`). The question that matters is whether a feature adds information *beyond* darkness; `rho_with_h_chr_mean`
   only says whether it duplicates darkness.
5. **Resolvability by scanner, and focus.**
   - Report boundary-feature AUC per scanner (3D Histech, Hamamatsu XR, S360, Aperio). The Sep-22 notes record that
     "about half of mitoses on whole-slide images lack clear hairy projections" and that "at 0.25 µm/px, spikes often
     appear only as a fuzzy edge" (from the sources those notes cite; not re-verified).
   - Add a focus measure, such as the variance of the Laplacian in the window. Blurred dark objects are a named
     false-positive class.
6. **A mask-free protrusion measure.** Use a white top-hat of H-OD with a disc of about 0.5–0.75 µm: the share of the
   object's OD carried by structures thinner than the disc. It needs no segmentation, so it cannot inherit mask
   failures (cf. F4/F6).
7. **Size relative to the neighbours.** Chromatin area divided by the median nucleus area within 101 px. Mitotic cells
   are the same size as their neighbours or larger; apoptotic bodies are smaller (Sep-22 notes).

---

## 4. Colour and pigment

**What the colour model does with pigment.**
- `rgb2hed` projects a neutral black absorber onto H at about 70% of its OD (computed 2026-10-02). A brown, blue-heavy
  absorber projects mostly onto D.
- So the H channel is contaminated by black or dark-brown pigment (`211`; incidental in `289`, `291`, `339` and `410`)
  much more than by melanin.

**A fixed colour rule will not transfer across scanners.** Measured on the darkest pixels of GT mitoses, the median
blue/red OD ratio per ROI is:

| scanner | range of ROI medians |
|---|---|
| Hamamatsu XR | 0.38–0.58 |
| 3D Histech | 0.45–0.88 |
| Aperio CS2 | 0.69–0.92 |
| Hamamatsu S360 | 0.80–0.84 |

**Saturation hides colour.**
- 64% of GT mitoses have at least one near-saturated pixel (mean RGB < 40) among their darkest pixels.
- So do 100% of `211`'s top-10 false positives.
- The earlier F9 audit found that a D/H ratio penalised saturated *mitoses* on `425`, `460` and `527`.
- Read colour on the mid-tone pixels.

**Options, in order:**
1. **A residual gate.** Estimate the H and E vectors per slide with Macenko's method [Macenko 2009]:
   - "find these particular stain vectors for each image based on the colors that are present";
   - the steps are an SVD of the OD values (dropping OD < 0.15), the plane of the two largest directions, and the
     robust extreme angles (1st and 99th percentiles).

   Two practical points:
   - The paper's algorithm returns two *unlabelled* vectors, so add a rule that assigns H, such as the one closer to
     Ruifrok's H.
   - The third vector is the residual, by cross product. Landini 2021 credits that step to their ImageJ plugin; the
     Ruifrok manuscript does not state it.

   Pixels whose OD is poorly explained by H + E (a large residual on that third axis) are clipped to the floor b in the
   search channel, with the same mechanics as the A2 clip.

   **Saturation guard.** Saturated pixels are near-neutral, so their residual is large: the F9 failure. Without a
   guard, this gate would hit the 64% of mitoses that have saturated pixels. Either decide "pigment" per object from its
   mid-tone pixels only, or exempt pixels with gray < 40 from the gate.

   Vahadane 2016 (sparse, non-negative stain separation) is the alternative. Only its abstract was read.
2. **Per-slide stain estimation as the general upgrade.** It stabilises every threshold that reads H across scanners.
   It is already in the per-iteration plan's backlog ("stain vectors").
3. **A per-ROI pigment vector as the third stain**, only if 1 and 2 are not enough.

**The test.** Measure how many of `211`'s top-10 false positives are removed (and those of `289`, `291`, `339` and
`410`), and how many mitoses are lost across all 49 ROIs. The per-iteration plan's gate is that the pigment veto must
reject 0 mitoses.

### 4.1 Floor-clip stress test across stains and scanners

**Status:** `spec'd`, **deferred until before shipping** (decided 2026-10-02).
- It does not block the per-feature analysis. The window floor ties `b_otsu` (P@10 +0.24, pool +0.04 mitoses per
  click), so the candidate pools barely depend on which floor ships.
- The clipped A3 pools are saved: `results/tm_floor_definition_49roi_6seed_pool.csv.gz`, 294 clicks, up to 100 each,
  with buckets.
- **Sources:** the plan an agent proposed after `tm_floor_definition_49roi_6seed.ipynb`, with the review corrections
  below. The descriptive checks are in `floor_robustness/`.

**Why.**
- The floor's effect is asymmetric:
  - at 0.5× `b_otsu`, about 40% of the P@10 gain is lost;
  - from 1.0× to 2.0× it is flat in the mean, but clicks worse than no clip rise from 20 to 58.
- The floor is estimated from each image, so on another lab's slides it lands wherever that slide's stain and scanner
  put it.
- **The question:** does the floor estimator stay at or above the knee when stain and scanner change?

**What theory predicts (checked in code, 2026-10-02).**
- The pipeline is invariant to H → a·H + c with a > 0:
  - the `TM_CCOEFF` score scales by a²;
  - the peak threshold is median + z·MAD of the score map;
  - the floor and the box refinement are multi-Otsu, with geometry in µm;
  - the template fill is the window's own mean;
  - the border band uses a −3e38 sentinel.
- So these perturbations should change nothing:
  - both stains scaled together (section thickness);
  - RGB gain and white balance (a constant added to each OD channel);
  - pure gamma (a scaling of OD).
- **What breaks the invariance:**
  - the `rgb2hed` clamp of negative values at 0 (3.9% of pixels in the median ROI, up to 25.6%);
  - RGB saturation and 8-bit rounding;
  - every change that is not an affine change of the measured H.

**Stage A: a screen of floors only, no searches.** Run first, minutes of compute.
- **Estimators:**
  - `b_otsu`;
  - `W128` (a window around the click, glass included);
  - a tissue-only window or tile floor;
  - TOTSU (tissue-only Otsu over the ROI).

  Optionally add a tissue quantile and an absolute floor as reference points. The product has no 2 mm² ROI, so the
  window or tile estimators are the ones that could ship.
- **Perturbations, at one realistic strength.** Only those the theory does not cover:
  - hematoxylin vs eosin ratio change (H ×0.6 and ×1.4, E ×1.5), synthesised in each slide's own stain basis:
    OD_new = OD + (k_H − 1)·c_H·h_slide + (k_E − 1)·c_E·e_slide, using the per-ROI Macenko vectors and keeping the
    residual;
  - black-point or contrast curves in RGB;
  - Gaussian blur (σ ≈ 1 px);
  - JPEG re-compression;
  - one stain gradient across the ROI (section thickness);
  - windows shifted onto glass, lumen or the ROI edge.

  Include one exact scaling (H × 2) as a null control.
- **Score.**
  - The equivalent level: the multiple of the original `b_otsu` that clips the same share of *tissue* pixels.
  - The overlap with the original clipped pixel set. This is the primary score for changes that alter pixel ranks.
  - Whether the refinement box changes, or refuses the click.
- **Trigger for stage B, written before stage A runs** (proposed): any estimator lands with its median below 0.75×,
  or any ROI below 0.5×, under a realistic perturbation. If every estimator stays at or above about 0.9×, stop. The
  existing ladder already shows that 1.0× to 2.0× is flat in the mean.

**Stage B: search, only if triggered.** About 1–1.5 h on 6 workers. Check RAM per worker in the smoke run.
- **Arms:**
  - no clip;
  - the floor at 1.0×;
  - the floor at 1.5×, the undershoot detector: if it significantly beats 1.0× under a perturbation, the floor landed
    below the knee there;
  - the product-feasible window or tile floor, under *every* selected perturbation.
- **Perturbations:** chosen from stage A by a rule fixed beforehand: lowest and highest landings, one gradient, blur,
  JPEG, and one null control.
- **Report three numbers per perturbation:**
  - clip vs no clip under the perturbation;
  - the clipped pipeline, perturbed vs original;
  - no clip, perturbed vs original.
- **Also report per-click harm:** clicks worse than no clip, and the worst click.
- **Decision.** The clip carries over if it still beats no clip with the CI clear of 0, and its benefit shrinks by no
  more than about a third of the current gain (1.1 pts at P@10, 0.9 at P@20, 0.5 pool mitoses per click).
  - This applies at realistic strength only; strong settings are descriptive.
  - Add a per-click guard. Example: clicks worse than no clip at most about 1.5× the unperturbed count; the user sets
    the number.
  - Correct losses for multiple testing. Null controls are harness gates.
- **Harness gate for the new runner:**
  - exact parity against the recorded no-clip and 1.0× runs;
  - H × 2, c = 0 must give bit-identical lists, since scaling by a power of two is exact in floating point;
  - a general a and c may differ only at tie or rank boundaries.

**Stage C: real cross-scanner data, later.**
- **First, the free version.** MIDOG++ breast spans three scanners with tumour type held fixed (Hamamatsu XR, S360,
  Aperio CS2). The held-out supply had 28 unused valid breast ROIs; check how they split by scanner. These are
  different cases, not re-scans, so this check cannot certify robustness.
- Re-scanned same-slide data only if that shows trouble, after checking which scanners it covers, its labels, its
  licence, and any case overlap with MIDOG++.

**Dropped from the proposed plan, and why.**
- Phase 0b (stain spread and the regression of clip gain): already done in `floor_robustness/`, and null (leakage
  ρ 0.03, nothing above |ρ| 0.26). The per-ROI Macenko spread also mixes tissue with stain.
- Phase 0's tests through RGB (gain, white balance, gamma): predicted null from the code; one exact gate is enough.
- Two strengths in the screen: one realistic strength.

**Pitfalls the review found.**
- **A no-op stain test.** Stain changes made through `rgb2hed`'s own fixed matrix test nothing: E there is invisible to
  H by construction, and H there is a pure scaling.
- **The clamp reaches the clipped pipeline.** It isn't immune: `b_otsu` and the box refinement are computed on the
  unclipped, clamped channel.
- **Fixed RGB constants.** The white-adjacent FP counts and TOTSU's glass mask use gray > 215. Compute them from the
  *original* image, or the metrics will move even when the lists are identical.
- **The mean cannot settle 1.0× vs 1.5×.** That trade-off is per-click harm.

**Out of scope:** a different scanner resolution (a whole-pipeline question), new tissue (the held-out question), later
rounds of the loop, and chromatin od.

---

## 5. Evaluation changes

- **Report per ROI.** Put the mitosis count, the ceiling min(1, (M−1)/K) and the tail count beside every P@K, and
  report ceiling-normalised precision.
- **Judge generator changes by pool content.** Count mitoses in the top 30 / 100 / 300 on the nine failing ROIs:
  `211`, `242`, `524`, `546`, `462`, `527`, `410`, `245`, `289`.
- **Judge loop changes** with the multi-round metrics in §1.1.
- **Validate on sparse held-out ROIs** before adopting anything. That is the population the product will see.

---

## Sources added in this file

Sources already cited in `low_precision_rois/2026-10-02-…md` are not repeated here.

All were read in session on 2026-10-02 unless marked otherwise.

- **[Hariharan 2012]** Hariharan B, Malik J, Ramanan D. Discriminative decorrelation for clustering and classification.
  ECCV 2012, LNCS 7575, pp. 459–472. doi:10.1007/978-3-642-33765-9_33. Read as the archived author PDF.
- **[Malisiewicz 2011]** Malisiewicz T, Gupta A, Efros AA. Ensemble of exemplar-SVMs for object detection and beyond.
  ICCV 2011, pp. 89–96. doi:10.1109/ICCV.2011.6126229.
- **[Manning 2008, ch. 9]** Manning CD, Raghavan P, Schütze H. *Introduction to Information Retrieval*, ch. 9,
  "Relevance feedback and query expansion" (eq. 9.3, Rocchio). Cambridge University Press (2008).
- **[Macenko 2009]** Macenko M, Niethammer M, Marron JS, et al. A method for normalizing histology slides for
  quantitative analysis. ISBI 2009, pp. 1107–1110. doi:10.1109/ISBI.2009.5193250.
- **[Vahadane 2016]** Vahadane A, et al. Structure-preserving color normalization and sparse stain separation for
  histological images. *IEEE TMI* 35:1962–1971 (2016). doi:10.1109/TMI.2016.2529665. **Abstract only** (paywalled).
- **[Landini 2021]** Landini G, Martinelli G, Piccinini F. Colour deconvolution: stain unmixing in histological imaging.
  *Bioinformatics* 37:1485–1487 (2021). doi:10.1093/bioinformatics/btaa847.
- **[Lewis 1995]** Lewis JP. Fast Normalized Cross-Correlation, the expanded Industrial Light & Magic version.
