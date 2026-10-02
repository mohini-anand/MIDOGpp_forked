# Why some ROIs still fail at P@10/P@20: the measured mechanism, the theory behind it, and how good this can get

Date: 2026-10-02.

**Question asked:** *"Look at my most recent 3-5 experiments. Identify the domains/ROIs that still have very low TP hit
rate within the top 10/20 … tell me why … ground your response in theory and logic … Do a thorough literature search …
how good can we get … Do not make anything up."*

**Status:** post hoc and descriptive. Every ROI-level number here was computed on the same 49 ROIs that the experiments
were run on, after the failing ROIs were known. Nothing is committed.

## Inputs

**Experiments read**, newest first:

| Experiment | Notebook | Arm used here |
|---|---|---|
| TM border + floor clip, 2026-10-01 | `production_hematoxylin_only/multiotsu_pad_arms/tm_border_floor_49roi_6seed.ipynb` | **A3** (clip + band), the adopted arm |
| Candidate-refinement re-rank selection, 2026-09-30 | `candrefine_rerank_mo_pad61_select.ipynb` | **C1 = REF**, an in-sample selection |
| Deep pool, uncapped `mo_pad61`, plus the solidity re-rank, 2026-09-30 | `deep_pool_mo_pad61_49roi_6seed.ipynb`, `shape_mo_pad61_deep_pool_49roi.ipynb` | — |
| Multi-Otsu pad arms and their domain slice, 2026-09-30 | — | — |

**New measurements made for this log**, all in `low_precision_rois/`:
- `roi_mechanism_measure.py` runs on every one of the 49 ROIs. It uses the click-free `baselines.nucleus_blobs` nucleus
  population, 6.2k–35k objects per ROI, and every GT mitosis and look-alike. At each object it computes:
  - `od51` (= `chromatin.chromatin_density`) and `od21`;
  - the colour of the darkest pixels.
- It does the same for every A3 top-10 detection.
- Outputs: `per_roi_mechanism.csv`, plus crop montages (`montage_*.jpg`) of GT mitoses, GT look-alikes and the A3 top 10
  for clicks 0 and 1 on both axes.
- **Reproducibility gaps.**
  - The per-point and per-blob outputs of `roi_mechanism_measure.py` (about 50 MB) were left in the session scratchpad.
    Re-run it with an output directory to regenerate them.
  - `roi_mechanism_analyse.py` reads its per-ROI precision input from the scratchpad path. That input is saved here as
    `per_roi_precision_inputs.csv`; point `pr = ...` at it.
  - The tail count, the partial correlations, the domain permutation test and the FP-repeat count ran as inline
    snippets, not saved scripts. Their results are in `per_roi_mechanism.csv` and in this log.

**Literature:** two independent agents fetched and read every source below in this session. Sources that could only be
read as an abstract are marked **[abstract]**. Anything that could not be verified is left out.

---

## TL;DR

1. **The list (ROIs, not domains).** Nine ROIs are below P@10 = 0.35 on both A3 axes:
   - `211`, `242` (canine lung);
   - `524`, `546`, `527` (melanoma);
   - `462`, `410` (soft-tissue sarcoma);
   - `245`, `289` (lymphoma).

   Five more sit at the edge: `509`, `201`, `339`, `404` and `364`.

   Domain levels are still not distinguishable from ROI noise (label-permutation p = 0.07–0.27 on A3 and C1).
2. **One click-free number explains most of the per-ROI spread.** It is the count of *unannotated* objects in the whole
   ROI that are darker (`od51`) than the ROI's 5th-darkest mitosis.
   - It is a median of **54** in the nine failing ROIs and **1** in the other 40.
   - Spearman ρ with A3 P@10 at `tm_score` is −0.75 [−0.87, −0.56]; holding the mitosis count fixed, the partial ρ is
     −0.70.
   - The `chromatin_od` values (−0.86, partial −0.83) are close to definitional, because that axis ranks by the same
     `od51` statistic. The `tm_score` numbers are the informative ones.
   - Mitosis count adds independently (partial ρ +0.60 at `tm_score`).
   - So there are two drivers, **a dark-object tail and low prevalence**, and the failing ROIs have one or both.
3. **Which ROIs fail is predicted without the click.** Rank every nucleus in the ROI by darkness alone, with no click and no
   template. That ranking's P@10 correlates with the pipeline's `tm_score` P@10 at ρ 0.79, and with C1 at 0.80 (0.89
   at `chromatin_od`, which reads the same statistic).
   - This is what theory predicts for `TM_CCOEFF`. It is an unnormalised zero-mean linear filter, so its output is a
     centre-minus-surround optical-density contrast (§3.1).
   - That matches the repo's finding that the template is mostly a mass detector.
4. **What the dark objects are differs by ROI.** These are by-eye readings of the montages and still need a pathologist:
   - black granular pigment (`211`);
   - small, smooth, very dark round nuclei beside alveolar or lumen white space (`242`);
   - smooth dark nuclei with perinuclear halos (`546`, `462`);
   - crowded, overlapping hyperchromatic tumour nuclei (`524`, `527`, `245`).

   Each class matches an imposter type the mitosis literature names (pigment, inflammatory cells, apoptotic cells,
   hyperchromatic nuclei, cells in close proximity). Pigment is outside the stain model
   altogether (§3.5).
5. **How good can it get.**
   - In the five worst ROIs the 300 highest-scoring candidates hold 7.1 mitoses on average, so even a perfect re-rank of
     the top 300 averages P@10 ≤ 0.71.
   - In ROIs with about 15 mitoses, P@20 ≤ 0.70 and P@30 ≤ 0.47 by arithmetic.
   - Trained state-of-the-art detectors lose about half their precision in imposter-rich regions: 0.805 → 0.400 in
     MIDOG 2025.
   - P@10 ≥ 0.8 on *every* ROI is not reachable by a darkness/mass ranker on these ROIs. The information that separates
     these mitoses from their imposters is orthogonal to darkness. The literature names it (protrusions or "hairy
     extensions", no nuclear membrane), and the repo has not yet measured it in a form that could carry it.

---

## 1. Which ROIs, by a stated rule

**Rule:** P@10 < 0.35 on *both* A3 axes. A3 is the adopted production-like arm. P@10 is the mean over 6 clicks.

- C1 is shown beside A3. It is an **in-sample** selection, so read its values as optimistic.
- "Top 30 / top 300" counts the mitoses among the highest-scoring candidates of the uncapped `mo_pad61` deep pool, as a
  mean per click.
- The **tail** is the number of unannotated nuclei in the whole ROI that are darker than the ROI's 5th-darkest mitosis
  (`od51`). It is click-free.

| ROI | domain | A3 P@10 tm / od | A3 P@20 tm / od | C1 P@10 | mitoses (unanimous) | tail | mitoses in top 30 / 300 |
|---|---|---|---|---|---|---|---|
| `211` | canine lung | 0.00 / 0.00 | 0.05 / 0.00 | 0.00 | 20 (70%) | **327** | 0.7 / 3.8 |
| `242` | canine lung | 0.08 / 0.10 | 0.07 / 0.10 | 0.17 | 27 (85%) | **60** | 1.2 / 5.8 |
| `524` | human melanoma | 0.13 / 0.18 | 0.10 / 0.12 | 0.22 | 15 (60%) | **54** | 2.5 / 7.0 |
| `546` | human melanoma | 0.15 / 0.10 | 0.12 / 0.07 | 0.13 | 28 (54%) | **55** | 3.0 / 10.7 |
| `462` | canine STS | 0.20 / 0.10 | 0.17 / 0.13 | 0.35 | 15 (87%) | **83** | 3.7 / 8.0 |
| `527` | human melanoma | 0.23 / 0.33 | 0.20 / 0.32 | 0.40 | 43 (67%) | 6 | 3.7 / 16.2 |
| `410` | canine STS | 0.25 / 0.27 | 0.14 / 0.24 | 0.30 | 17 (76%) | 9 | 3.8 / 6.2 |
| `245` | canine lymphoma | 0.30 / 0.13 | 0.32 / 0.12 | 0.37 | **90** (87%) | **47** | 7.7 / 24.5 |
| `289` | canine lymphoma | 0.32 / 0.33 | 0.22 / 0.25 | 0.48 | 18 (67%) | 4 | 3.8 / 9.7 |
| *edge:* `509` | human melanoma | 0.35 / 0.22 | 0.33 / 0.31 | 0.67 | 18 (100%) | 6 | 7.2 / 13.2 |
| *edge:* `201` | canine lung | 0.35 / 0.33 | 0.27 / 0.36 | 0.32 | 18 (**39%**) | 5 | 5.7 / 13.0 |
| *edge:* `339` | canine MCT | 0.35 / 0.45 | 0.28 / 0.30 | 0.57 | 15 (67%) | 5 | 6.5 / 10.0 |
| *edge:* `404` | human NET | 0.37 / 0.57 | 0.38 / 0.47 | 0.67 | 61 (84%) | 2 | 9.2 / 32.7 |
| *edge:* `364` | human NET | 0.40 / 0.48 | 0.22 / 0.27 | 0.40 | 23 (**30%**) | 4 | 5.5 / 15.7 |
| *reference:* `459` | canine STS | 1.00 / 0.90 | 0.98 / 0.87 | 0.98 | 131 (93%) | 1 | 28.8 / 81.3 |

**Reading the table**
- ROIs outside the five worst average 13.7 mitoses in the top 30. The other 40 ROIs' median tail is 1, and only one non-failing ROI (`235`,
  tail 50) is at or above 40.
- At P@20, the rule "< 0.30 on both A3 axes" picks the same set minus `527` and `245`, plus `364`.

**Domain is not the right unit**
- The failing nine are:
  - 3 of 7 melanoma ROIs;
  - 2 of 7 each in lung, STS and lymphoma;
  - none of the mast-cell, breast or neuroendocrine ROIs.
- A Kruskal–Wallis label-permutation test on the current arms gives p = 0.27 (A3 tm P@10), 0.10 (A3 od P@10), 0.11 and
  0.12 (A3 P@20), and 0.065 (C1 P@10). This repeats the 2026-09-30 finding that domain levels are not established.
- **The literature points the same way.**
  - MIDOG++'s trained baseline finds melanoma among the *easiest* domains, with in-domain F1 0.82. It finds
    neuroendocrine the hardest, at 0.58 [Aubreville 2023b, Table 4].
  - Our neuroendocrine ROIs are fine, plausibly in part because the repo's inclusion bar (≥ 15 mitoses) excludes the
    typical neuroendocrine ROI (untested). Its median is 4 mitoses per ROI [Aubreville 2023b, Table 2].

**Two caveats**
- **Selection.** The failing ROIs were chosen on their own poor P@10, so their features are descriptive and subject to
  regression to the mean.
- **C1 regressions.** C1 does not lift every ROI. It drops `233` (0.72 → 0.43), `499` (0.77 → 0.60) and `220`
  (0.67 → 0.53).

---

## 2. What the measurements show

### 2.1 The darkness tail and prevalence are separate drivers

Spearman ρ across the 49 ROIs. The tail count and the click-free ranking both use `od51`, the statistic `chromatin_od`
ranks by. So the `A3 od` column is close to definitional, and the evidence is in the `tm_score` and C1 columns.

| | A3 tm P@10 | A3 od P@10 | A3 tm P@20 | C1 P@10 |
|---|---|---|---|---|
| tail (unannotated darker than 5th-darkest mitosis) | −0.75 | −0.86 | −0.75 | −0.79 |
| partial, given mitosis count | −0.70 | −0.83 | −0.75 | −0.76 |
| mitosis count, partial given tail | +0.60 | +0.32 | +0.76 | +0.55 |
| P@10 of a click-free whole-ROI `od51` ranking of nuclei | +0.79 | +0.89 | +0.76 | +0.80 |

**Weaker or non-explanatory measures** (Spearman with A3 tm P@10):
- the share of top-10 FPs beside white space: −0.15;
- the share of mitoses with unanimous labels: +0.12;
- mitosis blob area: none.

**Darkness itself is not missing in the failing ROIs.**
- Whole-ROI `od21` AUC, mitosis versus other nuclei, has a median of 0.93 (range 0.88–0.99) in the failing nine,
  against 0.97 in the other 40.
- What differs is the **extreme tail**: hundreds of unannotated objects darker than the mitoses.
  - `211`: 2,162 objects darker than its *median* mitosis.
  - `242`: 1,535.
  - `410`: 1,155.
  - `245`: 633.
  - Median of the failing nine: 324. Median of the other 40: 28.
- The 99th percentile of ordinary-nucleus `od21`, as a multiple of the median mitosis, has a median of 1.03 in the
  failing nine, against 0.90 elsewhere. It reaches 1.21 / 1.17 / 1.09 / 1.06 in `211` / `410` / `242` / `546`.
- So in most failing ROIs, the top 1% of non-mitotic nuclei is darker than a typical mitosis.

**Annotation status of the top-10 FPs.**
- In good ROIs, most of the remaining top-10 FPs are *annotated* look-alikes (MIDOG++ hard negatives). The share
  correlates with P@10 at ρ +0.67.
- In the failing nine, they are almost all *unannotated*: 0–26% look-alike, and ≈ 0% look-alikes with any mitotic vote.
- **Caveat:** MIDOG++ annotated hard negatives non-exhaustively, aiming for about as many as mitoses [Aubreville 2023b].
  Low-count ROIs therefore have fewer look-alike labels for an FP to land on.
- The montages show what the unannotated FPs are (§4). By eye, they are mostly pigment, smooth dark round nuclei and
  crowded hyperchromatic nuclei, pending pathologist review.

### 2.2 What the measurements rule out

- **Label noise is not why these ROIs fail.** Their FPs are not contested look-alikes.
  - Two edge cases are the exception, where measured precision is understated: `201`, with 39% of mitoses unanimous, and
    `364`, with 30%.
  - Label noise does cap every ROI (§5).
- **Missing mitoses are not the cause.** `nucleus_blobs` puts a nucleus within 7.5 µm of 90–100% of GT mitoses in the
  failing nine (77–100% across all 49).
  - The deep pool covers about 93% of each ROI.
  - The mitoses are in the candidate set, but ranked deep: the median first-match rank in the five worst ROIs is 623,
    against 126 elsewhere.
- **Template geometry is not the lever for these ROIs.**
  - A ranking that never sees the click predicts the pipeline's per-ROI P@10 at ρ 0.79–0.89.
  - This is consistent with the earlier measured finding that template geometry is second-order once a re-rank follows
    (`2026-09-30` deep-pool and two-stage re-rank work), and that template gains are mostly order inside the top 30.
    It is a between-ROI correlation, so it does not show that no template edit could change a given ROI's pool.

---

## 3. Why: the mechanism, from first principles and the literature

### 3.1 What the score measures

OpenCV defines `TM_CCOEFF` as

    R(x,y) = Σ T′(x′,y′) · I′(x+x′, y+y′)

where T′ is the template minus its mean and I′ is the window minus its mean [OpenCV `imgproc.hpp`]. It has no
denominator; only `TM_CCOEFF_NORMED` divides by √(ΣT′²·ΣI′²).

**Derivation.** Because ΣT′ = 0, the window-mean term cancels: Σ T′·I′ = Σ T′·I. Lewis states the same identity for the
NCC numerator: "Since t′ has zero mean and thus zero sum the term f̄ Σ t′ … is also zero" [Lewis 1995b].

**Consequences**
- `TM_CCOEFF` is a **linear filter with a zero-mean kernel**. It ignores an additive offset, but its output scales
  linearly with local contrast. The repo measured exactly that: contrast 0.7 / 0.3 / 0.1 gives 0.700 / 0.300 / 0.101 of
  the full score (`midog_utils_full/chromatin.py` docstring).
- For a two-level template (value a on a core of n_c pixels, b on a surround of n_s, N = n_c + n_s), it reduces to
  the following (derived here; the same form as the repo's mass/shape decomposition):

      R = (a − b) · (n_c·n_s / N) · (mean_core(H-OD) − mean_surround(H-OD))

  So the score is *how much more hematoxylin OD sits in a nucleus-sized core than around it*.
  - The repo measured that this mass term carries about 4× what the template's shape term does: deleting mass costs
    12 pts; a flat box filter recovers 92–93%.
  - **Anything with a darker core, or a whiter surround, outranks a mitosis.**
- This is the `TM_CCOEFF` form of the failure Lewis describes for unnormalised correlation: "the correlation between the
  template and an exactly matching region in the image may be less than the correlation between the template and a
  bright spot" [Lewis 1995a]. Here the bright spot is a high-OD-contrast spot.
- `chromatin_od` (the mean of the darkest 10% of a 51 px window) is a darkness statistic by construction. So both ranking
  axes read the same physical quantity, and the ROIs where darkness fails are the ROIs where both fail. That is the
  ρ 0.79–0.89 in §2.1.

### 3.2 Why "darkest" is not "mitotic": the biology of hematoxylin

**Mitoses are dark for a reason other cells share**
- Hemalum colours DNA [Kiernan 2018, abstract].
- Mitotic chromatin "further condenses and individualizes to discrete chromosomes" [Antonin & Neumann 2016, abstract].
- Pyknosis also condenses chromatin. "Pyknosis is the result of chromatin condensation and this is the most
  characteristic feature of apoptosis". The apoptotic cell shows "dense purple nuclear chromatin fragments" [Elmore 2007].
- Mitotic figures' "hyperchromacity [is] the most salient feature. Moreover, many other cellular components can have a
  similar hyperchromatic appearance, such as apoptotic or necrotic nuclei, compressed nuclei, 'junk' particles and other
  artifacts" [Veta 2015].
- Darkness therefore enriches for mitoses, but it does not identify them.

**Named imposters** (all show up in our montages, §4)
- "apoptotic bodies, hyperchromatic nuclei, karyorrhectic debris, inflammatory cells, two cells in close proximity, and
  tissue artifacts"; "lymphocytes were most likely to resemble mitosis in prometaphase, and stromal cells were likely to
  mimic mitosis in metaphase" [Lin 2024].
- "apoptotic and necrotic cells, hyperchromatic nuclei, inflammatory cells, and ink for margin marking"
  [Aubreville 2026, preprint].

**What does separate them, per the literature, is not darkness**
- Mitoses "appear as hyperchromatic objects without nuclear membrane and hairy extensions of nuclear material"
  [Veta 2016].
- In AMIDA13, "The majority of the false positives produced by the methods with lower performance were dark objects that
  lack other characteristics of mitotic figures such as protrusions around the edges. On the other hand, many of the
  mitotic figures that had less intensive staining were not correctly detected" [Veta 2015].

**AMIDA13 had the same two failing kinds of case**
- "Case #9 contained many dark nuclei that are not mitotic figures, thus the methods that extensively rely on the
  staining intensity as a feature had very low performance for this case". That is our tail driver.
- Cases #1 and #7 "have very low mitotic activity, and the few ground truth mitoses have atypical appearance (lack of
  hyperchromacity)". That is our prevalence driver [Veta 2015].

### 3.3 Prevalence and the extreme tail: why few mitoses hurt so much

P@10 depends only on the top 10 of 6k–35k candidates, so it is set by the **upper tail** of the non-mitosis score
distribution, not by AUC.

- "a large change in the number of false positives can lead to a small change in the false positive rate used in ROC
  analysis" [Davis & Goadrich 2006].
- "Only PRC changes with the ratio of positives and negatives" [Saito & Rehmsmeier 2015].
- P@k "does not average well, since the total number of relevant documents for a query has a strong influence on
  precision at k" [Manning et al. 2008, §8.4].

**The arithmetic (a derivation, idealised).** Take an ROI with M evaluable mitoses among N other nuclei. For P@10 ≥ 0.5,
at most 5 of the N may score above the 5th-best mitosis.
- That is an operating point of FPR ≤ 5/N and TPR ≥ 5/M.
- With M = 14 and N = 15,000, that is TPR 0.36 at FPR 0.0003.

Under an equal-variance Gaussian model, the AUC required for P@10 = 0.5 and 0.8 is:

| evaluable mitoses M | AUC for P@10 = 0.5 | AUC for P@10 = 0.8 |
|---|---|---|
| 14 | 0.984 | 0.997 |
| 26 | 0.963 | 0.987 |
| 89 | 0.900 | 0.948 |
| 130 | 0.876 | 0.932 |

The same detector that is good enough in `459`, with 130 mitoses, is nowhere near good enough with 14. That is why P@10
tracks mitosis count (ρ +0.66) whatever the arm.

**The model is optimistic.** `245` and `289` have whole-ROI `od21` AUCs (0.94, 0.99) *above* their idealised
requirement, yet fail. Their non-mitosis tail is heavier than Gaussian: a handful of very dark objects. This is why §2.1's
tail count is the right statistic and AUC is the wrong one. It is also why the repo's earlier "od has AUC 0.84 but isn't
a top 10" result was not a contradiction.

**MIDOG 2022 notes the same effect on F1:** "especially for cases with a low true MFs count, the F1 score is strongly
influenced by the number of false positives" [Aubreville 2024].

### 3.4 Matched-filter theory: one click's template is the wrong filter in clutter

**What the theory requires**
- "The matched filter is the optimal filter for detecting a known waveform in stationary Gaussian noise" [Allen 2012].
  Correlating with the clicked cell itself is the **non-prewhitening matched filter** (NPW): t = Δḡᵀg [Li 2021].
- With correlated backgrounds, the optimal *linear* detector is the Hotelling observer, t = (K⁻¹Δḡ)ᵀg. It "is optimal
  among all linear observers in the sense that it maximizes the signal-to-noise ratio" [Li 2021]. It is the
  "prewhitening matched filter" [Poludniowski 2025].
- When the background is lumpy, "the ideal observer becomes nonlinear". The NPW model "predicts that performance will …
  drop precipitously with increasing lumpiness" [Rolland & Barrett 1992, abstract].
- When signal parameters are unknown (phase, orientation, size), the standard extension is a template bank with a
  maximised likelihood ratio. That is "in general not an optimal detection statistic" [Röver 2011].

**How tissue breaks those requirements**
- Template matching "is not invariant with respect to imaging scale, rotation, and perspective distortions"
  [Lewis 1995b].
- "Intraclass variability poses new problems for template matching" [Brunelli 2009, ch. 1].
- Phase matters especially. MIDOG 2022's largest category of mitoses missed on H&E was "prophase with resemblance to
  apoptosis" (31.23%) [Aubreville 2024].

**What that means here (a derivation, not measured).**
- The background in these ROIs is mostly *other nuclei*. Its covariance K has its largest variance along exactly the
  direction the template mostly measures: a nucleus-sized blob of hematoxylin mass.
- Prewhitening by K⁻¹ would *down-weight* that direction. It would put the weight on the components where mitoses differ
  from the clutter: outline, protrusions, texture.
- The NPW template does the opposite. That makes no difference where mitoses are the darkest objects (`459`), and it is
  fatal where they are not (`211`, `242`, `546`).
- So the theory and the measured tail statistic say the same thing. **In high-clutter ROIs the useful signal is
  orthogonal to mass**, and a correlation template built from one click cannot isolate it.

### 3.5 Pigment: absorbers the colour model cannot represent

**Colour deconvolution assumes three stains**
- It assumes Beer–Lambert absorbers ("OD for each channel is linear with the concentration of absorbing material") and
  up to three known stain vectors [Ruifrok & Johnston 2001].
- "Neutral (grey) dyes cannot be unmixed with this method" [Landini 2021, read through a summariser].
- Ruifrok & Johnston measured off-basis leakage: saturated DAB contributed "an average of 65%" of the measured
  hematoxylin IOD.

**What `rgb2hed` does with a neutral absorber (computed here)**
- skimage 0.24's `rgb2hed` maps an equal-OD (neutral black) absorber a to **H = 1.21a**, D = 0.88a, and E < 0, clipped
  to 0. Here a is the OD in each single channel; the comparison with hematoxylin is per unit vector norm, below.
- So a neutral absorber projects about 70% of its OD (by vector norm) onto H; a pure hematoxylin vector of the same
  norm projects 100%. Black pigment is, however, far more opaque than chromatin, so its *absolute* H value exceeds a
  mitosis's. In `211` the top-10 FPs' `od21` is 1.5× the median mitosis.
- Melanin absorbs "broad and monotonic extending from the UV to the near-IR" [Meredith 2006]. It registered in the
  deconvolved DAB channel of melanoma sections before immunostaining [Rexhepaj 2013].
- MIDOG++ names melanoma pigment an imposter source [Aubreville 2023b].
- Anthracosis was found histologically in 51.8% of one sample of dog lungs [Mohamed 2022, abstract]. **No source shows
  it in MIDOG++ ROIs.**

**Measured in `211`**
- The top-10 FPs' darkest pixels have a DAB-channel OD of 0.152, against 0.034 for the mitoses.
- Their blue/red OD ratio is 0.87, against 0.72.
- 52% of them have a DAB OD above twice the mitoses' 95th percentile and above 0.05.
- In the montage they are black granular clumps.

**Pigment explains only part of the list**
- It is a minority cause in `289` (11%) and `339` (12%), and is visible in `410` and `291`.
- In the failing melanoma ROIs (`524`, `546`, `527`) the colour test finds ~0%, and the montages show dark nuclei, not
  melanin.
- So pigment is the cause in `211` and incidental elsewhere.

### 3.6 White surround

**Mechanism (derived in the repo).** By the formula in §3.1, a white surround raises the score exactly as much as a
darker core does. The A2/A3 floor clip was the fix. It removed white-adjacent top-10 FPs (316 → 240) without hitting
lumen mitoses.

**Effect on the failing ROIs**
- After A3, white space is still the dominant FP context in `242`, where 53% of FPs are white-adjacent. This is
  alveolar or lumen space.
- It is also substantial in `289` (37%).
- Four of the seven mast-cell ROIs have 77–100% white-adjacent FPs but still pass, because they have enough mitoses.
- White matters only when it coincides with low prevalence.

---

## 4. Per-ROI diagnosis

**Source:** the `montage_*.jpg` files in this folder. These are by-eye morphology readings,
**not pathologist-confirmed**. The words "consistent with" mark interpretation.

| ROI | main driver | what the top-10 FPs are | what the mitoses look like |
|---|---|---|---|
| `211` | tail (327) + low count | Black granular pigment clumps (52% by the colour test), some at the ROI edge; consistent with anthracotic pigment | Moderately dark, irregular figures inside large cells; not the darkest objects |
| `242` | tail (60) + low count + white | Small, smooth, very dark round nuclei beside alveolar/lumen space (53% white-adjacent); consistent with lymphocytes, plasma cells or pyknotic cells | Irregular clumps in larger cells, often paler than the FPs |
| `524` | tail (54) + 14 evaluable mitoses + label noise (60% unanimous) | Overlapping and crowded hyperchromatic tumour nuclei (OD adds where nuclei overlap, by Beer–Lambert); one cell in a vessel; no melanin seen | Small dark clumps in the same crowded field, not darker than their neighbours |
| `546` | tail (55) + label noise (54% unanimous) | Smooth, round, dark nuclei, several with perinuclear halos | Spiky, irregular plates of similar darkness: separable by outline, not by OD |
| `462` | tail (83) + 14 evaluable mitoses | Condensed dark nuclei in clear halos or vacuoles (25% are annotated look-alikes) | Irregular bars |
| `527` | moderate tail (165 objects darker than its median mitosis) | Large, pleomorphic, hyperchromatic interphase nuclei | Irregular clumps of the same size and darkness |
| `410` | low count + 1,155 objects darker than its median mitosis | Dark nuclei in clear halos; nuclei beside golden-brown granular pigment (consistent with hemosiderin); hyperchromatic spindle nuclei | Several reach the top 10 |
| `245` | tail (47) despite 90 mitoses | Round dark nuclei with halos and tiny black specks, inside a uniformly dark lymphoid sheet | Subtle, small dark clumps scarcely darker than the sheet; contrast collapse |
| `289` | low count + white (37%) | Dark nuclei at lumen edges with red cells; brown-black pigment clusters (one annotated as a look-alike) | Small dark bars in a lymphoid sheet; several do reach the top 10 |
| `509` (edge) | order, not pool | Angulated, smudged dark nuclei in eosinophilic surroundings | Large, clear metaphase plates that are *paler* than the FPs; C1 lifts it to 0.67 |
| `201` (edge) | label noise (39% unanimous) + low count | 33% are annotated look-alikes; elongated dark nuclei | — |
| `404` (edge) | contested labels | Lymphocyte-like clusters, a blue-green artefact, 47% look-alikes (25% contested) | — |
| `364` (edge) | label noise (30% unanimous) | Crushed, elongated nuclei at the tissue edge | — |

### Summary of failure types

**A. Dark-object tail.**
- **ROIs:** `211`, `242`, `524`, `546`, `462`, `245`.
- **Why they fail:** dozens to hundreds of non-mitotic objects are darker than the ROI's few best mitoses.
- **What the dark objects are:** pigment, small dark round cells, apoptotic-looking halo nuclei, or crowded
  hyperchromatic nuclei.
- **What re-ranking can do:** no darkness or mass ranker can fix these, and re-ranking a 100- or 300-deep pool cannot
  reach 0.5 because the mitoses aren't there (§5).

**B. Prevalence-limited with a modest tail.**
- **ROIs:** `289`, `410`, `527`, and edge cases `509`, `339`.
- **Why they fail:** about 15–18 mitoses, so 5 of them must beat nearly everything.
- **What re-ranking can do:** this is where re-ranking helps; C1 lifts `509` to 0.67, `339` to 0.57 and `289` to 0.48.

**C. Label-limited.**
- **ROIs:** `201`, `364`, and part of `524` and `546`.
- **Why they fail:** measured precision understates the real precision when only 30–60% of mitoses are unanimous.

---

## 5. How good can this get? Four separate ceilings

**1. Labels.**
- In MIDOG++ "in approximately 20% of cases, the experts disagreed" (Table 3: 16.6–23.2% by domain)
  [Aubreville 2023b].
- Object-level Dice between expert pathologists is 0.57–0.67 [Veta 2016], and 0.483 between AMIDA13's two pathologists
  [Veta 2015].
- In AMIDA13, 29% [0.23, 0.36] of the best method's false positives were re-annotated as mitoses [Veta 2015].
- Expert H&E reads scored against PHH3-assisted truth reached F1 0.70–0.76 [Aubreville 2024].
- So no method's measured precision should be expected near 1.0. This ceiling applies everywhere and is not why the nine
  fail.

**2. Prevalence.**
- P@K ≤ min(1, (M−1)/K) [Manning 2008, §8.4, the same arithmetic].
- For `524` and `462`, with 14 evaluable mitoses: P@20 ≤ 0.70 and P@30 ≤ 0.47.
- In the 2026-09-30 measurement, P@30 ≥ 0.8 was impossible in lung, melanoma and STS.

**3. Pool (measured, `mo_pad61` deep pool).**
- A perfect re-rank of the top 30 / 100 / 300 by score gives pooled P@10 of 0.80 / 0.90 / 0.96.
- The five worst ROIs hold only 2.2 / 4.2 / 7.1 mitoses at those depths, so their P@10 is at most about
  0.22 / 0.42 / 0.71 on average, even with an oracle.
- The shape (solidity) re-rank measured on that pool added at most about half a mitosis there.

**4. Method class (literature; F1 is not P@10, but the ordering is informative).**

| System | Method type | Result | Source |
|---|---|---|---|
| ICPR 2012, IPAL | hand-crafted features + **trained** classifier | F1 0.718 | [Cireşan 2013, Table 2; Irshad 2013] |
| ICPR 2012 | CNN | F1 0.782 | same |
| AMIDA13, DTU | hand-crafted + SVM | F1 0.483 | [Veta 2015] |
| AMIDA13, IDSIA | CNN | F1 0.611 | [Veta 2015] |
| MIDOG++ RetinaNet, in-domain | trained detector | lung 0.66, STS 0.70, lymphosarcoma 0.79, melanoma 0.82, MCT 0.85, NET 0.58 | [Aubreville 2023b] |
| MIDOG 2025 participants | trained models | precision 0.805 (hotspot) → 0.614 (random) → **0.400** (imposter-rich); recall flat at 0.66–0.69 | [Aubreville 2026, preprint] |

**What the literature table means**
- Even trained state-of-the-art detectors lose precision, not recall, exactly where imposters are dense.
- The hand-crafted rows are already *trained* classifiers on 143+ features, and they still trail CNNs.
- A no-ML darkness or mass ranker in a low-prevalence, imposter-rich ROI should be expected well below those numbers.
- MIDOG++ ROIs are hotspots chosen for "high mitotic density" [Aubreville 2023b]. On random regions the product will
  see lower prevalence, so *more* ROIs will look like the nine.

**Bottom line**
- The 49-ROI pooled P@10 can still rise; the top-100 ceiling is 0.90.
- A per-ROI floor of 0.8 with low variance across ROIs is not attainable with darkness/mass ranking. For ROIs like
  `211`, `242`, `546` and `524` it would need a discriminator orthogonal to darkness with an AUC around 0.98+ over the
  whole ROI.
- The equal-variance model makes the requirements in §3.3 look *easier* than they are.

---

## 6. What follows (brief)

1. **Pigment is the one cheap, principled gate.**
   - Beer–Lambert says a neutral or brown absorber has an OD vector off the hematoxylin axis. The DAB channel or the
     blue/red OD ratio flags it (§3.5).
   - It is dominant only in `211`. It should be tested as a pre-registered arm, judged on whether it removes pigment
     without removing mitoses (mitoses' DAB-channel medians range 0.000–0.075 across ROIs and scanners).
2. **The rest needs features orthogonal to darkness.**
   - The literature names them: "hairy extensions", no nuclear membrane, "protrusions around the edges"
     [Veta 2016; Veta 2015]. Hotelling theory says these are the directions that matter when clutter shares the mass.
   - The repo's solidity measure (a dilated mask, AUC about 0.6) is not that cue. The pathologist-cue features have
     never been computed; `CANDIDATE_FEATURE_EXTRACTION_PROMPT.md` is where they would be.
   - Even done well, trained hand-crafted systems reached F1 0.48–0.72, so expect a modest gain. A trained classifier is
     the step change, subject to the commercial-licence constraint.
3. **Report per ROI, not only pooled.**
   - Put mitosis count, the tail count and the prevalence ceiling beside every P@K.
   - Judge generator changes by mitoses in the top 30/100 on the nine ROIs here.
   - The tail statistic is a diagnostic that needs GT. It is not something the product can compute.

---

## 7. Caveats

- **Post hoc.** The rule, the tail statistic and the montages were all chosen after the failing ROIs were known.
- **"Unannotated" ≠ "not mitotic".** Negatives are nuclei more than 7.5 µm from *any* annotation. Unannotated real
  mitoses would inflate the tail. MIDOG++ mitosis labels were screened with a low-threshold model to find misses
  [Aubreville 2023b], so this is probably small, but it is not zero [Veta 2015; Aubreville 2024].
- **The pigment test is a heuristic.** A looser version (blue/red ratio *and* DAB above the mitoses' 95th percentile)
  also flagged most FPs in some *good* ROIs, e.g. `301` and `459`. Only the `211` signal is strong enough to quote.
- **Morphology readings are by eye.** "Consistent with lymphocytes / pyknotic / hemosiderin / anthracotic" needs a
  pathologist.
- **The equal-variance model is idealised.** §3.3 shows it is optimistic on these data.

---

## Sources

All were read in this session. **[abstract]** means only the abstract was read.

**Datasets and challenges**
- **[Aubreville 2023b]** Aubreville M, Wilm F, Stathonikos N, et al. A comprehensive multi-domain dataset for mitotic
  figure detection. *Sci Data* 10:484 (2023). doi:10.1038/s41597-023-02327-4. Read with its 2024 correction,
  doi:10.1038/s41597-024-03548-x.
- **[Aubreville 2023a]** Aubreville M, et al. Mitosis domain generalization in histopathology images — The MIDOG
  challenge. *Med Image Anal* (2023). doi:10.1016/j.media.2022.102699. Read as arXiv 2204.03742 v1.
- **[Aubreville 2024]** Aubreville M, et al. Domain generalization across tumor types, laboratories, and species —
  insights from the 2022 edition of the Mitosis Domain Generalization Challenge. *Med Image Anal* (2024).
  doi:10.1016/j.media.2024.103155. Read as arXiv 2309.15589 v2.
- **[Aubreville 2026]** Aubreville M, et al. Mitosis Detection in the Wild: Multi-Tumor and Context-Aware Generalization
  in the MIDOG 2025 Challenge. arXiv 2606.07368 v1 (preprint).
- **[Veta 2015]** Veta M, van Diest PJ, Willems SM, et al. Assessment of algorithms for mitosis detection in breast
  cancer histopathology images. *Med Image Anal* 20:237–248 (2015). doi:10.1016/j.media.2014.11.010. Read as
  arXiv 1411.5825.
- **[Veta 2016]** Veta M, van Diest PJ, Jiwa M, Al-Janabi S, Pluim JPW. Mitosis counting in breast cancer: object-level
  interobserver agreement and comparison to an automatic method. *PLoS ONE* 11:e0161286 (2016).
- **[Cireşan 2013]** Cireşan DC, Giusti A, Gambardella LM, Schmidhuber J. Mitosis detection in breast cancer histology
  images with deep neural networks. MICCAI 2013, LNCS 8150. doi:10.1007/978-3-642-40763-5_51.
- **[Irshad 2013]** Irshad H. Automated mitosis detection in histopathology using morphological and multi-channel
  statistics features. *J Pathol Inform* 4:10 (2013). doi:10.4103/2153-3539.112695.
- **[Roux 2013]** Roux L, et al. Mitosis detection in breast cancer histological images: an ICPR 2012 contest.
  *J Pathol Inform* 4:8 (2013). doi:10.4103/2153-3539.112693.
- **[Lin 2024]** Lin S, et al. The 1000 Mitoses Project: a consensus-based international collaborative study on mitotic
  figures classification. *Int J Surg Pathol* 32:1449–1458 (2024). doi:10.1177/10668969241234321.

**Template matching and detection theory**
- **[OpenCV]** `modules/imgproc/include/opencv2/imgproc.hpp` (4.x), the `TemplateMatchModes` definitions.
- **[Lewis 1995a]** Lewis JP. Fast Template Matching. Vision Interface 95, pp. 120–123.
- **[Lewis 1995b]** Lewis JP. Fast Normalized Cross-Correlation, the expanded Industrial Light & Magic version
  (scribblethink.org).
- **[Brunelli 2009]** Brunelli R. *Template Matching Techniques in Computer Vision: Theory and Practice.* Wiley (2009).
  Only the publisher's chapter-1 sample was read.
- **[Allen 2012]** Allen B, et al. FINDCHIRP. *Phys Rev D* 85:122006 (2012). doi:10.1103/PhysRevD.85.122006.
- **[Röver 2011]** Röver C. A Student-t based filter for robust signal detection. *Phys Rev D* 84:122004 (2011).
  doi:10.1103/PhysRevD.84.122004.
- **[Li 2021]** Li K, Zhou W, Li H, Anastasio MA. Assessing the impact of deep neural network-based image denoising on
  binary signal detection tasks. *IEEE TMI* 40:2295–2305 (2021). doi:10.1109/TMI.2021.3076810.
- **[Poludniowski 2025]** Poludniowski G, Titternes R, Thor D. A practical approach to the spatial-domain calculation of
  nonprewhitening model observers in computed tomography. *Med Phys* 52:2106–2122 (2025). doi:10.1002/mp.17599.
- **[Rolland & Barrett 1992]** Rolland JP, Barrett HH. Effect of random background inhomogeneity on observer detection
  performance. *JOSA A* 9:649–658 (1992). doi:10.1364/josaa.9.000649. [abstract]

**Evaluation under imbalance**
- **[Davis & Goadrich 2006]** Davis J, Goadrich M. The relationship between precision-recall and ROC curves. ICML 2006.
  doi:10.1145/1143844.1143874.
- **[Saito & Rehmsmeier 2015]** Saito T, Rehmsmeier M. The precision-recall plot is more informative than the ROC plot
  when evaluating binary classifiers on imbalanced datasets. *PLoS ONE* 10:e0118432 (2015).
  doi:10.1371/journal.pone.0118432.
- **[Manning 2008]** Manning CD, Raghavan P, Schütze H. *Introduction to Information Retrieval*, §8.4. CUP (2008).

**Colour deconvolution and pigment**
- **[Ruifrok & Johnston 2001]** Ruifrok AC, Johnston DA. Quantification of histochemical staining by color
  deconvolution. *Anal Quant Cytol Histol* 23:291–299 (2001). Read as the author manuscript.
- **[Landini 2021]** Landini G, Martinelli G, Piccinini F. Colour deconvolution: stain unmixing in histological imaging.
  *Bioinformatics* 37:1485–1487 (2021). doi:10.1093/bioinformatics/btaa847. Read through a summariser.
- **[Meredith 2006]** Meredith P, et al. Towards structure–property–function relationships for eumelanin. *Soft Matter*
  2:37–44 (2006). doi:10.1039/b511922g.
- **[Rexhepaj 2013]** Rexhepaj E, et al. A texture based pattern recognition approach to distinguish melanoma from
  non-melanoma cells in histopathological tissue microarray sections. *PLoS ONE* 8:e62070 (2013). doi:10.1371/journal.pone.0062070.
- **[Mohamed 2022]** Mohamed RA, et al. *Am J Vet Res* 83(10) (2022). doi:10.2460/ajvr.22.05.0077. [abstract]

**Biology of nuclear darkness**
- **[Kiernan 2018]** Kiernan JA. *Biotech Histochem* 93:133–148 (2018). doi:10.1080/10520295.2017.1399466. [abstract]
- **[Antonin & Neumann 2016]** Antonin W, Neumann H. *Curr Opin Cell Biol* 40:15–22 (2016).
  doi:10.1016/j.ceb.2016.01.013. [abstract]
- **[Elmore 2007]** Elmore S. Apoptosis: a review of programmed cell death. *Toxicol Pathol* 35:495–516 (2007).
  doi:10.1080/01926230701320337.

**Not verified, so not used for any claim:** the imposter criteria in Donovan et al. 2021 (*Vet Pathol*, paywalled;
abstract only); hand-crafted-only F1 in Wang et al. 2014; Turin 1960; Barrett & Myers.
