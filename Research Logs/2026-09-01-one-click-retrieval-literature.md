# One-click mitosis retrieval: literature review and the fastest routes to adapt

Date: 2026-09-01
Question asked: *"my agents keep telling me native picture thresholding is better, but I want
template matching or something more sophisticated. Pathologist clicks a mitotic cell, the
algorithm finds similar structures in ~1 minute, the human only verifies."*
Follows: `2026-08-31-premise-test-results.md`, `2026-08-31-chromatin-density-rerank.md`
Status: review + recommendation. Nothing here had been implemented or measured beyond the
CPU timing table in §3 when this was written.

> **Superseded in part by `2026-09-01-click-ranking-experiment.md`.** §7's experiment has
> since been run. Two of this document's conclusions did not survive it:
> **§6 recommendation 4 is refuted** -- frozen pathology foundation-model embeddings with
> cosine to the click (`lunit_dino_vits16`, `kaiko_vits16`) are 3-42x *worse* than a
> seedless blob detector on reading depth, not a "ship it today" path; and **§7's own
> gate needed the decision metric it names but a harder bar than it sets** (the best
> seedless arm, not chromatin). What did survive, and strongly: the two-stage reframe, the
> claim that the click's information had never been measured in a representation that could
> carry it, and recommendation 5 -- adapt the encoder to the task. Read that log before
> acting on §6 here.

---

## 1. The reframe that reconciles your instinct with your own data

Your agents and you are each right about a *different stage*, and the argument only looks
like a disagreement because the pipeline fuses the two.

| stage | job | what your data already says |
|---|---|---|
| **A. Proposals** | enumerate every object worth showing | **Solved, seedlessly.** `nucleus_blobs` captures 96.6-100% of mitotic figures and 98-100% of look-alikes (`results/lookalike_auc_unbiased.csv`) in 7-11 s, with no click. The seeded correlation search never beat it on recall@budget at any z. |
| **B. Ranking** | order those candidates so the reader hits the mitoses first | **Wide open.** Chromatin density sits at AUC 0.693 (301) / 0.768 (246) against annotated look-alikes. That is the entire remaining problem. |

The premise test's verdict — "the click buys +0.05, less than the pipeline's deficit" — was
measured with the click acting through `cv2.matchTemplate(TM_CCOEFF_NORMED)`. Your own
chromatin log proves that statistic is invariant to `I -> aI + b`, i.e. blind to mean
intensity, which `results/morph_diag_bhattacharyya.csv` ranks as the **strongest**
mitotic-vs-ordinary-nucleus feature in all seven domains (0.925-3.703). So the measurement
is of *pixel template matching*, not of *the click*. **The information content of the click
has never been measured in a representation capable of carrying it.**

That is the review's organising claim: **keep the seedless proposal generator, and spend the
click on ranking.** Everything below is filed under stage A or stage B.

## 2. The latency budget is much softer than "one minute of compute"

Neither proposal generation nor candidate *embedding* depends on the click. Both are
functions of the ROI alone, so both belong in an ingest step that runs when the slide is
loaded (or overnight, or in the cloud). At click time the work is: embed one 64-96 px patch,
one matrix-vector cosine against ~20k cached vectors, sort. **Milliseconds.**

This matters for model selection: it removes "fast inference" as a constraint on the
*encoder*, and replaces it with "how long may ingest take, and on what hardware."

## 3. What this laptop can actually do (measured, 2026-09-01)

i7-8750H, 6 threads, no CUDA (Intel UHD 630; `torch.backends.mps` reports True but there is
no usable GPU backend here). torch 2.2.2 / torchvision 0.17.2, random weights, forward pass
only — **timing only, these say nothing about accuracy.** 20k candidates is the real
per-ROI count on 301.tiff.

| model @ crop | patches/s | 20k candidates | verdict |
|---|---:|---:|---|
| ResNet-18 @ 64 px | 456 | **44 s** fwd / **61 s** end-to-end | borderline click-time; see below |
| ResNet-18 @ 96 px | 256 | 78 s | ingest-time |
| ViT-S/16 @ 96 px | 133 | 150 s | ingest-time |
| ViT-S/16 @ 224 px | 26 | 13 min | ingest-time only |
| ResNet-50 @ 224 px | 14 | 23 min | ingest-time only |
| ViT-B/16 @ 224 px | 8.8 | 38 min | overnight |
| ViT-g/H FMs (UNI2-h, Virchow2, H-optimus) | ~0.3-0.6 est. | **10-25 h est.** | not on this machine |

Detector timing, same machine: `fcos_resnet50_fpn` at 1024 px tiles is 1.14 s/tile, ≈ **87 s
for a whole 5412x7215 ROI** at 30% overlap. So a real mitosis detector over the full ROI is
an ingest-cost item on CPU, not a click-time one.

The forward-pass column excludes cropping and normalisation, so the ResNet-18 row was
re-measured end-to-end on the real 301.tiff array: **0.9 s to cut 20k 64 px crops + 60.0 s to
normalise and embed = 61 s**. Cropping is free; the tensor conversion and normalisation cost
as much again as the convolutions. Treat the other rows' totals as ~1.4x optimistic for the
same reason. (Note: the system `python3` here has numpy 2.4.6 against a torch built for
numpy 1.x, so `torch.from_numpy` raises `RuntimeError: Numpy is not available` — the same
numpy<2 pin `requirements-midog-utils.txt` documents for cv2/tifffile applies to torch.)

Two consequences worth stating plainly:

* A **small CNN encoder is nearly click-time-viable with zero infrastructure** — 61 s
  measured end-to-end for 20k candidates on the laptop you have — and comfortably so once
  proposals come from a detector (a few hundred candidates) rather than blobs.
* **Pathology foundation models are a rented-GPU-at-ingest decision, not a latency
  decision.** 20k patches through ViT-L on one A100/4090 is 2-4 minutes. The cost of using
  UNI/Virchow2/H-optimus is an ops cost, not a user-facing one.

---

## 4. Stage A — proposals (seedless). Ranked by what you'd gain.

### A1. Use the released MIDOG++-trained detector instead of blobs — **MIT licensed, weights public**

* `jonas-amme/FCOS_Inference_CLI` — FCOS with ResNet-18 / ResNeXt-50 / ResNeXt-101
  backbones, **trained on MIDOG++**, F1 **0.737-0.753** on the official MIDOG 2022 test set,
  MIT license, `download_weights.py`, CLI with ROI and WSI modes.
  <https://github.com/jonas-amme/FCOS_Inference_CLI>
* The same `FCOS_x101.ckpt` ships as the **MIDOG 2025 Track 1 reference algorithm**
  (`DeepMicroscopy/MIDOG25_T1_reference_docker`, release v1.0.0), i.e. it is the challenge's
  own baseline. The MIDOG 2025 baselines paper reports the detection baseline at F1 0.6883
  on the final test set (0.7472 in hotspots, 0.4161 on "challenging" ROIs).

  **Leakage warning, load-bearing:** these weights were trained on MIDOG++, which is where
  001-506 come from. Any recall/precision number you compute on 301.tiff or 246.tiff with
  them is optimistic and not publishable. Use it as a *product* component; measure it only on
  held-out domains or the official test split.

  Why it still matters: a proposal list of a few hundred detections at F1 ≈ 0.75 replaces a
  20k-candidate blob list. That changes the reading burden by two orders of magnitude before
  the click is even involved, and it makes the click's job "re-rank 300 things," which is a
  far easier problem than "re-rank 20,000."

### A2. Current SOTA architecture, if you ever train your own

Every strong MIDOG 2025 entry is **two-stage: high-recall detector → patch classifier for
false-positive reduction.**

* *Mitosis Detection, Fast and Slow* (MDFS) — MIDOG21 F1 0.747 (joint 1st), MIDOG22 **0.764
  (1st)**. EUNet (EfficientNet-B0 encoder) segments candidates at reduced resolution, then a
  refinement classifier. **~1.53 s per 2 mm² ROI**, <3 min per WSI. This is the speed
  reference for the whole field. <https://arxiv.org/abs/2208.12587>
* *A bag of tricks for real-time MF detection* — single-stage RTMDet-S, F1 0.81 on the MIDOG
  2025 preliminary test set, explicitly targeting real-time deployment; hard-negative mining
  on necrosis/debris. <https://arxiv.org/abs/2508.19804>
* *YOLO11x proposals + ConvNeXt-Tiny classifier* — F1 0.882 fused; precision 0.762 → 0.839
  from the second stage alone. <https://arxiv.org/abs/2509.02627>
* *RF-DETR* and *attention-guided FP correction* are the other 2025 Track 1 entries
  (<https://arxiv.org/pdf/2509.02599>, <https://arxiv.org/pdf/2509.02598>); best reported
  Track 1 numbers cluster at F1 0.84 on validation / MIDOG++ test.

### A3. Generic nucleus detectors (only if you want candidates, not mitoses)

CellViT, StarDist, Cellpose, and the lighter **NuLite** (1.55-10.4x faster than CellViT)
give instance masks for *every* nucleus. Against your `nucleus_blobs` — which already
captures 96.6-100% of both classes in 7-11 s — these buy morphology features, not recall.
Low priority.

---

## 5. Stage B — click-conditioned ranking. This is where your product lives.

### B1. Frozen pathology-FM embedding + cosine to the clicked patch — **the recommended core**

The literature says the representation is the thing:

* **Beyond Classification: Pathology FMs as Detection Encoders for Mitotic Figures**
  (<https://arxiv.org/html/2607.28007v1>) — UNI, UNI2-h, Virchow, Virchow2, H-optimus-0/1,
  **frozen**, with Faster R-CNN / RetinaNet / Deformable DETR heads, on MIDOG++. Best frozen
  FM (H-optimus-0 + RetinaNet) **F1 0.7718** vs an end-to-end ResNet-50 baseline's 0.7917 and
  a *frozen* ResNet-50's 0.7498. Two things to take from it: frozen FM features are
  spatially resolved enough for mitosis-scale detection, and they do not beat a trained
  ResNet-50 — so the FM's value here is that it needs **no training**, not that it is better.
  Their FPN strides were [7,14,28,56,112] px, i.e. a ~50 px nucleus is 3-7 tokens; check this
  on your own crops rather than assuming it transfers.
* **Benchmarking FMs for Mitotic Figure Classification**
  (<https://arxiv.org/abs/2508.04441>, MELBA 2026:003) — **LoRA adaptation ≫ linear probing**,
  and LoRA-adapted FMs reach ~100% of full performance with **10% of the training data**.
  This is the upgrade path once your pathologists have produced verified clicks: their
  verifications are exactly the LoRA training set.
* **Retrieval precedent**: SMILY (<https://arxiv.org/pdf/1901.11112>) is the original
  "select a region, get visually similar ones back" system for pathology; Yottixel-style WSI
  retrieval with UNI/Virchow/GigaPath embeddings reaches top-5 F1 0.40-0.42 vs 0.27 for
  DenseNet (<https://www.nature.com/articles/s41598-025-88545-9>) — FM embeddings roughly
  1.5x a generic CNN at retrieval.
* **SimpleShot-style prototypes**: with >1 click, average the embeddings into a class
  prototype and rank by cosine. Standard in few-shot detection with frozen backbones
  (DE-ViT, DeFRCN lineage); it is three lines of numpy over a cached matrix, and it is how
  your tool naturally accumulates evidence as the pathologist verifies.

**Licensing — check before you benchmark, this is a product:**

| model | license | commercial? | size |
|---|---|---|---|
| **Midnight-12k** (kaiko-ai) | **MIT** | **yes** | ViT-g/14, 224 px, 1536-d |
| **H-optimus-0** (Bioptimus) | **Apache-2.0** | **yes** | 1.1 B ViT, 224 px @ 0.5 µm/px, 1536-d |
| DINOv3 (Meta) | DINOv3 license | yes, with terms + gated signup | ViT-S…H+, ConvNeXt |
| DINOv2 (Meta) | Apache-2.0 | yes | ViT-S…g |
| UNI / UNI2-h | research-only | **no** | ViT-L / ViT-h |
| Virchow2 | CC-BY-NC-ND-4.0 | **no** (explicitly prohibited) | ViT-H |
| Phikon-v2 | Owkin custom (check LICENSE.pdf) | verify | ViT-L |

For AnnotateDx the shortlist is therefore **Midnight-12k (MIT)**, **H-optimus-0
(Apache-2.0)**, and **DINOv2/v3** as the natural-image control. UNI and Virchow2 are
benchmark-only; do not build on them.

### B2. Deep-feature template matching — the minimal-change version of what you asked for

This keeps your architecture and swaps only the correlation surface: compute a dense feature
map over the ROI, correlate the clicked patch's feature vector against it, keep your robust-z
floor, your NMS, your whole evaluation harness.

* **OS2D** (<https://arxiv.org/pdf/2003.06800>, code `aosokin/os2d`, ECCV 2020) — one-stage
  one-shot detection by dense correlation of *learned* features plus a feed-forward geometric
  alignment, trained end-to-end. The direct deep analogue of `template_match.py`.
* DINOv2/v3 dense patch tokens are the modern, training-free way to get the same
  correspondence surface.

Honest assessment: this is the lowest-friction path *for your existing code* — one function
changes — but it re-enters the problem through detection, which your own data says the click
loses at. Its real value is as a **cheap falsifier**: if cosine-in-feature-space also fails
to beat chromatin, you have learned something about the click, not about `cv2`.

### B3. Exemplar-conditioned counting/detection — read, don't build

CounTR, LOCA, **DAVE** (detect-then-verify, **0.3 s/image**, 1.4 s for >300 objects), **GeCo
/ GeCo2** (single-stage, 0.4 s, +27% MAE over LOCA), CoDi, OCCAM (training-free). These take
1-3 exemplar boxes and find "more like this" — literally your interaction.

Why it is a read and not a build: they are trained on FSC147-style natural images where the
target class is visually *obvious* and the challenge is density. Mitosis-vs-mimic is a
fine-grained distinction between objects of identical scale, colour and texture; these models
will very plausibly return every nucleus. If you probe one, probe GeCo (single-stage, code
released) on one ROI, and measure it on the AUC bar in §7 before believing any count.

### B4. Product precedent — the systems that already do this UX

* **AnnotateAnyCell** (bioRxiv 2025.11.02.686114) — the closest published analogue to
  AnnotateDx: segmentation → contrastive cell embeddings → interactive UMAP where
  morphologically similar cells cluster → active learning picks what the expert reviews.
  Open source.
* **A deep active learning framework for mitotic figure detection with minimal manual
  annotation** (Liu et al., *Histopathology* 2025, doi 10.1111/his.15506) — the same loop
  scoped to your exact task. (Paywalled; not read here.)
* **NuClick** (<https://arxiv.org/abs/2005.14511>) — one click inside a nucleus yields a
  precise mask. Use it for the click→object-boundary step so the template/prototype is the
  *cell*, not a fixed 51 px box. It also gives your seed patch a consistent scale, which your
  own logs flagged as a confound (template size varies 25-51 px with the seed, and the
  robust-z null width varies with it).
* **TissueWand**, **Quick Annotator**, **SMILY** — earlier interactive-annotation precedent.

---

## 6. Recommendation, in order

1. **Split the pipeline in the code** into `propose()` and `rank(click, candidates)`. Every
   result you already have is a statement about one or the other; right now they are welded.
2. **Run the falsification experiment in §7 before building anything.** One afternoon.
3. **Swap proposals to the MIT-licensed MIDOG++ FCOS** for the product path (with the leakage
   caveat for any measurement).
4. **Ship stage B as: frozen encoder → cached embedding matrix per ROI → cosine to click.**
   Start with Midnight-12k or H-optimus-0 embeddings computed on a rented GPU at ingest; keep
   a ResNet-18 @ 64 px path as the offline/laptop fallback (44 s for 20k, measured).
5. **Once verifications accumulate, LoRA-adapt the encoder** on them — the benchmark paper's
   10%-of-data result is the strongest evidence in this review that your human-in-the-loop
   data has compounding value.
6. **Use the click for what seedless methods structurally cannot do**: per-case adaptation to
   *this* slide's stain, *this* tumour's morphology, and this pathologist's threshold. That
   is the defensible product claim, and it is a ranking claim, not a detection one.

## 7. The experiment that decides all of this

`lookalike_auc_unbiased.py` already does almost everything. It scores candidates on the
`nucleus_blobs` set (96.6-100% capture of both classes) and computes Mann-Whitney AUC with
Hanley-McNeil CIs. Add one statistic — **cosine similarity to the clicked patch's
embedding** — and one metric.

**Gate on read-to-50%, not on AUC.** The two annotated classes are 105-123 look-alikes
against ~20,000 unannotated candidates per ROI, so a reader working down the ranked list
meets unannotated candidates in overwhelming proportion. Read-to-50% is therefore dominated
by the *unannotated* contrast (chromatin: 0.935 / 0.962), not the look-alike one (0.693 /
0.768). An encoder that merely ties chromatin on look-alikes but reaches 0.99 against
unannotated would win decisively on reading burden — and an AUC-only rule would reject it.
Read-to-50% is also the metric every prior log used and the one that maps to pathologist
time.

So report three numbers per (ROI, seed, encoder), all off the same blob candidate set:

| quantity | incumbent bar, 301.tiff | incumbent bar, 246.tiff |
|---|---:|---:|
| **read-to-50% on the blob candidate set** (the decision metric) | *compute for chromatin — not yet measured on this set* | *ditto* |
| AUC vs annotated look-alikes (diagnostic) | 0.693 [0.635, 0.751] | 0.768 [0.707, 0.830] |
| AUC vs unannotated (diagnostic) | 0.935 [0.912, 0.958] | 0.962 [0.937, 0.987] |

The two AUCs do not decide anything; they explain *why* read-to-50% moved.

**Two protocol requirements, both from mistakes this project has already made once:**

1. **Drop the seed annotation from the positive set for its own seed.** Cosine-to-click
   scores the clicked cell at 1.0 by construction; chromatin has no seed. Leaving it in is
   the same "selection correlated with the feature under test" that
   `lookalike_auc_unbiased.py`'s docstring was written to eliminate, reintroduced in a new
   form.
2. **Gate on the worst of 5 seeds, not the median.** The premise test found read-to-50%
   varying 2.2x with which cell was clicked, and unreached on 2 of 5 seeds. A ranker that
   wins on median and fails on one click in five is exactly the product defect already
   identified.

Encoders to score, cheapest first: ResNet-18/ImageNet (control), DINOv2 ViT-S, Midnight-12k,
H-optimus-0.

**Decision rule, set in advance:** cosine-to-click must beat chromatin's read-to-50% on the
*worst* of 5 seeds, on both decision-grade ROIs. If it does not, the click does not work as a
ranker either, and the honest product is a seedless detector where the click means
"confirm/reject", not "find more like this". If it does, everything downstream is
engineering.

## 8. What I would not do

* Do not tune the correlation search further. Six operating-point moves on one ranking
  function have all come back neutral, for the reason the chromatin log identified.
* Do not benchmark UNI or Virchow2 as product candidates — you cannot ship them.
* Do not run a ViT-g foundation model on this laptop; the estimate is 10-25 h per ROI.
* Do not measure detection quality on 001-506 with the released FCOS weights and report it.

## Sources

MIDOG detectors: <https://github.com/jonas-amme/FCOS_Inference_CLI> ·
<https://github.com/DeepMicroscopy/MIDOG25_T1_reference_docker> ·
<https://arxiv.org/abs/2208.12587> · <https://arxiv.org/abs/2508.19804> ·
<https://arxiv.org/abs/2509.02627> · <https://arxiv.org/pdf/2509.02599> ·
<https://arxiv.org/html/2509.02630> · <https://link.springer.com/chapter/10.1007/978-3-032-25180-0_1>
Foundation models: <https://arxiv.org/html/2607.28007v1> · <https://arxiv.org/abs/2508.04441> ·
<https://huggingface.co/kaiko-ai/midnight> · <https://huggingface.co/bioptimus/H-optimus-0> ·
<https://huggingface.co/paige-ai/Virchow2> · <https://ai.meta.com/blog/dinov3-self-supervised-vision-model/> ·
<https://www.nature.com/articles/s41598-025-88545-9> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC11403354/>
One-shot / exemplar: <https://arxiv.org/pdf/2003.06800> · <https://github.com/jerpelhan/GeCo> ·
<https://arxiv.org/pdf/2409.18686> · <https://arxiv.org/html/2601.13871v1>
Interactive annotation: <https://www.biorxiv.org/content/10.1101/2025.11.02.686114v2.full> ·
<https://arxiv.org/abs/2005.14511> · <https://arxiv.org/pdf/1901.11112> ·
<https://onlinelibrary.wiley.com/doi/10.1111/his.15506> · <https://pmc.ncbi.nlm.nih.gov/articles/PMC7518350/>
Nucleus segmentation: <https://arxiv.org/pdf/2408.01797> · <https://github.com/TIO-IKIM/CellViT>
Mimickers: <https://journals.sagepub.com/doi/10.1177/0300985820980049>
