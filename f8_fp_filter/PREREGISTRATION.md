# Pre-registration F8: does a learned FP filter on the candidate list buy reading depth?

Date: 2026-09-08. **Written before any F8 code is run.** Status: PLAN — awaiting approval.

F1–F7 are taken (F2, F6, F7 in flight). This is the next free number.

---

## 0. Three corrections to the request, made before anything is built

The experiment as proposed is sound in shape. Three things in it do not survive contact with
what is on disk, and each is recorded here rather than silently adjusted.

### 0a. There is no pretrained MIDOG++ encoder to use

`README.md`, verbatim: *"Due to space restrictions, we can't make available the weights of the
trained models."* A filesystem search for `*.pth *.pt *.ckpt *.h5` across the repo returns
nothing. `training.py` builds a fastai RetinaNet from scratch; the only artefact of it in this
fork is the config.

What is actually available locally: `~/.cache/torch/hub/checkpoints/resnet18-f37072fd.pth`
(ImageNet) and `dinov2_vits14_pretrain.pth`.

**Substitution, registered:** *frozen ImageNet ResNet-18 as the encoder, with a linear head
trained leave-one-domain-out on MIDOG++ labels.* Rationale, in the order that matters here:

* **No leakage.** An ImageNet encoder has never seen a mitotic figure. Any MIDOG++-trained
  encoder would have to be held out by domain anyway, and we would be training it ourselves.
* **Licence.** `annotatedx-product-goal` records that AnnotateDx is commercial, so
  non-commercial pathology weights (UNI, Virchow2) are benchmark-only. ImageNet ResNet-18 is
  BSD-licensed via torchvision. DINOv2 is Apache-2.0 and therefore also clean, but see below.
* **Cost.** `annotatedx-cpu-only-hardware` records ResNet-18 at 61 s / 20k 64 px crops on this
  machine. The seed-0 pools total 397k candidates → **~20 min** for the whole embedding pass.
  DINOv2 ViT-S/14 at the same scale is 3–5 h, so it is registered as a **subset consistency
  axis** (3 ROIs), never the primary.

The phrase "lightweight pretrained MIDOG++ encoder" is therefore delivered as
**pretrained-generic encoder + MIDOG++-trained linear head, held out by domain.** That is the
nearest honest thing to the request.

### 0b. The 95%-recall target and the 3:1 target are mutually exclusive, and the baseline proves it

Both were measured **before** the run, from the committed ledger (§2). They are two readings of
the same step function taken at opposite ends of it:

* To reach 95% recall the reader is at a dilution of **8:1 to 175:1** (median 56:1).
* To hold dilution at 3:1 the reader stops at **K₃ = 122 candidates** (median) having found
  **46%** of the mitoses.

No filter reconciles them, because `K₃ ≤ 4m` is arithmetic: precision ≥ 0.25 requires
TP(K) ≥ K/4 and TP(K) ≤ m. So the 3:1 list can never be longer than 4× the mitoses present, and
"list length at 3:1" is not the free variable the request treats it as.

**Consequence for the design.** The two stated targets become two separate registered
questions, not one operating point:

* **Q-A (the 95% question):** hold recall ≥ 0.95 per ROI, minimise **reading depth**.
* **Q-B (the 3:1 question):** hold dilution ≤ 3:1, maximise **recall@K₃**. The quantity to
  improve is `recall@K₃` (baseline median 0.46, ceiling 1.00), *not* `K₃` — `K₃` alone is
  gameable, an empty list has perfect precision. 245.tiff already scores `K₃ = 0`.

### 0c. "Hit at least 95% recall on every ROI" is already true, for free

`deep_recall = 1.0` in all 14 seed-0 cells of `results/tm_recall_workload_cells.csv` — every
mitosis is in the pool. The baseline attains 100% recall on every ROI by reading the whole list.
The 95% bar is not a feasibility bar for either arm; it is a **depth** bar. Registered as such.

---

## 1. The question

With the search, the pool and the candidate coordinates held byte-identical, **how much reading
depth does a learned true-positive/false-positive filter remove at fixed recall, and how much
recall does it add at fixed dilution — when its threshold is chosen without seeing the test
domain?**

The last clause is the whole experiment. An oracle threshold is reported as an upper bound; the
deployable number is the transferred one.

---

## 2. Pre-run baseline (measured, not predicted)

Seed 0, all 14 ROIs, derived from `results/tm_recall_workload_tp_ledger.csv` +
`results/tm_recall_workload_cells.csv`. `K95` = candidates read to reach ⌈0.95·m⌉ mitoses.
`K₃` = deepest depth at which precision ≥ 0.25 (definition fixed in §7). `rec@K₃` = recall there.

| ROI | domain | m | K95 | dilution at K95 | K₃ | rec@K₃ | 4m (K₃ ceiling) |
|---|---|---:|---:|---:|---:|---:|---:|
| 201.tiff | canine lung | 17 | 666 | 38:1 | 24 | 0.353 | 68 |
| 402.tiff | neuroendocrine | 104 | 1,103 | 10:1 | 331 | 0.798 | 416 |
| 013.tiff | breast | 17 | 1,173 | 68:1 | 23 | 0.353 | 68 |
| 233.tiff | canine lung | 17 | 1,329 | 77:1 | 20 | 0.412 | 68 |
| 094.tiff | breast | 81 | 1,631 | 20:1 | 161 | 0.519 | 324 |
| 529.tiff | melanoma | 19 | 1,740 | 91:1 | 10 | 0.211 | 76 |
| 246.tiff | lymphosarcoma | 115 | 1,913 | 16:1 | 242 | 0.635 | 460 |
| 548.tiff | melanoma | 238 | 2,076 | 8:1 | 800 | 0.874 | 952 |
| 403.tiff | neuroendocrine | 52 | 4,286 | 85:1 | 84 | 0.442 | 208 |
| 460.tiff | soft tissue sarcoma | 35 | 5,380 | 157:1 | 49 | 0.486 | 140 |
| 300.tiff | mast cell | 180 | 6,371 | 36:1 | 380 | 0.539 | 720 |
| 459.tiff | soft tissue sarcoma | 130 | 8,854 | 70:1 | 187 | 0.369 | 520 |
| 301.tiff | mast cell | 217 | 9,434 | 45:1 | 479 | 0.553 | 868 |
| 245.tiff | lymphosarcoma | 89 | 14,989 | 175:1 | **0** | **0.000** | 356 |
| **median** | | | **1,994** | **56:1** | **122** | **0.464** | |
| **total** | | 1,111 | **60,945** | | **2,790** | | 5,244 |

Two facts to carry forward. **245.tiff misses the 3:1 bar by one candidate.** Its first true
positive sits at depth **5**, so its best precision at any depth is 1/5 = **0.200** — just under
0.25, which is the whole reason `K₃ = 0` there. Deleting any *one* of the four false positives
above it would make `K₃ ≥ 4`. `K₃ = 0` is therefore a knife-edge, not a wall, and §2b measures
what the far side of it costs. And **the total 95%-recall workload across 14 ROIs is 60,945
candidates**; at 1 s per verification that is 17 h, at 3 s 51 h. That total is the number Q-A
moves.

First-TP depth per ROI, for reference — 2 px on seven ROIs, never worse than 5:

| depth of 1st TP | 2 | 3 | 5 |
|---|---|---|---|
| ROIs | 233, 246, 300, 403, 460, 459, 548 | 013, 094, 201, 301, 529 | 245, 402 |

*(This corrects a draft of this section which read the 45th TP's depth, 3,136, as the first's.
The error is recorded rather than overwritten: it would have made P4's mechanism claim — that no
threshold could rescue 245.tiff — arithmetically void, since the candidates above a first TP
contain zero mitoses by definition and deleting them costs no retention at all.)*

## 2b. The effect size Q-B requires, computed before the run

Model the filter as deleting a fraction `f` of false positives **independently of TM rank**. The
`j`-th mitosis then sits at depth `j + (1−f)(d_j − j)`, and it is inside the 3:1 list iff

```
f  >=  1 - 3j / (d_j - j)
```

Solving per ROI for the `f` that **doubles** `recall@K₃` (and, for 245.tiff whose baseline is 0,
that reaches recall 0.10):

| ROI | baseline rec@K₃ | target | `f` required | `f` required for rec@K₃ = 0.95 |
|---|---:|---:|---:|---:|
| 201.tiff | 0.353 | 0.706 | **0.31** | 0.92 |
| 459.tiff | 0.369 | 0.738 | 0.76 | 0.96 |
| 529.tiff | 0.211 | 0.421 | 0.77 | 0.97 |
| 233.tiff | 0.412 | 0.824 | 0.87 | 0.96 |
| 245.tiff | 0.000 | 0.101 | 0.89 | 0.98 |
| 403.tiff | 0.442 | 0.885 | 0.90 | 0.97 |
| 013.tiff | 0.353 | 0.706 | 0.90 | 0.96 |
| 548.tiff | 0.874 | 1.000 | 0.91 | 0.63 |
| 094.tiff | 0.519 | 1.000 | 0.94 | 0.85 |
| 402.tiff | 0.798 | 1.000 | 0.94 | 0.70 |
| 246.tiff | 0.635 | 1.000 | 0.94 | 0.82 |
| 300.tiff | 0.539 | 1.000 | 0.95 | 0.92 |
| 301.tiff | 0.553 | 1.000 | 0.96 | 0.93 |
| 460.tiff | 0.486 | 0.971 | 0.98 | 0.98 |
| **median** | 0.464 | | **0.905** | **0.945** |

**Doubling `recall@K₃` requires deleting a median 90.5% of false positives while keeping
essentially every mitosis.** That is the bar, stated before the encoder is loaded. It is not a
bound — it is the requirement *under rank-independence*, and the deviation from it is itself the
measurement:

* If the head is **complementary** to `tm_score` — it deletes FPs that TM ranks highly — the
  observed `f` needed is **lower** than the table.
* If the head is **correlated** with `tm_score` — it deletes what TM already buried — the observed
  `f` is **higher**, and the filter is doing work the ranker had already done.

F8 reports observed `f` against required `f` per ROI, which separates those two cases directly.
201.tiff at `f = 0.31` is the one cell where the bar is plausibly reachable.

---

## 3. Sample and what is frozen

**Sample.** 14 ROIs of `images/extra_valid` (`MANIFEST.md`: 2 per tumour type, all 7 domains,
`n_mitotic ≥ 15`, seed pool ≥ 5), **seed index 0 only**, as requested. 14 cells.

**Frozen by construction, not by re-declaration.** No search is re-run. F8 reads
`results/tm_recall_workload_pool_z.npz`, whose arrays are score-descending so array order *is*
rank. That bakes in: padded ROI + `BORDER_REPLICATE`, `hematoxylin_od`, `cv2.TM_CCOEFF`,
largest-CC template tightening, `PEAK_MIN_DISTANCE = 7`, NMS radius 5.0 µm (the default),
5.0 px self-hit filter, `DEEP_FLOOR_Z = -1.5`, and the seed draw.

**The candidate list under test** is that pool thresholded at `CURRENT_Z = 1.0`, the operating
point `recall_workload_ledger.py` measures every saving against: 13,208–36,284 candidates per
ROI. The cutoff is **not a knob** — every depth reported here (K₃ ≤ 800, K95 ≤ 14,989) is
shallower than the shallowest pool, and high_z's prefix identity makes recall@K exactly
cutoff-invariant below `z_max`. Stated once, then fixed.

**Seed 0 is a real limitation.** D4's rule is worst-of-5 clicks, and `annotatedx-product-goal`
records per-seed variance as a product defect. A 1-seed result cannot speak to click stability.
Seeds 1–4 are already cached coordinates, so extending costs only the embedding pass (~80 min);
registered as a **conditional follow-up if the primary is positive**, not part of F8.

---

## 4. Gates, in order. Any failure stops the run.

**G0 — environment.** `torch.from_numpy` currently raises *"Numpy is not available"* on the
system python3 (numpy 2.4.6 against torch 2.2.2 built for numpy 1.x — the same pin
`requirements-midog-utils.txt` documents). A `numpy<2` venv is built and smoke-tested
(`torch.from_numpy(np.zeros((2,2),np.float32))`) **before any pipeline code is written**.

**G1 — reproduction.** Re-run `ev.bucket_detections` on the unpruned cached seed-0 pool and
assert the resulting TP ranks equal `tm_recall_workload_tp_ledger.csv`'s `rank` column
element-wise, per ROI. `match_radius_px` is read from the cells CSV, never recomputed from
`mpp`. The npz arrays are used in stored order and never re-sorted. `gt_eval` is rebuilt as
`ds.image_annotations(annotations, fn)` minus the row whose `ann_id` is that cell's
`seed_ann_id`. If this fails anywhere, every downstream number is void.

**G2 — instrument check, before any arm is scored.** F7 §5a is the precedent: its statistic was
disqualified pre-run for measuring the window rather than the object. The analogue here is that
the embedding encodes **scanner and stain identity** rather than mitosis identity — at ~0.3%
TP prevalence, a nuisance axis with unit variance buries a between-class term of order
p(1−p)d². Two probes are trained on the *same* frozen embeddings:

| probe | target | protocol |
|---|---|---|
| domain-ID | 7-way tumour type | **held out by ROI** — `GroupKFold` on `file_name`, so the probe is never scored on an ROI it was fitted on |
| TP/FP | binary | **the same 7 LODO folds the head uses**, so its number is directly comparable to §8's |

**Both probes are scored out-of-fold only.** A 512-d linear probe fitted and scored in-sample on
14 ROIs reaches ≈1.0 on either target by construction, which would make this gate unable to fire —
and firing early is its entire purpose. In-sample numbers are recorded beside the held-out ones
solely to show the size of that gap.

If domain-ID is near-perfect out-of-fold while TP/FP is near-chance, LODO transfer will fail and
F8 stops at §4 with that as its result. This costs minutes and is run before the arms.

---

## 5. The arms

Three. All three consume the identical candidate set; they differ only in the order (and
membership) of the list handed to the reader.

| arm | definition |
|---|---|
| **A — baseline** | rank by `tm_score` (i.e. `z`) descending. The shipped ranker, D5. |
| **B — filter then re-rank** | delete candidates whose head probability < τ; rank the survivors by `tm_score` descending. The requested variant. |
| **C — rank by probability** | no deletion; rank the full list by head probability descending. |

**Why C is in.** It is free once the embeddings exist and it bounds what B leaves on the table.
Neither dominates: B can win when the head has a good high-specificity point but a poor ordering,
C can win when the reverse holds. One extra arm, not a sweep.

**Encoder.** Frozen `torchvision.models.resnet18(weights=IMAGENET1K_V1)`, global-average-pooled
512-d penultimate features, `eval()`, no fine-tuning.

**Crop.** RGB, ImageNet normalisation. **Physically sized, not pixel sized**: `round(14.5 / mpp)`
px centred on `(cx, cy)`, resized to 64×64. `mpp` spans 0.226–0.253 across these ROIs, so a fixed
pixel crop is a 12% different field of view per scanner — exactly the nuisance axis G2 tests for.
The fixed-64-px variant is recorded alongside as a consistency axis, never as primary.
Candidates within the crop half-width of the ROI edge are read through
`tm.read_padded_patch`; if that returns `None` the candidate is **kept and never filtered**
(a border rule would be a second knob). Their count is reported per ROI.

**Head.** `sklearn.linear_model.LogisticRegression`, `class_weight='balanced'`, L2, on
standardised features. One knob (`C`) chosen inside the training fold by grouped CV over
training ROIs; never on the test domain. No MLP, no calibration stack — `prefers-direct-
measurement-over-method-machinery` applies, and a linear probe is the direct measurement of
"do these features separate the classes".

---

## 6. Labels — two definitions, deliberately

**Training labels are distance-based** and therefore order-independent: a candidate is positive
iff it lies within `match_radius_px` of a **category-1** annotation of that ROI, excluding the
seed's own `ann_id`. Greedy-match labels would make the training set depend on the ranking being
tested.

**Evaluation labels are greedy-match**, re-run **in each arm's own order**. This is not optional:
`Research Logs/2026-09-08-tp-fp-separability-audit.md` §2 (`C6`) found that scoring re-ranked
lists with TM-order labels biased `read_95` by up to 23% *against* the re-ranker, and that
`results/tp_fp_reading_depth.csv` carries the bias into a decision D5 cites. Arm A's order is the
matching order, so A is unaffected; B and C are re-matched.

**Category-2 neighbours are negatives.** Per that audit's `C1`, ~52% of the previously set-aside
"on annotation" population is within radius of a **look-alike** — a structure a pathologist
examined and *rejected*. They stay in the negative class, and `n_fp_on_lookalike` /
`n_fp_on_mitosis` / `n_fp_unannotated` are reported as separate columns so the composition is
visible rather than assumed.

---

## 7. Metrics — exact definitions, fixed now

Let the arm's list be `d₁ … d_n`, `TP(K)` the greedy-matched mitoses in the first `K`, `m` the
ROI's mitosis count (seed excluded).

* **`K95` — reading depth at 95% recall.** The smallest `K` with `TP(K) ≥ ⌈0.95·m⌉`, paired with
  a **`reached_95` boolean**. When an arm never reaches 95% — possible for arm B, where deletion
  can put it out of reach — the cell is **censored**: `K95` is null, `reached_95` is False, and
  the ROI is **excluded from every median**, with `n_censored` printed beside that median. It is
  *not* entered at `n`. Entering a sentinel into a median is the failure the
  `tm-recall-workload` audit found in that notebook's `*_median` columns, and an arm that fails
  on 2 ROIs must not be able to look better by having those 2 counted at full list length.
  A median over a different denominator than the arm it is compared to is reported as such.
* **`K₃` — the 3:1 list length.** The **largest** `K` with `TP(K)/K ≥ 0.25`. Last crossing, not
  first: precision along a ranked list is non-monotone, and the first crossing is decided by
  whether candidate #1 happens to be an FP (402.tiff: 331 vs 0; 403.tiff: 84 vs 9 — the same
  ROI under the two rules). Last crossing is the plain reading of *"3 FPs for every TP in the
  list"*. `K₃ = 0` is a value, not missing data.
* **`recall@K₃`** — always reported *with* `K₃`. This is the arbiter for Q-B.
* **`recall@K` at fixed budget** — primary `K = 250` (the repo's standing budget, D4/F1/F5);
  grid `K ∈ {50, 100, 250, 1000}` reported.
* **`coverage_frac`** at every operating point, per `evaluate.py`'s own rule that no un-budgeted
  recall number may be read without it. Filtering changes the retained geometry, so arm B's
  coverage is a different number from the pool's and must be computed, not inherited.
* **Retention** — `n_kept / n_candidates` and `mitoses_kept / m` per ROI for arm B.

**Pooled AUC is not the arbiter.** D5's reasoning about `od51` — best at `read_50`, worst tail,
therefore wrong for this product — applies unchanged. AUC is reported as a diagnostic; the depth
metrics decide.

---

## 8. Threshold selection: leave-one-domain-out, and an oracle for contrast

Seven folds. For each held-out domain: train the head on the 12 ROIs of the other six domains,
choose τ on those training ROIs as the **largest** threshold retaining ≥ 95% of *their* mitoses in
the worst training ROI, then apply that τ unchanged to the held-out domain's 2 ROIs. **Largest,
not smallest**: retention is monotone decreasing in τ, so the set of thresholds meeting the
retention bar is `[0, τ*]` and its smallest element is τ = 0 — a filter that deletes nothing and
turns arm B into arm A.

Reported as a pair, always:

| | how τ is chosen | what it means |
|---|---|---|
| **oracle** | per test ROI, the τ that retains exactly 95% of *its own* mitoses | upper bound; not deployable |
| **transferred** | LODO, as above | the deployable number, and the primary result |

**A transferred-τ failure is a result, not a bug.** The repo's LODO history is 3 failures in 35
rows (`tm_recall_workload_z_sweep`), and D5's own transfer table shows the outside cutoff sitting
above a domain's own on two domains. If arm B drops below 95% retention on some ROI under the
transferred τ, that is F8's honest primary finding about deployability.

**Training set choice — the fork, and the call made.** The alternative was to train on the 9
downloaded ROIs *outside* `extra_valid` (001, 002, 202, 350, 351, 405, 406, 505, 506) and keep
all 14 as a pure test set. LODO within the 14 is chosen instead: it measures the transfer risk
directly, uses ROIs whose seed pools and mitosis counts are known decision-grade (the 9 are not
all — 202.tiff explicitly fails the seed bar), and needs no new search runs. The cost is stated
in §11.

---

## 9. Registered predictions

Signed before any number is seen. Each is falsifiable against a table F8 will produce.

* **P1 (instrument).** The domain-ID probe reaches ≥ 0.95 pooled accuracy — the embedding is
  substantially scanner/stain — while the pooled TP/FP probe still reaches ROC-AUC ≥ 0.85.
  *Falsified if TP/FP AUC < 0.75, in which case F8 stops at G2.*
* **P2 (LODO transfer).** LODO TP/FP ROC-AUC has median ≥ 0.85 over the 7 folds but at least one
  fold falls below 0.75. In-domain handcrafted chromatin reached 0.92–0.997
  (`tp_fp_separability`); a LODO linear probe on a generic encoder should sit below that floor
  on its worst domain.
* **P3 (Q-A, the 95% question).** Arm B's median `K95` falls to **≤ 1,200** from the baseline
  1,994 (a ≥ 40% cut) under the *oracle* τ, but under the *transferred* τ **at least 2 of the 14
  ROIs are censored** (`reached_95` False), and the medians are reported over the surviving
  denominator per §7. 245.tiff and 459.tiff stay above `K95 = 4,000` in every arm.
* **P4 (Q-B, the 3:1 question).** The registered quantity is §2b's **observed `f` against required
  `f`**. Prediction: at a τ retaining ≥ 95% of mitoses, the head's observed FP-deletion fraction
  is **`f` ≤ 0.6 on at least 10 of the 14 ROIs** — below §2b's median requirement of 0.905 — so
  **median `recall@K₃` improves by less than the doubling §2b prices, landing in 0.50–0.65**
  against a baseline of 0.464. **201.tiff (`f` required = 0.31) is the one ROI predicted to
  double.** *Falsified if median `recall@K₃` ≥ 0.75, or if observed `f` exceeds required `f` on a
  majority of ROIs.*
* **P4b (245.tiff, the knife edge).** `K₃` goes **from 0 to non-zero** in arm B — it needs only
  one of the four FPs above its depth-5 first TP deleted — but `recall@K₃` stays **below 0.10**,
  because reaching its 5th mitosis (baseline depth 83) needs `f ≥ 0.89`. *Falsified if
  245.tiff's `recall@K₃` ≥ 0.10 at a τ retaining ≥ 95%.*
* **P5 (B vs C).** Arm C beats arm B on median `recall@250`, and arm B beats arm C on median
  `K95`. If one arm wins both, the registered reason for including C was wrong and that is
  recorded.
* **P6 (the null).** Explicitly allowed: this is the sixth re-ranking attempt in this repo after
  five nulls (F7 §2). A result of "no arm beats `tm_score` on the depth metrics at a transferable
  threshold" is a complete F8 and will be written up as one.

---

## 10. Time translation

Reader verification time is a **stated assumption, not a measurement** — no timing study exists
in this repo. Reported at **1 s and 3 s per candidate**, per ROI and summed over the 14:

* Baseline Q-A total: 60,945 candidates = **16.9 h @ 1 s**, 50.8 h @ 3 s.
* Baseline Q-B total: 2,790 candidates = **46 min @ 1 s**, 2.3 h @ 3 s.

The deliverable is Δ against these, per arm, with the assumption restated in the same sentence
as every hour figure. No single "time saved" number is quoted without the recall it was bought
at — Q-A's hours and Q-B's minutes are not the same reader doing the same job.

---

## 11. Outputs

```
f8_fp_filter/f8_embed.py          encoder pass → .cache_tail/f8_embeddings_seed0.npz   (~20 min)
f8_fp_filter/f8_filter_eval.py    G1/G2 gates, 7 LODO folds, 3 arms, all metrics
f8_fp_filter/f8_fp_filter.ipynb   tables + figures, reading only the CSVs below
results/f8_fp_filter_cells.csv        one row per (ROI, arm, τ-source): every §7 metric
results/f8_fp_filter_lodo.csv         one row per fold: AUC, τ, retention, failures
results/f8_fp_filter_probes.csv       G2: domain-ID vs TP/FP
results/f8_fp_filter_verification.csv G1 reproduction, per ROI
```

Nothing overwrites an F1–F7 artefact. The embedding cache is a pure cache — delete to regenerate.

---

## 12. Limitations, named before the run

1. **One seed.** §3. Cannot speak to click stability, which is a product defect metric here.
2. **Domain-transfer and ROI-transfer are confounded.** Each domain is 2 ROIs, so holding out a
   domain also holds out two specific ROIs. The same limitation the
   `tm-recall-workload` audit (F8 of that log) names and the repo already accepts.
3. **The encoder is generic, not pathology-pretrained.** A domain-appropriate encoder would
   likely do better; this measures the floor, not the ceiling, of the idea.
4. **`K₃` is a composition metric, not a stopping rule a reader could execute.** It requires
   knowing the labels. It answers "how long is the useful part of this list", not "when should
   the reader stop".
5. **1,111 mitoses over 14 ROIs, with 4 ROIs carrying 66% of them** (548, 301, 300, 246). Any
   pooled statistic is a statement about those four; every ROI-level test runs on ROI medians.

---

## 13. What would make this abandoned mid-run

* G1 fails → the cached coordinates do not reproduce the ledger; nothing downstream is valid.
* G2's TP/FP probe is at chance while domain-ID is perfect → the embedding is a scanner
  detector; report that and stop (P1 falsified).
* No τ exists that retains 95% of mitoses on the *training* domains without keeping >90% of the
  candidates → the filter has no operating point worth transferring; report the retention curve
  and stop.
