# What 100% sensitivity actually costs: the tail statistic

Date: 2026-09-01
Question asked: *"can PyHIST segment the cells, then a classifier pick the mitotic ones?
Use one mitotic cell as a proxy for the pathologist's click and find all the rest. We
definitely want to hit all the mitotic cells; some false positives are fine, but not
'flag every nucleus'."*
Code: `tail_sensitivity.py`, `tail_object_audit.py`, `miss_attribution.py`,
`fp_filter_probe.py`, `fp_filter_domain.py`, `watershed_split.py`
Results: `results/tail_sensitivity.csv`, `results/tail_object_audit.csv`,
`results/miss_attribution.csv`, `results/fp_filter_domain.csv`, `results/watershed_tail.csv`
Follows: `2026-09-01-one-click-retrieval-literature.md`, `2026-08-31-premise-test-results.md`

---

## Summary

1. **PyHIST is the wrong tool** — it is tissue-vs-glass tiling, which `tissue_mask` already
   does in one line. It has no instance segmentation. (§0)
2. **Every metric this project has used is a median; "hit all the mitoses" is a tail
   statistic.** They disagree, and the tail had never been measured. (§1)
3. As the detector stands, **100% sensitivity costs 24–97% of the candidate pool**, and two
   of seven ROIs cannot reach it at all. (§2)
4. **The cause is one defect, not two**: connected-component labelling fuses touching nuclei,
   so the merged object's centroid misses the annotation and its mean intensity is diluted.
   All eleven problem annotations sit in components that are among the *darkest* in the ROI —
   they are under-segmented, not faint. (§4)
5. **Splitting them with a watershed fixes the ceiling on 7 of 7 ROIs** and cuts the reading
   depth for 100% sensitivity by 6.6x on 246.tiff — while regressing 301.tiff and 245.tiff,
   which were already the two hardest. One untuned knob explains the regression. (§6)
6. A **22-feature logistic regression** false-positive filter rescues the single worst ROI
   (2.0x) and is a coin flip elsewhere — worth shipping, not an average-case win. (§5)
7. **Nothing here used the click.** Its role should be per-slide threshold-setting and
   online learning, not template matching — and that is still unmeasured. (§8, Step 3)

---

## 0. PyHIST: no, and you already have its useful half

PyHIST is a **tissue-vs-glass tile extractor** — Felzenszwalb graph segmentation on a
downsampled WSI, then it writes out tiles that contain tissue. It has no notion of a cell
or a nucleus; there is no instance segmentation in it at all. Its entire output is
functionally `baselines.tissue_mask`, which is already in this repo, is one `cv2` call, and
was validated at gray < 220 excluding **0 of 690** mitotic annotations across 14 ROIs.

So PyHIST is not a route to cell segmentation. `nucleus_blobs` (tiled Otsu → connected
components) is the cell-segmentation stage in this repo, and it runs in 7–11 s per ROI.

The instinct behind the question turns out to be **half right, and right about the half
this project had written off.** Going in, segmentation looked solved (96.6–100% capture) and
ranking looked like the whole problem. It is the other way round: §4 shows the segmentation
silently fuses touching nuclei, and §6 shows that fixing it is worth more than every ranking
change measured here put together. What PyHIST cannot supply — instance segmentation — is
exactly what was missing.

## 1. The reframe: every metric used so far is a median; the requirement is a tail

`recall@budget` and `read-to-50%` — the two metrics carrying every conclusion in this
project — are **median** statistics. "We definitely want to hit all the mitotic cells" is a
**tail** statistic: the reading depth at which the *worst-ranked* true mitosis appears.

These are not the same question and they do not have the same answer. AUC 0.986
(`blob_component_score` vs unannotated on 301.tiff) is compatible with one mitosis sitting
at rank 14,000 of 20,000 — which makes 100% sensitivity cost the whole slide while leaving
the AUC untouched. **Nobody had measured the tail.** `tail_sensitivity.py` measures it.

Protocol: the seedless `nucleus_blobs` candidate set (so the proposal ceiling is identical
across arms and only the *ordering* is under test), on the **7 ROIs with ≥ 14 annotated
mitoses** — which adds 300.tiff (181 mitoses) and 245.tiff (90), the second- and
fourth-densest available and both absent from the previous 7-ROI set. Sensitivity is
counted over the mitoses the proposal stage *captured*; the ones it missed are unreachable
at any depth and are reported separately as the ceiling.

## 2. With the detector as it stands, 100% sensitivity costs 24–97% of the pool

*(§6 changes this verdict. Everything in this section describes the current
`nucleus_blobs` connected-component detector, which §4 shows is defective.)*

Candidates a reader must work through, best ranker per ROI:

| ROI | domain | pool | mitoses | ceiling | **95%** | **98%** | **100%** |
|---|---|---:|---:|---:|---:|---:|---:|
| 301 | mast cell | 20,370 | 218 | 1.000 | **1,428** (7.0%) | 4,307 (21%) | 7,750 (**38%**) |
| 300 | mast cell | 17,376 | 181 | 1.000 | **1,157** (6.7%) | 2,510 (14%) | 4,246 (**24%**) |
| 246 | lymphosarcoma | 21,059 | 116 | **0.966** | **1,757** (8.3%) | 1,966 (9.3%) | 20,493 (**97%**) |
| 245 | lymphosarcoma | 34,910 | 90 | 1.000 | 13,378 (38%) | 15,959 (46%) | 17,473 (**50%**) |
| 201 | lung | 12,446 | 18 | **0.889** | **823** (6.6%) | 823 | 823 (6.6%) |
| 202 | lung | 13,440 | 16 | 1.000 | **2,304** (17%) | 2,304 | 2,304 (17%) |
| 405 | soft tissue | 7,577 | 14 | 1.000 | **1,102** (15%) | 1,102 | 1,102 (15%) |

**Reading 24–97% of every nucleus in the slide is the "flag every nucleus" outcome the
product exists to avoid.** At the same time, **95% costs only 6.7–8.3% of the pool on
301/300/246 — a 12–15x reduction, at 6–15 false positives per true positive.** The gap
between those two rows is the whole finding.

### How thin the evidence for that 7% is

Two qualifications, both visible in the table above and both load-bearing:

* **245.tiff needs 38% of the pool at 95%**, not 7%. It is the same tumour and scanner as
  246.tiff, so this is not a domain effect, and it is unexplained. The affordable operating
  point is *not* a single number.
* **The 95% column does not exist on 201, 202 and 405.** With 14–18 mitoses,
  ⌈0.95 × n⌉ = n, so their 95/98/100% columns are the same cell by construction and carry
  no information about the tail's shape.

So "95% is affordable" rests on **301, 300 and 246 — three ROIs, of which 301 and 300 share
tumour type and scanner. Effectively two domains.** Adding 300.tiff and 245.tiff did widen
the base the 2026-09-01 audit called too thin, and the base is still too thin. Both are
true.

### The cliff is one or two annotations

246.tiff: 98% of captured costs **1,966** candidates; 99% costs **17,321**. A single
mitosis sitting at rank ~17k multiplies the reading burden by **8.8x**. This is the defining
shape of the problem — the tail is set by a handful of individual objects, not by the bulk
of the distribution, so no amount of improving the *average* ranker moves it.

### The false positives are ordinary nuclei, not the pathologist's look-alikes

At 95% on 301.tiff the 1,428 candidates read contain **91 annotated look-alikes and 1,129
unannotated nuclei**. Across every ROI and level, look-alikes are 1–8% of the false-positive
mass. The "hard negative / mimicker" framing that has driven the texture and
chromatin work targets the small share; **the volume problem is ordinary nuclei**, and a
ranker that only separates mitoses from annotated mimics will not move the reading burden.

## 3. The committed chromatin ranking reverses at the operating point that matters

*(Measured on the pre-split connected-component candidate set. §6 changes that candidate
set, so this comparison — like the §5 filter — is owed a re-run on the split candidates.)*

Depth to 95% of captured, `blob_score` (the Otsu component's own mean hematoxylin — the
statistic `nucleus_blobs` already sorts by) against `chromatin_od` (the committed re-ranker):

| ROI | mitoses | `blob_score` | `chromatin_od` | winner |
|---|---:|---:|---:|---|
| 301 | 218 | **1,428** | 8,106 | blob, 5.7x |
| 300 | 181 | **1,157** | 5,096 | blob, 4.4x |
| 246 | 116 | **1,757** | 3,629 | blob, 2.1x |
| 245 | 90 | **13,481** | 15,458 | blob, 1.1x |
| 201 | 18 | 4,794 | **823** | chromatin, 5.8x |
| 202 | 16 | **2,304** | 5,409 | blob, 2.3x |
| 405 | 14 | 3,630 | **1,102** | chromatin, 3.3x |

**`blob_score` wins on all four ROIs with ≥ 90 mitoses; `chromatin_od` wins on two of the
three with ≤ 18.** Commit `7c3af93`'s 19.4x chromatin advantage is not contradicted — it is
a read-to-50% result, and read-to-50% is a median. The two statistics rank differently in
the *tail*, and the tail is what the product requirement is about. A rank-fusion of the two
lands between them everywhere and wins nowhere, consistent with `chromatin.rerank`'s own
docstring note that fusion never beat `od` alone on the median metric either.

**This is not a retraction of `7c3af93` and must not be read as one.** The 19.4x was
read-to-50% on the *seeded pipeline's* candidate set; this is depth-to-95% on the
*seedless blob* set. Different candidate set, different metric, different operating point —
neither result refutes the other, and the reversal has not been measured on the seeded set
at all. The scoped claim is: *on the blob candidate set, at 95% sensitivity, `blob_score`
ranks better on the ROIs with the most annotations.*

## 4. The ceiling and the tail are **one** defect: connected components fuse touching nuclei

Two ROIs cannot reach 100% at any depth — 246.tiff misses 4 of 116, 201.tiff misses 2 of 18.
A first attempt at attribution matched each miss to the **nearest component centroid** and
reported areas of 1–17 px, which read as "faint mitoses an intensity threshold cannot see."
**That attribution was wrong**, and the object-level audit caught it: those same four
annotations on 246.tiff sit in the **97th percentile of chromatin density** for the ROI.
They are among the darkest objects present. Nearest-centroid finds a stray speck because a
large fused component's centroid lies far from any nucleus inside it.

Asking the right question — *which connected component contains the annotated pixel*
(`miss_attribution.py`, `results/miss_attribution.csv`) — gives a single answer:

| ann | ROI | status | component at the click | chromatin pct | why |
|---|---|---|---:|---:|---|
| 6438 | 246 | missed | 1,173 px | 0.975 | passes the gate; **centroid > 30.2 px away** |
| 6446 | 246 | missed | 2,462 px | 0.975 | same |
| 6573 | 246 | missed | 1,734 px | 0.972 | same |
| 6603 | 246 | missed | 8,070 px | 0.938 | over `max_area=4000` |
| 4466 | 201 | missed | 4,520 px | 0.967 | over `max_area` |
| 4479 | 201 | missed | 5,156 px | 0.972 | over `max_area` |
| 6419 | 246 | matched at 28.2 px | 1,397 px | 0.971 | edge of the 30.2 px radius |
| 6609 | 246 | matched at 26.7 px | 2,148 px | 0.979 | edge of radius |
| **6626** | 246 | matched at 26.8 px, **rank 17,320** | 1,428 px | 0.848 | edge of radius |
| **6634** | 246 | matched at 27.7 px, **rank 20,492** | 1,461 px | 0.711 | edge of radius |
| 4491 | 201 | matched at 28.4 px | 3,408 px | 0.972 | edge of radius |

**Every one of the eleven problem annotations sits in a chromatin-dense component that
exists.** Three are rejected by the upper area gate; the rest are matched — or missed —
because `label(connectivity=2)` fused two or more touching nuclei into one component, so
its centroid lands *between* them and its mean hematoxylin is diluted by whatever it
merged with.

That last clause is the important one, because it explains the cliff as well as the
ceiling. 246.tiff's two deep-tail annotations, at ranks 17,320 and 20,492, are **both**
edge-of-radius matches (26.8 and 27.7 px against a 30.2 px radius) to fused components
whose diluted mean intensity is what put them in the bottom quintile. They are not badly
ranked mitoses. They are under-segmented ones.

**So the ceiling problem and the tail cliff are the same defect, and it is under-segmentation
— not faintness, not the ranker, and not the area filter except in three cases.** The
standard fix is a distance-transform watershed, which is ~15 lines of `skimage`. Whether it
works is measured in `watershed_split.py`.

## 5. The cheap supervised filter: rescues the worst ROI, a coin flip elsewhere

This is the "keep `nucleus_blobs`, then add a step that filters out false positives" idea,
measured. `fp_filter_probe.py` extracts **22 features per candidate** — the region
properties the component already has (area, perimeter, eccentricity, solidity, extent, axis
lengths, Euler number, Feret diameter, intensity mean/max/min, plus aspect, circularity and
intensity range) and window statistics on unclipped hematoxylin OD at 31/51/81 px, a 121 px
context value, and two contrasts. No embeddings, no CNN, no GPU. Training is seconds;
scoring 20k candidates is milliseconds.

Two protocols, because they disagree: **leave-one-ROI-out** (`loro`) and **leave-one-domain-out**
(`lodo`). The seven ROIs are domain pairs sharing tumour type *and* scanner, so `loro`
leaves 300.tiff in training when scoring 301.tiff. Only `lodo` is the deployment condition.
(405.tiff is the only soft-tissue ROI, so its `loro` and `lodo` columns are identical — a
free check that the split is doing what it claims, and they are.)

Depth to 95% of captured, leave-one-domain-out:

| ROI | `blob_score` | `chromatin_od51` | **`logreg_lodo`** | `hgb_lodo` | vs `blob_score` |
|---|---:|---:|---:|---:|---:|
| 245 | 13,481 | 15,458 | **6,785** | 9,358 | **0.50x** |
| 405 | 3,630 | 1,102 | **1,454** | 2,063 | **0.40x** |
| 202 | 2,304 | 5,409 | **1,198** | 2,952 | **0.52x** |
| 300 | 1,157 | 5,096 | **804** | 5,212 | **0.70x** |
| 201 | 4,794 | **823** | 5,016 | 1,260 | 1.05x |
| 301 | **1,428** | 8,106 | 1,738 | 4,965 | 1.22x |
| 246 | **1,757** | 3,629 | 2,555 | 5,471 | 1.45x |

The median ratio over all seven ROIs is 0.70x, but **that number should not be quoted.**
Three of the seven (201, 202, 405) have 14-18 mitoses, so their 95% column *is* their 100%
column and measures the single worst object rather than a quantile. Restricted to the four
ROIs where 95% resolves as a quantile:

| ROI | ratio | |
|---|---:|---|
| 245 | **0.50x** | win |
| 300 | **0.70x** | win |
| 301 | 1.22x | loss |
| 246 | 1.45x | loss |

**Two wins, two losses, and no median improvement.** The apparent 1.4x comes entirely from
the three low-count ROIs. The defensible claim is narrower and still useful: *the supervised
filter's one clear win is 245.tiff at 2.0x — the ROI on which every hand-made statistic is
worst — and elsewhere it is a coin flip.* It beats `chromatin_od51` on 5 of 7. It is worth
shipping because it is free at inference and it rescues the worst case, **not** because it
lifts the average.

Three things worth recording:

* **Logistic regression beats gradient boosting on 6 of 7 under `lodo`**, often heavily
  (301: 1,738 vs 4,965; 300: 804 vs 5,212). The boosted model is better within a domain and
  worse across one — it memorises the training domains' positives. The *simpler* model is
  the right choice here, which is convenient.
* **The domain leak is real but modest**: `loro` → `lodo` costs 0-20% on five ROIs and 49%
  on 246.tiff. Any figure quoted from a leave-one-ROI-out protocol on this dataset is
  optimistic by roughly a fifth.
* **It does not fix the tail.** On 246.tiff, `logreg_lodo` still needs 13,804 candidates to
  reach 99%. The under-segmentation defect in §4 is invisible to any re-ranker, because the
  object it should rank was never proposed as a separate candidate.

`fp_filter_probe.py` initially fixed the true-positive labels from a single greedy match in
`blob_score` order and reused them across arms, which silently penalised every other arm —
it reported `chromatin_od51` needing 12,124 candidates on 201.tiff where re-matching gives
823. `fp_filter_domain.py` re-matches per arm, and its `chromatin_od51` column now
reproduces `tail_sensitivity.csv` exactly on all seven ROIs. **Use
`results/fp_filter_domain.csv`, not `results/fp_filter_tail.csv`.**

## 6. Splitting merged nuclei fixes the ceiling on every ROI, and the cliff with it

`watershed_split.py` implements the §4 fix: connected components, then any component large
enough to hold two nuclei (≥ 2 × `min_area`) is re-segmented by a distance-transform
watershed **inside its own bounding box**. Restricting the work to plausibly-merged
foreground objects — rather than running the transform over the whole 5412 × 7215 ROI — took
the run from >10 min per ROI to 83–346 s, and leaves genuinely single nuclei untouched.

### The ceiling: 100% capture on all seven ROIs

| ROI | connected components | + `max_area` 10,000 only | **+ watershed** |
|---|---:|---:|---:|
| 246 | 0.9655 (112/116) | 0.9655 | **1.0000 (116/116)** |
| 201 | 0.8889 (16/18) | 0.9444 (17/18) | **1.0000 (18/18)** |
| other five | 1.0000 | 1.0000 | 1.0000 |

**Raising the upper area gate is not the fix.** On its own it recovers one annotation on
201.tiff and none on 246.tiff, exactly as §4's attribution predicts — three of the six misses
are over-large components, and the other three are fused components of *legal* size whose
centroid is displaced. Splitting recovers all six.

### The cliff: gone

Depth to reach **100% of captured** mitoses, and the pool it is read from:

| ROI | cc depth | ws depth | change | pool growth | ws depth as % of pool |
|---|---:|---:|---:|---:|---:|
| **246** | 20,493 *(at 0.966 ceiling)* | **3,082** *(at 1.000)* | **6.6x better** | 2.71x | 5.4% |
| **201** | 4,794 *(at 0.889)* | **3,087** *(at 1.000)* | **1.6x better** | 3.28x | 7.6% |
| **300** | 4,246 | **2,763** | 1.5x better | 1.86x | 8.6% |
| **202** | 2,304 | **1,859** | 1.2x better | 2.25x | 6.1% |
| 405 | 3,630 | **3,531** | 1.03x better | 1.66x | 28.0% |
| 245 | **17,473** | 21,500 | 1.23x worse | 1.74x | 35.5% |
| 301 | **7,750** | 14,252 | 1.84x worse | 1.85x | 37.7% |

246.tiff is the headline: the two annotations that sat at ranks 17,320 and 20,492 were fused
components, and once split they rank inside the first 3,082 — **while four more mitoses that
were previously unreachable at any depth are now captured.** The cliff between 98% (1,966)
and 99% (17,321) disappears entirely; the watershed arm runs 1,842 → 3,033 → 3,082.

### It is not free, and it is not tuned

The candidate pool grows **1.66–3.28x**, and on **301.tiff and 245.tiff the deep tail gets
worse** — 1.84x and 1.23x. That is the over-splitting failure mode: `MIN_PEAK_DIST = 7` is a
single untuned constant applied to every domain, and where nuclei are large or chromatin is
coarse it carves one nucleus into several, none of which carries the whole object's
intensity. Note the two losses are *tail* losses; at 95% on 245.tiff the split still wins
(13,481 → 11,514), and only 301.tiff loses at both levels.

So the correct statement is: **splitting is unambiguously right for the ceiling (7/7) and
right for the reading depth on 5 of 7 ROIs, with one untuned hyperparameter that plainly
needs a sweep before it is shipped.** A union of the two candidate sets is the obvious cheap
hedge and has not been measured.

## 7. Caveats

* **These are 2 mm² ROIs, not WSIs — but 2 mm² is the right scale.** The mitotic count is
  clinically defined per 2 mm², which is why MIDOG's ROIs are that size, so these
  measurements are at the clinical unit of work. The open whole-slide question is *hotspot
  selection* — which 2 mm² to read — which is a separate and better-studied problem, not a
  scaling-up of anything measured here.
* Sensitivity is over *captured* mitoses throughout, so the 95% column on 246.tiff is 95%
  of 96.6% = 91.8% of annotated mitoses. Read the ceiling column with it.
* The greedy matcher assigns one detection per annotation at the ROI's match radius
  (29.6–33.1 px); depth is therefore "rank of the detection that claimed this annotation",
  which is the reader's position in the list, as intended.
* `245.tiff` is an outlier on every arm — 34,910 candidates for 90 mitoses, and 38% of the
  pool even at 95%. It and 246.tiff are the same domain and scanner, so this is not a domain
  effect; it is unexplained and is the single best target for a follow-up diagnostic.

## 8. The plan

Ordered by measured value per unit of work. Nothing below needs a GPU, a foundation model,
a pretrained network, or a training run longer than seconds.

### Step 0 — the spec, restated on the evidence

The original goal — *hit every mitotic figure, at a false-positive count that still saves
the pathologist time* — **is achievable, but only after the detector is fixed.** With
`nucleus_blobs` as it stands it is not: the proposal stage caps at 88.9% and 96.6% on two
ROIs (§2), so no ranker can reach 100%. After splitting merged nuclei (§6) the ceiling is
**100% on all seven ROIs**, and 100% sensitivity costs:

* **1,859–3,087 candidates on 246, 300, 201, 202** — 5–9% of the candidate pool. Against
  the pool the reader faces *today*, that is the number the product turns on: on 246.tiff,
  **3,082 candidates for every mitosis, against 21,059 nuclei in the unsplit pool — a 6.8x
  saving at full sensitivity**;
* 3,531 on 405 (28% of a much smaller pool);
* 14,252 and 21,500 on 301 and 245 — 36–38%, i.e. still the "review every nucleus" outcome.

So the spec to commit to is **100% of annotated mitoses at ≤ ~3,000 candidates per 2 mm²
ROI, met on 4 of 7 ROIs today**, with 301 and 245 named as the open cases rather than
averaged away. Report the ceiling and `n_pool` beside every depth — a change that halves
`frac_of_pool` while doubling the pool has made the reader's job worse, and the watershed
grows the pool 1.7–3.3x.

Pre-commit the metric before any further tuning: **depth to 100% of annotated mitoses,
worst of N seeds, leave-one-domain-out.** Not the median — the premise test already found
reading depth varying 2.2x with which cell was clicked.

### Step 1 — ship the watershed split, then sweep its one knob (~1 day)

This is the highest-value change measured in this session and it is ~30 lines of `skimage`.
It fixes the ceiling on 7/7 and improves reading depth on 5/7, including a **6.6x**
improvement on 246.tiff *while capturing four mitoses that were previously unreachable*.

Two things to do before it ships:

1. **Sweep `MIN_PEAK_DIST`.** It is a single untuned constant (7 px) applied across four
   tumour domains, and it is why 301.tiff and 245.tiff regress. Sweep 5–15 per domain and
   pick by depth-to-100%, not by pool size.
2. **Measure the union of the split and unsplit candidate sets.** It cannot lower the
   ceiling and it bounds the regression on 301/245. Cheap; unmeasured.

Do **not** ship "raise `max_area`" as the fix. Measured, it recovers 1 of 6 misses (§6).

### Step 2 — add the cheap supervised filter, scoped honestly (~1 day)

The "keep the blobs, add a step that filters false positives" idea from the original
question, measured in §5. 22 morphometric features, **logistic regression** — *not* gradient
boosting, which wins within a domain and loses across one on 6 of 7 ROIs. Trains in seconds,
scores 20k candidates in milliseconds, adds no dependency.

Ship it for what it does: **it rescues the worst ROI (245.tiff, 2.0x) and is a coin flip on
the others.** That matters because 245 and 301 are exactly the two ROIs Step 1 does not fix,
so the two changes are complementary rather than redundant. It is not an average-case win
and must not be quoted as one.

**Retrain it on the post-watershed candidates**, not the current ones. The features change
when the components change, and `intensity_mean` — the single most useful one — is precisely
the quantity that merging corrupts.

### Step 3 — the click (currently the other agent's thread, `click_rank_*`)

Everything measured in this log is **seedless**. The click was not used at any point, so
nothing here is evidence for or against it. What it does establish is the bar and the
remaining job:

* The premise test's +0.05 was a **median** result measured through `TM_CCOEFF_NORMED`. The
  click has **never** been measured against the tail, which is what the requirement is about.
* The click's defensible role is not as a template. It is (a) a **per-slide threshold
  setter** — the clicked cell's filter score is a cut point calibrated to this slide's stain
  and this tumour's morphology, which no seedless method can do — and (b) the **first label
  in an online loop**, where each verify/reject refits the Step 2 logistic regression in
  milliseconds. That is the dynamic process the original question described, and Step 2's
  model is already the right shape for it.

### Step 4 — frozen encoder embeddings, only if Steps 1–3 leave a gap

Gated behind the decision rule in `2026-09-01-one-click-retrieval-literature.md` §7. Do not
start it before Steps 1–2, and do not start it on the strength of an AUC.

### The resulting algorithm, end to end

No GPU, no pretrained weights, ~1–6 min per 2 mm² ROI on this laptop:

1. `tissue_mask` — grayscale < 220, morphological close. (Not PyHIST.)
2. Tiled Otsu on the hematoxylin channel → binary foreground.
3. Connected components; **watershed-split any component ≥ 2 × `min_area`**; area gate.
4. Score every candidate with the logistic regression on 22 morphometric features
   (weights trained offline on other domains, shipped with the app).
5. Present the ranked list; cut at the depth the target sensitivity requires.
6. **The click** sets the cut point for *this* slide, and every verify/reject refits the
   model in the loop.

Steps 1–3 are the detector; 4 is the "classification layer" from the original question;
5–6 are the product.

### What not to do

* **PyHIST** (§0) — it is tissue tiling, and `tissue_mask` already is that.
* **Ship the watershed untuned.** It regresses 301.tiff by 1.84x at the tail.
* **More correlation-search tuning.** Seven operating-point experiments have now come back
  neutral or reversed.
* **Gradient boosting** on this feature set at this data scale — measured worse across
  domains on 6 of 7 ROIs.
* **Quote any leave-one-ROI-out figure as a generalisation estimate** — it is optimistic by
  ~20% against leave-one-domain-out on this dataset (§5).

## 9. Reproducing this

```
~/anaconda3/bin/python tail_sensitivity.py    # ceiling + depth, 3 unsupervised rankers
~/anaconda3/bin/python tail_object_audit.py   # the 8 deepest mitoses per ROI
~/anaconda3/bin/python miss_attribution.py    # why the 11 problem annotations fail
~/anaconda3/bin/python fp_filter_probe.py     # 22 features (caches to .cache_tail/)
~/anaconda3/bin/python fp_filter_domain.py    # leave-one-domain-out, per-arm re-matching
~/anaconda3/bin/python watershed_split.py     # does splitting fix the ceiling and the cliff
```

Requires the **anaconda** interpreter, not the system `python3`: `requirements-midog-utils.txt`
pins numpy < 2 and the system install is on 2.4.6, where `import cv2` fails.

`results/fp_filter_domain.csv` supersedes `results/fp_filter_tail.csv` (§5), and
`results/miss_attribution.csv` supersedes `results/tail_misses.csv`, whose `cause` column
was derived by a nearest-centroid rule that §4 shows is wrong. Both superseded files are
kept so the correction is auditable.
