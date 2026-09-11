# Audit of `tp_fp_separability.ipynb`: the separability is real, three of the conclusions built on it are not

Date: 2026-09-08. Script: `tp_fp_separability_audit.py` (eight checks, all re-derived from
`.cache_tail/tp_fp_candidate_features.csv`, the notebook's own input). Tables:
`results/tp_fp_separability_audit_*.csv`.

The notebook's central claim survives: measured against the matcher's own candidate pool rather
than against uniformly sampled nuclei, chromatin separates mitotic figures from false positives
at AUC 0.92–0.997 in all seven domains, and shape does not. Nothing below disputes that.

What does not survive is three things built on top of it, plus a set-aside that is half
mislabelled. Each is a defect in the **analysis**, not in the search — the candidate pool is
taken as given throughout, and it reproduces `results/f1_seed_sweep.csv` exactly as §1 claims.

---

## 1. The set-aside population is 52 % something else  (`C1`)

§3 removes 1,211 candidates it describes as *"false positives sitting on mitotic figures"*,
*"duplicate detections of true positives"*, and *"a reader who clicked one would not be wrong"*.
It builds them with

```python
gt = ds.image_annotations(anns, fn)          # <- pools category 1 AND category 2
gt = gt[gt.ann_id != int(seed_of[fn])]
d, _ = KDTree(gt[['cx','cy']].to_numpy()).query(...)
op['on_annotation'] = (op.bucket == UNANNB) & (op.dist_to_gt <= op.match_radius)
```

`dataset.image_annotations` applies no category filter, so `dist_to_gt` is the distance to the
nearest annotation **of either category**. Splitting it:

| | pooled | 094 | 201 | 246 | 301 | 402 | 459 | 548 |
|---|---|---|---|---|---|---|---|---|
| set aside as "on a mitosis" | 1,211 | 222 | 16 | 266 | 198 | 217 | 91 | 201 |
| actually within radius of a mitosis | 582 | 79 | 6 | 116 | 115 | 117 | 30 | 119 |
| within radius of a **look-alike** only | 629 | 143 | 10 | 150 | 83 | 100 | 61 | 82 |
| % that are look-alikes | **51.9 %** | 64 % | 63 % | 56 % | 42 % | 46 % | 67 % | 41 % |

Three consequences.

* §3's headline **"1.3 false positives per mitotic figure"** should be **0.65** (582 / 902).
* *"a reader who clicked one would not be wrong"* is false for 629 of them. A category-2
  annotation is a structure a pathologist examined and **rejected**.
* The 629 are second detections of pathologist-marked mimickers — precisely the population §10
  point 4 identifies as the one no feature in the set touches — and they are removed from the
  negative class before every AUC in the notebook is computed.

The quantitative effect is genuinely small: §7's own include-vs-exclude cell bounds it at
ΔAUC ≤ 0.003, and that bound is unaffected by the relabelling. The defect is that the notebook
is one AUC softer in the one place it is weakest, and states a mechanism for the wrong half of
the set.

**This definition is shared with `f1_seed_sweep.py:201-203`'s `n_dup_fp`**, which is why §3's
`assert t['subset ok'].all()` passes — the notebook is consistent with the repo and both are
consistent with a name that does not describe what is counted. Fixing it means a category filter
in both places, or renaming to `n_on_annotation_fp`.

§6's mechanism story needs re-splitting with it. Re-running the control on the two halves
separately (`C1`):

| | `od51` AUC vs mitosis-duplicates | `od51` AUC vs look-alike-duplicates |
|---|---|---|
| range over 7 ROIs | 0.705 – 0.874 | 0.783 – 0.947 |

so the "20–30 px off-centre, so the chromatin is weaker" mechanism §6 gives is the *weaker* of
the two effects being averaged, not the whole of it.

**What the same check clears.** The positives themselves are well centred — median 3.0–6.1 px
from their matched annotation, ≥ 92.6 % within 20 px in every ROI (`C3`). The mitotic class is
not contaminated with off-centre windows, so §6's finding does not turn round on the positives.

---

## 2. Every re-ranked reading depth is scored with the wrong labels  (`C6`)

§9's `depth()` sorts by `od51` or by the supervised score, then counts hits using `y2` — the
bucket labels `evaluate.greedy_match` assigned **walking the list in TM-score order**. Greedy
matching is order-dependent by construction (`evaluate.py`'s own docstring says so); under a
different ranking a different detection claims each annotation. A candidate that the new ranker
promotes above the TM-order winner is scored as a false positive while its mitosis is counted
found only further down.

Re-running the match in each ranker's own order, on the same pool:

| ROI | `od51` read_95, notebook | `od51` read_95, re-matched |
|---|---|---|
| 301.tiff | 9,585 | **7,345**  (−23 %) |
| 246.tiff | 3,391 | 2,819 |
| 548.tiff | 4,488 | 4,628 |

The bias is toward penalising every re-ranker, `z` alone excepted (its order *is* the matching
order). `results/tp_fp_reading_depth.csv` carries it, and `DECISIONS.md` D5 cites that file as
evidence.

---

## 3. `od_falloff` is never put through the metric the notebook says decides things  (`C6`)

§10 point 1 names `od_falloff` the strongest feature (AUC 0.939–0.993). §9 then measures reading
depth for `z`, `od51` and the supervised score — not for `od_falloff`. With the match re-run per
ranker:

| ranker | read_95 min | median | max |
|---|---|---|---|
| `od51` | 801 | 4,628 | 10,236 |
| `od_contrast` | 405 | 2,935 | 9,778 |
| `z` (shipped, D5) | 665 | 1,912 | 9,433 |
| `od_falloff` | 872 | **1,492** | **5,464** |
| mean of `z` and `od_falloff` per-ROI percentiles | 653 | **1,216** | 6,810 |

`od51` — the axis §9 chose as the feature baseline — is the **worst** of the four on this metric,
and it is the axis D5 already removed from the arbiter. `od_falloff` roughly halves `z`'s median
depth. No single axis wins everywhere: `od_falloff` is 3.4× *worse* than `od_contrast` on
201.tiff, and 459.tiff needs 5,464–10,236 whatever is used.

---

## 4. The supervised gain is quoted against a baseline D5 retired  (`C7`)

§10 point 7 reports the supervised leave-one-domain-out re-rank cutting `read_95` to
**0.13×–0.29×**. Those are ratios against `od51`. D5 states the production ranker is the
`TM_CCOEFF` score and that *"no experiment may declare `chromatin_od` its primary axis without
measuring, on that run's own data, that it beats `tm_score`"* — and on this run's own data it
does not (median read_95 4,628 vs 1,912).

Against the shipped ranker, with the match re-run:

| ROI | GB LODO read_95 | `z` read_95 | ratio |
|---|---|---|---|
| 301.tiff | 1,449 | 9,433 | 0.15 |
| 246.tiff | 753 | 1,912 | 0.39 |
| 459.tiff | 3,408 | 8,853 | 0.38 |
| 548.tiff | 734 | 2,075 | 0.35 |
| 402.tiff | 646 | 1,102 | 0.59 |
| 094.tiff | 1,451 | 1,630 | 0.89 |
| 201.tiff | 735 | 665 | 1.11 |
| median | **753** | **1,912** | **0.39** |

Still a real, transferring 2.5× — the notebook's finding stands, and looks better re-matched than
the 0.62 median its own `supervised / z` column reports. It is *not* the 3–8× the prose leads
with, and the prose picks the comparison that flatters it.

Two smaller notes on the same paragraph. 201.tiff has 17 mitoses, so `ceil(0.95 × 17) = 17` and
its `read_95` **is** `read_100` by construction; it is one of the two ROIs where the model
"loses", so "shorter on 5 of 7" is partly an n = 17 artefact — on an ROI §3 says carries no
weight. And §10 point 4's *"chromatin features drop to AUC 0.62–0.80"* against look-alikes omits
`od_falloff` on 459.tiff at **0.441** — the notebook's own best feature, reversed below chance.
The summary CSV has the honest range (0.441–0.760); the prose does not.

---

## 5. The k = 2 clustering result does not survive its own baseline  (`C5`)

§10 point 6: *"a 2-means split puts 99.0 % of the mitoses in one cluster holding 46.8 % of the
pool."* True. The comparison it is never given is the one-line alternative: keep the top 46.8 %
by per-ROI percentile on a single feature.

| | pooled recall at 46.8 % volume |
|---|---|
| k-means, k = 2, 8 features | 99.0 % |
| threshold on `od_falloff` | **99.2 %** |
| threshold on `z` | 98.6 % |
| threshold on `od_contrast` | 98.7 % |

Per ROI, the `od_falloff` threshold at each ROI's own k-means volume matches or beats k-means in
**7 of 7** (301: 98.6 vs 97.7; 459: 97.7 vs 96.2; the other five tie). The permutation null in §8
is careful and correct — it just nulls against chance, and the baseline that matters is not
chance, it is `sort()`. Clustering buys nothing here and costs a fitted, non-monotone,
initialisation-dependent object.

It is also fed a metric the notebook warned about. §4: *"any distance-based clustering would
weight the size axis several times over."* Inside the kept eight, `extent`–`solidity` sits at
ρ = **0.867** — above the 0.85 line the heatmap annotates as redundant — with
`tightened_size`–`area_frac_of_window` at 0.81 and `od51`–`od_contrast` at 0.80. Shape is counted
twice and size is counted twice in the k-means distance. §4's cell prints the 0.867 without
comment.

---

## 6. "Nothing below is limited by coverage" is not what was measured  (`C2`)

§1 prints `full_list_recall = 1.0` and concludes *"the search captures every mitotic figure, so
nothing below is limited by coverage"*. `evaluate.coverage_fraction` on the z ≥ 1.0 list:

| 201 | 548 | 094 | 459 | 301 | 402 | 246 |
|---|---|---|---|---|---|---|
| 0.896 | 0.939 | 0.942 | 0.965 | 0.974 | 0.982 | 0.992 |

89.6–99.2 % of every ROI lies within the match radius of some candidate. `evaluate.py`'s module
docstring names this the check to read *before* any un-budgeted recall number, for exactly this
reason. The separability analysis is unaffected — the positive class is complete either way, and
that is all §5 onward needs — but the sentence claims something the run did not establish.

A weaker instance of the same shape: §5's domain-predictability test (82.6 % raw → 16.4 %
ranked, chance 14.3 %). File name and domain are 1:1 across these seven ROIs, and the per-ROI
rank transform makes every feature's marginal uniform *within* each ROI by construction, so a
linear model on the ranked features is close to guaranteed to fail. The conclusion is right; the
test is much weaker evidence for it than it reads.

---

## 7. What the audit checked and found correct

Recorded so a future reader does not re-check them.

* `auc()`'s orientation. `scipy.stats.mannwhitneyu(a, b).statistic` returns U for the **first**
  sample in scipy 1.13.1 (the kernel's version), so `auc(pos, neg)` is signed the way the tables
  read. The legacy `min(U1, U2)` behaviour would have folded every AUC below 0.5.
* The permutation null in §8. It permutes labels **within** ROI, preserves each base rate, and
  correctly nulls the fact that `enrich_of_best` is a maximum over k clusters. This is the most
  careful construction in the notebook.
* No leakage in the leave-one-domain-out fit. The per-ROI rank transform is computed inside each
  ROI, so a held-out ROI's features never see training data, and ranking within the ROI you are
  deploying on is available at deployment.
* `depth()`'s quantile arithmetic (`ceil(q·m)`, then `+1`) is right and conservative — 207 of 217
  is 95.4 %, not 94.9 %.
* Candidate counts and template sizes reproduce `results/f1_seed_sweep.csv` on all seven ROIs, so
  §1's provenance gate does what it says.
* The mitotic class is well centred on its annotations (`C3`, above).

---

## 8. What this means for the shortlist question

The notebook is a **separability** result, and it is a good one. It is not yet evidence about a
deployable shortlist, for a reason that is structural rather than a defect: every `read_95` in it
— and every one in this audit — is an **oracle** quantity. It sorts the ROI and walks down until
95 % of that ROI's *known* mitoses have been found. At deployment the stopping point is the
unknown.

The repo already has the deployable form of the question on the score axis:
`results/tm_recall_workload_tolerance.csv` and `..._lodo.csv` fix a recall tolerance, derive the
cutoff on six domains and measure achieved recall on the seventh, 14 ROIs × 5 clicks. Its answer
at tolerance 0.95 is that a cutoff on `z` alone retains **3,591–18,311 candidates (median),
25,650 worst-case**. Thresholding the search score cannot produce a short list; that is settled.

Two candidate filters were measured here against that gap.

**Clustering — no.** §5 above: matched or beaten by a percentile threshold in 7 of 7 ROIs.

**Similarity to the clicked template — no, and the failure mode is the useful part** (`C8`).
Ranking candidates by Euclidean distance to the clicked cell's own feature vector, 30 clicks per
ROI:

| ROI | read_95 p10 | median | p90 | **worst click** | rec@1000 median | **rec@1000 worst click** |
|---|---|---|---|---|---|---|
| 402 | 708 | 1,030 | 1,398 | 6,083 | 94 % | 13 % |
| 246 | 1,478 | 1,756 | 2,164 | 14,118 | 88 % | 4 % |
| 094 | 1,909 | 2,122 | 2,760 | 13,816 | 79 % | 0 % |
| 548 | 1,679 | 2,186 | 3,845 | 9,239 | 79 % | 16 % |
| 301 | 2,880 | 5,484 | 6,939 | 21,914 | 82 % | 2 % |
| 459 | 6,448 | 8,698 | 12,305 | 20,842 | 47 % | 5 % |

The median is unremarkable — comparable to a single monotone feature — and the tail is
catastrophic: some clicks produce a ranking that needs most of the pool, and `rec@1000` collapses
to 0–16 %. The mechanism is that a **distance is symmetric**. A candidate that is *darker, more
sharply peaked, more mitotic-looking* than the clicked cell is pushed away from the prototype
exactly as hard as one that is paler. Click a faint or atypical mitosis and the ranking inverts.
Restricting the distance to the three chromatin features (`C8`, `proto_3chrom`) narrows the p90/p10
spread but does not fix the worst click.

(The prototype here is drawn from the pool's own mitotic candidates, not the real seed, which
passes `agreement_pool` and the largest-CC gates — so these are a lower bound. The symmetry
argument is structural, though, and the p10 column already fails to beat `od_falloff`.)

The conclusion is a shape constraint on any filter, and it is the one useful thing to carry
forward: **the ranking function has to be monotone in "more of the good thing", not a distance to
the click.** The click can set a *direction* or a *floor*; it cannot be a centre.

**Not measured, and the open question.** Whether the tolerance-and-LODO machinery of
`recall_workload_ledger.py` — which currently ranks on `tm_score` only, by design — produces a
usefully shorter list when run on `od_falloff`, on the `z`+`od_falloff` blend, or on the
supervised score. §3's median read_95 of 1,216–1,492 for the cheap blends and 753 for the
supervised model are the oracle numbers that would have to survive being converted into a rule
that picks K from ROI-observable quantities. Against a click-to-click SD of 260–2,744 on
`read_95` (`results/tm_ccoeff_headtohead_seed_variance.csv`), a single seed per ROI cannot
resolve the differences between the three cheap axes above — that comparison needs the 14 ROI ×
5 seed grid, not this notebook's 7 × 1.

---

# Part 2 — the corrected run, on 14 ROIs

Date: 2026-09-08, same day. Notebook: `tp_fp_separability_14roi.ipynb`. Extraction:
`tp_fp_feature_extract.py --images-dir images/extra_valid --select all` (52 min, 397,148
candidates). Tables: `results/tp_fp_14roi_*.csv`, `results/tp_fp_separability_14roi_summary.csv`.

All six defects above are fixed in the copy; the original notebook is untouched. Two further
changes were forced by the ROI count: leave-one-domain-out now holds out **both** ROIs of a domain
(with 2 per domain, holding out a file would be the sibling leak the 7-ROI version criticised in
`results/fp_filter_domain.csv`), and reading depth is measured on the **full** candidate pool,
duplicates included, because a reader at deployment sees them.

## The provenance gate

All 7 ROIs shared with the committed run reproduce it exactly. The notebook's gate checks the two columns `results/f1_seed_sweep.csv` actually carries (`n_pool`, `tightened_size`); `results/tp_fp_extract_summary_14roi.csv` against `results/tp_fp_extract_summary.csv` agrees on
`template_size`, `seed_ann_id`, `n_candidates`, `n_gt_mitotic`, `n_tp`, `n_lookalike`,
`full_list_recall` and `shape_fail_rate`. Adding seven ROIs to the loop changed nothing about the
original seven, which is what the per-image `default_rng([0, image_id])` seed stream promises.

An earlier bit-exact gate on the feature block *failed*, and the failure was in the gate: `score`
and `z` are `float32` in memory and the cache writes them at 7–8 significant digits, so
`float(np.float32(1.194488)) != float('1.194488')` by 5e-8. A md5 over CSV-parsed values compared
against in-memory float32 can never match. Recorded because it cost 13 minutes and looked, for
about a minute, like an order-dependence bug.

## What the corrections changed

| | 7 ROIs (as published) | 14 ROIs (corrected) |
|---|---|---|
| set aside as "on a mitosis" | 1,211 (**52 % were look-alikes**) | 913 — look-alikes kept in the negatives (900 of them) |
| duplicates per mitotic figure | 1.3 | **0.70** |
| effect of the set-aside rule on any AUC | ≤ 0.003 | ≤ 0.0022 |
| `coverage_frac` at z ≥ 1.0 | not measured | **0.743 – 0.998** |
| max \|ρ\| inside the kept 8 | 0.867, unremarked | 0.856, named (`extent`–`solidity`) |

## Three findings that did not survive the larger set

**1. `od_falloff`'s advantage over the shipped `z` is gone.** Median `read_95` **1,954.5 vs
1,993.5** — a dead heat, where 7 ROIs gave 1,492 vs 1,912. The click-to-click SD of `read_95` is
260–2,744 (`results/tm_ccoeff_headtohead_seed_variance.csv`); a 39-candidate gap is not a result.
What survives is the *worst case*: 9,487 vs 14,988.

**2. The separability floor was optimistic by ~0.07 AUC.** `245.tiff` (canine lymphosarcoma, the
second ROI of its domain) reads 0.854–0.905 where the 7-ROI minimum was 0.92, with an AP lift of
9–17x against the 37–277x that version reported.

**3. The clustering-vs-threshold margin narrowed.** Pooled, k = 2 k-means now edges the
matched-volume `od_falloff` threshold 98.8 % to 98.3 % (on 7 ROIs the threshold won outright, 7/7).
Per ROI the threshold still matches or beats it in **13 of 14**. The conclusion is unchanged in
direction and weaker in strength: the two are within a few tenths of a point, which is a reason not
to ship the fitted object, not a reason to call it beaten.

## The finding that only appears under the stated rule

The operating rule set for this work is **95 % mean recall across ROIs with a 90 % floor on every
individual ROI**. Under it, every ranker's binding constraint is the floor, and the ordering
inverts:

| ranker | median `read_95` | K for mean ≥ 95 % | K for every ROI ≥ 90 % |
|---|---|---|---|
| `od_contrast` | 3,074.0 | 3,232 | **5,655** |
| supervised LODO | **1,283.0** | **1,628** | 6,078 |
| `z+od_falloff` | 1,781.0 | 2,556 | 7,022 |
| `od_falloff` | 1,954.5 | 2,603 | 7,280 |
| `od51` | 3,709.5 | 4,637 | 9,895 |
| `z` (shipped, D5) | 1,993.5 | 5,357 | 12,739 |

Two claims come out of that table and only one is resolved.

**Resolved.** Against `z`, the ranker the repo ships, the cheapest rule costs **5,655 rather than
12,739** — 2.3x — from one within-image difference (`od51 − od_ctx`) needing no training and no
calibration.

**Not resolved.** That `od_contrast` beats the supervised leave-one-domain-out model, 5,655 to
6,078. **All six floors are set by the same ROI**, `245.tiff` — its `read_90` is the maximum for
every ranker (od_contrast 5,561, supervised 6,018, z+od_falloff 6,964, od_falloff 7,175, od51
9,879, z 12,554) — so the ordering is a ranking of one ROI's depth curve at one click.
`results/tm_ccoeff_headtohead_seed_variance.csv` puts that ROI's own click-to-click SD of `read_95`
at **2,153** candidates (5 seeds, mean 11,190, range 9,762–14,956). A 423-candidate margin is
**0.2 SD**; the 7,000-candidate margin over `z` is more than 3. `od_contrast` was also picked
best-of-six on the same 14 ROIs it is scored on.

What survives as a finding is structural rather than a winner: **a floor is priced by the hard ROI,
so median reading depth is the wrong statistic for choosing a ranker under this brief.** On the
median the supervised model is 2.4x ahead of `od_contrast`; under the floor it is 7 % behind.

At K = 5,655: mean recall 97.1 %, worst ROI 91.0 %, 4 of 14 ROIs below 95 %, and K is 16–43 % of
each ROI's pool. **The list is shorter; it is not short.**

## What is still owed

`K` above is an **oracle** — chosen knowing each ROI's mitosis count. The deployable form has to
pick `K` from ROI-observable quantities, fitted on other domains and validated on the held-out one.
`results/tm_recall_workload_tolerance.csv` and `..._lodo.csv` are exactly that shape on the
`tm_score` axis (14 ROIs × 5 seeds); running that machinery on `od_contrast` is the next piece of
work. It also supplies the seed dimension this notebook lacks — one click per ROI cannot resolve a
39-candidate median gap, and should not be asked to.
