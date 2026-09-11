# `chromatin_od` vs `tm_score`, and the stain-residual prune that the comparison suggested and the data refuses

Date: 2026-09-08.
Question asked: *"Haven't I already experimented with ranking on the chromatin component? Is there
a fair comparison between TM_CCOEFF and chromatin, and are the differences significant? From
whatever you find, how feasible is a prune on the stain-unmixing residual?"*

Status: **analysis + one new probe**. §1-2 are re-derived from `results/tm_ccoeff_high_z_deployable_cutoffs.csv`
and `results/tm_ccoeff_threshold_axis_sweep_largest_cc_high_z.csv`, both already in the repo.
§3-4 are new, measured by `residual_prune_probe.py` (written for this log), seed 0, 7 ROIs.
**Nothing here is decision-grade**: `DECISIONS.md` D5 requires the 14 `images/extra_valid/` ROIs
for any run comparing a chromatin axis, and `2026-09-04-largest-cc-z-headroom-audit.md` Tier 2 C
already flags this run's `chromatin_od` arm as outside that rule.

---

## 1. Yes, a fair comparison exists, and it is properly paired

`tm_threshold_axis_sweep_largest_cc_high_z` builds **one** deep pool per (ROI, seed) and ranks it
two ways -- `rank_key='score'` (`TM_CCOEFF`) and `rank_key='od'` (`chromatin.chromatin_density`,
mean of the darkest 10% of a 51 px hematoxylin-OD window). Identical candidates, identical ground
truth, identical match radius, re-matched per arm. Every comparison below is a within-ROI paired
difference at the deployable cutoff `z = 1.0`. 7 ROIs, seed 0.

## 2. The axes cross over, monotonically, and only the tail is resolved

| metric | chromatin better | median ratio chrom/tm | sign p |
|---|---|---:|---:|
| `read_50`  | 5/7 | **0.67x** | 0.453 |
| `read_80`  | 4/7 | 0.79x | 1.000 |
| `read_90`  | 3/7 | 1.10x | 1.000 |
| `read_95`  | 2/7 | 1.20x | 0.453 |
| `read_99`  | **0/7** | **1.55x** | **0.016** |
| `read_100` | **0/7** | **1.26x** | **0.016** |

recall@K tells the same story from the other side: chromatin loses at every budget except K = 250
(2/7 at K = 100, 4/7 at 250, 3/7 at 500, 2/7 at 1000 and 2000, 1/7 at 5000).

**Reading.** `chromatin_od` is a better *precision* signal and a worse *recall* signal. It reaches
the first half of the mitoses in a third fewer candidates and the last few in half again as many.
That is a coherent mechanism, not noise: absolute darkness is a strong marker of a *typical*
mitotic figure and an actively misleading one for a faint or telophase figure, which a darkness
ranking buries.

**Two caveats that bound the strength of the 0/7 rows.** (i) `binomtest(0, 7)` bottoms out at
2/2^7 = 0.0156, so p = 0.016 means "unanimous at n = 7" and nothing stronger. (ii) Seed 0 only;
`results/tm_ccoeff_headtohead_seed_variance.csv` puts the seed-to-seed SD of `read_95` at
260-2744 candidates, which is the same order as several of these deltas.

**Consequence for the stated goal.** "All the TPs from one template in the fewest candidates" is
the `read_100` column, and it is the one column where the two axes disagree unanimously and in
`tm_score`'s favour. Chromatin is the wrong axis for that target. D5 was right to keep it off the
arbiter, and this is a second, independent reason.

## 3. The residual prune: measured, and it does not pay

The crossover suggested a mechanism -- that unnormalised `TM_CCOEFF` promotes high-amplitude dark
*non-chromatin* material (pigment, RBC clusters, hemosiderin, ink), which a "fraction of unmixed OD
that is hematoxylin" criterion would delete without touching a mitosis. `residual_prune_probe.py`
tests it off the cached pool coordinates (`results/tm_recall_workload_pool_z.npz`), with the
deletion rate measured **only over candidates above that cell's depth-to-100%** -- the sole region
a prune can pay in -- and the cut placed at a TP quantile rather than at min(TP).

`h_frac` = hematoxylin share of clipped `rgb2hed` OD over the darkest 10% of a 31 px window.
Ordered by the separation actually available, `tp_med - fp_med`, not by the 1.00-retention column:

| ROI | domain | depth100 | gap | del @ 1.00 | @ 0.98 | @ 0.95 |
|---|---|---:|---:|---:|---:|---:|
| 548 | melanoma | 7,950 | **+0.062** | 2.4 % (179) | **26.0 % (1,983)** | 33.0 % |
| 094 | breast | 4,053 | **+0.061** | **22.7 % (883)** | 25.6 % (998) | 38.9 % |
| 301 | mast cell | 15,930 | +0.049 | 5.7 % (666) | 12.5 % (1,463) | 17.0 % |
| 402 | neuroendocrine | 5,458 | +0.048 | 9.7 % (513) | 15.4 % (814) | 39.7 % |
| 245 | lymphosarcoma | 16,227 | +0.003 | 0.3 % (40) | 2.0 % (232) | 5.1 % |
| 246 | lymphosarcoma | 6,269 | +0.003 | 1.0 % (63) | 1.2 % (71) | 1.7 % |
| 459 | soft tissue | 15,467 | **-0.039** | 0.03 % (4) | 0.3 % (33) | 0.6 % |

**Do not read the 1.00-retention column alone.** It is a min-order statistic over 79-236 TPs and one
atypical mitosis pins it: 548 goes 2.4 % -> 14.7 % -> 26.0 % on retention 1.00 -> 0.99 -> 0.98, an
11x swing from dropping two objects. 548 has the *widest* TP/FP separation of the seven and reads
as the worst ROI in that column. 245 has the same shape (0.3 % -> 2.0 %) on a gap 20x smaller.

**The verdict is the split, not a median across it.** Three groups, and the median of seven numbers
spanning 0.03 %-22.7 % summarises none of them:

* **Real but small signal** -- 548, 094, 301, 402, gap 0.048-0.062. At 98 % TP retention the rule
  deletes 12-26 % of the working region. Best case is 094: 883 candidates off a 4,053-deep read at
  *zero* TP loss, a 1.29x. On 301, the deepest of the four, 666 off 15,930 is 1.04x.
* **No signal** -- 245, 246, gap 0.003. 1-2 % at 98 % retention. 245 is the deepest read in the set.
* **Inverted sign** -- 459, gap **-0.039**: its mitoses are *less* hematoxylin-pure than the false
  positives above them, so the rule as written would have to be reversed there to delete anything.

That last row is the finding that closes the idea as posed. **A criterion whose sign flips between
domains is not a universal criterion**, and 459 is not a borderline case -- it is the only ROI where
`h_frac` also *rises* with rank depth (§4) and where `tp_shallow < tp_deep`. Whatever it is
measuring in soft-tissue sarcoma, it is not the quantity the rule assumes.

## 4. Why it fails, which is the part worth keeping

The proposed mechanism is **refuted**. Median `h_frac` of false positives by rank quartile
(q0 = best-ranked):

| ROI | q0 | q1 | q2 | q3 | TP shallow | TP deep |
|---|---:|---:|---:|---:|---:|---:|
| 548 | 0.847 | 0.829 | 0.821 | 0.816 | 0.904 | 0.875 |
| 301 | 0.532 | 0.514 | 0.502 | 0.489 | 0.602 | 0.531 |
| 094 | 0.472 | 0.435 | 0.405 | 0.393 | 0.488 | 0.473 |
| 402 | 0.924 | 0.906 | 0.892 | 0.886 | 0.962 | 0.938 |
| 245 | 0.286 | 0.287 | 0.289 | 0.289 | 0.294 | 0.289 |
| 246 | 0.377 | 0.376 | 0.378 | 0.380 | 0.379 | 0.383 |
| 459 | 0.515 | 0.530 | 0.534 | 0.534 | 0.457 | 0.522 |

`h_frac` **falls** from q0 to q3 on 4 of 7 ROIs, is flat on 2, and rises on 1. **`TM_CCOEFF` on the
hematoxylin-OD channel already sorts stain purity to the top of its own list.** Correlating against
a nucleus-shaped hematoxylin template excludes pigment and red cells before the ranking starts, so
there is no non-chromatin population left above the working depth for a residual criterion to
remove. What survives there is hematoxylin-positive nuclei -- precisely
`2026-09-03-fp-reduction-framing.md` §2's finding that ~95 % of the FP mass is ordinary nuclei
rather than mimickers, now confirmed on the `TM_CCOEFF` pool rather than the blob pool.

The weaker mechanism (b) has mild support: `tp_shallow > tp_deep` on 5 of 7, i.e. deep mitoses are
somewhat less hematoxylin-pure -- the *wrong* sign for the prune, since it puts the budget-setting
mitoses nearest the cut.

**A second, independent blow to "universal".** Median TP `h_frac` ranges from 0.291 (245) to 0.950
(402) -- a 3.3x spread across domains. No absolute threshold transfers; only a per-ROI relative cut
could, which is still click-free but is not the scanner-independent constant the idea promised.
`chromatin.py:126-128` already warned that OD is not calibrated across scanners; this is that
warning showing up in a ratio that was supposed to be immune to it.

## 5. What survives

* `tm_score` is the right axis for a 100 %-recall target, on this evidence and D5's.
* The **stain-purity** axis of the prune idea is **closed**. Do not spend more on it.
* ~~The **physical-size** axis (µm, from `mpp`) is untested~~ — **tested and closed, 2026-09-08.**
  F7 (`f7_size_prune/`, 14 ROIs x 5 seeds, reproduction gate 70/70 exact) measured it as a
  post-NMS deletion. At zero held-out true-positive loss the leave-one-domain-out bound keeps
  99.8 % of the pool and moves the median depth by 0 at every recall level; six of seven folds
  inherit the bound **[0.50, 99.85] µm**. The pre-registered a-priori bound of [4, 18] µm loses
  mitoses on **8 of 14** ROIs. Fitted within a slide it works (14/14 ROIs shorter, p = 0.0001), so
  mitotic size is learnable per-slide and does not transfer. See `DECISIONS.md` D6.
* The `read_50` crossover is real (chromatin reaches the first half 33 % faster) and is **not**
  being carried forward. Exploiting it means a two-axis reranker; F1-F5 and
  `2026-09-01-click-ranking-experiment.md` have each measured a reranker and each found it fragile
  out of domain, and nothing here changes that prior. Recorded so the next reader does not
  rediscover the crossover and re-run the experiment that its obvious exploitation implies.

## 5b. Side finding: the deep tail is **not** enriched in contested annotations

Computed on `results/tm_recall_workload_tp_ledger.csv` (14 ROIs x 5 seeds, every captured mitosis's
rank in its own deep pool), which closes the hypothesis `2026-09-03-fp-reduction-framing.md` §6b
raised from n = 32 and flagged as untested.

Per ROI, the **deepest 10 % of true positives are less unanimous on only 7 of 14 ROIs** (sign test
p = 1.0). Pooled unanimity by rank decile is flat — 0.733, 0.771, 0.762, 0.798, 0.755, 0.747, 0.762,
0.738, 0.731, 0.777 from shallowest to deepest, no trend, against an overall 0.757. The largest
single per-ROI gap is 460.tiff at −0.460, and 094.tiff runs the other way at +0.155.

**Two caveats belong with this null.** The 6,555 rows are not independent — the same annotation
recurs across up to 5 seeds — so the ROI-level sign test is the valid reading and the pooled decile
counts are inflated. And per `dataset.py`, for category 1 the only label multisets in MIDOG++ are
`(1,1)` x 8,917 and `(1,1,2)` x 3,020, so `unanimous` is exactly `n_votes == 2`: a proxy for "did
this need a third reader", not a measure of label quality.

Consequence: "100 % of annotated mitoses" does **not** disproportionately mean "100 % of objects the
pathologists argued about". The objects that set the reading budget are ordinary-agreement mitoses,
and the cost of the tail cannot be discounted on labelling grounds. `DECISIONS.md` "Still open"
item 3 is closed accordingly.

## 6. Reproducing

`residual_prune_probe.py` (seed 0, 7 ROIs, ~4 min). Writes
`results/residual_prune_probe.csv` and `results/residual_prune_mechanism.csv`. It re-matches
nothing -- it reads `results/tm_recall_workload_pool_z.npz` and
`results/tm_recall_workload_tp_ledger.csv`, both committed. §1-2 are pandas pivots over
`results/tm_ccoeff_high_z_deployable_cutoffs.csv`.
