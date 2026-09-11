# F7 results: the physical-size prune is a null, and where it bites it makes the read *longer*

Date: 2026-09-08. Design of record: [`PREREGISTRATION.md`](PREREGISTRATION.md), including the §5a
instrument amendment. Every number below comes from `results/` in this folder; nothing was
re-specified after the run.

---

## 0. Gate

**Reproduction gate: 70/70 cells exact.** Re-running `ev.bucket_detections` on the *unpruned*
cached pool reproduces `results/tm_recall_workload_tp_ledger.csv`'s ranks element-wise on every
cell, and the recomputed `n_gt_mitotic` matches the cells CSV 70/70. The fast matcher used for the
arms was separately asserted equal to `ev.bucket_detections` on 3 cells × 2 arms.

Seed draw is with replacement, so the 70 cells are **57 distinct clicks** (3-5 per ROI). Every
test below is at ROI level on ROI medians, never on 70 rows treated as independent.

## 1. The headline

| arm | bound µm | frac kept | ROIs losing ≥1 mitosis | Δdepth90 | Δdepth95 | Δdepth100 |
|---|---|---:|---:|---:|---:|---:|
| **lodo@1.00** (primary) | **[0.50, 99.85]** | **0.998** | 1/14 | **0** | **0** | **0** |
| lodo@0.99 | [2.82, 24.68] | 0.974 | 1/14 | −2 | 0 | −33 |
| lodo@0.98 | [3.62, 22.27] | 0.960 | 4/14 | +10 | +138 | −7 |
| lodo@0.95 | [3.96, 16.35] | 0.922 | 8/14 | **+550** | **+1,980** | +661 |
| **apriori** | [4, 18] | 0.918 | **8/14** | **+128** | **+1,115** | +48 |
| **insample@1.00** (ceiling) | [3.98, 20.11] | 0.933 | **0/14** | **−76** | **−142** | **−341** |

**The pre-committed primary contrast is inert.** The LODO bound at zero held-out TP loss keeps
99.8 % of the pool and moves the median depth by 0 at every recall level. It is not a prune; it is
a no-op that still manages to lose a mitosis on 402.tiff.

**Every LODO operating point that deletes enough to matter loses mitoses**, and at 0.95 retention
it makes the reader go **1,980 candidates deeper** to reach 95 %.

**Read the sign tests with the ties in view.** `lodo@1.00` reports "depth_95 shorter on 1/14,
p = 0.0018", which sounds like a strong negative; it is not. Twelve of the fourteen deltas are
exactly **0**, and the test counts a tie as a failure. The honest statement is that the arm is
tied on almost every ROI, not that it reliably hurts.

## 2. Why LODO is inert: there is no size interval that transfers

Per-fold bounds (hold out a whole domain, both ROIs; set the bound on the other six):

| held-out domain | n TP | bound from the others | that domain's own range | frac deleted |
|---|---:|---|---|---:|
| mast cell | 1,115 | [0.81, 99.85] | [0.50, 9.96] | 0.004 |
| lung | 163 | [0.50, 99.85] | [4.24, 22.36] | 0.000 |
| lymphosarcoma | 880 | [0.50, 99.85] | [1.81, 25.16] | 0.000 |
| soft tissue | 783 | [0.50, 99.85] | [1.43, 20.57] | 0.000 |
| breast | 452 | [0.50, 99.85] | [4.14, 33.66] | 0.000 |
| melanoma | 1,218 | [0.50, 99.85] | [2.82, 23.33] | 0.000 |
| neuroendocrine | 564 | [0.50, **33.66**] | [0.81, **99.85**] | **0.018** |

Six of seven folds inherit **[0.50, 99.85] µm** — a 200-fold interval that excludes nothing. Pooled
over all 70 cells the measurable TP diameters run min 0.50, p01 3.63, median 6.15, p99 21.34,
max 99.85 µm. **Two objects destroy the bound for every other domain**: a 0.50 µm "mitosis"
(segmentation noise) and a 99.85 µm one (a merged clump in 402.tiff). That is the min/max
order-statistic fragility §6 anticipated, and tightening to a quantile does not rescue it — the
retention curve never reaches a point that both saves and keeps every mitosis.

## 3. Why the a-priori bound is a net cost, and the mechanism that makes a prune *harmful*

`[4, 18] µm`, committed before any data was seen, keeps 91.8 % and **loses mitoses on 8 of 14
ROIs** (worst: 201.tiff, −0.118, i.e. 2 of 17). On four ROIs — 300, 201, 459, 548 — depth-to-100 %
becomes **undefined**, because the mitosis is gone.

But the more instructive number is that **depth gets longer where no mitosis was lost at all**.
233.tiff keeps full recall and its depth-to-95 % goes **1,328 → 17,341**, a 13× *increase*, from
deleting 8 % of the pool.

The mechanism is the one the design insisted on re-matching for: `evaluate.greedy_match` gives each
annotation to the **best-ranked** detection within the radius. Delete that detection and the
annotation is re-claimed by a **much deeper** one. The mitosis is still found, so recall is
unchanged — and the reader now has to go 13× further to reach it.

**This is invisible to any evaluation that reuses stored ranks after a deletion**, which is exactly
what a naive prune-then-report pipeline would do. It is the single most important methodological
finding of F7.

## 4. What the in-sample ceiling says

A bound fitted to each cell's *own* TPs at zero loss keeps 93.3 % and is shorter on **14/14 ROIs at
every level** (median −76 / −142 / −341, p = 0.0001). So the signal is real but small: 2-8 % off the
read, best case 013.tiff 366 → 260 (−29 %) and 460.tiff 10,622 → 8,345 (−21 %).

The distance between that and the inert LODO arm **is** the result: the size of a mitotic figure is
learnable within a slide and does not transfer between domains.

## 5. Instrument caveats, all pre-committed in §5a and reported regardless

* **`on_nucleus` = 0.22-0.71 per ROI** (013 0.25, 403 0.22, 459 0.71). Most TM candidates do not
  land on a segmented nucleus; their size is `NaN` and they are **kept by rule**. Only **78.9 %
  (5,175 / 6,555)** of true positives even have a measurable size. The knob can act on roughly half
  the pool at best, and that cap is part of the null.
* **Tile-invariance** (independent audit, all 70 cells): median |Δsize| between tile 512 and 1024 is
  0.0000 µm, but 9.6 % of candidates move >1 µm and **1.7 % of keep/drop decisions flip** (max
  4.9 %). Aggregates are stable; individual assignments are ~2 % unstable.
* **Residual merging**: 0-15.3 % of candidates sit in components above `nucleus_blobs`' validated
  4,000 px nuclear maximum (402 15.3 %, 094 14.4 %), and on four ROIs **5-11 % of the true positives
  themselves** are in merged components (201 10.6 %). 402.tiff's 99.85 µm "mitosis" is this.
* The bounding-box variant of the statistic is far worse — the independent audit measured TP
  retention under [4, 18] µm of median 0.822, **worst 0.396** — which is why §5 named equivalent
  diameter primary.

## 6. Verdict

**The physical-size prune does not work, on the axis it was proposed on.** At zero held-out TP loss
it is a no-op; at any bound tight enough to remove candidates it removes mitoses and lengthens the
read. The a-priori biological bound is a clear net cost. `DECISIONS.md` should record the size axis
as closed alongside the stain-purity axis
(`Research Logs/2026-09-08-chromatin-vs-tm-and-residual-prune.md` §5).

Two things survive and are worth more than the null:

1. **Re-matching after a deletion is mandatory.** §3's 13× lengthening on a full-recall ROI is the
   proof. Any future prune experiment that reuses stored ranks is measuring a fiction.
2. **The candidate list is not nucleus-centred** (`on_nucleus` 0.22-0.71). Every prune proposed so
   far — stain purity, physical size — has assumed the candidate *is* a segmentable object. Half of
   them are not, and that is a property of the `TM_CCOEFF` peak set worth understanding before a
   third such experiment is designed.

## 6b. Follow-up: the candidate list *is* nucleus-centred where it is read

§6 left one question open — why do only 22-71 % of `TM_CCOEFF` peaks land on a segmented nucleus?
Answered here from the cached score-descending arrays, no new segmentation
(`f7_on_nucleus_by_rank.csv`, `results/on_nucleus_by_rank.log`).

**`on_nucleus` is strongly rank-dependent. The head of the list is nucleus-centred; the tail is
not.** ROI medians over seeds:

| | top 100 | top 1,000 | q0 | q1 | q2 | q3 | whole pool | true positives |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| pooled median | **0.990** | **0.944** | 0.87 | 0.68 | 0.35 | 0.05 | **0.463** | **0.939** |

The head-minus-tail gap `q0 − q3` is **positive on 14 of 14 ROIs**, median **+0.696**, sign
p = 0.0001. Extremes: 013.tiff 0.820 → 0.009, 233.tiff 0.948 → 0.015, 529.tiff 0.880 → 0.031.

**The low overall figure is an artefact of the unread tail.** In the top 100 candidates, 99 % sit on
a segmented nucleus; in the top 1,000, 94 %; and 93.9 % of true positives do. The peak set is fine
where anyone looks at it.

**This is the same pattern for the third time, and it explains all three nulls at once.** Stain
purity (`2026-09-08-chromatin-vs-tm-and-residual-prune.md` §4), physical size (§1-2 above), and now
"is this even a nucleus" are each **prunable exactly where pruning is worthless** — deep in a list
the reader never reaches. It also disposes of the `on_nucleus` arm §5a deliberately held back: its
1.85x median true-positive enrichment over the pool is a *rank* artefact, and at the working depth
there is almost nothing left for it to delete.

**Two ROIs are different, and it is the segmenter, not the search.** 301.tiff has true-positive
`on_nucleus` of **0.097** against a pool rate of 0.245 — an enrichment of **0.40x**, the only
inversion in the set — and 403.tiff runs 0.269 against 0.142. These are the two ROIs with the
fewest measurable true positives (301: 21 of 217). `nucleus_blobs`' tiled Otsu fails on them —
plausibly the granules in mast cell tumour — which is a property of `baselines.py`, not of
`TM_CCOEFF`. Any future experiment that leans on that segmentation must report its per-domain
mitosis-miss rate first.

## 7. Files

`f7_size_features.py` (extraction + gate) → `results/f7_repro_gate.csv`,
`results/f7_candidate_sizes.npz`. `f7_analysis.py` (all arms) → `results/f7_prune_results.csv`,
`f7_retention_curve.csv`, `f7_lodo_fold_bounds.csv`, `f7_roi_*.csv`, `f7_on_nucleus_by_rank.csv`. `f7_prune_eval.py` holds the
prune, the re-match and the regimes. `results/f7_independent_audit.csv` is the separate
verification run. Reproduce with
`/Users/mohinianand/anaconda3/bin/python3.11 f7_size_prune/f7_analysis.py`.
