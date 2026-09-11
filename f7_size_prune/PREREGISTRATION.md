# Pre-registration F7: what does a physical-size prune on the candidate list buy?

Date: 2026-09-08. **Written before `f7_size_features.py` is run.** Status: PLAN.

F1-F6 are taken (F2 and F6 are in flight). This is the next free number.

**This is not F2's question.** F2 varies the *template's* `base_size` — the size of the thing the
matcher correlates with, which changes the candidate list itself. F7 changes nothing upstream: it
takes one fixed candidate list and deletes candidates whose *own measured object* is not a
plausible physical size for a mitotic figure. Template size and candidate size are different knobs
and this run holds F2's knob frozen at whatever the cached pool used.

---

## 1. The question

`Research Logs/2026-09-08-chromatin-vs-tm-and-residual-prune.md` closed the stain-purity axis of the
prune idea: `TM_CCOEFF` on `hematoxylin_od` already sorts stain purity to the top of its own list,
so a "this darkness is not chromatin" rule finds almost nothing above the working depth, and its
sign inverts on 459.tiff. That log's §5 leaves exactly one criterion open:

> The **physical-size** axis (µm, from `mpp`) is untested and is the remaining criterion that is
> universal by construction rather than by fitting. The TM candidate pool still has no size or
> tissue gate at all, unlike the blob path.

F7 measures it. **Q: with the search, the pool and the ranking held byte-identical, how much
reading burden does a post-NMS deletion by candidate physical size remove, and at what cost in
recall?**

## 2. Why this is a genuinely different lever from re-ranking

A prune is formally a re-ranker that sends what it rejects to the bottom, so no recall@K metric can
tell them apart. The reasons to prefer a prune here are epistemic, and they are the reasons this
experiment is worth running after five re-ranking nulls:

* **It uses information the ranker structurally cannot see.** `mpp` is TIFF metadata, not pixels.
  `cv2.matchTemplate` has no access to it; the same correlation value means a different physical
  size on a 0.226 µm/px scanner than on a 0.248 µm/px one.
* **The guarantee is one-sided and checkable.** "This rule lost 0 mitoses on the held-out domain" is
  one number. "This ordering is better" is a claim about the whole list.
* **It composes.** It is independent of whatever ranker ships later.

## 3. Design: one knob, and literally nothing else recomputed

**No search is re-run.** F7 reads the deep pool that `recall_workload_ledger.py` already committed:

| artefact | what it supplies |
|---|---|
| `results/tm_recall_workload_pool_z.npz` | `cx`, `cy`, `z` per (ROI, seed), **score-descending**; array order *is* the rank |
| `results/tm_recall_workload_cells.csv` | `seed_ann_id`, `mpp`, `match_radius_px`, `n_pool`, `n_gt_mitotic` |
| `results/tm_recall_workload_tp_ledger.csv` | each mitosis's rank in that pool — the reproduction gate |

So the padded ROI, `BORDER_REPLICATE`, `TM_CCOEFF`, `hematoxylin_od`, the largest-CC template
sizing, `PEAK_MIN_DISTANCE = 7`, the 5.0 µm NMS radius, the 5.0 px self-hit filter, `DEEP_FLOOR_Z
= -1.5` and the seed draw are all frozen **by construction, not by re-declaration** — they are
baked into the cached coordinates. The only new computation is a per-candidate size, and the only
new operation is a deletion.

**Sample.** All 14 ROIs of `images/extra_valid` × 5 seeds = 70 cells, per `DECISIONS.md` D5.
The ledger's seed draw is *with replacement across seed indices*, so some cells are duplicates
(201.tiff seeds 0/1 both draw ann 4457; 233.tiff seeds 0/1/4 all draw 5761). `n_distinct_clicks`
is reported per ROI and every ROI-level test runs on ROI medians, never on 70 rows treated as
independent.

## 4. The reproduction gate, which runs first

Re-run `ev.bucket_detections` on the **unpruned** cached pool for all 70 cells and assert the
resulting TP ranks equal `tm_recall_workload_tp_ledger.csv`'s `rank` column element-wise, per cell.
If that fails anywhere the run **stops**: it would mean the coordinates or the `gt_eval`
reconstruction diverge from what produced the ledger, and every downstream number would be void.

Two hazards handled explicitly:

* `match_radius_px` is read from the cells CSV, **not** recomputed from `mpp`, so no float
  difference can creep in.
* The npz arrays are used in stored order and never re-sorted — there is no `score` column to sort
  by, and array order *is* the rank.

`gt_eval` is rebuilt as `ds.image_annotations(annotations, fn)` minus the row whose `ann_id` equals
that cell's `seed_ann_id`, matching `recall_workload_ledger.py`'s construction exactly.

## 5. The one knob: candidate physical size

For each candidate at `(cx, cy)`:

1. `tm.read_padded_patch(gray_inverted, cx, cy, 51)` — the same channel and the same 51 px window
   the seed side uses. **`gray_inverted`, not `hematoxylin_od`** — `largest_cc_box` in
   `recall_workload_ledger.py` is applied to `gray_inv` via its `_check`, and F7 measures candidates
   on the same ruler as the seed. A reader would otherwise assume the hematoxylin channel.
2. `cv2.normalize(..., NORM_MINMAX)` to uint8, then Otsu — **mirroring `largest_cc_box` exactly.**
   The per-patch normalisation is kept: dropping it changes the segmentation, not just the gate.
3. `cv2.connectedComponentsWithStats(binary, connectivity=8)` — 8-connectivity matches
   `skimage.measure.label(..., connectivity=2)`.
4. **Component selection: the component containing the centre pixel**, falling back to the largest
   when the centre is background. *This deliberately differs from `largest_cc_box`.* On a seed
   patch the largest component *is* the annotated object, because the click is centred on it. On an
   arbitrary candidate in dense tissue the largest component in a 51 px window is frequently a
   neighbouring clump touching the border — which would systematically inflate measured size
   exactly in the dense ROIs (245, 246, 459) where the workload is worst. The largest-CC variant is
   recorded alongside as a consistency axis, never as the primary.
5. `size_px = max(component width, component height)` — the same statistic the seed side reduces
   its bbox to via `_odd_local(max(y1 - y0, x1 - x0))`. `size_um = size_px * mpp`.
   `equivalent diameter` (`2*sqrt(area/pi) * mpp`) is recorded as a secondary measure.

**Candidates whose size is undefined** (within 25 px of the ROI edge, so `read_padded_patch`
returns `None`; or a patch with no foreground) are **kept**, never pruned. Pruning them would
smuggle a border rule in as a second knob. Their count is reported per cell.

**No area, solidity or area-fraction gate.** `largest_cc_box`'s `min_area=50`,
`max_area_frac=0.85` and `min_solidity=0.5` are *not* applied: each would be an extra knob and F7
is a one-knob experiment.

## 5a. AMENDMENT (2026-09-08, before any result was seen): the instrument in §5 cannot measure

§5's statistic was checked before the run and **it does not measure physical size.** Recorded here
rather than silently replaced.

**What was found.** `cv2.normalize(..., NORM_MINMAX)` followed by Otsu forces a roughly 50/50 split
of *any* 51 px window regardless of its content — measured median foreground fraction 0.49 (013),
0.68 (245), 0.56 (459). The component under the candidate is therefore about half the window and
its bounding box pins to the window edge:

| ROI | window | `max(w, h)` **at the window ceiling** | `equiv_diam` median |
|---|---|---:|---:|
| 013 | 51 px (11.5 µm) | **64.5 %** | 7.91 µm |
| 013 | 101 px (22.9 µm) | **54.9 %** | 13.40 µm |
| 245 | 51 px (12.7 µm) | **64.2 %** | 9.75 µm |
| 245 | 101 px (25.1 µm) | **57.7 %** | 17.47 µm |
| 459 | 51 px (12.7 µm) | **54.2 %** | 9.41 µm |
| 459 | 101 px (25.1 µm) | **41.0 %** | 15.28 µm |

Two independent disqualifications. The bbox statistic saturates for the *majority* of candidates,
so an upper bound `HI` would have nothing to discriminate. And `equiv_diam` nearly **doubles when
the window doubles** (7.91 → 13.40 µm on 013) — a quantity that scales with the measuring
instrument is not a property of the object.

**The replacement.** Segment the whole ROI once with the segmentation `baselines.nucleus_blobs`
already uses — `tissue_mask` → `to_hematoxylin` → `_tiled_otsu_threshold(tile=512)` → connected
components — and give each candidate the size of the component its pixel falls in. There is no
per-candidate window, so there is no ceiling, and the threshold means "darker than this tissue's
nuclear threshold" rather than "darker than the median of these 51 px".

* **Primary statistic:** equivalent diameter, `2*sqrt(area/pi) * mpp`, in µm.
* Bbox max side and raw component area in px are recorded alongside.
* A candidate whose pixel lands on **no** component has size `NaN` and is **kept**, never pruned
  (§5's rule, unchanged). Its `on_nucleus` flag is recorded so the count is visible.
* `nucleus_blobs`' `min_area=80` / `max_area=4000` gate, its darkness ranking, and every
  solidity/area-fraction gate are **not** applied. Only the segmentation is reused. F7 stays a
  one-knob experiment.

**This choice was made on saturation and window-invariance grounds, with no sight of TP-vs-FP
separation.** No TP/FP contrast was computed on any statistic before the instrument was fixed.

**Two gates the replacement must pass before the run, pre-committed here:**

1. **Tile-invariance.** The size distribution must not move materially between `tile = 512` and
   `tile = 1024`. If it does, a window artefact has been traded for a tile artefact.
2. **Residual merging.** Report per ROI the fraction of candidates landing in a component larger
   than 4,000 px — `nucleus_blobs`' own validated nuclear maximum. If that fraction is large in the
   dense ROIs (245, 246, 459), adjacent nuclei are still merging under the global threshold, any
   upper bound is measuring clumps there, and a null in exactly those ROIs would be an artefact of
   the instrument rather than a result. **This must be reported with the results either way.**

**Both gates ran on 2026-09-08 before the extraction and both PASS.**

| ROI | tile | n components | on_nucleus | eqd p05 / med / p95 µm | frac area > 4,000 px |
|---|---:|---:|---:|---|---:|
| 013 | 512 | 85,930 | 0.170 | 1.80 / **7.75** / 15.27 | 0.033 |
| 013 | 1024 | 82,823 | 0.168 | 1.80 / **7.71** / 15.48 | 0.035 |
| 245 | 512 | 100,377 | 0.414 | 3.61 / **6.01** / 11.78 | 0.003 |
| 245 | 1024 | 100,602 | 0.413 | 3.58 / **5.99** / 11.76 | 0.002 |
| 246 | 512 | 166,756 | 0.586 | 2.41 / **7.17** / 18.03 | 0.054 |
| 246 | 1024 | 165,994 | 0.587 | 2.41 / **7.20** / 17.85 | 0.052 |
| 459 | 512 | 92,815 | 0.634 | 3.65 / **8.31** / 18.06 | 0.053 |
| 459 | 1024 | 92,759 | 0.635 | 3.64 / **8.34** / 18.12 | 0.054 |

*Gate 1 (tile-invariance)*: every median moves by < 0.5 % when the tile doubles, against the 70 %
move the window statistic showed for the same doubling. The measurement is a property of the
object. **Passes.**

*Gate 2 (residual merging)*: 0.3-5.4 % of candidates land in a component above 4,000 px, and the
dense ROIs are not the worst (245 is the *lowest* at 0.3 %). Merging is not dominating, so an
upper bound is not measuring clumps. **Passes.** The per-ROI fraction is reported with the results
regardless.

Medians of 6.0-8.3 µm with p95 of 11.8-18.1 µm are biologically plausible nuclear diameters, which
§5's instrument never produced.

*Coordinate convention, checked before the run*: the lookup is `lab[round(cy), round(cx)]`. If x
and y were transposed the annotated mitoses would show no enrichment over the pool. They do —
**013.tiff: GT mitoses land on a nucleus 0.889 of the time against 0.125 transposed** (pool 0.170,
segmentation foreground fraction 0.157); **459.tiff: 0.931 against 0.405 transposed** (pool 0.634).
True-positive *candidates* are on a nucleus 0.765 (013) and 0.938 (459) of the time, so the size
statistic is defined for most of the objects whose depth F7 measures. **Passes.**

**A third fact, not a gate, that bounds what F7 can find.** `on_nucleus` is only **0.170 on
013.tiff** and 0.41-0.63 on the others: most TM candidates do not land on a segmented nucleus at
all. Under §5's pre-registered rule those keep `NaN` and are never pruned, so on 013 the size knob
can act on 17 % of the pool at most. That rule is **not** being changed after seeing this — doing so
would be fitting — but it caps the achievable saving per ROI, and the cap must be quoted beside
every result. It is also a finding in its own right: the `TM_CCOEFF` candidate list is largely not
nucleus-centred.

**One consequence to carry forward.** The `on_nucleus` flag is *not* a "certainly not a mitosis"
criterion and must never be reported as one. `nucleus_blobs`' segmentation has a proposal ceiling
of 0.965 on 246.tiff (`2026-09-03-fp-reduction-framing.md` §3), i.e. it misses ~3.5 % of that ROI's
mitoses outright. Any candidate sitting on a mitosis the tiled Otsu failed to segment reads as
`on_nucleus = False`. If reported at all it is a separate, clearly-labelled arm with that floor on
its TP loss stated.

**Correction to the depth arithmetic.** The ledger's `rank` column is **1-based**
(`recall_workload_ledger.py` stores `ranks0 + 1`). The reproduction gate compares
`recomputed_index + 1 == ledger rank`. Earlier ad-hoc depth tables in this session and
`residual_prune_probe.py` added a further +1 and are therefore one too large (~0.01 % at a depth of
10,000; no conclusion changes). F7 uses the 1-based convention directly.

## 6. Where the bound comes from — three regimes, LODO is the headline

The prune is `keep if LO <= size_um <= HI`.

| regime | how LO/HI are set | what it answers |
|---|---|---|
| **a-priori** | **LO = 4, HI = 18 µm**, committed here | Does a bound chosen from biology alone, with no sight of the data, pay? |
| **in-sample** | quantiles of *this cell's own* TP sizes | The optimistic ceiling. Not deployable. |
| **LODO** | quantiles of TP sizes over the **other six domains** | The deployment condition, and the only one that speaks to "universal". |

The a-priori bound is a pre-committed guess the data may refute. Its basis: a condensed mitotic
chromatin mass runs roughly 5-15 µm across, widened to 4-18 µm for segmentation slack. It is **not**
derived from MIDOG's 50 px annotation box — that box is a fixed *pixel* size across scanners whose
`mpp` spans 0.226-0.248 in the cells CSV, so it carries no biological size statement.

**LODO holds out a whole domain (both its ROIs), not one ROI**, per `fp_filter_domain.py`'s finding
that the 14 come in domain pairs sharing tumour type *and* scanner, so leave-one-ROI-out measures
"generalises to another ROI of a tumour I have already seen". Seven folds.

Quantile levels: TP retention **1.00, 0.99, 0.98, 0.95**, so the trade is a curve rather than a
point. Reporting only the 1.00 level would repeat the residual probe's error — a min over ~100-240
TPs is an order statistic one atypical mitosis can pin.

## 7. What is measured, per (ROI, seed), before and after

Denominator is fixed at `n_gt_mitotic` from `gt_eval` in **both** arms, so a lost mitosis shows up
as a rising depth or a NaN rather than being hidden by a shrinking denominator.

* `n_pool` before / after → **candidate-list shortening**
* `deep_recall` before / after → **the ceiling**, reported as its own column so a lost TP is visible
* **depth to 90 %, 95 %, 100 %** of `n_gt_mitotic`, before and after (NaN when unreachable)
* the **paired per-cell delta** at each level — this is what the sign test runs on, not the two
  absolute tables

Ranking after the prune is the surviving candidates **in their existing order**. The pool is
already `TM_CCOEFF`-descending, so "prune then re-rank by TM_CCOEFF" is the identity on the
survivors; no re-sort is applied and none is needed.

## 8. Pre-committed reading

* **Primary contrast:** LODO bound at 1.00 retention, paired ROI-median delta in depth-to-100 %.
* The prune is worth shipping if it removes reading burden at **zero held-out TP loss** on a
  majority of the seven domains, with the ROI-level sign test reported alongside — and, per the
  residual probe's lesson, if the saving is not concentrated in the ROIs that were already cheap.
* **A sign flip between domains kills "universal"** exactly as it did for the residual criterion.
  That check is pre-committed: report LO and HI per LODO fold and whether any fold's bound would
  delete another fold's mitoses.
* n = 7 domains: `binomtest(0, 7)` bottoms out at p = 0.0156, so a unanimous result is the *most*
  significant outcome attainable. Report that floor rather than implying resolution the design
  cannot deliver.

## 9. Files

All under `f7_size_prune/`. No existing script is modified; `midog_utils` is imported, never edited.

* `f7_size_features.py` — reproduction gate, then per-candidate size extraction → `results/f7_candidate_sizes.npz`, `results/f7_repro_gate.csv`
* `f7_size_prune.ipynb` — the prune, the three regimes, the tables
* `results/` — all outputs
