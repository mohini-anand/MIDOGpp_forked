# Audit of `tightening_process_hem_vs_gray.ipynb`: both headline effects reproduce exactly and the run is clean, but the "extreme" example quotes a size the notebook's own figure contradicts, the stated mechanism fails on the four pairs where `gray_inverted`'s box is pinned to the window ceiling, and "the same click at `seed_index=0`" is true of a notebook it does not name and false on 3 of 7 ROIs of the one it does

**Scope.** Target: `production_hematoxylin_only/tightening_process_hem_vs_gray.ipynb` (untracked,
executed in place 2026-09-16 11:27, 11 cells, 5 code cells, 14 embedded figures). Artifacts:
`production_hematoxylin_only/tightening_process_summary.csv` (untracked, written 11:27:55) and the
14 `production_hematoxylin_only/tightening_*.png` it wrote. Primary sources:
`databases/MIDOG++.json` and the 14 `images/extra_valid/*.tiff`. Cited artifacts chased:
`results/precision_at_k_14roi_prodseed_chromatin_hembbox_{raw,per_roi}.csv` (untracked) and the
committed baseline `results/precision_at_k_14roi_prodseed_chromatin_halfpixfix_{raw,per_roi}.csv`.
Decisions read: `DECISIONS.md` D1, D5 (with its 2026-09-08 amendment), D7, D8, D9;
`D8_TEMPLATE_ANCHOR.md`. Prior audit reconciled:
`Research Logs/2026-09-16-hembbox-precision-at-k-audit.md`. Audit script:
`tightening_process_audit.py`; tables `results/tightening_process_audit_*.csv` (24 of them).
**Everything below is re-derived from `databases/MIDOG++.json`, the ROI pixels, and the notebook's
own persisted CSV — never from the notebook's printed output, and never by calling
`midog_utils.seed_selection` as the oracle.**

---

## Conflict of interest

`git log -1 --format='%an %ar' -- production_hematoxylin_only/tightening_process_hem_vs_gray.ipynb`
returns nothing: the notebook has never been committed. `git status --porcelain` lists the whole
`production_hematoxylin_only/` directory as `??`, and the three commits immediately preceding the
11:27 run are all `mohini-anand` earlier the same morning (`ad268f0` 10:23, `74fbea0` 10:16 —
the latter touching `midog_utils/find_and_suppress.py`, which the notebook imports for one
constant). So the work under audit and this audit are almost certainly the same hand, hours apart,
and only fresh context separates them. `git status --porcelain midog_utils/` is **clean**.

What limits it: every table here is produced by `tightening_process_audit.py` from
`databases/MIDOG++.json` and the ROI pixels, never from the notebook's output cells. The Otsu
tightening, the agreement/border pool, the RNG draw, the half-open-bbox centre formula and the
template cut are all re-implemented in plain `cv2`/`skimage`/`numpy`; `midog_utils.seed_selection`
is called only inside an explicitly-labelled **consistency** block that is not the Tier A oracle.
All 14 figures were decoded and read, not assumed.

What it cannot cover: the design premises in Part 4 — ROI as the exchangeable unit, `BASE_SIZE` = 51
as the Otsu window, the `min_area`/`max_area_frac`/`min_solidity` gate values themselves. Those go
to `premise-reviewer` and to a reader who did not make them, not to me.

---

## Part 0 — what reproduces

**The engineering is sound. Every number in the CSV, every number rendered into the 14 figures, and
every count in the closing note is the number the code in the working tree produces from these
pixels.** The findings in Part 1 are about a mis-paired quotation, a mechanism, a provenance label
and a plotting colour — not about the run.

| gate / check | result |
|---|---|
| **Execution coherence** | 11 cells, 5 code cells, `execution_count` **1…5, contiguous from 1**; 0 `output_type == 'error'`; 0 unrun code cells; last cell is the markdown closing note, and the last *code* cell (9) ran. Not vendored — imports `midog_utils`, reads `../images/extra_valid` and `../databases/MIDOG++.json`, names this dataset's ROI files. |
| **Provenance — modules** | All six imported modules (`seed_selection`, `template_match`, `channels`, `chromatin`, `dataset`, `find_and_suppress`) are **clean** in `git status` and every one has an mtime **before** the 11:27 run (latest: `find_and_suppress.py` 10:14). No module changed after the artifacts were written, so Tier A here is an **independent reproduction**, not consistency with a stale artifact. |
| **Provenance — artifacts** | The notebook, its CSV and its 14 PNGs are **untracked and never committed** — no commit provenance at all. mtimes are ordered correctly (modules → CSV 11:27:55 → notebook 11:27:57). No `.partial` file is read anywhere. |
| **Composition gate** | 13 checks, **13 pass**. 28 rows = 7 domains × 2 clicks × 2 channels; 28 distinct `(file_name, ann_id, channel)` keys; 0 duplicate rows; 7 distinct ROIs, 7 distinct domains, each domain present exactly 4×; 14 distinct `ann_id`s, matching my independent rebuild exactly; each channel exactly 14×; `base_size` null **iff** `accepted is False`; `reason` null **iff** `accepted is True`; every `base_size` odd and ≥ 5. |
| **Tier A — the CSV** | **112 independently-derived values compared, 0 divergences** (28 rows × `domain`, `accepted`, `base_size`, `reason`), rebuilt from `MIDOG++.json` + pixels. A further **28** values (`which`) are **reproduced by construction, not independently**: my script builds the same two literal labels in the same order the notebook does, so that column proves nothing and is counted separately. `domain` *is* independent — it comes from my own JSON parse, without `dataset.canonical_tumor`. `results/tightening_process_audit_tier_a_divergences.csv` is empty. |
| **Tier A — the selection** | The whole example set re-derived independently: 14 ROIs on disk, 7 domains, the alphabetically-first-file rule reproduces on 7/7 (013 over 094 included), and all 14 `(ann_id, cx, cy)` triples match the notebook's cell-5 table exactly. |
| **Tier B — figure-rendered numbers** | The 14 PNGs carry **140 title strings** (28 panel-rows × 5 panels) holding numbers and verdicts that exist in **no CSV**: the Otsu cut, the component label, area, solidity, the ACCEPTED/REFUSED verdict, the refusal reason, `base_size` and the template's N×N size. All 28 panel-rows were read off the renders (4 figures in full, 10 via title-strip crops) and **every one matches** `results/tightening_process_audit_figure_titles.csv`. The 23 "final template (N×N px)" panels were also checked against the actual crop shape: **23/23** are genuinely N×N, no silent border clipping. |
| **Refusal reasons — the notebook's own assertion cannot check these** | Cell 7's cross-check compares only accept/refuse, bbox and size/centre; `tighten_box_otsu` returns a bare `None`, so it cannot distinguish "click on background" from a downstream gate failure. Recomputed `labels[py, px]` for all 5 refusals: **5/5 are genuinely `center_label == 0`**, not a mislabelled area/solidity rejection. The paired half also holds: on 300/14581, 245/6274 and 013/254 `gray_inverted` accepts at that exact pixel (labels 2, 3, 3; areas 1126, 1258, 582 px). |
| **Config drift** | 7 constants checked three ways (notebook cell 1 ↔ module default ↔ `DECISIONS.md`): `OTSU_WINDOW = 51 == tm.BASE_SIZE`; `BORDER = 36 == FSConfig().patch_size // 2` (`tm.PATCH_SIZE = 73`); `min_area = 50`, `max_area_frac = 0.85`, `min_solidity = 0.5`, `method = "binary"`, `center_tolerance = 0` — all equal to `tighten_box_otsu`'s live signature defaults. **No stale override.** The repo's known `FSConfig` gaps (`tm_method = TM_CCOEFF_NORMED` vs D1's `TM_CCOEFF`; `nms_radius = None`) are both still present but **not applicable** — this notebook runs no template match and no NMS. |
| **Consistency vs the module** | `explore_tightening` reproduces `ss.tighten_box_otsu` on **28/28** accept/refuse decisions and `ss.tightened_template_box` on **28/28** `base_size` values — confirming both my reimplementation and the notebook's own printed "cross-check mismatches: 0". Labelled *consistency*, not independent. |
| **Channel locality** | Both channels verified **pure per-pixel** (`max_abs_diff = 0.0`, byte-identical `uint8` after the patch-local `NORM_MINMAX`), so a 51 px crop is a faithful recompute and the "same pixel, different Otsu cut" reading is a local statement, not an ROI-normalisation artifact. |
| **Production-gap checks the notebook skips** | `build_seed` additionally requires `_patch_readable(roi_shape, tx, ty, patch_size=73)` at the *recentred* point; **0 of 23** accepted examples would fail it, so no accepted example here is one production would have refused for that reason. Largest click→template offset: **14.30 px** (245.tiff ann 6243, `hematoxylin_od`). |
| **Cited sibling numbers** | The sibling's "13/14 ROIs shrink, median ~13 %" reproduces **exactly** from the two raw CSVs: 13 shrink, 1 grow, median −13.3 %. |
| **Cross-validation against two independent production runs** | This notebook's `gray_inverted` `base_size` equals the `half_pix_fix` baseline's recorded `base_size` on **7/7** ROIs (45, 51, 47, 33, 31, 37, 29), and on the **4** pairs both this notebook and the hem-bbox run measure on the same click, the gray→hem percentage is **identical to the decimal** (−43.1, −6.1, −10.8, −13.8). Two separately-executed pipelines agree exactly wherever they overlap. |
| **Claims with nothing persisted underneath** | Two: the 140 figure title strings (handled by Tier B above, from pixels) and the closing note's prose counts (handled by Tier A from the CSV + pixels). **Nothing in this notebook was un-checkable.** `cannot check` appears **zero** times in Part 2. |
| **Gates not applicable** | **Tier C: not run** — no divergence required it, and the whole notebook recomputes in 71 s under my script; the measurements that touch pixels are Part 0's Tier B rows and Findings T2-2/T2-4. **Step 4.2 (treatment/domain confound): N/A as a test, but see Part 4** — this notebook uses **one ROI per domain**, so ROI and stratum are perfectly collinear here and Step 4.1 and 4.2 are the same test. **Step 4.3 (recall triad): N/A** — no recall, precision or detection count is computed. **Step 4.4 (length-matched precision null): N/A** — no precision null is computed. **Step 4.6/4.7 (pre-registration, decision table): N/A** — `grep` over `Research Logs/*preregistration*` finds none for this notebook. **Step 4.10 (in-sample operating-point selection): N/A** — nothing is swept or chosen; the gate parameters are inherited from `tighten_box_otsu`'s defaults and the ROI/annotation selection rules are stated in advance and mechanical. **Step 4.11 (timing): N/A** — no timing claim is made. |

**Family.** `git log --diff-filter=A` returns nothing (the file is untracked), so **no family is
established by git**. Established instead from durable in-file evidence: the notebook's cell-0
markdown names `production_seed_precision_at_k_chromatin_hem_bbox.ipynb` as the thing it is a
"visual companion" to; that notebook sits in the same untracked directory and reads
`results/precision_at_k_14roi_prodseed_chromatin_halfpixfix_*.csv`, written by
`production_seed_precision_at_k/production_seed_precision_at_k_chromatin_half_pix_fix.ipynb`.
`grep -rln "tightening_process"` over every `.py`, `.ipynb` and `.md` returns only the notebook
itself, this audit's script, and the prior audit's out-of-scope note at line 751. So: **one
companion arm, no unreported sibling arms of the same measurement, and this headline is not a
best-of-N selection.**

---

## Part 1 — findings

No Tier 1. Neither of the notebook's two headline effects is overturned: `gray_inverted` accepts
13/14 and `hematoxylin_od` 10/14, and on all 10 pairs both channels accept, `hematoxylin_od`'s box
is strictly smaller. Everything below changes a number, a mechanism, a scope label or a rendering.

### Tier 2

#### T2-1 — "`300.tiff`#2 45->19px is the extreme" mis-pairs two different examples, and the mis-pairing manufactures a shrink larger than any real pair

The closing note reads: *"All 10 examples both channels accepted: `300.tiff`#2 45->19px is the
extreme."* `45` is not `300.tiff`#2's `gray_inverted` size. It is `300.tiff`**#1**'s — annotation
14581, whose `hematoxylin_od` run was **REFUSED** and therefore contributes no pair at all.

| example | ann_id | gray | hem | Δ | Δ% |
|---|---|---|---|---|---|
| `300.tiff` **#1** | 14581 | **45** | *refused* | — | — |
| `300.tiff` **#2** | 14575 | **33** | **19** | −14 | −42.4 % |
| `201.tiff` **#1** | 4457 | **51** | **29** | **−22** | **−43.1 %** |

The notebook's own figure contradicts its prose: `tightening_300_14575.png` is titled
`tightened bbox (base_size=33px)` on the top row and `(base_size=19px)` on the bottom.

Two consequences:
1. The stated pair `45 → 19` is a **−57.8 %** shrink — larger than any pair that actually exists in
   the run. The error invents the extreme it is quoting.
2. `300.tiff`#2 is not the extreme on either reading. `201.tiff`#1 is larger in absolute terms
   (−22 px vs −14 px) *and* in proportional terms (−43.1 % vs −42.4 %). The sibling notebook's own
   `SEED_COMPARE` table independently lists `201.tiff` at 51 → 29, −43.1 %.

**Caveat on my finding.** The −43.1 % vs −42.4 % margin is one 2 px step in an odd-integer size;
"one of the two largest" would be defensible, "the extreme" is not. Table:
`results/tightening_process_audit_pairs.csv`.

**CORRECTION TO BE MADE.** Target: cell 10 (closing note, effect 2) of
`tightening_process_hem_vs_gray.ipynb`. Change *"`300.tiff`#2 45->19px is the extreme"* to
*"`201.tiff`#1 51->29px (−22 px, −43 %) and `300.tiff`#2 33->19px (−14 px, −42 %) are the two
largest"*. Do not simply swap `45` for `33` and keep the word "extreme" — that leaves the wrong
example named.

---

#### T2-2 — The stated mechanism fails on exactly the pairs that drive the headline: on 4 of the 10, `gray_inverted`'s box is pinned to the 51 px Otsu window, not to any nuclear outline

The closing note explains effect 2 as: *"the `gray_inverted` component tends to trace the fuller,
softer nuclear outline, while the `hematoxylin_od` component is confined to the densest, most
heavily-stained chromatin."*

Recomputed geometry of all 23 accepted components
(`results/tightening_process_audit_component_geometry.csv`):

| channel | n accepted | median area / window | `base_size` at the 51 px ceiling | bbox touching the window edge | component ≥ 40 % of the window |
|---|---|---|---|---|---|
| `gray_inverted` | 13 | **0.322** | **4** | **5** | **6** |
| `hematoxylin_od` | 10 | **0.157** | 0 | 1 | **0** |

A `base_size` of 51 means the component's bounding box spans the entire Otsu window in at least one
axis — the value is **right-censored**: it reads "at least 51", not "51". Among the 10
both-accepted pairs, `gray_inverted` sits at that ceiling on **4** — 201#1, 201#2, 402#2 and
245#2 — and **those four are exactly the four largest absolute shrinks in the run**, with a clean
gap to the fifth:

| rank | pair | gray → hem | Δ px | Δ% | gray censored at the window? |
|---|---|---|---|---|---|
| 1 | 201#1 | 51 → 29 | **−22** | −43.1 % | **yes** |
| 2 | 201#2 | 51 → 31 | **−20** | −39.2 % | **yes** |
| 3 | 402#2 | 51 → 33 | **−18** | −35.3 % | **yes** |
| 4 | 245#2 | 51 → 35 | **−16** | −31.4 % | **yes** |
| 5 | 300#2 | 33 → 19 | −14 | −42.4 % | no |
| 6–10 | 529#2, 402#1, 529#1, 459#1, 459#2 | — | −8 … −2 | −19.5 … −4.9 % | no |

No censored pair falls outside the top four, and no uncensored pair falls inside it
(`results/tightening_process_audit_shrink_ranking.csv`). **The magnitude of effect 2 is carried
entirely by the pairs whose `gray_inverted` value is not a measurement.**

Reading it off the renders confirms what the numbers say: in `tightening_201_4457.png`'s top row
the accepted `gray_inverted` component (red, area 1592 px = **61.2 %** of the 2601 px window,
bbox 50 × 51) is a merged blob of adjacent nuclei and stroma filling the frame, sitting just under
the `max_area_frac = 0.85` gate — not "the fuller, softer nuclear outline" of the clicked nucleus.
On those four pairs the change from `gray_inverted` to `hematoxylin_od` is not "same object, tighter
boundary" but **a different object entirely**, and the quoted magnitudes (51→29, 51→31, 51→35,
51→33) are lower bounds on the `gray_inverted` side rather than measurements of it.

This does **not** overturn effect 2. The direction holds on all 10 pairs, and the six uncensored
pairs still shrink (−4.9 % to −42.4 %). It changes what the effect *is* on the pairs with the
largest numbers.

**Caveat on my finding.** "≥ 40 % of the window" and "bbox touches the edge" are my thresholds, not
the repo's; the censoring at `base_size == 51` is not a threshold choice — it is exactly what
`_odd(max(y1−y0, x1−x0))` can return given a 51 px patch.

**CORRECTION TO BE MADE.** Target: cell 10 (closing note, effect 2) of
`tightening_process_hem_vs_gray.ipynb`. Add, after the "fuller, softer nuclear outline" sentence:
*"On 4 of the 10 pairs (201#1, 201#2, 245#2, 402#2) the `gray_inverted` component's bounding box
spans the whole 51 px Otsu window, so its `base_size = 51` is censored at the window and the
quoted shrink there is a lower bound. On those the `gray_inverted` component covers 47–61 % of the
window and is a merged multi-nucleus blob rather than one nuclear outline."* Optionally add an
`area_frac = area / patch.size` column to `SUMMARY` so this is visible in the table, not only in
the figures.

---

#### T2-3 — "the first one is literally the same click `…_hem_bbox.ipynb` used at `seed_index=0`" is true of a notebook it does not name and false on 3 of 7 ROIs of the one it does

Cell 0 claims the `#1` pick is *"literally the same click
`production_seed_precision_at_k_chromatin_hem_bbox.ipynb` used for that ROI at `seed_index=0`"*.
The notebook draws `rng0 = default_rng([0, image_id]); idx0 = rng0.integers(len(pool))` with **no
retry** — cell 4 says so explicitly. `build_seed` / `draw_seed_with_retry` **do** retry: a candidate
the Otsu gate refuses is dropped and redrawn on the same stream.

Recomputed against both production runs' own raw CSVs
(`results/tightening_process_audit_seed_provenance.csv`):

| ROI | `#1` here | hem-bbox actually used | match? | half_pix_fix (gray) actually used | match? |
|---|---|---|---|---|---|
| 300.tiff | 14581 | **14499** (`n_retries = 1`) | ✗ | 14581 | ✓ |
| 201.tiff | 4457 | 4457 | ✓ | 4457 | ✓ |
| 245.tiff | 6274 | **6317** (`n_retries = 1`) | ✗ | 6274 | ✓ |
| 459.tiff | 22441 | 22441 | ✓ | 22441 | ✓ |
| 013.tiff | 254 | **249** (`n_retries = 1`) | ✗ | 254 | ✓ |
| 529.tiff | 24918 | 24918 | ✓ | 24918 | ✓ |
| 402.tiff | 20254 | 20254 | ✓ | 20254 | ✓ |
| | | | **4/7** | | **7/7** |

The three misses are exactly the three ROIs where `hematoxylin_od` refuses `#1` — the notebook's own
effect 1. What `#1` *is* is the first draw on that RNG stream, which is the click the
**gray-tightened `half_pix_fix` baseline** used at `seed_index=0` on **7/7** ROIs — and my
`gray_inverted` `base_size` matches that baseline's recorded `base_size` on **7/7** as well
(45, 51, 47, 33, 31, 37, 29), a clean external cross-validation of this notebook's whole
`gray_inverted` column.

**Caveat on my finding.** This is a provenance-label defect, not an evidentiary one: **neither
headline effect depends on it.** Both are within-example `gray` vs `hem` paired comparisons on
whatever click was drawn, and would be equally valid on an arbitrary set of 14 gated annotations.

**CORRECTION TO BE MADE.** Target: cell 0 (Scope paragraph) of
`tightening_process_hem_vs_gray.ipynb`. Replace *"the first one is literally the same click
`production_seed_precision_at_k_chromatin_hem_bbox.ipynb` used for that ROI at `seed_index=0`"*
with *"the first one is the first draw on that ROI's `seed_index=0` stream — the click the
gray-tightened `half_pix_fix` baseline used on all 7 ROIs, and the click
`…_chromatin_hem_bbox.ipynb` drew first on 7 but kept on only 4, since its `hematoxylin_od` gate
refused 300/14581, 245/6274 and 013/254 and redrew."*

---

#### T2-4 — `#2` is labelled `seed_index=1` but is drawn from a pool one row shorter, and matches a real `seed_index=1` draw on only 3 of 7 ROIs

`SUMMARY['which']` labels the second pick `#2 (seed_index=1)`, and cell 4 says the picks are
*"RNG-seeded exactly like `run_pipeline.py`/the precision notebooks (`[seed_index, image_id]`)"*.
The seed is indeed `default_rng([1, image_id])`, but the notebook draws it from `remaining` — the
pool with `#1` removed — so it calls `rng1.integers(len(pool) - 1)`, while
`run_pipeline.py:77` → `build_seed` calls `rng.integers(len(pool))` on the full border-filtered
pool.

Recomputed (same table): the annotation a real `seed_index=1` draw would pick matches `#2` on
**3/7** ROIs, and differs on 300 (14575 vs **14577**), 201 (4469 vs **4457**), 459 (22446 vs
**22445**) and 013 (243 vs **246**). On 201.tiff a real `seed_index=1` would have re-drawn
`#1`'s own annotation, 4457.

**Caveat on my finding.** Same as T2-3: the effects are within-example paired comparisons, so the
label is wrong but the evidence is not. The draw-without-replacement design is a *reasonable*
choice for a 2-per-ROI visual sample — it is the `seed_index=1` label on it that is not.

**CORRECTION TO BE MADE.** Target: cell 5 (`examples.append(... which='#2 (seed_index=1)')`) and
cell 4's markdown of `tightening_process_hem_vs_gray.ipynb`. Relabel `#2 (seed_index=1)` to
`#2 (second draw, rng [1, image_id] over the pool minus #1)`, and in cell 4 state that `#2` is
drawn **without replacement** and is therefore not the annotation `run_pipeline.py` would pick at
`seed_index=1`.

---

#### T2-5 — The closing note endorses, by citation, a mechanism the cited notebook's own data does not support

The last sentence of effect 2: *"This is the same pattern `…_hem_bbox.ipynb` measured at full scale
(13/14 ROIs shrink, median ~13 %) **and the mechanism proposed there for that run's precision drop:
a smaller, less distinctive template gives `TM_CCOEFF` less to key on.**"*

Three separate problems, in increasing order of weight:

1. **The two medians differ by ~12 points and are presented as one — but the runs agree exactly
   wherever they overlap, and the gap is pure composition.** This notebook's median over its 10
   both-accepted pairs is **−25.4 %**; the sibling's is **−13.3 %** over 14 production seeds. On
   the **4** pairs both runs measure on the same click (201/4457, 459/22441, 529/24918,
   402/20254), the percentages are **identical to the decimal** — −43.1, −6.1, −10.8, −13.8 in
   both (`results/tightening_process_audit_sibling_overlap.csv`). The gap is entirely which pairs
   each median is taken over:

   | subset | n | median |
   |---|---|---|
   | this notebook, `#1` draws only | 4 | **−12.3 %** |
   | sibling, same-seed ROIs only | 11 | **−12.8 %** |
   | sibling, all production seeds | 14 | **−13.3 %** |
   | this notebook, `#2` draws only | 6 | **−33.3 %** |
   | this notebook, all both-accepted pairs | 10 | **−25.4 %** |

   The `#2` draws appear in no production run at all. So retry is not the driver (the sibling's
   same-seed-only median is −12.8 %, essentially its all-14 value) and the two runs do not
   disagree — the medians are simply not comparable, for a stated reason rather than an asserted
   one. Reconciled this way the citation is *sound on the `#1` draws* and the sentence only needs
   to say which set each number covers.
2. **The sibling's 13/14 is not the paired comparison this notebook's 10/10 is.** Its single growing
   ROI is `013.tiff`, which is precisely a *different-click* ROI (254 → 249, `same_seed = False`).
   So this notebook's within-click 10/10 is the cleaner evidence of the two, and reconciling them
   *strengthens* effect 2 while showing the citation is loose.
3. **The cited mechanism is not supported by the cited notebook's own per-ROI data.** Recomputed
   independently (Spearman of Δ`base_size` against Δ`precision_at_K`, 2,000-permutation two-sided
   p, all 14 ROIs, `results/tightening_process_audit_cited_mechanism_spearman.csv`):

   | arm | K = 10 | K = 20 | K = 30 | K = 50 |
   |---|---|---|---|---|
   | `tm_score` | ρ = −0.153 (p 0.58) | −0.162 (0.58) | −0.306 (0.28) | −0.078 (0.78) |

   For "smaller template → worse `tm_score`" the correlation would have to be **positive** (both
   deltas negative together). **0 of 4** cells have that sign and **0 of 4** reach p < 0.05.

**Caveat on my finding.** n = 14; this is a weak null, not a refutation of the physics. And this
notebook is *reporting* a proposal made elsewhere, not making it — the defect is that it reports it
as settled ("the mechanism proposed there") in a sentence whose other clause ("the same pattern")
is presented as verified.

**Reconciliation.** This independently confirms Finding **T1-4** of
`Research Logs/2026-09-16-hembbox-precision-at-k-audit.md`, whose all-14 `tm_score` ρ values
(−0.153, −0.162, −0.307, −0.078) match mine to the third decimal. That audit's correction is
scoped to the sibling's cell 24; this notebook inherits the same sentence and needs the same fix.

**CORRECTION TO BE MADE.** Target: cell 10 (closing note, effect 2, last sentence) of
`tightening_process_hem_vs_gray.ipynb`. Replace *"…and the mechanism proposed there for that run's
precision drop: a smaller, less distinctive template gives `TM_CCOEFF` less to key on"* with
*"…measured at full scale over 14 production seeds (13/14 shrink, median −13.3 % per ROI). The two
runs agree to the decimal on all four pairs they share; the medians differ (−25.4 % here) only
because this notebook's `#2` draws appear in no production run — over the `#1` draws alone the
median here is −12.3 %. The `TM_CCOEFF`-distinctiveness mechanism proposed there is a hypothesis,
not a measured result: Δ`base_size` and Δ`precision` are uncorrelated and wrong-signed in that
run's own per-ROI data at every budget."*

---

### Tier 3

- **T3-1 — The click marker in panel 3 is black-on-black, and the closing note points the reader at
  it for exactly the case where it is invisible.** Cell 10 says the three hem refusals are *"visible
  in panel 2/3 of those figures as a click marker sitting outside every colored blob."* Panel 2's
  marker is `plot(marker='+', color='red', mew=2)` → strokes red on the binary's black background:
  visible, correct. Panel 3's is `plot(marker='+', color='white', ms=13, mew=2, mec='black')`. A
  `'+'` glyph has **no fill**, so `markerfacecolor='white'` is never drawn; the explicit
  `mec='black'` is the stroke, and `label2rgb(labels, bg_label=0)`'s background is **pure black
  (0, 0, 0)**. Verified by reading the render: `tightening_300_14581.png`'s bottom-row panel 3 shows
  six coloured components and **no marker at all**. Affects **5 of 5** refusals — every refusal in
  the run is "click on background", which is precisely the condition that makes the marker
  invisible. Table: `results/tightening_process_audit_figure_render_defects.csv`.
  **CORRECTION TO BE MADE:** cell 3, `plot_channel_row`, panel 3 — change
  `color='white', mec='black'` to `color='white', mec='white'` (or drop `mec` and add a contrasting
  outline, e.g. `path_effects.withStroke(linewidth=3, foreground='black')`), so the click stays
  visible on the `bg_label=0` background. Alternatively pass
  `label2rgb(..., bg_color=(0.25, 0.25, 0.25))`.

- **T3-2 — In all 14 figures the bottom row's three-line panel-3 title is occluded by the top row's
  axes.** The `components (label=N)` first line is hidden behind the top row's images; only
  `area=… solidity=…` and the verdict survive. Read off every figure. The information is not lost
  (it is in the CSV-adjacent recompute) but the bottom row is not the self-describing panel the top
  row is.
  **CORRECTION TO BE MADE:** cell 7 — add vertical headroom for the bottom row, e.g.
  `plt.subplots(2, 5, figsize=(20, 9.6))` plus `fig.tight_layout(h_pad=3.0)`, or move the
  `label=`/`area=`/`solidity=` line into an in-axes annotation instead of a third title line.

- **T3-3 — `area={r.area}px` renders as `1126.0px`.** `regionprops(...).area` returns a float under
  skimage 0.24, and the f-string carries no format spec, so every panel-3 title shows a spurious
  `.0` on an integer pixel count.
  **CORRECTION TO BE MADE:** cell 3, `plot_channel_row` — `area={int(r.area)}px`.

- **T3-4 — "Every number and pass/fail decision here is cross-checked against the real
  `tighten_box_otsu`/`tightened_template_box`" is overstated.** Cell 7's assertion covers
  accept/refuse, the bbox tuple, `base_size` and the centre. It does **not** cover the Otsu cut
  value, the component label, `area`, `solidity` or the refusal `reason` string — five of the eight
  quantities rendered into each panel-3/panel-2 title, and the ones the closing note's effect-1
  mechanism rests on. (`tighten_box_otsu` returns a bare `None`, so the `reason` strings are
  structurally uncheckable that way.) I verified all five independently and **all 28 rows are
  right** — this is a scope overstatement, not a defect.
  **CORRECTION TO BE MADE:** cell 0, last paragraph — change *"Every number and pass/fail decision
  here"* to *"Every accept/refuse decision, bounding box, `base_size` and template centre here"*.

- **T3-5 — `seed_selection.py:89-153` is off by one.** `tighten_box_otsu` spans lines **89–154** in
  the working-tree file. The mirroring claim itself reproduces: 28/28 accept/refuse and 28/28
  `base_size` agree with the real function, and I read the two implementations side by side —
  identical gate order, identical half-open bbox unpacking `(y0, y1, x0, x1)`, identical
  `patch.size` denominator, identical `_odd(..., minimum=5)` and identical `(x0 + x1 − 1) / 2`
  D8 centre.
  **CORRECTION TO BE MADE:** cell 2 — `seed_selection.py:89-154`, or drop the line numbers and cite
  the function name, which does not drift.

- **T3-6 — latent: the template crop uses raw negative-safe slicing.** `plot_channel_row` cuts
  `tpl = rgb[icy - bh : icy - bh + bs, icx - bh : icx - bh + bs]` with no bounds check. A negative
  start index would silently wrap to the far edge of the array rather than raise. It does not bind
  here — all 23 crops are exactly `base_size × base_size`, and `_patch_readable(..., 73)` holds on
  23/23 — but the recentred point can sit up to 25 px from the click while `BORDER` is only 36, so
  a click at the border margin with a 51 px template could in principle reach 50 px out and wrap.
  **CORRECTION TO BE MADE:** cell 3 — guard the template crop with
  `tm.read_padded_patch(rgb, cxr, cyr, bs)` (which returns `None` on an out-of-bounds read) instead
  of raw slicing, and render "template not readable at this centre" when it returns `None`.

---

## Part 2 — verdict per conclusion

| # | cell | claim (verbatim or condensed) | verdict |
|---|---|---|---|
| C1 | 0 | "7 domains -> 7 ROIs; the alphabetically-first file per domain, e.g. `013.tiff` over `094.tiff`" | **reproduces** — 14 ROIs on disk, 7 domains × 2, alphabetical rule holds 7/7 |
| C2 | 0 | "RNG-seeded `[0, image_id]` and `[1, image_id]` … the first one is literally the same click `…_hem_bbox.ipynb` used for that ROI at `seed_index=0`" | **does not reproduce** — 4/7, not 7/7; it is the `half_pix_fix` baseline's click on 7/7 (T2-3) |
| C3 | 0 | the five steps: 51 px `BASE_SIZE` window, `cv2.normalize` → `THRESH_BINARY+OTSU`, `label(connectivity=2)`, gates 50 / 0.85 / 0.5, `method="binary"`, `center_tolerance=0` | **reproduces** — all seven constants match the live module defaults |
| C4 | 0 | "Every number and pass/fail decision here is cross-checked against the real `tighten_box_otsu`/`tightened_template_box`" | **reproduces, overstated** — covers 3 of 8 rendered quantities; the other 5 are right but unasserted (T3-4) |
| C5 | 2 | "Mirrors `seed_selection.py:89-153` line for line" | **reproduces, overstated** — the mirroring is exact (28/28, 28/28); the span is 89–154 (T3-5) |
| C6 | 4 | "RNG-seeded exactly like `run_pipeline.py`/the precision notebooks (`[seed_index, image_id]`)"; `#2` labelled `seed_index=1` | **does not reproduce** for `#2` — drawn over `len(pool)−1`; matches a real `seed_index=1` on 3/7 (T2-4) |
| C7 | 5 | the domain→ROI map and the 14 `(ann_id, cx, cy)` rows | **reproduces** — all 14 rebuilt independently, exact |
| C8 | 7 | "cross-check mismatches …: 0" | **reproduces** — my own 28/28 + 28/28 consistency confirms it |
| C9 | 7 | the 14 figures: 140 panel titles, 23 template crops | **reproduces** — all 140 title strings verified against pixels; 23/23 crops genuinely N×N. Render defects T3-1/T3-2/T3-3 are separate |
| C10 | 9 | "accepted: 23/28"; gray 13/14, hem 10/14; the 28-row table | **reproduces** — 112 independently-derived values, 0 divergences |
| C11 | 10 | "Two distinct effects are visible across these 14 examples, not one" | **reproduces** — the accept/refuse effect and the size effect are separable and independently verified |
| C12 | 10 | "`hematoxylin_od` refuses outright more often. `gray_inverted` accepted 13/14; `hematoxylin_od` accepted only 10/14" | **reproduces, overstated** — the counts are exact; as a directional generalisation it rests on **3** discordant pairs, exact two-sided McNemar **p = 0.25**, and 3 of 7 ROIs carrying any difference (floor 2/2³ = 0.25). Supported claim: *"4 of 14 vs 1 of 14 in this sample; direction not resolvable at this n"* |
| C13 | 10 | three of four hem refusals: gray accepted the same click, hem called it background — "visible in panel 2/3 … as a click marker sitting outside every colored blob" | **reproduces, overstated** — 3-of-4 and the mechanism verified 5/5 from `labels[py, px]`; **panel 3's marker is invisible** (T3-1), so only panel 2 shows it |
| C14 | 10 | "in production (where a refusal triggers a redraw …) this means more retries" | **reproduces** — `build_seed`/`draw_seed_with_retry` do redraw; the sibling recorded `n_retries = 1` on exactly those 3 ROIs |
| C15 | 10 | "Where both channels accept, `hematoxylin_od` gives a smaller box — every time. All 10 examples both channels accepted" | **reproduces** — 10/10 strictly smaller, 0 tied, 0 larger. ROI-level exact sign-flip over **G = 6** ROIs: p = 0.031, at the floor |
| C16 | 10 | "`300.tiff`#2 45->19px is the extreme" | **does not reproduce** — `300.tiff`#2 is 33→19; the extreme is `201.tiff`#1, 51→29 (T2-1) |
| C17 | 10 | "`459.tiff`#1/#2 barely move (33->31, 41->39)" | **reproduces** — exact |
| C18 | 10 | "the `gray_inverted` component tends to trace the fuller, softer nuclear outline, while the `hematoxylin_od` component is confined to the densest … chromatin — a visibly different Otsu cut on the same tissue" | **reproduces, overstated** — the different Otsu cut is real (gray cuts 103–166, hem 28–104 on the same patch), but on 4/10 pairs gray's `base_size` is censored at the 51 px window and its component covers 47–61 % of it (T2-2) |
| C19 | 10 | "the same pattern `…_hem_bbox.ipynb` measured at full scale (13/14 ROIs shrink, median ~13 %) and the mechanism proposed there …" | **reproduces, overstated** — 13/14 and −13.3 % reproduce exactly, and the two runs agree **to the decimal on all 4 pairs they share**; the quoted medians nonetheless cover different click sets (−25.4 % here over 10 pairs, −12.3 % over the `#1` draws alone), the sibling's one grower is its one different-click ROI, and the cited mechanism is null and wrong-signed in the sibling's own data (T2-5) |
| C20 | 10 | "…not a free, behaviourally-neutral simplification — it changes both which clicks pass the gate and how large the resulting template is … both effects point the same direction" | **reproduces** — both effects verified, both in the same direction. Scope: 14 examples from 7 ROIs, 1 per domain; the *magnitudes* carry the T2-2 censoring caveat |

`cannot check`: **0 of 20.**

---

## Part 3 — what was re-run versus read

**Audit script.** `tightening_process_audit.py` at the repo root. Run start to finish, clean, under
`/Users/mohinianand/anaconda3/bin/python3` before this section was written: **exit 0, ~75 s, 24
tables** to `results/tightening_process_audit_*.csv`. It regenerates every table above. Nothing in
this log came from a heredoc that is not in that script, except the PIL crops used to read the
figure renders, which live in the scratchpad.

**Tier A (recomputed, 112 independently-derived values, 0 divergences; a further 28 `which`
labels reproduced by construction and counted apart).** `tightening_process_summary.csv`, rebuilt from
`databases/MIDOG++.json` + the 7 ROI TIFFs. The annotation frame, the `agreement_pool` split, the
`border_filter`, the `[seed_index, image_id]` draws, the Otsu tightening and the D8 centre formula
are all re-implemented in the script; `midog_utils.seed_selection` is called **only** in the
labelled consistency block. The composition, execution-coherence, provenance and config-drift gates
are all Tier A and all in the script.

**Tier B (recomputed from pixels, well inside the ~20 min budget — the whole script, Tier A and Tier B together, is 71 s).** 7 full ROI loads
(`tifffile`, 5412 × 7215 each), one full-ROI `rgb2hed` + one full-ROI `cvtColor` to prove channel
locality, then 28 × (51 px crop → channel → `NORM_MINMAX` → Otsu → `label` → `regionprops`) and 23
template crops. Chosen because the entire notebook is a claim about what happens inside that crop:
a wrong answer there falsifies everything, and nothing else in the notebook touches pixels. The
14 rendered PNGs were decoded and read — 4 in full (`300_14581`, `300_14575`, `201_4457`,
`013_243`) and the remaining 10 via title-strip crops — covering all 28 panel-rows.

**Tier C: not run.** No divergence required it, and my script recomputes the notebook's entire
substance in 71 s, so re-execution would add nothing. The measurements that touch pixels are Part 0's
Tier B rows and Findings T2-2 and T3-1.

**Read, not re-run.** `midog_utils/seed_selection.py` (all of `tighten_box_otsu`,
`tightened_template_box`, `build_seed`, `agreement_pool`, `border_filter`, `_patch_readable`),
`template_match.py` (`read_padded_patch`, `BASE_SIZE`, `PATCH_SIZE`, `build_augmentations`),
`channels.py`, `chromatin.py` (`hematoxylin_od`), `dataset.py` (`image_annotations`, `load_roi`),
`find_and_suppress.py` (`FSConfig`), `production_pipeline/run_pipeline.py` (`select_annotation`),
the sibling notebook's cells 1/3/5/17/18/24, `DECISIONS.md` D1/D5/D7/D8/D9 with
`grep -n '^## D\|^### Amendment'` run first (D5 carries a 2026-09-08 amendment; D8's body is a
pointer to `D8_TEMPLATE_ANCHOR.md`, which is the operative entry and which I read).

**What I looked for beyond Step 4's named modes.** (a) Whether the notebook's reimplementation
diverges from the real function in the *reason* strings, which its own assertion structurally
cannot check — it does not, 5/5 (Part 0). (b) Whether `hematoxylin_od` carries any global
normalisation that would make "the same pixel, a different Otsu cut" an ROI-level rather than a
local statement — it does not; unlike `to_hematoxylin`, which percentile-clips over the whole
image, `to_hematoxylin_od` is `rgb2hed(...)[:,:,0]`, pure per-pixel, verified numerically. (c)
Whether `base_size` can saturate at the Otsu window and censor the comparison — **it can and it
does, on 4 of 10 pairs** (T2-2); this is the finding none of modes 1–11 would have produced. (d)
Whether any displayed template crop is silently clipped at an ROI edge (numpy would return a
smaller array, not raise) — 23/23 are exactly N×N. (e) Whether `build_seed`'s extra
`_patch_readable` check at the recentred point would have refused any example the notebook accepts
— 0/23. (f) Whether `ds.image_annotations`' permissive `category_id=None` default caused the
known category-filter drift — it did not; cell 5 filters explicitly with `== ds.MITOTIC`. (g)
Whether the panel-3 and panel-4 markers land where the annotation database puts them — checked on
four renders against `cx`/`cy`. (h) Whether my own Tier A count was self-proving — it partly was:
28 of the 140 values are the `which` label my script constructs from the same literals the
notebook does, so the headline count is **112**, not 140. (i) Whether the −25.4 % / −13.3 %
median gap in T2-5 reflects disagreement between the two runs — it does not; they are **identical
to the decimal on all 4 pairs they share**, and the gap is composition (T2-5).

**Appendix facts that no longer hold, as of 2026-09-16.**
- **Notebook and figure counts have grown again.** 66 notebooks (appendix: 59) and **238** embedded
  PNGs (appendix: 217). Treat both as order-of-magnitude only.
- **`results/` is back to partially untracked.** 339 `results/*.csv`; `git status --porcelain`
  reports 52 untracked entries, 22 of which this audit wrote — so ~30 pre-existing untracked
  tables, against the appendix's "two days later all 233 were tracked". This is exactly the churn
  the appendix warns about, and exactly why Tier A is keyed to *persisted*, not *committed*.
- **`FSConfig.tm_method` still defaults to `cv2.TM_CCOEFF_NORMED` while D1 selects `TM_CCOEFF`**
  (`find_and_suppress.py:40`), and `nms_radius` still defaults to `None` (line 35). Both still
  true; neither applies to this notebook.
- **`dataset.image_annotations(annotations, file_name, category_id=None)`** still carries the
  permissive default (`dataset.py:122`). Still true.
- **`tm.BASE_SIZE = 51`, `tm.PATCH_SIZE = 73`** — both still true, and the Otsu window and the
  36 px border in this notebook are derived from them correctly.
- **Tier B costs are lower than the appendix's figures suggest for this class of work.** Not a
  find-and-suppress pass: one `load_roi` is ~1.0 s, one full-ROI `rgb2hed` is **4.2 s**, one
  full-ROI `cvtColor` is **0.26 s**, and one 51 px Otsu-tighten is milliseconds. The appendix's
  ~8 s single-seed figure is about the search path and does not bound a seed-gate audit.
- **`production_seed_precision_at_k_chromatin_hem_bbox.ipynb` is no longer at
  `production_seed_precision_at_k/`.** It was moved to `production_hematoxylin_only/` at 11:39 on
  2026-09-16, as the prior audit's own header records. Locate it by filename, not folder — the
  appendix's own lesson.

---

## Part 4 — premises this audit inherited

Each line is a **project decision I applied, not a fact I verified.**

| premise | source | what would falsify it | does a verdict lean on it? |
|---|---|---|---|
| **ROI is the exchangeable unit** | `DECISIONS.md` D5; `Research Logs/2026-09-04-f1-preregistration.md` §8 | Evidence that two annotations in the same ROI are as independent as two in different ROIs — e.g. within-ROI variance in `base_size` delta matching between-ROI variance. Here the two clicks in one ROI share the scanner, the stain batch and the tissue block, so I treated them as clustered. | **Yes.** It is the whole basis for reporting C15 at **G = 6 ROIs** (p = 0.031, at the floor) rather than **n = 10 examples** (p = 0.002). If the ROI is *not* the cluster, C15's evidence is ~16× stronger than I credited. |
| **The Otsu window is `BASE_SIZE` = 51 px, and `base_size` is capped by it** | `template_match.BASE_SIZE`; `seed_selection.tighten_box_otsu(otsu_window=tm.BASE_SIZE)`; D8 | A wider window returning `base_size > 51` on the four censored pairs, which would show the ceiling is an artifact of the window rather than of the tissue. | **Yes, T2-2 directly.** The censoring finding *is* this premise made visible; it is not a challenge to it. |
| **The gate values `min_area=50`, `max_area_frac=0.85`, `min_solidity=0.5`** | `seed_selection.tighten_box_otsu` signature defaults | Any of them being tuned on the same ROIs they are scored on. I did not verify this and could not — there is no held-out split anywhere in this repo. | **Yes, weakly.** They determine which of the 28 rows are refusals, so C12 and C13 rest on them. All 5 refusals are `center_label == 0`, which none of the three values can change, so C12/C13 are actually **robust** to them. |
| **D8: the template's size and centre both come from the accepted component, centred at `(x0+x1−1)/2`** | `DECISIONS.md` D8 → `D8_TEMPLATE_ANCHOR.md` (read through both superseded sections; the 2026-09-10 correction is operative) | A different anchor convention producing different `base_size` values. | **Yes.** The `base_size` column and every figure's panel 4/5 are D8 quantities. C15/C16/C18 all rest on it. |
| **Prior audits are not independent corroboration** | this agent file | — | **Applied.** I re-derived `Research Logs/2026-09-16-hembbox-precision-at-k-audit.md`'s T1-4 Spearman values myself before citing it; my −0.153 / −0.162 / −0.306 / −0.078 match its −0.153 / −0.162 / −0.307 / −0.078. |

**ROI-per-stratum count, measured.** `databases/MIDOG++.json` restricted to the 14 ROIs in
`images/extra_valid/`: **exactly 2 ROIs per `tumor_type`**, across 7 domains (300/301, 201/233,
245/246, 459/460, 013/094, 529/548, 402/403). But **this notebook uses one ROI per domain**, so
within its own 7-ROI analysis set the exchangeable unit is **perfectly collinear with the
stratum** — there is no within-domain replication at all, "cluster by ROI" and "cluster by domain"
are the same test, and Step 4.1 and Step 4.2 are one check, not two. The notebook does not claim
otherwise: cell 0 explicitly calls the one-per-domain rule *"a simple, stated rule, not a
significance claim"*. I note it because any reader tempted to read the 7 domains as 7 independent
replications of a *domain* effect cannot: they are 7 ROIs, each standing in for its domain, with
nothing to separate ROI variance from domain variance.

**Effective units, measured.** For C15 the honest denominator is **G = 6**, not 10 examples and not
7 ROIs: `013.tiff` contributes no both-accepted pair. All 6 move in the same direction, so exact
two-sided sign-flip over 2⁶ = 64 enumerated patterns gives **p = 0.031**, which *is* the floor —
"all six agree" is the strongest statement this design can make, and it cannot be strengthened
without more ROIs. A t interval on the 6 cluster means (5 df) puts the per-ROI mean change at
**−26.7 %, 95 % CI [−42.0 %, −11.4 %]**. For C12 the effective denominator is **G = 3** (only
300, 245 and 013 carry any gray-vs-hem difference), floor 2/2³ = **0.25**, which the observed data
sits exactly on.

**Which verdicts would change if a premise here turned out to be wrong.** Only the *strength*
verdicts, not the *reproduces* ones. If ROI is **not** the exchangeable unit, C15 moves from
"reproduces, at the resolution floor of a 6-cluster design" to "reproduces, p = 0.002" and C12's
overstatement caveat weakens — but 13/14, 10/14, 10/10 and 45≠33 are arithmetic and would not move.
If the 51 px Otsu window were wider, **T2-2 could dissolve entirely**: the four censored
`gray_inverted` values would become real measurements and the mechanism sentence in C18 might be
right after all. That is the single most load-bearing premise in this audit, and it is a one-line
experiment — re-run `tighten_box_otsu` at `otsu_window = 101` on the four affected clicks — which I
did not run because it is outside the notebook's stated scope. **It is worth a premise review.**
Every other verdict above survives any premise in this table being wrong.
