# Re-ranking `mo_pad61`'s 100-peak pool by per-candidate refinement features, against `chromatin_od`: plan

**Status:** Written 2026-09-30; revised the same day, and again after a review. The selection step (Section 6) has
run, as `production_hematoxylin_only/multiotsu_pad_arms/candrefine_rerank_mo_pad61_select.ipynb`, and picked
`CUE = C1 = REF`. The verdict step (Section 7) has not run.

## 0. Scope

- **The 49 ROIs on disk only:** `images/extra_valid/` (14) and `images/extra_valid/testing_set/` (35).
  - There is no held-out test set, so **every result is in-sample**.
  - The four-feature rule this plan builds on was found on these same ROIs. Every conclusion must say so.
- **`CUE` is chosen by leave-one-domain-out on the selection clicks** (Section 6).
- **The fresh verdict clicks** (Section 3) test robustness to a new click draw, not new candidates.
  - A median 88 % (IQR 80-94 %) of one click's pool lies within 7.5 µm of a candidate in the pools of the other 5
    selection clicks on its ROI.
  - For its true positives, the median is 100 %.
- **Baseline: `chromatin_od`**, by choice. The production default ranker is `tm_score` (D5); `tm` enters as the D5
  check (Section 7).
- **Adoption threshold: +3 points of P@10 for `REF` over `od`** (Section 7).
  - The threshold only decides the adoption verdict.
  - **Every contrast, at every K, for every combination, is reported in full whatever the verdict.**
- **"The pool" means the 100-peak pool:** 100 peaks before NMS (D9), so 83-100 candidates after NMS on the selection
  clicks. That is the set `chromatin_od` ranks.

## 1. The question

On the same candidates, does sorting `mo_pad61`'s 100-peak pool by features measured on each candidate's own
multi-Otsu refinement raise precision@10/20/30 over sorting it by `chromatin_od`? It is asked in two steps:
1. Does `REF`, the earlier four-feature rule plus TM rank, beat `od`?
2. Do the new cues add to `REF`? That is `CUE − REF`.

- All arms rank one identical set, so every difference is a ranking effect.
- P@30's ceiling is the pool's content, the same for every arm.

## 2. Held fixed: the `mo_pad61` search

Copy the steps from `MULTIOTSU_PAD_ARMS_PROMPT.md` § "Pipeline, per (click, arm)", for `mo_pad61` only:

| item | value |
|---|---|
| channel | `hematoxylin_od` (`rgb2hed` H) |
| refinement | `refine_template_box` at **commit `c767898`**: **`PAD_PX = 2`**, a fixed 2 px margin, not 0.5 µm |
| template | the refinement's box, padded to S = 61, the pad set to its own mean; centred on the box |
| blank | 51 x 51 at the click (D11) |
| matcher | `TM_CCOEFF` (D1) |
| peaks | `PEAK_MIN_DISTANCE` 7, `DEEP_FLOOR_Z` −1.5, `MAX_PEAKS` 100 before NMS (D9) |
| NMS radius = match radius | `ev.radius_px(mpp)`, 7.5 µm (D7). **Not** the 5.0 µm in older scripts |
| `od` | 51 px, darkest 10 %, unblanked channel, pad 25 |
| matching | `ev.bucket_detections`, the click's own annotation dropped, re-run on every arm's list (D6) |

**Copying the refinement.**
- Copy `refine_template_box` (cell 3) and its constants (cell 1) verbatim from
  `production_hematoxylin_only/multiotsu_pad_arms/multiotsu_template_refinement_14roi.ipynb` at `c767898`.
- **Do not** copy them from the pad-arms, gallery or shape notebooks. Those still carry `PAD_UM = 0.5`.

**Copying the template and the search.** Both notebooks need them.
- Copy them from `production_hematoxylin_only/multiotsu_pad_arms/multiotsu_pad_arms_49roi_6seed.ipynb`, for `mo_pad61`
  only:
  - the helpers in cell 11: `fill_outside`, `template_sha1`, `padded`, `read_window`, `arm_side`,
    `build_mo_templates`, `search_once`, `rank_by`, `score_top_k`, `template_flags`, `selfmatch_mask`,
    `selfmatch_counts`;
  - the per-click steps of the run loop in cell 15: build the template, blank, search, rank, match.
- Do not copy cell 15's per-run bookkeeping:
  - the near-seed, border-band and `ARM_RULE` columns;
  - `GALLERY`;
  - the `ACC` lookup. `bs`, `tcx` and `tcy` come from this notebook's own `refine_template_box` call.
- From its cell 1, copy only these constants: `BLANK_SIZE`, `HALF_BLANK`, `OD_PAD`, `MAX_SIDE`, `EDGE_READ_PAD`,
  `RING_PX`, `ARMS`, `PAD_ARMS`, `AXES`, `AXIS_KEY`, `BUDGETS`, `TEST_K`, `SELF_PX`, `SEED_INDICES`,
  `N_CLICKS_PER_ROI`, `N_BOOT`, `BOOT_SEED`, `SUBSET_DIRS`.
- Copy the imports these need:
  - `hashlib`;
  - from `midog_utils`: `chromatin as cm`, `evaluate as ev` and `production as prod`;
  - `from midog_utils.nms import nms_by_distance`.
- **Never copy** any of the following:
  - its `PAD_UM`, or its refinement constants;
  - `SMOKE_ROIS`, `SUFFIX`, `OUT_STEM`, `OUT` or `FIG_PATHS`. Its `OUT` paths are the recorded CSVs this plan reads,
    and its loop rewrites them;
  - its cell 4, the 0.5 µm refinement.
- The box centre is `(tcx, tcy) = (int(round(center_xy[0])), int(round(center_xy[1])))`, rounded once.

## 3. Two click sets on the 49 ROIs

| set | clicks | pools | role |
|---|---|---|---|
| **selection** | the 294 recorded pad-arms clicks | already saved; no new search | choose `CUE` by leave-one-domain-out (Section 6) |
| **verdict** | 294 fresh clicks, 6 per ROI | new searches | the statistics and the verdict |

**Selection pools.**
- The clicks are `results/precision_at_k_49roi_6seed_fillgeom_postD11_clicks.csv` (294 rows).
- Each click's pool is the first `n_detections_recorded` rows of its list in
  `results/deep_pool_mo_pad61_49roi_6seed_postD11_top1000.csv`, which carries `od`.
- `n_detections_recorded` is in that run's `_per_run.csv`.
- The deep-pool run's own gate already proved those rows are the recorded 100-peak pool.

**Verdict clicks.**
1. Build the annotation pool:
   - `mitotic = ds.image_annotations(annotations, file_name, category_id=1)`;
   - `pool, _ = ss.agreement_pool(mitotic)`;
   - then `ss.border_filter(pool, 36, rgb.shape)`.
2. Split it into **fresh** annotations (not among the 294 selection clicks) and **reused** ones (among them).
3. Shuffle the two parts:
   - sort each part by `ann_id`;
   - seed one generator, `rng = np.random.default_rng([20261001, image_id])`, with `image_id` from the clicks CSV;
   - permute the fresh part with `rng.permutation`, then the reused part with the same `rng`;
   - order the result fresh first, then reused.
4. Walk that order and run `refine_template_box` at each annotation centre. Keep the first 6 it accepts and log every
   refusal. A ROI uses reused annotations only when it has fewer than 6 accepted fresh ones.
5. The click is the annotation centre.

**Fresh supply, measured 2026-09-30** with the `c767898` refinement on all 49 ROIs:
- There are 1,754 fresh annotations, and 1,671 of them are accepted.
- Capped at 6 per ROI, that gives **271 fresh clicks, 23 short of 294**.
- **8 sparse ROIs are short:**

| ROI | domain | accepted fresh clicks |
|---|---|---|
| `201` | lung | 1 |
| `364` | neuroendocrine | 1 |
| `230` | lung | 2 |
| `524` | melanoma | 3 |
| `013` | breast | 4 |
| `339` | mast cell | 4 |
| `425` | soft tissue | 5 |
| `233` | lung | 5 |

- **So the verdict set has 271 fresh + 23 reused clicks.**
  - Every ROI keeps 6 clicks, so ROI weights match the selection set.
  - The sparse ROIs, where precision is hardest, are not under-weighted.
  - A reused click reproduces its selection-set pool exactly, because the pipeline is deterministic. So 23 of 294
    verdict clicks (7.8 %) overlap the set `CUE` was chosen on. The candidate overlap in Section 0 comes on top.
- An independent recount on 2026-09-30 reproduced every number above.
  - It also found all 294 selection clicks in the border-filtered agreement pool.
  - The 2 px refinement accepts all 294 recorded clicks.
- **The click populations differ.**
  - Fresh clicks pass only the refinement's own gate (1,671 of 1,754 accepted).
  - The selection clicks were drawn under the old Otsu joint gate.

**Pre-declared sensitivity: fresh clicks only (271).**
- Re-run every primary, secondary and D5 row with the 23 reused clicks dropped.
- A ROI's value is then its TP@K summed over its fresh clicks. The pooled denominator is K x the number of fresh
  clicks, the same convention the pad-arms run used for refusals.
- ROIs then keep unequal click counts (`201` and `364` keep 1), so the statistics use `pair_stats_clicks` (Section 7).
- Report it next to the primary result.
- **If the two disagree on any verdict, say so, and call that verdict not robust to the overlap.**

## 4. Refining each candidate

For every detection in a pool, run the same `refine_template_box`.
- **Inputs:** the detection's `(cx, cy)`, the **unblanked** channel and the ROI's `mpp`.
- **Matching location:** every arm is matched at the **original detection coordinates**, never at the refined box
  centre.

**Refusals.** The method refuses for four reasons:
- its 51 px window crosses the ROI edge;
- no piece centre lies within 3.75 µm;
- the selection covers less than 3.0 µm²;
- the selection covers more than 85 % of the window.

How they are handled:
- **Primary analysis:** a refused candidate, or one with an empty chromatin mask, is removed from every arm before
  ranking, `od` and `tm` included.
- **Unfiltered baseline:** print the unfiltered `chromatin_od` P@K beside it. Report `od` filtered minus unfiltered as
  a descriptive row.
- **Border effect:** refusals sit mostly at the edge, so removing them also drops border false positives from the
  baseline. Say so.
- **Tabulation:** count refusals by reason and by bucket (mitosis, look-alike, unlabelled).

## 5. Features, computed once per candidate

**No feature code exists yet.** These definitions are the spec for it.

F1 and F3-F9 are computed in the refinement's 51 px window centred on the candidate.
- **Inputs:**
  - `C` is the refinement's `class_map` (undilated; values 0, 1, 2), and `S` its `selection_mask` (dilated).
  - H, E and D come from one `rgb2hed(rgb.astype(np.float32) / 255.0)` call on the full ROI, cast to float32. That is
    the call `chromatin.hematoxylin_od` makes.
  - H must equal the `hem` window bit for bit.
- **Masks:**
  - `M_chr = S & (C == 2)` is the undilated chromatin.
  - `M_nuc` is the union of the components of `label(C >= 1, connectivity=1)` that contain a pixel of `M_chr`.
- **Direction:** "↑" means a larger value is more mitosis-like.
- **Box size:** none of the features reads the box size, so the 2 px margin affects the template only.

| id | name | definition | dir |
|---|---|---|---|
| F0 | `od` | as the pipeline computes it | ↑ |
| T | `tm_rank` | position in the score-ordered pool | low is better |
| F1 | `core_od` | mean H over `M_chr` | ↑ |
| F2 | `not_white` | −(share of the 61 px window with gray > 215); gray is `cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)`, the window cut from a `cv2.BORDER_REPLICATE`-padded copy | ↑ |
| F3 | `solidity` | `regionprops(S.astype(np.uint8))[0].solidity`, the selection as one region | ↓ |
| F4 | `condensed` | `|M_chr| / |M_nuc|` | ↑ |
| F5 | `texture` | std/mean of H over `M_chr` | set per fold (Section 6) |
| F6 | `rim_core` | mean H over the rim `M_nuc & ~E`, divided by mean H over `E`, where `E = cv2.erode(M_nuc.astype(np.uint8), DILATE_KERNEL, iterations=2) > 0` (default border). NaN if `E` or the rim is empty | ↓ |
| F7 | `spiky` | see the note below | ↑ |
| F8 | `e_over_h` | mean E / mean H over `M_chr` | ↓ |
| F9 | `d_over_h` | mean D / mean H over `M_chr` (catches pigment) | ↓ |

- **F7 `spiky`:**
  - Take the largest piece of `label(M_chr, connectivity=1)`. The feature is NaN below 30 px of area.
  - Get its closed outline with `cv2.findContours(piece.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)`.
    Resample it to 128 points at equal arc-length spacing, by linear interpolation along the closed polyline (the
    last-to-first segment included).
  - Take the distances from the outline points to the piece's centroid, minus their mean, and FFT them.
    - The contour points are (x, y). So the centroid is (x, y) = (mean column, mean row) of the piece's pixels.
  - The feature is the energy, `abs(np.fft.rfft(r)) ** 2`, at harmonics 6-64 divided by the energy at harmonics 1-64.
    It is NaN if that denominator is 0.
- **NaN handling:** a NaN feature takes the worst rank: rank with `pd.Series.rank(method='average',
  na_option='bottom')`. Do not use `scipy.stats.rankdata`, which returns all NaN when any input is NaN. Report each
  feature's NaN rate by bucket.
- **Out of scope for this round:** clear zone, nucleolus and paired bars.

## 6. Arms, and choosing `CUE` by leave-one-domain-out on the selection set

All arms rank the same accepted pool.
- **Combining:** equal-weight averages of within-pool ranks, from `pd.Series.rank(method='average',
  na_option='bottom')` in the declared direction (rank 1 = most mitosis-like).
- **Ties:** broken by `od` rank, then by `tm_rank`, using `kind='mergesort'`.

| arm | key | role |
|---|---|---|
| `od` | F0 | **baseline** |
| `tm` | score | the D5 check (Section 7) |
| `REF` | F0, F1, F2, F3, T | the earlier rule plus TM rank, re-measured on the refinement |
| `CUE` | chosen below | **the method** |

**`REF` re-implements the earlier rule.**
- The earlier rule measured its features differently:
  - on the single piece under (or nearest) the centre;
  - with no refusals;
  - with core darkness over the dilated piece.
- Here, F1 and F3 come from `refine_template_box`.
- As a descriptive check, print `REF`'s selection-set P@10. For comparison:
  - the earlier implementation gave about 0.64 on the deep pool's top 100 (a scratchpad analysis, not in the repo);
  - the recorded `mo_pad61` P@10 is 0.563 for `chromatin_od` and 0.542 for `tm_score`.

**On the selection set:**
1. **Visual check:** for `300.tiff`, `201.tiff`, `013.tiff` and `529.tiff`, save montages of the 10 highest- and 10
   lowest-scoring candidates per feature. These are the first ROIs of the mast cell, lung, breast and melanoma domains
   in clicks-CSV order. Change no feature definition here; the checkpoint in Section 9 takes the montages to the user.
2. **Labels:** re-match the accepted pool in `tm` order.
   - A candidate is a mitosis when its bucket is `human_correct_label`; every other bucket is negative.
   - Each (click, candidate) row counts once.
   - Every AUC uses these labels. A NaN value is left out of that feature's AUC.
3. **Score the 5 pre-declared combinations** on the re-matched lists:

| combination | terms |
|---|---|
| `C1` | `REF` |
| `C2` | `REF` + F8 + F9 (colour) |
| `C3` | `REF` + F4 + F5 (condensation) |
| `C4` | `REF` + F6 + F7 (edge) |
| `C5` | `REF` + F4-F9 (all) |

4. **Choose by leave-one-domain-out.** For each of the 7 domains, use only the other 6 domains' ROIs:
   - fix F5's direction by the sign of its AUC − 0.5 over all their candidates pooled; an AUC of exactly 0.5 gives ↑;
   - pick the combination with the highest mean over ROIs of `(P@10 + P@20) / 2`;
   - within 0.5 points (0.005 as a fraction) of the best, prefer fewer terms;
   - among equal term counts, take the higher mean; an exact tie goes to the lower `C` number.

   Record each fold's pick, and the held-out domain's `(P@10 + P@20) / 2` for that pick minus `REF`'s.
5. **Freeze `CUE`:**
   - Apply the same rule, F5's direction included, to all 49 ROIs.
   - `CUE` is that pick if at least 4 of the 7 folds also picked it; otherwise `CUE = REF`.
   - In a smoke run, the folds are the domains present, and the bar is more than half of them.
   - If `CUE = REF`, the secondary family and `CUE − tm` are skipped.
6. **Write the spec:** `production_hematoxylin_only/multiotsu_pad_arms/candrefine_rerank_spec.json`. It holds:
   - the chosen terms, directions and tie-break;
   - F5's direction from all 49 ROIs, whatever `CUE` is. The verdict set scores all 5 combinations with it;
   - each fold's pick and F5 direction;
   - the SHA-1 of every feature function's source, `hashlib.sha1(inspect.getsource(f).encode()).hexdigest()`;
   - the refusal policy, the K values, `BOOT_SEED` and the verdict `DRAW_SEED`.
7. **Report all of it, all descriptive; no verdicts come from this set:**
   - per-feature AUC per ROI;
   - every combination's P@10/20/30 and its deltas against `od` and `REF`;
   - each fold's pick and its held-out delta.

## 7. Verdict set: statistics and decision

**Statistics.** Copy cells 8 and 9 of `production_hematoxylin_only/multiotsu_pad_arms/multiotsu_pad_arms_49roi_6seed.ipynb`
verbatim, with the imports they use (`itertools`, `math`, `fractions.Fraction`, `scipy.stats as sps`). They hold:
- the window sweep's statistics;
- `pair_stats_clicks`, for ROIs with fewer than 6 clicks;
- `stats_for`, which picks between the two.

Keep both of their validations:
- the fill-geometry run's 10 recorded `delta_stats` rows must reproduce;
- `pair_stats_clicks` must equal `pair_stats` at 6 clicks per ROI.

Do not copy the statistics of `shape_vs_mass_template_terms_49roi_3seed.ipynb`, whose `pair_stats` divides by 3 clicks
per ROI. The bootstrap ROI order is the order in which `file_name` first appears in
`results/precision_at_k_49roi_6seed_fillgeom_postD11_per_run.csv`.
- **Unit:** the ROI (49). A ROI's value is its TP@K summed over its clicks.
- **Effect:** the pooled precision delta, in points (100 x `pooled_delta_precision`).
- **Uncertainty:** a cluster bootstrap CI (`N_BOOT` 10 000, `BOOT_SEED` 20260917), a t CI, W/L/T, `n_nonzero`,
  `min_attainable_p` and an exact sign-flip p.
- **Multiplicity:** Holm within each family.
- **D5 bar:** the CI excludes 0 **and** a majority of ROIs move in the delta's direction.

| family | tests | role |
|---|---|---|
| primary | `REF − od`, K = 10/20/30 | does the rule beat `od`? |
| secondary | `CUE − REF`, K = 10/20/30 | do the new cues add to it? Skipped if `CUE = REF` |
| D5 | `od − tm`, `REF − tm`, `CUE − tm`, K = 10/20/30 | the D5 check; `CUE − tm` skipped if `CUE = REF` |

**Always reported, whatever the verdict:**
- every row above, with its delta, both CIs, W/L/T, Holm p and D5 bar;
- `CUE − od`, descriptive;
- the unfiltered `chromatin_od` P@K, and `od` filtered minus unfiltered;
- the fresh-only sensitivity rerun (Section 3);
- all 5 combinations' P@K on the verdict set as well;
- a per-domain table and top-10 composition by bucket;
- ROIs with 20 or fewer mitoses, refusals, and milliseconds per candidate.

**Decision rule, fixed now:**
- **Definitions:**
  - **X beats Y at a given K** when the pooled delta of `X − Y` is positive, its Holm p is below 0.05 and the D5 bar
    holds.
  - **A significant loss** at a given K is a negative delta with Holm p below 0.05.
- **Adopt `REF` as `mo_pad61`'s ranking** only if all three hold:
  - it beats `od` at K = 10;
  - no K shows a significant loss against `od`;
  - the pooled P@10 gain over `od` (the point estimate) is at least **+3 points** (0.03 as a fraction).
- **Upgrade to `CUE`** only if all three hold:
  - `REF` is adopted;
  - `CUE` beats `REF` at K = 10;
  - no K shows a significant loss against `REF`.
- **D5 check:** at K = 10, the adopted ranker minus `tm` must have a positive pooled delta and meet the D5 bar.
  - If `CUE` fails it, fall back to `REF` and apply the check to `REF`.
  - If `REF` fails it, nothing is adopted.
- **Every verdict carries the in-sample label** (Section 0).
- **Secondary rows:** `CUE` was picked on candidates that overlap the verdict pools (Section 0), so their p-values are
  optimistic. Say so beside them.
- **Reading rules:** a CI that includes 0 means no difference was detected. No equivalence is claimed.
- **Production change:** adopting any ranker needs a DECISIONS.md entry first (D5).

## 8. Gates and flags

1. **Config:** assert the Section 2 constants against `midog_utils.production`. Also assert that `PAD_PX == 2` and
   that no `PAD_UM` name exists in the copied cells.
2. **Refinement tripwire.** Run the copied `refine_template_box` on the 294 recorded clicks. Compare `bs`, `tcx`,
   `tcy` and the `mo_pad61` `tpl_sha1` against the `condition == 'mo_pad61'` rows of
   `results/precision_at_k_49roi_6seed_mopad_postD11_templates.csv`, joined on `file_name` and `seed_index`.
   - **Always print the number of clicks that differ, and list them.**
   - **Halt only if more than 100 of the 294 differ.**
   - **Expected: 0.** The recorded `pad_px` is 2 on all 294 clicks, and `c767898` reports no box changes on these
     ROIs. A nonzero count below 100 does not halt, but it means the copy or its inputs changed. Explain it in the
     write-up.
   - **Differing clicks in the selection set:** re-search them with the new template, using Section 2's search code,
     instead of using the saved pool.
   - In a smoke run, this gate covers the smoke ROIs' clicks and only prints; the 100-click halt cannot fire there.
3. **Search parity,** in the verdict notebook before any click is drawn. Re-search the 12 recorded clicks on
   `337.tiff` and `509.tiff` whose boxes match.
   - On both axes, they must match the recorded `mo_pad61` rows: the top 30 (`score` cast to float32 first) and
     `tp_at_10/20/30`.
   - The whole post-NMS pool must also equal the first `n_detections_recorded` rows of `top1000`: `cx`, `cy`, `score`
     after a float32 cast, and `od` at `rtol = 1e-12`.
4. **`rgb2hed` parity:** H equals `hem` on every candidate window.
5. **Pool identity:** every arm's list is a permutation of the same accepted pool.
6. **Re-match:**
   - On the selection set, the unfiltered `od` list, re-matched, must give the recorded `mo_pad61` `chromatin_od`
     `tp_at_10/20/30` (`results/precision_at_k_49roi_6seed_mopad_postD11_per_run.csv`). A click re-searched under gate 2
     is compared with its new search instead.
   - On the verdict set, it must give the search's own TP@K.
7. **Verdict notebook,** full run only:
   - after the selection run, paste the spec file's SHA-1 into the verdict notebook as `SPEC_SHA1`;
   - assert that the file's SHA-1 equals it;
   - assert that each copied feature function's source SHA-1 equals the spec's.

**Flags (report, never halt):**
- refusals by bucket;
- candidates within 30 px of the ROI edge;
- pool sizes before and after refusals (the recorded range is 83-100);
- per verdict click, the pad-arms flags `box_outside_blank` and `mass_weight_nonpositive` (`template_flags`, cell 11);
- self-matches: non-TP detections within 2 px (Chebyshev) of the click or the box centre, per arm, in the top 10 and
  the top 30 (`selfmatch_counts`, cell 11).
  - If any arm has one in its top 10, drop the self-match detections, re-match, take the top K again, and rerun every
    tested contrast.
  - Report that rerun beside the primary result.
- `M_nuc` touching the edge of the 51 px window, with `|M_nuc|`. In dense tissue it can be a tissue-sized component,
  and then F4 and F6 measure tissue density.

**Supply check:** assert the verdict draw reproduces Section 3's numbers: 271 fresh clicks plus 23 reused, with the
same 8 ROIs short. A mismatch means the pool, the exclusion or the refinement changed. Full run only; a smoke run
prints its counts.

## 9. Files, runtime, environment

**Inputs,** every CSV read with `float_precision='round_trip'`:
- clicks: `results/precision_at_k_49roi_6seed_fillgeom_postD11_clicks.csv`;
- selection pools: `results/deep_pool_mo_pad61_49roi_6seed_postD11_{top1000,per_run}.csv`;
- tripwire and parity: `results/precision_at_k_49roi_6seed_mopad_postD11_{templates,per_run,top30}.csv`, the
  `mo_pad61` rows;
- statistics validation: `results/precision_at_k_49roi_6seed_fillgeom_postD11_{per_run,delta_stats}.csv`, read as
  `REF_RUN` and `REF_STATS`. `ROI_ORDER_ALL = list(pd.unique(REF_RUN['file_name']))`, as in pad-arms cell 2;
- ground truth: `images, annotations = ds.load_annotations('../../databases/MIDOG++.json')`, then
  `ds.image_annotations(annotations, file_name)`, with the click's own annotation dropped by `seed_ann_id`;
- ROI images: `{'original_14': '../../images/extra_valid', 'testing_35': '../../images/extra_valid/testing_set'}[subset] + '/' + file_name`,
  read with `ds.load_roi`.

**Files:**
- `production_hematoxylin_only/multiotsu_pad_arms/candrefine_rerank_mo_pad61_select.ipynb`: the selection set; writes
  the spec.
- `production_hematoxylin_only/multiotsu_pad_arms/candrefine_rerank_mo_pad61_verdict.ipynb`: fresh clicks and the
  statistics.
- `results/candrefine_rerank_mo_pad61_{select,verdict}_{clicks,pool,auc,combos,per_run,per_roi,delta_stats,refusals}.csv`.
  From the notebooks' directory, that is `../../results/`.

**Runtime:**
- Selection: no searches unless the tripwire finds differing clicks (expected 0). About 29k candidate refinements
  plus features, around 5-10 min.
- Verdict: 294 searches at about 3.3 s each (the recorded `mo_pad61` median), plus ROI loads and features, around
  25-30 min.

**Order of work:**
1. Run the selection notebook: smoke first, then the full run.
2. **Stop and report to the user:** the montages, the per-feature AUCs, the fold picks and the spec.
   - Change a feature definition only with the user's approval.
   - A change means re-running step 1.
3. Paste the spec's SHA-1 into the verdict notebook (gate 7). Run the verdict notebook: smoke first, then the full run.

**Environment:**
- Interpreter: `/Users/mohinianand/anaconda3/bin/python3`.
- Execute with `/Users/mohinianand/anaconda3/bin/jupyter nbconvert --to notebook --execute --inplace
  --ExecutePreprocessor.timeout=14400 <path>`. A bare `jupyter` resolves to a Python 3.8 install that cannot import
  `cv2`.
- Smoke runs:
  - use `CANDREFINE_SMOKE_ROIS=337.tiff,509.tiff,201.tiff`. `201.tiff` exercises the reused clicks and
    `pair_stats_clicks`;
  - write `_smoke`-suffixed outputs, a `_smoke` spec included. The verdict smoke reads the `_smoke` spec;
  - if the `_smoke` spec gives `CUE = REF`, set `CUE = C5` in the verdict smoke, so that every code path runs;
  - render the montages for whichever of the four Section 6 ROIs are present;
  - skip the supply check and the spec SHA-1 assert, which run on the full run only.
- Kill only your own PID, never by pattern.

**Rules:**
- If a check fails, or a value here contradicts the repo, stop and report. Do not adapt the design to make anything
  pass.
- Do not commit. Do not modify existing files.
- Docstrings and single-line `def`/call headers follow the repo convention.

## 10. Not in this plan

- a held-out test set;
- learned weights or a classifier;
- the clear-zone, nucleolus and paired-bar features;
- any change to the search, the template or the pool size;
- multi-iteration simulation;
- any production change without a DECISIONS.md entry.
