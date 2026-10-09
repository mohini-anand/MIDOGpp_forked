# Per-candidate feature extraction on the adopted A3 search's top 300 (`mo_pad61` + floor clip + border band): agent prompt

**Status:** spec only. Written 2026-10-01, amended twice on 2026-10-02, pruned on 2026-10-08. Nothing here has been
run. It replaces the feature part of `CANDIDATE_REFINEMENT_RERANK_PLAN.md` (Section 5). It does **not** replace that
plan's ranking, selection or verdict parts; those are deferred to a separate, later analysis that will read this
notebook's outputs.

**The 2026-10-08 revision.** The 2026-10-02 text, with its amendment tables, is in git at `0401ae3`.

| what changed | where |
|---|---|
| 31 feature columns dropped: exact functions of kept columns, near-duplicates, and columns with no counting criterion or false-positive mechanism behind them. 36 features remain | Section 5 |
| The nucleus mask `M_nuc` and its three columns dropped: `rim_core_nuc` failed by this spec's own measurements, `condensed_nuc` read size × surrounding tissue, `h_nuc_mean` had no defined property | Sections 4, 5 |
| Checks removed that only re-confirm copied code or hold by construction: the config asserts, the tripwire, channel and `od` parity, label consistency, the exact regression on the earlier features (and the five columns kept only for it), `difflib`, per-function and input hashes, the shuffle, blur and surround tests, the second implementation, the mask bounds report | Sections 2, 3, 7 |
| Float comparisons with the recorded runs (made on the Intel Mac) use an absolute tolerance of 1e-4; positions, ranks and labels must match exactly except for float near-ties | Sections 3, 7 |
| `directions.csv` gains a `role` column and is frozen before the smoke run; texture gets direction −1; circularity removed; `chr_roughness` gets 0 (it fails the lobulation test) | Sections 2, 5.0, 5.2, 8 |
| The AUC rules, every named column set, and the montage scope stated exactly | Sections 6, 8 |
| `chr_concavities_per_um`, `chr_roughness`, `chr_radial_hf` and `chr_edge_sharpness_per_um` constructions pinned down, with checked values | Sections 5.2, 5.8 |

## 0. What this is, and what it is not

**The job.** For every one of the 294 recorded `mo_pad61` clicks on the 49 ROIs:
1. take the 300 highest-scoring post-NMS candidates of a deep search under the adopted A3 settings (Section 3.1);
2. measure the features of Section 5 on each;
3. check that each feature measures what it claims to measure (Sections 2 and 6);
4. store every raw value.

**Not the job.** No ranking arms, no feature combinations, no P@K comparison between rankers, no selection rule, no
verdict, no production change. The downstream analysis decides all of that from your tables. If you find yourself
writing a rank average, stop.

**Two deliverables of equal weight:**
1. the feature tables (Section 8);
2. the evidence that each feature is valid, and the log of what you found wrong with this spec and with your own work
   (Sections 2 and 6).

A table of numbers without the second deliverable is a failed run.

## 1. Why this spec exists: what went wrong last time

The previous round (`candrefine_rerank_mo_pad61_select.ipynb`) implemented its feature spec faithfully, bit for bit,
and the spec was wrong in five places. Every one of them was visible in the notebook's own outputs. Nobody looked.

| feature | what it was meant to measure | what it measured | how it could have been caught |
|---|---|---|---|
| F7 `spiky` | hairy, spiky chromosome edges | inverse elongation: discs scored 0.99, 2:1 ellipses 0.007, because high-harmonic energy was divided by total outline energy and a disc's total energy is pixel noise | three synthetic shapes; the montage showed round nuclei at the top and metaphase plates at the bottom |
| F4 `condensed` | condensed vs diffuse chromatin | whether the nucleus mask ran off the 51 px window, which it did on 47 % of candidates | splitting its AUC by the edge-contact flag the notebook itself computed |
| F6 `rim_core` | a nuclear membrane rim | edge blur of whatever blob the mask was; a blurred disc with no membrane scores about 0.82 | one synthetic blurred disc |
| F5 `texture` | clumped chromatin texture | the coefficient of variation of intensity, identical for a smooth gradient and for the same pixels shuffled | one shuffle test |
| F9 `d_over_h` | brown pigment | nothing consistent: the chromatin mask is selected on strong hematoxylin, so pigment rarely enters it; on the darkest-stained ROIs the ratio rises with darkness and penalises the darkest mitoses | checking the feature on the melanoma ROIs it was built for; a per-ROI sign table |

Two further problems were structural:
- `od` and `core_od` had a rank correlation of 0.94, so one measurement held two votes. This spec keeps one column per
  measurement and reports overlaps in `redundancy.csv`.
- The direction of a feature was read from a pooled AUC, which mixes between-ROI differences into a within-pool
  statistic. F9 read 0.48 pooled and 0.58 within a click. Every AUC here is within click.

**The general failure:** the implementing agent treated the spec as truth and its own job as transcription. Treat this
spec as a hypothesis and your job as measurement. Section 2 says how.

**What the 2026-10-02 analysis adds.** Read its write-up,
`low_precision_rois/2026-10-02-low-precision-rois-mechanism-and-literature.md`, §2–4, before starting.
- **Both ranking axes are darkness detectors.** Per-ROI P@10 is predicted by one count: the unannotated nuclei in the
  whole ROI that are darker than the ROI's 5th-darkest mitosis. Spearman ρ −0.75 at `tm_score`.
- **Nine ROIs fail:** `211`, `242`, `524`, `546`, `462`, `527`, `410`, `245`, `289`. Their top-10 false positives are
  dark objects:
  - black pigment (`211`);
  - small, smooth, dark round nuclei beside lumen (`242`);
  - smooth dark nuclei with halos (`546`, `462`);
  - crowded, hyperchromatic nuclei (`524`, `527`, `245`).

**What that means for this run.**
- The question to ask of each feature is whether it carries information **beyond darkness**. Section 6's
  `auc_dark_stratified` measures that. Report it; do not select on it (Section 0).
- It must be measured where those failing ROIs' mitoses are, which is why the depth is 300.
- It must be measured on the false positives production will actually rank, which is why the candidates come from the
  A3 search (Section 3.1), not from the older unclipped list.
- The literature names the cues that separate mitoses from dark imposters (full map and sources in Section 5.0):
  - no nuclear membrane, and "clearly visible hairy extension of nuclear material" (the criteria of van Diest et al.
    1992 and Baak et al. 1993, as restated by Al-Janabi et al. 2013; summarised by Veta et al. 2016);
  - "short spikes or blunt projections of chromosomes on the surface of aggregated nuclear material", called the most
    reliable feature (VCGP 2021). Veta et al. 2015 found that most of the weaker algorithms' false positives were dark
    objects lacking "protrusions around the edges";
  - imposters are "Smooth contours, circular or oval, hyperchromatic bodies (single or multiple)" (VCGP 2021).

  Sections 5.2, 5.5 and 5.8 target them.

## 2. What you must check, and when

Keep a running log, `adversarial_log.csv`, with columns `item`, `concern`, `test`, `result`, `action`. Every defect you
find in this spec, and every definition you change, gets a row.

**Before any real data:**
1. **Reproduce every unit test and synthetic suite in Section 5.** They are the check that your code matches this
   text. A value that differs from a checked value means your implementation differs from the text: find out why
   before continuing. Where a suite gives an expected ordering without values, write the ordering down before running
   it. Store every result in `synthetic_tests.csv`. A failing test is reported as failing; never tune a constant or a
   tolerance until it passes.
2. **Write `directions.csv`** (Section 8). It records what each feature is supposed to show before any real-data value
   exists. It is frozen from here on: its SHA-1 goes into the spec JSON and is printed by both the smoke run and the
   full run, which must agree. Do not edit it after seeing real-data values; once you have seen them, a direction is
   no longer a hypothesis.

**The rule for disagreement with this spec.** There are two kinds of failure, treated differently:
- **A gate fails** (Section 7): the pipeline does not reproduce. Stop, report, do not adapt.
- **A definition is defective**: it cannot measure what it claims, as a synthetic test or a montage shows. Fix it, log
  it, keep the spec's version as a `spec_`-prefixed column with role `spec` (the corrected column keeps the original's
  `directions.csv` row; the `spec_` column gets a row with direction 0), and carry on. Do not stop the run over a
  definitional defect, and do not silently keep a definition you have shown to be wrong.

If anything in this spec contradicts the repo (`midog_utils`, `DECISIONS.md`), the repo wins, and you log it. One
exception: the floor clip. It was adopted by `tm_border_floor_arm_verdicts.ipynb` (2026-10-01) and is recorded as the
2026-10-02 amendment to D3 in `DECISIONS.md`, but it is not yet in `midog_utils.production`. Its absence from
`midog_utils` is not a contradiction to resolve in the repo's favour; the search code comes from the notebooks named in
Section 3.

**Smoke run.** Run everything on 4 ROIs (`CANDREFINE_SMOKE_ROIS=337.tiff,509.tiff,201.tiff,211.tiff`; `509` has a
brown-pigmented look-alike and `211` black pigment). Measure the seconds per search and the milliseconds per location,
and print the projection for the full run before starting it.

**After the full run:**
3. Produce the validity report (Section 6) and read it. For every feature, write one sentence: does the evidence say it
   measures what Section 5 claims? If not, say what it measures instead. Do not soften this.
4. Open every montage (Section 8) with the Read tool and look at it. For each, write one sentence on whether its
   highest and lowest tiles look like the property the feature claims. The F7 failure was visible in one glance.

## 3. Inputs, the candidate set, and what is held fixed

**The machine.** The recorded runs this spec compares against were made on the Intel Mac. On the current Apple M5 Max,
colour deconvolution differs in the last float bits (a different libm and BLAS), so `hem` differs by up to about 1e-6
on every ROI. Measured 2026-10-08 on the production pipeline: the same detections in the same order, with scores
differing by up to 1.2e-5 and `od` by up to 6e-8. So every float compared with a recorded value uses an **absolute
tolerance of 1e-4**, and positions, ranks and labels must match exactly except for float near-ties (gate B).

**Held fixed, copied verbatim** (nothing here is re-tuned):
- The 49 ROIs: `images/extra_valid/` (14) and `images/extra_valid/testing_set/` (35), 7 per domain.
- **The 294 recorded clicks and their refinement boxes.** Take the `arm == 'A0'` rows of
  `results/tm_border_floor_49roi_6seed_per_run.csv` with `rank_key == 'tm_score'`: `subset, domain, file_name,
  seed_index, seed_ann_id, click_cx, click_cy, tcx, tcy, bs, mpp, match_radius_px`. This is exactly how
  `tm_floor_definition_49roi_6seed.ipynb` cell 3 builds `CLK`. The search uses these recorded boxes; it does not
  re-refine the clicks.
- **Imports:** as in cell 1 of `tm_floor_definition_49roi_6seed.ipynb` (`ch`, `cm`, `ds`, `ev`, `prod`, `tm`,
  `nms_by_distance`), plus `from midog_utils import seed_selection as ss`, which `refine_template_box` uses.
- **`refine_template_box`**, verbatim from cell 3 of
  `production_hematoxylin_only/multiotsu_pad_arms/multiotsu_template_refinement_14roi.ipynb`, with the constants it
  reads from that notebook's cell 1 (`WINDOW`, `CLASSES`, `DILATE_KERNEL`, `DILATE_ITERATIONS`, `MIN_PIECE_UM2`,
  `MIN_BASE_UM`, `CLICK_RADIUS_UM`, `MIN_AREA_UM2`, `MAX_AREA_FRAC`, `PAD_PX`). Extract it at runtime and `exec` it:
  `s = ''.join(nb['cells'][3]['source'])`, then the `ast.FunctionDef` named `refine_template_box` in
  `ast.parse(s).body`, then `ast.get_source_segment(s, node)`. Gate A checks its SHA-1.
- **The search code**, verbatim from `production_hematoxylin_only/multiotsu_pad_arms/tm_floor_definition_49roi_6seed.ipynb`:
  every `def` in cell 5 (including `roi_floors`, `clip_floor`, `band_mask`, `build_mo_templates` and `detect`),
  `search_arm` from cell 11, and the constants in cell 1 that those functions read. Extract the def blocks with `ast`
  and `exec` them. The only new search code is `detect_deep` and `search_arm_deep` (Section 3.1).
- **From the earlier run**, verbatim from cell 8 of
  `production_hematoxylin_only/multiotsu_pad_arms/candrefine_rerank_mo_pad61_select.ipynb`: `refusal_class` (the
  refusal classes of Section 4) and `resample_closed` (the outline resampling of `chr_radial_hf`).
- **Matching:** `ev.bucket_detections` at `ev.radius_px(mpp)` (7.5 µm, D7), with the click's own annotation dropped. It
  matches greedily in rank order (`ev.greedy_match`), so the prefix of a list gets the labels the prefix alone would.
- **Channels:** `hem = ch.to_channel(rgb, 'hematoxylin_od')`, which is `rgb2hed(rgb.astype(np.float32) / 255.0)[:, :,
  0]` as float32. Compute `hed = rgb2hed(rgb.astype(np.float32) / 255.0).astype(np.float32)` once per ROI; E is
  `hed[:, :, 1]`, D is `hed[:, :, 2]`. `gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)`.
- **Per-ROI thresholds:** `fl = roi_floors(hem)`. The floor `b` and `t_nuc` are both `fl['b_otsu']`, the lower
  threshold of a 3-class multi-Otsu over the whole unclipped `hem`; `t_dark` is `fl['b_otsu_upper']`.

**The candidate set.** For each click, take the first 300 rows of its deep A3 list (Section 3.1), in descending score.
`cand_rank` is 0-based. That is exactly 88,200 candidate rows. The unclipped list had about 83,800 distinct locations;
report this run's count.

Do not drop any row for any reason. Refused candidates stay, with their refusal recorded and their features NaN. The
refusal policy is a downstream decision.

**`in_recorded_pool`** is `cand_rank < n_detections_recorded`: the 79 to 100 rows that production's A3 search ranks,
with `MAX_PEAKS` = 100 (D9). Gate B checks them against the recorded A3 pool.

**Why 300.**
- In the five worst ROIs (`211`, `242`, `462`, `524`, `546`), the top 100 holds about 4.2 mitoses per click and the top
  300 about 7.1. Those counts were measured on the *unclipped* deep pool; the clip adds about 1.5 pool mitoses per click.
- At depth 100, validity on exactly the ROIs these features are meant to help would rest on a handful of positives.

**Other inputs:**
- `databases/MIDOG++.json`, read as `_, annotations = ds.load_annotations(path)` (it returns `(images, annotations)`),
  then `ds.image_annotations(annotations, file_name)` with both categories.
- The recorded A3 run, for gate B:
  - `results/tm_floor_definition_49roi_6seed_per_run.csv`, `arm == 'A3'`: `threshold`, `n_detections`;
  - `results/tm_floor_definition_49roi_6seed_pool.csv.gz`, `arm == 'A3'`: the recorded pools.
- The earlier run's pool, for gate C: `results/candrefine_rerank_mo_pad61_select_pool.csv`.

Read every CSV with `float_precision='round_trip'`.

### 3.1 The candidate search: A3, with no peak cap, top 300 kept

**What A3 is.** The recorded `mo_pad61` search with two changes from `tm_border_floor_arm_verdicts.ipynb`
(2026-10-01). They have different status:
- **the floor clip** (arm A2) is the adopted production change, recorded as the 2026-10-02 amendment to D3;
- **the 30 px border band** (arm A1) is a benchmark-harness rule for ROI crops, not a production change. The product
  reads real pixels past any frame edge, so there is nothing for it to remove there. In the harness it makes 28 of
  2,782 annotated mitoses (1.0 %) unreachable.

Everything is exactly as `tm_floor_definition_49roi_6seed.ipynb` runs its `A3` arm:
- **The floor:** `b` (Section 3).
- **The search channel:** `clip_floor(hem, b)`, which is `np.maximum(hem, b)`. Its `channel_edge` is
  `padded(channel, EDGE_READ_PAD)`, as in cell 11.
- **The template:** the `mo_pad61` template from `build_mo_templates`, cut from the clipped channel at the recorded box.
- **The D11 blank:** the 51 px square around the click, which fills with `b`.
- **The score and peaks:** `TM_CCOEFF`; peak threshold median + `DEEP_FLOOR_Z`·MAD; `PEAK_MIN_DISTANCE`.
- **The 30 px border band:** `band_mask(hem.shape, PAD)` with `PAD` = 30, made invalid and its scores floored
  (`BAND_MODE = 'floored'`).
- **NMS** at `ev.radius_px(mpp)`, then ranking by score.
- **`od`:** computed on the unclipped channel padded by `OD_PAD`, as `detect` does. It is `cm.chromatin_density`: the
  mean of the darkest 10 % of the 51 px window around the detection (260 pixels). That is mostly the object's
  darkness, but on an object smaller than 260 px it includes surround, and a darker neighbour in the window can set it.

**Why the clip exists, and what it means here.**
- Under `TM_CCOEFF` the filled template is mostly a mass detector: the score grows with mean H in the template's kept
  box minus mean H in its filled surround, so dark objects beside lumen, glass or pale stroma were over-rewarded.
  Flooring the channel at `b` caps that reward. `b_otsu` is the whole lowest 3-class Otsu class, not glass alone: it
  clips 26.6–82.5 % of an ROI's pixels (mean 55 %).
- On A3 − A1 (the clip on the band-corrected harness): white-adjacent top-10 false positives fell 316 → 240, and P@10
  at `tm_score` rose 3.30 points.
- **Consequence for this spec.** The features read the unclipped channel (Section 4), so the surround columns
  (`h_ring1_mean` to `h_ring3_mean`, `ring_peak_prominence`, the `ring1_frac_*` shares, `e_ring23_mean`,
  `e_over_h_ring23`) carry the white-surround signal the clip removed from the score. That is correct for measurement,
  and not a reason to clip the features. Read every surround column beside its correlation with
  `white_frac_gray215_61` in `redundancy.csv`.

**The one change: no peak cap, then keep 300.**
- `detect_deep` is `detect` with exactly these changes:
  1. `tm.extract_peaks(..., UNCAPPED)` with `UNCAPPED = 10**7`.
  2. `keep = nms_by_distance(centers, scores, nms_radius)` over all peaks, and
     `keep_rec = nms_by_distance(centers[:prod.MAX_PEAKS], scores[:prod.MAX_PEAKS], nms_radius)`. Greedy NMS runs in
     score order, so production's 100-peak pool is the first `len(keep_rec)` rows of the deep list; gate B confirms it.
  3. Truncate to the first `DEPTH = 300` post-NMS rows **before** computing `od`, and compute `od` only on them. `od`
     costs about 6 ms a row, and an uncapped pool is about 15k.
  4. Return `n_peaks` (uncapped), `n_detections_deep` (the full post-NMS count), `n_detections_recorded`
     (`len(keep_rec)`) and `threshold`.
- `search_arm_deep` is `search_arm` with its `detect(...)` call replaced by `detect_deep(...)`. Nothing else changes.
- **Cost:** the unclipped deep pool took 45 min for the 294 searches on the Intel Mac. This Mac ran the production
  pipeline about 2.5 × faster. Use the smoke run's projection.

**There are no border candidates.** The band is invalid in the search, so every candidate has
`edge_px = min(cx, cy, w − 1 − cx, h − 1 − cy) ≥ 30`, and its 51 px and 61 px windows lie inside the ROI. There is no
padded re-run of any candidate. The search's own internal pads (the correlation pad, `OD_PAD`, the template read pad)
are part of the verbatim search and supply no feature pixels.

## 4. Per-candidate pipeline

**Refinement at the candidate.** Run the verbatim `refine_template_box(hem, cx, cy, mpp)` at the detection's own
coordinates, on the **unclipped, unblanked** `hem`. Every feature reads that channel too.
- The clip and the blank belong to the search: they decide which candidates exist.
- Production's own box refinement also reads the unclipped channel.
- The surround features need the real background pixels that the clip would flatten.

Record:
- `accepted`;
- `refusal`: `refusal_class(ref['reason'])`, one of `accepted`, `no piece near centre`, `selection too small`,
  `selection too large` (`roi border` cannot occur, Section 3.1). Also `empty chromatin mask`, which you set, with
  `accepted = False`, when an accepted refinement leaves `M_chr` empty, as the earlier run's `candidate_features` did;
- `mode`, `n_pieces` (`len(ref['selected'])`) and `base_size`.

Features depend only on the location, so compute them once per distinct `(file_name, cx, cy)` and join to clicks.

**Windows.** Each is centred on the detection `(cx, cy)`.
- **51 px:** the refinement window, `ref['hem_patch']`, where the masks live. `H_w` is its unclipped H, as float32;
  `rgb_w` is the same window of the raw `uint8` RGB.
- **61 px:** the window of `white_frac_gray215_61`.
- **101 px:** the window of the pair, context and neighbour columns. It can cross the ROI edge for candidates 30–50 px
  from it (about 1.6 % of rows). It is **truncated** to the pixels inside the ROI, never padded. `(y0, x0)` is its
  origin in ROI pixels, and `M_chr` sits in it at row and column offset `(cy − 25 − y0, cx − 25 − x0)`, which is 25
  only when the window is not truncated. `hem_101` is its unclipped H.

**Masks**, in the 51 px window:
- `M_chr`: `ref['selection_mask'] & (ref['class_map'] == 2)`, the undilated top multi-Otsu class inside the selection.
  Boundary features use this mask, because the refinement's 1 px dilation fills exactly the notches a boundary feature
  measures. Store `chr_touches_window`, true when `M_chr` reaches the window's edge; an outline cut by the window is
  not the object's outline.
- `M_sel`: `ref['selection_mask']`, the dilated selection, as the template uses it.

**`t_nuc` and `t_dark`** are the ROI-wide thresholds of Section 3. `t_nuc` is the search floor, so "above `t_nuc`" means
"untouched by the clip", and "below `t_nuc`" lumps glass, lumen, halos and pale cytoplasm together. About 45 % of a
typical ROI lies above it, cytoplasm included.

**Rings** for the surround features: `ring_k` for k = 1, 2, 3 is the band of pixels whose distance from `M_chr`
(`scipy.ndimage.distance_transform_edt(~M_chr) * mpp`) lies in `((k−1) µm, k µm]`, clipped to the 51 px window.
`M_chr` pixels are at distance 0, so no ring contains them.

**Clipped pixels.** A pixel is clipped when any RGB channel is 0 or 255.
- True clipping occurs on the 3D Histech ROIs: on `249`, `435`, `230` and `460`, 5.7–7.5 % of GT-mitosis `M_chr`
  pixels have G at exactly 0, with about twice as many at 0 as at 1 (a pile-up). `rgb2hed` clamps 0 to 1e-6, so G
  going from 1 to 0 moves H by −0.04 and E by +0.68: those pixels' colour and local contrast are artefacts.
- `M_chr_clip` = `M_chr` minus clipped pixels. The chromatin colour column (5.6) reads it, and is NaN with reason
  `clipped` when |`M_chr_clip`| < `NS_MIN_PX` = 10. The texture columns exclude clipped pixels as Section 5.4 says.
- **Do not exclude pixels merely because a channel is below 20.** On GT mitoses those are G < 20 dense blue-purple
  hematoxylin, mean RGB about (27–36, 9–16, 47–88), not saturation. Excluding them dropped 60 % of a mitosis's `M_chr`
  against 37 % of other candidates' on the recorded pool: a subset chosen by label and darkness.
- Section 5.8 (a)'s `chr_resid_frac` is the one exception and keeps [20, 235] on its own region, because black pigment
  is low in all channels. Its support is therefore label-dependent too (median `chr_ns_share` 0.76 for mitoses against
  0.44 for top-10 false positives on 7 ROIs), so its AUC is read beside `chr_ns_share`.

**Labels.** Match each click's 300-row list in score order with `ev.bucket_detections` against the ROI's ground truth
minus the click's own annotation. Store `bucket`, `matched_ann_id`, `is_mitosis`. Independently of matching, store
`nearest_mitosis_um` and `nearest_lookalike_um`, the distance to the nearest annotation of each category, excluding
the click's own; the downstream analysis needs labels that do not depend on greedy order.

## 5. Features

**Rules that apply to every feature:**
- Store the **raw measurement**. Never negate, rank, standardise or direction-encode a column. Directions are
  hypotheses and live in `directions.csv` (Section 8).
- Physical units where a quantity has one: µm, µm², optical density. Nothing may depend on `mpp` implicitly.
- Nothing may depend on the 51 px window's own intensity histogram except through `M_chr` and `M_sel`, which are the
  production refinement and are kept for that reason.
- NaN where undefined, with the reason in that family's `nan_reason_*` column (Section 8); never a sentinel number.
- Column names are `snake_case`, and every column is described in the notebook in one line with its unit.
- Every column has a role in `directions.csv`: `feature` (the 36 below), `diagnostic` (`chr_touches_window`,
  `chr_ns_share`, `nbr_n`), or `spec` (a definition you replaced, Section 2).

**The 36 features, by section:**
- 5.2: `chr_area_um2`, `chr_n_pieces`, `chr_minor_um`, `chr_elongation`, `chr_solidity`, `chr_roughness`,
  `chr_concavities_per_um`, `chr_hole_frac`, `chr_radial_hf`, `pair_sep_um`, `pair_area_ratio`, `pair_parallel_deg`,
  `pair_axis_deg`;
- 5.3: `h_chr_mean`, `h_ring1_mean`, `h_ring2_mean`, `h_ring3_mean`;
- 5.4: `h_chr_lap`, `h_chr_glcm_contrast`;
- 5.5: `ring_peak_prominence`, `ring1_frac_nuclear`, `ring1_frac_glass`, `ring1_frac_pale`, `white_frac_gray215_61`,
  `e_ring23_mean`, `e_over_h_ring23`;
- 5.6: `e_over_h_chr`, `brown_frac_surround`;
- 5.7: `dark_area_frac_101`, `nuc_area_frac_101`;
- 5.8: `chr_resid_frac`, `chr_tophat_edge`, `chr_tophat_core`, `chr_edge_sharpness_per_um`, `nbr_nuc_area_med_um2`,
  `chr_area_over_nbr`.

**Named column sets** used elsewhere in this spec:
- **boundary columns:** `chr_solidity`, `chr_roughness`, `chr_concavities_per_um`, `chr_radial_hf`, `chr_tophat_edge`,
  `chr_edge_sharpness_per_um`;
- **spike columns** (the lobulation test, 5.2): `chr_radial_hf`, `chr_roughness`, `chr_concavities_per_um`;
- **texture columns:** `h_chr_lap`, `h_chr_glcm_contrast`;
- **101 px columns:** the four `pair_*`, `dark_area_frac_101`, `nuc_area_frac_101`, `nbr_nuc_area_med_um2`,
  `chr_area_over_nbr`;
- **montage features** (Section 8): `chr_radial_hf`, `chr_tophat_edge`, `chr_roughness`, `chr_concavities_per_um`,
  `chr_hole_frac`, `ring_peak_prominence`, `chr_solidity`, `chr_edge_sharpness_per_um`, `h_chr_lap`,
  `h_chr_glcm_contrast`, `chr_tophat_core`, `ring1_frac_pale`, `e_over_h_ring23`, `chr_resid_frac`,
  `brown_frac_surround`. These are the features that claim a visible property, which a number alone cannot confirm.

### 5.0 What pathologists use, and which columns target it

Use this table as the `rationale` of `directions.csv`. A direction given here or in the section that defines a column
is the hypothesis to write down; a 0 stays 0 unless you log a reason before the smoke run. Quotes are verbatim from the
sources listed after the table.

| # | criterion | evidence | columns | direction |
|---|---|---|---|---|
| C1 | **Short spikes or rods on the surface of the chromatin**, the most reliable cue | VCGP 2021: "The most reliable microscopic feature to identify MF and AMF are the short spikes or blunt projections of chromosomes on the surface of aggregated nuclear material"; also "dark purple staining nuclear aggregates with an irregular contour (short rods or spikes of chromosomes, described by van Diest as hairy projections)". Meuten 2021: "short 'rods' of chromosomes protruding from the surface of aggregates of nuclear material". Al-Janabi 2013 (restating the van Diest and Baak criteria): "clearly visible hairy extension of nuclear material". Veta 2015, on weaker algorithms' false positives: "dark objects that lack other characteristics of mitotic figures such as protrusions around the edges". VCGP QRG 2022, under "Differentiate MF from Mitotic-Like-Figures (MLF)": "Surface contour: MLF are smooth, MF have projecting rods/spikes" | `chr_radial_hf`, `chr_tophat_edge`, `chr_concavities_per_um`; `chr_roughness` fails the lobulation test (5.2) | +1; `chr_roughness` 0 |
| C2 | **No nuclear membrane** | Al-Janabi 2013: "absence of the nuclear membrane". VCGP 2021: "lack of a nuclear membrane". Veta 2015: "hyperchromatic objects that lack a clear nuclear membrane" | `ring_peak_prominence`, `chr_hole_frac`; Section 5.5 says what each can see | 0 for both (5.5) |
| C3 | **Dark chromatin**, necessary but not sufficient | VCGP 2021: "dark purple staining nuclear aggregates". Veta 2015: "many other cellular components can have a similar hyperchromatic appearance" | `h_chr_mean` | +1; it is why `auc_dark_stratified` exists |
| C4 | **Configuration**: a clot or ring, a plate, separate clots | VCGP 2021: "Prometaphase: A dark cluster or ring of chromosomes with variably protruding chromosome rods and spikes"; "Metaphase: An equatorial plate (linear or ring-shaped aggregate with short rods or spikes)"; "Anaphase: Two separate aggregates, variable distances apart, linear with projections/spikes"; "Anaphase and telophase MF are counted as one". Al-Janabi 2013: "either clotted (beginning metaphase), in plane (metaphase/anaphase), or in separate clots (telophase)". Meuten 2021: identifying the phase is "not necessary" | `chr_elongation`, `chr_minor_um`, `chr_n_pieces`, `chr_hole_frac`, the `pair_*` columns (5.2) | 0 for all (below, and 5.2) |
| C5 | **Imposters are smooth, round or oval dark bodies, one or several** | VCGP 2021, mitotic-like figures: "Smooth contours, circular or oval, hyperchromatic bodies (single or multiple)". ICCR, in its own words (its method source is Elston & Ellis 1991): "hyperchromatic, karyorrhectic, or apoptotic nuclei are excluded". Cree 2021: "Strict criteria should be applied for accepting a mitotic figure in distinction from a pyknotic or apoptotic nucleus" | `chr_solidity`, `chr_edge_sharpness_per_um` | −1 for solidity; 0 for edge sharpness (5.8 c) |
| C6 | **Prophase is not counted on H&E** | VCGP 2021: "MF are counted from prometaphase through telophase … Prophase cannot be reliably identified in histological preparations". Lin 2024: dark chromosome aggregates "are most likely in the final stages of prophase, commonly referred to as prometaphase". Meyer 2005 recommends "that prophases be excluded from mitotic counts performed for cancer grading", though its own study counted them. MIDOG++ delegates its criteria to Donovan et al. 2021 and Ibrahim et al. 2022. MIDOG 2022's PHH3-assisted ground truth (not the H&E labels) did not annotate immunopositive cells lacking H&E morphology, "mostly early prophase MFs", and its largest subcategory of figures found only with PHH3 (within "borderline morphology", 87.58 %) was "prophase with resemblance to apoptosis" (31.23 %) | none directly; context for `ring_peak_prominence` (5.5) | — |
| C7 | **Named imposters**: immune cells, stromal cells, two cells in close proximity, karyorrhectic debris | Lin 2024, on which mimic each true mitotic figure resembled (most often the low-agreement ones), not on false positives: "karyorrhectic debris, immune cells, stromal cells, and two cells in close proximity"; "lymphocytes were most likely to resemble mitosis in prometaphase, and stromal cells were likely to mimic mitosis in metaphase" | `chr_area_um2`, `chr_area_over_nbr` (5.8 d); the lobulation test (5.2) | 0 |
| C8 | **Pigment** | Aubreville 2023 (MIDOG++): melanoma "pigment particles that contribute additional imposter structures ... although they have a different chromaticity (brown)" | `chr_resid_frac`, `brown_frac_surround` | −1 |
| C9 | **Cytoplasm and surround** | Ikpatt 2002, restating the criteria of Baak & Oort: "The cytoplasm of the mitotic cells was often larger during mitosis than in the resting cells". No verified source states a clear zone around mitotic figures in H&E. A clear halo is described around apoptotic hepatocytes in rodent liver (NTP atlas). VCGP 2021, mitotic-like figures: "MLF cytoplasm tends to be eosinophilic, however color is not reliable". VCGP QRG 2022: "Cytoplasm color (not reliable). MF tend to have more amphophilic/basophilic cytoplasm (increased RNA), MLF tend to have more eosinophilic cytoplasm" | lightness and H: the `h_ring*_mean` columns, the `ring1_frac_*` shares and `white_frac_gray215_61` (5.3, 5.5). Cytoplasm colour: `e_ring23_mean`, `e_over_h_ring23` (5.5) | 0 |
| C10 | **Context** | VCGP 2021: "A structure within a region of necrosis ... is more likely to be MLF"; "A structure within a region of high mitotic activity is more likely to be a MF". VCGP QRG 2022: "Context: In areas of poor cell preservation, MLF may be more likely" | none; `dark_area_frac_101` and `nuc_area_frac_101` measure crowding, the mechanism of the crowded failing ROIs, not this criterion | 0 for the crowding columns |

**Elongation.** A side-on metaphase plate is elongated, but so are spindle-cell and stromal nuclei, which Lin 2024 names
as metaphase mimics; `462` and `410` are sarcomas. Pathologists need elongation together with surface spikes, so
`chr_elongation` and `chr_minor_um` have direction 0.

**Texture and colour.** VCGP 2021's list of mitotic-like figures (item 5) describes hyperchromatic resting nuclei as
having "more heterochromatin with a more clumped chromatin pattern", heterochromatin "staining as dark granules". C5's
direction comes from the same list's description of imposters, so the texture columns (5.4) and `chr_tophat_core`
(5.8 b) get −1 on the same footing: coarser internal structure is more imposter-like. Record the confound in
`directions.csv`: a prometaphase clump has internal structure from its chromosome arms. No source gives a chromatin
colour criterion (VCGP's "color is not reliable" is about imposter *cytoplasm*), so `e_over_h_chr` has direction 0.
Pigment (C8) is the exception.

**Independence.** VCGP 2021, Meuten 2021, Donovan et al. 2021, MIDOG 2022 and MIDOG++ share authors; the VCGP QRG
2022 names no authors and cites the first three. MIDOG++ also includes van Diest and Veta, who co-authored Al-Janabi
2013 (van Diest), Veta 2015 (van Diest) and Veta 2016. So every C1 and C2 source shares an author with the dataset's
group; Lin 2024 (C6, C7) does not. Their agreement is one school's view, not independent corroboration.

**Sources** (read in full unless marked):
- **Al-Janabi 2013:** Al-Janabi S, …, van Diest PJ, Jiwa M. *PLoS ONE* 8:e82576.
  doi:10.1371/journal.pone.0082576. It restates, without quoting, the criteria "as defined before", citing **van Diest
  PJ et al. 1992**, *Hum Pathol* 23:603–607, and Baak et al. 1993; neither was read. van Diest is a co-author.
- **VCGP 2021:** Veterinary Cancer Guidelines and Protocols, "Morphologies of Mitotic Figures (MF)" v1.0 (Donovan,
  Bertram, Moore et al., May 2021), https://vcgp.org/documents/2021/05/morphologies-of-mitotic-figures.pdf. This is not
  Donovan et al. 2021, *Vet Pathol* 58:243, which was read as an abstract only and is cited for nothing.
- **Meuten 2021:** Meuten DJ, et al. *Vet Pathol* 58:766–794. doi:10.1177/03009858211013712. Full text.
- **Veta 2015:** *Med Image Anal* 20:237–248. doi:10.1016/j.media.2014.11.010. Read as arXiv:1411.5825.
  **Veta 2016:** *PLoS ONE* 11:e0161286.
- **VCGP QRG 2022:** Veterinary Cancer Guidelines and Protocols, "Mitotic Count and Mitotic Figure Morphology Quick
  Reference Guide (QRG)" v1.0,
  https://vcgp.org/documents/2022/03/mitotic-figure-mitotic-count-quick-reference-guide.pdf. Full text.
- **ICCR:** International Collaboration on Cancer Reporting, invasive breast carcinoma dataset, "Histological tumour
  grade" note. It cites Elston CW, Ellis IO 1991, *Histopathology* 19:403–410, which was not read.
- **Cree 2021:** Cree IA, et al. *Mod Pathol* 34:1651–1657. doi:10.1038/s41379-021-00825-7.
- **Lin 2024:** Lin S, et al. The 1000 Mitoses Project. *Int J Surg Pathol* 32:1449–1458.
  doi:10.1177/10668969241234321.
- **Meyer 2005:** Meyer JS, et al. *Mod Pathol* 18:1067–1078. doi:10.1038/modpathol.3800388.
- **Ikpatt 2002:** Ikpatt OF, Kuopio T, Collan Y. *Mod Pathol* 15:783–789. doi:10.1097/01.MP.0000021764.03552.BD.
  Its mitosis criteria are "those described by Baak and Oort", restated.
- **MIDOG 2022:** Aubreville M, et al. *Med Image Anal* 94:103155 (2024). doi:10.1016/j.media.2024.103155. Read as
  arXiv 2309.15589.
- **Aubreville 2023 (MIDOG++):** *Sci Data* 10:484. doi:10.1038/s41597-023-02327-4.
- **NTP atlas:** NTP Nonneoplastic Lesion Atlas, "Liver, Hepatocyte – Apoptosis". Rodent tissue.
- **Cisneros-Soberanis 2024:** *J Cell Biol* 223:e202403165. doi:10.1083/jcb.202403165. By serial block-face EM of
  cultured RPE1 cells, a single chromatid is "0.68 ± 0.04 μm" wide and a cohesed chromosome "1.23 ± 0.1 μm" (used in
  5.8 b). It bounds the physical width of a chromosome arm; it says nothing about what a scan resolves.

**Not verified, so used for nothing:** van Diest 1992's own text, including whether "two parallel, clearly separate
chromosome clots" count as one figure or two (open sources say both; VCGP 2021 counts anaphase and telophase as one
and notes "original reports indicated telophase should be counted as two"); Elston & Ellis 1991's own text; Donovan
et al. 2021 beyond its abstract; Ibrahim et al. 2022, *J Clin Pathol*; any published perimitotic clear zone in H&E.

**Not measured:** the nucleolus. A visible nucleolus marks a resting nucleus, but a nucleolus detector is out of scope
(Section 10).

### 5.1 Identity and pipeline columns

`subset, domain, file_name, image_id, mpp, seed_index, seed_ann_id, click_cx, click_cy, cand_rank, in_recorded_pool,
n_detections_recorded, cx, cy, score, od, dist_to_click_um, edge_px, accepted, refusal, mode, n_pieces, base_size,
bucket, matched_ann_id, is_mitosis, nearest_mitosis_um, nearest_lookalike_um, scanner, roi_group`.

- `dist_to_click_um` is the distance from `(cx, cy)` to `(click_cx, click_cy)`, times `mpp`.
- `n_pieces` is the number of refinement pieces the selection unioned, distinct from `chr_n_pieces` in Section 5.2.
- **`scanner`** is read from the `scanner` column of `images/extra_valid/MANIFEST.md` and
  `images/extra_valid/testing_set/MANIFEST.md`, mapped `3D Histech` → `3dhistech`, `Aperio CS2` → `aperio_cs2`,
  `Hamamatsu S360` → `hamamatsu_s360`, `Hamamatsu XR` → `hamamatsu_xr`. Assert that all 49 ROIs get one.
- **`roi_group`** is `failing9` for `211`, `242`, `524`, `546`, `462`, `527`, `410`, `245` and `289`, and `other40`
  otherwise. The nine are the ROIs below P@10 0.35 on both A3 axes (`low_precision_rois/`). The group was defined after
  seeing results, so every split by it is descriptive.

### 5.2 Geometry of the chromatin, on the undilated `M_chr`

**The piece.** Columns are computed on the largest 4-connected piece of `M_chr` (`skimage.measure.label(M_chr,
connectivity=1)`; a tie goes to the lower label) unless marked *whole*. They are NaN, with reason `small piece`, when
that piece has fewer than 30 px.

**The filled piece and its outline.** `chr_minor_um` through `chr_concavities_per_um`, and `chr_radial_hf`, read the
piece with its holes filled (`scipy.ndimage.binary_fill_holes`). Pathologists judge the outer surface contour (C1,
C5), and on a ring-shaped piece `regionprops` counts the inner boundary too: the vesicular nucleus of Section 5.5,
whose `M_chr` is its rim, scored solidity 0.21 and roughness 1.82 unfilled, the most "mitotic" values, against 0.99
and 1.00 filled. `chr_area_um2` and `chr_hole_frac` read the piece unfilled; the hole is `chr_hole_frac`'s job. The
**outline** is the longest contour of `cv2.findContours(filled.astype(np.uint8), cv2.RETR_EXTERNAL,
cv2.CHAIN_APPROX_NONE)`, as (x, y) points.

| column | definition | unit |
|---|---|---|
| `chr_area_um2` (*whole*) | pixel count of `M_chr` × mpp² | µm² |
| `chr_n_pieces` (*whole*) | 4-connected pieces of `M_chr` | count |
| `chr_minor_um` | `regionprops` `axis_minor_length` × mpp | µm |
| `chr_elongation` | `axis_major_length / axis_minor_length` | ratio |
| `chr_solidity` | `regionprops` `solidity` | ratio |
| `chr_roughness` | `regionprops(filled)[0].perimeter / regionprops(convex)[0].perimeter`, where `convex` is `regionprops(filled)[0].image_convex`. Both perimeters use skimage's estimator, so a disc scores 1.000 | ratio |
| `chr_concavities_per_um` | the number of convexity defects deeper than `CONCAVITY_MIN_UM` = 0.25 µm, divided by the outline length `cv2.arcLength(c, True)` × mpp. With `c` the outline as an `int32` (n, 1, 2) array: `d = cv2.convexityDefects(c, cv2.convexHull(c, returnPoints=False))`; a defect counts when `d[:, 0, 3] / 256 * mpp > 0.25` (the depth is fixed-point); `d is None` counts 0. If OpenCV raises (it can on an outline that touches itself), NaN with reason `convexity error`, and report the count. 0.25 µm is about 1 px, at the pixel staircase | per µm |
| `chr_hole_frac` | 1 − area / area with holes filled, on the largest piece. A ring-shaped `M_chr` scores high; Section 5.5 says why that matters | ratio |
| `chr_radial_hf` | `o = resample_closed(outline, 128)`; `r` = the distance of each point of `o` from the pixel centroid of the filled piece (x the mean column, y the mean row); `e = np.abs(np.fft.rfft(r − r.mean()))**2`; `e[6:65].sum() / (128 * r.mean()**2)`, so harmonics 6 to 64 inclusive. A `skimage` contour gives different values (disc 0.047 against 0.057 at r = 8), so use this one | ratio |

`chr_radial_hf` replaces the earlier `spiky`. Its denominator is the mean radius squared, not the total outline energy,
so a disc scores near zero instead of near one.

**Pair columns: anaphase and telophase (C4), read from the 101 px window, not from `M_chr`.**
- **Why not `M_chr`.** `M_chr` is restricted to the refinement's selection. In base mode that is the one piece under
  the detection; in fallback mode it is the pieces centred within 3.75 µm. So a pair seen through `M_chr` depends on
  where the detection landed: on a synthetic anaphase pair 4.5 µm apart, the two-piece separation of `M_chr` was 4.5 µm
  with the detection in the gap and undefined with it on a plate, and it was 0.9 µm on an interphase nucleus with a
  nucleolus.
- **Why 101 px.** A 51 px window centred on one plate does not contain a partner much more than 6 µm away.
- **Pieces `P`:** 4-connected components of `hem_101 > t_dark`, each at least `PAIR_MIN_UM2` = 1.0 µm² (the
  refinement's speck floor). Thresholding at the ROI-wide `t_dark` keeps the columns off the window's own histogram.
- **Anchor:** the piece of `P` that overlaps `M_chr` (placed as in Section 4) the most; ties go to the larger piece.
  All four columns NaN, reason `no anchor`, if no piece overlaps `M_chr`.
- **Partner:** among the other pieces with area ≥ `PAIR_MIN_RATIO` = 0.5 × the anchor's and centroid distance ≤
  `PAIR_MAX_UM` = 10 µm, the nearest by centroid. All four columns NaN, reason `no partner`, if none.
- **Columns:**
  - `pair_sep_um`: centroid distance;
  - `pair_area_ratio`: smaller area over larger;
  - `pair_parallel_deg`: the acute angle between the two pieces' major axes. With θ = `regionprops` `orientation`,
    a piece's major axis in (row, column) coordinates is `(cos θ, sin θ)`; the centroid line is
    `centroid_partner − centroid_anchor` in the same (row, column) order;
  - `pair_axis_deg`: the mean, over the two pieces, of the acute angle between the piece's major axis and the line
    joining the centroids. Separating plates seen side-on score near 90°.
  - The two angle columns alone are NaN, reason `round piece`, when either piece's major/minor axis ratio is below
    `PAIR_ELONG_MIN` = 1.2, where orientation is undefined.
- **What the columns cannot do.** A neighbouring dark nucleus of similar size is a partner too, so in crowded tissue
  most candidates have one. And side-by-side spindle nuclei in a fascicle give the same angles as separating plates
  (parallel near 0°, axis near 90°); `462` and `410` are sarcomas. So every pair column has direction 0; the
  hypothesis that `pair_axis_deg` near 90° marks anaphase goes in the rationale only.
- **Pair test, mandatory.** Build the synthetic anaphase windows of the construction below and run the verbatim
  refinement on them.
  - The pair columns must agree, within 0.25 µm and 5°, between a detection at the gap's midpoint and one at a plate's
    centre, at 4.5 and 7 µm centre-to-centre. At 9 µm only the plate-centred detection is tested: the verbatim
    refinement refuses the gap-centred one (`no piece centre within 3.75um of the click`), which is its own rule, not
    a pair-column failure.
  - **Oblique pairs.** The same pairs rotated by 30°, 60° and −30° about the canvas centre must give `pair_axis_deg`
    at least 85° and `pair_parallel_deg` at most 5°. An axis-aligned pair cannot catch a sign error in the orientation
    convention; at a 30° rotation a flipped sign gives about 30° instead of 90°.
  - The interphase nucleus with a nucleolus, the uniform dark disc and the metaphase plate must give `no partner`.
  - Checked values, from two independent implementations, identical for noise seeds 0 to 3: 4.5 / 7.0 / 9.0 µm, area
    ratio 1.00, parallel 0°, axis 90°, the same at both centrings; `no partner` on the other three. Rotated (one
    implementation): separation 4.49–4.50 / 7.02 µm, parallel 0°, axis 89.4–90.0°.
  - Construction: `mpp` 0.25; a 201 × 201 canvas at 0.03 OD plus Gaussian noise of σ 0.003 (any seed); then, drawn
    over it with constant values that replace the noise, a cell disc of radius 30 px at 0.045 OD and two plates, each
    7 × 29 px (`|x′ − x_k| ≤ 3`, `|y′| ≤ 14` in coordinates rotated about the centre) at 0.20 OD, centred 18, 28 or
    36 px apart along x′; then `scipy.ndimage.gaussian_filter` with σ 1 px, cast to float32. `t_dark` is the upper
    threshold of a 3-class multi-Otsu over the whole canvas. The other three windows are the membrane cases of
    Section 5.5, and a 9 × 33 px plate at 0.20 OD inside a clear disc of radius 22 px at 0.02 OD inside a cell disc of
    radius 34 px at 0.045 OD, drawn the same way.

**Synthetic suite, mandatory, for the outline columns.** Rasterise at two scales (radius r = 8 px and 16 px) each of
the shapes below, centred on a (6r + 1)² canvas, with x and y the centred coordinates, ρ = √(x² + y²) and
θ = atan2(y, x). The shape is `M_chr`.
- a disc, ρ ≤ r; a 2:1 ellipse with semi-axes r√2 (along x) and r/√2; a 4:1 ellipse with semi-axes 2r and r/2;
- 5-, 8- and 12-lobed stars, ρ ≤ r·(1 + 0.25·cos(nθ));
- a 2:1 ellipse with 12 small spikes: the 2:1 ellipse plus, at each angle φ = 2πk/12, the rectangle
  `e(φ) − 1 ≤ t ≤ e(φ) + 4`, `|o| ≤ 1`, where e(φ) is the ellipse's radius at φ and t, o are the along- and
  across-ray coordinates. The spikes are a fixed 2 × 5 px, not scaled with r;
- two touching discs, each of radius 0.9r, centred at x = ±0.9r;
- a disc with one deep notch: the disc minus the slot x ≥ 0, |y| ≤ max(1, r/8).

Checked values (one implementation, r = 8 / r = 16; for `chr_radial_hf` an independent implementation reproduced the
disc and 12-lobed star to the third decimal, and the 2026-10-08 implementation of this text reproduced every value):

| shape | `chr_radial_hf` | `chr_roughness` | `chr_concavities_per_um` |
|---|---|---|---|
| disc | 0.057 / 0.010 | 1.000 / 1.000 | 0.000 / 0.000 |
| 2:1 ellipse | 0.033 / 0.008 | 1.000 / 1.000 | 0.000 / 0.000 |
| 4:1 ellipse | 0.097 / 0.081 | 1.000 / 1.000 | 0.000 / 0.000 |
| 5-lobed star | 0.024 / 0.020 | 1.116 / 1.121 | 0.303 / 0.146 |
| 8-lobed star | 2.328 / 1.681 | 1.414 / 1.328 | 0.343 / 0.183 |
| 12-lobed star | 1.089 / 1.635 | 1.291 / 1.807 | 0.470 / 0.205 |
| spiked ellipse | 1.346 / 0.453 | 1.643 / 1.220 | 0.350 / 0.252 |
| touching discs | 0.180 / 0.158 | 1.180 / 1.215 | 0.088 / 0.043 |
| notched disc | 0.075 / 0.035 | 1.288 / 1.297 | 0.060 / 0.030 |

What the suite must show, and what the values say:
- **Discs and ellipses score alike and low** on all three columns, so none of them is secretly an elongation measure.
  `chr_radial_hf` has a partial elongation sensitivity (the 4:1 ellipse scores above the disc), because equal-arc
  sampling gives an elongated outline in-band energy even without pixels. Report it as such.
- **Spikes score above smooth shapes.** The 8- and 12-lobed stars and the spiked ellipse score clearly above every disc
  and ellipse on each column at each r. At a fixed relative amplitude a star's in-band outline energy barely depends on
  its lobe count, so do not expect `chr_radial_hf` to increase with it. The 5-lobed star's fundamental lies below the
  band, so it scores low on `chr_radial_hf`.
- **The lobulation test.** Two touching discs and the notched disc are smooth outlines, not spikes: they stand in for
  "two cells in close proximity" (Lin 2024) and for the crowded-nuclei false positives of Section 1. Every spike column
  must score both below the 8- and 12-lobed stars and the spiked ellipse, at each r. `chr_radial_hf` and
  `chr_concavities_per_um` pass. `chr_roughness` fails at r = 16 (notched disc 1.297 against spiked ellipse 1.220), so
  it is a lobulation measure, and its direction is 0. A column that fails gets direction 0 and is described as a
  lobulation measure in `directions.csv`.
- **Scale.** For the uniformly scaled shapes (disc, ellipses, stars, touching discs, notched disc), the two scales must
  agree within 10 % of the larger value. Do not widen the tolerance.
  - `chr_radial_hf` fails on every shape (the closest is the touching discs at 12 %): its low values are mostly
    pixel-staircase noise, which falls with r², and 8 or 12 lobes on an 8 px radius are 2–3 px wide, at the pixel
    pitch. Report those as failures caused by the resolution limit, with the radius below which they apply.
  - `chr_roughness` fails on the 12-lobed star at r = 8, for the same reason.
  - `chr_concavities_per_um` is exempt: it is a count per µm of outline, so with a fixed number of concavities it halves
    when r doubles, by construction.
  - The spiked ellipse is exempt, because its spikes keep a fixed pixel size and so shrink relative to the body as r
    grows, as real chromosome arms do on a larger figure. The consequence: fixed-size spikes score 3 × lower on
    `chr_radial_hf` on a body twice the size (1.346 → 0.453), so it is negatively size-confounded on its own cue.
    Section 6's `auc_area_tertiles` is the real-data check.
- **Elongation.** `chr_elongation` separates discs from ellipses, and no boundary column does.

### 5.3 Darkness and the surround profile

| column | definition | unit |
|---|---|---|
| `h_chr_mean` | mean H over `M_chr`, in float64; the earlier `core_od` | OD |
| `h_ring1_mean`, `h_ring2_mean`, `h_ring3_mean` | mean H over `ring_1`, `ring_2`, `ring_3` | OD |

Expect `h_chr_mean` and `od` to be nearly collinear (rank correlation 0.94 on the earlier pool); `od` stays as the
identity column and the stratifier of `auc_dark_stratified`. The chromatin-to-surround contrast
`h_chr_mean − h_ring3_mean` and the ring slope `(h_ring3_mean − h_ring1_mean) / 2` follow from these columns and are
not stored.

### 5.4 Chromatin texture

**The texture support:** `M_chr` eroded once by `scipy.ndimage.binary_erosion` with its default cross, minus every
pixel whose 3 × 3 neighbourhood contains a clipped pixel (Section 4). Both columns are NaN, reason `clipped`, when the
support has fewer than `NS_MIN_PX` = 10 pixels.

| column | definition | unit |
|---|---|---|
| `h_chr_lap` | std (ddof 0, in float64), over the texture support, of `cv2.Laplacian(H_w, cv2.CV_32F, ksize=1) / mpp**2`, divided by `h_chr_mean`. With ksize 1 the aperture reads only the 4-neighbours, which lie in `M_chr` for every support pixel, so no pixel outside the mask enters | 1/µm² |
| `h_chr_glcm_contrast` | `skimage.feature.graycomatrix` at distance 1, angles 0, π/4, π/2 and 3π/4, within `M_chr`'s bounding box (*whole* `M_chr`), on H / `h_chr_mean` quantised to `GLCM_LEVELS` = 64 equal bins over [0, `GLCM_REL_HI`] with `GLCM_REL_HI` = 2.0: level = `np.minimum(np.floor(H / h_chr_mean / GLCM_REL_HI * GLCM_LEVELS), GLCM_LEVELS − 1)` (rounding instead of flooring gives 30.19, not 23.86, on the checkerboard below). Pixels outside `M_chr` and clipped pixels get a 65th level, whose row and column are dropped from the matrix, so every pair touching them is excluded. `graycoprops(..., 'contrast')`, averaged over the four angles | ratio |

**Why relative bins.** Fixed OD bins made contrast scale with OD²: within-ROI ρ 0.60 with `h_chr_mean`, against
0.10–0.15 for relative bins. Within-object rank levels were rejected because they turn pixel noise on a flat body into
maximal texture (noise test below). Pixel noise that is fixed in OD spreads over more relative bins on a pale object.

Both texture columns have direction −1 (Section 5.0, texture).

**Noise test, mandatory.** On a 71 × 71 canvas at `mpp` 0.25, render a disc of radius 16 px at 0.20 OD on 0.03 OD
(`M_chr` = the disc) in four versions: flat; flat plus Gaussian noise of σ 0.003 everywhere; the noisy one blurred at
σ 1 px; and the granular disc of Section 5.8's reference construction (a 2 px checkerboard of +0.06 OD more than 3 px
inside the edge). Both texture columns must score the flat, noisy and noisy-blurred discs far below the checkerboard; a
column that does not is measuring noise or the blurred edge. Checked values, one implementation (flat / noisy /
noisy-blurred / checkerboard):
- `h_chr_glcm_contrast` 0.00 / 0.64 / 4.9 / 23.86 (the noisy-blurred value moves by about 0.05 with the noise draw).
  With 8 within-object rank levels it was 0.00 / 10.30 / 2.26 / 5.96, the noisy flat disc above the checkerboard, which
  is why rank levels were rejected;
- `h_chr_lap` 0.00 / 1.07 / 1.03 / 7.72.

**Clipped-pixel test, mandatory.** Render a uniform RGB disc of (40, 1, 100) on a (230, 200, 220) background, then set
G = 0 on a random 7 % of the disc's pixels (canvas 71 × 71, radius 16 px, `M_chr` = the disc). Both texture columns
must stay at their value on the unaltered disc, because the clipped pixels and their neighbourhoods are excluded.
Checked values, one implementation (`h_chr_lap` / `h_chr_glcm_contrast`): unaltered 0 / 0; with the clipped pixels
excluded 0 / 0; with them included 3.89 / 6.21. The exclusion costs support: 709 px fell to 383, because scattered
clipped pixels take their 3 × 3 neighbourhoods with them. On real data clipped pixels cluster in the densest core.

### 5.5 Surround: membrane, clear zone, white space, cytoplasm

| column | definition | unit |
|---|---|---|
| `ring_peak_prominence` | for k in {2, 3}: `h_ring_k_mean` minus the minimum of `h_ring_j_mean` over j < k; take the larger, divide by `h_chr_mean`, floor at 0. A rim that is darker than the band inside it scores above 0; a profile that only falls scores 0 | ratio |
| `ring1_frac_nuclear` | share of `ring_1` pixels above `t_nuc` | share |
| `ring1_frac_glass` | share of `ring_1` pixels at or below `t_nuc` and with `gray > 215` (the `WHITE_GRAY` of the floor-clip run): glass and lumen | share |
| `ring1_frac_pale` | share of `ring_1` pixels at or below `t_nuc` and with `gray ≤ 215`: pale tissue, cleared cytoplasm, halos. The three `ring1_frac_*` shares partition `ring_1` and sum to 1, so read them as two degrees of freedom | share |
| `white_frac_gray215_61` | share of the 61 px window with `gray > 215`, as `white_frac` in the floor-clip run computes it; the earlier `not_white` without the sign | share |
| `e_ring23_mean` | mean E (`hed[:, :, 1]`) over the cytoplasm ring: `ring_2` ∪ `ring_3` minus glass (`gray > 215`), minus pixels above `t_dark`, minus clipped pixels. NaN with reason `no cytoplasm` below `NS_MIN_PX` pixels | OD |
| `e_over_h_ring23` | mean E / mean H over the same cytoplasm ring | ratio |

All columns here have direction 0 (below, and C9).

**What the membrane columns can see (C2), checked 2026-10-02 with the verbatim refinement.** In an H&E interphase
nucleus the membrane is not a separate object that `M_chr` sits inside. Where it shows, it is a thin dark rim, and the
window's top multi-Otsu class takes the rim itself. Synthetic windows, `mpp` 0.25, blur σ 1 px:

| window | refinement mode | what `M_chr` becomes | `ring_peak_prominence` | `chr_hole_frac` |
|---|---|---|---|---|
| vesicular nucleus: 0.5 µm dark rim, pale interior | fallback | the rim | 0.06 | 0.78 |
| the same plus a nucleolus | fallback | rim + nucleolus | 0.00 | 0.80 |
| uniformly dark disc | base | the disc | 0.00 | 0.00 |
| 2 µm clump inside an intact rim | fallback | the clump | 0.15 | 0.00 |

- **`ring_peak_prominence`** scores highest on a dark envelope around an internal clump, but it is not specific to it.
  It gave 0.05 on an anaphase pair 4.5 µm apart with the detection on one plate (the partner plate falls in `ring_3`),
  and 0.03 on a metaphase plate whose clear zone is bounded by darker cytoplasm (Section 5.2's pair-test windows). No
  source describes prophase as a clump inside a dark rim: VCGP 2021 counts from prometaphase on, and its prometaphase
  is "A dark cluster or ring of chromosomes with variably protruding chromosome rods and spikes". Its direction is 0;
  the hypothesis that it marks an intact envelope (C2, C6) goes in the rationale only.
- **`chr_hole_frac`** (Section 5.2) is where a vesicular nucleus's membrane shows up, but a ring-shaped mitotic figure
  gives the same mask: VCGP 2021 counts "A dark cluster or ring of chromosomes with variably protruding chromosome rods
  and spikes" (prometaphase) and a "linear or ring-shaped aggregate" (metaphase) as mitotic figures. Its direction
  is 0.
- **On a uniformly dark nucleus no membrane column can work.** Pathologists separate those by their smooth outline
  (C5), which is the job of Sections 5.2 and 5.8 (b, c). A membrane column scoring them like mitoses is expected, not a
  defect.
- **Membrane test, mandatory.** Reproduce the four rows above within ±0.03 before real data. Construction, on the
  canvas of Section 5.2's pair test (201 × 201, 0.03 OD plus Gaussian noise of σ 0.003, any seed; shapes drawn over it
  with constant values that replace the noise; then blur σ 1 px and cast to float32), detection at the centre. With r
  the distance from a disc's centre, "the outer 2 px" of a disc of radius R means the pixels with R − 2 < r ≤ R; a
  distance-transform reading gives a thicker rim and moves the clump value to about 0.12.
  - vesicular: a disc of radius 24 px at 0.07 OD whose outer 2 px are 0.16 OD;
  - nucleolus: add a disc of radius 4 px at 0.16 OD, centred 3 px right and 2 px up of the centre;
  - dark disc: radius 18 px at 0.17 OD;
  - clump: a disc of radius 20 px at 0.07 OD whose outer 2 px are 0.15 OD, with a central disc of radius 8 px at
    0.20 OD.
  - Two independent implementations agree on all four rows under this reading (clump 0.148–0.150, vesicular
    0.058–0.061 / 0.78, across noise seeds; adding the noise after drawing moves them to about 0.145 and 0.066 / 0.79).

**Clear zone (C9).** No verified source says a mitosis sits in cleared cytoplasm. The H&E statement found is "The
cytoplasm of the mitotic cells was often larger during mitosis than in the resting cells" (Ikpatt 2002); a clear halo
is documented around apoptotic cells (rodent liver, NTP atlas), and this run's own false positives include smooth dark
nuclei with halos (Section 1). "Below `t_nuc`" lumps cleared cytoplasm, halos, lumen and glass together, which is why
`ring1_frac_glass` and `ring1_frac_pale` split them. The one contrast the ring shares are expected to carry: a dark
clump inside an interphase nucleus is surrounded by nucleoplasm above `t_nuc`, so its `ring1_frac_nuclear` should be
high.

### 5.6 Colour and pigment

RGB statistics are on the raw `uint8` image.

| column | definition | unit |
|---|---|---|
| `e_over_h_chr` | mean E / mean H over `M_chr_clip` (Section 4); NaN, reason `clipped`, below `NS_MIN_PX` pixels | ratio |
| `brown_frac_surround` | share of the 51 px window's pixels outside `M_sel` with `R > G > B`, `R − B > 40` and `gray < 180` | share |

The brown rule is a starting heuristic. If the top tiles of its montage (Section 8) are not pigment, calibrate it there,
log the change, and keep the original as `spec_brown_frac_surround`.

**Where pigment lands.**
- **Light brown pigment stays outside `M_chr`.** For lightly deposited melanin and the golden-brown granules
  (hemosiderin-like) seen on `410` and `291`, colour deconvolution assigns the OD mostly to D: for a brown OD direction
  ∝ (0.3, 0.5, 0.8) of unit norm, the H coefficient (`skimage.color.hed_from_rgb`) is about 0.05, against 1.0 for pure
  hematoxylin. These are deconvolution coefficients, not orthogonal projections. `brown_frac_surround` is the tool
  for it.
- **Dense melanin can enter `M_chr`.** Its H value grows linearly with density at that 0.05 coefficient, and in a dense
  deposit the blue channel saturates first, which pushes the measured OD vector towards neutral (coefficient about
  0.70). `chr_resid_frac` (5.8 a) scores it, as long as its rim keeps enough non-saturated pixels.
- **Neutral black pigment *is* the `M_chr` of its candidates.** A neutral (equal-OD) absorber of unit OD norm gets an H
  coefficient of about 0.70. On `211` the top-10 false positives are black granular clumps; their darkest pixels have
  D 0.152 against 0.034 for the mitoses (`rgb2hed` D averaged over the darkest 10 % by H of a 21 × 21 px core, as
  `colour_feats` in `low_precision_rois/roi_mechanism_measure.py` computes it). `chr_resid_frac` is the feature for
  this case.

**Pigment check, mandatory.** Report `chr_resid_frac`'s within-click AUC (Section 6) on each of `211`, `289` and `339`
(black or brown-black clusters seen near top-10 false positives). The montages of `chr_resid_frac` and
`brown_frac_surround` (Section 8) are the visual check.

### 5.7 Context and density

| column | definition | unit |
|---|---|---|
| `dark_area_frac_101` | share of the 101 px window's pixels above `t_dark` | share |
| `nuc_area_frac_101` | share of the 101 px window's pixels above `t_nuc` | share |

Both have direction 0 (C10). Store `t_nuc` and `t_dark` per ROI in `clicks.csv`.

### 5.8 In-mask pigment, thin protrusions, edge sharpness, size relative to neighbours

These four families target the dark false positives of the failing ROIs (Section 1). Every rule at the top of Section 5
applies. Every constant named here goes into the spec JSON.

**(a) `chr_resid_frac`, `chr_ns_share`: colour off the hematoxylin–eosin plane, on the chromatin and its 2 px rim.**
- **Per-pixel optical density:** `od = -np.log10(np.maximum(rgb_w.astype(np.float64), 1.0) / 255.0)`, with channels in
  R, G, B order.
- **The region:** `M_res` = `scipy.ndimage.binary_dilation(M_chr, structure=skimage.morphology.disk(RESID_DILATE_PX))`,
  with `RESID_DILATE_PX` = 2, on *whole* `M_chr`, within the 51 px window.
- **Non-saturated pixels:** `M_res_ns` = `M_res` ∧ all three channels in [`NS_LO`, `NS_HI`] = [20, 235].
  - Store `chr_ns_share` = |`M_res_ns`| / |`M_res`| (role `diagnostic`).
  - If |`M_res_ns`| < `NS_MIN_PX` = 10, then `chr_resid_frac` is NaN with reason `outside_20_235`.
- **Why the rim is included (measured 2026-10-02).** On `M_chr` alone, the black clumps on `211` are mostly saturated:
  44 % of `211`'s A3 top-10 false positives had fewer than 10 non-saturated pixels, so the feature was NaN on its own
  target, and 35 % of `289`'s mitoses were NaN too. The clumps' partly transmitting 2 px edge carries their colour.
  With the rim included, nothing was NaN on the seven ROIs checked.
- **The normal to the plane:** `h`, `e` = `skimage.color.rgb_from_hed[0]` and `[1]`, each divided by its norm. Then
  `n = np.cross(h, e)`, normalised, with its sign chosen so that `n[2] > 0`. That gives about (−0.332, −0.081, 0.940).
- **The feature:** `v` = the mean of `od` over `M_res_ns`, and `chr_resid_frac = (v @ n) / np.linalg.norm(v)`. It is
  signed and dimensionless.
- **What it measures.** Any mixture of the two stains, as `rgb2hed` defines them, lies in the plane spanned by `h` and
  `e` and scores 0. A neutral (black or grey) absorber scores 0.304, and a brown OD vector ∝ (0.3, 0.5, 0.8) scores
  0.618.
- **Unit test, mandatory.** Use uniform 21 × 21 patches rendered as `I = round(255 · 10**(−OD))` per channel, with
  the OD scaled so that every channel lies in [20, 235], and `M_chr` = the central 11 × 11 pixels, so that `M_res` stays
  inside the uniform patch. Expected values: pure `h` → 0.000 ± 0.02; 0.7·`h` + 0.3·`e` → 0.000 ± 0.02; neutral →
  0.304 ± 0.02; brown → 0.618 ± 0.02; a patch with every channel below 20 → NaN.
- **Fixed `rgb2hed` vectors.** Each scanner's real stain vectors differ from them (per-ROI Macenko H vectors are 5–15°
  away at the scanner median, `floor_robustness/macenko.csv`). That shifts every candidate within an ROI alike, so the
  within-click comparisons of Section 6 are unaffected, but pooled values are **not** comparable across scanners; say
  so wherever you report them. The mitoses' median baseline was about 0.18–0.21 on two Aperio ROIs (`300`, `339`),
  about 0.08 on `211`, and about 0.01 on `527` (Hamamatsu XR).
- **Real-data anchor (2026-10-02).** On A3 top-10 false positives against GT mitoses refined at the annotation centres
  (a sample, not this spec's population), the AUC that a false positive scores higher than a mitosis was 0.83 on `211`
  and 0.86 on `289`, both ROIs with pigment, and 0.20–0.49 on ROIs without pigment (`300`, `339`, `425`, `460`, `527`),
  where the *mitoses* score slightly higher. So expect a pigment-specific signal at the top of the range, not a monotone
  one. Do not tune `RESID_DILATE_PX` or `NS_LO` to change this.
- **Direction hypothesis:** −1 (pigment is higher).

**(b) `chr_tophat_edge`, `chr_tophat_core`: the share of the object's OD carried by structures thinner than about
1 µm, at the surface and inside.**
- **The thin-structure map.**
  - `r_th = max(1, round(TOPHAT_R_UM / mpp))` with `TOPHAT_R_UM` = 0.5, which gives 2 px at 0.23–0.25 µm/px.
  - `se = skimage.morphology.disk(r_th).astype(np.uint8)`.
  - `T = H_w − cv2.morphologyEx(H_w.astype(np.float32), cv2.MORPH_OPEN, se, borderType=cv2.BORDER_REPLICATE)`: a
    grey-level white top-hat, the high-OD structures narrower than the disc.
  - The disc is 5 px, about 1.25 µm, across. A single chromatid is "0.68 ± 0.04 μm" wide and a cohesed chromosome
    "1.23 ± 0.1 μm" (Cisneros-Soberanis 2024, by EM), so the threshold sits at the physical width of one chromosome
    arm. Whether a scan resolves an arm that thin is a separate question (5.8 c, and Section 6's per-scanner columns).
- **Regions, following the mask** (`M_chr` *whole*):
  - `R` = the window pixels with `distance_transform_edt(~M_chr) * mpp <= TOPHAT_RING_UM`, with `TOPHAT_RING_UM` = 1.0,
    so that fuzzy protrusions up to 1 µm outside the mask count. `M_chr` itself is included.
  - `M_core` = the pixels with `distance_transform_edt(M_chr) * mpp > TOPHAT_CORE_UM`, with `TOPHAT_CORE_UM` = 0.5:
    chromatin more than 0.5 µm inside its own edge.
  - `R_edge` = `R` ∧ ¬`M_core`.
- **The features:**
  - `chr_tophat_edge = T[R_edge].sum() / H_w[M_chr].sum()`: thin structure at the surface, the C1 cue;
  - `chr_tophat_core = T[M_core].sum() / H_w[M_chr].sum()`: thin structure inside, which granular interphase chromatin
    also has;
  - both NaN, reason `small piece`, when `M_chr` has fewer than 30 px, as in Section 5.2.
- **Why the region follows the mask.** A disc about the centroid follows the mask's size, not its shape: on an about
  4:1 bar with spikes at its tips, the tips lie outside it. And one undivided region mixes opposite cues: a disc with
  2 px internal granules scored above a disc with 2 px surface spikes. The split puts granules in `chr_tophat_core` and
  spikes in `chr_tophat_edge`.
- **Synthetic suite, mandatory**, expected ordering written first, checked values below. Render each shape as an H
  image, 0.20 OD on a 0.03 background, at radius 8 and 16 px:
  - a disc scores low on both columns;
  - the 12-spike shape (`w` = 2) scores clearly above the disc on `chr_tophat_edge` at both scales;
  - the same shape with `w` = 6 scores far below the `w` = 2 version on `chr_tophat_edge`;
  - the `w` = 2 shape blurred with Gaussian σ = 1 px: report the drop. Blur erasing the cue is the resolvability
    question Section 6 asks per scanner;
  - the granular disc scores clearly above the disc on `chr_tophat_core`, and no more than 0.01 above it on
    `chr_tophat_edge`;
  - the tip-spiked bar scores clearly above the plain bar on `chr_tophat_edge` at both scales;
  - **a known confound, not a pass:** the vesicular window of Section 5.5, refined verbatim, scores `chr_tophat_edge`
    0.34, above the 2 px spikes (0.256), because its `M_chr` is a thin continuous rim. A top-hat cannot tell a smooth
    thin membrane from discrete spikes. Record the value; `auc_hole_gt03` and `auc_hole_le03` (Section 6) read the effect on real
    data.
- **The scale rule does not apply.** These columns are exempt from Section 5.2's scale rule, because their width
  threshold is physical by design. At r = 8 the plain bar is 5 px wide, exactly the structuring element's width, so the
  opening removes only its tip corners (12 of its 85 px); a bar thinner than 5 px would be removed entirely.
- **Small bodies score high too.** The whole of a body narrower than about 1.25 µm counts as thin structure. So
  karyorrhectic fragments and other tiny round bodies, which are imposters (C5, C7), raise `chr_tophat_edge` as spikes
  do. Section 6's `auc_area_tertiles` reads this; write the confound into `directions.csv`.
- **Direction hypotheses:** `chr_tophat_edge` +1, with the small-body and vesicular-rim confounds above;
  `chr_tophat_core` −1 (Section 5.0, texture).

**(c) `chr_edge_sharpness_per_um`: how crisp the chromatin boundary is.**
- **The boundary `B`:** the hole-filled largest 4-connected piece of `M_chr` (Section 5.2), minus its `binary_erosion`
  with a 3 × 3 cross. NaN, reason `small piece`, below 30 px.
- **The gradient:** `gx = cv2.Sobel(H_w, cv2.CV_32F, 1, 0, ksize=3) / 8` and `gy` likewise; `g = np.hypot(gx, gy) / mpp`,
  in OD per µm.
- **The step:** `step = h_chr_mean − h_ring1_mean`, from Section 5.3. If `step` < `EDGE_MIN_STEP` = 0.005 OD, the
  feature is NaN with reason `no edge step`.
- **The feature:** `chr_edge_sharpness_per_um = g[B].mean() / step`, in 1/µm. It is roughly the inverse width of the
  edge.
- **Synthetic suite, mandatory.** A disc and a 2:1 ellipse, at radius 8 and 16 px, 0.20 OD on 0.03, each
  Gaussian-blurred at σ = 0, 1 and 2 px:
  - the value must fall strictly with σ;
  - the ellipse must score within 15 % of the disc at the same σ;
  - the two radii must score within 15 % of each other.
- **What it is for.** Blurred dark objects are a named false-positive class. The value is also the covariate for
  whether boundary features are resolvable on a given scan: in MIDOG 2022, 27.5 % of the mitoses found only with PHH3
  help were classed "out of focus" (Aubreville et al. 2024, Table 7).
- **Direction hypothesis:** 0. Blurred objects should score low, smooth pyknotic bodies high, and ragged mitotic edges
  in between, so no monotone direction is claimed.

**(d) `nbr_n`, `nbr_nuc_area_med_um2`, `chr_area_over_nbr`: chromatin size relative to the neighbouring nuclei.**
- **The window:** the 101 px window of Section 4.
- **The components:** 4-connected components of `hem_101 > t_nuc`.
- **Drop:**
  - every component touching the (truncated) window's edge;
  - every component intersecting the candidate's own `M_chr`, placed as in Section 4;
  - every component with an area outside [`NBR_MIN_UM2`, `NBR_MAX_UM2`] = [8, 150] µm². This drops fragments, and
    merged clusters or tissue sheets.
- **The outputs:**
  - `nbr_n` (role `diagnostic`) is the number of components kept;
  - `nbr_nuc_area_med_um2` is their median area, NaN with reason `few neighbours` when `nbr_n` < 3;
  - `chr_area_over_nbr = chr_area_um2 / nbr_nuc_area_med_um2`.
- **The risk in `t_nuc`.** `t_nuc` admits cytoplasm (Section 4). In dense tissue, neighbouring nuclei and their
  cytoplasm may merge into blobs above 150 µm² and be dropped, leaving `nbr_n` < 3 on exactly the crowded ROIs this
  column is for; and the ratio then compares `M_chr`'s chromatin with a neighbour's nucleus plus cytoplasm. The unit
  test below cannot see this.
- **Real-data check, mandatory.** For 10 random accepted candidates on each of `245` (lymphosarcoma) and `462`
  (sarcoma), show the 101 px window with the kept components outlined in one colour and the dropped ones in another,
  embedded in the notebook. If the kept components are not nuclei, the definition is defective (Section 2): fix it, for
  example with `t_dark` or different area bounds, and keep this version as `spec_` columns.
- **Unit test, mandatory.** A 101 × 101 field at `mpp` 0.25 holding 12 non-touching discs of radius 12 px (about
  28 µm²), one disc at the centre standing in for the candidate's `M_chr`, and two discs cut by the window edge. It must
  give `nbr_n` = 12 and a median equal to the rendered disc's pixel count × mpp², within 2 %.
- **Direction hypothesis:** 0. It should be +1 against apoptotic bodies and lymphocytes, which are smaller than their
  neighbours, but −1 against hyperchromatic interphase nuclei.

**Reference shapes and checked values (computed 2026-10-02 from the definitions above).** Reproduce these before
running on real data.

**Construction:**
- `mpp` = 0.25. The canvas is 51 px for r = 8 and 71 px for r = 16. Shapes are centred, at 0.20 OD on a 0.03
  background, as float32. The whole canvas is the window `H_w`.
- The disc is `x² + y² ≤ r²`. The ellipse has semi-axes r·√2 and r/√2.
- **The 12-spike shape** is the disc plus 12 rectangles at angles 2πk/12. Each rectangle covers `r − 1 ≤ t ≤ r + 4`
  along the ray and `|o| ≤ w/2` across it, with spike width `w` = 2. The "6 px lobes" variant uses `w` = 6.
- **The granular disc** is the disc plus 0.06 OD on the pixels more than 3 px inside its edge
  (`distance_transform_edt(disc) > 3`) where `(floor(x / 2) + floor(y / 2)) % 2 == 0`, with x and y the centred
  coordinates: a 2 px checkerboard. `M_chr` is the disc.
- **The bar** is `|x| ≤ r/4` and `|y| ≤ r`, about 4:1. **The tip-spiked bar** adds six rectangles built like the
  12-spike shape's (`w` = 2, `r − 1 ≤ t ≤ r + 4`), at angles π/2 + {−0.5, 0, +0.5} rad and −π/2 + {−0.5, 0, +0.5} rad,
  where `t = x·cos(a) + y·sin(a)`.
- Blur is `scipy.ndimage.gaussian_filter` with default settings. Ring 1 uses `scipy.ndimage.distance_transform_edt`.
- In every synthetic test, `M_chr` is the shape's own *unblurred* pixel set, even when the image is blurred.

**Checked values:**
- **`chr_resid_frac` unit test:** 0.001, 0.000, 0.304 and 0.615 for pure `h`, the mix, neutral and brown. The rendered
  pixels are (77, 70, 149), (106, 60, 165), (64, 64, 64) and (137, 90, 49); the ODs are 0.8·`h`,
  0.8·(0.7`h` + 0.3`e`), 0.6 per channel, and 0.9·(0.3, 0.5, 0.8).
- **The top-hat columns,** r = 8 / r = 16; two independent implementations reproduce every value:

  | shape | `chr_tophat_edge` | `chr_tophat_core` |
  |---|---|---|
  | disc | 0.035 / 0.000 | 0.000 / 0.000 |
  | 2 px spikes | 0.256 / 0.086 | 0.000 / 0.000 |
  | 6 px lobes | 0.014 / 0.030 | 0.000 / 0.000 |
  | 2 px spikes, σ = 1 | 0.118 / 0.043 | 0.001 / 0.000 |
  | granular disc | 0.033 / 0.000 | 0.056 / 0.092 |
  | bar | 0.120 / 0.034 | 0.000 / 0.000 |
  | tip-spiked bar | 0.313 / 0.146 | 0.000 / 0.000 |

- **`chr_edge_sharpness_per_um` for the disc at σ = 0 / 1 / 2:** 1.93 / 1.53 / 1.16 at r = 8, and 1.97 / 1.47 / 1.04 at
  r = 16. The ellipse is within 6 % of the disc.

## 6. The validity report (descriptive; no selection, no verdict)

**Rows.** One row per column with role `feature` or `spec` in `directions.csv`.

**The AUC.** Every AUC in this section is computed the same way. For a set of groups (clicks, or strata within clicks):
1. In each group, drop the rows where the column is NaN, then keep the group only if it still has at least one mitosis
   (`is_mitosis`) and one other row.
2. In each kept group, U = the number of (mitosis, other) pairs in which the mitosis has the larger value, plus 0.5 for
   each tie.
3. AUC = ΣU / Σ(n1·n0), summed over the kept groups, where n1 and n0 are each group's mitosis and other counts.

The AUC is the raw value of P(mitosis > other): never folded about 0.5, so a value below 0.5 means the column is lower
on mitoses.

**Columns of `validity.csv`**, computed on accepted rows with `in_recorded_pool`, the pool production ranks:
- `n_clicks`, `n1`, `n0`: the clicks kept and the rows they contribute.
- `auc_within_click`: the AUC with clicks as the groups.
- `auc_dark_stratified`: **the AUC beyond darkness.** Within each click, split its accepted rows into five strata by
  within-click rank of `od`: `pd.qcut(od.rank(method='first'), 5, labels=False)`. The AUC with the click-by-stratum
  cells as the groups. Also store `n_strata_used` and `n1_dark_stratified`. A value near 0.5 means the column adds
  nothing once the window's darkness is held roughly fixed, whatever its `auc_within_click`. Darkness columns
  themselves will sit nearer 0.5 here by construction; that is the point of the column, not a defect.
- `auc_area_tertiles`, for the boundary columns only (NaN elsewhere): as `auc_dark_stratified`, with three strata of
  `chr_area_um2` in place of five of `od`. Spikes on a small object approach the pixel pitch (Section 5.2, scale), and
  small bodies raise `chr_tophat_edge` (5.8 b), so a boundary column that separates only between size strata is
  measuring size.
- `auc_hole_gt03`, `auc_hole_le03`, for `chr_tophat_edge` only: `auc_within_click` on the rows with `chr_hole_frac`
  above 0.3 and at most 0.3. A ring-shaped `M_chr` (a vesicular nucleus's rim) raises `chr_tophat_edge` (5.8 b).
- `auc_wc_3dhistech`, `auc_wc_aperio_cs2`, `auc_wc_hamamatsu_xr`, `auc_wc_hamamatsu_s360`: `auc_within_click`
  restricted to each scanner's clicks, each with its click count. Scanner is confounded with tumour type for 42 of the
  49 ROIs, and S360 is 3 ROIs, so read these as scanner-and-tissue groups. For the boundary columns this is the "is the
  cue resolvable on this scan" check; report it beside each group's median `chr_edge_sharpness_per_um`.
- `auc_wc_failing9`, `auc_wc_other40`: `auc_within_click` by `roi_group`, each with its mitosis count. This shows
  whether a column carries signal where darkness fails. The grouping is post hoc, so it is descriptive.
- `nan_rate_mitosis`, `nan_rate_other`: the share of accepted mitosis and other rows where the column is NaN.
- `synthetic_pass`: pass, fail or n.a., from Section 5's suites.
- `direction_hypothesis`, copied from `directions.csv`, and `direction_measured`, the sign of
  `auc_within_click − 0.5`. Report agreement; do not change the hypothesis after the fact.

**`validity_depth300.csv`:** the same columns on all accepted rows (`cand_rank` < 300). The primary table stays on the
pool production ranks; the depth-300 table is where the failing ROIs have enough mitoses to read. Print the mitosis
count per `roi_group` for each table.

**`redundancy.csv`:** the matrix of mean within-click Spearman correlations (rows where both columns are non-NaN, clicks
with at least 10 such rows) between the role-`feature` columns and `od`. Print the pairs above 0.8 in absolute value.

**What you may conclude from this report:** whether each measurement is valid, and what it actually measures. **What
you may not conclude:** which features to use, how to weight them, or whether any ranker is better. Write the report
in that voice.

## 7. Gates (halt and report) and flags (report only)

The gates protect what this run cannot check any other way: that the copied refinement is the one at `c767898`, that
the new search glue reproduces the recorded A3 pool, and that the new candidate glue (channel, coordinates, mask)
reproduces the earlier run's measurements. Copied code is not re-checked against itself.

**Gates:**
- **A. The refinement.** The SHA-1 of the extracted `refine_template_box` segment + `"\n"` (Section 3) is
  `21c43a20db00dccf7ac3bab88aee423c345441b6`. Halt otherwise.
- **B. The search reproduces the recorded A3 pool.** For every click, the first `n_detections_recorded` deep rows
  against the recorded A3 pool (the `arm == 'A3'` rows of `pool.csv.gz` for its `(file_name, seed_index)`, ordered by
  their `rank` column):
  - exact: `n_detections_recorded` = the recorded `n_detections`, and `cx`, `cy`, rank order, `bucket` and
    `matched_ann_id` row by row;
  - within 1e-4 absolute: `score`, `od`, and `threshold` against the recorded `threshold`.

  A difference that a float shift below 1e-4 can explain is a near-tie: two scores within 1e-4 of each other swapping
  order (with their labels swapping with them), or a peak within 1e-4 of the threshold or of the 100th peak's score
  entering or leaving the pool. Report near-ties and continue. Any other difference halts.
- **C. The candidate measurements reproduce the earlier run.** On every row whose `(file_name, cx, cy)` exactly equals
  a location in `results/candrefine_rerank_mo_pad61_select_pool.csv` (about 15,700 distinct locations, 57 % of the
  recorded A3 pool's; print the count): `refusal` equals that run's `refusal`, and on accepted rows
  `|h_chr_mean − core_od|` ≤ 1e-4. A wiring error (the clipped or blanked channel, swapped coordinates, the dilated
  mask) mismatches nearly every row, while a float near-tie in the refinement's multi-Otsu flips a few. Print the
  mismatches; halt if more than 1 % of shared rows mismatch.

**Flags (report):**
- refusal counts by class and bucket;
- for every column, the NaN share by reason and by bucket, on accepted rows;
- the share of accepted rows with `chr_touches_window`, by bucket;
- the share of accepted rows with `chr_resid_frac` NaN (`outside_20_235`) and the median `chr_ns_share`, by bucket;
- the share with `nbr_n` < 3, by domain and by `roi_group`;
- the number of distinct locations, and per click `n_peaks`, `n_detections_deep` and the search seconds (also in
  `clicks.csv`).

## 8. Outputs and files

Notebook: `production_hematoxylin_only/multiotsu_pad_arms/candidate_features_a3_top300.ipynb`.

Tables, in `results/`, prefix `candidate_features_a3_top300_`:
- `candidates.csv`: 88,200 rows, one per (click, candidate): the Section 5.1 columns, every column with role
  `feature`, `diagnostic` or `spec`, and one reason column per family:

  | column | NaN reasons it can hold | columns it covers |
  |---|---|---|
  | `nan_reason_geom` | `small piece` | the 5.2 table except `chr_area_um2`, `chr_n_pieces` and `chr_concavities_per_um`; `chr_tophat_edge`, `chr_tophat_core`, `chr_edge_sharpness_per_um` |
  | `nan_reason_concav` | `small piece`, `convexity error` | `chr_concavities_per_um` |
  | `nan_reason_pair` | `no anchor`, `no partner` (all four), `round piece` (the two angles) | the `pair_*` columns |
  | `nan_reason_texture` | `clipped` | the texture columns |
  | `nan_reason_colour` | `clipped` | `e_over_h_chr` |
  | `nan_reason_cyto` | `no cytoplasm` | `e_ring23_mean`, `e_over_h_ring23` |
  | `nan_reason_resid` | `outside_20_235` | `chr_resid_frac` |
  | `nan_reason_edge` | `no edge step` | `chr_edge_sharpness_per_um` |
  | `nan_reason_nbr` | `few neighbours` | `nbr_nuc_area_med_um2`, `chr_area_over_nbr` |

  Refused rows have every feature NaN and every reason column empty; their `refusal` says why. A NaN this table does
  not explain is a bug or a defect: log it.
- `clicks.csv`: 294 rows. Per click: the click and its box (`bs`, `tcx`, `tcy`); the floor `b`; the search's
  `threshold`, `n_peaks`, `n_detections_deep`, `n_detections_recorded` and search seconds; `t_nuc` and `t_dark`; `mpp`
  and `radius_px`.
- `directions.csv`: one row per column with a role: `column`, `role`, `claims_to_measure` (one sentence a pathologist
  would accept), `direction_hypothesis` (+1 larger is more mitosis-like, −1, or 0 for no hypothesis), `rationale`,
  `confound`, `test`. Use the directions and rationale that Section 5.0 or the defining section gives. It is written by a
  cell that runs after the synthetic suites and before the smoke run, from a literal in the notebook (Section 2).
- `validity.csv`, `validity_depth300.csv`, `redundancy.csv` (Section 6), `synthetic_tests.csv` (Section 5).
- `adversarial_log.csv` (Section 2).

Also:
- **The spec JSON:** `production_hematoxylin_only/multiotsu_pad_arms/data/candidate_features_spec.json`, holding every
  constant (including Section 5.8's), the refinement SHA-1 and the `directions.csv` SHA-1, with sorted keys and no
  timestamps, so a rerun reproduces its SHA-1. The `data/` subfolder is where this folder's JSON and CSV outputs live,
  since `fee67b8`.
- **Montages, embedded in the notebook only.** Do not write image files into the repo; figures live in notebook outputs
  (`plt.show()`), by the user's 2026-10-02 rule. For your own look (Section 2, step 4), save temporary copies under
  `$TMPDIR` and open those with the Read tool.
  - **One figure per montage feature** (Section 5), 15 figures.
  - **Rows:** one per ROI: `300`, `201`, `013`, `529`, `425`, `211`, `546`, `524`. The pigment features use pigment
    ROIs instead: `chr_resid_frac` uses `211`, `289`, `339`; `brown_frac_surround` uses `509`, `410`, `291`.
  - **Tiles per row:** the ten highest values on the left and the ten lowest on the right, from the ROI's accepted
    depth-300 rows pooled over its 6 clicks, deduplicated by location, NaN excluded, picked greedily in value order and
    skipping any location within 7.5 µm of one already picked.
  - **Each tile:** the 61 px RGB window centred on the location, with `M_chr` outlined, a border coloured by bucket
    (mitosis green, look-alike orange, unannotated grey), and the value as a caption. Size the figure so that each tile
    shows at least 2 screen pixels per image pixel when you open it.
  - **Plus** the neighbour-component check of Section 5.8 (d).
  - Your one-sentence judgement per figure goes in the notebook's final markdown cell.

The final markdown cell holds, in this order: the per-feature validity verdicts in your own words; the montage
judgements; the adversarial log; what in this spec you changed and why; what you could not check.

Smoke runs write `_smoke`-suffixed copies of everything. The full run refuses to start if any output exists. Every
output is written at the end.

## 9. Environment, conventions, order of work

- Interpreter: `/Users/mohinianand/anaconda3/bin/python3`, the only one that imports `cv2` and `midog_utils`.
- Execute with `/Users/mohinianand/anaconda3/bin/jupyter nbconvert --to notebook --execute --inplace
  --ExecutePreprocessor.timeout=14400 <path>`. The notebook sits two levels below the repo root:
  `sys.path.insert(0, '../..')`, images at `../../images/...`, results at `../../results/`.
- Expected runtime: on the Intel Mac the deep search took about 45 minutes and refinement plus features were projected
  at 10 to 20 ms per location (about 84k locations); this Mac is faster. Use the smoke run's projection.
- Docstrings: the opening `"""` on its own line, a brief summary, one line per parameter with type and a short
  description, one line for the return. Function definitions and calls on a single line, never wrapped.
- Do not commit. Do not modify or delete any existing file. Kill only a PID you started.
- Order of work: the unit tests and synthetic suites (Section 2, step 1), then `directions.csv` (step 2), then the
  smoke run, then the full run, then the validity report and montages (steps 3 and 4), then the final markdown cell.
- Report back with, in this order:
  - both validity tables, with `auc_dark_stratified` and the `roi_group` columns shown for every feature;
  - the pairs above 0.8 in the redundancy matrix;
  - the synthetic failures;
  - the montage judgements;
  - the adversarial log.

  Lead with what is wrong, not with what reproduced.

## 10. Not in this spec

- Any ranking, combination, weighting, selection, statistics beyond the descriptive validity report, or verdict.
- Fresh clicks, held-out ROIs, or any search other than the one deep A3 search of Section 3.1.
- Changes to the production refinement, template, floor clip or matcher, or to the benchmark harness's border band.
  The only departures from production are lifting the peak cap to read deeper than its 100-peak pool, and the band,
  which is the harness's correction for ROI crops (Section 3.1).
- A nucleolus detector or a learned feature. If a cue in Section 5 turns out to need one, say so in the log and leave
  it out.
- **Stain-normalised channels.** No per-slide (Macenko) stain vectors anywhere in the features; Section 5.8 (a) uses
  the fixed `rgb2hed` vectors on purpose.
- The floor stress test and the multi-round feedback work (`FUTURE_IDEAS.md` §1 and §4.1).
