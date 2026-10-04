# Per-candidate feature extraction on the adopted A3 search's top 300 (`mo_pad61` + floor clip + border band), with built-in adversarial checking: agent prompt

**Status:** spec only, written 2026-10-01 and amended 2026-10-02. Nothing here has been run. It replaces the feature part of
`CANDIDATE_REFINEMENT_RERANK_PLAN.md` (Section 5). It does **not** replace that plan's ranking, selection or verdict
parts; those are deferred to a separate, later analysis that will read this notebook's outputs.

**The 2026-10-02 amendment.** It follows `low_precision_rois/2026-10-02-low-precision-rois-mechanism-and-literature.md`.

| what changed | where |
|---|---|
| Candidates now come from one new deep search under the adopted A3 settings (floor clip + 30 px border band), top 300 per click | Section 3.1 |
| Border candidates and the padded re-run removed; 101 px windows truncated, not padded | Sections 3.1, 4 |
| New identity columns | Section 5.1 |
| The pigment rationale corrected | Section 5.6 |
| Four new feature families | Section 5.8 |
| Darkness-stratified, per-scanner and per-ROI-group validity columns, and a depth-300 validity table | Section 6 |
| Montages embedded in the notebook, not written as PNGs; three failing ROIs added | Section 8 |
| How the refinement SHA-1 is computed | Section 3 |
| `211` added to the smoke ROIs | Section 2, step 7 |
| Gates 1–3, 6 and 7 re-pointed at the recorded A3 run; new flags | Section 7 |
| Stain normalisation ruled out | Section 10 |

Everything else is unchanged.

**The second 2026-10-02 amendment (review against the mitosis literature and the code).** It checked every feature
against published counting criteria and ran the verbatim refinement and the Section 5.8 reference shapes on synthetic
windows. Sources and the criterion map are in Section 5.0.

| what changed | where |
|---|---|
| Criterion → column map with verified sources, as the rationale for `directions.csv` | Section 5.0 |
| Primary sources added for the cues | Section 1 |
| Why the clip exists; clip (production) vs band (harness); the D3 carve-out | Sections 1, 2, 3.1, 10 |
| `t_nuc` ≡ `b_otsu` and `t_dark` ≡ `b_otsu_upper`, asserted | Sections 3, 4, 5.7 |
| Chromatin colour and texture exclude only truly clipped pixels (a channel at 0 or 255), not every pixel below 20; the F9 rationale corrected (the darkest chromatin pixels are dense blue-purple, not neutral); `chr_clip_share` added | Sections 4, 5.4, 5.6, 6, 7 |
| Pair columns read the 101 px window at `t_dark`, not `M_chr`, with an oblique test; `chr_two_piece_sep_um` kept as `spec_` | Section 5.2 |
| `chr_hole_frac` added; membrane columns re-described; membrane test cases | Sections 5.2, 5.5 |
| Elongation direction note; suite shapes pinned down; "increasing with lobe count" corrected; lobulation test; concavity depth in µm | Section 5.2 |
| `condensed_disc` removed (a size transform); `condensed_nuc` gets a surround test | Section 5.4 |
| GLCM on 64 bins of H / `h_chr_mean` over [0, 2]; gradient and Laplacian in µm units; a noise test for every texture column | Section 5.4 |
| Clear-zone, `ring_peak_prominence`, `rim_core_nuc`, `chr_hole_frac` and `pair_*` directions set to 0; ring shares partitioned into nuclear, glass and pale | Sections 5.2, 5.5 |
| Dense-melanin wording; H-coefficient wording | Section 5.6 |
| `chr_tophat_frac` split into `chr_tophat_edge` and `chr_tophat_core` on a mask-following region; old one kept as `spec_` | Section 5.8 (b) |
| The 101 px window's `M_chr` offset corrected for truncated windows | Sections 5.2, 5.8 (d) |
| Outline columns read the hole-filled piece; the vesicular rim recorded as a top-hat confound; `auc_split_hole` | Sections 5.2, 5.8 (b), 6 |
| Cytoplasm colour columns `e_ring23_mean`, `e_over_h_ring23` (VCGP QRG 2022), direction 0 | Sections 5.0, 5.5 |
| Boundary-column AUC within size tertiles | Section 6 |

An independent adversarial review then checked these edits against the code, real data (the recorded A3 pool, seeds
0–1) and the literature; its corrections are folded in.

## 0. What this is, and what it is not

**The job.** For every one of the 294 recorded `mo_pad61` clicks on the 49 ROIs:
1. take the 300 highest-scoring post-NMS candidates of a deep search under the adopted A3 settings (Section 3.1);
2. measure a defined set of per-candidate features on each;
3. validate that each feature measures what it claims to measure;
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

Two further problems were structural, not definitional:
- `od` and `core_od` had a rank correlation of 0.94, so one measurement held two votes;
- the direction of a feature was read from a pooled AUC, which mixes between-ROI differences into a within-pool
  statistic. F9 read 0.48 pooled and 0.58 within a click.

**The general failure:** the implementing agent treated the spec as truth and its own job as transcription. This spec
asks you to treat the spec as a hypothesis and your job as measurement. The protocol in Section 2 is how.

**What the 2026-10-02 analysis adds.** Read its write-up, `low_precision_rois/2026-10-02-low-precision-rois-mechanism-and-literature.md`, §2–4, before
starting.

**Both ranking axes are darkness detectors.**
- Per-ROI P@10 is predicted by one count: the unannotated nuclei in the whole ROI that are darker than the ROI's
  5th-darkest mitosis. Spearman ρ −0.75 at `tm_score`.
- Nine ROIs fail: `211`, `242`, `524`, `546`, `462`, `527`, `410`, `245`, `289`.
- Their top-10 false positives are dark objects:
  - black pigment (`211`);
  - small, smooth, dark round nuclei beside lumen (`242`);
  - smooth dark nuclei with halos (`546`, `462`);
  - crowded, hyperchromatic nuclei (`524`, `527`, `245`).

**What that means for this run.**
- The question to ask of each feature is whether it carries information **beyond darkness**. Section 6's
  `auc_dark_stratified` measures exactly that. Report it; do not select on it (Section 0).
- It must also be measured where those failing ROIs' mitoses are, which is why the depth is 300.
- It must be measured on the false positives production will actually rank. That is why the candidates come from the
  A3 search (the adopted floor clip, on the band-corrected benchmark harness; Section 3.1), not from the older
  unclipped list.
- The literature names the cues that separate mitoses from dark imposters (full map and sources in Section 5.0):
  - no nuclear membrane, and "clearly visible hairy extension of nuclear material" (the criteria of van Diest et al.
    1992 and Baak et al. 1993, as restated by Al-Janabi et al. 2013; summarised by Veta et al. 2016);
  - "short spikes or blunt projections of chromosomes on the surface of aggregated nuclear material", called the most
    reliable feature (VCGP 2021). Veta et al. 2015 found that most of the weaker algorithms' false positives were dark
    objects lacking "protrusions around the edges";
  - imposters are "Smooth contours, circular or oval, hyperchromatic bodies (single or multiple)" (VCGP 2021).

  Sections 5.2, 5.5 and 5.8 target them.

## 2. Your obligations: adversarial checking, of the spec and of yourself

Do these in order, and keep a running log (Section 8, `adversarial_log.csv`). Every row has: `item`, `concern`,
`test`, `result`, `action`, `evidence_path`.

**Before writing any feature code:**
1. For every feature in Section 5, write down, in the notebook, in your own words:
   - what it is supposed to measure, in one sentence a pathologist would accept;
   - two inputs that should score at opposite ends;
   - at least one confound: a different property that would produce the same value for the wrong reason;
   - the test that distinguishes the property from the confound, and what you expect it to show.
   Do this before running anything on real data, so that your expectations are on record.
2. Read Section 5 looking for definitions that would fail your own test. If you find one, do not implement it as
   written. Log it, implement the corrected definition, and keep the spec's version under a `spec_` prefix column so
   that nothing is hidden.
3. Check every constant and every formula for units, scale dependence and window dependence. A feature that changes
   when the ROI is rescaled, or when the 51 px window happens to contain more tissue, is measuring the window, not the
   candidate.

**While building:**
4. Every feature function gets a unit test on hand-built inputs with hand-derived expected values, written before the
   function.
5. Every geometric feature gets the synthetic shape suite in Section 5.2, with the expected ordering written first.
6. Every mask gets a bounds check: physical size, connectivity, containment, and the share of candidates on which it
   touches its bound. Report the share; above 10 % the mask is not doing its job and you say so.
7. Smoke-run on 4 ROIs (`CANDREFINE_SMOKE_ROIS=337.tiff,509.tiff,201.tiff,211.tiff`; `211` exercises the in-mask
   pigment feature of Section 5.8), then read your own notebook top to bottom
   as a reviewer would, before the full run.

**After the full run:**
8. Produce the validity report (Section 6) and read it. For every feature, write one sentence: does the evidence say it
   measures what Section 5 claims? If not, say what it measures instead. Do not soften this.
9. Open every montage with the Read tool and look at it. For every feature, write one sentence on whether its ten
   highest and ten lowest tiles look like the property it claims. The F7 failure was visible in one glance.
10. Re-implement five features of your choice a second time, in a deliberately different way (pure numpy where the
    first used skimage, or the reverse), on 300 random accepted candidates, and report the maximum discrepancy.

**The rule for disagreement with this spec.** There are two kinds of failure, treated differently:
- **A gate fails** (Section 7): something in the pipeline does not reproduce. Stop, report, do not adapt.
- **A definition is defective**: it cannot measure what it claims. Fix it, log it, keep the original under `spec_`,
  and carry on. Do not stop the run over a definitional defect, and do not silently keep a definition you have shown
  to be wrong.

If anything in this spec contradicts the repo (`midog_utils`, `DECISIONS.md`), the repo wins, and you log it.
One exception: the floor clip. It was adopted by `tm_border_floor_arm_verdicts.ipynb` (2026-10-01) and is recorded as
the 2026-10-02 amendment to D3 in `DECISIONS.md`, but it is not yet in `midog_utils.production`. Its absence from
`midog_utils` is not a contradiction to resolve in the repo's favour; the search code comes from the notebooks named in
Section 3.

## 3. Inputs, the candidate set, and what is held fixed

**Held fixed, copied verbatim** (nothing here is re-tuned):
- The 49 ROIs: `images/extra_valid/` (14) and `images/extra_valid/testing_set/` (35), 7 per domain.
- **The 294 recorded clicks and their refinement boxes.** Take the `arm == 'A0'` rows of
  `results/tm_border_floor_49roi_6seed_per_run.csv` with `rank_key == 'tm_score'`: `subset, domain, file_name,
  seed_index, seed_ann_id, click_cx, click_cy, tcx, tcy, bs, mpp, match_radius_px`. This is exactly how
  `tm_floor_definition_49roi_6seed.ipynb` cell 3 builds `CLK`. Assert that they are the same 294 clicks as
  `results/precision_at_k_49roi_6seed_fillgeom_postD11_clicks.csv`.
- **The search:** one new deep search under the adopted A3 settings, specified in Section 3.1. It is the only search
  in this run.
- `refine_template_box` and its constants, verbatim from
  `production_hematoxylin_only/multiotsu_pad_arms/multiotsu_template_refinement_14roi.ipynb` at commit `c767898`
  (`PAD_PX = 2`; no `PAD_UM`). Record its source SHA-1; it must equal `21c43a20db00dccf7ac3bab88aee423c345441b6`.
  - Compute it exactly as follows; this reproduces the value from the current notebook (checked 2026-10-02):
    - `s = ''.join(nb['cells'][3]['source'])`;
    - `node` = the `ast.FunctionDef` named `refine_template_box` in `ast.parse(s).body`;
    - `hashlib.sha1((ast.get_source_segment(s, node) + '\n').encode()).hexdigest()`.
- **The search code.** Take it verbatim from `production_hematoxylin_only/multiotsu_pad_arms/tm_floor_definition_49roi_6seed.ipynb`:
  - every `def` in cell 5 (including `rank_by` and `score_top_k`);
  - `search_arm` from cell 11;
  - the constants in cell 1 that those functions read.

  Extract the def blocks at runtime with `ast`, exactly as for `refine_template_box`, and `exec` them. Record each
  block's SHA-1 in the spec JSON.

  The only new code is `detect_deep` and `search_arm_deep` (Section 3.1). Assert with `difflib` that each differs from
  its verbatim original only in the lines Section 3.1 names.
- Matching: `ev.bucket_detections` at `ev.radius_px(mpp)` (7.5 µm, D7), with the click's own annotation dropped.
- Channels: `hem = ch.to_channel(rgb, 'hematoxylin_od')`; `hed = rgb2hed(rgb.astype(np.float32) / 255.0).astype(np.float32)`
  once per ROI; assert `hed[:, :, 0]` equals `hem` bit for bit (gate 4). E is `hed[:, :, 1]`, D is `hed[:, :, 2]`.
  `gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)`.

**The candidate set.** For each click, take the first 300 rows of its deep A3 list (Section 3.1), in descending score.
That is exactly 88,200 candidate rows. The unclipped list had about 83,800 distinct locations; report this run's
count.

Do not drop any row for any reason. Refused candidates stay, with their refusal recorded and their undefined features as
NaN. The refusal policy is a downstream decision.

**Two membership flags:**
- **`in_recorded_pool`.** The prefix of length `n_detections_recorded`, 79 to 100 rows. That is the pool production's
  A3 search ranks, with `MAX_PEAKS` = 100 (D9). Gate 3 asserts it equals the recorded A3 pool.
- **`in_top100`.** `cand_rank < 100`, the depth the 2026-10-01 version of this spec used.

**Why 300.**
- In the five worst ROIs (`211`, `242`, `462`, `524`, `546`), the top 100 holds about 4.2 mitoses per click and the top
  300 about 7.1. Those counts were measured on the *unclipped* deep pool; the clip adds about 1.5 pool mitoses per click.
- At depth 100, validity on exactly the ROIs these features are meant to help would rest on a handful of positives.

**Other inputs:**
- `databases/MIDOG++.json`, read as `_, annotations = ds.load_annotations(path)` (it returns `(images, annotations)`),
  then `ds.image_annotations(annotations, file_name)` with both categories.
- The clicks and boxes, `results/tm_border_floor_49roi_6seed_per_run.csv`.
- The floors, `results/tm_border_floor_49roi_6seed_floors.csv` (`b_otsu` and `b_otsu_upper` per ROI).
- The recorded A3 run, for gates 2–3:
  - `results/tm_floor_definition_49roi_6seed_per_run.csv`, `arm == 'A3'`: `tpl_sha1`, `threshold`, `n_detections`;
  - `results/tm_floor_definition_49roi_6seed_pool.csv.gz`, `arm == 'A3'`: the recorded pools.
- The regression reference for gate 7, `results/candrefine_rerank_mo_pad61_select_pool.csv`.

Read every CSV with `float_precision='round_trip'`. Record the SHA-256 prefix of every input.

### 3.1 The candidate search: A3, with no peak cap, top 300 kept

**What A3 is.** The recorded `mo_pad61` search with two changes from `tm_border_floor_arm_verdicts.ipynb`
(2026-10-01). They have different status:
- **the floor clip** (arm A2) is the adopted production change, recorded as the 2026-10-02 amendment to D3;
- **the 30 px border band** (arm A1) is a benchmark-harness rule for ROI crops, not a production change. The product
  reads real pixels past any frame edge, so there is nothing for it to remove there. In the harness it makes 28 of
  2,782 annotated mitoses (1.0 %) unreachable.

Everything is exactly as `tm_floor_definition_49roi_6seed.ipynb` runs its `A3` arm.
- **The floor.** Per ROI, `b = roi_floors(hem)['b_otsu']`: the lower 3-class multi-Otsu threshold of the whole
  unclipped `hem`. Assert it equals `floors.csv`.
- **The search channel:** `clip_floor(hem, b)`, which is `np.maximum(hem, b)`. Its `channel_edge` is
  `padded(channel, EDGE_READ_PAD)`, as in cell 11.
- **The template:** the `mo_pad61` template from `build_mo_templates`, cut from the clipped channel at the recorded box.
- **The D11 blank:** the 51 px square around the click, which fills with `b`.
- **The score and peaks:** `TM_CCOEFF`; peak threshold median + `DEEP_FLOOR_Z`·MAD; `PEAK_MIN_DISTANCE`.
- **The 30 px border band:** `band_mask(hem.shape, PAD)` with `PAD` = 30, made invalid and its scores floored
  (`BAND_MODE = 'floored'`).
- **NMS** at `ev.radius_px(mpp)`, then ranking by score.
- **`od`:** computed on the unclipped channel padded by `OD_PAD`, as `detect` does.

**Why the clip exists.**
- Under `TM_CCOEFF` the filled template is mostly a mass detector: besides a shape term, the score grows with mean H
  in the template's kept box minus mean H in its filled surround. For a uniform change, one unit of surround lightness
  counts as much as one unit of box darkness, so dark objects beside lumen, glass or pale stroma were over-rewarded,
  whatever they were.
- Flooring the channel at `b` caps that reward at `b`, and it also lowers the template's mass weight (median 0.040 →
  0.035). `b_otsu` is the whole lowest 3-class Otsu class, not glass alone: it clips 26.6–82.5 % of an ROI's pixels
  (mean 55 %), against a white share (gray > 215) of mean 9 % and median 4.9 %. A floor at 0.25 × `b_otsu` (glass plus
  the palest tissue) recovers only about a quarter of the gain, so the nuisance is pale surround broadly.
- On A3 − A1 (the clip on the band-corrected harness): white-adjacent top-10 false positives fell 316 → 240, other
  false positives 1019 → 998 (n.s.), mitoses left the top 10 at 8.6–10.3 % in every white-share bin up to 40 % (2.7 % above 40 %, n = 37), and
  P@10 at `tm_score` rose 3.30 points.
- **Consequence for this spec.** The features read the unclipped channel (Section 4), so the surround columns
  (`od_contrast`, the ring columns, the white shares) carry the same white-surround signal the clip removed from the
  score. That is correct for measurement, and not a reason to clip the features. But no validity column isolates it
  yet (a white-stratified AUC was proposed and not adopted). Read every surround column beside its correlation with
  `white_frac_gray215_61` in `redundancy.csv`; a correlation shows overlap but does not isolate the effect.

**The one change: no peak cap, then keep 300.**
- `detect_deep` is `detect` with exactly these changes:
  1. `tm.extract_peaks(..., UNCAPPED)` with `UNCAPPED = 10**7`, then assert `len(centers) < UNCAPPED`.
  2. `keep = nms_by_distance(centers, scores, nms_radius)` over all peaks. Also
     `keep_rec = nms_by_distance(centers[:prod.MAX_PEAKS], scores[:prod.MAX_PEAKS], nms_radius)`, and assert
     `np.array_equal(keep[:len(keep_rec)], keep_rec)`. Greedy NMS runs in score order, so the 100-peak pool is the
     prefix. `deep_pool_mo_pad61_49roi_6seed.ipynb`'s `search_deep` makes the same assertion.
  3. Truncate to the first `DEPTH = 300` post-NMS rows **before** computing `od`, and compute `od` only on them. `od`
     costs about 6 ms a row, and an uncapped pool is about 15k.
  4. Return `n_peaks` (uncapped), `n_detections_deep` (the full post-NMS count), `n_detections_recorded`
     (`len(keep_rec)`) and `threshold`.
- `search_arm_deep` is `search_arm` with its `detect(...)` call replaced by `detect_deep(...)`. Nothing else changes.
- **Cost:** about 45 min for the 294 searches; the unclipped deep pool took 45 min. Add the features' time to that.

**There are no border candidates.** The band is exactly the pixels with `edge_px < 30`, where
`edge_px = min(cx, cy, w − 1 − cx, h − 1 − cy)` (`assert_band_identity`). The band is invalid in the search, so every
candidate has `edge_px ≥ 30`.
- Its 51 px refinement window and 61 px white window therefore lie entirely inside the ROI.
- Assert, on every row:
  - `edge_px ≥ 30`;
  - no refinement refusal of class `roi border`.

  Either violation means the band was not applied. Halt.
- There is no padded re-run of any candidate and no `edge_read` column.
- The search's own internal pads are part of the verbatim A3 search and stay:
  - the correlation pad;
  - `OD_PAD`;
  - the template read pad.

  None of them supplies a candidate's feature pixels.

## 4. Per-candidate pipeline

**Refinement at the candidate.** Run the verbatim `refine_template_box(hem, cx, cy, mpp)` at the detection's own
coordinates, on the **unclipped, unblanked** `hem`. Every feature reads that channel too.
- The clip and the blank belong to the search: they decide which candidates exist.
- Production's own box refinement also reads the unclipped channel.
- The surround features need the real background pixels that the clip would flatten.

Record:
- `accepted`;
- the refusal class (`roi border`, which must never occur here, `no piece near centre`, `selection too small`,
  `selection too large`, `empty chromatin mask`);
- `mode`, `n_pieces`, `base_size`, and the box centre's offset from the detection.

Features depend only on the location, so compute them once per distinct `(file_name, cx, cy)` and join to clicks.

**Windows near the ROI edge.** With the band out, the 51 px and 61 px windows always lie inside the ROI (Section 3.1).
- The 101 px windows can still cross the edge, for candidates 30–50 px from it. They are **truncated** to the pixels
  inside the ROI, never padded.
- Store `win101_area_share`, the real pixels divided by 101². Expect it at 0.80 or above at one edge, and 0.64 or above
  at a corner.
- Every 101 px feature is computed on the truncated window. Report its validity split at `win101_area_share == 1`.

**Masks.** Three masks, all in the 51 px window, all stored by size in pixels and µm²:
- `M_chr`: the undilated top multi-Otsu class inside the selection, `selection_mask & (class_map == 2)`. Boundary
  features use this mask, because the refinement's 1 px dilation fills exactly the notches a boundary feature measures.
- `M_sel`: the dilated selection mask, as the template uses it. Kept for continuity with the earlier `solidity`.
- `M_nuc`: a **bounded nucleus mask** that does not depend on the window's own histogram. Default definition, which you
  must test and may replace with a logged reason:
  - `t_nuc` is one threshold per ROI, the lower threshold of `threshold_multiotsu(hem, classes=3)` computed over the
    whole ROI channel, so a dense window cannot lower it. Store `t_nuc` per ROI and the share of the ROI above it; a
    share outside 10 to 60 % is a flag.
  - **`t_nuc` is the search floor.** It is the same call as `roi_floors(hem)['b_otsu']` (Section 3.1). Assert
    `float(np.float32(t_nuc)) == b_otsu` from `floors.csv` on every ROI. So "above `t_nuc`" means "untouched by the
    clip" (the clip keeps `hem >= b`, this mask takes `hem > t_nuc`; they differ only at equality), and "below `t_nuc`"
    lumps glass, lumen, halos and pale cytoplasm together.
    Expect the share flag to fire on `246`, `293` and `201` (65–73 % of pixels above).
  - `c` is the centroid of `M_chr`, `r_eq` its equivalent radius, `R_nuc = r_eq + 2.5 µm`.
  - `M_nuc` is the union of the connected components (connectivity 1) of `hem_window > t_nuc` that intersect `M_chr`,
    intersected with the disc of radius `R_nuc` about `c`, unioned with `M_chr`. The centroid `c` may fall outside
    `M_chr` for a ring- or bar-shaped figure, which is why the rule is "intersects `M_chr`" and not "contains `c`".
  - Store `nuc_clipped` (the component reached the disc), `nuc_touches_window` (the mask reached the window edge), and
    `nuc_area_um2`. Section 6 requires the clipped share and every nucleus-based feature split by it.
- Ring regions for the surround features: `ring_k` for k = 1, 2, 3 is the band of pixels whose distance from `M_chr`
  (Euclidean distance transform of `~M_chr`) lies in `((k−1) µm, k µm]`, clipped to the window and excluding `M_chr`.
- **Dark pixels on the chromatin are mostly not saturated, but some are clipped.** Measured 2026-10-02:
  - The `M_chr` pixels with any channel below 20 are G < 20 with mean RGB about (27–36, 9–16, 47–88) on the GT mitoses
    of `425`, `527` and `300`: dense blue-purple hematoxylin. No mitosis pixel on the recorded pool has all three
    channels below 20. Among other candidates such pixels occur on 7 ROIs, mostly `211`'s black pigment (4.8 % of its
    candidates' `M_chr` pixels).
  - True clipping does occur. On the 3D Histech ROIs `249`, `435`, `230` and `460`, 5.7–7.5 % of GT-mitosis `M_chr`
    pixels have G at exactly 0, with about twice as many at 0 as at 1 (a pile-up); `425` has 1.3 % and `300` none.
    `rgb2hed` clamps 0 to 1e-6, so G going from 1 to 0 moves H by −0.04 and E by +0.68: those pixels' colour and local
    contrast are artefacts.
  - On the recorded A3 pool (seeds 0–1), excluding every pixel with a channel below 20 dropped 60 % of a mitosis's
    `M_chr` against 37 % of other candidates', and left fewer than 10 pixels for a support eroded by a 3 × 3 square on
    65 % of mitosis rows against 39 % of others. That is a subset chosen by label and darkness, so no chromatin colour
    or texture column of Sections 5.4 and 5.6 uses it.
  - `chr_clip_share`: the share of `M_chr` pixels with any RGB channel at 0 or 255. Store it. It exceeded 5 % on
    10.8 % of mitosis rows and 7.0 % of others on the recorded pool.
  - **Colour** columns on the chromatin read `M_chr_clip` = `M_chr` minus those clipped pixels, and are NaN with
    `nan_reason` `clipped` when |`M_chr_clip`| < `NS_MIN_PX` = 10.
  - **Texture** columns exclude clipped pixels too, as Section 5.4 says; store `tex_support_px`.
  - The darkness confound of both is read from `auc_dark_stratified` and `rho_with_h_chr_mean`; the clipping effect
    from `auc_split_chr_clip` (Section 6).
  - Section 5.8 (a)'s `chr_resid_frac` keeps [20, 235] on `M_res`, because black pigment is low in all channels. Its
    support is therefore label-dependent too (median `chr_ns_share` 0.76 for mitoses against 0.44 for top-10 false
    positives on 7 ROIs, third review's measurement), so report its AUC beside `chr_ns_share`.

**Labels.** Match each click's 300-row list in score order with `ev.bucket_detections` against the ROI's ground truth
minus the click's own annotation. Store `bucket`, `matched_ann_id`, `is_mitosis`. Gate 3 requires the recorded-pool
prefix's labels to equal the recorded A3 pool's. Independently of matching, store
`nearest_mitosis_um` and `nearest_lookalike_um`, the distance to the nearest annotation of each category, excluding
the click's own; the downstream analysis needs labels that do not depend on greedy order.

## 5. Features

**Rules that apply to every feature:**
- Store the **raw measurement**. Never negate, rank, standardise or direction-encode a column. Directions are
  hypotheses and live in `directions.csv` (Section 8), written before the full run.
- Physical units where a quantity has one: µm, µm², optical density. Nothing may depend on `mpp` implicitly.
- Nothing may depend on the 51 px window's own intensity histogram except through `M_chr` and `M_sel`, which are the
  production refinement and are kept for that reason.
- NaN where undefined, with the reason recorded in `nan_reason` for that family; never a sentinel number.
- Column names are `snake_case`, and every column is described in the notebook in one line with its unit.

### 5.0 What pathologists use, and which columns target it

Use this table as the `rationale` of `directions.csv`. A direction given here is the hypothesis to write down; a 0
stays 0 unless you log a reason before the full run. Quotes are verbatim from the sources listed after the table.

| # | criterion | evidence | columns | direction |
|---|---|---|---|---|
| C1 | **Short spikes or rods on the surface of the chromatin**, the most reliable cue | VCGP 2021: "The most reliable microscopic feature to identify MF and AMF are the short spikes or blunt projections of chromosomes on the surface of aggregated nuclear material"; also "dark purple staining nuclear aggregates with an irregular contour (short rods or spikes of chromosomes, described by van Diest as hairy projections)". Meuten 2021: "short 'rods' of chromosomes protruding from the surface of aggregates of nuclear material". Al-Janabi 2013 (restating the van Diest and Baak criteria): "clearly visible hairy extension of nuclear material". Veta 2015, on weaker algorithms' false positives: "dark objects that lack other characteristics of mitotic figures such as protrusions around the edges". VCGP QRG 2022, under "Differentiate MF from Mitotic-Like-Figures (MLF)": "Surface contour: MLF are smooth, MF have projecting rods/spikes" | `chr_radial_hf`, `chr_tophat_edge`; `chr_roughness` and `chr_concavities_per_um` only if they pass the lobulation test (5.2) | +1 |
| C2 | **No nuclear membrane** | Al-Janabi 2013: "absence of the nuclear membrane". VCGP 2021: "lack of a nuclear membrane". Veta 2015: "hyperchromatic objects that lack a clear nuclear membrane" | `rim_core_nuc`, `ring_peak_prominence`, `chr_hole_frac`; Section 5.5 says what each can see | per column, 5.5 |
| C3 | **Dark chromatin**, necessary but not sufficient | VCGP 2021: "dark purple staining nuclear aggregates". Veta 2015: "many other cellular components can have a similar hyperchromatic appearance" | Section 5.3 | +1; it is why `auc_dark_stratified` exists |
| C4 | **Configuration**: a clot or ring, a plate, separate clots | VCGP 2021: "Prometaphase: A dark cluster or ring of chromosomes with variably protruding chromosome rods and spikes"; "Metaphase: An equatorial plate (linear or ring-shaped aggregate with short rods or spikes)"; "Anaphase: Two separate aggregates, variable distances apart, linear with projections/spikes"; "Anaphase and telophase MF are counted as one". Al-Janabi 2013: "either clotted (beginning metaphase), in plane (metaphase/anaphase), or in separate clots (telophase)". Meuten 2021: identifying the phase is "not necessary" | `chr_elongation`, `chr_hole_frac`, the `pair_*` columns (5.2) | 0 for all (below, and 5.2) |
| C5 | **Imposters are smooth, round or oval dark bodies, one or several** | VCGP 2021, mitotic-like figures: "Smooth contours, circular or oval, hyperchromatic bodies (single or multiple)". ICCR, in its own words (its method source is Elston & Ellis 1991): "hyperchromatic, karyorrhectic, or apoptotic nuclei are excluded". Cree 2021: "Strict criteria should be applied for accepting a mitotic figure in distinction from a pyknotic or apoptotic nucleus" | `chr_solidity`, `chr_circularity`, `chr_edge_sharpness_per_um` | −1 for solidity and circularity; 0 for edge sharpness (5.8 c) |
| C6 | **Prophase is not counted on H&E** | VCGP 2021: "MF are counted from prometaphase through telophase … Prophase cannot be reliably identified in histological preparations". Lin 2024: dark chromosome aggregates "are most likely in the final stages of prophase, commonly referred to as prometaphase". Meyer 2005 recommends "that prophases be excluded from mitotic counts performed for cancer grading", though its own study counted them. MIDOG++ delegates its criteria to Donovan et al. 2021 and Ibrahim et al. 2022. MIDOG 2022's PHH3-assisted ground truth (not the H&E labels) did not annotate immunopositive cells lacking H&E morphology, "mostly early prophase MFs", and its largest subcategory of figures found only with PHH3 (within "borderline morphology", 87.58 %) was "prophase with resemblance to apoptosis" (31.23 %) | none directly; context for `ring_peak_prominence` (5.5) | — |
| C7 | **Named imposters**: immune cells, stromal cells, two cells in close proximity, karyorrhectic debris | Lin 2024, on which mimic each true mitotic figure resembled (most often the low-agreement ones), not on false positives: "karyorrhectic debris, immune cells, stromal cells, and two cells in close proximity"; "lymphocytes were most likely to resemble mitosis in prometaphase, and stromal cells were likely to mimic mitosis in metaphase" | `chr_area_over_nbr` (5.8 d); the lobulation test (5.2) | 0 |
| C8 | **Pigment** | Aubreville 2023 (MIDOG++): melanoma "pigment particles that contribute additional imposter structures ... although they have a different chromaticity (brown)" | Section 5.6; `chr_resid_frac` | −1 |
| C9 | **Cytoplasm and surround** | Ikpatt 2002, restating the criteria of Baak & Oort: "The cytoplasm of the mitotic cells was often larger during mitosis than in the resting cells". No verified source states a clear zone around mitotic figures in H&E. A clear halo is described around apoptotic hepatocytes in rodent liver (NTP atlas). VCGP 2021, mitotic-like figures: "MLF cytoplasm tends to be eosinophilic, however color is not reliable". VCGP QRG 2022: "Cytoplasm color (not reliable). MF tend to have more amphophilic/basophilic cytoplasm (increased RNA), MLF tend to have more eosinophilic cytoplasm" | lightness and H: the ring columns and white shares (5.5). Cytoplasm colour: `e_ring23_mean`, `e_over_h_ring23` (5.5) | 0 |
| C10 | **Context** | VCGP 2021: "A structure within a region of necrosis ... is more likely to be MLF"; "A structure within a region of high mitotic activity is more likely to be a MF". VCGP QRG 2022: "Context: In areas of poor cell preservation, MLF may be more likely" | none | not measured |

**Elongation.** A side-on metaphase plate is elongated, but so are spindle-cell and stromal nuclei, which Lin 2024 names
as metaphase mimics; `462` and `410` are sarcomas. Pathologists need elongation together with surface spikes, so
`chr_elongation` and `chr_eccentricity` have direction 0.

**Texture and colour.** VCGP 2021, in its list of mitotic-like figures (item 5), describes resting-cell chromatin as
having "more heterochromatin with a more clumped chromatin pattern", heterochromatin "staining as dark granules", but
gives no texture criterion for mitotic figures, so the texture columns (5.4) have direction 0. VCGP's "color is not reliable" is about mimic *cytoplasm*; no source
gives a chromatin colour criterion, so the chromatin colour columns (5.6) have direction 0 for lack of a source.
Pigment (C8) is the exception.

**Independence.** VCGP 2021, Meuten 2021, Donovan et al. 2021, MIDOG 2022 and MIDOG++ share authors; the VCGP QRG
2022 names no authors and cites the first three. MIDOG++ also includes van Diest and Veta, who co-authored Al-Janabi
2013 (van Diest), Veta 2015 (van Diest) and Veta 2016.
So every C1 and C2 source shares an author with the dataset's group; Lin 2024 (C6, C7) does not. Their agreement is
one school's view, not independent corroboration.

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

### 5.1 Identity and pipeline columns

`subset, domain, file_name, image_id, mpp, seed_index, seed_ann_id, click_cx, click_cy, cand_rank, in_recorded_pool,
n_detections_recorded, cx, cy, score, od, dist_to_click_um, edge_px, win101_area_share, accepted, refusal, mode, n_pieces,
base_size, box_offset_px, bucket, matched_ann_id, is_mitosis, nearest_mitosis_um, nearest_lookalike_um, in_top100,
scanner, roi_group`.

**`scanner`** comes from `images/extra_valid/MANIFEST.md` and `images/extra_valid/testing_set/MANIFEST.md`:

| scanner | ROIs |
|---|---|
| `3dhistech` | canine lung cancer, canine lymphosarcoma, canine soft tissue sarcoma |
| `aperio_cs2` | canine cutaneous mast cell tumor; breast `123`, `128` |
| `hamamatsu_s360` | breast `072`, `094`, `100` |
| `hamamatsu_xr` | human melanoma, human neuroendocrine tumor; breast `007`, `013` |

Assert that all 49 ROIs get a scanner, and that the table matches the manifests.

**`roi_group`** is `failing9` for `211`, `242`, `524`, `546`, `462`, `527`, `410`, `245` and `289`, and `other40`
otherwise.
- The nine are the ROIs below P@10 0.35 on both A3 axes (`low_precision_rois/`).
- The group was defined after seeing results, so every split by it is descriptive.

`n_pieces` is the number of refinement pieces the selection unioned (`len(ref['selected'])`), distinct from `chr_n_pieces` in Section 5.2.
`od` comes from `detect_deep`, which reads the unclipped channel. Recompute it with `cm.score_detections` on the
unclipped `hem`, for the first 2 clicks of every ROI, and assert equality at `rtol 1e-12` (gate 5).

### 5.2 Geometry of the chromatin, on the undilated `M_chr`

Computed on the largest 4-connected piece of `M_chr` unless marked *whole*. NaN below 30 px of area. The outline
columns (`chr_major_um` through `chr_concavities_per_um`, and `chr_radial_hf`) read that piece with its holes filled
(`scipy.ndimage.binary_fill_holes`). Pathologists judge the outer surface contour (C1, C5), and on a ring-shaped piece
`regionprops` counts the inner boundary too: the vesicular nucleus of Section 5.5, whose `M_chr` is its rim, scored
solidity 0.21, circularity 0.06 and roughness 1.82 unfilled, the most "mitotic" values, against 0.99, 0.95 and 1.00
filled. `chr_area_um2` and `chr_hole_frac` read the piece unfilled; the hole is `chr_hole_frac`'s job.

| column | definition | unit |
|---|---|---|
| `chr_area_um2` (*whole*) | pixel count × mpp² | µm² |
| `chr_n_pieces` (*whole*) | 4-connected pieces of `M_chr` | count |
| `spec_chr_two_piece_sep_um` | centroid distance of the two largest pieces of `M_chr`; NaN if one piece. Kept for comparison only; the pair columns below replace it | µm |
| `chr_major_um`, `chr_minor_um` | second-moment ellipse axes (`regionprops` major/minor axis length × mpp) | µm |
| `chr_elongation` | `chr_major_um / chr_minor_um` | ratio |
| `chr_eccentricity` | `regionprops` eccentricity | ratio |
| `chr_solidity` | area / convex-hull area | ratio |
| `chr_circularity` | 4π · area / perimeter² with `regionprops` perimeter | ratio |
| `chr_roughness` | perimeter / convex-hull perimeter | ratio |
| `chr_concavities_per_um` | `cv2.convexityDefects` deeper than `CONCAVITY_MIN_UM` = 0.25 µm (the returned depth is fixed-point: divide by 256 for pixels, then multiply by `mpp`), divided by the perimeter in µm | per µm |
| `chr_hole_frac` | 1 − area / area with holes filled (`scipy.ndimage.binary_fill_holes`), on the largest piece. A ring-shaped `M_chr` scores high; Section 5.5 says why that matters | ratio |
| `chr_radial_hf` | radial outline profile at 128 equal-arc points about the centroid, mean-centred; Σ over harmonics 6 to 64 of \|rfft\|² divided by (128 · mean radius²). The outline is the longest contour of `cv2.findContours(..., cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)`, closed and resampled at equal arc length by linear interpolation; the centroid is the region's pixel centroid. A `skimage` contour gives different values (disc 0.047 against 0.057 at r = 8), so use this one | ratio |
| `sel_solidity` | `regionprops(M_sel.astype(np.uint8))[0].solidity`, the earlier F3, unchanged | ratio |
| `spec_spiky` | the earlier F7, verbatim: high-harmonic energy over total outline energy; kept for comparison only | ratio |

`chr_radial_hf` replaces the earlier `spiky`. Its denominator is the mean radius squared, not the total outline
energy, so a disc scores near zero instead of near one. Keep the earlier definition as `spec_spiky` for comparison.

**Pair columns: anaphase and telophase (C4), read from a 101 px window, not from `M_chr`.**
- **Why not `M_chr`.** `M_chr` is restricted to the refinement's selection. In base mode that is the one piece under
  the detection; in fallback mode it is the pieces centred within 3.75 µm. So a pair seen through `M_chr` depends on
  where the detection landed. On a synthetic anaphase pair 4.5 µm apart, the verbatim refinement gives
  `spec_chr_two_piece_sep_um` 4.5 µm when the detection sits in the gap and NaN when it sits on a plate. It also gives
  0.9 µm on an interphase nucleus with a nucleolus.
- **Why 101 px.** A 51 px window centred on one plate does not contain a partner much more than 6 µm away. The 101 px
  window is the one Sections 5.3 and 5.8 (d) already use, truncated at the ROI edge.
- **Pieces `P`:** 4-connected components of `hem_101 > t_dark` (the ROI-wide upper threshold, Section 5.7), each at
  least `PAIR_MIN_UM2` = 1.0 µm² (the refinement's speck floor). Thresholding at the ROI-wide `t_dark` keeps the
  columns off the window's own histogram.
- **Anchor:** the piece of `P` that overlaps `M_chr` the most; ties go to the larger piece. NaN with `nan_reason`
  `no anchor` if no piece overlaps `M_chr`. Place `M_chr` in the 101 px window at row and column offset
  `(cy − 25 − y0, cx − 25 − x0)`, where `(y0, x0)` is the truncated window's origin in ROI pixels. That is 25 only when
  the window is not truncated; about 1.6 % of rows (`edge_px` 30–49) are truncated.
- **Partner:** among the other pieces with area ≥ `PAIR_MIN_RATIO` = 0.5 × the anchor's and centroid distance ≤
  `PAIR_MAX_UM` = 10 µm, the nearest by centroid. NaN with `nan_reason` `no partner` if none.
- **Columns:**
  - `pair_sep_um`: centroid distance;
  - `pair_area_ratio`: smaller area over larger;
  - `pair_parallel_deg`: the acute angle between the two pieces' major axes. With θ = `regionprops` `orientation`,
    a piece's major axis in (row, column) coordinates is `(cos θ, sin θ)`; the centroid line is
    `centroid_partner − centroid_anchor` in the same (row, column) order;
  - `pair_axis_deg`: the mean, over the two pieces, of the acute angle between the piece's major axis and the line
    joining the centroids. Separating plates seen side-on score near 90°.
  - The two angle columns are NaN (`nan_reason` `round piece`) when either piece's major/minor axis ratio is below
    `PAIR_ELONG_MIN` = 1.2, where orientation is undefined.
  - `pair_touches_window`: either piece touches the truncated window's edge.
- **What the columns cannot do.** A neighbouring dark nucleus of similar size is a partner too, so in crowded tissue
  most candidates have one. And side-by-side spindle nuclei in a fascicle give the same angles as separating plates
  (parallel near 0°, axis near 90°); `462` and `410` are sarcomas. So every pair column has direction 0; the
  hypothesis that `pair_axis_deg` near 90° marks anaphase goes in the rationale only.
- **Pair test, mandatory, written first.** Build the synthetic anaphase windows of the reference construction below
  and run the verbatim refinement on them.
  - The pair columns must agree, within 0.25 µm and 5°, between a detection at the gap's midpoint and one at a plate's
    centre, at 4.5 and 7 µm centre-to-centre. At 9 µm only the plate-centred detection is tested: the verbatim
    refinement refuses the gap-centred one (`no piece centre within 3.75um of the click`), which is its own rule, not
    a pair-column failure.
  - **Oblique pairs.** The same pairs rotated by 30°, 60° and −30° about the canvas centre must give `pair_axis_deg`
    at least 85° and `pair_parallel_deg` at most 5°. An axis-aligned pair cannot catch a sign error in the orientation
    convention; at a 30° rotation a flipped sign gives about 30° instead of 90°.
  - The interphase nucleus with a nucleolus, the uniform dark disc and the metaphase plate must give `no partner`.
  - Checked values, from two independent implementations of this text, identical for noise seeds 0 to 3: 4.5 / 7.0 /
    9.0 µm, area ratio 1.00, parallel 0°, axis 90°, the same at both centrings; `no partner` on the other three. Rotated
    (one implementation): separation 4.49–4.50 / 7.02 µm, parallel 0°, axis 89.4–90.0°.
  - Construction: `mpp` 0.25; a 201 × 201 canvas at 0.03 OD plus Gaussian noise of σ 0.003 (any seed); then, drawn
    over it with constant values that replace the noise, a cell disc of radius 30 px at 0.045 OD and two plates, each
    7 × 29 px (`|x′ − x_k| ≤ 3`, `|y′| ≤ 14` in coordinates rotated about the centre) at 0.20 OD, centred 18, 28 or
    36 px apart along x′; then `scipy.ndimage.gaussian_filter` with σ 1 px, cast to float32. `t_dark` is the upper
    threshold of a 3-class multi-Otsu over the whole canvas. The other three windows are the membrane cases of
    Section 5.5, and a 9 × 33 px plate at 0.20 OD inside a clear disc of radius 22 px at 0.02 OD inside a cell disc of
    radius 34 px at 0.045 OD, drawn the same way.

**Synthetic suite, mandatory, expected ordering written before running.** Rasterise at two scales (radius r = 8 px and
16 px) each of the shapes below, centred on a (6r + 1)² canvas, with x and y the centred coordinates, ρ = √(x² + y²)
and θ = atan2(y, x):
- a disc, ρ ≤ r; a 2:1 ellipse with semi-axes r√2 (along x) and r/√2; a 4:1 ellipse with semi-axes 2r and r/2;
- 5-, 8- and 12-lobed stars, ρ ≤ r·(1 + 0.25·cos(nθ));
- a 2:1 ellipse with 12 small spikes: the 2:1 ellipse plus, at each angle φ = 2πk/12, the rectangle
  `e(φ) − 1 ≤ t ≤ e(φ) + 4`, `|o| ≤ 1`, where e(φ) is the ellipse's radius at φ and t, o are the along- and
  across-ray coordinates. The spikes are a fixed 2 × 5 px, not scaled with r;
- two touching discs, each of radius 0.9r, centred at x = ±0.9r;
- a disc with one deep notch: the disc minus the slot x ≥ 0, |y| ≤ max(1, r/8).

For each boundary feature, the suite must show:
- discs and ellipses score alike and low on `chr_radial_hf`, `chr_roughness` and `chr_concavities_per_um`, so none of
  them is secretly an elongation measure. Expect `chr_radial_hf` to show a partial elongation sensitivity: the 4:1
  ellipse scores 0.097 / 0.081 against the disc's 0.057 / 0.010, and equal-arc sampling gives an elongated outline
  in-band energy even without pixels. Report it as such;
- the 8- and 12-lobed stars and the spiked ellipse score clearly above every disc and ellipse on `chr_radial_hf` at
  each r (the checked values below put the closest pair, the spiked ellipse against the 4:1 ellipse at r = 16, at
  5.6 ×; that ratio was read off the values, not set beforehand). Not "increasing with lobe count": at a fixed
  relative amplitude a star's in-band outline energy barely depends on its lobe count. By Parseval it is exactly the
  same (2.0) for any count inside the 6–64 band under equal-angle sampling, and about 1.4–1.6 for 6–24 lobes under
  this column's equal-arc sampling. Pixelation then reorders them at r = 8. The 5-lobed star's fundamental lies below
  the band, so it scores low;
- **the lobulation test.** Two touching discs and the notched disc are smooth outlines, not spikes: they stand in for
  "two cells in close proximity" (Lin 2024) and for the crowded-nuclei false positives of Section 1. Every column that
  claims to measure surface spikes must score both below the 8- and 12-lobed stars and the spiked ellipse. A column
  that fails is relabelled as a lobulation measure in `directions.csv` and the log, and kept. No expectation is
  recorded for `chr_roughness` or `chr_concavities_per_um`;
- **scale.** The two scales agree within 10 % of the larger value for the uniformly scaled shapes (disc, ellipses, stars, touching discs,
  notched disc). Expect `chr_radial_hf` to fail it on every shape (by the checked values below, the closest is the
  touching discs at 12 %): its low values are mostly pixel-staircase noise, which falls with r², and 8 or 12 lobes on an
  8 px radius are 2–3 px wide, at the pixel pitch. So at these sizes it is partly a size measure. Report those as failures caused by
  the resolution limit, with the radius below which it applies; do not widen the tolerance. The spiked ellipse is
  exempt, because its spikes keep a fixed pixel size and so shrink relative to the body as r grows, as real
  chromosome arms do on a larger figure. The consequence: fixed-size spikes score 3 × lower on a body twice the size
  (1.346 → 0.453), so `chr_radial_hf` is negatively size-confounded on its own cue. Section 6's `auc_area_tertiles` is
  the real-data check and is essential for it;
- `chr_elongation` separates discs from ellipses and nothing else does.

Checked `chr_radial_hf` values (one implementation, r = 8 / 16): disc 0.057 / 0.010; 2:1 ellipse 0.033 / 0.008;
4:1 ellipse 0.097 / 0.081; 5-, 8- and 12-lobed stars 0.024 / 0.020, 2.328 / 1.681 and 1.089 / 1.635; spiked ellipse
1.346 / 0.453; touching discs 0.180 / 0.158; notched disc 0.075 / 0.035. An independent implementation reproduced the
disc and 12-lobed star values to the third decimal.

Store the results in `synthetic_tests.csv`. Any feature that fails its own expected ordering is reported as failing
and still stored; do not quietly tune the harmonic range or a tolerance until it passes.

### 5.3 Darkness and intensity

| column | definition | unit |
|---|---|---|
| `h_chr_mean`, `h_chr_p90`, `h_chr_max` | H over `M_chr` | OD |
| `h_sel_mean` | H over `M_sel` | OD |
| `h_nuc_mean` | H over `M_nuc` | OD |
| `h_ring1_mean`, `h_ring2_mean`, `h_ring3_mean` | H over `ring_1`, `ring_2`, `ring_3` | OD |
| `h_win51_mean` | H over the whole 51 px window | OD |
| `h_win101_mean` | H over the 101 px window, truncated at the ROI edge (Section 4) | OD |
| `od_contrast` | `h_chr_mean − h_ring3_mean` | OD |
| `dark_frac_chr` | share of `M_chr` pixels with `gray < 40` | share |

`h_chr_mean` is the earlier `core_od`, unchanged in definition. `dark_frac_chr` (gray < 40) is a darkness covariate
for the F9 check in Section 5.6; it is not a saturation measure (Section 4). Expect `h_chr_mean` and `od` to be nearly collinear; Section 6 reports it.

### 5.4 Chromatin arrangement and texture

| column | definition | unit |
|---|---|---|
| `condensed_nuc` | `chr_area_um2 / nuc_area_um2` with the bounded `M_nuc` | ratio |
| `h_chr_cv` | std / mean of H over `M_chr`, ddof 0; the earlier `texture`, kept unchanged as a dispersion measure and for gate 7. It reads saturated pixels | ratio |
| `h_chr_grad` | mean, over the texture support, of the Sobel gradient magnitude of H (`cv2.Sobel`, ksize 3, each component divided by 8 and by `mpp`, so OD per µm), divided by `h_chr_mean`. The texture support is `M_chr` eroded once by `scipy.ndimage.binary_erosion` with its default cross, minus every pixel whose 3 × 3 neighbourhood contains a clipped pixel (Section 4); NaN with `nan_reason` `clipped` below `NS_MIN_PX`. The 3 × 3 kernel still reads diagonal neighbours outside the mask, so a flat disc scores 0.032, not 0 | 1/µm |
| `h_chr_lap` | std, over the same texture support, of the 3 × 3 Laplacian of H (`cv2.Laplacian`, ksize 1, divided by `mpp²`), divided by `h_chr_mean` | 1/µm² |
| `h_chr_glcm_contrast`, `h_chr_glcm_homog` | `graycomatrix` at distance 1, four angles averaged, within `M_chr`'s bounding box, on H / `h_chr_mean` quantised to `GLCM_LEVELS` = 64 equal bins over [0, `GLCM_REL_HI`] with `GLCM_REL_HI` = 2.0: level = `np.minimum(np.floor(H / h_chr_mean / GLCM_REL_HI * GLCM_LEVELS), GLCM_LEVELS − 1)` (rounding instead of flooring gives 30.19, not 23.86, on the checkerboard below); pairs with a pixel outside `M_chr`, or with a clipped pixel (Section 4), excluded (give those pixels a 65th level and drop its row and column from the matrix). Relative bins, like the normalisation of the two rows above: fixed OD bins made contrast scale with OD², within-ROI ρ 0.60 with `h_chr_mean` against 0.10–0.15 for relative bins (two review measurements, on 17 and 49 ROIs, seed 0). Within-object rank levels were rejected because they turn pixel noise on a flat body into maximal texture (noise test below). Pixel noise that is fixed in OD spreads over more relative bins on a pale object | ratio |

**`condensed_disc` is removed.** It was `chr_area_um2 / (π R_nuc²)` with `R_nuc = r_eq + 2.5 µm` and
`r_eq = √(area / π)`, which is `r_eq² / (r_eq + 2.5 µm)²`: a fixed increasing function of `chr_area_um2` (5, 10, 20 and
40 µm² give 0.11, 0.17, 0.25 and 0.35). It carried no information beyond size.

**`condensed_nuc` surround test, mandatory.** `M_nuc` is grown above `t_nuc`, which is the search floor (Section 4),
and about 45 % of a typical ROI lies above it. So around a mitosis, whose chromatin has no nucleus around it, `M_nuc`
takes in any above-floor cytoplasm within 2.5 µm. Fix `t_nuc` at 0.05 OD for this test. Render the 9 × 33 px plate of
Section 5.2's pair test, without the clear disc, inside a cell disc of radius 34 px whose cytoplasm is at 0.04 OD
(0.8 × `t_nuc`), and again at 0.075 OD (1.5 × `t_nuc`). Record `condensed_nuc` for both.
- Expect them to differ. The review's implementation gave 0.77 against 0.24; in above-floor cytoplasm `M_nuc` fills
  the bounding disc, so the column equals the removed `condensed_disc`.
- If they differ by more than 25 %, `directions.csv` describes the column as a surround measure as much as a
  condensation measure.
- On the recorded A3 pool (seeds 0–1, measured 2026-10-02), `nuc_clipped` was 79 %, far past Section 2 step 6's 10 %
  bound. That does not make `M_nuc` the disc: where `nuc_clipped` is true it is the disc ∩ the above-floor tissue
  connected to `M_chr`, equal to the whole disc on 0.2 % of those rows, with a median fill of 58 % of the disc
  (10th–90th percentile 44–79 %, third review's measurement). Report that distribution in the mask bounds check.

**Noise test, mandatory, written first.** On a 71 × 71 canvas at `mpp` 0.25, render a disc of radius 16 px at
0.20 OD on 0.03 OD (`M_chr` = the disc) in four versions: flat; flat plus Gaussian noise of σ 0.003 everywhere; the
noisy one blurred at σ 1 px; and the granular disc of Section 5.8's reference construction (a 2 px checkerboard of
+0.06 OD more than 3 px inside the edge). Every texture column must score the flat, noisy and noisy-blurred discs far
below the checkerboard; a column that does not is measuring noise or the blurred edge. Checked values, one
implementation (flat / noisy / noisy-blurred / checkerboard):
- `h_chr_glcm_contrast` 0.00 / 0.64 / 4.9 / 23.86 (relative bins; the noisy-blurred value moves by about 0.05 with the
  noise draw). With 8 within-object rank levels it was 0.00 / 10.30 / 2.26 / 5.96, the noisy flat disc above the
  checkerboard, which is why rank levels were rejected;
- `h_chr_grad` 0.032 / 0.064 / 0.113 / 0.359;
- `h_chr_lap` 0.00 / 1.07 / 1.03 / 7.72.

**Clipped-pixel test, mandatory.** Render a uniform RGB disc of (40, 1, 100) on a (230, 200, 220) background, then set
G = 0 on a random 7 % of the disc's pixels (canvas 71 × 71, radius 16 px, `M_chr` = the disc). Every texture column
must stay near its value on the unaltered disc, because the clipped pixels and their neighbourhoods are excluded.
Checked values, one implementation (`h_chr_grad` / `h_chr_lap` / `h_chr_glcm_contrast`): unaltered 0.036 / 0 / 0;
with the clipped pixels excluded 0.038 / 0 / 0; with them included 0.113 / 3.89 / 6.21. The tolerance, 10 %, was set
after seeing these values. The exclusion costs support: 709 px fell to 383, because scattered clipped pixels take
their 3 × 3 neighbourhoods with them. On real data clipped pixels cluster in the densest core, but read
`tex_support_px` beside the texture columns.

**Shuffle test, mandatory.** For 200 random accepted candidates, permute the H values inside `M_chr` and recompute.
`h_chr_cv` must be unchanged (it is a dispersion, and the test proves it). Every column that claims to be texture must
change by a median of at least 25 % in relative terms. Also test a smooth gradient against the same pixels shuffled:
the texture columns must differ, `h_chr_cv` must not.

### 5.5 Surround: membrane, clear zone, white space

| column | definition | unit |
|---|---|---|
| `rim_core_nuc` | mean H over the 2 px inner rim of `M_nuc` divided by mean H over its core (erode 2 iterations, 3×3 cross); NaN if either empty; the earlier F6 on the new mask | ratio |
| `ring_profile_slope` | least-squares slope of `(h_ring1_mean, h_ring2_mean, h_ring3_mean)` against distance | OD/µm |
| `ring_peak_prominence` | for k in {2, 3}: `h_ring_k_mean` minus the minimum of `h_ring_j_mean` over j < k; take the largest, divide by `h_chr_mean`, floor at 0. A rim that is darker than the band inside it scores above 0; a profile that only falls scores 0. It scores any band 1–3 µm outside `M_chr` that is darker than the band inside it: an envelope around a clump, but also a nearby partner plate or a cell edge; see "What the membrane columns can see" | ratio |
| `ring1_frac_nuclear` | share of `ring_1` pixels above `t_nuc` | share |
| `ring1_frac_glass` | share of `ring_1` pixels at or below `t_nuc` and with `gray > 215` (the `WHITE_GRAY` of the floor-clip run): glass and lumen | share |
| `ring1_frac_pale` | share of `ring_1` pixels at or below `t_nuc` and with `gray ≤ 215`: pale tissue, cleared cytoplasm, halos. The three `ring1_frac_*` shares partition `ring_1` and sum to 1 | share |
| `ring3_frac_nuclear` | share of `ring_3` pixels above `t_nuc` | share |
| `white_frac_gray215_61` | share of the 61 px window with `gray > 215`; the earlier `not_white` without the sign | share |
| `white_frac_od015_61` | share of the 61 px window whose per-channel OD sum, `−log10(R/255) − log10(G/255) − log10(B/255)` with each channel clipped to at least 1, is below 0.15 (about `gray > 227`; the `215` threshold is OD sum about 0.22) | share |
| `white_frac_od015_101` | the same on a 101 px window | share |
| `e_ring23_mean` | mean E (`hed[:, :, 1]`) over the cytoplasm ring: `ring_2` ∪ `ring_3` minus glass (`gray > 215`), minus pixels above `t_dark`, minus clipped pixels; NaN with `nan_reason` `no cytoplasm` below `NS_MIN_PX` | OD |
| `e_over_h_ring23` | mean E / mean H over the same cytoplasm ring | ratio |

**Blur test, mandatory.** Synthetic dark disc, Gaussian-blurred at σ = 1 and 2 px, no ring: `rim_core_nuc` will fall
below 1 from blur alone; record by how much. The same disc with a dark ring 2 px wide at radius `r + 1.5 µm`: record
the values. A feature that cannot tell the ring from the blur at σ = 2 is reported as unable to. Expect
`ring_peak_prominence` to pass and `rim_core_nuc` to fail; whichever way it comes out, say so.

**What the membrane columns can see (C2), checked 2026-10-02 with the verbatim refinement.** In an H&E interphase
nucleus the membrane is not a separate object that `M_chr` sits inside. Where it shows, it is a thin dark rim, and the
window's top multi-Otsu class takes the rim itself. Synthetic windows, `mpp` 0.25, blur σ 1 px:

| window | refinement mode | what `M_chr` becomes | `ring_peak_prominence` | `chr_hole_frac` |
|---|---|---|---|---|
| vesicular nucleus: 0.5 µm dark rim, pale interior | fallback | the rim | 0.06 | 0.78 |
| the same plus a nucleolus | fallback | rim + nucleolus | 0.00 | 0.80 |
| uniformly dark disc | base | the disc | 0.00 | 0.00 |
| 2 µm clump inside an intact rim | fallback | the clump | 0.15 | 0.00 |

- **`ring_peak_prominence`** scores highest on a dark envelope around an internal clump, but it is not specific to
  it. It gave 0.05 on an anaphase pair 4.5 µm apart with the detection on one plate (the partner plate falls in
  `ring_3`), and 0.03 on a metaphase plate whose clear zone is bounded by darker cytoplasm (Section 5.2's pair-test
  windows). No source describes prophase as a clump inside a dark rim: VCGP 2021 counts from prometaphase on, and its
  prometaphase is "A dark cluster or ring of chromosomes with variably protruding chromosome rods and spikes". Its
  direction is 0; the hypothesis that it marks an intact envelope (C2, C6) goes in the rationale only.
- **`chr_hole_frac`** (Section 5.2) is where a vesicular nucleus's membrane shows up, but a ring-shaped mitotic figure
  gives the same mask: VCGP 2021 counts "A dark cluster or ring of chromosomes with variably protruding chromosome rods
  and spikes" (prometaphase) and a "linear or ring-shaped aggregate" (metaphase) as mitotic figures. Its direction
  is 0.
- **`rim_core_nuc`** reads the rim against the core of `M_nuc`, so it can see a vesicular nucleus, but only where
  `M_nuc` is the nucleus. Where `nuc_clipped` is true (79 % of rows on the recorded A3 pool, seeds 0–1, measured
  2026-10-02) part of the rim is the edge of the 2.5 µm disc, in cytoplasm or a neighbour. Its direction is 0, because
  `directions.csv` holds one direction per column; the hypothesis "lower where `nuc_clipped` is false, because a
  darker rim is a membrane" goes in the rationale, and `auc_split_nuc_clipped` tests it.
- **On a uniformly dark nucleus no membrane column can work.** Pathologists separate those by their smooth outline
  (C5), which is the job of Sections 5.2 and 5.8 (b). A membrane column scoring them like mitoses is expected, not a
  defect.
- **Membrane test, mandatory, written first.** Reproduce the four rows above within ±0.03 before real data.
  Construction, on the canvas of Section 5.2's pair test (201 × 201, 0.03 OD plus Gaussian noise of σ 0.003, any seed;
  shapes drawn over it with constant values that replace the noise; then blur σ 1 px and cast to float32), detection
  at the centre. With r the distance from a disc's centre, "the outer 2 px" of a disc of radius R means the pixels with
  R − 2 < r ≤ R; a distance-transform reading gives a thicker rim and moves the clump value to about 0.12.
  - vesicular: a disc of radius 24 px at 0.07 OD whose outer 2 px are 0.16 OD;
  - nucleolus: add a disc of radius 4 px at 0.16 OD, centred 3 px right and 2 px up of the centre;
  - dark disc: radius 18 px at 0.17 OD;
  - clump: a disc of radius 20 px at 0.07 OD whose outer 2 px are 0.15 OD, with a central disc of radius 8 px at
    0.20 OD.
  - Two independent implementations agree on all four rows under this reading (clump 0.148–0.150, vesicular
    0.058–0.061 / 0.78, across noise seeds; adding the noise after drawing moves them to about 0.145 and 0.066 / 0.79).

**Clear-zone logic, corrected 2026-10-02 (C9).** The earlier text said a mitosis sits in cleared cytoplasm, so its
`ring1_frac_nuclear` should be low. That has no verified source.
- The H&E statement found is "The cytoplasm of the mitotic cells was often larger during mitosis than in the resting
  cells" (Ikpatt 2002). Rounding at mitosis is documented in cell biology, not as a counting criterion.
- A clear halo is documented around apoptotic cells (rodent liver, NTP atlas), and this run's own false positives
  include smooth dark nuclei with halos (Section 1).
- "Below `t_nuc`" is the clip floor's class (Section 4), so cleared cytoplasm, halos, lumen and glass fall together.
  `ring1_frac_glass` and `ring1_frac_pale` split glass and lumen from pale tissue.
- So `ring1_frac_nuclear`, `ring1_frac_glass`, `ring1_frac_pale` and `ring3_frac_nuclear` have direction 0. The one
  contrast they are expected to carry: a dark clump inside an interphase nucleus is surrounded by nucleoplasm above
  `t_nuc`, so its `ring1_frac_nuclear` should be high.

### 5.6 Colour and pigment

RGB statistics are on the raw `uint8` image. `hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)` once per ROI.

| column | definition | unit |
|---|---|---|
| `r_chr_mean`, `g_chr_mean`, `b_chr_mean` | RGB means over `M_chr_clip` | 0 to 255 |
| `blue_index_chr` | `(b_chr_mean − r_chr_mean) / 255` | ratio |
| `sat_chr_mean` | HSV saturation mean over `M_chr_clip` | 0 to 255 |
| `hue_chr_circmean` | circular mean of HSV hue over `M_chr_clip`. OpenCV's `uint8` hue is in half-degrees (0 to 179), so double it first | degrees |
| `e_over_h_chr` | mean E / mean H over `M_chr_clip` | ratio |
| `spec_e_over_h_chr` | mean E / mean H over all of `M_chr`: the earlier F8, unchanged, kept for gate 7 | ratio |
| `spec_d_over_h` | mean D / mean H over `M_chr`; the earlier F9, kept under its spec name for comparison | ratio |
| `d_mean_surround`, `d_p95_surround` | D over the 51 px window outside `M_sel` | OD |
| `brown_frac_surround` | share of window pixels outside `M_sel` with `R > G > B`, `R − B > 40` and `gray < 180` | share |
| `brown_frac_101` | the same on a 101 px window | share |
| `d_over_h_surround` | mean D / mean H over the window outside `M_sel` | ratio |

The chromatin colour columns (the first five rows) read `M_chr_clip` (Section 4): `M_chr` minus truly clipped pixels
(a channel at 0 or 255), NaN with `nan_reason` `clipped` below `NS_MIN_PX`. They do not drop pixels with a channel
below 20: those are dense blue-purple hematoxylin, and dropping them removed 60 % of a mitosis's chromatin against 37 %
of other candidates' (Section 4). Only the two `spec_` columns read all of `M_chr`, for comparison and regression.

The brown rule is a starting heuristic. Calibrate it on the melanoma montages if the top tiles are not pigment, log the
change, and keep the original under `spec_brown_frac_surround`.

**Where pigment lands, corrected 2026-10-02.** The earlier rationale ("pigment lives outside the hematoxylin-selected
mask") is true for only one kind of pigment.
- **Light brown pigment stays outside `M_chr`.** For lightly deposited melanin, and the golden-brown granules
  (hemosiderin-like) seen on `410` and `291`, the OD vector is blue-heavy, and colour deconvolution assigns it mostly
  to D: for a brown OD direction ∝ (0.3, 0.5, 0.8) of unit norm, the H coefficient (`skimage.color.hed_from_rgb`) is
  about 0.05, against 1.0 for pure hematoxylin. These are deconvolution coefficients, not orthogonal projections.
  Melanin's exact vector was not measured. So light brown pigment stays largely outside `M_chr`, and the surround
  features above are the right tools for it.
- **Dense melanin is the exception, for two reasons.** Under Beer–Lambert its OD direction is fixed, but its H value
  grows linearly with density, at that 0.05 coefficient. And in a dense deposit the blue channel saturates first,
  which pushes the measured OD vector towards the neutral case below (coefficient about 0.70). So a dense deposit looks
  brown-black and can enter `M_chr`. `chr_resid_frac` scores it either way, as long as its rim keeps enough
  non-saturated pixels.
- **Neutral black pigment *is* the `M_chr` of its candidates.** A neutral (equal-OD) absorber of unit OD norm gets an
  H coefficient of about 0.70. Black pigment is also far more opaque than chromatin.
  - On `211`, the top-10 false positives are black granular clumps.
  - Their darkest pixels have D 0.152 against 0.034 for the mitoses: `rgb2hed` D averaged over the darkest 10 % (by
    H) of a 21 × 21 px core, as `colour_feats` in `low_precision_rois/roi_mechanism_measure.py` computes it.
  - The feature for this case is `chr_resid_frac` in Section 5.8. It reads `M_chr` plus a 2 px rim, because the clumps'
    cores are saturated.
- **Why F9 (`spec_d_over_h`) failed.** It rose with darkness on the darkest-stained ROIs and so penalised the darkest
  *mitoses* (Section 1). Corrected 2026-10-02: those mitoses' darkest pixels are not neutral. Their green channel sits
  below 20 counts (92 % of the 2,658 refined GT mitoses have such pixels, on average half their `M_chr`, measured
  2026-10-02 against `NS_LO`/`NS_HI`), where one count is a large step in OD, while red and blue stay above it on
  average. The exact mechanism of the D/H rise was not
  isolated. Any colour feature on the chromatin must exclude truly clipped pixels and be read beside
  `rho_with_h_chr_mean`.

**Pigment checks, mandatory.**
- **Brown surround features:** run them on melanoma `509`, `524`, `527`, `546` and on `410`, `291`.
  - The 2026-10-02 montages showed no pigment among the top-10 false positives of `524`, `527` and `546`, so a weak
    result there is the expected outcome, not a failure.
  - `509` has a brown-pigmented look-alike, and `410` and `291` have golden-brown granules near candidates.
- **The in-mask black-pigment feature:** run it on `211`, plus `289` and `339` (black or brown-black clusters seen
  near top-10 false positives).

Verify that:
- the ten highest `brown_frac_surround` tiles are visibly brown pigment on the ROIs where pigment was seen (montage, your
  eyes);
- the ten highest `chr_resid_frac` accepted tiles on `211` are black pigment clumps (montage, your eyes), and report its
  within-click AUC on `211`, `289` and `339`;
- `brown_frac_surround` on melanoma ROIs is higher than on non-melanoma ROIs at the ROI median;
- `spec_d_over_h` rises with `h_chr_mean` on the darkest-stained ROIs (`425`, `460`, `201`) and falls elsewhere,
  reproducing the earlier finding, and report whether `dark_frac_chr` explains it: the Spearman of `spec_d_over_h`
  with `dark_frac_chr` within click.

### 5.7 Context and density

| column | definition | unit |
|---|---|---|
| `n_dark_comp_101` | 4-connected components above the ROI's upper multi-Otsu threshold `t_dark` in a 101 px window, each at least 1 µm² | count |
| `dark_area_frac_101` | share of the 101 px window above `t_dark` | share |
| `nuc_area_frac_101` | share above `t_nuc` | share |

`t_dark` and `t_nuc` are the two ROI-wide multi-Otsu thresholds; store both per ROI in `clicks.csv`. They are the
floor-clip run's `b_otsu_upper` and `b_otsu`, from the same call (Section 4). Assert
`float(np.float32(t_dark)) == b_otsu_upper` from `floors.csv` on every ROI.

### 5.8 Added 2026-10-02: in-mask pigment, thin protrusions, edge sharpness, size relative to neighbours

These four families target the dark false positives of the failing ROIs (Section 1). Every rule at the top of Section 5
applies.
- Every constant named here goes into the spec JSON.
- `H_w` is the unclipped 51 px `hem` window that the masks live in.
- `M_chr` is as in Section 4.
- `rgb_w` is the same 51 px window of the raw `uint8` RGB image. With the band out it always lies inside the ROI.

**(a) `chr_resid_frac`, `chr_ns_share`: colour off the hematoxylin–eosin plane, on the chromatin and its 2 px rim.**
- **Per-pixel optical density:** `od = -np.log10(np.maximum(rgb_w.astype(np.float64), 1.0) / 255.0)`, with channels in
  R, G, B order.
- **The region:** `M_res` = `scipy.ndimage.binary_dilation(M_chr, structure=skimage.morphology.disk(RESID_DILATE_PX))`,
  with `RESID_DILATE_PX` = 2. `M_chr` is *whole*; the region is clipped to the window.
- **Non-saturated pixels:** `M_res_ns` = `M_res` ∧ all three channels in [`NS_LO`, `NS_HI`] = [20, 235].
  - Store `chr_ns_share` = |`M_res_ns`| / |`M_res`|.
  - If |`M_res_ns`| < `NS_MIN_PX` = 10, then `chr_resid_frac` is NaN with `nan_reason` `outside_20_235`.
- **Why the rim is included (measured 2026-10-02).** On `M_chr` alone, the black clumps on `211` are mostly saturated.
  - 44 % of `211`'s A3 top-10 false positives had fewer than 10 non-saturated pixels, so the feature was NaN on its
    own target.
  - 35 % of `289`'s mitoses were NaN too.
  - The clumps' partly transmitting 2 px edge carries their colour. With the rim included, nothing was NaN on the
    seven ROIs checked.
- **The normal to the plane:** `h`, `e` = `skimage.color.rgb_from_hed[0]` and `[1]`, each divided by its norm. Then
  `n = np.cross(h, e)`, normalised, with its sign chosen so that `n[2] > 0`. That gives about (−0.332, −0.081, 0.940).
- **The feature:** `v` = the mean of `od` over `M_res_ns`, and `chr_resid_frac = (v @ n) / np.linalg.norm(v)`. It is
  signed and dimensionless.
- **What it measures.** Any mixture of the two stains, as `rgb2hed` defines them, lies in the plane spanned by `h` and
  `e` and scores 0. A neutral (black or grey) absorber scores 0.304, and a brown OD vector ∝ (0.3, 0.5, 0.8) scores
  0.618.
- **Unit test, written first.** Use uniform 21 × 21 patches rendered as `I = round(255 · 10**(−OD))` per channel, with
  the OD scaled so that every channel lies in [20, 235], and `M_chr` = the central 11 × 11 pixels, so that `M_res` stays
  inside the uniform patch. Expected values:
  - pure `h` → 0.000 ± 0.02;
  - 0.7·`h` + 0.3·`e` → 0.000 ± 0.02;
  - neutral → 0.304 ± 0.02;
  - brown → 0.618 ± 0.02;
  - a patch with every channel below 20 → NaN.
- **Why the fixed `rgb2hed` vectors are acceptable here.** Each scanner's real stain vectors differ from them:
  per-ROI Macenko H vectors are 5–15° away at the scanner median (`floor_robustness/macenko.csv`). That shifts every
  candidate within an ROI alike, so the within-click comparisons in Section 6 are unaffected. Pooled values are **not**
  comparable across scanners; say so wherever you report them.
  - The mitoses' median baseline was about 0.18–0.21 on two Aperio ROIs (`300`, `339`), about 0.08 on `211`, and about
    0.01 on `527` (Hamamatsu XR).
  - The OD formula takes 255 as white. A slide whose glass is darker adds a small neutral OD to every pixel, about
    0.02 at chromatin level against about 0.3 for pigment. That is a mild darkness confound, and
    `auc_dark_stratified` would expose it.
- **Real-data anchor (2026-10-02).** This compares A3 top-10 false positives with GT mitoses refined at the annotation
  centres. It is a sample, not this spec's candidate population.
  - The AUC that a false positive scores higher than a mitosis is 0.83 on `211` and 0.86 on `289`, both ROIs with
    pigment.
  - On ROIs without pigment it is 0.20–0.49 (`300`, `339`, `425`, `460`, `527`). There the *mitoses* score slightly
    higher.
  - So expect a pigment-specific signal at the top of the range, not a monotone one. Do not tune `RESID_DILATE_PX` or
    `NS_LO` to change this.
- **Direction hypothesis:** −1 (pigment is higher).

**(b) `chr_tophat_edge`, `chr_tophat_core`: the share of the object's OD carried by structures thinner than about
1 µm, at the surface and inside.**
- **The thin-structure map.**
  - `r_th = max(1, round(TOPHAT_R_UM / mpp))` with `TOPHAT_R_UM` = 0.5, which gives 2 px at 0.23–0.25 µm/px.
  - `se = skimage.morphology.disk(r_th).astype(np.uint8)`.
  - `T = H_w − cv2.morphologyEx(H_w.astype(np.float32), cv2.MORPH_OPEN, se, borderType=cv2.BORDER_REPLICATE)`. This is a
    grey-level white top-hat: the high-OD structures narrower than the disc.
  - The disc is 5 px, about 1.25 µm, across. A single chromatid is "0.68 ± 0.04 μm" wide and a cohesed chromosome
    "1.23 ± 0.1 μm" (Cisneros-Soberanis 2024, by EM), so the threshold sits at the physical width of one chromosome
    arm. Whether a scan resolves an arm that thin is a separate question (5.8 c, and Section 6's per-scanner split).
- **Regions, following the mask** (`M_chr` *whole*):
  - `R` = the window pixels with `distance_transform_edt(~M_chr) * mpp <= TOPHAT_RING_UM`, with `TOPHAT_RING_UM` = 1.0,
    so that fuzzy protrusions up to 1 µm outside the mask count. `M_chr` itself is included.
  - `M_core` = the pixels with `distance_transform_edt(M_chr) * mpp > TOPHAT_CORE_UM`, with `TOPHAT_CORE_UM` = 0.5:
    chromatin more than 0.5 µm inside its own edge.
  - `R_edge` = `R` ∧ ¬`M_core`.
- **The features:**
  - `chr_tophat_edge = T[R_edge].sum() / H_w[M_chr].sum()`: thin structure at the surface, the C1 cue;
  - `chr_tophat_core = T[M_core].sum() / H_w[M_chr].sum()`: thin structure inside, which granular interphase
    chromatin also has;
  - `spec_chr_tophat_frac = T[R_disc].sum() / H_w[M_chr].sum()`, the earlier definition, kept for comparison. `R_disc`
    is the disc of radius `r_eq + round(TOPHAT_RING_UM / mpp)` px about the centroid of `M_chr`, intersected with the
    window, with `r_eq = sqrt(|M_chr| / π)` px.
  - All three are NaN when `M_chr` has fewer than 30 px, as in Section 5.2.
  - Store `tophat_touches_window`, true when `R` reaches the window edge.
- **Why the region changed (measured on the reference construction below).**
  - A centroid disc follows the mask's size, not its shape. On an about 4:1 bar at r = 16 with three 2 px spikes at
    each tip, `spec_chr_tophat_frac` is 0.000 because the tips lie outside the disc; `chr_tophat_edge` is 0.146.
  - The earlier column mixed opposite cues. A disc with 2 px internal granules scored 0.092 at r = 16, above the 2 px
    spikes (0.086). The split puts that in `chr_tophat_core` (0.092), not `chr_tophat_edge` (0.000).
  - The new region depends on the mask's outline, but only through a 1 µm dilation and a 0.5 µm erosion, which smooth
    over a ragged edge.
- **Synthetic suite, expected ordering written first.** Render each shape as an H image, 0.20 OD on a 0.03 background,
  at radius 8 and 16 px:
  - a disc scores low on both columns;
  - the 12-spike shape (the reference construction at the end of this section, `w` = 2) scores clearly above the disc on
    `chr_tophat_edge` at both scales;
  - the same shape with `w` = 6 scores far below the `w` = 2 version on `chr_tophat_edge`;
  - the `w` = 2 shape blurred with Gaussian σ = 1 px: report the drop. Blur erasing the cue is the resolvability
    question Section 6 asks per scanner;
  - the granular disc scores clearly above the disc on `chr_tophat_core`, and no more than 0.01 above it on
    `chr_tophat_edge`;
  - the tip-spiked bar scores clearly above the plain bar on `chr_tophat_edge` at both scales;
  - **a known confound, not a pass:** the vesicular window of Section 5.5, refined verbatim, scores `chr_tophat_edge`
    0.34, above the 2 px spikes (0.256), because its `M_chr` is a thin continuous rim. A top-hat cannot tell a smooth
    thin membrane from discrete spikes. Record the value; `auc_split_hole` (Section 6) reads the effect on real
    data.
- **The scale rule does not apply.** These columns are exempt from Section 5.2's scale rule, because their width
  threshold is physical by design. At r = 8 the plain bar is 5 px wide, exactly the structuring element's width, so the
  opening removes only its tip corners (12 of its 85 px); a bar thinner than 5 px would be removed entirely.
- **Small bodies score high too.** The whole of a body narrower than about 1.25 µm counts as thin structure. So
  karyorrhectic fragments and other tiny round bodies, which are imposters (C5, C7), raise `chr_tophat_edge` as
  spikes do. Read the column with `chr_n_pieces` and `auc_area_tertiles`, and write this confound
  into `directions.csv`.
- **Direction hypotheses:** `chr_tophat_edge` +1, with the small-body confound above; `chr_tophat_core` 0
  (Section 5.0, texture); `spec_chr_tophat_frac` +1, as before.

**(c) `chr_edge_sharpness_per_um`: how crisp the chromatin boundary is.**
- **The boundary `B`:** the largest 4-connected piece of `M_chr`, minus its `binary_erosion` with a 3 × 3 cross.
- **The gradient:**
  - `gx = cv2.Sobel(H_w, cv2.CV_32F, 1, 0, ksize=3) / 8` and `gy` likewise;
  - `g = np.hypot(gx, gy) / mpp`, in OD per µm.
- **The step:** `step = h_chr_mean − h_ring1_mean`, from Section 5.3. If `step` < `EDGE_MIN_STEP` = 0.005 OD, the
  feature is NaN with `nan_reason` `no edge step`.
- **The feature:** `chr_edge_sharpness_per_um = g[B].mean() / step`, in 1/µm. It is roughly the inverse width of the
  edge.
- **Synthetic suite.** A disc and a 2:1 ellipse, at radius 8 and 16 px, 0.20 OD on 0.03, each Gaussian-blurred at
  σ = 0, 1 and 2 px:
  - the value must fall strictly with σ;
  - the ellipse must score within 15 % of the disc at the same σ;
  - the two radii must score within 15 % of each other.
- **What it is for.**
  - Blurred dark objects are a named false-positive class.
  - The value is also the covariate for whether boundary features are resolvable on a given scan. In MIDOG 2022,
    27.5 % of the mitoses found only with PHH3 help were classed "out of focus" (Aubreville et al. 2024, Table 7).
- **Direction hypothesis:** 0. Blurred objects should score low, smooth pyknotic bodies high, and ragged mitotic edges
  in between, so no monotone direction is claimed.

**(d) `nbr_n`, `nbr_nuc_area_med_um2`, `chr_area_over_nbr`: chromatin size relative to the neighbouring nuclei.**
- **The window:** the same 101 px window as `h_win101_mean`, centred on the detection and truncated at the ROI edge.
- **The components:** 4-connected components of `hem_101 > t_nuc`, with `t_nuc` from Section 4.
- **Drop:**
  - every component touching the (truncated) 101 px window's edge;
  - every component intersecting the candidate's own `M_chr`, placed at row and column offset
    `(cy − 25 − y0, cx − 25 − x0)` inside the 101 px window, where `(y0, x0)` is the truncated window's origin in ROI
    pixels (25 only when the window is not truncated);
  - every component with an area outside [`NBR_MIN_UM2`, `NBR_MAX_UM2`] = [8, 150] µm². This drops fragments, and
    merged clusters or tissue sheets.
- **The outputs:**
  - `nbr_n` is the number of components kept;
  - `nbr_nuc_area_med_um2` is their median area, NaN with `nan_reason` `few neighbours` when `nbr_n` < 3;
  - `chr_area_over_nbr = chr_area_um2 / nbr_nuc_area_med_um2`.
- **Calibrating the bounds.** The bounds are starting values. Check them on a lymphosarcoma montage (dense nuclei merge
  under `t_nuc`) and a sarcoma montage. Log any change, and keep the original under a `spec_` column.
- **Unit test.** A 101 × 101 field holding:
  - 12 non-touching discs of radius 12 px;
  - one disc at the centre, standing in for the candidate's `M_chr`;
  - two discs cut by the window edge.

  It must give `nbr_n` = 12 and a median equal to the rendered disc's pixel count × mpp², within 2 %.
- **Direction hypothesis:** 0. It should be +1 against apoptotic bodies and lymphocytes, which are smaller than their
  neighbours, but −1 against hyperchromatic interphase nuclei.

**Reference shapes and checked values (computed 2026-10-02 from the definitions above).** Reproduce these before
running on real data. A different value means your implementation differs from this text; find out why before
continuing.

**Construction:**
- `mpp` = 0.25.
- The canvas is 51 px for r = 8 and 71 px for r = 16.
- Shapes are centred, at 0.20 OD on a 0.03 background, as float32.
- The disc is `x² + y² ≤ r²`.
- The ellipse has semi-axes r·√2 and r/√2.
- **The 12-spike shape** is the disc plus 12 rectangles at angles 2πk/12. Each rectangle covers `r − 1 ≤ t ≤ r + 4`
  along the ray and `|o| ≤ w/2` across it, with spike width `w` = 2. The "6 px lobes" variant uses `w` = 6.
- **The granular disc** is the disc plus 0.06 OD on the pixels more than 3 px inside its edge
  (`distance_transform_edt(disc) > 3`) where `(floor(x / 2) + floor(y / 2)) % 2 == 0`, with x and y the centred
  coordinates: a 2 px checkerboard. `M_chr` is the disc.
- **The bar** is `|x| ≤ r/4` and `|y| ≤ r`, about 4:1. **The tip-spiked bar** adds six rectangles built like the
  12-spike shape's (`w` = 2, `r − 1 ≤ t ≤ r + 4`), at angles π/2 + {−0.5, 0, +0.5} rad and −π/2 + {−0.5, 0, +0.5} rad,
  where `t = x·cos(a) + y·sin(a)`.
- Blur is `scipy.ndimage.gaussian_filter` with default settings.
- Ring 1 uses `scipy.ndimage.distance_transform_edt`.
- In every synthetic test, `M_chr` is the shape's own *unblurred* pixel set, even when the image is blurred.
- The whole canvas is the window `H_w`.
- The neighbour unit test in (d) also uses `mpp` = 0.25. A radius-12 disc is about 28 µm², inside the bounds.

**Checked values:**
- **`chr_resid_frac` unit test:** 0.001, 0.000, 0.304 and 0.615 for pure `h`, the mix, neutral and brown.
  - The rendered pixels are (77, 70, 149), (106, 60, 165), (64, 64, 64) and (137, 90, 49).
  - The ODs are 0.8·`h`, 0.8·(0.7`h` + 0.3`e`), 0.6 per channel, and 0.9·(0.3, 0.5, 0.8).
- **The top-hat columns,** r = 8 / r = 16. The `spec_chr_tophat_frac` values of the first four rows are the values
  this spec already gave for `chr_tophat_frac`, except the 6 px lobes at r = 8: both independent implementations
  compute 0.0145, which rounds to 0.014, not the 0.015 given earlier. Two independent implementations reproduce every
  value in this table.

  | shape | `spec_chr_tophat_frac` | `chr_tophat_edge` | `chr_tophat_core` |
  |---|---|---|---|
  | disc | 0.035 / 0.000 | 0.035 / 0.000 | 0.000 / 0.000 |
  | 2 px spikes | 0.256 / 0.086 | 0.256 / 0.086 | 0.000 / 0.000 |
  | 6 px lobes | 0.014 / 0.030 | 0.014 / 0.030 | 0.000 / 0.000 |
  | 2 px spikes, σ = 1 | 0.118 / 0.041 | 0.118 / 0.043 | 0.001 / 0.000 |
  | granular disc | 0.089 / 0.092 | 0.033 / 0.000 | 0.056 / 0.092 |
  | bar | 0.120 / 0.000 | 0.120 / 0.034 | 0.000 / 0.000 |
  | tip-spiked bar | 0.193 / 0.000 | 0.313 / 0.146 | 0.000 / 0.000 |

- **`chr_edge_sharpness_per_um` for the disc at σ = 0 / 1 / 2:** 1.93 / 1.53 / 1.16 at r = 8, and 1.97 / 1.47 / 1.04 at
  r = 16. The ellipse is within 6 % of the disc.

**Flags for this section** (Section 7):
- by bucket: the share of accepted rows with `chr_resid_frac` NaN (`outside_20_235`), and the median `chr_ns_share`;
- the share with `tophat_touches_window`;
- by domain: the share with `nbr_n` < 3.

## 6. The validity report (descriptive; no selection, no verdict)

One row per feature in `validity.csv`, computed on accepted rows of the recorded A3 pool (`in_recorded_pool`), the
pool production ranks:

- `auc_within_click`: the mean over clicks of the per-click AUC that a mitosis has the larger value, over clicks with
  at least one mitosis and one non-mitosis; also `auc_within_click_pooled_u`, the ratio of summed U to summed
  `n1·n0` over clicks. Both unsigned, raw direction.
- `auc_pooled`: the pooled AUC over all rows, for contrast with the within-click value. Where the two differ in sign,
  say so; that is the F9 lesson.
- `rois_above_half`: ROIs whose per-ROI AUC is above 0.5, out of 49, and `domains_above_half` out of 7.
- `auc_split_nuc_clipped`, `auc_split_nuc_touches_window`, `auc_split_win101_full` (`win101_area_share == 1`):
  within-click AUC on the two halves of each flag. A nucleus-based feature whose overall AUC vanishes inside both halves is measuring the flag.
- `auc_split_hole` (`chr_hole_frac > 0.3`), for the boundary columns of Sections 5.2 and 5.8 (b, c): within-click
  AUC on the two halves. A ring-shaped `M_chr` (a vesicular nucleus's rim) raises `chr_tophat_edge`.
- `auc_split_chr_clip` (`chr_clip_share > 0.05`), for the texture (5.4) and chromatin colour (5.6) columns.
- `auc_area_tertiles`, for the boundary columns of Sections 5.2 and 5.8 (b, c) only: computed like
  `auc_dark_stratified` below, with three within-click strata of `chr_area_um2` in place of five of `od`. Spikes on a
  small object approach the pixel pitch (Section 5.2, scale), so a boundary column that separates only between size
  strata is measuring size.
- `rho_with_h_chr_mean`, `rho_with_chr_area`, `rho_with_od`: mean within-click Spearman with darkness, size and the
  pipeline's `od`. This exposes duplicates and proxies.
- `nan_rate_mitosis`, `nan_rate_other`.
- `synthetic_pass`: pass / fail / n.a. from Sections 5.2, 5.4, 5.5 and 5.8.
- `direction_hypothesis` copied from `directions.csv`, and `direction_measured` as the sign of
  `auc_within_click − 0.5`. Report agreement; do not change the hypothesis after the fact.
- **`auc_dark_stratified`: the AUC beyond darkness. The key column of the 2026-10-02 amendment.**
  - **Strata.** Within each click, split its evaluated rows into five strata by within-click rank of `od` (the
    pipeline's darkness): `pd.qcut(rank(od, method='first'), 5, labels=False)`. Rows with NaN `od` are left out, and
    their count is reported.
  - **Combining them.** In every stratum with at least one mitosis and one non-mitosis, among rows where the feature
    is not NaN, compute the Mann–Whitney U that a mitosis has the larger value (ties 0.5).
    `auc_dark_stratified = ΣU / Σ(n1·n0)` over all strata of all clicks.
  - Also store `n_strata_used` and `n1_dark_stratified`.
  - **Reading it.** A value near 0.5 means the feature adds nothing once darkness is held roughly fixed, whatever its
    `auc_within_click`. Darkness features themselves will sit nearer 0.5 here by construction. That is the point of the
    column, not a defect.
- **`auc_wc_3dhistech`, `auc_wc_aperio_cs2`, `auc_wc_hamamatsu_xr`, `auc_wc_hamamatsu_s360`:** `auc_within_click`
  restricted to each scanner's clicks, each with its click count.
  - Scanner is confounded with tumour type for 42 of the 49 ROIs, and S360 is 3 ROIs. Read these as scanner-and-tissue
    groups.
  - For the boundary features in Sections 5.2, 5.5 and 5.8 (b, c), this is the "is the cue resolvable on this scan"
    check. Report it beside each group's median `chr_edge_sharpness_per_um`.
- **`auc_wc_failing9`, `auc_wc_other40`:** `auc_within_click` by `roi_group`, with each group's total mitosis count.
  This answers whether a feature carries signal where darkness fails. The grouping is post hoc, so it is descriptive.

**A second table, `validity_depth300.csv`:** the same columns, computed on all accepted rows with `cand_rank` < 300. The
primary table stays on the pool production ranks; the depth-300 table is where the failing ROIs have enough mitoses to
read. Print, for each table, the mitosis count per `roi_group`.

Also write `redundancy.csv`: the full matrix of mean within-click Spearman correlations between all features, and
print the pairs above 0.8 in absolute value.

**What you may conclude from this report:** whether each measurement is valid, and what it actually measures. **What
you may not conclude:** which features to use, how to weight them, or whether any ranker is better. Write the report
in that voice.

## 7. Gates (halt and report) and flags (report only)

**Gates:**
1. **Config.**
   - `CHANNEL`, `TM_METHOD`, `PEAK_MIN_DISTANCE`, `DEEP_FLOOR_Z`, `MAX_PEAKS` and `OD_WINDOW` equal
     `midog_utils.production`.
   - `ev.MIDOG_RADIUS_UM == 7.5`, `PAD_PX == 2`, and no `PAD_UM` name exists.
   - The refinement's source SHA-1 is `21c43a20db00dccf7ac3bab88aee423c345441b6`.
   - The A3 constants hold: `TPL_ARM == 'mo_pad61'`, `S_TPL == 61`, `PAD == 30`, `BAND_MODE == 'floored'`,
     `BLANK_SIZE == 51` and `FLOOR_CLASSES == 3`.
   - The `difflib` check of Section 3 passes.
2. **Tripwire.** On all 294 recorded clicks:
   - the verbatim refinement gives the recorded `bs`, `tcx` and `tcy`;
   - the A3 template (from the clipped channel) has the recorded A3 `tpl_sha1` (`tm_floor_definition` per-run,
     `arm == 'A3'`);
   - `roi_floors` reproduces `b_otsu`.

   Print the count that differ, and halt above 0.
3. **Prefix identity.** Each click's first `n_detections_recorded` deep A3 rows equal the recorded A3 pool.
   - Compare `cx`, `cy`, `score` (float32), `od` (`rtol 1e-12`), `bucket` and `matched_ann_id`.
   - Also `threshold` and `n_detections_recorded` against the recorded `threshold` and `n_detections`.
   - Halt on any difference.
4. **Channel parity.** `rgb2hed` H equals `hem` bit for bit on every ROI.
5. **`od` parity.** Section 5.1.
6. **Label consistency.**
   - All 300 rows are matched in rank order with `ev.bucket_detections`, the click's own annotation dropped. Greedy
     matching in rank order makes the first `n_detections_recorded` rows' labels those of a prefix-only match, and
     gate 3 checks them against the recorded run.
   - Assert this directly on 5 clicks: match the prefix alone and compare.
   - (On the unclipped deep list the same property held on all 294 clicks, checked 2026-10-02.)
7. **Regression on unchanged features.** Features depend only on the location and the unclipped channel. So on every
   accepted row whose `(file_name, cx, cy)` exactly equals a location in
   `results/candrefine_rerank_mo_pad61_select_pool.csv`, these columns must equal that run's values exactly:
   - `h_chr_mean` = `core_od`;
   - `sel_solidity` = `solidity`;
   - `white_frac_gray215_61` = `−not_white`;
   - `h_chr_cv` = `texture`;
   - `spec_e_over_h_chr` = `e_over_h`;
   - `spec_d_over_h` = `d_over_h`;
   - `spec_spiky` = `spiky`;
   - the acceptance flag and the refusal class.

   About 57 % of the recorded A3 pool's distinct locations, roughly 15,700, are shared (measured 2026-10-02); print the
   count. A mismatch means the refinement or a mask changed; stop.
8. **Unit tests and synthetic suites** all pass, or every failure is reported as a failure in `synthetic_tests.csv`
   and in the log. A failing synthetic test is not a halt; hiding it is.

**Flags (report):** refusal counts by class and bucket; `nuc_clipped` and `nuc_touches_window` shares by bucket; the
ROI-level share above `t_nuc`; milliseconds per candidate; the number of distinct locations; any ROI where more than
5 % of accepted candidates have `chr_n_pieces > 1`; per click, `n_peaks` and `n_detections_deep`, plus the search time;
the minimum `win101_area_share`; the Section 5.8 flags; by bucket, the share of accepted rows whose chromatin colour
columns are NaN (`clipped`) and the median `chr_clip_share`, and the share of `pair_*` NaN by reason (`no anchor`,
`no partner`, `round piece`); the share with `pair_touches_window`; by bucket, the median `tex_support_px` and the
share of texture rows NaN (`clipped`); the share of `e_ring23_mean` NaN (`no cytoplasm`).

## 8. Outputs and files

Notebook: `production_hematoxylin_only/multiotsu_pad_arms/candidate_features_a3_top300.ipynb`.

Tables, in `results/`, prefix `candidate_features_a3_top300_`:
- `candidates.csv`: 88,200 rows, one per (click, candidate), every column of Section 5 plus Section 5.1.
- `locations.csv`: one row per distinct `(file_name, cx, cy)` with the location-only columns, so that a reader can see
  the features were computed once.
- `clicks.csv`: 294 rows. Per click:
  - the click and its template (`bs`, `tcx`, `tcy`, the A3 `tpl_sha1`);
  - the floor `b`;
  - the search's `threshold`, `n_peaks`, `n_detections_deep`, `n_detections_recorded` and search seconds;
  - `t_nuc`, `t_dark` and the ROI's share above each;
  - `mpp` and `radius_px`.
- `directions.csv`: one row per feature: `column`, `claims_to_measure`, `direction_hypothesis` (+1 larger is more
  mitosis-like, −1, or 0 for no hypothesis), `rationale`, `confound`, `test`. Where Section 5.0 or the section that
  defines a column gives a direction and a rationale, use them. It is written by a cell that runs
  before any feature is computed on real data, from a literal in the notebook; its SHA-1 goes into the spec JSON and
  is printed by both the smoke run and the full run, which must agree. Edit it after seeing real-data AUCs and the
  hypotheses are no longer hypotheses; if you must, log it.
- `validity.csv`, `validity_depth300.csv`, `redundancy.csv`, `synthetic_tests.csv` (Sections 6, 5.2, 5.4, 5.5, 5.8).
- `adversarial_log.csv` (Section 2).

Also:
- **The spec JSON:** `production_hematoxylin_only/multiotsu_pad_arms/data/candidate_features_spec.json`.
  - It holds every constant (including Section 5.8's), every feature function's source SHA-1 and every input's
    SHA-256 prefix, with sorted keys and no timestamps, so a rerun reproduces its SHA-1.
  - The `data/` subfolder is where this folder's JSON and CSV outputs live, since `fee67b8`.
- **Montages, embedded in the notebook only.** Do not write image files into the repo; figures live in notebook outputs
  (`plt.show()`), by the user's 2026-10-02 rule.
  - For your own visual checks (Section 2, step 9), save temporary copies under `$TMPDIR` and open those with the Read
    tool.
  - **Content.** Per feature, the ten highest and ten lowest distinct accepted locations at least 7.5 µm apart, framed
    by label, with `M_chr` and `M_nuc` outlined.
  - **ROIs `300`, `201`, `013`, `529` and `425`:** every feature.
  - **ROIs `211`, `546` and `524`** (three failing ROIs: black pigment; smooth dark round nuclei against spiky figures;
    crowded hyperchromatic nuclei): only the Section 5.2 boundary columns, Section 5.5, and Section 5.8.
  - Your one-sentence judgement per feature per montage goes in the notebook's final markdown cell.

The final markdown cell holds, in this order: the per-feature validity verdicts in your own words; the montage
judgements; the adversarial log; what in this spec you changed and why; what you could not check.

Smoke runs write `_smoke`-suffixed copies of everything. The full run refuses to start if any output exists. Every
output is written at the end.

## 9. Environment, conventions, order of work

- Interpreter: `/Users/mohinianand/anaconda3/bin/python3`, the only one that imports `cv2` and `midog_utils`.
- Execute with `/Users/mohinianand/anaconda3/bin/jupyter nbconvert --to notebook --execute --inplace
  --ExecutePreprocessor.timeout=14400 <path>`. The notebook sits two levels below the repo root:
  `sys.path.insert(0, '../..')`, images at `../../images/...`, results at `../../results/`.
- Expected runtime:
  - the deep A3 search, about 45 minutes;
  - about 84k distinct refinements plus features at roughly 10 to 20 ms each, so 15 to 30 minutes;
  - the tripwire, about 2.5 minutes. Measure milliseconds per location in the smoke run, and print the projection
  before the full run.
- Docstrings: the opening `"""` on its own line, a brief summary, one line per parameter with type and a short
  description, one line for the return. Function definitions and calls on a single line, never wrapped.
- Do not commit. Do not modify or delete any existing file. Kill only a PID you started.
- Order of work: Section 2 steps 1 to 3, then unit tests and synthetic suites, then the smoke run, then the
  self-review, then the full run, then Section 2 steps 8 to 10, then the final markdown cell.
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
