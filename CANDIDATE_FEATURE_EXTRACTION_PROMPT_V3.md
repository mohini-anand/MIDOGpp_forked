# Per-candidate feature extraction on the adopted A3 search's 100-peak pool (`mo_pad61` + floor clip + border band): agent prompt

**Status:** spec, **version 3**, written 2026-10-09. Nothing in it has been run.
- **What V3 is:** V2 (`CANDIDATE_FEATURE_EXTRACTION_PROMPT_V2.md`, kept unchanged) with two things applied:
  - the fix plan `CANDIDATE_FEATURE_EXTRACTION_FIX_PLAN.md`;
  - the user's two decisions of 2026-10-09.
- **Where the evidence is:** the fix plan holds the evidence behind every change. That is its verification of the V2 run (`d99e93a`) and of that run's audit (`a4a6c0f`).
- **What it replaces:** the feature part of `CANDIDATE_REFINEMENT_RERANK_PLAN.md` (Section 5), as before. It does **not** replace that plan's ranking, selection or verdict parts.

**This run is small (Section 0.1).** It runs no search and extracts no features. It:
- builds this run's tables from the V2 run's committed rows;
- adds one column;
- recomputes the validity report by the corrected method;
- draws one montage;
- writes the findings, in a new notebook.

Sections 3, 3.1, 4 and 5 still define the stage in full, for a future full run. Their first lines say so.

**The V3 revision (2026-10-09).**

| what changed | where |
|---|---|
| The candidate set is production's 100-peak pool: the 100 highest peaks, then NMS on those only (79–100 per click). Depth 300, `detect_deep` and every beyond-pool reading are removed (user) | Sections 0, 3, 3.1, 6 |
| The neighbour columns are dropped: `chr_area_over_nbr`, `nbr_n`, `nbr_nuc_area_med_um2`, and the V2 run's three `spec_` versions (user). Result: 24 features (23 from V2 plus `chr_frag`) and 7 diagnostics | Sections 5, 5.8 |
| One feature added: `chr_frag`, the fragmentation of `M_chr`. It carries the signal `chr_tophat_edge` was read as | Sections 5.2, 5.9 |
| The validity method: three readings (within click, at matched `od`, at matched `h_chr_mean`) with ROI-bootstrap CIs reported; the V2 sign rule unchanged; a `beyond_darkness` label. `auc_od` is V2's `auc_dark_stratified` under a new name (equal on every column). Dropped: `auc_area_tertiles`, `auc_hole_*`, `auc_absdev`, the scanner splits and the depth-300 table | Section 6 |
| What the columns measure, corrected from the V2 run's readings: `chr_tophat_edge`, `h_chr_lap`, `chr_area_um2`, `dark_area_frac_101`, `ring_peak_prominence` and others | Section 5.9 |
| This run is built from the committed V2 rows, in a new notebook, with expected values to reproduce | Sections 0.1, 2, 6.4, 8, 9 |
| The V2 run's and the audit's failures added to the record of what went wrong | Section 1 |
| Departures from the fix plan: `chr_frag` gets its own montage, whose row 9 (the most fragmented mitoses) no V2 montage shows; the plan (3.2) said no new montage. `chr_tophat_edge` stays a feature beside `chr_frag`: the plan's "don't carry both" is a downstream choice, recorded in 5.9, and this spec selects nothing. The plan's tie shuffle before forming strata is not needed, because neither covariate ties within a click | Sections 5.9, 6.1, 8 |

**Review corrections (2026-10-09).** A reviewer implemented Section 0.1's steps 2–3 and Section 6 from this text
alone. Every cell of 6.4 and every carried-over V2 value reproduced exactly. These were then corrected:

| what changed | where |
|---|---|
| The bootstrap matrix is `BOOT_W`. It was `W`, which shadowed the imported worker that the montage step reads | 6.2, 6.3, 6.4 |
| The `h_chr_lap` numerator. The 0.572 was the audit's `lap_sd`, which skips the clipped-pixel exclusion; the column's own numerator reads 0.541 | 5.9 |
| `chr_tophat_edge`'s single-piece CI, the `od` ranker CI and `chr_area_um2`'s mask-free controls, recomputed with this spec's bootstrap and rule. Two labels marked *by eye*, which had been called measured | 5.9, 6.1 |
| Pool numbers that were depth-300 or rounded loosely: 228 → 49 rows with `chr_hole_frac` > 0.3, 87 → 89 % neighbour NaN, 1.6 → 1.8 % truncated windows, and the per-ROI spread range | Sections 1, 4, 5.8 |
| Where the `validity.csv` columns this run does not recompute live | Section 8 |
| The fix plan's Section 1 rows 4, 6, 9, 10, 11, 17 and 18, which the audit errata copies | the fix plan |

**The V2 revision (2026-10-09).** It applies the 35-item plan in `CANDIDATE_FEATURE_EXTRACTION_EDIT_PLAN.md`. Two
reviewers agreed on that plan, and the user approved it. The plan holds the evidence behind each item.

| what changed | where |
|---|---|
| 8 feature columns dropped: `chr_roughness`, `chr_n_pieces`, `chr_minor_um`, `h_chr_glcm_contrast`, `e_over_h_chr` (with `M_chr_clip`), `ring1_frac_glass`, `e_ring23_mean`, `nuc_area_frac_101` | Sections 4, 5 |
| 5 columns become diagnostics: `h_ring1_mean`, `h_ring2_mean`, `h_ring3_mean`, `nbr_nuc_area_med_um2`, `chr_edge_sharpness_per_um` | Sections 5.3, 5.8 |
| 1 column added: `dark_concavities_per_um`, on the new mask `M_dark`, with its diagnostic `dark_touches_window`. Result: 24 features and 9 diagnostics | Sections 4, 5.2 |
| A mandatory end-to-end C1 test through the verbatim refinement, with one deciding cell | Section 5.2 |
| Directions: texture and `chr_tophat_core` from −1 to 0; the cytoplasm hue angle (which replaces the E/H ratio) gets −1 | Sections 5.0, 5.4, 5.5, 5.8 |
| `chr_solidity` added to the lobulation test. Lashen 2021 added for C1. The outline columns relabelled as measuring the dense core's ragged outline, and `chr_tophat_edge` as thin dark structure | Sections 5.0, 5.2, 5.8 |
| Definitions changed: the brown rule gains `G − B > 15`; `dark_area_frac_101` excludes the candidate's own selection; `ring1_frac_nuclear` renamed `ring1_frac_above_tnuc`; the top-hat columns' `small piece` rule is whole `M_chr` again | Sections 5.2, 5.5, 5.6, 5.7, 5.8 |
| New identity columns `matched_unanimous` and `matched_n_mitotic_votes` | Sections 4, 5.1 |
| Pinned: the `directions.csv` serialisation, the 101 px offset's rounding, and the brown-rule montage ROIs. Recorded: the 2026-10-09 look at real data, which unblinds the spike columns | Sections 2, 4, 6, 8 |
| Added at the user's request: measured orientation per feature (`orientation_recommended`, sign consistency across ROIs and domains, mitoses in the top 10 against the bottom 10 per click), montages sorted by the oriented value, and a findings document | Sections 2, 6, 8, 9 |

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
| Wording pinned: ROI ids, `is_mitosis`, `image_id`, the scanner table, `regionprops` input type, NaN reasons for empty rings and small texture supports, which reason column `chr_edge_sharpness_per_um` fills, the strata and `redundancy.csv` row sets, the count column names, `synthetic_pass`, and when `directions.csv` is written. Two readings chosen where the text allowed two: `brown_frac_surround` divides by the surround's pixels, not the whole window's; the neighbour unit test is a fully specified 3 × 3 grid with 8 neighbours (it had said 12, with no layout). A third reading, one `small piece` rule for every column, was reverted in V2 | Sections 2, 3, 4, 5.1, 5.2, 5.4, 5.6, 5.8, 6, 8 |

## 0. What this is, and what it is not

**The stage** (as a full run does it; this run re-uses V2's rows, Section 0.1). For every one of the 294 recorded `mo_pad61` clicks on the 49 ROIs:
1. take the candidates of production's A3 search: the 100 highest peaks, then NMS (Section 3.1);
2. measure the features of Section 5 on each;
3. check that each feature measures what it claims to measure (Sections 2 and 6);
4. store every raw value.

**Not the job.** No ranking arms, no feature combinations, no P@K comparison between rankers, no selection rule, no
verdict, no production change. The downstream analysis decides all of that from your tables. If you find yourself
writing a rank average, stop.

**Two deliverables of equal weight:**
1. the feature tables (Section 8);
2. the evidence for what each feature measures, and a record of what you found wrong with this spec and with your own work (Sections 2 and 6).

A table of numbers without the second deliverable is a failed run.

### 0.1 This run, step by step

**Why this run measures nothing anew.** The V2 run (`d99e93a`) already measured 23 of this spec's 24 features on
exactly this candidate set.
- **Its rows with `in_recorded_pool` are the 100-peak pool.** V2 computed them as `nms_by_distance(centers[:MAX_PEAKS], ...)` and asserted that they are its deep list's prefix.
- **Gate B matched them with the recorded A3 run:** 293 clicks exact, one adjacent-pixel near-tie.
- **Their features came from code that was checked.** V2's synthetic suites reproduced every checked value, and the audit reproduced every AUC.
- **The 24th feature already exists.** `chr_frag` is in the audit's committed per-location table, computed on V2's own `M_chr`.

So this run re-uses those rows. It writes new tables and a new findings document.

**Execute, in this order:**
1. **Set up.** Create the notebook `production_hematoxylin_only/multiotsu_pad_arms/candidate_features_a3_pool100.ipynb` (Section 9). It refuses to start if any of the three `results/candidate_features_a3_pool100_*.csv` tables or the findings document (Section 8) exists.
2. **Build the candidate table** (Section 8, `candidates.csv`), including the `chr_frag` join and its assert (Section 2.1).
3. **Produce the validity report** (Section 6). Compare it with the expected values (6.4) before going on.
4. **Draw the `chr_frag` montage** (Section 8), after the refinement SHA-1 assert (Section 2.1). Look at it.
5. **Write the findings document and the final markdown cell** (Section 8).
6. **Add the two text sections** to existing documents (Section 8).
7. **Report back** (Section 9).

**Do not execute in this run:**
- the search, the refinement or any feature code (the montage re-refines its own tiles only);
- the unit tests and synthetic suites;
- writing `directions.csv`;
- a smoke run;
- gates A–C and the flags of Section 7, except the one SHA-1 assert;
- the V2 run's 14 montages and its neighbour check;
- the spec JSON.

Every section that describes only a full run says so in its first line. "Mandatory" in Sections 4 and 5 means
mandatory for a full run.

**Runtime:** a few minutes. Reading `candidates.csv` (72 MB) and loading the montage ROIs take most of it. No worker
pool is needed.

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

**What that means for this stage.**
- **The question to ask of each feature is whether it carries information beyond darkness.** Section 6's two conditional readings measure that: at matched `od` and at matched `h_chr_mean`. Report them; do not select on them (Section 0).
- **It is measured on the 100-peak pool only** (user, 2026-10-09).
  - The failing ROIs' readings rest on 412 pool mitosis rows.
  - Their deeper mitoses are not measured.
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

**What went wrong in the V2 run (`d99e93a`) and in its audit (`a4a6c0f`).** Verified on 2026-10-09; the numbers are in
the fix plan.

| what | what was wrong | what is right |
|---|---|---|
| The headline numbers | They came from a depth-300 AUC, 54 % of whose pairs were pool mitoses against beyond-pool others, so it tracked search rank. Five columns changed label between it and the pool | Only the pool is analysed (Section 3) |
| "Beyond darkness" | It meant "within `od` quintiles". `od` holds dark mass, not core darkness: its darkest 260 px include surround on 36 % of cores, and dark neighbours | Two conditional readings, at matched `od` and at matched `h_chr_mean` (Section 6.1) |
| `chr_tophat_edge` | Read as failing "its claim and its direction" (x-1 at depth 300, 0.389), its high end by eye as "reticulated chromatin of pale nuclei, small fragments and pigment" | Null on single-piece `M_chr` (0.501). Its signal is fragmentation, which `chr_frag` measures (5.9) |
| `h_chr_lap` | Read by eye as "thinness plus texture" | Relative texture; its x-1 holds within area quintiles and on single-piece `M_chr` (5.9) |
| `dark_area_frac_101` | "Among equally dark candidates, mitoses have fewer dark neighbours" | Partly `od`'s own neighbour pixels: 0.407 at matched `od` against 0.436 at matched `h_chr_mean` (5.9) |
| The neighbour columns | NaN on 89 % of accepted pool rows; where defined, they tracked `chr_area_um2` (ρ 0.78) | Dropped (5.8 d) |
| The audit's darkness diagnosis | It said `od` strata "do not hold darkness fixed" and proposed `dark50` and `disc_h` | Quintiles of about 20 rows leave every darkness measure at 0.57–0.59 within its own strata. The replacements hold darkness no better. What matters is which covariate (6.1) |
| The audit's `chr_tophat_edge` statistic | 0.524 "within fragmentation quintiles" | A tie-order artefact. 72 % of rows are tied at 0, and `rank(method='first')` on rows in `cand_rank` order made those strata score strata |
| The audit's `e_h_angle` claim | "The pooled null hides strong ROI heterogeneity" | The SD of the per-ROI AUC is 0.09–0.16 on 21 of the 24 features (`brown_frac_surround` 0.02, `chr_frag` 0.08 and `pair_parallel_deg` 0.20 are outside it). `e_h_angle`, at 0.148, is not distinctive |

**The general failures:**
- a mechanism read from a montage and never measured;
- a stratifier built on a tie-heavy column;
- a pooled statistic whose pairs came mostly from two different populations.

Section 6 rules out the last two. Section 2 requires a measurement before any mechanism is stated.

## 2. What you must check, and when

Every deviation from this spec goes in the findings document's "What I changed" section, with its reason. There is no
separate log file in this run.

### 2.1 This run

1. **The `chr_frag` join is new glue, so check it.**
   - Read both files with `float_precision='round_trip'`, join on the float keys `(file_name, cx, cy)`, many-to-one.
   - Assert that every accepted row gets a non-NaN `chr_frag` and every refused row a NaN (checked on disk: 0 failures either way).
   - A join that silently NaNs rows through a float-key mismatch is the realistic failure.
2. **Reproduce the expected values (6.4).**
   - The validity table must reproduce Section 6.4's table to the printed decimals, including the labels.
   - The carried-over quantities must equal V2's committed values to float rounding.
   - A difference means your code differs from this text. Find out why before continuing.
   - Never tune a seed, a constant or a tolerance until it passes.
3. **The refinement SHA-1, before the montage re-runs the refinement.**
   - Assert `W.REFINE_SHA1 == W.REFINE_SHA1_EXPECTED` (`21c43a20db00dccf7ac3bab88aee423c345441b6`), with `W` the imported V2 worker.
   - The refinement is extracted from a notebook cell that could have changed since V2.
4. **Look at the montage** with the Read tool (Section 8). Write one sentence on each of:
   - the high end;
   - the low end;
   - what the most fragmented mitoses are.
5. **Read the validity table.** For each feature, write one sentence:
   - what it measures (5.9, or V2's `claims_to_measure` where 5.9 is silent);
   - what its three readings say.

   Do not soften either.

**Everything here is in-sample.**
- **The five unblinded columns.** `chr_concavities_per_um`, `dark_concavities_per_um`, `chr_solidity`, `chr_radial_hf` and `chr_tophat_edge` were measured on this pool before V2's `directions.csv` was frozen.
- **`chr_frag`** was found on this pool on 2026-10-09.
- **No unseen evidence remains.** Depth beyond the pool is out of scope, so none of these six has a test of its direction.
- **In-sample orientation is accepted.** The user accepts that orientation is chosen on these 49 ROIs; there is no held-out set (`images/extra_valid/testing_set/` is part of the 49).
- **The hypothesis record.** V2's frozen `directions.csv` stays the record of the hypotheses. Its frozen rows hash to `620d50f3…`, as V2's notebook asserts; the committed file hashes to `63c00a7b…` because it adds the three `spec_` amendment rows. Section 6 reports the measured orientation beside them.

**When the spec, the repo or the committed data disagree:**
- **The repo and the committed data win.** If this spec contradicts them (`midog_utils`, `DECISIONS.md`, the V2 tables), log it in "What I changed".
- **One exception: the floor clip.** It is adopted (the D3 amendment of 2026-10-02) but is not yet in `midog_utils.production`. Its absence there is not a contradiction to resolve in the repo's favour.
- **If a check of Section 2.1 fails,** stop and report. Do not adapt.

### 2.2 A future full run (not executed in this run)

A full run measures new data with the code of Sections 3–5. Before any real data it must do three things:
- **Reproduce every unit test and synthetic suite in Section 5,** including the `chr_frag` unit test (5.2).
  - V2's `results/candidate_features_a3_top300_synthetic_tests.csv` holds the values this code gave.
  - A value that differs from a checked value means the implementation differs from the text.
- **Copy V2's `directions.csv`.** Drop its rows for the six neighbour columns, and add a `chr_frag` row with direction 0 (post hoc, Section 5.9). Freeze it before the first smoke execution.
- **Do a smoke run** on 4 ROIs (`337`, `509`, `201`, `211`). Print the projected time of the full run before starting it.

It then runs the gates of Section 7.

## 3. Inputs, the candidate set, and what is held fixed

*Defines the stage for a full run. This run takes its rows from V2 (Section 0.1).*

**The machine.** The recorded runs this spec compares against were made on the Intel Mac. On the current Apple M5 Max,
colour deconvolution differs in the last float bits (a different libm and BLAS), so `hem` differs by up to about 1e-6
on every ROI. Measured 2026-10-08 on the production pipeline: the same detections in the same order, with scores
differing by up to 1.2e-5 and `od` by up to 6e-8. So every float compared with a recorded value uses an **absolute
tolerance of 1e-4**, and positions, ranks and labels must match exactly except for float near-ties (gate B, Section 7).

**Held fixed, copied verbatim** (nothing here is re-tuned):
- The 49 ROIs: `images/extra_valid/` (14) and `images/extra_valid/testing_set/` (35), 7 per domain. ROI ids in this
  spec drop the extension: `211` is the `file_name` `211.tiff`.
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
  and `exec` them. There is no new search code: `detect` and `search_arm` run as extracted (Section 3.1).
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

**The candidate set.** For each click: every post-NMS row of production's A3 search (Section 3.1), in descending
score. That is 79–100 rows per click and 29,025 rows in all. `cand_rank` is 0-based.

Do not drop any row for any reason. Refused candidates stay, with their refusal recorded and their features NaN. The
refusal policy is a downstream decision.

**Why only the pool (user, 2026-10-09).** The analysis is limited to the candidates production ranks. It has two
costs:
- **No unseen evidence:** the five unblinded columns and `chr_frag` have none (Section 2.1);
- **Little failing-ROI data:** the failing ROIs rest on 412 pool mitosis rows, and their deeper mitoses are not measured.

A deeper analysis comes back only if the user asks for it.

**Inputs for this run.** Read every CSV with `float_precision='round_trip'`, and `candidates.csv` also with
`low_memory=False`.

| file | what this run takes from it |
|---|---|
| `results/candidate_features_a3_top300_candidates.csv` | The rows with `in_recorded_pool`: identity, labels, features, diagnostics and NaN reasons |
| `results/candidate_features_a3_top300_audit_locations.csv.gz` | `chr_frag`, by `(file_name, cx, cy)` |
| `results/candidate_features_a3_top300_validity.csv` | The carried-over checks of Section 6.4 |
| `results/candidate_features_a3_top300_redundancy.csv` | The carried-over check of Section 6.4 |
| `results/candidate_features_a3_top300_directions.csv` | `direction_hypothesis` and `claims_to_measure` per feature |
| `production_hematoxylin_only/multiotsu_pad_arms/candidate_features_a3_top300.ipynb` | Cell 27's `tied_expectation` and `top_bottom`, extracted with `ast` (Section 6.2) |
| `production_hematoxylin_only/multiotsu_pad_arms/candidate_features_a3_top300_worker.py` | Imported as `W`, for the montage: `W.montage_windows` and `W.REFINE_SHA1`. Importing it extracts its verbatim code from three notebooks, which takes seconds |

The montage's ROI images: `images/extra_valid/<file_name>` when `subset` is `original_14`, and
`images/extra_valid/testing_set/<file_name>` when it is `testing_35`.

**Other inputs, for a full run:**
- `databases/MIDOG++.json`, read as `_, annotations = ds.load_annotations(path)` (it returns `(images, annotations)`),
  then `ds.image_annotations(annotations, file_name)` with both categories;
- for gate B (Section 7), the recorded A3 run:
  - `results/tm_floor_definition_49roi_6seed_per_run.csv`, `arm == 'A3'`: `threshold`, `n_detections`;
  - `results/tm_floor_definition_49roi_6seed_pool.csv.gz`, `arm == 'A3'`: the recorded pools.

### 3.1 The candidate search: A3, production's 100-peak pool

*Defines the stage for a full run. This run takes its rows from V2 (Section 0.1).*

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
- **The score and peaks:** `TM_CCOEFF`; peak threshold median + `DEEP_FLOOR_Z`·MAD; `PEAK_MIN_DISTANCE`; the 100
  highest peaks kept (`MAX_PEAKS`, D9).
- **The 30 px border band:** `band_mask(hem.shape, PAD)` with `PAD` = 30, made invalid and its scores floored
  (`BAND_MODE = 'floored'`).
- **NMS** at `ev.radius_px(mpp)` over those 100 peaks only, then ranking by score.
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
- **Consequence for this spec.** The features read the unclipped channel (Section 4), so the surround features
  (`ring_peak_prominence`, `ring1_frac_above_tnuc`, `ring1_frac_pale`, `e_h_angle_ring23_deg`) and the `h_ring*_mean`
  diagnostics carry the white-surround signal the clip removed from the score. That is correct for measurement, and
  not a reason to clip the features. Read every surround feature beside its correlation with `white_frac_gray215_61`
  in `redundancy.csv`.

**The pool: 100 peaks, then NMS.** This is `detect` from cell 5 exactly as written, with production's `MAX_PEAKS` = 100 (D9):
1. extract the score peaks above the peak threshold;
2. keep the 100 highest (`tm.extract_peaks(..., prod.MAX_PEAKS)`);
3. apply NMS at `ev.radius_px(mpp)` to those 100 only;
4. compute `od` on the rows NMS keeps.

**What it gives on these clicks:**
- Every one of the 294 clicks has at least 16,786 peaks above the threshold (V2's uncapped count), so the cap always binds. The effective threshold is the 100th peak's score.
- NMS leaves 79–100 candidates per click.

**This run** does not search. The V2 rows with `in_recorded_pool` are this pool (Section 0.1).

**Removed in V3:** `detect_deep`, `search_arm_deep`, the uncapped peaks, the prefix NMS and the 300-row depth.

**There are no border candidates.** The band is invalid in the search, so every candidate has
`edge_px = min(cx, cy, w − 1 − cx, h − 1 − cy) ≥ 30`, and its 51 px and 61 px windows lie inside the ROI. There is no
padded re-run of any candidate. The search's own internal pads (the correlation pad, `OD_PAD`, the template read pad)
are part of the verbatim search and supply no feature pixels.

## 4. Per-candidate pipeline

*Defines the stage for a full run. This run takes its rows from V2 (Section 0.1).*

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
- **101 px:** the window of the pair and context columns. It crosses the ROI edge for candidates 30–49 px
  from it (1.8 % of pool rows). It is **truncated** to the pixels inside the ROI, never padded. `(y0, x0)` is its
  origin in ROI pixels. With `ix = int(round(cx))` and `iy = int(round(cy))`, as `read_padded_patch` uses, the 51 px
  masks (`M_chr`, `M_sel`) sit in it at row and column offset `(iy − 25 − y0, ix − 25 − x0)`. That offset is 25 only
  when the window is not truncated. `hem_101` is its unclipped H.

**Masks**, in the 51 px window:
- `M_chr`: `ref['selection_mask'] & (ref['class_map'] == 2)`, the undilated top multi-Otsu class inside the selection.
  Boundary features use this mask, because the refinement's 1 px dilation fills exactly the notches a boundary feature
  measures. Store `chr_touches_window`, true when `M_chr` reaches the window's edge; an outline cut by the window is
  not the object's outline.
- `M_sel`: `ref['selection_mask']`, the dilated selection, as the template uses it.
- `M_dark`: a lower-threshold chromatin mask that keeps paler chromosome arms. Build it in three steps:
  1. take the 4-connected components of `H_w > t_dark` that overlap `M_chr`;
  2. keep the largest 4-connected piece, with a tie going to the lower label;
  3. fill its holes with `scipy.ndimage.binary_fill_holes`.

  Only `dark_concavities_per_um` (5.2) reads it. Store `dark_touches_window`, true when that piece reaches the window's
  edge.
  - **Why it exists.** On the recorded pool's mitoses, the window's upper Otsu threshold, which bounds `M_chr`, sits at
    a median of 0.72× the body's mean H (q10–q90 0.67–0.79), and `t_dark` at 0.44× (0.34–0.53). So an arm paler than
    about 0.7× the body leaves `M_chr`, but stays in `M_dark` down to about 0.45×.
  - **Its cost.** At the lower threshold, neighbours and debris fragments can merge into the piece. The window also
    cuts the piece more often than it cuts `M_chr`: on the recorded pool it cuts 21 % of other rows against 8 % of
    mitoses.

**`t_nuc` and `t_dark`** are the ROI-wide thresholds of Section 3. `t_nuc` is the search floor, so "above `t_nuc`" means
"untouched by the clip", and "below `t_nuc`" lumps glass, lumen, halos and pale cytoplasm together. About 45 % of a
typical ROI lies above it, cytoplasm included.

**Rings** for the surround features: `ring_k` for k = 1, 2, 3 is the band of pixels whose distance from `M_chr`
(`scipy.ndimage.distance_transform_edt(~M_chr) * mpp`) lies in `((k−1) µm, k µm]`, clipped to the 51 px window.
`M_chr` pixels are at distance 0, so no ring contains them. A large `M_chr` can leave a ring with no pixels inside the
window. Then every column that reads that ring is NaN with reason `empty ring` (`nan_reason_ring`, Section 8): the
ring's `h_ring_k_mean`, `ring_peak_prominence` if any of the three rings is empty, and the `ring1_frac_*` shares and
`chr_edge_sharpness_per_um` if `ring_1` is empty. (The cytoplasm columns of 5.5 keep their own `no cytoplasm` reason,
which an empty ring already triggers.)

**Clipped pixels.** A pixel is clipped when any RGB channel is 0 or 255.
- True clipping occurs on the 3D Histech ROIs: on `249`, `435`, `230` and `460`, 5.7–7.5 % of GT-mitosis `M_chr`
  pixels have G at exactly 0, with about twice as many at 0 as at 1 (a pile-up). `rgb2hed` clamps 0 to 1e-6, so G
  going from 1 to 0 moves H by −0.04 and E by +0.68: those pixels' colour and local contrast are artefacts.
- `h_chr_lap` excludes clipped pixels as Section 5.4 says, and the cytoplasm ring excludes them as Section 5.5 says.
  `NS_MIN_PX` = 10 is the minimum pixel count for both.
- **Do not exclude pixels merely because a channel is below 20.** On GT mitoses those are G < 20 dense blue-purple
  hematoxylin, mean RGB about (27–36, 9–16, 47–88), not saturation. Excluding them dropped 60 % of a mitosis's `M_chr`
  against 37 % of other candidates' on the recorded pool: a subset chosen by label and darkness.
- Section 5.8 (a)'s `chr_resid_frac` is the one exception and keeps [20, 235] on its own region, because black pigment
  is low in all channels. Its support is therefore label-dependent too (median `chr_ns_share` 0.76 for mitoses against
  0.44 for top-10 false positives on 7 ROIs), so its AUC is read beside `chr_ns_share`.

**Labels.** Match each click's list in score order with `ev.bucket_detections` against the ROI's ground truth
minus the click's own annotation. Store `bucket`, `matched_ann_id`, `is_mitosis`. `is_mitosis` is
`bucket == ev.HUMAN_CORRECT_LABEL`. Every other row, whether matched to a look-alike or unannotated, is an "other" row
in Section 6.
- **Expert agreement.** Store `matched_unanimous` and `matched_n_mitotic_votes`. Join them on
  `matched_ann_id = ann_id` from the `unanimous` and `n_mitotic_votes` columns that `ds.load_annotations` already
  returns, and leave them empty on unmatched rows. A quarter of MIDOG++ mitoses (3,020 of 11,937) are 2-to-1
  decisions. The downstream analysis computes any split by agreement; this spec adds no validity column for it.
- **Labels that don't depend on matching.** Independently of matching, store `nearest_mitosis_um` and
  `nearest_lookalike_um`: the distance to the nearest annotation of each category, excluding the click's own. The
  downstream analysis needs labels that do not depend on greedy order.

## 5. Features

*Defines the stage for a full run. Every "mandatory" test, check and montage, every permitted threshold change, every
spec-JSON entry and every smoke-run step in this section belongs to a full run;
the V2 run did them on this code, and this run executes none of them (Section 0.1). Section 5.9 records what the
columns measured on the V2 pool.*

**Rules that apply to every feature:**
- Store the **raw measurement**. Never negate, rank, standardise or direction-encode a column. Directions are
  hypotheses and live in `directions.csv` (Section 8).
- Physical units where a quantity has one: µm, µm², optical density. Nothing may depend on `mpp` implicitly.
- Nothing may depend on the 51 px window's own intensity histogram except through `M_chr` and `M_sel`, which are the
  production refinement and are kept for that reason. (`M_dark` uses the ROI-wide `t_dark`, not the window's
  histogram.)
- NaN where undefined, with the reason in that family's `nan_reason_*` column (Section 8); never a sentinel number.
- Column names are `snake_case`, and every column is described in the notebook in one line with its unit.
- Every column has a role in `directions.csv`: `feature` (the 24 below), `diagnostic` (the 7 below), or `spec` (a
  definition you replaced, Section 2).

**The 24 features, by section:**
- 5.2: `chr_area_um2`, `chr_elongation`, `chr_solidity`, `chr_concavities_per_um`, `dark_concavities_per_um`,
  `chr_hole_frac`, `chr_frag`, `chr_radial_hf`, `pair_sep_um`, `pair_area_ratio`, `pair_parallel_deg`, `pair_axis_deg`;
- 5.3: `h_chr_mean`;
- 5.4: `h_chr_lap`;
- 5.5: `ring_peak_prominence`, `ring1_frac_above_tnuc`, `ring1_frac_pale`, `white_frac_gray215_61`,
  `e_h_angle_ring23_deg`;
- 5.6: `brown_frac_surround`;
- 5.7: `dark_area_frac_101`;
- 5.8: `chr_resid_frac`, `chr_tophat_edge`, `chr_tophat_core`.

**The 7 diagnostics.** They are stored, and are not validity rows.
- `chr_touches_window` and `dark_touches_window` (Section 4);
- `h_ring1_mean`, `h_ring2_mean` and `h_ring3_mean` (5.3);
- `chr_ns_share` (5.8 a);
- `chr_edge_sharpness_per_um` (5.8 c).

**Named column sets** used elsewhere in this spec:
- **boundary columns:** `chr_solidity`, `chr_concavities_per_um`, `dark_concavities_per_um`, `chr_radial_hf`,
  `chr_tophat_edge`;
- **spike columns** (the mask-level lobulation test, 5.2): `chr_radial_hf` and `chr_concavities_per_um`. On a given
  mask, `dark_concavities_per_um` is the same computation as `chr_concavities_per_um`, so it is tested end to end
  instead;
- **texture column:** `h_chr_lap`;
- **101 px columns:** the four `pair_*` and `dark_area_frac_101`;
- **montage columns** (Section 8): `chr_radial_hf`, `chr_tophat_edge`, `chr_concavities_per_um`,
  `dark_concavities_per_um`, `chr_hole_frac`, `ring_peak_prominence`, `chr_solidity`, `chr_edge_sharpness_per_um`,
  `h_chr_lap`, `chr_tophat_core`, `ring1_frac_pale`, `e_h_angle_ring23_deg`, `chr_resid_frac`, `brown_frac_surround`, and in V3 `chr_frag`.
  These are the columns that claim a visible property, which a number alone cannot confirm. This run draws only
  `chr_frag`'s montage (Section 8); the V2 notebook holds the other 14.

### 5.0 What pathologists use, and which columns target it

Use this table as the `rationale` of `directions.csv`. A direction given here or in the section that defines a column
is the hypothesis to write down; a 0 stays 0 unless you log a reason before the smoke run. Quotes are verbatim from the
sources listed after the table.

| # | criterion | evidence | columns | direction |
|---|---|---|---|---|
| C1 | **Short spikes or rods on the surface of the chromatin**, the most reliable cue | VCGP 2021: "The most reliable microscopic feature to identify MF and AMF are the short spikes or blunt projections of chromosomes on the surface of aggregated nuclear material"; also "dark purple staining nuclear aggregates with an irregular contour (short rods or spikes of chromosomes, described by van Diest as hairy projections)". Meuten 2021: "short 'rods' of chromosomes protruding from the surface of aggregates of nuclear material". Al-Janabi 2013 (restating the van Diest and Baak criteria): "clearly visible hairy extension of nuclear material". Veta 2015, on weaker algorithms' false positives: "dark objects that lack other characteristics of mitotic figures such as protrusions around the edges". VCGP QRG 2022, under "Differentiate MF from Mitotic-Like-Figures (MLF)": "Surface contour: MLF are smooth, MF have projecting rods/spikes". Lashen 2021 (H&E breast cancer, light microscope against WSI): "Fifty-eight per cent of mitoses showed absent hairy-like projection in WSI; however, 89% retained their ragged nuclear border, which distinguished them from mimickers including apoptotic cells, lymphocytes and dark elongated hyperchromatic structures" | The outline columns `chr_concavities_per_um`, `dark_concavities_per_um` and `chr_radial_hf` measure irregularity of the dense core's outline, which is consistent with the ragged border that survives on WSI; that has not been shown to be the same cue. `chr_tophat_edge` measures thin dark structure within 1 µm of the core, which includes fine arms (5.2, 5.8 b). On the V2 pool, `chr_tophat_edge` measured `M_chr` fragmentation instead (5.9) | +1 for all four |
| C2 | **No nuclear membrane** | Al-Janabi 2013: "absence of the nuclear membrane". VCGP 2021: "lack of a nuclear membrane". Veta 2015: "hyperchromatic objects that lack a clear nuclear membrane" | No column measures C2 directly. `ring_peak_prominence` and `chr_hole_frac` see parts of it; Section 5.5 says what each can see | 0 for both (5.5) |
| C3 | **Dark chromatin**, necessary but not sufficient | VCGP 2021: "dark purple staining nuclear aggregates". Veta 2015: "many other cellular components can have a similar hyperchromatic appearance" | `h_chr_mean` | +1; it is why Section 6 reads every feature at matched darkness |
| C4 | **Configuration**: a clot or ring, a plate, separate clots | VCGP 2021: "Prometaphase: A dark cluster or ring of chromosomes with variably protruding chromosome rods and spikes"; "Metaphase: An equatorial plate (linear or ring-shaped aggregate with short rods or spikes)"; "Anaphase: Two separate aggregates, variable distances apart, linear with projections/spikes"; "Anaphase and telophase MF are counted as one". Al-Janabi 2013: "either clotted (beginning metaphase), in plane (metaphase/anaphase), or in separate clots (telophase)". Meuten 2021: identifying the phase is "not necessary" | `chr_elongation`, `chr_hole_frac`, the `pair_*` columns (5.2) | 0 for all (below, and 5.2) |
| C5 | **Imposters are smooth, round or oval dark bodies, one or several** | VCGP 2021, mitotic-like figures: "Smooth contours, circular or oval, hyperchromatic bodies (single or multiple)". ICCR, in its own words (its method source is Elston & Ellis 1991): "hyperchromatic, karyorrhectic, or apoptotic nuclei are excluded". Cree 2021: "Strict criteria should be applied for accepting a mitotic figure in distinction from a pyknotic or apoptotic nucleus" | `chr_solidity`. (`chr_edge_sharpness_per_um` is a diagnostic, 5.8 c) | −1; confound: touching or lobulated nuclei also score low (5.2) |
| C6 | **Prophase is not counted on H&E** | VCGP 2021: "MF are counted from prometaphase through telophase … Prophase cannot be reliably identified in histological preparations". Lin 2024: dark chromosome aggregates "are most likely in the final stages of prophase, commonly referred to as prometaphase". Meyer 2005 recommends "that prophases be excluded from mitotic counts performed for cancer grading", though its own study counted them. MIDOG++ delegates its criteria to Donovan et al. 2021 and Ibrahim et al. 2022. MIDOG 2022's PHH3-assisted ground truth (not the H&E labels) did not annotate immunopositive cells lacking H&E morphology, "mostly early prophase MFs", and its largest subcategory of figures found only with PHH3 (within "borderline morphology", 87.58 %) was "prophase with resemblance to apoptosis" (31.23 %) | none directly; context for `ring_peak_prominence` (5.5) | — |
| C7 | **Named imposters**: immune cells, stromal cells, two cells in close proximity, karyorrhectic debris | Lin 2024, on which mimic each true mitotic figure resembled (most often the low-agreement ones), not on false positives: "karyorrhectic debris, immune cells, stromal cells, and two cells in close proximity"; "lymphocytes were most likely to resemble mitosis in prometaphase, and stromal cells were likely to mimic mitosis in metaphase" | `chr_area_um2`; the lobulation test (5.2) | 0 |
| C8 | **Pigment** | Aubreville 2023 (MIDOG++): melanoma "pigment particles that contribute additional imposter structures ... although they have a different chromaticity (brown)" | `chr_resid_frac`, `brown_frac_surround` | −1 |
| C9 | **Cytoplasm and surround** | Ikpatt 2002, restating the criteria of Baak & Oort: "The cytoplasm of the mitotic cells was often larger during mitosis than in the resting cells". No verified source states a clear zone around mitotic figures in H&E. A clear halo is described around apoptotic hepatocytes in rodent liver (NTP atlas). VCGP 2021, mitotic-like figures: "MLF cytoplasm tends to be eosinophilic, however color is not reliable". VCGP QRG 2022: "Cytoplasm color (not reliable). MF tend to have more amphophilic/basophilic cytoplasm (increased RNA), MLF tend to have more eosinophilic cytoplasm". VCGP 2021, mitotic-like figures (item 4): "Eosinophilic cytoplasm can occur with denaturation of cytoplasmic proteins, loss of mRNA in degenerating, apoptotic or necrotic cells". NTP atlas (rodent liver): apoptotic cells have "condensed, deeply eosinophilic cytoplasm" | Lightness: `ring1_frac_above_tnuc`, `ring1_frac_pale` and `white_frac_gray215_61` (5.5); the `h_ring*_mean` profile is stored as a diagnostic (5.3). Cytoplasm colour: `e_h_angle_ring23_deg` (5.5) | 0 for lightness; −1 for `e_h_angle_ring23_deg`, with "color is not reliable" as its confound |
| C10 | **Context** | VCGP 2021: "A structure within a region of necrosis ... is more likely to be MLF"; "A structure within a region of high mitotic activity is more likely to be a MF". VCGP QRG 2022: "Context: In areas of poor cell preservation, MLF may be more likely" | none; necrosis context is not measured (Section 10). `dark_area_frac_101` measures crowding, the mechanism of the crowded failing ROIs, not this criterion | 0 |

**Elongation.** A side-on metaphase plate is elongated, but so are spindle-cell and stromal nuclei, which Lin 2024 names
as metaphase mimics; `462` and `410` are sarcomas. Pathologists need elongation together with surface spikes, so
`chr_elongation` has direction 0.

**Texture and colour.**
- **Texture: direction 0.** `h_chr_lap` (5.4) and `chr_tophat_core` (5.8 b) both get direction 0.
  - The only text behind a direction is one VCGP 2021 clause about resting nuclei (item 5.iv): "Cells in a resting
    phase have more heterochromatin with a more clumped chromatin pattern and a smaller nucleus". The same item
    describes pyknotic and karyorrhectic imposters only as "dark, hyperchromatic (darker) condensed nuclei", which
    gives no texture direction.
  - `h_chr_lap` depends on size: through the verbatim refinement, a flat blurred disc reads 1.74 / 1.51 / 1.26 /
    1.23 at r = 5 / 8 / 12 / 16 px.
  - Record two confounds in `directions.csv`. A prometaphase clump has internal structure from its chromosome arms. A
    pixel-noise floor may rise with darkness, depending on the noise model: with constant intensity noise
    `h_chr_lap` goes 2.1 → 3.6 over H 0.08–0.28, and with noise proportional to √I it stays flat.
  - Do not describe pyknotic nuclei as "homogeneous"; no source read for this spec says so.
- **Chromatin colour: no column.** No source gives a chromatin colour criterion; VCGP's "color is not reliable" is
  about imposter *cytoplasm*. Pigment (C8) is the exception.
- **Cytoplasm colour: direction −1** for `e_h_angle_ring23_deg`. VCGP 2021 item 4, the QRG 2022 and the NTP atlas
  state the tendency (C9). The reliability caveat is its confound.

**Independence.** VCGP 2021, Meuten 2021, Donovan et al. 2021, MIDOG 2022 and MIDOG++ share authors; the VCGP QRG
2022 names no authors and cites the first three. MIDOG++ also includes van Diest and Veta, who co-authored Al-Janabi
2013 (van Diest), Veta 2015 (van Diest) and Veta 2016. So every C1 and C2 source shares an author with the dataset's
group; Lin 2024 (C6, C7) does not. Their agreement is one school's view, not independent corroboration. Lashen 2021
(C1) is from the Rakha group, which wrote MIDOG++'s second criteria reference (Ibrahim et al. 2022). So it is
independent of the VCGP authors, but not of the dataset's criteria.

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
- **Lashen 2021:** Lashen A, Ibrahim A, Katayama A, Ball G, Mihai R, Toss M, Rakha E. Visual assessment of mitotic
  figures in breast cancer: a comparative study between light microscopy and whole slide images. *Histopathology*
  79:913–925. doi:10.1111/his.14543. Read as the abstract (Europe PMC).
- **Cisneros-Soberanis 2024:** *J Cell Biol* 223:e202403165. doi:10.1083/jcb.202403165. By serial block-face EM of
  cultured RPE1 cells, a single chromatid is "0.68 ± 0.04 μm" wide and a cohesed chromosome "1.23 ± 0.1 μm" (used in
  5.8 b). It bounds the physical width of a chromosome arm; it says nothing about what a scan resolves.

**Not verified, so used for nothing:** van Diest 1992's own text, including whether "two parallel, clearly separate
chromosome clots" count as one figure or two (open sources say both; VCGP 2021 counts anaphase and telophase as one
and notes "original reports indicated telophase should be counted as two"); Elston & Ellis 1991's own text; Donovan
et al. 2021 beyond its abstract; Ibrahim et al. 2022, *J Clin Pathol*; any published perimitotic clear zone in H&E.

**Not measured:**
- **The nucleolus.** A visible nucleolus marks a resting nucleus, but a nucleolus detector is out of scope (Section 10).
- **C2 directly.** See its row.
- **Necrosis context** (C10).

### 5.1 Identity and pipeline columns

`subset, domain, file_name, image_id, mpp, seed_index, seed_ann_id, click_cx, click_cy, cand_rank,
n_detections_recorded, cx, cy, score, od, dist_to_click_um, edge_px, accepted, refusal, mode, n_pieces, base_size,
bucket, matched_ann_id, is_mitosis, matched_unanimous, matched_n_mitotic_votes, nearest_mitosis_um,
nearest_lookalike_um, scanner, roi_group`.

- `dist_to_click_um` is the distance from `(cx, cy)` to `(click_cx, click_cy)`, times `mpp`.
- `n_pieces` is the number of refinement pieces the selection unioned.
- V2's `in_recorded_pool` is gone: every V3 row is in the 100-peak pool.
- `image_id` is the `id` of the entry in the `images` list returned by `ds.load_annotations` whose `file_name`
  matches.
- **`scanner`** is read from the `scanner` column of `images/extra_valid/MANIFEST.md` and
  `images/extra_valid/testing_set/MANIFEST.md`. Each file has more than one table, so use the table whose header has a
  `scanner` column. Map the values `3D Histech` → `3dhistech`, `Aperio CS2` → `aperio_cs2`,
  `Hamamatsu S360` → `hamamatsu_s360`, `Hamamatsu XR` → `hamamatsu_xr`. Assert that all 49 ROIs get one.
- **`roi_group`** is `failing9` for `211`, `242`, `524`, `546`, `462`, `527`, `410`, `245` and `289`, and `other40`
  otherwise. The nine are the ROIs below P@10 0.35 on both A3 axes (`low_precision_rois/`). The group was defined after
  seeing results, so every split by it is descriptive.

### 5.2 Geometry of the chromatin, on the undilated `M_chr`

**The piece.** Columns are computed on the largest 4-connected piece of `M_chr` (`skimage.measure.label(M_chr,
connectivity=1)`; a tie goes to the lower label) unless marked *whole*. They are NaN, with reason `small piece`, when
that piece, unfilled, has fewer than 30 px.
- **`small piece` has two rules in this spec.** The largest-piece rule above applies to the columns read from the
  largest piece (`chr_elongation`, `chr_solidity`, `chr_concavities_per_um`, `chr_hole_frac`, `chr_radial_hf`) and to
  `chr_edge_sharpness_per_um`. The two top-hat columns of 5.8 (b) read *whole* `M_chr`, and use "whole `M_chr` < 30 px"
  instead. `dark_concavities_per_um` applies the 30 px rule to its own `M_dark` piece.
- **Why the top-hat columns differ.** Clusters of small debris fragments pass the refinement with a whole `M_chr` of
  66–160 px but a largest piece of only 14–21 px. A largest-piece rule would make every such cluster NaN, which hides
  `chr_tophat_edge`'s highest "other" values (0.32–0.62), the small-body confound.
- **`regionprops`.** Every call takes `mask.astype(np.uint8)`: skimage 0.24 raises `TypeError` on a bool array.

**The filled piece and its outline.** `chr_elongation`, `chr_solidity`, `chr_concavities_per_um` and `chr_radial_hf`
read the piece with its holes filled (`scipy.ndimage.binary_fill_holes`). Pathologists judge the outer surface contour
(C1, C5), and on a ring-shaped piece `regionprops` counts the inner boundary too. The vesicular nucleus of Section 5.5,
whose `M_chr` is its rim, scored solidity 0.21 unfilled, the most "mitotic" value, against 0.99 filled.
`chr_area_um2` reads *whole* `M_chr` and `chr_hole_frac` reads the piece, both unfilled; the hole is `chr_hole_frac`'s
job. The **outline** is the longest contour of `cv2.findContours(filled.astype(np.uint8), cv2.RETR_EXTERNAL,
cv2.CHAIN_APPROX_NONE)`, as (x, y) points.

| column | definition | unit |
|---|---|---|
| `chr_area_um2` (*whole*) | pixel count of `M_chr` × mpp² | µm² |
| `chr_elongation` | `regionprops` `axis_major_length / axis_minor_length` | ratio |
| `chr_solidity` | `regionprops` `solidity` | ratio |
| `chr_concavities_per_um` | the number of convexity defects deeper than `CONCAVITY_MIN_UM` = 0.25 µm, divided by the outline length `cv2.arcLength(c, True)` × mpp. With `c` the outline as an `int32` (n, 1, 2) array: `d = cv2.convexityDefects(c, cv2.convexHull(c, returnPoints=False))`; a defect counts when `d[:, 0, 3] / 256 * mpp > 0.25` (the depth is fixed-point); `d is None` counts 0. If OpenCV raises (it can on an outline that touches itself), NaN with reason `convexity error`, and report the count. 0.25 µm is about 1 px, at the pixel staircase | per µm |
| `dark_concavities_per_um` | the same construction as `chr_concavities_per_um`, on the outline of the `M_dark` piece (Section 4), with the same `CONCAVITY_MIN_UM`. NaN, in `nan_reason_dark`: `small piece` when the piece has fewer than 30 px (empty included), or `convexity error`. A window cut is not a NaN; `dark_touches_window` records it | per µm |
| `chr_hole_frac` | 1 − area / area with holes filled, on the largest piece. A ring-shaped `M_chr` scores high; Section 5.5 says why that matters | ratio |
| `chr_frag` (*whole*) | 1 − (pixels of whole `M_chr`'s largest 4-connected piece) / (pixels of whole `M_chr`), with `skimage.measure.label(M_chr, connectivity=1)`. 0 when `M_chr` is one piece. Never NaN on an accepted row, because an empty `M_chr` is a refusal (Section 4), so it has no reason column. Added in V3 (5.9) | ratio |
| `chr_radial_hf` | `o = resample_closed(outline, 128)`; `r` = the distance of each point of `o` from the pixel centroid of the filled piece (x the mean column, y the mean row); `e = np.abs(np.fft.rfft(r − r.mean()))**2`; `e[6:65].sum() / (128 * r.mean()**2)`, so harmonics 6 to 64 inclusive. A `skimage` contour gives different values (disc 0.047 against 0.057 at r = 8), so use this one | ratio |

**`chr_frag` unit test (a full run).** On a 41 × 41 canvas:
- a 7 × 7 square alone gives 0;
- a 7 × 7 and a 5 × 5 square 3 px apart give 25 / 74 = 0.338;
- two 3 × 3 squares touching only at one corner give 0.5, because the pieces are 4-connected.

`chr_radial_hf` replaces the earlier `spiky`. Its denominator is the mean radius squared, not the total outline energy,
so a disc scores near zero instead of near one.

**What the outline columns claim, for `directions.csv`.**
- **`chr_concavities_per_um` and `chr_solidity`** measure the irregularity and the convexity of the dense core's outline.
- **`chr_radial_hf`** measures the high-harmonic outline energy of the dense core. Its confounds:
  - it fails the scale rule on every shape, because pixel-staircase noise falls with r²;
  - it reads like a disc when the arms are paler than about 0.7× the body (end-to-end test, below);
  - it overlaps `chr_solidity` (within-click ρ −0.67 on the V2 pool, `redundancy.csv`).
- **`dark_concavities_per_um`** measures the same irregularity on `M_dark`, which keeps arms down to about 0.45× the
  body.
  - **Rationale:** with synthetic arms at 0.7× and 0.5× body, `M_chr` reads 0–0.22 and this column 0.16–0.41.
  - **Confounds:** debris clusters merge into one ragged blob (0.20–0.23), and neighbours merge at `t_dark`. A
    window-cut outline reads about 0.02 per µm lower for both labels, and other rows are cut 2.5× as often. So compare
    it with `chr_concavities_per_um` only on rows where `dark_touches_window` is False.
  - **Status:** it was added after the 2026-10-09 look (Section 2).
- **All four** are consistent with the "ragged nuclear border" that Lashen 2021 found survives on WSI. That is a
  pathologist's visual judgement, and a defect count on a thresholded mask has not been shown to be the same cue.

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

Checked values (r = 8 / r = 16). For `chr_radial_hf`, an independent implementation reproduced the disc and the
12-lobed star to the third decimal, and the 2026-10-08 implementation of this text reproduced every value. The
`chr_solidity` values were computed on 2026-10-09 and reproduced by a second reviewer on six shapes.

| shape | `chr_radial_hf` | `chr_concavities_per_um` | `chr_solidity` |
|---|---|---|---|
| disc | 0.057 / 0.010 | 0.000 / 0.000 | 0.925 / 0.961 |
| 2:1 ellipse | 0.033 / 0.008 | 0.000 / 0.000 | 0.981 / 0.971 |
| 4:1 ellipse | 0.097 / 0.081 | 0.000 / 0.000 | 0.873 / 0.938 |
| 5-lobed star | 0.024 / 0.020 | 0.303 / 0.146 | 0.790 / 0.803 |
| 8-lobed star | 2.328 / 1.681 | 0.343 / 0.183 | 0.663 / 0.708 |
| 12-lobed star | 1.089 / 1.635 | 0.470 / 0.205 | 0.738 / 0.680 |
| spiked ellipse | 1.346 / 0.453 | 0.350 / 0.252 | 0.692 / 0.753 |
| touching discs | 0.180 / 0.158 | 0.088 / 0.043 | 0.831 / 0.861 |
| notched disc | 0.075 / 0.035 | 0.060 / 0.030 | 0.819 / 0.869 |

What the suite must show, and what the values say. Every result here is on ideal masks; the end-to-end test below
shows what survives the production mask. No synthetic body here is a ragged core, so the montages are the only check
that the columns see one.
- **Discs and ellipses score alike:** low on `chr_radial_hf` and `chr_concavities_per_um`, high on `chr_solidity`.
  So none of the columns is secretly an elongation measure.
  - `chr_radial_hf` has a partial elongation sensitivity (the 4:1 ellipse scores above the disc), because equal-arc
    sampling gives an elongated outline in-band energy even without pixels. Report it as such.
  - At r = 8, pixelation lowers the 4:1 ellipse's solidity to 0.873.
- **Spikes score beyond smooth shapes.** The 8- and 12-lobed stars and the spiked ellipse score clearly above every disc
  and ellipse on `chr_radial_hf` and `chr_concavities_per_um`, and clearly below them on `chr_solidity`, at each r.
  - At a fixed relative amplitude, a star's in-band outline energy barely depends on its lobe count, so do not expect
    `chr_radial_hf` to increase with it.
  - The 5-lobed star's fundamental lies below the band, so it scores low on `chr_radial_hf`.
- **The lobulation test.** Two touching discs and the notched disc are smooth outlines, not spikes. They stand in for
  "two cells in close proximity" (Lin 2024) and for the crowded-nuclei false positives of Section 1.
  - At each r, every spike column must score them below the 8- and 12-lobed stars and the spiked ellipse, and
    `chr_solidity` must score them above those three shapes.
  - All three pass. On `chr_solidity`, the lobulated shapes score 0.819–0.869 against at most 0.753 for the spiked
    shapes. Record `chr_solidity`'s confound: touching or lobulated nuclei score lower than a disc, in the mitotic
    direction.
  - A column that fails gets direction 0 and is described as a lobulation measure in `directions.csv`.
- **Scale.** For the uniformly scaled shapes (disc, ellipses, stars, touching discs, notched disc), the two scales must
  agree within 10 % of the larger value. Do not widen the tolerance.
  - **`chr_radial_hf` fails on every shape** (the closest is the touching discs at 12 %). Its low values are mostly
    pixel-staircase noise, which falls with r², and 8 or 12 lobes on an 8 px radius are 2–3 px wide, at the pixel
    pitch. Report those as failures caused by the resolution limit, with the radius below which they apply. This is
    the one predicted failure in this suite.
  - **`chr_solidity` passes** (largest difference 7.8 %, the 12-lobed star).
  - **`chr_concavities_per_um` is exempt.** It is a count per µm of outline, so with a fixed number of concavities it
    halves when r doubles, by construction.
  - **The spiked ellipse is exempt,** because its spikes keep a fixed pixel size and so shrink relative to the body as
    r grows, as real chromosome arms do on a larger figure. The consequence: fixed-size spikes score 3 × lower on
    `chr_radial_hf` on a body twice the size (1.346 → 0.453), so the column is negatively size-confounded on its own cue. V2's `auc_area_tertiles` was the real-data check; it moved
    no boundary column by more than 0.01, and V3 drops it.
- **Elongation.** `chr_elongation` separates discs from ellipses. Beyond the partial sensitivities noted, no boundary
  column does.

**End-to-end C1 test, mandatory.** The suite above hands the mask in directly. This test runs the verbatim refinement,
so it shows which arms survive into the masks the columns actually read.

- **Construction.**
  1. Start from a 201 × 201 canvas at 0.03 OD, with a cell disc (ρ ≤ 30 px) at 0.045 OD.
  2. Draw the shapes as constant values that replace what is below: the arms at the arm OD first, then the body at
     0.20 OD.
  3. Apply `scipy.ndimage.gaussian_filter` at σ.
  4. Add `np.random.default_rng(seed).normal(0, 0.003, (201, 201))`, then cast to float32.
  5. Run the verbatim `refine_template_box` at (100, 100), mpp 0.25.
  6. Build `M_chr` as in Section 4. Build `M_dark` as in Section 4 with `t_dark` **fixed at 0.088 OD**, which is
     0.44 × the body. Do not compute it from the canvas: a canvas-wide multi-Otsu moves with the drawn figure.
- **Shapes.**
  - **Spiked:** a body disc ρ ≤ r plus 8 rectangles at angles 2πk/8, covering r − 1 ≤ t ≤ r + L and |o| ≤ w/2, where
    t and o are the along- and across-ray coordinates. An arm is w + 1 px wide on axis-aligned rays and about w on
    diagonals.
  - **Touching discs:** radius 0.9r, centred at x = ±0.9r.
  - **Plain disc:** ρ ≤ r.
  - **Debris cluster:** a central disc of radius r_f, plus n_f − 1 more discs of radius r_f. Disc k sits at distance
    8 + 3·(k mod 2) px and angle 2πk/(n_f − 1) + 0.3, for (r_f, n_f) ∈ {(2.5, 8), (3.0, 8)}.
  - Every shape is drawn at 0.20 OD unless an arm OD is given.
- **The deciding cell** is the only cell that can move a direction: w = 4, arm OD 1.0 × body, σ = 1,
  r ∈ {8, 12}, L ∈ {4, 8}, seed 0.
  - **Pass:** at each (r, L), the spiked shape scores above the touching discs on `chr_radial_hf`,
    `chr_concavities_per_um`, `dark_concavities_per_um` and `chr_tophat_edge`, and below them on `chr_solidity`.
  - A column that fails gets direction 0, logged, as in the lobulation rule.
  - **Checked values, seed 0.** Every column passes every cell. The `M_chr` columns also pass on seeds 1–3, and
    `dark_concavities_per_um` on seed 1.

    | r | shape | L | `chr_radial_hf` | `chr_concavities_per_um` | `dark_concavities_per_um` | `chr_tophat_edge` | `chr_solidity` |
    |---|---|---|---|---|---|---|---|
    | 8 | touching | — | 0.142 | 0.097 | 0.091 | 0.015 | 0.824 |
    | 8 | spiked | 4 | 1.446 | 0.313 | 0.299 | 0.100 | 0.747 |
    | 8 | spiked | 8 | 3.791 | 0.189 | 0.183 | 0.136 | 0.576 |
    | 12 | touching | — | 0.136 | 0.062 | 0.060 | 0.008 | 0.870 |
    | 12 | spiked | 4 | 0.852 | 0.260 | 0.239 | 0.058 | 0.727 |
    | 12 | spiked | 8 | 2.406 | 0.163 | 0.157 | 0.101 | 0.526 |

- **The reporting cells** report resolvability and never move a direction:
  - arm OD 0.7× and 0.5× the body;
  - w = 2;
  - σ = 1.5;
  - the plain disc;
  - the two debris clusters.

  For each, report the share of arm pixels inside `M_chr` and every outline and top-hat column.
  - **Expected concavities, `M_chr` → `M_dark`** (seed 0, σ 1, w 4):

    | r | L | arms at 0.7× | arms at 0.5× |
    |---|---|---|---|
    | 8 | 4 | 0.075 → 0.323 | 0.000 → 0.402 |
    | 8 | 8 | 0.075 → 0.196 | 0.000 → 0.218 |
    | 12 | 4 | 0.000 → 0.256 | 0.000 → 0.303 |
    | 12 | 8 | 0.052 → 0.168 | 0.000 → 0.190 |

  - **Expected debris clusters** (seed 0, σ 1):
    - r_f 2.5: `chr_concavities_per_um` NaN (largest piece < 30 px), `dark_concavities_per_um` 0.227,
      `chr_tophat_edge` 0.330;
    - r_f 3.0: 0.237, 0.203 and 0.225.
  - **A loss of arms from `M_chr` in these cells is a documented property,** "the outline of the dense core", not a
    defect to fix.

### 5.3 Darkness and the surround profile

| column | definition | unit |
|---|---|---|
| `h_chr_mean` | mean H over `M_chr`, in float64; the earlier `core_od` | OD |
| `h_ring1_mean`, `h_ring2_mean`, `h_ring3_mean` (role `diagnostic`) | mean H over `ring_1`, `ring_2`, `ring_3` | OD |

Expect `h_chr_mean` and `od` to be nearly collinear (rank correlation 0.94 on the earlier pool); `od` stays as the
identity column and one of Section 6's two stratifiers.

The three ring means are stored as diagnostics, not features. `ring_peak_prominence` (5.5) is a function of them and
`h_chr_mean`, and it is the one surround-profile column tied to a hypothesis (C2). For the same reason, the
chromatin-to-surround contrast `h_chr_mean − h_ring3_mean` and the ring slope `(h_ring3_mean − h_ring1_mean) / 2` are
not stored.

### 5.4 Chromatin texture

**The texture support:** `M_chr` eroded once by `scipy.ndimage.binary_erosion` with its default cross, minus every
pixel whose 3 × 3 neighbourhood contains a clipped pixel (Section 4). `h_chr_lap` is NaN when the support has fewer
than `NS_MIN_PX` = 10 pixels. The reason is `small support` when the eroded `M_chr` alone has fewer than 10 pixels (a small or thin mask), and
`clipped` when the eroded `M_chr` has at least 10 and the clipped-pixel exclusion takes it below.

| column | definition | unit |
|---|---|---|
| `h_chr_lap` | std (ddof 0, in float64), over the texture support, of `cv2.Laplacian(H_w, cv2.CV_32F, ksize=1) / mpp**2`, divided by `h_chr_mean`. With ksize 1 the aperture reads only the 4-neighbours, which lie in `M_chr` for every support pixel, so no pixel outside the mask enters | 1/µm² |

**Direction 0** (Section 5.0, texture). A GLCM contrast column was dropped in V2 because it reads the blurred edge ramp
inside `M_chr`, and so measures size. Through the verbatim refinement, a flat blurred disc scored 10.0 / 5.1 / 2.4 / 1.9
at r = 5 / 8 / 12 / 16. Computed on the eroded support instead, it duplicated `h_chr_lap` (Spearman 0.90).
`h_chr_lap` depends on size too, more weakly (1.74 → 1.23 over the same radii). On the V2 pool its x-1 held within
area tertiles (V2's `auc_area_tertiles`, 0.397). V3 drops that reading; 5.9 gives the area-quintile value.

**Noise test, mandatory.** On a 71 × 71 canvas at `mpp` 0.25, render a disc of radius 16 px at 0.20 OD on 0.03 OD
(`M_chr` = the disc) in four versions: flat; flat plus Gaussian noise of σ 0.003 everywhere; the noisy one blurred at
σ 1 px; and the granular disc of Section 5.8's reference construction (a 2 px checkerboard of +0.06 OD more than 3 px
inside the edge). `h_chr_lap` must score the flat, noisy and noisy-blurred discs far below the checkerboard; if it
does not, it is measuring noise or the blurred edge. Checked values, one implementation (flat / noisy / noisy-blurred /
checkerboard): 0.00 / 1.07 / 1.03 / 7.72.

**Clipped-pixel test, mandatory.** Render a uniform RGB disc of (40, 1, 100) on a (230, 200, 220) background, then set
G = 0 on a random 7 % of the disc's pixels (canvas 71 × 71, radius 16 px, `M_chr` = the disc). `h_chr_lap` must stay
at its value on the unaltered disc, because the clipped pixels and their neighbourhoods are excluded. Checked values,
one implementation: unaltered 0; with the clipped pixels excluded 0; with them included 3.89. The exclusion costs support: 709 px fell to 383, because scattered
clipped pixels take their 3 × 3 neighbourhoods with them. On real data clipped pixels cluster in the densest core.

### 5.5 Surround: membrane, clear zone, white space, cytoplasm

| column | definition | unit |
|---|---|---|
| `ring_peak_prominence` | for k in {2, 3}: `h_ring_k_mean` minus the minimum of `h_ring_j_mean` over j < k; take the larger, divide by `h_chr_mean`, floor at 0. A rim that is darker than the band inside it scores above 0; a profile that only falls scores 0 | ratio |
| `ring1_frac_above_tnuc` | share of `ring_1` pixels above `t_nuc`: nucleus or dense cytoplasm, since `t_nuc` admits cytoplasm (Section 4) | share |
| `ring1_frac_pale` | share of `ring_1` pixels at or below `t_nuc` and with `gray ≤ 215` (the `WHITE_GRAY` of the floor-clip run): pale tissue, cleared cytoplasm, halos. The remaining share, glass and lumen (at or below `t_nuc` with `gray > 215`), is 1 − `ring1_frac_above_tnuc` − `ring1_frac_pale` and is not stored | share |
| `white_frac_gray215_61` | share of the 61 px window with `gray > 215`, as `white_frac` in the floor-clip run computes it; the earlier `not_white` without the sign | share |
| `e_h_angle_ring23_deg` | `np.degrees(np.arctan2(mean E, mean H))`, with E = `hed[:, :, 1]` and H = `hed[:, :, 0]`. The means are taken over the cytoplasm ring: `ring_2` ∪ `ring_3`, minus glass (`gray > 215`), minus pixels above `t_dark`, minus clipped pixels. NaN, reason `no cytoplasm`, when the ring has fewer than `NS_MIN_PX` pixels or when mean E and mean H are both 0. rgb2hed clamps stains at 0, so the angle lies in [0°, 90°]; higher is more eosinophilic | degrees |

**Directions.** `e_h_angle_ring23_deg` has direction −1 (C9). Its confounds:
- "color is not reliable" (VCGP 2021);
- a dark clump inside an interphase nucleus has nucleoplasm in its 1–3 µm ring, so it reads as basophilic, the
  mitosis-like direction;
- the fixed `rgb2hed` vectors make pooled values incomparable across scanners (5.8 a).

The angle ranks candidates exactly as the E/H ratio would, without dividing by a mean H that can approach 0 on
eosinophilic cytoplasm. Every other column here has direction 0 (below, and C9).

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
  (C5), which is the job of Sections 5.2 and 5.8 (b). A membrane column scoring them like mitoses is expected, not a
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
`ring1_frac_pale` separates pale tissue from glass and lumen. The one contrast the ring shares are expected to carry: a
dark clump inside an interphase nucleus is surrounded by nucleoplasm above `t_nuc`, so its `ring1_frac_above_tnuc`
should be high.

### 5.6 Colour and pigment

RGB statistics are on the raw `uint8` image.

| column | definition | unit |
|---|---|---|
| `brown_frac_surround` | the number of 51 px window pixels outside `M_sel` with `R > G > B`, `R − B > 40`, `G − B > 15` and `gray < 180`, divided by the number of window pixels outside `M_sel` | share |

There is no chromatin-colour column. No source gives a chromatin colour criterion, and pigment inside the mask is
`chr_resid_frac`'s job (5.8 a).

**Why the brown rule has four terms.** Without `G − B > 15`, the rule fires on red blood cells and pink stroma as often
as on pigment.

| ROI | what fires | share of ROI pixels, without the term | with it |
|---|---|---|---|
| `123` | red blood cells in a vessel | 0.280 % | 0.026 % |
| `128` | scattered pixels in pink stroma | 0.416 % | 0.041 % |
| `291` | golden-brown pigment granules | 0.417 % | 0.314 % |
| `410` | golden-brown pigment granules | 0.325 % | 0.189 % |

- Over the other ROIs, the maximum with the term is 0.041 % (`128`).
- The term was chosen after inspecting these ROIs, so its specificity is calibrated on the same data.
- It drops 42 % of `410`'s pigment pixels.
- The rule is a heuristic. If the top tiles of its montage are not pigment, you may change only its four
  thresholds. Log the change, and keep the original as `spec_brown_frac_surround`.

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
`brown_frac_surround` are the visual check; the V2 notebook holds them, and this run draws neither. The V2 run's pool values: 0.244, 0.326 and 0.353.

### 5.7 Context and density

| column | definition | unit |
|---|---|---|
| `dark_area_frac_101` | among the (truncated) 101 px window's pixels outside the candidate's own `M_sel` (placed as in Section 4), the share above `t_dark` | share |

Direction 0 (C10). The candidate's own selection is left out because it is 2.1–8.4 % of the window, which would put an
object-size term into a crowding column. Store `t_nuc` and `t_dark` per ROI in `clicks.csv`.

### 5.8 In-mask pigment, thin protrusions, edge sharpness

These three families target the dark false positives of the failing ROIs (Section 1). Edge sharpness (c) is stored as a
diagnostic, not a feature. Every rule at the top of Section 5 applies. Every constant named here goes into the spec
JSON.

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
- **Its relation to rgb2hed.** `n` equals the normalised third column of `skimage.color.hed_from_rgb` (checked
  numerically). So the column is rgb2hed's *unclamped* D of the mean OD, divided by that column's norm and by the OD's
  norm. This is why it works where the earlier `d_over_h` failed: the clamp, the support and the denominator all
  changed.
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
    arm. Whether a scan resolves an arm that thin is a separate question (5.8 c, and V2's per-scanner columns).
- **Regions, following the mask** (`M_chr` *whole*):
  - `R` = the window pixels with `distance_transform_edt(~M_chr) * mpp <= TOPHAT_RING_UM`, with `TOPHAT_RING_UM` = 1.0,
    so that fuzzy protrusions up to 1 µm outside the mask count. `M_chr` itself is included.
  - `M_core` = the pixels with `distance_transform_edt(M_chr) * mpp > TOPHAT_CORE_UM`, with `TOPHAT_CORE_UM` = 0.5:
    chromatin more than 0.5 µm inside its own edge.
  - `R_edge` = `R` ∧ ¬`M_core`.
- **The features:**
  - `chr_tophat_edge = T[R_edge].sum() / H_w[M_chr].sum()`: thin dark structure within 1 µm of the core's surface,
    as a share of the core's total H. Fine arms, debris fragments and thin rims all count;
  - `chr_tophat_core = T[M_core].sum() / H_w[M_chr].sum()`: thin structure inside, which granular interphase chromatin
    also has;
  - both NaN, reason `small piece`, when *whole* `M_chr` has fewer than 30 px. This is the whole-mask rule of
    Section 5.2, not the largest-piece rule.
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
    question V2's per-scanner columns ask;
  - the granular disc scores clearly above the disc on `chr_tophat_core`, and no more than 0.01 above it on
    `chr_tophat_edge`;
  - the tip-spiked bar scores clearly above the plain bar on `chr_tophat_edge` at both scales;
  - **a known confound, not a pass:** the vesicular window of Section 5.5, refined verbatim, scores `chr_tophat_edge`
    0.34, above the 2 px spikes (0.256), because its `M_chr` is a thin continuous rim. A top-hat cannot tell a smooth
    thin membrane from discrete spikes. Record the value. The effect cannot be read on real data: only 49 accepted
    pool rows have `chr_hole_frac` above 0.3 (228 at V2's depth 300), and no click has both classes there.
- **The scale rule does not apply.** These columns are exempt from Section 5.2's scale rule, because their width
  threshold is physical by design. At r = 8 the plain bar is 5 px wide, exactly the structuring element's width, so the
  opening removes only its tip corners (12 of its 85 px); a bar thinner than 5 px would be removed entirely.
- **Small bodies score high too.** The whole of a body narrower than about 1.25 µm counts as thin structure. So
  karyorrhectic fragments and other tiny round bodies, which are imposters (C5, C7), raise `chr_tophat_edge` as spikes
  do. Write the confound into `directions.csv`.
- **What the end-to-end test adds.** `chr_tophat_edge` is the only boundary column that reads the image beyond `M_chr`.
  When arms are paler than the body, it keeps a reduced signal (25–60 % of its value at arm OD 0.7×), where the
  outline columns on `M_chr` lose them.
- **Direction hypotheses:** `chr_tophat_edge` +1, with the small-body and vesicular-rim confounds above. Add to its
  `confound`: "on the recorded A3 pool its within-click AUC is below 0.5, but about 0.5 within od strata (seen before
  freezing, Section 2)". `chr_tophat_core` gets 0 (Section 5.0, texture).
- **Measured on the V2 pool (5.9).** `chr_tophat_edge` is null on single-piece `M_chr` (0.501); its signal is `M_chr`
  fragmentation, which `chr_frag` measures directly. It is not a thin-arm measure.

**(c) `chr_edge_sharpness_per_um` (role `diagnostic`): how crisp the chromatin boundary is.**
- **The boundary `B`:** the hole-filled largest 4-connected piece of `M_chr` (Section 5.2), minus its `binary_erosion`
  with a 3 × 3 cross. NaN, reason `small piece`, below 30 px.
- **The gradient:** `gx = cv2.Sobel(H_w, cv2.CV_32F, 1, 0, ksize=3) / 8` and `gy` likewise; `g = np.hypot(gx, gy) / mpp`,
  in OD per µm.
- **The step:** `step = h_chr_mean − h_ring1_mean`, from Section 5.3. If `step` < `EDGE_MIN_STEP` = 0.005 OD, the
  feature is NaN with reason `no edge step`.
- **Which reason column is filled.** Check in this order and stop at the first that applies: `small piece` goes in
  `nan_reason_geom`; `empty ring` (`ring_1` empty, so `step` is undefined) goes in `nan_reason_ring`; `no edge step`
  goes in `nan_reason_edge`. The other two stay empty.
- **The feature:** `chr_edge_sharpness_per_um = g[B].mean() / step`, in 1/µm. It is roughly the inverse width of the
  edge.
- **Synthetic suite, mandatory.** A disc and a 2:1 ellipse, at radius 8 and 16 px, 0.20 OD on 0.03, each
  Gaussian-blurred at σ = 0, 1 and 2 px:
  - the value must fall strictly with σ;
  - the ellipse must score within 15 % of the disc at the same σ;
  - the two radii must score within 15 % of each other.
- **What it is for.** It is the covariate for whether boundary features are resolvable on a given scan. In MIDOG
  2022, 27.5 % of the mitoses found only with PHH3 help were classed "out of focus" (Aubreville et al. 2024, Table 7).
  No source names blurred dark objects as a false-positive class, and no counting criterion stands behind the column,
  so it is a diagnostic.
  - V2's validity tables report each scanner group's median beside the per-scanner AUCs. V3 does not recompute them.
  - It stays in the montage list, because the montage checks that it measures edge crispness, whatever its role.
- **Its `directions.csv` row:** role `diagnostic`, direction 0.

**(d) Dropped in V3: the neighbour columns.**
- **What is gone:** `chr_area_over_nbr`, `nbr_n`, `nbr_nuc_area_med_um2`, and the V2 run's `spec_` versions at `t_nuc` (user, 2026-10-09).
- **Why:**
  - On the V2 pool they were NaN on 89 % of accepted rows.
  - Where defined, `chr_area_over_nbr` tracked `chr_area_um2`: within-click ρ 0.78, AUC 0.457 against 0.460.
  - Within a click the neighbour median barely varies, so a within-click rank of `chr_area_um2` already gives the normalisation.
- **The old text:** V2's Section 5.8 (d).

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

### 5.9 What the columns measured on the V2 pool (verified 2026-10-09)

These are the `measures` labels of Section 6. Each statement comes from the fix plan's measurements on this pool,
except the two parts marked *by eye*: those rest on V2's and the audit's montage readings, with no measurement behind
them (Section 1's rule). They change no `direction_hypothesis`. A column not listed keeps V2's `claims_to_measure`.

| column | measures |
|---|---|
| `chr_tophat_edge` | `M_chr` fragmentation, noisily. Null on single-piece `M_chr` (0.501 [0.466, 0.535]); its pool AUC comes from fragmented rows. Not a thin-arm measure. `chr_frag` measures the same thing directly |
| `chr_frag` | whether the window's top Otsu class inside the selection breaks into separate 4-connected pieces. The montage judgement goes in the findings document, not in this label |
| `h_chr_lap` | relative chromatin texture (Laplacian SD over core mean). Its x-1 holds within area quintiles (0.395), on single-piece `M_chr` (0.410) and at matched `h_chr_mean`. The numerator alone (`h_chr_lap` × `h_chr_mean`) is +1 marginally (0.541 [0.497, 0.585], 33 of 49 ROIs) and x-1 at matched `h_chr_mean` (0.414): mitoses have more absolute texture because they are darker, and less relative to their darkness |
| `chr_area_um2` | core area. x-1 marginally. +1 at matched core darkness: `h_chr_mean` 0.604 (`raw`). Two mask-free darkness measures point the same way, the window's darkest 50 px (0.566 [0.534, 0.599]) and a 1 µm disc at the H maximum (0.544 [0.510, 0.577]), but both are `non-monotone` by the ROI-share rule (0.65 and 0.61). No consistent relation at matched `od` (0.458), because `od` grows with core size on cores under 260 px (within-click ρ +0.10) |
| `chr_hole_frac` | holes in the largest piece. x-1 marginally; no consistent relation at matched darkness |
| `chr_radial_hf` | outline high-harmonic energy. `raw`, also at matched darkness. Overlaps `chr_solidity` (ρ −0.67) |
| `dark_area_frac_101` | dark area around the candidate (crowding). No marginal orientation; x-1 at matched darkness under both covariates. Its `od`-matched value is inflated: `od` counts neighbour pixels (within-click ρ with this column 0.23, against 0.11 for `h_chr_mean`) |
| `ring_peak_prominence` | a band 1–3 µm out that is darker than the band inside it (its definition). *By eye:* mostly a dark neighbour, or a darker nucleus around a clump, not an envelope |
| `chr_resid_frac` | chromatin colour off the H–E plane. Pigment is its high tail, and it is x-1 on non-pigment ROIs too. Not explained by the slide's glass level, the 2 px rim, or its support |
| `e_h_angle_ring23_deg` | the eosin of the 1–3 µm ring; *by eye*, red cells and stroma in the neighbourhood. No orientation. Its per-ROI spread (SD 0.148) is not distinctive: `h_chr_lap`, `pair_area_ratio` and `ring1_frac_pale` read 0.146–0.148 |
| `brown_frac_surround` | brown pixels beside the candidate (97.7 % of rows are 0). A pigment-ROI flag with almost no ranking signal even on its own ROIs: 291 reads 0.490, 410 reads 0.440 |

**More on `chr_frag`:**
- 72 % of pool rows have `chr_frag` = 0.
- Mitoses are fragmented on 18 % of rows, against 31 % of others.
- It is not the refinement unioning pieces: 86 % of fragmented rows have `n_pieces` = 1, and `n_pieces` itself carries no signal (AUC 0.491).
- V2 had dropped `chr_n_pieces`, which was a different count.
- The fix plan (3.2) recommends that the downstream carry `chr_frag` and not `chr_tophat_edge`, its noisier copy. This
  spec keeps both as measured columns and selects nothing (Section 0); the recommendation goes in the findings.
- Its unit test (5.2) was run on 2026-10-09 against the audit's own formula and gives 0, 0.338 and 0.5.

## 6. The validity report (descriptive; no selection, no verdict)

**Population.** The accepted rows of this run's `candidates.csv`: 28,883 rows from 294 clicks, 7,664 of them mitoses.

**Rows.** The 24 features of Section 5, in Section 5's order, with `chr_frag` after `chr_hole_frac`.

**The AUC.** Every AUC in this section is computed the same way, unchanged from V2. For a set of groups (clicks, or
strata within clicks):
1. In each group, drop the rows where the column is NaN, then keep the group only if it still has at least one mitosis
   (`is_mitosis`) and one other row.
2. In each kept group, U = the number of (mitosis, other) pairs in which the mitosis has the larger value, plus 0.5 for
   each tie.
3. AUC = ΣU / Σ(n1·n0), summed over the kept groups, where n1 and n0 are each group's mitosis and other counts.

The AUC is the raw value of P(mitosis > other): never folded about 0.5, so a value below 0.5 means the column is lower
on mitoses.

### 6.1 Three readings per feature

**The readings:**
- **`wc`:** the groups are the clicks.
- **`od`:** the groups are click × within-click `od` quintile.
- **`hmean`:** the groups are click × within-click `h_chr_mean` quintile.

**How the strata are formed:**
- Once per click, on the whole population, before any column's NaN rows are dropped.
- With `pd.qcut(v.rank(method='first'), 5, labels=False)`.
- Neither covariate has a tie within any click (checked), so the strata are exact.
- **Stratify on nothing else.** In particular, never stratify on a tie-heavy column with `rank(method='first')`. On rows in `cand_rank` order that turns the tied block into score strata, which is what produced the audit's 0.524.

**Why these two covariates.**
- They are the two darkness signals of the best existing rule: the 2026-09-30 rank average of `od`, core OD mean, not-white and solidity.
- They are tied as pool rankers: 0.770 [0.739, 0.799] against 0.761 [0.729, 0.794].
- **`od` holds dark mass, not core darkness.** Its darkest 260 px include surround on cores under 260 px (36 % of pool rows), and they include dark neighbours.
- So a beyond-darkness orientation that a rule can rely on has to hold under both.

**The strata's resolution.**
- Quintiles of about 20 rows leave each covariate separating within its own strata: `od` reads 0.571 within `od` quintiles, `h_chr_mean` 0.576 within its own. Print both.
- A darkness column's conditional AUC is read against that floor, not against 0.5.
- Caliper matching was tested. It moves no non-darkness column by more than 0.015. Do not add it.

### 6.2 Columns of `validity.csv`

**Identity:**
- `column`;
- `n_clicks`, `n1` and `n0` of the `wc` reading;
- `direction_hypothesis` and `claims_to_measure`, from V2's `directions.csv`. `chr_frag` has no V2 row: give it direction 0 (no hypothesis; it was found post hoc, which the findings say) and its 5.2 definition;
- `measures`: the 5.9 label, or `claims_to_measure` where 5.9 is silent.

**For each reading r ∈ {`wc`, `od`, `hmean`}:**
- **`auc_<r>`, `ci_lo_<r>`, `ci_hi_<r>`:** the ROI cluster bootstrap, as follows.
  - `rois` = the sorted list of the 49 `file_name`s.
  - `BOOT_W = np.random.default_rng(20261009).multinomial(49, np.full(49, 1 / 49), size=4000).astype(float)`, drawn once. The same `BOOT_W` serves every cell and the concavity lead (6.3); do not draw per column. Never name it `W`: `W` is the imported V2 worker, which the montage step reads after this one.
  - For each ROI, sum U and N = n1·n0 over its kept groups; an ROI with no kept group gets 0.
  - Bootstrap AUCs = `(BOOT_W @ U) / (BOOT_W @ N)`. The CI is `np.nanpercentile` at 2.5 and 97.5.
- **`roi_share_<r>`:** the share of ROIs with N > 0 whose AUC (their U / N) lies on the same side of 0.5 as `auc_<r>`. An ROI at exactly 0.5 does not agree.
- **`domains_<r>`:** the number of the 7 domains whose AUC (U and N summed over the domain's ROIs) lies on that side.
- **`orientation_<r>`:** the V2 rule unchanged.
  - `raw` when `auc_<r>` > 0.5, `roi_share_<r>` ≥ 2/3 and `domains_<r>` ≥ 5;
  - `x-1` when `auc_<r>` < 0.5 with the same two conditions;
  - otherwise `non-monotone`, including exactly 0.5.
  - The CI is reported and is not part of the rule. Across all 72 cells, adding "the CI excludes 0.5" changes no label.

**Plus:**
- **`beyond_darkness`:**
  - the shared label when `orientation_od` and `orientation_hmean` agree;
  - `entangled` when one is `raw` and the other `x-1`;
  - otherwise `od: <label> / h_chr_mean: <label>`.
- **`auc_wc_failing9`, `n1_failing9`, `auc_wc_other40`, `n1_other40`:** the `wc` AUC within each `roi_group`, as in V2.
- **Top 10 against bottom 10 per click,** oriented by `orientation_wc`; for `non-monotone`, by the side of `auc_wc`.
  - Extract V2's `top_bottom` and `tied_expectation` from cell 27 of `candidate_features_a3_top300.ipynb` with `ast`, as the V2 worker extracts code. Run them with `TB_K = 10` and `TB_MIN_ROWS = 20`, the population given the `click` key `file_name + '|' + seed_index.astype(str)`.
  - Store their ten outputs. The tie-averaged pair `top10_mitoses_tieavg` / `bottom10_mitoses_tieavg` is the one to report.

### 6.3 The other outputs of the report

**`redundancy.csv`.** The V2 definition, over the 24 features and `od`, on the same population.
- **The statistic:** V2's cell 29 code. Per click, `g[RED_COLS].astype(float).corr(method='spearman', min_periods=10)`, which ranks each pair on the rows where both are non-NaN; then the mean over clicks of each entry, skipping NaN (a pair with fewer than 10 such rows, or a column constant in the click). Do not rank whole columns first: that changes every pair with a NaN-bearing column and fails the 1e-9 check.
- **What to print:** the pairs above 0.8 in absolute value, and `chr_frag`'s three largest. Expected for `chr_frag`: `chr_tophat_edge` 0.35, `chr_solidity` −0.26, `chr_radial_hf` 0.25.

**The two concavity columns, on uncut rows.** Section 5.2 allows comparing `chr_concavities_per_um` with `dark_concavities_per_um` only on rows where `dark_touches_window` is False.
- **The rows:** accepted rows with `dark_touches_window == False` (23,836). The column loads as object dtype, so compare with `== False`, not `~`. Drop NaN per column, as every AUC does.
- Report both `wc` AUCs on those rows.
- Report the lead (dark minus chr), with its CI as the percentiles of the difference of the two columns' bootstrap AUCs under the same `BOOT_W`.

### 6.4 Expected values

The table below was computed independently on 2026-10-09. Its 72 AUC and CI cells also equal the audit's committed
`results/candidate_features_a3_top300_audit_stratified.csv` (rows `population == 'pool'`, controls `none`, `od` and
`h_chr_mean`), whose bootstrap is this `BOOT_W`; the uncut values are in `results/candidate_features_a3_top300_audit_checks.json`
(`uncut_pool_dark_minus_chr_concavities`).

**Your validity table must reproduce it to the printed decimals:**
- every AUC, CI, ROI share, domain count and label;
- `od` 0.571 and `h_chr_mean` 0.576 within their own quintiles;
- on uncut rows: `chr_concavities_per_um` 0.636, `dark_concavities_per_um` 0.656, lead +0.020 [+0.002, +0.039].

**These must equal V2's committed values to 1e-9,** on the 23 columns carried over from V2:
- in V2's `validity.csv`:
  - `auc_within_click` (= `auc_wc`);
  - `auc_dark_stratified` (= `auc_od`);
  - `roi_sign_agree_share` (= `roi_share_wc`);
  - `domain_sign_agree_n` (= `domains_wc`);
  - `orientation_recommended` (= `orientation_wc`);
  - `auc_wc_failing9` and `auc_wc_other40`;
  - all ten of `top_bottom`'s outputs;
- in V2's `redundancy.csv`: every entry between those columns and `od`.

Checked on 2026-10-09: the verification code reproduces every one of these exactly, using V2's own `top_bottom` (`auc_od` against `auc_dark_stratified`: largest difference 6e-17).

**`chr_frag`'s top 10 against bottom 10:** expect 851.5 / 585.0 mitoses with ties averaged. The `cand_rank` tie-break
gives 1,660 / 566: a tie straddles the top-10 cut on every click, because 72 % of rows are 0. Report the
tie-averaged pair.

**How to read the table:** each cell is AUC [95 % ROI CI], then ROI share / domains, then the cell's label.

| column | `wc` | `orientation_wc` | `od` | `orientation_od` | `hmean` | `orientation_hmean` | `beyond_darkness` |
|---|---|---|---|---|---|---|---|
| `chr_area_um2` | 0.460 [0.427, 0.493] 0.67/5 | x-1 | 0.458 [0.422, 0.494] 0.61/5 | non-monotone | 0.604 [0.574, 0.634] 0.84/7 | raw | od: non-monotone / h_chr_mean: raw |
| `chr_elongation` | 0.636 [0.599, 0.671] 0.78/7 | raw | 0.627 [0.588, 0.665] 0.76/7 | raw | 0.601 [0.561, 0.640] 0.67/6 | raw | raw |
| `chr_solidity` | 0.375 [0.340, 0.411] 0.88/7 | x-1 | 0.361 [0.326, 0.396] 0.86/7 | x-1 | 0.351 [0.318, 0.385] 0.94/7 | x-1 | x-1 |
| `chr_concavities_per_um` | 0.621 [0.595, 0.646] 0.90/7 | raw | 0.642 [0.619, 0.666] 0.96/7 | raw | 0.669 [0.647, 0.693] 0.98/7 | raw | raw |
| `dark_concavities_per_um` | 0.651 [0.623, 0.680] 0.94/7 | raw | 0.674 [0.646, 0.702] 1.00/7 | raw | 0.685 [0.658, 0.712] 0.98/7 | raw | raw |
| `chr_hole_frac` | 0.435 [0.409, 0.459] 0.67/6 | x-1 | 0.460 [0.434, 0.485] 0.61/6 | non-monotone | 0.493 [0.468, 0.516] 0.47/4 | non-monotone | non-monotone |
| `chr_frag` | 0.446 [0.426, 0.464] 0.82/7 | x-1 | 0.460 [0.440, 0.481] 0.73/6 | x-1 | 0.469 [0.449, 0.489] 0.69/5 | x-1 | x-1 |
| `chr_radial_hf` | 0.561 [0.530, 0.591] 0.82/6 | raw | 0.593 [0.562, 0.624] 0.84/7 | raw | 0.625 [0.598, 0.652] 0.92/7 | raw | raw |
| `pair_sep_um` | 0.581 [0.540, 0.619] 0.71/6 | raw | 0.624 [0.589, 0.658] 0.83/7 | raw | 0.624 [0.591, 0.656] 0.81/7 | raw | raw |
| `pair_area_ratio` | 0.542 [0.517, 0.564] 0.65/7 | non-monotone | 0.538 [0.512, 0.565] 0.56/7 | non-monotone | 0.526 [0.493, 0.555] 0.56/6 | non-monotone | non-monotone |
| `pair_parallel_deg` | 0.517 [0.489, 0.544] 0.46/4 | non-monotone | 0.532 [0.502, 0.558] 0.54/3 | non-monotone | 0.508 [0.472, 0.540] 0.50/3 | non-monotone | non-monotone |
| `pair_axis_deg` | 0.531 [0.508, 0.556] 0.69/6 | raw | 0.514 [0.481, 0.547] 0.65/5 | non-monotone | 0.537 [0.505, 0.569] 0.62/6 | non-monotone | non-monotone |
| `h_chr_mean` | 0.761 [0.729, 0.794] 0.96/7 | raw | 0.604 [0.569, 0.638] 0.78/7 | raw | 0.576 [0.557, 0.596] 0.92/7 | raw | raw |
| `h_chr_lap` | 0.397 [0.355, 0.439] 0.71/5 | x-1 | 0.422 [0.382, 0.460] 0.63/6 | non-monotone | 0.395 [0.357, 0.432] 0.71/6 | x-1 | od: non-monotone / h_chr_mean: x-1 |
| `ring_peak_prominence` | 0.418 [0.394, 0.445] 0.78/6 | x-1 | 0.405 [0.381, 0.432] 0.84/7 | x-1 | 0.417 [0.392, 0.446] 0.82/6 | x-1 | x-1 |
| `ring1_frac_above_tnuc` | 0.634 [0.594, 0.677] 0.80/7 | raw | 0.594 [0.552, 0.639] 0.76/6 | raw | 0.609 [0.565, 0.657] 0.80/7 | raw | raw |
| `ring1_frac_pale` | 0.422 [0.377, 0.466] 0.63/6 | non-monotone | 0.461 [0.416, 0.504] 0.55/5 | non-monotone | 0.446 [0.401, 0.491] 0.57/5 | non-monotone | non-monotone |
| `white_frac_gray215_61` | 0.363 [0.318, 0.407] 0.82/7 | x-1 | 0.347 [0.304, 0.390] 0.84/7 | x-1 | 0.354 [0.310, 0.401] 0.80/7 | x-1 | x-1 |
| `e_h_angle_ring23_deg` | 0.504 [0.463, 0.547] 0.49/3 | non-monotone | 0.509 [0.471, 0.549] 0.51/3 | non-monotone | 0.487 [0.448, 0.530] 0.57/4 | non-monotone | non-monotone |
| `brown_frac_surround` | 0.495 [0.490, 0.499] 0.29/6 | non-monotone | 0.495 [0.491, 0.498] 0.29/6 | non-monotone | 0.494 [0.489, 0.498] 0.27/6 | non-monotone | non-monotone |
| `dark_area_frac_101` | 0.466 [0.434, 0.496] 0.65/5 | non-monotone | 0.407 [0.374, 0.437] 0.86/7 | x-1 | 0.436 [0.403, 0.469] 0.73/7 | x-1 | x-1 |
| `chr_resid_frac` | 0.373 [0.336, 0.410] 0.80/7 | x-1 | 0.393 [0.362, 0.424] 0.78/7 | x-1 | 0.399 [0.370, 0.431] 0.78/7 | x-1 | x-1 |
| `chr_tophat_edge` | 0.460 [0.422, 0.497] 0.57/4 | non-monotone | 0.517 [0.478, 0.556] 0.57/5 | non-monotone | 0.472 [0.434, 0.509] 0.47/5 | non-monotone | non-monotone |
| `chr_tophat_core` | 0.481 [0.444, 0.517] 0.47/4 | non-monotone | 0.520 [0.482, 0.554] 0.57/6 | non-monotone | 0.507 [0.470, 0.542] 0.57/3 | non-monotone | non-monotone |

The `hmean` reading of `h_chr_mean` is the column within its own strata: the floor of Section 6.1. So its
`beyond_darkness` label (`raw`) is not a beyond-darkness finding; say so wherever it is reported.

**Reading it:**
- **Robust: the same direction at matched `od` and at matched core darkness.**
  - `raw`: `chr_concavities_per_um`, `dark_concavities_per_um`, `chr_elongation`, `chr_radial_hf`, `pair_sep_um`, `ring1_frac_above_tnuc`;
  - `x-1`: `chr_solidity`, `chr_frag`, `ring_peak_prominence`, `white_frac_gray215_61`, `chr_resid_frac`.
- **Depends on the covariate:** `chr_area_um2` and `h_chr_lap`.
- **Conditional only:** `dark_area_frac_101`.

**What you may conclude from this report:** whether each measurement is valid, what it actually measures, and which
orientation puts mitoses first (in-sample). **What you may not conclude:** which features to use, how to weight them,
or whether any ranker is better. Write the report in that voice.

## 7. Gates (halt and report) and flags (report only)

**This run has no gates** beyond the SHA-1 assert of Section 2.1. Nothing gated is recomputed:
- **The rows:** gate-checked by V2. Gates A, B and C held, with 0 of 17,735 shared locations mismatched in gate C.
- **Their numbers:** reproduced by the audit, with every AUC exact.
- **`chr_frag`:** computed by the audit on the run's own `M_chr`. Its top-hat check reproduced V2's column exactly from that mask.

The two other checks of Section 2.1 test this run's new code.

**A future full run (not executed in this run):**
- **A. The refinement.** The SHA-1 of the extracted `refine_template_box` segment + `"\n"` is
  `21c43a20db00dccf7ac3bab88aee423c345441b6`. Halt otherwise.
- **B. On these 49 ROIs, the pool reproduces the recorded A3 pool.**
  - **Reference:** the `arm == 'A3'` rows of `results/tm_floor_definition_49roi_6seed_pool.csv.gz`, by their `rank`.
  - **Exact:** `n_detections`, `cx`, `cy`, rank order, `bucket` and `matched_ann_id`.
  - **Within 1e-4:** `score`, `od` and `threshold`.
  - **Near-ties are reported, not halted on.** These are the cases a float shift below 1e-4 explains:
    - two scores within 1e-4 swapping order, with their labels;
    - a peak within 1e-4 of the 100th peak's score entering or leaving;
    - a peak moving to an adjacent pixel, with the score within 1e-4 and identical labels. The V2 run met one of these.
  - Anything else halts.
- **C. On these 49 ROIs, the candidate measurements reproduce this run's table.**
  - Compare on exactly shared `(file_name, cx, cy)`: `refusal` must be equal, and on accepted rows `h_chr_mean` must agree within 1e-4.
  - Halt above 1 % mismatches.
- **Flags:**
  - refusal counts by class and bucket;
  - the NaN share by reason and bucket for every column;
  - the `chr_touches_window` and `dark_touches_window` shares;
  - the `chr_resid_frac` NaN share and the median `chr_ns_share`;
  - per click, `n_peaks` and the search seconds.

## 8. Outputs and files

**Notebook:** `production_hematoxylin_only/multiotsu_pad_arms/candidate_features_a3_pool100.ipynb`.

**Tables,** in `results/`, prefix `candidate_features_a3_pool100_`:
- **`candidates.csv`: 29,025 rows,** one per (click, candidate) of the 100-peak pool. These are V2's rows with `in_recorded_pool`, in V2's order.
  - **Columns:** every V2 column except:
    - `in_recorded_pool`, which is now constant;
    - the eight neighbour columns: `chr_area_over_nbr`, `nbr_n`, `nbr_nuc_area_med_um2`, `spec_chr_area_over_nbr`, `spec_nbr_n`, `spec_nbr_nuc_area_med_um2`, `nan_reason_nbr` and `nan_reason_spec_nbr`.
  - **Plus `chr_frag`,** placed after `chr_hole_frac`.
  - **Assert** 29,025 rows, and on every click a row count equal to `n_detections_recorded`.
  - **NaN reasons:** V2's reason columns stay, minus the two neighbour ones. `chr_frag` is NaN only on refused rows, so it has no reason column.
- **`validity.csv`** (Section 6.2) and **`redundancy.csv`** (6.3).

**Not written by this run.** V2's committed versions stand:
- `clicks.csv`, which describes the same clicks;
- `directions.csv`, which is frozen;
- `synthetic_tests.csv`, `adversarial_log.csv` and the spec JSON;
- in V2's `validity.csv`, the columns this run does not recompute: `nan_rate_mitosis`, `nan_rate_other` and
  `synthetic_pass` for the 23 carried-over columns, all computed on the pool. `chr_frag` has no `synthetic_pass` (its
  unit test was run on 2026-10-09, outside any run, 5.9) and a NaN rate of 0 on accepted rows.

**The `chr_frag` montage,** embedded in the notebook only.
- **Where it lives:** do not write image files into the repo; figures live in notebook outputs, by the user's 2026-10-02 rule. Embed the mosaic as a JPEG with `IPython.display`, as V2 did (quality 88). For your own look, save a temporary PNG outside the repo and open it with the Read tool.
- **The orientation** is `orientation_wc` (expected `x-1`). The left ten tiles are the least fragmented, the end that should hold mitoses.
- **Rows 1–8:** one per ROI: `300`, `201`, `013`, `529`, `425`, `211`, `546`, `524`.
  - **The tiles** come from this run's `candidates.csv` (pool rows only, not V2's 300-deep table): the accepted rows, deduplicated by `(file_name, cx, cy)`. A location's border bucket is its majority bucket over its rows, ties to mitosis, then look-alike (V2's rule); its `chr_frag` is the same on every row.
  - **Picking:** greedily in oriented-value order, skipping any location within 7.5 µm of one already picked.
  - **Ties:** shuffle the deduplicated table once (`df.sample(frac=1, random_state=0)`), then, for each row of the montage, take its subset and sort it stably (`kind='stable'`) on the value alone. 72 % of rows are tied at 0, and V2's tie rule (by `cx`, then `cy`) would show the ROI's leftmost locations, not a sample.
  - **The 7.5 µm exclusion** is shared by the two ends of a row, as in V2: the low end skips locations near a high-end pick.
- **Row 9:** the 10 deduplicated locations with majority bucket mitosis and the highest `chr_frag` over all 49 ROIs, picked in the same shuffled stable order, with the 7.5 µm rule applied within an ROI. It has one end only.
  - **Why:** it shows whether a fragmentation penalty would demote one mitotic phase, for example anaphase or telophase pairs.
- **Each tile** is V2's tile:
  - the 61 px RGB window centred on the location, with `M_chr` outlined;
  - a border coloured by bucket: mitosis green, look-alike orange, unannotated grey;
  - the raw value and the sign as a caption;
  - at least 2 screen pixels per image pixel when you open it.
- **Re-running the refinement for the tiles:** `M_chr` comes from re-running the refinement at the tile with `W.montage_windows(dict(path=..., mpp=..., locations=[...]))`, after the SHA-1 assert. `montage_windows` returns `rgb61` and `m_chr51` per location. Both windows are centred on the location, so the 51 px mask sits at offset 5 in the 61 px tile.
- **Optional reuse:** you may extract `caption_strip`, `montage_tile` and `hstack_tiles` from cell 32 of V2's notebook with `ast`, with the constants they read: cell 32's first line (`TILE_GAP`, `CAP_H`, ...) and cell 31's `BUCKET_RGB` and `OUTLINE_RGB`. `BUCKET_RGB` is keyed by short names: map `bucket` with `{'human_correct_label': 'mitosis', 'human_rejected_label': 'look-alike', 'non_human_findings': 'unannotated'}` (V2's `BUCKET_SHORT`, cell 23) before taking the majority. Do not run cells 31–32, and do not reuse `pick_tiles` (V2's tie rule) or `montage_sign` (it reads the depth-300 table).

**The findings document:** `production_hematoxylin_only/multiotsu_pad_arms/candidate_features_a3_pool100_findings.md`.
- **It opens with:**
  - the scope: the 100-peak pool only, everything in-sample, no unseen evidence;
  - what changed from V2, briefly;
  - what is wrong:
    - the features that do not measure what V2 claimed;
    - the measured orientations that disagree with `direction_hypothesis`;
    - the `non-monotone` features;
    - the covariate-dependent features.
- **Then one section per feature:**
  - what it measures (5.9);
  - its three readings, with CIs and labels, and `beyond_darkness`;
  - the top-10 against bottom-10 counts;
  - failing9 against other40.
- **For `chr_frag`, also:** the montage judgement and its largest redundancy pairs.
- **Then:**
  - the uncut concavity comparison;
  - the redundancy pairs above 0.8;
  - the expected-value comparison;
  - "What I changed";
  - what you could not check.

**The final markdown cell** holds, in this order:
1. the per-feature verdicts in your own words;
2. the montage judgement;
3. "What I changed";
4. what you could not check.

**Two text sections in existing documents.** These are this run's only edits to existing files.
- **Top of `production_hematoxylin_only/multiotsu_pad_arms/candidate_features_a3_top300_findings.md`:** insert a dated section, "Superseded (<date>)".
  - It points to the V3 findings document.
  - It states that every depth-300 and beyond-pool number in the document is out of scope, and that the neighbour columns are dropped.
  - It lists the superseded statements:
    - **§1.1:**
      - the `chr_tophat_edge` row ("Fails its claim and its direction"; "Measured `x-1`: depth 300 AUC 0.389"): in the pool it is non-monotone and measures fragmentation;
      - the `h_chr_lap` row ("thinness plus texture");
      - the `chr_radial_hf` row ("does not hold deeper"): out of scope;
      - the `chr_area_over_nbr` row: dropped.
    - **§1.2:** the depth-300 column and the beyond-pool notes are out of scope.
    - **§1.3:** `chr_radial_hf` is `raw` in the pool, not non-monotone.
    - **§1.4:** item 1 (the prefix-NMS search) is replaced by the 100-peak search; item 3's neighbour columns are dropped.
    - **§3:** `chr_area_um2`, `h_chr_mean` ("beyond od", with the 0.571 floor), `dark_area_frac_101` (the `od` inflation), `chr_tophat_edge`, `h_chr_lap`, and the neighbour columns, as in Section 5.9.
    - **Throughout:** "beyond darkness" means "beyond `od`".
    - **§4:** "Reticulated `M_chr`" is now also a feature: `chr_frag` measures it (fragmented on 28 % of pool rows). Whether it is the right mask stays a refinement question.
- **End of `production_hematoxylin_only/multiotsu_pad_arms/candidate_features_a3_top300_audit.md`:** append "Errata (verified 2026-10-09)".
  - It holds `CANDIDATE_FEATURE_EXTRACTION_FIX_PLAN.md` Section 1's rows 4, 5, 6, 9, 10, 11, 17 and 18.
  - It adds the point the audit missed, from the paragraph under that table.
  - It notes that the audit's depth-300 and beyond-pool material is out of scope.
  - Copy the verdicts and numbers as they stand.

## 9. Environment, conventions, order of work

**Environment:**
- **Interpreter:** `/Users/mohinianand/anaconda3/bin/python3`, the only one that imports `cv2` and `midog_utils`.
- **Execute with** `/Users/mohinianand/anaconda3/bin/jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=3600 <path>`.
- **Checkout:** run in the main checkout, not a git worktree. The ROI TIFFs under `images/` are not tracked, so a
  worktree has no montage images unless each TIFF is linked in one by one.
- **Paths:** the notebook sits two levels below the repo root. Use `sys.path.insert(0, '../..')` and `sys.path.insert(0, '.')` (for the V2 worker); images are at `../../images/...` and results at `../../results/`.

**Code conventions:**
- **Copied code is not retyped.** Extract V2 helpers with `ast` from the cells named in Sections 6.2 and 8, and import the V2 worker as `W`.
- **Docstrings:** the opening `"""` on its own line, a brief summary, one line per parameter with type and a short description, and one line for the return.
- **Layout:** function definitions and calls go on a single line, never wrapped.

**Rules:**
- Do not commit.
- Do not modify or delete any existing file, except for the two text sections of Section 8.
- Kill only a PID you started.

**Order of work:** Section 0.1.

**Report back with, in this order:**
1. the findings document's path and its opening section;
2. the validity table: all 24 rows, with the three readings (AUC, CI, label) and `beyond_darkness`;
3. the expected-value comparison of Section 6.4, with any mismatch and its cause;
4. the montage judgement;
5. the redundancy pairs above 0.8, and `chr_frag`'s largest pairs;
6. "What I changed".

Lead with what is wrong, not with what reproduced.

## 10. Not in this spec

- Any ranking, combination, weighting, selection, statistics beyond the descriptive validity report, or verdict.
- Fresh clicks, held-out ROIs, or any search other than production's A3 search (Section 3.1).
- Changes to the production refinement, template, floor clip or matcher, or to the benchmark harness's border band.
  The only departure from production is the band, which is the harness's correction for ROI crops (Section 3.1).
- A nucleolus detector or a learned feature. If a cue in Section 5 turns out to need one, say so in the log and leave
  it out.
- The necrosis context of C10 (VCGP 2021: "A structure within a region of necrosis ... is more likely to be MLF"). No
  column measures it.
- A fine-arm detector beyond `chr_tophat_edge`, which measured fragmentation instead (5.9). A mask-level arm counter (an opening residue attached to the core) was
  tested on 2026-10-09 and rejected. It fired on a notched disc, and it carried no signal on the recorded pool.
- **Stain-normalised channels.** No per-slide (Macenko) stain vectors anywhere in the features; Section 5.8 (a) uses
  the fixed `rgb2hed` vectors on purpose.
- The floor stress test and the multi-round feedback work (`FUTURE_IDEAS.md` §1 and §4.1).
- Any candidate beyond the 100-peak pool, until the user asks for it (2026-10-09).
- The neighbour columns (5.8 d), and the audit's ROI-wide rebuild of them.
- **The audit's other proposed reworks, declined in the fix plan (its Section 3.4):**
  - a within-score-quintile AUC;
  - `dark50` or `disc_h` as the darkness control;
  - a ROI-heterogeneity statistic;
  - `tophat_fix`;
  - `eh_territory`;
  - a sub-pixel `chr_radial_hf`;
  - a sector rule for `ring_peak_prominence`;
  - a local-step `chr_edge_sharpness_per_um`;
  - ROI-median subtraction for `chr_resid_frac`.
- Recomputing any V2 feature in this run.
