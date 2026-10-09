# Candidate feature extraction (V2 run, d99e93a): verification of the audit and fix plan

**Status:** plan only, written 2026-10-09 and revised the same day with your decisions. Nothing in Section 3 has been
run. The numbers in Sections 1 and 2 are verification numbers I computed from the committed tables. They are
in-sample, like everything else on this run.

**Corrected 2026-10-09, in the review of V3** (`CANDIDATE_FEATURE_EXTRACTION_PROMPT_V3.md`, which applies this plan):
- Section 1: rows 4 (ρ), 9 (two CIs), 10 (the numerator), 11 (the SD range), 17 (both distance bases) and 18 (wording);
- Section 2: items 3 (the `od` CI) and 7 (the NaN share); Section 3.4's heterogeneity row.

Every corrected number uses V3's bootstrap (Section 6.2) and each column's own definition.

**Your decisions (2026-10-09):**
- **The candidate set is production's 100-peak pool only** (Section 0). Nothing beyond it is analysed for now: no depth-300 table and no beyond-pool readings.
- **The neighbour columns are dropped.**

**Inputs:**
- the run: `candidate_features_a3_top300.ipynb`, its findings document, and `results/candidate_features_a3_top300_*`;
- the audit: `candidate_features_a3_top300_audit.md` and `.py` (a4a6c0f), and its `results/..._audit_*` tables.

**Verification scripts:** session scratchpad `verify/v1.py`–`v12b.py`. They read `candidates.csv` and the audit's
committed `audit_locations.csv.gz`. Bootstrap: 49 ROIs, 4,000 draws, seed 20261009, as the audit used.

**Bottom line.**
- The run's code and arithmetic are sound, as the audit found.
- The pool `validity.csv` orientation labels are right.
  - A CI gate on top of the V2 sign rule changes none of them.
  - It changes none of the 72 labels in the three readings of Section 3.1 either.
- What went wrong is reading:
  - the headline numbers came from the depth-300 AUC, which is now out of scope;
  - "beyond darkness" rested on one covariate, `od`, which holds dark mass fixed, not core darkness;
  - three mechanisms were inferred by eye and misread.
- The audit gets several of these right. It also gets three things wrong:
  - its diagnosis of the `od` strata;
  - the statistic behind its `chr_tophat_edge` mechanism;
  - its `e_h_angle` heterogeneity claim.
- **The fix is small:**
  - a re-tabulation from the committed tables, with no image pass;
  - one added column (`chr_frag`) that already exists in the audit's locations file;
  - corrected text.
- **Not needed:** no search, refinement or feature re-run, and no gates. Nothing that is re-computed could fail them.

---

## 0. The candidate set: the 100-peak pool

**Definition.** For each of the 294 clicks, run the A3 search: the `mo_pad61` template on the floor-clipped channel,
with the 30 px border band.
1. Extract the score peaks above the peak threshold (median + `DEEP_FLOOR_Z` × MAD, with `DEEP_FLOOR_Z` = −1.5).
2. Keep the 100 highest-scoring peaks (`MAX_PEAKS` = 100, D9).
3. Apply NMS at `ev.radius_px(mpp)` (7.5 µm) to those 100 only.

That leaves 79–100 candidates per click, 29,025 rows in all: 28,883 accepted, 7,664 of them mitoses.

**What the threshold is in practice.** Every click has at least 16,786 peaks above the peak threshold, so the
100-peak cap always binds. The effective threshold is the 100th peak's score.

**No re-run.** This is the production search, so its output is already in the committed tables.
- In `candidates.csv` it is exactly the rows with `in_recorded_pool`. The run computed them as
  `nms_by_distance(centers[:MAX_PEAKS], ...)` and asserted that they are the deep list's prefix.
- Gate B matched them against the recorded A3 run: 293 clicks exact, 1 adjacent-pixel near-tie.
- On every click, the pool row count equals `n_detections_recorded`.

**Every future run of this stage:**
- use the verbatim A3 `detect` from `tm_floor_definition_49roi_6seed.ipynb` cell 5, with `MAX_PEAKS` = 100;
- drop `detect_deep`, the uncapped peaks, the prefix NMS and the top 300;
- drop the beyond-pool readings and the depth-300 table.

**Out of scope until you ask for it:**
- every `cand_rank` beyond the pool;
- `validity_depth300.csv`;
- `auc_wc_beyond_pool` and `auc_dark_beyond_pool`;
- the audit's depth-300 decomposition.

The committed files stay in git and are not used.

**Two costs of the narrower scope:**
- **The unblinded columns lose their unseen evidence.** `chr_concavities_per_um`, `dark_concavities_per_um`, `chr_solidity`, `chr_radial_hf` and `chr_tophat_edge` were looked at on this pool before `directions.csv` was frozen, and the beyond-pool rows were their only unseen evidence. Their pool labels stay in-sample.
- **failing9 has little data.** Its readings rest on 412 pool mitosis rows.

## 1. The audit, claim by claim

| # | audit claim | verdict | evidence |
|---|---|---|---|
| 1 | No implementation defect; every AUC reproduces | **Confirmed** | `reproduce_max_abs_diff` 0.0; `tophat_check` and `resid_check` 0.0 |
| 2 | The depth-300 pooled AUC is 54 % pool-mitosis × beyond-other pairs, so it tracks search rank | **Confirmed; now moot** | `h_chr_mean` reads 0.918 on those pairs, against 0.761 in the pool. Section 0 drops depth 300 altogether |
| 3 | Montage orientation inherits (2) | **Confirmed in principle, no consequence** | A `non-monotone` column is sorted by the side of its AUC. Every montage column has the same side in the pool and at depth 300, so no montage changes and none needs re-rendering |
| 4 | `od` quintiles "do not hold darkness fixed"; use `dark50` or `disc_h` | **Wrong diagnosis** | Quintile strata of about 20 rows leave every darkness measure separating within its own strata: `od` 0.571, `h_chr_mean` 0.576, `disc_h` 0.579, `dark50` 0.586. `h_chr_mean` reads 0.604 within `od` strata, 0.598 within `dark50` and 0.639 within `disc_h`, so neither replacement holds darkness better. `dark50` is nearly `h_chr_mean` (within-click ρ 0.92). Caliper matching (±5 percentile) moves no non-darkness column by more than 0.015 (the largest is `chr_area_um2` at matched `h_chr_mean`, 0.604 → 0.619), so the coarseness only matters for the darkness columns themselves. **The real problem is which covariate:** see Section 2, item 3 |
| 5 | "It printed the darkness leak and did not flag it" | **Wrong** | The findings say `h_chr_mean` "still carries darkness information beyond the 51 px window's od". The run did not carry that into the other columns' "beyond darkness" sentences; that part stands |
| 6 | `chr_area_um2` becomes +1 at matched darkness in the pool (and null beyond it) | **Partly right; stated too broadly** | At matched *core* darkness it is +1: `h_chr_mean` 0.604, `dark50` 0.566 [0.534, 0.599], `disc_h` 0.544 [0.510, 0.577]. At matched `od` it reads 0.458 and at matched score 0.430. Why they differ is in Section 2, item 4. The beyond-pool half is out of scope |
| 7 | `chr_area_over_nbr` and `chr_hole_frac` are mostly darkness proxies in the pool | **Confirmed** | `chr_hole_frac`: 0.460 within `od` (ROI share 0.61) and 0.493 within `h_chr_mean`, so neither is consistent. `chr_area_over_nbr`: dropped (Section 2, item 7) |
| 8 | `chr_radial_hf` is robust in the pool (and inconclusive beyond it) | **Pool half confirmed** | 0.593 / 0.625 at matched `od` / `h_chr_mean`. The beyond-pool half is out of scope |
| 9 | `chr_tophat_edge`'s x-1 is `M_chr` fragmentation plus noise: 0.524 within `chr_frag` quintiles | **Conclusion right; statistic wrong** | 72 % of pool rows have `chr_frag` = 0. `rank(method='first')` splits that tied block in row order, which is `cand_rank` order, so three of the five "fragmentation" strata are score strata (score strata alone give 0.518). With ties shuffled the AUC is 0.486; with the zero block as one stratum, 0.499. **The clean statistic:** on single-piece `M_chr` it is 0.501 [0.466, 0.535]. The "plus noise" part is not supported: smoothed, it reads 0.530 [0.493, 0.563] on single-piece rows |
| 10 | `h_chr_lap` measures texture, not thinness; the ratio is not a division artefact | **First half confirmed, second half wrong** | The x-1 holds within area quintiles (0.395) and on single-piece rows (0.410). (Within `chr_frag` quintiles it reads 0.421, but that stratifier has the tie-order artefact of row 9; with ties broken at random it reads 0.403.) But the numerator alone (`h_chr_lap` × `h_chr_mean`, on the column's own support) is +1 marginally (0.541 [0.497, 0.585], 33 of 49 ROIs), null at matched `disc_h` (0.482), and x-1 only under `h_chr_mean` (0.414), the control the audit itself calls over-correcting. So the ratio's x-1 does come through the division: mitoses have more absolute texture because they are darker, and less relative to their darkness. (The audit's `lap_sd` reads 0.572 / 0.501 / 0.431, but it is not this numerator: it skips the clipped-pixel exclusion, and it differs from it on 77 % of 3D Histech rows against 0.1 % of Aperio rows) |
| 11 | `e_h_angle_ring23_deg`'s pooled null "hides strong ROI heterogeneity" (SD 0.148 against a 0.042 null) | **Wrong as a distinction** | The SD of the per-ROI AUC lies in 0.094–0.156 for 21 of the 24 pool features (`h_chr_lap` 0.146, `white_frac` 0.139, `pair_sep` 0.156; median 0.132). The exceptions are `brown_frac_surround` 0.018 (ties at 0), `chr_frag` 0.075 and `pair_parallel_deg` 0.197. The within-click null understates the noise. An object-level null (4 px object clusters, which is a lower bound) gives about 0.07. `e_h_angle` sits with `h_chr_lap`, `pair_area_ratio` and `ring1_frac_pale` (0.146–0.148), so it is not distinctive. Its extremes sit on ROIs with 5–12 mitoses per click, and no label-free predictor of the per-ROI sign was shown. "No signal for a fixed orientation", the run's reading, stands |
| 12 | `chr_resid_frac`'s non-pigment x-1 is not a glass, rim or support artefact | **Confirmed** | glass I0 0.383; `M_chr` only 0.359. The two strata rows are not in the audit script, but they reproduce: `chr_ns_share` quintiles 0.358, support OD norm quintiles 0.370 |
| 13 | `dark_concavities_per_um` leads `chr_concavities_per_um` by +0.020 [+0.002, +0.039] on uncut pool rows | **Confirmed** | Same values |
| 14 | `M_dark` contains `M_chr` on about 99 % of rows | **Confirmed** | 1.1 % of others and 0.2 % of mitoses fail |
| 15 | Gate B was adapted after it halted | **Procedurally true; nothing to fix** | The 1 px move comes from a 1.2e-7 score difference, and the `od` shift follows from the window moving. That is "a difference a float shift below 1e-4 can explain", the spec's own clause. The candidates are unaffected |
| 16 | Gate C's "64.2 %" mixes populations | **Confirmed; nothing to fix** | 15,708 / 27,612 = 56.9 %. The label is in the notebook output only, not in the findings |
| 17 | "The six clicks per ROI re-find the same objects one pixel apart" | **Overstated** | Distances are from each pool location to the nearest pool location of the other click or clicks; "adjacent" is under 2 px, the same or a neighbouring pixel. Against one other click of the same ROI, averaged over the pairs, 55–59 % of a click's pool lies within 8–15 px and 16 % is adjacent. Against all five other clicks together, 85–86 % lies within 8–15 px and 46 % is adjacent. So most objects are re-found, but under half one pixel apart. The ROI is still the right resampling unit |
| 18 | "The script reproduces every number" | **True, with two numbers outside the script** | The script does not compute the two `chr_resid_frac` strata of row 12. Recomputed separately, both reproduce |
| 19 | Fixes: score-quintile AUC, `dark50`/`disc_h`, ROI heterogeneity statistic, `tophat_fix`, `eh_territory`, sub-pixel `radial_hf`, ROI-wide neighbours, ring sectors, local edge step, resid ROI-median | **Declined, or dropped by your decisions** | See Section 3.4 |

**What the audit missed.** `dark_area_frac_101`'s strong reading at matched darkness is partly `od`'s own construction.
- It reads 0.407 at matched `od`, 0.436 at matched `h_chr_mean` and 0.470 at matched score.
- `od` averages the window's darkest 260 px, so dark neighbours enter it: within click, ρ(`od`, `dark_area_frac_101`) is 0.23, against 0.11 for `h_chr_mean`.
- So the run's "among equally dark candidates, mitoses have fewer dark neighbours" is half `od` artefact. At matched core darkness the effect is smaller, but it is consistent (Section 5).

## 2. What was wrong in the run's execution and reading, after verification

1. **Headlines came from the depth-300 AUC.**
   - Five columns change label between depth 300 and the pool:
     - `chr_tophat_edge` (x-1 → non-monotone);
     - `chr_radial_hf` (non-monotone → raw);
     - `ring1_frac_pale` (x-1 → non-monotone);
     - `chr_area_over_nbr` (x-1 → non-monotone);
     - `pair_axis_deg` (non-monotone → raw).
   - With Section 0, only the pool labels exist.
2. **No uncertainty was reported.** Adding ROI-bootstrap CIs is right for reporting, but as a gate it changes no label. In all 72 cells (24 columns × the three readings of 3.1), the existing sign-consistency rule (≥ 2/3 of ROIs, ≥ 5/7 domains) already rejects every column a CI would reject. So the CI is reported beside each AUC and is not added as a rule.
3. **"Beyond darkness" means "beyond `od`".** `od` and `h_chr_mean` are tied as pool darkness rankers: 0.770 [0.739, 0.799] against 0.761 [0.729, 0.794]. Conditioning on `od` is decision-relevant, but `od` is not core darkness.
   - It counts surround on cores under 260 px, which is 36 % of pool rows.
   - It counts dark neighbours.
   - The existing best rule, from 2026-09-30, is a rank average of `od`, core OD mean, not-white and solidity. So a column's beyond-darkness orientation has to hold at matched `od` and at matched `h_chr_mean` (the core OD mean) before a rule can rely on it.
4. **`chr_area_um2`'s orientation depends on what is held fixed.**
   - **At equal dark mass (`od`):** no consistent relation. ρ(`od`, area) is −0.02 overall. On cores under 260 px it is +0.10, because a bigger core fills more of the darkest 260 px; on larger cores it is −0.07. ρ(`h_chr_mean`, area) is −0.37 and −0.32 on the same two groups.
   - **At equal core darkness:** mitoses are larger. ρ with area is −0.32 to −0.50, and the dark imposters (small round dark nuclei) are smaller.
   - **The reading for a ranker:** with `h_chr_mean` in it, area is +1. With only `od`, area adds nothing. "Mitoses have smaller cores" is a darkness effect.
5. **`chr_tophat_edge` was read as "reticulum, mitoses have less thin structure".**
   - It is null on single-piece `M_chr` (0.501).
   - Its signal is `M_chr` fragmentation: the share of the top Otsu class outside its largest 4-connected piece.
   - **That is not the refinement unioning pieces.** 86 % of fragmented rows have `n_pieces` = 1, and `n_pieces` itself carries nothing (0.491).
   - **Fragmentation is direct and stronger.** Mitoses are fragmented on 18 % of pool rows against 31 % of others. Measured directly, the AUC is 0.446 [0.426, 0.464], ROI share 0.82, and x-1 under both darkness covariates.
   - V2 had dropped `chr_n_pieces`, which was a different count.
6. **`h_chr_lap` was read as "thinness plus texture".** It is relative texture (Section 1, row 10). The x-1 holds at matched `h_chr_mean` (0.395). At matched `od` it is weaker and not consistent (0.422, ROI share 0.63).
7. **The neighbour columns are dropped (your decision).**
   - **What goes:** `chr_area_over_nbr`, `nbr_n`, `nbr_nuc_area_med_um2`, and the three `spec_` versions built at `t_nuc`.
   - **Why:**
     - they are NaN on 89 % of accepted pool rows;
     - on the rest, `chr_area_over_nbr` is `chr_area_um2` (ρ 0.78; it reads 0.457 against 0.460 marginally);
     - within a click the neighbour median barely varies, so a within-click rank of `chr_area_um2` already gives the normalisation.
   - **What happens to them:** they stay in the committed `candidates.csv` and appear in no corrected table and no downstream input. The audit's ROI-wide rebuild is not done.
8. **`ring_peak_prominence`, `e_h_angle_ring23_deg` and `chr_resid_frac`** were read correctly by the run. The audit's relabels for `ring_peak_prominence` and `chr_resid_frac` are right. Its `e_h_angle` heterogeneity point is not (Section 1, row 11).
9. **`brown_frac_surround`** is montage-correct, but carries almost no ranking signal even where it should. On 291 it reads 0.490, and on 410 0.440 with 5.5 mitoses per click. On 123 it fires on 29 % of pool rows, with tiny values. Its pooled 0.495 is ties.

## 3. The fixes

### 3.1 Corrected validity table: one script, no image pass

**Files.** The script is `production_hematoxylin_only/multiotsu_pad_arms/candidate_features_a3_top300_validity_fix.py`.
It writes `results/candidate_features_a3_top300_validity_corrected.csv`.

**Inputs:**
- `candidates.csv`, restricted to accepted rows with `in_recorded_pool` (Section 0). No other row is read into any statistic.
- `chr_frag` from `candidate_features_a3_top300_audit_locations.csv.gz`, joined on `(file_name, cx, cy)`. The audit computed it on the run's verbatim `M_chr`; its `tophat_check` reproduced the run's column exactly from that mask.

**Rows:** 24.
- The 23 role-`feature` columns left after dropping `chr_area_over_nbr`.
- `chr_frag`.

**Three readings per column,** each by the existing V2 AUC rule (within click, ties ½):
- **`pool`:** within click.
- **`pool_od`:** within click × within-click `od` quintile.
- **`pool_hmean`:** the same with `h_chr_mean` quintiles.
  - Strata are formed once per click, before any column's NaN rows are dropped.
  - Ties are broken at random with a fixed seed (shuffle, then `rank(method='first')`).
  - Only these two continuous covariates are ever stratified on.

**For each reading:**
- `auc`, with `ci_lo` and `ci_hi`: the ROI cluster bootstrap, 4,000 draws, seed 20261009, U and n1·n0 summed per ROI;
- `roi_share` and `domains`: the same definitions as V2;
- `orientation`: the V2 rule unchanged (`raw` / `x-1` / `non-monotone`). The CI is not part of it; verified, it would change none of the 72 labels.

**Plus:**
- `beyond_darkness`: the shared label when `pool_od` and `pool_hmean` agree. `entangled` when one is `raw` and the other `x-1`. Otherwise `od: <label> / h_chr_mean: <label>`.
- `measures`: the corrected one-line label (3.3).

**What stays where it is.** The committed `validity.csv` stays the source for these columns. All of them are pool-based and unaffected:
- top-10 / bottom-10;
- scanner;
- failing9 / other40;
- NaN rates;
- `synthetic_pass`.

**What is dropped:**
- `validity_depth300.csv` and every beyond-pool number (Section 0);
- `auc_area_tertiles`: no boundary column moved by more than 0.01;
- `auc_hole_*`: undefined, because no click has both classes above 0.3;
- `auc_absdev`.

**Runtime and output.** About a minute. One printout states each covariate's AUC within its own quintiles (`od` 0.571,
`h_chr_mean` 0.576). That number is the strata's resolution: a darkness column's beyond-darkness AUC is read against
it, not against 0.5.

### 3.2 One column added: `chr_frag`

**Definition.** 1 − (pixels of whole `M_chr`'s largest 4-connected piece) / (pixels of whole `M_chr`). 0 for a single
piece; NaN on refused rows.

**Its record:**
- Role: `feature`. Added post hoc, so it has no direction hypothesis; it gets no `directions.csv` row, because that file stays frozen.
- **It supersedes `chr_tophat_edge`** in the downstream. Don't carry both: `chr_tophat_edge` is a noisier copy (Section 2, item 5).
- **No new montage.** The high end of `chr_tophat_edge`'s montage already shows the fragmented masks.

### 3.3 Corrected labels (the `measures` column, and the findings text)

| column | measures |
|---|---|
| `chr_tophat_edge` | `M_chr` fragmentation, noisily. Null on single-piece masks. Not a thin-arm measure. Superseded by `chr_frag` |
| `chr_frag` | whether the window's top Otsu class inside the selection breaks into separate 4-connected pieces |
| `h_chr_lap` | relative chromatin texture (Laplacian SD over core mean). Lower on mitoses at matched core darkness |
| `chr_area_um2` | core area. x-1 marginally; +1 at matched core darkness; no relation at matched `od` (2.4) |
| `chr_hole_frac` | holes in the largest piece. x-1 marginally; no consistent relation at matched darkness |
| `chr_radial_hf` | outline high-harmonic energy. `raw`, also at matched darkness. Overlaps `chr_solidity` (ρ −0.67) |
| `dark_area_frac_101` | dark area around the candidate (crowding). No marginal orientation; x-1 at matched darkness under both covariates. Its `od` value is inflated (Section 1) |
| `ring_peak_prominence` | dark mass 1–3 µm out (a neighbour, or a darker nucleus around a clump), not an envelope |
| `chr_resid_frac` | chromatin colour off the H–E plane. Pigment is its high tail; it is x-1 on non-pigment ROIs too |
| `e_h_angle_ring23_deg` | the neighbourhood's eosin (red cells, stroma). No orientation |
| `brown_frac_surround` | brown pixels beside the candidate. A pigment-ROI flag with no ranking signal on its own ROIs |

Every other column keeps its V2 label.

### 3.4 Audit fixes declined, and why

| proposed | why not |
|---|---|
| Within-score-quintile AUC at depth 300 | Depth 300 is out of scope (Section 0) |
| `dark50` / `disc_h` as the darkness control | Neither holds darkness better than `od`. `dark50` is about `h_chr_mean`. The two ranker covariates are the decision-relevant ones. They stay cited for `chr_area_um2` only, as its mask-free support |
| ROI-heterogeneity statistic with an object-level null | 21 of the 24 columns are heterogeneous at the same scale (0.09–0.16; Section 1, row 11). A fixed-orientation ranker cannot use a per-ROI sign without a predictor for it |
| `tophat_fix` | Null on single-piece rows (0.501). `chr_frag` carries its signal directly |
| `eh_territory` | 0.463 [0.424, 0.506]; failing9 0.48; same ROI spread. A cell-territory build for a column VCGP calls "not reliable" |
| Sub-pixel `chr_radial_hf` | Its role is covered by `chr_solidity` and the concavity columns, which are stronger and pass the scale rule or are exempt from it |
| ROI-wide neighbour components | The neighbour columns are dropped (your decision, Section 2, item 7) |
| Sector rule for `ring_peak_prominence` | No failing-ROI imposter class is "clump in an envelope". C2 is unmeasurable on uniformly dark nuclei by the spec's own test |
| Local-step `chr_edge_sharpness_per_um` | A diagnostic with no consumer: it did not track the per-scanner AUCs |
| ROI-median subtraction and pigment flag for `chr_resid_frac` | Within-click AUCs are unaffected. That is a downstream pooled-model choice |

### 3.5 Text corrections

**The findings document.** Add a dated "Corrections (verified 2026-10-09)" section at its top.
- **Scope:** state it first. The analysis is limited to the 100-peak pool (Section 0); every depth-300 and beyond-pool number in the document is out of scope, and the neighbour columns are dropped.
- **Then point to the corrected table** and supersede these statements:
  - **§1.1:**
    - the `chr_tophat_edge` row ("x-1, 0.389"; "mitoses have less thin structure");
    - the `h_chr_lap` row ("thinness plus texture");
    - the `chr_radial_hf` row ("does not hold deeper"): out of scope;
    - the `chr_area_over_nbr` row: dropped.
  - **§1.2:** the depth-300 column and the beyond-pool notes are out of scope.
  - **§1.3:** `chr_radial_hf` is `raw` in the pool, not non-monotone.
  - **§1.4 item 1:** the prefix-NMS search is replaced by the 100-peak search for future runs (Section 0).
  - **§1.4 item 3:** the neighbour columns are dropped.
  - **§3 entries:**
    - `chr_area_um2`: 2.4;
    - `h_chr_mean`: "beyond od", with the 0.571 floor;
    - `dark_area_frac_101`: the `od` inflation;
    - `chr_tophat_edge` and `h_chr_lap`: as above;
    - the neighbour columns: dropped.
  - **Wording, throughout:** "beyond darkness" becomes "beyond `od`".
  - **§4:** "Reticulated `M_chr` ... not measured" becomes "measured as `chr_frag`, 28 % of pool rows".

**The audit document.** Add an "Errata (verified 2026-10-09)" section at its end. It holds:
- Section 1's rows 4, 5, 6, 9, 10, 11, 17 and 18;
- the missed `dark_area_frac_101` point;
- a note that its depth-300 and beyond-pool material is out of scope.

**Memory.**
- Update `candidate-features-v2-audit` and `candidate-feature-extraction-spec`.
- Their `chr_tophat_edge` (0.389), `od`-leak and `e_h_angle` statements are superseded. Done with this plan.

## 4. Where the floor clip is used

The candidates come from the floor-clipped search at the adopted level, 1.0 × `b_otsu` (the D3 amendment of
2026-10-02).
- The worker sets `b = fl['b_otsu']` and searches on `clip_floor(hem, b)`, which is `np.maximum(hem, b)`.
- The template is cut from that clipped channel, and D11's blank fills with `b`.
- No other level has been adopted: Stage B, which would test other levels, is deferred to before shipping.
- **The clip is not in `midog_utils.production` yet.** The code is copied verbatim from `tm_floor_definition_49roi_6seed.ipynb`.

The features are measured on the unclipped channel, by design. The D3 amendment states it: "the box refinement,
`chromatin_od` and every per-candidate measurement keep reading the unclipped one."
- That covers the refinement, `M_chr`, `M_dark`, every feature and `od`.
- The surround columns need the real background that the clip flattens.

The clip value does enter one feature path: as a threshold. `t_nuc` is the same number as `b`.
- `ring1_frac_above_tnuc` and `ring1_frac_pale` split ring 1 at `t_nuc`.
- So a pixel "above `t_nuc`" is one the clip leaves untouched.

## 5. Expected corrected table (verification run, 100-peak pool)

Each cell is AUC [95 % ROI CI], then ROI share / domains. Every label is the V2 sign rule.

| column | pool | label | at matched `od` | at matched `h_chr_mean` | beyond darkness |
|---|---|---|---|---|---|
| `chr_area_um2` | 0.460 [0.427, 0.493] 0.67/5 | x-1 | 0.458 [0.422, 0.494] 0.61/5 | 0.604 [0.574, 0.634] 0.84/7 | od: non-monotone / h: raw |
| `chr_elongation` | 0.636 [0.599, 0.671] 0.78/7 | raw | 0.627 0.76/7 | 0.601 0.67/6 | raw |
| `chr_solidity` | 0.375 [0.340, 0.411] 0.88/7 | x-1 | 0.361 0.86/7 | 0.351 0.94/7 | x-1 |
| `chr_concavities_per_um` | 0.621 [0.595, 0.646] 0.90/7 | raw | 0.642 0.96/7 | 0.669 0.98/7 | raw |
| `dark_concavities_per_um` | 0.651 [0.623, 0.680] 0.94/7 | raw | 0.674 1.00/7 | 0.685 0.98/7 | raw |
| `chr_hole_frac` | 0.435 [0.409, 0.459] 0.67/6 | x-1 | 0.460 0.61/6 | 0.493 0.47/4 | non-monotone |
| `chr_radial_hf` | 0.561 [0.530, 0.591] 0.82/6 | raw | 0.593 0.84/7 | 0.625 0.92/7 | raw |
| `pair_sep_um` | 0.581 [0.540, 0.619] 0.71/6 | raw | 0.624 0.83/7 | 0.624 0.81/7 | raw |
| `pair_area_ratio` | 0.542 0.65/7 | non-monotone | 0.538 0.56/7 | 0.526 0.56/6 | non-monotone |
| `pair_parallel_deg` | 0.517 0.46/4 | non-monotone | 0.532 0.54/3 | 0.508 0.50/3 | non-monotone |
| `pair_axis_deg` | 0.531 [0.508, 0.556] 0.69/6 | raw | 0.514 0.65/5 | 0.537 0.62/6 | non-monotone |
| `h_chr_mean` | 0.761 [0.729, 0.794] 0.96/7 | raw | 0.604 0.78/7 | 0.576 (self) | raw |
| `h_chr_lap` | 0.397 [0.355, 0.439] 0.71/5 | x-1 | 0.422 0.63/6 | 0.395 0.71/6 | od: non-monotone / h: x-1 |
| `ring_peak_prominence` | 0.418 [0.394, 0.445] 0.78/6 | x-1 | 0.405 0.84/7 | 0.417 0.82/6 | x-1 |
| `ring1_frac_above_tnuc` | 0.634 [0.594, 0.677] 0.80/7 | raw | 0.594 0.76/6 | 0.609 0.80/7 | raw |
| `ring1_frac_pale` | 0.422 0.63/6 | non-monotone | 0.461 0.55/5 | 0.446 0.57/5 | non-monotone |
| `white_frac_gray215_61` | 0.363 [0.318, 0.407] 0.82/7 | x-1 | 0.347 0.84/7 | 0.354 0.80/7 | x-1 |
| `e_h_angle_ring23_deg` | 0.504 [0.463, 0.547] 0.49/3 | non-monotone | 0.509 0.51/3 | 0.487 0.57/4 | non-monotone |
| `brown_frac_surround` | 0.495 0.29/6 | non-monotone | 0.495 | 0.494 | non-monotone |
| `dark_area_frac_101` | 0.466 [0.434, 0.496] 0.65/5 | non-monotone | 0.407 0.86/7 | 0.436 0.73/7 | x-1 |
| `chr_resid_frac` | 0.373 [0.336, 0.410] 0.80/7 | x-1 | 0.393 0.78/7 | 0.399 0.78/7 | x-1 |
| `chr_tophat_edge` | 0.460 [0.422, 0.497] 0.57/4 | non-monotone | 0.517 0.57/5 | 0.472 0.47/5 | non-monotone |
| `chr_tophat_core` | 0.481 0.47/4 | non-monotone | 0.520 | 0.507 | non-monotone |
| `chr_frag` (added) | 0.446 [0.426, 0.464] 0.82/7 | x-1 | 0.460 0.73/6 | 0.469 0.69/5 | x-1 |

**Reading it.**
- **Robust, in the same direction at matched `od` and at matched core darkness:**
  - `raw`: `chr_concavities_per_um`, `dark_concavities_per_um`, `chr_elongation`, `chr_radial_hf`, `pair_sep_um`, `ring1_frac_above_tnuc`;
  - `x-1`: `chr_solidity`, `ring_peak_prominence`, `white_frac_gray215_61`, `chr_resid_frac`, `chr_frag`.
- **Covariate-dependent:** `chr_area_um2` and `h_chr_lap`.
- **Conditional only:** `dark_area_frac_101`.

The executing agent's AUCs must match these to ±0.005. The CIs may differ slightly with how the bootstrap draw is
implemented.

## 6. Execution notes

- **Interpreter:** `/Users/mohinianand/anaconda3/bin/python3`.
- **Writes only:** the script and the corrected table (3.1), and the two text sections (3.5).
- **Leaves alone:** every committed table and the notebook.
- **Commits:** none until asked.
- **Order:**
  1. the script and table;
  2. compare with Section 5;
  3. the findings corrections, then the audit errata;
  4. report the table and any value outside ±0.005.
