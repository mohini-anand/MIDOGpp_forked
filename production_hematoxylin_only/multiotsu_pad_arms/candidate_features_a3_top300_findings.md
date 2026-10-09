# Per-candidate features on the A3 top 300: findings

**Run:** 2026-10-09, `candidate_features_a3_top300.ipynb`, on the spec `CANDIDATE_FEATURE_EXTRACTION_PROMPT_V2.md`.
- **Scope:** 49 ROIs, 294 recorded `mo_pad61` clicks, the first 300 post-NMS rows of a deep A3 search per click.
- **Size:** 88,200 candidate rows and 84,003 distinct locations. 99.7 % of rows are accepted.
- **Mitoses:** 11,152 mitosis rows, 7,705 of them in the recorded pool.
- **Outputs:**
  - tables in `results/candidate_features_a3_top300_*.csv`;
  - the spec JSON in `data/candidate_features_spec.json`;
  - the code in `candidate_features_a3_top300_worker.py`.
- **`directions.csv` SHA-1:** `620d50f307995984eed9a9c387a047b79abe5b94`. It was frozen before the first smoke run and printed by the smoke run and the full run. The three `spec_` amendment rows have their own SHA-1, `9b75f8c3…`.

**Everything below is in-sample on the 49 ROIs.**
- For `chr_concavities_per_um`, `chr_solidity`, `chr_radial_hf`, `chr_tophat_edge` and `dark_concavities_per_um`, the pool numbers were seen before freezing, on 2026-10-09. Their `validity.csv` rows are not tests.
- Their only unseen evidence is the *beyond-pool* AUC, on accepted rows with `in_recorded_pool` False (ranks of about 100–299).
- That is a deeper population, not a replication.

"Orientation" is `orientation_recommended`:
- `raw`: mitoses are higher;
- `x-1`: mitoses are lower;
- `non-monotone`: the sign is not consistent over at least 2/3 of ROIs and 5 of 7 domains.

AUC is P(mitosis > other) within a click, never folded.

---

## 1. What is wrong

### 1.1 Features that do not measure what they claim

**The montages suggest one cause behind most of the boundary family.** This is a by-eye reading of the montages, not a measured cause.
- On pale, open-chromatin interphase nuclei, `M_chr` (the 51 px window's top multi-Otsu class) looks like a branching network of chromatin clumps, not a dense core.
- These networks dominate the top or bottom tiles of the `chr_radial_hf`, `chr_solidity`, `chr_tophat_edge`, `chr_hole_frac` and `h_chr_lap` montages. This is clearest in `013`, `529`, `425` and `201`.
- So these columns appear to measure partly "is `M_chr` a reticulum, not a solid core", not only the property they name.
- It would be a property of the production refinement mask, not of any one column's arithmetic, so no column was redefined for it. A measured check (for example, `M_chr`'s fill share of its nucleus) was not run.

| column | claims | what the montage and tables show |
|---|---|---|
| `chr_tophat_edge` | fine arms / thin structure within 1 µm of the core (direction +1) | **Fails its claim and its direction.** The high end is reticulated chromatin of pale nuclei, small fragments and pigment. The low end is solid compact cores, including many mitoses. Measured `x-1`: depth 300 AUC 0.389, beyond pool 0.438, within size tertiles 0.384, so it is not only size. The vesicular-rim split (`auc_hole_gt03`) is undefined: only 228 accepted rows have `chr_hole_frac` > 0.3, holding 2 mitoses and no click with both classes. |
| `chr_radial_hf` | high-harmonic outline energy of the dense core (+1) | **The pool signal does not hold deeper.** The pool AUC 0.561 / 0.593 beyond darkness (unblinded) falls to 0.502 at depth 300 and 0.484 beyond the pool, so it is `non-monotone` at depth 300. The high end is branching chromatin networks, not spiky condensed figures. It also fails the scale rule on every shape, a failure the spec predicted. |
| `ring_peak_prominence` | a darker band 1–3 µm out: an envelope around a clump (direction 0) | The high end is a dark neighbour within 3 µm (crowding) or a clump inside a darker nucleus. The low end, exactly 0 on 43 % of rows, is isolated objects, including many mitoses and the 211 pigment. It measures nearby dark mass more than an envelope. Measured `x-1`, about 0.39–0.42. |
| `e_h_angle_ring23_deg` | eosinophilia of the cell's own cytoplasm (−1) | The high end is red blood cells and pink stroma or collagen around the candidate (`201`, `546`). The low end is bluish crowded surround. It reads the neighbourhood's eosin, not the cell's cytoplasm, and carries no signal: AUC 0.504 / 0.519, consistent in 3 of 7 domains. |
| `h_chr_lap` | fine chromatin texture (0) | The high end is small or thin masks and reticulated chromatin, where the eroded support still contains the edge ramp. The low end is large homogeneous dark cores, including mitoses. It separates strongly (`x-1`, 0.362 at depth 300, and 0.354 within size tertiles), but it measures thinness plus texture, not texture alone. |
| `chr_resid_frac` | pigment, colour off the H–E plane (−1) | **The top of the range is pigment, as claimed.** That is 211's black clumps, 289's brown-black cluster (0.39) and 339's two pigment tiles. But the signal is not pigment-specific: the AUC is 0.36 on the 40 ROIs without the pigment failure, and 0.165 on Aperio. The spec's anchor said mitoses score slightly *higher* on non-pigment ROIs, but it compared top-10 false positives with mitoses refined at their annotation centres, so the two are not like-for-like. Here, within click, mitoses score lower. Its support is label-dependent: median `chr_ns_share` is 0.58 for mitoses against 0.88 for unannotated rows. On mitoses it is measured mostly on the rim. |
| `chr_edge_sharpness_per_um` (diagnostic) | crispness of the core's boundary | The high end is chromatin clumps inside nuclei, where the step to ring 1 is small. It is not visibly crisp-edged objects, so treat it as a weak resolvability covariate. |
| `chr_area_over_nbr` | core area relative to neighbouring nuclei (0) | The spec's `t_nuc` definition is defective, and the definition was changed (Section 1.4). The corrected column is still NaN on 87–88 % of rows: the 101 px window plus the edge rule rarely leaves 3 whole nuclei. |

**These measure what they claim, on the montages:**
- `chr_concavities_per_um` and `dark_concavities_per_um`: ragged, notched outlines against smooth round cores. The high end is mostly irregular *interphase* nuclei, and mitoses are a minority even there.
- `chr_hole_frac`: rings and networks with holes against solid cores.
- `chr_solidity`, with the reticulum caveat above.
- `ring1_frac_pale`: dark cores in a pale non-glass rim, against clumps inside nucleoplasm.
- `brown_frac_surround`: golden-brown granules beside the candidate on `291` and `410`. On `123` the red-cell exclusion holds: its highest value is 0.016. No threshold was changed.

### 1.2 Measured orientation against `direction_hypothesis`

| column | hypothesis | pool | depth 300 | note |
|---|---|---|---|---|
| `chr_tophat_edge` | +1 | non-monotone (0.460) | **x-1** (0.389) | opposite to the hypothesis |
| `e_h_angle_ring23_deg` | −1 | non-monotone (0.504) | non-monotone (0.519) | no signal |
| `brown_frac_surround` | −1 | non-monotone (0.495) | non-monotone (0.497) | 97.7 % of rows are exactly 0; it is a pigment-ROI column, and pooled AUC cannot show it |
| `chr_radial_hf` | +1 | raw (0.561), **unblinded** | non-monotone (0.502) | does not hold beyond the pool |
| `chr_solidity` | −1 | x-1 (0.375), unblinded | x-1 (0.422) | agrees; beyond pool 0.410 |
| `chr_concavities_per_um` | +1 | raw (0.621), unblinded | raw (0.605) | agrees; beyond pool 0.603 |
| `dark_concavities_per_um` | +1 | raw (0.651), unblinded | raw (0.631) | agrees; beyond pool 0.627 |
| `h_chr_mean` | +1 | raw (0.761) | raw (0.835) | agrees |
| `chr_resid_frac` | −1 | x-1 (0.373) | x-1 (0.369) | agrees, but not only through pigment (1.1) |

**Columns with direction 0 that have a consistent measured orientation.** Read in-sample only:
- `chr_elongation`: raw, 0.64;
- `chr_hole_frac`: x-1, 0.39 at depth 300;
- `h_chr_lap`: x-1, 0.36;
- `ring_peak_prominence`: x-1, 0.39;
- `ring1_frac_above_tnuc`: raw, 0.64;
- `white_frac_gray215_61`: x-1, 0.40;
- `chr_area_um2`: x-1, 0.43;
- `pair_sep_um`: raw, 0.58, on the half of rows that have a partner;
- `ring1_frac_pale`: x-1 at depth 300 only;
- `chr_area_over_nbr`: x-1 at depth 300 only, on 13 % of rows.

### 1.3 `non-monotone` columns (depth 300) and which kind

`auc_absdev_within_click` is this run's added descriptor. It is the AUC of the distance from the click's median: below 0.5 means mitoses sit in the middle of the range; above 0.5 means they sit at the extremes.
- **No consistent relation, with AUC and absdev both about 0.5:**
  - `pair_area_ratio`, `pair_parallel_deg`, `pair_axis_deg` (the spec expected this: the cue is joint);
  - `brown_frac_surround`;
  - `chr_tophat_core`: 0.461 / 0.474, and the montage mixes networks with homogeneous cores at its low end.
- **Weak sign plus some mid-range concentration:**
  - `chr_radial_hf`: 0.502 / 0.469;
  - `dark_area_frac_101`: 0.464 / 0.464. Beyond darkness it reads 0.374, so mitoses have fewer dark neighbours than equally dark others.
  - `e_h_angle_ring23_deg`: 0.519 / 0.439. Its mitoses sit mid-range, away from red-cell and stroma surrounds.
- **`spec_chr_area_over_nbr`:** no consistent relation.

### 1.4 What I changed in the spec, and why

1. **The search reads only as deep as the top 300 needs (your instruction).**
   - The spec's `detect_deep` ran NMS over every uncapped peak and reported `n_detections_deep`. On some ROIs that is about 625k peaks per click (527, 364: plateau maxima where the clipped channel is flat). It took 18–19 GB per worker and crashed two full-run attempts.
   - The verbatim `nms_by_distance` now runs on a score-sorted prefix: 1,200 peaks, doubled until 300 are kept. No click needed more than 1,200.
   - Greedy NMS keeps a peak only against higher-scored ones, so the 300 rows are exactly those of the full NMS. On 527, all 1,800 candidate rows and every feature were identical to the all-peak run.
   - `n_detections_deep` is no longer computed; `clicks.csv` has `n_nms_prefix` instead.
2. **Gate B gained a third near-tie class.**
   - With only the spec's two example near-ties coded, the first complete full run halted on one row of 294 clicks: `242.tiff` click 1, rank 91.
   - The peak sits at (2279, 3772) here and at (2280, 3772) in the recorded run. This Mac scores those two pixels 1.2e-7 apart (one float32 ulp).
   - The score agrees within 1e-6, and the labels are equal. `od` differs by 2.4e-4 because its window moved one pixel.
   - That is "a difference that a float shift below 1e-4 can explain", the spec's own rule. The added class covers only adjacent-pixel (Chebyshev ≤ 1) peaks with the score within 1e-4 and identical labels, with every other row exact. A label change, a 3 px move or a score difference still halts.
   - Final gate B: 293 exact, 1 near-tie. Gate C: 17,735 shared locations, 0 mismatches.
3. **The neighbour columns use `t_dark`, not `t_nuc` (the spec's defect rule).**
   - On the mandated check (`245`, `462`), and on `300`, `211` and `524`, the kept `t_nuc` components were often cytoplasm and stroma strands, not nuclei. On dense `524` every nucleus merged into an edge-touching sheet and nothing was kept.
   - At `t_dark` the kept components are nuclei, sometimes only the dark part of a pale nucleus.
   - The `t_nuc` versions are kept as `spec_nbr_n`, `spec_nbr_nuc_area_med_um2` and `spec_chr_area_over_nbr`, role `spec`, direction 0, with their own reason column `nan_reason_spec_nbr`.
   - The frozen `directions.csv` rows of the corrected columns still say `t_nuc`; the amendment rows record the change.
   - The choice was made by eye and on a label-free count. No variant's AUC was looked at, although the smoke table had already printed the spec version's AUC on 17 mitoses.
4. **Added outputs (labelled as additions):**
   - `auc_absdev_within_click` (1.3);
   - `auc_wc_beyond_pool` / `auc_dark_beyond_pool` in `validity_depth300.csv`;
   - **tie-robust top/bottom counts.** The spec breaks ties by `cand_rank`. On tie-heavy columns that reports production's tm order, not the column. `brown_frac_surround` read 1,730 vs 409 mitoses (pool) under the spec's rule but 788 vs 749 with ties ordered at random; `chr_hole_frac` and `ring_peak_prominence` are similarly inflated. The tables below show the tie-averaged counts, with the spec's counts beside them where they differ.
5. **Readings where the text allowed more than one.** Each is logged:
   - the debris disc index, k = 0..n_f − 2, the only reading that reproduces the checked values;
   - a third pair-anchor tie-break by the lower label (the columns are symmetric, so it cannot change a value);
   - `M_dark`'s 30 px rule applies to the unfilled component;
   - `matched_ann_id` is −1 when unmatched;
   - "clearly above", "far below" and "low" in the suites were given explicit thresholds, written per row.
6. **Mechanics.**
   - **Montage tie order.** Montage picks break ties by `cx`, then `cy`. On tie-heavy columns, the tiles at the tied end (the zeros of `chr_hole_frac`, `ring_peak_prominence`, both concavity columns and `brown_frac_surround`) are therefore the leftmost locations in the ROI, not a random sample. The judgements describe the tied value correctly, but not a representative set of objects at it.
   - The code is in a worker module, as Stage A did.
   - Montage copies for my own look went to the session scratchpad, not `$TMPDIR`.
   - Montages are embedded as JPEG, not `plt.show()` PNG, to keep the notebook near 20 MB.
   - The full run was repeated after the tie-averaged counts were added. This session's earlier outputs were removed first; no pre-existing file was touched.

---

## 2. The validity tables

- **Columns:**
  - "beyond darkness" is `auc_dark_stratified`;
  - "within size" is `auc_area_tertiles`;
  - "ROIs / domains agreeing" is the share of ROIs, and the number of the 7 domains, on the same side of 0.5;
  - mitoses top 10 / bottom 10 are summed over clicks with ties ordered at random.
- **`validity.csv` population:** accepted rows with `in_recorded_pool`: 7,664 mitoses (412 in failing9) against 21,219 others.
- **`validity_depth300.csv` population:** all accepted rows: 11,089 mitoses (628 in failing9) against 76,815 others.
- **failing9:** its AUCs rest on 9 ROIs and a few hundred mitoses.

#### `validity.csv`

| column | hyp | orientation | AUC within click | beyond darkness | within size | ROIs / domains agreeing | mitoses top 10 / bottom 10 | failing9 / other40 | synthetic |
|---|---|---|---|---|---|---|---|---|---|
| `chr_area_um2` | 0 | x-1 | 0.460 | 0.458 |  | 0.67 / 5 | 809 / 617 | 0.428 / 0.463 | n.a. |
| `chr_elongation` | 0 | raw | 0.636 | 0.627 |  | 0.78 / 7 | 1267 / 626 | 0.731 / 0.627 | pass |
| `chr_solidity` † | −1 | x-1 | 0.375 | 0.361 | 0.371 | 0.88 / 7 | 960 / 362 | 0.292 / 0.383 | pass |
| `chr_concavities_per_um` † | +1 | raw | 0.621 | 0.642 | 0.622 | 0.90 / 7 | 1012 / 348 (spec rule 1012 / 324) | 0.686 / 0.614 | pass |
| `dark_concavities_per_um` † | +1 | raw | 0.651 | 0.674 | 0.649 | 0.94 / 7 | 1184 / 397 (spec rule 1186 / 348) | 0.712 / 0.645 | pass |
| `chr_hole_frac` | 0 | x-1 | 0.435 | 0.460 |  | 0.67 / 6 | 898 / 631 (spec rule 1488 / 632) | 0.481 / 0.430 | pass |
| `chr_radial_hf` † | +1 | raw | 0.561 | 0.593 | 0.571 | 0.82 / 6 | 897 / 565 | 0.604 / 0.557 | fail (predicted) |
| `pair_sep_um` | 0 | raw | 0.581 | 0.624 |  | 0.71 / 6 | 783 / 517 | 0.569 / 0.583 | pass |
| `pair_area_ratio` | 0 | non-monotone | 0.542 | 0.538 |  | 0.65 / 7 | 761 / 610 | 0.566 / 0.538 | pass |
| `pair_parallel_deg` | 0 | non-monotone | 0.517 | 0.532 |  | 0.46 / 4 | 655 / 595 | 0.557 / 0.510 | pass |
| `pair_axis_deg` | 0 | raw | 0.531 | 0.514 |  | 0.69 / 6 | 649 / 542 | 0.555 / 0.527 | pass |
| `h_chr_mean` | +1 | raw | 0.761 | 0.604 |  | 0.96 / 7 | 1673 / 255 | 0.699 / 0.767 | n.a. |
| `h_chr_lap` | 0 | x-1 | 0.397 | 0.422 | 0.397 | 0.71 / 5 | 905 / 455 | 0.369 / 0.400 | pass |
| `ring_peak_prominence` | 0 | x-1 | 0.418 | 0.405 |  | 0.78 / 6 | 893 / 474 (spec rule 1518 / 474) | 0.414 / 0.418 | pass |
| `ring1_frac_above_tnuc` | 0 | raw | 0.634 | 0.594 |  | 0.80 / 7 | 1075 / 402 | 0.694 / 0.629 | n.a. |
| `ring1_frac_pale` | 0 | non-monotone | 0.422 | 0.461 |  | 0.63 / 6 | 985 / 592 | 0.352 / 0.429 | n.a. |
| `white_frac_gray215_61` | 0 | x-1 | 0.363 | 0.347 |  | 0.82 / 7 | 1050 / 354 (spec rule 1181 / 354) | 0.338 / 0.365 | n.a. |
| `e_h_angle_ring23_deg` | −1 | non-monotone | 0.504 | 0.509 |  | 0.49 / 3 | 539 / 587 | 0.563 / 0.499 | n.a. |
| `brown_frac_surround` | −1 | non-monotone | 0.495 | 0.495 |  | 0.29 / 6 | 788 / 749 (spec rule 1730 / 409) | 0.494 / 0.495 | n.a. |
| `dark_area_frac_101` | 0 | non-monotone | 0.466 | 0.407 |  | 0.65 / 5 | 753 / 567 | 0.448 / 0.468 | n.a. |
| `chr_resid_frac` | −1 | x-1 | 0.373 | 0.393 |  | 0.80 / 7 | 1121 / 374 | 0.497 / 0.361 | pass |
| `chr_tophat_edge` † | +1 | non-monotone | 0.460 | 0.517 | 0.462 | 0.57 / 4 | 763 / 595 | 0.497 / 0.456 | pass |
| `chr_tophat_core` | 0 | non-monotone | 0.481 | 0.520 | 0.482 | 0.47 / 4 | 664 / 641 | 0.533 / 0.476 | pass |
| `chr_area_over_nbr` (t_dark) | 0 | non-monotone | 0.457 | 0.494 |  | 0.64 / 4 | 222 / 158 | 0.513 / 0.419 | pass |
| `spec_chr_area_over_nbr` | 0 | non-monotone | 0.465 | 0.481 |  | 0.62 / 5 | 107 / 99 | 0.446 / 0.470 | pass |
| `spec_nbr_n` | 0 | x-1 | 0.473 | 0.474 |  | 0.71 / 7 | 809 / 658 (spec rule 1512 / 532) | 0.450 / 0.475 | pass |
| `spec_nbr_nuc_area_med_um2` | 0 | non-monotone | 0.548 | 0.500 |  | 0.52 / 3 | 108 / 87 | 0.666 / 0.512 | pass |

† Unblinded on 2026-10-09: these rows were seen before `directions.csv` was frozen and are not tests.

#### `validity_depth300.csv`

| column | hyp | orientation | AUC within click | beyond darkness | within size | beyond pool (wc / dark) | ROIs / domains agreeing | mitoses top 10 / bottom 10 | failing9 / other40 | synthetic |
|---|---|---|---|---|---|---|---|---|---|---|
| `chr_area_um2` | 0 | x-1 | 0.429 | 0.448 |  | 0.352 / 0.391 | 0.73 / 7 | 449 / 218 | 0.403 / 0.431 | n.a. |
| `chr_elongation` | 0 | raw | 0.644 | 0.628 |  | 0.695 / 0.675 | 0.84 / 7 | 801 / 295 | 0.724 / 0.638 | pass |
| `chr_solidity` | −1 | x-1 | 0.422 | 0.395 | 0.412 | 0.410 / 0.399 | 0.76 / 7 | 440 / 110 | 0.326 / 0.430 | pass |
| `chr_concavities_per_um` | +1 | raw | 0.605 | 0.635 | 0.611 | 0.603 / 0.630 | 0.92 / 7 | 548 / 162 (spec rule 548 / 109) | 0.650 / 0.601 | pass |
| `dark_concavities_per_um` | +1 | raw | 0.631 | 0.670 | 0.635 | 0.627 / 0.660 | 0.96 / 7 | 602 / 153 (spec rule 603 / 77) | 0.697 / 0.626 | pass |
| `chr_hole_frac` | 0 | x-1 | 0.386 | 0.424 |  | 0.364 / 0.392 | 0.88 / 7 | 510 / 263 (spec rule 1485 / 263) | 0.427 / 0.383 | pass |
| `chr_radial_hf` | +1 | non-monotone | 0.502 | 0.546 | 0.517 | 0.484 / 0.511 | 0.61 / 4 | 323 / 301 | 0.564 / 0.497 | fail (predicted) |
| `pair_sep_um` | 0 | raw | 0.578 | 0.630 |  | 0.547 / 0.603 | 0.73 / 6 | 415 / 238 | 0.557 / 0.581 | pass |
| `pair_area_ratio` | 0 | non-monotone | 0.538 | 0.525 |  | 0.543 / 0.529 | 0.58 / 7 | 376 / 305 | 0.567 / 0.535 | pass |
| `pair_parallel_deg` | 0 | non-monotone | 0.527 | 0.530 |  | 0.519 / 0.524 | 0.54 / 5 | 422 / 372 | 0.551 / 0.524 | pass |
| `pair_axis_deg` | 0 | non-monotone | 0.534 | 0.523 |  | 0.546 / 0.534 | 0.62 / 6 | 456 / 285 | 0.539 / 0.533 | pass |
| `h_chr_mean` | +1 | raw | 0.835 | 0.677 |  | 0.784 / 0.684 | 1.00 / 7 | 1598 / 12 | 0.784 / 0.838 | n.a. |
| `h_chr_lap` | 0 | x-1 | 0.362 | 0.386 | 0.354 | 0.394 / 0.385 | 0.78 / 6 | 574 / 168 | 0.371 / 0.362 | pass |
| `ring_peak_prominence` | 0 | x-1 | 0.393 | 0.379 |  | 0.405 / 0.374 | 0.82 / 7 | 483 / 86 (spec rule 1517 / 86) | 0.397 / 0.393 | pass |
| `ring1_frac_above_tnuc` | 0 | raw | 0.641 | 0.594 |  | 0.602 / 0.585 | 0.90 / 7 | 716 / 137 | 0.684 / 0.638 | n.a. |
| `ring1_frac_pale` | 0 | x-1 | 0.411 | 0.462 |  | 0.445 / 0.469 | 0.78 / 6 | 627 / 237 | 0.356 / 0.415 | n.a. |
| `white_frac_gray215_61` | 0 | x-1 | 0.401 | 0.374 |  | 0.391 / 0.381 | 0.78 / 6 | 621 / 157 (spec rule 964 / 157) | 0.374 / 0.403 | n.a. |
| `e_h_angle_ring23_deg` | −1 | non-monotone | 0.519 | 0.525 |  | 0.508 / 0.523 | 0.59 / 3 | 129 / 244 | 0.588 / 0.514 | n.a. |
| `brown_frac_surround` | −1 | non-monotone | 0.497 | 0.496 |  | 0.499 / 0.499 | 0.49 / 7 | 373 / 345 (spec rule 1730 / 101) | 0.497 / 0.497 | n.a. |
| `dark_area_frac_101` | 0 | non-monotone | 0.464 | 0.374 |  | 0.464 / 0.379 | 0.59 / 6 | 305 / 200 | 0.426 / 0.467 | n.a. |
| `chr_resid_frac` | −1 | x-1 | 0.369 | 0.400 |  | 0.389 / 0.412 | 0.76 / 6 | 817 / 148 | 0.503 / 0.359 | pass |
| `chr_tophat_edge` | +1 | x-1 | 0.389 | 0.454 | 0.384 | 0.438 / 0.450 | 0.78 / 6 | 476 / 151 | 0.469 / 0.383 | pass |
| `chr_tophat_core` | 0 | non-monotone | 0.461 | 0.507 | 0.452 | 0.497 / 0.518 | 0.63 / 4 | 313 / 256 | 0.532 / 0.456 | pass |
| `chr_area_over_nbr` (t_dark) | 0 | x-1 | 0.399 | 0.454 |  | 0.356 / 0.402 | 0.69 / 5 | 235 / 89 | 0.413 / 0.391 | pass |
| `spec_chr_area_over_nbr` | 0 | non-monotone | 0.445 | 0.472 |  | 0.415 / 0.450 | 0.68 / 4 | 103 / 78 | 0.409 / 0.453 | pass |
| `spec_nbr_n` | 0 | x-1 | 0.474 | 0.473 |  | 0.481 / 0.475 | 0.76 / 7 | 390 / 271 (spec rule 1510 / 186) | 0.463 / 0.475 | pass |
| `spec_nbr_nuc_area_med_um2` | 0 | non-monotone | 0.510 | 0.477 |  | 0.457 / 0.428 | 0.50 / 3 | 104 / 116 | 0.645 / 0.479 | pass |

**Per scanner (pooled values are not comparable across scanners).**
- **Boundary columns are weakest on Aperio:** `chr_concavities_per_um` 0.565 at depth 300 against 0.61–0.64 elsewhere.
- **`chr_solidity` is strongest on the Hamamatsu scanners:** 0.369 on XR and 0.337 on S360, against 0.46 on 3D Histech and Aperio.
- **Median `chr_edge_sharpness_per_um`:** 3D Histech 1.38, Aperio 1.39, XR 1.32, S360 1.27 (depth 300). It does not track the per-scanner AUCs in any simple way.
- **Scanner and tissue are confounded:** S360 is 3 ROIs (18 clicks).

**Pigment check (`chr_resid_frac`, within-click AUC, pool / depth 300):**
- `211`: 0.244 / 0.319;
- `289`: 0.326 / 0.360;
- `339`: 0.353 / 0.297.

---

## 3. Per feature

Each entry gives the claim, the evidence, the pathology basis (spec Section 5.0) and the measured orientation (in-sample).

- **`chr_area_um2`:** core area.
  - Evidence: no synthetic test.
  - Measured: mitoses have smaller cores (x-1; 0.43 at depth 300, 0.35 beyond the pool). C7 cuts both ways.
- **`chr_elongation`:** major/minor axis ratio.
  - Evidence: the synthetic suite passes.
  - Measured: mitoses are more elongated (raw; 0.64, beyond the pool 0.70, strongest in failing9 at 0.72).
  - Basis: C4 (metaphase plates). The hypothesis stayed 0 because spindle and stromal nuclei are elongated too.
- **`chr_solidity`:** outline convexity. Unblinded.
  - Evidence: synthetic pass, including the lobulation and deciding cells.
  - Montage: the low end mixes reticulated pale nuclei with ragged figures.
  - Measured: x-1 holds beyond the pool (0.410).
  - Basis: C5.
- **`chr_concavities_per_um`:** defects deeper than 0.25 µm per µm of outline. Unblinded.
  - Evidence: synthetic pass.
  - Montage: notched against smooth outlines; the high end is mostly irregular interphase nuclei.
  - Measured: raw, holding beyond the pool (0.603) and within size (0.611).
  - Basis: C1 and Lashen's "ragged border".
- **`dark_concavities_per_um`:** the same on `M_dark`. Unblinded.
  - Evidence: synthetic pass, and it keeps pale arms that `M_chr` loses.
  - Montage: more mitoses at its high end than `chr_concavities_per_um`.
  - Measured: raw, the highest boundary AUC, holding beyond the pool (0.627) and beyond darkness (0.660).
  - **On uncut rows only**, as Section 5.2 requires for this comparison (`dark_touches_window` False). This was computed after the run, from `candidates.csv`.

    | population | `chr_concavities_per_um` | `dark_concavities_per_um` | lead |
    |---|---|---|---|
    | pool | 0.636 | 0.656 | +0.020 |
    | depth 300 | 0.613 | 0.630 | +0.017 |
    | beyond pool | 0.607 | 0.628 | +0.021 |

    On all rows the lead was +0.024 to +0.030, so the window cuts inflate it a little. 21.5 % of unannotated rows touch the window, against 8.3 % of mitoses. The uncut beyond-pool lead (+0.021) is the only unseen evidence. It matches the 2026-10-09 in-sample +0.020, but no CI was computed.
- **`chr_hole_frac`:** enclosed hole share.
  - Evidence: the membrane test passes.
  - Montage: rings and reticula against solid cores.
  - Measured: x-1, so mitoses have fewer holes (0.386).
  - Basis: C2 and C4. Ring-shaped mitoses exist, but are rare in this pool.
- **`chr_radial_hf`:** see 1.1.
  - Evidence: fails the scale rule, as predicted; passes the deciding cell.
  - Measured: non-monotone at depth 300. The unblinded pool signal does not hold deeper.
- **`pair_*`:** anaphase pairs at `t_dark` in the 101 px window.
  - Evidence: the pair test passes, including oblique pairs.
  - About half the rows have no partner.
  - Measured: `pair_sep_um` raw (0.58: mitoses' partners are farther away); the angle columns and the area ratio show no consistent relation. A joint cue cannot show in marginal AUCs, as the spec said.
- **`h_chr_mean`:** core darkness.
  - Measured: raw, 0.76 pool and 0.84 at depth 300, and 0.60 / 0.68 within `od` strata. It still carries darkness information beyond the 51 px window's `od`.
  - Spearman with `od` is 0.79.
- **`h_chr_lap`:** see 1.1.
  - Measured: strong x-1 (0.36), holding within size tertiles. The montage says thinness and reticulation drive the high end.
- **`ring_peak_prominence`:** see 1.1.
  - Measured: x-1. It measures nearby dark mass.
- **`ring1_frac_above_tnuc`:** share of ring 1 above `t_nuc`.
  - Measured: raw (0.64), so mitoses have a non-pale immediate surround.
  - Correlations: −0.917 with `ring1_frac_pale` (the only pair above 0.8) and −0.16 with `white_frac_gray215_61`.
- **`ring1_frac_pale`:** pale, non-glass ring 1.
  - Montage: halos and pale cytoplasm against nucleoplasm.
  - Measured: x-1 at depth 300 (0.41), but the montage shows mitoses with clear halos at the high end on `300`. It is ROI-dependent, with 6 of 7 domains agreeing.
- **`white_frac_gray215_61`:** glass and lumen share in the 61 px window.
  - Measured: x-1 (0.36 / 0.40). This is the white-surround reward the clip reduced; it still separates.
- **`e_h_angle_ring23_deg`:** see 1.1. No signal; it reads red cells and stroma.
- **`brown_frac_surround`:** brown granules beside the candidate.
  - Montage: correct on `291` and `410`, and the red-cell exclusion holds on `123`.
  - Measured: about 0.50 pooled, because 97.7 % of rows are 0. Its use is ROI-specific.
- **`dark_area_frac_101`:** crowding.
  - Measured: about 0.46 within click, but 0.37–0.41 beyond darkness. Among equally dark candidates, mitoses have fewer dark neighbours. It is non-monotone by the ROI/domain rule.
- **`chr_resid_frac`:** see 1.1. It is pigment at the top, but the signal extends beyond pigment.
- **`chr_tophat_edge`:** see 1.1. It fails its claim and its direction.
- **`chr_tophat_core`:** internal thin structure.
  - Montage: granular cores at the high end; networks and homogeneous cores both at the low end.
  - Measured: no consistent relation.
- **`chr_area_over_nbr` (t_dark):** 13 % coverage.
  - Measured: x-1 at depth 300 (0.40), so mitoses' cores are small relative to their neighbours' dark chromatin.
  - It is unreadable on sparse tissue: `nbr_n` < 3 on 96 % of sarcoma rows and 97 % of breast rows.

---

## 4. Redundancy, synthetic failures, and what could not be checked

**Redundancy:** mean within-click Spearman on the pool population.
- The only pair above 0.8 is `ring1_frac_above_tnuc` against `ring1_frac_pale` (−0.917). The two shares plus glass sum to 1.
- Next highest:
  - `chr_area_um2` against `chr_area_over_nbr`: 0.78;
  - `chr_radial_hf` against `chr_solidity`: −0.67;
  - `h_chr_lap` against `chr_tophat_edge`: 0.60.
- Every feature except `h_chr_mean` has |ρ| ≤ 0.23 with `od`.
- A pairwise matrix cannot show a column that is a function of several others.

**Synthetic suites:** 955 rows. Every checked value in the spec reproduced, including all 54 outline values, the deciding and reporting cells, the pair, membrane, noise, clipped-pixel, resid, top-hat, edge and neighbour tests.
- The one failure is `chr_radial_hf` on the scale rule, on all 8 uniformly scaled shapes; the spec predicted it.
- Two report-only values depend on the random 7 % draw: the included clipped value (3.97 against 3.89) and the support (709 → 368 against 383).

**What could not be checked:**
- **Hole split:** `auc_hole_gt03` is undefined, because no click has both classes with `chr_hole_frac` > 0.3.
- **`n_detections_deep`:** not computed (1.4).
- **Out-of-sample orientation:** none. All 49 ROIs are in-sample; this was accepted by you.
- **Window size for the neighbour columns:** left to you. A 101 px window rarely holds 3 whole nuclei; a wider window would be a new constant chosen on these data.
- **Pathologist review of montages:** not done. Every montage reading here is by eye and not pathologist-confirmed.
- **Reticulated `M_chr`:** whether it is the right mask for pale nuclei is a production-refinement question, outside this spec.
