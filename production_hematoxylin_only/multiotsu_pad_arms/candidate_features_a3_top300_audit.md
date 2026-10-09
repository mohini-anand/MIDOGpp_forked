# Audit of the V2 candidate-feature run (d99e93a)

**Audited:** 2026-10-09, against `CANDIDATE_FEATURE_EXTRACTION_PROMPT_V2.md` and the run's findings document
(`candidate_features_a3_top300_findings.md`).

**What was checked:**
- the worker module, every notebook cell and output, the findings document and the committed tables;
- every AUC, recomputed from `candidates.csv` with an independent implementation;
- the verbatim refinement, re-run through the committed worker on all 83,716 accepted locations, to measure the mechanisms
  the run had inferred from montages.

**Reproduce:** `candidate_features_a3_top300_audit.py`, about 70 s on 12 workers. It writes these files to `results/`:

| file | contents |
|---|---|
| `candidate_features_a3_top300_audit_locations.csv.gz` | per-location measurements |
| `candidate_features_a3_top300_audit_reproduce.csv` | the independent reproduction |
| `candidate_features_a3_top300_audit_decomposition.csv` | the depth-300 pair types |
| `candidate_features_a3_top300_audit_stratified.csv` | every AUC under every control, with CIs |
| `candidate_features_a3_top300_audit_mechanisms.csv` | the mechanism tables |
| `candidate_features_a3_top300_audit_checks.json` | the remaining checks |

**Status of the numbers:**
- Everything here is post hoc and in-sample.
- CIs are ROI cluster bootstraps (49 ROIs, 4,000 draws, seed 20261009). The six clicks per ROI re-find the same objects one
  pixel apart, so the ROI is the unit.
- Two re-definitions in Section 4, `tophat_fix` and `eh_territory`, were built after seeing the original columns'
  results. They are hypotheses for held-out data, not fixes.

## Bottom line

**Implementation defects:** none. Every AUC in both validity tables reproduces exactly.

**Two flaws in the validity method distort several conclusions.** Both came from the spec. The agent did not catch either,
although its own output showed the second.
- The pooled depth-300 AUC is dominated by pairs the search score already separates.
- The `od` strata do not hold darkness fixed.

**Verdicts that change under a proper darkness control:**
- **`chr_area_um2`:** "mitoses have smaller cores" becomes "larger cores" in the pool, and no relation beyond it.
- **`chr_area_over_nbr` and `chr_hole_frac`:** their pool signal is mostly darkness.
- **`chr_radial_hf`:** beyond the pool it is inconclusive, not failed.

**Misattributed mechanisms:**
- `chr_tophat_edge`'s x-1 is `M_chr` fragmentation plus pixel noise. The reticulum reading is right; "mitoses have less
  thin structure" is not.
- `h_chr_lap` measures texture. "Thinness" is not supported.

**Too strong:** `e_h_angle_ring23_deg` "no signal". The pooled null hides strong ROI heterogeneity.

**These survive every check:**
- both concavity columns, with dark over chr by +0.02 (CI above 0);
- solidity, elongation and `h_chr_mean`;
- `white_frac_gray215_61` and `ring1_frac_above_tnuc`;
- `chr_resid_frac` being not pigment-specific. My own alternative explanations for it failed.

## 1. What the run got right

**Implementation matches the spec.** The code follows the text, including the NaN-reason rules, for:
- geometry, `M_dark`, the rings and the texture support;
- the surround columns, the brown rule and `chr_resid_frac`;
- top-hat, edge sharpness, and the pair and neighbour columns.

The RGB and H windows are aligned to within 5.4e-8.

**Calculations reproduce:**
- `auc_within_click` (pool and depth 300) and `auc_wc_beyond_pool` for all 24 features match explicit pair counting exactly.
- `chr_resid_frac` and `chr_tophat_edge` reproduce exactly from the re-run.

**Gates and synthetic suites:**
- Gates A, B and C hold.
- The prefix-NMS shortcut gives identical output.
- Every checked synthetic value reproduces, except the scale failure of `chr_radial_hf`, which the spec predicted.

**Honest handling:**
- The unblinding is stated everywhere.
- The tie-averaged counts are a real catch: `brown_frac_surround` read 1,730/409 under the spec's tie-break and 788/749
  with ties ordered at random.
- The `t_nuc` neighbour defect was found by looking at real windows.
- The run called the spec's resid anchor not like-for-like.
- Its 228 rows with `chr_hole_frac` > 0.3 is right under the spec's `float_precision='round_trip'`: three values are
  0.30000000000000004.

**A concern of mine that turned out fine: `M_dark` does contain `M_chr`.**
- The window's upper Otsu sits at a median of 1.35 × `t_dark` for others and 1.60 × for mitoses; the 5th percentiles are
  1.09 and 1.24.
- `M_dark` fails to contain the filled `M_chr` piece on 1.1 % of others and 0.2 % of mitoses.

**Montage readings, checked by eye:**
- `chr_tophat_edge`'s high end is fragmented chromatin networks in pale open nuclei.
- `e_h_angle`'s high end is red cells and pink stroma.

## 2. What went wrong, by cause

### 2A. Spec-design flaws the agent followed without flagging

**(i) The pooled depth-300 AUC is mostly cross pairs.**
- Pairs of a pool mitosis with a beyond-pool other are 54 % of all pairs.
- A ranker over the top 300 does have to order them, but they largely repeat what the search score already separated.
- So the pooled AUC inflates any column that tracks rank. The depth-300 orientations and the montage sort signs inherit
  this.

| column | depth 300 (reported) | pool × pool (17 %) | beyond × beyond (23 %) | pool-mit × beyond-other (54 %) | beyond-mit × pool-other (7 %) | depth 300, within score quintiles |
|---|---|---|---|---|---|---|
| `chr_tophat_edge` | **0.389** | 0.460 | 0.438 | 0.325 | 0.570 | 0.474 |
| `h_chr_mean` | 0.835 | 0.761 | 0.784 | 0.918 | 0.517 | 0.741 |
| `chr_radial_hf` | 0.502 | 0.561 | 0.484 | 0.482 | 0.572 | 0.539 |
| `chr_hole_frac` | 0.386 | 0.435 | 0.364 | 0.377 | 0.417 | 0.411 |
| `h_chr_lap` | 0.362 | 0.397 | 0.394 | 0.327 | 0.453 | 0.409 |

- `chr_tophat_edge`'s cross terms have opposite signs (0.325 and 0.570). That is the signature of a column that tracks
  search rank.
- Its pool and beyond-pool values do not differ significantly: the bootstrap difference is [−0.010, +0.051].

**(ii) `od` quintiles do not hold darkness fixed.**
- The spec predicted darkness columns would sit "nearer 0.5" there. But `h_chr_mean` reads 0.604 in the pool and 0.677 at
  depth 300 inside `od` quintiles.
- Stratifying on `h_chr_mean` itself overcorrects: it is measured on the same mask as the geometry columns (ρ −0.48 with
  area).
- So the table below uses two darkness measures the mask does not set:
  - `dark50`: the mean of the window's darkest 50 px;
  - `disc`: the mean H in a 1 µm disc at the smoothed H maximum.

| column | pool, within click | pool, `dark50` / `disc` | beyond, within click | beyond, `dark50` / `disc` |
|---|---|---|---|---|
| `chr_area_um2` | 0.460 | **0.566 [0.534, 0.599] / 0.544 [0.510, 0.577]** | 0.352 | 0.493 / 0.473 (CIs include 0.5) |
| `chr_area_over_nbr` | 0.457 | 0.515 / 0.525 (CIs include 0.5) | 0.356 | 0.478 / 0.458 |
| `chr_hole_frac` | 0.435 | 0.472 / 0.487 [0.463, 0.511] | 0.364 | 0.411 / 0.436 |
| `chr_radial_hf` | 0.561 | 0.605 / 0.605 | 0.484 | 0.525 [0.489, 0.566] / 0.533 [0.498, 0.575] |
| `h_chr_lap` | 0.397 | 0.370 / 0.431 | 0.394 | 0.349 / 0.424 |
| `chr_tophat_edge` | 0.460 | 0.456 / 0.506 | 0.438 | 0.405 / 0.478 |

**What this changes:**
- **`chr_area_um2`:** the x-1 is a darkness effect. At matched darkness, mitoses' cores are larger in the pool (C7's
  direction against lymphocytes and apoptotic bodies) and show no difference beyond it. Its orientation depends on whether
  the downstream ranker also uses darkness.
- **`chr_area_over_nbr`:** its x-1 vanishes at matched darkness; it does not flip.
- **`chr_hole_frac`:** essentially a darkness proxy in the pool; it stays x-1 beyond.
- **`chr_radial_hf`:** the drop from pool to beyond is real (difference [0.040, 0.113]). But the beyond-pool 0.484 is
  partly darkness coupling: it reads 0.53–0.55 at matched darkness, with CIs touching 0.5. It is inconclusive beyond the
  pool, not failed, and robust in the pool (0.605).

**(iii) Montage orientation is taken from depth 300** (Section 8 of the spec), so it inherits (i).

### 2B. The agent's own judgement

**(a) It printed the darkness leak and did not flag it.** `h_chr_mean` read 0.60/0.68 within `od` strata, against the
spec's own "near 0.5" prediction. Every "beyond darkness" sentence in the findings rests on that column.

**(b) It led with the depth-300 numbers.** `chr_tophat_edge`'s headline 0.389 is not a property of the column within
either population. Within the pool it reads 0.460 [0.422, 0.497]; beyond it, 0.438 [0.397, 0.478].

**(c) `chr_tophat_edge`: the mechanism was inferred by eye. The measurement supports the reticulum, not the conclusion.**

Within-click Spearman at depth 300:

| measure | ρ with `chr_tophat_edge` |
|---|---|
| fragmentation (`chr_frag` = 1 − the largest piece's share of `M_chr`) | +0.40 |
| hole share | +0.35 |
| noise over `h_chr_mean` | +0.29 |
| skeleton branch density | +0.27 |
| `h_chr_mean` | −0.05 |

Pool AUC under controls:

| control | AUC |
|---|---|
| none | 0.460 |
| within `chr_frag` quintiles | 0.524 [0.488, 0.558] |
| H smoothed at σ 0.7 px | 0.487 |

- The x-1 comes from fragmentation and noise, not from mitoses lacking thin structure.
- The reticulum is a tail: the median `M_chr` is one 4-connected piece for both labels.
- Fragmentation carries signal itself: AUC 0.446 in the pool and 0.418 beyond, with mitoses less fragmented.
- "It does not detect fine arms" stands.

**(d) `h_chr_lap`: "thinness plus texture" was not tested, and does not hold.** The x-1 survives every control:

| control | pool | beyond |
|---|---|---|
| area quintiles | 0.395 | 0.370 |
| `chr_frag` | 0.421 | 0.405 |
| `dark50` | 0.370 | 0.349 |
| `disc` | 0.431 | 0.424 |

- The numerator alone (Laplacian SD) reads 0.572 raw and 0.431 at matched darkness, so the ratio is not a division
  artefact.
- It measures relative chromatin texture: smooth condensed chromatin against granular interphase chromatin.
- That is the VCGP resting-nucleus direction, encoded as −1 in V1.

**(e) `e_h_angle_ring23_deg`: "no signal" ignores heterogeneity.**
- The pooled AUC is null: 0.504 [0.463, 0.547].
- ROI AUCs range from 0.15 to 0.83, with an SD of 0.148. A within-click label permutation gives about 0.04; that null is
  understated, because the six clicks share objects, but the gap is large.
- failing9 reads 0.56 (0.59 at depth 300), against the hypothesis.
- The agent's diagnosis, that it reads the neighbourhood's eosin, is right.

**(f) No uncertainty anywhere.** `orientation_recommended` rests on sign counts alone. For example, `chr_tophat_edge` in
the pool is [0.422, 0.497].

**(g) Procedure and reporting.**
- **Gate B was adapted after it halted.** The spec says "stop, report, do not adapt".
  - The substance is defensible: the score difference is 1.5e-7, a float shift.
  - But the new class also exempts that row's `od` difference of 2.4e-4, which is above the 1e-4 tolerance.
- **Gate C's ratio mixes populations.** The notebook reports "17,735 shared (64.2 % of the recorded pool's 27,612)".
  - 15,708 of those are in the recorded pool: 56.9 %, exactly the spec's "about 15,700, 57 %".
  - The other 2,027 are ranks 100–299.
  - The gate itself is fine.

### 2C. Implementation defects

None found.

## 3. Alternative explanations I tested that failed (the run is right here)

**`chr_resid_frac` as a background or support artefact.** The glass is not 255: its G reads 215–230. But none of the
artefact controls moves the AUC:

| version | pool AUC |
|---|---|
| as committed | 0.373 |
| OD against each slide's glass I0 | 0.383 |
| within support-OD-norm quintiles | 0.370 |
| within `chr_ns_share` quintiles | 0.359 |
| `M_chr` only | 0.359 |

Its non-pigment x-1 is a real chromatin-colour difference.

**`h_chr_lap` as a signal-to-noise artefact:** refuted; see 2B (d).

## 4. Fixes

### Method

1. **At depth 300, take orientation from the within-score-quintile AUC,** and use it for the montage sort sign. Report it
   beside the pool × pool and beyond × beyond AUCs.
2. **Control darkness with a mask-independent measure:** `dark50`, the disc mean, or `od` together with one of them.
   Report the marginal and the conditional orientation side by side. A ranker that already uses darkness needs the
   conditional one.
3. **ROI cluster-bootstrap CIs on every AUC.**
   - `orientation_recommended` should also require the CI to exclude 0.5.
   - Add an ROI-heterogeneity statistic, with a null that permutes labels per object rather than per click.
4. **Restore a fragmentation diagnostic** (V2 dropped `chr_n_pieces`): `chr_frag`. Read the top-hat and radial columns
   against it.

### Features

| column | change |
|---|---|
| `chr_area_um2` | Record the orientation at matched darkness (+1 in the pool, null beyond) beside the marginal x-1. |
| `chr_area_over_nbr` | Record it as a darkness proxy; fix its coverage (below). |
| `chr_tophat_edge` | Reclassify as a fragmentation and noise proxy, not a C1 arm measure. `tophat_fix` (smoothed H, around the filled largest piece, over perimeter × darkness) reads 0.471 [0.434, 0.507] in the pool. *Exploratory.* A top-hat does not recover an arm cue at this resolution; `dark_concavities_per_um` is the C1 column that works. |
| `chr_radial_hf` | Keep it, labelled "robust in the pool, inconclusive beyond it at matched darkness". A sub-pixel outline may remove the scale failure (`skimage.measure.find_contours` on lightly smoothed H at the window's upper Otsu level). Untested. |
| `h_chr_lap` | Relabel as relative chromatin texture and drop the thinness caveat. |
| `chr_resid_frac` | Relabel as "chromatin colour off the H–E plane"; pigment is its tail. Flag pigment with an ROI-relative high quantile, and subtract the ROI median before any pooled use (scanner baselines are 0.01–0.21). |
| `e_h_angle_ring23_deg` | Restrict the ring to the candidate's own territory: pixels whose nearest `H > t_dark` component is its own. Exclude red cells. `eh_territory` reads 0.463 [0.424, 0.506] in the pool and 0.447 [0.401, 0.497] beyond; failing9 moves to 0.48. The ROI spread is unchanged (SD 0.146). *Exploratory.* |
| `chr_hole_frac` | Record it as a darkness proxy in the pool. |
| neighbour and pair columns | Label `H > t_dark` once per ROI and select by centroid distance (for example 15 µm). That replaces the 101 px window edge rule, which leaves 87 % of `chr_area_over_nbr` NaN. |
| `ring_peak_prominence` | For an envelope, require a dark ring-2/3 band in most angular sectors, inside the candidate's own `t_nuc` component, so one dark neighbour no longer scores. |
| `chr_edge_sharpness_per_um` | Normalise by a local step across the boundary (inner against outer 1 px band), not the ring-1 mean. Clumps inside nuclei have a tiny ring-1 step, so they score "sharp". |

## 5. Corrected verdicts (in-sample)

| verdict | columns |
|---|---|
| Holds | `chr_concavities_per_um`; `dark_concavities_per_um`, which leads chr on uncut rows by +0.020 [+0.002, +0.039] in the pool and +0.021 [+0.004, +0.040] beyond; `chr_solidity`; `chr_elongation`; `h_chr_mean`; `white_frac_gray215_61`; `ring1_frac_above_tnuc`; `pair_sep_um`; `h_chr_lap` (as texture); `chr_resid_frac` (as chromatin colour, plus pigment) |
| Depends on the darkness covariate | `chr_area_um2`: x-1 marginally, +1 at matched darkness in the pool |
| Robust in the pool, inconclusive beyond | `chr_radial_hf` |
| Mostly a darkness proxy | `chr_area_over_nbr`; `chr_hole_frac` in the pool |
| Fragmentation and noise proxy | `chr_tophat_edge` |
| Heterogeneous across ROIs, pooled null | `e_h_angle_ring23_deg` |
| No consistent relation | `chr_tophat_core`; `pair_area_ratio`; `pair_parallel_deg`; `pair_axis_deg`; `brown_frac_surround` (pooled) |
