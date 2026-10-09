# Edit plan for `CANDIDATE_FEATURE_EXTRACTION_PROMPT.md`

**Status:** a plan only, not applied, pending the user's approval. Drawn up 2026-10-09 by two reviewers:
- the main agent;
- a second, independent reviewer (subagent, premise-reviewer type). It wrote its own review before reading the
  main agent's, then two rounds of discussion followed, plus a third round on the spike metric (items 33–35).

All 35 items are agreed between the reviewers (items 33–35 are in the addendum at the end).
- **The baseline is the working tree.** It already holds the main agent's uncommitted wording edits of 2026-10-08
  (`git diff CANDIDATE_FEATURE_EXTRACTION_PROMPT.md`).
  - The reviewers keep those edits, except the `small piece` extension, which item 21 reverts.
- **Evidence for items 1–32 is synthetic, plus a pixel count over the 49 ROIs.** Items 33–35 add an in-sample look at real candidates (addendum). The numbers are quoted
  inline because the scratch scripts live in a session scratchpad.
- **The result:** 36 features → 23 (24 with item 33). The 8 diagnostics (9 with item 33) are `chr_touches_window`, `chr_ns_share`, `nbr_n`, the three
  `h_ring*_mean`, `nbr_nuc_area_med_um2` and `chr_edge_sharpness_per_um`.

## The answer to the review question

**Does the spec measure what pathologists use?** Partly.
- **The cue list is right.** It covers:
  - surface spikes (C1);
  - no membrane (C2);
  - dark chromatin (C3);
  - configuration (C4);
  - smooth round imposters (C5);
  - pigment (C8);
  - cytoplasm (C9).

  The second reviewer re-read VCGP 2021 and the QRG 2022 from the PDFs, and every quote the spec uses is verbatim.

**C1 (spikes) is conditional on how dark the arms are.**
- **The mechanism.** The outline columns read `M_chr`, the window's own upper multi-Otsu class.
  - Arms as dark as the body survive.
  - At 0.7× body OD, `M_chr` keeps a median of about 3 % of the arm pixels.
  - At 0.5×, it keeps none.
- **Why the paler case is likely the real one [inferred].** An arm is one chromosome deep, while the body is many.
- **External support.** Lashen et al. 2021 (*Histopathology*, n = 506, light microscope against WSI): "58% of mitoses
  showed absent hairy-like projection in WSI; however, 89% retained their ragged nuclear border".
  - So the outline columns measure the dense core's ragged border, which is the cue that survives on WSI.
  - `chr_tophat_edge` measures the fine arms. It reads up to 1 µm outside `M_chr`, keeps 25–60 % of its signal at
    0.7×, and is the only C1 column that degrades gracefully.

**C2 (no nuclear membrane) is effectively unmeasured,** and the nucleolus is out of scope. The resting-nucleus side
therefore rests on outline and texture alone.

## The edit list

Section numbers are the prompt's. The status of every item is **agreed**.

| # | section | change | rationale |
|---|---|---|---|
| 1 | 5.2, new subsection after the outline suite | **A mandatory end-to-end C1 test.** Spelled out below the table | The existing suites hand in ideal masks, so they never test what the refinement keeps |
| 2 | 5.0 C1 row; Sources; Independence | **Add Lashen 2021** with the quote above. `claims_to_measure` for the outline columns: "irregularity of the dense core's outline". For `chr_tophat_edge`: "fine projections within 1 µm of the core". Independence paragraph: Lashen is the Rakha group, MIDOG++'s criteria reference 26 (Ibrahim 2022). It is independent of the VCGP authors, not of the dataset's criteria | The one primary study of C1 at WSI scale |
| 3 | 5.2 "What the suite must show" | Change "`chr_radial_hf` and `chr_concavities_per_um` pass" to "pass on ideal masks; the end-to-end test (item 1) shows what survives the production mask". Add: "no synthetic body here is a ragged core; the montages are the only check of that" | The pass is only at mask level |
| 4 | 5.2 lobulation test | **Add `chr_solidity`** with checked values (r = 8 / 16): disc 0.925 / 0.961, touching discs 0.831 / 0.861, notched disc 0.819 / 0.869, 8-lobed star 0.663 / 0.708, 12-lobed star 0.738 / 0.680, spiked ellipse 0.692 / 0.753. Keep −1. Confound: "touching or lobulated nuclei score low, like mitoses" | It has a direction but was never tested for lobulation. It passes |
| 5 | 5, 5.2, 8 | Drop `chr_roughness` | It fails the spec's own lobulation test, has direction 0, and is a convex-hull comparison like solidity |
| 6 | 5, 5.2, 8 | Drop `chr_n_pieces` | It counts specks and detached spike tips: 5 on the tip-spiked bar, 6–8 beside a neighbour. It also depends on the refinement mode, the spec's own reason for moving the pair columns to the 101 px window. The refinement's dilation bridges gaps of up to 2 px. `n_pieces` stays |
| 7 | 5, 5.2 | Drop `chr_minor_um` | It equals √(4·A / (π·elongation)) exactly on ellipses and within 2–8 % on other shapes. A pairwise `redundancy.csv` cannot see that |
| 8 | 5, 5.4, 8 | Drop `h_chr_glcm_contrast`, its constants and its checked values | It reads the edge ramp, so it measures size. Blur alone raises it ×7 in the spec's own noise test (0.66 → 4.85 at r = 16). Through the verbatim refinement, a flat blurred disc reads 10.0 / 5.1 / 2.4 / 1.9 at r = 5 / 8 / 12 / 16. Computed on the eroded support it would duplicate `h_chr_lap` (Spearman 0.90) |
| 9 | 4, 5.6, 8 | Drop `e_over_h_chr`, `M_chr_clip`, `nan_reason_colour`, and the Section 4 sentence that refers to them | It has no counting criterion and no false-positive mechanism, which fails the spec's own pruning rule |
| 10 | 5.5 | Drop `ring1_frac_glass`. One line: the glass share is 1 minus the other two | The three shares sum to 1 |
| 11 | 5.5 | Drop `e_ring23_mean` | The hue cue lives in the angle (item 17) |
| 12 | 5.7 | Drop `nuc_area_frac_101` | `t_nuc` admits cytoplasm, so this reads roughly "not white". No named mechanism |
| 13 | 5.3, 5.5, 8 | `h_ring1_mean`, `h_ring2_mean` and `h_ring3_mean` become role `diagnostic`. `ring_peak_prominence` stays a feature. 5.3's last sentence says the ring means are stored as diagnostics because the prominence is a function of them | An exact functional dependency. This keeps the one column tied to C2 |
| 14 | 5.8(d) | `nbr_nuc_area_med_um2` becomes `diagnostic` | `chr_area_over_nbr` is the ratio over it |
| 15 | 5.0 texture paragraph; 5.4; 5.8(b) | `h_chr_lap` and `chr_tophat_core` move from −1 to direction 0. **Rationale:** the only text behind −1 is one VCGP clause about resting nuclei (item 5.iv). The lap column is measured to depend on size: through the refinement it falls 1.74 / 1.51 / 1.26 / 1.23 at r = 5 / 8 / 12 / 16. **Confound:** a pixel-noise floor may rise with darkness, depending on the noise model. With constant σ_I it goes 2.1 → 3.6 over H 0.08–0.28; with σ_I ∝ √I it stays flat. Do **not** call pyknotic nuclei "homogeneous": no fetched source says so | Weak basis, plus a size confound |
| 16 | 6 | Compute `auc_area_tertiles` for `h_chr_lap` and `chr_tophat_core` too | Both depend on size |
| 17 | 5.5 | `e_over_h_ring23` becomes `e_h_angle_ring23_deg` = degrees(atan2(mean E, mean H)) over the same cytoplasm ring. NaN, reason `no cytoplasm`, below `NS_MIN_PX` pixels or when mean E and mean H are both 0 | It ranks the same as the ratio. skimage 0.24 clamps the stains at 0, so mean H can approach 0 and the ratio would blow up |
| 18 | 5.0 C9; 5.5 | Item 17's column gets direction −1. **Confounds:** "color is not reliable" (VCGP); and a dark clump inside an interphase nucleus has nucleoplasm in its 1–3 µm ring, so it reads as mitosis-like | VCGP 2021 item 4 ("eosinophilic cytoplasm … in degenerating, apoptotic or necrotic cells"); QRG 2022 ("MF tend to have more amphophilic/basophilic cytoplasm … MLF tend to have more eosinophilic"); NTP atlas ("condensed, deeply eosinophilic cytoplasm") |
| 19 | 5.8(c) | Remove "Blurred dark objects are a named false-positive class" | Unsourced. Section 1's false-positive classes don't include it, and MIDOG 2022's out-of-focus figure is about mitoses |
| 20 | 5.8(c), 5, 6, 8 | `chr_edge_sharpness_per_um` becomes `diagnostic`: the covariate for whether a cue is resolvable on a scan. It stays in the montage list, and Section 6 still reports each scanner group's median beside the per-scanner AUCs | With item 19's sentence gone, no criterion or mechanism remains |
| 21 | 5.2 "The piece"; 5.8(b); 8; revision table | **Revert the 2026-10-08 extension.** `chr_tophat_edge` and `chr_tophat_core` are NaN `small piece` when *whole* `M_chr` < 30 px; the 5.2 columns keep the largest-piece rule. Rewrite three places: the 5.2 sentence beginning "This is the one `small piece` rule…" (it must name two rules), the 5.8(b) bullet, and the revision-table clause. The Section 8 `nan_reason_geom` row becomes: "`small piece` (whole `M_chr` < 30 px for `chr_tophat_edge` and `chr_tophat_core`; largest piece < 30 px for the rest)". No new flag: Section 7's NaN share by reason and bucket already gives the count | Clusters of 5–8 debris fragments pass the refinement with whole `M_chr` 66–160 px but a largest piece of 14–21 px. The extension makes them 100 % NaN and hides `chr_tophat_edge`'s highest "other" values (0.32–0.62), which are the small-body confound. On the earlier pool, whole `M_chr` was never below 30 px (minimum 34) |
| 22 | 5.7 | `dark_area_frac_101` excludes `M_sel` (placed as in Section 4) from both numerator and denominator | The candidate's own selection is 2.1–8.4 % of the window, an object-size term inside a crowding column |
| 23 | 4 Labels; 5.1; 8 | Store `matched_unanimous` and `matched_n_mitotic_votes`. Join them on `matched_ann_id = ann_id` from the `unanimous` and `n_mitotic_votes` columns that `ds.load_annotations` already returns; leave them empty when unmatched. No validity column | 25 % of MIDOG++ mitoses (3,020 of 11,937) are 2-to-1 decisions. The downstream analysis computes any split |
| 24 | 5.6; 8 montage list | **The `brown_frac_surround` rule** becomes `R > G > B`, `R − B > 40`, `G − B > 15`, `gray < 180`. **Montage rows:** `291`, `410` and `123`, replacing `509`; `123` checks the red-cell exclusion. "Calibrate it there" becomes "you may change only the four thresholds, logged, with the original kept as `spec_brown_frac_surround`". **Rationale text:** without the new term, the rule fires on red cells (`123`, 0.280 % of pixels) and pink stroma (`128`, 0.416 %) as often as on pigment (`291` 0.417 %, `410` 0.325 %). With it: 0.026 % / 0.041 % against 0.314 % / 0.189 %. The term was chosen after inspecting these ROIs, and it drops 42 % of `410`'s pigment pixels | Measured over all 49 ROIs and checked against crops. `509` has 96 qualifying pixels in the whole ROI, so its montage row would be blank |
| 25 | 5.8(a) | One line: "`n` is the normalised third column of `hed_from_rgb`, so the column is rgb2hed's unclamped D of the mean OD over the OD norm" | It explains why this works where `d_over_h` failed. Checked numerically |
| 26 | 5.5 | Rename `ring1_frac_nuclear` to `ring1_frac_above_tnuc` | `t_nuc` admits cytoplasm (about 45 % of an ROI) |
| 27 | 2, 8 | Pin the `directions.csv` serialisation: `to_csv(index=False, lineterminator='\n')`, UTF-8 | Otherwise "the bytes it serialises to" is not unique |
| 28 | 4 Windows | The 101 px offset uses `ix = int(round(cx))`, `iy = int(round(cy))` | Matches `read_padded_patch` |
| 29 | 5, 5.2, `directions.csv` | **Keep `chr_radial_hf`**, direction +1. `claims_to_measure`: "high-harmonic outline energy of the dense core". `confound`: fails the scale rule on every shape; disc-like when arms are ≤ 0.7× body OD; may overlap `chr_solidity` (read `redundancy.csv` and `auc_area_tertiles`) | It passes the deciding cell of item 1. The only duplication evidence is a designed synthetic family (ρ −0.92 with solidity) |
| 30 | 6 | One sentence: the `pair_*` columns carry the anaphase cue only jointly, so their marginal AUCs are expected near 0.5 | Prevents a misread |
| 31 | 5.0 C2 row; 10 | C2 row: "no column measures C2 directly". Section 10: "the necrosis context (VCGP item 6a) is not measured" | Honest scope |
| 32 | 0, 3.1, 5 header, 6, 8 | **Mechanical updates.** Change the feature count (36 → 23) and every named column set: boundary, spike (which becomes `chr_radial_hf`, `chr_concavities_per_um`), texture and 101 px. In the montage list, remove `chr_roughness` and `h_chr_glcm_contrast` and rename `e_over_h_ring23`. Section 3.1's surround list keeps only the surround *features*: `ring_peak_prominence`, `ring1_frac_above_tnuc`, `ring1_frac_pale`, `e_h_angle_ring23_deg` | Consistency. The main agent chose the shorter 3.1 option over adding diagnostics to `redundancy.csv`, which would add machinery |

### Item 1 in full: the end-to-end C1 test

**Construction:**
- a 201 × 201 canvas at 0.03 OD, with a cell disc of r 30 px at 0.045 OD;
- the arms drawn at the arm OD, then the body at 0.20 OD, as constant values;
- `scipy.ndimage.gaussian_filter` at σ;
- **then** `np.random.default_rng(0).normal(0, 0.003, (201, 201))` added, and the canvas cast to float32;
- verbatim `refine_template_box` at (100, 100), mpp 0.25, with `M_chr` = selection & (class_map == 2).

**Shapes:**
- **Spiked:** a body disc of radius r plus 8 rectangles at angles 2πk/8, covering r − 1 ≤ t ≤ r + L and |o| ≤ w/2.
- **Touching discs:** radius 0.9r, centred at x = ±0.9r.

**The deciding cell** is the only one that can move a direction: w = 4, arm OD 1.0 × body, σ = 1, r ∈ {8, 12},
L ∈ {4, 8}.
- **Pass:** at each (r, L), the spiked shape scores above the touching discs on `chr_radial_hf`,
  `chr_concavities_per_um` and `chr_tophat_edge`, and below them on `chr_solidity`.
- **Checked values, seed 0** (all four columns pass all four cells, and all 16 cells over seeds 0–3):

| r | shape | L | `chr_radial_hf` | `chr_concavities_per_um` | `chr_tophat_edge` | `chr_solidity` |
|---|---|---|---|---|---|---|
| 8 | touching | — | 0.142 | 0.097 | 0.015 | 0.824 |
| 8 | spiked | 4 | 1.446 | 0.313 | 0.100 | 0.747 |
| 8 | spiked | 8 | 3.791 | 0.189 | 0.136 | 0.576 |
| 12 | touching | — | 0.136 | 0.062 | 0.008 | 0.870 |
| 12 | spiked | 4 | 0.852 | 0.260 | 0.058 | 0.727 |
| 12 | spiked | 8 | 2.406 | 0.163 | 0.101 | 0.526 |

**The reporting cells** never move a direction:
- arm OD 0.7× and 0.5× body;
- w = 2;
- σ = 1.5;
- the plain disc;
- the 8-fragment debris cluster.

For each, report the spike-pixel survival in `M_chr` and every outline and top-hat column. A loss there is the
documented property "outline of the dense core", not a defect to fix.

What the reporting cells are expected to show, measured with the noise-first construction:

| arm OD | arm pixels kept in `M_chr` | `chr_concavities_per_um` | `chr_radial_hf` vs a disc | `chr_tophat_edge` |
|---|---|---|---|---|
| 1.0 | 38–94 % (w = 4, 6); 9–31 % (w = 2) | non-zero | 60–220× | full |
| 0.7 | median ~3 %, max 21 % | 0 on 20 of 24 shapes | median ~4× | keeps 25–60 % |
| 0.5 | ~0 | 0 | ~disc | weak |

## Where the two reviewers disagreed, and how it ended

| point | main agent | second reviewer | outcome |
|---|---|---|---|
| Why texture goes to direction 0 | Pyknotic/apoptotic imposters are homogeneous ("textbook") | Unsourced: VCGP says only "condensed", and Hou 2016 says "a few large clumps". Use the single-clause basis plus the measured size dependence | Reviewer's rationale (item 15) |
| Cytoplasm hue direction | 0 ("not reliable") | −1: VCGP item 4, QRG and NTP state the tendency, and texture had got −1 on less | −1 with confounds (item 18) |
| `chr_edge_sharpness_per_um` | Keep as a feature, direction 0 | No criterion left once the unsourced sentence goes, so diagnostic | Diagnostic (item 20) |
| `small piece` rule for the top-hat columns | Extended the largest-piece rule (2026-10-08 edit) | That makes debris clusters 100 % NaN and hides the confound | Reverted (item 21) |
| Expert-vote columns | Store the columns; decline a validity split | Store them, plus `auc_wc_unanimous` | Columns only (item 23) |
| `chr_radial_hf` | Keep and relabel | Drop, or keep and relabel | Keep and relabel (item 29) |
| Red cells and the brown rule | A real risk | "Drop the worry" (9 ROIs checked) | The main agent measured all 49 ROIs and found red cells (`123`) and stroma (`128`) at pigment rates. The reviewer conceded; `G − B > 15` added (item 24) |
| Deciding-cell numbers | Asked for one pinned construction | Round-1 numbers mixed two constructions | Recomputed under one; all pass (item 1) |

**The main agent accepted these corrections from the reviewer:**
- the dilation bridges 2 px gaps, not 1;
- "w + 1 px wide" holds only on axis-aligned spikes;
- the main agent's GLCM numbers were mask-level, about 2× the end-to-end effect;
- most surround columns are justified by measured false-positive mechanisms (white space, halos), not by C9.

**The reviewer withdrew three claims of its own:**
- that w = 3 spikes are physically realistic;
- that the texture noise floor rises with darkness (it depends on the noise model);
- the 1 µm² floor as a rescue for `chr_n_pieces`.

## Limits that remain

- **The arm-to-body OD ratio on real scans is unmeasured.** It decides how much of C1 the outline columns can see.
- **`G − B > 15` was chosen on the same ROIs** that show it working.
- **No real candidate was measured.** `redundancy.csv` and the montages in the full run are the real-data checks.

## Addendum, 2026-10-09: a spike metric (items 33–35, agreed after a third round)

**The question.** Is there a metric that sees chromosome arms paler than the body without responding to debris
clusters or thin rims?

**Who drew this up.** The main agent ran the test, and the second reviewer reviewed it adversarially in round 3: it
re-ran the synthetic deciding cell, and re-analysed the real-data rows without computing any new real-data variant.
The two agree on everything below; the one contested point (item 33, replace or add) ended at "add".

**What was tested.**
- **(a) The existing outline columns on `M_dark`.** `M_dark` is the 4-connected components of `H_w > t_dark` (the
  ROI-wide upper Otsu, `b_otsu_upper`) that overlap `M_chr`, taking the largest piece with holes filled.
- **(b) A new arm counter, defined before any data was seen.** Take the filled piece minus its opening with a 0.75 µm
  disc, and count the residue components that touch the core and reach at least 0.5 µm beyond it.

**How.**
- **Synthetic:** item 1's pinned construction, with `t_dark` fixed at 0.088 OD, plus debris, vesicular,
  touching-disc and neighbour cases.
- **Real:** every distinct location of the recorded A3 pool (27,612 locations; 28,883 accepted rows over 294 clicks,
  7,664 of them mitosis rows), run through the verbatim refinement.
  - AUC is P(mitosis > other) within a click.
  - CIs come from a 4,000-replicate bootstrap over ROIs.
  - The scripts are in the session scratchpad: `spike/`, and the reviewer's `reviewer2/r3_*.py`.

**Calibration on the 7,114 distinct mitosis locations.** As a fraction of the mitosis's mean `M_chr` H:
- the window's upper Otsu threshold, which bounds `M_chr`, sits at a median of 0.72× (q10–q90 0.67–0.79);
- the ROI-wide `t_dark` sits at 0.44× (0.34–0.53).

So `M_chr` drops any arm paler than about 0.7× the body's mean, and `M_dark` keeps arms down to about 0.45×. The
calibration used labelled mitosis locations.

**Real-data results (recorded A3 pool, in-sample):**

| column | within click [95 % CI] | od-stratified | failing9 |
|---|---|---|---|
| `h_chr_mean` (reference) | 0.761 [0.727, 0.794] | 0.604 | 0.699 |
| concavities on `M_chr` | 0.621 [0.595, 0.647] | 0.642 | 0.686 |
| concavities on `M_dark` | 0.651 [0.622, 0.682] | 0.674 | 0.712 |
| `chr_solidity` on `M_chr` (direction −1, so 1 − AUC) | 0.625 | 0.639 | 0.708 |
| solidity on `M_dark` (direction −1, so 1 − AUC) | 0.577 | 0.593 | 0.667 |
| `chr_radial_hf` on `M_chr` | 0.561 [0.530, 0.593] | 0.593 | 0.604 |
| `chr_radial_hf` on `M_dark` | 0.579 | 0.592 | 0.607 |
| `chr_tophat_edge` | 0.460 [0.423, 0.498] | 0.517 | 0.497 |
| arm counter, `arm_frac`, on `M_chr` | 0.493 | 0.535 | 0.566 |
| arm counter, `arm_frac`, on `M_dark` | 0.451 | 0.470 | 0.512 |

The failing9 figures rest on 9 ROIs; read them as rough.

**`M_dark` against `M_chr` for concavities: lead with the uncut rows.**
- On rows where `M_dark` does not touch the 51 px window, the gain is **+0.020 [0.002, 0.040]** (0.635 → 0.655),
  and 0.657 → 0.679 beyond darkness.
- On all rows it is +0.030 [0.013, 0.048]. About 0.013 of that comes from window cuts: a cut outline reads about 0.02
  per µm lower for both labels, and others are cut 2.5× as often (21 % against 8 %).
- The cut inflates the comparison between the masks, not `M_dark`'s own AUC (0.651 on all rows, 0.655 uncut).
- The `M_chr` columns are not inflated by touching. Their touching rows score near 0.5, which dilutes rather than
  inflates.
- Solidity is better on `M_chr` (+0.048 [0.020, 0.076]), and `chr_radial_hf` shows no difference.

**What it means:**
1. **Neither fine-arm metric we tested separates mitoses on this pool.** Both respond to arms on the synthetic shapes.
   - `chr_tophat_edge` is below 0.5 within click, but about 0.5 within od strata. So on this pool it carries about no
     information beyond darkness.
   - In synthetic cases it reads debris clusters at 0.21–0.57 and vesicular rims at 0.32–0.36, against 0.02–0.18 for
     spiked figures. That confound is a hypothesis for the real pool, not a finding; Section 6's hole-split and
     area-tertile AUCs will test it.
   - The arm counter fires on the notched disc (0.466 at r = 8), and on any outline tip sharper than 3 px.
   - This is consistent with Lashen 2021, where fine projections were lost on WSI in 58 % of mitoses. It doesn't show
     fine arms are unmeasurable, since only two definitions were tested.
2. **The outline-raggedness family carries signal beyond darkness.** That is concavities, solidity and, more weakly,
   `chr_radial_hf`: each separates at least as well once darkness is held fixed. This is *consistent with* Lashen's
   "ragged nuclear border" (89 % retained on WSI), which is a pathologist's visual judgement; a defect count on a
   thresholded mask has not been shown to be the same cue.
3. **`M_dark` lets the concavity count see arms down to about 0.45× body.** With synthetic arms at 0.5–0.7×, `M_chr`
   reads 0–0.22 and `M_dark` 0.16–0.41. Its confounds are merged debris clusters (0.20–0.23) and neighbours merging at
   `t_dark`.

**The edits:**
- **33. Add `dark_concavities_per_um` as a new feature column.** `chr_concavities_per_um` stays on `M_chr`, unchanged.
  The second reviewer argued for adding over replacing, and the main agent accepted: replacing would discard the only
  pre-look concavity hypothesis for a post-look choice whose clean gain has a lower bound near 0. Within-click Spearman
  between the two is 0.43, so they are not near-duplicates.
  - **Section 4 definition,** in the 51 px window, on the unclipped `H_w`:
    1. take the 4-connected components of `H_w > t_dark` that overlap `M_chr`, with `t_dark` = `fl['b_otsu_upper']`;
    2. keep the largest 4-connected piece, a tie going to the lower label;
    3. fill its holes with `scipy.ndimage.binary_fill_holes`.
  - **The value** is 5.2's concavity construction, with the same `CONCAVITY_MIN_UM`, on that piece's longest
    `RETR_EXTERNAL` contour as int32.
  - **`nan_reason_dark`** holds `small piece` (fewer than 30 px, including empty) or `convexity error`. A touch is not
    a NaN.
  - **New diagnostic `dark_touches_window`:** true when the *measured piece* reaches the window edge, matching
    `chr_touches_window`.
  - **Column sets.** It joins the boundary, spike and montage sets, and `auc_area_tertiles` applies to it.
  - **`directions.csv` row:**
    - **Direction:** +1.
    - **Rationale:** "`M_chr` is bounded at about 0.72× the body's mean H and `t_dark` at about 0.44×, so arms at
      0.45–0.7× body survive only here. With synthetic arms at 0.7× and 0.5× body, `M_chr` reads 0–0.22 and this
      column 0.16–0.41."
    - **Confound:** "debris clusters merge into one ragged blob (0.20–0.23), and neighbours merge at `t_dark`. A
      window-cut outline reads about 0.02 per µm lower for both labels, and others are cut 2.5× as often, so compare it
      with `chr_concavities_per_um` only on rows where `dark_touches_window` is False. Added after the 2026-10-09 look
      at the recorded A3 pool: its `in_recorded_pool` validity is in-sample (uncut 0.655 against 0.635), and rows with
      `in_recorded_pool` False are its only unseen evidence."
- **33b. Item 1's knock-on.**
  - In the end-to-end test, `t_dark` is fixed at 0.44 × body = **0.088 OD**, not computed from the canvas.
  - **The deciding cell for `dark_concavities_per_um`.** Measured under item 1's pinned construction, whose `M_chr`
    values reproduce item 1's table exactly. It passes all four cells on seeds 0 and 1:

    | r | shape | L | seed 0 (seed 1) |
    |---|---|---|---|
    | 8 | touching | — | 0.091 (0.092) |
    | 8 | spiked | 4 | 0.299 (0.296) |
    | 8 | spiked | 8 | 0.183 (0.183) |
    | 12 | touching | — | 0.060 (0.060) |
    | 12 | spiked | 4 | 0.239 (0.238) |
    | 12 | spiked | 8 | 0.157 (0.157) |

  - **Pale-arm reporting cells:** w 4, σ 1, arms at 0.7× and 0.5× body. Seed 0, `M_chr` → `M_dark`:

    | r | L | arms at 0.7× | arms at 0.5× |
    |---|---|---|---|
    | 8 | 4 | 0.075 → 0.323 | 0.000 → 0.402 |
    | 8 | 8 | 0.075 → 0.196 | 0.000 → 0.218 |
    | 12 | 4 | 0 → 0.256 | 0 → 0.303 |
    | 12 | 8 | 0.052 → 0.168 | 0 → 0.190 |

  - **Debris reporting cell** on `M_dark`: 0.227 (fragment radius 2.5 px, 8 fragments) and 0.203–0.222 (radius 3 px).
    It is NaN where the fragments stay separate at 0.088.
  - Reporting cells never move a direction.
- **34. Do not add the arm counter.**
- **35. `chr_tophat_edge` keeps direction +1.** Dropping it or moving it to 0 after seeing 0.460 would be pruning on
  the outcome.
  - **`claims_to_measure`:** "thin dark structure within 1 µm of the core, as a share of the core's total H: fine
    arms, debris fragments, thin rims".
  - **Confound,** added to 5.8(b)'s: "on the recorded A3 pool its within-click AUC is below 0.5, but about 0.5 within
    od strata (seen before freezing)".

**Caveat: what the look unblinded.** Section 2's log records all of the following.
- **When, and on what.** The 2026-10-09 look at the recorded A3 pool came before `directions.csv` was frozen. It
  covered the spec's primary validity population: accepted rows with `in_recorded_pool`.
- **Which columns.** For `chr_concavities_per_um`, `chr_solidity`, `chr_radial_hf`, `chr_tophat_edge` and
  `dark_concavities_per_um`, `validity.csv` is in-sample. So their `direction_measured` against
  `direction_hypothesis` is not a test. Rows with `in_recorded_pool` False (ranks 100–299) are the only unseen data.
- **What was and wasn't changed.** No direction was changed. `dark_concavities_per_um` was added after the look, and
  the choice of mask per column (concavities only, not solidity or `chr_radial_hf`) was informed by it.
- **What was fixed in advance.** `t_dark` and the arm-counter parameters were set before the run. The 0.44×
  calibration used labelled mitosis locations.

**Count.** With item 33, the plan gives **24 features** (23 + `dark_concavities_per_um`) and 9 diagnostics (the 8
above + `dark_touches_window`).

**Script defects found in round 3.** None changes a headline number:
- the touch flag in `real_pool.py` covered every overlapping `M_dark` component, not the measured piece (item 33
  defines it on the piece);
- the failing9 bootstrap drew all 49 ROIs, so the number of failing ROIs varied by replicate;
- the arm counter's attachment test used 4-adjacency, which is moot under item 34.
