# Saliency box refinement vs the production Otsu gate (2026-09-24)

This folder tests the design notes' clickless template refinement ("Mitosis annotation co-pilot pipeline improvement notes.pdf", Stage 1b, *Recommended implementation*) against this repo's production gate (D8, D12). Everything runs with `/Users/mohinianand/anaconda3/bin/python3` from this folder.

| script | what it does | output |
|---|---|---|
| `common.py` | the saliency rule as the notes specify it, production's gate, and production's search on a precomputed channel | |
| `parity.py` | gate: `common.search` must equal `run_production_pipeline`, and reproduce the bbox3way and fillgeom postD11 rows | PARITY OK on 300/128 |
| `geometry.py` | all 2,048 border-filtered agreement-pool annotations on 49 ROIs. Part A: the MIDOG 50 px box. Part B: simulated drawn boxes (+15/50/100 %, ±20 % jitter, 5 draws, two references) | `out/geometry_midog_box.csv`, `out/geometry_boxes.csv.gz` |
| `precision.py` | P@10/20/30, both rank keys: 294 joint clicks plus 87 refusal-stratum clicks, MIDOG-box arms and drawn-box arms | `out/precision_per_run.csv` |
| `precision_followup.py` | all 211 gray-refused clicks, and the notes' fill template (`sal_fill_51`) against the recorded `hem_mask_fill_1um` | `out/precision_followup_per_run.csv` |
| `analyze.py`, `analyze_followup.py` | ROI-level paired contrasts: exact sign-flip test, bootstrap CI over ROIs, Holm correction, W/L/T | `out/contrasts*.csv`, `out/analyze*.txt` |
| `box_vs_click.py` | D12's reversal test: the box-sized window against the 51 px window at the same clicks, with Holm and D5 | `out/box_vs_click.csv` |
| `contact_sheet.py` | visual check of what each method keeps | `out/contact_sheet.png` |

**The saliency rule as implemented:**
- 3-class Otsu (top class) on hematoxylin OD, over the box grown by 1 µm.
- Components below 1 µm² are dropped.
- Each component scores `mean OD × √area × exp(−d²/2σ²)`, with σ = 0.30 × the box diagonal.
- Kept: the winner, plus every component scoring at least 0.45 of the winner.
- `sal` uses that union as it is; `sal1` dilates it by 1 µm.

**How the template is cut:**
- **Square arms** (`sal_*`, `sal1_*`, `prod_*`): production's construction, a raw square of the kept region's bbox side, centred on the bbox.
- **`sal_fill_51`**: the notes' construction. A fixed 51 px window on the click keeps the `sal1` mask as real pixels and fills the rest with the mean of the filled pixels. This is `hem_mask_fill_1um` with only the keep mask swapped.

## Verification

Every check below passed.

- `default_51` (joint clicks, both axes) matches the fillgeom postD11 rows 588/588.
- `prod_51` (chromatin_od) matches the bbox3way postD11 `gray_bbox` rows 147/147, with the same `base_size`.
- `hem_mask_fill_1um`, re-derived, matches the recorded template SHA-1 and TP@10/20/30 on both axes for all 294 clicks.
- `gray_component` equals `tighten_box_otsu` on every pool annotation checked.
- The drawn-box arms equal the geometry pass's jitter-0 boxes.

## Results

### Which object, at the MIDOG box (2,048 annotations)

- **Refusals.** The gray gate refuses 212 (10.4 %). Saliency always returns a template (2 unreadable, 5 below the size floor).
- **Refusals are mostly multi-mass figures.** Saliency keeps ≥2 masses on 150 of the 213 clicks the gate does not refine (212 refused, 1 unreadable): anaphase/telophase clicks that land in the pale gap.
- **Where the gray gate accepts, both methods pick the same object on 1,830 of 1,835.** The 5 that differ are mostly gap clicks where gray grabbed the pale gap and saliency spans both plates. Saliency is right on 3–4 of the 5, judged by eye.
- **Size.** On those same objects, saliency's square is smaller on 97 % (median 0.77×; median side 27 px vs 39). `sal1`'s square is about production's size (35 px).

### Stability under how the box is drawn (Part B)

Each figure gets 15 boxes (3 paddings × 5 jitters).

| | production (D12 gray gate) | saliency (`sal` / `sal1`) |
|---|---|---|
| template-size range per figure, median | 20–28 px | 2 px |
| centre spread per figure, p90 | 1.3–2.2 µm | 0.1–0.2 µm |
| falls back to the drawn box | 9–21 % (the fallback is the box itself, up to 73 px) | 0–0.1 % |

This holds under both box references, including the gray-derived one, which favours production.

### Precision, 294 joint clicks, MIDOG box (tm_score, P@10)

| arm | P@10 | vs `prod_51` |
|---|---|---|
| `default_51` (no refinement) | 0.509 | **+3.8 pts**, Holm 0.001, 32/12/5 ROIs |
| `hem_mask_fill_1um` | 0.525 | +5.5 |
| `sal_fill_51` (notes' recipe) | 0.517 | +4.7 (−0.8 vs `hem_mask_fill_1um`, Holm 0.14) |
| `prod_51` (production) | 0.470 | — |
| `sal1_51` | 0.470 | −0.0 (null on every K and axis) |
| `sal_51` | 0.400 | **−7.0 pts**, Holm < 1e-4, 8/39/2 |

At chromatin_od the order is the same with smaller gaps: `sal_51` −3.0 pts, `sal1_51` −0.3, `default_51` +0.9. P@10 tracks template size almost monotonically. This matches the repo's mass-detector result (`shape_vs_mass_template_terms_49roi_3seed.ipynb`).

### Drawn boxes, 294 joint clicks, tm_score P@10

| box | `prod_box` | `sal1_box` | `sal_box` |
|---|---|---|---|
| +15 % (saliency ref) | 0.421 | **0.463** (+4.2, Holm 2e-4, 34/8/7) | 0.384 |
| +50 % (saliency ref) | 0.445 | 0.463 (+1.8, Holm 0.11) | 0.388 |
| +100 % (saliency ref) | 0.470 | 0.465 (null) | 0.396 |
| +15/50/100 % (gray ref) | 0.447/0.446/0.446 | 0.461/0.462/0.461 (+1.4–1.6, null) | 0.39 |

- **`sal1_box`'s precision win depends on the reference.** Its +4.2 on tight boxes holds only for saliency-derived boxes. Under the gray reference it is +1.4 to +1.6 and null. The stability result above holds under both references.
- **D12's own reversal test** (`box_vs_click.py`): the box-sized window against the 51 px window at the same clicks. D12 says a loss here should change the decision.
  - **tm_score:** `prod_box − prod_51` is −2.3 to −4.9 pts at K=10 on 5 of the 6 box conditions, under both references. D5 is met on all 5, and Holm < 0.05 on 4 (gray +100 %: Holm 0.13). Loose saliency-derived boxes are null.
  - **chromatin_od:** the losses are smaller and mixed. At K=10, −0.3 to −1.8 pts with Holm < 0.05 only on gray +15 %. More conditions clear Holm at K=20/30.
  - **`sal1_box − sal1_51`:** −0.5 to −0.9 pts; D5 is met on only 3 of the 36 tests.
  - **The boxes are simulated.** Whether D12 should keep the box as its window is for the user to decide.
- **The notes' "losing the click costs nothing" holds for saliency here.** It does not hold for the current D12 box path on tight boxes.

### The clicks production refuses (211 gray-refused clicks, 46 ROIs)

Production falls back to `default_51` on these clicks.

- **Fill template (`sal_fill_51` − `default_51`), tm_score:**
  - K=10: +2.3 pts (Holm 0.07, 22/11/13, D5 not met).
  - K=20: +2.2 (Holm 0.003, 29/6/11, D5 met).
  - K=30: +1.8 (Holm 1e-4, 30/5/11, D5 met).
  - At chromatin_od it is null.
- **This gain belongs to the fill, not to saliency choosing better.** On accepted clicks, fill alone adds about 1.7 over `default_51`. Saliency's contribution is that it makes a keep mask possible where both Otsu gates refuse.
- **Square saliency on the same clicks is worse than the fallback:** `sal_51` −2.2, `sal1_51` −1.2 at K=10.

## Conclusion

1. **Saliency is not more accurate than the gray gate at choosing the object.** They agree on 99.7 % of accepted clicks.
2. **Saliency is more reliable in two ways:**
   - it never refuses, and on refused gap clicks it handles anaphase/telophase;
   - its template does not depend on how loosely or off-centre the box is drawn.
3. **As a drop-in replacement for production's square crop, it loses.** The literal recipe (`sal`) costs 7 pts. The dilated version (`sal1`) ties. Shrinking the footprint is what costs precision, and saliency shrinks it more.
4. **It only helps inside the notes' own construction**, a fixed window with fill:
   - On accepted clicks it ties the repo's `hem_mask_fill_1um` (slightly worse).
   - It adds about 2 pts on the ~10 % of clicks the Otsu gates refuse.
   - Expected pooled gain over a hem-mask fill: roughly 0.2–0.4 pts.
5. **The largest effect measured is not saliency.** Production's own gray refinement is 3.8 pts below using no refinement at tm_score (Holm 0.001, 32/12/5, D5 met), and the fill adds about 1.7 on top of `default_51`.
   - This is scoped to the joint-gated reference clicks. The 2026-09-17 audit (F1) found the chromatin_od version of this cost shrinks about 30% on an ungated draw.
   - It bears on D8 at the product's default axis. It is a decision for the user, not a recommendation to revert.
6. **The notes' edge guard is uninformative on tight boxes.** The search extends 1 µm past the box, so it fires on 79 % of +15 % boxes. At the MIDOG box it fires on 6.2 %; the small-mask guard fires on 2.9 %.

## Caveats

- **Boxes are simulated** from a reference extent, which favours whichever method defines it, so both references were run.
- **Real pathologist boxes were not tested.** Neither were two mitoses in one box, or a box that clips a figure.
- **One jitter draw per click in the precision runs.**
- **The joint clicks exclude refusals by construction.** The refusal stratum covers that gap.
- **The fill arm was run at the MIDOG box only, not in drawn boxes.** Its box stability is inferred from the geometry pass.
- **Holm families:**
  - MIDOG-box family: 15 tests per click set × axis.
  - Drawn-box family: 36 tests.
  - Follow-up families: 12 and 9 tests.
