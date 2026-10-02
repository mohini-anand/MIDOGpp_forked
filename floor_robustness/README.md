# Floor robustness across scanners: descriptive checks

**Date:** 2026-10-02.

**Status:** post hoc and descriptive. Nothing is committed. Nothing here re-runs a template search.

This is a companion to `production_hematoxylin_only/multiotsu_pad_arms/tm_floor_definition_49roi_6seed.ipynb`. That
notebook's outputs are `results/tm_floor_definition_49roi_6seed_*.csv`.

## Scanner mapping

Taken from `images/extra_valid/MANIFEST.md` and `testing_set/MANIFEST.md`.

| scanner | tumour types |
|---|---|
| 3D Histech | canine lung, lymphoma, STS |
| Aperio CS2 | canine mast cell; breast `123` and `128` |
| Hamamatsu S360 | breast `072`, `094` and `100` |
| Hamamatsu XR | melanoma, neuroendocrine; breast `007` and `013` |

**Scanner and tumour type coincide for 42 of the 49 ROIs.** Only breast spans scanners. Read every "by scanner"
number as "by scanner and tissue group".

## Files

| file | what it holds |
|---|---|
| `floor_per_roi_by_scanner.csv` | Per ROI: absolute `b_otsu`; its ROI clip share; the white share; the TOTSU / W128 / SFILL floors as multiples of `b_otsu` (from the notebook's `floors.csv`); and the scanner. |
| `floor_ladder_delta_vs_A1_per_roi.csv` | Per ROI and arm, the change against no clip (`A1`): P@10 points and pool mitoses per click (from the notebook's `per_run.csv`). |
| `macenko_vectors.py` → `macenko.csv` | Per-ROI H and E stain vectors by Macenko 2009, Algorithm 1 (OD > 0.15; SVD plane through the origin; 1st/99th-percentile angles). Also their angles to `rgb2hed`'s fixed Ruifrok vectors, and `eosin_leak_ratio`, the H-channel response to one unit of the slide's own E vector relative to its H vector. |
| `eosin_leak.py` → `eosin_leak.csv` | The H level of eosin-rich non-nuclear tissue (top E quartile, below the upper multi-Otsu H threshold) relative to `b_otsu`, and the Spearman correlation of H with E in non-nuclear tissue. |

The per-scanner summary tables and correlations quoted in the chat answer were computed inline from these CSVs.

## Main numbers

**Floor values**
- Absolute `b_otsu` ranges from 0.015 to 0.083 OD across the 49 ROIs. The per-scanner medians are 0.017 (S360), 0.024
  (Aperio), 0.032 (3D Histech) and 0.044 (XR).
- A single fixed floor at the pooled median would land at 0.49–1.87× each ROI's `b_otsu` (5th–95th percentile).
  11 of 49 ROIs would sit below 0.75×.
- Tissue-only Otsu (TOTSU) against whole-ROI `b_otsu` rises with the ROI's white share (Spearman 0.58). It reaches
  1.26–1.36× in glass-rich ROIs (`401`, `343`).
- The seed-surround floor (SFILL) ranges from 0.30 to 1.55× `b_otsu`.

**Macenko vectors against the fixed Ruifrok vectors** (per-scanner medians)
- They differ by 5–15° for H and 13–30° for E.
- The implied eosin → H leakage ratio is 0.24–0.38 (per ROI, 0.02–0.54).
- **Caveat:** Macenko's E vector is the pinkest 1% extreme, not pure eosin, so this ratio is likely an upper estimate.
  The eosin-rich H level also rises with section thickness.

**What predicts each ROI's floor response:** nothing measured here.
- The 1.0× → 1.5× change and the 0.5× → 1.0× change show no correlation above |ρ| 0.26 with:
  - leakage;
  - white share;
  - absolute `b_otsu`;
  - the TOTSU ratio.
- The leakage ρ is 0.03 for both changes.
