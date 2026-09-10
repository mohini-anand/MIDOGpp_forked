# Pipeline debug visuals

A place to *look* at the pipeline at each stage, rather than read its numbers out of a CSV. Each
notebook here takes one question, runs the real pipeline on real ROIs, and saves a figure at every
stage into `figures/`. Intermediate data lands in `artifacts/` so the plotting cells can be re-run
without re-doing the expensive `matchTemplate` passes.

These notebooks are **diagnostics, not results**. They reproduce committed runs exactly (each one
gates on that before drawing anything) and then take them apart. Nothing here is sized to decide
anything — `DECISIONS.md` D5's bar for this class of question is a 5-seed sweep with the paired
delta clustered at the ROI, and a two-ROI single-seed walkthrough is not that.

## Notebooks

### `seed_refinement_variants.ipynb`

Why the three `precision_at_k` notebook families disagree. They differ in exactly one thing — how
the seed's search template is cut — and that difference is a **2×2**, not a list:

| | template centred on **the click** | template centred on **the component's bbox centre** |
|---|---|---|
| **ungated** — the *largest* Otsu component in the 51 px window, click never consulted | `ungated_click`<br>`precision_at_k_budgets_14roi_8aug.ipynb`'s inline `largest_cc_box` | `ungated_recentred`<br>**new — the empty cell, written here** |
| **gated** — the component the *click's own pixel* lands in; refuse and redraw if there is none | `gated_click`<br>`seed_selection.tightened_base_size`<br>**production today** (D8 amendment) | `gated_recentred`<br>`seed_selection.tightened_template_box`<br>D8 as originally written, superseded |

Two ROIs, because each isolates a different axis:

- **`245.tiff`** — the containment gate is **inert**: all four rules accept the same annotation
  (6274) at the same size (47 px), so the 2×2 collapses to two runs and the only free variable is a
  **4.30 px template shift**.
- **`403.tiff`** — the only ROI of 14 where the gate **changes the draw**: the ungated rule accepts
  a click sitting on a background pixel and borrows an unrelated component (ann 20334, base 33); the
  gate refuses and redraws onto a real, click-containing component (ann 20356, base 51).

14 stages, all four variants side by side: seed window and Otsu labelling → the four boxes → the
actual `matchTemplate` inputs and their pixel difference → the 8-augmentation bank → fused response,
its difference map and `best_aug` → peaks/NMS → the top 50 bucketed → rank strip → crop gallery →
which annotations were claimed → rank churn → the annotations that change hands.

Five gates run before any figure: the gate is inert on 245 (0) and live on 403 (1); all three
committed 8-aug runs reproduce exactly, **each at its own `SELF_HIT_RADIUS`** (2); that radius
mismatch between the committed notebooks is inert (3); and the frame the figures use is the frame
`evaluate_arms` scored (4).

**Configuration:** 8-augmentation bank (4 rotations × 2 flips × 1 scale), `TM_CCOEFF`,
`hematoxylin_od`, deep floor `z = -1.5`, NMS radius = match radius = 7.5 µm, `seed_index = 0` —
D1/D2/D3/D5/D7 held fixed so the seed rule is the only free variable. Roughly 4 minutes of compute
on a CPU-only machine.

## Running it elsewhere

Copy these paths, preserving their layout relative to the repo root:

```
midog_utils/                                        # the whole package directory
databases/MIDOG++.json                              # not MIDOG++.sqlite -- unused here
images/extra_valid/245.tiff                         # the two ROIs the notebook runs on
images/extra_valid/403.tiff
results/precision_at_k_14roi_8aug_per_roi.csv       # the three committed runs Gate 2 checks against
results/precision_at_k_14roi_gatedseed_8aug_per_roi.csv
results/precision_at_k_14roi_prodseed_8aug_per_roi.csv
pipeline_debug_visuals/                             # this folder
```

Open the notebook directly (Jupyter starts the kernel in the notebook's own directory, so the `../`
paths inside it resolve correctly), or from the repo root:

```
jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=3600 pipeline_debug_visuals/seed_refinement_variants.ipynb
```

## Environment

Use `requirements-midog-utils.txt` from the repo root, not `requirements.txt` (that one pins an
unrelated fastai/SlideRunner training stack this code never imports):

```
numpy==1.26.4
opencv-python==4.8.1
scikit-image==0.24.0
scikit-learn==1.5.1
scipy==1.13.1
tifffile==2023.4.12
pandas==2.2.2
```

Plus `matplotlib` for the figures.

**numpy must stay below 2.0** -- the pinned opencv/tifffile wheels are built against numpy's 1.x C
API and fail with `AttributeError: _ARRAY_API not found` under numpy 2.x. Install into a fresh
virtualenv/conda env from that file rather than into whatever numpy is already on the machine.
