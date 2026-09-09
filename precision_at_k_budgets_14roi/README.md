# Precision at fixed candidate budgets (K = 10, 20, 30, 50)

Two notebooks:

- `precision_at_k_budgets_14roi.ipynb` -- the full analysis, with verification checks.
- `precision_at_k_budgets_14roi_presentation.ipynb` -- the same computation and tables,
  with the checks and explanatory narrative stripped out.

Both compute precision per ROI at K in {10, 20, 30, 50} across the 14 ROIs of
`images/extra_valid`, one seed per ROI, plus a per-domain rollup. Outputs land in
`../results/` as `precision_at_k_14roi_*.csv`.

## Running it elsewhere

Copy these paths, preserving their layout relative to the repo root:

```
midog_utils/                       # the whole package directory
databases/MIDOG++.json             # not MIDOG++.sqlite -- unused here
images/extra_valid/*.tiff          # 013, 094, 201, 233, 245, 246, 300, 301,
                                    # 402, 403, 459, 460, 529, 548
precision_at_k_budgets_14roi/      # this folder
```

Open the notebook directly (Jupyter starts the kernel in the notebook's own directory,
so the `../` paths inside it resolve correctly).

## Environment

Use `requirements-midog-utils.txt` from the repo root, not `requirements.txt` (that one
pins an unrelated fastai/SlideRunner training stack this code never imports):

```
numpy==1.26.4
opencv-python==4.8.1
scikit-image==0.24.0
scikit-learn==1.5.1
scipy==1.13.1
tifffile==2023.4.12
pandas==2.2.2
```

**numpy must stay below 2.0** -- the pinned opencv/tifffile wheels are built against
numpy's 1.x C API and fail with `AttributeError: _ARRAY_API not found` under numpy 2.x.
Install into a fresh virtualenv/conda env from that file rather than into whatever numpy
is already on the machine.
