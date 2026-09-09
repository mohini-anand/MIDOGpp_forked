# Precision at fixed candidate budgets, ranked by chromatin density (K = 10, 20, 30, 50)

One notebook: `precision_at_k_budgets_14roi_chromatin.ipynb`.

Replicates `../precision_at_k_budgets_14roi/precision_at_k_budgets_14roi.ipynb` exactly through
NMS (same 14 ROIs, same seed, same search, same deep floor, same 7.5 um NMS) and changes only
the ranking step: instead of one `tm_score` arm, six arms rank the same post-NMS pool --
`tm_score` (D5's production ranker, included as the baseline) plus the five chromatin-density
axes `DECISIONS.md` D5 names as open (`chromatin_od`/`od51`, `od31`, `od_falloff`,
`mask_od_mean`, `od_contrast`).

Computes precision per (ROI, arm) at K in {10, 20, 30, 50} across the 14 ROIs of
`images/extra_valid`, one seed per ROI, plus a per-domain rollup per arm. Outputs land in
`../results/` as `precision_at_k_14roi_chromatin_*.csv`. A reproduction gate at the top of the
notebook (Gate 0) asserts the candidate pool matches the reference notebook's own
`../results/precision_at_k_14roi_per_roi.csv`, so any precision difference is attributable only
to ranking, not to a drifted search or NMS.

This is a **single-seed** run. It reports observations at precision@K, not a decision -- D5's
own bar for promoting a chromatin axis is a 5-seed sweep with the paired delta clustered at the
ROI, which this notebook does not attempt.

## Running it elsewhere

Copy these paths, preserving their layout relative to the repo root:

```
midog_utils/                             # the whole package directory
databases/MIDOG++.json                   # not MIDOG++.sqlite -- unused here
images/extra_valid/*.tiff                # 013, 094, 201, 233, 245, 246, 300, 301,
                                          # 402, 403, 459, 460, 529, 548
results/precision_at_k_14roi_per_roi.csv # the reference run this notebook's Gate 0 checks against
precision_at_k_budgets_14roi_chromatin/  # this folder
```

Open the notebook directly (Jupyter starts the kernel in the notebook's own directory, so the
`../` paths inside it resolve correctly).

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

**numpy must stay below 2.0** -- the pinned opencv/tifffile wheels are built against numpy's 1.x
C API and fail with `AttributeError: _ARRAY_API not found` under numpy 2.x. Install into a fresh
virtualenv/conda env from that file rather than into whatever numpy is already on the machine.
