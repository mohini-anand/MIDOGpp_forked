# Precision at fixed candidate budgets, ranked by `TM_CCOEFF_NORMED` (K = 10, 20, 30, 50)

One notebook: `precision_at_k_budgets_14roi_normed.ipynb`.

Replicates `../precision_at_k_budgets_14roi/precision_at_k_budgets_14roi.ipynb` exactly, with
one change: the template-matching method is `cv2.TM_CCOEFF_NORMED` instead of `cv2.TM_CCOEFF`,
used for both the response map and the ranking. Same 14 ROIs, same seed per ROI, same channel
(`hematoxylin_od`, unclipped -- D3), same deep-floor extraction (`z = -1.5` in per-ROI
robust-z units), same 7.5 um NMS radius (D7), same `tm_score` ranking key (D5). D2 (no
`tissue_mask`) is also unchanged.

`DECISIONS.md` D1 chose `TM_CCOEFF` over `TM_CCOEFF_NORMED` on a `read_50` comparison (`TM_CCOEFF`
won 49/49 decision-grade cells, median ratio 4.32x) but flagged that comparison as owed a
re-derivation at the newer `recall@K`/precision@K metric family, on the match score alone. This
notebook is that re-derivation, at the same precision@K budgets the reference notebook reports.

Computes precision per ROI at K in {10, 20, 30, 50} across the 14 ROIs of
`images/extra_valid`, one seed per ROI, plus a per-domain rollup, plus a head-to-head against
the reference notebook's already-on-disk `TM_CCOEFF` output (Table C). Outputs land in
`../results/` as `precision_at_k_14roi_normed_*.csv`, distinct from the reference notebook's
`precision_at_k_14roi_*.csv` so neither run overwrites the other; the head-to-head lands as
`precision_at_k_14roi_normed_vs_ccoeff_{per_roi,by_domain}.csv`.

**Verdict:** `TM_CCOEFF` beats `TM_CCOEFF_NORMED` on precision@K at every budget -- 28/28 domain
x budget cells, 54/56 ROI x budget cells (2 ties, 0 losses), pooled ratio 3.55x-4.42x rising
with K. This confirms D1's `read_50` finding at the precision@K metric family too; see Table C
in the notebook for the full breakdown.

**No pool-reproduction gate.** `../precision_at_k_budgets_14roi_chromatin/`'s Gate 0 checks that
its candidate pool exactly reproduces the `TM_CCOEFF` reference run's, because that notebook only
swaps the ranking step *after* an identical search and NMS. This notebook swaps the matching
method itself, so the response map, the deep-floor pool, its robust-z scale and the resulting
ranking are all expected to differ from `../results/precision_at_k_14roi_per_roi.csv` -- that
divergence is the subject of the comparison, not something to gate against. Table C instead
asserts the click itself (`seed_ann_id`, `base_size`, `n_gt_mitotic` -- all upstream of `METHOD`)
reproduced identically, so the precision comparison is attributable to the method alone.

**One extra check `TM_CCOEFF_NORMED` needs that `TM_CCOEFF` doesn't.** `TM_CCOEFF_NORMED` can
return NaN on a zero-variance window; `midog_utils.template_match.fused_response` converts that
to a large negative *finite* sentinel, which its sibling `robust_stats`'s `isfinite` filter does
not catch. This notebook counts and asserts zero such sentinel-contaminated pixels on every ROI
(verified: zero on all 14). Dormant in this run, but worth knowing if this method is reused on
noisier or more artifact-prone images -- `robust_stats` itself has no guard against it.

## Running it elsewhere

Copy these paths, preserving their layout relative to the repo root:

```
midog_utils/                       # the whole package directory
databases/MIDOG++.json             # not MIDOG++.sqlite -- unused here
images/extra_valid/*.tiff          # 013, 094, 201, 233, 245, 246, 300, 301,
                                    # 402, 403, 459, 460, 529, 548
precision_at_k_budgets_14roi_normed/  # this folder
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
