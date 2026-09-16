# `images/extra_valid/testing_set/` — the 35-ROI testing set

35 MIDOG++ ROIs, **5 per tumour type across all 7 domains**, assembled 2026-09-16 by
`download_testing_set.py`. None of them was on disk before that day: they are disjoint from the 14
ROIs in `images/extra_valid/` and from everything else in `images/`. Unlike the parent folder, these
are real downloads, not clones of `images/*.tiff`, and `images/` does not hold a copy.

## Inclusion criterion

The same bar as `images/extra_valid/MANIFEST.md` — an ROI is here only if **both** hold:

1. **`n_mitotic >= 15`** — `category_id == 1` annotations in `databases/MIDOG++.json`.
2. **seed pool >= 5** — `seed_selection.agreement_pool(mitotic)` then
   `seed_selection.border_filter(pool, 36, roi_shape)`: mitotic annotations that *every* rater who saw
   them called mitotic (2/2 or 3/3), at least 36 px from the ROI edge. This is stricter than "5
   unanimous annotations" — the border filter can only remove figures — and no ROI here relied on
   the contested 2-of-3 fallback tier.

Both were predicted from the annotations before download and re-measured on the real pixel shape
after it, together with the rest of `new_domain_roi_check.check_roi`'s battery: pixel dims equal the
JSON's width/height (the truncation guard), `roi_mpp` in 0.22-0.26, ROI area in 1.9-2.1 mm2, and
the tumour type. `6_tissue_mask_covers_mitotic` was recorded but not gated on (it tests
`baselines.tissue_mask`, not the ROI); see the `checks` column.

## Contents

| file | domain | scanner | split | n_mitotic | n_unanimous | seed pool | checks |
|---|---|---|---|---:|---:|---:|---|
| `322.tiff` | canine cutaneous mast cell tumor | Aperio CS2 | train | 23 | 23 | 23 | ALL |
| `337.tiff` | canine cutaneous mast cell tumor | Aperio CS2 | train | 30 | 13 | 13 | ALL |
| `339.tiff` | canine cutaneous mast cell tumor | Aperio CS2 | test | 15 | 10 | 10 | ALL |
| `343.tiff` | canine cutaneous mast cell tumor | Aperio CS2 | train | 37 | 28 | 28 | ALL |
| `344.tiff` | canine cutaneous mast cell tumor | Aperio CS2 | test | 44 | 30 | 30 | ALL |
| `211.tiff` | canine lung cancer | 3D Histech | train | 20 | 14 | 14 | ALL |
| `220.tiff` | canine lung cancer | 3D Histech | test | 77 | 51 | 51 | ALL |
| `230.tiff` | canine lung cancer | 3D Histech | test | 18 | 9 | 9 | ALL |
| `235.tiff` | canine lung cancer | 3D Histech | train | 22 | 20 | 20 | ALL |
| `242.tiff` | canine lung cancer | 3D Histech | train | 27 | 23 | 23 | ALL |
| `249.tiff` | canine lymphosarcoma | 3D Histech | test | 112 | 66 | 66 | ALL |
| `284.tiff` | canine lymphosarcoma | 3D Histech | train | 75 | 56 | 55 | ALL |
| `289.tiff` | canine lymphosarcoma | 3D Histech | train | 18 | 12 | 12 | ALL |
| `291.tiff` | canine lymphosarcoma | 3D Histech | test | 87 | 61 | 59 | ALL |
| `293.tiff` | canine lymphosarcoma | 3D Histech | train | 105 | 79 | 79 | ALL |
| `410.tiff` | canine soft tissue sarcoma | 3D Histech | train | 17 | 13 | 13 | ALL |
| `425.tiff` | canine soft tissue sarcoma | 3D Histech | test | 17 | 11 | 11 | ALL |
| `435.tiff` | canine soft tissue sarcoma | 3D Histech | train | 53 | 40 | 40 | ALL |
| `462.tiff` | canine soft tissue sarcoma | 3D Histech | train | 15 | 13 | 13 | ALL |
| `499.tiff` | canine soft tissue sarcoma | 3D Histech | train | 62 | 60 | 59 | ALL |
| `007.tiff` | human breast cancer | Hamamatsu XR | test | 22 | 17 | 17 | ALL |
| `072.tiff` | human breast cancer | Hamamatsu S360 | train | 49 | 32 | 32 | ALL |
| `100.tiff` | human breast cancer | Hamamatsu S360 | test | 44 | 22 | 22 | ALL |
| `123.tiff` | human breast cancer | Aperio CS2 | test | 36 | 27 | 27 | ALL |
| `128.tiff` | human breast cancer | Aperio CS2 | test | 30 | 28 | 28 | ALL |
| `509.tiff` | human melanoma | Hamamatsu XR | train | 18 | 18 | 18 | ALL |
| `524.tiff` | human melanoma | Hamamatsu XR | train | 15 | 9 | 9 | ALL |
| `527.tiff` | human melanoma | Hamamatsu XR | train | 43 | 29 | 29 | ALL |
| `546.tiff` | human melanoma | Hamamatsu XR | train | 28 | 15 | 15 | ALL |
| `549.tiff` | human melanoma | Hamamatsu XR | train | 28 | 22 | 22 | ALL |
| `353.tiff` | human neuroendocrine tumor | Hamamatsu XR | test | 63 | 50 | 47 | ALL |
| `364.tiff` | human neuroendocrine tumor | Hamamatsu XR | train | 23 | 7 | 7 | ALL |
| `400.tiff` | human neuroendocrine tumor | Hamamatsu XR | train | 104 | 95 | 91 | ALL |
| `401.tiff` | human neuroendocrine tumor | Hamamatsu XR | train | 19 | 14 | 14 | ALL |
| `404.tiff` | human neuroendocrine tumor | Hamamatsu XR | train | 61 | 51 | 49 | ALL |

`split` is the MIDOG++ cross-validation split from `datasets_xvalidation.csv`, recorded for
reference only; the draw did not filter on it.

## Draw rule

Uniformly at random within each domain, never by density — picking the densest candidates would bias
every depth metric optimistically. Reproducible as written:

```python
# i indexes the alphabetically sorted domains
candidates = sorted(valid ROIs of that domain not in images/ or images/extra_valid/, by file name)
order = candidates[np.random.default_rng([20260916, i]).permutation(len(candidates))]
# walk `order`, keeping each ROI that passes verification, until 5 are kept
```

The full ordering for every domain, with figshare ids, is in `results/testing_set_draw_order.csv`, written before the first
download. Every attempt, including any rejection, is in `results/testing_set_rois.csv`.

| domain | candidates | ranks tried | rejected |
|---|---:|---|---|
| canine cutaneous mast cell tumor | 23 | 0-4 | none |
| canine lung cancer | 19 | 0-4 | none |
| canine lymphosarcoma | 49 | 0-4 | none |
| canine soft tissue sarcoma | 22 | 0-4 | none |
| human breast cancer | 33 | 0-4 | none |
| human melanoma | 18 | 0-4 | none |
| human neuroendocrine tumor | 6 | 0-4 | none |

## Caveats

- **Human neuroendocrine tumor had only 6 valid candidates for 5 slots**, so this domain is
  close to the whole valid population not already in `images/extra_valid/`, not a sample from it:
  only `356.tiff` (rank 5: 18 mitotic, seed pool 11) was left out. It also holds the thinnest seed
  pool in the set, `364.tiff` at 7, two above the bar.
- The parent folder's caveat applies verbatim: `bbox_headroom_end_to_end.decision_grade_rois()` tests
  `n_mitotic >= 15` alone and never looks at seed pools, so the guarantee lives in how this folder
  was built, not in anything that re-checks it on read.
- Existing scripts that read `images/extra_valid/` use a non-recursive `*.tiff` listing, so this
  subfolder is invisible to them; the 14-ROI set is unchanged. A recursive glob (`rglob`, `**`)
  over `images/extra_valid/` would silently mix the two sets.
