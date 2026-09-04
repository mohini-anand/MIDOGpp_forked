# `images/extra_valid/` — the clean 14-ROI draw

14 MIDOG++ ROIs, **2 per tumour type across all 7 domains**, assembled
2026-09-04. Each file is a byte-identical copy of the same-named file in `images/` (made with APFS
`clonefile` via `cp -c`, so it shares blocks with the original and costs no extra disk). `images/`
itself is untouched and still holds all 23 downloaded ROIs.

This folder exists so that experiments can be pointed at a set that is balanced by domain and known
to be seedable, without re-deriving the filter each time.

## Inclusion criterion

An ROI is here if and only if **both** hold:

1. **`n_mitotic >= 15`** — the number of `category_id == 1` annotations for that image in
   `databases/MIDOG++.json`. This is the repo's pre-existing "decision-grade" bar,
   `tm_variant_sweep.DECISION_MIN_MITOTIC`.
2. **seed pool >= 5** — the pool being `seed_selection.agreement_pool(mitotic)` followed by
   `seed_selection.border_filter(pool, 36, roi_shape)`: the mitotic annotations that every rater
   who saw them called mitotic, restricted to those at least 36 px
   (`FSConfig.patch_size // 2`) from the ROI edge. This is the same pool
   `tm_variant_sweep.draw_seeds` draws from, and matches checks 4 and 5 of
   `new_domain_roi_check.py`. The Otsu foreground filter that `seed_selection.pick_seed` applies
   by default is deliberately **not** part of this bar, for the reason `draw_seeds` gives: it
   sizes a box these experiments do not size.

Together they are the conjunction `tm_variant_report.py` already applies as
`decision_grade & (n_seeds >= 5)`.

## Contents

| file | domain | scanner | n_mitotic | seed pool |
|---|---|---|---:|---:|
| `300.tiff` | canine cutaneous mast cell tumor | Aperio CS2 | 181 | 127 |
| `301.tiff` | canine cutaneous mast cell tumor | Aperio CS2 | 218 | 162 |
| `201.tiff` | canine lung cancer | 3D Histech | 18 | 7 |
| `233.tiff` | canine lung cancer | 3D Histech | 18 | 12 |
| `245.tiff` | canine lymphosarcoma | 3D Histech | 90 | 78 |
| `246.tiff` | canine lymphosarcoma | 3D Histech | 116 | 97 |
| `459.tiff` | canine soft tissue sarcoma | 3D Histech | 131 | 121 |
| `460.tiff` | canine soft tissue sarcoma | 3D Histech | 36 | 31 |
| `013.tiff` | human breast cancer | Hamamatsu XR | 18 | 11 |
| `094.tiff` | human breast cancer | Hamamatsu S360 | 82 | 49 |
| `529.tiff` | human melanoma | Hamamatsu XR | 20 | 13 |
| `548.tiff` | human melanoma | Hamamatsu XR | 239 | 184 |
| `402.tiff` | human neuroendocrine tumor | Hamamatsu XR | 105 | 65 |
| `403.tiff` | human neuroendocrine tumor | Hamamatsu XR | 53 | 36 |

All 14 are `train` split per `datasets_xvalidation.csv`, and none fell back to the
contested agreement tier (`agreement_pool` returned a non-empty unanimous set for every one).

## What was excluded, and why it matters

`202.tiff` is the **only** downloaded ROI that clears the mitosis bar but not the seed bar:
16 mitotic figures, of which just 4 are unanimous and survive the border filter. So the
two criteria are not interchangeable — there are 15 decision-grade ROIs on disk but only
14 valid ones, and `202.tiff` is the whole difference. A future reader comparing this folder
against a "decision-grade" count elsewhere in the repo should expect exactly that off-by-one.

## Provenance

Nine of these were already downloaded before 2026-09-04. Five were added that day — `013`, `233`,
`403`, `460`, `529` — one for each domain that had only a single valid ROI, drawn **uniformly at
random** rather than by picking the densest candidate, which would have biased every depth metric
optimistically. The rule, reproducible as written:

```python
# i indexes the alphabetically sorted domains lacking a second valid ROI
candidates = sorted(valid ROIs of that domain not yet in images/, by file name)
pick = candidates[np.random.default_rng([20260904, i]).integers(len(candidates))]
```

Each was verified after download against `new_domain_roi_check.py`'s battery — pixel dims vs the
JSON, `roi_mpp` in 0.22-0.26, ROI area in 1.9-2.1 mm2 — and the measured counts matched the
prediction made from the annotations alone in all five cases.

## Caveat

`bbox_headroom_end_to_end.decision_grade_rois("images/extra_valid")` returns 14 — but only
because the folder was pre-filtered. That function tests `n_mitotic >= 15` alone and never looks at
seed pools. If an ROI with many mitoses but a thin seed pool is later dropped in here, it will be
returned and this folder will quietly stop being the clean set. The guarantee lives in how the
folder was built, not in anything that re-checks it on read.
