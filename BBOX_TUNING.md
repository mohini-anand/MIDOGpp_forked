# Bounding-box tightening

MIDOG++ ground truth is a single pathologist **point click** per mitotic figure, not a
hand-drawn bounding box. Every `bbox` in `databases/MIDOG++.json` is a synthetic
`point ± 25` square — exactly 50x50px, verified for all 26,286 annotations
(`midog_utils/dataset.py`'s `check_invariants`). That box is rarely centred on the
actual chromatin mass and is almost always larger than the object itself, which makes
it a poor template for anything that matches on appearance (template search, a
classifier crop, etc.).

**Bbox tightening** replaces that synthetic box with the real extent of the nucleus
under the click: threshold a small crop around the point and take the connected
component the click lands in.

See **`bbox_tuning_walkthrough.ipynb`** for a fully worked, executed example — one
annotation walked through every preprocessing step, plus examples across all seven
MIDOG++ tumor-type domains.

## How it works

Implementation: `midog_utils/seed_selection.py`, function `tighten_box_otsu`.

1. Crop a small window around the click — `template_match.BASE_SIZE` (51px, the
   50x50 annotation box made odd so the click lands exactly on the centre pixel).
2. Convert the crop to a structural, "more object → higher value" channel —
   always `channels.to_gray_inverted` (grayscale, then `255 - gray`, since H&E stains
   chromatin dark on a light background) — independent of whichever channel a later
   template search might use.
3. Rescale the crop's own min/max to `0-255` `uint8` (`cv2.normalize`).
4. Otsu-threshold it (`cv2.threshold(..., cv2.THRESH_OTSU)`) — the split point is
   picked per-crop from its own histogram, no hand-tuned threshold.
5. Label connected components (`skimage.measure.label`, 8-connectivity). A different,
   non-touching nucleus elsewhere in the crop becomes its own label and is ignored.
6. Take the component under the click. If the click's pixel is background
   (`label == 0`), **return `None`** rather than guessing — see "What it refuses to
   do" below.
7. Sanity-check the accepted component (`skimage.measure.regionprops`):
   - `min_area=50` — rejects a degenerate sliver.
   - `max_area_frac=0.85` of the window's pixel count — rejects a component that has
     grown implausibly large (bridged into a neighbouring nucleus).
   - `min_solidity=0.5` (filled area / convex-hull area) — rejects a concave,
     dumbbell-shaped merge of two touching nuclei, while still accepting a genuinely
     elongated single mitotic figure (anaphase/telophase chromatin), since that stays
     convex.
8. Return the component's bounding box, `(y0, y1, x0, x1)`, in the crop's local
   coordinates — half-open, i.e. `x1`/`y1` are exclusive.

Any rejection above returns `None`. Nothing is patched over.

### API

```python
from midog_utils import dataset, channels, template_match as tm, seed_selection as ss

images, anns = dataset.load_annotations()
rgb = dataset.load_roi("images/002.tiff")
structural = channels.to_gray_inverted(rgb)

row = dataset.image_annotations(anns, "002.tiff", category_id=1).iloc[0]
patch = tm.read_padded_patch(structural, row.cx, row.cy, tm.BASE_SIZE)
box = ss.tighten_box_otsu(patch)  # (y0, y1, x0, x1) in patch-local coords, or None
```

Two higher-level entry points build on `tighten_box_otsu`:

- **`foreground_filter(df, structural_channel, ...)`** — keep only the rows of an
  annotation frame whose click lands inside (or, with `center_tolerance`, near) its
  own Otsu component. Used as one stage of `pick_seed`'s candidate-pool narrowing,
  alongside pathologist-agreement filtering.
- **`tightened_base_size(structural_channel, cx, cy, ...)`** — the tightened box's
  own (odd) longer side, for `experiment.FSConfig(base_size=...)`. Kept centred on the
  click rather than recentred to the component's centroid: only the *scale* is being
  corrected, not the location.

Both accept `method="binary"` (default, two-class Otsu) or `method="multiotsu"`
(three-class, keeps only the brightest class — tested but not adopted as default; see
`Research Logs/design_choices.md` §7), and `center_tolerance` (default `0`, exact-pixel
match; widening it trades a stricter "is this really the same object" guarantee for
recovering clicks that sit a few pixels outside their own mask — §7a/§7b).

## Why not the earlier approaches

Two methods were tried and superseded before landing on the one above — see
`Research Logs/2026-08-25-bbox-tuning-domain-mismatch.md` for the full experiment log.

- **A reference video-tracking algorithm** (`bbox tuning code reference/`), adapted
  for hydrogel-particle tracking: samples a single horizontal and vertical intensity
  profile through the crop's centre and finds threshold crossings. Reused via
  `explore_dataset.ipynb`'s standalone `tighten_box`/`finding_edges` functions (not
  part of `midog_utils`). Two failure modes drove it out:
  - **Non-convexity**: sampling one row and one column misses off-axis lobes on
    branching/irregular chromatin (a real mitotic figure's shape).
  - **Instability**: re-running detection with a wider padding window shifts the
    computed threshold enough to move the "stable" box by more than a few pixels on
    busy H&E tissue — the reference's own stability check almost never passed.
  - It also carried a dead retry condition (`x_tl_new - x_br_new < 2`, always true
    since the right edge is ≥ the left edge) that made the "retry with sharpening"
    path fire unconditionally — fixed in the adapted copy, but the two failure modes
    above were the reason it was replaced rather than kept.
- **A "largest component" fallback**, when the click's pixel isn't itself foreground:
  probed in `explore_dataset.ipynb`, not adopted. A quick fix for one specific case,
  but a weak heuristic generally — in a crowded crop it can grab an unrelated, larger
  structure instead of the annotated one.

The current 2D Otsu + connected-components approach (`midog_utils.seed_selection`)
fixes non-convexity directly (thresholds the whole 2D crop instead of two 1D
profiles) and needs no stability retry ladder. Full rationale:
`Research Logs/design_choices.md`, sections 1-2, 7, 7a, 7b.

## What it refuses to do

Two situations are excluded from the candidate pool rather than given a guessed box:

- **The click sits outside its own object's mask** — common when a pathologist's
  click lands right on a nucleus's edge rather than solidly inside it. `None`, not a
  fallback to the nearest or largest blob.
- **The accepted component fails a sanity check** — too small, implausibly large, or
  non-convex (a merge of two touching nuclei).

This matters most in densely-cellular domains: on canine lymphosarcoma (sheets of
small, uniformly dark, tightly packed lymphocyte nuclei), both this method and its
predecessors can fail to isolate a single nucleus — a domain property from local
contrast having nowhere to anchor, not something more threshold tuning fixes. See
`bbox_tuning_walkthrough.ipynb`'s per-domain grids and
`Research Logs/2026-08-25-bbox-tuning-domain-mismatch.md`, experiment 6.

## Environment

`midog_utils` imports none of `requirements.txt`'s fastai/SlideRunner training stack —
see `requirements-midog-utils.txt` for its (smaller) pinned environment. The pin that
matters most: **numpy < 2** (the pinned opencv wheel and `tifffile`'s LZW codec path
both need the numpy 1.x C API).
