# Bbox tightening on MIDOG++: domain mismatch findings

Date: 2026-08-25
Scope: exploratory, single-image/single-domain-sample tests (`explore_dataset.ipynb`)

## Dataset findings

- MIDOG++ `bbox` values are **not** real tight annotations. Ground truth is a single
  point click per object (SlideRunner `type=1` "spot", one row in
  `Annotations_coordinates` per annotation). Every `bbox` in `MIDOG++.json` is a
  synthetic **50x50px box centered on that point** (verified: `bbox = point ± 25`
  exactly, across samples).
- `Setup.ipynb`'s figshare id map is missing entries for `151.tiff`-`200.tiff`
  (50 human breast cancer images present in `MIDOG++.json` but not downloadable
  via the notebook as-is).

## Reference code (`bbox tuning code reference/`) — bug found

- `detect_and_refine_boundary`'s retry condition checks
  `x_tl_new - x_br_new < 2`. Since the right edge is always >= the left edge,
  this is always true, so the "retry with sharpening" path fires unconditionally
  — it never actually uses the plain (unsharpened) detection result. Comment
  above the line indicates intent was `x_br_new - x_tl_new < 2` (width check).
  Fixed in our adapted version.

## Experiment 1 — 1D cross-profile method, threshold tuning (`002.tiff`, 9 mitotic figs)

- Reference default `edge_relThresh` (0.15-0.25) barely tightens H&E crops —
  boxes stay ~48-50px on a side, because busy tissue keeps the crop's relative
  intensity range high even away from the target nucleus.
- Swept 0.4 / 0.5 / 0.6 on the same 9 objects; **0.5 chosen** as the working
  base threshold (meaningful, consistent tightening without over-shrinking).
  (`debug_pics/bbox_tuning_thresh{0.25_default,0.4,0.5,0.6}.png`)

## Experiment 2 — stability check never passes on H&E

- With base_thresh=0.5, **0/9** objects passed the reference's stability check
  (pixel_pad=0 tight pass vs. pixel_pad=5 re-check, tolerance 2.5px on all of
  x/y/w/h) — all fell back to the first successful attempt.
- Root cause (confirmed with measured diffs): the two passes use different-sized
  crop windows, so the relative threshold is recomputed from different min/max
  each time. In busy H&E tissue this shifts the detected edges by more than
  2.5px even on visually-good detections (e.g. ann 9: w_diff=4.7, otherwise
  clean; ann 6: h grew 49->58 because the wider pass-2 window pulled in a
  neighboring dark structure).

## Experiment 3 — non-convexity failure mode (`ann 7`)

- 1D method only samples one row + one column through the crop's exact center.
  `ann 7`'s chromatin is a branching/leaf-shaped mass; a lobe extending
  down-left of center was never sampled by either profile line, so the
  tightened box (35x29) cut through the object.
- Confirmed visually by overlaying the actual sampled row/column bands
  (`debug_pics/ann7_zoom.png`).

## Experiment 4 — 2D Otsu + connected components (fix for non-convexity)

- Replaced the 1D cross with: Otsu-threshold the whole crop, take the connected
  component under the annotated point.
- `ann 7`: 35x29 -> **47x41**, fully captures the lobe
  (`debug_pics/ann7_otsu_debug.png`, `otsu_cc_grid.png`).
- Ran on same 9 `002.tiff` objects — consistently tighter/more shape-accurate
  boxes than the 1D method, one failure (`ann 20`: no component at the center
  pixel at all).

## Experiment 5 — `ann 20` largest-component fallback (temporary probe)

- Cause: the annotated point sits right on the edge of the Otsu foreground
  blob, not inside it, so the center pixel itself is background.
- Quick patch: if center pixel isn't foreground, use the largest component in
  the crop instead. Gave a plausible box (31x25) for `ann 20`, but this is a
  weak heuristic in general (would grab an unrelated large structure in a
  crowded crop) — flagged as temporary, not adopted as a real fix.
  (`debug_pics/ann20_otsu_fallback.png`)

## Experiment 6 — cross-domain sweep (key finding)

Ran both methods (1D cross, Otsu+largest-CC-fallback) on up to 9 mitotic
figures each from one sample image per remaining domain: canine lung cancer,
canine lymphosarcoma, canine cutaneous mast cell tumor, human neuroendocrine
tumor, canine soft tissue sarcoma, human melanoma.
(`debug_pics/domain_<name>_{cross,otsu}.png`)

- **Both methods work reasonably** on lung cancer, mast cell tumor, soft
  tissue sarcoma, melanoma, neuroendocrine tumor — sparser/more isolated nuclei,
  comparable to `002.tiff`.
- **Both methods largely fail on canine lymphosarcoma.** 1D method: nearly
  every box stays ~49x49 (no signal to lock onto). Otsu+CC: several boxes
  latch onto tissue spanning multiple adjacent nuclei rather than one object
  (e.g. `ann 6245`).
- **Reframed diagnosis**: the failure driver isn't shape (non-convexity) alone
  — it's **cellular density**. Lymphosarcoma is sheets of small, uniformly dark,
  tightly packed lymphocyte nuclei with no local "isolated object vs. lighter
  background" contrast anywhere in the crop, which both an intensity-profile
  method and a per-crop-Otsu method structurally depend on. This is a domain
  property, not fixable by more threshold tuning.

## Literature findings (external search, 2026-08-25)

- Direct precedent: **OMG-Net** (*Communications Biology*, 2024) — SAM-based
  nucleus contour segmentation + small ResNet18 classifier, applied to
  pan-cancer mitotic figure detection in H&E specifically.
- **Aubreville et al.** (MIDOG++ authors) have a newer dataset,
  *Subphase-Labeled Mitotic Dataset* (Scientific Data, 2026 / bioRxiv 2025),
  extending MIDOG++ with segmentation masks — worth checking whether it covers
  our images before building more of our own tightening logic.
- Landscape of alternatives to the 1D/2D threshold approach, roughly tiered:
  - classical, no training: GrabCut, watershed, active contours (same local-
    contrast dependence as our Otsu approach — won't fix the lymphosarcoma case)
  - small pretrained instance-segmentation models: StarDist, CellPose,
    HoVer-Net/HoVer-UNet, NuLite — built specifically for touching/non-convex
    nuclei, likely needed for dense domains like lymphosarcoma
  - promptable foundation models: SAM, and H&E/microscopy-tuned variants
    (μSAM, CellSAM, SAMCell) — take our existing point/box as a prompt directly

## Open items

- Lymphosarcoma (and likely other densely-cellular domains) needs an
  instance-aware method, not more threshold tuning.
- Check whether Aubreville's subphase-labeled dataset already has masks for
  our MIDOG++ images before building further.
- `ann 20`-style "point sits outside its own object's mask" cases need a real
  fix, not the largest-component probe.
