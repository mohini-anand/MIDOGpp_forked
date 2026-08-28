# Design choices — next iteration

Decisions for the next `find_and_suppress` pass on MIDOG++, following the single-pass
run in `find_and_suppress_midog.ipynb` and the discussion that followed it. Supersedes
nothing yet run; documents intent before implementation.

## 1. Seed selection: pathologist agreement

- **Rule**: seed with a category-1 annotation all 3 raters voted mitotic. If none exist
  in the image, fall back to annotations 2/3 raters voted mitotic.
- **Why**: a seed the experts themselves disagreed on is plausibly a morphologically
  atypical example; using it as the *only* template risks building a poor template from
  the start.
- **Caveat**: unverified effect size — `recall_by_agreement` was degenerate in the last
  run's full-detection-list config (recall 1.0 for both unanimous and contested
  annotations everywhere), so agreement hasn't yet been shown to matter at a budget.
  Combined with the border-distance filter and the foreground filter below, the
  candidate pool can shrink to zero on sparse domains (e.g. 350.tiff has 4 mitotic
  annotations total) — needs an explicit "no candidates left" fallback, not a silent
  empty selection.

## 2. Template: bbox tightening (Otsu + connected component)

- **Method**: Otsu-threshold the padded crop, take the connected component containing
  the annotated center point, use its extent as the template. Scale is dropped —
  matching uses each tightened box's native size, not a fixed 51px canonical size.
- **Center-not-in-foreground handling**: when the click doesn't land inside a
  foreground component, exclude that annotation from the seed candidate pool rather
  than falling back to "largest component in the crop." Applied as a seed-selection
  filter, alongside the existing border-distance exclusion in `pick_seed`.
- **Caveat**: this resolves the "click sits outside its own object's mask" failure mode,
  not the separate "component spans multiple adjacent nuclei" failure mode seen on
  dense domains (documented on lymphosarcoma). A component-size/shape sanity check
  (comparable to `nucleus_blobs`' `min_area`/`max_area`) is still needed for *accepted*
  components, since an oversized multi-nucleus component can still pass the
  point-in-foreground check.
- **Caveat**: variable template sizes reintroduce template-size-dependent score bias —
  `TM_CCOEFF_NORMED`'s null variance shrinks as template pixel count grows, so scores
  aren't directly comparable across differently-sized templates. Relevant if scores are
  ever compared across seeds or images, not within one seed's own ranked list.
- **Padding**: add rotation-safety margin (~box diagonal, i.e. side x sqrt(2)) after the
  box is tightened, sized to the tightened box itself — not to the crop used to compute
  the tightening. Widening the tightening-input crop is the documented cause of unstable
  boxes in the reference method (edges shift with the recomputed min/max).
- **Component size/shape sanity check — implemented**, resolving the caveat above.
  `seed_selection.tighten_box_otsu` now also rejects the accepted component via
  `skimage.measure.regionprops` on three bounds: `min_area=50` (a degenerate sliver),
  `max_area_frac=0.85` of the 51px window's own pixel count (grown implausibly large),
  and `min_solidity=0.5` (filled area / convex-hull area — catches a concave,
  dumbbell-shaped merge of two touching nuclei). Same exclusion path as the
  centre-not-foreground case: the annotation is dropped from the seed candidate pool,
  not patched over.
  - *Why solidity, not eccentricity*: a real mitotic figure can legitimately be
    irregular or elongated (anaphase/telophase chromatin) — that reads as high
    eccentricity on a single genuine object. Solidity doesn't flag a convex ellipse
    like that; it only flags the non-convex bridge a merge produces.
  - *Calibration*: measured on 452 centre-pixel-accepted components across the 10
    downloaded ROIs' unanimous candidates (dense and sparse domains both — see
    `seed_selection_midog.ipynb`, section "2b"). The defaults exclude 5/452 (1.1%): 3
    degenerate slivers (area 4–17px²) and 2 non-convex merges (solidity 0.48–0.49) — a
    light backstop, not a narrowing of the population, matching `nucleus_blobs`'s own
    loose `min_area`/`max_area`. `max_area_frac` doesn't bind on this sample (largest
    observed component: 2102px², under the 2211px² cutoff) but is kept as the same
    kind of defensive ceiling.
  - *Still open*: calibrated on the actual `agreement_pool` population from these 10
    ROIs (unanimous tier, or the contested fallback on a flagged domain) -- not
    unanimous-tier only, as this bullet previously said. The 452-candidate count is
    exactly 448 unanimous-tier candidates across 9 domains plus 350.tiff's 4 contested-
    tier candidates (350.tiff has zero unanimous annotations, hence "flagged"); dropping
    350.tiff's 4 gives 448, confirming the original count already included the contested
    tier. Rechecking against the rest of the (non-downloaded) dataset is still open.

## 3. Search channel: hematoxylin and raw RGB variants

- Run both the hematoxylin channel (`FSConfig(channel='hematoxylin')`, implemented,
  not yet run) and raw RGB (3-channel `TM_CCOEFF_NORMED`) as variants alongside the
  current `gray_inverted` default.
- **Not a variant**: `gray_inverted` vs. non-inverted grayscale — mathematically
  identical under `TM_CCOEFF_NORMED`, which is invariant to a shared `255-x` transform
  applied to both template and search image. Not worth testing separately.
- **Caveat**: hematoxylin's plausible benefit is cross-scanner consistency, not
  mitosis-vs-lookalike specificity — the `nucleus_blobs` baseline already matches in
  hematoxylin space and still can't separate the two classes. Raw RGB risks the
  opposite problem: baking in scanner-specific stain-color variation directly into the
  match.

## 4. Augmentation footprint

- Run two additional variants alongside the current 24-augmentation default (12 angles x
  2 flips): no augmentation (1 template), and 90-degree rotations x flips (8 templates).
- Why: fusing more augmented maps via element-wise max inflates spurious high scores
  under noise (same mechanism as the rejected multi-scale max) — fewer augmentations
  should reduce that, at the cost of recall on off-angle mitoses.
- Compare variants using the same seed per image across all three configs, not a fresh
  random seed per variant, so the comparison isolates augmentation count.

## 5. Evaluation label renaming

- Rename detection buckets for clarity: `TP` -> `human_correct_label`,
  `FP_lookalike` -> `human_rejected_label`, `FP_unannotated` -> `non_human_findings`.
- Scope: the bucket labels themselves (`evaluate.py` bucket constants, `viz.py` colour
  mapping), not the underlying TP/FP counting logic, which is unchanged.

## 6. Simplified pipeline defaults for the next run — implemented

Decided for the next run specifically: the simplest configuration, least moving parts,
not a claim that it is the best-scoring one.

- **Score threshold**: `FSConfig.score_threshold` raised from 0.25 to **0.5**.
  `template_match.extract_peaks`'s own default raised to match. `evaluate.threshold_sweep`'s
  default sweep range trimmed to `(0.5, 0.6, 0.7, 0.8, 0.9)` — the two lower entries
  (0.25, 0.4) removed, since they described detections the search no longer reports at
  the new floor. **Why 0.25 was too low**: it let the detection list tile most of the
  ROI (`coverage_frac` 0.88–0.97 in the first single-pass run), the point past which
  full-list recall stops being evidence about the detector rather than about list
  length (`evaluate.py`'s own module docstring).
- **Augmentation footprint**: `FSConfig.n_angles`/`flips` defaults changed from
  `12`/`(False, True)` (24 templates) to `1`/`(False,)` (1 template) — i.e. the
  `no_augmentation` variant from section 4 is now the pipeline default, not just one of
  three compared options. `experiment.AUGMENTATION_VARIANTS`'s first entry renamed from
  `default_12angles_2flips` to `angles12_flips2` accordingly, since it is no longer the
  default.
- **Multi-scale**: already `scales=(1.0,)` — no change needed; single-scale was already
  the default per section 2's fusion-comparison caveat.
- **Caveat**: this is a deliberate simplification, not a validated improvement. The one
  measurement so far (`augmentation_variants` on 002.tiff ann 17) ranks all three the
  *opposite* of simplicity: `rot90_4angles_2flips` scores highest (`discrimination`
  0.068), the old 12-angle/2-flip default second (0.037), and the new `no_augmentation`
  default worst and negative (-0.034 — i.e. on this one seed it discriminates
  look-alikes from ordinary nuclei slightly *better* than it discriminates mitoses).
  Revisit once this runs across more seeds/domains rather than trusting a single-seed
  read either way.

## 7. Multi-level Otsu for bbox tightening — tested, NOT adopted as default

- **Method**: `seed_selection.tighten_box_otsu` gained a `method` parameter --
  `"binary"` (default, unchanged) or `"multiotsu"` (3-class
  `skimage.filters.threshold_multiotsu`, keeping only the top/brightest class as
  foreground before labelling connected components). Downstream checks (centre pixel,
  `min_area`/`max_area_frac`/`min_solidity`) unchanged either way. Threaded through
  unchanged as an opt-in argument on `foreground_filter`, `pick_seed`,
  `tightened_base_size`, and `experiment.run_one_image(bbox_method=...)` -- every
  caller still defaults to `"binary"`, so nothing about the pipeline's behaviour changed
  by adding this.
- **Why it was proposed**: binary Otsu under-separates on some annotations — e.g.
  350.tiff ann 18161 (`seed_selection_midog.ipynb`, "What the picked seed looks like"):
  the accepted component bridges into a neighbouring structure and the CC bbox comes out
  47x51, barely smaller than the 51px window itself, versus a clean 37px result on
  002.tiff ann 17 from the same method on the same window size.
- **Result: the motivating example itself doesn't clearly improve.** Re-measured under
  this section's own calibration (452 agreement-pool candidates, same 10 ROIs as
  section 2b): 350.tiff ann 18161 goes from a 51px binary box (solidity 0.770) to a
  **50px** multiotsu box (solidity **0.679**) — one pixel smaller, and *less* convex, not
  more. This is the single result that most directly argues against adopting multiotsu:
  the case the whole section was written to fix isn't actually fixed by it.
- **Aggregate numbers, as supporting evidence**:
  - *Exclusion/pool size*: of 452 agreement-pool candidates, binary's end-to-end
    survivor count is 447/452 (98.9% of centre-foreground candidates, 5 excluded by the
    size/shape check — unchanged from section 2b). Multiotsu's is 413/452 — a *lower*
    centre-foreground rate (415/452 vs 452/452) drives nearly all of the difference; the
    size/shape check itself barely moves (2 excluded vs 5). I.e. multiotsu's cost shows
    up almost entirely as "the click no longer lands in a foreground pixel at all," not
    as more merges caught.
  - *Per domain*: multiotsu is never less aggressive than binary in any of the 10 ROIs;
    the extra loss is concentrated in the three densest domains (245/246/300/301.tiff,
    roughly 6–11 points of survival-rate drop each); 002.tiff, 350.tiff, and 506.tiff are
    essentially unaffected.
  - *Size shift, paired correctly*: comparing only the 411 candidates **both** methods
    keep end-to-end (not the full, differently-sized survivor sets, which would
    understate the shrink by comparing populations multiotsu has already thinned) —
    median tightened size drops 36px -> 26px, and 93.4% of paired candidates get
    strictly smaller. Restricted further to binary boxes already near-saturated at the
    window edge (>=48px, the actual target failure mode): 94.7% of the ones multiotsu
    still keeps get meaningfully de-saturated (<=45px). Where multiotsu doesn't drop the
    candidate outright, the tightening effect is real.
  - *Thin-pool domain*: 350.tiff's 4 contested-tier candidates (see the section-2 fix
    above) all still survive under multiotsu — no new empty-pool risk there.
  - *Spot-check on large, already-solid real objects* (the check this section
    specifically asked for): of 187 candidates with a large (>=800px) and already-solid
    (solidity >= 0.5, i.e. binary already reads them as one clean object) binary
    component, multiotsu drops 20 (10.7%) to no-foreground-at-all, including several
    with high binary solidity (up to 0.96 — clearly not bridging cases). Splitting this
    group by whether binary's box was already near-saturated (the target failure mode,
    11.8% lost) or was already small and clean despite a large chromatin area (8.3%
    lost) gives comparable loss rates either way — multiotsu isn't losing only the
    suspect boxes, it's losing clean ones too, at a similar rate.
  - *Are those losses genuine fragmentation, or a brittle centre-pixel test?* For the 37
    candidates binary keeps but multiotsu drops outright, checked whether a multiotsu
    foreground pixel exists in a small neighbourhood of the click (not just its exact,
    possibly-off-by-a-few-px rounded pixel): within 1px, 19/37 (51%); within 2px, 26/37
    (70%); within 3px, 31/37 (84%). A majority of the "losses" are near-misses, not
    genuinely absent foreground — some real cost here is the centre-pixel test itself
    getting brittle under multiotsu's tighter boundary, not multiotsu failing to find
    the object at all. A worthwhile refinement *if this is revisited*: allow the centre
    check a small (2–3px) tolerance rather than requiring the exact rounded click pixel.
    Implemented and tested below (7a) rather than left as a suggestion.
  - *Near-flat guard*: never hit in this population (0/452) — confirmed theoretical
    only, as anticipated.
  - *End-to-end smoke test*: `run_one_image(fn, ..., bbox_method="multiotsu")` runs
    to completion on 246.tiff and 300.tiff with no exceptions; `tightened_base_size` and
    `n_det` both come out finite and non-pathological (e.g. 246.tiff: 51px/23 dets under
    binary vs 33px/330 dets under multiotsu, on whichever seed each method's own
    candidate pool happened to draw — the large `n_det` swing is the template-size ->
    `TM_CCOEFF_NORMED`-null-variance effect from section 2's own caveat, not evidence of
    a bug).
- **Conclusion**: not adopted as the default. The motivating case doesn't clean up, and
  the collateral loss on already-good candidates (10.7% of large/solid components,
  concentrated in dense domains) is a real cost the aggregate exclusion-rate number
  alone would have hidden. Kept as an explicit opt-in (`method="multiotsu"` /
  `bbox_method="multiotsu"`) for anyone who wants to revisit it — e.g. after adding the
  centre-pixel tolerance noted above, which the near-miss numbers suggest could recover
  a majority of the 37 net-new losses without touching the fixed cases.

### 7a. Centre-pixel tolerance follow-up — more favourable, still not adopted

Implemented `tighten_box_otsu(center_tolerance=...)`: instead of requiring the exact
rounded click pixel to be foreground, accept the foreground pixel *nearest* the click
within an L-inf square of this half-width (`_nearest_label_within`), then proceed with
that pixel's component exactly as before. Default `0` (exact-pixel-only, unchanged);
threaded through `foreground_filter`/`pick_seed`/`tightened_base_size` and
`experiment.run_one_image(bbox_center_tolerance=...)` the same way `method` was.

Re-ran this section's full calibration (452 candidates, 10 ROIs) at tolerance
0/1/2/3px, binary unchanged throughout:

- **The collateral-loss objection largely dissolves.** The large-solid spot-check
  (section 7's strongest argument against adopting multiotsu) drops from 20/187 (10.7%)
  lost at tolerance 0 to 7/187 (3.7%) at tol=1, 3/187 (1.6%) at tol=2, 2/187 (1.1%) at
  tol=3.
- **Overall pool retention reaches parity with binary, then exceeds it.** Same
  504-candidate agreement-pool base as the rest of this section (452 of which are
  binary's centre-foreground count; the full 504 is the right denominator here since
  tolerance changes exactly that centre-foreground step). End-to-end survivors: binary
  447/504 (88.7%); multiotsu at tol=0/1/2/3 = 413/504 (81.9%), 432/504 (85.7%),
  448/504 (88.9%), 456/504 (90.5%). Tol=2 reaches binary's rate; tol=3 exceeds it. Per
  domain, multiotsu at tol=3 matches or beats binary everywhere except 246.tiff (85.6%
  vs 89.7%) and ties it on 406.tiff (both 50%, n=4) — a reversal from tol=0, where it was
  worse or equal everywhere.
- **The tightening effect is unaffected by tolerance, as expected** (tolerance only
  changes *whether* a component is found, not *which* component, once found): paired
  median size still 37px -> 26px at every tolerance level, ~93.6-93.8% of paired
  candidates strictly smaller.
- **The motivating example is still not fixed.** 350.tiff ann 18161 is byte-for-byte
  identical across every tolerance (813px area, 0.679 solidity, 50px box) — tolerance
  only changes whether a *different* pixel is checked, and ann 18161 was already
  centre-foreground even at tolerance 0, so this result was never touched by the fix.
  This still stands as the section's clearest evidence that multiotsu solves a
  centre-pixel brittleness problem, not the bridging problem it was proposed for.
- **New concern this follow-up surfaced: recovered components aren't guaranteed to
  actually be the click's own object.** `_nearest_label_within` accepts *whichever*
  labelled pixel is nearest within the tolerance window — by construction, for every
  candidate tolerance newly recovers, the click's own exact pixel was background, so the
  accepted component's mask never contains the click by definition. Checked this
  directly for the 41 candidates tol=2 recovers that tol=0 missed: the accepted
  component's *bounding box* contains the click in 34/41 (83%) -- 7/41 (17%) accept a
  component whose bbox doesn't even reach the click, i.e. plausibly a neighbouring
  nucleus, not the annotated one. Even restricted to those 34, the component's centroid
  sits farther from the click than the component's own approximate radius
  (`sqrt(area/pi)`) in 14/41 (34%) of all 41 recovered candidates -- consistent with a
  partial, off-centre capture of an irregular figure rather than a clean single-object
  hit. `tightened_base_size` still centres the returned box on
  the click regardless, so for this subset the resulting template can be sized to a
  component the click sits at the *edge* of, not the middle of -- a real, if partial,
  version of the same "wrong-object" risk the check was designed to catch.
  Sanity-exclusion counts corroborate this: 2 (tol=0) -> 7 -> 8 -> 13 (tol=1/2/3), rising
  faster than binary's constant 5 -- more of what tolerance recovers is marginal enough
  to then fail `min_area`/`max_area_frac`/`min_solidity`.
- **If this is pursued further, prefer tol=2 over tol=3**: tol=2 captures nearly all the
  large-solid recovery (1.6% lost vs 1.1%) with markedly less sanity-exclusion inflation
  (8 vs 13) and less exposure to the off-click risk above (a smaller window is less
  likely to reach into a neighbouring structure).
- **Still entirely untested**: downstream matching quality. Nothing here measures
  whether the tighter multiotsu templates actually discriminate mitoses from ordinary
  nuclei better or worse than binary's -- only that the candidate pool and box sizes
  differ. The `n_det` swing seen in the end-to-end smoke test (246.tiff: 23 detections
  under binary vs 1,897 under multiotsu at tol=2/3, same rng draw, different seed and
  template size) shows the effect on the matcher's output is large; its direction
  (better or worse discrimination) is unknown. This is the same `augmentation_variants`/
  `fusion_variants`-style AUC comparison section 6's caveat calls for on the augmentation
  axis, not yet built for the bbox-method axis.
- **Conclusion**: materially more evidence than section 7's original test, but still not
  adopted as the default. Pool retention and the large-solid collateral loss both
  improve enough at tol=2 to flip from "net loss" to "roughly competitive, with a
  meaningfully tighter box" -- but the newly-found off-click risk (17% of recovered
  candidates, worse than that on the centroid-distance measure) and the completely
  untested downstream discrimination question are enough to keep this an explicit
  opt-in (`method="multiotsu", center_tolerance=2` / `bbox_center_tolerance=2`),
  not a default. A natural next step, not done here: reject a tolerance-recovered
  component outright when the click falls outside its bbox, and re-run this
  calibration to see whether that removes the off-click risk without giving back the
  recovery gains.

### 7b. Bounding-box identity guard on tolerance-recovered components — implemented

Added the guard 7a's conclusion named: `tighten_box_otsu` now rejects the accepted
component outright (returns `None`, same as any other "unusable as a seed" case) when
the click's own pixel falls outside that component's bounding box. This only bites the
tolerance-recovered path -- when the exact click pixel is itself foreground, it is
trivially inside its own component's bbox, so `center_tolerance=0` (every existing
caller's default) is completely unaffected. No new parameter; the guard is unconditional
whenever a component is found via `_nearest_label_within`.

Re-ran the same 504-candidate, 10-ROI calibration at tolerance 0/1/2/3 with the guard
active (binary is untouched throughout, as always):

- **Sanity-check exclusions stop inflating.** This is the clearest confirmation the guard
  is doing its job: 7a saw exclusions climb 2 (tol=0) -> 7 -> 8 -> 13 as tolerance grew,
  i.e. more of what tolerance recovered was marginal enough to then fail
  `min_area`/`max_area_frac`/`min_solidity`. With the guard, exclusions are flat: 2 -> 5
  -> 5 -> 5 across tol=0/1/2/3 -- matching binary's own count (5) exactly, and
  spot-checking confirms these are mostly *different* annotations from binary's 5
  (only one, 300.tiff ann 14616, overlaps) -- a coincidence in count, not the same cases,
  but a reassuring one: multiotsu-with-guard's residual sanity failures are no larger a
  population than binary's own.
- **Large-solid loss still drops sharply, just not quite as far as the unguarded
  version.** 20/187 (10.7%) at tol=0, unchanged (the guard never touches already-good
  candidates) -> 8/187 (4.3%) at tol=1 -> 4/187 (2.1%) at tol=2 and tol=3 (7a's unguarded
  numbers were 7/3/2 at the same tolerances -- a small number of what 7a counted as
  "recovered" were exactly the off-click false recoveries this guard now correctly
  excludes again).
- **Overall pool retention is essentially unchanged, landing at parity with binary rather
  than exceeding it.** Survivors: binary 447/504 (88.7%); multiotsu-with-guard at
  tol=0/1/2/3 = 413/504 (81.9%), 431/504 (85.5%), 444/504 (88.1%), 448/504 (88.9%). Tol=3
  now matches binary almost exactly (88.9% vs 88.7%) rather than exceeding it (7a's
  unguarded tol=3 hit 90.5%) -- consistent with some of that extra unguarded headroom
  having been the same false, neighbour-object recoveries the guard removes.
- **Recovery count drops a little, as expected.** Of the 89 candidates multiotsu misses
  at tol=0, the guard recovers 21/34/38 at tol=1/2/3 (7a's unguarded numbers: 24/41/54)
  -- the guard is turning away roughly a third of what tolerance alone would have
  accepted, concentrated exactly where 7a's identity check flagged risk: 246.tiff (dense,
  lymphosarcoma) stops gaining anything past tol=1 (83.5% at tol=1/2/3 alike, versus a
  climb to 85.6% unguarded) -- the extra "recoveries" tolerance found there past 1px were
  disproportionately off-click.
- **Everything else is unchanged, as expected**: the paired tightening effect (median
  37px -> 26px, ~93.6-93.7% smaller), the thin-pool domain (350.tiff's 4 candidates all
  still recover identically), and the motivating example (350.tiff ann 18161, still
  813px/0.679 solidity/50px at every tolerance) are all untouched by a guard that only
  ever turns an accept into a reject, never the reverse.
- **End-to-end smoke test**: `run_one_image(fn, ..., bbox_method="multiotsu",
  bbox_center_tolerance=2)` (and `=3`) still runs to completion on 246.tiff and 300.tiff
  with no exceptions; sizes and detection counts remain finite and non-pathological. On
  300.tiff specifically, tol=2 and tol=3 land on *different* seeds (ann 14581 vs 14580),
  different tightened sizes (33px vs 23px), and roughly a 2x swing in `n_det` (7,679 vs
  15,677) -- the aggregate parity between tol=2 and tol=3 below is a population-level
  statement, not a claim that either tolerance behaves identically on a given image.
- **Re-checked 7a's softer identity concern (click-to-centroid distance) on the
  post-guard recovered set, not just bbox-containment.** By construction every
  guard-recovered candidate now has `click_in_bbox = True` (34/34 at tol=2, confirmed
  directly against `tighten_box_otsu`'s own output rather than re-derived) -- the guard
  does fully close the failure it targets. But bbox-containment doesn't imply the click
  sits near the *middle* of what got captured: of those same 34, only 27 (79.4%) have
  the component's centroid within its own approximate radius (`sqrt(area/pi)`) of the
  click, versus 27/41 (65.9%) before the guard. So the guard measurably improves this too
  (fewer partial/off-centre captures survive, since guard-rejected components disproportionately
  had extreme centroid offsets) but does not eliminate it: roughly 1 in 5 guard-recovered
  candidates at tol=2 still have a click sitting outside the captured component's own
  radius -- e.g. an elongated neighbouring structure whose bbox reaches the click while
  its mass sits elsewhere would pass the bbox guard and still land here. Not addressed by
  this section; a stricter guard (reject on centroid distance, not just bbox) is the
  natural next check if this is pursued further.
- **Conclusion**: the guard is a clean, low-cost fix for the bbox-containment failure 7a
  surfaced -- that specific risk (17% of recovered candidates, tol=2) goes to 0% by
  construction, at the cost of giving back only a small fraction of tolerance's recovery
  gains (which turn out to have been disproportionately the same false recoveries), and
  it explains 7a's rising sanity-exclusion count as a symptom of the same underlying
  issue rather than a separate one. It does not fully resolve the softer, centroid-distance
  version of the same concern (still ~1 in 5 at tol=2, down from ~1 in 3). **Still not
  adopted as the default**, for the same two reasons that survive this fix unchanged: the
  section-7 motivating example (350.tiff ann 18161) is still not actually improved by any
  of this, and downstream matching quality (whether the tighter multiotsu templates
  discriminate mitoses from ordinary nuclei better or worse than binary's) has not been
  measured at any point across 7, 7a, or 7b -- only the seed-selection candidate pool and
  box geometry have been. Recommended setting if pursued further is still
  `center_tolerance=2` over `3`: population-level survival is close either way (444 vs
  448 of 504), but the end-to-end smoke test shows tol=3 can pick a different seed and a
  markedly different template size/detection count on the same image (300.tiff above),
  and a smaller window is inherently less likely to reach a neighbouring structure in the
  first place -- fewer chances for the still-open centroid-distance risk to bite, not a
  claim that tol=3 is measurably worse in aggregate.

## 8. No-bbox-tightening variant: agreement-only seed selection, native box size

Added `tighten_bbox: bool = True` to `seed_selection.pick_seed` and
`experiment.run_one_image`/`run_experiment` (threaded through identically to
`bbox_method`/`bbox_center_tolerance`). Default `True` is the unchanged pipeline
behaviour from sections 1-2/7/7a/7b; `False` is the variant this section adds.

- **What it does**: seed selection becomes pathologist agreement + the border filter
  only -- `foreground_filter` (the Otsu/CC step, in either `method`) is skipped
  entirely, so an annotation is never excluded for having its click land outside a
  foreground component. The template also keeps `cfg.base_size`'s native size
  (`FSConfig`'s own default, `template_match.BASE_SIZE` = 51px, or whatever the caller's
  `cfg` already specifies) instead of being resized to a tightened box --
  `tightened_base_size` is never called. `bbox_method`/`bbox_center_tolerance` become
  irrelevant in this mode (there is no Otsu/CC step left for them to configure).
- **Why this wasn't already a flag**: sections 1-2 built agreement and bbox-tightening
  as if they were both mandatory parts of one seed-selection pipeline; this section
  makes the second one optional, so section 1's original question ("does agreement
  alone, without any Otsu/CC filtering, still work as a seed-selection strategy") is
  actually answerable rather than only inferable from reading `foreground_filter`'s
  source.
- **Verification**:
  - `pick_seed(..., tighten_bbox=False)` directly, on 002.tiff/350.tiff/246.tiff:
    `n_after_foreground` correctly equals `n_after_border` in every case (no filtering
    ran). On 002.tiff and 350.tiff the selected seed is unchanged from `tighten_bbox=True`
    (both domains happened to have zero foreground-filter exclusions already, so the
    pool is identical either way); on 246.tiff (a dense domain where the foreground
    filter does exclude 10/97 candidates under `tighten_bbox=True`) the selected seed
    correctly differs, since the candidate pool itself is larger and differently
    ordered.
  - `run_one_image(..., tighten_bbox=False)` end-to-end on 002.tiff/350.tiff/246.tiff/
    300.tiff: `tightened_base_size` is 51px in every case (the untouched
    `FSConfig.base_size` default); passing a deliberately non-default `cfg.base_size=41`
    confirms it is respected, not silently overwritten. `n_det` stays finite and
    non-pathological throughout (7 to 3,302 across the four domains -- no crashes, no
    degenerate zero- or runaway-detection counts).
  - `run_experiment(..., tighten_bbox=False)` on a 2-domain batch: the flag (and
    `bbox_method`) correctly reach every row's metrics (`tighten_bbox=False`,
    `tightened_base_size=51` for both domains in the batch); `tighten_bbox=True` through
    the same entry point is unaffected.
  - Empty-pool edge case: 001.tiff (zero mitotic annotations at all) still raises
    `pick_seed`'s "no seed candidates left" error under `tighten_bbox=False`, correctly
    naming the "agreement" stage (not "foreground", which never ran) as the one that
    emptied the pool.
  - Full regression suite: `verify_fixes.py` (100% pass), the `plant_and_recover`
    coordinate gate, and a full re-execution of `seed_selection_midog.ipynb` (0 errors,
    unaffected since it never sets `tighten_bbox` and the new default is unchanged) all
    still pass after these changes.
- **Not yet done**: this section only confirms the variant runs correctly end-to-end and
  is wired consistently through the pipeline -- it does not compare `tighten_bbox=False`
  against `True` on discrimination or FP rate. Given section 7b's conclusion (tighter
  templates plausibly *raise* FP rate under `TM_CCOEFF_NORMED`'s template-size-dependent
  null variance, not lower it), the native-size variant this section adds is the more
  conservative default to actually run experiments with next, but that is a prediction
  from the documented mechanism, not a result measured here.

## 9. Score threshold 0.75 — tested, NOT adopted (0.5 retained)

Ran the section-6 simplified configuration again with `FSConfig.score_threshold` raised
from 0.5 to **0.75** and nothing else changed — same seven ROIs, same `seed=0`, same
binary-Otsu tightening, same three channels. Notebook:
`find_and_suppress_midog_simple_t075.ipynb`; results in `results/fs_t075_*.csv`.

- **The outcome is forced by the pipeline's structure, not discovered.** Three facts:
  `extract_peaks` applies the local-maximum test (grey dilation) *independently* of the
  score floor and only then intersects with `fused >= score_threshold`
  (`template_match.py:188`); `max_peaks` (250 000) never binds on these ROIs; and
  `nms_by_distance` is greedy in **descending** score order, so a peak scoring ≥ 0.75 can
  only be suppressed by a higher-scoring peak, which is also ≥ 0.75 and therefore also
  present in the 0.75 candidate set. Therefore
  `detections@0.75 = {d ∈ detections@0.5 : d.score ≥ 0.75}`, same ranks, same order.
  Consequently `recall@K` can only stay equal (when the surviving list still holds ≥ K
  detections, in which case the top-K is byte-identical) or fall (when it does not); it
  can never rise. The FROC curve is unchanged — the 0.75 curve is an initial segment of
  the 0.5 curve. **The threshold is an operating point on a fixed curve, not a change to
  the curve.**
- **Verified, not assumed.** The expected 0.75 result was computed from the *saved* 0.5
  detection file before the re-run, then compared against it: zero coordinate-set
  differences on all seven images, and observed `recall@K` equal to the predicted value
  on every row. No tie-ordering divergence appeared despite `np.argsort`'s default
  quicksort being unstable.
- **Measured effect.** 0/21 (image × channel) rows improved; 14 unchanged; 7 got worse,
  all by the list emptying completely. Mean `recall@K` fell for every channel
  (rgb 0.111 → 0.079, hematoxylin 0.078 → 0.036, gray_inverted 0.067 → 0.056). RGB total
  detections 25 592 → 978; mean `coverage_frac` 0.263 → 0.012.
- **Why three ROIs zero out.** 201.tiff / 246.tiff / 506.tiff have best *genuine*
  (non-self) match scores of 0.64 / 0.60 / 0.74 — all below 0.75, so the floor sits above
  everything the matcher found there and the whole detection list, FROC curve included,
  is empty. Under RGB that is 3/7 ROIs; the same three empty out under all three channels
  (506.tiff survives with 3 and 10 detections under gray_inverted/hematoxylin).
- **Precision and coverage are not evidence here, and precision does not even move in
  one direction.** `coverage_frac` falls by an order of magnitude everywhere, mechanically,
  because the list is shorter. `all_precision_mitotic` rises only on the two rows where a
  true positive survives the truncation (301.tiff 0.015 → 0.124; 002.tiff 0.002 → 0.048);
  on the other five it falls to 0 (genuinely 0/51 and 0/4 on 405.tiff and 350.tiff, and an
  undefined 0/0 on the three emptied ROIs). Truncating a signal-bearing ranked list higher
  raises precision only while the surviving prefix still contains a hit — cut above every
  hit and the same operation drives it to zero, which is the emptying collapse seen through
  a second metric. The one real gain is interpretability: at 0.5, `coverage_frac` was
  0.29–0.84 on four ROIs, high enough that full-list recall was partly a statement about
  tiling (see `evaluate.py`'s module docstring); at 0.75 it is ≤ 0.08 everywhere — but only
  on the rows where anything survives.
- **The threshold-free probe is unchanged by construction.** `experiment._probe_variants`
  reads the fused response directly and never touches `score_threshold`, so
  `auc_mitosis_vs_nucleus` and `discrimination` are provably identical at both floors —
  the cleanest available demonstration that 0.75 changed *selectivity*, not *signal*.
  `results/fs_simple_channel_probe.csv` is reused rather than re-run (~5 min of
  byte-identical output avoided).
- **Decision: keep 0.5.** A fixed absolute floor is the wrong instrument for this
  pipeline. `TM_CCOEFF_NORMED` scores are not comparable across ROIs — the tightened
  template size varies 35–51 px between images and the statistic's null variance depends
  on template size (section 7b) — so one number does not mean the same thing on 301.tiff
  as on 246.tiff. Shorter, more precise lists are better obtained by truncating at report
  time (by rank, or by score), which yields the identical detections while preserving the
  option to read further down on ROIs whose whole score distribution sits low. If a floor
  is wanted at all, a per-image quantile of that ROI's own peak distribution would be a
  consistent operating point; a constant is not.

## 10. Reference NMS ordering (keep the right-most box) — tested, NOT adopted

`bbox tuning code reference/bbox_tuning.py:483` (`nms_with_area`) resolves an overlapping
cluster by keeping whichever box sits **furthest right**: the frame is sorted by
`['humanMade', 'nms applied', 'x top left', 'y top left']` (bbox_tuning.py:785), the loop
takes `last` — the largest `x top left` — and suppresses the smaller-x members. The match
score is never consulted, and the column is dropped entirely at bbox_tuning.py:798.
`midog_utils/nms.py` instead walks in descending score order. Measured what that choice is
worth. Script: `nms_ordering_probe.py`; results in `results/fs_nms_ordering_metrics.csv`
and `results/fs_nms_ordering_tp_flips.csv`.

- **Method.** Same fused correlation map, same extracted peaks, same suppression radius
  (each image's evaluation match radius, 29.6–33.1 px), same self-hit removal; *only* the
  order the greedy loop walks in changes. The suppression criterion stays the distance
  test in both arms — changing the geometry too would confound ordering with criterion.
  Survivors are ranked by score in both arms so every metric stays well defined; this
  charges the reference rule only for picking the wrong cluster representative. RGB,
  `score_threshold=0.5`, `seed=0`, all seven ROIs.
- **Result: no measurable difference.** Mean `recall@K` 0.1114 (score-ordered) vs 0.1141
  (x-ordered). Six of seven images are identical on every metric and show **zero** TP
  flips in either direction. All movement is on 301.tiff, the densest ROI (218 mitotic in
  2 mm²): 0.304 → 0.323, from 6 flips (1 TP lost, 5 gained) out of 13,415 detections.
  Whether that is noise or a real geometric effect cannot be told from one dense image.
  There *is* a candidate mechanism: the minimum spacing between two MIDOG++ annotations is
  26.2 px (403.tiff; 245.tiff is 26.6), *below* the ~30 px NMS-and-match radius
  (`find_and_suppress.py`'s module docstring), so where two GT sit that close one detection
  can cover only one of them and which of the pair is covered depends on the tiling the
  suppression produces. Score-NMS tiles greedily outward from the global maximum; x-NMS
  tiles as a right-to-left sweep. Nothing says the former covers closely-spaced GT better.
  That the only ROI to move is the densest one is consistent with this — and with noise.
  Distinguishing them needs more seeds and more dense ROIs; **the decision below does not
  rest on this null.**
- **Why the effect is so small — NMS is nearly inert here.** Peaks before vs after
  suppression: 16236→13353, 129→127, 6→5, 4669→4479, 4030→3621, 2159→2067, 2009→1940.
  Only 3–18% of peaks are ever in a contested cluster, because `peak_min_distance=7`'s
  grey dilation has already thinned the map and the surviving local maxima above the 0.5
  floor are mostly further apart than the ~30 px NMS radius. An ordering rule can only
  matter where clusters exist, and here they barely do. Under the x rule 0.8–12% of
  survivors are not the highest-scoring peak within their own radius (vs 0–1.5% under the
  score rule), displaced by a median 11–27 px — real, but too rare to move a metric.
- **The cost is structural, not statistical, and is not what the above measures.** The
  reference does not merely choose survivors by position; it emits them **in x order with
  no score** (bbox_tuning.py:798). `recall@K`, FROC and `sens@Xfp` all consume a ranked
  list, so applying them to an x-ordered one is not a worse number but a different
  question — the entire evaluation of this pipeline would become undefined. Greedy NMS's
  correctness argument also depends on descending-score order: it is exactly the property
  that makes section 9's threshold result exact (a peak ≥ t can only be suppressed by a
  higher-scoring peak). Order by x and a low-scoring peak can delete a cluster's best
  peak purely for sitting further right.
- **Why the reference gets away with it.** Its consumer is a particle tracker that wants
  *a* box per object per frame, not a ranked candidate list; positional order is
  harmless, arguably convenient, there. This pipeline's contract is different.
- **Decision: keep score-ordered NMS.** Not because the positional rule was measured to
  hurt — it was not — but because it buys nothing measurable while removing the ranking
  every metric here is built on. Scope of the null: RGB only, one threshold, one seed per
  image; it says the ordering does not matter *at this operating point*, where NMS itself
  barely binds — not that it could never matter in a denser regime, and the 301.tiff row
  is the reason to say so explicitly.

## Known limitation none of this addresses

Mitotic figures are 0.09–1.1% of nuclei in a 2 mm² ROI; at an AUC of 0.94 (002.tiff,
`score_probe`), roughly 500 ordinary nuclei still outrank a typical mitosis on score
alone. That's a base-rate problem, not a template-quality problem — none of the three
changes above change the ratio. Expect incremental movement in `recall_at_k`, not a
transformation of it.
