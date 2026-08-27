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

## Known limitation none of this addresses

Mitotic figures are 0.09–1.1% of nuclei in a 2 mm² ROI; at an AUC of 0.94 (002.tiff,
`score_probe`), roughly 500 ordinary nuclei still outrank a typical mitosis on score
alone. That's a base-rate problem, not a template-quality problem — none of the three
changes above change the ratio. Expect incremental movement in `recall_at_k`, not a
transformation of it.
