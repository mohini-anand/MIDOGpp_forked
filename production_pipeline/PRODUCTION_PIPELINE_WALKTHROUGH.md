# The production pipeline, step by step

Written 2026-09-13. Scope: `midog_utils/production.py` and everything it calls — the
callable that turns one pathologist click into a ranked list of detections on that ROI.
This is the **current, working implementation**, not a notebook convention: it replaces
the by-hand reimplementation every `production_seed_precision_at_k*` notebook used to
carry (see `production.py`'s own module docstring).

**How this doc was verified, not just written.** Every function and line number below was
checked against the source at the time of writing (`git status` confirmed no uncommitted
edits landed between reading and writing). Beyond reading the code, the whole pipeline was
**executed**, twice, independently:

1. `production_pipeline/production_pipeline_demo.ipynb` was run end-to-end with
   `/Users/mohinianand/anaconda3/bin/jupyter nbconvert --execute --inplace` on both its
   ROIs (`245.tiff`, `403.tiff`) and both ranking axes. Zero errors. The notebook already
   contained saved outputs from an earlier run, and this run reproduced them
   **byte-for-byte** — the pipeline has no hidden nondeterminism.
2. A standalone script called `run_production_pipeline` directly (bypassing the notebook)
   and printed the full `info` dict the function returns — every field, not just the
   subset the demo notebook prints — to confirm each internal stage actually did what its
   docstring claims. Results are quoted inline below, at the stage they verify.

If you edit `production.py` or anything it imports, re-run the demo notebook and diff its
outputs — that reproducibility is itself a regression test.

**Updated 2026-09-17:** the seed's own match is now excluded by blanking its refined-template footprint (`base_size x base_size`, centred on `template_xy`) out of a copy of the search channel before correlation runs (`template_match.blank_seed_square`), replacing a 5.0 px post-NMS filter. See `SELF_HIT_MASKING_PLAN.md` and `DECISIONS.md` D11. Numbers below were re-verified after that change.

**Re-synced 2026-09-16** with the dead-code cleanup (`c066829`, `PRODUCTION_PIPELINE_CLEANUP.md`),
which left production outputs bit-identical. Line numbers below are from the 2026-09-13 code,
archived as `midog_utils_full/`.

---

## Relationship to the rest of the repo

- **`PIPELINE_BREAKDOWN.md`** (repo root, 2026-09-12) is the historical record of six
  gaps between "what `DECISIONS.md` says should happen" and "what the notebooks actually
  ran" that existed before `production.py` was written. This document does not replace
  it — `production.py` *is* that gap-closing work, and `PIPELINE_BREAKDOWN.md` is kept as
  the record of why each constant below is set the way it is (measurement provenance,
  caveats). All six of its punch-list items are resolved in the code described here; see
  "Punch-list resolution" near the end.
- **`midog_utils/compare.py`**'s `Arm`/`evaluate_arms` machinery, which every earlier
  notebook used for ranking, and `chromatin.rerank` were removed in `c066829`.
  `production.py` ranks with its own inline stable sort.
- **`midog_utils/FIND_AND_SUPPRESS_REFERENCE_DIFFS.md`** documents how `find_and_suppress`
  deliberately diverges from the original reference implementation it replaced (no image
  masking, scores kept for ranking, nothing random, etc.) — read it if you're asking "why
  doesn't this match the reference code."

---

## Verification summary

Driving `run_production_pipeline` directly and printing its full `info` dict on both test
ROIs gave:

| field | 245.tiff / tm_score | 403.tiff / tm_score | what it confirms |
|---|---|---|---|
| `n_peaks` | 100 | 100 | D9's `max_peaks=100` cap binds pre-NMS, as designed |
| `max_peaks_binding` | True | True | recorded, not invariant-checked (see below) |
| `n_blanked_px` | 2209 | 2601 | pixels of the refined template's own footprint (47^2 / 51^2) excluded before correlation |
| `n_after_nms` | 90 | 99 | NMS removed 10 / 1 of the 100 pre-NMS peaks |
| `n_detections` | 90 | 99 | final ranked-list size; nothing is removed after NMS |
| `deep_floor_median` / `mad` | -0.0074 / 0.2143 | -0.0071 / 0.7181 | the adaptive robust-z floor actually computed per-ROI, not a leftover fixed threshold |
| `score_threshold_used` | -0.329 | -1.084 | deliberately permissive — "deep floor" is not a strict filter |

No detection can be extracted from inside the blanked square, because its source pixels are
gone before correlation runs. The filter this replaced removed only points within 5.0 px, so
it missed the seed's own match whenever that match fell below the `max_peaks` cutoff or moved
under an augmented template bank (`FIND_AND_SUPPRESS_REFERENCE_DIFFS.md`).

Visual check: `viz.overlay` renders correctly on both ROIs — ground-truth boxes, a
white/black star at the click, a cyan `+` at the D8-recentred template centre (visually
indistinguishable from the click at this zoom, consistent with the printed 2-4 px
offsets), and detection circles colour-coded by evaluation bucket.

---

## Stage A — Load the data (`midog_utils/dataset.py`)

| Function | Location | What it does |
|---|---|---|
| `load_annotations(json_path)` | [dataset.py:55](../midog_utils/dataset.py#L55) | Reads `MIDOG++.json` into two tidy tables: `images` (one row per ROI) and `annotations` (one row per pathologist click, with `cx`/`cy` = box centre, `n_votes`/`n_mitotic_votes` = how many raters saw it and how many called it mitotic). |
| `load_roi(path)` | [dataset.py:155](../midog_utils/dataset.py#L155) | Reads one ROI TIFF into a full-resolution RGB `uint8` array. Handles both this dataset's TIFF layouts (flat and 7-level LZW pyramid) and strips the (uniform) alpha channel. |
| `roi_mpp(path)` | [dataset.py:171](../midog_utils/dataset.py#L171) | Reads the TIFF's resolution tags to get microns-per-pixel. Must honour `ResolutionUnit` — two of the four scanners in this dataset store centimetres, not inches; assuming inches would be silently wrong for those two. |
| `image_annotations(annotations, file_name)` | [dataset.py:141](../midog_utils/dataset.py#L141) | Filters the full annotations table down to one ROI's rows. |

The caller then does `gt_mitotic = gt[gt.category_id == ds.MITOTIC]` inline — a plain
pandas filter, not a function call — to get the pool `build_seed` draws from.

---

## Stage B — Pick and gate the seed (`midog_utils/seed_selection.py`)

Entry point: **`build_seed(gt_mitotic, structural_channel, rng, roi_shape)`**,
[seed_selection.py:513](../midog_utils/seed_selection.py#L513). One call replaces the
agreement/border/gate/draw/retry block every earlier experiment wrote inline.

```
build_seed(gt_mitotic, gray_inv, rng, roi_shape)              seed_selection.py:513
│
├─ agreement_pool(gt_mitotic)                                 seed_selection.py:48
│    Prefer annotations every rater called mitotic (unanimous). If an image has none,
│    fall back to the 2-of-3 contested tier and flag it — a seed the experts themselves
│    argued about risks building a poor template from an atypical example.
│
├─ border_filter(pool, border, roi_shape)                     seed_selection.py:69
│    Drop annotations too close to the ROI edge to ever read a full rotation-safe patch
│    (border = patch_size // 2 = 36 px).
│
└─ LOOP — draw without replacement, retry on rejection:
     │
     ├─ (draw one row: cx, cy = the raw click)
     │
     └─ tightened_template_box(structural_channel, cx, cy, ...)   seed_selection.py:431
          "the D8-current seed constructor: click gates, component recentres"
          │
          ├─ tm.read_padded_patch(channel, cx, cy, otsu_window=51)  template_match.py:102
          │    Cuts a 51x51 px window centred on the raw click. Returns None if the
          │    click is too close to the ROI edge to read the full window.
          │
          └─ tighten_box_otsu(patch, min_area=50, max_area_frac=0.85, min_solidity=0.5)
                                                                   seed_selection.py:104
               "the gate" — is the click actually sitting on real chromatin?
               ├─ cv2.threshold(u8, 0, 255, THRESH_BINARY + THRESH_OTSU)
               │    two-class Otsu split of the 51x51 window into foreground/background
               ├─ skimage.measure.label(binary, connectivity=2)
               │    connected-component labelling of the foreground
               ├─ center_label = labels[cy, cx]
               │    the component under the CLICK'S OWN rounded pixel — not "the
               │    largest component in the window", a materially different rule
               ├─ reject (return None) unless ALL of:
               │    - center_label != 0 (the click's pixel is foreground at all)
               │    - the click falls inside that component's own bounding box
               │    - region.area >= 50 px            (not Otsu noise)
               │    - region.area <= 0.85 * 51*51 px   (not an implausibly large blob)
               │    - region.solidity >= 0.5           (not two nuclei bridged together)
               └─ returns (y0, y1, x0, x1) — the accepted component's half-open bbox,
                    in patch-local coordinates
          │
          back in tightened_template_box:
          base_size = odd(max(y1-y0, x1-x0))          -- the template's native size
          half = otsu_window // 2                     -- 25, for the 51 px window
          ix, iy = round(cx), round(cy)                -- matches read_padded_patch's own rounding
          center_x = ix - half + (x0 + x1 - 1) / 2.0   <- D8's half-pixel correction
          center_y = iy - half + (y0 + y1 - 1) / 2.0      (a half-open bbox covers pixels
                                                            x0..x1-1, so (x1-1) is the last
                                                            real pixel — using x1 directly
                                                            would be half a pixel too far)
          returns (base_size, center_x, center_y)
     │
     └─ _patch_readable(roi_shape, center_x, center_y, patch_size=73)  seed_selection.py:501
          Re-checks that the RECENTRED point (which can sit closer to the edge than the
          original click did) still has room for a full 73x73 rotation-safe read. If not,
          this candidate is refused and the loop draws again — same as any other rejection.
│
└─ returns Seed(ann_id, click_xy, template_xy, base_size, recentred, offset_px, ...)
                                                               seed_selection.py:478 (dataclass)
     click_xy    = the raw pathologist click — the permanent ground-truth reference
     template_xy = the gated + recentred point — where the search template is actually cut
     These are DIFFERENT points (the template is always recentred). Nothing downstream may
     substitute one for the other: the search-channel blanking and NMS-adjacent radii
     reference template_xy; ground-truth exclusion and every match-radius computation reference
     click_xy.
```

After `build_seed` returns, the caller (not `build_seed` itself) excludes the seed's own
annotation from the evaluation set:

```python
gt_eval = gt[gt.ann_id != seed.ann_id].reset_index(drop=True)
```

**Why `tighten_box_otsu`'s gate is not "is the click on the biggest thing nearby":** it
reads the label under the click's own pixel (`labels[cy, cx]`), while a naive
largest-component rule would take `max(regions, key=area)`. The two coincide only when the
click's own component also happens to be the biggest thing in the 51 px window — true on
`245.tiff`'s test seed, false in general, and the reason `403.tiff` needed 1 retry in the
verification run (its first draw's click-component didn't pass the gate).

---

## Stage C — Run the search and rank detections (`midog_utils/production.py`)

Entry point: **`run_production_pipeline(rgb, seed, mpp, rank_key='tm_score', max_peaks=100)`**,
[production.py:44](../midog_utils/production.py#L44).

Module-level constants wiring in the decision record ([production.py:30-41](../midog_utils/production.py#L30)):

| Constant | Value | Decision |
|---|---|---|
| `CHANNEL` | `"hematoxylin_od"` | D3 — unclipped optical density |
| `TM_METHOD` | `cv2.TM_CCOEFF` | D1 |
| `PEAK_MIN_DISTANCE` | 7 | local-maxima thinning window |
| `DEEP_FLOOR_Z` | -1.5 | adaptive floor = median - 1.5 * MAD, deliberately permissive |
| `MAX_PEAKS` | 100 | D9 |
| `OD_WINDOW` | `tm.BASE_SIZE` = 51 | `chromatin_density`'s window for the opt-in axis |
| `AXES` | `{"tm_score": "score", "chromatin_od": "od"}` | the one shared rank-key registry |

```
run_production_pipeline(rgb, seed, mpp, rank_key)              production.py:44
│  "one gated seed in, ranked detections on the whole ROI out"
│
├─ ch.to_channel(rgb, "hematoxylin_od")                        channels.py:78
│    └─ to_hematoxylin_od(rgb)                                 channels.py:41
│         └─ cm.hematoxylin_od(rgb)                            chromatin.py:102
│              skimage.color.rgb2hed(rgb.astype(f32)/255)[:,:,0]
│              UNCLIPPED hematoxylin optical density (D3) — dark chromatin = high value,
│              never rescaled to 0-255 (rescaling would saturate exactly the densest,
│              most informative pixels at the ceiling)
│
├─ ev.radius_px(mpp)                                           evaluate.py:53
│    7.5 micrometres (MIDOG_RADIUS_UM, evaluate.py:50) / mpp -> pixels.
│    Computed PER IMAGE from that scanner's resolution (29.6-33.1 px across this
│    dataset's four scanners) — never a hardcoded pixel constant (D7).
│
├─ FSConfig(base_size=seed.base_size, scales=(1.0,), n_angles=1, flips=(False,),
│           peak_min_distance=7, max_peaks=100, nms_radius=...,
│           deep_floor_z=-1.5, border_pad=True,
│           tm_method=cv2.TM_CCOEFF)                     find_and_suppress.py:34 (dataclass)
│    Just a settings bundle — nothing computed yet.
│
├─ find_and_suppress(hem, seed.template_xy, cfg, nms_radius)   find_and_suppress.py:101
│  │  "cut a template, blank its own footprint from a copy, correlate that copy, suppress"
│  │
│  ├─ tm.read_padded_patch(hem, tx, ty, patch_size=73)         template_match.py:102
│  │    73x73 crop centred on the TEMPLATE point (seed.template_xy, not the click).
│  │    73 = ceil(51*sqrt(2)): large enough that the inner 51x51 template survives an
│  │    arbitrary rotation without needing fabricated (padded) pixels.
│  │
│  ├─ tm.build_augmentations(patch, base_size, scales=(1.0,), n_angles=1, flips=(False,))
│  │                                                            template_match.py:118
│  │    With production's settings this yields exactly ONE template — no rotation, no
│  │    flip, no scale variant. The augmentation-bank machinery exists for the general
│  │    case (up to 12 angles x 2 flips x 3 scales elsewhere in the repo) but is a no-op
│  │    at these settings today.
│  │
│  ├─ tm.blank_seed_square(hem, tx, ty, base_size)              template_match.py
│  │    Blanks a base_size x base_size copy of hem at the template's own footprint --
│  │    never the original hem, which production.py reuses unblanked for chromatin_od
│  │    ranking. Verified: n_blanked_px 2209 on 245.tiff, 2601 on 403.tiff.
│  │
│  ├─ [cfg.border_pad=True] cv2.copyMakeBorder(blanked_hem, pad, ..., BORDER_REPLICATE)
│  │    Pads the whole ROI by half the template size so correlation reaches every real
│  │    ROI pixel, including near the border, without an unreachable band. Note this
│  │    pads the BLANKED copy from the node above — every search stage below runs on
│  │    that copy. The unblanked hem survives untouched for chromatin_od ranking.
│  │
│  ├─ tm.fused_response(padded_blanked_hem, templates, scale_normalize=False, method=TM_CCOEFF)
│  │                                                            template_match.py:184
│  │    For each template (here: one), cv2.matchTemplate(image, template, TM_CCOEFF).
│  │    The raw map is indexed by template TOP-LEFT; it's shifted by (size-1)//2 so the
│  │    fused map is indexed by CENTRE instead — getting this wrong shifts every
│  │    detection by 10-25 px, which at a 30 px match radius looks exactly like "template
│  │    matching doesn't work." With a single template the "fuse across augmentations"
│  │    step (element-wise max) is structurally a no-op but runs the same code path.
│  │    Returns (fused score map, winning-augmentation-index-per-pixel, valid-pixel mask).
│  │
│  ├─ crop the padded result back down to the ROI's original shape
│  │
│  ├─ tm.robust_stats(fused, valid)                            template_match.py:66
│  │    Median and 1.4826*MAD of the score map on a stride-8 sample of valid pixels —
│  │    this ROI's own "typical background correlation," recomputed per run (verified:
│  │    245.tiff and 403.tiff got different median/MAD, confirming this isn't a stale
│  │    fixed value).
│  │
│  ├─ threshold = median + (-1.5) * MAD                        (cfg.deep_floor_z)
│  │    Required: find_and_suppress raises ValueError if cfg.deep_floor_z is None.
│  │    Deliberately permissive — the label "deep floor" means almost everything above
│  │    noise survives to the next stage; this is not the pipeline's real selectivity.
│  │
│  ├─ tm.extract_peaks(fused, valid, min_distance=7, threshold, max_peaks=100)
│  │                                                            template_match.py:250
│  │    cv2.dilate for cheap local-maxima thinning (k = 2*7+1 = 15 px window), keep
│  │    points >= threshold, sort by a GLOBAL key (score descending, then x, then y —
│  │    not just score, so exact score ties break the same way regardless of pool size),
│  │    keep the top 100. Verified: n_peaks was exactly 100 on both test ROIs — the
│  │    ~17,000+ candidate peaks a 39-megapixel ROI produces are cut down to 100 every
│  │    single time (D9's cap binds as designed).
│  │
│  ├─ nms_by_distance(centers, scores, nms_radius)             nms.py:22
│  │    Greedy suppression in descending-score order: a KD-tree radius query drops every
│  │    remaining candidate within ~30 px of each kept point. Verified: this stage alone
│  │    removed 10 of 100 candidates on 245.tiff and 1 of 100 on 403.tiff.
│  │
│  └─ returns (detections DataFrame[rank, cx, cy, score, angle, flip, scale], info dict)
│       90 rows on 245.tiff, 99 on 403.tiff (100 - NMS losses)
│
├─ info["max_peaks_binding"] = (n_peaks == max_peaks)          production.py:80
│    Recorded as a fact, not invariant-checked.
│
├─ inv.check_nms_radius(nms_radius, mpp, label=...)            invariants.py:177
│    Asserts the radius actually used equals evaluate.radius_px(mpp) recomputed fresh.
│
├─ [only if rank_key == "chromatin_od"]:
│    │
│    ├─ cv2.copyMakeBorder(hem, pad=25, ..., BORDER_REPLICATE)
│    │    Pads so a detection near the ROI border still gets a fully-readable window.
│    │
│    └─ cm.score_detections(shifted_detections, hem_padded, window=51)  chromatin.py:140
│         └─ cm.chromatin_density(hem_padded, cx, cy, window=51, frac=0.10)
│                                                                  chromatin.py:113
│              For each detection: cut a 51x51 window around it, take the mean of the
│              darkest (highest-value, in this "more object = higher value" convention)
│              10% of pixels in that window. This is "how dense is the chromatin
│              actually here" — an absolute-darkness statistic, independent of the
│              template-match correlation score.
│    detections = detections.assign(od = the computed values)
│
├─ info["od_computed"] = (rank_key == "chromatin_od")
│
├─ detections.sort_values(AXES[rank_key], ascending=False, kind="mergesort")
│    STABLE sort. Because `pool` arrives already tm_score-descending (from
│    find_and_suppress), any tie on the chromatin_od axis breaks in tm_score order —
│    tm_score is the implicit secondary sort key for every axis, chromatin_od included.
│
├─ detections = detections.assign(rank = 0..n-1)               renumber after the re-sort
├─ info["rank_key"] = rank_key
│
└─ returns (detections, info)
```

**What `find_and_suppress` deliberately does NOT do**, per
`midog_utils/FIND_AND_SUPPRESS_REFERENCE_DIFFS.md`: no masking of every existing annotation
before correlating (only the seed's own template footprint is blanked — narrower and
well-defined: TM_CCOEFF against a flat window computes to ~0, not an extreme value, so no
zero-variance-window artefact); no boundary refinement, stability filter,
or aspect-ratio/PSNR filters (evaluation here is point-based, so box geometry moves no
metric, and the stability check passes 0/9 objects on H&E anyway); nothing random
anywhere in the search (the reference implementation used unseeded `sample`/`random.sample`
calls this port replaced with caller-controlled seeding and score-based truncation).

---

## Stage D — Scoring and visualization (caller-owned, not inside `production.py`)

`run_production_pipeline` returns a ranked list and stops. Turning that into a
precision/recall number or a picture is the caller's job — shown here because the demo
notebook does both immediately afterward.

```
ev.bucket_detections(det, gt_eval, match_radius)               evaluate.py:108
│    gt_eval must already have the seed's own annotation removed and pool both
│    categories (mitotic + look-alike), so one detection can't claim a mitotic figure and
│    a look-alike simultaneously.
│
└─ ev.greedy_match(det_xy, gt_xy, radius)                       evaluate.py:57
     KD-tree radius query, processed BEST-FIRST in the ranked detections' own order:
     each detection claims the nearest still-unclaimed ground-truth point within
     `radius`; one-to-one (a claimed GT point can't be claimed twice). Order is
     deliberate — a top-K prefix's matches never depend on lower-ranked detections,
     per the module's own docstring.

     labels each detection: HUMAN_CORRECT_LABEL (claimed a mitotic figure),
     HUMAN_REJECTED_LABEL (claimed a look-alike), or NON_HUMAN_FINDINGS (claimed nothing)

viz.overlay(rgb, gt_eval, det, seed_xy=seed.click_xy, tpl_xy=seed.template_xy, ...)
                                                                 viz.py:19
     Draws: GT boxes (red = mitotic, yellow = look-alike), one circle per detection
     coloured by its bucket, a white/black star at the raw click, a cyan '+' at the D8
     template centre. `radius` is drawn at the true evaluation match radius, so what the
     picture shows is exactly what the scorer counts as a hit.

viz.draw_box(ax, box, ...)                                      viz.py:68
     Draws the base_size search-template footprint square around the template centre.
```

---

## Component index

| Module | Key functions | Path |
|---|---|---|
| Data loading | `load_annotations`, `load_roi`, `roi_mpp`, `image_annotations` | [midog_utils/dataset.py](../midog_utils/dataset.py) |
| Seed selection | `build_seed`, `agreement_pool`, `border_filter`, `tighten_box_otsu`, `tightened_template_box` | [midog_utils/seed_selection.py](../midog_utils/seed_selection.py) |
| Production entry point | `run_production_pipeline`, `AXES`, module constants | [midog_utils/production.py](../midog_utils/production.py) |
| Channel conversion | `to_channel`, `to_hematoxylin_od` | [midog_utils/channels.py](../midog_utils/channels.py) |
| Search core | `find_and_suppress`, `FSConfig` | [midog_utils/find_and_suppress.py](../midog_utils/find_and_suppress.py) |
| Template matching | `read_padded_patch`, `build_augmentations`, `fused_response`, `robust_stats`, `extract_peaks` | [midog_utils/template_match.py](../midog_utils/template_match.py) |
| Suppression | `nms_by_distance` | [midog_utils/nms.py](../midog_utils/nms.py) |
| Chromatin axis | `hematoxylin_od`, `chromatin_density`, `score_detections` | [midog_utils/chromatin.py](../midog_utils/chromatin.py) |
| Invariant checks | `check_nms_radius` | [midog_utils/invariants.py](../midog_utils/invariants.py) |
| Scoring | `radius_px`, `greedy_match`, `bucket_detections` | [midog_utils/evaluate.py](../midog_utils/evaluate.py) |
| Visualization | `overlay`, `draw_box` | [midog_utils/viz.py](../midog_utils/viz.py) |
| Usage example | end-to-end demo | [production_pipeline_demo.ipynb](production_pipeline_demo.ipynb) |

---

## Punch-list resolution (vs. `PIPELINE_BREAKDOWN.md`, 2026-09-12)

| Old gap | Status in `production.py` |
|---|---|
| `MAX_PEAKS = 2,000,000` never actually set to D9's 100 | Resolved — [production.py:35](../midog_utils/production.py#L35) hardcodes 100 |
| `FSConfig.tm_method` defaults to `TM_CCOEFF_NORMED`, not D1's `TM_CCOEFF` | Resolved — the default is `cv2.TM_CCOEFF` since `c066829`, and [production.py](../midog_utils/production.py) also passes `tm_method=TM_METHOD` explicitly |
| Stale "click-centred is production" notebook cell | Resolved — `pipeline_debug_visuals/seed_refinement_variants.ipynb` cell 0 carries a 2026-09-12 correction pointing to `D8_TEMPLATE_ANCHOR.md` |
| Per-notebook `AXES` dict, no shared registry | Resolved — [production.py:41](../midog_utils/production.py#L41) is the one shared dict |
| `check_no_cap` structurally vacuous at `max_peaks=100` | Resolved — `check_no_cap` removed in `c066829`; recorded as `info["max_peaks_binding"]` instead |
| D9/D5 caveats absent from user-facing docs | Resolved — the demo notebook's closing markdown cell states both explicitly |

## Things to know, not bugs to fix

1. **Resolved 2026-09-17: augmented template banks no longer leak the seed's own match.** The
   5.0 px post-NMS filter used to miss self-matches that moved under rotation/flip/scale
   (12.66 px on 013.tiff, 5.39 px on 233.tiff, with the harness's 8-template bank);
   blanking the template's own footprint before correlation covers both, verified. A
   narrower gap remains where a *different* nearby detection (not the seed's own match)
   sits within one match radius of the seed but outside the blanked square — see
   `SELF_HIT_MASKING_PLAN.md` sec 2, "Known cost, accepted".
2. **`OD_WINDOW = 51` px is the window `chromatin.py`'s own docstring argues against**
   (a 31 px window measured better on 6 of 7 ROIs, and at 51 px the statistic
   "substantially reads the neighbour rather than the object" for candidates within 25 px
   of each other). Not a bug — it faithfully reproduces the old notebooks' `od51` — but
   the opt-in `chromatin_od` axis ships the weaker of two already-measured window sizes.
3. **`inv.check_nms_radius` is currently tautological.** It compares `nms_radius`
   (computed in [production.py](../midog_utils/production.py) as `ev.radius_px(mpp)`)
   against the same `radius_px(mpp)` formula recomputed inside the check itself — it can
   only ever pass today. It's a guard against future drift (e.g. someone hardcoding a
   radius elsewhere), not live evidence of anything right now.
4. **`base_size` cannot exceed `patch_size` (73) only because `otsu_window` (51) bounds
   it.** `tightened_template_box` returns `odd(max(component height, component width))`
   from a 51 px window, so it can never exceed 51 in practice (observed: 47 and 51 on the
   two test ROIs). `build_augmentations` crops the rotated patch at
   `patch[c-half_base : c+half_base+1]` with `c = 36`; a `base_size` above 73 would make
   that slice run outside the array. Raising `otsu_window` without re-checking this bound
   would remove the safety margin silently.
