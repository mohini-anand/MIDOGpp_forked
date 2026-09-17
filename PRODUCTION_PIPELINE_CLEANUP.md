# Production pipeline cleanup — confirmed removals

**Date:** 2026-09-16. **Scope:** dead code reachable only from
`production_pipeline/run_pipeline.py` -> `midog_utils/{channels, dataset, evaluate,
production, seed_selection, viz}.py` -> `find_and_suppress.py`, `template_match.py`,
`nms.py`, `chromatin.py`, `invariants.py`. This file is the agreed-on punch list from a
conversation walking every file in that call graph in invocation order.

**Status: applied 2026-09-16 in `c066829`**, with production output verified bit-identical,
except `plant_and_recover`, which was kept. Line refs below are to the pre-cleanup code
(`c066829^`). Research code outside the production path that used removed code now breaks;
this was accepted.

**Standing assumption for this pass (not yet a repo decision):** the production ranker
is fixed to `chromatin_od`. This directly overrides [DECISIONS.md D5](DECISIONS.md#L276),
which currently states the opposite (`tm_score` is the production ranker, `chromatin_od`
is a reported second axis only). Nothing below acts on that assumption yet — the
`rank_key`/`AXES` machinery in `production.py` and the `--rank-key` flag in
`run_pipeline.py` are explicitly **not** on this list (see "Deferred" below) because the
ranking decision itself is still open.

**Verification status:** independently re-checked against the live codebase by a
separate agent on 2026-09-16, re-deriving reachability from scratch rather than trusting
this document's claims. Two items were caught as wrongly listed for removal
(`channels.to_gray_inverted`, kept; `template_match.method_sign`/`_SIGN`, removed only
together with `fused_response`'s sign flip) and corrected in place below; one
cosmetic misattribution (`LOOKALIKE`'s actual reader function) was also fixed. Everything
else in "Confirmed for removal" checked out with zero reachable callers.

---

## Confirmed for removal

### `midog_utils/__init__.py`
- Eager imports of `baselines` and `experiment` at [`__init__.py:14-17`](midog_utils/__init__.py#L14) — neither module is reachable from `run_pipeline.py`; both go away entirely (see "Whole modules" below).

### `midog_utils/dataset.py`
- `load_slide_metadata` — [dataset.py:90-102](midog_utils/dataset.py#L90)
- `canonical_scanner` and `SCANNER_ALIASES` — [dataset.py:30-31](midog_utils/dataset.py#L30), [:21](midog_utils/dataset.py#L21) (only consumer is `load_slide_metadata`)
- `check_invariants` — [dataset.py:105-119](midog_utils/dataset.py#L105)
- `points` — [dataset.py:129-139](midog_utils/dataset.py#L129)
- `roi_area_mm2` and `check_roi_scale` — [dataset.py:179-207](midog_utils/dataset.py#L179)
- `CATEGORY_NAMES` — [dataset.py:16](midog_utils/dataset.py#L16)

**Explicitly kept in this file:** `tumor_type` (column and its use elsewhere in the repo — only the throwaway selection in `run_pipeline.py` is unused, see "Deferred"), `LOOKALIKE` (read in `evaluate.py`'s `lookalike_attraction_rate` — load-bearing; `bucket_detections` itself only branches on `== MITOTIC`, corrected 2026-09-16), `MITOTIC`, `BOX_SIZE`, `UNANNOTATED_IMAGE_IDS`, `TUMOR_ALIASES`/`canonical_tumor`, `load_annotations`, `load_roi`, `roi_mpp`, `image_annotations`.

### `midog_utils/seed_selection.py`
- `foreground_filter` — [seed_selection.py:156-175](midog_utils/seed_selection.py#L156)
- `SeedInfo` and `pick_seed` — [seed_selection.py:178-224](midog_utils/seed_selection.py#L178) (their only caller is `experiment.py`'s `run_single_pass`, which is itself being removed — see "Whole modules")
- `tightened_base_size` — [seed_selection.py:227-253](midog_utils/seed_selection.py#L227)
- `build_seed`'s `recentre=False` branch — [seed_selection.py:374-376](midog_utils/seed_selection.py#L374) (production always draws with `recentre=True`, the D8 path)
- `tighten_box_otsu`'s `"multiotsu"` and `"headroom"` branches — [seed_selection.py:114-130](midog_utils/seed_selection.py#L114) (production always calls with `method="binary"`); `headroom_frac` param becomes dead as a direct consequence
- The `center_tolerance > 0` branch in `tighten_box_otsu` and `_nearest_label_within` — [seed_selection.py:66-86](midog_utils/seed_selection.py#L66), [:136-137](midog_utils/seed_selection.py#L136) (production always calls with `center_tolerance=0`)
- The `method`, `center_tolerance` and `headroom_frac` parameters of `tighten_box_otsu`, `tightened_template_box` and `build_seed`, and `build_seed`'s `recentre`. `Seed.recentred` is kept and is always `True`.

### `midog_utils/channels.py`
- `to_hematoxylin` — [channels.py:19-26](midog_utils/channels.py#L19)
- `to_rgb` — [channels.py:34-36](midog_utils/channels.py#L34)

**Explicitly kept:** `to_hematoxylin_od`, **and `to_gray_inverted`** — corrected 2026-09-16 after subagent verification caught that this was wrongly listed for removal. `production_pipeline/run_pipeline.py:104` calls `ch.to_channel(rgb, "gray_inverted")` directly, inside `run_pipeline_on_roi` itself, to build the structural channel `select_annotation`/`build_seed` gate the seed against (a different channel from `hematoxylin_od`, which is only used for the template-match search in `production.py`). Deleting it would raise `KeyError` on every invocation. The `to_channel`/`CHANNELS` dispatch therefore needs at least two live entries (`gray_inverted`, `hematoxylin_od`), not one.

### `midog_utils/find_and_suppress.py`
- `FSConfig.channel` field — [find_and_suppress.py:25](midog_utils/find_and_suppress.py#L25) (write-only: set by `production.py`, never read inside `find_and_suppress()`)
- `FSConfig.score_threshold` field — [find_and_suppress.py:32](midog_utils/find_and_suppress.py#L32) (only reachable when `deep_floor_z is None`; production always sets `deep_floor_z=-1.5`). `find_and_suppress` now raises `ValueError` when `deep_floor_z` is `None`.
- `FSConfig.max_detections` field — [find_and_suppress.py:38](midog_utils/find_and_suppress.py#L38) — confirmed 2026-09-16, structurally dead since `max_peaks=100` already bounds the candidate pool far below `max_detections=10**9`. **Conditional on also editing `find_and_suppress()` itself**: unlike `.channel`, this field *is* read, at [find_and_suppress.py:120](midog_utils/find_and_suppress.py#L120) (`centers, scores = centers[: cfg.max_detections], scores[: cfg.max_detections]`) — that truncation line must be deleted too (the post-self-hit-removal `centers`/`scores` just pass through unchanged), not just the field.

**Explicitly kept:** `scale_normalize` field and `template_match._robust_z` — dead at current settings but coupled to the augmentation-bank machinery, which is being preserved on purpose.

**Changed:** `FSConfig.tm_method` defaults to `cv2.TM_CCOEFF`.

### `midog_utils/template_match.py`
- `METHODS` dict — [template_match.py:22-29](midog_utils/template_match.py#L22)
- `method_sign`, `_SIGN` — [template_match.py:20](midog_utils/template_match.py#L20), [:32-40](midog_utils/template_match.py#L32) — **conditional on also editing `fused_response`** (added 2026-09-16): `fused_response` calls `method_sign(method)` unconditionally at [template_match.py:183](midog_utils/template_match.py#L183) to negate the two `cv2.TM_SQDIFF*` (distance) methods onto the same "higher = better" convention `TM_CCOEFF` already uses natively. D1 permanently pins `tm_method=cv2.TM_CCOEFF` (a similarity method, closed decision), so `sign` is always `1.0` in production — safe to delete `sign = method_sign(method)` and the `if sign < 0: res *= sign` block from `fused_response` itself as part of this same removal, then `method_sign`/`_SIGN` become genuinely unreachable. Unlike everything else on this list, this is a behavior-narrowing edit to a function we keep, not a pure deletion — its one real-world side effect is `tm_variant_sweep.py` (root-level, out of scope), which compares TM methods including the SQDIFF family and now crashes on `tm.METHODS`/`tm.method_sign`.

**Explicitly kept:** `build_augmentations`, `fused_response`'s multi-template max-fusion and `scale_normalize` branch, `robust_stats`, `extract_peaks`, `read_padded_patch`, `plant_and_recover` (listed here for removal, kept on 2026-09-16) — the augmentation-bank machinery, preserved even though production settings collapse it to one template.

**Changed:** the `method` default of `fused_response` and `plant_and_recover` is `cv2.TM_CCOEFF`.

### `midog_utils/chromatin.py`
- `rerank` — [chromatin.py:53-69](midog_utils/chromatin.py#L53) (module's own docstring already states nothing in the production pipeline calls it; confirmed zero callers)

### `midog_utils/invariants.py`
- `check_tissue_mask_covers_gt` — [invariants.py:43-84](midog_utils/invariants.py#L43)
- `check_min_separation` — [invariants.py:136-152](midog_utils/invariants.py#L136)
- `check_no_cap` — [invariants.py:21-40](midog_utils/invariants.py#L21) — confirmed 2026-09-16, paired with the `FSConfig.max_detections` removal above. **Conditional on also editing `production.py`**: its only reachable call is [production.py:56](midog_utils/production.py#L56) (`inv.check_no_cap(int(info["n_detections"]), caps=(cfg.max_detections,), label="production_pipeline")`) — that line must go too.

**Explicitly kept:** `check_distinct_seeds` — [invariants.py:87-121](midog_utils/invariants.py#L87) (still needed for ongoing multi-seed experiments), `check_nms_radius`, `InvariantError` (still the base class `check_nms_radius` raises).

### `midog_utils/evaluate.py`
**Audit complete 2026-09-16** — all 12 candidates from the file's line-by-line review now
decided.

Confirmed for removal:
- `froc` — [evaluate.py:117-124](midog_utils/evaluate.py#L117)
- `sensitivity_at_fp` — [evaluate.py:127-132](midog_utils/evaluate.py#L127) (only ever consumed `froc`'s output)
- `optimal_assignment_disagreement` — [evaluate.py:182-208](midog_utils/evaluate.py#L182) (validates the greedy-vs-optimal matching *strategy*, not TP/FP determination itself — that logic is `greedy_match`/`bucket_detections`, both staying)
- `full_list_breakdown` — [evaluate.py:211-232](midog_utils/evaluate.py#L211) ("for now" — simple to rewrite later if wanted back)
- `threshold_sweep` — [evaluate.py:263-280](midog_utils/evaluate.py#L263) (also hardcoded to the `"score"`/tm_score column, and models an absolute-cutoff selection mode the production pipeline never uses)

Confirmed to keep, unchanged: `coverage_fraction` ([:56-66](midog_utils/evaluate.py#L56)), `recall_at_k` ([:108-114](midog_utils/evaluate.py#L108)), `_found_within` ([:135-139](midog_utils/evaluate.py#L135)), `lookalike_attraction_rate` ([:142-154](midog_utils/evaluate.py#L142)), `topk_composition` ([:157-165](midog_utils/evaluate.py#L157)), `recall_by_agreement` ([:168-179](midog_utils/evaluate.py#L168)).

Confirmed to keep, **but requires an edit**: `evaluate_run` — [evaluate.py:235-260](midog_utils/evaluate.py#L235). It calls `froc` (:241), `sensitivity_at_fp` (:248), and `full_list_breakdown` (:255) internally, and its return tuple's last element is `froc`'s `(fp_mm2, sens)` curve — all of which are being removed above. Same pairing pattern as `method_sign`/`fused_response`. The edit:
- delete the `fp_mm2, sens = froc(...)` line
- delete the `metrics.update({f"sens@{k}fp_mm2": ...})` line
- delete the `metrics.update(full_list_breakdown(det_out, gt_out))` line
- change `return det_out, gt_out, metrics, (fp_mm2, sens)` to `return det_out, gt_out, metrics` — a real return-signature change (4-tuple to 3-tuple), not just an internal simplification, since any future caller of `evaluate_run` needs to stop unpacking a curve that no longer exists.

### `midog_utils/viz.py`
Reviewed 2026-09-16.

Confirmed for removal:
- `froc_plot` — [viz.py:93-106](midog_utils/viz.py#L93) — zero production callers, and its only data source (`evaluate.froc`) is already being removed, so there'd be nothing left upstream to feed it even if kept.
- `legend_handles` — [viz.py:80-90](midog_utils/viz.py#L80) — confirmed 2026-09-16. Zero production callers (the overlay figure currently ships with no legend); was a judgment call rather than a clear dead-code case (a correct, ready-made companion to `overlay()` if a legend is ever wanted), decided in favor of removal.

Confirmed to keep: `overlay` ([:18-61](midog_utils/viz.py#L18)) and `draw_box` ([:64-77](midog_utils/viz.py#L64)) — both called directly by `production_pipeline/run_pipeline.py`'s `visualize()`; `BUCKET_COLORS` ([:11-15](midog_utils/viz.py#L11)) — the color dict `overlay` reads, stays with it. Minor, repo-wide-unused nit not worth its own line item: `draw_box`'s `label` parameter is never passed a value anywhere, including in every notebook that defines its own copy.

### Whole modules removed entirely
- `midog_utils/baselines.py` (179 lines) — zero callers in the production graph
- `midog_utils/experiment.py` (591 lines) — zero callers in the production graph; also the only caller of `seed_selection.pick_seed`/`SeedInfo` above
- `midog_utils/compare.py` (259 lines) — `Arm`/`evaluate_arms`; not even imported by `midog_utils/__init__.py`, and the production walkthrough doc states directly that nothing in the pipeline calls it

### Explicitly preserved, untouched
- `midog_utils_full/` — the entire directory, kept as a historical/rollback reference in case some removed functionality needs to come back into production later. Not part of this cleanup pass.

---

## Deferred — explicitly not decided yet

- **`production.py`'s `rank_key`/`AXES` machinery** and **`run_pipeline.py`'s `--rank-key` CLI flag** — the chromatin_od-only assumption is not yet a settled decision; another re-ranker may be added, so the selectable-axis plumbing stays for now.
- **`run_pipeline.py`'s `tumor_type` selection** at [run_pipeline.py:97](production_pipeline/run_pipeline.py#L97) — resolved 2026-09-16: kept as-is, not a pending item. Confirmed unused *within this file* (read into `meta_ix`, never referenced again — not in the returned dict, not in `summarize()`'s output row, not printed anywhere), but kept deliberately for future domain-grouped reporting, matching how the `tumor_type` column is already used across the rest of the repo (e.g. the `*_by_domain.csv` result files).
