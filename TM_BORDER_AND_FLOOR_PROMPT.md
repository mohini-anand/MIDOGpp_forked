# Two fixes to the template-matching step, tested as arms on 49 ROIs x 6 clicks: agent prompt

**Status:** spec only, written 2026-10-01. Nothing here has been run. You are a fresh agent; everything you need is in
this file, the repo, and the files it names. Read `DECISIONS.md` entries D1, D3, D5, D6, D7, D9 and D11 before writing
code. Where this prompt and the repo disagree, the repo wins and you report the disagreement.

## 0. What this is

The production search (`midog_utils.production.run_production_pipeline`, D5) cuts a template from the clicked
mitosis, correlates it over the hematoxylin optical-density channel with `cv2.TM_CCOEFF`, and ranks the peaks. Two
mechanisms in that step produce false positives that have nothing to do with the cell:

1. **The ROI border.** The search channel is replicate-padded by half a template so the template can be centred on
   every pixel. Replication turns any dark pixel on the boundary row into a dark streak into the pad, and the
   unnormalised matcher, which is mostly a mass detector, scores the streak. Measured on the recorded `mo_pad61` run:
   the 2 px band at the edge is 0.14 % of ROI area and holds 0.94 % of detections, a 6.8x enrichment; 252 of those 273
   detections sit on the boundary row itself and 93 % are unlabelled. In the chromatin-od top 10 the band takes
   3.9 % of slots against 0.45 % at ranks 30 to 100.
2. **The white surround.** Under `TM_CCOEFF` the mean-filled template's response is exactly a shape term plus a mass
   term, and the mass term is proportional to the box mean minus the ring mean. A lighter ring raises the score as
   much as a darker box, and glass is the lightest ring there is. In the top 10 by score, 24 % of non-mitoses have
   more than 20 % white in their 61 px neighbourhood against 11 % of mitoses; in mast cell tumour it is 82 % against
   52 %. The D11 blank adds a local version: it fills a 51 px square at the seed with the channel minimum, which is
   white, so every candidate within 30 px of the seed sees a white patch in its ring.

**Your job** is to implement one fix for each, test them as arms against the recorded search on the same 294 clicks,
and report whether they remove the mechanism and what they do to precision. You do not change production and you do
not write a DECISIONS entry; you produce the evidence.

**Two constraints you must not break.** D1: the matcher stays `cv2.TM_CCOEFF`; no normalised variant. D3: no global
rescale, stretch or top clip of the channel. The floor clip below touches only the bottom of the OD range, which is
why it is admissible to test; say so in the write-up, and note that adopting it would need a dated amendment to D3.

## 1. The two fixes, exactly

Let `hem` be the ROI's `hematoxylin_od` channel, `S = 61` the `mo_pad61` template side, `pad = 30`.

**Fix 1, border band invalid (search-side).** Build the padded channel and the fused response exactly as production
does. Then, before `robust_stats` and before `extract_peaks`, set `valid` to False for every pixel within `pad` px of
the ROI edge. No scored position then reads a padded pixel. Consequences you must state: the deep-floor median and
MAD are computed over the interior only; a band detection can no longer suppress an interior neighbour in NMS, so
interior pools may gain a detection; the band holds 1.0 % of annotated mitoses, which become unreachable in this
benchmark. In the product there is no ROI border, so this is a benchmark correction. Chromatin od needs no change
under this fix, because no detection lies within its 25 px half-window of the edge.

**Fix 2, floor clip (channel-side).** Choose a per-ROI floor `b` and form `hem_c = np.maximum(hem, b)`. Use `hem_c`
for cutting the template, for the blank and for the search. Three details:
- The refinement box (`refine_template_box`) runs on the **unclipped** `hem`, as production does, so the box and its
  centre are unchanged; the tripwire in Section 5 asserts that.
- `blank_seed_square` fills with `channel.min()`. On `hem_c` that minimum is `b`, so the blank fill becomes the floor
  with no code change. Assert that the fill value equals `b` and that `n_blanked` is unchanged.
- Chromatin od stays on the **unclipped** channel with its production replicate pad, so the ranking axis is
  identical across arms and only the generator differs.

The primary floor is `b_otsu`, the lower threshold of `skimage.filters.threshold_multiotsu(hem, classes=3)` over the
whole ROI. The secondary floor is `b_med`, the ROI median of `hem`. Record both per ROI, the share of pixels each
clips, and the ratio of each to the click's own template fill value (`mean(outside)` of the recorded `mo_pad61`
template). If `b_otsu` clips more than 60 % of a ROI's pixels, flag it and say what that ROI looks like.

## 2. Arms

All arms use the recorded clicks, the recorded refinement box, the `mo_pad61` template geometry, `MAX_PEAKS = 100`
before NMS, NMS at `ev.radius_px(mpp)`, and both ranking axes (`tm_score`, `chromatin_od`).

| arm | border band | floor | role |
|---|---|---|---|
| `A0` | production (replicate pad, band valid) | none | the recorded search; must reproduce it exactly |
| `A1` | invalid | none | fix 1 alone |
| `A2` | production | `b_otsu` | fix 2 alone |
| `A3` | invalid | `b_otsu` | both; the candidate production change |
| `A2m` | production | `b_med` | floor sensitivity, secondary |

Runtime: about 3.3 s per search; 5 arms x 294 clicks is roughly 80 to 90 minutes plus ROI loads. Smoke first.

## 3. Pre-registration, written in the notebook before any search runs

**Mechanism checks** (each is a yes/no, reported per arm):
- M1. Detections with their centre on the boundary row or column (`edge_px == 0`): `A1` and `A3` have 0; `A0` has the
  recorded 252-row count on the 294 clicks.
- M2. Share of top-10 non-mitoses with more than 20 % white in their 61 px neighbourhood (`gray > 215`, computed on
  the unclipped RGB): falls from `A0` to `A2` and `A3`, with a ROI-cluster bootstrap CI on the paired difference.
- M3. Lumen-adjacent true positives: the number of mitoses in the top 30 with more than 20 % white, summed per ROI.
  Report the paired `A2 − A0` and `A3 − A0` with CIs. A significant loss here is a red flag that goes next to the
  primary result; it is not a stopping rule.
- M4. Near-seed leaks: non-TP detections whose centre is between 1 and 2 match radii from the click, per arm. The
  floor fill should reduce them from `A0` to `A2`.
- M5. Synthetic demonstrations, run before the real data (Section 4).

**Primary family**: `A3 − A0` at `tm_score`, K = 10, 20, 30. Holm within the family. The D5 bar is a bootstrap CI
excluding 0 and a majority of ROIs moving the delta's way.

**Secondary families**, each with its own Holm: the decomposition `A1 − A0`, `A2 − A0`, `A3 − A1`, `A3 − A2` at
`tm_score`; the same five contrasts at `chromatin_od`; `A2m − A2` at both axes.

**Also reported, descriptive**: mitoses in the whole pool (83 to 100 detections) per arm, since that is the generator
metric; per-domain tables with mast cell and lung called out; the per-click SD of TP@10; the worst click per ROI.

**Decision rule, fixed now.** Recommend `A3` for a DECISIONS entry only if all hold: `A3` beats `A0` at K = 10 on
`tm_score` under the D5 bar with a positive pooled delta; no K on either axis shows a significant loss; M1 holds;
M2 falls with a CI excluding 0. If M3 shows a significant loss, recommend `A1` alone and report `A2` as not adopted.
Every number is in-sample on the 49 ROIs; say so.

## 4. Adversarial checks, yours, logged as you go

Keep `adversarial_log.csv` with `item, concern, test, result, action`. Minimum contents:

**Before coding**, for each fix write what could make it fail or mislead, then test it:
- The clip raises the deep-floor median and shrinks the MAD, so the threshold moves. Report the threshold and the
  number of peaks before NMS per arm; if `A2` returns systematically fewer peaks, say what that does to the pool.
- In dense tissue the clip touches almost nothing, so `A2` equals `A0` there and the effect is confined to loose
  tissue. Report the share of pixels clipped per domain and check that the `A2 − A0` effect tracks it.
- A per-ROI floor is a per-scanner floor. State whether `b_otsu` is stable within a domain.
- Invalidating the band changes the robust statistics. Report median and MAD for `A0` and `A1` per ROI.
- The two fixes could interact: with the band invalid, the clip's remaining effect is only interior white. Check that
  `A3 − A1` is close to `A2 − A0` minus the band's share.

**Synthetic demonstrations, mandatory, before the real run.** On a flat synthetic channel at tissue level:
- A dark disc 30 px from a straight white half-plane versus the same disc in uniform tissue: show the recorded
  template's `TM_CCOEFF` peak is higher beside the white under `A0` and equal within 5 % under the clip. State the
  ratio.
- A dark disc bisected by the array edge with replicate padding: show the peak forms on the boundary row and
  exceeds the peak of the same disc placed 60 px inside; show the band invalidation removes it.
- A disc 20 px from a 51 px square filled with the channel minimum versus with `b`: show the response difference.
These are the proof that each fix addresses its mechanism, independent of whether precision moves.

**After the run:** read your own tables for anything that contradicts the mechanism story, and say so first.
Re-match three clicks by hand from the JSON to confirm the matching. Open the montages (Section 6) and write one
sentence per arm on what the top-10 false positives now look like on `300.tiff` (mast cell), `201.tiff` (lung) and
`546.tiff` (melanoma, a stuck ROI).

**The rule for disagreement.** A gate failure (Section 5) means stop and report; do not adapt the design to make it
pass. A finding that one of the two fixes is wrongly specified means implement what the mechanism requires, log the
change, and keep the specified version as an extra arm if it is cheap.

## 5. Gates (halt and report)

1. **Config.** `CHANNEL`, `TM_METHOD`, `PEAK_MIN_DISTANCE`, `DEEP_FLOOR_Z`, `MAX_PEAKS`, `OD_WINDOW` match
   `midog_utils.production`; `ev.MIDOG_RADIUS_UM == 7.5`; `PAD_PX == 2` in the copied refinement and no `PAD_UM`.
2. **Tripwire.** The copied `refine_template_box` on the 294 clicks reproduces the recorded `bs`, `tcx`, `tcy` and the
   `A0` template SHA-1 from `results/precision_at_k_49roi_6seed_mopad_postD11_templates.csv` (`condition == 'mo_pad61'`)
   on all 294. Halt on any difference.
3. **`A0` parity.** Every `A0` search reproduces the recorded pool: `cx`, `cy`, `score` (float32) and `od`
   (`rtol 1e-12`) equal the first `n_detections_recorded` rows of
   `results/deep_pool_mo_pad61_49roi_6seed_postD11_top1000.csv`, and `tp_at_10/20/30` on both axes equal
   `results/precision_at_k_49roi_6seed_mopad_postD11_per_run.csv`. This is the proof that the arms differ only by the
   fix.
4. **Clip identity.** `hem_c >= hem` everywhere and `hem_c == hem` wherever `hem >= b`; the blank fill equals `b`;
   `n_blanked == 51 ** 2`.
5. **Band identity.** In `A1` and `A3`, `valid` is False exactly within `pad` px of the edge, and no detection has
   `edge_px < pad`.
6. **Matching.** Re-run `ev.bucket_detections` on every arm's list (D6), click's own annotation dropped. Each list is
   a permutation of its own pool.
7. **Statistics.** Copy cells 8 and 9 of `production_hematoxylin_only/multiotsu_pad_arms/multiotsu_pad_arms_49roi_6seed.ipynb`
   verbatim; the fill-geometry run's 10 recorded `delta_stats` rows must reproduce, and `pair_stats_clicks` must equal
   `pair_stats` at 6 clicks per ROI. Bootstrap ROI order is the order of first appearance in
   `results/precision_at_k_49roi_6seed_fillgeom_postD11_per_run.csv`; `N_BOOT 10000`, `BOOT_SEED 20260917`.
8. **Self-matches.** `selfmatch_counts` on every list; if any arm has a self-match in its top 10, drop it, re-match,
   and rerun the contrasts, reporting both.

## 6. What to copy and from where

- `refine_template_box` and its constants: verbatim from
  `production_hematoxylin_only/multiotsu_pad_arms/multiotsu_template_refinement_14roi.ipynb` at commit `c767898`
  (source SHA-1 must be `21c43a20db00dccf7ac3bab88aee423c345441b6`).
- Template, search, ranking, matching and flag helpers: `fill_outside`, `template_sha1`, `padded`, `read_window`,
  `arm_side`, `build_mo_templates`, `search_once`, `rank_by`, `score_top_k`, `template_flags`, `selfmatch_mask`,
  `selfmatch_counts` from `multiotsu_pad_arms_49roi_6seed.ipynb` cell 11, and the per-click steps of its cell 15.
  `search_once` is where fix 1 goes (the `valid` mask) and the channel you pass it is where fix 2 goes.
- A worked example of the same copy-and-gate pattern, including the `research_pool` function and the parity
  asserts, is `production_hematoxylin_only/multiotsu_pad_arms/candrefine_rerank_mo_pad61_select.ipynb` cells 1 to 12.
- The white fraction: a 61 px window of `cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)` cut from a replicate-padded copy,
  share of pixels above 215.
- Ground truth: `ds.load_annotations('../../databases/MIDOG++.json')`, `ds.image_annotations(annotations, file_name)`,
  both categories, the click's `seed_ann_id` dropped. Images: `{'original_14': '../../images/extra_valid',
  'testing_35': '../../images/extra_valid/testing_set'}[subset]`, read with `ds.load_roi`, `mpp` from `ds.roi_mpp`.
- Montages: for `300.tiff`, `201.tiff` and `546.tiff`, click 0, the top 10 of each arm at `tm_score` as 61 px RGB
  tiles framed by bucket.

## 7. Outputs

Notebook: `production_hematoxylin_only/multiotsu_pad_arms/tm_border_floor_49roi_6seed.ipynb`.

`results/tm_border_floor_49roi_6seed_{per_run,top30,per_roi,delta_stats,floors,edge,white,nearseed,adversarial_log}.csv`:
- `per_run`: one row per click x arm x axis with `tp_at_10/20/30`, `n_peaks`, `n_detections`, threshold, median, MAD,
  `n_edge0`, `n_edge30`, mitoses in the pool, self-match counts.
- `top30`: the top 30 per click x arm x axis with `cx`, `cy`, `score`, `od`, `bucket`, `edge_px`, `white_frac_61`.
- `floors`: per ROI, `b_otsu`, `b_med`, share clipped by each, domain, scanner if known, median template fill value.
- `edge`, `white`, `nearseed`: the M1 to M4 tables with CIs.
- `delta_stats`: every contrast, both axes, with CI, Holm p, W/L/T and the D5 bar.
Montages beside the notebook as `tm_border_floor_montage_<roi>.png` [files since removed; the three montages are embedded in `multiotsu_pad_arms/tm_border_floor_49roi_6seed.ipynb`, cell 31]. Smoke runs write `_smoke` copies
(`TMFIX_SMOKE_ROIS=337.tiff,509.tiff,201.tiff`). The full run refuses to start if any output exists; all outputs are
written at the end.

The final markdown cell holds, in order: the decision-rule outcome; the M1 to M5 results; the primary and secondary
tables; the per-domain table; the adversarial log; what you changed in this prompt and why; what you could not check.

## 8. Environment and conventions

- Interpreter: `/Users/mohinianand/anaconda3/bin/python3`, the only one that imports `cv2` and `midog_utils`.
- Execute with `/Users/mohinianand/anaconda3/bin/jupyter nbconvert --to notebook --execute --inplace
  --ExecutePreprocessor.timeout=14400 <path>`. The notebook sits two levels below the repo root:
  `sys.path.insert(0, '../..')`.
- Docstrings: the opening `"""` on its own line, a brief summary, one line per parameter with type and a short
  description, one line for the return. Function definitions and calls on a single line, never wrapped.
- Do not commit. Do not modify or delete any existing file, `midog_utils` included; the fixes live in the notebook's
  copies of the search helpers. Kill only a PID you started.
- Report back with: the decision-rule outcome in one line; M1 to M5; the primary table; the lumen-adjacent TP
  result; the per-domain table; the adversarial log. Lead with anything that contradicts the mechanism story.

## 9. Not in this prompt

- Any change to `midog_utils`, `production.py` or `DECISIONS.md`.
- `TM_CCOEFF_NORMED`, a background-subtracted channel, a tissue-masked ring, rotations, or any other generator change.
  If the clip is insufficient, the tissue-masked ring via the shape-plus-mass decomposition is the next candidate;
  name it in the write-up, do not build it.
- Re-ranking, feature extraction, fresh clicks or held-out ROIs. The feature-extraction prompt
  (`CANDIDATE_FEATURE_EXTRACTION_PROMPT.md`) runs after this, on whichever search is adopted.
