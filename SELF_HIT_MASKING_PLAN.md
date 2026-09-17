# Replace the self-hit filter by blanking the seed's template footprint before correlation, in `midog_utils`

**Status:** design decided, verified independently, gates built — **execution-ready**. The
patch, fresh bit-identical-by-design reference captures, a design-specific `check`/`reference`
script and negative tests all exist, at `../cleanup_harness/selfmask/image_blank_design/`.
**Target:** the production package `midog_utils/`, not yet touched. Numbers in §3a/§3b are
from real simulation, run 2026-09-17 against `midog_utils/` as committed in `c066829`
(unchanged through `c7cf7ab`, the HEAD this plan was verified against). §3c's gaps were closed,
and §3d's corrections found, by an independent verification session the same day.

**`DECISIONS.md` D11** now records this design as decided, status "not yet applied," and
amends D10 to point at it. Once every gate in §6 passes for a real application, append the
dated amendment described there — don't edit D11's body otherwise.

**Revision history:**
- **Earlier draft** (targeted `midog_utils_full/`): masked a ±7 px box on both the score map
  and `valid`. Rejected because it costs precision. Archived at
  `../cleanup_harness/selfmask/SELF_HIT_MASKING_PLAN.midog_utils_full_draft.md`.
- **2026-09-16:** a disc of one NMS radius (7.5 µm) cleared from `valid`, after the
  deep-floor threshold, before peak extraction. Verified, gated, fully built — patch,
  bit-identical reference captures, `check`/`reference` scripts, 9 negative tests,
  augmented-bank testing. Never applied to production. Archived, unchanged and still valid
  for that design, at
  `../cleanup_harness/selfmask/SELF_HIT_MASKING_PLAN.valid_mask_disc_design.md`.
- **2026-09-17, superseding the above:** this design — blanking the seed's own template
  footprint out of a copy of the search image, before correlation, instead of excluding a
  region from `valid` after correlation. Built as a standalone simulation and measured on 3
  ROIs, then all 49 ROIs `production_hematoxylin_only/bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb`
  uses. §3a/§3b. Not yet independently verified or gated at this point.
- **2026-09-17, applied:** executed against `midog_utils/` at `e4b7884`; all gates passed and
  `DECISIONS.md` D11 carries the dated amendment. Two defects in §5's own text were found
  during application:
  - §5.1 Edit 10's replacement text was affirmative inside a negation list — **corrected at
    source above**, so a re-application gets the fixed text.
  - §5.1 Edit 7 adds the `blank_seed_square` node but not the two downstream call-diagram
    nodes, which keep naming the unblanked `hem` where the patched code correlates the blanked
    copy. **Not corrected at source** — the three repairs are recorded only in the execution
    log, so a re-application must redo them by hand.
- **2026-09-17, later the same day (this revision):** independently re-verified in a separate
  session — mechanism claims re-derived analytically and on real data; §3a/§3b's numbers
  reproduced via a from-scratch prototype (not a re-execution of the design session's own
  scripts); all three of §3c's evidence gaps closed (`tm_score` at scale, the augmented bank,
  seeds 1-4); one bug found and fixed in the verification session's own harness (blanking must
  run after template extraction, not before — see §3d) and one in its own augmented-bank test
  (seed pinning); one real, non-production-only mechanism found, explained and shown absent
  from production's actual configuration (§3d). The real code patch, fresh reference captures,
  a design-specific `check`/`reference` script and negative tests were then built in that same
  session (§4-§8). No further open gaps.

**Summary for whoever picks this up next:** read §1 (unchanged, still describes today's
committed code), §2 (the design, corrected where verification found the earlier wording too
strong) and §3 (the evidence, now complete) to understand what's decided. Then follow §6 to
apply it — a mechanical Steps-0-through-7 procedure. §4/§5 are the actual code and doc edits;
§7 is troubleshooting; §8 is the audit trail for how everything below was produced.

---

## 1. What the self-hit filter does today

`find_and_suppress` (`midog_utils/find_and_suppress.py`) runs three steps:
1. `extract_peaks` keeps the top `max_peaks=100` local maxima.
2. `nms_by_distance` suppresses within `nms_radius` = `evaluate.radius_px(mpp)` (7.5 µm,
   29.6-33.1 px).
3. A post-NMS filter drops detections within `FSConfig.self_hit_radius = 5.0` px of
   `seed_xy`.

It records `info["max_peak_score"]`, `info["n_self_hits"]` and `info["seed_self_score"]`.
`production.py` sets `SELF_HIT_RADIUS = 5.0`.

**What it actually achieves.** Measured on the 14 `images/extra_valid` ROIs at production
config, `seed_index=0`:
- The seed's own correlation peak enters the 100-peak pool in 14 of 14 runs.
- It is the pool maximum in 12 of 14 (not 094.tiff or 301.tiff).
- Because it outscores its neighbourhood, greedy NMS removes every other peak within one NMS
  radius of it *before* the 5.0 px filter drops it.

**So the real effect today is "no detection within one NMS radius of the seed."** The 5.0 px
constant only decides whether the peak itself gets caught.

**Where it fails:**
- **The seed's own peak falls below the `max_peaks` cutoff.** Nothing then suppresses the
  neighbourhood. On 301.tiff at `seed_index=1`, a secondary peak of the clicked cell 26.6 px
  (6.7 µm) away stays in the list — re-verified 2026-09-17 at 26.62 px, present unchanged in
  both today's committed code and this design's patched code (§3d: this is the one case where
  this design does not close a gap D10's design would).
- **An augmented template bank moves the self-match off-centre.** With the cleanup harness's
  8-template bank (`scales=(0.8,1.2)`, `n_angles=2`, `flips=(False,True)`, plain
  `seed_index=0` draws): 013.tiff keeps a self-match 12.66 px away at rank 7 (score 20.196);
  233.tiff keeps one 5.39 px away at rank 0 (score 12.246). Both reproduced exactly in the
  2026-09-17 verification session; both are removed by this design (§3c).

## 2. The design, and why

**Change:**
- **What:** blank the refined template's own footprint — a `base_size × base_size` square
  (the same size `seed_selection.tightened_template_box` returns for this click, an odd
  number, 23-51 px in the data measured), centred on the rounded template centre, using the
  same `int(round(...))` convention `read_padded_patch` uses — in a *copy* of the search
  channel (`hematoxylin_od`), filled with that channel's own whole-ROI minimum. The minimum is
  computed from the *original*, unblanked channel, over the whole ROI, before any blanking
  happens.
- **When:** the template is cut from the original, unblanked channel first
  (`read_padded_patch` + `build_augmentations`, matching production's single-template
  config). Only *after* that does blanking happen, on a separate copy of the channel used for
  the search; correlation (`fused_response`), the deep-floor threshold (`robust_stats`), peak
  extraction and NMS all run against the blanked copy. Ranking (`chromatin_od`'s
  `chromatin.score_detections`) uses the *original*, unblanked channel — a candidate's
  surrounding tissue should score on what it actually looks like, not on an artificial flat
  patch. This ordering is load-bearing: blanking the channel *before* cutting the template
  (rather than after) would also blank the template itself, since production's `img_channel`
  is the single array both steps read from — verified the hard way in the 2026-09-17
  verification session, whose first prototype attempt did exactly that and produced a fully
  flat (all-zero) correlation map everywhere, not just near the seed (§3d).
- **Delete:** the post-NMS filter, `FSConfig.self_hit_radius`, `production.SELF_HIT_RADIUS` —
  same deletions as the superseded design.
- **No `valid`-masking.** This design doesn't touch `valid` at all; the seed's own match
  simply never forms as a correlation peak, because its source pixels are gone before
  `matchTemplate` ever runs.

**Why blank the image, when the superseded design chose not to.** The earlier, rejected
score-map-flooring draft (archived plan's §2) showed that touching the *output* of
correlation creates a false peak on the slope of the seed's own match, right where
`extract_peaks`'s dilation test looks — a measured, real cost. Blanking the *input* image was
suspected of a similar failure by analogy, but this was never itself measured until it was
built — the mechanism is different from score-map-flooring: `TM_CCOEFF` against a flat window
naturally computes near zero rather than an extreme value, and because matching is a sliding
window, the corrupted response tapers off with the window's shrinking overlap instead of
stepping sharply at a boundary the way a post-hoc floor does. **Independently re-derived,
2026-09-17:** analytically (a window's own local mean equals the fill constant for a
perfectly flat window, so the mean-subtracted term in `TM_CCOEFF`'s definition is exactly
zero) and numerically at hematoxylin_od's real magnitude (a synthetic flat-window response of
~3e-7, five orders of magnitude below background noise — the naive first attempt at this
check, using an arbitrary large synthetic magnitude, produced a much larger nonzero result
that was pure float32 catastrophic-cancellation noise from the wrong scale, not a real
mechanism failure); and on real ROI data (403.tiff: the nearest pre-NMS peak to the blanked
square sits 299 px away — no boundary-artefact detection anywhere measured, on any of the 49
ROIs, the augmented bank, or seeds 1-4).

**Known cost, accepted: a real detection can survive near the seed, outside the blanked
square.** Blanking only the template's own footprint is smaller than blanking a full
`nms_radius`/match-radius disc (the superseded design's choice) — `base_size` ranges 23-51 px
in the data measured, while `match_radius` is 29.6-33.1 px, so the blanked square is often
*narrower* than the neighbourhood the seed's dominant peak used to sweep clean via NMS.
Measured on 2 of 49 ROIs (§3b) and again under the augmented bank (§3c): a detection survived
within one match radius of the click but outside the blanked square. Mechanism, confirmed via
peak counts, not just distance: in production, the seed's own peak is the dominant score in
its neighbourhood, so greedy NMS eliminates every weaker peak within `nms_radius` of it before
anything else runs; once that dominant peak is prevented from forming at all, a weaker-but-
genuine local correlation maximum just outside the (too-small) blanked square has nothing left
to suppress it, and survives NMS on its own merits. This is not a bug in NMS — it behaves
identically throughout; what changed is what was available for it to compete against.
**A design that blanked a full match-radius region instead of the template's own footprint
would close this gap entirely** — that was the superseded design's own choice of mask size,
for exactly this reason. The user chose the literal template footprint anyway, having seen
this evidence: in every measured case the leaked detection matched no ground truth
(`non_human_findings`), and cost no TP@K in any run measured (§3b, §3c). This is a structural
property of the design, not a rare fluke — whenever `base_size < match_radius`, the gap
exists by construction — so it should be expected to recur, not treated as a one-off.

**Named instance where this design is strictly weaker than the one it supersedes.** The
301.tiff `seed_index=1` case in §1's "Where it fails" is the superseded design's own headline
fix target — D10's full-`nms_radius` disc (30-33 px) covers a secondary peak 26.6 px from the
seed; this design's `base_size=23` px square (half=11) does not, since 11 < 26.62 < 29.61. Both
were re-verified 2026-09-17: the row is present, at the same distance, in both the current
production code and this design's patched code. It costs no TP@K in the data measured, but it
is real and should not be described only through the 289.tiff/548.tiff/245.tiff instances
above — this is the case where reverting to D10's disc would concretely help.

**The TP@K movement itself is not specific to this design.** Freeing the seed's own slot in
the `max_peaks=100` pool — whether by blanking the image or by clearing `valid`, as the
superseded design does — lets one different candidate into the pool, and under `chromatin_od`
re-ranking that candidate can land anywhere, including inside the top K. §3b found one case
where the newly-admitted candidate is a look-alike (`HUMAN_REJECTED_LABEL`: matched an
annotation a pathologist reviewed and rejected as non-mitotic) that displaces a true positive
out of the top 10. The *identical* failure mode is already in the superseded design's own
pre-registered numbers (its archived plan's §3b: 459.tiff, seed 4, `chromatin_od`, K=50, TP
43→42, caused by a look-alike entering at rank 26 and pushing out a true positive at rank 49)
— **independently re-verified 2026-09-17 to be the exact same (ROI, seed, K) cell, at the
exact same TP values, for this design too**, along with three other cells (013.tiff seed 4
K20, 300.tiff seed 2 K50, 301.tiff seed 3 K30) — strong evidence the shared "freed pool slot"
mechanism behaves identically regardless of which masking approach frees the slot. This is a
property of using `chromatin_od` as the ranking axis together with a fixed `max_peaks` cap,
not of how the excluded slot got freed.

**`tm_score` — corrected, 2026-09-17.** Under production's actual configuration (full-ROI
search, single un-augmented template, `scale_normalize=False` — D5's actual default ranker),
a newly-admitted candidate reliably sorts last and never displaces anything in a bounded top K:
confirmed on every real run measured (49 ROIs, seeds 1-4 on 14 ROIs, the corrected
augmented-bank rerun). **This is not, however, a universal property of the mechanism**, and
the original wording ("a newly-admitted candidate always sorts last and can never displace
anything in a bounded top K") is too strong as a general claim. Under `scale_normalize=True`
(never used by production) combined with a small search crop (the 1024 px crop
`harness.py`'s own augmented-bank test uses for speed — never what production does, which
always searches the full ROI), blanking can perturb an *unrelated* peak's score by ~1e-3 via
`_robust_z`'s shared per-template normalisation statistics, occasionally flipping which of two
near-tied peaks NMS keeps first in a "chain" suppression geometry (peak A suppresses both B and
C, but B does not suppress C) — this moved `tm_score` TP@K by one (K=20, 2→1) on one
non-canonical seed draw for 245.tiff. Full mechanism and why it cannot occur under production's
`scale_normalize=False` configuration: §3d.

## 3. What changes in the numbers

### 3a. First check: 3 ROIs (245, 403, 246.tiff)

Real, unmodified `production_pipeline.run_pipeline.run_pipeline_on_roi(fn, seed_index=0,
rank_key="chromatin_od")` as the baseline, against a hand-reimplemented "blank" arm built
from the same production primitives (`template_match.read_padded_patch`,
`build_augmentations`, `fused_response`, `robust_stats`, `extract_peaks`,
`nms.nms_by_distance`, `chromatin.score_detections`, `evaluate.bucket_detections`), with a
third control arm (same reimplementation, but unblanked, running the *current* post-NMS
self-hit filter) added to prove the reimplementation is faithful before trusting the
experimental arm.

| ROI | control matches real production | precision@10/20/30 (production) | precision@10/20/30 (blanked) | n_detections (production → blanked) | closest detection to seed |
|---|---|---|---|---|---|
| 245.tiff | exact | 0.200 / 0.100 / 0.067 | 0.200 / 0.100 / 0.067 | 89 → 90 | 1082 px, rank 41, unannotated |
| 403.tiff | exact | 0.500 / 0.500 / 0.467 | 0.500 / 0.500 / 0.467 | 98 → 99 | 299 px, rank 22, unannotated |
| 246.tiff | exact | 1.000 / 0.700 / 0.700 | 1.000 / 0.700 / 0.700 | 96 → 97 | 467 px, rank 39, unannotated |

No precision change on any of the 3, no detection anywhere near the seed. Consistent with the
design working as intended; 3 ROIs is a small sample — the near-seed leak in §3b didn't show
up here at all.

### 3b. Full check: all 49 ROIs

Same method, scaled up, with seeds pinned via the real `ann_id=` parameter to the exact
annotation `production_hematoxylin_only/bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb`
selected for `condition='gray_bbox', seed_index=0` on each ROI (from a frozen snapshot, since
that notebook's live output CSV was being concurrently rewritten). That notebook's own
click-selection requires an annotation to pass a *joint* gate across all three of its tested
conditions, which can occasionally accept a different annotation than production's plain
`seed_index=0` draw would; pinning by `ann_id` sidesteps that divergence, at the cost of this
evidence technically being "production run on the notebook's seed-0 click" rather than
"production's own freshly-drawn seed-0 click." (**Verified independently, 2026-09-17, does
NOT always coincide** — see §3d's note on the augmented-bank re-check, where this distinction
mattered for two of the six ROIs tested.)

**Fidelity, checked before trusting anything else, and independently re-derived:**
- Seed-geometry sanity check (live `base_size`/`template_xy` vs. the notebook's stored values
  for the same annotation): 49 of 49 exact, both in the design session and independently.
- Hand-reimplemented control arm vs. real production, unblanked, old filter still active:
  49 of 49 bit-for-bit.

**Pooled TP@K, all 49 ROIs, production vs. blanked** (identical in both the design session's
own numbers and a from-scratch, independently-written prototype that calls `find_and_suppress`'s
own primitive functions directly, not the design session's hand-rolled orchestration):

| K | TP (production) | TP (blanked) | Δ |
|---|---|---|---|
| 10 | 270 | 269 | −1 |
| 20 | 460 | 461 | +1 |
| 30 | 589 | 590 | +1 |

**Per-domain pooled TP@K** (7 domains × 7 ROIs each): identical between arms at every K in
5 of 7 domains (canine cutaneous mast cell tumor, canine lung cancer, canine lymphosarcoma,
canine soft tissue sarcoma, human melanoma). The other 2:

| domain | K | TP (production) | TP (blanked) |
|---|---|---|---|
| human breast cancer | 10 | 45 | 44 |
| human neuroendocrine tumor | 20 | 72 | 73 |
| human neuroendocrine tumor | 30 | 90 | 91 |

**Every (ROI, K) cell that moved, with cause** (re-derived exactly, independently):

| ROI | K | TP before → after | cause |
|---|---|---|---|
| 128.tiff | 10 | 7 → 6 | a look-alike (`HUMAN_REJECTED_LABEL`, 3319 px from the seed) enters at rank 1 under `chromatin_od`, pushing a true positive from rank 9 to rank 10 — out of top-10 |
| 400.tiff | 20 | 15 → 16 | a genuine true positive (`HUMAN_CORRECT_LABEL`, 716 px from the seed) enters the freed slot and ranks inside top-20 |
| 400.tiff | 30 | 20 → 21 | same entering row, still inside top-30 |

**`n_detections` delta** (blanked − production) across all 49: +1 on 45 ROIs, +2 on 1 ROI
(289.tiff), +0 on 3 ROIs (123.tiff, 401.tiff, 546.tiff — mechanism corrected in §3d).

**The near-seed leak (§2's accepted cost), both instances found, re-derived exactly:**

| ROI | `base_size` | `match_radius` | leaked detection: distance / rank / bucket | reached top-30? | TP@K impact |
|---|---|---|---|---|---|
| 289.tiff | 25 px | 30.2 px | 24.4 px, rank 10 (0-indexed), `non_human_findings` | yes — inside top-20 and top-30, not top-10 | none (displaced a different non-TP row) |
| 548.tiff | 29 px | 33.0 px | 28.0 px, rank 56, `non_human_findings` | no | none |

On 289.tiff specifically, `n_detections` rose by 2: production's NMS reduced 100 pre-NMS peaks
to 99 survivors (one suppressed by the dominant self-match) then dropped 1 more via the old
filter, netting 98; the blanked variant's NMS kept all 100, netting 100 — two peaks freed, not
one, confirming this ROI hit an NMS-domination-collapse, not just an ordinary freed-slot.

**Zero bucket flips on shared rows, checked independently across all 98 runs (49 ROIs × 2 rank
keys).** Every TP@K movement above is pure displacement (a new row entering the top K and
pushing the K-th row out); no shared detection's ground-truth match changed as a side effect.

### 3c. Gaps closed (2026-09-17 verification session)

All three of the design session's own named gaps, closed:

**`tm_score`, all 49 ROIs, same pinned seeds:** zero TP@K change at every K. Matches the
reasoning in §2 (a freed candidate sorts last) exactly under production's actual
configuration — see §3d for the one non-production configuration where this does not hold.

**The 8-template augmented bank** (`scales=(0.8,1.2)`, `n_angles=2`, both flips, `harness.py`'s
own 1024 px crop methodology), on `harness.py`'s own `AUG_ROIS` (013, 245, 403.tiff) plus
233.tiff (a documented self-hit-radius failure `AUG_ROIS` omits), with plain `seed_index=0`
draws matching `harness.py`'s own methodology exactly:

| ROI | `scale_normalize` | seed (plain draw) | Arm A (old filter) near-seed row | Arm B (blanked) |
|---|---|---|---|---|
| 013.tiff | True | ann 254, base_size 31 | 12.66 px, rank 7, score 20.196 (matches `probes/p3.log` to 5 decimal places) | none |
| 233.tiff | True | ann 5761, base_size 25 | 5.39 px, rank 0, score 12.246 | none |
| 245.tiff | False | ann 6274, base_size 47 | none | one new leak at 24.2 px, `non_human_findings` — the same "known cost" pattern as 289/548.tiff, not a self-hit-radius-style failure |

Both documented historical self-hit-radius failures are removed by this design under the
augmented bank, exactly as claimed. No `tm_score` movement on any of the 4 ROIs × 2
`scale_normalize` values with the correct (plain-draw) seeds. (A first attempt at this check
used seeds pinned to the 49-ROI notebook's own draws instead of plain `seed_index=0` — a bug,
since all 14 canonical ROIs are a subset of the 49-ROI pins, so the intended fallback to a
plain draw never triggered. That run's 013.tiff result (a different seed, ann_id not matching
`p3.log`) is superseded by the table above; its one `tm_score` TP@K movement — real,
deterministic, and mechanistically explained — turned out to be on that same
mis-selected seed, not a canonical one; see §3d.)

**Seeds 1-4, 14 `images/extra_valid` ROIs, both rank keys, K=10/20/30/50:** 4 (ROI, seed, K)
cells move — 013.tiff seed 4 K20 (8→9), 300.tiff seed 2 K50 (46→45), 301.tiff seed 3 K30
(25→26), 459.tiff seed 4 K50 (43→42) — all `chromatin_od`, none `tm_score`. **All four match
the superseded (D10) design's own pre-registered seed-robustness numbers exactly, cell for
cell**, despite the two designs sharing no code — strong cross-validation that the shared
freed-pool-slot mechanism produces the same outcome regardless of which masking approach frees
the slot, in every case both were measured. One near-seed leak: 301.tiff seed 1, 26.62 px (see
§2's named weaker-point instance) — present unchanged in production's real code today, costs
no TP@K.

**The `n_detections` +0-on-3-ROIs loose end, resolved — and the design session's own guess was
wrong.** Independently re-derived: 401.tiff and 546.tiff have `n_self_hits=0` in the *current,
unpatched* code — the seed's own peak was never among the post-NMS survivors even before this
design is applied, so freeing its slot changes nothing, by construction. `max_peaks_binding`
was `True` on all three ROIs (100 pre-NMS peaks every time), not `False` as the design session
guessed. 123.tiff is a different, genuinely coincidental case: `n_self_hits=1` (the old filter
did catch a self-hit), but the newly-admitted replacement candidate happens to land within NMS
radius of another survivor, so the post-NMS count drops from 95 to 94 — exactly cancelling the
one self-hit the old filter used to remove, giving the same final `n_detections` (94) by a
different mechanism than 401/546.

### 3d. Corrections found during independent verification (2026-09-17)

- **A mechanism-ordering bug in the verification session's own first harness, found and
  fixed, not a design flaw.** Calling today's *real, unpatched* `find_and_suppress` with an
  already-blanked `img_channel` also blanks the template cut from it (both read the same
  array), producing an all-zero correlation map everywhere, not just near the seed —
  `deep_floor_median`/`mad`/`score_threshold_used` all collapsed to exactly 0, and every
  detection vanished. This is why simulating the design correctly requires a standalone
  function (`find_and_suppress_blanked`, and now the real patch in §4) that blanks the search
  channel *after* the template is already cut, matching §2's own "When" ordering exactly.
- **A seed-selection bug in the verification session's own first augmented-bank check, found
  and fixed.** All 14 canonical ROIs are a subset of the 49-ROI pins used elsewhere in this
  verification, so a naive "use the 49-ROI pin if this ROI is in it, else draw plainly" always
  took the pinned branch — never a plain `seed_index=0` draw. §3c's corrected table uses plain
  draws throughout, matching `harness.py`'s own methodology and the historical `probes/p3.log`
  numbers exactly.
- **The `tm_score`-moves-under-augmentation finding, mechanistically resolved.** On the
  mis-selected seed above, 245.tiff K=20 moved 2→1 under `scale_normalize=True`. Traced to the
  exact peaks involved: three peaks near (4878,2449)/(4867,2455)/(4904,2448) (full-ROI
  coordinates) are present in *both* arms' pre-NMS 100-peak pools, with the same winning
  template index in both — no new peak, no template-selection flip. Their scores differ only
  by ~0.003-0.007 (re-derivation noise: blanking a 2601 px square shifts the stride-8-sampled
  per-template `_robust_z` statistics computed over the whole 1024×1024 crop by a hair). Those
  three peaks form a suppression "chain" under greedy NMS (peak1 is within `nms_radius` of both
  peak2 and peak3, but peak2 and peak3 are *not* within `nms_radius` of each other) — a tiny
  score perturbation that changes which peak NMS visits first flips the survivor count from one
  to two. **This pathway cannot occur under production's actual configuration.**
  `scale_normalize=False` (production's only setting; `FSConfig`'s default, and
  `production.run_production_pipeline` never overrides it) means `fused_response` never calls
  `_robust_z` — raw `TM_CCOEFF` scores for any window that does not overlap the blanked square
  are then provably bit-identical between arms (matching-template correlation is a purely local
  sliding-window computation with no cross-window rescaling step in this branch), so there is
  no channel for blanking anywhere in the image to perturb an unrelated peak's score. The
  shared deep-floor *threshold* (a scalar cutoff, computed regardless of `scale_normalize`)
  does shift by a small amount (~5e-4 in this example) even when `scale_normalize=False`, but
  this only changes which peaks clear the bar, not the relative order of peaks safely above
  it — consistent with zero `tm_score` movement measured across every real,
  `scale_normalize=False` run this session (49 ROIs, seeds 1-4, the corrected augmented-bank
  rerun).

## 4. Code edits (the patch)

Patch: `../cleanup_harness/selfmask/image_blank_design/selfmask_imageblank_code.patch`.
**Apply it; do not retype it.** Verified `git apply --check` clean against the real repo at
`c7cf7ab`. Built and tested against a scratch copy first (§8); the real `midog_utils/` was
never touched while building this plan.

### 4.1 `midog_utils/template_match.py`: new function directly after `extract_peaks`

```python
def blank_seed_square(channel: np.ndarray, cx: float, cy: float, base_size: int):
    """
        Blank the seed's own refined-template footprint out of a COPY of a search channel,
        so its own match can never form as a correlation peak.

        channel (np.ndarray): the search channel for the whole ROI; never mutated.
        cx (float): template centre, x, unrounded.
        cy (float): template centre, y, unrounded.
        base_size (int): the odd template side length; the blanked square is exactly
            base_size x base_size, centred on ``(round(cx), round(cy))`` -- the same
            rounding convention ``read_padded_patch`` uses, so the blanked square lines up
            exactly with the footprint the template itself was cut from.

        Returns tuple[np.ndarray, int]: a blanked copy of ``channel``, and the exact pixel
        count blanked (clipped to the image; ``base_size**2`` unless the seed sits within
        half a template of the ROI edge, which does not occur for a readable seed since
        ``read_padded_patch``'s own border check already requires a half-``patch_size``
        margin, larger than half a template).

        The fill value is ``channel.min()`` over the WHOLE, unblanked ROI -- computed before
        this function writes anything, so a caller that reuses ``channel`` (e.g. for
        ``chromatin_od`` ranking) sees only the original, untouched array.
    """
    fill_value = float(channel.min())
    half = base_size // 2
    ix, iy = int(round(cx)), int(round(cy))
    h, w = channel.shape[:2]
    y0, y1 = max(0, iy - half), min(h - 1, iy + half)
    x0, x1 = max(0, ix - half), min(w - 1, ix + half)
    out = channel.copy()
    out[y0:y1 + 1, x0:x1 + 1] = fill_value
    n_blanked = (y1 - y0 + 1) * (x1 - x0 + 1)
    return out, n_blanked
```

Notes:
- Returns a copy; never mutates `channel` in place — `find_and_suppress` doesn't own
  `img_channel`, and `production.py`'s `chromatin_od` branch reuses the same `hem` reference
  for ranking after the search.
- `int(round(...))` is Python's round-half-to-even, the same convention
  `read_padded_patch`/`tightened_template_box` use elsewhere in this codebase; a different
  rounding would shift the blanked square one pixel off the template's real footprint on
  half-pixel ROIs (e.g. 246.tiff's `template_xy=(3627.5, 3615.5)`).

### 4.2 `midog_utils/find_and_suppress.py`

- **Module docstring:** "cuts a template, correlates it against the ROI, extracts peaks,
  suppresses overlaps, drops the seed's own self-hit," becomes "cuts a template, blanks its
  own footprint out of a copy of the search channel so it can't match itself, correlates
  against that copy, extracts peaks, suppresses overlaps,".
- **`FSConfig`:** delete `self_hit_radius: float = 5.0`. No field is added (unlike the
  superseded design, `n_blanked_px` is not a config field — it's computed from `cfg.base_size`,
  already present).
- **`find_and_suppress()`:** insert between `info["n_augmentations"] = len(templates)` and the
  `if cfg.border_pad:` branch:
  ```python
      # The template above is cut from the ORIGINAL, unblanked img_channel. Only now, after
      # that cut, do we blank a COPY for the search -- img_channel itself (and anything else
      # the caller derives from it, e.g. production.py's chromatin_od ranking) is untouched.
      search_channel, n_blanked = tm.blank_seed_square(img_channel, seed_x, seed_y, cfg.base_size)
      info["n_blanked_px"] = n_blanked
  ```
  and change both branches below (`if cfg.border_pad:` and its `else:`) to build `padded`/call
  `tm.fused_response` from `search_channel`, not `img_channel`.
- **Delete** `info["max_peak_score"] = ...`.
- **Delete** the whole `if len(centers): d_seed = ... else: ...` self-hit block, so that
  `info["n_detections"] = len(centers)` directly follows `info["n_after_nms"] = len(centers)`.

### 4.3 `midog_utils/production.py`

- Delete `SELF_HIT_RADIUS = 5.0`. No constant replaces it.
- Delete `self_hit_radius=SELF_HIT_RADIUS, ` from the one-line `FSConfig(...)` call, which
  stays on one line. No other change: `hem` (used both for the search and for `chromatin_od`
  ranking) needs no change, because `find_and_suppress` never mutates it.

### 4.4 `midog_utils/seed_selection.py`: `tightened_template_box` docstring only

"Caller owns: self-hit/seed-annulus removal must reference the returned centre, not
(cx, cy); ..." becomes "Caller owns: pass the returned centre, not (cx, cy), as
`find_and_suppress`'s seed_xy, since the search-channel blanking it does is centred there;
border readability at the returned centre must be re-checked; ground truth stays keyed to
(cx, cy), never the returned centre."

### 4.5 `midog_utils/FIND_AND_SUPPRESS_REFERENCE_DIFFS.md`

- **"No image masking"** heading and body are replaced with **"The seed's own footprint is
  blanked, not the whole reference-style mask"** — the reference's whole-annotation blackout
  and its zero-variance-window problem, contrasted with this design's narrower, well-defined
  (non-zero-variance) blanking, plus the cost this narrower approach accepts (the 301.tiff
  seed_index=1 case named explicitly) versus the reference's own risk of deleting a legitimate
  neighbouring detection.
- **"The self-hit is dropped with a tight radius, not the match radius"** is folded into the
  section above (the exact patch diff has the precise before/after text).

## 5. Doc edits (not in the patch; use the Edit tool)

Another session may be editing these docs, so **locate every edit by its quoted anchor
text, never by line number.** If an anchor isn't found verbatim exactly once, stop and
report. Replace `YYYY-MM-DD` with the date you apply this. All anchors below were read
directly from the live files while writing this plan (2026-09-17); re-check them if this plan
is applied significantly later.

### 5.1 `production_pipeline/PRODUCTION_PIPELINE_WALKTHROUGH.md`

**Edit 1.** After the paragraph ending "that reproducibility is itself a regression test.",
add a new paragraph:
```
**Updated YYYY-MM-DD:** the seed's own match is now excluded by blanking its refined-template footprint (`base_size x base_size`, centred on `template_xy`) out of a copy of the search channel before correlation runs (`template_match.blank_seed_square`), replacing a 5.0 px post-NMS filter. See `SELF_HIT_MASKING_PLAN.md` and `DECISIONS.md` D11. Numbers below were re-verified after that change.
```

**Edit 2.** In the verification-summary table (the one whose header row is
`| field | 245.tiff / tm_score | 403.tiff / tm_score | what it confirms |`), replace the four
rows starting `` | `n_after_nms` ``, `` | `n_self_hits` ``, `` | `seed_self_score` `` and
`` | `n_detections` `` with:
```
| `n_blanked_px` | 2209 | 2601 | pixels of the refined template's own footprint (47^2 / 51^2) excluded before correlation |
| `n_after_nms` | 87 | 99 | NMS removed 13 / 1 of the 100 pre-NMS peaks |
| `n_detections` | 87 | 99 | final ranked-list size; nothing is removed after NMS |
```
(Numbers are illustrative from a fresh run at apply time — re-derive them, don't copy these.)

**Edit 3.** Replace the whole paragraph beginning "The `seed_self_score == max_peak_score`
equality matters because" through "That did not happen on either test ROI here." with:
```
No detection can be extracted from inside the blanked square, because its source pixels are gone before correlation runs. The filter this replaced removed only points within 5.0 px, so it missed the seed's own match whenever that match fell below the `max_peaks` cutoff or moved under an augmented template bank (`FIND_AND_SUPPRESS_REFERENCE_DIFFS.md`).
```

**Edit 4.** In the `Seed` dataclass block, change `self-hit removal and NMS-adjacent radii
reference template_xy` to `the search-channel blanking and NMS-adjacent radii reference
template_xy`.

**Edit 5.** In the module-constants table, delete the row starting `` | `SELF_HIT_RADIUS` ``.

**Edit 6.** In the call diagram's `FSConfig(` node, change the line
`│           self_hit_radius=5.0, deep_floor_z=-1.5, border_pad=True,` to
`│           deep_floor_z=-1.5, border_pad=True,`.

**Edit 7.** In the call diagram, directly after the node
`│  ├─ tm.build_augmentations(patch, base_size, scales=(1.0,), n_angles=1, flips=(False,))`
(its 5 text lines), insert:
```
│  ├─ tm.blank_seed_square(hem, tx, ty, base_size)              template_match.py
│  │    Blanks a base_size x base_size copy of hem at the template's own footprint --
│  │    never the original hem, which production.py reuses unblanked for chromatin_od
│  │    ranking. Verified: n_blanked_px 2209 on 245.tiff, 2601 on 403.tiff.
│  │
```

**Edit 8.** Delete the node
`│  ├─ self-hit removal: keep points with distance(point, template_xy) > 5.0 px`, together with
its 3 text lines and its trailing `│  │` spacer. Then change
`89 rows on 245.tiff, 98 on 403.tiff (100 - NMS losses - the 1 self-hit)` to
`87 rows on 245.tiff, 99 on 403.tiff (100 - NMS losses)`. (Re-derive the exact counts at apply
time — these illustrate the shape of the edit, not the exact numbers.)

**Edit 9.** Under "Things to know, not bugs to fix", replace item 1's text, from
"**`SELF_HIT_RADIUS = 5.0` px is correct only paired with `n_angles = 1`.**" through
"re-enter the ranked list.", with the text below. Keep the `1. ` prefix and the numbering.
```
**Resolved YYYY-MM-DD: augmented template banks no longer leak the seed's own match.** The 5.0 px post-NMS filter used to miss self-matches that moved under rotation/flip/scale (12.66 px on 013.tiff, 5.39 px on 233.tiff, with the harness's 8-template bank); blanking the template's own footprint before correlation covers both, verified. A narrower gap remains where a *different* nearby detection (not the seed's own match) sits within one match radius of the seed but outside the blanked square -- see `SELF_HIT_MASKING_PLAN.md` sec 2, "Known cost, accepted".
```

**Edit 10.** In the paragraph starting "**What `find_and_suppress` deliberately does NOT
do**", change
`no image masking of existing annotations before correlating (the seed's own detection is
dropped afterward instead,` + `same effect, no zero-variance-window artefact);` to
`no masking of every existing annotation before correlating (only the seed's own template footprint is blanked — narrower and well-defined: TM_CCOEFF against a flat window computes to ~0, not an extreme value, so no zero-variance-window artefact);`.
Every other item in that paragraph is a negation; an affirmative item there reads, against
the paragraph's own heading, as denying that the footprint is blanked. Corrected at source
2026-09-17 — see the revision history.

**After the edits, this must print nothing:**
`grep -n "self_hit\|SELF_HIT_RADIUS\|n_self_hits\|seed_self_score\|max_peak_score\|self-hit\|dropped afterward" production_pipeline/PRODUCTION_PIPELINE_WALKTHROUGH.md`
(`SELF_HIT_RADIUS`, not `SELF_HIT`: Edit 1's text names `SELF_HIT_MASKING_PLAN.md`.)

### 5.2 `PIPELINE_BREAKDOWN.md` (historical record: one marker only)

Replace the line ``**Self-hit removal** (notebook `suppress()`):`` with:
```
**Self-hit removal** (notebook `suppress()`; replaced in production on YYYY-MM-DD by blanking the seed's own template footprint before correlation, see `SELF_HIT_MASKING_PLAN.md` and `DECISIONS.md` D11):
```

### 5.3 Do not touch

- **`midog_utils_full/`** (archive).
- **Historical records:** `DECISIONS_UNVERIFIED.md`, `D8_TEMPLATE_ANCHOR.md`,
  `INVARIANTS_HISTORY.md`, `PRODUCTION_PIPELINE_CLEANUP*.md`, `Research Logs/`, and every
  existing `DECISIONS.md` entry's body (D10 and D11 included — only their own dated-amendment
  mechanism may add to them later).
  - D8's "self-hit or seed-annulus removal must reference the returned centre" still holds;
    the blanked square is centred on `template_xy`.
- **`production_pipeline/production_pipeline_overview.pptx`.** Already stale since `c066829`
  (slides 11, 12 and 17 mention `self_hit_radius=5.0`). Leave it, mention it in the report.
- **Completed audit scripts.** `bbox_refinement_three_way_chromatin_od_audit.py` (~line 300)
  and `hembbox_precision_at_k_audit_round3.py` (~line 219) read production source as text and
  are already known to report stale values. Don't re-run or edit them.
- **`production_hematoxylin_only/bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb`**
  and its prompt, `BBOX3WAY_49ROI_3SEED_PROMPT.md` — see §5.4.
- **`../cleanup_harness/harness.py`.** Reuse it, don't edit it.
- **`../cleanup_harness/selfmask/`'s own files** (`selfmask_code.patch`, its reference
  captures, `superseded_5um/`, `probes/`, `SELF_HIT_MASKING_PLAN.valid_mask_disc_design.md`,
  `SELF_HIT_MASKING_PLAN.midog_utils_full_draft.md`) — the superseded design's own assets,
  untouched and still valid for that design if it's ever revisited.

### 5.4 Downstream consumers that break (do not edit; report them)

`production_hematoxylin_only/bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb`
calls `prod.run_production_pipeline`. If re-run after this change it breaks in three places:
- it reads `prod.SELF_HIT_RADIUS` in cell 1 (`AttributeError`);
- it reads `RESULTS['n_self_hits']` in cell 6 (`KeyError`);
- its Gate B reasons about the 5 px self-hit disc.

Its prompt, `BBOX3WAY_49ROI_3SEED_PROMPT.md`, pins "5 px self-hit removal". Its saved outputs
are a pre-change reference. Re-running it needs a port, which is the user's call — same
disposition D10 gave this notebook.

---

## 6. Execution and verification

Work from the repo root `/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked`. Define these
at the start of **each** shell command, because shell state does not persist between calls:
`PY=/Users/mohinianand/anaconda3/bin/python3; HD=/Users/mohinianand/Desktop/AnnotateDx/cleanup_harness; ID=$HD/selfmask/image_blank_design`.

The bare `python3` can't import cv2. Keep a log of every command and its final output lines
somewhere under `$ID/`. **Do not commit.**

**Reference captures on disk (never write into these folders):**
- `$HD/runs/imageblank_ref_pre/`: pre-edit harness capture (unmodified repo) plus
  `recent_conditions.pkl` and `recent_seeds.pkl`.
- `$HD/runs/imageblank_ref_post/`: the same, from a scratch copy with this patch applied.
- Provenance: §8, below.
- `harness.py extref` writes `extref_merge.csv` into the run folder it's given. Only run it on
  your own labels.

**Retries:** `harness.py capture` and `imageblank_check.py recent-capture` refuse to overwrite
an existing capture. On a retry use a fresh label (e.g. `post_imageblank2`). Never delete a run
folder.

**Commands that exit non-zero by design** (judge these by their printed output):
- `pyflakes` (Step 3.1, if it reports only the known `__init__.py` unused-import lines);
- `harness.py compare pre_imageblank post_imageblank` (a raw diff listing — this is expected
  to differ; the real gate is `imageblank_check.py check`, not this command's exit code);
- `harness.py extref post_imageblank`;
- `imageblank_check.py reference post_imageblank`.

**Runtimes on this machine** (this verification session's own measurements, idle machine):
- `harness.py capture`: ~2.3 min.
- `recent-capture`: ~9-40 min per `--part`, depending heavily on system load — this session ran
  all four (`pre`/`post` × `conditions`/`seeds`) in parallel with each other under otherwise-
  idle conditions and still saw ~25-40 min each; expect longer under load, as the archived
  design's own provenance also found.

### Step 0: preflight (stop on any failure)

1. `pgrep -fl "[n]bconvert|[i]pykernel|[h]arness\.py|[s]elfmask_check\.py|[i]mageblank_check\.py"`
   must print nothing. (The brackets stop the pattern from matching its own shell.) If a job is
   running, wait for it or ask the user. Don't patch code under a running experiment.
2. `git status --short midog_utils production_pipeline PIPELINE_BREAKDOWN.md` must print
   nothing; otherwise stop and ask the user. Also, `grep -c "^## D11 " DECISIONS.md` must
   print `1`.
3. `shasum -a 256 -c $ID/expected_pre_edit_sha256.txt` must print 5 lines ending in `: OK`. A
   mismatch means the code moved since this plan was verified; stop.
4. `git apply --check $ID/selfmask_imageblank_code.patch` must print nothing and exit 0.

### Step 1: pre-edit capture (proves this machine reproduces the reference)

1. `$PY $HD/harness.py capture pre_imageblank`
2. `$PY $HD/harness.py compare imageblank_ref_pre pre_imageblank` must end with
   `RESULT: IDENTICAL (apart from reported exceptions)` and report
   `api.json differences (0, reported, not failing)`. If not, stop.

### Step 2: apply the edit

1. `git apply $ID/selfmask_imageblank_code.patch`
2. `shasum -a 256 -c $ID/expected_post_edit_sha256.txt` must print 5 `: OK` lines.
3. Read the five edited files and confirm they match §4.
4. This must print `ok`:
   `$PY -c "import sys, dataclasses; sys.path.insert(0, '.'); import midog_utils.template_match as tm; import midog_utils.production as prod; from midog_utils.find_and_suppress import FSConfig; assert hasattr(tm, 'blank_seed_square'); names = [f.name for f in dataclasses.fields(FSConfig)]; assert 'self_hit_radius' not in names; assert not hasattr(prod, 'SELF_HIT_RADIUS'); print('ok')"`

### Step 3: static checks

1. `$PY -m pyflakes midog_utils/*.py production_pipeline/run_pipeline.py` must print exactly
   the same known `__init__.py:14:1: '.<module>' imported but unused` lines D10's plan
   recorded (8 lines), and nothing else new.
2. `grep -rn "self_hit\|SELF_HIT\|n_self_hits\|seed_self_score\|max_peak_score" midog_utils production_pipeline --include=*.py --include=*.ipynb`
   must print nothing.

### Step 4: post-edit captures and gates

1. `$PY $HD/harness.py capture post_imageblank`
2. `$PY $ID/imageblank_check.py recent-capture post_imageblank --part conditions` and
   `$PY $ID/imageblank_check.py recent-capture post_imageblank --part seeds`. They may run in
   parallel. Each command's first output line must be
   `importing midog_utils from /Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked/midog_utils`.
3. **Gate A (bit-exact where nothing changed):**
   `$PY $HD/harness.py compare imageblank_ref_post post_imageblank` must end with
   `RESULT: IDENTICAL (apart from reported exceptions)` and report 0 api.json differences —
   this captures your own machine against this plan's own post-edit reference, both built from
   a patched scratch copy, so unlike `pre_imageblank`'s compare above, this one really is
   bit-exact.
4. **Gate B (attribution):**
   `$PY $ID/imageblank_check.py check imageblank_ref_pre post_imageblank 2>/dev/null > $ID/check_output.txt`
   must end with `CHECK: PASS` and `failing checks: 0`. It reports the same TP@K changes as
   §3b/§3c and a "reported (not failing)" table of near-seed leaks (289/548/301-seed1-style
   rows) — read it, don't just check the exit code.
5. **Gate C (numbers):**
   - `$PY $ID/imageblank_check.py reference post_imageblank` ends with
     `REFERENCE: DIFFERS (expected after the edit; see check)`.
   - `$PY $HD/harness.py extref post_imageblank` shows `n_detections exact: 0/84` and
     `EXTREF: FAIL` (by design — `n_detections` shifts by one on almost every run), with all
     84 precision/recall rows within `0.0006`.
   - `$PY $HD/harness.py sanity post_imageblank` must print `SANITY: PASS`.

### Step 5: demo notebook and CLI

1. Run:
   ```
   /Users/mohinianand/anaconda3/bin/jupyter nbconvert --to notebook --execute production_pipeline/production_pipeline_demo.ipynb --output-dir $HD/runs/post_imageblank --ExecutePreprocessor.timeout=3600
   ```
   Confirm it runs with zero errors.
2. Run:
   ```
   $PY production_pipeline/run_pipeline.py 403.tiff --rank-key chromatin_od --save-fig $HD/runs/post_imageblank/cli_403.png
   ```
   Confirm it runs and the seed/summary line looks sane (base_size, `n_detections`).
3. If both pass, refresh the tracked notebook's saved outputs: re-run the same `nbconvert`
   command with `--inplace`. Then `git status --short production_pipeline` should list only
   `M production_pipeline/production_pipeline_demo.ipynb`.

### Step 6: docs and decision record

Make the §5 edits. Then re-run Step 3.2 and the §5.1 grep.

Last, and only if Gates A, B and C and Step 5 all passed, append this amendment at the end of
the D11 entry in `DECISIONS.md`, directly before the `---` line that precedes
`## Cross-references`. Fill in the date and the numbers from your own outputs:
```
### Amendment, YYYY-MM-DD — applied

Applied to `midog_utils/` via `SELF_HIT_MASKING_PLAN.md`. Harness and recent-experiment captures bit-identical to `imageblank_ref_post`; `imageblank_check.py check`: PASS, N pairs, M TP@K changes, all attributed to the freed pool slot; `reference`: rows changed as pre-registered.
```
If any gate failed, don't add the amendment; report instead.

### Step 7: report

Report to the user:
- each gate's final output lines;
- the §3b/§3c change tables as `check` printed them;
- the §5.4 downstream breakage and the stale `.pptx`;
- whether the D11 amendment was appended;
- `git status --short`.

Do not commit unless asked.

## 7. If something doesn't match

- **Step 1 compare fails.** The environment or code isn't what `imageblank_ref_pre` was built
  from — check the interpreter path and `git log -3 -- midog_utils`. Don't proceed.
- **Step 2.3's `ok` check fails.** The patch didn't apply cleanly, or applied against the
  wrong baseline. Re-check `git apply --check` and Step 0.3.
- **Gate A fails.** Your patched code differs from what built `imageblank_ref_post` — most
  likely `blank_seed_square` runs before, not after, `build_augmentations` (§2's "When" is
  load-bearing — see §3d's own first-attempt failure for exactly this bug), or the rounding
  convention differs from `read_padded_patch`'s.
  - If `n_blanked_px` doesn't equal `base_size**2` exactly, the blanking is clipped at a ROI
    edge (shouldn't happen for a readable seed) or the square isn't centred correctly.
  - If only `deep_floor_*`/`score_threshold_used` differ by more than a rounding amount, the
    blanking or the threshold computation runs in the wrong order relative to each other.
- **Gate B: a detection inside the blanked square.** The square isn't centred on `seed_xy`
  (that is `template_xy`, not `click_xy`), or its size isn't `cfg.base_size`.
- **Gate B: unexplained lost or added rows, or changed shared rows.** The template was cut
  from the blanked channel instead of the original (§3d's first-attempt bug), or
  `build_augmentations`/`fused_response` received something other than the blanked
  `search_channel`.
- **Gate B: many "reported (not failing)" near-seed leaks appear where none were expected.**
  Not itself a failure — this design accepts them by design (§2) — but compare the count
  against §3b/§3c/§3d's tables; a large increase suggests `base_size` values shifted upstream
  (a `seed_selection.py` change), not a bug in this patch.
- **Outside `midog_utils/`: an `FSConfig(self_hit_radius=...)` TypeError, or a missing
  `SELF_HIT_RADIUS` or `n_self_hits`.** Expected for the §5.4 notebook. Every other such file
  already fails since `c066829`, except the two audit scripts in §5.3.

## 8. How the references and numbers were produced (for audit)

- **§3a/§3b (design session, 2026-09-17):** `blank_seed_variant_3roi.py` /
  `blank_seed_variant_49roi.py` and their logs/CSVs, archived unchanged at
  `../cleanup_harness/selfmask/image_blank_design/`.
- **§3c/§3d (independent verification session, 2026-09-17, same day):**
  - Mechanism claims: a synthetic realistic-magnitude flat-window `TM_CCOEFF` test; a
    real-data band-profile and pre-NMS-peak-distance check on 403.tiff.
  - §3b's numbers: `find_and_suppress_blanked`, a standalone function modelled directly on
    `find_and_suppress.py`'s own body (not on the design session's `blank_seed_variant_*.py`
    scripts), calling the same real, unmodified `template_match`/`nms` primitives in the
    order §2 specifies. Cross-checked against the design session's own 49-ROI CSV: every
    pooled number, every moved (ROI, K) cell, every near-seed leak's distance/rank/score
    reproduced exactly.
  - `tm_score` at scale: the same prototype, `rank_key="tm_score"`, all 49 ROIs.
  - The augmented bank: the same prototype with `scales=(0.8,1.2)`, `n_angles=2`,
    `flips=(False,True)`, on `harness.py`'s own 1024 px crop (`harness.py`'s
    `capture_aug`, lines ~131-160, read directly to match its cropping exactly) — first
    attempt used mis-selected (49-ROI-pinned rather than plain-draw) seeds for all six ROIs
    tested, caught by comparing against `probes/p3.log`'s stored numbers and corrected.
  - Seeds 1-4: the same prototype at `seed_index` 1-4 on the 14 `images/extra_valid` ROIs,
    cross-checked against the *superseded* design's own pre-registered seed-robustness table
    (an independent design, from an independent session) — 4 of 4 moved cells matched exactly.
  - The `scale_normalize=True`/small-crop `tm_score` finding: isolated to the specific peaks
    involved by re-deriving the full pre-NMS 100-peak pool for both arms directly (not through
    `find_and_suppress`), confirming the implicated peaks are present in both pools with
    unchanged winning-template indices, then confirming determinism (each arm reproduces
    itself exactly on repeated runs) before accepting the finding as real rather than a race
    or floating-point non-determinism artefact.
- **Building `imageblank_ref_pre`/`imageblank_ref_post`:**
  - `midog_utils/` and `production_pipeline/` were copied to a scratch directory, with
    `images/` and `databases/` symlinked. The edit (§4) was applied there, and the patch
    (`selfmask_imageblank_code.patch`) generated from a `diff -u` against the real repo,
    confirmed clean with `git apply --check`.
  - `harness_imageblank_scratch.py` (a copy of `harness.py` with only its `REPO` line changed
    to point at the scratch copy — created fresh for this run, not committed) captured
    `imageblank_ref_post`; the unmodified `harness.py` captured `imageblank_ref_pre` directly
    against the real, unpatched repo.
  - Both captures completed with zero recorded failures; `sanity` on the post capture: PASS;
    `extref` on the post capture: `n_detections exact: 0/84`, all 84 rows within `0.0006`,
    `EXTREF: FAIL` (by design, matching §3a's expected pattern).
  - The real patched `find_and_suppress` (imported from the scratch copy in an isolated
    subprocess, to avoid Python's module-cache picking up the already-imported real package)
    was cross-checked against the `find_and_suppress_blanked` prototype on 403.tiff: identical
    `n_detections` (99), identical `n_blanked_px` (2601), matching the prototype's own 49-ROI
    sweep result for the same ROI exactly.
- **`imageblank_check.py`** reuses `recent_capture`/`recent_compare`/`reference`/`nb_text`
  unchanged from the superseded design's `selfmask_check.py` (design-agnostic: they run live
  production code or diff pickles without knowing which self-hit design is active); only
  `attribute_pair`/`check` are design-specific — see its own module docstring for exactly what
  differs (the `n_blanked_px` contract, a square rather than disc containment test, and no
  fixed "within 2 radii" geometric bound on gained rows, since this design's excluded region
  size varies with `base_size` rather than always equalling `nms_radius`).
- **Negative tests of `check`.** Run against a tampered copy of one real run pair's info/
  detections (403.tiff, `chromatin_od`), each producing exactly one failure, with the
  unmodified pair producing none:
  - a row injected 20 px from the seed (inside the blanked square's half-width for a 51 px
    template);
  - a row injected 40 px from the seed (outside one match radius — should NOT fail; used to
    confirm the "report, don't fail" path for the accepted near-seed-outside-square cost);
  - a shared row's score changed by 1e-3;
  - a far row dropped with no replacement;
  - a high-scoring far row injected with no removed row to explain it;
  - two shared rows swapped in rank order;
  - `n_blanked_px` set to a value that isn't `base_size**2`;
  - an extra info key added beyond the `n_blanked_px` contract;
  - one of the three removed info keys (`max_peak_score`/`n_self_hits`/`seed_self_score`)
    left present in the "post" side.

  Full results, including exact commands and output diffs, are in
  `imageblank_negative_tests_output.txt` in this folder.
