# Edit plan: correcting the three `tm_threshold_axis_sweep` notebooks

Companion to `Research Logs/2026-09-03-tm-axis-sweep-audit.md` (12 findings, 3 tiers).
Reviewed by an independent subagent; 19 issues raised, all resolved (see "Review outcomes").
**Nothing here has been executed.** Awaiting go-ahead.

## Principles

* **P1 — Correct in place, keep the trail.** Wrong conclusions are rewritten *and* carry a
  Correction note saying what they used to say. Precedent: commit `bcacb55`.
* **P2 — Keep the `od` fix notebook-local.** Not because `midog_utils` is in flight
  (`chromatin.py` is tracked and clean — the original rationale was wrong), but because
  changing `score_detections`' silent-NaN-at-borders behaviour moves output for every other
  consumer of that function. That is a library change needing its own verification pass →
  follow-up F3.
* **P3 — Fix the artifact, not just the prose.** A CSV with a known-bad column is a landmine.
* **P4 — Corrections, not new experiments.** Where a question can't be answered by what was
  run, state the limitation; don't launch a bigger run.
* **P5 — Measured, not asserted.** These notebooks hold themselves to "check, don't assume".

## Determinism, ordering, baseline

Deterministic throughout (fixed `default_rng` streams, `kind="stable"`/`"mergesort"` sorts,
global `lexsort` in `extract_peaks`, order-insensitive NMS suppression). Residual risk is
environment drift since the original runs, so V1/V2 are **diff reports**, not asserts.

Order is forced by file reads: **v1 → v2 → lcc** (~220 s + 150 s + 370 s plus ROI loads).

### E0 · Baseline — must happen before any re-run ⚠ needs your decision
All three notebooks **and** all three CSVs are **untracked** (nothing is gitignored — they
were simply never `git add`ed). So "re-run and diff against the current file" has no baseline:
the persist cells overwrite unconditionally.

* **E0a (recommended, needs your OK).** `git add` exactly these six paths and commit as a
  pre-correction baseline. This also closes a pre-existing gap against the repo's own
  convention. **Never `git add -A` here** — the tree holds 5 modified `midog_utils` modules
  and ~40 untracked files from the unrelated bbox-threshold experiment; sweeping those in
  would entangle two experiments and make the correction diff unreadable.
* **E0b (unconditional).** Also copy the six files to the scratchpad. One `cp`; decouples
  the first re-run from the commit having landed.

---

# A. `tm_threshold_axis_sweep.ipynb` (v1) — 23 cells
No computational change: its CSV must come back byte-identical (gate V1).

| # | cell | change | finding |
|---|---|---|---|
| **A1** | 7 (code) | Bucket `pool` once and record `n_dup_fp` — detections within `radius` of some GT that bucket `non_human_findings` — into a **new side list**, *not* into `ctx`. (`compare.py:195 base = dict(context)` turns every `ctx` key into a CSV column; that is why v1/v2/lcc CSVs are 37/39/41 wide. A `ctx` entry would break V1 for an unrelated reason.) | 8 |
| **A2** | new, after 21 | Print `coverage_frac` pivoted by (domain, z) from `results`, plus A1's table. | 1 |
| **A3** | 12 (md) | Caveat the composition figure with the *correct* ranges: across the nine z it plots, `n_detections` spans **3,070–17,397** and `coverage_frac` **0.256–0.951**; `precision_full` spans **0.109%–4.87%** and is dominated by list length, not list quality. | 1 |
| **A4** | new, before 19 | Assert the real theorem: **for `arm=='tm_score'`**, `recall_at_budget` has zero spread across z in every (domain, budget) group where `n_detections >= budget` — verified 56 groups, 0 violations, max spread 0.0. **The arm filter is the content, not scoping**: the same statement on `chromatin_od` fails in **33 of 56** groups (spread up to 0.0968), because filtering by score removes entries from arbitrary positions in an od-ranked list. Use `np.ptp(...) < 1e-12`. Print the boundary case: only budget 5000, only z=3.0, 6 of 7 domains (breast cancer keeps all nine at n=5,711). | 10 |
| **A5** | 18 (md) | Say the `tm_score` curve is flat by construction at this budget — smallest pool anywhere is **3,070** (canine lung cancer, z=3.0), **31×** K=100 — so the only live comparison in the panel is `chromatin_od`. | 10 |
| **A6** | 21 (code) | Add `n_read95_finite` / `n_read95_unreachable`. (Additive, not `cp.summarise`: its default `by` includes `file_name`, so it would give n=1 per group, and its docstring is about aggregating over *seeds*, not domains.) Cover every `read_*` the cell might display — `read_99` and `read_100` contributors fall 6→1 across z. | 4 |
| **A7** | 20 (md) | Label the survivorship: contributors are 7/7/7/7/7/7/7/**5**/**4**, so `median_read_95` at z ≥ 2.5 is over an easier domain set. Drop "the single table to defend a `z` choice from". | 4 |
| **A8** | 22 (md) | Replace the deferral with the result: **at every budget where the pool never truncates (K ≤ 2000 here)** the floor does not change `recall_at_budget` at all on `tm_score`, and by ≤ 0.0096 on `chromatin_od` (4 of 7 domains; 0 in the rest). What the floor changes is candidate volume and ROI coverage. Add coverage + single-seed caveats. | 10, 1, 5 |
| **A9** | 2 (md) | Declare the limitation: one ROI/domain, `SEED_INDEX=0` only, against `compare.summarise`'s "Never report a single-seed number (Step 0 rule 5)". Cite `tm_ccoeff_headtohead_seed_variance.csv`: `read_95` seed SD 260–2,744 on these ROIs, and seed 0 sits on the **max** for 301.tiff, the **min** for 246 and 548. | 5 |
| **A10** | 0 (md) | Link the audit log. | — |
| **A11** | 3 (code) | Fix the ambiguous comment `# production median + 0.5*MAD floor`. **This lives in v1, not v2.** Point it at `tm_variant_sweep.FLOOR_Z = 0.5` / `premise_test.py`'s `med + 0.5*mad`, and note it is *not* `FSConfig.score_threshold = 0.5`, a raw `TM_CCOEFF_NORMED` score that coincidentally shares the number. Add a forward pointer to the later z=1.0 operating point. | new |

---

# B. `tm_threshold_axis_sweep_v2.ipynb` — 27 cells
**The only notebook whose CSV changes** (its `chromatin_od` rows only).

| # | cell | change | finding |
|---|---|---|---|
| **B1** ⚠ | 14 (code) | Fix the `od` NaN bug Fix 1 introduced: mirror lcc cell 19 — `OD_PAD = tm.BASE_SIZE // 2` (25, exactly sufficient since `score_detections` defaults `window=tm.BASE_SIZE`), `BORDER_REPLICATE`-pad `hem`, score at shifted coords, keep the naive `od` for the record, `assert n_nan_fixed == 0`. **Acceptance test:** `201.tiff` is the same seed and a byte-identical 51 px template as lcc, which already agrees on `n_detections` (19,873 at z=0.5; 16,529 at z=1.0) — so v2's `chromatin_od` `read_95` must land **exactly** on lcc's **814** (z=0.5) and **801** (z=1.0), up from 19,373 / 16,143. That catches a wrong shift sign or an off-by-one to the unit. | 3 |
| **B2** | 14 (code) | A1's `n_dup_fp`, side list not `ctx`. | 8 |
| **B3** | new, before 14 | Explain why `od` needs its own padding: the match's reachability margin (`PAD`, template-size-dependent) is not the chromatin window's margin (fixed 51 px), and `chromatin.rerank` already documents "nan (border) sorts last". | 3 |
| **B4** ⚠ | 22 (md) | Rewrite against the corrected table + Correction note. **After B1 this becomes self-contained**: the same cell's previous version said +18,720 and it now says ≈ **+161** — a v2-internal before/after, no cross-notebook footwork needed. Label it **"the cost of both fixes"** (padding *and* radius), not "radius cost" — v1→v2-fixed spans both. | 3 |
| **B5** | 0 (md) | Fix "stranded 2 of 217": 2 mitoses sit in the unreachable margin but only **1** (ann 14779) was missed as a result; ann 15022 was still claimed inside the 29.6 px radius. Accounting is **1 border + 5 NMS = 6**, matching 211/217 and 0.0276 × 217. | 6 |
| **B6** | 10 (code) | Measure it: cell 10 already rebuilds the old-radius NMS state — bucket it and print the missed `ann_id`s with `valid` at each click. Expect 6 missed, 1 with `valid=False`. | 6 + P5 |
| **B7** | 2 (md) | Disclose the full cost of overriding `check_nms_radius`: quote its documented defect, then give B2's measured consequence (111 → 281 duplicate-FPs on 301.tiff, ~1.1% of the FP pool). The decoupling is **not** reverted; only the cost statement was incomplete. Also note z=0.5 was the then-current floor and the operating point later moved to 1.0 (value unchanged — moving it would silently shift cell 23's snapshot). | 8 |
| **B8** | 17–19 | Coverage grid beside the recall grid; rewrite cell 19's verdict: full recall at every z ≤ 1.0 **on lists covering 88.1–99.8% of the ROI**. Credit the two fixes with the two named mechanisms (ann 14779; the 5 NMS anns), not with the aggregate recall number. | 1 |
| **B10** | 26 (md) | Drop `assert_floor_not_limiting` from the "verified" list — with `DEEP_FLOOR_Z = min(Z_LEVELS) − 0.5` it cannot fire; it is drift insurance, not evidence. Restate recall per B8; add A9's caveat. | 9, 1, 5 |
| **B11** | 0 (md) | Link the audit log. | — |
| **B12** | 17 (md) → new check | Cell 17 asserts the `tm_score`/`chromatin_od` full-list equality is "verified in the first notebook". v2 has no such check of its own — add v1 cell 11's three-liner (holds: 0/63). Same defect C9 fixes in lcc, one notebook earlier. | 11 |

*(B9 from the draft is withdrawn: v2 cell 3 has no comment. Its substance moved to A11.)*

---

# C. `tm_threshold_axis_sweep_largest_cc.ipynb` — 38 cells

| # | cell | change | finding |
|---|---|---|---|
| **C1** | 0 (md) **and 3 (comment)** | The "stale `CURRENT_Z`" wording appears twice. It is an **error of fact**, not just an anachronism: lcc says "v2's own stale `CURRENT_Z=0.5` *config-cell comment*", but v2's config cell has no comment — lcc is describing **v1's** comment and attributing it to v2. Also, v2 predates the z=1.0 decision (`tp_fp_score_distribution.ipynb` calls it "just decided on"). | new |
| **C2** | 15 (code) | Correct the coincidence denominator and **report both, labelled by draw rule**: cell 13's population uses v2's plain uniform rule → `n_coincide/n_pool` = **0.4015**; the main sweep's rule is `draw_seed_with_retry` → `n_coincide/largest_cc_ok` = **0.4165** (exact, not approximate: uniform-without-replacement-until-first-success is uniform over the passing subset). They agree to 0.015, and the populations coincide *measurably* — cite cell 19's own printed `seed matches v2's: 7/7`, `retries: 0/7`. Add the independence clause the existing `np.prod` already assumed. Conclusion ("7/7 is unsurprising") survives; the published 0.8699 was 2.2× optimistic. | 7 |
| **C3** | 19 (code) | A1's `n_dup_fp` (`det_out` already exists here), side list not `ctx`. | 8 |
| **C4** | new after 25; 26 | Coverage grid; rewrite the verdict per B8. lcc coverage at z ≤ 1.0 is 0.8957–0.9994. | 1 |
| **C5** | 27 (code) + 37 bullet 1 | Add `n_detections_ratio` and `coverage_delta` beside the recall delta. Rewrite bullet 1: the positive recall deltas at z ∈ {1.5…3.0} **are** the volume increase seen through a coverage-saturated metric, not a second, independent benefit — evidenced by `201.tiff`, where the template didn't change and all three deltas are exactly 0. (Soft-tissue sarcoma z=3.0: ×2.06, +0.229, +0.046; mast cell z=3.0: ×1.69, +0.170, +0.152.) | 2 |
| **C6** ⚠ | 28–29 | Depends on B1. Once v2's `chromatin_od` is clean, the "od-fix confound" caveat and the "`201.tiff` proves it directly … a swing of −15342" passage no longer describe reality — replace. Re-run cell 29 against the corrected CSV. Keep `201.tiff` as the null control, now stating the simpler fact: same seed, byte-identical template, all 13 shared columns identical **on both axes**. This *simplifies* deliverable (d) into a clean tightening-only comparison. | 3 |
| **C7** | 34 (md) | Strike **only** the contrastive clause "only deliverable (d)'s *read-depth* numbers carry the od-fix confound" — B1 removes that confound. **Keep the sentence before it** ("both its `chromatin_od` numbers and this arm's use the fixed (padded) od…"), which stays true: `v2_reference` came from `tp_fp_chromatin_*`'s already-fixed od, not from the v2 CSV. Then add the second contaminant: duplicates of TPs surviving the decoupled radius are high-scoring by construction, land inside the TP range, and do so asymmetrically between arms of different pool size. | 8 |
| **C8** | 37 bullet 3 | Replace "mixed, no consistent direction… no domain property predicts the sign" with: every `read_95_delta` (+1416 … −1547) is smaller than its own ROI's seed SD (260–2,744), and seed 0 sits on a distribution endpoint for 3 of 7 ROIs — the signs are noise. | 5 |
| **C9** | 24 (md) → new check | Re-measure axis independence rather than asserting it "since the same candidates just get re-sorted" (holds: 0/63). | 11 |
| **C10** | 2 (md) | A9's single-seed limitation. | 5 |
| **C11** | 37 "Net read" | Fold in C5/C8; drop the two-part "does not cost recall **and** does increase candidate volume" framing. | 2, 5 |
| **C12** | 0 (md) | Link the audit log. | — |
| **C13** | 37 bullet 2 | **Mechanism sentence only.** "A smaller template has a broader, less-selective correlation response (more of the map clears the same per-map z cut)" is contradicted at the deep floor in 2 of 6 tightened domains: `301.tiff` 25,449 → 25,447 (×1.000) despite 51→41 px, `246.tiff` ×1.017 — against `459.tiff` ×1.389 and `548.tiff` ×1.405. At low z the count is bounded by peak-spacing/NMS packing geometry, not the threshold. **Do not touch the bullet's numbers** — "1.02–1.37 in the 6 domains that actually tightened" is correct (measured 1.0163–1.3730). | 12 |

---

# E. Also
* **E1 — done already.** Three errors in the audit log itself, found by review: Finding 10's
  "3,424 / 34×" → **3,070 (canine lung cancer) / 31×**; the unqualified "at every budget"
  identity → budget-scoped; Finding 3's NaN range → 61–746 per (domain, z) cell.

# Not doing, and why

| Not doing | Why |
|---|---|
| `od` fix into `midog_utils/chromatin.py` | P2 as **restated** — changes border behaviour for every consumer; needs its own verification pass (→ F3). *(The draft's "midog_utils is in flight" reason was wrong: `chromatin.py` is tracked and clean, and these notebooks already depend on uncommitted `channels.py`/`compare.py`/`template_match.py`/`nms.py` changes.)* |
| Reverting v2's NMS decoupling | The argument is sound; de-dup and match radii answer different questions. Only the cost disclosure was incomplete (B7). |
| Re-running at 5 seeds | P4 — a new experiment (~1.5 h), not a correction → F1. |
| Re-running the two `tp_fp_*` notebooks | Verified: cell 9 in each reads only `arm=='tm_score'` rows and only `n_detections`, `full_list_recall`, `n_gt_mitotic`. No path exists by which the `od` fix touches a `tm_score` row — `sub` selects identical rows, `_rank` sorts by `'score'`, and `nan_rate` is already 0.0 on every `tm_score` row. |
| Changing `CURRENT_Z` in v1 or v2 | Each documents the floor in force when it ran; changing it silently moves the comparison snapshots. A11/B7 clarify instead. |
| Un-hardcoding lcc's `v2_reference` | All 14 transcribed values verified correct today → F2. |

# Follow-ups (not this pass)
* **F1.** Re-run lcc-vs-v2 over `SEED_INDEX` 0–4 so reading-depth deltas get an error bar.
  The only change that would let deliverable (d) support a domain-level claim.
* **F2.** Persist `tp_fp_chromatin_*`'s `sep_table` to `results/`; have lcc read it.
* **F3.** Decide whether `chromatin.score_detections`' silent border NaN should become an
  explicit error or an internal pad — the trap that produced Finding 3.

# Verification (V)
* **V0.** E0 baseline in place.
* **V1 — a gate, not a post-hoc check.** Re-run v1, diff its CSV against the baseline.
  v2 cells 21/23 read v1's CSV off disk for every before/after delta, so **if the diff is
  non-empty, stop and diagnose before editing v2** — otherwise B4's rewritten prose is
  written against moved numbers. This is also the safe interruption point.
* **V2.** Re-run v2; diff — expect `arm=='tm_score'` rows identical, only `chromatin_od`
  read-depth + `nan_rate` moved.
* **V3.** `nan_rate == 0` on all v2 rows.
* **V4.** v2's `201.tiff` `chromatin_od` `read_95` == **814** (z=0.5) and **801** (z=1.0),
  and its `tm_score` rows still match lcc's exactly — both arms, not just one.
* **V5.** Re-run the two `tp_fp_*` cross-checks; expect 63/63 unchanged.

# Review outcomes
19 issues raised. 5 blockers accepted (missing baseline; A4's assertion would have raised;
wrong min-pool figure; a false claim A8 would have written into v1; B9 targeting a
non-existent comment). 11 should-fixes accepted, incl. 3 edits the draft omitted entirely
(B12, C13, and C1's second location). 3 nits accepted. Independently re-verified before
acceptance: the 3,070 minimum, the K=5000 counterexample, the 56-group/0-violation vs
33-violation arm split, the two comment locations, `chromatin.py`'s git state, `ctx` column
leakage, and the deep-pool ratios behind C13.
