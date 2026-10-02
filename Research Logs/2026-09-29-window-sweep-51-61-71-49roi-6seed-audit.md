# Audit of `window_sweep_51_61_71_49roi_6seed.ipynb`: every number reproduces from the CSVs, and the headline arms reproduce from the pixels; the `tm_score` window gain is real and survives the blank D11's wording asks for, but is a disclosed post-hoc exception to the notebook's own rule; and Q2's stated reason for staying open, a self-match, moves the 71 px difference-in-differences by 0.17 of 8.30 points

**Scope.** Target: `production_hematoxylin_only/window_sweep_51_61_71_49roi_6seed.ipynb` (23 cells, 12 code, 2
figures). It writes `results/precision_at_k_49roi_6seed_windowsweep_postD11_{per_run,per_roi,delta_per_roi,delta_stats,top30}.csv`
and `production_hematoxylin_only/windowsweep_{contrasts,size}.png` [files since removed; the figures are embedded in `window_sweep_51_61_71_49roi_6seed.ipynb`]. It reads the committed parent run
`results/precision_at_k_49roi_6seed_fillgeom_postD11_{clicks,per_run,top30,delta_stats}.csv`. Documents read: its
pre-registration `WINDOW_SWEEP_PROMPT.md` (commit `480e2f3`), and `DECISIONS.md` D1, D3, D5 (with its 2026-09-08
amendment), D7, D9 and D11 (with both amendments, 2026-09-17 and 2026-09-22). Audit script:
`window_sweep_51_61_71_49roi_6seed_audit.py` (repo root). Tables: `results/window_sweep_51_61_71_49roi_6seed_audit_*.csv`.
**Everything below is re-derived from `databases/MIDOG++.json`, the TIFF pixels and resolution tags in
`images/extra_valid/` and `images/extra_valid/testing_set/`, the notebook's persisted `per_run` and `top30` tables,
and the parent run's `clicks`, `per_run` and `top30` tables.** No `midog_utils` import appears anywhere in the audit
script. The matcher, the sign-flip, Holm, Spearman, the template construction, the peak extraction, the NMS and the
`od` statistic are all written there, against `cv2`, `skimage`, `tifffile`, `numpy`, `pandas` and `scipy` directly.

## Conflict of interest

- **Signals I can see.** The notebook, its five CSVs and its two figures are **untracked** (`??`), so there is no git
  author or age to attach to them; `git log -1 -- <notebook>` returns nothing. The kernel ran 13:34–15:08 today
  (cell metadata, UTC 17:34–19:08); the file was last saved 15:34. Its pre-registration was committed 6 hours ago by
  `mohini-anand` with a Claude co-author line. A user-memory note, written by an assistant session at 19:39 UTC today,
  already summarises this notebook's result. Every imported module is committed and clean (`git status --porcelain
  midog_utils/` is empty). The work and this audit are probably the same project hand; only my fresh context
  separates them.
- **What limits the conflict.** No table is checked against itself:
  - Buckets come from my own greedy matcher run against the JSON, never from the notebook's `bucket` column.
  - Every p comes from my own integer convolution, and every Holm from my own step-down.
  - Tier B templates and searches come from my own implementation on the pixels.
  - The one thing I reuse by necessity is the bootstrap's numpy stream (`default_rng([20260917, stream])`), because
    a percentile interval reproduces bit for bit only from the same draws. The exact sign-flip and the t interval on
    48 df sit beside it as resampling-free checks.
- **What it cannot cover.** The design itself: the fixed 51 px blank, the 60 px near-seed rule, the joint-gate click
  draw, precision at K = 10, the D5 bar, and the ROI as the unit. Those go to `premise-reviewer` and to a reader who
  did not make them. Part 4 lists which of my verdicts lean on them.

## Part 0 — what reproduces

**The published statistics are the statistics the code and data in the working tree compute.** Provenance is clean.
The modules (`template_match.py` and `production.py` 2026-09-17, the rest 2026-09-14 to 09-16, all committed, all
unmodified) predate the pre-registration (12:45 today). The pre-registration predates the run (13:34–15:08), and the
run predates the artifacts (15:06). No `.partial` file exists. The artifacts are untracked, so they have no committed
provenance. Every row of them is re-derived below from its own inputs: re-matched against the JSON, re-summed and
re-tested. A subset is also regenerated from the pixels:

- templates for all 7 arms on all 294 clicks;
- the `fill_61` and `fill_71` searches on all 294 clicks;
- `default_61` and `default_71` on 39 clicks;
- `default_51`, `fill_51` and `fill_71_pinned` on 4 clicks.

The remaining search rows were checked for consistency and re-matched, not searched again. That matters for Q4
(claim 21), which rests on `fill_71_pinned`: its search was reproduced on only 8 lists, although its template and
fill value were reproduced on all 294 clicks. Tier A results are **independent reproductions** of the arithmetic, not
consistency checks against a stale artifact.

| check | tier | values compared | divergences |
|---|---|---:|---:|
| every top-30 detection re-matched against `MIDOG++.json` (bucket, `matched_ann_id`) | A | 246,960 | 0 |
| `per_run` TP@10/20/30 from my buckets | A | 12,348 | 0 |
| W = 51 parity from the artifacts: `default_51` / `fill_51` top-30 `cx, cy, bucket, ann, score (f32), od (exact)` and `per_run` TP / `n_detections` / `n_peaks` / `tpl_sha1`, against the parent's `default_51` / `hem_rect_fill` | A | 218,736 | 0 |
| `per_roi` sums and precision columns | A | 5,488 | 0 |
| `delta_per_roi`: the ROI delta and its six per-click deltas, 46 contrasts × 49 ROIs | A | 15,778 | 0 |
| `delta_stats`: 21 columns × 46 rows (pooled Δ, W/L/T, n_nonzero, exact p, min p, bootstrap CI and shares, t CI, Holm, floors, D5) | A | 966 | 0 |
| the notebook's own validation (cell 4): the parent's 10 `delta_stats` rows from the parent's `per_run` | A | 120 | 0 |
| near-seed columns (Chebyshev / Euclidean, all / non-TP), from `top30` + clicks | A | 16,464 | 0 |
| `tp_at_10_excl_border` (band dropped, re-matched afresh) and `n_border_band_top10` | A | 8,232 | 0 |
| one search per (click, arm): search-level columns identical on both axes | A | 14,406 | 0 |
| `od` identical at the same (click, `cx`, `cy`) across all 7 arms and both axes (pad 25 / 30 / 35), in groups | A | 38,635 | 0 |
| every printed table and prose number (pooled precision, near-seed, NMS, churn, mechanics, recorded baselines, re-seed stability, Q3, timing, self-match table) and the composition gate | A | 6,745 | 0 |
| B1: all 7 templates × 294 clicks rebuilt from the pixels (SHA-1, `fill_value` f32, `n_kept_px`, `n_filled_px`), plus mpp from the TIFF tags | B | 8,281 | 0 |
| B2: searches reproduced from the pixels (`fill_61` and `fill_71` on all 294 clicks, `default_61` / `default_71` on 39 clicks, all 7 arms on 4 clicks); per list TP@10/20/30, `n_detections` and near-seed; top-30 `cx, cy, bucket, ann, score (f32 exact), od (exact)` | B | 252,216 | 0 |
| **total** | | **845,375** | **0** |

Two printed values differ from the recomputation only by format, and neither is a divergence. The median keep region
is 612.5 px, printed `612` by `:.0f`. The printed Q3 table rounds as its format string says.

**Gates.**

- **Execution coherence: passes.** Code-cell counts run 1–12 contiguously, with no error cell and no unrun cell, and
  the last cell is markdown. One Tier 3 note (Part 1): cell 21 carries count 12 but started 208 s after cell 19
  finished, while every earlier cell started within 0.0 s of its predecessor. As its markdown says, it was executed on
  its own, and the count hides that session boundary. It is self-contained (it re-imports and reads only the saved
  CSVs), so nothing printed above it depends on it.
- **Provenance: clean** (above).
- **Composition: passes, 19 items, 0 failures.**
  - `per_run` holds 4116 rows = 49 × 6 × 7 × 2, with no duplicate key, and 294 rows per (arm, axis).
  - `top30` holds 123,480 rows, and every list ranks 0–29.
  - Click fields match the clicks CSV on all 294 × 6.
  - The domain label matches `MIDOG++.json`'s `tumor_type` on 49/49 ROIs.
  - **7 domains × 7 ROIs.** The unit (ROI) is not collinear with the stratum; the database itself holds 44–150
    annotated ROIs per tumour type.
  - The bootstrap ROI order is the parent's first-appearance order, as pre-registered.
- **Every claim had a persisted artifact under it** except one: the pilot's side of "matching the pilot" (cell 16 /
  cell 22). The pilot numbers exist only in the prompt and the commit message of `480e2f3`, and no pilot CSV exists.
  The claims computed live from pixels had Tier B checks: template SHA-1s, fill values, the wide arms' searches, the
  R_W identity and the self-match mechanism. The two figures had both steps of the figure protocol.
- **Tier C: not run.** Nothing in Tier A or B diverged, so there was nothing to explain. No sweep was re-run through
  the notebook; the measurements that touch pixels are Tier B §B1–§B5.

## Part 1 — findings

### Tier 1

**1. Q2 at `tm_score` is held "unresolved" because of a contamination that the notebook's own artifacts bound at
0.17 of 8.30 points. The same bound shows the prior "precision tracks footprint size" is robustly contradicted,
not "contradicted on contaminated arms".**

The notebook's reading (cell 22):

- the `tm_score` difference-in-differences "cannot be read at face value";
- "this is the contrast the artefact contaminates directly, and in its direction. It stays unresolved";
- "Q2 is not answered cleanly";
- on the prior: "the self-match is one contributor to that drop. This notebook does not separate it from the rest."

It can be separated, and the notebook already holds what is needed:

- Every self-match detection is `non_human_findings`: unmatched, so it claims no ground truth. Dropping it from a
  greedy best-first list therefore changes no other detection's bucket, and the top-30 list gives the new top 10
  exactly (§A6).
- The only channel the top 30 cannot see is NMS un-suppression. A self-match peak at the click could have suppressed
  a real neighbour within one radius, and any such neighbour must sit within 2 r of the click. `MIDOG++.json` has
  **no** other mitosis within r of any of the 294 clicks, and **12** clicks with one within 2 r.
- Tier B §B4 re-ran `default_61` and `default_71` with the self-match peak removed **before** the 100-peak cap and
  NMS. It covered 39 clicks: the 36 with a visible self-match or a mitosis within 2 r, and 3 more from the all-arm set.

| `tm_score`, K = 10 | as recorded | self-match dropped (top 30) | self-match peak removed pre-NMS (Tier B, exact) | worst case (every freed top-10 slot a TP, plus one TP per click with a mitosis within 2 r) |
|---|---|---|---|---|
| `default_61 − default_51` | −1.94 [−3.37, −0.58] | −1.94 | −1.94 (unchanged) | −1.43 [−2.86, −0.10] |
| `default_71 − default_51` | −6.29 [−8.71, −3.98] | −6.12 [−8.47, −3.88] | −6.12 [−8.47, −3.88] | −5.31 [−7.55, −3.16] |
| DiD at 61 | +3.54 [+2.31, +4.86], 31/6/12 | +3.54 | +3.54 (unchanged) | +3.03 [+1.84, +4.35], 30/9/10, D5 met |
| DiD at 71 | +8.30 [+6.22, +10.54], 41/2/6 | +8.13 [+6.09, +10.31], 41/2/6 | +8.13 [+6.09, +10.31], 41/2/6 | +7.31 [+5.37, +9.39], 40/5/4, D5 met |

At `default_71` / `tm_score`, 17 top-10 slots hold a self-match, 8 of them at rank 0, where the tool would return the
clicked cell itself as its first suggestion. They hid **5** true positives, not 17: the detection that moves up into
rank 10 is a TP 5 times. At `default_61` the 3 slots hid **0**. So the self-match raises the 71 px DiD by 0.17 points
and the 61 px DiD by nothing, where the notebook says "it raises both". The exact Tier B rerun (§B4) gives the same counts: +5 TP at `default_71` / `tm_score`, +3 at `default_71` / `chromatin_od`, and 0 at `default_61`. The NMS channel added nothing on any of the 39 clicks, including the 12 with a mitosis within 2 r.

The worst case leaves both DiDs positive with D5 met, and leaves `default_71 − default_51` negative with its CI
excluding 0. On `tm_score` it is a true bound:

- Ranking is by the NMS key, so a peak the self-match un-suppresses enters below the self-match's own rank.
- The 101st peak, admitted by a freed pool slot, ranks last.
- The one remaining channel is the blank's effect on windows that overlap it, and it shows in the near-seed counts.
  `default_71` has 2 near-seed true positives in its `tm_score` top 10, against `default_51`'s 3. That is one TP,
  0.03 points, inside the margin above.

**What is supported:** at `tm_score` the change is specific to the fill template. The DiD is +8.30: +8.13 exactly (§B4)
without the self-match, and at least +7.31 in the worst case. As the notebook correctly says, it is driven mostly by
the untouched window *losing* precision. The prior predicted a positive `default_W − default_51`; it is negative
beyond any contamination this design can produce.

Keeping the pre-registered flag on Q2 is the user's decision, recorded in the notebook. It stands, and under it Q2
stays labelled unresolved. My finding is narrower: **the evidential reason the notebook gives for Q2 being unreadable
does not hold.** The asymmetry is worth stating. The flag was lifted for Q1 because the artefact is absent from Q1's
arms, and kept for Q2 because the artefact is present in Q2's arms. But an artefact present at 0.17 of 8.30 points is
not a reason the contrast "cannot be read". The notebook proposes a follow-up, a blank covering the untouched
template's footprint "to measure the untouched window cleanly". That follow-up is not needed to answer Q2's sign or
size at `tm_score`.

*Caveats.* At `chromatin_od` the worst case is not a true bound. Ranking there is by `od`, not by the NMS key, so a
freed pool slot (the 101st peak by score) could enter the `od` top 10 anywhere in the ROI, and nothing in the
artifacts bounds that outside the 39 re-run clicks. I therefore claim no bound on the `chromatin_od` DiDs. They are
non-significant as recorded (+0.24 and +0.58), and the notebook claims nothing for them. This finding rests on the
ROI as the unit (Part 4).

### Tier 2

**2. The Q1 `tm_score` headline is a post-hoc exception to the notebook's own pre-registered rule. The exception is
well founded, but nothing durable records it.**

Cell 0 and the prompt fixed the rule before the run: "A rise above 0 under either definition on a wide arm makes every
`fill_61`/`fill_71` gain on that axis unresolved." It fired on both axes: `default_61` 3 and `default_71` 17 at
`tm_score`, and `default_71` 4 at `chromatin_od`. Under the pre-registration's own rule, Q1 at `tm_score` therefore
lands in *unresolved*, not in *D5 met*.

The notebook says so, and reports the user's post-run decision to lift the flag (cell 22). That is disclosure, not
concealment, which is why this is Tier 2. But "D5 met (flag lifted)" is not a confirmatory result. The lift is
recorded only in the notebook's prose and in a user-memory note, not in `DECISIONS.md` or any `Research Logs/`
entry: a grep of every `*.md` in the repo for `lifted`, `flag.*lift`, `window sweep`, `windowsweep` and
`window_sweep` finds only the prompt.

I checked whether the lift's rationale holds rather than taking the counts:

- **The near-seed top-10 detection *sets* match exactly, by (`cx`, `cy`), between each wide fill arm and `fill_51`.**
  - At `tm_score`: 294/294 clicks for `fill_61`, `fill_71` and `fill_71_pinned`.
  - At `chromatin_od`: 294/294, 294/294 and 293/294. The one difference (`353.tiff` s2, `fill_61`) is the same
    mitosis (ann 18280), detected 1 px lower.
- **A fill template cannot self-match.** At the click its keep pixels all fall on the blank, so
  `R = k (m_k − f_W)(min − μ_W) < 0`. Tier B §B3 measured the `fill_71` response at the click as negative on
  49/49 clicks, the largest being −0.118.
- The pre-registered target, μ_W inflation of the fill arms near the blank, left no trace in any fill arm's top 10.

So the rule's purpose was met on Q1's arms. What fired it was a different mechanism, the ring self-match (confirmed in
closed form, §B3), on other arms. The headline stands as a disclosed post-hoc exception on a sound basis. It belongs
in a dated `DECISIONS.md` note or a results log, so that a later reader never meets "D5 met" without "after the rule
fired" beside it.

The ROI-level result behind it is robust:

- Leave-one-ROI-out keeps `fill_61 − fill_51` within +1.46 to +1.74 (max Holm p 0.0023), and `fill_71 − fill_51`
  within +1.77 to +2.19 (max Holm p 0.0085).
- No single ROI carries more than 14% of the net; the top three carry 28% and 37%.
- All 6 click draws are net positive for both contrasts: +8, +8, +3, +16, +2, +10 TP, and +10, +9, +6, +9, +9, +16 TP.

**3. "At `chromatin_od` … costs 1.16 points at 71" (What this means) is stated as an effect. On the project's own bar
it is not established.**

The arithmetic reproduces: −1.16 [−2.14, −0.27], Holm 0.0385. The notebook's own Q1 section correctly says D5 is not
met, since 22/49 ROIs lose, short of the 25 required. The closing bullet drops that qualifier. On the pre-registered
bar the result is "CI excludes 0, D5 not met", which the project does not read as an established effect. Two
robustness checks inside the pre-registered design point the same way:

- Dropping `293.tiff` alone lifts its Holm p to 0.077.
- Excluding the edge band gives exact p 0.072 and a t CI of [−1.94, +0.11].

Two sensitivities *outside* the pre-registered design, shown for context only:

- Holm across all 8 tested rows, instead of the pre-registered 4 families of 2: 0.077 (`..._multiplicity.csv`).
- The domain as the unit: p = 0.094 (Part 4).

For balance: the loss is not one click draw. All six draws are net negative (−1, −5, −8, −7, −5, −8 TP). The supported
statement is that at `chromatin_od` widening the window does not help at either width, and a loss at 71 is suggested
but not established by D5. The `tm_score` Q1 contrasts pass every check above.

### Tier 3

4. **The 51 px blank departs from D11's wording, not from its rationale, and the `tm_score` gain survives the blank the wording asks for. Under the rule fixed below, this is resolved.**

D11 (applied 2026-09-17, kept 2026-09-22) blanks "a `base_size x base_size` square … the same footprint the template
itself was cut from". Its stated reason is that "the seed's own match simply has no source pixels left to correlate
from". The prompt fixed the blank at 51 for every arm, so the positions diverge three ways:

- D11's wording: blank = the template footprint, which is W.
- The prompt: 51 for every arm.
- The notebook: 51 for every arm.

For the *untouched* 61 / 71 px templates, 51 violates D11's rationale: the ring has source pixels, and it
self-matches (finding 1). For the *fill* templates it does not, because the only non-constant pixels are the keep
region, which sits inside the 51 px blank, and the response at the click is negative (§B3). D11's rationale is met at
51 for the fill arms; only its wording asks for W. "What this means" nevertheless frames the gain as a production
statement ("at `tm_score` (the production default)"). Tier B §B5 therefore measured the wording's version: `fill_61`
and `fill_71` on all 294 clicks with a 61 / 71 px blank filled with the channel minimum.

The response was patched as `R + corr(D, T − mean T)`, where D is the ring difference image. `TM_CCOEFF` is linear in
the image, so this is exact; §B5's own check against a full recompute on `509.tiff` agrees to float32 rounding (max
|Δ| 3.6 × 10⁻⁶ on a response scale of 4.9), with identical top-100 peaks in the same order. The contamination check counts near-seed non-TP detections out to the window's full reach:
60 px for 61, and **70 px** for 71 (35 + 35), not the notebook's 60.

**Decision rule, fixed before I read the final B5 table.** A two-ROI development pass had already shown no change,
so this is a rule written with a hint in hand, not a blind one. If the D11 `fill_71 − fill_51` and `fill_61 − fill_51` each keep
a bootstrap CI excluding 0, meet D5, and show 0 near-seed non-TPs at full reach, this is a Tier 3 scope observation,
resolved. If either loses its CI or D5, it becomes Tier 1 against the "production default" framing. (The rule concerns
the `tm_score` headline; at `chromatin_od` neither contrast met D5 as recorded. That clarification was added after
reading the table.)

| K = 10 | blank | near-seed reach | near-seed non-TP in top 10 | Δ, pts [boot CI] | W/L/T | exact p | D5 |
|---|---|---:|---:|---|---|---:|---|
| `tm_score` `fill_61 − fill_51` | 51 (as run, reproduced) | 60 | 0 | +1.60 [+0.82, +2.41] | 29/10/10 | 0.00058 | met |
| | 61 (D11 wording) | 60 | 0 | +1.60 [+0.82, +2.41] | 29/10/10 | 0.00058 | met |
| `tm_score` `fill_71 − fill_51` | 51 (as run, reproduced) | 60 | 0 | +2.01 [+0.75, +3.30] | 27/10/12 | 0.0044 | met |
| | 71 (D11 wording) | 70 | 0 | +2.01 [+0.75, +3.30] | 27/10/12 | 0.0044 | met |
| `chromatin_od` `fill_61 − fill_51` | 51 (as run, reproduced) | 60 | 0 | −0.58 [−1.53, +0.14] | 10/18/21 | 0.24 | not met |
| | 61 (D11 wording) | 60 | 0 | −0.58 [−1.53, +0.14] | 10/18/21 | 0.24 | not met |
| `chromatin_od` `fill_71 − fill_51` | 51 (as run, reproduced) | 60 | 0 | −1.16 [−2.14, −0.27] | 9/22/18 | 0.019 | not met |
| | 71 (D11 wording) | 70 | 0 | −1.19 [−2.24, −0.27] | 9/22/18 | 0.019 | not met |

The wording's blank changes TP@10 on **1 of 1176** fill lists (net −1 TP, on `chromatin_od` `fill_71`). Out to full
reach, the near-seed top-10 count is unchanged (3 at `tm_score`, 4 at `chromatin_od`), and none is a non-TP. Under
the rule fixed above, both `tm_score` contrasts keep their CI, D5 and a clean near-seed count. **This is a scope
observation, resolved: the headline is not an artefact of the 51 px blank.** One point of framing remains. The fill
template is not production's template, so "the production default" names the ranking axis, not a configuration
production runs.

5. **Q3's emphasis.** All Q3 numbers reproduce (22 printed values, 0 divergences). The prose leads with the
   between-ROI "larger-nucleus ROIs gain more (or lose less) … That is the pilot's direction." Between ROIs, size is
   largely species (Spearman −0.692 with canine), and at `tm_score` human ROIs gain +3.10 points against +1.19 for
   canine. The confound-free within-ROI association at `tm_score` has the **opposite sign**, −0.183
   [−0.316, −0.042]: within a ROI, clicks on larger nuclei gain less. The notebook prints both numbers and the
   confound beside them, and it concludes "no design decision follows from Q3". So this is emphasis, not a weakened
   claim, but the sign reversal deserves a sentence of its own.
6. **Cell 21's count masks a session boundary** (Part 0, execution gate).
7. **Figure 1's two panels use different x scales.** The data span −3.20 to +2.21 points at `chromatin_od` and
   −8.71 to +12.99 at `tm_score`, each auto-scaled. Every annotation matches the recomputation, but equal bar lengths
   across the two panels are not equal effects.
8. **Figure 2 over-plots at y = 0.** 18 of 49 points at `chromatin_od` and 12 of 49 at `tm_score` sit on the zero
   line, under ±0.25 px jitter. The series has all 49 points (28 canine, filled; 21 human, hollow), but the render
   cannot show them all.
9. **Timing.** The loop took 91.7 min (the per-search times sum to 84.5 min), 1.5–1.8× the prompt's 50–60 min
   budget.
   - The prompt's premise that `matchTemplate` gets monotonically faster with W (3.02 / 2.56 / 1.98 s) is not what
     the run measured: median per-search 2.27 / 2.52 / 1.96 s.
   - Each is a median over 588–882 searches on a machine shared with other processes, so it is a throughput under
     contention, not a clean benchmark.
   - My re-timing on `300.tiff` (3 repeats): load + `rgb2hed` 5.28 s; `matchTemplate` including the replicate pad 2.12 / 2.83 / 1.88 s at W = 51 / 61 / 71 (`..._retiming.csv`). That is slower at 61 than at 51, as the run found, so the prompt's monotone premise fails again.
   - No conclusion rests on timing.

**Observations that are not defects.**

- **In-sample selection (Step 4.10) is not material here.** Leave-one-ROI-out choice among all 9 configurations the two
  runs compared (4 at W = 51 from the parent, 5 here) picks `fill_71` on 49/49 folds. Its held-out gain over
  `default_51` equals the in-sample +3.81 points (p 3.5 × 10⁻⁶). Restricted to the fill window, the held-out gain
  over `fill_51` equals the in-sample +2.01. The untested `fill_61 − default_51` the notebook points to is +3.40
  [+2.28, +4.56], 36/8/5 (descriptive).
- **Product relevance (Step 4.8): no worst-seed regression at `tm_score`.**
  - `fill_71 − fill_51`: 85 clicks up, 47 down. The worst click loses 2 TP (4 such clicks); 21 clicks gain ≥ 2.
  - Lists with TP@10 = 0: 14 (`fill_51`), 12 (`fill_61`), 11 (`fill_71`).
  - The gain is fewer unannotated false positives (33.88% → 31.94% of top-10 slots). The look-alike share is flat
    (13.47% → 13.40%).
  - The untouched `default_71` goes the other way: 34 zero lists, and 44.25% unannotated false positives.
- **"Most of the gain has arrived by 61."** At `tm_score`, `fill_71 − fill_61` = +0.41 [−0.34, +1.16], p 0.34. At
  `chromatin_od` it is −0.58 [−1.05, −0.14], p 0.028. Both are descriptive and consistent with the prose.
- **The R_W identity (cell 0) holds numerically.** `Σ_keep (t − f_W)(I − μ_W)` and its decomposed form match cv2's
  `TM_CCOEFF` response to a relative 1.1 × 10⁻⁵ on 294 positions (§B3).
- **Config drift: none against production.**
  - `CHANNEL`, `TM_METHOD`, `PEAK_MIN_DISTANCE`, `DEEP_FLOOR_Z`, `MAX_PEAKS` and `OD_WINDOW` equal `production.py`.
  - The NMS radius equals the match radius, 7.5 µm (D7).
  - `FSConfig`'s `tm_method` default is now `cv2.TM_CCOEFF`, agreeing with D1.
  - The only drift is the blank (finding 4), and it is by design. Table: `..._config_drift.csv`.
- **Leakage: none.** The click's own annotation is removed from the ground truth, by `seed_ann_id`, on every list. My
  matcher does the same independently and reproduces all 246,960 bucket and annotation values. The self-match is the
  only self-hit, and it is finding 1.
- **Definitional drift: none found.**
  - `ds.image_annotations(annotations, file_name)` is called with `category_id=None`, which is correct here: the
    matcher needs both categories to bucket look-alikes, and TP is `category_id == 1`.
  - `n_gt_mitotic` counts category 1 only.
  - `n_nms_suppressed == n_peaks − n_detections` on 4116/4116 rows.
  - The notebook asserts its four invariants inline; none is persisted as a `*_verification.csv`.
- **Not applicable:**
  - recall triad (Step 4.3): no recall is claimed, and precision is at a fixed K with `n_detections` ≥ 79 on every
    list;
  - length-matched precision null (Step 4.4): every contrast compares arms at the same K on the same clicks;
  - dose confounded with domain (Step 4.2): there is no dose, and every ROI receives every arm;
  - method machinery (Step 4.9): none is layered over known labels.

## Part 2 — verdict per conclusion

| # | cell | claim (abridged) | verdict | tier |
|---|---|---|---|---|
| 1 | 0 | TM_CCOEFF on a fill template = Σ_keep (t − f_W)(I − μ_W); the first term is W-independent | reproduces (§B3) | — |
| 2 | 0 | `od` is identical for every pad by construction | reproduces (38,635 groups) | — |
| 3 | 0, 2 | recorded W = 51 run: 0 near-seed non-TP (28 / 16 slots); edge < 35 px = 4.70%, TP 0.246 vs 0.560; `MAX_PEAKS` 2352/2352 | reproduces | — |
| 4 | 1 | configuration equals `production.py` | reproduces | — |
| 5 | 4 | stats code reproduces the parent's 10 `delta_stats` rows; convolution = enumeration | reproduces (120 values) | — |
| 6 | 8, 22 | W = 51 parity 1176/1176; `od` bit-exact 35,280/35,280; asserts 1–4 on 2058 searches | reproduces (from artifacts, and from pixels §B2) | — |
| 7 | 10 | `MAX_PEAKS` binds on every row; NMS table; near-seed table; flag fires on both axes | reproduces | — |
| 8 | 10, 22 | mechanics: mass weight −2.3% / −4.4%, filled share +9.3% / +14.9%, count +56.3% / +122.7%; f_71 ≠ f_51 on 294/294; 7 distinct templates on 294/294; "the prompt's ~28% does not match" | reproduces (§B1 from pixels) | — |
| 9 | 12 | pooled precision table | reproduces | — |
| 10 | 14 | 8 tested + 38 descriptive `delta_stats` rows; re-seed stability | reproduces | — |
| 11 | 16, 22 | churn 69 / 132 clicks, overlap 5.49 / 5.57, net −34 / +59; pilot 9/12, 5.92, +4, −3/+7 | reproduces (this run's side); **cannot check** the pilot's side, because no pilot artifact exists | — |
| 12 | 21, 22 | the near-seed rise is only in `default_61` / `default_71`, every one 0–2 px from the click, a seed self-match; 24/17 and 13/4 lists | reproduces (closed form, §B3) | — |
| 13 | 22 | "it lowers them, and through them raises both DiDs" | reproduces, overstated: +0.17 pts at 71 and 0 at 61, from the top 30 (A6) and from the exact rerun (§B4) alike | 1 |
| 14 | 22 | Q1 `tm_score`: +1.60 / +2.01, Holm 0.0012 / 0.0044, D5 met (flag lifted) | reproduces, overstated: a post-hoc exception to the pre-registered rule (finding 2); it survives the D11-wording blank (finding 4) | 2 |
| 15 | 22 | the three facts behind the lift | reproduces; the near-seed *sets* are identical too (294/294 at `tm_score`) | — |
| 16 | 22 | Q1 `chromatin_od`: no difference at 61; −1.16 at 71, Holm 0.0385, D5 not met | reproduces | — |
| 17 | 22 | "the two axes disagree in sign" | reproduces at 71; at 61 the `chromatin_od` CI includes 0 | — |
| 18 | 22 | Q2 `tm_score` DiD "cannot be read at face value … stays unresolved" | reproduces as a rule outcome (the user's decision to keep the flag); the stated evidential reason does not hold: the self-match moves the DiDs by 0.17 and 0 pts (exact, §B4), and by at most 0.99 / 0.51 in the worst case | 1 |
| 19 | 22 | Q2 `chromatin_od`: no difference detected; no equivalence claimed | reproduces | — |
| 20 | 22 | "the stated prior is contradicted, on contaminated arms … not separated" | reproduces, overstated caution: `default_71 − default_51` is −6.12 (exact, §B4) without the self-match and −5.31 or lower in the worst case, CI excluding 0 | 1 |
| 21 | 22 | Q4: pinning changes almost nothing; the pinned arm keeps most of the window effect; the effect runs through μ_W | reproduces (descriptive; at `chromatin_od` the pinned arm keeps more than all of it, −1.46 vs −1.16) | — |
| 22 | 18, 22 | Q3 numbers (between / within, species split, recorded facts) | reproduces | — |
| 23 | 22 | Q3: "larger-nucleus ROIs gain more … the pilot's direction" | reproduces; the emphasis omits that the within-ROI sign reverses at `tm_score` (finding 5) | 3 |
| 24 | 22 | Q3: no design decision follows | reproduces | — |
| 25 | 22 | edge band: shares 6.19 → 6.56 → 7.01% and 3.20 → 3.30 → 3.47%; signs unchanged; `chromatin_od` loss borderline | reproduces | — |
| 26 | 22 | What this means: `tm_score` +1.6 to +2.0, D5 met, via μ_W, most of it by 61 | reproduces, overstated (as 14) | 2 |
| 27 | 22 | 0.5085 → 0.5265 → 0.5425; `default_51 → fill_61` untested | reproduces (+3.40 descriptively; LOO shows no selection optimism) | — |
| 28 | 22 | What this means: `chromatin_od` "costs 1.16 points at 71"; "a window choice made on one axis does not carry to the other" | reproduces, overstated: D5 not met (finding 3) | 2 |
| 29 | 22 | What this means: "Q2 is not answered cleanly … inflates the DiD" | reproduces, overstated (finding 1) | 1 |
| 30 | 22 | standing caveats (joint-gate draw, one channel, K = 10 only, no domain analysis) | reproduces | — |
| F1 | 19 | Figure 1: every K = 10 contrast with its CI, W/L/T and Holm | reproduces: series re-derived (22 rows), render read; unequal panel scales (Tier 3) | 3 |
| F2 | 19 | Figure 2: per-ROI `fill_71 − fill_51` against nucleus size, with the titles' rho and CIs | reproduces: series re-derived (49 + 49 points), render read; over-plotting at 0 (Tier 3) | 3 |

`cannot check` appears once, for the pilot's side of claim 11. One instance is not a reproducibility pattern.

## Part 3 — what was re-run versus read

- **Tier A (every claim with a persisted artifact under it):** A0–A14 in `window_sweep_51_61_71_49roi_6seed_audit.py`.
  - All statistics use plain numpy / pandas / scipy.
  - The exact sign-flip is an integer polynomial convolution (`signflip_p`). The bootstrap draws from the same
    documented numpy stream, which is necessary to reproduce a percentile CI. The domain-level null is a full 2⁷
    enumeration.
  - At G = 49 the cluster bootstrap is within its validity range, and I report the notebook's exact sign-flip and its
    t on 48 df beside it.
  - The attainable floor for the tested rows is 2^(1 − n_nonzero) ≤ 7.5 × 10⁻⁹, so no tested row is floor-blocked.
- **Tier B (19.8 min wall clock, 4 workers, on a machine shared with other processes):** I chose these checks
  because the notebook's only parity gate covers W = 51. The wide arms, which carry the headline, had four asserts
  and no reproduction.
  - §B1: all 7 templates for all 294 clicks, rebuilt from the pixels.
  - §B2: `fill_61` and `fill_71` searches on all 294 clicks, plus the other arms on the subsets above, reproduced from
    the pixels, top 30 bit for bit.
  - §B3: the R_W identity at 6 positions on one click per ROI (294 positions, max relative error 1.1 × 10⁻⁵), and the
    `default_71` response at the click against `Σ_ring T′(T − min)` on 39 re-run clicks (max relative error
    2.2 × 10⁻⁶; positive on 38/39).
  - §B4: `default_61` / `default_71` with the self-match peak removed before the cap and NMS, on 39 clicks.
  - §B5: the D11-wording blank for `fill_61` / `fill_71` on all 294 clicks, plus the shortcut's check against a full
    recompute.
- **Tier C:** not run (Part 0).
- **Not audited:** `window_sweep_51_61_71_summary.pptx` (repo root, untracked, saved 15:37, after the notebook). It is
  a downstream copy of these claims. It is worth checking that it carries the post-hoc caveat on claims 14 and 26 and
  the D5 qualifier on claim 28.
- **The script runs end to end.** The final clean run took 23.4 min and exited 0. It regenerates every table
  cited here: 33 files at `results/window_sweep_51_61_71_49roi_6seed_audit_*.csv`, including the ledger, each
  with an mtime inside that run.
- **Operational note.** Earlier runs of the audit script were stopped part-way through Tier B so that checks could be
  added; every number here comes from the final, uninterrupted run. Stopping the first used a pattern-matched
  `pkill -f multiprocessing.spawn`. That pattern can also have matched a pool worker of an unrelated long-running
  Python job on this machine (a `Tractography` scratch script, PID 82491). The job was still alive and spawning
  workers afterwards, but whoever owns it should check its output.
- **Beyond Step 4's named modes, I looked for:**
  - whether the Q1 gain is a near-seed or edge artefact: the near-seed sets are identical, and the edge band leaves
    every sign unchanged;
  - whether it is a rescue of catastrophic lists: no ROI carries more than 14% of the net, and zero-TP lists fall;
  - whether it depends on one click draw: all 6 draws are positive;
  - whether the gain comes from look-alikes or unannotated tissue: unannotated tissue;
  - whether it survives the blank D11's wording asks for (§B5);
  - whether the family definition carries the `chromatin_od` loss: it does (Holm-8 gives 0.077);
  - whether the self-match bound needs the NMS channel: 12 clicks have a mitosis within 2 r; §B4 closes it, with 0
    extra TP.
- **Family of arms.**
  - `git log --diff-filter=A` returns nothing for the notebook, because it is untracked.
  - The prompt's commit `480e2f3` is the only window-sweep commit.
  - `ls results/ | grep -i window` finds only this run's five CSVs.
  - A grep for `windowsweep` / `window_sweep` over `*.ipynb`, `*.py` and `*.md` finds only this notebook and its
    prompt.
  - The pilot (2 ROIs) is the only other arm, and it has no artifact.
  - The parent (`fill_geometry_precision_at_10_49roi_6seed.ipynb`) and its by-domain slice read the same clicks, but
    they are the stage before this one, not arms of it.

  **No family of window-sweep arms is established by these methods.** The selection breadth across the two runs is
  the 9 configurations in the LOO check.
- **Prior audits.** `ls "Research Logs"/*audit*` lists 20. None covers this notebook, its parent `fillgeom_postD11`
  run, or their artifacts, so there is nothing to reconcile, and every finding above is new rather than a
  confirmation or an overturning. `2026-09-23-annotation-frame-spec-review.md` cites the parent notebook as a
  reference, but it is a spec review, not an audit.
- **Appendix facts that no longer hold** (agent file dated 2026-09-12):
  - `FSConfig.tm_method` defaults to `cv2.TM_CCOEFF`, not `TM_CCOEFF_NORMED` (§A14).
  - `results/` holds 584 CSVs (was 233), 33 of them this audit's; excluding those, 5 are untracked (this run's five), where the appendix had every one tracked. The repo holds 85 notebooks (was 59) with 440 embedded PNGs (was 217) (script's final lines).
  - The appendix's timing concerns `to_gray_inverted` (0.45 s). A `hematoxylin_od` Tier B pays `rgb2hed` instead:
    5.28 s (load + `rgb2hed`, median of 3) on `300.tiff` alone, and 4–60 s per ROI under 4 concurrent workers on a machine shared with other jobs (the Tier B per-ROI lines).
  - Still true: `image_annotations`' `category_id=None` default, and `MIDOG_RADIUS_UM = 7.5`.

## Part 4 — premises this audit inherited

**ROIs per stratum, measured: 7 per domain in this sample (49 ROIs, 7 domains), and 44–150 annotated ROIs per tumour
type in `MIDOG++.json`.** Clustering by ROI and clustering by domain are different tests here; Step 4.1 and Step 4.2
are not collinear. Effective units: 37–43 of the 49 ROIs are nonzero on the tested contrasts at `tm_score`.
`211.tiff` has TP@10 = 0 on every arm and click, and carries no signal on any contrast.

**The ROI-premise sensitivity.** With the 7 domains as the units (exact sign-flip over the 7 domain sums, all 2⁷
patterns):

| K = 10 | domains + / − / 0 | domain sign-flip p | floor |
|---|---|---:|---:|
| `tm_score` `fill_61 − fill_51` | 6 / 1 / 0 (canine lymphosarcoma −1.19) | 0.0625 | 0.0156 |
| `tm_score` `fill_71 − fill_51` | 5 / 1 / 1 (canine lymphosarcoma −1.43, canine lung 0.00) | 0.0938 | 0.0313 |
| `tm_score` DiD 61 | 7 / 0 / 0 | 0.0156 | 0.0156 |
| `tm_score` DiD 71 | 7 / 0 / 0 | 0.0156 | 0.0156 |
| `tm_score` `default_71 − default_51` | 0 / 6 / 1 | 0.0313 | 0.0313 |
| `chromatin_od` `fill_71 − fill_51` | 1 / 5 / 1 | 0.0938 | 0.0313 |

This is not a per-domain analysis, which the prompt forbids and which I do not do. At G = 7, Q1 is not resolvable;
at the floor, the DiDs are resolved.

| premise | source | what would falsify it | leaned on by |
|---|---|---|---|
| the ROI is the exchangeable unit | `DECISIONS.md` D5; F5 §8 | domain-level (or scanner-level) heterogeneity making domains the effective unit; the sensitivity above shows the Q1 verdict depends on it | every Q1 verdict; findings 1 and 3 |
| the near-seed rule (60 px, non-TP, any wide arm) as the contamination instrument | `WINDOW_SWEEP_PROMPT.md` | a contamination that changes top-10 content without adding non-TP near-seed detections. I checked detection-set identity, not only counts, and widened the reach to 70 px for the 71 px blank in §B5 | findings 2 and 4 |
| the D5 bar (CI excludes 0, majority of ROIs) as the effect bar | D5 "What would change my mind" (written for axis promotion, reused here as a general bar) | none: it is a convention. A different bar would move the `chromatin_od` verdict, not the `tm_score` ones | findings 1, 2, 3 and 4 |
| NMS radius = match radius = 7.5 µm | D7 | a matching tolerance under which the wide templates' slightly displaced peaks (e.g. 1 px at `353.tiff`) match differently; not observed | all TP counts |
| precision@10 per click as the product metric, with the worst click beside it | D4; the product framing (reading burden per click) | a reader workflow in which K ≠ 10, or in which look-alikes weigh differently from unannotated false positives. The gain here is unannotated false positives, not look-alikes | findings 1–4 |
| the 294 joint-gate clicks represent a pathologist's click | parent run design (`precision_at_k_49roi_6seed_fillgeom_postD11_clicks.csv`: 59 of 294 clicks needed at least one redraw before both Otsu gates accepted) | a click distribution outside the gate's acceptance set (the own-draw arm, which the parent did not run) | every verdict's scope |
| the fixed 51 px blank as the tested configuration | `WINDOW_SWEEP_PROMPT.md` ("by decision"), against D11's wording | D11's wording applied to wide templates, measured in §B5 | findings 1 and 4 |
| two ranking axes are two outcomes, never pooled | the prompt; D5 | none needed; I follow it | all |
| prior audits as corroboration | — | none was leaned on, because none covers this material | — |

The premise that looks shakiest is the ROI as the exchangeable unit for a domain-structured product claim. It is
worth a `premise-reviewer` pass, together with the joint-gate click premise.

**If a premise here is wrong, these verdicts change:**

- If the domain rather than the ROI is the exchangeable unit, claims 14 and 26 (Q1 at `tm_score`) become "not
  resolvable at G = 7". So does claim 28's `chromatin_od` loss. Finding 1's `tm_score` DiDs (claims 18 and 29)
  survive at the floor, 7/7 domains, and so does the prior's contradiction (claim 20: 0 / 6 / 1 domains, p = 0.031,
  the floor).
- If the joint-gate clicks are not representative, every verdict keeps its arithmetic and loses its product scope.
- If the pre-registered near-seed rule is read literally, with no post-hoc lift allowed, claims 14 and 26 return to
  "unresolved". My check that the lift is well founded (finding 2) does not change that reading; it changes only what
  the reading costs.
