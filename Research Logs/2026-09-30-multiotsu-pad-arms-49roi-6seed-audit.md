# Audit of `multiotsu_pad_arms_49roi_6seed.ipynb`: every number reproduces, from the CSVs and from the pixels, and all 10 pre-registered verdicts stand; padding's `tm_score` gain over the bare box survives every check I could add, including the production (D11) blank, while the `chromatin_od` side is thinner, its `mo_pad61` null turning on one ROI and its `mo_pad51` "met" resting on one ROI of majority slack

**Scope.** Target: `production_hematoxylin_only/multiotsu_pad_arms_49roi_6seed.ipynb`. It has 29 cells (15 code) and 2
figures.

- **Writes:** `results/precision_at_k_49roi_6seed_mopad_postD11_{templates,per_run,top30,per_roi,delta_per_roi,delta_stats,selfmatch_sensitivity}.csv`,
  their `_smoke` twins, and `production_hematoxylin_only/mopad_{contrasts,templates}.png`.
- **Reads:** `results/precision_at_k_49roi_6seed_fillgeom_postD11_{clicks,per_run,delta_stats}.csv` (committed, `a9293d0`),
  `results/precision_at_k_49roi_6seed_windowsweep_postD11_{per_run,top30,delta_stats}.csv` (untracked), and
  `production_hematoxylin_only/multiotsurefine14_summary.csv` (committed, `154863e`).
- **Documents read:** its pre-registration (cell 0) and the prompt behind it (`MULTIOTSU_PAD_ARMS_PROMPT.md`,
  untracked); `DECISIONS.md` D5 (with its 2026-09-08 amendment), D7 and D9; D11 (with both amendments, 2026-09-17 and
  2026-09-22); and D8 in `D8_TEMPLATE_ANCHOR.md`.
- **Prior audit, read last:** `Research Logs/2026-09-29-window-sweep-51-61-71-49roi-6seed-audit.md`, which covers the
  reference artifacts.
- **Audit script:** `multiotsu_pad_arms_49roi_6seed_audit.py` (repo root). **Tables:**
  `results/multiotsu_pad_arms_49roi_6seed_audit_*.csv`.

**Family: what else was run before this was written up.** I established it from four sources: the git history
(`git log --oneline`; `git log --all -- '*multiotsu*'`), `grep -rl mopad` over the repo's `.md` / `.ipynb` / `.py`, the
`results/` stems the notebook reads, and file mtimes. The notebook itself is untracked, so
`git log --diff-filter=A -- <path>` returns nothing. It is the fourth step of one sequence on the same 294 clicks:

1. the fill-geometry parent run, which drew the clicks and ran 4 arms (`a9293d0`, 09-22);
2. the adoption of `refine_template_box` on 24 clicks (`154863e`, 09-29), 19 of which are among the 294;
3. the window sweep, 7 arms at W = 51 / 61 / 71 (prompt `480e2f3`; artifacts written 09-29 15:06);
4. this notebook.

Its four arms were fixed in the prompt before the run, and all four are reported with their three tested pad
contrasts. The headline is therefore not a best-of-N pick among pad variants. The references `fill_51` and `fill_61`
entered by matching side, not by outcome.

Excluded as not arms:

- the `_smoke` outputs: the same code on 2 ROIs, bit-identical to the full run;
- the `..._windowsweep_postD11_blankW_*` CSVs written at 13:06 today, during this audit: `default_61` and `default_71`
  only, a window-sweep follow-up.

**Everything below is re-derived from `databases/MIDOG++.json`, the TIFF pixels and resolution tags under
`images/extra_valid/`, the notebook's persisted `per_run`, `top30` and `templates` tables, the window sweep's
`per_run` and `top30`, and the fill-geometry `clicks`, `per_run` and `delta_stats`.** No `midog_utils` import appears
anywhere in the audit script. I wrote the following myself, against `json`, `tifffile`, `cv2`, `skimage`, `scipy`,
`numpy` and `pandas`:

- the matcher, and the ground-truth parse (TLBR box to centre);
- the exact sign-flip (integer DP), Holm, the t interval, and every flag;
- the multi-Otsu box, using `scipy.ndimage` for the morphology rather than `cv2` or `skimage.measure`;
- the template construction, the search (blank, replicate pad, `TM_CCOEFF`, robust floor, 15 px local maxima, top 100,
  greedy distance NMS), and `od`.

## Conflict of interest

- **Signals I can see.**
  - The notebook, its CSVs, its figures and its prompt are all **untracked** (`??`). `git log -1 -- <notebook>` returns
    nothing, so there is no author or age to attach. The file was last written 2026-09-29 23:51, and its artifacts
    23:46.
  - The imported modules are committed and clean (`git status --porcelain midog_utils/` is empty).
  - A user-memory note visible to this session already summarises the result ("padding the new box helps a lot, but
    padded it ties the old Otsu-rect fill; pad choice unanswered").
  - The repository was being worked on **while this audit ran**:
    - `c2b3c21` was committed at 12:30 today (a flowchart; it touches nothing audited here);
    - `window_sweep_51_61_71_49roi_6seed.ipynb` was re-saved at 13:08;
    - new `..._windowsweep_postD11_blankW_*` CSVs appeared at 13:06.

  The notebook and this audit are probably the same project hand; only my fresh context separates them.
- **What limits the conflict.** No table is checked against itself:
  - Buckets come from my matcher run against the JSON, never from the `bucket` column.
  - Every p comes from my own DP, and every Holm from my own step-down.
  - Tier B boxes, templates and searches come from my own implementation on the pixels.
  - The one thing I reuse by necessity is the bootstrap's numpy stream (`default_rng([20260917, stream])`), because a
    percentile interval reproduces bit for bit only from the same draws. Those columns are labelled
    **consistency**. The exact sign-flip and the t interval on 48 df sit beside them as resampling-free checks.
- **What it cannot cover.** The design itself goes to `premise-reviewer` and to a reader who did not make it, not to
  me: the joint-gate click draw, the fixed 51 px blank, precision@10, the D5 bar, the ROI as the unit, and the
  three-bucket nearest-GT matching. Part 4 lists which verdicts lean on each.

## Part 0 — what reproduces

**The published statistics are the statistics the code and data in the working tree compute.** Across the checks
below, **371,880 values compared, 0 divergences**: 333,004 Tier A, 38,682 Tier B, and 194 consistency (bootstrap
bounds on the same RNG stream). The findings in Part 1 are about inference and scope, not about the run.

| check | tier | values compared | divergences |
|---|---|---:|---:|
| composition gate (22 items; see below) | A | 22 | 0 |
| every top-30 detection of the 4 arms and the 3 references re-matched against `MIDOG++.json` (bucket, `matched_ann_id`) | A | 246,960 | 0 |
| TP@10/20/30 per list from my buckets (2352 + 1764 lists) | A | 12,348 | 0 |
| reference provenance: window-sweep `default_51` / `fill_51` against the **committed** fill-geometry `default_51` / `hem_rect_fill` (TP@10/20/30, `n_detections`, `n_peaks`, `tpl_sha1`, `n_blanked_px`), plus SHA against the clicks CSV | A | 9,408 | 0 |
| `per_roi` sums and click counts | A | 3,430 | 0 |
| pooled precision: the cell-19 table (7 arms × 2 axes × 4 metrics) and the six prose values | A | 62 | 0 |
| the notebook's own validation (cell 8): the 10 fill-geometry `delta_stats` rows, from the committed `per_run` | A | 110 | 0 |
| flags, `S`, `n_filled_px` and the Chebyshev offset, from the box geometry; the prose flag table | A | 8,245 | 0 |
| `delta_stats`, all 94 rows: pooled Δ, exact p, floor, t CI, W/L/T, `n_nonzero`, `n_same_direction`, D5; Holm on the 10 tested rows | A | 1,420 | 0 |
| `delta_stats` bootstrap bounds (same stream) and the 20-stream re-seed claim | consistency | 194 | 0 |
| `delta_per_roi` (2940 rows) | A | 2,940 | 0 |
| the cell-28 prose: the 10-test table; the rerun, mass-weight, floor and neighbour numbers; "only one sign changes"; `fill_61 − fill_51` against the sweep's saved row | A | 135 | 0 |
| self-match, near-seed (reach; 60 px Chebyshev and Euclidean), edge band, `near_seed_reach_px`, `MAX_PEAKS`, NMS suppression, from `top30` + clicks + JSON | A | 28,234 | 0 |
| `tp_at_10_excl_border` re-matched afresh from the top 30 (every list had ≥ 10 non-band detections in it) | A | 4,116 | 0 |
| box statistics, the 19-click tripwire, the 6 `box_ge_S` lists and SHAs | A | 392 | 0 |
| determinism: the smoke run's 337/509 rows and top 30 against the full run's | A | 15,168 | 0 |
| figure 1's plotted span; figure 2's eight picks | A | 10 | 0 |
| per-arm search-time medians | A | 4 | 0 |
| **B1**: multi-Otsu box re-implemented on **all 294** clicks (accepted, mode, `bs`, centre x/y); mpp from the TIFF tags | B | 1,764 | 0 |
| **B2**: template SHA-1 rebuilt from pixels, 4 arms × 294 and 3 references × 294; fill value, `m_keep` and mass-weight sign | B | 5,286 | 0 |
| **B3**: own search on the 17 stress clicks × 6 arms (4 + `fill_51` + `fill_61`) × 2 axes: top-30 `cx`, `cy`, float32 `score`, `od` (rtol 1e-9), matched ann; `n_detections`, TP@10, TP@30 | B | 31,212 | 0 |
| **B4**: the production (D11) blank on 7 ROIs × 6 clicks × 5 arms × 2 axes, TP@10 against the notebook blank | B | 420 | 0 changed |

(Rows are grouped from `results/multiotsu_pad_arms_49roi_6seed_audit_ledger.csv`, which lists all 46 checks.)

**Gates.**

- **Execution coherence: passes.**
  - The 15 code cells carry counts 1–15, contiguous from 1.
  - There is no error output and no unrun cell, and the last cell is markdown.
  - The two PNGs are the two cell-27 figures. (`..._audit_execution_gate.csv`)
- **Provenance: clean for what this notebook computes.**
  - The modules (`template_match.py`, `production.py` and `seed_selection.py` 2026-09-17; the rest 09-14 to 09-16)
    are committed, unmodified, and predate everything.
  - The prompt (09-29 21:55) predates the smoke run (22:29), which predates the full run's artifacts (23:46) and the
    notebook's last save (23:51).
  - No `.partial` file exists.
  - The notebook's own artifacts are untracked. Every row of them is re-derived above, so the Tier A results are
    **independent reproductions**, not consistency with a stale file.
- **Provenance of the three reference arms**, which are read, not computed. The window-sweep CSVs are untracked, and
  their notebook was re-saved after they were written: 09-29 21:53, and again at 13:08 **today, during this audit**.
  The CSVs themselves are unchanged since 09-29 15:06 (SHA-256 prefixes in `..._audit_provenance.csv`).
  - `default_51` and `fill_51` are bit-identical, row for row, to the **committed** fill-geometry run's `default_51`
    and `hem_rect_fill` (`a9293d0`) on all 294 clicks × 2 axes, 9,408 values. So those two references have committed
    provenance.
  - `fill_61` has none. It is backed by four things: the notebook's 12-click parity gate; my SHA rebuild of its
    template on all 294 clicks; my own search reproducing its top 30 exactly on the 17 stress clicks (both axes); and
    the prior audit's §B2, which reports reproducing it on all 294.
- **Composition: passes, 22 items, 0 failures.** (`..._audit_composition.csv`)
  - `per_run` holds 2352 rows = 294 × 4 × 2, with no duplicate key and 294 rows per (arm, axis).
  - `top30` holds 70,560 rows, and every list ranks 0–29.
  - `templates` holds 1176 rows, all accepted.
  - The other outputs hold `per_roi` 686, `delta_per_roi` 2940, `delta_stats` 94 (60 + 34 reruns) and
    `selfmatch_sensitivity` 20 rows.
  - The references hold 1764 rows and 52,920 top-30 rows.
  - Clicks match the clicks CSV on all 2352 rows.
  - **7 domains × 7 ROIs**, and the domain label equals `MIDOG++.json`'s `tumor_type` on 49/49. The database has 44–150
    ROIs per tumour type, so the ROI is not collinear with the stratum.
- **Claims with no persisted artifact under them.** Three kinds, each handled another way:
  - **Pixel-level claims** (the boxes, the templates, "no self-match anywhere in the full list", the mass-weight
    values): Tier B B1–B3.
  - **The `n_selfmatch_all = 0` claim.** It lives only as a column, because full lists are not persisted. I checked it
    from the pixels on the 17 clicks at risk (the 14 `box_outside_blank` and the 3 non-positive-mass-weight clicks),
    all 6 arms, and found 0.
  - **The two figures** had both steps of the figure protocol.
- **Tier C: not run.** Nothing in Tier A or B diverged, so there was nothing to explain. No sweep was re-run through
  the notebook. The measurements that touch pixels are Tier B §B1–§B4 and the re-timing.

## Part 1 — findings

**No Tier 1 and no Tier 2 finding.** Nothing I measured changes a conclusion or a number the notebook reports. What
follows are Tier 3 notes: fragility a reader should know about, one scope omission that I measured as immaterial, and
wording.

### Tier 3

**1. The `chromatin_od` `mo_pad61` null turns on one ROI.**

The notebook reads `chromatin_od` `mo_pad61 − mo_box` = +1.97 [−0.31, +3.98], Holm 0.0876, as "no gain is established".
That is correct under the pre-registered rule, and the verdict stands. But leave-one-ROI-out
(`..._audit_loo_roi.csv`) shows how little separates it from the other verdict:

- `293.tiff` (canine lymphosarcoma, 104 mitoses) contributes **−22 TP** over its 6 clicks.
- Without that ROI the contrast is **+2.78 [+1.28, +4.38], 27/10/11, D5 met**.
- It is the only one of the 490 leave-one-out refits across the 10 tested rows that changes a verdict.
- At the domain level the same ROI shows up as canine lymphosarcoma's −19 TP. The other six domains are +3 to +29
  (`..._audit_domain_level.csv`).
- `293.tiff` is also the most influential ROI for `chromatin_od` `mo_pad51 − mo_box` (−12 TP), where it flips nothing.

A reader should not take this null as evidence that `mo_pad61` fails at `chromatin_od`. All three pads move
`chromatin_od` the same way (+1.97 to +3.27), and `mo_pad61` misses the bar on one lymphosarcoma ROI. Dropping an
inconvenient ROI is not a licence to report the refit, and I do not.

This **extends** the prior audit's finding 3, where `293.tiff` alone moved the window sweep's `chromatin_od` `fill_71`
loss across Holm. In two notebooks now, one ROI steers `chromatin_od` verdicts.

**2. `chromatin_od` `mo_pad51 − mo_box` is the thinnest of the five "met" verdicts.**

The arithmetic reproduces exactly: +2.76 [+0.95, +4.56], 26/8/15, Holm 0.0105. How thin it is:

- **The majority clause has one ROI of slack.** 26 of 49 move up, and D5 needs 25.
- **Leave-one-ROI-out keeps it met on all 49 refits.** The lowest lower bound is +0.66 points.
- **Clustered at the domain, it is 6/7 positive**, with the one negative domain at −6 TP. The exact p is 0.047 and the
  Holm p over the pad family 0.094.
- **That is a resolution limit, not contrary evidence.** At G = 7 with a family of 3, the lowest attainable
  Holm-adjusted p is 3/64 = 0.047, and only a unanimous 7/7 reaches it; a 6/7 has exact p 3/64 before any
  adjustment. So at the domain level nothing short of unanimity resolves. `chromatin_od` `mo_ring5` and all three `tm_score` pad contrasts are 7/7 and reach the floor.

"Met under the ROI as the unit; not resolvable at the domain level unless unanimous." Part 4 names the premise.

**3. The notebook blanks 51 px at the click for every arm, not D11's footprint, and never says so. Measured on 420 lists, this changes nothing.**

D11, as applied 2026-09-17 and kept 2026-09-22, blanks "a `base_size x base_size` square centred on the rounded
template centre". `find_and_suppress` does exactly that. So in production the arms here would be blanked as follows:

| arm | production (D11) blank |
|---|---|
| `mo_box` | `bs × bs` (23–55 px) at the box centre |
| `mo_pad51` | 51 at the box centre |
| `mo_pad61` | 61 at the box centre |
| `mo_ring5` | `bs + 10` at the box centre |
| `fill_61` | 61 at the click |

The notebook's title says "post-D11". Its standing caveat ("the blank stays 51 × 51 on the click") discusses only box
overlap and does not say the convention departs from D11. The window sweep's own cell 25 did say so; this notebook
dropped it.

Tier B §B4 measured it. I re-ran all five affected arms with D11's blank on the first ROI of each domain in the
fixed ROI order: `300`, `201`, `245`, `459`, `013`, `529` and `402`, which is 42 clicks and 210 searches.

- **TP@10 changed on 0 of 420 (search, axis) lists.** Every tested contrast on the subset is identical under both
  blanks (`..._audit_tierb_prodblank_contrasts.csv`).
- **The top 30 is identical to the notebook's on 420 of 420 lists**, and `n_detections` is equal on 210 of 210
  searches.
- **The blanks do differ by construction.** For example, `mo_box`'s is `bs²` = 529–3025 px at the box centre, against
  2601 at the click. What keeps the lists identical is the 100-peak cap:
  - on 193 of 210 searches, the largest response within 80 px of the click is below the 100th peak's score (median
    ratio 0.72), so nothing near the click can enter the pool under either blank;
  - on the other 17 (maximum ratio 1.35), the top 30 is still identical.
- This **extends** the prior audit's finding 4 from the window sweep's fill arms to the four multi-Otsu arms.

*Caveat.* 7 of 49 ROIs. The mechanism, the near-click response sitting below the `MAX_PEAKS` cut, holds on 193 of the
210 searches and is not specific to these ROIs, but the other 42 ROIs were not re-run.

**4. "Which pad to adopt is not answered" is right by the pre-registration. At the production axis, the untested ordering is not ambiguous.**

The pad-vs-pad rows are descriptive by design, and the notebook is right not to verdict them. "The axes disagree on
which pad is best" is also accurate: `mo_pad61` leads at `tm_score` and `mo_ring5` at `chromatin_od`, and the same pad
wins all 49 leave-one-ROI-out folds on each axis (`..._audit_loo_selection.csv`).

D5, though, makes `tm_score` the production ranker and `chromatin_od` "a reported second axis". At `tm_score` the
descriptive ordering has CIs clear of 0:

- `mo_pad61 − mo_pad51` +1.70 [+0.68, +2.76], exact p 0.0023;
- `mo_ring5 − mo_pad61` −3.16 [−4.83, −1.60].

So the open question is narrower than the prose makes it: at the production axis, `mo_pad61` is the untested favourite.

**5. Wording: "the same holds for the three references" (self-match) covers only their top 30.**

Cell 17 prints that whole-list counts are not saved for the references. The cell-28 sentence follows a clause about
"any arm's full detection list", so it reads as if the references' full lists were checked. On the 17 stress clicks,
my own full-list searches of `fill_51` and `fill_61` find 0 self-matches, so nothing is wrong. The sentence claims one
level more than the notebook checked.

**6. Timing is one contended run.**

- **The notebook's figures.** "64.2 min for 1176 searches" and the per-arm medians (2.48 / 2.90 / 3.29 / 2.63 s)
  reproduce from `t_search_s`. They are a single execution on a machine the prompt itself describes as shared with
  another long-running job, with no repeats.
- **My re-timing.** Three repeats on `300.tiff` s0, at a 1-min load average of 3.53: `mo_box` median 2.21 s (2.01–2.60)
  and `mo_pad61` median 3.36 s (2.60–3.50) (`..._audit_retiming.csv`). That is the same ordering, with a ratio of 1.52
  against the notebook's 1.33.
- **How to read it.** The notebook presents these as a run log, not a cost claim, so I do not score it higher. As a
  per-search cost of adopting a 61 px pad, it is an order of magnitude, not a measurement.

**7. The pre-registration's timing is not durable.**

Cell 0 says it was written before the full run. The notebook and prompt are untracked, so no commit dates either
one. The prompt's mtime (21:55) precedes the smoke run (22:29) by 34 minutes.

Every design element in cell 0 matches the prompt, and each **[gap]** resolution is at least as conservative:

- the two-ended worst case;
- the union neighbourhood;
- `pair_stats_clicks`, never triggered because there were 0 refusals.

I found no departure between what cell 0 promises and what the code does. That covers the families, the tested set,
D5, the flags, the rerun definition and the contamination rule.

### Checked and found sound (so a reader need not repeat it)

- **Unit of analysis (mode 1).** Every test is an exact ROI-level sign-flip on click-summed integer deltas. The
  effective units are `n_nonzero`, 34–47 in the pad family, with floors 2⁻³³ to 2⁻⁴⁶. The cluster bootstrap is used at
  G = 49, above its small-G regime. The t interval on 48 df makes the same exclude-0 call as the bootstrap on all 10
  tested rows; the largest gap between the two intervals' bounds is 0.25 points (`..._audit_loo_roi.csv`).
- **Cross-family multiplicity.** Holm over all 10 tests instead of four families leaves every verdict unchanged. The
  weakest met row, `chromatin_od` `mo_pad51 − mo_box`, goes to 0.032 (`..._audit_cross_family_holm.csv`).
- **Product relevance (mode 8).** At `tm_score`, padding helps the product metric as well as the mean:
  - The worst seed per ROI rises +9.80, +11.22 and +6.73 points for `mo_pad51`, `mo_pad61` and `mo_ring5` over the
    bare box, p ≤ 2.5e-4.
  - No click loses more than 2 TP@10.
  - The bare box is the product regression: 26 clicks lose ≥ 3 TP against `default_51`.
- **The vs-incumbent nulls hide no tail.** `mo_pad51 − fill_51` and `mo_pad61 − fill_61` reproduce as nulls bounded
  to about ±1 point, as stated. The per-click and worst-seed reads (`..._audit_product_metrics.csv`) agree:
  - no click loses more than 2 TP@10 against the incumbent, on either axis;
  - the worst-seed precision per ROI moves −1.02 points (3/8/38, p 0.23) and +0.41 at `tm_score`, and −0.61 and −0.20
    at `chromatin_od`. All four are null.

  No regression is hiding behind the pooled nulls. The references entered by matching side (51 and 61), not by
  outcome, so the vs-incumbent family carries no selection of its own.
- **Precision null (mode 4), for context only.** The length-matched random-list precision@10, pooled over the 294
  clicks, is 0.0049. Every arm is about 100× above it; the contrasts between arms do not need it.
- **Recall triad (mode 3): not applicable.** No recall is claimed, and every list delivers K = 10: the fewest
  detections is 83, and 78 after the edge band.
- **Ties.** 0 of 2352 lists have a key tie across ranks 9/10 (`..._audit_ties.csv`).
- **The mass-weight mechanism.** The claim that "the response at the template's own position is positive" holds
  exactly, not approximately. With the kept box lying on the flat blank, `R = bs²·(m_k − f)·(min − Ī)`, whatever the
  pad covers. None of the three flagged clicks is `box_outside_blank`, and my full-list searches find no self-match on
  any of them.

## Part 2 — verdict per conclusion

| # | cell | claim | verdict | tier |
|---|---|---|---|---|
| 1 | 1, 28 | configuration matches `midog_utils.production` | reproduces (read `production.py`: channel, `TM_CCOEFF`, 7, −1.5, 100, 51) | — |
| 2 | 8, 28 | the stats code reproduces the 10 fill-geometry `delta_stats` rows; `pair_stats_clicks` ≡ `pair_stats` | reproduces (my own code, 110 values) | — |
| 3 | 6, 28 | 294/294 accepted, 0 refusals; tripwire 19/19; `pad_px` = 2 | reproduces (B1: box re-implemented on all 294) | — |
| 4 | 13, 28 | parity gate 48/48; `od` bit-exact 1440/1440 | reproduces in substance: `default_51` and `fill_51` ≡ committed fill-geometry on all 294 × 2; `fill_61` top 30 reproduced from pixels on 17 clicks | — |
| 5 | 15, 28 | asserts 1–4 on 1176 searches; 6 `box_ge_S` lists ≡ `mo_box` | reproduces | — |
| 6 | 21, 28 | `fill_61 − fill_51` reproduces the sweep's saved row | reproduces | — |
| 7 | 6, 28 | `bs` 23–55 (median 33), ≥ 51 on 6, base 259 / fallback 35, offset median 2.50 (max 9.66) | reproduces | — |
| 8 | 17, 25, 28 | contamination rule did not fire; 0 self-matches in top 10 / top 30 / full list; near-seed 0 | reproduces (full list checked from pixels on the 17 at-risk clicks); the reference clause claims one level more than was checked (finding 5) | 3 |
| 9 | 17, 28 | flag table; minimum mass weights 0.00053 / −0.00503 / 0.00445; the three `mo_pad61` clicks and their values | reproduces (Tier A geometry; B2 from pixels) | — |
| 10 | 28 | the mass-weight mechanism: response at own position positive; no self-match followed | reproduces (closed form; B3) | — |
| 11 | 23, 28 | flag-excluded reruns leave every D5 verdict unchanged; only one sign changes; the quoted rerun numbers | reproduces | — |
| 12 | 21, 28 | the 10-test table | reproduces (60 printed values; 1,420 underlying) | — |
| 13 | 21, 28 | no row floor-blocked (max floor 1.91e-06); 5 CIs exclude 0 on 20 re-seeds; closest bound +0.85 | reproduces (the re-seed part is consistency) | — |
| 14 | 28 | **headline:** at `tm_score` every pad beats the bare box, +5.00 to +8.16, Holm ≤ 5.1e-07, D5 met; 36–38/49 up | **reproduces**; survives LOO, domain clustering (7/7, at the floor), cross-family Holm and the D11 blank on the 7-ROI subset | — |
| 15 | 28 | pooled 0.4602 → 0.5102 / 0.5248 / 0.5418 | reproduces | — |
| 16 | 28 | at `chromatin_od`, `mo_pad51` and `mo_ring5` meet D5 | reproduces; `mo_pad51` has one ROI of majority slack and is not resolvable at the domain level unless unanimous (finding 2) | 3 |
| 17 | 28 | at `chromatin_od`, no gain established for `mo_pad61` | reproduces; the null turns on one ROI, `293.tiff` (finding 1) | 3 |
| 18 | 28 | edge band dropped: `mo_pad61 − mo_box` at `chromatin_od` +2.69 [+0.85, +4.52] | reproduces (descriptive) | — |
| 19 | 28 | "the `chromatin_od` results do not replicate the `tm_score` ones" | reproduces: in context it restates the pre-registered rule that one axis is not a replication of the other | — |
| 20 | 28 | vs incumbent: no difference detected on either axis; CIs about ±1; no equivalence claimed | reproduces; no per-click or worst-seed tail behind the nulls (checked and found sound) | — |
| 21 | 28 | the vs-incumbent family moves three things at once; nothing attributed | reproduces (correct design reading) | — |
| 22 | 21, 28 | descriptive table against `default_51` (12 values) | reproduces | — |
| 23 | 28 | the bare box is worse than `default_51`; much of the gain is recovering it; padded, it lands where the fill template lands | reproduces (descriptive; −4.83 of +6.46 at `mo_pad51`) | — |
| 24 | 28 | pad-vs-pad rows; "the axes disagree on which pad is best" | reproduces; at D5's production axis the untested ordering favours `mo_pad61` (finding 4) | 3 |
| 25 | 17, 28 | `MAX_PEAKS` binds on every row; NMS suppression 0.918–1.72; edge-band shares; band leaves every sign | reproduces | — |
| 26 | 15, 27, 28 | runtime 64.2 / 70.0 min; per-arm medians | reproduces as a record of one contended run (finding 6) | 3 |
| 27 | 28 | "Which pad to adopt is not answered here" | reproduces by the pre-registration (finding 4) | 3 |
| 28 | 28 | scope: 294 clicks drawn under the old Otsu joint gate | reproduces; a premise, Part 4 | — |
| 29 | 28 | caveat: blank 51 on the click; 14 `box_outside_blank`, no self-match, excluding them changes nothing | reproduces; omits that the blank is not D11's, which measured as immaterial (finding 3) | 3 |
| F1 | 27 | Figure 1: every K = 10 contrast, shared x-axis, tested solid / descriptive hollow, W/L/T and Holm | reproduces: series re-derived (30 rows, `..._audit_figure1_series.csv`); render read. All 30 annotations match; span −7.35 to +10.65 is inside the axis, so nothing is clipped; bars are the bootstrap CI as labelled | — |
| F2 | 27 | Figure 2: the four templates for 8 clicks with click, blank and box marked | reproduces: picks re-derived by the cell's own rule (8/8); every title and row label (`S`, fill, mass wt, flags, `bs`, offset) matches `..._audit_figure2_captions.csv`; templates rebuilt from pixels (B2); the box/blank geometry drawn agrees with the flags (e.g. `128.tiff` s4 box crossing the blank at offset 9) | — |

"cannot check" was not needed on any claim.

## Part 3 — what was re-run versus read

- **Audit script.** `multiotsu_pad_arms_49roi_6seed_audit.py` ran **start to finish, clean, in 20.3 min** under
  `/Users/mohinianand/anaconda3/bin/python3`, and regenerated every `results/multiotsu_pad_arms_49roi_6seed_audit_*.csv`
  table cited here. Every number in this log is printed by it or is in one of its tables.
- **Tier A, from persisted artifacts.** Covers everything in the Part 0 table marked A:
  - my matcher against the JSON;
  - my statistics: integer DP sign-flip, Holm, t on 48 df, D5;
  - my flag formulas from the templates' geometry;
  - my diagnostics from the top 30.
  - It also covers what the notebook did not run:
    - leave-one-ROI-out on the 10 tests;
    - domain-level sign-flip and t on 6 df;
    - Holm across all 10;
    - per-click and worst-seed product metrics;
    - leave-one-ROI-out selection among the pads and among all 7 arms;
    - the length-matched precision null;
    - K = 10 boundary ties.
- **Tier B, from the pixels.** 17.4 min of compute, 312 searches. The machine showed a load average of 5–12
  (`uptime`) during the runs, and 3.53 at the re-timing.
  - **B1**: my multi-Otsu box on all 294 clicks, plus mpp from the TIFF tags.
  - **B2**: all 7 templates (4 arms, 3 references) on all 294 clicks, with their SHA-1, fill and `m_keep`.
  - **B3**: my search on the 17 stress clicks. These are the 14 `box_outside_blank`, where the kept box reaches
    unblanked image and a self-match was most possible, and the 3 non-positive-mass-weight `mo_pad61` clicks, where
    the response at the template's own position is positive. It covered 6 arms × 2 axes, including full-list
    self-match counts.
  - **B4**: the D11 blank on 7 ROIs × 6 clicks × 5 arms.
  - **Re-timing**: `mo_box` against `mo_pad61`, 3 repeats each.

  I chose B3 because a self-match there would have been the single most damaging defect, and B4 because it is the
  only way the headline could be a convention artifact.
- **Tier C: not run** (no divergence to explain). No sweep was re-run through the notebook.
- **Read, not re-run.** `production.py`, `find_and_suppress.py`, `template_match.py`, `evaluate.py`, `nms.py`,
  `chromatin.py`, `dataset.py`, `channels.py` and `invariants.py`, for the implementation review; and every code cell
  of the notebook.
- **Implementation review.**
  - **Config.** No config drift: the notebook asserts equality with `production.py`, and I read it.
  - **Caps.** `MAX_PEAKS` = 100 binds everywhere, by D9's design. K = 10 is always delivered.
  - **Leakage.** The click's own annotation is removed from GT on every list.
  - **Matching.** Greedy one-to-one against both categories (the house three-bucket rule), with the radius per image
    from mpp (7.5 µm, D7).
  - **Pairing.** The same 294 clicks go to every arm, so the design is paired.
  - **Numerics.** No NaN in `od`, and no boundary ties.
  - **Definitional drift.** `ds.image_annotations(annotations, fn)` is called with `category_id=None` *deliberately*,
    so look-alikes are matchable (bucket `human_rejected_label`). `n_gt_mitotic` filters by category. "Self-match"
    means non-TP, which includes a look-alike match; that is the prompt's definition and is applied uniformly.
- **Beyond the named modes, I looked for:**
  - whether the headline could be an artifact of the non-D11 blank (B4: no);
  - whether the 100-peak cap could hide a near-click self-match below the top 30 (B3 full lists: no);
  - whether the references came from a different pipeline state (they are bit-identical to the committed run for 2
    of 3, and pixel-reproduced for the third);
  - whether one ROI or one domain carries any verdict (finding 1: yes, for the `chromatin_od` `mo_pad61` null);
  - whether the pad gain costs a worst seed (it does not);
  - whether the smoke run and the full run differ (bit-identical);
  - whether a rank-10 tie decides any list (none does).
- **Prior audits.** I read `2026-09-29-window-sweep-51-61-71-49roi-6seed-audit.md` last, after my own findings.
  - **Reconciliation of my findings:**
    - Finding 1 **extends** its finding 3 (`293.tiff` steering `chromatin_od`).
    - Finding 3 **extends** its finding 4 (the 51 px blank departs from D11's wording, measured immaterial) from the
      fill arms to the multi-Otsu arms.
    - My `fill_61` checks (SHA on all 294 clicks; bit-exact top 30 on 17; D11-blank TP@10 unchanged on 42 clicks ×
      2 axes) **confirm** its §B2 and §B5 on the one reference with no committed provenance.
  - **Its findings I did not reach:**
    - finding 1 (Q2's self-match bound on `default_61` / `default_71`) and finding 2 (Q1's flag lift): *did not look*.
      They concern arms this notebook does not use.
    - Its Tier 3 figure notes: *did not look*.
- **Facts in the agent file's appendix that no longer hold** (`..._audit_appendix_drift.csv`):
  - **Notebook, PNG and CSV counts.** The recount glob finds **86** notebooks, not 59, and **442** embedded PNGs, not
    217. `results/` holds 602 CSVs, of which **56 are untracked** (the audit's own tables
    excluded). The appendix's "all 233 tracked" no longer holds.
  - **`FSConfig.tm_method`** now defaults to `cv2.TM_CCOEFF`, not `TM_CCOEFF_NORMED`.
  - **`FSConfig.nms_radius = None`** now means "caller must supply", and `find_and_suppress` raises `ValueError` when
    it is unresolved. It no longer resolves to the image's match radius.
  - **Still true:** `image_annotations`' `category_id=None` default and `MIDOG_RADIUS_UM = 7.5`.
  - **Timing: neither confirmed nor contradicted.** My re-timed single search, excluding the ROI load, took a median 2.21 s
    (`mo_box`) under load. That is blank, pad, `matchTemplate`, dilation peaks, NMS and `od`. The appendix's ~8 s
    single-seed pass is a different composition: it includes `load_roi` and a `peak_local_max` call this pipeline does
    not use. The ~87 s FCOS figure was not re-measured.

## Part 4 — premises this audit inherited

| premise | source | what would falsify it | verdicts that lean on it |
|---|---|---|---|
| the ROI is the exchangeable unit | D5; F5 §8 | Measured: 7 ROIs per domain, so not collinear. If ROIs within a domain are not exchangeable (shared scanner or stain), the domain is the unit. At G = 7 with a family of 3, the Holm floor is 3/64 = 0.047, and only a unanimous 7/7 reaches it. Every `tm_score` pad contrast and `chromatin_od` `mo_ring5` are 7/7 and reach it. `chromatin_od` `mo_pad51 − mo_box` is 6/7 (Holm 0.094): not resolvable at that G, rather than refuted. | 14, 16 (finding 2) |
| D5's bar (a CI excluding 0 **and** a majority of ROIs) as the verdict rule | `DECISIONS.md` D5, "What would change my mind" | a rule without the majority clause, or with a two-sided equivalence margin, would re-label the `chromatin_od` rows. `mo_pad51` has one ROI of slack. | 14, 16, 17, 20 |
| the NMS radius is the match radius, 7.5 µm | D7 | a different radius changing which detections survive into the top 10 | all |
| `MAX_PEAKS` = 100 | D9 | a larger cap letting near-click peaks into the list. That would make the blank convention matter: finding 3 rests on the cap binding above the near-click response on 193 of 210 searches. | finding 3 |
| the three-bucket nearest-free-GT matching (a detection closer to a look-alike than to a mitosis is not a TP) | `evaluate.py`, house convention | a two-category matcher changing TP@10 differentially across template sizes | all |
| precision@10 per click, and worst seed per ROI, as the product metrics | D4; the product framing | a product that reads to K = 20–30, where the notebook's descriptive rows (not verdicts) apply | 14–20 |
| the 294 joint-gate clicks represent a pathologist's click | the fill-geometry click draw; stated in the notebook | a click population the multi-Otsu method would face in production, not both old Otsu gates. The notebook states this. | all |
| the fixed 51 px click-centred blank is representative of production | `MULTIOTSU_PAD_ARMS_PROMPT.md`, against D11 | tested on 7/49 ROIs (finding 3): 0/420 changes | 14, 16, 20 |
| the prior audit's findings | `Research Logs/2026-09-29-...-audit.md`, same project hand under the same premises | not independent corroboration. I re-derived what I relied on (the `fill_61` provenance) rather than citing it. | none |

**ROI-per-stratum count measured:** 7 ROIs in each of 7 tumour types (49); the database holds 44–150 per type.
**Effective units:** `n_nonzero` is 34–47 ROIs across the pad family and 21–35 across the vs-incumbent family.
`211.tiff` is 0 TP@10 on all 7 arms on both axes (the notebook's `per_roi`, reproduced above), so it carries no
permutation entropy in any contrast.

**Which verdicts would change if a premise were wrong.**

- If the domain rather than the ROI were the exchangeable unit, `chromatin_od` `mo_pad51 − mo_box` would go from
  "met" to "not resolvable at G = 7" (6/7 domains, below the unanimity the floor requires). Nothing else would
  change: every other met verdict is 7/7.
- If production kept many more than 100 peaks, the blank convention would need re-measuring before the `tm_score` pad
  gain could be carried over to the production pipeline.
- The headline, padding's `tm_score` gain over the bare box, leans on no premise in this list that I could find
  shaky. It is 7/7 domains, 36–38/49 ROIs, every leave-one-out refit, and unchanged under D11's blank on the subset I
  re-ran. What it does lean on is the click population, which neither this notebook nor this audit can test; that
  question is for `premise-reviewer`.
