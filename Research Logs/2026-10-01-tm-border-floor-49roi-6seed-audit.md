# Audit of `tm_border_floor_49roi_6seed.ipynb`: every number and the A3 recommendation reproduce, from the CSVs and from the pixels; but the write-up's lead "contradiction", that the clip's effect is a reshuffle and not a targeted removal of white-adjacent false positives, is contradicted by the notebook's own M2 and rests on counting one-pixel peak shifts as changes of object

**Scope.** Target: `production_hematoxylin_only/multiotsu_pad_arms/tm_border_floor_49roi_6seed.ipynb`. It has 35 cells
(18 code), no embedded figures, and 3 montage PNGs saved beside it.

- **Writes:** `results/tm_border_floor_49roi_6seed_{per_run,top30,per_roi,delta_stats,floors,edge,white,nearseed,adversarial_log}.csv`,
  their `_smoke` twins, and `tm_border_floor_montage_{300,201,546}.png` [files since removed; the three montages are embedded in `multiotsu_pad_arms/tm_border_floor_49roi_6seed.ipynb`, cell 31]. It also pickles the search to
  `$TMPDIR/tm_border_floor_49roi_6seed_search.pkl`. That pickle is the only home of the full post-NMS pools, so this
  audit re-matched them and cached a copy (see Part 3).
- **Reads:** the fill-geometry `clicks` / `per_run` / `delta_stats`, the deep-pool `per_run` / `top1000`, and the
  pad-arms `templates` / `per_run` / `top30` CSVs, plus `databases/MIDOG++.json`. Every SHA-256 prefix matches the one
  cell 2 printed.
- **Documents read:**
  - the pre-registration (cell 0) and the prompt behind it, `TM_BORDER_AND_FLOOR_PROMPT.md` (untracked);
  - `DECISIONS.md` D3 and D5 (with its 2026-09-08 amendment), D9, and D11 (with both amendments);
  - for D1, the cross-reference only. It now points to `DECISIONS_UNVERIFIED.md`, and nothing here depends on it.
- **Prior audit, read last:** `Research Logs/2026-09-30-multiotsu-pad-arms-49roi-6seed-audit.md`, which covers the run
  `A0` reproduces. No earlier audit covers this notebook.
- **Audit script:** `tm_border_floor_49roi_6seed_audit.py` (repo root). **Tables:**
  `results/tm_border_floor_49roi_6seed_audit_*.csv`.

**Family: what else was run before this was written up.** I established it as follows:

- **Git history.** The notebook, its outputs and its prompt are all untracked, so `git log --diff-filter=A` has
  nothing to say.
- **Repo search.** `grep -rl` for `tm_border_floor`, `np.maximum(hem`, `floor clip` and `b_otsu` over every `.md`,
  `.py` and `.ipynb` finds only the prompt and this notebook.
- **Artifact stems.** The only `results/` stems are `tm_border_floor_49roi_6seed_*` and its `_smoke_*` twin.

The smoke run is the same code on 3 ROIs (`337`, `201`, `509`). It is not an arm. Its `delta_stats` plan (54 rows:
families, contrasts, K) and its 15 adversarial concerns and tests are identical to the full run's, which is the only
durable evidence that the families and the decision rule were not changed after it.

The arms (`A0`–`A3`, `A1v`/`A3v`, `A2m`) and the floor (`b_otsu`) were all fixed by the prompt before any run. The
headline `A3 − A0` is therefore the one pre-registered candidate, not a best of N. The hypothesis itself was formed on
these same 49 ROIs' `A0` lists: the prompt's motivating 24 % / 11 % and 82 % / 52 % white shares come from the
recorded `mo_pad61` run.

**Everything below is re-derived from `results/tm_border_floor_49roi_6seed_{per_run,top30,per_roi,delta_stats,floors,edge,white,nearseed,adversarial_log}.csv`,
from the search pickle's pools, from `databases/MIDOG++.json`, and from the TIFF pixels under `images/extra_valid/`.**

I wrote the following myself in plain `numpy` / `pandas` / `scipy` / `cv2`:

- the greedy matcher, and the ground-truth parse (TLBR box to centre);
- the exact sign-flip (a dictionary over attainable sums, a different algorithm from the notebook's convolution), Holm,
  and the t interval;
- the template construction, the D11 blank, the replicate pad, the `TM_CCOEFF` correlation, the robust floor, the
  15 px local maxima, the top 100 with tie-break, the band floor, the greedy distance NMS, and `od`;
- the M5 synthetic images.

Two repo dependencies remain, both stated: `midog_utils.dataset.load_roi` / `roi_mpp` (TIFF reading) and
`midog_utils.channels.to_channel` (colour deconvolution).

## Conflict of interest

- **Signals I can see.**
  - The notebook, its CSVs, its montages and its prompt are all **untracked** (`??`). `git log -1 -- <notebook>`
    returns nothing, so there is no author or age to attach. The notebook was last written 2026-10-01 19:34, about
    2.5 h before this audit; its CSVs at 19:28.
  - `git status --porcelain midog_utils/` is empty. The last commit touching it is `9e22eb9` (D11, 2026-09-17), which
    predates every artifact.
  - A user-memory note visible to this session (`tm-border-floor-fix-plan`) already summarises this run's result,
    including the mechanism reading this audit partly overturns ("gain does NOT track white"). The work and its audit
    are very likely the same hand, separated only by my fresh context.
- **What limits the conflict.**
  - Every table was re-derived from the raw artifacts, never from the notebook's aggregate outputs. Its aggregates were
    only compared against.
  - The matcher, the statistics, the search and the synthetic demonstrations were re-implemented from their
    definitions, not imported. None of the notebook's helpers is called.
- **What it cannot cover.** The design choices: the 61 px template, the 51 px blank at the click, `b_otsu` as the
  floor, the 20 % / gray-215 "white-adjacent" proxy, round-1 single-click P@K as the product proxy, and the ROI as the
  unit. Those go to `premise-reviewer` and to a reader who did not make them. Part 4 lists them.

## Part 0 — what reproduces

**The published statistics are the statistics the code and data in the working tree compute.** Provenance is clean:

- the inputs' hashes match;
- `midog_utils` is clean and older than every artifact;
- the pickle's `PER_RUN` equals `per_run.csv` on all 181,104 numeric cells.

So the Tier A results below are **independent reproductions**, not consistency with a stale artifact. On top of that,
my own re-implementation of the search reproduces the notebook's lists from the pixels.

| check | tier | values compared | divergent |
|---|---|---|---|
| execution counts contiguous 1..18, no error output, no unrun code cell | coherence | 18 | 0 |
| input SHA-256 prefixes vs cell 2's printout | coherence | 9 | 0 |
| `adversarial_log.csv` fields vs the executed cell-33 printout | coherence | 60 | **1** (Tier 3 finding 6) |
| composition: 4116 rows = 49 × 6 × 7 × 2; unique keys; 7 domains × 7 ROIs (domain = JSON `tumor_type`); top30 = 4116 × 30 | A | 8 | 0 |
| every top-30 bucket and `matched_ann_id`, re-matched from the raw JSON with my own matcher | A | 246,960 | 0 |
| top-30 rank order and `edge_px` (from JSON width / height) | A | 242,844 | 0 |
| `per_run` TP@10/20/30 and the top-10/30 white, edge and near-seed counts, from the re-match | A | 49,392 | 0 |
| pickle `PER_RUN` vs `per_run.csv` | A | 181,104 | 0 |
| full-pool buckets, re-matched from JSON | A | 203,441 | 0 |
| `per_run` pool columns (n_det, pool mitoses and look-alikes, edge 0 / ≤ 2 / < 30, near-seed), from the pools | A | 14,406 | 0 |
| `delta_stats`: pooled, W/L/T, exact p, n_nonzero, bootstrap CI, t CI, D5 bar, Holm (54 rows × 12) | A | 648 | 0 |
| `white.csv` (M2 shares and CIs, M3 counts / WLT / p / CI / Holm) | A | 132 | 0 |
| `nearseed.csv` (M4), `edge.csv` (M1), pooled P@K vs cell 23 | A | 129 | 0 |
| `per_roi.csv` (sums, SD, worst click, pool mitoses, white counts) | A | 6,860 | 0 |
| cell-34 per-domain table (clip, white, P@10 × 5 arms, A3 − A0, W/L/T) | A | 63 | 0 |
| numbers stated only in cell 34's prose or in cell printouts (incl. the pickle-only flows and band-leak counts) | A | 73 | 3 (wording; findings 4, 5, 10) |
| `floors.csv` re-derived from pixels, 49 ROIs | B | 294 | 0 |
| independent search, 10 ROIs × 6 clicks × 5 arms × 2 axes: template SHA-1, fill, threshold, n_det, top-30 xy / score / od, TP@10/20/30 | B | 4,800 | 0 |
| `white_frac_61` recomputed from RGB (top 30, A0 and A3, `tm_score`) | B | 3,600 | 0 |
| M5 synthetic demonstrations, re-derived | B | 13 | 0 |
| D11-footprint blank (61 px at the template centre) vs the notebook's 51 px at the click, A0 and A3, 3 ROIs | B | 108 | 0 |

The rows are grouped from `results/tm_border_floor_49roi_6seed_audit_checks.csv`. Over all checks, 954,962 values were
compared and 4 diverged: the hand-patched log row and three wording items. The churn counts in findings 1 and 2 match
the notebook exactly under its own definition, so they are not counted as divergences.

**Gates, and which applied.**
- **Execution coherence.** Applied; passes on counts. The source was edited after the run (Tier 3 finding 6).
- **Provenance.** Applied; clean.
- **Composition.** Applied; passes.
- **Figure protocol.** Applied to the 3 montage PNGs. They are saved files, not embedded outputs.
  - **Step 1b:** every bucket letter and every rank is checked against `top30.csv`.
  - **Step 2:** all three renders were read. Each matches the CSV on order, buckets and the prose's description of
    which tile moves. One small wording point is under Tier 3.

**Claims with no persisted artifact underneath them**, and what I did for each:
- **The pool-level flows** (2,520 / 10,310 / 2,960 / 9,790, and 631 / 571 / 604 / 507), and the band-leak counts
  (68 / 67). These live only in the temp pickle, which I re-matched and cached.
- **The scratch leak test** ("189–301 interior peaks per click on 18 clicks"). **Cannot check:** `leak_test.py` was
  never saved. Its in-run successor (68 / 294 clicks) does reproduce from the pickle.
- **The thumbnail descriptions of 337, 529 and 401** ("cytoplasm, collagen and retraction spaces"). **Cannot check:**
  the thumbnails were never saved. The white shares quoted beside them (1 %, 55 %) reproduce from pixels.

## Part 1 — findings

### Tier 1

**1. Cell 34 says the clip's effect "is not a targeted removal of white-adjacent false positives". The notebook's own pre-registered M2 says it is, and the churn counts offered against M2 measure one-pixel peak displacement, not candidates changing.**

The sentence overturned is cell 34's own (contradiction 1):

> "But the clip also reshuffles the top 10. 1193 of 2940 top-10 slots change: 560 mitoses leave and 669 enter. That
> is not a targeted removal of white-adjacent false positives. The ROI-level effect is not concentrated in white
> tissue."

**The notebook already answers this.** M2, its pre-registered and ROI-clustered test, shows the white-adjacent share of
top-10 FPs falling from 0.239 to 0.194, with CI [−0.064, −0.029] clear of zero. An untargeted reshuffle would leave
that share near 0.24. So "not a targeted removal" contradicts the notebook's own inferential result. The churn counts
are the only thing set against it, and they do not measure what the prose says.

**What the counts measure.** 1,193 / 2,940 slots, with 560 mitoses leaving and 669 entering, reproduce exactly, but
only under **exact pixel-coordinate equality**. The clip changes both the image and the template, so almost every
surviving peak moves by a pixel. Re-counted with the same object allowed to move (`..._audit_churn_top10.csv`,
A0 → A3, `tm_score`):

| definition of "same slot" | slots changed | mitoses leave / enter | FPs leave / enter | white-adjacent FPs leave / enter |
|---|---|---|---|---|
| exact coordinates (the notebook) | 1,193 (41 %) | 560 / 669 | 633 / 524 | 159 / 77 |
| ≤ 1 px | 484 (16 %) | 159 / 268 | 325 / 216 | **103 / 21** |
| ≤ 3 px | 449 | 134 / 243 | 315 / 206 | 100 / 18 |
| mitoses by annotation id | — | **128 / 237** | — | — |

Counted by object, the descriptive split agrees with M2:

- **The rates.** White-adjacent FPs are 23.9 % of A0's top-10 FPs (322 / 1,347). They are 31.7 % of the FPs that
  leave (103 / 325), but only 9.7 % of those that enter (21 / 216).
- **The net.** 82 of the net −109 FPs are white-adjacent, against 27 non-white. That is 75 % of the reduction, from 24 %
  of the population. Cell 34 itself reports the 82, then discounts it on the strength of the jittered churn.
- **By domain** (`..._audit_white_fp_decomposition.csv`). In 6 of 7 domains the clip removes 17–42 % of white-adjacent
  top-10 FPs, against −2 % to +11 % of the others. Melanoma's 17 % is 1 of 6.
- **Mast cell tumour is the exception, and the notebook is right about it.** Only 12 of its 122 white-adjacent FPs go
  (10 %), and its non-white FPs rise from 26 to 29.

**What stands.**
- The ROI-level statement is right: ROI pixel white share does not predict the per-ROI gain (ρ +0.08, p 0.60).
- A better covariate gives the same weak answer. With A0's own white-adjacent top-10 FP count, ρ = +0.22 (p 0.13) for
  A2 − A0 and +0.12 (p 0.43) for A3 − A0 (`..._audit_white_covariates.csv`).
- `b_otsu` is not a glass floor: it clips more than 60 % of pixels on 20 ROIs.

**The overturned part.** The conclusion that the clip's effect is not a targeted removal of white-adjacent FPs, and
the framing under "Read this first" that the clip helps largely for a reason other than the white surround, are both
wrong at the level the mechanism acts on. At the detection level the clip does what the white-surround mechanism
predicts. That the gain fails to line up with ROI white share says only that a ROI-wide pixel share is a poor proxy
for which candidates sit beside glass. A cytoplasm and stroma contribution, the notebook's own interpretation, may
coexist with this, but it is not needed to explain where the false positives went.

**This does not touch the decision.** The A3 recommendation rests on the precision contrast, M1 and M2, and all of
them reproduce.

**Caveats on my own finding.**
- It leans on the notebook's white proxy (more than 20 % of a 61 px window at gray > 215). M2 uses the same proxy.
- The inferential weight is M2's: its ROI-cluster CI. The leave/enter split is a summed net decomposition with no
  inference of its own, and is not a causal attribution per candidate.
- The ≤ 1 px tolerance is my choice. The pattern is the same at 1, 2, 3, 5 and 10 px.

### Tier 2

**2. The same coordinate-equality definition drives the pool and `chromatin_od` story, and that conclusion does not survive it.**

Cell 34 says the clip "changes the pool set a lot": about 2,520 mitoses leave A0's pool and 2,960 enter A2's. It then
says "the extra pool mitoses therefore do not reach od's top 10 on net". Both rest on the same definition. The numbers
reproduce under exact coordinates: 2,524 / 10,309 / 2,962 / 9,794, and 632 / 570 / 604 / 507 (±1 on the last).
Re-counted by object (`..._audit_churn_pool.csv`, `..._audit_od_top10_flows.csv`):

- **Pool.** By annotation id, **168** mitoses leave A0's pool and **604** enter A2's. The net is +436, as stated. The
  gross churn is 15× smaller than "about 2,520".
- **The `od` top 10.** Counted by position with jitter allowed, new-to-pool mitoses reach A2's `od` top 10 net
  **+45**. That holds at ≤ 1 px (111 in, 66 out) and at ≤ 3 px (66 in, 21 out), and it does not depend on matching order.
  - The id count corroborates it: **47** mitoses whose annotation is absent from A0's pool are in A2's `od` top 10, and
    **2** in reverse. These ids compare `od`-order matches against `tm`-order pool matches, so small order effects are
    possible.
  - So the new mitoses **do** reach `od`'s top 10.
  - What sinks `od` (A2 − A0 = −10 TP at K = 10) is newly pooled high-`od` non-mitoses: net +84 at ≤ 1 px, +85 at
    ≤ 3 px. They displace about 55 incumbent mitoses.
- "This is not simply '`od` sees the pool set, not its order'" therefore has no support left. The measured flows are
  exactly a pool-set effect: the clip adds high-`od` non-mitoses as well as mitoses.

The precision numbers are unaffected. The tier is 2 because it changes stated numbers and a secondary inference.

**3. The headline bundles a benchmark-only component. The product-relevant clip effect is +3.20 to +3.30 points, and the band's `od` gain is benchmark-only too.**

The notebook's own premise (from the prompt) is that "the product has no ROI border", so fix 1 is a benchmark
correction. Under that premise:

- **The headline.** `A3 − A0` = +3.71 includes the band's +0.41 at `tm_score`. The clip's estimate for the product is
  `A2 − A0` = +3.20 (31/7/11, Holm 2.6e-5) or `A3 − A1` = +3.30 (32/6/11, Holm 8.2e-6). Both meet the D5 bar.
- **The band on `od`.** "The band earns its place on the second axis" (`od` `A1 − A0` = +2.28) is a benchmark
  statement as well. The `od` window (51 px, 25 px replicate pad) reads padded pixels at the ROI edge exactly as the
  search does. In a product that reads real pixels past the edge, that gain has nothing to remove.
- **The contradiction.** Cell 34's "its value in the product is only what the `od` axis gains here" contradicts its own
  premise.

The decision rule is pre-registered on `A3`, and this changes no verdict under it. A DECISIONS entry should quote the
clip estimate as the product effect, and should present the band as a benchmark-harness correction. That correction is
worth keeping for every future number on these ROIs: it removes 174 band entries, 76 % unlabelled, from A0's `od` top
10s.

### Tier 3

**4. "The four ROIs whose worst click finds 0 stay at four" hides a swap.**

- A0's four zero-worst-click ROIs are `211`, `242`, `462` and `524`. A3's are `211`, `242`, **`245`** and `524`.
- `462` is rescued; `245.tiff` click 1 goes from 1 TP@10 to 0.
- Per click, A3 − A0 is 106 better, **21 worse (every one by exactly 1 TP)** and 167 the same.
- The worst click per ROI improves on 22 ROIs and worsens on 2 (`245`, `499`).

The product reading (worst-seed) is still favourable. It is just not "unchanged at four".

**5. "`A2m` is indistinguishable from `A2` on every K and axis (all Holm p ≥ 0.21)".** Every Holm p in the
`median_floor` family is 1.00. The smallest *raw* p is 0.18 (`od`, K = 20), and 0.21 is the raw `tm` K = 20 value. The
conclusion stands; the number quoted is not a Holm p.

**6. The executed notebook and its log disagree on one row.**

The cell-19 source and `adversarial_log.csv` (row `m1_band_count_disagreement`) were edited by hand after the run, at
19:34 against 19:28. As a result:

- the executed cell-33 printout still says "the prompt's 273 matches neither exactly";
- the CSV and the cell-19 source say "its 273 is `edge_px` ≤ 2".

Cell 34 discloses the patch. All of 252, 273 and 264 reproduce (`edge_px` == 0, ≤ 2 and < 2 in the A0 pools). Counts
are contiguous and nothing else differs (59 of 60 fields match the printout). This is a coherence note, not a numbers
defect.

**7. The per-domain paragraph mixes contrasts.** "The largest gains are in lymphosarcoma (+6.9) and neuroendocrine
(+5.7) … smallest in melanoma (+1.7) and mast cell (+2.1)" quotes `A3 − A0`. The Spearman it illustrates is `A2 − A0`.
Under `A2 − A0` the ordering is lymphosarcoma +5.95, neuroendocrine +4.52 and lung +4.29, down to melanoma +0.95 and
mast cell +1.67. The narrative survives.

**8. D3's wording is stronger than "does not cover a bottom clip".**

- D3's "what would change my mind" clause admits only transforms that are "monotone and **do not clip**".
- The `hematoxylin_od` channel already carries a bottom clip at 0: skimage's `rgb2hed` clamps negatives. A median of
  3.9 % of pixels per ROI is exactly 0, with a range of 0–25.6 % over the 49 ROIs (`..._audit_tierb_floors.csv`,
  `zero_share`).

The notebook's conclusion, that adoption needs a dated amendment to D3, is the right one either way.

**9. Decision condition 3 (M1) cannot fail for the band arms.** Gate 5 asserts that no band-arm detection has
`edge_px` < 30, so `edge_px` == 0 is 0 by construction. A0's 252 is a reproduction check. The rule therefore had three
discriminating conditions, and the border mechanism's own demonstration, M5(b), failed (0.87×). The notebook reports
that failure plainly.

**10. Small wording.**
- "83–100 detections per click" is A0's range; A3's is 79–100. Immaterial at K ≤ 30.
- On `201.tiff`, "A1 reorders the same set" is not quite right. A1 drops A0's band entry (x = 26, unlabelled) and adds a
  new mitosis at rank 10.
- "Under the clip the false positives are … not white-adjacent ones": one of A2/A3's six FPs there (rank 7, white
  0.26) is white-adjacent, against two of six in A0.

### Checked and found sound (so a reader need not repeat it)

- **Unit of analysis (mode 1).** Every test is an exact ROI-level sign-flip on click-summed integer deltas.
  - The primary family has 41–45 nonzero ROIs, so its floor is at most 2⁻⁴⁰. The bootstrap runs at G = 49, above its
    small-G regime.
  - The t interval on 48 df agrees: `A3 − A0` K = 10 gives [+2.37, +5.04] against the bootstrap's [+2.48, +5.07].
  - **At the domain level** (7 units, 7 ROIs each, so not collinear with the ROI), `A3 − A0` at `tm_score` is positive in
    **7 of 7** domains at every K. The exact p is 1/64 = 0.0156, and Holm over the 3 K gives 0.047.
  - That is the attainable floor. The headline survives a domain-level null, but only by unanimity, and nothing short of
    it resolves at G = 7.
  - `od` `A3 − A0` is 5/1 (p 0.0625) at the domain level.
- **Pre-registration adherence (mode 6).**
  - The code does what cell 0 promises.
  - Cell 0 matches the prompt, except for the declared fix-1 change. That change is justified: the band leak is real on
    68 / 294 clicks and worth at most +0.09 pt.
  - The prompt is internally inconsistent: M3 is "not a stopping rule" in §3, yet "recommend A1 alone" in its decision
    rule. Cell 0 resolved it toward the decision rule, and the point is moot because M3 shows no loss.
  - The post-smoke M5 rows are labelled descriptive and change no pass/fail line.
- **Decision-table placement (mode 7).** Under its own rule the result lands in "recommend A3": c1 met, c2 no loss on 6
  tests, c3 holds, c4 CI [−0.064, −0.029]. I agree.
- **In-sample selection (mode 10).**
  - Nothing was swept. A3 and `b_otsu` were fixed in advance, and `A2m` is a sensitivity. Leave-one-ROI-out has no
    operating point to choose, so it does not apply.
  - The real limits are already stated in the notebook: one optimistic 49-ROI draw (D5), round 1 only, and a hypothesis
    formed on the same ROIs.
- **The independence of the search.** My search, written from the definitions, reproduces the notebook's lists from the
  pixels on 10 ROIs × 6 clicks × 5 arms × 2 axes:
  - template SHA-1, threshold, n_det, top-30 coordinates, float32 scores and `od`;
  - TP@10/20/30 on all 600 lists.

  The arms differ only by the fix, as gate 3 claims.
- **The D11 blank convention.** The notebook blanks 51 px at the click. D11's text blanks the template footprint at its
  centre, which is 61 px for `mo_pad61`. The prior audit measured this as immaterial for the pad arms. I extended that
  check to the clipped arm: on `300`, `201` and `291`, all 6 clicks, A0 and A3, the 61 px footprint blank changes
  TP@10/20/30 on 0 of 108 values, and the top 30 is identical on 36 of 36 lists.
- **Not applicable, with the reason.**
  - Recall triad (mode 3): no recall is claimed, and n_det ≥ 79 > K.
  - Length-matched precision null (mode 4): no claim is made against random.
  - Treatment confounded with domain (mode 2): every arm runs on every click.
  - Leave-one-ROI-out (mode 10): see above.
  - Timing (mode 11): the 53 min and 6.9 s per click are a run log, not a cost claim. 6.9 s reproduces from
    `t_click_s`; my re-timing is in Part 3.
- **Ties.** 4 adjacent exact score ties across all 4116 top-30 lists, none straddling K = 10 or 20.

## Part 2 — verdict per conclusion

F*n* in this table refers to finding *n* in Part 1.

| # | cell | conclusion | verdict | tier |
|---|---|---|---|---|
| 1 | 34, 23 | Decision: recommend `A3` for a DECISIONS entry; all four conditions hold; M3 no loss | **reproduces** (c3 is non-discriminating; F9) | — |
| 2 | 34, 23 | `A3 − A0` tm K=10 +3.71 (0.542 → 0.579), CI [+2.48, +5.07], 34/7/8, Holm 2.2e-7, D5 met | **reproduces**; also survives a domain-level null at its floor (7/7, Holm 0.047) | — |
| 3 | 34 | No significant loss at any K on either axis; tm +3.71/+2.99/+2.54, od +1.80/+0.97/+0.78 | **reproduces** | — |
| 4 | 34, 19 | M1 holds: 0 boundary detections on band arms; A0 = 252; `edge_px` < 30: 766 | **reproduces** (true by construction for band arms) | 3 |
| 5 | 34, 21 | M2 falls: 0.239 → 0.194, CI [−0.064, −0.029]; od flat | **reproduces** | — |
| 6 | 34 | Almost all of the gain is the clip (`A2 − A0` +3.20); band alone +0.41 n.s. | **reproduces** | — |
| 7 | 34 | The band earns its place on od (+2.28, Holm 3e-4) and M1; 174 od top-10 band entries, 76 % unlabelled | **reproduces, overstated**: benchmark-only by the notebook's own premise (F3) | 2 |
| 8 | 34 | Adoption needs a dated D3 amendment; D3 "does not cover a bottom clip" | **reproduces, overstated** wording (F8); conclusion right | 3 |
| 9 | 34 | Contradiction 1: the clip's gain does not track white | **reproduces** at ROI level (ρ +0.08) | — |
| 10 | 34 | … 82 of the −109 net FPs are white-adjacent | **reproduces** | — |
| 11 | 34 | … 1193/2940 slots change, 560 mitoses leave / 669 enter: "not a targeted removal of white-adjacent FPs" | **does not reproduce**. The counts reproduce under exact-pixel equality, but the quantity the prose names (mitoses leaving and entering the top 10) is 128 / 237, not 560 / 669: the count measures peak displacement. The inference contradicts the notebook's own M2 CI (F1) | **1** |
| 12 | 34, 15, 27 | Contradiction 2: the floor is not a glass floor (156/294 ≥ fill, median 1.03, 0.47–3.33; > 60 % on 20/49; fill ×1.30, max 3.40; mass weight ×0.89) | **reproduces** | — |
| 13 | 34 | Thumbnail descriptions of 337, 529, 401 | **cannot check**: thumbnails never saved; white shares reproduce | — |
| 14 | 34 | `A2m` indistinguishable from `A2` (all Holm p ≥ 0.21) | **reproduces**, number mislabelled (F5) | 3 |
| 15 | 34, 13 | Contradiction 3: M5(b) fails, 0.87×; constant-pad comparators 1.41× / 2.04× | **reproduces** (Tier B) | — |
| 16 | 34, 13, 21 | Contradiction 4: M4 untestable (2–3 near-seed per pool, 0 in any top 30); M5(c) 1.181 → 1.328, 1.173 | **reproduces** | — |
| 17 | 34, 13 | M5(a) PASS is trivial (t 0.0253 < b 0.0272); above-floor tissue 27 % / 36 % vs predicted 48 % | **reproduces** (Tier B) | — |
| 18 | 34 | M2 / M3 / M4 per-arm table (shares, CIs, 362 → 354, Holm 0.70, CI [−24, +4]) | **reproduces** | — |
| 19 | 34, 23 | Primary and secondary table (all 54 delta_stats rows) | **reproduces** (648 values) | — |
| 20 | 34, 27 | The two fixes add (gap +0.10); A3 − A1 +3.30 ≈ A2 − A0 +3.20 | **reproduces** | — |
| 21 | 34 | Clip and band split by axis (A2 − A0 od −0.34) | **reproduces** | — |
| 22 | 34 | "Not simply od-sees-pool-set": 2,520/10,310 leave, 2,960/9,790 enter; od top-10 631/571/604/507; "extra pool mitoses do not reach od's top 10 on net" | **does not reproduce**. The counts reproduce under exact-pixel equality, but the quantity the prose names (mitoses leaving and entering the pool) is 168 / 604, not about 2,520 / 2,960. New-to-pool mitoses do reach od's top 10, net +45 (F2) | 2 |
| 23 | 34, 25 | The generator grows: 24.67 → 26.21 per click, 44/4/1, p 4e-10; A2 26.15; A1 24.68 | **reproduces** | — |
| 24 | 34, 25 | Robustness: 291 +13 (12 %); +3.33 without it; SD 1.03 → 0.99; worst click 4.06 → 4.55; four zero-worst ROIs stay four | **reproduces**; "stay four" hides a swap (F4) | 3 |
| 25 | 34, 25 | Per-domain table; A3 − A0 positive in every domain; mast cell 82 %; od −0.5..+6.0 mostly band | **reproduces** (63 cells); mixed contrasts in the prose (F7) | 3 |
| 26 | 34, 31 | Montages (300, 201, 546): composition and moves as described | **reproduces** (CSV and render); two wording points (F10) | 3 |
| 27 | 34, 33 | Adversarial log (band leak 68/294 (67), ≤ 0.09 pt; thresholds; MAX_PEAKS binds; CV 0.14–0.44; robust-stat ratios; hem.min() 0; 28/2782; 17 band mitoses in 75 slots; self-matches 0; 252/273) | **reproduces**; one row hand-patched after the run (F6) | 3 |
| 28 | 34 | Scratch leak test: 189–301 interior peaks per click on 18 clicks | **cannot check**: `leak_test.py` not saved | — |
| 29 | 34 | A1v / A3v within 0.1 pt everywhere; the decision does not depend on the band variant | **reproduces** | — |
| 30 | 34 | Fix 1's product value is "only what the od axis gains here" | **reproduces, overstated**. Supported: by the notebook's own premise, fix 1 has no product value on either axis; it is a benchmark correction (F3) | 2 |
| 31 | 2, 10, 11, 17, 29 | 28/2782 band mitoses; stats validation; gates 3–6; sharing check; hand re-match | **reproduces** (and the re-match extended to all 123,480 rows) | — |
| 32 | 34 | Timing: 53 min, 6.9 s per click | **reproduces** as a single-run log | — |

Two "cannot check" verdicts (13, 28). Both are on side material whose evidence was never saved: a scratch script and
eyeballed thumbnails. Neither is a pattern.

## Part 3 — what was re-run versus read

**Highest tier reached: B, with full Tier A.** No sweep was re-run, so this is not Tier C. Nothing in Tier A or B
diverged in a way that reading could not explain.

**`tm_border_floor_49roi_6seed_audit.py` runs start to finish** under `/Users/mohinianand/anaconda3/bin/python3` and
regenerates every `..._audit_*.csv`. The reported run took 11.7 min, exit 0, with 954,962 values compared and 4
divergent.

**Tier A** (sections 1–8 of the script). Read and recomputed:
- every per-run, per-ROI, delta-stats, M1–M4, domain and printed number;
- all 123,480 top-30 rows, re-matched from the raw JSON;
- the 203,441 pool rows from the search pickle, also re-matched.

**The pickle is not durable.** It sits in the system temp directory and will be cleaned. On the audit run it existed
and matched `per_run.csv` exactly. The script caches what it derives from it:
- the re-matched pool as `results/tm_border_floor_49roi_6seed_audit_pool_from_pickle.csv` (about 22 MB);
- the band-leak table as `..._audit_leak_from_pickle.csv`.

The script falls back to those copies when the pickle is gone. A later re-run is then consistency with this audit's
copy, not independent of it.

**Tier B, about 11 min of compute.**
- **Floors on all 49 ROIs.** `b_otsu`, `b_med`, clip shares, white share, `hem.min()`, and the share of exact-zero
  pixels, from the TIFFs.
- **An independent search on 10 ROIs × 6 clicks × 5 arms × 2 axes.** The ROIs were chosen by damage:
  - `291`, the largest single-ROI gain;
  - `401` and `339`, the heaviest clips (82 %, 76 %);
  - `300`, the M5 and montage ROI;
  - `201`, the most boundary detections, and a montage;
  - `546`, the stuck melanoma montage;
  - `400`, a +7 ROI;
  - `499` and `100`, the two −2 ROIs;
  - `462`, an edge-heavy ROI.
- **`white_frac_61` from RGB** on the A0 and A3 top 30s of those 60 clicks.
- **M5 (a), (b), (b2/b3), (c)** with my own template and correlation.
- **The D11-footprint blank** on 3 ROIs.
- **Re-timing.** On this machine at a load average of about 4, one 61 px `TM_CCOEFF` correlation over a full ROI has a
  median of 1.19 s, and one robust floor + peaks + NMS pass 0.26 s, over 180 correlations and 300 passes
  (`..._audit_tierb_timing.csv`). Three correlations and seven passes come to about 5.4 s, against the notebook's
  6.9 s per click. The notebook's figure also includes template building and `od` scoring, so the two are
  consistent.

**What I looked for beyond the named modes (mode 12).**
- **What "changed" means.** This produced findings 1 and 2.
- **Whether the clip effect is one ROI or one domain.** It is neither: 7/7 domains, and +3.33 without `291`.
- **Whether matching or ties could manufacture a gain.** They do not: 0 / 123,480 bucket divergences, and no tie at
  K = 10 or 20.
- **Whether the 51 px blank convention interacts with the clipped blank fill.** It does not: 0 / 108.
- **Whether the od window, like the search, reads padded pixels.** It does, which is the basis of finding 3.
- **Whether the per-click distribution hides a regression.** It does not: no click loses more than 1 TP.
- **Whether the channel was already bottom-clipped.** It was, at 0.

**Appendix facts that no longer hold, or have drifted.**
- `FSConfig.tm_method` now defaults to `cv2.TM_CCOEFF`. The appendix says it still defaulted to `TM_CCOEFF_NORMED`.
- `FSConfig.nms_radius = None` now means "the caller must supply it" and raises if neither is set. It no longer means
  "use this image's match radius".
- The appendix's "~8 s per single-seed pass" (51 px, `peak_local_max`) is pessimistic for this path. A 61 px
  `TM_CCOEFF` correlation measured 1.19 s and a `cv2.dilate` peak pass 0.26 s, so one pass is about 1.5 s plus the
  ROI load. The appendix's `peak_local_max` at 3.84 s is the difference.
- The `image_annotations` signature (`category_id=None`) and the 7.5 µm radius still hold. Here the both-categories
  default is the intended D6 bucketing, not a defect: `TP` is `bucket == human_correct_label`.

**Prior audit reconciliation** (`2026-09-30-multiotsu-pad-arms-49roi-6seed-audit.md`, read after my own analysis).
- **Its finding 3**, the 51 px blank at the click rather than D11's footprint, measured as immaterial: **extended**. It
  is immaterial for the clipped arm too, 0 / 108.
- **Its domain-level resolution note** (Holm floor 3/64 at G = 7 with 3 tests): **confirmed** and applied here. The
  primary sits exactly on it.
- **Its findings 1, 2 and 4** (`od` nulls turning on `293.tiff`, `mo_pad51` slack, pad choice): **did not look**. They
  concern contrasts this notebook does not make.
- **Its timing note (finding 6)**: **agree** in kind. This notebook's timings are also a single contended run.
- None of my findings overturns it. Nothing there covered the churn definition, which is new here.

## Part 4 — premises this audit inherited

ROIs per tumour type in the 49, measured from `databases/MIDOG++.json`: **7 in each of the 7 domains**. The ROI and the
stratum are therefore **not** collinear, and the ROI-level and domain-level tests are two checks, not one. There are 49
units, and 41–45 carry signal in the primary family. At the domain level all 7 carry signal.

| premise | source | what would falsify it |
|---|---|---|
| The ROI is the exchangeable unit | D5, F5 §8 | ROIs within a domain so correlated that the domain is the real unit. Partially tested: the primary survives at the domain floor (7/7) |
| Matching at 7.5 µm, greedy, both categories bucketed, the seed dropped | D7, D6 | an arm-asymmetric change in TP under a different radius, or a different ground-truth parse |
| The D5 bar (bootstrap CI excluding 0 and a majority of ROIs) is the adoption bar | D5 | a held-out ROI set where `A3 − A0` fails it. None was run, and held-out ROIs are not in this sample |
| Unlabelled detections are false positives | MIDOG++ exhaustive-labelling assumption | a pathologist review finding the clip surfaces more unlabelled true mitoses than A0 (or fewer) |
| Round-1 single-click P@K is the product proxy | product memo, D4 / D5 | later rounds (depleted pools, blanked prior candidates) ranking A3 and A0 differently |
| The product has no ROI border | `TM_BORDER_AND_FLOOR_PROMPT.md` §1, a user statement | the product running per frame without context beyond its edge. Then the band fix and its `od` gain become product-relevant (finding 3 reverses) |
| "White-adjacent" = more than 20 % of a 61 px window at gray > 215 | the prompt, the notebook | pale cytoplasm or stroma above gray 215 misclassified as glass, or lumen below it. Finding 1 leans on this proxy |
| An object that moves ≤ 1 px is the same candidate | this audit | a 1 px shift changing which structure a candidate sits on. At 30 px match radius and 15 px peak spacing, implausible |
| The prior audits are corroboration | `Research Logs/` | they are by the same process under the same premises, not independent |

**Which verdicts would change if a premise were wrong.**
- If the "white-adjacent" proxy is miscalibrated, or the ≤ 1 px identity is wrong, finding 1 (Tier 1) weakens back
  toward the notebook's reading.
- If the product does have an ROI border, finding 3 reverses: the band becomes product-relevant.
- If the ROI is not exchangeable, or round 1 is not the product, or unlabelled detections are often real, the decision
  verdict (#1, #2) is what moves. Every arithmetic verdict stays as it is. The adoption bar is the premise I would
  send to `premise-reviewer` first. A 49-ROI in-sample draw that D5 itself calls optimistic is the only evidence
  behind a production change. `MIDOG++.json` lists 37–128 unused annotated ROIs per tumour type; whether they are
  downloadable and pass the `extra_valid` filters was not checked here.
