# Per-iteration precision vs. the precision that saves a pathologist time: agent prompt

**Input:** `NOTEBOOK`, a notebook that measured precision@K after some change on the 49-ROI MIDOG++ set. It must
have saved, for every click, its ranked candidate pools.

**Output:**
1. A new notebook that replays the click-review loop over several iterations with that notebook's best ranking
   method.
2. Per-page break-even precision targets, and whether the 49 ROIs reach them.
3. A short written verdict with an adversarial-check log.

**Treat this prompt as a hypothesis, not as truth.** Re-derive the cost model, check every number against the data,
and say plainly where the prompt or `NOTEBOOK` is wrong.

---

## 1. Background: the product loop being simulated

AnnotateDx is a click-to-verify tool. Each iteration works like this:
1. A pathologist draws a box on one mitosis.
2. The tool searches the ROI and shows the top K candidates.
3. The pathologist accepts the true ones and rejects the false ones.
4. Every shown candidate is hidden from later iterations.
5. The pathologist boxes another mitosis, and the loop repeats.

Two quantities matter at each iteration, or **page**:
- the page's precision;
- whether that precision **saves time** compared with finding the same mitoses by hand.

Precision@K after a single click, which is what `NOTEBOOK` reports, is only page 1.

## 2. The cost model

Re-derive it yourself. Check the units, the limits, and what changes with r.

**Per page:**

| symbol | meaning |
|---|---|
| `n` | candidates shown, which is K unless the pool ran out |
| `TP` | mitoses confirmed on the page |
| `p = TP / n` | the page's precision |
| `r` | ROI mitoses still unfound when the page is shown. The box just drawn counts as found. |
| `v` | seconds to accept or reject one candidate |
| `c` | seconds to draw the box |
| `L` | seconds of latency the user waits for |
| `T` | seconds for one complete manual pass over the ROI, assumed to find every mitosis |

**Regime A, the tool replaces search.** This is the primary regime.

Scanning a fraction f of the ROI finds f·r of the r remaining mitoses in expectation. The manual yield is therefore
`r / T` mitoses per second, and each mitosis costs `T / r` by hand. A page beats scanning when:

    TP / (n·v + c + L)  >  r / T
    ⇔  p  >  p* = r · (n·v + c + L) / (n · T)            (break-even precision)

    E = p / p* = TP · T / (r · (n·v + c + L))             (efficiency ratio; E > 1 saves time)
    seconds saved by the page = TP·T/r − (n·v + c + L) = (E − 1)·(n·v + c + L)

What follows from this:
- **p\* is proportional to r.** Sparse ROIs and late pages need *lower* precision to pay, because hand search is slowest
  there.
- **Targets:** E ≥ 1 is break-even. **E ≥ 2 is the target**, meaning the page is at least twice as fast as hand search.
- **Some targets are unreachable.** A page can be at most `p_max = min(1, r/n)` precise, so its best possible E is
  `E_max = p_max / p*`.
  - In dense ROIs `p*` is high: at the central setting with K = 10, `2·p* > 1` once r is above about 111. No method
    can reach E ≥ 2 there, because hand search is already fast.
  - Classify every page as unreachable at E ≥ 1, unreachable at E ≥ 2, or reachable.
  - **Report every reach rate twice:** over all pages, and over pages where the target was reachable. Judge the
    method on the second, so a dense ROI is not read as a method failure.

**Regime B, a full manual pass is mandatory anyway** (exhaustive annotation). Each tool-found mitosis then saves only
its decide-and-mark time inside that pass, `c_saved`:

    p*_B = (n·v + c + L) / (n · c_saved)

This does not depend on r. Report it, but it is secondary.

**Parameters.** None of these has been measured, so never report one setting alone.

| | values |
|---|---|
| `v` | 1.5, **2**, 3 s |
| `T` | 5, **10**, 20 min |
| `c` | 5 s |
| `L` | 2 s |
| `c_saved` | 4, 6, 10 s (regime B only) |

The central setting is v = 2 s and T = 10 min. Report every headline number at the central setting and across the
grid.

**Assumptions to state, and to challenge if you think they are wrong:**
- **T does not grow with ROI density.** If dense ROIs take longer to scan, their E is understated.
- **The pathologist judges perfectly.**
- **Mitoses are spread evenly, and the scan is systematic.**
- **Seed-finding is excluded from E.** A box drawn on a mitosis that is not yet found costs about `T / r` to find by
  hand and yields one mitosis at the manual rate. That is neutral, so it is left out of E. Report how often a later
  iteration's seed was already found by the tool: those boxes need no search.

## 3. Data

**First, find the best ranking method in `NOTEBOOK`.**
- Use the notebook's own pre-declared choice: its spec or selection rule.
- Then recompute pooled P@10/20/30 for every arm from its per-run CSV yourself.
- If the notebook's choice is not the highest-precision arm at K = 10, simulate both and say so.
- Always also simulate these comparison arms, when `NOTEBOOK` has them:
  - **the user's current ordering:** `chromatin_od` on the unfiltered pool (`od_unfiltered`);
  - **the D5 production default:** `tm_score` on the unfiltered pool (`tm_unfiltered`);
  - **the plain `od` arm** on `NOTEBOOK`'s filtered pool.
- **Comparing arms:** compare exact values. Ties are common at this granularity, and `idxmax` silently favours
  column order.
- **Do not re-run any search.** If the per-click pools are not saved, stop and report.

**Ground truth** comes from `databases/MIDOG++.json` via `ds.load_annotations` and `ds.image_annotations(annotations,
file_name)`. Category 1 is a mitosis and category 2 a look-alike; columns `ann_id`, `cx`, `cy`, `category_id`.
- **Pitfall:** per-run CSVs in this repo store `n_gt_mitotic` **excluding the clicked mitosis**. A ROI's total is the
  category-1 count from the JSON.
- **Pitfall:** the deep-pool `gt_rank` CSVs omit each click's own seed. Never use them as complete ground truth.

**Match radius:** `ev.radius_px(mpp)`, 7.5 µm (D7). The exclusion radius is the same value.

## 4. The simulation

For each ROI, take its recorded clicks (6 per ROI) and **20 random orders** of them, from
`np.random.default_rng([20261003, image_id])`. Use the same orders for every arm, so the comparison is paired.
Simulate K = 10, 20 and 30 as separate sessions, each with **6 iterations**.

Per order, start with `found = {}` (annotations confirmed), `claimed = {}` (annotations already matched to a shown
candidate, look-alikes included) and `shown = []` (shown candidate coordinates). Then for iteration `j` with click `s`:

1. **The box.** Add `s`'s own annotation to `found`. Note whether it was already in `found`.
2. **The pool.** Start from `s`'s saved pool: accepted candidates for filtered arms, every candidate for the
   unfiltered arms. Remove every candidate at distance ≤ the match radius from:
   - any coordinate in `shown`;
   - any annotation in `found` **other than `s`'s own**.

   `s`'s own search already blanked its seed, and `NOTEBOOK` scored it with the seed's annotation removed from ground
   truth. Keep that convention, so iteration 1 reproduces `NOTEBOOK` exactly.
3. **Rank** the remaining pool with the arm's method.
   - Single-key arms (`od`, `tm`) keep their order.
   - Rank-combination arms **recompute their ranks inside the remaining pool**, exactly as `NOTEBOOK`'s functions do,
     tie-break included.
   - **Sensitivity check:** also run the full-pool order with removed items skipped. Report the difference.
4. **The page.** Take the first `min(K, remaining)`.
   - If `r == 0`, the session ends.
   - If this click's pool is empty but `r > 0`, record a skipped page with `n = 0` and go on to the next click.
5. **Matching.** Match the page with `ev.bucket_detections(page, gt_eval, radius)`. Here `gt_eval` is the ROI's
   ground truth (both categories) minus every annotation in `found` (which includes `s`'s own) minus every
   annotation in `claimed`. At iteration 1 this is exactly ground truth minus the click, which is `NOTEBOOK`'s
   convention.
   - A page TP is a `HUMAN_CORRECT_LABEL` match.
   - Add every matched annotation to `claimed`, and the matched mitoses to `found`.
6. **Record** per page:
   - ROI, order, arm, K and iteration;
   - `n`, `TP`, `p`, and `r` before the page;
   - whether the seed was already found;
   - and, for every grid setting, `p*`, `E`, `E_max`, the reachability class, seconds saved and `p*_B`.

**Aggregation.**
- **Unit: the ROI.** Average over the 20 orders within a ROI first, then summarise across the 49 ROIs, so each ROI
  counts once.
- **Uncertainty:** ROI-cluster bootstrap CIs (`N_BOOT` 10 000, seed 20260917).
- **Density bins** use the ROI's total mitoses: 20 or fewer, 21-40, 41-80, more than 80. Also report by domain.

## 5. Gates: halt and report on any failure

1. **Pool reconstruction.** Each click's pool size and accepted count equal `NOTEBOOK`'s per-run values. Each ROI's
   clicks are distinct annotations.
2. **Ranking parity.** For every arm with stored positions, the reconstructed order of every click's pool equals
   `NOTEBOOK`'s stored positions exactly.
   - **Unfiltered arms:** they may have no stored positions; refused rows can be NaN. Rebuild them from the full pool
     **sorted by the score-order column first**, then a stable mergesort on the arm's key, NaN last.
   - **Their only check** is gate 3.
3. **TP parity.** Iteration 1, run on **every** click as a standalone page, gives TP@10/20/30 exactly equal to
   `NOTEBOOK`'s per-run values, for every simulated arm, unfiltered arms included.
4. **Ground-truth totals.** For every click, the JSON's category-1 count equals `n_gt_mitotic + 1`.
5. **Invariants.** These must hold throughout:
   - `found` only grows;
   - no annotation is counted twice;
   - `TP <= r`;
   - no candidate shown on a later page lies within the radius of one shown on an earlier page. Within one page, NMS
     already spaces candidates, so don't assert there;
   - `p*` and `E` are finite whenever a page is shown.
6. **Toy test.** Before the real run, check the simulator on a hand-built case of 2 clicks, 2 pools and 4 annotations
   with a known answer.

## 6. Outputs

**Per-page precision:** for each arm × K × iteration, by density bin and by domain. Mean over ROIs with CIs, plus the
share of zero-TP pages.

**Targets:** a table of `p*` and `2·p*` against r (5, 10, 20, 40, 80, 160) for K = 10/20/30, at the central setting
and across the grid. Mark any cell above 1 as unreachable. This is the answer to "what precision should we aim for".

**Reaching them?** At the central setting, for the best method and for the comparison arms:
- the share of pages with E ≥ 1 and with E ≥ 2, by iteration and density bin, over all pages **and** over reachable
  pages only;
- **page 1 per ROI**: the achieved p against `p*` and `2·p*`. List the ROIs below break-even, and those below target.

Then give the same headline shares across the grid, and the regime B shares.

**Time.** Seconds saved per ROI under each of these:
- all 6 pages reviewed;
- an oracle stop, at the first page with E < 1;
- an observable stop, after the first page with zero TP.

The pathologist cannot see r, so `p*` cannot drive a real stopping rule. A zero-TP page is something they do see.

Report minutes saved per ROI (median and IQR) and the share of ROIs with net savings. Then give the best method minus
production, as a paired ROI-level difference with a bootstrap CI.

**The two rankings within an iteration:** recomputed vs. full-pool order.

**Diagnostics:**
- under-filled pages (`n < K`) by iteration;
- the share of seeds already found;
- how late-page precision splits between depletion (r falling) and pool exhaustion.

**Files** (do not modify any existing file):
- the notebook `production_hematoxylin_only/multiotsu_pad_arms/iteration_breakeven_<NOTEBOOK stem>.ipynb`;
- `results/iteration_breakeven_<NOTEBOOK stem>_{pages,per_roi,targets,grid,summary}.csv`.

The last markdown cell holds the verdict and the adversarial-check log.

## 7. Adversarial checks: do these and log the outcome of each

1. **Re-derive the model.** Re-derive `p*` and E from scratch. Does `p*` really scale with r? Is anything
   double-counted, such as the box cost or seed finding? Is `(E − 1)·cost` really the seconds saved?
2. **Check `NOTEBOOK`'s own headline numbers** from its CSVs, e.g. its best-method-minus-`od` deltas. Note any
   mismatch.
3. **Selection bias.** Was the best method chosen on these same clicks? If so, every number for it is optimistic.
   Say so beside each headline.
4. **Why later pages fall.** Separate depletion (r falling), pool exhaustion (the saved pools hold about 100
   candidates) and genuinely worse ranking. In a product that hides shown candidates *before* the 100-peak cap,
   deeper candidates would refill the pool; the saved pools cannot model that.
5. **Matching across pages.** Can a mitosis or look-alike be claimed twice? Do the exclusions hide a mitosis the user
   never saw?
6. **Robustness of the verdict.** Does "reaching the target" survive the parameter grid, the ranking variant and the
   regime? If it hinges on one setting, say so.
7. **Sanity anchors.** Iteration-1 precision must equal `NOTEBOOK`'s P@K. Later pages should fall in sparse ROIs. If
   either fails, find out why before reporting anything.

## 8. Environment and conventions

- Interpreter: `/Users/mohinianand/anaconda3/bin/python3`. It is the only one that imports `cv2` and the repo's
  `midog_utils`. The notebook sits two levels below the repo root: `sys.path.insert(0, '../..')`.
- Execute with `/Users/mohinianand/anaconda3/bin/jupyter nbconvert --to notebook --execute --inplace
  --ExecutePreprocessor.timeout=14400 <path>`.
- No TIFFs are needed: everything comes from saved CSVs and the JSON. It should run in minutes.
- Docstrings: the opening `"""` on its own line, then a brief summary, one line per parameter (type and a short
  description) and one line for the return. Function definitions and calls go on a single line, never wrapped.
- Do not commit. Kill only a PID you started, never by pattern.

## 9. Your final message

Write these, in order:
1. **The answer** in 3-5 lines: what page precision should we aim for, and does the best method reach it on the
   49 ROIs, at the central setting and across the grid?
2. **The key tables:** targets, the E ≥ 1 / E ≥ 2 shares by iteration and density, and minutes saved.
3. **Caveats:** in-sample, unmeasured parameters, pool cap, regime.
4. **The adversarial-check log:** what you checked, what failed, and what you changed or flagged in this prompt or in
   `NOTEBOOK`.
