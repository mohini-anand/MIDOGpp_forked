# Which `cv2.TM_*` variant, on the hematoxylin channel?

**Status: complete, and independently audited.** All three experiments run; **4257 invariant
checks, 0 failed** -- and then three independent audits found a further **fifteen** defects the
gates did not catch, listed in section 8. One invalidated an entire method
(`TM_SQDIFF_NORMED`, now withdrawn); one showed the invariant meant to catch it was guarding
the wrong quantity; one invalidated a claim this log had presented as its cleanest sanity
check; and one violated a rule this log pre-committed to in its own section 4. All are
corrected in place above and marked where they occur. **Read section 8 before quoting
anything.**

The five surviving methods, the 29/29 paired `TM_CCOEFF` result, the rotation-augmentation
result and both click controls are unaffected. Sections 1-4 were written before
`tm_variant_sweep.py` produced a single grid row, so section 5 can be read as a test of them
rather than a description of them. **Not committed, though** -- an earlier draft said
"committed", and nothing in this experiment is in git: the script, this log and every CSV are
untracked, and the `midog_utils` changes are modified-uncommitted, so no result here ties to a
code version. That is defect M10 from `2026-09-01-premise-test-audit.md` repeating. File
mtimes do corroborate the mid-run sort fix disclosed in section 7 (`template_match.py` 19:47,
between the script at 19:18 and exp1's CSV at 20:16), and an independent audit found no
prediction that could only have been written after the results -- but mtimes are not
provenance. The only numbers below section 4 come from the gates
and from two pilot probes, both named as such.

---

## 0. Definitions -- every term this log uses

Written because an earlier draft used all of these without defining any of them.

**ROI** -- one downloaded MIDOG++ image, ~2 mm^2, 31-39 megapixels. 13 of them ran here.
**Domain** -- tumour type. 7 exist; each has 2 ROIs downloaded except human breast cancer.
**Seed / click** -- one mitotic annotation chosen to stand in for the pathologist's click. Its
51x51 px neighbourhood becomes the template.
**Cell** -- **one (ROI, seed) pair**, i.e. one run of one pipeline on one image from one click.
*Not* a biological cell. 13 ROIs x up to 5 seeds = **58 cells**; the 6 decision-grade ROIs give
**29**.
**Arm** -- one pipeline variant being scored: a method, a template source, and a ranking key.
`tm_ccoeff` = click-seeded, TM_CCOEFF, ranked by match score. `_od` suffix = same candidates,
ranked by chromatin density instead. `disc_*` = seeded from a synthetic disc (no click).
`lookalike_*` = seeded from a non-mitotic annotation. `blob_native` = no click at all.
**Pool** -- the ranked candidate list an arm produces for one cell, after thresholding, local-
maximum extraction and NMS. Typically 5000-20000 points on a 39 MP ROI.
**Candidate** -- one point in that list: a location the pipeline is asking the pathologist to
look at.
**`read_50` / `read_95` / `read_100`** -- **how many candidates you must read down the ranked
list to have seen 50% / 95% / 100% of that ROI's mitotic figures.** Lower is better. `NaN` when
that fraction is never reached at any depth. This is the headline cost metric: it is literally
"how much work does the pathologist do".
**`full_list_recall` / ceiling** -- fraction of that ROI's mitotic figures that appear anywhere
in the pool, at any depth. Answers "could this pipeline find them at all", separately from
"how deep do you have to read".
**`n_matched`** -- for one cell, the size of the *smallest* pool among the arms being compared.
Every arm is also emitted truncated to that length (`@matched`), so ceilings are compared at
equal list length rather than one arm being credited for simply emitting more candidates.
**`coverage_frac`** -- fraction of the ROI's *area* lying within one match radius (~30 px) of
some candidate, sampled on a stride-16 grid. A pool covering 95% of the image "finds" any
annotation by geometry, so a ceiling is only meaningful read next to this.
**Match radius** -- 7.5 um in pixels (29.6-33.1 px here); MIDOG's own scoring tolerance. A
candidate counts as finding an annotation if it lands within this distance.
**Decision grade** -- ROIs with >= 15 mitotic figures, the only ones any conclusion rests on:
201, 202, 245, 246, 300, 301. All canine, three domains.
**P1-P4** -- the four predictions registered in section 3 before the run.
**`match_pool`** -- the function that runs one match: template -> response map -> threshold ->
local maxima -> NMS -> ranked pool.
**Clamp gate** -- a check that no arm's ranking is mostly one block of tied scores, since tied
candidates are in arbitrary order.

## 1. Why

Every operating-point experiment in this project so far -- score threshold, template size,
search channel, augmentation count, NMS ordering, bbox tightening -- came back neutral. The
2026-08-31 chromatin log gave the reason: each was "a monotone operating-point move on one
fixed ranking function. None could have added information, so none did." And that ranking
function, `TM_CCOEFF_NORMED`, mean-centres and L2-normalises both template and window, so it
is invariant to `I -> aI + b`. A 30%-contrast copy of a patch scores 0.999999 against the
original. Meanwhile `results/morph_diag_bhattacharyya.csv` measures `mean_intensity` as the
strongest mitotic-vs-ordinary-nucleus feature in **all seven** domains (0.925-3.703).

The scoring function is invariant to the strongest measured discriminator in the dataset.

The response so far has been to bolt the missing signal back on afterwards, as a chromatin
re-rank of the candidates the matcher produced (median 19.4x less reading, 29/29 wins). This
experiment asks the question one level down instead: **is there a `cv2.matchTemplate` variant
that sees chromatin density natively?**

It has never been asked here. `grep -roh "TM_[A-Z_]*"` over every `.py`, `.md` and `.ipynb` in
the repo returns 61 hits, all `TM_CCOEFF_NORMED`. The other five appear zero times -- not in
`midog_utils/template_match.py`, not in any log, and not in the ancestor
`bbox tuning code reference/bbox_tuning.py`, which is CCOEFF_NORMED-only across all five of its
call sites.

## 2. The six variants are a 2x3 factorial over the two invariances that matter

Verified numerically rather than quoted (`W = aT + b` on a 51x51 random T; a value that does
not move between rows is an invariance):

| variant | `a=1,b=0` | `a=0.3` | `b=+0.5` | gain-invariant | offset-invariant |
|---|---|---|---|---|---|
| `TM_SQDIFF` | 3.7e-05 | 421.8 | 650.3 | no | no |
| `TM_SQDIFF_NORMED` | 4.3e-08 | 1.000 | 0.4187 | no | no |
| `TM_CCORR` | 860.9 | 258.3 | 1506.4 | no | no |
| `TM_CCORR_NORMED` | 1.000 | 1.000 | 0.9699 | **yes** | no |
| `TM_CCOEFF` | 220.1 | 66.03 | 220.1 | no | **yes** |
| `TM_CCOEFF_NORMED` | 1.000 | 1.000 | 1.000 | **yes** | **yes** |

So the two axes are not arbitrary. **Normed vs unnormed is the "can this method see chromatin
magnitude" axis**, and CCORR/CCOEFF/SQDIFF is how the local mean is treated.
`TM_CCOEFF_NORMED` -- the incumbent -- is the one cell of the table that discards both.

**`TM_SQDIFF_NORMED` is clamped to [0, 1].** OpenCV matches the hand formula wherever it is
below 1 (a=0.5 -> 0.500, a=0.8 -> 0.050, a=1.5 -> 0.1667) and returns exactly 1.0 wherever the
formula exceeds it (a=0.1, 0.3, 3.0). That is a saturating tie block at the *worst* end of the
ranking, which is exactly where a depth-to-100% statistic reads. It is gated, not assumed
benign -- see section 4.

## 3. Predictions, registered before the run

The load-bearing prior is `results/click_darkness.csv` (`click_darkness_probe.py`, un-logged
and un-audited, so treated as indicative). It made darkness **click-conditioned** --
`od51_sim = -|d - d_seed|` -- on a fixed candidate set, and lost to plain absolute darkness by
**30-345x on all 7 ROIs**. The mechanism its own docstring gives: darkness is monotone with
respect to the label, so a two-sided similarity "pushes down every candidate *darker* than
whichever mitosis the pathologist happened to click".

**P1. The two-sided variants lose.** `TM_SQDIFF` and `TM_SQDIFF_NORMED` are that same
two-sided function, `sum((T-I)^2)`, and are predicted worst. A pilot shows the mechanism
directly: over `nucleus_blobs` centroids, *negated* `TM_SQDIFF` correlates **negatively** with
local chromatin density (Spearman -0.54 on 301.tiff, -0.21 on 246.tiff). It actively ranks
dark candidates down, because a darker window has a larger `sum(I^2)` and so a larger penalty.

**P2. `TM_CCORR` is partly degenerate, and ROI-dependently so.** `sum(T*I)` on a non-negative
OD image is a click-shaped matched filter for darkness -- with a constant template it is
exactly a box filter. Pilot Spearman(click template, disc template) over real competing
candidates: **0.58 on 301.tiff, 0.88 on 246.tiff**; against a plain 51 px box filter of OD,
0.55 and 0.93. So on the denser ROI the click still steers the ordering and on 246.tiff it
barely does. Predicted: strong read-depth, but the disc control closes most of the gap.

**P3. `TM_CCOEFF` vs `TM_CCOEFF_NORMED` is the headline, and it is run first and alone
(`--exp1`).** It is the only single-variable contrast in the grid -- identical function, gain
normalisation the only difference -- so the delta isolates *does letting the matcher see
chromatin magnitude help?* with nothing else moving. Both are strongly click-specific in the
pilot (Spearman to the disc template 0.19-0.21 on 301.tiff, -0.05 on 246.tiff), so neither is
a disguised blob detector. Predicted: `TM_CCOEFF` better, by less than the 19.4x that the
explicit chromatin re-rank already delivers.

**P4. TM should beat the blob detector's *ceiling*.** `nucleus_blobs` caps at 0.966 on
246.tiff and 0.889 on 201.tiff, and `results/miss_attribution.csv` shows why: all eleven
problem annotations sit in the **97th percentile** of chromatin density for their ROI, so the
cause is `label(connectivity=2)` fusing touching nuclei, not faintness. TM peak extraction has
no connected-component step. Goal (1) -- find *every* mitotic figure -- is where TM can
genuinely win, if the ceiling is read at matched pool size.

**The bar all four are measured against.** `results/tm_vs_blob_depth.csv` (also un-logged and
un-audited) has `TM_CCOEFF_NORMED` losing to the *click-free* `blob_native` score on **7/7
ROIs at every quantile**: 4.3-216x at read-50, 1.9-12.9x at 95%, 1.0-5.7x at 100%.

**Caveat on the pilot Spearmans, which is why they only order the predictions.** Computed
first over 200k *uniform* ROI locations, `TM_CCORR` looked completely degenerate
(rho(click, disc) = 0.979). Restricted to the ~20k `nucleus_blobs` centroids that actually
compete for rank, it falls to 0.58. Uniform sampling is dominated by near-white background
(median OD 0.0144), where every template agrees, which inflates rho and says nothing about
ordering among candidates.

## 4. Design decisions that could otherwise be made after seeing the numbers

* **Channel: `chromatin.hematoxylin_od`, unclipped.** Not `channels.to_hematoxylin`, which
  min-max rescales between the ROI's 0.5/99.5 percentiles and clips, parking 0.5% of pixels on
  the 255 ceiling -- and those are the dense-chromatin pixels the unnormalised variants exist
  to read. Using it would delete the signal under test.
* **No threshold is ever compared across methods.** Each response map is reduced to its own
  median and MAD and the extraction floor is `med + 0.5*mad`. The robust z is monotone, so it
  reorders nothing; it only makes one floor mean the same search depth for methods whose ranges
  differ by eight orders of magnitude.
* **Template fixed at 51 px, no Otsu tightening**, so the only thing moving on the method axis
  is the method. `--exp1` keeps the shipped configuration (rgb, tightened) as a separate
  reference row.
* **Fusion is over rotations and flips at one scale only.** The unnormalised methods scale with
  template pixel count, so an element-wise max over a bank of *different sizes* is decided by
  size rather than fit. Same-size banks are safe.
* **`ceiling` is only read at matched pool size.** A deep pool tiles the ROI and reaches 1.0 by
  geometry -- the claim this repo already retracted once (coverage_frac 0.91-0.97;
  `random_in_tissue` reaching full-list recall 1.000 on 002 and 350). Every arm is therefore
  also emitted truncated to the smallest pool in its comparison.
* **The clamp threshold is pre-committed at `largest_tie_block / n_pool > 0.05`**, above which
  an arm's `read_100` is reported as arbitrary rather than as a measurement. Measured on
  candidate centroids, `TM_SQDIFF_NORMED` is at **0.199 on 301.tiff** and 0.0026 on 246.tiff,
  so the rule bites on real data and is not a formality.
* **The unreached rule is pre-committed.** A seed that never reaches a quantile is not imputed
  and not silently dropped; spread is computed over reached seeds and `n_unreached` is carried
  beside it. This is the defect the premise-test audit found, where a reported "28" was the
  median of {28, 28, 200}.
* **5 seeds, drawn without replacement, worst-of-5 as the headline** -- matching
  `tail_sensitivity.csv`, `click_darkness.csv` and `tm_vs_blob_depth.csv` so the numbers are
  comparable to the bar. ROIs delivering fewer than 5 distinct seeds are excluded from the
  headline table and kept in the CSV.

### Scope, and what it can actually decide

13 ROIs, all 7 domains, 2 per domain except human breast cancer (001.tiff has zero mitotic
figures and cannot be seeded). Only 301/300/246/245/201/202 carry n_mitotic >= 15, and 301/300
share tumour type *and* scanner -- so roughly **three effective domains arbitrate**, and the
rest are labelled `decision_grade=False` and drive no conclusion.

**And those six are narrower than "7 domains" suggests: all six are canine.** Three tumour
types over two scanners and two labs -- 3D Histech / VMU Vienna for 201/202/245/246, Aperio
CS2 / FU Berlin for 300/301. **No human ROI is decision grade**, and every human ROI in this
dataset is Hamamatsu XR / UMC Utrecht, so that scanner contributes *nothing* to any conclusion
here. This matters for a product aimed at human histopathology, and "blobs own 300 and 301"
rests on one tumour type on one scanner. Every conclusion in section 5d should be read with
that in front of it.

**Annotation crowding turns out not to constrain anything** -- though the argument originally
given here for that was wrong, see the retraction in section 5. Only 245.tiff has two mitoses
closer together than the match radius (26.6 px against 30.2 px), and both are in fact matched
in 22 of its cells. No arm is limited by annotation geometry on any ROI, so a variant that
stalls short of 1.0 has failed on its own merits.

### Gates, all passing before the grid ran

25 gates, 0 failures (`results/tm_variant_gates.csv`):

* `plant_and_recover` per method x 2 template banks (12): the fused map's peak lands within
  0 px of where the template was planted. Catches the per-map half-size offset bug, and any
  sign error -- un-negated `TM_SQDIFF` puts the planted patch at the map's *minimum*.
* `matchtemplate_precision` (6): OpenCV against float64 brute force at 200 random windows.
  Worst relative error 7.3e-05, and the **unnormalised methods are the more precise ones**
  (~6e-08 median) -- the DFT path is not a problem here.
* `z_is_monotone` (6): Spearman(sign x raw, robust-z) = 1.000000 exactly, so no read-depth is
  measuring a different ranking than the one named in the CSV.
* `regression_tie_ccoeff_normed` (1): the refactored `fused_response` reproduces the
  pre-refactor implementation's scores and its 2406-detection list **identically** on
  002.tiff. Without this, any exp1 delta could be an artefact of the refactor.

Per-ROI gates run inside the grid: `verify_shortcut` per method (coordinate-set equality for
the one-match-many-z economy), tissue-mask coverage, NMS radius, distinct seed streams, and
no-cap.

## 5. Results -- Experiment 1: `TM_CCOEFF` vs `TM_CCOEFF_NORMED`

13 ROIs x 5 seeds, `results/tm_ccoeff_headtohead.csv`, 8048 rows, **1058 invariant checks, 0
failed**. Decision-grade means n_mitotic >= 15, which is 201/202/245/246/300/301 -- 29 paired
seeds (202.tiff delivered 4 of 5, its seed pool being smaller than its mitotic count).

### P3 confirmed, and not marginally

Both arms saw the identical click on the identical ROI, so the paired delta carries no seed
variance at all:

| metric | pairs | `TM_CCOEFF` better | median delta | Wilcoxon p |
|---|---:|---:|---:|---:|
| `read_50` | 29 | **29** | -1338 candidates | < 1e-5 |
| `read_95` | 25 (4 unreached) | **25** | -3780 | < 1e-5 |
| `read_100` | 8 (21 unreached) | 8 | -1204 | 0.0078 -- **see the warning below** |
| `full_list_recall` | 29 | 10 better, 16 tied, 3 worse | +0.0042 mean | 0.033 |

> **The `read_100` row is conditioned on success and must not be quoted as a tail result.**
> Those 8 pairs are the cells where *both* arms reached full recall; `TM_CCOEFF` reached it in
> 11 of 29 cells on its own and `TM_CCOEFF_NORMED` in 8. Selecting the cells where the harder
> arm succeeded is exactly the defect the premise-test audit caught, where a reported "28" was
> the median of {28, 28, 200}. Goal (1) is settled by `full_list_recall` below, which is
> defined on all 29 pairs and conditions on nothing.

**Every single one of 29 paired seeds favours removing the gain normalisation.**

Those 29 cells are 6 ROIs x 5 seeds, so they are clustered, not independent, and the Wilcoxon
p-values above overstate their evidence. Collapsing to one value per ROI and running a sign
test on the 6 ROIs -- the conservative version -- gives: **`read_50` 6 ROIs better, 0 worse,
p = 0.031, which survives**; `full_list_recall` 3 better, 0 worse, 3 tied, **p = 0.25, which
does not**. So the depth result stands at ROI level and the ceiling result does not have an
independent statistical basis. Six ROIs is the floor of what a sign test can resolve; 0.031 is
the smallest p obtainable.

Worst-of-5
`read_50`, the number a pathologist actually pays:

| ROI | `TM_CCOEFF` | `TM_CCOEFF_NORMED` | **method gain** | `blob_native` | `pipeline_as_shipped` |
|---|---:|---:|---:|---:|---:|
| 201.tiff | **194** | 2256 | 11.6x | 81 | 4767 |
| 202.tiff | -- | -- | -- | -- | -- |
| 245.tiff | **2943** | 5997 | 2.0x | 945 | 6681 |
| 246.tiff | **156** | 534 | 3.4x | 194 | 1876 |
| 300.tiff | **371** | 3496 | 9.4x | 143 | 5538 |
| 301.tiff | **526** | 4110 | 7.8x | 168 | 6216 |

So the answer to the question this experiment was built to ask is **yes**: letting the matcher
see chromatin magnitude is worth **2.0-11.6x** the reading burden at 50% sensitivity. That is
the method effect, and it is clean -- same channel, same 51 px template, same floor, same NMS,
same seeds, only the `cv2.TM_*` constant differs.

**The `pipeline_as_shipped` column is not a method comparison and must not be quoted as one.**
It differs from `TM_CCOEFF` in three ways at once -- method, channel (`rgb` ->
`hematoxylin_od`) and template size (Otsu-tightened -> fixed 51 px) -- so its 2.3-24.6x spread
is a *pipeline* delta, not a method delta. How much of it is not the method is directly
visible: the shipped arm is worse than `TM_CCOEFF_NORMED` on **every** ROI (4767 vs 2256,
6681 vs 5997, 1876 vs 534, 5538 vs 3496, 6216 vs 4110), and those two share a method. The
channel and tightening changes therefore account for a large share of the shipped gap on their
own.

### The configuration to actually ship -- with the caveat exp2 later forced onto it

*(Convention for this table: `read_50` is the median of per-ROI medians over the five
five-seed decision-grade ROIs; the win counts are worst-of-5 per ROI. Other sections pool the
29 decision-grade cells instead, which is why numbers move slightly between tables -- e.g. the
median of the worst-of-5 values for `tm_ccoeff_od` is 139, not 135.)*

Not `TM_CCOEFF` ranked by its own score. **`TM_CCOEFF` as the candidate generator, ranked by
chromatin density** (`tm_ccoeff_od`) is better than both the incumbent and the click-free bar:

| arm | median `read_50` | wins vs `blob_native` | cv over seeds |
|---|---:|---:|---:|
| `tm_ccoeff_od` | **135** | **4 of 5 ROIs** | 0.098 |
| `blob_native` (no click) | 168 | -- | 0.008 |
| `blob_od` | 191 | 2 of 5 | 0.006 |
| `tm_ccoeff` | 202 | 1 of 5 | 0.348 |
| `tm_ccoeff_normed` (incumbent) | 1672 | 0 of 5 | 0.484 |

It beats `blob_native` at `read_50` on 201 (37 vs 81), 245 (595 vs 945), 246 (130 vs 194) and
300 (139 vs 143), losing only 301 (272 vs 168) -- and it is **3.5x more stable across seeds**
than ranking by the correlation score. This is the first configuration in this project where a
click-seeded pipeline beats the seedless blob detector on a matched candidate set;
`results/tm_vs_blob_depth.csv` had `TM_CCOEFF_NORMED` losing 4.3-216x.

It is also, notably, the same conclusion `7c3af93` reached from the other direction: the search
earns its place as the *candidate generator*, not as the ranker.

**But exp2 shows the method swap is not what buys this.** Ranked by chromatin,
`tm_ccoeff_od` and `tm_ccoeff_normed_od` are indistinguishable at `read_50` -- 8 of 29 paired
cells better, 3 tied, median delta **+1 candidate, p = 0.20**. The 135-vs-168 win over
`blob_native` belongs to "TM peaks as a generator, chromatin as a ranker", and would have been
available with the incumbent method. What the method swap does buy under chromatin ranking is
smaller and further out: `read_95` improves on 16 of 25 cells (median -91, p = 0.0125) and the
ceiling on 10 of 29 (p = 0.033, 11 full-recall cells against 8). See section 5b.

### Three qualifications, all of which matter

**1. The advantage is a head advantage, not a tail advantage.** At `read_95` the same gain
falls to **1.1-1.9x**, and at `read_100` only 8 of 29 cells reach full recall at all. The
click helps you clear the easy half of the mitoses quickly; it does very little for the last
few, which is where the reading cost actually lives (`read_95` medians are in the thousands).

**2. It still does not clear the click-free bar, though it nearly does** -- *true only of the
single-template configuration; see section 5e, where the augmented one clears it on 5 of 6
ROIs.*
`TM_CCOEFF` beats `blob_native` on only 1 of 5 ROIs at `read_50` (246.tiff, 156 vs 194) and is
2.4-3.1x behind on the other four. But the median-of-ROI-medians is **202 for `TM_CCOEFF`
against 168 for `blob_native`** -- against 4.3-216x behind for `TM_CCOEFF_NORMED` in
`results/tm_vs_blob_depth.csv`. The gap that motivated abandoning template matching as a
candidate generator has closed from two orders of magnitude to 1.2x.

And the combination wins outright: **TM candidates ranked by chromatin density
(`tm_ccoeff_od`) beat `blob_native` on 4 of 5 ROIs** -- 37 vs 81 (201), 595 vs 945 (245),
130 vs 194 (246), 139 vs 143 (300) -- losing only on 301 (272 vs 168). Median 135 vs 168.

**3. `TM_CCOEFF` is far more seed-sensitive than the thing it beats.** Median within-ROI
coefficient of variation of `read_50` across seeds:

| arm | cv (median) | cv (max) |
|---|---:|---:|
| `blob_native` | 0.008 | 0.079 |
| `tm_ccoeff_od` | 0.098 | 0.231 |
| `tm_ccoeff` | **0.348** | **0.746** |
| `tm_ccoeff_normed` | 0.484 | 0.762 |

The seedless arms are *nearly* deterministic -- their candidate list is byte-identical across
seeds, and the residual cv comes from the evaluation ground truth, which drops a different
annotation for each seed. (An earlier draft called them "deterministic by construction", which
this table's own non-zero cv contradicts.) So this is not a fair contest -- but for a
tool where the pathologist gets *one* click, a 35% coefficient of variation is the number that
decides whether the median is worth anything. Ranking the TM candidates by chromatin instead of
by correlation cuts it to 0.098 while also being faster to read, which is the same conclusion
`7c3af93` reached by a different route.

### P4 is confirmed on one ROI, and refuted as a general claim

**Goal (1) first, on the metric that conditions on nothing.** Cells (of 29 decision-grade
ROI x seed pairs) where the arm's pool contained *every* mitotic figure. Section 4 of this log
pre-committed to reading ceiling **at matched pool size**, because an unmatched pool reaches
1.0 partly by tiling the ROI -- so the matched column is the one that counts:

| arm | full pool | **matched pool (the pre-committed metric)** |
|---|---:|---:|
| `blob_native` | 19 / 29 | **15 / 29** |
| `tm_ccoeff` | 11 / 29 | **11 / 29** |
| `tm_ccoeff_normed` | 8 / 29 | **7 / 29** |

`TM_CCOEFF` improves the ceiling over `TM_CCOEFF_NORMED` (11 vs 7 matched cells, 10 better /
16 tied / 3 worse over 29 cells) -- but that survives only the clustered test, not the ROI-level
sign test (p = 0.25), and **the click-free blob detector still achieves goal (1) more often than
either**. P4 predicted TM would beat the blob ceiling generally; over the six decision-grade
ROIs it does not. What follows is the one ROI where it does.

> **Correction.** An earlier version of this log reported only the full-pool column (19 / 11 /
> 8) as the headline, which is the aggregation section 4 explicitly pre-committed *against*.
> The direction of the finding is unchanged, and the gap in blob's favour is narrower at
> matched pool, not wider. Found by independent audit. Blob's full pool is up to 2x TM's
> (34910 against ~17000 on 245.tiff), which is exactly why the rule existed.

Ceiling at **matched pool size** (every arm truncated to the same list length on that ROI):

| ROI | `tm_ccoeff` | `tm_ccoeff_normed` | `blob_native` | geometric cap | coverage (TM / blob) |
|---|---:|---:|---:|---:|---|
| 201.tiff | **1.0000** | 0.9412 | 0.8824 | 1.0000 | 0.83 / 0.76 |
| 202.tiff | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.84 / 0.75 |
| 245.tiff | 0.9888 | 0.9888 | 0.9888 | **0.9889** | 0.95 / 0.77 |
| 246.tiff | **0.9913** | 0.9826 | 0.9478 | 1.0000 | 0.95 / 0.85 |
| 300.tiff | 0.9889 | 0.9889 | **1.0000** | 1.0000 | 0.93 / 0.85 |
| 301.tiff | 0.9954 | 0.9862 | **1.0000** | 1.0000 | 0.95 / 0.88 |

TM beats the blob ceiling on the two ROIs where under-segmentation was documented to break it
-- 201.tiff (0.882) and 246.tiff (0.948 at matched pool) in `results/miss_attribution.csv`,
where `label(connectivity=2)` fuses touching nuclei and the merged centroid lands between them.
Template matching has no connected-component step and recovers those. But only **201.tiff
survives the coverage objection**: +0.12 ceiling for +0.07 coverage. On 246.tiff the gain is
+0.043 ceiling for +0.10 coverage, i.e. roughly proportional to the extra ROI the pool tiles,
so it is not independent evidence. **P4 is confirmed on 201.tiff and unresolved on 246.tiff**,
and on the two mast-cell ROIs, where blobs already reach 1.0, TM gives up 0.005-0.011.

> **RETRACTED -- the `geometric_ceiling` column is not a valid upper bound, and the 245.tiff
> "honesty check" built on it was wrong.** An earlier version of this section argued that all
> three arms landing on 0.9888 against a computed cap of 0.9889 showed nothing was being
> credited that geometry forbids. Two things are wrong with it. (i) The two numbers are
> different fractions on different denominators -- 88/89 = 0.988764 measured against `gt_eval`,
> and 89/90 = 0.988889 computed over all mitotic annotations including the seed -- so their
> agreement was a coincidence, not a match. (ii) More seriously, the bound itself is false.
> `geometric_ceiling` assumed that two mitoses closer together than the NMS radius cannot both
> be matched. They can: the detections need not sit on the ground-truth centres, so two peaks
> can be more than a radius apart from each other while each lies within a match radius of a
> different annotation. This is not hypothetical -- on 245.tiff **22 cells reach
> `full_list_recall` = 1.0000**, i.e. both members of the crowded pair (ann 6253/6254, 26.6 px
> apart, radius 30.2 px) were matched, and neither is ever the seed.
>
> `geometric_ceiling` and `n_mitotic_crowded` are **reported columns only** -- they are never
> used as a filter, a cap, or an input to any metric (`grep` confirms the only consumers are a
> log line and a display column), so no other number in this log is affected. But the claim
> "goal (1) is geometrically reachable on 12 of 13 ROIs and capped at 0.989 on 245.tiff",
> made in the Scope section, is unsupported and should be read as "no geometric obstruction was
> established either way". Found by independent audit, not by me.

**The caveat that has to travel with this table**: at equal list length TM still covers more of
the ROI than the blob detector does (0.83-0.95 against 0.75-0.88), because blob centroids
cluster. Equal `n` is not equal coverage, so part of TM's ceiling advantage is geometry rather
than detection -- which is why only the 201.tiff margin is claimed above.

> **Bug in the blob coverage column, found by independent audit.**
> `_method_arms` cached `blob_native@matched`'s coverage under a key that omitted the seed
> index, while the truncation depth `n_matched` varies between seeds -- so one seed's coverage
> was reported for all five. This is exactly the mis-keying `compare.evaluate_arms`' own
> docstring warns about, and it landed on the column used to adjudicate P4. Fixed in the code;
> the correct per-cell values are recomputed in
> `results/tm_variant_blob_matched_coverage_fix.csv` without re-running the grid, since the
> blob list is deterministic per ROI. **The error is at most 0.024 in this table's CSV and
> exactly 0.0000 on 201.tiff**, the only ROI carrying a P4 claim, so the conclusion is
> unaffected. In exp2's CSV, where the matched pool is shallower and varies more, it reaches
> 0.18 -- so do not quote blob coverage from `tm_variant_sweep.csv` without the fix file.

### The click is doing real work -- at the head

Ratio of the click-free disc template's depth to the real click's, same ROI, same method. 1.0
would mean the click is worthless:

| ROI | `ccoeff` @ read_50 | `ccoeff` @ read_95 | `ccoeff_normed` @ read_50 |
|---|---:|---:|---:|
| 201.tiff | 7.5x | 0.95x | 6.2x |
| 202.tiff | **46.6x** | 1.5x | 27.8x |
| 245.tiff | 1.6x | 1.1x | 1.5x |
| 246.tiff | 5.5x | 2.0x | 13.4x |
| 300.tiff | 7.6x | 1.3x | 4.6x |
| 301.tiff | 9.1x | 1.9x | 5.9x |

So P2's degeneracy worry does **not** apply to `TM_CCOEFF`: it is not a dark-blob detector
wearing a template. Properly **paired per (ROI, seed)** rather than taken as a ratio of two
medians, and counting only cells where both arms reached the quantile:

| method | paired median disc/click | cells | verdict |
|---|---:|---:|---|
| `TM_CCOEFF` | **6.89x** | 29/29 | click clearly helps |
| `TM_CCOEFF_NORMED` | 5.91x | 29/29 | click clearly helps |
| `TM_CCORR` | 3.12x | 29/29 | click helps |
| `TM_CCORR_NORMED` | 1.47x | **5/29** | not a measurement |
| ~~`TM_SQDIFF_NORMED`~~ | ~~0.71x~~ | ~~15/29~~ | withdrawn -- degenerate pool |
| `TM_SQDIFF` | **0.67x** | 22/29 | **click is worse than a disc**, p = 0.030 |

Three corrections to an earlier draft here. The headline ratio for `TM_CCOEFF` is **6.9x, not
8.1x** -- the 8.1x was a ratio of unpaired medians. The range "1.6-46.6x" was wrong twice
over: its maximum came from 202.tiff, which this log's own four-seed rule excludes from the
adjacent worst-of-5 table, so it was included in one table and excluded in the next; the
per-ROI paired range excluding it is **1.6-9.1x**. (`TM_SQDIFF_NORMED`'s ratio also flipped
sign when paired, but both its arms are built on the degenerate pool, so it is withdrawn rather
than reported.)

The `read_95` collapse still holds directionally, but note it is measured over unequal reached
subsets (201.tiff's click median is over 3 seeds against the disc's 5), which biases it in the
click's favour -- so "the click is worth nothing at the tail" is if anything understated.

Worth noting separately: the click-free `disc_ccoeff` (median `read_50` 1529) still beats the
click-seeded `tm_ccoeff_normed` (1672). Part of the win is the click and part is simply seeing
magnitude at all; those are separable and both real.

## 5b. Results -- Experiment 2: all six variants

13 ROIs x 5 seeds x 6 methods plus both controls. `results/tm_variant_sweep.csv`, 19968 rows,
**2600 invariant checks, 0 failed**, 56 min. Medians over the 29 decision-grade ROI x seed
cells, click-seeded, ranked by the method's own score:

| rank | method | `read_50` | `read_95` | ceiling | full-recall cells | click worth (disc/click) |
|---|---|---:|---:|---:|---:|---:|
| 1 | **`TM_CCOEFF`** | **189** | 6544 | **0.9913** | **11 / 29** | **6.9x** |
| 2 | `TM_CCORR` | 501 | 8788 | 0.9739 | 8 / 29 | 3.1x |
| 3 | `TM_CCOEFF_NORMED` *(incumbent)* | 1623 | 10462 | 0.9888 | 8 / 29 | 5.9x |
| 4 | `TM_CCORR_NORMED` | 3713 | 13870 | 0.9652 | 5 / 29 | 1.5x ‡ |
| 5 | `TM_SQDIFF` | 8149 † | 15493 †† | 0.8667 | 2 / 29 | **0.67x** |
| -- | ~~`TM_SQDIFF_NORMED`~~ | ~~2364~~ | ~~11509~~ | ~~0.9565~~ | ~~5 / 29~~ | ~~1.2x~~ |

Click worth is the **paired** median disc/click ratio (see the table further down), not a ratio
of medians. ‡ `TM_CCORR_NORMED`'s is over 5 of 29 cells and is not a measurement.

> **`TM_SQDIFF_NORMED` is WITHDRAWN from this table and from every other result in this log.**
> Its candidate pool is a tie-break artefact, not a measurement. OpenCV clamps
> `TM_SQDIFF_NORMED` to [0, 1], and on hematoxylin OD **83% of the response map saturates
> there** -- so after the sign flip the map's median is exactly -1.0 and its **MAD is exactly
> 0.0**. The extraction floor `med + 0.5*mad` therefore lands *on* the clamp for any z,
> admitting ~28 million tied local maxima; the `max_peaks = 250000` cap keeps the first 250k,
> and because `extract_peaks` breaks ties by x-coordinate the survivors are a strip down the
> left edge of the ROI. **97.4% of the retained peaks carry the identical score.** Measured on
> 002.tiff seed 0; the cap bound on **228 rows across 10 of 13 ROIs**, every one of them a
> `sqdiff_normed` variant.
>
> This also propagates: `n_matched` is a minimum over arms, and `sqdiff_normed` set that
> minimum in **20 of 58 cells**, so **every `@matched` number in exp2 is truncated to a length
> partly fixed by the artefact**. No `@matched` figure from `tm_variant_sweep.csv` is quoted in
> this log, and **exp1 and exp3 are entirely clean** -- zero capped rows, zero collapsed floors,
> and neither run includes the method. The headline ceiling comparisons in sections 5 and 5c
> come from those two runs.
>
> `match_pool` now raises `DegenerateMapError` on `mad <= 0` and on `n_peaks >= max_peaks`
> rather than returning a pool. Found by independent audit; see section 8, finding A.

† 3 of 29 cells never reached 50%; †† **23 of 29** never reached 95%, so that cell is not a
usable summary. Other medians in this table with material unreached counts, which section 4's
pre-committed rule says must travel with them: `tm_ccoeff` `read_95` 3/29 unreached,
`tm_ccorr_normed` `read_95` 13/29, `tm_sqdiff_normed` `read_95` 13/29, `disc_ccorr_normed`
`read_50` 24/29. The `read_50` column is fully reached for every arm except `TM_SQDIFF`
(3/29) and `TM_SQDIFF_NORMED` (1/29), so the ranking itself is sound.

### P1 confirmed, including the sign

`TM_SQDIFF` is last on every column of the five measurable methods. The registered prediction
named *both* two-sided variants, and only one of them can be evaluated at all: the other,
`TM_SQDIFF_NORMED`, is withdrawn above because its pool is degenerate. So **P1 is confirmed for
the two-sided method that produced a valid pool, and untested for the other** -- and the reason
the other is untestable is itself a property of the clamp, not an accident of this run.
The mechanism registered in section 3 shows up exactly as predicted: **its click-free disc control
reads 4750 candidates against the real click's 8149, so under `TM_SQDIFF` the pathologist's
click is actively worse than a blank disc** (0.58x). This is `od51_sim` again -- a two-sided
score penalising every candidate darker than whichever mitosis happened to be clicked -- now
reproduced through a completely different mechanism, `sum((T-I)^2)` rather than
`-|d - d_seed|`. `TM_SQDIFF` also has the worst ceiling by a wide margin (0.8667, full recall
in 2 of 29 cells), so it fails goal (1) as well as goal (2).

### P2 refined: `TM_CCORR` is not degenerate, but it is second best

The uniform-sampling pilot suggested `TM_CCORR` was a box filter in disguise (rho 0.979). It
is not: the click is worth **3.7x** against a shape- and magnitude-matched disc, and the disc
arm's ceiling collapses to 0.90 against the click's 0.974. But it does not beat `TM_CCOEFF`
anywhere, so "keep the local mean, drop the normalisation" loses to "drop the local mean, keep
the magnitude".

### The clamp gate fired -- and was still far too weak

> Superseded by the withdrawal above. The tie-fraction rule below did flag the right arm, but
> it flagged a *symptom* and let the disease through: a pool that is 97% one tie block is not
> "a tail reported as arbitrary", it is not a pool at all. The rule caught `disc_sqdiff_normed`
> and missed `tm_sqdiff_normed_od` entirely (`tie_frac` is computed on the *ranking* key, and
> the chromatin key has no ties), so 23 of the 38 capped rows passed every gate. The real check
> was the one that did not exist: assert the map has non-zero MAD and that `max_peaks` never
> binds. Both now exist.

`largest_tie_block / n_pool > 0.05` was fixed in section 4 before any grid row existed, and it
fires. **Across all 13 ROIs it flags 129 cells over six arms**, every one of them a
`sqdiff_normed` variant and no other method anywhere: `disc_sqdiff_normed` reaches **0.4954**
(351.tiff) and `tm_sqdiff_normed` **0.2291** (405.tiff). Those cells' `read_100` is not
reported as a measurement.

Restricted to the six decision-grade ROIs the picture is milder -- `disc_sqdiff_normed` peaks
at 0.4105 with 19 cells flagged, and the click-seeded `tm_sqdiff_normed` peaks at 0.0414 and is
flagged **zero** times, because there the extraction floor cuts most of the clamped block away
before it reaches the pool. **That is a property of the decision-grade ROIs, not of the arm**:
on the sparser ROIs the floor does not cut it away and the arm under test is flagged in 15
cells. An earlier draft of this section stated the decision-grade figures without the
restriction, which made `TM_SQDIFF_NORMED` look cleaner than it is.

### The finding that most changes what to do next

**Rank by chromatin density and the method choice almost vanishes.** Median `read_50` over the
same 29 cells, same pools, only the sort key changed:

| arm | `read_50` | `read_95` | ceiling | full-recall cells |
|---|---:|---:|---:|---:|
| ~~`tm_sqdiff_normed_od`~~ | ~~128~~ | ~~4654~~ | ~~0.9565~~ | ~~5 / 29~~ |
| `tm_ccoeff_od` | **130** | **4556** | **0.9913** | **11 / 29** |
| `tm_ccoeff_normed_od` | 130 | 5083 | 0.9888 | 8 / 29 |
| `tm_ccorr_normed_od` | 130 | 5394 | 0.9652 | 5 / 29 |
| `tm_ccorr_od` | 140 | 4418 | 0.9739 | 8 / 29 |
| `tm_sqdiff_od` | 142 | 7294 | 0.8667 | 2 / 29 |

(`tm_sqdiff_normed_od` is withdrawn with the rest of that method -- it appeared *best* in an
earlier draft of this table, on the artefact pool.) A 43x spread at `read_50` between the best
and worst measurable method (189 to 8149) collapses to **130 to 142** once chromatin does the
ranking. The variants differ enormously as *rankers* and barely
at all as *generators of the top of the list*. What survives is further out: the **ceiling**
still separates them (0.8667 to 0.9913, and 2 to 11 full-recall cells), and `read_95` mildly.

So the method choice is worth 8.6x if the TM score is the ranking key, and worth almost nothing
at `read_50` if it is not -- which is the same shape of result as `7c3af93`, arrived at from a
new direction.

### And the click's specificity does not survive the controls

Two controls, two different answers, and the pair is the interesting part:

*(medians; mitotic n=29 cells, look-alike n=6, disc n=29 -- see the warning below)*

| method | mitotic click | look-alike click | disc (no click) |
|---|---:|---:|---:|
| `TM_CCOEFF` | 189 | **180** | 1529 |
| `TM_CCOEFF_NORMED` | 1623 | **1068** | 7257 |
| `TM_CCORR` | **501** | 578 | 1868 |
| `TM_CCORR_NORMED` | **3713** | 5952 | 10832 |
| `TM_SQDIFF_NORMED` | **2364** | 5669 | 2724 |
| `TM_SQDIFF` | **8149** | 12436 | 4750 |

> **The table above is not a like-for-like comparison and should not be read as one.** The
> mitotic column pools 29 cells (5 seeds x 6 ROIs); the look-alike column is **6 cells** (the
> control was budgeted at one click per ROI). Comparing a 29-cell median against a 6-cell
> median is not a controlled contrast, and the apparent mitotic advantage on the other four
> methods is largely that artefact -- like-for-like, the mitotic click wins on only 2/6, 2/6,
> 2/6 and 3/6 ROIs.

The statement that survives the aggregation problem, checked per ROI:

> **On all 6 decision-grade ROIs, and for both `TM_CCOEFF` and `TM_CCOEFF_NORMED`, the
> look-alike click's `read_50` falls inside the range spanned by that ROI's five mitotic
> clicks** -- never below the best, never above the worst.

Concretely for `TM_CCOEFF`: 46 against a mitotic range of [40, 194] on 201.tiff; 19 against
[15, 130] on 202; 2457 against [1687, 2943] on 245; 145 against [92, 156] on 246; 216 against
[144, 371] on 300; 248 against [203, 526] on 301. A look-alike click is not distinguishable
from a mitotic click by this measure at all.

Against the disc the click is worth 6.9x, so the template *is* carrying real information -- but
that information is "a dark, nucleus-sized object of roughly this shape", not "a mitotic
figure". The mitotic *identity* contributes nothing measurable under `TM_CCOEFF`. (A uniform
disc has no texture, so this null cannot separate shape from texture; an earlier draft claiming
"the 8.1x is shape and texture" overstated what it isolates.)

That is consistent with, and sharper than, the premise test's +0.052 median recall@1000 for a
mitotic seed over a look-alike seed -- "consistent in sign everywhere, decisive only on
246.tiff". It also matches `2026-09-01-click-ranking-experiment.md`'s diagnosis from the other
end: only a task-trained (and leaky) backbone separated the two clicks; every frozen
representation, raw pixels included, did not.

**Caveat on the look-alike column**: it is one seed per ROI (6 cells) against the mitotic
column's 29, because the control was budgeted at one click per ROI. It is enough to refute a
large specificity effect and not enough to measure a small one.

## 5c. Results -- Experiment 3: rotation augmentation on the top two

`results/tm_variant_stage_b.csv`, 4176 rows, **574 checks, 0 failed**, 20 min. Same 29 cells,
`n_angles=4, flips=(False, True), scales=(1.0,)` -- 8 templates, all the same size, so
element-wise max fusion stays unbiased for the unnormalised methods (section 4).

| arm | `read_50` 1 tpl -> 8 tpl | cells better | p | ceiling 1 tpl -> 8 tpl | full-recall cells |
|---|---|---:|---:|---|---:|
| `tm_ccoeff` | 189 -> **136** | 24 / 29 | 0.0002 | 0.9913 -> **1.0000** | 11 -> **18** |
| `tm_ccorr` | 501 -> 346 | 25 / 29 | 0.0001 | 0.9739 -> 0.9778 | 8 -> 8 |
| `tm_ccoeff_od` | 130 -> 135 | 14 / 29 | 0.33 | 0.9913 -> **1.0000** | 11 -> **18** |

**This is the largest single move on goal (1) in the whole experiment.** Augmented
`TM_CCOEFF` reaches full recall in **18 of 29 cells at full pool against the blob detector's
19**, up from 11 -- and at the pre-committed **matched pool, 12 of 29 against 14**, up from 11.
The move is real on both aggregations; the gap to the blob detector is wider on the one that
counts. Mitotic figures have no canonical orientation, and it turns out that matters
far more for the *ceiling* than for the head of the list -- augmentation does nothing for
`read_50` once chromatin ranks (p = 0.33) but moves the ceiling on the same cells (p = 0.023).

It also mildly corrects an implication of the un-audited `results/tm_vs_blob_depth.csv`, where
rotation "helped on 24 of 42 cells, hurt on 18, never enough to change an outcome". Under
`TM_CCOEFF` it changes an outcome.

### TM and the blob detector fail on *disjoint* ROIs

Cells reaching full recall, per ROI, out of 5 seeds:

| ROI | aug `tm_ccoeff` full | `blob_native` full | **aug TM matched** | **blob matched** |
|---|---:|---:|---:|---:|
| 201.tiff | **3 / 5** | **0 / 5** | **3 / 5** | **0 / 5** |
| 246.tiff | **4 / 5** | **0 / 5** | **3 / 5** | **0 / 5** |
| 202.tiff | 4 / 4 | 4 / 4 | 3 / 4 | 4 / 4 |
| 245.tiff | 4 / 5 | 5 / 5 | 1 / 5 | 0 / 5 |
| 300.tiff | 1 / 5 | **5 / 5** | 1 / 5 | **5 / 5** |
| 301.tiff | 2 / 5 | **5 / 5** | 1 / 5 | **5 / 5** |

(202.tiff delivered 4 seeds, not 5 -- its column is out of 4.)

The two methods are complementary, and the split is mechanistic rather than incidental: TM
solves the two ROIs where `label(connectivity=2)` fuses touching nuclei and the blob detector
*never* reaches full recall at any depth, while the blob detector solves the two dense
mast-cell ROIs where TM's peak extraction misses one or two figures. **201.tiff, 246.tiff,
300.tiff and 301.tiff hold that pattern under both aggregations. 245.tiff does not** -- at
matched pool the blob detector drops from 5/5 to 0/5 there, so it is TM that is (barely) ahead,
and an earlier draft listing 245 as a blob win was reading the full-pool column only.

## 5d. What this changes

0. **Read section 5e first.** The best configuration measured here -- `TM_CCOEFF` with
   rot4 x flip2 augmentation -- beats the click-free blob detector on 5 of 6 decision-grade
   ROIs at `read_50` and ties it on ceiling. Several statements below were written before that
   comparison was run and are scoped to the single-template arms.
1. **If the TM score is the ranking key, switch `TM_CCOEFF_NORMED` to `TM_CCOEFF`.** 8.6x less
   reading at `read_50`, 29/29 paired seeds, a better ceiling, and a one-constant change to
   `FSConfig.tm_method`. There is no argument for keeping the normalisation.
2. **The better pipeline is still TM-as-generator plus chromatin-as-ranker**, but be careful
   what the method swap buys there: at `read_50` nothing (median +1 candidate, p = 0.20), and
   the ceiling and `read_95` advantages **do not survive a ROI-level sign test** (p = 0.25 and
   p = 0.22 over 6 ROIs) even though they clear a per-cell Wilcoxon. Treat them as suggestive.
3. **Turn rotation augmentation back on** for the unnormalised method: `n_angles=4,
   flips=(False, True), scales=(1.0,)`. It costs 8x the match time and takes full-recall cells
   from 11 to 18 of 29, which is the largest single move on goal (1) here. Keep it single-scale
   -- a multi-scale bank breaks unnormalised max fusion outright.
4. **Goal (1) is still not solved, but the route is now visible.** Augmented `TM_CCOEFF`
   reaches full recall in 12 of 29 matched-pool cells and `blob_native` in 14 -- **and they
   fail on disjoint ROIs**. TM owns 201 and 246, where the blob detector never gets there at
   any depth under either aggregation; blobs own 300 and 301. A union generator is the obvious
   next experiment, and `watershed_split.py` already exists for the blob half. Nothing here
   suggests either alone will close it.
5. **Stop budgeting on the assumption that the mitotic click is special.** Two independent
   controls now say the click buys nucleus-shaped priors, not mitotic identity.

## 5e. The comparison this log failed to run

Every "TM does not clear the click-free bar" statement above compares the **single-template**
TM arms against `blob_native`. Experiment 3 produced augmented TM numbers and they were only
ever compared *against single-template TM* -- a within-TM question -- because by that point
this log had already framed the blob detector as the incumbent. Running the comparison that
was actually available the whole time:

**Augmented `TM_CCOEFF` vs the click-free blob detector, `read_50`, median over 5 seeds:**

| ROI | aug TM (chromatin-ranked) | `blob_native` | ratio | aug TM ceiling | blob ceiling |
|---|---:|---:|---:|---:|---:|
| 201.tiff | **25** | 81 | **3.2x** | **1.0000** | 0.8824 |
| 202.tiff | **18.5** | 27 | 1.5x | 1.0000 | 1.0000 |
| 245.tiff | **494** | 945 | 1.9x | 1.0000 | 1.0000 |
| 246.tiff | **117** | 194 | 1.7x | **1.0000** | 0.9652 |
| 300.tiff | **135** | 143 | 1.1x | 0.9944 | 1.0000 |
| 301.tiff | 198 | **168** | 0.85x | 0.9954 | 1.0000 |

Median of ROI medians: **augmented `TM_CCOEFF` 125 by its own score, 126 chromatin-ranked,
against `blob_native`'s 155.5 and `blob_od`'s 170.** Full-recall cells 18/29 against 19/29;
median ceiling 1.0000 for both.

**So the click-seeded pipeline beats the click-free blob detector on 5 of 6 decision-grade
ROIs, and ties it on ceiling.** Sections 5 and 5b say the opposite ("still does not clear the
click-free bar", "loses to `blob_native` on 4 of 5 ROIs") -- those statements are true only of
the *single-template* configuration, which is not the configuration this log recommends.

Two further things this table settles that earlier sections got backwards:

* **Ranking augmented `TM_CCOEFF` by its own match score is as good as chromatin ranking**
  (125 vs 126), so the "the method advantage vanishes under chromatin ranking" finding in
  section 5b is a property of the *single-template* run, not of the method.
* **The seed-variance penalty largely disappears too**: augmented chromatin-ranked cv is
  0.147 against the single-template score-ranked 0.307.

Why it was missed is not a coding error. Experiment 3's scope was set to "does augmentation
help?" rather than "does augmentation change the answer to the question this whole log exists
to ask?", and that scoping followed from having already accepted the blob detector as the
standard. See section 9.

## 5f. Why TM loses on 300/301 -- neither loss is a matching failure

`results/tm_miss_300_301.csv`. Augmented `TM_CCOEFF`, all 5 seeds on both ROIs: **9 misses in
total**, out of 180 and 217 mitotic figures. For each, three separating questions were asked --
was there any local maximum near it before the floor; did the floor cut it; did NMS remove it.
Two causes, cleanly separated:

**Cause 1 -- one annotation outside the searchable region (4 of 9, all the same figure).**
`300.tiff` ann 14512 sits at (191, 4831) in a 6447x4835 image: **3 px from the bottom edge**. A
51 px template needs its centre 25 px clear of every edge, so the reachable region stops at
y = 4809 and this annotation is 22 px beyond it. `n_localmax_within_radius = 0` for every seed
-- there is no peak near it at any score, before any threshold. The template can only ever see
a truncated sliver of the object. `nucleus_blobs` reaches it because connected components have
no border margin. **This is a template border constraint, not a matching failure**, and it is
what the 300.tiff "blob win" mostly is.

**Cause 2 -- NMS removed a qualifying peak (5 of 9).** Every one had `above_floor = True` and at
least one local maximum inside the match radius, but after suppression the nearest survivor sat
at 30.1, 30.6, 30.8, 33.0 and 34.8 px -- just outside the 29.6 px scoring tolerance. These are
not faint objects: chromatin density 0.081-0.137, and **3 of the 5 unanimous**.

### The suppression radius sweep

`results/tm_nms_radius_sweep.csv`. Suppression radius varied, **scoring radius held at
29.6 px**, chromatin-ranked:

| ROI | NMS radius | pool | ceiling | missed | `read_50` |
|---|---:|---:|---:|---:|---:|
| 301 | 29.6 (current) | 16929 | 0.9954 | 4 | **198** |
| 301 | 24.0 | 23217 | **1.0000** | **0** | 266 |
| 301 | 20.0 | 29120 | **1.0000** | **0** | 319 |
| 301 | 12.0 | 42864 | **1.0000** | **0** | 457 |
| 300 | 29.6 (current) | 15363 | 0.9944 | 5 | **135** |
| 300 | 24.0 | 17595 | 0.9944 | 4 | 140 |
| 300 | 12.0 | 20332 | 0.9944 | 4 | 149 |

**On 301.tiff, dropping the suppression radius to 24 px reaches full recall on every seed** --
ceiling 1.0000, matching the blob detector -- for 37% more candidates and `read_50` 198 -> 266.
On 300.tiff it recovers only the one non-border miss and the ceiling stays pinned at 0.9944
(= 179/180, the border annotation), confirming that cause 1 is immune to this knob.

So of the two ROIs written up as blob wins: **301 is a tie once the suppression radius is
loosened**, and **300 is a tie except for a single annotation 3 px from the image edge**.
Neither is evidence that the blob detector recognises mitotic figures better.

**Do not read this as a recommendation to set the radius to 24.** The sweep covers only the two
ROIs where TM was losing -- a deliberately biased sample -- and it buys the tail by paying at
the head, which is the wrong trade on ROIs where TM already wins. The honest statement is that
the 300/301 deficit is a suppression-radius operating point plus a border margin, both
addressable, and neither a property of `TM_CCOEFF`.

## 6. The negative result section 6 pre-registered, and which half of it landed

Registered before the run: *"if the disc control matches the winning method, the finding is
that unnormalised TM on hematoxylin is a darkness detector the click does not steer."*

**That did not happen, and something adjacent did.** The disc control does not match
`TM_CCOEFF` -- the click is worth 8.1x at `read_50`, so the template is genuinely steering the
search and this is not a blob detector wearing a template. But the *look-alike* control does
match it (180 against 189), so what the click contributes is "a dark, textured, nucleus-sized
object", not "a mitotic figure". The click matters; its mitotic identity does not.

And the win is narrower than the exp1 headline suggests. Ranked by chromatin density -- which
is what any sensible pipeline would do, and what `7c3af93` already established -- all six
variants land between `read_50` 128 and 142, and `TM_CCOEFF` beats the incumbent by a median of
one candidate at p = 0.20. **The 8.6x is a fact about the TM score as a ranking key, not about
the candidate sets.** What survives the chromatin re-rank is the ceiling, and there the honest
number is that no variant solves goal (1) alone: 18 of 29 cells for the best configuration
against 19 for a detector that never sees a click.

**And under the configuration actually recommended, the click is worth close to nothing.**
The 8.1x is measured with the TM score as the ranking key. Rank by chromatin instead and the
click-free disc ties or beats the real click on 3 of 5 decision-grade ROIs -- `disc_ccoeff_od`
worst-of-5 reads 39 / 577 / 132 / 129 / 177 against `tm_ccoeff_od`'s 37 / 595 / 130 / 139 /
272. So in the pipeline this log recommends shipping, the pathologist's click is buying the
candidate *set* (where augmentation and the method both matter, section 5c) and almost nothing
about the *order*. That is the number the product should budget against, not the 8.1x.

Both of these belong in the log as plainly as the 8.6x does.

## 7. Reproducing

    python tm_variant_sweep.py --gates    # 25 gates, must print 0 failed
    python tm_variant_sweep.py --exp1     # ~26 min
    python tm_variant_sweep.py --exp2     # ~56 min
    python tm_variant_sweep.py --exp3 ccoeff ccorr   # ~20 min; NEEDS exp2's CSV -- see below
    python tm_variant_report.py tm_ccoeff_headtohead --pair tm_ccoeff tm_ccoeff_normed
    python tm_variant_report.py tm_variant_sweep

Section 5c's table is a *paired difference against exp2's single-template rows*, so
`--exp3` alone does not reproduce it: it writes `results/tm_variant_stage_b.csv` and has
nothing to difference against unless `results/tm_variant_sweep.csv` already exists. Run exp2
first, or the "189 -> 136, p = 0.0002" row cannot be regenerated from these commands.

Read `*_verification.csv` before any table, then `coverage_frac` and `n_pool`, then the tail.
A `ceiling` without its matched `n_pool` is not a result; a `read_*` without `n_unreached` and
`cv` is not a result.

### One defect found and fixed mid-run

The first exp1 pass failed `verify_shortcut` on 350.tiff: 2 of 8747 coordinates disagreed
between filtering the deep pool at z=2 and re-extracting at z=2. Cause was
`np.argsort(scores)[::-1]` -- quicksort is unstable, so two peaks with byte-identical scores
(0.5803979, 1 px apart) swapped depending on how many peaks were in the array. `extract_peaks`
now orders by a global key (score desc, then x, then y) and `nms_by_distance` uses a stable
sort, which makes the ordering independent of list length and the one-match-many-z economy
exact rather than exact-up-to-ties. That run was discarded and everything re-run on the fixed
code; the regression-tie gate confirms `TM_CCOEFF_NORMED`'s detection list is unchanged.

## 8. What independent audit found after the gates passed

The 4257 checks are all gates I wrote, and a gate only catches what its author anticipated.
Three independent audits -- one on code correctness, one re-deriving every number from the CSVs
without reading this log first, one on design validity -- found **fifteen** defects, including
one that invalidates an entire arm. Each is corrected in place above; this is the index.

**Invalidated a result outright:**

A. **`TM_SQDIFF_NORMED`'s pool is a tie-break artefact and the method is withdrawn.** OpenCV's
   [0,1] clamp saturates 83% of the map, collapsing the MAD to exactly zero, so the extraction
   floor sits on the clamp and admits ~28M tied maxima; the `max_peaks` cap then keeps a
   left-edge strip. 228 rows on 10 of 13 ROIs. It also set `n_matched` in 20 of 58 cells, so
   every `@matched` figure in `tm_variant_sweep.csv` is partly truncated by it. exp1 and exp3
   are clean and carry the headline ceiling numbers. In an earlier draft this method was
   reported as 4th of 6, and its chromatin-ranked variant as the *best* arm in the whole grid.

B. **The invariant that should have caught A was guarding the wrong quantity.**
   `invariants.check_no_cap` was applied to the post-NMS list length; the cap that actually
   binds is `max_peaks` *inside* `extract_peaks`. 228 capped rows passed all 4257 checks. This
   is defect M4 from `2026-09-01-premise-test-audit.md` -- "invariant #1 defeated in the code
   meant to enforce it" -- repeating in new code. `match_pool` now raises on both a zero MAD
   and a bound cap.

C. **`verify_shortcut` passes vacuously when the MAD is zero** (the cut equals the floor, so it
   compares against an unfiltered re-extraction), and runs only at seed 0 and one z level.

**Would have changed a conclusion:**

1. **Ceiling was reported at full pool, against section 4's own pre-committed rule.** The
   headline "18 of 29 against 19" is full-pool; matched-pool is 12 against 14. Direction
   unchanged, gap wider. (Section 5, 5c.)
2. **The 245.tiff "geometric ceiling" check was invalid** -- not a bound at all, and computed on
   a different denominator than what it was compared against. Retracted in section 5.
3. **The look-alike control compared a 29-cell median against a 6-cell median.** The corrected
   like-for-like statement is *stronger* than the original claim. (Section 5b.)
4. **Ceiling and `read_95` significance do not survive clustering.** 29 cells are 6 ROIs x 5
   seeds; at ROI level the depth result holds (p = 0.031) and those two do not. (Sections 5, 5d.)
5. **The disc ratio was a ratio of unpaired medians**: 8.1x becomes 6.9x paired, the "1.6-46.6x"
   range was built on a ROI the adjacent table excludes, and `TM_SQDIFF_NORMED`'s ratio flips
   sign when paired. (Section 5b.)

**Code bugs:**

6. **`blob_native@matched`'s coverage cache key omitted the seed index**, so one seed's coverage
   was reported for all five -- on the column used to adjudicate P4. Fixed; error is 0.0000 on
   the ROI that carries the claim. (Section 5.)
7. **`paired_delta` would silently mis-pair a `lookalike_*` arm**, which shares `seed_index=0`
   with a different `seed_ann_id`. Never invoked that way; now raises.

**Diagnostic columns and reporting:**

D. **Diagnostic columns are cross-assigned by method.** `got["method"] == m` also matches
   `disc_*` and `pipeline_as_shipped`, so those rows carry the *seeded* pool's `extract_floor`
   and `n_peaks`; `pipeline_as_shipped` additionally reports `channel = hematoxylin_od` and
   `base_size = 51` when it actually ran RGB with a per-seed tightened size. `map_median` /
   `mad_scale` are one method's values broadcast to all rows. No conclusion rests on these
   columns, but do not read them from the CSVs.
E. **`--exp3`'s emitted paired file answers a different question than exp3 asks** -- it is a
   ccoeff-vs-ccorr contrast *within* the augmented run, not augmented-vs-single. The rotation
   contrast reported in section 5c was computed separately and is correct; the file is not it.
F. **`tm_variant_report.py` gaps**: `spread_table` ignores `--all-rois`; the look-alike arms
   (`n_seeds == 1`) are filtered out of three of six tables; `paired_summary` reports n=29 for a
   test scipy runs on 13 non-zero pairs.

**Misstatements:**

8. **P1 was half-confirmed, stated as confirmed** -- `TM_SQDIFF_NORMED` ranks 4th of 6, not
   bottom. (Section 5b.)
9a. **Two docstring claims are false**: `disc_template` says only `TM_SQDIFF` depends on
   template scale (`TM_SQDIFF_NORMED` does too), and the module docstring's "normed vs unnormed
   *is* the can-it-see-magnitude axis" does not hold for `TM_SQDIFF_NORMED`, whose score for a
   window `a*T` is `(1-a)^2/|a|` -- minimised at `a = 1`, so neither invariant nor monotone in
   gain. Section 2's table was right; the prose generalising it was not.
9. **Scope, pre-registration wording, unreached counts, the clamp gate's undisclosed
   restriction to decision-grade ROIs, two wrong table cells, and "seedless arms are
   deterministic by construction"** (contradicted by this log's own cv table). All corrected.

**Cleared on inspection:** the pairing is genuine (`seed_ann_id` identical across all three
experiments on all 58 cells); the look-alike ground-truth difference is immaterial (the removed
annotation is 122-546 px from any mitosis against a 30 px radius); `TM_CCOEFF_NORMED`'s pools
are consistently *larger* than `TM_CCOEFF`'s, so the full-pool ceiling comparison was
conservative rather than flattering; and the `read_100` conditioning warning, the
`pipeline_as_shipped` confound warning and the coverage caveat were all correctly scoped.

**Not resolved:** seed 0 is the deepest of five draws under most methods on 3 of 6 ROIs
(binomial p = 0.099), and cannot be retested here because all three experiments share the same
draws. Re-drawing with different RNG streams is the check.


## 9. Where the framing came from

This log was written with the seedless blob detector positioned as the incumbent and template
matching as the challenger. That framing was not neutral and it cost a result (section 5e).
Its sources, traced:

1. **An unaudited file was given load-bearing status.** `results/tm_vs_blob_depth.csv` --
   the origin of the "4.3-216x" figure this log calls "the bar" -- is cited in **no** prior
   research log. It was produced by a script with no writeup and never went through the
   adversarial review every other result in this repo did. It was labelled "provisional" here
   and then used to frame the entire experiment anyway.
2. **The bar was built into the experimental design.** `blob_native` appears in every
   comparison as the standard TM must clear. An equally valid design would have treated it as
   one baseline among several and asked what the *click* buys, which is the actual product
   question.
3. **Worst-of-N systematically penalises the seeded method.** The blob detector's cv across
   seeds is 0.018 because it cannot see the click at all; TM's is 0.15-0.37 because it can.
   Reporting the worst seed as the headline is defensible for product risk, but it is
   structurally unfavourable to the method under test and that was never stated.
4. **Failures were reported asymmetrically.** TM missing 1 mitotic figure of 218 was written up
   as "blob wins"; the blob detector missing 2 of 18 on 201.tiff and 4 of 116 on 246.tiff was
   mentioned only where TM happened to win. Normalised, TM's worst per-ROI miss rate is 0.6%
   and the blob detector's is 11%.
5. **The union recommendation imported the weaker method's brittleness.** Combining TM with
   watershed-split blobs pulls in `min_area=80` / `max_area=4000` / `tile=512` / `gray_max=220`
   -- constants never tuned per domain, and `results/miss_attribution.csv` attributes 3 of 11
   documented misses directly to `max_area`. There is no ground truth at inference time to tune
   them against. The click, by contrast, supplies template size and stain intensity for free.

The correction is not "template matching wins". The blob detector is still ahead on 301.tiff
and on ceiling for both mast-cell ROIs. It is that the comparison was set up so that parity
read as defeat.
