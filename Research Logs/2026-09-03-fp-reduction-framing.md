# "Cover every mitosis" and "a short list" are one metric, not two

Date: 2026-09-03
Question asked: *"My postdoc wants the threshold raised so the candidate list covers the entire
mitotic population in any domain and any ROI. That takes thousands of candidates. My secondary
goal is to cut the false positives, and I don't think prior work has done anything about that.
What angle am I missing? My postdoc says maybe I should fix a small candidate budget and measure
how many mitoses fall inside it instead."*
Status: analysis only. No experiment was run to write this. Every number is either re-derived
by me from artifacts already in `results/` (*[verified]*), quoted from a prior log without
re-checking (*[unverified]*), or measured by me here from `databases/MIDOG++.json` (*[new]*).
Two nulls in §6 and one enrichment probe in §5 are new measurements. §3b reports an
adversarial audit of the `fcos_midogpp_r18` arm run on 2026-09-04, whose decisive finding I
re-derived independently with my own matcher.

---

## 1. The thing that makes this feel stuck

The two goals are being asked of one ranked list, and **on a ranked list "push true positives
up" and "remove false positives" are the same operation.** A ranked list has no false positives
in it until you cut it; you only ever choose a depth. Promoting one true positive past one
false positive *is* removing that false positive from every cut between them. So there is no
second lever to pull, and every experiment in this repo that tried to pull it — six
operating-point moves on the correlation search, the channel comparison, the NMS ordering, the
bbox tightening, the threshold sweeps — has come back neutral for a reason that is structural,
not empirical.

Four things, and as far as I can tell only four, make FP-removal a genuinely different lever
from ranking:

| lever | what makes it different | status in this repo |
|---|---|---|
| **(a) information the ranker cannot see** | a score computed from something outside the candidate patch | the click currently supplies none of it; reader rejections would (§4) |
| **(b) set-level criteria** | diversity / marginal-coverage selection, which is not expressible as a per-candidate score | named once in `2026-09-02-next-steps-plan.md` §4c, never run |
| **(c) scores that change while the reader works** | relevance feedback, online refit | **never run, on any path** |
| **(d) changing the unit of decision** | one rejection kills many candidates | **never considered; every metric here assumes 1 candidate = 1 decision** (§5) |

Everything else — better channel, better similarity, better threshold, better encoder — is
lever zero, i.e. ranking, and its ceiling is measured in §3.

## 2. Prior work *has* attacked the FP count. It attacked the wrong 5% of it.

The premise that "nobody has done anything about false positives" is not right, but the work
that was done targeted a small slice of the problem:

* `results/tail_sensitivity.csv` *[verified]* — at 100% sensitivity on 301.tiff the reader
  passes **105 annotated look-alikes and 7,427 unannotated nuclei**. Across every ROI and
  every sensitivity level, pathologist-marked mimickers are **1–8% of the false-positive mass**.
* So the chromatin work, the morphometric Bhattacharyya ranking, the texture plan, the
  look-alike AUC controls and the "hard negative" framing all target ~5% of the volume. The
  other 95% is ordinary nuclei. `2026-09-01-tail-sensitivity.md` §2 states this and the point
  has not propagated into any subsequent plan.
* The one direct FP-filter experiment that was run — a 22-feature logistic regression,
  leave-one-domain-out, `results/fp_filter_domain.csv` — is **2 wins / 2 losses** on the four
  ROIs where 95% sensitivity resolves as a quantile *[unverified, that log's §5]*. Worth
  shipping for the worst case; not a route to a 10x.

## 3. The bound: what the last 20% of sensitivity costs, at the ceiling

This is the number that should govern the decision, and it is sitting unplotted in
`results/click_rank_metrics.csv`. 301.tiff, 20,369 candidates, 217 mitoses, identical candidate
set across arms, 1 click, median over 5 seeds *[verified]*:

| ranker | read_50 | read_80 | **read_100** | 80% → 100% |
|---|---:|---:|---:|---:|
| `blob_native` (no click, free, 7-11 s) | 167 | 480 | **7,749** | **16.1x** |
| `chromatin` (the committed re-ranker) | 192 | 1,576 | 19,933 | 12.6x |
| `fcos_midogpp_r18_cl2n_snap` (**contaminated**) | 119 | 251 | **1,830** | **7.3x** |
| **`||e-mu||`, contaminated encoder, NO CLICK** *[new, §3b]* | **116** | 262 | **2,703** | **23.3x** |

`blob_native`'s 7,749 reproduces `results/tail_sensitivity.csv`'s independently computed 7,750
for `blob_score` on the same ROI. An adversarial audit (2026-09-04) confirmed these are genuinely
separate computations: `tail_sensitivity.py` re-reads the ROI, re-runs `nucleus_blobs`, and uses
the full GT (218) and full pool (20,370), where `click_rank_stage3.py` drops the seed annotation
and the snapped candidate (217 / 20,369). The snapped candidate sits at blob-rank 42, above depth
7,750, so every later depth shifts down by exactly one. *[verified]*

`fcos_midogpp_r18` is the ResNet-18 body of a detector **trained on MIDOG++, which contains
301.tiff with labels**. It is not a product candidate. It is the best case anything could do
here, and **even it needs ~1,830 candidates for the last mitosis on 301.tiff, against 251 for the
first 80%.**

Two consequences:

1. **"Cover the entire mitotic population at a small candidate count" is not reachable by any
   representation.** The gap between 251 and 1,830 is not a modelling deficit that a better
   encoder closes -- it is what a memorising oracle already costs.
2. On **246.tiff, `read_100` is NaN for every arm and every seed** *[verified]* -- the
   `nucleus_blobs` proposal ceiling is 0.965 (`full_list_recall` identical across all 23 arms),
   so 100% is undefined there before any ranking question arises. Raising a *matching* threshold
   cannot fix a *proposal* ceiling.

The shape underneath is a convex cost curve. `blob_score` on 301.tiff, marginal candidates read
per additional mitosis found *[verified, derived from `tail_sensitivity.csv`]*:

| band | extra candidates | extra mitoses | cost per mitosis |
|---|---:|---:|---:|
| 0 -> 50% | 167 | 109 | 1.5 |
| 50 -> 90% | 553 | 87 | 6.4 |
| 90 -> 95% | 708 | 11 | 64 |
| 95 -> 98% | 2,879 | 7 | 411 |
| 98 -> 100% | 3,443 | 4 | 861 |

**The pathologist is being asked to pay 861 candidate decisions for the last mitosis, which
moves a count of 218 by 0.46%.** No threshold setting makes that trade good.

## 3b. Audit of the one arm all of this leans on: the numbers are right, the story was wrong

An adversarial audit of `fcos_midogpp_r18` (2026-09-04) was commissioned because §3 rests on it.
Findings, most important first. I re-derived the decisive one myself with an independently
written greedy matcher; it reproduces exactly.

**The weights and the plumbing are clean** *[verified by audit]*. `FCOS_18.ckpt` is a genuine
Lightning checkpoint (epoch 77, 174 keys). `click_rank_embed.py:126-127` keeps only
`model.backbone.body.*` -- both detector heads and the FPN are excluded, so the "embedding" is a
512-d global-average-pooled ResNet-18 trunk feature, **not** a mitosis logit. The load is guarded
(`:132-133`): `missing == ['fc.weight','fc.bias']`, `unexpected == []`, so a silent random-init or
ImageNet fallback is impossible. Every tensor differs from `resnet18-f37072fd.pth` (mean relative
conv difference 0.057 at the stem rising to 0.340 at layer4) -- a fully fine-tuned trunk.
Evaluation is clean too: `gt_eval` is built once per (ROI, pass, seed) and handed to every arm;
the snapped candidate is removed from the scored list for all arms before scoring and its
annotation is out of `gt_eval`; candidate/embedding alignment (asserted in the code only by length
equality) was checked by re-cutting 15 crops from `images/301.tiff` -- byte-exact. Every number
quoted in §3 and §7 reproduces to the digit.

**But the headline result is reproducible with no click at all.** *[new, verified by me
independently]* Rank candidates by `||e_i - mu||`, the L2 distance of the raw fcos embedding from
that ROI's *mean* embedding. No query, no click, no labels -- and the ROI mean is exactly the
click-independent ingest-time quantity the log's own arms table (line 73) and
`click_rank_stage3.py:134-136` already argue is free:

| | 301 read_50 | 301 read_80 | 301 read_100 | 246 read_50 | 246 read_80 |
|---|---|---|---|---|---|
| **no click, `||e-mu||`** | **116 / 116 / 116** | 262 | 2,703 | **63 / 63 / 63** | 152 |
| 1 click, `cl2n_snap` | 114 / 119 / **220** | 251 | 1,830 | 60 / 64 / **80** | 120 |
| 3 clicks (prototype) *[audit]* | 113 / 123 / 135 | 223 | 1,678 | 59 / 61 / 69 | 116 |
| `blob_native` (the log's bar) | 167 | 480 | 7,749 | 193 | 577 |

(min / median / worst over 5 seeds where three numbers are given.)

Four things follow, and they change what the click-ranking log established:

1. **The pre-registered gate inverts.** That gate was worst-of-5-seeds `read_50`. The no-click
   ranker passes on 301 (116 on every seed against `blob_native`'s 167); the 1-click arm
   **fails** (worst seed 220). "Fails at one click, passes at three" is really "fails at one
   click, and at three clicks still does not beat doing nothing with the click" (3-click worst
   135 against no-click 116).
2. **The click's value is real but lives entirely in the tail.** At `read_50` it contributes
   nothing (119 vs 116) and costs stability (worst seed 220 vs a deterministic 116). At
   `read_80` it is 4% better on 301 and 21% on 246. At `read_100` it is **32% better** (1,830 vs
   2,703). **The project optimised `read_50` for four months -- the one depth at which the click
   is worthless.** This is the strongest reason yet to move to the §7 metric.
3. **The provenance control survives but argues the opposite of what the log claims.** The
   18-172x mitotic-vs-look-alike gap re-derives exactly, and the candidate set and `gt_eval` are
   genuinely identical across the three provenances (`click_rank_provenance.py:75-77`). But
   against the correct null the ordering is: **no click 116 < mitotic click 119 << look-alike
   click 1,294 << random click 11,256.** Query-dependence is not click *value*. The control shows
   a bad click destroys a good list, not that a good click builds one.
4. **The log inverted its own falsification test.** Lines 91-94 include this arm precisely so
   that "if even a representation built for this exact task does not make the click pay, the
   negative result is about the click and not about encoder choice." That is what happened. The
   log read it as the click paying.

**Two caveats, both load-bearing.** (i) The sign is one label-informed bit -- ranking by
*proximity* to the mean gives 20,149. "Rank nuclei by how unusual they are for this slide" is an
a-priori-natural direction rather than a fitted parameter, but it is a choice and should be
pre-registered, not chosen after seeing the answer. (ii) **This is specific to the contaminated
encoder.** In every frozen encoder the same click-free statistic is far *worse* than the click
(301 seed 0: `lunit_16um` click 1,067 vs no-click 5,656; `kaiko` 3,352 vs 8,992;
`resnet18_in1k` 5,640 vs 6,393). So "unusual = mitotic" holds only in a space already trained to
find mitoses. It is not a training-free shortcut; it is a statement about what that training
bought.

**The right characterisation** is therefore *not* "a contaminated upper bound on one-click
retrieval". It is **a valid upper bound on what a MIDOG++-contaminated mitosis-trained encoder
achieves on these two ROIs -- which is a click-free quantity.** §3's bound on the cost of the
last 20% stands, and is now broader than stated: it bounds *any* method built on this
representation, clicked or not.

**One defect in the log's own verification, worth recording.** Its header (line 13) cites
`coverage_frac` having `nunique()==1` across arms as evidence the candidate sets match. That is a
caching tautology -- `compare.py:188-193` caches coverage by `arm.coverage_key`, and
`click_rank_stage3.py` gives every arm the same key. `n_detections`, `full_list_recall`,
`n_gt_mitotic` and `n_gt_lookalike` do carry the weight; `coverage_frac` carries none. The same
header promises alignment checks that do not exist in either script (they were run during the
audit and pass).

## 4. The lever nobody has pulled: the reader's own rejections

Lever (c) above, and it is the only one of the four that is both cheap and matched to what the
product already does.

Today: one click → a static list of 20,000 → the reader works down it. By candidate 1,428 on
301.tiff the reader has generated **1,300+ labelled hard negatives on this exact slide, this
exact stain, this exact tumour** — and the system throws all of them away.

Why this is the right lever for *this* FP distribution specifically: §2 says the volume is
ordinary nuclei, and ordinary nuclei are highly redundant. A distribution whose mass sits in a
few dense modes is exactly the distribution where a handful of labelled negatives collapses a
large slice of the list. It is the opposite of the mimicker problem, where each negative is
individually hard and teaches you little.

It is also the one thing a click can supply that no seedless method and no pretrained encoder
structurally can. Note what the disc control already established
(`2026-09-02-next-steps-plan.md` §2 item 4, *[verified there]*): per-slide **stain** calibration
is available with **no click at all**, from the image's own OD percentiles. That removes one of
the two standing justifications for the click. Per-slide *decision boundary*, learned from this
reader's accept/reject stream, is the one that survives.

`2026-09-02-next-steps-plan.md` §4c calls negative clicks "the most under-exploited signal in
the repo" and then ranks the experiment fifth. It has still never been run on any path.

## 5. The unasked question: is the false-positive mass redundant?

Every metric in this project — `read_50`, `read_100`, `recall_at_budget`, depth-to-95% — counts
**candidates** and assumes one candidate costs one pathologist decision. Nobody has checked
whether that denominator is real.

The question, stated so it can be measured: **how many distinct visual modes do 301.tiff's 7,427
unannotated false positives fall into?** If the answer is tens rather than thousands, then the
same ranked list, presented as a grid of visually clustered thumbnails with reject-the-group,
becomes tens of decisions with **no change to any model**. That is the only candidate on this
page that plausibly delivers the order of magnitude being asked for, and it is testable from
artifacts already cached (`~/.cache/annotatedx_click_rank/`, per the click-ranking log).

`AnnotateAnyCell` (bioRxiv 2025.11.02.686114) is exactly this UX and is already cited in
`2026-09-01-one-click-retrieval-literature.md` §B4 — as an annotation tool, never connected to
the burden metric.

## 6. Two things I measured that close off plausible directions

### (a) There is no spatial prior to exploit. *[new]*

A natural idea — mitoses cluster in proliferative hotspots, so condition the candidate list on
distance from the click — is dead on this data. Measured from `databases/MIDOG++.json` over the
10 ROIs with ≥ 15 mitotic annotations:

| ROI | n | mean NN dist | expected under CSR | ratio |
|---|---:|---:|---:|---:|
| 094 | 82 | 313.4 | 340.8 | 0.920 |
| 201 | 18 | 673.9 | 671.3 | 1.004 |
| 202 | 16 | 707.1 | 712.0 | 0.993 |
| 245 | 90 | 260.0 | 300.2 | 0.866 |
| 246 | 116 | 250.5 | 264.4 | 0.947 |
| 300 | 181 | 191.9 | 207.5 | 0.925 |
| 301 | 218 | 188.6 | 189.1 | 0.998 |
| 402 | 105 | 330.6 | 304.9 | 1.084 |
| 459 | 131 | 234.0 | 248.8 | 0.940 |
| 548 | 239 | 209.1 | 202.1 | 1.035 |

Clark-Evans ratios of 0.87–1.08 — indistinguishable from a random point pattern. The
click-conditioned form is more direct: given a click on a random mitosis, the share of the
*other* mitoses within radius D against the share of ROI area inside that disc, averaged over
all possible clicks (301/246/300/245, D = 500–3,000 px): enrichment **0.96–1.09x** at every
radius on 301/246/300, and 1.31x only at 245.tiff's smallest radius where the count is tiny.

**Mechanism, so this reads as explained rather than surprising:** MIDOG++ ROIs were *selected*
as mitotic hotspots. Conditioning on hotspot membership removes the clustering, by construction.
Any spatial prior would have to operate at whole-slide hotspot-selection scale, which is a
different problem (`2026-09-01-tail-sensitivity.md` §7 says so).

### (b) A quarter of the target is a 2-of-3 majority vote, and the tail may be enriched in it. *[new, and only suggestive]*

`databases/MIDOG++.json` carries per-annotation votes. Over all **11,937** category-1 (mitotic)
annotations, **74.7% are unanimous**. Over the 605 in the seven tail-sensitivity ROIs, 78.2%.

Over the 32 deepest-ranked captured mitoses in `results/tail_object_audit.csv` (the 8 worst per
ROI, i.e. the objects that set the reading budget), the unanimous rate is **65.6%** — 0.444 on
300.tiff against a 0.724 base, 0.615 on 246.tiff against 0.862.

**Scope this tightly.** n = 32, top-8-per-ROI, no significance test — this is a cheap hypothesis
worth testing properly on the full rank distribution, not a finding. And carry `dataset.py`'s
own caveat: for category 1 the only label multisets in the dataset are `(1,1)` ×8,917 and
`(1,1,2)` ×3,020, so `unanimous` is *exactly* `n_votes == 2`, i.e. "did this need a third
reader". It is a difficulty proxy, not an independent measure of label quality.

If it holds, it changes the spec: "100% of annotated mitoses" partly means "100% of objects a
2-of-3 majority called mitotic", and the objects that cost 861 candidates each are
disproportionately the ones the pathologists themselves argued about.

## 7. The postdoc is right, and the strongest reason is stability

Fixing a budget K and measuring recall@K is the correct move, for four reasons, of which
convenience is the weakest:

1. **It is the only metric on which the two goals are commensurable.** Sensitivity and list
   length become one curve instead of two competing targets.
2. **It is bounded and stable; depth-to-target is neither.** Depth is hostage to individual
   objects: 246.tiff goes 98% → 1,966 candidates, 99% → 17,321 — **one annotation for 8.8x**
   *[verified]*. And it is undefined wherever the ceiling is below the target (246.tiff's
   `read_100` is NaN on every arm). recall@K is defined everywhere and moves smoothly.
3. **It is the decision the product actually makes.** The tool must choose a list length before
   the reader starts. Depth-to-100% is not a knob anyone can set.
4. It is already the plan of record — `2026-09-02-next-steps-plan.md` M6 adopts
   `recall_at_budget` as the standard and notes it "is already computed in every CSV and unused
   in every log". The curves exist and have never been plotted.

Here they are, from `results/click_rank_metrics.csv`, 1 click, identical candidate sets
*[verified]*:

**301.tiff** (217 mitoses, 20,369 candidates, ceiling 1.000) — median / worst of 5 seeds

| arm | @500 | @1000 | @2000 | @5000 |
|---|---|---|---|---|
| `blob_native` (no click) | 0.802 | 0.931 | 0.968 | 0.986 |
| `chromatin` | 0.631 | 0.760 | 0.843 | 0.908 |
| `fcos` (contaminated) | 0.945 / 0.922 | 0.986 / 0.977 | 1.000 / 0.995 | 1.000 / 1.000 |

**246.tiff** (115 mitoses, 21,058 candidates, ceiling **0.965**)

| arm | @500 | @1000 | @2000 | @5000 |
|---|---|---|---|---|
| `blob_native` (no click) | 0.765 | 0.861 | 0.948 | 0.948 |
| `chromatin` | 0.678 | 0.774 | 0.878 | 0.922 |
| `fcos` (contaminated) | 0.948 | 0.965 | 0.965 | 0.965 |

**At K = 500, a free seedless blob detector already returns 76–80% of the mitoses, and the
best-case contaminated oracle returns 95%** — which on 246.tiff is 98% of everything the
proposal stage made reachable. The entire remaining argument is about the last 5%, and §3 prices
it at 400–900 candidates each.

**The form to commit to, because plain median recall@K hides the concern that motivated the
postdoc's request in the first place:** *the smallest K such that recall@K ≥ target on the
**worst** ROI and the **worst** seed, reported per domain.* That is what "any domain, any ROI"
means operationally, and it is the gate this repo already adopted for other purposes
(`2026-09-01-one-click-retrieval-literature.md` §7 protocol requirement 2).

## 8. What is missing from the toolchain, as opposed to the science

The repo has a candidate generator, several rankers, and a metric. It has no **calibration
layer**, and the word "conformal" does not appear anywhere in it.

The question "how deep must the reader go on a *new* slide to be confident of 90% sensitivity"
is a distribution-free calibration question with a standard answer — conformal risk control /
Learn-then-Test / risk-controlling prediction sets — that converts a ranked list into a set with
a guarantee: *choose K such that recall ≥ 0.9 with probability ≥ 0.95 on an unseen ROI*, using
the other ROIs as the calibration set.

Two honest caveats, in order:

* **This does not reduce false positives.** If the ranker needs 7,749, conformal tells you 7,749
  with a confidence statement attached. Do not mistake it for an answer to the FP question.
* It is what makes the fixed-K framing *shippable* rather than arbitrary, it is the right home
  for the "click as per-slide threshold setter" idea (`2026-09-01-tail-sensitivity.md` §8 Step
  3a, still unmeasured), and per-domain calibration is exactly how "any domain, any ROI" stops
  being a hope and becomes a measured coverage claim.

## 9. If I were ordering the next month

1. **Restate the spec as recall@K, worst ROI, worst seed, per domain** (§7). One afternoon of
   plotting from CSVs that already exist. Everything below is unreadable without it.
2. **Measure FP redundancy** (§5). The only candidate that plausibly delivers 10x, testable from
   cached embeddings, and it changes the product rather than the model.
3. **Run one relevance-feedback experiment** (§4). Simulate a reader accepting/rejecting down
   the list, refit the 22-feature logistic regression every N decisions, and re-measure
   recall@K. The model already exists and trains in seconds.
4. **Test the tail-unanimity hypothesis properly** (§6b) before committing to any 100% target.
5. **Never run another representation experiment without a click-free null in the same
   representation** (§3b). The `||e-mu||` arm costs one line over cached embeddings and it is what
   separates "the click works" from "the encoder works". The click-ranking log's recommended
   follow-up -- retrain leaving these ROIs out -- must carry it, or it will re-measure the
   encoder and report it as the click.
6. **Do not** build a spatial prior (§6a), tune the correlation search further, or benchmark
   another frozen pathology foundation model (`2026-09-01-click-ranking-experiment.md` refutes
   that path at 2.8-37x worse than doing nothing with the click).

## 10. Reproducing the new measurements in §6

Both were computed inline from `databases/MIDOG++.json` and `results/tail_object_audit.csv`;
neither is scripted in the repo yet. The spatial null is a Clark-Evans nearest-neighbour ratio
(observed mean NN distance against `0.5 * sqrt(A/n)`) plus a Monte-Carlo disc-area comparison at
20,000 sample points per ROI. The unanimity probe is a join of `tail_object_audit.csv`'s
`ann_id` against the vote multisets in the JSON. If either is to be quoted, script it and stamp
the commit — `2026-09-02-next-steps-plan.md` M1.
