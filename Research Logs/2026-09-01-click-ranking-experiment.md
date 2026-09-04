# Does the click work as a *ranker*? The §7 experiment

Date: 2026-09-01
Plan: `2026-09-01-one-click-retrieval-literature.md` §7
Code: `click_rank_stage1.py`, `click_rank_embed.py`, `click_rank_stage3.py`
Results: `results/click_rank_metrics.csv`, `results/click_rank_auc.csv`
Intermediates: `~/.cache/annotatedx_click_rank/` (crops and embeddings, ~2.3 GB, regenerable)

Run: 2 ROIs x 5 seeds x 2 passes x 23 arms, plus a 3-provenance control. Adversarially
verified over three rounds by an independent re-implementation (own greedy matcher, no import of the scoring code): every
headline `read_50` reproduces exactly. Every arm in a comparison has an identical candidate
set -- `n_detections`, `coverage_frac` (0.92-0.94), `full_list_recall`, `n_gt_mitotic` and
`n_gt_lookalike` all have `nunique() == 1` across arms within each (ROI, pass, seed) -- and
the largest tie block anywhere is 6 rows (border-clamped duplicate crops). The `chromatin`
arm alone carries a 2.2% NaN ranking key (border windows), sorted last.

The 690 `compare.evaluate_arms` invariant records are **not** evidence of anything: neither
script sets `caps` or `nms_radius`, so every one of them is a vacuously-true `check_no_cap`.
The checks that carry weight are the ones in the paragraph above, plus the alignment and
self-hit checks below.

## The question

The premise test measured the pathologist's click through `TM_CCOEFF_NORMED`, which the
chromatin log proved is invariant to mean intensity -- the strongest measured discriminator
in this dataset. So it measured pixel template matching, not the click. This run asks the
question the premise test could not: **on one identical candidate set, does ranking by
similarity-to-the-click in a learned embedding beat ranking by chromatin density?**

Detection is not under test. `nucleus_blobs` generates the candidates for every arm, and it
captures 96.6-100% of both mitotic figures and annotated look-alikes. Only the sort key
differs between arms.

## Protocol

| | |
|---|---|
| ROIs | 301.tiff (lymphosarcoma, 217 mitotic / 107 look-alike captured), 246.tiff (mast cell, 111 / 123) -- the only two with enough of both classes |
| candidates | `nucleus_blobs`, 20,370 (301) and 21,059 (246), identical across arms within a seed |
| seeds | 5, `np.random.default_rng([seed_index, image_id])` into `seed_selection.pick_seed` -- the same convention as the premise test |
| passes | `1click` (the product interaction) and `3click` (SimpleShot prototype over 3 seeds, so "one click is too little signal" cannot explain a null) |
| crops | 32 um field, cut at each ROI's own mpp and resampled to 128 px, so a candidate is framed identically on ROIs of different resolution. 16 um variants are the exact centre crop |
| decision metric | `read_50` -- candidates a reader works through to reach 50% of the mitotic figures |
| gate | beat the best seedless arm on the **worst** of 5 seeds, on **both** ROIs |
| diagnostics | AUC vs annotated look-alikes, AUC vs unannotated candidates, per-seed spread |

Protocol requirements, each from a mistake this project already made once:

* `gt_eval` drops the seed's own annotation(s), and the *same* `gt_eval` goes to every arm,
  including the seedless ones.
* The clicked cell is dropped from the candidate list of every arm -- within
  `self_hit_radius`, plus the single nearest candidate whatever its distance -- so it
  cannot be returned as a discovery or occupy rank 1 at cosine 1.0.
* AUC labels come from one canonical bucketing of the candidate set rather than each arm's
  own ranking, because `bucket_detections` matches greedily in row order and would otherwise
  label the same candidate differently for different arms. That canonical order is
  `nucleus_blobs`'s output order -- which *is* the `blob_native` arm's ranking, since the
  detector returns score-sorted rows. The consequences are worked through under "AUC
  diagnostics" below; `read_50` is unaffected, being bucketed per arm.
* Read depth is reported per seed and gated on the worst one. The premise test found
  `read_50` varying 2.2x with which cell was clicked, and unreached on 2 of 5 seeds.

## Arms

Two seedless incumbents, and each encoder in three query constructions.

| arm | what it is |
|---|---|
| `chromatin` | mean of the darkest 10% of a 51 px window on the unclipped hematoxylin OD channel -- the incumbent ranker |
| `blob_native` | the Otsu component's own mean hematoxylin -- the other seedless ranker, and on this candidate set the harder bar |
| `cos_rawpixel_ccoeff` | cosine between mean-centred, L2-normalised grey-inverted crops. This *is* Pearson correlation, i.e. `TM_CCOEFF_NORMED` -- the incumbent representation, put on the same footing as the learned ones |
| `cos_<encoder>` | cosine to the click's embedding |
| `cos_<encoder>_cl2n` | SimpleShot CL2N: subtract the ROI's mean embedding, re-normalise, then cosine. The ROI mean is click-independent, so it is ingest-time information |
| `cos_<encoder>_cl2n_snap` | CL2N, with the query taken from the candidate the click lands on rather than a crop centred on the clicked pixel |

Encoders, all frozen, all run on CPU on the crop at its stored resolution (no resampling,
so no arm can win on interpolation):

| encoder | what it is | licence |
|---|---|---|
| `lunit_dino_vits16` | ViT-S/16, DINO self-supervised on TCGA pathology | Apache-2.0 (port) |
| `kaiko_vits16` | ViT-S/16, DINOv2 self-supervised on TCGA pathology | MIT |
| `resnet18_in1k` | ImageNet-supervised CNN -- the generic control | permissive |
| `fcos_midogpp_r18` | the ResNet-18 body of the MIT-licensed MIDOG++-trained FCOS detector (F1 0.737-0.753 on the MIDOG 2022 test set) | MIT |

DINOv2 ViT-S/14 was in the plan as an out-of-domain SSL control and was dropped: it is the
most expensive encoder here (patch-14 gives 81 tokens against 64) and the question it
answers -- do generic natural-image features work -- is already covered by
`resnet18_in1k`, with two in-domain encoders and a task-trained upper bound above it.

`fcos_midogpp_r18` is a **leaky upper bound**, not a product candidate: it was trained on
MIDOG++, which contains 301.tiff and 246.tiff with labels. It is here to answer a different
question -- if even a representation built for this exact task does not make the click pay,
the negative result is about the click and not about encoder choice.

Why 128 px rather than each checkpoint's native 224: what a ViT sees is set by
microns-per-token. These encoders were trained near 0.5 um/px, where a patch-16 token covers
~8 um. Feeding this 0.25 um/px crop at 224 px would make a token cover 2.3 um -- 3.5x too
zoomed -- against 4 um at 128 px, and costs 4x more. DINO and DINOv2 pretraining both used
96-98 px local crops, so 128 px is inside these models' experience.

## The answer, in one paragraph

**The click carries a great deal of rankable information — and no frozen off-the-shelf
encoder can extract it.** In the feature space of a network trained to detect mitotic
figures, clicking one mitotic figure cuts the reader's workload to **60-80 candidates on
246.tiff against the best seedless arm's 193** — a 2.4x reduction *on that ROI*; on
301.tiff a single click's worst seed loses to the seedless bar and only three clicks win it.
The same click on a pathologist-rejected look-alike costs 1,083-10,640, so it is the click's
*content* doing the work, not the encoder alone. In every frozen self-supervised or ImageNet
encoder tested, the same construction is **6x to 39x worse at the median than doing nothing
with the click at all** — that is one click, each encoder's best variant, against
`blob_native`. Allowing three clicks moves the band to 2.8-38.1x; across every frozen arm,
variant, pass and seed the tails run 1.5x to 65x. The
representation, not the interaction, was the bottleneck the whole time. The caveat is
load-bearing: the only encoder that works here was trained on MIDOG++, which contains both
of these ROIs, so its margin is an upper bound and not an estimate.

## read_50: candidates read to reach 50% of the mitotic figures

Min / median / **worst** over 5 seeds. `cl2n_snap` is each encoder's best variant.

**1 click**

| arm | 246.tiff | 301.tiff |
|---|---|---|
| `blob_native` — no click | 193 / 193 / **193** | 167 / 167 / **167** |
| `chromatin` — no click | 193 / 210 / **210** | 191 / 192 / **192** |
| **`fcos_midogpp_r18`** (leaky) | 60 / 64 / **80** | 114 / 119 / **220** |
| `lunit_dino_vits16_16um` | 1816 / 2047 / **6032** | 317 / 1067 / **3575** |
| `resnet18_in1k_16um` | 3284 / 4209 / **8863** | 962 / 2187 / **5399** |
| `kaiko_vits16` | 5662 / 7496 / **11701** | 3352 / 4618 / **7503** |
| `rawpixel_ccoeff` = `TM_CCOEFF_NORMED` | 1736 / 2404 / **4526** | 1788 / 4473 / **4866** |

**3 clicks** (SimpleShot prototype)

| arm | 246.tiff | 301.tiff |
|---|---|---|
| `blob_native` — no click | 192 / 192 / **192** | 166 / 166 / **166** |
| **`fcos_midogpp_r18`** (leaky) | 59 / 61 / **69** | 113 / 123 / **135** |
| `lunit_dino_vits16_16um` | 642 / 1156 / **1313** | 247 / 458 / **607** |
| `rawpixel_ccoeff` | 641 / 916 / **1149** | 1475 / 1825 / **2482** |

**The pre-registered gate — beat the best seedless arm on the worst of 5 seeds, on both
ROIs — fails at one click and passes at three.** At one click the task-trained arm wins
246.tiff outright (80 against 193) but its worst seed on 301.tiff is 220 against
`blob_native`'s 167. At three clicks it wins both on every seed (69 vs 192; 135 vs 166).
No other encoder comes within 3x of the bar in any configuration.

The premise test's primary metric agrees and more emphatically. recall@1000, median over
seeds, **all at one click**: **0.965 / 0.986** for the task-trained arm against
`blob_native`'s 0.861 / 0.931 and `chromatin`'s 0.774 / 0.760, while every frozen-encoder
arm lands between 0.061 and 0.479. At three clicks the same four quantities are 0.965 /
0.991, 0.867 / 0.930, 0.770 / 0.763, and 0.062-0.665.

## The control that decides what the click is worth

Same representation, same candidate set, same evaluation ground truth; **only the clicked
cell changes** (`results/click_rank_provenance.csv`). read_50, min-median-worst:

| encoder | click on a mitotic figure | on a rejected look-alike | on a random nucleus |
|---|---|---|---|
| `fcos_midogpp_r18`, 246.tiff | 60 / 64 / **80** | 1083 / 4086 / **10640** | 931 / 3227 / **19643** |
| `fcos_midogpp_r18`, 301.tiff | 114 / 119 / **220** | 430 / 1347 / **2284** | 3304 / 16865 / **19617** |
| `lunit_dino_vits16_16um`, 301.tiff | 317 / 1067 / **3575** | 1220 / 4700 / **5209** | 3844 / 7045 / **11898** |

Paired per seed, a click on a genuine mitotic figure is **18-172x** better than a click on a
look-alike on 246.tiff and **3.7-20x** better on 301.tiff, in the same space. The effect is
a property of that space and not of clicking as such: in `lunit_dino_vits16_16um` the same
contrast on 301.tiff is 0.3-14.9x, i.e. on one seed of five a look-alike click produces the
*better* list. This is the premise test's seed-provenance contrast — which measured
+0.05 recall through `TM_CCOEFF_NORMED` — re-run in a representation that can carry the
signal. **The click was never the weak part. The scoring function was.**

## AUC diagnostics, and why they were not the gate

Median over seeds; `tp_vs_lookalike` / `tp_vs_unannotated`:

| statistic | 246.tiff | 301.tiff |
|---|---|---|
| `chromatin` | 0.767 / 0.962 | 0.692 / 0.935 |
| `blob_native` | 0.742 / 0.972 | 0.686 / 0.986 |
| `fcos_midogpp_r18_cl2n_snap` | **0.953 / 0.999** | **0.869 / 0.996** |
| `lunit_dino_vits16_16um_cl2n_snap` | 0.624 / 0.784 | 0.665 / 0.881 |
| `kaiko_vits16_cl2n_snap` | 0.559 / 0.591 | 0.585 / 0.661 |
| `resnet18_in1k_cl2n_snap` | 0.501 / 0.551 | 0.563 / 0.614 |

One property of this table needs stating rather than burying. The AUC labels come from
bucketing the candidate set in `nucleus_blobs`'s own output order -- and that order *is*
`blob_native`'s ranking, since the detector returns rows sorted by its score. So the
positive-label set is fixed by one of the arms under comparison, and that matters: 57%
(246.tiff) and 49% (301.tiff) of mitotic annotations have more than one candidate inside the
match radius, so which candidate is called the true positive is genuinely arm-dependent for
about half the positives. The bias runs in favour of `blob_native` and **against every other arm in the
table, `chromatin` included** -- `chromatin` is seedless and is disadvantaged by the same
mechanism. The direction is exact rather than plausible: `greedy_match` walks detections in
row order and takes the nearest free ground truth, so the positive label goes to whichever
in-radius candidate `blob_native` ranks highest, and any other valid labelling can only
lower its AUC. The one loophole -- a candidate inside the radius of two annotations, where
"nearest free GT" could override "highest blob score" -- affects 0 of 21,059 candidates on
246.tiff and 2 of 20,370 on 301.tiff. The task-trained arm still posts 0.953/0.869 against
`blob_native`'s 0.742/0.686 under a labelling chosen to favour `blob_native`.

**None of this touches `read_50`.** `evaluate_arms` calls `bucket_detections` per arm on
that arm's own sorted frame; the canonical bucketing exists only to build the AUC label
masks. Do not "fix" this into a single shared bucketing -- it would silently score a
different ranking than the one reported.

Reading depth is set by the unannotated contrast — 20,000 unannotated candidates against
~110 look-alikes — which is exactly why the plan gated on read_50 and not on AUC. The
frozen encoders sit at 0.55-0.88 against unannotated where chromatin is 0.94-0.96, and that
gap is the whole story of their read_50 columns.

## What was actually wrong with the frozen representations

Measured on 301.tiff, seed 0's query, CL2N cosine over all 20,370 candidates. `concen` is
the norm of the mean unit embedding (1.0 would mean every crop points the same way);
`rho` columns are Spearman correlations of the cosine ranking against two *context* proxies
-- whole-crop brightness and the number of candidates within 150 px -- and against the
central cell's own chromatin density.

| encoder | concen | sd(cos) | sd(cl2n) | rho brightness | rho density | rho chromatin | AUC vs unann. |
|---|---:|---:|---:|---:|---:|---:|---:|
| `lunit_dino_vits16` 32 um | 0.904 | 0.076 | 0.200 | **0.33** | 0.20 | 0.09 | 0.688 |
| `kaiko_vits16` 32 um | 0.894 | 0.078 | 0.189 | **0.34** | 0.18 | 0.11 | 0.719 |
| `lunit_dino_vits16` 16 um | 0.860 | 0.075 | 0.176 | 0.22 | 0.08 | 0.20 | 0.876 |
| `resnet18_in1k` 16 um | 0.841 | 0.099 | 0.228 | 0.17 | 0.03 | 0.21 | 0.791 |
| `fcos_midogpp_r18` 32 um | 0.954 | 0.086 | **0.417** | **0.03** | **-0.04** | 0.24 | **0.997** |
| `chromatin` (for reference) | — | — | — | -0.20 | 0.04 | 1.00 | 0.935 |

Three things this says:

1. **Every embedding in an ROI points nearly the same way** (0.84-0.95). Raw cosine over
   20,000 crops of one tumour has a standard deviation of 0.08 around a mean near 0.85 --
   the discriminating signal is a small residual on top of a large "this is that tissue"
   component. That is why CL2N, which removes exactly that component, helps everywhere.
2. **Concentration is not the defect; what is left after centering is.** `fcos_midogpp_r18`
   is the *most* concentrated encoder here and has by far the largest centred spread
   (0.417 against 0.18-0.23). The frozen encoders do not lack dynamic range so much as
   spend it on the wrong thing.
3. **What they spend it on is context.** At 32 um the two pathology encoders' similarity
   ranks candidates mostly by whole-crop brightness (rho 0.33-0.34) and local cellularity
   (0.18-0.20), and barely by the central cell (0.09-0.11). Cutting to 16 um moves brightness
   0.33 -> 0.22 and chromatin 0.09 -> 0.20, which is the mechanism behind the field-of-view
   effect below. The task-trained backbone is the mirror image: **uncorrelated with context**
   (0.03 and -0.04) and tracking the cell, and it separates mitoses from unannotated nuclei
   at 0.997 -- above chromatin's own 0.935, so it is using more than chromatin density.

The honest scope of this. The readout under test is **frozen encoder -> global pooling ->
cosine**. A ViT-S/16 at 128 px pools an 8x8 token grid, in which a ~50 px nucleus is about
2x2 tokens: the cell is a small minority of what gets averaged, and the rest is the context
the correlations above expose. Taking the *central* tokens instead of the pooled vector, or
ROI-pooling at the cell's location, would discard that context directly rather than
crudely-by-cropping, and it is untested here. Given that cropping alone was worth 2-9x, it
is the cheapest remaining experiment before concluding anything about these encoders'
features in general. There is also a residual scale mismatch: these checkpoints were trained
near 0.5 um/px and 128 px input puts a patch-16 token at 4 um, still 2x finer than that.

## Four things that moved the numbers, in order of size

1. **What the representation is trained for.** 6-39x at one click, median, best variant
   per encoder. Nothing else comes close.
2. **Field of view.** Cutting the crop from 32 um to 16 um helps consistently, but the
   size depends on the pass and the first draft's "3-10x" was computed across mismatched
   variants. Like variant against like variant: at **one click**, median gains are 1.8-5.0x
   (`lunit_dino_vits16`) and 1.6-3.2x (`resnet18_in1k`), and worst-seed gains 1.05-2.41x
   and 0.95-2.06x -- on 246.tiff's worst seed, plain cosine at 16 um is very slightly
   *worse*. At **three clicks** the gains are systematically larger -- median gains span
   1.2-9.3x across the twelve (encoder x variant x ROI) cells, and every cell above 5x is a
   3-click cell while no 1-click cell exceeds 5.0x. Read 9.3x as the maximum of twelve
   cells over five correlated seeds, not as the effect size; the band is the honest
   summary. So: 2-3x at one click, 1.2-9.3x at three. A 32 um crop at this resolution is
   largely neighbouring tissue and the embedding is largely about it. `kaiko_vits16` was run only at 32 um, so it is
   reported at a configuration now known to be suboptimal -- but a 2-3x correction does not
   reach a ~30x deficit, so this is a fairness footnote, not a live hypothesis.
3. **CL2N centering** (SimpleShot). Consistent 10-40% improvement, free, click-independent.
4. **Number of clicks.** Three beats one everywhere, and — more importantly — it collapses
   the worst-seed spread: the task-trained arm goes from 114-220 to 113-135 on 301.tiff, and
   `lunit_dino_vits16_16um` from 317-3575 to 247-607. Read that stability claim with care:
   the five 3-click sets are rotations of the same five clicks, so any two share two of
   three members, and the tightness is partly prototype overlap rather than five independent
   draws. (Seed selection draws with replacement and can collide; here it did not -- all
   five `ann_id`s are distinct on both ROIs.)

## Per-seed instability is real and is a product problem

At one click on 301.tiff the task-trained arm's AUC against look-alikes ranges 0.519-0.910
across five seeds (0.831 / 0.519 / 0.893 / 0.869 / 0.910), and its uncentred variant reaches
**0.218 on seed 1 — an inverted ranking, far worse than chance.** The read_50 for that seed is 887 uncentred, 220 with CL2N and
snapping. One click in five produced a list twice as long as the other four. CL2N, snapping,
and a 3-click prototype each shrink this, and together they remove it on these two ROIs, but
a tool whose reading burden doubles depending on which mitosis the pathologist happened to
click needs this measured on more than two ROIs before it ships.

## What this does and does not establish

**Does not:** that any shippable model achieves these numbers. `fcos_midogpp_r18` is the
ResNet-18 body of a detector trained on MIDOG++, which contains 301.tiff and 246.tiff with
labels. Its features have seen these mitotic figures. **Every number in its rows is an
upper bound.** And the bound is weaker than "an inflated estimate": the entire pattern,
provenance control included, is reproducible by **pure memorisation**. If the backbone
memorised these ROIs' mitotic figures, a query taken from one of them lands in a tight
cluster containing the others while a look-alike query does not — which is exactly what was
observed, with no generalising mechanism. So not even the *direction* is established by this
experiment. The one piece of evidence against memorisation is external and indirect: the
same weights score F1 0.737-0.753 on the held-out MIDOG 2022 test set, so the representation
demonstrably generalises for detection. **That number is imported from the model's own
documentation and is the only figure in this log not re-derived from these artifacts.** And
it bounds less than it appears to: a detector can generalise to unseen domains *and* have
additionally memorised these 217 instances -- the two are not in tension, and held-out F1
constrains instance-specific overfitting on training data not at all. So it makes *pure*
memorisation unlikely while reducing the contamination of these particular numbers by
exactly zero.

**Does:** that the premise test's negative result cannot be blamed on the interaction
without also blaming `TM_CCOEFF_NORMED`. There exists a representation in which one click
reorders 20,000 candidates well enough to cut reading depth below every seedless arm — 2.4x
on 246.tiff at one click, and on both ROIs at three — and in which the ordering tracks
*what was clicked*: a look-alike click is 3.7-172x worse and a random nucleus 12.6-327x
worse than a click on a genuine mitotic figure. Whether any
*shippable* representation has that property is exactly what the next experiment must
establish.

**Also does:** rule out the cheap version of the product. Frozen pathology foundation-model
embeddings with cosine similarity — the architecture the literature review recommended as
the "ship it today" path — **do not work here at all**. `lunit_dino_vits16` and
`kaiko_vits16` are both TCGA-scale self-supervised pathology encoders, and at the median
both are 2.8-37x worse than a blob detector that ignores the click. §6 of the literature review is wrong on
this point and the recommendation has to change: the encoder must be adapted to the task,
not merely borrowed from the domain.

## What follows

1. **Re-run this with a mitosis-trained encoder that has never seen 301 or 246.** The
   MIDOG++ split is public and `jonas-amme/FCOS_Inference_CLI` is MIT-licensed with training
   code; a leave-these-two-ROIs-out retrain, or LoRA on a permissive backbone
   (Midnight-12k, H-optimus-0) against MIDOG++ minus these domains, converts every upper
   bound above into an estimate. This is now the single highest-value experiment.
2. **Tune the field of view** before anything else in the encoder. It was worth 2-3x at one
   click and up to 9x at three, and it was not in the plan.
3. **Default the product to a 3-click prototype**, not a single click. It is better on
   every ROI, every encoder and every metric measured, and it is what fixes the
   one-in-five bad click.
4. **Re-measure `kaiko_vits16` at 16 um** for fairness, not for hope: at its best variant
   it is 28-35x off the bar at the median (one click) and the field-of-view correction is
   worth 2-3x there. Note also that the
   best frozen arm anywhere -- `lunit_dino_vits16_16um_cl2n_snap`, 301.tiff, 3 clicks --
   still sits at 3.66x the seedless bar (607 against 166), so "no frozen encoder comes
   within 3x" is true but only just.
5. Extend beyond two ROIs before any claim about per-seed stability.
