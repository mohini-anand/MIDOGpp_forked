# Why `midog_utils/invariants.py`'s checks exist

`midog_utils/invariants.py` is assert-style structural checks that catch bug classes a
regression test cannot: `Research Logs/2026-08-31-next-steps-plan.md`, Step 0 rule 11,
notes that a test re-running seed 0 and asserting previously committed numbers cannot
catch the 2026-08-25 bug class *by construction* -- those bugs were present when the
baseline was committed, so such a test would lock them in. Each of the four checks below
is tied to a specific defect that actually occurred in this repo. This is the history;
the code itself now only names the check.

## 1. `check_no_cap`

A list length exactly equal to a configured maximum is almost never a coincidence.
`baselines.nucleus_blobs` used to default to `max_detections=20000` and silently bound on
the two densest ROIs, truncating to the *darkest* competitors and so biasing every rank
statistic derived from it.

(`baselines.py`'s own docstring still records this: "it used to default to...".)

## 2. `check_tissue_mask_covers_gt`

Already documented -- see `DECISIONS.md` D2 for the full account (506.tiff annotation
24249, 13 px from the ROI edge, click pixel gray 232 against the 220 cut).

## 3. `check_distinct_seeds`

An RNG rebuilt inside a loop (`np.random.default_rng(0)` per iteration) returns the same
draw every time, so a "5-seed sweep" silently becomes one seed measured five times.

The check is deliberately **not** the literal "distinct `seed_index` gives distinct
`seed_ann_id`" -- that statement is not an invariant of this experiment and would fail on
correct runs: seed selection draws uniformly *with replacement* from a pool of accepted
annotations, so a collision is expected, not a defect. The pools run 3-162 candidates and
5 draws collide with probability `1 - prod(1 - i/p)`, which is 6% at p=162 (301.tiff) and
certain at p<5 (350.tiff has 3 mitotic annotations after the seed is excluded).
`results/od_seed_sweep.csv`, produced by already-committed code using the exact RNG
construction this experiment uses, shows 4 of 7 ROIs drawing only 4 distinct annotations
across 5 seed indices -- correct behaviour that the literal check would wrongly report as
a bug.

What is asserted instead is the property the literal check was a proxy for: the RNG
streams themselves must differ (each `rng_spec` replayed and its first 8 draws compared),
and not every draw may be the same annotation when the pool offers a choice.

## 4. `check_nms_radius`

Already documented -- see `DECISIONS.md` D7 for the full account (a fixed 25.0 px radius
let two detections 26 px apart both survive inside one ground-truth object's 29.6-33.1 px
match radius, inflating the FROC's false-positive axis with duplicates).
