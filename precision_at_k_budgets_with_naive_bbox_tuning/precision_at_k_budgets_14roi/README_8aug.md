# Precision at fixed candidate budgets -- 8-augmentation templates

`precision_at_k_budgets_14roi_8aug.ipynb` is `precision_at_k_budgets_14roi.ipynb` (see
`README.md`) with exactly one change: the template bank matched against each ROI is 8
augmentations per seed click -- 4 rotations (0/90/180/270 degrees) x 2 flips x 1 scale
(`FSConfig(n_angles=4, flips=(False, True))`) -- fused into one response map by per-pixel
max (`tm.fused_response`), instead of the original's single template
(`FSConfig(n_angles=1, flips=(False,))`). This is the `rot90_4angles_2flips` variant of
`midog_utils.experiment.AUGMENTATION_VARIANTS`.

Everything else is identical except one downstream consequence of that config change: same
14 ROIs (`images/extra_valid`), same seed per ROI (`seed_index = 0`, same RNG stream, so
the same click is drawn on each ROI), same decisions D1/D2/D3/D5/D7, same
deep-floor-then-truncate shortcut for K in {10, 20, 30, 50}.

**The one downstream change:** `SELF_HIT_RADIUS` (the fixed radius `suppress()` uses to drop
the seed's own detection after NMS) is 10.0px here, vs. 5.0px in the original. Under a single
template the seed's self-correlation is provably the response map's global maximum, so NMS
alone already clears everything within one match radius of it -- 5.0px was pure headroom on
top of an empty annulus. Under 8 fused templates that guarantee doesn't automatically hold: a
rotated/flipped template can score marginally higher a few pixels off the exact click, which
left the seed's own detection sitting 6.08px away on 403.tiff, just past the original 5.0px
window. 10.0px covers that with margin. It was *not* widened all the way to the match radius
(~30px): `find_and_suppress.py` documents that the closest two real MIDOG++ annotations
anywhere in the dataset are 26.2px apart (403.tiff) / 26.6px (245.tiff) -- below the match
radius -- so a match-radius self-hit filter risks deleting a legitimate neighbouring detection,
not just the seed's own echo. See the notebook's `## The pipeline` cell for the full reasoning.

Outputs land in `../results/` as `precision_at_k_14roi_8aug_*.csv`, distinct from the
original's `precision_at_k_14roi_*.csv` so neither run overwrites the other.

## Runtime

matchTemplate cost scales with template count, but it is not the only per-ROI cost --
image load, channel conversion, peak extraction and NMS are all independent of template
count and dominate the wall-clock budget here. Across three full executions of this notebook
while developing it, total runtime ranged 169-267s for 14 ROIs on this CPU-only hardware (vs.
the original single-template notebook's one measured run at 193s) -- close to parity, well
short of an 8x multiplier, but with enough run-to-run variance on this machine that a
single-point ratio isn't worth quoting as fact.

## Running it elsewhere

Same environment and path layout as `README.md` -- see that file for the
`requirements-midog-utils.txt` pin (numpy must stay below 2.0) and the list of paths to
copy if running outside this checkout.
