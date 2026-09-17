# `find_and_suppress` vs. the reference implementation

`midog_utils/find_and_suppress.py` is the MIDOG++ port of the reference `find_and_suppress`
search. This is the record of every place it deliberately changes *behaviour* from
`bbox tuning code reference/bbox_tuning.py:757`, not just style. For what the module
actually does today, see its own module docstring, not this file.

## The seed's own footprint is blanked, not the whole reference-style mask

The reference blacks out every existing annotation before correlating
(`apply_mask_to_image(..., mode='black')`). On H&E that creates zero-variance windows --
where `TM_CCOEFF_NORMED` is undefined -- and hard edges that rotated templates correlate
against. Here, only the seed's own refined-template footprint (`base_size x base_size`,
`template_match.blank_seed_square`) is blanked, in a copy of the search channel, after the
template is cut from the original but before correlation runs -- `TM_CCOEFF` (D1, unlike
`_NORMED`) is well-defined against a flat window (it computes to ~0, not an extreme value:
the window's own local mean equals the constant, so the mean-subtracted term is exactly
zero), and because matching is a sliding window, the corrupted response tapers off with the
window's shrinking overlap rather than stepping sharply at a hard edge -- no boundary
artefact detection. The minimum spacing between two MIDOG++ annotations anywhere in the
dataset is 26.2 px (403.tiff; 245.tiff is 26.6), below the ~30 px match radius, so blanking
a region as large as a full match radius (the disc this design's predecessor used, `D10`'s
amendment) could delete a legitimate detection of a neighbouring ground-truth object; the
smaller, exact template footprint accepts a narrower version of the same cost instead (see
`SELF_HIT_MASKING_PLAN.md` sec 2's "Known cost, accepted", including the specific case,
301.tiff seed_index=1, where a real secondary match 26.6 px from the seed survives because
`base_size`'s half-width is smaller than that distance).

## Scores are kept and used to rank

The reference computes `match score` and then drops the column at
`bbox_tuning.py:798-801`, before an NMS that is ordered by `x top left` rather than by
score.

## Boundary refinement, the stability filter, and the aspect-ratio / PSNR filters are not run

Evaluation here is point-based, so box geometry moves no metric; and the stability check
passes 0/9 objects on H&E (see the bbox-tuning research log), so keeping it would delete
essentially every detection.

## Nothing is random

The reference picks its template with an unseeded `sample(frac=1)` and, when it has too
many matches, keeps `max_templates` of them with an unseeded `random.sample` --
discarding the best matches at random. Here the seed is chosen by the caller and
truncation is by score.
