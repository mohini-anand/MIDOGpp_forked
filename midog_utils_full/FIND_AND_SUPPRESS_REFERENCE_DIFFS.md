# `find_and_suppress` vs. the reference implementation

`midog_utils/find_and_suppress.py` is the MIDOG++ port of the reference `find_and_suppress`
search. This is the record of every place it deliberately changes *behaviour* from
`bbox tuning code reference/bbox_tuning.py:757`, not just style. For what the module
actually does today, see its own module docstring, not this file.

## No image masking

The reference blacks out every existing annotation before correlating
(`apply_mask_to_image(..., mode='black')`). On H&E that creates zero-variance windows --
where `TM_CCOEFF_NORMED` is undefined -- and hard edges that rotated templates correlate
against. The seed's own detection is dropped afterwards instead, which has the same effect
with no artefact. (Note the reference's `'none'` mode does not disable masking; it sets
masked pixels to 1.)

## The self-hit is dropped with a tight radius, not the match radius

The minimum spacing between two MIDOG++ annotations anywhere in the dataset is 26.2 px
(403.tiff; 245.tiff is 26.6), below the ~30 px match radius, so dropping everything within
a match radius of the seed could delete a legitimate detection of a neighbouring
ground-truth object.

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
