"""Assisted-annotation utilities written specifically for MIDOG++.

This package is a deliberate re-implementation of the `find_and_supress` pipeline in
`bbox tuning code reference/`, not a wrapper around it. That code was written for
bright, near-isotropic particles on a uniform dark background (hydrogel particle
tracking) and carries several bugs that would silently invalidate a histology
experiment -- most importantly a dead rotation search and a non-score-ranked NMS.
See `Research Logs/` for the full catalogue.

Nothing here imports fastai or SlideRunner, so it can be used without the pinned
training environment in `requirements.txt`.
"""

from . import (  # noqa: F401
    dataset, channels, template_match, nms, find_and_suppress,
    baselines, evaluate, viz, experiment, seed_selection,
)
