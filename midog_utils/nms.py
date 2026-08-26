"""Distance-based non-maximum suppression for point detections.

MIDOG++ ground truth is a point click, so suppression by centre distance is the
natural operation -- there are no real box extents to compute IoU over.

This mirrors `non_max_suppression_by_distance` in `utils/nms_WSI.py`, which is already
correctly score-ordered. It is re-implemented standalone only because importing that
module drags in `utils.detection_helper` and therefore fastai.

Note this is *not* what the reference `nms_with_area` does. That function walks the
DataFrame index and keeps the last element, so the survivor of a cluster is whichever
box happens to sit furthest right -- and the match score it computes is dropped
before NMS ever sees it.
"""

from __future__ import annotations

import numpy as np
from sklearn.neighbors import KDTree


def nms_by_distance(centers: np.ndarray, scores: np.ndarray, radius: float = 25.0) -> np.ndarray:
    """Greedy NMS in descending score order.

    Returns the indices to keep, ordered best-first, so downstream code can treat the
    result as a ranked detection list.
    """
    centers = np.asarray(centers, dtype=np.float64)
    scores = np.asarray(scores, dtype=np.float64)
    if len(centers) == 0:
        return np.zeros(0, dtype=int)

    tree = KDTree(centers)
    neighbours = tree.query_radius(centers, r=radius)

    order = np.argsort(scores)[::-1]
    suppressed = np.zeros(len(centers), dtype=bool)
    keep = []
    for idx in order:
        if suppressed[idx]:
            continue
        keep.append(idx)
        suppressed[neighbours[idx]] = True
    return np.asarray(keep, dtype=int)
