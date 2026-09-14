"""Distance-based non-maximum suppression for point detections."""

from __future__ import annotations

import numpy as np
from sklearn.neighbors import KDTree


def nms_by_distance(centers: np.ndarray, scores: np.ndarray, radius: float = 25.0) -> np.ndarray:
    """
        Greedy NMS in descending score order.

        centers (np.ndarray): (N, 2) point coordinates.
        scores (np.ndarray): (N,) scores, higher is better.
        radius (float): suppression radius.

        Returns np.ndarray: indices to keep, ordered best-first.
    """
    centers = np.asarray(centers, dtype=np.float64)
    scores = np.asarray(scores, dtype=np.float64)
    if len(centers) == 0:
        return np.zeros(0, dtype=int)

    tree = KDTree(centers)
    neighbours = tree.query_radius(centers, r=radius)

    order = np.argsort(-scores, kind="stable")  # stable so ties keep the caller's incoming order
    suppressed = np.zeros(len(centers), dtype=bool)
    keep = []
    for idx in order:
        if suppressed[idx]:
            continue
        keep.append(idx)
        suppressed[neighbours[idx]] = True
    return np.asarray(keep, dtype=int)
