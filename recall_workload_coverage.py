"""Precompute `coverage_frac` as a function of list length, for every (ROI, seed) cell.

Why this exists
---------------
`midog_utils/evaluate.py` opens by requiring `coverage_frac` before **any** un-budgeted
full-list number: *"a detection list long enough to tile the ROI answers 'is this annotation
within the match radius of some detection?' by geometry rather than by evidence."* Every
recall number in `tm_recall_workload_curve.ipynb`'s tolerance, count, z-sweep and LODO
sections is exactly that, and the first version of that notebook reported none of them beside
a coverage figure. This closes that gap.

Why it is not just a loop over `evaluate.coverage_fraction`
-----------------------------------------------------------
That function costs ~0.27 s per call (a KDTree query over ~127k probe points), and the notebook
needs it at 17,194 distinct (cell, cutoff) pairs -- Table 2's z sweep alone is 70 cells x 151
cutoffs = 10,570 of them. That is ~77 minutes, recomputed on every re-read. (An earlier version
of this docstring said 5.4 s and 840 pairs; both were wrong, in opposite directions, and the
conclusion survived only because the pair count was understated by more than the per-call cost
was overstated. The notebook prints its own coverage-call count at the end of every run, which
is where 17,194 comes from, so this figure cannot go stale unnoticed.)

The whole curve is obtainable from **one** query per cell instead. For a probe point `p`, let
`rank_min(p)` be the lowest rank among detections within the match radius of `p`. Because a
cutoff keeps exactly the pool's first `n` rows (the prefix identity), `p` is covered at list
length `n` **iff** `rank_min(p) <= n`. So

    coverage(n) = #{p : rank_min(p) <= n} / n_probes

for every `n` at once. `rank_min` comes from a single `query_radius` of the *detections*
against a KDTree of the *probes*, filled in reverse rank order so the lowest rank wins.

Stride, probe grid and radius are `evaluate.coverage_fraction`'s own, so the result is
bit-for-bit identical to calling it -- verified per cell against the real function at FOUR
points on each curve before anything is written: z = 1.0, 2.0, 3.0 and the full pool. A single
gate point would leave the shape of the curve between the ends unchecked, which is most of what
this file is for.

Output
------
``results/tm_recall_workload_coverage.npz`` -- per cell, ``<fn>|<seed>|cum``: the cumulative
count of covered probes at list length 0..n_pool (so ``cum[n] / n_probes`` is the coverage of
a list of length ``n``, exactly, with no grid), plus ``<fn>|<seed>|n_probes``.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from sklearn.neighbors import KDTree

from midog_utils import evaluate as ev

CELLS = 'results/tm_recall_workload_cells.csv'
POOLZ = 'results/tm_recall_workload_pool_z.npz'
OUT = 'results/tm_recall_workload_coverage.npz'
STRIDE = 16          # evaluate.coverage_fraction's default -- do not change without re-gating


def main():
    t0 = time.time()
    cells = pd.read_csv(CELLS, float_precision='round_trip').set_index(
        ['file_name', 'seed_index']).sort_index()
    pool = np.load(POOLZ)
    out = {}

    for a, b in cells.index:
        k = (str(a), int(b))
        t = time.time()
        z = pool[f'{k[0]}|{k[1]}|z']
        det = np.stack([pool[f'{k[0]}|{k[1]}|cx'], pool[f'{k[0]}|{k[1]}|cy']],
                       axis=1).astype(np.float64)
        h, w = int(cells.loc[k, 'roi_h']), int(cells.loc[k, 'roi_w'])
        radius = float(cells.loc[k, 'match_radius_px'])

        ys, xs = np.mgrid[0:h:STRIDE, 0:w:STRIDE]
        probes = np.stack([xs.ravel(), ys.ravel()], axis=1).astype(np.float64)
        neigh = KDTree(probes).query_radius(det, radius)     # det is already rank-ordered

        SENT = len(det) + 1
        rank_min = np.full(len(probes), SENT, dtype=np.int64)
        for i in range(len(neigh) - 1, -1, -1):   # reverse rank order: the lowest rank wins
            rank_min[neigh[i]] = i + 1

        # cum[n] = probes covered by a list of length n. cum[0] = 0 by construction.
        cum = np.zeros(len(det) + 1, dtype=np.int64)
        covered = rank_min[rank_min <= len(det)]
        cum[1:] = np.bincount(covered, minlength=len(det) + 1)[1:].cumsum()

        # Gate: identical to the function this replaces, at four points ALONG the curve --
        # three operating points plus the full pool, not just the current one.
        n1 = int(np.searchsorted(-z, -1.0, side='right'))
        for n in (n1, int(np.searchsorted(-z, -2.0, side='right')),
                  int(np.searchsorted(-z, -3.0, side='right')), len(det)):
            ref = ev.coverage_fraction(det[:n], (h, w), radius) if n else 0.0
            mine = cum[n] / len(probes)
            assert abs(ref - mine) < 1e-12, \
                f'{k}: curve {mine} != coverage_fraction {ref} at n={n}'
        assert cum[0] == 0 and bool(np.all(np.diff(cum) >= 0)), f'{k}: curve is not a cdf'
        mine = cum[n1] / len(probes)

        out[f'{k[0]}|{k[1]}|cum'] = cum.astype(np.int32)
        out[f'{k[0]}|{k[1]}|n_probes'] = np.int64(len(probes))
        print(f'[{k[0]} s{k[1]}] n_pool={len(det):6d} probes={len(probes):7d} '
              f'cov(z=1.0)={mine:.4f} [{time.time() - t:.1f}s]', flush=True)

    np.savez_compressed(OUT, **out)
    print(f'\ntotal {time.time() - t0:.0f}s | wrote {OUT} for {len(cells)} cells')
    print('Every cell gated against evaluate.coverage_fraction at z = 1.0, 2.0, 3.0 '
          'and the full pool.')


if __name__ == '__main__':
    main()
