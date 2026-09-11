"""Audit of `tp_fp_separability.ipynb` -- eight checks, each re-derived from the cache.

Reads `.cache_tail/tp_fp_candidate_features.csv` (the notebook's own input) and re-runs
the notebook's load-bearing steps against a definition or a baseline the notebook did not
use. Nothing here re-runs the search; the candidate pool is taken as given, so a
divergence is an analysis defect rather than a reproduction failure.

The three checks that changed a conclusion:

* **C1** `on_annotation` is built from `dataset.image_annotations(anns, fn)`, which pools
  **both** categories. 629 of the 1,211 rows the notebook sets aside as "false positives
  sitting on a mitotic figure" sit on a category-2 look-alike instead. Shared with
  `f1_seed_sweep.py`'s `n_dup_fp`, which is why the notebook's subset assertion passes.
* **C4** the notebook re-ranks by `od51` / the supervised score but keeps the bucket
  labels `evaluate.greedy_match` produced in **TM-score order**. Re-running the match in
  each ranker's own order is what `results/tp_fp_reading_depth.csv` should have stored.
* **C8/C9** the two candidate FP filters, measured against the baseline the notebook
  never ran: a one-line per-ROI percentile threshold on the single best feature.

Writes results/tp_fp_separability_audit_*.csv. Run with the notebook's kernel
(`~/anaconda3/bin/python3`), not the system 3.11 -- its cv2 is built against numpy 1.x.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.cluster import KMeans
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.neighbors import KDTree

from midog_utils import dataset as ds
from midog_utils import evaluate as ev

warnings.filterwarnings('ignore')
pd.set_option('display.width', 260)
pd.set_option('display.max_columns', 60)

Z_OP = 1.0
N_PROTO_DRAWS = 30
FEATURES = ['z', 'od51', 'od_contrast', 'od_falloff',
            'solidity', 'extent', 'tightened_size', 'area_frac_of_window']
CHROM = ['od51', 'od_contrast', 'od_falloff']
ALL = ['z', 'od31', 'od51', 'od81', 'od_ctx', 'od_contrast', 'od_falloff', 'area', 'solidity',
       'extent', 'eccentricity', 'perimeter', 'major_axis_length', 'minor_axis_length',
       'circularity', 'tightened_size', 'mask_od_mean', 'mask_od_max', 'area_frac_of_window']
TPB = ev.HUMAN_CORRECT_LABEL
UNANNB = ev.NON_HUMAN_FINDINGS


def load():
    cand = pd.read_csv('.cache_tail/tp_fp_candidate_features.csv')
    summ = pd.read_csv('results/tp_fp_extract_summary.csv').set_index('file_name')
    op = cand[cand.z >= Z_OP].copy().reset_index(drop=True)
    _, anns = ds.load_annotations()

    gt = {}
    for fn in summ.index:
        g = ds.image_annotations(anns, fn)
        gt[fn] = g[g.ann_id != int(summ.seed_ann_id[fn])].reset_index(drop=True)

    # the notebook's `dist_to_gt`, split by what the nearest annotation actually is
    for col in ('dist_all', 'dist_mit'):
        op[col] = np.nan
    for fn, s in op.groupby('file_name'):
        g, xy = gt[fn], s[['cx', 'cy']].to_numpy()
        op.loc[s.index, 'dist_all'] = KDTree(g[['cx', 'cy']].to_numpy()).query(xy, k=1)[0][:, 0]
        m = g[g.category_id == ds.MITOTIC]
        op.loc[s.index, 'dist_mit'] = KDTree(m[['cx', 'cy']].to_numpy()).query(xy, k=1)[0][:, 0]
    op['match_radius'] = op.file_name.map(summ.match_radius_px)
    op['on_annotation'] = (op.bucket == UNANNB) & (op.dist_all <= op.match_radius)   # notebook's
    op['on_mitosis'] = (op.bucket == UNANNB) & (op.dist_mit <= op.match_radius)      # as claimed
    op['on_look_only'] = op.on_annotation & ~op.on_mitosis
    return op, summ, gt


def cum_hits(xy, g, radius):
    """Cumulative distinct mitoses claimed, with the greedy match re-run in THIS order."""
    det_to_gt, _ = ev.greedy_match(xy, g[['cx', 'cy']].to_numpy(), radius)
    cls = g.category_id.to_numpy()
    return np.cumsum([(i >= 0 and cls[i] == ds.MITOTIC) for i in det_to_gt])


def read_at(cum, n_mit, q):
    ix = np.searchsorted(cum, int(np.ceil(q * n_mit)), side='left')
    return int(ix) + 1 if ix < len(cum) else np.nan


def rank_order(s, score):
    return s.iloc[np.argsort(-np.asarray(score), kind='mergesort')]


def _auc(a, b):
    a = np.asarray(a, float)[~np.isnan(np.asarray(a, float))]
    b = np.asarray(b, float)[~np.isnan(np.asarray(b, float))]
    if len(a) < 5 or len(b) < 5:
        return np.nan
    # scipy >= 1.7 returns U for the FIRST sample, so this is signed: > 0.5 = a scores higher
    return mannwhitneyu(a, b, alternative='two-sided').statistic / (len(a) * len(b))


def c1_duplicate_composition(op):
    print('\n### C1  what the 1,211 "FPs sitting on a mitotic figure" actually sit on')
    t = op.groupby('file_name').agg(mitotic=('bucket', lambda b: int((b == TPB).sum())),
                                    on_annotation=('on_annotation', 'sum'),
                                    on_mitosis=('on_mitosis', 'sum'),
                                    on_lookalike_only=('on_look_only', 'sum'))
    t['pct_lookalike'] = (100 * t.on_lookalike_only / t.on_annotation).round(1)
    t['dup_per_mitosis_claimed'] = (t.on_annotation / t.mitotic).round(2)
    t['dup_per_mitosis_true'] = (t.on_mitosis / t.mitotic).round(2)
    print(t.to_string())
    print(f'pooled {int(op.on_annotation.sum())} set aside, {int(op.on_mitosis.sum())} on a mitosis, '
          f'{int(op.on_look_only.sum())} on a look-alike '
          f'({100 * op.on_look_only.sum() / op.on_annotation.sum():.1f}%)')

    # sec 6's control, re-run on the two halves separately: the "20-30 px off-centre, so the
    # chromatin is weaker" mechanism is stated for a set that is half something else.
    print('\n    the notebook\'s sec-6 control, split -- AUC(mitotic vs each half):')
    rows = []
    for fn, s in op.groupby('file_name', sort=False):
        pos = s[s.bucket == TPB]
        r = {'file_name': fn}
        for half, sel in [('on_mitosis', s.on_mitosis), ('on_lookalike_only', s.on_look_only)]:
            neg = s[sel]
            r[f'n_{half}'] = len(neg)
            for f in ('z', 'od51', 'od_falloff'):
                r[f'{f}_vs_{half}'] = round(_auc(pos[f], neg[f]), 3)
        rows.append(r)
    ctrl = pd.DataFrame(rows)
    print(ctrl.to_string(index=False))
    return t.reset_index().merge(ctrl, on='file_name')


def c2_coverage(op, summ):
    print('\n### C2  coverage_fraction of the z >= 1.0 list (evaluate.py says read this first)')
    rows = []
    for fn, s in op.groupby('file_name', sort=False):
        r = summ.loc[fn]
        rows.append({'file_name': fn, 'n_cand': len(s),
                     'coverage_frac': round(ev.coverage_fraction(
                         s[['cx', 'cy']].to_numpy(), (int(r.roi_h), int(r.roi_w)),
                         float(r.match_radius_px), stride=16), 3)})
    t = pd.DataFrame(rows)
    print(t.to_string(index=False))
    return t


def c3_positive_centredness(op, summ):
    print('\n### C3  are the POSITIVES centred on their annotation? (control for C1\'s mechanism)')
    rows = []
    for fn, s in op.groupby('file_name', sort=False):
        d = s.loc[s.bucket == TPB, 'dist_mit']
        rows.append({'file_name': fn, 'n_mit': len(d), 'match_r': round(float(summ.match_radius_px[fn]), 1),
                     'median_px': round(d.median(), 1), 'q90_px': round(d.quantile(.9), 1),
                     'frac_over_20px': round(float((d > 20).mean()), 3)})
    t = pd.DataFrame(rows)
    print(t.to_string(index=False))
    return t


def c4_redundancy(work):
    print('\n### C4  correlated pairs still inside the kept 8 (the notebook\'s own line is 0.85)')
    sub = work[FEATURES].corr(method='spearman').abs()
    np.fill_diagonal(sub.values, 0)
    s = sub.stack().sort_values(ascending=False).drop_duplicates()
    t = s[s > 0.7].rename('abs_spearman').reset_index()
    t.columns = ['feature_a', 'feature_b', 'abs_spearman']
    print(t.to_string(index=False))
    return t


def c5_cluster_vs_threshold(work):
    print('\n### C5  k=2 clustering vs a per-ROI percentile threshold at MATCHED volume')
    rows = []
    for fn, s in work.groupby('file_name', sort=False):
        R = s[FEATURES].rank(pct=True)
        yv = (s.bucket == TPB).to_numpy().astype(int)
        lab = KMeans(n_clusters=2, n_init=10, random_state=0).fit_predict(R.to_numpy(float))
        n, t_ = np.bincount(lab), np.bincount(lab, weights=yv)
        b = int(np.argmax(t_ / n))
        keep = (lab == b).mean()
        r = {'file_name': fn, 'n_mit': int(yv.sum()), 'kmeans_keeps': round(keep, 3),
             'kmeans_recall': round(t_[b] / yv.sum(), 4)}
        for f in ['od_falloff', 'od_contrast', 'od51', 'z']:
            v = R[f].to_numpy()
            r[f'{f}_at_same_volume'] = round(float(yv[v >= np.quantile(v, 1 - keep)].sum() / yv.sum()), 4)
        rows.append(r)
    t = pd.DataFrame(rows)
    print(t.to_string(index=False))
    return t


def c6_rematched_depth(op, summ, gt):
    print('\n### C6  reading depth with the greedy match RE-RUN in each ranker\'s own order')
    rows = []
    for fn, s0 in op.groupby('file_name', sort=False):
        s = s0.reset_index(drop=True)
        g, radius = gt[fn], float(summ.match_radius_px[fn])
        n_mit = int((g.category_id == ds.MITOTIC).sum())
        R = s[FEATURES].rank(pct=True)
        axes = [('z', s.z.to_numpy()), ('od51', s.od51.to_numpy()),
                ('od_falloff', s.od_falloff.to_numpy()), ('od_contrast', s.od_contrast.to_numpy()),
                ('z+od_falloff_pct', (R.z + R.od_falloff).to_numpy() / 2)]
        for name, sc in axes:
            cum = cum_hits(rank_order(s, sc)[['cx', 'cy']].to_numpy(), g, radius)
            rows.append({'file_name': fn, 'n_mit': n_mit, 'n_cand': len(s), 'ranker': name,
                         'read_90': read_at(cum, n_mit, .90), 'read_95': read_at(cum, n_mit, .95),
                         'read_100': read_at(cum, n_mit, 1.0),
                         'rec_at_1000': round(cum[min(1000, len(cum)) - 1] / n_mit, 3),
                         'rec_at_2000': round(cum[min(2000, len(cum)) - 1] / n_mit, 3)})
    t = pd.DataFrame(rows)
    print(t.set_index(['file_name', 'ranker']).to_string())
    print('\nread_95 over the 7 ROIs:')
    print(t.groupby('ranker').read_95.agg(['min', 'median', 'max']).to_string())
    return t


def c7_supervised_vs_shipped(op, summ, gt):
    print('\n### C7  supervised LODO vs the SHIPPED ranker z (D5), match re-run per ranker')
    rk = op.groupby('file_name')[FEATURES].rank(pct=True)
    X, y = rk.to_numpy(float), (op.bucket == TPB).to_numpy().astype(int)
    roi = op.file_name.to_numpy()
    rows = []
    for fn in op.file_name.unique():
        te = roi == fn
        s = op[te].reset_index(drop=True)
        g, radius = gt[fn], float(summ.match_radius_px[fn])
        n_mit = int((g.category_id == ds.MITOTIC).sum())
        gb = HistGradientBoostingClassifier(max_iter=250, learning_rate=0.08, max_leaf_nodes=15,
                                            l2_regularization=1.0, class_weight='balanced',
                                            random_state=0).fit(X[~te], y[~te])
        c_gb = cum_hits(rank_order(s, gb.predict_proba(X[te])[:, 1])[['cx', 'cy']].to_numpy(), g, radius)
        c_z = cum_hits(rank_order(s, s.z.to_numpy())[['cx', 'cy']].to_numpy(), g, radius)
        rows.append({'file_name': fn, 'n_mit': n_mit,
                     'gb_read_95': read_at(c_gb, n_mit, .95), 'z_read_95': read_at(c_z, n_mit, .95),
                     'ratio_gb_over_z': round(read_at(c_gb, n_mit, .95) / read_at(c_z, n_mit, .95), 2),
                     'gb_read_100': read_at(c_gb, n_mit, 1.0), 'z_read_100': read_at(c_z, n_mit, 1.0),
                     'gb_rec_at_1000': round(c_gb[999] / n_mit, 3),
                     'z_rec_at_1000': round(c_z[999] / n_mit, 3)})
    t = pd.DataFrame(rows)
    print(t.to_string(index=False))
    return t


def c8_prototype_similarity(op, summ, gt):
    """Rank by distance to the clicked template in feature space, over many clicks.

    The prototype is drawn from the pool's own mitotic candidates rather than being the
    real seed (which passes `agreement_pool` and the largest-CC gates), so this is a
    lower bound on one-click similarity. The spread over draws is the point, not the
    median: a rule whose worst click is unusable is not deployable on an unseen domain.
    """
    print(f'\n### C8  distance-to-the-clicked-template as the ranker, {N_PROTO_DRAWS} clicks per ROI')
    rows = []
    for fn, s0 in op.groupby('file_name', sort=False):
        s = s0.reset_index(drop=True)
        g_all, radius = gt[fn], float(summ.match_radius_px[fn])
        R = {'proto_8feat': s[FEATURES].rank(pct=True).to_numpy(float),
             'proto_3chrom': s[CHROM].rank(pct=True).to_numpy(float)}
        mit = np.flatnonzero((s.bucket == TPB).to_numpy())
        picks = np.random.default_rng(0).choice(mit, size=min(N_PROTO_DRAWS, len(mit)), replace=False)
        for p in picks:
            g = g_all[g_all.ann_id != int(s.matched_ann_id.iloc[p])].reset_index(drop=True)
            n_mit = int((g.category_id == ds.MITOTIC).sum())
            keep = np.ones(len(s), bool)
            keep[p] = False
            for name, M in R.items():
                d = np.linalg.norm(M - M[p], axis=1)[keep]
                cum = cum_hits(rank_order(s[keep], -d)[['cx', 'cy']].to_numpy(), g, radius)
                rows.append({'file_name': fn, 'ranker': name, 'pick_row': int(p), 'n_mit': n_mit,
                             'read_95': read_at(cum, n_mit, .95),
                             'rec_at_1000': cum[min(1000, len(cum)) - 1] / n_mit})
    t = pd.DataFrame(rows)
    for m in t.ranker.unique():
        g = t[t.ranker == m].groupby('file_name')
        print(f'\n-- {m} --')
        print(pd.DataFrame({'n_mit': g.n_mit.median(),
                            'read95_p10': g.read_95.quantile(.1).round(),
                            'read95_med': g.read_95.median().round(),
                            'read95_p90': g.read_95.quantile(.9).round(),
                            'read95_worst_click': g.read_95.max(),
                            'rec1000_med': g.rec_at_1000.median().round(2),
                            'rec1000_worst_click': g.rec_at_1000.min().round(2)}).to_string())
    return t


def main():
    op, summ, gt = load()
    work = op[~op.on_annotation].copy().reset_index(drop=True)
    print(f'{len(op):,} candidates at z >= {Z_OP}, {int((op.bucket == TPB).sum())} mitotic, '
          f'{len(work):,} in the notebook\'s `work`')
    out = {
        'duplicate_composition': c1_duplicate_composition(op),
        'coverage': c2_coverage(op, summ),
        'positive_centredness': c3_positive_centredness(op, summ),
        'kept_feature_redundancy': c4_redundancy(work),
        'cluster_vs_threshold': c5_cluster_vs_threshold(work),
        'rematched_depth': c6_rematched_depth(op, summ, gt),
        'supervised_vs_shipped': c7_supervised_vs_shipped(op, summ, gt),
        'prototype_similarity': c8_prototype_similarity(op, summ, gt),
    }
    for name, df in out.items():
        df.to_csv(f'results/tp_fp_separability_audit_{name}.csv', index=False)
    print(f'\nwrote {len(out)} tables to results/tp_fp_separability_audit_*.csv')


if __name__ == '__main__':
    main()
