"""Independent re-verification of the chromatin_density ranker's claims -- the evidence for D5.

Re-derives, from `.cache_tail/tp_fp_candidate_features.csv` and `results/f1_seed_sweep.csv`,
every number the repo uses to justify `chromatin_od` as the primary ranking axis. Nothing here
trusts a published CSV: the AUCs use a rank-based Mann-Whitney written here rather than
`tp_fp_separability.ipynb`'s scipy call, and the reading depths use an independent
cumulative-hit loop. It reproduces `results/tp_fp_auc_2class.csv`,
`results/tp_fp_auc_vs_fp_human.csv` and `results/tp_fp_reading_depth.csv` cell for cell, and
F4 section 2's `chromatin_od` 0.6403 vs `tm_score` 0.6081 exactly.

Six sections, in the order the argument runs:

1. **The founding premise, retested under D1.** `chromatin.py`'s docstring and commit `7c3af93`
   justify the statistic as the signal `TM_CCOEFF_NORMED` is blind to. D1 retired that matcher.
   `TM_CCOEFF` scales *linearly* with contrast (ratio 0.700 / 0.300 / 0.101 at contrast
   0.7 / 0.3 / 0.1) where `TM_CCOEFF_NORMED` does not move (0.999 / 0.992 / 0.930).
2. **Which contrast the rationale came from.** `mean_intensity`'s Bhattacharyya distance is
   0.93-3.70 against *ordinary nuclei* and 0.16-0.33 against *look-alikes* in the three domains
   with real sample size. The commit quoted the first column.
3. **The head-to-head at F4/F5's declared primary configuration** (recall@250, z = 1.0),
   clustered at the ROI -- F5 section 8's own unit. Delta = +0.032, 95% CI [-0.047, +0.112],
   p = 0.36, positive on 3 of 7 ROIs.
4. **AUC against all non-mitotic and against pathologist-marked look-alikes.** `od51` is
   0.68-0.78 on the contrast that decides clinical value, no better than the search score.
5. **Reading depth.** `od51` wins at read_50 and has the longest tail of any candidate: median
   read_95 4,488 against `od31` 1,336, `od_falloff` 1,759, `mask_od_mean` 1,748.
6. **Why 51 px is the wrong window.** Between candidate pairs <= 25 px apart, `od51` correlates
   0.72-0.84 -- against `od31` at 0.44-0.63 and the search score at 0.45-0.64. The window is
   substantially reading the neighbour.

Two limits that bound every number here, both inherited from the sources:

* **The 7-ROI densest-per-domain draw**, via `experiment.select_domain_images`, which documents
  itself as optimistic. D5 requires the 14-ROI `images/extra_valid` set for any re-measurement.
* **One seed** (`SEED_INDEX = 0`) for sections 4-6; only section 3 has five.

Run with an interpreter whose numpy matches the compiled cv2 (numpy 1.26 here, not 2.x):

    /Users/mohinianand/anaconda3/bin/python verify_chromatin_ranker.py
"""
import os
import sys

import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy import stats

REPO = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO)
from midog_utils import dataset as ds

TPB, LOOKB, UNANNB = 'human_correct_label', 'human_rejected_label', 'non_human_findings'


def auc(pos, neg):
    """Tie-aware AUC from average ranks. Deliberately not scipy.mannwhitneyu."""
    pos = np.asarray(pos, float); neg = np.asarray(neg, float)
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    n1, n2 = len(pos), len(neg)
    if n1 < 5 or n2 < 5:
        return np.nan
    r = pd.Series(np.concatenate([pos, neg])).rank(method='average').to_numpy()
    return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n2)


def depth(s, key, fracs=(0.5, 0.9, 0.95, 1.0)):
    """Candidates a reader opens to reach each fraction of the ROI's mitotic figures."""
    hit = (s.sort_values(key, ascending=False, kind='mergesort').bucket == TPB).to_numpy()
    cum = np.cumsum(hit)
    return {f'read_{int(f*100)}': int(np.searchsorted(cum, int(np.ceil(f * hit.sum())))) + 1
            for f in fracs}


def load_work():
    """The tp_fp candidate pool at z >= 1.0, with duplicate-on-mitosis FPs set aside."""
    cand = pd.read_csv(f'{REPO}/.cache_tail/tp_fp_candidate_features.csv')
    summ = pd.read_csv(f'{REPO}/results/tp_fp_extract_summary.csv').set_index('file_name')
    op = cand[cand.z >= 1.0].copy().reset_index(drop=True)
    _, anns = ds.load_annotations()
    op['dist_to_gt'] = np.nan
    for fn, s in op.groupby('file_name'):
        gt = ds.image_annotations(anns, fn)
        gt = gt[gt.ann_id != int(summ.loc[fn, 'seed_ann_id'])]
        d, _ = cKDTree(gt[['cx', 'cy']].to_numpy()).query(s[['cx', 'cy']].to_numpy(), k=1)
        op.loc[s.index, 'dist_to_gt'] = d
    on_ann = (op.bucket == UNANNB) & (op.dist_to_gt <= op.file_name.map(summ.match_radius_px))
    return op[~on_ann].copy()


def contrast_invariance():
    """CCOEFF vs CCOEFF_NORMED under I -> mean + a*(I - mean). Closed form, no cv2."""
    rng = np.random.default_rng(0)
    yy, xx = np.mgrid[-25:26, -25:26]
    T = 0.02 + 0.20 * np.exp(-(yy**2 + xx**2) / (2 * 9.0**2))
    out = []
    for a in (1.0, 0.7, 0.5, 0.3, 0.1):
        I = T.mean() + a * (T - T.mean()) + rng.normal(0, 0.002, T.shape)
        t, i = T - T.mean(), I - I.mean()
        out.append(dict(contrast=a, ccoeff=float((t * i).sum()),
                        ccoeff_normed=float((t * i).sum() /
                                            (np.linalg.norm(t) * np.linalg.norm(i)))))
    d = pd.DataFrame(out)
    d['ccoeff_ratio'] = (d.ccoeff / d.ccoeff.iloc[0]).round(4)
    d['normed_ratio'] = (d.ccoeff_normed / d.ccoeff_normed.iloc[0]).round(4)
    return d


def head_to_head(z=1.0, budget=250):
    """chromatin_od vs tm_score on F1's own sweep, clustered at the ROI (F5's own unit)."""
    d = pd.read_csv(f'{REPO}/results/f1_seed_sweep.csv')
    s = d[(d.z == z) & (d.budget == budget)]
    p = s.pivot_table(index=['file_name', 'seed_index', 'arm_name'],
                      columns='axis', values='recall_at_budget')
    dl = p.chromatin_od - p.tm_score
    per = dl.groupby('file_name').mean()
    se = per.std(ddof=1) / np.sqrt(len(per))
    tcrit = stats.t.ppf(.975, len(per) - 1)
    return dict(n_cells=len(p), chromatin=round(p.chromatin_od.mean(), 4),
                tm=round(p.tm_score.mean(), 4), delta=round(dl.mean(), 4),
                roi_ci=(round(per.mean() - tcrit * se, 4), round(per.mean() + tcrit * se, 4)),
                p_roi=round(stats.ttest_1samp(per, 0).pvalue, 4),
                rois_positive=f'{(per > 0).sum()}/{len(per)}')


def neighbour_sharing(work, sep=25.0):
    """Correlation of a feature between candidate pairs whose windows overlap."""
    rows = []
    for fn, s in work.groupby('file_name', sort=False):
        pairs = cKDTree(s[['cx', 'cy']].to_numpy()).query_pairs(sep, output_type='ndarray')
        r = {'file': fn, 'n_pairs': len(pairs)}
        for f in ('od31', 'od51', 'od81', 'od_ctx', 'mask_od_mean', 'z'):
            v = s[f].to_numpy()
            r[f] = round(float(np.corrcoef(v[pairs[:, 0]], v[pairs[:, 1]])[0, 1]), 3)
        rows.append(r)
    return pd.DataFrame(rows)


if __name__ == '__main__':
    pd.set_option('display.width', 300); pd.set_option('display.max_columns', 60)
    print('--- 1. is the founding premise still true under D1 (TM_CCOEFF)? ---')
    print(contrast_invariance().to_string(index=False))
    print('\n--- 2. bhattacharyya: which contrast was mean_intensity measured on? ---')
    b = pd.read_csv(f'{REPO}/results/morph_diag_bhattacharyya.csv')
    print(b[['tumor_type', 'n_mitotic', 'n_lookalike',
             'mean_intensity__vs_nucleus', 'mean_intensity__vs_lookalike']].to_string(index=False))
    print('\n--- 3. chromatin_od vs tm_score at F4/F5 declared primary config ---')
    print(head_to_head())
    work = load_work()
    F = ['z', 'od31', 'od51', 'od81', 'od_contrast', 'od_falloff', 'mask_od_mean']
    print('\n--- 4. AUC, mitotic vs all non-mitotic / vs pathologist look-alike ---')
    for name, sel in (('2-class', lambda s: s.bucket != TPB),
                      ('vs look-alike', lambda s: s.bucket == LOOKB)):
        rows = [dict(file=fn, **{f: round(auc(s[f][s.bucket == TPB], s[f][sel(s)]), 3) for f in F})
                for fn, s in work.groupby('file_name', sort=False)]
        print(f'\n{name}:'); print(pd.DataFrame(rows).to_string(index=False))
    print('\n--- 5. reading depth per ranker ---')
    rows = [dict(file=fn, ranker=k, **depth(s, k))
            for fn, s in work.groupby('file_name', sort=False) for k in F]
    t = pd.DataFrame(rows)
    print(t.pivot(index='file', columns='ranker', values='read_95').to_string())
    print('\nmedian read_95:', t.groupby('ranker').read_95.median().round(0).to_dict())
    print('\n--- 6. neighbour window sharing (why 51 px is the wrong window) ---')
    print(neighbour_sharing(work).to_string(index=False))
