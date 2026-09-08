"""Analysis of `f2_base_size_dose.py`, implementing `Research Logs/2026-09-08-f2-preregistration.md`.

Written and committed BEFORE the run's results exist. Every parameter below is fixed by that
document. F1's Robustness table was computed ad hoc outside `f1_analysis.py` and used a test
the pre-registration did not specify (audit finding F1.2); this file exists so that cannot
happen again -- every number in the results log must come from here.

The primary is deliberately plain. F1 had to correlate its delta against an *observed* dose,
which forced sign-flip permutation over cells, an argument about null-cell anchoring, and a
duplicate-draw correction. With the dose assigned the contrast is within-cell, so the
estimator is a mean per dose with a cluster interval and nothing else.
"""
from __future__ import annotations
import itertools

import numpy as np
import pandas as pd

import f2_base_size_dose as run

# --- pre-registered constants (sections 3d, 6.1, 8); GATE 10 asserts they match the runner --
PRIMARY_AXIS, Z_PRIMARY, K_PRIMARY = 'tm_score', 1.0, 250
PRIMARY_RULE, REFERENCE = 'matched_z', 'b51'
DOSES = ('b19', 'b27', 'b33', 'b37', 'b43')          # the non-reference levels; |family| = 5
N_BOOT, BOOT_SEED, ALPHA = 10_000, 20260908, 0.05
F1_ROIS = ('094.tiff', '201.tiff', '246.tiff', '301.tiff', '402.tiff', '459.tiff', '548.tiff')
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 50)

assert (PRIMARY_AXIS, Z_PRIMARY, K_PRIMARY) == (run.PRIMARY_AXIS, run.Z_PRIMARY, run.K_PRIMARY)
assert len(DOSES) == run.HOLM_FAMILY, "the Holm family must be the five non-reference doses"
assert set(DOSES) | {REFERENCE} == {f'b{d}' for d in run.DOSES}


def load(suffix=''):
    cells = pd.read_csv(run.OUT_CELLS + suffix)
    assert bool(cells['run_complete'].all()), "run did not complete -- refusing to analyse"
    assert len(cells) == int(cells['n_cells_expected'].iloc[0]), \
        f"{len(cells)} cells, expected {cells['n_cells_expected'].iloc[0]}"
    df = pd.read_csv(run.OUT_MAIN + suffix)
    verif = pd.read_csv(run.OUT_VERIF + suffix)
    assert bool(verif['passed'].fillna(True).all()), "a gate failed -- read the verification CSV"
    return df, cells, verif


def cell_deltas(df, axis=PRIMARY_AXIS, rule=PRIMARY_RULE, z=Z_PRIMARY, k=K_PRIMARY,
                value='recall_at_budget'):
    """One row per (ROI, seed, dose) with the paired delta against that cell's own b51 arm.

    Pivots on (file_name, seed_index) x dose_tag -- the unique cell key. Putting anything
    else in the index makes pandas build a cartesian product padded with NaN, which inflates
    the row count and silently corrupts any n or mean computed from it (audit F1.7 / the
    `pivot_table(dropna=False)` bug F1 had to fix mid-flight).
    """
    s = df[(df['arm'].str.startswith(axis + '@')) & (df['rule'] == rule)
           & (df['budget'] == k)]
    s = s[s['z'].isna()] if rule == 'matched_n' else s[s['z'] == z]
    w = s.pivot(index=['file_name', 'seed_index'], columns='dose_tag', values=value)
    assert len(w) == s['file_name'].nunique() * s['seed_index'].nunique(), \
        f"cell_deltas: {len(w)} rows -- the pivot index is not the unique cell key"
    out = w[list(DOSES)].sub(w[REFERENCE], axis=0).stack().rename('delta').reset_index()
    meta = df[['file_name', 'seed_index', 'tumor_type']].drop_duplicates(
        ['file_name', 'seed_index'])
    return out.merge(meta, on=['file_name', 'seed_index'], validate='many_to_one')


def roi_means(d):
    """Mean-of-seeds within each ROI. The ROI is the design cluster, not the cell (section 8)."""
    return (d.groupby(['dose_tag', 'file_name', 'tumor_type'])['delta']
            .agg(mean='mean', n_seeds='size').reset_index())


def exact_signflip(x):
    """Exact randomization p over all 2^n sign assignments of the cluster means, two-sided.

    Exact rather than Monte-Carlo, so there is no permutation-count arbitrariness and no
    p floor. F1's results log published 0.0001 and 0.0012 for the same test because one was
    scipy's asymptotic p and the other a 10,000-draw permutation whose floor is 1/10001.
    """
    x = np.asarray(x, float)
    n = len(x)
    assert n <= 20, f"exact enumeration of 2^{n} is not affordable"
    signs = np.array(list(itertools.product([-1.0, 1.0], repeat=n)))
    null = signs @ x / n
    return float(np.mean(np.abs(null) >= abs(x.mean()) - 1e-15))


def cluster_boot(x, n_boot=N_BOOT, seed=BOOT_SEED, alpha=ALPHA):
    """Percentile CI from resampling the ROI means with replacement."""
    x = np.asarray(x, float)
    rng = np.random.default_rng(seed)
    draws = rng.choice(x, size=(n_boot, len(x)), replace=True).mean(axis=1)
    return float(np.quantile(draws, alpha / 2)), float(np.quantile(draws, 1 - alpha / 2))


def holm(pvals, alpha=ALPHA):
    """Holm-Bonferroni across the family, returned in the caller's order."""
    pvals = np.asarray(pvals, float)
    m = len(pvals)
    order = np.argsort(pvals)
    adj, running = np.empty(m), 0.0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * pvals[i])
        adj[i] = min(1.0, running)
    return adj, adj < alpha


def dose_table(d, label):
    """The pre-registered per-dose table: mean, CI, exact p, Holm-adjusted p, ROI signs."""
    rm = roi_means(d)
    rows = []
    for dose in DOSES:
        x = rm.loc[rm.dose_tag == dose, 'mean'].to_numpy()
        lo, hi = cluster_boot(x)
        rows.append(dict(dose=dose, base_size=int(dose[1:]),
                         area_ratio=round((int(dose[1:]) / 51) ** 2, 4), n_roi=len(x),
                         mean=x.mean(), ci_lo=lo, ci_hi=hi, p_exact=exact_signflip(x),
                         n_roi_neg=int((x < 0).sum()), n_roi_pos=int((x > 0).sum())))
    t = pd.DataFrame(rows)
    t['p_holm'], t['reject'] = holm(t['p_exact'].to_numpy())
    print(f"\n{'=' * 104}\n{label}\n{'=' * 104}")
    print(t.round(5).to_string(index=False))
    return t


def threshold_interval(t):
    """Section 6.2: an interval, never a point, and non-monotonicity reported as itself."""
    t = t.sort_values('base_size')
    rej = t.loc[t.reject, 'base_size'].tolist()
    keep = t.loc[~t.reject, 'base_size'].tolist()
    print(f"\n  doses rejecting 0 (Holm, alpha={ALPHA}): {rej or 'none'}")
    print(f"  doses not rejecting                  : {keep or 'none'}")
    if not rej:
        print("  -> CLEAN NEGATIVE: template size is not a decision-relevant lever at n=14 "
              "clusters.\n     F1's threshold does not replicate on the 14-ROI set.")
    elif not keep:
        print(f"  -> every dose down to {max(rej)} px still costs recall; the threshold is at "
              "or above 43 px\n     and this grid does not bracket it.")
    elif max(rej) < min(keep):
        print(f"  -> THRESHOLD INTERVAL: ({max(rej)}, {min(keep)}] px. Contiguous, so P1's "
              "monotone reading holds.")
    else:
        print(f"  -> NON-MONOTONE: rejection is not contiguous in dose order. The threshold "
              "framing does\n     not survive; the per-dose table above is the result.")
    return t


def main():
    df, cells, verif = load()
    print(f"loaded {len(df)} rows, {len(cells)} cells, {len(verif)} gate records; "
          f"all gates passed = {bool(verif['passed'].fillna(True).all())}")
    print(f"ROIs {df.file_name.nunique()}, seeds {df.seed_index.nunique()}, "
          f"doses {sorted(df.dose_tag.unique())}")

    # ---- PRIMARY (section 8) ------------------------------------------------------------
    d = cell_deltas(df)
    t = dose_table(d, f"PRIMARY (pre-registered): mean-of-{df.seed_index.nunique()} paired "
                      f"D(recall@{K_PRIMARY}) per ROI, z={Z_PRIMARY}, axis={PRIMARY_AXIS}, "
                      f"rule={PRIMARY_RULE}, Holm across {len(DOSES)}")
    threshold_interval(t)

    # ---- section 9: the two decontaminated readings, pre-committed -----------------------
    new7 = [f for f in df.file_name.unique() if f not in F1_ROIS]
    dose_table(d[d.file_name.isin(new7)],
               f"ROBUSTNESS (section 9a): the {len(new7)} ROIs F1 never ran -- underpowered "
               "by construction, reported anyway")
    dose_table(d[d.seed_index != 0], "ROBUSTNESS (section 9b): seed 0 excluded -- the draw "
                                     "F1 shares")

    # ---- section 8: effect size against each ROI's own within-arm seed SD ----------------
    print(f"\n{'=' * 104}\nEFFECT SIZE: |mean delta| against that ROI's own within-arm seed SD"
          f"\n{'=' * 104}")
    s = df[(df.arm == f'{PRIMARY_AXIS}@{REFERENCE}') & (df.rule == PRIMARY_RULE)
           & (df.z == Z_PRIMARY) & (df.budget == K_PRIMARY)]
    sd = s.groupby('file_name')['recall_at_budget'].std().rename('ref_seed_sd')
    rm = roi_means(d).merge(sd, on='file_name')
    rm['ratio'] = rm['mean'].abs() / rm['ref_seed_sd']
    print(rm.pivot(index=['file_name', 'tumor_type'], columns='dose_tag',
                   values='mean').round(4).join(sd.round(4)).to_string())
    print(f"\n  cells where |mean delta| exceeds that ROI's own seed SD: "
          f"{int((rm.ratio > 1).sum())} of {len(rm)}")

    # ---- section 6.1: the D4 reporting frame, descriptive, no p-value --------------------
    print(f"\n{'=' * 104}\nREPORTING FRAME (D4): recall@{K_PRIMARY} per domain, WORST ROI and "
          f"WORST seed, per dose -- descriptive\n{'=' * 104}")
    r = df[(df['arm'].str.startswith(PRIMARY_AXIS + '@')) & (df.rule == PRIMARY_RULE)
           & (df.z == Z_PRIMARY) & (df.budget == K_PRIMARY)]
    print(r.groupby(['tumor_type', 'dose_tag'])['recall_at_budget'].min().unstack()
          .round(4).to_string())

    # ---- S1 / S2 (section 6.4): the audit-F1.5 discriminator -----------------------------
    print(f"\n{'=' * 104}\nS1/S2 (unadjusted): does chromatin_od's dose curve converge to "
          f"tm_score's once the list length is matched?\n{'=' * 104}")
    curves = {}
    for axis, rule in [('tm_score', 'matched_z'), ('tm_score', 'matched_n'),
                       ('chromatin_od', 'matched_z'), ('chromatin_od', 'matched_n')]:
        rmx = roi_means(cell_deltas(df, axis=axis, rule=rule))
        curves[f'{axis}/{rule}'] = rmx.groupby('dose_tag')['mean'].mean()
    c = pd.DataFrame(curves).reindex(list(DOSES))
    print(c.round(5).to_string())
    gap_z = (c['chromatin_od/matched_z'] - c['tm_score/matched_z']).abs().mean()
    gap_n = (c['chromatin_od/matched_n'] - c['tm_score/matched_n']).abs().mean()
    print(f"\n  mean |chromatin_od - tm_score| gap: matched_z {gap_z:.5f} -> matched_n {gap_n:.5f}")
    print("  P2 predicts the gap SHRINKS. If it does not, F1's chromatin-specific reading "
          "stands and\n  the audit's list-length explanation is wrong -- a real outcome, "
          "reported either way.")

    # ---- S3 (section 3d): the axis head-to-head, consequence fixed in advance ------------
    print(f"\n{'=' * 104}\nS3 (unadjusted): chromatin_od - tm_score on the {REFERENCE} arm\n"
          f"{'=' * 104}")
    h = df[(df.dose_tag == REFERENCE) & (df.rule == PRIMARY_RULE) & (df.z == Z_PRIMARY)
           & (df.budget == K_PRIMARY)].pivot(index=['file_name', 'seed_index'],
                                             columns='axis', values='recall_at_budget')
    hd = (h['chromatin_od'] - h['tm_score']).groupby(level='file_name').mean()
    lo, hi = cluster_boot(hd.to_numpy())
    print(f"  mean {hd.mean():+.4f}  CI [{lo:+.4f}, {hi:+.4f}]  exact p "
          f"{exact_signflip(hd.to_numpy()):.4f}  positive on {int((hd > 0).sum())}/{len(hd)} ROIs")
    print("  Section 3d: if this replicates F5 section 2.1, both axes are reported with "
          "neither secondary.\n  The arbiter does NOT move and the Holm family stays at 5.")

    # ---- P3, P4, D2 (sections 4, 6.3, 6.5) ----------------------------------------------
    print(f"\n{'=' * 104}\nP3 volume / P4 censoring / D2 pool size beside every depth figure\n"
          f"{'=' * 104}")
    for v in ('n_detections', 'coverage_frac', 'read_95'):
        dv = cell_deltas(df, value=v)
        print(f"\n  {v} (delta vs {REFERENCE}), mean over ROI means:")
        print(roi_means(dv).groupby('dose_tag')['mean'].mean().reindex(list(DOSES))
              .round(4).to_string())
    nd = cell_deltas(df, value='n_detections')
    mono = (nd.pivot(index=['file_name', 'seed_index'], columns='dose_tag', values='delta')
            [list(DOSES)].diff(axis=1).iloc[:, 1:] <= 0).all(axis=1)
    print(f"\n  P3: n_detections rises monotonically as dose falls in {int(mono.sum())} of "
          f"{len(mono)} cells")
    cen = df[(df.rule == PRIMARY_RULE) & (df.z == Z_PRIMARY)].drop_duplicates(
        ['file_name', 'seed_index', 'arm'])
    print(f"  P4: read_95 censored (NaN) at z={Z_PRIMARY} in "
          f"{int(cen['read_95'].isna().sum())} of {len(cen)} (cell, arm) rows")

    # ---- section 6.6: the lcc residual, a diagnostic and not a test ----------------------
    print(f"\n{'=' * 104}\nSECTION 6.6 (diagnostic, not a test): lcc against the assigned "
          f"curve interpolated at its own size\n{'=' * 104}")
    wide = df[(df['arm'] == f'{PRIMARY_AXIS}@lcc') & (df.rule == PRIMARY_RULE)
              & (df.z == Z_PRIMARY) & (df.budget == K_PRIMARY)]
    ref = df[(df['arm'] == f'{PRIMARY_AXIS}@{REFERENCE}') & (df.rule == PRIMARY_RULE)
             & (df.z == Z_PRIMARY) & (df.budget == K_PRIMARY)]
    m = wide.merge(ref[['file_name', 'seed_index', 'recall_at_budget']],
                   on=['file_name', 'seed_index'], suffixes=('_lcc', '_ref'))
    m['delta_obs'] = m['recall_at_budget_lcc'] - m['recall_at_budget_ref']
    curve = d.groupby('dose_tag')['delta'].mean().reindex(list(DOSES))
    xs = np.array([int(t[1:]) for t in DOSES] + [51])
    ys = np.append(curve.to_numpy(), 0.0)
    m['delta_pred'] = np.interp(m['lcc_size'], xs, ys)
    m['residual'] = m['delta_obs'] - m['delta_pred']
    print(f"  lcc sizes {sorted(m.lcc_size.unique())}")
    print(f"  mean observed {m.delta_obs.mean():+.4f}   mean predicted by the curve "
          f"{m.delta_pred.mean():+.4f}   mean residual {m.residual.mean():+.4f}")
    print("  Centred on 0 means the rule does nothing beyond setting a size, which is what "
          "the notebook\n  claims. A systematic offset would contradict it.")


if __name__ == '__main__':
    main()
