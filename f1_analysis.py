"""Analysis of `f1_seed_sweep.py`, implementing `Research Logs/2026-09-04-f1-preregistration.md`.

Written and committed BEFORE the run's results were inspected. Every parameter below
(z=1.0, K=250, axis=chromatin_od, sign-flip permutation of Spearman rho, 10,000 draws,
two-sided, all 35 cells, seed-0-excluded as decisive) is fixed by that document.
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# --- pre-registered constants; changing any of these invalidates the pre-registration ---
Z_PRIMARY, K_PRIMARY, AXIS_PRIMARY = 1.0, 250, 'chromatin_od'
N_PERM, RNG_SEED = 10_000, 0
ARMS = ('base51', 'largest_cc')
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 40)


def load():
    cells = pd.read_csv('results/f1_seed_sweep_cells.csv')
    assert bool(cells['run_complete'].all()), "run did not complete -- refusing to analyse"
    assert len(cells) == int(cells['n_cells_expected'].iloc[0]), \
        f"{len(cells)} cells, expected {cells['n_cells_expected'].iloc[0]}"
    return pd.read_csv('results/f1_seed_sweep.csv'), cells


def paired(df, z, k, axis, value='recall_at_budget'):
    """One row per (domain, seed) with the paired delta. NaN is RIGHT-CENSORED, never dropped."""
    s = df[(df['z'] == z) & (df['budget'] == k) & (df['arm'] == axis)]
    w = s.pivot_table(index=['tumor_type', 'file_name', 'seed_index', 'tightened_size',
                            'area_ratio', 'null_cell'],
                      columns='arm_name', values=value, dropna=False).reset_index()
    w['delta'] = w['largest_cc'] - w['base51']
    return w


def signflip_spearman(x, d, n_perm=N_PERM, seed=RNG_SEED):
    """Exact randomization null for a paired design: each cell's delta sign is exchangeable.

    Each cell is its own stratum, so domain clustering is handled by construction rather
    than by a modelling assumption. Two-sided on |rho|.
    """
    x, d = np.asarray(x, float), np.asarray(d, float)
    rho_obs = spearmanr(x, d).statistic
    rng = np.random.default_rng(seed)
    flips = rng.choice([-1.0, 1.0], size=(n_perm, len(d)))
    null = np.array([spearmanr(x, d * f).statistic for f in flips])
    p = (np.sum(np.abs(null) >= abs(rho_obs) - 1e-12) + 1) / (n_perm + 1)
    return rho_obs, float(p), null


def report_primary(df, label, cells_filter=None):
    w = paired(df, Z_PRIMARY, K_PRIMARY, AXIS_PRIMARY)
    if cells_filter is not None:
        w = w[cells_filter(w)]
    assert w['delta'].notna().all(), "recall_at_budget should never be NaN"
    rho, p, _ = signflip_spearman(w['area_ratio'], w['delta'])
    print(f"\n{'='*94}\n{label}\n{'='*94}")
    print(w.sort_values('area_ratio')[['tumor_type', 'seed_index', 'tightened_size',
                                       'area_ratio', 'base51', 'largest_cc', 'delta']]
          .round(4).to_string(index=False))
    neg, pos, zero = (w.delta < 0).sum(), (w.delta > 0).sum(), (w.delta == 0).sum()
    print(f"\n  n={len(w)} cells ({zero} null)   signs: {neg} neg / {pos} pos / {zero} zero")
    print(f"  Spearman rho(area_ratio, delta) = {rho:+.4f}")
    print(f"  sign-flip permutation p (two-sided, {N_PERM} draws) = {p:.4f}")
    print(f"  mean delta {w.delta.mean():+.4f}   median {w.delta.median():+.4f}")
    return rho, p, w


def main():
    df, cells = load()
    print(f"loaded {len(df)} rows, {len(cells)} cells, run_complete={bool(cells.run_complete.all())}")
    print(f"gate: 0 retries anywhere = {not (cells.n_retries > 0).any()}; "
          f"null cells {int(cells.null_cell.sum())}/{len(cells)}")

    rho1, p1, w1 = report_primary(df, "PRIMARY (pre-registered): all 35 cells, "
                                      f"z={Z_PRIMARY}, K={K_PRIMARY}, axis={AXIS_PRIMARY}")
    rho2, p2, _ = report_primary(df, "DECISIVE (pre-registered): seed 0 EXCLUDED "
                                     "-- the cells that suggested the hypothesis",
                                 cells_filter=lambda w: w['seed_index'] != 0)

    print(f"\n{'='*94}\nSECONDARY 2: same test on the tm_score axis\n{'='*94}")
    wt = paired(df, Z_PRIMARY, K_PRIMARY, 'tm_score')
    rt, pt, _ = signflip_spearman(wt['area_ratio'], wt['delta'])
    wt0 = wt[wt.seed_index != 0]
    rt0, pt0, _ = signflip_spearman(wt0['area_ratio'], wt0['delta'])
    print(f"  all 35 cells : rho={rt:+.4f}, p={pt:.4f}")
    print(f"  seed-0 excl. : rho={rt0:+.4f}, p={pt0:.4f}")

    print(f"\n{'='*94}\nSECONDARY 3: WITHIN-arm seed SD under THIS configuration\n"
          f"(retires the audit's caveat that the only available SDs came from a different "
          f"pipeline)\n{'='*94}")
    s = df[(df.z == Z_PRIMARY) & (df.budget == K_PRIMARY) & (df.arm == AXIS_PRIMARY)]
    sd = (s.groupby(['tumor_type', 'arm_name'])['recall_at_budget']
          .agg(mean='mean', sd='std', min='min', max='max', n='size').reset_index())
    piv = sd.pivot(index='tumor_type', columns='arm_name', values=['mean', 'sd'])
    ad = w1.groupby('tumor_type')['delta'].agg(mean_delta='mean', sd_delta='std')
    print(piv.round(4).join(ad.round(4)).to_string())
    print(f"\n  pooled WITHIN-arm seed SD (unpaired): {sd['sd'].mean():.4f}")
    print(f"  pooled PAIRED delta SD              : {w1['delta'].std():.4f}")
    print(f"  |mean paired delta| / within-arm SD : {abs(w1.delta.mean())/sd['sd'].mean():.3f}")
    print("  These answer different questions: the first is 'is a single-seed number")
    print("  trustworthy', the second is 'is tightening's own effect resolvable'.")

    print(f"\n{'='*94}\nSECONDARY 5: read_95 right-censoring (NaN = deeper than the whole list)\n{'='*94}")
    r = paired(df, Z_PRIMARY, K_PRIMARY, AXIS_PRIMARY, value='read_95')
    print(f"  at z={Z_PRIMARY}: base51 NaN {int(r.base51.isna().sum())}, "
          f"largest_cc NaN {int(r.largest_cc.isna().sum())}")
    tot = cells[['read95_unreachable_base51', 'read95_unreachable_largest_cc']].sum()
    print(f"  across the whole z grid, unreachable (arm, z) cells: "
          f"base51 {int(tot.iloc[0])}, largest_cc {int(tot.iloc[1])}")
    print("  Never dropped: NaN ranks above every finite value (it means 'deeper than the list').")

    print(f"\n{'='*94}\nSECONDARY 6: paired volume / coverage deltas beside the quality delta\n{'='*94}")
    for v in ('n_detections', 'coverage_frac'):
        pv = paired(df, Z_PRIMARY, K_PRIMARY, AXIS_PRIMARY, value=v)
        g = pv.groupby('tumor_type')['delta'].agg(mean='mean', min='min', max='max')
        print(f"\n  {v}:"); print(g.round(4).to_string())

    print(f"\n{'='*94}\nSECONDARY 4 (EXPLORATORY, not a test): dose-response slope grid\n{'='*94}")
    rows = []
    for axis in df['arm'].unique():
        for z in sorted(df['z'].unique()):
            for k in sorted(df['budget'].unique()):
                w = paired(df, z, k, axis)
                if w['delta'].notna().sum() < 10 or w['delta'].nunique() < 3:
                    continue
                rows.append(dict(axis=axis, z=z, K=k,
                                 rho=round(spearmanr(w.area_ratio, w.delta).statistic, 3)))
    grid = pd.DataFrame(rows)
    for axis in grid['axis'].unique():
        print(f"\n  Spearman rho, axis={axis} (rows z, cols K) -- EXPLORATORY:")
        print(grid[grid.axis == axis].pivot(index='z', columns='K', values='rho').to_string())

    print(f"\n{'='*94}\nVERDICT against the pre-registered decision rule\n{'='*94}")
    if p1 < 0.05 and p2 < 0.05:
        v = ("H1 SUPPORTED: the effect is real and dose-dependent. lcc's per-domain signs "
             "were noise, but a monotone structure underlies them.")
    elif p1 < 0.05:
        v = ("NOT SUPPORTED: primary significant but the seed-0-excluded test is not. "
             "Per the pre-registration, the seed-0-excluded version is decisive -- treat the "
             "primary as generated by the cells that suggested the hypothesis.")
    else:
        v = ("NO RESOLVABLE DOSE-RESPONSE at n=35. Combined with the effect size against the "
             "within-arm seed SD, tightening is not a decision-relevant lever, and "
             "deliverable (d) should say so.")
    print(f"  primary p = {p1:.4f} | seed-0-excluded p = {p2:.4f}\n  -> {v}")


if __name__ == '__main__':
    main()
