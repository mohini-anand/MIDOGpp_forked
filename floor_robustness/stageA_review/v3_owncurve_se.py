"""Noise of an own-ROI-curve PL vs the pooled-curve PL: uniform landing c on all 49 ROIs, per-ROI losses from each ROI's own 6-click ladder."""
import numpy as np, pandas as pd
from scipy import stats
REPO = '/Users/mohinianand/Library/Mobile Documents/com~apple~CloudDocs/Desktop/AnnotateDx/MIDOGpp_forked'
pr = pd.read_csv(f'{REPO}/results/tm_floor_definition_49roi_6seed_per_run.csv')
tp = pr.pivot_table(index=['file_name', 'seed_index'], columns='arm', values='tp_at_10')
arms = ['A1', 'L025', 'L050', 'A3']; cs = np.array([0, .25, .5, 1.])
roi = tp.groupby(level=0)[arms].sum() / 60 * 100     # P@10 points per ROI
for c in (0.9, 0.8, 0.7):
    own = roi['A3'] - np.array([np.interp(c, cs, r) for r in roi.to_numpy()])
    print(f'uniform c={c}: pooled PL {own.mean():.3f} pt; own-ROI losses SD {own.std(ddof=1):.3f} pt -> SE of a 49-ROI mean {own.std(ddof=1) / 7:.3f} pt')
# does the per-ROI 0.5->1.0 response correlate with white share or clamp-related quantities?
fr = pd.read_csv(f'{REPO}/floor_robustness/floor_per_roi_by_scanner.csv').set_index('file_name')
step = (roi['A3'] - roi['L050'])
for col in ('white', 'b_otsu', 'TOTSU_med', 'W128_med'):
    r = stats.spearmanr(step, fr.loc[step.index, col])
    print(f'Spearman(0.5->1.0 gain, {col}) = {r.correlation:+.2f} (p {r.pvalue:.2f})')
