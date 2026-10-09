"""PL (plan v3 section 3.5) on the 7-ROI pilot dumps, plus checks of the PL construct against the ladder run's per-click data."""
import glob
import numpy as np, pandas as pd

REPO = '/Users/mohinianand/Library/Mobile Documents/com~apple~CloudDocs/Desktop/AnnotateDx/MIDOGpp_forked'
lad = pd.read_csv(f'{REPO}/results/tm_floor_definition_49roi_6seed_ladder.csv')
CS = np.array([0, 0.25, 0.5, 1.0]); PS = lad.set_index('arm').loc[['A1', 'L025', 'L050', 'A3'], 'P@10'].to_numpy()
print('ladder P@10 at c=0/.25/.5/1:', np.round(PS, 6))


def loss(c):
    """Plan v3: max(0, L(1)-L(c)) for c<=1, 0 for 1<c<=2; c>2 undefined in the plan -> NaN here so it shows up."""
    c = np.asarray(c, float)
    out = np.where(c <= 1, np.maximum(0, PS[-1] - np.interp(np.minimum(c, 1), CS, PS)), 0.0)
    return np.where(c > 2, np.nan, out) * 100      # in P@10 points


print(f'uniform 0.91 -> {loss(0.91):.3f} pt; 0.9 -> {loss(0.9):.3f}; 0.5 -> {loss(0.5):.3f}')

# ---- pilot PL (7 ROIs x 6 clicks = 42 landings)
D = pd.concat([pd.read_csv(f) for f in sorted(glob.glob('v3_dump_*.csv'))])
print(f'\npilot: {D.file_name.nunique()} ROIs, {D.groupby("var").size().iloc[0]} landings per variant -- NOT 49 ROIs')
order = ['UNPERT', 'GAMMA', 'JPG95', 'H06', 'H14', 'E15', 'CDN', 'CUP', 'GAIN09', 'GAIN11', 'WBR095', 'WBR105', 'WBB095', 'WBB105', 'BLUR', 'JPG', 'GRAD', 'SHIFT']
rows = []
for var in order:
    g = D[D['var'] == var]
    for est, col in (('b_otsu', 'c_b_otsu'), ('W128', 'c_W128')):
        if var == 'SHIFT' and est == 'b_otsu':
            continue
        c = g[col].to_numpy()
        pl = np.nanmean(loss(c)); n_nan = int(np.isnan(loss(c)).sum())
        roi_val = g.groupby('file_name')[col].median() if est == 'W128' else g.groupby('file_name')[col].first()
        over = (roi_val.median() > 1.5) or (roi_val > 2.0).any()
        rows.append(dict(var=var, est=est, PL=pl, PL_fires=pl >= 0.25, OVER=over, roi_median=roi_val.median(), roi_min=roi_val.min(), roi_max=roi_val.max(), share_lt_09=float((c < 0.9).mean()), n_gt2=n_nan))
T = pd.DataFrame(rows)
pd.set_option('display.width', 200)
print(T.to_string(index=False, float_format=lambda v: f'{v:.3f}'))

# ---- unperturbed W128 PL on all 294 clicks (from floors.csv; unperturbed c_eq == b/b_otsu to <3e-5)
fl = pd.read_csv(f'{REPO}/results/tm_floor_definition_49roi_6seed_floors.csv', float_precision='round_trip')
w = fl[fl.arm == 'W128'].sort_values(['file_name', 'seed_index'])
print(f'\nunperturbed W128 PL, all 49 ROIs x 6 clicks: {np.mean(loss(w.b_over_otsu)):.3f} pt  (share of clicks < 1: {(w.b_over_otsu < 1).mean():.2f}; this is also NULL2 PL for W128)')
tt = fl[fl.arm == 'TOTSU']; print(f'unperturbed TOTSU PL: {np.mean(loss(tt.b_over_otsu)):.3f}')

# ---- does PL predict a measured outcome? the W128 arm searched with per-click landings
pr = pd.read_csv(f'{REPO}/results/tm_floor_definition_49roi_6seed_per_run.csv')
tp = pr.pivot_table(index=['file_name', 'seed_index'], columns='arm', values='tp_at_10')
w = w.set_index(['file_name', 'seed_index'])
obs = (tp['W128'] - tp['A3']).mean() * 10    # P@10 points (tp/10 *100)
print(f'W128 - A3 observed: {obs:+.3f} pt; pooled-curve PL predicts {-np.mean(loss(w.b_over_otsu)):+.3f} pt (no overshoot credit)')
# own-curve prediction, per click, with the plan's no-credit rule and with full own curve
cs_all = np.array([0, 0.25, 0.5, 1.0, 1.5, 2.0]); arms = ['A1', 'L025', 'L050', 'A3', 'L150', 'L200']
cur = tp[arms].to_numpy(float)
cw = w.reindex(tp.index).b_over_otsu.to_numpy()
own_full = np.array([np.interp(c, cs_all, row) for c, row in zip(cw, cur)]) - tp['A3'].to_numpy()
own_nocredit = np.where(cw <= 1, own_full, 0.0)
print(f'own-click-curve prediction: full {own_full.mean() * 10:+.3f} pt; no-overshoot-credit {own_nocredit.mean() * 10:+.3f} pt')

# ---- heterogeneity: per-ROI 0.5->1.0 response (tp sum A3 - L050, per 6 clicks) and the pilot ROIs' own curves
roi = tp.groupby(level=0).sum()
step = (roi['A3'] - roi['L050'])
print(f'\nper-ROI TP@10 gain 0.5->1.0 (sum of 6 clicks): mean {step.mean():.2f}, sd {step.std():.2f}, range {step.min():.0f}..{step.max():.0f}')
pil = ['401.tiff', '343.tiff', '201.tiff', '529.tiff', '128.tiff', '211.tiff', '403.tiff']
print('pilot ROIs 0.5->1.0 gain:', step.reindex(pil).to_dict())
# own-ROI-curve PL vs pooled PL for the pilot H06/E15 b_otsu landings (no max clip; c<=1 interpolated on the ROI's own summed curve)
for var in ('H06', 'E15', 'CDN', 'GAIN09'):
    g = D[D['var'] == var].groupby('file_name').c_b_otsu.first()
    own = []
    for f, c in g.items():
        rc = roi.loc[f, arms].to_numpy(float) / 60 * 100     # P@10 points for the ROI
        own.append(rc[3] - np.interp(min(c, 1), cs_all, rc) if c <= 1 else 0.0)
    print(f'{var} b_otsu, 7 ROIs: pooled-curve PL {np.mean(loss(g.to_numpy())):.3f} pt; own-ROI-curve PL {np.mean(own):.3f} pt (per ROI {np.round(own, 2).tolist()})')
