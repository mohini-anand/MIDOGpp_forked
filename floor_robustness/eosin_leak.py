"""How much eosin-rich, non-nuclear tissue sits above the floor, per ROI (2026-10-02, descriptive)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from skimage.color import rgb2hed
from skimage.filters import threshold_multiotsu
REPO = Path('/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked')
sys.path.insert(0, str(REPO))
from midog_utils_full import dataset as ds

fl = pd.read_csv(REPO / 'results/tm_floor_definition_49roi_6seed_floors.csv')
b_otsu = fl[fl.arm == 'A3'].groupby('file_name').b.first()
rng = np.random.default_rng(0)
rows = []
for fn, b in b_otsu.items():
    p = REPO / 'images/extra_valid' / fn
    if not p.exists():
        p = REPO / 'images/extra_valid/testing_set' / fn
    rgb = ds.load_roi(p)
    gray = (0.299 * rgb[..., 0] + 0.587 * rgb[..., 1] + 0.114 * rgb[..., 2])
    idx = np.flatnonzero((gray <= 215).ravel())
    idx = rng.choice(idx, size=min(400_000, idx.size), replace=False)
    px = rgb.reshape(-1, 3)[idx].astype(np.float32) / 255.0
    hed = rgb2hed(px[None])[0]
    H, E = hed[:, 0], hed[:, 1]
    t_lo, t_hi = threshold_multiotsu(H, classes=3)
    nonnuc = H < t_hi
    e_hi = E > np.quantile(E[nonnuc], 0.75)
    sel = nonnuc & e_hi
    sel_lo = nonnuc & (E < np.quantile(E[nonnuc], 0.25))
    rows.append(dict(file_name=fn, b_otsu=b, h_eosin_rich_med=float(np.median(H[sel])), h_eosin_poor_med=float(np.median(H[sel_lo])),
                     share_eosin_rich_above_b=float((H[sel] > b).mean()), share_eosin_rich_above_1p5b=float((H[sel] > 1.5 * b).mean()),
                     share_eosin_poor_above_b=float((H[sel_lo] > b).mean()),
                     rho_HE_nonnuc=float(pd.Series(H[nonnuc]).corr(pd.Series(E[nonnuc]), method='spearman'))))
    print(fn, flush=True)
pd.DataFrame(rows).to_csv(Path(sys.argv[1]), index=False)
