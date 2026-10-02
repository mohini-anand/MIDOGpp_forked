"""Per-ROI Macenko stain vectors vs the fixed Ruifrok vectors rgb2hed uses, and the eosin->H leakage they imply."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from skimage.color import hed_from_rgb, rgb_from_hed
REPO = Path('/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked')
sys.path.insert(0, str(REPO))
from midog_utils_full import dataset as ds

BETA, ALPHA = 0.15, 1.0   # Macenko 2009 Algorithm 1 constants
fl = pd.read_csv(REPO / 'results/tm_floor_definition_49roi_6seed_floors.csv')
files = sorted(fl.file_name.unique())
rng = np.random.default_rng(0)
ref_h = rgb_from_hed[0] / np.linalg.norm(rgb_from_hed[0])
ref_e = rgb_from_hed[1] / np.linalg.norm(rgb_from_hed[1])


def ang(a, b):
    return float(np.degrees(np.arccos(np.clip(abs(a @ b) / np.linalg.norm(a) / np.linalg.norm(b), -1, 1))))


rows = []
for fn in files:
    p = REPO / 'images/extra_valid' / fn
    if not p.exists():
        p = REPO / 'images/extra_valid/testing_set' / fn
    rgb = ds.load_roi(p).reshape(-1, 3)
    px = rgb[rng.choice(len(rgb), 500_000, replace=False)].astype(np.float64)
    od = -np.log10((px + 1) / 256.0)
    od = od[(od > BETA).all(1)]
    _, _, vt = np.linalg.svd(od - 0 * od.mean(0), full_matrices=False)
    v = vt[:2].T                                  # 3x2 plane of the two largest directions
    if v[:, 0].sum() < 0:
        v[:, 0] *= -1
    if v[:, 1].sum() < 0:
        v[:, 1] *= -1
    proj = od @ v
    phi = np.arctan2(proj[:, 1], proj[:, 0])
    lo, hi = np.percentile(phi, [ALPHA, 100 - ALPHA])
    s1 = v @ np.array([np.cos(lo), np.sin(lo)])
    s2 = v @ np.array([np.cos(hi), np.sin(hi)])
    s1, s2 = s1 / np.linalg.norm(s1), s2 / np.linalg.norm(s2)
    h_est, e_est = (s1, s2) if s1[0] > s2[0] else (s2, s1)   # hematoxylin absorbs red more than eosin does
    # what rgb2hed's fixed matrix does with the slide's own pure stains: H-channel response per unit OD
    h_resp_h = float((h_est @ hed_from_rgb)[0])
    h_resp_e = float((e_est @ hed_from_rgb)[0])
    rows.append(dict(file_name=fn, h_est=np.round(h_est, 3).tolist(), e_est=np.round(e_est, 3).tolist(), ang_h_deg=ang(h_est, ref_h), ang_e_deg=ang(e_est, ref_e),
                     h_channel_per_unit_true_h=h_resp_h, h_channel_per_unit_true_e=h_resp_e, eosin_leak_ratio=h_resp_e / h_resp_h))
    print(fn, flush=True)
pd.DataFrame(rows).to_csv(Path(sys.argv[1]), index=False)
