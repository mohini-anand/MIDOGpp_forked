"""Shared helpers for the Stage A plan review (scratch only; reads the repo, writes nothing there)."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '1')
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import cv2
from skimage.filters import threshold_multiotsu

REPO = Path('/Users/mohinianand/Library/Mobile Documents/com~apple~CloudDocs/Desktop/AnnotateDx/MIDOGpp_forked')
sys.dont_write_bytecode = True
sys.path.insert(0, str(REPO))
from midog_utils import dataset as ds
from midog_utils import channels as ch

FLOORS = pd.read_csv(REPO / 'results/tm_floor_definition_49roi_6seed_floors.csv')
FILES = sorted(FLOORS.file_name.unique())
B_REC = FLOORS[FLOORS.arm == 'A3'].groupby('file_name').b.first()
W128_REC = FLOORS[FLOORS.arm == 'W128'].set_index(['file_name', 'seed_index']).b
PR = pd.read_csv(REPO / 'results/tm_border_floor_49roi_6seed_per_run.csv')
CLK = PR[(PR.arm == 'A0') & (PR.rank_key == 'tm_score')][['file_name', 'seed_index', 'click_cx', 'click_cy', 'tcx', 'tcy', 'bs', 'mpp']].sort_values(['file_name', 'seed_index'])


def path_of(fn):
    p = REPO / 'images/extra_valid' / fn
    return p if p.exists() else REPO / 'images/extra_valid/testing_set' / fn


def load(fn):
    p = path_of(fn)
    rgb = ds.load_roi(p)
    return rgb, ch.to_channel(rgb, 'hematoxylin_od'), ds.roi_mpp(p)


def roi_floors(hem):
    th = threshold_multiotsu(hem, classes=3)
    return dict(b_otsu=float(np.float32(th[0])), b_otsu_upper=float(np.float32(th[1])), b_med=float(np.float32(np.median(hem))))


def window_floor(hem, cx, cy, side):
    h, w = hem.shape
    half = side // 2
    ix, iy = int(round(cx)), int(round(cy))
    win = hem[max(0, iy - half):min(h, iy + half + 1), max(0, ix - half):min(w, ix + half + 1)]
    return float(np.float32(threshold_multiotsu(win, classes=3)[0])), int(win.size)


def window_slice(shape, cx, cy, side):
    h, w = shape
    half = side // 2
    ix, iy = int(round(cx)), int(round(cy))
    return slice(max(0, iy - half), min(h, iy + half + 1)), slice(max(0, ix - half), min(w, ix + half + 1))


def c_eq(ht, v0, b0, b):
    p = np.count_nonzero(ht < np.float32(b)) / ht.size
    return float(v0[min(int(round(p * ht.size)), ht.size - 1)]) / b0


def rres(ht, s0, b):
    r = np.maximum(ht - np.float32(b), 0, dtype=np.float32)
    tot = float(r.sum(dtype=np.float64))
    return float(r[s0].sum(dtype=np.float64)) / tot


def od(I):
    return -np.log10((I.astype(np.float64) + 1) / 256.0)


def to_rgb(o):
    return np.clip(np.rint(256.0 * 10.0 ** (-o) - 1), 0, 255).astype(np.uint8)


def macenko_basis(rgb, rng, beta=0.15, alpha=1.0, n=500_000):
    flat = rgb.reshape(-1, 3)
    px = flat[rng.choice(len(flat), n, replace=False)].astype(np.float64)
    o = -np.log10((px + 1) / 256.0)
    o = o[(o > beta).all(1)]
    _, _, vt = np.linalg.svd(o, full_matrices=False)
    v = vt[:2].T
    if v[:, 0].sum() < 0:
        v[:, 0] *= -1
    if v[:, 1].sum() < 0:
        v[:, 1] *= -1
    proj = o @ v
    phi = np.arctan2(proj[:, 1], proj[:, 0])
    lo, hi = np.percentile(phi, [alpha, 100 - alpha])
    s1 = v @ np.array([np.cos(lo), np.sin(lo)])
    s2 = v @ np.array([np.cos(hi), np.sin(hi)])
    s1, s2 = s1 / np.linalg.norm(s1), s2 / np.linalg.norm(s2)
    h, e = (s1, s2) if s1[0] > s2[0] else (s2, s1)
    return np.stack([h, e], 1)


def synth(rgb, M, k, rows=256, g=None):
    """OD' = (C*k) @ M.T + R in row blocks; g: optional per-row-block multiplier function (for GRAD) not used here."""
    pinvM = np.linalg.pinv(M)
    out = np.empty_like(rgb)
    k = np.asarray(k, float)
    for r0 in range(0, rgb.shape[0], rows):
        blk = rgb[r0:r0 + rows].reshape(-1, 3)
        O = od(blk)
        C = O @ pinvM.T
        R = O - C @ M.T
        O1 = (C * k) @ M.T + R
        out[r0:r0 + rows] = to_rgb(O1).reshape(rgb[r0:r0 + rows].shape)
    return out
