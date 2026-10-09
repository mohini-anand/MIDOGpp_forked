"""Per-click c_eq for b_otsu (ROI value) and W128 (own window) under every v2/v3 variant + unperturbed + SHIFT, one ROI -> CSV."""
import sys
import v2lib
from v2lib import *
from midog_utils import channels as _ch
from skimage.filters import threshold_multiotsu

_orig = _ch.to_channel
def _blocked(rgb, name, rows=1024):
    out = np.empty(rgb.shape[:2], np.float32)
    for r0 in range(0, rgb.shape[0], rows):
        out[r0:r0 + rows] = _orig(rgb[r0:r0 + rows], name)
    return out
ch.to_channel = _blocked

fn = sys.argv[1]
rgb, hem, mpp = load(fn)
gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY); T0 = gray <= 215
b0 = roi_floors(hem)['b_otsu']
v0 = np.sort(hem[T0])
M = macenko_basis(rgb, np.random.default_rng([20261008, FILES.index(fn)]))
side = int(round(128 / mpp)) | 1
clk = CLK[CLK.file_name == fn]
H, W = hem.shape
rows = []


def record(var, h1):
    ht = h1[T0]
    bo = float(np.float32(threshold_multiotsu(h1, classes=3)[0]))
    cb = c_eq(ht, v0, b0, bo)
    for r in clk.itertuples():
        bw = window_floor(h1, r.click_cx, r.click_cy, side)[0]
        rows.append(dict(file_name=fn, var=var, seed_index=int(r.seed_index), c_b_otsu=cb, c_W128=c_eq(ht, v0, b0, bw)))


record('UNPERT', hem)
# SHIFT per plan v3 section 3.3 (W128 only; b_otsu column = unperturbed 1.0)
d = int(round(64 / mpp))
dirs = [(d, 0), (-d, 0), (0, d), (0, -d), (d, d), (-d, d), (d, -d), (-d, -d)]
ht0 = hem[T0]
for r in clk.itertuples():
    best = None
    for dx, dy in dirs:
        cx = min(max(r.click_cx + dx, 0), W - 1); cy = min(max(r.click_cy + dy, 0), H - 1)
        sl = window_slice(hem.shape, cx, cy, side)
        g = float((gray[sl] > 215).mean())
        if best is None or g > best[0]:
            best = (g, cx, cy)
    bw = window_floor(hem, best[1], best[2], side)[0]
    rows.append(dict(file_name=fn, var='SHIFT', seed_index=int(r.seed_index), c_b_otsu=1.0, c_W128=c_eq(ht0, v0, b0, bw)))
for var, f in variants(rgb, M).items():
    record(var, ch.to_channel(f(), 'hematoxylin_od'))
pd.DataFrame(rows).to_csv(f'v3_dump_{fn}.csv', index=False)
print(fn, 'done', len(rows))
