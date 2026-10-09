"""Plan v2 variants on one ROI: per-ROI c_eq for b_otsu (single) and W128 (median of 6 clicks), R2 status, per-variant timing, peak RSS."""
import sys, time, resource
from v2lib import *
from skimage.filters import threshold_multiotsu

fn = sys.argv[1]
only = sys.argv[2].split(',') if len(sys.argv) > 2 else None
t_all = time.time()
t = time.time()
rgb, hem, mpp = load(fn)
gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY); T0 = gray <= 215
b0 = roi_floors(hem)['b_otsu']
hT = hem[T0]; v0 = np.sort(hT); del hT
M = macenko_basis(rgb, np.random.default_rng([20261008, FILES.index(fn)]))
side = int(round(128 / mpp)) | 1
clk = CLK[CLK.file_name == fn]
t_setup = time.time() - t
flat = rgb.reshape(-1, 3)
zero = (flat == 0).any(1)
print(f'{fn}: shape {rgb.shape[:2]}  b0 {b0:.5f}  ROI bin/b0 {(float(hem.max()) - float(hem.min())) / 256 / b0:.3f}  hem max {hem.max():.3f}  zero-channel px: {zero.mean():.5f} of ROI, {(zero & T0.reshape(-1)).sum() / T0.sum():.5f} of tissue  tissue clamp share {np.mean(hem[T0] == 0):.3f}  setup {t_setup:.1f}s')


def floors(h1):
    bo = float(np.float32(threshold_multiotsu(h1, classes=3)[0]))
    ws = [window_floor(h1, r.click_cx, r.click_cy, side)[0] for r in clk.itertuples()]
    return bo, ws


bo, ws = floors(hem)
ht0 = hem[T0]
w0 = np.median([c_eq(ht0, v0, b0, w) for w in ws])
print(f'  unperturbed: b_otsu {c_eq(ht0, v0, b0, bo):.3f}  W128 per-ROI median {w0:.3f}  (per click {np.round([c_eq(ht0, v0, b0, w) for w in ws], 3).tolist()})')
del ht0
V = variants(rgb, M)
for name, fn_v in V.items():
    if only and name not in only:
        continue
    t = time.time(); rgb1 = fn_v(); t_make = time.time() - t
    t = time.time(); hem1 = ch.to_channel(rgb1, 'hematoxylin_od'); t_hem = time.time() - t
    t = time.time(); ht = hem1[T0]; bo, ws = floors(hem1)
    cb = c_eq(ht, v0, b0, bo); cw_click = [c_eq(ht, v0, b0, w) for w in ws]; cw = float(np.median(cw_click)); t_est = time.time() - t
    flag = lambda c: 'R2' if (c < 0.5 or c > 2.0) else ('out' if (c < 0.9 or c > 1.5) else '')
    print(f'  {name:7s} b_otsu {cb:6.3f} {flag(cb):3s} | W128 {cw:6.3f} {flag(cw):3s} (clicks {min(cw_click):.3f}-{max(cw_click):.3f}) | make {t_make:.2f}s hem {t_hem:.2f}s est+c_eq {t_est:.2f}s')
    del rgb1, hem1, ht
print(f'total {time.time() - t_all:.1f}s  peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9:.2f} GB')
