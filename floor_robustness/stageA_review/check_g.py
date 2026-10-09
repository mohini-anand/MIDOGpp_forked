"""5g: the realistic perturbations (plus extras the plan omits) on one ROI; c_eq, Rres vs the 0.75x/0.9x refs, for b_otsu, TOTSU, W128, TW128."""
import sys, time
import numpy as np, cv2
from skimage.filters import threshold_multiotsu
from common import *

fn = sys.argv[1]
roi_index = FILES.index(fn)
rgb, hem, mpp = load(fn)
gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY); T0 = gray <= 215
b0 = roi_floors(hem)['b_otsu']
hT = hem[T0]; v0 = np.sort(hT); s0 = hT < np.float32(b0)
ref = {c: rres(hT, s0, float(np.float32(c * b0))) for c in (0.5, 0.75, 0.9)}
M = macenko_basis(rgb, np.random.default_rng([20261008, roi_index]))
side = int(round(128 / mpp)) | 1
clk = CLK[CLK.file_name == fn]
H, W = hem.shape
print(f'{fn}: b0 {b0:.5f}  Rres refs 0.5x {ref[0.5]:.5f} 0.75x {ref[0.75]:.5f} 0.9x {ref[0.9]:.5f}  clamp share (hem0==0) {np.mean(hem == 0):.3f}, tissue {np.mean(hT == 0):.3f}')


def synth_pos(rgb, M, k, rows=256):
    """Variant: scale only the non-negative part of each concentration."""
    pinvM = np.linalg.pinv(M); out = np.empty_like(rgb); k = np.asarray(k, float)
    for r0 in range(0, rgb.shape[0], rows):
        blk = rgb[r0:r0 + rows].reshape(-1, 3); O = od(blk); C = O @ pinvM.T; R = O - C @ M.T
        C1 = C + (k - 1) * np.maximum(C, 0)
        out[r0:r0 + rows] = to_rgb(C1 @ M.T + R).reshape(rgb[r0:r0 + rows].shape)
    return out


def grad(rgb, rows=256):
    out = np.empty_like(rgb)
    long_x = W >= H
    for r0 in range(0, H, rows):
        blk = rgb[r0:r0 + rows].astype(np.float64)
        O = -np.log10((blk + 1) / 256.0)
        if long_x:
            g = 0.8 + 0.4 * np.arange(W) / (W - 1); O = O * g[None, :, None]
        else:
            g = 0.8 + 0.4 * (np.arange(r0, r0 + blk.shape[0]) / (H - 1)); O = O * g[:, None, None]
        out[r0:r0 + rows] = to_rgb(O)
    return out


def lut(f):
    t = np.clip(np.rint(f(np.arange(256, dtype=np.float64))), 0, 255).astype(np.uint8)
    return cv2.LUT(rgb, t)


def jpeg(q, samp=cv2.IMWRITE_JPEG_SAMPLING_FACTOR_420):
    ok, buf = cv2.imencode('.jpg', cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, q, cv2.IMWRITE_JPEG_SAMPLING_FACTOR, samp])
    return cv2.cvtColor(cv2.imdecode(buf, cv2.IMREAD_COLOR), cv2.COLOR_BGR2RGB)


def gain(gr, gg, gb):
    g = np.array([gr, gg, gb])
    return np.clip(np.rint(rgb.astype(np.float32) * g), 0, 255).astype(np.uint8)


PERTS = {
    'H06': lambda: synth(rgb, M, [0.6, 1]),
    'H06pos': lambda: synth_pos(rgb, M, [0.6, 1]),
    'H14': lambda: synth(rgb, M, [1.4, 1]),
    'E15': lambda: synth(rgb, M, [1, 1.5]),
    'CDN': lambda: lut(lambda I: 12 + I * 231 / 255),
    'CUP': lambda: lut(lambda I: (I - 12) * 255 / 231),
    'BLUR': lambda: cv2.GaussianBlur(rgb, (0, 0), 1.0, borderType=cv2.BORDER_REFLECT_101),
    'JPG': lambda: jpeg(75),
    'JPG95': lambda: jpeg(95),
    'GRAD': lambda: grad(rgb),
    'WB_R0.95': lambda: gain(0.95, 1, 1),
    'WB_R1.05': lambda: gain(1.05, 1, 1),
    'GAIN0.9': lambda: gain(0.9, 0.9, 0.9),
}
only = sys.argv[2].split(',') if len(sys.argv) > 2 else list(PERTS)
print(f"{'pert':9s} {'est':6s} {'b/b0':>6s} {'c_eq':>6s} {'Rres':>8s} {'R-ref75':>8s} {'R-ref90':>8s}  T1<.75 T3>ref75 stop(c>=.9&R<=ref90)")
for name in only:
    t = time.time()
    rgb1 = PERTS[name]()
    hem1 = ch.to_channel(rgb1, 'hematoxylin_od'); gray1 = cv2.cvtColor(rgb1, cv2.COLOR_RGB2GRAY)
    ht = hem1[T0]
    est = {}
    est['b_otsu'] = [float(np.float32(threshold_multiotsu(hem1, classes=3)[0]))]
    est['TOTSU'] = [float(np.float32(threshold_multiotsu(hem1[gray1 <= 215], classes=3)[0]))]
    est['W128'], est['TW128'] = [], []
    for r in clk.itertuples():
        est['W128'].append(window_floor(hem1, r.click_cx, r.click_cy, side)[0])
        sl = window_slice(hem1.shape, r.click_cx, r.click_cy, side)
        tw = hem1[sl][gray1[sl] <= 215]
        est['TW128'].append(float(np.float32(threshold_multiotsu(tw, classes=3)[0])) if tw.size >= 0.05 * hem1[sl].size else np.nan)
    for e, bs in est.items():
        ce = np.median([c_eq(ht, v0, b0, b) for b in bs]); rr = np.median([rres(ht, s0, b) for b in bs])
        print(f"{name:9s} {e:6s} {np.median(bs) / b0:6.3f} {ce:6.3f} {rr:8.5f} {rr - ref[0.75]:+8.5f} {rr - ref[0.9]:+8.5f}  {str(ce < 0.75):5s}  {str(rr - ref[0.75] > 0):5s}  {str(ce >= 0.9 and rr - ref[0.9] <= 0)}")
    # how much of Rres is intrinsic to the perturbation: Rres at the quantile-matched floor (the b' that clips exactly the original share p0 on T0)
    p0 = s0.mean(); b_match = float(np.quantile(ht, p0))
    print(f"{name:9s} {'QMATCH':6s} {b_match / b0:6.3f} {c_eq(ht, v0, b0, b_match):6.3f} {rres(ht, s0, b_match):8.5f} {rres(ht, s0, b_match) - ref[0.75]:+8.5f} {rres(ht, s0, b_match) - ref[0.9]:+8.5f}   (floor placed at the original tissue share p0={p0:.3f}; {time.time() - t:.1f}s)")
