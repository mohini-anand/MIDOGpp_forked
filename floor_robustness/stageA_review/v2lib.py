"""Plan v2 synthesis and variants, as written in TM_FLOOR_STRESS_STAGE_A_PLAN.md v2 (scratch re-implementation for review)."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'
from common import *          # noqa: E402  (common imports numpy after the env is set)
import numpy as np, cv2       # noqa: E402
cv2.setNumThreads(1)
from skimage.color import hed_from_rgb  # noqa: E402

HCOL = hed_from_rgb[:, 0]


def synth_v2(rgb, M, kH, kE, rows=256, col_for_kH=0):
    """OD' = OD + (kH-1)*max(cH,0)*h + (kE-1)*max(cE,0)*e; pixels with any channel 0 copied. col_for_kH=1 simulates a wrong-column bug."""
    pinvM = np.linalg.pinv(M)
    out = np.empty_like(rgb)
    h, e = M[:, 0], M[:, 1]
    for r0 in range(0, rgb.shape[0], rows):
        blk = rgb[r0:r0 + rows]
        flat = blk.reshape(-1, 3)
        O = od(flat)
        C = O @ pinvM.T
        cpos = np.maximum(C, 0)
        O1 = O + (kH - 1) * cpos[:, [col_for_kH]] * h[None, :] + (kE - 1) * cpos[:, [1]] * e[None, :]
        o = to_rgb(O1)
        zero = (flat == 0).any(1)
        o[zero] = flat[zero]
        out[r0:r0 + rows] = o.reshape(blk.shape)
    return out


def grad_v2(rgb, rows=256):
    H, W = rgb.shape[:2]
    out = np.empty_like(rgb)
    long_x = W >= H
    for r0 in range(0, H, rows):
        blk = rgb[r0:r0 + rows]
        O = od(blk)
        if long_x:
            O = O * (0.8 + 0.4 * np.arange(W) / (W - 1))[None, :, None]
        else:
            O = O * (0.8 + 0.4 * np.arange(r0, r0 + blk.shape[0]) / (H - 1))[:, None, None]
        o = to_rgb(O)
        zero = (blk == 0).any(-1)
        o[zero] = blk[zero]
        out[r0:r0 + rows] = o
    return out


def lut_all(rgb, f):
    t = np.clip(np.rint(f(np.arange(256, dtype=np.float64))), 0, 255).astype(np.uint8)
    return cv2.LUT(rgb, t)


def lut_one(rgb, ch_i, f):
    t = np.clip(np.rint(f(np.arange(256, dtype=np.float64))), 0, 255).astype(np.uint8)
    out = rgb.copy()
    out[..., ch_i] = cv2.LUT(np.ascontiguousarray(rgb[..., ch_i]), t)
    return out


def jpeg(rgb, q):
    ok, buf = cv2.imencode('.jpg', cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, q, cv2.IMWRITE_JPEG_SAMPLING_FACTOR, cv2.IMWRITE_JPEG_SAMPLING_FACTOR_420])
    return cv2.cvtColor(cv2.imdecode(buf, cv2.IMREAD_COLOR), cv2.COLOR_BGR2RGB)


def variants(rgb, M):
    return {
        'GAMMA': lambda: lut_all(rgb, lambda I: 255 * (I / 255) ** 1.1),
        'JPG95': lambda: jpeg(rgb, 95),
        'H06': lambda: synth_v2(rgb, M, 0.6, 1.0),
        'H14': lambda: synth_v2(rgb, M, 1.4, 1.0),
        'E15': lambda: synth_v2(rgb, M, 1.0, 1.5),
        'CDN': lambda: lut_all(rgb, lambda I: 12 + I * 231 / 255),
        'CUP': lambda: lut_all(rgb, lambda I: (I - 12) * 255 / 231),
        'GAIN09': lambda: lut_all(rgb, lambda I: 0.9 * I),
        'GAIN11': lambda: lut_all(rgb, lambda I: 1.1 * I),
        'WBR095': lambda: lut_one(rgb, 0, lambda I: 0.95 * I),
        'WBR105': lambda: lut_one(rgb, 0, lambda I: 1.05 * I),
        'WBB095': lambda: lut_one(rgb, 2, lambda I: 0.95 * I),
        'WBB105': lambda: lut_one(rgb, 2, lambda I: 1.05 * I),
        'BLUR': lambda: cv2.GaussianBlur(rgb, (0, 0), 1.0, borderType=cv2.BORDER_REFLECT_101),
        'JPG': lambda: jpeg(rgb, 75),
        'GRAD': lambda: grad_v2(rgb),
    }
