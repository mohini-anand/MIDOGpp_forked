"""Per-channel tone-curve (scanner colour profile) perturbations the plan omits: OD of one channel scaled via RGB gamma. Also NULL2 Rres exactness."""
import sys
from common import *
import numpy as np, cv2
from skimage.filters import threshold_multiotsu

for fn in sys.argv[1:]:
    rgb, hem, mpp = load(fn)
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY); T0 = gray <= 215
    b0 = roi_floors(hem)['b_otsu']
    hT = hem[T0]; v0 = np.sort(hT); s0 = hT < np.float32(b0)
    ref = {c: rres(hT, s0, float(np.float32(c * b0))) for c in (0.75, 0.9)}
    side = int(round(128 / mpp)) | 1
    clk = CLK[CLK.file_name == fn]
    # NULL2 Rres / c_eq exactness for b_otsu and W128
    h2 = hem * np.float32(2); ht2 = h2[T0]
    b2 = float(np.float32(threshold_multiotsu(h2, classes=3)[0]))
    ok = [c_eq(ht2, v0, b0, b2) == c_eq(hT, v0, b0, b0), rres(ht2, s0, b2) == rres(hT, s0, b0)]
    for r in clk.itertuples():
        bw = window_floor(hem, r.click_cx, r.click_cy, side)[0]; bw2 = window_floor(h2, r.click_cx, r.click_cy, side)[0]
        ok += [bw2 == 2 * bw, c_eq(ht2, v0, b0, bw2) == c_eq(hT, v0, b0, bw), rres(ht2, s0, bw2) == rres(hT, s0, bw)]
    print(f'{fn}: NULL2 exact on all {len(ok)} checks: {all(ok)};  refs 0.75x {ref[0.75]:.5f} 0.9x {ref[0.9]:.5f}')
    del h2, ht2
    for ch_i, name in ((2, 'B'), (0, 'R')):
        for gam in (0.77, 1.3):
            t = np.clip(np.rint(255.0 * (np.arange(256) / 255.0) ** gam), 0, 255).astype(np.uint8)
            rgb1 = rgb.copy(); rgb1[..., ch_i] = cv2.LUT(rgb[..., ch_i], t)
            hem1 = ch.to_channel(rgb1, 'hematoxylin_od'); ht = hem1[T0]
            b1 = float(np.float32(threshold_multiotsu(hem1, classes=3)[0]))
            ws = [window_floor(hem1, r.click_cx, r.click_cy, side)[0] for r in clk.itertuples()]
            cw = np.median([c_eq(ht, v0, b0, w) for w in ws]); rw = np.median([rres(ht, s0, w) for w in ws])
            print(f'  OD_{name} x{gam}: b_otsu c_eq {c_eq(ht, v0, b0, b1):.3f} Rres-ref90 {rres(ht, s0, b1) - ref[0.9]:+.5f} | W128 c_eq {cw:.3f} Rres-ref90 {rw - ref[0.9]:+.5f}')
