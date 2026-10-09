"""T2 noise check: multi-Otsu bin width relative to the floor (ROI-wide, W128 windows, TW128 windows) and per-click c_eq jitter under near-null JPEG q95, plus SHIFT."""
import sys, time
from common import *
import numpy as np, cv2
from skimage.filters import threshold_multiotsu


def jpeg(rgb, q):
    ok, buf = cv2.imencode('.jpg', cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, q, cv2.IMWRITE_JPEG_SAMPLING_FACTOR, cv2.IMWRITE_JPEG_SAMPLING_FACTOR_420])
    return cv2.cvtColor(cv2.imdecode(buf, cv2.IMREAD_COLOR), cv2.COLOR_BGR2RGB)


def tw_floor(hem, gray, cx, cy, side):
    sl = window_slice(hem.shape, cx, cy, side)
    w = hem[sl]; m = gray[sl] <= 215
    if m.sum() < 0.05 * w.size:
        return np.nan, np.nan
    tw = w[m]
    b = float(np.float32(threshold_multiotsu(tw, classes=3)[0]))
    return b, (float(tw.max()) - float(tw.min())) / 256 / b


for fn in sys.argv[1:]:
    rgb, hem, mpp = load(fn)
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY); T0 = gray <= 215
    b0 = roi_floors(hem)['b_otsu']
    hT = hem[T0]; v0 = np.sort(hT)
    side = int(round(128 / mpp)) | 1
    roi_bin = (float(hem.max()) - float(hem.min())) / 256 / b0
    rgbq = jpeg(rgb, 95); hq = ch.to_channel(rgbq, 'hematoxylin_od'); gq = cv2.cvtColor(rgbq, cv2.COLOR_RGB2GRAY); htq = hq[T0]
    clk = CLK[CLK.file_name == fn]
    rows = []
    for r in clk.itertuples():
        sl = window_slice(hem.shape, r.click_cx, r.click_cy, side)
        w = hem[sl]
        bw = window_floor(hem, r.click_cx, r.click_cy, side)[0]
        bwq = window_floor(hq, r.click_cx, r.click_cy, side)[0]
        btw, tbin = tw_floor(hem, gray, r.click_cx, r.click_cy, side)
        btwq, _ = tw_floor(hq, gq, r.click_cx, r.click_cy, side)
        rows.append((r.seed_index, bw / b0, (float(w.max()) - float(w.min())) / 256 / bw, c_eq(hT, v0, b0, bw), c_eq(htq, v0, b0, bwq), btw / b0, tbin, c_eq(hT, v0, b0, btw), c_eq(htq, v0, b0, btwq), T0[sl].mean()))
    a = np.array(rows, float)
    print(f'\n{fn}: b0 {b0:.4f}  hem max {hem.max():.3f}  ROI bin/b0 {roi_bin:.3f}')
    print(' seed  W128/b0  Wbin/b  c_eq   c_eq(q95)  drop | TW128/b0  TWbin/b  c_eq   c_eq(q95)  drop | tissue')
    for x in a:
        print(f'  s{int(x[0])}  {x[1]:.3f}   {x[2]:.3f}   {x[3]:.3f}  {x[4]:.3f}    {x[3]-x[4]:+.3f} | {x[5]:.3f}    {x[6]:.3f}    {x[7]:.3f}  {x[8]:.3f}    {x[7]-x[8]:+.3f} | {x[9]:.2f}')
    print(f'  ROI medians: W128 c_eq {np.median(a[:,3]):.3f} -> q95 {np.median(a[:,4]):.3f};  TW128 {np.nanmedian(a[:,7]):.3f} -> {np.nanmedian(a[:,8]):.3f}')
    # SHIFT: 8 offsets at 64 um (unit directions incl. diagonals), choose highest glass share on the original gray
    off = 64 / mpp
    dirs = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    sh = []
    for r in clk.itertuples():
        best = None
        for dx, dy in dirs:
            n = np.hypot(dx, dy)
            cx, cy = r.click_cx + off * dx / n, r.click_cy + off * dy / n
            sl = window_slice(hem.shape, cx, cy, side)
            g = 1 - T0[sl].mean()
            if best is None or g > best[0]:
                best = (g, cx, cy, T0[sl].size / side ** 2)
        g, cx, cy, area = best
        bw = window_floor(hem, cx, cy, side)[0]
        btw, _ = tw_floor(hem, gray, cx, cy, side)
        sh.append((g, area, c_eq(hT, v0, b0, bw), c_eq(hT, v0, b0, btw) if np.isfinite(btw) else np.nan))
    sh = np.array(sh)
    print('  SHIFT per click (glass share, area share, W128 c_eq, TW128 c_eq):', np.round(sh, 3).tolist())
    print(f'  SHIFT ROI medians: W128 {np.median(sh[:,2]):.3f}  TW128 {np.nanmedian(sh[:,3]):.3f}')
