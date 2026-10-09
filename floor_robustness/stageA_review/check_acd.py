"""5a NULL2 exactness, 5c c_eq self-consistency + Rres(b0), 5d unperturbed dry run, on the ROIs given as argv."""
import sys, time
import numpy as np, cv2
from skimage.filters import threshold_multiotsu
from common import *

for fn in sys.argv[1:]:
    t = time.time()
    rgb, hem, mpp = load(fn)
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    T0 = gray <= 215
    b0 = roi_floors(hem)['b_otsu']
    print(f'\n##### {fn} shape {hem.shape} mpp {mpp:.4f} b0 {b0!r} recorded {B_REC[fn]!r} equal={b0 == B_REC[fn]}  (load+hem {time.time()-t:.1f}s)')
    # --- 5a NULL2 exactness
    th1 = threshold_multiotsu(hem, classes=3)
    hem2 = hem * np.float32(2)
    th2 = threshold_multiotsu(hem2, classes=3)
    print('5a ROI: th dtype', th1.dtype, 'th2 == 2*th1 bitwise:', np.array_equal(th2, 2 * th1), th1, th2)
    S1 = hem < np.float32(th1[0]); S2 = hem2 < np.float32(th2[0])
    print('5a ROI clipped sets identical:', np.array_equal(S1, S2))
    side = int(round(128 / mpp)) | 1
    clk = CLK[CLK.file_name == fn]
    for r in clk.itertuples():
        sl = window_slice(hem.shape, r.click_cx, r.click_cy, side)
        w1 = threshold_multiotsu(hem[sl], classes=3); w2 = threshold_multiotsu(hem2[sl], classes=3)
        bw, _ = window_floor(hem, r.click_cx, r.click_cy, side)
        same = np.array_equal(w2, 2 * w1) and np.array_equal(hem[sl] < np.float32(w1[0]), hem2[sl] < np.float32(w2[0]))
        rec = W128_REC.get((fn, r.seed_index), np.nan)
        # window bin width relative to the floor
        win = hem[sl]; bwid = (float(win.max()) - float(win.min())) / 256
        print(f'5a W128 s{r.seed_index}: x2 exact={same}  W128={bw:.6f} rec={rec:.6f} eq={bw == rec}  ratio_to_b0={bw / b0:.3f}  bin/b={bwid / bw:.3f}  max={win.max():.3f}')
    del hem2, S1, S2
    # --- 5c c_eq self-consistency
    hT = hem[T0]
    v0 = np.sort(hT)
    s0 = hT < np.float32(b0)
    Nt = hT.size
    print(f'N_t={Nt}  tissue share={Nt / hem.size:.3f}  unique tissue values={np.unique(v0).size}')
    for c in (0.25, 0.5, 0.75, 0.9, 1.0, 1.5, 2.0):
        b = float(np.float32(c * b0))
        ce = c_eq(hT, v0, b0, b)
        print(f'5c c={c:<5} c_eq={ce:.6f}  |c_eq-c|={abs(ce - c):.2e}  Rres={rres(hT, s0, b):.6f}')
    print('5c Rres(b0) exactly 0:', rres(hT, s0, b0) == 0.0)
    # c_eq(b0) range check
    ce1 = c_eq(hT, v0, b0, b0)
    print(f'5c c_eq(b0)={ce1!r} in [1,1.01): {1 <= ce1 < 1.01}')
    # Rres reference levels
    ref = {c: rres(hT, s0, float(np.float32(c * b0))) for c in (0.5, 0.75, 0.9)}
    print('5d Rres refs:', {k: round(v, 5) for k, v in ref.items()})
    # --- 5d dry run on the original image: TOTSU, W128, TW128
    btot = float(np.float32(threshold_multiotsu(hem[T0], classes=3)[0]))
    print(f'5d TOTSU b={btot:.6f} b/b0={btot / b0:.3f} c_eq={c_eq(hT, v0, b0, btot):.3f} Rres={rres(hT, s0, btot):.5f}')
    for r in clk.itertuples():
        bw, _ = window_floor(hem, r.click_cx, r.click_cy, side)
        sl = window_slice(hem.shape, r.click_cx, r.click_cy, side)
        tw = hem[sl][T0[sl]]
        btw = float(np.float32(threshold_multiotsu(tw, classes=3)[0])) if tw.size >= 0.05 * hem[sl].size else np.nan
        print(f'5d s{r.seed_index}: W128 b/b0={bw / b0:.3f} c_eq={c_eq(hT, v0, b0, bw):.3f} Rres={rres(hT, s0, bw):.5f} | TW128 b/b0={btw / b0:.3f} c_eq={c_eq(hT, v0, b0, btw):.3f} Rres={rres(hT, s0, btw):.5f}  tissue share of window={T0[sl].mean():.2f}')
    print(f'total {time.time() - t:.1f}s')
