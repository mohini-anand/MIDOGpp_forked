"""Is production hem on 1024-row blocks bit-identical to the whole-image call? Original + H06 + JPG + BLUR on one ROI; also odd block sizes. Peak RSS of the blocked path in a fresh process (argv[2] == 'rss')."""
import sys, time, resource
from v2lib import *


def hem_blocked(rgb, rows=1024):
    out = np.empty(rgb.shape[:2], np.float32)
    for r0 in range(0, rgb.shape[0], rows):
        out[r0:r0 + rows] = ch.to_channel(rgb[r0:r0 + rows], 'hematoxylin_od')
    return out


fn = sys.argv[1]
rgb = ds.load_roi(path_of(fn))
if len(sys.argv) > 2 and sys.argv[2] == 'rss':
    t = time.time(); h = hem_blocked(rgb); print(f'{fn} blocked hem {time.time() - t:.2f}s; peak RSS after load + blocked hem {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9:.2f} GB')
    sys.exit()
M = macenko_basis(rgb, np.random.default_rng([20261008, FILES.index(fn)]))
cases = {'original': rgb, 'H06': synth_v2(rgb, M, 0.6, 1.0), 'JPG': jpeg(rgb, 75), 'BLUR': cv2.GaussianBlur(rgb, (0, 0), 1.0)}
for name, im in cases.items():
    t = time.time(); whole = ch.to_channel(im, 'hematoxylin_od'); tw = time.time() - t
    for rows in (1024, 777, 1):
        if rows == 1 and name != 'original':
            continue
        t = time.time(); blk = hem_blocked(im, rows); tb = time.time() - t
        nd = int(np.count_nonzero(blk != whole))
        print(f'{fn} {name:8s} rows={rows:5d}: bit-identical={nd == 0} (differing px {nd}, max|d| {float(np.abs(blk - whole).max()):.2e})  whole {tw:.2f}s blocked {tb:.2f}s')
