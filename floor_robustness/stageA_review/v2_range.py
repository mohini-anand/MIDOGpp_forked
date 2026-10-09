"""Histogram range of b_otsu under synthesis: zero-channel pixels before/after, hem max, ROI bin width relative to b'."""
import sys
from v2lib import *
from skimage.filters import threshold_multiotsu
for fn in sys.argv[1:]:
    rgb, hem, mpp = load(fn)
    M = macenko_basis(rgb, np.random.default_rng([20261008, FILES.index(fn)]))
    V = variants(rgb, M)
    z0 = (rgb == 0).any(-1)
    b0 = float(np.float32(threshold_multiotsu(hem, classes=3)[0]))
    print(f'{fn}: original zero-channel px {int(z0.sum())}, hem max {hem.max():.3f}, bin/b {(hem.max()-hem.min())/256/b0:.3f}')
    for name in ('H06', 'H14', 'E15', 'GRAD', 'CUP', 'GAIN09'):
        r1 = V[name](); h1 = ch.to_channel(r1, 'hematoxylin_od'); z1 = (r1 == 0).any(-1)
        b1 = float(np.float32(threshold_multiotsu(h1, classes=3)[0]))
        print(f'  {name:6s} zero-channel px {int(z1.sum()):8d} (new {int((z1 & ~z0).sum()):8d})  hem max {h1.max():.3f}  bin/b\' {(h1.max()-h1.min())/256/b1:.3f}')
