"""Time one synthesis + rgb2hed on 201.tiff; thread env decided by the caller (set before numpy import here)."""
import os, sys
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1' if sys.argv[1] == '1' else '18'
import time
from common import *   # common sets env with setdefault before numpy; '0' path pre-imports numpy below
import numpy as np, cv2
if sys.argv[1] == '1':
    cv2.setNumThreads(1)
from threadpoolctl import threadpool_info
print('threads:', [(i['internal_api'], i['num_threads']) for i in threadpool_info()], 'cv2', cv2.getNumThreads())
rgb = ds.load_roi(path_of('201.tiff'))
M = macenko_basis(rgb, np.random.default_rng([20261008, FILES.index('201.tiff')]))
t = time.time(); rgb1 = synth(rgb, M, [0.6, 1]); ts = time.time() - t
t = time.time(); hem1 = ch.to_channel(rgb1, 'hematoxylin_od'); th = time.time() - t
t = time.time(); b = cv2.GaussianBlur(rgb, (0, 0), 1.0); tb = time.time() - t
print(f'synth {ts:.2f}s  rgb2hed {th:.2f}s  blur {tb:.3f}s')
