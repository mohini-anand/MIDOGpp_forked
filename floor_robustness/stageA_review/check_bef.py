"""5b identity + orthogonality, 5e peak RSS, 5f timing of one perturbation (synthesis vs rgb2hed), on one ROI (argv[1])."""
import sys, time, resource
import numpy as np, cv2
from skimage.filters import threshold_multiotsu
from common import *


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9  # bytes on macOS


fn = sys.argv[1]
roi_index = FILES.index(fn)
T = {}
t = time.time(); rgb = ds.load_roi(path_of(fn)); T['load'] = time.time() - t
print(f'{fn} shape {rgb.shape}  after load peak RSS {rss_gb():.2f} GB')
t = time.time(); hem = ch.to_channel(rgb, 'hematoxylin_od'); T['hem0'] = time.time() - t
print(f'after hem0 peak RSS {rss_gb():.2f} GB')
gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY); T0 = gray <= 215
b0 = roi_floors(hem)['b_otsu']
t = time.time(); hT = hem[T0]; v0 = np.sort(hT); s0 = hT < np.float32(b0); T['sort_v0'] = time.time() - t
del hT
print(f'after v0 sort peak RSS {rss_gb():.2f} GB  (N_t={v0.size})')
t = time.time(); M = macenko_basis(rgb, np.random.default_rng([20261008, roi_index])); T['macenko'] = time.time() - t
print('M (cols h, e):\n', np.round(M, 4), ' angle h-e deg', np.degrees(np.arccos(M[:, 0] @ M[:, 1])), ' cond', np.linalg.cond(M))

# --- 5b identity on the full image
t = time.time(); rgb_id = synth(rgb, M, [1, 1]); T['synth_ID'] = time.time() - t
print('5b ID bit-identical:', np.array_equal(rgb_id, rgb), ' n differing px:', int((rgb_id != rgb).any(-1).sum()))
del rgb_id
# identity with a deliberately WRONG basis (random unit vectors) -- does the gate still pass?
rng = np.random.default_rng(1)
Mbad = rng.normal(size=(3, 2)); Mbad /= np.linalg.norm(Mbad, axis=0)
print('5b ID with a random wrong basis bit-identical:', np.array_equal(synth(rgb[:512], Mbad, [1, 1]), rgb[:512]))
# orthogonality: residual vs stain plane on a 100k sample; the plan writes `M.T @ R` (R is N x 3 -> shape error); correct is R @ M
flat = rgb.reshape(-1, 3)
idx = np.random.default_rng(2).choice(len(flat), 100_000, replace=False)
O = od(flat[idx]); pinvM = np.linalg.pinv(M); C = O @ pinvM.T; R = O - C @ M.T
try:
    _ = M.T @ R
    print('M.T @ R ran (unexpected)')
except ValueError as e:
    print('plan literal `M.T @ R` raises:', e)
print(f'5b orthogonality max|R @ M| = {np.abs(R @ M).max():.3e}   max|R| = {np.abs(R).max():.3e}  max OD {O.max():.3f}')
print(f'   with the wrong basis: max|R @ Mbad| = {np.abs((O - (O @ np.linalg.pinv(Mbad).T) @ Mbad.T) @ Mbad).max():.3e}')
tis = T0.reshape(-1)[idx]
print(f'   residual share ||R||/||OD|| on tissue sample = {np.linalg.norm(R[tis]) / np.linalg.norm(O[tis]):.4f}')
print(f'   negative concentrations on tissue: cH<0 {np.mean(C[tis,0] < 0):.3f}, cE<0 {np.mean(C[tis,1] < 0):.3f}; p99 cH {np.percentile(C[tis,0],99):.3f} cE {np.percentile(C[tis,1],99):.3f}')

# --- 5f one perturbation end to end: H06
t = time.time(); rgb1 = synth(rgb, M, [0.6, 1]); T['synth_H06'] = time.time() - t
print(f'after synthesis peak RSS {rss_gb():.2f} GB')
t = time.time(); hem1 = ch.to_channel(rgb1, 'hematoxylin_od'); T['rgb2hed_H06'] = time.time() - t
print(f'after rgb2hed(perturbed) peak RSS {rss_gb():.2f} GB')
t = time.time(); gray1 = cv2.cvtColor(rgb1, cv2.COLOR_RGB2GRAY); ht = hem1[T0]; T['gather'] = time.time() - t
t = time.time(); b1 = float(np.float32(threshold_multiotsu(hem1, classes=3)[0])); T['otsu_roi'] = time.time() - t
t = time.time(); bt = float(np.float32(threshold_multiotsu(hem1[gray1 <= 215], classes=3)[0])); T['otsu_tissue'] = time.time() - t
t = time.time(); ce = c_eq(ht, v0, b0, b1); T['c_eq'] = time.time() - t
t = time.time(); rr = rres(ht, s0, b1); T['rres'] = time.time() - t
print(f'H06: b_otsu {b1:.5f} (b/b0 {b1 / b0:.3f}) c_eq {ce:.3f} Rres {rr:.5f}; TOTSU c_eq {c_eq(ht, v0, b0, bt):.3f}')
t = time.time(); _ = cv2.imdecode(cv2.imencode('.jpg', cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, 75, cv2.IMWRITE_JPEG_SAMPLING_FACTOR, cv2.IMWRITE_JPEG_SAMPLING_FACTOR_420])[1], cv2.IMREAD_COLOR); T['jpeg'] = time.time() - t
print('timings (s):', {k: round(v, 2) for k, v in T.items()})
print(f'final peak RSS {rss_gb():.2f} GB')
