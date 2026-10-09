"""Plan v2 gate 5 (direction) on one ROI: correct synthesis, an H/E label swap, a wrong k column, and the 1/6 factor."""
import sys, time
from v2lib import *

fn = sys.argv[1]
rgb, hem, mpp = load(fn)
gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY); T0 = gray <= 215
M = macenko_basis(rgb, np.random.default_rng([20261008, FILES.index(fn)]))
pinvM = np.linalg.pinv(M)
print(f'{fn}: h@Hcol {M[:,0] @ HCOL:.3f}  e@Hcol {M[:,1] @ HCOL:.3f}')
flat0 = rgb.reshape(-1, 3)
C = od(flat0) @ pinvM.T   # concentrations from the CORRECT basis, reused for the regressor where the gate would


def gate(rgb1, hem1, kvec, cols, Mreg, Creg, label):
    """Regress dhem on (k-1)*max(c,0)*(s@Hcol)/6 over the gate's pixel filter; Mreg/Creg are what the harness believes."""
    f1 = rgb1.reshape(-1, 3)
    sel = T0.reshape(-1) & (flat0 >= 20).all(1) & (flat0 <= 250).all(1) & (f1 >= 20).all(1) & (f1 <= 250).all(1) & (hem.reshape(-1) > 0.01) & (hem1.reshape(-1) > 0.01)
    x = np.zeros(sel.sum())
    for k, col in zip(kvec, cols):
        x += (k - 1) * np.maximum(Creg[sel, col], 0) * (Mreg[:, col] @ HCOL) / 6
    y = (hem1.reshape(-1) - hem.reshape(-1))[sel].astype(np.float64)
    slope = float((x @ y) / (x @ x))            # through the origin
    A = np.c_[x, np.ones_like(x)]; s2, icpt = np.linalg.lstsq(A, y, rcond=None)[0]
    r = np.corrcoef(x, y)[0, 1]
    print(f'  {label:42s} n={sel.sum():>9d}  slope(origin) {slope:6.3f}  slope(OLS) {s2:6.3f}  icpt {icpt:+.5f}  r {r:.3f}  pass={0.9 <= slope <= 1.1}')


for name, kH, kE in (('H06', 0.6, 1.0), ('H14', 1.4, 1.0), ('E15', 1.0, 1.5)):
    t = time.time()
    rgb1 = synth_v2(rgb, M, kH, kE); hem1 = ch.to_channel(rgb1, 'hematoxylin_od')
    print(f'{name}: (synthesis + hem {time.time() - t:.1f}s)')
    gate(rgb1, hem1, (kH, kE), (0, 1), M, C, 'correct synthesis, correct regressor')
    # bug A: basis labels swapped at source (macenko returns [e, h]); harness consistently uses the swapped labels
    Msw = M[:, ::-1].copy(); Csw = C[:, ::-1]
    rgbA = synth_v2(rgb, Msw, kH, kE); hemA = ch.to_channel(rgbA, 'hematoxylin_od')
    gate(rgbA, hemA, (kH, kE), (0, 1), Msw, Csw, 'BUG A: H/E labels swapped (self-consistent)')
    # bug B: the kH factor applied to the E concentration (wrong column), regressor as specified
    rgbB = synth_v2(rgb, M, kE, kH); hemB = ch.to_channel(rgbB, 'hematoxylin_od')
    gate(rgbB, hemB, (kH, kE), (0, 1), M, C, 'BUG B: k columns swapped (kH on e, kE on h)')
    # bug C: regressor without the 1/6 factor -> slope should be ~1/6
    gate(rgb1, hem1, [1 + 6 * (k - 1) for k in (kH, kE)], (0, 1), M, C, 'regressor without 1/6 (expect ~0.167)')
    zero_share = float(((flat0 == 0).any(1)).mean()); zero_tis = float(((flat0 == 0).any(1) & T0.reshape(-1)).sum() / T0.sum())
print(f'zero-channel pixels: {zero_share:.5f} of ROI, {zero_tis:.5f} of tissue')
