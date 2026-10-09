"""
    Stage A of the floor-clip stress test: the per-ROI worker for tm_floor_stress_stageA_49roi.ipynb.

    One task is one ROI. For every image variant of TM_FLOOR_STRESS_STAGE_A_PLAN.md section 3.2, the worker computes:
    - the four floor estimators of section 3.3;
    - each floor's equivalent level c_eq on the original ladder (section 3.4);
    - the template-box refinement, re-run on the variant;
    - the gate quantities of section 4.4.
    It writes nothing. The notebook writes every output.

    Spawn workers on macOS cannot unpickle functions defined in a notebook cell, so this code lives in a module.
    The BLAS/OpenMP thread limits are set here, before numpy is imported.
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'

import resource
import sys
import time
from pathlib import Path

import cv2
import numpy as np
from skimage.color import hed_from_rgb
from skimage.filters import threshold_multiotsu
from skimage.measure import label, regionprops
from threadpoolctl import threadpool_info, threadpool_limits

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from midog_utils import channels as ch  # noqa: E402
from midog_utils import dataset as ds  # noqa: E402
from midog_utils import seed_selection as ss  # noqa: E402
from midog_utils import template_match as tm  # noqa: E402

# the floor estimators (plan section 3.3)
FLOOR_CLASSES = 3              # lower threshold of a 3-class multi-Otsu, nbins 256 (skimage's default)
WHITE_GRAY = 215               # glass: gray > 215; tissue: gray <= 215
WINDOW_UM = 128.0              # W128 / TW128 window side
SHIFT_UM = 64.0                # SHIFT offset
TW_MIN_TISSUE = 0.05           # TW128 is NaN below this tissue share of its window
# the synthesis (plan section 3.2)
MACENKO_SEED = 20261008
MACENKO_N = 500_000
MACENKO_BETA = 0.15            # log10 OD, as floor_robustness/macenko_vectors.py
MACENKO_ALPHA = 1.0            # percent
SYNTH_BLOCK_ROWS = 256
HEM_BLOCK_ROWS = 1024          # production hematoxylin_od on row blocks (bit-identical, gated in the smoke run)
GRAD_RANGE = (0.8, 1.2)        # section-thickness gradient along the long axis, pale at index 0
# gates (plan section 4.4)
DIRECTION_SAMPLE = 2_000_000   # pixels on a fixed stride grid for the direction gate and the residual share
DIRECTION_RGB = (20, 250)      # every channel inside this range, before and after
DIRECTION_MIN_HEM = 0.01       # hem0 and hem' above this
LADDER_CHECK = (0.25, 0.5, 0.75, 0.9, 1.0, 1.5, 2.0)
HCOL = hed_from_rgb[:, 0]      # rgb2hed's H column
SHIFT_DIRS = (('E', 1, 0), ('W', -1, 0), ('S', 0, 1), ('N', 0, -1), ('SE', 1, 1), ('SW', -1, 1), ('NE', 1, -1), ('NW', -1, -1))
STAIN_K = {'H06': (0.6, 1.0), 'H14': (1.4, 1.0), 'E15': (1.0, 1.5)}
VARIANT_CLASS = {'UNPERT': 'reference', 'ID': 'gate', 'NULL2': 'gate', 'GAMMA': 'calibration', 'JPG95': 'calibration', 'H06': 'realistic', 'H14': 'realistic', 'E15': 'realistic', 'CDN': 'realistic', 'CUP': 'realistic', 'GAIN09': 'realistic', 'GAIN11': 'realistic', 'WBR095': 'realistic', 'WBR105': 'realistic', 'WBB095': 'realistic', 'WBB105': 'realistic', 'BLUR': 'realistic', 'JPG': 'realistic', 'GRAD': 'realistic', 'SHIFT': 'realistic'}
IMAGE_VARIANTS = ('GAMMA', 'JPG95', 'H06', 'H14', 'E15', 'CDN', 'CUP', 'GAIN09', 'GAIN11', 'WBR095', 'WBR105', 'WBB095', 'WBB105', 'BLUR', 'JPG', 'GRAD')

# ---------------------------------------------------------------------------------------------------------------------
# verbatim copies (gate 2): refine_template_box and its constants from tm_border_floor_49roi_6seed.ipynb cells 1 and 4,
# roi_floors from its cell 8, window_floor from tm_floor_definition_49roi_6seed.ipynb cell 9
# ---------------------------------------------------------------------------------------------------------------------
# refine_template_box's constants, copied verbatim from cell 1 of multiotsu_template_refinement_14roi.ipynb at c767898
WINDOW = tm.BASE_SIZE                                            # 51 px window around the click
CLASSES = 3                                                      # multi-Otsu classes; the top one is kept
DILATE_KERNEL = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))  # 3x3 cross
DILATE_ITERATIONS = 1                                            # 1 px
MIN_PIECE_UM2 = 1.0                                              # speck floor
MIN_BASE_UM = 5.0                                                # base case: the click's piece bbox, longer side
CLICK_RADIUS_UM = 3.75                                           # fallback: piece centre within this of the click
MIN_AREA_UM2 = 3.0                                               # selection must cover at least this
MAX_AREA_FRAC = 0.85                                             # ... and at most this fraction of the window
PAD_PX = 2                                                       # margin kept around the selection on every side
assert PAD_PX == 2 and 'PAD_UM' not in globals(), 'gate 1: the refinement margin must be the fixed 2 px of c767898'
REFINE_SHA1 = '21c43a20db00dccf7ac3bab88aee423c345441b6'         # inspect.getsource(refine_template_box) at c767898


def refine_template_box(hem_img, cx, cy, mpp):
    """
        Find the mitotic figure's chromatin around a click and return the template box.
        hem_img (np.ndarray): the full-ROI hematoxylin_od channel.
        cx, cy (float): click position in ROI pixel coordinates.
        mpp (float): this ROI's microns per pixel.
        Returns dict: accepted (bool), reason (str), base_size (int, template side in px), center_xy (ROI x, y of the template centre, on the selection's bounding box), plus intermediates for plotting.
    """
    hem_patch = tm.read_padded_patch(hem_img, cx, cy, WINDOW)
    result = dict(accepted=False, reason='too close to the ROI border to read the window', base_size=None, center_xy=None,
                  hem_patch=hem_patch, thresholds=None, class_map=None, labels=None, pieces=[], mode=None,
                  selected=[], selection_mask=None, bbox=None, radius_px=CLICK_RADIUS_UM / mpp, pad_px=PAD_PX)
    if hem_patch is None:
        return result

    patch = hem_patch.astype(np.float64)
    py, px = patch.shape[0] // 2, patch.shape[1] // 2
    thresholds = threshold_multiotsu(patch, classes=CLASSES)
    top_class = cv2.dilate((patch > thresholds[-1]).astype(np.uint8), DILATE_KERNEL, iterations=DILATE_ITERATIONS) > 0
    labels = label(top_class, connectivity=1)
    result.update(thresholds=[float(t) for t in thresholds], class_map=np.digitize(patch, thresholds), labels=labels)

    pieces = []
    for r in regionprops(labels):
        ys, xs = r.coords[:, 0], r.coords[:, 1]
        y0, x0, y1, x1 = r.bbox
        pieces.append(dict(label=int(r.label), area_um2=r.area * mpp ** 2, centre_um=float(np.hypot(ys.mean() - py, xs.mean() - px)) * mpp,
                           bbox_um=max(y1 - y0, x1 - x0) * mpp, centre_yx=(float(ys.mean()), float(xs.mean())), speck=r.area * mpp ** 2 < MIN_PIECE_UM2))
    result['pieces'] = pieces
    by_label = {p['label']: p for p in pieces}

    clicked = int(labels[py, px])
    if clicked and not by_label[clicked]['speck'] and by_label[clicked]['bbox_um'] >= MIN_BASE_UM:
        result['mode'] = 'base'
        selected = [clicked]
    else:
        result['mode'] = 'fallback'
        selected = [p['label'] for p in pieces if not p['speck'] and p['centre_um'] <= CLICK_RADIUS_UM]
    for p in pieces:
        p['selected'] = p['label'] in selected
    result['selected'] = selected
    if not selected:
        result['reason'] = f'no piece centre within {CLICK_RADIUS_UM}um of the click'
        return result

    mask = np.isin(labels, selected)
    result['selection_mask'] = mask
    if mask.sum() * mpp ** 2 < MIN_AREA_UM2:
        result['reason'] = f'selection {mask.sum() * mpp ** 2:.2f}um^2 < {MIN_AREA_UM2}um^2'
        return result
    if mask.sum() > MAX_AREA_FRAC * patch.size:
        result['reason'] = f'selection covers more than {MAX_AREA_FRAC:.0%} of the window'
        return result

    ys, xs = np.where(mask)
    y0, y1, x0, x1 = int(ys.min()), int(ys.max()) + 1, int(xs.min()), int(xs.max()) + 1
    half = WINDOW // 2
    ix, iy = int(round(cx)), int(round(cy))
    result.update(accepted=True, reason='accepted', bbox=(y0, y1, x0, x1), base_size=ss._odd(max(y1 - y0, x1 - x0), minimum=5) + 2 * result['pad_px'],
                  center_xy=(ix - half + (x0 + x1 - 1) / 2.0, iy - half + (y0 + y1 - 1) / 2.0))
    return result


def roi_floors(hem):
    """
        The two per-ROI floors over the whole unclipped channel: the lower 3-class multi-Otsu threshold and the median.

        hem (np.ndarray): this ROI's hematoxylin_od channel.

        Returns dict: b_otsu, b_otsu_upper (the second threshold, for reference) and b_med, as float32-representable floats.
    """
    th = threshold_multiotsu(hem, classes=FLOOR_CLASSES)
    return dict(b_otsu=float(np.float32(th[0])), b_otsu_upper=float(np.float32(th[1])), b_med=float(np.float32(np.median(hem))))


def window_floor(hem, cx, cy, side):
    """
        The lower 3-class multi-Otsu threshold over a square window centred on a click, truncated at the ROI edge.

        hem (np.ndarray): the unclipped channel.
        cx (float): click x.
        cy (float): click y.
        side (int): the odd window side in px.

        Returns tuple[float, int]: the float32-representable floor and the window's pixel count after truncation.
    """
    h, w = hem.shape
    half = side // 2
    ix, iy = int(round(cx)), int(round(cy))
    win = hem[max(0, iy - half):min(h, iy + half + 1), max(0, ix - half):min(w, ix + half + 1)]
    return float(np.float32(threshold_multiotsu(win, classes=FLOOR_CLASSES)[0])), int(win.size)


# ---------------------------------------------------------------------------------------------------------------------
# threads and memory
# ---------------------------------------------------------------------------------------------------------------------
_LIMITER = None


def init_worker():
    """
        Pool initialiser: pin every BLAS/OpenMP pool and OpenCV to one thread, and check that they are pinned.

        Returns None.
    """
    global _LIMITER
    _LIMITER = threadpool_limits(limits=1)
    cv2.setNumThreads(1)
    bad = [(m.get('internal_api'), m.get('num_threads')) for m in threadpool_info() if m.get('num_threads') != 1]
    assert not bad, f'gate 1: thread pools not pinned to 1: {bad}'


def thread_report():
    """
        The thread count of every loaded BLAS/OpenMP pool in this process.

        Returns str: 'api:threads' pairs joined by ';'.
    """
    return ';'.join(f"{m.get('internal_api')}:{m.get('num_threads')}" for m in threadpool_info())


def peak_rss_gb():
    """
        This process's peak resident set size so far.

        Returns float: GB (macOS reports ru_maxrss in bytes).
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9


# ---------------------------------------------------------------------------------------------------------------------
# channel, floors and the equivalent level
# ---------------------------------------------------------------------------------------------------------------------
def hem_of(rgb):
    """
        The production hematoxylin_od channel, computed on row blocks to bound memory.

        rgb (np.ndarray): H x W x 3 uint8 image.

        Returns np.ndarray: H x W float32, identical to ch.to_channel(rgb, 'hematoxylin_od') (gated in the smoke run).
    """
    out = np.empty(rgb.shape[:2], np.float32)
    for r0 in range(0, rgb.shape[0], HEM_BLOCK_ROWS):
        out[r0:r0 + HEM_BLOCK_ROWS] = ch.to_channel(rgb[r0:r0 + HEM_BLOCK_ROWS], 'hematoxylin_od')
    return out


def b_otsu_of(hem):
    """
        The ROI-wide floor: the same computation as roi_floors' b_otsu, without its median.

        hem (np.ndarray): a hematoxylin_od channel.

        Returns float: the float32-representable lower 3-class multi-Otsu threshold.
    """
    return float(np.float32(threshold_multiotsu(hem, classes=FLOOR_CLASSES)[0]))


def bin_width(values):
    """
        The width of one multi-Otsu histogram bin over these values (skimage bins the full range into 256).

        values (np.ndarray): the pixels the threshold was computed on.

        Returns float: (max - min) / 256.
    """
    return float((float(values.max()) - float(values.min())) / 256.0)


def window_slices(shape, cx, cy, side):
    """
        The click-centred window that window_floor reads, truncated at the ROI edge.

        shape (tuple): (h, w) of the channel.
        cx (float): window centre x.
        cy (float): window centre y.
        side (int): the odd window side in px.

        Returns tuple[slice, slice]: row and column slices.
    """
    h, w = shape
    half = side // 2
    ix, iy = int(round(cx)), int(round(cy))
    return slice(max(0, iy - half), min(h, iy + half + 1)), slice(max(0, ix - half), min(w, ix + half + 1))


def tissue_window_floor(hem, gray, cx, cy, side):
    """
        TW128: the lower multi-Otsu threshold over the tissue pixels of the click-centred window.

        hem (np.ndarray): the channel.
        gray (np.ndarray): the same image's gray (uint8); tissue is gray <= WHITE_GRAY.
        cx (float): click x.
        cy (float): click y.
        side (int): the odd window side in px.

        Returns tuple[float, float, float]: the floor (NaN below TW_MIN_TISSUE tissue), the window's tissue share, and the bin width.
    """
    sl = window_slices(hem.shape, cx, cy, side)
    tis = gray[sl] <= WHITE_GRAY
    share = float(tis.mean())
    if share < TW_MIN_TISSUE:
        return float('nan'), share, float('nan')
    vals = hem[sl][tis]
    return float(np.float32(threshold_multiotsu(vals, classes=FLOOR_CLASSES)[0])), share, bin_width(vals)


def c_eq(ht, v0, b0, b):
    """
        The equivalent level: the multiple of the original b_otsu that clips the same share of original tissue pixels.

        ht (np.ndarray): the variant's channel at the original tissue pixels (float32).
        v0 (np.ndarray): the original channel's tissue pixels, sorted (float32).
        b0 (float): the original ROI b_otsu.
        b (float): the floor on the variant.

        Returns float: c_eq, or NaN when b is NaN.
    """
    if not np.isfinite(b):
        return float('nan')
    p = np.count_nonzero(ht < np.float32(b)) / ht.size
    return float(v0[min(int(round(p * ht.size)), ht.size - 1)]) / b0


# ---------------------------------------------------------------------------------------------------------------------
# synthesis (plan section 3.2)
# ---------------------------------------------------------------------------------------------------------------------
def od10(rgb):
    """
        The synthesis model's optical density, -log10((I + 1) / 256), per channel.

        rgb (np.ndarray): uint8 pixels, any shape ending in 3.

        Returns np.ndarray: float64 OD of the same shape.
    """
    return -np.log10((rgb.astype(np.float64) + 1.0) / 256.0)


def rgb_from_od10(od):
    """
        Back from the synthesis OD to uint8 pixels.

        od (np.ndarray): float64 OD.

        Returns np.ndarray: uint8, clip(rint(256 * 10^-od - 1), 0, 255).
    """
    return np.clip(np.rint(256.0 * 10.0 ** (-od) - 1.0), 0, 255).astype(np.uint8)


def macenko_basis(rgb, rng):
    """
        The slide's own H and E stain vectors, Macenko 2009 Algorithm 1, as floor_robustness/macenko_vectors.py computes them.

        rgb (np.ndarray): H x W x 3 uint8 image.
        rng (np.random.Generator): the per-ROI generator for the pixel sample.

        Returns np.ndarray: 3 x 2 matrix [h e] of unit OD vectors (log10 units), h the one with the larger red component.
    """
    flat = rgb.reshape(-1, 3)
    px = flat[rng.choice(len(flat), MACENKO_N, replace=False)].astype(np.float64)
    od = -np.log10((px + 1) / 256.0)
    od = od[(od > MACENKO_BETA).all(1)]
    _, _, vt = np.linalg.svd(od, full_matrices=False)
    v = vt[:2].T
    if v[:, 0].sum() < 0:
        v[:, 0] *= -1
    if v[:, 1].sum() < 0:
        v[:, 1] *= -1
    proj = od @ v
    phi = np.arctan2(proj[:, 1], proj[:, 0])
    lo, hi = np.percentile(phi, [MACENKO_ALPHA, 100 - MACENKO_ALPHA])
    s1 = v @ np.array([np.cos(lo), np.sin(lo)])
    s2 = v @ np.array([np.cos(hi), np.sin(hi)])
    s1, s2 = s1 / np.linalg.norm(s1), s2 / np.linalg.norm(s2)
    h, e = (s1, s2) if s1[0] > s2[0] else (s2, s1)
    return np.stack([h, e], 1)


def synth_stain(rgb, M, kH, kE):
    """
        Rescale the positive part of the slide's H and E concentrations; keep the residual; copy pixels with any channel at 0.

        rgb (np.ndarray): H x W x 3 uint8 image.
        M (np.ndarray): 3 x 2 stain matrix [h e].
        kH (float): hematoxylin factor.
        kE (float): eosin factor.

        Returns np.ndarray: the perturbed uint8 image, OD' = OD + (kH - 1) max(cH, 0) h + (kE - 1) max(cE, 0) e.
    """
    pinv = np.linalg.pinv(M)
    h, e = M[:, 0], M[:, 1]
    out = np.empty_like(rgb)
    for r0 in range(0, rgb.shape[0], SYNTH_BLOCK_ROWS):
        blk = rgb[r0:r0 + SYNTH_BLOCK_ROWS]
        flat = blk.reshape(-1, 3)
        od = od10(flat)
        cpos = np.maximum(od @ pinv.T, 0)
        new = rgb_from_od10(od + (kH - 1.0) * cpos[:, [0]] * h[None, :] + (kE - 1.0) * cpos[:, [1]] * e[None, :])
        zero = (flat == 0).any(1)
        new[zero] = flat[zero]
        out[r0:r0 + SYNTH_BLOCK_ROWS] = new.reshape(blk.shape)
    return out


def synth_identity(rgb, M):
    """
        The ID gate: reconstruct OD fully as (C [1, 1]) M^T + R and re-encode; tests the encoding path only.

        rgb (np.ndarray): H x W x 3 uint8 image.
        M (np.ndarray): 3 x 2 stain matrix.

        Returns np.ndarray: the re-encoded uint8 image (must equal rgb bit for bit).
    """
    pinv = np.linalg.pinv(M)
    out = np.empty_like(rgb)
    for r0 in range(0, rgb.shape[0], SYNTH_BLOCK_ROWS):
        blk = rgb[r0:r0 + SYNTH_BLOCK_ROWS]
        od = od10(blk.reshape(-1, 3))
        C = od @ pinv.T
        R = od - C @ M.T
        out[r0:r0 + SYNTH_BLOCK_ROWS] = rgb_from_od10((C * np.array([1.0, 1.0])) @ M.T + R).reshape(blk.shape)
    return out


def synth_gradient(rgb):
    """
        Section-thickness gradient: OD scaled by g, linear from GRAD_RANGE[0] to GRAD_RANGE[1] along the long axis.

        rgb (np.ndarray): H x W x 3 uint8 image.

        Returns np.ndarray: the perturbed uint8 image; pixels with any channel at 0 are copied unchanged.
    """
    H, W = rgb.shape[:2]
    g0, g1 = GRAD_RANGE
    out = np.empty_like(rgb)
    for r0 in range(0, H, SYNTH_BLOCK_ROWS):
        blk = rgb[r0:r0 + SYNTH_BLOCK_ROWS]
        od = od10(blk)
        if W >= H:
            g = (g0 + (g1 - g0) * np.arange(W) / (W - 1))[None, :, None]
        else:
            g = (g0 + (g1 - g0) * np.arange(r0, r0 + blk.shape[0]) / (H - 1))[:, None, None]
        new = rgb_from_od10(od * g)
        zero = (blk == 0).any(-1)
        new[zero] = blk[zero]
        out[r0:r0 + SYNTH_BLOCK_ROWS] = new
    return out


def lut_all(rgb, f):
    """
        Apply one tone curve to every channel.

        rgb (np.ndarray): uint8 image.
        f (callable): maps float64 levels 0..255 to new levels.

        Returns np.ndarray: uint8 image through clip(rint(f(I)), 0, 255).
    """
    t = np.clip(np.rint(f(np.arange(256, dtype=np.float64))), 0, 255).astype(np.uint8)
    return cv2.LUT(rgb, t)


def lut_one(rgb, channel, f):
    """
        Apply a tone curve to one channel only.

        rgb (np.ndarray): uint8 image, RGB order.
        channel (int): 0 red, 1 green, 2 blue.
        f (callable): maps float64 levels 0..255 to new levels.

        Returns np.ndarray: a copy with that channel through clip(rint(f(I)), 0, 255).
    """
    t = np.clip(np.rint(f(np.arange(256, dtype=np.float64))), 0, 255).astype(np.uint8)
    out = rgb.copy()
    out[..., channel] = cv2.LUT(np.ascontiguousarray(rgb[..., channel]), t)
    return out


def jpeg(rgb, quality):
    """
        JPEG round trip at 4:2:0 through OpenCV.

        rgb (np.ndarray): uint8 image, RGB order.
        quality (int): JPEG quality.

        Returns np.ndarray: the decoded uint8 RGB image.
    """
    ok, buf = cv2.imencode('.jpg', cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, quality, cv2.IMWRITE_JPEG_SAMPLING_FACTOR, cv2.IMWRITE_JPEG_SAMPLING_FACTOR_420])
    assert ok, 'JPEG encode failed'
    return cv2.cvtColor(cv2.imdecode(buf, cv2.IMREAD_COLOR), cv2.COLOR_BGR2RGB)


def make_rgb_variant(var, rgb, M):
    """
        Build one RGB-level image variant (plan section 3.2).

        var (str): a name in IMAGE_VARIANTS.
        rgb (np.ndarray): the original uint8 image.
        M (np.ndarray): the ROI's 3 x 2 stain matrix.

        Returns np.ndarray: the perturbed uint8 image.
    """
    if var in STAIN_K:
        return synth_stain(rgb, M, *STAIN_K[var])
    makers = {
        'GAMMA': lambda: lut_all(rgb, lambda I: 255.0 * (I / 255.0) ** 1.1),
        'JPG95': lambda: jpeg(rgb, 95),
        'CDN': lambda: lut_all(rgb, lambda I: 12.0 + I * 231.0 / 255.0),
        'CUP': lambda: lut_all(rgb, lambda I: (I - 12.0) * 255.0 / 231.0),
        'GAIN09': lambda: lut_all(rgb, lambda I: 0.9 * I),
        'GAIN11': lambda: lut_all(rgb, lambda I: 1.1 * I),
        'WBR095': lambda: lut_one(rgb, 0, lambda I: 0.95 * I),
        'WBR105': lambda: lut_one(rgb, 0, lambda I: 1.05 * I),
        'WBB095': lambda: lut_one(rgb, 2, lambda I: 0.95 * I),
        'WBB105': lambda: lut_one(rgb, 2, lambda I: 1.05 * I),
        'BLUR': lambda: cv2.GaussianBlur(rgb, (0, 0), 1.0, borderType=cv2.BORDER_REFLECT_101),
        'JPG': lambda: jpeg(rgb, 75),
        'GRAD': lambda: synth_gradient(rgb),
    }
    return makers[var]()


# ---------------------------------------------------------------------------------------------------------------------
# the per-ROI task
# ---------------------------------------------------------------------------------------------------------------------
def shift_window_centre(gray0, cx, cy, side, d):
    """
        SHIFT: of the eight offsets of d px, the window centre whose truncated window holds the most glass.

        gray0 (np.ndarray): the original gray (uint8).
        cx (float): click x.
        cy (float): click y.
        side (int): the W128 window side in px.
        d (int): the offset in px.

        Returns tuple[float, float, str, float, float]: centre x, centre y, direction, glass share, and the window's area share of side^2.
    """
    H, W = gray0.shape
    best = None
    for name, sx, sy in SHIFT_DIRS:
        x = min(max(cx + sx * d, 0.0), W - 1.0)
        y = min(max(cy + sy * d, 0.0), H - 1.0)
        sl = window_slices(gray0.shape, x, y, side)
        win = gray0[sl]
        glass = float((win > WHITE_GRAY).mean())
        if best is None or glass > best[3]:
            best = (x, y, name, glass, win.size / side ** 2)
    return best


def estimator_rows(var, hem1, gray1, ht, v0, b0, clicks, side):
    """
        Every floor estimator on one variant, mapped to c_eq.

        var (str): the variant name.
        hem1 (np.ndarray): the variant's channel.
        gray1 (np.ndarray): the gray the variant's tissue-only estimators use.
        ht (np.ndarray): hem1 at the original tissue pixels.
        v0 (np.ndarray): sorted original tissue values.
        b0 (float): the original b_otsu.
        clicks (list): the ROI's click dicts.
        side (int): the W128 window side in px.

        Returns list[dict]: one row per (estimator, click); ROI-wide estimators carry seed_index -1.
    """
    rows = []
    b = b_otsu_of(hem1)
    rows.append(dict(var=var, est='b_otsu', seed_index=-1, b=b, c_eq=c_eq(ht, v0, b0, b), bin_width=bin_width(hem1), win_px=int(hem1.size), tissue_share=np.nan))
    tis = gray1 <= WHITE_GRAY
    vals = hem1[tis]
    b = float(np.float32(threshold_multiotsu(vals, classes=FLOOR_CLASSES)[0]))
    rows.append(dict(var=var, est='TOTSU', seed_index=-1, b=b, c_eq=c_eq(ht, v0, b0, b), bin_width=bin_width(vals), win_px=int(vals.size), tissue_share=float(tis.mean())))
    del tis, vals
    for r in clicks:
        b, n_px = window_floor(hem1, r['click_cx'], r['click_cy'], side)
        sl = window_slices(hem1.shape, r['click_cx'], r['click_cy'], side)
        rows.append(dict(var=var, est='W128', seed_index=r['seed_index'], b=b, c_eq=c_eq(ht, v0, b0, b), bin_width=bin_width(hem1[sl]), win_px=n_px, tissue_share=np.nan))
        b, share, bw = tissue_window_floor(hem1, gray1, r['click_cx'], r['click_cy'], side)
        rows.append(dict(var=var, est='TW128', seed_index=r['seed_index'], b=b, c_eq=c_eq(ht, v0, b0, b), bin_width=bw, win_px=n_px, tissue_share=share))
    return rows


def box_rows(var, hem1, clicks, mpp):
    """
        Re-run the template-box refinement on one variant and compare with the recorded boxes.

        var (str): the variant name.
        hem1 (np.ndarray): the variant's channel.
        clicks (list): the ROI's click dicts, with the recorded tcx, tcy, bs.
        mpp (float): microns per pixel.

        Returns list[dict]: one row per click: accepted, reason, bs, tcx, tcy, and whether the box equals the recorded one.
    """
    rows = []
    for r in clicks:
        res = refine_template_box(hem1, r['click_cx'], r['click_cy'], mpp)
        row = dict(var=var, seed_index=r['seed_index'], accepted=bool(res['accepted']), reason=res['reason'], bs=-1, tcx=-1, tcy=-1)
        if res['accepted']:
            row.update(bs=int(res['base_size']), tcx=int(round(res['center_xy'][0])), tcy=int(round(res['center_xy'][1])))
        row['box_matches'] = bool(row['accepted'] and (row['bs'], row['tcx'], row['tcy']) == (r['bs'], r['tcx'], r['tcy']))
        rows.append(row)
    return rows


def direction_slope(sample, rgb1, hem1, M, kH, kE):
    """
        Gate 5: the slope through the origin of the change in hem against the change the synthesis should cause.

        sample (dict): the stride-grid sample of the original (idx, rgb, hem, tissue, C).
        rgb1 (np.ndarray): the stain-perturbed image.
        hem1 (np.ndarray): its channel.
        M (np.ndarray): the 3 x 2 stain matrix the synthesis used.
        kH (float): hematoxylin factor.
        kE (float): eosin factor.

        Returns tuple[float, float, int]: the slope, the correlation, and the number of pixels used.
    """
    lo, hi = DIRECTION_RGB
    rgb1s = rgb1.reshape(-1, 3)[sample['idx']]
    hem1s = hem1.reshape(-1)[sample['idx']]
    rgb0s, hem0s = sample['rgb'], sample['hem']
    sel = sample['tissue'] & (rgb0s >= lo).all(1) & (rgb0s <= hi).all(1) & (rgb1s >= lo).all(1) & (rgb1s <= hi).all(1) & (hem0s > DIRECTION_MIN_HEM) & (hem1s > DIRECTION_MIN_HEM)
    C = np.maximum(sample['C'][sel], 0)
    x = (kH - 1.0) * C[:, 0] * float(M[:, 0] @ HCOL) / 6.0 + (kE - 1.0) * C[:, 1] * float(M[:, 1] @ HCOL) / 6.0
    y = hem1s[sel].astype(np.float64) - hem0s[sel].astype(np.float64)
    slope = float(x @ y / (x @ x)) if x.size else float('nan')
    r = float(np.corrcoef(x, y)[0, 1]) if x.size > 2 else float('nan')
    return slope, r, int(sel.sum())


def run_roi(task):
    """
        Stage A for one ROI: every variant, every estimator, every box, and the gate quantities.

        task (dict): file_name, path, mpp (recorded), roi_index, clicks (list of dicts), check_blocked (bool).

        Returns dict: rows (estimator landings), boxes, checks (one dict of per-ROI gate quantities), and the timing and memory record.
    """
    t_start = time.perf_counter()
    fn, clicks = task['file_name'], task['clicks']
    timing = {}
    t = time.perf_counter()
    rgb0 = ds.load_roi(task['path'])
    mpp = ds.roi_mpp(task['path'])
    assert mpp == task['mpp'], f'{fn}: mpp {mpp} != recorded {task["mpp"]}'
    hem0 = hem_of(rgb0)
    gray0 = cv2.cvtColor(rgb0, cv2.COLOR_RGB2GRAY)
    T0 = gray0 <= WHITE_GRAY
    fl = roi_floors(hem0)
    b0 = b_otsu_of(hem0)
    assert b0 == fl['b_otsu'], f'{fn}: b_otsu_of {b0} != roi_floors {fl["b_otsu"]}'
    ht0 = hem0[T0]
    v0 = np.sort(ht0)
    side = int(round(WINDOW_UM / mpp)) | 1
    d_shift = int(round(SHIFT_UM / mpp))
    M = macenko_basis(rgb0, np.random.default_rng([MACENKO_SEED, task['roi_index']]))
    timing['setup'] = time.perf_counter() - t
    h, w = hem0.shape
    checks = dict(file_name=fn, roi_index=task['roi_index'], height=h, width=w, mpp=mpp, n_tissue=int(T0.sum()), tissue_share=float(T0.mean()), hem0_min=float(hem0.min()), hem0_max=float(hem0.max()), b0=b0, b0_bin_width=bin_width(hem0), window_side=side, shift_px=d_shift)
    checks.update({f'h{i}': float(M[i, 0]) for i in range(3)})
    checks.update({f'e{i}': float(M[i, 1]) for i in range(3)})
    checks.update(h_resp=float(M[:, 0] @ HCOL), e_resp=float(M[:, 1] @ HCOL))
    checks['label_ok'] = bool(checks['h_resp'] > checks['e_resp'])
    zero = (rgb0 == 0).any(axis=2)
    checks.update(zero_px_share=float(zero.mean()), zero_px_tissue_share=float((zero & T0).sum() / T0.sum()))
    del zero
    for c in LADDER_CHECK:
        checks[f'selfcons_{c:g}'] = c_eq(ht0, v0, b0, float(np.float32(c * b0)))
    # the stride-grid sample for gate 5 and the residual share
    n_px = h * w
    idx = np.arange(0, n_px, max(1, n_px // DIRECTION_SAMPLE))
    s_rgb = rgb0.reshape(-1, 3)[idx]
    s_od = od10(s_rgb)
    s_C = s_od @ np.linalg.pinv(M).T
    s_R = s_od - s_C @ M.T
    s_tis = T0.reshape(-1)[idx]
    checks['resid_share'] = float(np.linalg.norm(s_R[s_tis]) / np.linalg.norm(s_od[s_tis]))
    sample = dict(idx=idx, rgb=s_rgb, hem=hem0.reshape(-1)[idx], tissue=s_tis, C=s_C)
    del s_od, s_R
    # ID: full reconstruction must re-encode bit for bit
    t = time.perf_counter()
    checks['id_n_diff'] = int(np.count_nonzero(synth_identity(rgb0, M) != rgb0))
    timing['ID'] = time.perf_counter() - t

    rows, boxes = [], []
    # the unperturbed reference, and SHIFT on it
    t = time.perf_counter()
    rows += estimator_rows('UNPERT', hem0, gray0, ht0, v0, b0, clicks, side)
    boxes += box_rows('UNPERT', hem0, clicks, mpp)
    for r in clicks:
        x, y, dname, glass, area = shift_window_centre(gray0, r['click_cx'], r['click_cy'], side, d_shift)
        b, n = window_floor(hem0, x, y, side)
        rows.append(dict(var='SHIFT', est='W128', seed_index=r['seed_index'], b=b, c_eq=c_eq(ht0, v0, b0, b), bin_width=bin_width(hem0[window_slices(hem0.shape, x, y, side)]), win_px=n, tissue_share=1.0 - glass, shift_dir=dname, shift_area_share=area))
    timing['UNPERT+SHIFT'] = time.perf_counter() - t
    # NULL2: the exact channel-level null
    t = time.perf_counter()
    hem1 = hem0 * np.float32(2.0)
    rows += estimator_rows('NULL2', hem1, gray0, hem1[T0], v0, b0, clicks, side)
    boxes += box_rows('NULL2', hem1, clicks, mpp)
    del hem1
    timing['NULL2'] = time.perf_counter() - t
    # every RGB-level variant
    for var in IMAGE_VARIANTS:
        t = time.perf_counter()
        rgb1 = make_rgb_variant(var, rgb0, M)
        hem1 = hem_of(rgb1)
        gray1 = cv2.cvtColor(rgb1, cv2.COLOR_RGB2GRAY)
        rows += estimator_rows(var, hem1, gray1, hem1[T0], v0, b0, clicks, side)
        boxes += box_rows(var, hem1, clicks, mpp)
        if var in STAIN_K:
            slope, r_xy, n_sel = direction_slope(sample, rgb1, hem1, M, *STAIN_K[var])
            checks.update({f'dir_slope_{var}': slope, f'dir_r_{var}': r_xy, f'dir_n_{var}': n_sel})
        if var == 'H06':
            checks['zero_px_after_H06'] = float((rgb1 == 0).any(axis=2).mean())
            h06_rgb = rgb1 if task.get('check_blocked') else None
        del rgb1, hem1, gray1
        timing[var] = time.perf_counter() - t
    checks['peak_rss_gb_work'] = peak_rss_gb()
    checks['threads'] = thread_report()
    # smoke only: the blocked channel must equal the whole-image production call (this check's own memory peak comes after the RSS record above)
    if task.get('check_blocked'):
        checks['blocked_n_diff_orig'] = int(np.count_nonzero(ch.to_channel(rgb0, 'hematoxylin_od') != hem0))
        checks['blocked_n_diff_H06'] = int(np.count_nonzero(ch.to_channel(h06_rgb, 'hematoxylin_od') != hem_of(h06_rgb)))
    for row in rows:
        row.update(file_name=fn, roi_index=task['roi_index'])
    for row in boxes:
        row['file_name'] = fn
    timing['total'] = time.perf_counter() - t_start
    return dict(file_name=fn, rows=rows, boxes=boxes, checks=checks, timing=timing, pid=os.getpid())
