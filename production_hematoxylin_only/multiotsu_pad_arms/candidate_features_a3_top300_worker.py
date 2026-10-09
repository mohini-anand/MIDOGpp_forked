"""
    Per-candidate feature extraction on the adopted A3 search's top 300: the worker module for candidate_features_a3_top300.ipynb.

    The spec is CANDIDATE_FEATURE_EXTRACTION_PROMPT_V2.md at the repo root. One task is one ROI: the deep A3 search of
    Section 3.1 on its 6 recorded clicks, the labels of Section 4, and the features of Section 5 once per distinct
    location. The worker writes nothing; the notebook writes every output.

    Spawn workers on macOS cannot unpickle functions defined in a notebook cell, so this code lives in a module (the
    Stage A pattern, tm_floor_stress_stageA_worker.py). The BLAS/OpenMP thread limits are set here, before numpy is
    imported.

    The copied code is not retyped. At import, the module extracts with ast and execs:
    - refine_template_box from cell 3 of multiotsu_template_refinement_14roi.ipynb, with the cell-1 constants it reads;
    - every def in cell 5 and search_arm from cell 11 of tm_floor_definition_49roi_6seed.ipynb, with the cell-1
      constants they read;
    - refusal_class and resample_closed from cell 8 of candrefine_rerank_mo_pad61_select.ipynb.
    The only new search code is detect_deep (written out below) and search_arm_deep (search_arm's own source with its
    one detect( call renamed).
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'

import ast
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import scipy.ndimage as ndi
from scipy.spatial import cKDTree
from skimage.color import rgb2hed, rgb_from_hed
from skimage.filters import threshold_multiotsu
from skimage.measure import label, regionprops
from skimage.morphology import disk
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from midog_utils import channels as ch  # noqa: E402
from midog_utils import chromatin as cm  # noqa: E402
from midog_utils import dataset as ds  # noqa: E402
from midog_utils import evaluate as ev  # noqa: E402
from midog_utils import invariants as inv  # noqa: E402
from midog_utils import production as prod  # noqa: E402
from midog_utils import seed_selection as ss  # noqa: E402,F401  (refine_template_box reads ss._odd)
from midog_utils import template_match as tm  # noqa: E402
from midog_utils.nms import nms_by_distance  # noqa: E402

# ---------------------------------------------------------------------------------------------------------------------
# verbatim code, extracted at import (Section 3)
# ---------------------------------------------------------------------------------------------------------------------
SRC_REFINE = HERE / 'multiotsu_template_refinement_14roi.ipynb'
SRC_SEARCH = HERE / 'tm_floor_definition_49roi_6seed.ipynb'
SRC_CANDREFINE = HERE / 'candrefine_rerank_mo_pad61_select.ipynb'
REFINE_SHA1_EXPECTED = '21c43a20db00dccf7ac3bab88aee423c345441b6'
REFINE_CONSTANTS_SPEC = ('WINDOW', 'CLASSES', 'DILATE_KERNEL', 'DILATE_ITERATIONS', 'MIN_PIECE_UM2', 'MIN_BASE_UM', 'CLICK_RADIUS_UM', 'MIN_AREA_UM2', 'MAX_AREA_FRAC', 'PAD_PX')


def _cell(path, idx):
    """
        One notebook cell's source.

        path (Path): the notebook.
        idx (int): the cell index.

        Returns str: the cell source.
    """
    with open(path) as fh:
        return ''.join(json.load(fh)['cells'][idx]['source'])


def _def_segments(src, names=None):
    """
        The source segments of a cell's top-level function definitions, in source order.

        src (str): the cell source.
        names (tuple or None): the names to keep; None keeps every def.

        Returns dict: name -> (segment, ast node).
    """
    out = {}
    for node in ast.parse(src).body:
        if isinstance(node, ast.FunctionDef) and (names is None or node.name in names):
            out[node.name] = (ast.get_source_segment(src, node), node)
    assert names is None or set(out) == set(names), f'missing defs: {set(names) - set(out)}'
    return out


def _loaded_names(node):
    """
        Every name an AST node reads.

        node (ast.AST): the node.

        Returns set: the names in Load context.
    """
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}


def _constant_segments(src, wanted):
    """
        The top-level single-name assignments of a cell that a set of names needs, with their own dependencies, in source order.

        src (str): the cell source.
        wanted (set): the names the copied functions read.

        Returns dict: name -> assignment source, for every assigned name in the closure of wanted.
    """
    assigns = {}
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            assigns[node.targets[0].id] = node
    need, todo = set(), set(wanted) & set(assigns)
    while todo:
        name = todo.pop()
        need.add(name)
        todo |= (_loaded_names(assigns[name].value) & set(assigns)) - need
    order = sorted(need, key=lambda n: assigns[n].lineno)
    return {n: ast.get_source_segment(src, assigns[n]) for n in order}


def _exec(segment, origin):
    """
        Execute one extracted segment into this module's namespace.

        segment (str): the source.
        origin (str): a label for tracebacks.

        Returns None.
    """
    exec(compile(segment, origin, 'exec'), globals())


# refine_template_box and the constants it reads (cells 1 and 3 of the refinement notebook)
_REFINE_CELL = _cell(SRC_REFINE, 3)
_REFINE_SEG, _REFINE_NODE = _def_segments(_REFINE_CELL, ('refine_template_box',))['refine_template_box']
REFINE_SHA1 = hashlib.sha1((_REFINE_SEG + '\n').encode()).hexdigest()
REFINE_CONSTANTS = _constant_segments(_cell(SRC_REFINE, 1), _loaded_names(_REFINE_NODE))
assert tuple(REFINE_CONSTANTS) == REFINE_CONSTANTS_SPEC, f'refine_template_box reads {tuple(REFINE_CONSTANTS)}, not the spec list'
for _name, _seg in REFINE_CONSTANTS.items():
    _exec(_seg, f'<{SRC_REFINE.name} cell 1: {_name}>')
_exec(_REFINE_SEG, f'<{SRC_REFINE.name} cell 3>')

# the search: every def in cell 5, search_arm from cell 11, and the cell-1 constants they read
_SEARCH_DEFS = _def_segments(_cell(SRC_SEARCH, 5))
_SEARCH_ARM_SEG, _SEARCH_ARM_NODE = _def_segments(_cell(SRC_SEARCH, 11), ('search_arm',))['search_arm']
_search_reads = set().union(*(_loaded_names(n) for _, n in _SEARCH_DEFS.values())) | _loaded_names(_SEARCH_ARM_NODE)
SEARCH_CONSTANTS = _constant_segments(_cell(SRC_SEARCH, 1), _search_reads | {'RANK_KEY'})   # RANK_KEY: the loop glue ranks by it, as cell 11's tp_rows does
for _name, _seg in SEARCH_CONSTANTS.items():
    _exec(_seg, f'<{SRC_SEARCH.name} cell 1: {_name}>')
for _name, (_seg, _) in _SEARCH_DEFS.items():
    _exec(_seg, f'<{SRC_SEARCH.name} cell 5: {_name}>')
_exec(_SEARCH_ARM_SEG, f'<{SRC_SEARCH.name} cell 11: search_arm>')
SEARCH_DEF_NAMES = tuple(_SEARCH_DEFS) + ('search_arm',)

# search_arm_deep: search_arm's own source, with its one detect( call renamed (Section 3.1)
assert _SEARCH_ARM_SEG.count('detect(') == 1 and _SEARCH_ARM_SEG.count('def search_arm(') == 1
SEARCH_ARM_DEEP_SRC = _SEARCH_ARM_SEG.replace('def search_arm(', 'def search_arm_deep(').replace('detect(', 'detect_deep(')
_exec(SEARCH_ARM_DEEP_SRC, '<search_arm_deep>')

# refusal_class and resample_closed from the earlier run (cell 8)
_CANDREFINE_DEFS = _def_segments(_cell(SRC_CANDREFINE, 8), ('refusal_class', 'resample_closed'))
for _name, (_seg, _) in _CANDREFINE_DEFS.items():
    _exec(_seg, f'<{SRC_CANDREFINE.name} cell 8: {_name}>')

assert BAND_MODE == 'floored' and RANK_KEY == 'tm_score' and TPL_ARM == 'mo_pad61' and S_TPL == 61 and PAD == 30  # noqa: F821
assert OD_PAD == prod.OD_WINDOW // 2  # noqa: F821

# ---------------------------------------------------------------------------------------------------------------------
# the deep search (Section 3.1)
# ---------------------------------------------------------------------------------------------------------------------
UNCAPPED = 10 ** 7
DEPTH = 300
NMS_PREFIX0 = 4 * DEPTH      # first prefix of the score-sorted peak list that NMS reads; doubled until DEPTH detections are kept


def detect_deep(fused, valid, band, band_mode, od_padded, od_pad, nms_radius):
    """
        detect with no peak cap, but NMS only as deep as the first DEPTH detections need: production's 100-peak pool as a checked prefix, and od on the DEPTH rows only.

        Greedy NMS keeps a peak only by comparing it with higher-scored peaks, so NMS on the first M peaks of the score-sorted list
        returns exactly the full list's kept peaks among those M. M starts at NMS_PREFIX0 and doubles until DEPTH are kept, so the
        DEPTH rows equal those of an NMS over every peak, without reading the rest of the list (the user's instruction, 2026-10-09).

        fused (np.ndarray): the fused map from fused_map.
        valid (np.ndarray): its reachability mask.
        band (np.ndarray): the border band mask.
        band_mode (str): 'production', 'invalid' or 'floored', as in detect.
        od_padded (np.ndarray): the UNCLIPPED, UNBLANKED hematoxylin_od channel padded by od_pad.
        od_pad (int): the od channel's replicate pad.
        nms_radius (float): this ROI's NMS radius in pixels.

        Returns tuple[pd.DataFrame, dict, np.ndarray]: the first DEPTH post-NMS detections (cx, cy, score, od), the search's counts and statistics (n_nms_prefix: the peaks NMS read), and the pre-NMS peak centres.
    """
    h, w = fused.shape
    if band_mode == 'production':
        f, v = fused, valid
    else:
        assert band_mode in ('invalid', 'floored'), band_mode
        v = valid & ~band
        f = fused if band_mode == 'invalid' else np.where(band, np.float32(tm._FLOOR), fused)
    med, mad = tm.robust_stats(f, v)
    threshold = med + prod.DEEP_FLOOR_Z * mad
    centers, scores = tm.extract_peaks(f, v, prod.PEAK_MIN_DISTANCE, threshold, UNCAPPED)
    peaks = centers.copy()
    peak_scores = np.asarray(scores, dtype=np.float64)
    m = min(NMS_PREFIX0, len(centers))
    while True:
        keep = nms_by_distance(centers[:m], scores[:m], nms_radius)
        if len(keep) >= DEPTH or m == len(centers):
            break
        m = min(2 * m, len(centers))
    keep_rec = nms_by_distance(centers[:prod.MAX_PEAKS], scores[:prod.MAX_PEAKS], nms_radius)
    assert np.array_equal(keep[:len(keep_rec)], keep_rec), 'production pool is not a prefix of the deep NMS list'
    assert len(keep) >= DEPTH, f'only {len(keep)} post-NMS detections'
    n_nms_prefix = m
    keep = keep[:DEPTH]
    centers, scores = centers[keep], scores[keep]
    det = pd.DataFrame({'rank': np.arange(len(centers)), 'cx': centers[:, 0], 'cy': centers[:, 1], 'score': scores})
    assert od_padded.shape[:2] == (h + 2 * od_pad, w + 2 * od_pad), f'the od channel is not padded by {od_pad}'
    shifted = det.assign(cx=det['cx'] + od_pad, cy=det['cy'] + od_pad)
    det = det.assign(od=cm.score_detections(shifted, od_padded, window=prod.OD_WINDOW)['od'].to_numpy())
    info = dict(n_peaks=len(peaks), n_nms_prefix=n_nms_prefix, n_detections_recorded=len(keep_rec), n_detections=len(det), threshold=float(threshold), median=float(med), mad=float(mad), score_peak_100=float(peak_scores[prod.MAX_PEAKS - 1]), score_peak_101=float(peak_scores[prod.MAX_PEAKS]) if len(peak_scores) > prod.MAX_PEAKS else np.nan)
    return det, info, peaks


# ---------------------------------------------------------------------------------------------------------------------
# feature constants (Sections 4 and 5); every one goes into the spec JSON
# ---------------------------------------------------------------------------------------------------------------------
MIN_PIECE_PX = 30            # 'small piece' below this many pixels (5.2)
CONCAVITY_MIN_UM = 0.25      # convexity defects deeper than this count (5.2)
RADIAL_POINTS = 128          # outline resampling of chr_radial_hf (5.2)
RADIAL_BAND = (6, 64)        # harmonics 6 to 64 inclusive (5.2)
PAIR_MIN_UM2 = 1.0           # pair pieces at least this large (5.2)
PAIR_MIN_RATIO = 0.5         # partner area at least this share of the anchor's (5.2)
PAIR_MAX_UM = 10.0           # partner centroid within this distance (5.2)
PAIR_ELONG_MIN = 1.2         # angle columns NaN below this major/minor ratio (5.2)
WIN101 = 101                 # pair, context and neighbour window (Section 4)
WIN61 = 61                   # white_frac_gray215_61's window
NS_MIN_PX = 10               # minimum support for texture, cytoplasm ring and chr_resid_frac (Sections 4, 5.4, 5.5, 5.8 a)
GLASS_GRAY = 215             # gray > 215 is glass or lumen (5.5), the floor-clip run's WHITE_GRAY
BROWN_R_MINUS_B = 40         # brown rule (5.6)
BROWN_G_MINUS_B = 15
BROWN_GRAY_MAX = 180
RESID_DILATE_PX = 2          # chr_resid_frac's rim (5.8 a)
NS_LO, NS_HI = 20, 235       # chr_resid_frac's non-saturated range (5.8 a)
TOPHAT_R_UM = 0.5            # top-hat disc radius (5.8 b)
TOPHAT_RING_UM = 1.0         # top-hat region reaches this far outside M_chr (5.8 b)
TOPHAT_CORE_UM = 0.5         # M_core: more than this inside M_chr's edge (5.8 b)
EDGE_MIN_STEP = 0.005        # chr_edge_sharpness_per_um needs at least this OD step (5.8 c)
NBR_MIN_UM2, NBR_MAX_UM2 = 8.0, 150.0   # neighbour component area bounds (5.8 d)
NBR_MIN_N = 3                # nbr_nuc_area_med_um2 needs at least this many neighbours (5.8 d)
CROSS = ndi.generate_binary_structure(2, 1)

_H = rgb_from_hed[0] / np.linalg.norm(rgb_from_hed[0])
_E = rgb_from_hed[1] / np.linalg.norm(rgb_from_hed[1])
_N = np.cross(_H, _E)
_N = _N / np.linalg.norm(_N)
RESID_NORMAL = _N if _N[2] > 0 else -_N

FEATURES = ('chr_area_um2', 'chr_elongation', 'chr_solidity', 'chr_concavities_per_um', 'dark_concavities_per_um', 'chr_hole_frac', 'chr_radial_hf', 'pair_sep_um', 'pair_area_ratio', 'pair_parallel_deg', 'pair_axis_deg', 'h_chr_mean', 'h_chr_lap', 'ring_peak_prominence', 'ring1_frac_above_tnuc', 'ring1_frac_pale', 'white_frac_gray215_61', 'e_h_angle_ring23_deg', 'brown_frac_surround', 'dark_area_frac_101', 'chr_resid_frac', 'chr_tophat_edge', 'chr_tophat_core', 'chr_area_over_nbr')
DIAGNOSTICS = ('chr_touches_window', 'dark_touches_window', 'h_ring1_mean', 'h_ring2_mean', 'h_ring3_mean', 'chr_ns_share', 'chr_edge_sharpness_per_um', 'nbr_n', 'nbr_nuc_area_med_um2')
REASONS = ('nan_reason_geom', 'nan_reason_concav', 'nan_reason_dark', 'nan_reason_pair', 'nan_reason_texture', 'nan_reason_cyto', 'nan_reason_resid', 'nan_reason_ring', 'nan_reason_edge', 'nan_reason_nbr', 'nan_reason_spec_nbr')
# spec_-prefixed columns of definitions replaced during the run (Section 2). 5.8 (d)'s neighbour components at t_nuc were not nuclei on the
# real-data check (cytoplasm and stroma fragments on 462 and 300; every nucleus merged into edge-touching sheets on 524), so the corrected
# columns use t_dark and the spec's t_nuc version is kept here, with its own reason column nan_reason_spec_nbr.
SPEC_COLUMNS = ('spec_chr_area_over_nbr', 'spec_nbr_n', 'spec_nbr_nuc_area_med_um2')
NBR_THRESHOLD = 't_dark'     # the corrected neighbour threshold (spec: t_nuc, kept as the spec_ columns)


def blank_row():
    """
        A feature row with every feature and diagnostic NaN and every reason empty, as a refused row has it.

        Returns dict: column -> value.
    """
    return dict({c: np.nan for c in FEATURES + DIAGNOSTICS + SPEC_COLUMNS}, **{r: '' for r in REASONS})


# ---------------------------------------------------------------------------------------------------------------------
# masks and pieces
# ---------------------------------------------------------------------------------------------------------------------
def largest_piece(mask):
    """
        The largest 4-connected piece of a mask, a tie going to the lower label.

        mask (np.ndarray): a boolean mask.

        Returns np.ndarray: the piece as a boolean mask (all False when the mask is empty).
    """
    lab = label(mask, connectivity=1)
    if lab.max() == 0:
        return np.zeros(mask.shape, bool)
    counts = np.bincount(lab.ravel())
    counts[0] = 0
    return lab == int(np.argmax(counts))


def touches_edge(mask):
    """
        Whether a mask reaches its array's outer edge.

        mask (np.ndarray): a boolean mask.

        Returns bool: True when any edge pixel is set.
    """
    return bool(mask[0].any() or mask[-1].any() or mask[:, 0].any() or mask[:, -1].any())


def outline_of(filled):
    """
        The longest external contour of a filled mask, as an int32 (n, 1, 2) array of (x, y) points.

        filled (np.ndarray): a boolean mask with holes filled.

        Returns np.ndarray: the contour.
    """
    contours, _ = cv2.findContours(filled.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    return max(contours, key=len).astype(np.int32)


def concavities_per_um(filled, mpp):
    """
        The number of convexity defects deeper than CONCAVITY_MIN_UM per µm of outline (5.2).

        filled (np.ndarray): the hole-filled piece.
        mpp (float): microns per pixel.

        Returns tuple[float, str]: the value (NaN on an OpenCV error) and the NaN reason ('' or 'convexity error').
    """
    c = outline_of(filled)
    try:
        d = cv2.convexityDefects(c, cv2.convexHull(c, returnPoints=False))
    except cv2.error:
        return np.nan, 'convexity error'
    n = 0 if d is None else int((d[:, 0, 3] / 256.0 * mpp > CONCAVITY_MIN_UM).sum())
    return n / (cv2.arcLength(c, True) * mpp), ''


def radial_hf(filled):
    """
        chr_radial_hf: the in-band (harmonics 6-64) energy of the mean-centred radius profile over 128 * mean radius squared (5.2).

        filled (np.ndarray): the hole-filled piece.

        Returns float: the ratio.
    """
    o = resample_closed(outline_of(filled)[:, 0, :], RADIAL_POINTS)  # noqa: F821  (verbatim, cell 8)
    ys, xs = np.nonzero(filled)
    r = np.hypot(o[:, 0] - xs.mean(), o[:, 1] - ys.mean())
    e = np.abs(np.fft.rfft(r - r.mean())) ** 2
    return float(e[RADIAL_BAND[0]:RADIAL_BAND[1] + 1].sum() / (RADIAL_POINTS * r.mean() ** 2))


# ---------------------------------------------------------------------------------------------------------------------
# feature families (Section 5); each returns a dict of its columns and reason columns
# ---------------------------------------------------------------------------------------------------------------------
def geometry_cols(m_chr, mpp):
    """
        Section 5.2's M_chr columns: area (whole), elongation, solidity, concavities, hole share, radial energy, window touch.

        m_chr (np.ndarray): the undilated chromatin mask.
        mpp (float): microns per pixel.

        Returns dict: the columns, nan_reason_geom and nan_reason_concav.
    """
    out = dict(chr_area_um2=float(m_chr.sum()) * mpp ** 2, chr_touches_window=touches_edge(m_chr), chr_elongation=np.nan, chr_solidity=np.nan, chr_concavities_per_um=np.nan, chr_hole_frac=np.nan, chr_radial_hf=np.nan, nan_reason_geom='', nan_reason_concav='')
    piece = largest_piece(m_chr)
    if piece.sum() < MIN_PIECE_PX:
        out.update(nan_reason_geom='small piece', nan_reason_concav='small piece')
        return out
    filled = ndi.binary_fill_holes(piece)
    rp = regionprops(filled.astype(np.uint8))[0]
    out.update(chr_elongation=float(rp.axis_major_length / rp.axis_minor_length) if rp.axis_minor_length > 0 else np.nan, chr_solidity=float(rp.solidity), chr_hole_frac=float(1.0 - piece.sum() / filled.sum()), chr_radial_hf=radial_hf(filled))
    out['chr_concavities_per_um'], out['nan_reason_concav'] = concavities_per_um(filled, mpp)
    return out


def dark_cols(h_w, m_chr, t_dark, mpp):
    """
        dark_concavities_per_um on M_dark: the largest 4-connected component of H_w > t_dark overlapping M_chr, holes filled (Section 4, 5.2).

        h_w (np.ndarray): the 51 px window's unclipped H.
        m_chr (np.ndarray): the undilated chromatin mask.
        t_dark (float): the ROI's upper 3-class multi-Otsu threshold.
        mpp (float): microns per pixel.

        Returns dict: dark_concavities_per_um, dark_touches_window and nan_reason_dark.
    """
    comp = label(h_w > t_dark, connectivity=1)
    ids = np.unique(comp[m_chr])
    ids = ids[ids > 0]
    if len(ids) == 0:
        return dict(dark_concavities_per_um=np.nan, dark_touches_window=np.nan, nan_reason_dark='small piece')
    sizes = np.bincount(comp.ravel())[ids]
    piece = comp == int(ids[int(np.argmax(sizes))])
    touch = touches_edge(piece)
    if piece.sum() < MIN_PIECE_PX:
        return dict(dark_concavities_per_um=np.nan, dark_touches_window=touch, nan_reason_dark='small piece')
    value, reason = concavities_per_um(ndi.binary_fill_holes(piece), mpp)
    return dict(dark_concavities_per_um=value, dark_touches_window=touch, nan_reason_dark=reason)


def ring_masks(m_chr, mpp):
    """
        ring_1, ring_2, ring_3: the window pixels whose distance from M_chr lies in ((k - 1) µm, k µm] (Section 4).

        m_chr (np.ndarray): the undilated chromatin mask.
        mpp (float): microns per pixel.

        Returns list[np.ndarray]: the three ring masks.
    """
    dist = ndi.distance_transform_edt(~m_chr) * mpp
    return [(dist > k - 1) & (dist <= k) for k in (1, 2, 3)]


def darkness_cols(h_w, m_chr, rings):
    """
        Section 5.3: h_chr_mean and the three ring means (diagnostics).

        h_w (np.ndarray): the 51 px window's unclipped H.
        m_chr (np.ndarray): the undilated chromatin mask.
        rings (list): ring_masks.

        Returns dict: h_chr_mean, h_ring1_mean, h_ring2_mean, h_ring3_mean.
    """
    h = h_w.astype(np.float64)
    out = dict(h_chr_mean=float(h[m_chr].mean()))
    for k, ring in enumerate(rings, start=1):
        out[f'h_ring{k}_mean'] = float(h[ring].mean()) if ring.any() else np.nan
    return out


def texture_cols(h_w, m_chr, clipped, mpp, h_chr_mean):
    """
        Section 5.4: h_chr_lap, the std of the 4-neighbour Laplacian over the eroded, clip-free support, over h_chr_mean.

        h_w (np.ndarray): the 51 px window's unclipped H, float32.
        m_chr (np.ndarray): the undilated chromatin mask.
        clipped (np.ndarray): pixels with any RGB channel at 0 or 255.
        mpp (float): microns per pixel.
        h_chr_mean (float): mean H over M_chr.

        Returns dict: h_chr_lap and nan_reason_texture.
    """
    eroded = ndi.binary_erosion(m_chr)
    support = eroded & ~ndi.binary_dilation(clipped, structure=np.ones((3, 3), bool))
    if eroded.sum() < NS_MIN_PX:
        return dict(h_chr_lap=np.nan, nan_reason_texture='small support')
    if support.sum() < NS_MIN_PX:
        return dict(h_chr_lap=np.nan, nan_reason_texture='clipped')
    lap = cv2.Laplacian(np.ascontiguousarray(h_w, dtype=np.float32), cv2.CV_32F, ksize=1) / mpp ** 2
    return dict(h_chr_lap=float(lap[support].astype(np.float64).std(ddof=0) / h_chr_mean), nan_reason_texture='')


def surround_cols(h_w, gray_w, e_w, clipped, rings, dk, t_nuc, t_dark):
    """
        Section 5.5 (except the 61 px white share): ring prominence, the ring-1 shares and the cytoplasm hue angle.

        h_w (np.ndarray): the 51 px window's unclipped H.
        gray_w (np.ndarray): the window's uint8 gray.
        e_w (np.ndarray): the window's rgb2hed E.
        clipped (np.ndarray): pixels with any RGB channel at 0 or 255.
        rings (list): ring_masks.
        dk (dict): darkness_cols' output.
        t_nuc (float): the ROI's lower 3-class multi-Otsu threshold.
        t_dark (float): the ROI's upper threshold.

        Returns dict: ring_peak_prominence, ring1_frac_above_tnuc, ring1_frac_pale, e_h_angle_ring23_deg, nan_reason_ring, nan_reason_cyto.
    """
    out = dict(ring_peak_prominence=np.nan, ring1_frac_above_tnuc=np.nan, ring1_frac_pale=np.nan, e_h_angle_ring23_deg=np.nan, nan_reason_ring='', nan_reason_cyto='')
    hr = [dk['h_ring1_mean'], dk['h_ring2_mean'], dk['h_ring3_mean']]
    if any(not ring.any() for ring in rings):
        out['nan_reason_ring'] = 'empty ring'
    else:
        p2 = hr[1] - hr[0]
        p3 = hr[2] - min(hr[0], hr[1])
        out['ring_peak_prominence'] = max(0.0, max(p2, p3) / dk['h_chr_mean'])
    r1 = rings[0]
    if r1.any():
        h1, g1 = h_w[r1], gray_w[r1]
        out.update(ring1_frac_above_tnuc=float((h1 > t_nuc).mean()), ring1_frac_pale=float(((h1 <= t_nuc) & (g1 <= GLASS_GRAY)).mean()))
    cyto = (rings[1] | rings[2]) & ~(gray_w > GLASS_GRAY) & ~(h_w > t_dark) & ~clipped
    if cyto.sum() < NS_MIN_PX:
        out['nan_reason_cyto'] = 'no cytoplasm'
        return out
    me, mh = float(e_w[cyto].astype(np.float64).mean()), float(h_w[cyto].astype(np.float64).mean())
    if me == 0 and mh == 0:
        out['nan_reason_cyto'] = 'no cytoplasm'
        return out
    out['e_h_angle_ring23_deg'] = float(np.degrees(np.arctan2(me, mh)))
    return out


def brown_col(rgb_w, gray_w, m_sel):
    """
        Section 5.6: brown_frac_surround, the share of window pixels outside M_sel passing the four-term brown rule.

        rgb_w (np.ndarray): the window's raw uint8 RGB.
        gray_w (np.ndarray): the window's uint8 gray.
        m_sel (np.ndarray): the dilated selection.

        Returns dict: brown_frac_surround.
    """
    r, g, b = (rgb_w[:, :, k].astype(np.int16) for k in range(3))
    brown = (r > g) & (g > b) & (r - b > BROWN_R_MINUS_B) & (g - b > BROWN_G_MINUS_B) & (gray_w < BROWN_GRAY_MAX)
    out_sel = ~m_sel
    return dict(brown_frac_surround=float(brown[out_sel].sum() / out_sel.sum()))


def resid_cols(rgb_w, m_chr):
    """
        Section 5.8 (a): chr_resid_frac, the mean OD's component off the H-E plane over its norm, on M_chr plus a 2 px rim, non-saturated pixels only.

        rgb_w (np.ndarray): the window's raw uint8 RGB.
        m_chr (np.ndarray): the undilated chromatin mask (whole).

        Returns dict: chr_resid_frac, chr_ns_share and nan_reason_resid.
    """
    m_res = ndi.binary_dilation(m_chr, structure=disk(RESID_DILATE_PX).astype(bool))
    ns = m_res & np.all((rgb_w >= NS_LO) & (rgb_w <= NS_HI), axis=2)
    out = dict(chr_ns_share=float(ns.sum() / m_res.sum()), chr_resid_frac=np.nan, nan_reason_resid='')
    if ns.sum() < NS_MIN_PX:
        out['nan_reason_resid'] = 'outside_20_235'
        return out
    od = -np.log10(np.maximum(rgb_w.astype(np.float64), 1.0) / 255.0)
    v = od[ns].mean(axis=0)
    out['chr_resid_frac'] = float((v @ RESID_NORMAL) / np.linalg.norm(v))
    return out


def tophat_cols(h_w, m_chr, mpp):
    """
        Section 5.8 (b): thin dark structure (a 0.5 µm white top-hat) within 1 µm of the core's surface and inside the core, each over the core's total H.

        h_w (np.ndarray): the window's unclipped H.
        m_chr (np.ndarray): the undilated chromatin mask (whole).
        mpp (float): microns per pixel.

        Returns dict: chr_tophat_edge, chr_tophat_core (NaN 'small piece' in nan_reason_geom when whole M_chr < 30 px).
    """
    if m_chr.sum() < MIN_PIECE_PX:
        return dict(chr_tophat_edge=np.nan, chr_tophat_core=np.nan, tophat_small=True)
    r_th = max(1, int(round(TOPHAT_R_UM / mpp)))
    se = disk(r_th).astype(np.uint8)
    hf = np.ascontiguousarray(h_w, dtype=np.float32)
    t = (hf - cv2.morphologyEx(hf, cv2.MORPH_OPEN, se, borderType=cv2.BORDER_REPLICATE)).astype(np.float64)
    region = ndi.distance_transform_edt(~m_chr) * mpp <= TOPHAT_RING_UM
    core = ndi.distance_transform_edt(m_chr) * mpp > TOPHAT_CORE_UM
    total = float(hf[m_chr].astype(np.float64).sum())
    return dict(chr_tophat_edge=float(t[region & ~core].sum() / total), chr_tophat_core=float(t[core].sum() / total), tophat_small=False)


def edge_col(h_w, m_chr, mpp, dk):
    """
        Section 5.8 (c): chr_edge_sharpness_per_um, the mean Sobel gradient on the filled largest piece's boundary over the chromatin-to-ring-1 step.

        h_w (np.ndarray): the window's unclipped H.
        m_chr (np.ndarray): the undilated chromatin mask.
        mpp (float): microns per pixel.
        dk (dict): darkness_cols' output.

        Returns dict: chr_edge_sharpness_per_um and which reason applies ('geom', 'ring', 'edge' or '').
    """
    piece = largest_piece(m_chr)
    if piece.sum() < MIN_PIECE_PX:
        return dict(chr_edge_sharpness_per_um=np.nan, edge_reason='geom')
    if not np.isfinite(dk['h_ring1_mean']):
        return dict(chr_edge_sharpness_per_um=np.nan, edge_reason='ring')
    step = dk['h_chr_mean'] - dk['h_ring1_mean']
    if step < EDGE_MIN_STEP:
        return dict(chr_edge_sharpness_per_um=np.nan, edge_reason='edge')
    filled = ndi.binary_fill_holes(piece)
    boundary = filled & ~ndi.binary_erosion(filled, structure=CROSS)
    hf = np.ascontiguousarray(h_w, dtype=np.float32)
    gx = cv2.Sobel(hf, cv2.CV_32F, 1, 0, ksize=3) / 8.0
    gy = cv2.Sobel(hf, cv2.CV_32F, 0, 1, ksize=3) / 8.0
    g = np.hypot(gx.astype(np.float64), gy.astype(np.float64)) / mpp
    return dict(chr_edge_sharpness_per_um=float(g[boundary].mean() / step), edge_reason='')


def _acute_deg(u, v):
    """
        The acute angle between two lines with direction vectors u and v.

        u (np.ndarray): a 2-vector.
        v (np.ndarray): a 2-vector.

        Returns float: degrees in [0, 90].
    """
    c = abs(float(np.dot(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v)))
    return float(np.degrees(np.arccos(min(1.0, c))))


def pair_cols(hem_101, m_chr_101, mpp, t_dark):
    """
        Section 5.2's pair columns from the 101 px window: the anchor piece under M_chr and its nearest similar-sized partner, at t_dark.

        hem_101 (np.ndarray): the (truncated) 101 px window's unclipped H.
        m_chr_101 (np.ndarray): M_chr placed in that window.
        mpp (float): microns per pixel.
        t_dark (float): the ROI's upper threshold.

        Returns dict: pair_sep_um, pair_area_ratio, pair_parallel_deg, pair_axis_deg, nan_reason_pair.
    """
    out = dict(pair_sep_um=np.nan, pair_area_ratio=np.nan, pair_parallel_deg=np.nan, pair_axis_deg=np.nan, nan_reason_pair='')
    lab = label(hem_101 > t_dark, connectivity=1)
    props = [p for p in regionprops(lab) if p.area * mpp ** 2 >= PAIR_MIN_UM2]
    if not props:
        out['nan_reason_pair'] = 'no anchor'
        return out
    overlap = {p.label: int(m_chr_101[p.coords[:, 0], p.coords[:, 1]].sum()) for p in props}
    cands = [p for p in props if overlap[p.label] > 0]
    if not cands:
        out['nan_reason_pair'] = 'no anchor'
        return out
    anchor = sorted(cands, key=lambda p: (-overlap[p.label], -p.area, p.label))[0]
    ca = np.asarray(anchor.centroid)
    partners = [(float(np.hypot(*(np.asarray(p.centroid) - ca))), p.label, p) for p in props if p.label != anchor.label and p.area >= PAIR_MIN_RATIO * anchor.area]
    partners = [t for t in partners if t[0] * mpp <= PAIR_MAX_UM]
    if not partners:
        out['nan_reason_pair'] = 'no partner'
        return out
    dist, _, partner = sorted(partners, key=lambda t: (t[0], t[1]))[0]
    out.update(pair_sep_um=dist * mpp, pair_area_ratio=float(min(anchor.area, partner.area) / max(anchor.area, partner.area)))
    ratios = [p.axis_major_length / p.axis_minor_length if p.axis_minor_length > 0 else np.inf for p in (anchor, partner)]
    if min(ratios) < PAIR_ELONG_MIN:
        out['nan_reason_pair'] = 'round piece'
        return out
    axes = [np.array([np.cos(p.orientation), np.sin(p.orientation)]) for p in (anchor, partner)]
    line = np.asarray(partner.centroid) - ca
    out.update(pair_parallel_deg=_acute_deg(axes[0], axes[1]), pair_axis_deg=0.5 * (_acute_deg(axes[0], line) + _acute_deg(axes[1], line)))
    return out


def neighbour_components(hem_101, m_chr_101, mpp, thr):
    """
        Section 5.8 (d)'s components of hem_101 > thr, split into kept neighbours and the three drop rules.

        hem_101 (np.ndarray): the (truncated) 101 px window's unclipped H.
        m_chr_101 (np.ndarray): M_chr placed in that window.
        mpp (float): microns per pixel.
        thr (float): the component threshold: t_dark for the corrected columns, t_nuc for the spec_ columns.

        Returns tuple[np.ndarray, dict]: the label image, and label -> 'kept', 'edge', 'own' or 'area'.
    """
    lab = label(hem_101 > thr, connectivity=1)
    edge_labels = set(np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))) - {0}
    own = set(np.unique(lab[m_chr_101])) - {0}
    area = np.bincount(lab.ravel()) * mpp ** 2
    status = {}
    for k in range(1, lab.max() + 1):
        if k in edge_labels:
            status[k] = 'edge'
        elif k in own:
            status[k] = 'own'
        elif not NBR_MIN_UM2 <= area[k] <= NBR_MAX_UM2:
            status[k] = 'area'
        else:
            status[k] = 'kept'
    return lab, status


def neighbour_cols(hem_101, m_chr_101, mpp, thr, chr_area_um2):
    """
        Section 5.8 (d): the neighbour count, their median area, and the candidate's chromatin area over it.

        hem_101 (np.ndarray): the (truncated) 101 px window's unclipped H.
        m_chr_101 (np.ndarray): M_chr placed in that window.
        mpp (float): microns per pixel.
        thr (float): the component threshold.
        chr_area_um2 (float): the candidate's whole M_chr area.

        Returns dict: nbr_n, nbr_nuc_area_med_um2, chr_area_over_nbr, nan_reason_nbr.
    """
    lab, status = neighbour_components(hem_101, m_chr_101, mpp, thr)
    kept = [k for k, s in status.items() if s == 'kept']
    out = dict(nbr_n=len(kept), nbr_nuc_area_med_um2=np.nan, chr_area_over_nbr=np.nan, nan_reason_nbr='')
    if len(kept) < NBR_MIN_N:
        out['nan_reason_nbr'] = 'few neighbours'
        return out
    area = np.bincount(lab.ravel()) * mpp ** 2
    med = float(np.median(area[kept]))
    out.update(nbr_nuc_area_med_um2=med, chr_area_over_nbr=chr_area_um2 / med)
    return out


def dark_area_col(hem_101, m_sel_101, t_dark):
    """
        Section 5.7: dark_area_frac_101, the share above t_dark among the 101 px window's pixels outside the candidate's own M_sel.

        hem_101 (np.ndarray): the (truncated) 101 px window's unclipped H.
        m_sel_101 (np.ndarray): M_sel placed in that window.
        t_dark (float): the ROI's upper threshold.

        Returns dict: dark_area_frac_101.
    """
    out_sel = ~m_sel_101
    return dict(dark_area_frac_101=float((hem_101[out_sel] > t_dark).sum() / out_sel.sum()))


def window_features(h_w, m_sel, m_chr, rgb_w, gray_w, e_w, hem_101, offset, mpp, t_nuc, t_dark):
    """
        Every feature and diagnostic of one accepted candidate from its windows and masks (Section 5), except the 61 px white share.

        h_w (np.ndarray): the 51 px window's unclipped H (float32).
        m_sel (np.ndarray): the dilated selection.
        m_chr (np.ndarray): the undilated chromatin mask, non-empty.
        rgb_w (np.ndarray): the 51 px window's raw uint8 RGB.
        gray_w (np.ndarray): the 51 px window's uint8 gray.
        e_w (np.ndarray): the 51 px window's rgb2hed E.
        hem_101 (np.ndarray): the (truncated) 101 px window's unclipped H.
        offset (tuple): (row, column) of the 51 px window's origin inside hem_101.
        mpp (float): microns per pixel.
        t_nuc (float): the ROI's lower threshold.
        t_dark (float): the ROI's upper threshold.

        Returns dict: the feature row (features, diagnostics, reasons).
    """
    row = blank_row()
    clipped = np.any((rgb_w == 0) | (rgb_w == 255), axis=2)
    geo = geometry_cols(m_chr, mpp)
    row.update(geo)
    row.update(dark_cols(h_w, m_chr, t_dark, mpp))
    rings = ring_masks(m_chr, mpp)
    dk = darkness_cols(h_w, m_chr, rings)
    row.update(dk)
    row.update(texture_cols(h_w, m_chr, clipped, mpp, dk['h_chr_mean']))
    row.update(surround_cols(h_w, gray_w, e_w, clipped, rings, dk, t_nuc, t_dark))
    row.update(brown_col(rgb_w, gray_w, m_sel))
    row.update(resid_cols(rgb_w, m_chr))
    th = tophat_cols(h_w, m_chr, mpp)
    if th.pop('tophat_small'):
        row['nan_reason_geom'] = 'small piece'
    row.update(th)
    ed = edge_col(h_w, m_chr, mpp, dk)
    reason = ed.pop('edge_reason')
    if reason == 'ring':
        row['nan_reason_ring'] = 'empty ring'
    elif reason == 'edge':
        row['nan_reason_edge'] = 'no edge step'
    row.update(ed)
    oy, ox = offset
    m_chr_101 = np.zeros(hem_101.shape, bool)
    m_sel_101 = np.zeros(hem_101.shape, bool)
    m_chr_101[oy:oy + m_chr.shape[0], ox:ox + m_chr.shape[1]] = m_chr
    m_sel_101[oy:oy + m_sel.shape[0], ox:ox + m_sel.shape[1]] = m_sel
    row.update(pair_cols(hem_101, m_chr_101, mpp, t_dark))
    row.update(dark_area_col(hem_101, m_sel_101, t_dark))
    row.update(neighbour_cols(hem_101, m_chr_101, mpp, t_dark, geo['chr_area_um2']))
    spec = neighbour_cols(hem_101, m_chr_101, mpp, t_nuc, geo['chr_area_um2'])
    row.update(spec_chr_area_over_nbr=spec['chr_area_over_nbr'], spec_nbr_n=spec['nbr_n'], spec_nbr_nuc_area_med_um2=spec['nbr_nuc_area_med_um2'], nan_reason_spec_nbr=spec['nan_reason_nbr'])
    return row


def window_101(img, cx, cy):
    """
        The 101 px window centred on a location, truncated at the ROI edge, and the 51 px window's origin inside it.

        img (np.ndarray): a full-ROI array.
        cx (float): x in ROI pixels.
        cy (float): y in ROI pixels.

        Returns tuple[np.ndarray, tuple, tuple]: the window, (y0, x0) its origin in ROI pixels, and (row, column) of the 51 px window inside it.
    """
    h, w = img.shape[:2]
    ix, iy = int(round(cx)), int(round(cy))
    half = WIN101 // 2
    y0, x0 = max(0, iy - half), max(0, ix - half)
    y1, x1 = min(h, iy + half + 1), min(w, ix + half + 1)
    oy, ox = iy - WINDOW // 2 - y0, ix - WINDOW // 2 - x0  # noqa: F821  (WINDOW is the refinement's 51)
    assert oy >= 0 and ox >= 0 and oy + WINDOW <= y1 - y0 and ox + WINDOW <= x1 - x0, 'the 51 px window is not inside the 101 px window'  # noqa: F821
    return img[y0:y1, x0:x1], (y0, x0), (oy, ox)


def location_features(hem, e_ch, rgb, gray, cx, cy, mpp, t_nuc, t_dark):
    """
        Refine one candidate at its detection coordinates on the unclipped channel and measure Section 5 on it.

        hem (np.ndarray): the ROI's unclipped, unblanked hematoxylin_od.
        e_ch (np.ndarray): the ROI's rgb2hed E channel.
        rgb (np.ndarray): the ROI's raw uint8 RGB.
        gray (np.ndarray): the ROI's uint8 gray.
        cx (float): detection x.
        cy (float): detection y.
        mpp (float): microns per pixel.
        t_nuc (float): the ROI's lower threshold.
        t_dark (float): the ROI's upper threshold.

        Returns tuple[dict, dict]: the row (refinement record, features, diagnostics, reasons), and the masks (empty when refused).
    """
    ref = refine_template_box(hem, cx, cy, mpp)  # noqa: F821  (verbatim)
    row = dict(accepted=bool(ref['accepted']), refusal=refusal_class(ref['reason']), mode=ref['mode'], n_pieces=len(ref['selected']), base_size=np.nan if ref['base_size'] is None else int(ref['base_size']))  # noqa: F821
    row.update(blank_row())
    if not ref['accepted']:
        return row, {}
    m_sel = ref['selection_mask']
    m_chr = m_sel & (ref['class_map'] == 2)
    if not m_chr.any():
        row.update(accepted=False, refusal='empty chromatin mask')
        return row, {}
    h_w = np.asarray(ref['hem_patch'], dtype=np.float32)
    rgb_w = tm.read_padded_patch(rgb, cx, cy, WINDOW)  # noqa: F821
    gray_w = tm.read_padded_patch(gray, cx, cy, WINDOW)  # noqa: F821
    e_w = tm.read_padded_patch(e_ch, cx, cy, WINDOW)  # noqa: F821
    hem_101, _, offset = window_101(hem, cx, cy)
    row.update(window_features(h_w, m_sel, m_chr, rgb_w, gray_w, e_w, hem_101, offset, mpp, t_nuc, t_dark))
    return row, dict(m_sel=m_sel, m_chr=m_chr)


# ---------------------------------------------------------------------------------------------------------------------
# one ROI (the spawn worker's task)
# ---------------------------------------------------------------------------------------------------------------------
def nearest_um(xy, ref_xy, mpp):
    """
        The distance from each point to the nearest reference point, in µm.

        xy (np.ndarray): (n, 2) points.
        ref_xy (np.ndarray): (m, 2) reference points.
        mpp (float): microns per pixel.

        Returns np.ndarray: (n,) distances, NaN when there is no reference point.
    """
    if len(ref_xy) == 0:
        return np.full(len(xy), np.nan)
    return cKDTree(ref_xy).query(xy, k=1)[0] * mpp


def run_roi(task):
    """
        The deep A3 search on one ROI's clicks, the labels of Section 4, and every feature once per distinct location.

        task (dict): file_name, path, mpp, clicks (list of dicts: seed_index, seed_ann_id, click_cx, click_cy, tcx, tcy, bs), gt (pd.DataFrame of the ROI's annotations).

        Returns dict: clicks (list of per-click dicts), cands (pd.DataFrame, 300 rows per click), feats (pd.DataFrame, one row per distinct location), timing.
    """
    cv2.setNumThreads(1)
    with threadpool_limits(1):
        t_roi = time.perf_counter()
        f, path = task['file_name'], task['path']
        rgb = ds.load_roi(path)
        mpp = ds.roi_mpp(path)
        assert mpp == task['mpp'], f'{f}: mpp {mpp} differs from the recorded {task["mpp"]}'
        hem = ch.to_channel(rgb, 'hematoxylin_od')
        fl = roi_floors(hem)  # noqa: F821
        b, t_nuc, t_dark = fl['b_otsu'], fl['b_otsu'], fl['b_otsu_upper']
        nms_radius = match_radius = ev.radius_px(mpp)
        inv.check_nms_radius(nms_radius, mpp, label='candidate_features_a3_top300')
        channel = clip_floor(hem, b)  # noqa: F821
        channel_edge = padded(channel, EDGE_READ_PAD)  # noqa: F821
        od_padded = padded(hem, OD_PAD)  # noqa: F821
        band = band_mask(hem.shape, PAD)  # noqa: F821
        assert_band_identity(band, PAD)  # noqa: F821
        gt = task['gt']
        click_rows, cand_frames = [], []
        for c in task['clicks']:
            t0 = time.perf_counter()
            det, info = search_arm_deep(channel, channel_edge, b, float(c['click_cx']), float(c['click_cy']), int(c['tcx']), int(c['tcy']), int(c['bs']), band, od_padded, nms_radius)  # noqa: F821
            t_s = time.perf_counter() - t0
            ranked = rank_by(det, RANK_KEY)  # noqa: F821
            assert np.array_equal(ranked[['cx', 'cy']].to_numpy(), det[['cx', 'cy']].to_numpy()), 'NMS order is not score order'
            gt_eval = gt[gt['ann_id'] != c['seed_ann_id']].reset_index(drop=True)
            out, _ = ev.bucket_detections(ranked, gt_eval, match_radius)
            xy = out[['cx', 'cy']].to_numpy(dtype=np.float64)
            mit = gt_eval[gt_eval['category_id'] == ds.MITOTIC][['cx', 'cy']].to_numpy(dtype=np.float64)
            lka = gt_eval[gt_eval['category_id'] == ds.LOOKALIKE][['cx', 'cy']].to_numpy(dtype=np.float64)
            h, w = hem.shape
            cand = pd.DataFrame(dict(file_name=f, seed_index=int(c['seed_index']), cand_rank=np.arange(len(out)), cx=out['cx'].to_numpy(), cy=out['cy'].to_numpy(), score=out['score'].to_numpy(), od=out['od'].to_numpy(), bucket=out['bucket'].to_numpy(), matched_ann_id=out['matched_ann_id'].to_numpy(), nearest_mitosis_um=nearest_um(xy, mit, mpp), nearest_lookalike_um=nearest_um(xy, lka, mpp), dist_to_click_um=np.hypot(xy[:, 0] - float(c['click_cx']), xy[:, 1] - float(c['click_cy'])) * mpp, edge_px=np.minimum(np.minimum(xy[:, 0], xy[:, 1]), np.minimum(w - 1 - xy[:, 0], h - 1 - xy[:, 1]))))
            cand['in_recorded_pool'] = cand['cand_rank'] < info['n_detections_recorded']
            cand['n_detections_recorded'] = info['n_detections_recorded']
            cand_frames.append(cand)
            click_rows.append(dict(file_name=f, seed_index=int(c['seed_index']), b=b, t_nuc=t_nuc, t_dark=t_dark, radius_px=match_radius, t_search_s=t_s, **{k: info[k] for k in ('threshold', 'median', 'mad', 'n_peaks', 'n_nms_prefix', 'n_detections_recorded', 'score_peak_100', 'score_peak_101', 'tpl_sha1', 'fill_value', 'm_keep', 'mass_weight', 'tpl_sd', 'degenerate')}))
        del channel, channel_edge, od_padded, band
        cands = pd.concat(cand_frames, ignore_index=True)
        assert (cands['edge_px'] >= PAD).all(), 'a candidate inside the 30 px band'  # noqa: F821
        e_ch = rgb2hed(rgb.astype(np.float32) / 255.0)[:, :, 1].astype(np.float32)
        gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
        wint = white_integral(rgb)  # noqa: F821
        locs = cands.drop_duplicates(['cx', 'cy'])[['cx', 'cy']].reset_index(drop=True)
        t0 = time.perf_counter()
        rows = []
        for x, y in locs.itertuples(index=False):
            r, _ = location_features(hem, e_ch, rgb, gray, float(x), float(y), mpp, t_nuc, t_dark)
            rows.append(r)
        t_feat = time.perf_counter() - t0
        feats = pd.concat([locs, pd.DataFrame(rows)], axis=1)
        wf = white_frac(wint, feats['cx'].to_numpy(), feats['cy'].to_numpy())  # noqa: F821
        feats['white_frac_gray215_61'] = np.where(feats['accepted'].to_numpy(bool), wf, np.nan)
        feats.insert(0, 'file_name', f)
        rss_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
        timing = dict(file_name=f, pid=os.getpid(), n_locations=len(locs), t_features_s=t_feat, ms_per_location=1000.0 * t_feat / max(1, len(locs)), t_roi_s=time.perf_counter() - t_roi, peak_rss_gb=rss_gb)
        return dict(clicks=click_rows, cands=cands, feats=feats, timing=timing)


# ---------------------------------------------------------------------------------------------------------------------
# montage and neighbour-check helpers (re-run the refinement at chosen locations)
# ---------------------------------------------------------------------------------------------------------------------
def montage_windows(task):
    """
        The 61 px RGB windows and the M_chr outlines at chosen locations of one ROI, for the montages.

        task (dict): path, mpp, locations (list of (cx, cy)).

        Returns dict: (cx, cy) -> dict(rgb61, m_chr51) (m_chr51 None when refused).
    """
    cv2.setNumThreads(1)
    with threadpool_limits(1):
        rgb = ds.load_roi(task['path'])
        hem = ch.to_channel(rgb, 'hematoxylin_od')
        out = {}
        for cx, cy in task['locations']:
            ref = refine_template_box(hem, cx, cy, task['mpp'])  # noqa: F821
            m_chr = ref['selection_mask'] & (ref['class_map'] == 2) if ref['accepted'] else None
            out[(cx, cy)] = dict(rgb61=tm.read_padded_patch(rgb, cx, cy, WIN61), m_chr51=m_chr)
        return out


def neighbour_check(task):
    """
        The 101 px windows of chosen candidates with Section 5.8 (d)'s component labels, for the real-data neighbour check.

        task (dict): path, mpp, t_nuc, t_dark, locations (list of (cx, cy)).

        Returns list[dict]: per location, rgb101, M_chr placed in the window, and for each of 'corrected' (t_dark) and 'spec' (t_nuc) the label image, the status map and the columns.
    """
    cv2.setNumThreads(1)
    with threadpool_limits(1):
        rgb = ds.load_roi(task['path'])
        hem = ch.to_channel(rgb, 'hematoxylin_od')
        out = []
        for cx, cy in task['locations']:
            ref = refine_template_box(hem, cx, cy, task['mpp'])  # noqa: F821
            m_chr = ref['selection_mask'] & (ref['class_map'] == 2)
            hem_101, _, (oy, ox) = window_101(hem, cx, cy)
            rgb_101, _, _ = window_101(rgb, cx, cy)
            m_chr_101 = np.zeros(hem_101.shape, bool)
            m_chr_101[oy:oy + WINDOW, ox:ox + WINDOW] = m_chr  # noqa: F821
            item = dict(cx=cx, cy=cy, rgb101=rgb_101, m_chr_101=m_chr_101)
            for name, thr in (('corrected', task['t_dark']), ('spec', task['t_nuc'])):
                lab, status = neighbour_components(hem_101, m_chr_101, task['mpp'], thr)
                item[name] = dict(lab=lab, status=status, cols=neighbour_cols(hem_101, m_chr_101, task['mpp'], thr, float(m_chr.sum()) * task['mpp'] ** 2))
            out.append(item)
        return out
