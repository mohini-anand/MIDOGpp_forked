"""
    Audit of the V2 candidate-feature run (d99e93a): reproduces its validity numbers and tests the mechanisms behind them.

    The write-up is candidate_features_a3_top300_audit.md, next to this file. Everything here is post hoc and in-sample.
    Run from anywhere with the anaconda3 interpreter; it takes about 3 minutes on 12 workers.

    Two passes:
    1. Per accepted location, the committed worker's verbatim refinement is re-run, and the audit's own measurements are
       taken: a mask-independent darkness, M_chr fragmentation, the M_dark superset check, de-noised and reworked top-hat
       columns, the resid column against the slide's own glass, and a cell-territory cytoplasm hue.
       Output: results/candidate_features_a3_top300_audit_locations.csv.gz.
    2. Statistics from candidates.csv and those measurements:
       - an independent reproduction of the committed AUCs;
       - the split of the depth-300 AUC into its pair types;
       - AUCs within strata of darkness, score, size and fragmentation, with ROI cluster-bootstrap CIs;
       - the mechanism tables and the remaining checks.
       Output: results/candidate_features_a3_top300_audit_{reproduce,decomposition,stratified,mechanisms}.csv and
       _audit_checks.json.

    Both passes read the committed candidates.csv, clicks.csv and validity tables of the run. Neither touches any of them.
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'

import json
import sys
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import scipy.ndimage as ndi
from skimage.color import rgb2hed
from skimage.measure import label
from skimage.morphology import disk, skeletonize

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _p in (str(REPO), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import candidate_features_a3_top300_worker as W  # noqa: E402
from midog_utils import channels as ch  # noqa: E402
from midog_utils import dataset as ds  # noqa: E402
from midog_utils import template_match as tm  # noqa: E402

RES = REPO / 'results'
STEM = 'candidate_features_a3_top300'
IN_CANDIDATES = RES / f'{STEM}_candidates.csv'
IN_CLICKS = RES / f'{STEM}_clicks.csv'
IN_VALIDITY = RES / f'{STEM}_validity.csv'
IN_VALIDITY_300 = RES / f'{STEM}_validity_depth300.csv'
IN_GATE_C_POOL = RES / 'candrefine_rerank_mo_pad61_select_pool.csv'
OUT_LOCATIONS = RES / f'{STEM}_audit_locations.csv.gz'
OUT_TABLES = {name: RES / f'{STEM}_audit_{name}.csv' for name in ('reproduce', 'decomposition', 'stratified', 'mechanisms')}
OUT_CHECKS = RES / f'{STEM}_audit_checks.json'
SUBSET_DIRS = {'original_14': REPO / 'images/extra_valid', 'testing_35': REPO / 'images/extra_valid/testing_set'}
N_WORKERS = 12

N_BOOT = 4000              # ROI cluster-bootstrap draws
BOOT_SEED = 20261009
N_PERM = 200               # within-click label permutations for the ROI-spread null
STRATA = 5                 # within-click quantile strata for every control
SMOOTH_SIGMA_PX = 0.7      # light de-noising of H for the top-hat variants
DARK_N_PX = 50             # dark50: mean H of the window's darkest 50 px
DISC_SEARCH_UM, DISC_R_UM = 2.0, 1.0   # disc_h: mean H in a 1 um disc at the smoothed-H maximum within 2 um of the centre
GLASS_MIN_PX = 1000        # fewer glass pixels than this and I0 falls back to 255
FAILING9 = ('211', '242', '524', '546', '462', '527', '410', '245', '289')
FEATURES = list(W.FEATURES)
EXTRAS = ['chr_frag', 'skel_branch_per_um', 'noise_h', 'lap_sd', 'tophat_edge_smooth', 'tophat_fix', 'resid_i0', 'resid_core_only', 'resid_vnorm', 'eh_territory', 'dark50', 'disc_h']


# ---------------------------------------------------------------------------------------------------------------------
# pass 1: per-location measurements
# ---------------------------------------------------------------------------------------------------------------------
def roi_i0(rgb):
    """
        The slide's glass intensity per channel, for optical density against the real background instead of 255.

        rgb (np.ndarray): the ROI's uint8 RGB.

        Returns tuple[np.ndarray, int]: the per-channel median over glass pixels (gray > 215 and every channel > 200), or 255s when fewer than GLASS_MIN_PX exist; and the glass pixel count.
    """
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    m = (gray > W.GLASS_GRAY) & np.all(rgb > 200, axis=2)
    if m.sum() < GLASS_MIN_PX:
        return np.full(3, 255.0), int(m.sum())
    return np.median(rgb[m].reshape(-1, 3).astype(np.float64), axis=0), int(m.sum())


def resid_of(od, support):
    """
        chr_resid_frac's statistic on one support: the mean OD's component along the H-E plane normal, over its norm.

        od (np.ndarray): per-pixel OD, shape (h, w, 3).
        support (np.ndarray): the pixels to average.

        Returns tuple[float, float]: the statistic and the mean OD's norm (both NaN below NS_MIN_PX pixels).
    """
    if support.sum() < W.NS_MIN_PX:
        return np.nan, np.nan
    v = od[support].mean(axis=0)
    return float(v @ W.RESID_NORMAL / np.linalg.norm(v)), float(np.linalg.norm(v))


def location_extras(hem, rgb, e_ch, gray, cx, cy, mpp, t_dark, i0):
    """
        The audit's measurements at one candidate, on the refinement the run used.

        hem (np.ndarray): the ROI's unclipped hematoxylin_od.
        rgb (np.ndarray): the ROI's uint8 RGB.
        e_ch (np.ndarray): the ROI's rgb2hed E channel.
        gray (np.ndarray): the ROI's uint8 gray.
        cx (float): detection x.
        cy (float): detection y.
        mpp (float): microns per pixel.
        t_dark (float): the ROI's upper 3-class multi-Otsu threshold.
        i0 (np.ndarray): the slide's glass intensity per channel.

        Returns dict or None: the measurements, or None when the run refused the candidate.
    """
    ref = W.refine_template_box(hem, cx, cy, mpp)
    if not ref['accepted']:
        return None
    m_chr = ref['selection_mask'] & (ref['class_map'] == 2)
    if not m_chr.any():
        return None
    h = np.asarray(ref['hem_patch'], np.float32)
    hd = h.astype(np.float64)
    rgb_w = tm.read_padded_patch(rgb, cx, cy, W.WINDOW)
    g_w = tm.read_padded_patch(gray, cx, cy, W.WINDOW)
    e_w = tm.read_padded_patch(e_ch, cx, cy, W.WINDOW)
    out = dict(align_err=float(np.abs(rgb2hed(rgb_w.astype(np.float32) / 255.0)[:, :, 0] - h).max()), t_win_upper=float(ref['thresholds'][1]))
    piece = W.largest_piece(m_chr)
    filled = ndi.binary_fill_holes(piece)
    # M_dark exactly as the worker builds it, and whether it contains the filled M_chr piece
    comp = label(h > t_dark, connectivity=1)
    own = np.unique(comp[m_chr])
    own = own[own > 0]
    mdark = ndi.binary_fill_holes(comp == int(own[int(np.argmax(np.bincount(comp.ravel())[own]))])) if len(own) else np.zeros_like(m_chr)
    out['filled_in_mdark'] = float((filled & mdark).sum() / max(filled.sum(), 1))
    # fragmentation of M_chr, and the skeleton's branch density on its largest piece
    sizes = np.bincount(label(m_chr, connectivity=1).ravel())[1:]
    out['mchr_n4'] = int((sizes > 0).sum())
    out['chr_frag'] = float(1.0 - sizes.max() / sizes.sum())
    sk = skeletonize(piece)
    nbr = ndi.convolve(sk.astype(np.uint8), np.ones((3, 3), np.uint8), mode='constant') - sk
    out['skel_branch_per_um'] = float(((nbr >= 3) & sk).sum() / (sk.sum() * mpp)) if sk.any() else np.nan
    # pixel noise: robust SD of H about its 3 x 3 median, over the window
    dev = hd - ndi.median_filter(hd, 3)
    out['noise_h'] = float(1.4826 * np.median(np.abs(dev - np.median(dev))))
    # h_chr_lap's numerator alone
    eroded = ndi.binary_erosion(m_chr)
    lap = cv2.Laplacian(np.ascontiguousarray(h), cv2.CV_32F, ksize=1) / mpp ** 2
    out['lap_sd'] = float(lap[eroded].astype(np.float64).std()) if eroded.sum() >= W.NS_MIN_PX else np.nan
    # top-hat: the committed construction (a check), the same on lightly smoothed H, and a reworked version
    se = disk(max(1, int(round(W.TOPHAT_R_UM / mpp)))).astype(np.uint8)
    hs = cv2.GaussianBlur(h, (0, 0), SMOOTH_SIGMA_PX)
    t_raw = (h - cv2.morphologyEx(h, cv2.MORPH_OPEN, se, borderType=cv2.BORDER_REPLICATE)).astype(np.float64)
    t_s = (hs - cv2.morphologyEx(hs, cv2.MORPH_OPEN, se, borderType=cv2.BORDER_REPLICATE)).astype(np.float64)
    r_edge = (ndi.distance_transform_edt(~m_chr) * mpp <= W.TOPHAT_RING_UM) & ~(ndi.distance_transform_edt(m_chr) * mpp > W.TOPHAT_CORE_UM)
    out['tophat_check'] = float(t_raw[r_edge].sum() / hd[m_chr].sum()) if m_chr.sum() >= W.MIN_PIECE_PX else np.nan
    out['tophat_edge_smooth'] = float(t_s[r_edge].sum() / hs.astype(np.float64)[m_chr].sum()) if m_chr.sum() >= W.MIN_PIECE_PX else np.nan
    out['tophat_fix'] = np.nan
    if piece.sum() >= W.MIN_PIECE_PX:
        reg = (ndi.distance_transform_edt(~filled) * mpp <= W.TOPHAT_RING_UM) & ~(ndi.distance_transform_edt(filled) * mpp > W.TOPHAT_CORE_UM)
        per_um = cv2.arcLength(W.outline_of(filled), True) * mpp
        out['tophat_fix'] = float(t_s[reg].sum() * mpp ** 2 / (per_um * hd[filled].mean()))
    # chr_resid_frac: the committed construction (a check), against the slide's glass, on M_chr alone, and its support's OD norm
    m_res = ndi.binary_dilation(m_chr, structure=disk(W.RESID_DILATE_PX).astype(bool))
    in_range = np.all((rgb_w >= W.NS_LO) & (rgb_w <= W.NS_HI), axis=2)
    od255 = -np.log10(np.maximum(rgb_w.astype(np.float64), 1.0) / 255.0)
    od_i0 = -np.log10(np.maximum(rgb_w.astype(np.float64), 1.0) / i0[None, None, :])
    out['resid_check'], out['resid_vnorm'] = resid_of(od255, m_res & in_range)
    out['resid_i0'], _ = resid_of(od_i0, m_res & in_range)
    out['resid_core_only'], _ = resid_of(od255, m_chr & in_range)
    # cytoplasm hue on the candidate's own cell territory: ring pixels whose nearest H > t_dark component is its own
    out['eh_territory'] = np.nan
    if comp.max() > 0 and len(own):
        _, (iy, ix) = ndi.distance_transform_edt(comp == 0, return_indices=True)
        rings = W.ring_masks(m_chr, mpp)
        clipped = np.any((rgb_w == 0) | (rgb_w == 255), axis=2)
        cyto = (rings[1] | rings[2]) & np.isin(comp[iy, ix], own) & (g_w <= W.GLASS_GRAY) & (h <= t_dark) & ~clipped
        if cyto.sum() >= W.NS_MIN_PX:
            out['eh_territory'] = float(np.degrees(np.arctan2(e_w[cyto].astype(np.float64).mean(), hd[cyto].mean())))
    # darkness covariates the refinement mask does not set
    out['dark50'] = float(np.sort(hd.ravel())[-DARK_N_PX:].mean())
    yy, xx = np.mgrid[0:W.WINDOW, 0:W.WINDOW] - W.WINDOW // 2
    s = cv2.GaussianBlur(hd, (0, 0), 1.0)
    iy0, ix0 = np.unravel_index(np.argmax(np.where(np.hypot(xx, yy) * mpp <= DISC_SEARCH_UM, s, -np.inf)), s.shape)
    out['disc_h'] = float(hd[np.hypot(xx + W.WINDOW // 2 - ix0, yy + W.WINDOW // 2 - iy0) * mpp <= DISC_R_UM].mean())
    return out


def measure_roi(task):
    """
        Pass 1 on one ROI: the audit's measurements at every accepted location.

        task (dict): file_name, path, mpp, t_dark, locs (list of (cx, cy)).

        Returns pd.DataFrame: one row per accepted location.
    """
    cv2.setNumThreads(1)
    rgb = ds.load_roi(task['path'])
    hem = ch.to_channel(rgb, 'hematoxylin_od')
    e_ch = rgb2hed(rgb.astype(np.float32) / 255.0)[:, :, 1].astype(np.float32)
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    i0, n_glass = roi_i0(rgb)
    rows = []
    for cx, cy in task['locs']:
        r = location_extras(hem, rgb, e_ch, gray, cx, cy, task['mpp'], task['t_dark'], i0)
        if r is not None:
            rows.append(dict(file_name=task['file_name'], cx=cx, cy=cy, t_dark=task['t_dark'], i0_r=i0[0], i0_g=i0[1], i0_b=i0[2], n_glass=n_glass, **r))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------------------------------
# pass 2: statistics
# ---------------------------------------------------------------------------------------------------------------------
def brute_auc(df, col):
    """
        The pooled within-click AUC by explicit pair counting, independent of the notebook's rank-based code.

        df (pd.DataFrame): rows with click, is_mitosis and the column.
        col (str): the column.

        Returns float: sum U over sum n1 n0, ties counted one half.
    """
    u = n = 0.0
    for _, g in df[df[col].notna()].groupby('click'):
        p = g.loc[g['is_mitosis'], col].to_numpy(float)
        q = np.sort(g.loc[~g['is_mitosis'], col].to_numpy(float))
        if len(p) and len(q):
            lt, le = np.searchsorted(q, p, 'left'), np.searchsorted(q, p, 'right')
            u += lt.sum() + 0.5 * (le - lt).sum()
            n += len(p) * len(q)
    return u / n if n else np.nan


def roi_u_n(df, col, keys, rois):
    """
        Per-ROI totals of Mann-Whitney U and n1 n0 over the groups in keys, NaN rows dropped, groups with both classes only.

        df (pd.DataFrame): rows with file_name, is_mitosis, the keys and the column.
        col (str): the column.
        keys (list): the grouping columns, click first.
        rois (np.ndarray): the ROI order of the output.

        Returns tuple[np.ndarray, np.ndarray]: U and n1 n0 per ROI.
    """
    ok = df[col].notna()
    d = df.loc[ok, keys + ['file_name']].assign(v=df.loc[ok, col].astype(float), pos=df.loc[ok, 'is_mitosis'].astype(int))
    d['rp'] = d.groupby(keys)['v'].rank(method='average') * d['pos']
    g = d.groupby(keys).agg(file_name=('file_name', 'first'), s1=('rp', 'sum'), n1=('pos', 'sum'), n=('pos', 'size'))
    g['n0'] = g['n'] - g['n1']
    g = g[(g['n1'] > 0) & (g['n0'] > 0)]
    g['U'] = g['s1'] - g['n1'] * (g['n1'] + 1) / 2.0
    g['N'] = g['n1'] * g['n0']
    t = g.groupby('file_name')[['U', 'N']].sum().reindex(rois, fill_value=0.0)
    return t['U'].to_numpy(), t['N'].to_numpy()


def boot(u, n, wts):
    """
        A pooled AUC and its ROI cluster-bootstrap distribution.

        u (np.ndarray): U per ROI.
        n (np.ndarray): n1 n0 per ROI.
        wts (np.ndarray): (N_BOOT, n_rois) multinomial ROI weights.

        Returns tuple[float, np.ndarray]: the AUC and the bootstrap AUCs (NaN where a draw has no pair).
    """
    with np.errstate(invalid='ignore', divide='ignore'):
        return float(u.sum() / n.sum()) if n.sum() else np.nan, (wts @ u) / (wts @ n)


def ci(bs):
    """
        The 95 % percentile interval of a bootstrap distribution.

        bs (np.ndarray): the bootstrap values.

        Returns tuple[float, float]: the 2.5th and 97.5th percentiles, NaNs ignored.
    """
    return float(np.nanpercentile(bs, 2.5)), float(np.nanpercentile(bs, 97.5))


def add_strata(df, col, name):
    """
        Within-click quantile strata of a column, formed once per population before any NaN is dropped.

        df (pd.DataFrame): the population, with click.
        col (str): the stratifying column.
        name (str): the stratum column to add.

        Returns None.
    """
    df[name] = df.groupby('click')[col].transform(lambda v: pd.qcut(v.rank(method='first'), STRATA, labels=False) if v.notna().sum() >= STRATA else np.nan)


def spearman_within_click(df, a, b, min_rows=10):
    """
        The mean over clicks of the within-click Spearman correlation of two columns.

        df (pd.DataFrame): rows with click.
        a (str): one column.
        b (str): the other.
        min_rows (int): clicks with fewer complete rows are skipped.

        Returns float: the mean correlation.
    """
    rs = [g[a].rank().corr(g[b].rank()) for _, g in df[['click', a, b]].dropna().groupby('click') if len(g) >= min_rows]
    return float(np.nanmean(rs))


def run_stats(C, L):
    """
        Pass 2: every audit table and check, from the candidate table and the pass-1 measurements.

        C (pd.DataFrame): candidates.csv.
        L (pd.DataFrame): the pass-1 measurements, one row per accepted location.

        Returns tuple[dict, dict]: the tables by name, and the checks.
    """
    A = C[C['accepted']].merge(L, on=['file_name', 'cx', 'cy'], how='left', validate='many_to_one')
    A['click'] = A['file_name'] + '|' + A['seed_index'].astype(str)
    assert A['dark50'].notna().all(), 'a pass-1 measurement is missing for an accepted row'
    rois = np.array(sorted(A['file_name'].unique()))
    wts = np.random.default_rng(BOOT_SEED).multinomial(len(rois), np.full(len(rois), 1.0 / len(rois)), size=N_BOOT).astype(float)
    pops = {'pool': A[A['in_recorded_pool']].copy(), 'beyond': A[~A['in_recorded_pool']].copy(), 'd300': A.copy()}
    controls = {'none': None, 'od': 'od', 'h_chr_mean': 'h_chr_mean', 'dark50': 'dark50', 'disc_h': 'disc_h', 'score': 'score', 'chr_area_um2': 'chr_area_um2', 'chr_frag': 'chr_frag'}
    for P in pops.values():
        for name, col in controls.items():
            if col is not None:
                add_strata(P, col, f'q_{name}')
    checks = {}

    # 1. independent reproduction of the committed AUCs
    V = pd.read_csv(IN_VALIDITY, float_precision='round_trip').set_index('column')
    V3 = pd.read_csv(IN_VALIDITY_300, float_precision='round_trip').set_index('column')
    rep = pd.DataFrame([dict(column=c, pool=brute_auc(pops['pool'], c), pool_committed=V.loc[c, 'auc_within_click'], d300=brute_auc(pops['d300'], c), d300_committed=V3.loc[c, 'auc_within_click'], beyond=brute_auc(pops['beyond'], c), beyond_committed=V3.loc[c, 'auc_wc_beyond_pool']) for c in FEATURES])
    rep['max_abs_diff'] = np.abs(rep[['pool', 'd300', 'beyond']].to_numpy() - rep[['pool_committed', 'd300_committed', 'beyond_committed']].to_numpy()).max(axis=1)
    checks['reproduce_max_abs_diff'] = float(rep['max_abs_diff'].max())
    checks['resid_check_max_abs_diff'] = float(np.nanmax(np.abs(A['resid_check'] - A['chr_resid_frac'])))
    checks['tophat_check_max_abs_diff'] = float(np.nanmax(np.abs(A['tophat_check'] - A['chr_tophat_edge'])))
    checks['rgb_h_window_align_max_abs_diff'] = float(A['align_err'].max())

    # 2. the depth-300 AUC split by pair type (pool or beyond for the mitosis, pool or beyond for the other)
    dec = []
    for c in FEATURES:
        row = dict(column=c, d300=V3.loc[c, 'auc_within_click'])
        parts = {}
        for name, mit_pool, oth_pool in (('pool_pool', True, True), ('beyond_beyond', False, False), ('poolmit_beyondoth', True, False), ('beyondmit_pooloth', False, True)):
            u = n = 0.0
            for _, g in A[A[c].notna()].groupby('click'):
                p = g.loc[g['is_mitosis'] & (g['in_recorded_pool'] == mit_pool), c].to_numpy(float)
                q = np.sort(g.loc[~g['is_mitosis'] & (g['in_recorded_pool'] == oth_pool), c].to_numpy(float))
                if len(p) and len(q):
                    lt, le = np.searchsorted(q, p, 'left'), np.searchsorted(q, p, 'right')
                    u += lt.sum() + 0.5 * (le - lt).sum()
                    n += len(p) * len(q)
            parts[name] = (u / n if n else np.nan, n)
        total = sum(v[1] for v in parts.values())
        for name, (a, n) in parts.items():
            row[f'auc_{name}'], row[f'pair_share_{name}'] = a, n / total
        dec.append(row)
    dec = pd.DataFrame(dec)

    # 3. AUCs within strata, with ROI cluster-bootstrap CIs
    strat, bs_keep = [], {}
    for pname, P in pops.items():
        for c in FEATURES + EXTRAS:
            for name in controls:
                keys = ['click'] if name == 'none' else ['click', f'q_{name}']
                u, n = roi_u_n(P, c, keys, rois)
                pt, bs = boot(u, n, wts)
                lo, hi = ci(bs)
                bs_keep[(pname, c, name)] = bs
                strat.append(dict(population=pname, column=c, control=name, auc=pt, ci_lo=lo, ci_hi=hi))
    strat = pd.DataFrame(strat)
    for c in ('chr_radial_hf', 'chr_tophat_edge'):
        checks[f'{c}_pool_minus_beyond_ci'] = ci(bs_keep[('pool', c, 'none')] - bs_keep[('beyond', c, 'none')])

    # 4. mechanism tables
    mech = []
    for a, b in (('chr_tophat_edge', 'chr_frag'), ('chr_tophat_edge', 'chr_hole_frac'), ('chr_tophat_edge', 'skel_branch_per_um'), ('chr_tophat_edge', 'noise_over_h'), ('chr_tophat_edge', 'chr_area_um2'), ('chr_tophat_edge', 'h_chr_mean'), ('chr_radial_hf', 'chr_frag'), ('chr_radial_hf', 'skel_branch_per_um'), ('chr_radial_hf', 'h_chr_mean'), ('h_chr_lap', 'noise_over_h'), ('h_chr_lap', 'chr_area_um2'), ('h_chr_lap', 'h_chr_mean'), ('chr_area_um2', 'h_chr_mean'), ('chr_area_um2', 'dark50'), ('chr_area_um2', 'disc_h'), ('chr_area_um2', 'od'), ('dark50', 'h_chr_mean'), ('disc_h', 'h_chr_mean'), ('chr_resid_frac', 'resid_vnorm'), ('chr_resid_frac', 'chr_ns_share')):
        D = A.assign(noise_over_h=A['noise_h'] / A['h_chr_mean'])
        mech.append(dict(kind='spearman_within_click_d300', a=a, b=b, value=spearman_within_click(D, a, b)))
    A['mdark_not_superset'] = A['filled_in_mdark'] < 0.99
    A['win_over_tdark'] = A['t_win_upper'] / A['t_dark']
    for lab, g in A.groupby('is_mitosis'):
        who = 'mitosis' if lab else 'other'
        mech.append(dict(kind='mdark_not_superset_share', a=who, b='', value=float(g['mdark_not_superset'].mean())))
        for q in (0.05, 0.5, 0.95):
            mech.append(dict(kind=f'window_upper_otsu_over_t_dark_q{int(q * 100):02d}', a=who, b='', value=float(g['win_over_tdark'].quantile(q))))
        for c in ('mchr_n4', 'chr_frag', 'skel_branch_per_um', 'resid_vnorm'):
            mech.append(dict(kind=f'median_{c}', a=who, b='', value=float(g[c].median())))
    for sc, g in A.groupby('scanner'):
        per_roi = g.groupby('file_name')[['i0_r', 'i0_g', 'i0_b']].first().median()
        for chn in ('i0_r', 'i0_g', 'i0_b'):
            mech.append(dict(kind='glass_i0_median_over_rois', a=sc, b=chn, value=float(per_roi[chn])))
    mech = pd.DataFrame(mech)

    # 5. remaining checks
    Q = {p: P[P['dark_touches_window'] == False] for p, P in pops.items()}  # noqa: E712
    for p, P in Q.items():
        ua, na = roi_u_n(P, 'chr_concavities_per_um', ['click'], rois)
        ub, nb = roi_u_n(P, 'dark_concavities_per_um', ['click'], rois)
        a, ba = boot(ua, na, wts)
        b, bb = boot(ub, nb, wts)
        checks[f'uncut_{p}_dark_minus_chr_concavities'] = dict(chr=a, dark=b, lead=b - a, ci=ci(bb - ba))
    rng = np.random.default_rng(BOOT_SEED)
    for c in ('e_h_angle_ring23_deg', 'eh_territory'):
        P = pops['pool']
        u, n = roi_u_n(P, c, ['click'], rois)
        ra = u[n > 0] / n[n > 0]
        g9 = np.isin([r.split('.')[0] for r in rois], FAILING9)
        null_sd = []
        for _ in range(N_PERM):
            Pp = P[['click', 'file_name', 'is_mitosis', c]].copy()
            Pp['is_mitosis'] = Pp.groupby('click')['is_mitosis'].transform(lambda v: rng.permutation(v.to_numpy()))
            up, np_ = roi_u_n(Pp, c, ['click'], rois)
            null_sd.append(np.std(up[np_ > 0] / np_[np_ > 0]))
        checks[f'roi_spread_pool_{c}'] = dict(median=float(np.median(ra)), q25=float(np.percentile(ra, 25)), q75=float(np.percentile(ra, 75)), min=float(ra.min()), max=float(ra.max()), sd=float(ra.std()), null_sd_p95_within_click_permutation=float(np.percentile(null_sd, 95)), share_above_half=float((ra > 0.5).mean()), failing9=float(u[g9].sum() / n[g9].sum()))
    old = pd.read_csv(IN_GATE_C_POOL, usecols=['file_name', 'cx', 'cy']).drop_duplicates()
    locs = C.groupby(['file_name', 'cx', 'cy'])['in_recorded_pool'].max().reset_index()
    shared = locs.merge(old, on=['file_name', 'cx', 'cy'])
    checks['gate_c'] = dict(shared_locations=len(shared), shared_in_recorded_pool=int(shared['in_recorded_pool'].sum()), recorded_pool_locations=int(locs['in_recorded_pool'].sum()), share_of_recorded_pool=float(shared['in_recorded_pool'].sum() / locs['in_recorded_pool'].sum()))
    checks['hole_frac_gt_0p3'] = dict(rows=int((A['chr_hole_frac'] > 0.3).sum()), mitoses=int(((A['chr_hole_frac'] > 0.3) & A['is_mitosis']).sum()))
    D = A.assign(has_partner=A['pair_sep_um'].notna().astype(float), pair_sep_none_far=A['pair_sep_um'].fillna(1e9))
    checks['pair_missingness_auc'] = {f'{c}_{p}': brute_auc(D[D['in_recorded_pool']] if p == 'pool' else D, c) for c in ('has_partner', 'pair_sep_none_far') for p in ('pool', 'd300')}
    return dict(reproduce=rep, decomposition=dec, stratified=strat, mechanisms=mech), checks


def main():
    """
        Run pass 1 in spawn workers, then pass 2, and write every audit output.

        Returns None.
    """
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    C = pd.read_csv(IN_CANDIDATES, float_precision='round_trip', low_memory=False)
    K = pd.read_csv(IN_CLICKS, float_precision='round_trip').groupby('file_name')['t_dark'].first()
    acc = C[C['accepted']].drop_duplicates(['file_name', 'cx', 'cy'])
    tasks = [dict(file_name=f, path=str(SUBSET_DIRS[g['subset'].iat[0]] / f), mpp=float(g['mpp'].iat[0]), t_dark=float(K[f]), locs=list(map(tuple, g[['cx', 'cy']].to_numpy(float)))) for f, g in acc.groupby('file_name')]
    with ProcessPoolExecutor(max_workers=N_WORKERS, mp_context=get_context('spawn')) as pool:
        L = pd.concat(list(pool.map(measure_roi, tasks)), ignore_index=True)
    L.to_csv(OUT_LOCATIONS, index=False)
    print(f'pass 1: {len(L)} accepted locations -> {OUT_LOCATIONS.name}', flush=True)
    tables, checks = run_stats(C, L)
    for name, df in tables.items():
        df.to_csv(OUT_TABLES[name], index=False)
    with open(OUT_CHECKS, 'w') as fh:
        json.dump(checks, fh, sort_keys=True, indent=1)
        fh.write('\n')
    print('pass 2:', ', '.join(p.name for p in list(OUT_TABLES.values()) + [OUT_CHECKS]))
    print(json.dumps(checks, sort_keys=True, indent=1))


if __name__ == '__main__':
    main()
