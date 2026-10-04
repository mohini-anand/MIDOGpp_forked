"""
    Independent audit of production_hematoxylin_only/window_sweep_51_61_71_49roi_6seed.ipynb.

    Every number in Research Logs/2026-09-29-window-sweep-51-61-71-49roi-6seed-audit.md is printed by this
    script, and every table it cites is written to results/window_sweep_51_61_71_49roi_6seed_audit_*.csv.

    Tier A re-derives the notebook's statistics from its persisted CSVs and from databases/MIDOG++.json with plain
    numpy / pandas / scipy (no midog_utils import anywhere in this file). Tier B re-derives templates and searches
    from the ROI pixels with cv2 / skimage / tifffile directly, in a multiprocessing pass over the 49 ROIs.

    Run: /Users/mohinianand/anaconda3/bin/python3 window_sweep_51_61_71_49roi_6seed_audit.py
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time
from multiprocessing import get_context

import numpy as np
import pandas as pd
from scipy import signal as spsig
from scipy import stats as sps

REPO = os.path.dirname(os.path.abspath(__file__))
NB_PATH = f'{REPO}/production_hematoxylin_only/window_sweep_51_61_71_49roi_6seed.ipynb'
WS = f'{REPO}/results/precision_at_k_49roi_6seed_windowsweep_postD11'
FG = f'{REPO}/results/precision_at_k_49roi_6seed_fillgeom_postD11'
OUT = f'{REPO}/results/window_sweep_51_61_71_49roi_6seed_audit'
DB = f'{REPO}/databases/MIDOG++.json'
IMG_DIRS = {'original_14': f'{REPO}/images/extra_valid', 'testing_35': f'{REPO}/images/extra_valid/testing_set'}

ARMS = ('default_51', 'fill_51', 'default_61', 'fill_61', 'default_71', 'fill_71', 'fill_71_pinned')
ARM_W = {a: int(a.split('_')[1]) for a in ARMS}
AXES = ('chromatin_od', 'tm_score')
AXIS_KEY = {'chromatin_od': 'od', 'tm_score': 'score'}
SEEDS = (0, 1, 2, 3, 4, 5)
BOOT_SEED, N_BOOT = 20260917, 10_000
TP = 'human_correct_label'
NEAR_PX, EDGE_PX, SELF_PX = 60, 35, 2
MITOTIC, LOOKALIKE = 1, 2
CONTRASTS = {
    'fill_61 - fill_51': {'fill_61': 1, 'fill_51': -1},
    'fill_71 - fill_51': {'fill_71': 1, 'fill_51': -1},
    '(fill_61 - fill_51) - (default_61 - default_51)': {'fill_61': 1, 'fill_51': -1, 'default_61': -1, 'default_51': 1},
    '(fill_71 - fill_51) - (default_71 - default_51)': {'fill_71': 1, 'fill_51': -1, 'default_71': -1, 'default_51': 1},
    'default_61 - default_51': {'default_61': 1, 'default_51': -1},
    'default_71 - default_51': {'default_71': 1, 'default_51': -1},
    'fill_71_pinned - fill_71': {'fill_71_pinned': 1, 'fill_71': -1},
    'fill_71_pinned - fill_51': {'fill_71_pinned': 1, 'fill_51': -1},
    'fill_51 - default_51': {'fill_51': 1, 'default_51': -1},
    'fill_61 - default_61': {'fill_61': 1, 'default_61': -1},
    'fill_71 - default_71': {'fill_71': 1, 'default_71': -1},
}
EXTRA_CONTRASTS = {'fill_71 - fill_61': {'fill_71': 1, 'fill_61': -1}}
TESTED = {'fill_61 - fill_51': 'window', 'fill_71 - fill_51': 'window', '(fill_61 - fill_51) - (default_61 - default_51)': 'fill_vs_raw', '(fill_71 - fill_51) - (default_71 - default_51)': 'fill_vs_raw'}
METRIC_K = {'tp_at_10': 10, 'tp_at_10_excl_border': 10, 'tp_at_20': 20, 'tp_at_30': 30}
TIER_B_CLICKS_ALL7 = [('301.tiff', 2), ('220.tiff', 2), ('549.tiff', 3), ('291.tiff', 0)]
N_WORKERS = 4

LEDGER = []   # every comparison: check, n_compared, n_diverge, note


def ledger(check, n_compared, n_diverge, note=''):
    """
        Record one comparison block in the audit ledger and echo it.

        check (str): what was compared.
        n_compared (int): values compared.
        n_diverge (int): values that diverged.
        note (str): tolerance or caveat.

        Returns None.
    """
    LEDGER.append(dict(check=check, n_compared=int(n_compared), n_diverge=int(n_diverge), note=note))
    print(f'  [ledger] {check}: {n_compared} compared, {n_diverge} diverge {("(" + note + ")") if note else ""}')


def section(title):
    """
        Print a section banner.

        title (str): the banner text.

        Returns None.
    """
    print('\n' + '=' * 110 + f'\n{title}\n' + '=' * 110)


# ---------------------------------------------------------------------------------------------------------------
# independent primitives (no midog_utils)
# ---------------------------------------------------------------------------------------------------------------

def load_gt():
    """
        Parse MIDOG++.json directly: images and annotation centres, both categories.

        Returns tuple[pd.DataFrame, pd.DataFrame]: images (file_name, width, height, tumor_type), annotations (ann_id, file_name, cx, cy, category_id).
    """
    raw = json.load(open(DB))
    alias = {'canine lymphoma': 'canine lymphosarcoma'}
    ims = pd.DataFrame([dict(image_id=i['id'], file_name=i['file_name'], width=i['width'], height=i['height'], tumor_type=alias.get(i['tumor_type'], i['tumor_type'])) for i in raw['images']])
    name = dict(zip(ims['image_id'], ims['file_name']))
    anns = pd.DataFrame([dict(ann_id=a['id'], file_name=name[a['image_id']], cx=(a['bbox'][0] + a['bbox'][2]) / 2.0, cy=(a['bbox'][1] + a['bbox'][3]) / 2.0, category_id=a['category_id']) for a in raw['annotations']])
    return ims, anns


def greedy_match(det_xy, gt_xy, radius):
    """
        Greedy best-detection-first one-to-one matching by centre distance (distance <= radius), nearest free GT wins.

        det_xy (np.ndarray): (N, 2) ranked detections.
        gt_xy (np.ndarray): (M, 2) ground truth.
        radius (float): match radius in px.

        Returns np.ndarray: per detection the GT index claimed, -1 if none.
    """
    out = np.full(len(det_xy), -1, dtype=int)
    if len(det_xy) == 0 or len(gt_xy) == 0:
        return out
    taken = np.zeros(len(gt_xy), bool)
    dist = np.hypot(det_xy[:, None, 0] - gt_xy[None, :, 0], det_xy[:, None, 1] - gt_xy[None, :, 1])
    for i in range(len(det_xy)):
        ok = (dist[i] <= radius) & ~taken
        if ok.any():
            cand = np.where(ok)[0]
            g = cand[np.argmin(dist[i, cand])]
            out[i] = g
            taken[g] = True
    return out


def bucket(det_xy, gt, radius):
    """
        Bucket ranked detections against a GT frame: human_correct_label / human_rejected_label / non_human_findings.

        det_xy (np.ndarray): (N, 2) ranked detections.
        gt (pd.DataFrame): ground truth with cx, cy, category_id, ann_id.
        radius (float): match radius in px.

        Returns tuple[np.ndarray, np.ndarray]: bucket labels and matched ann_id (-1 if none).
    """
    g = greedy_match(np.asarray(det_xy, float), gt[['cx', 'cy']].to_numpy(float), radius)
    cat, ann = gt['category_id'].to_numpy(), gt['ann_id'].to_numpy()
    labels = np.array([('non_human_findings' if j < 0 else (TP if cat[j] == MITOTIC else 'human_rejected_label')) for j in g], dtype=object)
    return labels, np.array([-1 if j < 0 else int(ann[j]) for j in g])


def signflip_p(d):
    """
        Exact two-sided sign-flip p over integer per-unit deltas, by integer polynomial convolution.

        d (array-like): one integer delta per unit; zeros carry no sign.

        Returns tuple[float, int]: p-value and the number of nonzero units.
    """
    d = np.asarray(d).astype(np.int64)
    nz = np.abs(d[d != 0])
    if len(nz) == 0:
        return 1.0, 0
    dist = np.array([1], dtype=np.int64)
    for v in nz:
        k = np.zeros(2 * v + 1, dtype=np.int64)
        k[0] = k[-1] = 1
        dist = np.convolve(dist, k)
    T = int(nz.sum())
    sums = np.arange(-T, T + 1)
    return float(int(dist[np.abs(sums) >= abs(int(d.sum()))].sum()) / 2 ** len(nz)), len(nz)


def holm(p):
    """
        Holm step-down adjusted p-values.

        p (array-like): raw p-values of one family.

        Returns np.ndarray: adjusted p-values in input order.
    """
    p = np.asarray(p, float)
    m = len(p)
    order = sorted(range(m), key=lambda i: (p[i], i))
    adj, run = np.empty(m), 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (m - rank) * p[i]))
        adj[i] = run
    return adj


def boot_idx(n, stream=0):
    """
        Bootstrap ROI resample indices, drawn from the same numpy stream the notebook documents (needed to reproduce its CI).

        n (int): number of units.
        stream (int): RNG stream under BOOT_SEED.

        Returns np.ndarray: (N_BOOT, n) indices.
    """
    return np.random.default_rng([BOOT_SEED, stream]).integers(n, size=(N_BOOT, n))


def contrast_stats(d, k, idx, n_clicks=6):
    """
        Every statistic of one per-ROI integer delta vector, recomputed independently.

        d (np.ndarray): integer per-ROI deltas in bootstrap ROI order.
        k (int): the budget K.
        idx (np.ndarray): bootstrap indices.
        n_clicks (int): clicks per ROI.

        Returns dict: the statistic set, column-named like the notebook's delta_stats.
    """
    d = np.asarray(d, np.int64)
    n = len(d)
    den = n_clicks * n * k
    pooled = d.sum() / den
    w, l, t = int((d > 0).sum()), int((d < 0).sum()), int((d == 0).sum())
    p, nnz = signflip_p(d)
    sums = d[idx].sum(axis=1)
    lo, hi = np.percentile(sums / den, [2.5, 97.5])
    x = d / (n_clicks * k)
    half = sps.t.ppf(0.975, n - 1) * x.std(ddof=1) / math.sqrt(n)
    same = w if pooled > 0 else (l if pooled < 0 else 0)
    ci_ex = bool(lo > 0 or hi < 0)
    return dict(n_roi=n, pooled_delta_precision=pooled, wins=w, losses=l, ties=t, n_nonzero=nnz, exact_p=p, min_attainable_p=2.0 ** (1 - nnz) if nnz else 1.0, ci_low=float(lo), ci_high=float(hi), boot_share_ge_0=float((sums >= 0).mean()), boot_share_le_0=float((sums <= 0).mean()), t_ci_low=float(x.mean() - half), t_ci_high=float(x.mean() + half), ci_excludes_0=ci_ex, t_ci_excludes_0=bool(x.mean() - half > 0 or x.mean() + half < 0), n_same_direction=same, majority_of_all_rois=bool(same > n / 2), d5_bar_met=bool(ci_ex and same > n / 2))


def spearman(x, y):
    """
        Spearman rho as the Pearson correlation of average ranks, with scipy's t-based p for the reported p only.

        x (array-like): first variable.
        y (array-like): second variable.

        Returns tuple[float, float]: rho and two-sided p; NaN when undefined.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan, np.nan
    rx, ry = sps.rankdata(x), sps.rankdata(y)
    rho = float(np.corrcoef(rx, ry)[0, 1])
    tt = rho * math.sqrt((len(x) - 2) / max(1e-300, 1 - rho ** 2))
    return rho, float(2 * sps.t.sf(abs(tt), len(x) - 2))


# ---------------------------------------------------------------------------------------------------------------
# Tier B worker: templates and searches from pixels, cv2 / skimage / tifffile only
# ---------------------------------------------------------------------------------------------------------------

def _roi_mpp(path):
    """
        Microns per pixel from the TIFF resolution tags, honouring ResolutionUnit.

        path (str): ROI TIFF path.

        Returns float: microns per pixel.
    """
    import tifffile
    with tifffile.TiffFile(path) as tf:
        page = tf.pages[0]
        num, den = page.tags['XResolution'].value
        unit = int(page.tags['ResolutionUnit'].value)
    return {2: 25400.0, 3: 10000.0}[unit] / (num / den)


def _load_hem(path):
    """
        Unclipped hematoxylin OD from the ROI's RGB, skimage rgb2hed directly.

        path (str): ROI TIFF path.

        Returns np.ndarray: float32 hematoxylin OD.
    """
    import tifffile
    from skimage.color import rgb2hed
    with tifffile.TiffFile(path) as tf:
        s = tf.series[0]
        lv = getattr(s, 'levels', None)
        arr = lv[0].asarray() if lv else s.asarray()
    rgb = np.ascontiguousarray(arr[:, :, :3])
    return rgb2hed(rgb.astype(np.float32) / 255.0)[:, :, 0].astype(np.float32)


def _window(img, cx, cy, W):
    """
        W x W window centred on round(cx), round(cy).

        img (np.ndarray): the channel.
        cx (float): x.
        cy (float): y.
        W (int): odd side.

        Returns np.ndarray: float32 copy.
    """
    h = W // 2
    ix, iy = int(round(cx)), int(round(cy))
    return np.array(img[iy - h:iy + h + 1, ix - h:ix + h + 1], dtype=np.float32, copy=True)


def _templates(hem, cx, cy, rect):
    """
        The seven arms' templates rebuilt from pixels, independently of the notebook's build_templates.

        hem (np.ndarray): hematoxylin OD.
        cx (float): click x.
        cy (float): click y.
        rect (tuple): (y0, y1, x0, x1) 51-local half-open.

        Returns dict: arm -> (template float32, fill value or nan, n_kept).
    """
    out = {}
    f51 = None
    for W in (51, 61, 71):
        win = _window(hem, cx, cy, W)
        off = (W - 51) // 2
        keep = np.zeros((W, W), bool)
        keep[rect[0] + off:rect[1] + off, rect[2] + off:rect[3] + off] = True
        out[f'default_{W}'] = (win.copy(), np.nan, W * W, keep)
        f = np.float32(win[~keep].astype(np.float64).mean())
        F = win.copy()
        F[~keep] = f
        out[f'fill_{W}'] = (F, f, int(keep.sum()), keep)
        if W == 51:
            f51 = f
        if W == 71:
            P = win.copy()
            P[~keep] = f51
            out['fill_71_pinned'] = (P, f51, int(keep.sum()), keep)
    return out


def _sha(T):
    """
        SHA-1 of a template's float32 bytes.

        T (np.ndarray): template.

        Returns str: hex digest.
    """
    return hashlib.sha1(np.ascontiguousarray(T, dtype=np.float32).tobytes()).hexdigest()


def _blank(hem, cx, cy, size, fill):
    """
        A copy of the channel with a size x size square centred on the rounded click set to fill.

        hem (np.ndarray): channel.
        cx (float): x.
        cy (float): y.
        size (int): odd side.
        fill (float): fill value.

        Returns np.ndarray: blanked copy.
    """
    h = size // 2
    ix, iy = int(round(cx)), int(round(cy))
    out = hem.copy()
    out[max(0, iy - h):iy + h + 1, max(0, ix - h):ix + h + 1] = fill
    return out


def _response(search, T):
    """
        Full-ROI TM_CCOEFF response at every ROI pixel, replicate-padded by (W - 1) // 2.

        search (np.ndarray): search channel.
        T (np.ndarray): W x W template.

        Returns np.ndarray: float32 response, same shape as search.
    """
    import cv2
    p = (T.shape[0] - 1) // 2
    sp = cv2.copyMakeBorder(search, p, p, p, p, borderType=cv2.BORDER_REPLICATE)
    return np.asarray(cv2.matchTemplate(np.ascontiguousarray(sp, dtype=np.float32), np.ascontiguousarray(T, dtype=np.float32), cv2.TM_CCOEFF), dtype=np.float32)


def _peaks(R, drop_near=None, max_peaks=100):
    """
        Local maxima (15 x 15 dilation), threshold median - 1.5 MAD on a stride-8 sample, best-first with an (x, y) tie-break, capped.

        R (np.ndarray): response.
        drop_near (tuple or None): (ix, iy): drop maxima within SELF_PX (Chebyshev) of it BEFORE the cap.
        max_peaks (int): the cap.

        Returns tuple[np.ndarray, np.ndarray]: (x, y) centres float64 and float32 scores.
    """
    import cv2
    s = R[::8, ::8].ravel()
    s = s[np.isfinite(s)]
    med = float(np.median(s))
    mad = float(1.4826 * np.median(np.abs(s - med)))
    thr = med + (-1.5) * mad
    dil = cv2.dilate(R, np.ones((15, 15), np.uint8))
    ys, xs = np.nonzero((R >= dil) & (R >= thr))
    sc = R[ys, xs]
    if drop_near is not None:
        keep = np.maximum(np.abs(xs - drop_near[0]), np.abs(ys - drop_near[1])) > SELF_PX
        ys, xs, sc = ys[keep], xs[keep], sc[keep]
    order = np.lexsort((ys, xs, -sc))[:max_peaks]
    return np.stack([xs[order], ys[order]], axis=1).astype(np.float64), sc[order]


def _nms(c, s, r):
    """
        Greedy distance NMS, descending score, stable on ties, suppression at distance <= r.

        c (np.ndarray): centres.
        s (np.ndarray): scores.
        r (float): radius.

        Returns np.ndarray: kept indices best-first.
    """
    order = np.argsort(-s.astype(np.float64), kind='stable')
    sup = np.zeros(len(c), bool)
    keep = []
    for i in order:
        if sup[i]:
            continue
        keep.append(i)
        sup |= np.hypot(c[:, 0] - c[i, 0], c[:, 1] - c[i, 1]) <= r
    return np.asarray(keep, int)


def _od(hem_pad25, c):
    """
        Mean of the darkest-10% (highest-OD) pixels of a 51 px window at each detection, on the unblanked channel padded by 25.

        hem_pad25 (np.ndarray): replicate-padded channel.
        c (np.ndarray): detection centres in unpadded coordinates.

        Returns np.ndarray: od per detection.
    """
    k = max(1, int(0.10 * 51 * 51))
    out = []
    for x, y in c.astype(int):
        flat = hem_pad25[y:y + 51, x:x + 51].ravel()
        out.append(float(np.partition(flat, -k)[-k:].mean()))
    return np.asarray(out)


def _search_to_lists(R, hem_pad25, gt_eval, r, drop_near=None):
    """
        Peaks -> NMS -> od -> ranking on both axes -> buckets.

        R (np.ndarray): response.
        hem_pad25 (np.ndarray): unblanked channel padded by 25.
        gt_eval (pd.DataFrame): GT without the seed annotation.
        r (float): NMS = match radius.
        drop_near (tuple or None): drop maxima near this pixel before the cap.

        Returns tuple[dict, int, int]: axis -> ranked frame (cx, cy, score, od, bucket, matched_ann_id), n_peaks, n_detections.
    """
    c, s = _peaks(R, drop_near)
    n_peaks = len(c)
    keep = _nms(c, s, r)
    c, s = c[keep], s[keep]
    det = pd.DataFrame({'cx': c[:, 0], 'cy': c[:, 1], 'score': s.astype(np.float32), 'od': _od(hem_pad25, c)})
    lists = {}
    for ax in AXES:
        ranked = det.sort_values(AXIS_KEY[ax], ascending=False, na_position='last', kind='mergesort').reset_index(drop=True)
        b, a = bucket(ranked[['cx', 'cy']].to_numpy(), gt_eval, r)
        lists[ax] = ranked.assign(rank=np.arange(len(ranked)), bucket=b, matched_ann_id=a)
    return lists, n_peaks, int(len(det))


def _local_delta(hem, cx, cy, T, blank_small, blank_big, fill):
    """
        Exact change in the TM_CCOEFF response when the blank grows from blank_small to blank_big: correlation of the
        zero-mean template with the difference image, which is nonzero only on the ring between the two squares.

        hem (np.ndarray): unblanked channel.
        cx (float): click x.
        cy (float): click y.
        T (np.ndarray): W x W template.
        blank_small (int): current blank side.
        blank_big (int): new blank side.
        fill (float): the blank's fill value.

        Returns tuple[np.ndarray, int, int]: delta response over the affected positions and its top-left ROI coordinate (x0, y0).
    """
    W = T.shape[0]
    ht, hb, hs = W // 2, blank_big // 2, blank_small // 2
    reach = hb + ht                                    # positions whose window touches the ring
    size = 2 * (reach + ht) + 1
    ix, iy = int(round(cx)), int(round(cy))
    D = np.zeros((size, size), np.float64)
    c0 = reach + ht
    big = hem[iy - hb:iy + hb + 1, ix - hb:ix + hb + 1].astype(np.float64)
    ring = fill - big
    ring[hb - hs:hb + hs + 1, hb - hs:hb + hs + 1] = 0.0
    D[c0 - hb:c0 + hb + 1, c0 - hb:c0 + hb + 1] = ring
    Tz = T.astype(np.float64) - T.astype(np.float64).mean()
    delta = spsig.correlate(D, Tz, mode='valid')     # (2 reach + 1)^2
    return delta, ix - reach, iy - reach


def tierb_roi(args):
    """
        Tier B for one ROI: rebuild all 7 templates per click, reproduce searches from pixels, and run the sensitivity variants.

        args (tuple): (file_name, subset, clicks records, gt frame, width, height, set of (seed) needing default reruns, set of seeds needing all 7 arms).

        Returns dict: templates rows, search summary rows, top-30 frames, timing.
    """
    import cv2
    cv2.setNumThreads(3)
    file_name, subset, clicks, gt, width, height, default_seeds, all7_seeds = args
    t0 = time.time()
    path = f'{IMG_DIRS[subset]}/{file_name}'
    hem = _load_hem(path)
    t_hem = time.time() - t0
    mpp = _roi_mpp(path)
    r = 7.5 / mpp
    hmin = float(hem.min())
    hp25 = cv2.copyMakeBorder(hem, 25, 25, 25, 25, borderType=cv2.BORDER_REPLICATE)
    tpl_rows, srch_rows, tops, identity_rows, selfmatch_rows = [], [], [], [], []
    t_mt = {51: [], 61: [], 71: []}
    for ck in clicks:
        cx, cy, s = float(ck['click_cx']), float(ck['click_cy']), int(ck['seed_index'])
        rect = (int(ck['rect_y0']), int(ck['rect_y1']), int(ck['rect_x0']), int(ck['rect_x1']))
        tpl = _templates(hem, cx, cy, rect)
        for arm, (T, f, nk, keep) in tpl.items():
            tpl_rows.append(dict(file_name=file_name, seed_index=s, condition=arm, sha=_sha(T), fill_value=float(f), n_kept=nk, n_filled=T.size - nk if arm.startswith('fill') else 0))
        gt_eval = gt[gt['ann_id'] != int(ck['seed_ann_id'])].reset_index(drop=True)
        search51 = _blank(hem, cx, cy, 51, hmin)
        arms_full = ['fill_61', 'fill_71']
        if s in default_seeds:
            arms_full += ['default_61', 'default_71']
        if s in all7_seeds:
            arms_full = list(ARMS)
        for arm in arms_full:
            T = tpl[arm][0]
            W = T.shape[0]
            t1 = time.time()
            R = _response(search51, T)
            t_mt[W].append(time.time() - t1)
            variants = [('repro_blank51', R, None, NEAR_PX)]
            if arm in ('fill_61', 'fill_71'):
                delta, x0, y0 = _local_delta(hem, cx, cy, T, 51, W, hmin)
                R2 = R.copy()
                ys0, xs0 = max(0, y0), max(0, x0)
                ys1, xs1 = min(height, y0 + delta.shape[0]), min(width, x0 + delta.shape[1])
                R2[ys0:ys1, xs0:xs1] = (R2[ys0:ys1, xs0:xs1].astype(np.float64) + delta[ys0 - y0:ys1 - y0, xs0 - x0:xs1 - x0]).astype(np.float32)
                variants.append((f'd11_blank{W}', R2, None, 2 * (W // 2)))   # a W window overlaps a W blank out to W//2 + W//2
                if W == 71 and s == int(clicks[0]['seed_index']):
                    # identity check: TM_CCOEFF of a fill template = sum_keep (t - f)(I - mu_W); evaluate at 5 peaks + the click
                    keepm = tpl[arm][3]
                    c_chk, _ = _peaks(R)
                    pts = [tuple(map(int, c_chk[j])) for j in range(5)] + [(int(round(cx)), int(round(cy)))]
                    hpad = cv2.copyMakeBorder(search51, 35, 35, 35, 35, borderType=cv2.BORDER_REPLICATE).astype(np.float64)
                    t_i = T.astype(np.float64)[keepm]
                    fW = float(tpl[arm][1])
                    m_k = t_i.mean()
                    for (px, py) in pts:
                        Iw = hpad[py:py + 71, px:px + 71]
                        muW = Iw.mean()
                        form1 = float(((t_i - fW) * (Iw[keepm] - muW)).sum())
                        form2 = float(((t_i - m_k) * Iw[keepm]).sum() + (m_k - fW) * keepm.sum() * (Iw[keepm].mean() - muW))
                        identity_rows.append(dict(file_name=file_name, seed_index=s, x=px, y=py, cv2_response=float(R[py, px]), formula_sum_keep=form1, formula_decomposed=form2))
            if arm in ('default_61', 'default_71'):
                variants.append(('selfmatch_removed', R, (int(round(cx)), int(round(cy))), NEAR_PX))
                if arm == 'default_71':
                    # closed form at the click: sum_ring T'(T - min), since the centre 51 is blank and sum T' = 0
                    Tz = T.astype(np.float64) - T.astype(np.float64).mean()
                    ringm = np.ones_like(T, bool)
                    ringm[10:61, 10:61] = False
                    closed = float((Tz[ringm] * (T.astype(np.float64)[ringm] - hmin)).sum())
                    selfmatch_rows.append(dict(file_name=file_name, seed_index=s, cv2_response_at_click=float(R[int(round(cy)), int(round(cx))]), closed_form_ring=closed))
            for vname, RR, drop, near_px in variants:
                lists, npk, ndet = _search_to_lists(RR, hp25, gt_eval, r, drop)
                for ax, L in lists.items():
                    hit = (L['bucket'].to_numpy() == TP)
                    cheb = np.maximum(np.abs(L['cx'].to_numpy()[:10] - cx), np.abs(L['cy'].to_numpy()[:10] - cy))
                    srch_rows.append(dict(file_name=file_name, seed_index=s, condition=arm, variant=vname, rank_key=ax, n_peaks=npk, n_detections=ndet, tp_at_10=int(hit[:10].sum()), tp_at_20=int(hit[:20].sum()), tp_at_30=int(hit[:30].sum()), near_px=near_px, n_near_seed_top10_cheb=int((cheb <= NEAR_PX).sum()), n_near_seed_fp_top10_cheb=int(((cheb <= NEAR_PX) & ~hit[:10]).sum()), n_near_seed_top10_cheb_reach=int((cheb <= near_px).sum()), n_near_seed_fp_top10_cheb_reach=int(((cheb <= near_px) & ~hit[:10]).sum()), n_selfmatch_top30=int((np.maximum(np.abs(L['cx'].to_numpy()[:30] - cx), np.abs(L['cy'].to_numpy()[:30] - cy)) <= SELF_PX).sum())))
                    tops.append(L.head(30).assign(file_name=file_name, seed_index=s, condition=arm, variant=vname, rank_key=ax))
            del R
    return dict(file_name=file_name, tpl=tpl_rows, srch=srch_rows, top=pd.concat(tops, ignore_index=True) if tops else pd.DataFrame(), identity=identity_rows, selfmatch=selfmatch_rows, t_hem=t_hem, t_mt=t_mt, t_total=time.time() - t0, mpp=mpp)


# ---------------------------------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------------------------------

def main():
    T_START = time.time()
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 60)
    pd.set_option('display.max_rows', 500)

    # ----------------------------------------------------------------------------------------------- A0
    section('A0  execution-coherence gate')
    nb = json.load(open(NB_PATH))
    code = [(i, c) for i, c in enumerate(nb['cells']) if c['cell_type'] == 'code']
    ecs = [c.get('execution_count') for _, c in code]
    errs = [i for i, c in code if any(o.get('output_type') == 'error' for o in c.get('outputs', []))]
    unrun = [i for i, c in code if c.get('execution_count') is None]
    starts = {i: c['metadata'].get('execution', {}).get('iopub.execute_input') for i, c in code}
    ends = {i: c['metadata'].get('execution', {}).get('shell.execute_reply') for i, c in code}
    print(f'code cells {len(code)}; execution counts {ecs}; contiguous from 1: {ecs == list(range(1, len(code) + 1))}; error cells {errs}; unrun {unrun}')
    gaps = []
    prev = None
    for i, _ in code:
        if prev is not None and starts[i] and ends[prev]:
            g = (pd.Timestamp(starts[i]) - pd.Timestamp(ends[prev])).total_seconds()
            gaps.append(dict(cell=i, gap_s_after_previous_code_cell=round(g, 1), start=starts[i]))
        prev = i
    gaps = pd.DataFrame(gaps)
    print(gaps.to_string(index=False))
    exec_gate = pd.DataFrame([dict(item='n_code_cells', value=len(code)), dict(item='execution_counts', value=str(ecs)), dict(item='contiguous_from_1', value=ecs == list(range(1, len(code) + 1))), dict(item='error_cells', value=str(errs)), dict(item='unrun_cells', value=str(unrun)), dict(item='last_cell_type', value=nb['cells'][-1]['cell_type']), dict(item='largest_gap_s', value=float(gaps['gap_s_after_previous_code_cell'].max())), dict(item='largest_gap_before_cell', value=int(gaps.loc[gaps['gap_s_after_previous_code_cell'].idxmax(), 'cell']))])
    exec_gate.to_csv(f'{OUT}_execution_gate.csv', index=False)

    # ----------------------------------------------------------------------------------------------- A1
    section('A1  composition gate')
    run = pd.read_csv(f'{WS}_per_run.csv', float_precision='round_trip')
    top = pd.read_csv(f'{WS}_top30.csv', float_precision='round_trip')
    per_roi = pd.read_csv(f'{WS}_per_roi.csv', float_precision='round_trip')
    dpr = pd.read_csv(f'{WS}_delta_per_roi.csv', float_precision='round_trip')
    dst = pd.read_csv(f'{WS}_delta_stats.csv', float_precision='round_trip')
    clicks = pd.read_csv(f'{FG}_clicks.csv', float_precision='round_trip')
    ref_run = pd.read_csv(f'{FG}_per_run.csv', float_precision='round_trip')
    ref_top = pd.read_csv(f'{FG}_top30.csv', float_precision='round_trip')
    ref_stats = pd.read_csv(f'{FG}_delta_stats.csv', float_precision='round_trip')
    ims, anns = load_gt()
    roi_order = list(pd.unique(ref_run['file_name']))
    ce = clicks.join(ims.set_index('file_name')[['width', 'height']], on='file_name')
    ce_min = np.minimum.reduce([ce['click_cx'], ce['click_cy'], ce['width'] - 1 - ce['click_cx'], ce['height'] - 1 - ce['click_cy']])
    print(f'nearest click to an ROI edge: {ce_min.min():.0f} px (>= 36 keeps the 71 px blank ring off the edge row)')
    assert ce_min.min() >= 36
    comp = []

    def comp_row(item, expected, observed):
        comp.append(dict(item=item, expected=str(expected), observed=str(observed), ok=str(expected) == str(observed)))

    comp_row('per_run rows = 49 x 6 x 7 x 2', 4116, len(run))
    comp_row('per_run duplicate (file, seed, condition, rank_key)', 0, int(run.duplicated(['file_name', 'seed_index', 'condition', 'rank_key']).sum()))
    comp_row('per_run rows per (condition, rank_key)', '{294}', set(run.groupby(['condition', 'rank_key']).size()))
    comp_row('per_run ROIs', 49, run['file_name'].nunique())
    comp_row('per_run seeds per ROI', '{6}', set(run.groupby('file_name')['seed_index'].nunique()))
    comp_row('per_run arms', sorted(ARMS), sorted(run['condition'].unique()))
    comp_row('top30 rows = 4116 x 30', 123480, len(top))
    g = top.groupby(['file_name', 'seed_index', 'condition', 'rank_key'])['rank']
    comp_row('top30 lists with ranks exactly 0..29', 4116, int(g.apply(lambda v: sorted(v) == list(range(30))).sum()))
    comp_row('per_roi rows = 49 x 7 x 2', 686, len(per_roi))
    comp_row('delta_per_roi rows = 46 x 49', 2254, len(dpr))
    comp_row('delta_stats rows', 46, len(dst))
    comp_row('delta_stats tested rows', 8, int((dst['role'] == 'tested').sum()))
    meta = run.drop_duplicates(['file_name', 'seed_index'])[['file_name', 'seed_index', 'click_cx', 'click_cy', 'seed_ann_id', 'mpp', 'domain', 'subset']].merge(clicks[['file_name', 'seed_index', 'click_cx', 'click_cy', 'seed_ann_id', 'mpp', 'domain', 'subset']], on=['file_name', 'seed_index'], suffixes=('', '_c'), how='outer', indicator=True)
    comp_row('per_run clicks found in clicks CSV (outer join both)', 294, int((meta['_merge'] == 'both').sum()))
    n_meta_diff = int(sum((meta[c] != meta[c + '_c']).sum() for c in ('click_cx', 'click_cy', 'seed_ann_id', 'mpp', 'domain', 'subset')))
    comp_row('per_run click fields != clicks CSV (6 fields x 294)', 0, n_meta_diff)
    dom = run.drop_duplicates('file_name').merge(ims[['file_name', 'tumor_type']], on='file_name', how='left')
    comp_row('per_run domain == MIDOG++.json tumor_type (49 ROIs)', 49, int((dom['domain'] == dom['tumor_type']).sum()))
    comp_row('domains in sample', 7, dom['domain'].nunique())
    comp_row('ROIs per domain in sample', '{7}', set(dom.groupby('domain').size()))
    ann_ims = ims[ims['file_name'].isin(anns['file_name'])]
    comp_row('tumour types in MIDOG++.json (annotated images)', 7, ann_ims['tumor_type'].nunique())
    comp_row('ROI order: parent per_run first appearance == this per_run first appearance', True, roi_order == list(pd.unique(run['file_name'])))
    comp = pd.DataFrame(comp)
    print(comp.to_string(index=False))
    print('\nMIDOG++.json annotated ROIs per tumour type (whole database):')
    print(ann_ims.groupby('tumor_type').size().to_string())
    comp.to_csv(f'{OUT}_composition.csv', index=False)
    ledger('composition gate items', len(comp), int((~comp['ok']).sum()))

    # ----------------------------------------------------------------------------------------------- A2
    section('A2  every top-30 list re-matched independently against MIDOG++.json; TP@10/20/30 from the buckets')
    gt_by_file = {f: g.reset_index(drop=True) for f, g in anns.groupby('file_name')}
    cmeta = clicks.set_index(['file_name', 'seed_index'])
    top = top.sort_values(['file_name', 'seed_index', 'condition', 'rank_key', 'rank']).reset_index(drop=True)
    my_b, my_a = np.empty(len(top), dtype=object), np.empty(len(top), dtype=np.int64)
    for (f, s, cond, ax), gi in top.groupby(['file_name', 'seed_index', 'condition', 'rank_key'], sort=False).groups.items():
        m = cmeta.loc[(f, s)]
        gt_eval = gt_by_file[f][gt_by_file[f]['ann_id'] != int(m['seed_ann_id'])]
        b, a = bucket(top.loc[gi, ['cx', 'cy']].to_numpy(float), gt_eval, 7.5 / float(m['mpp']))
        my_b[gi], my_a[gi] = b, a
    top['my_bucket'], top['my_ann'] = my_b, my_a
    ledger('top30 bucket re-matched (123,480 detections)', len(top), int((top['bucket'] != top['my_bucket']).sum()))
    ledger('top30 matched_ann_id re-matched', len(top), int((top['matched_ann_id'] != top['my_ann']).sum()))
    top['hit'] = top['my_bucket'] == TP
    tpk = top.groupby(['file_name', 'seed_index', 'condition', 'rank_key']).apply(lambda g: pd.Series({f'tp_at_{k}': int(g['hit'].to_numpy()[:k].sum()) for k in (10, 20, 30)})).reset_index()
    chk = run.merge(tpk, on=['file_name', 'seed_index', 'condition', 'rank_key'], suffixes=('', '_mine'))
    ledger('per_run tp_at_10/20/30 from re-matched top30', 3 * len(chk), int(sum((chk[f'tp_at_{k}'] != chk[f'tp_at_{k}_mine']).sum() for k in (10, 20, 30))))
    # W = 51 parity from the artifacts: this run's default_51 / fill_51 top30 vs the parent's default_51 / hem_rect_fill
    par = []
    for arm, cond in (('default_51', 'default_51'), ('fill_51', 'hem_rect_fill')):
        a = top[top['condition'] == arm].sort_values(['file_name', 'seed_index', 'rank_key', 'rank']).reset_index(drop=True)
        b = ref_top[ref_top['condition'] == cond].sort_values(['file_name', 'seed_index', 'rank_key', 'rank']).reset_index(drop=True)
        assert len(a) == len(b) == 17640
        for c in ('cx', 'cy', 'bucket', 'matched_ann_id'):
            par.append((f'{arm} {c}', len(a), int((a[c].to_numpy() != b[c].to_numpy()).sum())))
        par.append((f'{arm} score (float32)', len(a), int((a['score'].to_numpy().astype(np.float32) != b['score'].to_numpy().astype(np.float32)).sum())))
        par.append((f'{arm} od (exact)', len(a), int((a['od'].to_numpy() != b['od'].to_numpy()).sum())))
        ra = run[run['condition'] == arm].sort_values(['file_name', 'seed_index', 'rank_key']).reset_index(drop=True)
        rb = ref_run[ref_run['condition'] == cond].sort_values(['file_name', 'seed_index', 'rank_key']).reset_index(drop=True)
        for c in ('tp_at_10', 'tp_at_20', 'tp_at_30', 'n_detections', 'n_peaks', 'tpl_sha1'):
            par.append((f'{arm} per_run {c}', len(ra), int((ra[c].to_numpy() != rb[c].to_numpy()).sum())))
    for name, n, dv in par:
        ledger(f'W=51 parity from artifacts: {name}', n, dv)

    # ----------------------------------------------------------------------------------------------- A3
    section('A3  per-ROI sums and pooled precision')
    metrics = list(METRIC_K)
    roi = run.groupby(['file_name', 'condition', 'rank_key'])[metrics].sum()
    pr = per_roi.set_index(['file_name', 'condition', 'rank_key'])
    ledger('per_roi tp sums (686 x 4)', 686 * 4, int((roi.loc[pr.index, metrics].to_numpy() != pr[metrics].to_numpy()).sum()))
    ledger('per_roi precision columns (686 x 4)', 686 * 4, int(sum((~np.isclose(pr[m.replace('tp_at', 'precision_at')], pr[m] / (6 * k), rtol=1e-12, atol=0)).sum() for m, k in METRIC_K.items())))
    pooled = run.groupby(['rank_key', 'condition'])[metrics].sum()
    for m, k in METRIC_K.items():
        pooled[m] = pooled[m] / (6 * 49 * k)
    pooled.columns = [m.replace('tp_at', 'precision_at') for m in metrics]
    pooled = pooled.reset_index()
    printed_pooled = {('chromatin_od', 'default_51'): (0.5667, 0.5844, 0.5029, 0.4351), ('chromatin_od', 'fill_51'): (0.5680, 0.5861, 0.5073, 0.4437), ('chromatin_od', 'default_61'): (0.5585, 0.5779, 0.4915, 0.4230), ('chromatin_od', 'fill_61'): (0.5622, 0.5827, 0.5066, 0.4426), ('chromatin_od', 'default_71'): (0.5493, 0.5622, 0.4718, 0.4052), ('chromatin_od', 'fill_71'): (0.5565, 0.5769, 0.5051, 0.4422), ('chromatin_od', 'fill_71_pinned'): (0.5534, 0.5748, 0.5029, 0.4412),
                      ('tm_score', 'default_51'): (0.5085, 0.5065, 0.4361, 0.3906), ('tm_score', 'fill_51'): (0.5265, 0.5276, 0.4595, 0.4096), ('tm_score', 'default_61'): (0.4891, 0.4908, 0.4206, 0.3704), ('tm_score', 'fill_61'): (0.5425, 0.5439, 0.4724, 0.4188), ('tm_score', 'default_71'): (0.4456, 0.4456, 0.3854, 0.3422), ('tm_score', 'fill_71'): (0.5466, 0.5500, 0.4724, 0.4161), ('tm_score', 'fill_71_pinned'): (0.5439, 0.5480, 0.4730, 0.4161)}
    n_bad = 0
    for (ax, arm), vals in printed_pooled.items():
        row = pooled[(pooled['rank_key'] == ax) & (pooled['condition'] == arm)].iloc[0]
        for col, v in zip(('precision_at_10', 'precision_at_10_excl_border', 'precision_at_20', 'precision_at_30'), vals):
            n_bad += int(round(row[col], 4) != v)
    ledger('printed pooled precision table (cell 12), 4 d.p.', 56, n_bad)
    pooled.to_csv(f'{OUT}_pooled_precision.csv', index=False)
    print(pooled.to_string(index=False, float_format=lambda v: f'{v:.4f}'))

    # ----------------------------------------------------------------------------------------------- A4
    section('A4  contrasts: per-ROI deltas and every delta_stats column, recomputed')
    roi10 = {m: run.pivot_table(index=['rank_key', 'file_name'], columns='condition', values=m, aggfunc='sum') for m in metrics}
    seed_tab = {m: run.pivot_table(index=['rank_key', 'file_name', 'seed_index'], columns='condition', values=m, aggfunc='sum') for m in metrics}
    idx49 = boot_idx(49, 0)
    spec = [(ax, 'tp_at_10', c) for ax in AXES for c in CONTRASTS] + [(ax, m, c) for ax in AXES for m in metrics if m != 'tp_at_10' for c in TESTED]
    mine_rows, n_dpr_bad, n_dpr = [], 0, 0
    dpr_ix = dpr.set_index(['rank_key', 'metric', 'contrast', 'file_name'])
    for ax, m, cname in spec:
        coefs = CONTRASTS[cname]
        tab = roi10[m].loc[ax].reindex(roi_order)
        d = sum(c * tab[a] for a, c in coefs.items()).to_numpy().astype(np.int64)
        st_tab = seed_tab[m].loc[ax]
        for f, dv in zip(roi_order, d):
            rec = dpr_ix.loc[(ax, m, cname, f)]
            n_dpr += 7
            n_dpr_bad += int(rec['delta_tp_sum'] != dv)
            for s in SEEDS:
                n_dpr_bad += int(rec[f'delta_tp_s{s}'] != sum(c * st_tab.loc[(f, s), a] for a, c in coefs.items()))
        st = contrast_stats(d, METRIC_K[m], idx49)
        tested = m == 'tp_at_10' and cname in TESTED
        mine_rows.append(dict(rank_key=ax, metric=m, contrast=cname, role='tested' if tested else 'descriptive', holm_family=f'{ax} | {TESTED[cname]}' if tested else '', **st))
    ledger('delta_per_roi delta_tp_sum and 6 per-seed deltas (2254 x 7)', n_dpr, n_dpr_bad)
    mine = pd.DataFrame(mine_rows)
    mine['holm_p'] = np.nan
    mine['family_size'] = np.nan
    for fam, g in mine[mine['role'] == 'tested'].groupby('holm_family'):
        mine.loc[g.index, 'holm_p'] = holm(g['exact_p'].to_numpy())
        mine.loc[g.index, 'family_size'] = len(g)
    mine['holm_adjusted_floor'] = mine['min_attainable_p'] * mine['family_size']
    rec = dst.set_index(['rank_key', 'metric', 'contrast'])
    mi = mine.set_index(['rank_key', 'metric', 'contrast']).loc[rec.index]
    n_cmp, n_bad = 0, 0
    exact_cols = ('wins', 'losses', 'ties', 'n_nonzero', 'n_same_direction', 'ci_excludes_0', 't_ci_excludes_0', 'majority_of_all_rois', 'd5_bar_met', 'role')
    float_cols = ('pooled_delta_precision', 'exact_p', 'min_attainable_p', 'ci_low', 'ci_high', 'boot_share_ge_0', 'boot_share_le_0', 't_ci_low', 't_ci_high', 'holm_p', 'holm_adjusted_floor')
    for c in exact_cols:
        n_cmp += len(rec)
        n_bad += int((mi[c].astype(str).to_numpy() != rec[c].astype(str).to_numpy()).sum())
    for c in float_cols:
        a, b = mi[c].to_numpy(float), rec[c].to_numpy(float)
        both_nan = np.isnan(a) & np.isnan(b)
        n_cmp += len(rec)
        n_bad += int((~(both_nan | np.isclose(a, b, rtol=1e-9, atol=1e-15))).sum())
    ledger('delta_stats: 21 columns x 46 rows', n_cmp, n_bad, 'floats rtol 1e-9; bootstrap drawn from the documented numpy stream')
    tr = mine[mine['role'] == 'tested'].copy()
    tr['holm_all8'] = holm(tr['exact_p'].to_numpy())
    tr['bonferroni_all8'] = np.minimum(1.0, 8 * tr['exact_p'])
    print('multiplicity sensitivity -- Holm across all 8 tested rows instead of 4 families of 2:')
    print(tr[['rank_key', 'contrast', 'exact_p', 'holm_p', 'holm_all8', 'bonferroni_all8']].to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    tr[['rank_key', 'contrast', 'exact_p', 'holm_p', 'holm_all8', 'bonferroni_all8']].to_csv(f'{OUT}_multiplicity.csv', index=False)
    mine.to_csv(f'{OUT}_delta_stats_recomputed.csv', index=False)
    print(mine[mine['role'] == 'tested'][['rank_key', 'contrast', 'pooled_delta_precision', 'ci_low', 'ci_high', 't_ci_low', 't_ci_high', 'wins', 'losses', 'ties', 'n_nonzero', 'exact_p', 'holm_p', 'd5_bar_met']].to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    # re-seeded bootstrap stability of the 5 tested CIs that exclude 0
    stab = []
    for r_ in mine[(mine['role'] == 'tested') & mine['ci_excludes_0']].itertuples():
        tab = roi10['tp_at_10'].loc[r_.rank_key].reindex(roi_order)
        d = sum(c * tab[a] for a, c in CONTRASTS[r_.contrast].items()).to_numpy().astype(np.int64)
        bounds = []
        for sd in range(20):
            sums = d[boot_idx(49, 1 + sd)].sum(1) / 2940
            bounds.append(np.percentile(sums, [2.5, 97.5]))
        closest = min(abs(lo) if r_.pooled_delta_precision > 0 else abs(hi) for lo, hi in bounds)
        stab.append(dict(rank_key=r_.rank_key, contrast=r_.contrast, excludes_0_all_20=all(lo > 0 or hi < 0 for lo, hi in bounds), closest_bound_pts=round(100 * closest, 3)))
    stab = pd.DataFrame(stab)
    printed_stab = {('chromatin_od', 'fill_71 - fill_51'): 0.272, ('tm_score', 'fill_61 - fill_51'): 0.782, ('tm_score', 'fill_71 - fill_51'): 0.714, ('tm_score', '(fill_61 - fill_51) - (default_61 - default_51)'): 2.279, ('tm_score', '(fill_71 - fill_51) - (default_71 - default_51)'): 6.190}
    ledger('re-seeded bootstrap stability table (cell 14), 5 rows x 2', 10, int(sum((round(r_.closest_bound_pts, 3) != printed_stab[(r_.rank_key, r_.contrast)]) + (not r_.excludes_0_all_20) for r_ in stab.itertuples())))
    print(stab.to_string(index=False))

    # parent validation (cell 4): the parent's 10 delta_stats rows from the parent's per_run
    ref_roi10 = ref_run.pivot_table(index=['rank_key', 'file_name'], columns='condition', values='tp_at_10', aggfunc='sum')
    pv = []
    for r_ in ref_stats.itertuples():
        a, b = r_.pair.split(' - ')
        tab = ref_roi10.loc[r_.rank_key].reindex(roi_order)
        st = contrast_stats((tab[a] - tab[b]).to_numpy().astype(np.int64), int(r_.K), idx49)
        pv.append(dict(rank_key=r_.rank_key, pair=r_.pair, holm_family=r_.holm_family, **st))
    pv = pd.DataFrame(pv)
    pv['holm_p'] = np.nan
    for fam, g in pv.groupby('holm_family'):
        pv.loc[g.index, 'holm_p'] = holm(g['exact_p'].to_numpy())
    nb_ = 0
    for c in ('pooled_delta_precision', 'exact_p', 'ci_low', 'ci_high', 't_ci_low', 't_ci_high', 'holm_p'):
        nb_ += int((~np.isclose(pv[c].to_numpy(float), ref_stats[c].to_numpy(float), rtol=1e-9, atol=1e-15)).sum())
    for c in ('wins', 'losses', 'ties', 'n_nonzero', 'd5_bar_met'):
        nb_ += int((pv[c].astype(str).to_numpy() != ref_stats[c].astype(str).to_numpy()).sum())
    ledger("parent fill-geometry delta_stats (cell 4's validation), 10 rows x 12 cols", 120, nb_)

    # ----------------------------------------------------------------------------------------------- A5
    section('A5  near-seed and self-match diagnostics, from top30 + clicks')
    t = top.merge(clicks[['file_name', 'seed_index', 'click_cx', 'click_cy']], on=['file_name', 'seed_index'], validate='many_to_one')
    t['dx'], t['dy'] = (t['cx'] - t['click_cx']).abs(), (t['cy'] - t['click_cy']).abs()
    t['cheb'], t['euc'] = np.maximum(t['dx'], t['dy']), np.hypot(t['dx'], t['dy'])
    t10 = t[t['rank'] < 10]
    ns = t10.assign(cheb_all=t10['cheb'] <= NEAR_PX, cheb_fp=(t10['cheb'] <= NEAR_PX) & ~t10['hit'], euc_all=t10['euc'] <= NEAR_PX, euc_fp=(t10['euc'] <= NEAR_PX) & ~t10['hit']).groupby(['file_name', 'seed_index', 'condition', 'rank_key'])[['cheb_all', 'cheb_fp', 'euc_all', 'euc_fp']].sum().reset_index()
    chk = run.merge(ns, on=['file_name', 'seed_index', 'condition', 'rank_key'])
    nbad = sum(int((chk[a] != chk[b]).sum()) for a, b in (('n_near_seed_top10_cheb', 'cheb_all'), ('n_near_seed_fp_top10_cheb', 'cheb_fp'), ('n_near_seed_top10_euclid', 'euc_all'), ('n_near_seed_fp_top10_euclid', 'euc_fp')))
    ledger('per_run near-seed columns from top30 (4116 x 4)', 4 * len(chk), nbad)
    near_tab = ns.groupby(['rank_key', 'condition'])[['cheb_all', 'cheb_fp', 'euc_all', 'euc_fp']].sum().reset_index()
    printed_near = {('chromatin_od', 'default_51'): (4, 0, 2, 0), ('tm_score', 'default_51'): (3, 0, 2, 0), ('chromatin_od', 'fill_51'): (4, 0, 2, 0), ('tm_score', 'fill_51'): (3, 0, 2, 0), ('chromatin_od', 'default_61'): (4, 0, 2, 0), ('tm_score', 'default_61'): (5, 3, 4, 3), ('chromatin_od', 'fill_61'): (4, 0, 2, 0), ('tm_score', 'fill_61'): (3, 0, 2, 0), ('chromatin_od', 'default_71'): (8, 4, 6, 4), ('tm_score', 'default_71'): (19, 17, 18, 17), ('chromatin_od', 'fill_71'): (4, 0, 2, 0), ('tm_score', 'fill_71'): (3, 0, 2, 0), ('chromatin_od', 'fill_71_pinned'): (4, 0, 2, 0), ('tm_score', 'fill_71_pinned'): (3, 0, 2, 0)}
    nt = near_tab.set_index(['rank_key', 'condition'])[['cheb_all', 'cheb_fp', 'euc_all', 'euc_fp']]
    ledger('printed near-seed table (cell 10), 14 x 4', 56, int(sum(sum(int(x != y) for x, y in zip(nt.loc[key].to_numpy(), v)) for key, v in printed_near.items())))
    print(near_tab.to_string(index=False))
    # the non-TP near-seed detections: how far from the click, which bucket
    nfp = t10[(t10['cheb'] <= NEAR_PX) & ~t10['hit']]
    print('\nnear-seed non-TP detections: Chebyshev distance range and buckets per arm/axis')
    print(nfp.groupby(['condition', 'rank_key']).agg(n=('rank', 'size'), n_clicks=('seed_index', 'size'), cheb_min=('cheb', 'min'), cheb_max=('cheb', 'max'), buckets=('my_bucket', lambda b: ';'.join(sorted(set(b))))).to_string())
    sm = t[t['cheb'] <= SELF_PX]
    sm_tab = sm.groupby(['condition', 'rank_key']).agg(n_lists=('seed_index', 'size'), in_top10=('rank', lambda r: int((r < 10).sum())), at_rank0=('rank', lambda r: int((r == 0).sum())), rank_median=('rank', 'median'), buckets=('my_bucket', lambda b: ';'.join(sorted(set(b))))).reset_index()
    print('\nself-match (<= 2 px) detections in the top 30:')
    print(sm_tab.to_string(index=False))
    ledger('post-run cell 21 self-match table (default_61 tm 3/3; default_71 od 13/4, tm 24/17)', 6, int(sum([
        tuple(sm_tab[(sm_tab['condition'] == 'default_61') & (sm_tab['rank_key'] == 'tm_score')][['n_lists', 'in_top10']].iloc[0]) != (3, 3),
        tuple(sm_tab[(sm_tab['condition'] == 'default_71') & (sm_tab['rank_key'] == 'chromatin_od')][['n_lists', 'in_top10']].iloc[0]) != (13, 4),
        tuple(sm_tab[(sm_tab['condition'] == 'default_71') & (sm_tab['rank_key'] == 'tm_score')][['n_lists', 'in_top10']].iloc[0]) != (24, 17)])) * 2)
    sm_tab.to_csv(f'{OUT}_selfmatch_top30.csv', index=False)
    near_tab.to_csv(f'{OUT}_near_seed.csv', index=False)
    # identity of near-seed top-10 detections: fill_61 / fill_71 / pinned vs fill_51, per click and axis
    ident = []
    ns_sets = {k: set(zip(g['cx'], g['cy'])) for k, g in t10[t10['cheb'] <= NEAR_PX].groupby(['rank_key', 'condition', 'file_name', 'seed_index'])}
    for ax in AXES:
        for arm in ('fill_61', 'fill_71', 'fill_71_pinned', 'default_61', 'default_71'):
            same, n_nonempty = 0, 0
            for f in roi_order:
                for s in SEEDS:
                    a_ = ns_sets.get((ax, arm, f, s), set())
                    b_ = ns_sets.get((ax, 'fill_51' if arm.startswith('fill') else 'default_51', f, s), set())
                    same += int(a_ == b_)
                    n_nonempty += int(bool(a_ or b_))
                    if a_ != b_ and arm.startswith('fill'):
                        dif = t10[(t10['rank_key'] == ax) & (t10['file_name'] == f) & (t10['seed_index'] == s) & (t10['condition'].isin([arm, 'fill_51'])) & (t10['cheb'] <= NEAR_PX)]
                        print(f'  near-seed set differs: {ax} {arm} vs fill_51 at {f} s{s}:')
                        print(dif[['condition', 'rank', 'cx', 'cy', 'cheb', 'my_bucket', 'my_ann']].to_string(index=False))
            ident.append(dict(rank_key=ax, arm=arm, reference='fill_51' if arm.startswith('fill') else 'default_51', clicks_identical_near_seed_set=same, clicks=294, clicks_with_any_near_seed=n_nonempty))
    ident = pd.DataFrame(ident)
    print('\nnear-seed (Chebyshev <= 60) top-10 detection SETS, by exact (cx, cy), against the W = 51 arm of the same template kind:')
    print(ident.to_string(index=False))
    ident.to_csv(f'{OUT}_near_seed_identity.csv', index=False)

    # ----------------------------------------------------------------------------------------------- A6
    section('A6  self-match bound: drop <= 2 px detections from each top-30 list, re-take the top 10, recompute the contrasts')
    assert (sm['my_bucket'] == 'non_human_findings').all(), 'a self-match detection claimed a GT; dropping it would change later matches'
    t_ns = t[t['cheb'] > SELF_PX].copy()
    t_ns['r2'] = t_ns.groupby(['file_name', 'seed_index', 'condition', 'rank_key']).cumcount()
    corr = t_ns[t_ns['r2'] < 10].groupby(['file_name', 'seed_index', 'condition', 'rank_key'])['hit'].sum().rename('tp10_selfmatch_dropped').reset_index()
    runc = run.merge(corr, on=['file_name', 'seed_index', 'condition', 'rank_key'])
    runc['gain'] = runc['tp10_selfmatch_dropped'] - runc['tp_at_10']
    gg = runc.groupby(['condition', 'rank_key'])['gain'].agg(['sum', 'max'])
    print(gg[gg['sum'] != 0].to_string())
    # extra headroom from NMS un-suppression: a mitotic GT within 2 r of the click (the only TPs a peak at the click could have suppressed)
    rows2 = []
    for r_ in clicks.itertuples():
        g_ = gt_by_file[r_.file_name]
        g_ = g_[(g_['category_id'] == MITOTIC) & (g_['ann_id'] != r_.seed_ann_id)]
        dd = np.hypot(g_['cx'] - r_.click_cx, g_['cy'] - r_.click_cy)
        rr = 7.5 / r_.mpp
        rows2.append(dict(file_name=r_.file_name, seed_index=r_.seed_index, n_mitotic_within_r=int((dd <= rr).sum()), n_mitotic_within_2r=int((dd <= 2 * rr).sum()), nearest_mitotic_px=float(dd.min()) if len(dd) else np.nan))
    near_gt = pd.DataFrame(rows2)
    print(f"\nclicks with another mitotic GT within r: {int((near_gt['n_mitotic_within_r'] > 0).sum())}; within 2 r: {int((near_gt['n_mitotic_within_2r'] > 0).sum())} (of 294)")
    near_gt.to_csv(f'{OUT}_mitotic_near_click.csv', index=False)
    sb_rows = []
    for ax in AXES:
        piv = runc[runc['rank_key'] == ax].pivot_table(index='file_name', columns='condition', values='tp10_selfmatch_dropped', aggfunc='sum').reindex(roi_order)
        piv0 = roi10['tp_at_10'].loc[ax].reindex(roi_order)
        # worst case for the DiD: every top-10 self-match slot on a default arm was a TP displaced, plus one TP for every click with a mitotic GT within 2 r
        sm_top10 = sm[(sm['rank'] < 10) & (sm['rank_key'] == ax)].groupby(['condition', 'file_name']).size()
        ng = near_gt.groupby('file_name')['n_mitotic_within_2r'].apply(lambda v: int((v > 0).sum())).reindex(roi_order)
        worst = piv0.copy()
        for arm in ('default_61', 'default_71'):
            add = sm_top10.loc[arm].reindex(roi_order).fillna(0).astype(int) if arm in sm_top10.index.get_level_values(0) else pd.Series(0, index=roi_order)
            worst[arm] = piv0[arm] + add + ng
        for cname in ('default_61 - default_51', 'default_71 - default_51', '(fill_61 - fill_51) - (default_61 - default_51)', '(fill_71 - fill_51) - (default_71 - default_51)'):
            for label, tab in (('as recorded', piv0), ('self-match dropped from top30', piv), ('worst-case bound', worst)):
                d = sum(c * tab[a] for a, c in CONTRASTS[cname].items()).to_numpy().astype(np.int64)
                st = contrast_stats(d, 10, idx49)
                sb_rows.append(dict(rank_key=ax, contrast=cname, version=label, **{k: st[k] for k in ('pooled_delta_precision', 'ci_low', 'ci_high', 'wins', 'losses', 'ties', 'exact_p', 'd5_bar_met')}))
    sb = pd.DataFrame(sb_rows)
    print(sb.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    sb.to_csv(f'{OUT}_selfmatch_bound.csv', index=False)

    # ----------------------------------------------------------------------------------------------- A7
    section('A7  ROI-edge band: band share and tp_at_10_excl_border recomputed from top30 with fresh matching')
    t = t.join(ims.set_index('file_name')[['width', 'height']], on='file_name')
    t['edge'] = np.minimum.reduce([t['cx'], t['cy'], t['width'] - 1 - t['cx'], t['height'] - 1 - t['cy']])
    band10 = t[t['rank'] < 10].groupby(['file_name', 'seed_index', 'condition', 'rank_key'])['edge'].apply(lambda e: int((e < EDGE_PX).sum())).rename('band10').reset_index()
    xb = []
    short = 0
    for (f, s, cond, ax), g in t.groupby(['file_name', 'seed_index', 'condition', 'rank_key'], sort=False):
        kept = g[g['edge'] >= EDGE_PX]
        if len(kept) < 10:
            short += 1
        m = cmeta.loc[(f, s)]
        gt_eval = gt_by_file[f][gt_by_file[f]['ann_id'] != int(m['seed_ann_id'])]
        b, _ = bucket(kept[['cx', 'cy']].to_numpy(float)[:10], gt_eval, 7.5 / float(m['mpp']))
        xb.append(dict(file_name=f, seed_index=s, condition=cond, rank_key=ax, tpxb_mine=int((b == TP).sum())))
    xb = pd.DataFrame(xb).merge(band10, on=['file_name', 'seed_index', 'condition', 'rank_key'])
    chk = run.merge(xb, on=['file_name', 'seed_index', 'condition', 'rank_key'])
    ledger('per_run tp_at_10_excl_border re-matched from top30', len(chk), int((chk['tp_at_10_excl_border'] != chk['tpxb_mine']).sum()), f'{short} lists had <10 non-band rows in the top 30')
    ledger('per_run n_border_band_top10', len(chk), int((chk['n_border_band_top10'] != chk['band10']).sum()))
    edge_tab = chk.groupby(['rank_key', 'condition']).agg(band_slots=('band10', 'sum'), n=('band10', 'size')).reset_index()
    edge_tab['band_share_pct'] = 100 * edge_tab['band_slots'] / (10 * edge_tab['n'])
    print(edge_tab.to_string(index=False, float_format=lambda v: f'{v:.4f}'))
    # recorded parent baselines (cell 2 of the notebook)
    rt = ref_top[ref_top['rank'] < 10].merge(clicks[['file_name', 'seed_index', 'click_cx', 'click_cy']], on=['file_name', 'seed_index']).join(ims.set_index('file_name')[['width', 'height']], on='file_name')
    rtp = rt['bucket'] == TP
    rdx, rdy = (rt['cx'] - rt['click_cx']).abs(), (rt['cy'] - rt['click_cy']).abs()
    redge = np.minimum.reduce([rt['cx'], rt['cy'], rt['width'] - 1 - rt['cx'], rt['height'] - 1 - rt['cy']])
    base = dict(slots=len(rt), near_cheb=int((np.maximum(rdx, rdy) <= 60).sum()), near_cheb_fp=int(((np.maximum(rdx, rdy) <= 60) & ~rtp).sum()), near_euc=int((np.hypot(rdx, rdy) <= 60).sum()), near_euc_fp=int(((np.hypot(rdx, rdy) <= 60) & ~rtp).sum()), band25_pct=round(100 * (redge < 25).mean(), 2), tp_in_25=round(rtp[redge < 25].mean(), 3), tp_out_25=round(rtp[redge >= 25].mean(), 3), band35_pct=round(100 * (redge < 35).mean(), 2), tp_in_35=round(rtp[redge < 35].mean(), 3), tp_out_35=round(rtp[redge >= 35].mean(), 3), max_peaks_binding=int(ref_run['max_peaks_binding'].sum()))
    printed_base = dict(slots=23520, near_cheb=28, near_cheb_fp=0, near_euc=16, near_euc_fp=0, band25_pct=3.68, tp_in_25=0.113, tp_out_25=0.562, band35_pct=4.70, tp_in_35=0.246, tp_out_35=0.560, max_peaks_binding=2352)
    ledger('cell 2 recorded-run baselines', len(printed_base), int(sum(base[k] != v for k, v in printed_base.items())))
    print(base)

    # ----------------------------------------------------------------------------------------------- A8
    section('A8  churn fill_51 vs fill_71 (cell 16) and the pilot comparison')
    ks = {k: set(zip(g['cx'], g['cy'])) for k, g in t[t['rank'] < 10].groupby(['rank_key', 'condition', 'file_name', 'seed_index'])}
    tp10 = run.set_index(['rank_key', 'condition', 'file_name', 'seed_index'])['tp_at_10']
    ch_rows = [dict(rank_key=ax, file_name=f, seed_index=s, overlap=len(ks[(ax, 'fill_51', f, s)] & ks[(ax, 'fill_71', f, s)]), d=int(tp10[(ax, 'fill_71', f, s)] - tp10[(ax, 'fill_51', f, s)])) for ax in AXES for f in roi_order for s in SEEDS]
    churn = pd.DataFrame(ch_rows)
    out_churn = []
    for scope, fr in (('all', churn), ('pilot', churn[churn['file_name'].isin(['337.tiff', '509.tiff'])])):
        for ax in AXES:
            g = fr[fr['rank_key'] == ax]
            out_churn.append(dict(scope=scope, rank_key=ax, changed=int((g['d'] != 0).sum()), n=len(g), mean_overlap=round(g['overlap'].mean(), 2), net=int(g['d'].sum()), per_roi=dict(g.groupby('file_name')['d'].sum()) if scope == 'pilot' else ''))
    out_churn = pd.DataFrame(out_churn)
    print(out_churn.to_string(index=False))
    printed_churn = [(69, 5.49, -34), (132, 5.57, 59), (3, 5.00, 1), (9, 5.92, 4)]
    ledger('cell 16 churn numbers (4 rows x 3)', 12, int(sum(sum(int(a != b) for a, b in zip((r_.changed, r_.mean_overlap, r_.net), p_)) for r_, p_ in zip(out_churn.itertuples(), printed_churn))))
    out_churn.to_csv(f'{OUT}_churn.csv', index=False)

    # ----------------------------------------------------------------------------------------------- A9
    section('A9  template mechanics and the NMS / MAX_PEAKS table, from per_run')
    r1 = run[run['rank_key'] == 'chromatin_od'].copy()
    r2 = run[run['rank_key'] == 'tm_score'].set_index(['file_name', 'seed_index', 'condition'])
    same_cols = ['tpl_sha1', 'n_peaks', 'n_detections', 'fill_value', 'n_kept_px', 'n_filled_px', 't_search_s']
    r1i = r1.set_index(['file_name', 'seed_index', 'condition'])
    ledger('one search per (click, arm): search-level columns identical on both axes (2058 x 7)', 2058 * 7, int(sum(((r1i[c] != r2.loc[r1i.index, c]) & ~(r1i[c].isna() & r2.loc[r1i.index, c].isna())).sum() for c in same_cols)))
    nms = r1.groupby('condition').agg(max_peaks_binding=('max_peaks_binding', 'mean'), n_peaks_min=('n_peaks', 'min'), nms_mean=('n_nms_suppressed', 'mean'), nms_median=('n_nms_suppressed', 'median'), nms_min=('n_nms_suppressed', 'min'), nms_max=('n_nms_suppressed', 'max'), n_det_min=('n_detections', 'min')).reindex(ARMS)
    print(nms.to_string(float_format=lambda v: f'{v:.3g}'))
    printed_nms = {'default_51': (1, 100, 1.21, 0, 0, 18, 82), 'fill_51': (1, 100, 0.98, 0, 0, 15, 85), 'default_61': (1, 100, 1.49, 1, 0, 19, 81), 'fill_61': (1, 100, 1.01, 0, 0, 18, 82), 'default_71': (1, 100, 1.76, 1, 0, 21, 79), 'fill_71': (1, 100, 1.28, 1, 0, 17, 83), 'fill_71_pinned': (1, 100, 1.24, 1, 0, 16, 84)}
    ledger('cell 10 NMS table (7 x 7, 3 s.f.)', 49, int(sum(sum(int(float(f'{x:.3g}') != float(y)) for x, y in zip(nms.loc[a].to_numpy(float), v)) for a, v in printed_nms.items())))
    ledger('n_nms_suppressed == n_peaks - n_detections (4116)', len(run), int((run['n_nms_suppressed'] != run['n_peaks'] - run['n_detections']).sum()))
    fa = r1[r1['condition'].str.startswith('fill')].copy()
    fa['m_keep'] = fa['mass_weight'] + fa['fill_value']
    mk = fa.pivot_table(index=['file_name', 'seed_index'], columns='condition', values='m_keep')
    ledger('m_keep = mass_weight + fill_value identical across the 4 fill arms (294 x 3)', 294 * 3, int(sum((~np.isclose(mk[a], mk['fill_51'], rtol=1e-9, atol=1e-12)).sum() for a in ('fill_61', 'fill_71', 'fill_71_pinned'))))
    fv = fa.pivot_table(index=['file_name', 'seed_index'], columns='condition', values='fill_value')
    nk = fa.pivot_table(index=['file_name', 'seed_index'], columns='condition', values='n_kept_px')['fill_51']
    mw = fa.pivot_table(index=['file_name', 'seed_index'], columns='condition', values='mass_weight')
    mech = []
    for W in (51, 61, 71):
        mech.append(dict(W=W, fill_value_median=fv[f'fill_{W}'].median(), mass_weight_median=mw[f'fill_{W}'].median(), filled_share_median=((W * W - nk) / (W * W)).median()))
    mech = pd.DataFrame(mech)
    print(mech.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    mw_rows = []
    for W in (61, 71):
        c = mw[f'fill_{W}'] / mw['fill_51'] - 1
        fs = ((W * W - nk) / (W * W)) / ((51 * 51 - nk) / (51 * 51)) - 1
        nf = (W * W - nk) / (51 * 51 - nk) - 1
        mw_rows.append(dict(W=W, mass_weight_change_median_pct=round(100 * c.median(), 1), q25=round(100 * c.quantile(0.25), 1), q75=round(100 * c.quantile(0.75), 1), filled_share_change_median_pct=round(100 * fs.median(), 1), filled_count_change_median_pct=round(100 * nf.median(), 1), filled_share_abs_change_median_pts=round(100 * ((W * W - nk) / (W * W) - (51 * 51 - nk) / (51 * 51)).median(), 1)))
    mw_rows = pd.DataFrame(mw_rows)
    print(mw_rows.to_string(index=False))
    printed_mw = [(-2.3, -5.2, 0.1, 9.3, 56.3), (-4.4, -8.6, -0.6, 14.9, 122.7)]
    ledger('cell 10 mass-weight / filled-share lines (2 x 5)', 10, int(sum(sum(int(abs(a - b) > 1e-9) for a, b in zip((r_.mass_weight_change_median_pct, r_.q25, r_.q75, r_.filled_share_change_median_pct, r_.filled_count_change_median_pct), p_)) for r_, p_ in zip(mw_rows.itertuples(), printed_mw))))
    printed_mech = [(0.0281, 0.0734, 0.7645), (0.03049, 0.07162, 0.8354), (0.03184, 0.06945, 0.8785)]
    ledger('cell 10 template-mechanics medians (3 x 3, 4 s.f.)', 9, int(sum(sum(int(float(f'{a:.4g}') != b) for a, b in zip((r_.fill_value_median, r_.mass_weight_median, r_.filled_share_median), p_)) for r_, p_ in zip(mech.itertuples(), printed_mech))))
    ledger("median keep region printed '612 px' (value 612.5, printed with :.0f)", 1, int(f'{nk.median():.0f}' != '612'))
    shas = r1.pivot_table(index=['file_name', 'seed_index'], columns='condition', values='tpl_sha1', aggfunc='first')
    ledger('all 7 templates distinct per click (294)', 294, int((shas.nunique(axis=1) != 7).sum()))
    ledger('f_71 != f_51 on every click (294)', 294, int((fv['fill_71'] == fv['fill_51']).sum()))
    ledger('fill_71_pinned fill == f_51 in float32 (294)', 294, int((fv['fill_71_pinned'].astype(np.float32) != fv['fill_51'].astype(np.float32)).sum()))
    ck = clicks.set_index(['file_name', 'seed_index'])
    ledger('fill_51 n_kept_px == clicks n_kept_px_rect; fill_51 fill == clicks fill_value (294 x 2)', 588, int((nk != ck.loc[nk.index, 'n_kept_px_rect']).sum() + (fv['fill_51'].astype(np.float32) != ck.loc[fv.index, 'fill_value'].astype(np.float32)).sum()))
    mech.to_csv(f'{OUT}_template_mechanics.csv', index=False)

    # ----------------------------------------------------------------------------------------------- A10
    section('A10  Q3: nucleus size, recomputed (Spearman from ranks; bootstrap streams 100 + j and 200 + j as documented)')
    clicks['nucleus_px'] = np.maximum(clicks['rect_h'], clicks['rect_w'])
    rs = clicks.groupby('file_name', sort=False).agg(nucleus_px=('nucleus_px', 'median'), domain=('domain', 'first'))
    rs['canine'] = rs['domain'].str.startswith('canine').astype(int)
    q3 = []
    rho, p = spearman(rs['nucleus_px'], rs['canine'])
    q3.append(dict(item='Spearman(ROI nucleus px, canine) over 49', rho=rho, p=p, ci_low=np.nan, ci_high=np.nan, n=49))
    ref_roi = ref_run.pivot_table(index=['rank_key', 'file_name'], columns='condition', values='tp_at_10', aggfunc='sum')
    for ax in AXES:
        d = (ref_roi.loc[ax]['hem_rect_fill'] - ref_roi.loc[ax]['default_51'])
        rho, p = spearman(rs.loc[d.index, 'nucleus_px'], d)
        q3.append(dict(item=f'parent Spearman(nucleus, hem_rect_fill - default_51) {ax}', rho=rho, p=p, ci_low=np.nan, ci_high=np.nan, n=49))

    def cboot(groups, stream):
        idx = np.random.default_rng([BOOT_SEED, stream]).integers(len(groups), size=(N_BOOT, len(groups)))
        rh = []
        for row in idx:
            x = np.concatenate([groups[j][0] for j in row])
            y = np.concatenate([groups[j][1] for j in row])
            rh.append(spearman(x, y)[0])
        rh = np.asarray(rh)
        rh = rh[np.isfinite(rh)]
        return float(np.percentile(rh, 2.5)), float(np.percentile(rh, 97.5))

    q3_plot = {}
    for j, ax in enumerate(AXES):
        tab = roi10['tp_at_10'].loc[ax].reindex(roi_order)
        dl = (100 * (tab['fill_71'] - tab['fill_51']) / 60.0)
        x = rs.loc[roi_order, 'nucleus_px'].to_numpy()
        rho, p = spearman(x, dl.to_numpy())
        lo, hi = cboot([(np.array([xi]), np.array([yi])) for xi, yi in zip(x, dl.to_numpy())], 100 + j)
        q3.append(dict(item=f'between ROIs {ax}', rho=rho, p=p, ci_low=lo, ci_high=hi, n=49))
        q3_plot[ax] = pd.DataFrame({'file_name': roi_order, 'nucleus_px': x, 'delta_pts': dl.to_numpy(), 'canine': rs.loc[roi_order, 'canine'].to_numpy()})
        ck_ = clicks[['file_name', 'seed_index', 'nucleus_px']].copy()
        ck_['d'] = [int(tp10[(ax, 'fill_71', f, s)] - tp10[(ax, 'fill_51', f, s)]) for f, s in zip(ck_['file_name'], ck_['seed_index'])]
        for c in ('nucleus_px', 'd'):
            ck_[c + '_dm'] = ck_[c] - ck_.groupby('file_name')[c].transform('mean')
        rho, p = spearman(ck_['nucleus_px_dm'], ck_['d_dm'])
        groups = [(g['nucleus_px_dm'].to_numpy(), g['d_dm'].to_numpy()) for _, g in ck_.groupby('file_name', sort=False)]
        lo, hi = cboot(groups, 200 + j)
        q3.append(dict(item=f'within ROIs {ax}', rho=rho, p=p, ci_low=lo, ci_high=hi, n=294))
        for sp_, g in q3_plot[ax].groupby('canine'):
            q3.append(dict(item=f'{ax} {"canine" if sp_ else "human"} ROIs: mean fill_71 - fill_51 pts / median nucleus px', rho=g['delta_pts'].mean(), p=g['nucleus_px'].median(), ci_low=np.nan, ci_high=np.nan, n=len(g)))
    q3 = pd.DataFrame(q3)
    print(q3.to_string(index=False, float_format=lambda v: f'{v:+.4f}'))
    printed_q3 = [(-0.692, 3.7e-08, None, None), (+0.158, 0.28, None, None), (-0.322, 0.024, None, None), (+0.374, 0.008, 0.082, 0.614), (-9.7e-05, 0.999, -0.097, 0.097), None, None, (+0.496, 0.00029, 0.273, 0.663), (-0.183, 0.002, -0.316, -0.042), None, None]

    def tol(v):
        return 0.0005 if abs(v) >= 1e-3 else 0.5 * 10 ** (math.floor(math.log10(abs(v))) - 1)   # the notebook prints '+.3f' at >= 1e-3, else '.2g'

    n_q, n_qbad = 0, 0
    for ri, (r_, p_) in enumerate(zip(q3.itertuples(), printed_q3)):
        if p_ is None:
            continue
        for ci, (val, pv_) in enumerate(zip((r_.rho, r_.p, r_.ci_low, r_.ci_high), p_)):
            if pv_ is None:
                continue
            n_q += 1
            tl = (0.5 * 10 ** (math.floor(math.log10(abs(pv_))) - 1)) if (ri < 3 and ci == 1) else tol(pv_)   # the fact lines print p at '.2g'
            n_qbad += int(abs(val - pv_) > tl)
    ledger('Q3 printed rho / p / CI (cell 18)', n_q, n_qbad, 'at printed precision')
    sp_means = q3[q3['item'].str.contains('ROIs: mean')]
    printed_sp = [(-0.16, 33, 21), (-1.90, 25, 28), (3.10, 33, 21), (1.19, 25, 28)]   # q3 rows: per axis, human (canine = 0) first
    ledger('Q3 species split lines (4 x 3)', 12, int(sum(sum(int(a != b) for a, b in zip((round(r_.rho, 2), r_.p, r_.n), p_)) for r_, p_ in zip(sp_means.itertuples(), printed_sp))))
    q3.to_csv(f'{OUT}_q3.csv', index=False)

    # ----------------------------------------------------------------------------------------------- A11
    section('A11  beyond the notebook: LOO-ROI, worst click, zero lists, domain as unit, fill_71 - fill_61, od across pads')
    loo = []
    for ax in AXES:
        tab = roi10['tp_at_10'].loc[ax].reindex(roi_order)
        for cname in list(TESTED) + list(EXTRA_CONTRASTS):
            coefs = {**CONTRASTS, **EXTRA_CONTRASTS}[cname]
            d = sum(c * tab[a] for a, c in coefs.items()).to_numpy().astype(np.int64)
            vals, ps, holms = [], [], []
            fam = [k for k, v in TESTED.items() if cname in TESTED and v == TESTED[cname]]
            others = [sum(c * tab[a] for a, c in CONTRASTS[o].items()).to_numpy().astype(np.int64) for o in fam if o != cname]
            for i in range(49):
                dd = np.delete(d, i)
                vals.append(dd.sum() / (6 * 48 * 10))
                ps.append(signflip_p(dd)[0])
                if others:
                    holms.append(holm([ps[-1]] + [signflip_p(np.delete(o, i))[0] for o in others])[0])
            sgn = np.sign(d.sum()) if d.sum() else 1
            ordd = np.argsort(-sgn * d, kind='stable')
            imax = int(ordd[0])
            loo.append(dict(rank_key=ax, contrast=cname, pooled_pts=100 * d.sum() / 2940, loo_min_pts=100 * min(vals), loo_max_pts=100 * max(vals), loo_max_p=max(ps), loo_max_holm_p=max(holms) if holms else np.nan, loo_argmax_holm_roi=roi_order[int(np.argmax(holms))] if holms else '', largest_roi=roi_order[imax], largest_roi_delta_tp=int(d[imax]), largest_roi_share_of_net=float(d[imax] / d.sum()) if d.sum() else np.nan, top3_share_of_net=float(d[ordd[:3]].sum() / d.sum()) if d.sum() else np.nan))
    loo = pd.DataFrame(loo)
    print(loo.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    loo.to_csv(f'{OUT}_loo_roi.csv', index=False)
    # worst click and zero lists
    wc = []
    for ax in AXES:
        for cname in ('fill_61 - fill_51', 'fill_71 - fill_51'):
            coefs = CONTRASTS[cname]
            dc = sum(c * seed_tab['tp_at_10'].loc[ax][a] for a, c in coefs.items())
            wc.append(dict(rank_key=ax, contrast=cname, clicks_up=int((dc > 0).sum()), clicks_down=int((dc < 0).sum()), clicks_down_ge2=int((dc <= -2).sum()), clicks_up_ge2=int((dc >= 2).sum()), worst_click=int(dc.min()), best_click=int(dc.max()), worst_click_id=str(dc.idxmin())))
    wc = pd.DataFrame(wc)
    print(wc.to_string(index=False))
    zero = run.assign(zero=run['tp_at_10'] == 0).groupby(['rank_key', 'condition'])['zero'].sum().unstack('condition')[list(ARMS)]
    print('\nlists with TP@10 = 0 (of 294):')
    print(zero.to_string())
    # per-click worst TP@10 by arm (a reading-burden view: 10 - TP@10 false positives to read through)
    lowest = run.groupby(['rank_key', 'condition'])['tp_at_10'].describe()[['min', '25%', '50%', 'mean']]
    print(lowest.to_string())
    wc.to_csv(f'{OUT}_worst_click.csv', index=False)
    zero.reset_index().to_csv(f'{OUT}_zero_lists.csv', index=False)
    # domain as unit: G = 7 domain sums, exact sign flip by full enumeration (2^7)
    dom_rows = []
    roi_dom = run.drop_duplicates('file_name').set_index('file_name')['domain']
    for ax in AXES:
        tab = roi10['tp_at_10'].loc[ax].reindex(roi_order)
        for cname in list(TESTED) + ['default_61 - default_51', 'default_71 - default_51']:
            d = sum(c * tab[a] for a, c in CONTRASTS[cname].items())
            dd = d.groupby(roi_dom).sum()
            signs = np.array(np.meshgrid(*[[1, -1]] * 7)).reshape(7, -1).T
            null = np.abs(signs @ dd.to_numpy())
            p = float((null >= abs(dd.sum())).mean())
            dom_rows.append(dict(rank_key=ax, contrast=cname, domains_positive=int((dd > 0).sum()), domains_negative=int((dd < 0).sum()), domains_zero=int((dd == 0).sum()), domain_signflip_p=p, floor=2.0 ** (1 - int((dd != 0).sum())), **{k: round(100 * v / 420, 2) for k, v in dd.items()}))
    dom_rows = pd.DataFrame(dom_rows)
    print(dom_rows.to_string(index=False))
    dom_rows.to_csv(f'{OUT}_domain_as_unit.csv', index=False)
    # top-10 composition (what the gain is made of) and per-click-draw consistency of Q1
    comp10 = top[top['rank'] < 10].groupby(['rank_key', 'condition'])['my_bucket'].value_counts().unstack(fill_value=0)
    comp10 = (comp10 / 2940 * 100).round(2)
    print('top-10 composition, % of 2940 slots (human_correct = TP, human_rejected = look-alike, non_human = unannotated):')
    print(comp10.to_string())
    comp10.reset_index().to_csv(f'{OUT}_top10_composition.csv', index=False)
    by_seed = []
    for ax in AXES:
        st_ = seed_tab['tp_at_10'].loc[ax]
        for cname in ('fill_61 - fill_51', 'fill_71 - fill_51'):
            dc = sum(c * st_[a] for a, c in CONTRASTS[cname].items())
            g_ = dc.groupby(level='seed_index').sum()
            by_seed.append(dict(rank_key=ax, contrast=cname, **{f's{k}': int(v) for k, v in g_.items()}, draws_positive=int((g_ > 0).sum())))
    by_seed = pd.DataFrame(by_seed)
    print('net TP@10 change by click draw (seed_index), summed over 49 ROIs:')
    print(by_seed.to_string(index=False))
    by_seed.to_csv(f'{OUT}_by_click_draw.csv', index=False)
    # fill_71 - fill_61, descriptive
    ex = []
    for ax in AXES:
        tab = roi10['tp_at_10'].loc[ax].reindex(roi_order)
        d = (tab['fill_71'] - tab['fill_61']).to_numpy().astype(np.int64)
        ex.append(dict(rank_key=ax, contrast='fill_71 - fill_61', **contrast_stats(d, 10, idx49)))
    ex = pd.DataFrame(ex)
    print(ex[['rank_key', 'contrast', 'pooled_delta_precision', 'ci_low', 'ci_high', 'wins', 'losses', 'ties', 'exact_p']].to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    ex.to_csv(f'{OUT}_fill71_vs_fill61.csv', index=False)
    # od identical across pads / arms at the same (cx, cy) within a click
    g = top.groupby(['file_name', 'seed_index', 'cx', 'cy'])['od']
    nuniq = g.nunique()
    ledger('od identical for the same (click, cx, cy) across all 7 arms and both axes (groups)', len(nuniq), int((nuniq > 1).sum()), f'{int((g.size() > 1).sum())} groups carry >1 row')
    # t_search timing
    tim = r1.groupby('W')['t_search_s'].agg(['median', 'mean', 'sum', 'count'])
    print('\nper-search time (s) by W, from per_run (one row per search):')
    print(tim.to_string())
    ledger('printed per-search medians 2.27 / 2.52 / 1.96', 3, int(sum(round(tim.loc[W, 'median'], 2) != v for W, v in ((51, 2.27), (61, 2.52), (71, 1.96)))))
    print(f"sum of t_search_s over 2058 searches: {tim['sum'].sum() / 60:.1f} min (loop printed 91.7 min)")

    # ----------------------------------------------------------------------------------------------- A13
    section('A13  in-sample selection: leave-one-ROI-out choice of configuration (Step 4.10)')
    cfg_tab = {}
    for ax in AXES:
        par = ref_roi10_all = ref_run.pivot_table(index=['rank_key', 'file_name'], columns='condition', values='tp_at_10', aggfunc='sum').loc[ax].reindex(roi_order)
        cur = roi10['tp_at_10'].loc[ax].reindex(roi_order)
        assert (par['default_51'] == cur['default_51']).all() and (par['hem_rect_fill'] == cur['fill_51']).all()
        cfg_tab[ax] = pd.concat([par[['default_51', 'hem_rect_fill', 'hem_rect_fill_1um', 'hem_mask_fill_1um']], cur[['default_61', 'default_71', 'fill_61', 'fill_71', 'fill_71_pinned']]], axis=1)
    sel_rows = []
    for ax in AXES:
        C = cfg_tab[ax]
        for scope, cols, base in (('all 9 configurations across the two runs, vs default_51', list(C.columns), 'default_51'), ('fill window only {hem_rect_fill, fill_61, fill_71}, vs fill_51', ['hem_rect_fill', 'fill_61', 'fill_71'], 'hem_rect_fill')):
            ins_best = C[cols].sum().idxmax()
            ins_gain = (C[ins_best] - C[base]).sum() / 2940
            held, chosen = [], []
            for i, f in enumerate(roi_order):
                rest = C[cols].drop(index=f).sum()
                pick = rest.idxmax()
                chosen.append(pick)
                held.append(int(C.loc[f, pick] - C.loc[f, base]))
            held = np.asarray(held)
            p_, _ = signflip_p(held)
            sel_rows.append(dict(rank_key=ax, scope=scope, n_configs=len(cols), in_sample_best=ins_best, in_sample_gain_pts=100 * ins_gain, loo_gain_pts=100 * held.sum() / 2940, loo_exact_p=p_, loo_choices=';'.join(f'{k}:{v}' for k, v in pd.Series(chosen).value_counts().items())))
    sel = pd.DataFrame(sel_rows)
    print(sel.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    sel.to_csv(f'{OUT}_loo_selection.csv', index=False)
    # descriptive: the untested default_51 -> fill_61 step the notebook quotes
    tab = roi10['tp_at_10'].loc['tm_score'].reindex(roi_order)
    st = contrast_stats((tab['fill_61'] - tab['default_51']).to_numpy().astype(np.int64), 10, idx49)
    print(f"descriptive fill_61 - default_51 at tm_score: {100 * st['pooled_delta_precision']:+.2f} pts [{100 * st['ci_low']:+.2f}, {100 * st['ci_high']:+.2f}], W/L/T {st['wins']}/{st['losses']}/{st['ties']}, exact p {st['exact_p']:.3g}")

    # ----------------------------------------------------------------------------------------------- A14
    section('A14  configuration drift: notebook vs production.py vs FSConfig defaults vs decisions (source parsed as text)')
    import re
    prod_src = open(f'{REPO}/midog_utils/production.py').read()
    fs_src = open(f'{REPO}/midog_utils/find_and_suppress.py').read()
    ev_src = open(f'{REPO}/midog_utils/evaluate.py').read()

    def grab(src, pat):
        m = re.search(pat, src)
        return m.group(1).strip() if m else 'NOT FOUND'

    drift = pd.DataFrame([
        dict(setting='channel', notebook='hematoxylin_od', production_py=grab(prod_src, r'CHANNEL = "([^"]+)"'), fsconfig_default='(caller supplies channel)', decision='D3 hematoxylin_od'),
        dict(setting='tm_method', notebook='cv2.TM_CCOEFF (4)', production_py=grab(prod_src, r'TM_METHOD = (\S+)'), fsconfig_default=grab(fs_src, r'tm_method: int = (\S+)'), decision='D1 TM_CCOEFF'),
        dict(setting='peak_min_distance', notebook='7', production_py=grab(prod_src, r'PEAK_MIN_DISTANCE = (\S+)'), fsconfig_default=grab(fs_src, r'peak_min_distance: int = (\S+)'), decision='-'),
        dict(setting='deep_floor_z', notebook='-1.5', production_py=grab(prod_src, r'DEEP_FLOOR_Z = (\S+)'), fsconfig_default=grab(fs_src, r'deep_floor_z: Optional\[float\] = (\S+)'), decision='-'),
        dict(setting='max_peaks', notebook='100', production_py=grab(prod_src, r'MAX_PEAKS = (\S+)'), fsconfig_default=grab(fs_src, r'max_peaks: int = (\S+)'), decision='D9 100'),
        dict(setting='nms_radius', notebook='ev.radius_px(mpp) = 7.5 um', production_py='ev.radius_px(mpp)', fsconfig_default=grab(fs_src, r'nms_radius: float = (\S+)'), decision='D7 7.5 um = ' + grab(ev_src, r'MIDOG_RADIUS_UM = (\S+)')),
        dict(setting='od window', notebook='51', production_py=grab(prod_src, r'OD_WINDOW = (\S+)'), fsconfig_default='-', decision='-'),
        dict(setting='blanked square', notebook='51 for every arm, W = 51 / 61 / 71', production_py='seed.base_size (the template side)', fsconfig_default='cfg.base_size', decision='D11: the template footprint (base_size x base_size); diverges for W = 61 / 71'),
        dict(setting='border padding', notebook='replicate, (W - 1) // 2', production_py='border_pad=True', fsconfig_default=grab(fs_src, r'border_pad: bool = (\S+)'), decision='-'),
    ])
    print(drift.to_string(index=False))
    drift.to_csv(f'{OUT}_config_drift.csv', index=False)

    # ----------------------------------------------------------------------------------------------- A12
    section('A12  figures: the plotted quantities behind both figures')
    fig1 = mine[mine['metric'] == 'tp_at_10'].copy()
    fig1['label'] = [f"{100 * r_.pooled_delta_precision:+.2f}  W/L/T {r_.wins}/{r_.losses}/{r_.ties}" + (f"  Holm {r_.holm_p:.2g}" if r_.role == 'tested' else '') for r_ in fig1.itertuples()]
    print(fig1[['rank_key', 'contrast', 'label', 'ci_low', 'ci_high']].to_string(index=False))
    fig1.to_csv(f'{OUT}_figure1_series.csv', index=False)
    for ax in AXES:
        g_ = fig1[fig1['rank_key'] == ax]
        print(f"fig 1 {ax}: data span {100 * min(g_['ci_low'].min(), g_['pooled_delta_precision'].min()):+.2f} .. {100 * max(g_['ci_high'].max(), g_['pooled_delta_precision'].max()):+.2f} points over {len(g_)} rows")
    for ax in AXES:
        f2 = q3_plot[ax]
        print(f"fig 2 {ax}: {len(f2)} points ({int(f2['canine'].sum())} canine), y range {f2['delta_pts'].min():+.2f} .. {f2['delta_pts'].max():+.2f}, x range {f2['nucleus_px'].min()} .. {f2['nucleus_px'].max()}, points at y = 0: {int((f2['delta_pts'] == 0).sum())}")
        f2.assign(rank_key=ax).to_csv(f'{OUT}_figure2_series_{ax}.csv', index=False)

    # ----------------------------------------------------------------------------------------------- TIER B
    if '--skip-tierb' in sys.argv:
        led = pd.DataFrame(LEDGER)
        print(led.to_string(index=False))
        print(f"\nTOTAL (Tier A only): {led['n_compared'].sum():,} compared, {led['n_diverge'].sum():,} diverge")
        return
    section('TIER B  templates and searches from pixels, 49 ROIs, multiprocessing')
    need_default = set(map(tuple, sm[['file_name', 'seed_index']].drop_duplicates().to_numpy().tolist())) | set(map(tuple, near_gt[near_gt['n_mitotic_within_2r'] > 0][['file_name', 'seed_index']].to_numpy().tolist()))
    jobs = []
    tb_rois = roi_order if '--tierb-rois' not in sys.argv else [x for x in roi_order if x in sys.argv[sys.argv.index('--tierb-rois') + 1].split(',')]
    for f in tb_rois:
        cr = clicks[clicks['file_name'] == f].sort_values('seed_index')
        jobs.append((f, cr['subset'].iat[0], cr.to_dict('records'), gt_by_file[f], int(ims.set_index('file_name').loc[f, 'width']), int(ims.set_index('file_name').loc[f, 'height']), {s for (ff, s) in need_default if ff == f}, {s for (ff, s) in TIER_B_CLICKS_ALL7 if ff == f}))
    print(f'default_61/71 re-run on {len(need_default)} clicks; all 7 arms on {TIER_B_CLICKS_ALL7}')
    tb0 = time.time()
    with get_context('spawn').Pool(N_WORKERS) as pool:
        res = []
        for i, r_ in enumerate(pool.imap_unordered(tierb_roi, jobs)):
            res.append(r_)
            print(f'  [{i + 1}/49] {r_["file_name"]} {r_["t_total"]:.0f}s (hem {r_["t_hem"]:.1f}s)  elapsed {(time.time() - tb0) / 60:.1f} min', flush=True)
    TB_WALL = (time.time() - tb0) / 60
    b_order = [f for f in roi_order if f in set(tb_rois)]
    idxb = idx49 if len(b_order) == 49 else boot_idx(len(b_order), 0)
    tpl_b = pd.DataFrame([x for r_ in res for x in r_['tpl']])
    srch_b = pd.DataFrame([x for r_ in res for x in r_['srch']])
    top_b = pd.concat([r_['top'] for r_ in res], ignore_index=True)
    ident_b = pd.DataFrame([x for r_ in res for x in r_['identity']])
    selfm_b = pd.DataFrame([x for r_ in res for x in r_['selfmatch']])
    mt = {W: [x for r_ in res for x in r_['t_mt'][W]] for W in (51, 61, 71)}
    print(f'Tier B wall clock {TB_WALL:.1f} min; CPU-summed ROI time {sum(r_["t_total"] for r_ in res) / 60:.1f} min; matchTemplate median by W (s, 4 workers x 3 threads): ' + ', '.join(f'{W}: {np.median(v):.2f} (n={len(v)})' for W, v in mt.items() if v))
    ledger(f'Tier B: mpp from TIFF tags == per_run mpp ({len(res)})', len(res), int(sum(abs(r_['mpp'] - float(clicks[clicks['file_name'] == r_['file_name']]['mpp'].iat[0])) > 1e-15 for r_ in res)))
    # B1 templates
    rt_ = r1[r1['file_name'].isin(tb_rois)].set_index(['file_name', 'seed_index', 'condition'])
    tb = tpl_b.set_index(['file_name', 'seed_index', 'condition']).loc[rt_.index]
    ledger('B1 template SHA-1, all 7 arms x 294 clicks, rebuilt from pixels', len(tb), int((tb['sha'] != rt_['tpl_sha1']).sum()))
    fv_ok = ((tb['fill_value'].isna() & rt_['fill_value'].isna()) | (tb['fill_value'].astype(np.float32) == rt_['fill_value'].astype(np.float32)))
    ledger('B1 fill_value (float32), 7 x 294', len(tb), int((~fv_ok).sum()))
    ledger('B1 n_kept_px / n_filled_px, 7 x 294 x 2', 2 * len(tb), int((tb['n_kept'] != rt_['n_kept_px']).sum() + (tb['n_filled'] != rt_['n_filled_px']).sum()))
    # B2 reproduction of searches
    rep = srch_b[srch_b['variant'] == 'repro_blank51'].merge(run, on=['file_name', 'seed_index', 'condition', 'rank_key'], suffixes=('_b', ''))
    ledger('B2 search reproduction: n_detections, tp_at_10/20/30, near-seed (per list)', 6 * len(rep), int(sum((rep[c + '_b'] != rep[c]).sum() for c in ('n_detections', 'tp_at_10', 'tp_at_20', 'tp_at_30', 'n_near_seed_top10_cheb', 'n_near_seed_fp_top10_cheb'))), f'{len(rep)} lists: ' + ', '.join(f'{k}={v}' for k, v in rep.groupby('condition').size().items()))
    tb_rep = top_b[top_b['variant'] == 'repro_blank51'].sort_values(['file_name', 'seed_index', 'condition', 'rank_key', 'rank']).reset_index(drop=True)
    keys = tb_rep[['file_name', 'seed_index', 'condition', 'rank_key']].drop_duplicates()
    tp_rep = top.merge(keys, on=['file_name', 'seed_index', 'condition', 'rank_key']).sort_values(['file_name', 'seed_index', 'condition', 'rank_key', 'rank']).reset_index(drop=True)
    assert len(tp_rep) == len(tb_rep)
    for c in ('cx', 'cy', 'bucket', 'matched_ann_id'):
        ledger(f'B2 top-30 {c} reproduced from pixels', len(tb_rep), int((tb_rep[c].to_numpy() != tp_rep[c].to_numpy()).sum()))
    ledger('B2 top-30 score, float32 exact', len(tb_rep), int((tb_rep['score'].to_numpy(np.float32) != tp_rep['score'].to_numpy().astype(np.float32)).sum()))
    ledger('B2 top-30 od, exact', len(tb_rep), int((tb_rep['od'].to_numpy() != tp_rep['od'].to_numpy()).sum()))
    # B3 identity and closed form
    if len(ident_b):
        ident_b['rel_err_sum_keep'] = (ident_b['formula_sum_keep'] - ident_b['cv2_response']).abs() / ident_b['cv2_response'].abs().clip(lower=1e-9)
        ident_b['rel_err_decomposed'] = (ident_b['formula_decomposed'] - ident_b['cv2_response']).abs() / ident_b['cv2_response'].abs().clip(lower=1e-9)
        print(f"B3 R_W identity (fill_71, {len(ident_b)} positions): max relative error sum_keep form {ident_b['rel_err_sum_keep'].max():.2e}, decomposed form {ident_b['rel_err_decomposed'].max():.2e}; response at the click, min {ident_b.loc[ident_b.groupby(['file_name', 'seed_index']).tail(1).index, 'cv2_response'].max():+.4f} (max over clicks; must be < 0)")
        at_click = ident_b.loc[ident_b.groupby(['file_name', 'seed_index']).tail(1).index, 'cv2_response']
        print(f"B3 fill_71 response at the click negative on {int((at_click < 0).sum())}/{len(at_click)} clicks (max {at_click.max():+.4f})")
        ident_b.to_csv(f'{OUT}_tierb_identity.csv', index=False)
    if len(selfm_b):
        selfm_b['abs_err'] = (selfm_b['cv2_response_at_click'] - selfm_b['closed_form_ring']).abs()
        selfm_b['rel_err'] = selfm_b['abs_err'] / selfm_b['closed_form_ring'].abs().clip(lower=1e-9)
        print(f"B3 default_71 at the click: cv2 response vs sum_ring T'(T - min), {len(selfm_b)} clicks, max relative error {selfm_b['rel_err'].max():.2e}; response at the click positive on {int((selfm_b['closed_form_ring'] > 0).sum())}/{len(selfm_b)}")
        selfm_b.to_csv(f'{OUT}_tierb_selfmatch_closed_form.csv', index=False)
    # B4 exact self-match removal on default_61 / default_71
    sr = srch_b[srch_b['variant'] == 'selfmatch_removed'].merge(run[['file_name', 'seed_index', 'condition', 'rank_key', 'tp_at_10']], on=['file_name', 'seed_index', 'condition', 'rank_key'], suffixes=('_removed', ''))
    sr['gain'] = sr['tp_at_10_removed'] - sr['tp_at_10']
    print('\nB4 self-match peak removed BEFORE the 100-peak cap and NMS, default_61 / default_71, re-run clicks:')
    print(sr.groupby(['condition', 'rank_key']).agg(n_lists=('gain', 'size'), tp10_gain=('gain', 'sum'), lists_gain=('gain', lambda g: int((g > 0).sum())), lists_loss=('gain', lambda g: int((g < 0).sum())), selfmatch_left=('n_selfmatch_top30', 'sum'), near_fp_left=('n_near_seed_fp_top10_cheb', 'sum')).to_string())
    exact_rows = []
    for ax in AXES:
        piv0 = roi10['tp_at_10'].loc[ax].reindex(b_order).copy()
        g_ = sr[sr['rank_key'] == ax].groupby(['condition', 'file_name'])['gain'].sum()
        corr_tab = piv0.copy()
        for arm in ('default_61', 'default_71'):
            if arm in g_.index.get_level_values(0):
                corr_tab[arm] = piv0[arm] + g_.loc[arm].reindex(b_order).fillna(0).astype(int)
        for cname in ('default_61 - default_51', 'default_71 - default_51', '(fill_61 - fill_51) - (default_61 - default_51)', '(fill_71 - fill_51) - (default_71 - default_51)'):
            d = sum(c * corr_tab[a] for a, c in CONTRASTS[cname].items()).to_numpy().astype(np.int64)
            st = contrast_stats(d, 10, idxb)
            exact_rows.append(dict(rank_key=ax, contrast=cname, version='self-match peak removed pre-NMS (Tier B)', **{k: st[k] for k in ('pooled_delta_precision', 'ci_low', 'ci_high', 'wins', 'losses', 'ties', 'exact_p', 'd5_bar_met')}))
    exact_rows = pd.DataFrame(exact_rows)
    print(exact_rows.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    pd.concat([sb, exact_rows], ignore_index=True).to_csv(f'{OUT}_selfmatch_bound.csv', index=False)
    sr.to_csv(f'{OUT}_tierb_selfmatch_removed.csv', index=False)
    # B5 D11-conformant blank for the fill arms
    d11 = srch_b[srch_b['variant'].str.startswith('d11')]
    rep_f = srch_b[(srch_b['variant'] == 'repro_blank51') & srch_b['condition'].isin(['fill_61', 'fill_71'])]
    d11_rows = []
    for ax in AXES:
        f51 = roi10['tp_at_10'].loc[ax].reindex(b_order)['fill_51']
        for arm in ('fill_61', 'fill_71'):
            v = d11[(d11['rank_key'] == ax) & (d11['condition'] == arm)].groupby('file_name')['tp_at_10'].sum().reindex(b_order)
            v0 = rep_f[(rep_f['rank_key'] == ax) & (rep_f['condition'] == arm)].groupby('file_name')['tp_at_10'].sum().reindex(b_order)
            for label, tab in (('blank 51 (as run, reproduced)', v0), (f'blank {ARM_W[arm]} (D11: blank = template footprint)', v)):
                st = contrast_stats((tab - f51).to_numpy().astype(np.int64), 10, idxb)
                src = d11 if label.startswith('blank ' + str(ARM_W[arm])) else rep_f
                sel_ = src[(src['rank_key'] == ax) & (src['condition'] == arm)]
                d11_rows.append(dict(rank_key=ax, contrast=f'{arm} - fill_51', version=label, near_px_reach=int(sel_['near_px'].iat[0]), near_seed_all_top10_reach=int(sel_['n_near_seed_top10_cheb_reach'].sum()), near_seed_fp_top10_reach=int(sel_['n_near_seed_fp_top10_cheb_reach'].sum()), near_seed_fp_top10_60px=int(sel_['n_near_seed_fp_top10_cheb'].sum()), **{k: st[k] for k in ('pooled_delta_precision', 'ci_low', 'ci_high', 't_ci_low', 't_ci_high', 'wins', 'losses', 'ties', 'n_nonzero', 'exact_p', 'd5_bar_met')}))
    d11_rows = pd.DataFrame(d11_rows)
    print('\nB5 fill arms with a D11-conformant blank (blank side = W) vs the run as recorded (blank 51); fill_51 unchanged (already blank = 51):')
    print(d11_rows.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    d11_rows.to_csv(f'{OUT}_tierb_d11_blank.csv', index=False)
    # validation of the B5 shortcut on one click: patched response vs a full matchTemplate on the W-blanked channel
    import cv2
    cv2.setNumThreads(6)
    vr = clicks[clicks['file_name'] == '509.tiff'].iloc[0]
    vhem = _load_hem(f"{IMG_DIRS[vr['subset']]}/509.tiff")
    vmin = float(vhem.min())
    vt = _templates(vhem, vr['click_cx'], vr['click_cy'], (int(vr['rect_y0']), int(vr['rect_y1']), int(vr['rect_x0']), int(vr['rect_x1'])))
    for W in (61, 71):
        T = vt[f'fill_{W}'][0]
        R51 = _response(_blank(vhem, vr['click_cx'], vr['click_cy'], 51, vmin), T)
        Rfull = _response(_blank(vhem, vr['click_cx'], vr['click_cy'], W, vmin), T)
        dl, x0, y0 = _local_delta(vhem, vr['click_cx'], vr['click_cy'], T, 51, W, vmin)
        R2 = R51.copy()
        hh, ww = vhem.shape
        ys0, xs0, ys1, xs1 = max(0, y0), max(0, x0), min(hh, y0 + dl.shape[0]), min(ww, x0 + dl.shape[1])
        R2[ys0:ys1, xs0:xs1] = (R2[ys0:ys1, xs0:xs1].astype(np.float64) + dl[ys0 - y0:ys1 - y0, xs0 - x0:xs1 - x0]).astype(np.float32)
        dif = np.abs(R2.astype(np.float64) - Rfull)
        c1, _ = _peaks(R2)
        c2, _ = _peaks(Rfull)
        print(f'B5 shortcut check, 509.tiff s{int(vr["seed_index"])} fill_{W}: max |patched - full recompute| = {dif.max():.2e} (response scale {np.abs(Rfull).max():.2f}); top-100 peaks identical and in the same order: {np.array_equal(c1, c2)}')
    del vhem
    dd11 = d11.merge(rep_f, on=['file_name', 'seed_index', 'condition', 'rank_key'], suffixes=('_d11', '_51'))
    print(f"lists whose TP@10 changes under the D11 blank: {int((dd11['tp_at_10_d11'] != dd11['tp_at_10_51']).sum())}/{len(dd11)}; net TP@10 change {int((dd11['tp_at_10_d11'] - dd11['tp_at_10_51']).sum()):+d}")
    srch_b.to_csv(f'{OUT}_tierb_searches.csv', index=False)
    top_b.to_csv(f'{OUT}_tierb_top30.csv', index=False)

    # ----------------------------------------------------------------------------------------------- ledger
    section('LEDGER')
    led = pd.DataFrame(LEDGER)
    print(led.to_string(index=False))
    print(f"\nTOTAL: {led['n_compared'].sum():,} values compared, {led['n_diverge'].sum():,} divergences, across {len(led)} checks")
    led.to_csv(f'{OUT}_ledger.csv', index=False)
    # re-timing (Step 4.11), after Tier B so the pool is gone; three repeats per stage, medians
    import cv2
    cv2.setNumThreads(0)
    path = f"{IMG_DIRS['original_14']}/300.tiff"
    t_h, t_m = [], {51: [], 61: [], 71: []}
    for _ in range(3):
        t0 = time.time()
        hem = _load_hem(path)
        t_h.append(time.time() - t0)
        for W in (51, 61, 71):
            T = _window(hem, 1411.0, 3551.0, W)
            t0 = time.time()
            _response(hem, T)
            t_m[W].append(time.time() - t0)
    print(f"re-timing on 300.tiff (3 repeats, medians, default cv2 threads, other processes present): load + rgb2hed {np.median(t_h):.2f} s; matchTemplate incl. replicate pad " + ', '.join(f'W={W}: {np.median(v):.2f} s' for W, v in t_m.items()))
    pd.DataFrame([dict(stage='load_tiff_plus_rgb2hed', median_s=np.median(t_h), runs=str([round(x, 2) for x in t_h]))] + [dict(stage=f'matchTemplate_W{W}', median_s=np.median(v), runs=str([round(x, 2) for x in v])) for W, v in t_m.items()]).to_csv(f'{OUT}_retiming.csv', index=False)
    # repo recounts for the agent-file appendix
    import glob
    import subprocess
    n_csv = len(glob.glob(f'{REPO}/results/*.csv'))
    untracked = [l for l in subprocess.run(['git', '-C', REPO, 'status', '--porcelain', 'results/'], capture_output=True, text=True).stdout.splitlines() if l.startswith('??')]
    n_untracked_other = sum(1 for l in untracked if 'window_sweep_51_61_71_49roi_6seed_audit' not in l)
    nbs = glob.glob(f'{REPO}/**/*.ipynb', recursive=True)
    n_png = 0
    for pth in nbs:
        try:
            nbj = json.load(open(pth))
        except Exception:
            continue
        n_png += sum(1 for c in nbj['cells'] if c['cell_type'] == 'code' for o in c.get('outputs', []) if 'image/png' in o.get('data', {}))
    print(f'repo recount: results/*.csv {n_csv}; untracked results CSVs excluding this audit {n_untracked_other}; notebooks {len(nbs)}; embedded PNGs {n_png}')
    print(f'\naudit script total {(time.time() - T_START) / 60:.1f} min')


if __name__ == '__main__':
    main()
