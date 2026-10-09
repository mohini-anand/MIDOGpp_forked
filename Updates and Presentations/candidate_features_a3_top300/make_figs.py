"""
Figures for the candidate-features deck: the feature schematics, the example-candidate strips, and the two charts.

Reads the V2 candidate table, the audit's stratified AUC table, and ROI 128's image. Writes PNGs to ./figs/, next to this
file. The deck uploads them as assets; the PNGs are not tracked.

Run with the Miniconda 3.11 interpreter: ~/anaconda3/bin/python make_figs.py
"""
import os
from pathlib import Path

import cv2
import matplotlib
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, Polygon, Rectangle
from scipy.spatial import ConvexHull

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OUT = HERE / 'figs'
OUT.mkdir(exist_ok=True)

INK, SEC, PAPER = '#1C1B22', '#55526A', '#F6F3EE'
BLUE, ORANGE, GRAY = '#2A78D6', '#EB6834', '#9A98A6'
NUC, TISSUE, GLASS = '#3B2F7A', '#E7B4C6', '#FFFFFF'
plt.rcParams['font.family'] = 'Helvetica Neue'

# the six features on the deck: (column, orientation sign, short key). Sign +1 = high values look like a mitosis.
FEATS = [
    ('dark_concavities_per_um', +1, 'ragged'),
    ('chr_solidity', -1, 'convex'),
    ('chr_elongation', +1, 'elong'),
    ('h_chr_mean', +1, 'dark'),
    ('white_frac_gray215_61', -1, 'white'),
    ('ring1_frac_above_tnuc', +1, 'rim'),
]

# example-candidate strips: one ROI, 61 px windows shown at 2.6x, six rungs spaced along the oriented ranking
ROI = '128.tiff'
WIN, SCALE, GAP = 61, 2.6, 10
TILE = int(round(WIN * SCALE))
LADDER = [0.99, 0.90, 0.72, 0.45, 0.20, 0.03]
MIN_APART_PX = 30
TIE_SEED = 7


def tile(rgb, cx, cy, is_mit):
    """
        One candidate's window, upscaled, framed green with an M badge when it is an annotated mitosis.

        rgb (np.ndarray): the ROI as uint8 RGB.
        cx, cy (float): the candidate's position in ROI pixels.
        is_mit (bool): whether the candidate matched an annotated mitosis.

        Returns np.ndarray: the uint8 RGB tile.
    """
    ix, iy = int(round(cx)), int(round(cy))
    h = WIN // 2
    t = cv2.resize(rgb[iy - h:iy + h + 1, ix - h:ix + h + 1], (TILE, TILE), interpolation=cv2.INTER_CUBIC)
    col = (0, 131, 0) if is_mit else (154, 152, 166)
    cv2.rectangle(t, (0, 0), (TILE - 1, TILE - 1), col, 7 if is_mit else 2)
    if is_mit:
        cv2.rectangle(t, (0, 0), (36, 36), col, -1)
        cv2.putText(t, 'M', (7, 28), cv2.FONT_HERSHEY_SIMPLEX, 1.0, (255, 255, 255), 2, cv2.LINE_AA)
    return t


def ladder(df, col, sign):
    """
        The rungs of one feature's strip: locations at fixed percentiles of the oriented ranking, no two within 30 px.

        df (pd.DataFrame): distinct locations with cx, cy, the feature column and a random tie key.
        col (str): the feature column.
        sign (int): +1 when high values look like a mitosis, -1 when low values do.

        Returns list[pd.Series]: one row per rung, mitosis-like end first.
    """
    d = df.assign(v=df[col] * sign).dropna(subset=['v']).sort_values(['v', 'tie'], ascending=[False, True]).reset_index(drop=True)
    out, used = [], set()
    for p in LADDER:
        i = int(round((1 - p) * (len(d) - 1)))
        for j in [i] + [i + k * s for k in range(1, 40) for s in (1, -1)]:
            if 0 <= j < len(d) and j not in used:
                r = d.iloc[j]
                if all(np.hypot(r.cx - o.cx, r.cy - o.cy) >= MIN_APART_PX for o in out):
                    out.append(r)
                    used.add(j)
                    break
    return out


def canvas():
    """
        An empty 400 x 400 px figure with a square, axis-free drawing area.

        Returns tuple[plt.Figure, plt.Axes]: the figure and its axes.
    """
    fig = plt.figure(figsize=(4, 4), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-2.2, 2.2)
    ax.set_aspect('equal')
    ax.axis('off')
    return fig, ax


def blob(r0=1.0, amp=0.0, k=9, ph=0.0, extra=None):
    """
        A closed outline with a sinusoidal radius, as a polygon.

        r0 (float): the mean radius.
        amp (float): the amplitude of the lobes.
        k (int): the number of lobes.
        ph (float): the phase of the lobes in radians.
        extra (callable or None): an added radius term, a function of the angle.

        Returns tuple[np.ndarray, np.ndarray, np.ndarray]: the angles, the radii and the (n, 2) points.
    """
    th = np.linspace(0, 2 * np.pi, 360, endpoint=False)
    r = r0 + amp * np.sin(k * th + ph)
    if extra is not None:
        r = r + extra(th)
    return th, r, np.c_[r * np.cos(th), r * np.sin(th)]


def save(fig, name):
    """
        Save a schematic with a transparent background and close it.

        fig (plt.Figure): the figure.
        name (str): the file stem under figs/.

        Returns None.
    """
    fig.savefig(OUT / f'{name}.png', dpi=100, transparent=True)
    plt.close(fig)


def strips(P):
    """
        Write the six example-candidate strips of ROI 128.

        P (pd.DataFrame): the recorded-pool, accepted rows of the candidate table.

        Returns None.
    """
    rgb = cv2.cvtColor(cv2.imread(str(REPO / 'images/extra_valid/testing_set' / ROI)), cv2.COLOR_BGR2RGB)
    loc = P[P.file_name == ROI].groupby(['cx', 'cy']).agg(is_mitosis=('is_mitosis', 'max'), **{c: (c, 'first') for c, _, _ in FEATS}).reset_index()
    loc['tie'] = np.random.default_rng(TIE_SEED).random(len(loc))
    for col, sign, key in FEATS:
        rows = ladder(loc, col, sign)
        strip = np.full((TILE, TILE * 6 + GAP * 5, 3), 246, np.uint8)
        for i, r in enumerate(rows):
            strip[:, i * (TILE + GAP):i * (TILE + GAP) + TILE] = tile(rgb, r.cx, r.cy, bool(r.is_mitosis))
        cv2.imwrite(str(OUT / f'ladder_{key}.png'), cv2.cvtColor(strip, cv2.COLOR_RGB2BGR))
        print(key, 'mitosis per rung', [int(r.is_mitosis) for r in rows])


def schematics():
    """
        Write the six feature schematics.

        Returns None.
    """
    # ragged outline: a smooth core beside a notched one, the notches marked
    fig, ax = canvas()
    th, r, xy = blob(0.78, 0, 1)
    ax.add_patch(Polygon(xy + [-1.15, 0], closed=True, fc=NUC, ec='none'))
    th, r, xy = blob(0.88, 0.14, 8, 0.4)
    xy2 = xy + [1.15, 0]
    ax.add_patch(Polygon(xy2, closed=True, fc=NUC, ec='none'))
    for i in [i for i in range(len(th)) if r[i] <= r.min() + 1e-3 or (r[i] < r[i - 1] and r[i] < r[(i + 1) % len(th)])]:
        ax.add_patch(Circle(xy2[i], 0.095, fc=ORANGE, ec=PAPER, lw=1.5, zorder=5))
    save(fig, 'sch_ragged')

    # convexity: a lobed core inside its convex hull, the gap shaded
    fig, ax = canvas()
    th, r, xy = blob(0.95, 0.30, 3, 0.3, extra=lambda t: 0.12 * np.sin(7 * t))
    xy = xy * np.array([1.25, 0.95])
    hp = xy[ConvexHull(xy).vertices]
    ax.add_patch(Polygon(hp, closed=True, fc='none', ec=ORANGE, lw=3.2, ls=(0, (4, 3)), zorder=3))
    ax.add_patch(Polygon(hp, closed=True, fc=ORANGE, ec='none', alpha=0.28, zorder=1))
    ax.add_patch(Polygon(xy, closed=True, fc=NUC, ec='none', zorder=2))
    save(fig, 'sch_convex')

    # elongation: the long and the short axis of an ellipse
    fig, ax = canvas()
    ax.add_patch(Ellipse((0, 0), 3.4, 1.35, angle=-25, fc=NUC, ec='none'))
    a = np.deg2rad(-25)
    for ln, col, ang in ((1.7, BLUE, a), (0.675, ORANGE, a + np.pi / 2)):
        d = np.array([np.cos(ang), np.sin(ang)]) * ln
        ax.add_patch(FancyArrowPatch(-d, d, arrowstyle='<->', mutation_scale=22, lw=3.4, color=col, shrinkA=0, shrinkB=0, zorder=5))
    save(fig, 'sch_elong')

    # core darkness: a pale core beside a dark one
    fig, ax = canvas()
    for x, c in ((-1.12, '#9E94CC'), (1.12, '#1F1650')):
        th, r, xy = blob(0.82, 0.04, 6, 0.2)
        ax.add_patch(Polygon(xy + [x, 0], closed=True, fc=c, ec='none'))
    save(fig, 'sch_dark')

    # white surround: the 61 px window with tissue and a glass region around a nucleus
    fig, ax = canvas()
    ax.add_patch(Rectangle((-1.9, -1.9), 3.8, 3.8, fc=TISSUE, ec=INK, lw=2.4, zorder=1))
    glass = np.array([[-1.9, 1.9], [-0.2, 1.9], [0.15, 1.2], [0.7, 0.65], [1.0, -0.1], [1.9, -0.3], [1.9, 1.9]])
    ax.add_patch(Polygon(glass, closed=True, fc=GLASS, ec='none', zorder=2))
    ax.add_patch(Rectangle((-1.9, -1.9), 3.8, 3.8, fc='none', ec=INK, lw=2.4, zorder=6))
    th, r, xy = blob(0.62, 0.04, 6, 0.3)
    ax.add_patch(Polygon(xy + [-0.35, -0.25], closed=True, fc=NUC, ec='none', zorder=4))
    save(fig, 'sch_white')

    # non-pale rim: the first micron around the core, tissue segments against pale ones
    fig, ax = canvas()
    th, r, xy = blob(0.62, 0.04, 6, 0.3)
    n, ro, ri = 120, 1.28, 0.8
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    for i in range(n):
        a0, a1 = ang[i], ang[i + 1] if i < n - 1 else 2 * np.pi
        seg = np.array([[ri * np.cos(a0), ri * np.sin(a0)], [ro * np.cos(a0), ro * np.sin(a0)], [ro * np.cos(a1), ro * np.sin(a1)], [ri * np.cos(a1), ri * np.sin(a1)]])
        c = '#F2EEF3' if np.cos(a0 - 0.5) > 0.35 else TISSUE
        ax.add_patch(Polygon(seg, closed=True, fc=c, ec=c, lw=0.6, zorder=2))
    ax.add_patch(Polygon(xy, closed=True, fc=NUC, ec='none', zorder=4))
    ax.add_patch(Circle((0, 0), ro, fc='none', ec=GRAY, lw=1.6, ls=(0, (3, 3)), zorder=5))
    save(fig, 'sch_rim')


def od_auc(P):
    """
        od's pooled within-click AUC and its ROI cluster-bootstrap 95% interval, by the audit's estimator.

        P (pd.DataFrame): the recorded-pool, accepted rows, with click, file_name, is_mitosis and od.

        Returns tuple[float, tuple[float, float]]: the AUC and the 2.5th and 97.5th percentiles over 4,000 ROI draws.
    """
    d = P[['file_name', 'click', 'is_mitosis', 'od']].dropna()
    d['rp'] = d.groupby('click')['od'].rank(method='average') * d.is_mitosis.astype(int)
    g = d.groupby('click').agg(file_name=('file_name', 'first'), s1=('rp', 'sum'), n1=('is_mitosis', 'sum'), n=('is_mitosis', 'size'))
    g['n0'] = g.n - g.n1
    g = g[(g.n1 > 0) & (g.n0 > 0)]
    g['U'] = g.s1 - g.n1 * (g.n1 + 1) / 2
    g['N'] = g.n1 * g.n0
    t = g.groupby('file_name')[['U', 'N']].sum()
    w = np.random.default_rng(20261009).multinomial(len(t), np.ones(len(t)) / len(t), size=4000)
    b = (w @ t.U.to_numpy()) / (w @ t.N.to_numpy())
    return t.U.sum() / t.N.sum(), (np.percentile(b, 2.5), np.percentile(b, 97.5))


def chart_auc(P):
    """
        The dot plot of within-click AUC, as is and at matched darkness (dark50), oriented so above 0.5 is mitosis-like.

        P (pd.DataFrame): the recorded-pool, accepted rows, for the od reference.

        Returns None.
    """
    S = pd.read_csv(REPO / 'results/candidate_features_a3_top300_audit_stratified.csv')
    S = S[S.population == 'pool']

    def get(col, ctrl):
        r = S[(S.column == col) & (S.control == ctrl)].iloc[0]
        return r.auc, r.ci_lo, r.ci_hi

    auc, (lo, hi) = od_auc(P)
    print('od auc', auc, (lo, hi))
    rows = [('od (baseline)', None, +1), ('Ragged outline', 'dark_concavities_per_um', +1), ('Convexity', 'chr_solidity', -1), ('Elongation', 'chr_elongation', +1),
            ('Core darkness', 'h_chr_mean', +1), ('White surround', 'white_frac_gray215_61', -1), ('Non-pale rim', 'ring1_frac_above_tnuc', +1)]
    fig = plt.figure(figsize=(11.2, 6.6), dpi=100)
    ax = fig.add_axes([0.30, 0.17, 0.67, 0.80])
    ys = np.arange(len(rows))[::-1]
    for y, (_, col, sg) in zip(ys, rows):
        if col is None:
            ax.plot([lo, hi], [y, y], color=GRAY, lw=3, solid_capstyle='round', zorder=2)
            ax.plot([auc], [y], 'o', ms=13, color=GRAY, mec=PAPER, mew=2, zorder=4)
            continue
        a_raw, _, _ = get(col, 'none')
        a_m, l_m, h_m = get(col, 'dark50')
        if sg < 0:
            a_raw, a_m, l_m, h_m = 1 - a_raw, 1 - a_m, 1 - h_m, 1 - l_m
        ax.plot([l_m, h_m], [y, y], color=BLUE, lw=3, solid_capstyle='round', zorder=2)
        ax.plot([a_raw], [y], 'o', ms=12, mfc=PAPER, mec=ORANGE, mew=3, zorder=3)
        ax.plot([a_m], [y], 'o', ms=13, color=BLUE, mec=PAPER, mew=2, zorder=4)
    ax.axvline(0.5, color=SEC, lw=1.6, ls=(0, (4, 3)), zorder=1)
    ax.set_xlim(0.45, 0.85)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_yticks(ys)
    ax.set_yticklabels([n + (f'  ({"high" if s > 0 else "low"})' if c else '') for n, c, s in rows], fontsize=19, color=INK)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8])
    ax.tick_params(axis='x', labelsize=18, colors=SEC, length=0, pad=10)
    ax.tick_params(axis='y', length=0, pad=12)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.spines['bottom'].set_color('#CFCBC4')
    ax.grid(axis='x', color='#E4E0D8', lw=1)
    ax.set_axisbelow(True)
    ax.set_xlabel('How often a mitosis ranks above another candidate in its own click', fontsize=18, color=SEC, labelpad=14)
    ax.text(0.505, len(rows) - 0.55, 'no separation', fontsize=17, color=SEC, va='top')
    fig.patch.set_alpha(0)
    fig.savefig(OUT / 'chart_auc.png', dpi=200, transparent=True)
    plt.close(fig)


def chart_precision():
    """
        The bars of precision gain over od (refused candidates dropped), from the earlier re-rank run's printed table.

        Source: candrefine_rerank_mo_pad61_select.ipynb, the P@10/20/30 table: od, REF (C1) and REF plus F4-F9 (C5, F5 direction -1).

        Returns None.
    """
    od = np.array([0.5888, 0.5194, 0.4523])
    ref = np.array([0.6374, 0.5495, 0.4822])
    c5 = np.array([0.6177, 0.5219, 0.4576])
    fig = plt.figure(figsize=(11.2, 6.6), dpi=100)
    ax = fig.add_axes([0.09, 0.14, 0.89, 0.82])
    xs = np.arange(3)
    for off, gain, col in ((-0.2, (ref - od) * 100, BLUE), (0.2, (c5 - od) * 100, ORANGE)):
        ax.bar(xs + off, gain, 0.36, color=col, zorder=3)
        for x, v in zip(xs + off, gain):
            ax.text(x, v + 0.12, f'+{v:.1f}', ha='center', va='bottom', fontsize=21, color=INK, fontweight='bold')
    ax.axhline(0, color=INK, lw=1.8, zorder=4)
    ax.set_xticks(xs)
    ax.set_xticklabels(['P@10', 'P@20', 'P@30'], fontsize=22, color=INK, fontweight='bold')
    ax.tick_params(axis='x', length=0, pad=14)
    ax.set_ylim(0, 6.2)
    ax.set_yticks([0, 2, 4, 6])
    ax.set_yticklabels(['0', '+2', '+4', '+6'], fontsize=19, color=SEC)
    ax.tick_params(axis='y', length=0)
    ax.set_ylabel('Gain over od (points)', fontsize=19, color=SEC, labelpad=12)
    for s in ('top', 'right', 'left', 'bottom'):
        ax.spines[s].set_visible(False)
    ax.grid(axis='y', color='#E4E0D8', lw=1)
    ax.set_axisbelow(True)
    fig.patch.set_alpha(0)
    fig.savefig(OUT / 'chart_precision.png', dpi=200, transparent=True)
    plt.close(fig)


if __name__ == '__main__':
    C = pd.read_csv(REPO / 'results/candidate_features_a3_top300_candidates.csv', low_memory=False)
    P = C[C.in_recorded_pool & C.accepted].copy()
    P['click'] = P.file_name + '|' + P.seed_index.astype(str)
    strips(P)
    schematics()
    chart_auc(P)
    chart_precision()
    print(sorted(os.listdir(OUT)))
