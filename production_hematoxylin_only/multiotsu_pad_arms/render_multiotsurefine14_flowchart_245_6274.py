"""
    Render refine_template_box step by step for one click, 245.tiff / 6274, as a two-row flowchart
    from the 51 px input template to the refined template. Saved next to this script. The method is
    exec'd from `multiotsu_template_refinement_14roi.ipynb` (cells 1, 3, 4), not copied, and the
    result is asserted against that notebook's recorded box for this click.
"""

import json
import os

import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib.colors import to_rgb
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle
from skimage.measure import label

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)                                              # the notebook cells use paths relative to this folder

NB = 'multiotsu_template_refinement_14roi.ipynb'
FN, ANN_ID = '245.tiff', 6274
OUT = 'multiotsurefine14_flowchart_245_6274.png'

# Pull the method in verbatim: imports/constants (cell 1), refine_template_box (cell 3), context_crop (cell 4).
# The asserts below catch it if those cells move or the method changes.
g = {}
cells = json.load(open(NB))['cells']
for i in (1, 3, 4):
    exec(''.join(cells[i]['source']), g)
ds, ch, tm, plt = g['ds'], g['ch'], g['tm'], g['plt']
WINDOW, CONTEXT_WINDOW, CLASSES = g['WINDOW'], g['CONTEXT_WINDOW'], g['CLASSES']

images, annotations = ds.load_annotations(g['DB_PATH'])
row = ds.image_annotations(annotations, FN)
row = row[row['ann_id'] == ANN_ID]
assert len(row) == 1 and int(row.iloc[0]['category_id']) == ds.MITOTIC
cx, cy = float(row.iloc[0]['cx']), float(row.iloc[0]['cy'])
domain = images.set_index('file_name')['tumor_type'][FN]
rgb = ds.load_roi(f"{g['IMAGES_DIR']}/{FN}")
mpp = ds.roi_mpp(f"{g['IMAGES_DIR']}/{FN}")
hem = ch.to_channel(rgb, 'hematoxylin_od')
res = g['refine_template_box'](hem, cx, cy, mpp)

# Fidelity guards: the figure must show what the method actually does on this click.
assert res['accepted'] and res['mode'] == 'fallback' and res['selected'] == [3, 6], res['reason']
assert res['bbox'] == (10, 44, 12, 42) and res['base_size'] == 39 and res['center_xy'] == (6357.5, 536.5)

# ---------------------------------------------------------------- intermediates
patch = res['hem_patch'].astype(np.float64)
c = WINDOW // 2
thr = res['thresholds']
top_raw = patch > thr[-1]                                   # the method's own expression, before dilation
labels = res['labels']
dilated = labels > 0
added = dilated & ~top_raw
n_raw, n_dil = label(top_raw, connectivity=1).max(), labels.max()
pieces = res['pieces']
specks = [p for p in pieces if p['speck']]
real = [p for p in pieces if not p['speck']]
sel = res['selection_mask']
y0, y1, x0, x1 = res['bbox']
bs, pad = res['base_size'], res['pad_px']
long_side = max(y1 - y0, x1 - x0)

ix, iy = int(round(cx)), int(round(cy))                     # window origin = click - half
tcx, tcy = int(round(res['center_xy'][0])), int(round(res['center_xy'][1]))   # read_padded_patch's rounding
wx, wy = tcx - (ix - c), tcy - (iy - c)                     # template centre in window coordinates
bh = (bs - 1) // 2
rgb_win = tm.read_padded_patch(rgb, cx, cy, WINDOW)
template = tm.read_padded_patch(rgb, res['center_xy'][0], res['center_xy'][1], bs)
assert np.array_equal(template, rgb_win[wy - bh:wy + bh + 1, wx - bh:wx + bh + 1])
ctx = g['context_crop'](rgb, cx, cy)

# ---------------------------------------------------------------- style
INK, INK2, ARROW = '#0b0b0b', '#52514e', '#6b6a66'
CLICK, RADIUS, SEL, NOTSEL, BOX = 'red', 'yellow', 'lime', 'red', 'cyan'
ADDED = '#eb6834'
SPECK = '#6b6a66'
PIECE_COLOURS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#4a3aa7']
TITLE_FS, NOTE_FS = 9.5, 8

W, H = 21.5, 9.2
S = 2.6                                                     # image panel side, inches
fig = plt.figure(figsize=(W, H), facecolor='white')
IN = fig.dpi_scale_trans


def panel(x, y, w=S, h=S):
    """
        Add an image axes at a position given in inches.
        x, y (float): lower-left corner, inches from the figure's lower-left.
        w, h (float): width and height, inches.
        Returns matplotlib.axes.Axes: the axes, ticks off.
    """
    ax = fig.add_axes([x / W, y / H, w / W, h / H])
    ax.set_xticks([])
    ax.set_yticks([])
    return ax


def title(x, y, w, text):
    """
        Write a step label centred above a node.
        x, y (float): node's lower-left corner, inches.
        w (float): node width, inches; the label sits over its centre.
        text (str): first line is the step name (bold), the rest a short note.
        Returns None: draws on the figure.
    """
    head, *rest = text.split('\n')
    lift = 0.40 if y == ROW1 else 0.58
    fig.text(x + w / 2, y + S + lift, head, transform=IN, ha='center', va='bottom', fontsize=TITLE_FS, color=INK, fontweight='bold')
    if rest:
        fig.text(x + w / 2, y + S + lift - 0.05, '\n'.join(rest), transform=IN, ha='center', va='top', fontsize=NOTE_FS, color=INK2, linespacing=1.35)


def arrow(p0, p1, colour=ARROW, style='-|>', ls='-'):
    """
        Draw a straight arrow between two points given in inches.
        p0, p1 (tuple): start and end (x, y), inches.
        colour (str): arrow colour.
        style (str): matplotlib arrowstyle.
        ls (str): line style.
        Returns None: draws on the figure.
    """
    fig.patches.append(FancyArrowPatch(p0, p1, transform=IN, arrowstyle=style, mutation_scale=16, color=colour, lw=1.6, linestyle=ls, shrinkA=0, shrinkB=0))


def mark_click(ax, at=(c, c), colour=CLICK):
    """
        Mark the click with a plus.
        ax (matplotlib.axes.Axes): target axes.
        at (tuple): click (x, y) in the axes' pixel coordinates.
        colour (str): marker colour.
        Returns None: draws into ax.
    """
    ax.plot(*at, marker='+', color=colour, ms=12, mew=2)


def pieces_rgb(colour_of):
    """
        Paint the dilated pieces on black.
        colour_of (callable): piece dict -> colour string, or None to leave it black.
        Returns np.ndarray: WINDOW x WINDOW x 3 float image.
    """
    img = np.zeros((WINDOW, WINDOW, 3))
    for p in pieces:
        col = colour_of(p)
        if col is not None:
            img[labels == p['label']] = to_rgb(col)
    return img


# ---------------------------------------------------------------- layout
X0, G1 = 0.55, 0.95
ROW1, ROW2 = 5.4, 1.05
xs1 = [X0 + i * (S + G1) for i in range(6)]
mid1, mid2 = ROW1 + S / 2, ROW2 + S / 2

# 1. click + 51 px window
ax = panel(xs1[0], ROW1)
ax.imshow(ctx)
hc = CONTEXT_WINDOW // 2
ax.add_patch(Rectangle((hc - c - 0.5, hc - c - 0.5), WINDOW, WINDOW, fill=False, edgecolor='white', linestyle='--', linewidth=1.6))
mark_click(ax, (hc, hc))
title(xs1[0], ROW1, S, f'1. click on the figure\nred + = click, dashed = {WINDOW} px window')

# 2. input template
ax = panel(xs1[1], ROW1)
ax.imshow(rgb_win)
mark_click(ax)
title(xs1[1], ROW1, S, f'2. input template\n{WINDOW} × {WINDOW} px window')

# 3. hematoxylin OD
ax = panel(xs1[2], ROW1)
ax.imshow(patch, cmap='gray')
mark_click(ax)
title(xs1[2], ROW1, S, '3. hematoxylin OD channel\nbright = more hematoxylin')

# 4. multi-Otsu histogram
ax = fig.add_axes([(xs1[3] + 0.42) / W, (ROW1 + 0.42) / H, (S - 0.47) / W, (S - 0.47) / H])
bins = np.linspace(0, patch.max(), 41)
counts, edges = np.histogram(patch, bins=bins)
mids = (edges[:-1] + edges[1:]) / 2
cls = np.digitize(mids, thr)
face = np.array(['#2b2b2b', '#9a9a96', 'white'])[cls]
ax.bar(edges[:-1], counts, width=np.diff(edges), align='edge', color=face, edgecolor=INK2, linewidth=0.5)
for t in thr:
    ax.axvline(t, color=INK, linestyle='--', linewidth=1.2)
    ax.text(t, counts.max() * 1.04, f'{t:.3f}', ha='center', va='bottom', fontsize=7.5, color=INK, bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none'))
ax.set_ylim(0, counts.max() * 1.15)
ax.set_xlabel('hematoxylin OD', fontsize=7.5, color=INK2, labelpad=2)
ax.set_ylabel('pixels', fontsize=7.5, color=INK2, labelpad=2)
ax.tick_params(labelsize=6.5, colors=INK2, length=2)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
for s in ('left', 'bottom'):
    ax.spines[s].set_color('#b5b4ae')
title(xs1[3], ROW1, S, f'4. {CLASSES}-class multi-Otsu\n2 cuts split the window into 3 classes')

# 5. keep the top class
ax = panel(xs1[4], ROW1)
ax.imshow(res['class_map'], cmap='gray', vmin=0, vmax=CLASSES - 1)
mark_click(ax)
title(xs1[4], ROW1, S, f'5. keep the top class (white)\ndarkest chromatin: {int(top_raw.sum())} px, {n_raw} pieces')

# 6. 1 px dilation
ax = panel(xs1[5], ROW1)
img = np.zeros((WINDOW, WINDOW, 3))
img[top_raw] = 1.0
img[added] = to_rgb(ADDED)
ax.imshow(img)
mark_click(ax)
title(xs1[5], ROW1, S, f'6. dilate 1 px (3×3 cross)\norange = added, {n_raw} → {n_dil} pieces')

for a, b in zip(xs1[:-1], xs1[1:]):
    arrow((a + S + 0.1, mid1), (b - 0.1, mid1))

# ---- row 2
G2 = 1.0
DIA_W, DIA_H = 2.1, 1.75
x7 = X0
xdia = x7 + S + G2
x9 = xdia + DIA_W + G2
x10 = x9 + S + G2
x11 = x10 + S + G2
SO = S * bs / WINDOW                                        # output at the input template's pixel scale
x12 = x11 + S + G2

# wrap arrow: from under panel 6 to above panel 7's label
yw = ROW1 - 0.32
fig.lines.append(plt.Line2D([xs1[5] + S / 2, xs1[5] + S / 2, x7 + S / 2], [ROW1 - 0.08, yw, yw], transform=IN, color=ARROW, lw=1.6))
arrow((x7 + S / 2, yw), (x7 + S / 2, ROW2 + S + 0.86))

# 7. label pieces, drop specks
ax = panel(x7, ROW2)
colour_of = {p['label']: PIECE_COLOURS[i] for i, p in enumerate(real)}
ax.imshow(pieces_rgb(lambda p: SPECK if p['speck'] else colour_of[p['label']]))
for p in specks:
    ys, xs = np.where(labels == p['label'])
    ax.add_patch(Circle((xs.mean(), ys.mean()), 3.2, fill=False, edgecolor='white', linewidth=1))
mark_click(ax, colour='white')
ax.set_xlim(-0.5, WINDOW - 0.5)
ax.set_ylim(WINDOW - 0.5, -0.5)
title(x7, ROW2, S, f'7. label pieces (4-connected)\n{len(real)} pieces + {len(specks)} specks < {g["MIN_PIECE_UM2"]:.0f} µm² (grey, circled) dropped')

# 8. base-case check (decision diamond)
dcx, dcy = xdia + DIA_W / 2, mid2
fig.patches.append(Polygon([(dcx, dcy + DIA_H / 2), (dcx + DIA_W / 2, dcy), (dcx, dcy - DIA_H / 2), (dcx - DIA_W / 2, dcy)], closed=True, transform=IN, facecolor='#f4f3ef', edgecolor=INK2, linewidth=1.4))
fig.text(dcx, dcy, f'click on a piece\nwith bbox ≥ {g["MIN_BASE_UM"]:.0f} µm?', transform=IN, ha='center', va='center', fontsize=8.5, color=INK, linespacing=1.3)
title(xdia, ROW2, DIA_W, '8. base case?')
arrow((x7 + S + 0.1, mid2), (dcx - DIA_W / 2 - 0.05, mid2))
arrow((dcx + DIA_W / 2 + 0.05, mid2), (x9 - 0.1, mid2))
fig.text((dcx + DIA_W / 2 + x9) / 2, mid2 + 0.1, 'no', transform=IN, ha='center', va='bottom', fontsize=8.5, color=INK, fontweight='bold')
fig.text((dcx + DIA_W / 2 + x9) / 2, mid2 - 0.1, 'click is on\nbackground', transform=IN, ha='center', va='top', fontsize=7.5, color=INK2)
arrow((dcx, dcy - DIA_H / 2 - 0.05), (dcx, 0.62), colour='#b5b4ae', ls='--')
fig.text(dcx, 0.52, 'yes: keep the clicked piece alone\n(not taken here)', transform=IN, ha='center', va='top', fontsize=7.5, color='#8a8984', linespacing=1.3)

# 9. fallback: piece centres within the radius
ax = panel(x9, ROW2)
ax.imshow(pieces_rgb(lambda p: '#3a3a38' if p['speck'] else ('#2f7d2f' if p['selected'] else '#7a2d2d')))
for p in real:
    ax.contour(labels == p['label'], levels=[0.5], colors=SEL if p['selected'] else NOTSEL, linewidths=1.5 if p['selected'] else 1)
ax.add_patch(Circle((c, c), res['radius_px'], fill=False, edgecolor=RADIUS, linestyle='--', linewidth=1.5))
for p in real:
    py_, px_ = p['centre_yx']
    col = SEL if p['selected'] else NOTSEL
    ax.plot([c, px_], [c, py_], color=col, lw=0.8, alpha=0.9)
    ax.plot(px_, py_, 'o', color=col, ms=4, mec='black', mew=0.5)
    ax.text(float(np.clip(px_, 3, WINDOW - 4)), float(np.clip(py_ + (2.6 if py_ > c else -2.6), 2, WINDOW - 3)), f"{p['centre_um']:.2f}", color='white', fontsize=6.5, ha='center', va='center', bbox=dict(boxstyle='round,pad=0.12', fc='black', ec='none', alpha=0.6))
mark_click(ax, colour='white')
ax.set_xlim(-0.5, WINDOW - 0.5)
ax.set_ylim(WINDOW - 0.5, -0.5)
n_sel = len(res['selected'])
title(x9, ROW2, S, f'9. fallback: piece centre ≤ {g["CLICK_RADIUS_UM"]} µm\nyellow = radius, numbers = centre distance (µm)\nlime = selected ({n_sel}), red = not selected')

# 10. selection, area gates, bbox
ax = panel(x10, ROW2)
ax.imshow(patch, cmap='gray', vmin=0, vmax=patch.max() * 2.5)
over = np.zeros((WINDOW, WINDOW, 4))
over[sel] = (*to_rgb(SEL), 0.85)
ax.imshow(over)
ax.add_patch(Rectangle((x0 - 0.5, y0 - 0.5), x1 - x0, y1 - y0, fill=False, edgecolor='white', linestyle='--', linewidth=1.5))
ax.text(x0 + (x1 - x0) / 2, y1 + 0.2, f'{x1 - x0} px', color='white', fontsize=7, ha='center', va='top')
ax.text(x1 + 0.2, y0 + (y1 - y0) / 2, f'{y1 - y0} px', color='white', fontsize=7, ha='left', va='center', rotation=-90)
mark_click(ax)
area = sel.sum() * mpp ** 2
frac = sel.sum() / patch.size
title(x10, ROW2, S, f'10. selection + bounding box\n{area:.1f} µm² ≥ {g["MIN_AREA_UM2"]:.0f} µm² ✓   {frac:.0%} of window ≤ {g["MAX_AREA_FRAC"]:.0%} ✓\nbbox {y1 - y0} × {x1 - x0} px (dashed)')

# 11. square box centred on the bbox
ax = panel(x11, ROW2)
ax.imshow(rgb_win)
ax.contour(sel, levels=[0.5], colors=SEL, linewidths=1)
ax.add_patch(Rectangle((x0 - 0.5, y0 - 0.5), x1 - x0, y1 - y0, fill=False, edgecolor='white', linestyle='--', linewidth=1.2))
ax.add_patch(Rectangle((wx - bh - 0.5, wy - bh - 0.5), bs, bs, fill=False, edgecolor=BOX, linewidth=2))
ax.plot(wx, wy, 'o', color=SEL, ms=5, mec='black', mew=0.5)
mark_click(ax)
title(x11, ROW2, S, f'11. square box on the bbox centre\nside = odd({long_side}) + 2 × {pad} px ({g["PAD_UM"]} µm margin) = {bs} px\ncyan = box, lime dot = new centre')

arrow((x9 + S + 0.1, mid2), (x10 - 0.1, mid2))
arrow((x10 + S + 0.1, mid2), (x11 - 0.1, mid2))

# 12. output
ax = panel(x12, ROW2 + (S - SO) / 2, SO, SO)
ax.imshow(template)
for s in ax.spines.values():
    s.set_edgecolor(BOX)
    s.set_linewidth(2.5)
title(x12, ROW2, SO, f'12. refined template\n{bs} × {bs} px (same scale as 2)')
arrow((x11 + S + 0.1, mid2), (x12 - 0.12, mid2))

offset = np.hypot(res['center_xy'][0] - cx, res['center_xy'][1] - cy)
fig.text(W / 2, H - 0.22, f'Multi-Otsu template refinement, step by step:  {FN} / annotation {ANN_ID}  ({domain}, {mpp:.3f} µm/px)', transform=IN, ha='center', va='top', fontsize=13, color=INK, fontweight='bold')

fig.savefig(OUT, dpi=130, facecolor='white')
print(f'-> {OUT}; box {bs} px, centre offset {offset:.1f} px from the click')
