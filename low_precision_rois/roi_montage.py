"""Crop montages for the low-precision ROIs: GT mitoses vs A3 top-10 (tm_score, chromatin_od), click 0 and 1."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = Path('/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked')
sys.path.insert(0, str(REPO))
from midog_utils_full import dataset as ds

OUT = Path(sys.argv[1])
FILES = sys.argv[2].split(',')
C = 101
j = json.load(open(REPO / 'databases/MIDOG++.json'))
imgs = {i['id']: i for i in j['images']}
ann = pd.DataFrame(j['annotations'])
ann['file_name'] = ann.image_id.map(lambda i: imgs[i]['file_name'])
ann['cx'] = ann.bbox.map(lambda b: (b[0] + b[2]) / 2)
ann['cy'] = ann.bbox.map(lambda b: (b[1] + b[3]) / 2)
top = pd.read_csv(REPO / 'results/tm_border_floor_49roi_6seed_top30.csv')
top = top[(top.arm == 'A3') & (top['rank'] < 10)]
TAG = {'human_correct_label': 'MIT', 'human_rejected_label': 'LA', 'non_human_findings': '-'}


def crop(rgb, cx, cy):
    h = C // 2
    p = np.pad(rgb, ((h, h), (h, h), (0, 0)), mode='constant', constant_values=0)
    ix, iy = int(round(cx)) + h, int(round(cy)) + h
    return p[iy - h:iy + h + 1, ix - h:ix + h + 1]


for fn in FILES:
    p = REPO / 'images/extra_valid' / fn
    if not p.exists():
        p = REPO / 'images/extra_valid/testing_set' / fn
    rgb = ds.load_roi(p)
    a = ann[ann.file_name == fn]
    mit = a[a.category_id == 1].sample(n=min(10, (a.category_id == 1).sum()), random_state=0)
    la = a[a.category_id == 2].sample(n=min(10, (a.category_id == 2).sum()), random_state=0)
    rows = [('GT mitoses (random 10)', [(crop(rgb, x, y), f'v{sum(l == 1 for l in L)}/{len(L)}') for x, y, L in zip(mit.cx, mit.cy, mit['labels'])]),
            ('GT look-alikes (random 10)', [(crop(rgb, x, y), f'v{sum(l == 1 for l in L)}/{len(L)}') for x, y, L in zip(la.cx, la.cy, la['labels'])])]
    for s in (0, 1):
        for rk in ('tm_score', 'chromatin_od'):
            tt = top[(top.file_name == fn) & (top.seed_index == s) & (top.rank_key == rk)].sort_values('rank')
            rows.append((f'A3 {rk} click {s}', [(crop(rgb, x, y), f"{TAG[b]} e{int(e)} w{w:.2f}") for x, y, b, e, w in zip(tt.cx, tt.cy, tt.bucket, tt.edge_px, tt.white_frac_61)]))
    fig, axs = plt.subplots(len(rows), 10, figsize=(20, 2.2 * len(rows)))
    for r, (title, cells) in enumerate(rows):
        for c in range(10):
            ax = axs[r, c]
            ax.axis('off')
            if c < len(cells):
                ax.imshow(cells[c][0])
                ax.set_title(cells[c][1], fontsize=8)
        axs[r, 0].text(-0.1, 0.5, title, transform=axs[r, 0].transAxes, rotation=90, va='center', ha='right', fontsize=9)
    fig.suptitle(f'{fn} ({a.file_name.iat[0]})  crops {C}px', fontsize=12)
    plt.tight_layout()
    fig.savefig(OUT / f'montage_{fn.replace(".tiff", "")}.jpg', dpi=70)
    plt.close(fig)
    print(fn, flush=True)
