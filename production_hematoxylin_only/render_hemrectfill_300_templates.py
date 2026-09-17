"""
    Render the two hem_rect_fill templates for each 300.tiff click, beside the Otsu steps that
    produced the rectangle. One figure per click, saved next to this script. The templates are
    rebuilt from the pixels and checked against the recorded run: every SHA-1 must match the
    tpl_sha1 of `precision_at_k_49roi_3seed_chromatin_hemrectfill_per_run.csv`, so the panels
    show the arrays `hem_rect_fill_vs_default_chromatin_od_49roi_3seed.ipynb` actually searched with.
"""

from __future__ import annotations

import hashlib
import os
import sys

import cv2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from skimage.measure import label, regionprops

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from midog_utils import channels as ch
from midog_utils import dataset as ds
from midog_utils import seed_selection as ss
from midog_utils import template_match as tm

STEM = f'{ROOT}/results/precision_at_k_49roi_3seed_chromatin_hemrectfill'
FILE_NAME = '300.tiff'
ROI_PATH = f'{ROOT}/images/extra_valid/{FILE_NAME}'


def main():
    """
        Render one figure per 300.tiff click, after checking each template against the recorded run.

        Returns None: writes hemrectfill_300_<ann_id>.png beside this script.
    """
    click_table = pd.read_csv(f'{STEM}_clicks.csv', float_precision='round_trip')
    per_run = pd.read_csv(f'{STEM}_per_run.csv', float_precision='round_trip')
    rgb = ds.load_roi(ROI_PATH)
    hem = ch.to_channel(rgb, 'hematoxylin_od')

    for _, run in per_run[(per_run['file_name'] == FILE_NAME) & (per_run['condition'] == 'default_51')].sort_values('seed_index').iterrows():
        s, ann_id, cx, cy = int(run['seed_index']), int(run['seed_ann_id']), float(run['click_cx']), float(run['click_cy'])
        ix, iy = int(round(cx)), int(round(cy))

        # the cut the notebook's cut_templates makes
        patch = tm.read_padded_patch(hem, cx, cy, 73)
        assert patch is not None, f'{FILE_NAME} s{s}: the 73 px patch is not readable'
        templates, _ = tm.build_augmentations(patch, 51, (1.0,), 1, (False,))
        T = templates[0]
        W = tm.read_padded_patch(hem, cx, cy, 51)
        assert W is not None and np.array_equal(W, T), f'{FILE_NAME} s{s}: the Otsu window and the template are not the same pixels'
        bbox = ss.tighten_box_otsu(W)
        assert bbox is not None, f'{FILE_NAME} s{s}: tighten_box_otsu refused the window'
        y0, y1, x0, x1 = bbox
        outside = np.ones((51, 51), bool)
        outside[y0:y1, x0:x1] = False
        fill = np.float32(T[outside].astype(np.float64).mean())
        F = T.copy()
        F[outside] = fill

        # the arrays shown are the arrays the notebook ran: sha1 and rectangle must match its CSVs
        fill_run = per_run[(per_run['file_name'] == FILE_NAME) & (per_run['seed_index'] == s) & (per_run['condition'] == 'hem_rect_fill')].iloc[0]
        clk = click_table[(click_table['file_name'] == FILE_NAME) & (click_table['seed_index'] == s)].iloc[0]
        assert hashlib.sha1(T.tobytes()).hexdigest() == run['tpl_sha1'], f'{FILE_NAME} s{s}: default_51 template differs from the recorded one'
        assert hashlib.sha1(F.tobytes()).hexdigest() == fill_run['tpl_sha1'], f'{FILE_NAME} s{s}: hem_rect_fill template differs from the recorded one'
        assert (y0, y1, x0, x1) == (clk['rect_y0'], clk['rect_y1'], clk['rect_x0'], clk['rect_x1']) and int(clk['seed_ann_id']) == ann_id, f'{FILE_NAME} s{s}: rectangle differs from the recorded one'

        # tighten_box_otsu's own steps, kept for display, then checked against its result
        u8 = cv2.normalize(W.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
        cut, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        labels = label(binary, connectivity=2)
        comp = labels[25, 25]
        region = next(p for p in regionprops(labels) if p.label == comp)
        assert region.bbox == (y0, x0, y1, x1), f'{FILE_NAME} s{s}: displayed component is not the one tighten_box_otsu accepted'

        rect = dict(xy=(x0 - 0.5, y0 - 0.5), width=x1 - x0, height=y1 - y0, fill=False, edgecolor='red', linewidth=2)
        vmin, vmax = float(T.min()), float(T.max())

        fig, axes = plt.subplots(1, 5, figsize=(21, 5.2))
        axes[0].imshow(rgb[iy - 25: iy + 26, ix - 25: ix + 26])
        axes[0].add_patch(Rectangle(**rect))
        axes[0].set_title('1. RGB, 51 px window on the click', fontsize=10)

        axes[1].imshow(binary, cmap='gray')
        axes[1].set_title(f'2. hematoxylin_od Otsu binary (cut={cut:.0f}/255)', fontsize=10)

        comp_rgb = np.zeros((51, 51, 3))
        comp_rgb[labels > 0] = (0.55, 0.55, 0.55)  # other foreground components: grey
        comp_rgb[labels == comp] = (0.1, 0.35, 0.9)  # the accepted component: blue
        axes[2].imshow(comp_rgb)
        axes[2].contour(labels == comp, levels=[0.5], colors='lime', linewidths=2)
        axes[2].add_patch(Rectangle(**rect))
        axes[2].set_title(f'3. component under the click (blue, lime outline; other blobs grey)\narea={int(region.area)} px, solidity={region.solidity:.2f}; red = its bounding rectangle', fontsize=10)

        axes[3].imshow(T, cmap='gray', vmin=vmin, vmax=vmax)
        axes[3].add_patch(Rectangle(**rect))
        axes[3].set_title(f'4. default_51 template (51x51)\nsha1 {run["tpl_sha1"][:10]}', fontsize=10)

        axes[4].imshow(F, cmap='gray', vmin=vmin, vmax=vmax)
        axes[4].add_patch(Rectangle(**rect))
        axes[4].set_title(f'5. hem_rect_fill template (51x51)\noutside rectangle = {float(fill):.4f} ({int(outside.sum())} px, {outside.mean():.0%}); sha1 {fill_run["tpl_sha1"][:10]}', fontsize=10)

        for ax in axes:
            ax.plot(25, 25, marker='+', color='red' if ax is not axes[2] else 'white', ms=12, mew=2)
            ax.set_xticks([])
            ax.set_yticks([])
        fig.suptitle(f'{FILE_NAME} seed {s} (annotation {ann_id}, click {ix},{iy}): rectangle rows {y0}-{y1 - 1}, cols {x0}-{x1 - 1} of the 51 px window = {y1 - y0}x{x1 - x0} px.  '
                     f'Panels 2-5 are hematoxylin_od, bright = more hematoxylin; 4 and 5 share one grey scale.  '
                     f'TP@10/20/30: default_51 {tuple(int(run[f"tp_at_{k}"]) for k in (10, 20, 30))}, hem_rect_fill {tuple(int(fill_run[f"tp_at_{k}"]) for k in (10, 20, 30))}', fontsize=10.5)
        fig.tight_layout()
        out = f'{HERE}/hemrectfill_300_{ann_id}.png'
        fig.savefig(out, dpi=130, bbox_inches='tight')
        plt.close(fig)
        print(f'seed {s} ann {ann_id}: rect {(y0, y1, x0, x1)}, both templates match the recorded run -> {out}')


if __name__ == '__main__':
    main()
