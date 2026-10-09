"""Gate-3 dry check on a few ROIs: b0 bit-exact vs floors.csv (as float32), W128 vs recorded, refinement boxes reproduced from hem0."""
import sys, json, hashlib, inspect, time
from common import *
import numpy as np, cv2
from skimage.measure import label, regionprops
from midog_utils import template_match as tm
from midog_utils import seed_selection as ss

nb = json.load(open(REPO / 'production_hematoxylin_only/multiotsu_pad_arms/tm_border_floor_49roi_6seed.ipynb'))
src = ''.join(nb['cells'][4]['source'])
src = src[:src.index('\n\nassert hashlib')]
WINDOW = tm.BASE_SIZE; CLASSES = 3
DILATE_KERNEL = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3)); DILATE_ITERATIONS = 1
MIN_PIECE_UM2 = 1.0; MIN_BASE_UM = 5.0; CLICK_RADIUS_UM = 3.75; MIN_AREA_UM2 = 3.0; MAX_AREA_FRAC = 0.85; PAD_PX = 2
ns = dict(globals())
exec(src, ns)
refine_template_box = ns['refine_template_box']
print('refine source sha1 (exec-defined; inspect needs a file):', hashlib.sha1((src + '\n').encode()).hexdigest()[:12], 'vs 21c43a20db00')

for fn in sys.argv[1:]:
    rgb, hem, mpp = load(fn)
    b0 = roi_floors(hem)['b_otsu']
    side = int(round(128 / mpp)) | 1
    n_ok = 0; n_w = 0; maxdiff_w = 0
    for r in CLK[CLK.file_name == fn].itertuples():
        res = refine_template_box(hem, float(r.click_cx), float(r.click_cy), mpp)
        if res['accepted']:
            bs, tcx, tcy = int(res['base_size']), int(round(res['center_xy'][0])), int(round(res['center_xy'][1]))
            n_ok += (bs, tcx, tcy) == (int(r.bs), int(r.tcx), int(r.tcy))
        bw = window_floor(hem, r.click_cx, r.click_cy, side)[0]
        rec = W128_REC[(fn, r.seed_index)]
        n_w += np.float32(bw) == np.float32(rec)
        maxdiff_w = max(maxdiff_w, abs(bw - rec) / rec)
    print(f'{fn}: b0 float32-equal to floors.csv: {np.float32(b0) == np.float32(B_REC[fn])} (rel diff {abs(b0 - B_REC[fn]) / B_REC[fn]:.1e});  W128 float32-equal {n_w}/6 (max rel diff {maxdiff_w:.1e});  boxes reproduced {n_ok}/6;  hem0==0 share {np.mean(hem == 0):.3f}')
