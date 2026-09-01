"""Does the one-match-many-z shortcut hold at the low z levels? No -- not exactly.

`premise_test.py` derives every z level by filtering one deep candidate pool rather than
re-matching. Its own `verify_shortcut` checks that, but only at seed 0 and z in {1.0, 2.0}:
14 of 280 (ROI, seed, z) cells, and none of the z < 1.0 levels that carry the extended
frontier. This probe checks the missing levels -- z in {0.25, 0.5, 0.75}, two seeds, the
three decision-grade ROIs.

Result: 17 of 18 bit-identical, and one genuine violation. On 301.tiff seed 3 at z = 0.75 the
filtered and re-extracted lists differ by 1 detection out of 16,825. `np.argsort` is unstable
in both `extract_peaks` and `nms_by_distance`, so an exact-score tie pair inside the NMS
radius can be ordered one way in a large pool and the other way in a small one. The two
disputed candidates are 1.41 px apart against a 29.6 px match radius -- the same object -- so
no reported metric moves, but "exact, not an approximation" was false.

Run: PYTHONPATH=. python shortcut_lowz_probe.py
Writes results/premise_shortcut_lowz.csv.
"""
import numpy as np, pandas as pd
import premise_test as pt
from midog_utils import channels as ch, chromatin as cm, dataset as ds
from midog_utils import evaluate as ev, experiment as ex, find_and_suppress as fs
from midog_utils import seed_selection as ss

images, ann = ds.load_annotations(); sel = ex.select_domain_images(images, ann)
rows = []
for fn in ("301.tiff", "246.tiff", "201.tiff"):
    row = sel[sel.file_name == fn].iloc[0]; path = f"images/{fn}"
    rgb = ds.load_roi(path); radius = ev.radius_px(ds.roi_mpp(path))
    gray = ch.to_gray_inverted(rgb); hem = cm.hematoxylin_od(rgb)
    img = ch.to_channel(rgb, pt.CHANNEL); gt = ds.image_annotations(ann, fn)
    border = fs.FSConfig(channel=pt.CHANNEL).patch_size // 2
    for s in (0, 3):
        seed, _ = ss.pick_seed(gt[gt.category_id == ds.MITOTIC], gray,
                               np.random.default_rng([s, int(row["image_id"])]),
                               border, rgb.shape)
        base = ss.tightened_base_size(gray, seed["cx"], seed["cy"])
        cfg = fs.FSConfig(channel=pt.CHANNEL, base_size=base)
        pool, diag, maps = pt.match_pool(img, hem, (seed["cx"], seed["cy"]), base,
                                         radius, cfg, keep_map=True)
        for z in (0.25, 0.5, 0.75):
            cut = diag["med"] + z * diag["mad"]
            v = pt.verify_shortcut(maps[0], maps[1], radius, (seed["cx"], seed["cy"]),
                                   cfg, pool, cut)
            v.update({"file_name": fn, "seed_index": s, "z": z})
            rows.append(v); print(v, flush=True)
        del maps
d = pd.DataFrame(rows)
d.to_csv("results/premise_shortcut_lowz.csv", index=False)
print("\n=== identical:", int(d.identical.sum()), "/", len(d), "===")
