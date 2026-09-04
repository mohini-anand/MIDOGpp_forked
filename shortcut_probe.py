"""Is the one-match-many-z shortcut exact? No -- and the main run never checked where it isn't.

`premise_test.py` derives every z level by filtering one deep candidate pool rather than
re-matching, on the argument that greedy score-ordered NMS means a peak above t can only be
suppressed by a higher-scoring peak, which is also above t. Its own `verify_shortcut` tests
that, but only at seed 0 and z in {1.0, 2.0} -- 14 of 280 (ROI, seed, z) cells, and none of
the z < 1.0 levels the extended frontier rests on.

This probe widens the check: 7 ROIs at seed 0 and the three decision-grade ROIs at a second
seed, at z in {0.5, 0.75, 1.0, 2.0}. z = FLOOR_Z is deliberately excluded -- there the filter
is a no-op (the cut equals the extraction floor, so `n_detections == n_pool`), and such a
cell tests re-extraction determinism rather than the shortcut.

`argsort` is unstable in both `extract_peaks` and `nms_by_distance`, so an exact-score tie
pair inside the NMS radius can be ordered one way in a large pool and the other way in a
small one, flipping which of the two survives. That is observed, not hypothetical. A
disagreement counts as benign only if the two lists are the same length AND every disputed
point has a counterpart within `peak_min_distance` -- see `premise_test.verify_shortcut`.

Run: PYTHONPATH=. python shortcut_probe.py
Writes results/premise_shortcut_probe.csv.
"""
import numpy as np, pandas as pd
import premise_test as pt
from midog_utils import channels as ch, chromatin as cm, dataset as ds
from midog_utils import evaluate as ev, experiment as ex, find_and_suppress as fs
from midog_utils import seed_selection as ss

Z_CHECK = (0.5, 0.75, 1.0, 2.0)
DG = ("301.tiff", "246.tiff", "201.tiff")

def main():
    images, ann = ds.load_annotations()
    sel = ex.select_domain_images(images, ann)
    rows = []
    for _, row in sel.iterrows():
        fn = row["file_name"]; path = f"images/{fn}"
        rgb = ds.load_roi(path); radius = ev.radius_px(ds.roi_mpp(path))
        gray = ch.to_gray_inverted(rgb); hem = cm.hematoxylin_od(rgb)
        img = ch.to_channel(rgb, pt.CHANNEL); gt = ds.image_annotations(ann, fn)
        border = fs.FSConfig(channel=pt.CHANNEL).patch_size // 2
        seeds = (0, 3) if fn in DG else (0,)
        for s in seeds:
            seed, _ = ss.pick_seed(gt[gt["category_id"] == ds.MITOTIC], gray,
                                   np.random.default_rng([s, int(row["image_id"])]),
                                   border, rgb.shape)
            base = ss.tightened_base_size(gray, seed["cx"], seed["cy"])
            cfg = fs.FSConfig(channel=pt.CHANNEL, base_size=base)
            pool, diag, maps = pt.match_pool(img, hem, (seed["cx"], seed["cy"]), base,
                                             radius, cfg, keep_map=True)
            for z in Z_CHECK:
                cut = diag["med"] + z * diag["mad"]
                v = pt.verify_shortcut(maps[0], maps[1], radius, (seed["cx"], seed["cy"]),
                                       cfg, pool, cut)
                v.update({"file_name": fn, "seed_index": s, "z": z,
                          "cut": round(cut, 5), "match_radius_px": round(radius, 3),
                          "peak_min_distance": cfg.peak_min_distance})
                rows.append(v); print(v, flush=True)
            del maps
    d = pd.DataFrame(rows)
    d.to_csv("results/premise_shortcut_probe.csv", index=False)
    n_id, n_ben = int(d["identical"].sum()), int(d["benign"].sum())
    print(f"\n{n_id}/{len(d)} cells bit-identical; {n_ben}/{len(d)} benign "
          f"(same length, every disputed point within peak_min_distance)")
    if n_ben < len(d):
        print(d[~d["benign"]].to_string(index=False))
        raise SystemExit("shortcut diverges beyond a tie-break -- derived arms are unsafe")

if __name__ == "__main__":
    main()
