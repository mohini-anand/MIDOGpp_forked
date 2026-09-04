"""Why does the proposal stage miss six mitotic figures? Corrected attribution.

`tail_sensitivity.py`'s first attempt matched each missed annotation to the **nearest
component centroid** and read off its area. That is wrong, and the object audit exposed it:
the four misses on 246.tiff sit in the **97th percentile of chromatin density** -- they are
among the darkest objects in the ROI, not faint ones -- yet the nearest centroid carried
1-17 px. A large, sprawling merged component has a centroid far from any one of the nuclei
inside it, so nearest-centroid finds a stray speck instead.

This asks the right question: **which connected component contains the annotated pixel**,
and why did that component not survive the area gate. It also re-checks the two deep-tail
annotations on 246.tiff, whose matched detection sits 26.8 and 27.7 px away against a
30.2 px match radius -- i.e. at the very edge of what the greedy matcher will accept, which
is the signature of an unrelated neighbour claiming an annotation whose own object was
never proposed.

Writes results/miss_attribution.csv.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from skimage.measure import label, regionprops_table

from midog_utils import baselines as bl
from midog_utils import chromatin as cm
from midog_utils import dataset as ds
from midog_utils import evaluate as ev

# 246 has the four misses and the two edge-of-radius deep matches; 201 the other two misses.
ROIS = ["246.tiff", "201.tiff"]
MIN_AREA, MAX_AREA = 80, 4000


def main():
    _images, ann = ds.load_annotations()
    rows = []
    for fn in ROIS:
        path = f"images/{fn}"
        rgb = ds.load_roi(path)
        radius = ev.radius_px(ds.roi_mpp(path))
        gt = ds.image_annotations(ann, fn).reset_index(drop=True)
        hem_od = cm.hematoxylin_od(rgb)

        mask = bl.tissue_mask(rgb)
        h_chan = bl.to_hematoxylin(rgb)
        thr = bl._tiled_otsu_threshold(h_chan, mask, 512)
        binary = (h_chan > thr) & mask
        lab = label(binary, connectivity=2)
        props = pd.DataFrame(regionprops_table(lab, properties=("label", "area")))
        area_of = dict(zip(props["label"], props["area"]))

        blobs = bl.nucleus_blobs(rgb)
        det = cm.score_detections(blobs, hem_od).reset_index(drop=True)
        _det_out, gt_out = ev.bucket_detections(det, gt, radius)
        mit = gt_out[gt_out["category_id"] == ds.MITOTIC]

        for _, m in mit.iterrows():
            found = bool(m["found"])
            k = int(m["matched_rank"])
            dist = (float(np.hypot(det.at[k, "cx"] - m["cx"], det.at[k, "cy"] - m["cy"]))
                    if found else np.nan)
            # only the two edge-of-radius matches and the outright misses are of interest
            if found and dist < 0.8 * radius:
                continue
            y, x = int(round(m["cy"])), int(round(m["cx"]))
            lid = int(lab[y, x])
            # the component the annotated pixel falls in; 0 means it is background
            a = int(area_of.get(lid, 0)) if lid else 0
            # largest component touching a small disc around the click, in case the click
            # pixel itself sits in a gap between two touching nuclei
            y0, y1 = max(0, y - 10), min(lab.shape[0], y + 11)
            x0, x1 = max(0, x - 10), min(lab.shape[1], x + 11)
            near = np.unique(lab[y0:y1, x0:x1])
            near = [int(v) for v in near if v]
            a_near = max((area_of.get(v, 0) for v in near), default=0)
            reason = ("background_at_click" if lid == 0 and not near else
                      "merged_component_over_max_area" if a_near > MAX_AREA else
                      "fragment_under_min_area" if 0 < a_near < MIN_AREA else
                      "component_passes_gate" if a_near else "background")
            rows.append(dict(
                file_name=fn, ann_id=int(m["ann_id"]),
                status="edge_of_radius_match" if found else "missed",
                matched_rank=k if found else -1,
                ann_to_det_px=round(dist, 1) if found else np.nan,
                match_radius=round(radius, 1),
                comp_area_at_click=a, comp_area_near=a_near,
                od51_pct=round(float((det["od"] < cm.chromatin_density(
                    hem_od, float(m["cx"]), float(m["cy"]))).mean()), 4),
                reason=reason))
        print(f"{fn} done", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv("results/miss_attribution.csv", index=False)
    pd.set_option("display.width", 260)
    print(out.to_string(index=False))
    print("\nreason counts:")
    print(out["reason"].value_counts().to_string())


if __name__ == "__main__":
    main()
