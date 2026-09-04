"""Is the watershed arm's 100% capture real detection, or geometric carpeting?

`results/watershed_tail.csv` reports full-list recall of 1.0000 on all seven ROIs for the
watershed arm -- 653 of 653 annotated mitotic figures. The 2026-09-01 audit (M6) warned
that exactly this kind of number can be an artifact: a candidate list long enough to tile
the image will "find" every annotation regardless of whether it detected anything, because
every point in the ROI is within the match radius of *some* candidate.

The watershed pools are 12,611-60,609 candidates and the match radius is ~30 px, so the
discs sum to several times the ROI area. Whether they actually cover it depends on how much
they overlap, which is not something to reason about -- it is something to measure.

Two measurements per (ROI, arm):

* `coverage_frac` -- the fraction of the ROI within the match radius of some candidate,
  on a strided probe grid. If this is near 1.0 the recall number means nothing.
* the distribution of **annotation-to-matched-detection distance** for mitotic figures.
  Real detection puts a candidate centroid on the object, a few px out. Carpeting matches
  at the edge of the radius. The pre-split run already showed five annotations matched at
  26.7-28.4 px against a 30.2 px radius, which is the carpeting signature.

Writes results/blob_coverage.csv.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from skimage.measure import label, regionprops

from midog_utils import baselines as bl
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from watershed_split import MIN_AREA, gate_and_rank, split_components

ROIS = ["301.tiff", "300.tiff", "246.tiff", "245.tiff", "201.tiff", "202.tiff", "405.tiff"]


def main():
    _images, ann = ds.load_annotations()
    rows = []
    for fn in ROIS:
        t0 = time.time()
        path = f"images/{fn}"
        rgb = ds.load_roi(path)
        shape = rgb.shape[:2]
        radius = ev.radius_px(ds.roi_mpp(path))
        gt = ds.image_annotations(ann, fn)

        mask = bl.tissue_mask(rgb)
        h_chan = bl.to_hematoxylin(rgb)
        thr = bl._tiled_otsu_threshold(h_chan, mask, 512)
        binary = (h_chan > thr) & mask
        del thr, mask, rgb

        lab = label(binary, connectivity=2)
        cc = pd.DataFrame([(p.centroid[1], p.centroid[0], p.area, p.intensity_mean)
                           for p in regionprops(lab, intensity_image=h_chan)],
                          columns=["cx", "cy", "area", "score"])
        del lab
        arms = {"connected_components": gate_and_rank(cc),
                "watershed": gate_and_rank(split_components(binary, h_chan))}

        for arm, det in arms.items():
            xy = det[["cx", "cy"]].to_numpy()
            cov = ev.coverage_fraction(xy, shape, radius, stride=16)
            _d, gt_out = ev.bucket_detections(det.reset_index(drop=True),
                                              gt.reset_index(drop=True), radius)
            mit = gt_out[(gt_out["category_id"] == ds.MITOTIC) & gt_out["found"]]
            k = mit["matched_rank"].to_numpy()
            dist = np.hypot(det["cx"].to_numpy()[k] - mit["cx"].to_numpy(),
                            det["cy"].to_numpy()[k] - mit["cy"].to_numpy())
            rows.append(dict(
                file_name=fn, arm=arm, n_pool=len(det), radius=round(radius, 1),
                coverage_frac=round(float(cov), 4),
                n_captured=len(dist), n_gt=int((gt.category_id == ds.MITOTIC).sum()),
                dist_median=round(float(np.median(dist)), 1),
                dist_p90=round(float(np.percentile(dist, 90)), 1),
                dist_max=round(float(dist.max()), 1),
                frac_beyond_80pct_radius=round(float((dist > 0.8 * radius).mean()), 3)))
        print(f"{fn}: {time.time()-t0:.0f}s", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv("results/blob_coverage.csv", index=False)
    pd.set_option("display.width", 250)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
