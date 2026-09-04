"""How much of `nucleus_blobs`' output is duplicate at the evaluation's own resolution?

`baselines.nucleus_blobs` emits raw connected components with no suppression, while every
other arm in `premise_test.py` is suppressed at the match radius. Greedy matching is
one-to-one, so a second component inside one match radius of a first can never claim a
second ground-truth object -- it only consumes reading budget.

Two different quantities get confused here, so both are recorded:

* ``adjacency`` -- the fraction of blobs with *any* other blob inside one match radius.
* ``nms_removed`` -- the fraction actually dropped by `nms.nms_by_distance`, which is
  smaller, because a blob is removed only when the neighbour outranks it.

An earlier write-up quoted the second while describing the first.

Run: PYTHONPATH=. python blob_adjacency_probe.py
Writes results/blob_adjacency.csv.
"""
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from midog_utils import baselines as bl, dataset as ds, evaluate as ev, experiment as ex
from midog_utils.nms import nms_by_distance

def main():
    images, ann = ds.load_annotations()
    sel = ex.select_domain_images(images, ann)
    rows = []
    for _, row in sel.iterrows():
        fn = row["file_name"]; path = f"images/{fn}"
        rgb = ds.load_roi(path); radius = ev.radius_px(ds.roi_mpp(path))
        b = bl.nucleus_blobs(rgb); xy = b[["cx", "cy"]].to_numpy()
        tree = cKDTree(xy)
        adj = float(np.mean([len(tree.query_ball_point(p, radius)) > 1 for p in xy]))
        kept = nms_by_distance(xy, b["score"].to_numpy(), radius)
        rows.append({"file_name": fn, "n_blobs": len(b),
                     "match_radius_px": round(radius, 3),
                     "adjacency_frac": round(adj, 4),
                     "n_after_nms": len(kept),
                     "nms_removed_frac": round(1 - len(kept) / len(b), 4)})
        print(rows[-1], flush=True)
    d = pd.DataFrame(rows)
    d.to_csv("results/blob_adjacency.csv", index=False)
    print(f"\nadjacency {d.adjacency_frac.min():.3f}-{d.adjacency_frac.max():.3f}; "
          f"nms removal {d.nms_removed_frac.min():.3f}-{d.nms_removed_frac.max():.3f}")

if __name__ == "__main__":
    main()
