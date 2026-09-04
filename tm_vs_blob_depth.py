"""Head-to-head: how many candidates to find EVERY mitosis, TM score vs blob score.

The comparison that has never been run. Prior TM work reports `recall_at_k` (a budgeted
median) and `median_mitosis_rank`; the product requirement is the **tail** -- the rank of
the *worst* mitosis. This measures it for both rankers on the **identical**
`nucleus_blobs` candidate set, so candidate generation is held constant and only the
scoring function differs. That is the fair isolation: TM's own peak-picking is a separate
question and a multi-hour seeded run.

Arms:
* `blob_score`  -- mean hematoxylin of the segmented component. No click.
* `tm_single`   -- Pearson correlation of the candidate's 51 px patch against the clicked
                   cell's 51 px patch. This is exactly `TM_CCOEFF_NORMED`, evaluated at
                   candidate centroids instead of every pixel.
* `tm_rot24`    -- max over 12 rotations x 2 flips of that template, the rotation-invariant
                   version, since mitotic figures have no canonical orientation.

Channel is hematoxylin, which `fs_simple_channel_probe.csv` measures as TM's best
(median AUC 0.843 vs 0.796 rgb, 0.744 gray-inverted) -- i.e. TM is given its best shot.

Protocol: 5 seeds per ROI, the seed dropped from both the candidate list and the evaluation
ground truth (a click correlates with itself at 1.0 by construction), arms judged on the
**worst** seed. Writes results/tm_vs_blob_depth.csv.
"""

from __future__ import annotations

import time

import cv2
import numpy as np
import pandas as pd

from midog_utils import baselines as bl
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from fp_filter_domain import tail_rematched

S = 51                      # template_match.BASE_SIZE
N_SEEDS = 5
ANGLES = np.arange(0, 360, 30.0)   # 12 angles x 2 flips = 24 augmentations


def patch_matrix(chan: np.ndarray, cx, cy):
    """Stack an S x S patch per candidate; rows for border points are flagged invalid."""
    h, w = chan.shape
    half = S // 2
    n = len(cx)
    P = np.zeros((n, S * S), dtype=np.float32)
    ok = np.zeros(n, dtype=bool)
    xi = np.round(cx).astype(int)
    yi = np.round(cy).astype(int)
    for i in range(n):
        x, y = xi[i], yi[i]
        if x - half < 0 or y - half < 0 or x + half >= w or y + half >= h:
            continue
        P[i] = chan[y - half:y + half + 1, x - half:x + half + 1].ravel()
        ok[i] = True
    return P, ok


def zscore_rows(P: np.ndarray) -> np.ndarray:
    """Mean-centre and L2-normalise each row -- the CCOEFF_NORMED transform."""
    P = P - P.mean(axis=1, keepdims=True)
    nrm = np.linalg.norm(P, axis=1, keepdims=True)
    return P / np.where(nrm > 1e-9, nrm, 1.0)


def templates_from(patch2d: np.ndarray) -> np.ndarray:
    """The clicked patch, its flip, and 12 rotations of each, z-scored. (24, S*S)."""
    out = []
    c = S // 2
    for flip in (False, True):
        base = patch2d[:, ::-1] if flip else patch2d
        for a in ANGLES:
            if a == 0.0:
                r = base
            else:
                m = cv2.getRotationMatrix2D((float(c), float(c)), float(a), 1.0)
                r = cv2.warpAffine(np.ascontiguousarray(base), m, (S, S),
                                   flags=cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
            out.append(r.ravel())
    return zscore_rows(np.asarray(out, dtype=np.float32))


def main():
    df = pd.read_csv(".cache_tail/fp_filter_features.csv")
    _images, ann = ds.load_annotations()
    rows = []

    for fn, sub in df.groupby("file_name"):
        t0 = time.time()
        sub = sub.reset_index(drop=True)
        rgb = ds.load_roi(f"images/{fn}")
        chan = bl.to_hematoxylin(rgb)
        del rgb
        radius = ev.radius_px(ds.roi_mpp(f"images/{fn}"))
        gt_all = ds.image_annotations(ann, fn)

        P, ok = patch_matrix(chan, sub["cx"].to_numpy(), sub["cy"].to_numpy())
        Pz = zscore_rows(P)
        del P
        blob = sub["score"].to_numpy(float)

        pos = np.flatnonzero(sub["y"].to_numpy() == 1)
        rng = np.random.default_rng(0)
        seeds = rng.choice(pos, size=min(N_SEEDS, len(pos)), replace=False)

        for si, s in enumerate(seeds):
            if not ok[s]:
                continue
            ann_id = int(sub.at[s, "matched_ann_id"])
            gt = gt_all[gt_all["ann_id"] != ann_id]
            keep = np.ones(len(sub), bool); keep[s] = False
            cand = sub[keep].reset_index(drop=True)

            half = S // 2
            x, y = int(round(sub.at[s, "cx"])), int(round(sub.at[s, "cy"]))
            tpl = chan[y - half:y + half + 1, x - half:x + half + 1]
            T = templates_from(tpl)                     # (24, S*S)

            corr = Pz[keep] @ T.T                       # (n, 24) Pearson correlations
            valid = ok[keep]
            arms = {
                "blob_score": blob[keep],
                "tm_single": np.where(valid, corr[:, 0], -2.0),
                "tm_rot24": np.where(valid, corr.max(axis=1), -2.0),
            }
            for name, sc in arms.items():
                for r in tail_rematched(cand, sc.astype(float), gt, radius):
                    rows.append(dict(file_name=fn, seed=si, ranker=name, **r))
        print(f"{fn}: {time.time()-t0:.0f}s", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv("results/tm_vs_blob_depth.csv", index=False)
    pd.set_option("display.width", 250)
    for q in (0.95, 1.0):
        w = (out[out["quantile"] == q].groupby(["file_name", "ranker"])["depth"]
             .max().unstack())
        w["TM_best"] = w[["tm_single", "tm_rot24"]].min(axis=1)
        w["x_worse_than_blob"] = (w["TM_best"] / w["blob_score"]).round(2)
        print(f"\n=== WORST-of-5-seeds depth to {q:.0%} of captured ===")
        print(w.to_string())


if __name__ == "__main__":
    main()
