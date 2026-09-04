"""Does making darkness *click-conditioned* beat using darkness on its own?

The proposal: keep the seeded protocol, but replace `TM_CCOEFF_NORMED` with a similarity
built on the chromatin-darkness signal that the seedless work found useful -- i.e. rank a
candidate by how close its darkness is to the clicked cell's.

The concern this tests: darkness is **monotone** with respect to the label (mitoses are
darker than ordinary nuclei, Bhattacharyya 0.93-3.70 in all seven domains). A two-sided
similarity `-|d - d_seed|` discards that -- it pushes down every candidate *darker* than
whichever mitosis the pathologist happened to click, and the seed is only one draw from the
mitotic darkness distribution. So the click may be actively harmful here, in a way it is not
for a non-monotone representation like an embedding.

Arms, all on the identical cached `nucleus_blobs` candidate set so only the ordering differs:

* `blob_abs`, `od51_abs`  -- absolute darkness, no click. The incumbents.
* `od51_sim`             -- two-sided `-|d - d_seed|`. The literal proposal.
* `od51_onesided`        -- `-relu(d_seed - d)`, ties broken by darkness. Penalises being
                            *paler* than the click but not darker, so it keeps monotonicity
                            while still using the click to set the scale. The steel-man.
* `feat_sim`             -- negative Euclidean distance to the seed in the 22-feature space,
                            z-scored within the ROI. The richest cheap version of the idea.

Protocol, per the gates in `2026-09-01-one-click-retrieval-literature.md` §7:
the seed is **dropped from the evaluation ground truth** (a click scores itself perfectly by
construction), 5 seeds per ROI, and arms are judged on the **worst** seed, not the median.

Reads .cache_tail/fp_filter_features.csv. Writes results/click_darkness.csv.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from fp_filter_probe import FEATURES
from fp_filter_domain import tail_rematched

N_SEEDS = 5


def main():
    df = pd.read_csv(".cache_tail/fp_filter_features.csv")
    _images, ann = ds.load_annotations()
    rows = []

    for fn, sub in df.groupby("file_name"):
        sub = sub.reset_index(drop=True)
        gt_all = ds.image_annotations(ann, fn)
        radius = ev.radius_px(ds.roi_mpp(f"images/{fn}"))

        X = np.nan_to_num(sub[FEATURES].to_numpy(float), nan=0.0, posinf=0.0, neginf=0.0)
        Z = (X - X.mean(0)) / np.where(X.std(0) > 1e-9, X.std(0), 1.0)
        od = np.nan_to_num(sub["od51"].to_numpy(float), nan=-1e9)
        blob = sub["score"].to_numpy(float)

        pos = np.flatnonzero(sub["y"].to_numpy() == 1)
        if len(pos) < 2:
            continue
        rng = np.random.default_rng(0)
        seeds = rng.choice(pos, size=min(N_SEEDS, len(pos)), replace=False)

        for si, s in enumerate(seeds):
            # the clicked cell is removed from the evaluation ground truth and from the
            # candidate list, so no arm can score its own seed
            ann_id = int(sub.at[s, "matched_ann_id"])
            gt = gt_all[gt_all["ann_id"] != ann_id]
            keep = np.ones(len(sub), bool); keep[s] = False
            cand = sub[keep].reset_index(drop=True)

            d_seed = od[s]
            arms = {
                "blob_abs": blob[keep],
                "od51_abs": od[keep],
                "od51_sim": -np.abs(od[keep] - d_seed),
                "od51_onesided": -np.maximum(0.0, d_seed - od[keep]) * 1e6 + od[keep],
                "feat_sim": -np.linalg.norm(Z[keep] - Z[s], axis=1),
            }
            for name, sc in arms.items():
                for r in tail_rematched(cand, sc.astype(float), gt, radius):
                    rows.append(dict(file_name=fn, seed=si, seed_ann_id=ann_id,
                                     ranker=name, **r))
        print(f"{fn}: {len(seeds)} seeds done", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv("results/click_darkness.csv", index=False)
    pd.set_option("display.width", 250)
    for q in (0.95, 1.0):
        w = (out[out["quantile"] == q]
             .groupby(["file_name", "ranker"])["depth"].max().unstack())  # worst of 5 seeds
        print(f"\n=== WORST-of-{N_SEEDS}-seeds depth to {q:.0%} of captured ===")
        print(w.to_string())


if __name__ == "__main__":
    main()
