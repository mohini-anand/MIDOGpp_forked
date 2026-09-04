"""Does raising the Otsu foreground threshold help or hurt end-to-end detection?

`Research Logs/2026-09-03-bbox-threshold-sweep.md` found that raising `tighten_box_otsu`'s
effective foreground threshold (`method="headroom"`, `threshold = otsu_thresh +
headroom_frac*(255-otsu_thresh)`) un-merges a majority of the multi-nucleus connected
components that plague dense-cellular domains (lymphosarcoma, mast cell tumor) -- but that
analysis measured seed-selection *geometry* only (box size/shape/solidity) and explicitly
left "downstream template-matching/discrimination quality... completely untested." This
script closes that gap: an actual end-to-end `find_and_suppress` run, comparing
`bbox_method="binary"` (production default) against `bbox_method="headroom",
bbox_headroom_frac=0.15` (the sweep's recommended middle ground) and `bbox_method="multiotsu"`
(the repo's other existing opt-in) on real detection metrics, not box geometry.

Design, following `Research Logs/2026-09-02-next-steps-plan.md`'s established, already-audited
practice (its "0. Terms" and "6. Statistical power" sections):

- **Population**: the 10 "decision-grade" ROIs currently downloaded (`n_mitotic >= 15`,
  `tm_variant_sweep.DECISION_MIN_MITOTIC`) -- 094, 201, 202, 245, 246, 300, 301, 402, 459,
  548.tiff. Computed fresh each run from `dataset.load_annotations()` against `images/`, not
  hardcoded, so a changed set of downloaded ROIs is reflected automatically.
- **5 seeds per ROI**: `rng = np.random.default_rng(seed_index)` for `seed_index in range(5)`,
  a *fresh* generator per (ROI, seed_index, arm) call -- not one generator advanced across
  arms -- so each arm draws against its own filtered candidate pool from the same nominal
  draw index, not a state-drifted one.
- **`fs.FSConfig()` defaults, unmodified** (score_threshold=0.5, n_angles=1, flips=(False,),
  single-scale -- see the dataclass itself, `find_and_suppress.py`).
- **`run_one_image(..., run_baselines=False)`** for each of the 10 x 5 x 3 = 150 (ROI,
  seed_index, arm) triples. No baseline (nucleus-blob/random) comparison here -- only the
  three `bbox_method` arms matter for this question.
- **`recall_at_budget` at budgets (100, 250, 500)**, matching the Sep-2 plan's headline
  metric, derived via `evaluate.recall_at_k(det_out["bucket"], n_gt_mitotic_eval, k=budget)`
  on the already-bucketed `detections` frame `run_one_image` returns -- not a reimplementation;
  algebraically identical to `compare.evaluate_arms`'s own `tp_at_budget / n_mit` with
  `delivered = min(budget, len(list))`.
- **`seed_ann_id` recorded for every run.** `pick_seed`'s candidate pool depends on
  `bbox_method` (an annotation accepted under `binary`'s foreground filter can be excluded
  under `headroom`/`multiotsu`, or vice versa), so the same `rng` draw at the same (ROI,
  seed_index) can select a *different* underlying annotation across arms -- expected, not a
  bug, but it must be visible, not silently averaged away. See the Research Log for the
  seed-pool-mismatch analysis this enables.

Writes `results/bbox_headroom_end_to_end.csv` -- one row per (ROI, seed_index, arm) triple
(150 rows), the full metrics dict `run_one_image` returns plus the three `recall_at_budget_*`
/ `budget_delivered_*` columns. Does not aggregate, test, or analyse -- see
`Research Logs/2026-09-03-bbox-headroom-end-to-end.md` for the ROI-level median + sign-test
analysis built on top of this raw file.

Runtime: ~150 calls x ~5-7s each (CPU-only) -- about 15-20 minutes total.
"""
from __future__ import annotations

import time
import traceback

import numpy as np
import pandas as pd

from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from midog_utils import experiment as exp
from midog_utils import find_and_suppress as fs

DECISION_MIN_MITOTIC = 15  # matches tm_variant_sweep.DECISION_MIN_MITOTIC, applied by count

ARMS = [
    {"bbox_method": "binary", "bbox_headroom_frac": None},
    {"bbox_method": "headroom", "bbox_headroom_frac": 0.15},
    {"bbox_method": "multiotsu", "bbox_headroom_frac": None},
]

BUDGETS = (100, 250, 500)
N_SEEDS = 5

OUT_PATH = "results/bbox_headroom_end_to_end.csv"
ERR_PATH = "results/bbox_headroom_end_to_end_errors.csv"


def decision_grade_rois(images: pd.DataFrame, annotations: pd.DataFrame, images_dir="images"):
    """Every downloaded ROI with >= DECISION_MIN_MITOTIC category-1 annotations.

    Computed fresh, not hardcoded -- see this module's docstring.
    """
    from pathlib import Path

    on_disk = {p.name for p in Path(images_dir).glob("*.tiff")}
    counts = (
        annotations[annotations["category_id"] == ds.MITOTIC]
        .groupby("file_name").size().rename("n_mitotic")
    )
    df = images[images["file_name"].isin(on_disk)].merge(
        counts, left_on="file_name", right_index=True, how="left"
    )
    df["n_mitotic"] = df["n_mitotic"].fillna(0).astype(int)
    dg = df[df["n_mitotic"] >= DECISION_MIN_MITOTIC].sort_values("file_name")
    return dg["file_name"].tolist()


def run(out_path=OUT_PATH, err_path=ERR_PATH):
    images, annotations = ds.load_annotations()
    rois = decision_grade_rois(images, annotations)
    print(f"Decision-grade ROIs ({len(rois)}): {rois}")

    rows = []
    errors = []
    t_start = time.time()
    n_total = len(rois) * N_SEEDS * len(ARMS)
    n_done = 0

    for roi in rois:
        for seed_idx in range(N_SEEDS):
            for arm in ARMS:
                # Fresh generator per (roi, seed_idx, arm): each arm draws against its own
                # filtered pool from the same nominal draw index, not a state-drifted one.
                rng = np.random.default_rng(seed_idx)
                t0 = time.time()
                try:
                    out = exp.run_one_image(
                        roi, annotations, cfg=fs.FSConfig(), rng=rng, run_baselines=False,
                        bbox_method=arm["bbox_method"],
                        bbox_headroom_frac=arm["bbox_headroom_frac"],
                    )
                    dt = time.time() - t0
                    metrics = dict(out["metrics"][0])
                    det_out = out["detections"]
                    n_mit = metrics["n_gt_mitotic_eval"]
                    bucket = det_out["bucket"].to_numpy()
                    for b in BUDGETS:
                        metrics[f"recall_at_budget_{b}"] = ev.recall_at_k(bucket, n_mit, k=b)
                        metrics[f"budget_delivered_{b}"] = min(b, len(det_out))
                    metrics["seed_index"] = seed_idx
                    metrics["roi"] = roi
                    metrics["arm"] = arm["bbox_method"]
                    metrics["wall_s"] = round(dt, 3)
                    rows.append(metrics)
                except Exception as e:
                    dt = time.time() - t0
                    err = {
                        "roi": roi, "seed_index": seed_idx, "arm": arm["bbox_method"],
                        "bbox_headroom_frac": arm["bbox_headroom_frac"],
                        "error": str(e), "wall_s": round(dt, 3),
                    }
                    errors.append(err)
                    print("ERROR:", err)
                    traceback.print_exc()
                n_done += 1
                elapsed = time.time() - t_start
                print(f"[{n_done}/{n_total}] {roi} seed={seed_idx} arm={arm['bbox_method']} "
                      f"dt={dt:.2f}s elapsed={elapsed:.1f}s")

    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    print(f"Saved {len(df)} rows to {out_path}")
    print(f"Total wall clock: {time.time() - t_start:.1f}s")
    if errors:
        print(f"\n{len(errors)} ERRORS:")
        for e in errors:
            print(e)
        pd.DataFrame(errors).to_csv(err_path, index=False)
    return df


if __name__ == "__main__":
    run()
