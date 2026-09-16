"""
    Script form of `production_pipeline_demo.ipynb`: one clicked annotation on one ROI in,
    ranked candidate detections out, no notebook required.

    Every actual pipeline step below is a direct call into `midog_utils` -- this file only
    wires them in the order the demo notebook does. The one piece of logic that isn't already
    a midog_utils function is `select_annotation`: deciding *which* ground-truth annotation on
    an ROI becomes the click. Everything downstream of that (gating/recentring the template,
    searching, ranking, bucketing against ground truth) is unchanged production code:
    `seed_selection.build_seed` and `production.run_production_pipeline`.

    Run directly for a quick look -- see the example invocations just below this docstring.
"""

# Example invocations (anaconda python required -- see the repo's own notes on why):
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_pipeline.py 403.tiff
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_pipeline.py 245.tiff
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_pipeline.py 403.tiff --rank-key chromatin_od
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_pipeline.py 245.tiff --ann-id 6274        # pin the click to one specific annotation
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_pipeline.py 403.tiff --seed-index 1 --top-n 20
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_pipeline.py 403.tiff --save-fig /tmp/403_overlay.png

from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from midog_utils import channels as ch
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from midog_utils import production as prod
from midog_utils import seed_selection as ss
from midog_utils import viz

DEFAULT_DB_PATH = ROOT / "databases" / "MIDOG++.json"
DEFAULT_IMAGES_DIR = ROOT / "images" / "extra_valid"


def select_annotation(gt_mitotic: pd.DataFrame, gray_inv: np.ndarray, roi_shape: tuple, image_id: int, *, ann_id: int | None = None, seed_index: int = 0) -> ss.Seed:
    """
        Pick which ground-truth mitotic figure on this ROI becomes the click, then gate
        and recentre it into a template via `seed_selection.build_seed`.

        gt_mitotic (pd.DataFrame): this ROI's mitotic ground-truth rows, the draw pool.
        gray_inv (np.ndarray): the ROI's inverted-grayscale structural channel.
        roi_shape (tuple): the ROI array's (height, width, channels) shape.
        image_id (int): this ROI's MIDOG++ image id; seeds the random draw.
        ann_id (int or None): pin the click to this annotation id; None draws randomly.
        seed_index (int): rng draw index used when ann_id is None.

        Returns ss.Seed: the gated/recentred template -- click point, template point, size,
        and retry count.

        `ann_id=None` (default) reproduces the demo notebook's convention: draw one via
        rejection sampling seeded on `[seed_index, image_id]`. Passing `ann_id` instead pins
        the click to that specific annotation -- the click-to-verify case, where a point was
        already chosen -- by narrowing `build_seed`'s draw pool to that one row; `build_seed`
        still runs its own agreement/border/gate checks against it and raises if that
        annotation doesn't survive them (e.g. too close to the ROI edge, or its component
        fails the Otsu gate), the same as a rejected candidate in the random-draw path.
    """
    pool = gt_mitotic
    if ann_id is not None:
        pool = gt_mitotic[gt_mitotic["ann_id"] == ann_id].reset_index(drop=True)
        if pool.empty:
            raise ValueError(f"ann_id={ann_id} is not a mitotic annotation on this ROI")

    rng = np.random.default_rng([seed_index, image_id]) # random number generator
    return ss.build_seed(pool, gray_inv, rng, roi_shape)


def run_pipeline_on_roi(fn: str, *, ann_id: int | None = None, seed_index: int = 0, rank_key: str = "tm_score", images_dir: Path | str = DEFAULT_IMAGES_DIR, db_path: Path | str = DEFAULT_DB_PATH) -> dict:
    """
        Run the full production pipeline on one ROI for one annotation: load the ROI,
        pick the click, search, rank, and score against ground truth.

        fn (str): ROI filename, e.g. "403.tiff".
        ann_id (int or None): pin the click to this annotation id; None draws one via seed_index.
        seed_index (int): rng draw index used when ann_id is None.
        rank_key (str): which axis ranks detections, "tm_score" or "chromatin_od".
        images_dir (Path or str): directory holding the ROI TIFFs.
        db_path (Path or str): path to the MIDOG++.json annotations file.

        Returns dict: the ROI, mpp, seed, match radius, ranked detections, ground truth, and
        pipeline info for this run.
    """
    images, annotations = ds.load_annotations(str(db_path))
    meta_ix = images.set_index("file_name")[["image_id", "tumor_type"]]
    image_id = int(meta_ix.loc[fn, "image_id"])

    rgb = ds.load_roi(f"{images_dir}/{fn}")
    mpp = ds.roi_mpp(f"{images_dir}/{fn}")
    roi_shape = rgb.shape
    match_radius = ev.radius_px(mpp) # currently redundant, already returned by info as well.
    gray_inv = ch.to_channel(rgb, "gray_inverted")

    gt = ds.image_annotations(annotations, fn)
    gt_mitotic = gt[gt["category_id"] == ds.MITOTIC]
    seed = select_annotation(gt_mitotic, gray_inv, roi_shape, image_id, ann_id=ann_id, seed_index=seed_index)

    # Ground truth stays keyed to the click; only the template moved (D8 item 3).
    gt_eval = gt[gt["ann_id"] != seed.ann_id].reset_index(drop=True)

    detections, info = prod.run_production_pipeline(rgb, seed, mpp, rank_key=rank_key)
    det_out, gt_out = ev.bucket_detections(detections, gt_eval, match_radius)

    return dict(fn=fn, rank_key=rank_key, rgb=rgb, mpp=mpp, seed=seed, match_radius=match_radius, detections=det_out, gt_eval=gt_out, info=info)


def summarize(result: dict, top_n: int = 10) -> dict:
    """
        Compute the precision/recall-at-top_n summary row for one pipeline result,
        matching the demo notebook's summary table.

        result (dict): one `run_pipeline_on_roi` output.
        top_n (int): how many top-ranked detections to score precision/recall over.

        Returns dict: one summary row -- roi, rank_key, peak/detection counts,
        precision_at_top_n, and recall_at_top_n.
    """
    det = result["detections"]
    n_gt = int((result["gt_eval"]["category_id"] == ds.MITOTIC).sum())
    top = det.head(top_n)
    tp = int((top["bucket"] == ev.HUMAN_CORRECT_LABEL).sum())
    return dict(roi=result["fn"], rank_key=result["rank_key"], n_peaks=result["info"]["n_peaks"], max_peaks_binding=result["info"]["max_peaks_binding"], n_detections=result["info"]["n_detections"], **{f"precision_at_{top_n}": round(tp / min(top_n, len(det)), 3) if len(det) else float("nan"), f"recall_at_{top_n}": round(tp / n_gt, 3) if n_gt else float("nan")})


def visualize(result: dict, ax=None, downsample: int = 4, top_n: int = 10):
    """
        Draw the demo notebook's overlay: ground truth, ranked detections, and the D8
        template box, on one ROI.

        result (dict): one `run_pipeline_on_roi` output.
        ax (matplotlib Axes or None): axes to draw on; None creates a new figure.
        downsample (int): display downsampling factor for the overlay image.
        top_n (int): how many top-ranked detections to draw.

        Returns matplotlib Axes: the axes the overlay was drawn on.
    """
    import matplotlib.pyplot as plt

    seed = result["seed"]
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))

    viz.overlay(result["rgb"], result["gt_eval"], result["detections"], seed_xy=seed.click_xy, tpl_xy=seed.template_xy, ax=ax, downsample=downsample, top_n=top_n, radius=result["match_radius"], title=f"{result['fn']} -- {result['rank_key']} (top {top_n})")

    # base_size is always odd (template_match._odd), so this box's pixel centre is exactly
    # seed.template_xy -- seed_selection.tightened_template_box's own half-pixel convention.
    half = (seed.base_size - 1) // 2
    tx, ty = seed.template_xy
    box = (ty - half, ty - half + seed.base_size, tx - half, tx - half + seed.base_size)
    viz.draw_box(ax, box, "#00e5ff", lw=1.2, scale=1.0 / downsample)
    return ax


def main():
    parser = argparse.ArgumentParser(description=textwrap.dedent(__doc__ or ""), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("roi", help="ROI filename, e.g. 403.tiff")
    parser.add_argument("--ann-id", type=int, default=None, help="pin the click to this ground-truth annotation id (default: random draw via --seed-index, matching the demo notebook)")
    parser.add_argument("--seed-index", type=int, default=0, help="rng draw index used when --ann-id is not given")
    parser.add_argument("--rank-key", default="tm_score", choices=sorted(prod.AXES))
    parser.add_argument("--images-dir", default=str(DEFAULT_IMAGES_DIR))
    parser.add_argument("--db-path", default=str(DEFAULT_DB_PATH))
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--save-fig", default=None, help="path to save the overlay figure to, instead of plt.show()")
    args = parser.parse_args()

    result = run_pipeline_on_roi(args.roi, ann_id=args.ann_id, seed_index=args.seed_index, rank_key=args.rank_key, images_dir=args.images_dir, db_path=args.db_path)

    seed = result["seed"]
    print(f"{args.roi}: seed ann {seed.ann_id}, base_size {seed.base_size}, click-to-template offset {seed.offset_px:.2f}px, n_retries {seed.n_retries}")
    print(pd.Series(summarize(result, top_n=args.top_n)).to_frame().T.to_string(index=False))

    ax = visualize(result, top_n=args.top_n)
    if args.save_fig:
        ax.figure.savefig(args.save_fig, dpi=150, bbox_inches="tight")
        print(f"saved figure to {args.save_fig}")
    else:
        import matplotlib.pyplot as plt
        plt.show()


if __name__ == "__main__":
    main()
