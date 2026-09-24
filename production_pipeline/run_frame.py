"""
    Script form of the annotation-frame path: one ROI and one click in, the annotation frame out --
    the click, the template box it was refined to, and every ranked prediction -- scored against
    MIDOG++ ground truth when the click is a MIDOG annotation. The frame-path counterpart of
    `run_pipeline.py`, which stays as it is.

    Every step is a direct call into `midog_utils`: `dataset` loads the ROI, and `annotation_frame`
    builds the input frame (`click_frame`, `human_seed_frame` or `draw_seed_rows`), runs it
    (`run_roi`) and scores it (`evaluate`). This file only chooses which builder to use.

    Run directly for a quick look -- see the example invocations just below this docstring.
"""

# Example invocations (anaconda python required -- see the repo's own notes on why):
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff                     # validation draw, seed 0
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --seed-index 2
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --ann-id 14499       # pin the click to one MIDOG annotation
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --click 1411 3551    # a raw click: predictions only, no scoring
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --rank-key chromatin_od --refinement-channel hematoxylin_od
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 322.tiff --save-csv /tmp/322_frame.csv

from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from midog_utils import annotation_frame as af
from midog_utils import dataset as ds
from midog_utils import production as prod

DEFAULT_DB_PATH = ROOT / "databases" / "MIDOG++.json"
DEFAULT_IMAGES_DIR = ROOT / "images" / "extra_valid"
BUDGETS = (10, 20, 30)


def roi_path(images_dir: Path | str, fn: str) -> Path:
    """
        Find an ROI's TIFF in images_dir, or in its testing_set subfolder where the 35 non-original ROIs live.

        images_dir (Path or str): directory holding the ROI TIFFs.
        fn (str): ROI filename, e.g. "300.tiff".

        Returns Path: the TIFF's path. Raises FileNotFoundError when neither location has it.
    """
    for path in (Path(images_dir) / fn, Path(images_dir) / "testing_set" / fn):
        if path.exists():
            return path
    raise FileNotFoundError(f"{fn} is in neither {images_dir} nor its testing_set subfolder")


def build_seed_frame(fn: str, rgb, gt: pd.DataFrame, image_id: int, *, click: tuple | None = None, click_id: int = 0, edge_origin: bool = False, ann_id: int | None = None, seed_index: int = 0) -> tuple[pd.DataFrame, bool]:
    """
        Build run_roi's one-row input frame from a raw click, a pinned MIDOG annotation, or the validation draw.

        fn (str): ROI filename.
        rgb (np.ndarray): the ROI's RGB image; the validation draw reads it.
        gt (pd.DataFrame): this ROI's annotations, all categories.
        image_id (int): this ROI's MIDOG++ image id; seeds the validation draw.
        click (tuple or None): a raw (x, y) click; takes precedence over ann_id and seed_index.
        click_id (int): the raw click's object_id.
        edge_origin (bool): the raw click uses pixel-edge coordinates.
        ann_id (int or None): pin the click to this mitotic annotation.
        seed_index (int): validation-draw seed, used when neither click nor ann_id is given.

        Returns tuple[pd.DataFrame, bool]: the input frame, and whether its object_id is a MIDOG ann_id that evaluate can score.
    """
    if click is not None:
        return af.click_frame(fn, click[0], click[1], click_id, edge_origin=edge_origin), False
    if ann_id is not None:
        rows = gt[(gt["ann_id"] == ann_id) & (gt["category_id"] == ds.MITOTIC)]
        if rows.empty:
            raise ValueError(f"ann_id={ann_id} is not a mitotic annotation on {fn}")
        return af.human_seed_frame(rows.iloc[0], fn), True
    # Draws are sequential -- each seed excludes the clicks earlier seeds accepted -- so seed k is the last of range(k + 1).
    return af.draw_seed_rows(rgb, gt, image_id, seed_indices=range(seed_index + 1))[-1][0], True


def run_frame_on_roi(fn: str, *, click: tuple | None = None, click_id: int = 0, edge_origin: bool = False, ann_id: int | None = None, seed_index: int = 0, rank_key: str = "tm_score", refinement_channel: str | None = "gray_inverted", images_dir: Path | str = DEFAULT_IMAGES_DIR, db_path: Path | str = DEFAULT_DB_PATH) -> dict:
    """
        Run the annotation-frame path on one ROI for one click: load, build the input frame, run_roi, and evaluate when the click is a MIDOG annotation.

        fn (str): ROI filename, e.g. "300.tiff".
        click (tuple or None): a raw (x, y) click.
        click_id (int): the raw click's object_id.
        edge_origin (bool): the raw click uses pixel-edge coordinates.
        ann_id (int or None): pin the click to this mitotic annotation.
        seed_index (int): validation-draw seed, used when neither click nor ann_id is given.
        rank_key (str): which axis ranks the predictions, "tm_score" or "chromatin_od".
        refinement_channel (str or None): the Otsu gate's channel, or None for the untightened 51 px template.
        images_dir (Path or str): directory holding the ROI TIFFs.
        db_path (Path or str): path to the MIDOG++.json annotations file.

        Returns dict: fn, the output frame and info, and evaluate's scored predictions, gt_out and metrics (None for a raw click).
    """
    images, annotations = ds.load_annotations(str(db_path))
    image_id = int(images.set_index("file_name").loc[fn, "image_id"])

    path = roi_path(images_dir, fn)
    rgb = ds.load_roi(path)
    mpp = ds.roi_mpp(path)
    gt = ds.image_annotations(annotations, fn)

    seed_frame, scorable = build_seed_frame(fn, rgb, gt, image_id, click=click, click_id=click_id, edge_origin=edge_origin, ann_id=ann_id, seed_index=seed_index)
    frame, info = af.run_roi(rgb, seed_frame, mpp, rank_key=rank_key, refinement_channel=refinement_channel)
    scored, gt_out, metrics = af.evaluate(frame, gt, mpp, rgb.shape, budgets=BUDGETS) if scorable else (None, None, None)
    return dict(fn=fn, frame=frame, info=info, scored=scored, gt_out=gt_out, metrics=metrics)


def summarize(result: dict) -> dict:
    """
        One summary row for a run: the seed and template, the detection count, and tp/precision at each budget when scored.

        result (dict): one `run_frame_on_roi` output.

        Returns dict: the summary row; precision_at_k = tp_at_k / k, as the reference driver computes it.
    """
    info, metrics = result["info"], result["metrics"]
    row = dict(roi=result["fn"], seed_object_id=info["seed_object_id"], refinement_channel=info["refinement_channel"], refined=info["refined"], refusal_reason=info["refusal_reason"], base_size=info["base_size"], rank_key=info["rank_key"], n_detections=info["n_detections"])
    if metrics is not None:
        for k in BUDGETS:
            row[f"tp_at_{k}"] = metrics[f"tp_at_{k}"]
            row[f"precision_at_{k}"] = round(metrics[f"tp_at_{k}"] / k, 3)
    return row


def main():
    parser = argparse.ArgumentParser(description=textwrap.dedent(__doc__ or ""), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("roi", help="ROI filename, e.g. 300.tiff")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--click", type=float, nargs=2, metavar=("X", "Y"), default=None, help="a raw click in full-resolution ROI pixels; predictions only, no ground-truth scoring")
    source.add_argument("--ann-id", type=int, default=None, help="pin the click to this MIDOG mitotic annotation")
    parser.add_argument("--click-id", type=int, default=0, help="object_id for a --click row")
    parser.add_argument("--edge-origin", action="store_true", help="--click uses pixel-edge coordinates (COCO / canvas); 0.5 is subtracted")
    parser.add_argument("--seed-index", type=int, default=0, help="validation-draw seed, used when neither --click nor --ann-id is given")
    parser.add_argument("--rank-key", default="tm_score", choices=sorted(prod.AXES))
    parser.add_argument("--refinement-channel", default="gray_inverted", choices=("gray_inverted", "hematoxylin_od", "none"), help="the Otsu gate's channel; none = the untightened 51 px template")
    parser.add_argument("--images-dir", default=str(DEFAULT_IMAGES_DIR))
    parser.add_argument("--db-path", default=str(DEFAULT_DB_PATH))
    parser.add_argument("--show", type=int, default=5, help="how many predictions to print after the seed rows")
    parser.add_argument("--save-csv", default=None, help="path to write the full output frame to (read back with float_precision='round_trip')")
    args = parser.parse_args()

    channel = None if args.refinement_channel == "none" else args.refinement_channel
    click = tuple(args.click) if args.click is not None else None
    result = run_frame_on_roi(args.roi, click=click, click_id=args.click_id, edge_origin=args.edge_origin, ann_id=args.ann_id, seed_index=args.seed_index, rank_key=args.rank_key, refinement_channel=channel, images_dir=args.images_dir, db_path=args.db_path)

    info = result["info"]
    print(f"{args.roi}: seed {info['seed_object_id']} at {info['click_center']}, template {info['base_size']} px at {info['template_center']}, offset {info['offset_px']:.2f}px, refined={info['refined']}, refusal_reason={info['refusal_reason']}")
    print(pd.Series(summarize(result)).to_frame().T.to_string(index=False))
    if result["metrics"] is None:
        print("raw click: not scored (its object_id is not a MIDOG ann_id)")

    frame = result["frame"]
    n_seed_rows = int((frame["class"] != af.PREDICTION).sum())
    print(frame.head(n_seed_rows + args.show).drop(columns="current_time").to_string())
    if args.save_csv:
        frame.to_csv(args.save_csv, index=False)
        print(f"saved {len(frame)} rows to {args.save_csv}")


if __name__ == "__main__":
    main()
