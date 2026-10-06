"""
    The pipeline's single entry point: one ROI and one human seed in -- a click, a drawn box, or a
    MIDOG++ annotation, chosen with --source -- and the annotation frame out: the seed row, the
    template box it was refined to, every ranked prediction, and each prediction's box refined
    the same way (unless --no-refine-predictions), scored against MIDOG++ ground truth when the
    seed is a MIDOG annotation, and plotted on request.

    Every step is a direct call into `midog_utils`: `dataset` loads the ROI, and `annotation_frame`
    builds the input frame (`click_frame`, `box_frame`, `human_seed_frame` or `draw_seed_rows`),
    runs it (`run_roi`) and scores it (`evaluate`). This file only chooses which builder to use.
    A drawn box's centre is the click and its size the Otsu gate's window and fallback template;
    a click or MIDOG seed always uses 51 px. Each prediction is refined with the seed's gate channel
    in a 51 px window (its own box) and stored as a PREDICTION_REFINED row, whatever the human later
    decides. `run_pipeline.py` is deprecated in favour of this file.

    Every run is one iteration of a session (`annotation_frame.run_iteration`). Without
    --session-csv it is the first, iteration 0, and blanks nothing. With it, every box the earlier
    iterations stored (their refined seeds and shown predictions) is blanked in the search, the new
    seed's box may not contain a stored box's centre (unless --allow-overlap), and --save-session
    writes the session with this iteration's seed and first --k shown predictions appended.

    Run directly for a quick look -- see the example invocations just below this docstring.
"""

# Example invocations (anaconda python required -- see the repo's own notes on why):
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --source midog                        # validation draw, seed 0
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --source midog --seed-index 2
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --source midog --ann-id 14499          # pin the seed to one MIDOG annotation
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --source click --click 1411 3551       # a raw click: predictions only, no scoring
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --source box --box 1391 3531 41 41     # a drawn box, pixel-edge x y w h: predictions only
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --source midog --ann-id 14499 --rank-key chromatin_od --refinement-channel hematoxylin_od
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --source midog --ann-id 14499 --save-fig /tmp/300_frame.png
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 322.tiff --source midog --save-csv /tmp/322_frame.csv
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --source midog --no-refine-predictions   # predictions without PREDICTION_REFINED rows
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --source midog --ann-id 14499 --save-session /tmp/300_session.csv                                # iteration 0 of a session
#   /Users/mohinianand/anaconda3/bin/python3 production_pipeline/run_frame.py 300.tiff --source midog --ann-id 14578 --session-csv /tmp/300_session.csv --save-session /tmp/300_session.csv   # iteration 1: blanks iteration 0's stored boxes

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
SOURCE_OPTIONS = {"click": ("click", "object_id", "edge_origin"), "box": ("box", "object_id"), "midog": ("ann_id", "seed_index")}  # the seed arguments each source takes
REQUIRED_OPTION = {"click": "click", "box": "box"}


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


def check_seed_args(source: str, *, click: tuple | None = None, box: tuple | None = None, object_id: int | None = None, edge_origin: bool | None = None, ann_id: int | None = None, seed_index: int | None = None) -> None:
    """
        Check a seed source and its arguments before any loading: the source is known, it gets only its own arguments and
        every required one, and a box is a valid drawn box. None means "not given".

        source (str): "click", "box" or "midog".
        click (tuple or None): a raw (x, y) click; source "click" only, required there.
        box (tuple or None): a drawn (x, y, w, h) box in pixel-edge coordinates; source "box" only, required there.
        object_id (int or None): the click's or box's object_id; sources "click" and "box".
        edge_origin (bool or None): the click uses pixel-edge coordinates; source "click" only.
        ann_id (int or None): pin the seed to this mitotic annotation; source "midog" only.
        seed_index (int or None): validation-draw seed, when ann_id is not given; source "midog" only.

        Returns None. Raises ValueError naming the first problem.
    """
    if source not in SOURCE_OPTIONS:
        raise ValueError(f"source must be one of {list(SOURCE_OPTIONS)}, got {source!r}")
    given = {k for k, v in dict(click=click, box=box, object_id=object_id, edge_origin=edge_origin, ann_id=ann_id, seed_index=seed_index).items() if v is not None}
    foreign = sorted(given - set(SOURCE_OPTIONS[source]))
    if foreign:
        raise ValueError(f"source {source!r} does not take {', '.join('--' + k.replace('_', '-') for k in foreign)}")
    required = REQUIRED_OPTION.get(source)
    if required and required not in given:
        raise ValueError(f"source {source!r} needs --{required}")
    if ann_id is not None and seed_index is not None:
        raise ValueError("--ann-id and --seed-index are exclusive: a pinned annotation is not drawn")
    if seed_index is not None and seed_index < 0:
        raise ValueError(f"--seed-index must be >= 0, got {seed_index}")
    if box is not None:
        af.box_frame("", *box, 0)


def build_seed_frame(fn: str, rgb, gt: pd.DataFrame, image_id: int, source: str, *, click: tuple | None = None, box: tuple | None = None, object_id: int | None = None, edge_origin: bool | None = None, ann_id: int | None = None, seed_index: int | None = None) -> tuple[pd.DataFrame, bool]:
    """
        Build run_roi's one-row input frame for a seed source: a raw click, a drawn box, or a MIDOG annotation (pinned, or
        the validation draw).

        fn (str): ROI filename.
        rgb (np.ndarray): the ROI's RGB image; the validation draw reads it.
        gt (pd.DataFrame): this ROI's annotations, all categories.
        image_id (int): this ROI's MIDOG++ image id; seeds the validation draw.
        source (str): "click", "box" or "midog".
        click (tuple or None): a raw (x, y) click, pixel-centre coordinates unless edge_origin.
        box (tuple or None): a drawn (x, y, w, h) box, pixel-edge coordinates.
        object_id (int or None): the click's or box's object_id; 0 when None.
        edge_origin (bool or None): the click uses pixel-edge coordinates.
        ann_id (int or None): pin the seed to this mitotic annotation.
        seed_index (int or None): validation-draw seed, used when ann_id is None; 0 when None.

        Returns tuple[pd.DataFrame, bool]: the input frame, and whether its object_id is a MIDOG ann_id that evaluate can score.
        Raises ValueError as check_seed_args does, or when ann_id is not a mitotic annotation on this ROI.
    """
    check_seed_args(source, click=click, box=box, object_id=object_id, edge_origin=edge_origin, ann_id=ann_id, seed_index=seed_index)
    if source == "click":
        return af.click_frame(fn, click[0], click[1], object_id or 0, edge_origin=bool(edge_origin)), False
    if source == "box":
        return af.box_frame(fn, *box, object_id or 0), False
    if ann_id is not None:
        rows = gt[(gt["ann_id"] == ann_id) & (gt["category_id"] == ds.MITOTIC)]
        if rows.empty:
            raise ValueError(f"ann_id={ann_id} is not a mitotic annotation on {fn}")
        return af.human_seed_frame(rows.iloc[0], fn), True
    # Draws are sequential -- each seed excludes the clicks earlier seeds accepted -- so seed k is the last of range(k + 1).
    return af.draw_seed_rows(rgb, gt, image_id, seed_indices=range((seed_index or 0) + 1))[-1][0], True


def read_session(path: Path | str) -> pd.DataFrame:
    """
        Read a session CSV written by --save-session, with round_trip float parsing, and validate it.

        path (Path or str): the CSV.

        Returns pd.DataFrame: the session (`annotation_frame.check_session`).
    """
    return af.check_session(pd.read_csv(path, float_precision="round_trip"))


def run_frame_on_roi(fn: str, source: str, *, click: tuple | None = None, box: tuple | None = None, object_id: int | None = None, edge_origin: bool | None = None, ann_id: int | None = None, seed_index: int | None = None, rank_key: str = "tm_score", refinement_channel: str | None = "gray_inverted", refine_predictions: bool = True, session: pd.DataFrame | None = None, k: int = af.DEFAULT_K, check_overlap: bool = True, images_dir: Path | str = DEFAULT_IMAGES_DIR, db_path: Path | str = DEFAULT_DB_PATH, keep_rgb: bool = False) -> dict:
    """
        Run one iteration of the annotation-frame path on one ROI for one seed: load, build the input frame, run_iteration
        (run_roi with the session's stored boxes blanked), and evaluate when the seed is a MIDOG annotation and this is
        iteration 0. Later iterations are not scored: evaluate would read run_roi's raw list, not the page shown, against
        ground truth that still holds the mitoses found earlier.

        fn (str): ROI filename, e.g. "300.tiff".
        source (str): "click", "box" or "midog".
        click (tuple or None): a raw (x, y) click.
        box (tuple or None): a drawn (x, y, w, h) box, pixel-edge coordinates.
        object_id (int or None): the click's or box's object_id; the session's next free id when None.
        edge_origin (bool or None): the click uses pixel-edge coordinates.
        ann_id (int or None): pin the seed to this mitotic annotation.
        seed_index (int or None): validation-draw seed, used when ann_id is None.
        rank_key (str): which axis ranks the predictions, "tm_score" or "chromatin_od".
        refinement_channel (str or None): the Otsu gate's channel, or None for no tightening: the template is the human's box (51 px for clicks and MIDOG).
        refine_predictions (bool): also refine every prediction's box the way the seed's was and append it as a PREDICTION_REFINED row.
        session (pd.DataFrame or None): the earlier iterations' stored rows (`read_session`); None for iteration 0.
        k (int): how many predictions to show and store in the session.
        check_overlap (bool): raise when the new seed's box contains the centre of a stored box.
        images_dir (Path or str): directory holding the ROI TIFFs.
        db_path (Path or str): path to the MIDOG++.json annotations file.
        keep_rgb (bool): keep the ROI image in the result, for visualize (about 93 MB on 300.tiff).

        Returns dict: fn, source, the output frame (every survivor) and info, the new session, evaluate's scored predictions, gt_out and metrics (None unless scorable at iteration 0), and rgb when keep_rgb.
    """
    check_seed_args(source, click=click, box=box, object_id=object_id, edge_origin=edge_origin, ann_id=ann_id, seed_index=seed_index)
    images, annotations = ds.load_annotations(str(db_path))
    image_id = int(images.set_index("file_name").loc[fn, "image_id"])

    path = roi_path(images_dir, fn)
    rgb = ds.load_roi(path)
    mpp = ds.roi_mpp(path)
    gt = ds.image_annotations(annotations, fn)

    if source in ("click", "box") and object_id is None:
        object_id = af.next_object_id(session)
    seed_frame, scorable = build_seed_frame(fn, rgb, gt, image_id, source, click=click, box=box, object_id=object_id, edge_origin=edge_origin, ann_id=ann_id, seed_index=seed_index)
    new_session, frame, info = af.run_iteration(rgb, session, mpp, seed_frame, k=k, check_overlap=check_overlap, rank_key=rank_key, refinement_channel=refinement_channel, refine_predictions=refine_predictions)
    scored, gt_out, metrics = af.evaluate(frame, gt, mpp, rgb.shape, budgets=BUDGETS) if scorable and info["iteration"] == 0 else (None, None, None)
    result = dict(fn=fn, source=source, frame=frame, info=info, session=new_session, scored=scored, gt_out=gt_out, metrics=metrics)
    if keep_rgb:
        result["rgb"] = rgb
    return result


def summarize(result: dict) -> dict:
    """
        One summary row for a run: the iteration, source, seed and template, the detection count, how many prediction boxes
        the gate refined (when refined), the page shown, and tp/precision at each budget when scored.

        result (dict): one `run_frame_on_roi` output.

        Returns dict: the summary row; precision_at_k = tp_at_k / k, as the reference driver computes it.
    """
    info, metrics = result["info"], result["metrics"]
    row = dict(roi=result["fn"], iteration=info["iteration"], source=result["source"], seed_object_id=info["seed_object_id"], refinement_channel=info["refinement_channel"], refined=info["refined"], refusal_reason=info["refusal_reason"], base_size=info["base_size"], rank_key=info["rank_key"], n_detections=info["n_detections"])
    if "n_predictions_refined" in info:
        row["n_predictions_refined"] = info["n_predictions_refined"]
    row.update(n_blanked_boxes=len(info["session_blanks"]), n_repeats_in_page=info["n_repeats_in_page"], n_shown=info["n_shown"])
    if metrics is not None:
        for k in BUDGETS:
            row[f"tp_at_{k}"] = metrics[f"tp_at_{k}"]
            row[f"precision_at_{k}"] = round(metrics[f"tp_at_{k}"] / k, 3)
    return row


def visualize(result: dict, ax=None, downsample: int = 4, top_n: int = 10):
    """
        Draw one run on its ROI: ground truth (when scored), the top-ranked predictions (from iteration 1 on, the page
        shown, repeats skipped), the seed (star) and template centre (plus), the gate window (white, when a gate ran), the
        template footprint (cyan), the boxes earlier iterations stored and this run blanked (orange) and, when the frame has
        them, the drawn predictions' refined boxes (magenta when the gate refined them, grey for the fallback).

        result (dict): one `run_frame_on_roi` output made with keep_rgb=True.
        ax (matplotlib Axes or None): axes for the overview; None creates a figure with the overview and a full-resolution zoom on the seed.
        downsample (int): display downsampling factor for the overlay image.
        top_n (int): how many top-ranked predictions to draw.

        Returns matplotlib Axes: the axes the overlay was drawn on.
    """
    import matplotlib.pyplot as plt

    from midog_utils import viz

    if result.get("rgb") is None:
        raise ValueError("visualize needs the ROI image: run run_frame_on_roi with keep_rgb=True")
    frame, info, rgb = result["frame"], result["info"], result["rgb"]
    seed_id = info["seed_object_id"]
    mitotic = frame[(frame["class"] == af.MITOTIC_GT) & (frame["object_id"] == seed_id)]
    seed = mitotic[mitotic["made_by"] == af.HUMAN].iloc[0]
    refined = mitotic[mitotic["made_by"] == af.ALGORITHM].iloc[0]

    page = info["iteration"] > 0  # from iteration 1 on, draw the page shown: repeats are skipped, so it is not the raw top
    if page:
        rows = frame[(frame["class"] == af.PREDICTION) & frame["object_id"].isin(info["shown_object_ids"])].sort_values("object_id", kind="stable")
    else:
        rows = (result["scored"] if result["scored"] is not None else frame[frame["class"] == af.PREDICTION]).sort_values("object_id", kind="stable")
    dets = pd.DataFrame({"cx": af.box_center(rows["x_tl"], rows["width"]).to_numpy(), "cy": af.box_center(rows["y_tl"], rows["height"]).to_numpy()})
    if "bucket" in rows:
        dets["bucket"] = rows["bucket"].to_numpy()
    gt = result["gt_out"] if result["gt_out"] is not None else pd.DataFrame({"cx": [], "cy": [], "category_id": []})

    how = "refined" if info["refined"] else (info["refusal_reason"] or "untightened")
    title = f"{result['fn']} -- iteration {info['iteration']}, {result['source']} seed {seed_id}, {info['base_size']} px template ({how}), {info['rank_key']} ({'page' if page else 'top'} {top_n})"
    if result["metrics"] is not None:
        title += f", tp@10 = {result['metrics']['tp_at_10']}"
    zoom = None
    if ax is None:
        _, (ax, zoom) = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw=dict(width_ratios=(3, 1)))
    viz.overlay(rgb, gt, dets, seed_xy=info["click_center"], tpl_xy=info["template_center"], ax=ax, downsample=downsample, top_n=top_n, radius=info["match_radius_px"], title=title)

    # Half-open (y0, y1, x0, x1) pixel boxes: the gate window on the rounded seed centre, as read_padded_patch reads it, and
    # the template footprint from the refined row (x_tl is the centre of its first pixel, width = base_size - 1).
    cx, cy = info["click_center"]
    ix, iy, half = int(round(cx)), int(round(cy)), af.gate_window_side(seed["width"], seed["height"]) // 2
    window = (iy - half, iy + half + 1, ix - half, ix + half + 1)
    footprint = (int(refined["y_tl"]), int(refined["y_tl"] + refined["height"]) + 1, int(refined["x_tl"]), int(refined["x_tl"] + refined["width"]) + 1)

    boxes = [(footprint, "#00e5ff", 1.2)]
    if info["refinement_requested"]:
        boxes.insert(0, (window, "white", 1.0))
    for box, color, lw in boxes:
        viz.draw_box(ax, box, color, lw=lw, scale=1.0 / downsample)
    for sq in info["session_blanks"].itertuples(index=False):
        bx, by, h = int(round(sq.cx)), int(round(sq.cy)), int(sq.side) // 2  # blank_squares' footprint
        viz.draw_box(ax, (by - h, by + h + 1, bx - h, bx + h + 1), "#ff9800", lw=0.9, scale=1.0 / downsample)
    drawn = info["shown_object_ids"][:top_n] if page else range(top_n)
    refined_preds = frame[(frame["class"] == af.PREDICTION_REFINED) & frame["object_id"].isin(drawn)]
    for pred in refined_preds.itertuples(index=False):
        box = (int(pred.y_tl), int(pred.y_tl + pred.height) + 1, int(pred.x_tl), int(pred.x_tl + pred.width) + 1)
        viz.draw_box(ax, box, "#ff40ff" if pred.refined else "#bdbdbd", lw=0.9, scale=1.0 / downsample)
    if zoom is not None:
        r = af.MAX_WINDOW_PX
        y0, y1, x0, x1 = max(iy - r, 0), min(iy + r + 1, rgb.shape[0]), max(ix - r, 0), min(ix + r + 1, rgb.shape[1])
        zoom.imshow(rgb[y0:y1, x0:x1], extent=(x0 - 0.5, x1 - 0.5, y1 - 0.5, y0 - 0.5))
        zoom.plot(cx, cy, marker="*", markersize=16, color="white", markeredgecolor="black", markeredgewidth=0.8)
        zoom.plot(*info["template_center"], marker="P", markersize=10, color="#00e5ff", markeredgecolor="black", markeredgewidth=0.7)
        for box, color, lw in boxes:
            viz.draw_box(zoom, box, color, lw=lw + 0.6)
        gate = f"gate window {window[1] - window[0]} px (white), " if info["refinement_requested"] else "no gate, "
        zoom.set_title(f"{gate}template {info['base_size']} px (cyan)", fontsize=10)
        zoom.axis("off")
    return ax


def main():
    parser = argparse.ArgumentParser(description=textwrap.dedent(__doc__ or ""), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("roi", help="ROI filename, e.g. 300.tiff")
    parser.add_argument("--source", required=True, choices=tuple(SOURCE_OPTIONS), help="where the seed comes from: a raw click, a drawn box, or a MIDOG++ annotation")
    parser.add_argument("--click", type=float, nargs=2, metavar=("X", "Y"), default=None, help="--source click: the click in full-resolution ROI pixels, pixel-centre coordinates unless --edge-origin; predictions only, no scoring")
    parser.add_argument("--box", type=float, nargs=4, metavar=("X", "Y", "W", "H"), default=None, help="--source box: the drawn box's left, top, width and height in pixel-edge coordinates (COCO / canvas); 18 to 73 px on its longer side; predictions only, no scoring")
    parser.add_argument("--object-id", type=int, default=None, help="--source click or box: the seed row's object_id (default: the session's next free id, 0 without one)")
    parser.add_argument("--edge-origin", action="store_true", default=None, help="--source click: --click uses pixel-edge coordinates (COCO / canvas); 0.5 is subtracted")
    parser.add_argument("--ann-id", type=int, default=None, help="--source midog: pin the seed to this MIDOG mitotic annotation")
    parser.add_argument("--seed-index", type=int, default=None, help="--source midog without --ann-id: the validation draw's seed (default 0)")
    parser.add_argument("--rank-key", default="tm_score", choices=sorted(prod.AXES))
    parser.add_argument("--refinement-channel", default="gray_inverted", choices=("gray_inverted", "hematoxylin_od", "none"), help="the Otsu gate's channel; none = no tightening: the template is the human's box (51 px for clicks and MIDOG)")
    parser.add_argument("--no-refine-predictions", action="store_true", help="skip refining each prediction's box the way the seed's was (no PREDICTION_REFINED rows)")
    parser.add_argument("--session-csv", default=None, help="the session so far (a --save-session file): its stored boxes are blanked and this run is its next iteration; without it this is iteration 0")
    parser.add_argument("--save-session", default=None, help="path to write the session with this iteration appended (may be the --session-csv file)")
    parser.add_argument("--k", type=int, default=af.DEFAULT_K, help="predictions shown, and stored in the session, per iteration")
    parser.add_argument("--allow-overlap", action="store_true", help="seed even when the new box contains the centre of a box the session stored")
    parser.add_argument("--images-dir", default=str(DEFAULT_IMAGES_DIR))
    parser.add_argument("--db-path", default=str(DEFAULT_DB_PATH))
    parser.add_argument("--show", type=int, default=5, help="how many predictions to print after the seed rows")
    parser.add_argument("--save-csv", default=None, help="path to write the full output frame to (read back with float_precision='round_trip')")
    parser.add_argument("--top-n", type=int, default=10, help="how many top-ranked predictions to draw")
    parser.add_argument("--save-fig", default=None, help="path to save the overlay figure to")
    parser.add_argument("--plot", action="store_true", help="show the overlay figure")
    args = parser.parse_args()

    seed_args = dict(click=None if args.click is None else tuple(args.click), box=None if args.box is None else tuple(args.box), object_id=args.object_id, edge_origin=args.edge_origin, ann_id=args.ann_id, seed_index=args.seed_index)
    try:
        check_seed_args(args.source, **seed_args)
    except ValueError as exc:
        parser.error(str(exc))
    if args.k < 1:
        parser.error(f"--k must be at least 1, got {args.k}")
    if args.save_fig and not args.plot:
        import matplotlib
        matplotlib.use("Agg")  # save without opening a window

    channel = None if args.refinement_channel == "none" else args.refinement_channel
    session = read_session(args.session_csv) if args.session_csv else None
    result = run_frame_on_roi(args.roi, args.source, **seed_args, rank_key=args.rank_key, refinement_channel=channel, refine_predictions=not args.no_refine_predictions, session=session, k=args.k, check_overlap=not args.allow_overlap, images_dir=args.images_dir, db_path=args.db_path, keep_rgb=bool(args.save_fig or args.plot))

    info = result["info"]
    print(f"{args.roi}: iteration {info['iteration']}, {args.source} seed {info['seed_object_id']} at {info['click_center']}, template {info['base_size']} px at {info['template_center']}, offset {info['offset_px']:.2f}px, refined={info['refined']}, refusal_reason={info['refusal_reason']}")
    if info["iteration"] > 0:
        print(f"blanked {info['n_extra_blanks']} boxes stored by iterations 0..{info['iteration'] - 1} ({info['n_extra_blanked_px']} px); {info['n_repeats_in_page']} repeats skipped from the first {info['k']} (all repeats: {info['repeat_object_ids'] or 'none'})")
    print(f"page shown: {info['n_shown']} predictions, object_ids {info['shown_object_ids']}")
    print(pd.Series(summarize(result)).to_frame().T.to_string(index=False))
    if result["metrics"] is None:
        why = "multi-round scoring is deferred: evaluate would read the raw list, not the page, against ground truth that still holds earlier finds" if info["iteration"] > 0 else "its object_id is not a MIDOG ann_id"
        print(f"{args.source} seed: not scored ({why})")

    frame = result["frame"]
    n_seed_rows = int((frame["class"] == af.MITOTIC_GT).sum())
    print(frame.head(n_seed_rows + args.show).drop(columns="current_time").to_string())
    if "prediction_refinement" in info:
        table = info["prediction_refinement"]
        reasons = table["refusal_reason"].value_counts().to_dict()
        gate = f"the gate refined {info['n_predictions_refined']} of {len(table)} boxes; fallbacks {reasons or 'none'}" if info["refinement_requested"] else f"no gate (--refinement-channel none), all {len(table)} boxes are the predictions' own"
        print(f"refined predictions: {gate}")
        print(frame[frame["class"] == af.PREDICTION_REFINED].head(args.show).drop(columns="current_time").to_string())
    if args.save_csv:
        frame.to_csv(args.save_csv, index=False)
        print(f"saved {len(frame)} rows to {args.save_csv}")
    if args.save_session:
        result["session"].to_csv(args.save_session, index=False)
        print(f"saved the session ({info['iteration'] + 1} iterations, {len(result['session'])} rows) to {args.save_session}")

    if args.save_fig or args.plot:
        ax = visualize(result, top_n=args.top_n)
        if args.save_fig:
            ax.figure.savefig(args.save_fig, dpi=150, bbox_inches="tight")
            print(f"saved figure to {args.save_fig}")
        if args.plot:
            import matplotlib.pyplot as plt
            plt.show()


if __name__ == "__main__":
    main()
