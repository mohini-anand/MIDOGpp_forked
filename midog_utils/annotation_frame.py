"""
    The annotation frame: one DataFrame carrying every box of a click-to-verify session -- the
    human's click, the template box the algorithm refined it to, every prediction, and, on request,
    each prediction refined the same way. `run_roi` takes a frame holding the click and returns a
    new frame with the refined seed and the predictions appended; `evaluate` scores that output
    against MIDOG++ ground truth;
    `click_frame` and `box_frame` build the input row for a raw click or a drawn box, and
    `draw_seed_rows` draws the human-seed rows the validation path feeds `run_roi`. The seed row's
    box sets the Otsu gate's window and the fallback template (`gate_window_side`): 51 px for every
    click and MIDOG row, the drawn box's own size for a box. `run_iteration` chains runs into a
    session: each iteration's seed, refined seed and first k shown predictions are stored with an
    `iteration` number, and every later search blanks those stored boxes. A layer on top of the
    frozen modules, which do every numerical step. Spec: `ANNOTATION_DATAFRAME_SPEC.md`.
"""

from __future__ import annotations

import math
import time
from datetime import datetime, timezone
from types import SimpleNamespace

import numpy as np
import pandas as pd

from . import channels as ch
from . import dataset as ds
from . import evaluate as ev
from . import production as prod
from . import seed_selection as ss
from . import template_match as tm

COLUMNS = ("file_name", "x_tl", "y_tl", "width", "height", "class", "object_id", "made_by", "refined", "current_time", "score", "od")
DTYPES = {"file_name": object, "x_tl": np.float64, "y_tl": np.float64, "width": np.float64, "height": np.float64, "class": object, "object_id": np.int64, "made_by": object, "refined": bool, "current_time": object, "score": np.float64, "od": np.float64}
COORD_COLUMNS = ("x_tl", "y_tl", "width", "height")

MITOTIC_GT = "MITOTIC_GT"
PREDICTION = "PREDICTION"
PREDICTION_REFINED = "PREDICTION_REFINED"  # a prediction's box after the seed's template choice; object_id = that prediction's rank
HUMAN = "HUMAN"
ALGORITHM = "ALGORITHM"

GATE_REFUSED = "gate_refused"                # the Otsu gate refused the click
TEMPLATE_UNREADABLE = "template_unreadable"  # the gate accepted, but the refined template is within 36 px of an edge
NEAR_EDGE = "near_edge"                      # a prediction centre within 36 px of an edge, not gated (a seed there raises instead)

DEFAULT_BASE_SIZE = tm.BASE_SIZE   # 51: the untightened template, and the fallback, for clicks and MIDOG rows (gate_window_side(50, 50))
OTSU_WINDOW = tm.BASE_SIZE         # 51: the validation draw's gate window (_joint_spec)
PREDICTION_BOX = float(ds.BOX_SIZE)  # 50.0: a prediction row's width and height
CLICK_BOX = float(ds.BOX_SIZE)       # 50.0: a raw click's row width and height, MIDOG's own box size; fixes a click's 51 px window and fallback
MAX_WINDOW_PX = tm.PATCH_SIZE      # 73: the largest gate window; a larger human box raises
MIN_WINDOW_PX = 19                 # the smallest gate window, i.e. a drawn box of 18 px or more (DECISIONS.md D12)
REFINEMENT_CHANNELS = (None, "gray_inverted", "hematoxylin_od")

ITERATION = "iteration"                        # the session column: 0 for the first iteration, then 1, 2, ...
SESSION_COLUMNS = COLUMNS + (ITERATION,)
SESSION_KEY = [ITERATION, "class", "object_id", "made_by", "refined"]  # a session's row key (one file_name per session)
DEFAULT_K = 10                                 # predictions shown, and stored, per iteration (user, 2026-10-06)
JOINT_CONDITIONS = {"default_51": None, "gray_bbox": "gray_inverted", "hem_bbox": "hematoxylin_od"}  # the driver's names and order


def now_utc() -> str:
    """
        The frame's timestamp: ISO-8601 UTC with microseconds.

        Returns str: the current time.
    """
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


def box_center(x_tl, width):
    """
        A box's continuous centre on one axis, ``x_tl + width / 2``, never rounded.

        x_tl (float or pd.Series): top-left coordinate, pixel-centre origin.
        width (float or pd.Series): extent, ``x2 - x1``.

        Returns float or pd.Series: the centre.
    """
    return x_tl + width / 2


def footprint_box(base_size: int, tx: float, ty: float) -> tuple[float, float, float, float]:
    """
        The frame box of the ``base_size`` pixel footprint centred on ``(round(tx), round(ty))``. Its width is ``base_size - 1``
        (so ``base_size = width + 1``): the box's edges are the centres of the footprint's outermost pixels. Rounds half to
        even, as ``read_padded_patch`` does.

        base_size (int): the footprint's side in pixels; must be odd.
        tx (float): footprint centre x, unrounded.
        ty (float): footprint centre y, unrounded.

        Returns tuple[float, float, float, float]: (x_tl, y_tl, width, height), with width = height = base_size - 1.
    """
    if base_size % 2 != 1:
        raise ValueError(f"base_size must be odd for width = base_size - 1 to hold, got {base_size}")
    half = base_size // 2
    return float(int(round(tx)) - half), float(int(round(ty)) - half), float(base_size - 1), float(base_size - 1)


def gate_window_side(width, height) -> int:
    """
        The Otsu gate's window, and the fallback template's size, for a human row's box: the pixel count of its longer side
        (extent + 1), raised to odd. 51 for every click and MIDOG row (width 50).

        width (float): the box's extent, ``x2 - x1``, so it spans width + 1 pixels.
        height (float): the box's extent on y.

        Returns int: the odd window side, MIN_WINDOW_PX..MAX_WINDOW_PX.
        Raises ValueError when width or height is not finite or negative, or the window falls outside those limits.
    """
    width, height = float(width), float(height)
    if not (math.isfinite(width) and math.isfinite(height) and width >= 0 and height >= 0):
        raise ValueError(f"a human box's width and height must be finite and >= 0, got ({width}, {height})")
    n_px = int(round(max(width, height))) + 1  # round half to even
    side = n_px if n_px % 2 else n_px + 1  # odd: read_padded_patch and blank_seed_square act on 2 * (n // 2) + 1 pixels
    if side > MAX_WINDOW_PX:
        raise ValueError(f"a human box of {n_px} px needs a {side} px gate window, above the {MAX_WINDOW_PX} px limit (tm.PATCH_SIZE)")
    if side < MIN_WINDOW_PX:
        raise ValueError(f"a human box of {n_px} px gives a {side} px gate window, below the {MIN_WINDOW_PX} px limit; draw a box of {MIN_WINDOW_PX - 1} px or more")
    assert side % 2 == 1 and side <= tm.PATCH_SIZE
    return side


def _snap_axis(x0: float, w: float) -> tuple[float, float]:
    """
        One axis of a drawn box snapped to whole pixels: ``round(w)`` pixels raised to odd, centred on the pixel containing the
        drawn centre. ``floor(x0 + w / 2)`` is round half up in the frame's centre coordinates, a deliberate carve-out from the
        spec's ban (§1): it picks the pixel containing an edge-coordinate point, has no tie, and does not depend on the box's
        position, and the snapped centre is a whole pixel, so the frozen modules' half-to-even rounding never meets a tie.

        x0 (float): the box's left (or top) edge, pixel-edge coordinates.
        w (float): the box's width (or height) in pixels.

        Returns tuple[float, float]: (x_tl, width) in the frame's pixel-centre coordinates, width = the odd pixel count - 1.
        Raises ValueError when x0 or w is not finite, or w rounds below one pixel.
    """
    x0, w = float(x0), float(w)
    if not (math.isfinite(x0) and math.isfinite(w)):
        raise ValueError(f"a drawn box's origin and size must be finite, got ({x0}, {w})")
    n = int(round(w))
    if n < 1:
        raise ValueError(f"a drawn box must be at least one pixel wide, got {w}")
    n_odd = n if n % 2 else n + 1
    c = math.floor(x0 + w / 2)
    return float(c - n_odd // 2), float(n_odd - 1)


def _block(n: int, **columns) -> pd.DataFrame:
    """
        A frame block of ``n`` rows in the schema's column order and dtypes.

        n (int): the number of rows.
        columns (object or sequence): one value, repeated n times, or one n-long sequence per schema column.

        Returns pd.DataFrame: the block.
    """
    data = {c: pd.Series(columns[c] if np.ndim(columns[c]) else [columns[c]] * n, dtype=DTYPES[c]) for c in COLUMNS}
    return pd.DataFrame(data)


def human_seed_frame(ann: pd.Series, file_name: str, current_time: str = None) -> pd.DataFrame:
    """
        The one-row human-seed frame for a MIDOG++ annotation. ``x_tl = cx - w / 2`` recovers MIDOG's ``x1`` bit-exactly.

        ann (pd.Series): one ``dataset.load_annotations`` row (ann_id, cx, cy, w, h).
        file_name (str): the ROI's file name.
        current_time (str or None): the click time; now when None.

        Returns pd.DataFrame: one MITOTIC_GT, HUMAN, unrefined row whose object_id is the ann_id.
    """
    cx, cy, w, h = float(ann["cx"]), float(ann["cy"]), float(ann["w"]), float(ann["h"])
    return _block(1, file_name=file_name, x_tl=cx - w / 2, y_tl=cy - h / 2, width=w, height=h, **{"class": MITOTIC_GT}, object_id=int(ann["ann_id"]), made_by=HUMAN, refined=False, current_time=current_time or now_utc(), score=np.nan, od=np.nan)


def click_frame(file_name: str, cx: float, cy: float, object_id: int, *, edge_origin: bool = False, current_time: str = None) -> pd.DataFrame:
    """
        The one-row human-seed frame for a raw click: a CLICK_BOX-wide box centred on it, so ``x_tl + width / 2`` gives the
        click back, and the gate window and fallback template are 51 px, as for a MIDOG row. The production counterpart of
        ``human_seed_frame``, which needs a MIDOG annotation row.

        file_name (str): the ROI's file name.
        cx (float): click x in full-resolution ROI pixels.
        cy (float): click y in full-resolution ROI pixels.
        object_id (int): the click's id, unique among the ROI's HUMAN MITOTIC_GT rows.
        edge_origin (bool): True when the click uses pixel-edge coordinates (pixel k spans [k, k + 1), as in COCO and canvas APIs); 0.5 is subtracted to reach the frame's pixel-centre origin.
        current_time (str or None): the click time; now when None.

        Returns pd.DataFrame: one MITOTIC_GT, HUMAN, unrefined row.
    """
    cx, cy = float(cx), float(cy)
    if edge_origin:
        cx, cy = cx - 0.5, cy - 0.5
    half = CLICK_BOX / 2
    return _block(1, file_name=file_name, x_tl=cx - half, y_tl=cy - half, width=CLICK_BOX, height=CLICK_BOX, **{"class": MITOTIC_GT}, object_id=int(object_id), made_by=HUMAN, refined=False, current_time=current_time or now_utc(), score=np.nan, od=np.nan)


def box_frame(file_name: str, x0: float, y0: float, w: float, h: float, object_id: int, *, current_time: str = None) -> pd.DataFrame:
    """
        The one-row human-seed frame for a box the human drew. Each side is snapped to whole pixels -- ``round`` of its size,
        raised to odd, centred on the pixel containing the drawn centre -- so the stored row is exactly the gate window (for a
        square box) and the fallback template. Sub-pixel coordinates, and an even size's difference from the next odd one,
        are not kept. For drawn boxes only: a MIDOG box is already in pixel-centre coordinates (``human_seed_frame``).

        file_name (str): the ROI's file name.
        x0 (float): the box's left edge, pixel-edge coordinates (pixel k spans [k, k + 1), as in COCO and canvas APIs).
        y0 (float): the box's top edge, pixel-edge coordinates.
        w (float): the box's width in pixels.
        h (float): the box's height in pixels.
        object_id (int): the box's id, unique among the ROI's HUMAN MITOTIC_GT rows.
        current_time (str or None): the drawing time; now when None.

        Returns pd.DataFrame: one MITOTIC_GT, HUMAN, unrefined row.
        Raises ValueError on a non-finite or empty box, or one whose gate window falls outside 19..73 px (``gate_window_side``).
    """
    x_tl, width = _snap_axis(x0, w)
    y_tl, height = _snap_axis(y0, h)
    try:
        gate_window_side(width, height)
    except ValueError as exc:
        raise ValueError(f"drawn box {w} x {h} px snaps to {int(width) + 1} x {int(height) + 1} px: {exc}") from None
    return _block(1, file_name=file_name, x_tl=x_tl, y_tl=y_tl, width=width, height=height, **{"class": MITOTIC_GT}, object_id=int(object_id), made_by=HUMAN, refined=False, current_time=current_time or now_utc(), score=np.nan, od=np.nan)


def _select_seed(frame: pd.DataFrame, seed_object_id) -> pd.Series:
    """
        The one HUMAN MITOTIC_GT row, with ``seed_object_id`` when given.

        frame (pd.DataFrame): an annotation frame.
        seed_object_id (int or None): the seed's object_id, or None when only one such row exists.

        Returns pd.Series: the seed row. Raises ValueError on zero or several matches.
    """
    rows = frame[(frame["class"] == MITOTIC_GT) & (frame["made_by"] == HUMAN)]
    if seed_object_id is not None:
        rows = rows[rows["object_id"] == seed_object_id]
    if len(rows) != 1:
        which = "" if seed_object_id is None else f" with object_id {seed_object_id}"
        raise ValueError(f"expected exactly one HUMAN MITOTIC_GT seed row{which}, found {len(rows)}")
    return rows.iloc[0]


def _input_rows(annotations_df: pd.DataFrame) -> pd.DataFrame:
    """
        A validated copy of ``run_roi``'s input rows, coordinate columns cast to float64 (lossless) and nothing else changed.

        annotations_df (pd.DataFrame): the caller's frame; never mutated.

        Returns pd.DataFrame: the copy, in the schema's column order.
        Raises ValueError on missing or extra columns, a row that is not HUMAN MITOTIC_GT, a duplicated object_id, or a row
        whose box gives no valid gate window (``gate_window_side``).
    """
    missing = [c for c in COLUMNS if c not in annotations_df.columns]
    extra = [c for c in annotations_df.columns if c not in COLUMNS]
    if missing or extra:
        raise ValueError(f"input frame columns must be exactly {list(COLUMNS)}; missing {missing}, extra {extra}")
    frame = annotations_df.loc[:, list(COLUMNS)].copy().reset_index(drop=True)
    frame = frame.astype({c: np.float64 for c in COORD_COLUMNS})
    not_seed = ~((frame["class"] == MITOTIC_GT) & (frame["made_by"] == HUMAN))
    if not_seed.any():
        raise ValueError(f"{int(not_seed.sum())} input rows are not HUMAN MITOTIC_GT; run_roi takes only the human's clicks, never a previous run's output")
    if frame["object_id"].duplicated().any():
        raise ValueError(f"input object_ids must be unique; duplicated: {sorted(frame.loc[frame['object_id'].duplicated(), 'object_id'].unique().tolist())}")
    for object_id, width, height in zip(frame["object_id"], frame["width"], frame["height"]):
        try:
            gate_window_side(width, height)
        except ValueError as exc:
            raise ValueError(f"input row object_id {int(object_id)}: {exc}") from None
    return frame


def _check_schema(frame: pd.DataFrame) -> None:
    """
        Raise unless a frame has the schema's columns, in order, with the schema's dtypes.

        frame (pd.DataFrame): the frame to check.

        Returns None. Raises ValueError naming the first mismatch.
    """
    if tuple(frame.columns) != COLUMNS:
        raise ValueError(f"frame columns {list(frame.columns)} != {list(COLUMNS)}")
    for c in COLUMNS:
        if frame[c].dtype != np.dtype(DTYPES[c]):
            raise ValueError(f"frame column {c!r} has dtype {frame[c].dtype}, expected {np.dtype(DTYPES[c])}; an input column breaks the schema")


def _template_choice(rgb: np.ndarray, cx: float, cy: float, refinement_channel, window: int, gate_image: np.ndarray = None):
    """
        The template spec for a readable seed: the accepted Otsu component's, or the human's box (``window`` px) centred on the
        seed.

        rgb (np.ndarray): the ROI's RGB image.
        cx (float): seed centre x, unrounded.
        cy (float): seed centre y, unrounded.
        refinement_channel (str or None): the gate's channel, or None for no refinement.
        window (int): the gate window and fallback size, ``gate_window_side`` of the seed row; 51 for clicks and MIDOG rows.
        gate_image (np.ndarray or None): the whole ROI in refinement_channel, when the caller already has it; converted here when None.

        Returns tuple: (base_size, tx, ty, refined, refusal_reason).
        Raises RuntimeError when the window is not readable at the seed centre, so a refusal always means the gate refused.
    """
    if refinement_channel is None:
        return window, cx, cy, False, None
    if not ss._patch_readable(rgb.shape, cx, cy, window):
        raise RuntimeError(f"the {window} px gate window at ({cx}, {cy}) is not readable; run_roi's edge check should have raised")
    if gate_image is None:
        gate_image = ch.to_channel(rgb, refinement_channel)
    got = ss.tightened_template_box(gate_image, cx, cy, otsu_window=window)
    if got is None:
        return window, cx, cy, False, GATE_REFUSED
    base_size, tx, ty = int(got[0]), float(got[1]), float(got[2])
    if not ss._patch_readable(rgb.shape, tx, ty, tm.PATCH_SIZE):
        return window, cx, cy, False, TEMPLATE_UNREADABLE
    return base_size, tx, ty, True, None


def _prediction_rows(detections: pd.DataFrame, file_name: str, current_time: str) -> pd.DataFrame:
    """
        One PREDICTION row per production detection, in production's order, with object_id = rank.

        detections (pd.DataFrame): ``run_production_pipeline``'s ranked detections, with od.
        file_name (str): the ROI's file name.
        current_time (str): the prediction stage's stamp.

        Returns pd.DataFrame: the prediction block.
    """
    n = len(detections)
    rank = detections["rank"].to_numpy(dtype=np.int64)
    if not np.array_equal(rank, np.arange(n)):
        raise RuntimeError("production's rank is not its row position; object_id = rank would not be the prediction's position")
    cx, cy = detections["cx"].to_numpy(dtype=np.float64), detections["cy"].to_numpy(dtype=np.float64)
    return _block(n, file_name=file_name, x_tl=cx - PREDICTION_BOX / 2, y_tl=cy - PREDICTION_BOX / 2, width=np.full(n, PREDICTION_BOX), height=np.full(n, PREDICTION_BOX), **{"class": PREDICTION}, object_id=rank, made_by=ALGORITHM, refined=False, current_time=current_time, score=detections["score"].to_numpy(dtype=np.float64), od=detections["od"].to_numpy(dtype=np.float64))


def _refined_prediction_rows(rgb: np.ndarray, gate_image, detections: pd.DataFrame, refinement_channel, file_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
        The seed's template choice (``_template_choice``) applied at every prediction's centre, with the prediction's own box as
        the window: 51 px (user, 2026-10-06). One PREDICTION_REFINED row per detection, in production's order, object_id = rank;
        on a refusal the row is the fallback, the prediction's own box, with refined False. A prediction centre within 36 px of
        an edge is not gated (NEAR_EDGE), the predicate that makes a seed raise.

        rgb (np.ndarray): the ROI's RGB image.
        gate_image (np.ndarray or None): the whole ROI in refinement_channel; None when refinement_channel is None.
        detections (pd.DataFrame): ``run_production_pipeline``'s ranked detections.
        refinement_channel (str or None): the seed's gate channel, or None for no refinement.
        file_name (str): the ROI's file name.

        Returns tuple[pd.DataFrame, pd.DataFrame]: the PREDICTION_REFINED block, and one row per prediction with object_id,
        base_size, the unrounded template centre (template_cx, template_cy), offset_px from the prediction centre, refined and
        refusal_reason.
    """
    window = gate_window_side(PREDICTION_BOX, PREDICTION_BOX)  # 51: the prediction row's own box
    cx, cy = detections["cx"].to_numpy(dtype=np.float64), detections["cy"].to_numpy(dtype=np.float64)
    choices = []
    for x, y in zip(cx, cy):
        if gate_image is not None and not ss._patch_readable(rgb.shape, x, y, tm.PATCH_SIZE):
            choices.append((window, x, y, False, NEAR_EDGE))
        else:
            choices.append(_template_choice(rgb, x, y, refinement_channel, window, gate_image=gate_image))
    table = pd.DataFrame(choices, columns=["base_size", "template_cx", "template_cy", "refined", "refusal_reason"]).astype({"base_size": np.int64, "template_cx": np.float64, "template_cy": np.float64, "refined": bool})
    table.insert(0, "object_id", detections["rank"].to_numpy(dtype=np.int64))
    table.insert(4, "offset_px", np.hypot(table["template_cx"].to_numpy() - cx, table["template_cy"].to_numpy() - cy))
    boxes = np.array([footprint_box(b, tx, ty) for b, tx, ty in zip(table["base_size"], table["template_cx"], table["template_cy"])], dtype=np.float64).reshape(-1, 4)
    rows = _block(len(table), file_name=file_name, x_tl=boxes[:, 0], y_tl=boxes[:, 1], width=boxes[:, 2], height=boxes[:, 3], **{"class": PREDICTION_REFINED}, object_id=table["object_id"].to_numpy(), made_by=ALGORITHM, refined=table["refined"].to_numpy(), current_time=now_utc(), score=np.nan, od=np.nan)
    return rows, table


def run_roi(rgb: np.ndarray, annotations_df: pd.DataFrame, mpp: float, *, seed_object_id=None, rank_key: str = "tm_score", refinement_channel="gray_inverted", refine_predictions: bool = False, extra_blanks=None) -> tuple[pd.DataFrame, dict]:
    """
        One human click or box in, the refined seed and every post-NMS prediction appended -- and, with refine_predictions,
        every prediction refined the same way. Never draws a seed and never re-sorts. The seed row's box sets the gate window
        and the fallback template (``gate_window_side``). Only the search sees extra_blanks: the gates, the template and od
        read the unblanked image.

        rgb (np.ndarray): the ROI's full-resolution RGB image (`dataset.load_roi`).
        annotations_df (pd.DataFrame): the frame's input rows, all HUMAN MITOTIC_GT of one ROI; never mutated.
        mpp (float): the ROI's microns per pixel (`dataset.roi_mpp`).
        seed_object_id (int or None): which input row is the seed; None when there is exactly one.
        rank_key (str): "tm_score" or "chromatin_od" -- the sort order only; score and od are both always computed.
        refinement_channel (str or None): the Otsu gate's channel, "gray_inverted" or "hematoxylin_od", or None for no tightening: the template is the human's box (51 px for clicks and MIDOG).
        refine_predictions (bool): also apply the seed's template choice at every prediction (51 px window) and append one PREDICTION_REFINED row each; False leaves the frame and info exactly as before.
        extra_blanks (list[tuple[float, float, int]] or None): (cx, cy, side) squares blanked in production's search copy beside the seed's own (`find_and_suppress.blank_squares`), e.g. the boxes earlier session iterations stored (`run_iteration`); None leaves the frame and info exactly as before.

        Returns tuple[pd.DataFrame, dict]: the input rows, the refined seed, the predictions and (refine_predictions) the refined predictions, and the run's info.
        Raises ValueError on a bad argument or input frame (including a human box outside the 19..73 px window limits), or a
        seed centre within 36 px of an edge.
    """
    t_start = time.time()
    if rank_key not in prod.AXES:
        raise ValueError(f"rank_key must be one of {sorted(prod.AXES)}, got {rank_key!r}")
    if refinement_channel not in REFINEMENT_CHANNELS:
        raise ValueError(f"refinement_channel must be one of {list(REFINEMENT_CHANNELS)}, got {refinement_channel!r}")
    frame = _input_rows(annotations_df)
    seed = _select_seed(frame, seed_object_id)
    file_name = seed["file_name"]
    if not (frame["file_name"] == file_name).all():
        raise ValueError(f"every input row must share the seed's file_name {file_name!r}; found {sorted(frame['file_name'].unique().tolist())}")
    window = gate_window_side(seed["width"], seed["height"])

    cx, cy = box_center(float(seed["x_tl"]), float(seed["width"])), box_center(float(seed["y_tl"]), float(seed["height"]))
    if not ss._patch_readable(rgb.shape, cx, cy, tm.PATCH_SIZE):
        raise ValueError(f"seed centre ({cx}, {cy}) is within {tm.PATCH_SIZE // 2} px of the ROI edge; seeds that close to the edge are not supported")

    t0 = time.time()
    gate_image = None if refinement_channel is None else ch.to_channel(rgb, refinement_channel)
    base_size, tx, ty, refined, refusal_reason = _template_choice(rgb, cx, cy, refinement_channel, window, gate_image=gate_image)
    t_refine = time.time() - t0
    refine_time = now_utc()
    if not refine_predictions:
        gate_image = None  # free it before the search, as before; only the prediction refinement reads it again

    detections, prod_info = prod.run_production_pipeline(rgb, SimpleNamespace(base_size=base_size, template_xy=(tx, ty)), mpp, rank_key=rank_key, always_od=(rank_key != "chromatin_od"), extra_blanks=extra_blanks)
    predict_time = now_utc()

    x_tl, y_tl, width, height = footprint_box(base_size, tx, ty)
    refined_row = _block(1, file_name=file_name, x_tl=x_tl, y_tl=y_tl, width=width, height=height, **{"class": MITOTIC_GT}, object_id=int(seed["object_id"]), made_by=ALGORITHM, refined=refined, current_time=refine_time, score=np.nan, od=np.nan)
    blocks = [frame, refined_row, _prediction_rows(detections, file_name, predict_time)]
    if refine_predictions:
        t0 = time.time()
        refined_preds, pred_table = _refined_prediction_rows(rgb, gate_image, detections, refinement_channel, file_name)
        t_refine_preds = time.time() - t0
        blocks.append(refined_preds)
    out = pd.concat(blocks, ignore_index=True)
    _check_schema(out)

    added = dict(file_name=file_name, mpp=mpp, match_radius_px=ev.radius_px(mpp), refinement_channel=refinement_channel, refinement_requested=refinement_channel is not None, channel=prod.CHANNEL, tm_method=prod.TM_METHOD, peak_min_distance=prod.PEAK_MIN_DISTANCE, deep_floor_z=prod.DEEP_FLOOR_Z, max_peaks=prod.MAX_PEAKS, patch_size=tm.PATCH_SIZE, base_size=base_size, seed_object_id=int(seed["object_id"]), click_center=(cx, cy), refined=refined, refusal_reason=refusal_reason, template_center=(tx, ty), offset_px=float(np.hypot(tx - cx, ty - cy)), t_refine_s=round(t_refine, 2))
    if refine_predictions:
        added.update(prediction_refinement=pred_table, n_predictions_refined=int(pred_table["refined"].sum()), t_refine_predictions_s=round(t_refine_preds, 2))
    clash = sorted(set(added) & set(prod_info))
    if clash:
        raise RuntimeError(f"run_roi's info keys would overwrite production's: {clash}")
    info = dict(prod_info)
    info.update(added)
    info["t_total_s"] = round(time.time() - t_start, 2)
    return out, info


def evaluate(result_df: pd.DataFrame, gt_df: pd.DataFrame, mpp: float, roi_shape, *, seed_object_id=None, budgets=(10, 20, 30)) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
        Score a ``run_roi`` frame against MIDOG++ ground truth, with the seed's own annotation excluded. MIDOG only: the seed
        row's object_id must be its ann_id.

        result_df (pd.DataFrame): ``run_roi``'s output frame.
        gt_df (pd.DataFrame): the ROI's full ground truth, mitotic and look-alike, in ``dataset.load_annotations`` format.
        mpp (float): the ROI's microns per pixel.
        roi_shape (tuple): the ROI's (h, w) or (h, w, c).
        seed_object_id (int or None): chooses the seed row as ``run_roi`` does; never used to filter ground truth.
        budgets (tuple[int]): the K values for tp_at_K.

        Returns tuple[pd.DataFrame, pd.DataFrame, dict]: the bucketed predictions sorted by object_id, ``evaluate_run``'s gt_out, and the metrics plus tp_at_K and seed_object_id.
        Raises ValueError when the seed row's object_id is not one mitotic ann_id at its centre, or the prediction ids are not 0..n-1.
    """
    seed = _select_seed(result_df, seed_object_id)
    seed_id = int(seed["object_id"])
    cx, cy = box_center(float(seed["x_tl"]), float(seed["width"])), box_center(float(seed["y_tl"]), float(seed["height"]))
    own = gt_df[gt_df["ann_id"] == seed_id]
    if len(own) != 1 or int(own["category_id"].iat[0]) != ds.MITOTIC or float(own["cx"].iat[0]) != cx or float(own["cy"].iat[0]) != cy:
        raise ValueError(f"seed object_id {seed_id} is not one mitotic ann_id at the seed's centre ({cx}, {cy}); evaluate needs the seed row's object_id to be its MIDOG ann_id")
    gt_eval = gt_df[gt_df["ann_id"] != seed_id].reset_index(drop=True)

    cands = result_df[result_df["class"] == PREDICTION].sort_values("object_id", kind="stable").reset_index(drop=True)
    if not np.array_equal(cands["object_id"].to_numpy(), np.arange(len(cands))):
        raise ValueError("prediction object_ids are not the contiguous 0..n-1 run_roi emits; gt_out's matched_rank would not be an object_id")
    detections = pd.DataFrame({"rank": cands["object_id"].to_numpy(), "cx": box_center(cands["x_tl"], cands["width"]).to_numpy(), "cy": box_center(cands["y_tl"], cands["height"]).to_numpy(), "score": cands["score"].to_numpy(), "od": cands["od"].to_numpy()})
    h, w = roi_shape[0], roi_shape[1]
    area_mm2 = h * w * mpp * mpp / 1e6  # midog_utils_full.dataset.roi_area_mm2's exact expression; mpp ** 2 differs in the last bit
    det_out, gt_out, metrics = ev.evaluate_run(detections, gt_eval, ev.radius_px(mpp), area_mm2, roi_shape)

    hit = (det_out["bucket"] == ev.HUMAN_CORRECT_LABEL).to_numpy()
    for k in budgets:
        metrics[f"tp_at_{k}"] = int(hit[:k].sum())
    metrics["seed_object_id"] = seed_id
    scored = cands.assign(bucket=det_out["bucket"].to_numpy(), matched_ann_id=det_out["matched_ann_id"].to_numpy(), matched_category=det_out["matched_category"].to_numpy())
    return scored, gt_out, metrics


SESSION_KINDS = {(MITOTIC_GT, HUMAN), (MITOTIC_GT, ALGORITHM), (PREDICTION, ALGORITHM), (PREDICTION_REFINED, ALGORITHM)}  # the (class, made_by) pairs a session stores


def empty_session() -> pd.DataFrame:
    """
        A session with no iterations yet, in SESSION_COLUMNS order and dtypes.

        Returns pd.DataFrame: zero rows.
    """
    return _block(0, **{c: [] for c in COLUMNS}).assign(**{ITERATION: pd.Series([], dtype=np.int64)})


def check_session(session: pd.DataFrame) -> pd.DataFrame:
    """
        A validated copy of a session -- every completed iteration's stored rows -- in SESSION_COLUMNS order, coordinates cast
        to float64 and iteration to int64 (both lossless) and nothing else changed.

        session (pd.DataFrame): the session; never mutated.

        Returns pd.DataFrame: the copy.
        Raises ValueError on missing or extra columns, a dtype outside the schema, more than one file_name, iterations that are
        not 0..T-1, a (class, made_by) pair a session never stores, a duplicated row key (SESSION_KEY), a human object_id used
        twice, an iteration without exactly one human seed and its refined seed, refined predictions that are not one per
        stored prediction, or a stored algorithm box that is not a square of odd pixel side.
    """
    missing = [c for c in SESSION_COLUMNS if c not in session.columns]
    extra = [c for c in session.columns if c not in SESSION_COLUMNS]
    if missing or extra:
        raise ValueError(f"session columns must be exactly {list(SESSION_COLUMNS)}; missing {missing}, extra {extra}")
    if session.empty:
        return empty_session()
    frame = session.loc[:, list(SESSION_COLUMNS)].copy().reset_index(drop=True)
    frame = frame.astype({**{c: np.float64 for c in COORD_COLUMNS}, ITERATION: np.int64})
    for c in COLUMNS:
        if frame[c].dtype != np.dtype(DTYPES[c]):
            raise ValueError(f"session column {c!r} has dtype {frame[c].dtype}, expected {np.dtype(DTYPES[c])}")

    names = frame["file_name"].unique()
    if len(names) != 1:
        raise ValueError(f"a session holds one ROI; found file_names {sorted(names.tolist())}")
    iterations = np.sort(frame[ITERATION].unique())
    if not np.array_equal(iterations, np.arange(len(iterations))):
        raise ValueError(f"session iterations must be 0..T-1; found {iterations.tolist()}")
    kinds = set(zip(frame["class"], frame["made_by"])) - SESSION_KINDS
    if kinds:
        raise ValueError(f"a session never stores (class, made_by) {sorted(kinds)}")
    if frame.duplicated(SESSION_KEY).any():
        raise ValueError(f"{int(frame.duplicated(SESSION_KEY).sum())} session rows repeat a row key {SESSION_KEY}")
    seeds = frame["class"] == MITOTIC_GT
    human = frame[seeds & (frame["made_by"] == HUMAN)]
    if human["object_id"].duplicated().any():
        raise ValueError(f"human object_ids must be unique across the session; repeated: {sorted(human.loc[human['object_id'].duplicated(), 'object_id'].unique().tolist())}")
    for t, it in frame.groupby(ITERATION, sort=True):
        h, a = it[(it["class"] == MITOTIC_GT) & (it["made_by"] == HUMAN)], it[(it["class"] == MITOTIC_GT) & (it["made_by"] == ALGORITHM)]
        if len(h) != 1 or len(a) != 1 or int(h["object_id"].iat[0]) != int(a["object_id"].iat[0]):
            raise ValueError(f"iteration {t} must hold one human seed and its refined seed (same object_id); found {len(h)} and {len(a)}")
        pred_ids = sorted(it.loc[it["class"] == PREDICTION, "object_id"].tolist())
        refined_ids = sorted(it.loc[it["class"] == PREDICTION_REFINED, "object_id"].tolist())
        if refined_ids and refined_ids != pred_ids:
            raise ValueError(f"iteration {t}: refined predictions {refined_ids} are not one per stored prediction {pred_ids}")
    stored = frame[frame["made_by"] == ALGORITHM]
    side = stored["width"] + 1
    square = (stored["width"] == stored["height"]) & (side == np.round(side)) & (side % 2 == 1)
    if not square.all():
        raise ValueError(f"{int((~square).sum())} stored algorithm boxes are not squares of odd pixel side (width + 1)")
    return frame


def stored_blanks(session: pd.DataFrame) -> pd.DataFrame:
    """
        The squares the next iteration blanks: each iteration's refined seed, and each stored prediction's refined box (its
        PREDICTION_REFINED row: the gate's box when it refined, the prediction's own 51 px box on fallback), or its PREDICTION
        row when that iteration stored no refined predictions. Human rows are never blanked: as D11 does for the current seed,
        the refined seed's footprint is. A square is its row's footprint: centre ``x_tl + width / 2``, side ``width + 1``.

        session (pd.DataFrame): a validated session (``check_session``).

        Returns pd.DataFrame: one row per square -- iteration, class, object_id, cx, cy, side -- in session order.
    """
    with_refined = set(session.loc[session["class"] == PREDICTION_REFINED, ITERATION].tolist())
    pick = ((session["class"] == MITOTIC_GT) & (session["made_by"] == ALGORITHM)) | (session["class"] == PREDICTION_REFINED) | ((session["class"] == PREDICTION) & ~session[ITERATION].isin(with_refined))
    rows = session[pick]
    return pd.DataFrame({ITERATION: rows[ITERATION].to_numpy(dtype=np.int64), "class": rows["class"].to_numpy(dtype=object), "object_id": rows["object_id"].to_numpy(dtype=np.int64), "cx": box_center(rows["x_tl"], rows["width"]).to_numpy(dtype=np.float64), "cy": box_center(rows["y_tl"], rows["height"]).to_numpy(dtype=np.float64), "side": (rows["width"] + 1).to_numpy(dtype=np.float64).astype(np.int64)})


def next_object_id(session) -> int:
    """
        The smallest object_id above every human seed's in a session: the default id for the next click or box.

        session (pd.DataFrame or None): the session; None or empty for the first iteration.

        Returns int: 0 for no session, else the largest human object_id + 1.
    """
    if session is None or len(session) == 0:
        return 0
    human = session[(session["class"] == MITOTIC_GT) & (session["made_by"] == HUMAN)]
    return int(human["object_id"].max()) + 1


def _check_seed_overlap(session: pd.DataFrame, seed: pd.Series) -> None:
    """
        Raise if the new seed's box contains the centre of any box the session stored: a human seed, a refined seed, a
        prediction or a refined prediction. Bounds are inclusive, ``x_tl <= cx <= x_tl + width`` and likewise on y; a click's
        box is its 50 px CLICK_BOX. Nearness alone is allowed: a box beside a stored one may hold another mitosis.

        session (pd.DataFrame): a validated, non-empty session.
        seed (pd.Series): the new seed row.

        Returns None. Raises ValueError naming the first stored boxes whose centres fall inside.
    """
    cx, cy = box_center(session["x_tl"], session["width"]), box_center(session["y_tl"], session["height"])
    x0, y0, w, h = float(seed["x_tl"]), float(seed["y_tl"]), float(seed["width"]), float(seed["height"])
    inside = (cx >= x0) & (cx <= x0 + w) & (cy >= y0) & (cy <= y0 + h)
    if inside.any():
        hits = session[inside].head(3)
        named = "; ".join(f"iteration {int(r[ITERATION])} {r['class']}/{r['made_by']} object_id {int(r['object_id'])} at ({cx[i]}, {cy[i]})" for i, r in hits.iterrows())
        raise ValueError(f"the new seed's box (x_tl {x0}, y_tl {y0}, {w} x {h}) contains the centres of {int(inside.sum())} stored boxes: {named}; pass check_overlap=False (run_frame.py --allow-overlap) to seed there anyway")


def prepare_iteration(session, seed_df: pd.DataFrame, *, check_overlap: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, int]:
    """
        The input wrapper: validate a session and the human's new seed, number the iteration, and list the squares the search
        must blank. The seed must be on the session's ROI, with an object_id no stored human row uses, and (check_overlap) a
        box containing no stored box's centre.

        session (pd.DataFrame or None): every completed iteration's stored rows (SESSION_COLUMNS); None or empty for the first iteration.
        seed_df (pd.DataFrame): the new seed, one HUMAN MITOTIC_GT row in COLUMNS (``click_frame`` or ``box_frame``).
        check_overlap (bool): raise when the seed's box contains a stored box's centre.

        Returns tuple: the validated session, the validated one-row seed frame, the squares to blank (``stored_blanks``; none for the first iteration), and the iteration number.
        Raises ValueError on an invalid session or seed, another ROI's seed, a used object_id, or (check_overlap) an overlap.
    """
    session = empty_session() if session is None else check_session(session)
    seed_frame = _input_rows(seed_df)
    if len(seed_frame) != 1:
        raise ValueError(f"the new seed must be exactly one row, got {len(seed_frame)}")
    seed = seed_frame.iloc[0]
    if session.empty:
        return session, seed_frame, stored_blanks(session), 0
    if seed["file_name"] != session["file_name"].iat[0]:
        raise ValueError(f"the new seed is on {seed['file_name']!r}, the session on {session['file_name'].iat[0]!r}")
    human = session[(session["class"] == MITOTIC_GT) & (session["made_by"] == HUMAN)]
    if int(seed["object_id"]) in set(human["object_id"].tolist()):
        raise ValueError(f"object_id {int(seed['object_id'])} is already a human seed in this session; the next free id is {next_object_id(session)}")
    if check_overlap:
        _check_seed_overlap(session, seed)
    return session, seed_frame, stored_blanks(session), int(session[ITERATION].max()) + 1


def merge_iteration(session: pd.DataFrame, frame: pd.DataFrame, info: dict, blanks: pd.DataFrame, iteration: int, *, k: int = DEFAULT_K) -> tuple[pd.DataFrame, dict]:
    """
        The output wrapper: append one ``run_roi`` frame to the session as iteration ``iteration``. It stores the human seed,
        its refined seed, and the first k predictions in rank order that are not repeats, with their refined rows; nothing
        beyond the page is stored. A repeat is a prediction whose centre is within the match radius
        (``info["match_radius_px"]``) of a blanked square's centre -- a leak beside a blanked object, or the uncovered part of
        a split figure -- and is skipped, never shown.

        session (pd.DataFrame): the validated session the run started from.
        frame (pd.DataFrame): ``run_roi``'s output for this iteration's seed.
        info (dict): ``run_roi``'s info for that run.
        blanks (pd.DataFrame): the squares the run blanked (``stored_blanks``).
        iteration (int): this iteration's number.
        k (int): how many predictions to show and store.

        Returns tuple[pd.DataFrame, dict]: the new session, and the info keys this step adds: iteration, k, n_shown,
        shown_object_ids, repeat_object_ids (every prediction flagged, any rank), n_repeats_in_page (repeats among the raw
        first k) and session_blanks (the blanked squares).
    """
    preds = frame[frame["class"] == PREDICTION].sort_values("object_id", kind="stable")
    px, py = box_center(preds["x_tl"], preds["width"]).to_numpy(), box_center(preds["y_tl"], preds["height"]).to_numpy()
    bx, by = blanks["cx"].to_numpy(dtype=np.float64), blanks["cy"].to_numpy(dtype=np.float64)
    repeat = np.zeros(len(preds), dtype=bool)
    if len(preds) and len(blanks):
        repeat = (np.hypot(px[:, None] - bx[None, :], py[:, None] - by[None, :]) <= info["match_radius_px"]).any(axis=1)
    shown = preds[~repeat].head(k)
    shown_ids = shown["object_id"].tolist()
    seed_rows = frame[frame["class"] == MITOTIC_GT]
    refined = frame[(frame["class"] == PREDICTION_REFINED) & frame["object_id"].isin(shown_ids)].sort_values("object_id", kind="stable")
    new = pd.concat([seed_rows, shown, refined], ignore_index=True).assign(**{ITERATION: np.int64(iteration)})
    merged = new if session.empty else pd.concat([session, new], ignore_index=True)
    added = dict(iteration=int(iteration), k=int(k), n_shown=len(shown_ids), shown_object_ids=[int(i) for i in shown_ids], repeat_object_ids=[int(i) for i in preds["object_id"].to_numpy()[repeat]], n_repeats_in_page=int(repeat[:k].sum()), session_blanks=blanks)
    return check_session(merged), added


def run_iteration(rgb: np.ndarray, session, mpp: float, seed_df: pd.DataFrame, *, k: int = DEFAULT_K, check_overlap: bool = True, rank_key: str = "tm_score", refinement_channel="gray_inverted", refine_predictions: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
        One iteration of a click-to-verify session: ``prepare_iteration``, then ``run_roi`` with every stored box blanked in the
        search (after the template is cut), then ``merge_iteration``. The first iteration blanks nothing, so its run is
        exactly ``run_roi``'s.

        rgb (np.ndarray): the ROI's full-resolution RGB image.
        session (pd.DataFrame or None): every completed iteration's stored rows; None or empty for the first iteration.
        mpp (float): the ROI's microns per pixel.
        seed_df (pd.DataFrame): the human's new seed, one HUMAN MITOTIC_GT row (``click_frame`` or ``box_frame``).
        k (int): how many predictions to show and store; DEFAULT_K (10).
        check_overlap (bool): raise when the new seed's box contains the centre of a stored box.
        rank_key (str): "tm_score" or "chromatin_od".
        refinement_channel (str or None): the Otsu gate's channel, or None for no tightening.
        refine_predictions (bool): store each shown prediction's refined box, which later iterations blank; False blanks the prediction's own 51 px box instead.

        Returns tuple[pd.DataFrame, pd.DataFrame, dict]: the new session; ``run_roi``'s full frame for this iteration (every survivor, for scoring and plots); and its info plus ``merge_iteration``'s keys.
        Raises ValueError on a bad k, session or seed, as prepare_iteration and run_roi do.
    """
    if isinstance(k, bool) or not isinstance(k, (int, np.integer)) or k < 1:
        raise ValueError(f"k must be a positive integer, got {k!r}")
    session, seed_frame, blanks, iteration = prepare_iteration(session, seed_df, check_overlap=check_overlap)
    squares = list(zip(blanks["cx"].tolist(), blanks["cy"].tolist(), blanks["side"].tolist())) if len(blanks) else None
    frame, info = run_roi(rgb, seed_frame, mpp, rank_key=rank_key, refinement_channel=refinement_channel, refine_predictions=refine_predictions, extra_blanks=squares)
    new_session, added = merge_iteration(session, frame, info, blanks, iteration, k=int(k))
    clash = sorted(set(added) & set(info))
    if clash:
        raise RuntimeError(f"merge_iteration's info keys would overwrite run_roi's: {clash}")
    info.update(added)
    return new_session, frame, info


def _joint_spec(channel, hem: np.ndarray, cx: float, cy: float):
    """
        One condition's template spec for one click, as the driver's ``template_spec`` computes it.

        channel (np.ndarray or None): the condition's gate channel, or None for the 51 px default.
        hem (np.ndarray): the ROI's hematoxylin_od channel, for the readability check.
        cx (float): click x.
        cy (float): click y.

        Returns tuple[int, float, float] or None: (base_size, tx, ty), or None if the gate refuses or the template is unreadable.
    """
    got = (DEFAULT_BASE_SIZE, cx, cy) if channel is None else ss.tightened_template_box(channel, cx, cy, otsu_window=OTSU_WINDOW)
    if got is None or tm.read_padded_patch(hem, got[1], got[2], tm.PATCH_SIZE) is None:
        return None
    return int(got[0]), float(got[1]), float(got[2])


def draw_seed_rows(rgb: np.ndarray, gt: pd.DataFrame, image_id: int, *, seed_indices=(0, 1, 2)) -> list[tuple[pd.DataFrame, dict]]:
    """
        The validation draw: a port of ``rerun_bbox3way_postD11.select_clicks``'s joint gate. A click is accepted only if every
        condition in JOINT_CONDITIONS can build a template from it; each seed excludes the clicks earlier seeds accepted.

        rgb (np.ndarray): the ROI's RGB image.
        gt (pd.DataFrame): the ROI's annotations, all categories, as ``dataset.image_annotations`` returns them.
        image_id (int): the ROI's image id; seeds each draw's RNG with [seed_index, image_id].
        seed_indices (iterable[int]): the seeds to draw, in order.

        Returns list[tuple[pd.DataFrame, dict]]: per seed, the one-row human-seed frame and its draw_info.
        Raises ValueError when gt is not one image, fewer clicks pass the joint gate than seeds are asked for, or a pool empties.
    """
    seed_indices = list(seed_indices)
    file_names, image_ids = gt["file_name"].unique(), gt["image_id"].unique()
    if len(file_names) != 1 or len(image_ids) != 1 or int(image_ids[0]) != int(image_id):
        raise ValueError(f"gt must hold one file_name and one image_id equal to image_id={image_id}; found {list(file_names)} and {list(image_ids)}")
    file_name = file_names[0]

    channels = {name: ch.to_channel(rgb, name) for name in ("gray_inverted", "hematoxylin_od")}
    pool, flagged = ss.agreement_pool(gt[gt["category_id"] == ds.MITOTIC])
    n_agreement_pool = len(pool)
    pool = ss.border_filter(pool, tm.PATCH_SIZE // 2, rgb.shape)
    gate_specs = {int(row["ann_id"]): {condition: _joint_spec(None if name is None else channels[name], channels["hematoxylin_od"], float(row["cx"]), float(row["cy"])) for condition, name in JOINT_CONDITIONS.items()} for _, row in pool.iterrows()}
    del channels
    n_joint_valid = sum(all(spec is not None for spec in specs.values()) for specs in gate_specs.values())
    if n_joint_valid < len(seed_indices):
        raise ValueError(f"{file_name}: only {n_joint_valid} pool annotations pass the joint gate, fewer than the {len(seed_indices)} distinct clicks needed")

    out, accepted = [], []
    for seed_index in seed_indices:
        seed_pool = pool[~pool["ann_id"].isin(accepted)] if accepted else pool
        rng = np.random.default_rng([int(seed_index), int(image_id)])
        working, refusals = seed_pool.copy(), []
        while True:
            if len(working) == 0:
                raise ValueError(f"{file_name} seed {seed_index}: seed pool exhausted -- no annotation passed every condition")
            idx = int(rng.integers(len(working)))
            row = working.iloc[idx]
            specs = gate_specs[int(row["ann_id"])]
            if all(spec is not None for spec in specs.values()):
                break
            refusals.append(f"{int(row['ann_id'])}:{'+'.join(c for c, spec in specs.items() if spec is None)}")
            working = working.drop(working.index[idx])
        accepted.append(int(row["ann_id"]))
        draw_info = dict(seed_index=seed_index, seed_ann_id=int(row["ann_id"]), n_agreement_pool=n_agreement_pool, n_after_border=len(pool), n_joint_valid=n_joint_valid, n_seed_pool=len(seed_pool), n_retries=len(refusals), refused_draws=";".join(refusals), contested_seed_tier=bool(flagged), specs=dict(specs))
        out.append((human_seed_frame(row, file_name), draw_info))
    return out
