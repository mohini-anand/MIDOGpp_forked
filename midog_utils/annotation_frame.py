"""
    The annotation frame: one DataFrame carrying every box of a click-to-verify session -- the
    human's click, the template box the algorithm refined it to, and every prediction. `run_roi`
    takes a frame holding the click and returns a new frame with the refined seed and the
    predictions appended; `evaluate` scores that output against MIDOG++ ground truth;
    `click_frame` and `box_frame` build the input row for a raw click or a drawn box, and
    `draw_seed_rows` draws the human-seed rows the validation path feeds `run_roi`. The seed row's
    box sets the Otsu gate's window and the fallback template (`gate_window_side`): 51 px for every
    click and MIDOG row, the drawn box's own size for a box. A layer on top of the frozen modules,
    which do every numerical step. Spec: `ANNOTATION_DATAFRAME_SPEC.md`.
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
HUMAN = "HUMAN"
ALGORITHM = "ALGORITHM"

GATE_REFUSED = "gate_refused"                # the Otsu gate refused the click
TEMPLATE_UNREADABLE = "template_unreadable"  # the gate accepted, but the refined template is within 36 px of an edge

DEFAULT_BASE_SIZE = tm.BASE_SIZE   # 51: the untightened template, and the fallback, for clicks and MIDOG rows (gate_window_side(50, 50))
OTSU_WINDOW = tm.BASE_SIZE         # 51: the validation draw's gate window (_joint_spec)
PREDICTION_BOX = float(ds.BOX_SIZE)  # 50.0: a prediction row's width and height
CLICK_BOX = float(ds.BOX_SIZE)       # 50.0: a raw click's row width and height, MIDOG's own box size; fixes a click's 51 px window and fallback
MAX_WINDOW_PX = tm.PATCH_SIZE      # 73: the largest gate window; a larger human box raises
MIN_WINDOW_PX = 19                 # the smallest gate window, i.e. a drawn box of 18 px or more (DECISIONS.md D12)
REFINEMENT_CHANNELS = (None, "gray_inverted", "hematoxylin_od")
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


def _template_choice(rgb: np.ndarray, cx: float, cy: float, refinement_channel, window: int):
    """
        The template spec for a readable seed: the accepted Otsu component's, or the human's box (``window`` px) centred on the
        seed.

        rgb (np.ndarray): the ROI's RGB image.
        cx (float): seed centre x, unrounded.
        cy (float): seed centre y, unrounded.
        refinement_channel (str or None): the gate's channel, or None for no refinement.
        window (int): the gate window and fallback size, ``gate_window_side`` of the seed row; 51 for clicks and MIDOG rows.

        Returns tuple: (base_size, tx, ty, refined, refusal_reason).
        Raises RuntimeError when the window is not readable at the seed centre, so a refusal always means the gate refused.
    """
    if refinement_channel is None:
        return window, cx, cy, False, None
    if not ss._patch_readable(rgb.shape, cx, cy, window):
        raise RuntimeError(f"the {window} px gate window at ({cx}, {cy}) is not readable; run_roi's edge check should have raised")
    got = ss.tightened_template_box(ch.to_channel(rgb, refinement_channel), cx, cy, otsu_window=window)
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


def run_roi(rgb: np.ndarray, annotations_df: pd.DataFrame, mpp: float, *, seed_object_id=None, rank_key: str = "tm_score", refinement_channel="gray_inverted") -> tuple[pd.DataFrame, dict]:
    """
        One human click or box in, the refined seed and every post-NMS prediction appended. Never draws a seed and never
        re-sorts. The seed row's box sets the gate window and the fallback template (``gate_window_side``).

        rgb (np.ndarray): the ROI's full-resolution RGB image (`dataset.load_roi`).
        annotations_df (pd.DataFrame): the frame's input rows, all HUMAN MITOTIC_GT of one ROI; never mutated.
        mpp (float): the ROI's microns per pixel (`dataset.roi_mpp`).
        seed_object_id (int or None): which input row is the seed; None when there is exactly one.
        rank_key (str): "tm_score" or "chromatin_od" -- the sort order only; score and od are both always computed.
        refinement_channel (str or None): the Otsu gate's channel, "gray_inverted" or "hematoxylin_od", or None for no tightening: the template is the human's box (51 px for clicks and MIDOG).

        Returns tuple[pd.DataFrame, dict]: the input rows, the refined seed and the predictions, and the run's info.
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
    base_size, tx, ty, refined, refusal_reason = _template_choice(rgb, cx, cy, refinement_channel, window)
    t_refine = time.time() - t0
    refine_time = now_utc()

    detections, prod_info = prod.run_production_pipeline(rgb, SimpleNamespace(base_size=base_size, template_xy=(tx, ty)), mpp, rank_key=rank_key, always_od=(rank_key != "chromatin_od"))
    predict_time = now_utc()

    x_tl, y_tl, width, height = footprint_box(base_size, tx, ty)
    refined_row = _block(1, file_name=file_name, x_tl=x_tl, y_tl=y_tl, width=width, height=height, **{"class": MITOTIC_GT}, object_id=int(seed["object_id"]), made_by=ALGORITHM, refined=refined, current_time=refine_time, score=np.nan, od=np.nan)
    out = pd.concat([frame, refined_row, _prediction_rows(detections, file_name, predict_time)], ignore_index=True)
    _check_schema(out)

    added = dict(file_name=file_name, mpp=mpp, match_radius_px=ev.radius_px(mpp), refinement_channel=refinement_channel, refinement_requested=refinement_channel is not None, channel=prod.CHANNEL, tm_method=prod.TM_METHOD, peak_min_distance=prod.PEAK_MIN_DISTANCE, deep_floor_z=prod.DEEP_FLOOR_Z, max_peaks=prod.MAX_PEAKS, patch_size=tm.PATCH_SIZE, base_size=base_size, seed_object_id=int(seed["object_id"]), click_center=(cx, cy), refined=refined, refusal_reason=refusal_reason, template_center=(tx, ty), offset_px=float(np.hypot(tx - cx, ty - cy)), t_refine_s=round(t_refine, 2))
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
