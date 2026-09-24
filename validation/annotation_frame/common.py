"""
    Shared setup for the annotation-frame validation harness: pins the working directory and
    import path to the worktree, loads the tracked references, and holds the exact comparator and
    the driver-formula rebuild of one R1 row. Every harness script imports this first.
"""

from __future__ import annotations

import os
import pickle
import sys
from pathlib import Path

WORKTREE = Path(__file__).resolve().parents[2]
os.chdir(WORKTREE)
sys.path.insert(0, ".")
sys.path.insert(0, str(WORKTREE))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import midog_utils  # noqa: E402
from midog_utils import channels as ch  # noqa: E402,F401
from midog_utils import dataset as ds  # noqa: E402
from midog_utils import evaluate as ev  # noqa: E402
from midog_utils import production as prod  # noqa: E402
from midog_utils import seed_selection as ss  # noqa: E402,F401
from midog_utils import template_match as tm  # noqa: E402,F401
import rerun_bbox3way_postD11 as rr  # noqa: E402

for _mod in (midog_utils, prod, rr):
    assert str(Path(_mod.__file__).resolve()).startswith(str(WORKTREE) + os.sep), f"{_mod.__name__} imported from {_mod.__file__}, outside the worktree {WORKTREE}"

OUT = WORKTREE / "validation" / "annotation_frame" / "out"
STEP0 = OUT / "step0"

R1_PATH = "results/precision_at_k_49roi_3seed_chromatin_bbox3way_postD11_per_run.csv"
R2_PATH = "results/precision_at_k_49roi_3seed_chromatin_bbox3way_postD11_top30.csv"
R3_PATH = "results/precision_at_k_49roi_3seed_chromatin_bbox3way_per_run.csv"
R4_PER_RUN_PATH = "results/precision_at_k_49roi_6seed_fillgeom_postD11_per_run.csv"
R4_TOP30_PATH = "results/precision_at_k_49roi_6seed_fillgeom_postD11_top30.csv"
R4_CLICKS_PATH = "results/precision_at_k_49roi_6seed_fillgeom_postD11_clicks.csv"

CONDITION_CHANNEL = dict(rr.CONDITIONS)  # default_51 -> None, gray_bbox -> gray_inverted, hem_bbox -> hematoxylin_od
TOP_N = max(rr.BUDGETS)

# V9's definition, frozen before the edit: the first original_14 ROI of each of the first five domains in build_rois order,
# a 1024 px crop of the full-ROI hematoxylin channel around the reference seed-0 click, clamped inside the image.
V9_CROP = 1024
V9_N_ROIS = 5


def read_ref(path):
    """
        Read a tracked reference CSV the way every check must.

        path (str): repo-relative CSV path.

        Returns pd.DataFrame: the CSV, parsed with round_trip float precision.
    """
    return pd.read_csv(path, float_precision="round_trip")


def v9_rois(rois):
    """
        The five original_14 ROIs V9 crops: the first original_14 ROI of each of the first five domains, in build_rois order.

        rois (pd.DataFrame): rr.build_rois()'s ROI table.

        Returns list[str]: five file names.
    """
    orig = rois[rois["subset"] == "original_14"]
    firsts = orig.drop_duplicates("domain", keep="first")
    return firsts["file_name"].tolist()[:V9_N_ROIS]


def v9_crop_box(shape, cx, cy):
    """
        The integer top-left of V9's 1024 px crop around a click, clamped inside the image.

        shape (tuple): the ROI's array shape.
        cx (float): click x, integer-valued.
        cy (float): click y, integer-valued.

        Returns tuple[int, int]: (x0, y0) of the crop in full-ROI pixels.
    """
    assert float(cx).is_integer() and float(cy).is_integer(), f"V9 expects an integer reference click, got ({cx}, {cy})"
    h, w = shape[:2]
    x0 = min(max(int(cx) - V9_CROP // 2, 0), w - V9_CROP)
    y0 = min(max(int(cy) - V9_CROP // 2, 0), h - V9_CROP)
    return x0, y0


def v9_config(mpp):
    """
        The augmented-bank FSConfig V9 runs, exactly as the hand-off prompt names it.

        mpp (float): the ROI's microns per pixel.

        Returns FSConfig: the 8-template config.
    """
    from midog_utils.find_and_suppress import FSConfig
    return FSConfig(scales=(0.8, 1.2), n_angles=2, flips=(False, True), deep_floor_z=-1.5, nms_radius=ev.radius_px(mpp), max_peaks=100)


def save_pickle(obj, path):
    """
        Write a pickle atomically, so an interrupted run never leaves a truncated file.

        obj: anything picklable.
        path (Path): destination.

        Returns None.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "wb") as fh:
        pickle.dump(obj, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, path)


def load_pickle(path):
    """
        Read a pickle written by save_pickle.

        path (Path): source.

        Returns object: the unpickled object.
    """
    with open(path, "rb") as fh:
        return pickle.load(fh)


def _kind(series):
    """
        The comparison kind of a column: 'num' for int/float, 'bool', 'str' for object, else its dtype string.

        series (pd.Series): the column.

        Returns str: the kind.
    """
    k = series.dtype.kind
    if k in "iuf":
        return "num"
    if k == "b":
        return "bool"
    if k == "O":
        return "str"
    return str(series.dtype)


def compare_frames(mine, ref, columns, *, float32_cols=(), label="", max_detail=10):
    """
        Exact, positional, cell-by-cell comparison of two frames: float32 on the named score columns, NaN equals NaN, strings after fillna("").

        mine (pd.DataFrame): the frame under test.
        ref (pd.DataFrame): the reference frame.
        columns (list[str]): the columns to compare; every other column is excluded by name.
        float32_cols (tuple[str]): columns compared after casting both sides to float32.
        label (str): identifies the comparison in the detail lines.
        max_detail (int): how many diverging cells to describe.

        Returns dict: n_rows, n_columns, n_values, n_diverged and a detail list.
    """
    mine = mine.reset_index(drop=True)
    ref = ref.reset_index(drop=True)
    detail = []
    if len(mine) != len(ref):
        return dict(label=label, n_rows=max(len(mine), len(ref)), n_columns=len(columns), n_values=max(len(mine), len(ref)) * len(columns), n_diverged=max(len(mine), len(ref)) * len(columns), detail=[f"{label}: row count {len(mine)} != {len(ref)}"])
    n_div = 0
    for c in columns:
        if c not in mine.columns or c not in ref.columns:
            n_div += len(mine)
            detail.append(f"{label}: column {c!r} missing from {'mine' if c not in mine.columns else 'ref'}")
            continue
        a, b = mine[c], ref[c]
        if c in float32_cols:
            av, bv = a.to_numpy(dtype=np.float64).astype(np.float32), b.to_numpy(dtype=np.float64).astype(np.float32)
            ok = (av == bv) | (np.isnan(av) & np.isnan(bv))
        else:
            ka, kb = _kind(a), _kind(b)
            if ka != kb:
                ok = np.zeros(len(a), dtype=bool)
                detail.append(f"{label}: column {c!r} kind {ka} ({a.dtype}) != {kb} ({b.dtype})")
            elif ka == "num":
                av, bv = a.to_numpy(dtype=np.float64), b.to_numpy(dtype=np.float64)
                ok = (av == bv) | (np.isnan(av) & np.isnan(bv))
                if a.dtype.kind in "iu" or b.dtype.kind in "iu":
                    ok &= a.to_numpy() == b.to_numpy()  # ints compared as ints too, not only through float64
            elif ka == "bool":
                ok = a.to_numpy() == b.to_numpy()
            elif ka == "str":
                ok = np.array([x == y and type(x) is type(y) for x, y in zip(a.fillna(""), b.fillna(""))], dtype=bool)
            else:
                ok = (a.to_numpy() == b.to_numpy())
        bad = np.flatnonzero(~ok)
        n_div += len(bad)
        for i in bad[:max(0, max_detail - len(detail))]:
            detail.append(f"{label}: row {i} column {c!r}: {a.iat[i]!r} != {b.iat[i]!r}")
    return dict(label=label, n_rows=len(mine), n_columns=len(columns), n_values=len(mine) * len(columns), n_diverged=n_div, detail=detail)


def same_value(a, b):
    """
        Exact equality of two info values: same type and equal, with NaN equal to NaN and tuples compared element-wise.

        a: first value.
        b: second value.

        Returns bool: True when the two are identical.
    """
    if isinstance(a, tuple) and isinstance(b, tuple):
        return len(a) == len(b) and all(same_value(x, y) for x, y in zip(a, b))
    if type(a) is not type(b):
        return False
    if isinstance(a, float) and np.isnan(a) and np.isnan(b):
        return True
    return bool(a == b)


def compare_info(mine, ref, *, exclude=(), label=""):
    """
        Key-for-key, value-for-value comparison of two info dicts, excluding every t_* key and the named keys.

        mine (dict): the info under test.
        ref (dict): the reference info.
        exclude (tuple[str]): further keys to leave out.
        label (str): identifies the comparison in the detail lines.

        Returns dict: n_values, n_diverged and a detail list.
    """
    keys_a = {k for k in mine if not k.startswith("t_") and k not in exclude}
    keys_b = {k for k in ref if not k.startswith("t_") and k not in exclude}
    detail = [f"{label}: key {k!r} only in {'mine' if k in keys_a else 'ref'}" for k in sorted(keys_a ^ keys_b)]
    n_div = len(keys_a ^ keys_b)
    for k in sorted(keys_a & keys_b):
        if not same_value(mine[k], ref[k]):
            n_div += 1
            detail.append(f"{label}: info[{k!r}] {mine[k]!r} != {ref[k]!r}")
    return dict(label=label, n_values=len(keys_a | keys_b), n_diverged=n_div, detail=detail)


def driver_row(roi, seed_index, ann, specs, refusals, n_seed_pool, digests, condition, detections, info, gt, pool, flagged, n_agreement, n_joint_valid, mpp):
    """
        Rebuild one R1 row and its R2 rows from raw production output, with the driver's own formulas (rerun_bbox3way_postD11.py:230-252).

        roi (pd.Series): the ROI's build_rois row.
        seed_index (int): the draw's seed index.
        ann (pd.Series): the accepted annotation.
        specs (dict): condition -> (base_size, tx, ty).
        refusals (list[str]): the draw's refusal entries.
        n_seed_pool (int): the draw's pool size.
        digests (dict): condition -> template SHA-1.
        condition (str): the run's condition.
        detections (pd.DataFrame): production's ranked survivors.
        info (dict): production's info.
        gt (pd.DataFrame): the ROI's ground truth.
        pool (pd.DataFrame): the border-filtered pool.
        flagged (bool): whether the pool is the contested tier.
        n_agreement (int): agreement-tier pool size.
        n_joint_valid (int): joint-valid pool size.
        mpp (float): microns per pixel.

        Returns tuple[dict, pd.DataFrame]: the per-run row (t_pipeline_s left out) and the top-30 frame.
    """
    fn = roi["file_name"]
    match_radius = ev.radius_px(mpp)
    click_xy = (float(ann["cx"]), float(ann["cy"]))
    seed_ann_id = int(ann["ann_id"])
    gt_eval = gt[gt["ann_id"] != seed_ann_id].reset_index(drop=True)
    n_gt = int((gt_eval["category_id"] == ds.MITOTIC).sum())
    nearest_other = float(np.hypot(gt_eval["cx"] - click_xy[0], gt_eval["cy"] - click_xy[1]).min())
    base_size, tx, ty = specs[condition]
    offset = float(np.hypot(tx - click_xy[0], ty - click_xy[1]))
    det_out, tp = rr.score_top_k(detections, gt_eval, match_radius)
    d_click = np.hypot(det_out["cx"] - click_xy[0], det_out["cy"] - click_xy[1]).to_numpy()
    d_tpl = np.hypot(det_out["cx"] - tx, det_out["cy"] - ty).to_numpy()
    near = np.flatnonzero(d_click <= match_radius)
    near_detail = ";".join(f"rank{int(det_out['rank'].iat[i])}:{det_out['bucket'].iat[i]}:ann{int(det_out['matched_ann_id'].iat[i])}:d_click={d_click[i]:.1f}:d_tpl={d_tpl[i]:.1f}" for i in near)
    od = det_out["od"].to_numpy()
    row = dict(subset=roi["subset"], file_name=fn, domain=roi["domain"], image_id=int(roi["image_id"]), seed_index=seed_index, condition=condition, seed_ann_id=seed_ann_id, click_cx=click_xy[0], click_cy=click_xy[1], base_size=int(base_size), tpl_cx=tx, tpl_cy=ty, tpl_offset_px=round(offset, 3), tpl_sha1=digests[condition], n_agreement_pool=n_agreement, n_after_border=len(pool), n_joint_valid=n_joint_valid, n_seed_pool=n_seed_pool, n_retries=len(refusals), refused_draws=";".join(refusals), contested_seed_tier=bool(flagged), n_gt_mitotic=n_gt, nearest_other_ann_px=round(nearest_other, 2), n_peaks=int(info["n_peaks"]), max_peaks_binding=bool(info["max_peaks_binding"]), n_detections=int(info["n_detections"]), n_blanked_px=int(info["n_blanked_px"]), n_near_click=len(near), near_click_detail=near_detail, od_nan=int(np.isnan(od).sum()), od_descending=bool(np.all(np.diff(od) <= 0)), match_radius_px=match_radius, mpp=mpp, **{f"tp_at_{k}": tp[k] for k in rr.BUDGETS})
    top = det_out.head(TOP_N).assign(subset=roi["subset"], file_name=fn, seed_index=seed_index, condition=condition)
    return row, top


def finish_per_run(rows, rois):
    """
        Add precision_at_k and sort rows the way the driver's main() does before writing R1.

        rows (list[dict]): per-run rows.
        rois (pd.DataFrame): build_rois()'s ROI table, in its order.

        Returns pd.DataFrame: the per-run frame.
    """
    results = pd.DataFrame(rows)
    for k in rr.BUDGETS:
        results[f"precision_at_{k}"] = results[f"tp_at_{k}"] / k
    roi_order = rois["file_name"].tolist()
    results["_roi_o"] = results["file_name"].map({fn: i for i, fn in enumerate(roi_order)})
    results["_cond_o"] = results["condition"].map({c: i for i, c in enumerate(rr.CONDITIONS)})
    return results.sort_values(["_roi_o", "seed_index", "_cond_o"]).drop(columns=["_roi_o", "_cond_o"]).reset_index(drop=True)


def summarize(results):
    """
        Collapse a list of comparison dicts into totals.

        results (list[dict]): compare_frames / compare_info outputs.

        Returns dict: n_values, n_diverged and the first detail lines.
    """
    return dict(n_values=int(sum(r["n_values"] for r in results)), n_diverged=int(sum(r["n_diverged"] for r in results)), detail=[d for r in results for d in r["detail"]][:20])
