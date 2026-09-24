"""
    The checks with their own named inputs (spec section 7.2): V5 (gate fallback and template_unreadable), V6 (half-integer
    clicks), V7 (edges and inputs), V8 (schema, int64 input), V9 (augmented bank, after the edit vs Step 0) and V1's
    contested-tier path. The smoke tier runs V8 on one all-int ROI and V9 on one crop; the full tier runs all of them.
    Every candidate click is re-verified before use. Usage: named_checks.py --tier smoke|full
"""

from __future__ import annotations

import argparse
import json
import time
from types import SimpleNamespace

import numpy as np
import pandas as pd

import common as C
from common import ch, ds, prod, rr, ss, tm
import checks
from checks import Tally
from midog_utils import annotation_frame as af
from midog_utils.find_and_suppress import find_and_suppress

RAW = json.loads(open("databases/MIDOG++.json").read())
V5_CANDIDATES = [("462.tiff", 22673, "gray_inverted"), ("529.tiff", 24934, "gray_inverted"), ("300.tiff", 14581, "hematoxylin_od"), ("245.tiff", 6274, "hematoxylin_od")]
V5_UNREADABLE = ("300.tiff", 6410.0, 250.0)
V6_CANDIDATES = [("072.tiff", 1835), ("007.tiff", 144)]
V8_ALL_INT = ["211.tiff", "364.tiff", "529.tiff"]
CONTESTED_CANDIDATES = ["201.tiff", "233.tiff", "013.tiff"]
OTHER = {"gray_inverted": "hematoxylin_od", "hematoxylin_od": "gray_inverted"}
RANK_KEYS = ("chromatin_od", "tm_score")


class Rois:
    """
        Loads each ROI once per process and keeps the last few.
    """

    def __init__(self):
        self.rois, self.annotations = rr.build_rois()
        self.cache = {}

    def get(self, fn):
        """
            The ROI's (roi row, rgb, mpp, gt).

            fn (str): file name.

            Returns tuple: (pd.Series, np.ndarray, float, pd.DataFrame).
        """
        if fn not in self.cache:
            if len(self.cache) >= 2:
                self.cache.pop(next(iter(self.cache)))
            roi = self.rois[self.rois["file_name"] == fn].iloc[0]
            self.cache[fn] = (roi, ds.load_roi(roi["path"]), ds.roi_mpp(roi["path"]), ds.image_annotations(self.annotations, fn))
        return self.cache[fn]


def synthetic_seed(fn, cx, cy, object_id=1):
    """
        A one-row human-seed frame for a synthetic click, width 50 on the half-pixel grid.

        fn (str): file name.
        cx (float): click x.
        cy (float): click y.
        object_id (int): the row's id.

        Returns pd.DataFrame: the frame.
    """
    return pd.DataFrame({"file_name": pd.Series([fn], dtype=object), "x_tl": [cx - 25.0], "y_tl": [cy - 25.0], "width": [50.0], "height": [50.0], "class": pd.Series(["MITOTIC_GT"], dtype=object), "object_id": np.array([object_id], dtype=np.int64), "made_by": pd.Series(["HUMAN"], dtype=object), "refined": np.array([False]), "current_time": pd.Series([af.now_utc()], dtype=object), "score": [np.nan], "od": [np.nan]})


def preds(frame):
    """
        A frame's prediction rows.

        frame (pd.DataFrame): run_roi's output.

        Returns pd.DataFrame: the PREDICTION rows.
    """
    return frame[frame["class"] == "PREDICTION"]


def gate(rgb, channel, cx, cy):
    """
        The gate's verdict on a click: the accepted spec or None, and whether the refined template is readable.

        rgb (np.ndarray): the ROI.
        channel (str): gate channel.
        cx (float): click x.
        cy (float): click y.

        Returns tuple: (spec or None, readable or None).
    """
    got = ss.tightened_template_box(ch.to_channel(rgb, channel), cx, cy, otsu_window=tm.BASE_SIZE)
    return got, (None if got is None else ss._patch_readable(rgb.shape, got[1], got[2], tm.PATCH_SIZE))


def expect_raise(tally, name, fn_call, what):
    """
        Record whether a call raises, with its exception type, message and elapsed time.

        tally (Tally): where to record.
        name (str): the check.
        fn_call (callable): the call.
        what (str): what the case is.

        Returns dict: the case record.
    """
    t0 = time.time()
    try:
        fn_call()
        raised, err = False, None
    except Exception as exc:  # noqa: BLE001
        raised, err = True, f"{type(exc).__name__}: {exc}"
    rec = dict(case=what, raised=raised, error=err, seconds=round(time.time() - t0, 3))
    tally.flag(name, raised, f"{what}: did not raise")
    return rec


def v5(tally, R, log):
    """
        V5: a refused click falls back to the default_51 template exactly; the other gate accepts; template_unreadable.

        tally (Tally): where to record.
        R (Rois): the ROI loader.
        log (dict): per-case records.

        Returns None.
    """
    cases = []
    for fn, ann_id, refused in V5_CANDIDATES:
        roi, rgb, mpp, gt = R.get(fn)
        ann = gt[gt["ann_id"] == ann_id].iloc[0]
        cx, cy = float(ann["cx"]), float(ann["cy"])
        g_ref, _ = gate(rgb, refused, cx, cy)
        g_oth, oth_readable = gate(rgb, OTHER[refused], cx, cy)
        verified = g_ref is None and g_oth is not None and oth_readable and ss._patch_readable(rgb.shape, cx, cy, tm.PATCH_SIZE)
        rec = dict(file_name=fn, ann_id=ann_id, refused_channel=refused, click=(cx, cy), refused_gate=g_ref, other_gate=g_oth, other_readable=oth_readable, verified=bool(verified))
        tally.flag("V5 candidate re-verified", verified, f"{fn} ann {ann_id}: refused gate {g_ref}, other gate {g_oth} readable {oth_readable}")
        if not verified:
            cases.append(rec)
            continue
        seed = af.human_seed_frame(ann, fn)
        for rank_key in RANK_KEYS:
            f_none, i_none = af.run_roi(rgb, seed, mpp, rank_key=rank_key, refinement_channel=None)
            f_ref, i_ref = af.run_roi(rgb, seed, mpp, rank_key=rank_key, refinement_channel=refused)
            f_oth, i_oth = af.run_roi(rgb, seed, mpp, rank_key=rank_key, refinement_channel=OTHER[refused])
            lbl = f"{fn} ann {ann_id} {rank_key}"
            tally.flag("V5 gate_refused: refined False and exact reason", i_ref["refined"] is False and i_ref["refusal_reason"] == "gate_refused" and bool(f_ref["refined"].iat[1]) is False, f"{lbl}: refined {i_ref['refined']!r} reason {i_ref['refusal_reason']!r}", runs=1)
            tally.add("V5 gate_refused: predictions identical to refinement_channel=None", C.compare_frames(preds(f_ref), preds(f_none), checks.COMPARED_FRAME_COLUMNS, label=f"{lbl} fallback"), runs=1)
            tally.flag("V5 extra: fallback template is default_51 at the click", i_ref["base_size"] == 51 and C.same_value(i_ref["template_center"], (cx, cy)) and f_ref.iloc[1][["x_tl", "y_tl", "width", "height"]].equals(f_none.iloc[1][["x_tl", "y_tl", "width", "height"]]) and i_none["refusal_reason"] is None and i_none["refined"] is False and i_ref["refinement_requested"] is True and i_none["refinement_requested"] is False, f"{lbl}: fallback template {i_ref['base_size']} {i_ref['template_center']}")
            tally.flag("V5 other channel: refined True", i_oth["refined"] is True and i_oth["refusal_reason"] is None and bool(f_oth["refined"].iat[1]) is True, f"{lbl}: other channel refined {i_oth['refined']!r}", runs=1)
            rec[f"n_pred_{rank_key}"] = int(len(preds(f_ref)))
        cases.append(rec)

    fn, cx, cy = V5_UNREADABLE
    roi, rgb, mpp, gt = R.get(fn)
    gates = {c: gate(rgb, c, cx, cy) for c in ("gray_inverted", "hematoxylin_od")}
    verified = ss._patch_readable(rgb.shape, cx, cy, tm.PATCH_SIZE) and all(g is not None and r is False for g, r in gates.values())
    rec = dict(file_name=fn, click=(cx, cy), width=int(rgb.shape[1]), x_from_right=int(rgb.shape[1] - cx), gates={c: dict(spec=g, readable=r) for c, (g, r) in gates.items()}, verified=bool(verified))
    tally.flag("V5 template_unreadable click re-verified", verified, f"{fn} ({cx}, {cy}): gates {gates}")
    if verified:
        seed = synthetic_seed(fn, cx, cy)
        for rank_key in RANK_KEYS:
            f_none, _ = af.run_roi(rgb, seed, mpp, rank_key=rank_key, refinement_channel=None)
            for channel in ("gray_inverted", "hematoxylin_od"):
                f_c, i_c = af.run_roi(rgb, seed, mpp, rank_key=rank_key, refinement_channel=channel)
                lbl = f"{fn} ({cx}, {cy}) {channel} {rank_key}"
                tally.flag("V5 template_unreadable: refined False and exact reason", i_c["refined"] is False and i_c["refusal_reason"] == "template_unreadable" and bool(f_c["refined"].iat[1]) is False, f"{lbl}: refined {i_c['refined']!r} reason {i_c['refusal_reason']!r}", runs=1)
                tally.add("V5 template_unreadable: predictions identical to refinement_channel=None", C.compare_frames(preds(f_c), preds(f_none), checks.COMPARED_FRAME_COLUMNS, label=f"{lbl} unreadable"), runs=1)
    cases.append(rec)
    log["V5"] = cases


def v6(tally, R, log):
    """
        V6: half-integer clicks -- the human row reproduces the JSON, the refined row rounds half to even, the blank size, and
        run_roi's predictions against production called directly with a template_spec seed.

        tally (Tally): where to record.
        R (Rois): the ROI loader.
        log (dict): per-case records.

        Returns None.
    """
    cases = []
    for fn, ann_id in V6_CANDIDATES:
        roi, rgb, mpp, gt = R.get(fn)
        ann = gt[gt["ann_id"] == ann_id].iloc[0]
        cx, cy = float(ann["cx"]), float(ann["cy"])
        channels = {name: ch.to_channel(rgb, name) for name in ("gray_inverted", "hematoxylin_od")}
        gates = {c: gate(rgb, c, cx, cy) for c in ("gray_inverted", "hematoxylin_od")}
        half = (cx % 1 == 0.5) or (cy % 1 == 0.5)
        verified = half and all(g is not None for g, _ in gates.values())
        rec = dict(file_name=fn, ann_id=ann_id, click=(cx, cy), gates={c: dict(spec=g, readable=r) for c, (g, r) in gates.items()}, half_integer=half, verified=bool(verified), x_half_even=int(round(cx)), x_half_up=int(np.floor(cx + 0.5)), y_half_even=int(round(cy)), y_half_up=int(np.floor(cy + 0.5)))
        tally.flag("V6 candidate re-verified", verified, f"{fn} ann {ann_id}: half {half} gates {gates}")
        if not verified:
            cases.append(rec)
            continue
        raw = next(a for a in RAW["annotations"] if a["id"] == ann_id)
        seed = af.human_seed_frame(ann, fn)
        tally.flag("V6 human row x_tl/y_tl == JSON x1/y1", seed["x_tl"].iat[0] == float(raw["bbox"][0]) and seed["y_tl"].iat[0] == float(raw["bbox"][1]), f"{fn} ann {ann_id}: ({seed['x_tl'].iat[0]}, {seed['y_tl'].iat[0]}) != {raw['bbox'][:2]}", runs=1, rows=1)
        rec["json_bbox"] = raw["bbox"]
        rec["templates"] = {}
        for condition, channel in C.CONDITION_CHANNEL.items():
            spec = rr.template_spec(condition, channels, cx, cy)
            for rank_key in RANK_KEYS:
                frame, info = af.run_roi(rgb, seed, mpp, rank_key=rank_key, refinement_channel=channel)
                lbl = f"{fn} ann {ann_id} {condition} {rank_key}"
                b, tx, ty = info["base_size"], *info["template_center"]
                tally.flag("V6 template matches template_spec", spec is not None and C.same_value((b, tx, ty), spec), f"{lbl}: run_roi template {(b, tx, ty)} != template_spec {spec}")
                if spec is None:
                    continue
                ref_row = frame.iloc[1]
                tally.flag("V6 refined row rounds half to even", ref_row["x_tl"] == float(int(round(tx)) - b // 2) and ref_row["y_tl"] == float(int(round(ty)) - b // 2) and ref_row["width"] == float(b - 1), f"{lbl}: refined ({ref_row['x_tl']}, {ref_row['y_tl']}) for template ({b}, {tx}, {ty})", runs=1, rows=1)
                tally.flag("V6 n_blanked_px == base_size ** 2", info["n_blanked_px"] == b ** 2, f"{lbl}: n_blanked_px {info['n_blanked_px']} != {b ** 2}")
                det, _ = prod.run_production_pipeline(rgb, SimpleNamespace(base_size=spec[0], template_xy=(spec[1], spec[2])), mpp, rank_key=rank_key, always_od=(rank_key != "chromatin_od"))
                tally.add("V6 predictions vs production with a template_spec seed", C.compare_frames(preds(frame), checks.convert_detections(det, fn), checks.COMPARED_FRAME_COLUMNS, label=f"{lbl} V6"), runs=1)
                rec["templates"][f"{condition}/{rank_key}"] = dict(base_size=b, tx=tx, ty=ty, refined_x_tl=float(ref_row["x_tl"]), refined_y_tl=float(ref_row["y_tl"]), tx_half_up_differs=int(round(tx)) != int(np.floor(tx + 0.5)), ty_half_up_differs=int(round(ty)) != int(np.floor(ty + 0.5)))
        del channels
        cases.append(rec)
    log["V6"] = cases


def v7(tally, R, log):
    """
        V7: edges (x = 35 / 36 and w - 36 / w - 37, plus the same on y as an extra) and every input-validation raise.

        tally (Tally): where to record.
        R (Rois): the ROI loader.
        log (dict): per-case records.

        Returns None.
    """
    fn = "300.tiff"
    roi, rgb, mpp, gt = R.get(fn)
    h, w = rgb.shape[:2]
    cases = []
    ymid, xmid = float(h // 2), float(w // 2)
    for name, (cx, cy), should_raise in [("V7 edges", (35.0, ymid), True), ("V7 edges", (36.0, ymid), False), ("V7 edges", (float(w - 36), ymid), True), ("V7 edges", (float(w - 37), ymid), False), ("V7 extra: edges on y", (xmid, 35.0), True), ("V7 extra: edges on y", (xmid, 36.0), False), ("V7 extra: edges on y", (xmid, float(h - 36)), True), ("V7 extra: edges on y", (xmid, float(h - 37)), False)]:
        what = f"{fn} click ({cx}, {cy}) of ({w}, {h})"
        if should_raise:
            cases.append(expect_raise(tally, name, lambda: af.run_roi(rgb, synthetic_seed(fn, cx, cy), mpp), what))
        else:
            t0 = time.time()
            try:
                frame, info = af.run_roi(rgb, synthetic_seed(fn, cx, cy), mpp)
                ok, err = True, None
                extra = dict(refined=info["refined"], refusal_reason=info["refusal_reason"], n_rows=len(frame))
            except Exception as exc:  # noqa: BLE001
                ok, err, extra = False, f"{type(exc).__name__}: {exc}", {}
            tally.flag(name, ok, f"{what}: raised {err}")
            cases.append(dict(case=what, raised=not ok, error=err, seconds=round(time.time() - t0, 3), **extra))

    ann = gt[gt["ann_id"] == 14499].iloc[0]
    seed = af.human_seed_frame(ann, fn)
    frame, info = af.run_roi(rgb, seed, mpp, refinement_channel=None)
    two = pd.concat([seed, af.human_seed_frame(gt[gt["ann_id"] == 14578].iloc[0], fn)], ignore_index=True)
    dup = pd.concat([seed, seed], ignore_index=True)
    not_human = seed.assign(made_by="ALGORITHM")
    not_gt = seed.assign(**{"class": "PREDICTION"})
    other_file = pd.concat([seed, af.human_seed_frame(gt[gt["ann_id"] == 14578].iloc[0], "301.tiff")], ignore_index=True)
    for name, call, what in [
        ("V7 zero or two seed rows raise", lambda: af.run_roi(rgb, seed.iloc[0:0], mpp), "zero input rows"),
        ("V7 zero or two seed rows raise", lambda: af.run_roi(rgb, seed, mpp, seed_object_id=999999), "seed_object_id matching no row"),
        ("V7 zero or two seed rows raise", lambda: af.run_roi(rgb, two, mpp), "two HUMAN MITOTIC_GT rows, seed_object_id None"),
        ("V7 non-HUMAN-MITOTIC_GT input rows raise", lambda: af.run_roi(rgb, frame, mpp, seed_object_id=14499), "run_roi's own output fed back in"),
        ("V7 non-HUMAN-MITOTIC_GT input rows raise", lambda: af.run_roi(rgb, pd.concat([seed, frame.iloc[1:2]], ignore_index=True), mpp, seed_object_id=14499), "input plus the ALGORITHM refined row"),
        ("V7 non-HUMAN-MITOTIC_GT input rows raise", lambda: af.run_roi(rgb, not_human, mpp), "a MITOTIC_GT row made_by ALGORITHM"),
        ("V7 non-HUMAN-MITOTIC_GT input rows raise", lambda: af.run_roi(rgb, not_gt, mpp), "a HUMAN row of class PREDICTION"),
        ("V7 unknown rank_key / refinement_channel raise", lambda: af.run_roi(rgb, seed, mpp, rank_key="od"), "rank_key='od'"),
        ("V7 unknown rank_key / refinement_channel raise", lambda: af.run_roi(rgb, seed, mpp, refinement_channel="gray"), "refinement_channel='gray'"),
        ("V7 duplicate object_ids raise", lambda: af.run_roi(rgb, dup, mpp, seed_object_id=14499), "two rows with object_id 14499"),
        ("V7 extra: rows from two file_names raise", lambda: af.run_roi(rgb, other_file, mpp, seed_object_id=14499), "a second row from another file"),
        ("V7 extra: missing column raises", lambda: af.run_roi(rgb, seed.drop(columns="od"), mpp), "input without the od column"),
    ]:
        cases.append(expect_raise(tally, name, call, what))

    ok_scored, _, ok_metrics = af.evaluate(frame, gt, mpp, rgb.shape)
    tally.flag("V7 evaluate raises on a renumbered seed", ok_metrics["seed_object_id"] == 14499, "positive control: the valid frame did not evaluate")
    look = int(gt.loc[gt["category_id"] == ds.LOOKALIKE, "ann_id"].iloc[0])
    other_mit = int(gt.loc[(gt["category_id"] == ds.MITOTIC) & (gt["ann_id"] != 14499), "ann_id"].iloc[0])
    for new_id, what in [(1, "object_id 1 (an ann_id of image 1)"), (look, f"object_id {look}, a look-alike in {fn}"), (other_mit, f"object_id {other_mit}, a mitotic ann_id in {fn} at another centre"), (999999, "object_id 999999, no such ann_id")]:
        bad = frame.copy()
        bad.loc[(bad["made_by"] == "HUMAN"), "object_id"] = new_id
        cases.append(expect_raise(tally, "V7 evaluate raises on a renumbered seed", lambda: af.evaluate(bad, gt, mpp, rgb.shape), f"evaluate with the seed row renumbered to {what}"))
    log["V7"] = cases


def v8(tally, R, log, names):
    """
        V8: schema on a frame built from the raw JSON of an all-int image, with int64 coordinate columns.

        tally (Tally): where to record.
        R (Rois): the ROI loader.
        log (dict): per-case records.
        names (list[str]): the all-int ROIs to use.

        Returns None.
    """
    cases = []
    for fn in names:
        roi, rgb, mpp, gt = R.get(fn)
        image_id = int(roi["image_id"])
        anns = [a for a in RAW["annotations"] if a["image_id"] == image_id]
        all_int = all(type(v) is int for a in anns for v in a["bbox"])
        tally.flag("V8 image is all-int in the raw JSON", all_int, f"{fn}: not all-int")
        pick = next(a for a in anns if a["category_id"] == ds.MITOTIC and ss._patch_readable(rgb.shape, (a["bbox"][0] + a["bbox"][2]) / 2, (a["bbox"][1] + a["bbox"][3]) / 2, tm.PATCH_SIZE))
        x1, y1, x2, y2 = pick["bbox"]
        inp = pd.DataFrame({"file_name": [fn], "x_tl": [x1], "y_tl": [y1], "width": [x2 - x1], "height": [y2 - y1], "class": ["MITOTIC_GT"], "object_id": [pick["id"]], "made_by": ["HUMAN"], "refined": [False], "current_time": [af.now_utc()], "score": [np.nan], "od": [np.nan]})
        tally.flag("V8 input coordinate columns are int64", all(str(inp[c].dtype) == "int64" for c in ("x_tl", "y_tl", "width", "height")), f"{fn}: input dtypes {inp.dtypes.to_dict()}")
        before = inp.copy(deep=True)
        frame, info = af.run_roi(rgb, inp, mpp)
        lbl = f"{fn} ann {pick['id']}"
        checks.schema_checks(tally, "V8 schema", frame, info, before, inp, info["refined"], lbl)
        tally.runs("V8 schema")
        row = frame.iloc[0]
        cast_ok = all(str(frame[c].dtype) == "float64" and type(row[c]) is np.float64 and row[c] == float(v) for c, v in (("x_tl", x1), ("y_tl", y1), ("width", x2 - x1), ("height", y2 - y1)))
        tally.flag("V8 int64 input comes back float64 with the same values", cast_ok and str(inp["x_tl"].dtype) == "int64", f"{lbl}: row {row[['x_tl', 'y_tl', 'width', 'height']].to_dict()}", runs=1, rows=1)
        cases.append(dict(file_name=fn, ann_id=pick["id"], bbox=pick["bbox"], all_int=all_int, n_rows=len(frame), refined=info["refined"], refusal_reason=info["refusal_reason"], dtypes={c: str(t) for c, t in frame.dtypes.items()}))
    log["V8"] = cases


def v9(tally, R, log, n):
    """
        V9: the augmented-bank find_and_suppress on the frozen 1024 px crops, after the edit, against Step 0's capture.

        tally (Tally): where to record.
        R (Rois): the ROI loader.
        log (dict): per-case records.
        n (int): how many of the five crops.

        Returns None.
    """
    cases = []
    for fn in C.v9_rois(R.rois)[:n]:
        s0 = C.load_pickle(C.STEP0 / f"{fn}.pkl")
        pre = s0["v9_pre"]
        roi, rgb, mpp, gt = R.get(fn)
        ann = s0["draws"][0][1]
        x0, y0 = C.v9_crop_box(rgb.shape, float(ann["cx"]), float(ann["cy"]))
        hem = ch.to_channel(rgb, "hematoxylin_od")
        crop = np.ascontiguousarray(hem[y0:y0 + C.V9_CROP, x0:x0 + C.V9_CROP])
        del hem
        click = (float(ann["cx"]) - x0, float(ann["cy"]) - y0)
        det, info = find_and_suppress(crop, click, C.v9_config(mpp))
        tally.flag("V9 crop definition unchanged", (x0, y0, click) == (pre["x0"], pre["y0"], pre["click_in_crop"]), f"{fn}: crop {(x0, y0, click)} != {(pre['x0'], pre['y0'], pre['click_in_crop'])}")
        tally.flag("V9 detections identical", list(det.columns) == list(pre["detections"].columns), f"{fn}: columns differ")
        tally.add("V9 detections identical", C.compare_frames(det, pre["detections"], list(pre["detections"].columns), label=f"{fn} V9"), runs=1)
        tally.add("V9 info identical (non-t_*)", C.compare_info(info, pre["info"], label=f"{fn} V9 info"), runs=1)
        cases.append(dict(file_name=fn, x0=x0, y0=y0, click_in_crop=click, n_detections=len(det), n_augmentations=info["n_augmentations"]))
    log["V9"] = cases


def v1_contested(tally, R, log):
    """
        V1's contested-tier path: with an ROI's unanimous mitotic rows removed from a copy of the table, draw_seed_rows must
        equal rr.select_clicks on that table, with flagged True.

        tally (Tally): where to record.
        R (Rois): the ROI loader.
        log (dict): per-case records.

        Returns None.
    """
    rec = dict(tried=[])
    for fn in CONTESTED_CANDIDATES:
        roi = R.rois[R.rois["file_name"] == fn].iloc[0]
        a = R.annotations
        drop = (a["file_name"] == fn) & (a["category_id"] == ds.MITOTIC) & (a["n_mitotic_votes"] == a["n_votes"])
        modified = a[~drop].copy()
        try:
            rgb, mpp, gt, pool, flagged, n_agreement, n_joint_valid, draws = rr.select_clicks(roi, modified)
        except ValueError as exc:
            rec["tried"].append(dict(file_name=fn, error=str(exc)))
            continue
        gt_mod = ds.image_annotations(modified, fn)
        mine = af.draw_seed_rows(rgb, gt_mod, int(roi["image_id"]))
        s0 = dict(draws=draws, n_agreement=n_agreement, pool=pool, n_joint_valid=n_joint_valid, flagged=flagged)
        sub = Tally()
        for (frame, di), (seed_index, ann, specs, refusals, n_seed_pool, _) in zip(mine, draws):
            for field, m, r in [("seed_index", di["seed_index"], seed_index), ("seed_ann_id", di["seed_ann_id"], int(ann["ann_id"])), ("refused_draws", di["refused_draws"], ";".join(refusals)), ("n_retries", di["n_retries"], len(refusals)), ("n_seed_pool", di["n_seed_pool"], n_seed_pool), ("n_agreement_pool", di["n_agreement_pool"], n_agreement), ("n_after_border", di["n_after_border"], len(pool)), ("n_joint_valid", di["n_joint_valid"], n_joint_valid), ("flagged", di["contested_seed_tier"], bool(flagged))]:
                sub.flag("V1 contested tier vs select_clicks", C.same_value(m, r), f"{fn} s{seed_index} {field}: {m!r} != {r!r}")
            sub.flag("V1 contested tier vs select_clicks", checks.same_specs(di["specs"], specs), f"{fn} s{seed_index} specs")
            sub.runs("V1 contested tier vs select_clicks")
        sub.flag("V1 contested tier: flagged True", bool(flagged) is True and all(di["contested_seed_tier"] is True for _, di in mine), f"{fn}: flagged {flagged}")
        for name, e in sub.d.items():
            t = tally._entry(name)
            for k in ("n_runs", "n_rows", "n_values", "n_diverged"):
                t[k] += e[k]
            t["detail"].extend(e["detail"])
        rec.update(file_name=fn, n_removed_unanimous=int(drop.sum()), n_agreement_pool=n_agreement, n_after_border=len(pool), n_joint_valid=n_joint_valid, flagged=bool(flagged), clicks=[di["seed_ann_id"] for _, di in mine], refused=[di["refused_draws"] for _, di in mine])
        del rgb
        break
    tally.flag("V1 contested tier: an ROI with >= 3 joint-valid contested clicks found", "file_name" in rec, f"none of {CONTESTED_CANDIDATES} qualified")
    log["V1_contested"] = rec


def main():
    """
        Run the named-input checks for a tier and write out/<tier>/named.json.

        Returns None.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--tier", required=True, choices=("smoke", "full"))
    parser.add_argument("--only", default="", help="comma-separated subset of V1c,V5,V6,V7,V8,V9")
    args = parser.parse_args()
    only = {s.strip() for s in args.only.split(",") if s.strip()}
    smoke = args.tier == "smoke"
    R = Rois()
    tally, log = Tally(), {}
    t0 = time.time()
    steps = [("V5", lambda: v5(tally, R, log)), ("V6", lambda: v6(tally, R, log)), ("V7", lambda: v7(tally, R, log)), ("V8", lambda: v8(tally, R, log, V8_ALL_INT[:1] if smoke else V8_ALL_INT)), ("V9", lambda: v9(tally, R, log, 1 if smoke else C.V9_N_ROIS)), ("V1c", lambda: v1_contested(tally, R, log))]
    for key, step in steps:
        if only and key not in only:
            continue
        t = time.time()
        step()
        print(f"{key} done [{time.time() - t:.0f}s]", flush=True)
    out = dict(tier=args.tier, checks=tally.d, cases=log, t_total_s=round(time.time() - t0, 1))
    suffix = "" if not only else "_" + "_".join(sorted(only))
    path = C.OUT / args.tier / f"named{suffix}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1, default=str))
    for name in sorted(tally.d):
        e = tally.d[name]
        print(f"  {name}: runs={e['n_runs']} rows={e['n_rows']} values={e['n_values']} diverged={e['n_diverged']}")
        for line in e["detail"][:5]:
            print(f"      {line}")
    print(f"-> {path}")
    if tally.n_diverged:
        raise SystemExit(f"stopped: {tally.n_diverged} divergences")


if __name__ == "__main__":
    main()
