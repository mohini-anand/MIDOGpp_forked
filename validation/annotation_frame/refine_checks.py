"""
    Checks for the refined predictions (run_roi's refine_predictions; ANNOTATION_DATAFRAME_SPEC.md §3 "Refined
    predictions"). For each seed and configuration: the flag-on frame starts with the flag-off frame and adds exactly one
    PREDICTION_REFINED row per prediction; the row key stays unique; fallbacks are their prediction's own box; near_edge is
    the set of prediction centres within 36 px of an edge, counted by hand; evaluate is unchanged; and the seed path, run on
    a click at each of the first 10 predictions, gives the identical box. Writes out/refine/refine.json. Usage:
    refine_checks.py
"""

from __future__ import annotations

import json
import time

import common as C
from common import np, pd, tm
from checks import Tally
from midog_utils import annotation_frame as af
from named_checks import Rois

ROIS = ("300.tiff", "013.tiff", "245.tiff")
CHANNELS = ("gray_inverted", "hematoxylin_od", None)
RANK_KEYS = ("tm_score", "chromatin_od")
NEW_INFO_KEYS = ("prediction_refinement", "n_predictions_refined", "t_refine_predictions_s")
ROW_KEY = ["file_name", "class", "object_id", "made_by", "refined"]
VALUE_COLUMNS = [c for c in af.COLUMNS if c != "current_time"]
SEED_PATH_TOP = 10


def near_edge_by_hand(shape, cx, cy):
    """
        Whether each centre is within tm.PATCH_SIZE // 2 px of an edge, from the rounded centre, without seed_selection.

        shape (tuple): the ROI's array shape.
        cx (np.ndarray): centre x.
        cy (np.ndarray): centre y.

        Returns np.ndarray: one bool per centre.
    """
    half, (h, w) = tm.PATCH_SIZE // 2, shape[:2]
    ix, iy = np.rint(cx).astype(int), np.rint(cy).astype(int)
    return (ix < half) | (iy < half) | (ix + half >= w) | (iy + half >= h)


def seeds_for(roi, rgb, gt):
    """
        The seeds checked on one ROI: the validation draw's seed 0, and on 300.tiff a 41 px drawn box.

        roi (pd.Series): the ROI's row.
        rgb (np.ndarray): the ROI's image.
        gt (pd.DataFrame): the ROI's annotations.

        Returns dict: seed name -> one-row input frame.
    """
    fn = roi["file_name"]
    seeds = {"midog0": af.draw_seed_rows(rgb, gt, int(roi["image_id"]), seed_indices=(0,))[0][0]}
    if fn == "300.tiff":
        seeds["box41"] = af.box_frame(fn, 1391, 3531, 41, 41, 7)
    return seeds


def check_run(tally, log, rgb, gt, mpp, seed, channel, rank_key, scorable, label):
    """
        One seed and configuration, run with the flag off and on.

        tally (Tally): where to record.
        log (list): per-run records.
        rgb (np.ndarray): the ROI's image.
        gt (pd.DataFrame): the ROI's annotations.
        mpp (float): the ROI's microns per pixel.
        seed (pd.DataFrame): the one-row input frame.
        channel (str or None): the refinement channel.
        rank_key (str): the rank key.
        scorable (bool): the seed is a MIDOG annotation, so evaluate can score the frame.
        label (str): identifies the run in detail lines.

        Returns tuple: the PREDICTION_REFINED rows, info["prediction_refinement"], the prediction centres (cx, cy) and the near-edge mask.
    """
    f_off, i_off = af.run_roi(rgb, seed, mpp, rank_key=rank_key, refinement_channel=channel)
    f_on, i_on = af.run_roi(rgb, seed, mpp, rank_key=rank_key, refinement_channel=channel, refine_predictions=True)
    n0, n = len(f_off), i_on["n_detections"]
    tally.flag("R1 flag off adds no info keys", not any(k in i_off for k in NEW_INFO_KEYS), f"{label}: flag-off info has {[k for k in NEW_INFO_KEYS if k in i_off]}", runs=1)
    tally.add("R1 flag-on head equals the flag-off frame", C.compare_frames(f_on.iloc[:n0].reset_index(drop=True), f_off, VALUE_COLUMNS, label=label), runs=1)
    tally.add("R1 flag-on info equals flag-off info apart from the new keys", C.compare_info(i_on, i_off, exclude=NEW_INFO_KEYS, label=label), runs=1)

    ref, table = f_on.iloc[n0:].reset_index(drop=True), i_on["prediction_refinement"]
    preds = f_on[f_on["class"] == af.PREDICTION].reset_index(drop=True)
    layout = len(f_on) == n0 + n and len(table) == n and (ref["class"] == af.PREDICTION_REFINED).all() and (ref["made_by"] == af.ALGORITHM).all() and np.array_equal(ref["object_id"].to_numpy(), np.arange(n)) and np.array_equal(table["object_id"].to_numpy(), np.arange(n)) and ref["score"].isna().all() and ref["od"].isna().all()
    tally.flag("R2 one PREDICTION_REFINED row per prediction, object_id = rank, score and od NaN", bool(layout), f"{label}: {len(f_on)} rows, {n0} + {n} expected", runs=1, rows=n)
    tally.flag("R2 refined flags match info", np.array_equal(ref["refined"].to_numpy(), table["refined"].to_numpy()) and i_on["n_predictions_refined"] == int(table["refined"].sum()), f"{label}: refined flags vs info", runs=1)
    tally.flag("R3 row key unique", not f_on.duplicated(ROW_KEY).any(), f"{label}: {int(f_on.duplicated(ROW_KEY).sum())} duplicated row keys", runs=1)

    fallback = ~table["refined"].to_numpy()
    tally.add("R4 a fallback row is its prediction's box", C.compare_frames(ref.loc[fallback, list(af.COORD_COLUMNS)].reset_index(drop=True), preds.loc[fallback, list(af.COORD_COLUMNS)].reset_index(drop=True), list(af.COORD_COLUMNS), label=label), runs=1)
    cx, cy = af.box_center(preds["x_tl"], preds["width"]).to_numpy(), af.box_center(preds["y_tl"], preds["height"]).to_numpy()
    near = near_edge_by_hand(rgb.shape, cx, cy) & (channel is not None)
    tally.flag("R4 near_edge is exactly the predictions within 36 px of an edge", np.array_equal(near, (table["refusal_reason"] == af.NEAR_EDGE).to_numpy()), f"{label}: near_edge set differs", runs=1)
    reasons_ok = table.loc[~fallback, "refusal_reason"].isna().all() and (table.loc[fallback, "refusal_reason"].isin([af.GATE_REFUSED, af.TEMPLATE_UNREADABLE, af.NEAR_EDGE]).all() if channel else table["refusal_reason"].isna().all())
    tally.flag("R4 refusal reasons: none when refined, a defined reason on every gated fallback", bool(reasons_ok), f"{label}: reasons {table['refusal_reason'].value_counts(dropna=False).to_dict()}", runs=1)
    tally.flag("R4 offset_px from the prediction centre", np.array_equal(table["offset_px"].to_numpy(), np.hypot(table["template_cx"].to_numpy() - cx, table["template_cy"].to_numpy() - cy)), f"{label}: offset_px", runs=1)

    if scorable:
        s_off, _, m_off = af.evaluate(f_off, gt, mpp, rgb.shape)
        s_on, _, m_on = af.evaluate(f_on, gt, mpp, rgb.shape)
        tally.flag("R5 evaluate ignores the refined rows", m_off == m_on and s_off.drop(columns="current_time").equals(s_on.drop(columns="current_time")), f"{label}: evaluate differs", runs=1)
    log.append(dict(run=label, n_detections=n, n_predictions_refined=i_on["n_predictions_refined"], reasons={str(k): int(v) for k, v in table["refusal_reason"].value_counts().items()}, t_refine_predictions_s=i_on["t_refine_predictions_s"]))
    return ref, table, cx, cy, near


def check_seed_path(tally, rgb, mpp, fn, channel, ref, table, cx, cy, near, label):
    """
        The seed's own path, run on a click at each of the first predictions, must give the identical refined box.

        tally (Tally): where to record.
        rgb (np.ndarray): the ROI's image.
        mpp (float): the ROI's microns per pixel.
        fn (str): the ROI's file name.
        channel (str): the refinement channel.
        ref (pd.DataFrame): the PREDICTION_REFINED rows.
        table (pd.DataFrame): info["prediction_refinement"].
        cx (np.ndarray): prediction centres, x.
        cy (np.ndarray): prediction centres, y.
        near (np.ndarray): which predictions are near an edge.
        label (str): identifies the run in detail lines.

        Returns None.
    """
    for k in range(min(SEED_PATH_TOP, len(ref))):
        click = af.click_frame(fn, cx[k], cy[k], 1000 + k)
        if near[k]:
            try:
                af.run_roi(rgb, click, mpp, refinement_channel=channel)
                tally.flag("R6 a near-edge prediction raises as a seed", False, f"{label}: prediction {k} did not raise")
            except ValueError:
                tally.flag("R6 a near-edge prediction raises as a seed", True, "")
            continue
        frame, info = af.run_roi(rgb, click, mpp, refinement_channel=channel)
        seed_row, row = frame.iloc[1], ref.iloc[k]
        same = all(float(seed_row[c]) == float(row[c]) for c in af.COORD_COLUMNS) and bool(seed_row["refined"]) == bool(row["refined"]) and info["refusal_reason"] == table["refusal_reason"].iat[k] and info["template_center"] == (table["template_cx"].iat[k], table["template_cy"].iat[k]) and info["base_size"] == table["base_size"].iat[k]
        tally.flag("R6 the seed path at a prediction's centre gives its refined box", same, f"{label}: prediction {k}: seed {info['base_size']} at {info['template_center']} vs {table.iloc[k].to_dict()}", runs=1)


def check_empty(tally, rgb):
    """
        No detections: an empty, schema-valid block and table.

        tally (Tally): where to record.
        rgb (np.ndarray): any ROI image.

        Returns None.
    """
    empty = pd.DataFrame({"rank": np.array([], dtype=np.int64), "cx": np.array([], dtype=np.float64), "cy": np.array([], dtype=np.float64)})
    rows, table = af._refined_prediction_rows(rgb, None, empty, None, "x.tiff")
    ok = len(rows) == 0 and len(table) == 0 and tuple(rows.columns) == af.COLUMNS and all(rows[c].dtype == np.dtype(af.DTYPES[c]) for c in af.COLUMNS)
    tally.flag("R7 no detections gives an empty, schema-valid block", ok, f"{len(rows)} rows, dtypes {rows.dtypes.to_dict()}")


def main():
    """
        Run every refined-prediction check and write out/refine/refine.json.

        Returns None.
    """
    t0 = time.time()
    tally, log = Tally(), []
    R = Rois()
    for fn in ROIS:
        roi, rgb, mpp, gt = R.get(fn)
        for sname, seed in seeds_for(roi, rgb, gt).items():
            for channel in CHANNELS:
                for rank_key in RANK_KEYS:
                    label = f"{fn} {sname} {channel} {rank_key}"
                    ref, table, cx, cy, near = check_run(tally, log, rgb, gt, mpp, seed, channel, rank_key, sname == "midog0", label)
                    if channel is not None and rank_key == "tm_score" and sname == "midog0":
                        check_seed_path(tally, rgb, mpp, fn, channel, ref, table, cx, cy, near, label)
                    print(f"{label} done [{time.time() - t0:.0f}s]", flush=True)
    check_empty(tally, R.get(ROIS[-1])[1])
    out = dict(checks=tally.d, runs=log, t_total_s=round(time.time() - t0, 1))
    path = C.OUT / "refine" / "refine.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1, default=str))
    for name in sorted(tally.d):
        e = tally.d[name]
        print(f"  {name}: runs={e['n_runs']} values={e['n_values']} diverged={e['n_diverged']}")
        for line in e["detail"][:5]:
            print(f"      {line}")
    print(f"-> {path}")
    if tally.n_diverged:
        raise SystemExit(f"stopped: {tally.n_diverged} divergences")


if __name__ == "__main__":
    main()
