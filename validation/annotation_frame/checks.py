"""
    Per-ROI validation of the annotation frame (spec section 7.2): V1 (the draw), V2's re-runs of production, V3 (the adapter
    against B1 and R1/R2), V4 (rank-key independence, B2 and R4), V8's per-frame schema checks and V10 (evaluation), on one
    ROI against its Step 0 capture and the tracked references. check_roi records every comparison in a Tally; the caller
    writes the ROI's results and stops on any divergence.
"""

from __future__ import annotations

import time
from datetime import datetime

import numpy as np
import pandas as pd

import common as C
from common import ch, ds, ev, prod, rr, ss
from midog_utils import annotation_frame as af

# Spec section 2, written out here rather than read from the module under test.
FRAME_COLUMNS = ("file_name", "x_tl", "y_tl", "width", "height", "class", "object_id", "made_by", "refined", "current_time", "score", "od")
FRAME_DTYPES = {"file_name": "object", "x_tl": "float64", "y_tl": "float64", "width": "float64", "height": "float64", "class": "object", "object_id": "int64", "made_by": "object", "refined": "bool", "current_time": "object", "score": "float64", "od": "float64"}
ROW_KEY = ["file_name", "class", "object_id", "made_by", "refined"]
COMPARED_FRAME_COLUMNS = [c for c in FRAME_COLUMNS if c != "current_time"]

# Spec section 5: the keys run_roi adds to production's info.
ADDED_INFO_KEYS = ("file_name", "mpp", "match_radius_px", "refinement_channel", "refinement_requested", "channel", "tm_method", "peak_min_distance", "deep_floor_z", "max_peaks", "patch_size", "base_size", "seed_object_id", "click_center", "refined", "refusal_reason", "template_center", "offset_px", "t_refine_s", "t_total_s")

R2_SKIPPED = ("angle", "flip", "scale")
R4_TOP_COLS = ["rank", "cx", "cy", "score", "od", "bucket", "matched_ann_id", "matched_category"]
R4_RUN_COLS = ["tp_at_10", "tp_at_20", "tp_at_30", "n_peaks", "n_detections", "n_blanked_px"]
DRAW_COLS = ["seed_ann_id", "n_agreement_pool", "n_after_border", "n_joint_valid", "n_seed_pool", "n_retries", "refused_draws", "contested_seed_tier"]


class Tally:
    """
        Per-check counts: runs, rows, values compared and values diverged, plus the first divergence details.
    """

    def __init__(self):
        self.d = {}

    def _entry(self, name):
        return self.d.setdefault(name, dict(n_runs=0, n_rows=0, n_values=0, n_diverged=0, detail=[]))

    def add(self, name, res, runs=0, rows=None):
        """
            Fold a compare_frames / compare_info result into a check's counts.

            name (str): the check.
            res (dict): the comparison result.
            runs (int): runs this comparison covers.
            rows (int or None): rows it covers; res['n_rows'] when None.

            Returns None.
        """
        e = self._entry(name)
        e["n_runs"] += runs
        e["n_rows"] += res.get("n_rows", 0) if rows is None else rows
        e["n_values"] += res["n_values"]
        e["n_diverged"] += res["n_diverged"]
        e["detail"].extend(res["detail"][:max(0, 20 - len(e["detail"]))])

    def flag(self, name, ok, what, runs=0, rows=0):
        """
            Record one boolean comparison as one value.

            name (str): the check.
            ok (bool): whether it held.
            what (str): the detail line if it did not.
            runs (int): runs it covers.
            rows (int): rows it covers.

            Returns None.
        """
        e = self._entry(name)
        e["n_runs"] += runs
        e["n_rows"] += rows
        e["n_values"] += 1
        if not ok:
            e["n_diverged"] += 1
            if len(e["detail"]) < 20:
                e["detail"].append(what)

    def runs(self, name, n=1):
        """
            Count runs for a check without adding values.

            name (str): the check.
            n (int): runs to add.

            Returns None.
        """
        self._entry(name)["n_runs"] += n

    @property
    def n_diverged(self):
        return sum(e["n_diverged"] for e in self.d.values())


def merge_tallies(tallies):
    """
        Sum several Tally dicts check by check.

        tallies (list[dict]): Tally.d dicts.

        Returns dict: the summed counts per check.
    """
    out = {}
    for t in tallies:
        for name, e in t.items():
            o = out.setdefault(name, dict(n_runs=0, n_rows=0, n_values=0, n_diverged=0, detail=[]))
            for k in ("n_runs", "n_rows", "n_values", "n_diverged"):
                o[k] += e[k]
            o["detail"].extend(e["detail"][:max(0, 20 - len(o["detail"]))])
    return out


def convert_detections(detections, file_name):
    """
        Production detections converted to frame prediction rows by spec section 2, written independently of the module.

        detections (pd.DataFrame): run_production_pipeline's output.
        file_name (str): the ROI's file name.

        Returns pd.DataFrame: prediction rows in frame columns (current_time left out; od NaN when production computed none).
    """
    n = len(detections)
    cx, cy = detections["cx"].to_numpy(dtype=np.float64), detections["cy"].to_numpy(dtype=np.float64)
    od = detections["od"].to_numpy(dtype=np.float64) if "od" in detections else np.full(n, np.nan)
    return pd.DataFrame({"file_name": pd.Series([file_name] * n, dtype=object), "x_tl": cx - 25.0, "y_tl": cy - 25.0, "width": np.full(n, 50.0), "height": np.full(n, 50.0), "class": pd.Series(["PREDICTION"] * n, dtype=object), "object_id": detections["rank"].to_numpy(dtype=np.int64), "made_by": pd.Series(["ALGORITHM"] * n, dtype=object), "refined": np.zeros(n, dtype=bool), "score": detections["score"].to_numpy().astype(np.float64), "od": od})


def expected_human_row(ann, file_name):
    """
        The human-seed row spec section 2 prescribes for a load_annotations row, written independently of the module.

        ann (pd.Series): the annotation.
        file_name (str): the ROI's file name.

        Returns dict: column -> value (current_time left out).
    """
    cx, cy, w, h = float(ann["cx"]), float(ann["cy"]), float(ann["w"]), float(ann["h"])
    return {"file_name": file_name, "x_tl": cx - w / 2, "y_tl": cy - h / 2, "width": w, "height": h, "class": "MITOTIC_GT", "object_id": int(ann["ann_id"]), "made_by": "HUMAN", "refined": False, "score": float("nan"), "od": float("nan")}


def row_matches(row, expected):
    """
        Whether a frame row equals an expected dict exactly, NaN equal to NaN.

        row (pd.Series): the frame row.
        expected (dict): column -> expected value.

        Returns list[str]: the columns that differ.
    """
    bad = []
    for c, v in expected.items():
        a = row[c]
        if isinstance(v, float) and np.isnan(v):
            if not (isinstance(a, (float, np.floating)) and np.isnan(a)):
                bad.append(c)
        elif not (a == v):
            bad.append(c)
    return bad


def schema_checks(tally, name, frame, info, seed_before, seed_after, expected_refined, label):
    """
        V8's per-frame checks, plus V3's row count: columns and dtypes, row key, prediction ids, stamps, the caller's frame
        unmutated, and the row layout (input row, refined seed, predictions).

        tally (Tally): where to record.
        name (str): the check name.
        frame (pd.DataFrame): run_roi's output.
        info (dict): run_roi's info.
        seed_before (pd.DataFrame): a deep copy of the input frame taken before the call.
        seed_after (pd.DataFrame): the caller's input frame after the call.
        expected_refined (bool): the refined flag this run must carry.
        label (str): identifies the run.

        Returns None.
    """
    tally.flag(name, tuple(frame.columns) == FRAME_COLUMNS, f"{label}: columns {list(frame.columns)}")
    tally.flag(name, {c: str(frame[c].dtype) for c in FRAME_COLUMNS if c in frame} == FRAME_DTYPES, f"{label}: dtypes {frame.dtypes.to_dict()}")
    tally.flag(name, not frame.duplicated(ROW_KEY).any(), f"{label}: duplicated row key")
    preds = frame[frame["class"] == "PREDICTION"]
    n_in = len(seed_before)
    tally.flag(name, np.array_equal(preds["object_id"].to_numpy(), np.arange(len(preds))), f"{label}: prediction object_ids not 0..n-1")
    t_refine = datetime.fromisoformat(frame["current_time"].iat[n_in])
    t_pred = [datetime.fromisoformat(t) for t in preds["current_time"].unique()]
    tally.flag(name, len(t_pred) <= 1 and all(t_refine <= t for t in t_pred), f"{label}: refinement stamp {t_refine} later than prediction stamp {t_pred}")
    tally.flag(name, seed_before.equals(seed_after) and seed_before.dtypes.equals(seed_after.dtypes) and seed_before.index.equals(seed_after.index), f"{label}: caller's frame mutated")
    layout = (len(frame) == n_in + 1 + info["n_detections"]) and (frame["made_by"].iloc[:n_in] == "HUMAN").all() and frame["made_by"].iat[n_in] == "ALGORITHM" and frame["class"].iat[n_in] == "MITOTIC_GT" and (frame["class"].iloc[n_in + 1:] == "PREDICTION").all()
    tally.flag(name, bool(layout), f"{label}: row layout / count {len(frame)} != {n_in} + 1 + {info['n_detections']}")
    copied = frame.iloc[:n_in].reset_index(drop=True)
    want = seed_before.astype({c: "float64" for c in ("x_tl", "y_tl", "width", "height")}).reset_index(drop=True)[list(FRAME_COLUMNS)]
    tally.flag(name, copied.equals(want), f"{label}: input rows not copied through unchanged")
    ref_row = frame.iloc[n_in]
    tally.flag(name, bool(ref_row["refined"] == expected_refined) and type(ref_row["refined"]) in (bool, np.bool_) and int(ref_row["object_id"]) == int(seed_before["object_id"].iat[0]), f"{label}: refined row flag/id {ref_row['refined']!r}/{ref_row['object_id']!r}")


def refined_row_check(tally, frame, base_size, tpl_cx, tpl_cy, label):
    """
        V3/V6: the refined row's box equals (int(round(tx)) - b // 2, int(round(ty)) - b // 2, b - 1, b - 1).

        tally (Tally): where to record.
        frame (pd.DataFrame): run_roi's output.
        base_size (int): the reference base_size.
        tpl_cx (float): the reference template centre x.
        tpl_cy (float): the reference template centre y.
        label (str): identifies the run.

        Returns bool: whether round-half-even and round-half-up differ on either axis for this template centre.
    """
    b = int(base_size)
    want = (float(int(round(tpl_cx)) - b // 2), float(int(round(tpl_cy)) - b // 2), float(b - 1), float(b - 1))
    row = frame[(frame["made_by"] == "ALGORITHM") & (frame["class"] == "MITOTIC_GT")]
    got = tuple(float(v) for v in row[["x_tl", "y_tl", "width", "height"]].iloc[0]) if len(row) == 1 else None
    for i, c in enumerate(("x_tl", "y_tl", "width", "height")):
        tally.flag("V3 refined row vs R1", got is not None and got[i] == want[i] and type(row[c].iat[0]) is np.float64, f"{label}: refined {c} {None if got is None else got[i]!r} != {want[i]!r}", rows=1 if i == 0 else 0)
    return any(int(round(v)) != int(np.floor(v + 0.5)) for v in (tpl_cx, tpl_cy))


def r1_row_from_frame(roi, draw_info, condition, info, scored, gt_out, metrics, hem):
    """
        One R1 row recomputed from the frame, info, draw_info and evaluate, with the driver's diagnostic formulas.

        roi (pd.Series): the ROI's build_rois row.
        draw_info (dict): the seed's draw_info.
        condition (str): the run's condition name.
        info (dict): run_roi's info.
        scored (pd.DataFrame): evaluate's candidate rows.
        gt_out (pd.DataFrame): evaluate's gt_out.
        metrics (dict): evaluate's metrics.
        hem (np.ndarray): the ROI's hematoxylin channel, for tpl_sha1.

        Returns dict: the row, t_pipeline_s left out.
    """
    click = info["click_center"]
    tx, ty = info["template_center"]
    cx = (scored["x_tl"] + scored["width"] / 2).to_numpy()
    cy = (scored["y_tl"] + scored["height"] / 2).to_numpy()
    d_click, d_tpl = np.hypot(cx - click[0], cy - click[1]), np.hypot(cx - tx, cy - ty)
    near = np.flatnonzero(d_click <= info["match_radius_px"])
    near_detail = ";".join(f"rank{int(scored['object_id'].iat[i])}:{scored['bucket'].iat[i]}:ann{int(scored['matched_ann_id'].iat[i])}:d_click={d_click[i]:.1f}:d_tpl={d_tpl[i]:.1f}" for i in near)
    od = scored["od"].to_numpy()
    nearest_other = float(np.hypot(gt_out["cx"] - click[0], gt_out["cy"] - click[1]).min())
    row = dict(subset=roi["subset"], file_name=info["file_name"], domain=roi["domain"], image_id=int(roi["image_id"]), seed_index=draw_info["seed_index"], condition=condition, seed_ann_id=draw_info["seed_ann_id"], click_cx=click[0], click_cy=click[1], base_size=info["base_size"], tpl_cx=tx, tpl_cy=ty, tpl_offset_px=round(info["offset_px"], 3), tpl_sha1=rr.template_digest(hem, info["base_size"], *info["template_center"]), n_agreement_pool=draw_info["n_agreement_pool"], n_after_border=draw_info["n_after_border"], n_joint_valid=draw_info["n_joint_valid"], n_seed_pool=draw_info["n_seed_pool"], n_retries=draw_info["n_retries"], refused_draws=draw_info["refused_draws"], contested_seed_tier=draw_info["contested_seed_tier"], n_gt_mitotic=metrics["n_gt_mitotic_eval"], nearest_other_ann_px=round(nearest_other, 2), n_peaks=info["n_peaks"], max_peaks_binding=info["max_peaks_binding"], n_detections=info["n_detections"], n_blanked_px=info["n_blanked_px"], n_near_click=len(near), near_click_detail=near_detail, od_nan=int(np.isnan(od).sum()), od_descending=bool(np.all(np.diff(od) <= 0)), match_radius_px=info["match_radius_px"], mpp=info["mpp"])
    row.update({f"tp_at_{k}": metrics[f"tp_at_{k}"] for k in rr.BUDGETS})
    return row


def top_from_scored(scored, n=C.TOP_N):
    """
        The first n candidates of evaluate's output in the driver's top-30 column names (rank, cx, cy, score, od, buckets).

        scored (pd.DataFrame): evaluate's candidate rows, sorted by object_id.
        n (int): how many.

        Returns pd.DataFrame: the top-n frame.
    """
    s = scored.head(n)
    return pd.DataFrame({"rank": s["object_id"].to_numpy(), "cx": (s["x_tl"] + s["width"] / 2).to_numpy(), "cy": (s["y_tl"] + s["height"] / 2).to_numpy(), "score": s["score"].to_numpy(), "od": s["od"].to_numpy(), "bucket": s["bucket"].to_numpy(), "matched_ann_id": s["matched_ann_id"].to_numpy(), "matched_category": s["matched_category"].to_numpy()})


def driver_seed(ann, specs, refusals, flagged, n_agreement, pool, condition):
    """
        The ss.Seed rr.run_roi builds for one (click, condition) (rerun_bbox3way_postD11.py:233-241).

        ann (pd.Series): the accepted annotation.
        specs (dict): condition -> spec.
        refusals (list[str]): the draw's refusals.
        flagged (bool): contested tier.
        n_agreement (int): agreement pool size.
        pool (pd.DataFrame): border-filtered pool.
        condition (str): the condition.

        Returns ss.Seed: the seed.
    """
    click_xy = (float(ann["cx"]), float(ann["cy"]))
    base_size, tx, ty = specs[condition]
    offset = float(np.hypot(tx - click_xy[0], ty - click_xy[1]))
    return ss.Seed(ann_id=int(ann["ann_id"]), click_xy=click_xy, template_xy=(tx, ty), base_size=base_size, recentred=rr.CONDITIONS[condition] is not None, offset_px=offset, n_retries=len(refusals), agreement_flagged=bool(flagged), n_agreement_pool=n_agreement, n_after_border=len(pool))


def same_specs(a, b):
    """
        Two condition -> spec dicts equal in key order and in every value's type and bits.

        a (dict): first.
        b (dict): second.

        Returns bool: True when identical.
    """
    return list(a) == list(b) and all((a[k] is None and b[k] is None) or (a[k] is not None and b[k] is not None and C.same_value(tuple(a[k]), tuple(b[k]))) for k in a)


def check_draws(tally, roi, s0, draws3, draws6, r4_clicks):
    """
        V1's draw checks on one ROI: draw_seed_rows against the captured rr.select_clicks, and seeds 0-5 against R4's clicks CSV.

        tally (Tally): where to record.
        roi (pd.Series): the ROI's build_rois row.
        s0 (dict): the ROI's Step 0 capture.
        draws3 (list): draw_seed_rows with the default seeds.
        draws6 (list): draw_seed_rows with seed_indices=range(6).
        r4_clicks (pd.DataFrame): R4's clicks CSV.

        Returns None.
    """
    fn = roi["file_name"]
    name = "V1 draw vs select_clicks"
    for (frame, di), (seed_index, ann, specs, refusals, n_seed_pool, _) in zip(draws3, s0["draws"]):
        lbl = f"{fn} s{seed_index}"
        pairs = [("seed_index", di["seed_index"], seed_index), ("seed_ann_id", di["seed_ann_id"], int(ann["ann_id"])), ("refused_draws", di["refused_draws"], ";".join(refusals)), ("n_retries", di["n_retries"], len(refusals)), ("n_seed_pool", di["n_seed_pool"], n_seed_pool), ("n_agreement_pool", di["n_agreement_pool"], s0["n_agreement"]), ("n_after_border", di["n_after_border"], len(s0["pool"])), ("n_joint_valid", di["n_joint_valid"], s0["n_joint_valid"]), ("flagged", di["contested_seed_tier"], bool(s0["flagged"]))]
        for field, mine, ref in pairs:
            tally.flag(name, C.same_value(mine, ref), f"{lbl} {field}: {mine!r} != {ref!r}")
        tally.flag(name, same_specs(di["specs"], specs), f"{lbl} specs: {di['specs']!r} != {specs!r}")
        tally.runs(name)
        tally.d[name]["n_rows"] += 1
        bad = row_matches(frame.iloc[0], expected_human_row(ann, fn))
        tally.flag("V1 human-seed row", len(frame) == 1 and not bad and tuple(frame.columns) == FRAME_COLUMNS, f"{lbl}: human row differs on {bad}", runs=1, rows=1)
    tally.flag(name, len(draws3) == len(s0["draws"]), f"{fn}: {len(draws3)} draws != {len(s0['draws'])}")

    mine6 = pd.DataFrame([{**{c: di[c] for c in DRAW_COLS}, "seed_index": di["seed_index"], "click_cx": float(f["x_tl"].iat[0] + f["width"].iat[0] / 2), "click_cy": float(f["y_tl"].iat[0] + f["height"].iat[0] / 2)} for f, di in draws6])
    ref6 = r4_clicks[r4_clicks["file_name"] == roi["file_name"]].sort_values("seed_index", kind="stable").reset_index(drop=True)
    tally.add("V1 seeds 0-5 vs R4 clicks", C.compare_frames(mine6, ref6, ["seed_index"] + DRAW_COLS, label=f"{fn} R4 clicks"), runs=len(mine6))
    tally.add("V1 extra: seeds 0-5 click centre vs R4 clicks", C.compare_frames(mine6, ref6, ["click_cx", "click_cy"], label=f"{fn} R4 click centre"), runs=len(mine6))
    same3 = all(same_specs(a[1]["specs"], b[1]["specs"]) and all(C.same_value(a[1][k], b[1][k]) for k in DRAW_COLS + ["seed_index"]) and a[0].drop(columns="current_time").equals(b[0].drop(columns="current_time")) for a, b in zip(draws3, draws6[:3]))
    tally.flag("V1 extra: range(6) prefix equals default draw", same3, f"{fn}: seeds 0-2 of range(6) differ from the default draw")


def check_roi(roi, annotations, refs, s0):
    """
        Run V1-V4, V8's per-frame checks and V10 on one ROI.

        roi (pd.Series): the ROI's build_rois row.
        annotations (pd.DataFrame): the full annotation table.
        refs (dict): the reference frames (r1, r4_per_run, r4_top, r4_clicks).
        s0 (dict): the ROI's Step 0 capture.

        Returns dict: the tally and the recomputed R1 rows, R2 frames and selection rows.
    """
    t_start = time.time()
    fn, image_id = roi["file_name"], int(roi["image_id"])
    tally = Tally()
    rgb = ds.load_roi(roi["path"])
    mpp = ds.roi_mpp(roi["path"])
    tally.flag("setup: mpp and shape match Step 0", C.same_value(mpp, s0["mpp"]) and tuple(rgb.shape) == tuple(s0["shape"]), f"{fn}: mpp/shape differ from Step 0")
    gt = ds.image_annotations(annotations, fn)
    tally.flag("setup: gt matches Step 0", gt.equals(s0["gt"]), f"{fn}: gt differs from Step 0")

    draws3 = af.draw_seed_rows(rgb, gt, image_id)
    draws6 = af.draw_seed_rows(rgb, gt, image_id, seed_indices=range(6))
    check_draws(tally, roi, s0, draws3, draws6, refs["r4_clicks"])
    hem = ch.to_channel(rgb, "hematoxylin_od")

    r1_rows, r2_frames, sel_rows, r4_tops, n_half_differ = [], [], [], [], 0
    area = rgb.shape[0] * rgb.shape[1] * mpp * mpp / 1e6
    r1 = refs["r1"]
    for (seed_frame, di), (seed_index, ann, specs, refusals, n_seed_pool, _) in zip(draws3, s0["draws"]):
        seed_ann_id = di["seed_ann_id"]
        gt_eval = gt[gt["ann_id"] != seed_ann_id].reset_index(drop=True)
        for condition, channel in C.CONDITION_CHANNEL.items():
            lbl = f"{fn} s{seed_index} {condition}"
            b1_det, b1_info = s0["runs"][(seed_index, condition, "chromatin_od")]
            ref_row = r1[(r1["file_name"] == fn) & (r1["seed_index"] == seed_index) & (r1["condition"] == condition)]
            assert len(ref_row) == 1, lbl
            ref_row = ref_row.iloc[0]

            out = {}
            for rank_key in ("chromatin_od", "tm_score"):
                before = seed_frame.copy(deep=True)
                frame, info = af.run_roi(rgb, seed_frame, mpp, rank_key=rank_key, refinement_channel=channel)
                schema_checks(tally, "V8 per-frame schema", frame, info, before, seed_frame, channel is not None, f"{lbl} {rank_key}")
                tally.runs("V8 per-frame schema")
                scored, gt_out, metrics = af.evaluate(frame, gt, mpp, rgb.shape)
                out[rank_key] = (frame, info, scored, gt_out, metrics)

            frame_c, info_c, scored_c, gt_out_c, metrics_c = out["chromatin_od"]
            frame_t, info_t, scored_t, gt_out_t, metrics_t = out["tm_score"]

            # V3: predictions vs B1 converted; production info vs B1 info; refined row vs R1; R1/R2 rows.
            preds_c = frame_c[frame_c["class"] == "PREDICTION"]
            tally.add("V3 predictions vs B1", C.compare_frames(preds_c, convert_detections(b1_det, fn), COMPARED_FRAME_COLUMNS, label=f"{lbl} V3"), runs=1)
            prod_part = {k: v for k, v in info_c.items() if k not in ADDED_INFO_KEYS}
            tally.add("V3 production info vs B1", C.compare_info(prod_part, b1_info, label=f"{lbl} V3 info"), runs=1)
            tally.flag("V3 row count 2 + n_detections", len(frame_c) == 2 + b1_info["n_detections"], f"{lbl}: {len(frame_c)} rows", runs=1)
            n_half_differ += refined_row_check(tally, frame_c, ref_row["base_size"], ref_row["tpl_cx"], ref_row["tpl_cy"], lbl)
            tally.runs("V3 refined row vs R1")
            added_ok = [info_c["refinement_channel"] == channel, info_c["refinement_requested"] == (channel is not None), info_c["refined"] == (channel is not None), info_c["refusal_reason"] is None, C.same_value(info_c["match_radius_px"], ev.radius_px(mpp)), info_c["channel"] == prod.CHANNEL, info_c["tm_method"] == prod.TM_METHOD, info_c["max_peaks"] == prod.MAX_PEAKS, info_c["seed_object_id"] == seed_ann_id, C.same_value(tuple(info_c["template_center"]), tuple(specs[condition][1:])), info_c["base_size"] == specs[condition][0], C.same_value(info_c["click_center"], (float(ann["cx"]), float(ann["cy"])))]
            tally.flag("V3 extra: run_roi's added info", all(added_ok) and all(k in info_c for k in ADDED_INFO_KEYS), f"{lbl}: added info {added_ok}", runs=1)
            r1_rows.append(r1_row_from_frame(roi, di, condition, info_c, scored_c, gt_out_c, metrics_c, hem))
            r2_frames.append(top_from_scored(scored_c).assign(subset=roi["subset"], file_name=fn, seed_index=seed_index, condition=condition))
            sel_rows.append(dict(file_name=fn, seed_index=seed_index, condition=condition, seed_ann_id=di["seed_ann_id"], base_size=info_c["base_size"], tpl_cx=info_c["template_center"][0], tpl_cy=info_c["template_center"][1], tpl_offset_px=round(info_c["offset_px"], 3), tpl_sha1=rr.template_digest(hem, info_c["base_size"], *info_c["template_center"]), n_agreement_pool=di["n_agreement_pool"], n_after_border=di["n_after_border"], n_joint_valid=di["n_joint_valid"], n_seed_pool=di["n_seed_pool"], n_retries=di["n_retries"], refused_draws=di["refused_draws"]))

            # V10: evaluate vs evaluate_run on B1's detections directly; tp_at_k vs R1.
            _, _, direct = ev.evaluate_run(b1_det, gt_eval, ev.radius_px(mpp), area, rgb.shape)
            tally.add("V10 metrics vs evaluate_run on B1", C.compare_info({k: metrics_c[k] for k in direct if k in metrics_c}, direct, label=f"{lbl} V10"), runs=1)
            for k in rr.BUDGETS:
                tally.flag("V10 tp_at_k vs R1", C.same_value(metrics_c[f"tp_at_{k}"], int(ref_row[f"tp_at_{k}"])), f"{lbl}: tp_at_{k} {metrics_c[f'tp_at_{k}']} != {ref_row[f'tp_at_{k}']}")
            tally.runs("V10 tp_at_k vs R1")
            tally.flag("V10 extra: seed_object_id", metrics_c["seed_object_id"] == seed_ann_id, f"{lbl}: metrics seed_object_id {metrics_c['seed_object_id']}")

            # V4: same coordinate set and identical score/od at each coordinate under both rank keys.
            pc = frame_c[frame_c["class"] == "PREDICTION"][["x_tl", "y_tl", "score", "od"]]
            pt = frame_t[frame_t["class"] == "PREDICTION"][["x_tl", "y_tl", "score", "od"]]
            m = pc.merge(pt, on=["x_tl", "y_tl"], how="outer", suffixes=("_c", "_t"), indicator=True)
            same_set = len(pc) == len(pt) == len(m) and (m["_merge"] == "both").all() and not pc.duplicated(["x_tl", "y_tl"]).any()
            tally.flag("V4 coordinate set, chromatin_od vs tm_score", bool(same_set), f"{lbl}: coordinate sets differ ({len(pc)} vs {len(pt)}, merged {len(m)})", runs=1, rows=len(pc))
            both = m[m["_merge"] == "both"]
            tally.add("V4 score/od at each coordinate", C.compare_frames(both[["score_t", "od_t"]].set_axis(["score", "od"], axis=1), both[["score_c", "od_c"]].set_axis(["score", "od"], axis=1), ["score", "od"], label=f"{lbl} V4"), runs=1)
            tally.flag("V4 extra: tm_score info matches chromatin_od info apart from ranking", all(C.same_value(info_t[k], info_c[k]) for k in info_c if not k.startswith("t_") and k not in ("rank_key",)), f"{lbl}: tm_score info differs beyond rank_key")

            if roi["subset"] == "original_14":
                b2_det, b2_info = s0["runs"][(seed_index, condition, "tm_score")]
                preds_t = frame_t[frame_t["class"] == "PREDICTION"]
                tally.add("V4 tm_score predictions vs B2", C.compare_frames(preds_t, convert_detections(b2_det, fn), [c for c in COMPARED_FRAME_COLUMNS if c != "od"], label=f"{lbl} V4 B2"), runs=1)

                # V2: production re-run, default always_od and always_od=True, against B2.
                seed = driver_seed(ann, specs, refusals, s0["flagged"], s0["n_agreement"], s0["pool"], condition)
                d0, i0 = prod.run_production_pipeline(rgb, seed, mpp, rank_key="tm_score")
                tally.flag("V2 default always_od vs B2", list(d0.columns) == list(b2_det.columns), f"{lbl}: columns {list(d0.columns)}")
                tally.add("V2 default always_od vs B2", C.compare_frames(d0, b2_det, list(b2_det.columns), label=f"{lbl} V2 default"), runs=1)
                tally.add("V2 default always_od vs B2 (info)", C.compare_info(i0, b2_info, label=f"{lbl} V2 default info"), runs=1)
                d1, i1 = prod.run_production_pipeline(rgb, seed, mpp, rank_key="tm_score", always_od=True)
                tally.flag("V2 always_od=True vs B2", list(d1.columns) == list(b2_det.columns) + ["od"], f"{lbl}: columns {list(d1.columns)}")
                tally.add("V2 always_od=True vs B2", C.compare_frames(d1, b2_det, list(b2_det.columns), label=f"{lbl} V2 always_od"), runs=1)
                tally.add("V2 always_od=True vs B2 (info)", C.compare_info(i1, b2_info, exclude=("od_computed",), label=f"{lbl} V2 always_od info"), runs=1)
                tally.flag("V2 always_od=True vs B2 (info)", i1["od_computed"] is True and b2_info["od_computed"] is False, f"{lbl}: od_computed {i1['od_computed']!r} / B2 {b2_info['od_computed']!r}")
                mo = d1[["cx", "cy", "od"]].merge(b1_det[["cx", "cy", "od"]], on=["cx", "cy"], how="left", suffixes=("", "_b1"))
                tally.add("V2 extra: always_od od vs B1 od at each coordinate", C.compare_frames(mo[["od"]], mo[["od_b1"]].set_axis(["od"], axis=1), ["od"], label=f"{lbl} V2 od vs B1"), runs=1)
                tally.add("V2 extra: always_od od vs run_roi tm_score od", C.compare_frames(d1[["od"]], preds_t[["od"]], ["od"], label=f"{lbl} V2 od vs run_roi"), runs=1)

            if condition == "default_51":
                r4t = refs["r4_top"]
                ref_top = r4t[(r4t["file_name"] == fn) & (r4t["seed_index"] == seed_index) & (r4t["condition"] == "default_51") & (r4t["rank_key"] == "tm_score")]
                mine_top = top_from_scored(scored_t)
                tally.add("V4 default_51 tm_score top30 vs R4", C.compare_frames(mine_top, ref_top, R4_TOP_COLS, float32_cols=("score",), label=f"{lbl} V4 R4 top"), runs=1)
                r4r = refs["r4_per_run"]
                ref_run = r4r[(r4r["file_name"] == fn) & (r4r["seed_index"] == seed_index) & (r4r["condition"] == "default_51") & (r4r["rank_key"] == "tm_score")]
                mine_run = pd.DataFrame([dict(tp_at_10=metrics_t["tp_at_10"], tp_at_20=metrics_t["tp_at_20"], tp_at_30=metrics_t["tp_at_30"], n_peaks=info_t["n_peaks"], n_detections=info_t["n_detections"], n_blanked_px=info_t["n_blanked_px"])])
                tally.add("V4 default_51 tm_score per-run vs R4", C.compare_frames(mine_run, ref_run, R4_RUN_COLS, label=f"{lbl} V4 R4 run"), runs=1)
                r4_tops.append(mine_top.assign(file_name=fn, seed_index=seed_index))

    del rgb, hem
    return dict(file_name=fn, tally=tally.d, r1_rows=r1_rows, r2_frames=r2_frames, sel_rows=sel_rows, r4_tops=r4_tops, n_half_even_vs_half_up_differ=int(n_half_differ), t_check_s=round(time.time() - t_start, 1))
