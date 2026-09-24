"""
    Checks for the drawn-box seed (DECISIONS.md D12, SINGLE_ENTRY_POINT_PLAN.md V-C): gate_window_side's vectors,
    box_frame's snapping, the gate running on the box's own window, the fallback to the human's box, the 36 px edge, and the
    input-row size check. Expected values come from pixel counts worked out here, never from af.gate_window_side. Writes
    out/box/box.json. Usage: box_checks.py
"""

from __future__ import annotations

import json
import time
from decimal import ROUND_FLOOR, ROUND_HALF_EVEN, Decimal

import common as C
from common import ch, ds, np, pd, ss, tm
from checks import Tally
from midog_utils import annotation_frame as af
from named_checks import Rois

FN = "300.tiff"


def hand_pixels(w):
    """
        A drawn side's stored pixel count, worked out independently: w rounded half to even, raised to odd.

        w (float): the drawn side in pixels.

        Returns int: the odd pixel count.
    """
    n = int(Decimal(repr(float(w))).quantize(Decimal(1), rounding=ROUND_HALF_EVEN))
    return n if n % 2 else n + 1


def hand_centre(x0, w):
    """
        The pixel containing a drawn box's centre, worked out independently.

        x0 (float): the box's left edge, pixel-edge coordinates.
        w (float): the drawn side.

        Returns int: the pixel index.
    """
    return int((Decimal(repr(float(x0))) + Decimal(repr(float(w))) / 2).to_integral_value(rounding=ROUND_FLOOR))


def raises(call):
    """
        Whether a call raises, and what.

        call (callable): the call.

        Returns tuple[bool, str or None]: (raised, "Type: message").
    """
    try:
        call()
        return False, None
    except Exception as exc:  # noqa: BLE001
        return True, f"{type(exc).__name__}: {exc}"


def geometry(row):
    """
        A frame row's box.

        row (pd.Series): one frame row.

        Returns tuple[float, float, float, float]: (x_tl, y_tl, width, height).
    """
    return float(row["x_tl"]), float(row["y_tl"]), float(row["width"]), float(row["height"])


def check_window_vectors(tally, log):
    """
        1. gate_window_side on frame widths (extent = pixel count - 1), with the expected window counted by hand. The window
        comes from the longer side, so a narrow box is valid when its longer side is 18 px or more.

        tally (Tally): where to record.
        log (dict): per-case records.

        Returns None.
    """
    name = "C1 gate_window_side vectors"
    cases = []
    for (w, h), want in [((50.0, 50.0), 51), ((51.0, 51.0), 53), ((71.0, 71.0), 73), ((72.0, 72.0), 73), ((17.0, 17.0), 19), ((18.0, 18.0), 19), ((20.0, 50.0), 51), ((50.0, 20.0), 51), ((50.0, 16.0), 51), ((16.0, 72.0), 73), ((16.0, 17.0), 19), ((40.0, 40.0), 41), ((49.5, 49.5), 51), ((50.5, 50.5), 51), ((51.5, 51.5), 53)]:
        got = af.gate_window_side(w, h)
        tally.flag(name, got == want and type(got) is int, f"({w}, {h}) -> {got!r}, expected {want}")
        cases.append(dict(width=w, height=h, window=got))
    for w, h in [(73.0, 73.0), (16.0, 16.0), (80.0, 50.0), (15.0, 16.0), (float("nan"), 50.0), (50.0, float("nan")), (-1.0, 50.0), (float("inf"), 50.0)]:
        r, err = raises(lambda: af.gate_window_side(w, h))
        tally.flag(name, r and err.startswith("ValueError"), f"({w}, {h}) did not raise ValueError: {err}")
        cases.append(dict(width=w, height=h, error=err))
    log["C1"] = cases


def check_box_frame(tally, log):
    """
        2. box_frame: stored geometry for fixed boxes, the raises, shift behaviour, click_frame equivalence at integer
        origins, and non-square boxes.

        tally (Tally): where to record.
        log (dict): per-case records.

        Returns None.
    """
    cases = []
    for box, want in [((0, 0, 51, 51), (0.0, 0.0, 50.0, 50.0)), ((0, 0, 52, 52), (0.0, 0.0, 52.0, 52.0)), ((0, 0, 40, 40), (0.0, 0.0, 40.0, 40.0)), ((0, 0, 73, 73), (0.0, 0.0, 72.0, 72.0)), ((0, 0, 18, 18), (0.0, 0.0, 18.0, 18.0)), ((10, 20, 30, 60), (10.0, 20.0, 30.0, 60.0)), ((10, 20, 31, 61), (10.0, 20.0, 30.0, 60.0)), ((11.5, 11.5, 51, 51), (12.0, 12.0, 50.0, 50.0))]:
        f = af.box_frame(FN, *box, 7)
        row = f.iloc[0]
        ok_schema = raises(lambda: af._check_schema(f))[0] is False and len(f) == 1
        ok_cols = row["class"] == af.MITOTIC_GT and row["made_by"] == af.HUMAN and not bool(row["refined"]) and int(row["object_id"]) == 7 and np.isnan(row["score"]) and np.isnan(row["od"]) and row["file_name"] == FN
        tally.flag("C2 box_frame geometry", geometry(row) == want, f"{box} -> {geometry(row)}, expected {want}")
        tally.flag("C2 box_frame row schema", ok_schema and ok_cols, f"{box}: schema {ok_schema}, columns {ok_cols}")
        cases.append(dict(box=box, stored=geometry(row)))
    for box in [(0, 0, 74, 74), (0, 0, 17, 17), (0, 0, 0, 0), (0, 0, 0.4, 30), (0, 0, 30, 80), (0, 0, float("nan"), 30), (float("inf"), 0, 30, 30)]:
        r, err = raises(lambda: af.box_frame(FN, *box, 0))
        tally.flag("C2 box_frame raises", r and err.startswith("ValueError"), f"{box} did not raise ValueError: {err}")
        cases.append(dict(box=box, error=err))

    # Shift behaviour: the side depends only on w, the stored box moves with the drawn one pixel for pixel, and the
    # stored centre is the pixel containing the drawn centre, for every origin fraction and both origin parities.
    n_shift = 0
    for w in np.arange(18.0, 73.5, 0.5):
        widths = set()
        for k in (100, 101):
            for frac in (0.0, 0.25, 0.5, 0.75):
                x0 = k + frac
                a = af.box_frame(FN, x0, x0, w, w, 0).iloc[0]
                b = af.box_frame(FN, x0 + 1, x0 + 1, w, w, 0).iloc[0]
                widths.add(float(a["width"]))
                centre = float(a["x_tl"]) + float(a["width"]) / 2
                ok = float(a["width"]) == hand_pixels(w) - 1 and centre == hand_centre(x0, w) and float(b["x_tl"]) == float(a["x_tl"]) + 1 and float(a["y_tl"]) == float(a["x_tl"]) and abs(centre - (x0 + w / 2 - 0.5)) <= 0.5
                tally.flag("C2 box_frame shift behaviour", ok, f"w {w} x0 {x0}: width {a['width']}, x_tl {a['x_tl']}, shifted x_tl {b['x_tl']}")
                n_shift += 1
        tally.flag("C2 box_frame shift behaviour", len(widths) == 1, f"w {w}: widths {sorted(widths)} vary with the origin")

    # Integer origins: a 51 px box is its centre pixel's click, row for row.
    n_eq = 0
    for x0, y0 in [(0, 0), (1390, 3530), (1391, 3531), (6000, 17), (37, 4000)]:
        a = af.box_frame(FN, x0, y0, 51, 51, 3, current_time="t")
        b = af.click_frame(FN, x0 + 25, y0 + 25, 3, current_time="t")
        tally.flag("C2 box_frame == click_frame at integer origins", a.equals(b), f"({x0}, {y0}): {a.iloc[0].to_dict()} != {b.iloc[0].to_dict()}")
        n_eq += 1
    # Not claimed at other origins: the box snaps to the pixel containing its centre, the click keeps its half pixel.
    a, b = af.box_frame(FN, 11.5, 11.5, 51, 51, 0).iloc[0], af.click_frame(FN, 36.5, 36.5, 0).iloc[0]
    log["C2"] = dict(cases=cases, n_shift_cases=n_shift, n_click_equivalence=n_eq, non_integer_origin_example=dict(box_centre=float(a["x_tl"]) + 25, click_centre=float(b["x_tl"]) + 25))


def gate_cases(rgb, gt):
    """
        Boxes on 300.tiff whose own window changes the gate's answer against 51 px: accepted boxes and refused boxes, each
        verified before use.

        rgb (np.ndarray): the ROI.
        gt (pd.DataFrame): its annotations.

        Returns tuple[list[dict], list[dict]]: (accepted cases, refused cases).
    """
    gray, hem = ch.to_channel(rgb, "gray_inverted"), ch.to_channel(rgb, "hematoxylin_od")
    channels = {"gray_inverted": gray, "hematoxylin_od": hem}
    mit = gt[gt["category_id"] == ds.MITOTIC]
    mit = mit[[ss._patch_readable(rgb.shape, cx, cy, tm.PATCH_SIZE) for cx, cy in zip(mit["cx"], mit["cy"])]]
    # Each spec: (label, channel, drawn w, drawn h, shift). The drawn centre is at the annotation's pixel centre plus shift,
    # shift in [-0.5, 0.5), so it stays inside that pixel; integer annotation centres only.
    specs = [("even-pixel box 40 px", "gray_inverted", 40, 40, 0.0), ("box 31 px", "gray_inverted", 31, 31, 0.0), ("box 61 px", "gray_inverted", 61, 61, 0.0), ("non-square 30 x 44 px", "gray_inverted", 30, 44, 0.0), ("half-integer origin, 45 px", "gray_inverted", 45, 45, -0.5), ("hem box 35 px", "hematoxylin_od", 35, 35, 0.0), ("box 72 px, origin fraction .25", "gray_inverted", 72, 72, -0.25)]
    accepted, refused, used = [], [], set()
    for label, channel, w, h, shift in specs:
        # Prefer an annotation not used yet whose template size changes with the window; else one whose centre does.
        found = {}
        for _, a in mit.iterrows():
            cx, cy = int(a["cx"]), int(a["cy"])
            if (float(cx), float(cy)) != (float(a["cx"]), float(a["cy"])) or int(a["ann_id"]) in used:
                continue
            x0, y0 = cx + 0.5 - w / 2 + shift, cy + 0.5 - h / 2
            if hand_centre(x0, w) != cx or hand_centre(y0, h) != cy:
                continue
            side = max(hand_pixels(w), hand_pixels(h))
            got = ss.tightened_template_box(channels[channel], float(cx), float(cy), otsu_window=side)
            g51 = ss.tightened_template_box(channels[channel], float(cx), float(cy), otsu_window=51)
            if side != 51 and got is not None and got != g51 and ss._patch_readable(rgb.shape, got[1], got[2], tm.PATCH_SIZE):
                kind = "size" if g51 is None or got[0] != g51[0] else "centre"
                found.setdefault(kind, dict(label=label, channel=channel, box=(x0, y0, w, h), centre=(cx, cy), side=side, expect=got, at_51=g51, ann_id=int(a["ann_id"]), differs_in=kind))
                if kind == "size":
                    break
        pick = found.get("size") or found.get("centre")
        if pick:
            accepted.append(pick)
            used.add(pick["ann_id"])
    # Refusals: the square and the non-square box at points the gate refuses at their own window (background near
    # 14499, searched outward), verified before use.
    rng = np.random.default_rng(0)
    for label, w, h in [("refused square 41 px", 41, 41), ("refused non-square 25 x 57 px", 25, 57)]:
        side = max(hand_pixels(w), hand_pixels(h))
        for _ in range(500):
            cx, cy = int(rng.integers(200, rgb.shape[1] - 200)), int(rng.integers(200, rgb.shape[0] - 200))
            if ss.tightened_template_box(gray, float(cx), float(cy), otsu_window=side) is None:
                x0, y0 = cx + 0.5 - w / 2, cy + 0.5 - h / 2
                assert hand_centre(x0, w) == cx and hand_centre(y0, h) == cy
                refused.append(dict(label=label, channel="gray_inverted", box=(x0, y0, w, h), centre=(cx, cy), side=side))
                break
    return accepted, refused


def check_gate_and_fallback(tally, log, R):
    """
        3-4. The gate runs on the box's own window (accepted boxes), and a refused box falls back to the human's box.

        tally (Tally): where to record.
        log (dict): per-case records.
        R (Rois): the ROI loader.

        Returns None.
    """
    roi, rgb, mpp, gt = R.get(FN)
    accepted, refused = gate_cases(rgb, gt)
    tally.flag("C3 at least 4 verified accepted boxes off 51", len(accepted) >= 4, f"only {len(accepted)} accepted cases found")
    tally.flag("C4 verified refused boxes found (square and non-square)", len(refused) == 2, f"only {len(refused)} refused cases found")
    cases = []
    for c in accepted:
        seed = af.box_frame(FN, *c["box"], 1)
        frame, info = af.run_roi(rgb, seed, mpp, refinement_channel=c["channel"])
        b, tx, ty = c["expect"]
        ref = frame[(frame["made_by"] == af.ALGORITHM) & (frame["class"] == af.MITOTIC_GT)].iloc[0]
        ok = info["refined"] is True and info["refusal_reason"] is None and info["base_size"] == b and C.same_value(tuple(info["template_center"]), (float(tx), float(ty))) and info["base_size"] <= c["side"] and info["n_blanked_px"] == info["base_size"] ** 2 and geometry(ref) == af.footprint_box(b, tx, ty) and C.same_value(tuple(info["click_center"]), (float(c["centre"][0]), float(c["centre"][1])))
        tally.flag("C3 the gate uses the box's window", ok, runs=1, what=f"{c['label']}: info base {info['base_size']} at {info['template_center']} refined {info['refined']} {info['refusal_reason']}; expected {c['expect']} (51 px gives {c['at_51']}); blanked {info['n_blanked_px']}; refined row {geometry(ref)}")
        cases.append(dict(label=c["label"], channel=c["channel"], box=c["box"], ann_id=c["ann_id"], differs_in=c["differs_in"], window=c["side"], base_size=info["base_size"], template_center=info["template_center"], at_51=c["at_51"], n_detections=info["n_detections"]))

    for c in refused:
        seed = af.box_frame(FN, *c["box"], 1)
        frame, info = af.run_roi(rgb, seed, mpp, refinement_channel=c["channel"])
        f_none, i_none = af.run_roi(rgb, seed, mpp, refinement_channel=None)
        human = frame[frame["made_by"] == af.HUMAN].iloc[0]
        ref = frame[(frame["made_by"] == af.ALGORITHM) & (frame["class"] == af.MITOTIC_GT)].iloc[0]
        cx, cy = c["centre"]
        half = c["side"] // 2
        enclosing = (float(cx - half), float(cy - half), float(c["side"] - 1), float(c["side"] - 1))
        square = c["box"][2] == c["box"][3]
        ok = info["refined"] is False and info["refusal_reason"] == af.GATE_REFUSED and info["base_size"] == c["side"] and C.same_value(tuple(info["template_center"]), (float(cx), float(cy))) and geometry(ref) == enclosing and (geometry(ref) == geometry(human) if square else geometry(ref) != geometry(human))
        tally.flag("C4 a refused box falls back to the human's box", ok, runs=1, what=f"{c['label']}: base {info['base_size']} at {info['template_center']} {info['refusal_reason']}; refined row {geometry(ref)}, human {geometry(human)}, enclosing {enclosing}")
        tally.add("C4 fallback predictions identical to refinement_channel=None", C.compare_frames(frame[frame["class"] == af.PREDICTION], f_none[f_none["class"] == af.PREDICTION], ["x_tl", "y_tl", "width", "height", "class", "object_id", "made_by", "refined", "score", "od"], label=c["label"]), runs=1)
        tally.flag("C4 refinement_channel=None uses the human's box", i_none["base_size"] == c["side"] and i_none["refusal_reason"] is None, f"{c['label']}: None gave base {i_none['base_size']}")
        cases.append(dict(label=c["label"], box=c["box"], window=c["side"], base_size=info["base_size"], refusal_reason=info["refusal_reason"], refined_row=geometry(ref), human_row=geometry(human)))

    # template_unreadable with a box window: search boxes centred on V5's unreadable click; record the outcome either way.
    cx, cy = 6410, 250
    unreadable = None
    for channel in ("gray_inverted", "hematoxylin_od"):
        for w in (41, 45, 55, 61, 71):
            got = ss.tightened_template_box(ch.to_channel(rgb, channel), float(cx), float(cy), otsu_window=w)
            if got is not None and not ss._patch_readable(rgb.shape, got[1], got[2], tm.PATCH_SIZE):
                unreadable = (channel, w)
                break
        if unreadable:
            break
    rec = dict(click=(cx, cy), found=unreadable)
    if unreadable:
        channel, w = unreadable
        seed = af.box_frame(FN, cx + 0.5 - w / 2, cy + 0.5 - w / 2, w, w, 1)
        frame, info = af.run_roi(rgb, seed, mpp, refinement_channel=channel)
        ok = info["refusal_reason"] == af.TEMPLATE_UNREADABLE and info["base_size"] == w and C.same_value(tuple(info["template_center"]), (float(cx), float(cy)))
        tally.flag("C4 template_unreadable box falls back to the human's box", ok, runs=1, what=f"({cx}, {cy}) {channel} {w}: {info['refusal_reason']} base {info['base_size']}")
        rec.update(channel=channel, window=w, refusal_reason=info["refusal_reason"], base_size=info["base_size"])
    log["C3"] = cases[:len(accepted)]
    log["C4"] = dict(refused=cases[len(accepted):], template_unreadable=rec)


def check_edges(tally, log, R):
    """
        5. A 73 px box centred exactly 36 px from an edge runs, and the window guard never fires on valid input; the guard
        itself raises on an unreadable window.

        tally (Tally): where to record.
        log (dict): per-case records.
        R (Rois): the ROI loader.

        Returns None.
    """
    roi, rgb, mpp, gt = R.get(FN)
    h, w = rgb.shape[:2]
    cases = []
    for label, box in [("left edge", (0.0, float(h // 2 - 36), 73.0, 73.0)), ("bottom edge", (float(w // 2 - 36), float(h - 73), 73.0, 73.0))]:
        seed = af.box_frame(FN, *box, 1)
        r, err = raises(lambda: af.run_roi(rgb, seed, mpp))
        centre = (float(seed["x_tl"].iat[0]) + 36, float(seed["y_tl"].iat[0]) + 36)
        tally.flag("C5 a 73 px box 36 px from an edge runs", not r and (centre[0] == 36.0 or centre[1] == float(h - 37)), runs=1, what=f"{label} {box} centre {centre}: {err}")
        cases.append(dict(label=label, box=box, centre=centre, error=err))
    r, err = raises(lambda: af._template_choice(rgb, 30.0, 2000.0, "gray_inverted", 73))
    tally.flag("C5 the window guard raises on an unreadable window", r and err.startswith("RuntimeError"), f"guard did not raise: {err}")
    cases.append(dict(label="guard at x = 30, window 73", error=err))
    log["C5"] = cases


def check_input_rows(tally, log, R):
    """
        6. run_roi raises on a seed or non-seed input row whose box gives no valid window, whichever row is the seed.

        tally (Tally): where to record.
        log (dict): per-case records.
        R (Rois): the ROI loader.

        Returns None.
    """
    roi, rgb, mpp, gt = R.get(FN)
    good = af.human_seed_frame(gt[gt["ann_id"] == 14499].iloc[0], FN)
    cases = []
    bad_seed = good.assign(width=80.0, height=80.0)
    r, err = raises(lambda: af.run_roi(rgb, bad_seed, mpp))
    tally.flag("C6 input rows outside the window limits raise", r and err.startswith("ValueError") and "gate window" in err, f"width-80 seed: {err}")
    cases.append(dict(case="seed width 80", error=err))
    other = af.human_seed_frame(gt[gt["ann_id"] == 14578].iloc[0], FN)
    tally.flag("C6 positive control: a valid two-row input passes the input checks", not raises(lambda: af._input_rows(pd.concat([good, other], ignore_index=True)))[0], "valid two-row input raised")
    for what, bad in [("width 80", other.assign(width=80.0)), ("width 16 (window 17)", other.assign(width=16.0, height=16.0)), ("NaN width", other.assign(width=np.nan))]:
        two = pd.concat([good, bad], ignore_index=True)
        for seed_id in (14499, 14578):
            r, err = raises(lambda: af.run_roi(rgb, two, mpp, seed_object_id=seed_id))
            tally.flag("C6 input rows outside the window limits raise", r and err.startswith("ValueError") and "object_id 14578" in err, f"non-seed {what}, seed {seed_id}: {err}")
            cases.append(dict(case=f"row 14578 {what}, seed {seed_id}", error=err))
    log["C6"] = cases


def main():
    """
        Run every box check and write out/box/box.json.

        Returns None.
    """
    t0 = time.time()
    tally, log = Tally(), {}
    R = Rois()
    for key, step in [("C1", lambda: check_window_vectors(tally, log)), ("C2", lambda: check_box_frame(tally, log)), ("C6", lambda: check_input_rows(tally, log, R)), ("C5", lambda: check_edges(tally, log, R)), ("C3-C4", lambda: check_gate_and_fallback(tally, log, R))]:
        t = time.time()
        step()
        print(f"{key} done [{time.time() - t:.0f}s]", flush=True)
    out = dict(checks=tally.d, cases=log, t_total_s=round(time.time() - t0, 1))
    path = C.OUT / "box" / "box.json"
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
