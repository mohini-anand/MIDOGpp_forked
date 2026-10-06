"""
    Checks for multi-iteration sessions (annotation_frame.run_iteration; ANNOTATION_DATAFRAME_SPEC.md §3 "Sessions" and
    DECISIONS.md D14): iteration 0 is exactly run_roi's run; blank_squares paints blank_seed_square's footprint; the template
    is cut before the blanks and the search sees them; later iterations blank every stored box, store k non-repeat
    predictions and nothing else; the session survives a CSV round trip; the overlap guard is inclusive and box-based; bad
    sessions and seeds raise. Expected values are worked out here, never read back from the code under test. Writes
    out/session/session.json. Usage: session_checks.py
"""

from __future__ import annotations

import json
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace

import common as C
from common import np, pd, prod, tm
from checks import Tally
from midog_utils import annotation_frame as af
from midog_utils import find_and_suppress as fs
from named_checks import Rois

K = 10
SESSION_ROIS = ("300.tiff", "245.tiff")
N_ITERATIONS = 3
MERGE_KEYS = ("iteration", "k", "n_shown", "shown_object_ids", "repeat_object_ids", "n_repeats_in_page", "session_blanks")
EXTRA_BLANK_KEYS = ("n_extra_blanks", "n_extra_blanked_px")
VALUE_COLUMNS = [c for c in af.COLUMNS if c != "current_time"]
SESSION_VALUE_COLUMNS = VALUE_COLUMNS + [af.ITERATION]


def centres(rows):
    """
        Box centres of frame rows, x_tl + width / 2.

        rows (pd.DataFrame): frame or session rows.

        Returns tuple[np.ndarray, np.ndarray]: cx, cy.
    """
    return (rows["x_tl"] + rows["width"] / 2).to_numpy(dtype=np.float64), (rows["y_tl"] + rows["height"] / 2).to_numpy(dtype=np.float64)


def expected_page(frame, k):
    """
        Iteration 0's stored rows, by hand: the two seed rows, the first k predictions and their refined rows.

        frame (pd.DataFrame): run_roi's output.
        k (int): the page size.

        Returns pd.DataFrame: the rows, iteration 0.
    """
    seeds = frame[frame["class"] == af.MITOTIC_GT]
    preds = frame[(frame["class"] == af.PREDICTION) & (frame["object_id"] < k)]
    refined = frame[(frame["class"] == af.PREDICTION_REFINED) & (frame["object_id"] < k)]
    return pd.concat([seeds, preds, refined], ignore_index=True).assign(iteration=0)


def check_iteration_zero(tally, R):
    """
        S1: with no session, run_iteration's run is run_roi's exactly, it adds no production info key, and it stores the first
        k predictions; None and an empty session give the same session.

        tally (Tally): where to record.
        R (Rois): the ROI loader.

        Returns None.
    """
    cases = [("300.tiff", "midog0", "gray_inverted"), ("300.tiff", "midog0", None), ("300.tiff", "box41", "gray_inverted"), ("245.tiff", "midog0", "hematoxylin_od")]
    for fn, sname, channel in cases:
        roi, rgb, mpp, gt = R.get(fn)
        seed = af.draw_seed_rows(rgb, gt, int(roi["image_id"]), seed_indices=(0,))[0][0] if sname == "midog0" else af.box_frame(fn, 1391, 3531, 41, 41, 7)
        label = f"{fn} {sname} {channel}"
        f_ref, i_ref = af.run_roi(rgb, seed, mpp, refinement_channel=channel, refine_predictions=True)
        session, frame, info = af.run_iteration(rgb, None, mpp, seed, k=K, refinement_channel=channel)
        tally.add("S1 iteration 0's frame equals run_roi's", C.compare_frames(frame, f_ref, VALUE_COLUMNS, label=label), runs=1)
        tally.add("S1 iteration 0's info equals run_roi's apart from the merge keys", C.compare_info(info, i_ref, exclude=MERGE_KEYS + ("prediction_refinement",), label=label), runs=1)
        tally.flag("S1 iteration 0's prediction table equals run_roi's", info["prediction_refinement"].equals(i_ref["prediction_refinement"]), f"{label}: prediction_refinement differs", runs=1)
        tally.flag("S1 iteration 0 adds no production key", not any(k in info for k in EXTRA_BLANK_KEYS), f"{label}: {[k for k in EXTRA_BLANK_KEYS if k in info]}", runs=1)
        tally.add("S1 iteration 0 stores the seed rows and the first k predictions with their refined rows", C.compare_frames(session, expected_page(f_ref, K), SESSION_VALUE_COLUMNS, label=label), runs=1)
        merge_ok = info["iteration"] == 0 and info["n_shown"] == K and info["shown_object_ids"] == list(range(K)) and info["repeat_object_ids"] == [] and info["n_repeats_in_page"] == 0 and info["session_blanks"].empty
        tally.flag("S1 iteration 0's merge keys", merge_ok, f"{label}: {[(k, info[k]) for k in MERGE_KEYS if k != 'session_blanks']}", runs=1)
        if fn == "300.tiff" and sname == "midog0" and channel == "gray_inverted":
            s_empty, f_empty, _ = af.run_iteration(rgb, af.empty_session(), mpp, seed, k=K, refinement_channel=channel)
            tally.add("S1 an empty session is iteration 0, as None is", C.compare_frames(s_empty, session, SESSION_VALUE_COLUMNS, label=label), runs=1)
        print(f"S1 {label} done", flush=True)


def check_blank_squares(tally):
    """
        S2: blank_squares paints exactly blank_seed_square's footprint (same rounding and clipping), counts the union once, and
        rejects bad squares.

        tally (Tally): where to record.

        Returns None.
    """
    rng = np.random.default_rng(0)
    channel = rng.random((120, 160)).astype(np.float32) + 1.0
    fill = float(channel.min())
    for cx, cy, side in [(40.0, 50.0, 51), (40.5, 50.5, 51), (41.5, 49.5, 23), (2.0, 3.0, 21), (158.5, 119.0, 33), (80.0, 60.0, 1), (79.5, 60.5, 73)]:
        ref, n_ref = tm.blank_seed_square(channel, cx, cy, side)
        mine = channel.copy()
        n_sq, n_px = fs.blank_squares(mine, [(cx, cy, side)], fill)
        tally.flag("S2 one square equals blank_seed_square's footprint", np.array_equal(mine, ref) and n_sq == 1 and n_px == n_ref, f"({cx}, {cy}, {side}): {n_px} px vs {n_ref}")
    squares = [(30.0, 30.0, 21), (40.0, 30.0, 21), (100.0, 80.0, 11)]
    mine = channel.copy()
    n_sq, n_px = fs.blank_squares(mine, squares, fill)
    hand = np.zeros(channel.shape, bool)
    for cx, cy, side in squares:
        h = side // 2
        hand[int(cy) - h:int(cy) + h + 1, int(cx) - h:int(cx) + h + 1] = True
    tally.flag("S2 overlapping squares: union painted and counted once", n_sq == 3 and n_px == int(hand.sum()) == 21 * 31 + 121 and np.array_equal(mine == fill, hand | (channel == fill)) and np.array_equal(mine[~hand], channel[~hand]), f"{n_px} px, hand {int(hand.sum())}")
    mine = channel.copy()
    tally.flag("S2 no squares changes nothing", fs.blank_squares(mine, [], fill) == (0, 0) and np.array_equal(mine, channel), "empty list changed the channel")
    for bad in [(10.0, 10.0, 20), (10.0, 10.0, 0), (10.0, 10.0, 7.5), (float("nan"), 10.0, 21), (10.0, float("inf"), 21)]:
        try:
            fs.blank_squares(channel.copy(), [bad], fill)
            tally.flag("S2 a bad square raises", False, f"{bad} did not raise")
        except ValueError:
            tally.flag("S2 a bad square raises", True, "")


def check_cut_before_blank(tally, diag, R):
    """
        S3: production's template is cut from the unblanked channel even when a blank covers it, the correlation sees every
        extra square at the fill value, od reads the unblanked channel, and an empty list changes nothing but adds the keys.

        tally (Tally): where to record.
        diag (dict): logged diagnostics.
        R (Rois): the ROI loader.

        Returns None.
    """
    fn = "300.tiff"
    roi, rgb, mpp, gt = R.get(fn)
    seed = SimpleNamespace(base_size=49, template_xy=(1413.0, 3550.0))
    squares = [(1413.0, 3550.0, 73), (1450.0, 3520.0, 51), (3000.5, 2000.5, 31)]  # the first covers the whole template and its 73 px patch
    captured = {}
    build, fused = fs.tm.build_augmentations, fs.tm.fused_response

    def spy_build(patch, *a, **kw):
        captured.setdefault("patch", []).append(patch.copy())
        return build(patch, *a, **kw)

    def spy_fused(image, *a, **kw):
        captured.setdefault("image", []).append(image)
        return fused(image, *a, **kw)

    fs.tm.build_augmentations, fs.tm.fused_response = spy_build, spy_fused
    try:
        d0, i0 = prod.run_production_pipeline(rgb, seed, mpp, always_od=True)
        d1, i1 = prod.run_production_pipeline(rgb, seed, mpp, always_od=True, extra_blanks=squares)
        de, ie = prod.run_production_pipeline(rgb, seed, mpp, always_od=True, extra_blanks=[])
    finally:
        fs.tm.build_augmentations, fs.tm.fused_response = build, fused
    tally.flag("S3 the template is cut before the blanks", np.array_equal(captured["patch"][0], captured["patch"][1]), "the template patch changed under a covering blank", runs=2)
    hem = C.ch.to_channel(rgb, prod.CHANNEL)
    fill, pad, img = float(hem.min()), i1["pad_px"], captured["image"][1]
    ok, hand = True, np.zeros(hem.shape, bool)
    for cx, cy, side in squares:
        h, ix, iy = side // 2, int(round(cx)), int(round(cy))
        hand[iy - h:iy + h + 1, ix - h:ix + h + 1] = True
        ok &= bool((img[pad + iy - h:pad + iy + h + 1, pad + ix - h:pad + ix + h + 1] == fill).all())
    tally.flag("S3 the correlation sees every extra square at the fill value", ok, "an extra square is not at the fill value in the searched image", runs=1)
    tally.flag("S3 the extra blank counts", i1["n_extra_blanks"] == 3 and i1["n_extra_blanked_px"] == int(hand.sum()), f"{i1.get('n_extra_blanks')} squares, {i1.get('n_extra_blanked_px')} px vs {int(hand.sum())}", runs=1)
    tally.flag("S3 the seed's own blank is unchanged", i1["n_blanked_px"] == i0["n_blanked_px"] == 49 ** 2, f"n_blanked_px {i1['n_blanked_px']} vs {i0['n_blanked_px']}", runs=1)
    hits = np.zeros(len(d1), bool)
    for cx, cy, side in squares:
        hits |= (np.abs(d1["cx"].to_numpy() - cx) <= side // 2) & (np.abs(d1["cy"].to_numpy() - cy) <= side // 2)
    diag["s3_detections_inside_a_blank"] = int(hits.sum())  # logged, not judged: nothing guarantees a blank holds no peak
    merged = d1.merge(d0, on=["cx", "cy"], suffixes=("", "_ref"))
    tally.flag("S3 od reads the unblanked channel", len(merged) > 50 and np.array_equal(merged["od"].to_numpy(), merged["od_ref"].to_numpy()), f"{len(merged)} shared detections, od differs on {int((merged['od'] != merged['od_ref']).sum())}", runs=1, rows=len(merged))
    tally.add("S3 an empty list leaves the detections unchanged", C.compare_frames(de, d0, list(d0.columns), label="empty list"), runs=1)
    tally.add("S3 an empty list adds only the zero counts to info", C.compare_info(ie, i0, exclude=EXTRA_BLANK_KEYS, label="empty list"), runs=1)
    tally.flag("S3 an empty list's counts are zero", ie.get("n_extra_blanks") == 0 and ie.get("n_extra_blanked_px") == 0, f"{[ie.get(k) for k in EXTRA_BLANK_KEYS]}", runs=1)
    print("S3 done", flush=True)


def check_session_run(tally, log, R, fn, tmp):
    """
        S4: N_ITERATIONS iterations on one ROI with the validation draw's seeds 0..N-1: each blanks exactly the stored boxes,
        stores k non-repeat predictions and nothing else, and matches its own run_roi call; then a CSV round trip, and the
        last iteration re-run from the CSV session. Logs, without judging, how many of each page's predictions the same seed
        would have repeated without blanking.

        tally (Tally): where to record.
        log (list): per-iteration records.
        R (Rois): the ROI loader.
        fn (str): the ROI.
        tmp (Path): a scratch directory for the CSV.

        Returns None.
    """
    roi, rgb, mpp, gt = R.get(fn)
    seeds = [s for s, _ in af.draw_seed_rows(rgb, gt, int(roi["image_id"]), seed_indices=range(N_ITERATIONS))]
    radius = C.ev.radius_px(mpp)
    session = None
    for t, seed in enumerate(seeds):
        label = f"{fn} iteration {t}"
        before = af.empty_session() if session is None else session
        session, frame, info = af.run_iteration(rgb, session, mpp, seed, k=K)
        new = session[session[af.ITERATION] == t].reset_index(drop=True)

        stored = before[before["made_by"] == af.ALGORITHM]
        hand = stored[(stored["class"] == af.MITOTIC_GT) | (stored["class"] == af.PREDICTION_REFINED)]
        bx, by = centres(hand)
        sides = (hand["width"] + 1).to_numpy()
        blanks = info["session_blanks"]
        same_blanks = len(blanks) == len(hand) == t * (1 + K) and np.array_equal(blanks["cx"].to_numpy(), bx) and np.array_equal(blanks["cy"].to_numpy(), by) and np.array_equal(blanks["side"].to_numpy(), sides)
        tally.flag("S4 the blanks are every earlier refined seed and stored refined prediction", same_blanks, f"{label}: {len(blanks)} blanks, {len(hand)} by hand", runs=1)
        if t:
            tally.flag("S4 production blanked them all", info["n_extra_blanks"] == len(hand), f"{label}: n_extra_blanks {info['n_extra_blanks']}", runs=1)

        f_ref, i_ref = af.run_roi(rgb, seed, mpp, refine_predictions=True, extra_blanks=list(zip(bx, by, sides)) if t else None)
        tally.add("S4 the frame equals run_roi with the stored boxes blanked", C.compare_frames(frame, f_ref, VALUE_COLUMNS, label=label), runs=1)

        preds = f_ref[f_ref["class"] == af.PREDICTION]
        px, py = centres(preds)
        repeat = (np.hypot(px[:, None] - bx[None, :], py[:, None] - by[None, :]) <= radius).any(axis=1) if t else np.zeros(len(preds), bool)
        page = preds[~repeat].head(K)
        expected = pd.concat([f_ref[f_ref["class"] == af.MITOTIC_GT], page, f_ref[(f_ref["class"] == af.PREDICTION_REFINED) & f_ref["object_id"].isin(page["object_id"])]], ignore_index=True).assign(iteration=t)
        tally.add("S4 the iteration stores the seed rows and the first k non-repeats with their refined rows", C.compare_frames(new, expected, SESSION_VALUE_COLUMNS, label=label), runs=1)
        tally.flag("S4 the repeat keys", info["repeat_object_ids"] == preds.loc[repeat, "object_id"].tolist() and info["n_repeats_in_page"] == int(repeat[:K].sum()) and info["n_shown"] == len(page) == K, f"{label}: repeats {info['repeat_object_ids']} vs {preds.loc[repeat, 'object_id'].tolist()}", runs=1)
        if t:
            sx, sy = centres(new[new["class"] == af.PREDICTION])
            tally.flag("S4 no stored prediction is within the match radius of a blank", not (np.hypot(sx[:, None] - bx[None, :], sy[:, None] - by[None, :]) <= radius).any(), f"{label}: a stored prediction sits on a blank", runs=1)
        tally.add("S4 earlier iterations are untouched", C.compare_frames(session[session[af.ITERATION] < t], before, SESSION_VALUE_COLUMNS + ["current_time"], label=label), runs=1)

        record = dict(roi=fn, iteration=t, seed_ann_id=int(seed["object_id"].iat[0]), n_blanks=len(blanks), n_extra_blanked_px=info.get("n_extra_blanked_px", 0), n_repeats_in_page=info["n_repeats_in_page"], repeat_object_ids=info["repeat_object_ids"], t_total_s=info["t_total_s"])
        if t:
            f_plain, _ = af.run_roi(rgb, seed, mpp)
            qx, qy = centres(f_plain[f_plain["class"] == af.PREDICTION].head(K))
            record["unblanked_first_k_within_radius_of_a_stored_box"] = int((np.hypot(qx[:, None] - bx[None, :], qy[:, None] - by[None, :]) <= radius).any(axis=1).sum())
        log.append(record)
        print(f"S4 {label} done: {record}", flush=True)

    path = tmp / f"{fn}_session.csv"
    session.to_csv(path, index=False)
    back = af.check_session(pd.read_csv(path, float_precision="round_trip"))
    tally.add("S4 the session survives a CSV round trip", C.compare_frames(back, session, list(af.SESSION_COLUMNS), label=fn), runs=1)
    tally.flag("S4 the round trip keeps the dtypes", (back.dtypes == session.dtypes).all(), f"{fn}: {back.dtypes.to_dict()}", runs=1)
    prior = session[session[af.ITERATION] < N_ITERATIONS - 1]
    prior_csv = back[back[af.ITERATION] < N_ITERATIONS - 1]
    s_mem, f_mem, i_mem = af.run_iteration(rgb, prior, mpp, seeds[-1], k=K)
    s_csv, f_csv, i_csv = af.run_iteration(rgb, prior_csv, mpp, seeds[-1], k=K)
    tally.add("S4 the last iteration from the CSV session equals the in-memory one", C.compare_frames(s_csv, s_mem, SESSION_VALUE_COLUMNS, label=fn), runs=2)
    tally.add("S4 the last iteration's info from the CSV session equals the in-memory one", C.compare_info(i_csv, i_mem, exclude=("prediction_refinement", "session_blanks"), label=fn), runs=1)
    tally.flag("S4 and its blanks", i_csv["session_blanks"].equals(i_mem["session_blanks"]), f"{fn}: session_blanks differ", runs=1)


def check_guard_and_inputs(tally, R):
    """
        S5: the overlap guard raises exactly when the new box contains a stored centre (bounds inclusive) and can be switched
        off; used ids, another ROI's seed, bad sessions and a bad k raise; next_object_id. S6: without refined predictions
        the stored predictions' own 51 px boxes are blanked.

        tally (Tally): where to record.
        R (Rois): the ROI loader.

        Returns None.
    """
    fn = "245.tiff"

    def raises(label, fn_call):
        try:
            fn_call()
            tally.flag(label, False, "did not raise")
        except ValueError:
            tally.flag(label, True, "")

    def passes(label, fn_call):
        try:
            fn_call()
            tally.flag(label, True, "")
        except ValueError as exc:
            tally.flag(label, False, f"raised: {exc}")

    # A hand-built session whose only stored centres are a seed at (3000, 3000) and one fallback prediction at (1000, 1000),
    # so each case's expectation is the inclusive box test on that one centre, worked out here.
    now = af.now_utc()
    human = af.click_frame(fn, 3000.0, 3000.0, 5, current_time=now)
    pred = af._block(1, file_name=fn, x_tl=975.0, y_tl=975.0, width=50.0, height=50.0, **{"class": af.PREDICTION}, object_id=0, made_by=af.ALGORITHM, refined=False, current_time=now, score=1.0, od=0.1)
    synthetic = af.check_session(pd.concat([human, human.assign(made_by=af.ALGORITHM), pred, pred.assign(**{"class": af.PREDICTION_REFINED, "score": np.nan, "od": np.nan})], ignore_index=True).assign(iteration=0))
    for dx, dy in [(0.0, 0.0), (24.0, -24.0), (25.0, 0.0), (-25.0, 25.0), (25.5, 0.0), (0.0, -25.5), (40.0, 40.0)]:
        click = af.click_frame(fn, 1000.0 + dx, 1000.0 + dy, 6)
        inside = abs(dx) <= 25.0 and abs(dy) <= 25.0
        name = "S5 a box containing a stored centre raises (bounds inclusive)" if inside else "S5 a box beside a stored centre passes"
        (raises if inside else passes)(name, lambda click=click: af.prepare_iteration(synthetic, click))
    box = af.box_frame(fn, 960.0, 990.0, 73, 21, 6)  # pixel-edge box: centre pixel (996, 1000), 73 x 21 px, so x 960..1032 and y 990..1010
    raises("S5 a drawn box containing a stored centre raises", lambda: af.prepare_iteration(synthetic, box))
    passes("S5 check_overlap=False lets it through", lambda: af.prepare_iteration(synthetic, box, check_overlap=False))
    passes("S5 a drawn box beside it passes", lambda: af.prepare_iteration(synthetic, af.box_frame(fn, 960.0, 1026.0, 73, 21, 6)))

    roi, rgb, mpp, gt = R.get(fn)
    seed0 = af.draw_seed_rows(rgb, gt, int(roi["image_id"]), seed_indices=(0,))[0][0]
    session, _, _ = af.run_iteration(rgb, None, mpp, seed0, k=K)
    sx, sy = centres(session[session["class"] == af.PREDICTION])
    margin, (h, w) = tm.PATCH_SIZE // 2, rgb.shape[:2]
    first = int(np.flatnonzero((sx >= margin) & (sy >= margin) & (sx < w - margin) & (sy < h - margin))[0])  # a seed needs a 36 px margin
    pcx, pcy = float(sx[first]), float(sy[first])
    next_id = af.next_object_id(session)
    raises("S5 a click on a stored prediction raises", lambda: af.prepare_iteration(session, af.click_frame(fn, pcx, pcy, next_id)))
    s_over, f_over, i_over = af.run_iteration(rgb, session, mpp, af.click_frame(fn, pcx, pcy, next_id), k=K, check_overlap=False)
    tally.flag("S5 a seed on a blanked prediction still runs, iteration 1", i_over["iteration"] == 1 and i_over["n_extra_blanks"] == 1 + K, f"iteration {i_over['iteration']}", runs=1)

    raises("S5 a used object_id raises", lambda: af.prepare_iteration(session, af.click_frame(fn, 3000.5, 3000.5, int(seed0["object_id"].iat[0]))))
    raises("S5 another ROI's seed raises", lambda: af.prepare_iteration(session, af.click_frame("300.tiff", 3000.5, 3000.5, next_id)))
    raises("S5 a two-row seed raises", lambda: af.prepare_iteration(session, pd.concat([af.click_frame(fn, 3000.5, 3000.5, next_id), af.click_frame(fn, 3500.5, 3000.5, next_id + 1)], ignore_index=True)))
    raises("S5 k = 0 raises", lambda: af.run_iteration(rgb, session, mpp, af.click_frame(fn, 3000.5, 3000.5, next_id), k=0))
    bad = {
        "no iteration column": session.drop(columns=af.ITERATION),
        "an extra column": session.assign(verdict="ACCEPTED"),
        "iterations 0 and 2": pd.concat([session, s_over[s_over[af.ITERATION] == 1].assign(iteration=2)], ignore_index=True),
        "two file_names": pd.concat([session, s_over[s_over[af.ITERATION] == 1].assign(file_name="300.tiff")], ignore_index=True),
        "a second human seed in iteration 0": pd.concat([session, s_over[(s_over[af.ITERATION] == 1) & (s_over["made_by"] == af.HUMAN)].assign(iteration=0)], ignore_index=True),
        "a refined prediction without its prediction": session.drop(index=session.index[session["class"] == af.PREDICTION][0]),
        "a non-square stored box": session.assign(width=np.where(session["class"] == af.PREDICTION, 48.0, session["width"])),
        "an even stored side": session.assign(width=np.where(session["class"] == af.PREDICTION, 51.0, session["width"]), height=np.where(session["class"] == af.PREDICTION, 51.0, session["height"])),
        "a duplicated row": pd.concat([session, session.iloc[[5]]], ignore_index=True),
    }
    for what, frame in bad.items():
        raises(f"S5 a bad session raises: {what}", lambda frame=frame: af.check_session(frame))
    tally.flag("S5 next_object_id", af.next_object_id(None) == 0 and af.next_object_id(af.empty_session()) == 0 and next_id == int(seed0["object_id"].iat[0]) + 1 and af.next_object_id(s_over) == next_id + 1, f"next ids {af.next_object_id(None)}, {next_id}, {af.next_object_id(s_over)}")

    s_plain, _, _ = af.run_iteration(rgb, None, mpp, seed0, k=K, refine_predictions=False)
    blanks = af.stored_blanks(s_plain)
    preds = s_plain[s_plain["class"] == af.PREDICTION]
    px, py = centres(preds)
    ok = (s_plain["class"] != af.PREDICTION_REFINED).all() and len(blanks) == 1 + K and (blanks["side"].to_numpy()[1:] == 51).all() and np.array_equal(blanks["cx"].to_numpy()[1:], px) and np.array_equal(blanks["cy"].to_numpy()[1:], py)
    tally.flag("S6 without refined predictions their own 51 px boxes are blanked", bool(ok), f"{len(blanks)} blanks, sides {blanks['side'].tolist()}", runs=1)
    print("S5, S6 done", flush=True)


def main():
    """
        Run every session check and write out/session/session.json.

        Returns None.
    """
    t0 = time.time()
    tally, log, diag = Tally(), [], {}
    R = Rois()
    check_blank_squares(tally)
    check_iteration_zero(tally, R)
    check_cut_before_blank(tally, diag, R)
    with tempfile.TemporaryDirectory() as tmp:
        for fn in SESSION_ROIS:
            check_session_run(tally, log, R, fn, Path(tmp))
    check_guard_and_inputs(tally, R)
    out = dict(checks=tally.d, iterations=log, diagnostics=diag, t_total_s=round(time.time() - t0, 1))
    path = C.OUT / "session" / "session.json"
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
