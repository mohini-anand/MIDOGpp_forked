"""
    Step 2: the module's pieces on one ROI before any long run -- draw_seed_rows against rr.select_clicks (via the Step 0
    capture), the section 2 conversions, and run_roi for each refinement setting against B1. Usage: step2_pieces.py [--roi F]
"""

from __future__ import annotations

import argparse
import json

import numpy as np

import common as C
from common import ds, rr
import checks
from checks import Tally
from midog_utils import annotation_frame as af


def conversions(tally):
    """
        The section 1/2 conversions: footprint_box's half-to-even rounding and +1, the odd-size guard, and the human row's
        exact recovery of MIDOG x1/y1 on every annotation in the JSON.

        tally (Tally): where to record.

        Returns dict: what was checked.
    """
    cases = [((51, 100.0, 100.0), (75.0, 75.0, 50.0, 50.0)), ((51, 100.5, 100.5), (75.0, 75.0, 50.0, 50.0)), ((51, 101.0, 101.5), (76.0, 77.0, 50.0, 50.0)), ((51, 101.5, 102.5), (77.0, 77.0, 50.0, 50.0)), ((27, 6414.5, 257.0), (6401.0, 244.0, 26.0, 26.0)), ((33, 37.0, 281.0), (21.0, 265.0, 32.0, 32.0)), ((5, 2.5, 3.5), (0.0, 2.0, 4.0, 4.0))]
    for (b, tx, ty), want in cases:
        got = af.footprint_box(b, tx, ty)
        tally.flag("step2 footprint_box", got == want and all(type(v) is float for v in got), f"footprint_box{(b, tx, ty)} = {got} != {want}")
    for b in (50, 52):
        try:
            af.footprint_box(b, 100.0, 100.0)
            ok = False
        except ValueError:
            ok = True
        tally.flag("step2 footprint_box rejects even sizes", ok, f"footprint_box({b}) did not raise")

    raw = json.loads(open(C.WORKTREE / "databases/MIDOG++.json").read())
    x1 = {a["id"]: a["bbox"] for a in raw["annotations"]}
    _, ann = ds.load_annotations("databases/MIDOG++.json")
    frames = [af.human_seed_frame(r, r["file_name"], "t") for _, r in ann.iterrows()]
    bad = 0
    for f, (_, r) in zip(frames, ann.iterrows()):
        bb = x1[int(r["ann_id"])]
        row = f.iloc[0]
        ok = row["x_tl"] == float(bb[0]) and row["y_tl"] == float(bb[1]) and row["x_tl"] + row["width"] / 2 == float(r["cx"]) and row["y_tl"] + row["height"] / 2 == float(r["cy"]) and row["width"] == 50.0 and row["height"] == 50.0
        bad += not ok
    tally.d.setdefault("step2 human row == MIDOG x1/y1 and centre round-trips", dict(n_runs=0, n_rows=len(frames), n_values=4 * len(frames), n_diverged=bad, detail=[]))
    return dict(n_annotations=len(frames), n_bad=bad)


def main():
    """
        Run the step 2 checks on one ROI and print them.

        Returns None.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--roi", default="300.tiff")
    args = parser.parse_args()
    rois, annotations = rr.build_rois()
    roi = rois[rois["file_name"] == args.roi].iloc[0]
    s0 = C.load_pickle(C.STEP0 / f"{args.roi}.pkl")
    tally = Tally()
    rec = conversions(tally)

    rgb, mpp = ds.load_roi(roi["path"]), ds.roi_mpp(roi["path"])
    gt = ds.image_annotations(annotations, args.roi)
    draws3 = af.draw_seed_rows(rgb, gt, int(roi["image_id"]))
    draws6 = af.draw_seed_rows(rgb, gt, int(roi["image_id"]), seed_indices=range(6))
    checks.check_draws(tally, roi, s0, draws3, draws6, C.read_ref(C.R4_CLICKS_PATH))

    seed_frame, di = draws3[0]
    for condition, channel in C.CONDITION_CHANNEL.items():
        b1_det, b1_info = s0["runs"][(0, condition, "chromatin_od")]
        frame, info = af.run_roi(rgb, seed_frame, mpp, rank_key="chromatin_od", refinement_channel=channel)
        lbl = f"{args.roi} s0 {condition}"
        tally.add("step2 run_roi predictions vs B1", C.compare_frames(frame[frame["class"] == "PREDICTION"], checks.convert_detections(b1_det, args.roi), checks.COMPARED_FRAME_COLUMNS, label=lbl), runs=1)
        tally.add("step2 run_roi production info vs B1", C.compare_info({k: v for k, v in info.items() if k not in checks.ADDED_INFO_KEYS}, b1_info, label=lbl), runs=1)
        b, tx, ty = di["specs"][condition]
        tally.flag("step2 run_roi refined row", tuple(frame.iloc[1][["x_tl", "y_tl", "width", "height"]]) == af.footprint_box(b, tx, ty) and info["refined"] == (channel is not None), f"{lbl}: refined row {tuple(frame.iloc[1][['x_tl', 'y_tl', 'width', 'height']])}")
        print(f"{lbl}: {len(frame)} rows, refined={info['refined']}, template=({info['base_size']}, {info['template_center']}), info keys added={sorted(set(info) - set(b1_info))}")
        print(frame.head(4).to_string())
        print(frame.dtypes.to_dict())

    for name in sorted(tally.d):
        e = tally.d[name]
        print(f"  {name}: runs={e['n_runs']} rows={e['n_rows']} values={e['n_values']} diverged={e['n_diverged']}")
        for line in e["detail"][:5]:
            print(f"      {line}")
    (C.OUT / "step2").mkdir(parents=True, exist_ok=True)
    (C.OUT / "step2" / f"pieces_{args.roi}.json").write_text(json.dumps(dict(checks=tally.d, conversions=rec), indent=1, default=str))


if __name__ == "__main__":
    main()
