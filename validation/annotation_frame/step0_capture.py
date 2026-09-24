"""
    Step 0 capture, from the unedited code: B1 (rank_key="chromatin_od", all 441 runs), B2 (rank_key="tm_score", the 126 runs of
    the 14 original_14 ROIs), each ROI's rr.select_clicks draw, and V9's augmented-bank outputs. Every run's seed is built exactly
    as rr.run_roi builds it (rerun_bbox3way_postD11.py:232-243). One pickle per ROI under out/step0/; an ROI whose pickle exists
    is skipped, so an interrupted pass resumes. Usage: step0_capture.py [--shard I/N] [--rois A,B] [--v9-tag pre]
"""

from __future__ import annotations

import argparse
import gc
import time

import numpy as np

import common as C
from common import ds, prod, rr, ss, ch


def capture_v9(rgb, mpp, ann):
    """
        V9's augmented-bank find_and_suppress on a 1024 px crop of the full-ROI hematoxylin channel around one click.

        rgb (np.ndarray): the ROI's RGB image.
        mpp (float): microns per pixel.
        ann (pd.Series): the reference seed-0 click.

        Returns dict: crop origin, click in crop, detections and info.
    """
    from midog_utils.find_and_suppress import find_and_suppress
    cx, cy = float(ann["cx"]), float(ann["cy"])
    x0, y0 = C.v9_crop_box(rgb.shape, cx, cy)
    hem = ch.to_channel(rgb, "hematoxylin_od")
    crop = np.ascontiguousarray(hem[y0:y0 + C.V9_CROP, x0:x0 + C.V9_CROP])
    del hem
    click_in_crop = (cx - x0, cy - y0)
    detections, info = find_and_suppress(crop, click_in_crop, C.v9_config(mpp))
    return dict(x0=x0, y0=y0, click_in_crop=click_in_crop, ann_id=int(ann["ann_id"]), detections=detections, info=info)


def capture_roi(roi, annotations, v9_names):
    """
        Capture one ROI's draw, B1 runs, B2 runs (original_14 only) and, for V9's ROIs, the augmented-bank output.

        roi (pd.Series): one build_rois row.
        annotations (pd.DataFrame): the full annotation table.
        v9_names (list[str]): the five V9 ROIs.

        Returns dict: everything captured for the ROI.
    """
    t0 = time.time()
    rgb, mpp, gt, pool, flagged, n_agreement, n_joint_valid, draws = rr.select_clicks(roi, annotations)
    runs = {}
    rank_keys = ["chromatin_od"] + (["tm_score"] if roi["subset"] == "original_14" else [])
    for seed_index, ann, specs, refusals, n_seed_pool, digests in draws:
        click_xy = (float(ann["cx"]), float(ann["cy"]))
        seed_ann_id = int(ann["ann_id"])
        for condition, (base_size, tx, ty) in specs.items():
            offset = float(np.hypot(tx - click_xy[0], ty - click_xy[1]))
            seed = ss.Seed(ann_id=seed_ann_id, click_xy=click_xy, template_xy=(tx, ty), base_size=base_size, recentred=rr.CONDITIONS[condition] is not None, offset_px=offset, n_retries=len(refusals), agreement_flagged=bool(flagged), n_agreement_pool=n_agreement, n_after_border=len(pool))
            for rank_key in rank_keys:
                detections, info = prod.run_production_pipeline(rgb, seed, mpp, rank_key=rank_key)
                runs[(seed_index, condition, rank_key)] = (detections, info)
    v9 = capture_v9(rgb, mpp, draws[0][1]) if roi["file_name"] in v9_names else None
    out = dict(file_name=roi["file_name"], subset=roi["subset"], shape=rgb.shape, mpp=mpp, gt=gt, pool=pool, flagged=flagged, n_agreement=n_agreement, n_joint_valid=n_joint_valid, draws=draws, runs=runs, v9_pre=v9, t_capture_s=round(time.time() - t0, 1))
    del rgb
    gc.collect()
    return out


def main():
    """
        Capture every requested ROI that has no pickle yet.

        Returns None.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", default="0/1", help="I/N: capture ROIs whose build_rois position is I mod N")
    parser.add_argument("--rois", default="", help="comma-separated file names to restrict to")
    args = parser.parse_args()
    shard_i, shard_n = (int(v) for v in args.shard.split("/"))

    rois, annotations = rr.build_rois()
    v9_names = C.v9_rois(rois)
    print(f"V9 ROIs: {v9_names}", flush=True)
    wanted = [s.strip() for s in args.rois.split(",") if s.strip()]
    t_run = time.time()
    for i, roi in rois.iterrows():
        if i % shard_n != shard_i or (wanted and roi["file_name"] not in wanted):
            continue
        path = C.STEP0 / f"{roi['file_name']}.pkl"
        if path.exists():
            print(f"[{i + 1}/49 {roi['file_name']}] exists, skipped", flush=True)
            continue
        got = capture_roi(roi, annotations, v9_names)
        C.save_pickle(got, path)
        print(f"[{i + 1}/49 {roi['file_name']}] {roi['subset']}: {len(got['runs'])} runs, clicks {[int(d[1]['ann_id']) for d in got['draws']]}, v9={'yes' if got['v9_pre'] else 'no'} [{got['t_capture_s']:.0f}s; {(time.time() - t_run) / 60:.1f} min]", flush=True)


if __name__ == "__main__":
    main()
