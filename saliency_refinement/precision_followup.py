"""
    Follow-up to precision.py, on the two questions it left open.

    1. The clicks where the saliency rule and production's gray gate actually pick differently: EVERY gray-refused
       annotation in out/geometry_midog_box.csv (production falls back to the untightened 51 px box there) plus the
       gray-accepted annotations whose saliency winner lies outside the gray component. Arms: default_51 (= production's
       fallback), prod_51, sal_51, sal1_51 as in precision.py, plus sal_fill_51.
    2. The notes' own template construction: a fixed 51 px window on the click, the saliency keep mask (dilated 1 um) kept
       as real pixels and everything else filled with the mean of the filled pixels (sal_fill_51) -- the repo's
       hem_mask_fill_1um arm with only the keep mask swapped. Run on the 294 joint clicks too, next to a re-derivation of
       hem_mask_fill_1um that must reproduce the recorded fillgeom postD11 rows (template SHA-1 and TP@10/20/30, both axes)
       before anything else runs. No feathering (the notes measured it as 0.000 +/- 0.005).

    Writes out/precision_followup_per_run.csv. Usage: precision_followup.py [--workers N] [--rois A,B]
"""

import argparse
import hashlib
from multiprocessing import Pool

import cv2
import numpy as np
import pandas as pd
from skimage.measure import label, regionprops

import common as C
import precision as P
from midog_utils.nms import nms_by_distance

W = C.tm.BASE_SIZE           # 51
OFF = C.HALF_PATCH - W // 2  # 11: the 51 px window's origin inside the 73 px patch


def fill_outside(T: np.ndarray, keep: np.ndarray):
    """
        The fillgeom notebook's fill: every pixel outside keep set to the float64 mean of those pixels, cast to float32.

        T (np.ndarray): 51x51 float32 window.
        keep (np.ndarray): 51x51 boolean keep mask.

        Returns tuple[np.ndarray, int]: the filled template and the number of filled pixels (0 = untouched window).
    """
    outside = ~keep
    if not outside.any():
        return np.ascontiguousarray(T, dtype=np.float32), 0
    F = T.copy()
    F[outside] = np.float32(T[outside].astype(np.float64).mean())
    return np.ascontiguousarray(F, dtype=np.float32), int(outside.sum())


def hem_mask_fill_1um(T: np.ndarray, mpp: float):
    """
        The recorded hem_mask_fill_1um template: the hem window's own Otsu component under the click, dilated by 1 um.

        T (np.ndarray): 51x51 hematoxylin OD window on the click.
        mpp (float): microns per pixel.

        Returns np.ndarray or None: the template, or None when the hem gate refuses.
    """
    if C.ss.tighten_box_otsu(T) is None:
        return None
    u8 = cv2.normalize(T.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    _, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    labels = label(binary, connectivity=2)
    comp = labels == labels[W // 2, W // 2]
    m = int(round(1.0 / mpp))
    se = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * m + 1, 2 * m + 1))
    return fill_outside(T, cv2.dilate(comp.astype(np.uint8), se).astype(bool))[0]


def sal_fill_51(hem: np.ndarray, ann: pd.Series, mpp: float):
    """
        The notes' template: saliency keep mask (dilated 1 um) inside the MIDOG box, fixed 51 px window, constant fill.

        hem (np.ndarray): the ROI's hematoxylin OD channel.
        ann (pd.Series): the click's annotation row.
        mpp (float): microns per pixel.

        Returns tuple[np.ndarray, int]: the template (the untouched window when the box holds no foreground) and filled px.
    """
    T = C.tm.read_padded_patch(hem, float(ann.cx), float(ann.cy), W)
    _, sal = C.saliency_spec(hem, *C.midog_box(ann), mpp, W, dilate=True)
    if sal["mask_1um"] is None:
        return np.ascontiguousarray(T, dtype=np.float32), 0
    return fill_outside(T, sal["mask_1um"][OFF:OFF + W, OFF:OFF + W])


def search_template(roi: dict, template: np.ndarray, cx: float, cy: float):
    """
        The fillgeom notebook's search_once for a 51 px template on the click, post-D11: blank the 51 px footprint, pad,
        correlate, threshold, peaks, NMS, od on the unblanked channel.

        roi (dict): load_roi output.
        template (np.ndarray): 51x51 float32 template.
        cx (float): click x.
        cy (float): click y.

        Returns tuple[pd.DataFrame, dict]: detections in score order with od, and the counts.
    """
    hem = roi["hem"]
    h, w = hem.shape
    search, n_blanked = C.tm.blank_seed_square(hem, cx, cy, W)
    pad = C.OD_PAD
    padded = cv2.copyMakeBorder(search, pad, pad, pad, pad, borderType=cv2.BORDER_REPLICATE)
    fused_p, _, valid_p = C.tm.fused_response(padded, [template], False, method=C.prod.TM_METHOD)
    fused, valid = fused_p[pad:pad + h, pad:pad + w], valid_p[pad:pad + h, pad:pad + w]
    assert valid.all()
    med, mad = C.tm.robust_stats(fused, valid)
    centers, scores = C.tm.extract_peaks(fused, valid, C.prod.PEAK_MIN_DISTANCE, med + C.prod.DEEP_FLOOR_Z * mad, C.prod.MAX_PEAKS)
    n_peaks = len(centers)
    keep = nms_by_distance(centers, scores, C.ev.radius_px(roi["mpp"]))
    centers, scores = centers[keep], scores[keep]
    det = pd.DataFrame({"rank": np.arange(len(centers)), "cx": centers[:, 0], "cy": centers[:, 1], "score": scores})
    shifted = det.assign(cx=det["cx"] + pad, cy=det["cy"] + pad)
    det = det.assign(od=C.cm.score_detections(shifted, roi["od_padded"], window=C.prod.OD_WINDOW)["od"].to_numpy())
    return det, dict(n_peaks=n_peaks, n_detections=len(det), n_blanked_px=int(n_blanked))


def sha1(t: np.ndarray) -> str:
    """
        SHA-1 of a template's float32 bytes.

        t (np.ndarray): the template.

        Returns str: the hex digest.
    """
    return hashlib.sha1(np.ascontiguousarray(t, dtype=np.float32).tobytes()).hexdigest()


def roi_rows(args):
    """
        Every arm for this ROI's follow-up clicks; asserts hem_mask_fill_1um parity on its joint clicks first.

        args (tuple): (roi_table row as dict, click rows, recorded fillgeom per_run rows for this ROI, recorded click rows).

        Returns list[dict]: per_run rows.
    """
    meta, clicks, ref_run, ref_clicks = args
    roi = C.load_roi(meta["path"])
    gt, mpp = roi["gt"], roi["mpp"]
    rows, cache = [], {}

    def emit(c, arm, det, info, extra):
        for rank_key in P.RANK_KEYS:
            sc = C.score(C.rank_by(det, rank_key), gt, int(c["ann_id"]), (float(c["cx"]), float(c["cy"])), mpp)
            rows.append(dict(file_name=meta["file_name"], domain=meta["domain"], subset=meta["subset"], image_id=int(meta["image_id"]), click_set=c["click_set"], seed_index=int(c["seed_index"]), seed_ann_id=int(c["ann_id"]), click_cx=float(c["cx"]), click_cy=float(c["cy"]), mpp=mpp, rank_key=rank_key, arm=arm, **extra, **info, **sc))

    for c in clicks:
        ann = gt[gt.ann_id == c["ann_id"]].iloc[0]
        cx, cy = float(ann.cx), float(ann.cy)
        T = C.tm.read_padded_patch(roi["hem"], cx, cy, W)
        if c["click_set"] == "joint":
            hmf = hem_mask_fill_1um(T, mpp)
            rc = ref_clicks[ref_clicks.seed_ann_id == c["ann_id"]].iloc[0]
            assert sha1(hmf) == rc.sha_hem_mask_fill_1um, (meta["file_name"], c["ann_id"], "hem_mask_fill_1um template differs from the recorded one")
            det, info = search_template(roi, hmf, cx, cy)
            for rank_key in P.RANK_KEYS:
                sc = C.score(C.rank_by(det, rank_key), gt, int(c["ann_id"]), (cx, cy), mpp)
                rr = ref_run[(ref_run.seed_ann_id == c["ann_id"]) & (ref_run.condition == "hem_mask_fill_1um") & (ref_run.rank_key == rank_key)].iloc[0]
                assert all(sc[f"tp_at_{k}"] == int(rr[f"tp_at_{k}"]) for k in C.BUDGETS) and info["n_peaks"] == int(rr.n_peaks), (meta["file_name"], c["ann_id"], rank_key, "hem_mask_fill_1um TP differs from the recorded run")
            emit(c, "hem_mask_fill_1um", det, info, dict(base_size=W, status="refined", n_filled=int((hmf != T).sum())))
        else:
            for s in P.arm_specs(roi, meta, ann):
                if not s["arm"].endswith("_51"):
                    continue
                key = (int(s["base_size"]), int(round(s["tx"])), int(round(s["ty"])))
                if key not in cache:
                    cache[key] = C.search(roi, *key[:1], float(s["tx"]), float(s["ty"]))
                emit(c, s["arm"], *cache[key], dict(base_size=int(s["base_size"]), status=s["status"], n_filled=0))
        sf, n_filled = sal_fill_51(roi["hem"], ann, mpp)
        det, info = search_template(roi, sf, cx, cy)
        emit(c, "sal_fill_51", det, info, dict(base_size=W, status="refined" if n_filled else "untouched", n_filled=n_filled))
    print(meta["file_name"], len(clicks), "clicks", flush=True)
    return rows


def main():
    """
        Build the click sets, run every ROI with a process pool, write the per-run CSV.

        Returns None.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--rois", default=None)
    args = ap.parse_args()
    rois = C.roi_table()
    if args.rois:
        rois = rois[rois.file_name.isin(args.rois.split(","))]
    geo = pd.read_csv(C.OUT_DIR / "geometry_midog_box.csv")
    pick = geo[((geo.prod_status != "refined") | (geo.best_in_gray < 0.5)) & (geo.sal_status != "unreadable")]
    diff = pd.DataFrame(dict(file_name=pick.file_name, click_set=np.where(pick.prod_status != "refined", "gray_refused", "different_object"), seed_index=-1, ann_id=pick.ann_id, cx=pick.cx, cy=pick.cy))
    jc = pd.read_csv(C.RESULTS / "precision_at_k_49roi_6seed_fillgeom_postD11_clicks.csv")
    joint = pd.DataFrame(dict(file_name=jc.file_name, click_set="joint", seed_index=jc.seed_index, ann_id=jc.seed_ann_id, cx=jc.click_cx, cy=jc.click_cy))
    clicks = pd.concat([joint, diff], ignore_index=True)
    ref_run = pd.read_csv(C.RESULTS / "precision_at_k_49roi_6seed_fillgeom_postD11_per_run.csv")
    jobs = [(m, clicks[clicks.file_name == m["file_name"]].to_dict("records"), ref_run[ref_run.file_name == m["file_name"]], jc[jc.file_name == m["file_name"]]) for m in rois.to_dict("records")]
    jobs = sorted([j for j in jobs if j[1]], key=lambda j: -len(j[1]))
    with Pool(args.workers) as p:
        results = p.map(roi_rows, jobs, chunksize=1)
    out = pd.DataFrame([r for rs in results for r in rs])
    suffix = "" if not args.rois else "_subset"
    out.to_csv(C.OUT_DIR / f"precision_followup_per_run{suffix}.csv", index=False)
    print(f"{len(out)} rows, {out.groupby(['file_name', 'seed_ann_id']).ngroups} clicks -> {C.OUT_DIR}")


if __name__ == "__main__":
    main()
