"""
    Geometry pass, no matching: how the saliency rule and production's gray Otsu gate behave on every border-filtered
    agreement-pool mitotic annotation of the 49 validation ROIs.

    Part A, the MIDOG box (every MIDOG++ annotation is a 50 x 50 box on the click, so this is production's own click path):
    refusals, which object each method picks, template sizes, the notes' two guards.

    Part B, simulated drawn boxes (the notes' Stage 1b protocol): a reference extent padded by 15 / 50 / 100 % and jittered
    by up to +/-20 % of the box size, JITTERS draws each, sides clipped to the D12 limits (18-72 px). Every method sees the
    same box. The reference is either the saliency mask's bbox or the gray component's bbox at the MIDOG box; each favours
    its own method, so both are run.

    Writes out/geometry_midog_box.csv (one row per annotation) and out/geometry_boxes.csv.gz (one row per annotation x
    reference x padding x jitter). Usage: geometry.py [--workers N] [--rois A,B]
"""

import argparse
import math
from multiprocessing import Pool

import numpy as np
import pandas as pd

import common as C

PADS = (0.15, 0.50, 1.00)
JITTER = 0.20
JITTERS = 5
MIN_SIDE, MAX_SIDE = 18.0, 72.0


def frac_inside(a: np.ndarray, b: np.ndarray) -> float:
    """
        Share of a's pixels that are also in b.

        a (np.ndarray): boolean mask.
        b (np.ndarray): boolean mask.

        Returns float: |a & b| / |a|, nan when a is empty.
    """
    n = int(a.sum())
    return float((a & b).sum() / n) if n else float("nan")


def square_iou(base_size: int, tx: float, ty: float, ref) -> float:
    """
        IoU of a template's pixel footprint with a reference pixel bbox, both in image coordinates.

        base_size (int): template side.
        tx (float): template centre x.
        ty (float): template centre y.
        ref (tuple): (y0, y1, x0, x1) half-open image-coordinate bbox.

        Returns float: the IoU.
    """
    h = base_size // 2
    ix, iy = int(round(tx)), int(round(ty))
    a = (iy - h, iy + h + 1, ix - h, ix + h + 1)
    iy0, iy1, ix0, ix1 = max(a[0], ref[0]), min(a[1], ref[1]), max(a[2], ref[2]), min(a[3], ref[3])
    inter = max(0, iy1 - iy0) * max(0, ix1 - ix0)
    union = base_size * base_size + (ref[1] - ref[0]) * (ref[3] - ref[2]) - inter
    return inter / union


def simulated_box(ref, pad: float, rng):
    """
        One simulated drawn box around a reference extent, in the edge coordinates box_frame takes.

        ref (tuple): (y0, y1, x0, x1) half-open image-coordinate bbox.
        pad (float): fractional padding of each side.
        rng (np.random.Generator): jitter stream.

        Returns tuple: (x0, y0, w, h, clipped) edge-coordinate box and whether a side was clipped to 18-72 px.
    """
    y0, y1, x0, x1 = ref
    w, h = (x1 - x0) * (1 + pad), (y1 - y0) * (1 + pad)
    cx, cy = (x0 + x1) / 2.0 + rng.uniform(-JITTER, JITTER) * w, (y0 + y1) / 2.0 + rng.uniform(-JITTER, JITTER) * h
    wc, hc = min(max(w, MIN_SIDE), MAX_SIDE), min(max(h, MIN_SIDE), MAX_SIDE)
    return cx - wc / 2, cy - hc / 2, wc, hc, bool(wc != w or hc != h)


def roi_rows(args):
    """
        Parts A and B for one ROI.

        args (tuple): (roi table row as dict,).

        Returns tuple[list[dict], list[dict]]: Part A rows and Part B rows.
    """
    (meta,) = args
    roi = C.load_roi(meta["path"])
    hem, gray, mpp, shape = roi["hem"], roi["gray"], roi["mpp"], roi["hem"].shape
    gt = roi["gt"]
    pool, flagged = C.ss.agreement_pool(gt[gt.category_id == C.ds.MITOTIC])
    pool = C.ss.border_filter(pool, C.HALF_PATCH, shape)
    rows_a, rows_b = [], []
    for _, a in pool.iterrows():
        cx, cy, ann = float(a.cx), float(a.cy), int(a.ann_id)
        ix, iy = int(round(cx)), int(round(cy))
        base = dict(file_name=meta["file_name"], domain=meta["domain"], subset=meta["subset"], image_id=int(meta["image_id"]), ann_id=ann, cx=cx, cy=cy, mpp=mpp, contested_tier=bool(flagged))

        # Part A: the MIDOG box
        x_tl, y_tl, w, h = C.midog_box(a)
        window = C.af.gate_window_side(w, h)
        prodA = C.production_spec(gray, cx, cy, window, shape)
        hemA = C.ss.tightened_template_box(hem, cx, cy, otsu_window=51)
        hem_ok = hemA is not None and C.readable(shape, hemA[1], hemA[2])
        gpatch = C.tm.read_padded_patch(gray, cx, cy, 51)
        gcomp = C.gray_component(gpatch)
        salA, sal = C.saliency_spec(hem, x_tl, y_tl, w, h, mpp, window, dilate=False)
        salA1, _ = C.saliency_spec(hem, x_tl, y_tl, w, h, mpp, window, dilate=True)
        row = dict(base, window=window, prod_status=prodA["status"], prod_base=prodA["base_size"], prod_tx=prodA["tx"], prod_ty=prodA["ty"], hem_gate_ok=hem_ok, joint_gate_ok=bool(prodA["status"] == "refined" and hem_ok), sal_status=salA["status"], sal_base=salA["base_size"], sal_tx=salA["tx"], sal_ty=salA["ty"], sal1_status=salA1["status"], sal1_base=salA1["base_size"], sal1_tx=salA1["tx"], sal1_ty=salA1["ty"], sal_n_comp=sal["n_comp"], sal_n_above_floor=sal["n_above_floor"], sal_n_kept=sal["n_kept"], sal_below_floor=sal["below_floor"], sal_guard_small=sal["guard_small"], sal_guard_edge=sal["guard_edge"], sal_mask_px=int(sal["mask"].sum()) if sal["mask"] is not None else 0)
        if gcomp is not None and sal["mask"] is not None:
            # gray_component works on the 51 px window; lift it into the 73 px patch the saliency mask lives in
            gmask = np.zeros((C.PATCH, C.PATCH), bool)
            o = C.HALF_PATCH - 25
            gmask[o:o + 51, o:o + 51] = gcomp[0]
            row.update(best_in_gray=frac_inside(sal["best"], gmask), kept_in_gray=frac_inside(sal["mask"], gmask), gray_in_kept=frac_inside(gmask, sal["mask"]))
        row["centre_dist_um"] = float(math.hypot(salA["tx"] - prodA["tx"], salA["ty"] - prodA["ty"]) * mpp)
        rows_a.append(row)

        # Part B: simulated drawn boxes from each reference extent
        refs = {}
        if sal["mask"] is not None:
            y0, y1, x0, x1 = C.mask_bbox(sal["mask"])
            refs["sal"] = (iy - C.HALF_PATCH + y0, iy - C.HALF_PATCH + y1, ix - C.HALF_PATCH + x0, ix - C.HALF_PATCH + x1)
        if gcomp is not None:
            y0, y1, x0, x1 = gcomp[1]
            refs["gray"] = (iy - 25 + y0, iy - 25 + y1, ix - 25 + x0, ix - 25 + x1)
        for ref_name, ref in refs.items():
            ref_side = max(ref[1] - ref[0], ref[3] - ref[2])
            ref_c = ((ref[2] + ref[3] - 1) / 2.0, (ref[0] + ref[1] - 1) / 2.0)
            for pi, pad in enumerate(PADS):
                rng = np.random.default_rng([int(meta["image_id"]), ann, pi, 0 if ref_name == "sal" else 1])
                for j in range(JITTERS):
                    bx0, by0, bw, bh, clipped = simulated_box(ref, pad, rng)
                    fr = C.af.box_frame(meta["file_name"], bx0, by0, bw, bh, 0).iloc[0]
                    x_tl, y_tl, width, height = float(fr.x_tl), float(fr.y_tl), float(fr.width), float(fr.height)
                    bcx, bcy = C.af.box_center(x_tl, width), C.af.box_center(y_tl, height)
                    if not C.readable(shape, bcx, bcy):
                        continue
                    window = C.af.gate_window_side(width, height)
                    specs = {"prod": C.production_spec(gray, bcx, bcy, window, shape)}
                    specs["sal"], salB = C.saliency_spec(hem, x_tl, y_tl, width, height, mpp, window, dilate=False)
                    specs["sal1"], _ = C.saliency_spec(hem, x_tl, y_tl, width, height, mpp, window, dilate=True)
                    for method, s in specs.items():
                        rows_b.append(dict(base, ref=ref_name, pad=pad, jitter=j, clipped=clipped, box_w=width + 1, box_h=height + 1, window=window, ref_side=ref_side, method=method, status=s["status"], base_size=s["base_size"], tx=s["tx"], ty=s["ty"], iou_ref=square_iou(s["base_size"], s["tx"], s["ty"], ref), centre_err_um=float(math.hypot(s["tx"] - ref_c[0], s["ty"] - ref_c[1]) * mpp), guard_small=salB["guard_small"] if method != "prod" else None, guard_edge=salB["guard_edge"] if method != "prod" else None))
    print(meta["file_name"], len(rows_a), "annotations", flush=True)
    return rows_a, rows_b


def main():
    """
        Run both parts over the ROIs with a process pool and write the two CSVs.

        Returns None.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--rois", default=None)
    args = ap.parse_args()
    rois = C.roi_table()
    if args.rois:
        rois = rois[rois.file_name.isin(args.rois.split(","))]
    C.OUT_DIR.mkdir(exist_ok=True)
    with Pool(args.workers) as p:
        results = p.map(roi_rows, [(r,) for r in rois.to_dict("records")], chunksize=1)
    a = pd.DataFrame([r for ra, _ in results for r in ra])
    b = pd.DataFrame([r for _, rb in results for r in rb])
    suffix = "" if not args.rois else "_subset"
    a.to_csv(C.OUT_DIR / f"geometry_midog_box{suffix}.csv", index=False)
    b.to_csv(C.OUT_DIR / f"geometry_boxes{suffix}.csv.gz", index=False)  # ~33 MB uncompressed
    print(f"{len(a)} annotations, {len(b)} box rows -> {C.OUT_DIR}")


if __name__ == "__main__":
    main()
