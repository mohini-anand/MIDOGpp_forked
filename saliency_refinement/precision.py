"""
    Precision@10/20/30 through production's own search for the saliency rule against production's gray Otsu gate, on
    the same seed and the same box for every arm.

    Click sets:
      joint    the 294 post-D11 reference clicks (49 ROIs x 6 seeds, results/..._fillgeom_postD11_clicks.csv), drawn under a
               gate that required gray AND hem-under-click to accept -- so they exclude the clicks production refuses.
      refusal  up to 2 clicks per ROI that gate excluded (gray or hem-under-click refused), from out/geometry_midog_box.csv,
               drawn with default_rng([4242, image_id]). In production a gray refusal falls back to the untightened box.

    Arms, all cut as production cuts a template (a raw hematoxylin OD square) and searched by production's find_and_suppress:
      MIDOG box (50 x 50 on the click, production's click path):  default_51, prod_51 (gray gate), sal_51, sal1_51 (1 um)
      simulated drawn box, jitter draw 0 of geometry.py's stream, reference "sal" and "gray", padding 15/50/100 %:
               prod_box (D12: gray gate in the box's window, fallback the box), sal_box, sal1_box

    Writes out/precision_per_run.csv, one row per click x arm x rank key. Usage: precision.py [--workers N] [--rois A,B]
"""

import argparse
from multiprocessing import Pool

import numpy as np
import pandas as pd

import common as C
from geometry import PADS, simulated_box

REFUSAL_PER_ROI = 2
RANK_KEYS = ("tm_score", "chromatin_od")


def click_sets(rois: pd.DataFrame) -> pd.DataFrame:
    """
        The joint and refusal click sets.

        rois (pd.DataFrame): roi_table rows to include.

        Returns pd.DataFrame: file_name, click_set, seed_index, ann_id, cx, cy.
    """
    joint = pd.read_csv(C.RESULTS / "precision_at_k_49roi_6seed_fillgeom_postD11_clicks.csv")
    joint = joint[joint.file_name.isin(rois.file_name)]
    out = [pd.DataFrame(dict(file_name=joint.file_name, click_set="joint", seed_index=joint.seed_index, ann_id=joint.seed_ann_id, cx=joint.click_cx, cy=joint.click_cy))]
    geo = pd.read_csv(C.OUT_DIR / "geometry_midog_box.csv")
    for fn, g in geo[(~geo.joint_gate_ok) & geo.file_name.isin(rois.file_name)].groupby("file_name", sort=True):
        rng = np.random.default_rng([4242, int(g.image_id.iat[0])])
        pick = g.iloc[rng.permutation(len(g))[:REFUSAL_PER_ROI]]
        out.append(pd.DataFrame(dict(file_name=fn, click_set="refusal", seed_index=-1, ann_id=pick.ann_id, cx=pick.cx, cy=pick.cy)))
    return pd.concat(out, ignore_index=True)


def arm_specs(roi: dict, meta: dict, ann: pd.Series) -> list[dict]:
    """
        Every arm's template spec for one click.

        roi (dict): load_roi output.
        meta (dict): the ROI's roi_table row.
        ann (pd.Series): the click's annotation row.

        Returns list[dict]: one dict per arm with arm, ref, pad, base_size, tx, ty, status, window.
    """
    hem, gray, mpp, shape = roi["hem"], roi["gray"], roi["mpp"], roi["hem"].shape
    cx, cy, ann_id = float(ann.cx), float(ann.cy), int(ann.ann_id)
    ix, iy = int(round(cx)), int(round(cy))
    x_tl, y_tl, w, h = C.midog_box(ann)
    specs = [dict(arm="default_51", ref="", pad=np.nan, base_size=51, tx=cx, ty=cy, status="untightened", window=51)]
    specs.append(dict(arm="prod_51", ref="", pad=np.nan, window=51, **C.production_spec(gray, cx, cy, 51, shape)))
    s, sal = C.saliency_spec(hem, x_tl, y_tl, w, h, mpp, 51, dilate=False)
    specs.append(dict(arm="sal_51", ref="", pad=np.nan, window=51, **s))
    s1, _ = C.saliency_spec(hem, x_tl, y_tl, w, h, mpp, 51, dilate=True)
    specs.append(dict(arm="sal1_51", ref="", pad=np.nan, window=51, **s1))

    refs = {}
    if sal["mask"] is not None:
        y0, y1, x0, x1 = C.mask_bbox(sal["mask"])
        refs["sal"] = (iy - C.HALF_PATCH + y0, iy - C.HALF_PATCH + y1, ix - C.HALF_PATCH + x0, ix - C.HALF_PATCH + x1)
    gcomp = C.gray_component(C.tm.read_padded_patch(gray, cx, cy, 51))
    if gcomp is not None:
        y0, y1, x0, x1 = gcomp[1]
        refs["gray"] = (iy - 25 + y0, iy - 25 + y1, ix - 25 + x0, ix - 25 + x1)
    for ref_name, ref in refs.items():
        for pi, pad in enumerate(PADS):
            rng = np.random.default_rng([int(meta["image_id"]), ann_id, pi, 0 if ref_name == "sal" else 1])
            bx0, by0, bw, bh, _ = simulated_box(ref, pad, rng)  # jitter draw 0, the same box geometry.py records as jitter 0
            fr = C.af.box_frame(meta["file_name"], bx0, by0, bw, bh, 0).iloc[0]
            bx_tl, by_tl, bwidth, bheight = float(fr.x_tl), float(fr.y_tl), float(fr.width), float(fr.height)
            bcx, bcy = C.af.box_center(bx_tl, bwidth), C.af.box_center(by_tl, bheight)
            if not C.readable(shape, bcx, bcy):
                continue
            window = C.af.gate_window_side(bwidth, bheight)
            specs.append(dict(arm="prod_box", ref=ref_name, pad=pad, window=window, **C.production_spec(gray, bcx, bcy, window, shape)))
            sb, _ = C.saliency_spec(hem, bx_tl, by_tl, bwidth, bheight, mpp, window, dilate=False)
            specs.append(dict(arm="sal_box", ref=ref_name, pad=pad, window=window, **sb))
            sb1, _ = C.saliency_spec(hem, bx_tl, by_tl, bwidth, bheight, mpp, window, dilate=True)
            specs.append(dict(arm="sal1_box", ref=ref_name, pad=pad, window=window, **sb1))
    return specs


def roi_rows(args):
    """
        Every click x arm x rank key on one ROI, searching each distinct template once.

        args (tuple): (roi_table row as dict, that ROI's click rows as a list of dicts).

        Returns list[dict]: per_run rows.
    """
    meta, clicks = args
    roi = C.load_roi(meta["path"])
    gt = roi["gt"]
    cache, rows = {}, []
    for c in clicks:
        ann = gt[gt.ann_id == c["ann_id"]].iloc[0]
        assert (float(ann.cx), float(ann.cy)) == (float(c["cx"]), float(c["cy"]))
        for s in arm_specs(roi, meta, ann):
            key = (int(s["base_size"]), int(round(s["tx"])), int(round(s["ty"])))
            if key not in cache:
                cache[key] = C.search(roi, int(s["base_size"]), float(s["tx"]), float(s["ty"]))
            det, info = cache[key]
            for rank_key in RANK_KEYS:
                sc = C.score(C.rank_by(det, rank_key), gt, int(c["ann_id"]), (float(c["cx"]), float(c["cy"])), roi["mpp"])
                rows.append(dict(file_name=meta["file_name"], domain=meta["domain"], subset=meta["subset"], image_id=int(meta["image_id"]), click_set=c["click_set"], seed_index=int(c["seed_index"]), seed_ann_id=int(c["ann_id"]), click_cx=float(c["cx"]), click_cy=float(c["cy"]), mpp=roi["mpp"], rank_key=rank_key, template_key=f"{key[0]}@{key[1]},{key[2]}", n_peaks=int(info["n_peaks"]), n_detections=int(info["n_detections"]), n_blanked_px=int(info["n_blanked_px"]), **s, **sc))
    print(meta["file_name"], len(clicks), "clicks", len(cache), "searches", flush=True)
    return rows


def main():
    """
        Run every ROI with a process pool and write the per-run CSV.

        Returns None.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--rois", default=None)
    args = ap.parse_args()
    rois = C.roi_table()
    if args.rois:
        rois = rois[rois.file_name.isin(args.rois.split(","))]
    clicks = click_sets(rois)
    jobs = [(m, clicks[clicks.file_name == m["file_name"]].to_dict("records")) for m in rois.to_dict("records")]
    jobs.sort(key=lambda j: -len(j[1]))
    with Pool(args.workers) as p:
        results = p.map(roi_rows, jobs, chunksize=1)
    out = pd.DataFrame([r for rs in results for r in rs])
    suffix = "" if not args.rois else "_subset"
    out.to_csv(C.OUT_DIR / f"precision_per_run{suffix}.csv", index=False)
    print(f"{len(out)} rows, {out.groupby(['file_name', 'click_set', 'seed_ann_id']).ngroups} clicks -> {C.OUT_DIR}")


if __name__ == "__main__":
    main()
