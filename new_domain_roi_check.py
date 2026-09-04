"""Pixel-level viability verification for newly downloaded MIDOG++ ROIs."""
import csv, gc, os, sys, traceback
from pathlib import Path
import numpy as np
from midog_utils import dataset, seed_selection, baselines, invariants

CANDS = [("459", "canine soft tissue sarcoma"), ("094", "human breast cancer"),
         ("548", "human melanoma"), ("402", "human neuroendocrine tumor")]
OUT = Path("results/new_domain_rois.csv")
FIELDS = ["file_name", "image_id", "domain", "json_tumor_type", "scanner", "lab", "species",
          "width_px", "height_px", "mpp", "roi_area_mm2", "n_mitotic", "n_lookalike",
          "n_unanimous_mitotic", "seed_pool_size", "agreement_flagged",
          "tm_n_excluded_all", "tm_n_excluded_mitotic", "file_size_bytes", "checks_passed", "notes"]

images, anns = dataset.load_annotations()
dataset.check_invariants(anns)
meta = dataset.load_slide_metadata()
FIDS = {"459": "40283401", "094": "40282360", "548": "40283659", "402": "40283125"}

writer = None
fh = OUT.open("w", newline="")
writer = csv.DictWriter(fh, fieldnames=FIELDS); writer.writeheader(); fh.flush()

for nnn, want_domain in CANDS:
    fn = f"{nnn}.tiff"; path = Path("images") / fn
    checks, notes = [], []
    try:
        im = images[images["file_name"] == fn].iloc[0]
        a = anns[anns["file_name"] == fn]
        mit = a[a["category_id"] == 1]; look = a[a["category_id"] == 2]

        # check 7: tumor_type matches, and file was not already present pre-task
        c7 = (im["tumor_type"] == want_domain)
        checks.append(("7_domain_matches", c7))
        if not c7: notes.append(f"json tumor_type={im['tumor_type']!r} != {want_domain!r}")

        # check 1
        rgb = dataset.load_roi(path)
        c1 = (rgb.ndim == 3 and rgb.shape[2] == 3 and rgb.dtype == np.uint8)
        checks.append(("1_load_roi_hxwx3_uint8", c1))
        # truncation guard: pixel dims must match this image's own JSON dims
        cdim = (rgb.shape[0] == im["height"] and rgb.shape[1] == im["width"])
        checks.append(("dims_match_json", cdim))
        if not cdim: notes.append(f"pixels {rgb.shape[:2]} != json ({im['height']},{im['width']})")

        # check 2
        mpp = dataset.roi_mpp(path)
        c2 = 0.22 <= mpp <= 0.26
        checks.append(("2_mpp_sane", c2))
        if not c2: notes.append(f"mpp={mpp:.4f} outside 0.22-0.26")

        # check 3
        try:
            area = dataset.check_roi_scale(path, rgb.shape); c3 = True
        except AssertionError as e:
            area = dataset.roi_area_mm2(path, rgb.shape); c3 = False; notes.append(str(e))
        checks.append(("3_roi_scale_1.9_2.1mm2", c3))

        # check 4
        c4 = len(mit) >= 15
        checks.append(("4_n_mitotic_ge_15", c4))

        # check 5 -- on real pixel shape
        pool, flagged = seed_selection.agreement_pool(mit)
        seeds = seed_selection.border_filter(pool, 36, rgb.shape)
        c5 = len(seeds) >= 5
        checks.append(("5_seed_pool_ge_5", c5))
        if flagged: notes.append("agreement_pool fell back to contested tier")

        # check 6 -- full annotation frame so look-alike exclusions are visible
        mask = baselines.tissue_mask(rgb)
        res = invariants.check_tissue_mask_covers_gt(mask, a, label=fn)
        c6 = bool(res["passed"])
        checks.append(("6_tissue_mask_covers_mitotic", c6))
        if res["n_excluded"]:
            notes.append(f"mask excludes {res['n_excluded']} ann(s) "
                         f"({res['n_excluded_strict']} mitotic): {res['excluded_ann_ids']}")

        mrow = meta[meta["image_id"] == int(nnn)]
        row = dict(
            file_name=fn, image_id=int(nnn), domain=want_domain, json_tumor_type=im["tumor_type"],
            scanner=mrow["Scanner"].iloc[0] if len(mrow) else "",
            lab=mrow["Origin"].iloc[0] if len(mrow) else "",
            species=mrow["Species"].iloc[0] if len(mrow) else "",
            width_px=rgb.shape[1], height_px=rgb.shape[0], mpp=round(float(mpp), 5),
            roi_area_mm2=round(float(area), 4), n_mitotic=len(mit), n_lookalike=len(look),
            n_unanimous_mitotic=int((mit["n_mitotic_votes"] == mit["n_votes"]).sum()),
            seed_pool_size=len(seeds), agreement_flagged=bool(flagged),
            tm_n_excluded_all=int(res["n_excluded"]), tm_n_excluded_mitotic=int(res["n_excluded_strict"]),
            file_size_bytes=path.stat().st_size,
            checks_passed="ALL" if all(v for _, v in checks) else
                          "FAILED:" + ",".join(k for k, v in checks if not v),
            notes="; ".join(notes))
        writer.writerow(row); fh.flush()
        print(f"--- {fn} [{want_domain}] ---")
        for k, v in checks: print(f"    {'PASS' if v else 'FAIL'}  {k}")
        for k in FIELDS: print(f"    {k}: {row[k]}")
        print()
        del rgb, mask; gc.collect()
    except Exception:
        print(f"--- {fn}: EXCEPTION ---"); traceback.print_exc()
        writer.writerow({"file_name": fn, "image_id": int(nnn), "domain": want_domain,
                         "checks_passed": "EXCEPTION", "notes": traceback.format_exc()[-300:]})
        fh.flush(); gc.collect()
fh.close()
print("wrote", OUT)
