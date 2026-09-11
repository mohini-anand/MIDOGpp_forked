"""F7 stage 1: reproduction gate, then per-candidate physical size.

See `f7_size_prune/PREREGISTRATION.md`, including the 2026-09-08 instrument amendment (§5a)
that replaced the window-local Otsu with the tiled-Otsu segmentation `nucleus_blobs` uses.

Nothing upstream of the prune is recomputed. The deep pool comes from
`results/tm_recall_workload_pool_z.npz`, which `recall_workload_ledger.py` committed, so the
padded ROI, TM_CCOEFF, hematoxylin_od, largest-CC template sizing, the 5.0 um NMS radius, the
self-hit filter, DEEP_FLOOR_Z and the seed draw are frozen by construction rather than by
re-declaration. The only new computation is a per-candidate size; the only new operation is a
deletion, applied downstream in the notebook.

Outputs (all under f7_size_prune/results/):
  f7_repro_gate.csv        one row per cell: does re-matching the UNPRUNED cached pool
                           reproduce the committed TP ledger exactly?
  f7_candidate_sizes.npz   per cell, every candidate's component-equivalent diameter in um,
                           bbox max side in um, component area in px, and on_nucleus flag,
                           all in pool (score-descending) order.

Run from the repository root:  python f7_size_prune/f7_size_features.py
"""

from __future__ import annotations

import gc
import os
import sys
import time

import cv2
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from skimage.color import rgb2hed

from midog_utils import baselines as bl
from midog_utils import dataset as ds
from midog_utils import evaluate as ev

IMAGES_DIR = "images/extra_valid"
SEEDS = tuple(range(5))
OTSU_TILE = 512                     # baselines.nucleus_blobs' own default

POOL_NPZ = "results/tm_recall_workload_pool_z.npz"
CELLS = "results/tm_recall_workload_cells.csv"
LEDGER = "results/tm_recall_workload_tp_ledger.csv"

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
OUT_GATE = os.path.join(OUT_DIR, "f7_repro_gate.csv")
OUT_SIZES = os.path.join(OUT_DIR, "f7_candidate_sizes.npz")


def hematoxylin_banded(rgb: np.ndarray, band: int = 512) -> np.ndarray:
    """`baselines.to_hematoxylin`, computed in row bands to bound peak memory.

    Identical output, different memory profile. `skimage.color.rgb2hed` promotes to float64
    internally, so calling it on a whole 39-megapixel ROI allocates ~0.9 GB transiently on top
    of everything else -- which OOM-killed a first attempt at this run on 245.tiff. Banding
    caps the transient at `band` rows. The 0.5/99.5 percentile rescale is taken over the
    **whole** image afterwards, exactly as `baselines.to_hematoxylin` does, so the result is
    the same array.
    """
    h, w = rgb.shape[:2]
    out = np.empty((h, w), dtype=np.float32)
    for y0 in range(0, h, band):
        y1 = min(y0 + band, h)
        out[y0:y1] = rgb2hed(rgb[y0:y1].astype(np.float32) / 255.0)[:, :, 0].astype(np.float32)
    lo, hi = np.percentile(out, [0.5, 99.5])
    if hi <= lo:
        return np.zeros_like(out)
    np.subtract(out, lo, out=out)
    np.divide(out, hi - lo, out=out)
    np.clip(out, 0.0, 1.0, out=out)
    np.multiply(out, 255.0, out=out)
    return out


def segment_roi(rgb: np.ndarray, tile: int = OTSU_TILE):
    """Label every nucleus in the ROI, using `nucleus_blobs`' own segmentation.

    `baselines.nucleus_blobs` is `tissue_mask` -> `to_hematoxylin` -> `_tiled_otsu_threshold`
    -> connected components. F7 reuses the segmentation and *not* the `min_area`/`max_area`
    gate or the darkness ranking, because those are `nucleus_blobs`' knobs and F7 has one.

    Why this and not an Otsu over a 51 px window around each candidate (the original section
    5): `cv2.normalize(NORM_MINMAX)` followed by Otsu forces a roughly 50/50 split of *any*
    window regardless of its content, so the component under the candidate is about half the
    window and its bbox pins to the window edge. That is a measurement of the window, not of
    the object. See PREREGISTRATION.md section 5a for the numbers.
    """
    mask = bl.tissue_mask(rgb)
    h_chan = hematoxylin_banded(rgb)
    thr = bl._tiled_otsu_threshold(h_chan, mask, tile)
    binary = ((h_chan > thr) & mask).astype(np.uint8)
    del h_chan, thr, mask
    gc.collect()
    n, lab, stats, _ = cv2.connectedComponentsWithStats(binary, connectivity=8)
    del binary
    gc.collect()
    return lab, stats[:, cv2.CC_STAT_AREA].astype(np.float64), \
        stats[:, cv2.CC_STAT_WIDTH].astype(np.float64), \
        stats[:, cv2.CC_STAT_HEIGHT].astype(np.float64), n - 1


def main():
    t_all = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    images, annotations = ds.load_annotations()
    ds.check_invariants(annotations)
    meta = images.set_index("file_name")[["image_id", "tumor_type"]]

    pool_npz = np.load(POOL_NPZ)
    cells = pd.read_csv(CELLS)
    ledger = pd.read_csv(LEDGER)

    files = sorted(f for f in os.listdir(IMAGES_DIR) if f.endswith(".tiff"))
    print(f"{len(files)} ROIs x {len(SEEDS)} seeds", flush=True)

    gate_rows, size_out = [], {}

    for fn in files:
        t_roi = time.time()
        rgb = ds.load_roi(f"{IMAGES_DIR}/{fn}")
        lab, area_px, w_px, h_px, n_comp = segment_roi(rgb)
        H, W = lab.shape
        del rgb
        gc.collect()
        gt_all = ds.image_annotations(annotations, fn)
        print(f"{fn}: segmented, {n_comp} components, {time.time()-t_roi:.0f}s", flush=True)

        for si in SEEDS:
            sel = cells[(cells.file_name == fn) & (cells.seed_index == si)]
            assert len(sel) == 1, f"{fn} s{si}: expected one cell row, got {len(sel)}"
            row = sel.iloc[0]
            seed_ann_id = int(row["seed_ann_id"])
            match_radius = float(row["match_radius_px"])   # from the CSV, never recomputed
            mpp = float(row["mpp"])

            cx = pool_npz[f"{fn}|{si}|cx"]
            cy = pool_npz[f"{fn}|{si}|cy"]
            assert len(cx) == int(row["n_pool"]), \
                f"{fn} s{si}: npz has {len(cx)} candidates, cells CSV says {row['n_pool']}"

            # ---- reproduction gate: re-match the UNPRUNED pool, in stored order -------------
            # The npz is score-descending and carries no `score` column: array order IS the
            # rank, so it is used as-is and never re-sorted.
            gt_eval = gt_all[gt_all["ann_id"] != seed_ann_id].reset_index(drop=True)
            det, _ = ev.bucket_detections(pd.DataFrame({"cx": cx, "cy": cy}),
                                          gt_eval, match_radius)
            tp = det[det["bucket"] == ev.HUMAN_CORRECT_LABEL]
            # `rank` in the ledger is 1-based -- recall_workload_ledger.py stores ranks0 + 1,
            # "the list LENGTH a reader would have to be given to receive it".
            got = np.sort(tp.index.to_numpy()) + 1
            want = np.sort(ledger[(ledger.file_name == fn)
                                  & (ledger.seed_index == si)]["rank"].to_numpy())
            exact = bool(len(got) == len(want) and np.array_equal(got, want))

            n_mit = int((gt_eval["category_id"] == ds.MITOTIC).sum())
            gate_rows.append(dict(
                file_name=fn, tumor_type=meta.loc[fn, "tumor_type"], seed_index=si,
                seed_ann_id=seed_ann_id, mpp=mpp, match_radius_px=match_radius,
                n_pool=len(cx), n_gt_mitotic=n_mit, n_gt_mitotic_csv=int(row["n_gt_mitotic"]),
                n_tp_recomputed=len(got), n_tp_ledger=len(want), ranks_exact=exact,
                n_components=n_comp))

            # ---- the one knob: per-candidate physical size ----------------------------------
            ix = np.round(cx).astype(np.int64).clip(0, W - 1)
            iy = np.round(cy).astype(np.int64).clip(0, H - 1)
            ids = lab[iy, ix]
            on = ids > 0
            eqd = np.full(len(cx), np.nan, dtype=np.float32)
            side = np.full(len(cx), np.nan, dtype=np.float32)
            apx = np.full(len(cx), np.nan, dtype=np.float32)
            a = area_px[ids[on]]
            eqd[on] = (2.0 * np.sqrt(a / np.pi) * mpp).astype(np.float32)
            side[on] = (np.maximum(w_px[ids[on]], h_px[ids[on]]) * mpp).astype(np.float32)
            apx[on] = a.astype(np.float32)

            size_out[f"{fn}|{si}|eqd_um"] = eqd
            size_out[f"{fn}|{si}|side_um"] = side
            size_out[f"{fn}|{si}|area_px"] = apx
            size_out[f"{fn}|{si}|on_nucleus"] = on

            gate_rows[-1]["frac_on_nucleus"] = float(on.mean())
            gate_rows[-1]["n_size_undefined"] = int((~on).sum())
            print(f"  {fn} s{si}: {len(cx)} cand, gate exact={exact} "
                  f"({len(got)}/{len(want)} TP), on_nucleus={on.mean():.3f}", flush=True)

        del lab, area_px, w_px, h_px
        gc.collect()
        print(f"{fn} done in {time.time()-t_roi:.0f}s", flush=True)

    gate = pd.DataFrame(gate_rows)
    gate.to_csv(OUT_GATE, index=False)
    np.savez_compressed(OUT_SIZES, **size_out)
    print(f"\nwrote {OUT_GATE} and {OUT_SIZES} in {time.time()-t_all:.0f}s")
    print(f"reproduction gate: {int(gate.ranks_exact.sum())}/{len(gate)} cells exact")
    print(f"n_gt_mitotic agrees with cells CSV: "
          f"{int((gate.n_gt_mitotic == gate.n_gt_mitotic_csv).sum())}/{len(gate)}")
    if not bool(gate.ranks_exact.all()):
        print("GATE FAILED -- downstream numbers are void", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
