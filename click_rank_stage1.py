"""Step 1 of the click-ranking experiment: candidates, crops, seeds.

`Research Logs/2026-09-01-one-click-retrieval-literature.md` §7. The question is whether
ranking the *same* seedless candidate set by similarity-to-the-click beats ranking it by
chromatin density -- i.e. whether the pathologist's click carries information that
`TM_CCOEFF_NORMED` (mean-intensity invariant, see the 2026-08-31 chromatin log) could not
carry.

Split into three scripts because the two things this needs do not live in one interpreter:
`midog_utils` needs numpy<2 for cv2/tifffile (the anaconda install), and the torch wheel on
this machine is built against numpy 1.x but the system python has numpy 2.x. Stage 1 (here)
and stage 3 run under anaconda python; stage 2 (`click_rank_embed.py`) runs under the system
python with torch. The handoff is .npy/.csv in the scratchpad.

Crops are cut at a fixed *micron* field (32 um), not a fixed pixel size, so a candidate is
framed identically on ROIs with different mpp. The 16 um variant is the exact centre crop of
the 32 um one, so the scale control costs no extra extraction and no extra storage.
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd

from midog_utils import baselines as bl
from midog_utils import channels as ch
from midog_utils import chromatin as cm
from midog_utils import dataset as ds
from midog_utils import evaluate as ev
from midog_utils import find_and_suppress as fs
from midog_utils import invariants as inv
from midog_utils import seed_selection as ss

ROIS = ("301.tiff", "246.tiff")   # the only two ROIs with enough of both classes
N_SEEDS = 5
CROP_UM = 32.0                    # field of view of the crop handed to an encoder
CROP_PX = 128                     # every crop is resampled to this, so all ROIs match
CHANNEL = "rgb"
OUT = os.environ.get("CLICK_RANK_DIR", os.path.expanduser("~/.cache/annotatedx_click_rank"))


def cut_crops(rgb: np.ndarray, cx: np.ndarray, cy: np.ndarray, half: int) -> tuple:
    """Crops of side ``2*half``, edge-clamped. Returns (crops, shifted_flags)."""
    h, w, _ = rgb.shape
    n = len(cx)
    out = np.empty((n, 2 * half, 2 * half, 3), np.uint8)
    shifted = np.zeros(n, bool)
    for i in range(n):
        x0 = int(round(cx[i])) - half
        y0 = int(round(cy[i])) - half
        x1, y1 = x0 + 2 * half, y0 + 2 * half
        sx = min(max(x0, 0), max(w - 2 * half, 0))
        sy = min(max(y0, 0), max(h - 2 * half, 0))
        shifted[i] = (sx != x0) or (sy != y0)
        out[i] = rgb[sy:sy + 2 * half, sx:sx + 2 * half, :3]
    return out, shifted


def main():
    os.makedirs(OUT, exist_ok=True)
    images, ann = ds.load_annotations()
    border = fs.FSConfig(channel=CHANNEL).patch_size // 2
    checks = []

    for fn in ROIS:
        t0 = time.time()
        path = f"images/{fn}"
        image_id = int(images.loc[images["file_name"] == fn, "image_id"].iloc[0])
        rgb = ds.load_roi(path)
        mpp = ds.roi_mpp(path)
        radius = ev.radius_px(mpp)
        gray = ch.to_gray_inverted(rgb)
        hem = cm.hematoxylin_od(rgb)
        gt = ds.image_annotations(ann, fn)
        gt_mit = gt[gt["category_id"] == ds.MITOTIC]

        mask = bl.tissue_mask(rgb)
        checks.append(inv.check_tissue_mask_covers_gt(mask, gt, label=fn))

        blobs = bl.nucleus_blobs(rgb, mask)
        blobs = cm.score_detections(blobs, hem)      # adds the `od` column
        half_native = int(round(CROP_UM / mpp / 2))  # 32 um half-width at this ROI's mpp
        crops, shifted = cut_crops(rgb, blobs["cx"].to_numpy(), blobs["cy"].to_numpy(),
                                   half_native)
        # one resample to a common pixel grid, so encoders see identical geometry per ROI
        if crops.shape[1] != CROP_PX:
            import cv2
            crops = np.stack([cv2.resize(c, (CROP_PX, CROP_PX),
                                         interpolation=cv2.INTER_AREA) for c in crops])
        blobs = blobs.assign(border_shifted=shifted)

        seed_rows, seed_crops = [], []
        for s in range(N_SEEDS):
            rng = np.random.default_rng([s, image_id])
            seed, info = ss.pick_seed(gt_mit, gray, rng, border, rgb.shape)
            sc, sshift = cut_crops(rgb, np.array([seed["cx"]]), np.array([seed["cy"]]),
                                   half_native)
            if sc.shape[1] != CROP_PX:
                import cv2
                sc = np.stack([cv2.resize(c, (CROP_PX, CROP_PX),
                                          interpolation=cv2.INTER_AREA) for c in sc])
            seed_crops.append(sc[0])
            seed_rows.append(dict(seed_index=s, ann_id=int(seed["ann_id"]),
                                  cx=float(seed["cx"]), cy=float(seed["cy"]),
                                  border_shifted=bool(sshift[0]),
                                  pool_after_foreground=int(info.n_after_foreground)))

        blobs.to_csv(f"{OUT}/cands_{fn}.csv", index=False)
        gt.to_csv(f"{OUT}/gt_{fn}.csv", index=False)
        pd.DataFrame(seed_rows).to_csv(f"{OUT}/seeds_{fn}.csv", index=False)
        np.save(f"{OUT}/crops_{fn}.npy", crops)
        np.save(f"{OUT}/seedcrops_{fn}.npy", np.stack(seed_crops))
        json.dump({"file_name": fn, "image_id": image_id, "mpp": mpp, "radius": radius,
                   "roi_shape": list(rgb.shape), "n_candidates": int(len(blobs)),
                   "crop_um": CROP_UM, "crop_px": CROP_PX,
                   "half_native_px": half_native, "n_border_shifted": int(shifted.sum())},
                  open(f"{OUT}/meta_{fn}.json", "w"), indent=1)
        print(f"[{fn}] {len(blobs)} candidates, crops {crops.shape}, "
              f"{shifted.sum()} border-shifted, {time.time() - t0:.0f}s", flush=True)

    bad = [c for c in checks if not c.get("passed", True)]
    print(f"invariant checks: {len(checks)} run, {len(bad)} failed")
    for c in bad:
        print("  FAIL", c)


if __name__ == "__main__":
    sys.exit(main())
