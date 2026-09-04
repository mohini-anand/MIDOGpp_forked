"""Does raising the Otsu foreground threshold un-merge multi-nucleus components?

`Research Logs/2026-08-25-bbox-tuning-domain-mismatch.md` (Experiment 6) found that on
dense-cellular domains, `tighten_box_otsu`'s connected component under the click can span
several adjacent nuclei rather than one -- e.g. `245.tiff` ann 6245: the accepted binary-Otsu
component is the ENTIRE 51x51 window (area_frac=0.58, solidity=0.72), yet still passes every
sanity check (`min_area=50`, `max_area_frac=0.85`, `min_solidity=0.5`). That log's own
conclusion was "a domain property, not fixable by more threshold tuning" -- this script is a
deliberate attempt to test/falsify that conclusion at scale, not an oversight of it.

`design_choices.md` section 7 already tested `method="multiotsu"` (3-class Otsu, keep only the
brightest class) but calibrated against a DIFFERENT motivating example (350.tiff ann 18161, a
"bridges into ONE neighbour" case multiotsu did not fix) and measured aggregate pool-survival /
paired-size-shift, not specifically "does this resolve a component spanning multiple nuclei."
That calibration's own scripts no longer exist in the repo (checked via git history), so its
452-504-candidate numbers are trusted as written, not re-derived here.

Three arms, one shared diagnostic pipeline (`threshold_diag.py`, individually validated against
production `seed_selection.tighten_box_otsu` before this script was written -- frac=0 and
multiotsu both reproduce production byte-for-byte on spot checks including ann 6245 and ann 20):

  (a) binary   -- production default, `cv2.THRESH_OTSU`'s own split.
  (b) multiotsu -- production's existing opt-in, `skimage.filters.threshold_multiotsu` top class.
  (c) headroom sweep -- NEW: threshold = otsu + frac*(255-otsu) for frac in
      {0.05, 0.10, 0.15, 0.20, 0.30}. Not a flat 0-255 pixel delta (would behave inconsistently
      across crops with different histogram shapes -- the same problem the superseded 1D method
      had with a flat relative threshold on busy tissue) -- a fraction of the *headroom* between
      Otsu's own split and the crop's max, which is always 255 since `cv2.normalize` runs before
      thresholding in every arm.

Population: full `agreement_pool` + `border_filter` (tm.PATCH_SIZE//2 = 36px) candidates across
every currently-downloaded ROI (18 files, all 7 domains) -- not a demo per-domain sample, and not
just the "10 ROI" set section 2b/7 used (which can no longer be exactly reconstructed), so
numbers here are NOT directly comparable to design_choices.md's own percentages, only internally
comparable across this script's own arms.

Split (on binary's own diagnostic, per candidate, `otsu_window` pinned to `tm.BASE_SIZE`=51px
everywhere):
  - no_foreground     : click's pixel isn't foreground at all (Experiment-5 territory -- the
                        centre-check/fallback question -- explicitly OUT OF SCOPE here).
  - min_area_rejected : component too small (also out of scope -- unrelated to merging).
  - flagged_rejected  : component fails max_area_frac or min_solidity -- i.e. ALREADY discarded
                        by production as a probable merge. This is the population the earlier,
                        naive version of this analysis (anchored only to *accepted* output) would
                        have missed entirely -- exactly the failure mode a reviewing subagent
                        caught before this script was run.
  - flagged_accepted  : accepted, but area_frac>=0.5 or longer bbox side>=45px -- passes every
                        sanity check yet still looks suspiciously large (ann 6245's bucket).
  - control           : accepted, area_frac<0.3 -- already clean, single-object candidates.
  - gap               : accepted, everything else -- reported as a population share, not
                        otherwise analysed.

Every cross-arm metric is paired (binary vs. one arm at a time, not a many-way intersection --
a 7-way intersection over a suspicious population is close to a guaranteed null, per the same
reviewer).

Writes results/bbox_threshold_sweep.csv (one row per candidate x arm) and prints per-domain and
pooled summary tables. Does NOT measure downstream template-matching/discrimination quality --
same explicit gap `design_choices.md` section 7b already admits it never closed.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from threshold_diag import diagnose_binary, diagnose_multiotsu, diagnose_headroom, OTSU_WINDOW  # noqa: E402

from midog_utils import dataset, channels, seed_selection as ss, template_match as tm  # noqa: E402

FRACS = [0.05, 0.10, 0.15, 0.20, 0.30]
ARMS = ["multiotsu"] + [f"frac_{f:.2f}" for f in FRACS]

MERGE_AREA_FRAC = 0.5
MERGE_LONGER_SIDE = 45
CONTROL_AREA_FRAC = 0.3


def classify(diag_bin: dict) -> str:
    if not diag_bin["center_foreground"]:
        return "no_foreground"
    if not diag_bin["accepted"]:
        if diag_bin["reject_reason"] == "min_area":
            return "min_area_rejected"
        return "flagged_rejected"  # max_area_frac or min_solidity
    if diag_bin["area_frac"] >= MERGE_AREA_FRAC or diag_bin["longer_side"] >= MERGE_LONGER_SIDE:
        return "flagged_accepted"
    if diag_bin["area_frac"] < CONTROL_AREA_FRAC:
        return "control"
    return "gap"


def looks_merged(diag: dict) -> bool:
    return diag["accepted"] and (diag["area_frac"] >= MERGE_AREA_FRAC or diag["longer_side"] >= MERGE_LONGER_SIDE)


def run_arm(name: str, patch: np.ndarray) -> dict:
    if name == "multiotsu":
        return diagnose_multiotsu(patch)
    frac = float(name.split("_")[1])
    return diagnose_headroom(patch, frac)


def main():
    t0 = time.time()
    images, anns = dataset.load_annotations()
    files = sorted(p.name for p in Path("images").glob("*.tiff"))

    rows = []
    for fn in files:
        mit = dataset.image_annotations(anns, fn, category_id=1)
        if len(mit) == 0:
            continue
        meta = images[images.file_name == fn].iloc[0]
        pool, agreement_flagged = ss.agreement_pool(mit)
        pool = ss.border_filter(pool, tm.PATCH_SIZE // 2, (meta.height, meta.width))
        if len(pool) == 0:
            continue

        rgb = dataset.load_roi(f"images/{fn}")
        structural = channels.to_gray_inverted(rgb)

        for _, row in pool.iterrows():
            cx, cy = float(row["cx"]), float(row["cy"])
            patch = tm.read_padded_patch(structural, cx, cy, OTSU_WINDOW)
            if patch is None:
                continue  # shouldn't happen post border_filter, but stay defensive

            diag_bin = diagnose_binary(patch)
            bucket = classify(diag_bin)

            rec = dict(
                file_name=fn, ann_id=int(row["ann_id"]), tumor_type=meta.tumor_type,
                bucket=bucket,
                bin_accepted=diag_bin["accepted"], bin_reject_reason=diag_bin["reject_reason"],
                bin_area=diag_bin["area"], bin_solidity=diag_bin["solidity"],
                bin_area_frac=diag_bin["area_frac"], bin_longer_side=diag_bin["longer_side"],
                bin_threshold=diag_bin["threshold_used"],
            )
            for arm in ARMS:
                d = run_arm(arm, patch)
                rec[f"{arm}_accepted"] = d["accepted"]
                rec[f"{arm}_reject_reason"] = d["reject_reason"]
                rec[f"{arm}_area"] = d["area"]
                rec[f"{arm}_solidity"] = d["solidity"]
                rec[f"{arm}_area_frac"] = d["area_frac"]
                rec[f"{arm}_longer_side"] = d["longer_side"]
                rec[f"{arm}_looks_merged"] = looks_merged(d) if d["accepted"] else None
            rows.append(rec)

        print(f"{fn} done ({len(pool)} candidates), {time.time()-t0:.0f}s elapsed", flush=True)

    df = pd.DataFrame(rows)
    Path("results").mkdir(exist_ok=True)
    df.to_csv("results/bbox_threshold_sweep.csv", index=False)
    print(f"\nwrote results/bbox_threshold_sweep.csv, {len(df)} rows, {time.time()-t0:.0f}s total")

    print("\n=== population, pooled ===")
    print(df["bucket"].value_counts().to_string())

    print("\n=== population, per domain ===")
    print(pd.crosstab(df["tumor_type"], df["bucket"]).to_string())


if __name__ == "__main__":
    main()
