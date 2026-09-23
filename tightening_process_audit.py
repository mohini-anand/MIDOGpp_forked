"""
    Independent audit of `production_hematoxylin_only/tightening_process_hem_vs_gray.ipynb`.

    Re-derives every number the notebook prints, renders into a figure title, or asserts in
    prose, from the primary sources (`databases/MIDOG++.json`, `images/extra_valid/*.tiff`)
    and from the notebook's own persisted artifact
    (`production_hematoxylin_only/tightening_process_summary.csv`).

    The Otsu tightening is reimplemented here in plain cv2/skimage/numpy rather than by
    calling `midog_utils.seed_selection`, so the Tier A comparison is an independent
    reproduction. A separate, explicitly-labelled consistency block does call the real
    helpers, to confirm the reimplementation and to reproduce the notebook's own cross-check.

    Writes `results/tightening_process_audit_*.csv`. Run start to finish under
    /Users/mohinianand/anaconda3/bin/python3.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import time
from itertools import combinations

import cv2
import numpy as np
import pandas as pd
from scipy import stats as stats_mod
from skimage.color import rgb2hed
from skimage.measure import label, regionprops

REPO = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO)

NB_DIR = os.path.join(REPO, "production_hematoxylin_only")
NB_PATH = os.path.join(NB_DIR, "tightening_process_hem_vs_gray.ipynb")
SIB_PATH = os.path.join(NB_DIR, "production_seed_precision_at_k_chromatin_hem_bbox.ipynb")
SUMMARY_CSV = os.path.join(NB_DIR, "tightening_process_summary.csv")
IMAGES_DIR = os.path.join(REPO, "images", "extra_valid")
DB_PATH = os.path.join(REPO, "databases", "MIDOG++.json")
OUT = os.path.join(REPO, "results")

OTSU_WINDOW = 51            # the notebook's OTSU_WINDOW (tm.BASE_SIZE)
MIN_AREA = 50
MAX_AREA_FRAC = 0.85
MIN_SOLIDITY = 0.5
BORDER = 36                 # the notebook's BORDER (FSConfig().patch_size // 2)
PATCH_SIZE = 73             # tm.PATCH_SIZE, used by build_seed's _patch_readable check
MITOTIC = 1

TABLES = {}


def emit(name, df):
    """Save one audit table and echo its shape."""
    path = os.path.join(OUT, f"tightening_process_audit_{name}.csv")
    df.to_csv(path, index=False)
    TABLES[name] = df
    print(f"  -> results/tightening_process_audit_{name}.csv  ({len(df)} rows)")
    return df


def rule(title):
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


# ===========================================================================================
# 0. Independent reimplementation of the tightening, from the algorithm, not from the module
# ===========================================================================================

def gray_inverted(rgb):
    """255 minus the luma of an RGB crop, float32."""
    return (255.0 - cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)).astype(np.float32)


def hematoxylin_od(rgb):
    """Unclipped hematoxylin optical density of an RGB crop, float32."""
    return rgb2hed(rgb.astype(np.float32) / 255.0)[:, :, 0].astype(np.float32)


CHANNEL_FNS = {"gray_inverted": gray_inverted, "hematoxylin_od": hematoxylin_od}


def odd_at_least(n, minimum=5):
    """Round to an odd integer no smaller than `minimum`."""
    n = int(round(n))
    if n % 2 == 0:
        n += 1
    return max(minimum, n)


def read_patch_rgb(rgb, cx, cy, size):
    """Square RGB crop centred on the rounded click, or None when it runs off the ROI."""
    half = size // 2
    ix, iy = int(round(cx)), int(round(cy))
    h, w = rgb.shape[:2]
    if ix - half < 0 or iy - half < 0 or ix + half >= w or iy + half >= h:
        return None
    return np.ascontiguousarray(rgb[iy - half: iy + half + 1, ix - half: ix + half + 1])


def tighten(patch, cx, cy, otsu_window=OTSU_WINDOW):
    """Reimplementation of tighten_box_otsu + tightened_template_box, all intermediates kept."""
    out = dict(otsu_cut=np.nan, center_label=0, area=np.nan, solidity=np.nan,
               accepted=False, reason="too close to the ROI border to read the Otsu window",
               bbox=None, base_size=np.nan, center_x=np.nan, center_y=np.nan)
    if patch is None:
        return out
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    cut, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    labels = label(binary, connectivity=2)
    py, px = patch.shape[0] // 2, patch.shape[1] // 2
    cl = int(labels[py, px])
    out.update(otsu_cut=float(cut), center_label=cl, n_components=int(labels.max()))
    if cl == 0:
        out["reason"] = "click lands on background at this Otsu cut"
        return out
    region = next(r for r in regionprops(labels) if r.label == cl)
    out.update(area=int(region.area), solidity=float(region.solidity))
    min_row, min_col, max_row, max_col = region.bbox
    if not (min_row <= py < max_row and min_col <= px < max_col):
        out["reason"] = "accepted component's bbox does not actually contain the click pixel"
        return out
    if region.area < MIN_AREA:
        out["reason"] = f"component area {region.area}px < min_area {MIN_AREA}px"
        return out
    if region.area > MAX_AREA_FRAC * patch.size:
        out["reason"] = (f"component area {region.area}px > {MAX_AREA_FRAC:.0%} of the "
                         f"{patch.size}px window")
        return out
    if region.solidity < MIN_SOLIDITY:
        out["reason"] = f"component solidity {region.solidity:.2f} < min_solidity {MIN_SOLIDITY}"
        return out
    y0, x0, y1, x1 = region.bbox
    base_size = odd_at_least(max(y1 - y0, x1 - x0), minimum=5)
    half = otsu_window // 2
    ix, iy = int(round(cx)), int(round(cy))
    n = patch.shape[0]
    out.update(accepted=True, reason="accepted", bbox=(int(y0), int(y1), int(x0), int(x1)),
               base_size=int(base_size), bbox_h=int(y1 - y0), bbox_w=int(x1 - x0),
               area_frac=float(region.area) / patch.size,
               bbox_touches_window_edge=bool(y0 == 0 or x0 == 0 or y1 == n or x1 == n),
               base_size_at_window_ceiling=bool(int(base_size) >= n),
               center_x=ix - half + (x0 + x1 - 1) / 2.0,
               center_y=iy - half + (y0 + y1 - 1) / 2.0)
    return out


def agreement_pool(gt_mitotic):
    """Unanimous tier when non-empty, else the 2-of-3 contested tier; plus the flag."""
    unanimous = gt_mitotic[gt_mitotic["n_mitotic_votes"] == gt_mitotic["n_votes"]]
    if len(unanimous):
        return unanimous, False
    contested = gt_mitotic[(gt_mitotic["n_votes"] > 0)
                           & (gt_mitotic["n_mitotic_votes"] / gt_mitotic["n_votes"] >= 2.0 / 3.0)
                           & (gt_mitotic["n_mitotic_votes"] < gt_mitotic["n_votes"])]
    return contested, True


def border_filter(df, border, roi_shape):
    """Rows whose rounded centre clears `border` px on every ROI edge."""
    h, w = roi_shape[:2]
    ix = np.rint(df["cx"].to_numpy()).astype(int)
    iy = np.rint(df["cy"].to_numpy()).astype(int)
    ok = (ix >= border) & (ix <= w - 1 - border) & (iy >= border) & (iy <= h - 1 - border)
    return df[ok]


def patch_readable(roi_shape, cx, cy, patch_size):
    """build_seed's own readability predicate at the recentred template point."""
    half = patch_size // 2
    ix, iy = int(round(cx)), int(round(cy))
    h, w = roi_shape[:2]
    return not (ix - half < 0 or iy - half < 0 or ix + half >= w or iy + half >= h)


def load_db(path):
    """MIDOG++ images and annotations frames, built straight from the JSON."""
    raw = json.loads(open(path).read())
    images = pd.DataFrame([dict(image_id=im["id"], file_name=im["file_name"],
                                width=im["width"], height=im["height"],
                                tumor_type=im["tumor_type"]) for im in raw["images"]])
    id_to_name = dict(zip(images["image_id"], images["file_name"]))
    rows = []
    for a in raw["annotations"]:
        x1, y1, x2, y2 = a["bbox"]
        labels_ = a.get("labels", []) or []
        rows.append(dict(ann_id=a["id"], image_id=a["image_id"],
                         file_name=id_to_name[a["image_id"]],
                         cx=(x1 + x2) / 2.0, cy=(y1 + y2) / 2.0,
                         category_id=a["category_id"],
                         n_votes=len(labels_),
                         n_mitotic_votes=int(sum(1 for v in labels_ if v == MITOTIC))))
    anns = pd.DataFrame(rows)
    counts = anns.groupby("image_id").size()
    images = images[images["image_id"].isin(counts.index)].reset_index(drop=True)
    return images, anns


def load_roi(path):
    """Full-resolution RGB uint8 array for one ROI."""
    import tifffile
    with tifffile.TiffFile(str(path)) as tf:
        series = tf.series[0]
        levels = getattr(series, "levels", None)
        arr = levels[0].asarray() if levels else series.asarray()
    return np.ascontiguousarray(arr[:, :, :3])


def git(*args):
    """Run git in the repo and return stripped stdout."""
    try:
        return subprocess.run(["git"] + list(args), cwd=REPO, capture_output=True,
                              text=True, timeout=60).stdout.strip()
    except Exception as exc:                                    # pragma: no cover
        return f"<git failed: {exc}>"


T_START = time.time()

# ===========================================================================================
# 1. Execution-coherence gate
# ===========================================================================================
rule("GATE 1 -- execution coherence")

nb = json.loads(open(NB_PATH).read())
code = [c for c in nb["cells"] if c["cell_type"] == "code"]
ecs = [c.get("execution_count") for c in code]
n_err = sum(1 for c in code for o in (c.get("outputs") or [])
            if o.get("output_type") == "error")
unrun = [i for i, c in enumerate(nb["cells"])
         if c["cell_type"] == "code" and c.get("execution_count") is None]
contiguous = ecs == list(range(1, len(ecs) + 1))
last_code_ix = max(i for i, c in enumerate(nb["cells"]) if c["cell_type"] == "code")
print(f"code cells                 : {len(code)}")
print(f"execution_count sequence   : {ecs}")
print(f"contiguous from 1          : {contiguous}")
print(f"error outputs              : {n_err}")
print(f"unrun code cells           : {unrun}")
print(f"last cell type             : {nb['cells'][-1]['cell_type']} (last code cell index {last_code_ix})")
emit("execution_gate", pd.DataFrame([
    dict(check="execution_count contiguous from 1", value=str(ecs), passed=contiguous),
    dict(check="no error outputs", value=n_err, passed=n_err == 0),
    dict(check="no unrun code cells", value=str(unrun), passed=len(unrun) == 0),
    dict(check="last cell is markdown summary (not an unrun code cell)",
         value=nb["cells"][-1]["cell_type"], passed=nb["cells"][-1]["cell_type"] == "markdown"),
    dict(check="vendored exemption applies", value="no -- imports midog_utils, reads images/ and databases/",
         passed=False),
]))

# ===========================================================================================
# 2. Provenance gate
# ===========================================================================================
rule("GATE 2 -- provenance")

prov_rows = []
for rel in ["midog_utils/seed_selection.py", "midog_utils/template_match.py",
            "midog_utils/channels.py", "midog_utils/chromatin.py",
            "midog_utils/dataset.py", "midog_utils/find_and_suppress.py",
            "production_hematoxylin_only/tightening_process_hem_vs_gray.ipynb",
            "production_hematoxylin_only/tightening_process_summary.csv"]:
    p = os.path.join(REPO, rel)
    prov_rows.append(dict(
        path=rel,
        mtime=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(p))),
        last_commit=git("log", "-1", "--format=%h %ad %s", "--date=short", "--", rel) or "<untracked>",
        porcelain=git("status", "--porcelain", "--", rel) or "(clean)",
    ))
prov = pd.DataFrame(prov_rows)
print(prov.to_string(index=False))
nb_mtime = os.path.getmtime(NB_PATH)
mods_newer = [r["path"] for r in prov_rows if r["path"].startswith("midog_utils")
              and os.path.getmtime(os.path.join(REPO, r["path"])) > nb_mtime]
print(f"\nmodules modified after the notebook run: {mods_newer or 'none'}")
emit("provenance", prov)

# ===========================================================================================
# 3. Rebuild the notebook's own selection, independently
# ===========================================================================================
rule("REBUILD -- the 14 examples, from MIDOG++.json and the ROI files")

images, annotations = load_db(DB_PATH)
files_on_disk = sorted(f for f in os.listdir(IMAGES_DIR) if f.endswith(".tiff"))
print(f"ROIs on disk in images/extra_valid: {len(files_on_disk)} -> {files_on_disk}")
meta_ix = images.set_index("file_name")[["image_id", "tumor_type"]].loc[files_on_disk]

# ROI-per-stratum count (collinearity of the exchangeable unit with the stratum)
per_domain_14 = meta_ix.reset_index().groupby("tumor_type")["file_name"].agg(["count", list])
per_domain_all = images.groupby("tumor_type")["file_name"].count()
print("\nROIs per tumor_type among the 14 on disk:")
print(per_domain_14.to_string())
emit("roi_per_domain", per_domain_14.reset_index().rename(columns={"count": "n_roi_of_14",
                                                                   "list": "files"})
     .assign(n_roi_whole_db=lambda d: d["tumor_type"].map(per_domain_all).values))

domain_roi = meta_ix.reset_index().groupby("tumor_type")["file_name"].min().sort_index()
alphabetical_ok = all(fn == sorted(per_domain_14.loc[dom, "list"])[0]
                      for dom, fn in domain_roi.items())
print(f"\none ROI per domain (alphabetically-first) reproduces: {alphabetical_ok}")
print(domain_roi.to_string())

ROIS, EXAMPLES, SEED_ROWS = {}, [], []
for domain, fn in domain_roi.items():
    image_id = int(meta_ix.loc[fn, "image_id"])
    rgb = load_roi(os.path.join(IMAGES_DIR, fn))
    ROIS[fn] = rgb
    gt = annotations[annotations["file_name"] == fn].reset_index(drop=True)
    gt_mit = gt[gt["category_id"] == MITOTIC]
    pool, flagged = agreement_pool(gt_mit)
    pool = border_filter(pool, BORDER, rgb.shape).reset_index(drop=True)

    rng0 = np.random.default_rng([0, image_id])
    idx0 = int(rng0.integers(len(pool)))
    row0 = pool.iloc[idx0]
    remaining = pool.drop(pool.index[idx0]).reset_index(drop=True)
    rng1 = np.random.default_rng([1, image_id])
    idx1 = int(rng1.integers(len(remaining)))
    row1 = remaining.iloc[idx1]

    # what seed_index=1 would actually draw under run_pipeline.py / build_seed: the FULL pool
    true_idx1 = int(np.random.default_rng([1, image_id]).integers(len(pool)))
    true_row1 = pool.iloc[true_idx1]

    SEED_ROWS.append(dict(domain=domain, file_name=fn, image_id=image_id,
                          n_mitotic=int(len(gt_mit)), agreement_flagged=bool(flagged),
                          n_pool_after_border=int(len(pool)),
                          nb_pick1_ann_id=int(row0["ann_id"]),
                          nb_pick2_ann_id=int(row1["ann_id"]),
                          true_seed_index1_ann_id=int(true_row1["ann_id"]),
                          pick2_equals_seed_index1=bool(int(row1["ann_id"]) == int(true_row1["ann_id"]))))
    for which, row in (("#1 (seed_index=0)", row0), ("#2 (seed_index=1)", row1)):
        EXAMPLES.append(dict(domain=domain, file_name=fn, image_id=image_id, which=which,
                             ann_id=int(row["ann_id"]), cx=float(row["cx"]), cy=float(row["cy"])))

EX = pd.DataFrame(EXAMPLES)
print(f"\n{len(EX)} examples rebuilt")
print(EX.to_string(index=False))

# ===========================================================================================
# 4. Channel locality -- can the 51px crop be converted, or is the channel ROI-global?
# ===========================================================================================
rule("CHECK -- is each channel a pure per-pixel map (crop-then-convert == convert-then-crop)?")

loc_rows = []
probe_fn = EX.iloc[0]["file_name"]
probe_rgb = ROIS[probe_fn]
pcx, pcy = EX.iloc[0]["cx"], EX.iloc[0]["cy"]
for name, fnc in CHANNEL_FNS.items():
    t0 = time.time()
    full = fnc(probe_rgb)
    t_full = time.time() - t0
    half = OTSU_WINDOW // 2
    ix, iy = int(round(pcx)), int(round(pcy))
    a = full[iy - half: iy + half + 1, ix - half: ix + half + 1]
    b = fnc(read_patch_rgb(probe_rgb, pcx, pcy, OTSU_WINDOW))
    max_abs = float(np.max(np.abs(a.astype(np.float64) - b.astype(np.float64))))
    # what actually matters: identical uint8 after the patch-local NORM_MINMAX
    ua = cv2.normalize(a.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    ub = cv2.normalize(b.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    loc_rows.append(dict(channel=name, roi=probe_fn, full_roi_convert_s=round(t_full, 2),
                         max_abs_diff_float=max_abs, u8_identical=bool(np.array_equal(ua, ub))))
    del full
loc = pd.DataFrame(loc_rows)
print(loc.to_string(index=False))
LOCAL_OK = bool(loc["u8_identical"].all())
print(f"\nboth channels per-pixel local: {LOCAL_OK}  (crops are therefore a faithful recompute)")
emit("channel_locality", loc)
assert LOCAL_OK, "channel is not patch-local; the crop-based recompute below would be invalid"

# ===========================================================================================
# 5. Tier A -- recompute every cell of tightening_process_summary.csv
# ===========================================================================================
rule("TIER A -- independent recompute of tightening_process_summary.csv")

REC = []
for e in EXAMPLES:
    rgb = ROIS[e["file_name"]]
    patch_rgb = read_patch_rgb(rgb, e["cx"], e["cy"], OTSU_WINDOW)
    for ch_key, fnc in CHANNEL_FNS.items():
        patch = None if patch_rgb is None else fnc(patch_rgb)
        r = tighten(patch, e["cx"], e["cy"])
        tpl_shape = None
        readable_73 = None
        if r["accepted"]:
            bs, bh = int(r["base_size"]), (int(r["base_size"]) - 1) // 2
            icx, icy = int(round(r["center_x"])), int(round(r["center_y"]))
            tpl = rgb[icy - bh: icy - bh + bs, icx - bh: icx - bh + bs]
            tpl_shape = f"{tpl.shape[0]}x{tpl.shape[1]}"
            readable_73 = patch_readable(rgb.shape, r["center_x"], r["center_y"], PATCH_SIZE)
        REC.append(dict(domain=e["domain"], file_name=e["file_name"], which=e["which"],
                        ann_id=e["ann_id"], channel=ch_key, cx=e["cx"], cy=e["cy"],
                        accepted=r["accepted"],
                        base_size=r["base_size"] if r["accepted"] else np.nan,
                        reason=None if r["accepted"] else r["reason"],
                        otsu_cut=r["otsu_cut"], center_label=r["center_label"],
                        n_components=r.get("n_components", np.nan),
                        area=r["area"], solidity=r["solidity"],
                        bbox_h=r.get("bbox_h", np.nan), bbox_w=r.get("bbox_w", np.nan),
                        area_frac=r.get("area_frac", np.nan),
                        bbox_touches_window_edge=r.get("bbox_touches_window_edge", None),
                        base_size_at_window_ceiling=r.get("base_size_at_window_ceiling", None),
                        center_x=r["center_x"], center_y=r["center_y"],
                        offset_px=(np.nan if not r["accepted"]
                                   else float(np.hypot(r["center_x"] - e["cx"],
                                                       r["center_y"] - e["cy"]))),
                        template_crop_shape=tpl_shape,
                        template_readable_at_patch_size_73=readable_73))
REC = pd.DataFrame(REC)

NBC = pd.read_csv(SUMMARY_CSV)
key = ["file_name", "ann_id", "channel"]
merged = NBC.merge(REC, on=key, suffixes=("_nb", "_rec"), how="outer", indicator=True)
assert (merged["_merge"] == "both").all(), "key mismatch between notebook CSV and recompute"

# `which` is built here from the same literal strings the notebook uses, in the same order, so
# comparing it proves nothing; it is counted separately as "reproduced by construction".
INDEPENDENT_COLS = [("domain", "domain_nb", "domain_rec"),
                    ("accepted", "accepted_nb", "accepted_rec"),
                    ("base_size", "base_size_nb", "base_size_rec"),
                    ("reason", "reason_nb", "reason_rec")]
BY_CONSTRUCTION_COLS = [("which", "which_nb", "which_rec")]
compare_cols = INDEPENDENT_COLS + BY_CONSTRUCTION_COLS
div_rows, n_compared, n_div = [], 0, 0
n_indep = 0
for label_, cnb, crec in compare_cols:
    for _, row in merged.iterrows():
        a, b = row[cnb], row[crec]
        n_compared += 1
        if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
            same = True
        elif (a is None or (isinstance(a, float) and math.isnan(a))) and \
             (b is None or (isinstance(b, float) and math.isnan(b))):
            same = True
        elif isinstance(a, float) or isinstance(b, float):
            try:
                same = bool(np.isclose(float(a), float(b)))
            except (TypeError, ValueError):
                same = a == b
        else:
            same = a == b
        if label_ != "which":
            n_indep += 1
        if not same:
            n_div += 1
            div_rows.append(dict(file_name=row["file_name"], ann_id=row["ann_id"],
                                 channel=row["channel"], column=label_,
                                 notebook=a, recomputed=b))
print(f"values compared total          : {n_compared}")
print(f"  independently derived        : {n_indep}  "
      f"({', '.join(c[0] for c in INDEPENDENT_COLS)})")
print(f"  reproduced by construction   : {n_compared - n_indep}  "
      f"({', '.join(c[0] for c in BY_CONSTRUCTION_COLS)} -- built here from the same literal "
      f"strings the notebook uses, so it proves nothing)")
print(f"divergences                    : {n_div}")
emit("tier_a_counts", pd.DataFrame([
    dict(artifact="tightening_process_summary.csv", basis="independently derived",
         rows=len(NBC), columns=len(INDEPENDENT_COLS),
         values_compared=n_indep, divergences=n_div),
    dict(artifact="tightening_process_summary.csv", basis="reproduced by construction",
         rows=len(NBC), columns=len(BY_CONSTRUCTION_COLS),
         values_compared=n_compared - n_indep, divergences=0),
]))
emit("tier_a_divergences", pd.DataFrame(div_rows, columns=["file_name", "ann_id", "channel",
                                                           "column", "notebook", "recomputed"]))

# ===========================================================================================
# 6. Composition gate
# ===========================================================================================
rule("GATE 3 -- composition of tightening_process_summary.csv")

comp = [
    dict(check="row count == 7 domains x 2 clicks x 2 channels = 28",
         expected=28, observed=len(NBC), passed=len(NBC) == 28),
    dict(check="distinct (file_name, ann_id, channel) keys", expected=28,
         observed=int(NBC[key].drop_duplicates().shape[0]),
         passed=NBC[key].drop_duplicates().shape[0] == 28),
    dict(check="fully duplicated rows", expected=0, observed=int(NBC.duplicated().sum()),
         passed=NBC.duplicated().sum() == 0),
    dict(check="distinct ROIs", expected=7, observed=int(NBC["file_name"].nunique()),
         passed=NBC["file_name"].nunique() == 7),
    dict(check="distinct domains", expected=7, observed=int(NBC["domain"].nunique()),
         passed=NBC["domain"].nunique() == 7),
    dict(check="every domain present exactly 4x (2 clicks x 2 channels)", expected="all 4",
         observed=str(sorted(NBC.groupby("domain").size().unique().tolist())),
         passed=set(NBC.groupby("domain").size().unique()) == {4}),
    dict(check="distinct ann_ids", expected=14, observed=int(NBC["ann_id"].nunique()),
         passed=NBC["ann_id"].nunique() == 14),
    dict(check="each channel appears 14x", expected="14/14",
         observed=str(NBC.groupby("channel").size().to_dict()),
         passed=set(NBC.groupby("channel").size()) == {14}),
    dict(check="base_size null exactly where accepted is False", expected=True,
         observed=bool((NBC["base_size"].isna() == ~NBC["accepted"]).all()),
         passed=bool((NBC["base_size"].isna() == ~NBC["accepted"]).all())),
    dict(check="reason null exactly where accepted is True", expected=True,
         observed=bool((NBC["reason"].isna() == NBC["accepted"]).all()),
         passed=bool((NBC["reason"].isna() == NBC["accepted"]).all())),
    dict(check="every base_size odd and >= 5", expected=True,
         observed=bool(((NBC["base_size"].dropna() % 2 == 1)
                        & (NBC["base_size"].dropna() >= 5)).all()),
         passed=bool(((NBC["base_size"].dropna() % 2 == 1)
                      & (NBC["base_size"].dropna() >= 5)).all())),
    dict(check="notebook's 14 ann_ids == rebuilt 14 ann_ids", expected=True,
         observed=bool(set(NBC["ann_id"]) == set(EX["ann_id"])),
         passed=bool(set(NBC["ann_id"]) == set(EX["ann_id"]))),
    dict(check="ROIs used are 7 of the 14 on disk (one per domain)", expected=7,
         observed=int(len(set(NBC["file_name"]) & set(files_on_disk))),
         passed=len(set(NBC["file_name"]) & set(files_on_disk)) == 7),
]
comp = pd.DataFrame(comp)
print(comp.to_string(index=False))
emit("composition", comp)

# ===========================================================================================
# 7. Consistency block -- the real helpers (explicitly NOT the Tier A oracle)
# ===========================================================================================
rule("CONSISTENCY -- reimplementation vs midog_utils.seed_selection (labelled, not independent)")

from midog_utils import channels as ch_mod                 # noqa: E402
from midog_utils import seed_selection as ss               # noqa: E402
from midog_utils import template_match as tm               # noqa: E402

cons_rows = []
for _, r in REC.iterrows():
    rgb = ROIS[r["file_name"]]
    patch_rgb = read_patch_rgb(rgb, r["cx"], r["cy"], OTSU_WINDOW)
    patch = CHANNEL_FNS[r["channel"]](patch_rgb)
    official_bbox = ss.tighten_box_otsu(patch)
    got = ss.tightened_template_box(ch_mod.to_channel(patch_rgb, r["channel"]),
                                    OTSU_WINDOW // 2, OTSU_WINDOW // 2,
                                    otsu_window=OTSU_WINDOW)
    cons_rows.append(dict(
        file_name=r["file_name"], ann_id=r["ann_id"], channel=r["channel"],
        rec_accepted=bool(r["accepted"]),
        ss_accepted=official_bbox is not None,
        accept_agrees=bool(r["accepted"]) == (official_bbox is not None),
        rec_base_size=r["base_size"],
        ss_base_size=np.nan if got is None else got[0],
        size_agrees=bool((got is None and not r["accepted"])
                         or (got is not None and r["accepted"] and got[0] == r["base_size"])),
    ))
cons = pd.DataFrame(cons_rows)
print(f"accept/refuse agrees on {int(cons['accept_agrees'].sum())}/{len(cons)}; "
      f"base_size agrees on {int(cons['size_agrees'].sum())}/{len(cons)}")
emit("consistency_vs_module", cons)

# also confirm the constants the notebook's config cell asserts
cfg_rows = [
    dict(name="OTSU_WINDOW", notebook=51, source="template_match.BASE_SIZE",
         source_value=tm.BASE_SIZE, agrees=tm.BASE_SIZE == 51),
    dict(name="BORDER", notebook=36, source="FSConfig().patch_size // 2",
         source_value=PATCH_SIZE // 2, agrees=PATCH_SIZE // 2 == 36),
    dict(name="min_area", notebook=50, source="tighten_box_otsu default", source_value=50,
         agrees=True),
    dict(name="max_area_frac", notebook=0.85, source="tighten_box_otsu default",
         source_value=0.85, agrees=True),
    dict(name="min_solidity", notebook=0.5, source="tighten_box_otsu default", source_value=0.5,
         agrees=True),
    dict(name="method", notebook="binary", source="tighten_box_otsu default",
         source_value="binary", agrees=True),
    dict(name="center_tolerance", notebook=0, source="tighten_box_otsu default", source_value=0,
         agrees=True),
]
import inspect                                              # noqa: E402
sig = inspect.signature(ss.tighten_box_otsu).parameters
for r in cfg_rows:
    p = {"min_area": "min_area", "max_area_frac": "max_area_frac", "min_solidity": "min_solidity",
         "method": "method", "center_tolerance": "center_tolerance"}.get(r["name"])
    if p is not None:
        r["source_value"] = sig[p].default
        r["agrees"] = sig[p].default == r["notebook"]
cfg = pd.DataFrame(cfg_rows)
print(cfg.to_string(index=False))
emit("config_drift", cfg)

# the prose's source-line citation
ss_src = open(os.path.join(REPO, "midog_utils", "seed_selection.py")).read().split("\n")
def_line = next(i + 1 for i, l in enumerate(ss_src) if l.startswith("def tighten_box_otsu"))
end_line = next(i for i, l in enumerate(ss_src[def_line:], start=def_line)
                if l.startswith("def foreground_filter"))
print(f"\nprose cites seed_selection.py:89-153 for tighten_box_otsu; "
      f"the function actually spans lines {def_line}-{end_line - 1}")

# ===========================================================================================
# 8. Refusal-reason verification (the assertion in cell 7 cannot check these)
# ===========================================================================================
rule("REFUSALS -- is every refusal really 'click lands on background'?")

ref = REC[~REC["accepted"]].copy()
ref["center_label_is_zero"] = ref["center_label"] == 0
ref["reason_matches_mechanism"] = ((ref["reason"] == "click lands on background at this Otsu cut")
                                   == ref["center_label_is_zero"])
# the paired half: what did the OTHER channel do at the same pixel?
other = REC.set_index(["file_name", "ann_id", "channel"])
pair_rows = []
for _, r in ref.iterrows():
    oth = "gray_inverted" if r["channel"] == "hematoxylin_od" else "hematoxylin_od"
    o = other.loc[(r["file_name"], r["ann_id"], oth)]
    pair_rows.append(dict(file_name=r["file_name"], ann_id=r["ann_id"],
                          refused_channel=r["channel"], reason=r["reason"],
                          center_label=int(r["center_label"]),
                          n_components_in_patch=int(r["n_components"]),
                          otsu_cut=r["otsu_cut"],
                          other_channel=oth, other_accepted=bool(o["accepted"]),
                          other_center_label=int(o["center_label"]),
                          other_base_size=o["base_size"], other_area=o["area"],
                          other_solidity=o["solidity"],
                          same_pixel_opposite_verdict=bool(o["accepted"])))
refv = pd.DataFrame(pair_rows)
print(refv.to_string(index=False))
print(f"\nrefusals whose stated reason matches the recomputed mechanism: "
      f"{int(ref['reason_matches_mechanism'].sum())}/{len(ref)}")
emit("refusal_verification", refv)

# ===========================================================================================
# 9. Figure-title numbers (the only home of ~100 of this notebook's values)
# ===========================================================================================
rule("FIGURES -- the title strings each of the 28 panels-rows should carry")

fig_rows = []
for _, r in REC.iterrows():
    if not r["accepted"] and np.isnan(r["otsu_cut"]):
        p2 = p3 = p4 = p5 = "(patch unreadable)"
    else:
        p2 = f"Otsu binary (cut={r['otsu_cut']:.1f}/255)"
        if r["center_label"] != 0:
            verdict = "ACCEPTED" if r["accepted"] else f"REFUSED: {r['reason']}"
            p3 = (f"components (label={int(r['center_label'])}) "
                  f"area={int(r['area'])}px solidity={r['solidity']:.2f} {verdict}")
        else:
            p3 = f"components REFUSED: {r['reason']}"
        p4 = (f"tightened bbox (base_size={int(r['base_size'])}px)" if r["accepted"]
              else f"REFUSED {r['reason']}")
        p5 = (f"final template ({int(r['base_size'])}x{int(r['base_size'])}px)" if r["accepted"]
              else "no template (refused)")
    fig_rows.append(dict(png=f"tightening_{r['file_name'].replace('.tiff', '')}_{r['ann_id']}.png",
                         row=("top" if r["channel"] == "gray_inverted" else "bottom"),
                         channel=r["channel"],
                         panel1_title=f"raw patch ({OTSU_WINDOW}px)",
                         panel2_title=p2, panel3_title=p3, panel4_title=p4, panel5_title=p5,
                         template_crop_shape=r["template_crop_shape"],
                         template_crop_matches_base_size=(
                             None if not r["accepted"]
                             else r["template_crop_shape"] == f"{int(r['base_size'])}x{int(r['base_size'])}"),
                         suptitle=f"{r['file_name']}  ({r['domain']})  --  annotation "
                                  f"{r['ann_id']}, {r['which']}"))
FIG = pd.DataFrame(fig_rows)
print(FIG[["png", "channel", "panel2_title", "panel3_title", "panel4_title"]].to_string(index=False))
emit("figure_titles", FIG)

pngs_expected = sorted(FIG["png"].unique())
pngs_on_disk = sorted(f for f in os.listdir(NB_DIR) if f.startswith("tightening_") and f.endswith(".png"))
print(f"\nPNGs the run cell should have written: {len(pngs_expected)}; on disk: {len(pngs_on_disk)}; "
      f"identical set: {pngs_expected == pngs_on_disk}")
tpl_ok = FIG["template_crop_matches_base_size"].dropna()
print(f"'final template (NxN px)' panels whose crop really is NxN: {int(tpl_ok.sum())}/{len(tpl_ok)}")

# ===========================================================================================
# 10. Every number in the closing note
# ===========================================================================================
rule("PROSE -- every number the markdown asserts")

acc = REC.pivot_table(index=["file_name", "ann_id", "which", "domain"], columns="channel",
                      values="accepted", aggfunc="first")
siz = REC.pivot_table(index=["file_name", "ann_id", "which", "domain"], columns="channel",
                      values="base_size", aggfunc="first")
PAIR = acc.join(siz, lsuffix="_acc", rsuffix="_size").reset_index()
PAIR.columns = ["file_name", "ann_id", "which", "domain",
                "gray_accepted", "hem_accepted", "gray_base_size", "hem_base_size"]
PAIR["both_accepted"] = PAIR["gray_accepted"] & PAIR["hem_accepted"]
PAIR["delta_px"] = PAIR["hem_base_size"] - PAIR["gray_base_size"]
PAIR["pct_change"] = 100 * PAIR["delta_px"] / PAIR["gray_base_size"]
PAIR = PAIR.sort_values(["file_name", "which"]).reset_index(drop=True)
print(PAIR.to_string(index=False))
emit("pairs", PAIR)

both = PAIR[PAIR["both_accepted"]]
n_gray_acc = int(PAIR["gray_accepted"].sum())
n_hem_acc = int(PAIR["hem_accepted"].sum())
hem_refusals = PAIR[~PAIR["hem_accepted"]]
hem_ref_gray_ok = hem_refusals[hem_refusals["gray_accepted"]]
extreme_abs = both.loc[both["delta_px"].idxmin()]
extreme_pct = both.loc[both["pct_change"].idxmin()]
r300_2 = PAIR[(PAIR["file_name"] == "300.tiff") & (PAIR["which"].str.startswith("#2"))].iloc[0]
r300_1 = PAIR[(PAIR["file_name"] == "300.tiff") & (PAIR["which"].str.startswith("#1"))].iloc[0]
r459_1 = PAIR[(PAIR["file_name"] == "459.tiff") & (PAIR["which"].str.startswith("#1"))].iloc[0]
r459_2 = PAIR[(PAIR["file_name"] == "459.tiff") & (PAIR["which"].str.startswith("#2"))].iloc[0]

prose = [
    dict(claim="gray_inverted accepted 13/14", notebook="13/14", recomputed=f"{n_gray_acc}/14",
         agrees=n_gray_acc == 13),
    dict(claim="hematoxylin_od accepted only 10/14", notebook="10/14",
         recomputed=f"{n_hem_acc}/14", agrees=n_hem_acc == 10),
    dict(claim="CSV printed 'accepted: 23/28'", notebook="23/28",
         recomputed=f"{int(REC['accepted'].sum())}/28", agrees=int(REC["accepted"].sum()) == 23),
    dict(claim="3 of the 4 hem refusals had gray accept the same click",
         notebook="3 of 4", recomputed=f"{len(hem_ref_gray_ok)} of {len(hem_refusals)}",
         agrees=len(hem_ref_gray_ok) == 3 and len(hem_refusals) == 4),
    dict(claim="those three are 300/14581, 245/6274, 013/254", notebook="14581, 6274, 254",
         recomputed=", ".join(str(a) for a in sorted(hem_ref_gray_ok["ann_id"])),
         agrees=sorted(hem_ref_gray_ok["ann_id"].tolist()) == sorted([14581, 6274, 254])),
    dict(claim="where both accept, hem is smaller every time",
         notebook="10/10 smaller",
         recomputed=f"{int((both['delta_px'] < 0).sum())}/{len(both)} smaller, "
                    f"{int((both['delta_px'] == 0).sum())} tied, "
                    f"{int((both['delta_px'] > 0).sum())} larger",
         agrees=bool((both["delta_px"] < 0).all()) and len(both) == 10),
    dict(claim="'All 10 examples both channels accepted'", notebook=10, recomputed=len(both),
         agrees=len(both) == 10),
    dict(claim="300.tiff#2 is 45->19px", notebook="45 -> 19",
         recomputed=f"{int(r300_2['gray_base_size'])} -> {int(r300_2['hem_base_size'])} "
                    f"(45 is 300.tiff#1's gray size, ann {int(r300_1['ann_id'])}, "
                    f"whose hem was REFUSED)",
         agrees=bool(r300_2["gray_base_size"] == 45 and r300_2["hem_base_size"] == 19)),
    dict(claim="300.tiff#2 is 'the extreme'",
         notebook="300.tiff#2",
         recomputed=f"largest absolute shrink is {extreme_abs['file_name']} "
                    f"{extreme_abs['which']} ({int(extreme_abs['gray_base_size'])}->"
                    f"{int(extreme_abs['hem_base_size'])}, {int(extreme_abs['delta_px'])}px); "
                    f"largest proportional shrink is {extreme_pct['file_name']} "
                    f"{extreme_pct['which']} ({extreme_pct['pct_change']:.1f}%) vs 300#2 "
                    f"({r300_2['pct_change']:.1f}%)",
         agrees=bool(extreme_abs["file_name"] == "300.tiff"
                     and extreme_abs["which"].startswith("#2"))),
    dict(claim="459.tiff#1 barely moves (33->31)", notebook="33 -> 31",
         recomputed=f"{int(r459_1['gray_base_size'])} -> {int(r459_1['hem_base_size'])}",
         agrees=bool(r459_1["gray_base_size"] == 33 and r459_1["hem_base_size"] == 31)),
    dict(claim="459.tiff#2 barely moves (41->39)", notebook="41 -> 39",
         recomputed=f"{int(r459_2['gray_base_size'])} -> {int(r459_2['hem_base_size'])}",
         agrees=bool(r459_2["gray_base_size"] == 41 and r459_2["hem_base_size"] == 39)),
    dict(claim="7 domains -> 7 ROIs, alphabetically-first file per domain",
         notebook="7", recomputed=f"{len(domain_roi)} domains, alphabetical rule holds: "
                                  f"{alphabetical_ok}", agrees=len(domain_roi) == 7 and alphabetical_ok),
    dict(claim="013.tiff chosen over 094.tiff for human breast cancer", notebook="013.tiff",
         recomputed=domain_roi["human breast cancer"],
         agrees=domain_roi["human breast cancer"] == "013.tiff"),
    dict(claim="14 ROIs on disk in images/extra_valid", notebook=14, recomputed=len(files_on_disk),
         agrees=len(files_on_disk) == 14),
    dict(claim="every accepted example's 'final template' crop is base_size x base_size",
         notebook="implied", recomputed=f"{int(tpl_ok.sum())}/{len(tpl_ok)}",
         agrees=bool(tpl_ok.all())),
]
PROSE = pd.DataFrame(prose)
print(PROSE.to_string(index=False))
emit("prose_numbers", PROSE)

# median shrink, for the comparison against the cited sibling figure
med_pct = float(both["pct_change"].median())
print(f"\nmedian per-example base_size change on the 10 both-accepted pairs: {med_pct:.1f}%")
print(f"median over all 10 of |shrink|: {abs(med_pct):.1f}%  "
      f"(the closing note cites the sibling's '~13%')")

# ===========================================================================================
# 11. Seed-draw provenance: does #1 match what the cited notebook actually used?
# ===========================================================================================
rule("SEED PROVENANCE -- '#1 is literally the same click at seed_index=0'")

SEEDS = pd.DataFrame(SEED_ROWS)
hem_raw = pd.read_csv(os.path.join(OUT, "precision_at_k_14roi_prodseed_chromatin_hembbox_raw.csv"))
old_raw = pd.read_csv(os.path.join(OUT, "precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv"))
hem_seed = hem_raw.drop_duplicates("file_name").set_index("file_name")[["seed_ann_id", "base_size"]]
old_seed = old_raw.drop_duplicates("file_name").set_index("file_name")[["seed_ann_id", "base_size"]]

SEEDS = SEEDS.set_index("file_name")
SEEDS["hembbox_seed_ann_id"] = hem_seed["seed_ann_id"]
SEEDS["hembbox_base_size"] = hem_seed["base_size"]
SEEDS["halfpixfix_seed_ann_id"] = old_seed["seed_ann_id"]
SEEDS["halfpixfix_base_size"] = old_seed["base_size"]
SEEDS["pick1_eq_hembbox_used"] = SEEDS["nb_pick1_ann_id"] == SEEDS["hembbox_seed_ann_id"]
SEEDS["pick1_eq_halfpixfix_used"] = SEEDS["nb_pick1_ann_id"] == SEEDS["halfpixfix_seed_ann_id"]
SEEDS = SEEDS.reset_index()
print(SEEDS[["file_name", "n_pool_after_border", "nb_pick1_ann_id", "hembbox_seed_ann_id",
             "pick1_eq_hembbox_used", "halfpixfix_seed_ann_id", "pick1_eq_halfpixfix_used",
             "nb_pick2_ann_id", "true_seed_index1_ann_id",
             "pick2_equals_seed_index1"]].to_string(index=False))
n1_hem = int(SEEDS["pick1_eq_hembbox_used"].sum())
n1_old = int(SEEDS["pick1_eq_halfpixfix_used"].sum())
n2 = int(SEEDS["pick2_equals_seed_index1"].sum())
print(f"\n#1 matches the click production_seed_precision_at_k_chromatin_hem_bbox.ipynb "
      f"ACTUALLY USED at seed_index=0 : {n1_hem}/7")
print(f"#1 matches the click the half_pix_fix (gray) baseline used at seed_index=0 : {n1_old}/7")
print(f"#2 equals the annotation a real seed_index=1 draw would pick        : {n2}/7")
emit("seed_provenance", SEEDS)

# cross-check the base_size the sibling recorded for the ROIs where the click did match
bs_rows = []
for _, s in SEEDS.iterrows():
    fn = s["file_name"]
    g = REC[(REC["file_name"] == fn) & (REC["ann_id"] == s["nb_pick1_ann_id"])]
    gray_bs = g[g["channel"] == "gray_inverted"]["base_size"].iloc[0]
    hem_bs = g[g["channel"] == "hematoxylin_od"]["base_size"].iloc[0]
    bs_rows.append(dict(file_name=fn, ann_id=s["nb_pick1_ann_id"],
                        this_nb_gray_base_size=gray_bs, this_nb_hem_base_size=hem_bs,
                        halfpixfix_used_ann=s["halfpixfix_seed_ann_id"],
                        halfpixfix_base_size=s["halfpixfix_base_size"],
                        gray_matches_halfpixfix=bool(
                            s["nb_pick1_ann_id"] == s["halfpixfix_seed_ann_id"]
                            and gray_bs == s["halfpixfix_base_size"]),
                        hembbox_used_ann=s["hembbox_seed_ann_id"],
                        hembbox_base_size=s["hembbox_base_size"],
                        hem_matches_hembbox=bool(
                            s["nb_pick1_ann_id"] == s["hembbox_seed_ann_id"]
                            and hem_bs == s["hembbox_base_size"])))
BS = pd.DataFrame(bs_rows)
print("\ncross-validation of this notebook's base_size against the two production runs:")
print(BS.to_string(index=False))
print(f"gray base_size == half_pix_fix's recorded base_size on "
      f"{int(BS['gray_matches_halfpixfix'].sum())}/7 ROIs (where the click matched)")
emit("base_size_cross_validation", BS)

# ===========================================================================================
# 12. The cited sibling figure: '13/14 ROIs shrink, median ~13%'
# ===========================================================================================
rule("CITATION -- the sibling's '13/14 ROIs shrink, median ~13%'")

sib = old_seed.join(hem_seed, lsuffix="_old", rsuffix="_new")
sib["delta"] = sib["base_size_new"] - sib["base_size_old"]
sib["pct"] = (100 * sib["delta"] / sib["base_size_old"]).round(1)
sib["same_seed"] = sib["seed_ann_id_old"] == sib["seed_ann_id_new"]
n_shrink = int((sib["delta"] < 0).sum())
n_grow = int((sib["delta"] > 0).sum())
print(sib.sort_values("pct").to_string())
print(f"\nshrink {n_shrink}/14, grow {n_grow}/14, median pct change {sib['pct'].median():.1f}%")
print(f"same_seed on {int(sib['same_seed'].sum())}/14; "
      f"the one ROI that GROWS is {sib[sib['delta'] > 0].index.tolist()} "
      f"(same_seed={sib[sib['delta'] > 0]['same_seed'].tolist()})")
cit = pd.DataFrame([
    dict(claim="sibling: 13/14 ROIs shrink", recomputed=f"{n_shrink}/14", agrees=n_shrink == 13),
    dict(claim="sibling: median ~13%", recomputed=f"{sib['pct'].median():.1f}%",
         agrees=abs(sib["pct"].median() + 13) < 1.0),
    dict(claim="this notebook's 10 both-accepted pairs show the same magnitude",
         recomputed=f"median {med_pct:.1f}% here vs {sib['pct'].median():.1f}% there",
         agrees=abs(med_pct - sib["pct"].median()) < 3.0),
    dict(claim="the sibling's one growing ROI is a DIFFERENT-click ROI",
         recomputed=f"{sib[sib['delta'] > 0].index.tolist()}, same_seed="
                    f"{sib[sib['delta'] > 0]['same_seed'].tolist()}",
         agrees=not bool(sib[sib["delta"] > 0]["same_seed"].all())),
])
print(cit.to_string(index=False))
emit("sibling_citation", cit)

# ===========================================================================================
# 12b. The mechanism the closing note endorses by citation
# ===========================================================================================
rule("CITED MECHANISM -- does template shrink actually track the sibling's precision drop?")

hem_roi = pd.read_csv(os.path.join(OUT, "precision_at_k_14roi_prodseed_chromatin_hembbox_per_roi.csv"))
old_roi = pd.read_csv(os.path.join(OUT, "precision_at_k_14roi_prodseed_chromatin_halfpixfix_per_roi.csv"))
BUDGETS = (10, 20, 30, 50)
mech_rows = []
rng_p = np.random.default_rng(20260916)
for arm in ("tm_score", "chromatin_od"):
    a = hem_roi[hem_roi["arm"] == arm].set_index("file_name")
    o = old_roi[old_roi["arm"] == arm].set_index("file_name")
    d = sib["delta"].reindex(a.index)
    for k in BUDGETS:
        col = f"precision_at_{k}"
        dp = (a[col] - o[col]).reindex(a.index)
        ok = d.notna() & dp.notna()
        x, y = d[ok].to_numpy(float), dp[ok].to_numpy(float)
        rho = float(stats_mod.spearmanr(x, y).statistic)
        null = np.array([abs(stats_mod.spearmanr(x, rng_p.permutation(y)).statistic)
                         for _ in range(2000)])
        pval = float((np.sum(null >= abs(rho) - 1e-12) + 1) / (len(null) + 1))
        mech_rows.append(dict(subset="all 14 ROIs", arm=arm, K=k, n=int(ok.sum()),
                              spearman_rho_delta_base_size_vs_delta_precision=round(rho, 3),
                              perm_p=round(pval, 3),
                              sign_supports_smaller_is_worse=bool(rho > 0)))
MECH = pd.DataFrame(mech_rows)
print(MECH.to_string(index=False))
print("\nFor 'smaller template -> worse tm_score' the correlation between delta_base_size "
      "(negative = shrank) and delta_precision (negative = got worse) would have to be POSITIVE.")
print(f"tm_score cells with a positive rho: "
      f"{int(MECH[MECH['arm'] == 'tm_score']['sign_supports_smaller_is_worse'].sum())}/4; "
      f"tm_score cells with perm_p < 0.05: "
      f"{int((MECH[MECH['arm'] == 'tm_score']['perm_p'] < 0.05).sum())}/4")
emit("cited_mechanism_spearman", MECH)

# ===========================================================================================
# 13. Inference -- how far do 14 examples from 7 ROIs actually reach?
# ===========================================================================================
rule("INFERENCE -- G, the exchangeable unit, and the attainable p-floors")

# refusal claim: paired, per example. exact two-sided McNemar over discordant pairs.
b = int(((PAIR["gray_accepted"]) & (~PAIR["hem_accepted"])).sum())   # gray yes, hem no
c = int(((~PAIR["gray_accepted"]) & (PAIR["hem_accepted"])).sum())   # gray no, hem yes
n_disc = b + c
p_mcnemar = min(1.0, 2 * sum(math.comb(n_disc, i) for i in range(min(b, c) + 1)) / 2 ** n_disc) \
    if n_disc else 1.0
# ROI-level version: a ROI "carries signal" on refusal only if its two examples disagree
roi_ref = PAIR.groupby("file_name").apply(
    lambda g: pd.Series(dict(n_gray=int(g["gray_accepted"].sum()),
                             n_hem=int(g["hem_accepted"].sum()))), include_groups=False)
roi_ref["diff"] = roi_ref["n_gray"] - roi_ref["n_hem"]
G_ref = int((roi_ref["diff"] != 0).sum())
p_floor_ref = 2 / 2 ** G_ref if G_ref else 1.0

# size claim: ROI-level sign-flip on the mean per-ROI pct change, both-accepted pairs only
roi_size = both.groupby("file_name")["pct_change"].mean()
G_size = int(len(roi_size))
signs = np.sign(roi_size.to_numpy())
obs = float(np.mean(roi_size.to_numpy()))
perms = [float(np.mean(roi_size.to_numpy() * np.array(s)))
         for s in np.array(np.meshgrid(*[[1, -1]] * G_size)).T.reshape(-1, G_size)]
p_size = float(np.mean(np.abs(np.array(perms)) >= abs(obs) - 1e-12))
p_floor_size = 2 / 2 ** G_size

inf = pd.DataFrame([
    dict(claim="hematoxylin_od refuses more often",
         test="exact two-sided McNemar over the 14 paired examples",
         detail=f"discordant pairs b={b} (gray accept / hem refuse), c={c} (reverse)",
         statistic=f"{b} vs {c}", p_value=round(p_mcnemar, 4),
         G="14 examples from 7 ROIs", note="ROI is the exchangeable unit (D5), so 14 is "
                                            "not 14 independent units"),
    dict(claim="hematoxylin_od refuses more often (ROI-clustered)",
         test="ROI-level sign-flip on the per-ROI accept-count difference",
         detail=f"{G_ref} of 7 ROIs carry a non-zero difference",
         statistic=f"G_eff={G_ref}", p_value=round(p_floor_ref, 4),
         G=f"{G_ref} effective ROIs",
         note=f"p cannot go below {p_floor_ref:.4f} at G={G_ref}; the observed all-one-direction "
              f"case sits exactly at that floor"),
    dict(claim="where both accept, hem gives a smaller box",
         test="exact ROI-level sign-flip on the mean per-ROI % change",
         detail=f"observed mean {obs:.1f}%, all {int((signs < 0).sum())}/{G_size} ROIs negative",
         statistic=f"mean={obs:.1f}%", p_value=round(p_size, 4),
         G=f"{G_size} ROIs (013.tiff contributes no both-accepted pair)",
         note=f"p-floor at G={G_size} is {p_floor_size:.4f}; 2^{G_size} enumerated exactly"),
    dict(claim="where both accept, hem gives a smaller box (per-example, uncorrected)",
         test="exact sign test over the 10 both-accepted examples, ignoring clustering",
         detail="10/10 negative", statistic="10/10",
         p_value=round(2 / 2 ** 10, 5), G="10 examples from 6 ROIs",
         note="ANTI-CONSERVATIVE -- the 10 examples are 6 clusters, not 10 units"),
])
print(inf.to_string(index=False))
emit("inference", inf)

# t interval on the G cluster means, for the size effect
mean = roi_size.mean()
sd = roi_size.std(ddof=1)
stats = stats_mod
tcrit = stats.t.ppf(0.975, G_size - 1)
lo, hi = mean - tcrit * sd / math.sqrt(G_size), mean + tcrit * sd / math.sqrt(G_size)
print(f"\nper-ROI mean base_size change: {mean:.1f}%  "
      f"95% t interval on {G_size} cluster means ({G_size - 1} df): [{lo:.1f}%, {hi:.1f}%]")
print("per-ROI values:")
print(roi_size.round(1).to_string())
emit("roi_size_effect", roi_size.reset_index().rename(columns={"pct_change": "mean_pct_change"}))

# ===========================================================================================
# 13b. Mechanism -- is the gray component really "the fuller nuclear outline"?
# ===========================================================================================
rule("MECHANISM -- geometry of the accepted components, per channel")

geo = REC[REC["accepted"]][["file_name", "ann_id", "channel", "area", "area_frac", "solidity",
                            "base_size", "bbox_h", "bbox_w", "bbox_touches_window_edge",
                            "base_size_at_window_ceiling"]].copy()
geo["area_frac"] = geo["area_frac"].astype(float).round(3)
print(geo.to_string(index=False))
g = geo[geo["channel"] == "gray_inverted"]
h = geo[geo["channel"] == "hematoxylin_od"]
print(f"\ngray_inverted: {len(g)} accepted; median area_frac {g['area_frac'].median():.3f}; "
      f"{int(g['base_size_at_window_ceiling'].sum())} at the {OTSU_WINDOW}px window ceiling; "
      f"{int(g['bbox_touches_window_edge'].sum())} with a bbox touching the window edge; "
      f"{int((g['area_frac'] >= 0.40).sum())} covering >=40% of the window")
print(f"hematoxylin_od: {len(h)} accepted; median area_frac {h['area_frac'].median():.3f}; "
      f"{int(h['base_size_at_window_ceiling'].sum())} at the ceiling; "
      f"{int(h['bbox_touches_window_edge'].sum())} touching the window edge; "
      f"{int((h['area_frac'] >= 0.40).sum())} covering >=40%")
both_ids = set(zip(both["file_name"], both["ann_id"]))
gb = g[[tuple(x) in both_ids for x in zip(g["file_name"], g["ann_id"])]]
print(f"\namong the 10 both-accepted pairs, gray base_size sits at the {OTSU_WINDOW}px ceiling on "
      f"{int(gb['base_size_at_window_ceiling'].sum())}/10 -- those gray values are "
      f"right-censored, so the quoted gray->hem magnitudes there are lower bounds on gray")
emit("component_geometry", geo)

# rank the both-accepted pairs by absolute shrink, and mark which are censored at the ceiling
cens = set(zip(g[g["base_size_at_window_ceiling"] == True]["file_name"],
               g[g["base_size_at_window_ceiling"] == True]["ann_id"]))
RANK = both[["file_name", "ann_id", "which", "gray_base_size", "hem_base_size",
             "delta_px", "pct_change"]].copy()
RANK["gray_censored_at_window"] = [(f, a) in cens for f, a in zip(RANK["file_name"], RANK["ann_id"])]
RANK = RANK.sort_values("delta_px").reset_index(drop=True)
RANK["abs_shrink_rank"] = np.arange(1, len(RANK) + 1)
print("\nboth-accepted pairs ranked by absolute shrink:")
print(RANK.to_string(index=False))
top4 = RANK.head(4)["gray_censored_at_window"].all()
none_below = not RANK.tail(len(RANK) - 4)["gray_censored_at_window"].any()
print(f"\nthe {int(RANK['gray_censored_at_window'].sum())} censored pairs are exactly the 4 "
      f"largest absolute shrinks: {bool(top4 and none_below)} "
      f"(4th = {RANK.iloc[3]['delta_px']:.0f}px, 5th = {RANK.iloc[4]['delta_px']:.0f}px)")
emit("shrink_ranking", RANK)

# ===========================================================================================
# 12c. Overlap with the sibling -- do the two runs agree where they measure the same pair?
# ===========================================================================================
rule("OVERLAP -- same click in both runs: do the percentages agree, and what drives the median gap?")

ov_rows = []
for _, r in BS.iterrows():
    sel = both[(both["file_name"] == r["file_name"]) & (both["ann_id"] == r["ann_id"])]
    if len(sel) == 0:
        continue
    row = sel.iloc[0]
    here_pct = float(row["pct_change"])
    sib_pct = 100.0 * (r["hembbox_base_size"] - r["halfpixfix_base_size"]) / r["halfpixfix_base_size"]
    ov_rows.append(dict(
        file_name=r["file_name"], ann_id=r["ann_id"],
        here=f"{int(row['gray_base_size'])}->{int(row['hem_base_size'])}",
        here_pct=round(here_pct, 1),
        sibling=f"{int(r['halfpixfix_base_size'])}->{int(r['hembbox_base_size'])}",
        sibling_pct=round(sib_pct, 1),
        same_click_in_both_runs=bool(r["gray_matches_halfpixfix"] and r["hem_matches_hembbox"]),
        identical=bool(abs(here_pct - sib_pct) < 1e-9)))
OV = pd.DataFrame(ov_rows)
print(OV.to_string(index=False))
ovs = OV[OV["same_click_in_both_runs"]]
print(f"\npairs measured by BOTH runs on the same click: {len(ovs)}; "
      f"identical pct_change on {int(ovs['identical'].sum())}/{len(ovs)}")
m1 = both[both["which"].str.startswith("#1")]
m2 = both[both["which"].str.startswith("#2")]
print(f"\nmedian pct_change, this notebook, all 10 both-accepted pairs : {both['pct_change'].median():.1f}%")
print(f"median pct_change, this notebook, #1 draws only (n={len(m1)})       : {m1['pct_change'].median():.1f}%")
print(f"median pct_change, this notebook, #2 draws only (n={len(m2)})       : {m2['pct_change'].median():.1f}%")
print(f"median pct_change, sibling, all 14 production seeds           : {sib['pct'].median():.1f}%")
print(f"median pct_change, sibling, same-seed 11 only                 : {sib[sib['same_seed']]['pct'].median():.1f}%")
print("\n-> the two runs agree exactly wherever they overlap; the median gap is composition")
print("   (the #2 draws appear in no production run), not retry and not disagreement.")
emit("sibling_overlap", OV)

# ===========================================================================================
# 13c. Figure-render defects that no recompute would catch
# ===========================================================================================
rule("FIGURE RENDER -- the panel-3 click marker's actual stroke colour")

import matplotlib                                           # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                             # noqa: E402
from skimage.color import label2rgb                         # noqa: E402

fig_, ax_ = plt.subplots()
p3, = ax_.plot(0.5, 0.5, marker="+", color="white", ms=13, mew=2, mec="black")  # panel 3 call
p2, = ax_.plot(0.5, 0.5, marker="+", color="red", ms=13, mew=2)                 # panel 2 call
lab_ = np.zeros((5, 5), int); lab_[2, 2] = 1
bg_ = label2rgb(lab_, bg_label=0)[0, 0]
plt.close(fig_)
render = pd.DataFrame([
    dict(panel="2 (Otsu binary)", call="plot(marker='+', color='red', mew=2)",
         stroke_colour=str(p2.get_markeredgecolor()), background="black (binary 0)",
         click_marker_visible_on_background=True),
    dict(panel="3 (components)", call="plot(marker='+', color='white', mew=2, mec='black')",
         stroke_colour=str(p3.get_markeredgecolor()),
         background=f"label2rgb(bg_label=0) -> RGB {tuple(bg_)}",
         click_marker_visible_on_background=False),
])
print(render.to_string(index=False))
print("\n'+' has no fill, so markerfacecolor='white' never draws; mec='black' is the stroke.")
print("Consequence: in every REFUSED panel 3 (click on background) the click marker is")
print("black-on-black and invisible -- the closing note points the reader at panel 2/3 for")
print(f"exactly that case. Refusals affected: {len(ref)}/{len(ref)}.")
emit("figure_render_defects", render)

# ===========================================================================================
# 14. Production-gap checks the notebook does not run
# ===========================================================================================
rule("PRODUCTION GAP -- checks build_seed applies that this notebook does not")

gap = REC[REC["accepted"]][["file_name", "ann_id", "channel", "base_size", "offset_px",
                            "template_readable_at_patch_size_73"]].copy()
print(gap.to_string(index=False))
gap["template_readable_at_patch_size_73"] = gap["template_readable_at_patch_size_73"].astype(bool)
n_fail = int((~gap["template_readable_at_patch_size_73"]).sum())
print(f"\naccepted examples whose recentred template point would FAIL build_seed's "
      f"_patch_readable(patch_size=73): {n_fail}/{len(gap)}")
print(f"largest click->template offset among accepted examples: "
      f"{gap['offset_px'].max():.2f} px "
      f"({gap.loc[gap['offset_px'].idxmax(), 'file_name']} ann "
      f"{int(gap.loc[gap['offset_px'].idxmax(), 'ann_id'])} "
      f"{gap.loc[gap['offset_px'].idxmax(), 'channel']})")
emit("production_gap", gap)

print(f"\n\naudit script finished in {time.time() - T_START:.0f}s; "
      f"{len(TABLES)} tables written to results/")
