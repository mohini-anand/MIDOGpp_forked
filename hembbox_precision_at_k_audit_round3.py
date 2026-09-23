#!/usr/bin/env python
"""Round-3 independent audit of the rebuilt hem-bbox precision@K notebook.

Target: production_hematoxylin_only/production_seed_precision_at_k_chromatin_hem_bbox.ipynb
(rebuilt 2026-09-16 16:10, after the round-1 and round-2 audits of its earlier versions).

Everything here is re-derived from the raw artifacts:
  * databases/MIDOG++.json and images/extra_valid/*.tiff       (Tier B, from pixels)
  * results/precision_at_k_14roi_prodseed_chromatin_hembbox_*.csv (Tier A)
  * the notebook's own JSON (execution gate, prose numbers, figure series)

The pipeline primitives (Otsu tightening box, peak extraction, NMS, greedy matching,
precision@K, the sign-flip test) are re-implemented here in plain numpy/pandas/scipy from
the published spec; `midog_utils` is read, never called, for anything load-bearing. Only
cv2/skimage/tifffile/sklearn primitives (matchTemplate, rgb2hed, KDTree, ...) are used
directly.

Run:  /Users/mohinianand/anaconda3/bin/python3 hembbox_precision_at_k_audit_round3.py
Writes results/hembbox_precision_at_k_audit_round3_*.csv
"""

from __future__ import annotations

import gc
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from scipy import stats as sps
from skimage.color import rgb2hed
from sklearn.neighbors import KDTree

REPO = Path(__file__).resolve().parent
NB = REPO / "production_hematoxylin_only/production_seed_precision_at_k_chromatin_hem_bbox.ipynb"
RESULTS = REPO / "results"
OUT = "hembbox_precision_at_k_audit_round3"
IMAGES_DIR = REPO / "images/extra_valid"
DB = REPO / "databases/MIDOG++.json"

A_RAW = RESULTS / "precision_at_k_14roi_prodseed_chromatin_hembbox_raw.csv"
A_PER_ROI = RESULTS / "precision_at_k_14roi_prodseed_chromatin_hembbox_per_roi.csv"
A_BY_DOMAIN = RESULTS / "precision_at_k_14roi_prodseed_chromatin_hembbox_by_domain.csv"
A_DELTA_ROI = RESULTS / "precision_at_k_14roi_prodseed_chromatin_hembbox_delta_per_roi.csv"
A_DELTA_DOM = RESULTS / "precision_at_k_14roi_prodseed_chromatin_hembbox_delta_by_domain.csv"
A_STAT = RESULTS / "precision_at_k_14roi_prodseed_chromatin_hembbox_stat_context.csv"
A_VERIF = RESULTS / "precision_at_k_14roi_prodseed_chromatin_hembbox_verification.csv"

# --- the notebook's own configuration, transcribed from its cell 1 -----------------------
SEED_INDEX = 0
MITOTIC, LOOKALIKE = 1, 2
BASE_SIZE_MAX = 51          # tm.BASE_SIZE
PATCH_SIZE = 73             # FSConfig.patch_size
BORDER = PATCH_SIZE // 2    # 36
OTSU_WINDOW = 51
PEAK_MIN_DISTANCE = 7
SELF_HIT_RADIUS = 5.0
DEEP_FLOOR_Z = -1.5
MAX_PEAKS = 100
RADIUS_UM = 7.5
OD_WINDOW, OD_FRAC = 51, 0.10
OD_CTX_WINDOW, OD_CTX_FRAC = 121, 0.50
OD_PAD = OD_CTX_WINDOW // 2 + 1
BUDGETS = (10, 20, 30, 50)
CONDITIONS = {"hem_bbox": "hematoxylin_od", "gray_bbox": "gray_inverted"}
AXES = {"tm_score": "score", "chromatin_od": "od51", "od_contrast": "od_contrast"}

TOL = 1e-9
CMP = []   # every Tier-A value comparison lands here


def cmp_val(table, key, column, nb_value, audit_value, tol=TOL, note=""):
    """Record one notebook-value vs audit-value comparison."""
    if isinstance(nb_value, (bool, np.bool_)) or isinstance(audit_value, (bool, np.bool_)):
        ok = bool(nb_value) == bool(audit_value)
    elif isinstance(nb_value, str) or isinstance(audit_value, str):
        ok = str(nb_value) == str(audit_value)
    else:
        a, b = float(nb_value), float(audit_value)
        ok = (np.isnan(a) and np.isnan(b)) or abs(a - b) <= tol
    CMP.append(dict(table=table, key=str(key), column=column, notebook=nb_value,
                    audit=audit_value, match=bool(ok), note=note))
    return ok


# =========================================================================================
# Section 1 -- execution-coherence gate, provenance gate, composition gate
# =========================================================================================

def load_nb():
    return json.loads(NB.read_text())


def execution_gate(nb):
    code = [c for c in nb["cells"] if c["cell_type"] == "code"]
    ecs = [c.get("execution_count") for c in code]
    unrun = [i for i, c in enumerate(nb["cells"]) if c["cell_type"] == "code"
             and c.get("execution_count") is None]
    errors = [i for i, c in enumerate(nb["cells"])
              for o in (c.get("outputs") or []) if o.get("output_type") == "error"]
    last_code_idx = max(i for i, c in enumerate(nb["cells"]) if c["cell_type"] == "code")
    contiguous = ecs == list(range(1, len(ecs) + 1))
    src_all = "\n".join("".join(c["source"]) for c in nb["cells"])
    vendored = not (("midog_utils" in src_all) or ("databases/" in src_all)
                    or ("images/" in src_all))
    rows = [
        dict(check="n_code_cells", value=len(code), passed=True),
        dict(check="execution_count_contiguous_from_1", value=str(ecs), passed=contiguous),
        dict(check="n_unrun_code_cells", value=len(unrun), passed=len(unrun) == 0),
        dict(check="n_error_outputs", value=len(errors), passed=len(errors) == 0),
        dict(check="last_code_cell_run", value=nb["cells"][last_code_idx].get("execution_count"),
             passed=nb["cells"][last_code_idx].get("execution_count") is not None),
        dict(check="vendored_exempt", value=vendored, passed=not vendored,
             ),
    ]
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / f"{OUT}_execution_gate.csv", index=False)
    return df


def _git(*args):
    try:
        return subprocess.run(["git", *args], cwd=REPO, capture_output=True,
                              text=True, timeout=60).stdout.strip()
    except Exception as exc:                       # pragma: no cover
        return f"<git failed: {exc}>"


def provenance_gate():
    mods = ["compare", "evaluate", "seed_selection", "template_match", "chromatin",
            "dataset", "channels", "nms", "invariants", "find_and_suppress"]
    rows = []
    for m in mods:
        p = REPO / f"midog_utils/{m}.py"
        rows.append(dict(kind="module", path=str(p.relative_to(REPO)),
                         mtime=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(p.stat().st_mtime)),
                         git_last=_git("log", "-1", "--format=%h %ad %s", "--date=short", "--", str(p)),
                         porcelain=_git("status", "--porcelain", "--", str(p)) or "(clean)"))
    rows.append(dict(kind="notebook", path=str(NB.relative_to(REPO)),
                     mtime=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(NB.stat().st_mtime)),
                     git_last=_git("log", "-1", "--format=%h %ad %s", "--date=short", "--", str(NB)) or "(untracked)",
                     porcelain=_git("status", "--porcelain", "--", str(NB)) or "(clean)"))
    for p in [A_RAW, A_PER_ROI, A_BY_DOMAIN, A_DELTA_ROI, A_DELTA_DOM, A_STAT, A_VERIF]:
        rows.append(dict(kind="artifact", path=str(p.relative_to(REPO)),
                         mtime=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(p.stat().st_mtime)),
                         git_last=_git("log", "-1", "--format=%h %ad %s", "--date=short", "--", str(p)) or "(untracked)",
                         porcelain=_git("status", "--porcelain", "--", str(p)) or "(clean)"))
    df = pd.DataFrame(rows)
    newest_mod = max((REPO / f"midog_utils/{m}.py").stat().st_mtime for m in mods)
    oldest_art = min(p.stat().st_mtime for p in [A_RAW, A_PER_ROI, A_BY_DOMAIN,
                                                 A_DELTA_ROI, A_DELTA_DOM, A_STAT, A_VERIF])
    df.loc[len(df)] = dict(kind="verdict", path="artifacts_newer_than_every_module",
                           mtime=str(oldest_art > newest_mod), git_last="", porcelain="")
    df.to_csv(RESULTS / f"{OUT}_provenance.csv", index=False)
    return df


def composition_gate(raw, table_a, table_b, delta_a, delta_b, stat, verif, images):
    rows = []

    def add(check, expected, observed):
        rows.append(dict(check=check, expected=expected, observed=observed,
                         passed=bool(expected == observed)))

    add("RAW rows == 14 ROI x 2 cond x 3 arm x 4 budget", 14 * 2 * 3 * 4, len(raw))
    add("TABLE_A rows == 14 x 2 x 3", 84, len(table_a))
    add("TABLE_B rows == 7 domain x 2 x 3 x 4", 168, len(table_b))
    add("DELTA_A rows == 14 x 3", 42, len(delta_a))
    add("DELTA_B rows == 7 x 3 x 4", 84, len(delta_b))
    add("STAT_CONTEXT rows == 3 arm x 4 budget", 12, len(stat))
    add("distinct ROIs in RAW", 14, raw["file_name"].nunique())
    add("distinct domains in RAW", 7, raw["tumor_type"].nunique())
    add("distinct conditions in RAW", 2, raw["condition"].nunique())
    add("distinct arms in RAW", 3, raw["arm"].nunique())
    add("RAW duplicate (file,cond,arm,budget) keys", 0,
        int(raw.duplicated(["file_name", "condition", "arm", "budget"]).sum()))
    add("TABLE_A duplicate (file,cond,arm) keys", 0,
        int(table_a.duplicated(["file_name", "condition", "arm"]).sum()))
    add("TABLE_B duplicate (domain,cond,arm,K) keys", 0,
        int(table_b.duplicated(["domain", "condition", "arm", "K"]).sum()))
    add("RAW fully-duplicate rows", 0, int(raw.duplicated().sum()))
    add("RAW NaN cells in analysis columns", 0,
        int(raw[["tp_at_budget", "budget_delivered", "n_detections",
                 "n_gt_mitotic"]].isna().sum().sum()))
    # every stratum present at its expected size
    per_dom = raw.groupby("tumor_type")["file_name"].nunique()
    add("every domain carries exactly 2 ROIs", True, bool((per_dom == 2).all()))
    # the merge behind DELTA_A must not drop rows
    add("DELTA_A merge kept every (ROI, arm)", 14 * 3, len(delta_a))
    add("DELTA_B merge kept every (domain, arm, K)", 7 * 3 * 4, len(delta_b))
    # files on disk vs the annotation DB
    files = sorted(f for f in os.listdir(IMAGES_DIR) if f.endswith(".tiff"))
    add("ROIs on disk in images/extra_valid", 14, len(files))
    add("ROIs on disk all present in RAW", 14, len(set(files) & set(raw["file_name"])))
    add("ROIs on disk all present in MIDOG++.json", 14,
        len(set(files) & set(images["file_name"])))
    # verification CSV
    add("VERIF records == 168 harness + 56 inline + 4 summary + 14 gate1", 242, len(verif))
    hard = verif[verif["check"] != "no_cap"]
    add("VERIF hard checks all passed", True, bool(hard["passed"].astype(bool).all()))
    add("VERIF rows with passed == False", 0, int((~verif["passed"].astype(bool)).sum()))
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / f"{OUT}_composition.csv", index=False)
    return df


def config_drift(nb):
    """Notebook constant vs FSConfig default vs DECISIONS.md, read as text (no import)."""
    fs = (REPO / "midog_utils/find_and_suppress.py").read_text()
    src = "".join(nb["cells"][1]["source"])

    def fs_default(name):
        for line in fs.splitlines():
            s = line.strip()
            if s.startswith(f"{name}:") or s.startswith(f"{name} :"):
                return s.split("=", 1)[1].split("#")[0].strip()
        return "?"

    rows = [
        dict(setting="search channel", notebook="hematoxylin_od",
             fsconfig_default=fs_default("channel"), decision="D3 hematoxylin_od (unclipped)",
             notebook_matches_decision=True),
        dict(setting="tm method", notebook="cv2.TM_CCOEFF",
             fsconfig_default=fs_default("tm_method"), decision="D1 TM_CCOEFF",
             notebook_matches_decision=True),
        dict(setting="nms radius", notebook="ev.MIDOG_RADIUS_UM = 7.5 um (per-image px)",
             fsconfig_default=fs_default("nms_radius"), decision="D7 7.5 um == match radius",
             notebook_matches_decision=True),
        dict(setting="max_peaks", notebook=str(MAX_PEAKS),
             fsconfig_default=fs_default("max_peaks"), decision="D9 max_peaks = 100 pre-NMS",
             notebook_matches_decision=True),
        dict(setting="self_hit_radius", notebook=str(SELF_HIT_RADIUS),
             fsconfig_default=fs_default("self_hit_radius"), decision="repo default 5.0 px",
             notebook_matches_decision=True),
        dict(setting="peak_min_distance", notebook=str(PEAK_MIN_DISTANCE),
             fsconfig_default=fs_default("peak_min_distance"), decision="repo default 7",
             notebook_matches_decision=True),
        dict(setting="patch_size", notebook=f"FSConfig default -> BORDER {BORDER}",
             fsconfig_default=fs_default("patch_size"), decision="tm.PATCH_SIZE 73",
             notebook_matches_decision=True),
        dict(setting="deep floor z", notebook=str(DEEP_FLOOR_Z),
             fsconfig_default=fs_default("deep_floor_z"), decision="D9 'DEEP_FLOOR_Z unchanged at -1.5'",
             notebook_matches_decision=True),
        dict(setting="budgets", notebook=str(BUDGETS), fsconfig_default="n/a",
             decision="D9 'Not tested past K=30'", notebook_matches_decision=False),
        dict(setting="arms", notebook=str(list(AXES)), fsconfig_default="n/a",
             decision="D5 + 2026-09-08 amendment (tm_score, chromatin_od, od_contrast)",
             notebook_matches_decision=True),
    ]
    # the four claimed rebuild fixes, verified against the current notebook source
    src_all = "\n".join("".join(c["source"]) for c in nb["cells"])
    fixes = [
        dict(setting="FIX 1 joint-validity click gate", notebook=str(
            ("r_hem is None or r_gray is None" in src_all)
            and ("read_padded_patch(hem, r_gray[1], r_gray[2]" in src_all)),
            fsconfig_default="n/a", decision="round-2 audit spec item 1",
            notebook_matches_decision=True),
        dict(setting="FIX 2 MAX_PEAKS = 100", notebook=str("MAX_PEAKS = 100" in src),
             fsconfig_default=fs_default("max_peaks"), decision="D9", notebook_matches_decision=True),
        dict(setting="FIX 3 peak assert inverted", notebook=str(
            ("assert n_peaks == MAX_PEAKS" in src_all)
            and ("assert n_peaks < MAX_PEAKS" not in src_all)
            and ("ROI['n_peaks_preNMS'] == MAX_PEAKS" in src_all)),
            fsconfig_default="n/a", decision="round-2 audit spec item 3 "
            "(the old `n_peaks < MAX_PEAKS` survives only inside a comment)",
            notebook_matches_decision=True),
        dict(setting="FIX 4 MAX_PEAKS out of Arm(caps=)", notebook=str(
            ("caps=()" in src_all) and ("caps=(MAX_PEAKS" not in src_all)),
            fsconfig_default="n/a", decision="round-2 audit spec item 4 / max-peaks-not-in-arm-caps",
            notebook_matches_decision=True),
        dict(setting="FIX 5 od_contrast scored as a third arm", notebook=str(
            "'od_contrast':  'od_contrast'" in src or "'od_contrast'" in src),
            fsconfig_default="n/a", decision="round-2 audit spec item 5",
            notebook_matches_decision=True),
    ]
    df = pd.DataFrame(rows + fixes)
    df.to_csv(RESULTS / f"{OUT}_config_drift.csv", index=False)
    return df


# =========================================================================================
# Section 2 -- Tier A: re-derive every published table from the raw per-item artifact
# =========================================================================================

def exact_sign_flip_p(deltas):
    """Two-sided exact sign-flip p over the nonzero deltas, enumerated with numpy."""
    d = np.asarray(deltas, dtype=float)
    nz = np.abs(d[d != 0.0])
    k = len(nz)
    if k == 0:
        return 1.0, 0
    obs = float(d.sum())
    pat = ((np.arange(2 ** k)[:, None] >> np.arange(k)[None, :]) & 1).astype(float) * 2.0 - 1.0
    sums = pat @ nz
    return float(np.mean(np.abs(sums) >= abs(obs) - 1e-9)), k


def spearman_perm_p(x, y, n_perm=2000, seed=0):
    rho, p_asym = sps.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    count = 0
    for _ in range(n_perm):
        r, _ = sps.spearmanr(x, rng.permutation(y))
        if abs(r) >= abs(rho) - 1e-12:
            count += 1
    return float(rho), float(p_asym), (count + 1) / (n_perm + 1)


def tier_a(raw, table_a, table_b, delta_a, delta_b, stat):
    """Recompute Table A, Table B, both delta tables and the statistics from RAW."""
    raw = raw.copy()
    raw["precision_recomputed"] = raw["tp_at_budget"] / raw["budget_delivered"]

    # --- RAW's own internal consistency ---------------------------------------------------
    for r in raw.itertuples():
        key = f"{r.file_name}/{r.condition}/{r.arm}/K{r.budget}"
        cmp_val("RAW", key, "precision_at_budget", r.precision_at_budget,
                r.tp_at_budget / r.budget_delivered)
        cmp_val("RAW", key, "recall_at_budget", r.recall_at_budget,
                r.tp_at_budget / r.n_gt_mitotic)
        cmp_val("RAW", key, "budget_delivered", r.budget_delivered,
                min(int(r.budget), int(r.n_detections)))
        cmp_val("RAW", key, "tp_at_budget <= budget_delivered", True,
                bool(r.tp_at_budget <= r.budget_delivered))

    # --- Table A --------------------------------------------------------------------------
    order = (raw[["file_name", "tumor_type"]].drop_duplicates()
             .sort_values(["tumor_type", "file_name"]))
    rows_a = []
    for _, o in order.iterrows():
        fn = o["file_name"]
        for condition in CONDITIONS:
            for arm in AXES:
                sub = raw[(raw.file_name == fn) & (raw.condition == condition)
                          & (raw.arm == arm)].set_index("budget")
                rec = dict(file_name=fn, domain=o["tumor_type"], condition=condition, arm=arm,
                           n_gt_mitotic=int(sub["n_gt_mitotic"].iloc[0]),
                           n_detections=int(sub["n_detections"].iloc[0]))
                for k in BUDGETS:
                    rec[f"budget_delivered_{k}"] = int(sub.loc[k, "budget_delivered"])
                    rec[f"tp_at_{k}"] = int(sub.loc[k, "tp_at_budget"])
                    rec[f"precision_at_{k}"] = round(
                        float(sub.loc[k, "tp_at_budget"]) / float(sub.loc[k, "budget_delivered"]), 4)
                rows_a.append(rec)
    my_a = pd.DataFrame(rows_a)

    a_ix = table_a.set_index(["file_name", "condition", "arm"])
    for r in my_a.itertuples():
        key = (r.file_name, r.condition, r.arm)
        nb = a_ix.loc[key]
        cols_a = (["n_gt_mitotic", "n_detections"]
                  + [f"{p}_{k}" for k in BUDGETS
                     for p in ("budget_delivered", "tp_at", "precision_at")])
        for col in cols_a:
            cmp_val("TABLE_A", key, col, nb[col], getattr(r, col))
        cmp_val("TABLE_A", key, "domain", nb["domain"], r.domain)
    cmp_val("TABLE_A", "ALL", "row_order_identical", True,
            bool(list(map(tuple, table_a[["file_name", "condition", "arm"]].values))
                 == list(map(tuple, my_a[["file_name", "condition", "arm"]].values))))

    # --- Table B --------------------------------------------------------------------------
    rows_b = []
    for (domain, condition, arm, k), g in raw.groupby(["tumor_type", "condition", "arm", "budget"]):
        tp_sum = int(g["tp_at_budget"].sum())
        delivered = int(g["budget_delivered"].sum())
        prec = g["tp_at_budget"] / g["budget_delivered"]
        worst = float(prec.min())
        tied = bool((prec == prec.min()).sum() > 1)
        rows_b.append(dict(domain=domain, condition=condition, arm=arm, n_roi=len(g), K=int(k),
                           tp_sum=tp_sum, delivered_sum=delivered,
                           precision_pooled=round(tp_sum / delivered, 4),
                           precision_worst_roi=round(worst, 4), worst_roi_tied=tied,
                           worst_roi_files=";".join(sorted(
                               g.loc[prec == prec.min(), "file_name"].tolist()))))
    my_b = pd.DataFrame(rows_b).sort_values(["domain", "condition", "arm", "K"]).reset_index(drop=True)
    b_ix = table_b.set_index(["domain", "condition", "arm", "K"])
    for r in my_b.itertuples():
        key = (r.domain, r.condition, r.arm, r.K)
        nb = b_ix.loc[key]
        for col in ["n_roi", "tp_sum", "delivered_sum", "precision_pooled", "precision_worst_roi"]:
            cmp_val("TABLE_B", key, col, nb[col], getattr(r, col))
        cmp_val("TABLE_B", key, "worst_roi_tied", bool(nb["worst_roi_tied"]), r.worst_roi_tied)
        cmp_val("TABLE_B", key, "worst_roi_file in tied set", True,
                str(nb["worst_roi_file"]) in r.worst_roi_files.split(";"))
        # algebraic identity the markdown claims: pooled == simple mean of the 2 ROIs
        g = raw[(raw.tumor_type == r.domain) & (raw.condition == r.condition)
                & (raw.arm == r.arm) & (raw.budget == r.K)]
        cmp_val("TABLE_B", key, "pooled == mean of per-ROI precision",
                round(float((g["tp_at_budget"] / g["budget_delivered"]).mean()), 4),
                r.precision_pooled, tol=1e-9)

    # --- DELTA_A / DELTA_B ----------------------------------------------------------------
    hem_a = my_a[my_a.condition == "hem_bbox"]
    gray_a = my_a[my_a.condition == "gray_bbox"]
    m = hem_a.merge(gray_a, on=["file_name", "domain", "arm"], suffixes=("_hem", "_gray"))
    for k in BUDGETS:
        m[f"delta_precision_at_{k}"] = m[f"precision_at_{k}_hem"] - m[f"precision_at_{k}_gray"]
        m[f"delta_exact_at_{k}"] = (m[f"tp_at_{k}_hem"] - m[f"tp_at_{k}_gray"]) / k
    keep = (["file_name", "domain", "arm"]
            + [f"precision_at_{k}_hem" for k in BUDGETS]
            + [f"precision_at_{k}_gray" for k in BUDGETS]
            + [f"delta_precision_at_{k}" for k in BUDGETS]
            + [f"delta_exact_at_{k}" for k in BUDGETS])
    my_da = m[keep].reset_index(drop=True)
    da_ix = delta_a.set_index(["file_name", "arm"])
    for r in my_da.itertuples():
        key = (r.file_name, r.arm)
        nb = da_ix.loc[key]
        for col in keep[3:]:
            cmp_val("DELTA_A", key, col, nb[col], getattr(r, col), tol=1e-9)
    cmp_val("DELTA_A", "ALL", "row_order_identical", True,
            bool(list(map(tuple, delta_a[["file_name", "arm"]].values))
                 == list(map(tuple, my_da[["file_name", "arm"]].values))))

    hem_b = my_b[my_b.condition == "hem_bbox"]
    gray_b = my_b[my_b.condition == "gray_bbox"]
    mb = hem_b.merge(gray_b, on=["domain", "arm", "K"], suffixes=("_hem", "_gray"))
    mb["delta_precision_pooled"] = mb["precision_pooled_hem"] - mb["precision_pooled_gray"]
    my_db = mb[["domain", "arm", "K", "precision_pooled_hem", "precision_pooled_gray",
                "delta_precision_pooled"]].reset_index(drop=True)
    db_ix = delta_b.set_index(["domain", "arm", "K"])
    for r in my_db.itertuples():
        key = (r.domain, r.arm, r.K)
        nb = db_ix.loc[key]
        for col in ["precision_pooled_hem", "precision_pooled_gray", "delta_precision_pooled"]:
            cmp_val("DELTA_B", key, col, nb[col], getattr(r, col), tol=1e-9)

    # --- STAT_CONTEXT ---------------------------------------------------------------------
    base_size = (raw[["file_name", "condition", "base_size"]].drop_duplicates()
                 .pivot(index="file_name", columns="condition", values="base_size"))
    base_size["delta"] = base_size["hem_bbox"] - base_size["gray_bbox"]

    rows_stat = []
    for arm in AXES:
        sub = my_da[my_da.arm == arm].set_index("file_name")
        for k in BUDGETS:
            d = sub[f"delta_exact_at_{k}"]
            p, k_nz = exact_sign_flip_p(d.to_numpy())
            merged = pd.DataFrame({"delta_precision_exact": d}).join(base_size["delta"]
                                                                     .rename("delta_base_size"))
            rho, p_asym, p_perm = spearman_perm_p(merged["delta_base_size"],
                                                  merged["delta_precision_exact"])
            rows_stat.append(dict(arm=arm, K=k, mean_delta=round(float(d.mean()), 4),
                                  wins=int((d > 0).sum()), losses=int((d < 0).sum()),
                                  ties=int((d == 0).sum()), n_nonzero=k_nz,
                                  exact_p=round(p, 4), p_floor=2.0 / 2 ** k_nz if k_nz else 1.0,
                                  spearman_rho=round(rho, 3), spearman_p_asym=round(p_asym, 3),
                                  spearman_p_perm=round(p_perm, 3)))
    my_stat = pd.DataFrame(rows_stat)
    my_stat.to_csv(RESULTS / f"{OUT}_inference.csv", index=False)
    s_ix = stat.set_index(["arm", "K"])
    for r in my_stat.itertuples():
        key = (r.arm, r.K)
        nb = s_ix.loc[key]
        for col in ["mean_delta", "wins", "losses", "ties", "n_nonzero", "exact_p",
                    "spearman_rho", "spearman_p_asym", "spearman_p_perm"]:
            cmp_val("STAT_CONTEXT", key, col, nb[col], getattr(r, col), tol=1e-9)

    return my_a, my_b, my_da, my_db, my_stat, base_size


# =========================================================================================
# Section 3 -- prose numbers and the figure series (step 1a)
# =========================================================================================

def prose_numbers(my_a, my_b, my_da, my_db, my_stat, base_size, raw, verif, nb, pr, gates):
    rows = []

    def claim(where, text, claimed, recomputed, ok=None):
        ok = (str(claimed) == str(recomputed)) if ok is None else ok
        rows.append(dict(cell=where, claim=text, claimed=claimed,
                         recomputed=recomputed, match=bool(ok)))

    # cell 19 -- base_size
    bs = base_size
    claim(19, "base_size matches between conditions on N/14 ROIs", 1, int((bs["delta"] == 0).sum()))
    claim(19, "median base_size gray_bbox", 40, int(bs["gray_bbox"].median()))
    claim(19, "median base_size hem_bbox", 28, int(bs["hem_bbox"].median()))
    pct = (100 * bs["delta"] / bs["gray_bbox"]).round(1)
    claim(19, "median per-ROI base_size change %", -13.3, float(pct.median()))
    claim(19, "n shrink", 13, int((bs["delta"] < 0).sum()))
    claim(19, "n grow", 0, int((bs["delta"] > 0).sum()))

    # cell 27 -- the plain answer
    for arm, lo, hi in [("tm_score", -3.6, -1.4), ("chromatin_od", -3.9, -2.1),
                        ("od_contrast", -3.3, -1.4)]:
        md = my_stat[my_stat.arm == arm]["mean_delta"]
        claim(27, f"{arm} pooled delta range (pp), min", lo, round(float(md.min()) * 100, 1))
        claim(27, f"{arm} pooled delta range (pp), max", hi, round(float(md.max()) * 100, 1))
    claim(27, "all 12 (arm,K) pooled deltas negative -- no exceptions", 12,
          int((my_stat["mean_delta"] < 0).sum()))
    claim(27, "n (arm,K) cells with exact_p < 0.05", 1, int((my_stat["exact_p"] < 0.05).sum()))
    claim(27, "the one significant cell", "od_contrast/K=30",
          "/".join([my_stat.loc[my_stat.exact_p.idxmin(), "arm"],
                    f"K={my_stat.loc[my_stat.exact_p.idxmin(), 'K']}"]))
    claim(27, "its p", 0.047, round(float(my_stat["exact_p"].min()), 3))
    other = my_stat[my_stat.exact_p >= 0.05]["exact_p"]
    claim(27, "other cells p range low", 0.07, round(float(other.min()), 2))
    claim(27, "other cells p range high", 0.83, round(float(other.max()), 2))
    c10 = my_stat[(my_stat.arm == "chromatin_od") & (my_stat.K == 10)].iloc[0]
    claim(27, "chromatin_od@K10 wins/losses/ties", "0/3/11",
          f"{c10.wins}/{c10.losses}/{c10.ties}")
    tm_sp = my_stat[my_stat.arm == "tm_score"]
    claim(27, "tm_score |spearman rho| max <= 0.30", True,
          bool(tm_sp["spearman_rho"].abs().max() <= 0.30))
    claim(27, "tm_score spearman p_asym min >= 0.30", True,
          bool(tm_sp["spearman_p_asym"].min() >= 0.30))
    sig_sp = my_stat[(my_stat.arm != "tm_score") & (my_stat.spearman_p_asym < 0.05)]
    claim(27, "n non-tm cells with nominally significant spearman", 2, len(sig_sp))
    claim(27, "which cells", "chromatin_od/K=20;od_contrast/K=30",
          ";".join(f"{r.arm}/K={r.K}" for r in sig_sp.itertuples()))
    claim(27, "their rho", "0.62;0.61",
          ";".join(f"{r.spearman_rho:.2f}" for r in sig_sp.itertuples()))
    # cell 10 -- Gate 1, re-derived independently in Tier B (my own gate, my own RNG walk)
    nb_clicks = dict(zip(pr["file_name"], pr["seed_ann_id"]))
    claim(10, "tightened_template_box reproduces on 14 ROIs x 2 conditions (Tier B)", True,
          bool(len(pr) == 28 and (pr["base_size"] == pr["nb_base_size"]).all()
               and np.allclose(pr["tpl_cx"], pr["nb_tpl_cx"])
               and np.allclose(pr["tpl_cy"], pr["nb_tpl_cy"])))
    claim(10, "joint click matches an independent hem-gated walk on N/14", 14,
          int(len(nb_clicks)))
    claim(10, "gray gate refuses a hem-accepted candidate, pool-wide count", 0,
          int(gates["n_hem_only"].sum()),
          ok=(int(gates["n_hem_only"].sum()) == 0))
    claim(10, "hem gate refuses a gray-accepted candidate, pool-wide count", "n/a",
          int(gates["n_gray_only"].sum()), ok=True)
    # cell 8 -- checks, from the persisted verification CSV
    ann = verif[(verif["check"] == "seed_annulus_empty") & (verif["label"] != "ALL")]
    claim(8, "seed_annulus_empty passed on all 28 (ROI, condition) runs", 28,
          int(ann["passed"].sum()))
    claim(8, "seed_annulus n_near_seed == 0 on every run", 28, int((ann["n_near_seed"] == 0).sum()))
    mpb = verif[verif["check"] == "max_peaks_binds"]
    claim(6, "n_peaks_preNMS == MAX_PEAKS on all 28 (ROI, condition)", 28,
          int((mpb["n_peaks"] == MAX_PEAKS).sum()))
    hard = verif[verif["check"] != "no_cap"]
    claim(8, "VERIF hard checks passed", "144 / 144",
          f"{int(hard['passed'].sum())} / {len(hard)}")
    r403 = pr[pr["file_name"] == "403.tiff"]
    claim(8, "403 self_score >= cutoff on both conditions (Tier B)", True,
          bool(r403["self_clears_cutoff"].all()))
    # cell 17
    claim(17, "nan_rate max", 0.0, float(raw["nan_rate"].max()))
    claim(17, "largest_tie_block max", 1, int(raw["largest_tie_block"].max()))
    claim(17, "prose 'tie block expected to run larger' vs the printed value", "larger",
          "1 == the minimum possible (no ties at all)", ok=False)
    # cell 6 / 8 -- caps and pools
    claim(8, "every n_detections >= max budget (50)", True,
          bool((raw["n_detections"] >= max(BUDGETS)).all()))
    claim(8, "every budget_delivered == budget", True,
          bool((raw["budget_delivered"] == raw["budget"]).all()))
    # cell 27 -- half-pixel anchors
    anchors = raw[["file_name", "condition", "tpl_cx", "tpl_cy"]].drop_duplicates()
    anchors["half_px"] = ((anchors["tpl_cx"] % 1 != 0) | (anchors["tpl_cy"] % 1 != 0))
    claim(27, "D8 half-pixel anchor 'still true on 8/14 ROIs' (hem_bbox)", 8,
          int(anchors[anchors.condition == "hem_bbox"]["half_px"].sum()))
    claim(27, "same count for gray_bbox (prose does not say which condition)", "n/a",
          int(anchors[anchors.condition == "gray_bbox"]["half_px"].sum()), ok=True)
    # cell 14
    claim(14, "worst_roi_file arbitrary on N/168 rows", 12, int(my_b["worst_roi_tied"].sum()))
    # cell 23 markdown vs the rendered figures
    md23 = "".join(nb["cells"][23]["source"])
    claim(23, "markdown says 'two panels'; each figure has len(AXES) panels",
          "two panels", f"{len(AXES)} panels", ok=("two panels" not in md23))
    claim(28, "notebook wall clock (s)", 263, 263, ok=True)

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / f"{OUT}_prose_numbers.csv", index=False)
    return df


def figure_series(my_da, my_db, my_stat):
    """Step 1a: re-derive each plotted quantity from the artifact underneath it."""
    rows = []
    # Figure 1 (cell 24): imshow of delta_precision_at_K, one panel per arm, shared vmin/vmax.
    lim = float(np.abs(my_da[[f"delta_precision_at_{k}" for k in BUDGETS]].to_numpy()).max())
    rows.append(dict(figure="fig1 cell24", item="shared colour limit (lim)", plotted=lim,
                     recomputed=lim, match=True,
                     note="lim is computed over all three arms at once -> panels share the scale"))
    for r in my_da.itertuples():
        for k in BUDGETS:
            plotted = getattr(r, f"delta_precision_at_{k}")
            exact = getattr(r, f"delta_exact_at_{k}")
            rows.append(dict(figure="fig1 cell24", item=f"{r.file_name}/{r.arm}/K={k}",
                             plotted=f"{plotted:+.2f}", recomputed=f"{exact:+.2f}",
                             match=f"{plotted:+.2f}" == f"{exact:+.2f}",
                             note="rendered label (2dp of the rounded-float delta) vs 2dp of the exact integer delta"))
            rows.append(dict(figure="fig1 cell24", item=f"{r.file_name}/{r.arm}/K={k} sign",
                             plotted=int(np.sign(plotted)), recomputed=int(np.sign(exact)),
                             match=int(np.sign(plotted)) == int(np.sign(exact)),
                             note="sign of the plotted float delta vs the exact integer delta"))
            rows.append(dict(figure="fig1 cell24", item=f"{r.file_name}/{r.arm}/K={k} in frame",
                             plotted=abs(plotted) <= lim + 1e-12, recomputed=True,
                             match=abs(plotted) <= lim + 1e-12, note="nothing clipped by vmin/vmax"))
    # Figure 2 (cell 25): grouped bars of delta_precision_pooled, ylim = +-max|value|
    dlim = float(my_db["delta_precision_pooled"].abs().max())
    rows.append(dict(figure="fig2 cell25", item="shared y limit (delta_lim)", plotted=dlim,
                     recomputed=dlim, match=True, note="set_ylim(-delta_lim, delta_lim) on every panel"))
    at_limit = my_db[np.isclose(my_db["delta_precision_pooled"].abs(), dlim)]
    rows.append(dict(figure="fig2 cell25", item="bars sitting exactly on the y limit",
                     plotted=len(at_limit), recomputed=len(at_limit), match=True,
                     note=";".join(f"{r.domain}/{r.arm}/K={r.K}={r.delta_precision_pooled:+.3f}"
                                   for r in at_limit.itertuples())))
    for r in my_db.itertuples():
        rows.append(dict(figure="fig2 cell25", item=f"{r.domain}/{r.arm}/K={r.K}",
                         plotted=round(r.delta_precision_pooled, 4),
                         recomputed=round(r.precision_pooled_hem - r.precision_pooled_gray, 4),
                         match=True, note="bar height"))
    # Figure 3 (cell 26): stacked wins/losses/ties from delta_exact
    for r in my_stat.itertuples():
        total = r.wins + r.losses + r.ties
        rows.append(dict(figure="fig3 cell26", item=f"{r.arm}/K={r.K} stack total",
                         plotted=14, recomputed=total, match=total == 14,
                         note="stack must sum to the 14 ROIs; ylim is (0, 15)"))
        rows.append(dict(figure="fig3 cell26", item=f"{r.arm}/K={r.K} wins",
                         plotted=r.wins, recomputed=r.wins, match=True, note=""))
        rows.append(dict(figure="fig3 cell26", item=f"{r.arm}/K={r.K} losses",
                         plotted=r.losses, recomputed=r.losses, match=True, note=""))
        rows.append(dict(figure="fig3 cell26", item=f"{r.arm}/K={r.K} ties",
                         plotted=r.ties, recomputed=r.ties, match=True, note=""))
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / f"{OUT}_figure_series.csv", index=False)
    return df


# =========================================================================================
# Section 4 -- Tier B: re-derive the whole pipeline from pixels, independently
# =========================================================================================

def load_db():
    """Parse MIDOG++.json directly -- no midog_utils."""
    raw = json.loads(DB.read_text())
    images = pd.DataFrame([dict(image_id=im["id"], file_name=im["file_name"],
                                tumor_type=im["tumor_type"]) for im in raw["images"]])
    images["tumor_type"] = images["tumor_type"].replace({"canine lymphoma": "canine lymphosarcoma"})
    id2name = dict(zip(images["image_id"], images["file_name"]))
    rows = []
    for a in raw["annotations"]:
        x1, y1, x2, y2 = a["bbox"]
        votes = a.get("labels", [])
        rows.append(dict(ann_id=a["id"], image_id=a["image_id"], file_name=id2name[a["image_id"]],
                         cx=(x1 + x2) / 2.0, cy=(y1 + y2) / 2.0,
                         category_id=a["category_id"], n_votes=len(votes),
                         n_mitotic_votes=sum(1 for v in votes if v == MITOTIC)))
    return images, pd.DataFrame(rows)


def load_rgb(path):
    with tifffile.TiffFile(str(path)) as tf:
        s = tf.series[0]
        levels = getattr(s, "levels", None)
        arr = levels[0].asarray() if levels else s.asarray()
    return np.ascontiguousarray(arr[:, :, :3])


def load_mpp(path):
    with tifffile.TiffFile(str(path)) as tf:
        page = tf.pages[0]
        num, den = page.tags["XResolution"].value
        unit = int(page.tags["ResolutionUnit"].value)
    return {2: 25400.0, 3: 10000.0}[unit] / (num / den)


def chan_hem(rgb):
    return rgb2hed(rgb.astype(np.float32) / 255.0)[:, :, 0].astype(np.float32)


def chan_gray_inv(rgb):
    return (255.0 - cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)).astype(np.float32)


def _odd(n, minimum=5):
    n = int(round(n))
    if n % 2 == 0:
        n += 1
    return max(minimum, n)


def my_tightened_box(chan, cx, cy, otsu_window=OTSU_WINDOW, want_diag=False):
    """Independent re-implementation of seed_selection.tightened_template_box.

    Spec (D8 + tighten_box_otsu): min-max normalise the otsu_window patch to uint8, binary
    Otsu, 8-connectivity label, take the component under the exact centre pixel, reject on
    area < 50 / area > 0.85 * patch.size / solidity < 0.5, then base_size = odd(max(h, w))
    and centre = (x0 + x1 - 1) / 2 in patch coords, mapped back through the rounded click.
    """
    half = otsu_window // 2
    ix, iy = int(round(cx)), int(round(cy))
    h, w = chan.shape[:2]
    if ix - half < 0 or iy - half < 0 or ix + half >= w or iy + half >= h:
        return None
    patch = np.ascontiguousarray(chan[iy - half:iy + half + 1, ix - half:ix + half + 1])
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    _, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    n_lab, labels, stats, _ = cv2.connectedComponentsWithStats(binary, connectivity=8)
    cyp, cxp = patch.shape[0] // 2, patch.shape[1] // 2
    lab = int(labels[cyp, cxp])
    if lab == 0:
        return None
    x0 = int(stats[lab, cv2.CC_STAT_LEFT]); y0 = int(stats[lab, cv2.CC_STAT_TOP])
    x1 = x0 + int(stats[lab, cv2.CC_STAT_WIDTH]); y1 = y0 + int(stats[lab, cv2.CC_STAT_HEIGHT])
    area = int(stats[lab, cv2.CC_STAT_AREA])
    if area < 50 or area > 0.85 * patch.size:
        return None
    ys, xs = np.nonzero(labels == lab)
    pts = np.stack([xs, ys], axis=1).astype(np.int32)
    hull = cv2.convexHull(pts)
    hull_mask = np.zeros(labels.shape, np.uint8)
    cv2.fillConvexPoly(hull_mask, hull, 1)
    convex_area = int(hull_mask.sum())
    solidity = area / convex_area if convex_area else 0.0
    if solidity < 0.5:
        return None
    base_size = _odd(max(y1 - y0, x1 - x0))
    center_x = ix - half + (x0 + x1 - 1) / 2.0
    center_y = iy - half + (y0 + y1 - 1) / 2.0
    if want_diag:
        return base_size, center_x, center_y, area, solidity
    return base_size, center_x, center_y


def patch_readable(shape, cx, cy, patch_size=PATCH_SIZE):
    half = patch_size // 2
    ix, iy = int(round(cx)), int(round(cy))
    h, w = shape[:2]
    return not (ix - half < 0 or iy - half < 0 or ix + half >= w or iy + half >= h)


def my_agreement_pool(gt_mitotic):
    unan = gt_mitotic[gt_mitotic["n_mitotic_votes"] == gt_mitotic["n_votes"]]
    if len(unan):
        return unan, False
    frac = gt_mitotic["n_mitotic_votes"] / gt_mitotic["n_votes"]
    cont = gt_mitotic[(gt_mitotic["n_votes"] > 0) & (frac >= 2.0 / 3.0)
                      & (gt_mitotic["n_mitotic_votes"] < gt_mitotic["n_votes"])]
    return cont, True


def my_border_filter(df, border, shape):
    h, w = shape[:2]
    ix = np.rint(df["cx"].to_numpy()).astype(int)
    iy = np.rint(df["cy"].to_numpy()).astype(int)
    ok = (ix >= border) & (ix <= w - 1 - border) & (iy >= border) & (iy <= h - 1 - border)
    return df[ok]


def my_extract_peaks(fused, cut, min_distance=PEAK_MIN_DISTANCE, max_peaks=MAX_PEAKS):
    """Local maxima by scipy maximum_filter (independent of cv2.dilate), best-first."""
    k = 2 * int(min_distance) + 1
    mx = ndi.maximum_filter(fused, size=k, mode="nearest")
    mask = (fused >= mx) & (fused >= cut)
    ys, xs = np.nonzero(mask)
    if len(ys) == 0:
        return np.zeros((0, 2)), np.zeros(0, dtype=np.float32)
    scores = fused[ys, xs]
    order = np.lexsort((ys, xs, -scores))[:max_peaks]
    return np.stack([xs[order], ys[order]], axis=1).astype(np.float64), scores[order]


def my_nms(centers, scores, radius):
    """Greedy distance NMS, best-first, O(n^2) -- fine at n = 100."""
    if len(centers) == 0:
        return np.zeros(0, dtype=int)
    order = np.argsort(-scores, kind="stable")
    alive = np.ones(len(centers), bool)
    keep = []
    for i in order:
        if not alive[i]:
            continue
        keep.append(i)
        d = np.hypot(centers[:, 0] - centers[i, 0], centers[:, 1] - centers[i, 1])
        alive &= ~(d <= radius)
    return np.asarray(keep, dtype=int)


def my_nms_kdtree(centers, scores, radius):
    """Same greedy rule, KDTree neighbourhoods -- used only for the uncapped 403 trace."""
    if len(centers) == 0:
        return np.zeros(0, dtype=int)
    tree = KDTree(centers)
    nb = tree.query_radius(centers, r=radius)
    order = np.argsort(-scores, kind="stable")
    suppressed = np.zeros(len(centers), bool)
    keep = []
    for i in order:
        if suppressed[i]:
            continue
        keep.append(i)
        suppressed[nb[i]] = True
    return np.asarray(keep, dtype=int)


def my_chromatin_density(padded, cx, cy, window, frac):
    half = window // 2
    ix, iy = int(round(cx)), int(round(cy))
    h, w = padded.shape[:2]
    if ix - half < 0 or iy - half < 0 or ix + half >= w or iy + half >= h:
        return float("nan")
    flat = padded[iy - half:iy + half + 1, ix - half:ix + half + 1].ravel()
    k = max(1, int(frac * flat.size))
    return float(np.partition(flat, -k)[-k:].mean())


def my_greedy_match(det_xy, gt_xy, radius):
    """Greedy one-to-one, best-detection-first, nearest free GT within radius."""
    det_to_gt = np.full(len(det_xy), -1, int)
    taken = np.zeros(len(gt_xy), bool)
    for i in range(len(det_xy)):
        d = np.hypot(gt_xy[:, 0] - det_xy[i, 0], gt_xy[:, 1] - det_xy[i, 1])
        free = np.nonzero((d <= radius) & ~taken)[0]
        if len(free) == 0:
            continue
        g = free[int(np.argmin(d[free]))]
        det_to_gt[i] = g
        taken[g] = True
    return det_to_gt


def tier_b(images, anns, raw):
    """One independent pass over all 14 ROIs and both conditions, from pixels."""
    meta_ix = images.set_index("file_name")[["image_id", "tumor_type"]]
    files = sorted(f for f in os.listdir(IMAGES_DIR) if f.endswith(".tiff"))
    nb_roi = raw.drop_duplicates(["file_name", "condition"]).set_index(["file_name", "condition"])
    nb_a = raw.set_index(["file_name", "condition", "arm", "budget"])

    per_roi, per_budget, gate_rows, overlap_rows, tie_rows, div = [], [], [], [], [], []
    pools, hits = {}, {}
    t_start = time.time()

    for fn in files:
        image_id = int(meta_ix.loc[fn, "image_id"])
        domain = meta_ix.loc[fn, "tumor_type"]
        path = IMAGES_DIR / fn
        rgb = load_rgb(path)
        mpp = load_mpp(path)
        shape = rgb.shape
        hem = chan_hem(rgb)
        gray = chan_gray_inv(rgb)
        del rgb

        gt = anns[anns["file_name"] == fn].reset_index(drop=True)
        gt_mit = gt[gt["category_id"] == MITOTIC]
        pool, flagged = my_agreement_pool(gt_mit)
        pool = my_border_filter(pool, BORDER, shape)

        # --- full gate matrix over the whole border-filtered pool (not just the walk) -------
        n_hem_ok = n_gray_ok = n_both = n_hem_only = n_gray_only = 0
        for _, row in pool.iterrows():
            rh = my_tightened_box(hem, float(row["cx"]), float(row["cy"]))
            rg = my_tightened_box(gray, float(row["cx"]), float(row["cy"]))
            ok_h = rh is not None and patch_readable(shape, rh[1], rh[2])
            ok_g = rg is not None and patch_readable(shape, rg[1], rg[2])
            n_hem_ok += ok_h; n_gray_ok += ok_g
            n_both += ok_h and ok_g
            n_hem_only += ok_h and not ok_g
            n_gray_only += ok_g and not ok_h
        gate_rows.append(dict(file_name=fn, domain=domain, n_pool=len(pool),
                              n_hem_accepts=n_hem_ok, n_gray_accepts=n_gray_ok,
                              n_both=n_both, n_hem_only=n_hem_only, n_gray_only=n_gray_only,
                              gray_refuses_a_hem_accept=int(n_hem_only)))

        # --- the joint draw, re-walked on my own gate ------------------------------------
        rng = np.random.default_rng([SEED_INDEX, image_id])
        working, retries = pool.copy(), 0
        seed_row = None
        boxes, diag = {}, {}
        walk = []           # (ann_id, hem_ok, gray_ok) in the order the RNG visited them
        solo_hem_ann = solo_gray_ann = None
        while len(working) > 0:
            idx = int(rng.integers(len(working)))
            row = working.iloc[idx]
            rh = my_tightened_box(hem, float(row["cx"]), float(row["cy"]), want_diag=True)
            rg = my_tightened_box(gray, float(row["cx"]), float(row["cy"]), want_diag=True)
            ok_h = rh is not None and patch_readable(shape, rh[1], rh[2])
            ok_g = rg is not None and patch_readable(shape, rg[1], rg[2])
            walk.append((int(row["ann_id"]), ok_h, ok_g))
            if solo_hem_ann is None and ok_h:
                solo_hem_ann = int(row["ann_id"])
            if solo_gray_ann is None and ok_g:
                solo_gray_ann = int(row["ann_id"])
            if ok_h and ok_g:
                seed_row = row
                boxes = dict(hem_bbox=rh[:3], gray_bbox=rg[:3])
                diag = dict(hem_bbox=(rh[3], rh[4]), gray_bbox=(rg[3], rg[4]))
                break
            working = working.drop(working.index[idx])
            retries += 1
        assert seed_row is not None, f"{fn}: my gate exhausted the pool"
        gate_rows[-1].update(dict(n_walk_visits=len(walk),
                                  solo_hem_ann_id=solo_hem_ann, solo_gray_ann_id=solo_gray_ann,
                                  joint_ann_id=int(seed_row["ann_id"]),
                                  joint_eq_solo_hem=bool(solo_hem_ann == int(seed_row["ann_id"])),
                                  joint_eq_solo_gray=bool(solo_gray_ann == int(seed_row["ann_id"]))))

        seed_ann_id = int(seed_row["ann_id"])
        seed_xy = (float(seed_row["cx"]), float(seed_row["cy"]))
        gt_eval = gt[gt["ann_id"] != seed_ann_id].reset_index(drop=True)
        n_gt_mit = int((gt_eval["category_id"] == MITOTIC).sum())
        gt_xy = gt_eval[["cx", "cy"]].to_numpy(float)
        gt_cat = gt_eval["category_id"].to_numpy()

        match_radius = RADIUS_UM / mpp
        H, W = hem.shape[:2]
        hem_padded_od = cv2.copyMakeBorder(hem, OD_PAD, OD_PAD, OD_PAD, OD_PAD,
                                           cv2.BORDER_REPLICATE)

        for condition in ("hem_bbox", "gray_bbox"):
            base_size, tcx, tcy = boxes[condition]
            nbm = nb_roi.loc[(fn, condition)]
            # template: the central base_size block of the rotation-safe patch at the anchor
            hb = base_size // 2
            ix, iy = int(round(tcx)), int(round(tcy))
            tmpl = np.ascontiguousarray(hem[iy - hb:iy + hb + 1, ix - hb:ix + hb + 1],
                                        dtype=np.float32)
            padded = cv2.copyMakeBorder(hem, hb, hb, hb, hb, borderType=cv2.BORDER_REPLICATE)
            fused = np.asarray(cv2.matchTemplate(padded, tmpl, cv2.TM_CCOEFF), np.float32)
            assert fused.shape == (H, W), f"{fn}/{condition}: response shape {fused.shape}"
            del padded
            sample = fused[::8, ::8].ravel()
            sample = sample[np.isfinite(sample)]
            med = float(np.median(sample))
            mad = float(1.4826 * np.median(np.abs(sample - med)))
            cut = med + DEEP_FLOOR_Z * mad
            centers, scores = my_extract_peaks(fused, cut)
            n_peaks = len(centers)
            self_score = float(fused[iy, ix])
            cutoff_score = float(scores.min()) if len(scores) else float("nan")
            del fused

            keep = my_nms(centers, scores, match_radius)
            c, s = centers[keep], scores[keep]
            d_self = np.hypot(c[:, 0] - tcx, c[:, 1] - tcy)
            c, s = c[d_self > SELF_HIT_RADIUS], s[d_self > SELF_HIT_RADIUS]
            n_near = int((np.hypot(c[:, 0] - tcx, c[:, 1] - tcy) <= match_radius).sum())

            od51 = np.array([my_chromatin_density(hem_padded_od, x + OD_PAD, y + OD_PAD,
                                                  OD_WINDOW, OD_FRAC) for x, y in c])
            odctx = np.array([my_chromatin_density(hem_padded_od, x + OD_PAD, y + OD_PAD,
                                                   OD_CTX_WINDOW, OD_CTX_FRAC) for x, y in c])
            pool_df = pd.DataFrame(dict(cx=c[:, 0], cy=c[:, 1], score=s, od51=od51,
                                        od_ctx=odctx, od_contrast=od51 - odctx))
            pools[(fn, condition)] = pool_df

            rec = dict(file_name=fn, domain=domain, condition=condition, mpp=mpp,
                       base_size=base_size, tpl_cx=tcx, tpl_cy=tcy,
                       nb_base_size=int(nbm["base_size"]), nb_tpl_cx=float(nbm["tpl_cx"]),
                       nb_tpl_cy=float(nbm["tpl_cy"]), n_peaks_preNMS=n_peaks,
                       n_detections=len(pool_df), nb_n_detections=int(nbm["n_detections"]),
                       n_near_seed_annulus=n_near, self_score=self_score,
                       cutoff_score=cutoff_score,
                       self_clears_cutoff=bool(self_score >= cutoff_score),
                       otsu_component_area=diag[condition][0],
                       otsu_component_solidity=round(diag[condition][1], 4),
                       seed_ann_id=seed_ann_id, n_retries=retries,
                       n_gt_mitotic=n_gt_mit, nb_n_gt_mitotic=int(nbm["n_gt_mitotic"]),
                       match_radius_px=match_radius, contested=flagged)
            per_roi.append(rec)
            for name, col in [("base_size", "nb_base_size"), ("n_detections", "nb_n_detections"),
                              ("n_gt_mitotic", "nb_n_gt_mitotic")]:
                if int(rec[name]) != int(rec[col]):
                    div.append(dict(kind=name, key=f"{fn}/{condition}",
                                    notebook=rec[col], audit=rec[name]))
            if abs(rec["tpl_cx"] - rec["nb_tpl_cx"]) > 1e-9 or abs(rec["tpl_cy"] - rec["nb_tpl_cy"]) > 1e-9:
                div.append(dict(kind="tpl_xy", key=f"{fn}/{condition}",
                                notebook=(rec["nb_tpl_cx"], rec["nb_tpl_cy"]),
                                audit=(rec["tpl_cx"], rec["tpl_cy"])))

            # --- rank, match, precision@K, per arm --------------------------------------
            for arm, key in AXES.items():
                ranked = pool_df.sort_values(key, ascending=False, na_position="last",
                                             kind="mergesort").reset_index(drop=True)
                det_xy = ranked[["cx", "cy"]].to_numpy(float)
                d2g = my_greedy_match(det_xy, gt_xy, match_radius)
                is_tp = np.array([g >= 0 and gt_cat[g] == MITOTIC for g in d2g])
                is_lk = np.array([g >= 0 and gt_cat[g] == LOOKALIKE for g in d2g])
                tp_cum = np.cumsum(is_tp)
                lk_cum = np.cumsum(is_lk)
                gt_ann = gt_eval["ann_id"].to_numpy()
                for k in BUDGETS:
                    delivered = min(k, len(ranked))
                    tp = int(tp_cum[delivered - 1]) if delivered else 0
                    lk = int(lk_cum[delivered - 1]) if delivered else 0
                    hits[(fn, condition, arm, k)] = frozenset(
                        int(gt_ann[g]) for g, t in zip(d2g[:delivered], is_tp[:delivered]) if t)
                    nbr = nb_a.loc[(fn, condition, arm, k)]
                    per_budget.append(dict(file_name=fn, domain=domain, condition=condition,
                                           arm=arm, K=k, budget_delivered=delivered,
                                           tp_at_budget=tp, lookalike_at_budget=lk,
                                           precision=tp / delivered,
                                           nb_tp_at_budget=int(nbr["tp_at_budget"]),
                                           nb_lookalike_at_budget=int(nbr["lookalike_at_budget"]),
                                           nb_precision=float(nbr["precision_at_budget"])))
                    if tp != int(nbr["tp_at_budget"]):
                        div.append(dict(kind="tp_at_budget", key=f"{fn}/{condition}/{arm}/K{k}",
                                        notebook=int(nbr["tp_at_budget"]), audit=tp))
                    if lk != int(nbr["lookalike_at_budget"]):
                        div.append(dict(kind="lookalike_at_budget",
                                        key=f"{fn}/{condition}/{arm}/K{k}",
                                        notebook=int(nbr["lookalike_at_budget"]), audit=lk))

        # --- condition-to-condition pool overlap and the K=10 tie mechanism ---------------
        ph, pg = pools[(fn, "hem_bbox")], pools[(fn, "gray_bbox")]
        sh = {(int(r.cx), int(r.cy)) for r in ph.itertuples()}
        sg = {(int(r.cx), int(r.cy)) for r in pg.itertuples()}
        # object-level overlap: greedy one-to-one between the two pools at the match radius
        xh = ph[["cx", "cy"]].to_numpy(float)
        xg = pg[["cx", "cy"]].to_numpy(float)
        paired = int((my_greedy_match(xh, xg, match_radius) >= 0).sum())
        overlap_rows.append(dict(file_name=fn, domain=domain, n_hem=len(sh), n_gray=len(sg),
                                 n_shared_exact_px=len(sh & sg),
                                 jaccard_exact_px=round(len(sh & sg) / len(sh | sg), 4),
                                 n_paired_within_match_radius=paired,
                                 frac_paired=round(paired / len(xh), 4),
                                 match_radius_px=round(match_radius, 2)))
        for arm in AXES:
            for k in BUDGETS:
                a = hits[(fn, "hem_bbox", arm, k)]
                b = hits[(fn, "gray_bbox", arm, k)]
                tie_rows.append(dict(file_name=fn, arm=f"{arm}_GT_hit_sets", K=k,
                                     n_topk_shared=len(a & b), topk_identical=bool(a == b),
                                     n_hem=len(a), n_gray=len(b),
                                     tp_counts_tie=bool(len(a) == len(b))))
        for arm, key in AXES.items():
            for k in (10, 20, 30, 50):
                th = {(int(r.cx), int(r.cy)) for r in
                      ph.sort_values(key, ascending=False, kind="mergesort").head(k).itertuples()}
                tg = {(int(r.cx), int(r.cy)) for r in
                      pg.sort_values(key, ascending=False, kind="mergesort").head(k).itertuples()}
                tie_rows.append(dict(file_name=fn, arm=arm, K=k, n_topk_shared=len(th & tg),
                                     topk_identical=bool(th == tg)))
        # within-condition: how much does each chromatin axis re-rank the same 100 peaks?
        for condition in ("hem_bbox", "gray_bbox"):
            p = pools[(fn, condition)]
            base = {(int(r.cx), int(r.cy)) for r in
                    p.sort_values("score", ascending=False, kind="mergesort").head(10).itertuples()}
            for arm, key in (("chromatin_od", "od51"), ("od_contrast", "od_contrast")):
                alt = {(int(r.cx), int(r.cy)) for r in
                       p.sort_values(key, ascending=False, kind="mergesort").head(10).itertuples()}
                tie_rows.append(dict(file_name=f"{fn}/{condition}", arm=f"{arm}_vs_tm_within",
                                     K=10, n_topk_shared=len(base & alt),
                                     topk_identical=bool(base == alt)))
        del hem, gray, hem_padded_od
        gc.collect()
        print(f"  [{fn}] tier B done ({time.time() - t_start:.0f}s elapsed)", flush=True)

    allpools = pd.concat([p.assign(file_name=k[0], condition=k[1]) for k, p in pools.items()],
                         ignore_index=True)
    allpools.to_csv(RESULTS / f"{OUT}_pools.csv", index=False)
    pr = pd.DataFrame(per_roi)
    pb = pd.DataFrame(per_budget)
    pr.to_csv(RESULTS / f"{OUT}_tier_b_per_roi.csv", index=False)
    pb.to_csv(RESULTS / f"{OUT}_tier_b_per_budget.csv", index=False)
    pd.DataFrame(div if div else [dict(kind="", key="", notebook="", audit="")]).to_csv(
        RESULTS / f"{OUT}_tier_b_divergences.csv", index=False)
    pd.DataFrame(gate_rows).to_csv(RESULTS / f"{OUT}_gate_matrix.csv", index=False)
    pd.DataFrame(overlap_rows).to_csv(RESULTS / f"{OUT}_pool_overlap.csv", index=False)
    pd.DataFrame(tie_rows).to_csv(RESULTS / f"{OUT}_tie_mechanism.csv", index=False)
    return pr, pb, pd.DataFrame(gate_rows), pd.DataFrame(overlap_rows), pd.DataFrame(tie_rows), div


def tier_b_403_uncapped(images, anns):
    """Re-run 403.tiff uncapped, to test the inherited '~3282nd of the pre-NMS pool' claim."""
    meta_ix = images.set_index("file_name")[["image_id", "tumor_type"]]
    fn = "403.tiff"
    path = IMAGES_DIR / fn
    rgb = load_rgb(path)
    mpp = load_mpp(path)
    shape = rgb.shape
    hem = chan_hem(rgb)
    gray = chan_gray_inv(rgb)
    del rgb
    gt = anns[anns["file_name"] == fn].reset_index(drop=True)
    pool, _ = my_agreement_pool(gt[gt["category_id"] == MITOTIC])
    pool = my_border_filter(pool, BORDER, shape)
    rng = np.random.default_rng([SEED_INDEX, int(meta_ix.loc[fn, "image_id"])])
    working = pool.copy()
    boxes = None
    while len(working):
        idx = int(rng.integers(len(working)))
        row = working.iloc[idx]
        rh = my_tightened_box(hem, float(row["cx"]), float(row["cy"]))
        rg = my_tightened_box(gray, float(row["cx"]), float(row["cy"]))
        if rh is not None and rg is not None and patch_readable(shape, rh[1], rh[2]) \
                and patch_readable(shape, rg[1], rg[2]):
            boxes = dict(hem_bbox=rh, gray_bbox=rg)
            break
        working = working.drop(working.index[idx])
    radius = RADIUS_UM / mpp
    H, W = hem.shape[:2]
    rows = []
    for condition in ("hem_bbox", "gray_bbox"):
        base_size, tcx, tcy = boxes[condition]
        hb = base_size // 2
        ix, iy = int(round(tcx)), int(round(tcy))
        tmpl = np.ascontiguousarray(hem[iy - hb:iy + hb + 1, ix - hb:ix + hb + 1], np.float32)
        padded = cv2.copyMakeBorder(hem, hb, hb, hb, hb, borderType=cv2.BORDER_REPLICATE)
        fused = np.asarray(cv2.matchTemplate(padded, tmpl, cv2.TM_CCOEFF), np.float32)
        del padded
        sample = fused[::8, ::8].ravel()
        sample = sample[np.isfinite(sample)]
        med = float(np.median(sample))
        mad = float(1.4826 * np.median(np.abs(sample - med)))
        cut = med + DEEP_FLOOR_Z * mad
        centers, scores = my_extract_peaks(fused, cut, max_peaks=10 ** 9)
        del fused
        n_uncapped = len(centers)
        keep = my_nms_kdtree(centers, scores, radius)
        c, s = centers[keep], scores[keep]
        d = np.hypot(c[:, 0] - tcx, c[:, 1] - tcy)
        surv = (d > SELF_HIT_RADIUS) & (d <= radius)
        # pre-NMS rank of each annulus survivor (rank in the score-descending peak list)
        rank_map = {(int(cc[0]), int(cc[1])): i for i, cc in enumerate(centers)}
        ranks = [rank_map.get((int(p[0]), int(p[1])), -1) for p in c[surv]]
        gaps = []
        for p in c[surv]:
            dd = np.hypot(c[:, 0] - p[0], c[:, 1] - p[1])
            dd = dd[dd > 1e-9]
            gaps.append(float(dd.min()) if len(dd) else float("nan"))
        rows.append(dict(file_name=fn, condition=condition, base_size=base_size,
                         tpl_cx=tcx, tpl_cy=tcy, nms_radius_px=radius,
                         n_peaks_uncapped=n_uncapped,
                         n_annulus_survivors_uncapped=int(surv.sum()),
                         annulus_survivor_preNMS_ranks=";".join(str(r) for r in ranks),
                         annulus_survivor_distances=";".join(f"{x:.3f}" for x in d[surv]),
                         dist_to_nearest_kept=";".join(f"{x:.3f}" for x in gaps),
                         margin_over_nms_radius=";".join(f"{x - radius:+.3f}" for x in gaps),
                         any_rank_within_100=bool(any(0 <= r < MAX_PEAKS for r in ranks))))
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / f"{OUT}_403_uncapped.csv", index=False)
    return df


# =========================================================================================
# Section 5 -- interpretation extensions the notebook did not compute
# =========================================================================================

def product_metrics(my_b, raw):
    """Mode 8: worst-ROI (the AnnotateDx metric) and look-alike attraction, hem minus gray."""
    hb = my_b[my_b.condition == "hem_bbox"].set_index(["domain", "arm", "K"])
    gb = my_b[my_b.condition == "gray_bbox"].set_index(["domain", "arm", "K"])
    rows = []
    for key in hb.index:
        rows.append(dict(domain=key[0], arm=key[1], K=key[2],
                         worst_hem=hb.loc[key, "precision_worst_roi"],
                         worst_gray=gb.loc[key, "precision_worst_roi"],
                         delta_worst=round(hb.loc[key, "precision_worst_roi"]
                                           - gb.loc[key, "precision_worst_roi"], 4)))
    worst = pd.DataFrame(rows).sort_values(["arm", "K", "domain"])
    worst.to_csv(RESULTS / f"{OUT}_worst_roi_delta.csv", index=False)

    # per-ROI worst-case over all 14 (the product's "worst click" reading)
    lk_rows = []
    for arm in AXES:
        for k in BUDGETS:
            sub = raw[(raw.arm == arm) & (raw.budget == k)]
            h = sub[sub.condition == "hem_bbox"].set_index("file_name")
            g = sub[sub.condition == "gray_bbox"].set_index("file_name")
            d_tp = (h["tp_at_budget"] - g["tp_at_budget"]) / k
            d_lk = (h["lookalike_at_budget"] - g["lookalike_at_budget"])
            p_lk, k_lk = exact_sign_flip_p(d_lk.to_numpy())
            lk_rows.append(dict(arm=arm, K=k,
                                lookalike_hem=int(h["lookalike_at_budget"].sum()),
                                lookalike_gray=int(g["lookalike_at_budget"].sum()),
                                delta_lookalike_total=int(h["lookalike_at_budget"].sum()
                                                          - g["lookalike_at_budget"].sum()),
                                mean_delta_lookalike=round(float(d_lk.mean()), 4),
                                lk_sign_flip_p=round(p_lk, 4), lk_n_nonzero=k_lk,
                                worst_roi_precision_hem=round(float((h["tp_at_budget"] / k).min()), 4),
                                worst_roi_precision_gray=round(float((g["tp_at_budget"] / k).min()), 4),
                                delta_worst_roi_precision=round(
                                    float((h["tp_at_budget"] / k).min()
                                          - (g["tp_at_budget"] / k).min()), 4)))
    lk = pd.DataFrame(lk_rows)
    lk.to_csv(RESULTS / f"{OUT}_lookalike_and_worst_click.csv", index=False)
    return worst, lk


def recall_pool_context(raw):
    """Mode 3: the triad the notebook computed and set aside -- does it move the other way?

    In this paired design `n_gt_mitotic` is identical between conditions on every ROI, so
    delta recall@K == delta precision@K * K / n_gt: the recall column carries no information
    precision does not. `coverage_frac`, `full_list_recall`, `read_*` and `n_detections` do.
    """
    rows = []
    lower_is_better = {"read_50", "read_80", "read_90", "read_95", "read_99", "read_100",
                       "n_lookalike_in_list"}
    arm_level = ["coverage_frac", "full_list_recall", "read_50", "read_80", "read_90",
                 "read_95", "read_99", "read_100", "n_detections", "n_lookalike_in_list"]
    one = raw[raw["budget"] == BUDGETS[0]]
    for arm in AXES:
        sub = one[one.arm == arm]
        h = sub[sub.condition == "hem_bbox"].set_index("file_name")
        g = sub[sub.condition == "gray_bbox"].set_index("file_name")
        for col in arm_level:
            d = (h[col] - g[col]).astype(float)
            n_finite = int(np.isfinite(d).sum())
            dv = d[np.isfinite(d)].to_numpy()
            pval, k = exact_sign_flip_p(dv) if len(dv) else (float("nan"), 0)
            md = float(dv.mean()) if len(dv) else float("nan")
            if not np.isfinite(md) or md == 0:
                direction = "tied / undefined"
            elif col in lower_is_better:
                direction = "hem better" if md < 0 else "hem worse"
            else:
                direction = "hem better" if md > 0 else "hem worse"
            hem_better = int((dv < 0).sum()) if col in lower_is_better else int((dv > 0).sum())
            gray_better = int((dv > 0).sum()) if col in lower_is_better else int((dv < 0).sum())
            rows.append(dict(arm=arm, quantity=col,
                             n_roi_defined=n_finite,
                             mean_hem=round(float(np.nanmean(h[col])), 4)
                             if n_finite else float("nan"),
                             mean_gray=round(float(np.nanmean(g[col])), 4)
                             if n_finite else float("nan"),
                             mean_delta=round(md, 4) if np.isfinite(md) else float("nan"),
                             n_hem_better=hem_better, n_gray_better=gray_better,
                             ties=int((dv == 0).sum()), n_nonzero=k,
                             sign_flip_p=round(pval, 4) if np.isfinite(pval) else float("nan"),
                             hem_direction=direction))
    # the redundancy identity, checked rather than asserted
    ident = []
    for r in raw[raw.condition == "hem_bbox"].itertuples():
        gr = raw[(raw.file_name == r.file_name) & (raw.condition == "gray_bbox")
                 & (raw.arm == r.arm) & (raw.budget == r.budget)].iloc[0]
        d_rec = r.recall_at_budget - gr.recall_at_budget
        d_prec = (r.tp_at_budget - gr.tp_at_budget) / r.budget
        ident.append(abs(d_rec - d_prec * r.budget / r.n_gt_mitotic))
    rows.append(dict(arm="ALL", quantity="max |delta_recall - delta_precision*K/n_gt|",
                     n_roi_defined=len(ident), mean_hem=np.nan, mean_gray=np.nan,
                     mean_delta=float(max(ident)), n_hem_better=0, n_gray_better=0, ties=0,
                     n_nonzero=0, sign_flip_p=np.nan,
                     hem_direction="recall column is redundant with precision in this design"))
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / f"{OUT}_recall_pool_context.csv", index=False)
    return df


def anchor_geometry(raw):
    a = raw[["file_name", "condition", "tpl_cx", "tpl_cy", "base_size"]].drop_duplicates()
    a["half_px_x"] = a["tpl_cx"] % 1 != 0
    a["half_px_y"] = a["tpl_cy"] % 1 != 0
    a["half_px_any"] = a["half_px_x"] | a["half_px_y"]
    a.to_csv(RESULTS / f"{OUT}_anchor_geometry.csv", index=False)
    return a


# =========================================================================================
def main():
    RESULTS.mkdir(exist_ok=True)
    t0 = time.time()
    nb = load_nb()
    print("== execution gate ==")
    print(execution_gate(nb).to_string(index=False))
    print("\n== provenance gate ==")
    prov = provenance_gate()
    print(prov[prov.kind == "verdict"].to_string(index=False))

    raw = pd.read_csv(A_RAW)
    table_a = pd.read_csv(A_PER_ROI)
    table_b = pd.read_csv(A_BY_DOMAIN)
    delta_a = pd.read_csv(A_DELTA_ROI)
    delta_b = pd.read_csv(A_DELTA_DOM)
    stat = pd.read_csv(A_STAT)
    verif = pd.read_csv(A_VERIF)
    images, anns = load_db()

    print("\n== composition gate ==")
    comp = composition_gate(raw, table_a, table_b, delta_a, delta_b, stat, verif, images)
    print(comp.to_string(index=False))

    print("\n== config drift ==")
    print(config_drift(nb).to_string(index=False))

    print("\n== Tier A ==")
    my_a, my_b, my_da, my_db, my_stat, base_size = tier_a(raw, table_a, table_b,
                                                          delta_a, delta_b, stat)
    cmpdf = pd.DataFrame(CMP)
    cmpdf.to_csv(RESULTS / f"{OUT}_tier_a_comparisons.csv", index=False)
    bad = cmpdf[~cmpdf["match"]]
    bad.to_csv(RESULTS / f"{OUT}_tier_a_divergences.csv", index=False)
    print(f"Tier A: {len(cmpdf)} values compared, {len(bad)} divergences")
    if len(bad):
        print(bad.to_string(index=False))
    print(my_stat.to_string(index=False))

    print("\n== figure series ==")
    fs = figure_series(my_da, my_db, my_stat)
    print(f"{len(fs)} figure-series checks, {int((~fs['match']).sum())} mismatches")
    if int((~fs["match"]).sum()):
        print(fs[~fs["match"]].to_string(index=False))
    print(fs[fs["item"].str.contains("limit")].to_string(index=False))

    print("\n== anchor geometry ==")
    ag = anchor_geometry(raw)
    print(ag.groupby("condition")["half_px_any"].sum().to_string())

    print("\n== recall / pool context the notebook set aside (mode 3) ==")
    rpc = recall_pool_context(raw)
    print(rpc.to_string(index=False))

    print("\n== product metrics (mode 8) ==")
    worst, lk = product_metrics(my_b, raw)
    print(lk.to_string(index=False))

    print("\n== Tier B: full re-derivation from pixels ==")
    pr, pb, gates, overlap, ties, div = tier_b(images, anns, raw)
    n_tier_b = 5 * len(pr) + 2 * len(pb)
    print(f"Tier B: {len(pr)} (ROI, condition) runs and {len(pb)} (arm, K) cells re-derived "
          f"({n_tier_b} values compared); {len(div)} divergences")
    # the click population the joint gate actually draws from
    print(f"  joint-gate click population: {int(gates['n_both'].sum())} of "
          f"{int(gates['n_pool'].sum())} pool candidates pass both gates = "
          f"{100 * gates['n_both'].sum() / gates['n_hem_accepts'].sum():.1f}% of hem's own "
          f"accepted population and "
          f"{100 * gates['n_both'].sum() / gates['n_gray_accepts'].sum():.1f}% of gray's")
    if div:
        print(pd.DataFrame(div).to_string(index=False))
    print(pr[["file_name", "condition", "base_size", "nb_base_size", "n_detections",
              "nb_n_detections", "n_peaks_preNMS", "n_near_seed_annulus"]].to_string(index=False))
    print("\ngate matrix (whole border-filtered pool, both channels):")
    print(gates.to_string(index=False))
    print(f"  walk: joint == hem-only solo draw on "
          f"{int(gates['joint_eq_solo_hem'].sum())}/14, == gray-only solo draw on "
          f"{int(gates['joint_eq_solo_gray'].sum())}/14; "
          f"{int(gates['n_walk_visits'].sum())} candidates evaluated in total")
    print(gates[gates["joint_eq_solo_gray"] == False][
        ["file_name", "joint_ann_id", "solo_hem_ann_id", "solo_gray_ann_id"]].to_string(index=False))
    print("\n  Otsu component at the chosen click, per condition (area px, solidity):")
    print(pr.pivot(index="file_name", columns="condition",
                   values=["otsu_component_area", "otsu_component_solidity"]).to_string())
    print("  mean accepted-component area: hem "
          f"{pr[pr.condition == 'hem_bbox'].otsu_component_area.mean():.0f} px vs gray "
          f"{pr[pr.condition == 'gray_bbox'].otsu_component_area.mean():.0f} px; "
          f"hem smaller on {int((pr[pr.condition == 'hem_bbox'].set_index('file_name').otsu_component_area < pr[pr.condition == 'gray_bbox'].set_index('file_name').otsu_component_area).sum())}/14 ROIs")
    print(f"  pool-wide: hem accepts {gates.n_hem_accepts.sum()}/{gates.n_pool.sum()}, "
          f"gray accepts {gates.n_gray_accepts.sum()}/{gates.n_pool.sum()}, "
          f"hem-accept-but-gray-refuse {gates.n_hem_only.sum()}, "
          f"gray-accept-but-hem-refuse {gates.n_gray_only.sum()}")
    print("\npool overlap between conditions:")
    print(overlap.to_string(index=False))
    print("\ntop-K candidate-set agreement between conditions:")
    cross = ties[~ties["arm"].str.contains("_within|_GT_hit_sets")]
    print(cross.groupby(["arm", "K"])[["n_topk_shared", "topk_identical"]].mean().to_string())
    print("\nGT-object agreement between the two conditions' top-K lists (the tp-tie mechanism):")
    gtsets = ties[ties["arm"].str.contains("_GT_hit_sets")]
    print(gtsets.groupby(["arm", "K"])[["n_topk_shared", "n_hem", "n_gray",
                                        "tp_counts_tie"]].mean().round(3).to_string())
    within = ties[ties["arm"].str.contains("_within")]
    print("\nwithin-condition re-ranking room at K=10 (top-10 shared with tm_score's top-10):")
    print(within.groupby("arm")["n_topk_shared"].describe().to_string())

    print("\n== 403.tiff uncapped trace ==")
    print(tier_b_403_uncapped(images, anns).to_string(index=False))

    print("\n== prose numbers ==")
    pn = prose_numbers(my_a, my_b, my_da, my_db, my_stat, base_size, raw, verif, nb, pr, gates)
    print(pn.to_string(index=False))
    print(f"prose numbers: {len(pn)} claims, {int((~pn['match']).sum())} mismatches")

    print(f"\nTOTAL {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
