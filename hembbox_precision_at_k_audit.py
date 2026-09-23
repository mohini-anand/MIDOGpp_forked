"""Independent audit of `production_seed_precision_at_k/production_seed_precision_at_k_chromatin_hem_bbox.ipynb`.

Regenerates every table in `Research Logs/2026-09-16-hembbox-precision-at-k-audit.md`.
Run with the anaconda interpreter:

    /Users/mohinianand/anaconda3/bin/python3 hembbox_precision_at_k_audit.py

Tier A work reads the persisted CSVs and `databases/MIDOG++.json` directly. Tier B work
re-runs the full 14-ROI pipeline from pixels with a locally written matcher, NMS, greedy
matcher and precision@K, so no midog_utils helper is its own oracle.
"""

from __future__ import annotations

import ast
import itertools
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import cast

import cv2
import numpy as np
import pandas as pd
import tifffile
from scipy import stats
from skimage.color import rgb2hed
from skimage.measure import label, regionprops

REPO = Path(__file__).resolve().parent
RES = REPO / "results"


def find_notebook(stem: str) -> Path:
    """Locate a notebook by filename anywhere in the repo -- this repo reorganises folders mid-run."""
    hits = [q for q in REPO.rglob(stem) if ".ipynb_checkpoints" not in str(q)]
    if not hits:
        raise FileNotFoundError(f"{stem} is not anywhere under {REPO}")
    return sorted(hits)[0]


NB = find_notebook("production_seed_precision_at_k_chromatin_hem_bbox.ipynb")
BASE_NB = find_notebook("production_seed_precision_at_k_chromatin_half_pix_fix.ipynb")
IMAGES_DIR = REPO / "images/extra_valid"
OUT = "hembbox_precision_at_k_audit"

NEW_RAW = RES / "precision_at_k_14roi_prodseed_chromatin_hembbox_raw.csv"
NEW_PER_ROI = RES / "precision_at_k_14roi_prodseed_chromatin_hembbox_per_roi.csv"
NEW_BY_DOMAIN = RES / "precision_at_k_14roi_prodseed_chromatin_hembbox_by_domain.csv"
NEW_VERIF = RES / "precision_at_k_14roi_prodseed_chromatin_hembbox_verification.csv"
NEW_DELTA_A = RES / "precision_at_k_14roi_prodseed_chromatin_hembbox_vs_halfpixfix_per_roi.csv"
NEW_DELTA_B = RES / "precision_at_k_14roi_prodseed_chromatin_hembbox_vs_halfpixfix_by_domain.csv"

OLD_RAW = RES / "precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv"
OLD_PER_ROI = RES / "precision_at_k_14roi_prodseed_chromatin_halfpixfix_per_roi.csv"
OLD_BY_DOMAIN = RES / "precision_at_k_14roi_prodseed_chromatin_halfpixfix_by_domain.csv"

BUDGETS = (10, 20, 30, 50)
AXES = {"tm_score": "score", "chromatin_od": "od51"}
SEED_INDEX = 0
PATCH_SIZE = 73
BORDER = PATCH_SIZE // 2
OTSU_WINDOW = 51
OD_WINDOW = 51
OD_PAD = OD_WINDOW // 2 + 1
OD_FRAC = 0.10
PEAK_MIN_DISTANCE = 7
SELF_HIT_RADIUS = 5.0
DEEP_FLOOR_Z = -1.5
MATCH_RADIUS_UM = 7.5
MITOTIC, LOOKALIKE = 1, 2

TABLES: dict[str, pd.DataFrame] = {}


def save(name: str, df: pd.DataFrame) -> pd.DataFrame:
    """Persist one audit table and echo it."""
    path = RES / f"{OUT}_{name}.csv"
    df.to_csv(path, index=False)
    TABLES[name] = df
    print(f"\n--- {name} -> {path.name} ({len(df)} rows)")
    with pd.option_context("display.width", 200, "display.max_columns", 40):
        print(df.to_string(index=False))
    return df


def banner(text: str) -> None:
    """Print a section header."""
    print("\n" + "=" * 100)
    print(text)
    print("=" * 100)


# ---------------------------------------------------------------------------------------
# Part 0a -- provenance gate
# ---------------------------------------------------------------------------------------

def git(*args: str) -> str:
    """Run a git command in the repo and return stripped stdout."""
    return subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True).stdout.strip()


def strip_ast(src: str) -> str:
    """AST dump of a module with every docstring removed, so comments/docstrings do not count."""
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            body = node.body
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                node.body = body[1:] or [ast.Pass()]
    return ast.dump(ast.fix_missing_locations(tree))


def provenance_gate() -> pd.DataFrame:
    """Compare artifact mtimes and module code against the commit the baseline was written at."""
    baseline_commit = git("rev-list", "-1", "--before=2026-09-10 15:51", "HEAD")
    mods = ["channels", "chromatin", "dataset", "evaluate", "seed_selection",
            "template_match", "compare", "nms", "invariants", "find_and_suppress"]
    rows = []
    for m in mods:
        old = git("show", f"{baseline_commit}:midog_utils/{m}.py")
        new = (REPO / f"midog_utils/{m}.py").read_text()
        rows.append(dict(
            item=f"midog_utils/{m}.py",
            mtime=time.strftime("%Y-%m-%d %H:%M", time.localtime((REPO / f"midog_utils/{m}.py").stat().st_mtime)),
            last_commit=git("log", "-1", "--format=%h %ad", "--date=short", "--", f"midog_utils/{m}.py"),
            uncommitted=bool(git("status", "--porcelain", "--", f"midog_utils/{m}.py")),
            code_identical_to_baseline_run=(strip_ast(old) == strip_ast(new)),
        ))
    for p in [NB, NEW_RAW, NEW_PER_ROI, NEW_BY_DOMAIN, NEW_VERIF, NEW_DELTA_A, NEW_DELTA_B,
              OLD_RAW, OLD_PER_ROI, OLD_BY_DOMAIN]:
        tracked = bool(git("ls-files", "--", str(p.relative_to(REPO))))
        rows.append(dict(
            item=str(p.relative_to(REPO)),
            mtime=time.strftime("%Y-%m-%d %H:%M", time.localtime(p.stat().st_mtime)),
            last_commit=git("log", "-1", "--format=%h %ad", "--date=short", "--", str(p.relative_to(REPO))) or "UNTRACKED",
            uncommitted=(not tracked),
            code_identical_to_baseline_run=np.nan,
        ))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------
# Part 0b -- execution coherence
# ---------------------------------------------------------------------------------------

def execution_gate() -> pd.DataFrame:
    """Check execution_count contiguity, error outputs and unrun cells for the target notebook."""
    nb = json.loads(NB.read_text())
    code = [c for c in nb["cells"] if c["cell_type"] == "code"]
    counts = [c.get("execution_count") for c in code]
    contiguous = counts == list(range(1, len(code) + 1))
    n_err = sum(1 for c in code for o in (c.get("outputs") or []) if o.get("output_type") == "error")
    unrun = [i for i, c in enumerate(nb["cells"]) if c["cell_type"] == "code" and c.get("execution_count") is None]
    n_png = sum(1 for c in code for o in (c.get("outputs") or []) if "image/png" in o.get("data", {}))
    return pd.DataFrame([dict(
        n_cells=len(nb["cells"]), n_code_cells=len(code),
        execution_counts=str(counts), contiguous_from_1=contiguous,
        n_error_outputs=n_err, n_unrun_code_cells=len(unrun),
        last_cell_is_run=(nb["cells"][-1].get("execution_count") is not None),
        n_embedded_pngs=n_png,
    )])


# ---------------------------------------------------------------------------------------
# Part 0c -- composition gate
# ---------------------------------------------------------------------------------------

def composition_gate(raw, per_roi, by_domain, delta_a, delta_b, old_per_roi) -> pd.DataFrame:
    """Row counts, key uniqueness, stratum coverage and merge losses, against the design."""
    files = sorted(f for f in os.listdir(IMAGES_DIR) if f.endswith(".tiff"))
    checks = []

    def add(name, expected, got):
        checks.append(dict(check=name, expected=expected, got=got, passed=(expected == got)))

    add("raw rows = 14 ROI x 2 arms x 4 budgets", 112, len(raw))
    add("per_roi rows = 14 ROI x 2 arms", 28, len(per_roi))
    add("by_domain rows = 7 domains x 2 arms x 4 budgets", 56, len(by_domain))
    add("delta_a rows (inner merge, no loss)", 28, len(delta_a))
    add("delta_b rows (inner merge, no loss)", 56, len(delta_b))
    add("distinct ROIs on disk", 14, len(files))
    add("distinct ROIs in raw", 14, raw["file_name"].nunique())
    add("distinct ROIs in per_roi", 14, per_roi["file_name"].nunique())
    add("ROI set in raw == ROI set on disk", True, sorted(raw["file_name"].unique()) == files)
    add("distinct domains in raw", 7, raw["tumor_type"].nunique())
    add("every domain has exactly 2 ROIs", True,
        bool((raw.groupby("tumor_type")["file_name"].nunique() == 2).all()))
    add("raw duplicate (file,arm,budget) keys", 0, int(raw.duplicated(["file_name", "arm", "budget"]).sum()))
    add("per_roi duplicate (file,arm) keys", 0, int(per_roi.duplicated(["file_name", "arm"]).sum()))
    add("raw fully duplicated rows", 0, int(raw.duplicated().sum()))
    add("arms in raw", sorted(AXES), sorted(raw["arm"].unique()))
    add("budgets in raw", list(BUDGETS), sorted(raw["budget"].unique()))
    add("baseline per_roi ROIs available for merge", 14, old_per_roi["file_name"].nunique())
    add("baseline arms include both audited arms", True,
        set(AXES) <= set(old_per_roi["arm"].unique()))
    add("delta_a domain_new agrees with baseline domain", True,
        bool((delta_a["domain"].notna()).all()))
    add("no NaN in raw tp_at_budget", 0, int(raw["tp_at_budget"].isna().sum()))
    add("budget_delivered == budget everywhere", True, bool((raw["budget_delivered"] == raw["budget"]).all()))
    add("n_detections >= max budget everywhere", True, bool((raw["n_detections"] >= max(BUDGETS)).all()))
    return pd.DataFrame(checks)


# ---------------------------------------------------------------------------------------
# Part 1 -- Tier A: recompute Table A, Table B, the deltas, the win counts, the summary
# ---------------------------------------------------------------------------------------

def recompute_table_a(raw, per_roi) -> tuple[pd.DataFrame, int]:
    """Rebuild Table A from RAW and diff it cell by cell against the published per-ROI CSV."""
    rows = []
    for fn in sorted(raw["file_name"].unique()):
        for arm in AXES:
            sub = raw[(raw["file_name"] == fn) & (raw["arm"] == arm)].set_index("budget")
            r: dict[str, object] = dict(file_name=fn, arm=arm)
            for k in BUDGETS:
                r[f"budget_delivered_{k}"] = int(sub.loc[k, "budget_delivered"])
                r[f"tp_at_{k}"] = int(sub.loc[k, "tp_at_budget"])
                r[f"precision_at_{k}"] = round(float(sub.loc[k, "tp_at_budget"]) / float(sub.loc[k, "budget_delivered"]), 4)
            rows.append(r)
    mine = pd.DataFrame(rows)
    merged = mine.merge(per_roi, on=["file_name", "arm"], suffixes=("_audit", "_nb"))
    cols = [f"{p}_{k}" for k in BUDGETS for p in ("budget_delivered", "tp_at", "precision_at")]
    diffs, n_cmp = [], 0
    for c in cols:
        a, b = merged[f"{c}_audit"].to_numpy(float), merged[f"{c}_nb"].to_numpy(float)
        n_cmp += len(a)
        bad = np.where(~np.isclose(a, b, atol=1e-9))[0]
        for i in bad:
            diffs.append(dict(file_name=merged.loc[i, "file_name"], arm=merged.loc[i, "arm"],
                              column=c, audit=a[i], notebook=b[i]))
    print(f"Table A: {n_cmp} values compared, {len(diffs)} divergences")
    return pd.DataFrame(diffs), n_cmp


def recompute_table_b(raw, by_domain) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Rebuild Table B (pooled and worst-ROI) from RAW and diff against the published CSV."""
    raw = raw.copy()
    raw["p"] = raw["tp_at_budget"] / raw["budget_delivered"]
    rows = []
    for (dom, arm, k), g in raw.groupby(["tumor_type", "arm", "budget"]):
        tp, dl = int(g["tp_at_budget"].sum()), int(g["budget_delivered"].sum())
        worst = g.loc[g["p"].idxmin()]
        # tie-aware: is the worst-ROI label unique?
        n_at_min = int((g["p"] == g["p"].min()).sum())
        rows.append(dict(domain=dom, arm=arm, K=int(k), n_roi=len(g), tp_sum=tp, delivered_sum=dl,
                         precision_pooled=round(tp / dl, 4),
                         precision_worst_roi=round(float(worst["p"]), 4),
                         worst_roi_file=str(worst["file_name"]),
                         worst_roi_is_tied=bool(n_at_min > 1)))
    mine = pd.DataFrame(rows)
    m = mine.merge(by_domain, on=["domain", "arm", "K"], suffixes=("_audit", "_nb"))
    diffs, n_cmp = [], 0
    for c in ["n_roi", "tp_sum", "delivered_sum", "precision_pooled", "precision_worst_roi"]:
        a, b = m[f"{c}_audit"].to_numpy(float), m[f"{c}_nb"].to_numpy(float)
        n_cmp += len(a)
        for i in np.where(~np.isclose(a, b, atol=1e-9))[0]:
            diffs.append(dict(domain=m.loc[i, "domain"], arm=m.loc[i, "arm"], K=m.loc[i, "K"],
                              column=c, audit=a[i], notebook=b[i]))
    n_cmp += len(m)
    for i in range(len(m)):
        if m.loc[i, "worst_roi_file_audit"] != m.loc[i, "worst_roi_file_nb"]:
            diffs.append(dict(domain=m.loc[i, "domain"], arm=m.loc[i, "arm"], K=m.loc[i, "K"],
                              column="worst_roi_file", audit=m.loc[i, "worst_roi_file_audit"],
                              notebook=m.loc[i, "worst_roi_file_nb"]))
    print(f"Table B: {n_cmp} values compared, {len(diffs)} divergences")
    return mine, pd.DataFrame(diffs), n_cmp


def recompute_deltas(per_roi, old_per_roi, by_domain, old_by_domain, delta_a, delta_b, same_seed_rois):
    """Rebuild the new-minus-old per-ROI and per-domain deltas and diff against the published CSVs."""
    old_a = old_per_roi[old_per_roi["arm"].isin(AXES)]
    m = per_roi.merge(old_a, on=["file_name", "arm"], suffixes=("_new", "_old"))
    for k in BUDGETS:
        m[f"delta_precision_at_{k}"] = m[f"precision_at_{k}_new"] - m[f"precision_at_{k}_old"]
    m["same_seed"] = m["file_name"].isin(same_seed_rois)

    d = m.merge(delta_a, on=["file_name", "arm"], suffixes=("_audit", "_nb"))
    diffs, n_cmp = [], 0
    for k in BUDGETS:
        for pre in ("precision_at_%d_new", "precision_at_%d_old", "delta_precision_at_%d"):
            c = pre % k
            a, b = d[f"{c}_audit"].to_numpy(float), d[f"{c}_nb"].to_numpy(float)
            n_cmp += len(a)
            for i in np.where(~np.isclose(a, b, atol=1e-9))[0]:
                diffs.append(dict(unit=f"{d.loc[i,'file_name']}/{d.loc[i,'arm']}", column=c,
                                  audit=a[i], notebook=b[i]))
    a, b = d["same_seed_audit"].to_numpy(), d["same_seed_nb"].to_numpy()
    n_cmp += len(a)
    for i in np.where(a != b)[0]:
        diffs.append(dict(unit=f"{d.loc[i,'file_name']}/{d.loc[i,'arm']}", column="same_seed",
                          audit=a[i], notebook=b[i]))

    old_b = old_by_domain[old_by_domain["arm"].isin(AXES)]
    mb = by_domain.merge(old_b, on=["domain", "arm", "K"], suffixes=("_new", "_old"))
    mb["delta_precision_pooled"] = mb["precision_pooled_new"] - mb["precision_pooled_old"]
    db = mb.merge(delta_b, on=["domain", "arm", "K"], suffixes=("_audit", "_nb"))
    for c in ["precision_pooled_new", "precision_pooled_old", "delta_precision_pooled"]:
        a, b = db[f"{c}_audit"].to_numpy(float), db[f"{c}_nb"].to_numpy(float)
        n_cmp += len(a)
        for i in np.where(~np.isclose(a, b, atol=1e-9))[0]:
            diffs.append(dict(unit=f"{db.loc[i,'domain']}/{db.loc[i,'arm']}/K{db.loc[i,'K']}",
                              column=c, audit=a[i], notebook=b[i]))
    print(f"Deltas: {n_cmp} values compared, {len(diffs)} divergences")
    return m, mb, pd.DataFrame(diffs), n_cmp


def seed_compare(new_raw, old_raw) -> pd.DataFrame:
    """Per-ROI seed identity and base_size, new run vs the half_pix_fix baseline."""
    new = new_raw.drop_duplicates("file_name").set_index("file_name")[["tumor_type", "seed_ann_id", "base_size"]]
    old = old_raw.drop_duplicates("file_name").set_index("file_name")[["seed_ann_id", "base_size"]]
    old.columns = ["old_seed_ann_id", "old_base_size"]
    out = new.join(old).reset_index()
    out["same_seed"] = out["seed_ann_id"] == out["old_seed_ann_id"]
    out["delta_base_size"] = out["base_size"] - out["old_base_size"]
    out["pct_change_base_size"] = (100 * out["delta_base_size"] / out["old_base_size"]).round(1)
    return out


# ---------------------------------------------------------------------------------------
# Part 2 -- inference: ROI-level exact tests, effective units, leave-one-ROI-out
# ---------------------------------------------------------------------------------------

def exact_sign_flip(d: np.ndarray) -> tuple[float, int]:
    """Two-sided exact sign-flip p for the mean of paired deltas; returns (p, n_nonzero)."""
    nz = d[d != 0]
    k = len(nz)
    if k == 0:
        return float("nan"), 0
    obs = abs(nz.sum())
    hits = 0
    for signs in itertools.product((1, -1), repeat=k):
        if abs(float(np.dot(signs, nz))) >= obs - 1e-12:
            hits += 1
    return hits / 2 ** k, k


def inference_table(delta_a_audit) -> pd.DataFrame:
    """Per (arm, K): mean delta, exact sign-flip p, exact sign-test p, t CI, LOO range."""
    rows = []
    for subset_name, mask in [("same_seed (n=11)", delta_a_audit["same_seed"]),
                              ("all ROIs (n=14)", pd.Series(True, index=delta_a_audit.index))]:
        for arm in AXES:
            sub = delta_a_audit[mask & (delta_a_audit["arm"] == arm)]
            for k in BUDGETS:
                d = sub[f"delta_precision_at_{k}"].to_numpy(float)
                g = len(d)
                wins, losses, ties = int((d > 0).sum()), int((d < 0).sum()), int((d == 0).sum())
                p_flip, k_eff = exact_sign_flip(d)
                p_sign = (stats.binomtest(wins, wins + losses, 0.5).pvalue
                          if (wins + losses) else float("nan"))
                sd = d.std(ddof=1) if g > 1 else float("nan")
                se = sd / np.sqrt(g) if g > 1 else float("nan")
                tcrit = stats.t.ppf(0.975, g - 1) if g > 1 else float("nan")
                loo = np.array([np.delete(d, i).mean() for i in range(g)])
                rows.append(dict(
                    subset=subset_name, arm=arm, K=k, n_roi=g,
                    mean_delta=round(float(d.mean()), 4),
                    wins=wins, losses=losses, ties=ties,
                    n_effective_units=k_eff,
                    p_floor_sign_flip=(2.0 ** (1 - k_eff) if k_eff else float("nan")),
                    p_exact_sign_flip=(round(p_flip, 4) if k_eff else float("nan")),
                    p_exact_sign_test=(round(float(p_sign), 4) if (wins + losses) else float("nan")),
                    ci95_lo=round(float(d.mean() - tcrit * se), 4),
                    ci95_hi=round(float(d.mean() + tcrit * se), 4),
                    loo_mean_min=round(float(loo.min()), 4),
                    loo_mean_max=round(float(loo.max()), 4),
                    loo_most_influential_roi=str(sub["file_name"].to_numpy()[int(np.argmax(np.abs(loo - d.mean())))]),
                    loo_sign_ever_flips=bool((np.sign(loo) != np.sign(d.mean())).any() and d.mean() != 0),
                ))
    return pd.DataFrame(rows)


def win_counts(delta_a_audit) -> pd.DataFrame:
    """Recompute the win/loss/tie counts the third figure plots (never persisted by the notebook)."""
    same = delta_a_audit[delta_a_audit["same_seed"]]
    rows = []
    for arm in AXES:
        sub = same[same["arm"] == arm]
        for k in BUDGETS:
            d = sub[f"delta_precision_at_{k}"]
            rows.append(dict(arm=arm, K=k, n_roi=len(d), wins=int((d > 0).sum()),
                             losses=int((d < 0).sum()), ties=int((d == 0).sum())))
    return pd.DataFrame(rows)


def mechanism_test(delta_a_audit, seedcmp, n_perm=2000) -> pd.DataFrame:
    """Does the base_size shrink explain the precision delta? Spearman + permutation p."""
    rng = np.random.default_rng(0)
    rows = []
    for subset_name, files in [("all ROIs", seedcmp["file_name"].tolist()),
                               ("same_seed only", seedcmp.loc[seedcmp["same_seed"], "file_name"].tolist())]:
        s = seedcmp[seedcmp["file_name"].isin(files)].set_index("file_name")
        for arm in AXES:
            sub = delta_a_audit[(delta_a_audit["arm"] == arm) & (delta_a_audit["file_name"].isin(files))]
            for k in BUDGETS:
                x = s.loc[sub["file_name"], "delta_base_size"].to_numpy(float)
                y = sub[f"delta_precision_at_{k}"].to_numpy(float)
                rho = stats.spearmanr(x, y).statistic
                rx, ry = stats.rankdata(x), stats.rankdata(y)
                null = np.array([stats.spearmanr(rx, rng.permutation(ry)).statistic for _ in range(n_perm)])
                p = float((np.abs(null) >= abs(rho) - 1e-12).mean())
                rows.append(dict(subset=subset_name, arm=arm, K=k, n=len(x),
                                 spearman_rho=round(float(rho), 4), p_permutation_2000=round(p, 4)))
    return pd.DataFrame(rows)


def recall_pool_triad(new_raw, old_raw) -> pd.DataFrame:
    """The context columns the notebook computed and set aside: pool size, coverage, tail recall."""
    cols = ["n_detections", "coverage_frac", "full_list_recall", "read_95", "read_100",
            "n_lookalike_in_list", "largest_tie_block", "nan_rate", "n_gt_mitotic"]
    new = new_raw[new_raw["budget"] == 10][["file_name", "arm"] + cols]
    old = old_raw[(old_raw["budget"] == 10) & (old_raw["arm"].isin(AXES))][["file_name", "arm"] + cols]
    m = new.merge(old, on=["file_name", "arm"], suffixes=("_new", "_old"))
    for c in cols:
        m[f"delta_{c}"] = m[f"{c}_new"] - m[f"{c}_old"]
    return m


def recall_pool_summary(triad) -> pd.DataFrame:
    """Per-arm summary of the pool/recall context, with an exact sign-flip on full_list_recall."""
    rows = []
    for arm in AXES:
        sub = triad[triad["arm"] == arm]
        for c in ["n_detections", "coverage_frac", "full_list_recall", "read_95", "read_100",
                  "n_lookalike_in_list"]:
            d = sub[f"delta_{c}"].to_numpy(float)
            d = np.nan_to_num(d, nan=0.0) if c.startswith("read") else d
            p, keff = exact_sign_flip(d)
            rows.append(dict(arm=arm, column=c,
                             mean_new=round(float(sub[f"{c}_new"].mean()), 4),
                             mean_old=round(float(sub[f"{c}_old"].mean()), 4),
                             mean_delta=round(float(np.nanmean(d)), 4),
                             n_worse=int((d < 0).sum()) if c in ("full_list_recall", "coverage_frac")
                             else int((d > 0).sum()),
                             n_effective_units=keff,
                             p_exact_sign_flip=round(p, 4) if keff else float("nan")))
    return pd.DataFrame(rows)


def family_annulus_check() -> pd.DataFrame:
    """Was this really the first notebook in the family to violate `seed_annulus_empty`?"""
    rows = []
    for p in sorted(RES.glob("precision_at_k_14roi_prodseed*_verification.csv")):
        v = pd.read_csv(p)
        sa = v[v["check"] == "seed_annulus_empty"]
        per_roi = sa[sa["label"] != "ALL"]
        rows.append(dict(verification_csv=p.name, n_rows=len(v),
                         has_seed_annulus_check=bool(len(sa)),
                         n_per_roi_rows=len(per_roi),
                         n_violations=int((~per_roi["passed"].astype(bool)).sum()) if len(per_roi) else 0,
                         violating_labels=";".join(per_roi.loc[~per_roi["passed"].astype(bool), "label"].astype(str))))
    return pd.DataFrame(rows)


def config_drift() -> pd.DataFrame:
    """Notebook config vs FSConfig dataclass defaults vs the operative DECISIONS.md entry."""
    sys.path.insert(0, str(REPO))
    from midog_utils import find_and_suppress as fs
    from midog_utils import template_match as tm
    from midog_utils import evaluate as ev
    d = fs.FSConfig()
    rows = [
        dict(setting="channel (search)", notebook="hematoxylin_od", fsconfig_default=d.channel,
             decisions="D3: hematoxylin_od unclipped", agrees_with_decision=True),
        dict(setting="tm_method", notebook="cv2.TM_CCOEFF (=4)", fsconfig_default=f"{d.tm_method} (TM_CCOEFF_NORMED=5)",
             decisions="D1: TM_CCOEFF", agrees_with_decision=True),
        dict(setting="nms_radius", notebook=f"radius_px(mpp, {ev.MIDOG_RADIUS_UM})",
             fsconfig_default=str(d.nms_radius), decisions="D7: = match radius, 7.5 um",
             agrees_with_decision=True),
        dict(setting="match_radius", notebook=f"radius_px(mpp, {ev.MIDOG_RADIUS_UM})",
             fsconfig_default="n/a", decisions="D7: 7.5 um", agrees_with_decision=True),
        dict(setting="self_hit_radius", notebook="5.0 px (fixed)", fsconfig_default=str(d.self_hit_radius),
             decisions="not decided; fixed-px while match radius is mpp-scaled", agrees_with_decision=True),
        dict(setting="peak_min_distance", notebook=str(PEAK_MIN_DISTANCE), fsconfig_default=str(d.peak_min_distance),
             decisions="not decided", agrees_with_decision=True),
        dict(setting="MAX_PEAKS (pre-NMS)", notebook="2_000_000", fsconfig_default=str(d.max_peaks),
             decisions="D9 (2026-09-12): 100", agrees_with_decision=False),
        dict(setting="deep floor z", notebook=str(DEEP_FLOOR_Z), fsconfig_default=str(d.deep_floor_z),
             decisions="D9: unchanged at -1.5", agrees_with_decision=True),
        dict(setting="base_size source", notebook="tightened_template_box (Otsu on hematoxylin_od)",
             fsconfig_default=str(d.base_size), decisions="D8: accepted component, centre (x0+x1-1)/2",
             agrees_with_decision=True),
        dict(setting="otsu_window", notebook=str(OTSU_WINDOW), fsconfig_default=str(tm.BASE_SIZE),
             decisions="D8 default", agrees_with_decision=True),
        dict(setting="patch_size", notebook=str(PATCH_SIZE), fsconfig_default=str(d.patch_size),
             decisions="tm.PATCH_SIZE", agrees_with_decision=True),
        dict(setting="Arm(caps=...)", notebook="(MAX_PEAKS,) -- a PRE-NMS cap checked against a POST-NMS length",
             fsconfig_default="n/a", decisions="invariants.check_no_cap contract", agrees_with_decision=False),
    ]
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------
# Part 3 -- Tier B: full 14-ROI re-run from pixels, written locally
# ---------------------------------------------------------------------------------------

def load_json_annotations() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse MIDOG++.json here rather than through midog_utils, so the audit owns its ground truth."""
    raw = json.loads((REPO / "databases/MIDOG++.json").read_text())
    images = pd.DataFrame([{"image_id": im["id"], "file_name": im["file_name"],
                            "tumor_type": im["tumor_type"]} for im in raw["images"]])
    name = dict(zip(images["image_id"], images["file_name"]))
    rows = []
    for a in raw["annotations"]:
        x1, y1, x2, y2 = a["bbox"]
        votes = a.get("labels", [])
        rows.append(dict(ann_id=a["id"], image_id=a["image_id"], file_name=name[a["image_id"]],
                         cx=(x1 + x2) / 2.0, cy=(y1 + y2) / 2.0, category_id=a["category_id"],
                         n_votes=len(votes), n_mitotic_votes=sum(1 for v in votes if v == MITOTIC)))
    return images, pd.DataFrame(rows)


def roi_mpp(path: Path) -> float:
    """Microns per pixel from the TIFF resolution tags, honouring ResolutionUnit."""
    unit_um = {2: 25400.0, 3: 10000.0}
    with tifffile.TiffFile(str(path)) as tf:
        page = cast(tifffile.TiffPage, tf.pages[0])
        num, den = page.tags["XResolution"].value
        unit = int(page.tags["ResolutionUnit"].value)
    return unit_um[unit] / (num / den)


def read_patch(img, cx, cy, size):
    """Square patch centred on the rounded pixel, or None when it runs off the edge."""
    half = size // 2
    ix, iy = int(round(cx)), int(round(cy))
    h, w = img.shape[:2]
    if ix - half < 0 or iy - half < 0 or ix + half >= w or iy + half >= h:
        return None
    return np.ascontiguousarray(img[iy - half: iy + half + 1, ix - half: ix + half + 1])


def odd(n, minimum=5):
    """Round to an odd integer no smaller than `minimum`."""
    n = int(round(n))
    if n % 2 == 0:
        n += 1
    return max(minimum, n)


def tighten_box(patch):
    """Locally written `tighten_box_otsu`: binary Otsu, centre-containment gate, sanity gates."""
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)  # type: ignore[call-overload]
    _, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    labels = cast(np.ndarray, label(binary, connectivity=2))
    cy, cx = patch.shape[0] // 2, patch.shape[1] // 2
    lab = labels[cy, cx]
    if lab == 0:
        return None
    region = next(p for p in regionprops(labels) if p.label == lab)
    y0, x0, y1, x1 = region.bbox
    if not (y0 <= cy < y1 and x0 <= cx < x1):
        return None
    if region.area < 50 or region.area > 0.85 * patch.size or region.solidity < 0.5:
        return None
    return int(y0), int(y1), int(x0), int(x1)


def tightened_template(chan, cx, cy, window=OTSU_WINDOW):
    """Locally written `tightened_template_box`: (base_size, centre_x, centre_y) or None."""
    patch = read_patch(chan, cx, cy, window)
    if patch is None:
        return None
    box = tighten_box(patch)
    if box is None:
        return None
    y0, y1, x0, x1 = box
    half = window // 2
    ix, iy = int(round(cx)), int(round(cy))
    return odd(max(y1 - y0, x1 - x0)), ix - half + (x0 + x1 - 1) / 2.0, iy - half + (y0 + y1 - 1) / 2.0


def nms_local(centers, scores, radius):
    """Greedy distance NMS, best-first, ties keeping incoming order."""
    from sklearn.neighbors import KDTree
    tree = KDTree(centers)
    nb = tree.query_radius(centers, r=radius)
    order = np.argsort(-scores, kind="stable")
    dead = np.zeros(len(centers), bool)
    keep = []
    for i in order:
        if dead[i]:
            continue
        keep.append(i)
        dead[nb[i]] = True
    return np.asarray(keep, int)


def greedy_precision(det_xy, gt_xy, gt_cat, radius, budgets):
    """Locally written greedy one-to-one matcher; returns tp counts at each budget and full-list TP."""
    from sklearn.neighbors import KDTree
    claimed = np.full(len(gt_xy), -1, int)
    hit = np.zeros(len(det_xy), int)          # 1 = mitotic GT, 2 = look-alike, 0 = nothing
    nb = KDTree(gt_xy).query_radius(det_xy, r=radius) if len(gt_xy) else [[]] * len(det_xy)
    for i, cands in enumerate(nb):
        free = [g for g in cands if claimed[g] == -1]
        if not free:
            continue
        d = np.hypot(gt_xy[free, 0] - det_xy[i, 0], gt_xy[free, 1] - det_xy[i, 1])
        g = free[int(np.argmin(d))]
        claimed[g] = i
        hit[i] = 1 if gt_cat[g] == MITOTIC else 2
    tp_cum = np.cumsum(hit == 1)
    return {k: int(tp_cum[min(k, len(tp_cum)) - 1]) for k in budgets}, int(tp_cum[-1]) if len(tp_cum) else 0


def tier_b_run(images, anns) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Re-run the notebook's pipeline over all 14 ROIs from pixels, with locally written code."""
    meta_ix = images.set_index("file_name")[["image_id", "tumor_type"]]
    files = sorted(f for f in os.listdir(IMAGES_DIR) if f.endswith(".tiff"))
    rows, annulus_rows = [], []
    for fn in files:
        t0 = time.time()
        image_id = int(meta_ix.loc[fn, "image_id"])
        domain = meta_ix.loc[fn, "tumor_type"]
        path = IMAGES_DIR / fn
        bgr = cv2.imread(str(path), cv2.IMREAD_COLOR)
        rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
        mpp = roi_mpp(path)
        radius = MATCH_RADIUS_UM / mpp
        hem = rgb2hed(rgb.astype(np.float32) / 255.0)[:, :, 0].astype(np.float32)
        H, W = hem.shape[:2]
        del bgr, rgb

        gt = anns[anns["file_name"] == fn].reset_index(drop=True)
        mit = gt[gt["category_id"] == MITOTIC]
        unanimous = mit[mit["n_mitotic_votes"] == mit["n_votes"]]
        pool = unanimous if len(unanimous) else mit[
            (mit["n_votes"] > 0) & (mit["n_mitotic_votes"] / mit["n_votes"] >= 2 / 3)
            & (mit["n_mitotic_votes"] < mit["n_votes"])]
        ix, iy = np.rint(pool["cx"]).astype(int), np.rint(pool["cy"]).astype(int)
        pool = pool[(ix >= BORDER) & (ix <= W - 1 - BORDER) & (iy >= BORDER) & (iy <= H - 1 - BORDER)]

        rng = np.random.default_rng([SEED_INDEX, image_id])
        working, retries = pool.copy(), 0
        seed = spec = None
        while len(working):
            j = int(rng.integers(len(working)))
            row = working.iloc[j]
            got = tightened_template(hem, float(row["cx"]), float(row["cy"]))
            if got is not None and read_patch(hem, got[1], got[2], PATCH_SIZE) is not None:
                seed, spec = row, got
                break
            working = working.drop(working.index[j])
            retries += 1
        if seed is None or spec is None:
            raise RuntimeError(f"{fn}: no valid seed found among {len(pool)} candidates")
        base_size, tx, ty = spec

        patch = read_patch(hem, tx, ty, PATCH_SIZE)
        assert patch is not None, f"{fn}: seed patch at ({tx}, {ty}) exceeds ROI bounds"
        c = PATCH_SIZE // 2
        hb = base_size // 2
        tmpl = np.ascontiguousarray(patch[c - hb: c + hb + 1, c - hb: c + hb + 1], dtype=np.float32)
        pad = (tmpl.shape[0] - 1) // 2
        padded = cv2.copyMakeBorder(hem, pad, pad, pad, pad, borderType=cv2.BORDER_REPLICATE)
        fused = np.asarray(cv2.matchTemplate(np.ascontiguousarray(padded, np.float32), tmpl,
                                             cv2.TM_CCOEFF), dtype=np.float32)
        assert fused.shape == (H, W), f"{fn}: response shape {fused.shape} != ROI {(H, W)}"
        del padded

        sample = fused[::8, ::8].ravel()
        sample = sample[np.isfinite(sample)]
        med = float(np.median(sample))
        mad = float(1.4826 * np.median(np.abs(sample - med)))
        cut = med + DEEP_FLOOR_Z * mad
        k = 2 * PEAK_MIN_DISTANCE + 1
        dil = np.asarray(cv2.dilate(fused, np.ones((k, k), np.uint8)), dtype=np.float32)
        mask = (fused >= dil) & (fused >= cut)
        ys, xs = np.nonzero(mask)
        sc = fused[ys, xs]
        order = np.lexsort((ys, xs, -sc))
        centers = np.stack([xs[order], ys[order]], 1).astype(np.float64)
        sc = sc[order]
        del fused, dil, mask

        keep = nms_local(centers, sc, radius)
        cc, ss = centers[keep], sc[keep]
        ok = np.hypot(cc[:, 0] - tx, cc[:, 1] - ty) > SELF_HIT_RADIUS
        cc, ss = cc[ok], ss[ok]

        hem_pad = cv2.copyMakeBorder(hem, OD_PAD, OD_PAD, OD_PAD, OD_PAD, cv2.BORDER_REPLICATE)
        nk = max(1, int(OD_FRAC * OD_WINDOW * OD_WINDOW))
        od = np.empty(len(cc))
        for i, (x, y) in enumerate(zip(cc[:, 0] + OD_PAD, cc[:, 1] + OD_PAD)):
            patch_i = read_patch(hem_pad, x, y, OD_WINDOW)
            assert patch_i is not None, f"{fn}: OD window at ({x}, {y}) exceeds padded bounds"
            od[i] = np.partition(patch_i.ravel(), -nk)[-nk:].mean()
        del hem, hem_pad

        gt_eval = gt[gt["ann_id"] != int(seed["ann_id"])].reset_index(drop=True)
        gt_xy = gt_eval[["cx", "cy"]].to_numpy(float)
        gt_cat = gt_eval["category_id"].to_numpy()
        n_gt = int((gt_cat == MITOTIC).sum())

        # 403 annulus diagnostic: any survivor inside one match radius of the template centre
        d_tpl = np.hypot(cc[:, 0] - tx, cc[:, 1] - ty)
        for i in np.where(d_tpl <= radius)[0]:
            annulus_rows.append(dict(
                file_name=fn, base_size=base_size, match_radius_px=round(radius, 3),
                cand_x=cc[i, 0], cand_y=cc[i, 1],
                dist_to_template_centre=round(float(d_tpl[i]), 3),
                dist_to_click=round(float(np.hypot(cc[i, 0] - seed["cx"], cc[i, 1] - seed["cy"])), 3),
                inside_click_radius=bool(np.hypot(cc[i, 0] - seed["cx"], cc[i, 1] - seed["cy"]) <= radius),
                rank_tm_score=int(np.argsort(np.argsort(-ss, kind="stable"), kind="stable")[i]) + 1,
                rank_chromatin_od=int(np.argsort(np.argsort(-od, kind="stable"), kind="stable")[i]) + 1,
                n_detections=len(cc)))

        for arm, key in AXES.items():
            vals = ss if key == "score" else od
            idx = np.argsort(-vals, kind="stable")
            tps, full_tp = greedy_precision(cc[idx], gt_xy, gt_cat, radius, BUDGETS)
            r = dict(file_name=fn, domain=domain, arm=arm, base_size=base_size,
                     seed_ann_id=int(seed["ann_id"]), n_retries=retries,
                     n_gt_mitotic=n_gt, n_detections=len(cc),
                     match_radius_px=round(radius, 3), full_list_recall=round(full_tp / n_gt, 6))
            for kk in BUDGETS:
                r[f"tp_at_{kk}"] = tps[kk]
                r[f"precision_at_{kk}"] = round(tps[kk] / kk, 4)
            rows.append(r)
        print(f"  [tierB] {fn} base={base_size} n_det={len(cc)} n_gt={n_gt} [{time.time()-t0:.0f}s]", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(annulus_rows)


def subset_sensitivity(delta_a_audit) -> pd.DataFrame:
    """How much of the headline depends on dropping the 3 different-seed ROIs."""
    rows = []
    for arm in AXES:
        sub = delta_a_audit[delta_a_audit["arm"] == arm]
        for k in BUDGETS:
            same = sub.loc[sub["same_seed"], f"delta_precision_at_{k}"]
            diff = sub.loc[~sub["same_seed"], f"delta_precision_at_{k}"]
            allr = sub[f"delta_precision_at_{k}"]
            rows.append(dict(arm=arm, K=k,
                             mean_same_seed_n11=round(float(same.mean()), 4),
                             mean_diff_seed_n3=round(float(diff.mean()), 4),
                             mean_all_rois_n14=round(float(allr.mean()), 4),
                             reported_by_notebook="same_seed n=11",
                             sign_flips_between_subsets=bool(np.sign(same.mean()) != np.sign(allr.mean())
                                                            and allr.mean() != 0 and same.mean() != 0)))
    return pd.DataFrame(rows)


def figure_series(delta_a_audit, delta_b_audit, wc) -> pd.DataFrame:
    """Re-derive the exact series each of the three figures plots, from the persisted deltas."""
    dom_col = "domain_new" if "domain_new" in delta_a_audit.columns else "domain"
    order = delta_a_audit.drop_duplicates("file_name").sort_values([dom_col, "file_name"])["file_name"].tolist()
    rows = []
    for arm in AXES:
        sub = delta_a_audit[delta_a_audit["arm"] == arm].set_index("file_name").loc[order]
        for fn in order:
            for k in BUDGETS:
                rows.append(dict(figure="fig1_roi_delta_heatmap", arm=arm, row=fn, col=f"K={k}",
                                 value=round(float(sub.loc[fn, f"delta_precision_at_{k}"]), 4),
                                 annotation="starred" if not bool(sub.loc[fn, "same_seed"]) else ""))
    for arm in AXES:
        sub = delta_b_audit[delta_b_audit["arm"] == arm]
        for dom in sorted(sub["domain"].unique()):
            for k in BUDGETS:
                v = sub[(sub["domain"] == dom) & (sub["K"] == k)]["delta_precision_pooled"].iloc[0]
                rows.append(dict(figure="fig2_domain_delta_bars", arm=arm, row=dom, col=f"K={k}",
                                 value=round(float(v), 4), annotation=""))
    for _, r in wc.iterrows():
        rows.append(dict(figure="fig3_win_counts", arm=r["arm"], row="wins", col=f"K={r['K']}",
                         value=int(r["wins"]), annotation=f"ties={int(r['ties'])} (not plotted)"))
        rows.append(dict(figure="fig3_win_counts", arm=r["arm"], row="losses", col=f"K={r['K']}",
                         value=int(r["losses"]), annotation="stacked on wins"))
    return pd.DataFrame(rows)


def roi_403_selfpeak(images, anns) -> pd.DataFrame:
    """Trace the one seed-annulus violation: was the self-peak outscored and then NMS-suppressed?"""
    fn = "403.tiff"
    meta = images.set_index("file_name")
    image_id = int(meta.loc[fn, "image_id"])
    path = IMAGES_DIR / fn
    bgr = cv2.imread(str(path), cv2.IMREAD_COLOR)
    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
    mpp = roi_mpp(path)
    radius = MATCH_RADIUS_UM / mpp
    hem = rgb2hed(rgb.astype(np.float32) / 255.0)[:, :, 0].astype(np.float32)
    H, W = hem.shape[:2]
    del bgr, rgb
    gt = anns[anns["file_name"] == fn].reset_index(drop=True)
    mit = gt[gt["category_id"] == MITOTIC]
    pool = mit[mit["n_mitotic_votes"] == mit["n_votes"]]
    ix, iy = np.rint(pool["cx"]).astype(int), np.rint(pool["cy"]).astype(int)
    pool = pool[(ix >= BORDER) & (ix <= W - 1 - BORDER) & (iy >= BORDER) & (iy <= H - 1 - BORDER)]
    rng = np.random.default_rng([SEED_INDEX, image_id])
    working = pool.copy()
    seed = spec = None
    while len(working):
        j = int(rng.integers(len(working)))
        row = working.iloc[j]
        got = tightened_template(hem, float(row["cx"]), float(row["cy"]))
        if got is not None and read_patch(hem, got[1], got[2], PATCH_SIZE) is not None:
            seed, spec = row, got
            break
        working = working.drop(working.index[j])
    if seed is None or spec is None:
        raise RuntimeError(f"{fn}: no valid seed found among {len(pool)} candidates")
    base_size, tx, ty = spec
    patch = read_patch(hem, tx, ty, PATCH_SIZE)
    assert patch is not None, f"{fn}: seed patch at ({tx}, {ty}) exceeds ROI bounds"
    c, hb = PATCH_SIZE // 2, base_size // 2
    tmpl = np.ascontiguousarray(patch[c - hb: c + hb + 1, c - hb: c + hb + 1], dtype=np.float32)
    pad = (tmpl.shape[0] - 1) // 2
    padded = cv2.copyMakeBorder(hem, pad, pad, pad, pad, borderType=cv2.BORDER_REPLICATE)
    fused = np.asarray(cv2.matchTemplate(np.ascontiguousarray(padded, np.float32), tmpl, cv2.TM_CCOEFF), np.float32)
    sample = fused[::8, ::8].ravel()
    med = float(np.median(sample))
    mad = float(1.4826 * np.median(np.abs(sample - med)))
    cut = med + DEEP_FLOOR_Z * mad
    k = 2 * PEAK_MIN_DISTANCE + 1
    dil = np.asarray(cv2.dilate(fused, np.ones((k, k), np.uint8)), dtype=np.float32)
    mask = (fused >= dil) & (fused >= cut)
    ys, xs = np.nonzero(mask)
    sc = fused[ys, xs]
    order = np.lexsort((ys, xs, -sc))
    centers = np.stack([xs[order], ys[order]], 1).astype(np.float64)
    sc = sc[order]
    d_pre = np.hypot(centers[:, 0] - tx, centers[:, 1] - ty)
    i_self = int(np.argmin(d_pre))
    keep = nms_local(centers, sc, radius)
    survivor = None
    for i in keep:
        if d_pre[i] <= radius and d_pre[i] > SELF_HIT_RADIUS:
            survivor = i
            break
    return pd.DataFrame([dict(
        file_name=fn, base_size=base_size, match_radius_px=round(radius, 3),
        global_max_score=round(float(sc[0]), 2),
        global_max_dist_to_template_centre=round(float(d_pre[0]), 3),
        nearest_pre_nms_peak_dist=round(float(d_pre[i_self]), 3),
        nearest_pre_nms_peak_score=round(float(sc[i_self]), 2),
        self_peak_inside_self_hit_radius=bool(d_pre[i_self] <= SELF_HIT_RADIUS),
        self_peak_survived_nms=bool(i_self in set(keep.tolist())),
        survivor_dist_to_template_centre=round(float(d_pre[survivor]), 3) if survivor is not None else np.nan,
        survivor_score=round(float(sc[survivor]), 2) if survivor is not None else np.nan,
        survivor_outscores_self_peak=bool(sc[survivor] > sc[i_self]) if survivor is not None else False,
        survivor_to_self_peak_px=round(float(np.hypot(centers[survivor, 0] - centers[i_self, 0],
                                                      centers[survivor, 1] - centers[i_self, 1])), 3)
        if survivor is not None else np.nan,
        survivor_within_nms_radius_of_self_peak=bool(
            np.hypot(centers[survivor, 0] - centers[i_self, 0],
                     centers[survivor, 1] - centers[i_self, 1]) <= radius) if survivor is not None else False,
        margin_past_nms_radius_px=round(float(np.hypot(centers[survivor, 0] - centers[i_self, 0],
                                                       centers[survivor, 1] - centers[i_self, 1]) - radius), 3)
        if survivor is not None else np.nan,
        template_centre_xy=f"({tx}, {ty})",
        self_peak_xy=f"({centers[i_self, 0]}, {centers[i_self, 1]})",
        click_xy=f"({seed['cx']}, {seed['cy']})",
        self_hit_radius_px=SELF_HIT_RADIUS,
        gap_between_self_hit_and_match_radius_px=round(float(radius - SELF_HIT_RADIUS), 3),
    )])



def dropped_arms(old_raw) -> pd.DataFrame:
    """Pooled precision@K of every arm the baseline scored, showing which this notebook dropped."""
    rows = []
    for (arm, k), g in old_raw.groupby(["arm", "budget"]):
        rows.append(dict(arm=arm, K=int(k),
                         precision_pooled_14roi=round(float(g["tp_at_budget"].sum() / g["budget_delivered"].sum()), 4),
                         carried_into_hembbox_notebook=(arm in AXES)))
    out = pd.DataFrame(rows).pivot_table(index=["arm", "carried_into_hembbox_notebook"],
                                         columns="K", values="precision_pooled_14roi").reset_index()
    out.columns = [f"P@{c}" if isinstance(c, (int, np.integer)) else c for c in out.columns]
    return out.sort_values("P@10", ascending=False)


def prose_numbers(seedcmp, wc, annulus) -> pd.DataFrame:
    """Every number quoted in the notebook's prose that is not printed by any cell."""
    rows = []

    def add(claim, quoted, recomputed):
        rows.append(dict(prose_claim=claim, notebook_says=str(quoted), audit_finds=str(recomputed),
                         agrees=(str(quoted) == str(recomputed))))

    add("median old base_size", 39, int(seedcmp["old_base_size"].median()))
    add("median new base_size", 28, int(seedcmp["base_size"].median()))
    add("median per-ROI pct change", -13.3, round(float(seedcmp["pct_change_base_size"].median()), 1))
    add("n shrink", 13, int((seedcmp["delta_base_size"] < 0).sum()))
    add("n grow", 1, int((seedcmp["delta_base_size"] > 0).sum()))
    add("201.tiff base_size old->new", "51->29",
        f"{int(seedcmp.set_index('file_name').loc['201.tiff','old_base_size'])}->"
        f"{int(seedcmp.set_index('file_name').loc['201.tiff','base_size'])}")
    add("245.tiff base_size old->new", "47->23",
        f"{int(seedcmp.set_index('file_name').loc['245.tiff','old_base_size'])}->"
        f"{int(seedcmp.set_index('file_name').loc['245.tiff','base_size'])}")
    add("403.tiff base_size old->new", "51->31",
        f"{int(seedcmp.set_index('file_name').loc['403.tiff','old_base_size'])}->"
        f"{int(seedcmp.set_index('file_name').loc['403.tiff','base_size'])}")
    tm = wc[wc["arm"] == "tm_score"]
    add("tm_score same-seed wins range across K", "2-3",
        f"{int(tm['wins'].min())}-{int(tm['wins'].max())}")
    add("tm_score same-seed losses range across K", "6-8",
        f"{int(tm['losses'].min())}-{int(tm['losses'].max())}")
    ch = wc[wc["arm"] == "chromatin_od"]
    add("chromatin_od same-seed ties range across K", "7-10",
        f"{int(ch['ties'].min())}-{int(ch['ties'].max())}")
    a = annulus.iloc[0]
    add("403 annulus candidate distance to click (px)", 33.29, round(float(a["dist_to_click"]), 2))
    add("403 match_radius (px)", 33.06, round(float(a["match_radius_px"]), 2))
    add("403 annulus candidate rank, tm_score (ordinal)", "3281st", f"{int(a['rank_tm_score'])}th")
    add("403 annulus candidate rank, chromatin_od (ordinal)", "393rd", f"{int(a['rank_chromatin_od'])}th")
    add("403 pool size", 17626, int(a["n_detections"]))
    add("403 annulus candidate reaches K<=50", False, bool(int(a["rank_tm_score"]) <= 50 or int(a["rank_chromatin_od"]) <= 50))
    return pd.DataFrame(rows)


def main() -> None:
    """Run every gate, every Tier A recomputation and the Tier B re-run, and write the tables."""
    t_start = time.time()
    pd.set_option("display.width", 220)

    banner("PROVENANCE GATE")
    save("provenance", provenance_gate())

    banner("EXECUTION-COHERENCE GATE")
    save("execution_gate", execution_gate())

    new_raw = pd.read_csv(NEW_RAW)
    new_a = pd.read_csv(NEW_PER_ROI)
    new_b = pd.read_csv(NEW_BY_DOMAIN)
    d_a = pd.read_csv(NEW_DELTA_A)
    d_b = pd.read_csv(NEW_DELTA_B)
    old_raw = pd.read_csv(OLD_RAW)
    old_a = pd.read_csv(OLD_PER_ROI)
    old_b = pd.read_csv(OLD_BY_DOMAIN)

    banner("COMPOSITION GATE")
    save("composition", composition_gate(new_raw, new_a, new_b, d_a, d_b, old_a))

    banner("TIER A -- TABLE A / TABLE B / DELTAS")
    diffs_a, n_a = recompute_table_a(new_raw, new_a)
    mine_b, diffs_b, n_b = recompute_table_b(new_raw, new_b)
    save("table_b_recomputed", mine_b)

    seedcmp = seed_compare(new_raw, old_raw)
    save("seed_compare", seedcmp)
    same_rois = seedcmp.loc[seedcmp["same_seed"], "file_name"].tolist()

    delta_a_audit, delta_b_audit, diffs_d, n_d = recompute_deltas(
        new_a, old_a, new_b, old_b, d_a, d_b, same_rois)
    all_diffs = pd.concat([diffs_a, diffs_b, diffs_d], ignore_index=True) if \
        (len(diffs_a) or len(diffs_b) or len(diffs_d)) else pd.DataFrame(
            columns=["unit", "column", "audit", "notebook"])
    save("tier_a_divergences", all_diffs)
    print(f"\nTIER A TOTAL: {n_a + n_b + n_d} values compared, {len(all_diffs)} divergences")
    save("tier_a_counts", pd.DataFrame([dict(block="Table A", values=n_a, divergences=len(diffs_a)),
                                        dict(block="Table B", values=n_b, divergences=len(diffs_b)),
                                        dict(block="Deltas", values=n_d, divergences=len(diffs_d)),
                                        dict(block="TOTAL", values=n_a + n_b + n_d,
                                             divergences=len(all_diffs))]))

    banner("CLOSING-SUMMARY MEANS")
    wc = win_counts(delta_a_audit)
    save("win_counts_recomputed", wc)
    summ = []
    for arm in AXES:
        sub = delta_a_audit[delta_a_audit["same_seed"] & (delta_a_audit["arm"] == arm)]
        r: dict[str, object] = dict(arm=arm, n_roi=len(sub))
        for k in BUDGETS:
            r[f"mean_delta_K{k}"] = round(float(sub[f"delta_precision_at_{k}"].mean()), 4)
        summ.append(r)
    save("closing_summary_means", pd.DataFrame(summ))

    banner("INFERENCE -- ROI-LEVEL EXACT TESTS")
    save("inference", inference_table(delta_a_audit))

    banner("SUBSET SENSITIVITY -- same-seed vs all 14 ROIs")
    save("subset_sensitivity", subset_sensitivity(delta_a_audit))

    banner("FIGURE SERIES -- re-derived from the persisted deltas")
    save("figure_series", figure_series(delta_a_audit, delta_b_audit, wc))

    banner("MECHANISM -- base_size shrink vs precision delta")
    save("mechanism_spearman", mechanism_test(delta_a_audit, seedcmp))

    banner("RECALL / POOL CONTEXT (the triad the notebook set aside)")
    triad = recall_pool_triad(new_raw, old_raw)
    save("recall_pool_per_roi", triad)
    save("recall_pool_summary", recall_pool_summary(triad))

    banner("FAMILY seed_annulus_empty HISTORY")
    save("family_annulus", family_annulus_check())

    banner("ARMS DROPPED FROM THE BASELINE'S SIX")
    save("dropped_arms", dropped_arms(old_raw))

    banner("CONFIG DRIFT")
    save("config_drift", config_drift())

    banner("TIER B -- FULL 14-ROI RE-RUN FROM PIXELS")
    images, anns = load_json_annotations()
    tb, annulus = tier_b_run(images, anns)
    save("tier_b_per_roi", tb)
    save("tier_b_annulus", annulus)
    save("tier_b_403_selfpeak", roi_403_selfpeak(images, anns))
    save("prose_numbers", prose_numbers(seedcmp, wc, annulus))

    m = tb.merge(new_a, on=["file_name", "arm"], suffixes=("_tierb", "_nb"))
    tb_diffs, n_tb = [], 0
    for c in ["n_gt_mitotic", "n_detections"] + [f"tp_at_{k}" for k in BUDGETS] + \
             [f"precision_at_{k}" for k in BUDGETS]:
        a, b = m[f"{c}_tierb"].to_numpy(float), m[f"{c}_nb"].to_numpy(float)
        n_tb += len(a)
        for i in np.where(~np.isclose(a, b, atol=1e-9))[0]:
            tb_diffs.append(dict(unit=f"{m.loc[i,'file_name']}/{m.loc[i,'arm']}", column=c,
                                 tier_b=a[i], notebook=b[i]))
    seed_m = tb.drop_duplicates("file_name").set_index("file_name")[["seed_ann_id", "base_size"]].join(
        seedcmp.set_index("file_name")[["seed_ann_id", "base_size"]], rsuffix="_nb")
    for fn, sr in seed_m.iterrows():
        n_tb += 2
        if int(sr["seed_ann_id"]) != int(sr["seed_ann_id_nb"]):
            tb_diffs.append(dict(unit=fn, column="seed_ann_id", tier_b=sr["seed_ann_id"], notebook=sr["seed_ann_id_nb"]))
        if int(sr["base_size"]) != int(sr["base_size_nb"]):
            tb_diffs.append(dict(unit=fn, column="base_size", tier_b=sr["base_size"], notebook=sr["base_size_nb"]))
    save("tier_b_divergences", pd.DataFrame(tb_diffs) if tb_diffs else
         pd.DataFrame(columns=["unit", "column", "tier_b", "notebook"]))
    print(f"\nTIER B: {n_tb} values compared, {len(tb_diffs)} divergences")

    banner(f"DONE in {time.time() - t_start:.0f}s -- {len(TABLES)} tables written to results/{OUT}_*.csv")


if __name__ == "__main__":
    main()
