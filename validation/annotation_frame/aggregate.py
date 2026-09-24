"""
    Combine a tier's per-ROI results: sum the per-check tallies, compare the recomputed R1 rows and R2 top-30 rows against the
    references for the covered ROIs, run rr.compare_selection against R3, and (with --git) V2's repository checks.
    Usage: aggregate.py --tier smoke|full [--git BASE]
"""

from __future__ import annotations

import argparse
import json
import subprocess

import pandas as pd

import common as C
from common import rr
import checks


def git(*args):
    """
        Run git in the worktree.

        args (str): git arguments.

        Returns str: stdout.
    """
    return subprocess.run(["git", *args], cwd=C.WORKTREE, capture_output=True, text=True, check=True).stdout


def git_checks(base):
    """
        V2's repository half: the diff against the base touches only production.py, git status lists only that and the new
        module under midog_utils/, and every other tracked file is byte-identical to the base.

        base (str): the base commit.

        Returns dict: the diff, status and the check results.
    """
    diff_names = git("diff", "--name-only", base, "--", "midog_utils/").split()
    status = git("status", "--short", "midog_utils/").splitlines()
    others = git("diff", "--name-only", base, "--", ".", ":(exclude)midog_utils/production.py").split()
    all_status = git("status", "--short", "--untracked-files=all").splitlines()
    return dict(diff_names=diff_names, status_midog_utils=status, other_tracked_changed=others, status_all=all_status, diff_production=git("diff", base, "--", "midog_utils/production.py"), ok=diff_names == ["midog_utils/production.py"] and sorted(status) == sorted([" M midog_utils/production.py", "?? midog_utils/annotation_frame.py"]) and others == [])


def main():
    """
        Aggregate a tier and write out/<tier>/summary.json.

        Returns None.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--tier", required=True, choices=("smoke", "full"))
    parser.add_argument("--git", default="", help="base commit for V2's repository checks")
    args = parser.parse_args()

    rois, _ = rr.build_rois()
    roi_dir = C.OUT / args.tier / "rois"
    results = [C.load_pickle(roi_dir / f"{fn}.pkl") for fn in rois["file_name"] if (roi_dir / f"{fn}.pkl").exists()]
    covered = [r["file_name"] for r in results]
    tally = checks.merge_tallies([r["tally"] for r in results])

    r1 = C.read_ref(C.R1_PATH)
    per_run = C.finish_per_run([row for r in results for row in r["r1_rows"]], rois)
    r1_sub = r1[r1["file_name"].isin(covered)].reset_index(drop=True)
    r1_cols = [c for c in r1.columns if not c.startswith("t_")]
    key = ["file_name", "seed_index", "condition"]
    assert list(per_run.columns) == r1_cols, f"recomputed R1 columns differ:\n{list(per_run.columns)}\n{r1_cols}"
    assert per_run[key].equals(r1_sub[key]), "R1 key order differs"
    res = C.compare_frames(per_run, r1_sub, r1_cols, label="V3 R1 recompute")
    tally["V3 R1 recomputed from frame/info/draw_info/evaluate"] = dict(n_runs=len(per_run), n_rows=res["n_rows"], n_values=res["n_values"], n_diverged=res["n_diverged"], detail=res["detail"], columns=len(r1_cols))

    r2 = C.read_ref(C.R2_PATH)
    const = (r2["angle"] == 0.0).all() and (r2["flip"] == False).all() and (r2["scale"] == 1.0).all()  # noqa: E712
    tally["V3 R2 angle/flip/scale are constants"] = dict(n_runs=0, n_rows=len(r2), n_values=3 * len(r2), n_diverged=0 if const else int((r2["angle"] != 0.0).sum() + (r2["flip"] != False).sum() + (r2["scale"] != 1.0).sum()), detail=[] if const else ["R2 angle/flip/scale not constant"])  # noqa: E712
    top30 = pd.concat([f for r in results for f in r["r2_frames"]], ignore_index=True)
    r2_sub = r2[r2["file_name"].isin(covered)].reset_index(drop=True)
    r2_cols = [c for c in r2.columns if c not in checks.R2_SKIPPED]
    key2 = ["file_name", "seed_index", "condition", "rank"]
    assert top30[key2].equals(r2_sub[key2]), "R2 key order differs"
    res2 = C.compare_frames(top30, r2_sub, r2_cols, float32_cols=("score",), label="V3 R2 recompute")
    tally["V3 R2 recomputed from frame + evaluate"] = dict(n_runs=int(top30.groupby(["file_name", "seed_index", "condition"]).ngroups), n_rows=res2["n_rows"], n_values=res2["n_values"], n_diverged=res2["n_diverged"], detail=res2["detail"], columns=len(r2_cols))

    r3 = C.read_ref(C.R3_PATH)
    sel = pd.DataFrame([row for r in results for row in r["sel_rows"]])
    cmp = rr.compare_selection(sel, r3)
    tally["V1 compare_selection vs R3"] = dict(n_runs=len(cmp), n_rows=len(cmp), n_values=len(cmp) * len(rr.SELECTION_COLS), n_diverged=int((~cmp["passed"]).sum()), detail=cmp.loc[~cmp["passed"], "detail"].tolist()[:20], n_passed=int(cmp["passed"].sum()))

    n_half = sum(r["n_half_even_vs_half_up_differ"] for r in results)
    summary = dict(tier=args.tier, n_rois=len(covered), rois=covered, n_half_even_vs_half_up_differ=n_half, checks=tally)
    if args.git:
        summary["v2_git"] = git_checks(args.git)
    if args.tier == "full":
        want_runs = {"V3 predictions vs B1": 441, "V3 production info vs B1": 441, "V3 row count 2 + n_detections": 441, "V3 refined row vs R1": 441, "V10 metrics vs evaluate_run on B1": 441, "V10 tp_at_k vs R1": 441, "V4 coordinate set, chromatin_od vs tm_score": 441, "V4 score/od at each coordinate": 441, "V4 tm_score predictions vs B2": 126, "V4 default_51 tm_score top30 vs R4": 147, "V4 default_51 tm_score per-run vs R4": 147, "V2 default always_od vs B2": 126, "V2 default always_od vs B2 (info)": 126, "V2 always_od=True vs B2": 126, "V2 always_od=True vs B2 (info)": 126, "V8 per-frame schema": 882, "V1 draw vs select_clicks": 147, "V1 seeds 0-5 vs R4 clicks": 294}
        completeness = {"n_rois == 49": len(covered) == 49, "R1 recompute rows == 441": len(per_run) == 441, "R2 recompute rows == 13230": len(top30) == 13230, "compare_selection passed == 441": tally["V1 compare_selection vs R3"]["n_passed"] == 441 and len(cmp) == 441}
        completeness.update({f"{name} runs == {n}": tally.get(name, {}).get("n_runs") == n for name, n in want_runs.items()})
        summary["completeness"] = completeness
        print("completeness:", "all ok" if all(completeness.values()) else {k: v for k, v in completeness.items() if not v})
    (C.OUT / args.tier).mkdir(parents=True, exist_ok=True)
    (C.OUT / args.tier / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
    print(f"{args.tier}: {len(covered)} ROIs; runs with round-half-even != round-half-up on a template axis: {n_half}")
    for name in sorted(tally):
        e = tally[name]
        print(f"  {name}: runs={e['n_runs']} rows={e['n_rows']} values={e['n_values']} diverged={e['n_diverged']}")
        for line in e["detail"][:5]:
            print(f"      {line}")
    if args.git:
        g = summary["v2_git"]
        print(f"  V2 git: ok={g['ok']} diff_names={g['diff_names']} status={g['status_midog_utils']} other_tracked_changed={g['other_tracked_changed']}")
    print(f"-> {C.OUT / args.tier / 'summary.json'}")


if __name__ == "__main__":
    main()
