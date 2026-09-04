"""Summary metrics for bbox_threshold_sweep.py's per-candidate CSV.

Pairwise (binary vs. one arm at a time, never a many-way intersection -- a 7-way
intersection over a suspicious population would be close to a guaranteed null).
"""
from __future__ import annotations

import pandas as pd

ARMS = ["multiotsu", "frac_0.05", "frac_0.10", "frac_0.15", "frac_0.20", "frac_0.30"]

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 20)

df = pd.read_csv("results/bbox_threshold_sweep.csv")


def flagged_accepted_metrics(sub: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for arm in ARMS:
        acc = sub[f"{arm}_accepted"]
        merged = sub[f"{arm}_looks_merged"]
        n = len(sub)
        resolved = int((acc & (merged == False)).sum())  # noqa: E712
        still_merged = int((acc & (merged == True)).sum())  # noqa: E712
        excluded = int((~acc).sum())
        paired = sub[acc]  # binary already accepted (this is flagged_accepted by definition)
        sol_delta = (paired[f"{arm}_solidity"] - paired["bin_solidity"]).median() if len(paired) else float("nan")
        rows.append(dict(arm=arm, n=n, resolved=resolved, resolved_pct=round(100*resolved/n, 1),
                         still_merged=still_merged, still_merged_pct=round(100*still_merged/n, 1),
                         excluded=excluded, excluded_pct=round(100*excluded/n, 1),
                         median_solidity_delta=round(sol_delta, 3) if pd.notna(sol_delta) else None))
    return pd.DataFrame(rows)


def flagged_rejected_metrics(sub: pd.DataFrame) -> pd.DataFrame:
    rows = []
    n = len(sub)
    for arm in ARMS:
        acc = sub[f"{arm}_accepted"]
        merged = sub[f"{arm}_looks_merged"]
        recovered = int((acc & (merged == False)).sum())  # noqa: E712
        accepted_still_merged = int((acc & (merged == True)).sum())  # noqa: E712
        still_rejected = int((~acc).sum())
        rows.append(dict(arm=arm, n=n, recovered=recovered, accepted_still_merged=accepted_still_merged,
                         still_rejected=still_rejected))
    return pd.DataFrame(rows)


def control_metrics(sub: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for arm in ARMS:
        acc = sub[f"{arm}_accepted"]
        n = len(sub)
        excluded = int((~acc).sum())
        paired = sub[acc]
        shrink = None
        if len(paired):
            rel = (paired["bin_longer_side"] - paired[f"{arm}_longer_side"]) / paired["bin_longer_side"]
            shrink = round(100 * rel.median(), 1)
        rows.append(dict(arm=arm, n=n, collateral_excluded=excluded,
                         collateral_excluded_pct=round(100*excluded/n, 1),
                         median_longer_side_shrink_pct=shrink))
    return pd.DataFrame(rows)


print("=" * 100)
print("POPULATION (pooled)")
print("=" * 100)
print(df["bucket"].value_counts().to_string())
print()

print("=" * 100)
print(f"FLAGGED_ACCEPTED (n={len(df[df.bucket=='flagged_accepted'])}) -- accepted by binary Otsu, but "
      "area_frac>=0.5 or longer_side>=45px. Does raising the threshold resolve these?")
print("=" * 100)
fa = df[df.bucket == "flagged_accepted"]
print(flagged_accepted_metrics(fa).to_string(index=False))
print()

print("--- flagged_accepted, per domain (dense domains only: lymphosarcoma, mast cell tumor) ---")
for dom in ["canine lymphosarcoma", "canine cutaneous mast cell tumor"]:
    sub = fa[fa.tumor_type == dom]
    print(f"\n{dom} (n={len(sub)}):")
    print(flagged_accepted_metrics(sub).to_string(index=False))
print()

print("=" * 100)
print(f"FLAGGED_REJECTED (n={len(df[df.bucket=='flagged_rejected'])}) -- binary Otsu ALREADY rejects these "
      "(max_area_frac or min_solidity). Does raising the threshold recover a clean box?")
print("=" * 100)
fr = df[df.bucket == "flagged_rejected"]
if len(fr):
    print(flagged_rejected_metrics(fr).to_string(index=False))
    print("\n(individual rows, n is small -- inspect directly)")
    print(fr[["file_name", "ann_id", "bin_area_frac", "bin_longer_side", "bin_solidity", "bin_reject_reason"]].to_string(index=False))
else:
    print("(empty)")
print()

print("=" * 100)
print(f"CONTROL (n={len(df[df.bucket=='control'])}) -- already clean, area_frac<0.3. "
      "Collateral damage check: does raising the threshold needlessly exclude/shrink these?")
print("=" * 100)
ctrl = df[df.bucket == "control"]
print(control_metrics(ctrl).to_string(index=False))
print()

print("--- control, per domain (dense domains only) ---")
for dom in ["canine lymphosarcoma", "canine cutaneous mast cell tumor"]:
    sub = ctrl[ctrl.tumor_type == dom]
    print(f"\n{dom} (n={len(sub)}):")
    print(control_metrics(sub).to_string(index=False))

print()
print("=" * 100)
print("GAP ZONE population share (0.3<=area_frac<0.5 and longer_side<45, not otherwise analysed)")
print("=" * 100)
print(f"{len(df[df.bucket=='gap'])} / {len(df)} = {100*len(df[df.bucket=='gap'])/len(df):.1f}% of all candidates")
