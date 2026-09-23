"""
    Round-2 audit of `production_seed_precision_at_k_chromatin_hem_bbox.ipynb`.

    Round 1 (`Research Logs/2026-09-16-hembbox-precision-at-k-audit.md`) raised 16
    "CORRECTION TO BE MADE" blocks plus a "Correction dependencies" section. A later pass
    edited the notebook and re-executed it (13:24). This script verifies, for every one of
    those 16 items and for the dependency section's requirements, whether the correction was
    applied, applied correctly and applied completely -- and re-derives every number the
    corrected notebook now asserts, from the raw CSVs and from the ROI pixels.

    Writes results/hembbox_precision_at_k_audit_round2_*.csv. Run under
    /Users/mohinianand/anaconda3/bin/python3.
"""

from __future__ import annotations

import itertools
import json
import os
import re
import subprocess
import time

import cv2
import numpy as np
import pandas as pd
from scipy import stats as sps
from skimage.color import rgb2hed
from skimage.measure import label, regionprops
from sklearn.neighbors import KDTree

REPO = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(REPO, 'results')
OUT = os.path.join(RES, 'hembbox_precision_at_k_audit_round2_{}.csv')

NB = next(p for p in
          (os.path.join(dp, f) for dp, _, fs in os.walk(REPO) for f in fs)
          if p.endswith('production_seed_precision_at_k_chromatin_hem_bbox.ipynb')
          and '.audit_tmp' not in p)

NEW_RAW = os.path.join(RES, 'precision_at_k_14roi_prodseed_chromatin_hembbox_raw.csv')
NEW_PER_ROI = os.path.join(RES, 'precision_at_k_14roi_prodseed_chromatin_hembbox_per_roi.csv')
NEW_BY_DOMAIN = os.path.join(RES, 'precision_at_k_14roi_prodseed_chromatin_hembbox_by_domain.csv')
NEW_VERIF = os.path.join(RES, 'precision_at_k_14roi_prodseed_chromatin_hembbox_verification.csv')
NEW_STAT = os.path.join(RES, 'precision_at_k_14roi_prodseed_chromatin_hembbox_stat_context.csv')
NEW_DELTA_A = os.path.join(RES, 'precision_at_k_14roi_prodseed_chromatin_hembbox_vs_halfpixfix_per_roi.csv')
NEW_DELTA_B = os.path.join(RES, 'precision_at_k_14roi_prodseed_chromatin_hembbox_vs_halfpixfix_by_domain.csv')
OLD_RAW_P = os.path.join(RES, 'precision_at_k_14roi_prodseed_chromatin_halfpixfix_raw.csv')
OLD_PER_ROI = os.path.join(RES, 'precision_at_k_14roi_prodseed_chromatin_halfpixfix_per_roi.csv')

R1_TIER_B = os.path.join(RES, 'hembbox_precision_at_k_audit_tier_b_per_roi.csv')
R1_SPEARMAN = os.path.join(RES, 'hembbox_precision_at_k_audit_mechanism_spearman.csv')
R1_INFERENCE = os.path.join(RES, 'hembbox_precision_at_k_audit_inference.csv')
R1_RECALL = os.path.join(RES, 'hembbox_precision_at_k_audit_recall_pool_summary.csv')

IMAGES_DIR = os.path.join(REPO, 'images', 'extra_valid')
DB = os.path.join(REPO, 'databases', 'MIDOG++.json')

BUDGETS = (10, 20, 30, 50)
AXES = {'tm_score': 'score', 'chromatin_od': 'od51'}
MITOTIC, LOOKALIKE = 1, 2
BASE_SIZE, PATCH_SIZE = 51, 73
SELF_HIT_RADIUS = 5.0
DEEP_FLOOR_Z = -1.5
PEAK_MIN_DISTANCE = 7
OD_WINDOW, OD_PAD = 51, 26
OD_FRAC = 0.10          # chromatin.DEFAULT_FRAC -- the darkest 10 % of the window
MIDOG_RADIUS_UM = 7.5


def cells(nb_path):
    """All cells of a notebook as (index, cell_type, source, cell) tuples."""
    nb = json.load(open(nb_path))
    return [(i, c['cell_type'], ''.join(c['source']), c) for i, c in enumerate(nb['cells'])]


NBC = cells(NB)
SRC_ALL = '\n'.join(s for _, _, s, _ in NBC)
CODE_ALL = '\n'.join(s for _, t, s, _ in NBC if t == 'code')
MD_ALL = '\n'.join(s for _, t, s, _ in NBC if t == 'markdown')


# ---------------------------------------------------------------------------------------
# 1. Execution-coherence gate
# ---------------------------------------------------------------------------------------
def execution_gate():
    """Contiguity of execution_count, error outputs, unrun cells, embedded figures."""
    code = [(i, c) for i, t, _, c in NBC if t == 'code']
    ecs = [c.get('execution_count') for _, c in code]
    contiguous = ecs == list(range(1, len(ecs) + 1))
    n_err = sum(1 for _, c in code for o in c.get('outputs', []) if o.get('output_type') == 'error')
    n_unrun = sum(1 for e in ecs if e is None)
    n_png = sum(1 for _, c in code for o in c.get('outputs', [])
                if 'image/png' in o.get('data', {}))
    rows = [
        dict(check='n_cells', value=len(NBC), note='was 26 at round 1'),
        dict(check='n_code_cells', value=len(code), note='was 15 at round 1'),
        dict(check='execution_count_contiguous_from_1', value=bool(contiguous), note=str(ecs)),
        dict(check='n_error_outputs', value=n_err, note=''),
        dict(check='n_unrun_code_cells', value=n_unrun, note=''),
        dict(check='last_cell_run', value=bool(ecs[-1] is not None), note=''),
        dict(check='n_embedded_pngs', value=n_png, note='was 3 at round 1'),
    ]
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------
# 2. Correction-application matrix -- the 16 blocks + the dependency section
# ---------------------------------------------------------------------------------------
def cell_src(i):
    return NBC[i][2]


def cell_flat(i):
    """Cell source with comment markers and line wrapping removed, so a substring test on
    prose is not defeated by where the author happened to break the line."""
    lines = [ln.lstrip().lstrip('#').strip() for ln in cell_src(i).split('\n')]
    return ' '.join(' '.join(lines).split())


def correction_matrix():
    """One row per correction item: what was required, what the notebook now contains."""
    r = []

    def add(item, requirement, observed, verdict, evidence):
        r.append(dict(item=item, requirement=requirement, observed=observed,
                      verdict=verdict, evidence=evidence))

    # --- dependency section ---------------------------------------------------------
    max_peaks = re.search(r'^MAX_PEAKS\s*=\s*(\S+)', cell_src(1), re.M)
    max_peaks = max_peaks.group(1) if max_peaks else 'ABSENT'
    add('DEP-1 combined rerun (T1-3+T2-2+T3-3+T2-3 together, one rerun)',
        'all four pipeline changes land together and are rerun once over 14 ROIs',
        f'MAX_PEAKS={max_peaks}; arms={list(AXES)}; joint-seed helper present='
        f'{"joint" in CODE_ALL.lower()}; caps=(MAX_PEAKS,) present={"caps=(MAX_PEAKS,)" in CODE_ALL}',
        'NOT-APPLIED', 'cell 1 config + cell 5 run_roi')

    add('DEP-2 K=50 delivery risk under a 100-peak cap',
        'after the rerun, check budget_delivered==50 on every ROI and arm',
        'risk never incurred -- the cap was never changed; budget_delivered==50 on 112/112 '
        'rows, but at MAX_PEAKS=2_000_000, so this is no evidence the cap is safe',
        'NOT-APPLICABLE', 'raw CSV; MAX_PEAKS unchanged')

    verif = pd.read_csv(NEW_VERIF)
    annulus_viol = verif[(verif['check'] == 'seed_annulus_empty') &
                         (verif['label'] != 'ALL') & (~verif['passed'].astype(bool))]
    add('DEP-3 T1-1/T3-2 conditional on the cap',
        'if the 403 event vanishes, drop the mechanism paragraph and the rank fix',
        f'event still occurs on {len(annulus_viol)} ROI(s): '
        f'{list(annulus_viol["label"])} -- so the "still occurs" branch applies',
        'CONDITION-UNCHANGED', 'verification CSV')

    # --- T1-1 -----------------------------------------------------------------------
    c8 = cell_flat(8)
    c26 = cell_flat(26)
    c7 = cell_flat(7)
    old_mech = 'a nearby, higher-variance neighbour can outscore it'
    add('T1-1.1 cell 8 -- delete the 3-step mechanism, replace with the measured one',
        'replacement must name the 0.5px D8 half-pixel anchor / integer-peak gap as what '
        'opens the annulus (anchor geometry), not the self_hit/match radius gap',
        f'old mechanism removed={old_mech not in c8}; '
        f'mentions half-pixel/anchor geometry={bool(re.search(r"half.?pix|anchor", c8, re.I))}; '
        f'names self_hit/match-radius gap as "the real cause"='
        f'{"The real cause is geometric" in c8}',
        'APPLIED-INCORRECTLY', 'cell 8 comment')
    add('T1-1.2 cell 24->26 -- delete "31px template ... no longer correlates most strongly"',
        'delete outright, no softened replacement',
        f'phrase present={"no longer correlates most strongly" in c26}',
        'APPLIED-CORRECTLY', 'closing summary')
    add('T1-1.3 cell 24->26 -- delete "consistent with ... this run\'s one seed_annulus_empty '
        'violation"',
        'remove the causal linkage entirely',
        f'phrase present={"one seed_annulus_empty violation" in c26}',
        'APPLIED-CORRECTLY', 'closing summary')
    add('T1-1.4 cell 24->26 -- replace "a genuine, if rare, consequence of the channel change"',
        'replacement must attribute the event to D8 half-pixel template-anchor geometry',
        f'old phrase present={"genuine, if rare, consequence" in c26}; '
        f'replacement names the self-hit/match-radius gap instead of the anchor='
        f'{"gap between the fixed 5px self-hit radius" in c26}; '
        f'markdown cell 7 still calls the event "informative about what the channel change '
        f'does to template distinctiveness"={"informative about what the channel" in c7}',
        'APPLIED-INCORRECTLY', 'closing summary + markdown cell 7 (untouched)')

    # --- T1-2 -----------------------------------------------------------------------
    add('T1-2.1 rerun the exact sign-flip on the post-rerun jointly-gated dataset',
        'new p-values over however many ROIs got a jointly-valid seed (up to 14)',
        'sign-flip IS computed in new cell 21, but on the un-jointly-gated 11/3 split; '
        'p-values are therefore the pre-rerun ones recomputed on identical data',
        'APPLIED-IN-FALLBACK', 'cell 21; precondition T1-3 not met')
    add('T1-2.2 rewrite cell 24 tm_score/chromatin_od claims from the new test',
        'supported-claim sentence form; drop "consistently worse ... 4-6 pp"',
        f'"consistently worse" present={"consistently worse" in c26}; '
        f'"not resolvable at this design" present={"not resolvable at this design" in c26}; '
        f'sentence still quotes 11-ROI subset counts={"6-8 of 11 ROIs" in c26}',
        'APPLIED-IN-FALLBACK', 'closing summary')
    add('T1-2.3 note the single-click caveat in the closing summary',
        'joint gate removes the mismatch confound but does not replicate the estimate',
        f'single-seed caveat present={"Single-seed" in c26}; '
        f'small-G caveat present={"Small G" in c26}',
        'APPLIED-CORRECTLY', 'closing summary caveats')

    # --- T1-3 -----------------------------------------------------------------------
    has_joint = bool(re.search(r'joint', CODE_ALL, re.I))
    both_channels = ('gray_inverted' in CODE_ALL and
                     re.search(r'to_channel\([^)]*gray_inverted', CODE_ALL) is not None)
    add('T1-3.2 joint-validity seed draw (both channels accept the same candidate)',
        'notebook-local helper requiring tightened_template_box AND _patch_readable to '
        'succeed under BOTH hematoxylin_od and gray_inverted',
        f'any joint helper in code={has_joint}; gray_inverted ever converted in code='
        f'{both_channels}; _patch_readable referenced={"_patch_readable" in CODE_ALL}; '
        f'_check gates on `hem` only=True',
        'NOT-APPLIED', 'cell 5 `_check` closure')
    wlt_same_only = "DELTA_A[DELTA_A['same_seed']]" in cell_src(25)
    add('T1-3.3 rerun 14 ROIs, recompute Table A/B/deltas/win-counts/figures; drop the '
        'starred different-seed annotation from Figure 1',
        'no different-seed subset should remain',
        f'star annotation still in figure code={"different seed_ann_id" in cell_src(23)}; '
        f'same_seed column still built={"same_seed" in cell_src(19)}; '
        f'win/loss/tie still restricted to same-seed={wlt_same_only}',
        'NOT-APPLIED', 'cells 19/23/25')
    add('T1-3.4 report any ROI excluded because the joint pool is exhausted',
        'explicit exclusion reporting',
        'moot -- no joint gate exists, so no exclusion can arise',
        'NOT-APPLICABLE', '-')
    add('T1-3.5 rewrite the closing summary to one number per arm per budget over N ROIs; '
        'remove the "more severe confound" exclusion rationale',
        'exclusion rationale removed entirely',
        f'"more severe confound" still in markdown cell 22='
        f'{"the more severe confound" in cell_flat(22)}; closing summary now reports BOTH '
        f'the 11-ROI and the all-14 view={"Over all 14 ROIs" in c26}',
        'APPLIED-IN-FALLBACK', 'cells 22/26')

    # --- T1-4 -----------------------------------------------------------------------
    add('T1-4.1 delete "a smaller template is a less distinctive one ... consistent with '
        'the precision drop"',
        'delete, do not soften',
        f'old causal claim present={"less distinctive one for unnormalized" in c26}; '
        f'new text says the story "has no support in this run\'s own data"='
        f'{"has no support in" in c26}',
        'APPLIED-CORRECTLY', 'closing summary')
    add('T1-4.2 do not replace with a softened causal claim',
        'state only that no mechanism is established',
        f'"isn\'t established as the *reason*" present={"established as the *reason*" in c26}',
        'APPLIED-CORRECTLY', 'closing summary')

    # --- T2-1 (user declined) --------------------------------------------------------
    quotes_t21 = all(t in c26 for t in ['read_95', 'coverage_frac', 'pool size'])
    computed_t21 = any(t in CODE_ALL for t in ['read_95', 'coverage_frac', 'full_list_recall'])
    add('T2-1 (none -- user explicitly declined; recall/pool/reading-depth out of scope)',
        'nothing should be added',
        f'closing summary now quotes the declined T2-1 numbers={quotes_t21}; '
        f'any cell computes read_95/coverage_frac/full_list_recall deltas={computed_t21}',
        'ADDED-THOUGH-DECLINED', 'closing summary para 4; no code computes it')

    # --- T2-2 -----------------------------------------------------------------------
    add('T2-2.1 set MAX_PEAKS = 100 (D9)',
        'MAX_PEAKS = 100 in the config cell',
        f'MAX_PEAKS = {max_peaks}',
        'NOT-APPLIED', 'cell 1')
    add('T2-2.2 rerun both arms under the corrected cap',
        'recompute Table A/B/deltas/summary/figures at the cap',
        'notebook re-executed but at the superseded cap; every number identical to the '
        'pre-correction run',
        'NOT-APPLIED', 'cell 1 + raw CSV')
    add('T2-2.3 drop the "no coverage at the D9 production configuration" caveat',
        'only after the cap is applied',
        f'caveat retained and expanded in the closing summary='
        f'{"Known scope gaps this run does not resolve" in c26}',
        'CORRECTLY-RETAINED', 'closing summary (correct, since 2.1/2.2 not done)')

    # --- T2-3 -----------------------------------------------------------------------
    add('T2-3.1 add od_contrast as a third scored arm',
        'AXES carries three arms, rerun over 14 ROIs',
        f'AXES in config = {list(AXES)}; od_contrast appears only in markdown '
        f'({sum(1 for _, t, s, _ in NBC if t == "markdown" and "od_contrast" in s)} cells)',
        'NOT-APPLIED', 'cell 1')
    add('T2-3.2 recompute Table A/B/deltas/summary/figures with the third arm',
        '3 arms in every output table',
        f'arms in raw CSV = {sorted(pd.read_csv(NEW_RAW)["arm"].unique())}',
        'NOT-APPLIED', 'raw CSV')
    add('T2-3.3 report od_contrast delta beside its per-candidate compute cost',
        'measured margin vs one extra chromatin_density call at window 121',
        'no measurement; closing summary states the gap in prose only',
        'NOT-APPLIED', 'closing summary')

    # --- T2-4 -----------------------------------------------------------------------
    add('T2-4 (subsumed by T1-2; sentence shape "moves on X of N ROIs")',
        'fill the shape from the new data',
        f'"only 1 of 11 same-seed ROIs moves at all" present={"1 of 11 same-seed ROIs moves" in c26}; '
        f'"unresolvable at this design, not confirmed" present='
        f'{"unresolvable at this design, not confirmed" in c26}',
        'APPLIED-IN-FALLBACK', 'closing summary')

    # --- T2-5 -----------------------------------------------------------------------
    c9 = cell_flat(9)
    add('T2-5 markdown cell 9 -- remove "doubles as a fidelity check on the channel change"',
        'remove the claim and its justification, no softened replacement',
        f'claim now explicitly retracted={"It does not" in c9}; '
        f'the phrase survives only inside the retraction='
        f'{c9.count("doubles as a fidelity check") == 1}',
        'APPLIED-CORRECTLY', 'markdown cell 9')

    # --- T3-1 -----------------------------------------------------------------------
    add('T3-1 "2-3 wins" -> "1-3 wins" (fallback if the rerun is skipped)',
        'fallback text edit',
        f'"1-3 wins vs. 6-8 losses" present={"1-3 wins vs. 6-8 losses" in c26}; '
        f'"2-3 wins" present={"2-3 wins" in c26}',
        'APPLIED-CORRECTLY', 'closing summary')

    # --- T3-2 -----------------------------------------------------------------------
    add('T3-2 rank off-by-one (0-based rank reported as an ordinal)',
        'recompute the rank fresh under the new cap; report as ordinal',
        f'cell 8 now states both 0-indexed and ordinal='
        f'{"0-indexed rank is 3281" in c8 and "3282nd and 394th by ordinal count" in c8}',
        'APPLIED-CORRECTLY', 'cell 8 (numbers re-verified from pixels below)')

    # --- T3-3 -----------------------------------------------------------------------
    verif_caps = pd.read_csv(NEW_VERIF)
    n_no_cap = int((verif_caps['check'] == 'no_cap').sum())
    add('T3-3 drop MAX_PEAKS from Arm(caps=(MAX_PEAKS,)) in the verification pass',
        'must land together with T2-2',
        f'caps=(MAX_PEAKS,) still in cell 5={"caps=(MAX_PEAKS,)" in cell_src(5)}; '
        f'no_cap rows still in verification CSV={n_no_cap}',
        'NOT-APPLIED', 'cell 5 + verification CSV')

    # --- T3-4 -----------------------------------------------------------------------
    add('T3-4 add a cell printing raw[nan_rate].max() and raw[largest_tie_block].max()',
        'one new cell after the raw table',
        f'"nan_rate" anywhere in the notebook={"nan_rate" in SRC_ALL}; '
        f'"largest_tie_block" anywhere={"largest_tie_block" in SRC_ALL}',
        'NOT-APPLIED', 'whole notebook')

    # --- T3-5 -----------------------------------------------------------------------
    c25 = cell_src(25)
    ties_stacked = "for part in ('wins', 'losses', 'ties')" in c25
    add('T3-5 add ties as a visible third segment in the win-count bars',
        'bar height reflects all 11 ROIs',
        f'ties stacked={ties_stacked}; '
        f'ylim set to n+1={"set_ylim(0, len(SAME_SEED_ROIS) + 1)" in c25}',
        'APPLIED-CORRECTLY', 'cell 25')

    # --- T3-6 -----------------------------------------------------------------------
    c23 = cell_src(23)
    add('T3-6 use one shared colour scale across both heatmap panels',
        'same vmin/vmax on both panels instead of a per-panel lim',
        f'per-panel lim still computed inside the loop='
        f'{"lim = max(float(np.abs(mat).max()), 1e-6)" in c23}; '
        f'markdown caveat substituted instead={"Caveat (T3-6, audit)" in cell_flat(22)}',
        'APPLIED-INCORRECTLY', 'cells 22/23')

    # --- T3-7 -----------------------------------------------------------------------
    c14 = cell_src(14)
    tb = pd.read_csv(NEW_BY_DOMAIN)
    add('T3-7 flag or join tied worst_roi_file instead of silently keeping the first row',
        'joined list of tied file_names, or a "tied" marker',
        f'worst_roi_tied column added={"worst_roi_tied" in c14}; '
        f'worst_roi_file still the arbitrary single pick={"worst_roi_file=str(worst" in c14}; '
        f'rows flagged tied in the by_domain CSV='
        f'{int(tb["worst_roi_tied"].sum()) if "worst_roi_tied" in tb.columns else "COLUMN ABSENT"}',
        'APPLIED-CORRECTLY', 'cell 14 + by_domain CSV')

    return pd.DataFrame(r)


# ---------------------------------------------------------------------------------------
# 3. Tier A -- statistics re-derived from the raw CSVs
# ---------------------------------------------------------------------------------------
def exact_sign_flip_p(deltas):
    """Two-sided exact sign-flip p over the nonzero deltas, and the effective unit count."""
    d = np.asarray(deltas, dtype=float)
    nz = d[d != 0]
    if len(nz) == 0:
        return 1.0, 0
    obs = abs(d.sum())
    a = np.abs(nz)
    tot = cnt = 0
    for signs in itertools.product([1, -1], repeat=len(a)):
        tot += 1
        if abs(float(np.dot(signs, a))) >= obs - 1e-9:
            cnt += 1
    return cnt / tot, len(nz)


def load_frames():
    """New raw/per-ROI/baseline tables plus the same-seed partition, re-derived here."""
    new_raw = pd.read_csv(NEW_RAW)
    old_raw = pd.read_csv(OLD_RAW_P)
    new_a = pd.read_csv(NEW_PER_ROI)
    old_a = pd.read_csv(OLD_PER_ROI)
    old_a = old_a[old_a['arm'].isin(AXES)].reset_index(drop=True)
    new_seed = new_raw.drop_duplicates('file_name').set_index('file_name')[['seed_ann_id', 'base_size']]
    old_seed = old_raw.drop_duplicates('file_name').set_index('file_name')[['seed_ann_id', 'base_size']]
    same = (new_seed['seed_ann_id'] == old_seed['seed_ann_id'])
    return new_raw, old_raw, new_a, old_a, new_seed, old_seed, same


def tier_a_inference():
    """Sign-flip / win-loss-tie / t-CI over ROI deltas, same-seed-11 and all-14."""
    new_raw, old_raw, new_a, old_a, new_seed, old_seed, same = load_frames()
    m = new_a.merge(old_a, on=['file_name', 'arm'], suffixes=('_new', '_old'))
    rows = []
    for arm in AXES:
        for k in BUDGETS:
            m[f'd{k}'] = m[f'precision_at_{k}_new'] - m[f'precision_at_{k}_old']
        sub = m[m['arm'] == arm].set_index('file_name')
        for label, idx in (('same_seed_11', same[same].index), ('all_14', same.index)):
            for k in BUDGETS:
                d = sub.loc[idx, f'd{k}'].to_numpy()
                p, keff = exact_sign_flip_p(d)
                w, l, t = int((d > 0).sum()), int((d < 0).sum()), int((d == 0).sum())
                n = len(d)
                se = d.std(ddof=1) / np.sqrt(n)
                tcrit = sps.t.ppf(0.975, n - 1)
                rows.append(dict(subset=label, arm=arm, K=k, n_roi=n,
                                 mean_delta=round(float(d.mean()), 4),
                                 wins=w, losses=l, ties=t,
                                 n_effective_units=keff,
                                 p_floor=round(min(1.0, 2 / 2 ** keff), 4) if keff else 1.0,
                                 p_exact_sign_flip=round(p, 4),
                                 ci95_lo=round(float(d.mean() - tcrit * se), 4),
                                 ci95_hi=round(float(d.mean() + tcrit * se), 4)))
    return pd.DataFrame(rows)


def tie_zero_check():
    """Cell 20 claims a same-tp tie always subtracts to exactly 0.0. Test it per (arm, K)."""
    new_raw, old_raw, new_a, old_a, *_ = load_frames()
    m = new_a.merge(old_a, on=['file_name', 'arm'], suffixes=('_new', '_old'))
    rows = []
    for arm in AXES:
        for k in BUDGETS:
            sub = m[m['arm'] == arm]
            d = (sub[f'precision_at_{k}_new'] - sub[f'precision_at_{k}_old'])
            d.index = sub['file_name']
            ntp = new_raw[(new_raw['arm'] == arm) & (new_raw['budget'] == k)].set_index('file_name')['tp_at_budget']
            otp = old_raw[(old_raw['arm'] == arm) & (old_raw['budget'] == k)].set_index('file_name')['tp_at_budget']
            eq = (ntp - otp == 0)
            eq_tp = int(eq.sum())
            rows.append(dict(arm=arm, K=k, n_delta_exactly_zero=int((d == 0).sum()),
                             n_equal_tp=eq_tp,
                             agrees=bool(int((d == 0).sum()) == eq_tp),
                             max_abs_residual_on_equal_tp=float(
                                 np.max(np.abs(d.loc[eq[eq].index].to_numpy()))
                                 if eq_tp else 0.0)))
    return pd.DataFrame(rows)


def spearman_four_ways():
    """Spearman(delta_base_size, delta_precision) under both delta definitions,
    both p-value methods, and both ROI subsets -- to adjudicate cell 20's claim."""
    new_raw, old_raw, new_a, old_a, new_seed, old_seed, same = load_frames()
    dbs = (new_seed['base_size'] - old_seed['base_size']).rename('dbs')
    m = new_a.merge(old_a, on=['file_name', 'arm'], suffixes=('_new', '_old'))
    rng = np.random.default_rng(0)
    rows = []
    for arm in AXES:
        for k in BUDGETS:
            sub = m[m['arm'] == arm].set_index('file_name')
            d_round = (sub[f'precision_at_{k}_new'] - sub[f'precision_at_{k}_old'])
            ntp = new_raw[(new_raw['arm'] == arm) & (new_raw['budget'] == k)].set_index('file_name')['tp_at_budget']
            otp = old_raw[(old_raw['arm'] == arm) & (old_raw['budget'] == k)].set_index('file_name')['tp_at_budget']
            d_exact = (ntp - otp) / k
            for dname, dser in (('rounded_precision_cols', d_round), ('integer_tp_diff', d_exact)):
                for sname, idx in (('all_14', same.index), ('same_seed_11', same[same].index)):
                    x = dbs.loc[idx].to_numpy(dtype=float)
                    y = dser.loc[idx].to_numpy(dtype=float)
                    rho, p_asym = sps.spearmanr(x, y)
                    obs = abs(rho)
                    cnt = 0
                    for _ in range(2000):
                        r2, _ = sps.spearmanr(x, rng.permutation(y))
                        cnt += abs(r2) >= obs - 1e-12
                    rows.append(dict(delta_def=dname, subset=sname, arm=arm, K=k, n=len(x),
                                     rho=round(float(rho), 3),
                                     p_asymptotic=round(float(p_asym), 3),
                                     p_permutation_2000=round(cnt / 2000, 3)))
    return pd.DataFrame(rows)



def rounding_artifact():
    """Cell 20 blames IEEE float subtraction for the round-1 Spearman divergence. Test the
    competing explanation: Table A stores `round(tp/K, 4)`, so two ROIs with the SAME integer
    tp difference can land on different rounded deltas, which is what splits Spearman's ranks."""
    new_raw, old_raw, new_a, old_a, *_ = load_frames()
    m = new_a.merge(old_a, on=['file_name', 'arm'], suffixes=('_new', '_old'))
    rows = []
    for arm in AXES:
        for k in BUDGETS:
            sub = m[m['arm'] == arm].set_index('file_name')
            d_round = (sub[f'precision_at_{k}_new'] - sub[f'precision_at_{k}_old'])
            ntp = new_raw[(new_raw['arm'] == arm) & (new_raw['budget'] == k)].set_index('file_name')['tp_at_budget']
            otp = old_raw[(old_raw['arm'] == arm) & (old_raw['budget'] == k)].set_index('file_name')['tp_at_budget']
            dtp = (ntp - otp).reindex(d_round.index)
            g = pd.DataFrame({'dtp': dtp, 'bitexact': d_round, 'q4': d_round.round(9)})
            split_bits = g.groupby('dtp')['bitexact'].nunique()
            split_q4 = g.groupby('dtp')['q4'].nunique()
            spreads = g.groupby('dtp')['bitexact'].agg(lambda v: float(v.max() - v.min()))
            worst = spreads.abs().idxmax()
            rows.append(dict(arm=arm, K=k,
                             n_distinct_tp_diffs=int(g['dtp'].nunique()),
                             tp_diffs_split_by_ieee_noise=int((split_bits > 1).sum()),
                             tp_diffs_split_by_4dp_quantisation=int((split_q4 > 1).sum()),
                             worst_tp_diff=int(worst),
                             max_spread_within_one_tp_diff=float(spreads.abs().max()),
                             example=str(sorted(set(
                                 repr(v) for v in g[g['dtp'] == worst]['bitexact'])))[:120]))
    return pd.DataFrame(rows)


def stat_context_check():
    """Every cell of the notebook's own stat_context CSV, rebuilt from the raw tables."""
    nb_stat = pd.read_csv(NEW_STAT)
    mine = tier_a_inference()
    sp = spearman_four_ways()
    sp = sp[(sp['delta_def'] == 'integer_tp_diff') & (sp['subset'] == 'all_14')]
    rows = []
    for _, o in nb_stat.iterrows():
        arm, k = o['arm'], int(o['K'])
        s11 = mine[(mine['subset'] == 'same_seed_11') & (mine['arm'] == arm) & (mine['K'] == k)].iloc[0]
        a14 = mine[(mine['subset'] == 'all_14') & (mine['arm'] == arm) & (mine['K'] == k)].iloc[0]
        spr = sp[(sp['arm'] == arm) & (sp['K'] == k)].iloc[0]
        checks = [
            ('mean_delta_same_seed_11', o['mean_delta_same_seed_11'], s11['mean_delta']),
            ('wlt_same_seed_11', o['wlt_same_seed_11'],
             f"{s11['wins']}/{s11['losses']}/{s11['ties']}"),
            ('exact_p_same_seed_11', o['exact_p_same_seed_11'], s11['p_exact_sign_flip']),
            ('mean_delta_all_14', o['mean_delta_all_14'], a14['mean_delta']),
            ('wlt_all_14', o['wlt_all_14'], f"{a14['wins']}/{a14['losses']}/{a14['ties']}"),
            ('exact_p_all_14', o['exact_p_all_14'], a14['p_exact_sign_flip']),
            ('spearman_rho_vs_delta_base_size', o['spearman_rho_vs_delta_base_size'], spr['rho']),
            ('spearman_p', o['spearman_p'], spr['p_asymptotic']),
        ]
        for name, got, exp in checks:
            ok = (str(got) == str(exp)) if isinstance(exp, str) else \
                bool(np.isclose(float(got), float(exp), atol=5e-4))
            rows.append(dict(arm=arm, K=k, column=name, notebook=got, audit=exp, agrees=ok))
    return pd.DataFrame(rows)


def recall_pool_context():
    """The declined-T2-1 numbers cell 26 now quotes: recompute them from the raw CSVs."""
    new_raw, old_raw, *_ = load_frames()
    rows = []
    specs = [('n_detections', 'tm_score'), ('coverage_frac', 'tm_score'),
             ('full_list_recall', 'tm_score'), ('read_95', 'tm_score'),
             ('read_95', 'chromatin_od')]
    for col, arm in specs:
        n = new_raw[(new_raw['arm'] == arm) & (new_raw['budget'] == 10)].set_index('file_name')[col]
        o = old_raw[(old_raw['arm'] == arm) & (old_raw['budget'] == 10)].set_index('file_name')[col]
        d = (n - o).fillna(0.0)
        p, keff = exact_sign_flip_p(d.to_numpy())
        rows.append(dict(column=col, arm=arm, mean_new=round(float(n.mean()), 4),
                         mean_old=round(float(o.mean()), 4),
                         mean_delta=round(float(d.mean()), 4),
                         n_better=int((d > 0).sum()), n_worse=int((d < 0).sum()),
                         n_effective_units=keff, p_exact_sign_flip=round(p, 4)))
    n = new_raw[(new_raw['arm'] == 'tm_score') & (new_raw['budget'] == 10)].set_index('file_name')['full_list_recall']
    o = old_raw[(old_raw['arm'] == 'tm_score') & (old_raw['budget'] == 10)].set_index('file_name')['full_list_recall']
    reach1 = sorted(n[(n >= 1.0) & (o < 1.0)].index)
    rows.append(dict(column='full_list_recall_reaching_1.000', arm='tm_score',
                     mean_new=np.nan, mean_old=np.nan, mean_delta=np.nan,
                     n_better=len(reach1), n_worse=0, n_effective_units=len(reach1),
                     p_exact_sign_flip=np.nan))
    df = pd.DataFrame(rows)
    df.attrs['reach1'] = reach1
    return df, reach1


def prose_numbers(inf, recall, reach1, spear):
    """Every number the corrected notebook's prose asserts, against my own recomputation."""
    new_raw, old_raw, new_a, old_a, new_seed, old_seed, same = load_frames()
    c26 = cell_flat(26)
    c8 = cell_flat(8)

    def g(subset, arm, k, col):
        r = inf[(inf['subset'] == subset) & (inf['arm'] == arm) & (inf['K'] == k)]
        return r.iloc[0][col]

    rows = []

    def add(claim, says, finds, ok, where):
        rows.append(dict(claim=claim, notebook_says=says, audit_finds=finds,
                         agrees=ok, where=where))

    p20 = g('same_seed_11', 'tm_score', 20, 'p_exact_sign_flip')
    add('tm_score resolvable at K=20 only, p=0.039', '0.039', round(float(p20), 4),
        bool(np.isclose(p20, 0.039, atol=6e-4)), 'cell 26')
    others = {k: float(g('same_seed_11', 'tm_score', k, 'p_exact_sign_flip'))
              for k in (10, 30, 50)}
    add('K=10/30/50 do not clear p<0.05', 'none clear 0.05', others,
        all(v >= 0.05 for v in others.values()), 'cell 26')
    wl = {k: (int(g('same_seed_11', 'tm_score', k, 'wins')),
              int(g('same_seed_11', 'tm_score', k, 'losses'))) for k in BUDGETS}
    add('1-3 wins vs 6-8 losses for tm_score across the four budgets', '1-3 / 6-8', wl,
        min(w for w, _ in wl.values()) == 1 and max(w for w, _ in wl.values()) == 3
        and min(l for _, l in wl.values()) == 6 and max(l for _, l in wl.values()) == 8,
        'cell 26')
    add('tm_score falls on 6-8 of 11 ROIs at every budget', '6-8 of 11',
        {k: int(g('same_seed_11', 'tm_score', k, 'losses')) for k in BUDGETS},
        all(6 <= int(g('same_seed_11', 'tm_score', k, 'losses')) <= 8 for k in BUDGETS),
        'cell 26')
    m10 = float(g('all_14', 'tm_score', 10, 'mean_delta'))
    add('all-14 K=10 mean delta +0.021, 6 better / 6 worse / 2 tied, p=0.75',
        '+0.021, 6/6/2, p=0.75',
        f"{m10:+.4f}, {int(g('all_14','tm_score',10,'wins'))}/"
        f"{int(g('all_14','tm_score',10,'losses'))}/{int(g('all_14','tm_score',10,'ties'))}, "
        f"p={float(g('all_14','tm_score',10,'p_exact_sign_flip')):.4f}",
        bool(np.isclose(m10, 0.021, atol=6e-4)
             and (int(g('all_14', 'tm_score', 10, 'wins')),
                  int(g('all_14', 'tm_score', 10, 'losses')),
                  int(g('all_14', 'tm_score', 10, 'ties'))) == (6, 6, 2)
             and np.isclose(float(g('all_14', 'tm_score', 10, 'p_exact_sign_flip')), 0.75, atol=6e-3)),
        'cell 26')
    pen = [abs(float(g('all_14', 'tm_score', k, 'mean_delta'))) for k in (20, 30, 50)]
    add('over all 14 ROIs the tm_score penalty shrinks to roughly 1-2pp',
        '1-2 pp', [round(100 * p, 2) for p in pen],
        all(0.005 <= p <= 0.025 for p in pen), 'cell 26')
    add('chromatin_od: only 1 of 11 same-seed ROIs moves at K=10', '1 of 11',
        11 - int(g('same_seed_11', 'chromatin_od', 10, 'ties')),
        11 - int(g('same_seed_11', 'chromatin_od', 10, 'ties')) == 1, 'cell 26')

    # declined-T2-1 numbers
    rr = recall.set_index(['column', 'arm'])
    add('pool size +399 candidates/ROI, 10/14 larger, p=0.007', '+399, 10/14, p=0.007',
        f"{rr.loc[('n_detections','tm_score'),'mean_delta']:+.1f}, "
        f"{int(rr.loc[('n_detections','tm_score'),'n_better'])}/14, "
        f"p={rr.loc[('n_detections','tm_score'),'p_exact_sign_flip']:.4f}",
        bool(abs(rr.loc[('n_detections', 'tm_score'), 'mean_delta'] - 399) < 1.5
             and int(rr.loc[('n_detections', 'tm_score'), 'n_better']) == 10
             and abs(rr.loc[('n_detections', 'tm_score'), 'p_exact_sign_flip'] - 0.007) < 6e-4),
        'cell 26 (declined T2-1)')
    add('coverage_frac +0.007, p=0.006', '+0.007, p=0.006',
        f"{rr.loc[('coverage_frac','tm_score'),'mean_delta']:+.4f}, "
        f"p={rr.loc[('coverage_frac','tm_score'),'p_exact_sign_flip']:.4f}",
        bool(abs(rr.loc[('coverage_frac', 'tm_score'), 'mean_delta'] - 0.007) < 6e-4
             and abs(rr.loc[('coverage_frac', 'tm_score'), 'p_exact_sign_flip'] - 0.006) < 6e-4),
        'cell 26 (declined T2-1)')
    add('read_95 for tm_score falls by 755 candidates, p=0.027', '-755, p=0.027',
        f"{rr.loc[('read_95','tm_score'),'mean_delta']:+.1f}, "
        f"p={rr.loc[('read_95','tm_score'),'p_exact_sign_flip']:.4f}",
        bool(abs(rr.loc[('read_95', 'tm_score'), 'mean_delta'] + 755) < 1.5
             and abs(rr.loc[('read_95', 'tm_score'), 'p_exact_sign_flip'] - 0.027) < 6e-4),
        'cell 26 (declined T2-1)')
    add('300.tiff and 301.tiff reach full-list recall 1.000', '300, 301', reach1,
        sorted(reach1) == ['300.tiff', '301.tiff'], 'cell 26 (declined T2-1)')

    # mechanism / Spearman claims
    sp_all = spear[(spear['delta_def'] == 'integer_tp_diff') & (spear['subset'] == 'all_14')]
    sp_11 = spear[(spear['delta_def'] == 'integer_tp_diff') & (spear['subset'] == 'same_seed_11')]
    worst11 = sp_11.loc[sp_11['p_asymptotic'].idxmin()]
    add('"every [Spearman] coefficient is small and non-significant"',
        'none significant; all small',
        f"computed subset (all_14): max |rho|={sp_all['rho'].abs().max():.3f}, min p_asym="
        f"{sp_all['p_asymptotic'].min():.3f}. UNCOMPUTED subset (same_seed_11, 8 of the "
        f"audit's 16 cells): {worst11['arm']}/K={worst11['K']} rho={worst11['rho']}, "
        f"p_asym={worst11['p_asymptotic']:.3f}, p_perm={worst11['p_permutation_2000']:.3f}",
        bool(sp_all['rho'].abs().max() < 0.3 and sp_11['p_permutation_2000'].min() >= 0.05),
        'cell 26')
    add('"on the all-14 set the sign runs the wrong way"',
        'negative rho on all-14 tm_score',
        dict(zip(sp_all[sp_all['arm'] == 'tm_score']['K'],
                 sp_all[sp_all['arm'] == 'tm_score']['rho'])),
        bool((sp_all[sp_all['arm'] == 'tm_score']['rho'] <= 0).all()), 'cell 26')
    d245 = int(new_seed.loc['245.tiff', 'base_size']), int(old_seed.loc['245.tiff', 'base_size'])
    add('245.tiff shrank 47 -> 23 px', '47 -> 23', f'{d245[1]} -> {d245[0]}',
        d245 == (23, 47), 'cell 26')
    dbs = (new_seed['base_size'] - old_seed['base_size'])
    add('base_size shrinks on 13 of 14 ROIs', '13 shrink, 1 grows',
        f'{int((dbs < 0).sum())} shrink, {int((dbs > 0).sum())} grow, '
        f'{int((dbs == 0).sum())} unchanged',
        bool(int((dbs < 0).sum()) == 13 and int((dbs > 0).sum()) == 1), 'cells 17/26')
    add('seed_ann_id matches the baseline on 11/14 ROIs', '11/14',
        f'{int(same.sum())}/14', int(same.sum()) == 11, 'cell 18')
    add('MAX_PEAKS = 2_000_000, D9 supersedes it with 100',
        'notebook still at 2_000_000',
        re.search(r'^MAX_PEAKS\s*=\s*(\S+)', cell_src(1), re.M).group(1),
        True, 'cell 26 scope-gap paragraph (self-reported, and true)')
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------
# 4. Tier B -- recomputed from ROI pixels, no midog_utils in the path
# ---------------------------------------------------------------------------------------
def load_db():
    """Images and annotations from MIDOG++.json, read straight from the JSON (no midog_utils).
    bbox is TLBR, so the click is ((x1+x2)/2, (y1+y2)/2); votes come from `labels`."""
    d = json.load(open(DB))
    images = pd.DataFrame([dict(image_id=im['id'], file_name=im['file_name'],
                                tumor_type=im['tumor_type']) for im in d['images']])
    rows = []
    for a in d['annotations']:
        x1, y1, x2, y2 = a['bbox']
        votes = a.get('labels', [])
        rows.append(dict(ann_id=a['id'], image_id=a['image_id'],
                         cx=(x1 + x2) / 2.0, cy=(y1 + y2) / 2.0,
                         category_id=a['category_id'], n_votes=len(votes),
                         n_mitotic_votes=sum(1 for v in votes if v == MITOTIC)))
    return images, pd.DataFrame(rows), None


def my_hematoxylin_od(rgb):
    """Unclipped hematoxylin optical density, float32 -- the `hematoxylin_od` channel."""
    return rgb2hed(rgb.astype(np.float32) / 255.0)[:, :, 0].astype(np.float32)


def my_read_patch(img, cx, cy, size):
    """Square patch centred on the rounded (cx, cy), or None at the border."""
    half = size // 2
    ix, iy = int(round(cx)), int(round(cy))
    h, w = img.shape[:2]
    if ix - half < 0 or iy - half < 0 or ix + half >= w or iy + half >= h:
        return None
    return np.ascontiguousarray(img[iy - half:iy + half + 1, ix - half:ix + half + 1])


def my_tighten(patch):
    """Otsu gate + centre-containment + the three sanity gates. Returns (y0,y1,x0,x1)."""
    u8 = cv2.normalize(patch.astype(np.float32), None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    _, binary = cv2.threshold(u8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    labels = label(binary, connectivity=2)
    cy, cx = patch.shape[0] // 2, patch.shape[1] // 2
    cl = labels[cy, cx]
    if cl == 0:
        return None
    region = next(p for p in regionprops(labels) if p.label == cl)
    y0, x0, y1, x1 = region.bbox
    if not (y0 <= cy < y1 and x0 <= cx < x1):
        return None
    if region.area < 50 or region.area > 0.85 * patch.size or region.solidity < 0.5:
        return None
    return int(y0), int(y1), int(x0), int(x1)


def my_tightened_template_box(chan, cx, cy, otsu_window=BASE_SIZE):
    """D8 constructor: odd base_size and the bbox pixel centre (x0+x1-1)/2."""
    patch = my_read_patch(chan, cx, cy, otsu_window)
    if patch is None:
        return None
    bb = my_tighten(patch)
    if bb is None:
        return None
    y0, y1, x0, x1 = bb
    n = int(round(max(y1 - y0, x1 - x0)))
    if n % 2 == 0:
        n += 1
    base = max(5, n)
    half = otsu_window // 2
    ix, iy = int(round(cx)), int(round(cy))
    return base, ix - half + (x0 + x1 - 1) / 2.0, iy - half + (y0 + y1 - 1) / 2.0


def my_nms(centers, scores, radius):
    """Greedy NMS in descending score order, stable on ties."""
    tree = KDTree(centers)
    nb = tree.query_radius(centers, r=radius)
    order = np.argsort(-scores, kind='stable')
    sup = np.zeros(len(centers), dtype=bool)
    keep = []
    for i in order:
        if sup[i]:
            continue
        keep.append(i)
        sup[nb[i]] = True
    return np.asarray(keep, dtype=int)


def my_greedy_match(det_xy, gt_xy, radius):
    """One-to-one greedy match, best-first; returns each detection's claimed GT index."""
    det_to_gt = np.full(len(det_xy), -1, dtype=int)
    gt_to_det = np.full(len(gt_xy), -1, dtype=int)
    if not len(det_xy) or not len(gt_xy):
        return det_to_gt
    tree = KDTree(gt_xy)
    nb = tree.query_radius(det_xy, r=radius)
    for i, cands in enumerate(nb):
        free = [g for g in cands if gt_to_det[g] == -1]
        if not free:
            continue
        d = np.hypot(gt_xy[free, 0] - det_xy[i, 0], gt_xy[free, 1] - det_xy[i, 1])
        g = free[int(np.argmin(d))]
        det_to_gt[i] = g
        gt_to_det[g] = i
    return det_to_gt


def my_chromatin_density(chan, cx, cy, window=OD_WINDOW, frac=OD_FRAC):
    """Mean of the darkest `frac` of pixels in a window-sized box."""
    p = my_read_patch(chan, cx, cy, window)
    if p is None:
        return float('nan')
    f = p.ravel()
    k = max(1, int(frac * f.size))
    return float(np.partition(f, -k)[-k:].mean())


def roi_mpp(path):
    """Microns per pixel from the TIFF's own resolution tags (unit 2 = inch, 3 = cm)."""
    import tifffile
    with tifffile.TiffFile(str(path)) as tf:
        page = tf.pages[0]
        num, den = page.tags['XResolution'].value
        unit = int(page.tags['ResolutionUnit'].value)
    return (10000.0 if unit == 3 else 25400.0) / (num / den)


def load_rgb(path):
    """Full-resolution RGB uint8 array for one ROI, alpha stripped."""
    import tifffile
    with tifffile.TiffFile(str(path)) as tf:
        series = tf.series[0]
        levels = getattr(series, 'levels', None)
        arr = levels[0].asarray() if levels else series.asarray()
    return np.ascontiguousarray(arr[:, :, :3])


def tier_b(files, trace_403=True):
    """Full from-pixels rerun of the notebook's pipeline for `files`. No midog_utils."""
    images, ann, _ = load_db()
    meta = images.set_index('file_name')[['image_id', 'tumor_type']]
    rows, trace = [], []
    for fn in files:
        t0 = time.time()
        path = os.path.join(IMAGES_DIR, fn)
        image_id = int(meta.loc[fn, 'image_id'])
        rgb = load_rgb(path)
        mpp = roi_mpp(path)
        H, W = rgb.shape[:2]
        match_radius = MIDOG_RADIUS_UM / mpp
        nms_radius = match_radius
        hem = my_hematoxylin_od(rgb)
        del rgb

        gt = ann[ann['image_id'] == image_id].reset_index(drop=True)
        mit = gt[gt['category_id'] == MITOTIC]
        unan = mit[mit['n_mitotic_votes'] == mit['n_votes']]
        pool = unan if len(unan) else mit[(mit['n_votes'] > 0) &
                                          (mit['n_mitotic_votes'] / mit['n_votes'] >= 2 / 3) &
                                          (mit['n_mitotic_votes'] < mit['n_votes'])]
        border = PATCH_SIZE // 2
        ix = np.rint(pool['cx'].to_numpy()).astype(int)
        iy = np.rint(pool['cy'].to_numpy()).astype(int)
        pool = pool[(ix >= border) & (ix <= W - 1 - border) &
                    (iy >= border) & (iy <= H - 1 - border)]

        rng = np.random.default_rng([0, image_id])
        working, retries = pool.copy(), 0
        seed = spec = None
        while len(working):
            i = int(rng.integers(len(working)))
            row = working.iloc[i]
            got = my_tightened_template_box(hem, float(row['cx']), float(row['cy']))
            ok = got is not None and my_read_patch(hem, got[1], got[2], PATCH_SIZE) is not None
            if ok:
                seed, spec = row, got
                break
            working = working.drop(working.index[i])
            retries += 1
        base_size, tx, ty = spec
        seed_ann = int(seed['ann_id'])
        gt_eval = gt[gt['ann_id'] != seed_ann].reset_index(drop=True)
        n_gt = int((gt_eval['category_id'] == MITOTIC).sum())

        patch = my_read_patch(hem, tx, ty, PATCH_SIZE)
        c = patch.shape[0] // 2
        hb = base_size // 2
        tmpl = np.ascontiguousarray(patch[c - hb:c + hb + 1, c - hb:c + hb + 1], dtype=np.float32)
        PAD = (tmpl.shape[0] - 1) // 2
        hp = cv2.copyMakeBorder(hem, PAD, PAD, PAD, PAD, borderType=cv2.BORDER_REPLICATE)
        res = np.asarray(cv2.matchTemplate(np.ascontiguousarray(hp, dtype=np.float32),
                                           tmpl, cv2.TM_CCOEFF), dtype=np.float32)
        del hp
        fused = res[0:H, 0:W]   # matchTemplate on the PAD-replicated ROI is exactly HxW
        sample = fused[::8, ::8].ravel()
        sample = sample[np.isfinite(sample)]
        med = float(np.median(sample))
        mad = float(1.4826 * np.median(np.abs(sample - med)))
        cut = med + DEEP_FLOOR_Z * mad
        k = 2 * PEAK_MIN_DISTANCE + 1
        dil = cv2.dilate(fused, np.ones((k, k), np.uint8))
        mask = (fused >= dil) & (fused >= cut)
        ys, xs = np.nonzero(mask)
        sc = fused[ys, xs]
        order = np.lexsort((ys, xs, -sc))
        centers = np.stack([xs[order], ys[order]], axis=1).astype(np.float64)
        scores = sc[order]
        n_pre_nms = len(centers)

        if trace_403 and fn == '403.tiff':
            gmax = float(fused.max())
            gy, gx = np.unravel_index(int(np.argmax(fused)), fused.shape)
            d = np.hypot(centers[:, 0] - tx, centers[:, 1] - ty)
            j = int(np.argmin(d))
            trace.append(dict(file_name=fn, base_size=base_size, n_pre_nms_peaks=n_pre_nms,
                              global_max_score=round(gmax, 3),
                              global_max_xy=f'({gx},{gy})',
                              global_max_dist_to_anchor=round(
                                  float(np.hypot(gx - tx, gy - ty)), 3),
                              template_anchor=f'({tx},{ty})',
                              anchor_x_is_half_integer=bool(abs(tx % 1 - 0.5) < 1e-9),
                              anchor_y_is_half_integer=bool(abs(ty % 1 - 0.5) < 1e-9),
                              nearest_peak_xy=f'({centers[j,0]:.0f},{centers[j,1]:.0f})',
                              nearest_peak_dist=round(float(d[j]), 3),
                              nearest_peak_score=round(float(scores[j]), 3),
                              match_radius_px=round(match_radius, 3)))

        keep = my_nms(centers, scores, nms_radius)
        cc, ss = centers[keep], scores[keep]
        if trace_403 and fn == '403.tiff':
            dd = np.hypot(cc[:, 0] - tx, cc[:, 1] - ty)
            jj = int(np.argmin(dd))
            trace[-1].update(self_peak_survived_nms=bool(dd[jj] <= 1.0),
                             self_peak_score_post_nms=round(float(ss[jj]), 3),
                             self_peak_inside_self_hit_radius=bool(dd[jj] <= SELF_HIT_RADIUS))
        ok = np.hypot(cc[:, 0] - tx, cc[:, 1] - ty) > SELF_HIT_RADIUS
        cc, ss = cc[ok], ss[ok]
        pooldf = pd.DataFrame({'cx': cc[:, 0], 'cy': cc[:, 1], 'score': ss})

        hem_pad = cv2.copyMakeBorder(hem, OD_PAD, OD_PAD, OD_PAD, OD_PAD, cv2.BORDER_REPLICATE)
        pooldf['od51'] = [my_chromatin_density(hem_pad, x + OD_PAD, y + OD_PAD)
                          for x, y in zip(pooldf['cx'], pooldf['cy'])]
        del hem_pad, hem, fused, dil, res

        d_seed = np.hypot(pooldf['cx'] - tx, pooldf['cy'] - ty)
        n_near = int((d_seed <= match_radius).sum())

        gt_xy = gt_eval[['cx', 'cy']].to_numpy()
        gt_cls = gt_eval['category_id'].to_numpy()
        for arm, key in AXES.items():
            ranked = pooldf.sort_values(key, ascending=False, na_position='last',
                                        kind='mergesort').reset_index(drop=True)
            d2g = my_greedy_match(ranked[['cx', 'cy']].to_numpy(), gt_xy, match_radius)
            tp = np.cumsum([(g >= 0 and gt_cls[g] == MITOTIC) for g in d2g])
            rec = dict(file_name=fn, arm=arm, base_size=base_size, seed_ann_id=seed_ann,
                       n_retries=retries, n_gt_mitotic=n_gt, n_detections=len(ranked),
                       n_pre_nms_peaks=n_pre_nms, n_near_seed_annulus=n_near,
                       match_radius_px=round(match_radius, 3),
                       anchor_x=tx, anchor_y=ty,
                       anchor_half_pixel=bool(abs(tx % 1 - 0.5) < 1e-9 or abs(ty % 1 - 0.5) < 1e-9),
                       t_s=round(time.time() - t0, 1))
            for b in BUDGETS:
                rec[f'tp_at_{b}'] = int(tp[b - 1])
                rec[f'precision_at_{b}'] = round(float(tp[b - 1]) / b, 6)
            rows.append(rec)

            if trace_403 and fn == '403.tiff' and n_near:
                cand = pooldf[d_seed <= match_radius].iloc[0]
                r0 = int(ranked.index[(ranked['cx'] == cand['cx']) &
                                      (ranked['cy'] == cand['cy'])][0])
                trace[-1][f'annulus_rank0_{arm}'] = r0
                trace[-1][f'annulus_ordinal_{arm}'] = r0 + 1
                trace[-1]['annulus_xy'] = f"({cand['cx']:.0f},{cand['cy']:.0f})"
                trace[-1]['annulus_score'] = round(float(cand['score']), 3)
                trace[-1]['annulus_dist_to_anchor'] = round(
                    float(np.hypot(cand['cx'] - tx, cand['cy'] - ty)), 3)
                trace[-1]['annulus_dist_to_click'] = round(
                    float(np.hypot(cand['cx'] - float(seed['cx']),
                                   cand['cy'] - float(seed['cy']))), 3)
                sp = trace[-1]['nearest_peak_xy'].strip('()').split(',')
                trace[-1]['annulus_dist_to_self_peak'] = round(
                    float(np.hypot(cand['cx'] - float(sp[0]), cand['cy'] - float(sp[1]))), 3)
                trace[-1]['n_detections'] = len(pooldf)
        print(f'  tier B {fn}: base={base_size} pre_nms={n_pre_nms} pool={len(pooldf)} '
              f'[{time.time()-t0:.0f}s]', flush=True)
    return pd.DataFrame(rows), pd.DataFrame(trace)


def tier_b_divergences(tb, new_raw):
    """Compare every Tier B value against the notebook's own 13:24 raw CSV."""
    piv = new_raw.pivot_table(index=['file_name', 'arm'], columns='budget',
                              values=['tp_at_budget', 'budget_delivered'])
    piv.columns = [f'{a}_{b}' for a, b in piv.columns]
    piv = piv.reset_index()
    base = new_raw.drop_duplicates('file_name').set_index('file_name')
    rows = []
    for _, r in tb.iterrows():
        fn, arm = r['file_name'], r['arm']
        p = piv[(piv['file_name'] == fn) & (piv['arm'] == arm)].iloc[0]
        for c in ('base_size', 'seed_ann_id', 'n_gt_mitotic', 'n_detections'):
            rows.append(dict(file_name=fn, arm=arm, quantity=c, audit=r[c],
                             notebook=base.loc[fn, c],
                             agrees=bool(int(r[c]) == int(base.loc[fn, c]))))
        for b in BUDGETS:
            rows.append(dict(file_name=fn, arm=arm, quantity=f'tp_at_{b}',
                             audit=r[f'tp_at_{b}'], notebook=int(p[f'tp_at_budget_{b}']),
                             agrees=bool(int(r[f'tp_at_{b}']) == int(p[f'tp_at_budget_{b}']))))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------
# 5. Config drift and provenance
# ---------------------------------------------------------------------------------------

def cap_feasibility(tb):
    """Would T2-2 (`MAX_PEAKS = 100`) as specified actually run? `extract_peaks` truncates with
    `order[:max_peaks]`, so `n_peaks` becomes exactly 100 wherever more peaks clear the floor --
    and cell 5 asserts `n_peaks < MAX_PEAKS`. Checked against the measured pre-NMS peak counts."""
    per_roi = tb.drop_duplicates('file_name')[['file_name', 'n_pre_nms_peaks']]
    rows = []
    for _, r in per_roi.iterrows():
        n_at_cap = min(int(r['n_pre_nms_peaks']), 100)
        rows.append(dict(file_name=r['file_name'],
                         n_pre_nms_peaks_uncapped=int(r['n_pre_nms_peaks']),
                         n_peaks_at_MAX_PEAKS_100=n_at_cap,
                         cell5_assert_n_peaks_lt_MAX_PEAKS=bool(n_at_cap < 100),
                         check_no_cap_would_raise=bool(n_at_cap == 100)))
    return pd.DataFrame(rows)


def config_drift():
    """Three-way: notebook config cell vs FSConfig defaults vs DECISIONS.md."""
    cfg = cell_src(1)

    def val(name):
        m = re.search(rf'^{name}\s*=\s*([^\n#]+)', cfg, re.M)
        return m.group(1).strip() if m else 'ABSENT'
    import sys
    sys.path.insert(0, REPO)
    from midog_utils import find_and_suppress as fs
    d = fs.FSConfig()
    rows = [
        dict(setting='MAX_PEAKS', notebook=val('MAX_PEAKS'), dataclass_default='n/a',
             decisions='D9 (2026-09-12): 100', agrees=val('MAX_PEAKS') == '100'),
        dict(setting='CHANNEL', notebook=val('CHANNEL'), dataclass_default=repr(d.channel),
             decisions='D3 hematoxylin_od', agrees='hematoxylin_od' in val('CHANNEL')),
        dict(setting='METHOD', notebook=val('METHOD'), dataclass_default=str(d.tm_method),
             decisions='D1 TM_CCOEFF', agrees='TM_CCOEFF' in val('METHOD')
             and 'NORMED' not in val('METHOD')),
        dict(setting='PEAK_MIN_DISTANCE', notebook=val('PEAK_MIN_DISTANCE'),
             dataclass_default=str(d.peak_min_distance), decisions='repo default 7',
             agrees=val('PEAK_MIN_DISTANCE') == str(d.peak_min_distance)),
        dict(setting='SELF_HIT_RADIUS', notebook=val('SELF_HIT_RADIUS'),
             dataclass_default=str(d.self_hit_radius), decisions='FSConfig 5.0 fixed px',
             agrees=float(val('SELF_HIT_RADIUS')) == float(d.self_hit_radius)),
        dict(setting='NMS_RADIUS_UM', notebook=val('NMS_RADIUS_UM'), dataclass_default='n/a',
             decisions='D7 7.5 um', agrees='MIDOG_RADIUS_UM' in val('NMS_RADIUS_UM')),
        dict(setting='MATCH_RADIUS_UM', notebook=val('MATCH_RADIUS_UM'), dataclass_default='n/a',
             decisions='D7 7.5 um', agrees='MIDOG_RADIUS_UM' in val('MATCH_RADIUS_UM')),
        dict(setting='DEEP_FLOOR_Z', notebook=val('DEEP_FLOOR_Z'),
             dataclass_default=str(getattr(d, 'deep_floor_z', 'n/a')),
             decisions='repo convention -1.5', agrees=val('DEEP_FLOOR_Z') == '-1.5'),
        dict(setting='OD_WINDOW', notebook=val('OD_WINDOW'), dataclass_default='n/a',
             decisions='D5 amendment: od51', agrees='BASE_SIZE' in val('OD_WINDOW')),
        dict(setting='AXES (arm set)', notebook=str(list(AXES)), dataclass_default='n/a',
             decisions="D5 amendment lists od_contrast as a candidate axis",
             agrees='od_contrast' in AXES),
        dict(setting='SEED_INDEX', notebook=val('SEED_INDEX'), dataclass_default='n/a',
             decisions='one click per ROI', agrees=val('SEED_INDEX') == '0'),
    ]
    return pd.DataFrame(rows)


def provenance():
    """git status / last-commit for the notebook, its artifacts and the modules it imports."""
    def sh(cmd):
        return subprocess.run(cmd, shell=True, capture_output=True, text=True,
                              cwd=REPO).stdout.strip()
    targets = [os.path.relpath(NB, REPO)] + [
        os.path.relpath(p, REPO) for p in
        [NEW_RAW, NEW_PER_ROI, NEW_BY_DOMAIN, NEW_VERIF, NEW_STAT, NEW_DELTA_A, NEW_DELTA_B,
         OLD_RAW_P, OLD_PER_ROI]] + [
        f'midog_utils/{m}.py' for m in
        ('channels', 'chromatin', 'compare', 'dataset', 'evaluate', 'find_and_suppress',
         'invariants', 'nms', 'seed_selection', 'template_match')]
    rows = []
    for t in targets:
        full = os.path.join(REPO, t)
        rows.append(dict(path=t,
                         mtime=time.strftime('%Y-%m-%d %H:%M:%S',
                                             time.localtime(os.path.getmtime(full)))
                         if os.path.exists(full) else 'MISSING',
                         git_status=sh(f'git status --porcelain -- "{t}"') or 'clean',
                         last_commit=sh(f'git log -1 --format="%h %ad %an" --date=short -- "{t}"')
                         or 'NEVER COMMITTED'))
    return pd.DataFrame(rows)


def figure_series():
    """Re-derive the plotted series of each figure from the delta tables (step 1a)."""
    da = pd.read_csv(NEW_DELTA_A)
    db = pd.read_csv(NEW_DELTA_B)
    new_raw, old_raw, new_a, old_a, new_seed, old_seed, same = load_frames()
    m = new_a.merge(old_a, on=['file_name', 'arm'], suffixes=('_new', '_old'))
    rows = []
    for _, r in da.iterrows():
        ref = m[(m['file_name'] == r['file_name']) & (m['arm'] == r['arm'])].iloc[0]
        for k in BUDGETS:
            mine = ref[f'precision_at_{k}_new'] - ref[f'precision_at_{k}_old']
            rows.append(dict(figure='fig1_roi_heatmap', key=f"{r['file_name']}/{r['arm']}/K{k}",
                             notebook=round(float(r[f'delta_precision_at_{k}']), 6),
                             audit=round(float(mine), 6),
                             agrees=bool(np.isclose(r[f'delta_precision_at_{k}'], mine, atol=1e-9))))
        rows.append(dict(figure='fig1_star_annotation', key=r['file_name'],
                         notebook=bool(r['same_seed']),
                         audit=bool(same.loc[r['file_name']]),
                         agrees=bool(bool(r['same_seed']) == bool(same.loc[r['file_name']]))))
    tb_new = pd.read_csv(NEW_BY_DOMAIN)
    for _, r in db.iterrows():
        n = tb_new[(tb_new['domain'] == r['domain']) & (tb_new['arm'] == r['arm']) &
                   (tb_new['K'] == r['K'])].iloc[0]['precision_pooled']
        rows.append(dict(figure='fig2_domain_bars',
                         key=f"{r['domain']}/{r['arm']}/K{r['K']}",
                         notebook=round(float(r['delta_precision_pooled']), 6),
                         audit=round(float(n - r['precision_pooled_old']), 6),
                         agrees=bool(np.isclose(r['delta_precision_pooled'],
                                                n - r['precision_pooled_old'], atol=1e-9))))
    return pd.DataFrame(rows)


def composition():
    """Shape of every artifact the corrected run wrote, against the design it implies."""
    raw = pd.read_csv(NEW_RAW)
    a = pd.read_csv(NEW_PER_ROI)
    b = pd.read_csv(NEW_BY_DOMAIN)
    v = pd.read_csv(NEW_VERIF)
    s = pd.read_csv(NEW_STAT)
    da = pd.read_csv(NEW_DELTA_A)
    db = pd.read_csv(NEW_DELTA_B)
    disk = sorted(f for f in os.listdir(IMAGES_DIR) if f.endswith('.tiff'))
    checks = [
        ('raw rows = 14 ROI x n_arms x 4 budgets', 14 * len(AXES) * 4, len(raw)),
        ('raw distinct arms', len(AXES), raw['arm'].nunique()),
        ('raw distinct ROIs', 14, raw['file_name'].nunique()),
        ('raw ROIs match disk', len(disk), int(raw['file_name'].isin(disk).all()) * 14),
        ('per_roi rows', 14 * len(AXES), len(a)),
        ('by_domain rows', 7 * len(AXES) * 4, len(b)),
        ('delta per_roi rows (inner merge loses nothing)', 14 * len(AXES), len(da)),
        ('delta by_domain rows', 7 * len(AXES) * 4, len(db)),
        ('domains x 2 ROIs each', 7, raw.groupby('tumor_type')['file_name'].nunique().eq(2).sum()),
        ('duplicate (file,arm,budget) keys', 0,
         int(raw.duplicated(['file_name', 'arm', 'budget']).sum())),
        ('duplicate raw rows', 0, int(raw.duplicated().sum())),
        ('NaN tp_at_budget', 0, int(raw['tp_at_budget'].isna().sum())),
        ('budget_delivered == budget on all rows', len(raw),
         int((raw['budget_delivered'] == raw['budget']).sum())),
        ('n_detections >= 50 on all rows', len(raw), int((raw['n_detections'] >= 50).sum())),
        ('verification records', 87, len(v)),
        ('verification no_cap rows (T3-3 should have removed these)', 0,
         int((v['check'] == 'no_cap').sum())),
        ('stat_context rows = n_arms x 4 budgets', len(AXES) * 4, len(s)),
        ('seed_annulus_empty per-ROI violations', 0,
         int((v[(v['check'] == 'seed_annulus_empty') & (v['label'] != 'ALL')]['passed']
              == False).sum())),
    ]
    return pd.DataFrame([dict(check=c, expected=e, got=g, passed=bool(e == g))
                         for c, e, g in checks])


def main():
    t0 = time.time()
    print('== execution gate'); eg = execution_gate(); eg.to_csv(OUT.format('execution_gate'), index=False)
    print(eg.to_string(index=False))

    print('== correction matrix'); cm = correction_matrix()
    cm.to_csv(OUT.format('correction_matrix'), index=False)
    print(cm[['item', 'verdict']].to_string(index=False))

    print('== composition'); comp = composition(); comp.to_csv(OUT.format('composition'), index=False)
    print(comp.to_string(index=False))

    print('== provenance'); pv = provenance(); pv.to_csv(OUT.format('provenance'), index=False)

    print('== config drift'); cd = config_drift(); cd.to_csv(OUT.format('config_drift'), index=False)
    print(cd.to_string(index=False))

    print('== tier A inference'); inf = tier_a_inference()
    inf.to_csv(OUT.format('inference'), index=False)
    print(inf.to_string(index=False))

    print('== tie/zero check'); tz = tie_zero_check(); tz.to_csv(OUT.format('tie_zero'), index=False)
    print(tz.to_string(index=False))

    print('== rounding artifact'); ra = rounding_artifact()
    ra.to_csv(OUT.format('rounding_artifact'), index=False)
    print(ra.to_string(index=False))

    print('== spearman four ways'); sp = spearman_four_ways()
    sp.to_csv(OUT.format('spearman'), index=False)
    print(sp[sp['p_asymptotic'] < 0.10].to_string(index=False))

    print('== stat_context check'); sc = stat_context_check()
    sc.to_csv(OUT.format('stat_context_check'), index=False)
    print('stat_context cells compared:', len(sc), 'divergences:', int((~sc['agrees']).sum()))
    if (~sc['agrees']).any():
        print(sc[~sc['agrees']].to_string(index=False))

    print('== recall/pool context (declined T2-1 numbers)')
    rc, reach1 = recall_pool_context(); rc.to_csv(OUT.format('recall_pool'), index=False)
    print(rc.to_string(index=False)); print('reach 1.000:', reach1)

    print('== prose numbers'); pn = prose_numbers(inf, rc, reach1, sp)
    pn.to_csv(OUT.format('prose_numbers'), index=False)
    print(pn.to_string(index=False))

    print('== figure series'); fs_ = figure_series(); fs_.to_csv(OUT.format('figure_series'), index=False)
    print('figure values compared:', len(fs_), 'divergences:', int((~fs_['agrees']).sum()))

    files = sorted(f for f in os.listdir(IMAGES_DIR) if f.endswith('.tiff'))
    print(f'== tier B from pixels, {len(files)} ROIs')
    tb, tr = tier_b(files)
    tb.to_csv(OUT.format('tier_b_per_roi'), index=False)
    tr.to_csv(OUT.format('tier_b_403_trace'), index=False)
    div = tier_b_divergences(tb, pd.read_csv(NEW_RAW))
    div.to_csv(OUT.format('tier_b_divergences'), index=False)
    print('tier B values compared:', len(div), 'divergences:', int((~div['agrees']).sum()))
    if (~div['agrees']).any():
        print(div[~div['agrees']].to_string(index=False))
    print(tr.T.to_string())

    cf = cap_feasibility(tb); cf.to_csv(OUT.format('cap_feasibility'), index=False)
    print(cf.to_string(index=False))
    print(f"ROIs where `assert n_peaks < MAX_PEAKS` would FAIL at MAX_PEAKS=100: "
          f"{int((~cf['cell5_assert_n_peaks_lt_MAX_PEAKS']).sum())}/{len(cf)}")

    anchors = tb.drop_duplicates('file_name')[
        ['file_name', 'anchor_x', 'anchor_y', 'anchor_half_pixel', 'n_near_seed_annulus']]
    anchors.to_csv(OUT.format('anchor_geometry'), index=False)
    print(anchors.to_string(index=False))

    print(f'\nDONE in {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
