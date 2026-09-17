"""
    Round-2 audit of production_hematoxylin_only/bbox_refinement_three_way_chromatin_od.ipynb (re-executed 19:18-19:24).

    Checks the nine fixes and the two new cells (TIES, contested SENSITIVITY) against independent recomputation:
    TP counts re-matched from databases/MIDOG++.json, template centres and refusals from round 1's pixel re-implementation
    (results/bbox_refinement_three_way_chromatin_od_audit_tier_b_*.csv, persisted), sign-flip p in exact rational
    arithmetic. Also checks that the notebook still runs against the working-tree midog_utils, which changed during
    and after the re-execution.

    Helpers (matcher, sign flip, Holm, JSON/TIFF readers, D8 gate) are imported from round 1's script, which imports
    nothing from midog_utils. midog_utils is imported only in the provenance section, as the object under test.
    Writes results/bbox_refinement_three_way_chromatin_od_audit_round2_*.csv only.
"""

from __future__ import annotations

import importlib.util
import json
import math
import os
import re
import subprocess
import sys
import time
from fractions import Fraction

import numpy as np
import pandas as pd

REPO = '/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked'
sys.path.insert(0, REPO)
import bbox_refinement_three_way_chromatin_od_audit as r1  # noqa: E402  (round 1 helpers; defines only)

NB_PATH = r1.NB_PATH
STEM = r1.STEM
R1 = f'{REPO}/results/bbox_refinement_three_way_chromatin_od_audit'
OUT = f'{REPO}/results/bbox_refinement_three_way_chromatin_od_audit_round2'
CONDITIONS, PAIRS, BUDGETS = r1.CONDITIONS, r1.PAIRS, r1.BUDGETS
RULES = ('all_mitoses', 'contested_excluded', 'contested_as_fp')
LEDGER = []


def log(msg):
    """
        Print with flush.

        msg (str): message.

        Returns None.
    """
    print(msg, flush=True)


def check(table, key, column, nb_value, audit_value, tol=0.0):
    """
        Record one notebook-vs-audit comparison.

        table (str): source table or claim.
        key (str): row key.
        column (str): column.
        nb_value (object): notebook value.
        audit_value (object): audit value.
        tol (float): absolute tolerance for numbers.

        Returns bool: agreement.
    """
    try:
        ok = bool(abs(float(nb_value) - float(audit_value)) <= tol + 1e-12)
    except (TypeError, ValueError):
        ok = str(nb_value) == str(audit_value)
    LEDGER.append(dict(table=table, key=key, column=column, notebook=str(nb_value), audit=str(audit_value), match=ok))
    return ok


def save(df, name):
    """
        Write one round-2 table.

        df (pd.DataFrame): table.
        name (str): suffix.

        Returns None.
    """
    df.to_csv(f'{OUT}_{name}.csv', index=False)
    log(f'  -> results/{os.path.basename(OUT)}_{name}.csv ({len(df)} rows)')


def exact_sign_flip(deltas):
    """
        Exact two-sided sign-flip p for rational per-ROI deltas, in integer arithmetic (no float tolerance).

        deltas (list[Fraction]): one delta per ROI.

        Returns tuple[Fraction, int]: exact p and the number of nonzero deltas.
    """
    nz = [abs(d) for d in deltas if d != 0]
    n = len(nz)
    if n == 0:
        return Fraction(1), 0
    lcm = 1
    for d in list(nz) + [sum(deltas, Fraction(0))]:
        lcm = lcm * d.denominator // math.gcd(lcm, d.denominator)
    ints = np.array([int(d * lcm) for d in nz], dtype=np.int64)
    obs = abs(int(sum(deltas, Fraction(0)) * lcm))
    bits = (np.arange(2 ** n)[:, None] >> np.arange(n)) & 1
    sums = np.abs(((2 * bits - 1) * ints).sum(axis=1))
    return Fraction(int((sums >= obs).sum()), 2 ** n), n


def holm_exact(pvals):
    """
        Holm step-down adjusted p, computed on exact fractions.

        pvals (list[Fraction]): raw p-values.

        Returns list[Fraction]: adjusted p in input order.
    """
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    out, running = [None] * m, Fraction(0)
    for rank, i in enumerate(order):
        running = max(running, min(Fraction(1), (m - rank) * pvals[i]))
        out[i] = running
    return out


# --------------------------------------------------------------------------------------------------

def execution_and_provenance():
    """
        Execution gate, artifact/module provenance, and whether the notebook still runs on the working tree.

        Returns tuple[pd.DataFrame, pd.DataFrame, dict]: execution gate, provenance, and the run-now result.
    """
    nb = json.load(open(NB_PATH))
    rows = []
    for i, c in enumerate(nb['cells']):
        if c['cell_type'] == 'code':
            ex = c.get('metadata', {}).get('execution', {})
            rows.append(dict(cell=i, execution_count=c.get('execution_count'), has_error=any(o['output_type'] == 'error' for o in c.get('outputs', [])), started_utc=ex.get('iopub.execute_input'), idle_utc=ex.get('iopub.status.idle')))
    ex = pd.DataFrame(rows)
    ex['contiguous_from_1'] = ex['execution_count'].tolist() == list(range(1, len(ex) + 1))

    files = [NB_PATH] + [f'{STEM}_{s}.csv' for s in ('per_roi', 'summary', 'by_domain', 'delta_per_roi', 'delta_stats', 'top30', 'verification', 'contested_sensitivity')] + [f'{REPO}/production_hematoxylin_only/bbox3way_{s}.png' for s in ('pooled_precision', 'roi_delta_heatmap', 'win_counts')] + [r1.REF_RAW] + [f'{REPO}/midog_utils/{m}.py' for m in ('channels', 'dataset', 'evaluate', 'production', 'seed_selection', 'template_match', 'find_and_suppress', 'nms', 'chromatin', 'invariants', '__init__')]
    prov = []
    for f in files:
        rel = os.path.relpath(f, REPO)
        st = subprocess.run(['git', '-C', REPO, 'status', '--porcelain', '--', rel], capture_output=True, text=True).stdout.strip()
        prov.append(dict(file=rel, mtime=time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(f))), git_status=st[:2] if st else 'clean'))
    prov = pd.DataFrame(prov)

    # names the notebook uses from midog_utils that the working tree no longer defines
    src = '\n'.join(''.join(c['source']) for c in nb['cells'] if c['cell_type'] == 'code')
    used = sorted(set(re.findall(r'\b(ds|ss|ev|tm|ch|prod)\.([A-Za-z_]\w*)', src)))
    probe = 'import sys, json; sys.path.insert(0, ".."); from midog_utils import channels as ch, dataset as ds, evaluate as ev, production as prod, seed_selection as ss, template_match as tm; print(json.dumps({a + "." + n: hasattr(eval(a), n) for a, n in ' + repr(used) + '}))'
    res = subprocess.run([sys.executable, '-c', probe], capture_output=True, text=True, cwd=f'{REPO}/production_hematoxylin_only')
    present = json.loads(res.stdout.strip().splitlines()[-1])
    missing = [k for k, v in present.items() if not v]
    # reproduce the first statements of cell 4 exactly as written (cells 1 and 3 define only; nothing is written before the failure point)
    cell1 = ''.join(nb['cells'][1]['source'])
    cell4_head = '\n'.join(''.join(nb['cells'][4]['source']).split('\n')[:2])
    runner = f'import os\nos.chdir({repr(REPO + "/production_hematoxylin_only")})\n' + cell1.replace("print(CONFIG_DRIFT.to_string(index=False))", "") + '\n' + cell4_head + '\nprint("CELL4_HEAD_OK")'
    run = subprocess.run([sys.executable, '-c', runner], capture_output=True, text=True, cwd=f'{REPO}/production_hematoxylin_only')
    tail = (run.stdout + run.stderr).strip().splitlines()[-1] if (run.stdout + run.stderr).strip() else ''
    return ex, prov, dict(missing_names=missing, cell4_head_runs=('CELL4_HEAD_OK' in run.stdout), cell4_error=tail)


def module_equivalence(images, ann, tmeta, files):
    """
        HEAD vs working-tree seed_selection.tightened_template_box on every border-filtered unanimous candidate, both channels.

        images (pd.DataFrame): image table.
        ann (pd.DataFrame): annotations.
        tmeta (dict): per-file TIFF metadata.
        files (list[str]): ROI files.

        Returns pd.DataFrame: per ROI and channel, candidates compared and disagreements.
    """
    head_src = subprocess.run(['git', '-C', REPO, 'show', 'HEAD:midog_utils/seed_selection.py'], capture_output=True, text=True).stdout
    import midog_utils  # noqa: F401  (package context for the relative import in the HEAD source)
    from midog_utils import seed_selection as ss_now
    spec = importlib.util.spec_from_loader('midog_utils._seed_selection_head', loader=None)
    ss_head = importlib.util.module_from_spec(spec)
    ss_head.__package__ = 'midog_utils'
    sys.modules[spec.name] = ss_head
    exec(compile(head_src, 'HEAD:midog_utils/seed_selection.py', 'exec'), ss_head.__dict__)
    rows = []
    meta = images.set_index('file_name')
    for fn in files:
        rgb = r1.load_rgb(f'{r1.IMG_DIR}/{fn}')
        chans = dict(gray_inverted=r1.to_gray_inv(rgb), hematoxylin_od=r1.to_hem_od(rgb))
        shape = rgb.shape
        del rgb
        g = ann[(ann['image_id'] == int(meta.loc[fn, 'image_id'])) & (ann['category_id'] == r1.MITOTIC)]
        g = g[g['n_mitotic_votes'] == g['n_votes']]
        ix, iy = np.rint(g['cx'].to_numpy()).astype(int), np.rint(g['cy'].to_numpy()).astype(int)
        g = g[(ix >= 36) & (ix <= shape[1] - 37) & (iy >= 36) & (iy <= shape[0] - 37)]
        for name, chan in chans.items():
            n_diff = n_acc = 0
            for cx, cy in zip(g['cx'], g['cy']):
                a = ss_head.tightened_template_box(chan, float(cx), float(cy), otsu_window=51)
                b = ss_now.tightened_template_box(chan, float(cx), float(cy), otsu_window=51)
                n_diff += int(a != b)
                n_acc += int(b is not None)
            rows.append(dict(file_name=fn, channel=name, candidates=len(g), accepted_now=n_acc, head_vs_worktree_disagreements=n_diff))
        log(f'  [{fn}] module equivalence done')
    eq = pd.DataFrame(rows)
    from midog_utils import production as prod_now
    top = pd.read_csv(f'{STEM}_top30.csv')
    pr = pd.read_csv(f'{STEM}_per_roi.csv')
    cols = ['file_name', 'condition', 'tpl_cx', 'tpl_cy', 'base_size']
    geo = pd.read_csv(f'{R1}_tier_b_pipeline_runs.csv')[cols].set_index(['file_name', 'condition'])
    run_rows = []
    for fn in ('245.tiff', '403.tiff'):
        rgb = r1.load_rgb(f'{r1.IMG_DIR}/{fn}')
        mpp = tmeta[fn]['mpp_x']
        for cond in CONDITIONS:
            g = geo.loc[(fn, cond)]
            seed = ss_now.Seed(ann_id=int(pr[(pr['file_name'] == fn)]['seed_ann_id'].iloc[0]), click_xy=(0.0, 0.0), template_xy=(float(g['tpl_cx']), float(g['tpl_cy'])), base_size=int(g['base_size']), recentred=cond != 'default_51', offset_px=0.0, n_retries=0, agreement_flagged=False, n_agreement_pool=0, n_after_border=0)
            det, info = prod_now.run_production_pipeline(rgb, seed, mpp, rank_key='chromatin_od')
            nb = top[(top['file_name'] == fn) & (top['condition'] == cond)].sort_values('rank')
            run_rows.append(dict(file_name=fn, condition=cond, top30_identical=bool(np.array_equal(det[['cx', 'cy']].to_numpy()[:30], nb[['cx', 'cy']].to_numpy())), n_detections_now=int(info['n_detections']), n_detections_notebook=int(pr[(pr['file_name'] == fn) & (pr['condition'] == cond)]['n_detections'].iloc[0])))
        del rgb
    return eq, pd.DataFrame(run_rows)


def main():
    """
        Run every round-2 check and write the tables.

        Returns None.
    """
    t0 = time.time()
    images, ann = r1.load_db(r1.DB_PATH)
    files = sorted(f for f in os.listdir(r1.IMG_DIR) if f.endswith('.tiff'))
    meta = images.set_index('file_name')
    tmeta = {fn: r1.tiff_meta(f'{r1.IMG_DIR}/{fn}') for fn in files}
    roi_order = meta.loc[files].reset_index().sort_values(['tumor_type', 'file_name'])['file_name'].tolist()

    log('== execution, provenance, runnability ==')
    ex, prov, runnow = execution_and_provenance()
    save(ex, 'execution_gate')
    save(prov, 'provenance')
    save(pd.DataFrame([dict(missing_midog_utils_names_used_by_notebook=';'.join(runnow['missing_names']), cell4_first_two_statements_run=runnow['cell4_head_runs'], error=runnow['cell4_error'])]), 'runnability')
    log(f"  contiguous {bool(ex['contiguous_from_1'].iloc[0])}, errors {int(ex['has_error'].sum())}; missing names {runnow['missing_names']}; cell 4 head runs now: {runnow['cell4_head_runs']} ({runnow['cell4_error']})")

    per_roi = pd.read_csv(f'{STEM}_per_roi.csv')
    summary = pd.read_csv(f'{STEM}_summary.csv')
    by_domain = pd.read_csv(f'{STEM}_by_domain.csv')
    dpr = pd.read_csv(f'{STEM}_delta_per_roi.csv')
    dstats = pd.read_csv(f'{STEM}_delta_stats.csv')
    top = pd.read_csv(f'{STEM}_top30.csv')
    verif = pd.read_csv(f'{STEM}_verification.csv')
    sens = pd.read_csv(f'{STEM}_contested_sensitivity.csv')

    log('== value identity with the 17:40 run (round-1 persisted tables) ==')
    r1det = pd.read_csv(f'{R1}_top30_rematch_detections.csv')
    r1lists = pd.read_csv(f'{R1}_tier_b_full_lists_chromatin_od.csv')
    r1ledger = pd.read_csv(f'{R1}_comparisons.csv')
    m = top.merge(r1det, on=['file_name', 'condition', 'rank'], how='outer', indicator=True, suffixes=('', '_r1'))
    check('identity', 'top30', 'rows present in both runs', 1260, int((m['_merge'] == 'both').sum()))
    for c, c1 in (('cx', 'cx_r1'), ('cy', 'cy_r1'), ('od', 'od_r1'), ('bucket', 'nb_bucket'), ('matched_ann_id', 'nb_matched_ann')):
        eq = (m[c] == m[c1]) if c in ('bucket', 'matched_ann_id', 'cx', 'cy') else np.isclose(m[c], m[c1], rtol=0, atol=0)
        check('identity', 'top30 vs 17:40 top30', c, 1260, int(eq.sum()))
    pix = r1lists[r1lists['rank'] < 30].merge(top, on=['file_name', 'condition', 'rank'], suffixes=('_pix', ''))
    check('identity', 'top30 vs round-1 pixel re-implementation', 'identical (cx, cy)', 1260, int(((pix['cx'] == pix['cx_pix']) & (pix['cy'] == pix['cy_pix'])).sum()))
    old_nb = r1ledger[r1ledger['table'].isin(['per_roi', 'summary', 'by_domain', 'delta_per_roi'])]
    log(f'  round-1 ledger rows for per_roi/summary/by_domain/delta_per_roi: {len(old_nb)} (all matched in round 1)')
    check('identity', 'verification', 'rows / passed', '267/267', f"{len(verif)}/{int(verif['passed'].sum())}")

    log('== TP re-match from the JSON, all derived tables ==')
    tp_rows = []
    for (fn, cond), g in top.groupby(['file_name', 'condition'], sort=False):
        pr = per_roi[(per_roi['file_name'] == fn) & (per_roi['condition'] == cond)].iloc[0]
        gt = r1.build_gt(ann, int(meta.loc[fn, 'image_id']), int(pr['seed_ann_id']))
        radius = r1.RADIUS_UM / tmeta[fn]['mpp_x']
        g = g.sort_values('rank')
        det_xy = g[['cx', 'cy']].to_numpy(float)
        idx, _ = r1.greedy(det_xy, gt[['cx', 'cy']].to_numpy(float), radius)
        cat = gt['category_id'].to_numpy()
        contested = ((gt['category_id'] == r1.MITOTIC) & (gt['n_mitotic_votes'] < gt['n_votes'])).to_numpy()
        hit = np.array([i >= 0 and cat[i] == r1.MITOTIC for i in idx])
        chit = np.array([i >= 0 and cat[i] == r1.MITOTIC and contested[i] for i in idx])
        matched = [int(gt['ann_id'].iloc[i]) if i >= 0 else -1 for i in idx]
        check('top30', f'{fn}/{cond}', 'bucket agreement (30)', 30, int(sum((b == 'human_correct_label') == h for b, h in zip(g['bucket'], hit))))
        for k in BUDGETS:
            tp_rows.append(dict(file_name=fn, condition=cond, K=k, tp=int(hit[:k].sum()), contested_hits=int(chit[:k].sum()), matched_mitoses=frozenset(a for a, h in zip(matched[:k], hit[:k]) if h), coords=frozenset(map(tuple, det_xy[:k]))))
            check('per_roi', f'{fn}/{cond}', f'tp_at_{k}', int(pr[f'tp_at_{k}']), int(hit[:k].sum()))
            check('per_roi', f'{fn}/{cond}', f'precision_at_{k}', pr[f'precision_at_{k}'], hit[:k].sum() / k, tol=1e-12)
    tp = pd.DataFrame(tp_rows)
    T = tp.set_index(['file_name', 'condition', 'K'])

    for _, r in summary.iterrows():
        s = tp[(tp['condition'] == r['condition']) & (tp['K'] == r['K'])]
        p = s['tp'] / r['K']
        key = f"{r['condition']}/K{r['K']}"
        check('summary', key, 'tp_sum', r['tp_sum'], int(s['tp'].sum()))
        check('summary', key, 'precision_pooled', r['precision_pooled'], round(s['tp'].sum() / (14 * r['K']), 4), tol=1e-12)
        check('summary', key, 'worst_roi_precision', r['worst_roi_precision'], round(p.min(), 4), tol=1e-12)
        check('summary', key, 'n_roi_at_worst', r['n_roi_at_worst'], int((p == p.min()).sum()))
        check('summary', key, 'best_roi_precision', r['best_roi_precision'], round(p.max(), 4), tol=1e-12)
    dom = tp.merge(images[['file_name', 'tumor_type']], on='file_name')
    for _, r in by_domain.iterrows():
        s = dom[(dom['condition'] == r['condition']) & (dom['K'] == r['K']) & (dom['tumor_type'] == r['domain'])]
        check('by_domain', f"{r['domain']}/{r['condition']}/K{r['K']}", 'precision_pooled', r['precision_pooled'], round(s['tp'].sum() / (2 * r['K']), 4), tol=1e-12)
        check('by_domain', f"{r['domain']}/{r['condition']}/K{r['K']}", 'tp_sum', r['tp_sum'], int(s['tp'].sum()))
    for _, r in dpr.iterrows():
        a, b = r['pair'].split(' - ')
        d = int(T.loc[(r['file_name'], a, r['K']), 'tp'] - T.loc[(r['file_name'], b, r['K']), 'tp'])
        check('delta_per_roi', f"{r['file_name']}/{r['pair']}/K{r['K']}", 'delta_tp', r['delta_tp'], d)

    # delta stats incl. the new columns, exact
    st_rows = []
    for a, b in PAIRS:
        for k in BUDGETS:
            d = [Fraction(int(T.loc[(fn, a, k), 'tp'] - T.loc[(fn, b, k), 'tp']), k) for fn in roi_order]
            p, n = exact_sign_flip(d)
            st_rows.append(dict(pair=f'{a} - {b}', K=k, mean=float(sum(d, Fraction(0)) / 14), wins=sum(x > 0 for x in d), losses=sum(x < 0 for x in d), ties=sum(x == 0 for x in d), n_nonzero=n, p=p, min_p=Fraction(2, 2 ** n) if n else Fraction(1)))
    hp = holm_exact([r['p'] for r in st_rows])
    for r, h in zip(st_rows, hp):
        r['holm'] = h
    stx = pd.DataFrame(st_rows)
    for _, r in dstats.iterrows():
        a = stx[(stx['pair'] == r['pair']) & (stx['K'] == r['K'])].iloc[0]
        key = f"{r['pair']}/K{r['K']}"
        check('delta_stats', key, 'mean_delta_precision', r['mean_delta_precision'], round(a['mean'], 4), tol=1e-12)
        for c in ('wins', 'losses', 'ties', 'n_nonzero'):
            check('delta_stats', key, c, r[c], a[c])
        check('delta_stats', key, 'exact_p (6 dp)', r['exact_p'], round(float(a['p']), 6), tol=1e-12)
        check('delta_stats', key, 'min_attainable_p (6 dp)', r['min_attainable_p'], round(float(a['min_p']), 6), tol=1e-12)
        check('delta_stats', key, 'holm_p (6 dp)', r['holm_p'], round(float(a['holm']), 6), tol=1e-12)
    save(stx.assign(p=stx['p'].astype(str), min_p=stx['min_p'].astype(str), holm=stx['holm'].astype(str), p_float=[float(x) for x in stx['p']], holm_float=[float(x) for x in stx['holm']]), 'delta_stats_exact')

    log('== contested sensitivity, independent definitions, exact rational test ==')
    sens_rows, float_vs_exact = [], []
    for rule in RULES:
        block = []
        for a, b in PAIRS:
            for k in BUDGETS:
                def prec(fn, c):
                    """
                        Precision of one run at K under this rule, as an exact fraction.

                        fn (str): ROI.
                        c (str): condition.

                        Returns Fraction: precision.
                    """
                    t, h = int(T.loc[(fn, c, k), 'tp']), int(T.loc[(fn, c, k), 'contested_hits'])
                    return Fraction(t, k) if rule == 'all_mitoses' else (Fraction(t - h, k - h) if rule == 'contested_excluded' else Fraction(t - h, k))
                d = [prec(fn, a) - prec(fn, b) for fn in roi_order]
                p, n = exact_sign_flip(d)
                fl = np.array([float(x) for x in d])
                signs = (2.0 * ((np.arange(2 ** n)[:, None] >> np.arange(n)) & 1) - 1.0) if n else None
                p_float = float(np.mean(np.abs(signs @ np.abs(fl[np.abs(fl) > 1e-12])) >= abs(fl.sum()) - 1e-9)) if n else 1.0  # the notebook's float-tolerant rule
                float_vs_exact.append(dict(rule=rule, pair=f'{a} - {b}', K=k, exact=str(p), exact_float=float(p), float_tolerant=p_float, equal=abs(p_float - float(p)) < 1e-15))
                block.append(dict(rule=rule, pair=f'{a} - {b}', K=k, mean=float(sum(d, Fraction(0)) / 14), wins=sum(x > 0 for x in d), losses=sum(x < 0 for x in d), ties=sum(x == 0 for x in d), n_nonzero=n, p=p))
        for r, h in zip(block, holm_exact([r['p'] for r in block])):
            r['holm'] = h
        sens_rows.extend(block)
    sx = pd.DataFrame(sens_rows)
    for _, r in sens.iterrows():
        a = sx[(sx['rule'] == r['rule']) & (sx['pair'] == r['pair']) & (sx['K'] == r['K'])].iloc[0]
        key = f"{r['rule']}/{r['pair']}/K{r['K']}"
        check('contested_sensitivity', key, 'mean_delta_precision', r['mean_delta_precision'], round(a['mean'], 4), tol=1e-12)
        check('contested_sensitivity', key, 'wins', r['wins'], a['wins'])
        check('contested_sensitivity', key, 'losses', r['losses'], a['losses'])
        check('contested_sensitivity', key, 'exact_p (6 dp)', r['exact_p'], round(float(a['p']), 6), tol=1e-12)
        check('contested_sensitivity', key, 'holm_p (6 dp)', r['holm_p'], round(float(a['holm']), 6), tol=1e-12)
    save(sx.assign(p=sx['p'].astype(str), holm=sx['holm'].astype(str), p_float=[float(x) for x in sx['p']], holm_float=[float(x) for x in sx['holm']]), 'contested_sensitivity_exact')
    r1rules = pd.read_csv(f'{R1}_delta_stats_by_scoring_rule.csv').set_index(['rule', 'pair', 'K'])
    for mine, theirs in (('contested_neutral', 'contested_excluded'), ('contested_mitoses_dropped', 'contested_as_fp')):
        for a, b in PAIRS:
            for k in BUDGETS:
                x = r1rules.loc[(mine, f'{a} - {b}', k)]
                y = sx[(sx['rule'] == theirs) & (sx['pair'] == f'{a} - {b}') & (sx['K'] == k)].iloc[0]
                check('rule_definitions', f'{theirs} vs round-1 {mine}/{a} - {b}/K{k}', 'exact_p', round(float(x['exact_p']), 9), round(float(y['p']), 9), tol=1e-9)
                check('rule_definitions', f'{theirs} vs round-1 {mine}/{a} - {b}/K{k}', 'mean_delta', round(float(x['mean_delta_precision']), 9), round(y['mean'], 9), tol=1e-9)
    fve = pd.DataFrame(float_vs_exact)
    save(fve, 'sign_flip_float_vs_exact')
    share = tp.groupby(['condition', 'K'])[['contested_hits', 'tp']].sum()
    share['share'] = share['contested_hits'] / share['tp']
    save(share.reset_index(), 'contested_share')
    min_den = int((tp['K'] - tp['contested_hits']).min())
    log(f"  contested share {share['share'].min():.3f}-{share['share'].max():.3f}; smallest contested_excluded denominator {min_den}; float-tolerant == exact on {int(fve['equal'].sum())}/{len(fve)}")

    log('== TIES recomputed ==')
    runs1 = pd.read_csv(f'{R1}_tier_b_pipeline_runs.csv').set_index(['file_name', 'condition'])  # pixel-derived template geometry
    base = per_roi.set_index(['file_name', 'condition'])['base_size']
    tie_rows = []
    for a, b in PAIRS:
        for k in BUDGETS:
            n_diff, shared_xy, both_all, both_diff = 0, 0, 0, []
            for fn in roi_order:
                same = runs1.loc[(fn, a), 'base_size'] == runs1.loc[(fn, b), 'base_size'] and (round(runs1.loc[(fn, a), 'tpl_cx']), round(runs1.loc[(fn, a), 'tpl_cy'])) == (round(runs1.loc[(fn, b), 'tpl_cx']), round(runs1.loc[(fn, b), 'tpl_cy']))
                same_px = same  # identical pixels: same channel, same odd size, same rounded centre
                n_both = len(T.loc[(fn, a, k), 'matched_mitoses'] & T.loc[(fn, b, k), 'matched_mitoses'])
                both_all += n_both
                if not same_px:
                    n_diff += 1
                    shared_xy += len(T.loc[(fn, a, k), 'coords'] & T.loc[(fn, b, k), 'coords'])
                    both_diff.append(n_both)
            gap = float(np.mean([abs(int(base.loc[(fn, a)]) - int(base.loc[(fn, b)])) for fn in roi_order]))
            ties = int(sum(T.loc[(fn, a, k), 'tp'] == T.loc[(fn, b, k), 'tp'] for fn in roi_order))
            tie_rows.append(dict(pair=f'{a} - {b}', K=k, mean_abs_size_gap_px=round(gap, 1), ties=ties, rois_with_different_template=n_diff, shared_topk_coordinates_where_templates_differ=shared_xy, mean_mitoses_found_by_both_all_14=round(both_all / 14, 2), mean_mitoses_found_by_both_where_templates_differ=round(float(np.mean(both_diff)), 2)))
    ties_df = pd.DataFrame(tie_rows)
    save(ties_df, 'ties_recomputed')
    nb = json.load(open(NB_PATH))
    ties_txt = ''.join(''.join(o['data']['text/plain']) for o in nb['cells'][19]['outputs'] if 'data' in o)
    for line in ties_txt.strip().splitlines()[1:]:
        parts = line.split()
        pair, k = f'{parts[1]} - {parts[3]}', int(parts[4])
        r = ties_df[(ties_df['pair'] == pair) & (ties_df['K'] == k)].iloc[0]
        for j, c in enumerate(('mean_abs_size_gap_px', 'ties', 'rois_with_different_template', 'shared_topk_coordinates_where_templates_differ', 'mean_mitoses_found_by_both_all_14')):
            check('TIES', f'{pair}/K{k}', c, float(parts[5 + j]), r[c], tol=1e-9)

    log('== refusal attribution vs round-1 pixel walk ==')
    draws = pd.read_csv(f'{R1}_tier_b_seed_draws.csv')
    geo_txt = ''.join(''.join(o['data']['text/plain']) for o in nb['cells'][10]['outputs'] if 'data' in o)
    stream = ''.join(''.join(o.get('text', '')) for o in nb['cells'][10]['outputs'] if o['output_type'] == 'stream')
    nb_refused = {m.group(1): m.group(2) for m in re.finditer(r'^(\d{3}\.tiff)\s+\S.*?\s(\d+)\s+(\S+:\S+)\s', geo_txt, flags=re.M)}
    nb_refused = {m.group(1): m.group(3) for m in re.finditer(r'^(\d{3}\.tiff)\s+.+?\s+(\d)\s+(\d+:[a-z_0-9+]+)\s+', geo_txt, flags=re.M)}
    ref_rows = []
    gray_differs_audit = []
    for fn in files:
        d = draws[(draws['file_name'] == fn) & (draws['draw'] != 'summary') & (draws['joint_accept'].astype(str) == 'False')]
        entries = []
        for _, r in d.iterrows():
            refusing = []
            if str(r['gray_template_readable']) != 'True':
                refusing.append('gray_bbox')
            if str(r['hem_template_readable']) != 'True':
                refusing.append('hem_bbox')
            entries.append(f"{int(r['ann_id'])}:{'+'.join(refusing)}")
        audit = ';'.join(entries)
        if any('gray_bbox' not in e.split(':')[1] for e in entries):
            gray_differs_audit.append(fn)
        ref_rows.append(dict(file_name=fn, refused_draws_audit=audit, refused_draws_notebook=nb_refused.get(fn, '')))
        check('refusals', fn, 'refused_draws', nb_refused.get(fn, ''), audit)
    save(pd.DataFrame(ref_rows), 'refusal_attribution')
    nb_gray = re.search(r"accepted by the gray gate on (\[.*?\])", stream).group(1)
    nb_default = re.search(r"refused at least one draw on (\[.*?\])", stream).group(1)
    check('refusals', 'ALL', 'production gray_bbox would click differently', nb_gray, str(sorted(gray_differs_audit)))
    check('refusals', 'ALL', 'default_51 alone would click differently', nb_default, str(sorted(draws[(draws['draw'] == 'summary') & draws['joint_accept'].str.contains('n_retries=1')]['file_name'].tolist())))

    log('== closing-summary claims ==')
    stx_i = stx.set_index(['pair', 'K'])
    sx_i = sx.set_index(['rule', 'pair', 'K'])
    tie_i = ties_df.set_index(['pair', 'K'])
    pooled = {(c, k): tp[(tp['condition'] == c) & (tp['K'] == k)]['tp'].sum() / (14 * k) for c in CONDITIONS for k in BUDGETS}
    claims = []

    def claim(text, stated, recomputed, holds):
        """
            Record one closing-summary claim.

            text (str): the claim.
            stated (str): as written.
            recomputed (str): as recomputed.
            holds (bool): verdict.

            Returns None.
        """
        claims.append(dict(claim=text, stated=stated, recomputed=recomputed, holds=bool(holds)))

    fmt = lambda x, n: f'{float(x):.{n}f}'
    claim('pooled table 9 values (3 dp)', '0.671 0.657 0.636 / 0.632 0.604 0.564 / 0.567 0.550 0.524', ' / '.join(' '.join(fmt(pooled[(c, k)], 3) for c in CONDITIONS) for k in BUDGETS), ' / '.join(' '.join(fmt(pooled[(c, k)], 3) for c in CONDITIONS) for k in BUDGETS) == '0.671 0.657 0.636 / 0.632 0.604 0.564 / 0.567 0.550 0.524')
    hd = [stx_i.loc[('hem_bbox - default_51', k)] for k in BUDGETS]
    claim('hem-default size 3.6/6.8/4.3', '3.6/6.8/4.3', '/'.join(fmt(-100 * r['mean'], 1) for r in hd), '/'.join(fmt(-100 * r['mean'], 1) for r in hd) == '3.6/6.8/4.3')
    claim('hem-default losses 5/11/10, wins 0/1/1', '5/11/10; 0/1/1', f"{'/'.join(str(r['losses']) for r in hd)}; {'/'.join(str(r['wins']) for r in hd)}", [r['losses'] for r in hd] == [5, 11, 10] and [r['wins'] for r in hd] == [0, 1, 1])
    win013 = sorted({fn for fn in files for k in BUDGETS if T.loc[(fn, 'hem_bbox', k), 'tp'] > T.loc[(fn, 'default_51', k), 'tp']})
    claim('the one hem-default win is 013', '013.tiff', str(win013), win013 == ['013.tiff'])
    claim('hem-default exact p 0.0625/0.0029/0.0068', '0.0625/0.0029/0.0068', '/'.join(str(r['p']) for r in hd), [fmt(hd[0]['p'], 4), fmt(hd[1]['p'], 4), fmt(hd[2]['p'], 4)] == ['0.0625', '0.0029', '0.0068'])
    claim('hem-default Holm 0.44/0.026/0.055', '0.44/0.026/0.055', '/'.join(fmt(r['holm'], 6) for r in hd), [fmt(hd[0]['holm'], 2), fmt(hd[1]['holm'], 3), fmt(hd[2]['holm'], 3)] == ['0.44', '0.026', '0.055'])
    claim('K=10 hem-default: only 5 ROIs differ, min p 0.0625', '5; 0.0625', f"{hd[0]['n_nonzero']}; {hd[0]['min_p']}", hd[0]['n_nonzero'] == 5 and hd[0]['min_p'] == Fraction(1, 16))
    claim('contested share 22.5-25.5% of top-K TPs', '22.5-25.5%', f"{100 * share['share'].min():.1f}-{100 * share['share'].max():.1f}%", (round(100 * share['share'].min(), 1), round(100 * share['share'].max(), 1)) == (22.5, 25.5))
    ce = [sx_i.loc[('contested_excluded', 'hem_bbox - default_51', k)] for k in BUDGETS]
    cf = [sx_i.loc[('contested_as_fp', 'hem_bbox - default_51', k)] for k in BUDGETS]
    claim('contested_excluded keeps it: Holm 0.029 at K=20 and K=30', '0.029/0.029', f"{float(ce[1]['holm']):.6f}/{float(ce[2]['holm']):.6f}", fmt(ce[1]['holm'], 3) == '0.029' and fmt(ce[2]['holm'], 3) == '0.029')
    claim('contested_as_fp does not: p 0.064/0.037, Holm 0.52/0.33', '0.064/0.037; 0.52/0.33', f"{float(cf[1]['p']):.6f}/{float(cf[2]['p']):.6f}; {float(cf[1]['holm']):.6f}/{float(cf[2]['holm']):.6f}", [fmt(cf[1]['p'], 3), fmt(cf[2]['p'], 3), fmt(cf[1]['holm'], 2), fmt(cf[2]['holm'], 2)] == ['0.064', '0.037', '0.52', '0.33'])
    neg = all(sx_i.loc[(rule, 'hem_bbox - default_51', k), 'mean'] < 0 for rule in RULES for k in BUDGETS)
    claim('hem-default direction negative at every K under all three rules', 'negative 9/9', f'{neg}', neg)
    resolves = {rule: [k for k in BUDGETS if sx_i.loc[(rule, 'hem_bbox - default_51', k), 'holm'] <= Fraction(1, 20)] for rule in RULES}
    claim("heading: hem-default 'resolves, but only under this TP convention'", 'resolves only under all_mitoses', f'Holm <= 0.05 at K: {resolves}', len(resolves['contested_excluded']) == 0)
    gd = [stx_i.loc[('gray_bbox - default_51', k)] for k in BUDGETS]
    claim('gray-default -1.4/-2.9/-1.7; ties 10/10/6; losses 3/4/5; wins 1/0/3; p 0.625/0.125/0.297', 'as stated', '; '.join(['/'.join(fmt(100 * r['mean'], 1) for r in gd), '/'.join(str(r['ties']) for r in gd), '/'.join(str(r['losses']) for r in gd), '/'.join(str(r['wins']) for r in gd), '/'.join(fmt(r['p'], 3) for r in gd)]), ['/'.join(fmt(100 * r['mean'], 1) for r in gd), [r['ties'] for r in gd], [r['losses'] for r in gd], [r['wins'] for r in gd], [fmt(r['p'], 3) for r in gd]] == ['-1.4/-2.9/-1.7', [10, 10, 6], [3, 4, 5], [1, 0, 3], ['0.625', '0.125', '0.297']])
    claim('gray-default: 4 non-tied at K=10 and 20, min attainable p 0.125 (cannot resolve)', '4, 0.125', f"{gd[0]['n_nonzero']}/{gd[1]['n_nonzero']}; {gd[0]['min_p']}/{gd[1]['min_p']}", gd[0]['n_nonzero'] == gd[1]['n_nonzero'] == 4 and gd[0]['min_p'] == Fraction(1, 8))
    lose = sorted({fn for fn in files for k in BUDGETS if T.loc[(fn, 'gray_bbox', k), 'tp'] < T.loc[(fn, 'default_51', k), 'tp']})
    win = sorted({fn for fn in files for k in BUDGETS if T.loc[(fn, 'gray_bbox', k), 'tp'] > T.loc[(fn, 'default_51', k), 'tp']})
    claim('gray-default loses at some K on 246,402,459,529,548; wins on 094,233,301', '246,402,459,529,548 | 094,233,301', f'{lose} | {win}', lose == ['246.tiff', '402.tiff', '459.tiff', '529.tiff', '548.tiff'] and win == ['094.tiff', '233.tiff', '301.tiff'])
    hg = [stx_i.loc[('hem_bbox - gray_bbox', k)] for k in BUDGETS]
    claim('hem-gray -2.1/-3.9/-2.6; W/L/T 0/3/11, 3/9/2, 3/7/4; p 0.25/0.090/0.098; Holm 0.75/0.54/0.54', 'as stated', '; '.join(['/'.join(fmt(-100 * r['mean'], 1) for r in hg), ', '.join(f"{r['wins']}/{r['losses']}/{r['ties']}" for r in hg), '/'.join(fmt(r['p'], 4) for r in hg), '/'.join(fmt(r['holm'], 4) for r in hg)]), ['/'.join(fmt(-100 * r['mean'], 1) for r in hg), [(r['wins'], r['losses'], r['ties']) for r in hg], [fmt(hg[0]['p'], 2), fmt(hg[1]['p'], 3), fmt(hg[2]['p'], 3)], [fmt(r['holm'], 2) for r in hg]] == ['2.1/3.9/2.6', [(0, 3, 11), (3, 9, 2), (3, 7, 4)], ['0.25', '0.090', '0.098'], ['0.75', '0.54', '0.54']])
    k10 = [tie_i.loc[(f'{a} - {b}', 10)] for a, b in PAIRS]
    claim('K=10 ties 9-11 of 14; mean gap 9.9-21.7 px', '9-11; 9.9-21.7', f"{min(r['ties'] for r in k10)}-{max(r['ties'] for r in k10)}; {min(r['mean_abs_size_gap_px'] for r in k10)}-{max(r['mean_abs_size_gap_px'] for r in k10)}", (min(r['ties'] for r in k10), max(r['ties'] for r in k10)) == (9, 11))
    claim('K=10, where templates differ: 1-2 shared coordinates in total', '1-2', '/'.join(str(r['shared_topk_coordinates_where_templates_differ']) for r in k10), all(1 <= r['shared_topk_coordinates_where_templates_differ'] <= 2 for r in k10))
    claim('K=10, where templates differ: lists find 5.2-6.1 of the same mitoses per ROI', '5.2-6.1 (scoped to ROIs where templates differ)', f"all 14 ROIs: {'/'.join(str(r['mean_mitoses_found_by_both_all_14']) for r in k10)}; ROIs where templates differ: {'/'.join(str(r['mean_mitoses_found_by_both_where_templates_differ']) for r in k10)}", max(r['mean_mitoses_found_by_both_where_templates_differ'] for r in k10) < 6.15)
    k20 = [tie_i.loc[(f'{a} - {b}', 20)] for a, b in PAIRS]
    claim('K=20: gray-default (11.9 px) 10 ties; hem-default (21.7) and hem-gray (9.9, smallest) 2 each', '10/2/2; 11.9/21.7/9.9', f"{'/'.join(str(r['ties']) for r in k20)}; {'/'.join(str(r['mean_abs_size_gap_px']) for r in k20)}", [r['ties'] for r in k20] == [10, 2, 2] and [r['mean_abs_size_gap_px'] for r in k20] == [11.9, 21.7, 9.9])
    claim('caveat: default_51 alone differs on 013,245,300,403; production gray on 013,245,300', '013,245,300,403 | 013,245,300', f"{nb_default} | {sorted(gray_differs_audit)}", sorted(gray_differs_audit) == ['013.tiff', '245.tiff', '300.tiff'])
    save(pd.DataFrame(claims), 'closing_summary_claims')

    log('== figure labels (Figure 3 now .3g) ==')
    fig = [dict(panel=f'{a} vs {b}', K=k, p_label=f"p={round(float(stx_i.loc[(f'{a} - {b}', k), 'p']), 6):.3g}") for a, b in PAIRS for k in BUDGETS]
    save(pd.DataFrame(fig), 'figure3_p_labels_expected')

    if '--skip-module-equivalence' not in sys.argv:
        log('== HEAD vs working-tree seed_selection (pixels) ==')
        eqv, runs_now = module_equivalence(images, ann, tmeta, files)
        save(eqv, 'module_equivalence')
        save(runs_now, 'production_worktree_spotcheck')
        for _, r in runs_now.iterrows():
            check('worktree_spotcheck', f"{r['file_name']}/{r['condition']}", 'top30 identical with current production.py', True, r['top30_identical'])
            check('worktree_spotcheck', f"{r['file_name']}/{r['condition']}", 'n_detections', r['n_detections_notebook'], r['n_detections_now'])
        log(f"  candidates x channels {int(eqv['candidates'].sum())}, disagreements {int(eqv['head_vs_worktree_disagreements'].sum())}")

    led = pd.DataFrame(LEDGER)
    save(led, 'comparisons')
    cnt = led.groupby('table')['match'].agg(compared='size', divergences=lambda s: int((~s).sum())).reset_index()
    save(cnt, 'comparison_counts')
    log(cnt.to_string(index=False))
    cl = pd.DataFrame(claims)
    log(f"closing-summary claims: {int(cl['holds'].sum())}/{len(cl)} hold")
    log(f'total {len(led)} compared, {int((~led["match"]).sum())} divergences; {time.time() - t0:.0f}s')


if __name__ == '__main__':
    main()
