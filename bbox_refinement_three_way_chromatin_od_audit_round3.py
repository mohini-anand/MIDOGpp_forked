"""
    Round-3 audit of production_hematoxylin_only/bbox_refinement_three_way_chromatin_od.ipynb (re-executed 20:13-20:17
    against the finished-but-uncommitted midog_utils dead-code cleanup).

    Checks the four round-2 fixes, every claim in the rewritten closing summary, and whether the notebook runs against
    the cleaned midog_utils:
      - execution gate, module mtimes vs kernel start, git status of the cleanup;
      - Tier C: a sibling copy executed with only its output paths redirected to a temp dir, every CSV compared
        cell-by-cell as text and every printed output / figure compared against the notebook's own;
      - F1: HEAD dataset.check_invariants vs the inlined asserts, on the real table and on mutated copies;
      - Tier A: TP re-matched from databases/MIDOG++.json, every derived table, exact rational sign-flip and Holm,
        the contested rules, the TIES table (all nine values of the renamed column), refusal attribution;
      - value identity of every round-2-compared value between the 19:23 and 20:16 runs;
      - Tier B: the production seed draw (production_pipeline/run_pipeline.select_annotation) on all 14 ROIs.

    Helpers come from round 1's script (imports nothing from midog_utils) and round 2's (exact sign flip, exact Holm).
    midog_utils is imported only as the object under test (runnability, F1 base table, production seed draw).
    Writes results/bbox_refinement_three_way_chromatin_od_audit_round3_*.csv only; the Tier C copy is deleted.
"""

from __future__ import annotations

import base64
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from fractions import Fraction
from types import SimpleNamespace

import numpy as np
import pandas as pd
from PIL import Image

REPO = '/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked'
sys.path.insert(0, REPO)
import bbox_refinement_three_way_chromatin_od_audit as r1  # noqa: E402  (round-1 helpers; defines only)
import bbox_refinement_three_way_chromatin_od_audit_round2 as r2  # noqa: E402  (round-2 exact tests; defines only)

PY = '/Users/mohinianand/anaconda3/bin/python3'
JUPYTER = '/Users/mohinianand/anaconda3/bin/jupyter'
NB_PATH = r1.NB_PATH
NB_DIR = os.path.dirname(NB_PATH)
STEM = r1.STEM
R1, R2 = r2.R1, r2.OUT
OUT = f'{REPO}/results/bbox_refinement_three_way_chromatin_od_audit_round3'
CONDITIONS, PAIRS, BUDGETS, RULES = r1.CONDITIONS, r1.PAIRS, r1.BUDGETS, r2.RULES
CSV_SUFFIXES = ('per_roi', 'summary', 'by_domain', 'delta_per_roi', 'delta_stats', 'top30', 'verification', 'contested_sensitivity')
FIG_SUFFIXES = ('pooled_precision', 'roi_delta_heatmap', 'win_counts')
PRE_CLEANUP = 'af8bd0d'      # last commit before the dead-code cleanup (still defines dataset.check_invariants)
CLEANUP = 'c066829'          # the cleanup commit, made 20:30:49 while this audit was running
LEDGER, CLAIMS = [], []


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
        Write one round-3 table.

        df (pd.DataFrame): table.
        name (str): suffix.

        Returns None.
    """
    df.to_csv(f'{OUT}_{name}.csv', index=False)
    log(f'  -> results/{os.path.basename(OUT)}_{name}.csv ({len(df)} rows)')


def norm(text):
    """
        Collapse whitespace, so wrapped markdown can be searched for a sentence.

        text (str): text.

        Returns str: normalised text.
    """
    return re.sub(r'\s+', ' ', text).strip()


def claim(cid, quote, summary_text, stated, recomputed, holds, source):
    """
        Record one closing-summary claim; the quote must appear verbatim (whitespace-normalised) in cell 28.

        cid (str): claim id.
        quote (str): the words of the claim as written in the notebook.
        summary_text (str): normalised cell-28 text.
        stated (str): the number(s) as stated.
        recomputed (str): the number(s) as recomputed.
        holds (bool): verdict.
        source (str): what the recomputation rests on.

        Returns None.
    """
    CLAIMS.append(dict(id=cid, quote=quote, quote_in_cell_28=norm(quote) in summary_text, stated=stated, recomputed=recomputed, holds=bool(holds), source=source))


def fmt(x, n):
    """
        Fixed-point string.

        x (object): number.
        n (int): decimals.

        Returns str: formatted.
    """
    return f'{float(x):.{n}f}'


# --------------------------------------------------------------------------------------------------
# execution, provenance, runnability
# --------------------------------------------------------------------------------------------------

def module_snapshot():
    """
        mtime and git status of every midog_utils module (present or deleted) and the production entry script.

        Returns pd.DataFrame: one row per file.
    """
    status = subprocess.run(['git', '-C', REPO, 'status', '--porcelain', '--', 'midog_utils', 'production_pipeline/run_pipeline.py'], capture_output=True, text=True).stdout
    st = {line[3:]: line[:2].strip() for line in status.splitlines()}
    names = sorted({f'midog_utils/{f}' for f in os.listdir(f'{REPO}/midog_utils') if f.endswith('.py')} | {k for k in st if k.endswith('.py')} | {'production_pipeline/run_pipeline.py'})
    rows = []
    for rel in names:
        p = f'{REPO}/{rel}'
        rows.append(dict(file=rel, exists=os.path.exists(p), mtime=time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(p))) if os.path.exists(p) else '', git_status=st.get(rel, 'clean')))
    return pd.DataFrame(rows)


def execution_gate(nb):
    """
        Execution counts, errors, and kernel timestamps per code cell.

        nb (dict): notebook JSON.

        Returns pd.DataFrame: one row per code cell.
    """
    rows = []
    for i, c in enumerate(nb['cells']):
        if c['cell_type'] == 'code':
            ex = c.get('metadata', {}).get('execution', {})
            rows.append(dict(cell=i, execution_count=c.get('execution_count'), has_error=any(o['output_type'] == 'error' for o in c.get('outputs', [])), started_utc=ex.get('iopub.execute_input'), idle_utc=ex.get('iopub.status.idle')))
    df = pd.DataFrame(rows)
    df['contiguous_from_1'] = df['execution_count'].tolist() == list(range(1, len(df) + 1))
    return df


def utc_to_local(stamp):
    """
        Kernel ISO UTC timestamp to local epoch seconds.

        stamp (str): e.g. 2026-09-17T00:13:11.789950Z.

        Returns float: epoch seconds.
    """
    return pd.Timestamp(stamp).timestamp()


def runnability_now(nb):
    """
        Does the notebook's use of midog_utils still resolve against the working tree, right now.

        nb (dict): notebook JSON.

        Returns tuple[pd.DataFrame, pd.DataFrame]: the run-now summary, and every backticked identifier in the prose with where it resolves.
    """
    src = '\n'.join(''.join(c['source']) for c in nb['cells'] if c['cell_type'] == 'code')
    used = set(re.findall(r'\b(ds|ss|ev|tm|ch|prod)\.([A-Za-z_]\w*)', src))
    cfg_line = next(line for line in ''.join(nb['cells'][1]['source']).splitlines() if line.startswith('REFERENCE_CONFIG = dict('))
    used |= {('prod', k) for k in re.findall(r'(\w+)=', cfg_line)}  # read through getattr(prod, k), invisible to the regex above
    used = sorted(used)
    probe = 'import sys, json; sys.path.insert(0, ".."); from midog_utils import channels as ch, dataset as ds, evaluate as ev, production as prod, seed_selection as ss, template_match as tm; print(json.dumps({a + "." + n: hasattr(eval(a), n) for a, n in ' + repr(used) + '}))'
    res = subprocess.run([PY, '-c', probe], capture_output=True, text=True, cwd=NB_DIR)
    present = json.loads(res.stdout.strip().splitlines()[-1])
    missing = [k for k, v in present.items() if not v]

    cell1 = ''.join(nb['cells'][1]['source'])
    lines4 = ''.join(nb['cells'][4]['source']).split('\n')
    second_assert = [i for i, line in enumerate(lines4) if line.startswith('assert')][1]  # the two inlined invariant asserts come first
    head4 = '\n'.join(lines4[:second_assert + 1])
    runner = f'import os\nos.chdir({repr(NB_DIR)})\n' + cell1 + '\n' + head4 + '\nprint("CELL4_HEAD_OK", len(annotations))'
    run = subprocess.run([PY, '-c', runner], capture_output=True, text=True, cwd=NB_DIR)
    tail = (run.stdout + run.stderr).strip().splitlines()[-1] if (run.stdout + run.stderr).strip() else ''
    summary = pd.DataFrame([dict(midog_utils_names_used=len(used), missing_names=';'.join(missing), cell1_and_cell4_through_asserts_run=('CELL4_HEAD_OK' in run.stdout), last_line=tail)])

    # every backticked identifier in markdown and code comments, and where it resolves in the working tree
    prose = '\n'.join(''.join(c['source']) for c in nb['cells'] if c['cell_type'] == 'markdown')
    prose += '\n' + '\n'.join(line.split('#', 1)[1] for line in src.splitlines() if '#' in line)
    tree = {}
    for d in ('midog_utils', 'production_pipeline'):
        for f in os.listdir(f'{REPO}/{d}'):
            if f.endswith('.py'):
                tree[f'{d}/{f}'] = open(f'{REPO}/{d}/{f}').read()
    csv_cols = set()
    for s in CSV_SUFFIXES:
        csv_cols |= set(pd.read_csv(f'{STEM}_{s}.csv', nrows=1).columns)
    rows = []
    tokens = set(re.findall(r'`([^`\n]+)`', prose)) | {m for m in re.findall(r'\b((?:dataset|evaluate|production|seed_selection|template_match|channels|chromatin|invariants|find_and_suppress)\.[A-Za-z_]\w*)', prose)}
    for tok in sorted(tokens):
        t = tok.strip().rstrip('()')
        if re.search(r'\.(py|md|csv|ipynb|tiff)$', t):
            cands = [f'{REPO}/{t}', f'{NB_DIR}/{t}', f'{REPO}/Research Logs/{os.path.basename(t)}', f'{REPO}/midog_utils/{t}', f'{r1.IMG_DIR}/{t}']
            where = 'file exists' if any(os.path.exists(p) for p in cands) else 'FILE MISSING'
        elif re.fullmatch(r'[A-Za-z_][\w.]*', t):
            name = t.split('.')[-1]
            hits = [f for f, text in tree.items() if re.search(rf'\b{re.escape(name)}\b', text)]
            if hits:
                where = 'defined/used in ' + ','.join(sorted(hits)[:3])
            elif name in csv_cols:
                where = 'notebook CSV column'
            elif re.search(rf'\b{re.escape(name)}\b', src):
                where = 'notebook code'
            else:
                where = 'UNRESOLVED'
        else:
            continue
        rows.append(dict(token=tok, resolves=where))
    return summary, pd.DataFrame(rows)


def tier_c(nb):
    """
        Execute a sibling copy of the notebook with only OUT_STEM and FIG_STEM redirected, and compare everything it produces.

        nb (dict): the notebook as persisted.

        Returns tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]: run record, CSV identity, printed-output identity, figure-file identity.
    """
    tmp = tempfile.mkdtemp(prefix='bbox3way_round3_tierc_')
    copy_path = f'{NB_DIR}/.audit_tmp_bbox3way_round3.ipynb'
    cp = json.loads(json.dumps(nb))
    src1 = ''.join(cp['cells'][1]['source'])
    old_out, old_fig = "OUT_STEM = '../results/precision_at_k_14roi_prodseed_chromatin_bbox3way'", "FIG_STEM = 'bbox3way'"
    assert src1.count(old_out) == 1 and src1.count(old_fig) == 1, 'cell 1 output paths changed; update the redirect'
    src1 = src1.replace(old_out, f"OUT_STEM = '{tmp}/precision_at_k_14roi_prodseed_chromatin_bbox3way'").replace(old_fig, f"FIG_STEM = '{tmp}/bbox3way'")
    cp['cells'][1]['source'] = src1
    code = '\n'.join(''.join(c['source']) for c in cp['cells'] if c['cell_type'] == 'code')
    writes = re.findall(r'\.(to_csv|savefig)\(([^,)]*)', code)
    assert all(arg.strip().startswith(('OUT_', "f'{FIG_STEM}")) for _, arg in writes), f'a write not routed through OUT_*/FIG_STEM: {writes}'
    for c in cp['cells']:
        if c['cell_type'] == 'code':
            c['outputs'], c['execution_count'] = [], None
    t0 = time.time()
    try:
        json.dump(cp, open(copy_path, 'w'))
        proc = subprocess.run([JUPYTER, 'nbconvert', '--to', 'notebook', '--execute', '--ExecutePreprocessor.timeout=21600', '--output-dir', tmp, '--output', 'executed.ipynb', copy_path], capture_output=True, text=True, cwd=NB_DIR)
    finally:
        if os.path.exists(copy_path):
            os.remove(copy_path)
    secs = time.time() - t0
    run = pd.DataFrame([dict(started=time.strftime('%H:%M:%S', time.localtime(t0)), seconds=round(secs), returncode=proc.returncode, copy_deleted=not os.path.exists(copy_path), stderr_tail=proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else '')])
    if proc.returncode != 0:
        log(proc.stderr[-3000:])
        return run, pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    ex = json.load(open(f'{tmp}/executed.ipynb'))
    run['execution_counts_contiguous'] = [c.get('execution_count') for c in ex['cells'] if c['cell_type'] == 'code'] == list(range(1, 18))
    run['errors'] = sum(o['output_type'] == 'error' for c in ex['cells'] for o in c.get('outputs', []))

    csv_rows = []
    for s in CSV_SUFFIXES:
        a = pd.read_csv(f'{STEM}_{s}.csv', dtype=str, keep_default_na=False)
        b = pd.read_csv(f'{tmp}/precision_at_k_14roi_prodseed_chromatin_bbox3way_{s}.csv', dtype=str, keep_default_na=False)
        same_shape = a.shape == b.shape and list(a.columns) == list(b.columns)
        diffs = int((a.values != b.values).sum()) if same_shape else -1
        csv_rows.append(dict(csv=s, rows=len(a), columns=a.shape[1], cells_compared=a.size if same_shape else 0, text_divergences=diffs))
        check('tier_c_csv', s, 'cells textually identical', a.size, a.size - diffs if same_shape else -1)

    def mask(text):
        """
            Mask wall-clock timings and the redirected output path.

            text (str): printed output.

            Returns str: masked text.
        """
        text = text.replace(tmp, '../results')
        text = re.sub(r'\[\s*\d+(\.\d+)?s\]', '[T]', text)
        return re.sub(r'\bin \d+s\b', 'in Ts', text)

    def images(cell):
        """
            Decoded PNG outputs of one cell.

            cell (dict): notebook cell.

            Returns list[np.ndarray]: RGBA arrays.
        """
        return [np.array(Image.open(io.BytesIO(base64.b64decode(o['data']['image/png'])))) for o in cell.get('outputs', []) if 'image/png' in o.get('data', {})]

    def text(cell):
        """
            Concatenated stream and text/plain outputs of one cell.

            cell (dict): notebook cell.

            Returns str: text.
        """
        return ''.join(''.join(o.get('text', '')) if o['output_type'] == 'stream' else ''.join(o.get('data', {}).get('text/plain', '')) if o['output_type'] == 'execute_result' else '' for o in cell.get('outputs', []))

    out_rows = []
    for i, (c_nb, c_ex) in enumerate(zip(nb['cells'], ex['cells'])):
        if c_nb['cell_type'] != 'code':
            continue
        t_nb, t_ex = mask(text(c_nb)), mask(text(c_ex))
        im_nb, im_ex = images(c_nb), images(c_ex)
        img_same = len(im_nb) == len(im_ex) and all(x.shape == y.shape and np.array_equal(x, y) for x, y in zip(im_nb, im_ex))
        out_rows.append(dict(cell=i, text_chars=len(t_nb), text_identical_after_masking_timings=t_nb == t_ex, n_images=len(im_nb), images_pixel_identical=img_same))
        check('tier_c_outputs', f'cell {i}', 'printed text (timings masked)', True, t_nb == t_ex)
        if im_nb or im_ex:
            check('tier_c_outputs', f'cell {i}', 'embedded figure pixels', True, img_same)
    fig_rows = []
    for s in FIG_SUFFIXES:
        a = np.array(Image.open(f'{NB_DIR}/bbox3way_{s}.png'))
        b = np.array(Image.open(f'{tmp}/bbox3way_{s}.png'))
        same = a.shape == b.shape and np.array_equal(a, b)
        fig_rows.append(dict(figure=f'bbox3way_{s}.png', shape=str(a.shape), pixel_identical=same))
        check('tier_c_figures', s, 'saved PNG pixels', True, same)
    shutil.rmtree(tmp, ignore_errors=True)
    return run, pd.DataFrame(csv_rows), pd.DataFrame(out_rows), pd.DataFrame(fig_rows)


# --------------------------------------------------------------------------------------------------
# F1: the inlined invariants
# --------------------------------------------------------------------------------------------------

def f1_equivalence(nb, ann_json):
    """
        HEAD dataset.check_invariants against the two asserts now inlined in cell 4, on the real table and on mutations.

        nb (dict): notebook JSON.
        ann_json (pd.DataFrame): annotations parsed from the JSON by round 1 (independent w/h/category).

        Returns tuple[pd.DataFrame, pd.DataFrame]: per-case outcomes, and the constants check.
    """
    head_src = subprocess.run(['git', '-C', REPO, 'show', f'{PRE_CLEANUP}:midog_utils/dataset.py'], capture_output=True, text=True).stdout
    now_src = open(f'{REPO}/midog_utils/dataset.py').read()
    consts = {}
    for name in ('BOX_SIZE', 'MITOTIC', 'LOOKALIKE'):
        h = int(re.search(rf'^{name} = (\d+)', head_src, flags=re.M).group(1))
        n = int(re.search(rf'^{name} = (\d+)', now_src, flags=re.M).group(1))
        consts[name] = dict(head=h, working_tree=n, equal=h == n)
    fn_src = re.search(r'^def check_invariants\(.*?(?=^def )', head_src, flags=re.M | re.S).group(0)
    ns = dict(pd=pd, BOX_SIZE=consts['BOX_SIZE']['head'], MITOTIC=consts['MITOTIC']['head'], LOOKALIKE=consts['LOOKALIKE']['head'])
    exec(compile(fn_src, f'{PRE_CLEANUP}:midog_utils/dataset.py', 'exec'), ns)
    head_check = ns['check_invariants']
    asserts = [line for line in ''.join(nb['cells'][4]['source']).split('\n') if line.startswith('assert')][:2]  # the two inlined invariants; later asserts belong to the sweep
    assert all('annotations[' in x for x in asserts), asserts
    ds_ns = SimpleNamespace(**{k: v['working_tree'] for k, v in consts.items()})

    def inline_check(annotations):
        """
            Run the notebook's inlined asserts verbatim.

            annotations (pd.DataFrame): table.

            Returns None. Raises AssertionError.
        """
        exec('\n'.join(asserts), dict(annotations=annotations, ds=ds_ns))

    sys.path.insert(0, REPO)
    from midog_utils import dataset as ds_now
    _, base = ds_now.load_annotations(r1.DB_PATH)
    j = ann_json.set_index('ann_id').loc[base['ann_id']]
    indep = dict(rows=len(base), w_equal_json=bool((base['w'].to_numpy() == j['w'].to_numpy()).all()), h_equal_json=bool((base['h'].to_numpy() == j['h'].to_numpy()).all()), category_equal_json=bool((base['category_id'].to_numpy() == j['category_id'].to_numpy()).all()))

    cases = {'real table': base}
    for label, col, val in (('one w = 49', 'w', 49), ('one h = 51', 'h', 51), ('one w = NaN', 'w', np.nan), ('one category = 3', 'category_id', 3), ('one category = 0', 'category_id', 0)):
        m = base.copy()
        if isinstance(val, float) and np.isnan(val):
            m[col] = m[col].astype(float)
        m.loc[m.index[0], col] = val
        cases[label] = m
    cases['look-alikes only'] = base[base['category_id'] == ds_ns.LOOKALIKE]
    cases['empty table'] = base.iloc[:0]
    rows = []
    for label, df in cases.items():
        out = {}
        for who, fn in (('head', head_check), ('inline', inline_check)):
            try:
                fn(df)
                out[who] = 'passes'
            except AssertionError:
                out[who] = 'raises'
        rows.append(dict(case=label, pre_cleanup_check_invariants=out['head'], inlined_asserts=out['inline'], same=out['head'] == out['inline']))
        check('f1_invariants', label, 'af8bd0d check_invariants vs inlined outcome', out['head'], out['inline'])
    const_df = pd.DataFrame([dict(constant=k, **v) for k, v in consts.items()] + [dict(constant=f'working-tree load_annotations vs JSON: {k}', head=None, working_tree=None, equal=v) for k, v in indep.items() if k != 'rows'])
    return pd.DataFrame(rows), const_df


# --------------------------------------------------------------------------------------------------
# Tier B: the production seed draw
# --------------------------------------------------------------------------------------------------

def production_seed_draws(files, per_roi):
    """
        Production's own click (select_annotation -> build_seed) on gray_inverted and on hematoxylin_od, per ROI.

        files (list[str]): ROI files.
        per_roi (pd.DataFrame): the notebook's per-ROI table (joint click).

        Returns pd.DataFrame: per ROI, production gray / hem clicks next to the joint click and round 1's pixel walk.
    """
    sys.path.insert(0, f'{REPO}/production_pipeline')
    import run_pipeline as rp
    from midog_utils import channels as ch
    from midog_utils import dataset as ds
    images, annotations = ds.load_annotations(r1.DB_PATH)
    draws = pd.read_csv(f'{R1}_tier_b_seed_draws.csv')
    summ = draws[draws['draw'] == 'summary'].set_index('file_name')
    rows = []
    for fn in files:
        rgb = ds.load_roi(f'{r1.IMG_DIR}/{fn}')
        gt = ds.image_annotations(annotations, fn)
        gtm = gt[gt['category_id'] == ds.MITOTIC]
        image_id = int(images.loc[images['file_name'] == fn, 'image_id'].iloc[0])
        got = {}
        for name in ('gray_inverted', 'hematoxylin_od'):
            seed = rp.select_annotation(gtm, ch.to_channel(rgb, name), rgb.shape, image_id, seed_index=0)
            got[name] = int(seed.ann_id)
        del rgb
        joint = int(per_roi.loc[per_roi['file_name'] == fn, 'seed_ann_id'].iloc[0])
        r1_gray = int(re.search(r'=(\d+)', summ.loc[fn, 'gray_reason']).group(1))
        r1_hem = int(re.search(r'=(\d+)', summ.loc[fn, 'hem_reason']).group(1))
        rows.append(dict(file_name=fn, joint_click=joint, production_gray_click=got['gray_inverted'], production_hem_click=got['hematoxylin_od'], round1_pixel_gray_click=r1_gray, round1_pixel_hem_click=r1_hem, production_gray_differs_from_joint=got['gray_inverted'] != joint))
        check('production_seed', fn, 'production gray click == round-1 pixel walk', r1_gray, got['gray_inverted'])
        check('production_seed', fn, 'production hem click == round-1 pixel walk', r1_hem, got['hematoxylin_od'])
        log(f'  [{fn}] joint {joint}, production gray {got["gray_inverted"]}, production hem {got["hematoxylin_od"]}')
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------

def main():
    """
        Run every round-3 check and write the tables.

        Returns None.
    """
    t0 = time.time()
    nb = json.load(open(NB_PATH))
    images, ann = r1.load_db(r1.DB_PATH)
    files = sorted(f for f in os.listdir(r1.IMG_DIR) if f.endswith('.tiff'))
    meta = images.set_index('file_name')
    tmeta = {fn: r1.tiff_meta(f'{r1.IMG_DIR}/{fn}') for fn in files}
    roi_order = meta.loc[files].reset_index().sort_values(['tumor_type', 'file_name'])['file_name'].tolist()
    summary_text = norm(''.join(nb['cells'][28]['source']))

    log('== execution gate, module provenance ==')
    snap0 = module_snapshot()
    save(snap0, 'module_snapshot')
    ex = execution_gate(nb)
    save(ex, 'execution_gate')
    kernel_start = utc_to_local(ex['started_utc'].iloc[0])
    kernel_end = utc_to_local(ex['idle_utc'].iloc[-1])
    latest_module = max(os.path.getmtime(f'{REPO}/{f}') for f in snap0.loc[snap0['exists'], 'file'])
    prov = pd.DataFrame([dict(kernel_start=time.strftime('%H:%M:%S', time.localtime(kernel_start)), kernel_end=time.strftime('%H:%M:%S', time.localtime(kernel_end)), latest_module_mtime=time.strftime('%H:%M:%S', time.localtime(latest_module)), all_modules_older_than_kernel_start=latest_module < kernel_start, contiguous_1_to_17=bool(ex['contiguous_from_1'].iloc[0]) and len(ex) == 17, errors=int(ex['has_error'].sum()), notebook_mtime=time.strftime('%H:%M:%S', time.localtime(os.path.getmtime(NB_PATH))))])
    git = lambda *a: subprocess.run(['git', '-C', REPO, *a], capture_output=True, text=True)
    prov['head'] = git('log', '-1', '--format=%h').stdout.strip()
    prov['cleanup_commit_time'] = git('log', '-1', '--format=%ad', '--date=format:%H:%M:%S', CLEANUP).stdout.strip()
    prov['cleanup_parent_is_pre_cleanup'] = git('rev-parse', '--short', f'{CLEANUP}^').stdout.strip() == PRE_CLEANUP
    prov['worktree_midog_utils_equals_cleanup_commit'] = git('diff', '--quiet', CLEANUP, '--', 'midog_utils', 'production_pipeline').returncode == 0
    prov['files_deleted_by_cleanup'] = ';'.join(git('diff', '--name-only', '--diff-filter=D', PRE_CLEANUP, CLEANUP, '--', 'midog_utils').stdout.split())
    prov['notebook_and_csvs_tracked'] = git('ls-files', '--error-unmatch', NB_PATH).returncode == 0
    for s in CSV_SUFFIXES:
        prov[f'{s}_mtime'] = time.strftime('%H:%M:%S', time.localtime(os.path.getmtime(f'{STEM}_{s}.csv')))
    save(prov, 'kernel_vs_modules')
    log(prov.T.to_string())

    log('== runnability now ==')
    runnow, idents = runnability_now(nb)
    save(runnow, 'runnability_now')
    save(idents, 'prose_identifiers')
    log(runnow.T.to_string())
    log(f"  prose identifiers unresolved: {idents[idents['resolves'].str.contains('UNRESOLVED|MISSING')]['token'].tolist()}")

    per_roi = pd.read_csv(f'{STEM}_per_roi.csv')
    summary = pd.read_csv(f'{STEM}_summary.csv')
    by_domain = pd.read_csv(f'{STEM}_by_domain.csv')
    dpr = pd.read_csv(f'{STEM}_delta_per_roi.csv')
    dstats = pd.read_csv(f'{STEM}_delta_stats.csv')
    top = pd.read_csv(f'{STEM}_top30.csv')
    verif = pd.read_csv(f'{STEM}_verification.csv')
    sens = pd.read_csv(f'{STEM}_contested_sensitivity.csv')

    log('== F1: pre-cleanup check_invariants vs inlined asserts ==')
    f1, consts = f1_equivalence(nb, ann)
    save(f1, 'f1_invariant_equivalence')
    save(consts, 'f1_constants')
    log(f1.to_string(index=False))

    log('== value identity with the 17:40 run and round 1 pixel lists ==')
    r1det = pd.read_csv(f'{R1}_top30_rematch_detections.csv')
    r1lists = pd.read_csv(f'{R1}_tier_b_full_lists_chromatin_od.csv')
    m = top.merge(r1det, on=['file_name', 'condition', 'rank'], how='outer', indicator=True, suffixes=('', '_r1'))
    check('identity', 'top30', 'rows present in both runs', 1260, int((m['_merge'] == 'both').sum()))
    for c, c1 in (('cx', 'cx_r1'), ('cy', 'cy_r1'), ('od', 'od_r1'), ('bucket', 'nb_bucket'), ('matched_ann_id', 'nb_matched_ann')):
        check('identity', 'top30 vs 17:40 top30', c, 1260, int((m[c] == m[c1]).sum()))
    pix = r1lists[r1lists['rank'] < 30].merge(top, on=['file_name', 'condition', 'rank'], suffixes=('_pix', ''))
    check('identity', 'top30 vs round-1 pixel re-implementation', 'identical (cx, cy)', 1260, int(((pix['cx'] == pix['cx_pix']) & (pix['cy'] == pix['cy_pix'])).sum()))
    check('identity', 'verification', 'rows / passed', '267/267', f"{len(verif)}/{int(verif['passed'].sum())}")

    log('== TP re-match from the JSON; TP-convention alternatives ==')
    tp_rows, alt_rows, both_class, click_rows = [], [], 0, []
    runs1 = pd.read_csv(f'{R1}_tier_b_pipeline_runs.csv').set_index(['file_name', 'condition'])
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
        rules, _, _ = r1.score_list(det_xy, gt, radius)
        for k in BUDGETS:
            alt_rows.append(dict(file_name=fn, condition=cond, K=k, repo=rules['repo_mixed_greedy'][k], hungarian_mitoses=rules['mitotic_only_max_matching'][k], mitoses_only_greedy=rules['mitotic_only_greedy'][k], radius_minus_1px=rules['radius_minus_1px'][k], radius_plus_1px=rules['radius_plus_1px'][k]))
        d_all = np.hypot(det_xy[:, None, 0] - gt['cx'].to_numpy()[None, :], det_xy[:, None, 1] - gt['cy'].to_numpy()[None, :]) <= radius
        both_class += int(((d_all & (cat == r1.MITOTIC)[None, :]).any(axis=1) & (d_all & (cat == 2)[None, :]).any(axis=1)).sum())
    tp = pd.DataFrame(tp_rows)
    T = tp.set_index(['file_name', 'condition', 'K'])
    alts = pd.DataFrame(alt_rows)
    save(alts, 'tp_convention_alternatives')
    for fn in files:
        a = ann[ann['image_id'] == int(meta.loc[fn, 'image_id'])]
        sid = int(per_roi.loc[per_roi['file_name'] == fn, 'seed_ann_id'].iloc[0])
        s = a[a['ann_id'] == sid].iloc[0]
        o = a[a['ann_id'] != sid]
        radius = r1.RADIUS_UM / tmeta[fn]['mpp_x']
        d = np.hypot(o['cx'] - s['cx'], o['cy'] - s['cy'])
        runs = r1lists[r1lists['file_name'] == fn]
        dc = np.hypot(runs['cx'] - s['cx'], runs['cy'] - s['cy'])
        click_rows.append(dict(file_name=fn, seed_ann_id=sid, click_cx=s['cx'], click_cy=s['cy'], radius_px=radius, nearest_other_annotation_px=float(d.min()), full_list_candidates_within_r_of_click=int((dc <= radius).sum()), full_list_candidates=len(runs)))
    clicks = pd.DataFrame(click_rows).set_index('file_name')
    save(clicks.reset_index(), 'click_neighbourhood')

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
        check('delta_per_roi', f"{r['file_name']}/{r['pair']}/K{r['K']}", 'delta_tp', r['delta_tp'], int(T.loc[(r['file_name'], a, r['K']), 'tp'] - T.loc[(r['file_name'], b, r['K']), 'tp']))

    st_rows = []
    for a, b in PAIRS:
        for k in BUDGETS:
            d = [Fraction(int(T.loc[(fn, a, k), 'tp'] - T.loc[(fn, b, k), 'tp']), k) for fn in roi_order]
            p, n = r2.exact_sign_flip(d)
            st_rows.append(dict(pair=f'{a} - {b}', K=k, mean=float(sum(d, Fraction(0)) / 14), wins=sum(x > 0 for x in d), losses=sum(x < 0 for x in d), ties=sum(x == 0 for x in d), n_nonzero=n, p=p, min_p=Fraction(2, 2 ** n) if n else Fraction(1)))
    for r, h in zip(st_rows, r2.holm_exact([r['p'] for r in st_rows])):
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

    log('== contested sensitivity ==')
    sens_rows = []
    for rule in RULES:
        block = []
        for a, b in PAIRS:
            for k in BUDGETS:
                def prec(fn, c, rule=rule, k=k):
                    """
                        Precision of one run at K under a rule, as an exact fraction.

                        fn (str): ROI.
                        c (str): condition.
                        rule (str): scoring rule.
                        k (int): budget.

                        Returns Fraction: precision.
                    """
                    t, h = int(T.loc[(fn, c, k), 'tp']), int(T.loc[(fn, c, k), 'contested_hits'])
                    return Fraction(t, k) if rule == 'all_mitoses' else (Fraction(t - h, k - h) if rule == 'contested_excluded' else Fraction(t - h, k))
                d = [prec(fn, a) - prec(fn, b) for fn in roi_order]
                p, n = r2.exact_sign_flip(d)
                block.append(dict(rule=rule, pair=f'{a} - {b}', K=k, mean=float(sum(d, Fraction(0)) / 14), wins=sum(x > 0 for x in d), losses=sum(x < 0 for x in d), ties=sum(x == 0 for x in d), n_nonzero=n, p=p))
        for r, h in zip(block, r2.holm_exact([r['p'] for r in block])):
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
    share = tp.groupby(['condition', 'K'])[['contested_hits', 'tp']].sum()
    share['share'] = share['contested_hits'] / share['tp']

    log('== TIES (renamed column, all nine values) ==')
    base = per_roi.set_index(['file_name', 'condition'])['base_size']
    same_px = {}
    for fn in roi_order:
        for a, b in PAIRS:
            same_px[(fn, a, b)] = bool(runs1.loc[(fn, a), 'base_size'] == runs1.loc[(fn, b), 'base_size'] and (round(runs1.loc[(fn, a), 'tpl_cx']), round(runs1.loc[(fn, a), 'tpl_cy'])) == (round(runs1.loc[(fn, b), 'tpl_cx']), round(runs1.loc[(fn, b), 'tpl_cy'])))
    tie_rows = []
    for a, b in PAIRS:
        for k in BUDGETS:
            differ = [fn for fn in roi_order if not same_px[(fn, a, b)]]
            shared_xy = sum(len(T.loc[(fn, a, k), 'coords'] & T.loc[(fn, b, k), 'coords']) for fn in differ)
            both = [len(T.loc[(fn, a, k), 'matched_mitoses'] & T.loc[(fn, b, k), 'matched_mitoses']) for fn in differ]
            gap = float(np.mean([abs(int(base.loc[(fn, a)]) - int(base.loc[(fn, b)])) for fn in roi_order]))
            ties = int(sum(T.loc[(fn, a, k), 'tp'] == T.loc[(fn, b, k), 'tp'] for fn in roi_order))
            tie_rows.append(dict(pair=f'{a} - {b}', K=k, mean_abs_size_gap_px=round(gap, 1), ties=ties, rois_with_different_template=len(differ), shared_topk_coordinates_where_templates_differ=shared_xy, mean_mitoses_found_by_both_where_templates_differ=round(float(np.mean(both)), 2), identical_template_rois=';'.join(fn for fn in roi_order if same_px[(fn, a, b)])))
    ties_df = pd.DataFrame(tie_rows)
    save(ties_df, 'ties_recomputed')
    header = ''.join(''.join(o['data']['text/plain']) for o in nb['cells'][19]['outputs'] if 'data' in o).strip().splitlines()
    check('TIES', 'header', 'last column name', 'mean_mitoses_found_by_both_where_templates_differ', header[0].split()[-1])
    for line in header[1:]:
        parts = line.split()
        pair, k = f'{parts[1]} - {parts[3]}', int(parts[4])
        r = ties_df[(ties_df['pair'] == pair) & (ties_df['K'] == k)].iloc[0]
        for j, c in enumerate(('mean_abs_size_gap_px', 'ties', 'rois_with_different_template', 'shared_topk_coordinates_where_templates_differ', 'mean_mitoses_found_by_both_where_templates_differ')):
            check('TIES', f'{pair}/K{k}', c, float(parts[5 + j]), r[c], tol=1e-9)
    cell19 = ''.join(nb['cells'][19]['source'])
    check('TIES', 'code', 'denominator is n_differ, and n_differ > 0 for every pair', True, 'shared_mitoses / n_differ' in cell19 and ties_df['rois_with_different_template'].min() > 0)

    log('== refusal attribution ==')
    draws = pd.read_csv(f'{R1}_tier_b_seed_draws.csv')
    geo_txt = ''.join(''.join(o['data']['text/plain']) for o in nb['cells'][10]['outputs'] if 'data' in o)
    stream10 = ''.join(''.join(o.get('text', '')) for o in nb['cells'][10]['outputs'] if o['output_type'] == 'stream')
    nb_refused = {mm.group(1): mm.group(3) for mm in re.finditer(r'^(\d{3}\.tiff)\s+.+?\s+(\d)\s+(\d+:[a-z_0-9+]+)\s+', geo_txt, flags=re.M)}
    gray_differs, refusal_sets, ref_rows = [], [], []
    for fn in files:
        d = draws[(draws['file_name'] == fn) & (draws['draw'] != 'summary') & (draws['joint_accept'].astype(str) == 'False')]
        entries = []
        for _, r in d.iterrows():
            refusing = [c for c, col in (('gray_bbox', 'gray_template_readable'), ('hem_bbox', 'hem_template_readable')) if str(r[col]) != 'True']
            refusal_sets.append(dict(file_name=fn, ann_id=int(r['ann_id']), refusing='+'.join(refusing), click_readable=str(r['click_readable'])))
            entries.append(f"{int(r['ann_id'])}:{'+'.join(refusing)}")
        if any('gray_bbox' not in e.split(':')[1] for e in entries):
            gray_differs.append(fn)
        ref_rows.append(dict(file_name=fn, refused_draws_audit=';'.join(entries), refused_draws_notebook=nb_refused.get(fn, '')))
        check('refusals', fn, 'refused_draws', nb_refused.get(fn, ''), ';'.join(entries))
    save(pd.DataFrame(ref_rows), 'refusal_attribution')
    refusal_sets = pd.DataFrame(refusal_sets)
    nb_gray = re.search(r'accepted by the gray gate on (\[.*?\])', stream10).group(1)
    nb_default = re.search(r'refused at least one draw on (\[.*?\])', stream10).group(1)
    default_differs = sorted(draws[(draws['draw'] == 'summary') & draws['joint_accept'].str.contains('n_retries=1')]['file_name'].tolist())
    check('refusals', 'ALL', 'production gray_bbox would click differently', nb_gray, str(sorted(gray_differs)))
    check('refusals', 'ALL', 'default_51 alone would click differently', nb_default, str(default_differs))

    if '--skip-tier-b' not in sys.argv:
        log('== Tier B: production seed draw on all 14 ROIs ==')
        prod_seeds = production_seed_draws(files, per_roi)
        save(prod_seeds, 'production_seed_draws')
        prod_gray_differs = sorted(prod_seeds.loc[prod_seeds['production_gray_differs_from_joint'], 'file_name'])
        check('production_seed', 'ALL', 'ROIs where production gray click != joint click', nb_gray, str(prod_gray_differs))
    else:
        prod_gray_differs = None

    log('== identity 19:23 -> 20:16 over every value round 2 compared ==')
    led3 = pd.DataFrame(LEDGER)
    led2 = pd.read_csv(f'{R2}_comparisons.csv', dtype=str, keep_default_na=False)
    keys = ['table', 'key', 'column']
    jn = led2.merge(led3.astype({'notebook': str}), on=keys, suffixes=('_1923', '_2016'))
    jn['identical'] = jn['notebook_1923'] == jn['notebook_2016']
    save(jn[keys + ['notebook_1923', 'notebook_2016', 'identical']], 'identity_1923_vs_2016')
    not_carried = led2.merge(led3[keys], on=keys, how='left', indicator=True)
    not_carried = not_carried[not_carried['_merge'] == 'left_only'].groupby(['table', 'column']).size().reset_index(name='rows')
    save(not_carried, 'round2_checks_not_repeated')
    log(f"  round-2 values re-found in round 3: {len(jn)}, identical: {int(jn['identical'].sum())}")
    log(not_carried.to_string(index=False))
    for _, r in jn[~jn['identical']].iterrows():
        log(f"  CHANGED {r['table']}/{r['key']}/{r['column']}: {r['notebook_1923']} -> {r['notebook_2016']}")

    log('== closing summary, every claim ==')
    stx_i, sx_i, tie_i = stx.set_index(['pair', 'K']), sx.set_index(['rule', 'pair', 'K']), ties_df.set_index(['pair', 'K'])
    pooled = {(c, k): Fraction(int(tp[(tp['condition'] == c) & (tp['K'] == k)]['tp'].sum()), 14 * k) for c in CONDITIONS for k in BUDGETS}
    prec = {(fn, c, k): Fraction(int(T.loc[(fn, c, k), 'tp']), k) for fn in files for c in CONDITIONS for k in BUDGETS}
    click_xy = {fn: (clicks.loc[fn, 'click_cx'], clicks.loc[fn, 'click_cy']) for fn in files}
    ref_ids = eval(re.search(r'REF_SEED_ANN_IDS = (\{.*?\})', ''.join(nb['cells'][1]['source'])).group(1))
    r1_accept = draws[draws['draw'] == 'summary'].set_index('file_name')['ann_id'].astype(int).to_dict()

    vv = verif.groupby('check')['passed'].agg(['sum', 'size'])
    claim('C01', 'the D9 cap bound on 42/42 runs', summary_text, '42/42', f"verification {vv.loc['max_peaks_binds', 'sum']}/{vv.loc['max_peaks_binds', 'size']}; round-1 pixel n_peaks==100 on {int((runs1['n_peaks'] == 100).sum())}/42", vv.loc['max_peaks_binds', 'sum'] == 42 and int((runs1['n_peaks'] == 100).sum()) == 42, 'verification CSV; round-1 pixel runs')
    nd = per_roi.set_index(['file_name', 'condition'])['n_detections']
    claim('C02', 'every K was fully delivered', summary_text, 'n_detections >= 30 on 42/42', f"min n_detections {int(nd.min())}; equal to round-1 pixel n_detections on {int((nd == runs1['n_detections'].reindex(nd.index)).sum())}/42", nd.min() >= 30 and int((nd == runs1['n_detections'].reindex(nd.index)).sum()) == 42, 'per_roi CSV; round-1 pixel runs')
    claim('C03', 'no candidate survived within one match radius of the click (42/42)', summary_text, '0 on 42/42', f"full post-NMS lists (round-1 pixels): {int(clicks['full_list_candidates_within_r_of_click'].sum())} within r of the JSON click over {int(clicks['full_list_candidates'].sum())} candidates", int(clicks['full_list_candidates_within_r_of_click'].sum()) == 0, 'JSON click, TIFF mpp, round-1 full lists')
    joint = per_roi.groupby('file_name')['seed_ann_id'].first().to_dict()
    claim('C04', "the 14 clicks match the reference notebook's", summary_text, '14/14', f"joint == REF_SEED_ANN_IDS on {sum(joint[f] == ref_ids[f] for f in files)}/14; == round-1 pixel walk on {sum(joint[f] == r1_accept[f] for f in files)}/14", all(joint[f] == ref_ids[f] == r1_accept[f] for f in files), 'per_roi CSV; cell 1 dict; round-1 pixel walk')
    ref = pd.read_csv(r1.REF_RAW)
    ref = ref[(ref['arm'] == 'chromatin_od') & ref['budget'].isin(BUDGETS)]
    g2 = 0
    for (fn, cond), g in ref.groupby(['file_name', 'condition']):
        rr = g.set_index('budget')
        ok = int(rr['base_size'].iloc[0]) == int(runs1.loc[(fn, cond), 'base_size']) and float(rr['tpl_cx'].iloc[0]) == float(runs1.loc[(fn, cond), 'tpl_cx']) and float(rr['tpl_cy'].iloc[0]) == float(runs1.loc[(fn, cond), 'tpl_cy']) and int(rr['n_detections'].iloc[0]) == int(nd.loc[(fn, cond)]) and all(int(rr.loc[k, 'tp_at_budget']) == int(T.loc[(fn, cond, k), 'tp']) for k in BUDGETS)
        g2 += int(ok)
    claim('C05', "both Otsu conditions reproduce its `chromatin_od` numbers on 28/28 rows", summary_text, '28/28', f'{g2}/{ref.groupby(["file_name", "condition"]).ngroups}', g2 == 28, 'reference raw CSV vs round-1 pixel geometry and JSON-rematched TP')
    same_gd = [fn for fn in roi_order if same_px[(fn, 'gray_bbox', 'default_51')]]
    same_lists = [fn for fn in same_gd if T.loc[(fn, 'gray_bbox', 30), 'coords'] == T.loc[(fn, 'default_51', 30), 'coords'] and all(T.loc[(fn, 'gray_bbox', k), 'tp'] == T.loc[(fn, 'default_51', k), 'tp'] for k in BUDGETS)]
    claim('C06', "`default_51` returns the identical top-30 list wherever its template is the same pixels as `gray_bbox`'s", summary_text, 'identical wherever same pixels', f"same-pixel ROIs {same_gd}; identical top-30 and TP on {same_lists} (the notebook's Gate 3 tested only offset==0: 245)", same_gd == same_lists and len(same_gd) > 0, 'round-1 pixel geometry; top30 CSV')
    r1_src = open(f'{REPO}/bbox_refinement_three_way_chromatin_od_audit.py').read()
    r1_counts = pd.read_csv(f'{R1}_comparison_counts.csv')
    r1_imports_mu = bool(re.search(r'^\s*(from|import)\s+midog_utils', r1_src, flags=re.M))
    claim('C07', 're-derived all 42 top-30 lists and every TP count from the pixels and the annotation JSON, without `midog_utils`, and found 0 divergences', summary_text, '42 lists, 0 divergences, no midog_utils', f"round-1 imports midog_utils: {r1_imports_mu}; pixel top-30 identical {int(runs1['top30_identical_to_notebook'].sum())}/42; round-1 divergences {int(r1_counts['divergences'].sum())}; 20:16 top-30 == round-1 pixel lists {int(((pix['cx'] == pix['cx_pix']) & (pix['cy'] == pix['cy_pix'])).sum())}/1260", (not r1_imports_mu) and int(runs1['top30_identical_to_notebook'].sum()) == 42 and int(r1_counts['divergences'].sum()) == 0 and int(((pix['cx'] == pix['cx_pix']) & (pix['cy'] == pix['cy_pix'])).sum()) == 1260, 'round-1 script and tables; identity chain to 20:16')
    rad = [r1.RADIUS_UM / tmeta[fn]['mpp_x'] for fn in files]
    claim('C08', '7.5 um (29.6-33.1 px per ROI)', summary_text, '29.6-33.1', f'{min(rad):.1f}-{max(rad):.1f}', (fmt(min(rad), 1), fmt(max(rad), 1)) == ('29.6', '33.1'), 'TIFF resolution tags')
    alt_ok = {c: int((alts[c] == alts['repo']).sum()) for c in ('hungarian_mitoses', 'mitoses_only_greedy', 'radius_minus_1px', 'radius_plus_1px')}
    claim('C09', 'TP counts were identical in all 126 (run, K) cells under each of', summary_text, '126/126 under Hungarian, mitoses only, radius +-1 px', str(alt_ok), all(v == 126 for v in alt_ok.values()) and len(alts) == 126, '20:16 top30 CSV rematched from JSON under each rule')
    claim('C10', 'No top-30 detection had both a mitosis and a look-alike within the radius', summary_text, '0', str(both_class), both_class == 0, '20:16 top30 CSV vs JSON')
    claim('C11', 'no other annotation lies within one radius of any click', summary_text, 'none', f"min over clicks of (nearest other annotation px - radius px) = {float((clicks['nearest_other_annotation_px'] - clicks['radius_px']).min()):.1f}", bool((clicks['nearest_other_annotation_px'] > clicks['radius_px']).all()), 'JSON')
    table = ' / '.join(' '.join(fmt(pooled[(c, k)], 3) for c in CONDITIONS) for k in BUDGETS)
    claim('C12', '| 10 | 0.671 | 0.657 | 0.636 |', summary_text, '0.671 0.657 0.636 / 0.632 0.604 0.564 / 0.567 0.550 0.524', table, table == '0.671 0.657 0.636 / 0.632 0.604 0.564 / 0.567 0.550 0.524', 'JSON-rematched TP')
    order_ok = all(pooled[('default_51', k)] > pooled[('gray_bbox', k)] > pooled[('hem_bbox', k)] for k in BUDGETS)
    claim('C13', 'pooled precision falls in the same order at every budget: `default_51` > `gray_bbox` > `hem_bbox`', summary_text, 'strict order at K=10/20/30; neither Otsu beats default', str(order_ok), order_ok, 'JSON-rematched TP')
    resolved = {f'{a} - {b}': [k for k in BUDGETS if stx_i.loc[(f'{a} - {b}', k), 'holm'] <= Fraction(1, 20)] for a, b in PAIRS}
    claim('C14', '`hem_bbox` vs. `default_51` is the one contrast that resolves.', summary_text, 'only hem-default has Holm <= 0.05', str(resolved), resolved['hem_bbox - default_51'] and not resolved['gray_bbox - default_51'] and not resolved['hem_bbox - gray_bbox'], 'exact sign flip, exact Holm over 9')
    by_rule = {rule: [k for k in BUDGETS if sx_i.loc[(rule, 'hem_bbox - default_51', k), 'holm'] <= Fraction(1, 20)] for rule in RULES}
    claim('C15', 'It survives excluding contested mitoses, but not counting them as false positives.', summary_text, 'contested_excluded resolves at some K; contested_as_fp at none', str(by_rule), bool(by_rule['contested_excluded']) and not by_rule['contested_as_fp'], 'exact, Holm within each rule')
    hd = [stx_i.loc[('hem_bbox - default_51', k)] for k in BUDGETS]
    s16 = '/'.join(fmt(-100 * r['mean'], 1) for r in hd)
    claim('C16', 'costs 3.6 / 6.8 / 4.3 points of precision at K = 10 / 20 / 30', summary_text, '3.6/6.8/4.3', s16, s16 == '3.6/6.8/4.3', 'JSON-rematched TP')
    win_hd = sorted({fn for fn in files for k in BUDGETS if T.loc[(fn, 'hem_bbox', k), 'tp'] > T.loc[(fn, 'default_51', k), 'tp']})
    claim('C17', 'it loses on 5 / 11 / 10 ROIs and wins on 0 / 1 / 1 (the one win is `013.tiff`)', summary_text, '5/11/10; 0/1/1; 013', f"{[r['losses'] for r in hd]}; {[r['wins'] for r in hd]}; {win_hd}", [r['losses'] for r in hd] == [5, 11, 10] and [r['wins'] for r in hd] == [0, 1, 1] and win_hd == ['013.tiff'], 'JSON-rematched TP')
    claim('C18', 'exact sign-flip p = 0.0625 / 0.0029 / 0.0068', summary_text, '0.0625/0.0029/0.0068', '/'.join(str(r['p']) for r in hd), [fmt(hd[0]['p'], 4), fmt(hd[1]['p'], 4), fmt(hd[2]['p'], 4)] == ['0.0625', '0.0029', '0.0068'], 'exact rational enumeration')
    claim('C19', 'that is 0.44 / 0.026 / 0.055, so K = 20 is resolved and K = 30 just misses', summary_text, '0.44/0.026/0.055; K20 <= 0.05 < K30', '/'.join(str(r['holm']) for r in hd), [fmt(hd[0]['holm'], 2), fmt(hd[1]['holm'], 3), fmt(hd[2]['holm'], 3)] == ['0.44', '0.026', '0.055'] and hd[1]['holm'] <= Fraction(1, 20) < hd[2]['holm'], 'exact Holm over 9')
    claim('C20', 'At K = 10 only 5 ROIs differ, so the smallest p the test could return is 0.0625', summary_text, '5; 0.0625', f"{hd[0]['n_nonzero']}; {hd[0]['min_p']}", hd[0]['n_nonzero'] == 5 and hd[0]['min_p'] == Fraction(1, 16), 'exact')
    claim('C21', 'about a quarter (22.5-25.5%) of top-K TPs are contested mitoses', summary_text, '22.5-25.5%', f"{100 * share['share'].min():.1f}-{100 * share['share'].max():.1f}% (per condition x K)", (fmt(100 * share['share'].min(), 1), fmt(100 * share['share'].max(), 1)) == ('22.5', '25.5'), 'JSON votes')
    ce = [sx_i.loc[('contested_excluded', 'hem_bbox - default_51', k)] for k in BUDGETS]
    claim('C22', 'keeps the result (Holm p 0.029 at K = 20 and at K = 30)', summary_text, '0.029/0.029', f"{ce[1]['holm']}/{ce[2]['holm']}", fmt(ce[1]['holm'], 3) == '0.029' and fmt(ce[2]['holm'], 3) == '0.029', 'exact')
    cf = [sx_i.loc[('contested_as_fp', 'hem_bbox - default_51', k)] for k in BUDGETS]
    claim('C23', 'Counting them as false positives does not (p = 0.064 / 0.037, Holm 0.52 / 0.33)', summary_text, '0.064/0.037; 0.52/0.33', f"{float(cf[1]['p']):.6f}/{float(cf[2]['p']):.6f}; {float(cf[1]['holm']):.6f}/{float(cf[2]['holm']):.6f}", [fmt(cf[1]['p'], 3), fmt(cf[2]['p'], 3), fmt(cf[1]['holm'], 2), fmt(cf[2]['holm'], 2)] == ['0.064', '0.037', '0.52', '0.33'], 'exact')
    neg = all(sx_i.loc[(rule, 'hem_bbox - default_51', k), 'mean'] < 0 for rule in RULES for k in BUDGETS)
    claim('C24', 'The direction stays negative at every K under all three rules', summary_text, '9/9 negative', str(neg), neg, 'exact')
    gd = [stx_i.loc[('gray_bbox - default_51', k)] for k in BUDGETS]
    got25 = ['/'.join(fmt(100 * r['mean'], 1) for r in gd), [r['ties'] for r in gd], [r['losses'] for r in gd], [r['wins'] for r in gd], [fmt(r['p'], 3) for r in gd]]
    claim('C25', 'ties on 10 / 10 / 6 ROIs, losses on 3 / 4 / 5, wins on 1 / 0 / 3 (p = 0.625 / 0.125 / 0.297)', summary_text, '-1.4/-2.9/-1.7; 10/10/6; 3/4/5; 1/0/3; 0.625/0.125/0.297', str(got25), got25 == ['-1.4/-2.9/-1.7', [10, 10, 6], [3, 4, 5], [1, 0, 3], ['0.625', '0.125', '0.297']], 'exact')
    claim('C26', 'With only 4 non-tied ROIs at K = 10 and 20, the smallest attainable p is 0.125', summary_text, '4; 0.125', f"{gd[0]['n_nonzero']}/{gd[1]['n_nonzero']}; {gd[0]['min_p']}/{gd[1]['min_p']}", gd[0]['n_nonzero'] == gd[1]['n_nonzero'] == 4 and gd[0]['min_p'] == gd[1]['min_p'] == Fraction(1, 8), 'exact')
    lose = sorted({fn for fn in files for k in BUDGETS if T.loc[(fn, 'gray_bbox', k), 'tp'] < T.loc[(fn, 'default_51', k), 'tp']})
    win = sorted({fn for fn in files for k in BUDGETS if T.loc[(fn, 'gray_bbox', k), 'tp'] > T.loc[(fn, 'default_51', k), 'tp']})
    claim('C27', 'it loses at some K on `246`, `402`, `459`, `529` and `548`, and wins at some K on `094`, `233` and `301`', summary_text, '246,402,459,529,548 | 094,233,301', f'{lose} | {win}', lose == ['246.tiff', '402.tiff', '459.tiff', '529.tiff', '548.tiff'] and win == ['094.tiff', '233.tiff', '301.tiff'], 'JSON-rematched TP')
    off = {(fn, c): float(np.hypot(runs1.loc[(fn, c), 'tpl_cx'] - click_xy[fn][0], runs1.loc[(fn, c), 'tpl_cy'] - click_xy[fn][1])) for fn in files for c in CONDITIONS}
    claim('C28', "On `245.tiff` (0 px offset) and `201.tiff` (0.5 px, which `read_padded_patch` rounds to the same pixel) the template is identical to `default_51`'s, and so is the list", summary_text, '245 0 px, 201 0.5 px; identical template and list', f"offsets 245 {off[('245.tiff', 'gray_bbox')]:.2f}, 201 {off[('201.tiff', 'gray_bbox')]:.2f}; same-pixel ROIs {same_gd}; identical lists {same_lists}", fmt(off[('245.tiff', 'gray_bbox')], 2) == '0.00' and fmt(off[('201.tiff', 'gray_bbox')], 2) == '0.50' and same_gd == same_lists == ['201.tiff', '245.tiff'], 'round-1 pixel geometry; python round() as read_padded_patch')
    ceiling = [fn for fn in roi_order if int(runs1.loc[(fn, 'gray_bbox'), 'base_size']) == 51]
    shift = [fn for fn in ceiling if not same_px[(fn, 'gray_bbox', 'default_51')]]
    s403 = len(T.loc[('403.tiff', 'gray_bbox', 30), 'coords'] & T.loc[('403.tiff', 'default_51', 30), 'coords'])
    tp403 = all(T.loc[('403.tiff', 'gray_bbox', k), 'tp'] == T.loc[('403.tiff', 'default_51', k), 'tp'] for k in BUDGETS)
    claim('C29', '`403.tiff` is the only real same-size centring shift (2 px): its top-30 list changes completely, but its TP counts at K = 10/20/30 don\'t', summary_text, 'only 403; 2 px; 0 shared of 30; TP equal', f"ceiling {ceiling}; shifted {shift}; offset {off[('403.tiff', 'gray_bbox')]:.2f}; shared top-30 coords {s403}; TP equal {tp403}", shift == ['403.tiff'] and fmt(off[('403.tiff', 'gray_bbox')], 1) == '2.0' and s403 == 0 and tp403, 'round-1 pixel geometry; top30 CSV; JSON TP')
    ref_nb = json.load(open(r1.REF_NB_PATH))
    ref_out = ''.join(''.join(o.get('text', '')) + ''.join(o.get('data', {}).get('text/plain', '')) for o in ref_nb['cells'][24]['outputs'])
    ref_rows = {int(mm.group(1)): mm.groups() for mm in re.finditer(r'^\s*\d+\s+chromatin_od\s+(\d+)\s+(-?[\d.]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)', ref_out, flags=re.M)}
    hg = [stx_i.loc[('hem_bbox - gray_bbox', k)] for k in BUDGETS]
    ref_match = all(ref_rows[k][1] == fmt(r['mean'], 4) and (int(ref_rows[k][2]), int(ref_rows[k][3]), int(ref_rows[k][4])) == (r['wins'], r['losses'], r['ties']) and ref_rows[k][6] == fmt(r['p'], 4) for k, r in zip(BUDGETS, hg))
    got30 = ['/'.join(fmt(-100 * r['mean'], 1) for r in hg), [(r['wins'], r['losses'], r['ties']) for r in hg], [fmt(hg[0]['p'], 2), fmt(hg[1]['p'], 3), fmt(hg[2]['p'], 3)], [fmt(r['holm'], 2) for r in hg]]
    claim('C30', '`hem_bbox` vs. `gray_bbox` reproduces the reference notebook exactly.** -2.1 / -3.9 / -2.6 points; wins/losses/ties 0/3/11, 3/9/2 and 3/7/4; p = 0.25 / 0.090 / 0.098 (Holm 0.75 / 0.54 / 0.54)', summary_text, 'as stated; equal to reference notebook cell 24 chromatin_od rows', f'{got30}; reference rows equal: {ref_match}', got30 == ['2.1/3.9/2.6', [(0, 3, 11), (3, 9, 2), (3, 7, 4)], ['0.25', '0.090', '0.098'], ['0.75', '0.54', '0.54']] and ref_match, 'exact; reference notebook printed output')
    worst = {(c, k): sorted(fn for fn in files if prec[(fn, c, k)] == min(prec[(f, c, k)] for f in files)) for c in CONDITIONS for k in BUDGETS}
    expected = {key: ['245.tiff'] for key in worst}
    expected[('gray_bbox', 30)] = ['245.tiff', '529.tiff']
    claim('C31', '`245.tiff` is the worst ROI under every condition at every K (tied with `529.tiff` for `gray_bbox` at K = 30)', summary_text, '245 everywhere; +529 at gray K=30 only', str({f'{c}/K{k}': v for (c, k), v in worst.items()}), worst == expected, 'JSON-rematched TP')
    p245 = {(c, k): prec[('245.tiff', c, k)] for c in CONDITIONS for k in (20, 30)}
    claim('C32', '`hem_bbox` pulls it from 0.20 down to 0.10 at K = 20, and from 0.23 to 0.17 at K = 30. `default_51` and `gray_bbox` are identical there by construction', summary_text, '0.20->0.10; 0.23->0.17; default == gray', f"K20 {fmt(p245[('default_51', 20)], 2)}->{fmt(p245[('hem_bbox', 20)], 2)}; K30 {fmt(p245[('default_51', 30)], 2)}->{fmt(p245[('hem_bbox', 30)], 2)}; default==gray {all(prec[('245.tiff', 'default_51', k)] == prec[('245.tiff', 'gray_bbox', k)] for k in BUDGETS)}", (fmt(p245[('default_51', 20)], 2), fmt(p245[('hem_bbox', 20)], 2), fmt(p245[('default_51', 30)], 2), fmt(p245[('hem_bbox', 30)], 2)) == ('0.20', '0.10', '0.23', '0.17') and '245.tiff' in same_gd, 'JSON-rematched TP; pixel geometry')
    med_size = [float(np.median([int(runs1.loc[(fn, c), 'base_size']) for fn in files])) for c in CONDITIONS]
    med_off = [float(np.median([off[(fn, c)] for fn in files])) for c in ('gray_bbox', 'hem_bbox')]
    claim('C33', 'median template size (51 -> 40 -> 28 px)', summary_text, '51/40/28; offsets 2.5 gray, 3.5 hem', f'{med_size}; {[round(x, 3) for x in med_off]}', med_size == [51.0, 40.0, 28.0] and (fmt(med_off[0], 1), fmt(med_off[1], 1)) == ('2.5', '3.5') and 'median offset 2.5 px gray, 3.5 px hem' in summary_text, 'round-1 pixel geometry; JSON click')
    k10 = [tie_i.loc[(f'{a} - {b}', 10)] for a, b in PAIRS]
    mm10 = [r['mean_mitoses_found_by_both_where_templates_differ'] for r in k10]
    claim('C34', 'all three pairs tie on 9-11 of 14 ROIs, whatever the size gap between their templates (mean 9.9-21.7 px)', summary_text, '9-11; 9.9-21.7', f"{min(r['ties'] for r in k10)}-{max(r['ties'] for r in k10)}; {min(r['mean_abs_size_gap_px'] for r in k10)}-{max(r['mean_abs_size_gap_px'] for r in k10)}", (min(r['ties'] for r in k10), max(r['ties'] for r in k10), min(r['mean_abs_size_gap_px'] for r in k10), max(r['mean_abs_size_gap_px'] for r in k10)) == (9, 11, 9.9, 21.7), 'JSON TP; per_roi base_size')
    claim('C35', 'the top-10 lists share only 1-2 exact candidate coordinates in total, yet find 5.2-6.8 of the same mitoses per ROI on average', summary_text, '1-2; 5.2-6.8', f"{[r['shared_topk_coordinates_where_templates_differ'] for r in k10]}; {mm10}", all(1 <= r['shared_topk_coordinates_where_templates_differ'] <= 2 for r in k10) and (fmt(min(mm10), 1), fmt(max(mm10), 1)) == ('5.2', '6.8'), 'top30 CSV; JSON matching; pixel geometry for same-template')
    ref32 = norm(''.join(ref_nb['cells'][32]['source']))
    ref31 = ''.join(''.join(o.get('text', '')) + ''.join(o.get('data', {}).get('text/plain', '')) for o in ref_nb['cells'][31]['outputs'])
    ref_shared = float(re.search(r'chromatin_od\s+[\d.]+\s+[\d.]+\s+([\d.]+)', ref31).group(1))
    claim('C36', 'which is the reading the reference notebook now gives for its own K = 10 ties', summary_text, 'reference reads K=10 ties as converging on the same mitoses', f"reference cell 32 says 'recover the same mitotic objects': {'recover the same mitotic objects' in ref32}; its chromatin_od hem/gray mean_shared {ref_shared} vs this notebook's hem-gray K=10 {tie_i.loc[('hem_bbox - gray_bbox', 10), 'mean_mitoses_found_by_both_where_templates_differ']}", 'recover the same mitotic objects' in ref32 and abs(ref_shared - tie_i.loc[('hem_bbox - gray_bbox', 10), 'mean_mitoses_found_by_both_where_templates_differ']) < 0.01, 'reference notebook text and printed output')
    k20 = [tie_i.loc[(f'{a} - {b}', 20)] for a, b in PAIRS]
    claim('C37', '`gray_bbox - default_51` (mean gap 11.9 px) ties on 10 ROIs, while `hem_bbox - default_51` (21.7 px) and `hem_bbox - gray_bbox` (9.9 px, the smallest gap) each tie on only 2', summary_text, '10/2/2; 11.9/21.7/9.9', f"{[r['ties'] for r in k20]}; {[r['mean_abs_size_gap_px'] for r in k20]}", [r['ties'] for r in k20] == [10, 2, 2] and [r['mean_abs_size_gap_px'] for r in k20] == [11.9, 21.7, 9.9], 'JSON TP; per_roi base_size')
    first_refusals = refusal_sets.groupby('file_name').first()
    only_hem_first = sorted(first_refusals[first_refusals['refusing'] == 'hem_bbox'].index)
    claim('C38', 'Production `gray_bbox` (gray gate alone) would also have clicked differently on `013`, `245` and `300`, where only the `hematoxylin_od` gate refused the first draw. That includes `245`, the worst ROI.', summary_text, 'default alone: 013,245,300,403; gray: 013,245,300; only hem refused first draw there', f"default {default_differs}; gray (attribution) {sorted(gray_differs)}; gray (production select_annotation) {prod_gray_differs}; first draw refused only by hem on {only_hem_first}", default_differs == ['013.tiff', '245.tiff', '300.tiff', '403.tiff'] and sorted(gray_differs) == only_hem_first == ['013.tiff', '245.tiff', '300.tiff'] and prod_gray_differs in (None, ['013.tiff', '245.tiff', '300.tiff']), 'round-1 pixel walk; production seed draw (Tier B)')
    f4_ok = bool(len(refusal_sets)) and refusal_sets['refusing'].ne('').all() and refusal_sets['click_readable'].eq('True').all()
    claim('C39', "On a click either gate refuses, at least one tightened condition doesn't exist, so the three-way comparison isn't defined there.", summary_text, 'every refused draw: >= 1 tightened condition refused, default_51 never refuses', f"{len(refusal_sets)} refused draws; refusing sets {sorted(refusal_sets['refusing'].unique())}; click readable on all: {refusal_sets['click_readable'].eq('True').all()}", f4_ok and set(refusal_sets['refusing']) <= {'gray_bbox', 'hem_bbox', 'gray_bbox+hem_bbox'}, 'round-1 pixel walk')
    steps = all(abs(round(r['precision_pooled'] * 2 * r['K']) / (2 * r['K']) - r['precision_pooled']) < 5e-5 for _, r in by_domain.iterrows())
    claim('C40', 'Domain-level numbers move in steps of 1/(2K)', summary_text, 'multiples of 1/(2K)', str(steps), steps and images.set_index('file_name').loc[files].groupby('tumor_type').size().eq(2).all(), 'by_domain CSV; JSON tumour types (2 per domain)')
    d5 = re.search(r'^## D5.*?(?=^## D6)', open(f'{REPO}/DECISIONS.md').read(), flags=re.M | re.S).group(0)
    d9 = re.search(r'^## D9.*?(?=^## |\Z)', open(f'{REPO}/DECISIONS.md').read(), flags=re.M | re.S).group(0)
    cites = dict(d5_five_seeds_14_rois='swept over 5 seeds on the 14 ROIs' in norm(d5), d5_clustered_at_roi='clustered at the ROI' in norm(d5), d9_tm_score_arm='one seed each, `tm_score` arm' in norm(d9), seed_index_0="SEED_INDEX = 0" in ''.join(nb['cells'][1]['source']))
    claim('C41', "D5's bar (5 seeds x 14 ROIs, ROI-clustered paired delta) is not met", summary_text, 'D5 amendment names 5 seeds x 14 ROIs, clustered at ROI; D9 validated on tm_score; seed_index 0', str(cites), all(cites.values()) and "D9's `MAX_PEAKS = 100` was also validated with `tm_score` ranking" in summary_text, 'DECISIONS.md text; cell 1')
    cl = pd.DataFrame(CLAIMS)
    save(cl, 'closing_summary_claims')
    log(cl[['id', 'quote_in_cell_28', 'holds', 'recomputed']].to_string(index=False, max_colwidth=110))

    if '--skip-tier-c' not in sys.argv:
        log('== Tier C: execute a redirected sibling copy against the current midog_utils ==')
        run, csv_id, out_id, fig_id = tier_c(nb)
        save(run, 'tier_c_run')
        if len(csv_id):
            save(csv_id, 'tier_c_csv_identity')
            save(out_id, 'tier_c_output_identity')
            save(fig_id, 'tier_c_figure_identity')
        log(run.T.to_string())

    snap1 = module_snapshot()
    unchanged = snap0[['file', 'exists', 'mtime']].equals(snap1[['file', 'exists', 'mtime']])
    save(pd.DataFrame([dict(modules_unchanged_during_audit=unchanged, audit_started=time.strftime('%H:%M:%S', time.localtime(t0)), audit_finished=time.strftime('%H:%M:%S'))]), 'module_stability')
    led = pd.DataFrame(LEDGER)
    save(led, 'comparisons')
    cnt = led.groupby('table')['match'].agg(compared='size', divergences=lambda s: int((~s).sum())).reset_index()
    save(cnt, 'comparison_counts')
    log(cnt.to_string(index=False))
    log(f"closing-summary claims: {int(cl['holds'].sum())}/{len(cl)} hold; quotes found in cell 28: {int(cl['quote_in_cell_28'].sum())}/{len(cl)}")
    log(f"total {len(led)} compared, {int((~led['match']).sum())} divergences; modules unchanged during audit: {unchanged}; {time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
