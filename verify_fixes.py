"""Post-run verification of every fix. Each check prints PASS/FAIL with evidence."""
import json, sys
import numpy as np, pandas as pd
from sklearn.neighbors import KDTree
sys.path.insert(0, '.')
from midog_utils import dataset as ds, evaluate as ev, experiment as ex, find_and_suppress as fs

fails = []
def check(name, ok, evidence=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {evidence}" if evidence else ""))
    if not ok: fails.append(name)

m   = pd.read_csv('results/fs_metrics.csv')
det = pd.read_csv('results/fs_detections.csv')
FS  = m[m.method == 'find_and_suppress'].set_index('file_name')
BL  = m[m.method == 'nucleus_blobs'].set_index('file_name')
RN  = m[m.method == 'random_in_tissue'].set_index('file_name')
images, ann = ds.load_annotations('databases/MIDOG++.json')

print("\n===== E1  recall_by_agreement is budgeted =====")
check("E1 columns exist", {'recall_unanimous_at_k','recall_contested_at_k'} <= set(m.columns))
sub = FS[['recall_unanimous_at_k','recall_contested_at_k']]
check("E1 no longer pinned at 1.0", not bool((sub.fillna(0) == 1.0).all().all()),
      f"values:\n{sub.to_string()}")
check("E1 full-list versions retained", {'recall_unanimous','recall_contested'} <= set(m.columns))

print("\n===== E2  blob cap removed / budget honestly reported =====")
check("E2 no blob list sits exactly at the old 20000 cap", not (BL.n_detections_total == 20000).any(),
      f"blob counts: {BL.n_detections_total.tolist()}")
check("E2 budget_delivered recorded", 'budget_delivered' in m.columns)
check("E2 budget_delivered == min(matcher_budget, blob count)",
      bool((BL.budget_delivered == np.minimum(BL.matcher_budget, BL.n_detections_total)).all()),
      BL[['matcher_budget','budget_delivered','n_detections_total']].to_string())

print("\n===== E3  tissue mask =====")
check("E3 tissue_fraction recorded", 'tissue_fraction' in m.columns)
check("E3 tissue fraction plausible (>=0.6, was 0.38-0.73 under Otsu)",
      bool((FS.tissue_fraction >= 0.6).all()), FS.tissue_fraction.to_dict())
# every mitotic GT must now lie inside the mask -> blob detector can at least see them
import cv2
from midog_utils import baselines as bl
worst = []
for fn in FS.index:
    rgb = ds.load_roi(f'images/{fn}'); mask = bl.tissue_mask(rgb)
    a = ann[(ann.file_name == fn) & (ann.category_id == ds.MITOTIC)]
    ys = np.rint(a.cy).astype(int); xs = np.rint(a.cx).astype(int)
    worst.append((fn, int((~mask[ys, xs]).sum()), len(a))); del rgb, mask
check("E3 zero mitotic GT excluded by the mask (was 31/218 on 301.tiff)",
      all(w[1] == 0 for w in worst), str(worst))

print("\n===== E4  NMS radius derived from the image =====")
check("E4 nms_radius recorded and == match radius",
      bool(np.allclose(FS.nms_radius, FS.match_radius_px, atol=0.05)),
      FS[['nms_radius','match_radius_px','mpp']].to_string())
bad = []
for fn in FS.index:
    d = det[det.file_name == fn][['cx','cy']].to_numpy()
    r = float(FS.loc[fn,'nms_radius'])
    dist, _ = KDTree(d).query(d, k=2)
    bad.append((fn, round(float(dist[:,1].min()),2), r))
check("E4 no surviving pair closer than the NMS radius",
      all(b[1] >= b[2] - 1e-6 for b in bad), str(bad))

print("\n===== E5  coverage_frac =====")
check("E5 column present for every method", m.coverage_frac.notna().all())
fn = FS.index[0]; r = float(FS.loc[fn,'match_radius_px'])
im = images[images.file_name == fn].iloc[0]
indep = ev.coverage_fraction(det[det.file_name == fn][['cx','cy']].to_numpy(),
                             (int(im.height), int(im.width)), r)
check("E5 matches an independent recomputation",
      abs(indep - float(FS.loc[fn,'coverage_frac'])) < 1e-3,
      f"{fn}: stored {FS.loc[fn,'coverage_frac']} vs recomputed {indep:.4f}")
print(f"       coverage by method: fs={FS.coverage_frac.tolist()}")

print("\n===== E6  independent per-image seeds =====")
seeds = FS.seed_ann_id.astype(int)
quantiles = {}
for fn in FS.index:
    a = ann[(ann.file_name == fn) & (ann.category_id == ds.MITOTIC)].reset_index(drop=True)
    im = images[images.file_name == fn].iloc[0]
    b = fs.FSConfig().patch_size // 2
    ix, iy = np.rint(a.cx).astype(int), np.rint(a.cy).astype(int)
    ok = a[(ix >= b) & (ix <= im.width-1-b) & (iy >= b) & (iy <= im.height-1-b)].reset_index(drop=True)
    pos = int(ok.index[ok.ann_id == seeds[fn]][0]); quantiles[fn] = round(pos/len(ok), 3)
qs = list(quantiles.values())
check("E6 seed positions are no longer one repeated quantile (was 0.75-0.85 for all 7)",
      (max(qs) - min(qs)) > 0.3, str(quantiles))

print("\n===== E7  optimal-assignment reporting =====")
need = {'greedy_matches','optimal_matches','greedy_only','optimal_only','n_detections_differing'}
d = pd.DataFrame({'cx':[1.,-3.],'cy':[0.,0.]}); g = pd.DataFrame({'cx':[0.,4.],'cy':[0.,0.]})
res = ev.optimal_assignment_disagreement(d, g, 5.0)
check("E7 new fields present", need <= set(res))
check("E7 counts the unmatched-in-one case the old //2 floored away",
      res['greedy_matches']==1 and res['optimal_matches']==2 and res['optimal_only']==2, str(res))

print("\n===== E8  max score after self-hit removal =====")
check("E8 both fields present", {'max_peak_score','max_detection_score'} <= set(m.columns))
check("E8 max_peak_score is the seed self-hit",
      bool(np.allclose(FS.max_peak_score, FS.seed_self_score)))
check("E8 max_detection_score is strictly lower and informative",
      bool((FS.max_detection_score < FS.max_peak_score).all()),
      FS[['max_peak_score','max_detection_score']].round(4).to_string())

print("\n===== E9  threshold_metrics deleted =====")
check("E9 function gone", not hasattr(ev, 'threshold_metrics'))
check("E9 no orphan columns", not ({'fp_composition_lookalike'} & set(m.columns)))

print("\n===== E10 pick_seed predicate matches the patch reader =====")
from midog_utils import template_match as tm
img = np.zeros((200,200), np.float32); mism = 0
for c in np.r_[np.arange(30.,40.,.05), np.arange(158.,170.,.05)]:
    df = pd.DataFrame({'cx':[c],'cy':[100.]})
    try: ex.pick_seed(df, np.random.default_rng(0), 36, (200,200)); e = True
    except ValueError: e = False
    if e != (tm.read_padded_patch(img, c, 100., 73) is not None): mism += 1
check("E10 zero predicate mismatches over the border range", mism == 0, f"{mism} mismatches")

print("\n===== E11 rank statistics =====")
pr = pd.read_csv('results/fs_score_probe.csv')
check("E11 mean and median both reported", {'mean_mitosis_rank','median_mitosis_rank'} <= set(pr.columns))
# AUC is stored rounded to 3 dp, worth +-0.0005*n_nuclei rank units, plus +-0.5 from
# int(round()) on the mean itself. Verified exact at full precision separately, below.
_tol = 0.0005*pr.n_nuclei + 1.0
check("E11 mean == (1-AUC)*n_nuclei, within CSV rounding",
      bool((abs(pr.mean_mitosis_rank - (1-pr.auc_mitosis_vs_nucleus)*pr.n_nuclei) <= _tol).all()),
      pr[['file_name','auc_mitosis_vs_nucleus','n_nuclei','mean_mitosis_rank','median_mitosis_rank']].to_string())
_a,_m,_md = ex._rank_stats(np.random.default_rng(0).normal(1,1,200),
                           np.random.default_rng(1).normal(0,1,7000))
check("E11 mean == (1-AUC)*n exactly at full precision",
      abs(_m - (1-_a)*7000) < 1e-9, f"mean={_m:.9f} vs {(1-_a)*7000:.9f}")
check("E11 mean and median are genuinely different numbers",
      bool((pr.mean_mitosis_rank > 1.5*pr.median_mitosis_rank).all()),
      f"mean {pr.mean_mitosis_rank.tolist()} vs median {pr.median_mitosis_rank.tolist()}")
check("E11 nucleus population excludes annotated objects",
      bool((pr.n_nuclei < pr.n_nuclei_all_blobs).all()),
      pr[['n_nuclei_all_blobs','n_nuclei']].to_string())

print("\n===== E12 select_domain_images returns real rows =====")
sel = ex.select_domain_images(images, ann)
check("E12 one row per domain", len(sel) == sel.tumor_type.nunique() == 7)
check("E12 every row exists in the source frame",
      all(((images.file_name==r.file_name)&(images.image_id==r.image_id)).any() for r in sel.itertuples()))

print("\n===== E14 fusion variants regenerated =====")
fv = pd.read_csv('results/fs_fusion_variants.csv')
check("E14 all four variants on both images", len(fv) == 8 and fv.variant.nunique() == 4)
check("E14 discrimination column present", 'discrimination' in fv.columns)
check("E14 shared nucleus population per image (identical n_nuclei across variants)",
      bool(fv.groupby('file_name').n_nuclei.nunique().eq(1).all()),
      fv.groupby('file_name').n_nuclei.unique().to_dict())

print("\n===== internal consistency (unchanged invariants) =====")
check("recall_at_k == topk_tp / n_gt_mitotic_eval",
      bool(np.allclose(FS.recall_at_k, FS.topk_tp/FS.n_gt_mitotic_eval)))
check("recall_at_k == mitotic_recall_at_k (one-to-one matching)",
      bool(np.allclose(FS.recall_at_k, FS.mitotic_recall_at_k)))
check("topk buckets sum to k", bool(((FS.topk_tp+FS.topk_fp_lookalike+FS.topk_fp_unannotated) == FS.k).all()))
check("all buckets sum to the detection count",
      bool(((FS.all_tp+FS.all_fp_lookalike+FS.all_fp_unannotated) == FS.all_n_detections).all()))
check("n_after_nms == detections + self-hits", bool((FS.n_after_nms == FS.n_detections_total + FS.n_self_hits).all()))
check("exactly one self-hit per run", bool((FS.n_self_hits == 1).all()))
check("detection CSV length matches n_detections_total",
      bool((det.groupby('file_name').size().reindex(FS.index) == FS.n_detections_total).all()))
check("every GT claimed at most once (one-to-one)",
      bool((det[det.matched_ann_id >= 0].groupby(['file_name','matched_ann_id']).size().max()) == 1))
check("seed never appears in its own eval GT",
      all(FS.loc[fn,'seed_ann_id'] not in set(det[(det.file_name==fn)].matched_ann_id) for fn in FS.index))
check("ROI areas all ~2 mm^2", bool(FS.roi_area_mm2.between(1.9,2.1).all()))
check("24 augmentations (not 72)", bool((FS.n_augmentations == 24).all()))

print("\n" + "="*60)
print("FAILURES:", fails if fails else "none")
print("="*60)
