"""v2 of the residual-prune feasibility probe. Three fixes over v1:

1. Cut at TP *quantiles* (retention 1.00/0.99/0.98/0.95), not at min(TP) -- a min over
   ~200 TPs is an order statistic pinned by one atypical mitosis, not a bound.
2. Deletion measured over the candidates ABOVE that cell's depth-to-100%, which is the
   only region a prune can pay in. The top-5000 band of v1 is not the working depth on
   any ROI that matters.
3. Discriminates the two mechanisms that both explain chromatin_od's 0/7 at read_100:
   (a) dark non-chromatin material floats to the top  -> h_frac should FALL with rank
       improving, i.e. the shallowest FPs are the least hematoxylin-pure;
   (b) deep mitoses are simply less dense           -> h_frac of deep TPs < shallow TPs.
"""
import numpy as np, pandas as pd, cv2, os
from skimage.color import rgb2hed

W, DARK_FRAC, SEED = 31, 0.10, 0
ROIS = ["548.tiff", "301.tiff", "094.tiff", "245.tiff", "402.tiff", "459.tiff", "246.tiff"]

pool = np.load("results/tm_recall_workload_pool_z.npz")
led = pd.read_csv("results/tm_recall_workload_tp_ledger.csv")

def hfrac(P):
    n, h, w, _ = P.shape
    hed = np.clip(rgb2hed(P.reshape(-1, 1, 3).astype(np.float32) / 255.0), 0, None).reshape(n, h*w, 3)
    tot = hed.sum(2); k = max(1, int(DARK_FRAC*h*w))
    idx = np.argpartition(tot, -k, axis=1)[:, -k:]
    sel = hed[np.arange(n)[:, None], idx, :]
    s = sel.sum((1, 2))
    return np.where(s > 0, sel[:, :, 0].sum(1)/np.maximum(s, 1e-9), np.nan)

rows, mech = [], []
for fn in ROIS:
    img = cv2.imread(f"images/extra_valid/{fn}", cv2.IMREAD_COLOR)[:, :, ::-1]
    H, Wd = img.shape[:2]; r = W//2
    cx = pool[f"{fn}|{SEED}|cx"].astype(int); cy = pool[f"{fn}|{SEED}|cy"].astype(int)
    tp_rank = np.sort(led[(led.file_name == fn) & (led.seed_index == SEED)]["rank"].to_numpy())
    depth100 = int(tp_rank.max()) + 1                      # the working region for "all TPs"
    # every candidate above depth100, subsampled if huge, plus every TP
    band = np.arange(depth100)
    if len(band) > 12000:
        rng = np.random.default_rng(0); band = np.sort(rng.choice(band, 12000, replace=False))
    want = np.union1d(band, tp_rank)
    ok = (cx[want] >= r) & (cy[want] >= r) & (cx[want] < Wd-r) & (cy[want] < H-r)
    want = want[ok]
    hf = hfrac(np.stack([img[cy[i]-r:cy[i]+r+1, cx[i]-r:cx[i]+r+1] for i in want]))
    is_tp = np.isin(want, tp_rank)
    tp_hf, fp_hf, fp_rank = hf[is_tp], hf[~is_tp], want[~is_tp]
    rec = dict(roi=fn, depth100=depth100, n_tp=int(is_tp.sum()), n_fp_band=int((~is_tp).sum()),
               tp_med=np.nanmedian(tp_hf), fp_med=np.nanmedian(fp_hf))
    for ret in (1.00, 0.99, 0.98, 0.95):
        cut = np.nanquantile(tp_hf, 1-ret)
        rec[f"del@{ret:.2f}"] = float(np.mean(fp_hf < cut))
    rows.append(rec)
    # mechanism (a): does h_frac rise with depth among FPs?  (b): deep TPs vs shallow TPs
    q = pd.qcut(fp_rank, 4, labels=False)
    mech.append(dict(roi=fn, **{f"fp_q{i}": float(np.nanmedian(fp_hf[q == i])) for i in range(4)},
                     tp_shallow=float(np.nanmedian(tp_hf[want[is_tp] < np.median(tp_rank)])),
                     tp_deep=float(np.nanmedian(tp_hf[want[is_tp] >= np.median(tp_rank)]))))
    print(rows[-1], flush=True); del img

d = pd.DataFrame(rows); m = pd.DataFrame(mech)
out = os.path.dirname(__file__)
d.to_csv(f"{out}/residual_probe2.csv", index=False); m.to_csv(f"{out}/residual_mech.csv", index=False)
print("\n=== deletion rate among candidates above depth-to-100%, at fixed TP retention ===")
print(d.round(4).to_string(index=False))
print("\n=== mechanism: median h_frac of FPs by rank quartile (q0 = best-ranked), and TPs by depth ===")
print(m.round(4).to_string(index=False))
