"""Per-ROI mechanism measurements for the low-precision ROI question (2026-10-02).

For every one of the 49 ROIs:
  * the whole-ROI nucleus population from `baselines.nucleus_blobs` (no click, no TM);
  * od51 (= chromatin.chromatin_density, darkest 10% of hematoxylin OD in 51 px) and od21 at
    every blob, every GT mitosis and look-alike, and every A3 top-10 detection;
  * colour of the darkest pixels (per-channel RGB OD, blue/red OD ratio, rgb2hed E and D);
  * blob area.
Writes points + per-ROI summaries to the scratchpad. Read-only on the repo.
"""
import sys, json, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from skimage.color import rgb2hed

REPO = Path('/Users/mohinianand/Desktop/AnnotateDx/MIDOGpp_forked')
sys.path.insert(0, str(REPO))
from midog_utils_full import dataset as ds
from midog_utils_full.baselines import nucleus_blobs
from midog_utils_full.chromatin import hematoxylin_od

OUT = Path(sys.argv[1])
ONLY = sys.argv[2].split(',') if len(sys.argv) > 2 else None

j = json.load(open(REPO / 'databases/MIDOG++.json'))
imgs = {i['id']: i for i in j['images']}
ann = pd.DataFrame(j['annotations'])
ann['file_name'] = ann.image_id.map(lambda i: imgs[i]['file_name'])
ann['cx'] = ann.bbox.map(lambda b: (b[0] + b[2]) / 2)
ann['cy'] = ann.bbox.map(lambda b: (b[1] + b[3]) / 2)
ann['n_mit_votes'] = ann['labels'].map(lambda l: sum(1 for x in l if x == 1))

top = pd.read_csv(REPO / 'results/tm_border_floor_49roi_6seed_top30.csv')
top = top[(top.arm == 'A3') & (top['rank'] < 10)]


def find_path(fn):
    for d in ('images/extra_valid', 'images/extra_valid/testing_set'):
        p = REPO / d / fn
        if p.exists():
            return p
    raise FileNotFoundError(fn)


def darkest_mean(img, cx, cy, win, frac=0.10):
    h = win // 2
    H, W = img.shape[:2]
    ix, iy = int(round(cx)), int(round(cy))
    p = img[max(0, iy - h):min(H, iy + h + 1), max(0, ix - h):min(W, ix + h + 1)].ravel()
    k = max(1, int(frac * p.size))
    return float(np.partition(p, -k)[-k:].mean())


def colour_feats(rgb, hem, cx, cy, win=21, frac=0.10):
    """RGB optical density and HED E/D of the darkest-by-hematoxylin `frac` of a win x win core."""
    h = win // 2
    H, W = hem.shape
    ix, iy = int(round(cx)), int(round(cy))
    sl = (slice(max(0, iy - h), min(H, iy + h + 1)), slice(max(0, ix - h), min(W, ix + h + 1)))
    hp = hem[sl].ravel()
    k = max(1, int(frac * hp.size))
    idx = np.argpartition(hp, -k)[-k:]
    px = rgb[sl].reshape(-1, 3)[idx].astype(np.float32)
    od = -np.log10(np.maximum(px, 1.0) / 255.0)
    hed = rgb2hed(px[None, :, :] / 255.0)[0]
    odm = od.mean(0)
    return dict(odR=odm[0], odG=odm[1], odB=odm[2], br=odm[2] / max(odm[0], 1e-6),
                E=float(hed[:, 1].mean()), D=float(hed[:, 2].mean()), Hc=float(hed[:, 0].mean()),
                gray_min=float(px.mean(1).min()))


def white_frac(gray, cx, cy, win=61, thr=215):
    h = win // 2
    H, W = gray.shape
    ix, iy = int(round(cx)), int(round(cy))
    p = gray[max(0, iy - h):min(H, iy + h + 1), max(0, ix - h):min(W, ix + h + 1)]
    return float((p > thr).mean())


files = sorted(top.file_name.unique())
if ONLY:
    files = [f for f in files if f in ONLY]
pts_all, roi_rows = [], []
for fn in files:
    t0 = time.time()
    path = find_path(fn)
    rgb = ds.load_roi(path)
    mpp = ds.roi_mpp(path)
    r = 7.5 / mpp
    hem = hematoxylin_od(rgb)
    gray = (0.299 * rgb[..., 0] + 0.587 * rgb[..., 1] + 0.114 * rgb[..., 2]).astype(np.float32)
    blobs = nucleus_blobs(rgb)
    a = ann[ann.file_name == fn]
    mit = a[a.category_id == 1]
    la = a[a.category_id == 2]
    # label blobs: nearest annotation within r
    lab = np.array(['neg'] * len(blobs), dtype=object)
    bxy = blobs[['cx', 'cy']].to_numpy()
    if len(la):
        d, _ = cKDTree(la[['cx', 'cy']].to_numpy()).query(bxy)
        lab[d <= r] = 'la'
    d_m, i_m = cKDTree(mit[['cx', 'cy']].to_numpy()).query(bxy)
    lab[d_m <= r] = 'mit'
    blobs = blobs.assign(label=lab, mit_id=np.where(d_m <= r, mit.id.to_numpy()[i_m], -1))
    # mitosis coverage by blobs
    dm2, _ = cKDTree(bxy).query(mit[['cx', 'cy']].to_numpy())
    cov = float((dm2 <= r).mean())
    btree = cKDTree(bxy)

    def pts(df, kind, extra=None):
        rows = []
        for k, (cx, cy) in enumerate(zip(df.cx.to_numpy(), df.cy.to_numpy())):
            row = dict(file_name=fn, kind=kind, cx=cx, cy=cy,
                       od51=darkest_mean(hem, cx, cy, 51), od21=darkest_mean(hem, cx, cy, 21),
                       white61=white_frac(gray, cx, cy))
            row.update(colour_feats(rgb, hem, cx, cy))
            dd, ii = btree.query([cx, cy])
            row['blob_area'] = float(blobs.area.iat[ii]) if dd <= r else np.nan
            if extra is not None:
                for c in extra:
                    row[c] = df[c].iat[k]
            rows.append(row)
        return rows

    rows = []
    rows += pts(mit.assign(ann_id=mit.id, n_mit_votes=mit.n_mit_votes), 'gt_mit', ['ann_id', 'n_mit_votes'])
    rows += pts(la.assign(ann_id=la.id, n_mit_votes=la.n_mit_votes), 'gt_la', ['ann_id', 'n_mit_votes'])
    tt = top[top.file_name == fn]
    rows += pts(tt.assign(), 'top10', ['rank', 'seed_index', 'rank_key', 'bucket', 'matched_ann_id', 'edge_px', 'white_frac_61', 'score'])
    P = pd.DataFrame(rows)
    pts_all.append(P)
    # blob population: od51/od21 only (colour on a 3000-blob random subsample for speed)
    blobs['od51'] = [darkest_mean(hem, x, y, 51) for x, y in bxy]
    blobs['od21'] = [darkest_mean(hem, x, y, 21) for x, y in bxy]
    rng = np.random.default_rng(0)
    sub = rng.choice(len(blobs), size=min(3000, len(blobs)), replace=False)
    cf = pd.DataFrame([colour_feats(rgb, hem, bxy[i, 0], bxy[i, 1]) for i in sub])
    cf['white61'] = [white_frac(gray, bxy[i, 0], bxy[i, 1]) for i in sub]
    cf['label'] = blobs.label.to_numpy()[sub]
    cf['od51'] = blobs.od51.to_numpy()[sub]
    cf['area'] = blobs.area.to_numpy()[sub]
    cf['file_name'] = fn
    blobs.assign(file_name=fn).to_parquet(OUT / f'blobs_{fn}.parquet') if False else blobs.assign(file_name=fn).to_csv(OUT / f'blobs_{fn}.csv.gz', index=False)
    cf.to_csv(OUT / f'blobcolour_{fn}.csv.gz', index=False)
    roi_rows.append(dict(file_name=fn, mpp=mpp, n_blobs=len(blobs), n_mit=len(mit), n_la=len(la),
                         blob_cov_mit=cov, n_blob_mit=int((lab == 'mit').sum()), n_blob_la=int((lab == 'la').sum())))
    print(f'{fn}: {len(blobs)} blobs, cov {cov:.2f}, {time.time() - t0:.0f}s', flush=True)
    del rgb, hem, gray

pd.concat(pts_all).to_csv(OUT / 'points.csv', index=False)
pd.DataFrame(roi_rows).to_csv(OUT / 'roi_meta.csv', index=False)
print('done')
