"""Summarise the per-ROI mechanism measurements (reads mech/, writes mech/summary.csv)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm, mannwhitneyu

SP = Path('/private/tmp/claude-501/-Users-mohinianand-Desktop-AnnotateDx-MIDOGpp-forked/42eaa83a-927f-4eb9-93e3-6ebb100b4d9d/scratchpad')
M = Path(sys.argv[1]) if len(sys.argv) > 1 else SP / 'mech'
pr = pd.read_csv(SP / 'per_roi_recent2.csv', index_col=0)
pts = pd.read_csv(M / 'points.csv')
meta = pd.read_csv(M / 'roi_meta.csv').set_index('file_name')


def auc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    return mannwhitneyu(pos, neg).statistic / (len(pos) * len(neg))


rows = []
for fn in meta.index:
    b = pd.read_csv(M / f'blobs_{fn}.csv.gz')
    bc = pd.read_csv(M / f'blobcolour_{fn}.csv.gz')
    P = pts[pts.file_name == fn]
    mit = P[P.kind == 'gt_mit']
    t10 = P[P.kind == 'top10']
    fp = t10[t10.bucket != 'human_correct_label']
    neg = b[b.label == 'neg']
    mb = b[b.label == 'mit']
    # one value per mitosis: its darkest matched blob
    mb1 = mb.sort_values('od51', ascending=False).drop_duplicates('mit_id')
    n_mit = meta.at[fn, 'n_mit']
    N = len(b)
    # whole-ROI od51 ranking of blobs: unique mitoses in the first 10 / 20 blobs
    srt = b.sort_values('od51', ascending=False)
    def uniq_hits(k):
        s = srt.head(k)
        return s[s.label == 'mit'].mit_id.nunique()
    med_mit = mb1.od51.median()
    # idealised equal-variance requirement for E[TP in top 10] = 5 with M = n_mit - 1 positives among N
    Mpos = n_mit - 1
    t = norm.isf(5 / N)
    dreq = t - norm.isf(min(5 / Mpos, 0.999))
    # colour: pigment-like = darkest-pixel blue/red OD ratio and DAB above the ROI's mitosis 95th pct
    br95, d95 = mit.br.quantile(0.95), mit.D.quantile(0.95)
    pig = (fp.br > br95) & (fp.D > d95)
    rows.append(dict(
        file_name=fn, N_blobs=N, n_mit=n_mit, blob_cov=meta.at[fn, 'blob_cov_mit'],
        prev_per_1000=1000 * n_mit / N,
        od51_auc_mit_vs_neg=auc(mb1.od51, neg.od51), od21_auc=auc(mb.sort_values('od21', ascending=False).drop_duplicates('mit_id').od21, neg.od21),
        od51_auc_mit_vs_la=auc(mb1.od51, b[b.label == 'la'].od51),
        mit_od_pctile=float((b.od51.to_numpy()[:, None] < mb1.od51.to_numpy()[None, :]).mean(0).mean()) if len(mb1) else np.nan,
        n_neg_darker_than_median_mit=int((neg.od51 > med_mit).sum()),
        roi_od_P10=uniq_hits(10) / 10, roi_od_P20=uniq_hits(20) / 20,
        d_req=dreq, auc_req=norm.cdf(dreq / np.sqrt(2)),
        fp_n=len(fp), fp_white=float((fp.white_frac_61 > 0.20).mean()), fp_edge=float((fp.edge_px < 35).mean()),
        fp_la=float((fp.bucket == 'human_rejected_label').mean()),
        fp_la_contested=float(((fp.bucket == 'human_rejected_label') & fp.matched_ann_id.isin(P[(P.kind == 'gt_la') & (P.n_mit_votes > 0)].ann_id)).mean()),
        fp_pigment=float(pig.mean()),
        fp_od21_over_mit=fp.od21.median() / mit.od21.median(),
        fp_area_over_mit=fp.blob_area.median() / mit.blob_area.median(),
        mit_area_med=mit.blob_area.median(),
        mit_od21_pct_of_negblobs=float((neg.od21.to_numpy() < mit.od21.median()).mean()),
        neg_od21_p99_over_mit_med=neg.od21.quantile(0.99) / mit.od21.median(),
        mit_br_med=mit.br.median(), fp_br_med=fp.br.median(), mit_D_med=mit.D.median(), fp_D_med=fp.D.median(),
        blob_D_p99=bc.D.quantile(0.99),
    ))
S = pd.DataFrame(rows).set_index('file_name')
S = pr[['domain', 'A3_tm_P10', 'A3_ch_P10', 'A3_tm_P20', 'A3_ch_P20', 'C1_P10', 'C1_P20', 'mit_unan', 'la_per_mit', 'white_share']].join(S)
S.to_csv(M / 'summary.csv')
pd.set_option('display.width', 300); pd.set_option('display.max_columns', 60)
S = S.sort_values('A3_tm_P10')
print(S[['domain', 'A3_tm_P10', 'A3_ch_P10', 'C1_P10', 'n_mit', 'N_blobs', 'prev_per_1000', 'blob_cov', 'od51_auc_mit_vs_neg', 'od21_auc', 'od51_auc_mit_vs_la', 'mit_od_pctile', 'n_neg_darker_than_median_mit', 'roi_od_P10', 'auc_req']].round(3).to_string())
print(S[['A3_tm_P10', 'fp_n', 'fp_white', 'fp_edge', 'fp_la', 'fp_la_contested', 'fp_pigment', 'fp_od21_over_mit', 'fp_area_over_mit', 'mit_area_med', 'neg_od21_p99_over_mit_med', 'mit_br_med', 'fp_br_med', 'mit_D_med', 'fp_D_med', 'mit_unan']].round(3).to_string())
num = S.select_dtypes('number')
print(num.corr('spearman')[['A3_tm_P10', 'A3_ch_P10', 'C1_P10']].round(2).to_string())
