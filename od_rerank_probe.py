"""SUPERSEDED exploration -- kept for provenance.

This probe chose the statistic (mean of the darkest 10% of the window, over the Otsu
component alternative). It was run against `channels.to_hematoxylin`, which saturates dense
chromatin at its 255 ceiling; the pipeline now uses `chromatin.hematoxylin_od` (unclipped).
The *choice* it drove is unchanged and was re-confirmed post-fix, but the AUCs printed here
are pre-fix and are NOT the authoritative numbers -- see
`Research Logs/2026-08-31-chromatin-density-rerank.md` and `results/od_workload_ab.csv`.

Can a single post-hoc scalar -- mean hematoxylin OD under each detection -- re-rank
the existing detection list to cut FPs while keeping TPs?

Two variants:
  comp = mean OD over the tighten_box_otsu component  (needs the object gate)
  disc = mean OD over a fixed 21px central disc       (no gate, no TP cost)
Compared against the pipeline's own correlation score as a ranker.
"""
import numpy as np, pandas as pd
from midog_utils import dataset as ds, channels as ch, seed_selection as sel, template_match as tm

R = 10
yy, xx = np.mgrid[-R:R+1, -R:R+1]
DISC = (yy**2 + xx**2) <= R*R

def auc(pos, neg):
    if not len(pos) or not len(neg): return np.nan
    return float((pos[:, None] > neg[None, :]).mean() + 0.5*(pos[:, None] == neg[None, :]).mean())

det = pd.read_csv("results/fs_simple_detections.csv")
met = pd.read_csv("results/fs_simple_metrics.csv")
ngt = dict(zip(met.file_name, met.n_gt_mitotic_eval))

out, curves = [], []
for fn in ["301.tiff", "405.tiff", "002.tiff", "506.tiff", "350.tiff"]:
    hem = ch.to_hematoxylin(ds.load_roi(f"images/{fn}"))
    d = det[det.file_name == fn]
    tp_all = d[d.bucket == "human_correct_label"]
    un_all = d[d.bucket == "non_human_findings"]
    un = un_all.sample(min(2000, len(un_all)), random_state=0)
    sub = pd.concat([tp_all, d[d.bucket == "human_rejected_label"], un])

    comp, disc, keep = [], [], []
    for cx, cy in zip(sub.cx.values, sub.cy.values):
        p = tm.read_padded_patch(hem, cx, cy, tm.BASE_SIZE)
        if p is None:
            comp.append(np.nan); disc.append(np.nan); keep.append(False); continue
        c = p.shape[0] // 2
        disc.append(float(p[c-R:c+R+1, c-R:c+R+1][DISC].mean()))
        b = sel.tighten_box_otsu(p)
        if b is None:
            comp.append(np.nan); keep.append(False)
        else:
            y0, y1, x0, x1 = b
            comp.append(float(p[y0:y1, x0:x1].mean())); keep.append(True)
    sub = sub.assign(od_comp=comp, od_disc=disc, gated=keep)

    T = sub[sub.bucket == "human_correct_label"]; U = sub[sub.bucket == "non_human_findings"]
    frac_un = len(un) / len(un_all)
    out.append(dict(image=fn, n_mit=ngt[fn], n_det=len(d), n_tp=len(tp_all), n_unann=len(un_all),
        auc_ncc_score=round(auc(T.score.values, U.score.values), 3),
        auc_od_disc=round(auc(T.od_disc.dropna().values, U.od_disc.dropna().values), 3),
        auc_od_comp=round(auc(T.od_comp.dropna().values, U.od_comp.dropna().values), 3)))

    # workload curve: threshold od_disc at each FP percentile, no object gate
    for q in (50, 75, 90, 95):
        t = np.nanpercentile(U.od_disc.values, q)
        tp_ret = float(np.nanmean(T.od_disc.values >= t))
        fp_kept = float(np.nanmean(U.od_disc.values >= t))
        curves.append(dict(image=fn, cut=f"FP p{q}",
            tp_retained=round(tp_ret, 3),
            mitotic_recall_end_to_end=round(tp_ret * len(tp_all) / ngt[fn], 3),
            candidates_left=int(round(fp_kept*len(un_all) + tp_ret*len(tp_all))),
            was=len(d),
            per_mitosis_found=round((fp_kept*len(un_all)) / max(tp_ret*len(tp_all), 1), 1)))

pd.set_option("display.width", 250)
print(pd.DataFrame(out).to_string(index=False))
print()
print(pd.DataFrame(curves).to_string(index=False))
