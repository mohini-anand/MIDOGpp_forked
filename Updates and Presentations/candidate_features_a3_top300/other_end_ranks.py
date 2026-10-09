"""
Where the mitoses at a feature's "other end" sit in the search's own ranking.

For each of the deck's six features, a click's recorded pool is ordered by the feature in the mitosis-like direction. The
"other end" is the bottom 20 % of that order (tie-averaged percentile). The table reports, for the mitoses there, their
rank in the click's pool by search score and by od (1 = top), and the share in the score top 10.

Post hoc, in-sample, on the V2 candidate table (294 clicks, 7,664 pool mitoses). Prints a table; writes nothing.

Run with the Miniconda 3.11 interpreter: ~/anaconda3/bin/python other_end_ranks.py
"""
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OTHER_END = 0.20
FEATS = [('Ragged outline', 'dark_concavities_per_um', +1), ('Convexity', 'chr_solidity', -1), ('Elongation', 'chr_elongation', +1),
         ('Core darkness', 'h_chr_mean', +1), ('White surround', 'white_frac_gray215_61', -1), ('Non-pale rim', 'ring1_frac_above_tnuc', +1)]


def other_end_row(P, name, col, sign):
    """
        One feature's row: the mitoses in the bottom 20 % of each click's oriented ranking and where search puts them.

        P (pd.DataFrame): the recorded-pool, accepted rows, with click, score_rank, od_rank and is_mitosis.
        name (str): the plain feature name.
        col (str): the feature column.
        sign (int): +1 when high values look like a mitosis, -1 when low values do.

        Returns dict: the row of the table.
    """
    pct = (P[col] * sign).groupby(P['click']).rank(pct=True, method='average')
    m = P[P.is_mitosis & P[col].notna()]
    oth = m[pct.loc[m.index] <= OTHER_END]
    return {'feature': name, 'mitoses at the other end': len(oth) / len(m), 'median score rank (those)': oth.score_rank.median(),
            'share in score top 10 (those)': (oth.score_rank <= 10).mean(), 'median od rank (those)': oth.od_rank.median(),
            'median score rank (all mitoses)': m.score_rank.median(), 'share in score top 10 (all mitoses)': (m.score_rank <= 10).mean()}


if __name__ == '__main__':
    use = ['file_name', 'seed_index', 'in_recorded_pool', 'accepted', 'is_mitosis', 'od', 'score'] + [c for _, c, _ in FEATS]
    C = pd.read_csv(REPO / 'results/candidate_features_a3_top300_candidates.csv', usecols=use, low_memory=False)
    P = C[C.in_recorded_pool & C.accepted].copy()
    P['click'] = P.file_name + '|' + P.seed_index.astype(str)
    P['od_rank'] = P.groupby('click').od.rank(ascending=False, method='first')
    P['score_rank'] = P.groupby('click').score.rank(ascending=False, method='first')
    print(len(P), 'pool rows,', int(P.is_mitosis.sum()), 'mitoses')
    print(pd.DataFrame([other_end_row(P, *f) for f in FEATS]).round(3).to_string(index=False))
