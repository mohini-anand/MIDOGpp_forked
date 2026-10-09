# Candidate features: what helps (deck, 2026-10-09)

A six-slide deck on the V2 candidate-feature run (`production_hematoxylin_only/multiotsu_pad_arms/candidate_features_a3_top300*`)
and the earlier re-rank run (`candrefine_rerank_mo_pad61_select.ipynb`). The deck is a claude.ai artifact (private until
shared): https://claude.ai/artifact/NLcRKHzs1d6kcupwxovqbv

| file | what it is |
|---|---|
| `deck/project/` | the deck's source: `deck.json` and one HTML file per slide. Images are referenced by uploaded-asset id, so it renders only inside the artifact |
| `make_figs.py` | writes the 14 figures in `figs/` (schematics, ROI 128 example strips, the two charts). `figs/` is not tracked; rerun the script |
| `other_end_ranks.py` | prints the table of where mitoses at a feature's "other end" sit in the search ranking (not in the deck) |

Everything is post hoc and in-sample on the 49 ROIs. The od reference AUC (0.77) and the other-end table are computed by these
scripts; the AUC-at-matched-darkness values come from `results/candidate_features_a3_top300_audit_stratified.csv` and the
precision bars from the earlier notebook's printed P@K table.

The V2 features have not been run through a ranker, so the deck makes no claim about their precision over od.
