# Does raising the Otsu threshold un-merge multi-nucleus bbox-tightening components?

Date: 2026-09-03
Scope: `bbox_threshold_sweep.py` / `bbox_threshold_summary.py` / `threshold_diag.py`, 933 candidates
across all 18 currently-downloaded ROIs. Design reviewed and revised twice by an independent
subagent before running (see conversation; not reproduced here). Deliberate attempt to test/falsify
`2026-08-25-bbox-tuning-domain-mismatch.md`'s "domain property, not fixable by more threshold
tuning" conclusion at scale, not an oversight of it.

## What Experiment 5 and Experiment 6 actually were (they are not the same claim)

Re-read from `explore_dataset.ipynb` cells 15-24 and re-run against the current codebase, since the
2026-08-25 log's prose conflates two distinct mechanisms under "bbox tightening — dense domains":

- **Experiment 5** (`002.tiff` ann 20): the click's own pixel is *background* — it sits just
  outside its nucleus's Otsu mask. The question is a fallback policy question (grab the largest
  component in the crop when the click misses, vs. exclude). Re-run: in this one example the
  fallback happens to land on the correct, adjacent object (`debug_pics/ann20_otsu_fallback.png`)
  — the "would grab an unrelated structure" concern in the log is a *stated theoretical risk*, not
  something this example demonstrates. This mechanism is unrelated to thresholding and is **not**
  what this sweep tests.
- **Experiment 6** (cross-domain sweep, e.g. `245.tiff` ann 6245): the click's pixel *is*
  foreground, no fallback involved — but its connected component is the entire 51x51 window
  (area_frac=0.58, solidity=0.722), and **passes every current sanity check**
  (`min_area=50`, `max_area_frac=0.85`, `min_solidity=0.5`) despite visibly spanning a dense sheet
  of touching lymphocyte nuclei, not one object. This is what "neighbouring structures" actually
  refers to, and it's what this sweep tests.

## Design

Three arms, one shared, production-validated diagnostic pipeline (`threshold_diag.py` — `frac=0`
and `method="multiotsu"` both reproduce `seed_selection.tighten_box_otsu` byte-for-byte on spot
checks, ann 6245/20/10 included):

- **binary** — production default (`cv2.THRESH_OTSU`).
- **multiotsu** — production's existing opt-in (`skimage.filters.threshold_multiotsu`, 3-class, top
  class only).
- **headroom sweep** (new) — `threshold = otsu + frac*(255-otsu)` for `frac` in
  `{0.05, 0.10, 0.15, 0.20, 0.30}`. Not a flat pixel delta (would behave inconsistently across
  crops with different histogram shapes, the same problem the superseded 1D method had with a
  fixed relative threshold on busy tissue) — a fraction of the headroom between Otsu's own split
  and the crop's max (always 255, since `cv2.normalize` runs before thresholding in every arm).

Population: full `agreement_pool` + `border_filter` (`tm.PATCH_SIZE//2`=36px) across all 18
downloaded ROIs — 933 candidates, all 7 domains. Not the same 10-ROI set section 2b/7 used (that
population can no longer be exactly reconstructed — its scripts are gone from the repo), so these
percentages are internally comparable across this sweep's own arms, not directly against
`design_choices.md`'s numbers.

Split on binary's own diagnostic (`otsu_window` pinned to 51px everywhere):

| bucket | definition | n |
|---|---|---|
| `no_foreground` | click's pixel isn't foreground (Experiment-5 territory, out of scope) | 84 |
| `min_area_rejected` | component too small (unrelated to merging) | 4 |
| `flagged_rejected` | **already** rejected on `max_area_frac`/`min_solidity` | 3 |
| `flagged_accepted` | accepted, but `area_frac>=0.5` or `longer_side>=45px` | 206 |
| `control` | accepted, `area_frac<0.3` (clean) | 525 |
| `gap` | accepted, everything else | 111 |

A first draft of this design anchored "merged" only to binary's *accepted* output — a subagent
review caught that this structurally excludes the worst offenders (candidates `tighten_box_otsu`
already rejects for looking like a merge), which is why `flagged_rejected` exists as its own
population. In practice `min_solidity=0.5` almost never fires (n=3 across the whole dataset) — a
tightly-packed cluster of round nuclei stays convex enough to pass solidity even when it clearly
spans multiple objects, so the sanity checks are a weak filter for *this* failure mode specifically;
nearly all of it lives in `flagged_accepted`, not `flagged_rejected`.

## Results

**`flagged_accepted` (n=206, 93 lymphosarcoma + 50 mast-cell-tumor = 69% of it) — does raising the
threshold shrink these below the merge flag?**

| arm | resolved | still merged | excluded |
|---|---|---|---|
| frac=0.05 | 19.9% | 77.7% | 2.4% |
| frac=0.10 | 38.3% | 58.3% | 3.4% |
| frac=0.15 | 51.0% | 43.7% | 5.3% |
| frac=0.20 | 65.0% | 26.7% | 8.3% |
| frac=0.30 | 75.7% | 11.2% | 13.1% |
| multiotsu | 78.6% | 7.8% | 13.6% |

Monotone in `frac`, as expected. frac=0.30 lands within ~3 points of multiotsu on every metric —
raising the threshold most of the way to multiotsu's jump gets nearly all of multiotsu's effect.

**`control` (n=525, already clean) — collateral cost of the same intervention:**

| arm | collateral excluded (pooled) | ...in lymphosarcoma (n=51) |
|---|---|---|
| frac=0.05 | 1.0% | 2.0% |
| frac=0.15 | 3.2% | 11.8% |
| frac=0.30 | 5.1% | 13.7% |
| multiotsu | 4.0% | 13.7% |

Cost is real and concentrated exactly where the benefit is (dense domains) — frac=0.30's collateral
exclusion rate in lymphosarcoma matches multiotsu's exactly (13.7%), despite resolving slightly
fewer flagged cases (80.6% vs. 83.9%). There is no frac setting that gets multiotsu's benefit for
free; smaller fracs buy a smaller, genuinely dialable version of the same trade, not a free lunch.

**`flagged_rejected` (n=3, already-excluded-from-the-pool cases) — does raising the threshold
recover a usable box?** All 3 get a plausible-looking box at multiotsu / frac>=0.20 (visually
audited, `debug_pics/bbox_threshold_audit_recovered.png`). n is too small to generalize from, but
directionally positive and costs nothing extra to note.

## Visual audit (the geometry-only metrics above cannot show this on their own)

8 `flagged_accepted` "resolved" cases, sampled across lymphosarcoma/soft-tissue-sarcoma/
neuroendocrine/melanoma (`debug_pics/bbox_threshold_audit_resolved.png`): in every one, the binary
box visibly spans the annotated nucleus plus a neighbour or open tissue gap, and the multiotsu box
tightens cleanly onto the single densest, most plausible nucleus at the click. No "resolved-by-
geometry-but-visibly-wrong" case turned up in this sample (fragment, wrong neighbour captured) —
the geometric shrink tracked a real visual improvement every time sampled.

4 "still merged" cases (`debug_pics/bbox_threshold_audit_still_merged.png`), including the
`design_choices.md` §7 motivating example (350.tiff ann 18161): these do **not** look like
multiple touching nuclei that resisted separation. They look like single, genuinely large or
irregular chromatin masses (plausible anaphase/telophase figures) that the area/longer-side flag
mistakenly caught. **This means the ~79% multiotsu "resolved" rate above should not be read as
"79% were real merges that got fixed"** — the flag is a size proxy, not a merge detector; some
non-trivial share of `flagged_accepted` is not a merge at all, and those are exactly the cases
raising the threshold correctly leaves alone (a large object that's genuinely one object shouldn't
shrink to match a threshold-based expectation).

4 "excluded" cases (`debug_pics/bbox_threshold_audit_excluded.png`): genuinely dense, low-contrast
fields where multiotsu finds no valid component at all — real collateral loss, not spurious.

## Conclusion

The 2026-08-25 log's "domain property, not fixable by more threshold tuning" is **not fully
right** as a blanket claim: raising the threshold demonstrably separates a majority of the
window-saturating `flagged_accepted` population, visually confirmed as genuine single-nucleus
tightening in an 8-case sample with zero visibly-wrong resolutions. But it's not simply "more
tuning fixes it" either — the collateral cost is concentrated in exactly the same dense domains as
the benefit (up to ~14% of already-clean candidates lost in lymphosarcoma at multiotsu-level
aggressiveness), and a real share of the "still merged" population is not a merge at all, just a
legitimately large object the geometric flag can't distinguish from one.

**Not decided here**: whether to change the production default. This sweep only measures seed-
selection geometry, the same scope `design_choices.md` §7 stayed within — downstream
template-matching/discrimination quality (does a tightened-via-multiotsu template actually find
mitoses better or worse) remains completely untested, same gap §7b already named. If a change is
wanted, `frac=0.15` (headroom sweep, new in this analysis) is a candidate middle ground: resolves
about half of `flagged_accepted` (60.2% in mast-cell-tumor, 50.5% in lymphosarcoma) at roughly a
third of multiotsu's collateral cost in the domains that matter (11.8% vs. 13.7% lymphosarcoma
control exclusion) — an explicit opt-in worth prototyping end-to-end (matching a real
`find_and_suppress` run, not just geometry) before considering it further, not a recommendation to
adopt as the default.

Raw data: `results/bbox_threshold_sweep.csv` (933 rows x per-arm columns). Reproducible via
`python bbox_threshold_sweep.py && python bbox_threshold_summary.py`.
