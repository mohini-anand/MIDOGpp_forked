# `human_neuroendocrine_followup/` — an isolated, higher-n test of one domain result

## Why this folder exists

`production_hematoxylin_only/hem_rect_fill_vs_default_chromatin_od_by_domain.ipynb` sliced the
49-ROI `hem_rect_fill` vs `default_51` result by domain (7 domains x 7 ROIs) and found one domain,
`human neuroendocrine tumor`, with a small, broad (non-outlier), consistently negative direction —
but n = 7 per domain cannot resolve it: a back-of-envelope power calculation put the ROI count
needed at roughly 11-33 depending on K and how conservatively the multiple-comparison correction is
read. This folder pulls more `human neuroendocrine tumor` ROIs and re-runs the same two-arm test at
that higher n, **isolated from the main repo**: nothing here is written to `images/`,
`images/extra_valid/`, or `results/` at the repo root, and the analysis notebook lives in this
folder too, alongside its own `results/`.

## Inclusion criterion — relaxed, and why

The repo's usual "decision-grade" bar (`images/extra_valid/MANIFEST.md`) is `n_mitotic >= 15` and a
unanimous, border-filtered seed pool `>= 5`. At that bar there is exactly **one** unused valid ROI
left in the whole domain (`356.tiff`) — nowhere near enough. This folder instead requires only what
the experiment itself needs: **at least 3 mitotic annotations that are unanimous
(`agreement_pool`), survive the 36 px border filter, and pass `tighten_box_otsu`** — the same gate
`hem_rect_fill`'s template cut already applies, so nothing here is a click the main experiment
wouldn't also accept. Even that bar exhausts the domain's unanimous-tier candidates at 13 new ROIs
(`pull_rois.py`'s `FIDS`), so a further 7 come from `agreement_pool`'s contested 2-of-3 tier
(`FIDS_CONTESTED`) — tried only because the unanimous tier is exhausted, and every such ROI is
tagged `contested_2of3` throughout and reported separately, never silently merged with the rest.

## Contents

**22 ROIs used in the analysis** (7 existing + 15 new):

| tier | file(s) | n |
|---|---|---|
| `existing_7` — already used in the 49-ROI domain analysis, cloned here (not re-downloaded); the notebook reuses their exact recorded clicks, not redrawn | `353` `364` `400` `401` `402` `403` `404` | 7 |
| `unanimous_new` — new, unanimous-tier, `n_gate_valid >= 3` | `351` `356` `363` `365` `374` `376` `389` `394` `399` | 9 |
| `contested_2of3` — new, contested-tier backfill, `n_gate_valid >= 3` | `350` `361` `373` `379` `383` `392` | 6 |

**5 ROIs downloaded but excluded** (`n_gate_valid < 3`, cannot draw 3 distinct clicks): `352.tiff`
(2), `362.tiff` (2), `390.tiff` (1), `391.tiff` (2), `371.tiff` (1, contested tier). Left on disk in
`images/` for transparency; excluded by the notebook's own ROI list, not by deletion. Full counts
for every ROI tried (`n_mitotic`, `n_unanimous_pool`, `agreement_flagged`, `n_gate_valid`, `mpp`,
basic checks) are in `results/pull_verification.csv`.

**Basic checks** (image loads as HxWx3 uint8, pixel dims match the JSON, domain is
`human neuroendocrine tumor`, `mpp` in 0.22-0.26): **ALL 27/27 downloaded/cloned ROIs pass.**

## How this was built

`pull_rois.py`, run with the anaconda interpreter from this folder
(`/Users/mohinianand/anaconda3/bin/python3 pull_rois.py`). Existing ROIs are cloned via APFS
clonefile (`cp -c`, no extra disk); new ones stream from figshare via the file ids in
`../Setup.ipynb`'s download table, the same source `download_testing_set.py` used. That script and
`new_domain_roi_check.py` were **not** reused directly: both call
`midog_utils.dataset.load_slide_metadata` and `dataset.check_roi_scale`, neither of which exists in
the current `midog_utils` package (`production-cleanup-committed` — those two top-level scripts are
among the files that break by design post-cleanup). `pull_rois.py` re-implements the subset of
checks that still apply, from primitives `midog_utils` does still expose, plus the click-gate check
those older scripts never had a reason to compute.

## Caveat: iCloud eviction

This repo lives under iCloud Desktop & Documents sync, and the disk this ran on was within ~10-15
GiB of full. Both together mean a just-downloaded file can be evicted to an on-demand cloud
placeholder almost immediately — `ls`/`stat` still report the right apparent size, but a read can
hang or return a short/empty buffer while the real bytes are re-fetched. `pull_rois.py`'s
`materialize()` forces a `brctl download` and polls actual resident bytes (`st_blocks`) before every
read, with retries; this is why re-running the script is safe (already-verified files are skipped,
not re-fetched) but a first run may print `[retry ...]` lines. Anything reading these TIFFs later
(the analysis notebook included) should expect the same and may want the same guard if it starts
failing with `TimeoutError` or `tifffile.TiffFileError: not a TIFF file b''` — neither means the
data is bad; both have reproduced and resolved on retry here.
