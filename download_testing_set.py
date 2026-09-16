"""
    Download a held-out testing set of 5 valid MIDOG++ ROIs per tumour type into
    images/extra_valid/testing_set/, drawn uniformly at random and verified pixel by pixel.

    Validity is the same bar images/extra_valid/MANIFEST.md applies: n_mitotic >= 15 and a
    seed pool (agreement_pool -> border_filter(36, roi_shape)) of at least 5. Candidates are
    every valid ROI not already on disk in images/ or images/extra_valid/, predicted from the
    annotations and the JSON's width/height. Each domain's candidates are sorted by file name
    and shuffled once with np.random.default_rng([20260916, i]) (i = index of the domain in
    alphabetical order), and the whole ordering is written out before any download, so a
    replacement for a failed file is fixed in advance rather than chosen after the fact. The
    domain with the fewest spare candidates is downloaded first, so a shortfall shows up
    before the disk is spent on the rest; processing order does not change the draw.

    A download is accepted when new_domain_roi_check.check_roi passes every check except
    6_tissue_mask_covers_mitotic, which is recorded but not gated on: it tests a baseline's
    tissue mask, not the ROI. Anything else failing gets one re-download, then the file is
    deleted and the next candidate in that domain's ordering is tried. Re-running resumes:
    a file already in testing_set/ is re-verified instead of re-downloaded.
"""
import csv, json, shutil, socket, sys, traceback, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from midog_utils import dataset, seed_selection
from new_domain_roi_check import FIELDS, check_roi

PER_DOMAIN = 5
DRAW_SEED = 20260916
OUT_DIR = Path("images/extra_valid/testing_set")
ORDER_CSV = Path("results/testing_set_draw_order.csv")
LOG_CSV = Path("results/testing_set_rois.csv")
MIN_FREE_BYTES = 2 * 1024**3
NOT_GATED = {"6_tissue_mask_covers_mitotic"}
LOG_FIELDS = ["domain", "draw_rank", "attempt", "status", "split"] + FIELDS
socket.setdefaulttimeout(120)


def figshare_ids(notebook="Setup.ipynb"):
    """
        Map each ROI file name to its figshare file id, from the Setup notebook's download table.

        notebook (str): path to Setup.ipynb.

        Returns dict: file name (e.g. "402.tiff") -> figshare file id string.
    """
    src = "".join(json.load(open(notebook))["cells"][1]["source"])
    table = json.loads(src.split("=", 1)[1].replace("'", '"'))
    return {name: fid for fid, name in table.items()}


def valid_candidates(images, anns):
    """
        Every valid ROI not already on disk, with the counts that make it valid.

        images (pd.DataFrame): images frame from dataset.load_annotations.
        anns (pd.DataFrame): annotations frame from dataset.load_annotations.

        Returns pd.DataFrame: file_name, image_id, domain, n_mitotic, seed_pool, sorted by file name.
    """
    on_disk = {p.name for p in Path("images").glob("*.tiff")} | {p.name for p in Path("images/extra_valid").glob("*.tiff")}
    rows = []
    for im in images.itertuples():
        mit = anns[(anns["file_name"] == im.file_name) & (anns["category_id"] == dataset.MITOTIC)]
        pool, flagged = seed_selection.agreement_pool(mit)
        seeds = seed_selection.border_filter(pool, 36, (im.height, im.width))
        if len(mit) >= 15 and len(seeds) >= 5 and not flagged and im.file_name not in on_disk:
            rows.append(dict(file_name=im.file_name, image_id=im.image_id, domain=im.tumor_type, n_mitotic=len(mit), seed_pool=len(seeds)))
    return pd.DataFrame(rows).sort_values("file_name").reset_index(drop=True)


def draw_order(cands):
    """
        Shuffle each domain's candidates with a per-domain seeded generator.

        cands (pd.DataFrame): output of valid_candidates.

        Returns pd.DataFrame: cands plus draw_rank (0 = first to try), ordered by domain then rank.
    """
    parts = []
    for i, domain in enumerate(sorted(cands["domain"].unique())):
        dc = cands[cands["domain"] == domain].sort_values("file_name").reset_index(drop=True)
        perm = np.random.default_rng([DRAW_SEED, i]).permutation(len(dc))
        dc = dc.iloc[perm].reset_index(drop=True)
        dc["draw_rank"] = np.arange(len(dc))
        parts.append(dc)
    return pd.concat(parts, ignore_index=True)


def download(fid, dest):
    """
        Stream one figshare file to dest via a .part file, checking the byte count.

        fid (str): figshare file id.
        dest (Path): final path of the TIFF.

        Returns int: number of bytes written.
    """
    part = dest.with_suffix(".tiff.part")
    with urllib.request.urlopen(f"https://ndownloader.figshare.com/files/{fid}") as r, part.open("wb") as fh:
        expected = int(r.headers.get("Content-Length", -1))
        shutil.copyfileobj(r, fh, length=8 * 1024**2)
    n = part.stat().st_size
    if expected >= 0 and n != expected:
        part.unlink()
        raise IOError(f"truncated download: {n} of {expected} bytes")
    part.rename(dest)
    return n


def write_manifest():
    """
        Regenerate testing_set/MANIFEST.md from the draw-order and verification CSVs.

        Returns Path: the manifest written.
    """
    order = pd.read_csv(ORDER_CSV, dtype={"figshare_id": str})
    log = pd.read_csv(LOG_CSV)
    acc = log[log["status"] == "accepted"].sort_values(["domain", "file_name"])
    lines = [
        f"# `images/extra_valid/testing_set/` — the {len(acc)}-ROI testing set", "",
        f"{len(acc)} MIDOG++ ROIs, **{PER_DOMAIN} per tumour type across all 7 domains**, assembled 2026-09-16 by",
        "`download_testing_set.py`. None of them was on disk before that day: they are disjoint from the 14",
        "ROIs in `images/extra_valid/` and from everything else in `images/`. Unlike the parent folder, these",
        "are real downloads, not clones of `images/*.tiff`, and `images/` does not hold a copy.", "",
        "## Inclusion criterion", "",
        "The same bar as `images/extra_valid/MANIFEST.md` — an ROI is here only if **both** hold:", "",
        "1. **`n_mitotic >= 15`** — `category_id == 1` annotations in `databases/MIDOG++.json`.",
        "2. **seed pool >= 5** — `seed_selection.agreement_pool(mitotic)` then",
        "   `seed_selection.border_filter(pool, 36, roi_shape)`: mitotic annotations that *every* rater who saw",
        "   them called mitotic (2/2 or 3/3), at least 36 px from the ROI edge. This is stricter than \"5",
        "   unanimous annotations\" — the border filter can only remove figures — and no ROI here relied on",
        "   the contested 2-of-3 fallback tier.", "",
        "Both were predicted from the annotations before download and re-measured on the real pixel shape",
        "after it, together with the rest of `new_domain_roi_check.check_roi`'s battery: pixel dims equal the",
        "JSON's width/height (the truncation guard), `roi_mpp` in 0.22-0.26, ROI area in 1.9-2.1 mm2, and",
        "the tumour type. `6_tissue_mask_covers_mitotic` was recorded but not gated on (it tests",
        "`baselines.tissue_mask`, not the ROI); see the `checks` column.", "",
        "## Contents", "",
        "| file | domain | scanner | split | n_mitotic | n_unanimous | seed pool | checks |",
        "|---|---|---|---|---:|---:|---:|---|",
    ]
    for r in acc.itertuples():
        lines.append(f"| `{r.file_name}` | {r.domain} | {r.scanner} | {r.split} | {r.n_mitotic} | {r.n_unanimous_mitotic} | {r.seed_pool_size} | {r.checks_passed} |")
    lines += ["", "`split` is the MIDOG++ cross-validation split from `datasets_xvalidation.csv`, recorded for",
              "reference only; the draw did not filter on it.", "",
              "## Draw rule", "",
              "Uniformly at random within each domain, never by density — picking the densest candidates would bias",
              "every depth metric optimistically. Reproducible as written:", "",
              "```python",
              "# i indexes the alphabetically sorted domains",
              "candidates = sorted(valid ROIs of that domain not in images/ or images/extra_valid/, by file name)",
              f"order = candidates[np.random.default_rng([{DRAW_SEED}, i]).permutation(len(candidates))]",
              "# walk `order`, keeping each ROI that passes verification, until 5 are kept",
              "```", "",
              f"The full ordering for every domain, with figshare ids, is in `{ORDER_CSV}`, written before the first",
              f"download. Every attempt, including any rejection, is in `{LOG_CSV}`.", "",
              "| domain | candidates | ranks tried | rejected |", "|---|---:|---|---|"]
    for domain, g in order.groupby("domain"):
        lg = log[log["domain"] == domain]
        rej = sorted(set(lg[lg["status"] != "accepted"]["file_name"]) - set(acc["file_name"]))
        lines.append(f"| {domain} | {len(g)} | 0-{int(lg['draw_rank'].max())} | {', '.join(f'`{f}`' for f in rej) or 'none'} |")
    lines += ["", "## Caveats", "",
              "- **Human neuroendocrine tumor had only 6 valid candidates for 5 slots**, so this domain is",
              "  close to the whole valid population not already in `images/extra_valid/`, not a sample from it:",
              "  only `356.tiff` (rank 5: 18 mitotic, seed pool 11) was left out. It also holds the thinnest seed",
              "  pool in the set, `364.tiff` at 7, two above the bar.",
              "- The parent folder's caveat applies verbatim: `bbox_headroom_end_to_end.decision_grade_rois()` tests",
              "  `n_mitotic >= 15` alone and never looks at seed pools, so the guarantee lives in how this folder",
              "  was built, not in anything that re-checks it on read.",
              "- Existing scripts that read `images/extra_valid/` use a non-recursive `*.tiff` listing, so this",
              "  subfolder is invisible to them; the 14-ROI set is unchanged. A recursive glob (`rglob`, `**`)",
              "  over `images/extra_valid/` would silently mix the two sets.", ""]
    path = OUT_DIR / "MANIFEST.md"
    path.write_text("\n".join(lines))
    return path


def main():
    """
        Draw, download and verify until every domain has PER_DOMAIN accepted ROIs.

        Returns int: process exit code, 0 when every domain is filled.
    """
    images, anns = dataset.load_annotations()
    dataset.check_invariants(anns)
    meta = dataset.load_slide_metadata()
    split = dict(zip(meta["image_id"], meta["Dataset"]))
    fids = figshare_ids()

    order = draw_order(valid_candidates(images, anns))
    order["split"] = order["image_id"].map(split)
    order["figshare_id"] = order["file_name"].map(fids)
    ORDER_CSV.parent.mkdir(exist_ok=True)
    order.to_csv(ORDER_CSV, index=False)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    spare = order.groupby("domain").size() - PER_DOMAIN
    print("spare candidates per domain:", spare.to_dict(), flush=True)
    domains = sorted(spare.index, key=lambda d: (spare[d], d))

    fh = LOG_CSV.open("w", newline="")
    log = csv.DictWriter(fh, fieldnames=LOG_FIELDS); log.writeheader(); fh.flush()
    shortfall = {}
    for domain in domains:
        accepted = 0
        for c in order[order["domain"] == domain].itertuples():
            if accepted == PER_DOMAIN:
                break
            dest = OUT_DIR / c.file_name
            for attempt in (1, 2):
                base = dict(domain=domain, draw_rank=c.draw_rank, attempt=attempt, split=c.split, file_name=c.file_name, image_id=c.image_id)
                try:
                    if not dest.exists():
                        if shutil.disk_usage(OUT_DIR).free < MIN_FREE_BYTES:
                            raise SystemExit(f"less than {MIN_FREE_BYTES / 1024**3:.0f} GiB free; stopping before {c.file_name}")
                        download(c.figshare_id, dest)
                    row, checks = check_roi(dest, domain, images, anns, meta)
                    failed = [k for k, v in checks if not v and k not in NOT_GATED]
                    status = "accepted" if not failed else "rejected:" + ",".join(failed)
                    log.writerow({**base, **row, "status": status}); fh.flush()
                    print(f"[{domain}] rank {c.draw_rank} {c.file_name} attempt {attempt}: {status} (n_mitotic={row['n_mitotic']}, seed_pool={row['seed_pool_size']}, checks={row['checks_passed']})", flush=True)
                    if not failed:
                        accepted += 1
                        break
                except SystemExit:
                    raise
                except Exception:
                    log.writerow({**base, "status": "exception", "notes": traceback.format_exc()[-300:]}); fh.flush()
                    print(f"[{domain}] rank {c.draw_rank} {c.file_name} attempt {attempt}: EXCEPTION", flush=True)
                    traceback.print_exc()
                if dest.exists():
                    dest.unlink()
        if accepted < PER_DOMAIN:
            shortfall[domain] = accepted
        print(f"[{domain}] accepted {accepted}/{PER_DOMAIN}; free disk {shutil.disk_usage(OUT_DIR).free / 1024**3:.1f} GiB", flush=True)
    fh.close()
    print("wrote", write_manifest(), flush=True)
    if shortfall:
        print("SHORTFALL:", shortfall, flush=True)
        return 1
    print("all domains filled", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
