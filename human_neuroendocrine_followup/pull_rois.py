"""
    Pull the human_neuroendocrine_tumor ROIs needed for an isolated, higher-n follow-up of the
    domain-level hem_rect_fill vs default_51 result (see
    ../production_hematoxylin_only/hem_rect_fill_vs_default_chromatin_od_by_domain.ipynb). Kept in
    this folder, separate from images/, images/extra_valid/ and results/, so this exploratory pull
    never touches the main repo's data.

    Relaxes the repo's usual "decision-grade" bar (n_mitotic >= 15, seed pool >= 5, from
    images/extra_valid/MANIFEST.md) to just "at least 3 unanimous, border-filtered,
    Otsu-tightenable mitotic annotations" -- enough to draw 3 distinct clicks. At the strict bar
    there is exactly one unused valid ROI left in the whole domain (356.tiff); the relaxed bar
    admits 13.

    new_domain_roi_check.check_roi and download_testing_set.py are not reused: both call
    midog_utils.dataset.load_slide_metadata and dataset.check_roi_scale, neither of which exists in
    the current midog_utils package (production-cleanup-committed; those two top-level scripts are
    among the files that break by design post-cleanup). This script re-implements the subset of
    checks that still apply, from primitives midog_utils does still expose.
"""
import shutil
import socket
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, '..')
from midog_utils import channels as ch, dataset as ds, seed_selection as ss, template_match as tm

socket.setdefaulttimeout(120)
OUT_DIR = Path('images')
DOMAIN = 'human neuroendocrine tumor'
MIN_FREE_BYTES = 2 * 1024 ** 3

# the 7 ROIs already used in the 49-ROI domain analysis -- cloned (not re-downloaded), and their
# recorded clicks are reused verbatim by the notebook, not redrawn
EXISTING = {
    '353.tiff': '../images/extra_valid/testing_set/353.tiff',
    '364.tiff': '../images/extra_valid/testing_set/364.tiff',
    '400.tiff': '../images/extra_valid/testing_set/400.tiff',
    '401.tiff': '../images/extra_valid/testing_set/401.tiff',
    '402.tiff': '../images/extra_valid/402.tiff',
    '403.tiff': '../images/extra_valid/403.tiff',
    '404.tiff': '../images/extra_valid/testing_set/404.tiff',
}
# already downloaded elsewhere in the repo (images/), not part of the 49-ROI set -- just needs cloning
ALREADY_ON_DISK = {'351.tiff': '../images/351.tiff'}
# figshare file ids (../Setup.ipynb's download table) for the rest -- none previously downloaded anywhere in the repo
FIDS = {
    '352.tiff': '40282987', '356.tiff': '40283002', '362.tiff': '40283083', '363.tiff': '40283032',
    '365.tiff': '40283026', '374.tiff': '40283050', '376.tiff': '40283053', '389.tiff': '40283080',
    '390.tiff': '40283182', '391.tiff': '40283098', '394.tiff': '40283236', '399.tiff': '40283122',
}
# backfill: the unanimous-tier candidate pool above is exhausted (no unflagged image in the domain has
# a border-filtered pool >= 3 left after the 13 above), so these fall back to agreement_pool's contested
# 2-of-3 tier (flagged=True) -- kept separate and reported separately, never silently merged into FIDS
ALREADY_ON_DISK_CONTESTED = {'350.tiff': '../images/350.tiff'}
FIDS_CONTESTED = {
    '392.tiff': '40283200', '361.tiff': '40283020', '379.tiff': '40283146',
    '383.tiff': '40283167', '371.tiff': '40283113', '373.tiff': '40283131',
}


def clone(src, dest):
    """
        A byte-identical copy of src at dest via APFS clonefile (cp -c), costing no extra disk.

        src (str): source path.
        dest (Path): destination path.
    """
    subprocess.run(['cp', '-c', src, str(dest)], check=True)


def download(fid, dest):
    """
        Stream one figshare file to dest via a .part file, checking the byte count.

        fid (str): figshare file id.
        dest (Path): final path of the TIFF.

        Returns int: number of bytes written.
    """
    part = dest.with_suffix('.tiff.part')
    with urllib.request.urlopen(f'https://ndownloader.figshare.com/files/{fid}') as r, part.open('wb') as fh:
        expected = int(r.headers.get('Content-Length', -1))
        shutil.copyfileobj(r, fh, length=8 * 1024 ** 2)
    n = part.stat().st_size
    if expected >= 0 and n != expected:
        part.unlink()
        raise IOError(f'truncated download: {n} of {expected} bytes')
    part.rename(dest)
    return n


def verify(path, image_id, images, anns):
    """
        The subset of new_domain_roi_check.check_roi's checks that still apply against the current
        midog_utils package, plus the actual click-gate pool size this experiment needs.

        path (Path): the cloned/downloaded TIFF.
        image_id (int): the MIDOG++ image id this file name maps to.
        images (pd.DataFrame): images frame from dataset.load_annotations.
        anns (pd.DataFrame): annotations frame from dataset.load_annotations.

        Returns dict: file_name, n_mitotic, n_unanimous_pool, agreement_flagged, n_gate_valid (unanimous, border-filtered, Otsu-tightenable), mpp, checks_passed.
    """
    im = images[images['image_id'] == image_id].iloc[0]
    rgb = ds.load_roi(path)
    checks = {
        'loads_hxwx3_uint8': rgb.ndim == 3 and rgb.shape[2] == 3 and rgb.dtype == np.uint8,
        'dims_match_json': rgb.shape[0] == im['height'] and rgb.shape[1] == im['width'],
        'domain_matches': im['tumor_type'] == DOMAIN,
    }
    mpp = ds.roi_mpp(path)
    checks['mpp_sane_0.22_0.26'] = 0.22 <= mpp <= 0.26

    mit = anns[(anns['image_id'] == image_id) & (anns['category_id'] == ds.MITOTIC)]
    pool, flagged = ss.agreement_pool(mit)
    pool = ss.border_filter(pool, 36, rgb.shape)
    hem = ch.to_channel(rgb, 'hematoxylin_od')
    n_gate_valid = 0
    for _, row in pool.iterrows():
        cx, cy = float(row['cx']), float(row['cy'])
        if tm.read_padded_patch(hem, cx, cy, tm.PATCH_SIZE) is None:
            continue
        window = tm.read_padded_patch(hem, cx, cy, tm.BASE_SIZE)
        if window is not None and ss.tighten_box_otsu(window) is not None:
            n_gate_valid += 1
    del rgb, hem

    return dict(file_name=path.name, image_id=image_id, n_mitotic=len(mit), n_unanimous_pool=len(pool),
                agreement_flagged=bool(flagged), n_gate_valid=n_gate_valid, mpp=round(float(mpp), 5),
                checks_passed='ALL' if all(checks.values()) else 'FAILED:' + ','.join(k for k, v in checks.items() if not v))


def materialize(path, timeout=90):
    """
        Force a full local fetch of a possibly cloud-evicted file, and wait until it is fully resident.

        This repo lives under iCloud Desktop & Documents sync (confirmed via `brctl status`) on a
        near-full disk, which evicts a just-written file to an on-demand placeholder almost
        immediately -- `ls`/`ds.load_roi` see the right apparent size but a read can hang or return a
        short/empty buffer while the real bytes are being re-fetched. `du`/`st_blocks` sees the actual
        resident bytes, which is what this polls.

        path (Path): the file to materialize.
        timeout (float): seconds to wait for the fetch to finish before giving up (does not raise).
    """
    apparent = path.stat().st_size
    if apparent == 0:
        return
    subprocess.run(['brctl', 'download', str(path)], check=False, capture_output=True)
    t0 = time.time()
    while time.time() - t0 < timeout:
        if path.stat().st_blocks * 512 >= apparent * 0.99:
            return
        time.sleep(1)


def verify_retrying(path, image_id, images, anns, attempts=4):
    """
        verify(), materializing the file first and retrying on transient file I/O errors.

        path, image_id, images, anns: passed through to verify().
        attempts (int): total tries before giving up.

        Returns dict: verify()'s return value.
    """
    assert attempts >= 1
    for attempt in range(1, attempts + 1):
        materialize(path)
        try:
            return verify(path, image_id, images, anns)
        except Exception as e:
            if attempt == attempts:
                raise
            print(f"    [retry {attempt}/{attempts - 1}] {path.name}: {e!r}; re-materializing", flush=True)
            time.sleep(5 * attempt)
    raise AssertionError('unreachable')


def main():
    """
        Clone the 7 existing ROIs and download the 13 new ones, verify each, and write a report.

        Returns pd.DataFrame: one verification row per ROI.
    """
    images, anns = ds.load_annotations('../databases/MIDOG++.json')
    OUT_DIR.mkdir(exist_ok=True)

    rows = []
    for fn, src in {**EXISTING, **ALREADY_ON_DISK}.items():
        dest = OUT_DIR / fn
        if not dest.exists():
            clone(src, dest)
        image_id = int(images.loc[images['file_name'] == fn, 'image_id'].iloc[0])
        tier = 'existing_7' if fn in EXISTING else 'unanimous_new'
        rows.append(dict(tier=tier, **verify_retrying(dest, image_id, images, anns)))
        print(f"[clone]    {fn}: {rows[-1]}", flush=True)

    for fn, fid in FIDS.items():
        dest = OUT_DIR / fn
        if not dest.exists():
            if shutil.disk_usage(OUT_DIR).free < MIN_FREE_BYTES:
                raise SystemExit(f'less than {MIN_FREE_BYTES / 1024 ** 3:.0f} GiB free; stopping before {fn}')
            download(fid, dest)
        image_id = int(images.loc[images['file_name'] == fn, 'image_id'].iloc[0])
        rows.append(dict(tier='unanimous_new', **verify_retrying(dest, image_id, images, anns)))
        print(f"[download] {fn}: {rows[-1]}", flush=True)

    # backfill from the contested (2-of-3) tier, tried only because the unanimous tier is exhausted
    # (no unflagged candidate above n_gate_valid=2 remains in the domain) -- reported, never silently merged
    contested_rows = []
    for fn, src in ALREADY_ON_DISK_CONTESTED.items():
        dest = OUT_DIR / fn
        if not dest.exists():
            clone(src, dest)
        image_id = int(images.loc[images['file_name'] == fn, 'image_id'].iloc[0])
        contested_rows.append(dict(tier='contested_2of3', **verify_retrying(dest, image_id, images, anns)))
        print(f"[clone, contested tier]    {fn}: {contested_rows[-1]}", flush=True)
    for fn, fid in FIDS_CONTESTED.items():
        dest = OUT_DIR / fn
        if not dest.exists():
            if shutil.disk_usage(OUT_DIR).free < MIN_FREE_BYTES:
                raise SystemExit(f'less than {MIN_FREE_BYTES / 1024 ** 3:.0f} GiB free; stopping before {fn}')
            download(fid, dest)
        image_id = int(images.loc[images['file_name'] == fn, 'image_id'].iloc[0])
        contested_rows.append(dict(tier='contested_2of3', **verify_retrying(dest, image_id, images, anns)))
        print(f"[download, contested tier] {fn}: {contested_rows[-1]}", flush=True)

    report = pd.DataFrame(rows + contested_rows)
    Path('results').mkdir(exist_ok=True)
    report.to_csv('results/pull_verification.csv', index=False)
    print()
    print(report.to_string(index=False))

    bad = report[report['checks_passed'] != 'ALL']
    short = report[report['n_gate_valid'] < 3]
    print(f"\n{len(report)} ROIs verified (by tier: {report['tier'].value_counts().to_dict()}); "
          f"{len(bad)} failed a basic check; {len(short)} have fewer than 3 gate-valid annotations")
    if len(bad):
        print('FAILED:', bad[['file_name', 'checks_passed']].to_string(index=False))
    if len(short):
        print('SHORT:', short[['file_name', 'tier', 'n_gate_valid']].to_string(index=False))
    return report


if __name__ == '__main__':
    main()
