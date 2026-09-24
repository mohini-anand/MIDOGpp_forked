"""
    Post-edit per-ROI pass: V1-V4, V8's per-frame checks and V10 via checks.check_roi, one ROI at a time, each ROI's results
    written to out/<tier>/rois/<file>.pkl before the next starts (an ROI whose pickle exists is skipped). Stops on the first
    ROI with a divergence. Usage: run_tier.py --tier smoke|full [--rois A,B] [--shard I/N]
"""

from __future__ import annotations

import argparse
import sys
import time

import common as C
from common import rr
import checks


def main():
    """
        Check every requested ROI; stop on the first divergence.

        Returns None.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--tier", required=True, choices=("smoke", "full"))
    parser.add_argument("--rois", default="", help="comma-separated file names to restrict to")
    parser.add_argument("--shard", default="0/1", help="I/N: ROIs whose build_rois position is I mod N")
    args = parser.parse_args()
    shard_i, shard_n = (int(v) for v in args.shard.split("/"))

    rois, annotations = rr.build_rois()
    refs = dict(r1=C.read_ref(C.R1_PATH), r4_per_run=C.read_ref(C.R4_PER_RUN_PATH), r4_top=C.read_ref(C.R4_TOP30_PATH), r4_clicks=C.read_ref(C.R4_CLICKS_PATH))
    wanted = [s.strip() for s in args.rois.split(",") if s.strip()]
    out_dir = C.OUT / args.tier / "rois"
    t_run = time.time()
    for i, roi in rois.iterrows():
        if i % shard_n != shard_i or (wanted and roi["file_name"] not in wanted):
            continue
        path = out_dir / f"{roi['file_name']}.pkl"
        if path.exists():
            print(f"[{i + 1}/49 {roi['file_name']}] exists, skipped", flush=True)
            continue
        s0 = C.load_pickle(C.STEP0 / f"{roi['file_name']}.pkl")
        res = checks.check_roi(roi, annotations, refs, s0)
        C.save_pickle(res, path)
        n_div = sum(e["n_diverged"] for e in res["tally"].values())
        print(f"[{i + 1}/49 {roi['file_name']}] {roi['subset']}: {sum(e['n_values'] for e in res['tally'].values())} values, {n_div} diverged [{res['t_check_s']:.0f}s; {(time.time() - t_run) / 60:.1f} min]", flush=True)
        if n_div:
            for name, e in res["tally"].items():
                for line in e["detail"]:
                    print(f"  DIVERGENCE {name}: {line}", flush=True)
            sys.exit(f"stopped: {n_div} divergences on {roi['file_name']}")


if __name__ == "__main__":
    main()
