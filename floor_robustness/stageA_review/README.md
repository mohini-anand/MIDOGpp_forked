# Stage A plan review: scripts and outputs

**Date:** 2026-10-08. A reviewer agent wrote these while checking `TM_FLOOR_STRESS_STAGE_A_PLAN.md` (repo root) over
three rounds. They were copied here from the session scratchpad unchanged, so the numbers the plan cites in §3.6 and §6
can be re-run.

**Status:** scratch review code, not the Stage A runner.
- The pilot uses **7 hand-picked extreme ROIs**: clamp-heavy (128), pigment (211), glass-rich (401) and the largest
  (529), plus 201, 343 and 403. **These are not 49-ROI estimates.**
- Nothing here runs a template search.
- Nothing is committed.

## How to run

- **Directory.** Run from inside this folder:
  - `v3_dump.py` and `v3_pl.py` read and write `v3_dump_*.csv` in the working directory;
  - `v2_variants_blocked.py` runs `v2_variants.py` from it.
- **Interpreter.** Use `~/anaconda3/bin/python`.
- **Repo path.** `common.py`, `v3_pl.py` and `v3_owncurve_se.py` hard-code the repo path
  `/Users/mohinianand/Library/Mobile Documents/com~apple~CloudDocs/Desktop/AnnotateDx/MIDOGpp_forked`.
- **Threads.** Pin them to 1 before numpy is imported, e.g.
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ~/anaconda3/bin/python v2_variants.py 128.tiff`. Unpinned, one synthesis
  took about 950 s instead of 1.44 s when 3 processes ran at once.
- **Helpers.** `common.py` has the shared ones: loading, `roi_floors`, `window_floor`, the v1 synthesis and `c_eq`.
  `v2lib.py` re-implements the v2/v3 synthesis and the 16 variants.

## Files by round

### Round 1 (plan v1)

These use the v1 synthesis: full concentrations, with zero-channel pixels synthesised.

| script | output | what it checks |
|---|---|---|
| `check_acd.py` | `check_acd.out` | NULL2 bit-exactness; `c_eq` self-consistency and Rres(b0) = 0; the unperturbed dry run (201, 343) |
| `check_bef.py` | `check_bef_529.out` | the identity round trip, including with a random wrong basis; the orthogonality check; peak RSS (5.08 GB on 529); timing |
| `check_g.py` | `check_g_{201,343,529}.out`, `check_g_{128,401}.tiff.out` | every realistic perturbation plus GAIN, WB and JPG95; `c_eq` and Rres against the 0.75×/0.9× references, for `b_otsu`, TOTSU, W128 and TW128; QMATCH rows. Timings in the 201/343/529 files are contaminated by the thread stall |
| `check_t2.py` | `check_t2.out` | multi-Otsu bin widths; per-click jitter under JPEG q95; SHIFT (211, 401, 128, 435, 509) |
| `check_box.py` | `check_box.out` | gate 3: `b0` and W128 against the recorded values, and boxes reproduced (8 ROIs) |
| `check_chgamma.py` | `check_chgamma.out` | per-channel tone curves; NULL2 Rres exactness |
| `time_synth.py` | — | synthesis against `rgb2hed`, timed pinned and unpinned |

### Round 2 (plan v2)

| script | output | what it checks |
|---|---|---|
| `v2_direction.py` | `v2_direction_201.out` | the direction gate: slopes ≈ 1.01; without the 1/6 factor 0.168; with the k column swapped it fails; an [e, h] label swap still passes |
| `v2_variants.py` | `v2_variants_{128,211,403}.tiff.out`, `v2_h06e15_{201,343,401}.tiff.out` | per-ROI `c_eq` for all v2 variants |
| `v2_blocked_hem.py` | `v2_blocked_hem.out` | production `hem` on row blocks is bit-identical to the whole-image call (529: original, H06, JPG, BLUR) |
| `v2_variants_blocked.py` | `v2_variants_blocked_529.out` | full-worker peak RSS with blocked `hem` (2.46 GB) and runtime (22.4 s) |
| `v2_range.py` | `v2_range.out` | zero-channel pixel counts, `hem` max, and bin width over b′ under synthesis |

### Round 3 (plan v3, the PL rule)

| script | output | what it checks |
|---|---|---|
| `v3_dump.py` | `v3_dump_<ROI>.csv`, `v3_dump.log` | per-click `c_eq` for `b_otsu` and W128 under every variant, plus unperturbed and SHIFT (7 ROIs) |
| `v3_pl.py` | `v3_pl.out` | PL and OVER on the pilot. Also: unperturbed W128 PL on all 49 ROIs (0.056); PL against the searched W128 arm (predicted −0.056, observed +0.24); own-curve PL |
| `v3_owncurve_se.py` | `v3_owncurve_se.out` | SE of a per-ROI-curve PL; correlation of the per-ROI response with landing covariates |
| — | `v3_contrib.out` | per-ROI contributions to PL |
