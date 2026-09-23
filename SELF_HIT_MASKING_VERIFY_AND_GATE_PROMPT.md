# Task: verify the image-blanking self-hit design, then build its gates

## Goal

`SELF_HIT_MASKING_PLAN.md` (repo root) describes a design decision that is not yet
execution-ready: replace the current post-NMS self-hit filter by blanking the seed's refined
template footprint out of a copy of the search channel, before correlation runs. Its §2 and
§3 are a real design and real pre-registered numbers, from real simulation this session —
but nothing has been independently re-verified, and none of the gate infrastructure a
patch-applying agent would need (patch, bit-identical reference captures, a `check` script,
negative tests, augmented-bank and multi-seed testing) exists yet. §4 of that document lists
exactly what's missing.

**Your job has two phases, in order, and stops hard between them:**

1. **Verify.** Independently confirm §2's mechanistic claims and §3's numbers are actually
   true, and close the specific evidence gaps §3c names (not vague "more testing" — those
   three named gaps specifically). If anything doesn't hold up, **stop, do not proceed to
   phase 2, and report exactly what diverged.** You do not have authority to redesign the
   mechanism if verification fails — that's a decision for the user, not you.
2. **Build the gates.** Only once phase 1 fully passes: build the real code patch (against a
   scratch copy, never applied to the real `midog_utils/`), fresh bit-identical reference
   captures, a `check`/`reference`-equivalent script with negative tests, and rewrite
   `SELF_HIT_MASKING_PLAN.md` into the same execution-ready shape the superseded design had —
   so that whoever applies this patch to production next can run a short, mechanical
   Steps-0-through-7 procedure instead of re-deriving any of this. Also reconcile
   `DECISIONS.md`.

**You do not apply the patch to the real `midog_utils/` at any point.** That's a separate,
later step for whoever picks up the execution-ready plan you leave behind. Your deliverable
is a *ready* repository, not a *changed* production pipeline.

**Stay narrowly on this one feature.** If you notice something else that looks wrong or
stale in the repo while you're in here, do not fix it, refactor it, or even flag it unless it
directly blocks this task. Report it at the very end, one line, and move on.

## Trust nothing you haven't personally re-derived

This applies to every document you're about to read, not just the plan.

- **The plan's numbers, the archived design's numbers, and any `DECISIONS.md` entry's
  numbers are claims, not facts, until you've reproduced them yourself against real code and
  real data.** A number being written down in this repo, cited with a specific figure, or
  matched by a stored CSV is not evidence it's correct — it's evidence someone (possibly
  you, in an earlier attempt; possibly a previous session) computed it once. Recompute the
  ones you're relying on.
- **The user already reviewed this design in conversation and chose to proceed with it.**
  That is not evidence it's correct either. It reflects a judgment call made with the
  evidence available *at the time*, some of which (§3c's three gaps) was explicitly known to
  be incomplete. Your job is to find out whether that judgment still holds once the gaps are
  closed — not to produce a report that confirms a decision that's already been made. If
  phase 1 turns up a real problem, say so plainly; a report that finds nothing wrong is only
  worth as much as the effort you put into trying to find something wrong.
- **Don't just re-run the existing simulation scripts and treat a matching result as
  independent confirmation.** If `blank_seed_variant_49roi.py` has a subtle bug — in the
  fill-value computation, the match-radius/rank convention, the bucket definitions it reads
  from `evaluate.py`, how it pins seeds via `ann_id` — re-running that same script will
  reproduce the same wrong answer and you'll have "confirmed" nothing. Read the scripts in
  `cleanup_harness/selfmask/image_blank_design/` line by line, adversarially, before you
  trust their output: check that `fill_value` really is computed from the whole ROI on the
  *unblanked* channel before any mutation; check the rank indexing (0-indexed vs 1-indexed)
  matches what the plan's prose claims when it says "rank 10"; check `bucket_detections`'
  bucket definitions actually mean what §3's tables say they mean; check that pinning by
  `ann_id` really does route through the same gates a normal draw would. For the core
  pooled-TP@K numbers, write independent code that computes the same quantity a different
  way (e.g. call the real pipeline functions directly rather than re-deriving the search by
  hand a second time in the same style as the original script) so a shared bug can't hide
  behind agreement between two copies of the same approach.
- **Where a claim rests on a repo document *other* than the plan** (a `DECISIONS.md` entry's
  stated rationale, a docstring's description of what a function does, a comment asserting a
  channel or constant), and that claim is load-bearing for something you're verifying, check
  it against the actual running code rather than citing the document as if citing it settles
  the question. You don't need to re-litigate settled, unrelated decisions (D1-D9) that this
  task doesn't depend on — but anything this specific verification leans on, verify directly.

## Read first, fully, before writing anything

- `SELF_HIT_MASKING_PLAN.md` — the current design and its evidence. This is what you're
  verifying and then building gates for.
- `cleanup_harness/selfmask/SELF_HIT_MASKING_PLAN.valid_mask_disc_design.md` — the
  *superseded* design (masking `valid`, not blanking the image). Read this as a template for
  the **shape** your gate infrastructure should take: its §4 (patch format), §6 (execution
  procedure), §7 (troubleshooting), §8 (audit trail) are the right level of rigor to match.
  Do not reuse its content — the mechanism is different — and do not modify this file or
  anything under `cleanup_harness/selfmask/` that isn't yours to create (see "Don't touch"
  below).
- `cleanup_harness/selfmask/image_blank_design/` — the evidence behind the current plan's
  §3: both simulation scripts, their logs, the 49-ROI per-ROI CSV, and the notebook-seed
  snapshot they pinned from. This is real code you can re-run, not just numbers to trust.
- `DECISIONS.md` **D10** — records the *superseded* design as the decision, status "not yet
  applied." It is currently stale (describes a design this repo isn't pursuing) and
  reconciling it is part of your phase-2 job.
- `midog_utils/find_and_suppress.py`, `template_match.py`, `production.py`,
  `seed_selection.py`, `FIND_AND_SUPPRESS_REFERENCE_DIFFS.md` — today's real, committed code.
  Confirm what you read matches what the plan describes; don't assume the plan's quotes are
  still accurate.
- `cleanup_harness/harness.py` and `cleanup_harness/selfmask/selfmask_check.py` — reusable
  as-is for a lot of your phase-2 machinery (capture/compare/sanity/extref already work
  against whatever code is in `midog_utils/`; `selfmask_check.py`'s `recent-capture`/
  `recent-compare`/`reference`/`nb-text` subcommands are design-agnostic — only its `check`
  subcommand's attribution logic is specific to the superseded design's info keys and would
  need a design-specific equivalent, not a copy).

## Phase 1: verify

Confirm, independently — re-run code, don't just re-read CSVs — that:

1. **The mechanism claims in §2 hold.** Specifically: `TM_CCOEFF` against a flat (blanked)
   window computes near zero, not an extreme value; the corrupted correlation response
   tapers with window overlap rather than stepping sharply at the blanked square's edge (no
   artefact peak forms there); NMS behaves identically before and after — what changes is
   only what's available for it to suppress against.
2. **§3a and §3b's numbers reproduce.** Re-run both simulations (or write your own from
   scratch and cross-check against the archived scripts — either is fine, but at least one of
   your checks should not just re-execute the exact same script unmodified). Confirm: the
   3-ROI and 49-ROI pooled TP@K tables, the seed-geometry sanity check, the control-arm
   fidelity check, the 3 specific (ROI, K) cells that moved and their stated cause, the 2
   near-seed leak instances and their stated mechanism (including the `n_peaks`/`n_after_nms`
   evidence for 289.tiff's NMS-domination-collapse claim specifically — don't take "2 extra
   peaks survived" on faith, verify it).
3. **Close §3c's three named gaps** (this is required, not optional, and doubles as gate
   evidence for phase 2):
   - **`tm_score`**, all 49 ROIs, same pinned seeds: confirm zero TP@K change, matching the
     reasoning in §2 (a freed candidate always sorts last under `tm_score`).
   - **The 8-template augmented bank** (`scales=(0.8,1.2)`, `n_angles=2`, both flips — same
     bank the superseded design and the original self-hit-radius failures used), on however
     many ROIs the harness's existing augmented-bank machinery covers: does the near-seed
     leak get worse under augmentation (more instances, larger distances), the way the
     *original* self-hit-radius filter's real failures did? This is the single highest-value
     unclosed gap — augmentation is where the filter this design replaces actually failed
     historically.
   - **Seeds 1-4**, not just seed 0 (or its `ann_id`-pinned equivalent), on at least the 14
     `images/extra_valid` ROIs: does the near-seed leak rate or the TP@K movement rate change
     at other seeds?
   Also resolve the loose end in §3c: explain the 3 ROIs where `n_detections` didn't move at
   all (check whether `max_peaks_binding` was false for those runs, which would make the
   blanking a no-op on pool size by construction).
4. **Pass/fail:** if every one of the above reproduces within the plan's stated numbers (or
   you can explain a discrepancy as measurement noise, e.g. environment/package-version
   drift — check `cv2`/`numpy`/`pandas`/`skimage` versions against what the archived design's
   provenance recorded), proceed to phase 2. If `tm_score` shows any TP@K movement, if the
   augmented bank produces a leak that reaches a top-K list or costs a TP, if a correlation
   artefact peak actually appears anywhere, or if any pooled/per-cell number in §3 doesn't
   reproduce — **stop. Write up exactly what you found and why it contradicts the plan.
   Do not proceed to phase 2, do not edit `DECISIONS.md`, do not build a patch.** Report this
   as your final deliverable instead.

## Phase 2: build the gates (only after phase 1 fully passes)

Work in `cleanup_harness/selfmask/image_blank_design/` for everything new. Do not touch
`cleanup_harness/selfmask/selfmask_code.patch`, its reference captures, or anything else
belonging to the superseded design.

1. **Write the real patch** against `midog_utils/` (5-ish files, likely
   `template_match.py`, `find_and_suppress.py`, `production.py`, `seed_selection.py`'s
   docstring, `FIND_AND_SUPPRESS_REFERENCE_DIFFS.md` — confirm the exact set by writing the
   edit, not by assuming it matches the superseded design's file list). Key constraints from
   §2 of the plan, not optional:
   - The function that blanks the search channel must return a **copy**, never mutate the
     caller's array in place (`find_and_suppress` doesn't own `img_channel`).
   - The template must be cut from the **original, unblanked** channel before any blanking.
   - Blanking happens **before** `fused_response`, not after the threshold (this design runs
     earlier in the pipeline than the superseded one did).
   - `chromatin_od` ranking scores against the **original, unblanked** channel.
   - Add an info key recording what was blanked (e.g. `n_blanked_px`), matching the
     superseded design's `n_seed_masked_px` in spirit.
   - Build and test this patch against a **scratch copy** of the repo (mirror the archived
     design's §8 methodology: copy `midog_utils/` + `production_pipeline/` to a scratch
     directory, symlink `images/`/`databases/`, apply the patch there, run `harness.py`
     — pointed at the copy — from there). The real `midog_utils/` is never touched.
2. **Fresh reference captures**, built from that scratch copy: `harness.py capture`/
   `compare` sanity, plus recent-experiment-style captures covering everything phase 1
   measured (49-ROI chromatin_od and tm_score, the augmented bank, seeds 1-4). These become
   the new `selfmask_ref_pre`/`selfmask_ref_post`-equivalent this design's own gates check
   against — pick label names that don't collide with the superseded design's `runs/`
   entries.
3. **A `check`/`reference`-equivalent script**, with the same job the superseded design's
   `check` did: for every run pair, verify the info-key changes are exactly what's expected,
   every other field is identical, no detection lies inside the blanked square (and, per this
   design's known cost, explicitly allow — and report, don't fail on — a detection within one
   match radius but outside the blanked square, since that's the accepted tradeoff, not a
   bug), every lost/added row is explained. **Write negative tests** proving this script
   actually catches tampering — mirror the superseded design's 9 negative tests
   (`REFERENCE_PROVENANCE.md` in the archived folder lists them) adapted to this design's
   info keys and invariants.
4. **Rewrite `SELF_HIT_MASKING_PLAN.md`** (the live file, in place) once all of the above
   exists and passes: keep §1-3 (updated with whatever phase 1 closed — the plan's §3c gaps
   should shrink or disappear), and replace §4 with real §4 (code edits)/§5 (doc edits,
   mirroring the superseded design's §5 structure and its "locate every edit by quoted
   anchor text" discipline)/§6 (execution and verification, a Steps-0-through-7 procedure a
   future agent can run mechanically)/§7 (troubleshooting)/§8 (audit trail) — matching the
   superseded design's shape closely enough that someone who's read one can navigate the
   other.
5. **`DECISIONS.md`.** Add a new entry (the next letter after D10) recording this design as
   decided, status "not yet applied" — mirroring D10's own structure and level of evidence
   citation. Do not rewrite D10's body (its own rule); instead append a short note to D10
   pointing at the new entry, the same way D10 itself handles amendments. Get the entry
   number right by checking what the highest existing `## D` heading actually is at the time
   you write this, not by assuming D11.

## Don't

- Don't apply anything to the real `midog_utils/` or `production_pipeline/`.
- Don't touch `midog_utils_full/`, `DECISIONS_UNVERIFIED.md`, `D8_TEMPLATE_ANCHOR.md`,
  `INVARIANTS_HISTORY.md`, `PRODUCTION_PIPELINE_CLEANUP*.md`, `Research Logs/`, any existing
  `DECISIONS.md` entry's body, or anything under `cleanup_harness/selfmask/` that belongs to
  the superseded design or its `superseded_5um/` subfolder.
- Don't re-run or edit `bbox_refinement_three_way_chromatin_od_audit.py` (~line 300) or
  `hembbox_precision_at_k_audit_round3.py` (~line 219) — they read production source as text
  and are known to already report stale values; not your problem to fix.
- Don't touch `production_hematoxylin_only/bbox_refinement_three_way_chromatin_od_49roi_3seed.ipynb`
  or its prompt — it will break under this design the same way it breaks under the
  superseded one (reads `prod.SELF_HIT_RADIUS`/`n_self_hits`). Report it, don't fix it.
- Don't second-guess or redesign the mechanism (literal template footprint, whole-ROI-minimum
  fill, blank before correlation) unless phase 1 verification actually falsifies it — that's
  the stop condition above, not an invitation to pick a different design you think is better.
- Don't commit anything unless asked.
- Use `/Users/mohinianand/anaconda3/bin/python3` for anything touching `cv2` — the bare
  `python3` can't import it.
- Before starting phase 2 specifically, check no other notebook/pipeline job is running
  (`pgrep -fl "[n]bconvert|[i]pykernel|[h]arness\.py|[s]elfmask_check\.py"` should print
  nothing) and that `git status --short midog_utils production_pipeline` is clean, same
  preconditions the superseded design's Step 0 required — building against a moving baseline
  produces gates that don't mean anything.

## Report back

- Phase 1: which claims you re-derived and how, whether each matched, and the three
  previously-open gaps' results in full (tm_score, augmented bank, seeds 1-4) — even if
  everything passed, show the actual numbers, not just "confirmed."
- If you stopped after phase 1: exactly what didn't reproduce, your best explanation, and
  explicitly that you did not proceed further.
- If you completed phase 2: what the patch touches, where the new reference captures and
  checker live, the negative-test results, confirmation `midog_utils/` is unmodified
  (`git status --short midog_utils`), and a one-line pointer to the rewritten plan's new §6
  so the next agent knows where to start.
- One line, at most, for anything unrelated you noticed but didn't touch.
