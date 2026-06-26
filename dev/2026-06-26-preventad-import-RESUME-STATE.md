# PREVENT-AD MEG import — RESUME STATE (handoff)

**Snapshot taken:** 2026-06-26 ~11:00.
**Status:** PAUSED at 66/111 subjects. Batch process NOT running. Brainstorm GUI is
running interactively (user working in it). Resume is safe and deterministic.

---

## 1. Where we are

- **66 / 111 FreeSurfer subjects fully imported** (64 via batch + `sub-MTL0002`,
  `sub-MTL0005` done earlier). Protocol has 67 subjects (66 MTL + `Group_analysis`).
- **0 failures.** Paused cleanly: `sub-MTL0350` (#66) was the last full DONE;
  `sub-MTL0353` (#67) had just started and was **deleted** during the pause (so no
  half-imported subject remains).
- **45 subjects remain** (exact list — name-based resume will import exactly these):

```
sub-MTL0353 sub-MTL0363 sub-MTL0375 sub-MTL0376 sub-MTL0380 sub-MTL0387 sub-MTL0391
sub-MTL0394 sub-MTL0395 sub-MTL0404 sub-MTL0414 sub-MTL0425 sub-MTL0432 sub-MTL0442
sub-MTL0453 sub-MTL0457 sub-MTL0458 sub-MTL0461 sub-MTL0462 sub-MTL0473 sub-MTL0488
sub-MTL0492 sub-MTL0504 sub-MTL0509 sub-MTL0520 sub-MTL0521 sub-MTL0534 sub-MTL0537
sub-MTL0550 sub-MTL0552 sub-MTL0572 sub-MTL0574 sub-MTL0589 sub-MTL0598 sub-MTL0599
sub-MTL0602 sub-MTL0606 sub-MTL0608 sub-MTL0612 sub-MTL0632 sub-MTL0648 sub-MTL0649
sub-MTL0653 sub-MTL0669 sub-MTL0703
```
(17 of the 128 PET/whole-cohort subjects have no FreeSurfer and are out of scope; 111
is the ico5-importable set.)

## 2. Key locations

| Item | Path | Size |
|---|---|---|
| Brainstorm dev repo (branch `development`) | `/Users/diellorbasha/workspace/research/code/brainstorm3` | 1.6 GB |
| Import script (single subj) | `…/brainstorm3/dev/preventad_import.m` | — |
| Import batch wrapper | `…/brainstorm3/dev/preventad_import_batch.m` | — |
| Custom process (manifold) | `…/brainstorm3/toolbox/process/functions/process_tess_manifold.m` | — |
| Protocol (imported data) | `/Volumes/SpikeData-2/workspace/library/datasets/brainstorm_db/preventad` (STUDIES=`/data`, SUBJECTS=`/anat`) | 150 GB |
| BIDS source (raw .ds + FreeSurfer) | `/Volumes/SpikeData-2/workspace/library/datasets/preventad/meg` | 300 GB |
| Batch progress log | `…/brainstorm3/dev/preventad_import_batch.log` | — |
| nxr-compute plugin (managed, NOT in repo) | `~/.brainstorm/plugins/nxr-compute` | 42 MB |

## 3. Pipeline (what each subject gets) — see dev/2026-06-25-preventad-import-design.md

Import BIDS (icosphere **ico5**, 20484 v) → remove/refine headpoints → CTF continuous
→ **resample 2400→1200 Hz** → notch(60-300) + highpass 0.3 → PSD snapshot → delete
intermediates → artifact: detect cardiac(ECG)+blink(VEOG)+saccade(HEOG), remove
cardiac∩blink, SSP cardiac/blink/saccade → **bad-segment detect → rename bad_** →
noisecov(task-noise) → overlapping-spheres headmodel → **manifold (gauge=trivial)** →
**dSPM UNCONSTRAINED** (`SourceOrient 'free'`) → **Dirac dSPM** (process_inverse_dirac,
400 modes/hemi, tau 0.5) → power maps (standard dSPM only). Imports into EXISTING
`preventad` protocol (never deletes it). Per subject ≈ 2 GB, ≈ 7-10 min.

## 4. HOW TO RESUME (local)

1. If the interactive Brainstorm is running (it is now), stop it first so the headless
   run is the SOLE DB writer:
   - in the MCP/interactive MATLAB: `db_save; brainstorm stop`
2. Launch headless (snapshots already disabled in the batch → no figure-render hang):
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3/dev
MLB=/Applications/MATLAB_R2023b.app/bin/matlab
ROOT=/Users/diellorbasha/workspace/research/code/brainstorm3
BIDS=/Volumes/SpikeData-2/workspace/library/datasets/preventad/meg
LOG=$ROOT/dev/preventad_import_batch.log
nohup "$MLB" -nosplash -batch "addpath('$ROOT'); brainstorm nogui; addpath('$ROOT/dev'); preventad_import_batch('$BIDS','$LOG')" >> "$ROOT/dev/preventad_batch_console.log" 2>&1 &
```
   It skips the 66 done subjects by name and imports the 45 remaining (~6 h).

## 5. HOW TO RESUME (HPC cluster) — see cluster notes

Copy (rsync): the repo (#2), the protocol (#3 — or skip it and import-remaining-into-a
-fresh-protocol then merge), the nxr plugin (#4), the BIDS source (#5). Then:
- ⚑ **Rebuild `nxr_compute.mexa64`** from current source (bundled Linux MEX is
  2026-06-01, **predates the Dirac/DEC/facets code** the import needs → would fail).
  Match MATLAB **R2023b**.
- Register the copied protocol at its new path (internal DB refs are relative →
  relocatable; only the global registration is absolute).
- Run under `xvfb-run` (headless Linux has no display; Brainstorm starts an OpenGL
  engine at launch). Same `preventad_import_batch(clusterBids, clusterLog)` call.

## 6. nxr-compute build context (for the upcoming adaptation work)

- The import depends on RECENT nxr features: `tess_manifold` (facets + DEC operators,
  Embedded schema v2 incl face.area) and `bst_dirac`/`process_inverse_dirac` (Dirac
  eigenmodes). Dirac landed after 2026-06-10.
- Plugin holds multi-platform MEX in `~/.brainstorm/plugins/nxr-compute/nxr-compute-mex-r2023b/`:
  `nxr_compute.mexmaca64` (macOS arm64, **live 2026-06-23**), `.mexa64` (Linux,
  **stale 2026-06-01**), `.mexw64` (Windows, stale 2026-06-01).
- ⚑ **Managed-plugin stale-binary trap**: after rebuilding nxr, you must copy the new
  MEX into the `~/.brainstorm/plugins/...-r2023b/` folder or Brainstorm loads the old
  one. (Same trap applies per-platform on the cluster.)

## 7. Known issues / gotchas

- **Snapshot hang (fixed):** headless 3D contact-sheet render hung on `sub-MTL0032`
  in run 1. Fix: `preventad_import(..., DoSnapshots)` 3rd arg; batch passes `false`.
  Lost only cosmetic QC thumbnails; data unaffected.
- **`sub-MTL0032` consistency check — PENDING:** its run-1 process was killed mid-final
  step (after data writes), so a few late files (likely the `Group_analysis` power-map
  links) may be on disk but not in `protocol.mat`. After the batch finishes: `bst_get`
  0032 vs disk; if mismatched, `db_reload_database`. Non-destructive; data is intact
  (manifold + 4 kernels confirmed on disk).
- **Resume skips by NAME** → any subject left half-imported is wrongly skipped. Always
  DELETE a partial subject before resuming (we did this for `sub-MTL0353`). A future
  hardening would be to make the batch skip only on a completeness check (e.g. presence
  of the Dirac kernel) and delete+redo incompletes automatically.
- **One DB writer at a time:** never run the headless batch and an interactive
  Brainstorm (or two cluster jobs) against the same protocol simultaneously.

## 8. Post-completion verification (when all 111 done)

1. `ls anat | grep -c sub-MTL` == 111.
2. Spot-check a few subjects: cortex 20484 (ico5), manifold node, unconstrained dSPM +
   Dirac dSPM kernels (both rest runs), SSP cardiac/blink/saccade, bad_* events.
3. Run the `sub-MTL0032` index check (item 7).
4. Then PET import (see dev/2026-06-26-preventad-pet-design.md — drafted, on hold).
