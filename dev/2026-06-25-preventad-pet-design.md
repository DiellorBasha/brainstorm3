# PREVENT-AD PET import — register to existing MEG anatomy

**Date:** 2026-06-25 (rev 2026-06-26)   **Status:** IMPLEMENTED + validated (4D base)

## 2026-06-26 REWORK — 4D registered base (current design)

`dev/preventad_pet_import.m` now produces, per tracer, the **4D dynamic PET volume
registered + resliced to the subject's T1** as the BASE PET data — no averaging, no
SUVR, no surface projection. Pipeline per tracer:
`import_mri (4D)` → `mri_realign(..., 'ignore')` (keep all frames) →
`mri_coregister(..., 'spm', 1)` (4D estimate + frame-wise mri_reslice) → label the
result **`PET <tracer>`** (like the original BIDS file) → delete the raw import.
Realign + coregister run on in-memory STRUCTS, so no intermediate DB nodes are
created (only the raw import, which is then removed). Validated on sub-MTL0005:
`PET 18FNAV4694` [256x256x256x6] + `PET 18Fflortaucipir` [256x256x256x4] on the T1
grid; integral conserved across all frames.
⚑ Required a fix to `toolbox/anatomy/mri_coregister.m` (4D path): SPM's own reslice
of the "other" frames silently dropped ~58% of every non-first frame; now
SPM-estimate + mri_reslice (committed b81b4f84). SUVR / masking / surface projection
are NO LONGER part of the import — they become downstream steps on the 4D base.

---

## Original design (2026-06-25, SUVR-based — SUPERSEDED by the 4D rework above)

## Goal

Batch-import the PREVENT-AD PET volumes and register them to the SAME subjects /
FreeSurfer ico5 anatomy already imported by the MEG pipeline, into the existing
`preventad` protocol. Per-subject MEG ↔ amyloid ↔ tau fusion on a shared mesh.

## Dataset

`/Volumes/SpikeData-2/workspace/library/datasets/preventad/pet` — BIDS, 128 subjects,
2 tracers each (`ses-01/pet/`): `trc-18FNAV4694` (amyloid), `trc-18Fflortaucipir`
(tau). Dynamic: 6 frames × 5 min. **111 of 128 match the MEG FreeSurfer cohort**
(only matched subjects are processed; PET registers to the MEG anatomy).

## Why this is feasible / clean

- `import_mri(iSubject, file, [], 0, 0, 'PET ...')` imports a volume INTO an
  existing subject by index → PET attaches to the MEG anatomy. No anatomy re-import.
- Subject matching is trivial (identical BIDS ids). PET = ses-01, MEG = ses-02, same
  subject anatomy.
- Surface-projected SUVR lands on the same ico5 cortex as the MEG/Dirac source maps.
- Fast vs MEG: no headmodel/inverse — realign + coregister + SUVR(×2) ≈ a few min/subj.

## Per-tracer pipeline (from tutorial_pet_processing, Basha 2025)

1. `import_mri(iSubject, petFile, [], 0, 0, ['PET ' tracer])`  (Comment must contain
   'PET' → volType=PET; carries the trc tag)
2. `mri_realign(imp, 'spm_realign', 0, 'mean')`   (6 frames → mean)
3. `mri_coregister(agg, T1, 'spm', 1)`            (to the subject's FreeSurfer T1 + reslice)
4. `pet_process(coreg, 'ASEG', 'Cerebellum', 'Brainmask', 1, 1)`  (SUVR rescale +
   brainmask + project to cortex surface)

NB: PET steps are DIRECT anatomy functions (no `process_` wrappers exist for
realign/coregister/pet_process; only `import_mri` has one) — fully scriptable, just
not `bst_process('CallProcess')` like the MEG batch.

## Scripts (dev/)

- `preventad_pet_import.m`  — one subject: resolve existing subject + T1 + ASEG atlas,
  discover ses-*/pet tracers, run the 4-step chain per tracer. Idempotent (skips a
  tracer whose SUVR anatomy — comment contains tracer tag + `rescaled` — already
  exists). Opts: AtlasName/RefROI/MaskROI/ApplyMask/DoProject/DoSUVR/Aggregation/CoregMethod.
- `preventad_pet_import_batch.m` — discover PET subjects, process only those present
  in the protocol (skip unmatched = no anatomy), per-subject try/catch, timestamped
  log, resumable. Mirrors `preventad_import_batch`.

## Defaults / decisions

- SUVR reference ROI = `Cerebellum` (GUI default; standard for amyloid & tau). Single
  ref for both tracers; can be overridden per need.
- Atlas Comment = `ASEG`; mask = `Brainmask`; project to surface = on.
- Runs in the SAME `preventad` protocol; never deletes/recreates it.

## VERIFY before running (needs MATLAB; do tomorrow, batch idle)

1. ⚑ ASEG `*_volatlas` exists on the imported subjects — the MEG import used
   `process_import_bids` (icosphere); confirm it set isVolumeAtlas so the ASEG atlas
   exists (pet_process needs it). If absent: import aseg.mgz per subject or re-run with
   isVolumeAtlas. Also confirm the atlas Comment is exactly `ASEG` and that
   `Cerebellum`/`Brainmask` are valid ROI labels in it.
2. Static lint (`check_matlab_code`) both scripts.
3. Single-subject validation on `sub-MTL0002` (already has anatomy) before batch.
4. Headless launch: same recipe as MEG — `db_save; brainstorm stop` in the MCP
   session first, then `matlab -batch "addpath(ROOT); brainstorm nogui; addpath(dev);
   preventad_pet_import_batch(...)"`.
