# Boundary-Based Registration (BBR) in Brainstorm — design

**Date:** 2026-06-26  **Branch:** `experimental/pet`  **Status:** Approved

## Goal
Native replication of FreeSurfer `bbregister` (Greve & Fischl 2009, *NeuroImage* 48:63-72)
to robustly coregister low-contrast PET (esp. flortaucipir/tau) to the subject's anatomy.
Roots the registration on the GM/WM surface instead of weak global PET intensity similarity
— fixing the tau PET→T1 coregistration failures (sub-MTL0020/0039) that SPM's MI coreg
produced. Also sets up future differential/spectral-geometry work (manifold/operators/eigen):
the boundary-normal sampling rides on the same surface machinery.

## Algorithm (faithful to Greve BBR)
A 6-DOF RIGID transform `T` (PET→anat) maximizing GM/WM contrast at the white surface:
1. White-surface vertices `v` (SCS) + outward vertex normals `n` (`tess_normals`).
2. For candidate `T`: sample PET just OUTSIDE (`v + d·n`, GM side) and INSIDE (`v − d·n`, WM
   side) — `cs_convert` SCS→world→voxel (composing `T`) + `interp3` of the PET cube.
3. Per-vertex percent contrast `Q = 100·(I_gm − I_wm)/(½(I_gm+I_wm))`.
4. Cost `= mean(1 + tanh(M·(Q − Q0)))`, minimized over `T` (sharper boundary → lower cost).
   Sampling step `d ≈ 1-2 mm`; defaults `M`, `Q0` per Greve.

## Contrast-sign handling (PET-specific)
PET GM/WM contrast direction varies by tracer (amyloid GM>WM; tau GM≈WM). Estimate the sign
of mean(Q) at the initial pose and orient the cost accordingly (so "correct alignment" is
always the cost minimum), rather than hard-coding the T1/BOLD convention.

## Components
- **`toolbox/anatomy/mri_bbregister.m`** (new): `[Transf, cost] = mri_bbregister(sMriPet,
  WhiteSurf, initTransf, Opts)`. Returns the refined PET→anat 4×4 rigid transform (world
  space). Pure geometry + optimization; reuses `tess_normals`, `cs_convert`, `interp3`,
  `fminsearch` (6 params: 3 rot + 3 trans; all base MATLAB). Opts: `.SampleDist` (1.5mm),
  `.SubSample` (use ≤ N vertices), `.MultiStart` (n).
- **`toolbox/anatomy/mri_coregister.m`**: add `case 'bbr'` — global SPM estimate for the
  init, then `mri_bbregister` refine, applied via the existing `vox2ras` reslice path
  (output identical in form to the `'spm'` method).
- **`dev/preventad_pet_import.m`**: switch the coreg method to `'bbr'` (SPM init → BBR refine).

## Initialization (the main risk — local optimizer)
SPM coreg (global, MI) provides the init; BBR refines within its capture range. The observed
tau drift was a moderate shift (not a gross flip), so SPM's estimate is within range.
**Multi-start fallback:** if the BBR cost stays high, retry from small perturbations of the
init and keep the best — guards against a bad SPM init on the worst scans.

## Validation
Re-register the two tau outliers with BBR; recompute observed-vs-PETSurfer regional pattern
(no PVC). Expect 0020 to jump from −0.13 toward the ~0.9 its amyloid achieves, and 0039 from
0.48 upward. Also confirm the good subjects (amyloid 0002) are unchanged/improved. Cheap,
decisive — reuses the existing observed-vs-PETSurfer comparison.

## Files
- `toolbox/anatomy/mri_bbregister.m` — new core.
- `toolbox/anatomy/mri_coregister.m` — `'bbr'` method.
- `dev/preventad_pet_import.m` — use `'bbr'`.
- validation: reuse the inline observed-vs-PETSurfer check (no new file).

## Out of scope (future)
Manifold/operator/eigen boundary features (keep `mri_bbregister` factored so the sampling
step can be swapped later); cross-tracer same-session shared-init; surface-based PET-MR QC UI.
