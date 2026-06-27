# PET GTM (Rousset) regional PVC — design

**Date:** 2026-06-26  **Branch:** `experimental/pet`  **Status:** Approved (approach A)

## Goal
Add a Geometric Transfer Matrix (GTM, Rousset 1998) regional partial-volume-correction
option to the Brainstorm PET pipeline, alongside the existing voxelwise Müller-Gärtner
(MG, via PETPVE12). GTM corrects **every** region (cortical + subcortical + WM + CSF +
cerebellum), addressing MG's GM-only weakness, and is the field-standard regional method
(PETSurfer's GTM is the same algorithm). Native implementation — no PETPVE12 job / atlas
export / XLS round-trip, and no NIfTI geometry pitfalls (operates on Brainstorm cubes).

## Algorithm (Rousset GTM)
For a parcellation into R regions with masks `1_Rj`, scanner PSF `h` (auto-FWHM from
metadata), and observed PET:
- GTM matrix `W (R×R)`: `W(i,j) = mean over region i of (h ⊗ 1_Rj)` — the fraction of
  region j's unit activity that, after PSF blur, lands in region i.
- Observed regional means `m(i) = mean of PET over region i`.
- True (corrected) regional activities: `t = W \ m`  (T = GTM⁻¹·m).

## Component: `pet_gtm.m` (new)
`[MriFileGtm, errMsg, regTable] = pet_gtm(PetFile, fwhm, gtmOpts)`
1. Load the static 3D PET (Brainstorm DB).
2. Build a complete parcellation `L` on the PET/T1 grid: **Desikan-Killiany** (cortical
   1000+/2000+) ∪ **ASEG** (subcortical, WM, CSF, cerebellum, brainstem, <1000). Unlabeled
   voxels (air/skull/extracerebral) → one "rest" region so the partition is complete.
3. FWHM: `pet_scanner_fwhm` if not supplied (HRRT → 2.5 mm).
4. Regions = unique labels with ≥ `minVox` (default 50); smaller regions folded into "rest"
   (logged). Build `W` (separable Gaussian conv per region mask), `m`, solve `t = W \ m`;
   if `cond(W)` is large, fall back to `pinv` and warn. Report `cond(W)`.
5. Output a **piecewise-constant corrected volume** (each region's voxels = `t_i`) saved as
   a Brainstorm anatomy node tagged `_gtmpvc`, geometry inherited from the input PET (so it
   stays voxel-aligned — same grid, no NIfTI round-trip). Plus `regTable`
   (regionId, name, nvox, observed, corrected).

## Integration
`pet_pvc` gains `pvcOpts.method` (`'mg'` default | `'gtm'`); for `'gtm'` it delegates to
`pet_gtm` and returns. Downstream SUVR (`mri_rescale`) + surface projection
(`mri_interp_vol2tess`) are unchanged — they consume the corrected volume as before. Tree
menu / `panel_process_pet` can expose the method choice later (not in this piece).

## Validation
Extend / reuse `compare_pvc_3way`: with our GTM volume, **OURS-GTM vs PETSurfer-GTM becomes
a SAME-method comparison** (was cross-method) → expect agreement near the PETSurfer internal
ceiling (~0.9) for shared cortical/subcortical regions. Spot-check `cond(W)` and that
cerebellum-normalized cortical SUVR is sensible (~1–2). Run on `sub-MTL0002` both tracers
first, then a few cohort subjects.

## Known limitations (recorded in dev/2026-06-26-pet-pvc-limitations.md)
- Region set = Desikan+ASEG; **no extracerebral/off-target regions** (PETSurfer `gtmseg`
  adds skull/meninges/extra) → tau off-target still unmodeled. Future: `gtmseg`-style atlas.
- Classic Rousset GTM (not symmetric sGTM); piecewise-constant output (no within-region
  detail — inherent to regional GTM).

## Files
- `toolbox/anatomy/pet_gtm.m` — new.
- `toolbox/anatomy/pet_pvc.m` — `method` dispatch to `pet_gtm`.
- `dev/benchmarks/compare_pvc_3way.m` — allow our-GTM as the "ours" input (optional flag).
