# PET PVC — known limitations / weaknesses (to address)

Branch `experimental/pet`. Recorded 2026-06-26 from the 3-way PVC comparison
(`dev/benchmarks/compare_pvc_3way.m`) of our PETPVE12 Müller-Gärtner (MG) against
PETSurfer MGX (voxelwise MG) and PETSurfer GTM (Rousset regional), on
`sub-MTL0002` / `18FNAV4694`, regional SUVR over FreeSurfer ROIs, cerebellum-normalized.

## Quantitative status — COHORT (10 subjects x 2 tracers, batch_pvc_3way)
POOLED cortical ROIs (N=1360, each ROI cerebellum-normalized):
| comparison | type | pooled r | pooled CCC |
|---|---|---|---|
| OURS vs MGX | cross-pipeline, **same method (MG)** | **0.90** | **0.84** |
| MGX vs GTM | same pipeline (PETSurfer), cross-method | 0.97 | 0.97 |
| OURS vs GTM | cross-pipeline, cross-method | 0.85 | 0.80 |

Per-subject cortex r (mean +/- SD): OURS-MGX amyloid **0.76 +/- 0.09** (robust),
tau **0.60 +/- 0.41** (variable); MGX-GTM ceiling 0.88-0.90 +/- 0.05. So at the cohort
level our MG tracks PETSurfer's MG well (pooled r=0.90, CCC=0.84) — the single-subject
"compression" (sub-MTL0002 CCC 0.26) was NOT representative.

**New finding — tau variability:** flortaucipir is much more variable than amyloid, with
2/10 clear outliers (sub-MTL0020 r=-0.45, sub-MTL0039 r=0.36 vs PETSurfer MG). Suspects:
flortaucipir off-target binding (choroid plexus/basal ganglia/meninges) interacting with
SPM-Segment tissue maps, and/or coregistration on those specific tau scans. → *investigate
the 2 outliers individually; off-target compartments argue for GTM (Baker 2017).*

## Weaknesses (ranked)
1. **GM-only correction (Müller-Gärtner).** MG corrects only the GM compartment and
   assigns WM/CSF constant activity. WM / CSF / subcortical-GM ROIs therefore diverge
   from GTM **by design** (they scatter widely in the comparison). No valid PVC for
   subcortical structures. → *Add a regional GTM / symmetric-GTM (Rousset/Greve) and/or
   RBV option so all tissue classes are corrected; GTM is effectively mandatory for tau
   (flortaucipir off-target in choroid plexus / basal ganglia / meninges).*
2. **Dynamic-range compression (low CCC despite good r).** Our cortical SUVR clusters in
   a narrow band (~1.3–1.8) vs PETSurfer MG's wider spread (~0.7–2.5): r=0.73 but
   CCC=0.26 — a slope/scale mismatch, i.e. we under-recover spill-out variation. The
   synthetic benchmark shows MG at the correct FWHM recovers RC≈1.0, so the compression
   is **not** the core MG math — suspects: SPM-Segment tissue maps vs FreeSurfer gtmseg,
   the `gmThresh` masking, static (mean-of-frames) vs PETSurfer per-frame, WM/CSF
   constant-estimation (`wmcsfMethod`). → *Investigate segmentation + MG parameters; add a
   slope/Bland-Altman calibration check to the comparison.*
3. **Segmentation source mismatch.** We use SPM12 Segment (c1/c2/c3); PETSurfer uses
   FreeSurfer `gtmseg`. Different GM/WM boundaries cap the cross-pipeline same-method
   agreement at 0.73 (vs 0.89 internal). → *Option to derive tissue maps from the
   subject's FreeSurfer segmentation already in Brainstorm (ASEG/ribbon) for consistency.*
4. **No regional / surface-native PVC.** Only voxelwise MG in volume; no GTM regional
   table, no surface SGTM. → *Roadmap: GTM regional + the locked volume→surface step.*
5. **Interim static image.** PVC currently runs on a plain `mri_aggregate` mean; the
   proper metadata-aware, duration-weighted, window-validated static-image step is pending.

## Already fixed this session (not weaknesses)
- PSF FWHM now auto-derived from scanner metadata (`pet_scanner_fwhm`; HRRT→2.5mm),
  generic 6mm only as fallback (was hard-coded 6mm → 15–19% over-correction on HRRT).
- Geometry: resliced volumes kept a stale NIfTI vox2ras → PVC output misaligned. Fixed at
  the source in `mri_reslice` (inherit reference geometry) + `pet_pvc` guard. Cortex
  agreement vs PETSurfer GTM went 0.03 → 0.61 after the fix.
