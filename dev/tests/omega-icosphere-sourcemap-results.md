# Source-mapping regression test — icosphere `tess_downsize` (OMEGA / sub-0002)

Verifies that the new FreeSurfer/MNE-style **icosphere** cortex downsampling (see
`dev/ico-downsize.md`) does not break any source-mapping method. Run live against MATLAB R2023b
on `feature/eigenmode-analysis`, dev Brainstorm install.

- **Test script:** `dev/tests/test_omega_icosphere_sourcemap.m`
- **Run:** `results = test_omega_icosphere_sourcemap();`  → `results.pass = 1`
- **Protocol:** `TutorialOmega_IcoTest` (deleted + recreated on every run; safe to leave or delete).

## What it does

Reproduces the cortex-dependent half of `tutorial_omega` (anatomy → forward model → inverse) on a
single OMEGA subject, twice, changing **only the cortex downsampling method**:

| Arm | Subject | Downsampling | Cortex |
|-----|---------|--------------|--------|
| Test | `sub-0002-test` | **icosphere** (new), ico5 | 20,484 vert |
| Baseline | `sub-0002-test-rp` | reducepatch (legacy), 15000 | 15,002 vert |

Both arms use the **same** FreeSurfer recon (`derivatives/freesurfer/sub-0002/ses-mri/anat`), the
**same** resting MEG run (`sub-0002 task-rest`), and the **same** empty-room noise
(`sub-emptyroom task-noise`). Pipeline per arm: link raw → CTF continuous → notch (60/120/180/240/300)
→ high-pass 0.3 Hz → **noise covariance** (empty room) + **data covariance** (resting, for LCMV) →
**overlapping-spheres** head model → four inverses: **MNE, dSPM, sLORETA, LCMV** (all fixed
orientation, MEG, kernel-only). Each inverse kernel is applied to a 5-second block of the real
resting recording to confirm a finite, non-degenerate source map.

## Scope (what this test is / is not)

- **IS** a numerical-integrity test of everything that consumes the cortex geometry produced by
  `tess_downsize`: vertices/faces/normals → overlapping-spheres gain → inverse kernels → source maps.
- **NOT** a localization-accuracy test. The single-subject manual import uses default (MNI-derived)
  fiducials — no digitized head coils in the OMEGA headshape file — so MEG/MRI co-registration is
  approximate. It is **identical for both arms**, so it does not affect the comparison. Steps that do
  not depend on the cortex (head-point refine, SSP, resting power maps) are intentionally omitted.

## Results (OVERALL: PASS)

| Check (per arm) | icosphere `sub-0002-test` | reducepatch `sub-0002-test-rp` |
|---|---|---|
| Cortex vertices | 20,484 | 15,002 |
| Connected components | 2 | 2 |
| Vertex normals finite | yes | yes |
| **2-manifold** | **yes** | **no** (expected — the defect icosphere fixes) |
| Head model gain (modeled rows finite) | 300×61,452, finite | 300×45,006, finite |
| MNE kernel | 20,484×270, finite | 15,002×270, finite |
| dSPM kernel | 20,484×270, finite | 15,002×270, finite |
| sLORETA kernel | 20,484×270, finite | 15,002×270, finite |
| LCMV kernel | 20,484×270, finite | 15,002×270, finite |

End-to-end reconstruction (mean per-vertex temporal std over a 5 s block) is **finite and
non-degenerate** for every method, and the same order of magnitude across the two meshes:

| Method | icosphere | reducepatch |
|---|---|---|
| MNE     | 6.87e-11 | 8.97e-11 |
| dSPM    | 1.51e+00 | 1.49e+00 |
| sLORETA | 1.41e-09 | 1.55e-09 |
| LCMV    | 8.54e-01 | 8.44e-01 |

## Interpretation

- **No source-mapping method is broken by the new `tess_downsize`.** On the icosphere ico5 cortex,
  the overlapping-spheres head model and all four inverse operators compute cleanly (finite,
  correctly sized: rows = 20,484 sources, cols = 270 MEG channels), and produce finite,
  non-degenerate source maps on real data.
- **Equivalent to the legacy path.** Reconstruction magnitudes match the reducepatch baseline
  per method — the icosphere cortex behaves like the existing cortex for source estimation.
- **Bonus quality:** the icosphere cortex is 2-manifold; the reducepatch baseline is not. This is the
  expected, motivating difference, not a regression — so the baseline's `manifold=no` is reported but
  does not count against it.

## Notes on test reliability (checks corrected during bring-up)

Three first-run "failures" were traced to the test harness, not the surfaces (each confirmed by
evidence before fixing):

1. **Head model `finite=no`** — `all(isfinite(Gain(:)))` included 4 non-modeled rows (ECG/EOG/SysClock)
   that overlapping-spheres leaves as `NaN` by design. Fixed: check finiteness over modeled rows only
   (rows not entirely `NaN` = the 270 MEG + 26 MEG REF). Both arms were affected → never a surface issue.
2. **Reconstruction `skipped(iTimesBl)`** — `in_bst(...)`'s 5th argument (`RemoveBaseline`) is a string
   (`'all'`/`'time'`/`'no'`), not numeric; passing `0` left `iTimesBl` unassigned in `in_fread.m`.
   Fixed: pass `'no'`.
3. **reducepatch `manifold=no` flagged as fail** — reducepatch is inherently non-manifold; the manifold
   criterion now gates only the icosphere arm.
