# Eigenmode source-mapping validation

- Protocol: `omega_tutorial_test`
- Date: 2026-06-02 12:13
- Method under test: eigenmode (GBF) dSPM, prior=log

## Level 1 - Resolution metrics (REQUIRED)

- Study 7 (`@rawsub-0002_ses-01_task-rest_run-01_meg_notch_high`), base HM `sub-0002/@rawsub-0002_ses-01_task-rest_run-01_meg_notch_high/headmodel_surf_os_meg.mat`.
- 270 good MEG channels; 600 eigenmodes available.

| Method | median LocError (mm) | median SpatialDisp (mm) | depth-bias slope (amp/mm) |
|--------|----------------------|--------------------------|----------------------------|
| eigenmode-dSPM | 21.58 | 45.59 | +4.4511e-03 |
| wMNE | 13.24 | 34.93 | +4.5475e-03 |
| dSPM | 10.02 | 35.57 | +5.3205e-03 |
| sLORETA | 0.00 | 37.84 | +6.1167e-03 |

_Depth proxy = distance of each source from the cortex centroid (mm); slope = LSQ fit of peak-PSF amplitude (normalized) vs that distance. A larger-magnitude slope indicates stronger dependence of resolution amplitude on this depth proxy (i.e. more depth bias); slopes near 0 indicate depth-uniform resolution. (Cortex-centroid distance is a crude proxy; treat magnitudes comparatively across methods, not absolutely.)_

## Level 2 - Ground-truth simulation (BEST EFFORT)

- Self-contained simulation: 12 focal seeds spread across the cortex; localization error averaged (median) over seeds.

| SNR (amp) | eigenmode-dSPM median LocError (mm) | standard dSPM median LocError (mm) |
|-----------|--------------------------------------|--------------------------------------|
| 3 | 27.91 | 11.67 |
| 10 | 20.51 | 9.96 |

_Ground truth = leadfield forward of a single cortical vertex + colored sensor noise; error = distance from reconstructed peak to seed, median over seeds. (The GUI simulate processes were not driven headlessly; this self-contained forward simulation provides the same ground-truth comparison.)_

## Level 3 - OMEGA GBF-vs-dSPM (REQUIRED) + phantom (BEST EFFORT)

### OMEGA: spatial correlation of |source| maps (eigenmode-dSPM vs standard dSPM)

| Subject study | nVert | spatial corr (peak time) | spatial corr (time-avg) |
|---------------|-------|--------------------------|--------------------------|
| @rawsub-0002_ses-01_task-rest_run-01_meg_notch_high | 20484 | 0.036 | 0.104 |

### Phantom localization

SKIPPED: no phantom protocol loaded (no protocol with name containing "phantom").

