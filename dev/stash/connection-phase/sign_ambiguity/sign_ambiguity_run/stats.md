# Sign-ambiguity continuity test

- Surface: `Subject01/tess_cortex_mid_low.mat`
- Kernel: `link|Subject01/S01_AEF_20131218_01_notch/results_MN_MEG_KERNEL_260605_0111.mat|Subject01/S01_AEF_20131218_01_notch/data_deviant_average_260603_0914.mat`

## Component 1 (data-free smoothness)

| metric | sulcal | crown |
|---|---|---|
| normal var dn (median rad) | 0.4441 | 0.5336 |
| Fiedler var df (median rad) | 0.0690 | 0.0766 |

- singularity energy fraction of df^2: 0.013
- decoupling: dn sulci/crown ratio = 0.83, df sulci/crown ratio = 0.90

## Component 2 (real-data phase continuity)

- peak time: 91 ms, conditioned pairs: 65
- constrained sign-flip rate: 0.923

| frame | median phase disc (rad) | disc rate (>pi/2) |
|---|---|---|
| Fiedler | 1.369 | 0.477 |
| tess_tangents | 1.096 | 0.369 |
| global-xyz | 1.298 | 0.446 |
