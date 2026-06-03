# Eigenmode Source-Mapping Accuracy Benchmark - Report

Anatomies: auditory, neuromag. Methods: wmne, dspm, sloreta, eig_mne_log, eig_dspm_log. Regimes: focal, patch, distributed. SNR (dB): [2 4 6 10 20]. K (total): [600 1200 2000].

## K-sweep (focal, eig\_mne\_log)

| K | median LocError (mm) |
|---|---|
| 600 | 18.53 |
| 1200 | 19.81 |
| 2000 | 18.15 |

Plateau-K (focal): **600 total modes** (improvement <= 1.0 mm beyond this).

## Competitiveness (paired Wilcoxon, eig vs standard)

| regime | eig | vs | median diff (mm) | p |
|---|---|---|---|---|
| focal | eig_mne_log | wmne | +0.31 | 0.08938 |
| focal | eig_mne_log | dspm | +7.16 | 7.233e-09 |
| focal | eig_mne_log | sloreta | +16.80 | 1.203e-24 |
| focal | eig_dspm_log | wmne | +2.96 | 0.0007728 |
| focal | eig_dspm_log | dspm | +8.25 | 3.196e-11 |
| focal | eig_dspm_log | sloreta | +17.73 | 2.689e-25 |
| patch | eig_mne_log | wmne | +0.00 | 0.04284 |
| patch | eig_mne_log | dspm | +4.03 | 9.722e-06 |
| patch | eig_mne_log | sloreta | +10.42 | 6.521e-23 |
| patch | eig_dspm_log | wmne | +4.51 | 6.393e-11 |
| patch | eig_dspm_log | dspm | +6.91 | 6.096e-12 |
| patch | eig_dspm_log | sloreta | +16.12 | 3.601e-24 |
| distributed | eig_mne_log | wmne | +0.83 | 0.004806 |
| distributed | eig_mne_log | dspm | +6.80 | 1.07e-07 |
| distributed | eig_mne_log | sloreta | +13.11 | 5.4e-19 |
| distributed | eig_dspm_log | wmne | +4.20 | 3.579e-07 |
| distributed | eig_dspm_log | dspm | +11.68 | 1.267e-11 |
| distributed | eig_dspm_log | sloreta | +17.26 | 5.168e-22 |

_Eigenmode is "competitive" where median diff is small and p is not significant, or where the diff favours eig._

## Figures

![f1_distribution](figures/f1_distribution.png)

![f2_snr_sweep](figures/f2_snr_sweep.png)

![f4_per_regime](figures/f4_per_regime.png)

![f5_ksweep](figures/f5_ksweep.png)

![f3_cortex_reconstruction](figures/f3_cortex_reconstruction.png)

