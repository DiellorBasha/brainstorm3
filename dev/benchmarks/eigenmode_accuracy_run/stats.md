# Eigenmode accuracy benchmark - statistics

## Median LocError (mm) by method x regime

| regime | method | median | IQR | n |
|---|---|---|---|---|
| focal | wmne | 15.25 | 16.35 | 150 |
| focal | dspm | 9.34 | 9.34 | 150 |
| focal | sloreta | 0.00 | 0.00 | 150 |
| focal | eig_mne_log | 18.15 | 18.60 | 150 |
| focal | eig_dspm_log | 18.41 | 19.71 | 150 |
| patch | wmne | 13.79 | 11.83 | 150 |
| patch | dspm | 12.20 | 8.74 | 150 |
| patch | sloreta | 4.65 | 4.78 | 150 |
| patch | eig_mne_log | 15.05 | 19.81 | 150 |
| patch | eig_dspm_log | 20.77 | 24.41 | 150 |
| distributed | wmne | 20.05 | 14.26 | 150 |
| distributed | dspm | 16.38 | 14.50 | 150 |
| distributed | sloreta | 10.94 | 8.62 | 150 |
| distributed | eig_mne_log | 23.10 | 26.29 | 150 |
| distributed | eig_dspm_log | 26.63 | 28.08 | 150 |

## Paired Wilcoxon (eig vs standard), LocError

| regime | eig | vs | median diff (mm) | n | p |
|---|---|---|---|---|---|
| focal | eig_mne_log | wmne | +0.31 | 150 | 0.08938 |
| focal | eig_mne_log | dspm | +7.16 | 150 | 7.233e-09 |
| focal | eig_mne_log | sloreta | +16.80 | 150 | 1.203e-24 |
| focal | eig_dspm_log | wmne | +2.96 | 150 | 0.0007728 |
| focal | eig_dspm_log | dspm | +8.25 | 150 | 3.196e-11 |
| focal | eig_dspm_log | sloreta | +17.73 | 150 | 2.689e-25 |
| patch | eig_mne_log | wmne | +0.00 | 150 | 0.04284 |
| patch | eig_mne_log | dspm | +4.03 | 150 | 9.722e-06 |
| patch | eig_mne_log | sloreta | +10.42 | 150 | 6.521e-23 |
| patch | eig_dspm_log | wmne | +4.51 | 150 | 6.393e-11 |
| patch | eig_dspm_log | dspm | +6.91 | 150 | 6.096e-12 |
| patch | eig_dspm_log | sloreta | +16.12 | 150 | 3.601e-24 |
| distributed | eig_mne_log | wmne | +0.83 | 150 | 0.004806 |
| distributed | eig_mne_log | dspm | +6.80 | 150 | 1.07e-07 |
| distributed | eig_mne_log | sloreta | +13.11 | 150 | 5.4e-19 |
| distributed | eig_dspm_log | wmne | +4.20 | 150 | 3.579e-07 |
| distributed | eig_dspm_log | dspm | +11.68 | 150 | 1.267e-11 |
| distributed | eig_dspm_log | sloreta | +17.26 | 150 | 5.168e-22 |
