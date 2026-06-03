# Eigenmode accuracy benchmark - statistics

## Median LocError (mm) by method x regime

| regime | method | median | IQR | n |
|---|---|---|---|---|
| focal | wmne | 12.88 | 14.49 | 150 |
| focal | dspm | 10.71 | 8.10 | 150 |
| focal | sloreta | 0.00 | 0.00 | 150 |
| focal | eig_mne_log | 17.13 | 19.13 | 150 |
| focal | eig_dspm_log | 17.40 | 19.70 | 150 |
| patch | wmne | 12.11 | 12.04 | 150 |
| patch | dspm | 13.02 | 8.99 | 150 |
| patch | sloreta | 5.71 | 8.97 | 150 |
| patch | eig_mne_log | 13.60 | 14.68 | 150 |
| patch | eig_dspm_log | 15.67 | 16.41 | 150 |
| distributed | wmne | 17.94 | 12.43 | 150 |
| distributed | dspm | 15.93 | 11.97 | 150 |
| distributed | sloreta | 10.15 | 8.14 | 150 |
| distributed | eig_mne_log | 20.60 | 17.56 | 150 |
| distributed | eig_dspm_log | 23.21 | 18.47 | 150 |

## Paired Wilcoxon (eig vs standard), LocError

| regime | eig | vs | median diff (mm) | n | p |
|---|---|---|---|---|---|
| focal | eig_mne_log | wmne | +1.96 | 150 | 7.742e-06 |
| focal | eig_mne_log | dspm | +7.13 | 150 | 1.054e-08 |
| focal | eig_mne_log | sloreta | +15.63 | 150 | 2.168e-25 |
| focal | eig_dspm_log | wmne | +1.75 | 150 | 3.366e-07 |
| focal | eig_dspm_log | dspm | +6.94 | 150 | 6.208e-10 |
| focal | eig_dspm_log | sloreta | +16.69 | 150 | 7.296e-26 |
| patch | eig_mne_log | wmne | +0.30 | 150 | 0.008055 |
| patch | eig_mne_log | dspm | +1.57 | 150 | 0.05047 |
| patch | eig_mne_log | sloreta | +6.28 | 150 | 1.616e-16 |
| patch | eig_dspm_log | wmne | +2.00 | 150 | 2.476e-07 |
| patch | eig_dspm_log | dspm | +3.70 | 150 | 0.000969 |
| patch | eig_dspm_log | sloreta | +8.83 | 150 | 4.233e-19 |
| distributed | eig_mne_log | wmne | +0.08 | 150 | 0.1938 |
| distributed | eig_mne_log | dspm | +7.04 | 150 | 0.0003401 |
| distributed | eig_mne_log | sloreta | +9.59 | 150 | 2.259e-15 |
| distributed | eig_dspm_log | wmne | +2.51 | 150 | 0.0002458 |
| distributed | eig_dspm_log | dspm | +8.18 | 150 | 1.909e-06 |
| distributed | eig_dspm_log | sloreta | +10.70 | 150 | 1.606e-19 |
