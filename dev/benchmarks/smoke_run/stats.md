# Eigenmode accuracy benchmark - statistics

## Median LocError (mm) by method x regime

| regime | method | median | IQR | n |
|---|---|---|---|---|
| focal | wmne | 23.91 | 9.95 | 4 |
| focal | dspm | 9.44 | 1.05 | 4 |
| focal | sloreta | 4.54 | 9.14 | 4 |
| focal | eig_mne_log | 98.43 | 18.65 | 4 |
| focal | eig_dspm_log | 98.15 | 38.95 | 4 |

## Paired Wilcoxon (eig vs standard), LocError

| regime | eig | vs | median diff (mm) | n | p |
|---|---|---|---|---|---|
| focal | eig_mne_log | wmne | +80.40 | 4 | 0.06789 |
| focal | eig_mne_log | dspm | +87.70 | 4 | 0.06789 |
| focal | eig_mne_log | sloreta | +93.89 | 4 | 0.06789 |
| focal | eig_dspm_log | wmne | +74.25 | 4 | 0.06789 |
| focal | eig_dspm_log | dspm | +88.85 | 4 | 0.06789 |
| focal | eig_dspm_log | sloreta | +93.62 | 4 | 0.06789 |
