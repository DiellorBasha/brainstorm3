# ESI benchmark framework (v1) — design

**Date:** 2026-06-02
**Author:** Diellor Basha (design captured with Claude)
**Status:** Design — pending review before implementation plan
**Repo:** `research/code/brainstorm3` (branch off `development`)

---

## 1. Motivation

We need a **principled, reusable instrument** to benchmark, validate, and compare
source-estimation methods (the eigenmode/GBF method under test vs established
inverses) with statistical results — *before* further method development. The
existing `tutorial_eigenmodes_validation.m` is a one-off diagnostic prototype
(single subject, single time point, no statistics, wrong regime, limited
comparators). It revealed the regime problem but cannot *prove* anything.

This framework replaces the prototype's role as the benchmark of record while
**keeping the prototype alongside** (not deleted).

### Why our earlier benchmark was inconclusive (settled)
We tested focal, near-noiseless, single-snapshot MEG — the worst regime for a
smoothness prior. GBF's reported wins come from distributed sources under colored
noise at moderate SNR, scored on fidelity metrics. The framework must control
source structure, noise/SNR, and metrics so comparisons are fair and meaningful.

---

## 2. Scope (v1)

### In scope
- Ground truth from **Brainstorm-native simulation** (`process_simulate_*` /
  direct forward) — full control over source structure, SNR, and repetitions.
- **MEG only**, on the **OMEGA** anatomy (2 subjects).
- Static **spatial** regimes: focal point, small patch, distributed (smooth) patch.
- Comparator panel: native **wMNE, dSPM, sLORETA, LCMV** + our **eigenmode
  variants** (mne/dspm/sloreta × prior log/flat/power); **eLORETA only if
  FieldTrip is installed**.
- **Descriptive statistics**: medians, IQR, bootstrap 95% CIs, paired
  difference distributions. (Formal significance tests are a later layer.)
- Metric suite: localization error, AUC, NRMSE, correlation, spatial dispersion.
- Reusable, seeded, reproducible, with committed report artifacts.

### Explicitly deferred (architecture built to extend)
Significance/hypothesis tests · temporal/propagating sources + time-resolved
metrics (the hook for the L2/L3 spatiotemporal-prior work) · EEG modality ·
eLORETA/MxNE without FieldTrip · GBF-dataset replication · phantom data.

### Comparator availability (confirmed in code)
`bst_inverse_linear_2018` provides `minnorm`={amplitude(wMNE), dspm2018, sloreta},
`gls`, `lcmv`. **eLORETA** exists only via FieldTrip (`process_ft_sourceanalysis`).
**MxNE** is not available natively → deferred.

---

## 3. Architecture — composable library + thin driver

Five focused units (each independently testable) + a thin orchestrating driver.

```
 bst_benchmark_sources   (ground-truth generator: regime → GT vertex map(s) + time course + seed)
        │
        ▼
 bst_benchmark_simulate  (GT → sensor data; colored noise from noise cov @ target SNR)
        │
        ▼
 bst_benchmark_inverse   (comparator panel → vertex source estimate per method, matched configs)
        │
        ▼
 bst_benchmark_metrics   (GT vs estimate → LocError, AUC, NRMSE, correlation, spatial dispersion)
        │
        ▼
 bst_benchmark_report    (aggregate over realizations × regime × SNR × method →
                          medians/IQR/bootstrap CIs + paired diffs → CSV + md + PNG)

 tutorial_benchmark_esi.m  (driver: defines the sweep, calls the above, writes the report)
```

### Files

| File | Action | Responsibility |
|---|---|---|
| `toolbox/math/bst_benchmark_sources.m` | Create | Parameterized GT generator (focal / patch / distributed), seeded, N locations. |
| `toolbox/math/bst_benchmark_simulate.m` | Create | Forward-project GT + add colored noise at target SNR (reuses the noise-cov eigendecomposition used by `process_simulate_recordings`). |
| `toolbox/inverse/bst_benchmark_inverse.m` | Create | Run the comparator panel with matched configs; return vertex estimates. |
| `toolbox/inverse/bst_benchmark_metrics.m` | Create | Metric suite (extends `bst_resolution_metrics`). |
| `toolbox/math/bst_benchmark_report.m` | Create | Aggregate → descriptive stats + paired diffs → CSV/md/PNG. |
| `toolbox/script/tutorial_benchmark_esi.m` | Create | Thin driver / sweep definition. |
| `tutorial_eigenmodes_validation.m` | Keep | Prototype stays alongside (not superseded). |
| `dev/tests/test_benchmark_*_pure.m`, `test_benchmark_esi_e2e.m` | Create | Unit + smoke tests. |

---

## 4. Component contracts

### 4.1 `bst_benchmark_sources`
Given a cortex surface (+ atlas/`VertConn` for geodesic radius) and a regime spec,
produce N seeded ground-truth realizations:
- **focal**: single seed vertex active.
- **patch**: vertices within a geodesic radius `r` of the seed, uniform amplitude.
- **distributed**: smooth Gaussian profile on the surface (amplitude falls off with
  geodesic distance, scale `σ`) — the GBF-favorable regime.
Returns, per realization: GT vertex amplitude map `g [nVert×1]`, a scalar time
course `c [1×nTime]` (simple Gaussian-windowed burst), the GT source matrix
`G = g·c [nVert×nTime]`, and the seed vertex index (for LocError). RNG seeded.

### 4.2 `bst_benchmark_simulate`
`F = L · G(:,t)` over good channels via the head-model `Gain`, then add **colored
noise** drawn from the study noise covariance `C` via eigendecomposition
(`V·√D·randn`, the `get_noise_signals` model), scaled so the achieved sensor SNR
matches the target in dB: `noise_power = signal_power / 10^(SNR/10)` (GBF's
definition). Returns simulated sensor data + the noise cov used. M noise draws per
GT realization.

### 4.3 `bst_benchmark_inverse`
Run each method on the simulated data with matched inputs (same good channels,
noise cov, eval time = peak of the time course):
- **Native** via `bst_inverse_linear_2018` on the base head model: wMNE
  (`minnorm`/`amplitude`), dSPM (`dspm2018`), sLORETA; **LCMV** (`lcmv`, needs a
  data covariance computed from the simulated data).
- **Eigenmode** via `bst_inverse_eigenmodes` on the composed head model: mne/dspm/
  sloreta × prior {log, flat, power}; vertex estimate `Φ·M̃`.
- **eLORETA** via `process_ft_sourceanalysis` **iff FieldTrip is on the path**;
  otherwise logged as skipped.
Returns a vertex estimate per method. `λ`/SNR tuned per method for fairness.

### 4.4 `bst_benchmark_metrics`
GT map vs estimate (at eval time):
- **LocError** (mm): seed→estimate-peak distance; top-K Hungarian-matched for
  multi-source patches.
- **AUC**: ROC for detecting GT-active vertices (`|g|>thr`) from estimate magnitude.
- **NRMSE**: norm-matched RMSE / GT range (GBF definition).
- **Correlation**: Pearson(GT, estimate).
- **Spatial dispersion** (mm): from `bst_resolution_metrics`.

### 4.5 `bst_benchmark_report`
Long-format table (realization, subject, regime, SNR, method, metric, value) →
per (regime × SNR × method): **median, IQR, bootstrap 95% CI**; plus the
**paired per-realization difference** (eigenmode − each comparator) with its CI
(the input the later significance layer will consume). Writes CSV + markdown
summary tables + PNG figures (metric vs SNR per method per regime; paired-diff
distributions) to `dev/benchmarks/<date>/`. Seeded.

---

## 5. Evaluation design (v1)

- **Regimes:** focal · small patch (`r≈5–10 mm`) · distributed (`σ≈15–25 mm`).
- **SNR sweep:** {0, 3, 6, 10} dB (colored noise from the OMEGA noise cov).
- **Repetitions:** N≈20–40 seed locations spread across cortex × M≈5 noise draws
  each → distributions and CIs (not point estimates).
- **Subjects:** OMEGA sub-0002 and sub-0003 (MEG).
- **Eval:** at the time course peak (static-spatial v1).

---

## 6. Statistics (v1 — descriptive)

Per (regime × SNR × method): median + IQR + bootstrap 95% CI for each metric.
Per (regime × SNR): the paired per-realization difference (eigenmode − comparator)
distribution + its bootstrap CI. **No p-values in v1** — but the paired structure
is produced so the significance layer (Wilcoxon signed-rank / permutation, with
multiple-comparison control) drops in without re-running simulations.

---

## 7. Error handling

- Missing head model / eigenmodes / noise cov → clear error, skip that
  subject/regime with a logged reason (never crash the whole sweep).
- FieldTrip absent → eLORETA skipped, logged.
- LCMV with a rank-deficient data cov → regularize (Brainstorm convention) or skip
  with a logged note.
- Degenerate metric inputs (all-zero estimate) → metric returns a defined sentinel,
  not NaN propagation.

## 8. Testing

- **Pure:** `test_benchmark_sources_pure` (focal→one active vertex at seed; patch→
  active set within radius; distributed→monotone falloff with distance; shapes,
  seeding reproducibility). `test_benchmark_metrics_pure` (perfect estimate → LE 0,
  corr 1, AUC 1, NRMSE 0; a known-shifted estimate → LE = the known mm distance).
  `test_benchmark_simulate_pure` (achieved SNR ≈ target within tolerance; sample
  noise covariance ≈ input cov).
- **e2e smoke:** `test_benchmark_esi_e2e` — tiny sweep (1 regime × 1 SNR × 2 reps ×
  2 methods) on OMEGA → report artifacts produced, table shapes correct. Skips
  cleanly without a suitable protocol.

## 9. Reproducibility

Seeded `rng` (seed recorded in the report). Results to `dev/benchmarks/<date>/`
(CSV + markdown + PNG), re-runnable end-to-end via the MATLAB MCP. The prototype
`tutorial_eigenmodes_validation.m` is retained unchanged.

## 10. Open implementation details (resolved during planning)

- Geodesic radius/profile: exact `VertConn`/distance method for patch extent.
- Bootstrap CI: number of resamples; per-metric handling of bounded metrics (AUC).
- LCMV data-covariance estimation window on simulated single-burst data.
- Whether `bst_benchmark_simulate` calls `process_simulate_recordings` directly for
  uniform patches vs the direct forward+noise path for smooth profiles (the direct
  path is preferred for arbitrary vertex profiles; confirm it reproduces the
  process's noise model exactly).
