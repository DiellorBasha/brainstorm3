# Eigenmode Source-Mapping Accuracy Benchmark — Design

**Date:** 2026-06-02
**Author:** Diellor Basha (with Claude)
**Status:** Design — pending review before implementation plan

## Goal

Demonstrate that the LBO-eigenmode source-mapping method (GBF-style: `eig_mne/log`,
`eig_dspm/log`) achieves **localization accuracy competitive with the standard MNE
family** (`wMNE`, `dSPM`, `sLORETA`) against ground truth, across focal, patch, and
distributed source regimes and a range of SNRs. Deliverables are **statistics** (tables
+ significance tests) and **MATLAB figures** illustrating the comparison.

A secondary, measured question folds in: **does increasing eigenmode bandwidth (K)
recover focal-source accuracy, or does the spectral prior cap it?** This is answered
empirically via a K-sweep rather than asserted.

## Scope decisions (locked)

| Decision | Choice | Rationale |
|---|---|---|
| Primary claim | Competitive accuracy vs ground truth (head-to-head) | User goal |
| Data track | **Synthetic sources on real cortex** only | Surface-native ground truth on the *actual* cortical Laplacian; sources sit on the surface every method solves on → fair, no depth floor. Exercises the geometry-aware prior that is the method's defining strength. |
| Phantom (Elekta) | **Dropped** | Its 32 dipoles lie on two perpendicular intersecting planes (volume source space); no single 2-manifold surface contains them, so the surface-based eigenmode method cannot run on it without an artificial floor. Its only irreplaceable value (real fields from known dipoles) is inherently volumetric. |
| Anatomies | Auditory (CTF) + Neuromag (Elekta) | Two real head models, two MEG vendor systems, real noise covariance |
| Methods (5) | `wMNE`, `dSPM`, `sLORETA`, `eig_mne/log`, `eig_dspm/log` | All imaging methods producing full cortical maps → scored identically (apples-to-apples) |
| Regimes (3) | focal, patch, distributed | focal = eigenmode weak point, distributed = eigenmode strong point; "competitive across all three" is the claim |
| SNR sweep | `[2 4 6 10 20]` dB | Span low→high; shows accuracy holds as noise rises |
| Replicates | ~15 per (anatomy, regime, SNR, K) | Distributions + paired significance tests |
| Eigenmode bandwidth | **Sweep K ∈ {300, 600, 1000} per hemisphere** | Turns "more bandwidth fixes focal" into a measured curve; compute 1000 once, subset for free. Note data rank ceiling: Auditory MEG rank ≈ 272, Neuromag ≈ 205 — modes beyond the leadfield's effective bandwidth are prior-dominated, so the curve is expected to plateau. |

## Reused primitives (already implemented, tested in `dev/tests/`)

| Function | Role | Verified by |
|---|---|---|
| `bst_benchmark_sources(Surface, regime, ...)` → `.GT, .Sources, .SeedVertex` | Ground-truth generator (focal/patch/distributed) | `test_benchmark_sources_pure.m` |
| `bst_benchmark_simulate(L, Sources, NoiseCov, 'SNR', s, 'Seed', k)` → `.F` | Forward-project + add noise to target SNR | `test_benchmark_simulate_pure.m` |
| `bst_benchmark_inverse(F, baseHmFile, ncFile, chFile, goodMask, SNR)` → `Est.<method>` | Runs all methods on a base surface head model + cortex eigenmodes; returns `[nVert × nTime]` per method | `test_benchmark_inverse_e2e.m` |
| `bst_benchmark_metrics(gt, est, GridLoc, seedVertex)` → `.LocError, .Correlation, .NRMSE, .AUC, .SpatialDispersion` | Metric suite | `test_benchmark_metrics_pure.m` |

**Only extension required:** parameterize the eigenmode arms of `bst_benchmark_inverse`
by **K (modes/hemisphere)** so the same simulated data can be inverted at K ∈ {300,600,1000}
without recomputing the eigendecomposition (compute 1000 once on the cortex, subset).

## New modules (plain scripts under `dev/benchmarks/`)

1. **`bench_config.m`** — one config struct: anatomies (protocol/subject), methods, regimes,
   SNR grid, K grid, replicate count, master RNG seed, output dir. Plus a `smoke` preset
   (1 anatomy, focal only, 2 SNRs, 2 reps, K={300}) for ~1–2 min pipeline validation.
2. **`bench_run.m`** — the synthetic driver. For each anatomy: load base surface head model
   gain + cortex + noise covariance; compute 1000 eigenmodes/hemisphere **once**. Then loop
   regime × SNR × replicate: `sources → simulate → inverse(all methods; eig arms at each K)
   → metrics`. Emit one tidy row per (anatomy, regime, SNR, replicate, method, K).
3. **`bench_stats.m`** — aggregate rows → median/IQR per (method, regime, SNR, K); paired
   **Wilcoxon signed-rank** (each eigenmode arm vs each standard method, on matched
   replicates). Write `stats.csv` + a markdown summary table.
4. **`bench_figures.m`** — produce the 5 figures below (PNG + `.fig`).
5. **`benchmark_eigenmodes.m`** — top-level entry: `config → bench_run → bench_stats →
   bench_figures → REPORT.md`.

## Results schema (`synthetic.csv`)

`anatomy, regime, snr_db, replicate, method, K, locerror_mm, correlation, nrmse, auc, spatial_dispersion_mm`

`K` is `NaN`/`full` for the three standard methods (no mode basis).

## Statistics

- Per (method, regime, SNR, K): **median + IQR** of each metric over replicates.
- **Paired Wilcoxon signed-rank** test: each eigenmode arm vs each standard method, matched
  on (anatomy, regime, SNR, replicate). Reported as p-value + median paired difference, so
  "competitive" is supported by *failure to reject* / small effect size, not just overlapping
  boxes. Eigenmode arm uses its plateau-K (see below).
- **Plateau-K selection:** from the K-sweep on the focal regime, pick the smallest K beyond
  which median focal LocError no longer improves meaningfully (≤ ~1 mm). Main comparison
  figures (1–4) use this K; the choice is justified by Figure 5.

## Figures

| # | Figure | Content |
|---|---|---|
| 1 | **Distribution** | Box/violin of `locerror_mm` per method (at plateau-K), pooled over anatomies/SNR; the headline accuracy comparison |
| 2 | **SNR sweep** | `locerror_mm` (mean ± error bar) vs SNR per method; one panel per regime |
| 3 | **Example reconstructions on cortex** | Ground-truth seed vs each method's estimate rendered on the real cortical surface for a representative focal case; the qualitative "where did it land" panel (Brainstorm surface rendering) |
| 4 | **Per-regime breakdown** | Grouped bars: median `locerror_mm` by regime × method, IQR whiskers |
| 5 | **K-sweep curve** | `locerror_mm` vs K ∈ {300,600,1000} for the eigenmode arms (focal regime emphasis), with standard-method medians as horizontal reference lines; shows whether focal accuracy keeps improving or plateaus near the data rank |

## Outputs

All under `dev/benchmarks/eigenmode_accuracy_<YYYYMMDD>/`:
`synthetic.csv`, `stats.csv`, `figures/*.png` (+ `.fig`), `REPORT.md` (tables + figure
references + plateau-K finding + competitiveness verdict).

## Error handling & reproducibility

- Each (condition, method, K) wrapped in try/catch → logged `NaN` row; never aborts the run.
- Surface manifold gating via existing `tess_manifold` (eigenmodes require a 2-manifold).
- Deterministic per-replicate RNG seed derived from the master seed; same config → same numbers.
- Smoke preset must pass before a full run.

## Runtime estimate

Eigendecomposition (1000 modes/hemi) is the one-time cost (~1 min/surface). All inversions
are precomputed linear kernels (matrix multiplies). Full grid ≈ 2 anatomies × 3 regimes ×
5 SNR × 15 reps × (3 standard + 2 eig × 3 K) = 9 method-configs ≈ 4,050 inversions →
minutes, not hours.

## Risks / open questions

- **Prior may cap focal resolution regardless of K.** The `log` prior downweights high-λ
  modes. If Figure 5 shows focal LocError plateauing well above the standard methods, the
  spectral prior — not the mode count — is the binding constraint. Follow-up (out of scope
  here): compare `log` vs `flat`/tempered prior on the focal regime.
- **Data-rank ceiling.** Modes beyond MEG rank (~272 Auditory, ~205 Neuromag) are
  prior-constrained, not data-constrained; the K-sweep is expected to plateau near there.
- **`bst_benchmark_inverse` K-parameterization** is the one new code path; must be unit-tested
  (subset-of-modes inverse equals full recompute at that K).

## Out of scope

- Real evoked-data validation (Auditory/Neuromag *measured* responses vs expected location) —
  no exact ground truth; excluded from an accuracy claim.
- Phantom / volume source spaces.
- Vector eigenmodes, dispersion/wavelet (dynamics) analyses.
- Prior-family comparison beyond the focal follow-up noted above.
