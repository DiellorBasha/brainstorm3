# Eigenmode Dispersion Analysis (wave vs diffusion) — M3 increment 2 Design

**Date:** 2026-05-27
**Branch:** `feature/eigenmode-dispersion` (off `development`, which has M1 transform + M2 noise floor + M3-inc1 wavelet tensor)
**Status:** Approved design, ready for implementation plan

## Background

The eigenmode pipeline now produces stationary joint `(λ,ω)` power spectra `P(k,f)` — from M1's `process_fft` of the coefficients, or M2's cleaned-PSD / SNR spectra — and (M3 inc 1) the time-resolved complex wavelet tensor. The proposals (spatiotemporal filter banks §3.2) argue that *dynamical regime* shows up as a characteristic signature in the `(λ,ω)` plane:

- **Traveling waves** concentrate energy on a **dispersion curve** `ω = c√λ` (peak frequency rises with `√λ`); the slope gives propagation speed `c`.
- **Heat diffusion** fills a **wedge**: each mode `k` has a Lorentzian temporal spectrum of half-width `∝ λ_k`, so spectral *bandwidth* rises with `λ` (not a single ridge).

This increment estimates those signatures from the stationary spectrum and discriminates the regime. It runs on the small `P(k,f)` (K×nFreq), so it sidesteps the deferred wavelet-tensor memory problem.

See [[eigenmode-meg-project]].

## Goal (this increment)

Given a stationary `(λ,ω)` power spectrum and the eigenvalues `λ_k`, fit a **wave** model and a **diffusion** model and report which regime fits better, with the associated physical parameter:

1. **Wave fit** → propagation speed `c`, `R²_wave`.
2. **Diffusion fit** → diffusivity `α`, `R²_diff`.
3. **Regime** = the better-fitting model + a discrimination margin.

## Non-goals (YAGNI — deferred)

- Physiological speed-limit mask `ω ≤ c_max√λ` + admissible-power fraction (deferred by decision).
- Time-varying / instantaneous dispersion from the full wavelet tensor (this increment is stationary).
- Full 2-D power-template fitting / Bayesian model comparison (the feature-based discrimination below is the chosen, robust realization).
- Phase coherence, cross-frequency coupling, optical flow, vector work.

## Method — feature-based two-model discrimination

For each eigenmode `k`, reduce its spectrum `P(k,·)` to two scalar features, then test how each scales with `λ_k`:

- **Peak frequency** `f*(k) = argmax_f P(k,f)`.
- **Spectral bandwidth** `w(k)` = power-weighted standard deviation of frequency:
  `fbar(k)=Σ_f P·f / Σ_f P`, `w(k)=sqrt(Σ_f P·(f−fbar)² / Σ_f P)`.
- **Per-mode weight** `p_k = Σ_f P(k,f)` (total power) — used to down-weight modes with little/no signal in the fits.

**Wave model** (`ω = c√λ` ⇒ `f* = (c/2π)·√λ`): weighted through-origin fit of `f*(k)` on `√λ_k`:
`a = Σ_k p_k·√λ_k·f*_k / Σ_k p_k·λ_k`, propagation speed **`c = 2π·a`**; `R²_wave` = weighted squared correlation of `f*` with `√λ`.

**Diffusion model** (Lorentzian half-width `∝ λ` ⇒ `w ∝ λ`): weighted through-origin fit of `w(k)` on `λ_k`:
`b = Σ_k p_k·λ_k·w_k / Σ_k p_k·λ_k²`, diffusivity **`α = 2π·b`** (proportional); `R²_diff` = weighted squared correlation of `w` with `λ`.

**Regime:** `wave` if `R²_wave ≥ R²_diff`, else `diffusion`; discrimination margin `ΔR² = |R²_wave − R²_diff|`.

Rationale: a wave packs power into a peak that *moves up* with `√λ` while keeping bandwidth roughly fixed; diffusion keeps the peak near DC while *broadening* with `λ`. The two features (peak location vs width) are exactly the distinguishing signatures, and weighted `R²` gives a symmetric, robust goodness for the comparison.

**Units:** `λ_k` are LBO eigenvalues on a metre-scale mesh (Brainstorm surfaces are in metres), so `√λ` is in rad/m and `c = 2π·a` is in **m/s**. `α` is reported in the eigenvalue-consistent units (proportional to diffusivity). Documented as a dependency on the surface being in metres.

## Component contracts

### `toolbox/math/bst_eigenmodes_dispersion.m` (new — pure)

```matlab
Out = bst_eigenmodes_dispersion(P, lambdas, Freqs, varargin)
```

- **Inputs**
  - `P` `[K × nFreq]`: non-negative `(λ,ω)` power per eigenmode per frequency.
  - `lambdas` `[K × 1]`: eigenvalues (rad²/m²), aligned with `P`'s rows.
  - `Freqs` `[1 × nFreq]`: frequencies (Hz).
- **Options** (name-value): `'MinPowerFrac'` (default `0`): drop modes whose total power is below this fraction of the max per-mode power before fitting (0 = keep all, power-weighting handles it).
- **Computation:** per-mode `p_k`, `f*(k)`, `w(k)`; the two weighted through-origin fits and weighted `R²` above; only modes with `lambdas>0` and finite features are used.
- **Outputs** (`Out` struct):
  - `.PeakFreq` `[K×1]`, `.Bandwidth` `[K×1]`, `.Weights` `[K×1]`
  - `.c` (m/s, wave speed), `.alpha` (diffusivity), `.R2wave`, `.R2diff`
  - `.Regime` (`'wave'` | `'diffusion'`), `.Margin` (`|R²_wave − R²_diff|`)
- **Invariants** (unit-tested): a synthetic **wave** spectrum (Gaussian peak at `f*=c₀√λ/2π`, fixed width) → `Regime='wave'`, `R²_wave>R²_diff`, `c≈c₀` (±15%); a synthetic **diffusion** spectrum (per-mode Lorentzian half-width `α₀λ`) → `Regime='diffusion'`, `R²_diff>R²_wave`. Pure — no file I/O.

### `toolbox/process/functions/process_eigenmodes_dispersion.m` (new — thin)

- `GetDescription`: Comment `'Eigenmode dispersion (wave vs diffusion)'`, Category `'Custom'`, SubGroup `'Sources'`, Index `336.8` (after the wavelet tensor 336.7, verified free), InputTypes `{'timefreq'}`, OutputTypes `{'matrix'}`, nInputs 1, nMinFiles 1. Option `minpowerfrac` (value, default 0); `label_info` (label).
- `Run`:
  1. Load the timefreq (`in_bst_timefreq`): `TF`, `Freqs`, `Measure`, `SurfaceFile`. Convert to power **per the file's `Measure`** (so already-power inputs aren't squared again): `'none'` (complex, e.g. the M3 tensor or complex FFT) → `|TF|²`; `'magnitude'` → `TF.^2`; `'power'` (e.g. M2 cleaned-PSD/SNR, or PSD-unit FFT) → `TF`. Then average over time → `P [K × nFreq]` (handles both stationary `[K×1×nFreq]` spectra and the `[K×nTime×nFreq]` wavelet tensor).
  2. `λ_k` from `in_tess_eigenmodes(SurfaceFile).Values(1:K)` (error if absent or vertex/row mismatch).
  3. `Out = bst_eigenmodes_dispersion(P, lambdas, Freqs, 'MinPowerFrac', minpowerfrac)`.
  4. `bst_report('Info', …)` with `Regime`, `c`, `α`, `R²_wave`, `R²_diff`, margin.
  5. Save a small `matrix` file (`matrix_eigendispersion`): `Value = [PeakFreq'; Bandwidth']` `[2 × K]`, `Time = 1:K` (mode index), `Description = {'PeakFreq (Hz)','Bandwidth (Hz)'}`, `SurfaceFile`, Comment carrying the regime/c/α/R² summary, history.

## Verification (this increment)

1. **Pure unit test** `dev/tests/test_eigenmodes_dispersion_pure.m` (DB-free): the synthetic wave and diffusion spectra above — assert correct `Regime`, the winning `R²` ordering, and wave-speed recovery within tolerance.
2. **checkcode** clean on both files (Brainstorm idioms excepted).
3. **Optional manual smoke**: build a `(λ,ω)` spectrum on a real (downsampled/epoched) recording — e.g. `process_eigenmodes_transform` → `process_fft` (or M2's cleaned PSD) — run `process_eigenmodes_dispersion`, and confirm it reports a regime + finite `c`/`α`/`R²`. (No fixed expected regime for real resting data; this is a runs-clean check.)

## Risks

- **Eigenvalue units / `c` interpretation:** `c` in m/s only holds if the surface is in metres (Brainstorm convention) and `lambdas` are the metric LBO eigenvalues. Documented; the regime decision and `R²` are unit-independent.
- **Ridge robustness:** flat/noisy high-`k` modes can distort the fits; mitigated by power-weighting and the optional `MinPowerFrac` cut. (A future increment can restrict to M2's reliable-mode set `K*(f)`.)
- **DC / removed modes:** modes with `λ_k ≤ 0` (e.g. removed DC) are excluded from the fits.
- **Diffusion bandwidth on a finite grid:** a Lorentzian's std is formally infinite; on the finite `Freqs` grid the power-weighted width is finite and still grows monotonically with `λ`, which is what the discrimination needs.
- **Input type:** expects a `(λ,ω)` power spectrum whose rows are eigenmodes (M1 FFT / M2 PSD / M3 tensor). A non-eigenmode timefreq would mismatch `λ_k` — guarded by the vertex/row-count check against `SurfaceFile`.
