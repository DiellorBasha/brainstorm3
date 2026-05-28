# Eigenmode Complex Wavelet Tensor — Milestone 3 (increment 1) Design

**Date:** 2026-05-27
**Branch:** `feature/eigenmode-wavelet` (off `development`, which already has M1 transform + M2 noise floor)
**Status:** Approved design, ready for implementation plan

## Background

M1 produces raw eigenmode coefficient time series `θₖ(t)` (`matrix_eigentransform`) and a stationary `(λ,ω)` spectrum via FFT. M2 adds an empty-room noise floor. The next step in the proposals (spatiotemporal filter banks §4) is the **complex spatiotemporal wavelet tensor** `Wₖ(s,t)` — the time-resolved `(λ, ω, t)` decomposition that carries both amplitude and phase, and is the data structure the later analyses (dispersion, inter-eigenmode phase coherence, cross-frequency coupling, instantaneous dispersion) build on.

Brainstorm already has a validated complex Morlet CWT (`toolbox/timefreq/morlet_transform.m`) and `process_timefreq` accepts `matrix` inputs. So the tensor is computed by Morlet-transforming the eigenmode coefficients — the same "reuse an existing transform on the coefficients" pattern M1 used with `process_fft`. This increment wraps that in a pure, unit-tested core plus a thin process, deferring all analyses.

See [[eigenmode-transform-first]] and [[eigenmode-meg-project]].

## Goal (this increment, scoped narrowly)

Compute and store the **complex** wavelet tensor `Wₖ(s,t)` of the eigenmode coefficients:

1. A pure function that Morlet-transforms `[K × nTime]` coefficients into a complex `[K × nTime × nFreq]` tensor (amplitude `|W|`, phase `arg(W)`).
2. A thin process that applies it to a `matrix_eigentransform` file and saves a complex, λ-labeled timefreq file, viewable in Brainstorm's standard time-frequency viewer.

## Non-goals (YAGNI — deferred to later M3 increments)

- Dispersion analysis (`ω = c√λ` ridge fit, empirical propagation speed).
- Wave-vs-diffusion (dispersion-curve vs wedge) discrimination.
- Physiological speed-limit mask.
- Inter-eigenmode phase coherence, cross-frequency coupling, instantaneous (time-varying) dispersion.
- Optical-flow wave detection.
- Any new Morlet math (we wrap the existing `morlet_transform`, not reimplement it).

## The math / mechanics

For each eigenmode coefficient series `θₖ(t)`, the complex Morlet CWT gives

```
Wₖ(f,t) = ∫ θₖ(τ) ψ*_{f}(τ − t) dτ ∈ ℂ
```

implemented by `morlet_transform(Coeffs, t, Freqs, Fc, FWHM_tc, 'n')` with `'n'` = un-squared (complex) coefficients. **`morlet_transform` permutes internally and returns `[K × nTime × nFreq]` directly** (its own header comment claiming `[K × nFreq × nTime]` is wrong), so **no extra permute is applied**. We apply `conj()` so phase follows the standard positive-rotation analytic-signal convention (`morlet_transform`'s native convention is the conjugate); `|conj(W)| = |W|`, so amplitude is unaffected. `|W|` is amplitude, `arg(W)` is phase. Squaring (`'y'`) would discard phase — we explicitly do not, since phase is the point.

## Component contracts

### `toolbox/math/bst_eigenmodes_wavelet.m` (new — pure)

```matlab
[W, Freqs] = bst_eigenmodes_wavelet(Coeffs, sfreq, Freqs, varargin)
```

- **Inputs**
  - `Coeffs` `[K × nTime]`: eigenmode coefficient time series.
  - `sfreq`: sampling rate (Hz).
  - `Freqs`: `[1 × nFreq]` frequencies (Hz). If empty `[]`, a default **log-spaced** grid is built: `logspace(log10(2), log10(min(100, 0.4*sfreq)), 40)`.
- **Options** (name-value): `'MorletFc'` (default 1), `'MorletFwhmTc'` (default 3) — Brainstorm Morlet parameters.
- **Computation**
  - `t = (0:nTime-1) / sfreq;`
  - if `Freqs` empty → build default grid (above).
  - `W = conj(morlet_transform(double(Coeffs), t, Freqs, MorletFc, MorletFwhmTc, 'n'));`  % already `[K × nTime × nFreq]` (morlet permutes internally); conj → standard phase
- **Outputs**: complex `W` `[K × nTime × nFreq]`; the `Freqs` actually used (row vector).
- **Invariants** (unit-tested): a single mode carrying a pure sinusoid at `f₀` (others zero) yields `|W|` peaking at that `(mode, f₀)` cell; phase at `f₀` advances at ≈ `2πf₀` rad/s; other modes' amplitude ≈ 0; `W` is complex; empty `Freqs` yields a valid `[K × nTime × 40]` tensor.

Pure — no file I/O.

### `toolbox/process/functions/process_eigenmodes_wavelet.m` (new — thin)

- `GetDescription`: Comment `'Eigenmode wavelet tensor (complex)'`, Category `'Custom'`, SubGroup `'Sources'`, Index `336.7` (after denoise 336.6, verified free), InputTypes `{'matrix'}`, OutputTypes `{'timefreq'}`, nInputs 1, nMinFiles 1. Options:
  - `flo` (value, default `2`, Hz), `fhi` (value, default `0`, Hz; `0` ⇒ auto `min(100, 0.4*sfreq)`), `nfreqs` (value, default `40`).
  - `morletfc` (value, default `1`), `morletfwhmtc` (value, default `3`).
  - `label_info` (label): notes the output is the complex `(λ,ω,t)` tensor (amplitude + phase); input is an eigenmode-coefficient matrix.
- `Run`:
  1. Load the matrix file (`in_bst_matrix` / `in_bst_data` for the `matrix` type); `Coeffs = Value` `[K × nTime]`, `Time`, `Description` (mode labels), `sfreq = 1/(Time(2)-Time(1))`.
  2. Build `Freqs`: `fhiEff = fhi>0 ? fhi : min(100, 0.4*sfreq)`; `Freqs = logspace(log10(flo), log10(fhiEff), nfreqs)`.
  3. `[W, Freqs] = bst_eigenmodes_wavelet(Coeffs, sfreq, Freqs, 'MorletFc',morletfc, 'MorletFwhmTc',morletfwhmtc)`.
  4. Save a timefreq file (`db_template('timefreqmat')`): `TF = W` (complex `[K × nTime × nFreq]`), `Time`, `Freqs`, `RowNames = Description`, `Measure = 'none'` (complex), `Method = 'morlet'`, `DataType = 'matrix'`, `SurfaceFile` (carried from the matrix file if present), `Comment`, history. `db_add_data`. Prefix `timefreq_eigenwavelet`.
  5. Report `K`, `nFreq`, freq range via `bst_report`.

## Verification (this increment)

1. **Pure unit test** `dev/tests/test_eigenmodes_wavelet_pure.m` (DB-free): the single-mode single-frequency synthetic above — assert tensor shape `[K × nTime × nFreq]`, complexity, amplitude peak at the correct `(mode, f₀)`, negligible leakage to other modes, phase-advance rate ≈ `2πf₀/sfreq` in the central (edge-effect-free) region, and that empty `Freqs` produces the default 40-frequency grid.
2. **checkcode** clean on both files (Brainstorm idioms excepted).
3. **Optional manual smoke** on the `EigenSmoke` substrate: run `process_eigenmodes_wavelet` on the existing `matrix_eigentransform` coefficients → confirm a finite complex `[K × nTime × nFreq]` tensor; open it in the TF viewer and confirm a time-resolved `(λ,ω)` amplitude image (e.g. an alpha band that waxes/wanes over time). Save a figure.

## Risks

- **Output orientation (resolved during implementation):** `morlet_transform` returns `[K × nTime × nFreq]` (it permutes internally; its header comment is wrong), so the code does NOT permute again and applies `conj()` for standard phase rotation. The unit test's shape + peak-cell assertions guard this.
- **Complex timefreq handling:** `Measure='none'` keeps complex TF; confirm `db_add_data` + the TF viewer accept a complex `matrix`-type timefreq (the optional smoke validates display). If the viewer needs a real measure for display, magnitude view is the fallback; the stored tensor stays complex.
- **Memory:** `W` is complex `[K × nTime × nFreq]` — for `K=200`, `nTime=60 s × sfreq`, `nFreq=40` this is large. The default `nfreqs=40` and the user's choice of block length bound it; document that long blocks × many modes × many freqs are memory-heavy (a sparse/els strategy is out of scope here).
- **Branch base:** depends on M1's `matrix_eigentransform` output (already on `development`).
