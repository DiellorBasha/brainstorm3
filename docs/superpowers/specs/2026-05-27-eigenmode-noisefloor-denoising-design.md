# Eigenmode Noise-Floor Denoising — Milestone 2 Design

**Date:** 2026-05-27
**Branch:** `feature/eigenmode-noisefloor` (based on `feature/eigenmode-transform-spectrum`, which is not yet merged)
**Status:** Approved design, ready for implementation plan

## Background

Milestone 1 added the unregularized spatial transform `A = pinv(L·Φ)` (`bst_eigenmodes_transform`) and a process that produces raw eigenmode coefficient time series `θₖ(t)` (`process_eigenmodes_transform`), FFT-able into a joint (λ,ω) spectrum. The real-data smoke confirmed it works but also confirmed the motivating gap: the raw coefficients have **no noise reference**. High-k modes are noise-amplified and line noise (60 Hz) dominates, with no principled way to say which (k,f) cells carry signal.

Brainstorm already estimates an empty-room **noise covariance** `Σ_n` (a time-averaged *spatial* second-order statistic) for inverse whitening. This milestone extends that idea to the spatial-mode × temporal-frequency plane: treat the empty-room **recording** symmetrically with the data — push it through the **same** transform and FFT it — to obtain a frequency-resolved noise floor `N(k,f)`.

This is the transform-first philosophy applied to noise: you transform the data, you transform the noise, you compare. See [[eigenmode-transform-first]].

## Goal

Given a data recording (with a surface head model and precomputed eigenmodes) and a paired empty-room recording, compute the joint (λ,ω) noise floor and use it to:

1. **Diagnose** — an SNR-resolved spectrum `SNR(k,f) = P_data/N` over the (λ,ω) plane, plus a data-adaptive reliable-mode cutoff `K*(f)`.
2. **Denoise (spectrum)** — power spectral subtraction `P̂(k,f) = max(P_data − αN, 0)` → a cleaned (λ,ω) power spectrum.

Optionally (off by default): **Wiener-gain the signal** — `G = max(P_data−N,0)/P_data` applied to the data coefficients to yield cleaned time-domain coefficients.

## The math (and why this shape)

The data and empty-room recordings contain **different realizations** of the same noise process. For a (λ,ω) cell, `Θ_data = S + N_d` and `Θ_noise = N_e` with `N_d ⟂ N_e`.

- **Complex subtraction is wrong:** `Θ_data − Θ_noise = S + (N_d − N_e)` has `Var(N_d−N_e) = 2N` — it *doubles* the noise.
- **Power subtraction is correct (in expectation):** with averaging, the cross-term vanishes, so `E[|Θ_data|²] = |S|² + N` and `E[|Θ_noise|²] = N`. Hence `P_data − N → |S|²`. This is classic spectral subtraction; clamp at ≥ 0 because estimates are noisy.
- **Averaging is required:** a single-window `|FFT|²` is a high-variance (χ²) estimate. Both spectra must be **Welch-averaged PSDs** in density units (power/Hz) so different data/empty-room durations are comparable.
- **Wiener** is the MMSE sibling: `Θ_clean = G·Θ_data`, `G = max(P_data−N,0)/P_data ∈ [0,1]`, preserves phase and inverts back to a time series.

### Kernel handling (critical)

The transform kernel `A = pinv(L̃)` is built from the **data's** head model, eigenmodes, and good-channel set. The empty room has no head model, so we do **not** build a separate kernel for it. Instead we apply the **data's** kernel to the empty-room sensor data on the **common good-channel set** of the two recordings:

```
θ_data(t)  = A · d_data(iCommon, t)
θ_noise(t) = A · d_noise(iCommon, t)      % same A, same channels
```

If the channel sets differ, intersect them and rebuild `A` on the common set so both projections use identical rows.

## Architecture & data flow

```
data recording (head model Φ, A=pinv(LΦ)) ──┐
empty-room recording ────────────────────────┤  apply SAME kernel A to both
                                              ▼    (common good channels)
              θ_data [K×T]   θ_noise [K×Tn]
                    │              │
            Welch PSD        Welch PSD          (density units, same settings)
                    ▼              ▼
              P_data(k,f)      N(k,f)
                    └──────┬───────┘
                           ▼   bst_eigenmodes_noisefloor (pure)
        SNR(k,f) ; CleanPSD=max(P−αN,0) ; K*(f) ; [Wiener gain G(k,f)]
                           ▼
   timefreq_eigensnr (SNR over λ,ω) + timefreq_eigencleanpsd (cleaned spectrum)
                           + reliable-mode cutoff reported
              [optional] matrix_eigentransform_clean (Wiener'd coefficients)
```

## Component contracts

### `toolbox/math/bst_eigenmodes_noisefloor.m` (new — pure)

```matlab
[Out] = bst_eigenmodes_noisefloor(Pdata, Nnoise, varargin)
```

- **Inputs**
  - `Pdata` `[K × nFreq]`: data PSD (power/Hz) per eigenmode per frequency.
  - `Nnoise` `[K × nFreq]`: empty-room PSD, same `K`, same frequency grid.
- **Options** (name-value): `'Alpha'` (over-subtraction, default `1`), `'Floor'` (spectral floor fraction β of N, default `0`), `'SnrThresh'` (linear SNR for the reliable-mode cutoff, default `1`).
- **Outputs** (`Out` struct):
  - `.SNR`        `[K × nFreq]` = `Pdata ./ max(Nnoise, eps)`.
  - `.CleanPSD`   `[K × nFreq]` = `max(Pdata − Alpha.*Nnoise, Floor.*Nnoise)`.
  - `.Gain`       `[K × nFreq]` = `max(Pdata − Nnoise, 0) ./ max(Pdata, eps)` (Wiener gain ∈ [0,1]).
  - `.Kstar`      `[1 × nFreq]` = largest `k` with `SNR(k,f) ≥ SnrThresh` (0 if none).
- **Invariants** (unit-tested): with `Pdata = S + N` (`S,N ≥ 0`), `CleanPSD ≈ S` when `Alpha=1,Floor=0`; `SNR = (S+N)/N`; `Gain ∈ [0,1]` and `Gain = 0` wherever `Pdata ≤ Nnoise`; `Kstar` monotone-correct on a synthetic SNR ramp.

No file I/O, no Welch (the caller supplies PSDs) — pure combine math.

### `toolbox/process/functions/process_eigenmodes_denoise.m` (new)

- `GetDescription`: Comment `'Eigenmode noise-floor denoising'`, Category `'Custom'`, SubGroup `'Sources'`, Index `336.6` (just after the transform 336.5), InputTypes `{'data','raw'}`, nInputs 1, nMinFiles 1. Options:
  - `noisefile` (`filename` selector): the empty-room recording (raw or imported).
  - `nmodes` (value, default `0` = auto = `min(nch, available)`).
  - `alpha` (value, default `1`), `snrthresh` (value, default `1`), `floor` (value, default `0`).
  - `dowiener` (checkbox, default `0`): also output Wiener-cleaned coefficients.
  - `label_info` (label): notes power (not complex) subtraction + Welch averaging.
- `Run`:
  1. Resolve the data study → head model → `SurfaceFile`; load eigenmodes; build constrained gain (`in_bst_headmodel(...,1)`); good MEG/EEG channels (mirror `process_eigenmodes_transform`).
  2. Load the empty-room recording; compute the **common** good-channel set; build `A = bst_eigenmodes_transform(Gain(iCommon,:), Phi)`.
  3. `θ_data = A·d_data(iCommon,:)`, `θ_noise = A·d_noise(iCommon,:)`.
  4. Welch PSD of each (Brainstorm's PSD machinery, e.g. `bst_psd`), identical settings, density units, on the common frequency grid → `P_data`, `N`.
  5. `Out = bst_eigenmodes_noisefloor(P_data, N, 'Alpha',alpha, 'Floor',floor, 'SnrThresh',snrthresh)`.
  6. Save `timefreq_eigensnr` (`Out.SNR`, rows = modes labeled by λₖ, freqs) and `timefreq_eigencleanpsd` (`Out.CleanPSD`); report `K*(f)` summary via `bst_report`.
  7. If `dowiener`: apply `Out.Gain` to the full-block FFT of each `θ_data` row, iFFT (Hermitian) → cleaned `θ`; save `matrix_eigentransform_clean`.

## Verification (Milestone 2)

1. **Pure unit test** `dev/tests/test_eigenmodes_noisefloor_pure.m` (DB-free): synthetic `S`, `N`; `Pdata=S+N`; assert `CleanPSD≈S` (Alpha=1,Floor=0), `SNR=(S+N)./N`, `Gain∈[0,1]` and `=0` where `Pdata≤N`, `Floor` clamp respected for Alpha>1, and `Kstar` correct on an SNR ramp.
2. **checkcode** clean on both new files (Brainstorm idioms excepted).
3. **Optional manual OMEGA smoke** (reuse the `EigenSmoke` protocol / sub-0002): run denoise with the empty-room recording; confirm the 60 Hz line shows **low** SNR (it is present in both data and noise floor → correctly *not* flagged as signal), the alpha band shows SNR > 1 at low modes, and `K*(f)` falls with frequency. Save a figure.

## Non-goals (YAGNI — deferred)

- Per-mode (non-frequency-resolved) noise floor.
- The analytic covariance-propagation route (`Σ_n(f)` through `A`) — empirical only here.
- The broader heat/lowpass/bandpass filter library on coefficients (separate milestone; `bst_eigenmodes_filter` already has the transfer functions for vertex data).
- Fully reconstructed denoised **source-map time series** beyond the optional Wiener coefficient output.
- Over-subtraction/musical-noise refinements beyond α and a spectral floor.
- Vector / connection-Laplacian / Hodge work.

## Risks

- **Channel mismatch** data vs empty room → handled by intersecting good channels and rebuilding `A` on the common set.
- **PSD comparability** (different durations/settings) → enforce identical Welch settings and density (power/Hz) units for both.
- **Estimate variance / musical noise** from subtraction → mitigated by Welch averaging, the spectral floor β, and α; documented as a known limitation.
- **Empty-room availability** → the process errors clearly if no empty-room file is provided; the analytic covariance route remains a future fallback.
- **Branch base:** depends on unmerged M1 (`bst_eigenmodes_transform`); M2 branches off M1 and should merge after (or together with) it.
