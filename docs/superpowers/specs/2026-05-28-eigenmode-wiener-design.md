# Eigenmode Wiener Filter (noise-floor spectral gain) — Design

**Date:** 2026-05-28
**Branch:** `feature/eigenmode-wiener` (off `development`, which has M1 transform + M2 noise floor + coefficient filter library + M3 wavelet tensor + dispersion)
**Status:** Approved design, ready for implementation plan

## Background

The transform-first pipeline produces raw eigenmode coefficients `θₖ(t)`, then applies optional filters. The spatial filter library (`bst_eigenmodes_filter_gain` + `process_eigenmodes_coeffsfilter`) applies a **frequency-flat** per-mode gain `h(λₖ)` by plain time-series multiplication. This increment adds the first **spectral** (frequency-dependent) filter: a Wiener gain `G(k,f)` estimated from an empty-room noise floor.

M2 already computes the Wiener gain — `bst_eigenmodes_noisefloor` returns `Gain(k,f) = max(Pdata−Nnoise,0)/Pdata`, the per-`(k,f)` gain in `[0,1]` — but the M2 process (`process_eigenmodes_denoise`) computes it and **discards** it, saving only SNR and CleanPSD. So this increment is about **applying** that gain to data, not deriving it.

Because the gain varies with frequency, it cannot be applied by multiplying the time series (as the spatial filters do); it must be applied in the frequency domain. Brainstorm's `bst_bandpass_fft` gives the exact reusable idiom: mirror the signal → `fft` → multiply by a zero-phase magnitude response → `real(ifft(…))` → trim.

See [[eigenmode-meg-project]], [[eigenmode-transform-first]].

## Goal (this increment)

A self-contained 2-input process (Files A = data, Files B = empty-room) that, from the **same** coefficients and PSDs, computes the Wiener gain `G(k,f)` and applies it to the data's eigenmode-coefficient time series, producing a denoised coefficient time series (and, optionally, the gain spectrum and a vertex reconstruction).

## Non-goals (YAGNI — deferred)

- Time-varying / wavelet-tensor Wiener (the gain here is stationary, from Welch-averaged PSDs).
- Consuming a separately-saved gain product (the composable design was considered and rejected in favor of correctness-by-construction).
- New filter types in the spatial gain library (`bst_eigenmodes_filter_gain` is untouched).
- Vector / connection-Laplacian work.

## Architecture & data flow

```
process_eigenmodes_wiener  (A = data, B = empty-room)
  └─ S = process_eigenmodes_denoise('GetCoeffsAndPSD', sProcess, sInputsA, sInputsB, nModes, WinLen)
        S = { Coeffs Θ(t) [K×nTime], Pdata [K×nFreq], Nnoise [K×nFreq], Freqs [1×nFreq],
              lambdas [K×1], Phi [nVert×K], SurfaceFile, K, Time, sfreq, Info }   (or [] after bst_report error)
  └─ G = bst_eigenmodes_noisefloor(S.Pdata, S.Nnoise, 'Alpha',alpha, 'GainFloor',gainfloor).Gain   [K×nFreq]
  └─ Θ̂ = bst_eigenmodes_wiener(S.Coeffs, S.sfreq, G, S.Freqs, 'Mirror',domirror)   [K×nTime], real
  └─ save matrix_eigenwiener (Θ̂)   [+ optional timefreq_eigenwienergain, + optional results_eigenwiener]
```

The Wiener gain is computed from the same coefficients/PSDs it filters — no cross-file consistency footgun.

## Component contracts

### (a) `toolbox/math/bst_eigenmodes_noisefloor.m` (extend — behavior-preserving)

Make the `Gain` field honor over-subtraction and a floor:

```matlab
Out.Gain = max( max(Pdata - Alpha.*Nnoise, 0) ./ max(Pdata, eps), GainFloor );
```

- `Alpha` already exists (default `1`, used by CleanPSD). New option `GainFloor` (default `0`).
- Validation: error if `GainFloor` not in `[0,1]`; error if `Alpha < 0`.
- **Back-compatible:** defaults `Alpha=1, GainFloor=0` reproduce the current `max(Pdata−Nnoise,0)./max(Pdata,eps)` exactly. `SNR`, `CleanPSD`, `Kstar` are unchanged.
- Update the header to document `GainFloor` and the new `Gain` formula.
- **Guard:** the existing `dev/tests/test_eigenmodes_noisefloor_pure.m` must still pass unchanged.

### (b) `process_eigenmodes_denoise('GetCoeffsAndPSD', …)` (extract — shared helper)

Extract M2's `Run` setup into a subfunction exposed via the `eval(macro_method)` dispatch:

```matlab
S = process_eigenmodes_denoise('GetCoeffsAndPSD', sProcess, sInputsA, sInputsB, nModesOpt, WinLen)
```

Moves out of `Run`, verbatim in behavior: head-model / surface / eigenmode lookup and surface-type check; common good-channel intersection (by name) between data and empty-room; kernel `A = pinv(L·Φ) = bst_eigenmodes_transform(Gain(iCommonA,:), Phi)`; coefficients `thD = A·dataF`, `thN = A·noiseF`; Welch PSDs `bst_psd(…)` of both with the common-grid check. Returns struct `S` with fields `Coeffs` (=`thD` `[K×nTime]`), `Pdata`, `Nnoise` `[K×nFreq]`, `Freqs`, `lambdas` `[K×1]`, `Phi` `[nVert×K]`, `SurfaceFile`, `K`, `Time` (`DA.Time`), `sfreq` (`1/(DA.Time(2)−DA.Time(1))`), and `Info` (transform info). On any failure it calls `bst_report('Error', sProcess, sInputsA, …)` and returns `[]`.

`process_eigenmodes_denoise('Run', …)` is refactored to call the helper, then do its own PSD-domain combine/save — **behavior-preserving** (same output files as today). The existing `test_process_eigenmodes_denoise_options.m` must still pass; the `Run` behavior is re-verified by the live smoke (no automated process-level test exists for M2's `Run`).

### (c) `toolbox/math/bst_eigenmodes_wiener.m` (new — pure)

```matlab
Filtered = bst_eigenmodes_wiener(Coeffs, sfreq, Gain, GainFreqs, varargin)
```

- **Inputs:** `Coeffs` `[K×nTime]` real coefficient time series; `sfreq` (Hz); `Gain` `[K×nGainFreq]` in `[0,1]`; `GainFreqs` `[1×nGainFreq]` (Hz, the Welch grid, ascending).
- **Options:** `'Mirror'` (default `true`) — reflect the signal at both ends before the FFT and trim after, to suppress circular-convolution edge artifacts (as in `bst_bandpass_fft`).
- **Computation:** with FFT length `N` (after optional mirroring), build the two-sided frequency vector and fold it to `[0, Nyquist]`; for each mode `k`, interpolate `Gain(k,·)` from `GainFreqs` onto the folded grid (`interp1` linear, endpoints held flat beyond `GainFreqs`, then clamp to `[0,1]`). Because the gain is interpolated on the **folded** (absolute) frequency, the resulting response `H_k` is Hermitian-symmetric, so `Filtered = real(ifft(fft(Coeffs,[],2) .* H, [], 2))` is real and **zero-phase**. Trim the mirror.
- **Output:** `Filtered` `[K×nTime]`, real, same size as `Coeffs`.
- **Errors:** size mismatch between `size(Gain,1)` and `size(Coeffs,1)`; `numel(GainFreqs) ~= size(Gain,2)`.
- Pure — no file I/O.

### (d) `toolbox/process/functions/process_eigenmodes_wiener.m` (new)

- `GetDescription`: Comment `'Eigenmode Wiener filter'`, Category `'Custom'`, SubGroup `'Sources'`, Index `336.95` (after coeffsfilter 336.9; verified free), InputTypes `{'data','raw'}`, OutputTypes `{'matrix'}`, `nInputs=2` (A=data, B=empty-room), `nMinFiles=1`. Options:
  - `nmodes` (value, `0` = auto = `min(nCh, nModes)`).
  - `noisewin` (value, `2 s` — Welch window).
  - `alpha` (value, `1` — over-subtraction, `≥1`).
  - `gainfloor` (value, `0` — gain floor `Gmin∈[0,1]`).
  - `domirror` (checkbox, default `1`).
  - `dorecon` (checkbox, default `0`).
  - `savegain` (checkbox, default `0`).
  - `label_info` (label).
- `FormatComment`: reflect `alpha`/`gainfloor` (e.g. `Eigenmode Wiener filter (a=1.0, Gmin=0.10)`).
- `Run(sProcess, sInputsA, sInputsB)`:
  1. Error if `sInputsB` empty.
  2. `S = process_eigenmodes_denoise('GetCoeffsAndPSD', sProcess, sInputsA, sInputsB, nmodes, noisewin)`; return if `[]`.
  3. `NF = bst_eigenmodes_noisefloor(S.Pdata, S.Nnoise, 'Alpha',alpha, 'GainFloor',gainfloor)`; `G = NF.Gain`.
  4. `Filtered = bst_eigenmodes_wiener(S.Coeffs, S.sfreq, G, S.Freqs, 'Mirror',domirror)`.
  5. Save `matrix_eigenwiener` (`db_template('matrixmat')`: `Value=Filtered`, `Time=S.Time`, `Description` = per-mode `Mode k (lam=…)` row names, `SurfaceFile`, `nAvg=1`, history). Append to `OutputFiles`.
  6. If `savegain`: save `G` as `timefreq_eigenwienergain` (`[K×1×nFreq]`, `Freqs=S.Freqs`, `Measure='power'` — matches M2's unitless-ratio convention so it renders in the TF viewer, mode row names, `SurfaceFile`, history) — `db_add_data` sidecar, **not** in `OutputFiles`.
  7. If `dorecon`: `Q = S.Phi*Filtered` `[nVert×nTime]`; save `results_eigenwiener` (`ImageGridAmp=Q`, `ImagingKernel=[]`, `nComponents=1`, `Time=S.Time`, `SurfaceFile`, `HeadModelType='surface'`, `ColormapType='source'`, history) — `db_add_data` sidecar, **not** in `OutputFiles`.

`OutputFiles` stays homogeneous (`matrix` only); the gain spectrum and the vertex reconstruction are DB sidecars (mirrors the coeffsfilter recon-sidecar decision).

## Verification (this increment)

1. **`dev/tests/test_eigenmodes_wiener_pure.m`** (DB-free):
   - **Identity:** `Gain ≡ 1` ⇒ output ≈ input (tight tol).
   - **Null:** `Gain ≡ 0` ⇒ output ≈ 0.
   - **Selectivity:** two-tone signal (tones `fA`, `fB`); gain `1` around `fA`, `0` around `fB` ⇒ `fA` preserved, `fB` removed (compare per-tone power before/after).
   - **Zero-phase:** a single cosine in a pass band stays in phase (no time shift).
   - **Shape/realness:** output is real and `size == size(Coeffs)`.
   - **Errors:** row mismatch (`size(Gain,1)≠size(Coeffs,1)`) and `numel(GainFreqs)≠size(Gain,2)` both error.
2. **`dev/tests/test_eigenmodes_noisefloor_pure.m`** (extend): `Alpha=1,GainFloor=0` reproduces the prior `Gain`; `GainFloor=g` ⇒ `min(Gain(:)) ≥ g`; `Alpha>1` ⇒ elementwise-smaller (or equal) gain; `GainFloor` outside `[0,1]` errors; existing assertions unchanged.
3. **`dev/tests/test_process_eigenmodes_wiener_options.m`** (DB-free): Category/SubGroup/Index/InputTypes/OutputTypes/`nInputs` + all options present with documented defaults.
4. **`checkcode`** clean on all touched files (Brainstorm idioms excepted).
5. **Optional live smoke** (throwaway protocol, imported data + empty-room with a surface head model + eigenmodes): run with defaults and with `alpha>1`/`gainfloor>0`; confirm the filtered coefficients have reduced low-SNR/high-frequency power vs. input, `savegain` writes a `[0,1]` gain spectrum, and `dorecon` yields a finite vertex map. Tear down the protocol. Report `DONE_WITH_CONCERNS` if substrate setup proves intractable — the pure tests carry the math.

## Risks

- **Circular-convolution wrap-around** from full-signal FFT → mitigated by edge-mirroring (default on; same as `bst_bandpass_fft`).
- **M2 `Run` refactor has no automated process-level guard** (only the pure noise-floor test + the options test are automated) → mitigated by extracting the setup verbatim (helper returns exactly M2's prior intermediates) and re-running the live smoke for M2's denoise path.
- **Coarse gain grid:** the Welch PSD `df` is coarser than the coefficient FFT `df`; interpolation smooths the response. Standard and acceptable for a stationary Wiener gain.
- **Zero-phase / magnitude-only:** the Wiener gain is a real magnitude, so signal phase is preserved by construction (documented). No attempt to estimate or alter phase.
- **Short recordings:** if `nTime` is too small for mirroring, the apply function falls back to no mirror (documented) rather than erroring.
