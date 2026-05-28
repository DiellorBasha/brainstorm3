# Eigenmode Coefficient Filter Library — M2/M3 increment Design

**Date:** 2026-05-28
**Branch:** `feature/eigenmode-coeffsfilter` (off `development`, which has M1 transform + M2 noise floor + M3 wavelet tensor + dispersion)
**Status:** Approved design, ready for implementation plan

## Background

The transform-first design (M1/M2) is: produce raw eigenmode coefficients, inspect the spectrum, then apply a filter. This increment builds the **"then filter" half** as a per-mode spatial-spectral filter on the coefficient time series. Filtering coefficients is trivial — multiply `θₖ(t)` by a per-mode gain `h(λₖ)` (no project/reconstruct needed, since the data is already in the eigenmode basis).

The existing `toolbox/math/bst_eigenmodes_filter.m` already implements these transfer functions, but for **vertex** fields (project → multiply `h` → reconstruct), with types `lowpass/highpass/bandpass/heat/inverse_heat/custom`. To avoid two copies of the formulas, the transfer-function computation is **extracted into a shared pure function** that both the vertex filter and the new coefficient filter call (DRY refactor, confirmed).

The Wiener-from-noise-floor filter (M2's per-`(k,f)` gain, frequency-domain) is a separate, later increment.

See [[eigenmode-meg-project]], [[eigenmode-transform-first]].

## Goal (this increment)

1. Extract `bst_eigenmodes_filter_gain(lambdas, FilterType, …)` → `[K×1]` gain `h(λₖ)` — the single source of the transfer functions — and add a `tikhonov` type `h = 1/(1+β·λ)`.
2. Refactor `bst_eigenmodes_filter` to compute its gain via that function (behavior-preserving; guarded by its existing test).
3. Add `process_eigenmodes_coeffsfilter` (Index 336.9): apply `h(λₖ)` to an eigenmode-coefficient `matrix` file → filtered coefficients (+ optional vertex reconstruction).

## Non-goals (YAGNI — deferred)

- Wiener-from-noise-floor filter (per-`(k,f)`, frequency-domain, needs the M2 noise floor) — separate increment.
- Time-varying / wavelet-tensor filtering.
- New filter types beyond the existing set + `tikhonov`.
- Vector / connection-Laplacian work.

## Component contracts

### `toolbox/math/bst_eigenmodes_filter_gain.m` (new — pure)

```matlab
h = bst_eigenmodes_filter_gain(lambdas, FilterType, varargin)
```

- **Input:** `lambdas` `[K×1]` eigenvalues (defines `K = numel(lambdas)`); `FilterType` one of `lowpass`/`highpass`/`bandpass`/`heat`/`inverse_heat`/`tikhonov`/`custom`.
- **Options** (name-value, same defaults as today): `CutoffMode` (50), `ModeRange` ([20 80]), `DiffusionTime` (0.01), `MaxGain` (10), `TransferFn` ([]), and new `RegBeta` (1, for `tikhonov`).
- **Gains** (verbatim from the current `bst_eigenmodes_filter` transfer block, plus tikhonov):
  - `lowpass`: `h(1:min(CutoffMode,K))=1`, else 0.
  - `highpass`: `h(max(1,min(CutoffMode,K)):end)=1`.
  - `bandpass`: `h(max(1,ModeRange(1)):min(K,ModeRange(2)))=1`.
  - `heat`: `exp(-DiffusionTime*lambdas)` (error if `DiffusionTime<=0`).
  - `inverse_heat`: `min(exp(DiffusionTime*lambdas), MaxGain)` (error if `DiffusionTime<=0`).
  - `tikhonov`: `1 ./ (1 + RegBeta*lambdas)` (error if `RegBeta<0`).
  - `custom`: `TransferFn(lambdas)` (error if not a handle or wrong length); `h=h(:)`.
  - otherwise: error listing the valid types.
- **Output:** `h` `[K×1]`. Pure, no I/O.
- **Unit-tested:** lowpass/bandpass masks exact; `heat` → identity as `t→0`, strong suppression of high-λ as `t` large; `tikhonov` decreasing in λ and `=1` at `λ=0`; `inverse_heat` clamped at `MaxGain`; `custom` passthrough; unknown type errors.

### `toolbox/math/bst_eigenmodes_filter.m` (refactor — behavior-preserving)

Replace the inline option-parse + `BUILD TRANSFER FUNCTION` block with a single call:
```matlab
h = bst_eigenmodes_filter_gain(lambdas, FilterType, varargin{:});
```
Keep the `Phi`/`lambdas` extraction, the Data/MassMatrix validation, and the apply step `Filtered = Phi * (h .* (Phi' * (MassMatrix * Data)))`. Update the header's type list to include `tikhonov`. The existing `dev/tests/test_eigenmodes_filter_pure.m` must still pass unchanged (it guards the refactor).

### `toolbox/process/functions/process_eigenmodes_coeffsfilter.m` (new)

- `GetDescription`: Comment `'Eigenmode coefficient filter'`, Category `'Custom'`, SubGroup `'Sources'`, Index `336.9` (verified free), InputTypes `{'matrix'}`, OutputTypes `{'matrix'}`, nInputs 1, nMinFiles 1. Options:
  - `filtertype` (radio_linelabel: `lowpass`/`highpass`/`bandpass`/`heat`/`inverse_heat`/`tikhonov`, default `heat`).
  - `cutoffmode` (value, 50), `moderange` (value pair, [20 80] — or two value options `kmin`/`kmax`), `diffusiontime` (value, 0.01 s), `regbeta` (value, 1).
  - `dorecon` (checkbox, default 0): also reconstruct vertex sources `Q = Φ·θ_filt`.
  - `label_info` (label).
- `Run`:
  1. `M = in_bst_matrix(file)`; `Coeffs = M.Value` `[K×nTime]`, `Time`, `Description`, `SurfaceFile`.
  2. Error if no `SurfaceFile`; `Em = in_tess_eigenmodes(SurfaceFile)` (error if absent / fewer than K eigenvalues); `lambdas = Em.Values(1:K)`.
  3. `h = bst_eigenmodes_filter_gain(lambdas, filtertype, <opts from GUI>)`; `Coeffs_filt = h .* Coeffs`.
  4. Save filtered coefficients as `matrix_eigenfilt` (`db_template('matrixmat')`: `Value=Coeffs_filt`, `Time`, `Description`, `SurfaceFile`, Comment naming the filter+params, history).
  5. If `dorecon`: `Phi = Em.Vectors(:,1:K)`; `Q = Phi*Coeffs_filt` `[nVert×nTime]`; save `results_eigenfilt` (`ImageGridAmp=Q`, `ImagingKernel=[]`, `nComponents=1`, `Time`, `SurfaceFile`, `HeadModelType='surface'`, history).

## Verification (this increment)

1. **Pure unit test** `dev/tests/test_eigenmodes_filter_gain_pure.m` (DB-free): assert each gain type (masks, heat limits, tikhonov monotonicity + `h(λ=0)=1`, inverse_heat clamp, custom passthrough, unknown-type error).
2. **Refactor guard:** the existing `test_eigenmodes_filter_pure` (vertex filter) still passes unchanged.
3. **Process options test** `dev/tests/test_process_eigenmodes_coeffsfilter_options.m` (DB-free): SubGroup/Index/InputTypes/OutputTypes + options present.
4. **Optional live smoke:** a throwaway protocol with a sphere surface + eigenmodes + a synthetic coefficient matrix → run the process with `heat` → confirm filtered coefficients have suppressed high modes relative to low (and `dorecon` produces a finite vertex map). Tear down the protocol. (As in the dispersion increment; report DONE_WITH_CONCERNS if substrate setup proves intractable — the pure tests cover the math.)

## Risks

- **Refactor behavior change:** the extracted gain must reproduce the current formulas byte-for-byte; `test_eigenmodes_filter_pure` is the guard (its lowpass + heat-limit assertions must still pass). The option defaults must match exactly (`CutoffMode=50`, `ModeRange=[20 80]`, `DiffusionTime=0.01`, `MaxGain=10`).
- **Option threading:** `bst_eigenmodes_filter` now passes `varargin{:}` straight to the gain function — the option names/defaults must stay identical so existing callers (`process_eigenmodes_filter`) are unaffected.
- **λ alignment:** the coefficient matrix's `K` rows must align with `Em.Values(1:K)` (guarded by the eigenvalue-count check); modes with `λ≤0` are passed through by heat/tikhonov harmlessly (`exp(0)=1`, `1/(1+0)=1`).
- **`dorecon` size:** `Q=[nVert×nTime]` can be large for long recordings; it is off by default and the primary output is the (small) filtered coefficient matrix.
