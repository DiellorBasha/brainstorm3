# Eigenmode Spatial Transform + Spectrum — Milestone 1 Design

**Date:** 2026-05-27
**Branch:** `feature/eigenmode-transform-spectrum`
**Status:** Approved design, ready for implementation plan

## Background

The repo already has a scalar LBO eigenmode pipeline: `tess_eigenmodes`
(compute Φ), `in_/out_tess_eigenmodes` (storage), `bst_eigenmodes_project` /
`bst_eigenmodes_filter` (pure spatial-spectral math), and
`bst_inverse_eigenmodes` + `process_eigenmodes_inverse` (a **regularized**
eigenmode-space inverse: lead-field compression, noise whitening, a `λ^-α`
source prior, and MNE/dSPM/sLORETA normalization, all baked into one kernel).

The intended scientific architecture is **transform-first**, mirroring how
temporal signal analysis works in Brainstorm:

| Temporal (existing) | Spatial (this work) |
|---|---|
| Raw time series `d(t)` | Raw eigenmode coefficients `θₖ(t)` |
| FFT → spectrum (no filtering) | sensor→eigenmode transform, then FFT → joint (mode, freq) spectrum |
| *Look* at the spectrum, find the bands | *Look* at where power concentrates over spatial scale & time |
| *Then* choose a filter (bandpass/lowpass/notch) | *Then* choose a filter (heat/lowpass; noise-floor subtraction) |

Consequence: the **default product is the raw, unregularized, possibly-noisy
eigenmode coefficients**. All regularization (Tikhonov ≈ heat/lowpass filter,
dSPM/sLORETA, empty-room noise-floor denoising) becomes an *optional, separable*
second step — a future filter library applied to `θ`, not something baked into
the transform. The existing `bst_inverse_eigenmodes` therefore becomes *one
future entry* in that filter library, not the default path. It is **left
untouched** by this milestone.

## Goal (Milestone 1, scoped narrowly)

Implement the **pure, unregularized spatial eigenmode transform** of MEG/EEG
recordings, then FFT and visualize the resulting eigenmode spectrum using
Brainstorm's existing spectral tools:

1. Build the composite transform kernel `A = pinv(L̃)`, `L̃ = L·Φ` — no
   regularization, via SVD.
2. Apply it to recordings → eigenmode time series `Θ̂ = A·D` `[K × nTime]`,
   saved as a Brainstorm matrix file.
3. FFT of `Θ̂` via Brainstorm's existing `process_fft` (which accepts `matrix`
   inputs) → joint (eigenmode, frequency) spectrum.
4. Visualize that spectrum — the `(λₖ, ω)` plane — in Brainstorm's standard
   spectrum viewer.

## Non-goals (YAGNI — explicitly deferred)

- Empty-room noise-floor estimation / subtraction (denoising).
- The regularization filter library (Tikhonov/heat/lowpass; dSPM/sLORETA-as-filters).
- Any agreement / convergence testing against standard MNE/dSPM/sLORETA.
- Vector / connection-Laplacian / Hodge work.
- Refactoring `bst_inverse_eigenmodes` or `process_eigenmodes_inverse`.
- A bespoke `(λ, ω)` figure beyond what the standard spectrum viewer renders
  (an optional `imagesc` demo lives only in the test script).

## The math

The constrained (fixed-orientation) lead field is `L ∈ ℝ[nch × nVert]`
(`HeadModel.Gain` loaded with `ApplyOrient = 1`). The eigenmode matrix is
`Φ ∈ ℝ[nVert × K]`. Compress once:

```
L̃ = L · Φ ∈ ℝ[nch × K]      % column k = sensor topography of eigenmode k
```

The recordings relate to eigenmode coefficients by `D = L̃·Θ + noise`. With **no
regularization**, recover `Θ` by the Moore–Penrose pseudoinverse, computed via
SVD (`L̃ = U S Vᵀ`):

```
A = L̃⁺ = V S⁻¹ Uᵀ ∈ ℝ[K × nch]
Θ̂ = A · D ∈ ℝ[K × nTime]
```

Optional reconstruction back to vertices (raw, unfiltered): `Q = Φ·Θ̂`.

**Why SVD, not the normal equations.** The closed form flips with the
mode/channel ratio: the left-inverse `(L̃ᵀL̃)⁻¹L̃ᵀ` is valid only for `K ≤ nch`
(overdetermined), the right-inverse `L̃ᵀ(L̃L̃ᵀ)⁻¹` only for `K ≥ nch`
(underdetermined) — each form's normal matrix is singular in the other regime.
SVD-based `pinv` is correct in **both** regimes and avoids squaring the
condition number (forming `L̃ᵀL̃` would), which matters because `L̃` is
genuinely ill-conditioned: high-`λ` modes are nearly invisible to the sensors,
so their singular values are tiny.

**Regime / identifiability.** Default `K = min(nch, K_available)` makes the
transform well-determined (overdetermined LSQ). `K` is exposed as a parameter so
the underdetermined regime can be explored, with the understanding that modes
beyond `rank(L̃) ≈ nch` are min-norm fill-in, not data-determined.

**Noise is expected, not a bug.** Because `pinv` inverts the tiny singular
values of high-`k` modes, their coefficient time series are noise-amplified.
That is the point: the Milestone-1 spectrum will show clean banded structure at
low `k` and a noise floor at high `k`, and *that picture* is the diagnostic for
where to place a spatial filter later. MATLAB `pinv` zeroes only
*numerically*-zero singular values (tol = `max(size)·eps·max(s)`), so the
default is "unregularized but rank-safe" — honest noise amplification, no `Inf`.

## Architecture & data flow

```
Cortex surface (Φ on file)  ─┐
HeadModel.Gain (ApplyOrient=1)├─► bst_eigenmodes_transform ─► A=pinv(L̃) [K×nch]
                              ┘        (pure, SVD)            + Info (s, rank, cond)
                                                                   │
Recordings D [nch×nTime] ──────────────────────────────────────────┤ Θ̂ = A·D
                                                                   ▼
                          matrix_eigentransform  [K × nTime]  (Brainstorm matrix file)
                                                                   │
                                              process_fft (existing, accepts 'matrix')
                                                                   ▼
                          timefreq spectrum [K × nFreq]  →  standard spectrum viewer
                                                              (image = (λₖ, ω) plane)
```

## Component contracts

### `toolbox/math/bst_eigenmodes_transform.m` (new — pure)

Pairs symmetrically with `bst_eigenmodes_project` (which transforms a *known
cortical field* onto modes); this transforms *sensor data* onto modes through
the lead field. Pure linear algebra, no file I/O.

```matlab
[Kernel, Info] = bst_eigenmodes_transform(Gain, Phi, varargin)
```

- **Inputs**
  - `Gain` `[nch × nVert]`: constrained (fixed-orientation), channel-selected lead field.
  - `Phi`  `[nVert × K]`: eigenmode matrix (caller truncates to K).
- **Options**
  - `'Tol'` (default `[]` → MATLAB `pinv` default `max(size(L̃))·eps·max(s)`):
    singular-value floor (rank guard). Raising it is the crudest possible
    spatial filter — but Milestone 1 uses the default only.
- **Computation**
  - `L_tilde = Gain * Phi;`
  - `[U,S,V] = svd(L_tilde, 'econ'); s = diag(S);`
  - rank-safe inverse: `sinv = 1./s; sinv(s <= tol) = 0;`
  - `Kernel = V * diag(sinv) * U';`  % `[K × nch]`
- **Outputs**
  - `Kernel` `[K × nch]`: the transform `A = pinv(L̃)`.
  - `Info` struct: `.CompressedLF` (`L̃`), `.SingularValues`, `.Rank`,
    `.ConditionNumber` (`s(1)/s(rank)`), `.Tol`, `.nModes` (K).
- **Invariants** (the unit test asserts these): `Kernel*L_tilde ≈ Iₖ` for
  `K ≤ nch`; `L_tilde*Kernel ≈ I_nch` for `K ≥ nch`.

No whitening, prior, SNR, or normalization — by design.

### `toolbox/process/functions/process_eigenmodes_transform.m` (new)

The "first step" process the user runs. Reuses the file-resolution pattern of
`process_eigenmodes_inverse.m:138-217` (study → head model → surface →
`in_tess_eigenmodes` → vertex-count check → channel selection via
`good_channel`, MEG then EEG fallback), but **drops** the method/prior/SNR/noise
options.

- `GetDescription`
  - `Comment = 'Eigenmode transform (spatial FFT)'`, `Category = 'Custom'`,
    `SubGroup = 'Sources'`, `Index = 338` (immediately before the regularized
    inverse at 339).
  - `InputTypes = {'data','raw'}`, `OutputTypes = {'data','raw'}`,
    `nInputs = 1`, `nMinFiles = 1`.
  - Options:
    - `nmodes` (value, default `0` = auto = `min(nch, available)`).
    - `dorecon` (checkbox, default `0`): also save raw vertex reconstruction
      `Q = Φ·Θ̂` as a results file.
    - `label_info`: notes this is the unregularized transform; coefficients are
      raw (noisy at high modes by design); requires precomputed eigenmodes.
- `Run`
  - Resolve `HeadModelFile`, `SurfaceFile`, `Eigenmodes`, `GoodChannel`
    (identical to the inverse process).
  - `K = ` requested or `min(nGoodChan, Eigenmodes.nModes)`; `Phi =
    Eigenmodes.Vectors(:,1:K)`; `lambdas = Eigenmodes.Values(1:K)`.
  - Load constrained gain: `HM = in_bst_headmodel(HeadModelFile, 1)`; select
    good channels → `Gain` `[nGoodChan × nVert]`.
  - `[Kernel, Info] = bst_eigenmodes_transform(Gain, Phi);` report
    `Info.ConditionNumber`, `Info.Rank` via `bst_report('Info', ...)`.
  - **Imported data**: `Θ̂ = Kernel * D(iGood,:)` → save
    `matrix_eigentransform` (`db_template('matrixmat')`): `.Value = Θ̂`,
    `.Time = DataMat.Time`, `.Description = {'Mode k (lam=…)'}`,
    `.SurfaceFile`, `bst_history` entries. **This file is the FFT input.**
    If `dorecon`, also save `results_eigentransform` with
    `ImagingKernel = Φ·Kernel` `[nVert × nGoodChan]`, `nComponents = 1`,
    standard results metadata (mirror `process_eigenmodes_inverse.m:392-426`).
  - **Raw data**: cannot precompute coefficients (same constraint as the
    inverse process, `:264-271`); save the kernel `Φ·Kernel` as a kernel-only
    `results_eigentransform` and warn that coefficients require imported data.
  - File prefixes: `matrix_eigentransform`, `results_eigentransform`.

### FFT + visualization (no new code)

- Run the existing **`process_fft`** on the `matrix_eigentransform` file
  (`process_fft.InputTypes` already includes `'matrix'`,
  `process_fft.m:37`) → a Brainstorm timefreq/spectrum file `[K × nFreq]`.
- Open it in the **standard spectrum viewer**; in image mode the rows are the
  eigenmodes (ordered by `λₖ`, i.e. spatial frequency) and the columns are
  temporal frequency `ω` — directly the joint `(λ, ω)` plane.
- The test script may additionally `imagesc(freqs, lambdas, 10*log10(power))`
  as a quick standalone demo, but no production visualization code is added.

## Verification (Milestone 1 — minimal, matches repo norms)

Heavy/comparison testing is explicitly deferred. Two lightweight checks only:

1. **Pure unit test** `dev/tests/test_eigenmodes_transform_pure.m` (DB-free,
   synthetic, in the style of `test_eigenmodes_project_pure.m`):
   - Random `Gain` `[nch × nVert]`, M-orthonormal `Phi`, `K < nch`.
   - Noise-free `D = (Gain·Phi)·c0` → `Kernel·D` recovers `c0` to ~1e-9.
   - Left-inverse invariant `Kernel·L̃ ≈ Iₖ` (K < nch); right-inverse invariant
     `L̃·Kernel ≈ I_nch` for a second case with `K > nch`.
   - `Info.ConditionNumber` finite and ≥ 1; `Info.Rank == min(nch, K)`.
2. **checkcode / M-Lint** clean (Brainstorm idioms aside) on both new files.

Optional, manual (not required to close the milestone): an end-to-end smoke on
OMEGA sub-0002 reusing `test_omega_icosphere_sourcemap.m`'s setup — compute
eigenmodes, run the transform process on a short imported block, run
`process_fft`, eyeball the `(λ, ω)` image.

## Risks

- **`K > nch` confusion.** Users may request more modes than channels and read
  the (min-norm, not data-determined) high modes as signal. Mitigated by the
  `min(nch, …)` default and the `label_info` note; fully addressed later by the
  noise-floor step.
- **Plumbing duplication** with `process_eigenmodes_inverse` (~80 lines of
  file/channel resolution). Accepted for Milestone 1 to keep the change additive
  and the working inverse untouched; a shared resolver helper is a noted future
  cleanup, not part of this milestone.
- **Matrix-file FFT specifics.** `process_fft` on a `matrix` file is the assumed
  path; if its options need a non-obvious setting for mode-row spectra, that is
  resolved at implementation time (it changes only how the FFT step is invoked,
  not the transform contract).
