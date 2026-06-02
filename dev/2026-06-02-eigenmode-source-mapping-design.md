# Eigenmode-based source mapping (GBF) — forward/inverse split

**Date:** 2026-06-02
**Author:** Diellor Basha (design captured with Claude)
**Status:** Design — pending review before implementation plan
**Repo:** `research/code/brainstorm3` (branch `development`)

---

## 1. Motivation

We want to replicate the Geometric Basis Functions (GBF) method
(`research/code/GBFs`) for MEG/EEG source mapping inside Brainstorm. GBF solves
the inverse problem in a compact Laplace–Beltrami eigenmode (LBO) basis rather
than in vertex space:

1. Compute `K` LBO eigenmodes `(Φ, Λ)` on the cortical surface.
2. **Build the eigenmode leadfield** `L̃ = L·Φ` (`nChannels × K`) by composing the
   standard vertex-space leadfield `L` with the eigenmode matrix `Φ`.
3. Define a spectral prior covariance on the mode coefficients.
4. Solve the MAP estimate for the coefficients `θ̂`.
5. Reconstruct vertex sources `x̂ = Φ·θ̂`.

### What is wrong with the current Brainstorm code

The current integration has three different "sensor → mode coefficient" paths,
and the two that are actually wired into the viewer / spectrum analyses are *not*
GBF:

- **`bst_inverse_eigenmodes.m`** fuses the forward composition (`L̃ = L·Φ`) and
  the regularized inverse into one function. Mixing the forward and inverse makes
  the prior/regularization step awkward and non-idiomatic.
- **Harmonic path** (`bst_eigenmodes_harmonic` = `pinv(iW·L·Φ)·iW`) is
  unregularized — GBF with `β→0` and no prior, the ill-conditioned regime GBF was
  designed to avoid. This is what the time-series viewer currently uses.
- **Projection path** (`bst_eigenmodes_project` / `_modekernel` = `Φ'·M·K_std`)
  projects an *already reconstructed* standard inverse onto the eigenbasis. The
  hard, ill-posed vertex inverse has already happened; the compression buys no
  conditioning benefit. This is a legitimate **eigenspectrum analysis** tool but
  is not source mapping.

### Design principle

Brainstorm already separates the **forward** model (`headmodel_*.mat`) from the
**inverse** (`results_*.mat`). GBF's `L̃ = L·Φ` is *strictly a forward solution*.
We therefore split the method along Brainstorm's existing seam:

```
FORWARD  (toolbox/forward/)   base leadfield  →  composed eigenmode leadfield  L̃ = L·Φ
INVERSE  (toolbox/inverse/)   composed L̃  +  data cleaning  +  spectral prior R  →  M̃
DISPLAY                        coefficients θ = M̃·d   and   cortex sources Φ·M̃
```

This makes the inverse a clean, standard-conventions Brainstorm inverse whose
spatial prior is the eigenvalue spectrum, and lets the forward be reused with any
physics method (overlapping spheres, OpenMEEG BEM, DUNEuro FEM).

---

## 2. Scope

### In scope

- A standalone forward composer that turns any base head model into a composed
  eigenmode head model (`Gain = L·Φ`).
- A dedicated mode-space inverse that consumes the composed head model, applies
  the full Brainstorm data-cleaning pipeline (bad channels, SSP/ICA projectors,
  noise whitening), and solves the MAP estimate with the eigenvalue spectral prior
  as the source covariance `R`.
- Standard output nodes: a coefficient matrix and a kernel-only cortex results
  node.
- A three-level validation harness (resolution matrix, ground-truth simulation,
  real data) comparing against Brainstorm defaults.

### Out of scope (unchanged)

- **Eigenspectrum methods** — `bst_eigenmodes_project/_filter/_spectrum/_wiener/
  _dispersion/_wavelet/_modekernel` and their processes. They analyze fields
  already mapped to the surface; a separate concern that stays as-is.

### Retired (clean rebuild)

- `bst_eigenmodes_harmonic.m`, `bst_eigenmodes_transform.m`, the bespoke Harmonic
  results node, and `view_eigenmodes_timeseries.m`. The unregularized "harmonic"
  behavior survives as a **flat-`R` option** of the new inverse. Coefficients become
  a standard matrix node (viewable normally, consumable by the eigenspectrum tools).

---

## 3. Architecture & data flow

```
   [ base head model ]        headmodel_*.mat  (Gain = K, full BEM/OS/FEM leadfield)
   [ surface eigenmodes ]     tess_*.mat       (Φ, λ_k  via tess_eigenmodes)
              │
              ▼   STAGE 1 — FORWARD  (toolbox/forward/)
   bst_eigenmode_leadfield(baseHeadModel, Eigenmodes)
        L  = constrained(K)             [nChan × nVert]   (apply surface-normal orientation)
        L̃ = L · Φ                       [nChan × K]
              │
              ▼   writes a new node
   headmodel_eigenmode_*.mat   Gain = L̃,  stores λ_k + SurfaceFile + nModes + isEigenmode
              │
              ▼   STAGE 2 — INVERSE  (toolbox/inverse/)
   bst_inverse_eigenmodes(composedHM, NoiseCov, Projector, ChannelFlag, ...)
        clean : drop bad channels, apply SSP projectors, whiten by C (bst_whitener)
        prior : R = diag(σ²(λ_k))        ← 2026 log prior (replaces depth weighting)
        solve : M̃ = R L̃ᵀ (L̃ R L̃ᵀ + λ C)⁻¹   [K × nChan]   (MNE / dSPM / sLORETA)
              │
              ▼   STAGE 3 — OUTPUTS
   ┌─ coefficients  matrix_*.mat   θ = M̃·d   [K × nTime]   → feeds eigenspectrum tools
   └─ cortex sources results_*.mat ImagingKernel = Φ·M̃   [nVert × nChan]  → normal 3D display
```

### Files

| File | Status | Role |
|---|---|---|
| `toolbox/forward/bst_eigenmode_leadfield.m` | **new** | Engine: base head model + `Φ` → composed head model struct (`Gain = L·Φ`). |
| `toolbox/process/functions/process_eigenmode_leadfield.m` | **new** | Process wrapper → writes the `headmodel_eigenmode_*.mat` node. |
| `toolbox/inverse/bst_inverse_eigenmodes.m` | **rewritten** | Consumes the composed leadfield (never builds `L·Φ`); data cleaning + spectral-prior `R` + regularized solve; returns `M̃`. |
| `toolbox/process/functions/process_eigenmodes_inverse.m` | **rewritten** | Runs the inverse; emits coefficient + cortex nodes. |
| `bst_eigenmodes_harmonic.m`, `bst_eigenmodes_transform.m`, `view_eigenmodes_timeseries.m` | **retired** | Replaced by flat-`R` option + standard nodes. |

### Composed head-model representation

Saved as a real `headmodel_*.mat` node so "forward feeds inverse" holds literally:

- `Gain = L̃` `[nChan × K]`
- `HeadModelType = 'surface'` plus marker field `isEigenmode = 1`
- `nModes = K`, `Eigenvalues = λ_k`
- `SurfaceFile` = original cortex (source of `Φ` for reconstruction)
- `GridLoc`/`GridOrient` left empty (not needed; the cortex results node carries the
  real surface). **`Φ` is not duplicated into the head model** — reloaded from the
  surface (`in_tess_eigenmodes`) only at the `Φ·M̃` reconstruction step.

---

## 4. The mechanics

### 4.1 Spectral prior `R` from eigenvalues

The spectral prior **is** the source covariance `R` in the standard inverse
`J = R·Lᵀ·(L·R·Lᵀ + λ·C)⁻¹·d`. In mode space it is diagonal and indexed by mode
`k` (not by vertex), so it **replaces depth weighting** — no `GridLoc` needed.

Default = **2026 log prior**. GBF's penalty is `Λ·diag(−1/log λ_k)`, i.e.
`Σ⁻¹ ∝ −1/log λ_k`, so `R = diag(−log λ_k)` up to scale.

**Eigenvalue scale — match GBF exactly (millimetre mesh).** `−log λ_k` is positive
only for `λ_k < 1`. GBF achieves this by building the mesh in **millimetres**.
Brainstorm surfaces are in **metres**, but the cotan stiffness is scale-invariant and
only the mass matrix scales with area, so converting is exact and leaves the
eigenvectors unchanged:

```
λ_mm = λ_m · 1e−6        (m⁻² → mm⁻²; coordinate scale 1e3 ⇒ eigenvalue scale 1e−6)
R_k  ∝ −log(λ_mm,k)      (DC mode λ_1≈0 swapped to λ_2; clamp λ_mm ∈ (0,1))
```

(`R` is the source *covariance*. GBF writes the penalty as the *precision*
`Σ⁻¹ = −1/log λ`, so the covariance is `R = Σ = −log λ` — positive and **decreasing**
in `λ`: smooth low-`λ` modes get more prior variance, fine modes less.)

The millimetre scale adds a large constant offset (`−log(1e−6) ≈ 13.8`) to every
`R_k`, so the high modes are **gently** rolled off rather than annihilated — the edge
mode keeps ≈0.25–0.8× the DC variance. `R` is normalized so `max(R)=1`; absolute
scale is absorbed by the global regularizer. This is intentionally **scale-dependent**
(the millimetre scale is physical) — that is the GBF design.

> **Why not a dimensionless `λ_ref` normalization?** An earlier draft normalized by
> `λ_ref = λ_{K+1}` to be unit-independent. Validation (Level 1 below) showed this
> drives the highest *retained* mode's variance to ~1e−7 — it deletes the modes that
> localize, giving ~2× worse localization than dSPM. GBF's raw millimetre scaling
> avoids this. **Lesson: reproduce GBF's scaling literally; do not re-derive it.**

**Prior options:** `log` (default) · `flat` (`R = I`) · `power` (`λ^{−α}`, legacy).
The retired **harmonic** behavior = `flat` prior **with regularization disabled**
(`λ→0`, i.e. `pinv(L̃)`), exposed as an explicit "unregularized" switch rather than
a prior choice.

### 4.2 Regularization

Brainstorm SNR convention: `λ = trace(L̃ R L̃ᵀ) / (nChan · SNR²)`, default `SNR = 3`.
This is GBF's `Λ`, data-scaled so the SNR knob means the same thing regardless of
leadfield scale.

### 4.3 Data cleaning (reused, matching the standard inverse)

1. **Bad channels** — `good_channel` / `ChannelFlag` selection.
2. **SSP/ICA projectors** — `ChannelMat.Projector`, applied as the 2018 inverse does.
3. **Whitening** — `bst_whitener(NoiseCovMat, …)`, the same call the 2018 inverse
   makes (`NoiseMethod` reg/shrink/diag/none); whiten `L̃` and data.
   Order matches 2018 (projector → whiten).

### 4.4 Outputs (standard nodes, no bespoke viewer)

- **Coefficients** — `matrix_*.mat`, `Value = θ = M̃·d` `[K×nTime]`, row labels carry
  `λ_k`, `SurfaceFile` set → consumable by the untouched eigenspectrum tools.
- **Cortex sources** — `results_*.mat`, **kernel-only**, `ImagingKernel = Φ·M̃`
  `[nVert×nChan]`, `nComponents=1`, `Function='eigenmode_{mne,dspm,sloreta}'`,
  `HeadModelFile`=composed, `SurfaceFile`=original → standard 3D cortex display.
- **Raw inputs** — kernel-only results node (coefficients require import, as today).

### 4.5 Failure handling

- Reject non-surface head models.
- Vertex-count mismatch between head model and eigenmodes (the surface-repair
  gotcha) → clear "recompute head model" message.
- Missing eigenmodes → "Compute eigenmodes first".
- Missing noise covariance → warn + identity whitening (discouraged for
  dSPM/sLORETA).

---

## 5. Testing & validation

Two parallel tracks: **(A)** implementation correctness, **(B)** scientific
benchmarking against Brainstorm defaults.

### Track A — Implementation correctness (`dev/tests/`)

Fast, deterministic, no DB (pure) + a couple of e2e, following the existing
`test_*_pure.m` / `test_*_e2e.m` convention.

| Test | Asserts |
|---|---|
| `test_eigenmode_leadfield_pure` | `L̃ = L·Φ` exactly on synthetic `L,Φ`; constrained-orientation extraction matches `bst_gain_orient`; `K` clamping; shapes. |
| `test_inverse_eigenmodes_prior_pure` | `R = diag(σ²(λ_k))` for log/flat/power; `λ̃∈(0,1)`; ratio-preservation invariant; `λ_1≈0` handling. |
| `test_inverse_eigenmodes_clean_pure` | Whitening matches `C`; SSP projector folding annihilates the projected-out subspace (`M̃·P_out≈0`); bad-channel drop changes dims correctly. |
| `test_inverse_eigenmodes_norm_pure` | MNE/dSPM/sLORETA normalization formulas vs analytic on synthetic. |
| `test_inverse_eigenmodes_harmonic_limit` | flat-`R` + regularization disabled (`λ→0`), `K=nChan` reproduces the retired harmonic `pinv` kernel within tol. |
| `test_eigenmode_roundtrip_e2e` | Plant a mode-space source, forward through `L̃`, reconstruct `Φ·M̃·d ≈ truth` at high SNR; nodes created with right shapes. |

### Track B — Scientific validation vs Brainstorm defaults

All levels compare **eigenmode-MAP vs wMNE / dSPM / sLORETA** on the same head
model, cortex, and noise covariance.

**Level 1 — Resolution-matrix analysis (analytic, no data).**
New helper `bst_resolution_metrics`: `Res = Kernel · L_full`, then per-vertex
point-spread metrics — localization error (peak displacement, mm), spatial
dispersion (focality), overall gain (depth bias), binned by cortical depth. Tests
the core GBF claim: better conditioning + smoothness prior ⇒ tighter PSFs and less
depth bias than wMNE.

**Level 2 — Simulation with ground truth (Brainstorm-native).**
`process_simulate_sources` → `process_simulate_recordings`: plant known scout
patches, forward-project, add sensor noise at set SNR; reconstruct all methods.
Metrics: DLE (mm), AUC, spatial dispersion, time-course correlation, RMSE. Sweep
depth (superficial vs deep — the depth-bias test), SNR (≈1/3/6/10 dB), source
extent, `#modes K`, prior type. Simulation runs on the **real OMEGA subjects'**
head models / sensor layouts / noise covariances so ground truth sits on realistic
geometry. Mirrors GBF's `demo_simulation`.

**Level 3 — Real data.**
- **Elekta phantom** (`tutorial_phantom_elekta`): physical dipoles at known mm
  positions → **absolute localization error** vs dSPM/LCMV. Hard real-data ground
  truth.
- **OMEGA tutorial (2 subjects)**: a focused qualitative comparison of the
  **GBF method vs dSPM** localization on real resting data (no known generator →
  consistency/plausibility, not accuracy).

### Harness & acceptance

A committed `toolbox/script/tutorial_eigenmodes_validation.m` runs Levels 1–3 and
writes a results markdown + figures (same pattern as the existing
`dev/tests/omega-icosphere-sourcemap-results.md`), runnable via the MATLAB MCP.

**"Validated" =**

- All Track-A tests pass; flat-`R` reproduces the harmonic kernel.
- Level 1: eigenmode-MAP PSF spatial dispersion ≤ wMNE, with reduced depth bias.
- Level 2: localization error / AUC competitive-or-better than dSPM/sLORETA across
  the SNR sweep, with less depth bias than wMNE.
- Level 3: phantom localization error not worse than dSPM; OMEGA maps localize
  plausibly relative to dSPM.

### Out of scope / future

- Cross-language cross-check against GBF's Python pipeline on identical geometry
  (guards against porting errors) — optional follow-up.

---

## 6. Open implementation details (resolved during planning)

- Exact `λ_ref` selection (`λ_{K+1}` vs `λ_K·(1+ε)`) and `ε`.
- Whether projector application reuses a shared helper or replicates the 2018
  inverse inline.
- Default `K` policy (clamp to `nChannels`, or allow `K > nChannels` leaning on the
  prior — GBF's regime). Note: GBF keeps `K=300 > nChannels` and relies on the
  prior; current Brainstorm clamps to `min(nChannels, available)`.
- dSPM/sLORETA normalization in mode space (carry over from current
  `bst_inverse_eigenmodes`, validated by `test_inverse_eigenmodes_norm_pure`).
