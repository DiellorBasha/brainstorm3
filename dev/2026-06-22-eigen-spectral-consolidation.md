# Eigen Spectral Consolidation — Plan & Preservation Notes

**Date:** 2026-06-22
**Author:** Diellor Basha (with Claude)
**Status:** Plan approved (decisions captured below); execution phased, deletions gated on review.

This document does two jobs:

- **Part I** — the consolidation/deprecation plan for retiring the legacy scalar
  Laplace–Beltrami (LBO) "eigenmode" path and folding its useful capability into the
  canonical **file-based** Dirac method.
- **Part II** — a rigorous preservation note (math + architecture) for every spectral
  function we deprecate, so the capability can be reborn as a single **operator-agnostic
  eigen-spectral umbrella** keyed on an `eigen_*.mat` file.

---

## North-star decisions (2026-06-22)

1. **One operator-agnostic eigen-spectral family.** Spectral operations (project,
   transform, spectrum fft/psd/freq, wavelet, dispersion, denoise/wiener/noisefloor,
   filter, prior) become functions of an **`eigen_` file**, which already carries the
   eigenvalue (λ) axis, the operator `Variant` (Dirac / Connection Laplacian /
   Laplace–Beltrami / …) in its `Provenance`, and `Phi`. They stop caring *which*
   operator produced the modes. `bst_dirac_eigenmodes_filter` is the existing template
   for the vector/quaternion basis pattern.
2. **Deprecate the legacy scalar-LBO spectral functions**, but **preserve their math
   and architecture** (Part II) for the future umbrella. Do not lose the logic.
3. **Drop the "eigenmode" label** in user-facing surfaces and naming; everything is
   "the file-based method" (manifold / operator / eigen nodes).
4. **Connection-Laplacian phase axis must be `eigen_`-node based**, not embedded in the
   surface file. It becomes a member of the same eigen-spectral family.
5. **Experimental inverse spurs stay, just isolated** from the canonical path (documented
   here, conceptually out of the main flow). `bst_inverse_dirac` is the canonical inverse.

---

# Part 0 — Incremental execution plan (supersedes the Tier ordering as the *how*)

The Tiers below (Part I) remain the **destination**. This section is the **route**: small,
independently-shippable, inline steps, each leaving the tree green, easiest first. The big-bang
refactor is explicitly rejected. Implement one increment at a time, validate, commit, move on.

## Guiding idea — give the new DB types first-class in/out, like every native Brainstorm type

`results_`, `timefreq_`, `headmodel_` etc. each have: a `db_template` schema, a `db_add_*`
(save+register) path, a `bst_get('…File')` resolver, and a dedicated **`in_bst_*` loader**.
Our new types (`eigen_`, `manifold_`, `operator_`) have the schema, `db_add_*`, and the
`bst_get` resolver — but **no dedicated loader**. Every consumer hand-rolls
`load(file_fullpath(...))` (≈14 sites: `bst_dirac`, `bst_dirac_filter`, `tess_eigen`,
`tess_manifold`, `process_dirac_filter`, `process_vortex_track`, `view_helmholtz`,
`view_eigenmode_spectrum`, `view_wavelet_designer`, `panel_wavelet_designer`). Centralizing
this is the foundation that later schema/validation/migration work hangs off.

**Out side:** `db_add_eigen/operator/manifold` already ARE the registered-save path (the
Brainstorm idiom — there is no `out_bst_results`; saving = `bst_save` + `db_add`). So do NOT
create `out_bst_*`; `db_add_*` is canonical. The legacy embedded pattern
(`in/out_tess_eigenmodes`, `in/out_tess_conn_eigenmodes` storing data *inside* the surface
`.mat`) is exactly what the node-file + `db_add` + `in_bst_*` pattern replaces.

## Increment ladder (do in order; each is a small, self-contained commit)

**Inc 1 — Add dedicated loaders (purely additive, zero risk).**
Create in `toolbox/io/`:
- `in_bst_eigen(EigenFile, fields...)`
- `in_bst_operator(OperatorFile, fields...)`
- `in_bst_manifold(ManifoldFile, fields...)`

Convention (mirror `in_bst_headmodel`/`in_bst_results`): accept a relative OR full path;
`file_fullpath` internally; if a struct is passed, return it unchanged (passthrough); optional
trailing field list → `load(f, fields{:})`; light schema check against `db_template` (presence
of `Variant`/`Phi`/`Lambda`/`K` for eigen; `Operator`/`Mass` for operator; the 5 geometry
groups for manifold). No callers changed yet — ship and test in isolation.

**Inc 2 — Route ad-hoc loads through the loaders (mechanical, one file per commit).**
Replace each `load(file_fullpath(X.OperatorFile))` → `in_bst_operator(X.OperatorFile)`, etc.
Order by simplicity: `tess_eigen:181`, `bst_dirac:235/281`, `bst_dirac_filter:71`,
`process_dirac_filter:122`, then the GUI/views (`view_helmholtz`, `view_eigenmode_spectrum`,
`view_wavelet_designer`, `panel_wavelet_designer`), then `tess_manifold:96`. Each diff is tiny
and independently verifiable.

**Inc 3 — Fold the "find a child node matching criteria" duplication into the loaders.**
`bst_dirac.local_find_dirac_eigen` and `tess_eigen.local_find_eigen` are near-duplicates
(scan surface `.Eigen`, match `Variant`/`K`/`Tau`). Provide a find mode, e.g.
`in_bst_eigen('Find', SurfaceFile, Variant, K, Tau)` returning `[EigenMat, EigenFile]`, and
have both call sites use it. Removes duplicated DB-scan logic; single source of truth.

**Inc 4 — Trivial dead-weight sweep (independent, can slot anywhere).**
Delete genuinely-dead `bst_face_eigenmode_leadfield` (0 production callers) and the dangling
no-definition test/benchmark stubs (`tess_dirac_eigenmodes`, `bst_eigenmodes_field`,
`bst_dirac_eigenmode_field`). No live code touched.

**Inc 5 — Connection-phase onto the `eigen_` node (medium; the first real "integrate" step).**
Route `bst_conn_phase`/`view_connection_phase` through `tess_eigen('Connection Laplacian')` +
`in_bst_eigen`, with a parity test vs the current surface-embedded `ConnEigenmodes`. Then retire
`in/out_tess_conn_eigenmodes`, `tess_conn_eigenmodes`, `bst_conn_eigenmodes_ensure`. (Untangle
`bst_conn_eigenmodes_ensure.m:52`'s call to `bst_eigenmodes_ensure` here.)

**Inc 6 — Build the operator-agnostic spectral umbrella** (`bst_eigen_spectral_*`, Part II) so
Tier-2 deprecations have a landing place. Largest unit; only after the I/O foundation is solid.

**Inc 7 — Tier 1 retirements**, then **Inc 8 — Tier 2 deprecations**, per Part I, once the
umbrella covers the capability and the math notes (Part II) are committed.

Recommended starting point: **Inc 1** (additive loaders) — or **Inc 4** if you want a quick
visible win first.

---

# Part I — Consolidation Plan (destination / tiers)

## The three storage mechanisms (root of the mess)

| Mechanism | Storage | Status |
|---|---|---|
| **Dirac / file-based nodes** | `eigen_*.mat`, `operator_*.mat`, `manifold_*.mat` via `tess_eigen`/`db_add_eigen` | ✅ Canonical |
| **Scalar LBO** | embedded in surface `.mat` (`Eigenmodes` field) via `in/out_tess_eigenmodes` | ❌ Legacy — retire |
| **Connection-Laplacian** | embedded in surface `.mat` (`ConnEigenmodes` field) via `in/out_tess_conn_eigenmodes` | ⚠️ Active research, wrong storage → migrate to `eigen_` node |

## KEEP — canonical Dirac path (no action)

`bst_dirac`, `bst_inverse_dirac`, `panel_inverse_dirac`, `process_inverse_dirac`,
`process_dirac_filter`, `bst_dirac_filter`, `bst_dirac_eigenmodes_filter`, `tess_eigen`,
`tess_operators`, `tess_manifold`, `db_add_{eigen,operator,manifold}`, `bst_eigs_smallest`,
`nxr_safe_create`, `bst_dirac_helmholtz{,_face}`, `bst_vortex_*`, `process_vortex_track`,
`view_manifold*`, `view_helmholtz`/`panel_helmholtz`, `panel_wavelet_designer`,
`view_leadfield_vectors`.

**The three "eigen" viewers are Dirac, despite their names** and are wired into the live
GUI (tree right-click, `bst_figures`, `figure_timeseries`/`figure_image`):
`view_eigenmodes` (vector-field viewer; explicitly rejects an LBO file),
`view_eigenmode_spectrum`, `view_eigen_timeseries`. Keep all three. (Rename away from the
"eigenmode" label per decision #3 as a later cosmetic pass.)

## DO NOT DELETE — shared / load-bearing despite legacy appearance

- **`bst_eigenmodes_filter_gain`** — pure transfer-function kernel **reused by the Dirac
  filter** (`bst_dirac_eigenmodes_filter.m:101`). Keep; becomes the gain hub of the umbrella.
- **`manifold_ft` / `manifold_ift`** — operator-agnostic primitives; foundation of the umbrella.
- **Connection-Laplacian axis** (`tess_conn_eigenmodes`, `bst_conn_eigenmodes_ensure`,
  `in/out_tess_conn_eigenmodes`, `bst_conn_phase`, `view_connection_phase`) — live
  research, migrate to `eigen_` node (see Integration task), do not delete.
- **`tess_tangents` / `bst_tangent_face2vertex`** — self-marked deprecated but still
  load-bearing for `bst_face_leadfield` (constrained/loose), `bst_wavefront_track`,
  `tess_nxr_populate`. Keep until the face leadfield constrained path migrates to `tess_manifold`.
- **`bst_helmet_eigenmodes`** — zero production callers; sensor-side tool. Orphaned but
  harmless; keep as a research utility (isolated).

## Tier 0 — dead, safe to delete now

- `bst_face_eigenmode_leadfield` (0 production callers)
- `bst_eigenmode_oscillator_inverse` (0 production callers) → but see decision #5: it is an
  experimental spur; *isolate rather than delete* unless confirmed truly dead. (It is also
  orphaned, so deletion is defensible; keep per #5 to be safe — move/document instead.)
- `bst_eigenmode_poisson_sharpen` (0 production callers) → same caveat as above.
- Stale no-definition references (delete the dangling test/benchmark stubs only):
  `tess_dirac_eigenmodes`, `bst_eigenmodes_field`, `bst_dirac_eigenmode_field`.

> Per decision #5 ("keep all, just isolate"), the orphaned experimental inverses are
> *retained and isolated*, not deleted. Tier 0 therefore reduces to: the genuinely dead
> `bst_face_eigenmode_leadfield` and the dangling no-def test/benchmark stubs.

## Tier 1 — retire superseded core (Dirac equivalent already exists)

Excise the dormant scalar-LBO branches inside live files first:
- `panel_headmodel` — `isEigenSpace` path (~488–491, 521, 865–877): `in_tess_eigenmodes` +
  `bst_eigenmode_leadfield`.
- `process_inverse_2018` — `'eigenmode'` branch (~773, 787–789): `bst_inverse_eigenmodes` +
  `bst_eigenmode_reconstruct`.
- `panel_inverse_2018` — permanently-hidden Eigen-dSPM radio (~173–184).

Then delete (all reachable only from each other / their process wrappers / the excised branches):
- Forward/inverse: `bst_eigenmode_leadfield`, `process_eigenmode_leadfield`,
  `bst_inverse_eigenmodes`, `bst_eigenmode_reconstruct`, `bst_eigenmode_prior`,
  `process_eigenmodes_inverse`.
- Compute/IO/GUI for the *scalar* surface-embedded basis: `tess_eigenmodes`,
  `in_tess_eigenmodes`, `out_tess_eigenmodes`, `bst_eigenmodes_ensure`,
  `panel_eigenmodes_compute`, `process_eigenmodes`, `process_eigenmodes_view`.
  - **Untangle first:** `bst_conn_eigenmodes_ensure.m:52` calls `bst_eigenmodes_ensure`
    only to derive a default mode count — replace with a local count before deleting.
- Spatial filter duplicate of the Dirac filter: `bst_eigenmodes_filter`,
  `process_eigenmodes_filter`, `process_eigenmodes_coeffsfilter`.

## Tier 2 — deprecate spectral-analysis family, preserve math (Part II)

Deprecate (after Part II notes are committed): `bst_eigenmodes_project`,
`bst_eigenmodes_transform`, `bst_eigenmodes_modekernel`, `bst_eigenmodes_wavelet`,
`bst_eigenmodes_wiener`, `bst_eigenmodes_dispersion`, `bst_eigenmodes_noisefloor`, and
their process wrappers `process_eigenmodes_{transform,denoise,dispersion,fft,freq,psd,
spectrum,wavelet,wiener}`, plus `tutorial_eigenmodes_validation` and the eigenmode arm of
`tutorial_benchmark_esi`.

These return as members of the operator-agnostic umbrella (Part II synthesis). Keep
`manifold_ft/ift` and `bst_eigenmodes_filter_gain` (the umbrella's primitives).

## Integration task — connection-phase onto the `eigen_` node

`tess_eigen('Connection Laplacian')` + `tess_operators('Connection Laplacian')` already
exist and produce the per-hemisphere complex Hermitian eigenbasis. Route
`bst_conn_phase` / `view_connection_phase` through a loaded `eigen_` node (via `db_add_eigen`)
instead of `in/out_tess_conn_eigenmodes`; then retire `tess_conn_eigenmodes`,
`bst_conn_eigenmodes_ensure`, `in/out_tess_conn_eigenmodes`. **Verify** the file-based
`'Connection Laplacian'` modes are bit-compatible with what `bst_conn_phase` expects
(`ConnEig.Vectors [nV x nModes]` complex) before retiring the surface-embedded path.

## Isolation task (decision #5)

Document the experimental inverse spurs as **off the canonical path** (no deletion):
`bst_eigenmode_analytic_inverse`, `bst_eigenmode_oscillator_inverse`,
`bst_eigenmode_poisson_sharpen`, `bst_cwt_fiedler_pipeline` (+ `bst_eigenmode_cwt_inverse`),
`bst_wavefront_track`, `bst_benchmark_inverse`, `bst_helmet_eigenmodes`. The canonical
inverse is `bst_inverse_dirac`.

## Suggested execution order

1. Commit this document (Part II is the deletion gate).
2. Tier 0 (genuinely dead + dangling stubs).
3. Tier 1 dormant-branch excision → then Tier 1 deletions (untangle the conn-ensure count first).
4. Tier 2 deprecations.
5. Integration (connection-phase → `eigen_` node) with a parity test.
6. Cosmetic: drop the "eigenmode" label from kept Dirac viewers/processes.

---

# Part II — Preservation Notes (math + architecture)

> Full mathematical and numerical reference for every deprecated spectral function, plus
> the proposed operator-agnostic umbrella design. **Do not delete the source functions
> until the new family is validated against these formulas.**

## Core primitives (already operator-agnostic — KEEP)

### `manifold_ft` — forward manifold Fourier transform
- Purpose: project a vertex field onto an M-orthonormal eigenbasis.
- I/O: `Phi [nV×K]`, `M [nV×nV]`, `U [nV×nT]` → `C [K×nT]`.
- Math: `C = Phi' * (M * U)`, with `Phi' M Phi = I_K`.
- Notes: left-multiplies M once (sparse×dense); assumes Phi pre-selected/ordered. Fully general.

### `manifold_ift` — inverse manifold Fourier transform
- I/O: `Phi [nV×K]`, `C [K×nT]` → `U [nV×nT]`.
- Math: `U = Phi * C` (no mass needed). `ift∘ft` = M-orthogonal projector onto span(Phi).

## Transforms / kernels (GENERALIZE)

### `bst_eigenmodes_project`
- Project scalar data, optionally reconstruct over a canonical mode range.
- Math: `Coeffs = manifold_ft(Phi, M, Data)`; `Recon = manifold_ift(Phi(:,iSel), Coeffs(iSel,:))`,
  `iSel = Order(k1:k2)` for permutation-agnostic canonical (ascending-λ) ordering.
- Keep: `.Order` field for reorder; clamp ModeRange.

### `bst_eigenmodes_transform` — unregularized sensor→mode transform
- Math: `L_tilde = Gain*Phi` `[nCh×K]`; `[U,S,V]=svd(L_tilde,'econ')`;
  `Tol = max(size)·eps(max(s))`; `Kernel = V·diag(1./s | s>Tol)·U'` `[K×nCh]`;
  `Theta = Kernel*Data`.
- Why SVD not normal equations: correct for both K≤nCh and K≥nCh; avoids squaring the
  (already large) condition number of Gain*Phi. Unregularized but rank-safe.

### `bst_eigenmodes_modekernel`
- Fold `Phi'·M` into an imaging kernel: `Projector = Phi(:,1:K)'·M`;
  `ModeKernel = Projector·ImagingKernel` (or bare `Projector` if kernel empty).
- Precomputes the projector once (vs. per-time-point in project). Fully general.

## Filtering / spectral analysis

### `bst_eigenmodes_filter_gain` — per-mode transfer function `h(λ)` (KEEP — umbrella hub)
- Types & math:
  - lowpass: `h_k = 1[k ≤ CutoffMode]`; highpass: `1[k ≥ CutoffMode]`;
    bandpass: `1[k1 ≤ k ≤ k2]` (index masks — assume monotone λ ordering).
  - heat: `h_k = exp(-t·λ_k)`; inverse_heat: `min(exp(+t·λ_k), MaxGain)`;
    tikhonov: `1/(1+β·λ_k)`; custom: `TransferFn(λ)`.
  - kernel-based types delegate to `bst_eigfilter_kernel`/`bst_eigfilter_evaluate`.

### `bst_eigenmodes_filter` — vertex-field spectral filter (GENERALIZE)
- Math: `Filtered = Phi · diag(h) · (Phi'·M) · Data`. Three-step project→scale→reconstruct,
  zero-phase (real, time-independent gains).

### `bst_eigenmodes_wavelet` — complex Morlet CWT of coefficients (REUSE)
- Math: `t=(0:nT-1)/sfreq`; `W = conj(morlet_transform(Coeffs, t, Freqs, Fc, FwhmTc, 'n'))`
  → complex `[K×nT×nFreq]`. Default `Freqs = logspace(log10(2), log10(min(100,0.4·sfreq)), 40)`.
  `|W|`=amplitude, `arg(W)`=phase (conj → ω>0 convention). Operator-independent.

### `bst_eigenmodes_dispersion` — wave vs diffusion (GENERALIZE w/ caution)
- Per mode: weight `p_k=Σ_f P(k,f)`; peak `f*_k=argmax_f P`; width
  `w_k=sqrt(Σ_f P f²/p_k − (Σ_f P f/p_k)²)`.
- Wave fit (through origin, power-weighted): `f*_k = a·sqrt(λ_k)`, `c = 2π·a`,
  `R2wave = wR²(sqrt(λ), f*, p)`.
- Diffusion fit: `w_k = b·λ_k`, `α = 2π·b`, `R2diff = wR²(λ, w, p)`.
- Decision: `Regime = 'wave' if R2wave ≥ R2diff else 'diffusion'`; `Margin=|R2wave−R2diff|`.
- Operator note: `c [m/s]`, `α [m²/s]` assume LBO on a metric surface; R²/Regime are
  unit-free and general.

### `bst_eigenmodes_noisefloor` — joint (λ,ω) SNR + spectral subtraction (REUSE)
- `SNR = Pdata./Nnoise`;
  `CleanPSD = max(Pdata − α·Nnoise, Floor·Nnoise)`;
  `Gain = max( max(Pdata − α·Nnoise, 0)./Pdata, GainFloor )` (≡ CleanPSD/Pdata clamped to [GainFloor,1]);
  `K*(f) = max{ k : SNR(k,f) ≥ SnrThresh }`.
- Power-domain (not complex) because data/empty-room noise are independent realizations.

### `bst_eigenmodes_wiener` — frequency-domain Wiener filtering of coefficients (REUSE)
- Optional mirror `[fliplr, x, fliplr]`; folded freq grid `fvec` folded at Nyquist;
  `H(k,f)=clamp(interp1(GainFreqs, Gain(k,:), fvec, 'linear', edge-hold), 0, 1)`;
  `Y = real(ifft(fft(x).*H))`; trim mirror. Zero-phase (H real, Hermitian via folding).

### `bst_eigenmode_prior` — diagonal source-covariance prior `R(λ)` (GENERALIZE w/ caution)
- flat: `R=1`; power: `R_k ∝ λ_k^{-α}` (via eigfilter 'power'); log (GBF mm-scale):
  `λ_mm = λ_m·1e-6`, `R_k = −log(λ_mm)`. DC guard: if `λ(1) ≤ 1e-12·max(λ)` set `λ(1)=λ(2)`.
  Normalize `R/=max(R)`. The 'log' mm rescaling is LBO/metric-specific.

## Reference template — `bst_dirac_eigenmodes_filter` (the vector pattern)

Already operator-agnostic in structure. Differences vs scalar that the umbrella must
encapsulate per `Variant`:

| Aspect | Scalar LBO | Dirac vector |
|---|---|---|
| Embedding | scalar ∈ ℝ | 3-vector → pure-imaginary quaternion `[0,x,y,z]` |
| Mass | `M [nV×nV]` | `B_h = kron(M_h, I4) [4nVh×4nVh]` per hemi |
| Eigvecs | `Phi [nV×K]` | `Phi_h [4nVh×K]` per hemi (4-fold multiplet) |
| Project | `C = Phi'·M·U` | `C = Phi'·B·embed(J)` |
| Scale | real diag `h(λ)` | real diag `h(λ)` (uniform across multiplet) |
| Extra | — | chirality/helicity projector (optional) |
| Reconstruct | `Phi·C` | `Phi·C` then take imaginary parts → `(x,y,z)` |

Shared & already general: eigenvalue/gain logic (`bst_eigenmodes_filter_gain`),
per-hemisphere cell processing, cell-based multi-component storage.

## Proposed umbrella (keyed on `eigen_` file)

`eigen_*.mat` carries: `Phi {1×2}`, `Lambda {1×2}`, `Variant`, `K`, `Order` (optional),
`GlobalVertices {1×2}`, `GlobalFaces {1×2}` (face variants), `Provenance`.

Three top-level functions, each dispatching internally on `Eigen.Variant` and hiding the
embedding/extraction:

1. **`bst_eigen_spectral_transform`** — subsumes `project` + `transform` + `modekernel`.
2. **`bst_eigen_spectral_filter`** — subsumes `filter` (project → `bst_eigenmodes_filter_gain` → reconstruct).
3. **`bst_eigen_spectral_analysis`** — subsumes `wavelet` + `dispersion` + `noisefloor` + `wiener` + `prior`.

Principles: single gain hub; per-hemisphere cells everywhere; embedding/extraction behind
the dispatcher; operator-specific scales (e.g. mm log-prior) encapsulated/recorded in
`Provenance`; `manifold_ft/ift` as the scalar core.
