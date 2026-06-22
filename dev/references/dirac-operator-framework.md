# The Dirac-Operator Framework for Vector MEG/EEG Source Mapping and Flow Analysis

**Status:** Reference / implementation guide
**Scope:** How the Dirac operator and its eigenmodes are computed, used to map MEG/EEG
sources as cortical **vector fields**, and then used to **analyse** those recovered fields
(Helmholtz–Hodge decomposition, vortices, chirality).
**Code root:** `toolbox/` of this Brainstorm fork. File:line references are given throughout.

This document is the connective tissue between the conceptual proposals
([`vector_eigenmode_meg_proposal.md`](vector_eigenmode_meg_proposal.md),
[`inverse-methods-as-spectral-filters.md`](inverse-methods-as-spectral-filters.md),
[`face_based_source_model.md`](face_based_source_model.md),
[`Wang2026_Geometry_Aware_Brain_Dynamics.md`](Wang2026_Geometry_Aware_Brain_Dynamics.md))
and the implementation that now lives in the toolbox.

---

## Contents

1. [Why a Dirac operator](#1-why-a-dirac-operator)
2. [Mathematical setting: vectors as quaternion fields](#2-mathematical-setting-vectors-as-quaternion-fields)
3. [Step 1 — The operator node (`tess_operators`)](#3-step-1--the-operator-node-tess_operators)
4. [Step 2 — The eigenmodes (`tess_eigen`)](#4-step-2--the-eigenmodes-tess_eigen)
5. [Step 3 — Forward: leadfield → mode basis (`bst_dirac` Transform)](#5-step-3--forward-leadfield--mode-basis-bst_dirac-transform)
6. [Step 4 — Inverse: whitened MNE in the mode basis (`bst_inverse_dirac`)](#6-step-4--inverse-whitened-mne-in-the-mode-basis-bst_inverse_dirac)
7. [Step 5 — Reconstruct & view the source vector field](#7-step-5--reconstruct--view-the-source-vector-field)
8. [Spectral filtering in the Dirac eigenbasis](#8-spectral-filtering-in-the-dirac-eigenbasis)
9. [Analysing the recovered field — the first-order Dirac](#9-analysing-the-recovered-field--the-first-order-dirac)
10. [Benchmarks & validation figures](#10-benchmarks--validation-figures)
11. [Provenance, caveats, open issues](#11-provenance-caveats-open-issues)

---

## 1. Why a Dirac operator

Standard distributed MEG/EEG inverse methods either (a) constrain each source to the local
surface **normal** (one scalar per vertex), or (b) solve an **unconstrained** 3-component
dipole per vertex with no geometric coupling between neighbours. Both are problematic for
*vector* analysis of cortical currents:

- The cortical normal flips by ~180° across every sulcal wall, so a smooth tangential
  current acquires spurious π phase jumps in the normal-constrained representation.
- The unconstrained representation has no notion of a **consistent frame**, so "the same"
  flow on two sides of a fold is represented by unrelated 3-vectors.

The **Dirac operator** on the cortical surface provides a geometry-aware basis for tangential
vector fields: its eigenmodes are smooth quaternion-valued fields that transport correctly
across folds (the gauge is built in). Expanding the leadfield in this basis yields an inverse
that is (i) a curvature-aware spectral prior and (ii) returns a genuine cortical **vector
field**, which can then be decomposed with the *first-order* Dirac into irrotational /
solenoidal / harmonic flow for vortex and wave analysis.

Two distinct uses of "Dirac" appear in the code and must not be conflated:

| Object | Operator | Used for | Dimensions (per hemisphere) |
|---|---|---|---|
| **Squared Dirac** `A` | `OperatorMat.Operator{h}` | the *eigenbasis* for source mapping | `[4nVₕ × 4nVₕ]` |
| **First-order Dirac** `D` | `OperatorMat.FirstOrder.Intrinsic{h}` | *differential analysis* (div/curl) of a field | `[4nFₕ × 4nVₕ]` |

The eigenbasis (Steps 1–5) uses the squared operator; the flow analysis (Section 9) uses the
first-order operator.

---

## 2. Mathematical setting: vectors as quaternion fields

A tangential current at a vertex is a 3-vector `(Jx, Jy, Jz)`. The Dirac framework embeds it
as a **pure-imaginary quaternion** `(w, i, j, k) = (0, Jx, Jy, Jz)`. A field over `nV`
vertices therefore becomes a `4nV`-vector with the layout (rows interleaved per vertex):

```
ψ(1:4:end) = w  = 0          % real/scalar slot stays zero for a current field
ψ(2:4:end) = i  = Jx
ψ(3:4:end) = j  = Jy
ψ(4:4:end) = k  = Jz
```

This embedding is used identically in the forward transform
(`bst_dirac.m:136-140`), the inverse, the spectral filter
(`bst_dirac_eigenmodes_filter.m:94-98`), and the flow analysis
(`bst_dirac_helmholtz.m`). Recovering a 3-vector is the reverse: drop the `w` slot and read
`(i, j, k)`.

Quaternion **multiplication** is what makes this more than bookkeeping: the squared Dirac
operator's `4×4` blocks are quaternion left-multiplications by edge vectors, so the operator
mixes the three spatial components according to the surface geometry (the cotan Laplacian is
exactly the `w`–`w` scalar part; the off-diagonal blocks carry the tangent-frame transport,
`tess_operators.m:298-323`).

---

## 3. Step 1 — The operator node (`tess_operators`)

**File:** `toolbox/anatomy/tess_operators.m`. **Backend:** the `nxr-compute` plugin
(geometry-central). **Entry:** `tess_operators(SurfaceFile, Variant, 'Tau', τ)`.

`tess_operators` builds and caches a derived-anatomy **operator node** (`operator_*.mat`)
nested under a cortical surface. Three variants exist (`tess_operators.m:14-20`):

| Variant | `A = Operator{h}` | `B = Mass{h}` |
|---|---|---|
| `'Laplace-Beltrami'` | cotan Laplacian `[nVₕ × nVₕ]` | Galerkin mass `[nVₕ × nVₕ]` |
| `'Connection Laplacian'` | Levi-Civita connection Laplacian `[nVₕ × nVₕ]`, Hermitian | Galerkin mass |
| `'Dirac'` | squared Dirac `[4nVₕ × 4nVₕ]` | `kron(Mass, I₄)` `[4nVₕ × 4nVₕ]` |

**The squared Dirac** (`tess_operators.m:219`):

```
A = (1 − τ)·(D²_int / s_L)  +  τ·(E / s_E)            [4nVₕ × 4nVₕ]
```

- `D²_int = D'·M_F·D` — the *intrinsic* (immersion-based) squared Dirac, built from the
  first-order operator `D` and the face mass `M_F = kron(face-area, I₄)`
  (`tess_operators.m:298-323`). Its scalar (`w`) part **is** the cotan Laplacian.
- `E` — the *extrinsic* (Gauss-map) squared Dirac from nxr (`nxr_compute('operators',h,'dirac',1)`).
- `s_L = λ_max(D²_int, B)`, `s_E = λ_max(E, B)` — per-block largest generalized eigenvalues
  (`local_lambda_max`, `tess_operators.m:281-295`); they co-normalize the two blocks so that
- `τ ∈ [0,1]` is a **dimensionless** mixing weight (default `0.5`), portable across mesh
  resolution and units. `τ=0` = purely intrinsic; `τ=1` = purely extrinsic.

**Stored fields** (`tess_operators.m:256-266`):

```
.Variant         'Dirac'
.Operator        {A_L, A_R}                 squared Dirac, [4nVₕ × 4nVₕ]
.Mass            {B_L, B_R}                 kron(mass, I₄)
.GlobalVertices  {vH_L, vH_R}               global vertex indices of each hemisphere
.FirstOrder      .Intrinsic{1×2}            first-order D, [4nFₕ × 4nVₕ]  (← used in Section 9)
                 .Extrinsic{1×2}
.FaceMass        {W_F,L, W_F,R}             kron(face-area, I₄), [4nFₕ × 4nFₕ]
.Provenance      Backend='nxr', Tau, DiracScale={[s_L s_E]_L, [s_L s_E]_R}, ...
```

Everything is **block-diagonal by hemisphere** (left = index 1, right = 2); the operator and
eigensolve never couple the hemispheres. `GlobalVertices{h}` maps a hemisphere's local indices
to global vertex numbers.

> **Note (mass matrices).** The Dirac node's `Mass{h}` is the `4nVₕ` quaternion mass
> `kron(Mass_h, I₄)`. The *scalar* `nVₕ×nVₕ` mass used by the Poisson solves in Section 9
> comes from the separate **Laplace-Beltrami** node (`LBO.Mass{h}`).

---

## 4. Step 2 — The eigenmodes (`tess_eigen`)

**File:** `toolbox/anatomy/tess_eigen.m`. **Entry:** `tess_eigen(SurfaceFile, Variant, 'K', K, 'Tau', τ)`.
This is the find-or-load-or-compute wrapper: it locates a cached `eigen_*.mat` node (reusing
it if the variant, `K`, and `τ` match), else solves the generalized eigenproblem on the
operator node (creating the operator first if needed — `tess_eigen.m:150-164`).

**The generalized eigenproblem** per hemisphere is
`A φ = λ B φ`, smallest `λ` first (`bst_eigs_smallest`, `tess_eigen.m:237`). For the scalar
LBO this is routine. For the **Dirac** operator there is a structural subtlety:

### 4.1 The 4-fold quartet structure

Because the field is quaternion-valued, every distinct geometric eigenvalue `λ` comes with a
**4-fold multiplet** (the four quaternion units span a degenerate eigenspace). `eigs` returns
a rank-deficient, non-`B`-orthonormal spanning set across each multiplet, so `tess_eigen`:

1. **over-fetches** (`nRequest = K + ~30%`, `tess_eigen.m:219-221`) to capture whole multiplets;
2. runs a **Rayleigh–Ritz degenerate recovery** (`local_ritz_basis`, `tess_eigen.m:391-429`):
   form the Gram `G = U'BU`, keep its well-conditioned directions, build a `B`-orthonormal
   span `W`, diagonalize `Lr = W'AW`, and return genuine eigenpairs `Φ = W·V_r`,
   `λ = diag(D_r)`.

The result is a `B`-orthonormal basis (`Φ' B Φ = I`) with `λ` stored **once per multiplet**
(no repeated eigenvalues in `Lambda`).

### 4.2 Stored fields (`tess_eigen.m:145-154`)

```
.Variant         'Dirac'
.OperatorFile    file_short path to the operator_ node solved
.Phi             {Φ_L, Φ_R}     Dirac eigenvectors, [4nVₕ × K], B-orthonormal
.Lambda          {λ_L, λ_R}     [K × 1] ascending eigenvalues
.K               number of modes stored (nested: a cached K′≥K is truncated)
.GlobalVertices  {vH_L, vH_R}
.Provenance      Ortho='Rayleigh-Ritz' (Dirac/Conn) | 'B-orthonormal' (LBO), Tau, K, ...
```

The eigenvalues `λ_k ≥ 0` order the modes from **coarse (small λ, smooth, global)** to
**fine (large λ, oscillatory, local)** — this is the spatial-frequency axis that the inverse
prior and all spectral filters act on.

---

## 5. Step 3 — Forward: leadfield → mode basis (`bst_dirac` Transform)

**File:** `toolbox/forward/bst_dirac.m`. **Entry:** `CompHM = bst_dirac(HeadModel, 'Transform', 'nModes', K, 'Tau', τ)`.

Takes a standard **unconstrained** surface head model (`Gain [nCh × 3nV]`) and projects it
into the Dirac eigenbasis. Per hemisphere (`bst_dirac.m:136-142`):

1. Embed each channel's `3nVₕ` gain row as a quaternion field `Ψ [4nVₕ × nCh]` (Section 2).
2. Project onto the eigenmodes using the `B`-inner product:

```
L̃_h = Ψ_h' · (B_h · Φ_h)          [nCh × K]          % "mode leadfield"
```

i.e. `L̃_h(c,k) = ⟨gain_c, φ_k⟩_B` — channel `c`'s response to eigenmode `k`. The two
hemispheres are concatenated into the composed head model (`bst_dirac.m:147-166`):

```
CompHM.Gain            [nCh × 2K]                % stacked [L̃_L , L̃_R]
CompHM.Eigenvalues     [2K × 1]   λ_k per mode
CompHM.ModeHemisphere  [2K × 1]   1 or 2
CompHM.isDiracEigenmode = 1
CompHM.DiracEigenFile  path to the eigen_ node (so Reconstruct can reload Φ)
CompHM.DiracTau = τ
```

**Reconstruct** (`bst_dirac.m:170-202`) is the inverse map: given mode coefficients/rows
`c [m × 2K]`, it computes `R = Φ·cₕ` per hemisphere and reads the vector part to return a
per-vertex field `J [m × 3nV]`. Because `Φ` is `B`-orthonormal with `B = kron(mass, I₄)`,
`‖c‖₂ = ‖J‖_B` (the mode-coefficient norm equals the mass-weighted current norm) — the basis
is an isometry between mode space and the cortical `L²` of currents.

---

## 6. Step 4 — Inverse: whitened MNE in the mode basis (`bst_inverse_dirac`)

**File:** `toolbox/inverse/bst_inverse_dirac.m`. This is a **whitened minimum-norm** solver
that operates in the Dirac mode basis, built and validated stage by stage. Stage 1 was checked
**bit-identical** to Brainstorm's `bst_inverse_linear_2018` whitener. The five stages:

| Stage | Operation | Key lines |
|---|---|---|
| 1. **Whiten** | per-modality noise-cov regularization (Hamäläinen mean-eigenvalue ridge), inverse whitener `iW`; `Gw = iW·G_mode` | `:125-193` |
| 2. **Observe** | SVD of the whitened mode leadfield `Gw = U·S·Vᵀ`; `S` = **observability** singular values; `λ = SNR / mean(S²)` | `:195-219` |
| 3. **Filter** | Wiener/Tikhonov window `g(s) = λs/(λs²+1)`; mode kernel `K_mode = V·diag(g)·Uᵀ` | `:221-232` |
| 4. **Reconstruct** | lift the right singular vectors to the cortex `W_res = Φ·V`; vertex kernel `K_vtx = W_res·diag(g)·Uᵀ` | `:234-248` |
| 5. **Normalize** | amplitude (MN) / dSPM (per-source noise norm) / sLORETA (3×3 resolution norm); fold the whitener back: `ImagingKernel = K_norm·iW` | `:250-290` |

The **three-axis** way to read this (the conceptual framework that organizes the method):

- **Noise axis** — Stage 1 whitening, geometry-blind, set by the noise covariance.
- **Observability axis** — Stage 2 SVD, the ~tens of degrees of freedom the sensor array can
  actually see (basis-invariant; it is MEG physics, not a choice of basis).
- **Geometry/prior axis** — the eigenvalue spectrum `λ` (and any prior `R(λ)`), the
  curvature-aware spatial prior unique to the eigenmode approach.

A Matérn/`(κ²+λ)^{-ν}` prior can be layered on the geometry axis, but its gains are marginal
because the observability ceiling dominates (see `bench_dirac_matern_prior.png`, Section 10).

**Output / results node** (`bst_inverse_dirac.m:295-312`, written by
`process_inverse_dirac.m:171-200`): a **shared kernel** results file with

```
.ImagingKernel      [3nV × nChan]    vertex kernel on raw data  → J(t) = ImagingKernel · M(t)
.ImagingKernelMode  [2K × nChan]     mode-coefficient kernel    → c(t) = ImagingKernelMode · M(t)
.nComponents = 3                     unconstrained (vector) source
.Eigenvalues, .ModeHemisphere, .DiracEigenFile, .DiracTau
.Comment            e.g. 'Dirac: dSPM: MEG'
```

**Interactive entry:** the process **"Compute sources: Dirac eigenmodes"**
(`process_inverse_dirac.m`) exposes `measure` (Current density / dSPM / sLORETA), `snr`,
`nmodes` (default 400/hemisphere), `tau` (default 0.5), `noisereg` (0.1), `sensortypes`
(default MEG). Use MEG-only (≈274 ch on CTF) — the observable subspace is what matters.

---

## 7. Step 5 — Reconstruct & view the source vector field

The shared kernel is applied to recordings to give the cortical **vector field**:

```
J(t) = ImagingKernel · M(t)              [3nV × nTime]
per-vertex current at v:  J([3v-2, 3v-1, 3v], t)
```

Because the result keeps the **mode kernel** too, the same estimate can be viewed in three
linked ways from the results node's tree popup:

- **On the cortex** — "Display on cortex" shows the 3-component field as a vector quiver +
  norm colormap (the native unconstrained-source display).
- **Eigenvalue image / time series** (`view_eigen_timeseries`) — the mode-coefficient
  continuum `c_k(t) = ImagingKernelMode · M(t)` as a `(λ × time)` image or stacked traces:
  the spectral view between the sensors and the cortex. Travelling waves appear as tilted
  stripes; diffusion as a downward drift of the energy centroid.
- **Eigenspectrum (single time)** (`view_eigenmode_spectrum`) — the instantaneous
  `|c_k|²` spectrum over `λ`.

See `analyze_alpha_dirac_spectrogram.png` and `analyze_alpha_dirac_scale.png` (Section 10) for
the alpha-band scale-spectrogram and the diffusion/standing-wave test built on `c_k(t)`.

---

## 8. Spectral filtering in the Dirac eigenbasis

Once the field lives in the mode basis, **filtering is multiplication by a transfer function
`g(λ)`** — a band-limit / smoothing / sharpening of the *spatial* spectrum.

### 8.1 The kernel library (`bst_eigfilter_kernel`)

**File:** `toolbox/math/eigfilter/bst_eigfilter_kernel.m`. A string-dispatched registry of
analytic kernels `g(λ)`; `bst_eigfilter_kernel('info', name)` returns a metadata struct
(`display`, `params`, `bandpass`, `priorAdmissible`).

| Kernel | `g(λ)` | Role |
|---|---|---|
| `heat` | `exp(−tλ)` | low-pass / diffusion (smoothing at scale ~√t) |
| `inverse_heat` | `min(exp(+tλ), g_max)` | sharpening (capped) |
| `tikhonov` | `1/(1+βλ)` | low-pass regularization |
| `power` | `λ^{−α}` | smoothness / `1/f` prior |
| `matern` | `(κ²+λ)^{−ν}` | curvature-aware GP/SPDE prior |
| `mexhat` | `(tλ)·exp(−tλ)` | **band-pass** (Mexican hat) — not a valid prior |
| `dog` | `exp(−t₁λ) − exp(−t₂λ)` | **band-pass** (difference of heat kernels), `t₁<t₂` enforced |

(The GUI sections currently expose `{mexhat, dog, heat, inverse_heat, tikhonov}` via
`panel_eigenfilter_design('Kernels')`.) The same library feeds both the **analysis filters** and
the inverse **prior** (`bst_eigenmode_prior`), per the design in
[`2026-06-02-eigfilter-library-design.md`](../2026-06-02-eigfilter-library-design.md) and the
unifying view in
[`inverse-methods-as-spectral-filters.md`](inverse-methods-as-spectral-filters.md).

### 8.2 The vector filter (`bst_dirac_eigenmodes_filter`)

**File:** `toolbox/math/bst_dirac_eigenmodes_filter.m`.
Signature `[JFilt, h, Coeffs] = bst_dirac_eigenmodes_filter(EigenMat, MassCell, J, FilterType, ...)`.
Per hemisphere it does exactly the project → scale → reconstruct cycle:

```
ψ  = embed(J)                       % pure-imaginary quaternion field
c  = Φ' · B · ψ                     % project onto the Dirac eigenmodes   [K × nT]
c  ← diag(g(λ)) · c                 % scale by the transfer function
JFilt = imag( Φ · c )              % reconstruct, drop the w slot
```

Because every member of a 4-fold multiplet shares one `λ` (and one `g(λ)`), the scale filter
acts cleanly on whole multiplets. An optional **chirality** projector
`P_±(n̂) = (I ∓ i·R_n̂)/2` (`'Chirality'` option) splits the field into right/left circular
polarizations about an axis — used by the Wavelet Designer's helicity control
(`demo_dirac_wavelet_chirality.png`).

### 8.3 Where the filter is used in the GUI

- **Wavelet Designer** — a *localized* atom: vertex δ + direction (in the local manifold
  frame) + kernel + chirality + spectrum tiling.
- **Helmholtz "Smoothing"** — the same filter, folded into the Helmholtz view (Section 9.4)
  as an active-frame low-pass so the vortex-core detection can be tuned by spatial scale.
  *(The earlier standalone "Spatial filter" tool has been folded into this.)*

---

## 9. Analysing the recovered field — the first-order Dirac

This is the second half of the framework: given the recovered vector field `J`, use the
**first-order** intrinsic Dirac `D` to read its differential structure.

### 9.1 Helmholtz–Hodge decomposition (`bst_dirac_helmholtz`)

**File:** `toolbox/math/bst_dirac_helmholtz.m`. On a surface, a tangential field splits as

```
J = ∇φ  +  ∇⊥ψ  +  h          (irrotational + solenoidal + harmonic)
```

The pipeline, per hemisphere (`Prepare` caches the operators + a Cholesky factor; `Frame`
decomposes one time frame on demand):

1. **First-order Dirac applied to the embedded field** (`bst_dirac_helmholtz.m:105-110`):
   `q = D·ψ` over faces. Reading the quaternion output:

   ```
   ω (vorticity)   = q(1:4:end)                       % the w-part, per face
   ∇·J (divergence) = Σ q(2:4:end·) · n̂_face          % imaginary part · face normal
   ```

   > **Important convention (validated empirically, this fork):** the **w-part is the
   > vorticity** and the **imaginary·n̂ is the divergence** — not the reverse. This was checked
   > against synthetic pure-gradient and pure-rotational fields; with the labels swapped, a
   > pure gradient (which is curl-free) excited the wrong channel. Getting this backwards is
   > exactly what made the vortex cores fail to track the swirls in an earlier iteration.

2. **Poisson solves on the scalar cotan Laplace–Beltrami** (`LBO.Operator{h}`, `LBO.Mass{h}`,
   pinned mean-zero, cached factor):
   `K ψ = M ω`  (stream function from vorticity);   `K φ = M (∇·J)`  (potential from divergence).

3. **Component vector fields** via a per-face FEM gradient: `Virr = ∇φ`, `Vsol = n̂ × ∇ψ`,
   `Vharm = J − Virr − Vsol` (exact residual). A harmonic-energy fraction is reported.

4. **Singular points** (component-aware): vortex **cores** = local extrema of `ψ` (signed by
   `ω`); **sources/sinks** = local extrema of `φ` (signed by `∇·J`)
   (`FindCores`). A vortex centre is *any* extremum of ψ (a max and a min are a vortex /
   antivortex pair of opposite handedness), not only the global maximum — see Section 9.3.

`bst_helmholtz_curl.png` and `bst_helmholtz_divergence.png` show the curl/vorticity and
divergence scalar fields of a real source map; `bst_tangential_sourcesink.png` the
source/sink structure.

### 9.2 A note on operator consistency (why the Dirac div/curl is kept)

An operator-consistent **discrete-Hodge** variant was prototyped — div/curl from the FEM
gradient's own adjoint (`A = G'·diag(2·area)·G`), which gives cleaner component vector fields
(harmonic residual ~76% → ~27% on real data). But its cores/markers were less interpretable
than the **Dirac-based** div/curl, which detect the vortex cores well even though the
FEM-gradient component vectors leak more energy into the harmonic residual. The shipped code
therefore uses the **Dirac** div/curl (with the corrected w↔vorticity convention). The
consistent-projection prototype is preserved in git history for future revisiting.

### 9.3 Why there are many cores

`FindCores` marks every **local** extremum of ψ over its 1-ring (not the single global max/min).
A multi-vortex flow has a ψ landscape with many hilltops and basins — each is a vortex — so
multiple maxima/minima is correct, and forced by Morse / Poincaré–Hopf
(`#max − #saddle + #min = χ`). The raw count over-detects on noisy frames (every wrinkle in ψ
is a tiny extremum), which is why smoothing (Section 9.4) and a magnitude gate are provided.

### 9.4 The interactive view (`view_helmholtz` / `panel_helmholtz`)

Launched from the Dirac source **results node** popup ("Helmholtz / vorticity (Dirac)"). It
opens the **native** unconstrained-source figure and adds a control panel:

- **Component** radio — Total `|J|` / Irrotational `∇φ` / Solenoidal `∇⊥ψ` / Harmonic `h`.
  Switching changes the **scalar colormap** (signed `stat2` for φ/ψ; one-sided for `|J|`/`|h|`)
  while the **total-field quiver stays on** for every component.
- **Smoothing** (Dirac eigenmode kernel + scale) — band-limits the active frame *before*
  decomposing, so the cores follow the chosen spatial scale (66 → 34 cores in one real frame).
- **Marker threshold** — drops cores/sources with weak `|ω|`/`|div|` (34 → 2).
- **Show vectors / Show singular points** toggles; the quiver respects the Surface panel's
  **Data threshold** (Amplitude) slider rather than a hardcoded gate.

It is **active-frame, on-demand**: each frame is filtered → decomposed → cored as the time
cursor moves (the cotan factor is cached). Whole-series / batch decomposition is intended as a
separate `process_*` development.

The flow phenomenology this enables is shown in `dirac_alpha_vortex_pair.png`,
`dirac_alpha_vortex_timecourse.png`, `dirac_vortex_spin_vs_vorticity.png`,
`dirac_helmholtz_vortexscale.png`, and `dirac_superiorparietal_vortex*.png`.

---

## 10. Benchmarks & validation figures

All figures are in `dev/benchmarks/`; descriptions below are grounded in the generating
scripts where one exists (`compare_*.m`, `bench_*.m`, `analyze_*.m`, `demo_*.m`) and on the
filename otherwise.

**Eigenmode source mapping vs standard methods**
- `compare_dirac_vs_mne.png`, `compare_dirac_vs_mne_MMN.png` (`compare_dirac_vs_mne.m`) —
  Dirac-dSPM vs standard dSPM on real evoked data (M100, auditory MMN): spatial-map
  correlation, peak-vertex distance, Dice overlap, peak timecourse correlation. Dirac matches
  MNE/dSPM/sLORETA on real M100/MMN.
- `dirac_inverse_stage4_localization.png`, `dirac_inverse_stage5_dspm.png`
  (`bench_dirac_source_mapping.m`) — staged localization accuracy and dSPM maps.
- `dirac_leadfield_eigvec_analysis.png`, `compare_leadfield_*.png`
  (`compare_mne_dirac_leadfield.m`) — leadfield singular values and eigenvector agreement
  between the standard and Dirac bases.
- `eigenmode_accuracy_run/REPORT.md` + `figures/f1…f5` — the synthetic-on-real-cortex harness
  (`eigenmode_accuracy_benchmark_design.md`): LocError / Correlation / NRMSE / AUC across
  anatomy × regime × SNR × K. This harness **caught the `ModeIndices` reconstruction bug**
  (naive `Phi(:,1:K)` scrambled estimates → ~80 mm error); `eig` ties wMNE on distributed
  sources, with focal accuracy capped by the observability ceiling.

**Forward / observability / sensor geometry**
- `bench_dirac_forward_observability.png`, `dirac_noise_whitened_observability.png` —
  observability `S²` vs eigenvalue `λ` (raw and noise-whitened): how the sensor array filters
  spatial frequencies.
- `bench_sensor_dirac.png`, `bench_sensor_dirac_vsh.png` (`bench_sensor_dirac.m`) — helmet
  (sensor-side) Laplacian spectrum and VSH internal/external separation.

**Prior design**
- `bench_dirac_matern_prior.png` (`bench_dirac_matern_prior.m`) — Matérn prior helps
  *extended* sources (patch Dice / centroid) but not *focal* localization; quantifies the
  marginal prior gains under the observability ceiling.
- `demo_mne_dirac_filters.png` — filter gain `g(λ)` and the vertex point-spread (delta
  response) for several kernels.

**Alpha-band scale analysis**
- `analyze_alpha_dirac_spectrogram.png`, `analyze_alpha_dirac_scale.png`
  (`analyze_alpha_dirac.m`) — the `(λ × time)` scale-spectrogram of `c_k(t)`, the
  λ-centroid diffusion test, and the characteristic spatial scale over a posterior-alpha burst.

**First-order Dirac / flow analysis**
- `demo_dirac_firstorder_flux.png`, `demo_dirac_heat_transport.png` (frame-transport: Dirac
  heat flips the ambient field ~180° across a fold while a naïve scalar heat does not),
  `demo_dirac_wavelet_chirality.png` (helicity split of a vector wavelet).
- `bst_helmholtz_curl.png`, `bst_helmholtz_divergence.png`, `bst_tangential_sourcesink.png`,
  `dirac_rotation_rate.png`, `dirac_basis_ambient_vs_normal.png`,
  `dirac_timeresolve_chirality_div.png`.
- Vortex phenomenology: `dirac_alpha_vortex_pair*.png`, `dirac_alpha_vortex_timecourse.png`,
  `dirac_alpha_vortex_flat.png`, `dirac_vortex_spin_vs_vorticity.png`,
  `dirac_vortex_scale_spectrum.png`, `dirac_helmholtz_vortexscale.png`,
  `dirac_superiorparietal_vortex*.png`, `bst_vortex_*` (per-session top/zoom, advection
  montage, scale-isolated).

**GUI** — `filterdesigner_panel.png`, `wavelet_seed_marker.png`, `view_manifold_*.png`,
`manifold_v2_*.png`, `explore_alpha_*.png` document the viewers and control panels.

---

## 11. Provenance, caveats, open issues

- **Backend.** Every operator ultimately comes from the `nxr-compute` plugin (geometry-central).
  `tess_operators` is the *only* entry point that calls `nxr_compute('create', …)`; the exact
  cotan/connection/Dirac construction inside nxr is C++ and not visible from the toolbox.
  Use `bst_canonical_cortex` meshes for `create` (hand-built meshes can segfault nxr).
- **Hemisphere blocks.** All operators/eigenmodes are block-diagonal by hemisphere; a "mode"
  lives on one hemisphere only. The forward/inverse stack the two hemispheres' `K` modes into
  `2K` and carry `ModeHemisphere` throughout.
- **`τ` and caching.** Eigen nodes are reused only when `Variant`, `K`, **and** `τ` match;
  changing `τ` recomputes. A cached `K′ ≥ K` is truncated (the basis is nested in ascending λ).
- **Observability ceiling.** The honest ceiling on accuracy is the ~tens of observable DOF
  (Stage 2 SVD), set by MEG physics; richer priors (Matérn) and more modes mostly refine the
  null-space, not the data-space, fit.
- **Div/curl convention.** w-part = vorticity, imag·n̂ = divergence (Section 9.1). This is the
  single most error-prone point; the synthetic gradient/rotational test in
  `dev/tests/test_dirac_helmholtz.m` guards it (`irrotational is divergence-dominated`,
  `solenoidal is curl-dominated`).
- **Core over-detection.** Local-extrema detection is scale-sensitive; the smoothing scale +
  magnitude gate are the intended controls. A topological-persistence detector is a possible
  future upgrade.
- **Batch flow analysis.** The Helmholtz view is active-frame only; a whole-series/batch
  `process_dirac_helmholtz` is a planned separate development.

### Related reference documents (same folder)

- [`vector_eigenmode_meg_proposal.md`](vector_eigenmode_meg_proposal.md) — the motivating
  proposal (trivial connections, vector heat method, Hodge decomposition, joint (λ, ω) analysis).
- [`inverse-methods-as-spectral-filters.md`](inverse-methods-as-spectral-filters.md) — MNE /
  LORETA / dSPM / sLORETA / LCMV as transfer functions `T_k = σ_k²/(σ_k² + α r_k)` in eigenmode
  space.
- [`face_based_source_model.md`](face_based_source_model.md) — face-flux (2-form) source model.
- [`spatiotemporal_filters_proposal.md`](spatiotemporal_filters_proposal.md) — joint (λ, ω)
  spatiotemporal filter banks.
- [`Wang2026_Geometry_Aware_Brain_Dynamics.md`](Wang2026_Geometry_Aware_Brain_Dynamics.md) —
  the GBF empirical validation of eigenmode source imaging (log prior).
- [`diffusionnet-cortical-aplications.md`](diffusionnet-cortical-aplications.md) — learned heat
  diffusion on surfaces (geometry-aware, mesh-agnostic).

---

*Generated 2026-06-16. File:line references are to this fork's `toolbox/`. Figure descriptions
are grounded in the generating scripts where available and in filenames otherwise — regenerate
the scripts in `dev/benchmarks/` for authoritative captions.*
