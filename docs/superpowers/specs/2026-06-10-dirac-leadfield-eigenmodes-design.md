# Dirac Eigenbasis + Unconstrained-Leadfield Transform — Design

**Date:** 2026-06-10
**Status:** Design (pending review)
**Repos:** `brainstorm3` (MATLAB). Depends on the `nxr-compute` `operators(...,'dirac',τ)` operator (already implemented on `main`).

## Goal

Give the **unconstrained (vector) leadfield** a principled eigenmode transform by
expanding it in the **Dirac eigenbasis** instead of the scalar Laplace–Beltrami
basis. The scalar LBO basis (`Φ [nV×K]`) can only expand scalar fields; the
cortical leadfield is a per-vertex 3-vector field, so it needs a vector/spinor
eigenbasis. The relative-Dirac family `L(τ) = (1−τ)(cotanL⊗I₄) + τ·E` provides a
curvature-aware, quaternionic basis whose imaginary part lives in ambient ℝ³ —
exactly the space of leadfield vectors.

## Key insight — reuse the existing eigenmode inverse

Brainstorm already has the scalar mode-space pipeline. This work swaps only the
**basis** and the **leadfield composition**; the downstream inverse is reused
unchanged:

| Scalar (exists) | Dirac (new) |
|---|---|
| `tess_eigenmodes` → `Φ [nV×K]` | **Phase A:** Dirac eigenbasis `Φ_D [4V×K]`/hemi |
| `bst_eigenmode_leadfield` → `L̃ = L·Φ` | **Phase B:** quaternionic compose `L̃ = Φ_Dᵀ B Ψ` |
| `bst_inverse_eigenmodes` (+ `bst_eigenmode_prior`) | **reused unchanged** |

## Decisions (resolved during brainstorming)

- **Scope:** both pieces, decomposed — Phase A (eigenbasis) then Phase B
  (leadfield transform). Each becomes its own implementation plan.
- **Eigensolve location:** MATLAB `eigs(L, B, K, 'smallestabs')`, matching
  `tess_eigenmodes` (operator matrices from nxr; nxr's own `solve` is unused in
  Brainstorm). No new nxr command.
- **Embedding:** ambient-direct — each leadfield 3-vector `g=(gx,gy,gz)` (world
  coords) → pure-imaginary quaternion `ψ=[0,gx,gy,gz]` (`w=0`). Matches the
  relative-Dirac operator's own convention (built from ambient vertex normals).
- **Default `τ = 0.5`** (elliptic/PSD, curvature-aware), user-overridable.
- **Default `K = 400`** (≈100 quaternionic 4-fold multiplets/hemisphere),
  user-overridable.
- **Storage:** the Dirac eigenbasis is stored on the surface `TessMat` as a
  per-hemisphere field (analogous to the planned `Gauge.Eigen`), not a separate
  eigenmodes file.

## Prerequisite

The installed MEX (`~/.brainstorm/plugins/nxr-compute/...-r2023b/nxr_compute.mexmaca64`)
must expose `nxr_compute('operators', h, 'dirac', τ)`. The dirac operator is on
nxr `main` but the currently-installed binary predates it (it was built for the
facets work). **Verify, and if absent build + install** (`scripts/build.sh
Release` → copy into the plugin dir, backing up the current binary) before
Phase A. This is a build step, not a design fork.

---

## Phase A — Dirac eigenbasis (`tess_dirac_eigenmodes`)

New writer, sibling to `tess_eigenmodes`, self-contained (Brainstorm I/O; no
shared helper):

1. `TessMat = in_tess_bst(SurfaceFile)`; `TessFile = file_fullpath(...)`.
2. Cache-return: if `TessMat.DiracEigen` present (`1×2`, matching `Tau`/`nModes`)
   and not `ForceRecompute`, return it.
3. Require nxr-compute; split with `tess_hemisplit` (atlas L/R, **never**
   `conncomp`; require a Structures atlas with L/R labels, else error).
4. **Per hemisphere** (build the local submesh, `nxr_compute('create',Vloc,Floc)`):
   - `L = nxr_compute('operators', h, 'dirac', Tau)` → `[4Vₕ×4Vₕ]` real symmetric.
   - `Mg = nxr_compute('operators', h, 'mass', 'galerkin')` → `[Vₕ×Vₕ]`;
     `B = kron(Mg, speye(4))` → `[4Vₕ×4Vₕ]`.
   - `[Phi, D] = eigs(L, B, K, 'smallestabs')`; `lam = diag(D)`; sort ascending.
   - **B-orthonormalize** so `PhiᵀBPhi = I` (`normalizeEigenmodes`-style:
     `Phi = Phi / chol-or-sqrtm of (PhiᵀBPhi)`; within a 4-fold multiplet, any
     B-orthonormal basis of the eigenspace is acceptable — projection energy is
     gauge-invariant, per the nxr Dirac spec).
   - `nxr_compute('destroy', h)`.
5. **Store** `TessMat.DiracEigen` as a `1×2` per-hemisphere struct array
   ((1)=L,(2)=R), each element:
   - `Vectors [4Vₕ×K]` (vertex-interleaved `4v+c`, order `[w,x,y,z]`),
   - `Values [K×1]` (ascending), `nModes`, `Order`,
   - `Tau`, `GlobalVertices`, `Hemisphere`,
   - `Provenance` (`Backend='nxr'`, `NxrVersion`, `Tau`, `K`, `ComputeDate`).
   `bst_history` + `bst_save(TessFile, ..., 'v7')` unless `NoSave`.

**Options:** `Tau` (0.5), `K` (400), `NoSave`, `ForceRecompute`.

**Phase A tests** (real 20484 cortex): `DiracEigen` is `1×2`; `Vectors` is
`[4Vₕ×K]`; `Values` ascending and `≥ −1e-9`; `ΦᵀBΦ ≈ I` (`<1e-9`); eigenvalues
cluster in 4-fold multiplets; `Tau`/provenance present; cache-return + `NoSave`
behave; save isolation via `onCleanup`.

---

## Phase B — unconstrained-leadfield transform (`bst_dirac_eigenmode_leadfield`)

Composer, sibling to `bst_eigenmode_leadfield`. Input: an **unconstrained**
surface head model (`Gain [nCh × 3·nSrc]`, `nComponents=3`) whose source space is
the cortex carrying `DiracEigen`.

1. Load the head model `Gain` and the cortex `TessMat.DiracEigen` (compute it via
   Phase A if absent).
2. Guard: head model is surface + unconstrained (3 components/source); source
   count matches the cortex vertex count.
3. **Per hemisphere** (local vertices `vH = DiracEigen(hh).GlobalVertices`):
   - Build `Ψₕ [4Vₕ × nCh]`: for local vertex `vloc` ↔ global source
     `s = vH(vloc)`, with gain columns `Gs = Gain(:, 3*(s-1)+(1:3))` (`[nCh×3]`):
     `Ψₕ(4(vloc−1)+1, :) = 0` (w); rows `+2,+3,+4 = Gs(:,1)', Gs(:,2)', Gs(:,3)'`.
   - `Bₕ = kron(Mgₕ, speye(4))` (re-pull `Mg` from nxr, or recompute via Phase A
     provenance); `L̃ₕ = Ψₕᵀ · Bₕ · DiracEigen(hh).Vectors` → `[nCh × K]`
     (channels × modes — the head-model `Gain` convention; this is the transpose
     of the projection coefficients `Φ_Dᵀ B Ψ`).
4. **Stack** hemispheres horizontally: `Gain_tilde = [L̃_L, L̃_R]` (`[nCh × 2K]`),
   `Eigenvalues = [λ_L; λ_R]` (`[2K × 1]`), `nModes = 2K`, plus a per-mode
   hemisphere/`GlobalVertices` map for reconstruction.
5. Write a **composed eigenmode head model** with `Gain = Gain_tilde`
   (`[nCh × 2K]`), `.Eigenvalues`, `.nModes` — the exact struct
   `bst_inverse_eigenmodes` consumes (`L_all = HM.Gain` is `[nAllCh × K]`).

**Phase B tests:** `Gain_tilde` is `[nCh × 2K]`; projection identity
`L̃ₕ == Ψₕᵀ Bₕ Vectorsₕ`; in-span round-trip (a field built as `Φ_D·c` recovers `c`
within tol via `c = Vectorsᵀ B Ψ`); `bst_inverse_eigenmodes('SolvePure',
Gain_tilde, Eigenvalues, …)` returns a `[2K × nCh]` kernel without error.

---

## Reconstruction / inverse reuse

`bst_inverse_eigenmodes` (unchanged) returns a mode-space kernel
`M̃ [2K × nCh]`. Per-vertex source currents are recovered per hemisphere as
`J = imag( DiracEigen(hh).Vectors · (M̃_block · data) )` → a **3-vector
(unconstrained) estimate per vertex** directly (the `w` part is dropped). This
maps to `ImageGridAmp`/`nComponents=3` naturally — a property the scalar pipeline
lacks. Wiring the reconstruction into the results path is downstream of this
design (the mode-space inverse + spectral prior `bst_eigenmode_prior` are
reused as-is).

## Edge cases / notes

- **4-fold multiplets:** not canonicalized; only the B-orthonormal eigenspace is
  used. Truncating `K` mid-multiplet is acceptable (round `K` to a multiple of 4
  if a clean cutoff is desired — default 400 already is).
- **τ provenance:** `DiracEigen.Provenance.Tau` must match the value Phase B
  assumes; Phase B reads `Tau` from the stored basis (does not re-pick it).
- **Hemisphere split:** atlas-only via `tess_hemisplit`; error without a
  Structures L/R atlas (never the geometric fallback).
- **Mass `B` in Phase B:** prefer recomputing `Mg` from nxr for the same submesh
  rather than storing the (large) `B` on disk; `Vectors` already carries the
  basis.

## Out of scope

- Changes to `bst_inverse_eigenmodes` / `bst_eigenmode_prior` (reused as-is).
- A new nxr eigensolve command (MATLAB `eigs` per the existing convention).
- Cross-surface quaternion-gauge canonicalization (Liu §4.3) — single fixed
  cortex only.
- The results/visualization wiring of the unconstrained reconstruction (separate).
- Deleting the older scalar eigenmode path (coexists).

## Decomposition

- **Plan 1 (Phase A):** `tess_dirac_eigenmodes` + tests (+ MEX-dirac
  build/install prerequisite check).
- **Plan 2 (Phase B):** `bst_dirac_eigenmode_leadfield` + tests, validated end to
  end through `bst_inverse_eigenmodes('SolvePure', …)`.
