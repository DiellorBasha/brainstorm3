# Dirac-Connectome: a vector connectome operator — design

**Date:** 2026-07-02
**Status:** design (approved sketch + 2 decisions; spec for review)
**Part of:** the Dynamics-portal / eigen-operator program. Builds on the Dirac atoms, the mode-kernel
atoms, and the scalar LB-Connectome operator.

---

## 1. Motivation

We have four operator cells; three are built:

| | over geometry (surface) | over fibers (connectome) |
|---|---|---|
| scalar | Laplace-Beltrami | LB-Connectome |
| 3-vector | **Dirac** | **Dirac-Connectome (this)** |

The scalar pair differs only in the graph the Laplacian is built on (mesh cotangent vs the fiber
connectome, `K_LBO + γ(D−W)`). The **Dirac-Connectome** is that substitution promoted to 3-D vectors: a
current field that spreads **along white-matter fibers** rather than over the cortical surface, with
orientation preserved.

Ambient-flat is the key enabler: the surface Dirac's transport is already ambient-flat (direction ~constant
in ℝ³), and its quaternion/curvature coupling comes from the *surface's* differential geometry — which the
connectome graph does not have. So the natural (and well-defined) vector-connectome operator is the
**ambient, component-wise** one: apply the scalar connectome operator independently to each ambient
component. Ambient coordinates are also the *right* choice here because two fiber-connected vertices share
no tangent frame — so this is a Dirac-style (ambient), not a connection-style (tangent), object.

## 2. Decisions (approved)

- **Base operator = LB-Connectome** (the geometry+fiber blend `K_LBO + γ(D−W)` behind the scalar
  "Connectomic" atoms), lifted to 3-D. (Not the pure normalized graph Laplacian.)
- **Scope = full source-space:** atoms (impulse + Apply), the frame scalogram, and the differential
  (div/curl/Helmholtz) stack — all on the fiber-spread vector field. **No sensor forward** (fibers are not
  a forward model). Source-space analysis basis, like scalar LB-Connectome.

## 3. Architecture

### 3.1 The lift (construction) — reuse the scalar eigenbasis, don't re-solve

Take the scalar **LB-Connectome** eigenbasis (single-block whole-brain): `Φ_s [nV × K]`, `Λ_s [K]`,
mass `M_s [nV × nV]` (verified: LB-Connectome is single-block, `GlobalVertices{1}=all`, `{2}` empty).
**Lift** each scalar mode `φ_k` to three quaternion modes in the imaginary slots (`w=0`), in the
interleaved `[w,x,y,z]` layout (rows `(v−1)*4 + {1,2,3,4}`):

```
ψ_{k,x} : imag-x slot = φ_k   ψ_{k,y} : imag-y = φ_k   ψ_{k,z} : imag-z = φ_k
```

- Eigenbasis `Φ_q [4nV × 3K]`, eigenvalues `Λ_q = [Λ_s Λ_s Λ_s]` (each `λ_k` 3-fold), quaternion mass
  `M_q = kron(M_s, I₄)` (verified: the Dirac mass is exactly `kron(M_scalar, I₄)`).
- Orthonormal by construction: `⟨ψ_{k,d}, ψ_{k',d'}⟩_{M_q} = δ_{dd'} ⟨φ_k,φ_{k'}⟩_{M_s} = δ_{dd'}δ_{kk'}`.

Filtering a 3-vector source on `Φ_q` with `g(λ)` = applying the scalar connectome filter to each ambient
component (direction-preserved). This *is* `L_conn ⊗ I₃`, realized by a lift rather than a 3×-larger
eigensolve.

### 3.2 New `Dirac-Connectome` variant (operator + eigen + registry)

- **`tess_operators`** — add a `Dirac-Connectome` variant: reuse `local_build_connectome_operator`'s
  LB-Connectome scalar operator (`K = K_LBO + γ(D−W)`, mass `M_s`), then present it as a vector operator
  by supplying the quaternion mass `M_q = kron(M_s, I₄)` and the field metadata. Single-block whole-brain.
- **`tess_eigen`** — add a `Dirac-Connectome` path: compute (or load) the scalar LB-Connectome eigenbasis
  `{Φ_s, Λ_s}` (reuse the existing solve), then lift → store `Φ_q [4nV × 3K]`, `Λ_q`, `GlobalVertices`,
  `Provenance` (Backend='lift', BaseVariant='LB-Connectome'). Find-or-create like every other variant.
- **`bst_nxr_registry`** — connectome operators are non-nxr (no `Registry.Primary`), so `field_type` comes
  from the **layout fallback** in `bst_eigenfilter('Fiber')` (`size(Φ,1)/nV = 4 → quaternion`), which
  already handles a missing registry. (Optionally add a map entry for documentation.)

### 3.3 Maximal reuse via the quaternion layout

Because `Dirac-Connectome` uses the `[4n]` quaternion layout (`w=0`), `Fiber` → C=4 (`quaternion`), and
`manifold_quat_imag`, `i_dirac_recon`, `i_atom_realise_core`, `bst_eigenfilter('Atom')`, and the quaternion
branch of `bst_eigenwavelet('Scalogram')` all work **unchanged**. The impulse (Design) preview works via
the generic realiser as soon as the eigenbasis exists.

### 3.4 Apply as a vector operator (no mode-kernel, no sensor)

The just-shipped free-projection uses the *inverse's* Dirac basis (`ImagingKernelMode`); the connectome
basis ≠ the inverse's basis, so `Dirac-Connectome` uses the **reconstruct-then-project** vector path:
reconstruct `J` from the source (`GetResultsValues`) → embed → project onto the connectome-Dirac `ax` →
filter `g(λ)` → reconstruct (`i_dirac_recon`) → cortex magnitude + quivers. In `i_atom_apply`, admit
`Dirac-Connectome` into the vector branch alongside `Dirac`, but:
- do **not** take the mode-kernel sub-branch (that stays gated on `variant=='Dirac' && i_is_dirac_dspm`);
- `i_dirac_leadfield` returns `[]` for a non-`Dirac` vector variant → `Dfilt=[]` → **no sensor overlay**
  (cortex only). (The existing `i_dirac_forward`/`i_dirac_forward_modes` already return `[]` for empty
  `Leig`.)

### 3.5 Panel wiring

- Operator dropdown / `opVariants`: add a 5th entry `'Dirac-Connectome'` (label e.g. "Dirac (connectome)").
- `i_gate_mask`: extend to 5 operators. Unconstrained (`nComponents==3`) admits `Dirac-Connectome`
  (vector); constrained (scalar) disables it (`[LBO, LB-Connectome, Connection-L, Dirac, Dirac-Connectome]`
  → constrained `[1 1 1 0 0]`, unconstrained `[1 1 0 1 1]`).
- `i_atom_axes`: the canonical branch already routes any non-Dirac-dSPM Variant through
  `bst_eigen('Axes', …, 'Variant', variant)`, so `Dirac-Connectome` resolves via the standard find-or-create
  once tess_eigen supports it (no dSPM/mode-kernel branch — Dirac-Connectome is never a dSPM basis).

### 3.6 Helmholtz / differential on the fiber-spread field

The scope includes running the div/curl/Helmholtz stack on the **connectome-Dirac-filtered** vector field
(fiber-mediated sources/sinks/vortices). The atom Apply produces `V3 [nV×3]` (the fiber-spread current);
feed that to `process_helmholtz('Compute', V3(:), Cov)` (the existing engine; `Cov` is the surface's
Covariant operator, geometry-based — the *differential* is still taken on the surface, of a field that was
*spread over fibers*). Surface a "Helmholtz on the filtered field" action (or extend the differential
overlay to operate on the atom-filtered `V3` rather than only the raw source). Reuse `process_helmholtz`,
`bst_divergence`/`bst_curl`, `view_helmholtz` unchanged.

## 4. Scope & edge cases

- **In:** `Dirac-Connectome` operator + eigenbasis (lifted LB-Connectome); impulse + Apply (cortex
  magnitude + quivers, fiber-spread); frame scalogram; Helmholtz/differential on the filtered field.
- **Out:** sensor forward (no fiber leadfield); the pure-fiber (Connectome-Laplacian) vector variant
  (deferred); fiber-orientation-aware coupling (a *non-trivial* connection along tracts using tractography
  endpoint frames — a research increment; v1 is the ambient/trivial connection).
- **Edge:** whole-brain single-block (the connectome couples hemispheres) — the atom machinery's hemi
  loop handles a single block (`{2}` empty); a fiber-spread impulse legitimately crosses hemispheres.
  Mode count is `3K` (3× the scalar K) — keep the default K modest (e.g. 200 scalar → 600 vector) or expose
  it. Requires the surface to have a connectome (own or HCP-registered via `tess_connectome`); if absent,
  the variant is unavailable (the operator dropdown item disables / errors cleanly).

## 5. Testing

- **Headless:**
  - **Lift orthonormality:** `Φ_qᵀ M_q Φ_q ≈ I_{3K}`; `Λ_q` = `Λ_s` repeated 3×.
  - **Component-wise equivalence:** filtering an embedded 3-vector source on `Φ_q` with `g(λ)` equals
    applying the scalar LB-Connectome filter to each ambient component independently (to `~1e-10`).
  - **Fiber vs geodesic:** a Dirac-Connectome diffusion impulse at a seed spreads to fiber-connected
    vertices that are geodesically distant (its peak-spread set differs from the surface-Dirac impulse's).
  - `Fiber(ax)` → C=4 / `quaternion` for the lifted basis (via layout fallback).
  - Scalogram runs on the connectome-Dirac coefficients (per-vertex `W`, full `√λ` range).
- **Live (controller, MCP):** operator dropdown shows "Dirac (connectome)" for an unconstrained source;
  Create atom → fiber-spread cortex magnitude + quivers (visibly different from surface-Dirac spread);
  Apply on the real source; scalogram; run Helmholtz on the filtered field → fiber-mediated div/curl map.
  Screenshots.

## 6. Risks / notes

- **Whole-brain single-block vs the hemi-loop:** the Dirac machinery loops `for h=1:numel(ax.Phi)`;
  confirm single-block (`{2}` empty) is handled everywhere the Dirac-Connectome `ax` flows (it is for the
  scalar LB-Connectome path; verify for the quaternion decode/recon).
- **`3K` mode count + memory:** `Φ_q` is `[4nV × 3K]` (~4× rows, 3× cols vs the scalar). Keep K modest.
- **Helmholtz-on-filtered wiring** is the least-templated piece (§3.6) — the differential currently runs on
  the raw per-frame source; operating it on the atom-filtered `V3` is a new path. Keep it a distinct action
  so the raw-source differential is unaffected.
- **Registry/field_type:** rely on the layout fallback (non-nxr operator); do not fabricate an nxr id.
- Keep all existing operators (scalar, Dirac, mode-kernel) byte-unchanged; `Dirac-Connectome` is additive.
