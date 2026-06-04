# Design: Connection Laplacian in Brainstorm (`tess_connection_laplacian`)

- **Date:** 2026-06-04
- **Author:** Diellor Basha (with Claude)
- **Status:** Design approved — ready for implementation plan
- **Branch:** feature/connection-laplacian (off `development`)

## 1. North star (context, not this milestone)

The long-term objective is a **phase-bearing spectral basis** on the cortical
surface: the eigenmodes of the **connection Laplacian** carry phase (unlike the
scalar Laplace–Beltrami eigenmodes), so the argument of an eigenmode's complex
tangent field can mark angular location. The premise — a hypothesis to be
tested, not assumed — is that the first non-trivial eigenmode of the connection
Laplacian, read out against a globally consistent reference frame, completes one
cycle over the hemisphere, analogous to how the first DFT basis function's phase
marks time. See `dev/vector_eigenmode_meg_proposal.md` and the project memory
`connection-laplacian-phase`.

This document covers **only Milestone 1**: integrating the connection Laplacian
**operator** into Brainstorm. Eigenmodes, phase readout, registration alignment,
and storage are explicitly deferred to later milestones.

## 2. Key conceptual decisions (from brainstorming)

These were settled before design and constrain the implementation:

1. **The operator is the canonical Levi-Civita vertex connection Laplacian**
   (what nxr/geometry-central assembles from `(V, F)` alone). We do **not** build
   it from the trivial connection. The trivial connection collapses the operator
   toward the real scalar Laplacian in its trivializing gauge (real eigenvectors,
   phase only 0/π except at singularities), destroying the distributed phase
   winding that is the point. Levi-Civita is gauge-covariant: its **eigenvalues
   are intrinsic**, and eigenvector **phase is gauge-dependent**.

2. **The operator depends on the mesh alone**, not on any tangent frame.
   nxr's per-vertex frames are local/arbitrary (anchored to each vertex's first
   outgoing halfedge) and are *not* globally consistent — they cannot serve as a
   reference field. The smooth reference frame (Brainstorm's trivial-connection
   `TangentFrame`) enters **only at phase readout**, as a gauge transform
   `z_i → e^{-iθ_i} z_i`. This is why the operator's signature is a pure `(V, F)`
   function: it is the mathematically honest interface.

3. **Cross-subject comparison is a readout-layer concern (next phase).** It is
   achieved by expressing the intrinsic eigenmode against the FreeSurfer-pole
   reference frame. Deferred subtleties: vertex correspondence (reg sphere),
   global phase ambiguity per eigenmode (gauge-fix at a pole), sign/ordering/
   degeneracy, and the per-face → per-vertex reference-frame transfer.

4. **Dedicated function, not a branch of `tess_laplacian`.** The scalar cotan
   Laplacian returns real PSD `[L, M]`; the connection Laplacian returns a complex
   Hermitian operator with different metadata, has no MATLAB fallback, and shares
   no computational body. A single function with output meaning toggled by a flag
   would have an opaque contract. `tess_tangents` (also nxr-only,
   trivial-connection-based) is the sibling precedent. The only real duplication —
   the nxr-loaded guard — is factored into a shared helper.

5. **Complex format.** Per user preference, the operator is returned as a native
   MATLAB complex sparse matrix (MATLAB has first-class complex support); the
   real2N packing is not exposed.

## 3. nxr-compute capabilities relied upon (verified)

- `nxr.manifold.operator.connectionLaplacian(mctx, opts)` →
  `assembleConnectionLaplacian`. Options: `domain` (`vertex`/`face`/`edge`),
  `nSym` (default 1), `regularization` (default 1e-8), `format`
  (`real2N`/`complex`). C++ mirrors geometry-central's
  `computeVertexConnectionLaplacian` (Levi-Civita transport raised to `nSym`).
- For `format='complex'`, the MEX returns the struct fields `K_real` and `K_imag`
  as two **real** sparse matrices (`marshal.h:327` `connectionLaplacianToStruct`).
  MATLAB recombines: `K = K_real + 1i*K_imag`.
- `domain` and `format` are returned as **numeric enum codes** (doubles), not
  strings: domain Vertex=0/Face=1/Edge=2; format Real2N=0/Complex=1.
- `nxr.manifold.context(V, F)` exposes `ctx.M` = geometry-central's
  `vertexLumpedMassMatrix` (diagonal, `area/3` per vertex) = the lumped vertex
  mass. Faces are passed **1-based**; the MEX marshalling subtracts 1 (already
  proven to machine precision by the `tess_laplacian` parity test).
- The connection Laplacian assembly requires a clean 2-manifold; nxr throws
  `nxr:nonManifold` otherwise.

## 4. Interface

```matlab
[K, M, Info] = tess_connection_laplacian(Vertices, Faces, varargin)
```

**Inputs**
- `Vertices` : `[nV x 3]` double, vertex positions.
- `Faces`    : `[nF x 3]` triangle indices, **1-based** (Brainstorm convention).

**Options (name-value)**

| Option | Default | Meaning |
|---|---|---|
| `nSym` | `1` | Bundle symmetry. `1` = true vector field (carries phase); `2` = line field; `4` = cross field. |
| `Domain` | `'vertex'` | `'vertex'` (target) \| `'face'` \| `'edge'`. |
| `Regularization` | `1e-8` | ε·I diagonal added by nxr for strict positive-definiteness. |
| `CheckManifold` | `false` | Optional `tess_manifold` pre-check for a friendlier error message. |

**Outputs**
- `K` : `[N x N]` **complex Hermitian** sparse connection Laplacian
  (`N = nV` for the vertex domain), assembled as `K_real + 1i*K_imag`.
- `M` : `[N x N]` real diagonal **lumped vertex mass** (`area/3` per vertex), the
  companion mass for the eventual generalized eigenproblem `K φ = λ M φ`. In
  complex format no block duplication is required. From `ctx.M`.
- `Info` : struct with fields `nSym`, `Domain` (string), `Regularization`,
  `baseDim` (= `N`), `Format` (= `'complex'`), `Backend` (= `'nxr'`).

## 5. Backend & error handling

nxr-only, mirroring `tess_tangents` (no MATLAB fallback — the operator requires
geometry-central's Levi-Civita transport vectors).

- The nxr-loaded check is factored out of `tess_laplacian` into a shared helper
  (e.g. `nxr_is_loaded`) that both functions call. This is the only
  `tess_laplacian` change — a behaviour-preserving refactor.
- If nxr is not loaded: raise a clear error pointing to
  `bst_plugin('Install','nxr-compute')` (no silent degradation).
- If `CheckManifold` is true: run `tess_manifold(..., 'Repair', 0)` first and
  raise a friendly error on failure; otherwise let nxr surface `nxr:nonManifold`.

## 6. Correctness checks (verified, not assumed)

- **Hermitian:** `max(max(abs(K - K')))` ≈ 0 (conjugate transpose).
- **Positive (semi-)definite:** smallest eigenvalues ≥ 0 (≈ `Regularization`).
- **Faces convention:** 1-based passed straight through (proven by existing
  parity test); no off-by-one / degenerate operator.
- **Enum mapping:** numeric `domain`/`format` codes mapped back to strings in
  `Info`.

## 7. Tests

Function-style scripts under `dev/tests/` (repo idiom; run via the MATLAB MCP,
not `runtests`). On `tess_sphere(642)` unless noted:

1. **Smoke** (`test_connection_laplacian_smoke.m`): assemble `K`; assert
   `N x N`, complex, sparse, Hermitian; `M` diagonal and positive; `Info` fields
   correct (including string `Domain`/`Format`).
2. **Spectral sanity** (`test_connection_laplacian_spectrum.m`): solve the
   smallest `k` modes via `nxr_compute('solve', K, M, k)`; assert eigenvalues
   real and ≥ 0; the lowest non-trivial mode is a smooth field (qualitative —
   full spectral validation belongs to the eigenmode milestone).
3. **nSym variants** (`test_connection_laplacian_nsym.m`): `nSym = 1, 2, 4` each
   assemble and remain Hermitian.
4. **Backend guard** (`test_connection_laplacian_guard.m`): with nxr unloaded,
   an informative error is raised.

## 8. Files

| Path | Change |
|---|---|
| `toolbox/anatomy/tess_connection_laplacian.m` | New: the operator function. |
| `toolbox/anatomy/tess_laplacian.m` | Refactor: extract the nxr-loaded guard into the shared helper. |
| (shared helper) | New small function for the nxr-loaded check, called by both. |
| `dev/tests/test_connection_laplacian_*.m` | New: the four tests above. |

## 9. Out of scope (later milestones)

- Eigenmode computation and the complex spectral basis.
- Phase readout and the smooth reference frame (gauge transform).
- Registration alignment and cross-subject comparison.
- Per-vertex trivial-connection `TangentFrame` (face → vertex transfer).
- Persistence on the surface file, GUI/process wrappers.
- Vector-valued source-map analysis built on the operator.
