# Design: Connection Eigenmodes as a tess-file axis (`ConnEigenmodes`)

- **Date:** 2026-06-04
- **Author:** Diellor Basha (with Claude)
- **Status:** Design approved — ready for implementation plan
- **Branch:** feature/conn-eigenmodes (off `development`)
- **Builds on:** `tess_connection_laplacian` (Milestone 1, merged) — see
  `dev/connection_laplacian_integration.md`.

## 1. Goal

Add `ConnEigenmodes`, a **canonical, intrinsic spectral axis** on the Brainstorm
surface file holding the **eigenmodes of the connection Laplacian** — a
vector-field (n-RoSy, `nSym=1`) basis that carries phase. It is a sibling to the
existing scalar `TessMat.Eigenmodes` axis and is stored/loaded/ensured through a
parallel quartet of functions. The two axes are independent and serve different
purposes: `Eigenmodes` for scalar functions, `ConnEigenmodes` for tangent vector
fields.

This is **Milestone 2**. Phase readout, gauge-fixing, and registration alignment
remain a later milestone (§9). See project memory `connection-laplacian-phase`.

## 2. Decisions (from brainstorming)

1. **Naming:** field `TessMat.ConnEigenmodes`; functions `tess_conn_eigenmodes`,
   `in_tess_conn_eigenmodes`, `out_tess_conn_eigenmodes`,
   `bst_conn_eigenmodes_ensure`. Parallel to the scalar quartet.
2. **Default count matches the scalar axis.** `bst_conn_eigenmodes_ensure`
   derives the per-component count from the surface's scalar `Eigenmodes`
   (`round(nModes / nComponents)`), ensuring the scalar axis exists first
   (`bst_eigenmodes_ensure`). An explicit `nModesPerHemi` argument overrides.
3. **Store the operator `K`** in the struct (`ConnLaplacian`), mirroring how the
   scalar struct stores `Laplacian` — reusable downstream (readout, filtering).
4. **`nSym = 1`** (true vector field — the phase-bearing case).
5. **Per-connected-component solve**, mirroring `tess_eigenmodes`. The connection
   Laplacian is block-diagonal across hemispheres; each block is solved
   independently and tagged with `Component`/`CompRank`.
6. **No DC removal.** The connection Laplacian has no zero mode (no globally
   consistent parallel vector field on a curved closed surface; smallest
   eigenvalue is bounded away from 0), so `nRemoved = 0` always.
7. **Raw storage, no gauge-fixing.** Eigenvectors are stored exactly as the
   solver returns them, in nxr's intrinsic per-vertex frames. Phase readout and
   registration are deferred (§9).
8. **No repair** (mirror the scalar `ensure` no-repair policy; non-manifold →
   error). nxr-only (the operator already requires it).

## 3. Architecture — a parallel quartet

| Function | Mirrors | Responsibility |
|---|---|---|
| `tess_conn_eigenmodes(V, F, ...)` | `tess_eigenmodes` | Assemble `[K, M]` via `tess_connection_laplacian(V,F,'nSym',1)`; solve per component with MATLAB `eigs`; package the struct. |
| `out_tess_conn_eigenmodes(SurfaceFile, ConnEig, V, F[, isInteractive])` | `out_tess_eigenmodes` | Write `TessMat.ConnEigenmodes` (vectors complex `single`; operators `double`); add `ComputeDate` + history; single `bst_save`. |
| `in_tess_conn_eigenmodes(SurfaceFile)` | `in_tess_eigenmodes` | Read it back via `in_tess_bst(SurfaceFile, 0)`; cast vectors → `double` complex; backfill `Order`/`Component`/`CompRank`/`nComponents`; return `[ConnEig, isComputed]`. |
| `bst_conn_eigenmodes_ensure(SurfaceFile[, nModesPerHemi])` | `bst_eigenmodes_ensure` | Reuse if present; else derive the count from the scalar axis (ensuring it first), compute, store. |

The scalar functions and `tess_connection_laplacian` are **not modified**.

### 3.1 Solver
For each connected component with vertex index set `idx`:
```matlab
Kc = K(idx, idx);   Mc = M(idx, idx);       % exact block (operator is block-diagonal)
[Phi_c, Lam_c] = eigs(Kc, Mc, k, 'smallestabs');   % complex Hermitian generalized
```
`'smallestabs'` targets the smoothest (lowest) modes — the spectrum is strictly
positive, so the smallest-magnitude modes are the smallest eigenvalues. nxr's own
`solve` is real-only and cannot be used (see M1 spec §3). Eigenvalues are real
(imaginary part negligible); take `real(...)`. Modes are placed back into the
whole-mesh `Vectors` at rows `idx`, columns for that component (block structure,
zeros elsewhere), tagged with `Component`/`CompRank`. `Order` is the global
ascending sort of `Values`.

### 3.2 Match-scalar count derivation (`bst_conn_eigenmodes_ensure`)
```matlab
% Reuse if present
[ConnEig, isComputed] = in_tess_conn_eigenmodes(SurfaceFile);
if isComputed && ~isempty(ConnEig), Eig = ConnEig; return; end
% Derive count from the scalar axis (compute its default if absent)
if nargin < 2 || isempty(nModesPerHemi)
    sEig = bst_eigenmodes_ensure(SurfaceFile);          % scalar axis (1000/hemi default if absent)
    nModesPerHemi = max(1, round(sEig.nModes / sEig.nComponents));
end
% Guard manifoldness (no repair), compute, store
... tess_conn_eigenmodes(V, F, 'nModes', nModesPerHemi) ... out_tess_conn_eigenmodes(...)
```

## 4. The `ConnEigenmodes` struct

| Field | Type | Notes |
|---|---|---|
| `Vectors` | `[nV × nModes]` **complex single** | block-structured; each mode nonzero only on its component |
| `Values` | `[nModes × 1]` real double | eigenvalues (ascending within component), > 0 |
| `nModes` | scalar | total modes across components |
| `Component` | `[nModes × 1]` | connected-component id per mode |
| `CompRank` | `[nModes × 1]` | within-component rank |
| `Order` | `[nModes × 1]` | global ascending order of `Values` |
| `nComponents` | scalar | components solved |
| `MassMatrix` | sparse double | lumped vertex mass (real diagonal), basis is M-orthonormal |
| `ConnLaplacian` | sparse **complex** double | the Hermitian operator `K` (reuse downstream) |
| `OperatorType` | char | `'Connection-LeviCivita'` |
| `nSym` | scalar | `1` |
| `Regularization` | scalar | `1e-8` |
| `Sigma` | scalar/char | eigs shift used (`'smallestabs'`) |
| `Tolerance` | scalar | eigs tolerance |
| `nRemoved` | scalar | `0` (no DC mode in the connection bundle) |
| `ComputeTime` | scalar | seconds |
| `ComputeDate` | char | set by `out_tess_conn_eigenmodes` |

### 4.1 Storage notes
- **Complex `single` dense vectors** save correctly under `bst_save(..., 'v7')`.
- **Sparse must be `double`** — MATLAB has no single sparse — so `MassMatrix` and
  the complex `ConnLaplacian` stay double (negligible vs `Vectors`).
- `in_tess_conn_eigenmodes` casts `Vectors` back to `double` complex on load
  (mirrors the scalar loader casting single→double).

## 5. Correctness / sanity (verified by tests, not assumed)
- `Vectors` is `nV × nModes`, complex, and block-structured (a mode's nonzeros
  lie only on its component's vertices).
- `Values` real and > 0 (no zero mode).
- M-orthonormality within a component: `Phi_c' * Mc * Phi_c ≈ I` (Hermitian inner
  product; `'` is conjugate transpose).
- `Order` indexes `Values` in ascending order.
- I/O round-trip preserves vectors (single→double complex) and the complex
  operator.

## 6. Tests

Function-style scripts under `dev/tests/` (run via the MATLAB MCP). They resolve a
real 20484-vertex cortex from the loaded protocol by vertex count (same resolver
idiom as the M1 tests; on this machine: TutorialAuditory / Subject01), and SKIP
cleanly if none is present. To stay fast, the compute/round-trip tests use a small
`nModes` (e.g. 20), not the full default.

1. **Compute** (`test_conn_eigenmodes_compute.m`): `tess_conn_eigenmodes(V,F,'nModes',20)`
   → assert struct fields/types; `Vectors` complex `nV×nModes`; `Values` real & >0;
   per-component block structure; `Order` sorted ascending; M-orthonormality per
   component.
2. **I/O round-trip** (`test_conn_eigenmodes_roundtrip.m`): compute on a temp copy
   of the surface, `out_tess_conn_eigenmodes`, then `in_tess_conn_eigenmodes`;
   assert `isComputed` toggles false→true, `Vectors` returns as `double` complex
   and matches (single round-trip tolerance), `ConnLaplacian` preserved. Work on a
   temp copy so the DB surface is not mutated (mirror `test_tess_tangents.m`).
3. **Ensure** (`test_conn_eigenmodes_ensure.m`): on a temp copy, with an explicit
   small `nModesPerHemi`, `bst_conn_eigenmodes_ensure` computes + stores, and a
   second call reuses (idempotent, no recompute). Separately assert the
   match-scalar derivation: with the scalar axis present at a known per-component
   count, the no-arg ensure requests `round(nModes/nComponents)` connection modes.

## 7. Files

| Path | Change |
|---|---|
| `toolbox/anatomy/tess_conn_eigenmodes.m` | New: compute the struct (per-component eigs). |
| `toolbox/io/out_tess_conn_eigenmodes.m` | New: write `TessMat.ConnEigenmodes`. |
| `toolbox/io/in_tess_conn_eigenmodes.m` | New: read it back, backfill metadata. |
| `toolbox/anatomy/bst_conn_eigenmodes_ensure.m` | New: reuse/compute-default (match scalar). |
| `dev/tests/test_conn_eigenmodes_*.m` | New: the three tests above. |

No changes to `tess_eigenmodes` / `in_/out_tess_eigenmodes` / `bst_eigenmodes_ensure`
/ `tess_connection_laplacian`.

## 8. Out of scope (later milestones)
- Phase readout (gauge transform into a smooth global frame) and gauge-fixing the
  per-eigenmode global phase.
- The per-vertex trivial-connection reference frame (face→vertex transfer) and
  registration alignment for cross-subject phase comparison.
- GUI/process wrappers, viewers.
- Using the basis for vector-valued source-map analysis.
