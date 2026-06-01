# Design: `tess_tangents.m` — globally consistent tangent frame field

- **Date:** 2026-06-01
- **Author:** Diellor Basha (with Claude)
- **Status:** Design — pending review
- **Branch:** feature/tess-tangents (off `development`)

## 1. Goal

Compute a **globally consistent tangent frame field** on a cortical surface using
nxr-compute's trivial connection, with singularities pinned at the FreeSurfer
registration-sphere poles (north + south) of each hemisphere, and store it on the
surface for later vector-field analysis of source maps. New function
`tess_tangents.m`, in the spirit of `tess_normals.m`.

A trivial connection produces a direction field that is smooth (parallel) everywhere
except at prescribed singularities, whose indices must sum to the surface's Euler
characteristic (Gauss-Bonnet). Pinning the two singularities at the registration-sphere
poles yields a frame field with a consistent, anatomically reproducible structure
across subjects (a "globe-like" frame aligned to the FreeSurfer spherical coordinate
system).

## 2. Scope

**In scope**
- `toolbox/anatomy/tess_tangents.m`: compute a **per-face** tangent frame `(U, V)` and
  store it on the surface file; return it.
- **Per-hemisphere** solve (see §3): each hemisphere is an independent genus-0 sphere.
- Singularities at the FreeSurfer registration-sphere poles (max/min z per hemisphere).
- nxr-only (SPM-style guarded install, no MATLAB fallback).
- An integration test on a real FreeSurfer-registered cortex + a missing-registration
  guard test.

**Out of scope (deferred)**
- Per-**vertex** frame. nxr's trivial connection is per-face; transferring it to
  vertices (via the Hodge star / DEC machinery) is future work. The storage format
  reserves a `Domain` field for this.
- Using the frame to express/analyze source-map vector fields (the downstream consumer).
- Any GUI/`process_*.m` wrapper.
- Multi-platform nxr binaries (covered by the nxr plugin's own roadmap).

## 3. Background & guiding principle

### 3.1 nxr trivial connection (per-face)
`nxr.manifold.interpolate.trivial(mctx, singVerts, singValues)`:
- `mctx = nxr.manifold.context(V, F)` (V `[nV×3]` double, F `[nF×3]`, 1-based; the mex
  marshalling converts to 0-based internally — pass 1-based faces unchanged).
- `singVerts`: 1-based vertex indices of singularities; `singValues`: their indices,
  which must sum to the Euler characteristic χ (Gauss-Bonnet).
- Returns a struct: `directionVectors [nF×3]` (e1), `orthogonalVectors [nF×3]` (e2,
  rotated 90° in the face tangent plane), `connections [nE×1]`, `eulerCharacteristic`,
  `gaussBonnetSatisfied` (true iff |Σindices − χ| < 1e-3).
- **Per-FACE** output. nxr-only; no MATLAB equivalent.

### 3.2 Brainstorm surface data model
- `TessMat.Reg.Sphere.Vertices [nV×3]`: per-vertex coordinates on the FreeSurfer
  registration sphere (populated by `tess_addsphere` from `?h.sphere.reg`; a unit
  sphere in FreeSurfer RAS). Poles = max-z / min-z vertices.
- `tess_hemisplit(TessMat)` → `[rH, lH, isConnected, ...]`: per-hemisphere vertex index
  lists (uses the 'Structures' atlas; falls back to a y-coordinate split).
- Per-vertex/-face arrays are cached in the surface `.mat` and saved via
  `bst_save(SurfaceFile, TessMat, 'v7')` — exactly how `VertNormals` is computed and
  cached in `in_tess_bst`.

### 3.3 Guiding principle: per-hemisphere always
In Brainstorm/FreeSurfer the left and right hemispheres are **disconnected, separate
objects**. Every nxr-compute operation that involves a solve (trivial connection,
eigenmodes, geodesics, …) must be applied **per hemisphere**: each hemisphere is its
own closed genus-0 sphere with its own Euler characteristic (χ = 2). Solving on the
concatenated 2-component cortex risks per-component Gauss-Bonnet imbalance (a global
Σindices = 4 can hide one hemisphere summing to 1 and the other to 3). This design
solves each hemisphere independently; the principle generalizes to future nxr work.

## 4. Interface

`toolbox/anatomy/tess_tangents.m`:

```matlab
% USAGE:  [U, V] = tess_tangents(SurfaceFile)              % compute, store, return
%         [U, V] = tess_tangents(SurfaceFile, 'NoSave', 1) % compute + return only
%
% INPUT:
%    - SurfaceFile : Relative path to a Brainstorm cortex surface with a FreeSurfer
%                    registration sphere (TessMat.Reg.Sphere.Vertices).
% OPTIONS:
%    - NoSave : (logical) If true, do not write TangentFrame back to the file.
%               Default: false (store).
% OUTPUT:
%    - U : [nF×3] per-face first tangent direction  (e1)
%    - V : [nF×3] per-face orthogonal tangent       (e2), with U×V aligned to the
%          face normal.
```

SurfaceFile-based (not pure `Vertices/Faces` like `tess_normals`) because the
computation inherently needs `Reg.Sphere` and the hemisphere labeling, which are
surface-file-level. The per-hemisphere nxr solve is factored into a local helper
(`solve_hemisphere`) so each hemisphere is computed and reasoned about independently.

## 5. Algorithm (per-hemisphere pipeline)

1. Resolve + load the surface: `TessMat = in_tess_bst(SurfaceFile)`.
2. **Require registration:** if `Reg.Sphere.Vertices` is absent/empty → error
   `tess_tangents:noRegSphere` ("requires a FreeSurfer registration sphere; run
   tess_addsphere / import with surface registration").
3. **Ensure nxr:** `bst_plugin('Install','nxr-compute')`; error `tess_tangents:nxrUnavailable`
   if it cannot be installed/loaded.
4. **Hemisphere split:** `[rH, lH, isConnected] = tess_hemisplit(TessMat)`. If
   `isConnected` (hemispheres connected / not cleanly separable, e.g. the y-split
   fallback) → error `tess_tangents:connectedHemispheres` (the per-hemisphere solve
   requires disconnected hemispheres).
5. Initialize `U = zeros(nF,3)`, `V = zeros(nF,3)`, `assigned = false(nF,1)`.
6. For each hemisphere `vH ∈ {lH, rH}` (skip if empty), call `solve_hemisphere`:
   - `isVH = false(nV,1); isVH(vH) = true;`
   - `fMask = all(isVH(Faces), 2);` (faces wholly within the hemisphere; valid because
     hemispheres are disconnected components). `assert(~any(assigned & fMask))`.
   - Build local submesh: `map = zeros(nV,1); map(vH) = 1:numel(vH); Fh = map(Faces(fMask,:)); Vh = Vertices(vH,:);`
   - Poles from the registration sphere restricted to the hemisphere:
     `sph = Reg.Sphere.Vertices(vH,:); [~,iN] = max(sph(:,3)); [~,iS] = min(sph(:,3));`
     (local 1-based indices into `vH`).
   - `mctx = nxr.manifold.context(Vh, Fh);`
     `r = nxr.manifold.interpolate.trivial(mctx, [iN; iS], [1; 1]);`
   - `assert(r.gaussBonnetSatisfied)` → else error `tess_tangents:gaussBonnet`.
   - Scatter back: `U(fMask,:) = r.directionVectors; V(fMask,:) = r.orthogonalVectors; assigned(fMask)=true;`
   - Record singularity metadata: global vertex indices `vH([iN iS])`, indices `[1 1]`,
     hemisphere tag.
7. **Coverage check:** `assert(all(assigned))` — every face assigned to exactly one
   hemisphere (error `tess_tangents:unassignedFaces` otherwise; would indicate a
   cross-hemisphere face or a vertex missing from both hemisphere sets).
8. Orientation: ensure `cross(U,V)` points along the face normal; flip `V` per face if
   the dot product with the face normal is negative (so the frame is right-handed and
   consistent with the surface orientation).
9. Store (unless `NoSave`) and return.

## 6. Storage format

A struct field on the surface `.mat`, saved via `bst_save(SurfaceFile, TessMat, 'v7')`:

```matlab
TessMat.TangentFrame = struct( ...
   'Domain',        'face', ...                 % 'face' now; 'vertex' later (Hodge/DEC)
   'U',             single(U), ...              % [nF×3] e1 (directionVectors)
   'V',             single(V), ...              % [nF×3] e2 (orthogonalVectors)
   'Singularities', struct('Vertices', singVerts, ...   % [4×1] global vertex indices
                           'Indices',  singIdx, ...     % [4×1] (all +1)
                           'Hemisphere', {hemiTag}), ...% {'L','L','R','R'} etc.
   'Method',        'nxr trivial-connection (FreeSurfer poles)');
% Append a History line per Brainstorm convention.
```

A struct (rather than a bare `[nF×6]`) so the deferred `Domain='vertex'` variant and
the singularity provenance slot in without a format change. `single` matches the
precision convention of other cached per-element surface arrays (`Curvature`, `SulciMap`).

## 7. Dependency & error handling

- **nxr-only:** guarded `bst_plugin('Install','nxr-compute')` at the top; no MATLAB
  fallback (the trivial connection has no built-in equivalent). Mirrors the
  `bst_normalize_mni`/SPM pattern.
- **Clean errors** (each with a `tess_tangents:*` identifier): missing `Reg.Sphere`
  (`noRegSphere`); nxr unavailable (`nxrUnavailable`); connected/not-cleanly-separable
  hemispheres (`connectedHemispheres`, from `tess_hemisplit`'s `isConnected` — the
  per-hemisphere solve requires disconnected hemispheres); a failed Gauss-Bonnet check
  (`gaussBonnet`); a face assigned to both hemispheres (`overlap`); unassigned faces
  after the per-hemisphere pass (`unassignedFaces`).

## 8. Testing

Function-style scripts under `dev/tests/`, run by direct invocation (the repo idiom;
`runtests` does not accept these), via the MATLAB MCP.

1. **Integration test** (`test_tess_tangents.m`) on a real FreeSurfer-registered cortex
   (available from the BIDS/OMEGA import work; pick a subject whose cortex has
   `Reg.Sphere`). Assertions:
   - Per hemisphere, the nxr solve reports `gaussBonnetSatisfied`.
   - `U`, `V` are `[nF×3]`; every face assigned.
   - Orthonormality: `|U|≈1`, `|V|≈1`, `|U·V|≈0` per face; `U×V` aligned to the face
     normal (positive dot product).
   - Singularities are the four poles (two per hemisphere), indices all `+1`.
   - Determinism: two runs give identical `U`, `V`.
   - Storage round-trip: after a default run, reloading the surface yields
     `TangentFrame` with matching `U`, `V`.
2. **Guard test:** on a surface without `Reg.Sphere` (e.g. a bare `tess_sphere` saved as
   a surface), `tess_tangents` errors with `tess_tangents:noRegSphere`.

> Test-fixture note (to resolve in the plan): confirm a registered cortex is available
> in the local DB; if not, construct a synthetic two-hemisphere fixture (two offset
> spheres + a `Reg.Sphere` + a 'Structures' atlas or y-split) so the test is
> self-contained.

## 9. Out of scope / future

- **Per-vertex frame** via Hodge star / DEC transfer of the per-face field (the reserved
  `Domain` slot). This is the natural next step for per-vertex source-map analysis.
- **Source-map analysis** that expresses vector fields in this frame.
- GUI/process integration.

## 10. Decisions log (from brainstorming)

- **Frame domain:** per-FACE now (nxr-native); per-vertex deferred to a Hodge/DEC
  transfer step.
- **Solve structure:** per-hemisphere always — hemispheres are disconnected genus-0
  spheres; each solved independently (χ = 2, two poles, indices `[+1, +1]`). General
  principle for all nxr solves in Brainstorm.
- **Singularities:** FreeSurfer registration-sphere poles (max/min z) per hemisphere.
- **Interface:** SurfaceFile-based; computes, stores `TessMat.TangentFrame`, returns
  `[U, V]`; `'NoSave'` option to skip storage.
- **Dependency:** nxr-only, SPM-style guarded install, no MATLAB fallback.
- **Docs/tests:** design+plan in `dev/`; tests in `dev/tests/`.
