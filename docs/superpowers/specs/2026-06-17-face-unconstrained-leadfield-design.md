# Full-unconstrained face leadfield + head model — Design

**Date:** 2026-06-17
**Status:** design (awaiting user review)
**Repo:** brainstorm3 (`development`)

---

## 1. Motivation

The face-based pipeline (Dirac/Helmholtz, and the upcoming face-Dirac eigenbasis inverse) needs a forward model whose sources are FULL-3D UNCONSTRAINED dipoles at face centroids — exact parity with the vertex leadfield (raw Cartesian point dipoles), just relocated to faces. Per the standing directive: *all inverse methods operate on the full 3D vector; never a tangential/constrained inverse.*

**Prior work is misaligned.** `bst_face_leadfield` / `bst_face_headmodel` (committed `84710660`, `3fdeb3b6`) produce only **constrained** (1 col/face along n̂) or **"loose"** (3 cols in a mixed `{n̂·A_f, U_f, V_f}` frame) leadfields, and derive the frame from **`tess_tangents`** (deprecated; memory `dirac-quaternion-frame-rotation`). Neither is clean full-unconstrained.

**The rigorous model is simpler:** the overlapping-spheres Sarvas solver (`bst_meg_sph`) is position-agnostic, so the raw gain at a face centroid is already `[nCh × 3]` Cartesian (x,y,z) — three independent orientation columns. Stacking over faces gives `[nCh × 3F]` directly, with no frame, no projection, no `tess_tangents`.

## 2. Design

Add an `'unconstrained'` Mode to both functions (and make it the default for the rigorous path). Constrained/loose modes are retained untouched for back-compat.

### `bst_face_leadfield(SurfaceFile, Channel, Param, 'Mode','unconstrained')`
- **Geometry from `tess_manifold` (NOT `tess_tangents`, which is deprecated).** `tess_manifold(SurfaceFile)` returns the canonical nxr geometry backbone; assemble the global per-face arrays by scattering each hemisphere's data through its `GlobalFaces` map:
  - centroids `x_f = M.Embedded(hh).face.centroid` `[F×3]` (verified == barycentric mean, exact)
  - face normals `n̂_f = M.Embedded(hh).face.normal` `[F×3]` (gauge-consistent, **unit**, orientation verified vs geometry) — used ONLY for `GridOrient` display, never for projection
  - areas `A_f = ½|e₁×e₂|` `[F×1]` (pure geometry; `tess_manifold` carries no face-area field) — metadata only (`GridAreas`), not applied to the gain
  - This also makes the leadfield's face indexing/orientation match the manifold the downstream face-Dirac inverse will use.
- **Leadfield:** block over faces, `L_face(:, 3f-2:3f) = bst_meg_sph(x_f', Channel, Param)` — i.e. the raw `[nCh × 3F]` Sarvas gain. **No area weighting** (parity with the vertex model; area/depth weighting belongs in the inverse's source prior).
- Returns `FaceGeom` with `.Centroids .Normals .Areas` (no `.U/.V` in this mode).

### `bst_face_headmodel(BaseHeadModelFile, 'Mode','unconstrained')`
- Calls `bst_face_leadfield(...,'Mode','unconstrained')`.
- Saves a standard surface head model: `Gain [nCh×3F]`, `GridLoc = centroids`, `GridOrient = n̂_f` (geometry normals), `HeadModelType='surface'`, `MEGMethod = BaseHM.MEGMethod`, `nComponents = 3`, `isFaceBased = 1`, `GridAreas = A_f`. `GridU/GridV = []` (no frame in unconstrained mode). Comment tag `Face-unconstrained`.
- Reuses the base os_meg `Param` + MEG channel mapping (no sphere re-fit), exactly as the existing function does.

## 3. Correctness contract (tests)

`dev/tests/test_face_leadfield_unconstrained.m`:

| Check | Assertion |
|---|---|
| Shape | `L_face` is `[nMEG × 3F]`, all finite |
| Pure Sarvas | `L_face(:,3f-2:3f)` equals `bst_meg_sph(x_f(f,:)', …)` for a sampled face (it IS that — exactness, no frame) |
| **Constrained consistency** | the existing `'constrained'` column `f` equals `L_unconstr(:,3f-2:3f) · n̂_f · A_f` (ties the new mode to the validated constrained one) |
| Geometry from `tess_manifold` | unconstrained path obtains centroids/normals from `tess_manifold` (deprecated `tess_tangents` is never called); centroids == barycentric mean (exact), normals unit |
| Observability ceiling | the whitened face leadfield's effective rank ≈ the vertex leadfield's (~tens of DOF — the MEG physics ceiling, memory `bst_inverse_dirac`); a sanity bound, not bit-equality |
| Head model struct | saved HM has `nComponents==3`, `isFaceBased==1`, `size(Gain,2)==3*nF`, `GridOrient` unit normals |

## 4. Validation benchmark

`dev/benchmarks/bench_face_leadfield.m`: on the S01 head model,
1. Build the unconstrained face head model.
2. **End-to-end parity:** run `bst_inverse_dirac` (MEG-only, dspm2018) with the FACE head model on the real alpha/M100 frame, and compare the resulting unconstrained source map to the VERTEX `bst_inverse_dirac` map (interpolated for display). Expect the same large-scale localization (the observable subspace is basis-invariant ~tens of DOF, memory `bst_inverse_dirac`) — the face map is finer-resolution but co-located.
3. Report: gain shape, whitened-rank (face vs vertex), peak-localization agreement; save a side-by-side PNG.

This is the empirical input to the face-Dirac eigenbasis INVERSE (next phase).

## 5. Out of scope (next phases)

- Face-Dirac eigenbasis inverse (reuses this leadfield + the face Dirac/`gradFace`/`lapFace` eigenmodes).
- EEG/OpenMEEG face forward (OS-MEG only here; OpenMEEG is also position-agnostic but deferred).
- Retiring/refactoring the constrained/loose `tess_tangents` modes (left intact; only a new mode is added).

## 6. Self-review

- **Placeholders:** none — the construction is a direct call to the existing position-agnostic `bst_meg_sph`; geometry is the standard centroid/normal/area already used in `bst_dirac_helmholtz_face`.
- **Consistency:** the `constrained-consistency` test pins the new mode to the already-validated constrained leadfield, so we're not asserting correctness in a vacuum.
- **Scope:** one new mode in two existing functions + test + benchmark. No change to the forward solver, no new nxr work.
- **Directive parity:** full-3D unconstrained (3 Cartesian cols/face), no tangential projection, no `tess_tangents` — matches the vertex model exactly.
