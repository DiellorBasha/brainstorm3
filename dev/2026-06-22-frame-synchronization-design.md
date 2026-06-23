# Frame Synchronization: canonical tangent frame for operator eigenmodes

**Date:** 2026-06-22
**Status:** design approved (architecture decisions made); implementation pending review of this doc.

## Problem

Vector/complex operator eigenmodes are coefficients in a per-vertex tangent frame:
- **Connection Laplacian** (complex): `field(v) = real(U[v,k])·e1(v) + imag(U[v,k])·e2(v)`.
- Decoding a complex eigenmode to a 3-D tangent vector (for the wavelet atom display, 3-D
  seeding, or any geometric analysis) REQUIRES the exact frame the operator was built in.

**The risk (why this matters):** today no frame is recorded with the operator or its
eigenmodes, and the one frame that *is* stored — in the manifold file — can be a DIFFERENT
frame than the operator uses. Interpreting eigenmodes with the wrong frame is silently wrong.

## Findings (investigation 2026-06-22)

**nxr (`~/workspace/research/code/nxr-compute`, C++/CMake/MEX):**
- Canonical frame = geometry-central `vertexTangentBasis`: `e1 = basisX` (first-outgoing-
  halfedge axis), `e2 = n × e1`, `normal = vertexNormals`. DETERMINISTIC (no RNG, stable
  ordering) → two calls on the same mesh+embedding give bitwise-identical frames.
  (`src/vertex_frames.cpp:16`; nxr test `test_facets.cpp:50` asserts grid==vertexGrid 1e-12.)
- The manifold `Embedded.vertex.grid` PACKS the frame as a complex `[nV×3]`:
  `grid = e1 + i·e2` (so `e1=real(grid)`, `e2=imag(grid)`; verified orthonormal, |e1|=|e2|=1,
  e1⊥e2⊥n, right-handed). Row norm √2 is the packing, not a bug.
- The Connection-Laplacian operator decodes in the RAW (levi-civita) `vertexTangentBasis`
  and already computes `frameE1/frameE2` internally (`src/connection_laplacian.cpp:249`),
  but the `nxr_compute('operators','laplacian','connection')` path returns ONLY the matrix —
  the frame is not exposed.
- The **manifold defaults to the `'trivial'` gauge**, whose grid is COMBED
  (`rotation .* grid`, `src/facets.cpp` GaugeFacet::grid) — a DIFFERENT per-vertex frame than
  the operator's raw frame. ← the concrete desync.
- **Dirac uses only vertex normals; quaternions are in AMBIENT world coordinates**
  (`src/dirac_operator.cpp:33` uses `vf.normals` only). It is frame-FREE → no tangent-frame
  desync; the existing full-quaternion Dirac wavelet is already consistent.

**Brainstorm:**
- `tess_operators` stores `Operator`, `Mass`, `GlobalVertices/Faces`, `Provenance` — NO frame
  (`operatormat` schema has no frame field).
- `tess_manifold` stores `Embedded/Intrinsic/Extrinsic/Gauge` per hemi; default gauge trivial.
- `eigen_` nodes store complex `Phi` with no attached frame.

## Decisions

1. **Frame home: manifold + operator, validated.** Store the explicit raw frame
   `(e1,e2,normal)` in BOTH the manifold file and the operator file, each captured from the
   same deterministic nxr routine, with a test asserting `manifold.frame == operator.frame ==
   nxr.vertexFrames` bitwise. Consumers read the manifold (the single reference); the operator
   copy keeps the decoding basis next to the eigenmodes.
2. **Canonical convention: raw levi-civita.** The operator-aligned frame is the raw
   geometry-central `vertexTangentBasis`. The `'trivial'` gauge combing becomes a display-only,
   SEPARATELY-stored field (`Gauge.vertex.rotation`, already present). No operator rebuild.
3. **Dirac is ambient/frame-free** — documented; needs no frame. (Connection is the only
   frame-dependent vector operator for now.)

## Implementation plan (phased)

**Phase A — nxr exposure (C++).** Add a first-class `nxr_compute('vertexFrames', h)` (and
`'faceFrames'`) MEX command returning `e1,e2,normal` `[nV×3]` from `geometry::vertexFrames`
(routine already exists). Pure exposure, no math change. Rebuild the MEX (CMake) and deploy
into `~/.brainstorm/plugins/nxr-compute/...` (⚑ stale-binary trap — overwrite the loaded
binary). Add/extend an nxr C++ test asserting `vertexFrames == real/imag(vertexGrid)`.

**Phase B — schema.** Add a `Frame` field to `operatormat` (and confirm the manifold carries
explicit `e1/e2/normal`, not only the packed grid). `Frame` = 1×2 cell of
`struct('e1',[nV×3],'e2',[nV×3],'normal',[nV×3])` per hemi (empty for ambient/scalar variants).
db migration if needed (⚑ struct-array widening: rebuild template, don't in-place patch).

**Phase C — capture at build time.** `tess_operators`: for frame-dependent variants
(Connection Laplacian), fetch `nxr_compute('vertexFrames', h)` in the same hemi loop that
builds the operator and store it in `OperatorMat.Frame{hh}`. `tess_manifold`: store the explicit
raw frame in `Embedded.vertex` (e1/e2/normal) alongside the gauge grid.

**Phase D — validation (the "trusted" guarantee).** A test that, for a real surface, builds the
manifold + operator and asserts the stored frames are bitwise-equal to each other and to a fresh
`nxr_compute('vertexFrames')` call. Document determinism as the guarantee.

**Phase E — consumers.** A frame accessor (`tess_frame`-style getter) returning the canonical
`(e1,e2,normal)` for a surface/variant. Connection wavelet/display: decode complex eigenmode
`z → real(z)·e1 + imag(z)·e2`, add U(1) phase orientation control + glyph in
`panel_eigenwavelet_options`. Document Dirac as ambient (no frame).

## Non-goals / notes
- The neighbor-angle "ambient smoothness" probe is NOT a valid frame check (connection modes are
  COVARIANTLY smooth, ambiently rotating). Validation is structural (determinism + bitwise eq).
- No change to operator math or the trivial-gauge connection; raw frame is canonical.
- Dirac wavelet (full quaternion, merged) is unaffected.
