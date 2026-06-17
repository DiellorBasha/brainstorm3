# Intrinsic face Dirac (nxr-compute) + intrinsic face-domain Helmholtz (Brainstorm)

**Date:** 2026-06-17
**Status:** approved (design)
**Spans two repos:** Phase A in `nxr-compute` (C++); Phase B in `brainstorm3` (MATLAB).
**Builds on / motivated by:** the face-domain Helmholtz prototype, which proved the only
face Dirac nxr exposes (`diracFaceD`) is EXTRINSIC (Gauss-map) and therefore inconsistent
with the intrinsic cotan Poisson the vertex pipeline uses — a planted pure-skew-gradient
field recovered at corr 0.36 / HarmFrac 4790. An apples-to-apples intrinsic face Helmholtz
needs an INTRINSIC face Dirac, which does not exist yet.

## Background (the operators, established from nxr source/tests)

- Vertex extrinsic `matrix`/`diracD` [4F×4V]: block from **face-normal** differences (Gauss map).
- Vertex intrinsic `matrixIntrinsic`/`diracIntrinsicD` [4F×4V]: same block structure but from
  **vertex positions** `(P_r − P_q)/(2A_f)` (the immersion edge vector). Its `W_F`-Galerkin
  square's scalar block == the cotan Laplacian (Crane spin-connection property). The vertex
  Helmholtz pipeline uses THIS + the cotan LBO → self-consistent.
- Face extrinsic `matrixFace`/`diracFaceD` [4V×4F]: per vertex `v`, over cyclically ordered
  incident faces, block from **face-normal** differences `−leftMulImag(N_{f_{k+1}} − N_{f_{k-1}})/(2 Ã_v)`,
  `Ã_v` = barycentric vertex dual area. `D̃ᵀ W_V D̃ == diracFace(1)` (extrinsic, curvature-coupled).

## Decisions (settled with user)

- **Dual immersion = barycentric face centroids** (robust; consistent with the barycentric
  `W_V`; circumcenters are fragile on cortical triangles).
- **Correctness bar = functional + self-consistent square** (NOT forced equality to the
  circumcentric DEC 2-form Laplacian): (i) first-order (kills constant face fields), (ii) its
  `W_V`-Galerkin square scalar block `K̃ᵢₙₜ` is a valid PSD Laplacian used as-is in the dual
  Poisson, (iii) a planted dual-skew-gradient field recovers with HarmFrac→0, recovered
  potential corr ≈ 1. Report (informational) the relative distance `‖K̃ᵢₙₜ − DEC-2form-Lap‖`.
- **One spec, two phases, executed A → build/release → B.**

## Phase A — `diracFaceIntrinsicD` in nxr-compute (C++)

The intrinsic face Dirac `D̃ᵢₙₜ : ℍ^F → ℍ^V` [4V×4F] mirrors `matrixFace` but replaces the
incident face NORMALS with the incident face CENTROIDS (the dual immersion):

    for each interior vertex v, Ã_v = vertexDualArea(v); incident faces f_0..f_{d-1} (cyclic),
        centroids C_k = (1/3) Σ positions of f_k's vertices;
        s = -1/(2 Ã_v);
        block(v, f_k) += s · leftMulImag( C_{(k+1) mod d} − C_{(k-1+d) mod d} );

Closed-mesh v1 (same boundary guard as `matrixFace`). Geometry-only (vertex positions + dual
areas), cached.

**Files (nxr-compute):**
- `src/dirac_operator.cpp`: add `matrixFaceIntrinsic(Manifold&)`; declare in
  `include/nxr/compute.h` (`ops::dirac::matrixFaceIntrinsic`).
- `include/nxr/facets.h` + `src/facets.cpp`: `OperatorsFacet::diracFaceIntrinsicD()`
  (const-ref, cached via a new `OperatorId::DiracFaceIntrinsicD`).
- `bindings/mex/src/nxr_compute_mex.cpp`: add `'diracFaceIntrinsicD'` to the `operators`
  dispatch + the usage doc block + the error string list.
- `bindings/wasm/src/nxr_compute_wasm.cpp` + `bindings/wasm/js/index.d.ts`/`index.mjs`:
  expose `diracFaceIntrinsicD` (parity with `diracFaceD`).

**Tests (nxr-compute):**
- C++ `test/test_dirac_operator.cpp`: `D̃ᵢₙₜ` size `[4V×4F]`; `D̃ᵢₙₜ·Uf ≈ 0` (kills 4 constant
  face quaternions); `K̃ᵢₙₜ = (D̃ᵢₙₜᵀ W_V D̃ᵢₙₜ)` scalar block is symmetric PSD and kills the
  constant; print `‖K̃ᵢₙₜ − Δ₂(DEC)‖/‖Δ₂‖` (informational).
- MATLAB `bindings/mex/test/test_dirac_first_order.m`: size + kills-constants + cache stability,
  mirroring the existing `diracFaceD` checks.

**Build/validate:** local `bash scripts/build.sh Release && ./build/Release/test_dirac_operator`;
then run the MATLAB test against the fresh mex. Release via the (now hang-tolerant) `publish-mex`
workflow on a `v*` tag; install the new `nxr_compute.mexmaca64` into the plugin folder (back up
the current one per the `.bak` convention).

## Phase B — intrinsic face Helmholtz in Brainstorm (MATLAB)

1. **`tess_operators.m`**: in the per-hemisphere Dirac build, also pull
   `nxr_compute('operators', h, 'diracFaceIntrinsicD')` and store it (`FirstOrderFaceInt{hh}`,
   [4V×4F]) alongside the existing `FirstOrderInt`/`FirstOrderExt`; bump the operator schema/
   provenance. (Requires the Phase-A mex installed.)
2. **`bst_dirac_helmholtz_face.m`** (rewrite the prototype): use `D̃ᵢₙₜ`; form the intrinsic
   face Laplacian `K̃ᵢₙₜ` = scalar (w-w) block of `D̃ᵢₙₜᵀ W_V D̃ᵢₙₜ`; the dual Poisson
   `K̃ᵢₙₜ ψ̃ = (mass) · ω̃` solves the stream/potential on FACES; reconstruct the solenoidal /
   irrotational component fields via the `D̃ᵢₙₜ`-consistent dual skew-gradient (derived so that
   the planted-field test passes — the residual derivation piece).
3. **`bench_dirac_face_helmholtz.m`**: re-run the vertex-vs-face comparison; the face decomposition
   is now intrinsic and apples-to-apples.

## Testing (Brainstorm, Phase B)

- `dev/tests/test_dirac_helmholtz_face.m` (rewrite): the planted pure-skew-gradient field now
  recovers with **HarmFrac ≈ 0** and recovered potential **corr ≈ 1** with the planted one (the
  gate that the intrinsic operator + dual Poisson + dual reconstruction are consistent);
  reconstruction exact `Virr+Vsol+Vharm == Jf`; shapes; HarmFrac ∈ [0,1].
- Regression: the existing vertex `bst_dirac_helmholtz` suite + all other suites unaffected.

## Risks / open derivation

- **Dual reconstruction operator** (face-potential → component field) is the one piece still to
  be derived; the planted-field HarmFrac→0 test is the gate. If the barycentric construction
  does not give clean recovery, fall back to the circumcentric immersion (Q2-B) for `K̃ᵢₙₜ`.
- **Build coupling:** Phase B is blocked until the Phase-A mex is built + installed.
- **C++ validation only via build** (no MATLAB-MCP compile): local CMake build + tests.

## Non-goals

The face-based LEADFIELD and the face-Dirac eigenbasis INVERSE (the phase after this); open-mesh
boundary handling; circumcentric construction unless barycentric fails the planted-field gate.
