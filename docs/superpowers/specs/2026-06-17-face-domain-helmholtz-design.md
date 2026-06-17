# Face-domain Helmholtz via the dual Dirac (D̃) — prototype

**Date:** 2026-06-17
**Status:** approved (design)
**Builds on:** bst_dirac_helmholtz (vertex pipeline), tess_operators (nxr handle build), bst_vortex_persistence.

## Motivation

The current Dirac/Helmholtz pipeline carries the source field on **vertices** (3D per
vertex) and uses the vertex Dirac `D = diracIntrinsicD` ([4F×4V]). nxr also exposes the
**dual face Dirac** `D̃ = diracFaceD` ([4V×4F]), which acts on a field living on **faces**.
A face is a flat triangle with one canonical tangent plane/normal, so a per-face vector has
an unambiguous tangential/normal split — more principled than the averaged vertex normal.
This prototype runs a **face-domain Helmholtz decomposition** on the existing solution
(interpolated to faces) and compares it to the vertex pipeline, establishing the operator
path a future **face-based leadfield + face-Dirac eigenbasis inverse** will reuse.

The inverse stays **full unconstrained 3D** throughout — the face change moves the DOF
domain (vertices → face centroids), it does NOT project to the tangent plane.

## Decisions (settled with user)

- Extract curl/div with the **dual face-Dirac `D̃`** (apples-to-apples with the current
  Dirac pipeline; exercises the exact operator the face inverse will use).
- Build a **reusable `bst_dirac_helmholtz_face`** (the face inverse will reuse it), not
  throwaway bench code.
- Divergence in this dual path uses **vertex normals** (a known asymmetry); the prototype
  quantifies its effect rather than hiding it.

## A. Components

1. **`toolbox/math/bst_dirac_helmholtz_face.m`** (new) — face-domain analog of
   `bst_dirac_helmholtz`, `Prepare`/`Frame` dispatch:
   - `Prepare(DiracFaceOp, LBO, Surf)`: per hemisphere, obtain `D̃` [4Vh×4Fh] (see B on
     sourcing it), and cache: vertex cotan stiffness `K` + mass `M` + Cholesky of the pinned
     `K` (reuse LBO — no new Laplacian), `Surf.VertNormals(vH,:)`, local faces `Floc`, the
     per-face FEM gradient operators `Gx/Gy/Gz` (same construction as `bst_dirac_helmholtz`
     `Prepare`), and `GlobalVertices vH`.
   - `Frame(Op, Jf)`: `Jf` = per-face ambient field as `[3*nF x 1]` (Jx,Jy,Jz stacked) or
     `[nF x 3]`. Per hemisphere: embed pure-imaginary quaternion `qF` [4Fh]
     (`qF(2:4:end)=Jx; 3:4:=Jy; 4:4:=Jz`); `qV = D̃ * qF` [4Vh]; **vorticity** `omV = qV(1:4:end)`
     (w-part, normal-free); **divergence** `dvV = imag(qV)·N_vertex` per vertex; Poisson
     `psi = K\(M·omV)`, `phi = K\(M·dvV)` (mean-zero, pinned — reuse `i_poisson` /
     `bst_dirac_helmholtz('PoissonSolve',...)`); component fields per face
     `gphi=[Gx Gy Gz]·phi`, `gpsi=...`, `Vsol = Nf × gpsi`, `Virr = gphi`,
     `Vharm = Jf - Virr - Vsol`; HarmFrac from the face-area-weighted residual energy.
     Cores/Sources via `bst_vortex_persistence` on psi/phi (per hemisphere, with the existing
     `.hemi` tagging convention). Returns `Ht` with the same field names as the vertex
     `Frame`: `Curl/Div/Psi/Phi` (per vertex), `Vtot/Virr/Vsol/Vharm` (per face),
     `Fmag/Hmag`, `HarmFrac`, `Cores/Sources`.

2. **`dev/benchmarks/bench_dirac_face_helmholtz.m`** (new) — comparison driver:
   load the S01 alpha dSPM frame (kernel via `bst_inverse_dirac`, like the existing dirac
   benchmarks); compute the **vertex** `Ht` with `bst_dirac_helmholtz('Frame', Op, Jt)`;
   barycentric-interpolate `Jt` (3D/vertex) → face centroids `Jf` (mean of the 3 vertex
   vectors per face); compute the **face** `Ht` with `bst_dirac_helmholtz_face`; compare and
   visualize; save a PNG.

## B. Sourcing D̃

`tess_operators` builds per-hemisphere nxr handles and already pulls `diracIntrinsicD` /
`diracD`; it does NOT currently store `diracFaceD`. For the prototype, `Prepare` builds its
own per-hemisphere handle (`nxr_safe_create(Vloc, Floc)`), calls
`nxr_compute('operators', h, 'diracFaceD')`, then `nxr_compute('destroy', h)` — mirroring the
`tess_operators` loop. (Promoting `diracFaceD` into the stored Dirac operator is deferred to
the face-inverse phase.)

## C. Convention validation (critical — not assumed)

`D̃` is the dual operator; its quaternion convention (which part is vorticity vs divergence,
and sign) is NOT assumed to match the vertex `D`. Validate on the real cortex with synthetic
per-face fields of known content:
- a **pure-rotational** face field (e.g. `Nf × (c - centroid)` about a core) → expect nonzero
  vorticity dominating divergence;
- a **pure-gradient/divergent** face field (e.g. radial `(centroid - c)`) → expect divergence
  dominating vorticity.
Confirm the w-part is vorticity and `imag·n` is divergence (and the sign) before trusting the
decomposition; a Poisson round-trip (`K psi = M omega` recovers `psi`) is also checked.

## D. Comparison metrics + visualization (vertex vs face, same frame)

- harmonic-energy fraction (face vs vertex);
- #vortices/#antivortices and #sources/#sinks **per hemisphere** (persistence-gated);
- top-core vertex/location agreement;
- curl/div magnitude distributions (min/median/max);
- side-by-side cortex figure: ψ colormap + gated cores, vertex pipeline vs face pipeline.
Saved to `dev/benchmarks/bench_dirac_face_helmholtz.png`.

## E. Testing

- **`dev/tests/test_dirac_helmholtz_face.m`**: the C convention validation (rotational/
  gradient face fields), Poisson round-trip, finite + shape checks
  (`Curl/Div` [nV×1], component fields [nF×3]), `HarmFrac ∈ [0,1]`, and `Cores/Sources` are
  struct arrays carrying `.persistence`/`.hemi`.
- **Integration**: the benchmark runs end-to-end on the real frame and emits the figure;
  sanity asserts (finite, harmonic fraction in range, cores nonempty).

## Non-goals (next phase)

The face leadfield (gain for 3D dipoles at face centroids) and the face-Dirac eigenbasis
inverse; promoting `diracFaceD` into the stored Dirac operator; the face-native tangential
Hodge variant (option B). This prototype is analysis-only on the existing vertex solution.
