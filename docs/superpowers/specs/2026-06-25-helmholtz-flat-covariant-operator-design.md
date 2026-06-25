# Helmholtz decomposition on the flat covariant operator — design

**Goal:** Rebuild `bst_helmholtz`'s **decomposition** so it uses the correct, self-consistent
differential operator — the flat covariant ("Ambient") operator from nxr-compute — instead of the
inconsistent intrinsic-Dirac-plus-mismatched-FEM-gradient mix it uses today. The decomposition becomes a
machine-exact Helmholtz–Hodge of the full 3-component source current, with the normal degree of freedom
represented honestly.

**Status:** design, ready to plan.

**Author:** Diellor Basha, 2026-06-25

---

## 1. Motivation

`bst_helmholtz`'s vertex decomposition stitches three *mismatched* discretizations: the quaternionic
**intrinsic Dirac** (`FirstOrder.Intrinsic`, a repurposed spin-transformation operator) for divergence/curl,
a **separate FEM gradient** for reconstruction, and the **cotan Laplacian** for the Poisson solve. Because
the operators are not a consistent family, the decomposition is not a true Hodge decomposition. Verified on a
real alpha frame:

- **Round-trip residual ≈ 0.68** (`J ≠ grad φ + n×grad ψ`) — two-thirds of the field is lost.
- **Sign inverted**: the operator outputs `−divergence` and `−vorticity` (each correlates −1.0 with the FEM
  reference), and that sign rides through the Poisson solve, so a source reads as a sink.
- **`HarmFrac` is meaningless**: it reports ~4e-8 while the actual residual by its own formula is ~0.68; on a
  genus-0 cortex the true cohomological harmonic part is identically zero, so the leftover it bins as
  "harmonic" is really the normal component plus reconstruction error.

The fix is to compute the decomposition from a **single consistent operator** so divergence, the Laplacian,
and reconstruction cannot drift. nxr-compute already provides exactly this — the flat covariant ("Ambient")
operator — and the same file's **face-domain** path already proves the consistent approach (machine-exact,
`HarmFrac` 5e-24). We bring the vertex decomposition to that standard.

## 2. The operator (verified)

The flat covariant operator is nxr's **Ambient** covariant operator:
`gauge('levi-civita', coupling='ambient')` → `covariantLaplacian`, with the first-order covariant gradient
`covariantGradient` (`G : 3N→3E`, flat full-frame transport `P_ij = Fⱼᵀ Fᵢ`) and edge weights `⋆`. It is the
operator we identified and tested as correct for divergence/curl analysis. Empirically established during
brainstorming:

- **No curvature artifact:** a constant ambient field gives `div = curl =` machine zero, *even on folds*
  (the curvature-energy `D_N` operator fails this — it manufactures spurious sources on every fold).
- **Carries the shape operator but is flat:** its per-edge transport includes the normal↔tangent tilt, yet
  holonomy is ~1e-15 — so it transports the 3-D vector correctly across folds without injecting a curvature
  bias. This is exactly what "first-order differential analysis" wants (flatness is the virtue); the
  curvature-energy flavor belongs to the spectral/eigenbasis role, not here.
- **div/curl match the validated reference** (`divFEM = Gx·Jx+Gy·Jy+Gz·Jz`, `omegaFEM`) at corr ±1.0.
- **Consistent Hodge round-trips to machine zero** (gate 1: 2.5e-14 on a field built in the
  `{grad, n×grad}` span) — versus the current path's 0.68.
- **Normal component participates:** `‖div_full − div_tangential-only‖/‖div_full‖ = 0.81`.

Why **not** the curvature-energy `Δ₃ + D_N`: it injects div/curl proportional to (curvature × field), i.e.
a folding artifact, and would smear the geometry into the field analysis — the opposite of what
differentiation needs. `D_N`'s home is a curvature-aware *eigenbasis*, which already exists separately
(`tess_eigen('Dirac')`) and is out of scope here.

Why **not** a tangential-only operator (connection Laplacian / the current intrinsic path): it discards the
3-D rotation/scaling of the current; the commitment is the full 3-component ambient field.

## 3. Architecture — two layers

Mirrors the existing operator pattern (operators are built in `tess_operators` from nxr-compute and fetched
by consumers via the surface file).

### 3.1 Layer 1 — `tess_operators`: the `'Covariant'` operator node

A new operator variant **`'Covariant'`** alongside `'Laplace-Beltrami'`, `'Connection Laplacian'`, `'Dirac'`,
`'Dirac-Face'`, `'Hodge-Face'`. Built per hemisphere from nxr:

- the covariant gradient `G` (`covariantGradient`),
- the covariant Laplacian (`covariantLaplacian`, with its cached factorization for the Poisson solve),
- the edge weights `⋆`,
- the per-vertex frames/normals needed to (a) contract the *scalar* surface divergence and vorticity out of
  `G·J` and (b) split off the normal amplitude.

Stored as a derived-anatomy node on the surface file with provenance + schema version, behind the existing
recompute-on-stale gate, and fetched via `bst_get_operator_node(SurfaceFile, 'Covariant')`.

The precise contraction that yields the scalar `div`/`curl` from nxr's vector `G` is pinned during this
layer's verification (gates 1–2 are the arbiter); if `G` alone does not cleanly yield the scalar pair, the
node also pulls nxr's scalar gradient (`Gradient3D`) + cotan Laplacian as the matched scalar family. Either
way, divergence, the Laplacian, and reconstruction all derive from one gradient and cannot drift.

### 3.2 Layer 2 — `bst_helmholtz` consumes the node

`bst_helmholtz('Prepare', …)` fetches the `'Covariant'` node via the surface file (instead of
`Dirac.FirstOrder.Intrinsic`). `'Frame'` computes, all from the single `G`:

- **`Div`, `Curl`** — flat-covariant scalar divergence / vorticity (full-3-D; the feature fields for
  sources and vortices).
- **`Phi`, `Psi`** — scalar potentials from the consistent Poisson (`K = gᵀ⋆g`, cached factorization).
- **`Virr = grad φ`, `Vsol = n×grad ψ`** — reconstruction; round-trips the tangential field exactly.
- **`Jnormal = J·n`** — the normal/radial amplitude, the third scalar DOF.
- **sign** applied per the gate-4 calibration (so a source reads as a source).

The complete 3-DOF decomposition is `J = grad φ + n×grad ψ + (J·n)·n` — two tangential potentials plus the
normal scalar, the full 3-component current with nothing discarded. (Rationale: `φ`/`ψ` carry exactly the 2
tangential DOF; the normal DOF cannot be represented by a scalar potential's gradient, so `Jnormal` is
required for completeness — and is where a radial source on a flat patch lives, since such a source has zero
surface divergence.)

## 4. Output contract

The `Ht` struct keeps the field names the consumers read — `Div`, `Curl`, `Phi`, `Psi`, `Fmag`, `Vtot`,
`Virr`, `Vsol`, `Cores`, `Sources` — now with **correct, sign-calibrated** values. Changes:

- **Add `Jnormal`** (the normal-amplitude scalar field).
- **Replace the dishonest fields**: `Vharm`/`Hmag` (the mislabeled "harmonic" bin = normal component +
  reconstruction error) are removed in favor of `Jnormal` + a true residual that is now ~0. `HarmFrac`
  becomes a genuine ~0 number, or is dropped (on genus-0 it is meaningless).

Field names that consumers depend on are preserved; only `Jnormal` is added and the two dishonest fields
retired.

## 5. Scope boundaries

**Untouched:**
- The **`Dirac` eigenbasis** (`tess_eigen('Dirac')` → the curvature-aware τ-mixed Dirac) and the
  **scale-axis smoothing** (`bst_eigenfilter`). These are correct and curvature-aware by design; the `Dirac`
  node keeps serving the eigenbasis — it simply stops being the decomposition operator.
- The **face-domain** coupled-variational Hodge (already machine-exact).

**Out of scope (follow-on spec):** the **source-feature detector** — sources (`Div`/`Jnormal` extrema),
vortices (`Curl` extrema), saddles (`φ`/`ψ` separatrices). It is a separate effort that begins only after
this operator foundation is merged and validated, and sits on the now-correct decomposition.

## 6. Phasing and validation

**Phase A — Layer 1 (`tess_operators` `'Covariant'` node).** Ships only when all five gates pass; gates 1–3
are already empirically green:
1. **Exact round-trip Hodge** — machine zero (shown: 2.5e-14, vs the current 0.68).
2. **Matches validated div/curl** — corr ±1.0 with `divFEM`/`omegaFEM` (shown).
3. **No curvature artifact** — constant ambient field → div/curl machine zero, even on folds (shown).
4. **Correct sign** — a synthetic point source (`∇·J > 0`) yields `φ` with the agreed-sign extremum there;
   the sign is fixed and documented once, here.
5. **Complete 3-DOF** — `‖J − (grad φ + n×grad ψ) − (J·n)·n‖` is machine zero.

**Phase B — Layer 2 (`bst_helmholtz` consumption).** Emits the corrected fields + `Jnormal` from the node.
**Consumer regression gate:** `bst_operators`, `bst_divergence`, `bst_curl`, `view_helmholtz`,
`process_helmholtz_events`, `process_vortex_track` all run unchanged (same API) and produce sensible,
sign-calibrated output, verified on the alpha-vortex frame against the prior qualitative results (vortex
cores still found, sources sensible, sign now physical).

## 7. Risks / notes

- **Scalar-from-vector-gradient contraction**: the one piece pinned during Phase-A verification (gates 1–2),
  not assumed. nxr's scalar gradient + cotan Laplacian is the fallback matched family if needed.
- **nxr handle / gauge call**: Layer 1 must obtain the nxr handle for the surface and call the Ambient gauge;
  follow `tess_operators`' existing nxr-handle pattern (and the canonical-mesh / `nxr_safe_create` guard).
- **Consumer field dependence**: confirm which `Ht` fields each of the 6 consumers reads before retiring
  `Vharm`/`HarmFrac`, so nothing breaks silently.
- **Sign convention**: must be calibrated on a synthetic source and documented; the current code's silent
  inversion is the cautionary precedent.
- **Boundary**: the cortex hemisphere has a medial-wall boundary; the round-trip gate uses in-span synthetic
  fields (boundary-robust) and the consumer gate uses real data — both already exercised in brainstorming.
