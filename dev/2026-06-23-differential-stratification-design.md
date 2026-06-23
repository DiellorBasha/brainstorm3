# Differential module: field-type stratification + shared Poisson/Cholesky infrastructure

**Date:** 2026-06-23
**Author:** Diellor Basha (design assisted by Claude)
**Status:** Design — pending implementation plan

## 1. Motivation

`toolbox/differential/` currently holds `bst_gradient`, `bst_divergence`, `bst_curl`
(2-D tangent operators via the manifold DEC primitives), `bst_helmholtz` (3-D ambient
Hodge decomposition via the quaternion Dirac root), and the `bst_operators` orchestrator.
Two gaps:

1. The **differential side is not stratified** the way the eigen/operator side already is.
   `tess_operators` builds **Laplace–Beltrami** (ℝ, scalar), **Connection Laplacian**
   (ℂ Hermitian, tangent, with a per-vertex frame `e1/e2/n`), and **Dirac / Dirac-Face /
   Hodge-Face** (quaternion, 3-D), and `bst_eigen` dispatches on them. The differential
   verbs have no equivalent `switch` over field type, and no guard against ill-posed
   operations (e.g. curl of a scalar).

2. **Poisson/Cholesky is duplicated and inconsistent.** `bst_helmholtz` factorizes the
   pinned cotan stiffness once in `Prepare` and reuses it per frame (good). `bst_operators`'
   `i_poisson` **re-factorizes per time column** (a real performance bug). There is no
   shared factor cache and no `bst_poisson`.

This design introduces a field-type stratification across the differential verbs, a lazy
factor-caching helper (`tess_cholesky`), and a stratified Poisson solver (`bst_poisson`),
then refactors the existing engines onto that shared infrastructure.

## 2. Theory anchor (the stratification)

Under grad/curl/div sits the exterior derivative `d` on `k`-forms with `d² = 0`. On a 2-D
surface the de Rham complex is `Λ⁰ →d→ Λ¹ →d→ Λ²` (scalar → tangent vector → scalar):
there is **no vector-valued curl** — the curl of a tangent field is the scalar perpendicular
vorticity, and divergence is also a scalar. The surface Hodge star `⋆` on 1-forms is a 90°
rotation (multiplication by `i`), which is why tangent vectors are naturally complex and the
connection Laplacian is Hermitian.

`d² = 0` is metric-free, so `curl grad = 0` and `div curl = 0` hold at any curvature — the
Helmholtz/Hodge decomposition is always safe. But differentiating a *vector* requires
parallel transport, and curvature is the failure of second derivatives to commute. On a
surface the entire Riemann tensor is the Gauss curvature `K`, and the **Weitzenböck identity**
`Δ_Hodge = Δ_∇ + Ric`, with `Ric = K·g`, means the Hodge Laplacian (`δd + dδ`) and the
connection/Bochner Laplacian (`∇*∇`) **differ by `K` outright**. This is the
scalar-vs-tangent-vector face of the LBO-vs-connection-Laplacian split.

A current `J` is **ℝ³-valued**, not a tangent field: it is a section of the trivial ambient
bundle restricted to the surface, with a tangential and a normal part. Differentiating it
couples the two through the second fundamental form (Gauss/Weingarten):
`∇^amb_X Y = ∇^surf_X Y + II(X,Y)N` and `∇^amb_X N = −S(X)` (shape operator). This extrinsic
coupling is invisible to the LBO (metric only) and to the connection Laplacian (tangent only);
it is the `D_N` term in the extrinsic Dirac, `D² = Δ₃ + D_N`. A concrete instance: the surface
divergence of `J` picks up a mean-curvature term coupling its normal component to `H`
(see §6).

### The three strata

| Stratum | Field representation | Row layout | Operator family | Laplacian |
|---|---|---|---|---|
| **scalar** Λ⁰ | per-vertex scalar | `[nV × nT]` real | Laplace–Beltrami (cotan `K`, mass `M`) | `Δ = M⁻¹K` |
| **tangent** Λ¹ | tangent vector | DEC: per-face ambient `[3nF]`; connection: complex `[nV]` | de Rham `d/δ` (DEC) | Hodge `δd+dδ` **and** connection `∇*∇` (differ by `K`, Weitzenböck) |
| **ambient** | ℝ³ vector (tangential + normal) | per-vertex `[3nV]` | quaternion Dirac `D` (`D²=Δ₃+D_N`) | `D²` |

## 3. Operation matrix (router + guards)

`bst_operators` is the router. It infers (or is told) the stratum and dispatches
`Method → engine`, rejecting invalid `(Method × stratum)` pairs with a clear error.

| Method ↓ \ stratum → | scalar | tangent | ambient |
|---|---|---|---|
| `gradient`   | → tangent (`♯d₀`)           | ✗ covariant derivative (out of scope) | ✗ |
| `divergence` | ✗                           | → scalar (DEC adjoint)               | → scalar (`imag·n` **+ mean-curvature `−2H(J·N)`**) |
| `curl`       | ✗                           | → scalar vorticity (DEC)             | → scalar vorticity (Dirac `w`-part) |
| `laplacian`  | → scalar (LBO)              | → tangent (Hodge **or** connection)  | → ambient (`D²`) |
| `poisson`    | → scalar                    | → tangent                            | → scalar potentials (via helmholtz) |
| `helmholtz`  | ✗                           | (= div + curl)                       | → full Hodge decomposition |

✗ cells raise `bst_operators:badFieldType`, e.g. `"curl is undefined for a scalar field
(stratum=scalar); did you mean gradient?"`. This is the spurious-operation guard.

## 4. Field-type detection

`bst_operators` gains `OPTIONS.FieldType ∈ {'auto','scalar','tangent','ambient'}`
(default `'auto'`). Auto-detection from the row count and realness, given the surface's
`nV` and `nF`:

- `rows == nV` and real → `scalar`
- `rows == nV` and complex → `tangent` (connection representation)
- `rows == 3*nF` → `tangent` (DEC per-face ambient representation)
- `rows == 3*nV` → `ambient`

The one ambiguous case is `nV == nF` (degenerate, essentially never on a real cortex); auto
then errors and demands an explicit `FieldType`. An explicit `FieldType` always overrides
auto-detection.

## 5. Solver infrastructure (new)

### 5.1 `tess_cholesky` (in `toolbox/anatomy`)

Lazy **find-or-compute-or-load** factor helper, mirroring `bst_get_operator_node`. The factor
is the factorization of *one specific assembled operator matrix*, so it is **persisted on the
`operator_` node** (next to the matrix it factorizes), not on the geometry-only `manifold_`
node. Each `operator_` node gains a `Cholesky{hh}` sub-struct.

- **Storage:** `Cholesky{hh} = struct('L',L,'p',p,'free',free,'n',n)`, from
  `[L,~,p] = chol(A(free,free),'lower','vector')`. The opaque `decomposition` object is **not**
  serialized (unreliable across MATLAB versions); `L` (sparse lower factor) and `p`
  (permutation vector) are plain data that save/reload cleanly and reconstruct the two
  triangular solves. `free` records the pinned-index convention; `n` the full dimension.
- **Resolution order:** in-session memory cache (persistent map keyed by operator-node
  identity + pin convention) → factor stored on the `operator_` node on disk → else compute,
  persist to the node, and populate the memory cache.
- **API (provisional):**
  - `dF = tess_cholesky(OperatorNode, hh, PinSpec)` → factor struct for hemisphere `hh`.
  - `x  = tess_cholesky('solve', dF, rhs)` → permuted forward/back-substitution, pinned
    entries zero (recentering is the caller's / `bst_poisson`'s concern).
- `PinSpec` selects the nullspace handling (e.g. pin vertex 1 for the scalar LBO; the
  `[1; nFh+1]` pin for the coupled face system). Different pins under the same operator are
  cached under different keys.

### 5.2 `bst_poisson` (in `toolbox/differential`)

Stratified solver for `L φ = f`, picking the operator by stratum:

- **scalar:** `K φ = M f`, where `K` is the cotan stiffness and `M` the Galerkin mass; project
  the RHS to the mean-zero subspace, pinned solve via `tess_cholesky`, recenter.
- **tangent:** connection-Laplacian solve (complex) when requested; otherwise the de Rham
  route handled by the DEC verbs.
- **ambient:** the two scalar potential solves used inside the Hodge decomposition (stream `ψ`
  from vorticity, potential `φ` from divergence) — i.e. `bst_poisson` is the building block
  `bst_helmholtz` calls, not a separate ambient operator.

The nullspace handling (mean-zero projection → pinned solve → recenter) lives **once** in
`bst_poisson`. This replaces:

- `bst_operators` `i_poisson` (the per-column re-factorization bug — now one cached factor),
- `bst_helmholtz` private `i_poisson` and its `Prepare`-time `decomposition(...)` calls.

## 6. Verb engines (stratified; ambient = wrap the Dirac engine)

The ambient branch reuses the quantities **already** computed by `bst_helmholtz` — no new
physics, per the "wrap the Dirac engine" decision.

- **`bst_gradient`** — `scalar → tangent`, unchanged DEC `grad f = ♯ d₀ f`. Guards non-scalar
  input.
- **`bst_divergence`**:
  - `tangent (DEC face field) → scalar` — unchanged stable adjoint route
    `div V = −⋆₀⁻¹ d₀' ♯' W_F V` (avoids the unstable `⋆₁⁻¹`).
  - `ambient (3nV) → scalar` — the Hodge-decomposition divergence (`imag·n`) **plus the
    explicit mean-curvature term**: `div_Σ(J) = div_Σ(J_tan) − 2H·(J·N)`. `bst_helmholtz`
    uses the **intrinsic** Dirac root, which does not carry this extrinsic coupling; we add
    the `−2H·(J·N)` term explicitly from the shape operator `S` (mean curvature
    `H = ½ tr S`). Chosen over swapping to the extrinsic Dirac root because it is visible and
    independently testable: it must **vanish on a flat patch** and grow with curvature.
- **`bst_curl`**:
  - `tangent → scalar vorticity` — unchanged `curl V = −div(N × V)`.
  - `ambient → scalar vorticity` — the Dirac `w`-part from `bst_helmholtz`.
- **`bst_laplacian`** (new thin engine, or a `bst_operators` branch): `scalar → LBO`;
  `tangent → Hodge or connection` (option, Weitzenböck-documented); `ambient → D²`.
- **`bst_helmholtz`** — refactored onto `bst_poisson` / `tess_cholesky`; behavior bit-identical
  (HarmFrac and all components unchanged within numerical tolerance).

## 7. File-level changes

New:
- `toolbox/anatomy/tess_cholesky.m`
- `toolbox/differential/bst_poisson.m`
- (optional) `toolbox/differential/bst_laplacian.m` — or fold into `bst_operators`.

Modified:
- `toolbox/differential/bst_divergence.m` — add ambient branch + mean-curvature term.
- `toolbox/differential/bst_curl.m` — add ambient branch.
- `toolbox/differential/bst_gradient.m` — add scalar guard.
- `toolbox/differential/bst_helmholtz.m` — refactor onto `bst_poisson` / `tess_cholesky`.
- `toolbox/differential/bst_operators.m` — router, guard table, `FieldType`, drop `i_poisson`.
- `operatormat` DB template (`db_template`) — add `Cholesky` sub-field.
- Docs: Weitzenböck note (Hodge vs connection Laplacian) on the tangent stratum.

## 8. Sequencing (each step independently testable)

1. **`tess_cholesky` + operator_-node persistence** — round-trip solve identical to backslash;
   factor reuse verified (second call hits memory/disk, no refactor); persisted factor reloads
   across a fresh load of the node.
2. **`bst_poisson`** — swap both `i_poisson` sites; results bit-identical to current
   (`bst_operators` `poisson`/`laplacian`, and the `bst_helmholtz` ψ/φ solves).
3. **Refactor `bst_helmholtz`** onto the shared solver — HarmFrac and components unchanged.
4. **Stratified `bst_divergence` / `bst_curl` ambient branch + mean-curvature term** — validate
   on a synthetic ambient field; the `−2H(J·N)` term vanishes on a flat patch and tracks
   curvature on a folded patch.
5. **`bst_operators` router + guard table + `FieldType`** — validate the full operation matrix,
   including every ✗ cell raising the guard error.
6. **Docs** — Weitzenböck note on the tangent stratum.

## 9. Out of scope

- Covariant derivative / Jacobian of a vector field (`gradient` of tangent/ambient).
- Standalone first-class extrinsic grad/div/curl operators (the ambient branch wraps the
  validated Dirac engine instead).
- Replacing the stable DEC tangent div/curl with a connection-Laplacian reimplementation
  (tangent differential ops stay de Rham; the connection Laplacian remains the smoothing/eigen
  Laplacian for the tangent stratum).
- Time-resolved batch driving (stays in `process_` wrappers that loop `bst_operators`).

## 10. Validation

- Unit: `tess_cholesky` solve vs backslash; cached vs fresh factor identical.
- Regression: `bst_operators` `poisson`/`laplacian` and `bst_helmholtz` outputs bit-identical
  before/after the `bst_poisson` refactor.
- Synthetic: mean-curvature divergence term — zero on a flat patch, grows with curvature.
- Guard: each ✗ cell in the operation matrix raises `bst_operators:badFieldType`.
- Identity checks: `curl(grad f) ≈ 0`, `div(curl V) ≈ 0` (metric-free, must hold at any
  curvature).
