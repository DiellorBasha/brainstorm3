# nxr-compute v0.2.0 upgrade + operator-registry adoption

**Date:** 2026-06-26   **Branch:** `feature/nxr-v020-registry`   **Status:** IMPLEMENTED + live-validated 2026-06-26

## Implementation status (2026-06-26)

Implemented on `feature/nxr-v020-registry` and live-validated against the managed
v0.2.0 install on the `preventad` protocol (ico5 canonical cortex). All six
operator variants stamp the correct registry primary id
(`laplaceBeltrami`, `leviCivitaConnectionLaplacian`, `relativeDirac`,
`relativeFaceDirac`, `faceLaplacianGreenGauss`, `flatCovariantLaplacian`).
Validation scripts (all green): `dev/tests/test_nxr_v020_smoke.m`,
`test_bst_nxr_registry.m`, `test_operatormat_registry_field.m`,
`test_tess_operators_registry.m`, `test_tess_eigen_registry_check.m`,
`test_nxr_v020_registry_e2e.m`.

Notes from validation:
- ⚑ The published v0.2.0 `nxr_compute('version')` still returns the string
  `"nxr-compute 0.1.0"` (hardcoded, not bumped upstream). Functionally v0.2.0 —
  `operatorInfo`/`fieldInfo` (PR #5) are present. Cosmetic; bump in a future nxr release.
- ⚑ Pre-existing (NOT this project): `tess_operators` `'Hodge-Face'` (the
  `TessMat.VertNormals(vH(Floc(:,1)),:)` line) crashes with a cryptic index error
  when the surface has empty `VertNormals` (e.g. the high-res `tess_cortex_*_high`
  surfaces). Tests use `bst_canonical_cortex(20484)` (ico5, has VertNormals). Worth
  a defensive guard in a separate change.

## Goal

Upgrade the Brainstorm `nxr-compute` plugin from the manually-copied v0.1.0
build to a **managed install of the published v0.2.0 release**, and adopt the
new **operator registry** (controlled-vocabulary metadata) so our operator DB
nodes carry, and our consumers can reason from, each operator's type/role.

User-selected scope: **Upgrade + adopt registry** (not a package migration).

## Key finding — v0.2.0 is additive, not a breaking rewrite

Verified directly against the published v0.2.0 source
(`bindings/mex/src/nxr_compute_mex.cpp`, `src/operator_registry.cpp`):

- The release advertises a new `+nxr` MATLAB package
  (`nxr.manifold.context`, `nxr.manifold.operator.*`, …). That package is a
  **thin shim over the same flat `nxr_compute('cmd', …)` dispatcher** we
  already use — no new compute, and it does NOT expose all the operators we
  rely on (our face-Dirac stack is dispatcher-only). We keep the flat
  dispatcher; **no consumer call sites change calling convention.**
- Every command + operator literal our consumers pass is present **unchanged**
  in v0.2.0: `create`, `destroy`, `operators` (`laplacian/cotan`,
  `mass/galerkin`, `mass/lumped`, `laplacian/connection`, `dirac`, `diracD`,
  `diracIntrinsicD`, `diracFaceIntrinsicD`, `diracFace`, `gradFace`, `lapFace`,
  `dec`, `hodge/h0|h1|h2`), `gauge/levi-civita`, `vertexFrames`, `facets`,
  `heat`, `tracePath`, `version`. v0.2.0 even adds `diracFaceD`, `gradient3D`,
  `connectionGradient`.
- The MEX entry point is still `nxr_compute` (binary names
  `nxr_compute.mex{maca64,w64,a64}` unchanged) — so `TestFile` does not change.

**The genuinely new surface** (PR #5, "Operator registry +
controlled-vocabulary metadata"):

- `nxr_compute('operatorInfo', id)` → struct fields: `id, label, bundle,
  holonomy, order, role, field_type, domain, singular, gauge, coupling,
  natural_mass, graded (double 0/1), tau_presets, status, notes, squares_to,
  square_of, relation, input_field, output_field`.
- `nxr_compute('fieldInfo', id)` → `id, label, domain, bundle, field_type,
  n_form, representation, gauge, nSym (double), notes`.

Registry operator ids (controlled vocabulary):
`laplaceBeltrami, graphLaplacian, faceLaplacianGreenGauss, faceLaplacian2Form,
leviCivitaConnectionLaplacian, trivialConnectionLaplacian,
leviCivitaConnectionGradient, flatCovariantLaplacian, productCovariantLaplacian,
covariantGradient, faceGradient, extrinsicWeitzenbockLaplacian, intrinsicDirac,
extrinsicDirac, relativeDirac, intrinsicFaceDirac, extrinsicFaceDirac,
relativeFaceDirac, massLumped, massGalerkin, d0, d1, hodge0, hodge1, hodge2,
hodge1inv`.

---

## Part A — Plugin upgrade (manual copy → managed v0.2.0)

`toolbox/core/bst_plugin.m`, the `nxr-compute` PlugDesc block (~lines 192–221):

- `PlugDesc(end).Version` : `'0.1.0'` → `'0.2.0'`.
- `nxrRel` base URL and both track zip names → v0.2.0:
  - `https://github.com/neurodynamics-xr/nxr-compute/releases/download/v0.2.0/`
  - `nxr-compute-mex-r2023b-v0.2.0.zip` (R2023b+ / macOS Apple Silicon)
  - `nxr-compute-mex-r2023a-v0.2.0.zip` (R2023a; Windows + Linux)
- `TestFile` per OS: **unchanged** (`nxr_compute.mexmaca64` / `.mexw64` /
  `.mexa64`).
- `CompiledStatus=1`, `AutoUpdate=0`, `AutoLoad=0`: **unchanged**
  (install-on-demand). Keep the existing R2023b-only-for-mac comment, retargeted
  to v0.2.0.

**Migrate off the manual copy.** Current install is a hand-placed
`~/.brainstorm/plugins/nxr-compute/nxr-compute-mex-r2023b/` + hand-made
`plugin.mat`. Steps (interactive, MATLAB):

1. `bst_plugin('Uninstall', 'nxr-compute')` (or delete the dir) and **verify the
   directory is gone** — the documented "managed-plugin stale-binary trap": if
   the old binary lingers, BST loads it instead of the download.
2. `bst_plugin('Install', 'nxr-compute')` → downloads + unzips the v0.2.0
   release, writes a proper managed `plugin.mat`. The v0.2.0 binary is a
   **superset** of the local build (contains all custom Dirac/face operators).
3. Confirm the extracted folder (`nxr-compute-mex-r2023b/`) is on path and the
   `TestFile` resolves.

**Regression gate (must pass before Part B lands):** on the canonical cortex
(`bst_canonical_cortex(20484)` + `nxr_safe_create`), confirm `nxr_compute('version')`
== v0.2.0 and that each command/operator above returns a finite-sized result:
`create` → all `operators` variants + `gauge/levi-civita` + `vertexFrames` +
`facets` + `heat` + `tracePath` → `destroy`. This is the backward-compat proof.

---

## Part B — Registry adoption

### B1. New helper `toolbox/anatomy/bst_nxr_registry.m`

Single source of truth, two jobs:

1. Guarded wrappers:
   - `meta = bst_nxr_registry('operator', id)` → `operatorInfo` struct, or `[]`
     on any error (pre-registry binary, unknown id). Never throws.
   - `meta = bst_nxr_registry('field', id)` → `fieldInfo` struct, or `[]`.
2. The **Variant → primary registry id** map (defined here, nowhere else):

   | Brainstorm Variant     | primary registry id              | components recorded |
   |------------------------|----------------------------------|---------------------|
   | Laplace-Beltrami       | `laplaceBeltrami`                | massGalerkin        |
   | Connection Laplacian   | `leviCivitaConnectionLaplacian`  | massGalerkin        |
   | Dirac                  | `relativeDirac`                  | intrinsicDirac, extrinsicDirac, massGalerkin |
   | Dirac-Face             | `relativeFaceDirac`              | intrinsicFaceDirac, extrinsicFaceDirac, massLumped |
   | Hodge-Face             | `faceLaplacianGreenGauss`        | faceGradient        |
   | Covariant              | `flatCovariantLaplacian`         | laplaceBeltrami, massGalerkin |

   Decision (per user): primary id = the **assembled** operator `A`
   (`relativeDirac`/`relativeFaceDirac` for the τ-blends), and the blend
   **components** are also recorded so the node documents what it was built from.

   Helper accessors:
   - `id = bst_nxr_registry('idForVariant', Variant)`
   - `ids = bst_nxr_registry('componentsForVariant', Variant)`

### B2. `operatormat` DB template

`toolbox/db/db_template.m`, `case 'operatormat'`: add a `Registry` field
(sibling of `Provenance`), default `[]`. Holds `struct('Primary', <operatorInfo
struct>, 'Components', <1xN struct array of operatorInfo>)`. Additive: existing
nodes load with `Registry=[]`.

### B3. `tess_operators.m`

After the per-hemisphere build loop, once `Variant` and `Tau` are known:

- `pid = bst_nxr_registry('idForVariant', Variant)`;
  `OperatorMat.Registry.Primary = bst_nxr_registry('operator', pid)`.
- `OperatorMat.Registry.Components` = `operatorInfo` for each component id.
- All guarded — if the registry is unavailable, `Registry` stays `[]` and the
  node is built exactly as today.
- **Optional name-validation:** before each `nxr_compute('operators', …)` family
  call, the corresponding registry id is known; if `operatorInfo(id)` returns
  empty on a binary that DOES support the registry (i.e. version >= v0.2.0 but id
  missing), warn once with a clear message. Do not hard-fail (keeps older
  binaries working).

### B4. Consume the metadata where it removes a hardcoded assumption

Surgical only — do NOT rewrite the `Variant` switch:

- In `tess_eigen` / `bst_dirac`, the "is this a 4-component (quaternion)
  operator?" / "vertex vs face domain?" decisions currently re-derive from
  `Variant`. Where a node has `Registry.Primary`, read `graded` /
  `field_type` (`quaternion`) / `domain` (`vertex`|`face`) instead, falling back
  to the existing `Variant` switch when `Registry` is `[]`.
- Target the one or two concrete spots only; leave the rest untouched.

---

## Part C — Validation & docs

- `dev/tests/` (or `dev/benchmarks/`) smoke script: every Variant end-to-end on
  the canonical cortex against the managed v0.2.0 binary; assert
  `OperatorMat.Registry.Primary.id` matches the expected primary id and the
  operator/mass shapes are unchanged vs a pre-upgrade baseline.
- `operatorInfo`/`fieldInfo` round-trip check for each primary id.
- Update memory: `nxr-compute-plugin` (managed v0.2.0, registry), and a note on
  the registry in `nxr-bundle-surface-fields`.
- This design doc is the implementation gate.

## Out of scope (YAGNI)

- No migration to the `nxr.manifold.*` package (cosmetic, partial, no gain).
- No rewrite of the `Variant` dispatch switch in `tess_operators`.
- No changes to WASM / Node / CLI bindings (upstream repo, separate concern).

## Risks / watch-items

- **Stale-binary trap** (Part A step 1): the #1 failure mode — old manual MEX
  shadows the download. Verify removal explicitly.
- **Mac Intel:** v0.2.0 ships no prebuilt Intel-mac MEX (Apple-Silicon only on
  the R2023b track). Not our machine, but note for the cluster/fleet.
- **Cluster (Linux) parity:** the managed install fetches the Linux MEX from the
  same release — resolves the earlier "rebuild nxr Linux MEX" TODO for v0.2.0,
  provided the release's `mexa64` loads on the cluster's glibc.
