# Eigenmode Function Suite — Reorganization Design

**Date:** 2026-05-27
**Branch:** `feature/eigenmode-reorg`
**Status:** Approved design, ready for implementation plan

## Background

The eigenmode-analysis feature added seven functions, all currently in
`toolbox/anatomy/`:

| Function | Role |
|---|---|
| `tess_laplacian` | Cotangent Laplace–Beltrami operator + mass matrix (pure mesh operator) |
| `tess_eigenmodes` | Solve `L φ = λ M φ` for the LBO eigenmodes |
| `tess_eigenmodes_load` | Read precomputed eigenmodes from a surface file |
| `tess_eigenmodes_save` | Write eigenmodes back into a surface file |
| `tess_eigenmodes_project` | Project a scalar field onto the eigenmode basis |
| `tess_eigenmodes_filter` | Spatial-spectral filtering via eigenmodes |
| `tess_eigenmodes_leadfield` | Eigenmode-space inverse imaging kernel |

Five of these are misplaced in `toolbox/anatomy/`: two are file I/O, two are
spectral math, and one is an inverse solver. They also use ad-hoc file access
(`who('-file')`, raw `load`, double `bst_save`) instead of Brainstorm's io and
database conventions.

Separately, the cortex is now downsampled with **icosphere** (FreeSurfer/MNE-style
per-hemisphere) downsampling (`tess_downsize.m`), which produces a clean
2-manifold mesh **without changing vertex identity**. This makes automatic mesh
repair both unnecessary and undesirable.

## Goals

1. Relocate the five misplaced functions to their correct Brainstorm folders and
   rename them to match Brainstorm's prefix conventions.
2. Replace ad-hoc file access with Brainstorm io / database / forward-inverse
   conventions.
3. Make `tess_eigenmodes` / `tess_laplacian` **expect a manifold surface**:
   check manifoldness before computing, and **never repair automatically**.
4. Keep all numerical behavior identical (this is a reorganization, not a
   re-derivation of the math).

## Non-goals (YAGNI)

- Migrating `bst_inverse_eigenmodes` to the full `bst_inverse_linear_2018`
  `OPTIONS`-struct interface (future work).
- Storing the mass matrix `M` on disk (recomputed on demand, as today).
- Hardening every head-model/eigenmode vertex-count guard across the codebase.
- Remapping atlases after a repair (atlas is still cleared on the rare opt-in
  repair path, with a warning).

## Phase A — Verification of the unchanged core (completed)

Both compute functions were reviewed against the icosphere output and need **no
numerical changes**:

- **`tess_laplacian`** is a purely local, per-face assembly (`tess_laplacian.m:157-182`).
  It is agnostic to global topology: two disconnected hemispheres yield a
  block-diagonal `L`/`M`; boundary edges (reducepatch path) get the natural
  Neumann boundary condition; degenerate triangles are guarded by
  `denom(denom==0)=eps` (`:133`). The icosphere guarantees the two properties it
  needs — no zero-area faces and no isolated vertices — so the mass-matrix
  diagonal is strictly positive and `L φ = λ M φ` is well-posed.
- **`tess_eigenmodes`** DC-mode removal, M-orthonormalization, and small-negative
  eigenvalue clamping are all correct for two closed hemispheres.

The only behavioral change to these two functions is the manifold policy below.

## Decisions

1. **Storage:** eigenmodes remain **embedded** as a field in the cortex surface
   `.mat`. `in_tess_bst` does a full `load()` and preserves the `Eigenmodes`
   field untouched (`in_tess_bst.m:48`), so this is fully compatible.
2. **Inverse home:** the leadfield function moves to `toolbox/inverse/` (it
   outputs an `ImagingKernel`, the direct analogue of `bst_inverse_linear_2018`).
3. **Manifold policy:** check, never auto-repair. Repair is strictly opt-in.

## Target layout

| Current (in `toolbox/anatomy/`) | New location | New name |
|---|---|---|
| `tess_laplacian` | `toolbox/anatomy/` (stays) | unchanged |
| `tess_eigenmodes` | `toolbox/anatomy/` (stays) | unchanged |
| `tess_eigenmodes_load` | `toolbox/io/` | `in_tess_eigenmodes` |
| `tess_eigenmodes_save` | `toolbox/io/` | `out_tess_eigenmodes` |
| `tess_eigenmodes_project` | `toolbox/math/` | `bst_eigenmodes_project` |
| `tess_eigenmodes_filter` | `toolbox/math/` | `bst_eigenmodes_filter` |
| `tess_eigenmodes_leadfield` | `toolbox/inverse/` | `bst_inverse_eigenmodes` |

`tess_*` is kept only for the two genuine mesh operators; io functions get
`in_`/`out_`; math/inverse operators get `bst_`.

## Component contracts

### `tess_anatomy/tess_eigenmodes.m` (modified — manifold gate)

Replace the `FixMesh`-default-true auto-repair (`tess_eigenmodes.m:113-116`) with
a **check-only gate**. Rename option `FixMesh` → `Repair`, **default `false`**.

```matlab
% Validate ONLY — never auto-repair. Repair changes vertex/edge counts,
% which desyncs head models, lead fields, and atlases built on this surface.
[~, ~, isManifold, report] = tess_manifold(Vertices, Faces, 'Repair', 0, 'Verbose', Verbose);
if ~isManifold
    if Repair    % explicit opt-in only
        [Vertices, Faces] = tess_manifold(Vertices, Faces, 'Repair', 1, 'Verbose', Verbose);
    else
        error(['Surface is not a clean 2-manifold (%s). Re-mesh with icosphere ' ...
               'downsampling, or pass ''Repair'',true to attempt a risky repair ' ...
               'that changes the vertex count.'], strjoin(report.summary, '; '));
    end
end
```

Rationale: the no-repair default is what keeps the eigenmodes vertex-aligned with
the head model and lead field that were built on the same surface.

### `toolbox/anatomy/tess_laplacian.m` (doc + optional check)

- Document that it **assumes a clean 2-manifold input** (e.g. icosphere cortex).
- Add an optional `'CheckManifold'` flag, **default `false`**, that warns if the
  input is non-manifold. Default off because `tess_laplacian` is recomputed on
  every `project`/`filter` call and `tess_eigenmodes` is the authoritative gate;
  a full manifold scan in that inner loop would be wasteful.

### `toolbox/io/in_tess_eigenmodes.m` (was `tess_eigenmodes_load`)

```matlab
[Eigenmodes, isComputed] = in_tess_eigenmodes(SurfaceFile)
```

- Resolve path via `file_fullpath`; load via `in_tess_bst(SurfaceFile, 0)`
  (`isComputeMissing = 0` so a frequently-called read does not recompute
  curvature/normals — those recomputes are gated on `isComputeMissing` at
  `in_tess_bst.m:128,135,142,149`).
- Extract `.Eigenmodes`; return `isComputed = false` when the field is absent.
- Convert `Vectors` single → double for computation.

### `toolbox/io/out_tess_eigenmodes.m` (was `tess_eigenmodes_save`)

```matlab
out_tess_eigenmodes(SurfaceFile, Eigenmodes, Vertices, Faces, isInteractive)
```

- **Single** `bst_save` (drop the double save at `tess_eigenmodes_save.m:126,132`);
  `bst_history('add', ...)` applied before that single save.
- Store `Vectors` as `single`, `Values` as double column vector (as today).
- Repaired-mesh handling (update `Vertices`/`Faces`, recompute `VertConn`/
  `VertNormals`, clear `Curvature`/`SulciMap`, clear `Atlas` with a warning) is
  preserved but now only ever fires on the rare opt-in repair path.
- The database refresh (`db_reload_subjects`) stays in the process layer, which
  already calls it after `Compute`.

### `toolbox/math/bst_eigenmodes_project.m` (was `tess_eigenmodes_project`)

Made **pure** — no file I/O:

```matlab
[Coeffs, Reconstructed] = bst_eigenmodes_project(Eigenmodes, Data, MassMatrix, 'ModeRange', [k1 k2])
```

- `Eigenmodes`: struct from `in_tess_eigenmodes` (uses `.Vectors`, `.nModes`).
- `MassMatrix`: `[nV x nV]` sparse `M` (computed once by the caller via
  `tess_laplacian`).
- Math is identical: `Coeffs = Phi' * (M * Data)`; optional reconstruction over
  `ModeRange`.

### `toolbox/math/bst_eigenmodes_filter.m` (was `tess_eigenmodes_filter`)

Made **pure** — no file I/O:

```matlab
Filtered = bst_eigenmodes_filter(Eigenmodes, Data, MassMatrix, FilterType, 'CutoffMode', N, ...)
```

- Same filter types and transfer functions (`lowpass`, `highpass`, `bandpass`,
  `heat`, `inverse_heat`, `custom`).
- Math is identical: `Filtered = Phi * (h .* (Phi' * (M * Data)))`.

### `toolbox/inverse/bst_inverse_eigenmodes.m` (was `tess_eigenmodes_leadfield`)

```matlab
[Results, errMsg] = bst_inverse_eigenmodes(HeadModelFile, SurfaceFile, NoiseCovFile, ...)
```

- Algorithm preserved verbatim (lead-field compression `L̃ = L·Φ`, whitening,
  source prior, MNE/dSPM/sLORETA in eigenmode space).
- Only io swap: the `tess_eigenmodes_load` call (`tess_eigenmodes_leadfield.m:172`)
  becomes `in_tess_eigenmodes`. Head model already loads via `in_bst_headmodel`
  and noise cov via `load(file_fullpath(...))` (canonical) — unchanged.

## Process-plugin updates (5 files, the only call sites)

| Plugin | Edit |
|---|---|
| `process_eigenmodes.m` | `:159` load → `in_tess_eigenmodes`; `:226` save → `out_tess_eigenmodes`; manifold-gate + repair opt-in (below) |
| `process_eigenmodes_filter.m` | `:178` load → `in_tess_eigenmodes`; `:231/235` → `bst_eigenmodes_filter` (load `M` once via `tess_laplacian`, pass in) |
| `process_eigenmodes_inverse.m` | `:158` load → `in_tess_eigenmodes`; `:223` → `bst_inverse_eigenmodes` |
| `process_eigenmodes_spectrum.m` | `:151` load → `in_tess_eigenmodes`; `:171` → `bst_eigenmodes_project` (pass `Eigenmodes` + `M`) |
| `process_eigenmodes_view.m` | `:136` load → `in_tess_eigenmodes` |

Plus `SEE ALSO` comment updates in the moved files and the plugins.

### Manifold gate in `process_eigenmodes.m`

- The `fixmesh` checkbox flips meaning and default:
  *"Attempt repair if surface is non-manifold (risky — changes vertex count)"*,
  **default unchecked (`Value = 0`)**.
- `Compute` already loads `Vertices`/`Faces`; add an up-front
  `tess_manifold(..., 'Repair', 0)` check:
  - **Non-manifold + repair off** → `bst_report('Error', ...)` and abort
    (pipeline-safe; never silently reshapes the mesh).
  - **`ComputeInteractive`** → a `java_dialog` warns the surface is non-manifold
    and asks whether to attempt the risky repair; *No/Cancel* aborts, *Yes* sets
    `Repair = true`.
- Pass `Repair` (from checkbox/dialog) into `tess_eigenmodes` instead of `FixMesh`.

### Mass-matrix sharing (perf note)

Because `project`/`filter` are now pure, the process layer computes `M` once via
`tess_laplacian(Vertices, Faces, 'MassType', Eigenmodes.MassType)` and reuses it.
This fixes a real inefficiency: `process_eigenmodes_filter.m:231` currently
recomputes `M` on every frequency iteration.

## Verification plan

1. **M-Lint** (`checkcode`) on every moved/renamed file — no new warnings.
2. **Numerical sanity test** on a synthetic icosphere mesh (`tess_sphere`),
   no database required:
   - `L` symmetric and positive semidefinite (`min(eig) >= -tol`);
   - `M` diagonal strictly positive;
   - eigenvalues ascending and `>= 0`;
   - `UᵀMU ≈ I` (M-orthonormality);
   - correct DC-mode count.
   This empirically closes the Phase-A "works on icosphere" check.
3. **Manifold-policy test:** non-manifold input **errors by default**, mutates
   nothing; only repairs when `Repair = true`.
4. **io round-trip:** `out_tess_eigenmodes` → `in_tess_eigenmodes` returns an
   equivalent struct (single/double handled).
5. **Repo-wide grep:** zero remaining references to the five old function names
   across `toolbox/`, `dev/`, and `toolbox/script/`.
6. *Optional, if a protocol is available:* re-run
   `dev/tests/test_omega_icosphere_sourcemap.m`.

## Risks

- **Behavioral change for existing callers** relying on auto-repair: the default
  now errors on non-manifold input. Mitigated by the icosphere workflow (clean by
  construction) and the explicit `Repair` opt-in. This is intentional.
- **Renamed functions** break any out-of-tree scripts referencing the old names.
  Mitigated by the repo-wide grep (step 5); only the five process plugins call
  them in-tree.
