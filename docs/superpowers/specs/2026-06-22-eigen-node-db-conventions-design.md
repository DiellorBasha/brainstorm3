# Design: DB-convention-correct find / overwrite for manifold·operator·eigen nodes

**Date:** 2026-06-22
**Status:** Approved (design); spec under review before implementation.
**Branch:** `refactor/file-based-dirac-consolidation` (continues Increments 1–2).
**Relates to:** `dev/2026-06-22-eigen-spectral-consolidation.md` (Increment 3).

## Goal

Make the derived-anatomy node types (`manifold_`, `operator_`, `eigen_`) behave like
first-class Brainstorm DB citizens for *finding* and *overwriting*, reusing the existing
DB machinery rather than hand-rolling logic in `tess_*`. Concretely:

1. **Consolidate** the duplicated find-logic (`tess_eigen.local_find_eigen` /
   `local_find_operator`, `bst_dirac.local_find_dirac_eigen`).
2. **Fix** three caveats: (a) operator-find ignores Tau; (b) duplicate operators
   accumulate with no GC; (c) a raw `load(f,'Variant')` bypasses the loaders.
3. **Add** an interactive **Overwrite / Cancel** confirmation when a compute would
   replace an existing node.

## Conventions this design follows (established from the codebase)

| Concern | Function family | Reference |
|---|---|---|
| Find by criteria | `bst_get` (cache-only) | `bst_get('SurfaceFileByType')`; `in_bst_results` refines by loading |
| Load a file | `in_bst_*` (pure read) | `in_bst_eigen/operator/manifold` (Increments 1–2) |
| Append a node | `db_add_*` | `db_add_eigen/operator/manifold` |
| Replace with prompt | `db_set_*` (`ReplaceFile=[ask]`) | `db_set_noisecov` |
| Delete a node | `node_delete` (`isUserConfirm` → `java_dialog`) | `node_delete.m` |
| Schema + migration | `db_template` + `db_update` | v5.04 added the Manifold/Eigen/Operator lists |

`bst_get` never touches the filesystem, so for it to resolve a full spec the
discriminating params must live in the cache entry (decision below).

## Decisions (agreed)

- **Dialog:** **Overwrite / Cancel** (two-way). Cancel ⇒ reuse the existing node; Overwrite
  ⇒ delete the existing node and recompute fresh.
- **Trigger:** **exact-spec match.** A non-matching request (different nModes/Tau/Gauge) creates a
  new node *without* prompting.
- **Operator dependents:** **warn + cascade delete.** Overwriting an operator also deletes
  the eigen nodes that reference it; the confirm message states the count.
- **Find mechanism:** **store the spec in the cache entry** so `bst_get` resolves fully from
  cache (no file loading during a find).

## Current state (verified)

`db_template('surface')` child entries (lines 36–39):
```
Manifold : struct('FileName','Comment')                       % + Gauge
Eigen    : struct('FileName','Comment','Variant')             % + nModes, Tau, OperatorFile
Operator : struct('FileName','Comment','Variant')             % + Tau
```
`db_update` is at **v5.06**; `NormalizeSurfaceArray` backfills missing *template* fields
(empty). `db_add_eigen` sets `newEntry.{FileName,Comment,Variant}` only.

---

## Part 3a — Cache the spec; let `bst_get` find it

**Field naming.** All new fields use clear, descriptive names. The ambiguous `K` is
**renamed to `nModes`** (modes per hemisphere) — which also matches the existing public
option name (`bst_dirac(..., 'nModes', …)`, `bst_eigenmodes_ensure`). This rename applies to
**both** the new cache entry **and** the on-disk `EigenMat` field (see "Field rename" below).
`Tau` (the documented relative-Dirac curvature-mixing weight) and `Lambda` (eigenvalues) are
kept — they are domain-standard, documented names, not opaque letters.

**Schema (`db_template('surface')`)** — final entry shapes (Eigen carries `OperatorFile`
so the operator→eigen cascade lookup is also cache-only, see below):
```
Eigen    : struct('FileName','Comment','Variant','nModes','Tau','OperatorFile')
Operator : struct('FileName','Comment','Variant','Tau')
Manifold : struct('FileName','Comment','Gauge')
```

**`db_add_*` populate the new fields** (single source of truth = the saved struct):
- `db_add_eigen`:    `newEntry.nModes = EigenMat.nModes`; `newEntry.Tau = EigenMat.Provenance.Tau` (or `[]` for non-Dirac variants); `newEntry.OperatorFile = EigenMat.OperatorFile`.
- `db_add_operator`: `newEntry.Tau = OperatorMat.Provenance.Tau` (or `[]`).
- `db_add_manifold`: `newEntry.Gauge = ManifoldMat.Provenance.Gauge` (or `''`).

**Field rename `EigenMat.K` → `EigenMat.nModes`.** Update `db_template('eigenmat')`,
`tess_eigen`, `bst_dirac`, and any other consumers. Existing `eigen_*.mat` files on disk
carry the legacy `.K`; `in_bst_eigen` gains a backward-compat shim (if `.K` present and
`.nModes` absent, set `.nModes = .K`) — exactly the legacy-field handling pattern of
`in_bst_headmodel`. So old files keep loading without a file rewrite.

**Migration `db_update` v5.07 — value-backfill.** After `NormalizeSurfaceArray` adds the new
(empty) fields, a v5.07 step walks every subject's surfaces and, for each Eigen/Operator/
Manifold entry whose new field is empty, loads the file once via `in_bst_*` and fills
`nModes`/`Tau`/`Gauge`/`OperatorFile` from the struct / `Provenance`. Stale entries (file missing) are left empty
and skipped (cannot match a spec → effectively absent, which is correct).

**New `bst_get` by-spec cases** (mirror `SurfaceFileByType`; cache-only; return the
existing 4-tuple shape `[sSubject, iSubject, iSurface, iNode]`, empties if none):
- `bst_get('EigenFileForSurface',    SurfaceFile, Variant, nModes, Tau)` — match: `Variant` equal, entry `nModes >= nModes` (requested), and (Dirac/Dirac-Face/Hodge-Face) `Tau` equal; Tau ignored for Laplace-Beltrami / Connection Laplacian. Among matches, return the one with the smallest sufficient `nModes`.
- `bst_get('OperatorFileForSurface', SurfaceFile, Variant, Tau)` — match: `Variant` equal and (Dirac-type) `Tau` equal.
- `bst_get('ManifoldFileForSurface', SurfaceFile, Gauge)` — match: `Gauge` equal.

Reverse dependency lookup for the cascade is cache-only because the Eigen entry now carries
`OperatorFile` (above): `bst_get('EigenFilesForOperator', SurfaceFile, OperatorFile)` scans the
surface's Eigen entries and returns those whose `OperatorFile` matches. (Backfilled by the same
v5.07 step.)

## Part 3b — Consolidate finders via the new `bst_get` cases

- `tess_eigen`: replace `local_find_eigen` with `bst_get('EigenFileForSurface', …)` +
  `in_bst_eigen` load; replace `local_find_operator` with
  `bst_get('OperatorFileForSurface', …, Tau)` + `in_bst_operator` load (now **Tau-aware** —
  caveat 1 — and no raw `load(f,'Variant')` — caveat 3).
- `bst_dirac`: replace `local_find_dirac_eigen` with `bst_get('EigenFileForSurface', …)`.
- Delete the now-dead local helpers. Behavior preserved except the Tau fix.

## Part 3c — Make the nodes deletable (`node_delete`)

`node_delete.m` has no case for `manifold`/`operator`/`eigen` (gap: they're not deletable via
the tree today). Add cases that: confirm (when `isUserConfirm`), `file_delete` the node file,
remove the entry from `sSubject.Surface(iSurface).<Type>`, persist (`bst_set` +
`db_save`), and refresh the tree (`panel_protocols('UpdateNode', …)`). The exact tree node
*type string* (from `node_create_subject`) must be confirmed during implementation.

## Part 3d — Interactive Overwrite / Cancel

- Add `'Interactive', 0|1` (default **0**) to `tess_eigen` / `tess_operators` /
  `tess_manifold`. The tree "Compute …" menu callbacks pass `'Interactive', 1`
  (`tree_callbacks.m:1177/1183-1185/1191-1193`). All programmatic callers keep the default →
  **no prompt, current behavior** (automation-safe).
- When `Interactive=1` and the relevant `bst_get(...ForSurface...)` returns an exact-spec
  match:
  - `java_dialog('confirm', msg, title)` → **Overwrite / Cancel**.
  - **Cancel** ⇒ return the existing node (reuse), no compute.
  - **Overwrite** ⇒ delete the matched node via the 3c `node_delete` path, then
    recompute + `db_add_*` (so no duplicate accumulates — caveat 2 on the interactive path).
- **Operator cascade:** before overwriting an operator, `bst_get('EigenFilesForOperator', …)`
  lists dependents; the confirm message appends *"This will also remove N dependent
  eigenbasis node(s)."*; Overwrite deletes them via `node_delete` first.
- The prompt covers only the *primary* artifact the user asked to compute. Dependencies
  (e.g. the operator a `tess_eigen` compute needs) are found-or-created silently.

## Data integrity / edge cases

- **Nested basis vs overwrite:** a request for `nModes=400` matching an existing `nModes=500`
  Dirac/Tau=0.5 node counts as a match (reuse on Cancel). The confirm message shows the
  existing node's actual `nModes`/Tau so Overwrite (which would replace 500 with 400) is an
  informed choice.
- **Non-Dirac Tau:** Laplace-Beltrami / Connection Laplacian store `Tau=[]`; their find ignores
  Tau.
- **Stale entries** (file missing): never match a spec; `node_delete`/overwrite tolerate a
  missing file (remove the entry regardless).
- **Reuse picks the smallest sufficient `nModes`** so a tiny request doesn't load a huge basis
  when a right-sized one exists.

## Validation

- Extend `dev/tests/test_in_bst_nodes_pure.m` (or a sibling) — but the by-spec finders need a
  registered protocol, so add a live-protocol smoke (in MATLAB MCP): create two Dirac eigen
  nodes at different Tau, assert `bst_get('EigenFileForSurface', …, Tau)` discriminates;
  assert operator Tau discrimination; assert `EigenFilesForOperator` returns dependents.
- `checkcode` clean on every changed file.
- Migration: load the TutorialAuditory protocol (pre-v5.07) and confirm entries backfill
  nModes/Tau/Gauge/OperatorFile without error; re-running is a no-op.
- Overwrite flow: GUI path validated by code review + a non-interactive equivalence check
  (Interactive=0 unchanged); a scripted Overwrite=Yes path test deletes old + creates fresh.

## Out of scope

- The legacy scalar-LBO retirement (Increments 4+).
- The `bst_eigen_spectral_*` umbrella.
- Changing `tess_operators`' default (non-interactive) "always recompute" behavior.

## Sub-commit sequence

1. **3a** — schema fields + `db_add_*` populate + `db_update` v5.07 migration + `bst_get` by-spec cases (+ `EigenFilesForOperator`).
2. **3b** — refactor `tess_eigen` / `bst_dirac` finders onto `bst_get`; delete dead local helpers (Tau fix + caveat 3).
3. **3c** — `node_delete` support for `manifold`/`operator`/`eigen`.
4. **3d** — `Interactive` option + Overwrite/Cancel + operator cascade; menu callbacks pass `Interactive=1`.
