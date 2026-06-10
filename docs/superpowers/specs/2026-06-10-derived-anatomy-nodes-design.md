# Derived-Anatomy Nodes — Surface-Child `manifold_`/`eigen_`/`operator_` Files — Design

**Date:** 2026-06-10
**Status:** Design (pending review)
**Repo:** `brainstorm3` (MATLAB DB + tree) — and, only if a generic Java node type can't be reused, the `bst-java` fork (`/Users/diellorbasha/workspace/research/code/bst-java`).

## Goal

Replace the tangle of ad-hoc structural fields on the cortical `TessMat` with a
principled set of **derived-anatomy nodes** — separate files nested as children of
a surface in the database tree, mirroring how a study owns nested
Channel/HeadModel/Result nodes. Three node types:

- **`manifold_`** — the lightweight geometry backbone (Topology, Geometry, Gauge,
  Embedding from nxr-compute), assembled by `tess_manifold`.
- **`eigen_`** — eigenbases (`Scalar`=Laplace-Beltrami, `Vector2`=connection
  Laplacian, `Vector3`=Dirac), variant by label; assembled by `tess_eigen`.
- **`operator_`** — operators (Mass, DEC, Laplacian, Dirac), variant by label;
  assembled by `tess_operators`.

Path references (a head model already references its surface by path) and lazy
loading fall out naturally; context menus attach to the specific node type, so
"View tangent basis"/"View eigenmodes" appear only on surfaces that actually have
the relevant node — not on plain surfaces (head mask, `mid_`, …).

**This spec covers SP1: the node infrastructure (new DB data types + tree nodes),
proven end-to-end with the `manifold_` node.** `eigen`/`operator` reuse the same
pattern in later sub-projects.

## Background / problem (verified)

The cortex `TessMat` currently carries overlapping ad-hoc fields — `Eigenmodes`
(LBO), `ConnEigenmodes`, `DiracEigen`, `Operators`, `Coordinates` — plus
`nxr`/`TangentFrame` written elsewhere, from ~8 writers
(`tess_eigenmodes`/`tess_conn_eigenmodes`/`tess_dirac_eigenmodes`/`tess_laplacian`/
`tess_connection_laplacian`/`tess_nxr_populate`/`tess_frame`/`tess_tangents`) and
`bst_*_ensure` + `in_/out_tess_*eigenmodes` IO. Eigenbases are large
(`[4V×K]`≈130 MB/hemi at K=400), so storing everything in TessMat bloats the
surface file. No principled Topology/Geometry/Gauge/Embedding backbone exists.

## Decisions (resolved during brainstorming)

- **Storage = dedicated file-type nodes** nested under the surface (option B),
  using Brainstorm's study→children pattern (filename-keyed entries; parent lists
  its children; nested in the tree).
- **Geometry backbone = a `manifold_` node** (NOT default TessMat fields) — only
  the surfaces that need it get one; added as an analysis result on demand.
- **Tree placement = nested children of the surface.**
- **Java node type = reuse a generic existing type** with a custom icon if
  possible; only build in the `bst-java` fork if forced (never push/merge to
  upstream).
- **Breaking changes are accepted this round** — downstream consumers
  (`tess_frame`/`tess_tangents`/viewers/etc.) are NOT migrated here; the workflow
  is being redesigned and migration lands in SP4.
- **Rename** the existing `tess_manifold` (2-manifold validator/repair) →
  `tess_repair`; do not fix its callers this round.

## Decomposition (each its own spec→plan)

- **SP1 (this spec):** derived-node DB+tree infrastructure + the `manifold_` node
  end-to-end (`tess_manifold` assembler). The backbone + the proven pattern.
- **SP2:** `eigen_` node (`tess_eigen`, variants Scalar/Vector2/Vector3) + migrate
  the eigenmode leadfield/inverse consumers to read the eigen node by path.
- **SP3:** `operator_` node (`tess_operators`, variants Mass/DEC/Laplacian/Dirac).
- **SP4:** retire old TessMat fields + old writers + `in_/out_tess_*eigenmodes`;
  redesign/migrate the frame/tangent consumers; migration for existing files.

---

## SP1 design

### A. Rename the validator
`toolbox/anatomy/tess_manifold.m` (2-manifold validate/repair) → `tess_repair.m`
(function renamed; help/identifiers updated). Its 7 callers
(`tess_laplacian`/`tess_connection_laplacian`/`tess_eigenmodes`/
`tess_conn_eigenmodes`/`bst_eigenmodes_ensure`/`bst_conn_eigenmodes_ensure`/
`process_eigenmodes`) will break — **left broken intentionally** (redesigned in
SP4). This frees `tess_manifold` for the assembler.

### B. DB data types (`toolbox/db/db_template.m`)
- Extend `db_template('surface')` with typed child lists (default empty), mirroring
  `sStudy.Channel`/`.HeadModel`/`.Result`:
  ```
  .Manifold : struct('FileName',{},'Comment',{})
  .Eigen    : struct('FileName',{},'Comment',{},'Variant',{})   % wired now; populated SP2
  .Operator : struct('FileName',{},'Comment',{},'Variant',{})   % wired now; populated SP3
  ```
- New `db_template('manifoldmat')` (file content):
  ```
  Comment, ParentSurface,            % back-reference (path)
  Topology, Embedded, Intrinsic, Extrinsic, Gauge,  % nxr facet bundle, 1x2 per-hemi
  Provenance                          % nxr version, gauge, date
  ```
  (i.e. the nxr `facets` bundle = "Topology / Geometry(=Intrinsic+Extrinsic) /
  Gauge / Embedding(=Embedded)".)
- Files named `manifold_*.mat` in the subject anat folder (next to the surface);
  filename is the DB key; `ParentSurface` stored inside for resolution.

### C. Registration + resolution
- `toolbox/db/db_add_manifold.m` — `db_add_manifold(iSubject, ParentSurfaceFile,
  ManifoldStruct, Comment)`: save `manifold_*.mat`, append an entry to the parent
  surface's `.Manifold` list in `sSubject`, `bst_set('ProtocolSubjects', …)`.
- `bst_get` hooks: `bst_get('ManifoldFile', FileName)` → `[sSubject, iSubject,
  iSurface, iManifold]`; a helper to list a surface's manifold node(s). Parent
  surface resolved from the manifold entry's position (it lives under a surface)
  and/or the file's `ParentSurface`.

### D. Frontend (tree)
- `toolbox/tree/node_create_subject.m`: in the surface loop, after `nodeSurface`
  is created, iterate `sSubject.Surface(iSurface).Manifold` (and `.Eigen`/
  `.Operator`) and `nodeSurface.add(childNode)` for each — nested children.
- **Java node type:** first implementation step is to enumerate the existing
  `BstJava` node types and pick a reusable generic one (custom icon + comment +
  leaf-with-context-menu). If none fits, add the type in the `bst-java` fork
  (`dev`/feature branch; never upstream) and rebuild. Document which path was taken.
- `toolbox/tree/tree_callbacks.m`: surface context menu **"Compute manifold"**
  (→ `tess_manifold`); `manifold_` node menu **"View"** (a viewer hook — minimal in
  SP1, e.g. print/inspect; rich viz later) and **"Delete"** (file_delete + remove
  the `.Manifold` entry + reload).

### E. Assembler `tess_manifold(SurfaceFile, ...)` (new)
Pulls the nxr facet bundle for the surface and writes a `manifold_` node:
1. `in_tess_bst`, guards (nxr-compute installed; Structures L/R atlas;
   `Reg.Sphere` for the trivial gauge) — same guards as the retired `tess_frame`.
2. `tess_hemisplit` (atlas L/R, never `conncomp`); per hemisphere
   `nxr_compute('facets', h, 'trivial', opts)` (FreeSurfer-pole singularities),
   attach scatter maps → 1x2 per-hemi `Topology/Embedded/Intrinsic/Extrinsic/Gauge`.
3. Build the `manifoldmat` struct (+ `ParentSurface`, `Provenance`); `db_add_manifold`.
Options: `Gauge` (trivial), `NoSave`, `ForceRecompute`.

### F. Validation / tests
- `db_template('surface')` returns the new child-list fields; `db_template('manifoldmat')`
  returns the schema (DB-free unit checks).
- `db_add_manifold` on a live protocol creates a `manifold_*.mat`, the surface's
  `.Manifold` list gains an entry, and `bst_get('ManifoldFile', …)` resolves it.
- `tess_manifold` on the real 20484 cortex produces a node whose
  `Embedded.vertex.grid` is orthonormal and `Gauge` Gauss-Bonnet holds (reuse the
  facet checks), and the node appears nested under the cortex in the tree
  (visual + programmatic `sSubject.Surface(iCortex).Manifold` non-empty).
- Tree rendering / icon / menus: visual verification in the GUI.

## Edge cases / risks

- **Java node type** is the main unknown — mitigated by reusing a generic type
  first; the `bst-java` fork is the fallback.
- **`bst_set('ProtocolSubjects')` round-trip** must preserve the new `.Manifold`
  child list across save/reload (`db_save`/`db_reload`); covered by a test.
- **Default-anatomy / shared surfaces:** a `manifold_` node attaches to the
  specific surface file; if surfaces are shared via default anat, the node lives
  with the surface's subject (same as other anat). Confirm during implementation.

## Out of scope (SP1)

- `eigen_`/`operator_` **assemblers** and their `*mat` schemas (SP2/SP3; the child
  lists + nesting are wired now).
- **Consumer migration** (`tess_frame`/`tess_tangents`/viewers/leadfield/inverse) —
  SP4. Breaking changes accepted now.
- Rich node visualizations (eigenbasis browser, frame/topology viewers) — later.
- Retiring old TessMat fields / functions / IO — SP4.
