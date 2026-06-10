# Derived-Anatomy Nodes — SP1 (infrastructure + `manifold_` node) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Establish Brainstorm's first **derived-anatomy node** mechanism — files nested as children of a surface in the DB tree — proven end-to-end with the `manifold_` node (the nxr geometry backbone), assembled by `tess_manifold`. `eigen`/`operator` reuse this pattern in later sub-projects.

**Architecture:** Mirror the study→children DB pattern: `db_template('surface')` gains typed child lists (`.Manifold`/`.Eigen`/`.Operator`); a `manifold_*.mat` file (schema `db_template('manifoldmat')`) is registered under its parent surface by `db_add_manifold`; `node_create_subject` nests it; `tree_callbacks` gives it menus; `bst_get('ManifoldFile', …)` resolves it. The tree node uses a dedicated `'manifold'` type string — rendered with a default icon by the unmodified Java `BstNode` (graceful unknown-type fallback verified), behavior dispatched MATLAB-side in `tree_callbacks`. No Java build in SP1 (a custom icon is optional later work in the `bst-java` fork — `/Users/diellorbasha/workspace/research/code/bst-java`, never push upstream).

**Tech Stack:** MATLAB (Brainstorm DB/tree), nxr-compute MEX (`facets`), MATLAB MCP for live DB/tree verification.

**Prereqs verified:** nxr `facets` installed; `BstNode` unknown-type fallback graceful; `tree_callbacks` dispatches on `lower(nodeType)`; `bst_get` at `toolbox/core/bst_get.m`.

**Reference patterns (read before each task):**
- `toolbox/db/db_add_surface.m` — db_add pattern (db_template, file_short, `bst_get/bst_set('ProtocolSubjects')`, `panel_protocols('UpdateNode')`).
- `toolbox/db/db_template.m` — `'surface'` (32), `'study'` child arrays (371+), how `*mat` schemas are defined.
- `toolbox/tree/node_create_subject.m` — surface loop (110-131), `CreateNode` (144).
- `toolbox/tree/tree_callbacks.m` — surface context menu (search `case 'cortex'`/`'surface'`) and node dispatch (`switch lower(nodeType)`).
- `toolbox/anatomy/tess_frame.m` — the facet-bundle assembly (`in_tess_bst` → atlas guard → `tess_hemisplit` → per-hemi `nxr_compute('facets',…)` → scatter maps) to reuse in `tess_manifold`.
- `toolbox/io/file_gettype.m` — registering a new `manifold_` file type.

**Session discipline:** run via the MATLAB MCP; only ever `clear <named functions>`, never bare `clear`. DB-mutating verifications must restore the protocol (delete test nodes via `file_delete` + `db_reload_studies`/`db_reload_database('current')`).

---

## Task 1: Free the name + de-risk the node type

**Files:** rename `toolbox/anatomy/tess_manifold.m` → `toolbox/anatomy/tess_repair.m`.

- [ ] **Step 1: Rename the validator** (accept broken callers; not fixed this round)
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git mv toolbox/anatomy/tess_manifold.m toolbox/anatomy/tess_repair.m
```
Edit the file: change the function signature `function [...] = tess_manifold(` → `tess_repair(`, the header `% TESS_MANIFOLD:` → `% TESS_REPAIR:`, and any `error('tess_manifold:...')` identifiers → `tess_repair:`. Do NOT touch its 7 callers (left broken intentionally; SP4 redesigns them).

- [ ] **Step 2: De-risk — confirm a dedicated `'manifold'` BstNode type renders without Java** (live)
Via the MATLAB MCP, build a throwaway tree node with the new type string and confirm no error + it renders:
```matlab
if ~brainstorm('status'); brainstorm nogui; end
import org.brainstorm.tree.BstNode;
n = BstNode('manifold', 'TEST manifold node', 'manifold_test.mat', 1, 1);
fprintf('type=%s comment=%s  (constructed OK)\n', char(n.getType()), char(n.getComment()));
```
Expected: constructs without error (icon defaults gracefully). If it throws, STOP — a Java icon branch is needed (escalate to the `bst-java` fork path) before proceeding.

- [ ] **Step 3: Commit**
```bash
git add -A
git commit -m "refactor(tess): rename tess_manifold (2-manifold validator) -> tess_repair; free name for the manifold-node assembler"
```

---

## Task 2: DB data types (`db_template` + `file_gettype`)

**Files:** `toolbox/db/db_template.m`, `toolbox/io/file_gettype.m`.

- [ ] **Step 1: Add child lists to the `surface` template**
In `db_template.m`, replace the `case 'surface'` template with one that adds the three child lists (default empty struct arrays, study-pattern):
```matlab
    case 'surface'
        template = struct('Comment',     '', ...
                          'FileName',    '', ...
                          'SurfaceType', '', ...
                          'Manifold',    struct('FileName',{},'Comment',{}), ...
                          'Eigen',       struct('FileName',{},'Comment',{},'Variant',{}), ...
                          'Operator',    struct('FileName',{},'Comment',{},'Variant',{}));
```

- [ ] **Step 2: Add the `manifoldmat` file-content template**
In `db_template.m`, add (near the other `*mat` cases, e.g. after `'surfacemat'`):
```matlab
    case 'manifoldmat'
        template = struct(...
              'Comment',       '', ...
              'ParentSurface', '', ...   % relative path to the parent surface
              'Topology',      [], ...   % 1x2 per-hemisphere nxr facet structs
              'Embedded',      [], ...
              'Intrinsic',     [], ...
              'Extrinsic',     [], ...
              'Gauge',         [], ...
              'Provenance',    []);
```

- [ ] **Step 3: Register the `manifold_` file type**
In `toolbox/io/file_gettype.m`, add a branch so a filename beginning `manifold_` returns type `'manifold'` (read the function; mirror how e.g. `headmodel`/`tess`/`cortex` prefixes are matched — typically a `~isempty(strfind(fileName,'manifold_'))` / `strncmp` test returning `'manifold'`).

- [ ] **Step 4: Verify (DB-free)**
```matlab
rehash;
s = db_template('surface');  assert(all(isfield(s,{'Manifold','Eigen','Operator'})), 'surface child lists');
m = db_template('manifoldmat'); assert(all(isfield(m,{'Topology','Embedded','Gauge','ParentSurface'})), 'manifoldmat');
assert(strcmpi(file_gettype('manifold_test_250101.mat'),'manifold'), 'file_gettype manifold');
disp('TASK2 OK');
```

- [ ] **Step 5: Commit** (`db_template.m`, `file_gettype.m`) — `feat(db): manifold/eigen/operator surface child lists + manifoldmat schema + manifold file type`.

---

## Task 3: Registration (`db_add_manifold`) + resolution (`bst_get('ManifoldFile')`)

**Files:** create `toolbox/db/db_add_manifold.m`; modify `toolbox/core/bst_get.m`.

- [ ] **Step 1: `db_add_manifold`** — mirror `db_add_surface`, but append to the parent surface's `.Manifold` list.
```matlab
function iManifold = db_add_manifold(iSubject, ParentSurfaceFile, ManifoldMat, Comment)
% DB_ADD_MANIFOLD: Save a manifold_*.mat under a surface and register it as a child node.
% Authors: Diellor Basha, 2026
    ProtocolSubjects = bst_get('ProtocolSubjects');
    if (iSubject == 0), sSubject = ProtocolSubjects.DefaultSubject; else, sSubject = ProtocolSubjects.Subject(iSubject); end
    iSurface = find(strcmpi({sSubject.Surface.FileName}, file_short(ParentSurfaceFile)), 1);
    if isempty(iSurface), error('db_add_manifold:noSurface', 'Parent surface not found: %s', ParentSurfaceFile); end
    % Save the file next to the surface (anat folder)
    anatDir = bst_fileparts(file_fullpath(ParentSurfaceFile));
    OutputFile = file_unique(bst_fullfile(anatDir, ['manifold_' bst_get('FileNameTemplate') '.mat']));  % or build a manifold_*.mat name
    ManifoldMat.ParentSurface = file_short(ParentSurfaceFile);
    if isempty(Comment), Comment = 'Manifold'; end
    ManifoldMat.Comment = Comment;
    bst_save(OutputFile, ManifoldMat, 'v7');
    % Append child entry
    sEntry = db_template('Surface'); sEntry = struct('FileName', file_short(OutputFile), 'Comment', Comment); %#ok<NASGU>
    e.FileName = file_short(OutputFile); e.Comment = Comment;
    sSubject.Surface(iSurface).Manifold(end+1) = e;
    iManifold = numel(sSubject.Surface(iSurface).Manifold);
    if (iSubject == 0), ProtocolSubjects.DefaultSubject = sSubject; else, ProtocolSubjects.Subject(iSubject) = sSubject; end
    bst_set('ProtocolSubjects', ProtocolSubjects);
    panel_protocols('UpdateNode', 'Subject', iSubject);
    db_save();
end
```
NOTE for the implementer: read `db_add_surface` for the exact `file_unique`/naming helper and the correct `bst_get('ProtocolSubjects')` round-trip; adapt the manifold filename construction to the project convention (a `manifold_<tag>.mat` in the anat dir). Keep the child-append + `bst_set` + `UpdateNode` shape above.

- [ ] **Step 2: `bst_get('ManifoldFile', FileName)`** — add a case in `bst_get.m` that scans subjects' `Surface(*).Manifold` for the filename and returns `[sSubject, iSubject, iSurface, iManifold]`. Mirror the existing `'SurfaceFile'` case (read it first).

- [ ] **Step 3: Verify on a live protocol** (create + resolve + cleanup)
```matlab
% pick a subject with a cortex surface; make a tiny manifoldmat; add; resolve; delete.
[sSubj,iSubj] = bst_get('Subject', 1);  % adjust
surfFile = sSubj.Surface(sSubj.iCortex).FileName;
mm = db_template('manifoldmat'); mm.Topology = 1; % stub content for the DB test
iM = db_add_manifold(iSubj, surfFile, mm, 'TEST manifold');
[s2,is2,isurf,im] = bst_get('ManifoldFile', s2dummy); % use the created file name
% assert the surface's .Manifold has the entry, then cleanup:
% file_delete(file_fullpath(<created>),1); reload; remove entry; db_save;
```
(Implementer: finalize the resolve+cleanup; the assertion is that `.Manifold` gains an entry and `bst_get('ManifoldFile',…)` returns it; restore the protocol afterward.)

- [ ] **Step 4: Commit** — `feat(db): db_add_manifold + bst_get ManifoldFile resolution`.

---

## Task 4: The assembler `tess_manifold`

**Files:** create `toolbox/anatomy/tess_manifold.m`; test `dev/tests/test_tess_manifold.m`.

- [ ] **Step 1: Write `tess_manifold(SurfaceFile, ...)`** — reuse the retired `tess_frame` facet logic, but write a node:
  - `in_tess_bst`; guards (nxr-compute installed; Structures L/R atlas; `Reg.Sphere` for trivial gauge).
  - `tess_hemisplit` (atlas L/R, never `conncomp`); per hemisphere build submesh, `nxr_compute('create')` → `nxr_compute('facets', h, 'trivial', struct('singVerts',poles,'singValues',[1;1]))` → `destroy`; attach `GlobalVertices/GlobalFaces/Hemisphere`.
  - Assemble `ManifoldMat` = `db_template('manifoldmat')` filled with 1×2 `Topology/Embedded/Intrinsic/Extrinsic/Gauge` + `Provenance` (Backend='nxr', NxrVersion, Gauge, ComputeDate).
  - `iSubject` via `bst_get('SurfaceFile', SurfaceFile)`; `db_add_manifold(iSubject, SurfaceFile, ManifoldMat, Comment)` unless `NoSave`.
  - Options: `Gauge` (trivial), `NoSave`, `ForceRecompute` (skip if the surface already has a manifold node and not ForceRecompute).
  Return the `ManifoldMat` struct.

- [ ] **Step 2: Test `test_tess_manifold.m`** (real 20484 cortex; backup/restore the protocol; `onCleanup`):
  - `tess_manifold(cortexFile, 'ForceRecompute', 1)` creates a `manifold_` node; the cortex's `.Manifold` list is non-empty; the saved file has 1×2 `Embedded/Topology/Gauge`.
  - `Embedded(hh).vertex.grid` orthonormal; `Gauge(hh)` Gauss-Bonnet `sum(singularity.indices)==2` (reuse facet checks).
  - Cleanup: delete the created node + restore.

- [ ] **Step 3: Run / commit** — `feat(tess-manifold): assemble nxr geometry backbone as a manifold_ node`.

---

## Task 5: Tree nesting + context menus

**Files:** `toolbox/tree/node_create_subject.m`, `toolbox/tree/tree_callbacks.m`.

- [ ] **Step 1: Nest manifold children under the surface** — in `node_create_subject.m`, inside the surface loop (after `nodeSubject.add(nodeSurface)`), iterate `sSubject.Surface(iSurface).Manifold` and add child nodes:
```matlab
        for iM = 1:numel(sSubject.Surface(iSurface).Manifold)
            chN = CreateNode('manifold', char(sSubject.Surface(iSurface).Manifold(iM).Comment), ...
                char(sSubject.Surface(iSurface).Manifold(iM).FileName), iM, iSubject, iSearch);
            if ~isempty(chN), nodeSurface.add(chN); end
        end
```
(CreateNode returns `[isCreated,node]`; adapt to capture both. Mirror the surrounding usage.)

- [ ] **Step 2: Surface "Compute manifold" menu + manifold-node menu** — in `tree_callbacks.m`:
  - In the surface (`'cortex'`/`'surface'`) context-menu section, add `gui_component('MenuItem', jPopup, [], 'Compute manifold', [], [], @(h,ev)bst_call(@tess_manifold, filenameRelative))`.
  - Add a `case 'manifold'` in the node context-menu `switch lower(nodeType)` with "View" (minimal: `@(h,ev)disp(in_bst(...))` or a stub inspector) and "Delete" (`file_delete` + remove the `.Manifold` entry + `db_save` + `UpdateNode`). Mirror an existing leaf node's delete handling.

- [ ] **Step 3: Verify in the GUI** (visual + programmatic): run `tess_manifold` on the cortex, confirm a `manifold_` node appears nested under the cortex surface in the tree, right-click shows View/Delete, and "Compute manifold" appears on the surface. Programmatic: `sSubject.Surface(iCortex).Manifold` non-empty after Compute.

- [ ] **Step 4: Commit** — `feat(tree): nest manifold nodes under surface + compute/view/delete menus`.

---

## Final verification

- [ ] End-to-end: on the real cortex, "Compute manifold" → a nested `manifold_` node appears; its file holds the nxr facet bundle; View/Delete work; `db_template`/`file_gettype`/`bst_get('ManifoldFile')` all consistent. Restore the protocol (remove test node) when done.
- [ ] Then **superpowers:finishing-a-development-branch**, and proceed to SP2 (`eigen_` node + consumer migration).

---

## Notes for the implementer
- **No Java build in SP1.** The `'manifold'` node type renders with a default icon; behavior is MATLAB-side in `tree_callbacks`. A dedicated icon is optional later work in the `bst-java` fork (branch off `dev`; never push upstream).
- **Breaking changes accepted:** the renamed `tess_repair`'s 7 callers stay broken (SP4). Do not migrate `tess_frame`/`tess_tangents`/viewers here.
- **Each core-file edit: read the surrounding code first** (`db_template`, `bst_get`, `node_create_subject`, `tree_callbacks`, `file_gettype`, `db_add_surface`) — the code blocks above are the additions; match the file's exact local idioms.
- **DB-mutating tests must restore the protocol** (delete created nodes; `db_reload_database('current')`), per the project's DB-deletion rule (`file_delete` + reload, never raw `.mat` deletion).
