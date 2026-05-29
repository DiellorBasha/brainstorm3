# Eigenmodes R2: Per-Hemisphere Computation + Registered-Result Viewer — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make eigenmode computation per-hemisphere, re-architect the viewer onto a registered Source result (fixing custom-max + giving a diverging default + uniting hemispheres), and consolidate the compute dialog with auto-open.

**Architecture:** `tess_eigenmodes` solves each connected component independently and tags every mode with `Component`/`CompRank`. The viewer builds a paired display matrix (column k = each component's rank-k mode) as a real Brainstorm Source result displayed via the standard surface path, so colormap UI works natively; ←/→ step modes via a custom key handler driving `panel_time('SetCurrentTime')`, and the result is auto-removed on figure close. A custom panel collects compute options in one dialog; compute auto-opens the viewer.

**Tech Stack:** MATLAB, Brainstorm (`tess_eigenmodes`/`tess_laplacian`/`tess_vertconn`, `conncomp`/`graph`, `in_/out_tess_eigenmodes`, `view_surface_data`, `db_add_data`, `bst_get('AnalysisIntraStudy')`, `panel_time`, `gui_show_dialog`/`gui_river`/`gui_component`, `bst_colormaps`), MATLAB MCP for tests + static analysis.

---

## Spec

`docs/superpowers/specs/2026-05-28-eigenmode-perhemisphere-viewer-r2-design.md`

## File Structure

- `toolbox/anatomy/tess_eigenmodes.m` — per-component solve; output `Component`/`CompRank`; still return whole-mesh `L`,`M`.
- `toolbox/io/out_tess_eigenmodes.m` — persist `Component`/`CompRank`.
- `toolbox/io/in_tess_eigenmodes.m` — backward-compat defaults for files without the metadata.
- `toolbox/gui/view_eigenmodes.m` — rewrite: `BuildPairedGrid` pure helper + `ViewFigure` (registered result, diverging colormap, legend, custom stepping, auto-clean).
- `toolbox/gui/panel_eigenmodes_compute.m` — **new** custom modal panel (modes + mass-type dropdown).
- `toolbox/process/functions/process_eigenmodes.m` — `ComputeInteractive`: use the panel, auto-open the viewer, drop the success msgbox.
- `dev/tests/test_eigenmodes_perhemisphere.m` — **new** (per-component computation).
- `dev/tests/test_view_eigenmodes_pure.m` — replace with `BuildPairedGrid` tests (the R1 GetModeDisplay/StepMode helpers are superseded by the registered-result viewer).
- `dev/tests/test_io_eigenmodes_roundtrip.m` — extend for the new metadata.

## Conventions

Tests follow the repo idiom: a function whose name equals the file, using `assert`, printing `ALL TESTS PASSED`. Run via MATLAB MCP `evaluate_matlab_code` calling the function name (NOT `run_matlab_test_file`). Prefix interactive runs with `rng('default')` if the session is in legacy RNG mode. GUI code (Tasks 4–5) is verified by `check_matlab_code` + interactive validation (Task 7). Commit after each task. Work on branch `feature/eigenmode-context-menus`.

---

### Task 1: Per-component computation in `tess_eigenmodes`

**Files:**
- Modify: `toolbox/anatomy/tess_eigenmodes.m`
- Test: `dev/tests/test_eigenmodes_perhemisphere.m` (new)

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigenmodes_perhemisphere.m`:

```matlab
function test_eigenmodes_perhemisphere
% Verify tess_eigenmodes solves per connected component and tags modes with
% Component/CompRank, keeps each mode localized to its component, removes one
% DC mode per component, and stays globally M-orthonormal.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% ----- Two disjoint octahedra (each a closed 2-manifold) -----
V1 = [ 1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
F1 = [5 1 3; 5 3 2; 5 2 4; 5 4 1; 6 3 1; 6 2 3; 6 4 2; 6 1 4];
V  = [V1; bsxfun(@plus, V1, [10 0 0])];   % second component offset in x
F  = [F1; F1 + 6];
nV = size(V,1);

nModes = 3;
[Eig, ~, M] = tess_eigenmodes(V, F, 'nModes', nModes, 'Verbose', 0);

% Two components, nModes each
assert(size(Eig.Vectors,2) == 2*nModes, 'Expected nModes per component (2 components).');
assert(isequal(Eig.Component(:), [1;1;1;2;2;2]), 'Component labels wrong.');
assert(isequal(Eig.CompRank(:),  [1;2;3;1;2;3]), 'CompRank wrong.');
assert(Eig.nRemoved == 2, 'Expected one DC mode removed per component.');

% Localization: component-1 modes zero on component-2 vertices and vice versa
c1 = (Eig.Component==1);  c2 = (Eig.Component==2);
assert(max(max(abs(Eig.Vectors(7:12, c1)))) < 1e-8, 'Comp-1 modes must vanish on comp-2 verts.');
assert(max(max(abs(Eig.Vectors(1:6,  c2)))) < 1e-8, 'Comp-2 modes must vanish on comp-1 verts.');

% Global M-orthonormality (block-diagonal M => per-component orthonormality)
G = Eig.Vectors' * M * Eig.Vectors;
assert(norm(G - eye(size(G)), 'fro') < 1e-6, 'Modes not M-orthonormal.');

% ----- Single-component mesh still works (one octahedron) -----
[EigS, ~, ~] = tess_eigenmodes(V1, F1, 'nModes', 2, 'Verbose', 0);
assert(all(EigS.Component == 1), 'Single component must be labeled 1.');
assert(isequal(EigS.CompRank(:), (1:numel(EigS.Values))'), 'Single-component CompRank must be 1..K.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run via MATLAB MCP `evaluate_matlab_code`: `cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); test_eigenmodes_perhemisphere`
Expected: FAIL — current `tess_eigenmodes` returns no `Component`/`CompRank` (and solves whole-mesh, so the count/localization asserts fail).

- [ ] **Step 3: Replace the solve body with a per-component loop**

In `toolbox/anatomy/tess_eigenmodes.m`, **replace everything from `%% ===== STEP 2: ASSEMBLE OPERATORS =====` (line ~148) through the end of `%% ===== STEP 7: PACKAGE OUTPUT =====` (the `end` at line ~274 that closes the main function)** with the following. (Keep STEP 1 the manifold check and `nV`/`maxModes` lines above it unchanged; keep the `mOrthonormalize` helper below unchanged.)

```matlab
%% ===== STEP 2: WHOLE-MESH OPERATORS (for return + orthonormality reference) =====
if Verbose
    fprintf('BST> tess_eigenmodes: Assembling L and M (%s mass)...\n', MassType);
end
[L, M] = tess_laplacian(Vertices, Faces, 'MassType', MassType);

%% ===== STEP 3: SPLIT INTO CONNECTED COMPONENTS =====
% Each disconnected component (e.g. a cortical hemisphere) is a separate
% manifold; its Laplace-Beltrami operator must be solved independently. Solving
% the union as one block-diagonal operator interleaves components by eigenvalue
% and, for near-degenerate eigenvalues (near-symmetric hemispheres), lets eigs
% return meaningless cross-component mixtures.
compId = conncomp(graph(tess_vertconn(Vertices, Faces)));   % [1 x nV]
nComp  = max(compId);
if Verbose
    fprintf('BST> tess_eigenmodes: %d connected component(s).\n', nComp);
end

%% ===== STEP 4: SOLVE PER COMPONENT AND ASSEMBLE =====
VectorsAll = zeros(nV, 0);
ValuesAll  = zeros(0, 1);
Component  = zeros(0, 1);
CompRank   = zeros(0, 1);
nRemoved   = 0;
for c = 1:nComp
    vIdx = find(compId == c);
    nvc  = numel(vIdx);
    % Sub-mesh: faces fully inside this component, vertex indices remapped to 1..nvc
    faceMask = all(ismember(Faces, vIdx), 2);
    remap = zeros(nV, 1);
    remap(vIdx) = 1:nvc;
    Fsub = remap(Faces(faceMask, :));
    Vsub = Vertices(vIdx, :);
    % Cap modes to this component's size
    nModesC = min(nModes, nvc - 2);
    [Uc, lambdasC, nRemC] = SolveComponent(Vsub, Fsub, nModesC, MassType, ...
        RemoveDC, Sigma, Tolerance, MaxIter, Verbose);
    % Zero-pad into full vertex space and append
    Ufull = zeros(nV, size(Uc, 2));
    Ufull(vIdx, :) = Uc;
    VectorsAll = [VectorsAll, Ufull];                              %#ok<AGROW>
    ValuesAll  = [ValuesAll;  lambdasC(:)];                        %#ok<AGROW>
    Component  = [Component;  c * ones(numel(lambdasC), 1)];       %#ok<AGROW>
    CompRank   = [CompRank;   (1:numel(lambdasC))'];               %#ok<AGROW>
    nRemoved   = nRemoved + nRemC;
end

%% ===== STEP 5: PACKAGE OUTPUT =====
totalTime = toc(tStart);
Eigenmodes = struct();
Eigenmodes.Vectors     = VectorsAll;
Eigenmodes.Values      = ValuesAll;
Eigenmodes.nModes      = numel(ValuesAll);
Eigenmodes.Component   = Component;
Eigenmodes.CompRank    = CompRank;
Eigenmodes.nComponents = nComp;
Eigenmodes.MassType    = MassType;
Eigenmodes.Sigma       = Sigma;
Eigenmodes.Tolerance   = Tolerance;
Eigenmodes.nRemoved    = nRemoved;
Eigenmodes.ComputeTime = totalTime;

if Verbose
    fprintf('BST> tess_eigenmodes: Done. %d modes across %d component(s), total %.1f s.\n', ...
        Eigenmodes.nModes, nComp, totalTime);
end

end


%% ===== HELPER: SOLVE ONE COMPONENT =====
function [U, lambdas, nRemoved] = SolveComponent(Vertices, Faces, nModes, MassType, RemoveDC, Sigma, Tolerance, MaxIter, Verbose)
    [L, M] = tess_laplacian(Vertices, Faces, 'MassType', MassType);
    nV = size(Vertices, 1);
    nExtra   = 5;
    nRequest = min(nModes + nExtra, nV - 1);
    opts = struct('tol', Tolerance, 'maxit', MaxIter, 'disp', 0);
    try
        [U, D] = eigs(L, M, nRequest, Sigma, opts);
    catch
        [U, D] = eigs(L, M, nRequest, 'smallestabs', opts);
    end
    lambdas = real(diag(D));
    [lambdas, sortIdx] = sort(lambdas);
    U = real(U(:, sortIdx));
    % Remove DC modes (one per closed component)
    nRemoved = 0;
    if RemoveDC
        dcThreshold = max(abs(Sigma) * 10, 1e-6);
        dcMask = abs(lambdas) < dcThreshold;
        nRemoved = sum(dcMask);
        lambdas = lambdas(~dcMask);
        U = U(:, ~dcMask);
    end
    % Trim to requested count
    if numel(lambdas) > nModes
        lambdas = lambdas(1:nModes);
        U = U(:, 1:nModes);
    end
    % M-orthonormalize and clamp tiny negatives
    U = mOrthonormalize(U, M);
    lambdas(lambdas < 0) = 0;
end
```

Notes for the implementer: `tStart = tic;` already exists above STEP 1 — keep it. The `nModes` reduction-vs-`maxModes` block above STEP 2 stays (it caps the whole-mesh request; per-component capping is additionally applied via `nModesC`). The function's output signature `[Eigenmodes, L, M, Vertices, Faces]` is unchanged (`L`,`M` are the whole-mesh operators from STEP 2).

- [ ] **Step 4: Run test to verify it passes**

Run: `cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); test_eigenmodes_perhemisphere`
Expected: PASS — `ALL TESTS PASSED`.

- [ ] **Step 5: Static-analysis check**

MATLAB MCP `check_matlab_code` on `toolbox/anatomy/tess_eigenmodes.m`. Expected: no errors (AGROW is suppressed; pre-existing notes acceptable).

- [ ] **Step 6: Commit**

```bash
git add toolbox/anatomy/tess_eigenmodes.m dev/tests/test_eigenmodes_perhemisphere.m
git commit -m "tess_eigenmodes: solve per connected component, tag Component/CompRank

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Persist + load `Component`/`CompRank`

**Files:**
- Modify: `toolbox/io/out_tess_eigenmodes.m`, `toolbox/io/in_tess_eigenmodes.m`
- Test: `dev/tests/test_io_eigenmodes_roundtrip.m` (extend)

- [ ] **Step 1: Add assertions to the round-trip test**

Open `dev/tests/test_io_eigenmodes_roundtrip.m`. After the existing round-trip assertions (where it has loaded the eigenmodes back into a variable — call it `EigLoaded` to match the file; if the file uses a different variable name, adapt), add:

```matlab
% ----- New per-hemisphere metadata round-trips -----
assert(isfield(EigLoaded, 'Component') && isfield(EigLoaded, 'CompRank'), ...
    'Component/CompRank must persist through save/load.');
assert(isequal(EigLoaded.Component(:), Eig.Component(:)), 'Component not preserved.');
assert(isequal(EigLoaded.CompRank(:),  Eig.CompRank(:)),  'CompRank not preserved.');

% ----- Backward compatibility: a struct WITHOUT the metadata loads with defaults -----
% Simulate an old surface file by stripping the metadata before saving.
TessTmp = load(SurfaceFileFull);                  % SurfaceFileFull = the temp surface used above
TessTmp.Eigenmodes = rmfield(TessTmp.Eigenmodes, {'Component', 'CompRank'});
bst_save(SurfaceFileFull, TessTmp, 'v7');
[EigOld, isOld] = in_tess_eigenmodes(SurfaceFileFull);
assert(isOld, 'Old-format eigenmodes must still load.');
nOld = size(EigOld.Vectors, 2);
assert(isequal(EigOld.Component(:), ones(nOld,1)), 'Missing Component must default to a single component.');
assert(isequal(EigOld.CompRank(:),  (1:nOld)'),    'Missing CompRank must default to 1..K.');
```

(If the existing test does not already build `Eig` with `Component`/`CompRank`, the simplest path is to compute `Eig` via `tess_eigenmodes` on the two-octahedra mesh from Task 1 so the metadata exists; the implementer should make the variable names consistent with the existing test header.)

- [ ] **Step 2: Run test to verify it fails**

Run: `cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); test_io_eigenmodes_roundtrip`
Expected: FAIL — `out_tess_eigenmodes` does not store the fields and `in_tess_eigenmodes` provides no defaults.

- [ ] **Step 3: Persist the metadata in `out_tess_eigenmodes`**

In `toolbox/io/out_tess_eigenmodes.m`, in the `EigenmodesStore` block, after `EigenmodesStore.nModes = Eigenmodes.nModes;` add:

```matlab
    if isfield(Eigenmodes, 'Component'),   EigenmodesStore.Component   = Eigenmodes.Component(:);   end
    if isfield(Eigenmodes, 'CompRank'),    EigenmodesStore.CompRank    = Eigenmodes.CompRank(:);    end
    if isfield(Eigenmodes, 'nComponents'), EigenmodesStore.nComponents = Eigenmodes.nComponents;    end
```

- [ ] **Step 4: Default the metadata in `in_tess_eigenmodes`**

In `toolbox/io/in_tess_eigenmodes.m`, after the block that doubles `Eigenmodes.Vectors` (just before `isComputed = true;`), add:

```matlab
% Backward compatibility: older files have no per-hemisphere metadata.
% Treat them as a single component (rank = column index).
nK = size(Eigenmodes.Vectors, 2);
if ~isfield(Eigenmodes, 'Component') || isempty(Eigenmodes.Component)
    Eigenmodes.Component = ones(nK, 1);
end
if ~isfield(Eigenmodes, 'CompRank') || isempty(Eigenmodes.CompRank)
    Eigenmodes.CompRank = (1:nK)';
end
if ~isfield(Eigenmodes, 'nComponents') || isempty(Eigenmodes.nComponents)
    Eigenmodes.nComponents = max(Eigenmodes.Component);
end
```

- [ ] **Step 5: Run test to verify it passes**

Run: `cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); test_io_eigenmodes_roundtrip`
Expected: PASS — `ALL TESTS PASSED`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/io/out_tess_eigenmodes.m toolbox/io/in_tess_eigenmodes.m dev/tests/test_io_eigenmodes_roundtrip.m
git commit -m "io eigenmodes: persist + default Component/CompRank metadata

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: `BuildPairedGrid` pure helper (viewer math)

**Files:**
- Modify: `toolbox/gui/view_eigenmodes.m` (replace the file's helpers with the dispatch + `BuildPairedGrid`)
- Test: `dev/tests/test_view_eigenmodes_pure.m` (replace contents)

- [ ] **Step 1: Replace the pure-helper test**

Overwrite `dev/tests/test_view_eigenmodes_pure.m` with:

```matlab
function test_view_eigenmodes_pure
% Verify the paired-grid builder: column k holds each component's rank-k mode,
% summed (disjoint support), so a single display column shows both hemispheres.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% Fabricate a 2-component Eigenmodes struct: comp 1 on verts 1:4, comp 2 on 5:8
nV = 8;
Eig = struct();
Eig.Vectors = zeros(nV, 4);
Eig.Vectors(1:4, 1) = [1;2;3;4];     % comp 1, rank 1
Eig.Vectors(1:4, 2) = [5;6;7;8];     % comp 1, rank 2
Eig.Vectors(5:8, 3) = [9;10;11;12];  % comp 2, rank 1
Eig.Vectors(5:8, 4) = [13;14;15;16]; % comp 2, rank 2
Eig.Values    = [1;2;1;2];
Eig.Component = [1;1;2;2];
Eig.CompRank  = [1;2;1;2];

[Grid, K, Info] = view_eigenmodes('BuildPairedGrid', Eig);
assert(K == 2, 'K must equal the max within-component rank.');
assert(isequal(size(Grid), [nV, 2]), 'Grid must be [nV x K].');
% Column 1 = comp1-rank1 (verts 1:4) + comp2-rank1 (verts 5:8)
assert(isequal(Grid(:,1), [1;2;3;4;9;10;11;12]), 'Paired column 1 wrong.');
% Column 2 = comp1-rank2 + comp2-rank2
assert(isequal(Grid(:,2), [5;6;7;8;13;14;15;16]), 'Paired column 2 wrong.');
assert(Info.K == 2 && numel(Info.Values) == 4, 'Info fields wrong.');

% Single-component fallback (no metadata): each column maps to its own rank
Eig2 = struct('Vectors', magic(4), 'Values', (1:4)');
[Grid2, K2] = view_eigenmodes('BuildPairedGrid', Eig2);
assert(K2 == 4 && isequal(Grid2, magic(4)), 'Single-component grid must pass modes through.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); test_view_eigenmodes_pure`
Expected: FAIL — `view_eigenmodes('BuildPairedGrid', ...)` not implemented (current file dispatches `GetModeDisplay`/`StepMode`).

- [ ] **Step 3: Replace `view_eigenmodes.m` header + helpers**

Overwrite `toolbox/gui/view_eigenmodes.m` with the dispatch + `BuildPairedGrid` below. (`ViewFigure` is added in Task 4; the test only exercises `BuildPairedGrid`.)

```matlab
function varargout = view_eigenmodes(varargin)
% VIEW_EIGENMODES: Browse Laplace-Beltrami eigenmodes on a surface.
%
% USAGE:  hFig = view_eigenmodes(SurfaceFile)
%         [Grid, K, Info] = view_eigenmodes('BuildPairedGrid', Eigenmodes)
%
% The viewer displays each component's rank-k mode together (mode k shows both
% hemispheres) as a registered Brainstorm Source result, so the standard colormap
% UI applies; Left/Right arrows step modes.
%
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

if (nargin >= 1) && ischar(varargin{1}) && strcmp(varargin{1}, 'BuildPairedGrid')
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: paired display grid (column k = each component's CompRank==k mode) =====
function [Grid, K, Info] = BuildPairedGrid(Eig)
    nV = size(Eig.Vectors, 1);
    nK = size(Eig.Vectors, 2);
    if isfield(Eig, 'CompRank') && ~isempty(Eig.CompRank)
        CompRank  = Eig.CompRank(:);
        Component = Eig.Component(:);
    else
        CompRank  = (1:nK)';
        Component = ones(nK, 1);
    end
    K = max(CompRank);
    Grid = zeros(nV, K);
    for k = 1:K
        cols = find(CompRank == k);
        if ~isempty(cols)
            Grid(:, k) = sum(Eig.Vectors(:, cols), 2);   % disjoint support across components
        end
    end
    Info = struct('K', K, 'Component', Component, 'CompRank', CompRank, 'Values', Eig.Values(:));
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); test_view_eigenmodes_pure`
Expected: PASS — `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/view_eigenmodes.m dev/tests/test_view_eigenmodes_pure.m
git commit -m "view_eigenmodes: add BuildPairedGrid helper (mode k = each component's rank-k)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Viewer `ViewFigure` — registered result, diverging colormap, stepping, auto-clean

**Files:**
- Modify: `toolbox/gui/view_eigenmodes.m` (append `ViewFigure` + nested callbacks)

GUI code — verified by `check_matlab_code` here and interactively in Task 7. The auto-clean-on-close lifecycle is the known risk (spec §B): if it proves unreliable interactively, fall back to leaving the named node (delete the `DeleteFcn` wiring) and report it.

- [ ] **Step 1: Append `ViewFigure`**

Append to `toolbox/gui/view_eigenmodes.m`:

```matlab
%% ===== GUI: display the paired modes as a transient registered Source result =====
function hFig = ViewFigure(SurfaceFile, ~)
    hFig = [];
    % Load eigenmodes
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed || isempty(Eig) || ~isfield(Eig, 'Vectors') || isempty(Eig.Vectors) || ~isfield(Eig, 'Values')
        bst_error(['No eigenmodes found on this surface.' 10 ...
                   'Right-click the cortex and run "Compute eigenmodes" first.'], 'View eigenmodes', 0);
        return;
    end
    % Resolve subject + intra study
    [sSubject, iSubject] = bst_get('SurfaceFile', SurfaceFile); %#ok<ASGLU>
    [sStudy, iStudy] = bst_get('AnalysisIntraStudy', iSubject);
    if isempty(iStudy)
        bst_error('Could not find the intra-subject study.', 'View eigenmodes', 0);
        return;
    end
    % Build paired display grid (mode k shows every component's rank-k mode)
    [Grid, K, Info] = BuildPairedGrid(Eig);

    % Build a Source result (modes as the "time" axis)
    ResMat = db_template('resultsmat');
    ResMat.ImageGridAmp  = Grid;
    ResMat.ImagingKernel = [];
    ResMat.nComponents   = 1;
    ResMat.Time          = 1:K;
    ResMat.SurfaceFile   = SurfaceFile;
    ResMat.HeadModelType = 'surface';
    ResMat.nAvg          = 1;
    ResMat.Leff          = 1;
    ResMat.ColormapType  = 'stat2';   % diverging, non-absolute (signed +/- lobes) without touching 'source'
    ResMat.Comment       = sprintf('Eigenmode viewer (%d modes/component, %d component(s))', K, Eig.nComponents);
    ResMat = bst_history('add', ResMat, 'eigenmodes_view', 'Transient eigenmode viewer result');

    % Save to the intra study and register
    StudyDir   = bst_fileparts(file_fullpath(sStudy.FileName));
    OutputFile = bst_process('GetNewFilename', StudyDir, 'results_eigenview');
    bst_save(OutputFile, ResMat, 'v6');
    db_add_data(iStudy, OutputFile, ResMat);

    % Display via the standard surface-data path (colormap UI works natively)
    hFig = view_surface_data(SurfaceFile, file_short(OutputFile));
    if isempty(hFig)
        bst_error('Could not open the surface figure.', 'View eigenmodes', 0);
        return;
    end
    set(hFig, 'Name', ['Eigenmodes: ' SurfaceFile]);

    % Current mode (closure state), starting at mode 1
    curMode = 1;
    % Bottom-left legend
    hLabel = uicontrol('Style', 'text', 'String', '...', 'Units', 'Pixels', ...
        'Position', [6 0 560 20], 'HorizontalAlignment', 'left', ...
        'FontUnits', 'points', 'FontSize', bst_get('FigFont'), ...
        'ForegroundColor', [.9 .9 .9], 'BackgroundColor', [0 0 0], 'Parent', hFig);
    % Custom keyboard stepping (drives the global time = mode index) + legend
    KeyPressFcn_bak = get(hFig, 'KeyPressFcn');
    set(hFig, 'KeyPressFcn', @KeyPress_Callback);
    % Auto-remove the transient result when the figure is destroyed
    set(hFig, 'DeleteFcn', @(h,e) CleanupResult());
    % Initial position + legend
    SetMode(1);

    % ===== NESTED: move to mode k and refresh the legend =====
    function SetMode(k)
        curMode = min(max(round(k), 1), K);
        panel_time('SetCurrentTime', curMode);   % Time vector is 1:K, so time==mode index
        % Per-component eigenvalue(s) for this mode (rank == curMode)
        lv = Info.Values(Info.CompRank == curMode);
        if numel(lv) >= 2
            lamStr = sprintf('lambda = [%.4g, %.4g]', lv(1), lv(2));
        elseif ~isempty(lv)
            lamStr = sprintf('lambda = %.4g', lv(1));
        else
            lamStr = 'lambda = n/a';
        end
        set(hLabel, 'String', sprintf('Mode %d / %d     %s', curMode, K, lamStr));
    end

    % ===== NESTED: keyboard navigation =====
    function KeyPress_Callback(h, keyEvent)
        switch (keyEvent.Key)
            case 'leftarrow',  SetMode(curMode - 1);
            case 'rightarrow', SetMode(curMode + 1);
            case 'pageup',     SetMode(curMode + 10);
            case 'pagedown',   SetMode(curMode - 10);
            case 'h'
                java_dialog('msgbox', ...
                    ['Eigenmode viewer shortcuts:' 10 10 ...
                     '   Left / Right arrow  :  previous / next mode' 10 ...
                     '   Page Up / Page Down :  +/- 10 modes' 10 ...
                     '   H                   :  this help'], 'Eigenmode viewer');
            otherwise
                if ~isempty(KeyPressFcn_bak)
                    KeyPressFcn_bak(h, keyEvent);
                end
        end
    end

    % ===== NESTED: delete the transient result on figure close =====
    function CleanupResult()
        try
            file_delete(file_fullpath(OutputFile), 1);
            db_reload_studies(iStudy);
        catch
            % Non-fatal: leave the node if cleanup fails (user can delete it)
        end
    end
end
```

- [ ] **Step 2: Static-analysis check**

MATLAB MCP `check_matlab_code` on `toolbox/gui/view_eigenmodes.m`. Expected: no errors (nested shared vars `curMode`/`hLabel`/`KeyPressFcn_bak`/`OutputFile`/`iStudy` and the `%#ok<ASGLU>` are expected).

- [ ] **Step 3: No-regression on the pure helper**

Run: `cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); test_view_eigenmodes_pure`
Expected: PASS (appending `ViewFigure` must not break `BuildPairedGrid`).

- [ ] **Step 4: Commit**

```bash
git add toolbox/gui/view_eigenmodes.m
git commit -m "view_eigenmodes: registered-result viewer (paired modes, diverging colormap, stepping, auto-clean)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Custom compute panel + `ComputeInteractive` (consolidated dialog + auto-open)

**Files:**
- Create: `toolbox/gui/panel_eigenmodes_compute.m`
- Modify: `toolbox/process/functions/process_eigenmodes.m` (`ComputeInteractive`)

GUI — `check_matlab_code` + the existing options test (no regression) + interactive (Task 7).

- [ ] **Step 1: Create the panel**

Create `toolbox/gui/panel_eigenmodes_compute.m`:

```matlab
function varargout = panel_eigenmodes_compute(varargin)
% PANEL_EIGENMODES_COMPUTE: Options for computing Laplace-Beltrami eigenmodes (GUI).
%
% USAGE:  bstPanel = panel_eigenmodes_compute('CreatePanel')
%                s = panel_eigenmodes_compute('GetPanelContents')

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

eval(macro_method);
end


%% ===== CREATE PANEL =====
function [bstPanelNew, panelName] = CreatePanel() %#ok<DEFNU>
    panelName = 'EigenmodesCompute';
    import java.awt.*;
    import javax.swing.*;
    MassList = {'barycentric', 'voronoi', 'galerkin'};

    jPanelNew = gui_river([5,5], [10,15,12,10], 'Compute eigenmodes');
        gui_component('label', jPanelNew, '', 'Number of eigenmodes per hemisphere: ', [], [], [], []);
        jTextN = gui_component('text', jPanelNew, 'tab', '300', [], [], [], []);
        gui_component('label', jPanelNew, 'br', 'Mass matrix type: ', [], [], [], []);
        jComboMass = gui_component('combobox', jPanelNew, 'tab', [], {MassList}, [], [], []);
        jComboMass.setSelectedIndex(0);   % barycentric
    % Validation buttons
    gui_component('button', jPanelNew, 'br right', 'Cancel', [], [], @ButtonCancel_Callback, []);
    gui_component('button', jPanelNew, [], 'OK', [], [], @ButtonOk_Callback, []);

    bst_mutex('create', panelName);
    ctrl = struct('jTextN', jTextN, 'jComboMass', jComboMass);
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);

    function ButtonCancel_Callback(varargin)
        gui_hide(panelName);
    end
    function ButtonOk_Callback(varargin)
        bst_mutex('release', panelName);
    end
end


%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents() %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'EigenmodesCompute');
    nModes = str2double(char(ctrl.jTextN.getText()));
    if isnan(nModes) || (nModes < 1)
        error('Number of eigenmodes must be a positive integer.');
    end
    s.nModes   = round(nModes);
    s.MassType = char(ctrl.jComboMass.getSelectedItem());
end
```

- [ ] **Step 2: Static-analysis check**

MATLAB MCP `check_matlab_code` on `toolbox/gui/panel_eigenmodes_compute.m`. Expected: no errors.

- [ ] **Step 3: Replace `ComputeInteractive` to use the panel + auto-open viewer**

In `toolbox/process/functions/process_eigenmodes.m`, replace the entire `ComputeInteractive` subfunction body with:

```matlab
%% ===== COMPUTE INTERACTIVE (called from the cortex right-click menu) =====
function ComputeInteractive(iSubject, SurfaceFile) %#ok<DEFNU>
    % One consolidated dialog: number of modes + mass type
    opts = gui_show_dialog('Compute eigenmodes', @panel_eigenmodes_compute);
    if isempty(opts)
        return;   % cancelled
    end
    nModes   = opts.nModes;
    MassType = opts.MassType;
    % Confirm overwrite if eigenmodes already exist
    [~, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if isComputed
        isOverwrite = java_dialog('confirm', ...
            'Eigenmodes already exist for this surface. Overwrite?', 'Compute eigenmodes');
        if ~isOverwrite
            return;
        end
    end
    % Compute (RemoveDC=true, Repair=false, isInteractive=true)
    errMsg = Compute(SurfaceFile, nModes, MassType, true, false, true);
    if ~isempty(errMsg)
        bst_error(errMsg, 'Compute eigenmodes', 0);
        return;
    end
    % Reflect the new eigenmodes in the database
    db_reload_subjects(iSubject);
    % Visual confirmation: open the viewer on the freshly-computed modes
    view_eigenmodes(SurfaceFile);
end
```

- [ ] **Step 4: Static-analysis check**

MATLAB MCP `check_matlab_code` on `toolbox/process/functions/process_eigenmodes.m`. Expected: no errors.

- [ ] **Step 5: No-regression on the process options test**

Run: `cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); test_process_eigenmodes_options`
Expected: PASS (process description/options unchanged).

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/panel_eigenmodes_compute.m toolbox/process/functions/process_eigenmodes.m
git commit -m "Compute eigenmodes: single custom-panel dialog + auto-open viewer

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: Downstream regression

**Files:** none (verification only)

- [ ] **Step 1: Run the eigenmode test suite**

Run each via MATLAB MCP `evaluate_matlab_code` (prefix `cd('/Users/diellorbasha/workspace/research/code/brainstorm3');`), expecting `ALL TESTS PASSED` for each:

```
test_eigenmodes_perhemisphere
test_view_eigenmodes_pure
test_io_eigenmodes_roundtrip
test_eigenmodes_manifold_gate
test_process_eigenmodes_options
test_eigenmodes_transform_pure
test_eigenmodes_project_pure
test_eigenmodes_filter_pure
test_eigenmodes_filter_gain_pure
test_eigenmodes_noisefloor_pure
test_eigenmodes_dispersion_pure
test_eigenmodes_wavelet_pure
test_eigenmodes_wiener_pure
test_process_eigenmodes_transform_options
test_process_eigenmodes_denoise_options
test_process_eigenmodes_dispersion_options
test_process_eigenmodes_wavelet_options
test_process_eigenmodes_wiener_options
test_process_eigenmodes_coeffsfilter_options
```

Expected: all PASS. If any fails, STOP and report — the per-component storage change likely affected that consumer; fix before proceeding.

- [ ] **Step 2: Commit (only if a fix was needed)**

If a downstream fix was required, commit it with a message describing the consumer and the adjustment. If everything passed, no commit (this task is a gate).

---

### Task 7: Interactive validation (user) + report

GUI behavior can't be unit-tested; this is the user validation checkpoint.

- [ ] **Step 1: Launch / refresh Brainstorm** from this checkout (`rehash; clear functions` if menus/code don't refresh). Use a subject with a manifold ico cortex (per the project notes).

- [ ] **Step 2: Compute** — right-click cortex → *Compute eigenmodes*. Confirm the **single dialog** with a numeric field + **mass-type dropdown**, OK/Cancel. On OK, after the progress bar the **viewer opens automatically** (no success popup). Re-run → overwrite prompt appears.

- [ ] **Step 3: United stepping** — in the viewer, ←/→ step modes; confirm **both hemispheres update together** (mode k on left and right), PgUp/PgDn jump ±10, H help, legend shows `Mode k / K` + eigenvalues.

- [ ] **Step 4: Colormap** — confirm **diverging default** (signed +/− lobes, not magnitude); right-click colorbar → *Colormap: stat2* → **Maximum → Custom** opens the dialog and applies (#4 fixed).

- [ ] **Step 5: Auto-clean** — close the viewer figure; confirm the transient `results_eigenview*` node is removed from the subject's (Analysis) folder (no clutter). If it lingers or errors, note it — fallback is to keep the named node.

- [ ] **Step 6: Report** to the user: dialog/auto-open, united stepping, custom-max, diverging default, auto-clean outcome, and any anomalies (especially the auto-clean lifecycle).

---

## Self-Review

**Spec coverage:**
- #1 consolidated dialog → Task 5 (panel). ✓
- #2 auto-open viewer, no msgbox → Task 5 (`ComputeInteractive`). ✓
- #3 per-hemisphere computation → Task 1; metadata persistence → Task 2; united viz (paired grid) → Tasks 3–4. ✓
- #4 custom-max works → Task 4 (registered result). Validated Task 7 Step 4. ✓
- #5 diverging default → Task 4 (`ColormapType='stat2'`). ✓
- Backward compatibility → Task 2 (defaults). Downstream regression → Task 6. ✓

**Placeholder scan:** every code step has complete code; every run step states the exact MCP call + expected result. The auto-clean risk is a real implementation with a documented fallback, not a placeholder.

**Type/name consistency:** `view_eigenmodes('BuildPairedGrid', Eig)` → `[Grid, K, Info]` used identically in Task 3 test, Task 3 impl, and Task 4 `ViewFigure`. `Eigenmodes` fields `Component`/`CompRank`/`nComponents` are produced in Task 1, persisted/defaulted in Task 2, and consumed in Tasks 3–4. `panel_eigenmodes_compute('GetPanelContents')` returns `.nModes`/`.MassType`, consumed by `ComputeInteractive` in Task 5. `Compute(SurfaceFile, nModes, MassType, true, false, true)` matches the existing core signature. `ResMat.Time = 1:K` aligns with `panel_time('SetCurrentTime', curMode)` (time == mode index).

**Known residual risks (validated interactively, per spec):** the auto-clean `DeleteFcn` lifecycle; `panel_time('SetCurrentTime')` couples to global time if other time figures are open (standard Brainstorm behavior); `'stat2'` is a shared colormap type (adjusting it affects stat displays) — acceptable per spec, dedicated type deferred.
