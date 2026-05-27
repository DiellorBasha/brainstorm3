# Eigenmode Function Suite Reorganization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Relocate the five misplaced `tess_eigenmodes_*` functions out of `toolbox/anatomy/` into their correct Brainstorm folders (io, math, inverse), switch them to Brainstorm io/database conventions, and make `tess_eigenmodes`/`tess_laplacian` validate manifoldness without ever auto-repairing.

**Architecture:** Eigenmodes stay embedded as a field in the cortex surface `.mat` file; io goes through `in_tess_bst`/`bst_save`. The compute functions check manifoldness and error on non-manifold input unless repair is explicitly opted in (repair changes vertex counts and would desync head models/lead fields). The spectral math functions (`project`, `filter`) become pure (data + mass matrix in, data out); the process layer does the file io and computes the mass matrix once.

**Tech Stack:** MATLAB (R2023b), Brainstorm toolbox. Tests are MATLAB scripts in `dev/tests/` run via the MATLAB MCP (`run_matlab_file`). Each test bootstraps `brainstorm nogui` for path + `GlobalData`.

**Spec:** `docs/superpowers/specs/2026-05-27-eigenmode-functions-reorganization-design.md`

---

## Conventions used throughout

- **Running a test:** invoke the MATLAB MCP tool `run_matlab_file` with the absolute path to the test script. A passing test prints a final `ALL TESTS PASSED: <name>` line and raises no error. A failing test raises a MATLAB error (assertion or undefined function).
- **M-Lint:** invoke the MATLAB MCP tool `check_matlab_code` on the file. Expected: no new errors (pre-existing Brainstorm idiom warnings such as unused `varargout` are acceptable).
- **`tess_sphere(N)`** returns `[Vertices, Faces]` for a closed icosphere; valid `N` values include 42, 162, 642, 2562 (`toolbox/anatomy/tess_sphere.m:37`).
- **Test bootstrap block** (top of every test script):
  ```matlab
  thisDir  = fileparts(mfilename('fullpath'));
  repoRoot = fileparts(fileparts(thisDir));   % dev/tests -> repo root
  addpath(repoRoot);
  if ~brainstorm('status')
      brainstorm nogui
  end
  ```

---

## Task 1: Manifold gate — `tess_eigenmodes` + `process_eigenmodes`

Make the compute path validate manifoldness and never auto-repair. Repair is opt-in via a new `Repair` option (library) and an unchecked-by-default checkbox + interactive dialog (process).

**Files:**
- Modify: `toolbox/anatomy/tess_eigenmodes.m` (option parse + Step 1 + docstring)
- Modify: `toolbox/process/functions/process_eigenmodes.m` (GetDescription checkbox, Run, Compute, ComputeInteractive)
- Test: `dev/tests/test_eigenmodes_manifold_gate.m` (new)
- Test: `dev/tests/test_process_eigenmodes_options.m` (new)

- [ ] **Step 1: Write the failing library test**

Create `dev/tests/test_eigenmodes_manifold_gate.m`:
```matlab
function test_eigenmodes_manifold_gate
% Verify tess_eigenmodes validates manifoldness and never auto-repairs.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% ---- Clean manifold mesh ----
[V, F] = tess_sphere(642);
nV0 = size(V, 1);

% Case 1: clean mesh computes; vertex count preserved.
[Em, ~, ~, Vout] = tess_eigenmodes(V, F, 'nModes', 20, 'Verbose', 0);
assert(Em.nModes == 20, 'Expected 20 modes, got %d', Em.nModes);
assert(size(Vout,1) == nV0, 'Clean mesh vertex count changed (%d -> %d)', nV0, size(Vout,1));
fprintf('PASSED: clean manifold mesh computes, vertices preserved.\n');

% ---- Non-manifold mesh: extra face creates an edge shared by 3 faces ----
vExtra = find(~ismember(1:nV0, F(1,:)), 1);
Fnm = [F; F(1,1), F(1,2), vExtra];

% Case 2: non-manifold errors by DEFAULT.
threw = false;
try
    tess_eigenmodes(V, Fnm, 'nModes', 20, 'Verbose', 0);
catch ME
    threw = true;
    assert(~isempty(regexpi(ME.message, 'manifold')), ...
        'Expected a manifold error, got: %s', ME.message);
end
assert(threw, 'Non-manifold input did NOT error by default.');
fprintf('PASSED: non-manifold errors by default.\n');

% Case 3: explicit Repair=true succeeds.
Em3 = tess_eigenmodes(V, Fnm, 'nModes', 10, 'Repair', true, 'Verbose', 0);
assert(Em3.nModes >= 1, 'Repair path produced no modes.');
fprintf('PASSED: Repair=true repairs and computes.\n');

fprintf('ALL TESTS PASSED: test_eigenmodes_manifold_gate\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run via MATLAB MCP `run_matlab_file`: `dev/tests/test_eigenmodes_manifold_gate.m`
Expected: FAIL at Case 2 — current `FixMesh` default auto-repairs, so no error is thrown (`assert(threw, ...)` fails).

- [ ] **Step 3: Rename the `FixMesh` option to `Repair` (default false) in `tess_eigenmodes.m`**

In the PARSE INPUTS block (`tess_eigenmodes.m:88-108`), change:
```matlab
FixMesh   = true;
```
to:
```matlab
Repair    = false;
```
and change the switch case:
```matlab
        case 'fixmesh',   FixMesh   = varargin{i+1};
```
to:
```matlab
        case 'repair',    Repair    = varargin{i+1};
```

- [ ] **Step 4: Replace the Step-1 repair block with a check-only gate**

Replace `tess_eigenmodes.m:112-116`:
```matlab
%% ===== STEP 1: VALIDATE AND REPAIR MESH =====
if FixMesh
    % tess_manifold validates first and repairs only if defective (no-op on a clean mesh).
    [Vertices, Faces, isManifold] = tess_manifold(Vertices, Faces, 'Repair', 1, 'Verbose', Verbose); %#ok<ASGLU>
end
```
with:
```matlab
%% ===== STEP 1: VALIDATE MESH (NEVER AUTO-REPAIR) =====
% Repair changes vertex/edge counts, which desyncs head models, lead fields,
% and atlases built on this surface. So validate only; repair is opt-in.
[~, ~, isManifold, report] = tess_manifold(Vertices, Faces, 'Repair', 0, 'Verbose', Verbose);
if ~isManifold
    if Repair
        if Verbose
            fprintf('BST> tess_eigenmodes: Surface non-manifold; attempting repair (vertex count may change)...\n');
        end
        [Vertices, Faces] = tess_manifold(Vertices, Faces, 'Repair', 1, 'Verbose', Verbose);
    else
        if isfield(report, 'summary') && ~isempty(report.summary)
            defects = strjoin(report.summary, '; ');
        else
            defects = 'non-manifold mesh';
        end
        error(['tess_eigenmodes: Surface is not a clean 2-manifold (%s). Re-mesh with ' ...
               'icosphere downsampling, or pass ''Repair'',true to attempt a risky repair ' ...
               'that changes the vertex count.'], defects);
    end
end
```

- [ ] **Step 5: Update the `tess_eigenmodes` docstring**

In the OPTIONS section of the header comment (`tess_eigenmodes.m:34`), replace the `FixMesh` line:
```matlab
%     FixMesh    : Repair non-manifold defects (tess_manifold) before assembly (default: true)
```
with:
```matlab
%     Repair     : Attempt to repair non-manifold defects before assembly
%                  (default: false). When false, a non-manifold surface raises
%                  an error. Repair changes vertex/edge counts (risky).
```

- [ ] **Step 6: Run the library test to verify it passes**

Run via MATLAB MCP `run_matlab_file`: `dev/tests/test_eigenmodes_manifold_gate.m`
Expected: PASS — prints `ALL TESTS PASSED: test_eigenmodes_manifold_gate`.

- [ ] **Step 7: Write the failing process-options test**

Create `dev/tests/test_process_eigenmodes_options.m`:
```matlab
function test_process_eigenmodes_options
% Verify the Compute-eigenmodes process exposes an opt-in, unchecked repair box.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes('GetDescription');
assert(isfield(sProcess.options, 'repair'), 'Missing "repair" option.');
assert(~isfield(sProcess.options, 'fixmesh'), 'Old "fixmesh" option still present.');
assert(sProcess.options.repair.Value == 0, 'Repair must default to 0 (no auto-repair).');
assert(~isempty(regexpi(sProcess.options.repair.Comment, 'risky')), ...
    'Repair comment must warn that repair is risky.');
fprintf('ALL TESTS PASSED: test_process_eigenmodes_options\n');
end
```

- [ ] **Step 8: Run the process-options test to verify it fails**

Run via MATLAB MCP `run_matlab_file`: `dev/tests/test_process_eigenmodes_options.m`
Expected: FAIL — `sProcess.options.repair` does not exist yet (`isfield(...,'repair')` is false).

- [ ] **Step 9: Flip the checkbox in `process_eigenmodes.m` GetDescription**

Replace `process_eigenmodes.m:96-99`:
```matlab
    % === FIX MESH ===
    sProcess.options.fixmesh.Comment = 'Repair non-manifold defects before computing';
    sProcess.options.fixmesh.Type    = 'checkbox';
    sProcess.options.fixmesh.Value   = 1;
```
with:
```matlab
    % === REPAIR (RISKY, OPT-IN) ===
    sProcess.options.repair.Comment = 'Attempt repair if surface is non-manifold (risky: changes vertex count)';
    sProcess.options.repair.Type    = 'checkbox';
    sProcess.options.repair.Value   = 0;
```

- [ ] **Step 10: Update the `Run` option read and `Compute` call**

In `process_eigenmodes.m`, replace `:129`:
```matlab
    FixMesh   = sProcess.options.fixmesh.Value;
```
with:
```matlab
    Repair    = sProcess.options.repair.Value;
```
and replace the `Compute` invocation at `:171`:
```matlab
    errMsg = Compute(SurfaceFile, nModes, MassType, RemoveDC, FixMesh);
```
with:
```matlab
    errMsg = Compute(SurfaceFile, nModes, MassType, RemoveDC, Repair, false);
```

- [ ] **Step 11: Update the `Compute` signature, add the manifold pre-check, pass `Repair`**

Replace the `Compute` signature at `process_eigenmodes.m:186`:
```matlab
function errMsg = Compute(SurfaceFile, nModes, MassType, RemoveDC, FixMesh)
```
with:
```matlab
function errMsg = Compute(SurfaceFile, nModes, MassType, RemoveDC, Repair, isInteractive)
```
Immediately after the `Vertices`/`Faces` are loaded (`process_eigenmodes.m:201-202`, the lines `Vertices = double(TessMat.Vertices);` and `Faces = double(TessMat.Faces);`), insert:
```matlab

    % ===== MANIFOLD CHECK (never auto-repair) =====
    [~, ~, isManifold, report] = tess_manifold(Vertices, Faces, 'Repair', 0, 'Verbose', 0);
    if ~isManifold && ~Repair
        if isfield(report, 'summary') && ~isempty(report.summary)
            defects = strjoin(report.summary, '; ');
        else
            defects = 'non-manifold mesh';
        end
        if isInteractive
            resp = java_dialog('question', ...
                ['Surface is not a clean 2-manifold (' defects ').' 10 10 ...
                 'Attempt repair? This changes the vertex count and may desync the ' ...
                 'head model, lead field, and atlas built on this surface.'], ...
                'Compute eigenmodes', [], {'Repair', 'Cancel'}, 'Cancel');
            if ~strcmpi(resp, 'Repair')
                errMsg = 'Cancelled: surface is non-manifold and repair was declined.';
                return;
            end
            Repair = true;
        else
            errMsg = ['Surface is not a clean 2-manifold (' defects '). ' ...
                      'Re-mesh with icosphere downsampling, or enable the repair option.'];
            return;
        end
    end
```
Then in the `tess_eigenmodes` call (`process_eigenmodes.m:210-215`), replace the option line:
```matlab
            'FixMesh',   FixMesh, ...
```
with:
```matlab
            'Repair',    Repair, ...
```

- [ ] **Step 12: Update the `ComputeInteractive` call to `Compute`**

Replace `process_eigenmodes.m:278`:
```matlab
    errMsg = Compute(SurfaceFile, nModes, MassType, true, true);
```
with:
```matlab
    errMsg = Compute(SurfaceFile, nModes, MassType, true, false, true);
```
(RemoveDC=true, Repair=false so the dialog decides, isInteractive=true.)

- [ ] **Step 13: Run the process-options test to verify it passes**

Run via MATLAB MCP `run_matlab_file`: `dev/tests/test_process_eigenmodes_options.m`
Expected: PASS — prints `ALL TESTS PASSED: test_process_eigenmodes_options`.

- [ ] **Step 14: M-Lint the two modified source files**

Run MATLAB MCP `check_matlab_code` on `toolbox/anatomy/tess_eigenmodes.m` and `toolbox/process/functions/process_eigenmodes.m`.
Expected: no new errors.

- [ ] **Step 15: Commit**
```bash
git add toolbox/anatomy/tess_eigenmodes.m toolbox/process/functions/process_eigenmodes.m \
        dev/tests/test_eigenmodes_manifold_gate.m dev/tests/test_process_eigenmodes_options.m
git commit -m "Eigenmodes: validate manifoldness, never auto-repair (opt-in only)"
```

---

## Task 2: `tess_laplacian` — assume-manifold doc + optional `CheckManifold`

Document that the operator assumes a clean 2-manifold, and add an off-by-default `CheckManifold` flag that warns on non-manifold input.

**Files:**
- Modify: `toolbox/anatomy/tess_laplacian.m` (option parse + check + docstring)
- Test: `dev/tests/test_laplacian_ico.m` (new)

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_laplacian_ico.m`:
```matlab
function test_laplacian_ico
% Verify cotangent Laplacian on an icosphere, and the CheckManifold warning.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

[V, F] = tess_sphere(642);
[L, M] = tess_laplacian(V, F);   % default barycentric mass

% L symmetric
assert(norm(L - L', 1) < 1e-9, 'L is not symmetric.');
% Barycentric M is diagonal and strictly positive
md = full(diag(M));
assert(all(md > 0), 'M has non-positive diagonal entries.');
assert(nnz(M - spdiags(md, 0, numel(md), numel(md))) == 0, 'Barycentric M is not diagonal.');
% L positive semidefinite
lmin = eigs(L, 1, 'smallestreal');
assert(lmin > -1e-6, 'L is not PSD (min eig %.3e).', lmin);
% Constant function is in the null space: row sums ~ 0
assert(max(abs(sum(L, 2))) < 1e-6, 'L row sums are not ~0.');
fprintf('PASSED: L symmetric/PSD, M positive diagonal, constant null space.\n');

% CheckManifold warns on non-manifold input
vExtra = find(~ismember(1:size(V,1), F(1,:)), 1);
Fnm = [F; F(1,1), F(1,2), vExtra];
lastwarn('');
tess_laplacian(V, Fnm, 'CheckManifold', true);
[w, ~] = lastwarn();
assert(~isempty(w), 'CheckManifold=true did not warn on non-manifold input.');
fprintf('PASSED: CheckManifold warns on non-manifold input.\n');

fprintf('ALL TESTS PASSED: test_laplacian_ico\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run via MATLAB MCP `run_matlab_file`: `dev/tests/test_laplacian_ico.m`
Expected: FAIL at the CheckManifold block — the flag is ignored today, so no warning is issued (`assert(~isempty(w), ...)` fails). (The numerical asserts pass; they characterize existing behavior.)

- [ ] **Step 3: Add the `CheckManifold` option to `tess_laplacian.m`**

Replace the PARSE INPUTS block (`tess_laplacian.m:80-87`):
```matlab
MassType = 'barycentric';
Symmetrize = true;
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'masstype',    MassType = lower(varargin{i+1});
        case 'symmetrize',  Symmetrize = varargin{i+1};
    end
end
```
with:
```matlab
MassType = 'barycentric';
Symmetrize = true;
CheckManifold = false;
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'masstype',      MassType = lower(varargin{i+1});
        case 'symmetrize',    Symmetrize = varargin{i+1};
        case 'checkmanifold', CheckManifold = varargin{i+1};
    end
end
```

- [ ] **Step 4: Add the optional manifold check after input validation**

After the input-validation block (immediately after `Faces = double(Faces);` at `tess_laplacian.m:94`), insert:
```matlab

% Optional manifold check (OFF by default). This operator is recomputed on
% every project/filter call, and tess_eigenmodes is the authoritative gate,
% so a full manifold scan in the inner loop would be wasteful.
if CheckManifold
    [~, ~, isManifold] = tess_manifold(Vertices, Faces, 'Repair', 0, 'Verbose', 0);
    if ~isManifold
        warning('tess_laplacian:NonManifold', ...
            ['Input mesh is not a clean 2-manifold; the cotangent Laplacian may be ' ...
             'inaccurate. Validate with tess_manifold or re-mesh with icosphere downsampling.']);
    end
end
```

- [ ] **Step 5: Update the docstring**

In the header DESCRIPTION (`tess_laplacian.m:7`, after the first description line), add a sentence and an option entry. Under OPTIONS (`tess_laplacian.m:38`, after the `Symmetrize` line), add:
```matlab
%     CheckManifold : (logical) Warn if the input mesh is not a clean
%                     2-manifold. Default: false. This function ASSUMES a
%                     clean 2-manifold input (e.g. icosphere-downsampled
%                     cortex); validate upstream with tess_manifold.
```

- [ ] **Step 6: Run the test to verify it passes**

Run via MATLAB MCP `run_matlab_file`: `dev/tests/test_laplacian_ico.m`
Expected: PASS — prints `ALL TESTS PASSED: test_laplacian_ico`.

- [ ] **Step 7: M-Lint**

Run MATLAB MCP `check_matlab_code` on `toolbox/anatomy/tess_laplacian.m`. Expected: no new errors.

- [ ] **Step 8: Commit**
```bash
git add toolbox/anatomy/tess_laplacian.m dev/tests/test_laplacian_ico.m
git commit -m "tess_laplacian: document manifold assumption, add optional CheckManifold warning"
```

---

## Task 3: io pair — `in_tess_eigenmodes` + `out_tess_eigenmodes`

Relocate load/save to `toolbox/io/`, rewrite them on Brainstorm io (`in_tess_bst`, single `bst_save`), and update every caller of the old names.

**Files:**
- Move + rewrite: `toolbox/anatomy/tess_eigenmodes_load.m` → `toolbox/io/in_tess_eigenmodes.m`
- Move + rewrite: `toolbox/anatomy/tess_eigenmodes_save.m` → `toolbox/io/out_tess_eigenmodes.m`
- Modify (caller, load): `process_eigenmodes.m:159`, `process_eigenmodes_filter.m:178`, `process_eigenmodes_inverse.m:158`, `process_eigenmodes_spectrum.m:151`, `process_eigenmodes_view.m:136`, `toolbox/anatomy/tess_eigenmodes_leadfield.m:172`
- Modify (caller, save): `process_eigenmodes.m:226`
- Test: `dev/tests/test_io_eigenmodes_roundtrip.m` (new)

- [ ] **Step 1: Write the failing round-trip test**

Create `dev/tests/test_io_eigenmodes_roundtrip.m`:
```matlab
function test_io_eigenmodes_roundtrip
% Verify out_tess_eigenmodes -> in_tess_eigenmodes round-trips an embedded struct.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

[V, F] = tess_sphere(162);
nV = size(V, 1);

Em = struct();
Em.Vectors     = randn(nV, 15);
Em.Values      = sort(rand(15, 1));
Em.nModes      = 15;
Em.MassType    = 'barycentric';
Em.Sigma       = -1e-8;
Em.Tolerance   = 1e-10;
Em.nRemoved    = 2;
Em.ComputeTime = 1.23;

tmpDir = tempname; mkdir(tmpDir);
cleanup = onCleanup(@() rmdir(tmpDir, 's'));

% Surface WITH eigenmodes (file must be named like a Brainstorm tess file).
SurfFile = fullfile(tmpDir, 'tess_cortex_test.mat');
bst_save(SurfFile, struct('Vertices', V, 'Faces', F, 'Comment', 'test cortex'), 'v7');
out_tess_eigenmodes(SurfFile, Em, V, F, false);

[Em2, isComputed] = in_tess_eigenmodes(SurfFile);
assert(isComputed, 'in_tess_eigenmodes reported not computed.');
assert(Em2.nModes == 15, 'nModes mismatch (%d).', Em2.nModes);
assert(isa(Em2.Vectors, 'double'), 'Vectors not converted to double.');
assert(isequal(size(Em2.Vectors), [nV, 15]), 'Vectors size mismatch.');
assert(max(abs(Em2.Values - Em.Values)) < 1e-12, 'Values changed on round-trip.');
assert(max(abs(Em2.Vectors(:) - Em.Vectors(:))) < 1e-5, 'Vectors exceeded single-precision round-trip error.');
fprintf('PASSED: round-trip preserves eigenmodes.\n');

% Surface WITHOUT eigenmodes -> isComputed false.
SurfFile2 = fullfile(tmpDir, 'tess_cortex_empty.mat');
bst_save(SurfFile2, struct('Vertices', V, 'Faces', F, 'Comment', 'empty'), 'v7');
[~, isComp2] = in_tess_eigenmodes(SurfFile2);
assert(~isComp2, 'Empty surface wrongly reported computed.');
fprintf('PASSED: missing eigenmodes report isComputed=false.\n');

fprintf('ALL TESTS PASSED: test_io_eigenmodes_roundtrip\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run via MATLAB MCP `run_matlab_file`: `dev/tests/test_io_eigenmodes_roundtrip.m`
Expected: FAIL — `out_tess_eigenmodes`/`in_tess_eigenmodes` are undefined (functions not yet at the new path).

- [ ] **Step 3: Move and rewrite the loader**
```bash
git mv toolbox/anatomy/tess_eigenmodes_load.m toolbox/io/in_tess_eigenmodes.m
```
Replace the entire body below the license header (i.e. everything from the `Eigenmodes = [];` line to the final `end`) so the function reads:
```matlab
function [Eigenmodes, isComputed] = in_tess_eigenmodes(SurfaceFile)
% IN_TESS_EIGENMODES: Load precomputed Laplace-Beltrami eigenmodes from a surface file.
%
% USAGE:  [Eigenmodes, isComputed] = in_tess_eigenmodes(SurfaceFile)
%
% DESCRIPTION:
%     Loads the embedded Eigenmodes field from a Brainstorm surface file via
%     the canonical loader in_tess_bst. Returns [] and false if eigenmodes
%     have not been computed for this surface.
%
% SEE ALSO: out_tess_eigenmodes, tess_eigenmodes, in_tess_bst
```
(keep the existing license block), then the implementation:
```matlab
Eigenmodes = [];
isComputed = false;

% Load via the canonical Brainstorm loader. isComputeMissing=0 so a
% frequently-called read does not recompute curvature/normals (those recomputes
% are gated on isComputeMissing at in_tess_bst.m:128,135,142,149).
TessMat = in_tess_bst(SurfaceFile, 0);

if ~isfield(TessMat, 'Eigenmodes') || isempty(TessMat.Eigenmodes)
    return;
end

Eigenmodes = TessMat.Eigenmodes;
if isfield(Eigenmodes, 'Vectors') && isa(Eigenmodes.Vectors, 'single')
    Eigenmodes.Vectors = double(Eigenmodes.Vectors);
end
isComputed = true;
end
```

- [ ] **Step 4: Move and rewrite the saver (single `bst_save`)**
```bash
git mv toolbox/anatomy/tess_eigenmodes_save.m toolbox/io/out_tess_eigenmodes.m
```
Replace the function declaration and SEE ALSO so the header reads:
```matlab
function TessMat = out_tess_eigenmodes(SurfaceFile, Eigenmodes, Vertices, Faces, isInteractive)
% OUT_TESS_EIGENMODES: Save precomputed eigenmodes into a Brainstorm surface file.
%
% USAGE:  TessMat = out_tess_eigenmodes(SurfaceFile, Eigenmodes, Vertices, Faces)
%         TessMat = out_tess_eigenmodes(SurfaceFile, Eigenmodes, Vertices, Faces, isInteractive)
%
% SEE ALSO: in_tess_eigenmodes, tess_eigenmodes
```
(keep the existing license block). Then replace the implementation body (everything from `if nargin < 5` to the final `end`) with:
```matlab
if nargin < 5
    isInteractive = true;
end

%% ===== VALIDATE INPUTS =====
if ~isstruct(Eigenmodes) || ~isfield(Eigenmodes, 'Vectors') || ~isfield(Eigenmodes, 'Values')
    error('Eigenmodes must be a structure with Vectors and Values fields.');
end
if size(Eigenmodes.Vectors, 1) ~= size(Vertices, 1)
    error('Eigenmodes.Vectors (%d rows) must match Vertices (%d rows).', ...
        size(Eigenmodes.Vectors, 1), size(Vertices, 1));
end

%% ===== LOAD SURFACE FILE =====
SurfaceFileFull = file_fullpath(SurfaceFile);
if ~file_exist(SurfaceFileFull)
    error('Surface file not found: %s', SurfaceFileFull);
end
TessMat = load(SurfaceFileFull);

%% ===== UPDATE MESH IF REPAIRED =====
if (size(Vertices, 1) ~= size(TessMat.Vertices, 1)) || (size(Faces, 1) ~= size(TessMat.Faces, 1))
    if isInteractive
        fprintf('BST> Mesh was repaired (%d->%d V, %d->%d F). Updating surface.\n', ...
            size(TessMat.Vertices, 1), size(Vertices, 1), size(TessMat.Faces, 1), size(Faces, 1));
    end
    TessMat.Vertices    = Vertices;
    TessMat.Faces       = Faces;
    TessMat.VertConn    = tess_vertconn(Vertices, Faces);
    TessMat.VertNormals = tess_normals(Vertices, Faces, TessMat.VertConn);
    TessMat.Curvature   = [];
    TessMat.SulciMap    = [];
    TessMat.tess2mri_interp = [];
    if isfield(TessMat, 'Atlas') && ~isempty(TessMat.Atlas)
        warning('BST:AtlasCleared', 'Atlas cleared because mesh vertices changed during repair.');
        TessMat.Atlas = [];
    end
end

%% ===== STORE EIGENMODES (single precision vectors) =====
EigenmodesStore = struct();
EigenmodesStore.Vectors     = single(Eigenmodes.Vectors);
EigenmodesStore.Values      = Eigenmodes.Values(:);
EigenmodesStore.nModes      = Eigenmodes.nModes;
EigenmodesStore.MassType    = Eigenmodes.MassType;
EigenmodesStore.Sigma       = Eigenmodes.Sigma;
EigenmodesStore.Tolerance   = Eigenmodes.Tolerance;
EigenmodesStore.nRemoved    = Eigenmodes.nRemoved;
EigenmodesStore.ComputeTime = Eigenmodes.ComputeTime;
EigenmodesStore.ComputeDate = datestr(now, 'yyyy-mm-dd HH:MM:SS');
TessMat.Eigenmodes = EigenmodesStore;

%% ===== HISTORY + SINGLE SAVE =====
TessMat = bst_history('add', TessMat, 'eigenmodes', ...
    sprintf('Computed %d Laplace-Beltrami eigenmodes (%s mass, sigma=%.1e)', ...
        Eigenmodes.nModes, Eigenmodes.MassType, Eigenmodes.Sigma));
bst_save(SurfaceFileFull, TessMat, 'v7');

if isInteractive
    fprintf('BST> Saved %d eigenmodes to: %s\n', Eigenmodes.nModes, SurfaceFile);
end
end
```

- [ ] **Step 5: Update all load call sites to `in_tess_eigenmodes`**

In each of these files, replace `tess_eigenmodes_load(` with `in_tess_eigenmodes(` on the indicated line:
- `toolbox/process/functions/process_eigenmodes.m:159` (`[~, isComputed] = tess_eigenmodes_load(SurfaceFile);`)
- `toolbox/process/functions/process_eigenmodes_filter.m:178`
- `toolbox/process/functions/process_eigenmodes_inverse.m:158`
- `toolbox/process/functions/process_eigenmodes_spectrum.m:151`
- `toolbox/process/functions/process_eigenmodes_view.m:136`
- `toolbox/anatomy/tess_eigenmodes_leadfield.m:172`

- [ ] **Step 6: Update the save call site to `out_tess_eigenmodes`**

In `toolbox/process/functions/process_eigenmodes.m:226`, replace:
```matlab
        tess_eigenmodes_save(SurfaceFile, Eigenmodes, Vertices, Faces, true);
```
with:
```matlab
        out_tess_eigenmodes(SurfaceFile, Eigenmodes, Vertices, Faces, true);
```

- [ ] **Step 7: Run the round-trip test to verify it passes**

Run via MATLAB MCP `run_matlab_file`: `dev/tests/test_io_eigenmodes_roundtrip.m`
Expected: PASS — prints `ALL TESTS PASSED: test_io_eigenmodes_roundtrip`.

- [ ] **Step 8: Verify the old load/save names are gone**
```bash
grep -rn "tess_eigenmodes_load\|tess_eigenmodes_save" toolbox/
```
Expected: no output (zero matches).

- [ ] **Step 9: M-Lint the two new io files**

Run MATLAB MCP `check_matlab_code` on `toolbox/io/in_tess_eigenmodes.m` and `toolbox/io/out_tess_eigenmodes.m`. Expected: no new errors.

- [ ] **Step 10: Commit**
```bash
git add -A toolbox/io/in_tess_eigenmodes.m toolbox/io/out_tess_eigenmodes.m \
       toolbox/process/functions/process_eigenmodes.m \
       toolbox/process/functions/process_eigenmodes_filter.m \
       toolbox/process/functions/process_eigenmodes_inverse.m \
       toolbox/process/functions/process_eigenmodes_spectrum.m \
       toolbox/process/functions/process_eigenmodes_view.m \
       toolbox/anatomy/tess_eigenmodes_leadfield.m \
       dev/tests/test_io_eigenmodes_roundtrip.m
git commit -m "Move eigenmode load/save to toolbox/io as in_/out_tess_eigenmodes"
```

---

## Task 4: `bst_eigenmodes_project` (pure) + spectrum plugin

Relocate `project` to `toolbox/math/`, make it pure (Eigenmodes + Data + MassMatrix in), and update the spectrum plugin to load the mass matrix once and pass it in.

**Files:**
- Move + rewrite: `toolbox/anatomy/tess_eigenmodes_project.m` → `toolbox/math/bst_eigenmodes_project.m`
- Modify: `toolbox/process/functions/process_eigenmodes_spectrum.m` (compute M; call new name)
- Test: `dev/tests/test_eigenmodes_project_pure.m` (new)

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigenmodes_project_pure.m`:
```matlab
function test_eigenmodes_project_pure
% Verify the pure projection recovers M-orthonormal coefficients exactly.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

rng(0);
n = 50; k = 8;
d = 0.5 + rand(n, 1);
M = spdiags(d, 0, n, n);

% Build an M-orthonormal basis: Phi = A * inv(chol(A'MA)) -> Phi'M Phi = I.
A = randn(n, k);
G = A' * (M * A);
R = chol((G + G') / 2);
Phi = A / R;
assert(norm(Phi' * (M * Phi) - eye(k), 'fro') < 1e-9, 'setup: Phi not M-orthonormal');

Em = struct('Vectors', Phi, 'Values', (1:k)', 'nModes', k, 'MassType', 'barycentric');

% Known coefficients -> data -> recover exactly (proves M weighting is applied).
c0 = randn(k, 3);
Data = Phi * c0;
[C, Recon] = bst_eigenmodes_project(Em, Data, M);
assert(max(abs(C(:) - c0(:))) < 1e-9, 'Projection did not recover coefficients (M weighting bug?).');
assert(max(abs(Recon(:) - Data(:))) < 1e-9, 'Full reconstruction not exact.');

% Parseval in the M-norm.
uMnorm2 = sum(sum(Data .* (M * Data)));
assert(abs(uMnorm2 - sum(C(:).^2)) < 1e-7, 'Parseval identity failed.');

% ModeRange reconstruction uses only the selected modes.
[~, ReconLP] = bst_eigenmodes_project(Em, Data, M, 'ModeRange', [1 3]);
expected = Phi(:, 1:3) * C(1:3, :);
assert(max(abs(ReconLP(:) - expected(:))) < 1e-9, 'ModeRange reconstruction wrong.');

fprintf('ALL TESTS PASSED: test_eigenmodes_project_pure\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run via MATLAB MCP `run_matlab_file`: `dev/tests/test_eigenmodes_project_pure.m`
Expected: FAIL — `bst_eigenmodes_project` is undefined.

- [ ] **Step 3: Move and rewrite `project` as a pure function**
```bash
git mv toolbox/anatomy/tess_eigenmodes_project.m toolbox/math/bst_eigenmodes_project.m
```
Replace the function declaration + USAGE + SEE ALSO so the header reads:
```matlab
function [Coeffs, Reconstructed] = bst_eigenmodes_project(Eigenmodes, Data, MassMatrix, varargin)
% BST_EIGENMODES_PROJECT: Project scalar data onto a precomputed eigenmode basis.
%
% USAGE:  Coeffs = bst_eigenmodes_project(Eigenmodes, Data, MassMatrix)
%         Coeffs = bst_eigenmodes_project(Eigenmodes, Data, MassMatrix, 'ModeRange', [1 100])
%         [Coeffs, Reconstructed] = bst_eigenmodes_project(Eigenmodes, Data, MassMatrix, ...)
%
% INPUTS:
%     Eigenmodes : struct from in_tess_eigenmodes (uses .Vectors)
%     Data       : [nVertices x nTime] scalar field(s) on the mesh
%     MassMatrix : [nVertices x nVertices] sparse mass matrix (from tess_laplacian)
%
% OPTIONS:
%     ModeRange  : [1 x 2] mode index range for reconstruction (default: all)
%
% SEE ALSO: tess_eigenmodes, in_tess_eigenmodes, bst_eigenmodes_filter, tess_laplacian
```
(keep the existing license block). Then replace the implementation body (everything from the first `%% ===== PARSE INPUTS =====` to the final `end`) with:
```matlab
%% ===== PARSE INPUTS =====
ModeRange = [];
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'moderange', ModeRange = varargin{i+1};
    end
end

Phi    = double(Eigenmodes.Vectors);   % [nV x nModes]
nV     = size(Phi, 1);
nModes = size(Phi, 2);

%% ===== VALIDATE =====
Data = double(Data);
if size(Data, 1) ~= nV
    error('Data has %d rows but eigenmodes have %d vertices.', size(Data, 1), nV);
end
if (size(MassMatrix, 1) ~= nV) || (size(MassMatrix, 2) ~= nV)
    error('MassMatrix must be %dx%d.', nV, nV);
end

%% ===== MODE RANGE =====
if isempty(ModeRange)
    ModeRange = [1, nModes];
end
ModeRange(1) = max(1, ModeRange(1));
ModeRange(2) = min(nModes, ModeRange(2));
iModes = ModeRange(1):ModeRange(2);
if isempty(iModes)
    error('Empty mode range [%d, %d] (have %d modes).', ModeRange(1), ModeRange(2), nModes);
end

%% ===== PROJECT: c_k = phi_k' * M * u =====
Coeffs = Phi' * (MassMatrix * Data);   % [nModes x nTime]

%% ===== RECONSTRUCT (if requested) =====
if nargout >= 2
    Reconstructed = Phi(:, iModes) * Coeffs(iModes, :);   % [nV x nTime]
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run via MATLAB MCP `run_matlab_file`: `dev/tests/test_eigenmodes_project_pure.m`
Expected: PASS — prints `ALL TESTS PASSED: test_eigenmodes_project_pure`.

- [ ] **Step 5: Update the spectrum plugin to compute M and call the new name**

In `toolbox/process/functions/process_eigenmodes_spectrum.m`, after the vertex-count check (immediately after the `if nV_data ~= nV_eigen ... end` block ending at `:165`), insert:
```matlab

        % Mass matrix for the M-weighted projection (computed once per surface).
        sSurf = in_tess_bst(SurfaceFile, 0);
        [~, MassMatrix] = tess_laplacian(sSurf.Vertices, sSurf.Faces, 'MassType', Eigenmodes.MassType);
```
Then replace the projection call at `:171`:
```matlab
        Coeffs = tess_eigenmodes_project(SurfaceFile, double(Sources));
```
with:
```matlab
        Coeffs = bst_eigenmodes_project(Eigenmodes, double(Sources), MassMatrix);
```

- [ ] **Step 6: M-Lint and verify no stale references**

Run MATLAB MCP `check_matlab_code` on `toolbox/math/bst_eigenmodes_project.m` and `toolbox/process/functions/process_eigenmodes_spectrum.m`. Expected: no new errors.
```bash
grep -rn "tess_eigenmodes_project" toolbox/
```
Expected: no output.

- [ ] **Step 7: Commit**
```bash
git add -A toolbox/math/bst_eigenmodes_project.m \
       toolbox/process/functions/process_eigenmodes_spectrum.m \
       dev/tests/test_eigenmodes_project_pure.m
git commit -m "Move eigenmode projection to toolbox/math as pure bst_eigenmodes_project"
```

---

## Task 5: `bst_eigenmodes_filter` (pure) + filter plugin

Relocate `filter` to `toolbox/math/`, make it pure, and update the filter plugin to compute the mass matrix once (fixing the per-frequency recompute).

**Files:**
- Move + rewrite: `toolbox/anatomy/tess_eigenmodes_filter.m` → `toolbox/math/bst_eigenmodes_filter.m`
- Modify: `toolbox/process/functions/process_eigenmodes_filter.m` (compute M once; call new name)
- Test: `dev/tests/test_eigenmodes_filter_pure.m` (new)

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigenmodes_filter_pure.m`:
```matlab
function test_eigenmodes_filter_pure
% Verify the pure spectral filter (lowpass + heat kernel limits).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

rng(1);
n = 40; k = 10;
d = 0.5 + rand(n, 1);
M = spdiags(d, 0, n, n);
A = randn(n, k);
G = A' * (M * A);
R = chol((G + G') / 2);
Phi = A / R;

lambdas = (0:k-1)'.^2;   % eigenvalues; lambda(1)=0 is the DC mode
Em = struct('Vectors', Phi, 'Values', lambdas, 'nModes', k, 'MassType', 'barycentric');

c0 = randn(k, 1);
Data = Phi * c0;

% Lowpass cutoff 4 keeps modes 1..4 only.
LP = bst_eigenmodes_filter(Em, Data, M, 'lowpass', 'CutoffMode', 4);
assert(max(abs(LP - Phi(:, 1:4) * c0(1:4))) < 1e-9, 'lowpass mismatch.');

% Heat kernel with t -> 0 is the identity.
H0 = bst_eigenmodes_filter(Em, Data, M, 'heat', 'DiffusionTime', 1e-12);
assert(max(abs(H0 - Data)) < 1e-6, 'heat t->0 is not identity.');

% Heat kernel with large t isolates the DC mode (lambda=0 -> gain 1; rest -> 0).
HT = bst_eigenmodes_filter(Em, Data, M, 'heat', 'DiffusionTime', 1e6);
assert(max(abs(HT - Phi(:, 1) * c0(1))) < 1e-6, 'heat large t did not isolate DC.');

fprintf('ALL TESTS PASSED: test_eigenmodes_filter_pure\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run via MATLAB MCP `run_matlab_file`: `dev/tests/test_eigenmodes_filter_pure.m`
Expected: FAIL — `bst_eigenmodes_filter` is undefined.

- [ ] **Step 3: Move and rewrite `filter` as a pure function**
```bash
git mv toolbox/anatomy/tess_eigenmodes_filter.m toolbox/math/bst_eigenmodes_filter.m
```
Replace the function declaration + USAGE + SEE ALSO so the header reads:
```matlab
function Filtered = bst_eigenmodes_filter(Eigenmodes, Data, MassMatrix, FilterType, varargin)
% BST_EIGENMODES_FILTER: Spatial spectral filtering of cortical data via eigenmodes.
%
% USAGE:  Filtered = bst_eigenmodes_filter(Eigenmodes, Data, MassMatrix, 'lowpass',  'CutoffMode', 50)
%         Filtered = bst_eigenmodes_filter(Eigenmodes, Data, MassMatrix, 'bandpass', 'ModeRange', [20 80])
%         Filtered = bst_eigenmodes_filter(Eigenmodes, Data, MassMatrix, 'heat',     'DiffusionTime', 0.01)
%
% INPUTS:
%     Eigenmodes : struct from in_tess_eigenmodes (uses .Vectors, .Values)
%     Data       : [nVertices x nTime] scalar field(s) on the mesh
%     MassMatrix : [nVertices x nVertices] sparse mass matrix (from tess_laplacian)
%     FilterType : 'lowpass','highpass','bandpass','heat','inverse_heat','custom'
%
% OPTIONS:
%     CutoffMode, ModeRange, DiffusionTime, MaxGain, TransferFn (see code)
%
% SEE ALSO: bst_eigenmodes_project, in_tess_eigenmodes, tess_eigenmodes, tess_laplacian
```
(keep the existing license block). Then replace the implementation body (everything from the first `%% ===== PARSE INPUTS =====` to the final `end`) with:
```matlab
%% ===== PARSE INPUTS =====
CutoffMode    = 50;
ModeRange     = [20, 80];
DiffusionTime = 0.01;
MaxGain       = 10;
TransferFn    = [];
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'cutoffmode',    CutoffMode    = varargin{i+1};
        case 'moderange',     ModeRange     = varargin{i+1};
        case 'diffusiontime', DiffusionTime = varargin{i+1};
        case 'maxgain',       MaxGain       = varargin{i+1};
        case 'transferfn',    TransferFn    = varargin{i+1};
    end
end

Phi     = double(Eigenmodes.Vectors);   % [nV x nModes]
lambdas = double(Eigenmodes.Values(:));  % [nModes x 1]
nV      = size(Phi, 1);
nModes  = size(Phi, 2);

%% ===== VALIDATE =====
Data = double(Data);
if size(Data, 1) ~= nV
    error('Data has %d rows but eigenmodes have %d vertices.', size(Data, 1), nV);
end
if (size(MassMatrix, 1) ~= nV) || (size(MassMatrix, 2) ~= nV)
    error('MassMatrix must be %dx%d.', nV, nV);
end

%% ===== BUILD TRANSFER FUNCTION =====
h = zeros(nModes, 1);
switch lower(FilterType)
    case 'lowpass'
        c = min(CutoffMode, nModes);
        h(1:c) = 1;
    case 'highpass'
        c = max(1, min(CutoffMode, nModes));
        h(c:end) = 1;
    case 'bandpass'
        k1 = max(1, ModeRange(1));
        k2 = min(nModes, ModeRange(2));
        h(k1:k2) = 1;
    case 'heat'
        if DiffusionTime <= 0
            error('DiffusionTime must be positive (got %g).', DiffusionTime);
        end
        h = exp(-DiffusionTime * lambdas);
    case 'inverse_heat'
        if DiffusionTime <= 0
            error('DiffusionTime must be positive (got %g).', DiffusionTime);
        end
        h = min(exp(DiffusionTime * lambdas), MaxGain);
    case 'custom'
        if isempty(TransferFn) || ~isa(TransferFn, 'function_handle')
            error('Custom filter requires a TransferFn option (function handle).');
        end
        h = TransferFn(lambdas);
        if numel(h) ~= nModes
            error('TransferFn must return a vector of length %d (got %d).', nModes, numel(h));
        end
        h = h(:);
    otherwise
        error('Unknown filter type: %s. Use lowpass, highpass, bandpass, heat, inverse_heat, or custom.', FilterType);
end

%% ===== APPLY: u_filtered = Phi * (h .* (Phi' * M * u)) =====
Coeffs   = Phi' * (MassMatrix * Data);
Coeffs   = bsxfun(@times, h, Coeffs);
Filtered = Phi * Coeffs;
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run via MATLAB MCP `run_matlab_file`: `dev/tests/test_eigenmodes_filter_pure.m`
Expected: PASS — prints `ALL TESTS PASSED: test_eigenmodes_filter_pure`.

- [ ] **Step 5: Update the filter plugin to compute M once and call the new name**

In `toolbox/process/functions/process_eigenmodes_filter.m`, after the vertex-count check (immediately after the `if nV_data ~= nV_eigen ... end` block ending at `:198`), insert:
```matlab

    % Mass matrix for the M-weighted filter (computed once, reused across freqs).
    sSurf = in_tess_bst(SurfaceFile, 0);
    [~, MassMatrix] = tess_laplacian(sSurf.Vertices, sSurf.Faces, 'MassType', Eigenmodes.MassType);
```
Then replace the frequency-loop call (`:231-232`):
```matlab
            sInput.A(:, :, iFreq) = tess_eigenmodes_filter(SurfaceFile, ...
                sInput.A(:, :, iFreq), FilterType, filterArgs{:});
```
with:
```matlab
            sInput.A(:, :, iFreq) = bst_eigenmodes_filter(Eigenmodes, ...
                sInput.A(:, :, iFreq), MassMatrix, FilterType, filterArgs{:});
```
and replace the single-band call (`:235`):
```matlab
        sInput.A = tess_eigenmodes_filter(SurfaceFile, sInput.A, FilterType, filterArgs{:});
```
with:
```matlab
        sInput.A = bst_eigenmodes_filter(Eigenmodes, sInput.A, MassMatrix, FilterType, filterArgs{:});
```

- [ ] **Step 6: M-Lint and verify no stale references**

Run MATLAB MCP `check_matlab_code` on `toolbox/math/bst_eigenmodes_filter.m` and `toolbox/process/functions/process_eigenmodes_filter.m`. Expected: no new errors.
```bash
grep -rn "tess_eigenmodes_filter" toolbox/
```
Expected: no output.

- [ ] **Step 7: Commit**
```bash
git add -A toolbox/math/bst_eigenmodes_filter.m \
       toolbox/process/functions/process_eigenmodes_filter.m \
       dev/tests/test_eigenmodes_filter_pure.m
git commit -m "Move eigenmode filter to toolbox/math as pure bst_eigenmodes_filter"
```

---

## Task 6: `bst_inverse_eigenmodes` (relocate to toolbox/inverse)

Relocate the eigenmode-space inverse solver. The algorithm is unchanged; the eigenmode load was already switched to `in_tess_eigenmodes` in Task 3.

**Files:**
- Move + rename: `toolbox/anatomy/tess_eigenmodes_leadfield.m` → `toolbox/inverse/bst_inverse_eigenmodes.m`
- Modify: `toolbox/process/functions/process_eigenmodes_inverse.m:223` (call new name)

> **No standalone unit test:** this function needs a head model, channel file, and noise covariance, so it is validated by M-Lint, the no-stale-references grep, and the optional OMEGA regression test (`dev/tests/test_omega_icosphere_sourcemap.m`) in Task 7.

- [ ] **Step 1: Move and rename the function**
```bash
git mv toolbox/anatomy/tess_eigenmodes_leadfield.m toolbox/inverse/bst_inverse_eigenmodes.m
```
Replace the function declaration (`:1`):
```matlab
function [Results, errMsg] = tess_eigenmodes_leadfield(HeadModelFile, SurfaceFile, NoiseCovFile, varargin)
% TESS_EIGENMODES_LEADFIELD: Eigenmode-space source mapping via compressed lead field.
```
with:
```matlab
function [Results, errMsg] = bst_inverse_eigenmodes(HeadModelFile, SurfaceFile, NoiseCovFile, varargin)
% BST_INVERSE_EIGENMODES: Eigenmode-space source mapping via compressed lead field.
```

- [ ] **Step 2: Update the SEE ALSO line**

In the header (`tess_eigenmodes_leadfield.m:98`), replace:
```matlab
% SEE ALSO: tess_eigenmodes_load, tess_eigenmodes_project, process_eigenmodes_inverse
```
with:
```matlab
% SEE ALSO: in_tess_eigenmodes, bst_eigenmodes_project, process_eigenmodes_inverse
```

- [ ] **Step 3: Update the caller in the inverse plugin**

In `toolbox/process/functions/process_eigenmodes_inverse.m:223`, replace:
```matlab
    [InvResults, errMsg] = tess_eigenmodes_leadfield(HeadModelFile, SurfaceFile, NoiseCovFile, ...
```
with:
```matlab
    [InvResults, errMsg] = bst_inverse_eigenmodes(HeadModelFile, SurfaceFile, NoiseCovFile, ...
```

- [ ] **Step 4: M-Lint and verify no stale references**

Run MATLAB MCP `check_matlab_code` on `toolbox/inverse/bst_inverse_eigenmodes.m` and `toolbox/process/functions/process_eigenmodes_inverse.m`. Expected: no new errors.
```bash
grep -rn "tess_eigenmodes_leadfield" toolbox/
```
Expected: no output.

- [ ] **Step 5: Commit**
```bash
git add -A toolbox/inverse/bst_inverse_eigenmodes.m \
       toolbox/process/functions/process_eigenmodes_inverse.m
git commit -m "Move eigenmode-space inverse to toolbox/inverse as bst_inverse_eigenmodes"
```

---

## Task 7: Final verification sweep

Confirm no stale references remain anywhere, fix lingering SEE ALSO comments, M-Lint all touched files, and run the full numerical + structural test suite together.

**Files:**
- Modify (comments only, as needed): `toolbox/process/functions/process_eigenmodes*.m` SEE ALSO lines
- Modify: `toolbox/anatomy/tess_eigenmodes.m` SEE ALSO if it names moved functions

- [ ] **Step 1: Repo-wide grep for all five old names**
```bash
grep -rn "tess_eigenmodes_load\|tess_eigenmodes_save\|tess_eigenmodes_project\|tess_eigenmodes_filter\|tess_eigenmodes_leadfield" toolbox/ dev/
```
Expected: **no output**. If any remain (most likely `SEE ALSO` comment lines in `process_eigenmodes*.m` headers), update them to the new names:
- `tess_eigenmodes_load` → `in_tess_eigenmodes`
- `tess_eigenmodes_save` → `out_tess_eigenmodes`
- `tess_eigenmodes_project` → `bst_eigenmodes_project`
- `tess_eigenmodes_filter` → `bst_eigenmodes_filter`
- `tess_eigenmodes_leadfield` → `bst_inverse_eigenmodes`

- [ ] **Step 2: Re-run the grep to confirm zero matches**
```bash
grep -rn "tess_eigenmodes_load\|tess_eigenmodes_save\|tess_eigenmodes_project\|tess_eigenmodes_filter\|tess_eigenmodes_leadfield" toolbox/ dev/
```
Expected: no output.

- [ ] **Step 3: Run the full new test suite**

Run each via MATLAB MCP `run_matlab_file`, confirming each prints its `ALL TESTS PASSED` line:
- `dev/tests/test_eigenmodes_manifold_gate.m`
- `dev/tests/test_process_eigenmodes_options.m`
- `dev/tests/test_laplacian_ico.m`
- `dev/tests/test_io_eigenmodes_roundtrip.m`
- `dev/tests/test_eigenmodes_project_pure.m`
- `dev/tests/test_eigenmodes_filter_pure.m`

- [ ] **Step 4: (Optional) Run the OMEGA icosphere regression test**

If a `TutorialOmega` dataset is available locally, run `dev/tests/test_omega_icosphere_sourcemap.m` via the MATLAB MCP to confirm the full cortex → head model → eigenmode-inverse chain still works end-to-end. Skip if the dataset is not present.

- [ ] **Step 5: Commit any comment fixes**
```bash
git add -A
git commit -m "Update SEE ALSO references to relocated eigenmode functions"
```
(If Step 1 produced no output and no edits were needed, skip this commit.)

- [ ] **Step 6: Final review**

Confirm the branch is clean and the seven functions are in their final homes:
```bash
git status
ls toolbox/anatomy/tess_laplacian.m toolbox/anatomy/tess_eigenmodes.m \
   toolbox/io/in_tess_eigenmodes.m toolbox/io/out_tess_eigenmodes.m \
   toolbox/math/bst_eigenmodes_project.m toolbox/math/bst_eigenmodes_filter.m \
   toolbox/inverse/bst_inverse_eigenmodes.m
```
Expected: all seven paths exist; working tree clean.

---

## Self-review notes (for the implementer)

- **Order matters.** Task 3 (io) must precede Tasks 4–6 so that `in_tess_eigenmodes` exists before the plugins/inverse reference it. Within Task 3, the load-rename also touches `tess_eigenmodes_leadfield.m:172` even though that file is renamed later — that is intentional.
- **Transient state.** After Task 1, `process_eigenmodes` still calls the old `tess_eigenmodes_load`/`_save` names (renamed in Task 3); that is fine because those files still exist until Task 3.
- **`tess_manifold` report.** The gate uses `report.summary`; the code guards `isfield(report,'summary')` so it degrades gracefully if the field name differs. Verify against `toolbox/anatomy/tess_manifold.m` during Task 1 and adjust the field name if needed.
- **`java_dialog('question', ...)`** returns the chosen button label as a string (empty on cancel); see `toolbox/sensors/channel_fixunits.m:86` for the argument order.
