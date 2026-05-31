# nxr-compute Brainstorm Plugin — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Package the existing nxr-compute MATLAB MEX binding as a native, SPM-style Brainstorm plugin (macOS first) and prove it end-to-end by making `tess_laplacian` nxr-first with a MATLAB fallback.

**Architecture:** nxr-compute already ships a `+nxr` wrapper + `nxr_compute.mex*` dispatcher. We (1) package wrapper+binary as a plugin zip, (2) register a per-OS `PlugDesc` in `bst_plugin.m`, and (3) wire `tess_laplacian` to use `nxr.manifold.*` when the plugin is loaded, falling back to the existing MATLAB code otherwise. The Brainstorm side targets the stable `nxr.manifold.*` wrapper, never the raw mex dispatcher. No code is added to the nxr-compute repo except a packaging script — nxr stays application-agnostic.

**Tech Stack:** MATLAB R2023b, Brainstorm (`bst_plugin.m`), nxr-compute C++/MEX (geometry-central + Eigen + Spectra), bash packaging, MATLAB MCP for test execution.

**Two repos are touched.** Each task states its repo and working directory:
- `BST = /Users/diellorbasha/workspace/research/code/brainstorm3` (branch `feature/nxr-compute-plugin`)
- `NXR = /Users/diellorbasha/workspace/research/code/nxr-compute`

**Key facts already verified (do not re-derive):**
- The mex `mxToFaceBuffer` **subtracts 1** (1-based MATLAB → 0-based C). So pass Brainstorm's native **1-based `Faces` directly** to `nxr.manifold.context`. Do NOT subtract 1. (The `+nxr/.../context.m` docstring saying "zero-based" is a doc bug; the marshalling governs.)
- `bst_plugin('Install', Name, isInteractive)` installs **and** loads, returns `[isOk, errMsg, PlugDesc]`, and downloads `URLzip` over http. `bst_plugin('Load', Name)` / `bst_plugin('Unload', Name)` toggle path membership. `bst_plugin('GetInstalled', Name)` returns a `PlugDesc` with field `.isLoaded`.
- A plugin is discovered offline when `<bst_get('UserPluginsDir')>/nxr-compute/` contains the `TestFile`. So wiring tests stage the locally-built plugin there and `Load` it — no GitHub release required. The network `Install` path is validated later, post-publish (Task 8).
- Test idiom (see `dev/tests/test_laplacian_ico.m`): function-form script, `addpath(repoRoot); if ~brainstorm('status'), brainstorm nogui; end`, build mesh with `[V,F]=tess_sphere(642)`, `assert(...)`, print `PASSED`.
- The mex serves **Voronoi mass only** (`assembleManifoldOperators` takes exactly `(V,F)`); stiffness `L` is mass-independent and available for all variants.

---

## File Structure

**nxr-compute repo (`NXR`):**
- Create: `NXR/scripts/package-plugin.sh` — builds (if needed) and assembles `dist/plugin/nxr-compute/` (wrapper + mex + README) and a versioned zip.

**Brainstorm repo (`BST`):**
- Modify: `BST/toolbox/core/bst_plugin.m` — add the `nxr-compute` `PlugDesc` entry in `GetSupported()`.
- Modify: `BST/toolbox/anatomy/tess_laplacian.m` — extract `local_mass_matlab`, add `nxr_is_loaded` + `resolve_backend` helpers and the `Backend` option + nxr path.
- Create: `BST/dev/tests/test_nxr_plugin_lifecycle.m` — stage + load + compute + unload.
- Create: `BST/dev/tests/test_laplacian_nxr_parity.m` — nxr vs MATLAB `L`/`M` parity.
- Create: `BST/dev/tests/test_laplacian_backend_select.m` — auto/forced backend selection + fallback.

---

## Task 1: Packaging script (NXR repo)

**Repo:** NXR · **cwd:** `/Users/diellorbasha/workspace/research/code/nxr-compute`

**Files:**
- Create: `NXR/scripts/package-plugin.sh`

- [ ] **Step 1: Write the packaging script**

Create `NXR/scripts/package-plugin.sh`:

```bash
#!/usr/bin/env bash
# Package the nxr-compute MATLAB plugin (wrapper + mex) for Brainstorm.
# Usage: scripts/package-plugin.sh [version] [os-tag]
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSION="${1:-dev}"
OS_TAG="${2:-mac}"
STAGE_ROOT="$REPO/dist/plugin"
STAGE="$STAGE_ROOT/nxr-compute"
ZIP="$REPO/dist/nxr-compute-${VERSION}-${OS_TAG}.zip"

# 1. Ensure the mex is built
if ! ls "$REPO"/build/Release/nxr_compute.mex* >/dev/null 2>&1; then
  echo "No mex found; building (Release)..."
  bash "$REPO/scripts/build.sh" Release
fi

# 2. Assemble staging layout: +nxr, mex, README at the root
rm -rf "$STAGE"
mkdir -p "$STAGE"
cp -R "$REPO/bindings/mex/matlab/+nxr" "$STAGE/"
cp "$REPO/bindings/mex/matlab/README.md" "$STAGE/"
cp "$REPO"/build/Release/nxr_compute.mex* "$STAGE/"

# 3. Zip (folder-rooted so it extracts to nxr-compute/)
rm -f "$ZIP"
( cd "$STAGE_ROOT" && zip -rq "$ZIP" "nxr-compute" )

echo "Packaged: $ZIP"
echo "Staged:   $STAGE"
unzip -l "$ZIP"
```

- [ ] **Step 2: Make it executable and run it**

Run:
```bash
chmod +x /Users/diellorbasha/workspace/research/code/nxr-compute/scripts/package-plugin.sh
/Users/diellorbasha/workspace/research/code/nxr-compute/scripts/package-plugin.sh dev mac
```
Expected: prints `Packaged: .../dist/nxr-compute-dev-mac.zip` and an `unzip -l` listing that includes `nxr-compute/+nxr/`, `nxr-compute/nxr_compute.mexmaca64`, and `nxr-compute/README.md`.

- [ ] **Step 3: Verify the staged layout has the wrapper and binary side by side**

Run:
```bash
ls /Users/diellorbasha/workspace/research/code/nxr-compute/dist/plugin/nxr-compute
test -d /Users/diellorbasha/workspace/research/code/nxr-compute/dist/plugin/nxr-compute/+nxr && \
test -f /Users/diellorbasha/workspace/research/code/nxr-compute/dist/plugin/nxr-compute/nxr_compute.mexmaca64 && \
echo "STAGE OK"
```
Expected: lists `+nxr  README.md  nxr_compute.mexmaca64` and prints `STAGE OK`.

- [ ] **Step 4: Commit (NXR repo)**

```bash
cd /Users/diellorbasha/workspace/research/code/nxr-compute
git add scripts/package-plugin.sh
git commit -m "Add Brainstorm plugin packaging script (wrapper + mex zip)"
```
> Note: `dist/` is build output — confirm it is git-ignored (it is, per `.gitignore`); do not commit the zip/staging.

---

## Task 2: Register the `PlugDesc` entry (BST repo)

**Repo:** BST · **cwd:** `/Users/diellorbasha/workspace/research/code/brainstorm3`

**Files:**
- Modify: `BST/toolbox/core/bst_plugin.m` (inside `GetSupported()`)
- Test: `BST/dev/tests/test_nxr_plugin_registered.m`

- [ ] **Step 1: Write the failing test**

Create `BST/dev/tests/test_nxr_plugin_registered.m`:

```matlab
function test_nxr_plugin_registered
% Verify nxr-compute is registered in bst_plugin GetSupported with expected fields.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

PlugDesc = bst_plugin('GetSupported', 'nxr-compute');
assert(~isempty(PlugDesc) && strcmp(PlugDesc.Name, 'nxr-compute'), ...
    'nxr-compute not found in GetSupported.');
assert(strcmp(PlugDesc.Category, 'Anatomy'), 'Unexpected Category: %s', PlugDesc.Category);
assert(PlugDesc.AutoLoad == 0, 'nxr-compute must be AutoLoad=0 (install-on-demand).');
assert(~isempty(PlugDesc.URLinfo), 'URLinfo must be set.');
% On Apple Silicon the TestFile must point at the mac arm binary
if strcmp(bst_get('OsType'), 'mac64arm')
    assert(strcmp(PlugDesc.TestFile, 'nxr_compute.mexmaca64'), ...
        'TestFile must be nxr_compute.mexmaca64 on mac64arm, got: %s', PlugDesc.TestFile);
end
fprintf('ALL TESTS PASSED: test_nxr_plugin_registered\n');
end
```

- [ ] **Step 2: Run it to confirm it fails**

Run (MATLAB MCP `run_matlab_test_file`): `dev/tests/test_nxr_plugin_registered.m`
Expected: FAIL — `nxr-compute not found in GetSupported` (entry not added yet).

- [ ] **Step 3: Add the `PlugDesc` entry**

In `BST/toolbox/core/bst_plugin.m`, inside `GetSupported()`, add this block alongside the other Anatomy-category plugins (e.g. near the `iso2mesh`/`brain2mesh` entries). Insert after an existing `PlugDesc(end+1) = GetStruct(...)` block:

```matlab
    % === ANATOMY: NXR-COMPUTE (geometry compute backend) ===
    PlugDesc(end+1)              = GetStruct('nxr-compute');
    PlugDesc(end).Version        = 'github-master';
    PlugDesc(end).Category       = 'Anatomy';
    PlugDesc(end).AutoUpdate     = 0;
    PlugDesc(end).AutoLoad       = 0;            % SPM-style install-on-demand
    PlugDesc(end).URLinfo        = 'https://github.com/neurodynamics-xr/nxr-compute';
    PlugDesc(end).ReadmeFile     = 'README.md';
    PlugDesc(end).CompiledStatus = 1;            % native code, download-only
    switch bst_get('OsType')
        case 'mac64arm'
            PlugDesc(end).URLzip   = 'https://github.com/neurodynamics-xr/nxr-compute/releases/download/plugin-dev/nxr-compute-dev-mac.zip';
            PlugDesc(end).TestFile = 'nxr_compute.mexmaca64';
        % 'linux64' / 'win64' / 'mac64' arms added when those binaries exist
    end
```

> The `URLzip` is the eventual release asset (set for real in Task 8). It is unused by the offline wiring tests, which `Load` a locally-staged copy. `Version='github-master'` keeps `AutoUpdate=0` from forcing version comparisons.

- [ ] **Step 4: Run the test to confirm it passes**

Run (MATLAB MCP `run_matlab_test_file`): `dev/tests/test_nxr_plugin_registered.m`
Expected: PASS — `ALL TESTS PASSED: test_nxr_plugin_registered`.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/core/bst_plugin.m dev/tests/test_nxr_plugin_registered.m
git commit -m "nxr-compute plugin: register PlugDesc in bst_plugin (mac arm)"
```

---

## Task 3: Plugin lifecycle test — offline stage + load + compute + unload (BST repo)

**Repo:** BST · **Prereq:** Task 1 staging dir exists (`NXR/dist/plugin/nxr-compute`).

**Files:**
- Create: `BST/dev/tests/test_nxr_plugin_lifecycle.m`

- [ ] **Step 1: Write the lifecycle test**

Create `BST/dev/tests/test_nxr_plugin_lifecycle.m`:

```matlab
function test_nxr_plugin_lifecycle
% Stage the locally-built nxr-compute plugin, load it, run a compute, unload.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));   % .../brainstorm3
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% Locate the staged plugin built by NXR/scripts/package-plugin.sh
codeDir  = fileparts(repoRoot);             % .../research/code
nxrStage = fullfile(codeDir, 'nxr-compute', 'dist', 'plugin', 'nxr-compute');
assert(isfolder(nxrStage), ...
    'Staged plugin not found at %s. Run nxr-compute/scripts/package-plugin.sh first.', nxrStage);

% Stage into the Brainstorm user plugins dir (offline "install")
destDir = fullfile(bst_get('UserPluginsDir'), 'nxr-compute');
if isfolder(destDir)
    file_delete(destDir, 1, 3);
end
mkdir(destDir);
copyfile(fullfile(nxrStage, '*'), destDir);

% Load and verify it is reported as loaded
[isOk, errMsg] = bst_plugin('Load', 'nxr-compute');
assert(isOk, 'bst_plugin Load failed: %s', errMsg);
PlugDesc = bst_plugin('GetInstalled', 'nxr-compute');
assert(~isempty(PlugDesc) && PlugDesc.isLoaded, 'Plugin not reported as loaded.');
fprintf('PASSED: nxr-compute staged and loaded.\n');

% Prove the binary actually computes (context builds K and M)
[V, F] = tess_sphere(162);
ctx = nxr.manifold.context(V, F);
assert(isfield(ctx, 'K') && size(ctx.K, 1) == size(V, 1), 'nxr stiffness K wrong size.');
assert(isfield(ctx, 'M') && size(ctx.M, 1) == size(V, 1), 'nxr mass M wrong size.');
fprintf('PASSED: nxr.manifold.context returns K/M (%dx%d).\n', size(ctx.K,1), size(ctx.K,2));

% Unload restores state
bst_plugin('Unload', 'nxr-compute');
PlugDesc2 = bst_plugin('GetInstalled', 'nxr-compute');
assert(isempty(PlugDesc2) || ~PlugDesc2.isLoaded, 'Plugin still loaded after Unload.');
fprintf('ALL TESTS PASSED: test_nxr_plugin_lifecycle\n');
end
```

- [ ] **Step 2: Run it**

Run (MATLAB MCP `run_matlab_test_file`): `dev/tests/test_nxr_plugin_lifecycle.m`
Expected: PASS — three `PASSED` lines and `ALL TESTS PASSED: test_nxr_plugin_lifecycle`.
If it fails with a faces/manifold error, that is the convention landmine — confirm `tess_sphere` faces are passed 1-based unchanged (they must be) before proceeding.

- [ ] **Step 3: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add dev/tests/test_nxr_plugin_lifecycle.m
git commit -m "nxr-compute plugin: offline lifecycle test (stage/load/compute/unload)"
```

---

## Task 4: Refactor `tess_laplacian` — extract `local_mass_matlab` (no behavior change)

**Repo:** BST · This is a pure refactor; the existing `test_laplacian_ico.m` is the guard.

**Files:**
- Modify: `BST/toolbox/anatomy/tess_laplacian.m`

- [ ] **Step 1: Confirm the guard test passes before refactor**

Run (MATLAB MCP `run_matlab_test_file`): `dev/tests/test_laplacian_ico.m`
Expected: PASS — `ALL TESTS PASSED: test_laplacian_ico`.

- [ ] **Step 2: Add the `local_mass_matlab` local function**

At the end of `BST/toolbox/anatomy/tess_laplacian.m`, before the final `end`, add this local function (verbatim extraction of the existing mass logic; recomputes the geometric primitives it needs so it is self-contained):

```matlab
function M = local_mass_matlab(Vertices, Faces, MassType)
% Pure-MATLAB mass matrix for a triangle mesh (barycentric/voronoi/galerkin).
nV = size(Vertices, 1);
Faces = double(Faces);
V1 = Vertices(Faces(:,1), :);
V2 = Vertices(Faces(:,2), :);
V3 = Vertices(Faces(:,3), :);
e1 = V2 - V1;  e2 = V3 - V1;  e3 = V3 - V2;
dot12   = sum(e1 .* e2, 2);
cross12 = cross(e1, e2, 2);
dblArea = sqrt(sum(cross12.^2, 2));
dot_at2 = sum((-e1) .* e3, 2);
dot_at3 = sum((-e2) .* (-e3), 2);
denom = dblArea;  denom(denom == 0) = eps;
areas = 0.5 * dblArea;

switch MassType
    case 'barycentric'
        vertAreas = accumarray(Faces(:), repmat(areas, 3, 1) / 3, [nV, 1]);
        M = sparse(1:nV, 1:nV, vertAreas, nV, nV);

    case 'voronoi'
        l23sq = sum(e3.^2, 2);
        l13sq = sum(e2.^2, 2);
        l12sq = sum(e1.^2, 2);
        obtuse1 = (dot12 < 0);
        obtuse2 = (dot_at2 < 0);
        obtuse3 = (dot_at3 < 0);
        anyObtuse = obtuse1 | obtuse2 | obtuse3;
        raw_cot1 = dot12 ./ denom;
        raw_cot2 = dot_at2 ./ denom;
        raw_cot3 = dot_at3 ./ denom;
        areaV1 = (1/8) * (raw_cot3 .* l12sq + raw_cot2 .* l13sq);
        areaV2 = (1/8) * (raw_cot3 .* l12sq + raw_cot1 .* l23sq);
        areaV3 = (1/8) * (raw_cot2 .* l13sq + raw_cot1 .* l23sq);
        areaV1(obtuse1) = areas(obtuse1) / 2;
        areaV2(obtuse2) = areas(obtuse2) / 2;
        areaV3(obtuse3) = areas(obtuse3) / 2;
        otherObtuse1 = anyObtuse & ~obtuse1;
        otherObtuse2 = anyObtuse & ~obtuse2;
        otherObtuse3 = anyObtuse & ~obtuse3;
        areaV1(otherObtuse1) = areas(otherObtuse1) / 4;
        areaV2(otherObtuse2) = areas(otherObtuse2) / 4;
        areaV3(otherObtuse3) = areas(otherObtuse3) / 4;
        vertAreas = accumarray(Faces(:), [areaV1; areaV2; areaV3], [nV, 1]);
        M = sparse(1:nV, 1:nV, vertAreas, nV, nV);

    case 'galerkin'
        ii_diag = Faces(:);
        jj_diag = Faces(:);
        ww_diag = repmat(areas / 6, 3, 1);
        ii_off = [Faces(:,1); Faces(:,2); Faces(:,3); Faces(:,2); Faces(:,3); Faces(:,1)];
        jj_off = [Faces(:,2); Faces(:,3); Faces(:,1); Faces(:,1); Faces(:,2); Faces(:,3)];
        ww_off = repmat(areas / 12, 6, 1);
        M = sparse([ii_diag; ii_off], [jj_diag; jj_off], [ww_diag; ww_off], nV, nV);
        M = (M + M') / 2;

    otherwise
        error('Unknown MassType: %s. Use ''barycentric'', ''voronoi'', or ''galerkin''.', MassType);
end
end
```

- [ ] **Step 3: Replace the inline mass block with a call to the helper**

In `BST/toolbox/anatomy/tess_laplacian.m`, replace the entire `%% ===== MASS MATRIX =====` section (the `areas = 0.5 * dblArea;` line through the closing `end` of the `switch MassType`, i.e. the original mass `switch`) with:

```matlab
%% ===== MASS MATRIX =====
M = local_mass_matlab(Vertices, Faces, MassType);
```

> Leave the cotangent-Laplacian assembly and the `Symmetrize` handling for `L` untouched. The `galerkin` `Symmetrize` is now applied inside the helper (was `(M+M')/2`), matching prior behavior.

- [ ] **Step 4: Run the guard test to confirm no behavior change**

Run (MATLAB MCP `run_matlab_test_file`): `dev/tests/test_laplacian_ico.m`
Expected: PASS — `ALL TESTS PASSED: test_laplacian_ico`.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/anatomy/tess_laplacian.m
git commit -m "tess_laplacian: extract local_mass_matlab helper (no behavior change)"
```

---

## Task 5: Add nxr backend + `Backend` option to `tess_laplacian`

**Repo:** BST

**Files:**
- Modify: `BST/toolbox/anatomy/tess_laplacian.m`

- [ ] **Step 1: Add the `Backend` option to the input parser**

In the `%% ===== PARSE INPUTS =====` block, add a default and a case. Change:

```matlab
MassType = 'barycentric';
Symmetrize = true;
CheckManifold = false;
```
to:
```matlab
MassType = 'barycentric';
Symmetrize = true;
CheckManifold = false;
Backend = 'auto';          % 'auto' | 'nxr' | 'matlab'
```
and in the `switch lower(varargin{i})` block add:
```matlab
        case 'backend',       Backend = lower(varargin{i+1});
```

- [ ] **Step 2: Add the two local helpers**

At the end of `BST/toolbox/anatomy/tess_laplacian.m` (before the final `end`), add:

```matlab
function tf = nxr_is_loaded()
% True only if the nxr-compute plugin is currently loaded (cheap, in-memory).
tf = false;
try
    PlugDesc = bst_plugin('GetInstalled', 'nxr-compute');
    tf = ~isempty(PlugDesc) && isfield(PlugDesc, 'isLoaded') && PlugDesc.isLoaded;
catch
    tf = false;   % bst_plugin unavailable / any error -> MATLAB path
end
end

function backend = resolve_backend(req)
% Map the requested backend to an effective one ('nxr' or 'matlab').
switch req
    case 'matlab'
        backend = 'matlab';
    case 'nxr'
        if ~nxr_is_loaded()
            error('tess_laplacian:nxrNotLoaded', ...
                ['Backend ''nxr'' requested but the nxr-compute plugin is not loaded. ' ...
                 'Install/load it via bst_plugin(''Install'',''nxr-compute'').']);
        end
        backend = 'nxr';
    otherwise   % 'auto'
        if nxr_is_loaded()
            backend = 'nxr';
        else
            backend = 'matlab';
        end
end
end
```

- [ ] **Step 3: Insert the backend branch before the MATLAB assembly**

Immediately after the input-validation block (after `Faces = double(Faces);` and the optional `CheckManifold` block, i.e. right before `%% ===== COTANGENT LAPLACIAN =====`), insert:

```matlab
%% ===== BACKEND SELECTION =====
% tess_laplacian is a hot-loop operator: it uses nxr only if already loaded
% and never auto-installs here. Installation is a one-time action (Plugins
% menu, or a higher-level pipeline calling bst_plugin('Install','nxr-compute')).
backend = resolve_backend(Backend);
if strcmp(backend, 'nxr')
    try
        ctx = nxr.manifold.context(Vertices, Faces);   % faces passed 1-based (marshal subtracts 1)
        L = ctx.K;
        if strcmp(MassType, 'voronoi')
            M = ctx.M;                                  % nxr Voronoi mass
        else
            M = local_mass_matlab(Vertices, Faces, MassType);
        end
        if Symmetrize
            L = (L + L') / 2;
        end
        return;
    catch ME
        if strcmp(Backend, 'nxr')   % explicit request: surface the error
            rethrow(ME);
        end
        % 'auto': degrade gracefully to the MATLAB implementation
        warning('tess_laplacian:nxrFallback', ...
            'nxr backend failed (%s); falling back to MATLAB.', ME.message);
    end
end
% Fall through to the pure-MATLAB implementation below.
```

> The existing MATLAB cotangent assembly and `M = local_mass_matlab(...)` (from Task 4) remain as the fall-through path.

- [ ] **Step 4: Run the guard test (must still pass whether or not nxr is loaded)**

Run (MATLAB MCP `run_matlab_test_file`): `dev/tests/test_laplacian_ico.m`
Expected: PASS — `ALL TESTS PASSED: test_laplacian_ico` (it does not load nxr, so it exercises the MATLAB path through `'auto'`).

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/anatomy/tess_laplacian.m
git commit -m "tess_laplacian: nxr-first backend with MATLAB fallback (Backend option)"
```

---

## Task 6: Parity test — nxr vs MATLAB `L` and Voronoi `M` (the pass/fail gate)

**Repo:** BST · **Prereq:** staged plugin (Task 1).

**Files:**
- Create: `BST/dev/tests/test_laplacian_nxr_parity.m`

- [ ] **Step 1: Write the parity test**

Create `BST/dev/tests/test_laplacian_nxr_parity.m`:

```matlab
function test_laplacian_nxr_parity
% nxr vs MATLAB parity for tess_laplacian: stiffness L (all mass types) and
% Voronoi mass M. Stages + loads the locally-built nxr-compute plugin.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% Stage + load the plugin (same offline install as the lifecycle test)
codeDir  = fileparts(repoRoot);
nxrStage = fullfile(codeDir, 'nxr-compute', 'dist', 'plugin', 'nxr-compute');
assert(isfolder(nxrStage), 'Run nxr-compute/scripts/package-plugin.sh first (missing %s).', nxrStage);
destDir = fullfile(bst_get('UserPluginsDir'), 'nxr-compute');
if isfolder(destDir), file_delete(destDir, 1, 3); end
mkdir(destDir);
copyfile(fullfile(nxrStage, '*'), destDir);
[isOk, errMsg] = bst_plugin('Load', 'nxr-compute');
assert(isOk, 'Load failed: %s', errMsg);
cleanup = onCleanup(@() bst_plugin('Unload', 'nxr-compute'));

[V, F] = tess_sphere(642);

% --- Stiffness L parity, for every mass type (L is mass-independent) ---
tolL = 1e-6;
for mt = {'barycentric', 'voronoi', 'galerkin'}
    Ln = tess_laplacian(V, F, 'MassType', mt{1}, 'Backend', 'nxr');
    Lm = tess_laplacian(V, F, 'MassType', mt{1}, 'Backend', 'matlab');
    dL = full(max(abs(Ln(:) - Lm(:))));
    fprintf('L parity (%s): max|dL| = %.3e\n', mt{1}, dL);
    assert(dL < tolL, 'Stiffness L parity failed for %s (max|dL|=%.3e).', mt{1}, dL);
end

% --- Voronoi mass M parity (the only variant nxr serves) ---
% Two implementations may differ slightly in obtuse-triangle handling, so use
% a relative tolerance against the total mesh area.
[~, Mn] = tess_laplacian(V, F, 'MassType', 'voronoi', 'Backend', 'nxr');
[~, Mm] = tess_laplacian(V, F, 'MassType', 'voronoi', 'Backend', 'matlab');
totalArea = full(sum(diag(Mm)));
dM = full(max(abs(diag(Mn) - diag(Mm))));
relM = dM / (totalArea / size(V,1));   % relative to mean vertex area
fprintf('M parity (voronoi): max|dM| = %.3e, relative = %.3e\n', dM, relM);
assert(relM < 1e-3, 'Voronoi mass parity failed (relative=%.3e). Investigate obtuse handling.', relM);

% Sanity on the nxr stiffness itself
Ln = tess_laplacian(V, F, 'MassType', 'voronoi', 'Backend', 'nxr');
assert(norm(Ln - Ln', 1) < 1e-9, 'nxr L not symmetric.');
assert(max(abs(sum(Ln, 2))) < 1e-6, 'nxr L row sums not ~0.');
fprintf('ALL TESTS PASSED: test_laplacian_nxr_parity\n');
end
```

- [ ] **Step 2: Run it**

Run (MATLAB MCP `run_matlab_test_file`): `dev/tests/test_laplacian_nxr_parity.m`
Expected: PASS — per-variant `L parity` lines, an `M parity` line, and `ALL TESTS PASSED`.
If `L parity` fails by a large margin: suspect the faces convention or a sign flip — print `full([diag(Ln(1:3,1:3)), diag(Lm(1:3,1:3))])` to compare. If `M parity` fails only slightly: record the observed `relM` and, if the scheme genuinely differs, relax the tolerance with a comment; if large, restrict the nxr `M` path.

- [ ] **Step 3: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add dev/tests/test_laplacian_nxr_parity.m
git commit -m "tess_laplacian: nxr-vs-MATLAB parity test (L all variants, Voronoi M)"
```

---

## Task 7: Backend-selection + fallback test

**Repo:** BST · **Prereq:** staged plugin (Task 1).

**Files:**
- Create: `BST/dev/tests/test_laplacian_backend_select.m`

- [ ] **Step 1: Write the selection/fallback test**

Create `BST/dev/tests/test_laplacian_backend_select.m`:

```matlab
function test_laplacian_backend_select
% Verify backend selection: 'auto' tracks load state; 'nxr' errors when not
% loaded; 'matlab' always works; auto==nxr when loaded, auto==matlab when not.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

codeDir  = fileparts(repoRoot);
nxrStage = fullfile(codeDir, 'nxr-compute', 'dist', 'plugin', 'nxr-compute');
assert(isfolder(nxrStage), 'Run nxr-compute/scripts/package-plugin.sh first (missing %s).', nxrStage);
destDir = fullfile(bst_get('UserPluginsDir'), 'nxr-compute');
if isfolder(destDir), file_delete(destDir, 1, 3); end
mkdir(destDir);
copyfile(fullfile(nxrStage, '*'), destDir);

[V, F] = tess_sphere(162);

% --- With nxr UNLOADED ---
bst_plugin('Unload', 'nxr-compute');
La = tess_laplacian(V, F, 'Backend', 'auto');     % must use MATLAB
Lm = tess_laplacian(V, F, 'Backend', 'matlab');
assert(isequal(La, Lm), 'auto (unloaded) did not equal MATLAB backend.');
threw = false;
try
    tess_laplacian(V, F, 'Backend', 'nxr');        % must error when not loaded
catch ME
    threw = strcmp(ME.identifier, 'tess_laplacian:nxrNotLoaded');
end
assert(threw, 'Backend ''nxr'' did not raise nxrNotLoaded when unloaded.');
fprintf('PASSED: unloaded -> auto uses MATLAB; explicit nxr errors.\n');

% --- With nxr LOADED ---
[isOk, errMsg] = bst_plugin('Load', 'nxr-compute');
assert(isOk, 'Load failed: %s', errMsg);
cleanup = onCleanup(@() bst_plugin('Unload', 'nxr-compute'));
La2 = tess_laplacian(V, F, 'Backend', 'auto');
Ln  = tess_laplacian(V, F, 'Backend', 'nxr');
assert(isequal(La2, Ln), 'auto (loaded) did not equal nxr backend.');
fprintf('PASSED: loaded -> auto uses nxr.\n');

fprintf('ALL TESTS PASSED: test_laplacian_backend_select\n');
end
```

- [ ] **Step 2: Run it**

Run (MATLAB MCP `run_matlab_test_file`): `dev/tests/test_laplacian_backend_select.m`
Expected: PASS — two `PASSED` lines and `ALL TESTS PASSED`.

- [ ] **Step 3: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add dev/tests/test_laplacian_backend_select.m
git commit -m "tess_laplacian: backend selection + graceful fallback test"
```

---

## Task 8: Publish the macOS release and wire the real `URLzip` (OUTWARD-FACING — get user go-ahead)

**Repo:** NXR (release) + BST (URL) · **Requires the user's explicit confirmation** before publishing, since this pushes a public release asset.

**Files:**
- Modify: `BST/toolbox/core/bst_plugin.m` (confirm/adjust the `URLzip` + `Version`)

- [ ] **Step 1: Create the release and upload the zip (after confirmation)**

```bash
cd /Users/diellorbasha/workspace/research/code/nxr-compute
scripts/package-plugin.sh 1.0.0 mac
gh release create plugin-v1.0.0 \
  dist/nxr-compute-1.0.0-mac.zip \
  --repo neurodynamics-xr/nxr-compute \
  --title "Brainstorm plugin 1.0.0 (macOS)" \
  --notes "macOS (Apple Silicon) MEX + +nxr wrapper for the Brainstorm nxr-compute plugin."
```
Expected: `gh` prints the release URL; the asset `nxr-compute-1.0.0-mac.zip` is attached.

- [ ] **Step 2: Set the real `URLzip` and `Version` in the PlugDesc**

In `BST/toolbox/core/bst_plugin.m`, update the `mac64arm` arm to the published asset and bump `Version`:

```matlab
    PlugDesc(end).Version        = '1.0.0';
    ...
        case 'mac64arm'
            PlugDesc(end).URLzip   = 'https://github.com/neurodynamics-xr/nxr-compute/releases/download/plugin-v1.0.0/nxr-compute-1.0.0-mac.zip';
            PlugDesc(end).TestFile = 'nxr_compute.mexmaca64';
```

- [ ] **Step 3: Validate the real network install end-to-end**

Run (MATLAB MCP `evaluate_matlab_code`):
```matlab
% Remove any locally-staged copy so this exercises the real download
d = fullfile(bst_get('UserPluginsDir'), 'nxr-compute');
if isfolder(d), file_delete(d, 1, 3); end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute', 1);
assert(isOk, 'Network install failed: %s', errMsg);
v = nxr.manifold.context(tess_sphere(42)); %#ok<NASGU>
disp('NETWORK INSTALL OK');
```
Expected: downloads the asset, loads it, prints `NETWORK INSTALL OK`.

- [ ] **Step 4: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/core/bst_plugin.m
git commit -m "nxr-compute plugin: point URLzip at published macOS release 1.0.0"
```

---

## Task 9: Final integration — run the whole suite

**Repo:** BST

- [ ] **Step 1: Run all nxr-related tests in sequence**

Run (MATLAB MCP `run_matlab_test_file`), each expecting `ALL TESTS PASSED`:
- `dev/tests/test_nxr_plugin_registered.m`
- `dev/tests/test_nxr_plugin_lifecycle.m`
- `dev/tests/test_laplacian_ico.m`            (regression: existing behavior intact)
- `dev/tests/test_laplacian_nxr_parity.m`
- `dev/tests/test_laplacian_backend_select.m`

- [ ] **Step 2: Confirm the branch is clean and push (only if the user asks)**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git status
git log --oneline development..HEAD
```
Expected: working tree clean; commits from Tasks 2–8 on `feature/nxr-compute-plugin`. Push only on explicit request.

---

## Self-Review

**Spec coverage:**
- Deliverable 1 (package, mac) → Task 1 (+ Task 8 publish). ✓
- Deliverable 2 (`PlugDesc` entry) → Task 2. ✓
- Deliverable 3 (`tess_laplacian` nxr-first + fallback; L all variants, Voronoi M, others MATLAB) → Tasks 4–5; verified by Task 6. ✓
- Deliverable 4 (tests: lifecycle, parity gate, fallback) → Tasks 3, 6, 7; suite in Task 9. ✓
- Faces 0/1-based reconciliation → resolved (pass 1-based), guarded by Tasks 3 & 6. ✓
- Sign/symmetry/manifold reconciliation → Task 5 (`Symmetrize`, try/catch fallback) + Task 6 sanity asserts. ✓
- Multi-platform deferred → documented in spec; per-OS `switch` laid down in Task 2; CI not built. ✓
- `bst_nxr` bridge / other capabilities / GUI processes → explicitly out of scope; no tasks. ✓

**Placeholder scan:** No `TBD`/`TODO`/"handle edge cases". The `plugin-dev`/`1.0.0` URLs are concrete values set in Task 2 (dev) and finalized in Task 8 (release). ✓

**Type/name consistency:** `nxr_is_loaded`, `resolve_backend`, `local_mass_matlab`, and the `Backend` option are defined in Tasks 4–5 and used consistently in Tasks 6–7. `bst_plugin` command strings (`GetSupported`, `GetInstalled`, `Load`, `Unload`, `Install`) match the verified API. Test mesh helper `tess_sphere(N)` used consistently. ✓
