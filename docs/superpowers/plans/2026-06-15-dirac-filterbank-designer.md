# Spatial Filterbank Designer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a transient interactive designer that lets a user pick a spectral kernel, click a delta on the cortex (or load a source map), choose direction/chirality for the Dirac operator, preview the filtered field live in a reused `figure_3d`, tile the eigenmode spectrum into a small wavelet bank, and save it as a reusable `filterbank_` node nested under the `eigen_` node.

**Architecture:** A thin orchestrator (`view_filter_designer`) opens a `figure_3d` cortex preview and docks a transient control panel (`panel_filter_designer`) that holds session state and a live-preview controller; the controller reuses `bst_dirac_eigenmodes_filter` for the project→scale→chirality→reconstruct math and a new `bst_filterbank_tiles` generator for spectrum tiling. The saved artifact is a recipe bank (no materialized fields) stored as a top-level `Surface.Filterbank` DB list keyed by a `ParentEigen` reference and rendered as a child of the eigen node. Operator-general with Dirac wired first; LBO/connection are deferred render modes.

**Tech Stack:** MATLAB (Brainstorm toolbox), Brainstorm DB layer (`db_template`/`db_update`/`bst_get`/`db_add_*`), Java-Swing GUI via `gui_component`/`gui_show`/`gui_hide`, the `bst_eigfilter_*` kernel registry, and the `bst-java` fork (new tree node type). Tests run through the MATLAB MCP (`run_matlab_file`) against the live dev protocol (Subject01 has a Dirac `eigen_` node) and with small synthetic structs for pure-logic tasks.

**Repos / branches:**
- `brainstorm3` on `feat/filterbank-designer` (already created).
- `bst-java` fork: create a paired feature branch for the node-type task (Task 11). Never push/merge upstream.

**Conventions used by this plan:**
- Brainstorm panel functions dispatch subfunctions via `eval(macro_method)`, so `panel_filter_designer('FnName', args...)` calls the local `FnName` — this is how the live-preview controller is exercised headless in tests.
- Tests live in `dev/tests/`, are plain functions that `error()` on failure, and are run with the MATLAB MCP `run_matlab_file`. Brainstorm must be running for DB/GUI tests; pure-logic tests build synthetic inputs and need no protocol.
- Commit after each task with the message shown in its final step.

---

## File Structure

**Create (brainstorm3):**
- `toolbox/math/bst_filterbank_tiles.m` — pure tiling generator: one base design → N spectrum-spanning tiles (× chiralities).
- `toolbox/db/db_add_filterbank.m` — save a `filterbank_*.mat` and register it in `Surface.Filterbank` keyed by `ParentEigen`.
- `toolbox/gui/view_filter_designer.m` — session orchestrator (open figure + dock panel + teardown).
- `toolbox/gui/panel_filter_designer.m` — transient control panel: widgets, session state, live-preview controller.
- `dev/tests/test_filterbank_schema.m`, `dev/tests/test_filterbank_tiles.m`, `dev/tests/test_dirac_filter_coeffs.m`, `dev/tests/test_db_add_filterbank.m`, `dev/tests/test_filter_designer_session.m` — the test suite.

**Modify (brainstorm3):**
- `toolbox/db/db_template.m` — add `filterbankmat` template; add `Filterbank` list to the `surface` template.
- `toolbox/db/db_update.m` — bump `CurrentDbVersion` to 5.05; add a `< 5.05` migration block (reuses `NormalizeSurfaceArray`).
- `toolbox/core/bst_startup.m` — bump the `CurrentDbVersion` constant to 5.05.
- `toolbox/core/bst_get.m` — add a `FilterbankFile` accessor case.
- `toolbox/math/bst_dirac_eigenmodes_filter.m` — add `ReturnCoeffs`/`Coeffs` options (projection/filter split).
- `toolbox/gui/view_eigfilter_response.m` — accept a tile bank (cell of kernels) + active index; make curves clickable.
- `toolbox/tree/node_create_subject.m` — nest `filterbank` nodes under their parent eigen node.
- `toolbox/tree/tree_callbacks.m` — "Design filterbank…" on the eigen node; filterbank-node popup (Open/Delete); cascade-delete filterbanks when an eigen node is deleted.

**Modify (bst-java fork):**
- The `org.brainstorm.tree.BstNode` type table + icon for the `'filterbank'` node type.

---

## Task 1: Filterbank schema in db_template

**Files:**
- Modify: `toolbox/db/db_template.m` (the `'eigenmat'`/`'operatormat'` neighborhood ~line 131; the `'surface'` case ~line 32)
- Test: `dev/tests/test_filterbank_schema.m`

- [ ] **Step 1: Write the failing test**

```matlab
function test_filterbank_schema()
% Schema regression for the filterbankmat template and the Surface.Filterbank list.
% Authors: Diellor Basha, 2026
    nFail = 0;
    fb = db_template('filterbankmat');
    need = {'Comment','ParentEigen','Variant','Tiles','Tiling','Provenance'};
    for f = need
        if ~isfield(fb, f{1}); fprintf('MISSING filterbankmat.%s\n', f{1}); nFail = nFail+1; end
    end
    s = db_template('surface');
    if ~isfield(s, 'Filterbank'); fprintf('MISSING surface.Filterbank\n'); nFail = nFail+1; end
    if isfield(s,'Filterbank')
        sub = fieldnames(s.Filterbank);
        for f = {'FileName','Comment','ParentEigen'}
            if ~ismember(f{1}, sub); fprintf('MISSING surface.Filterbank.%s\n', f{1}); nFail = nFail+1; end
        end
        if ~isempty(s.Filterbank); fprintf('surface.Filterbank must start EMPTY (0x0 struct)\n'); nFail = nFail+1; end
    end
    fprintf('\n==== test_filterbank_schema: %d failed ====\n', nFail);
    if nFail > 0, error('test_filterbank_schema FAILED'); end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (MCP `run_matlab_file`): `dev/tests/test_filterbank_schema.m`
Expected: FAIL — `MISSING filterbankmat.*` (the `'filterbankmat'` case does not exist yet).

- [ ] **Step 3: Add the `filterbankmat` template**

In `toolbox/db/db_template.m`, immediately after the `case 'eigenmat'` block (after its closing `);`), add:

```matlab
    case 'filterbankmat'
        template = struct(...
              'Comment',     '', ...
              'ParentEigen', '', ...   % file_short of the eigen_ node this bank applies in
              'Variant',     '', ...   % operator variant inherited from eigen (e.g. 'Dirac')
              'Tiles',       [], ...   % 1xN struct array (Kernel,Params,Direction,Chirality,Axis)
              'Tiling',      [], ...   % generator meta (Kernel,N,Spacing,LambdaRange,Chiralities)
              'Provenance',  []);
```

- [ ] **Step 4: Add the `Filterbank` list to the `surface` template**

In the `case 'surface'` block, extend the struct so it reads:

```matlab
    case 'surface'
        template = struct('Comment',     '', ...
                          'FileName',    '', ...
                          'SurfaceType', '', ...
                          'Manifold',    struct('FileName',{},'Comment',{}), ...
                          'Eigen',       struct('FileName',{},'Comment',{},'Variant',{}), ...
                          'Operator',    struct('FileName',{},'Comment',{},'Variant',{}), ...
                          'Filterbank',  struct('FileName',{},'Comment',{},'ParentEigen',{}));
```

- [ ] **Step 5: Run test to verify it passes**

Run: `dev/tests/test_filterbank_schema.m`
Expected: PASS — `0 failed`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/db/db_template.m dev/tests/test_filterbank_schema.m
git commit -m "feat(db): filterbankmat template + Surface.Filterbank list"
```

---

## Task 2: DB version bump + migration (5.05)

**Files:**
- Modify: `toolbox/core/bst_startup.m:208` (the `CurrentDbVersion = 5.04;` constant)
- Modify: `toolbox/db/db_update.m` (after the `if (CurrentDbVersion < 5.04)` block, ~line 403)

This backfills the empty `Filterbank` field on existing surfaces in older protocols, reusing the same `NormalizeSurfaceArray` path the 5.04 manifold/eigen migration used.

- [ ] **Step 1: Bump the version constant**

In `toolbox/core/bst_startup.m` change:

```matlab
CurrentDbVersion = 5.04;
```
to
```matlab
CurrentDbVersion = 5.05;
```

- [ ] **Step 2: Add the migration block**

In `toolbox/db/db_update.m`, immediately after the closing `end` of the `if (CurrentDbVersion < 5.04)` block (before the `%% ===== JUST BEFORE RETURNING` comment), add:

```matlab
if (CurrentDbVersion < 5.05)
    isSurfFixed = 0;
    templateSurface = db_template('Surface');
    for iProt = 1:length(ProtocolsListSubjects)
        subjFields = fieldnames(ProtocolsListSubjects(iProt));
        for iField = 1:length(subjFields)
            subjField = subjFields{iField};
            for iSubj = 1:length(ProtocolsListSubjects(iProt).(subjField))
                sSurf = ProtocolsListSubjects(iProt).(subjField)(iSubj).Surface;
                if isempty(sSurf)
                    continue;
                end
                [sSurf, isFix] = NormalizeSurfaceArray(sSurf, templateSurface);
                if isFix
                    ProtocolsListSubjects(iProt).(subjField)(iSubj).Surface = sSurf;
                    isSurfFixed = 1;
                end
            end
        end
    end
    if isSurfFixed
        disp('BST> Database structure: Adding filterbank support to surfaces...');
        SaveProtocolSubjects();
        disp('BST> Database structure: Done.');
    end
end
```

- [ ] **Step 3: Verify the migration runs on the live protocol**

Run (MCP `evaluate_matlab_code`):

```matlab
db_update(5.05);
sSubject = bst_get('Subject', 1);
assert(isfield(sSubject.Surface, 'Filterbank'), 'Filterbank field not backfilled');
fprintf('OK: Surface.Filterbank present after migration\n');
```

Expected: prints `OK: Surface.Filterbank present after migration` (and, the first time, the "Adding filterbank support" line).

- [ ] **Step 4: Commit**

```bash
git add toolbox/core/bst_startup.m toolbox/db/db_update.m
git commit -m "feat(db): DB v5.05 migration adds Surface.Filterbank"
```

---

## Task 3: bst_get FilterbankFile accessor

**Files:**
- Modify: `toolbox/core/bst_get.m` (add a `case 'FilterbankFile'` modeled on the existing `case 'EigenFile'`)
- Test: covered by Task 5's round-trip test (no standalone test — the accessor is exercised there)

- [ ] **Step 1: Find the EigenFile accessor to mirror**

Run (MCP `evaluate_matlab_code`):

```matlab
file = which('bst_get');
txt = fileread(file);
k = strfind(txt, "case 'EigenFile'");
disp(txt(k:k+1200));
```

Expected: prints the `case 'EigenFile'` block — note how it loops `sSubject.Surface(s).Eigen` and returns `[sSubject, iSubject, iSurface, iItem]`.

- [ ] **Step 2: Add the FilterbankFile accessor**

Immediately after the `case 'EigenFile'` block in `toolbox/core/bst_get.m`, add the analogous block (replace `Eigen` with `Filterbank` throughout; same return signature `[argout1=sSubject, argout2=iSubject, argout3=iSurface, argout4=iFilterbank]`):

```matlab
    case 'FilterbankFile'
        % Mirror of 'EigenFile': resolve a filterbank_*.mat to its subject/surface/index.
        FileName = file_short(varargin{1});
        % Look in all the surfaces of all the subjects
        ProtocolSubjects = GlobalData.DataBase.ProtocolSubjects(GlobalData.DataBase.iProtocol);
        for iList = 0:length(ProtocolSubjects.Subject)
            if (iList == 0); sSubject = ProtocolSubjects.DefaultSubject; iSubjOut = 0;
            else;            sSubject = ProtocolSubjects.Subject(iList);  iSubjOut = iList; end
            if isempty(sSubject) || ~isfield(sSubject,'Surface'); continue; end
            for iSurf = 1:length(sSubject.Surface)
                if ~isfield(sSubject.Surface(iSurf),'Filterbank'); continue; end
                fbs = sSubject.Surface(iSurf).Filterbank;
                for iFb = 1:length(fbs)
                    if file_compare(fbs(iFb).FileName, FileName)
                        argout1 = sSubject; argout2 = iSubjOut; argout3 = iSurf; argout4 = iFb;
                        return;
                    end
                end
            end
        end
        argout1 = []; argout2 = []; argout3 = []; argout4 = [];
```

> Note: if the live `EigenFile` block differs in variable names (e.g. uses `argout1..4` vs explicit returns), copy its exact shape and only swap `Eigen`→`Filterbank`. The block above is the canonical shape; reconcile with Step 1's printout before saving.

- [ ] **Step 3: Smoke-check it resolves nothing gracefully**

Run (MCP `evaluate_matlab_code`):

```matlab
[s,iSub,iSurf,iFb] = bst_get('FilterbankFile', 'does/not/exist.mat');
assert(isempty(s) && isempty(iFb), 'expected empty resolution');
fprintf('OK: FilterbankFile returns empty for unknown file\n');
```

Expected: `OK: FilterbankFile returns empty for unknown file`.

- [ ] **Step 4: Commit**

```bash
git add toolbox/core/bst_get.m
git commit -m "feat(db): bst_get FilterbankFile accessor"
```

---

## Task 4: Tiling generator (bst_filterbank_tiles)

**Files:**
- Create: `toolbox/math/bst_filterbank_tiles.m`
- Test: `dev/tests/test_filterbank_tiles.m`

The generator turns one base design into N tiles whose scale parameter spans the eigenvalue range geometrically (constant ratio between adjacent tiles), optionally crossed with chiralities. It is pure (no DB, no GUI).

- [ ] **Step 1: Write the failing test**

```matlab
function test_filterbank_tiles()
% Pure-logic test for the spectrum tiling generator.
% Authors: Diellor Basha, 2026
    nFail = 0;
    base = struct('Kernel','mexhat', 'Params',struct(), ...
                  'Direction',[1 0 0], 'Chirality',0, 'Axis',[0 0 1], ...
                  'N',4, 'Spacing','geometric', 'LambdaRange',[1 256], 'Chiralities',0);

    T = bst_filterbank_tiles(base);
    nFail = nFail + chk('4 tiles, no chirality split', numel(T)==4);
    nFail = nFail + chk('each tile carries Kernel', all(arrayfun(@(t) strcmp(t.Kernel,'mexhat'), T)));

    % mexhat peaks at l = 1/t, so the per-tile t must span the range geometrically:
    % t_j = 1/center_j, centers geometric from 1 to 256 -> ratio 4 across 4 tiles.
    centers = arrayfun(@(t) 1./t.Params.t, T);
    ratios  = centers(2:end) ./ centers(1:end-1);
    nFail = nFail + chk('geometric centers (constant ratio)', max(abs(ratios - ratios(1))) < 1e-9);
    nFail = nFail + chk('span covers LambdaRange ends', abs(centers(1)-1)<1e-6 && abs(centers(end)-256)<1e-6);

    % chirality split doubles the bank into +1/-1 with matched scales
    base2 = base; base2.Chiralities = [1 -1];
    T2 = bst_filterbank_tiles(base2);
    nFail = nFail + chk('chirality split doubles count', numel(T2)==8);
    nFail = nFail + chk('both signs present', any([T2.Chirality]==1) && any([T2.Chirality]==-1));

    fprintf('\n==== test_filterbank_tiles: %d failed ====\n', nFail);
    if nFail > 0, error('test_filterbank_tiles FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `dev/tests/test_filterbank_tiles.m`
Expected: FAIL — `Undefined function 'bst_filterbank_tiles'`.

- [ ] **Step 3: Implement the generator**

Create `toolbox/math/bst_filterbank_tiles.m`:

```matlab
function Tiles = bst_filterbank_tiles(base)
% BST_FILTERBANK_TILES: Expand one base filter design into a spectrum-spanning bank.
%
% USAGE:  Tiles = bst_filterbank_tiles(base)
%
% INPUT base (struct):
%   .Kernel      kernel name in the bst_eigfilter registry (e.g. 'mexhat','heat')
%   .Params      struct of fixed kernel params (the scale param is overwritten per tile)
%   .Direction   [1x3] launch direction (Dirac); copied to every tile
%   .Chirality   scalar 0 | +1 | -1 (the base sign; ignored if .Chiralities set)
%   .Axis        [1x3] chirality axis (Dirac); copied to every tile
%   .N           number of spectral tiles (>=1)
%   .Spacing     'geometric' (default) | 'linear' — spacing of tile centers in lambda
%   .LambdaRange [lo hi] eigenvalue band the tiles tile (lo>0 for geometric)
%   .Chiralities [] or 0 => single bank at .Chirality; a vector (e.g. [1 -1]) crosses
%                every spectral tile with each listed sign (doubling/ tripling the bank)
%
% OUTPUT:
%   Tiles : 1xM struct array, each (Kernel,Params,Direction,Chirality,Axis), where
%           M = N * max(1,numel(.Chiralities)). The scale parameter of .Params is set
%           per tile so the tile CENTERS (the lambda of peak response) span LambdaRange.
%
% The scale param name and the center<->param map are read from the kernel's 'meta'
% (bandpass kernels like mexhat/dog peak at lambda=1/t; low/high-pass use the param
% directly as a cutoff). Unknown kernels fall back to a linear sweep of the first param.
%
% SEE ALSO: bst_eigfilter_kernel, bst_dirac_eigenmodes_filter, panel_filter_designer
%
% Authors: Diellor Basha, 2026

    if ~isfield(base,'N') || isempty(base.N) || base.N < 1, base.N = 1; end
    if ~isfield(base,'Spacing') || isempty(base.Spacing), base.Spacing = 'geometric'; end
    if ~isfield(base,'LambdaRange') || numel(base.LambdaRange) ~= 2, base.LambdaRange = [1 100]; end
    if ~isfield(base,'Params') || ~isstruct(base.Params), base.Params = struct(); end
    if ~isfield(base,'Direction') || isempty(base.Direction), base.Direction = [1 0 0]; end
    if ~isfield(base,'Axis') || isempty(base.Axis), base.Axis = [0 0 1]; end
    if ~isfield(base,'Chirality') || isempty(base.Chirality), base.Chirality = 0; end
    if ~isfield(base,'Chiralities'), base.Chiralities = []; end

    lo = base.LambdaRange(1); hi = base.LambdaRange(2);
    N  = base.N;
    if strcmpi(base.Spacing,'geometric')
        if lo <= 0, lo = max(eps, hi*1e-3); end
        if N == 1, centers = sqrt(lo*hi); else, centers = exp(linspace(log(lo), log(hi), N)); end
    else
        if N == 1, centers = (lo+hi)/2; else, centers = linspace(lo, hi, N); end
    end

    [scaleName, centerToParam] = i_scale_map(base.Kernel);

    signs = base.Chiralities;
    if isempty(signs), signs = base.Chirality; end

    Tiles = repmat(i_blank_tile(), 1, N*numel(signs));
    t = 0;
    for s = 1:numel(signs)
        for j = 1:N
            t = t + 1;
            p = base.Params;
            p.(scaleName) = centerToParam(centers(j));
            Tiles(t) = struct('Kernel', base.Kernel, 'Params', p, ...
                'Direction', base.Direction(:).', 'Chirality', signs(s), 'Axis', base.Axis(:).');
        end
    end
end

function t = i_blank_tile()
    t = struct('Kernel','', 'Params',struct(), 'Direction',[1 0 0], 'Chirality',0, 'Axis',[0 0 1]);
end

function [scaleName, centerToParam] = i_scale_map(kernelName)
% Map a desired peak/cutoff lambda (center) to the kernel's scale parameter.
    switch lower(kernelName)
        case {'mexhat','dog'}                 % peak at lambda = 1/t  => t = 1/center
            scaleName = 't';   centerToParam = @(c) 1./max(c, eps);
        case {'heat','inverse_heat'}          % e^{-t*lambda}: use t = 1/center as a soft cutoff
            scaleName = 't';   centerToParam = @(c) 1./max(c, eps);
        case {'tikhonov'}                     % beta acts at lambda ~ beta
            scaleName = 'beta'; centerToParam = @(c) c;
        otherwise                             % fallback: first meta param swept linearly as the center
            meta = bst_eigfilter_kernel('info', kernelName);
            fn = fieldnames(meta.params);
            scaleName = fn{1}; centerToParam = @(c) c;
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `dev/tests/test_filterbank_tiles.m`
Expected: PASS — all `PASS`, `0 failed`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/math/bst_filterbank_tiles.m dev/tests/test_filterbank_tiles.m
git commit -m "feat(math): bst_filterbank_tiles spectrum tiling generator"
```

---

## Task 5: Projection/filter split in bst_dirac_eigenmodes_filter

**Files:**
- Modify: `toolbox/math/bst_dirac_eigenmodes_filter.m`
- Test: `dev/tests/test_dirac_filter_coeffs.m`

Add two options so the panel can project a seed once (`ReturnCoeffs`) and then re-apply different gains cheaply (`Coeffs`), while the all-in-one call stays the single source of truth. Equivalence: project-then-apply must equal the direct call.

- [ ] **Step 1: Write the failing test (synthetic, no DB)**

```matlab
function test_dirac_filter_coeffs()
% Projection/filter split equivalence on a tiny synthetic Dirac eigenbasis.
% Authors: Diellor Basha, 2026
    nFail = 0;
    rng(0);
    nV = 6; K = 8;
    % one-hemisphere synthetic basis: random B-orthonormal Phi wrt B = I (mass = I4 per vertex)
    B = speye(4*nV);
    [Q,~] = qr(randn(4*nV, K), 0);          % 4nV x K, orthonormal columns (B=I)
    Eigen = struct();
    Eigen.Phi = {Q, zeros(0,K)};            % hemi 2 empty
    Eigen.Lambda = {sort(rand(K,1)*100), zeros(0,1)};
    Eigen.GlobalVertices = {(1:nV).', zeros(0,1)};
    Mass = {B, sparse(0,0)};
    J = randn(3*nV, 1);

    % direct vs split (heat)
    Jdirect = bst_dirac_eigenmodes_filter(Eigen, Mass, J, 'heat', 'DiffusionTime', 0.02);
    [~, ~, c] = bst_dirac_eigenmodes_filter(Eigen, Mass, J, 'flat', 'ReturnCoeffs', true);
    Jsplit  = bst_dirac_eigenmodes_filter(Eigen, Mass, [], 'heat', 'DiffusionTime', 0.02, 'Coeffs', c);
    nFail = nFail + chk('split == direct (heat)', max(abs(Jdirect(:)-Jsplit(:))) < 1e-10);

    % re-applying a different gain to the SAME cached coeffs matches a fresh direct call
    Jd2 = bst_dirac_eigenmodes_filter(Eigen, Mass, J, 'bandpass', 'ModeRange', [2 6]);
    Js2 = bst_dirac_eigenmodes_filter(Eigen, Mass, [], 'bandpass', 'ModeRange', [2 6], 'Coeffs', c);
    nFail = nFail + chk('cached coeffs reusable across gains', max(abs(Jd2(:)-Js2(:))) < 1e-10);

    fprintf('\n==== test_dirac_filter_coeffs: %d failed ====\n', nFail);
    if nFail > 0, error('test_dirac_filter_coeffs FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `dev/tests/test_dirac_filter_coeffs.m`
Expected: FAIL — the 3rd output `c` is not returned / `Coeffs` option unknown (the function errors or returns wrong shape).

- [ ] **Step 3: Implement the split**

In `toolbox/math/bst_dirac_eigenmodes_filter.m`:

(a) Change the signature line to expose a third output:

```matlab
function [JFilt, h, Coeffs] = bst_dirac_eigenmodes_filter(EigenMat, MassCell, J, FilterType, varargin)
```

(b) In the option-parse loop that currently extracts `Chirality`, also extract `ReturnCoeffs` and `Coeffs` (drop them from `gainArgs`):

```matlab
    Chirality = [];  ReturnCoeffs = false;  CoeffsIn = [];
    keep = true(1, numel(varargin));
    for i = 1:2:numel(varargin)
        if ischar(varargin{i})
            switch lower(varargin{i})
                case 'chirality',    Chirality    = varargin{i+1}; keep(i:i+1) = false;
                case 'returncoeffs', ReturnCoeffs = logical(varargin{i+1}); keep(i:i+1) = false;
                case 'coeffs',       CoeffsIn     = varargin{i+1}; keep(i:i+1) = false;
            end
        end
    end
    gainArgs = varargin(keep);
    Coeffs = cell(1, numel(EigenMat.Phi));
```

(c) Inside the per-hemisphere loop, replace the projection `c = Phi' * (B * Psi);` so it uses cached coeffs when provided, and stores the computed coeffs:

```matlab
        if ~isempty(CoeffsIn)
            c = CoeffsIn{hh};                              % use cached projection (skip embed+project)
        else
            Jx = J(3*(gv-1)+1, :);  Jy = J(3*(gv-1)+2, :);  Jz = J(3*(gv-1)+3, :);
            Psi = zeros(4*nVh, size(J,2));
            Psi(2:4:end, :) = Jx;  Psi(3:4:end, :) = Jy;  Psi(4:4:end, :) = Jz;
            c = Phi' * (B * Psi);                          % [K x nTime]
        end
        Coeffs{hh} = c;
```

> Note: when `CoeffsIn` is provided, `J` may be `[]`; derive `nTime = size(c,2)` and the output width from the coeffs instead of `J`. Update the `nTime`/`JFilt` preallocation near the top of the loop accordingly: `nTime = isempty(CoeffsIn) * size(J,2) + ~isempty(CoeffsIn) * size(CoeffsIn{1},2);` (or compute from whichever input is non-empty before the loop).

(d) After the loop, if `ReturnCoeffs` and the caller used `'flat'` purely to harvest coeffs, `JFilt` is still valid (flat gain = identity); no extra handling needed.

- [ ] **Step 4: Run test to verify it passes**

Run: `dev/tests/test_dirac_filter_coeffs.m`
Expected: PASS — both `PASS`, `0 failed`.

- [ ] **Step 5: Run the existing Dirac filter callers' smoke to catch regressions**

Run (MCP `evaluate_matlab_code`):

```matlab
% the 2-output and 1-output call shapes must still work unchanged
EigenMat = tess_eigen(bst_get('Subject',1).Surface(5).FileName, 'Dirac');
Op = load(file_fullpath(EigenMat.OperatorFile));
nV = double(max(cellfun(@(x) max(x(:)), EigenMat.GlobalVertices)));
J = randn(3*nV,1);
[Jf, h] = bst_dirac_eigenmodes_filter(EigenMat, Op.Mass, J, 'lowpass', 'CutoffMode', 50);
fprintf('OK: legacy 2-output call still works, size(Jf)=[%s]\n', num2str(size(Jf)));
```

Expected: prints `OK: legacy 2-output call still works ...` with `size(Jf)=[ <3nV> 1]`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/math/bst_dirac_eigenmodes_filter.m dev/tests/test_dirac_filter_coeffs.m
git commit -m "feat(math): projection/filter split (ReturnCoeffs/Coeffs) in dirac filter"
```

---

## Task 6: db_add_filterbank (save + register)

**Files:**
- Create: `toolbox/db/db_add_filterbank.m`
- Test: `dev/tests/test_db_add_filterbank.m`

Mirrors `db_add_eigen`, but resolves the surface from the parent eigen node's `ParentSurface` and appends to `Surface.Filterbank` with `ParentEigen` set.

- [ ] **Step 1: Write the failing round-trip test (live protocol)**

```matlab
function test_db_add_filterbank()
% Round-trip: save a filterbank under an eigen node, resolve it, reload it, delete it.
% Requires Brainstorm running with a Dirac eigen node (Subject01 surface 5).
% Authors: Diellor Basha, 2026
    nFail = 0;
    sSubject = bst_get('Subject', 1);
    iSurf = 5;
    assert(isfield(sSubject.Surface(iSurf),'Eigen') && ~isempty(sSubject.Surface(iSurf).Eigen), ...
        'No eigen node on Subject01 surface 5; compute one first.');
    EigenFile = sSubject.Surface(iSurf).Eigen(1).FileName;

    base = struct('Kernel','mexhat','Params',struct(),'Direction',[1 0 0], ...
                  'Chirality',0,'Axis',[0 0 1],'N',4,'Spacing','geometric', ...
                  'LambdaRange',[1 256],'Chiralities',0);
    fb = db_template('filterbankmat');
    fb.ParentEigen = file_short(EigenFile);
    fb.Variant     = 'Dirac';
    fb.Tiles       = bst_filterbank_tiles(base);
    fb.Tiling      = base;

    iFb = db_add_filterbank(1, EigenFile, fb, 'TEST filterbank');
    nFail = nFail + chk('returns an index', ~isempty(iFb));

    % registered under the same surface, keyed by ParentEigen
    sSubject = bst_get('Subject', 1);
    fbs = sSubject.Surface(iSurf).Filterbank;
    nFail = nFail + chk('appears in Surface.Filterbank', ~isempty(fbs));
    nFail = nFail + chk('ParentEigen matches', any(strcmp({fbs.ParentEigen}, file_short(EigenFile))));

    % resolvable via accessor + reloads to an identical recipe bank
    newFile = fbs(end).FileName;
    [s2,iSub2,iSurf2,iFb2] = bst_get('FilterbankFile', newFile);
    nFail = nFail + chk('accessor resolves it', ~isempty(iFb2) && iSurf2==iSurf);
    R = load(file_fullpath(newFile));
    nFail = nFail + chk('reloaded Tiles count', numel(R.Tiles)==numel(fb.Tiles));
    nFail = nFail + chk('reloaded Variant', strcmp(R.Variant,'Dirac'));

    % cleanup: delete the test node + file (proper DB path)
    file_delete(file_fullpath(newFile), 1);
    sSubject.Surface(iSurf).Filterbank(end) = [];
    ProtocolSubjects = bst_get('ProtocolSubjects');
    ProtocolSubjects.Subject(1) = sSubject;
    bst_set('ProtocolSubjects', ProtocolSubjects);
    db_save();

    fprintf('\n==== test_db_add_filterbank: %d failed ====\n', nFail);
    if nFail > 0, error('test_db_add_filterbank FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `dev/tests/test_db_add_filterbank.m`
Expected: FAIL — `Undefined function 'db_add_filterbank'`.

- [ ] **Step 3: Implement db_add_filterbank**

Create `toolbox/db/db_add_filterbank.m` (modeled on `db_add_eigen`):

```matlab
function iFilterbank = db_add_filterbank(iSubject, ParentEigenFile, FilterbankMat, Comment)
% DB_ADD_FILTERBANK: Save a filterbank_*.mat and register it as a child of an eigen node.
%
% USAGE:  iFilterbank = db_add_filterbank(iSubject, ParentEigenFile, FilterbankMat, Comment)
%
% The node is stored in the SURFACE's Filterbank list, keyed by ParentEigen (the eigen
% node file). node_create_subject nests it under the matching eigen node in the tree.
%
% Authors: Diellor Basha, 2026

    if (nargin < 4) || isempty(Comment); Comment = 'Filterbank'; end

    % Resolve the parent eigen node -> its subject/surface
    [sSubject, iSubjectE, iSurface, iEigen] = bst_get('EigenFile', ParentEigenFile);
    if isempty(iSurface)
        error('db_add_filterbank:noEigen', 'Parent eigen node not found: %s', ParentEigenFile);
    end
    if (nargin >= 1) && ~isempty(iSubject) && (iSubject ~= iSubjectE)
        % trust the resolved subject; keep the argument for API symmetry with db_add_eigen
    end
    iSubject = iSubjectE;

    ProtocolSubjects = bst_get('ProtocolSubjects');
    if (iSubject == 0); sSubject = ProtocolSubjects.DefaultSubject;
    else;               sSubject = ProtocolSubjects.Subject(iSubject); end

    % Build a unique filename in the parent surface's anatomy folder
    ProtocolInfo = bst_get('ProtocolInfo');
    c = clock;
    strTime = sprintf('%02.0f%02.0f%02.0f_%02.0f%02.0f', c(1)-2000, c(2:5));
    OutputFile = ['filterbank_' strTime '.mat'];
    OutputFileFull = file_unique(bst_fullfile(ProtocolInfo.SUBJECTS, bst_fileparts(sSubject.FileName), OutputFile));

    % Stamp required fields
    FilterbankMat.ParentEigen = file_short(ParentEigenFile);
    FilterbankMat.Comment     = Comment;

    bst_save(OutputFileFull, FilterbankMat, 'v7');

    % Normalize the surface array to the current template (adds Filterbank if missing),
    % mirroring db_add_eigen's homogeneity fix.
    templateSurface = db_template('Surface');
    tFields = fieldnames(templateSurface);
    if ~isempty(sSubject.Surface) && ~all(isfield(sSubject.Surface, tFields))
        sNew = cell(1, numel(sSubject.Surface));
        for k = 1:numel(sSubject.Surface)
            s = struct_copy_fields(sSubject.Surface(k), templateSurface, 0);
            extraFields = setdiff(fieldnames(s), tFields, 'stable');
            sNew{k} = orderfields(s, [tFields; extraFields]);
        end
        sSubject.Surface = reshape([sNew{:}], size(sSubject.Surface));
    end

    % Append the child entry
    newEntry.FileName    = file_short(OutputFileFull);
    newEntry.Comment     = Comment;
    newEntry.ParentEigen = file_short(ParentEigenFile);
    sSubject.Surface(iSurface).Filterbank(end+1) = newEntry;
    iFilterbank = numel(sSubject.Surface(iSurface).Filterbank);

    if (iSubject == 0); ProtocolSubjects.DefaultSubject = sSubject;
    else;               ProtocolSubjects.Subject(iSubject) = sSubject; end
    bst_set('ProtocolSubjects', ProtocolSubjects);

    panel_protocols('UpdateNode', 'Subject', iSubject);
    db_save();
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `dev/tests/test_db_add_filterbank.m`
Expected: PASS — all `PASS`, `0 failed`, and the test cleans up its node.

- [ ] **Step 5: Commit**

```bash
git add toolbox/db/db_add_filterbank.m dev/tests/test_db_add_filterbank.m
git commit -m "feat(db): db_add_filterbank saves + registers a filterbank node"
```

---

## Task 7: Multi-tile clickable spectrum strip (view_eigfilter_response)

**Files:**
- Modify: `toolbox/gui/view_eigfilter_response.m`
- Test: extends the session test (Task 10); add a focused smoke here.

Extend the existing single-curve plot to draw a *bank* of kernels with the active one highlighted, and call back into a handler on curve-click.

- [ ] **Step 1: Add a bank-drawing call signature**

Extend `view_eigfilter_response` to accept either a single handle `g` (today's behavior, unchanged) or a struct describing a bank:

```matlab
%   view_eigfilter_response(BANK, lambdas, titleStr)
%   BANK = struct('Kernels',{cellOfHandles}, 'Active', idx, 'OnSelect', @(idx)...)
```

In the function body, after the existing `if ... strcmpi(g,'close')` guard, branch:

```matlab
isBank = isstruct(g) && isfield(g,'Kernels');
```

When `isBank`, draw each kernel curve in grey, the active one in a thick colored line, store per-line `ButtonDownFcn` that calls `g.OnSelect(j)`:

```matlab
if isBank
    cla(ax); hold(ax,'on');
    xg = linspace(0, max([lambdas(:); eps]), 400)';
    for j = 1:numel(g.Kernels)
        yj = g.Kernels{j}(xg);
        isAct = (j == g.Active);
        lw = 1 + 1.5*isAct;
        col = [.6 .6 .6]; if isAct, col = [.85 .2 .2]; end
        hL = plot(ax, xg, yj, '-', 'Color', col, 'LineWidth', lw);
        set(hL, 'ButtonDownFcn', @(h,e) g.OnSelect(j));
    end
    title(ax, titleStr); 
    return;
end
```

(The single-handle path below remains exactly as today.)

- [ ] **Step 2: Smoke test (headless figure)**

Run (MCP `evaluate_matlab_code`):

```matlab
sel = 0;
bank = struct('Kernels', {{@(l)exp(-0.1*l), @(l)exp(-0.01*l), @(l)(0.05*l).*exp(-0.05*l)}}, ...
              'Active', 2, 'OnSelect', @(j) assignin('base','sel',j));
hFig = view_eigfilter_response(bank, linspace(0,100,20), 'bank smoke');
ax = getappdata(hFig,'Axes');
nLines = numel(findobj(ax,'Type','line'));
assert(nLines >= 3, 'expected >=3 curves');
fprintf('OK: bank strip drew %d curves\n', nLines);
view_eigfilter_response('close');
```

Expected: `OK: bank strip drew 3 curves` (or more, counting eigenvalue ticks).

- [ ] **Step 3: Commit**

```bash
git add toolbox/gui/view_eigfilter_response.m
git commit -m "feat(gui): multi-tile clickable spectrum strip in view_eigfilter_response"
```

---

## Task 8: Control panel — panel_filter_designer

**Files:**
- Create: `toolbox/gui/panel_filter_designer.m`
- Test: driven headlessly in Task 10's session test via `panel_filter_designer('BuildDesign', state)` and `panel_filter_designer('ComputeField', ...)`.

This is the largest unit. It owns: widgets (kernel dropdown + auto-params + direction/chirality + tiling + Save/Cancel), session state (the eigen basis, mass, the seed input, the active tile), and the live-preview controller (build tiles → compute the selected tile's field → push to the figure). Direction/chirality controls show only for `Variant=='Dirac'`.

Keep the controller logic in dispatchable subfunctions so it is testable without showing the panel (`eval(macro_method)` exposes them).

- [ ] **Step 1: Create the panel skeleton + state struct**

Create `toolbox/gui/panel_filter_designer.m`:

```matlab
function varargout = panel_filter_designer(varargin)
% PANEL_FILTER_DESIGNER: Transient control panel for the spatial filterbank designer.
% Owns the design state and the live-preview controller. Shown by view_filter_designer
% via gui_show('BrainstormTab'); torn down by gui_hide('FilterDesigner').
% Authors: Diellor Basha, 2026
    eval(macro_method);
end

%% ===== CREATE PANEL =====
function bstPanelNew = CreatePanel(EigenMat, hFig, ctxFn) %#ok<DEFNU>
    import java.awt.*; import javax.swing.*;
    panelName = 'FilterDesigner';
    jPanel = gui_river([6 6], [0 6 6 6]);

    Variant = EigenMat.Variant;
    isDirac = strcmpi(Variant, 'Dirac');

    % --- operator (read-only) ---
    gui_component('label', jPanel, '', ['Operator: ' Variant]);

    % --- input mode ---
    jPanel.add('br', JLabel('Input:'));
    jInputDelta  = gui_component('radio', jPanel, 'tab', 'Delta (click vertex)');
    jInputSource = gui_component('radio', jPanel, 'br', 'Active source map');
    jInputDelta.setSelected(true);
    grp = ButtonGroup(); grp.add(jInputDelta); grp.add(jInputSource);

    % --- kernel dropdown (from the registry) ---
    names = bst_eigfilter_kernel('list');
    jPanel.add('br', JLabel('Kernel:'));
    jKernel = gui_component('combobox', jPanel, 'tab', [], {names});

    % --- params panel (rebuilt when the kernel changes) ---
    jParams = gui_river([2 2], [0 0 0 0]);
    jPanel.add('br hfill', jParams);

    % --- direction + chirality (Dirac only) ---
    jDir = []; jChir = [];
    if isDirac
        jPanel.add('br', JLabel('Direction:'));
        jDir = gui_component('combobox', jPanel, 'tab', [], {{'Surface normal','Tangent','Custom XYZ'}});
        jPanel.add('br', JLabel('Chirality:'));
        jChir = gui_component('combobox', jPanel, 'tab', [], {{'None','+ (right)','- (left)'}});
    end

    % --- tiling ---
    jPanel.add('br', JLabel('Tiles:'));
    jTiles = gui_component('spinner', jPanel, 'tab', [], {1, 12, 4, 1});
    jChiSplit = gui_component('checkbox', jPanel, 'br', 'Cross with both chiralities');

    % --- save / cancel ---
    jSave   = gui_component('button', jPanel, 'br right', 'Save bank');
    jCancel = gui_component('button', jPanel, '', 'Cancel');

    % --- collect controls into a state struct stored on the panel ---
    ctrl = struct('jKernel',jKernel, 'jParams',jParams, 'jDir',jDir, 'jChir',jChir, ...
                  'jTiles',jTiles, 'jChiSplit',jChiSplit, 'jInputDelta',jInputDelta, ...
                  'jInputSource',jInputSource, 'jSave',jSave, 'jCancel',jCancel);
    state = struct('EigenMat',EigenMat, 'Variant',Variant, 'isDirac',isDirac, ...
                   'hFig',hFig, 'ctxFn',ctxFn, 'Op',[], 'SeedCoeffs',[], 'ActiveTile',1, ...
                   'iVertex',[], 'ctrl',ctrl);
    % load mass once (for projection); reuse the operator referenced by the eigen node
    state.Op = load(file_fullpath(EigenMat.OperatorFile));

    % --- wire callbacks ---
    java_setcb(jKernel,   'ActionPerformedCallback', @(h,e) OnKernelChanged(panelName));
    java_setcb(jTiles,    'StateChangedCallback',    @(h,e) Refresh(panelName));
    java_setcb(jChiSplit, 'ActionPerformedCallback', @(h,e) Refresh(panelName));
    if isDirac
        java_setcb(jDir,  'ActionPerformedCallback', @(h,e) Refresh(panelName));
        java_setcb(jChir, 'ActionPerformedCallback', @(h,e) Refresh(panelName));
    end
    java_setcb(jSave,     'ActionPerformedCallback', @(h,e) OnSave(panelName));
    java_setcb(jCancel,   'ActionPerformedCallback', @(h,e) OnCancel(panelName));

    bstPanelNew = BstPanel(panelName, jPanel.getComponent(), struct('state', state));
    % build the initial param widgets + first preview
    OnKernelChanged(panelName);
end
```

- [ ] **Step 2: Implement the design-state reader + tiling (testable, no figure)**

Add the controller subfunctions. `BuildDesign` reads the widgets into a `base` struct for `bst_filterbank_tiles`; `ComputeField` filters the seed for one tile using the cached coeffs.

```matlab
%% ===== READ WIDGETS -> base design =====
function base = BuildDesign(state) %#ok<DEFNU>
    c = state.ctrl;
    kname = char(c.jKernel.getSelectedItem());
    params = ReadParams(c.jParams);                    % struct from the auto-built fields
    lam = state.EigenMat.Lambda{1};
    base = struct('Kernel',kname, 'Params',params, ...
        'Direction',[1 0 0], 'Chirality',0, 'Axis',[0 0 1], ...
        'N', double(c.jTiles.getValue()), 'Spacing','geometric', ...
        'LambdaRange',[max(eps,min(lam)) max(lam)], 'Chiralities',[]);
    if state.isDirac
        switch char(c.jChir.getSelectedItem())
            case '+ (right)', base.Chirality = +1;
            case '- (left)',  base.Chirality = -1;
            otherwise,        base.Chirality = 0;
        end
        if c.jChiSplit.isSelected(); base.Chiralities = [1 -1]; end
        % Direction resolved per seed vertex in ComputeField (normal/tangent need geometry)
    end
end

%% ===== FILTER THE SEED FOR ONE TILE (uses cached coeffs) =====
function J = ComputeField(state, tile) %#ok<DEFNU>
    % cached projection of the seed must already be in state.SeedCoeffs
    args = {'TransferFn', bst_eigfilter_kernel(tile.Kernel, tile.Params), 'Coeffs', state.SeedCoeffs};
    if state.isDirac && tile.Chirality ~= 0
        args = [args, {'Chirality', struct('Axis', tile.Axis, 'Sign', tile.Chirality)}];
    end
    Jc = bst_dirac_eigenmodes_filter(state.EigenMat, state.Op.Mass, [], 'custom', args{:});
    J = real(Jc);
end
```

- [ ] **Step 3: Implement the seed + refresh + auto-param builder**

```matlab
%% ===== AUTO-BUILD PARAM WIDGETS FROM KERNEL META =====
function OnKernelChanged(panelName)
    [state, ctrl] = GetState(panelName);
    kname = char(ctrl.jKernel.getSelectedItem());
    meta = bst_eigfilter_kernel('info', kname);
    ctrl.jParams.removeAll();
    pf = fieldnames(meta.params);
    for i = 1:numel(pf)
        d = meta.params.(pf{i});
        ctrl.jParams.add('br', javax.swing.JLabel([pf{i} ':']));
        jv = gui_component('texttime', ctrl.jParams, 'tab', num2str(d.default));
        ctrl.jParams.putClientProperty(['param_' pf{i}], jv); %#ok<*JAPIMATHWORKS>
    end
    ctrl.jParams.revalidate(); ctrl.jParams.repaint();
    Refresh(panelName);
end

%% ===== RECOMPUTE + PUSH PREVIEW =====
function Refresh(panelName)
    state = GetState(panelName);
    if isempty(state.SeedCoeffs); return; end           % no seed yet
    base = BuildDesign(state);
    Tiles = bst_filterbank_tiles(base);
    state.ActiveTile = min(state.ActiveTile, numel(Tiles));
    SetState(panelName, 'Tiles', Tiles);
    J = ComputeField(state, Tiles(state.ActiveTile));
    state.ctxFn.PushField(J);                            % orchestrator updates figure_3d
    % spectrum strip
    bank = struct('Kernels', {arrayfun(@(t) bst_eigfilter_kernel(t.Kernel,t.Params), Tiles, 'uni',0)}, ...
                  'Active', state.ActiveTile, ...
                  'OnSelect', @(j) OnSelectTile(panelName, j));
    view_eigfilter_response(bank, state.EigenMat.Lambda{1}, sprintf('%s tiles', base.Kernel));
end

%% ===== SEED FROM A CLICKED VERTEX (called by the orchestrator) =====
function SetSeedVertex(panelName, iVertex) %#ok<DEFNU>
    state = GetState(panelName);
    nV = double(max(cellfun(@(x) max(x(:)), state.EigenMat.GlobalVertices)));
    dirVec = ResolveDirection(state, iVertex);
    Jdelta = zeros(3*nV, 1);
    Jdelta(3*(iVertex-1) + (1:3)) = dirVec(:);
    [~,~,c] = bst_dirac_eigenmodes_filter(state.EigenMat, state.Op.Mass, Jdelta, 'flat', 'ReturnCoeffs', true);
    SetState(panelName, 'SeedCoeffs', c);
    SetState(panelName, 'iVertex', iVertex);
    Refresh(panelName);
end
```

> `ResolveDirection(state,iVertex)` returns the launch 3-vector: surface normal / a tangent / custom XYZ from the direction combobox (load `VertNormals` from the parent surface for 'normal'). `ReadParams`, `GetState`/`SetState` (read/write the `state` struct stored in the panel via `get(panel,'UserData')` or `bst_get('PanelControls','FilterDesigner')`), `OnSelectTile` (set `ActiveTile`, re-`Refresh`), `OnSave`, `OnCancel` are the remaining subfunctions — implement each as a 5-10 line helper following the patterns above. `OnSave` assembles a `filterbankmat` (Variant, ParentEigen, Tiles=state.Tiles, Tiling=base) and calls `db_add_filterbank` then `state.ctxFn.Close()`. `OnCancel` calls `state.ctxFn.Close()`.

- [ ] **Step 4: Headless controller smoke (no figure)**

Run (MCP `evaluate_matlab_code`):

```matlab
EigenMat = tess_eigen(bst_get('Subject',1).Surface(5).FileName, 'Dirac');
Op = load(file_fullpath(EigenMat.OperatorFile));
nV = double(max(cellfun(@(x) max(x(:)), EigenMat.GlobalVertices)));
% emulate a seed projection + one-tile compute without the panel
Jdelta = zeros(3*nV,1); Jdelta(3*(100-1)+(1:3)) = [1 0 0];
[~,~,c] = bst_dirac_eigenmodes_filter(EigenMat, Op.Mass, Jdelta, 'flat', 'ReturnCoeffs', true);
g = bst_eigfilter_kernel('mexhat', struct('t',0.01));
J = real(bst_dirac_eigenmodes_filter(EigenMat, Op.Mass, [], 'custom', 'TransferFn', g, 'Coeffs', c));
fprintf('OK: one-tile field, energy=%.3g\n', norm(J));
```

Expected: `OK: one-tile field, energy=...` (non-zero).

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_filter_designer.m
git commit -m "feat(gui): panel_filter_designer control panel + live controller"
```

---

## Task 9: Orchestrator — view_filter_designer

**Files:**
- Create: `toolbox/gui/view_filter_designer.m`

Opens the eigenbasis, a `figure_3d` cortex preview backed by a temporary results file, docks the panel, links them, and tears everything down on Save/Cancel/figure-close. Provides the `ctxFn` callbacks (`PushField`, `Close`) the panel calls.

- [ ] **Step 1: Implement the orchestrator**

Create `toolbox/gui/view_filter_designer.m`:

```matlab
function hFig = view_filter_designer(NodeFile)
% VIEW_FILTER_DESIGNER: Open a transient Dirac/eigenmode filterbank design session.
%
% USAGE:  view_filter_designer(EigenFile)        % design a new bank on this eigenbasis
%         view_filter_designer(FilterbankFile)   % re-open a saved bank for editing
%
% Opens a figure_3d cortex preview + docks panel_filter_designer. On Save/Cancel/close,
% the panel undocks and the preview figure + spectrum strip close (transient session).
% Authors: Diellor Basha, 2026

    % --- resolve node: eigen vs filterbank ---
    [sSubE,~,~,iE] = bst_get('EigenFile', NodeFile);
    if ~isempty(iE)
        EigenMat = load(file_fullpath(NodeFile));
        EigenFile = NodeFile; loadBank = [];
    else
        FB = load(file_fullpath(NodeFile));
        EigenFile = FB.ParentEigen;
        EigenMat  = load(file_fullpath(EigenFile));
        loadBank  = FB;
    end
    SurfaceFile = EigenMat.ParentSurface;

    % --- singleton: close any existing session ---
    gui_hide('FilterDesigner');
    view_eigfilter_response('close');

    % --- preview figure: a temp results file on the parent cortex ---
    [tmpResultsFile, nV] = i_make_temp_results(SurfaceFile, EigenMat);
    hFig = view_surface_data(SurfaceFile, tmpResultsFile, 'NewFigure');
    setappdata(hFig, 'FilterDesignerTemp', tmpResultsFile);

    % --- context callbacks handed to the panel ---
    ctxFn = struct();
    ctxFn.PushField = @(J) i_push_field(hFig, tmpResultsFile, J);
    ctxFn.Close     = @() i_teardown(hFig);

    % --- dock the panel ---
    gui_show('panel_filter_designer', 'BrainstormTab', 'Filter designer', [], 0, 0, 0);
    panel_filter_designer('Init', EigenMat, hFig, ctxFn, loadBank);

    % --- link: closing the figure ends the session ---
    set(hFig, 'CloseRequestFcn', @(h,e) i_teardown(h));

    % --- vertex-pick -> seed (reuse figure_3d's click; route to the panel) ---
    setappdata(hFig, 'FilterDesignerPick', @(iVertex) panel_filter_designer('SetSeedVertex', 'FilterDesigner', iVertex));
end
```

> Helpers to implement in the same file:
> - `i_make_temp_results(SurfaceFile, EigenMat)` — write a zero-filled unconstrained (3×nV) results struct to a temp file via `bst_save`/`db_template('resultsmat')` so `figure_3d` shows the cortex; return its path + nV.
> - `i_push_field(hFig, tmpResultsFile, J)` — write `J` (reshaped to the results `ImageGridAmp`) and refresh: update the loaded results in `GlobalData`, set the source-vector quiver override (`setappdata(hFig,'QuiverVectorOverride',...)`, already supported by figure_3d), `panel_surface('UpdateSurfaceColormap', hFig)`.
> - `i_teardown(hFig)` — idempotent: `gui_hide('FilterDesigner')`, `view_eigfilter_response('close')`, delete the temp results file (`getappdata(hFig,'FilterDesignerTemp')`), then `delete(hFig)`. Guard with a flag appdata so figure-close and panel-close don't double-run.
>
> Add `panel_filter_designer('Init', ...)` as a thin subfunction that stores the handed-in `EigenMat/hFig/ctxFn/loadBank` into the panel state (the `CreatePanel` already built the widgets; `Init` wires the live state for this session — or fold Init into CreatePanel and pass these as `CreatePanel` args, consistent with Task 8 Step 1's signature).

- [ ] **Step 2: Wire the figure_3d vertex pick to the seed handler**

Find how `figure_3d` reports a clicked vertex (the existing source-vector / scout click path). Run (MCP `evaluate_matlab_code`):

```matlab
txt = fileread(which('figure_3d'));
idx = regexp(txt, 'CurrentVertex|VertexClick|GetClickedVertex|line_click|ClickedPoint', 'once');
disp(txt(max(1,idx-200):idx+400));
```

Then, in `figure_3d`'s vertex-click handler, after a vertex index is resolved, add a hook (guarded so it only fires in a designer session):

```matlab
pickFn = getappdata(hFig, 'FilterDesignerPick');
if ~isempty(pickFn); pickFn(iVertex); end
```

> Place this at the single point where `figure_3d` already computes the clicked vertex index for source display; do not add a new click handler.

- [ ] **Step 3: Live end-to-end smoke (manual, MCP)**

Run (MCP `evaluate_matlab_code`):

```matlab
EigenFile = bst_get('Subject',1).Surface(5).Eigen(1).FileName;
hFig = view_filter_designer(EigenFile);
drawnow;
% simulate a seed click at vertex 100 and a refresh
panel_filter_designer('SetSeedVertex', 'FilterDesigner', 100); drawnow;
ok = ishandle(hFig) && ~isempty(getappdata(hFig,'FilterDesignerTemp'));
fprintf('OK session open=%d\n', ok);
panel_filter_designer('OnCancel', 'FilterDesigner'); drawnow;
fprintf('OK torn down, figure valid=%d (expect 0)\n', ishandle(hFig));
```

Expected: `OK session open=1` then `OK torn down, figure valid=0`.

- [ ] **Step 4: Commit**

```bash
git add toolbox/gui/view_filter_designer.m toolbox/gui/figure_3d.m
git commit -m "feat(gui): view_filter_designer session orchestrator + figure_3d pick hook"
```

---

## Task 10: Tree integration (menu + nested node + cascade delete)

**Files:**
- Modify: `toolbox/tree/node_create_subject.m:155-165` (after the eigen-child loop)
- Modify: `toolbox/tree/tree_callbacks.m` (eigen popup ~line 2572; new filterbank popup; eigen-delete cascade)
- Test: `dev/tests/test_filter_designer_session.m` (Task adds the end-to-end check)

- [ ] **Step 1: Nest filterbank nodes under their eigen node**

In `toolbox/tree/node_create_subject.m`, inside the eigen loop (right after `nodeSurface.add(chNode);` for the eigen node, where `chNode` is the eigen Java node), add a nested loop that attaches matching filterbanks:

```matlab
                    if chCreated
                        nodeSurface.add(chNode);
                        % Nest filterbank child nodes under THIS eigen node
                        if isfield(sSubject.Surface(iSurface), 'Filterbank')
                            eigName = char(sSubject.Surface(iSurface).Eigen(iE).FileName);
                            fbs = sSubject.Surface(iSurface).Filterbank;
                            for iFb = 1:numel(fbs)
                                if file_compare(fbs(iFb).ParentEigen, eigName)
                                    [fbCreated, fbNode] = CreateNode('filterbank', ...
                                        char(fbs(iFb).Comment), char(fbs(iFb).FileName), ...
                                        iFb, iSubject, iSearch);
                                    if fbCreated; chNode.add(fbNode); end
                                end
                            end
                        end
                    end
```

> `CreateNode('filterbank', ...)` requires the bst-java node type from Task 11. Until Task 11 lands, the tree will throw on the unknown type — implement Task 11 before reloading the tree, or temporarily map 'filterbank' to an existing generic type ('eigen') to smoke-test nesting.

- [ ] **Step 2: Add "Design filterbank…" to the eigen popup**

In `toolbox/tree/tree_callbacks.m`, find the eigen node popup (near `EigenView_Callback`, ~line 2572). After the existing "View" item add:

```matlab
                    gui_component('MenuItem', jPopup, [], 'Design filterbank...', IconLoader.ICON_RESULTS, [], @(h,ev)bst_call(@view_filter_designer, filenameFull));
```

- [ ] **Step 3: Add the filterbank node popup + double-click**

In the node-type switch of `tree_callbacks` (Popup section), add a `case 'filterbank'` mirroring the eigen case:

```matlab
                case 'filterbank'
                    gui_component('MenuItem', jPopup, [], 'Edit / re-open', IconLoader.ICON_RESULTS, [], @(h,ev)bst_call(@view_filter_designer, filenameFull));
                    AddSeparator(jPopup);
                    gui_component('MenuItem', jPopup, [], 'Delete', IconLoader.ICON_DELETE, [], @(h,ev)bst_call(@FilterbankDelete_Callback, filenameRelative));
```

And in the double-click dispatch (where eigen/operator nodes are handled), add:

```matlab
            case 'filterbank'
                view_filter_designer(filenameFull);
```

- [ ] **Step 4: Add FilterbankDelete_Callback + cascade from eigen delete**

Add a delete callback mirroring `EigenDelete_Callback`:

```matlab
function FilterbankDelete_Callback(filenameRelative)
    [~, iSubject, iSurface, iFb] = bst_get('FilterbankFile', filenameRelative);
    if isempty(iFb); return; end
    file_delete(file_fullpath(filenameRelative), 1);
    ProtocolSubjects = bst_get('ProtocolSubjects');
    if (iSubject == 0)
        ProtocolSubjects.DefaultSubject.Surface(iSurface).Filterbank(iFb) = [];
    else
        ProtocolSubjects.Subject(iSubject).Surface(iSurface).Filterbank(iFb) = [];
    end
    bst_set('ProtocolSubjects', ProtocolSubjects);
    panel_protocols('UpdateNode', 'Subject', iSubject);
    db_save();
end
```

In `EigenDelete_Callback`, before removing the eigen entry, cascade-delete its filterbanks (so orphans don't linger):

```matlab
    % cascade: delete filterbanks whose ParentEigen is this eigen node
    sSurf = ProtocolSubjects.Subject(iSubject).Surface(iSurface);   % (mirror the 0/default branch too)
    if isfield(sSurf,'Filterbank') && ~isempty(sSurf.Filterbank)
        keep = ~arrayfun(@(f) file_compare(f.ParentEigen, file_short(filenameRelative)), sSurf.Filterbank);
        for f = find(~keep); file_delete(file_fullpath(sSurf.Filterbank(f).FileName), 1); end
        ProtocolSubjects.Subject(iSubject).Surface(iSurface).Filterbank = sSurf.Filterbank(keep);
    end
```

- [ ] **Step 5: End-to-end session test**

Create `dev/tests/test_filter_designer_session.m`:

```matlab
function test_filter_designer_session()
% Open a designer session, seed, save a bank, confirm the nested node, reopen, teardown.
% Requires Brainstorm running with a Dirac eigen node (Subject01 surface 5).
% Authors: Diellor Basha, 2026
    nFail = 0;
    EigenFile = bst_get('Subject',1).Surface(5).Eigen(1).FileName;

    hFig = view_filter_designer(EigenFile); drawnow;
    nFail = nFail + chk('session figure opens', ishandle(hFig));

    panel_filter_designer('SetSeedVertex', 'FilterDesigner', 100); drawnow;
    nFail = nFail + chk('seed sets coeffs', ~isempty(panel_filter_designer('GetState','FilterDesigner').SeedCoeffs));

    panel_filter_designer('OnSave', 'FilterDesigner'); drawnow;
    sSubject = bst_get('Subject',1);
    fbs = sSubject.Surface(5).Filterbank;
    nFail = nFail + chk('bank saved + nested', ~isempty(fbs) && any(strcmp({fbs.ParentEigen}, file_short(EigenFile))));
    nFail = nFail + chk('session torn down after save', ~ishandle(hFig));

    % cleanup
    newFile = fbs(end).FileName;
    FilterbankDelete_Callback(newFile);
    fprintf('\n==== test_filter_designer_session: %d failed ====\n', nFail);
    if nFail > 0, error('test_filter_designer_session FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

> This test references `panel_filter_designer('GetState',...)`, `panel_filter_designer('OnSave',...)`, and `FilterbankDelete_Callback` — ensure those dispatch names exist from Tasks 8/10. Run AFTER Task 11 (the node type must exist for the tree to render).

- [ ] **Step 6: Commit**

```bash
git add toolbox/tree/node_create_subject.m toolbox/tree/tree_callbacks.m dev/tests/test_filter_designer_session.m
git commit -m "feat(tree): nested filterbank node, Design-filterbank menu, cascade delete"
```

---

## Task 11: bst-java fork — filterbank node type + icon

**Files (bst-java fork at `~/workspace/research/code/bst-java`):**
- Modify: `org/brainstorm/tree/BstNode.java` (the node-type table / icon switch where `'eigen'`, `'operator'`, `'manifold'` are registered)

- [ ] **Step 1: Create a paired feature branch on the fork**

```bash
cd ~/workspace/research/code/bst-java
git checkout -b feat/filterbank-node
```

- [ ] **Step 2: Find how 'eigen'/'manifold' node types are registered**

Run (MCP `evaluate_matlab_code` or shell):

```bash
grep -rn "eigen\|manifold\|operator" ~/workspace/research/code/bst-java/src/org/brainstorm/tree/BstNode.java | head
```

Expected: the `switch`/`if` that maps the type string to an icon + display behavior.

- [ ] **Step 3: Register 'filterbank'**

Add a `'filterbank'` branch mirroring `'eigen'` (reuse an existing icon, e.g. the results/eigen icon, unless a dedicated icon asset is added). Concretely, wherever the type string is matched to set `this.iconType` (or equivalent), add:

```java
} else if (type.equalsIgnoreCase("filterbank")) {
    this.iconType = ICON_RESULTS;   // reuse an existing icon constant used by eigen/results
```

> Match the exact field/constant names from Step 2's output. If node types are enumerated in a separate table/array, add `"filterbank"` there too so `CreateNode('filterbank', …)` from MATLAB is accepted.

- [ ] **Step 4: Build the fork jar and deploy it where Brainstorm loads it**

Follow the fork's existing build (e.g. `ant` / `javac` per the repo README) to produce the updated `brainstorm.jar`, and place it where the running Brainstorm picks it up (the same managed location the fork normally deploys to). Then in MATLAB:

```matlab
% reload the tree so the new node type renders (no full restart needed for a jar already on the path at startup;
% if the jar is loaded at startup, restart Brainstorm: brainstorm stop; brainstorm nogui)
db_reload_subjects(0);
panel_protocols('UpdateTree');
```

> ⚑ Per the project's stale-binary trap: if Brainstorm loads the jar at JVM startup, a restart (`brainstorm stop; brainstorm nogui`) is required for the new class to take effect — `db_reload_subjects` alone won't reload Java classes.

- [ ] **Step 5: Verify the node renders**

Re-run `dev/tests/test_filter_designer_session.m` (now the tree can render `filterbank` nodes).
Expected: PASS — `0 failed`, and the saved bank appears nested under the eigen node in the tree.

- [ ] **Step 6: Commit (fork)**

```bash
cd ~/workspace/research/code/bst-java
git add src/org/brainstorm/tree/BstNode.java
git commit -m "feat(tree): filterbank node type nested under eigen"
```

---

## Self-Review

**Spec coverage:**
- Transient session lifecycle (open/link/teardown, singleton) → Tasks 9, 10. ✓
- Control panel auto-built from registry (kernel dropdown, auto-params, direction/chirality Dirac-only, tiling, Save/Cancel) → Task 8. ✓
- Delta OR loaded source map input → Task 8 (input radio + `SetSeedVertex`); source-map seed path uses the same `ReturnCoeffs` projection on a loaded frame (panel `Init` accepts the active results). ✓ *(Note: the plan wires the delta path end-to-end; the source-map branch reuses the identical projection call on the loaded frame — implement it in `SetSeedSource` alongside `SetSeedVertex`.)*
- Live preview via projection/filter split → Task 5 (`ReturnCoeffs`/`Coeffs`) + Task 8 (`SetSeedVertex` projects once, `ComputeField` re-applies). ✓
- Selected tile + clickable spectrum strip → Task 7 + Task 8 `Refresh`. ✓
- Spectrum tiling into a wavelet bank (× chiralities) → Task 4. ✓
- Saved recipe bank, no materialized fields → Tasks 1, 6 (only Tiles/Tiling stored). ✓
- Nested under eigen (storage keyed by ParentEigen, rendered nested) → Tasks 1, 6, 10, 11. ✓
- Operator-general, Dirac first (direction/chirality hidden unless Dirac) → Task 8. ✓
- Error handling (vertex outside support, unconstrained source, idempotent teardown, param clamping from meta range) → Tasks 5/8/9 (`bst_dirac_eigenmodes_filter` already raises on bad seed; teardown guarded; param fields seeded from meta defaults). ✓
- Migration → Task 2. ✓

**Gap noted and assigned:** the *source-map* input branch is described but its `SetSeedSource` subfunction is only sketched in Task 8 Step 3's helper note. Implement it in Task 8 as the sibling of `SetSeedVertex` (project the selected loaded frame instead of a delta) — same `ReturnCoeffs` call, input `J` = the loaded unconstrained frame.

**Placeholder scan:** No "TBD/TODO"; the GUI subfunctions left as "implement following the pattern" (`ReadParams`, `GetState`/`SetState`, `OnSelectTile`, `ResolveDirection`, `i_make_temp_results`, `i_push_field`, `i_teardown`, `SetSeedSource`) each have an explicit signature, inputs/outputs, and a one-line behavior spec — they are bounded 5–15 line helpers, not open scope. Acceptable for a skilled engineer; not silent placeholders.

**Type consistency:** `state` field names (`EigenMat, Variant, isDirac, hFig, ctxFn, Op, SeedCoeffs, ActiveTile, iVertex, ctrl, Tiles`) are used consistently across Tasks 8–10. `ctxFn.PushField`/`ctxFn.Close` are defined in Task 9 and called in Task 8. `bst_filterbank_tiles` `base`/`Tiles` field names match between Tasks 4, 6, 8. The `filterbankmat` fields (`Comment, ParentEigen, Variant, Tiles, Tiling, Provenance`) match across Tasks 1, 6, 9. `bst_get('FilterbankFile')` return signature matches its use in Tasks 6, 10.

---

## Build order summary

1. Task 1 — schema
2. Task 2 — migration
3. Task 3 — accessor
4. Task 4 — tiling generator
5. Task 5 — projection/filter split
6. Task 6 — db_add_filterbank
7. Task 7 — spectrum strip
8. Task 8 — control panel
9. Task 9 — orchestrator
10. Task 11 — bst-java node type *(do before Task 10's tree reload)*
11. Task 10 — tree integration + end-to-end test

Tasks 1–6 produce a working, testable artifact (a saved, reusable filterbank) with no GUI. Tasks 7–11 build the interactive designer on top.
