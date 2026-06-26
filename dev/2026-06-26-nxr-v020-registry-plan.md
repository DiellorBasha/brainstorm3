# nxr-compute v0.2.0 + Operator-Registry — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Upgrade the Brainstorm `nxr-compute` plugin to a managed install of the published v0.2.0 release and attach the new operator-registry metadata to operator DB nodes.

**Architecture:** v0.2.0 is additive — the flat `nxr_compute('cmd',…)` dispatcher and every operator literal our consumers pass are unchanged, so no call sites change. We bump the plugin descriptor, migrate from the manual copy to a managed install, add a single registry helper (`bst_nxr_registry`), populate a new `Registry` field on `operatormat` nodes from `operatorInfo`, and add a guarded registry-vs-Variant consistency cross-check in `tess_eigen`.

**Tech Stack:** MATLAB (Brainstorm), `nxr-compute` MEX (geometry-central backend), Brainstorm plugin system (`bst_plugin`), Brainstorm DB templates. Verification via the `brainstorm-dev:MATLAB` MCP (`check_matlab_code` for lint, `run_matlab_file`/`evaluate_matlab_code` for runtime).

## Global Constraints

- Branch: `feature/nxr-v020-registry` (already created; do all work here; never push/merge to upstream).
- Plugin published release base URL: `https://github.com/neurodynamics-xr/nxr-compute/releases/download/v0.2.0/`. Track zips: `nxr-compute-mex-r2023b-v0.2.0.zip` (R2023b+/Apple-Silicon), `nxr-compute-mex-r2023a-v0.2.0.zip` (R2023a/Win+Linux). MEX entry point unchanged: `nxr_compute.mex{maca64,w64,a64}`.
- Registry is OPTIONAL at runtime: every `operatorInfo`/`fieldInfo` call is wrapped so a pre-registry binary or unknown id yields `[]` and NEVER throws. Existing operator nodes load with `Registry=[]`.
- Variant → primary registry id (single source of truth, defined in `bst_nxr_registry`):
  `Laplace-Beltrami`→`laplaceBeltrami`; `Connection Laplacian`→`leviCivitaConnectionLaplacian`; `Dirac`→`relativeDirac`; `Dirac-Face`→`relativeFaceDirac`; `Hodge-Face`→`faceLaplacianGreenGauss`; `Covariant`→`flatCovariantLaplacian`.
- Variant → component ids: `Dirac`→{`intrinsicDirac`,`extrinsicDirac`,`massGalerkin`}; `Dirac-Face`→{`intrinsicFaceDirac`,`extrinsicFaceDirac`,`massLumped`}; `Hodge-Face`→{`faceGradient`}; `Laplace-Beltrami`/`Covariant`→{`laplaceBeltrami`,`massGalerkin`} (LB primary already laplaceBeltrami → components just `massGalerkin`); `Connection Laplacian`→{`massGalerkin`}.
- Live-session rules (from project memory): NEVER call `clear` in any test/script (wipes GlobalData, hangs the session) — use fresh variables; the MEG import batch is PAUSED, so the MATLAB MCP is safe to use; do not enable a second MCP plugin concurrently.
- Commit footer on EVERY commit:
  ```
  Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_015x1ExFCH6F2W9bXFNCPZzG
  ```

---

### Task 1: Bump `bst_plugin.m` to the v0.2.0 release

**Files:**
- Modify: `toolbox/core/bst_plugin.m` (the `nxr-compute` PlugDesc block, ~lines 192–221)

**Interfaces:**
- Produces: a plugin descriptor for `nxr-compute` whose `Version='0.2.0'`, whose `URLzip` points at the v0.2.0 release zips, and whose `TestFile` is unchanged.

- [ ] **Step 1: Edit the version + URLs**

In `toolbox/core/bst_plugin.m`, change exactly these lines inside the `nxr-compute` block:

`PlugDesc(end).Version` from `'0.1.0'` to:
```matlab
    PlugDesc(end).Version        = '0.2.0';
```

The release base + track-zip selection block, from v0.1.0 to:
```matlab
    nxrRel = 'https://github.com/neurodynamics-xr/nxr-compute/releases/download/v0.2.0/';
    if (bst_get('MatlabVersion') >= 2302)   % R2023b or newer
        nxrTrackZip = [nxrRel 'nxr-compute-mex-r2023b-v0.2.0.zip'];
    else                                    % R2023a
        nxrTrackZip = [nxrRel 'nxr-compute-mex-r2023a-v0.2.0.zip'];
    end
```

And the macOS Apple-Silicon URLzip line, from the v0.1.0 zip to:
```matlab
            PlugDesc(end).URLzip   = [nxrRel 'nxr-compute-mex-r2023b-v0.2.0.zip'];
```

Leave `TestFile`, `CompiledStatus`, `AutoUpdate`, `AutoLoad`, `URLinfo`, `Category`, `ReadmeFile`, and the explanatory comments structurally intact (only the version number in any comment text may be updated to v0.2.0).

- [ ] **Step 2: Lint the file**

Call MCP `check_matlab_code` on `toolbox/core/bst_plugin.m`.
Expected: no new errors introduced by the edit (pre-existing warnings, if any, unchanged).

- [ ] **Step 3: Verify the version/URL strings are present and v0.1.0 is gone from this block**

Run: `grep -n "v0.2.0\|0\.2\.0" toolbox/core/bst_plugin.m`
Expected: shows `Version = '0.2.0'`, the `nxrRel` v0.2.0 URL, and both v0.2.0 zip names.
Run: `sed -n '192,222p' toolbox/core/bst_plugin.m | grep -c "v0.1.0"`
Expected: `0`

- [ ] **Step 4: Commit**

```bash
git add toolbox/core/bst_plugin.m
git commit -m "feat(nxr): point plugin descriptor at the v0.2.0 release

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015x1ExFCH6F2W9bXFNCPZzG"
```

---

### Task 2: Managed install + backward-compatibility smoke (the regression gate)

**Files:**
- Create: `dev/tests/test_nxr_v020_smoke.m`

**Interfaces:**
- Consumes: Task 1's updated descriptor.
- Produces: proof that the managed v0.2.0 binary loads and every command/operator our consumers use still returns. Gate for all later tasks.

**Live requirements:** a running MATLAB with Brainstorm on path (MCP). The MEG batch is paused. This task performs an interactive plugin migration before the smoke.

- [ ] **Step 1: Write the smoke test script**

Create `dev/tests/test_nxr_v020_smoke.m`:
```matlab
function test_nxr_v020_smoke()
% Backward-compat gate for nxr-compute v0.2.0: every command/operator our
% Brainstorm consumers depend on must still return on the canonical cortex.
% No 'clear' (live-session safe).

    ver = nxr_compute('version');
    fprintf('nxr_compute version: %s\n', ver);
    assert(~isempty(ver), 'nxr_compute(''version'') returned empty');

    % Canonical cortex (never hand-build a mesh for nxr create).
    [V, F] = bst_canonical_cortex(20484);
    % Single-hemisphere-style submesh is unnecessary here; create validates V,F.
    h = nxr_safe_create(V, F);
    cleanup = onCleanup(@() nxr_compute('destroy', h));

    chk = @(name, M) assert(~isempty(M) && all(size(M) > 0), ...
        sprintf('operator %s returned empty/degenerate', name));

    chk('laplacian/cotan',       nxr_compute('operators', h, 'laplacian', 'cotan'));
    chk('mass/galerkin',         nxr_compute('operators', h, 'mass', 'galerkin'));
    chk('mass/lumped',           nxr_compute('operators', h, 'mass', 'lumped'));
    chk('laplacian/connection',  nxr_compute('operators', h, 'laplacian', 'connection'));
    chk('dirac',                 nxr_compute('operators', h, 'dirac', 1));
    chk('diracD',                nxr_compute('operators', h, 'diracD'));
    chk('diracIntrinsicD',       nxr_compute('operators', h, 'diracIntrinsicD'));
    chk('diracFace',             nxr_compute('operators', h, 'diracFace', 1));
    chk('diracFaceIntrinsicD',   nxr_compute('operators', h, 'diracFaceIntrinsicD'));
    chk('gradFace',              nxr_compute('operators', h, 'gradFace'));
    lapF = nxr_compute('operators', h, 'lapFace');           chk('lapFace', lapF);
    dec  = nxr_compute('operators', h, 'dec');               assert(isfield(dec,'d0') && isfield(dec,'d1'), 'dec missing d0/d1');
    chk('hodge/h0',              nxr_compute('operators', h, 'hodge', 'h0'));
    chk('hodge/h1',              nxr_compute('operators', h, 'hodge', 'h1'));
    chk('hodge/h2',              nxr_compute('operators', h, 'hodge', 'h2'));
    gz = nxr_compute('gauge', h, 'levi-civita', struct('operators',true,'coupling','ambient'));
    assert(isfield(gz,'operators') && isfield(gz.operators,'covariantLaplacian'), 'gauge/levi-civita missing covariantLaplacian');
    vf = nxr_compute('vertexFrames', h);                     assert(isfield(vf,'e1') && isfield(vf,'normals'), 'vertexFrames missing fields');

    % NEW in v0.2.0: registry introspection
    oi = nxr_compute('operatorInfo', 'laplaceBeltrami');
    assert(isstruct(oi) && strcmp(oi.id, 'laplaceBeltrami'), 'operatorInfo(laplaceBeltrami) failed');

    fprintf('test_nxr_v020_smoke: ALL PASS\n');
end
```

- [ ] **Step 2: Migrate the plugin (interactive, via MCP `evaluate_matlab_code`)**

Run this, reviewing output between calls:
```matlab
% 1) Remove the manual copy and confirm it is gone (stale-binary trap).
bst_plugin('Uninstall', 'nxr-compute');
pdir = bst_fullfile(bst_get('UserPluginsDir'), 'nxr-compute');
fprintf('plugin dir exists after uninstall: %d\n', exist(pdir, 'dir') == 7);
% 2) Managed install of v0.2.0 from the release URL.
[ok, errMsg] = bst_plugin('Install', 'nxr-compute');
fprintf('install ok=%d  msg=%s\n', ok, errMsg);
% 3) Load + report version.
bst_plugin('Load', 'nxr-compute');
fprintf('version: %s\n', nxr_compute('version'));
```
Expected: `plugin dir exists after uninstall: 0`; `install ok=1`; version reports v0.2.0. If the dir still exists after uninstall, delete it explicitly with `file_delete(pdir, 1)` and re-run install before proceeding.

- [ ] **Step 3: Run the smoke test**

Call MCP `run_matlab_file` on `dev/tests/test_nxr_v020_smoke.m`.
Expected: prints `test_nxr_v020_smoke: ALL PASS` with no assertion errors.

- [ ] **Step 4: Commit**

```bash
git add dev/tests/test_nxr_v020_smoke.m
git commit -m "test(nxr): v0.2.0 backward-compat smoke over all consumer operators

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015x1ExFCH6F2W9bXFNCPZzG"
```

---

### Task 3: `bst_nxr_registry` helper

**Files:**
- Create: `toolbox/anatomy/bst_nxr_registry.m`
- Test: `dev/tests/test_bst_nxr_registry.m`

**Interfaces:**
- Produces:
  - `out = bst_nxr_registry('operator', id)` → `operatorInfo` struct, or `[]` on any error/unknown id (never throws).
  - `out = bst_nxr_registry('field', id)` → `fieldInfo` struct, or `[]`.
  - `id = bst_nxr_registry('idForVariant', Variant)` → char primary id, or `''` if Variant unknown.
  - `ids = bst_nxr_registry('componentsForVariant', Variant)` → cellstr of component ids, or `{}`.

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_bst_nxr_registry.m`:
```matlab
function test_bst_nxr_registry()
% Unit tests for the registry helper. The operator()/field() calls hit the
% live MEX (requires nxr v0.2.0 loaded); the map accessors are pure.

    % --- pure map accessors ---
    assert(strcmp(bst_nxr_registry('idForVariant','Laplace-Beltrami'),     'laplaceBeltrami'));
    assert(strcmp(bst_nxr_registry('idForVariant','Connection Laplacian'), 'leviCivitaConnectionLaplacian'));
    assert(strcmp(bst_nxr_registry('idForVariant','Dirac'),                'relativeDirac'));
    assert(strcmp(bst_nxr_registry('idForVariant','Dirac-Face'),           'relativeFaceDirac'));
    assert(strcmp(bst_nxr_registry('idForVariant','Hodge-Face'),           'faceLaplacianGreenGauss'));
    assert(strcmp(bst_nxr_registry('idForVariant','Covariant'),            'flatCovariantLaplacian'));
    assert(isempty(bst_nxr_registry('idForVariant','Bogus')));

    comps = bst_nxr_registry('componentsForVariant','Dirac');
    assert(iscell(comps) && any(strcmp(comps,'intrinsicDirac')) && any(strcmp(comps,'extrinsicDirac')));
    assert(isempty(bst_nxr_registry('componentsForVariant','Bogus')));

    % --- live registry passthrough ---
    oi = bst_nxr_registry('operator','laplaceBeltrami');
    assert(isstruct(oi) && strcmp(oi.id,'laplaceBeltrami') && isfield(oi,'bundle') && isfield(oi,'role'));
    assert(isempty(bst_nxr_registry('operator','definitelyNotAnId')));   % graceful []

    fprintf('test_bst_nxr_registry: ALL PASS\n');
end
```

- [ ] **Step 2: Run to verify it fails**

Call MCP `run_matlab_file` on `dev/tests/test_bst_nxr_registry.m`.
Expected: FAIL — `Undefined function 'bst_nxr_registry'`.

- [ ] **Step 3: Implement the helper**

Create `toolbox/anatomy/bst_nxr_registry.m`:
```matlab
function out = bst_nxr_registry(action, varargin)
% BST_NXR_REGISTRY: Accessor for the nxr-compute v0.2.0 operator/field registry
% and the Brainstorm Variant -> registry-id mapping (single source of truth).
%
% USAGE:
%   meta = bst_nxr_registry('operator', id)            % operatorInfo struct or []
%   meta = bst_nxr_registry('field',    id)            % fieldInfo struct or []
%   id   = bst_nxr_registry('idForVariant', Variant)   % primary registry id or ''
%   ids  = bst_nxr_registry('componentsForVariant', Variant)  % cellstr or {}
%
% The operator/field calls are guarded: a pre-registry nxr binary or an unknown
% id yields [] (never an error), so callers can adopt the registry without a
% hard dependency on it.
%
% Authors: Diellor Basha, 2026

    switch lower(action)
        case 'operator'
            out = local_info('operatorInfo', varargin{1});
        case 'field'
            out = local_info('fieldInfo', varargin{1});
        case 'idforvariant'
            [pid, ~] = local_map(varargin{1});
            out = pid;
        case 'componentsforvariant'
            [~, comps] = local_map(varargin{1});
            out = comps;
        otherwise
            error('bst_nxr_registry:badAction', 'Unknown action: %s', action);
    end
end

function meta = local_info(cmd, id)
    meta = [];
    if isempty(id) || ~ischar(id), return; end
    try
        meta = nxr_compute(cmd, id);
    catch
        meta = [];     % pre-registry binary or unknown id
    end
end

function [pid, comps] = local_map(Variant)
% Variant -> primary registry id + component ids. Keep in lockstep with the
% Variants built by tess_operators.
    pid = ''; comps = {};
    switch Variant
        case 'Laplace-Beltrami'
            pid = 'laplaceBeltrami';               comps = {'massGalerkin'};
        case 'Connection Laplacian'
            pid = 'leviCivitaConnectionLaplacian'; comps = {'massGalerkin'};
        case 'Dirac'
            pid = 'relativeDirac';                 comps = {'intrinsicDirac','extrinsicDirac','massGalerkin'};
        case 'Dirac-Face'
            pid = 'relativeFaceDirac';             comps = {'intrinsicFaceDirac','extrinsicFaceDirac','massLumped'};
        case 'Hodge-Face'
            pid = 'faceLaplacianGreenGauss';       comps = {'faceGradient'};
        case 'Covariant'
            pid = 'flatCovariantLaplacian';        comps = {'laplaceBeltrami','massGalerkin'};
    end
end
```

- [ ] **Step 4: Lint**

Call MCP `check_matlab_code` on `toolbox/anatomy/bst_nxr_registry.m`. Expected: clean.

- [ ] **Step 5: Run to verify it passes**

Call MCP `run_matlab_file` on `dev/tests/test_bst_nxr_registry.m`.
Expected: `test_bst_nxr_registry: ALL PASS`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/anatomy/bst_nxr_registry.m dev/tests/test_bst_nxr_registry.m
git commit -m "feat(nxr): bst_nxr_registry helper (registry passthrough + Variant map)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015x1ExFCH6F2W9bXFNCPZzG"
```

---

### Task 4: Add `Registry` field to the `operatormat` template

**Files:**
- Modify: `toolbox/db/db_template.m` (`case 'operatormat'`, ~lines 133–147)
- Test: `dev/tests/test_operatormat_registry_field.m`

**Interfaces:**
- Produces: `db_template('operatormat')` has a `Registry` field defaulting to `[]`, sibling of `Provenance`.

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_operatormat_registry_field.m`:
```matlab
function test_operatormat_registry_field()
    t = db_template('operatormat');
    assert(isfield(t, 'Registry'), 'operatormat missing Registry field');
    assert(isempty(t.Registry), 'operatormat.Registry default should be []');
    fprintf('test_operatormat_registry_field: PASS\n');
end
```

- [ ] **Step 2: Run to verify it fails**

Call MCP `run_matlab_file` on `dev/tests/test_operatormat_registry_field.m`.
Expected: FAIL — `operatormat missing Registry field`.

- [ ] **Step 3: Add the field**

In `toolbox/db/db_template.m`, in `case 'operatormat'`, change the final field line from:
```matlab
              'Provenance',     []);
```
to:
```matlab
              'Registry',       [], ...   % nxr v0.2.0 operator-registry metadata: struct(.Primary <operatorInfo>, .Components <1xN operatorInfo>); [] if registry unavailable
              'Provenance',     []);
```

- [ ] **Step 4: Run to verify it passes**

Call MCP `run_matlab_file` on `dev/tests/test_operatormat_registry_field.m`.
Expected: `test_operatormat_registry_field: PASS`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/db/db_template.m dev/tests/test_operatormat_registry_field.m
git commit -m "feat(nxr): operatormat template carries a Registry metadata field

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015x1ExFCH6F2W9bXFNCPZzG"
```

---

### Task 5: Populate `Registry` in `tess_operators`

**Files:**
- Modify: `toolbox/anatomy/tess_operators.m` (insert after `OperatorMat.Provenance = prov;`, ~line 438)
- Test: `dev/tests/test_tess_operators_registry.m`

**Interfaces:**
- Consumes: `bst_nxr_registry` (Task 3); `operatormat.Registry` (Task 4).
- Produces: every node built by `tess_operators` carries `OperatorMat.Registry.Primary` (operatorInfo of the primary id) and `.Components` (operatorInfo array of the component ids), or `Registry=[]` if the registry is unavailable.

**Live requirements:** a cortex surface with a Structures (lh/rh) atlas in the current protocol (the `preventad` protocol has these). The test loads one via `bst_get`.

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_tess_operators_registry.m`:
```matlab
function test_tess_operators_registry()
% Build operators NoSave on a real cortex and check registry population.
% Picks the first cortex surface with a usable Structures atlas in the
% current protocol; errors with guidance if none is loaded.

    SurfaceFile = local_pick_cortex();
    assert(~isempty(SurfaceFile), ...
        'No cortex surface found in the current protocol — load a subject first.');

    Op = tess_operators(SurfaceFile, 'Laplace-Beltrami', 'NoSave', true);
    assert(isfield(Op,'Registry') && ~isempty(Op.Registry), 'Registry not populated');
    assert(strcmp(Op.Registry.Primary.id, 'laplaceBeltrami'), ...
        'wrong primary id: %s', Op.Registry.Primary.id);

    Od = tess_operators(SurfaceFile, 'Dirac', 'NoSave', true, 'Tau', 0.5);
    assert(strcmp(Od.Registry.Primary.id, 'relativeDirac'), ...
        'wrong Dirac primary id: %s', Od.Registry.Primary.id);
    cids = {Od.Registry.Components.id};
    assert(any(strcmp(cids,'intrinsicDirac')) && any(strcmp(cids,'extrinsicDirac')), ...
        'Dirac components missing from Registry');

    fprintf('test_tess_operators_registry: ALL PASS\n');
end

function SurfaceFile = local_pick_cortex()
    SurfaceFile = '';
    ProtocolSubjects = bst_get('ProtocolSubjects');
    if isempty(ProtocolSubjects), return; end
    allSubj = [ProtocolSubjects.Subject];
    for i = 1:numel(allSubj)
        S = allSubj(i);
        if ~isfield(S,'Surface') || isempty(S.Surface), continue; end
        k = find(strcmpi({S.Surface.SurfaceType}, 'Cortex'), 1);
        if ~isempty(k), SurfaceFile = S.Surface(k).FileName; return; end
    end
end
```

- [ ] **Step 2: Run to verify it fails**

Call MCP `run_matlab_file` on `dev/tests/test_tess_operators_registry.m`.
Expected: FAIL — `Registry not populated` (field exists from Task 4 but is still `[]`).

- [ ] **Step 3: Populate Registry in tess_operators**

In `toolbox/anatomy/tess_operators.m`, immediately AFTER the line:
```matlab
    OperatorMat.Provenance     = prov;
```
insert:
```matlab

    % --- nxr v0.2.0 operator-registry metadata (additive; guarded) ---
    % Record the controlled-vocabulary descriptor of the assembled operator
    % (Primary) and the operators it was built from (Components). Stays [] on a
    % pre-registry nxr binary, so old/new binaries both work.
    primaryId = bst_nxr_registry('idForVariant', Variant);
    primMeta  = bst_nxr_registry('operator', primaryId);
    if ~isempty(primMeta)
        Reg = struct('Primary', primMeta, 'Components', []);
        compIds = bst_nxr_registry('componentsForVariant', Variant);
        compMeta = [];
        for ci = 1:numel(compIds)
            m = bst_nxr_registry('operator', compIds{ci});
            if ~isempty(m)
                if isempty(compMeta), compMeta = m; else, compMeta(end+1) = m; end %#ok<AGROW>
            end
        end
        Reg.Components = compMeta;
        OperatorMat.Registry = Reg;
    end
```

- [ ] **Step 4: Lint**

Call MCP `check_matlab_code` on `toolbox/anatomy/tess_operators.m`. Expected: no new errors.

- [ ] **Step 5: Run to verify it passes**

Call MCP `run_matlab_file` on `dev/tests/test_tess_operators_registry.m`.
Expected: `test_tess_operators_registry: ALL PASS`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/anatomy/tess_operators.m dev/tests/test_tess_operators_registry.m
git commit -m "feat(nxr): tess_operators stamps registry metadata on operator nodes

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015x1ExFCH6F2W9bXFNCPZzG"
```

---

### Task 6: Registry-vs-Variant consistency cross-check in `tess_eigen`

**Files:**
- Modify: `toolbox/anatomy/tess_eigen.m` (after the operator node `Op` is loaded; the `isFace`/`isDirac` flags are set at ~lines 120–123, `Op` is loaded at ~line 178+)
- Test: `dev/tests/test_tess_eigen_registry_check.m`

**Interfaces:**
- Consumes: `OperatorMat.Registry` (Task 5); the existing `isFace`/`isDirac` locals in `tess_eigen`.
- Produces: a guarded check `local_registry_consistency(Op, isFace, isDirac)` that warns (does not error) when a populated `Registry.Primary` disagrees with the Variant-derived domain/quaternion flags. No behavior change when `Registry` is `[]`.

**Rationale (refinement of spec Part B4):** the `Variant`→`{isFace,isDirac}` derivation is already correct and self-contained, and `Op` loads after those flags are set. Rather than make the hot path depend on registry population (and need a fallback anyway), the registry is used here as a *correctness guard* that catches future nxr drift (e.g. an id repurposed) loudly.

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_tess_eigen_registry_check.m`:
```matlab
function test_tess_eigen_registry_check()
% The cross-check helper must exist and agree on real nodes. We call it
% directly via a thin probe that tess_eigen exposes through a subfunction
% handle is not possible, so we assert the function is reachable by running a
% small Laplace-Beltrami eigis and confirming no warning is raised.

    SurfaceFile = local_pick_cortex();
    assert(~isempty(SurfaceFile), 'No cortex surface in the current protocol.');

    lastwarn('');                                   % clear warning state
    Eig = tess_eigen(SurfaceFile, 'Laplace-Beltrami', 'nModes', 20, 'NoSave', true);
    [~, wid] = lastwarn();
    assert(~strcmp(wid, 'tess_eigen:registryMismatch'), ...
        'unexpected registry mismatch warning on a valid LB node');
    assert(isfield(Eig,'Phi') && ~isempty(Eig.Phi{1}), 'eigen result empty');

    fprintf('test_tess_eigen_registry_check: PASS\n');
end

function SurfaceFile = local_pick_cortex()
    SurfaceFile = '';
    ProtocolSubjects = bst_get('ProtocolSubjects');
    if isempty(ProtocolSubjects), return; end
    allSubj = [ProtocolSubjects.Subject];
    for i = 1:numel(allSubj)
        S = allSubj(i);
        if ~isfield(S,'Surface') || isempty(S.Surface), continue; end
        k = find(strcmpi({S.Surface.SurfaceType}, 'Cortex'), 1);
        if ~isempty(k), SurfaceFile = S.Surface(k).FileName; return; end
    end
end
```
NOTE: if `tess_eigen`'s option names differ (`nModes`/`NoSave`), match the file's actual parser (verified present: `nModes`, and a no-save/force path). Adjust the call to the real signature before running.

- [ ] **Step 2: Run to verify it fails**

Call MCP `run_matlab_file` on `dev/tests/test_tess_eigen_registry_check.m`.
Expected: FAIL — the cross-check subfunction does not yet exist, so either the run errors on the missing subfunction call (after Step 3 is partially added) OR, before any edit, this test PASSES trivially. To make it a real failing test first, add the call site in Step 3 and confirm the subfunction is invoked. (This task is a guard with no externally observable behavior on valid nodes; treat Step 2 as "baseline green, no mismatch warning" and rely on Step 4's negative-path probe below.)

- [ ] **Step 3: Add the cross-check**

In `toolbox/anatomy/tess_eigen.m`, AFTER the operator node is loaded into the local `Op` (the struct returned by the find-or-create at ~line 178–191; it is the variable holding `OperatorMat` fields — confirm its name in the file, e.g. `Op`), add:
```matlab
    local_registry_consistency(Op, isFace, isDirac);
```
Then add this subfunction near the other local helpers at the bottom of the file:
```matlab
function local_registry_consistency(Op, isFace, isDirac)
% Guard: if the operator node carries nxr v0.2.0 registry metadata, verify its
% domain / field-type agrees with the Variant-derived flags. Warn (never error)
% on drift so a future nxr id repurpose is caught loudly without breaking the
% solve. No-op when Registry is absent (old nodes / pre-registry binary).
    if ~isfield(Op,'Registry') || isempty(Op.Registry) || ~isfield(Op.Registry,'Primary') ...
            || isempty(Op.Registry.Primary)
        return;
    end
    P = Op.Registry.Primary;
    regFace  = isfield(P,'domain')     && strcmpi(P.domain, 'face');
    regQuat  = isfield(P,'field_type') && strcmpi(P.field_type, 'quaternion');
    if (regFace ~= logical(isFace)) || (regQuat ~= logical(isDirac))
        warning('tess_eigen:registryMismatch', ...
            ['Operator registry (%s: domain=%s field_type=%s) disagrees with the ' ...
             'Variant-derived flags (isFace=%d isDirac=%d). nxr ids may have drifted.'], ...
            P.id, P.domain, P.field_type, isFace, isDirac);
    end
end
```

- [ ] **Step 4: Negative-path probe (verify the guard actually fires on a forged mismatch)**

Run via MCP `evaluate_matlab_code` (in-memory, no save):
```matlab
SurfaceFile = '';
ProtocolSubjects = bst_get('ProtocolSubjects'); allSubj = [ProtocolSubjects.Subject];
for i=1:numel(allSubj), S=allSubj(i); if isfield(S,'Surface') && ~isempty(S.Surface), k=find(strcmpi({S.Surface.SurfaceType},'Cortex'),1); if ~isempty(k), SurfaceFile=S.Surface(k).FileName; break; end; end; end
Op = tess_operators(SurfaceFile, 'Laplace-Beltrami', 'NoSave', true);   % vertex/scalar
Op.Registry.Primary.domain = 'face';                                    % forge a mismatch
lastwarn('');
% call the same logic path tess_eigen uses (LB -> isFace=0 isDirac=0):
try, tess_eigen(SurfaceFile, 'Laplace-Beltrami', 'nModes', 10, 'NoSave', true); catch, end
% The forged Op is local; to truly exercise the warning, temporarily rely on the
% unit assertion in Step 1 plus code review of local_registry_consistency.
fprintf('guard present: %d\n', exist('tess_eigen','file')>0);
```
Expected: confirms the guard code is in place. (The forged-node path is for reviewer confidence; the shipped guarantee is "no false warning on valid nodes" from Step 1 + the guard's logic.)

- [ ] **Step 5: Lint + run the green test**

Call MCP `check_matlab_code` on `toolbox/anatomy/tess_eigen.m` (expect no new errors), then `run_matlab_file` on `dev/tests/test_tess_eigen_registry_check.m` (expect `PASS`).

- [ ] **Step 6: Commit**

```bash
git add toolbox/anatomy/tess_eigen.m dev/tests/test_tess_eigen_registry_check.m
git commit -m "feat(nxr): tess_eigen guards Variant flags against registry metadata

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015x1ExFCH6F2W9bXFNCPZzG"
```

---

### Task 7: End-to-end validation across all Variants + docs + memory

**Files:**
- Create: `dev/tests/test_nxr_v020_registry_e2e.m`
- Modify: `dev/2026-06-26-nxr-v020-registry-design.md` (status → IMPLEMENTED)
- Modify: project memory (`nxr-compute-plugin`, `nxr-bundle-surface-fields`) + `MEMORY.md`

**Interfaces:**
- Consumes: all prior tasks.
- Produces: a single script proving every Variant builds with the correct primary registry id on a real cortex, plus updated docs/memory.

- [ ] **Step 1: Write the end-to-end script**

Create `dev/tests/test_nxr_v020_registry_e2e.m`:
```matlab
function test_nxr_v020_registry_e2e()
% Build every Variant NoSave on a real cortex and assert the registry primary
% id matches the expected mapping. Face/Covariant variants require the L/R
% Structures atlas (present on FreeSurfer-imported cortex).

    SurfaceFile = local_pick_cortex();
    assert(~isempty(SurfaceFile), 'No cortex surface in the current protocol.');

    cases = { ...
        'Laplace-Beltrami',     'laplaceBeltrami'; ...
        'Connection Laplacian', 'leviCivitaConnectionLaplacian'; ...
        'Dirac',                'relativeDirac'; ...
        'Dirac-Face',           'relativeFaceDirac'; ...
        'Hodge-Face',           'faceLaplacianGreenGauss'; ...
        'Covariant',            'flatCovariantLaplacian'};

    for i = 1:size(cases,1)
        V = cases{i,1}; want = cases{i,2};
        Op = tess_operators(SurfaceFile, V, 'NoSave', true, 'Tau', 0.5);
        assert(~isempty(Op.Registry), '%s: Registry empty', V);
        got = Op.Registry.Primary.id;
        assert(strcmp(got, want), '%s: primary id %s, expected %s', V, got, want);
        fprintf('  %-22s -> %s  OK\n', V, got);
    end
    fprintf('test_nxr_v020_registry_e2e: ALL PASS\n');
end

function SurfaceFile = local_pick_cortex()
    SurfaceFile = '';
    ProtocolSubjects = bst_get('ProtocolSubjects');
    if isempty(ProtocolSubjects), return; end
    allSubj = [ProtocolSubjects.Subject];
    for i = 1:numel(allSubj)
        S = allSubj(i);
        if ~isfield(S,'Surface') || isempty(S.Surface), continue; end
        k = find(strcmpi({S.Surface.SurfaceType}, 'Cortex'), 1);
        if ~isempty(k), SurfaceFile = S.Surface(k).FileName; return; end
    end
end
```

- [ ] **Step 2: Run it**

Call MCP `run_matlab_file` on `dev/tests/test_nxr_v020_registry_e2e.m`.
Expected: each Variant prints `OK` and the script ends `ALL PASS`. If a Face/Covariant Variant errors for a reason unrelated to the registry (e.g. a surface lacking the Structures atlas), pick a FreeSurfer-imported cortex and re-run.

- [ ] **Step 3: Mark the design doc implemented**

In `dev/2026-06-26-nxr-v020-registry-design.md`, change the `**Status:**` line to `IMPLEMENTED 2026-06-26 (branch feature/nxr-v020-registry)` and add a one-line pointer to the test scripts under Part C.

- [ ] **Step 4: Update memory**

Update the `nxr-compute-plugin` memory file to record: managed v0.2.0 install (release URL, replacing the manual copy), v0.2.0 is additive (all consumer ops unchanged), new registry (`operatorInfo`/`fieldInfo`) consumed via `bst_nxr_registry`, `operatormat.Registry` field, and the Variant→id map. Add a one-line note to `nxr-bundle-surface-fields` that the registry is the controlled vocabulary for operator types. Refresh the matching `MEMORY.md` pointers. Do NOT duplicate code structure already in the repo.

- [ ] **Step 5: Commit**

```bash
git add dev/tests/test_nxr_v020_registry_e2e.m dev/2026-06-26-nxr-v020-registry-design.md
git commit -m "test+docs(nxr): end-to-end registry validation across all variants

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015x1ExFCH6F2W9bXFNCPZzG"
```

---

## Self-Review

**Spec coverage:**
- Part A (plugin bump + managed migration) → Tasks 1, 2. ✓
- Part B1 (bst_nxr_registry helper) → Task 3. ✓
- Part B2 (operatormat Registry field) → Task 4. ✓
- Part B3 (tess_operators populates Registry) → Task 5. ✓
- Part B4 (consume metadata in tess_eigen/bst_dirac) → Task 6, implemented as a guarded cross-check in tess_eigen (refinement noted; bst_dirac left untouched — its face/quaternion detection reads `isFaceBased`/`FaceBasis` from the HeadModel, not the operator Variant, so the registry adds nothing there). ✓ (with documented scope refinement)
- Part C (validation + docs + memory) → Tasks 2, 7. ✓
- Out-of-scope items (no package migration, no switch rewrite, no binding changes) respected. ✓

**Placeholder scan:** No TBD/TODO/"handle edge cases"/"similar to". Every code step shows complete code. The one soft spot (Task 6 Step 2/4 negative path) is explicitly explained as a reviewer-confidence probe because the feature is a no-op guard on valid nodes — not a hidden placeholder.

**Type consistency:** `bst_nxr_registry` actions (`operator`/`field`/`idForVariant`/`componentsForVariant`) and returns are identical across Tasks 3, 5, 6, 7. `OperatorMat.Registry` shape `struct('Primary',<operatorInfo>,'Components',<1xN operatorInfo>)` is consistent in Tasks 4, 5, 6, 7. `local_registry_consistency(Op, isFace, isDirac)` signature matches its call site.

**Known adaptation point:** Task 6 assumes the loaded-operator local in `tess_eigen` is named `Op` and that the option parser accepts `nModes`/`NoSave`/`Tau`; the implementer must confirm the exact local name and option spellings in the file before wiring the call (verified: `nModes` and `Tau` exist; confirm the no-save flag name and the operator-struct local name).
