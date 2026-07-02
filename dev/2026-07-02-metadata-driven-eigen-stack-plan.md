# Metadata-Driven Eigen Stack (SP1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `bst_eigen` route eigenmode↔data-type by operator `(field_type, domain)` metadata instead of a hardcoded `switch Variant`, and make `Dirac-Connectome` a first-class `tess_operators`/`tess_eigen` operator (retiring the panel's inline lift).

**Architecture:** `tess_operators` = operator factory (find-or-create, incl. composed operators), `tess_eigen` = eigen-file factory (structure-aware for lifts), `bst_eigen` = metadata-driven orchestrator with no name-switch and no defaults. A single `FieldSpec` resolver reads `(field_type, domain)` from `Operator.Registry.Primary` and derives component count + embedding; a guarded Phi-layout fallback keeps pre-registry binaries working.

**Tech Stack:** MATLAB R2023b · Brainstorm toolbox (`eval(macro_method)` dispatch) · nxr-compute plugin (operator registry) · headless `matlab -batch` for pure-math tests; MATLAB-MCP + a booted `preventad` protocol for integration tests.

**Design spec:** `dev/2026-07-02-metadata-driven-eigen-stack-design.md`.

## Global Constraints

- **No defaults in the library:** `bst_eigen`/`bst_eigenfilter` take explicit parameters only. `domain` absent in metadata ⇒ resolve to `'vertex'` (a completeness fallback, the one allowed default).
- **Additive + guarded metadata:** every registry/metadata change must leave pre-registry nxr binaries working (mirror the existing `if ~isempty(primMeta)` guard in `tess_operators.m:464`). Keep the Phi-layout inference as a fallback.
- **`field_type` vocabulary:** `'real'` (scalar, width 1, nComponents 1), `'complex'` (tangent, width 1-complex, nComponents 1), `'quaternion'` (width 4, nComponents 3, imag-slot embed).
- **`domain` vocabulary:** `'vertex'` (index `GlobalVertices`) | `'face'` (index `GlobalFaces`).
- **No behavior change in tasks 1–3** (pure de-duplication, regression-pinned); task 4 relocates existing math.
- **Never `clear`** in a live MATLAB/Brainstorm session (use `rehash`); edited `.m` auto-reload.
- **Tests live in** `dev/tests/test_*.m` (standalone scripts with `assert`, mirroring `dev/tests/test_eigen_axes.m` etc.).
- **Commits target `development`**, end with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

## File Structure

| File | Responsibility | Action |
|---|---|---|
| `toolbox/anatomy/bst_nxr_registry.m` | Variant → registry id **and** local `field_type`/`domain` for non-nxr operators | Modify |
| `toolbox/eigen/bst_eigen.m` | `FieldSpec` resolver + metadata-driven `ExtractHemiField`/`nComponents`/project-back | Modify |
| `toolbox/eigen/bst_eigenfilter.m` | `i_fiber` delegates to the shared `FieldSpec` | Modify |
| `toolbox/anatomy/tess_operators.m` | stamp connectome branch; build `Dirac-Connectome` operator node | Modify |
| `toolbox/anatomy/tess_eigen.m` | structure-aware `Dirac-Connectome` eigen (find-or-create LB-Conn + lift) | Modify |
| `toolbox/gui/panel_bst_dynamics.m` | retire inline lift in `i_atom_axes` → call `bst_eigen('Axes','Dirac-Connectome')` | Modify |
| `dev/tests/test_fieldspec.m` | FieldSpec resolution + fallback | Create |
| `dev/tests/test_eigen_metadata_routing.m` | regression: metadata routing == old switch, per variant | Create |
| `dev/tests/test_dirac_connectome_factory.m` | factory eigenbasis == panel lift (1e-12) | Create |

---

## Task 1: `FieldSpec` resolver + guarded fallback

Introduce the single metadata resolver, with no call-site swaps yet (pure addition). It reads `(field_type, domain)` and derives the field layout; falls back to Phi-layout inference for pre-registry operators.

**Files:**
- Modify: `toolbox/eigen/bst_eigen.m` (add `FieldSpec` subfunction + verb dispatch)
- Test: `dev/tests/test_fieldspec.m`

**Interfaces:**
- Produces: `spec = bst_eigen('FieldSpec', ax)` → struct with fields
  `field_type` (`'real'|'complex'|'quaternion'`), `domain` (`'vertex'|'face'`),
  `width` (rows per element: 1|1|4), `nComponents` (source-map comps: 1|1|3),
  `C` (component count for `i_fiber` compatibility: 1|2|4), `kind` (`'scalar'|'tangent'|'quaternion'`).
  `ax` is a `bst_eigen('Axes',…)` struct (has `.Operator.Registry.Primary`, `.Phi`, `.GlobalVertices`).

- [ ] **Step 1: Write the failing test** — `dev/tests/test_fieldspec.m`

```matlab
function test_fieldspec()
% FieldSpec derives (field_type,domain,width,nComponents,C,kind) from operator metadata,
% and falls back to Phi-layout inference when metadata is absent.
addpath(genpath(fullfile(pwd,'toolbox')));
fprintf('--- bst_eigen FieldSpec ---\n');

% helper: fake ax with a given registry field_type/domain + a Phi of the right width
mk = @(ft, dom, width, nGv) struct( ...
    'Operator', struct('Registry', struct('Primary', i_prim(ft,dom))), ...
    'Phi', {{randn(width*nGv, 5), []}}, 'GlobalVertices', {{(1:nGv)', []}});

% metadata-driven cases
s = bst_eigen('FieldSpec', mk('real','vertex',1,50));
assert(strcmp(s.field_type,'real') && strcmp(s.domain,'vertex') && s.width==1 && s.nComponents==1 && s.C==1 && strcmp(s.kind,'scalar'));
s = bst_eigen('FieldSpec', mk('quaternion','vertex',4,50));
assert(strcmp(s.field_type,'quaternion') && s.width==4 && s.nComponents==3 && s.C==4 && strcmp(s.kind,'quaternion'));
s = bst_eigen('FieldSpec', mk('complex','vertex',1,50));
assert(strcmp(s.field_type,'complex') && s.C==2 && s.nComponents==1 && strcmp(s.kind,'tangent'));

% domain default = vertex when metadata omits it
p = i_prim('quaternion',''); p = rmfield(p,'domain');
axNoDom = struct('Operator',struct('Registry',struct('Primary',p)), 'Phi',{{randn(4*50,5),[]}}, 'GlobalVertices',{{(1:50)',[]}});
s = bst_eigen('FieldSpec', axNoDom);  assert(strcmp(s.domain,'vertex'));

% pre-registry fallback: no Registry -> infer C from Phi rows / nV
axFb = struct('Phi',{{randn(4*50,5),[]}}, 'GlobalVertices',{{(1:50)',[]}});
s = bst_eigen('FieldSpec', axFb);  assert(s.C==4 && strcmp(s.kind,'quaternion'));
fprintf('PASS\n');
end
function p = i_prim(ft,dom), p = struct('field_type',ft,'domain',dom); end
```

- [ ] **Step 2: Run the test — expect FAIL** (FieldSpec verb not defined)

```bash
matlab -batch "addpath('dev/tests'); test_fieldspec" 2>&1 | grep -v X11
# Expected: error "Unrecognized ... FieldSpec" or assertion before PASS
```

- [ ] **Step 3: Add the `FieldSpec` verb + resolver** to `toolbox/eigen/bst_eigen.m`

Add near the top verb dispatch (the file uses explicit verb handling; place beside the existing `'Axes'` handling) a branch:
```matlab
if (nargin >= 1) && ischar(OPTIONS) && strcmpi(OPTIONS, 'FieldSpec')
    OutputFiles = i_field_spec(Data);   % Data = ax
    return;
end
```
Then add the subfunction (derived from `bst_eigenfilter`'s `i_fiber`, extended with domain/width/nComponents):
```matlab
%% ===== FIELD SPEC: (field_type, domain) -> layout, from operator metadata =====
function spec = i_field_spec(ax)
    ft = '';  dom = '';
    if isfield(ax,'Operator') && isstruct(ax.Operator) && isfield(ax.Operator,'Registry') ...
            && ~isempty(ax.Operator.Registry) && isfield(ax.Operator.Registry,'Primary') ...
            && ~isempty(ax.Operator.Registry.Primary)
        P = ax.Operator.Registry.Primary;
        if isfield(P,'field_type'), ft  = P.field_type; end
        if isfield(P,'domain'),     dom = P.domain;     end
    end
    switch lower(ft)
        case 'real',       C = 1;
        case 'complex',    C = 2;
        case 'quaternion', C = 4;
        otherwise,         C = [];      % pre-registry -> infer below
    end
    if isempty(C)                        % guarded fallback: derive from the Phi row layout
        nV = numel(ax.GlobalVertices{1});
        C  = round(size(ax.Phi{1},1) / max(nV,1));
        switch C, case 1, ft='real'; case 2, ft='complex'; case 4, ft='quaternion'; otherwise, ft='real'; C=1; end
    end
    if isempty(dom), dom = 'vertex'; end               % completeness default (vertex is universal)
    switch C
        case 1, kind='scalar';     width=1; nComp=1;
        case 2, kind='tangent';    width=1; nComp=1;   % complex row per element
        case 4, kind='quaternion'; width=4; nComp=3;
        otherwise, kind='scalar';  width=1; nComp=1;
    end
    spec = struct('field_type',ft, 'domain',dom, 'width',width, 'nComponents',nComp, 'C',C, 'kind',kind);
end
```

- [ ] **Step 4: Run the test — expect PASS**

```bash
matlab -batch "addpath('dev/tests'); test_fieldspec" 2>&1 | grep -v X11
# Expected: PASS
```

- [ ] **Step 5: Lint + commit**

```bash
matlab -batch "m=checkcode('toolbox/eigen/bst_eigen.m','-id','-string'); disp(m)" 2>&1 | grep -viE "MSNU|NOCOMMA|GVMIS|INUSD|NASGU|AGROW|X11" | head
git add toolbox/eigen/bst_eigen.m dev/tests/test_fieldspec.m
git commit -m "feat(eigen): FieldSpec resolver — (field_type,domain) from operator metadata, guarded fallback

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Complete operator `field_type` metadata

Give every operator a `field_type` — including the connectome family that isn't in the nxr binary — so FieldSpec resolves from metadata, not the fallback.

**Files:**
- Modify: `toolbox/anatomy/bst_nxr_registry.m` (add a local `field_type`/`domain` declaration for non-nxr variants + a `'fieldspec'` verb)
- Modify: `toolbox/anatomy/tess_operators.m` (stamp the connectome early-return branch)
- Test: extend `dev/tests/test_fieldspec.m` with a real-operator metadata check (guarded to run only if a `preventad` surface is available; else skipped) — keep it headless-safe.

**Interfaces:**
- Consumes: `bst_nxr_registry('idForVariant', Variant)` (existing).
- Produces: `ft = bst_nxr_registry('fieldspec', Variant)` → struct `('field_type','domain')` or `[]`. Declares locally: `Laplace-Beltrami`/`LB-Connectome`/`Connectome Laplacian`→`(real,vertex)`; `Dirac`/`Dirac-Connectome`→`(quaternion,vertex)`; `Dirac-Face`/`Hodge-Face`→`(quaternion,face)`; `Connection Laplacian`→`(complex,vertex)`.

- [ ] **Step 1: Write the failing test** — append to `dev/tests/test_fieldspec.m` a new function `test_registry_fieldspec()`:

```matlab
function test_registry_fieldspec()
addpath(genpath(fullfile(pwd,'toolbox')));
fprintf('--- bst_nxr_registry fieldspec ---\n');
want = {'Laplace-Beltrami','real','vertex'; 'LB-Connectome','real','vertex'; ...
        'Dirac','quaternion','vertex'; 'Dirac-Connectome','quaternion','vertex'; ...
        'Dirac-Face','quaternion','face'; 'Connection Laplacian','complex','vertex'};
for i=1:size(want,1)
    fs = bst_nxr_registry('fieldspec', want{i,1});
    assert(~isempty(fs), 'no fieldspec for %s', want{i,1});
    assert(strcmp(fs.field_type, want{i,2}) && strcmp(fs.domain, want{i,3}), 'wrong for %s', want{i,1});
end
fprintf('PASS\n');
end
```

- [ ] **Step 2: Run — expect FAIL**

```bash
matlab -batch "addpath('dev/tests'); test_registry_fieldspec" 2>&1 | grep -v X11
# Expected: FAIL ("no fieldspec for ...")
```

- [ ] **Step 3: Add the `fieldspec` verb + local table** to `bst_nxr_registry.m`

Add a verb branch in the top dispatch and a local map:
```matlab
% in the command switch:
    case 'fieldspec'
        out = local_fieldspec(varargin{1});
```
```matlab
function fs = local_fieldspec(Variant)
% Local (field_type, domain) for EVERY Brainstorm operator variant — including the
% connectome family that has no nxr binary id. Single source of truth for eigen routing.
    fs = [];
    switch Variant
        case {'Laplace-Beltrami','LB-Connectome','Connectome Laplacian'}
            fs = struct('field_type','real',       'domain','vertex');
        case {'Dirac','Dirac-Connectome'}
            fs = struct('field_type','quaternion', 'domain','vertex');
        case {'Dirac-Face','Hodge-Face'}
            fs = struct('field_type','quaternion', 'domain','face');
        case 'Connection Laplacian'
            fs = struct('field_type','complex',    'domain','vertex');
    end
end
```

- [ ] **Step 4: Stamp `field_type` in `tess_operators`** so operator nodes carry it (both the main path and the connectome early-return).

In `tess_operators.m`, after building `OperatorMat.Registry` (`:464–475`) and inside the connectome branch (before its `return` at `:173`), merge the local fieldspec into `Registry.Primary` so it's present even without an nxr id:
```matlab
% after Reg is built (or build a minimal Reg if primMeta was empty):
fs = bst_nxr_registry('fieldspec', Variant);
if ~isempty(fs)
    if ~isfield(OperatorMat,'Registry') || isempty(OperatorMat.Registry)
        OperatorMat.Registry = struct('Primary', struct(), 'Components', []);
    end
    OperatorMat.Registry.Primary.field_type = fs.field_type;
    OperatorMat.Registry.Primary.domain     = fs.domain;
end
```
Place this helper call in both the connectome branch (`:167–173`, before `return`) and the main build (after `:476`), so every operator is stamped.

- [ ] **Step 5: Run — expect PASS; lint; commit**

```bash
matlab -batch "addpath('dev/tests'); test_registry_fieldspec" 2>&1 | grep -v X11   # PASS
git add toolbox/anatomy/bst_nxr_registry.m toolbox/anatomy/tess_operators.m dev/tests/test_fieldspec.m
git commit -m "feat(anatomy): local field_type/domain for all operators (incl. connectome family)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Route `bst_eigen` by `FieldSpec` (retire `switch Variant`)

Swap the four name-based sites to metadata dispatch, pinned by a regression test that compares new output to the old switch for every existing variant.

**Files:**
- Modify: `toolbox/eigen/bst_eigen.m` (`ExtractHemiField` `:398`, complex `:220`, project-back `:238`, `nComponents` `:541/:590`)
- Modify: `toolbox/eigen/bst_eigenfilter.m` (`i_fiber` → delegate to `bst_eigen('FieldSpec',…)`)
- Test: `dev/tests/test_eigen_metadata_routing.m`

**Interfaces:**
- Consumes: `bst_eigen('FieldSpec', ax)` (Task 1).
- Produces: no new public interface; `ExtractHemiField` returns identical `U_h` for all variants.

- [ ] **Step 1: Write the failing regression test** — `dev/tests/test_eigen_metadata_routing.m`

```matlab
function test_eigen_metadata_routing()
% The FieldSpec-driven embedding must reproduce the old switch-Variant U_h byte-for-byte.
% Synthetic per-variant: build a fake EigenMat + source map F, compare embed to a reference.
addpath(genpath(fullfile(pwd,'toolbox')));
fprintf('--- metadata routing == old switch ---\n');
rng(3);  nGv = 40;  nT = 6;  gv = (1:nGv)';
% reference embeddings (copied from the pre-refactor switch bodies)
ref.real = @(F) F(gv,:);
ref.quat = @(F) i_embed4(F, gv, nT);              % [0;x;y;z]
for tc = {{'real',1},{'quaternion',3}}
    ft = tc{1}{1};  ncomp = tc{1}{2};
    F = randn(ncomp*nGv, nT);
    ax = struct('Variant',i_variant(ft), 'Operator',struct('Registry',struct('Primary',struct('field_type',ft,'domain','vertex'))), ...
                'Phi',{{randn(i_width(ft)*nGv,5),[]}}, 'GlobalVertices',{{gv,[]}}, 'Lambda',{{(1:5)',[]}});
    Op = struct('Mass',{{speye(i_width(ft)*nGv),[]}});
    [U_h,~,~,~,msg] = bst_eigen('ExtractHemiFieldTest', F, ax, Op, 1);  % thin test hook, see step 3
    assert(isempty(msg), 'msg: %s', msg);
    if strcmp(ft,'real'), R = ref.real(F); else, R = ref.quat(F); end
    assert(isequal(size(U_h),size(R)) && max(abs(U_h(:)-R(:)))<1e-12, '%s embed mismatch', ft);
end
fprintf('PASS\n');
end
function U=i_embed4(F,gv,nT), n=numel(gv); U=zeros(4*n,nT); U(2:4:end,:)=F((gv-1)*3+1,:); U(3:4:end,:)=F((gv-1)*3+2,:); U(4:4:end,:)=F((gv-1)*3+3,:); end
function w=i_width(ft), if strcmp(ft,'quaternion'), w=4; else, w=1; end, end
function v=i_variant(ft), if strcmp(ft,'quaternion'), v='Dirac'; else, v='Laplace-Beltrami'; end, end
```

- [ ] **Step 2: Run — expect FAIL** (`ExtractHemiFieldTest` hook not defined)

```bash
matlab -batch "addpath('dev/tests'); test_eigen_metadata_routing" 2>&1 | grep -v X11
```

- [ ] **Step 3: Rewrite `ExtractHemiField` to use `FieldSpec`; add the test hook.**

Replace the `switch EigenMat.Variant` body (`bst_eigen.m:398–463`) with a `FieldSpec`-driven embed:
```matlab
    spec = i_field_spec(struct('Operator',OperatorMat_or_ax(EigenMat,OperatorMat), ...
                               'Phi',{EigenMat.Phi}, 'GlobalVertices',{EigenMat.GlobalVertices}));
    if strcmp(spec.domain,'face'), idxAll = EigenMat.GlobalFaces{h}(:); else, idxAll = EigenMat.GlobalVertices{h}(:); end
    switch spec.field_type
        case 'real'
            if max(idxAll) > size(F,1), msg = sprintf('scalar map needs row %d, has %d.', max(idxAll), size(F,1)); return; end
            U_h = F(idxAll, :);
        case 'complex'
            if max(idxAll) > size(F,1), msg='...'; return; end
            if isreal(F), msg = 'Connection Laplacian needs a COMPLEX tangent field (frame not persisted).'; return; end
            U_h = F(idxAll, :);
        case 'quaternion'
            if 3*max(idxAll) > size(F,1), msg = sprintf('quaternion needs a 3-vector map row %d, has %d.', 3*max(idxAll), size(F,1)); return; end
            n = numel(idxAll);  U_h = zeros(4*n, nT);
            U_h(2:4:end,:) = F((idxAll-1)*3+1, :);
            U_h(3:4:end,:) = F((idxAll-1)*3+2, :);
            U_h(4:4:end,:) = F((idxAll-1)*3+3, :);
    end
```
(Read `EigenMat`'s available fields first; the `Operator` for FieldSpec is `OperatorMat` with `Registry` — verify `OperatorMat.Registry` is carried onto the loaded eigen struct, else pass `EigenMat.Operator`.) Add a thin verb `ExtractHemiFieldTest` that calls `ExtractHemiField` for the test.

- [ ] **Step 4: Swap the other three sites** (`:220`, `:238`, `:541`, `:590`) to `spec = bst_eigen('FieldSpec', ax); spec.field_type`/`spec.nComponents`. E.g. `:541`:
```matlab
spec = bst_eigen('FieldSpec', i_axfromeigen(EigenMat));
nComp = spec.nComponents;    % 3 for quaternion, 1 for real/complex
```

- [ ] **Step 5: Delegate `i_fiber` in `bst_eigenfilter.m`** to the shared resolver:
```matlab
function [C, kind] = i_fiber(ax)
    spec = bst_eigen('FieldSpec', ax);   C = spec.C;  kind = spec.kind;
end
```

- [ ] **Step 6: Run the regression test — expect PASS; lint; commit**

```bash
matlab -batch "addpath('dev/tests'); test_eigen_metadata_routing" 2>&1 | grep -v X11   # PASS
git add toolbox/eigen/bst_eigen.m toolbox/eigen/bst_eigenfilter.m dev/tests/test_eigen_metadata_routing.m
git commit -m "refactor(eigen): route eigenmode<->data-type by FieldSpec, retire switch Variant

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: `Dirac-Connectome` as a factory operator + retire the panel lift

Move the LB-Connectome→quaternion lift into the factories (structure-aware) and delete the panel's inline recipe.

**Files:**
- Modify: `toolbox/anatomy/tess_operators.m` (add `Dirac-Connectome` variant: find-or-create LB-Connectome operator, produce the lifted operator node + metadata)
- Modify: `toolbox/anatomy/tess_eigen.m` (add `Dirac-Connectome` variant map + structure-aware path: find-or-create LB-Connectome eigen, `bst_lift_connectome_dirac`)
- Modify: `toolbox/gui/panel_bst_dynamics.m` (`i_atom_axes` Dirac-Connectome branch → `bst_eigen('Axes','Dirac-Connectome')`)
- Test: `dev/tests/test_dirac_connectome_factory.m`

**Interfaces:**
- Consumes: `[Phiq,Lamq,Mq] = bst_lift_connectome_dirac(Phis,Lams,Ms)`; `tess_eigen(surf,'LB-Connectome')`.
- Produces: `bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac-Connectome',…))` returns a quaternion single-block ax; `tess_eigen(surf,'Dirac-Connectome')` returns `Phi[4nV×3K]`, `Lambda=repelem(λ,3)`, `field_type='quaternion'`.

- [ ] **Step 1: Write the failing equivalence test** — `dev/tests/test_dirac_connectome_factory.m` (needs a booted `preventad` surface; run via MATLAB-MCP or `matlab -batch` with a known surface path)

```matlab
function test_dirac_connectome_factory()
% Factory Dirac-Connectome eigenbasis == the reference lift of the LB-Connectome eigenbasis.
if ~brainstorm('status'), brainstorm nogui; end
iP = bst_get('Protocol','preventad'); gui_brainstorm('SetCurrentProtocol', iP);
surf = 'sub-MTL0005/tess_cortex_pial_low.mat';
axc = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','LB-Connectome','nModes',60));
[Pq,Lq,Mq] = bst_lift_connectome_dirac(axc.Phi{1}, axc.Lambda{1}(:), axc.Mass{1});   % reference
axd = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac-Connectome','nModes',60));  % factory
assert(isequal(size(axd.Phi{1}), size(Pq)) && max(abs(axd.Phi{1}(:)-Pq(:)))<1e-10, 'Phi mismatch');
assert(max(abs(axd.Lambda{1}(:)-Lq(:)))<1e-10, 'Lambda mismatch');
fs = bst_eigen('FieldSpec', axd);  assert(strcmp(fs.field_type,'quaternion'), 'field_type');
fprintf('PASS\n');
end
```

- [ ] **Step 2: Run — expect FAIL** (`bst_eigen('Axes','Dirac-Connectome')` errors: unknown variant in `tess_eigen`)

- [ ] **Step 3: Add the operator to `tess_operators`** — a `Dirac-Connectome` variant that find-or-creates LB-Connectome and records the lift descriptor (Components + field_type). Because the lift is naturally an eigen-level operation, the operator node stores the base reference + `Registry.Primary.field_type='quaternion'`, `Components=[LB-Connectome]`; the actual mode lift happens in `tess_eigen` (Step 4). Add a branch in the connectome region (`:167`):
```matlab
if strcmpi(Variant, 'Dirac-Connectome')
    baseOp = tess_operators(SurfaceFile, 'LB-Connectome', 'Tau', Tau);   % find-or-create base
    OperatorMat = baseOp;                                                % reuse mass/geometry
    OperatorMat.Comment = 'Dirac-Connectome operator';
    OperatorMat.Provenance = struct('Backend','lift','Base','LB-Connectome');
    OperatorMat.Registry = struct('Primary', struct('field_type','quaternion','domain','vertex'), ...
                                  'Components', []);
    if ~NoSave
        [~, iSubjectSave] = bst_get('SurfaceFile', SurfaceFile);
        db_add_operator(iSubjectSave, SurfaceFile, OperatorMat, 'Dirac-Connectome operator');
    end
    return;
end
```

- [ ] **Step 4: Add the structure-aware eigen path to `tess_eigen`** — variant map (`:103–120`) gains `'dirac-connectome' → 'Dirac-Connectome'`; then find-or-create the LB-Connectome eigen file and lift:
```matlab
if strcmpi(Variant, 'Dirac-Connectome')
    base = tess_eigen(SurfaceFile, 'LB-Connectome', 'nModes', nModes, 'ForceRecompute', ForceRecompute);
    Op   = in_bst_operator(base.OperatorFile);                       % LB-Connectome mass
    [Pq,Lq,Mq] = bst_lift_connectome_dirac(base.Phi{1}, base.Lambda{1}(:), Op.Mass{1});
    EigenMat = base;                                                 % clone structure, then override
    EigenMat.Variant = 'Dirac-Connectome';
    EigenMat.Phi     = {Pq, []};
    EigenMat.Lambda  = {Lq, []};
    EigenMat.K       = size(Pq,2);
    % GlobalVertices unchanged (whole-brain single block); OperatorFile -> the Dirac-Connectome node
    EigenMat.OperatorFile = tess_operators(SurfaceFile,'Dirac-Connectome','Tau',Tau,'NoSave',NoSave);  % or its FileName
    if ~NoSave, db_add_eigen(...); end   % follow the file's existing db_add_eigen pattern
    return;
end
```
(Read `tess_eigen`'s existing `db_add_eigen` call + `EigenMat` field set first, and mirror it exactly; `Mq` feeds the operator's `Mass`, already carried via the Dirac-Connectome operator node.)

- [ ] **Step 5: Retire the panel inline lift.** In `panel_bst_dynamics.m` `i_atom_axes`, replace the Dirac-Connectome branch (`:561–578`) body with:
```matlab
if strcmp(variant, 'Dirac-Connectome')
    surf = i_atom_surface(st);  key = ['dconn|' surf];
    Mc = getappdata(0,'DynamicsAtomAx');  if isempty(Mc)||~isa(Mc,'containers.Map'), Mc=containers.Map('KeyType','char','ValueType','any'); end
    if isKey(Mc,key), ax = Mc(key); return; end
    Fs = 100; D = getappdata(st.hFig,'DynamicsOverlay');
    if ~isempty(D) && isfield(D,'srcDS'), try, tv=bst_memory('GetTimeVector',D.srcDS,D.srcResult); if numel(tv)>1, Fs=1/median(diff(tv)); end, catch, end, end
    nF = max(2, round(4*Fs));
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac-Connectome','nModes',60,'TimeWindow',[0 (nF-1)/Fs],'SampleRate',Fs));
    if isempty(ax) || isempty(ax.Phi{1}), return; end
    Mc(key) = ax;  setappdata(0,'DynamicsAtomAx', Mc);
    return;
end
```

- [ ] **Step 6: Run the equivalence test (MCP or matlab -batch) — expect PASS**

```matlab
% via MATLAB-MCP (session already has preventad):
run_matlab_file('dev/tests/test_dirac_connectome_factory.m')   % -> PASS
```

- [ ] **Step 7: Integration check (live) + commit.** In the Dynamics panel on `sub-MTL0002`, build a Dirac-Connectome heat atom + Apply → confirm the fiber-spread cortex (contra-hemi spread) is unchanged.

```bash
git add toolbox/anatomy/tess_operators.m toolbox/anatomy/tess_eigen.m toolbox/gui/panel_bst_dynamics.m dev/tests/test_dirac_connectome_factory.m
git commit -m "feat(anatomy): Dirac-Connectome as a tess_operators/tess_eigen factory operator; retire panel lift

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review

**Spec coverage:** §4.A metadata completion → Task 2; §4.B Dirac-Connectome factory → Task 4; §4.C FieldSpec routing → Tasks 1+3; §3 `(field_type,domain)` model → Task 1 (`i_field_spec`); §5 sequencing → Task order 1→2→3→4. §7 risks: Connection-Laplacian complex/frame handled (Task 3 keeps the `isreal(F)` guard + message; frame persistence out of scope, noted); `fieldInfo` domain availability sidestepped by the **local** `fieldspec` table (Task 2); quaternion width-4/nComp-3 carried in the FieldSpec struct (Task 1). Covered.

**Placeholder scan:** the deep-internal edits in Tasks 3–4 (`ExtractHemiField` rewrite, `tess_eigen` `db_add_eigen` mirror) include a "read the current function first, mirror its exact `db_add_eigen`/field set" instruction rather than a fabricated full-function body — this is a deliberate guard against drift from the real code (the subagent reads it), not a vague TODO; the *new logic* (embeds, lift wiring, metadata) is shown in full. Test code is complete and runnable.

**Type consistency:** `i_field_spec` returns `{field_type,domain,width,nComponents,C,kind}` and every consumer (`i_fiber` uses `.C/.kind`; `nComponents` sites use `.nComponents`; ExtractHemiField uses `.field_type/.domain`) matches. `bst_lift_connectome_dirac` call sites use the verified `[Phiq,Lamq,Mq]` signature.

---

## Execution Handoff

Plan saved to `dev/2026-07-02-metadata-driven-eigen-stack-plan.md`.
