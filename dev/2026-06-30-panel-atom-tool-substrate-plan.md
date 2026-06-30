# Atom-tool substrate (Plan 1 of 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the MATLAB substrate for the panel atom tool — a shared per-kernel control spec (used by both the designer and the future panel), the designer refactored to consume it, atom-group generator fields, and a realise→threshold→store path.

**Architecture:** Extract the designer's per-kernel "which params + value→kernel mapping" into a shared `bst_eigfilter_controls`; refactor `view_atom_designer` to consume it; add generator fields (`KernelName/KernelParams/Threshold`) to the `atomgroup` template; add `bst_dynamics('AtomFromKernel', …)` that realises via `bst_eigenfilter('Atom')`, thresholds via `Levelset`, and returns a populated atom group. No GUI changes beyond the designer refactor; the panel GUI is Plan 2.

**Tech Stack:** MATLAB R2023b, Brainstorm dev fork. Engine: `bst_eigen('Axes')`, `bst_eigenfilter('Atom')`, `bst_dynamics` (Levelset/NewGroup/Save/Load), `bst_eigfilter_kernel`.

## Global Constraints

- MATLAB R2023b; Brainstorm dev fork; **no new dependencies**.
- Tests are MATLAB assertion scripts in `dev/tests/`, run in the **live Brainstorm session** by the controller (`rehash; run(...)`). Per repo policy, **implementer subagents do static checks only** (`check_matlab_code` + grep + commit); they write the test files but do NOT run live MATLAB. The controller runs the consolidated live pass.
- The shared control spec must produce, for every existing kernel, the **same kernel-param struct** the designer's `i_phys2kernel` produced before the refactor (regression — designer behavior unchanged).
- `bst_eigfilter_controls` lives in `toolbox/eigen/eigfilter/` but must NOT match `bst_eigfilter_design_*.m` (so the registry's `dir` discovery does not list it as a kernel).
- Atom = sparse reference: the realised field is never stored; the atom group carries `{KernelName, KernelParams, Threshold}` + the materialized Scout (`region`/`radius`) + Event (`times`).
- `lint` every edited `.m` (MCP `check_matlab_code`); pre-existing idiom warnings acceptable, no new structural ones.
- Commit after each task; end every message with: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`

---

### Task 1: `bst_eigfilter_controls` — shared per-kernel control spec

**Files:**
- Create: `toolbox/eigen/eigfilter/bst_eigfilter_controls.m`
- Test: `dev/tests/test_eigfilter_controls.m`

**Interfaces:**
- Produces:
  - `S = bst_eigfilter_controls('Sliders', kernelName, bounds)` → `1×3` struct array, each `struct('label','Freq (Hz)','lo',0,'hi',50,'def',10,'fmt','%.1f')` or `struct('label','', ...)` for a disabled row. `bounds = struct('scaleMinMM',a,'scaleMaxMM',b,'rateMinMM2',c,'rateMaxMM2',d)` supplies the spectrum-dependent ranges.
  - `kp = bst_eigfilter_controls('ToKernel', kernelName, vals, lmax)` → kernel param struct; `vals = [s1 s2 s3]` (the three sliders' values).
- Consumes: nothing (pure data + mapping).

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigfilter_controls.m`:

```matlab
% test_eigfilter_controls - Sliders config + ToKernel mapping match the designer's prior behavior
b = struct('scaleMinMM',7,'scaleMaxMM',95,'rateMinMM2',49,'rateMaxMM2',9025);
lmax = 40;
% --- Sliders: labels per kernel ---
g = bst_eigfilter_controls('Sliders','gabor',b);
assert(strcmp(g(1).label,'Scale (mm)') && strcmp(g(2).label,'Freq (Hz)') && strcmp(g(3).label,'BW (Hz)'), 'gabor rows');
r = bst_eigfilter_controls('Sliders','resonator',b);
assert(strcmp(r(1).label,'Freq (Hz)') && strcmp(r(2).label,'Q') && isempty(r(3).label), 'resonator rows');
d = bst_eigfilter_controls('Sliders','diffusion',b);
assert(strcmp(d(1).label,'Rate (mm^2/s)') && d(1).lo==49 && d(1).hi==9025, 'diffusion rate range from bounds');
% --- ToKernel: matches the designer's i_phys2kernel formulas ---
kp = bst_eigfilter_controls('ToKernel','gabor',[30 12 2],lmax);
assert(abs(kp.k0 - 2*pi/0.03)<1e-9 && kp.f0==12 && kp.sf==2 && kp.lmax==lmax, 'gabor ToKernel');
kp = bst_eigfilter_controls('ToKernel','wave',[0 4 0],lmax);
assert(abs(kp.alpha - 4*sqrt(lmax)/2)<1e-12, 'wave alpha');
kp = bst_eigfilter_controls('ToKernel','dampedwave',[0 4 0.5],lmax);
assert(abs(kp.beta - 1/0.5)<1e-12, 'dampedwave beta');
kp = bst_eigfilter_controls('ToKernel','diffusion',[400 0 0],lmax);
assert(abs(kp.tau - max((400/1e6)*lmax,eps))<1e-15, 'diffusion tau');
kp = bst_eigfilter_controls('ToKernel','resonator',[10 6 0],lmax);
assert(kp.f0==10 && kp.Q==6, 'resonator');
kp = bst_eigfilter_controls('ToKernel','stmatern',[50 1.5 0],lmax);
assert(abs(kp.kappa - 2*pi/0.05)<1e-9 && kp.nu==1.5, 'stmatern');
disp('OK');
```

- [ ] **Step 2: Run to verify it fails** — `Undefined function 'bst_eigfilter_controls'`.

- [ ] **Step 3: Create the file**

`toolbox/eigen/eigfilter/bst_eigfilter_controls.m`:

```matlab
function out = bst_eigfilter_controls(action, kernelName, varargin)
% BST_EIGFILTER_CONTROLS: shared per-kernel GUI control spec for the atom designer AND panel_bst_dynamics.
%   S  = bst_eigfilter_controls('Sliders', name, bounds)        -> 1x3 struct(label,unit,lo,hi,def,fmt)
%   kp = bst_eigfilter_controls('ToKernel', name, vals, lmax)   -> kernel param struct (vals=[s1 s2 s3])
% Single source of truth for which params a kernel exposes, their ranges, and how slider values map to
% kernel params. Spectrum-dependent ranges (Scale/Rate) come from bounds (scaleMinMM/scaleMaxMM/rate*).
% Authors: Diellor Basha, 2026
switch lower(action)
    case 'sliders',  out = i_sliders(lower(kernelName), varargin{1});
    case 'tokernel', out = i_tokernel(lower(kernelName), varargin{1}, varargin{2});
    otherwise, error('bst_eigfilter_controls: unknown action ''%s''.', action);
end
end

function S = i_sliders(k, b)
    e   = i_row('',[],[],[],'');
    scD = round((b.scaleMinMM + b.scaleMaxMM)/2);               % spectrum-derived defaults
    rtD = round(((b.scaleMinMM + b.scaleMaxMM)/2)^2);
    sc  = i_row('Scale (mm)', b.scaleMinMM, b.scaleMaxMM, scD, '%.0f');
    switch k
        case 'diffusion',   S = [i_row('Rate (mm^2/s)', b.rateMinMM2, b.rateMaxMM2, rtD, '%.0f'), e, e];
        case {'heat','mexhat','diffgauss','flat','ideal','inverse_heat','log','matern','power','tikhonov'}
            S = [sc, e, e];
        case 'wave',        S = [e, i_row('Speed (m/s)',0.1,10,1,'%.2g'), e];
        case 'kleingordon', S = [e, i_row('Speed (m/s)',0.1,10,1,'%.2g'), e];
        case 'dampedwave',  S = [e, i_row('Speed (m/s)',0.1,10,1,'%.2g'), i_row('Decay (s)',0.05,2,0.5,'%.2g')];
        case 'gabor',       S = [sc, i_row('Freq (Hz)',0,50,10,'%.1f'), i_row('BW (Hz)',0.5,20,2,'%.1f')];
        case 'travwave',    S = [i_row('Speed (m/s)',0.05,3,1,'%.2g'), i_row('RidgeW (Hz)',0.5,20,2,'%.1f'), e];
        case 'resonator',   S = [i_row('Freq (Hz)',0,50,10,'%.1f'), i_row('Q',1,30,6,'%.1f'), e];
        case 'stmatern',    S = [i_row('Corr (mm)', b.scaleMinMM, b.scaleMaxMM, scD, '%.0f'), i_row('nu',0.5,4,1.5,'%.1f'), e];
        otherwise,          S = [sc, e, e];
    end
end

function kp = i_tokernel(k, v, lmax)
    s1 = v(1); s2 = v(2); s3 = v(3);
    kp = struct('lmax', lmax);
    switch k
        case 'wave',        kp.alpha = s2 * sqrt(lmax) / 2;
        case 'kleingordon', kp.alpha = s2 * sqrt(lmax) / 2;  kp.mu = 0.1*lmax;
        case 'dampedwave',  kp.alpha = s2 * sqrt(lmax) / 2;  kp.beta = 1/max(s3,eps);
        case 'diffusion',   kp.tau   = max((s1/1e6) * lmax, eps);
        case 'heat',        lamS = (2*pi/(s1/1000))^2;  kp.t = log(2)/max(lamS,eps);
        case 'mexhat',      lamS = (2*pi/(s1/1000))^2;  kp.t = 1/max(lamS,eps);
        case 'gabor',       kp.k0 = 2*pi/max(s1/1000,eps);  kp.f0 = s2;  kp.sf = max(s3,eps);
        case 'travwave',    kp.c  = s1;  kp.width = max(s2,eps);
        case 'resonator',   kp.f0 = s1;  kp.Q = max(s2,eps);
        case 'stmatern',    kp.kappa = 2*pi/max(s1/1000,eps);  kp.nu = max(s2,eps);
    end
end

function r = i_row(label, lo, hi, def, fmt)
    r = struct('label',label, 'lo',lo, 'hi',hi, 'def',def, 'fmt',fmt);
end
```

- [ ] **Step 4: Run to verify it passes** — `OK`.

- [ ] **Step 5: Lint + commit**

`check_matlab_code` on the new file (no structural warnings).
```bash
git add toolbox/eigen/eigfilter/bst_eigfilter_controls.m dev/tests/test_eigfilter_controls.m
git commit -m "feat(eigfilter): bst_eigfilter_controls shared per-kernel control spec (Sliders + ToKernel)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Refactor `view_atom_designer` to consume `bst_eigfilter_controls`

**Files:**
- Modify: `toolbox/gui/view_atom_designer.m` (`i_config_sliders` body → call `bst_eigfilter_controls('Sliders', …)`; `i_phys2kernel` body → call `bst_eigfilter_controls('ToKernel', …)`)
- Test: `dev/tests/test_designer_uses_controls.m`

**Interfaces:**
- Consumes: `bst_eigfilter_controls('Sliders'/'ToKernel', …)` from Task 1.
- Produces: no signature change; the designer's `i_config_sliders`/`i_phys2kernel` now delegate.

- [ ] **Step 1: Write the failing test**

`dev/tests/test_designer_uses_controls.m` — asserts the designer file delegates (the regression that designer *behavior* is unchanged is covered by Task 1's ToKernel tests, which encode the exact prior formulas):

```matlab
% test_designer_uses_controls - the designer delegates its control spec/mapping to the shared module
src = fileread('/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox/gui/view_atom_designer.m');
assert(contains(src, "bst_eigfilter_controls('Sliders'"), 'i_config_sliders delegates to shared Sliders');
assert(contains(src, "bst_eigfilter_controls('ToKernel'"), 'i_phys2kernel delegates to shared ToKernel');
% and the per-kernel switch tables are gone from the designer (now centralized)
assert(~contains(src, "case 'travwave',    kp.c"), 'i_phys2kernel switch removed');
disp('OK');
```

- [ ] **Step 2: Run to verify it fails** — assertion fails (designer still has its own switches).

- [ ] **Step 3: Replace `i_config_sliders` body**

In `view_atom_designer.m`, replace the whole `i_config_sliders(k)` function (currently lines ~156-199) with:

```matlab
    function i_config_sliders(k)
        % Per-kernel slider config comes from the shared control spec (designer + panel single source).
        b = struct('scaleMinMM',scaleMinMM, 'scaleMaxMM',scaleMaxMM, 'rateMinMM2',rateMinMM2, 'rateMaxMM2',rateMaxMM2);
        S = bst_eigfilter_controls('Sliders', k, b);
        H = {hScale,hScaleL,hScaleV; hSpeed,hSpeedL,hSpeedV; hDecay,hDecayL,hDecayV};
        for ii = 1:3
            i_setrow(H{ii,1}, H{ii,2}, H{ii,3}, S(ii).label, S(ii).lo, S(ii).hi, S(ii).def, S(ii).fmt);
        end
    end
```

- [ ] **Step 4: Replace `i_phys2kernel` body**

Replace the whole `i_phys2kernel()` function (currently lines ~345-359) with:

```matlab
    function kp = i_phys2kernel()
        kp = bst_eigfilter_controls('ToKernel', kernel, [get(hScale,'Value'), get(hSpeed,'Value'), get(hDecay,'Value')], lmax);
    end
```

**Note on vestigial state:** after this delegation, `i_config_sliders` no longer reads `pScaleMM`/`pRate` (it uses the shared spec's `def`, which equals their init mid-range values — so the designer's on-switch defaults are unchanged), and `i_phys2kernel` no longer reads `pSpeed`/`pDecay`. Those four state vars (and the `pScaleMM`/`pRate` clamps in `OperatorChanged`) become set-but-unused. Leave them as-is this task (harmless; behavior unchanged) — `check_matlab_code` may emit "value assigned might be unused" for them, which is acceptable; a follow-up can remove them.

- [ ] **Step 5: Run to verify it passes** — `OK`. Lint `view_atom_designer.m` (the only new warnings should be the vestigial-state notes above; `i_setrow`/`i_setslider` still used).

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/view_atom_designer.m dev/tests/test_designer_uses_controls.m
git commit -m "refactor(gui): atom designer consumes shared bst_eigfilter_controls (DRY with future panel)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Atom-group generator fields

**Files:**
- Modify: `toolbox/db/db_template.m` (`atomgroup` case — add `KernelName`, `KernelParams`, `Threshold`)
- Test: `dev/tests/test_atomgroup_generator.m`

**Interfaces:**
- Produces: `db_template('atomgroup')` carries `KernelName` (''), `KernelParams` ([]), `Threshold` ([]). `bst_dynamics('NewGroup')` inherits them; `Save`/`Load` round-trip them (whole-struct save).

- [ ] **Step 1: Write the failing test**

`dev/tests/test_atomgroup_generator.m`:

```matlab
% test_atomgroup_generator - atom group carries the kernel generator, round-trips through Save/Load
G = bst_dynamics('NewGroup','t');
assert(isfield(G,'KernelName') && isfield(G,'KernelParams') && isfield(G,'Threshold'), 'generator fields present');
G.KernelName = 'gabor';  G.KernelParams = struct('f0',10,'k0',209,'sf',2);  G.Threshold = 0.5;
T = bst_dynamics('New','gen'); T.SurfaceFile = 'x'; T = bst_dynamics('AddGroup', T, G);
tmp = [tempname '.mat'];  bst_dynamics('Save', tmp, T);  T2 = bst_dynamics('Load', tmp);  delete(tmp);
g2 = T2.Groups(1);
assert(strcmp(g2.KernelName,'gabor') && g2.KernelParams.f0==10 && g2.Threshold==0.5, 'generator round-trips');
disp('OK');
```

- [ ] **Step 2: Run to verify it fails** — `isfield` assert fails (fields absent).

- [ ] **Step 3: Add the fields to the `atomgroup` template**

In `toolbox/db/db_template.m`, locate `case 'atomgroup'`, and add three fields immediately before the final `'SurfaceFile', '');` line of that struct:

```matlab
            'KernelName',   '', ...     % GENERATOR: eigfilter kernel name (atom = thresholded localized filter)
            'KernelParams', [], ...     % GENERATOR: kernel param struct (bst_eigfilter_controls('ToKernel'))
            'Threshold',    [], ...     % GENERATOR: level-set threshold (fraction of peak) used for Scout/Event
            'SurfaceFile', '');         % provenance: cortex (space source)
```

(i.e. insert the three lines, keeping the existing `'SurfaceFile', '');` as the struct's last field.)

- [ ] **Step 4: Run to verify it passes** — `OK`.

- [ ] **Step 5: Lint + commit**
```bash
git add toolbox/*/db_template.m dev/tests/test_atomgroup_generator.m
git commit -m "feat(dynamics): atomgroup carries the kernel generator (KernelName/KernelParams/Threshold)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: `bst_dynamics('AtomFromKernel', …)` — realise → threshold → store path

**Files:**
- Modify: `toolbox/dynamics/bst_dynamics.m` (add `AtomFromKernel` verb)
- Test: `dev/tests/test_atom_from_kernel.m`

**Interfaces:**
- Produces: `G = bst_dynamics('AtomFromKernel', ax, kernelName, kernelParams, seed, thr)` → a populated `atomgroup`: `vertices=seed`, `region={scoutVertices}`, `times=[tEvtStart; tEvtEnd]`, `type='extended'`, `SurfaceFile=ax.SurfaceFile`, plus the generator fields. `ax` from `bst_eigen('Axes')`; `thr` = fraction of peak (default 0.5).
- Consumes: `bst_eigenfilter('Atom')` (engine), `Levelset`/`NewGroup` (this module).

- [ ] **Step 1: Write the failing test (synthetic eigenbasis)**

`dev/tests/test_atom_from_kernel.m`:

```matlab
% test_atom_from_kernel - realise+threshold an atom into a populated group; threshold monotonicity
nV=60; K=20; nT=64; Fs=100;
[Q,~]=qr(reshape(cos(1:(nV*K)),nV,K),0); Phi=Q; Lam=(linspace(0,5,K)').^2; M=speye(nV);
ax=struct('nT',nT,'NFFT',nT,'Fs',Fs,'SurfaceFile','synthetic'); ax.Phi={Phi}; ax.Lambda={Lam}; ax.Mass={M};
ax.GlobalVertices={(1:nV)'}; ax.tlag=(0:nT-1)/Fs; ax.omega=(0:nT-1)*(Fs/nT);
kp = struct('lmax',max(Lam),'tau',0.3);
G = bst_dynamics('AtomFromKernel', ax, 'diffusion', kp, 13, 0.5);
assert(G.vertices==13 && strcmp(G.KernelName,'diffusion') && G.Threshold==0.5, 'generator fields');
assert(~isempty(G.region{1}) && numel(G.times)==2 && G.times(2)>=G.times(1), 'Scout + Event populated');
% threshold monotonicity: higher threshold -> Scout is a subset
G9 = bst_dynamics('AtomFromKernel', ax, 'diffusion', kp, 13, 0.9);
assert(all(ismember(G9.region{1}, G.region{1})), 'higher threshold -> subset Scout');
disp('OK');
```

- [ ] **Step 2: Run to verify it fails** — unknown command `AtomFromKernel`.

- [ ] **Step 3: Add the verb**

Add to `toolbox/dynamics/bst_dynamics.m` (a new local function, reached via `eval(macro_method)`):

```matlab
function G = AtomFromKernel(ax, kernelName, kernelParams, seed, thr) %#ok<DEFNU>
    % Realise an eigfilter atom from its generator and threshold it into a populated atom group:
    % Scout = spatial level set (region), Event = temporal level set (times). Atom = thresholded filter.
    if (nargin < 6) || isempty(thr), thr = 0.5; end
    [W, gv] = bst_eigenfilter('Atom', ax, kernelName, kernelParams, seed);
    LS = Levelset(W, gv, thr);
    G = NewGroup(sprintf('%s @vtx%d', kernelName, seed));
    G.vertices = seed;
    if isfield(ax,'SurfaceFile'), G.SurfaceFile = ax.SurfaceFile; end
    G.region   = {LS.scoutVertices(:)'};
    G.times    = [ax.tlag(LS.eventSamples(1)); ax.tlag(LS.eventSamples(end))];
    G.type     = 'extended';
    G.KernelName = kernelName;  G.KernelParams = kernelParams;  G.Threshold = thr;
end
```

- [ ] **Step 4: Run to verify it passes** — `OK`.

- [ ] **Step 5: Lint + commit**
```bash
git add toolbox/dynamics/bst_dynamics.m dev/tests/test_atom_from_kernel.m
git commit -m "feat(dynamics): AtomFromKernel realises+thresholds an atom into a Scout+Event group

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Done criteria

- `bst_eigfilter_controls('Sliders'/'ToKernel', …)` covers all kernels; ToKernel reproduces the designer's prior `i_phys2kernel` formulas (Task 1 tests).
- The designer delegates to it; no per-kernel switch left in `view_atom_designer.m`; designer behavior unchanged.
- `db_template('atomgroup')` carries `KernelName/KernelParams/Threshold`, round-tripping through Save/Load.
- `bst_dynamics('AtomFromKernel', …)` realises+thresholds an atom into a populated Scout+Event group with monotone thresholding.
- All `dev/tests/test_{eigfilter_controls,designer_uses_controls,atomgroup_generator,atom_from_kernel}.m` print `OK`; lint clean.

## Follow-on (Plan 2 — not this plan)

The panel Java/Swing GUI: replace the Navigator with filters-as-tools + contextual params (rendering `bst_eigfilter_controls('Sliders')` via `gui_component`), click-to-localize, live preview via `bst_eigenfilter('Atom')`, and Store via `bst_dynamics('AtomFromKernel')`; retire `bst_geodesic_tool` (the `heat` preset). This substrate is its foundation.
