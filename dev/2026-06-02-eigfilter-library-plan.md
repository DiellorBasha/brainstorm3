# bst_eigfilter Library Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `bst_eigfilter`, a GSP-style library of named analytic spectral kernels `g(λ)` for LBO eigenmodes, and unify the two existing `g(λ)` switches (`bst_eigenmode_prior`, `bst_eigenmodes_filter_gain`) onto it.

**Architecture:** One factory file per kernel under `toolbox/math/eigfilter/` (`bst_eigfilter_design_<name>.m`, each returns a function handle `@(λ)…`, GSP-faithful), a separate registry `bst_eigfilter_kernel.m` (string→factory + `list`/`info`, auto-discovering), plus `bst_eigfilter_evaluate.m` and `bst_eigfilter_compose.m`. The prior and filter-gain functions delegate to the library; their existing pure tests are kept as numeric parity guards.

**Tech Stack:** MATLAB, Brainstorm toolbox conventions, MATLAB MCP for running tests.

---

## Background references (read before starting)

- Spec: `dev/2026-06-02-eigfilter-library-design.md`
- Existing kernels to port / unify:
  - `toolbox/math/bst_eigenmode_prior.m` — prior `R = g(λ)` for `{flat, power, log}` (DC swap, `log` mm-scaling, `max(R)=1` normalization all live here and STAY here).
  - `toolbox/math/bst_eigenmodes_filter_gain.m` — gain `h = g(λ)` for `{lowpass, highpass, bandpass, heat, inverse_heat, tikhonov, custom}`.
- Parity-guard tests (must pass UNCHANGED after refactor):
  - `dev/tests/test_eigenmode_prior_pure.m`
  - `dev/tests/test_eigenmodes_filter_gain_pure.m`
- Test idiom: function script; `addpath(repoRoot)`; `if ~brainstorm('status'); brainstorm nogui; end`; `assert`; final `disp('ALL TESTS PASSED')`.

**Running a test (MATLAB MCP):** load `mcp__plugin_brainstorm-dev_MATLAB__run_matlab_file` via ToolSearch, call it on the absolute `.m` path. Success prints `ALL TESTS PASSED`. A MATLAB+Brainstorm session is already running.

**Branch:** before Task 1, create `feature/eigfilter-library` off `development`:
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git checkout development && git checkout -b feature/eigfilter-library
```

---

## File structure

| File | Responsibility |
|------|----------------|
| `toolbox/math/eigfilter/bst_eigfilter_design_<name>.m` (10 files) | Pure analytic factory: returns `@(λ)…` handle (or cell-of-handles for vector scale); answers `'meta'` |
| `toolbox/math/eigfilter/bst_eigfilter_kernel.m` | Registry: `(name,params)`→factory, `'list'`, `'info'` (auto-discovers by scanning `bst_eigfilter_design_*.m`) |
| `toolbox/math/eigfilter/bst_eigfilter_evaluate.m` | `h = evaluate(g, lambdas)`; handles single handle and cell bank |
| `toolbox/math/eigfilter/bst_eigfilter_compose.m` | `g = compose(g1,g2,…)`; pointwise product (serial filtering) |
| `toolbox/math/eigfilter/Contents.m` | Module index |
| `toolbox/math/bst_eigenmode_prior.m` (modify) | Delegate raw shape to library; keep DC/scaling/normalize/admissibility |
| `toolbox/math/bst_eigenmodes_filter_gain.m` (modify) | Delegate analytic types; index masks unchanged |
| `dev/tests/test_eigfilter_pure.m` | Factory + registry + evaluate + compose tests |

A standard Brainstorm header (the `@====` license block + `Authors: Diellor Basha, 2026`) should top each new file; the plan shows the functional body. Match the header of `bst_eigenmodes_filter_gain.m`.

---

## Task 1: Core infrastructure + reference kernel (heat)

**Files:**
- Create: `toolbox/math/eigfilter/bst_eigfilter_design_heat.m`, `bst_eigfilter_evaluate.m`, `bst_eigfilter_compose.m`, `bst_eigfilter_kernel.m`
- Test: `dev/tests/test_eigfilter_pure.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigfilter_pure.m`:

```matlab
function test_eigfilter_pure
% Verify the bst_eigfilter library: factories return handles, vector scale -> bank,
% 'meta' specs, registry list/info/dispatch/error, evaluate (handle + bank), compose.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

lam = [0; 1; 4; 9; 16; 25];

% ---- heat factory: handle, analytic value, t->0 ~ 1 ----
g = bst_eigfilter_design_heat(struct('t', 0.1));
assert(isa(g, 'function_handle'), 'factory must return a handle.');
assert(max(abs(g(lam) - exp(-0.1*lam))) < 1e-12, 'heat value wrong.');
g0 = bst_eigfilter_design_heat(struct('t', 1e-12));
assert(max(abs(g0(lam) - 1)) < 1e-6, 'heat t->0 must be ~1.');

% ---- vector scale -> cell bank ----
gb = bst_eigfilter_design_heat(struct('t', [0.1 0.5 1.0]));
assert(iscell(gb) && numel(gb) == 3, 'vector scale must return a 3-cell bank.');
assert(max(abs(gb{2}(lam) - exp(-0.5*lam))) < 1e-12, 'bank element wrong.');

% ---- meta ----
m = bst_eigfilter_design_heat('meta');
assert(strcmp(m.name,'heat') && isfield(m,'params') && isfield(m,'priorAdmissible'), 'meta malformed.');

% ---- evaluate: handle -> [K x 1], bank -> [K x Nf] ----
h1 = bst_eigfilter_evaluate(g, lam);
assert(isequal(size(h1), [numel(lam) 1]), 'evaluate handle shape.');
hB = bst_eigfilter_evaluate(gb, lam);
assert(isequal(size(hB), [numel(lam) 3]), 'evaluate bank shape.');

% ---- compose: product of two handles ----
g2 = bst_eigfilter_design_heat(struct('t', 0.2));
gc = bst_eigfilter_compose(g, g2);
assert(max(abs(gc(lam) - g(lam).*g2(lam))) < 1e-12, 'compose must be pointwise product.');

% ---- registry: dispatch, list, info, unknown error ----
gr = bst_eigfilter_kernel('heat', struct('t', 0.1));
assert(max(abs(gr(lam) - exp(-0.1*lam))) < 1e-12, 'registry dispatch wrong.');
names = bst_eigfilter_kernel('list');
assert(any(strcmp(names, 'heat')), 'list must include heat.');
mi = bst_eigfilter_kernel('info', 'heat');
assert(strcmp(mi.name, 'heat'), 'info must return heat meta.');
threw = false;
try, bst_eigfilter_kernel('no_such_kernel'); catch, threw = true; end
assert(threw, 'unknown kernel must error.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test, confirm it FAILS**

Run `dev/tests/test_eigfilter_pure.m`. Expected: FAIL (`Undefined function 'bst_eigfilter_design_heat'`).

- [ ] **Step 3: Create the four infrastructure files**

`toolbox/math/eigfilter/bst_eigfilter_design_heat.m` (body):
```matlab
function out = bst_eigfilter_design_heat(params)
% BST_EIGFILTER_DESIGN_HEAT: Heat / diffusion low-pass kernel g(l) = exp(-t*l).
% USAGE:  g = bst_eigfilter_design_heat(struct('t',0.01))   -> handle (cell if t is a vector)
%         m = bst_eigfilter_design_heat('meta')              -> metadata
% Optional params.lmax normalizes the spectrum: exp(-t*l/lmax).
if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','heat', 'display','Heat / diffusion (low-pass)', ...
        'params', struct('t', struct('default',0.01,'range',[0 Inf])), ...
        'bandpass', false, 'priorAdmissible', true);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'t') || isempty(params.t); params.t = 0.01; end
if ~isfield(params,'lmax'); params.lmax = []; end
sc = 1; if ~isempty(params.lmax) && params.lmax > 0; sc = 1/params.lmax; end
t = params.t;
if numel(t) > 1
    out = cell(numel(t),1);
    for ii = 1:numel(t); out{ii} = i_make(t(ii), sc); end
else
    out = i_make(t, sc);
end
end
function g = i_make(t, sc)
if t < 0; error('bst_eigfilter_design_heat: t must be >= 0 (got %g).', t); end
g = @(l) exp(-t * sc * double(l(:)));
end
```

`toolbox/math/eigfilter/bst_eigfilter_evaluate.m`:
```matlab
function h = bst_eigfilter_evaluate(g, lambdas)
% BST_EIGFILTER_EVALUATE: Evaluate a kernel handle (or cell bank) on eigenvalues.
%   h = [K x 1] for a single handle; h = [K x Nf] for a cell array of Nf handles.
lambdas = double(lambdas(:));
if iscell(g)
    h = zeros(numel(lambdas), numel(g));
    for ii = 1:numel(g); h(:,ii) = g{ii}(lambdas); end
else
    h = g(lambdas); h = h(:);
end
end
```

`toolbox/math/eigfilter/bst_eigfilter_compose.m`:
```matlab
function g = bst_eigfilter_compose(varargin)
% BST_EIGFILTER_COMPOSE: Serial composition of kernels = pointwise product of gains.
%   g = bst_eigfilter_compose(g1, g2, ...)  ->  @(l) g1(l).*g2(l).*...
gs = varargin;
for ii = 1:numel(gs)
    if ~isa(gs{ii}, 'function_handle')
        error('bst_eigfilter_compose: all arguments must be function handles.');
    end
end
g = @(l) i_prod(gs, l);
end
function y = i_prod(gs, l)
l = double(l(:));
y = ones(size(l));
for ii = 1:numel(gs); y = y .* gs{ii}(l); end
end
```

`toolbox/math/eigfilter/bst_eigfilter_kernel.m`:
```matlab
function out = bst_eigfilter_kernel(name, params)
% BST_EIGFILTER_KERNEL: Registry mapping a kernel name string to its factory.
% USAGE:  g     = bst_eigfilter_kernel(name, params)   -> handle from bst_eigfilter_design_<name>
%         names = bst_eigfilter_kernel('list')          -> available kernel names
%         meta  = bst_eigfilter_kernel('info', name)    -> metadata for one kernel
if nargin < 1; error('bst_eigfilter_kernel: name required.'); end
switch lower(name)
    case 'list'
        d = dir(fullfile(fileparts(mfilename('fullpath')), 'bst_eigfilter_design_*.m'));
        pre = numel('bst_eigfilter_design_');
        out = cellfun(@(nm) nm(pre+1:end-2), {d.name}, 'UniformOutput', false);
        return;
    case 'info'
        out = feval(['bst_eigfilter_design_' lower(params)], 'meta');
        return;
end
if nargin < 2; params = struct(); end
fn = ['bst_eigfilter_design_' lower(name)];
if ~exist(fn, 'file')
    error('bst_eigfilter_kernel: unknown kernel ''%s''.', name);
end
out = feval(fn, params);
end
```

- [ ] **Step 4: Run test, confirm it PASSES**

Run `dev/tests/test_eigfilter_pure.m`. Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Lint** `/lint-matlab` on the four new files (no new errors).

- [ ] **Step 6: Commit**
```bash
git add toolbox/math/eigfilter/bst_eigfilter_design_heat.m toolbox/math/eigfilter/bst_eigfilter_evaluate.m toolbox/math/eigfilter/bst_eigfilter_compose.m toolbox/math/eigfilter/bst_eigfilter_kernel.m dev/tests/test_eigfilter_pure.m
git commit -m "feat(eigfilter): core library (registry, evaluate, compose) + heat kernel

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Port the existing kernels (flat, power, log, inverse_heat, tikhonov, ideal)

**Files:** Create `bst_eigfilter_design_{flat,power,log,inverse_heat,tikhonov,ideal}.m`; extend `dev/tests/test_eigfilter_pure.m`.

- [ ] **Step 1: Extend the test (write failing assertions first)**

Append before `disp('ALL TESTS PASSED')` in `dev/tests/test_eigfilter_pure.m`:
```matlab
lp = [1; 2; 4; 8; 16];   % strictly positive for power/log
% flat
gf = bst_eigfilter_design_flat();
assert(all(bst_eigfilter_evaluate(gf, lp) == 1), 'flat must be ones.');
% power: lambda^-alpha, decreasing
gp = bst_eigfilter_design_power(struct('alpha',1));
assert(max(abs(gp(lp) - lp.^(-1))) < 1e-12, 'power value wrong.');
% log: -log(lambda), needs lambda in (0,1) for positivity (prior scales it)
gl = bst_eigfilter_design_log();
assert(max(abs(gl([0.1;0.5]) - (-log([0.1;0.5])))) < 1e-12, 'log value wrong.');
% inverse_heat: clamped
gih = bst_eigfilter_design_inverse_heat(struct('t',0.1,'maxgain',5));
assert(max(gih(lp)) <= 5 + 1e-12, 'inverse_heat must clamp at maxgain.');
% tikhonov
gt = bst_eigfilter_design_tikhonov(struct('beta',2));
assert(max(abs(gt(lp) - 1./(1+2*lp))) < 1e-12, 'tikhonov value wrong.');
% ideal band [2 8] inclusive
gi = bst_eigfilter_design_ideal(struct('band',[2 8]));
assert(isequal(gi(lp), double(lp>=2 & lp<=8)), 'ideal mask wrong.');
% registry sees all of them
nm = bst_eigfilter_kernel('list');
for k = {'flat','power','log','inverse_heat','tikhonov','ideal'}
    assert(any(strcmp(nm, k{1})), sprintf('list missing %s.', k{1}));
end
```

- [ ] **Step 2: Run test, confirm new assertions FAIL** (`Undefined function 'bst_eigfilter_design_flat'`).

- [ ] **Step 3: Create the six factories**

`bst_eigfilter_design_flat.m`:
```matlab
function out = bst_eigfilter_design_flat(params) %#ok<INUSD>
% BST_EIGFILTER_DESIGN_FLAT: g(l) = 1 (no spectral preference).
if nargin >= 1 && ischar(params) && strcmpi(params,'meta')
    out = struct('name','flat','display','Flat (no spectral prior)', ...
        'params', struct(), 'bandpass', false, 'priorAdmissible', true);
    return;
end
out = @(l) ones(numel(double(l(:))), 1);
end
```

`bst_eigfilter_design_power.m`:
```matlab
function out = bst_eigfilter_design_power(params)
% BST_EIGFILTER_DESIGN_POWER: g(l) = l.^(-alpha). Caller must pass l > 0.
if nargin >= 1 && ischar(params) && strcmpi(params,'meta')
    out = struct('name','power','display','Power law l^-alpha', ...
        'params', struct('alpha', struct('default',1,'range',[0 Inf])), ...
        'bandpass', false, 'priorAdmissible', true);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'alpha') || isempty(params.alpha); params.alpha = 1; end
a = params.alpha;
out = @(l) double(l(:)).^(-a);
end
```

`bst_eigfilter_design_log.m`:
```matlab
function out = bst_eigfilter_design_log(params) %#ok<INUSD>
% BST_EIGFILTER_DESIGN_LOG: raw g(l) = -log(l). Meaningful only for l in (0,1);
% the eigenmode prior performs the millimetre rescaling and DC handling.
if nargin >= 1 && ischar(params) && strcmpi(params,'meta')
    out = struct('name','log','display','Logarithmic (GBF 2026)', ...
        'params', struct(), 'bandpass', false, 'priorAdmissible', true);
    return;
end
out = @(l) -log(double(l(:)));
end
```

`bst_eigfilter_design_inverse_heat.m`:
```matlab
function out = bst_eigfilter_design_inverse_heat(params)
% BST_EIGFILTER_DESIGN_INVERSE_HEAT: sharpening g(l) = min(exp(+t*l), maxgain).
if nargin >= 1 && ischar(params) && strcmpi(params,'meta')
    out = struct('name','inverse_heat','display','Inverse heat (sharpening)', ...
        'params', struct('t', struct('default',0.01,'range',[0 Inf]), ...
                         'maxgain', struct('default',10,'range',[1 Inf])), ...
        'bandpass', false, 'priorAdmissible', true);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'t') || isempty(params.t); params.t = 0.01; end
if ~isfield(params,'maxgain') || isempty(params.maxgain); params.maxgain = 10; end
if params.t < 0; error('bst_eigfilter_design_inverse_heat: t must be >= 0.'); end
t = params.t; mg = params.maxgain;
out = @(l) min(exp(t * double(l(:))), mg);
end
```

`bst_eigfilter_design_tikhonov.m`:
```matlab
function out = bst_eigfilter_design_tikhonov(params)
% BST_EIGFILTER_DESIGN_TIKHONOV: low-pass g(l) = 1/(1+beta*l).
if nargin >= 1 && ischar(params) && strcmpi(params,'meta')
    out = struct('name','tikhonov','display','Tikhonov / membrane (low-pass)', ...
        'params', struct('beta', struct('default',1,'range',[0 Inf])), ...
        'bandpass', false, 'priorAdmissible', true);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'beta') || isempty(params.beta); params.beta = 1; end
if params.beta < 0; error('bst_eigfilter_design_tikhonov: beta must be >= 0.'); end
b = params.beta;
out = @(l) 1 ./ (1 + b * double(l(:)));
end
```

`bst_eigfilter_design_ideal.m`:
```matlab
function out = bst_eigfilter_design_ideal(params)
% BST_EIGFILTER_DESIGN_IDEAL: brick-wall indicator g(l) = 1[lo <= l <= hi].
% params.band = [lo hi] (default [0 Inf]).
if nargin >= 1 && ischar(params) && strcmpi(params,'meta')
    out = struct('name','ideal','display','Ideal (brick-wall) band', ...
        'params', struct('band', struct('default',[0 Inf])), ...
        'bandpass', true, 'priorAdmissible', true);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'band') || isempty(params.band); params.band = [0 Inf]; end
lo = params.band(1); hi = params.band(2);
out = @(l) double(double(l(:)) >= lo & double(l(:)) <= hi);
end
```

- [ ] **Step 4: Run test, confirm PASS** (`ALL TESTS PASSED`).

- [ ] **Step 5: Lint** the six files.

- [ ] **Step 6: Commit**
```bash
git add toolbox/math/eigfilter/bst_eigfilter_design_flat.m toolbox/math/eigfilter/bst_eigfilter_design_power.m toolbox/math/eigfilter/bst_eigfilter_design_log.m toolbox/math/eigfilter/bst_eigfilter_design_inverse_heat.m toolbox/math/eigfilter/bst_eigfilter_design_tikhonov.m toolbox/math/eigfilter/bst_eigfilter_design_ideal.m dev/tests/test_eigfilter_pure.m
git commit -m "feat(eigfilter): port flat/power/log/inverse_heat/tikhonov/ideal kernels

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: New kernels (matern, mexhat, dog)

**Files:** Create `bst_eigfilter_design_{matern,mexhat,dog}.m`; extend `dev/tests/test_eigfilter_pure.m`.

- [ ] **Step 1: Extend the test (failing first)**

Append before `disp('ALL TESTS PASSED')`:
```matlab
lpp = [1; 2; 4; 8; 16];
% matern: (kappa^2 + l)^-nu, decreasing, positive
gm = bst_eigfilter_design_matern(struct('kappa',1,'nu',1.5));
assert(max(abs(gm(lpp) - (1 + lpp).^(-1.5))) < 1e-12, 'matern value wrong.');
assert(all(diff(gm(lpp)) < 0), 'matern must be decreasing.');
% mexhat band-pass: zero at 0, peak interior, decay; vector t -> bank
gh = bst_eigfilter_design_mexhat(struct('t',0.1));
assert(abs(gh(0)) < 1e-12, 'mexhat must be 0 at lambda=0.');
ll = (0:0.01:50)'; v = gh(ll);
assert(v(1) < max(v) && v(end) < max(v), 'mexhat must peak in the interior.');
ghb = bst_eigfilter_design_mexhat(struct('t',[0.05 0.1 0.2]));
assert(iscell(ghb) && numel(ghb)==3, 'mexhat vector t must return a bank.');
mh = bst_eigfilter_design_mexhat('meta');
assert(mh.priorAdmissible == false, 'mexhat must be flagged not prior-admissible.');
% dog band-pass: non-negative for t1<t2, zero at 0
gd = bst_eigfilter_design_dog(struct('t1',0.1,'t2',0.4));
assert(abs(gd(0)) < 1e-12, 'dog must be 0 at lambda=0.');
assert(all(gd(lpp) >= -1e-12), 'dog must be non-negative for t1<t2.');
md = bst_eigfilter_design_dog('meta');
assert(md.priorAdmissible == false, 'dog must be flagged not prior-admissible.');
```

- [ ] **Step 2: Run, confirm FAIL.**

- [ ] **Step 3: Create the three factories**

`bst_eigfilter_design_matern.m`:
```matlab
function out = bst_eigfilter_design_matern(params)
% BST_EIGFILTER_DESIGN_MATERN: SPDE/Matern prior g(l) = (kappa^2 + l)^(-nu).
if nargin >= 1 && ischar(params) && strcmpi(params,'meta')
    out = struct('name','matern','display','Matern / SPDE (low-pass)', ...
        'params', struct('kappa', struct('default',1,'range',[0 Inf]), ...
                         'nu',    struct('default',1,'range',[0 Inf])), ...
        'bandpass', false, 'priorAdmissible', true);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'kappa') || isempty(params.kappa); params.kappa = 1; end
if ~isfield(params,'nu') || isempty(params.nu); params.nu = 1; end
k = params.kappa; nu = params.nu;
out = @(l) (k^2 + double(l(:))).^(-nu);
end
```

`bst_eigfilter_design_mexhat.m`:
```matlab
function out = bst_eigfilter_design_mexhat(params)
% BST_EIGFILTER_DESIGN_MEXHAT: Mexican-hat band-pass g(l) = (t*l).*exp(-t*l).
% Vector t returns a cell-array filterbank. Zero at l=0 (not prior-admissible).
if nargin >= 1 && ischar(params) && strcmpi(params,'meta')
    out = struct('name','mexhat','display','Mexican hat (band-pass)', ...
        'params', struct('t', struct('default',0.01,'range',[0 Inf])), ...
        'bandpass', true, 'priorAdmissible', false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'t') || isempty(params.t); params.t = 0.01; end
t = params.t;
if numel(t) > 1
    out = cell(numel(t),1);
    for ii = 1:numel(t); out{ii} = i_mh(t(ii)); end
else
    out = i_mh(t);
end
end
function g = i_mh(t)
if t < 0; error('bst_eigfilter_design_mexhat: t must be >= 0.'); end
g = @(l) (t * double(l(:))) .* exp(-t * double(l(:)));
end
```

`bst_eigfilter_design_dog.m`:
```matlab
function out = bst_eigfilter_design_dog(params)
% BST_EIGFILTER_DESIGN_DOG: difference-of-Gaussians band-pass
% g(l) = exp(-t1*l) - exp(-t2*l). Non-negative when t1 < t2; zero at l=0.
if nargin >= 1 && ischar(params) && strcmpi(params,'meta')
    out = struct('name','dog','display','Difference of Gaussians (band-pass)', ...
        'params', struct('t1', struct('default',0.01,'range',[0 Inf]), ...
                         't2', struct('default',0.04,'range',[0 Inf])), ...
        'bandpass', true, 'priorAdmissible', false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'t1') || isempty(params.t1); params.t1 = 0.01; end
if ~isfield(params,'t2') || isempty(params.t2); params.t2 = 0.04; end
if params.t1 >= params.t2
    error('bst_eigfilter_design_dog: require t1 < t2 (got t1=%g, t2=%g).', params.t1, params.t2);
end
t1 = params.t1; t2 = params.t2;
out = @(l) exp(-t1 * double(l(:))) - exp(-t2 * double(l(:)));
end
```

- [ ] **Step 4: Run, confirm PASS.**
- [ ] **Step 5: Lint** the three files.
- [ ] **Step 6: Commit**
```bash
git add toolbox/math/eigfilter/bst_eigfilter_design_matern.m toolbox/math/eigfilter/bst_eigfilter_design_mexhat.m toolbox/math/eigfilter/bst_eigfilter_design_dog.m dev/tests/test_eigfilter_pure.m
git commit -m "feat(eigfilter): add matern, mexhat, dog kernels

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Unify bst_eigenmode_prior onto the library

**Files:** Modify `toolbox/math/bst_eigenmode_prior.m`. Parity guard: `dev/tests/test_eigenmode_prior_pure.m` (unchanged). New admissibility assertions added to that test.

- [ ] **Step 1: Read the current file** `toolbox/math/bst_eigenmode_prior.m`. It computes (after DC swap): `flat`→`ones`; `power`→`lamK.^(-alpha)` with `lamK = max(lam(1:K), max(lam(1:K))*1e-12)`; `log`→`-log(lamMM)` with `lamMM = min(max(lam(1:K)*1e-6, eps), 1-1e-12)`; then `R = R/max(R)` (except flat returns early). Preserve every one of these numerics.

- [ ] **Step 2: Refactor the `switch` to delegate the raw shape, keep scaling/normalize/admissibility**

Replace the `switch lower(priorType) … end` block (the `flat`/`power`/`log` cases) with:
```matlab
switch lower(priorType)
    case 'flat'
        g = bst_eigfilter_kernel('flat');
        R = bst_eigfilter_evaluate(g, lam(1:K));      % ones
        % (no normalization needed; already all ones)
    case 'power'
        lamK = lam(1:K);
        lamK = max(lamK, max(lamK) * 1e-12);
        g = bst_eigfilter_kernel('power', struct('alpha', alpha));
        R = bst_eigfilter_evaluate(g, lamK);
    case 'log'
        M2MM2 = 1e-6;
        lamMM = lam(1:K) * M2MM2;
        lamMM = min(max(lamMM, eps), 1 - 1e-12);
        g = bst_eigfilter_kernel('log');
        R = bst_eigfilter_evaluate(g, lamMM);
    otherwise
        error('bst_eigenmode_prior:UnknownPrior', 'Unknown priorType: %s', priorType);
end

% Admissibility: a prior is a covariance (finite, non-negative). Guard against
% a caller selecting an analysis-only kernel (e.g. mexhat/dog) as a prior.
minfo = bst_eigfilter_kernel('info', priorType);
if (isfield(minfo,'priorAdmissible') && ~minfo.priorAdmissible) ...
        || any(~isfinite(R)) || any(R < 0)
    error('bst_eigenmode_prior:Inadmissible', ...
        'Kernel ''%s'' is not admissible as a prior (covariance must be finite and >= 0).', priorType);
end

% Normalize to max 1 (flat is already ones)
R = R / max(R);
```
Keep the existing function signature, the leading `lambdas = double(lambdas(:))`, the `K = max(1, min(K, nAvail))` clamp, the input validation, and the DC-mode swap exactly as they are. Only the per-type body + the new admissibility/normalize tail change. (Note: `flat` previously returned early; now it flows through `R = R/max(R)` which leaves ones unchanged — numerically identical.)

Update the header comment to note delegation to `bst_eigfilter_kernel`.

- [ ] **Step 3: Run the parity guard** `dev/tests/test_eigenmode_prior_pure.m`. Expected: `ALL TESTS PASSED` (flat/power/log numerics unchanged: shapes, positivity, `max==1`, monotonicity, DC swap, gentle log rolloff).

- [ ] **Step 4: Add an admissibility assertion to the prior test**

Append before `disp('ALL TESTS PASSED')` in `dev/tests/test_eigenmode_prior_pure.m`:
```matlab
% Analysis-only kernels must be rejected as priors
threw = false;
try, bst_eigenmode_prior(lambdas, K, 'mexhat', 0); catch, threw = true; end
assert(threw, 'mexhat must be rejected as a prior.');
threw = false;
try, bst_eigenmode_prior(lambdas, K, 'dog', 0); catch, threw = true; end
assert(threw, 'dog must be rejected as a prior.');
```
Run again. Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Lint** `bst_eigenmode_prior.m`.

- [ ] **Step 6: Commit**
```bash
git add toolbox/math/bst_eigenmode_prior.m dev/tests/test_eigenmode_prior_pure.m
git commit -m "refactor(eigfilter): bst_eigenmode_prior delegates to the kernel library

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Unify bst_eigenmodes_filter_gain onto the library

**Files:** Modify `toolbox/math/bst_eigenmodes_filter_gain.m`. Parity guard: `dev/tests/test_eigenmodes_filter_gain_pure.m` (unchanged).

- [ ] **Step 1: Delegate only the analytic types; keep masks and custom in place**

In the `switch lower(FilterType)` block, replace the `heat`, `inverse_heat`, and `tikhonov` case bodies (keeping their existing parameter-validation guards) with delegations; leave `lowpass`/`highpass`/`bandpass`/`custom`/`otherwise` exactly as they are.

```matlab
    case 'heat'
        if DiffusionTime <= 0
            error('DiffusionTime must be positive (got %g).', DiffusionTime);
        end
        g = bst_eigfilter_kernel('heat', struct('t', DiffusionTime));
        h = bst_eigfilter_evaluate(g, lambdas);
    case 'inverse_heat'
        if DiffusionTime <= 0
            error('DiffusionTime must be positive (got %g).', DiffusionTime);
        end
        g = bst_eigfilter_kernel('inverse_heat', struct('t', DiffusionTime, 'maxgain', MaxGain));
        h = bst_eigfilter_evaluate(g, lambdas);
    case 'tikhonov'
        if RegBeta < 0
            error('RegBeta must be non-negative (got %g).', RegBeta);
        end
        g = bst_eigfilter_kernel('tikhonov', struct('beta', RegBeta));
        h = bst_eigfilter_evaluate(g, lambdas);
```
These produce `exp(-DiffusionTime*lambdas)`, `min(exp(DiffusionTime*lambdas), MaxGain)`, and `1./(1+RegBeta*lambdas)` respectively — numerically identical to the current inline code. Update the SEE ALSO / header to note delegation.

- [ ] **Step 2: Run the parity guard** `dev/tests/test_eigenmodes_filter_gain_pure.m`. Expected: `ALL TESTS PASSED` (all seven types incl. the index masks and `custom` unchanged).

- [ ] **Step 3: Lint** `bst_eigenmodes_filter_gain.m`.

- [ ] **Step 4: Commit**
```bash
git add toolbox/math/bst_eigenmodes_filter_gain.m
git commit -m "refactor(eigfilter): bst_eigenmodes_filter_gain delegates analytic kernels

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: Contents.m, path verification, integration check

**Files:** Create `toolbox/math/eigfilter/Contents.m`. Verify path + run full suite.

- [ ] **Step 1: Verify the new subdirectory is on the Brainstorm path**

Run via MATLAB MCP (throwaway, not committed):
```matlab
if ~brainstorm('status'); brainstorm nogui; end
disp(exist('bst_eigfilter_kernel','file'));   % expect 2
disp(numel(bst_eigfilter_kernel('list')));    % expect 10
```
Expected: `2` and `10`. If `exist` returns `0`, Brainstorm did not add `toolbox/math/eigfilter/` to the path. In that case inspect how `toolbox/math` subfolders are added at startup (search `brainstorm.m` / `bst_set` path setup for `genpath`/`addpath` of the toolbox tree) and report the finding; the fix (e.g. ensuring the folder is included) becomes an added step. Do NOT hand-wave — confirm `list` returns all 10 kernel names: flat, power, log, heat, inverse_heat, tikhonov, ideal, matern, mexhat, dog.

- [ ] **Step 2: Create `toolbox/math/eigfilter/Contents.m`**
```matlab
% EIGFILTER  Spectral-filter kernels for Laplace-Beltrami eigenmodes (bst_eigfilter)
%
% Filters along the eigenvalue (spatial-frequency) axis, the spatial-spectral
% counterpart of bst_freqfilter (temporal frequency). Each kernel is an analytic
% factory returning a function handle g = @(lambda) ...; a vector-valued scale
% parameter returns a cell-array filterbank.
%
% Registry / apply:
%   bst_eigfilter_kernel    - name string -> factory; 'list' / 'info'
%   bst_eigfilter_evaluate  - evaluate a handle (or bank) on eigenvalues
%   bst_eigfilter_compose   - serial composition (pointwise product of gains)
%
% Kernel factories (bst_eigfilter_design_<name>):
%   flat          - g = 1
%   power         - g = lambda^-alpha
%   log           - g = -log(lambda)            (GBF 2026; prior rescales)
%   heat          - g = exp(-t*lambda)          (low-pass / diffusion)
%   inverse_heat  - g = min(exp(+t*lambda), maxgain)   (sharpening)
%   tikhonov      - g = 1/(1+beta*lambda)       (low-pass / membrane)
%   ideal         - g = 1[lo <= lambda <= hi]   (brick-wall band)
%   matern        - g = (kappa^2 + lambda)^-nu  (SPDE / Gaussian field)
%   mexhat        - g = (t*lambda).*exp(-t*lambda)     (band-pass; bank-capable)
%   dog           - g = exp(-t1*lambda) - exp(-t2*lambda)  (band-pass)
%
% Consumed by: bst_eigenmode_prior (prior R), bst_eigenmodes_filter_gain (analysis).
%
% Authors: Diellor Basha, 2026
```

- [ ] **Step 3: Run the whole eigfilter-related suite** via MATLAB MCP, each must print `ALL TESTS PASSED`:
  - `dev/tests/test_eigfilter_pure.m`
  - `dev/tests/test_eigenmode_prior_pure.m`
  - `dev/tests/test_eigenmodes_filter_gain_pure.m`
  - `dev/tests/test_eigenmodes_filter_pure.m` (consumer of filter_gain — regression)
  - `dev/tests/test_inverse_eigenmodes_pure.m` (consumer of the prior — regression)

- [ ] **Step 4: Commit**
```bash
git add toolbox/math/eigfilter/Contents.m
git commit -m "docs(eigfilter): module Contents.m index

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final verification

- [ ] All five tests in Task 6 Step 3 print `ALL TESTS PASSED`.
- [ ] `bst_eigfilter_kernel('list')` returns all 10 kernel names.
- [ ] `git diff development -- toolbox/math/bst_eigenmodes_filter.m toolbox/process/functions/process_eigenmodes_coeffsfilter.m` is **empty** (consumers untouched; behavior preserved purely via the parity guards).
- [ ] Use `superpowers:finishing-a-development-branch`.
