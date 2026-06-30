# Joint-spectral (js) kernels Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add four electrophysiology-motivated joint-spectral kernels (`gabor`, `travwave`, `resonator`, `stmatern`) to the eigfilter library and wire their parameters into the atom designer.

**Architecture:** Each kernel is a self-contained `bst_eigfilter_design_*.m` returning `meta` or a `g(l,w)` handle; `w` is the 0..Fs Hz grid, folded internally to signed frequencies so the realised atom is real. The designer's three sliders are relabeled/re-ranged per kernel and read directly by `i_phys2kernel`. No engine changes — kernels ride the Phase-1 `bst_eigenfilter('Atom')` js branch.

**Tech Stack:** MATLAB R2023b, Brainstorm dev fork. Registry `bst_eigfilter_kernel` auto-discovers `bst_eigfilter_design_*.m`. Realiser `bst_eigenfilter('Atom', ax, name, params, seed)`; axes `bst_eigen('Axes', OPTIONS)`.

## Global Constraints

- MATLAB R2023b; Brainstorm dev fork; **no new dependencies**.
- Tests are MATLAB assertion scripts in `dev/tests/`, run in the **live Brainstorm session** via the brainstorm-dev MATLAB MCP (`rehash; run('<path>')`). A test passes iff it prints `OK` with no error.
- Each design file follows the registry convention: `function out = bst_eigfilter_design_NAME(params)`; on `params=='meta'` return the meta struct (`name, display, params, domain, separable, bandpass, priorAdmissible`); else parse params (default any missing, incl. `lmax`) and return `out = @(l,w) i_eval(...)` calling file-local functions.
- **Realness idiom** (local `i_signed` in each js design file): `w=double(w(:)'); Fs=numel(w)*(w(2)-w(1)); ws = w - Fs.*(w >= Fs/2);` — define `g` even/Hermitian in `ws`.
- `l` is a `[K×1]` column, `w` a `[1×N]` row; `g → [K×N]`. Spatial frequency `k=√λ` (rad/m); `mm`-calibration `k = 2π/(mm/1000)`.
- `lint` every edited `.m` (MCP `check_matlab_code`); pre-existing globals/`%#ok`/comma-idiom warnings acceptable, new structural ones are not.
- Commit after each task; end every message with: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
- Tests use the live protocol surface `sub-MTL0002/tess_cortex_pial_low.mat` where a real cortex is needed; synthetic `ax` otherwise.

---

### Task 1: `gabor` kernel — joint scale×frequency packet

**Files:**
- Create: `toolbox/eigen/eigfilter/bst_eigfilter_design_gabor.m`
- Test: `dev/tests/test_js_gabor.m`

**Interfaces:**
- Produces: `g = bst_eigfilter_design_gabor(struct('f0',Hz,'k0',rad/m,'sf',Hz,'lmax',λmax))` → `@(l,w)→[K×N]`; `meta.domain='js'`, `meta.bandpass=true`.
- Consumes: `bst_eigenfilter('Atom')`, `manifold_ft/ift` (existing).

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_js_gabor.m`:

```matlab
% test_js_gabor - gabor js atom: real, and temporal spectrum peaks at f0
nV=60; K=20; nT=100; Fs=100;
[Q,~]=qr(reshape(cos(1:(nV*K)),nV,K),0); Phi=Q; Lam=(linspace(0,5,K)').^2; M=speye(nV);
ax=struct('nT',nT,'NFFT',nT,'Fs',Fs); ax.Phi={Phi}; ax.Lambda={Lam}; ax.Mass={M};
ax.GlobalVertices={(1:nV)'}; ax.tlag=(0:nT-1)/Fs; ax.omega=(0:nT-1)*(Fs/nT);
f0=12; k0=2*pi/0.03; sf=2;                                   % 12 Hz, ~30 mm scale
[W,gv]=bst_eigenfilter('Atom', ax, 'gabor', struct('f0',f0,'k0',k0,'sf',sf,'lmax',max(Lam)), 13);
assert(isequal(size(W),[nV nT]), 'shape');
assert(all(isfinite(W(:))), 'finite');
% temporal spectrum of the most active vertex peaks near f0
[~,iv]=max(sum(W.^2,2)); P=abs(fft(W(iv,:))).^2; fb=(0:nT-1)*(Fs/nT);
half=2:floor(nT/2); [~,ip]=max(P(half)); fpk=fb(half(ip));
assert(abs(fpk-f0) <= Fs/nT*1.5, sprintf('peak %.1f near f0=%.1f', fpk, f0));
disp('OK');
```

- [ ] **Step 2: Run to verify it fails**

Run (MCP): `rehash; try, run('<repo>/dev/tests/test_js_gabor.m'); catch e, fprintf('FAIL: %s\n',e.message); end`
Expected: FAIL — `unknown kernel 'gabor'`.

- [ ] **Step 3: Create the kernel**

`toolbox/eigen/eigfilter/bst_eigfilter_design_gabor.m`:

```matlab
function out = bst_eigfilter_design_gabor(params)
% BST_EIGFILTER_DESIGN_GABOR: joint-spectral Gabor packet g(l,w) localised at spatial wavenumber k0
% (rad/m) and temporal frequency f0 (Hz, bandwidth sf). Hermitian in temporal frequency -> real atom.
% Authors: Diellor Basha, 2026
if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','gabor', 'display','Gabor packet (scale x freq)', ...
        'params', struct('f0',struct('default',10,'range',[0 Inf]), ...
                         'k0',struct('default',209,'range',[0 Inf]), ...
                         'sf',struct('default',2,'range',[0 Inf])), ...
        'domain','js', 'separable',true, 'bandpass',true, 'priorAdmissible',false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'f0') || isempty(params.f0); params.f0 = 10;  end
if ~isfield(params,'k0') || isempty(params.k0); params.k0 = 209; end
if ~isfield(params,'sf') || isempty(params.sf); params.sf = 2;   end
out = @(l,w) i_gabor(l, w, params.k0, max(params.k0/3,eps), params.f0, max(params.sf,eps));
end
function G = i_gabor(l, w, k0, sk, f0, sf)
    ws = i_signed(w);
    sl = exp(-((sqrt(double(l(:)))-k0).^2)/(2*sk^2));                          % [K x 1]
    gt = exp(-((ws-f0).^2)/(2*sf^2)) + exp(-((ws+f0).^2)/(2*sf^2));            % [1 x N] Hermitian
    G  = sl * gt;
end
function ws = i_signed(w)
    w = double(w(:)');  Fs = numel(w)*(w(2)-w(1));  ws = w - Fs.*(w >= Fs/2);
end
```

- [ ] **Step 4: Run to verify it passes**

Run: `rehash; run('<repo>/dev/tests/test_js_gabor.m')` → `OK`. Also assert real: add nothing — `W` is already real via the realiser's `real(ifft)`; the Hermitian `gt` makes `max|imag|` ~1e-16 before the `real()`.

- [ ] **Step 5: Lint + commit**

`check_matlab_code` on the new file (no structural warnings).
```bash
git add toolbox/eigen/eigfilter/bst_eigfilter_design_gabor.m dev/tests/test_js_gabor.m
git commit -m "feat(eigfilter): gabor js kernel (joint scale x frequency packet)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: `travwave` kernel — dispersive traveling wave

**Files:**
- Create: `toolbox/eigen/eigfilter/bst_eigfilter_design_travwave.m`
- Test: `dev/tests/test_js_travwave.m`

**Interfaces:**
- Produces: `g = bst_eigfilter_design_travwave(struct('c',m/s,'width',Hz,'lmax',λmax))`; `meta.domain='js'`, `separable=false`, `bandpass=true`. Ridge `fr(λ)=c·√λ/(2π)` Hz.

- [ ] **Step 1: Write the failing test**

`dev/tests/test_js_travwave.m`:

```matlab
% test_js_travwave - traveling wave: real; peak temporal frequency grows with sqrt(lambda)
nV=60; K=20; nT=100; Fs=100;
[Q,~]=qr(reshape(cos(1:(nV*K)),nV,K),0); Phi=Q; Lam=(linspace(0.2,6,K)').^2; M=speye(nV);
ax=struct('nT',nT,'NFFT',nT,'Fs',Fs); ax.Phi={Phi}; ax.Lambda={Lam}; ax.Mass={M};
ax.GlobalVertices={(1:nV)'}; ax.tlag=(0:nT-1)/Fs; ax.omega=(0:nT-1)*(Fs/nT);
c=4; width=2;
g = bst_eigfilter_design_travwave(struct('c',c,'width',width,'lmax',max(Lam)));
G = g(Lam, ax.omega);                                            % [K x N]
assert(all(isfinite(G(:))), 'finite');
% ridge: per-mode peak frequency should increase with sqrt(lambda)
fb=(0:nT-1)*(Fs/nT); half=1:floor(nT/2);
[~,iLo]=min(Lam); [~,iHi]=max(Lam);
[~,pLo]=max(abs(G(iLo,half))); [~,pHi]=max(abs(G(iHi,half)));
assert(fb(half(pHi)) > fb(half(pLo)), 'higher lambda -> higher ridge frequency');
% realised atom is real
W = bst_eigenfilter('Atom', ax, 'travwave', struct('c',c,'width',width,'lmax',max(Lam)), 13);
assert(isequal(size(W),[nV nT]) && all(isfinite(W(:))), 'atom shape/finite');
disp('OK');
```

- [ ] **Step 2: Run to verify it fails** — `unknown kernel 'travwave'`.

- [ ] **Step 3: Create the kernel**

`bst_eigfilter_design_travwave.m`:

```matlab
function out = bst_eigfilter_design_travwave(params)
% BST_EIGFILTER_DESIGN_TRAVWAVE: joint-spectral traveling wave -- energy on the dispersion ridge
% f = c*sqrt(lambda)/(2*pi) Hz (phase speed c m/s, ridge width Hz). Non-separable; even in freq -> real.
% Authors: Diellor Basha, 2026
if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','travwave', 'display','Traveling wave (speed c)', ...
        'params', struct('c',struct('default',1,'range',[0 Inf]), ...
                         'width',struct('default',2,'range',[0 Inf])), ...
        'domain','js', 'separable',false, 'bandpass',true, 'priorAdmissible',false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'c')     || isempty(params.c);     params.c = 1;     end
if ~isfield(params,'width') || isempty(params.width); params.width = 2; end
out = @(l,w) i_trav(l, w, params.c, max(params.width,eps));
end
function G = i_trav(l, w, c, sg)
    ws = i_signed(w);
    fr = c .* sqrt(double(l(:))) ./ (2*pi);                                   % [K x 1] ridge freq (Hz)
    G  = exp(-((abs(ws) - fr).^2) / (2*sg^2));                                % [K x N], even in ws
end
function ws = i_signed(w)
    w = double(w(:)');  Fs = numel(w)*(w(2)-w(1));  ws = w - Fs.*(w >= Fs/2);
end
```

- [ ] **Step 4: Run to verify it passes** — `OK`.

- [ ] **Step 5: Lint + commit**
```bash
git add toolbox/eigen/eigfilter/bst_eigfilter_design_travwave.m dev/tests/test_js_travwave.m
git commit -m "feat(eigfilter): travwave js kernel (dispersive traveling wave)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: `resonator` kernel — damped harmonic oscillator (Lorentzian)

**Files:**
- Create: `toolbox/eigen/eigfilter/bst_eigfilter_design_resonator.m`
- Test: `dev/tests/test_js_resonator.m`

**Interfaces:**
- Produces: `g = bst_eigfilter_design_resonator(struct('f0',Hz,'Q',q,'lmax',λmax))`; `meta.domain='js'`, `bandpass=true`. Atom = decaying oscillation at `f0`; higher `Q` decays slower.

- [ ] **Step 1: Write the failing test**

`dev/tests/test_js_resonator.m`:

```matlab
% test_js_resonator - real; oscillates at f0; higher Q -> longer-lasting envelope
nV=60; K=20; nT=200; Fs=100;
[Qr,~]=qr(reshape(cos(1:(nV*K)),nV,K),0); Phi=Qr; Lam=(linspace(0,5,K)').^2; M=speye(nV);
ax=struct('nT',nT,'NFFT',nT,'Fs',Fs); ax.Phi={Phi}; ax.Lambda={Lam}; ax.Mass={M};
ax.GlobalVertices={(1:nV)'}; ax.tlag=(0:nT-1)/Fs; ax.omega=(0:nT-1)*(Fs/nT);
f0=10;
Wlo = bst_eigenfilter('Atom', ax, 'resonator', struct('f0',f0,'Q',3, 'lmax',max(Lam)), 13);
Whi = bst_eigenfilter('Atom', ax, 'resonator', struct('f0',f0,'Q',12,'lmax',max(Lam)), 13);
assert(isequal(size(Wlo),[nV nT]) && all(isfinite(Wlo(:))), 'shape/finite');
[~,iv]=max(sum(Whi.^2,2)); fb=(0:nT-1)*(Fs/nT); half=2:floor(nT/2);
P=abs(fft(Whi(iv,:))).^2;
[~,ip]=max(P(half)); assert(abs(fb(half(ip))-f0) <= Fs/nT*2, 'oscillates at f0');
% envelope persistence: late-window energy fraction larger for high Q
lateE = @(W) sum(sum(W(:,round(0.6*nT):end).^2)) / max(sum(W(:).^2),eps);
assert(lateE(Whi) > lateE(Wlo), 'higher Q -> more late-window energy');
disp('OK');
```

- [ ] **Step 2: Run to verify it fails** — `unknown kernel 'resonator'`.

- [ ] **Step 3: Create the kernel**

`bst_eigfilter_design_resonator.m`:

```matlab
function out = bst_eigfilter_design_resonator(params)
% BST_EIGFILTER_DESIGN_RESONATOR: joint-spectral damped harmonic oscillator (Lorentzian), peak at f0 Hz
% with quality Q. Hermitian in signed frequency -> real decaying-oscillation atom. lambda-independent.
% Authors: Diellor Basha, 2026
if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','resonator', 'display','Resonator (f0, Q)', ...
        'params', struct('f0',struct('default',10,'range',[0 Inf]), ...
                         'Q', struct('default',6, 'range',[0 Inf])), ...
        'domain','js', 'separable',true, 'bandpass',true, 'priorAdmissible',false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'f0') || isempty(params.f0); params.f0 = 10; end
if ~isfield(params,'Q')  || isempty(params.Q);  params.Q  = 6;  end
out = @(l,w) i_reson(l, w, params.f0, max(params.Q,eps));
end
function G = i_reson(l, w, f0, Q)
    ws = i_signed(w);
    H  = (f0.^2) ./ (f0.^2 - ws.^2 + 1i.*ws.*(f0/Q));                         % [1 x N] Hermitian
    G  = ones(numel(l),1) * H;                                               % [K x N], lambda-independent
end
function ws = i_signed(w)
    w = double(w(:)');  Fs = numel(w)*(w(2)-w(1));  ws = w - Fs.*(w >= Fs/2);
end
```

- [ ] **Step 4: Run to verify it passes** — `OK`.

- [ ] **Step 5: Lint + commit**
```bash
git add toolbox/eigen/eigfilter/bst_eigfilter_design_resonator.m dev/tests/test_js_resonator.m
git commit -m "feat(eigfilter): resonator js kernel (Lorentzian, f0/Q)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: `stmatern` kernel — spatiotemporal 1/f (Whittle–Matérn)

**Files:**
- Create: `toolbox/eigen/eigfilter/bst_eigfilter_design_stmatern.m`
- Test: `dev/tests/test_js_stmatern.m`

**Interfaces:**
- Produces: `g = bst_eigfilter_design_stmatern(struct('kappa',rad/m,'nu',ν,'lmax',λmax))`; `meta.domain='js'`, `bandpass=false`, `priorAdmissible=true`. Space–time speed `v` fixed at `1.0` m/s internally.

- [ ] **Step 1: Write the failing test**

`dev/tests/test_js_stmatern.m`:

```matlab
% test_js_stmatern - real, low-pass: positive-ish, temporal PSD monotone-decreasing (1/f)
nV=60; K=20; nT=128; Fs=100;
[Q,~]=qr(reshape(cos(1:(nV*K)),nV,K),0); Phi=Q; Lam=(linspace(0,5,K)').^2; M=speye(nV);
ax=struct('nT',nT,'NFFT',nT,'Fs',Fs); ax.Phi={Phi}; ax.Lambda={Lam}; ax.Mass={M};
ax.GlobalVertices={(1:nV)'}; ax.tlag=(0:nT-1)/Fs; ax.omega=(0:nT-1)*(Fs/nT);
g = bst_eigfilter_design_stmatern(struct('kappa',2*pi/0.05,'nu',1.5,'lmax',max(Lam)));
G = g(Lam, ax.omega);
assert(all(isfinite(G(:))) && all(G(:)>=0), 'finite, nonneg (power spectrum)');
% per-mode temporal spectrum decreases over the lower half (1/f)
half=1:floor(nT/2); row=abs(G(1,half));
assert(row(1) >= row(end), 'temporal spectrum decays with frequency');
W = bst_eigenfilter('Atom', ax, 'stmatern', struct('kappa',2*pi/0.05,'nu',1.5,'lmax',max(Lam)), 13);
assert(isequal(size(W),[nV nT]) && all(isfinite(W(:))), 'atom shape/finite');
disp('OK');
```

- [ ] **Step 2: Run to verify it fails** — `unknown kernel 'stmatern'`.

- [ ] **Step 3: Create the kernel**

`bst_eigfilter_design_stmatern.m`:

```matlab
function out = bst_eigfilter_design_stmatern(params)
% BST_EIGFILTER_DESIGN_STMATERN: spatiotemporal Whittle-Matern spectral density g = (kappa^2 + lambda +
% (w/v)^2)^(-nu) -- the 1/f aperiodic background / a spatiotemporal prior. v (space-time speed, m/s) fixed.
% Authors: Diellor Basha, 2026
V_SPEED = 1.0;   % m/s, fixed space-time coupling
if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','stmatern', 'display','Spatiotemporal 1/f (Whittle-Matern)', ...
        'params', struct('kappa',struct('default',126,'range',[0 Inf]), ...
                         'nu',   struct('default',1.5,'range',[0 Inf])), ...
        'domain','js', 'separable',false, 'bandpass',false, 'priorAdmissible',true);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'kappa') || isempty(params.kappa); params.kappa = 126; end
if ~isfield(params,'nu')    || isempty(params.nu);    params.nu    = 1.5; end
out = @(l,w) i_stm(l, w, params.kappa, params.nu, V_SPEED);
end
function G = i_stm(l, w, kappa, nu, v)
    ws = i_signed(w);
    base = (kappa.^2 + double(l(:))) * ones(1,numel(ws)) + ones(numel(l),1) * (ws./v).^2;   % [K x N]
    G = base .^ (-nu);
end
function ws = i_signed(w)
    w = double(w(:)');  Fs = numel(w)*(w(2)-w(1));  ws = w - Fs.*(w >= Fs/2);
end
```

- [ ] **Step 4: Run to verify it passes** — `OK`.

- [ ] **Step 5: Lint + commit**
```bash
git add toolbox/eigen/eigfilter/bst_eigfilter_design_stmatern.m dev/tests/test_js_stmatern.m
git commit -m "feat(eigfilter): stmatern js kernel (spatiotemporal 1/f / Whittle-Matern)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Designer control wiring (per-kernel three-slider mapping)

**Files:**
- Modify: `toolbox/gui/view_atom_designer.m` — capture all three slider label handles; replace `i_spatial_mode`/`ApplySpatial`/`i_spatial_str` with per-kernel `i_config_sliders`; rewrite `i_phys2kernel` to read sliders directly and map the 4 new kernels; simplify `ParamChanged`/`SyncControls`/`Status`.
- Test: `dev/tests/test_designer_controls.m`

**Interfaces:**
- Consumes: the 4 kernels (Tasks 1–4) via `bst_eigfilter_kernel('info', name)`.
- Produces: `i_config_sliders(kernel)` (configures the 3 sliders) and `i_phys2kernel()` returning the right param struct per kernel from current slider values.

- [ ] **Step 1: Write the failing test**

`dev/tests/test_designer_controls.m` — exercises the pure mapping with a stand-in slider state (the designer's `i_phys2kernel` will be refactored to read named handles; here we test the mapping logic in isolation by replicating it):

```matlab
% test_designer_controls - the 4 js kernels are discoverable, in the dynamic group, and realisable
ks = bst_eigfilter_kernel('list');
for n = {'gabor','travwave','resonator','stmatern'}
    assert(any(strcmp(ks, n{1})), sprintf('%s registered', n{1}));
    m = bst_eigfilter_kernel('info', n{1});
    assert(isfield(m,'domain') && strcmpi(m.domain,'js'), sprintf('%s domain js', n{1}));
end
% bandpass flags (gabor/travwave/resonator = wavelets; stmatern = low-pass)
assert(bst_eigfilter_kernel('info','gabor').bandpass==1, 'gabor bandpass');
assert(bst_eigfilter_kernel('info','stmatern').bandpass==0, 'stmatern lowpass');
disp('OK');
```

(The live slider relabel/range behaviour is GUI state; it is verified by the manual reopen check in Step 6, not by a headless assert.)

- [ ] **Step 2: Run to verify it passes already for registration** — this asserts Tasks 1–4 shipped; run `rehash; run('<repo>/dev/tests/test_designer_controls.m')` → `OK`. (If a kernel is missing it FAILs.) This is the guard that the kernels are wired into the registry before the GUI work.

- [ ] **Step 3: Capture the slider label handles**

In `view_atom_designer.m`, change the Speed/Decay slider creation (currently `[hSpeed,hSpeedV] = ...` / `[hDecay,hDecayV] = ...`) to capture the label handle:

```matlab
    [hScale,hScaleV,hScaleL] = i_slider(hP,'Scale', ys(1), lw,sw,vw, scaleMinMM, scaleMaxMM, pScaleMM, @ParamChanged);
    [hSpeed,hSpeedV,hSpeedL] = i_slider(hP,'Speed', ys(2), lw,sw,vw, 0.1, 10,  pSpeed, @ParamChanged);
    [hDecay,hDecayV,hDecayL] = i_slider(hP,'Decay', ys(3), lw,sw,vw, 0.05, 2,  pDecay, @ParamChanged);
```

- [ ] **Step 4: Replace `ApplySpatial`/`i_spatial_mode`/`i_spatial_str` with `i_config_sliders`**

Delete the functions `ApplySpatial`, `i_spatial_mode` (local, bottom), and `i_spatial_str`. Replace the body of `SyncControls` and add `i_config_sliders` + `i_setrow`:

```matlab
    function SyncControls()
        i_config_sliders(kernel);
        Status();
    end
    function i_config_sliders(k)
        % Per-kernel labels/units/ranges for the three sliders; values read live by i_phys2kernel.
        % Each row: {label, lo, hi, default, fmt} or [] to disable.
        sc = {scaleMinMM, scaleMaxMM};
        switch lower(k)
            case 'diffusion'
                i_setrow(hScale,hScaleL,hScaleV,'Rate (mm^2/s)',rateMinMM2,rateMaxMM2,pRate,'%.0f');
                i_setrow(hSpeed,hSpeedL,hSpeedV,[],0,1,0,'');  i_setrow(hDecay,hDecayL,hDecayV,[],0,1,0,'');
            case {'heat','mexhat','diffgauss','flat','ideal','inverse_heat','log','matern','power','tikhonov'}
                i_setrow(hScale,hScaleL,hScaleV,'Scale (mm)',sc{1},sc{2},pScaleMM,'%.0f');
                i_setrow(hSpeed,hSpeedL,hSpeedV,[],0,1,0,'');  i_setrow(hDecay,hDecayL,hDecayV,[],0,1,0,'');
            case 'wave'
                i_setrow(hScale,hScaleL,hScaleV,[],0,1,0,'');
                i_setrow(hSpeed,hSpeedL,hSpeedV,'Speed (m/s)',0.1,10,pSpeed,'%.2g');
                i_setrow(hDecay,hDecayL,hDecayV,[],0,1,0,'');
            case 'kleingordon'
                i_setrow(hScale,hScaleL,hScaleV,[],0,1,0,'');
                i_setrow(hSpeed,hSpeedL,hSpeedV,'Speed (m/s)',0.1,10,pSpeed,'%.2g');
                i_setrow(hDecay,hDecayL,hDecayV,[],0,1,0,'');
            case 'dampedwave'
                i_setrow(hScale,hScaleL,hScaleV,[],0,1,0,'');
                i_setrow(hSpeed,hSpeedL,hSpeedV,'Speed (m/s)',0.1,10,pSpeed,'%.2g');
                i_setrow(hDecay,hDecayL,hDecayV,'Decay (s)',0.05,2,pDecay,'%.2g');
            case 'gabor'
                i_setrow(hScale,hScaleL,hScaleV,'Scale (mm)',sc{1},sc{2},pScaleMM,'%.0f');
                i_setrow(hSpeed,hSpeedL,hSpeedV,'Freq (Hz)',0,50,10,'%.1f');
                i_setrow(hDecay,hDecayL,hDecayV,'BW (Hz)',0.5,20,2,'%.1f');
            case 'travwave'
                i_setrow(hScale,hScaleL,hScaleV,'Speed (m/s)',0.05,3,1,'%.2g');
                i_setrow(hSpeed,hSpeedL,hSpeedV,'RidgeW (Hz)',0.5,20,2,'%.1f');
                i_setrow(hDecay,hDecayL,hDecayV,[],0,1,0,'');
            case 'resonator'
                i_setrow(hScale,hScaleL,hScaleV,'Freq (Hz)',0,50,10,'%.1f');
                i_setrow(hSpeed,hSpeedL,hSpeedV,'Q',1,30,6,'%.1f');
                i_setrow(hDecay,hDecayL,hDecayV,[],0,1,0,'');
            case 'stmatern'
                i_setrow(hScale,hScaleL,hScaleV,'Corr (mm)',sc{1},sc{2},pScaleMM,'%.0f');
                i_setrow(hSpeed,hSpeedL,hSpeedV,'nu',0.5,4,1.5,'%.1f');
                i_setrow(hDecay,hDecayL,hDecayV,[],0,1,0,'');
            otherwise
                i_setrow(hScale,hScaleL,hScaleV,'Scale (mm)',sc{1},sc{2},pScaleMM,'%.0f');
                i_setrow(hSpeed,hSpeedL,hSpeedV,[],0,1,0,'');  i_setrow(hDecay,hDecayL,hDecayV,[],0,1,0,'');
        end
    end
    function i_setrow(hs, hl, hv, label, lo, hi, val, fmt)
        if isempty(label)
            set([hs hv], 'Enable','off');  set(hl,'String','');  set(hv,'String','--');
        else
            set(hl,'String',label);  i_setslider(hs, lo, hi, val);
            set([hs hv],'Enable','on');  set(hv,'String',sprintf(fmt, val));
        end
    end
```

- [ ] **Step 5: Rewrite `ParamChanged`, `i_phys2kernel`, `Status`**

`ParamChanged` (read each enabled slider's value into its readout, then regen):

```matlab
    function ParamChanged(~,~)
        i_refresh_readout(hScale,hScaleV);  i_refresh_readout(hSpeed,hSpeedV);  i_refresh_readout(hDecay,hDecayV);
        Regen();
    end
    function i_refresh_readout(hs, hv)
        if strcmpi(get(hs,'Enable'),'on'), set(hv,'String',sprintf('%.3g', get(hs,'Value'))); end
    end
```

`i_phys2kernel` (read slider values by their per-kernel role; `s1/s2/s3` are sliders 1/2/3):

```matlab
    function kp = i_phys2kernel()
        kp = struct('lmax', lmax);  s1=get(hScale,'Value'); s2=get(hSpeed,'Value'); s3=get(hDecay,'Value');
        switch lower(kernel)
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
```

`Status` (drop `i_spatial_str`; show the active kernel + first slider label/value):

```matlab
    function Status()
        s1 = get(hScale,'Value');  l1 = get(hScaleL,'String');
        set(hLabel,'String',sprintf('[%s] %s @ vtx %d | %s=%.3g | norm: %s | arrows=time', ...
            state, kernel, seedVtx, l1, s1, normMode));
    end
```

Also remove the now-unused `i_en2b` call from `Status` (handled above) and leave `i_en`, `i_setslider` (still used). Remove the `pSpeed`/`pDecay`/`pScaleMM`/`pRate` *reads* that referenced the deleted morphing, but KEEP those four as init defaults (they seed `i_config_sliders`); `OperatorChanged`'s clamps on `pScaleMM`/`pRate` stay (they bound the Scale/Rate defaults to the new spectrum).

- [ ] **Step 6: Verify — registration test + live reopen**

Run `rehash; run('<repo>/dev/tests/test_designer_controls.m')` → `OK`.
Manual (live GUI) check to record in the commit message: reopen the designer, select each js kernel, confirm the sliders relabel (gabor→Scale/Freq/BW; resonator→Freq/Q; travwave→Speed/RidgeW; stmatern→Corr/ν) and the atom updates. (Headless GUI launch is not asserted; the realiser path is covered by Tasks 1–4.)

- [ ] **Step 7: Lint + commit**

`check_matlab_code` on `view_atom_designer.m` (no new structural warnings).
```bash
git add toolbox/gui/view_atom_designer.m dev/tests/test_designer_controls.m
git commit -m "feat(gui): wire js kernels (gabor/travwave/resonator/stmatern) into the atom designer

Per-kernel 3-slider mapping (i_config_sliders + i_phys2kernel); Freq/Q/speed/nu controls.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Done criteria

- `bst_eigfilter_kernel('list')` includes `gabor, travwave, resonator, stmatern`, all `domain='js'`.
- Each realises a **real**, finite atom via `bst_eigenfilter('Atom')` and passes its defining-property test (this is the first real exercise of the js branch).
- The designer shows them in the dynamic group with relabeled sliders; `i_phys2kernel` maps slider values → kernel params.
- All `dev/tests/test_js_*.m` + `test_designer_controls.m` print `OK`; lint clean.
