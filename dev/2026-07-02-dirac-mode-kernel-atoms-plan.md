# Dirac mode-kernel atoms Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** For a Dirac-dSPM source, filter atoms directly on the inverse's eigen-coefficients (`c = ImagingKernelMode·recordings`) in the inverse's own `DiracEigenFile` basis, instead of reconstructing the field and re-projecting onto a 60-mode canonical basis.

**Architecture:** Reuse the proven `view_eigen_timeseries`/`view_eigenmode_spectrum` machinery: get `c` from the mode kernel (free), filter `c_filt=g(λ)·c`, and reconstruct only for the views. The atom's axes become the inverse's `DiracEigenFile` basis. Scalar/non-dSPM paths fall back to today's reconstruct-then-project. The mode-kernel path yields the amplitude current (direction-identical to dSPM); a per-vertex dSPM scale toggles the cortex magnitude.

**Tech Stack:** MATLAB, Brainstorm; `bst_memory('GetRecordingsValues'/'GetResultsValues')`, `in_bst_results`/`in_bst_eigen`/`in_bst_operator`, `bst_eigenwavelet('Scalogram'/'JTVAtoms')`, `bst_dirac`, the dimensional-atoms decode/quiver display.

## Global Constraints

- **Mode coefficients:** `c = double(R.ImagingKernelMode) * double(bst_memory('GetRecordingsValues', D.srcDS, R.GoodChannel, iWin))` → `[nMode×nWin]`, split per hemi by `R.ModeHemisphere`, each hemi ordered ascending `λ` (matches `view_eigen_timeseries`: `ordL=find(hemi==1)` sorted by `Eigenvalues`, then `ordR`).
- **Session basis (Dirac-dSPM):** `E=in_bst_eigen(R.DiracEigenFile)` → `Phi{1×2}[4Vh×Kh]`, `GlobalVertices{1×2}`; `O=in_bst_operator(E.OperatorFile)` → `Mass{1×2}`. `ax.Lambda{h}` = `Eigenvalues` rows for hemi `h`, ascending. Interleaved quaternion `[w,x,y,z]`; physical current = imag rows `2:4:end/3:4:end/4:4:end`, `w=0`.
- **Verified invariant (anchor):** `SIR(v)·reconstruct(Φ, ImagingKernelMode)_v == ImagingKernel_v` per vertex; per-vertex `cos(reconstruct, ImagingKernel) = 1`. `SIR(v)=‖ImagingKernel_v‖/‖reconstruct_v‖` is data-independent.
- **Measure:** filtering is in the amplitude mode domain. Cortex magnitude = amplitude (default) or `SIR(v)·amplitude` (dSPM toggle). Direction/quiver/differential/sensor are measure-invariant.
- **Fallback (byte-unchanged):** scalar operators (`Laplace-Beltrami`/`LB-Connectome`) and non-Dirac-dSPM sources keep the current reconstruct-then-project path.
- **A Dirac-dSPM source** = `~isempty(R.ImagingKernelMode) && ~isempty(R.DiracEigenFile)`.
- No `clear` while Brainstorm runs (use `rehash`); prefer `matlab -batch` for headless tests (Apple-silicon GUI drops `GlobalData`); reuse — don't re-roll the recordings loader or the basis load.
- Commit trailers (development): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>` + `Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL`.

---

### Task 1: Detect Dirac-dSPM + mode coefficients

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m` (add locals `i_is_dirac_dspm`, `i_mode_coeffs`); Test `dev/test_mode_coeffs.m`.

**Interfaces:**
- Produces `tf = panel_bst_dynamics('i_is_dirac_dspm', D)` — true if the linked source has `ImagingKernelMode` + `DiracEigenFile`.
- Produces `[cCell, ax_meta] = panel_bst_dynamics('i_mode_coeffs', st, D, iWin)` — `cCell{1×2}` per-hemi `[Kh×nWin]` (asc-λ order); `ax_meta` = struct `.Eigenvalues{1×2}` (asc-λ λ per hemi), `.DiracEigenFile`. Cached on `getappdata(0,'DynamicsModeCoeffCache')` keyed `DiracEigenFile|srcResult|iWin`.

- [ ] **Step 1: Failing test** (`dev/test_mode_coeffs.m`)

```matlab
function tests = test_mode_coeffs
tests = functiontests(localfunctions);
end
function tc = i_launch(~)
    % helper: launch a Dirac-dSPM Dynamics session, return st + D
end
function test_shape_and_split(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st), 'launch a Dirac dSPM Dynamics session first');
    D  = getappdata(st.hFig,'DynamicsOverlay');
    assertTrue(t, panel_bst_dynamics('i_is_dirac_dspm', D));
    iWin = 1:5;
    [cCell, meta] = panel_bst_dynamics('i_mode_coeffs', st, D, iWin);
    verifyEqual(t, numel(cCell), 2);
    verifyEqual(t, size(cCell{1},2), 5);
    % per-hemi mode counts match ModeHemisphere
    src = panel_bst_dynamics('i_src_resultfile', D);
    R = in_bst_results(src,0,'ModeHemisphere');
    verifyEqual(t, size(cCell{1},1), sum(R.ModeHemisphere==1));
    verifyEqual(t, size(cCell{2},1), sum(R.ModeHemisphere==2));
    % each hemi ascending lambda
    verifyTrue(t, issorted(meta.Eigenvalues{1}));
end
```

- [ ] **Step 2: Run (headless)** `matlab -batch "addpath(genpath('/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox')); brainstorm nogui; <launch a Dirac dSPM Dynamics session>; runtests('dev/test_mode_coeffs.m')"`. Expected FAIL (`i_mode_coeffs` undefined). (The controller supplies the launched session for the live gate; the implementer transcribes the code.)

- [ ] **Step 3: Implement**

```matlab
% True when the linked source is a Dirac-dSPM inverse (carries the mode kernel + its eigenbasis file).
function tf = i_is_dirac_dspm(D) %#ok<DEFNU>
    tf = false;  src = i_src_resultfile(D);  if isempty(src), return; end
    try
        R = in_bst_results(src, 0, 'ImagingKernelMode', 'DiracEigenFile');
        tf = isfield(R,'ImagingKernelMode') && ~isempty(R.ImagingKernelMode) ...
          && isfield(R,'DiracEigenFile')   && ~isempty(R.DiracEigenFile);
    catch %#ok<CTCH>
    end
end

% Mode coefficients c = ImagingKernelMode * recordings(GoodChannel, iWin), split per hemisphere and
% ordered ascending-lambda (mirrors view_eigen_timeseries). Free vs projecting the reconstructed field.
function [cCell, meta] = i_mode_coeffs(st, D, iWin) %#ok<DEFNU>
    cCell = {[],[]};  meta = struct('Eigenvalues',{{[],[]}},'DiracEigenFile','');
    src = i_src_resultfile(D);  if isempty(src), return; end
    R = in_bst_results(src, 0, 'ImagingKernelMode','Eigenvalues','ModeHemisphere','GoodChannel','DiracEigenFile');
    if isempty(R.ImagingKernelMode), return; end
    key = sprintf('%s|%s|%d-%d', R.DiracEigenFile, src, iWin(1), iWin(end));
    M = getappdata(0,'DynamicsModeCoeffCache');
    if ~isempty(M) && isstruct(M) && strcmp(M.key,key), cCell = M.cCell; meta = M.meta; return; end
    gc = R.GoodChannel;  if isempty(gc), gc = 1:size(R.ImagingKernelMode,2); end
    d  = double(bst_memory('GetRecordingsValues', D.srcDS, gc, iWin));   % [nGoodChan x nWin]
    cAll = double(R.ImagingKernelMode) * d;                             % [nMode x nWin]
    lam  = double(R.Eigenvalues(:));  hemi = double(R.ModeHemisphere(:));
    for h = 1:2
        ord = find(hemi==h);  [ls, s] = sort(lam(ord),'ascend');  ord = ord(s);
        cCell{h} = cAll(ord,:);  meta.Eigenvalues{h} = ls;
    end
    meta.DiracEigenFile = R.DiracEigenFile;
    setappdata(0,'DynamicsModeCoeffCache', struct('key',key,'cCell',{cCell},'meta',meta));
end
```

- [ ] **Step 4: Run** — Expected PASS (2 hemis, right counts, asc λ).
- [ ] **Step 5: Commit** — `feat(dynamics): Dirac-dSPM mode coefficients c=KernelMode*recordings (free projection)`.

---

### Task 2: Session Dirac basis in `i_atom_axes`

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m` (`i_atom_axes`, line ~506); Test `dev/test_session_basis.m`.

**Interfaces:**
- Consumes `i_is_dirac_dspm` (Task 1).
- Produces: `ax = i_atom_axes(st,'Dirac')` returns the inverse's basis when the source is Dirac-dSPM: `ax.Phi/GlobalVertices` from `in_bst_eigen(R.DiracEigenFile)`, `ax.Mass/Operator` from `in_bst_operator(E.OperatorFile)`, `ax.Lambda{h}` = ascending `Eigenvalues` for hemi `h`, `ax.SurfaceFile`, plus the existing `ax.nT/tlag/omega/NFFT` fields. Canonical `bst_eigen('Axes')` otherwise.

- [ ] **Step 1: Failing test** (`dev/test_session_basis.m`)

```matlab
function tests = test_session_basis
tests = functiontests(localfunctions);
end
function test_uses_inverse_basis(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st));
    ax = panel_bst_dynamics('i_atom_axes', st, 'Dirac');
    D  = getappdata(st.hFig,'DynamicsOverlay');
    src = panel_bst_dynamics('i_src_resultfile', D);
    R = in_bst_results(src,0,'Eigenvalues','ModeHemisphere','DiracEigenFile');
    verifyEqual(t, size(ax.Lambda{1},1), sum(R.ModeHemisphere==1));   % inverse's mode count, not 60
    verifyGreaterThan(t, size(ax.Lambda{1},1), 60);
    verifyEqual(t, size(ax.Phi{1},1), 4*numel(ax.GlobalVertices{1})); % quaternion layout
    verifyTrue(t, issorted(ax.Lambda{1}));
end
```

- [ ] **Step 2: Run** — Expected FAIL (returns 60-mode canonical ax).

- [ ] **Step 3: Implement** — at the top of `i_atom_axes(st, variant)`, before the canonical `bst_eigen('Axes')` build, add the Dirac-dSPM branch:

```matlab
    % Dirac-dSPM: use the INVERSE's own eigenbasis (DiracEigenFile) so atoms filter the mode kernel's
    % coefficients directly, at the inverse's full mode count (not a fresh 60-mode canonical basis).
    if strcmp(variant, 'Dirac')
        D = getappdata(st.hFig, 'DynamicsOverlay');
        if ~isempty(D) && i_is_dirac_dspm(D)
            surf = i_atom_surface(st);
            key  = ['dspm|' variant '|' surf];
            Mc = getappdata(0,'DynamicsAtomAx');  if isempty(Mc)||~isa(Mc,'containers.Map'), Mc=containers.Map('KeyType','char','ValueType','any'); end
            if isKey(Mc,key), ax = Mc(key); return; end
            src = i_src_resultfile(D);
            R = in_bst_results(src, 0, 'DiracEigenFile','Eigenvalues','ModeHemisphere','SurfaceFile');
            E = in_bst_eigen(R.DiracEigenFile);
            O = in_bst_operator(E.OperatorFile);
            if isfield(E,'Phi') && ~isempty(E.Phi) && isfield(O,'Mass') && numel(O.Mass)==2
                lam = double(R.Eigenvalues(:));  hemi = double(R.ModeHemisphere(:));
                ax = struct('Variant','Dirac','SurfaceFile',R.SurfaceFile, ...
                            'Phi',{E.Phi},'GlobalVertices',{E.GlobalVertices},'Mass',{O.Mass},'Operator',O);
                for h = 1:2
                    ord = find(hemi==h);  ls = sort(lam(ord),'ascend');  ax.Lambda{h} = ls;
                end
                Fs = 100; try, tv = bst_memory('GetTimeVector', D.srcDS, D.srcResult); if numel(tv)>1, Fs=1/median(diff(tv)); end, catch, end %#ok<CTCH>
                nF = max(2, round(4*Fs));  ax.nT = nF;  ax.tlag = (0:nF-1)/Fs;  ax.omega = (0:nF-1)*(Fs/nF);  ax.NFFT = nF;
                Mc(key) = ax;  setappdata(0,'DynamicsAtomAx', Mc);
                return;
            end
            % basis missing/malformed -> fall through to the canonical path
        end
    end
```

(The existing canonical body — `ax = bst_eigen('Axes', struct(...,'nModes',60,...))` — remains as the fallback below.)

- [ ] **Step 4: Run** — Expected PASS (inverse's mode count, quaternion layout, asc λ).
- [ ] **Step 5: Commit** — `feat(dynamics): atom axes use the inverse DiracEigenFile basis for a Dirac-dSPM source`.

---

### Task 3: Per-vertex dSPM scale + the lossless-projection anchor

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m` (add `i_dspm_scale`); Test `dev/test_dspm_scale.m`.

**Interfaces:**
- Consumes `i_atom_axes` (Task 2).
- Produces `sir = panel_bst_dynamics('i_dspm_scale', st, D)` — `[nV×1]` per-vertex positive scale `‖ImagingKernel_v‖/‖reconstruct(ImagingKernelMode)_v‖`; `[]` for an amplitude-measure source (all ratios ≈1) or if `ImagingKernel` absent. Cached keyed `DiracEigenFile|srcResult`.

- [ ] **Step 1: Failing test** (`dev/test_dspm_scale.m`) — this test IS the anchor:

```matlab
function tests = test_dspm_scale
tests = functiontests(localfunctions);
end
function test_anchor_direction_and_scale(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st));
    D  = getappdata(st.hFig,'DynamicsOverlay');
    src = panel_bst_dynamics('i_src_resultfile', D);
    R = in_bst_results(src,0,'ImagingKernel','ImagingKernelMode','Eigenvalues','ModeHemisphere','DiracEigenFile');
    E = in_bst_eigen(R.DiracEigenFile);
    lam=double(R.Eigenvalues(:)); hemi=double(R.ModeHemisphere(:));
    nV=0; for h=1:2, nV=max(nV,max(E.GlobalVertices{h}(:))); end
    Km=double(R.ImagingKernelMode); Krec=zeros(3*nV,size(Km,2));
    for h=1:2
        ord=find(hemi==h);[~,s]=sort(lam(ord),'ascend');ord=ord(s);
        Ph=double(E.Phi{h}); gv=E.GlobalVertices{h}(:); Uf=Ph*Km(ord,:);
        Krec((gv-1)*3+1,:)=Uf(2:4:end,:);Krec((gv-1)*3+2,:)=Uf(3:4:end,:);Krec((gv-1)*3+3,:)=Uf(4:4:end,:);
    end
    Kv=double(R.ImagingKernel);
    % direction identical per vertex
    cosv=zeros(nV,1);
    for v=1:nV, a=Krec((v-1)*3+(1:3),:); b=Kv((v-1)*3+(1:3),:); na=norm(a,'fro'); nb=norm(b,'fro');
        if na>0&&nb>0, cosv(v)=sum(a(:).*b(:))/(na*nb); end, end
    verifyGreaterThan(t, median(cosv(cosv~=0)), 1-1e-6);
    % i_dspm_scale reproduces the per-vertex ratio: SIR .* recon == ImagingKernel
    sir = panel_bst_dynamics('i_dspm_scale', st, D);
    Kfix = Krec .* reshape(repmat(sir(:)',3,1),[],1);
    verifyLessThan(t, norm(Kfix(:)-Kv(:))/norm(Kv(:)), 1e-6);
end
```

- [ ] **Step 2: Run** — Expected FAIL (`i_dspm_scale` undefined).

- [ ] **Step 3: Implement**

```matlab
% Per-vertex dSPM/sLORETA scale SIR(v) = ||ImagingKernel_v|| / ||reconstruct(ImagingKernelMode)_v||.
% Data-independent (from the stored kernels); [] when the source is amplitude-measure (ratios ~1) or has
% no vertex kernel. Lets the cortex magnitude toggle between amplitude and the source's dSPM measure.
function sir = i_dspm_scale(st, D) %#ok<DEFNU>
    sir = [];  src = i_src_resultfile(D);  if isempty(src), return; end
    R = in_bst_results(src, 0, 'ImagingKernel','ImagingKernelMode','Eigenvalues','ModeHemisphere','DiracEigenFile');
    if isempty(R.ImagingKernel) || isempty(R.ImagingKernelMode), return; end
    key = sprintf('%s|%s|sir', R.DiracEigenFile, src);
    Mc = getappdata(0,'DynamicsDspmScale');
    if ~isempty(Mc) && isstruct(Mc) && strcmp(Mc.key,key), sir = Mc.sir; return; end
    E = in_bst_eigen(R.DiracEigenFile);  lam=double(R.Eigenvalues(:)); hemi=double(R.ModeHemisphere(:));
    nV=0; for h=1:2, nV=max(nV,max(E.GlobalVertices{h}(:))); end
    Km=double(R.ImagingKernelMode); Krec=zeros(3*nV,size(Km,2));
    for h=1:2
        ord=find(hemi==h);[~,s]=sort(lam(ord),'ascend');ord=ord(s);
        Ph=double(E.Phi{h}); gv=E.GlobalVertices{h}(:); Uf=Ph*Km(ord,:);
        Krec((gv-1)*3+1,:)=Uf(2:4:end,:);Krec((gv-1)*3+2,:)=Uf(3:4:end,:);Krec((gv-1)*3+3,:)=Uf(4:4:end,:);
    end
    Kv=double(R.ImagingKernel);  sir=ones(nV,1);
    for v=1:nV
        a=Krec((v-1)*3+(1:3),:); b=Kv((v-1)*3+(1:3),:);  na=norm(a,'fro');
        if na>0, sir(v)=norm(b,'fro')/na; end
    end
    if max(abs(sir-1)) < 1e-6, sir = []; end     % amplitude source: no renorm needed
    setappdata(0,'DynamicsDspmScale', struct('key',key,'sir',sir));
end
```

- [ ] **Step 4: Run** — Expected PASS (cos≈1; `SIR·recon==ImagingKernel` to 1e-6).
- [ ] **Step 5: Commit** — `feat(dynamics): recover per-vertex dSPM scale + lossless-projection anchor test`.

---

### Task 4: Reconstruct + sensor-forward from coefficients

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m` (add `i_dirac_recon`, `i_dirac_forward_modes`); Test `dev/test_dirac_recon.m`.

**Interfaces:**
- Consumes `i_atom_axes` (Task 2), the D-era `i_dirac_leadfield`.
- Produces:
  - `[V3, mag] = panel_bst_dynamics('i_dirac_recon', ax, cCell)` — reconstruct per hemi `Uf=Φ_h·cCell{h}`, extract imag → full-surface `V3 [nV×3×nT]` and per-vertex magnitude `mag [nV×nT]` (amplitude).
  - `Dfilt = panel_bst_dynamics('i_dirac_forward_modes', ax, Leig, cCell)` — stack `cCell` L-then-R and `Dfilt = Leig * cstack`; `[]` if `Leig` empty.

- [ ] **Step 1: Failing test** (`dev/test_dirac_recon.m`) — reconstruct unfiltered mode c and check it equals the amplitude field direction (per §4):

```matlab
function tests = test_dirac_recon
tests = functiontests(localfunctions);
end
function test_recon_matches_amplitude_direction(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st));
    D  = getappdata(st.hFig,'DynamicsOverlay');
    ax = panel_bst_dynamics('i_atom_axes', st, 'Dirac');
    iWin = 1:4;
    [cCell,~] = panel_bst_dynamics('i_mode_coeffs', st, D, iWin);
    [V3, mag] = panel_bst_dynamics('i_dirac_recon', ax, cCell);
    nV = size(V3,1);  verifyEqual(t, size(V3,3), 4);  verifyEqual(t, size(mag), [nV 4]);
    % direction equals GetResultsValues (dSPM) per vertex at frame 1
    Jgt = double(bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], iWin, 0)); % [3nV x nWin]
    v1 = reshape(V3(:,:,1)',[],1);                 % [3nV x 1] recon frame 1
    g1 = Jgt(:,1);
    cang = dot(v1,g1)/(norm(v1)*norm(g1)+eps);
    verifyGreaterThan(t, cang, 0.999);             % same overall direction (amplitude vs dSPM = per-vtx scale)
end
```

- [ ] **Step 2: Run** — Expected FAIL (`i_dirac_recon` undefined).

- [ ] **Step 3: Implement**

```matlab
% Reconstruct a per-hemi mode-coefficient set to the full-surface ambient 3-vector (imag quaternion slots)
% + per-vertex magnitude (amplitude current). cCell{h} = [Kh x nT].
function [V3, mag] = i_dirac_recon(ax, cCell) %#ok<DEFNU>
    nV=0; for h=1:numel(ax.GlobalVertices), nV=max(nV,max(ax.GlobalVertices{h}(:))); end
    nT = 0; for h=1:numel(cCell), if ~isempty(cCell{h}), nT=size(cCell{h},2); break; end, end
    V3 = zeros(nV,3,nT);
    for h=1:numel(ax.Phi)
        Ph=ax.Phi{h}; if isempty(Ph)||isempty(cCell{h}), continue; end
        gv=ax.GlobalVertices{h}(:);  Uf = Ph * cCell{h};                  % [4Vh x nT]
        V3(gv,1,:)=reshape(Uf(2:4:end,:),numel(gv),1,nT);
        V3(gv,2,:)=reshape(Uf(3:4:end,:),numel(gv),1,nT);
        V3(gv,3,:)=reshape(Uf(4:4:end,:),numel(gv),1,nT);
    end
    mag = squeeze(sqrt(sum(V3.^2,2)));  if nT==1, mag=mag(:); end
end

% Sensor forward from mode coefficients: stack cCell L-then-R and Dfilt = Leig * cstack. [] if no Leig.
function Dfilt = i_dirac_forward_modes(ax, Leig, cCell) %#ok<DEFNU>
    Dfilt = [];  if isempty(Leig), return; end
    cstack = [];  for h=1:numel(cCell), cstack = [cstack; cCell{h}]; end %#ok<AGROW>
    Dfilt = Leig * cstack;
end
```

- [ ] **Step 4: Run** — Expected PASS (shape + direction cos>0.999).
- [ ] **Step 5: Commit** — `feat(dynamics): reconstruct + sensor-forward Dirac atoms from mode coefficients`.

---

### Task 5: Wire the Dirac Apply branch to the mode kernel

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m` (`i_atom_apply` Dirac branch, ~line 913; `i_dirac_leadfield` to accept the session basis). Test: live (controller) — the Apply cortex + sensor path.

**Interfaces:** Consumes Tasks 1–4. The Dirac Apply branch, for a Dirac-dSPM source, sources `cCell` from `i_mode_coeffs`, filters `c_filt{h}=g(ax.Lambda{h}).*cCell{h}`, reconstructs cortex via `i_dirac_recon`, forwards sensors via `i_dirac_forward_modes`. Non-dSPM Dirac keeps the current `i_dirac_forward(J)` path.

- [ ] **Step 1:** In the `i_atom_apply` Dirac branch, after the static-kernel guard and `Leig = i_dirac_leadfield(st, ax)`, branch on `i_is_dirac_dspm(D)`:

```matlab
        if i_is_dirac_dspm(D)
            [cCell,~] = i_mode_coeffs(st, D, iWin);
            g = bst_eigfilter_kernel(kernel, kp);
            cf = cCell;  for h=1:numel(cf), if ~isempty(cf{h}), cf{h} = g(ax.Lambda{h}(:)) .* cf{h}; end, end
            [V3, mag] = i_dirac_recon(ax, cf);                 % amplitude current
            sir = i_dspm_scale(st, D);
            magShow = mag;
            if strcmpi(i_field(st,'atomMeasure','amplitude'),'dspm') && ~isempty(sir), magShow = mag .* sir(:); end
            nV = size(mag,1);
            if ~isempty(st.hFig) && ishandle(st.hFig)
                view_dynamics('SetFilteredField', st.hFig, magShow, (1:nV)', iWin, false);
                Vq = zeros(nV,3);  Vq(:,:) = V3(:,:,max(1,round(size(V3,3)/2)));  % mid-window frame for quivers
                setappdata(st.hFig,'QuiverVectorOverride', Vq);
                try, figure_3d('SetShowSourceVectors', st.hFig, D.iTess, 1); catch, end %#ok<CTCH>
            end
            Dfilt = i_dirac_forward_modes(ax, Leig, cf);
            % ... existing sensor-overlay block using Dfilt (unchanged) ...
            return;
        end
```

(Keep the existing non-dSPM Dirac path — `GetResultsValues` + `i_dirac_forward(ax,Leig,J,g)` — below this branch.)

- [ ] **Step 2:** Update `i_dirac_leadfield` so `L_eig`'s Tau/nModes come from the session `ax` (already the inverse basis for a dSPM source) and it asserts `CompHM.Eigenvalues` matches `[ax.Lambda{:}]` (reuse the existing D guard; the count is now the inverse's, e.g. 400/hemi).

- [ ] **Step 3: Live verify (controller, MCP):** Dirac-dSPM session → Apply ON → cortex filtered magnitude (amplitude) + quivers; sensor overlay draws; toggling the measure to dSPM rescales the cortex magnitude but not the quivers; scalar operator still uses the old path. Confirm `getappdata(hFig,'QuiverVectorOverride')` non-empty and the info string reports the sensor overlay. Screenshot.
- [ ] **Step 4: Commit** — `feat(dynamics): Dirac Apply filters the inverse mode coefficients (cortex+sensor)`.

---

### Task 6: Scalogram + Localize from the mode kernel

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m` (`i_apply_projection` Dirac-dSPM branch feeding `OnAnalyzeWindow`/`OnLocalizeBands`). Test `dev/test_scalogram_modes.m`.

**Interfaces:** Consumes Tasks 1–2. For a Dirac-dSPM source, `i_apply_projection` returns `C = cCell` (mode coefficients per hemi) + `gvAll`, so `bst_eigenwavelet('Scalogram', ax, gCell, C)` runs on the inverse's coefficients (no reconstruct-then-project).

- [ ] **Step 1: Failing test** (`dev/test_scalogram_modes.m`)

```matlab
function tests = test_scalogram_modes
tests = functiontests(localfunctions);
end
function test_scalogram_from_mode_c(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st));
    D  = getappdata(st.hFig,'DynamicsOverlay');
    ax = panel_bst_dynamics('i_atom_axes', st, 'Dirac');
    nV = 0; for h=1:2, nV=max(nV,max(ax.GlobalVertices{h}(:))); end
    [C, gvAll] = panel_bst_dynamics('i_apply_projection', st, ax, D, 1:8, nV);
    verifyEqual(t, numel(C), 2);
    verifyEqual(t, size(C{1},1), size(ax.Lambda{1},1));   % coefficients = inverse mode count, per hemi
    verifyEqual(t, size(C{1},2), 8);
end
```

- [ ] **Step 2: Run** — Expected FAIL (current `i_apply_projection` returns scalar-magnitude projection sized `[nModes(60) x nWin]`).

- [ ] **Step 3: Implement** — at the top of `i_apply_projection`, add the Dirac-dSPM branch (before the existing scalar-magnitude body):

```matlab
    if strcmp(i_atom_op(st),'Dirac') && i_is_dirac_dspm(D)
        [cCell,~] = i_mode_coeffs(st, D, iWin);
        C = cCell;  gvAll = [];
        for h=1:numel(ax.GlobalVertices), gvAll=[gvAll; ax.GlobalVertices{h}(:)]; end %#ok<AGROW>
        return;
    end
```

- [ ] **Step 4: Run** — Expected PASS. Also live-check (controller): Analyze-window scalogram now spans the inverse's full `√λ` range (not the old 60-mode cap); Localize-bands produces band atoms.
- [ ] **Step 5: Commit** — `feat(dynamics): scalogram + localize consume the Dirac mode coefficients`.

---

### Task 7: Cortex-magnitude measure toggle (amplitude | dSPM)

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m` (Atom/Preview section: a measure toggle; `i_select_atom_load`/`OnApply` refresh). Test `dev/test_measure_toggle.m` (pure) + live.

**Interfaces:** Consumes Task 3 (`i_dspm_scale`), Task 5 (Apply reads `i_field(st,'atomMeasure','amplitude')`). Produces `OnSetMeasure(name)` storing `st.atomMeasure ∈ {'amplitude','dspm'}` then `i_atom_preview()`; the control is shown only for a Dirac-dSPM session (`i_is_dirac_dspm` true and `i_dspm_scale` non-empty).

- [ ] **Step 1: Failing test** (`dev/test_measure_toggle.m`)

```matlab
function tests = test_measure_toggle
tests = functiontests(localfunctions);
end
function test_measure_default_and_set(t)
    verifyEqual(t, panel_bst_dynamics('i_measure_default'), 'amplitude');
    verifyTrue(t, ismember('dspm', panel_bst_dynamics('i_measure_options')));
    verifyTrue(t, ismember('amplitude', panel_bst_dynamics('i_measure_options')));
end
```

- [ ] **Step 2: Run** — Expected FAIL (`i_measure_default`/`i_measure_options` undefined).
- [ ] **Step 3: Implement** the pure helpers + a two-item toggle (radio/combobox) in the Atom section built only when `i_is_dirac_dspm(D) && ~isempty(i_dspm_scale(st,D))`; `OnSetMeasure(name)` sets `st.atomMeasure` and calls `i_atom_preview()`:

```matlab
function d = i_measure_default(), d = 'amplitude'; end %#ok<DEFNU>
function o = i_measure_options(), o = {'amplitude','dspm'}; end %#ok<DEFNU>
function OnSetMeasure(name) %#ok<DEFNU>
    st = getappdata(0,'DynamicsTarget');  if isempty(st), return; end
    if ~ismember(name, i_measure_options()), return; end
    st.atomMeasure = name;  setappdata(0,'DynamicsTarget', st);  i_atom_preview();
end
```

- [ ] **Step 4: Run** — Expected PASS. Live-check (controller): toggling amplitude↔dSPM rescales the cortex magnitude (per-vertex `SIR`) while quivers/sensor/scalogram are unchanged.
- [ ] **Step 5: Commit** — `feat(dynamics): amplitude|dSPM cortex-magnitude toggle for filtered Dirac atoms`.

---

## Self-Review

**Spec coverage:** mode coefficients → Task 1; session basis → Task 2; dSPM scale + anchor → Task 3; reconstruct/forward → Task 4; Apply wiring (cortex+sensor) → Task 5; scalogram+localize → Task 6; measure toggle → Task 7. Fallback (scalar/non-dSPM) preserved in Tasks 2/5/6 branches. Design impulse "unify on inverse basis" is automatic via Task 2 (the existing `bst_eigenfilter('Atom')` uses whatever `ax` it's given).

**Placeholder scan:** none — every step has concrete code and the verified loader/ordering/anchor. The one test helper stub (`i_launch` in Task 1) is replaced by the controller launching the session for the live gate; the assertion tests are complete.

**Type consistency:** `i_mode_coeffs → [cCell{1×2}[Kh×nWin], meta]` used verbatim in Tasks 4/5/6; `i_atom_axes` Dirac-dSPM ax fields (`Phi/GlobalVertices/Mass/Operator/Lambda`) consistent Tasks 3/4/5/6; `i_dspm_scale → [nV×1]` used in Tasks 5/7; `i_dirac_recon → [V3 [nV×3×nT], mag [nV×nT]]` and `i_dirac_forward_modes → Dfilt` consistent Task 5; measure key `st.atomMeasure ∈ {'amplitude','dspm'}` consistent Tasks 5/7; imag slots `2:4:end/3:4:end/4:4:end` consistent with `RowMap`/Task-4 across the plan.
