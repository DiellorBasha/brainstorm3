# Dirac-Connectome (vector connectome operator) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A `Dirac-Connectome` vector operator whose atoms spread a 3-D current field over white-matter fibers (not the cortical surface), built by lifting the scalar LB-Connectome eigenbasis into the quaternion imaginary slots.

**Architecture:** Lift the scalar LB-Connectome eigenbasis `{Φ_s,Λ_s,M_s}` to a quaternion basis `Φ_q [4nV×3K]`, `Λ_q=Λ_s×3`, `M_q=kron(M_s,I₄)` — i.e. `L_conn ⊗ I₃`, no 3× eigensolve. The `[4n]` quaternion layout reuses the whole Dirac atom/scalogram/decode stack. Apply is reconstruct-then-project (no mode-kernel, no sensor); source-space only. Additive throughout.

**Tech Stack:** MATLAB, Brainstorm; `bst_eigen('Axes')`/`bst_eigenfilter`(`Fiber`/`Atom`/`RowMap`)/`bst_eigenwavelet`(`Scalogram`), `manifold_ft`/`manifold_ift`/`manifold_quat_imag`, `process_helmholtz`.

## Global Constraints

- **Lift:** given scalar LB-Connectome `Φ_s [nV×K]`, `Λ_s [K]`, `M_s [nV×nV]` (single-block whole-brain, `GV{1}=all`, `{2}` empty), produce `Φ_q [4nV×3K]` where mode `k`, axis `d∈{1,2,3}` puts `φ_k` at rows `(v-1)*4 + (d+1)` (imag `x,y,z`; `w=0`); `Λ_q` groups the three axes per mode; `M_q = kron(M_s, I₄)`. Orthonormal: `Φ_qᵀ M_q Φ_q = I_{3K}`.
- **Interleaved quaternion** `[w,x,y,z]`; physical current = imag rows; `Fiber` → C=4 (`quaternion`) via the layout fallback (connectome operators are non-nxr; no `Registry`).
- **`i_gate_mask` order** = `{Laplace-Beltrami, LB-Connectome, Connection Laplacian, Dirac, Dirac-Connectome}`; constrained (`nComponents==1`) → `[1 1 1 0 0]`; unconstrained (`==3`) → `[1 1 0 1 1]`.
- **Apply for Dirac-Connectome** = reconstruct `J` from the source → embed → project onto the lifted `ax` → filter `g(λ)` → reconstruct → cortex magnitude + quivers. NO mode-kernel (that stays `variant=='Dirac' && i_is_dirac_dspm`); NO sensor (`i_dirac_leadfield` returns `[]` for non-`Dirac` variants → `Dfilt=[]`).
- **Helmholtz** on the fiber-spread field: `process_helmholtz('Compute', V3col, Cov)` where `V3col` is the atom-filtered current `[3nV×1]` and `Cov` = the surface Covariant operator (the differential is on the surface, of a field spread over fibers).
- Whole-brain single-block; a fiber impulse legitimately crosses hemispheres. No `clear`; prefer `matlab -batch` for headless tests; commit trailers (development): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>` + `Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL`.

---

### Task 1: Lift helper `bst_lift_connectome_dirac`

**Files:** Create `toolbox/eigen/bst_lift_connectome_dirac.m`; Test `dev/test_lift_connectome_dirac.m`.

**Interfaces:**
- Produces `[Phiq, Lamq, Mq] = bst_lift_connectome_dirac(Phis, Lams, Ms)` — `Phis [n×K]`, `Lams [K×1]`, `Ms [n×n]` → `Phiq [4n×3K]`, `Lamq [3K×1]`, `Mq [4n×4n]` (`kron(Ms,I₄)`). Column order: for `k=1..K`, three columns (x,y,z) → `Lamq = repelem(Lams,3)` grouping, or interleave; the test pins the exact convention.

- [ ] **Step 1: Write the failing test** (`dev/test_lift_connectome_dirac.m`)

```matlab
function tests = test_lift_connectome_dirac
tests = functiontests(localfunctions);
end
function test_lift_shape_orthonormal(t)
    n = 6; K = 4;  rng(0);
    Ms = speye(n);                                  % identity mass -> Phis just needs orthonormal cols
    [Q,~] = qr(randn(n,K),0);  Phis = Q;  Lams = (1:K)';
    [Phiq, Lamq, Mq] = bst_lift_connectome_dirac(Phis, Lams, Ms);
    verifyEqual(t, size(Phiq), [4*n 3*K]);
    verifyEqual(t, size(Mq), [4*n 4*n]);
    verifyEqual(t, numel(Lamq), 3*K);
    % orthonormal in Mq
    G = Phiq' * Mq * Phiq;
    verifyLessThan(t, max(abs(G(:) - reshape(eye(3*K),[],1))), 1e-10);
    % each lifted column has w=0 and lives in exactly one imag axis
    for c = 1:size(Phiq,2)
        col = reshape(Phiq(:,c), 4, n);            % [w;x;y;z] x n
        verifyEqual(t, col(1,:), zeros(1,n), 'AbsTol', 0);      % w = 0
        nzAx = find(any(col(2:4,:) ~= 0, 2));
        verifyEqual(t, numel(nzAx), 1);            % exactly one of x/y/z active
    end
    % Lamq contains each Lams value 3 times
    verifyEqual(t, sort(Lamq), sort(repelem(Lams,3)));
end
function test_kron_mass(t)
    n=3; Ms=sprand(n,n,0.5); Ms=(Ms+Ms')/2+n*speye(n);
    [~,~,Mq]=bst_lift_connectome_dirac(eye(n),(1:n)',Ms);
    verifyLessThan(t, max(abs(full(Mq)-kron(full(Ms),eye(4))),[],'all'), 1e-12);
end
```

- [ ] **Step 2: Run (headless)** `matlab -batch "addpath(genpath('.../toolbox')); runtests('dev/test_lift_connectome_dirac.m')"`. Expected FAIL (function undefined).

- [ ] **Step 3: Implement** `toolbox/eigen/bst_lift_connectome_dirac.m` (Brainstorm header + body):

```matlab
function [Phiq, Lamq, Mq] = bst_lift_connectome_dirac(Phis, Lams, Ms)
% BST_LIFT_CONNECTOME_DIRAC: Lift a scalar eigenbasis to an ambient quaternion (Dirac) basis.
%   [Phiq,Lamq,Mq] = bst_lift_connectome_dirac(Phis,Lams,Ms)
% Each scalar mode phi_k becomes three quaternion modes with phi_k in the imaginary x/y/z slots
% (w=0), eigenvalue lam_k each; mass Mq = kron(Ms,I4). Realizes L_conn (x) I3 (ambient, component-wise).
    Phis = double(Phis);  [n, K] = size(Phis);
    Phiq = zeros(4*n, 3*K);
    Lamq = zeros(3*K, 1);
    for k = 1:K
        for d = 1:3                                  % d=1,2,3 -> imag x,y,z at rows (v-1)*4 + (d+1)
            c = (k-1)*3 + d;
            Phiq((0:n-1)*4 + (d+1), c) = Phis(:, k);
            Lamq(c) = Lams(k);
        end
    end
    Mq = kron(Ms, speye(4));
end
```

- [ ] **Step 4: Run** — Expected 2/2 PASS.
- [ ] **Step 5: Commit** — `feat(eigen): bst_lift_connectome_dirac -- lift a scalar eigenbasis to an ambient quaternion basis`.

---

### Task 2: `Dirac-Connectome` axes in the panel + `RowMap` support

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m` (`i_atom_axes`), `toolbox/eigen/bst_eigenfilter.m` (`RowMap`); Test `dev/test_dirac_connectome_axes.m`.

**Interfaces:**
- Consumes `bst_lift_connectome_dirac` (Task 1), `bst_eigen('Axes', …, 'LB-Connectome')`.
- Produces `ax = i_atom_axes(st, 'Dirac-Connectome')` — lifted quaternion basis: `ax.Phi{1}[4nV×3K]`, `ax.Lambda{1}[3K]`, `ax.Mass{1}=kron(M_s,I₄)`, `ax.GlobalVertices` (same as LB-Connectome, single-block), `ax.Variant='Dirac-Connectome'`, window fields. `RowMap` treats `Dirac-Connectome` like `Dirac`.

- [ ] **Step 1: Failing test** (`dev/test_dirac_connectome_axes.m`)

```matlab
function tests = test_dirac_connectome_axes
tests = functiontests(localfunctions);
end
function test_axes_and_rowmap(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st));   % controller launches a Dirac-dSPM session
    ax = panel_bst_dynamics('i_atom_axes', st, 'Dirac-Connectome');
    axs = bst_eigen('Axes', struct('SurfaceFile',ax.SurfaceFile,'Variant','LB-Connectome','nModes',60,'TimeWindow',[0 .04],'SampleRate',100));
    nV = numel(axs.GlobalVertices{1});  K = size(axs.Phi{1},2);
    verifyEqual(t, size(ax.Phi{1}), [4*nV 3*K]);
    verifyEqual(t, numel(ax.Lambda{1}), 3*K);
    [C,kind] = bst_eigenfilter('Fiber', ax);
    verifyEqual(t, C, 4);  verifyEqual(t, kind, 'quaternion');
    % RowMap maps a 3-vector source into the quaternion imag slots for this variant
    F = zeros(3*nV,1);  [srcRows, dstRows, nrows, msg] = bst_eigenfilter('RowMap', F, ax, 1);
    verifyEmpty(t, msg);  verifyEqual(t, nrows, 4*nV);
end
```

- [ ] **Step 2: Run** — FAIL (`i_atom_axes` has no Dirac-Connectome branch; `RowMap` errors on the variant).

- [ ] **Step 3a:** In `bst_eigenfilter`'s `RowMap`, add `'Dirac-Connectome'` to the Dirac case (the `case {'Dirac','Dirac-Face','Hodge-Face'}` line → add `'Dirac-Connectome'`; it uses `EigenMat.GlobalVertices{h}` and the `(idx-1)*3+[1,2,3]`→`(0:n-1)*4+[2,3,4]` mapping, which is correct for the lifted basis).

- [ ] **Step 3b:** In `i_atom_axes`, add a `Dirac-Connectome` branch (before the canonical `bst_eigen('Axes')` body, after the Dirac-dSPM branch):

```matlab
    if strcmp(variant, 'Dirac-Connectome')
        surf = i_atom_surface(st);  key = ['dconn|' surf];
        Mc = getappdata(0,'DynamicsAtomAx');  if isempty(Mc)||~isa(Mc,'containers.Map'), Mc=containers.Map('KeyType','char','ValueType','any'); end
        if isKey(Mc,key), ax = Mc(key); return; end
        Fs = 100; D = getappdata(st.hFig,'DynamicsOverlay');
        if ~isempty(D) && isfield(D,'srcDS') && isfield(D,'srcResult')
            try, tv = bst_memory('GetTimeVector', D.srcDS, D.srcResult); if numel(tv)>1, Fs=1/median(diff(tv)); end, catch, end %#ok<CTCH>
        end
        nF = max(2, round(4*Fs));
        axs = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','LB-Connectome','nModes',60,'TimeWindow',[0 (nF-1)/Fs],'SampleRate',Fs));
        if isempty(axs) || isempty(axs.Phi{1}), return; end
        [Phiq,Lamq,Mq] = bst_lift_connectome_dirac(axs.Phi{1}, axs.Lambda{1}(:), axs.Mass{1});
        ax = struct('Variant','Dirac-Connectome','SurfaceFile',surf, ...
                    'Phi',{{Phiq,[]}}, 'GlobalVertices',{axs.GlobalVertices}, 'Mass',{{Mq,[]}}, 'Lambda',{{Lamq,[]}});
        ax.nT = nF;  ax.tlag = (0:nF-1)/Fs;  ax.omega = (0:nF-1)*(Fs/nF);  ax.NFFT = nF;
        Mc(key) = ax;  setappdata(0,'DynamicsAtomAx', Mc);
        return;
    end
```

- [ ] **Step 4: Run** — Expected PASS (shape, Fiber quaternion, RowMap ok). Also live: `bst_eigenfilter('Atom', ax, 'heat', kp, seed)` realises (impulse spreads).
- [ ] **Step 5: Commit** — `feat(dynamics): Dirac-Connectome axes (lifted LB-Connectome) + RowMap support`.

---

### Task 3: Operator dropdown + `i_gate_mask` (5th operator)

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m` (`opDefs`, `i_gate_mask`); Test `dev/test_gate_mask_dconn.m`.

**Interfaces:** `opVariants` gains `'Dirac-Connectome'` (label "Dirac (connectome)"); `i_gate_mask(nComp)` returns a 5-logical.

- [ ] **Step 1: Failing test** (`dev/test_gate_mask_dconn.m`)

```matlab
function tests = test_gate_mask_dconn
tests = functiontests(localfunctions);
end
function test_masks(t)
    verifyEqual(t, panel_bst_dynamics('i_gate_mask', 1), logical([1 1 1 0 0]));  % constrained
    verifyEqual(t, panel_bst_dynamics('i_gate_mask', 3), logical([1 1 0 1 1]));  % unconstrained (Dirac + Dirac-Connectome)
    verifyEqual(t, panel_bst_dynamics('i_gate_mask', []), logical([1 1 1 1 1]));
end
```

- [ ] **Step 2: Run** — FAIL (masks are 4-element).
- [ ] **Step 3:** Update `opDefs` (add `'Dirac (connectome)','Dirac-Connectome'` as the 5th row) and `i_gate_mask`:

```matlab
function m = i_gate_mask(nComponents) %#ok<DEFNU>
    % Order: {Laplace-Beltrami, LB-Connectome, Connection Laplacian, Dirac, Dirac-Connectome}
    m = true(1,5);
    if     isequal(nComponents, 1), m = logical([1 1 1 0 0]);   % constrained: scalar bases + tangent
    elseif isequal(nComponents, 3), m = logical([1 1 0 1 1]);   % unconstrained: scalar-norm + Dirac + Dirac-Connectome
    end
end
```

- [ ] **Step 4: Run** — 1/1 PASS. (Update `dev/test_operator_gate.m` expectations to 5-element in the same commit.)
- [ ] **Step 5: Commit** — `feat(dynamics): add Dirac-Connectome operator + 5-way compat gate`.

---

### Task 4: `i_vector_coeffs` + Apply (cortex + quivers, no sensor)

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m` (`i_vector_coeffs` new; `i_atom_apply`; `i_dirac_leadfield`). Test: live (controller).

**Interfaces:**
- Produces `cCell = i_vector_coeffs(st, ax, D, iWin)` — reconstruct `J=GetResultsValues[3nV×nWin]`, embed into quaternion source layout, project per block: `cCell{h}=manifold_ft(ax.Phi{h}, ax.Mass{h}, U_h)`. (General reconstruct-then-project for a vector `ax`; analog of `i_mode_coeffs` but from the field.)
- `i_atom_apply`: admit `Dirac-Connectome` into the vector branch; use `i_vector_coeffs` → filter → `i_dirac_recon` → cortex magnitude + `QuiverVectorOverride`; `Dfilt` empty → no sensor.
- `i_dirac_leadfield`: return `[]` unless `variant=='Dirac'` (Dirac-Connectome has no leadfield).

- [ ] **Step 1:** Add `i_vector_coeffs` (embed `J` rows `(v-1)*3+{1,2,3}` → quaternion imag `(v-1)*4+{2,3,4}`, `w=0`, per block via `RowMap`, then `manifold_ft`), mirroring the embed in `bst_eigenfilter('Atom')`/`RowMap`.

- [ ] **Step 2:** In `i_atom_apply`, change the vector gate from `if strcmp(variant,'Dirac')` to `if any(strcmp(variant,{'Dirac','Dirac-Connectome'}))`; inside, keep the mode-kernel sub-branch gated on `strcmp(variant,'Dirac') && i_is_dirac_dspm(D)`; add an `else` (covers Dirac-Connectome AND non-dSPM Dirac): `cCell=i_vector_coeffs(st,ax,D,iWin); cf=g(λ)·cCell; [V3,mag]=i_dirac_recon(ax,cf);` paint magnitude + set `QuiverVectorOverride=V3(:,:,mid)`; `Leig=i_dirac_leadfield(st,ax)` (→ [] for Dirac-Connectome) so `i_dirac_forward_modes` yields `Dfilt=[]` → no sensor overlay.

- [ ] **Step 3:** In `i_dirac_leadfield`, add at the top: `if ~strcmp(ax.Variant,'Dirac'), Leig=[]; return; end` (only the surface Dirac has a physical eigen-leadfield).

- [ ] **Step 4: Live verify (controller, MCP):** unconstrained Dirac-dSPM session → select "Dirac (connectome)" → Create atom (static heat) → Apply → cortex magnitude + quivers that spread over **fibers** (visibly different extent from surface-Dirac; crosses hemispheres); info reports no sensor overlay. Screenshot. Confirm `getappdata(hFig,'QuiverVectorOverride')` non-empty, `Dfilt` empty path taken.
- [ ] **Step 5: Commit** — `feat(dynamics): Dirac-Connectome Apply -- fiber-spread cortex field + quivers (no sensor)`.

---

### Task 5: Scalogram + Localize for Dirac-Connectome

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m` (`i_apply_projection`; the `OnAnalyzeWindow`/`OnLocalizeBands` gates). Test `dev/test_scalogram_dconn.m`.

**Interfaces:** `i_apply_projection` returns the vector coefficients for `Dirac-Connectome`; the Analyze/Localize gates admit it.

- [ ] **Step 1: Failing test** (`dev/test_scalogram_dconn.m`)

```matlab
function tests = test_scalogram_dconn
tests = functiontests(localfunctions);
end
function test_scalogram(t)
    st = getappdata(0,'DynamicsTarget');  assert(~isempty(st));
    D  = getappdata(st.hFig,'DynamicsOverlay');
    ax = panel_bst_dynamics('i_atom_axes', st, 'Dirac-Connectome');
    nV = numel(ax.GlobalVertices{1});
    [C, ~] = panel_bst_dynamics('i_apply_projection', st, ax, D, 1:8, nV);
    verifyEqual(t, size(C{1},1), size(ax.Lambda{1},1));   % 3K vector coeffs
    lmax = max(ax.Lambda{1}(:)); N=4;
    gC = cell(1,N); for m=1:N, gC{m}=bst_eigfilter_design_itersine(struct('member',m,'Nf',N,'lmax',lmax)); end
    scal = bst_eigenwavelet('Scalogram', ax, gC, C);
    verifyEqual(t, size(scal.W,1), nV);                  % per-vertex magnitude (not 4nV)
    verifyEqual(t, size(scal.energy,3), N);
end
```

- [ ] **Step 2: Run** — FAIL (`i_apply_projection` returns the scalar-magnitude projection; gates exclude Dirac-Connectome).
- [ ] **Step 3a:** In `i_apply_projection`, add before the scalar body: `if strcmp(i_atom_op(st),'Dirac-Connectome'), C = i_vector_coeffs(st, ax, D, iWin); gvAll = [ax.GlobalVertices{1}(:)]; return; end`.
- [ ] **Step 3b:** In `OnAnalyzeWindow` and `OnLocalizeBands`, extend the scalar-only gate to admit `Dirac-Connectome`: `&& ~strcmp(variant,'Dirac-Connectome')` on the bail condition (alongside the existing Dirac-dSPM admission).
- [ ] **Step 4: Run** — PASS. Live: Analyze on a Dirac-Connectome atom → scalogram spanning `√λ`; Localize → band atoms.
- [ ] **Step 5: Commit** — `feat(dynamics): scalogram + localize for Dirac-Connectome`.

---

### Task 6: Helmholtz on the fiber-spread field

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m` (a `OnHelmholtzFiltered` action + wiring). Test: live.

**Interfaces:** an action that runs `process_helmholtz('Compute', V3col, Cov)` on the current Dirac-Connectome-filtered field and displays/reports div/curl.

- [ ] **Step 1:** Add `OnHelmholtzFiltered()`: rebuild the filtered field for the cursor frame (`i_vector_coeffs` → `g(λ)` → `i_dirac_recon` → `V3`), take the mid-frame `V3col=[3nV×1]`, load `Cov=getappdata(st.hFig,'DynamicsOverlay').Cov`, `Ht=process_helmholtz('Compute', V3col, Cov)`, and paint/report the divergence (or curl) magnitude via the existing `view_helmholtz`/differential-overlay path. Gate it to `Dirac-Connectome` (and Dirac) atoms.
- [ ] **Step 2:** Add a toolbar/menu entry "Helmholtz (filtered)" that calls it (near the Analyze/Localize buttons).
- [ ] **Step 3: Live verify (controller):** Dirac-Connectome atom → Helmholtz (filtered) → a div/curl map computed from the fiber-spread field (differs from the raw-source Helmholtz); no error; `Ht.Div`/`Ht.Curl` finite. Screenshot.
- [ ] **Step 4: Commit** — `feat(dynamics): Helmholtz/differential on the Dirac-Connectome-filtered field`.

---

## Self-Review

**Spec coverage:** lift → Task 1; variant axes + RowMap → Task 2; operator UI + gate → Task 3; Apply (fiber-spread cortex+quivers, no sensor) → Task 4; scalogram + localize → Task 5; Helmholtz on filtered field → Task 6. Impulse (Design) works after Tasks 1–3 (generic realiser on the lifted ax). Out-of-scope (sensor, pure-fiber variant, fiber-orientation coupling) correctly absent.

**Placeholder scan:** none — lift code, ax branch, gate masks, RowMap edit, apply/projection branches, and the Helmholtz call are all concrete; headless run command spelled out.

**Type consistency:** `bst_lift_connectome_dirac(Phis,Lams,Ms)→[Phiq,Lamq,Mq]` used verbatim in Task 2; the lifted `ax` fields (`Phi{1}[4nV×3K]`, `Lambda{1}[3K]`, `Mass{1}`, `GlobalVertices`) consistent Tasks 2/4/5; `i_vector_coeffs(st,ax,D,iWin)→cCell{1×2}` consistent Tasks 4/5/6; `i_gate_mask`→5-logical consistent Task 3; quaternion imag layout (`(v-1)*4+{2,3,4}`, `manifold_quat_imag`) consistent with `RowMap`/`i_dirac_recon` across the plan; `process_helmholtz('Compute', V3col, Cov)` per its signature.
