# Dynamics Frame Core (Sub-project B) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the Dynamics atom bank into a real frame — live frame-bounds/coverage readout, a one-click Design-tight-frame (itersine) generate, and cached windowed-source projection that makes the selected atom's Apply instant.

**Architecture:** Add `itersine` as a first-class eigfilter registry kernel so generated tight-frame members are ordinary atoms. Compute bounds by reusing `bst_eigenwavelet('Evaluate'/'Bounds')` over an ad-hoc frame built from the bank's kernels. Render coverage by reusing+extending `view_eigfilter_response` (its existing bank mode). Add a docked **Frame** GUI section (numbers + N spinner + Design-tight-frame button + Show-coverage toggle) to `panel_bst_dynamics.m`, and refactor `i_atom_apply` into project-once + apply-gain with a cached projection.

**Tech Stack:** MATLAB (Brainstorm fork). Verification via the brainstorm-dev MATLAB MCP (`check_matlab_code`, and controller-run live/headless) + `git`/`grep`.

## Global Constraints

- Files touched: `toolbox/eigen/eigfilter/bst_eigfilter_design_itersine.m` (new), `toolbox/gui/view_eigfilter_response.m`, `toolbox/gui/panel_bst_dynamics.m`, `toolbox/gui/panel_eigenfilter_design.m` (one-line exclusion), and new headless test files under `dev/`. No other files.
- **Behavior-preserving for existing features:** the `view_eigfilter_response` extension must be backward-compatible (new `Coverage`/`Bounds` fields optional; existing single-kernel and bank callers unaffected). The `i_atom_apply` refactor must produce output identical to the current path.
- **Static-only for implementers:** implementers run `git`, `grep`, `check_matlab_code`, and WRITE headless test `.m` files but DO NOT run/evaluate MATLAB or launch Brainstorm (it disturbs the live session). The controller runs the consolidated live pass (headless tests + GUI smoke) at the end.
- Kernel-registry convention: `bst_eigfilter_design_<name>(params)` returns a handle `@(l)` (or a cell bank for a vector param); `bst_eigfilter_design_<name>('meta')` returns `struct('name','display','params',...,'bandpass',...,'priorAdmissible',...)`. `bst_eigfilter_kernel('list')` auto-discovers by filename; `('info',name)` returns the meta.
- Itersine member formula (verbatim from `bst_eigenwavelet('Design','itersine')`): `overlap=2; scale=lmax/(Nf-overlap+1)*overlap; kf=@(x) sin(0.5*pi*(cos(pi*x)).^2).*(x>=-0.5 & x<=0.5); g=@(l) kf(l/scale-(member-overlap/2)/overlap)./sqrt(overlap).*sqrt(2)`.
- Frame bounds: `S(λ)=Σ_m|g_m(λ)|²`, `A=min S`, `B=max S`, tightness `=B/A` (✓ when `|B/A-1|<0.05`).
- Cache invalidation key for the Apply projection: `(srcResult, iWin, operator)` — NOT seed or kernel params.
- Commit after each task: `feat(dynamics): <summary>` (or `refactor(dynamics):` for Task 5) with the standard session trailers.
- Branch: `development`. Do not push.

---

### Task 1: `itersine` tight-frame member registry kernel

**Files:**
- Create: `toolbox/eigen/eigfilter/bst_eigfilter_design_itersine.m`
- Create (test): `dev/test_itersine_kernel.m`

**Interfaces:**
- Consumes: nothing.
- Produces: `bst_eigfilter_design_itersine(struct('member',ii,'Nf',N,'lmax',L))` → handle `@(l)`; `bst_eigfilter_design_itersine('meta')` → metadata. Discoverable via `bst_eigfilter_kernel('list')` and buildable via `bst_eigfilter_kernel('itersine', kp)`.

- [ ] **Step 1: Write the failing headless test**

Create `dev/test_itersine_kernel.m`:
```matlab
function tests = test_itersine_kernel
tests = functiontests(localfunctions);
end

function test_meta(tc)
m = bst_eigfilter_design_itersine('meta');
verifyEqual(tc, m.name, 'itersine');
verifyTrue(tc, m.bandpass);
end

function test_single_member_handle(tc)
g = bst_eigfilter_design_itersine(struct('member',2,'Nf',6,'lmax',10));
verifyTrue(tc, isa(g,'function_handle'));
y = g(linspace(0,10,50)');
verifyEqual(tc, numel(y), 50);
verifyGreaterThanOrEqual(tc, min(y), -1e-12);      % itersine window is non-negative
end

function test_tight_frame_property(tc)
% Sum of squares of all Nf members is ~constant over the interior -> tight (B/A ~ 1)
Nf = 6; lmax = 10;
lam = linspace(0.1*lmax, 0.9*lmax, 400)';          % interior (edges roll off)
S = zeros(size(lam));
for ii = 1:Nf
    g = bst_eigfilter_design_itersine(struct('member',ii,'Nf',Nf,'lmax',lmax));
    S = S + g(lam).^2;
end
A = min(S); B = max(S);
verifyLessThan(tc, B/A, 1.05);                     % tight frame on the interior
end
```

- [ ] **Step 2: (implementer) confirm the test references only the not-yet-created function**

The test calls `bst_eigfilter_design_itersine`, which does not exist yet — it will error until Step 3. Do NOT run MATLAB (controller runs it). Record that the RED state is "undefined function bst_eigfilter_design_itersine".

- [ ] **Step 3: Create the kernel**

Create `toolbox/eigen/eigfilter/bst_eigfilter_design_itersine.m`:
```matlab
function out = bst_eigfilter_design_itersine(params)
% BST_EIGFILTER_DESIGN_ITERSINE: One member of an itersine (half-cosine) TIGHT frame.
% The Nf members' squared responses sum to ~constant over the interior (a tight frame),
% so a bank of these is a robust frame (B/A -> 1). Used by the Dynamics "Design tight
% frame" generate; parameters (member/Nf/lmax) are set programmatically, not slider-tuned.
% USAGE:  g = bst_eigfilter_design_itersine(struct('member',ii,'Nf',N,'lmax',L))  -> handle
%         m = bst_eigfilter_design_itersine('meta')                                -> metadata
%
% Authors: Diellor Basha, 2026

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@

if nargin >= 1 && ischar(params) && strcmpi(params,'meta')
    out = struct('name','itersine','display','Itersine (tight-frame member)', ...
        'params', struct('member',struct('default',1,'range',[1 Inf]), ...
                         'Nf',struct('default',6,'range',[2 Inf]), ...
                         'lmax',struct('default',1,'range',[0 Inf])), ...
        'bandpass', true, 'priorAdmissible', false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
member = 1;  if isfield(params,'member') && ~isempty(params.member); member = params.member; end
Nf     = 6;  if isfield(params,'Nf')     && ~isempty(params.Nf);     Nf     = params.Nf;     end
lmax   = 1;  if isfield(params,'lmax')   && ~isempty(params.lmax);   lmax   = params.lmax;   end
if ~(lmax > 0); error('bst_eigfilter_design_itersine: lmax must be > 0.'); end
overlap = 2;
scale   = lmax / (Nf - overlap + 1) * overlap;
kf      = @(x) sin(0.5*pi*(cos(pi*x)).^2) .* (x >= -0.5 & x <= 0.5);
out     = @(l) kf(double(l(:))/scale - (member - overlap/2)/overlap) ./ sqrt(overlap) .* sqrt(2);
end
```

- [ ] **Step 4: Static check**

Run `check_matlab_code` on `bst_eigfilter_design_itersine.m` and `dev/test_itersine_kernel.m`. Expected: no errors.

- [ ] **Step 5: Grep — registry discovery**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
ls toolbox/eigen/eigfilter/bst_eigfilter_design_itersine.m && echo "kernel file present"
```
Expected: the file is present (so `bst_eigfilter_kernel('list')` will include `itersine`).

- [ ] **Step 6: Commit**

```bash
git add toolbox/eigen/eigfilter/bst_eigfilter_design_itersine.m dev/test_itersine_kernel.m
git commit -m "feat(dynamics): itersine tight-frame member kernel + headless test

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

*(Controller runs `test_itersine_kernel` live in the consolidated pass.)*

---

### Task 2: extend `view_eigfilter_response` with Coverage + Bounds (bank mode)

**Files:**
- Modify: `toolbox/gui/view_eigfilter_response.m` (bank-mode branch, currently ~lines 51-72)

**Interfaces:**
- Consumes: nothing new.
- Produces: `view_eigfilter_response(g, lambdas, titleStr)` bank mode now honors two OPTIONAL fields on `g`: `g.Coverage` (bool → overlay `S(λ)=Σ_m g_m(λ)²` as a bold line + a faint mean reference) and `g.Bounds` (struct `A`,`B`,`tightness` → annotate the title). Existing callers (no `Coverage`/`Bounds`) are unchanged.

- [ ] **Step 1: Add the Coverage/Bounds overlay to the bank-mode branch**

In `view_eigfilter_response.m`, the bank-mode branch (`if isstruct(g) && isfield(g,'Kernels')`) currently ends by setting `xlim`/labels/title. Insert the coverage overlay just before `xlim(ax, [0 lmax]);` (after the per-member plotting loop), and extend the title. Replace this block:
```matlab
    xlim(ax, [0 lmax]);
    xlabel(ax, '\lambda (eigenvalue)'); ylabel(ax, 'g(\lambda)'); grid(ax, 'on');
    if (nargin >= 3) && ~isempty(titleStr); title(ax, titleStr, 'Interpreter','none'); end
    return;
```
with:
```matlab
    % --- optional coverage overlay: S(lambda) = sum_m g_m(lambda)^2 (frame tightness) ---
    if isfield(g,'Coverage') && ~isempty(g.Coverage) && g.Coverage
        S = zeros(size(xg));
        for j = 1:nT; yj = g.Kernels{j}(xg); S = S + real(yj(:)).^2; end
        plot(ax, xg, S, '-', 'Color',[.1 .1 .8], 'LineWidth', 2.5);
        plot(ax, [0 lmax], mean(S)*[1 1], ':', 'Color',[.1 .1 .8 ]);   % faint "flat" reference
    end
    xlim(ax, [0 lmax]);
    xlabel(ax, '\lambda (eigenvalue)'); ylabel(ax, 'g(\lambda) / S(\lambda)'); grid(ax, 'on');
    ttl = '';  if (nargin >= 3) && ~isempty(titleStr); ttl = titleStr; end
    if isfield(g,'Bounds') && ~isempty(g.Bounds) && isstruct(g.Bounds)
        chk = ''; if abs(g.Bounds.tightness - 1) < 0.05; chk = ' OK'; end
        ttl = sprintf('%s   A=%.3g  B=%.3g  B/A=%.3g%s', ttl, g.Bounds.A, g.Bounds.B, g.Bounds.tightness, chk);
    end
    if ~isempty(ttl); title(ax, ttl, 'Interpreter','none'); end
    return;
```

- [ ] **Step 2: Static check**

Run `check_matlab_code` on `view_eigfilter_response.m`. Expected: no errors.

- [ ] **Step 3: Backward-compat grep**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
grep -n "isfield(g,'Coverage')\|isfield(g,'Bounds')" toolbox/gui/view_eigfilter_response.m
```
Expected: both guards present (so callers without those fields are unaffected). Confirm the single-kernel branch (non-struct `g`, below the bank branch) is untouched.

- [ ] **Step 4: Commit**

```bash
git add toolbox/gui/view_eigfilter_response.m
git commit -m "feat(dynamics): view_eigfilter_response bank mode gains Coverage + Bounds overlay

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 3: panel — `i_frame_response`, Frame GUI section, live refresh

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (`CreatePanel` + new subfunctions + wire into existing callbacks)

**Interfaces:**
- Consumes: `bst_eigenwavelet('Evaluate'/'Bounds')`; `bst_eigfilter_kernel(name, kp)`; the existing `i_atom_axes(st, variant)`, `i_atom_op(st)`, `i_cs()`, `bst_eigfilter_controls('ToKernel', ...)`.
- Produces: `i_frame_response(st, ax)` → `struct('lam',..,'S',..,'A',..,'B',..,'tightness',..,'nMembers',..,'gCell',{..})`; `i_frame_refresh()`; a `jFrameA/jFrameB/jFrameT` label set, a `jFrameN` spinner, and a `jFrameShow` toggle stored in the panel `ctrl` struct.

- [ ] **Step 1: Add the bounds-compute subfunction**

Add to `panel_bst_dynamics.m`:
```matlab
% Frame response over the CURRENT operator's atoms: coverage S(lambda) + bounds A/B/tightness.
function fr = i_frame_response(st, ax)
    fr = struct('lam',[], 'S',[], 'A',NaN, 'B',NaN, 'tightness',NaN, 'nMembers',0, 'gCell',{{}});
    if isempty(ax) || ~isfield(ax,'Lambda') || isempty(ax.Lambda), return; end
    variant = i_atom_op(st);
    lamAll = ax.Lambda{1}(:);
    if numel(ax.Lambda) > 1 && ~isempty(ax.Lambda{2}), lamAll = [lamAll; ax.Lambda{2}(:)]; end
    lmax = max(lamAll);  lminPos = min(lamAll(lamAll > 1e-9));  if isempty(lminPos), lminPos = 0; end
    % gather g(lambda) for every atom on the current operator
    gCell = {};
    for i = 1:numel(st.T.Groups)
        G = st.T.Groups(i);
        op = 'Laplace-Beltrami'; if isfield(G,'Operator') && ~isempty(G.Operator), op = G.Operator; end
        if ~strcmp(op, variant), continue; end
        kp = G.KernelParams;  if ~isstruct(kp), kp = struct(); end
        if ~isfield(kp,'lmax') || isempty(kp.lmax), kp.lmax = lmax; end
        try, g = bst_eigfilter_kernel(G.KernelName, kp); catch, continue; end
        if iscell(g), gCell = [gCell, g(:)']; else, gCell{end+1} = g; end %#ok<AGROW>
    end
    fr.nMembers = numel(gCell);  fr.gCell = gCell;
    if fr.nMembers == 0, return; end
    lam = linspace(max(lminPos,eps), lmax, 400)';
    S = zeros(size(lam));
    for m = 1:numel(gCell), y = gCell{m}(lam); S = S + real(y(:)).^2; end
    fr.lam = lam;  fr.S = S;  fr.A = min(S);  fr.B = max(S);
    if fr.nMembers >= 2 && fr.A > 0, fr.tightness = fr.B / fr.A; end   % undefined for a single band
end
```

- [ ] **Step 2: Add the live-refresh subfunction**

```matlab
% Recompute the frame response and update the Frame labels (+ the coverage view if shown).
function i_frame_refresh()
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    if ~isfield(ctrl,'jFrameA'), return; end
    ax = [];  try, ax = i_atom_axes(st, i_atom_op(st)); catch, end
    fr = i_frame_response(st, ax);
    if fr.nMembers == 0
        ctrl.jFrameA.setText('A —');  ctrl.jFrameB.setText('B —');  ctrl.jFrameT.setText('B/A —');
    else
        ctrl.jFrameA.setText(sprintf('A %.3g', fr.A));
        ctrl.jFrameB.setText(sprintf('B %.3g', fr.B));
        if isnan(fr.tightness), ctrl.jFrameT.setText('B/A —');
        else
            chk = ''; if abs(fr.tightness-1) < 0.05, chk = ' ✓'; end
            ctrl.jFrameT.setText(sprintf('B/A %.3g%s', fr.tightness, chk));
        end
    end
    if isfield(ctrl,'jFrameShow') && ~isempty(ctrl.jFrameShow) && ctrl.jFrameShow.isSelected() && ~isempty(ax)
        bnd = struct('A',fr.A,'B',fr.B,'tightness',fr.tightness);
        gstruct = struct('Kernels',{fr.gCell}, 'Active', max(1,i_field(st,'curAtom',1)), ...
            'OnSelect', @(j)bst_call(@()SetSelectedAtom(j)), 'Coverage', true, 'Bounds', bnd);
        lamMark = ax.Lambda{1}(:);
        try, view_eigfilter_response(gstruct, lamMark, sprintf('Frame coverage (%s)', i_atom_op(st))); catch, end
    end
end
```

- [ ] **Step 3: Build the Frame GUI section in `CreatePanel`**

In `CreatePanel`, after the `jAtom` section is added to `jPanelMain` (SOUTH), add a Frame section. Insert before `jPanelNew.add(jPanelMain, BorderLayout.CENTER);`:
```matlab
    % Frame section (SOUTH, below Atom): bounds readout + tight-frame generate + coverage toggle
    jFrame = gui_river([0 0], [0 2 0 2], 'Frame');
    jFrameA = gui_component('label', jFrame, [], 'A —');
    jFrameB = gui_component('label', jFrame, 'tab', 'B —');
    jFrameT = gui_component('label', jFrame, 'tab', 'B/A —');
    jFrame.add('br', javax.swing.JLabel('N'));
    jFrameN = gui_component('spinner', jFrame, 'tab', []);
    jFrameN.setModel(javax.swing.SpinnerNumberModel(int32(6), int32(2), int32(24), int32(1)));
    gui_component('button', jFrame, 'tab', 'Design tight frame', [], 'Replace the bank with an itersine tight frame of N members', @(h,e)bst_call(@OnDesignFrame));
    jFrameShow = gui_component('checkbox', jFrame, 'br', 'Show coverage', [], 'Open/close the frame coverage response view', @(h,e)bst_call(@OnFrameShow));
    jPanelMain.add(jFrame, BorderLayout.SOUTH);
```
Note: the Atom section currently uses `BorderLayout.SOUTH` of `jPanelMain`. Wrap the existing Atom section and this Frame section into a vertical container so both fit SOUTH — create `jSouth = gui_component('Panel'); jSouth.setLayout(BoxLayout(jSouth, BoxLayout.Y_AXIS)); jSouth.add(jAtom); jSouth.add(jFrame); jPanelMain.add(jSouth, BorderLayout.SOUTH);` (replacing the current `jPanelMain.add(jAtom, BorderLayout.SOUTH);`).

Add the four handles to the returned `BstPanel` struct: `'jFrameA',jFrameA, 'jFrameB',jFrameB, 'jFrameT',jFrameT, 'jFrameN',jFrameN, 'jFrameShow',jFrameShow`.

- [ ] **Step 4: Add the Show-coverage toggle handler**

```matlab
% Show-coverage toggle: open/refresh or close the frame coverage response view.
function OnFrameShow() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end %#ok<ASGLU>
    if isfield(ctrl,'jFrameShow') && ~isempty(ctrl.jFrameShow) && ctrl.jFrameShow.isSelected()
        i_frame_refresh();
    else
        try, view_eigfilter_response('close'); catch, end %#ok<CTCH>
    end
end
```

- [ ] **Step 5: Wire `i_frame_refresh()` into the bank-mutating callbacks**

Add a call to `i_frame_refresh();` at the end of each of these existing subfunctions (just before they return): `OnCreateAtom`, `i_atom_writeback`, `AtomsListValueChanged_Callback`, `OnSetOperator`, `AtomDeleteGroup`, `i_select_atom_load`. (These are every path that adds/edits/selects/deletes an atom or changes the operator.)

- [ ] **Step 6: Static check + grep**

Run `check_matlab_code` on `panel_bst_dynamics.m`. Expected: no undefined-reference errors.
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
grep -c "i_frame_refresh" toolbox/gui/panel_bst_dynamics.m   # expect >=7 (def + >=6 call sites)
grep -c "jFrameA\|jFrameShow" toolbox/gui/panel_bst_dynamics.m # expect >0
```

- [ ] **Step 7: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "feat(dynamics): live frame bounds/coverage (Frame section + i_frame_response)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 4: `OnDesignFrame` generate + exclude itersine from the hand-pick combobox

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (add `OnDesignFrame`)
- Modify: `toolbox/gui/panel_eigenfilter_design.m` (`AtomKernels` excludes `itersine`)

**Interfaces:**
- Consumes: Task 1's `itersine` kernel; Task 3's `jFrameN` spinner + `i_frame_refresh`; existing `i_atom_axes`, `i_default_atom`, `bst_dynamics('AddGroup')`, `UpdateAtomList`, `SetSelectedAtom`.
- Produces: `OnDesignFrame()` replaces the bank with N itersine members on the current operator.

- [ ] **Step 1: Exclude itersine from the hand-pick kernel list**

In `panel_eigenfilter_design.m`, `AtomKernels()` starts `keys = bst_eigfilter_kernel('list');`. Insert immediately after:
```matlab
    keys = keys(~strcmpi(keys, 'itersine'));   % itersine is generate-only (Design tight frame), not hand-picked
```
(So itersine never appears in the Filter combobox, but remains a valid realisable kernel.)

- [ ] **Step 2: Add `OnDesignFrame`**

Add to `panel_bst_dynamics.m`:
```matlab
% Design tight frame: replace the bank with N itersine members spanning the current operator's spectrum.
function OnDesignFrame() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    op = i_atom_op(st);
    ax = i_atom_axes(st, op);  if isempty(ax), return; end
    lamAll = ax.Lambda{1}(:);
    if numel(ax.Lambda) > 1 && ~isempty(ax.Lambda{2}), lamAll = [lamAll; ax.Lambda{2}(:)]; end
    lmax = max(lamAll);  lminPos = min(lamAll(lamAll > 1e-9));  if isempty(lminPos), lminPos = eps; end
    N = 6;  if isfield(ctrl,'jFrameN') && ~isempty(ctrl.jFrameN), N = double(ctrl.jFrameN.getValue()); end
    N = max(2, round(N));
    if ~isempty(st.T.Groups)
        if ~java_dialog('confirm', sprintf('Replace the current %d atom(s) with a %d-member itersine tight frame?', numel(st.T.Groups), N), 'Design tight frame')
            return;
        end
        st.T.Groups(:) = [];  st.T.nGroups = 0;
    end
    seed = ax.GlobalVertices{1}(1);
    for ii = 1:N
        kp = struct('member',ii, 'Nf',N, 'lmin',lminPos, 'lmax',lmax, 'vals',[]);
        G  = i_default_atom('itersine', kp, seed, ax.SurfaceFile, sprintf('itersine %d/%d', ii, N), op);
        st.T = bst_dynamics('AddGroup', st.T, G);
    end
    setappdata(0, 'DynamicsTarget', st);
    UpdateAtomList();
    SetSelectedAtom(1);              % selects member 1, loads it, and triggers i_frame_refresh
    bst_progress('text', sprintf('Designed %d-member itersine tight frame on %s', N, op));
end
```

- [ ] **Step 3: Static check + grep**

Run `check_matlab_code` on both files. Expected: no errors.
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
grep -c "OnDesignFrame" toolbox/gui/panel_bst_dynamics.m         # def + the CreatePanel button wiring (>=2)
grep -n "strcmpi(keys, 'itersine')" toolbox/gui/panel_eigenfilter_design.m  # the exclusion line
```

- [ ] **Step 4: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m toolbox/gui/panel_eigenfilter_design.m
git commit -m "feat(dynamics): Design-tight-frame generate (itersine bank); hide itersine from hand-pick

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 5: cached-projection instant Apply

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (`i_atom_apply` refactor + a cache helper)

**Interfaces:**
- Consumes: existing `i_atom_apply` context (`D.srcDS`/`D.srcResult`, `i_cursor_window`, `i_paintable_scalar`, `i_overlay_nv`, `i_atom_axes`, `bst_eigfilter_controls('ToKernel')`, `manifold_ft`/`manifold_ift`).
- Produces: `i_apply_projection(st, ax, D, iWin, nV)` → per-hemi cached modal coefficients `C` (cell) for `(srcResult, iWin, operator)`; `i_atom_apply` reuses it so slider-drag re-synthesis skips the fetch+projection.

- [ ] **Step 1: Add the projection cache helper**

```matlab
% Windowed-source modal coefficients C{h} = Phi{h}' * B{h} * F(gv{h},:), cached by
% (srcResult, iWin, operator). Seed and kernel params do NOT invalidate it (Apply filters the
% whole field). Scalar operators only (F reduced to per-vertex magnitude).
function [C, gvAll] = i_apply_projection(st, ax, D, iWin, nV)
    C = {};  gvAll = [];
    key = sprintf('%s|%d-%d|%s', D.srcResult, iWin(1), iWin(end), i_atom_op(st));
    M = getappdata(0, 'DynamicsApplyCache');
    if ~isempty(M) && isstruct(M) && isfield(M,'key') && strcmp(M.key, key)
        C = M.C;  gvAll = M.gvAll;  return;
    end
    F = double(bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], iWin, 0));
    Fr = i_paintable_scalar(F, nV);                       % scalar per-vertex magnitude field
    C = cell(1, numel(ax.Phi));  gvAll = [];
    for h = 1:numel(ax.Phi)
        if isempty(ax.Phi{h}), continue; end
        gv = ax.GlobalVertices{h}(:);
        C{h} = manifold_ft(ax.Phi{h}, ax.Mass{h}, Fr(gv,:));
        gvAll = [gvAll; gv]; %#ok<AGROW>
    end
    setappdata(0, 'DynamicsApplyCache', struct('key',key, 'C',{C}, 'gvAll',gvAll));
end
```

- [ ] **Step 2: Refactor `i_atom_apply` to use the cache for STATIC kernels (dynamic falls back)**

Replace the reconstruct+filter block in `i_atom_apply` (the `Fr = i_paintable_scalar(F, nV)` → `Ffilt = i_atom_filter_field(...)` → reduce region) with a domain-gated flow. Static kernels use the cached fast path; dynamic (ts/js) kernels keep the existing `i_atom_filter_field` (joint time-vertex) path, because the cached spatial-gain apply is only correct for a time-invariant `g(λ)`:
```matlab
    % --- domain-gated apply: static g(lambda) uses the cached projection (instant on param drag);
    %     dynamic ts/js kernels keep the joint time-vertex path (cache does not apply). ---
    iWin = i_cursor_window(D.srcDS, D.srcResult, 4);
    if isempty(iWin), ctrl.jAtomInfo.setText('Apply: no recording window'); return; end
    meta = bst_eigfilter_kernel('info', kernel);
    dom  = 'static';  if isfield(meta,'domain') && ~isempty(meta.domain), dom = meta.domain; end
    bst_progress('start', 'Atom', 'Filtering the real source...');
    if strcmpi(dom, 'static')
        [C, ~] = i_apply_projection(st, ax, D, iWin, nV);         % cached by (srcResult,iWin,operator)
        if isempty(C), bst_progress('stop'); ctrl.jAtomInfo.setText('Apply: projection failed'); return; end
        g = bst_eigfilter_kernel(kernel, kp);
        Ffilt = zeros(nV, numel(iWin));
        for h = 1:numel(ax.Phi)
            if isempty(ax.Phi{h}) || isempty(C{h}), continue; end
            gv = ax.GlobalVertices{h}(:);  hgain = g(ax.Lambda{h}(:));
            Ffilt(gv,:) = manifold_ift(ax.Phi{h}, hgain(:) .* C{h});
        end
    else
        F  = double(bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], iWin, 0));
        Fr = i_paintable_scalar(F, nV);
        Ffilt = i_atom_filter_field(Fr, ax, variant, kernel, kp);   % dynamic: JTVAnalysis (unchanged)
    end
    bst_progress('stop');
```
Keep the existing scalar-operator guard ABOVE this block unchanged (Dirac/vector still bail with the "scalar-only for now" message), and keep the subsequent `if isempty(Ffilt) ... i_paintable_scalar(Ffilt,nV) ... normalize ... view_dynamics('SetFilteredField', ...)` lines. `variant`/`kernel`/`kp`/`nV`/`ax` are already computed earlier in `i_atom_apply` — reuse them. `i_atom_filter_field` remains used (dynamic branch), not orphaned.

- [ ] **Step 3: Invalidate the cache on operator/session change**

In `OnSetOperator` and `SetTarget`, add `setappdata(0, 'DynamicsApplyCache', []);` so a stale projection from a previous operator/session is never reused. (The window part of the key already guards cursor moves.)

- [ ] **Step 4: Add the cached-path equivalence headless test**

Create `dev/test_apply_cache_equiv.m` — asserts the cached apply-gain equals the direct `bst_eigenfilter('Analysis')` result on a synthetic scalar field + eigenbasis. Because this needs a real eigenbasis, mark it a CONTROLLER-run test (it loads an operator). Minimal form:
```matlab
function tests = test_apply_cache_equiv
tests = functiontests(localfunctions);
end
function test_cached_equals_analysis(tc)
% Controller runs this with a surface on the path. Build a small LB eigenbasis, a random
% scalar field, and verify C-then-gain == bst_eigenfilter('Analysis').
ax = bst_eigen('Axes', struct('SurfaceFile', getenv('BST_TEST_SURF'), 'Variant','Laplace-Beltrami', 'nModes',40));
nV = 0; for h=1:numel(ax.GlobalVertices), nV = max(nV, max(ax.GlobalVertices{h}(:))); end
F  = randn(nV, 5);
kp = struct('t', 0.02, 'lmax', max(ax.Lambda{1}(:)));
g  = bst_eigfilter_kernel('heat', kp);
% cached path
Ffilt = zeros(nV,5);
for h=1:numel(ax.Phi)
    gv = ax.GlobalVertices{h}(:); C = manifold_ft(ax.Phi{h}, ax.Mass{h}, F(gv,:));
    hg = g(ax.Lambda{h}(:)); Ffilt(gv,:) = manifold_ift(ax.Phi{h}, hg(:).*C);
end
% direct Analysis
EigenMat = struct('Phi',{ax.Phi},'Lambda',{ax.Lambda},'Variant','Laplace-Beltrami','GlobalVertices',{ax.GlobalVertices});
OperatorMat = struct('Mass',{ax.Mass});
[Fana,~,err] = bst_eigenfilter('Analysis', F, EigenMat, OperatorMat, 'heat', kp);
verifyEqual(tc, err, 0);
verifyLessThan(tc, max(abs(Ffilt(:)-Fana(:))), 1e-9);
end
```

- [ ] **Step 5: Static check + grep**

Run `check_matlab_code` on `panel_bst_dynamics.m` and `dev/test_apply_cache_equiv.m`. Expected: no errors.
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
grep -c "i_apply_projection\|DynamicsApplyCache" toolbox/gui/panel_bst_dynamics.m  # def + uses + 2 invalidations
```

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m dev/test_apply_cache_equiv.m
git commit -m "refactor(dynamics): cache windowed-source projection for instant Apply re-synthesis

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 6: controller live consolidated pass (headless tests + GUI smoke)

**Files:** none (verification only; fixes go to the relevant file if a regression surfaces).

- [ ] **Step 1: Run the headless tests** in the established MATLAB session (controller): `test_itersine_kernel` (tight-frame property), and `test_apply_cache_equiv` (set `BST_TEST_SURF` to sub-MTL0002's cortex; cached == Analysis).

- [ ] **Step 2: GUI smoke** (controller, via MCP, mirroring sub-project A's Task 6): open the sub-MTL0002 Dirac session; `OnCreateAtom`; `OnSetOperator('Laplace-Beltrami')`; set N=6; `OnDesignFrame` → 6 itersine atoms, `i_frame_refresh` labels show tightness ✓; `jFrameShow` on → coverage view renders flat `S(λ)` (screenshot). Add a mismatched hand kernel → tightness degrades (screenshot). Apply mode: drag a slider → instant re-paint from cache, no progress bar (screenshot). Cleanup: close session, do not save.

- [ ] **Step 3: Record** results in the ledger; fix any regression in-place and re-verify.

---

## Acceptance criteria (whole plan)

- `itersine` is a discoverable registry kernel whose N-member set is a tight frame (B/A→1 on the interior); it does NOT appear in the hand-pick Filter combobox.
- The Frame section shows live A/B/tightness; the coverage view (reused `view_eigfilter_response`) renders `S(λ)` + per-member tiles with click-to-select.
- "Design tight frame" replaces the bank with N itersine members (confirm if non-empty); tightness reads ✓.
- Apply re-synthesis on a slider drag is instant (cached projection); the cached result is numerically identical to `bst_eigenfilter('Analysis')`.
- `check_matlab_code` clean; existing `view_eigfilter_response` callers unaffected; controller live pass green.
