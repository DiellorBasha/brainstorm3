# Dirac Eigenmode Vector-Field Viewer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the deprecated scalar-LBO `view_eigenmodes` with a standalone Dirac eigenmode **vector-field** viewer launched from the `eigen_` DB node, and fully retire the `panel_eigenmodes` display-lever subsystem.

**Architecture:** Rewrite `toolbox/gui/view_eigenmodes.m` as a self-contained viewer (modeled on `view_leadfield_vectors`): load an `eigen_*.mat` node, dispatch on `Variant` (Dirac implemented; others `bst_error`), reconstruct each eigenvector's ambient 3-vector from its quaternion vector part, and draw a `quiver3` field on a translucent cortex with single-mode keyboard stepping. Then delete `panel_eigenmodes.m` and unhook its couplings in `gui_initialize`, `bst_figures`, and `panel_surface`; rewire the eigen-node menu; and delete the orphaned lever/old-viewer tests.

**Tech Stack:** MATLAB R2023b, Brainstorm dev fork. Live tests run inside a running Brainstorm session via the MATLAB MCP. Headless unit tests use the `eval(macro_method)` subfunction-dispatch pattern (`view_eigenmodes('Name', args)`).

**Spec:** `docs/superpowers/specs/2026-06-11-dirac-eigenmode-vector-viewer-design.md`

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `toolbox/gui/view_eigenmodes.m` | Standalone Dirac vector-field viewer + pure `ReconstructModeField` | **Rewrite** |
| `dev/tests/test_eigenmode_vector_field.m` | Headless unit tests for `ReconstructModeField` | **Create** |
| `dev/tests/test_eigenmode_vector_viewer.m` | Live-figure integration test | **Create** |
| `toolbox/tree/tree_callbacks.m` | Eigen-node "View" → new viewer; remove cortex "View eigenmodes" | **Modify** |
| `toolbox/process/functions/process_eigenmodes.m` | Drop the legacy auto-view call | **Modify** |
| `toolbox/gui/panel_eigenmodes.m` | The display lever | **Delete** |
| `toolbox/gui/gui_initialize.m` | Lever tab registration | **Modify** |
| `toolbox/core/bst_figures.m` | `FireModesChanged` + lever UpdatePanel hooks | **Modify** |
| `toolbox/gui/panel_surface.m` | Lever `ApplyToColumn`/`IsActive` display hooks | **Modify** |
| 13 `dev/tests/test_eigenmode_lever_*`, `test_eigenmode_viewer_*`, `test_view_eigenmodes_pure`, `test_eigfilter_design_*` | Tests of the retired lever/old viewer | **Delete** |

**Kept untouched:** `panel_eigenmodes_compute.m` (independent legacy compute dialog), `process_eigenmodes` legacy compute + its menu item, the entire `in_tess_eigenmodes` ecosystem, and the `toolbox/math/eigfilter/` library.

---

## Task 1: Rewrite view_eigenmodes.m as the standalone Dirac vector viewer

**Files:**
- Rewrite: `toolbox/gui/view_eigenmodes.m`
- Test: `dev/tests/test_eigenmode_vector_field.m`

**Context:** The current `view_eigenmodes(SurfaceFile)` loads deprecated surface-stored scalar modes and drives the `panel_eigenmodes` lever. The new contract is `view_eigenmodes(EigenFile)` — an `eigen_*.mat` node path. `EigenMat` (`db_template('eigenmat')`) has `Phi{1,2}` per hemisphere; for Dirac each `Phi{hh}(:,k)` is `[4·nVh × 1]` (quaternion `w,i,j,k` per vertex). The ambient 3-vector is rows `2:4` of each 4-block (drop the `w` slot at `1:4:end`), scattered into the full vertex grid via `GlobalVertices{hh}` — exactly as `bst_dirac/local_reconstruct` does. Dirac eigenvectors are real.

**CRITICAL — axes-reset trap:** `quiver3` is a high-level plot function; on the `'Axes3D'`-tagged axes (default `NextPlot='replace'`) it runs `newplot` and **resets the axes — wiping the tag and deleting the cortex patch**. `view_leadfield_vectors` avoids this with `hold on` after `view_surface`. This viewer **must** call `hold(hAxes,'on')` before any `quiver3`. (See memory: brainstorm-axes-nextplot-trap.)

- [ ] **Step 1: Write the failing headless test**

Create `dev/tests/test_eigenmode_vector_field.m`:

```matlab
function test_eigenmode_vector_field()
% TEST_EIGENMODE_VECTOR_FIELD  Headless regression for the Dirac eigenmode
% vector-field reconstruction (view_eigenmodes pure core). Requires Brainstorm
% on path so view_eigenmodes('ReconstructModeField', ...) dispatches.
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;

    % Synthetic 2-hemi EigenMat: nVh=2 per hemi (Phi{hh} is [4*2 x K]), K=3.
    % Global vertices: hemi L -> [1;2], hemi R -> [3;4]; nVert=4.
    % Per mode k, quaternion [w i j k] per vertex:
    %   L v1 = [0 k 0 0] (x=k), L v2 = [0 0 k 0] (y=k),
    %   R v3 = [0 0 0 k] (z=k), R v4 = [0 k k k] (k,k,k).
    EM.Variant = 'Dirac';
    EM.GlobalVertices = {[1;2], [3;4]};
    PhiL = zeros(8,3); PhiR = zeros(8,3);
    for k = 1:3
        PhiL(1:4,k) = [0;k;0;0];
        PhiL(5:8,k) = [0;0;k;0];
        PhiR(1:4,k) = [0;0;0;k];
        PhiR(5:8,k) = [0;k;k;k];
    end
    EM.Phi = {PhiL, PhiR};
    EM.Lambda = {(1:3)', (1:3)'};

    V3 = view_eigenmodes('ReconstructModeField', EM, 2, 4);
    [nPass,nFail] = chk('v1 = (2,0,0)', isequal(V3(1,:),[2 0 0]), nPass,nFail);
    [nPass,nFail] = chk('v2 = (0,2,0)', isequal(V3(2,:),[0 2 0]), nPass,nFail);
    [nPass,nFail] = chk('v3 = (0,0,2)', isequal(V3(3,:),[0 0 2]), nPass,nFail);
    [nPass,nFail] = chk('v4 = (2,2,2)', isequal(V3(4,:),[2 2 2]), nPass,nFail);

    % off-support vertices stay zero (nVert larger than mapped indices)
    V3b = view_eigenmodes('ReconstructModeField', EM, 1, 6);
    [nPass,nFail] = chk('off-support vertices zero', ...
        isequal(V3b(5,:),[0 0 0]) && isequal(V3b(6,:),[0 0 0]), nPass,nFail);

    % w slot is dropped: a large w must not change the vector part
    EM2 = EM; EM2.Phi{1}(1,2) = 99;   % w of L v1, mode 2
    V3c = view_eigenmodes('ReconstructModeField', EM2, 2, 4);
    [nPass,nFail] = chk('w slot ignored', isequal(V3c(1,:),[2 0 0]), nPass,nFail);

    % out-of-range mode errors
    err = false; try, view_eigenmodes('ReconstructModeField', EM, 9, 4); catch, err = true; end
    [nPass,nFail] = chk('out-of-range mode errors', err, nPass,nFail);

    fprintf('\n==== test_eigenmode_vector_field: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_eigenmode_vector_field: %d test(s) FAILED.', nFail); end
end

function [p,f] = chk(name, c, p, f)
    if c, fprintf('  PASS  %s\n', name); p=p+1; else, fprintf('  FAIL  %s\n', name); f=f+1; end
end
```

- [ ] **Step 2: Run the test to verify it fails**

In the MATLAB MCP session: `test_eigenmode_vector_field`
Expected: FAIL — the old `view_eigenmodes` dispatch list does not include `'ReconstructModeField'`, so it routes to `ViewFigure('ReconstructModeField')` and errors (or returns wrong type).

- [ ] **Step 3: Rewrite `toolbox/gui/view_eigenmodes.m`**

Replace the **entire file** with:

```matlab
function varargout = view_eigenmodes(varargin)
% VIEW_EIGENMODES: Standalone vector-field viewer for Dirac eigenmodes.
%
% USAGE:  hFig = view_eigenmodes(EigenFile)
%         V3   = view_eigenmodes('ReconstructModeField', EigenMat, k, nVert)
%
% Loads an eigen_ DB node (db_template('eigenmat')) and renders each Dirac
% eigenvector as an ambient 3D quiver field on the parent cortex. Cycle modes
% with the keyboard (Left/Right = -/+1, PgUp/PgDn = +/-10), like
% view_leadfield_vectors cycles sensor channels.
%
% Only the Dirac variant is implemented; LBO / Connection Laplacian raise a
% bst_error (added later as variant branches). Dirac eigenvectors are real
% [4*nVh x K] quaternion fields; the per-vertex 3-vector is the quaternion
% vector part (rows 2:4 of each 4-block; the w slot 1:4:end is dropped) --
% exactly as bst_dirac/local_reconstruct extracts the ambient field.
%
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
%
% Authors: Diellor Basha, 2026

if (nargin >= 1) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'ReconstructModeField'}))
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: per-vertex ambient 3-vector for Dirac eigenmode k =====
function V3 = ReconstructModeField(EigenMat, k, nVert)
% V3 [nVert x 3]: zeros off-support; quaternion vector part (i,j,k)->(x,y,z),
% w slot dropped; scattered to global vertices via EigenMat.GlobalVertices.
    if ~isfield(EigenMat,'Phi') || isempty(EigenMat.Phi) || numel(EigenMat.Phi) ~= 2
        error('view_eigenmodes:badEigen', 'EigenMat.Phi must be a 1x2 per-hemisphere cell.');
    end
    V3 = zeros(nVert, 3);
    for hh = 1:2
        vH  = EigenMat.GlobalVertices{hh}(:);
        Phi = EigenMat.Phi{hh};
        if (k < 1) || (k > size(Phi,2))
            error('view_eigenmodes:badMode', 'Mode index %d out of range 1..%d.', k, size(Phi,2));
        end
        if size(Phi,1) ~= 4*numel(vH)
            error('view_eigenmodes:shapeMismatch', ...
                'Hemisphere %d: Phi has %d rows, expected 4*nV=%d.', hh, size(Phi,1), 4*numel(vH));
        end
        col = double(Phi(:, k));
        V3(vH,1) = col(2:4:end);   % quaternion i -> x
        V3(vH,2) = col(3:4:end);   % quaternion j -> y
        V3(vH,3) = col(4:4:end);   % quaternion k -> z
    end
end


%% ===== GUI: standalone quiver viewer =====
function hFig = ViewFigure(EigenFile)
    hFig = [];
    % --- load + validate eigen node ---
    EigenFull = file_fullpath(EigenFile);
    if ~file_exist(EigenFull)
        bst_error('Eigen file not found.', 'View eigenmodes', 0);
        return;
    end
    EigenMat = load(EigenFull);
    if ~isfield(EigenMat,'Variant') || isempty(EigenMat.Variant)
        bst_error('Eigen file has no Variant field.', 'View eigenmodes', 0);
        return;
    end
    % --- variant dispatch (Dirac implemented; others deferred) ---
    switch lower(EigenMat.Variant)
        case 'dirac'
            % implemented below
        case {'laplace-beltrami','connection laplacian'}
            bst_error(sprintf(['Vector viewer currently supports Dirac eigenmodes only.' 10 ...
                'This is a "%s" node.'], EigenMat.Variant), 'View eigenmodes', 0);
            return;
        otherwise
            bst_error(sprintf('Unknown eigen variant: %s', EigenMat.Variant), 'View eigenmodes', 0);
            return;
    end
    if ~isfield(EigenMat,'Phi') || isempty(EigenMat.Phi) || numel(EigenMat.Phi) ~= 2 ...
            || isempty(EigenMat.Phi{1}) || isempty(EigenMat.Phi{2})
        bst_error('Dirac eigen node has empty Phi.', 'View eigenmodes', 0);
        return;
    end

    Surface  = EigenMat.ParentSurface;
    TessMat  = in_tess_bst(Surface);
    Vertices = TessMat.Vertices;
    nVert    = size(Vertices, 1);
    K        = min(size(EigenMat.Phi{1},2), size(EigenMat.Phi{2},2));
    Tau      = NaN;
    if isfield(EigenMat,'Provenance') && isstruct(EigenMat.Provenance) ...
            && isfield(EigenMat.Provenance,'Tau') && ~isempty(EigenMat.Provenance.Tau)
        Tau = EigenMat.Provenance.Tau;
    end

    % --- display cortex (translucent gray) ---
    hFig = view_surface(Surface, 0.5, [0.5 0.5 0.5], 'NewFigure');
    if isempty(hFig)
        bst_error('Could not open the surface figure.', 'View eigenmodes', 0);
        return;
    end
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    % CRITICAL: hold the axes so quiver3 (high-level) does not run newplot and
    % reset the 'Axes3D' axes (which would delete the cortex patch).
    hold(hAxes, 'on');
    figure_3d('SetStandardView', hFig, 'left');
    set(hFig, 'Name', ['Eigenmodes: ' EigenMat.Variant ' | ' Surface]);

    % --- state (closure vars) ---
    iMode              = 1;
    quiverSize         = 1;
    quiverWidth        = 1;
    thresholdAmplitude = 1;     % fraction of cumulative norm kept (1 = all)
    thresholdBalance   = 0;     % 0 = keep small (<=), 1 = keep large (>)
    useNormalize       = false;

    % --- legend ---
    hLabel = uicontrol('Style','text','String','...','Units','Pixels', ...
        'Position',[6 1 1600 35],'HorizontalAlignment','left', ...
        'FontUnits','points','FontSize',bst_get('FigFont'), ...
        'ForegroundColor',[.3 1 .3],'BackgroundColor',[0 0 0],'Parent',hFig);

    % --- keyboard ---
    KeyPressFcn_bak = get(hFig, 'KeyPressFcn');
    set(hFig, 'KeyPressFcn', @KeyPress_Callback);

    DrawArrows();

    % ===== NESTED: draw the current mode's field =====
    function DrawArrows()
        delete(findobj(hAxes, '-depth', 1, 'Tag', 'eigArrows'));
        V3 = ReconstructModeField(EigenMat, iMode, nVert);
        if useNormalize
            nv = sqrt(sum(V3.^2, 2));
            nz = nv > eps;
            V3(nz,:) = V3(nz,:) ./ nv(nz);
        end
        % cumulative-norm amplitude gate (like view_leadfield_vectors)
        normV = sqrt(sum(V3.^2, 2));
        [sv, ind] = sort(normV, 'ascend');
        cdf = cumsum(sv);
        if cdf(end) > 0, cdf = cdf / cdf(end); end
        if thresholdBalance == 0
            keep = find(cdf <= thresholdAmplitude);
        else
            keep = find(cdf > thresholdAmplitude);
        end
        Vre = Vertices(ind, :);
        Dre = zeros(numel(ind), 3);
        Dre(keep, :) = V3(ind(keep), :);
        quiver3(Vre(:,1), Vre(:,2), Vre(:,3), Dre(:,1), Dre(:,2), Dre(:,3), quiverSize, ...
            'Parent', hAxes, 'LineWidth', quiverWidth, 'Color', [.3 1 .3], 'Tag', 'eigArrows');
        % legend
        lamL = EigenMat.Lambda{1}(iMode);
        lamR = EigenMat.Lambda{2}(iMode);
        tauStr = ''; if ~isnan(Tau), tauStr = sprintf(' | tau=%.3g', Tau); end
        normStr = ''; if useNormalize, normStr = ' | unit'; end
        set(hLabel, 'String', sprintf(['Mode %d / %d   |   lambdaL=%.4g, lambdaR=%.4g%s%s   ' ...
            '[arrows: %d | H for help]'], iMode, K, lamL, lamR, tauStr, normStr, numel(keep)));
    end

    % ===== NESTED: keyboard navigation =====
    function KeyPress_Callback(h, keyEvent)
        switch (keyEvent.Key)
            case 'leftarrow'
                if     ismember('shift',   keyEvent.Modifier), quiverSize  = quiverSize  / 1.2;
                elseif ismember('control', keyEvent.Modifier), quiverWidth = quiverWidth / 1.2;
                elseif ismember('alt',     keyEvent.Modifier), thresholdAmplitude = thresholdAmplitude - 0.01;
                else,  iMode = iMode - 1; end
            case 'rightarrow'
                if     ismember('shift',   keyEvent.Modifier), quiverSize  = quiverSize  * 1.2;
                elseif ismember('control', keyEvent.Modifier), quiverWidth = quiverWidth * 1.2;
                elseif ismember('alt',     keyEvent.Modifier), thresholdAmplitude = thresholdAmplitude + 0.01;
                else,  iMode = iMode + 1; end
            case 'uparrow'
                if     ismember('shift',   keyEvent.Modifier), quiverSize  = quiverSize  * 1.2;
                elseif ismember('control', keyEvent.Modifier), quiverWidth = quiverWidth * 1.2;
                elseif ismember('alt',     keyEvent.Modifier), thresholdAmplitude = thresholdAmplitude + 0.01;
                else,  return; end
            case 'downarrow'
                if     ismember('shift',   keyEvent.Modifier), quiverSize  = quiverSize  / 1.2;
                elseif ismember('control', keyEvent.Modifier), quiverWidth = quiverWidth / 1.2;
                elseif ismember('alt',     keyEvent.Modifier), thresholdAmplitude = thresholdAmplitude - 0.01;
                else,  return; end
            case 'pageup',   iMode = iMode + 10;
            case 'pagedown', iMode = iMode - 10;
            case 'n',        useNormalize = ~useNormalize;
            case 'return'
                if ismember('alt', keyEvent.Modifier), thresholdBalance = ~thresholdBalance; else, return; end
            case 'h'
                java_dialog('msgbox', ['<HTML><TABLE>' ...
                    '<TR><TD><B>Left/Right</B></TD><TD>Previous/next mode</TD></TR>' ...
                    '<TR><TD><B>PgUp/PgDn</B></TD><TD>+/- 10 modes</TD></TR>' ...
                    '<TR><TD><B>Shift+Left/Right</B></TD><TD>Arrow length -/+</TD></TR>' ...
                    '<TR><TD><B>Control+Left/Right</B></TD><TD>Arrow width -/+</TD></TR>' ...
                    '<TR><TD><B>Alt+Left/Right</B></TD><TD>Amplitude threshold -/+</TD></TR>' ...
                    '<TR><TD><B>Alt+Enter</B></TD><TD>Toggle threshold direction</TD></TR>' ...
                    '<TR><TD><B>N</B></TD><TD>Toggle unit-normalized arrows</TD></TR>' ...
                    '<TR><TD><B>0-9</B></TD><TD>Change view</TD></TR>' ...
                    '</TABLE>'], 'Keyboard shortcuts', [], 0);
                return;
            otherwise
                if ~isempty(KeyPressFcn_bak), KeyPressFcn_bak(h, keyEvent); end
                return;
        end
        if iMode < 1, iMode = K; end
        if iMode > K, iMode = 1; end
        if thresholdAmplitude < 0, thresholdAmplitude = 0; end
        if thresholdAmplitude > 1, thresholdAmplitude = 1; end
        DrawArrows();
    end
end
```

- [ ] **Step 4: Run the headless test to verify it passes**

In the MATLAB MCP session: `rehash; test_eigenmode_vector_field`
Expected: `==== test_eigenmode_vector_field: 7 passed, 0 failed ====`

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/view_eigenmodes.m dev/tests/test_eigenmode_vector_field.m
git commit -m "feat(eigen): rewrite view_eigenmodes as standalone Dirac vector-field viewer

New contract view_eigenmodes(EigenFile): loads an eigen_ node, dispatches on
Variant (Dirac implemented; LBO/Connection -> bst_error), reconstructs each
eigenvector's ambient 3-vector from the quaternion vector part, and draws a
quiver3 field on a translucent cortex with single-mode keyboard stepping.
Holds the Axes3D (NextPlot='add') so quiver3 does not reset the cortex axes.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 2: Wire the eigen_ node to the viewer; remove the deprecated cortex menu + auto-view

**Files:**
- Modify: `toolbox/tree/tree_callbacks.m:4081-4102` (`EigenView_Callback`)
- Modify: `toolbox/tree/tree_callbacks.m:1185` (cortex "View eigenmodes" menu item)
- Modify: `toolbox/process/functions/process_eigenmodes.m:287-289` (legacy auto-view)

**Context:** The eigen-node popup already has a "View" item calling `EigenView_Callback(filenameFull)` (`tree_callbacks.m:2579`), but the callback is a stub that prints field names. Repoint it at the new viewer. Remove the deprecated cortex-node "View eigenmodes" item (it used the surface-stored scalar path). Decouple the legacy `process_eigenmodes` from the viewer (its old `view_eigenmodes(SurfaceFile)` call is invalid under the new contract).

- [ ] **Step 1: Replace the `EigenView_Callback` body**

In `toolbox/tree/tree_callbacks.m`, replace the function at ~4081-4102:

```matlab
function EigenView_Callback(filenameFull)
    if ~file_exist(filenameFull)
        bst_error('Eigen file not found.', 'View eigen', 0);
        return;
    end
    view_eigenmodes(filenameFull);
end
```

- [ ] **Step 2: Remove the deprecated cortex "View eigenmodes" menu item**

In `toolbox/tree/tree_callbacks.m`, delete line ~1185:

```matlab
                    gui_component('MenuItem', jPopup, [], 'View eigenmodes', IconLoader.ICON_RESULTS, [], @(h,ev)bst_call(@view_eigenmodes, filenameRelative));
```

Leave the adjacent "Compute eigenmodes (legacy)" (~1183) and "View connection phase" (~1186) items intact.

- [ ] **Step 3: Remove the legacy auto-view in process_eigenmodes**

In `toolbox/process/functions/process_eigenmodes.m`, delete the two lines at ~287-289:

```matlab
    % Visual confirmation: open the viewer on the freshly-computed modes
    view_eigenmodes(SurfaceFile);
```

(Keep the preceding `db_reload_subjects(iSubject);`. The `end` that closes the function stays.)

- [ ] **Step 4: Verify no stale references to the old viewer contract**

Run:
```bash
grep -rn "view_eigenmodes(" toolbox/ | grep -v "view_eigenmodes.m"
```
Expected: only the new callsite `EigenView_Callback` → `view_eigenmodes(filenameFull)`. No `view_eigenmodes(SurfaceFile)` / `view_eigenmodes(filenameRelative)` callers remain.

- [ ] **Step 5: Commit**

```bash
git add toolbox/tree/tree_callbacks.m toolbox/process/functions/process_eigenmodes.m
git commit -m "feat(eigen): launch the Dirac vector viewer from the eigen_ node

EigenView_Callback now opens view_eigenmodes(EigenFile). Remove the deprecated
cortex-node 'View eigenmodes' item (surface-stored scalar path) and the legacy
process_eigenmodes auto-view call that used the old contract.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 3: Retire the panel_eigenmodes display-lever subsystem

**Files:**
- Delete: `toolbox/gui/panel_eigenmodes.m`
- Modify: `toolbox/gui/gui_initialize.m:53`
- Modify: `toolbox/core/bst_figures.m` (lines ~890-892, ~997-1040, ~1163-1165)
- Modify: `toolbox/gui/panel_surface.m` (lines ~1815-1817, ~1857-1862)

**Context:** The lever exists only to drive the old scalar viewer. Its only `GlobalData.UserModes` co-writer is `panel_eigenmodes.m` itself, and `FireModesChanged` is only invoked from `panel_eigenmodes.m:828`, so deleting the panel makes `FireModesChanged` dead. The `panel_surface` hooks (`ApplyToColumn`, `IsActive`) revert to plain display. `panel_eigenmodes_compute.m` is a **separate** legacy compute dialog (no calls into the lever) and is **kept**.

- [ ] **Step 1: Delete the lever panel**

```bash
git rm toolbox/gui/panel_eigenmodes.m
```

- [ ] **Step 2: Remove the lever tab registration**

In `toolbox/gui/gui_initialize.m`, delete line 53:

```matlab
gui_show('panel_eigenmodes', 'BrainstormTab', 'tools');
```

- [ ] **Step 3: Remove the `DeleteFigure` UpdatePanel hook in bst_figures.m**

At ~890-892, delete:

```matlab
        if gui_brainstorm('isTabVisible', 'EigenModes')
            panel_eigenmodes('UpdatePanel');
        end
```

- [ ] **Step 4: Remove the `SetCurrentFigure` UpdatePanel hook in bst_figures.m**

At ~1163-1165, delete:

```matlab
                if gui_brainstorm('isTabVisible', 'EigenModes')
                    panel_eigenmodes('UpdatePanel', hFig);
                end
```

- [ ] **Step 5: Delete the `FireModesChanged` function in bst_figures.m**

Delete the whole block at ~995-1040 (header comment through the function's closing `end`):

```matlab
%% ===== FIRE MODES CHANGED =====
% Repaint visible 3D source figures after the eigenmode lever changes.
function FireModesChanged() %#ok<DEFNU>
    ...
end
```

(The next function `FireCurrentFreqChanged` and its `%% ===== FIRE CURRENT FREQUENCY CHANGED =====` banner remain.)

- [ ] **Step 6: Remove the `ApplyToColumn` hook in panel_surface.m**

At ~1813-1822, replace:

```matlab
                TessInfo(iTess).Data = bst_memory('GetResultsValues', iDS, iResult, [], 'CurrentTimeIndex');
                % Eigenmode scale lever: live, non-destructive band-limited
                % reconstruction of the displayed column (no-op when the lever is
                % inactive or set for a different surface).
                TessInfo(iTess).Data = panel_eigenmodes('ApplyToColumn', ...
                    TessInfo(iTess).SurfaceFile, TessInfo(iTess).Data);
                if isempty(TessInfo(iTess).Data)
                    isOk = 0;
                    return;
                end
```

with:

```matlab
                TessInfo(iTess).Data = bst_memory('GetResultsValues', iDS, iResult, [], 'CurrentTimeIndex');
                if isempty(TessInfo(iTess).Data)
                    isOk = 0;
                    return;
                end
```

- [ ] **Step 7: Remove the `IsActive` DataMinMax override in panel_surface.m**

At ~1857-1862, delete:

```matlab
                % Eigenmode lever: when active, scale the colormap to the FILTERED
                % column's range -- the raw DataMinMax would compress a narrow-band
                % reconstruction to a near-uniform map.
                if ~isempty(TessInfo(iTess).Data) && panel_eigenmodes('IsActive', TessInfo(iTess).SurfaceFile)
                    TessInfo(iTess).DataMinMax = [min(TessInfo(iTess).Data(:)), max(TessInfo(iTess).Data(:))];
                end
```

- [ ] **Step 8: Verify no dangling references**

Run:
```bash
grep -rn "panel_eigenmodes(" toolbox/ ; echo "---" ; grep -rn "FireModesChanged\|UserModes\|ModesChangedCallback\|EigenView'" toolbox/
```
Expected: first grep empty (note: `panel_eigenmodes_compute` does NOT match `panel_eigenmodes(`); second grep empty. (`getappdata(...,'EigenView')` is gone with `FireModesChanged`.)

- [ ] **Step 9: Smoke-test Brainstorm startup**

In the MATLAB MCP session: `brainstorm stop; brainstorm start`
Expected: starts cleanly, no error about `panel_eigenmodes`, no "EigenModes" tab in the Tools area.

- [ ] **Step 10: Commit**

```bash
git add -A toolbox/gui/panel_eigenmodes.m toolbox/gui/gui_initialize.m toolbox/core/bst_figures.m toolbox/gui/panel_surface.m
git commit -m "refactor(eigen): retire the panel_eigenmodes display-lever subsystem

Delete panel_eigenmodes.m and unhook it: drop the Tools-tab registration
(gui_initialize), the two isTabVisible('EigenModes') UpdatePanel hooks and the
now-dead FireModesChanged (bst_figures), and the ApplyToColumn/IsActive display
hooks (panel_surface). panel_eigenmodes_compute (legacy compute dialog) and the
in_tess_eigenmodes ecosystem are untouched.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 4: Delete the orphaned lever / old-viewer tests

**Files (delete):**
- `dev/tests/test_eigenmode_lever_e2e.m`
- `dev/tests/test_eigenmode_lever_integration.m`
- `dev/tests/test_eigenmode_lever_lifecycle.m`
- `dev/tests/test_eigenmode_lever_paired.m`
- `dev/tests/test_eigenmode_lever_panel.m`
- `dev/tests/test_eigenmode_lever_state.m`
- `dev/tests/test_eigenmode_lever_weights.m`
- `dev/tests/test_eigenmode_panel_centerwidth.m`
- `dev/tests/test_eigenmode_viewer_e2e.m`
- `dev/tests/test_eigenmode_viewer_synth.m`
- `dev/tests/test_view_eigenmodes_pure.m`
- `dev/tests/test_eigfilter_design_smoke.m`
- `dev/tests/test_eigfilter_design_pure.m`

**Context:** Every file above exercises the deleted lever (`panel_eigenmodes('ResetState'/'GetWeights'/'KernelPairedWeights'/'GetDisplayColumn'/...)`) or the old `view_eigenmodes` contract (`ModesChangedCallback`, paired-grid synthesis). The `test_eigfilter_design_*` pair tests the lever's **kernel-design** API (heat/dog band weights), not the separate living `toolbox/math/eigfilter/` library, so they are part of the lever retirement.

- [ ] **Step 1: Confirm each file targets only retired code**

Run:
```bash
for f in test_eigenmode_lever_e2e test_eigenmode_lever_integration test_eigenmode_lever_lifecycle \
         test_eigenmode_lever_paired test_eigenmode_lever_panel test_eigenmode_lever_state \
         test_eigenmode_lever_weights test_eigenmode_panel_centerwidth test_eigenmode_viewer_e2e \
         test_eigenmode_viewer_synth test_view_eigenmodes_pure test_eigfilter_design_smoke \
         test_eigfilter_design_pure; do
  echo "== $f =="; grep -l "panel_eigenmodes\|ModesChangedCallback\|BuildPairedGrid\|SynthColumn" dev/tests/$f.m;
done
```
Expected: each filename echoes and matches (every file references retired symbols).

- [ ] **Step 2: Delete the files**

```bash
git rm dev/tests/test_eigenmode_lever_e2e.m dev/tests/test_eigenmode_lever_integration.m \
       dev/tests/test_eigenmode_lever_lifecycle.m dev/tests/test_eigenmode_lever_paired.m \
       dev/tests/test_eigenmode_lever_panel.m dev/tests/test_eigenmode_lever_state.m \
       dev/tests/test_eigenmode_lever_weights.m dev/tests/test_eigenmode_panel_centerwidth.m \
       dev/tests/test_eigenmode_viewer_e2e.m dev/tests/test_eigenmode_viewer_synth.m \
       dev/tests/test_view_eigenmodes_pure.m dev/tests/test_eigfilter_design_smoke.m \
       dev/tests/test_eigfilter_design_pure.m
```

- [ ] **Step 3: Verify no remaining test references the lever/old contract**

Run:
```bash
grep -rln "panel_eigenmodes\|ModesChangedCallback\|BuildPairedGrid\|SynthColumn" dev/tests/
```
Expected: empty.

- [ ] **Step 4: Commit**

```bash
git commit -m "test(eigen): delete orphaned panel_eigenmodes lever + old-viewer tests

Remove the 13 tests that exercised the retired display lever (incl. its kernel-
design API) and the old surface-based view_eigenmodes contract. The living
toolbox/math/eigfilter library is unaffected.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 5: Live integration test for the vector viewer

**Files:**
- Create: `dev/tests/test_eigenmode_vector_viewer.m`

**Context:** Requires a running Brainstorm session with a loaded protocol. Builds the fixture from the canonical cortex (`bst_canonical_cortex`) and a freshly computed Dirac eigen node (`tess_eigen(...,'Dirac',...)`), resolves the new node's filename by diffing the surface's `Eigen` list (the pattern used in `test_tess_eigen.m`), opens the viewer, and asserts the cortex survives, arrows draw, mode stepping updates the field, and a non-Dirac node is rejected.

- [ ] **Step 1: Write the integration test**

Create `dev/tests/test_eigenmode_vector_viewer.m`:

```matlab
function test_eigenmode_vector_viewer()
% TEST_EIGENMODE_VECTOR_VIEWER  Live-figure regression for the Dirac eigenmode
% vector viewer. Requires Brainstorm running with a loaded protocol. Builds a
% canonical cortex + Dirac eigen node, opens view_eigenmodes, and asserts the
% cortex survives the quiver draw, arrows update on mode step, and a non-Dirac
% node is rejected. Created eigen/operator nodes are left in place (cheap).
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;

    % --- fixture: canonical cortex + Dirac eigen node ---
    SurfaceFile = bst_canonical_cortex(20484);
    preEig  = local_eigen_names(SurfaceFile);
    tess_eigen(SurfaceFile, 'Dirac', 'K', 40, 'Tau', 0.5);   % saves + registers
    postEig = local_eigen_names(SurfaceFile);
    newEig  = setdiff(postEig, preEig);
    assert(~isempty(newEig), 'No new Dirac eigen node was created.');
    EigenFile = newEig{1};

    % --- open the viewer ---
    close(findobj(0, 'type', 'figure', 'Tag', '3DViz'));
    hFig = view_eigenmodes(EigenFile);
    drawnow;
    [nPass,nFail] = chk('viewer returns a figure', ~isempty(hFig) && ishandle(hFig), nPass,nFail);

    hAx3D = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    hPatch = findobj(hAx3D, 'Type', 'patch');
    hQ = findobj(hFig, 'Tag', 'eigArrows');
    [nPass,nFail] = chk('Axes3D survives quiver draw', ~isempty(hAx3D), nPass,nFail);
    [nPass,nFail] = chk('cortex patch survives quiver draw', ~isempty(hPatch), nPass,nFail);
    [nPass,nFail] = chk('arrows drawn', ~isempty(hQ) && numel(get(hQ(1),'UData')) > 0, nPass,nFail);

    % --- mode stepping updates the field ---
    U1 = get(findobj(hFig,'Tag','eigArrows'),'UData');
    KeyOnFig(hFig, 'rightarrow');  drawnow;
    U2 = get(findobj(hFig,'Tag','eigArrows'),'UData');
    [nPass,nFail] = chk('mode step changes the field', ~isequal(U1, U2), nPass,nFail);

    % --- quiver-size key changes arrow length (UData scales) ---
    Ua = get(findobj(hFig,'Tag','eigArrows'),'UData');
    KeyOnFig(hFig, 'rightarrow', {'shift'});  drawnow;   % Shift+Right = longer
    Ub = get(findobj(hFig,'Tag','eigArrows'),'UData');
    [nPass,nFail] = chk('quiver-size key rescales arrows', ~isequal(Ua, Ub), nPass,nFail);

    close(hFig);

    % --- non-Dirac node is rejected ---
    preL  = local_eigen_names(SurfaceFile);
    tess_eigen(SurfaceFile, 'Laplace-Beltrami', 'K', 40);
    postL = local_eigen_names(SurfaceFile);
    lboNew = setdiff(postL, preL);
    rejected = false;
    if ~isempty(lboNew)
        try
            hbad = view_eigenmodes(lboNew{1});
            rejected = isempty(hbad);   % bst_error path returns [] without a figure
            if ~isempty(hbad) && ishandle(hbad), close(hbad); end
        catch
            rejected = true;
        end
    end
    [nPass,nFail] = chk('LBO node rejected (Dirac-only)', rejected, nPass,nFail);

    fprintf('\n==== test_eigenmode_vector_viewer: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_eigenmode_vector_viewer: %d test(s) FAILED.', nFail); end
end

function names = local_eigen_names(SurfaceFile)
    names = {};
    [sSubject, ~, iSurface] = bst_get('SurfaceFile', SurfaceFile);
    if ~isempty(sSubject) && ~isempty(iSurface) ...
            && isfield(sSubject.Surface(iSurface), 'Eigen') ...
            && ~isempty(sSubject.Surface(iSurface).Eigen)
        names = {sSubject.Surface(iSurface).Eigen.FileName};
    end
end

function KeyOnFig(hFig, keyName, modifier)
    if nargin < 3, modifier = {}; end
    ev.Key = keyName; ev.Character = ''; ev.Modifier = modifier;
    cb = get(hFig, 'KeyPressFcn');
    cb(hFig, ev);
end

function [p,f] = chk(name, c, p, f)
    if c, fprintf('  PASS  %s\n', name); p=p+1; else, fprintf('  FAIL  %s\n', name); f=f+1; end
end
```

- [ ] **Step 2: Run the integration test**

In the MATLAB MCP session (Brainstorm running, a protocol with anatomy loaded): `test_eigenmode_vector_viewer`
Expected: `==== test_eigenmode_vector_viewer: 7 passed, 0 failed ====`

- [ ] **Step 3: Re-run the headless test (regression guard)**

`test_eigenmode_vector_field`
Expected: `7 passed, 0 failed`.

- [ ] **Step 4: Final retirement verification**

```bash
grep -rn "panel_eigenmodes(" toolbox/ ; echo "--- next two should be empty ---" ; \
grep -rn "FireModesChanged\|ModesChangedCallback" toolbox/ ; \
grep -rln "panel_eigenmodes\|ModesChangedCallback\|BuildPairedGrid" dev/tests/
```
Expected: all empty.

- [ ] **Step 5: Commit**

```bash
git add dev/tests/test_eigenmode_vector_viewer.m
git commit -m "test(eigen): live integration test for the Dirac vector viewer

Builds a canonical-cortex Dirac eigen node, opens view_eigenmodes, asserts the
cortex/Axes3D survive the quiver draw, mode stepping and the size key update the
field, and a non-Dirac (LBO) node is rejected by the variant dispatch.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Self-Review notes

- **Spec coverage:** §4 viewer rewrite → Task 1; §4 wiring (eigen node / cortex menu / process_eigenmodes) → Task 2; §4 lever retirement → Task 3; §8.1 headless → Task 1; §8.2 live → Task 5; §8.3 retirement greps → Tasks 3 & 5. Orphaned-test removal (implied by retirement) → Task 4.
- **Type/name consistency:** `ReconstructModeField(EigenMat, k, nVert) -> [nVert×3]` is defined in Task 1 and called identically in the Task 1 test and the Task 5 fixture math. `EigenFile` (node path) is the viewer's single argument throughout. Quiver tag `'eigArrows'` is consistent across the viewer and both tests.
- **Quaternion mapping** matches `bst_dirac/local_reconstruct` exactly: x=`2:4:end`, y=`3:4:end`, z=`4:4:end`, w=`1:4:end` dropped.
- **Axes-reset trap** explicitly handled (`hold(hAxes,'on')` before `quiver3`).
