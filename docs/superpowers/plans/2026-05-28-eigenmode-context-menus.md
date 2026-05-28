# Eigenmode Context Menus (Compute + Viewer) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add right-click cortex → *Compute eigenmodes* and → *View eigenmodes* (a leadfield-style transient viewer that steps modes with ←/→ and updates the cortex colormap live).

**Architecture:** A new `view_eigenmodes.m` mirrors `view_leadfield_sensitivity.m`: open a transient `view_surface` figure, register a `source` colormap, write each mode's per-vertex values into `TessInfo(1).Data`, and refresh via `panel_surface('UpdateSurfaceColormap')`; a `KeyPressFcn` steps the mode index. The viewer's index/data math lives in pure, headlessly-testable subfunctions reached through the standard Brainstorm macro dispatch. Compute reuses the **existing** `process_eigenmodes('Compute', …)` core via a thin new `ComputeInteractive` method. Two menu items are wired into the cortex branch of `tree_callbacks.m`.

**Tech Stack:** MATLAB, Brainstorm GUI (`tree_callbacks`, `view_surface`, `panel_surface`, `bst_colormaps`, `java_dialog`, `gui_component`), MATLAB MCP for static analysis + interactive validation.

---

## Reconciliation with the spec (read first)

The spec assumed the `Run` body of `process_eigenmodes.m` needed extracting into a shared `Compute` core. **That core already exists** — `process_eigenmodes.m:186` is `Compute(SurfaceFile, nModes, MassType, RemoveDC, Repair, isInteractive)`, already called by `Run` and already prompting interactively on non-manifold input. So:
- **No `Run` refactor is performed** (it would be a no-op).
- The spec's "compute core round-trip" test is **dropped as redundant** — the compute path is already covered by `dev/tests/test_eigenmodes_manifold_gate.m`, `test_io_eigenmodes_roundtrip.m`, and `test_process_eigenmodes_options.m`. Compute is validated end-to-end by the user in Task 5.

Everything else in the spec stands.

## File Structure

- **Create** `toolbox/gui/view_eigenmodes.m` — the viewer. Macro-dispatch header; pure subfunctions `GetModeDisplay`, `StepMode`; GUI subfunction `ViewFigure` with nested `UpdateMode` / `KeyPress_Callback`. One file, one responsibility (browse modes on a surface).
- **Create** `dev/tests/test_view_eigenmodes_pure.m` — headless tests for `GetModeDisplay` + `StepMode`.
- **Modify** `toolbox/process/functions/process_eigenmodes.m` — add `ComputeInteractive(iSubject, SurfaceFile)` (thin wrapper over existing `Compute`).
- **Modify** `toolbox/tree/tree_callbacks.m` — add two cortex-node menu items after the `Display` item (~line 1172).

## Conventions for this plan

- Tests follow the repo idiom: a function whose name equals the file, using `assert`, printing `ALL TESTS PASSED` at the end. Run via the MATLAB MCP `evaluate_matlab_code` calling the function name (NOT `run_matlab_test_file`).
- GUI code (Tasks 2–4) cannot be unit-tested; each such task is verified by (a) MATLAB MCP `check_matlab_code` (M-Lint, catches syntax/undefined-symbol errors headlessly) and (b) the interactive validation in Task 5.
- Commit after each task.

---

### Task 1: Viewer pure helpers (`GetModeDisplay`, `StepMode`) — TDD

**Files:**
- Create: `toolbox/gui/view_eigenmodes.m`
- Test: `dev/tests/test_view_eigenmodes_pure.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_view_eigenmodes_pure.m`:

```matlab
function test_view_eigenmodes_pure
% Verify the pure viewer helpers reached through the macro dispatch:
% GetModeDisplay (data column, symmetric color limits, eigenvalue/wavelength,
% label, index clamping) and StepMode (clamped stepping).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% ----- Fabricate a small Eigenmodes struct (as in_tess_eigenmodes returns) -----
rng(3);
nVert = 50; K = 8;
Eig = struct();
Eig.Vectors  = randn(nVert, K);
Eig.Values   = sort(abs(randn(K,1))) + 0.1;   % ascending, > 0
Eig.nModes   = K;
Eig.MassType = 'barycentric';

% ----- GetModeDisplay: basic correctness -----
d = view_eigenmodes('GetModeDisplay', Eig, 3);
assert(isequal(d.Data, Eig.Vectors(:,3)), 'Data must be the requested mode column.');
assert(d.iMode == 3 && d.nModes == K, 'iMode/nModes mismatch.');
m = max(abs(Eig.Vectors(:,3)));
assert(isequal(d.CLim, [-m, m]), 'CLim must be symmetric [-max|v|, +max|v|].');
assert(d.CLim(1) < d.CLim(2), 'CLim must be a non-degenerate increasing range.');
assert(abs(d.Lambda - Eig.Values(3)) < 1e-12, 'Lambda must be Values(iMode).');
assert(abs(d.Wavelength - 2*pi/sqrt(Eig.Values(3))) < 1e-9, 'Wavelength = 2*pi/sqrt(lambda).');
assert(ischar(d.Label) && ~isempty(d.Label), 'Label must be a non-empty string.');

% ----- GetModeDisplay: index clamping -----
dLo = view_eigenmodes('GetModeDisplay', Eig, -5);
assert(dLo.iMode == 1, 'iMode below 1 must clamp to 1.');
dHi = view_eigenmodes('GetModeDisplay', Eig, K+99);
assert(dHi.iMode == K, 'iMode above K must clamp to K.');

% ----- GetModeDisplay: degenerate (all-zero) mode guard -----
EigZ = Eig; EigZ.Vectors(:,2) = 0;
dz = view_eigenmodes('GetModeDisplay', EigZ, 2);
assert(dz.CLim(1) < dz.CLim(2), 'All-zero mode must still yield a non-degenerate CLim.');

% ----- GetModeDisplay: non-positive eigenvalue -> wavelength n/a -----
EigN = Eig; EigN.Values(1) = 0;
dn = view_eigenmodes('GetModeDisplay', EigN, 1);
assert(isnan(dn.Wavelength), 'lambda<=0 must give NaN wavelength.');
assert(~isempty(strfind(dn.Label, 'n/a')), 'Label must show n/a for lambda<=0.');

% ----- StepMode: stepping + clamping -----
assert(view_eigenmodes('StepMode', 3, +1, K) == 4, 'StepMode +1 failed.');
assert(view_eigenmodes('StepMode', 3, -1, K) == 2, 'StepMode -1 failed.');
assert(view_eigenmodes('StepMode', 1, -1, K) == 1, 'StepMode must clamp at 1.');
assert(view_eigenmodes('StepMode', K, +1, K) == K, 'StepMode must clamp at K.');
assert(view_eigenmodes('StepMode', K-2, +10, K) == K, 'StepMode +10 must clamp at K.');
assert(view_eigenmodes('StepMode', 3, -10, K) == 1, 'StepMode -10 must clamp at 1.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run via MATLAB MCP `evaluate_matlab_code`: `test_view_eigenmodes_pure`
Expected: FAIL — `Undefined function or variable 'view_eigenmodes'` (file does not exist yet).

- [ ] **Step 3: Create `view_eigenmodes.m` with dispatch + the two pure helpers**

Create `toolbox/gui/view_eigenmodes.m`:

```matlab
function varargout = view_eigenmodes(varargin)
% VIEW_EIGENMODES: Interactively browse Laplace-Beltrami eigenmodes on a surface.
%
% USAGE:  hFig  = view_eigenmodes(SurfaceFile)              % open viewer at mode 1
%         hFig  = view_eigenmodes(SurfaceFile, iMode)       % open viewer at mode iMode
%         disp  = view_eigenmodes('GetModeDisplay', Eig, iMode)     % pure helper (testable)
%         iNext = view_eigenmodes('StepMode', iMode, delta, nModes) % pure helper (testable)
%
% Modeled on view_leadfield_sensitivity.m: the Left/Right arrows step the
% displayed mode, the cortex colormap updates live, and a legend reports the
% mode index, eigenvalue, and approximate spatial wavelength. The figure is
% transient (no database node is created).
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

% ===== METHOD DISPATCH (headless-testable pure helpers) =====
if (nargin >= 1) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'GetModeDisplay', 'StepMode'}))
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
% ===== MAIN ENTRY: open the viewer =====
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: per-mode display package =====
function d = GetModeDisplay(Eig, iMode)
    nModes = size(Eig.Vectors, 2);
    iMode  = min(max(round(iMode), 1), max(nModes, 1));   % clamp to [1, nModes]
    data   = double(Eig.Vectors(:, iMode));
    m      = max(abs(data));
    if ~(m > 0)        % degenerate (all-zero) mode: avoid a zero-width colormap range
        m = 1e-30;
    end
    lambda = Eig.Values(iMode);
    if lambda > 0
        wavelength = 2*pi / sqrt(lambda);
        wlStr      = sprintf('%.3g', wavelength);
    else
        wavelength = NaN;
        wlStr      = 'n/a';
    end
    d = struct();
    d.iMode      = iMode;
    d.nModes     = nModes;
    d.Data       = data;
    d.CLim       = [-m, m];          % symmetric: eigenmodes are signed (+/- lobes)
    d.Lambda     = lambda;
    d.Wavelength = wavelength;
    d.Label      = sprintf('Mode %d / %d     lambda = %.4g     wavelength ~ %s', ...
                           iMode, nModes, lambda, wlStr);
end


%% ===== PURE: clamped mode stepping =====
function iNext = StepMode(iMode, delta, nModes)
    iNext = min(max(round(iMode) + round(delta), 1), max(nModes, 1));
end
```

(The `ViewFigure` subfunction is added in Task 2; the dispatch above only routes the two pure helpers, so this file is runnable for the test now.)

- [ ] **Step 4: Run test to verify it passes**

Run via MATLAB MCP `evaluate_matlab_code`: `test_view_eigenmodes_pure`
Expected: PASS — prints `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/view_eigenmodes.m dev/tests/test_view_eigenmodes_pure.m
git commit -m "Add view_eigenmodes pure helpers (GetModeDisplay, StepMode) + tests

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Viewer figure (`ViewFigure` + nested callbacks)

**Files:**
- Modify: `toolbox/gui/view_eigenmodes.m` (append the `ViewFigure` subfunction)

GUI code — no unit test. Verified by M-Lint here and by interactive validation in Task 5.

- [ ] **Step 1: Append the `ViewFigure` subfunction**

Add to the end of `toolbox/gui/view_eigenmodes.m` (after `StepMode`):

```matlab
%% ===== GUI: open the surface figure and wire up mode browsing =====
function hFig = ViewFigure(SurfaceFile, iMode)
    hFig = [];
    if (nargin < 2) || isempty(iMode)
        iMode = 1;
    end
    % Load eigenmodes embedded in the surface file
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed || isempty(Eig) || ~isfield(Eig, 'Vectors') || isempty(Eig.Vectors)
        bst_error(['No eigenmodes found on this surface.' 10 ...
                   'Right-click the cortex and run "Compute eigenmodes" first.'], ...
                   'View eigenmodes', 0);
        return;
    end
    nModes = size(Eig.Vectors, 2);
    iMode  = StepMode(iMode, 0, nModes);   % clamp into range

    % Open a transient surface figure (no database node)
    [hFig, iDS, iFig] = view_surface(SurfaceFile, 0, [], 'NewFigure', 0); %#ok<ASGLU>
    if isempty(hFig)
        bst_error('Could not open the surface figure.', 'View eigenmodes', 0);
        return;
    end
    % Register the source colormap so the overlay is colormapped (+ colorbar)
    bst_colormaps('AddColormapToFigure', hFig, 'source');
    % Figure name + default view
    set(hFig, 'Name', ['Eigenmodes: ' SurfaceFile]);
    figure_3d('SetStandardView', hFig, 'left');
    % Legend (bottom-left), mirrors view_leadfield_sensitivity
    hLabel = uicontrol('Style', 'text', 'String', '...', 'Units', 'Pixels', ...
        'Position', [6 0 520 20], 'HorizontalAlignment', 'left', ...
        'FontUnits', 'points', 'FontSize', bst_get('FigFont'), ...
        'ForegroundColor', [.9 .9 .9], 'BackgroundColor', [0 0 0], 'Parent', hFig);
    % Hook the keyboard callback (preserve the original for unhandled keys)
    KeyPressFcn_bak = get(hFig, 'KeyPressFcn');
    set(hFig, 'KeyPressFcn', @KeyPress_Callback);
    % Initial render
    UpdateMode();

    % ===== NESTED: render the current mode onto the surface =====
    function UpdateMode()
        d = GetModeDisplay(Eig, iMode);
        TessInfo = getappdata(hFig, 'Surface');
        TessInfo(1).DataSource.Type     = 'Source';
        TessInfo(1).DataSource.FileName = SurfaceFile;   % real DB file -> satisfies CLim guard
        TessInfo(1).Data                = d.Data;
        TessInfo(1).DataMinMax          = d.CLim;
        TessInfo(1).DataLimitValue      = d.CLim;
        setappdata(hFig, 'Surface', TessInfo);
        panel_surface('UpdateSurfaceColormap', hFig);
        set(hLabel, 'String', d.Label);
    end

    % ===== NESTED: keyboard navigation =====
    function KeyPress_Callback(h, keyEvent)
        switch (keyEvent.Key)
            case 'leftarrow'
                iMode = StepMode(iMode, -1, nModes);   UpdateMode();
            case 'rightarrow'
                iMode = StepMode(iMode, +1, nModes);   UpdateMode();
            case 'pageup'
                iMode = StepMode(iMode, +10, nModes);  UpdateMode();
            case 'pagedown'
                iMode = StepMode(iMode, -10, nModes);  UpdateMode();
            case 'h'
                java_dialog('msgbox', ...
                    ['Eigenmode viewer shortcuts:' 10 10 ...
                     '   Left / Right arrow  :  previous / next mode' 10 ...
                     '   Page Up / Page Down :  +/- 10 modes' 10 ...
                     '   H                   :  this help'], 'Eigenmode viewer');
            otherwise
                if ~isempty(KeyPressFcn_bak)
                    KeyPressFcn_bak(h, keyEvent);
                end
        end
    end
end
```

- [ ] **Step 2: Static-analysis check (headless)**

Run MATLAB MCP `check_matlab_code` on `toolbox/gui/view_eigenmodes.m`.
Expected: no errors/undefined symbols. (Brainstorm idioms like the unused `iDS`/`iFig` are suppressed with `%#ok<ASGLU>`; nested-function shared variables `iMode`, `hLabel`, `KeyPressFcn_bak` are expected.)

- [ ] **Step 3: Re-run the pure-helper test (no regression)**

Run via MATLAB MCP `evaluate_matlab_code`: `test_view_eigenmodes_pure`
Expected: PASS — `ALL TESTS PASSED` (dispatch still routes the pure helpers; appending `ViewFigure` must not break them).

- [ ] **Step 4: Commit**

```bash
git add toolbox/gui/view_eigenmodes.m
git commit -m "Add view_eigenmodes figure: arrow-stepped mode browser on the cortex

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: `ComputeInteractive` method in `process_eigenmodes.m`

**Files:**
- Modify: `toolbox/process/functions/process_eigenmodes.m` (add a subfunction; dispatch already exists via `eval(macro_method)` at line 43)

No unit test (interactive dialogs). The underlying `Compute` core is already covered by existing tests; this wrapper is validated interactively in Task 5.

- [ ] **Step 1: Add the `ComputeInteractive` subfunction**

Insert into `toolbox/process/functions/process_eigenmodes.m`, immediately **after the `Compute` subfunction** (it ends at line 261 with its closing `end`):

```matlab

%% ===== COMPUTE INTERACTIVE (called from the cortex right-click menu) =====
function ComputeInteractive(iSubject, SurfaceFile) %#ok<DEFNU>
    % Ask number of eigenmodes
    [res, isCancel] = java_dialog('input', 'Number of eigenmodes:', 'Compute eigenmodes', [], '300');
    if isCancel || isempty(res)
        return;
    end
    nModes = str2double(res);
    if isnan(nModes) || (nModes < 1)
        bst_error('Number of eigenmodes must be a positive integer.', 'Compute eigenmodes', 0);
        return;
    end
    nModes = round(nModes);
    % Ask mass matrix type (combo returns the selected string, [] if cancelled)
    MassType = java_dialog('combo', 'Mass matrix type:', 'Compute eigenmodes', [], ...
                           {'barycentric', 'voronoi', 'galerkin'});
    if isempty(MassType)
        return;
    end
    % Confirm overwrite if eigenmodes already exist on this surface
    [~, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if isComputed
        isOverwrite = java_dialog('confirm', ...
            'Eigenmodes already exist for this surface. Overwrite?', 'Compute eigenmodes');
        if ~isOverwrite
            return;
        end
    end
    % Compute: RemoveDC=true, Repair=false, isInteractive=true (manifold prompt handled in Compute)
    errMsg = Compute(SurfaceFile, nModes, MassType, true, false, true);
    if ~isempty(errMsg)
        bst_error(errMsg, 'Compute eigenmodes', 0);
        return;
    end
    % Reload the subject so the database reflects the new eigenmodes
    db_reload_subjects(iSubject);
    % Confirmation
    java_dialog('msgbox', ...
        sprintf('Computed %d eigenmodes (%s mass) on:\n%s', nModes, MassType, SurfaceFile), ...
        'Compute eigenmodes');
end
```

- [ ] **Step 2: Static-analysis check (headless)**

Run MATLAB MCP `check_matlab_code` on `toolbox/process/functions/process_eigenmodes.m`.
Expected: no errors/undefined symbols. (`Compute`, `in_tess_eigenmodes`, `java_dialog`, `bst_error`, `db_reload_subjects` are all in scope/on path.)

- [ ] **Step 3: Confirm existing process tests still pass (no regression)**

Run via MATLAB MCP `evaluate_matlab_code`: `test_process_eigenmodes_options`
Expected: PASS — adding a subfunction must not change the process description/options.

- [ ] **Step 4: Commit**

```bash
git add toolbox/process/functions/process_eigenmodes.m
git commit -m "process_eigenmodes: add ComputeInteractive (dialog wrapper over Compute)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Wire the two menu items into the cortex node

**Files:**
- Modify: `toolbox/tree/tree_callbacks.m` (surface-node block, immediately after the `Display` item at ~line 1172)

GUI code — no unit test. Verified by M-Lint here and by interactive validation in Task 5.

- [ ] **Step 1: Insert the menu items**

In `tree_callbacks.m`, in the block `case {'scalp', 'cortex', 'outerskull', 'innerskull', 'other'}`, find the `% === DISPLAY ===` group that ends at the `Display` menu item (~line 1172):

```matlab
                else
                    gui_component('MenuItem', jPopup, [], 'Display', IconLoader.ICON_DISPLAY, [], @(h,ev)view_surface(filenameRelative));
                end
```

Immediately **after** that closing `end` (before `% === SET SURFACE TYPE ===`), insert:

```matlab

                % === EIGENMODES (cortex only) ===
                if strcmpi(nodeType, 'cortex')
                    AddSeparator(jPopup);
                    if ~bst_get('ReadOnly')
                        gui_component('MenuItem', jPopup, [], 'Compute eigenmodes', IconLoader.ICON_SURFACE_CORTEX, [], @(h,ev)bst_call(@process_eigenmodes, 'ComputeInteractive', iSubject, filenameRelative));
                    end
                    gui_component('MenuItem', jPopup, [], 'View eigenmodes', IconLoader.ICON_RESULTS, [], @(h,ev)bst_call(@view_eigenmodes, filenameRelative));
                end
```

`iSubject` (line 1164) and `filenameRelative` (used at line 1171) are already in scope here. `AddSeparator`, `gui_component`, `IconLoader`, `bst_call` are all available in `tree_callbacks.m`.

- [ ] **Step 2: Static-analysis check (headless)**

Run MATLAB MCP `check_matlab_code` on `toolbox/tree/tree_callbacks.m`.
Expected: no NEW errors/undefined symbols introduced by the insertion (the file is large; confirm the added lines are clean and braces/`end`s balance).

- [ ] **Step 3: Commit**

```bash
git add toolbox/tree/tree_callbacks.m
git commit -m "tree_callbacks: add Compute/View eigenmodes to the cortex right-click menu

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Interactive validation (user) + report

GUI behavior cannot be unit-tested; this is the user-driven validation checkpoint described in the spec. Provide the exact click-path, surface findings, and report.

- [ ] **Step 1: Launch Brainstorm from this checkout**

Start Brainstorm (e.g. via `/brainstorm-start`) on a subject whose cortex is a manifold ico surface with a matching head model (per the eigenmode project notes, an icosphere/ico5 cortex avoids the OMEGA drift). If no suitable substrate exists, build one from the OMEGA sub-0002 harness (`dev/tests/test_omega_icosphere_sourcemap.m`).

- [ ] **Step 2: Validate Compute**

Right-click the cortex surface → **Compute eigenmodes** → enter a mode count (e.g. 100 for a quick run) → pick a mass type → confirm. Expect a progress bar, then the confirmation dialog. Re-open the menu and run again to confirm the **overwrite** prompt appears.

- [ ] **Step 3: Validate View + arrow navigation**

Right-click the cortex → **View eigenmodes**. Expect a surface figure titled `Eigenmodes: …` with mode 1 shown and the bottom-left legend reporting mode index + eigenvalue + wavelength. Press **→/←** to step modes (colormap updates live; higher modes show finer spatial structure), **PgUp/PgDn** to jump ±10, **H** for the help popup, and confirm rotate/zoom still work (unhandled keys fall through).

- [ ] **Step 4: Validate the error path**

On a cortex **without** eigenmodes, choose **View eigenmodes** → expect the friendly "Run Compute eigenmodes first" error.

- [ ] **Step 5: Check the bipolar colormap appearance (the one visual unknown)**

Confirm modes render as **signed** +/- lobes (diverging colors), not single-polarity magnitude. If they appear absolute/single-polarity, the `source` colormap is in absolute mode — fix by either toggling "Absolute values" off in the figure's colormap right-click menu, or, to make it the default, adjust the viewer to register/use a non-absolute diverging colormap. Capture which fix was applied.

- [ ] **Step 6: Report to the user**

Report: which substrate was used, mode count / mass type, that Compute + overwrite + View + arrows + PgUp/PgDn + H + error-path all behaved, and the colormap outcome (Step 5). Stop for the user to validate by interacting themselves, per the agreed loop.

---

## Self-Review

**Spec coverage:**
- Menu hooks (cortex: Compute + View, after Display, View works read-only) → Task 4. ✓
- `ComputeInteractive` (# modes + mass type dialog, RemoveDC on / Repair off, overwrite confirm, reuse `Compute`) → Task 3. ✓
- `view_eigenmodes.m` transient viewer (load+error, surface overlay, ←/→ + PgUp/PgDn + H, symmetric bipolar CLim, eigenvalue/wavelength legend, no DB node) → Tasks 1–2, validated Task 5. ✓
- Pure helpers `GetModeDisplay` / `StepMode` headlessly tested → Task 1. ✓
- Error handling (no modes, non-manifold via Compute, dialog cancel, clamping, all-zero guard) → Tasks 1 (clamp/zero), 2 (no-modes), 3 (cancel; non-manifold delegated to existing Compute). ✓
- Spec's compute-core test → intentionally dropped (reconciliation section); compute already covered by existing tests. ✓

**Placeholder scan:** No TBD/TODO; every code step contains complete code; every run step states the exact MCP call and expected result. The Step-5 colormap item is a real interactive check with two concrete remedies, not a placeholder.

**Type/name consistency:** `view_eigenmodes('GetModeDisplay'|'StepMode', …)` dispatch matches the subfunction names and the test calls; `GetModeDisplay` returns `.Data/.CLim/.iMode/.nModes/.Lambda/.Wavelength/.Label` (used identically in the test and in `UpdateMode`); `StepMode(iMode, delta, nModes)` signature matches all call sites; `Compute(SurfaceFile, nModes, MassType, RemoveDC, Repair, isInteractive)` matches the existing subfunction (process_eigenmodes.m:186) and the `ComputeInteractive` call.
