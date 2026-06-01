# view_tangents Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `view_tangents.m`, an interactive 3D viewer that overlays the per-face tangent frame field stored by `tess_tangents` (`TessMat.TangentFrame`) on a Brainstorm surface figure, and wire it into the cortex node's right-click menu.

**Architecture:** Standalone GUI function modeled on `toolbox/gui/view_leadfield_vectors.m`. It opens a `view_surface` figure, `hold on`, and draws the frame with three `quiver3` sets (U and V in one color, the face normal in another) plus `plot3` markers at the singularity poles. A hijacked `KeyPressFcn` gives interactive density/length/width controls and toggles, falling through to the saved callback for unhandled keys. Compute stays in `tess_tangents`; this file is display-only and auto-computes the frame if it is missing.

**Tech Stack:** MATLAB, Brainstorm toolbox (`view_surface`, `in_tess_bst`, `file_fullpath`, `tess_normals`, `tess_tangents`, `gui_component`, `bst_call`), MATLAB built-in `quiver3`/`plot3`. Tests run via the MATLAB MCP `evaluate_matlab_code` tool (function-style tests invoked directly, NOT `runtests`).

---

## Reference details (verified in the codebase)

- `view_surface(SurfaceFile, SurfAlpha, SurfColor, TargetFigure)` — pass `'NewFigure'` as the 4th arg to create a new figure; returns `hFig`. The 3D axes carry `Tag = 'Axes3D'`.
- `tess_normals(Vertices, Faces)` returns `[VertNormals, FaceNormals]`, both unit-normalized.
- `tess_tangents(SurfaceFile)` computes, stores `TessMat.TangentFrame` (fields: `Domain='face'`, `U`/`V` as `single` `[nF x 3]`, `Singularities` with `.Vertices`/`.Indices`/`.Hemisphere`), and returns `[U, V]`.
- `dev/tests/test_tess_tangents.m` contains `find_registered_cortex()` — finds the smallest cortex in the DB that has `Reg.Sphere.Vertices`. Reuse it (repeated verbatim in the test task below).
- Cortex node popup is built in `toolbox/tree/tree_callbacks.m`, case `{'scalp','cortex',...}` near line 1162; the initial "Display" menu item is at lines 1168–1172; `filenameRelative` and `IconLoader` are in scope there.

## File Structure

- **Create** `toolbox/gui/view_tangents.m` — the viewer. Public entry `hFig = view_tangents(SurfaceFile, varargin)`. Internals: nested `DrawArrows()` and `KeyPress_Callback()`, subfunctions `ArrowSubsample(nF, nArrows)` and `ArrowField(...)` (pure geometry).
- **Modify** `toolbox/tree/tree_callbacks.m` — one `gui_component('MenuItem', ...)` in the cortex popup.
- **Create** `dev/tests/test_view_tangents.m` — integration test (auto-compute, render, geometry-via-readback, interaction).

---

## Task 1: view_tangents core — load, auto-compute, render frame + singularities

**Files:**
- Create: `toolbox/gui/view_tangents.m`
- Create: `dev/tests/test_view_tangents.m`

- [ ] **Step 1: Write the failing test** (`dev/tests/test_view_tangents.m`)

```matlab
function test_view_tangents
% Integration test for view_tangents: auto-compute path, render, frame
% geometry (via the drawn quiver objects), and interactive callbacks.
% Runs on a temp COPY so the DB surface is not mutated.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% --- Find a low-res FreeSurfer-registered cortex in the DB ---
srcFile = find_registered_cortex();
assert(~isempty(srcFile), 'No FreeSurfer-registered cortex found in the DB to test against.');
fprintf('Source cortex: %s\n', srcFile);

% --- Work on a temp copy with NO stored frame (exercise auto-compute) ---
tmpFile = fullfile(tempdir, 'tess_cortex_viewtan_test.mat');
copyfile(file_fullpath(srcFile), tmpFile);
cleanup = onCleanup(@() cleanupFcn(tmpFile));
T = load(tmpFile);
if isfield(T, 'TangentFrame')
    T = rmfield(T, 'TangentFrame');
    save(tmpFile, '-struct', 'T');
end

% --- Render ---
hFig = view_tangents(tmpFile, 'MaxArrows', 500);
assert(ishandle(hFig), 'view_tangents did not return a valid figure handle.');

% --- Auto-compute fired: frame now stored ---
Tafter = load(tmpFile, 'TangentFrame');
assert(isfield(Tafter, 'TangentFrame') && ~isempty(Tafter.TangentFrame), ...
    'Auto-compute did not store TangentFrame on the file.');

% --- Three quiver sets with equal counts; U/V share a color != normal ---
hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
hU = findobj(hAxes, 'Tag', 'tangentU');
hV = findobj(hAxes, 'Tag', 'tangentV');
hN = findobj(hAxes, 'Tag', 'tangentN');
assert(~isempty(hU) && ~isempty(hV) && ~isempty(hN), 'Missing tangentU/tangentV/tangentN quiver objects.');
nU = numel(hU.UData);  nV = numel(hV.UData);  nN = numel(hN.UData);
assert(nU == nV && nV == nN, 'U/V/N arrow counts differ.');
assert(nU <= 500, 'Drew more arrows (%d) than requested (500).', nU);
assert(isequal(hU.Color, hV.Color), 'U and V must share a color.');
assert(~isequal(hU.Color, hN.Color), 'Normal color must differ from U/V.');

% --- Frame arrows are equal length (a frame field has no magnitude) ---
lenU = sqrt(hU.UData(:).^2 + hU.VData(:).^2 + hU.WData(:).^2);
assert((max(lenU) - min(lenU)) < 1e-6, 'U arrows are not equal length.');

% --- Singularity markers match the stored pole count ---
hS = findobj(hAxes, 'Tag', 'tangentSing');
assert(~isempty(hS), 'No singularity markers drawn.');
assert(numel(hS.XData) == numel(Tafter.TangentFrame.Singularities.Vertices), ...
    'Singularity marker count does not match stored Singularities.Vertices.');

close(hFig);
fprintf('ALL TESTS PASSED: test_view_tangents\n');
end


function cleanupFcn(tmpFile)
    if ~isempty(findobj(0, 'Type', 'figure'))
        close all force;
    end
    if exist(tmpFile, 'file')
        delete(tmpFile);
    end
end


function SurfaceFile = find_registered_cortex()
% Return a low-res cortex FileName that has Reg.Sphere.Vertices, or '' if none.
SurfaceFile = '';
best = inf;
sSubjects = bst_get('ProtocolSubjects');
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex'), continue; end
        try
            T = load(file_fullpath(surf(iF).FileName), 'Reg', 'Vertices');
        catch
            continue;
        end
        if isfield(T,'Reg') && isstruct(T.Reg) && isfield(T.Reg,'Sphere') ...
           && isfield(T.Reg.Sphere,'Vertices') && ~isempty(T.Reg.Sphere.Vertices)
            n = size(T.Vertices, 1);
            if n < best
                best = n;
                SurfaceFile = surf(iF).FileName;
            end
        end
    end
end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run (MATLAB MCP `evaluate_matlab_code`):
```matlab
run('dev/tests/test_view_tangents.m')
```
Expected: FAIL — `Unrecognized function or variable 'view_tangents'` (function does not exist yet).

- [ ] **Step 3: Write the implementation** (`toolbox/gui/view_tangents.m`)

```matlab
function hFig = view_tangents(SurfaceFile, varargin)
% VIEW_TANGENTS: Display the per-face tangent frame field stored by tess_tangents.
%
% USAGE:  hFig = view_tangents(SurfaceFile)
%         hFig = view_tangents(SurfaceFile, 'MaxArrows', 3000)
%
% DESCRIPTION:
%     Overlays the per-face tangent frame field (TessMat.TangentFrame, computed
%     by tess_tangents) on a Brainstorm surface figure. U and V (the in-plane
%     tangent cross) are drawn in one color; the face normal in another. The
%     four registration-pole singularities are marked. If the surface has no
%     stored frame, it is computed and stored via tess_tangents first.
%
%     Keyboard (figure focused):
%       Left/Right          fewer / more arrows
%       Shift + Up/Down     arrow length
%       Ctrl  + Up/Down     line width
%       N                   toggle normal arrows
%       P                   toggle singularity markers
%       H                   help
%
% INPUT:
%     - SurfaceFile : Brainstorm cortex surface (relative or full path).
% OPTIONS:
%     - MaxArrows : initial number of face frames drawn (default 2000).
% OUTPUT:
%     - hFig : handle to the 3D figure.
%
% SEE ALSO: tess_tangents, view_leadfield_vectors, tess_normals

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

%% ===== PARSE INPUTS =====
MaxArrows = 2000;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'maxarrows', MaxArrows = varargin{i+1};
    end
end

%% ===== LOAD SURFACE =====
TessMat = in_tess_bst(SurfaceFile);

%% ===== ENSURE TANGENT FRAME (auto-compute if missing) =====
if ~isfield(TessMat, 'TangentFrame') || isempty(TessMat.TangentFrame)
    tess_tangents(SurfaceFile);          % compute + store (errors propagate)
    TessMat = in_tess_bst(SurfaceFile);  % reload to get Singularities
end
TF = TessMat.TangentFrame;
if ~strcmpi(TF.Domain, 'face')
    error('view_tangents:unsupportedDomain', ...
        'view_tangents supports per-face frames only (Domain=''face''); got ''%s''.', TF.Domain);
end

Vtx = TessMat.Vertices;
Fcs = double(TessMat.Faces);
nF  = size(Fcs, 1);
U   = double(TF.U);
V   = double(TF.V);

%% ===== DISPLAY GEOMETRY (plain MATLAB) =====
[~, FaceNormals] = tess_normals(Vtx, Fcs);
Centroids = (Vtx(Fcs(:,1),:) + Vtx(Fcs(:,2),:) + Vtx(Fcs(:,3),:)) / 3;
e1 = sqrt(sum((Vtx(Fcs(:,2),:) - Vtx(Fcs(:,1),:)).^2, 2));
e2 = sqrt(sum((Vtx(Fcs(:,3),:) - Vtx(Fcs(:,2),:)).^2, 2));
e3 = sqrt(sum((Vtx(Fcs(:,1),:) - Vtx(Fcs(:,3),:)).^2, 2));
meanEdge = mean([e1; e2; e3]);
SingXYZ  = Vtx(TF.Singularities.Vertices, :);

%% ===== OPEN SURFACE FIGURE =====
hFig = view_surface(SurfaceFile, 0.5, [.5 .5 .5], 'NewFigure');
if isempty(hFig)
    error('view_tangents:noFigure', 'Could not open the surface figure.');
end
set(hFig, 'Name', ['Tangent basis: ' SurfaceFile]);
hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
hold(hAxes, 'on');

%% ===== DISPLAY STATE =====
nArrows     = min(MaxArrows, nF);
quiverSize  = 1;
quiverWidth = 1;
showNormals = true;
showSing    = true;
colTangent  = [1 1 0];   % yellow : U and V
colNormal   = [1 0 1];   % magenta: face normal
colSing     = [1 0 0];   % red    : singularities

% Hijack keyboard callback (fall through to original for unhandled keys)
KeyPressFcn_bak = get(hFig, 'KeyPressFcn');
set(hFig, 'KeyPressFcn', @KeyPress_Callback);

% Legend / status label
hLabel = uicontrol('Style', 'text', 'String', '...', 'Units', 'Pixels', ...
    'Position', [6 4 1000 18], 'HorizontalAlignment', 'left', ...
    'FontUnits', 'points', 'FontSize', bst_get('FigFont'), ...
    'ForegroundColor', [1 1 1], 'BackgroundColor', [0 0 0], 'Parent', hFig);

DrawArrows();

%% =================================================================================
%% ===== DRAW =====
    function DrawArrows()
        delete(findobj(hAxes, '-depth', 1, '-regexp', 'Tag', '^tangent'));
        idx = ArrowSubsample(nF, nArrows);
        [B, Uvec, Vvec, Nvec] = ArrowField(Centroids, FaceNormals, U, V, idx, ...
            quiverSize * meanEdge, 0.25 * meanEdge);
        % U and V (same color) — autoscale OFF (scale arg 0) for equal lengths
        quiver3(B(:,1), B(:,2), B(:,3), Uvec(:,1), Uvec(:,2), Uvec(:,3), 0, ...
            'Parent', hAxes, 'Color', colTangent, 'LineWidth', quiverWidth, 'Tag', 'tangentU');
        quiver3(B(:,1), B(:,2), B(:,3), Vvec(:,1), Vvec(:,2), Vvec(:,3), 0, ...
            'Parent', hAxes, 'Color', colTangent, 'LineWidth', quiverWidth, 'Tag', 'tangentV');
        % Face normal (different color)
        if showNormals
            quiver3(B(:,1), B(:,2), B(:,3), Nvec(:,1), Nvec(:,2), Nvec(:,3), 0, ...
                'Parent', hAxes, 'Color', colNormal, 'LineWidth', quiverWidth, 'Tag', 'tangentN');
        end
        % Singularity poles
        if showSing && ~isempty(SingXYZ)
            plot3(SingXYZ(:,1), SingXYZ(:,2), SingXYZ(:,3), 'o', ...
                'Parent', hAxes, 'MarkerFaceColor', colSing, 'MarkerEdgeColor', [.2 .2 .2], ...
                'MarkerSize', 10, 'LineStyle', 'none', 'Tag', 'tangentSing');
        end
        set(hLabel, 'String', sprintf( ...
            'U,V (yellow) - normal (magenta) - %d/%d faces shown - H for help', numel(idx), nF));
    end

%% ===== KEYBOARD CALLBACK =====
    function KeyPress_Callback(h, ev)
        switch ev.Key
            case 'rightarrow'
                nArrows = min(nF, ceil(nArrows * 1.5));
            case 'leftarrow'
                nArrows = max(min(50, nF), floor(nArrows / 1.5));
            case 'uparrow'
                if ismember('shift', ev.Modifier)
                    quiverSize = quiverSize * 1.2;
                elseif ismember('control', ev.Modifier)
                    quiverWidth = quiverWidth * 1.2;
                else
                    KeyPressFcn_bak(h, ev);  return;
                end
            case 'downarrow'
                if ismember('shift', ev.Modifier)
                    quiverSize = quiverSize / 1.2;
                elseif ismember('control', ev.Modifier)
                    quiverWidth = max(0.5, quiverWidth / 1.2);
                else
                    KeyPressFcn_bak(h, ev);  return;
                end
            case 'n'
                showNormals = ~showNormals;
            case 'p'
                showSing = ~showSing;
            case 'h'
                java_dialog('msgbox', ['<HTML><TABLE>' ...
                    '<TR><TD><B>Left / Right</B></TD><TD>Fewer / more arrows</TD></TR>' ...
                    '<TR><TD><B>Shift + Up/Down</B></TD><TD>Arrow length</TD></TR>' ...
                    '<TR><TD><B>Ctrl + Up/Down</B></TD><TD>Line width</TD></TR>' ...
                    '<TR><TD><B>N</B></TD><TD>Toggle normal arrows</TD></TR>' ...
                    '<TR><TD><B>P</B></TD><TD>Toggle singularity markers</TD></TR>' ...
                    '</TABLE>'], 'Tangent basis shortcuts', [], 0);
                return;
            otherwise
                KeyPressFcn_bak(h, ev);  return;
        end
        DrawArrows();
    end
end


%% ========================================================================
function idx = ArrowSubsample(nF, nArrows)
% Deterministic uniform stride over face index (reproducible for tests).
nArrows = max(1, min(nArrows, nF));
idx = unique(round(linspace(1, nF, nArrows)));
end


%% ========================================================================
function [B, Uvec, Vvec, Nvec] = ArrowField(Centroids, FaceNormals, U, V, idx, len, offset)
% Pure geometry: equal-length unit arrows, bases offset off the surface along
% the face normal. No figure handles — testable in isolation.
unit = @(X) X ./ max(sqrt(sum(X.^2, 2)), eps);
n = unit(FaceNormals(idx, :));
B = Centroids(idx, :) + offset .* n;
Uvec = unit(U(idx, :)) .* len;
Vvec = unit(V(idx, :)) .* len;
Nvec = n .* len;
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run (MATLAB MCP `evaluate_matlab_code`):
```matlab
run('dev/tests/test_view_tangents.m')
```
Expected: PASS — prints `ALL TESTS PASSED: test_view_tangents`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/view_tangents.m dev/tests/test_view_tangents.m
git commit -m "view_tangents: per-face tangent frame viewer (render + auto-compute)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Interactive density and toggle controls

The implementation in Task 1 already contains `KeyPress_Callback`. This task adds the assertions that lock its behavior in. (If Task 1's callback is correct, Step 3 is a no-op; keep it for the engineer reading tasks out of order.)

**Files:**
- Modify: `dev/tests/test_view_tangents.m` (extend before the final `close(hFig)`)
- Modify (if needed): `toolbox/gui/view_tangents.m`

- [ ] **Step 1: Add the failing interaction assertions**

In `dev/tests/test_view_tangents.m`, insert the following block immediately before the `close(hFig);` line in `test_view_tangents`:

```matlab
% --- Interaction: right arrow increases density ---
kp = get(hFig, 'KeyPressFcn');
evMore = struct('Key', 'rightarrow', 'Modifier', {{}});
kp(hFig, evMore);
hUmore = findobj(hAxes, 'Tag', 'tangentU');
assert(numel(hUmore.UData) > nU, 'Right arrow did not increase arrow density.');

% --- Interaction: N toggles the normal arrows off ---
evN = struct('Key', 'n', 'Modifier', {{}});
kp(hFig, evN);
assert(isempty(findobj(hAxes, 'Tag', 'tangentN')), 'N did not toggle the normal arrows off.');

% --- Interaction: P toggles the singularity markers off ---
evP = struct('Key', 'p', 'Modifier', {{}});
kp(hFig, evP);
assert(isempty(findobj(hAxes, 'Tag', 'tangentSing')), 'P did not toggle the singularity markers off.');
```

- [ ] **Step 2: Run the test to verify the new assertions are exercised**

Run (MATLAB MCP `evaluate_matlab_code`):
```matlab
run('dev/tests/test_view_tangents.m')
```
Expected: PASS if Task 1's callback is correct. If it FAILs (e.g. density did not change, or a toggle did not remove the tagged objects), proceed to Step 3.

- [ ] **Step 3: Fix the callback if a new assertion failed**

Only if Step 2 failed: in `toolbox/gui/view_tangents.m`, verify `KeyPress_Callback` handles `'rightarrow'` (increments `nArrows`), `'n'` (flips `showNormals`), and `'p'` (flips `showSing`), and that each handled key ends by calling `DrawArrows()`. The exact correct body is shown in Task 1 Step 3 — match it.

- [ ] **Step 4: Re-run to confirm pass**

Run (MATLAB MCP `evaluate_matlab_code`):
```matlab
run('dev/tests/test_view_tangents.m')
```
Expected: PASS — prints `ALL TESTS PASSED: test_view_tangents`.

- [ ] **Step 5: Commit**

```bash
git add dev/tests/test_view_tangents.m toolbox/gui/view_tangents.m
git commit -m "view_tangents: lock interactive density/toggle controls with tests

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Tree menu entry — cortex node "Display tangent basis"

**Files:**
- Modify: `toolbox/tree/tree_callbacks.m` (cortex popup, after the initial "Display" item at lines 1168–1172)

- [ ] **Step 1: Add the menu item**

In `toolbox/tree/tree_callbacks.m`, find the cortex popup "Display" block (the `if strcmpi(nodeType, 'other') ... else ... view_surface(filenameRelative) ... end` at lines 1168–1172). Immediately AFTER that `end` (before the `% === SET SURFACE TYPE ===` comment), insert:

```matlab
                % === DISPLAY TANGENT BASIS (cortex only) ===
                if strcmpi(nodeType, 'cortex') && (length(bstNodes) == 1)
                    gui_component('MenuItem', jPopup, [], 'Display tangent basis', IconLoader.ICON_DISPLAY, [], @(h,ev)bst_call(@view_tangents, filenameRelative));
                end
```

- [ ] **Step 2: Verify it parses and the callback target is on the path**

Run (MATLAB MCP `evaluate_matlab_code`):
```matlab
addpath(genpath(pwd));
err = [];
try
    pcode('toolbox/tree/tree_callbacks.m', '-inplace');  % parse-check; harmless
catch e
    err = e;
end
delete('toolbox/tree/tree_callbacks.p');
assert(isempty(err), 'tree_callbacks.m failed to parse: %s', getReport(err));
assert(exist('view_tangents', 'file') == 2, 'view_tangents is not on the MATLAB path.');
fprintf('PASSED: tree_callbacks parses and view_tangents resolves.\n');
```
Expected: PASS — prints `PASSED: tree_callbacks parses and view_tangents resolves.` (If your MATLAB lacks `pcode` write permission, substitute `checkcode('toolbox/tree/tree_callbacks.m')` and assert it reports no errors.)

- [ ] **Step 3: Confirm the menu line is present**

Run:
```bash
grep -n "Display tangent basis" toolbox/tree/tree_callbacks.m
```
Expected: one match showing the `gui_component(... @(h,ev)bst_call(@view_tangents, filenameRelative))` line.

- [ ] **Step 4: Commit**

```bash
git add toolbox/tree/tree_callbacks.m
git commit -m "view_tangents: add 'Display tangent basis' to cortex node menu

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review

**1. Spec coverage** (against `dev/view_tangents_design.md`):
- Interface (`hFig = view_tangents(SurfaceFile[, 'MaxArrows', N])`) → Task 1 Step 3. ✓
- Data flow: load, auto-compute if missing, reload, Domain validation → Task 1 Step 3. ✓
- Geometry: `tess_normals` face normals, centroids, meanEdge, normal offset, equal-length arrows → Task 1 Step 3 (`ArrowField`) + asserted in Task 1 Step 1 (equal-length). ✓
- Rendering: three `quiver3` sets, U/V yellow, normal magenta, Tags `tangentU/V/N` → Task 1. ✓
- Singularity markers (`tangentSing`, count from `Singularities.Vertices`) → Task 1. ✓
- Interaction: density (Left/Right), length (Shift+Up/Down), width (Ctrl+Up/Down), N/P toggles, H help, fall-through → Task 1 Step 3, asserted in Task 2. ✓
- Tree entry "Display tangent basis", cortex-only, `bst_call` → Task 3. ✓
- Error handling: `unsupportedDomain`; auto-compute errors propagate (not caught) → Task 1 Step 3. ✓
- Testing: integration (auto-compute + render + geometry readback) and interaction → Tasks 1–2. ✓

**2. Placeholder scan:** No "TBD"/"TODO"/"handle edge cases"/"similar to Task N" remain. `find_registered_cortex` and `cleanupFcn` are written out in full in Task 1.

**3. Type/name consistency:** Tags `tangentU`/`tangentV`/`tangentN`/`tangentSing` are identical across implementation and tests. Subfunction names `ArrowSubsample`/`ArrowField` and their signatures match their call sites. State variables (`nArrows`, `quiverSize`, `quiverWidth`, `showNormals`, `showSing`) are consistent between `DrawArrows` and `KeyPress_Callback`. `MaxArrows` option name matches the test call.

**Note on the spec's "isolated unit test of the geometry helper":** MATLAB cannot call a subfunction from an external test file. The geometry (`ArrowField`) is therefore verified through the rendered quiver objects (equal arrow length in Task 1; deterministic counts via `ArrowSubsample` implied by the `<= MaxArrows` and density-change assertions). This tests the same logic without an artificial public entry point — a faithful refinement of the spec, recorded here.
