# M3 Plan C — connection-phase cortex viewer (core)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A cortex viewer `view_connection_phase` that shows a connection-Laplacian eigenmode: the FS-gauge phase as a cyclic surface colormap (the between-subject location coordinate), the gauge-independent 3D field as glyphs (the reference-free within-subject winding), and the field singularities as markers — with keyboard controls to browse modes and toggle layers.

**Architecture:** Reuse the proven Brainstorm patterns mapped from `view_tangents` (open opaque surface via `view_surface` + `panel_surface` opacity/edges; glyph overlays via headless `quiver3` tagged for redraw; lollipop markers; `KeyPressFcn` with fallthrough) and `view_eigenmodes` (register a transient Source result with `ImageGridAmp` + `view_surface_data` + `panel_surface('UpdateSurfaceData')` for the scalar colormap). The data pipeline chains M2/M3: `bst_conn_eigenmodes_ensure` → nxr `vertexFrame` → `bst_tangent_face2vertex` (FS frame) → `bst_conn_phase`.

**Tech Stack:** MATLAB, Brainstorm GUI (`view_surface`, `view_surface_data`, `panel_surface`, `bst_colormaps`, `bst_figures`, `bst_progress`), nxr `+nxr` wrappers, the Plan A/B functions.

**Reference spec:** `dev/connection_phase_readout_integration.md` (§3.4, §4).

**Depends on:** Plan A (`nxr.manifold.measure.vertexFrame`, installed), Plan B (`bst_tangent_face2vertex`, `bst_conn_phase`), M2 (`bst_conn_eigenmodes_ensure`).

**Scope (this plan = the core viewer).** Deferred to a follow-on "viewer layers" plan: the scalar-LBO-eigenmode colormap driver (reusing `panel_eigenmodes`), the extrinsic FS-frame and intrinsic-frame glyph layers, iso-phase contour/stripe renderings, streamlines, a Java control panel, and the `view_tangents`/`view_eigenmodes` layer-module refactor. The within-subject *colormap* (intrinsic phase) is represented here by the field glyphs (reference-free); a cyclic scalar colormap uses the FS gauge (between-subject), per spec §2.

**GUI verification note.** GUI layers are smoke-tested (figure handle valid; tagged layer objects exist with expected counts; key callbacks change state) — they cannot fully self-verify visually. After each task, in addition to the smoke test, the viewer is launched in the live session and a screenshot confirmed by the controller.

---

## File Structure

| Path | Responsibility |
|---|---|
| `toolbox/gui/view_connection_phase.m` | The viewer: data pipeline, surface + phase colormap, glyph/marker overlays, keyboard controls. |
| `dev/tests/test_view_connection_phase.m` | Smoke test: figure opens with surface data + glyph + marker layers; keyboard mode-step works. |

Single viewer file (the glyph-geometry helper `ArrowField` is copied locally from `view_tangents`, ~10 lines — local duplication is the established Brainstorm idiom and avoids touching the working `view_tangents`).

---

## Task 1: Viewer shell + data pipeline + phase colormap

**Files:** Create `toolbox/gui/view_connection_phase.m`; Test `dev/tests/test_view_connection_phase.m`.

- [ ] **Step 1: Write the failing smoke test**

Create `dev/tests/test_view_connection_phase.m`:

```matlab
function test_view_connection_phase
% Smoke test: view_connection_phase opens a cortex figure carrying the connection
% phase as surface data plus a field-glyph layer, singularity markers, and a
% working keyboard mode-step. Works on a temp COPY (the viewer caches ConnEigenmodes).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

srcFile = find_cortex_20484V();
if isempty(srcFile)
    fprintf('SKIP: no 20484-vertex cortex in the current protocol.\n');
    return;
end
tmpFile = fullfile(tempdir, 'tess_cortex_connphaseview.mat');
copyfile(file_fullpath(srcFile), tmpFile);
cleanup = onCleanup(@() cleanupFig(tmpFile));

% Small mode count for a fast open.
hFig = view_connection_phase(tmpFile, 'nModes', 12, 'MaxArrows', 400);
assert(ishandle(hFig), 'view_connection_phase did not return a valid figure handle.');

hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
assert(~isempty(hAxes), 'No Axes3D in the figure.');

% Surface data is registered (the phase colormap drives the patch).
TessInfo = getappdata(hFig, 'Surface');
assert(~isempty(TessInfo) && ~isempty(TessInfo(1).DataSource.FileName), ...
    'Figure should carry surface data (the phase result).');

% Field-glyph layer + singularity markers exist (Task 2/Task 3 add these; for the
% shell-only first pass this asserts the appdata state object is present).
st = getappdata(hFig, 'ConnPhase');
assert(~isempty(st) && isfield(st, 'Rank') && st.Rank == 1, 'ConnPhase state must start at Rank 1.');
assert(isfield(st, 'nModes') && st.nModes >= 2, 'ConnPhase state must record the mode count.');

fprintf('PASSED: view_connection_phase opens with phase surface data (Rank %d of %d modes).\n', st.Rank, st.nModes);
fprintf('ALL TESTS PASSED: test_view_connection_phase\n');
end


function cleanupFig(tmpFile)
    hs = findobj(0, 'Type', 'figure');
    for h = hs(:)'
        if ~isempty(getappdata(h, 'ConnPhase')), close(h); end
    end
    if exist(tmpFile, 'file'), delete(tmpFile); end
end


function SurfaceFile = find_cortex_20484V()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects)
    return;
end
allSubj = [sSubjects.Subject];
fallback = '';
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex')
            continue;
        end
        try
            T = load(file_fullpath(surf(iF).FileName), 'Vertices', 'Reg');
        catch
            continue;
        end
        if size(T.Vertices, 1) ~= 20484
            continue;
        end
        hasReg = isfield(T,'Reg') && isstruct(T.Reg) && isfield(T.Reg,'Sphere') ...
                 && isfield(T.Reg.Sphere,'Vertices') && ~isempty(T.Reg.Sphere.Vertices);
        if hasReg
            SurfaceFile = surf(iF).FileName;
            return;
        elseif isempty(fallback)
            fallback = surf(iF).FileName;
        end
    end
end
if isempty(SurfaceFile)
    SurfaceFile = fallback;
end
end
```

- [ ] **Step 2: Run the test to verify it FAILS**

Run `dev/tests/test_view_connection_phase.m` via the MATLAB MCP.
Expected: FAIL — `Unrecognized function or variable 'view_connection_phase'`.

- [ ] **Step 3: Write the viewer (shell + pipeline + phase colormap)**

Create `toolbox/gui/view_connection_phase.m`:

```matlab
function hFig = view_connection_phase(SurfaceFile, varargin)
% VIEW_CONNECTION_PHASE: Connection-Laplacian eigenmode viewer (phase + field).
%
% USAGE:  hFig = view_connection_phase(SurfaceFile)
%         hFig = view_connection_phase(SurfaceFile, 'nModes', 50, 'MaxArrows', 2000)
%
% DESCRIPTION:
%     Opens an opaque cortex figure showing a connection-Laplacian eigenmode:
%       - SURFACE COLORMAP: the FreeSurfer-gauge phase (cyclic) -- the between-
%         subject location coordinate (arg of the field in the per-vertex
%         trivial-connection frame).
%       - FIELD GLYPHS: the gauge-independent 3D tangent field (reference-free,
%         shows the intrinsic winding).  [Task 2]
%       - SINGULARITY MARKERS: where the field vanishes.                [Task 2]
%     Keyboard: m / M previous / next mode rank; g toggle glyphs; p toggle
%     singularities; Left / Right fewer / more glyphs.                  [Task 3]
%
%     Data pipeline: bst_conn_eigenmodes_ensure -> nxr vertexFrame ->
%     bst_tangent_face2vertex (FS frame) -> bst_conn_phase.
%
% Requires the nxr-compute plugin and a FreeSurfer-registered cortex.
%
% SEE ALSO: bst_conn_phase, bst_conn_eigenmodes_ensure, tess_tangents, view_eigenmodes

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
nModes    = 50;
MaxArrows = 2000;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'nmodes',    nModes = varargin{i+1};
        case 'maxarrows', MaxArrows = varargin{i+1};
    end
end

bst_progress('start', 'Connection phase viewer', 'Preparing eigenmodes and frames...');

%% ===== DATA PIPELINE =====
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
if ~isOk
    bst_progress('stop');
    error('view_connection_phase:nxr', 'nxr-compute required: %s', errMsg);
end
bst_plugin('Load', 'nxr-compute');

ConnEig = bst_conn_eigenmodes_ensure(SurfaceFile, nModes);
TessMat = in_tess_bst(SurfaceFile);
Vtx = TessMat.Vertices;
Fcs = double(TessMat.Faces);
Nv  = TessMat.VertNormals;

mctx   = nxr.manifold.context(Vtx, Fcs);
vFrame = nxr.manifold.measure.vertexFrame(mctx);
[Uf, ~]  = tess_tangents(SurfaceFile, 'NoSave', 1);
[Uv, Vv] = bst_tangent_face2vertex(Fcs, Uf, Nv);
FsFrame  = struct('e1', Uv, 'e2', Vv);

% Phase/field/singularities for the Fiedler (Rank 1), FS gauge.
R = bst_conn_phase(ConnEig, vFrame, 'Rank', 1, 'FsFrame', FsFrame, 'nSing', 2);

%% ===== SURFACE FIGURE + PHASE COLORMAP =====
% Register a transient Source result whose ImageGridAmp is the phase (mirrors
% view_eigenmodes). Off-support / undefined-phase vertices are set to 0.
phase = R.Phase;
phase(~isfinite(phase)) = 0;
[sStudy, iStudy] = GetAnyStudy(SurfaceFile);
ResMat = db_template('resultsmat');
ResMat.ImageGridAmp  = [phase, phase];     % 2 samples => valid static result
ResMat.ImagingKernel = [];
ResMat.nComponents   = 1;
ResMat.Time          = [0 1];
ResMat.SurfaceFile   = SurfaceFile;
ResMat.HeadModelType = 'surface';
ResMat.nAvg          = 1;
ResMat.Leff          = 1;
ResMat.ColormapType  = 'connphase';        % dedicated cyclic colormap type
ResMat.Comment       = sprintf('Connection phase (Fiedler, FS gauge)');
ResMat = bst_history('add', ResMat, 'connphase_view', 'Transient connection-phase viewer result');
StudyDir   = bst_fileparts(file_fullpath(sStudy.FileName));
OutputFile = bst_process('GetNewFilename', StudyDir, 'results_connphase');
bst_save(OutputFile, ResMat, 'v6');
db_add_data(iStudy, OutputFile, ResMat);

% Dedicated cyclic colormap (full-range signed phase in [-pi, pi]).
EnsureCyclicColormap();

hFig = view_surface_data(SurfaceFile, file_short(OutputFile));
hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
panel_surface('SetSurfaceSmooth', hFig, 1, 0, 0);   % opaque, unsmoothed
panel_surface('SetSurfaceEdges',  hFig, 1, 1);      % wireframe on

%% ===== STORE STATE (for Task 2/3 layers + redraw) =====
st = struct();
st.SurfaceFile = SurfaceFile;
st.ResultsFile = file_short(OutputFile);
st.ConnEig     = ConnEig;
st.vFrame      = vFrame;
st.FsFrame     = FsFrame;
st.R           = R;
st.Rank        = 1;
st.nModes      = max(ConnEig.CompRank(:));
st.MaxArrows   = MaxArrows;
st.Vtx         = Vtx;
st.Fcs         = Fcs;
st.showGlyphs  = true;
st.showSing    = true;
setappdata(hFig, 'ConnPhase', st);

bst_progress('stop');
end


%% ===== ENSURE A CYCLIC COLORMAP TYPE EXISTS =====
function EnsureCyclicColormap()
    % Register/refresh a dedicated 'connphase' colormap: a cyclic HSV map over the
    % full signed range [-pi, pi]. Using a DEDICATED type (not 'source') so this
    % viewer never alters other source displays.
    N = 256;
    CMap = hsv2rgb([linspace(0, 1, N)', ones(N, 1), ones(N, 1)]);
    sColormap = bst_colormaps('GetColormap', 'connphase');
    if isempty(sColormap)
        % Seed from 'source' then override; if the type is unknown, GetColormap
        % returns [] -> build a minimal struct compatible with bst_colormaps.
        sColormap = bst_colormaps('GetColormap', 'source');
    end
    sColormap.Name             = 'Connection phase (cyclic)';
    sColormap.CMap             = CMap;
    sColormap.isAbsoluteValues = 0;       % signed: phase wraps across 0
    sColormap.DisplayColorbar  = 1;
    sColormap.MaxMode          = 'custom';
    sColormap.MinValue         = -pi;
    sColormap.MaxValue         = pi;
    bst_colormaps('SetColormap', 'connphase', sColormap);
end


%% ===== FIND ANY STUDY TO HOST THE TRANSIENT RESULT =====
function [sStudy, iStudy] = GetAnyStudy(SurfaceFile)
    % Host the transient result in the surface's subject default study (mirrors
    % how view_eigenmodes attaches a transient result to the DB).
    sSubject = bst_get('SurfaceFile', SurfaceFile);
    [sStudy, iStudy] = bst_get('DefaultStudy', sSubject.iSubject);
    if isempty(iStudy)
        [sStudy, iStudy] = bst_get('Study');   % current study fallback
    end
end
```

- [ ] **Step 4: Run the test to verify it PASSES**

Run `dev/tests/test_view_connection_phase.m` via the MATLAB MCP.
Expected: `ALL TESTS PASSED: test_view_connection_phase`. Paste the output.

The two API points most likely to need a small live fix (adjust and report what you changed, do not change the math):
- **Cyclic colormap registration** (`EnsureCyclicColormap`): if `bst_colormaps('GetColormap','connphase')`/`SetColormap` on an unknown type errors, use `bst_colormaps('CreateColormap','connphase')` or `NewCustomColormap` per the actual `bst_colormaps` API, then set `.CMap`. The goal: a dedicated cyclic type so other source figures are untouched.
- **Host study** (`GetAnyStudy`): if `DefaultStudy` is empty for this subject, fall back to any valid study so `db_add_data` succeeds (mirror `view_eigenmodes`'s study lookup if needed).

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/gui/view_connection_phase.m dev/tests/test_view_connection_phase.m
git commit -m "feat(conn-phase): view_connection_phase shell + phase colormap"
```

---

## Task 2: Field-glyph + singularity-marker overlays

**Files:** Modify `toolbox/gui/view_connection_phase.m`; extend `dev/tests/test_view_connection_phase.m`.

- [ ] **Step 1: Extend the test (add layer assertions)**

In `dev/tests/test_view_connection_phase.m`, immediately before the final two `fprintf` lines, insert:

```matlab
% Field-glyph layer present, count capped by MaxArrows.
hG = findobj(hAxes, 'Tag', 'connPhaseField');
assert(~isempty(hG), 'Missing connPhaseField glyph layer.');
nG = numel(hG.UData);
assert(nG <= 400, 'Glyph count (%d) must respect MaxArrows (400).', nG);
assert(strcmpi(hG.ShowArrowHead, 'off'), 'Field glyphs should be headless.');
% Singularity markers present (2 per component).
hS = findobj(hAxes, 'Tag', 'connPhaseSing');
assert(~isempty(hS), 'Missing connPhaseSing markers.');
assert(numel(hS.XData) == numel(st.R.Singularities), 'Singularity marker count mismatch.');
```

- [ ] **Step 2: Run it to verify it FAILS**

Run the test. Expected: FAIL at `Missing connPhaseField glyph layer` (the shell from Task 1 draws no overlays yet).

- [ ] **Step 3: Add the draw helpers + call them**

In `view_connection_phase.m`, immediately before `setappdata(hFig, 'ConnPhase', st);`, add a call to draw the overlays:

```matlab
DrawOverlays(hFig, hAxes, st);
```

Wait — `st` must exist before drawing. Reorder so `st` is built first, then `setappdata`, then `DrawOverlays(hFig, hAxes, st)`. Concretely, replace the `%% ===== STORE STATE ... setappdata(hFig, 'ConnPhase', st);` block's end with:

```matlab
setappdata(hFig, 'ConnPhase', st);
DrawOverlays(hFig, hAxes, st);
```

Then add these subfunctions (after `EnsureCyclicColormap`):

```matlab
%% ===== DRAW OVERLAYS (field glyphs + singularity markers) =====
function DrawOverlays(hFig, hAxes, st)
    delete(findobj(hAxes, '-depth', 1, '-regexp', 'Tag', '^connPhase(Field|Sing)'));
    Vtx = st.Vtx;  Fcs = st.Fcs;  R = st.R;

    % --- Field glyphs (per-vertex, subsampled to MaxArrows on the support) ---
    if st.showGlyphs
        supp = find(any(R.Field ~= 0, 2));
        if ~isempty(supp)
            step = max(1, ceil(numel(supp) / st.MaxArrows));
            idx  = supp(1:step:end);
            % Per-vertex glyph length from local edge scale.
            meanEdge = MeanEdgeLength(Vtx, Fcs);
            len = 0.8 * meanEdge;
            n   = st.vFrame.normals(idx, :);
            unitN = n ./ max(sqrt(sum(n.^2,2)), eps);
            w   = R.Field(idx, :);
            wlen = sqrt(sum(w.^2, 2));
            dirw = w ./ max(wlen, eps);
            B    = Vtx(idx, :) + (0.1 * meanEdge) .* unitN;   % lift off surface
            Vec  = dirw .* len;
            quiver3(B(:,1), B(:,2), B(:,3), Vec(:,1), Vec(:,2), Vec(:,3), 0, ...
                'Parent', hAxes, 'Color', [0.95 0.85 0.1], 'LineWidth', 1.0, ...
                'ShowArrowHead', 'off', 'Tag', 'connPhaseField');
        end
    end

    % --- Singularity markers ---
    if st.showSing && ~isempty(R.Singularities)
        P = Vtx(R.Singularities, :);
        plot3(P(:,1), P(:,2), P(:,3), 'o', 'Parent', hAxes, ...
            'MarkerFaceColor', [0.9 0.1 0.1], 'MarkerEdgeColor', [.2 .2 .2], ...
            'MarkerSize', 10, 'LineStyle', 'none', 'Tag', 'connPhaseSing');
    end
end


%% ===== MEAN EDGE LENGTH =====
function L = MeanEdgeLength(Vtx, Fcs)
    e1 = Vtx(Fcs(:,2),:) - Vtx(Fcs(:,1),:);
    e2 = Vtx(Fcs(:,3),:) - Vtx(Fcs(:,2),:);
    e3 = Vtx(Fcs(:,1),:) - Vtx(Fcs(:,3),:);
    L = mean([sqrt(sum(e1.^2,2)); sqrt(sum(e2.^2,2)); sqrt(sum(e3.^2,2))]);
end
```

- [ ] **Step 4: Run the test to verify it PASSES**

Run the test. Expected: `ALL TESTS PASSED`. Paste the output.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/gui/view_connection_phase.m dev/tests/test_view_connection_phase.m
git commit -m "feat(conn-phase): field-glyph + singularity-marker overlays"
```

---

## Task 3: Keyboard controls (mode browse + layer toggles)

**Files:** Modify `toolbox/gui/view_connection_phase.m`; extend `dev/tests/test_view_connection_phase.m`.

- [ ] **Step 1: Extend the test (mode-step + toggle)**

In `dev/tests/test_view_connection_phase.m`, immediately before the final two `fprintf` lines, insert:

```matlab
% Keyboard: 'M' advances the mode rank and recomputes the phase/field.
kp = get(hFig, 'KeyPressFcn');
assert(~isempty(kp), 'Viewer must install a KeyPressFcn.');
kp(hFig, struct('Key','m','Modifier',{{'shift'}}));   % next mode (Shift+m == 'M')
st2 = getappdata(hFig, 'ConnPhase');
assert(st2.Rank == 2, 'Shift+m should advance to Rank 2 (got %d).', st2.Rank);
% Toggle glyphs off.
kp(hFig, struct('Key','g','Modifier',{{}}));
assert(isempty(findobj(hAxes, 'Tag', 'connPhaseField')), 'g should toggle glyphs off.');
```

- [ ] **Step 2: Run it to verify it FAILS**

Run the test. Expected: FAIL — either no `KeyPressFcn` installed, or `Rank` does not advance.

- [ ] **Step 3: Install the key handler + recompute-on-mode-change**

In `view_connection_phase.m`, immediately after `DrawOverlays(hFig, hAxes, st);` (end of the main function), add:

```matlab
set(hFig, 'KeyPressFcn', @(h,ev) KeyPress_Callback(h, ev));
```

Then add these subfunctions:

```matlab
%% ===== KEYBOARD =====
function KeyPress_Callback(hFig, ev)
    st = getappdata(hFig, 'ConnPhase');
    if isempty(st), return; end
    switch ev.Key
        case 'm'
            if ismember('shift', ev.Modifier)
                st.Rank = min(st.nModes, st.Rank + 1);
            else
                st.Rank = max(1, st.Rank - 1);
            end
            st = Recompute(hFig, st);
        case 'g'
            st.showGlyphs = ~st.showGlyphs;
        case 'p'
            st.showSing = ~st.showSing;
        case 'rightarrow'
            st.MaxArrows = ceil(st.MaxArrows * 1.5);
        case 'leftarrow'
            st.MaxArrows = max(50, floor(st.MaxArrows / 1.5));
        otherwise
            return;
    end
    setappdata(hFig, 'ConnPhase', st);
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    DrawOverlays(hFig, hAxes, st);
end


%% ===== RECOMPUTE PHASE/FIELD FOR A NEW MODE RANK =====
function st = Recompute(hFig, st)
    global GlobalData;
    st.R = bst_conn_phase(st.ConnEig, st.vFrame, 'Rank', st.Rank, 'FsFrame', st.FsFrame, 'nSing', 2);
    % Update the surface phase scalar in place (mirrors view_eigenmodes).
    phase = st.R.Phase;
    phase(~isfinite(phase)) = 0;
    [iDS, iResult] = bst_memory('GetDataSetResult', st.ResultsFile);
    GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp = [phase, phase];
    panel_surface('UpdateSurfaceData', hFig);
end
```

- [ ] **Step 4: Run the test to verify it PASSES**

Run the test. Expected: `ALL TESTS PASSED`. Paste the output.

If `bst_memory('GetDataSetResult', ...)` does not resolve the transient result handle, use the `[iDS,iResult]` returned at registration time (store them in `st` during Task 1) — adjust and report.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/gui/view_connection_phase.m dev/tests/test_view_connection_phase.m
git commit -m "feat(conn-phase): keyboard mode-browse + layer toggles"
```

---

## Notes for the implementer

- **Run via the MATLAB MCP**; never bare `clear`. Resolve the cortex by vertex count; SKIP if absent. Work on the temp copy (the viewer caches `ConnEigenmodes`).
- **GUI tasks are smoke-tested** (handle valid, tagged layer objects exist with expected counts, key callbacks change state). The controller additionally launches the viewer live and confirms a screenshot after each task.
- **Two known live-iteration points** (Task 1 Step 4): the `bst_colormaps` cyclic-type registration and the host-study lookup. Use the real API; report any adjustment. Do not alter the phase math.
- **Deferred to a follow-on plan:** scalar-eigenmode colormap driver, extrinsic/intrinsic frame glyph layers, iso-phase contours, stripes, streamlines, Java control panel, and the `view_tangents`/`view_eigenmodes` layer-module refactor.
```
