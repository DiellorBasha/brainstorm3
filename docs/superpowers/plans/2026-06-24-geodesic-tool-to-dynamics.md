# Geodesic Area Tool → Dynamics Suite Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the custom geodesic Area tool out of `panel_scout` into a dynamics-owned `bst_geodesic_tool` that emits a **transient** heat-distance disk (seed + radius + vertices = the source-axis point+extent primitive), and rewire the atom panel's Capture to use it instead of a Scout.

**Architecture:** A new I/O-light `toolbox/dynamics/bst_geodesic_tool.m` (verb-dispatched) holds the tool state and reuses the unchanged `tess_scout_area` engine (stays in `toolbox/anatomy`); it draws a transient overlay, never a scout. `figure_3d` gets two small parallel branches routing cortical click/scroll to the tool when its pick-flag is set. The Area tool is excised from `panel_scout` (the geodesic Line tool stays). `panel_bst_dynamics` Capture reads `bst_geodesic_tool('GetState')`.

**Tech Stack:** MATLAB, Brainstorm GUI; `eval(macro_method)` dispatch; `tess_scout_area` (heat-distance engine); `select3d` (clicked-vertex resolution); tested headless with Brainstorm live.

## Global Constraints

- Phase 2 of `docs/superpowers/specs/2026-06-24-atom-tensor-architecture-analysis.md`; spec is `docs/superpowers/specs/2026-06-24-geodesic-tool-to-dynamics-design.md`.
- The tool emits a **transient disk** (overlay patch tag `'GeodesicToolDisk'`), **never** a `CreateScout` — no atlas entry, no bulk-scout-delete hazard.
- `tess_scout_area` is **unchanged** and **stays in `toolbox/anatomy`** (used by both the new tool and the scout geodesic **Line** tool).
- The geodesic **Line** tool stays in `panel_scout`; only the **Area** tool is removed.
- Constants mirror today's Area tool: initial radius `R0 = 0.003` m, scroll step `STEP = 0.003` m; scroll up (`VerticalScrollCount < 0`) grows.
- Removing `panel_scout('AreaToolScroll')`/`IsAreaToolActive` and editing the `figure_3d` scroll branch must land in the **same task** (Task 2), or the deleted calls error.
- Mutual exclusion: activating the dynamics Region tool clears `isSelectingCorticalSpot` on the figures (and the scout tools already clear their own modes), so there is never more than one cortical-pick mode at once.
- New module is verb-dispatched via `eval(macro_method)`, lives in `toolbox/dynamics/` (auto on path).
- Do not start implementation on `development`; this work is on branch `feature/geodesic-tool-to-dynamics` (spec already committed there).
- Tests run headless in MATLAB with Brainstorm live (`brainstorm nogui`, TutorialAuditory). Do not `clear`/restart Brainstorm or close all figures; `rehash` and re-run. The Subject01 cortex surface backs the geometry tests.

---

### Task 1: `bst_geodesic_tool` module (transient heat-disk)

**Files:**
- Create: `toolbox/dynamics/bst_geodesic_tool.m`
- Create/Test: `dev/test_bst_geodesic_tool.m`

**Interfaces:**
- Consumes: `tess_scout_area(SurfaceFile, seed, R [, phi])` → `[vertices, phi]` (anatomy); `in_tess_bst`; `bst_figures`; `select3d`; `file_gettype`.
- Produces:
  - `bst_geodesic_tool('Seed', SurfaceFile, vi)` — compute the disk around vertex `vi`, cache state. Headless.
  - `bst_geodesic_tool('Grow', scrollCount)` — re-threshold the cached `phi` at a new radius (`R = max(STEP, R − scrollCount·STEP)`), update cached `vertices`. Headless.
  - `st = bst_geodesic_tool('GetState')` — `struct('seed','phi','radius','vertices','pos','SurfaceFile'[, 'hFig'])` or `[]`.
  - `bst_geodesic_tool('Draw', hFig)` / `('Clear', hFig)` — transient overlay patch `'GeodesicToolDisk'`.
  - `bst_geodesic_tool('Toggle', onoff)` / `('IsActive')` — enter/exit dynamics pick mode (sets `isDynamicsGeodesicPick` appdata).
  - `bst_geodesic_tool('OnClick', hFig)` / `('OnScroll', scrollCount)` — figure_3d entry points (Task 2 wires them).

- [ ] **Step 1: Write the failing test**

Create `dev/test_bst_geodesic_tool.m`:

```matlab
function test_bst_geodesic_tool()
% TEST_BST_GEODESIC_TOOL: the dynamics heat-disk tool (seed/grow/state + draw/clear).
%
% USAGE:  test_bst_geodesic_tool   % Brainstorm running, TutorialAuditory loaded
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};  pass = true;

    % cortex surface for the geometry
    sSubject = bst_get('Subject', 'Subject01');
    if isempty(sSubject) || isempty(sSubject.iCortex)
        fprintf('SKIPPED (no Subject01 cortex)\n');  fprintf('\n==== SUITE: %s ====\n', PF{pass+1});  return;
    end
    SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
    Surf = in_tess_bst(SurfaceFile, 0);
    vi = round(size(Surf.Vertices,1)/3);

    % ---------- T1: Seed -> GetState ----------
    bst_geodesic_tool('Seed', SurfaceFile, vi);
    st = bst_geodesic_tool('GetState');
    n0 = 0;  if ~isempty(st), n0 = numel(st.vertices); end
    ok1 = ~isempty(st) && (st.seed==vi) && (n0>0) && isequal(size(st.pos),[1 3]) ...
       && strcmp(st.SurfaceFile, SurfaceFile) && isequal(st.pos, Surf.Vertices(vi,:));
    fprintf('T1 seed: verts=%d seed=%d pos=%d => %s\n', n0, (st.seed==vi), isequal(st.pos,Surf.Vertices(vi,:)), PF{ok1+1});
    pass = pass && ok1;

    % ---------- T2: Grow (up = -1 grows; down = +N shrinks; radius floors) ----------
    bst_geodesic_tool('Grow', -1);  s2 = bst_geodesic_tool('GetState');  nUp = numel(s2.vertices);
    bst_geodesic_tool('Grow', +10); s3 = bst_geodesic_tool('GetState');  nDn = numel(s3.vertices);
    ok2 = (nUp >= n0) && (nDn <= nUp) && (s3.radius >= 0.003 - 1e-12) && (s2.radius > s3.radius);
    fprintf('T2 grow: n0=%d up=%d down=%d radiusFloor=%d => %s\n', n0, nUp, nDn, (s3.radius>=0.003-1e-12), PF{ok2+1});
    pass = pass && ok2;

    % ---------- T3: Draw/Clear overlay on a figure ----------
    hFig = view_surface(SurfaceFile);  drawnow;
    bst_geodesic_tool('Seed', SurfaceFile, vi);
    bst_geodesic_tool('Draw', hFig);   drawnow;
    nPatch = numel(findobj(hFig, 'Tag', 'GeodesicToolDisk'));
    bst_geodesic_tool('Clear', hFig);  drawnow;
    nAfter = numel(findobj(hFig, 'Tag', 'GeodesicToolDisk'));
    ok3 = (nPatch==1) && (nAfter==0);
    fprintf('T3 draw/clear: patch=%d cleared=%d => %s\n', nPatch, (nAfter==0), PF{ok3+1});
    pass = pass && ok3;
    if ishandle(hFig), close(hFig); end

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (MATLAB, Brainstorm live):
```matlab
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/dev'); rehash; test_bst_geodesic_tool
```
Expected: FAIL/error — `bst_geodesic_tool` is undefined (`Unrecognized function or variable`).

- [ ] **Step 3: Create `bst_geodesic_tool.m`**

Create `toolbox/dynamics/bst_geodesic_tool.m`:

```matlab
function varargout = bst_geodesic_tool( varargin )
% BST_GEODESIC_TOOL: interactive heat-distance disk on the cortex (the dynamics source-axis
% point+extent primitive). A seed vertex (center) + a geodesic radius (extent) define a disk
% via the heat-distance field; the tool draws it as a TRANSIENT overlay (never a scout) and the
% atom panel snapshots it into an atom. Reuses the tess_scout_area engine (toolbox/anatomy).
%
% USAGE:
%   bst_geodesic_tool('Toggle', onoff)        % enter/exit the cortical-pick mode
%   tf  = bst_geodesic_tool('IsActive')
%   bst_geodesic_tool('OnClick', hFig)        % figure_3d: clicked-vertex -> Seed + Draw
%   ok  = bst_geodesic_tool('OnScroll', n)    % figure_3d: grow/shrink + Draw; returns consumed
%   bst_geodesic_tool('Seed', SurfaceFile, vi)% headless: compute the disk around vi
%   bst_geodesic_tool('Grow', scrollCount)    % headless: re-threshold at a new radius
%   st  = bst_geodesic_tool('GetState')       % struct(seed,phi,radius,vertices,pos,SurfaceFile[,hFig]) | []
%   bst_geodesic_tool('Draw', hFig) / ('Clear', hFig)
%
% SEE ALSO: tess_scout_area, panel_bst_dynamics
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

eval(macro_method);
end


%% ===== persistent state cache (one active disk) =====
function out = Cache(in)
    persistent C;
    if (nargin >= 1), C = in; end
    out = C;
end


%% ===== TOGGLE: enter/exit the dynamics cortical-pick mode =====
function Toggle(onoff) %#ok<DEFNU>
    hFigures = bst_figures('GetFiguresForScouts');
    if isempty(hFigures), hFigures = bst_figures('GetFigureWithSurfaces'); end
    if onoff
        if isempty(hFigures)
            java_dialog('warning', 'Open a 3D cortex figure first.', 'Region tool');  return;
        end
        for hFig = hFigures(:)'
            setappdata(hFig, 'isDynamicsGeodesicPick', 1);
            setappdata(hFig, 'isSelectingCorticalSpot', 0);   % mutual exclusion with scout pick
            set(hFig, 'Pointer', 'cross');
        end
        Cache([]);                                            % fresh session
        SurfaceFile = i_tool_surface(hFigures);
        if ~isempty(SurfaceFile)
            bst_progress('start', 'Region tool', 'Pre-factorizing the heat operator...');
            try, tess_scout_area('prewarm', SurfaceFile); catch, end %#ok<CTCH>
            bst_progress('stop');
        end
    else
        for hFig = hFigures(:)'
            setappdata(hFig, 'isDynamicsGeodesicPick', 0);
            set(hFig, 'Pointer', 'arrow');
            Clear(hFig);
        end
    end
end


%% ===== IS ACTIVE =====
function tf = IsActive() %#ok<DEFNU>
    tf = false;
    hFigures = bst_figures('GetFiguresForScouts');
    for hFig = hFigures(:)'
        if isappdata(hFig, 'isDynamicsGeodesicPick') && getappdata(hFig, 'isDynamicsGeodesicPick')
            tf = true;  return;
        end
    end
end


%% ===== SEED: compute the geodesic disk around vi (headless core) =====
function Seed(SurfaceFile, vi) %#ok<DEFNU>
    if isempty(SurfaceFile) || isempty(vi), return; end
    R0 = 0.003;     % initial geodesic radius [m] = 3 mm
    isProg = ~bst_progress('isVisible');
    if isProg, bst_progress('start', 'Region tool', 'Computing geodesic distance...'); end
    [verts, phi] = tess_scout_area(SurfaceFile, vi, R0);
    if isProg, bst_progress('stop'); end
    Surf = in_tess_bst(SurfaceFile, 0);
    Cache(struct('seed',double(vi), 'phi',phi, 'radius',R0, 'vertices',verts, ...
                 'pos',Surf.Vertices(vi,:), 'SurfaceFile',SurfaceFile));
end


%% ===== GROW: re-threshold cached phi at a new radius (headless core) =====
function Grow(scrollCount) %#ok<DEFNU>
    c = Cache();
    if isempty(c), return; end
    STEP = 0.003;   % 3 mm per scroll tick
    R = max(STEP, c.radius - double(scrollCount) * STEP);   % scroll up (<0) grows
    c.vertices = tess_scout_area(c.SurfaceFile, c.seed, R, c.phi);   % reuse cached distance
    c.radius = R;
    Cache(c);
end


%% ===== GET STATE =====
function st = GetState() %#ok<DEFNU>
    st = Cache();
end


%% ===== ON CLICK (figure_3d): resolve clicked vertex -> seed + draw =====
function OnClick(hFig) %#ok<DEFNU>
    TessInfo = getappdata(hFig, 'Surface');  iTess = getappdata(hFig, 'iSurface');
    if isempty(iTess) || isempty(TessInfo) || isempty(TessInfo(iTess).hPatch), return; end
    if strcmpi(TessInfo(iTess).Name, 'Anatomy'), return; end
    [~, vout, vi] = select3d(TessInfo(iTess).hPatch);
    if isempty(vout) || isempty(vi), return; end
    Seed(TessInfo(iTess).SurfaceFile, vi);
    c = Cache();  c.hFig = hFig;  Cache(c);    % remember the figure for OnScroll redraws
    Draw(hFig);
end


%% ===== ON SCROLL (figure_3d): grow/shrink + redraw; returns consumed flag =====
function handled = OnScroll(scrollCount) %#ok<DEFNU>
    handled = false;
    if ~IsActive(), return; end
    c = Cache();
    if isempty(c) || ~isfield(c,'hFig') || ~ishandle(c.hFig), return; end
    Grow(scrollCount);
    Draw(c.hFig);
    handled = true;
end


%% ===== DRAW the transient disk overlay =====
function Draw(hFig) %#ok<DEFNU>
    c = Cache();
    if isempty(c) || isempty(hFig) || ~ishandle(hFig), return; end
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');  if isempty(hAxes), return; end
    hAxes = hAxes(1);  set(hAxes, 'NextPlot', 'add');     % low-level: avoid the axes-reset trap
    Clear(hFig);
    Surf = in_tess_bst(c.SurfaceFile, 0);
    inR = false(size(Surf.Vertices,1), 1);  inR(c.vertices) = true;
    fIn = all(inR(Surf.Faces), 2);
    if ~any(fIn), return; end
    patch('Faces', Surf.Faces(fIn,:), 'Vertices', Surf.Vertices, 'Parent', hAxes, ...
        'FaceColor', [0.2 0.7 1.0], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'Tag', 'GeodesicToolDisk');
end


%% ===== CLEAR overlay =====
function Clear(hFig) %#ok<DEFNU>
    if isempty(hFig) || ~ishandle(hFig), return; end
    delete(findobj(hFig, 'Tag', 'GeodesicToolDisk'));
end


%% ===== cortex surface backing the figures (for prewarm) =====
function SurfaceFile = i_tool_surface(hFigures)
    SurfaceFile = '';
    for hFig = hFigures(:)'
        TessInfo = getappdata(hFig, 'Surface');
        if isempty(TessInfo), continue; end
        for i = 1:numel(TessInfo)
            sf = TessInfo(i).SurfaceFile;
            if ~isempty(sf) && strcmpi(file_gettype(sf), 'cortex'), SurfaceFile = sf;  return; end
        end
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```matlab
rehash; test_bst_geodesic_tool
```
Expected: `T1 seed: ... => PASS`, `T2 grow: ... => PASS`, `T3 draw/clear: patch=1 cleared=1 => PASS`, `==== SUITE: PASS ====`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/dynamics/bst_geodesic_tool.m dev/test_bst_geodesic_tool.m
git commit -m "feat(dynamics): bst_geodesic_tool — transient heat-disk source-axis primitive"
```

---

### Task 2: Re-home the figure interaction + excise the Area tool from Scouts

**Files:**
- Modify: `toolbox/gui/figure_3d.m` (click branch ~624; scroll branch ~943-947)
- Modify: `toolbox/gui/panel_scout.m` (remove the Area-tool button, functions, dispatch, controls entry, mutual-exclusion refs)
- Test: `dev/test_geodesic_move.m` (Scout panel loads, Line tool intact, dynamics scroll routed)

**Interfaces:**
- Consumes: `bst_geodesic_tool('OnClick'/'OnScroll'/'IsActive')` (Task 1).
- Produces: figure_3d routes cortical click → `bst_geodesic_tool('OnClick')` when `isDynamicsGeodesicPick` set, and scroll → `bst_geodesic_tool('OnScroll')`; `panel_scout` no longer defines/uses the Area tool.

- [ ] **Step 1: Write the failing test**

Create `dev/test_geodesic_move.m`:

```matlab
function test_geodesic_move()
% TEST_GEODESIC_MOVE: the Area tool is gone from Scouts; the dynamics tool drives the scroll.
%
% USAGE:  test_geodesic_move   % Brainstorm running
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};  pass = true;

    % T1: panel_scout no longer exposes the Area tool verbs (removed), Line tool stays
    haveArea = ~isempty(which('panel_scout'));
    areaGone = true;
    try, panel_scout('IsAreaToolActive');  areaGone = false; catch, areaGone = true; end %#ok<CTCH>
    lineStays = true;
    try, panel_scout('IsGeodesicToolActive'); catch, lineStays = false; end %#ok<CTCH>
    ok1 = haveArea && areaGone && lineStays;
    fprintf('T1 scout surgery: areaVerbGone=%d lineToolStays=%d => %s\n', areaGone, lineStays, PF{ok1+1});
    pass = pass && ok1;

    % T2: figure_3d source no longer references the removed scout Area verbs, and references the dynamics tool
    src = fileread(which('figure_3d'));
    noScoutArea = isempty(strfind(src, 'AreaToolScroll')) && isempty(strfind(src, 'IsAreaToolActive')); %#ok<STREMP>
    hasDynTool  = ~isempty(strfind(src, 'bst_geodesic_tool')); %#ok<STREMP>
    ok2 = noScoutArea && hasDynTool;
    fprintf('T2 figure_3d rewire: noScoutAreaRefs=%d hasDynTool=%d => %s\n', noScoutArea, hasDynTool, PF{ok2+1});
    pass = pass && ok2;

    % T3: with the tool OFF, OnScroll passes through (returns false) so the wheel still zooms
    bst_geodesic_tool('Toggle', 0);                 % ensure inactive
    passthrough = bst_geodesic_tool('OnScroll', -1);
    ok3 = islogical(passthrough) && (passthrough == false);
    fprintf('T3 scroll passthrough (tool off): consumed=%d => %s\n', passthrough, PF{ok3+1});
    pass = pass && ok3;

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```matlab
rehash; test_geodesic_move
```
Expected: FAIL — T1 `panel_scout('IsAreaToolActive')` still succeeds (Area tool not yet removed) so `areaGone=0`; T2 `figure_3d` still contains `AreaToolScroll`/`IsAreaToolActive` so `noScoutAreaRefs=0`.

- [ ] **Step 3: Edit `figure_3d.m` — click branch**

In `toolbox/gui/figure_3d.m`, the cortical-spot click branch (~line 622-624) reads:
```matlab
        % === SELECTING CORTICAL SCOUTS (New scout or geodesic Area tool) ===
        elseif isSelectingCorticalSpot
            panel_scout('CreateScoutMouse', hFig);
```
Replace with (dynamics pick takes precedence):
```matlab
        % === DYNAMICS GEODESIC REGION TOOL (heat-disk pick) ===
        elseif getappdata(hFig, 'isDynamicsGeodesicPick')
            bst_geodesic_tool('OnClick', hFig);
        % === SELECTING CORTICAL SCOUTS (New scout, geodesic Line tool) ===
        elseif isSelectingCorticalSpot
            panel_scout('CreateScoutMouse', hFig);
```

- [ ] **Step 4: Edit `figure_3d.m` — scroll branch**

In `toolbox/gui/figure_3d.m`, the scroll branch (~lines 939-947) reads:
```matlab
    % Scout geodesic Area tool: while the tool is active, the wheel grows/shrinks the SELECTED
    % scout (heat-distance disk) instead of zooming -- gated on the toggle being pressed and a
    % scout being selected (not on a keyboard modifier, which is unreliable once the scout list
    % has focus). With the tool off, or with no scout selected, the wheel zooms as usual.
    if ~isempty(event) && panel_scout('IsAreaToolActive')
        if panel_scout('AreaToolScroll', double(event.VerticalScrollCount))
            return;
        end
    end
```
Replace with:
```matlab
    % Dynamics geodesic Region tool: while active, the wheel grows/shrinks the heat-disk overlay
    % instead of zooming. OnScroll returns true only when the tool is active and a disk is seeded;
    % otherwise the wheel zooms as usual.
    if ~isempty(event) && bst_geodesic_tool('OnScroll', double(event.VerticalScrollCount))
        return;
    end
```

- [ ] **Step 5: Remove the Area tool from `panel_scout.m`**

Make these deletions (the geodesic **Line** tool and `i_get_scout_surface` stay):

1. **Button (lines ~180-187):** delete the `jToggleArea` block (the `% Geodesic area tool` comment through `jToggleArea.setVerticalAlignment(...)`), i.e. from the line `% Geodesic area tool (own row, below the swell buttons)` through the `jToggleArea.setVerticalAlignment(javax.swing.SwingConstants.CENTER);` line. Leave the following `% Geodesic line tool` block intact, and change its `gui_component('toggle', jPanelScoutOptions, '', ...)` first positional arg from `''` to `'br'` so the Line button starts its own row now that the Area button (which had `'br'`) is gone.
2. **Controls struct (line ~242):** delete the line `'jToggleArea',           jToggleArea, ...`.
3. **Tool-release block (lines ~1914-1917):** this block deselects the Area toggle on some action; delete the `if isfield(ctrl,'jToggleArea') && ~isempty(ctrl.jToggleArea)` / `ctrl.jToggleArea.setSelected(0);` / `end` (the Area-specific deselect). Keep any sibling `jToggleGeo` deselect.
4. **`CreateScoutMouse` Area dispatch (lines ~2982-2986):** delete the block
   ```matlab
       % ===== GEODESIC AREA TOOL =====
       % ... comment ...
       if ~isVolumeAtlas && IsAreaToolActive()
           CreateScoutArea(vi, TessInfo(iTess).SurfaceFile);
           return;
       end
   ```
   Keep the `% ===== GEODESIC LINE TOOL =====` block that follows.
5. **Functions:** delete `IsAreaToolActive` (~3058-3064), `AreaToolToggle` (~3069-3111), `AreaCache` (~3113-3120), `CreateScoutArea` (~3122-3144), `AreaToolScroll` (~3146-3193). Keep `i_get_scout_surface` (~3196-3213) — the Line tool uses it.
6. **GeodesicToolToggle mutual-exclusion (line ~3254):** delete the line
   `if isfield(ctrl,'jToggleArea') && ~isempty(ctrl.jToggleArea), ctrl.jToggleArea.setSelected(0); end`
   (no Area toggle to deselect anymore).

After editing, confirm no `jToggleArea`, `AreaToolToggle`, `IsAreaToolActive`, `AreaToolScroll`, `AreaCache`, or `CreateScoutArea` token remains in `panel_scout.m`:
```bash
grep -n "jToggleArea\|AreaToolToggle\|IsAreaToolActive\|AreaToolScroll\|AreaCache\|CreateScoutArea" toolbox/gui/panel_scout.m
```
Expected: no output.

- [ ] **Step 6: Run test to verify it passes**

Run:
```matlab
rehash; test_geodesic_move
```
Expected: `T1 scout surgery: areaVerbGone=1 lineToolStays=1 => PASS`, `T2 figure_3d rewire: noScoutAreaRefs=1 hasDynTool=1 => PASS`, `T3 scroll route: ... => PASS`, `==== SUITE: PASS ====`.

- [ ] **Step 7: Smoke-check the Scout panel still loads + Line tool toggles**

Run:
```matlab
gui_show('panel_scout', 'BrainstormTab', 'tools');  drawnow;
ctrl = bst_get('PanelControls', 'Scout');
fprintf('scout panel loaded=%d, hasLineToggle=%d, hasAreaToggle=%d\n', ~isempty(ctrl), isfield(ctrl,'jToggleGeo'), isfield(ctrl,'jToggleArea'));
```
Expected: `scout panel loaded=1, hasLineToggle=1, hasAreaToggle=0`.

- [ ] **Step 8: Commit**

```bash
git add toolbox/gui/figure_3d.m toolbox/gui/panel_scout.m dev/test_geodesic_move.m
git commit -m "refactor(scout): move geodesic Area tool interaction to bst_geodesic_tool; figure_3d routes to it"
```

---

### Task 3: Rewire atom Capture to the dynamics tool + Region toggle

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (`OnCaptureRegion` ~line 394; add a Region-tool toggle in the Record row + controls struct)
- Test: `dev/test_dynamics_atoms.m` (add a Capture-via-tool integration check; T1–T7 stay green)

**Interfaces:**
- Consumes: `bst_geodesic_tool('GetState'/'Toggle')` (Task 1); `bst_dynamics('AttachRegion', …)`; `file_compare`.
- Produces: `OnCaptureRegion` snapshots `bst_geodesic_tool('GetState')` into the active atom occurrence; a Region-tool toggle button activates the dynamics tool from the panel.

- [ ] **Step 1: Write the failing test**

In `dev/test_dynamics_atoms.m`, after the T7 block (before the `% cleanup` block), add T8:

```matlab
    % T8: OnCaptureRegion reads bst_geodesic_tool('GetState') (Scout-decoupled capture)
    st = getappdata(0,'DynamicsTarget');  Tt = st.T;
    SurfT8 = in_tess_bst(Tt.SurfaceFile, 0);  seed8 = round(size(SurfT8.Vertices,1)/4);
    bst_geodesic_tool('Seed', Tt.SurfaceFile, seed8);          % seed the dynamics tool (no scout)
    gs8 = bst_geodesic_tool('GetState');
    % select occurrence 1 of the trough phase child as the active atom
    iWinNode8 = find(arrayfun(@(k) strcmp(st.nodeInfo(k).kind,'window'), 1:numel(st.nodeInfo)), 1);
    ctrl = bst_get('PanelControls','Dynamics');
    ctrl.jTree.setSelectionPath(javax.swing.tree.TreePath(st.nodeList{iWinNode8}.getPath()));  drawnow;
    % find a list row that is a trough occurrence and select it
    model = ctrl.jListOccur.getModel();  rowSel = -1;
    for r = 0:(model.getSize()-1)
        if ~isempty(strfind(char(model.getElementAt(r)), 'trough')), rowSel = r; break; end
    end
    nG8 = 0;
    if (rowSel >= 0)
        ctrl.jListOccur.setSelectedIndex(rowSel);  drawnow;
        st = getappdata(0,'DynamicsTarget');  row = ctrl.jListOccur.getSelectedIndex()+1;
        g8 = st.occMap(row,1);  o8 = st.occMap(row,2);
        panel_bst_dynamics('OnCaptureRegion');  drawnow;
        st = getappdata(0,'DynamicsTarget');
        nG8 = numel(st.T.Groups(g8).region);  reg8 = st.T.Groups(g8).region{o8};
        capOK = ~isempty(reg8) && isequal(double(reg8), double(gs8.vertices)) && (st.T.Groups(g8).vertices(o8)==seed8);
    else
        capOK = false;
    end
    ok8 = ~isempty(gs8) && (rowSel>=0) && capOK;
    fprintf('T8 capture-via-tool: seeded=%d rowSel=%d regionWritten=%d => %s\n', ~isempty(gs8), (rowSel>=0), capOK, PF{ok8+1});
    pass = pass && ok8;
    bst_geodesic_tool('Clear', hFig);
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```matlab
rehash; test_dynamics_atoms
```
Expected: FAIL at T8 — `OnCaptureRegion` still calls `panel_scout('GetSelectedScouts')` (no scout selected → it warns and returns without writing a region), so `capOK=0`.

- [ ] **Step 3: Rewire `OnCaptureRegion`**

In `toolbox/gui/panel_bst_dynamics.m`, `OnCaptureRegion`, the current scout read + guards are:
```matlab
    % the geodesic region = the currently selected Scout (grown with the Scout "Area" tool)
    [sScout, ~, sSurf] = panel_scout('GetSelectedScouts');
    if isempty(sScout) || isempty(sScout(1).Vertices)
        java_dialog('warning', 'Grow a region with the Scout "Area" tool first.', 'Capture region');  return;
    end
    sScout = sScout(1);
    if ~isempty(sSurf) && ~isempty(SurfaceFile) && ~file_compare(sSurf.FileName, SurfaceFile)
        java_dialog('warning', 'The selected region is on a different surface than the atoms.', 'Capture region');  return;
    end
    if ~isempty(sScout.Seed), seed = double(sScout.Seed(1)); else, seed = double(sScout.Vertices(1)); end
    Surf = getappdata(st.hFig, 'DynamicsSurf');
    if isempty(Surf), Surf = in_tess_bst(SurfaceFile, 0);  setappdata(st.hFig, 'DynamicsSurf', Surf); end
    pos  = Surf.Vertices(seed, :);
    hemi = 1 + (pos(2) < 0);                                       % SCS Y>0 = left
    st.T.Groups(g) = bst_dynamics('AttachRegion', st.T.Groups(g), o, sScout.Vertices, seed, pos, hemi);
```
Replace that block (from the `% the geodesic region` comment through the `AttachRegion` line) with:
```matlab
    % the geodesic region = the dynamics Region tool's current heat-disk (no scout)
    gs = bst_geodesic_tool('GetState');
    if isempty(gs) || isempty(gs.vertices)
        java_dialog('warning', 'Seed a region with the Region tool first.', 'Capture region');  return;
    end
    if ~isempty(gs.SurfaceFile) && ~isempty(SurfaceFile) && ~file_compare(gs.SurfaceFile, SurfaceFile)
        java_dialog('warning', 'The region is on a different surface than the atoms.', 'Capture region');  return;
    end
    seed = double(gs.seed);
    pos  = gs.pos;
    hemi = 1 + (pos(2) < 0);                                       % SCS Y>0 = left
    st.T.Groups(g) = bst_dynamics('AttachRegion', st.T.Groups(g), o, gs.vertices, seed, pos, hemi);
```

- [ ] **Step 4: Add the Region-tool toggle to the Record row**

In `toolbox/gui/panel_bst_dynamics.m`, `CreatePanel`, the Record row currently has the two buttons:
```matlab
    gui_component('button', jRec, 'tab hfill', 'Record at cursor', [], 'Store the shaped field''s extrema at the cursor time as atoms', @(h,e)bst_call(@OnRecord));
    gui_component('button', jRec, 'br hfill', 'Capture region -> active atom', [], 'Snapshot the selected Scout''s vertices into the selected atom (localizes a time-only marker)', @(h,e)bst_call(@OnCaptureRegion));
```
Insert a Region-tool toggle before the Capture button and update the Capture tooltip:
```matlab
    gui_component('button', jRec, 'tab hfill', 'Record at cursor', [], 'Store the shaped field''s extrema at the cursor time as atoms', @(h,e)bst_call(@OnRecord));
    jRegionTool = gui_component('toggle', jRec, 'br', 'Region tool', [], 'Heat-disk tool: click a cortex vertex to seed a region, scroll to grow/shrink it', @(h,e)bst_call(@()bst_geodesic_tool('Toggle', ctrl_region_state())));
    gui_component('button', jRec, 'tab hfill', 'Capture region -> active atom', [], 'Snapshot the Region tool''s heat-disk into the selected atom (localizes a time-only marker)', @(h,e)bst_call(@OnCaptureRegion));
```
Add `jRegionTool` to the `BstPanel(...)` controls struct (extend the final struct argument list, after `'jPhaseItems',jPhaseItems`):
```matlab
        'jSpaceParams',jSpaceParams, 'jSpacePot',jSpacePot, 'jSpaceStr',jSpaceStr, 'jPeaks',jPeaks, 'jPhaseItems',jPhaseItems, 'jRegionTool',jRegionTool));
```
Add this small helper (place after `OnCaptureRegion`) so the toggle passes its pressed state:
```matlab
% Region-tool toggle state (1 when pressed) -> bst_geodesic_tool('Toggle', state)
function s = ctrl_region_state()
    s = 0;
    ctrl = bst_get('PanelControls', 'Dynamics');
    if ~isempty(ctrl) && isfield(ctrl,'jRegionTool') && ~isempty(ctrl.jRegionTool)
        s = double(ctrl.jRegionTool.isSelected());
    end
end
```

- [ ] **Step 5: Run test to verify it passes**

Run:
```matlab
rehash; test_dynamics_atoms
```
Expected: `T8 capture-via-tool: seeded=1 rowSel=1 regionWritten=1 => PASS`, and T1–T7 still PASS, ending `==== SUITE: PASS ====`.

- [ ] **Step 6: Confirm the geodesic-tool suite still passes**

Run:
```matlab
test_bst_geodesic_tool
test_geodesic_move
```
Expected: both end `==== SUITE: PASS ====`.

- [ ] **Step 7: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m dev/test_dynamics_atoms.m
git commit -m "feat(dynamics): atom Capture reads bst_geodesic_tool; Region-tool toggle replaces Scout dependency"
```

---

## Self-Review

**1. Spec coverage:**
- New `bst_geodesic_tool` (transient disk, reuses `tess_scout_area`, never a scout) → Task 1.
- figure_3d two parallel branches (click + scroll) → Task 2 Steps 3-4.
- Excise Area tool from `panel_scout`; keep Line tool + `i_get_scout_surface` + `tess_scout_area` → Task 2 Step 5.
- `AreaToolScroll` removal coupled with the figure_3d scroll edit (same task) → Task 2.
- Capture rewired off `GetSelectedScouts` → Task 3 Step 3.
- Region-tool toggle in the panel → Task 3 Step 4.
- Mutual exclusion (`Toggle` clears `isSelectingCorticalSpot`) → Task 1 `Toggle`.
- Constants R0/STEP = 0.003 → Task 1 `Seed`/`Grow`.
- Tests: tool unit (Task 1), scout-surgery + figure rewire (Task 2), Capture integration + T1–T7 green (Task 3).
- Out of scope (correctly deferred): source navigator block; renaming `tess_scout_area`; moving the Line tool; soft/wavelet disk.

**2. Placeholder scan:** none — every code step is complete; every run step has an exact command + expected line.

**3. Type consistency:** `GetState` returns `struct(seed,phi,radius,vertices,pos,SurfaceFile[,hFig])` — the same fields are read in `OnCaptureRegion` (Task 3: `gs.vertices/seed/pos/SurfaceFile`) and the tests (`st.vertices/seed/pos/SurfaceFile`). Verb names (`Toggle/IsActive/Seed/Grow/GetState/Draw/Clear/OnClick/OnScroll`) match between the module, figure_3d, the panel, and the tests. `isDynamicsGeodesicPick` appdata key is identical in `Toggle`/`IsActive`/`OnClick`-gate (figure_3d) and `test_geodesic_move`. Tag `'GeodesicToolDisk'` matches between `Draw`, `Clear`, and the Task 1 test.
