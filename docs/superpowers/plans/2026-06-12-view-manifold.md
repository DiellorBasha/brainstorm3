# view_manifold Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `view_manifold`, a consolidated viewer that renders a manifold node's canonical per-vertex tangent frame on the cortex; replace `view_tangents`; repoint `view_connection_phase` and the Manifold node menu to the manifold; deprecate (keep) `tess_tangents`.

**Architecture:** `view_manifold(ManifoldFile)` loads the `manifold_*.mat` node, derives the per-vertex orthonormal frame from its complex `Embedded.vertex.grid` and `Gauge.vertex.rotation` (`U=real(grid.*rot)`, `V=imag(grid.*rot)`, `N=cross(U,V)`), and draws headless `quiver3` glyphs on an opaque cortex — modeled on `view_tangents`. A `'View'` dispatch (default `'frame'`) leaves room for future view modes.

**Tech Stack:** MATLAB R2023b, Brainstorm dev fork. Live tests run via the MATLAB MCP against a running session with the registered manifold node. Headless tests use the `eval(macro_method)` dispatch (`view_manifold('DeriveVertexFrame', ...)`).

**Spec:** `docs/superpowers/specs/2026-06-12-view-manifold-design.md`

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `toolbox/gui/view_manifold.m` | Consolidated manifold viewer + pure `DeriveVertexFrame` | **Create** |
| `dev/tests/test_manifold_frame.m` | Headless unit test for `DeriveVertexFrame` | **Create** |
| `dev/tests/test_view_manifold.m` | Live-figure integration test | **Create** |
| `toolbox/tree/tree_callbacks.m` | Manifold node → single "View manifold" item | **Modify** |
| `toolbox/gui/view_connection_phase.m` | Repoint frame source to the manifold (`tess_frame`) | **Modify** |
| `toolbox/anatomy/tess_tangents.m` | `@deprecated` header note | **Modify** |
| `toolbox/gui/view_tangents.m` | Replaced by `view_manifold` | **Delete** |
| `dev/tests/test_view_tangents.m` | Tests deleted viewer | **Delete** |

**Kept untouched:** `tess_tangents.m` body (still serves `bst_wavefront_track`, `tess_nxr_populate`), `tess_frame.m`, `tess_manifold.m`.

---

## Task 1: view_manifold core + headless frame test

**Files:**
- Create: `toolbox/gui/view_manifold.m`
- Create: `dev/tests/test_manifold_frame.m`

**Context:** A manifold node (`db_template('manifoldmat')`) stores 1×2 per-hemisphere `Embedded` and `Gauge` structs. `Embedded(hh).vertex.grid` is `[nVh x 3] complex` (encodes `U + iV`); `Gauge(hh).vertex.rotation` is `[nVh x 1] complex unit`. The frame is `cRot = grid.*rot; U=real(cRot); V=imag(cRot); N=cross(U,V)` — verified on the canonical cortex: `U,V` unit, orthogonal, ⟂ `N`, `cross(U,V)=normal` to 1e-16. `Embedded(hh).GlobalVertices` scatters hemisphere-local rows into the full vertex grid; `Embedded(hh).vertex.position` are the glyph anchors; `Gauge(hh).singularity.vertices` are local pole indices.

**CRITICAL — axes-reset trap:** `quiver3` on the `'Axes3D'` axes runs `newplot` and resets it (deleting the cortex). Call `hold(hAxes,'on')` before any `quiver3` (as `view_tangents` does).

- [ ] **Step 1: Write the failing headless test**

Create `dev/tests/test_manifold_frame.m`:

```matlab
function test_manifold_frame()
% TEST_MANIFOLD_FRAME  Headless regression for view_manifold's DeriveVertexFrame.
% Requires Brainstorm on path so view_manifold('DeriveVertexFrame', ...) dispatches.
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;

    % 2 hemispheres, nVh=2 each. Global vertices: L->[1;2], R->[3;4]; nVert=4.
    % Per vertex, grid = e1 + 1i*e2 for chosen orthonormal (e1,e2); rot = 1.
    %   L v1: e1=[1 0 0] e2=[0 1 0] -> N=[0 0 1]
    %   L v2: e1=[0 1 0] e2=[0 0 1] -> N=[1 0 0]
    %   R v3: e1=[0 0 1] e2=[1 0 0] -> N=[0 1 0]
    %   R v4: e1=[1 0 0] e2=[0 0 1] -> N=[0 -1 0]
    Emb(1).GlobalVertices = [1;2];
    Emb(1).vertex.grid     = [1 0 0; 0 1 0] + 1i*[0 1 0; 0 0 1];
    Emb(1).vertex.position = [10 0 0; 11 0 0];
    Emb(2).GlobalVertices = [3;4];
    Emb(2).vertex.grid     = [0 0 1; 1 0 0] + 1i*[1 0 0; 0 0 1];
    Emb(2).vertex.position = [0 10 0; 0 11 0];
    Ga(1).vertex.rotation = [1;1];
    Ga(1).singularity.vertices = uint32(1);          % local -> global 1
    Ga(2).vertex.rotation = [1;1];
    Ga(2).singularity.vertices = uint32([]);

    G = view_manifold('DeriveVertexFrame', Emb, Ga, 4);
    [nPass,nFail] = chk('U v1', isequal(G.U(1,:),[1 0 0]), nPass,nFail);
    [nPass,nFail] = chk('V v1', isequal(G.V(1,:),[0 1 0]), nPass,nFail);
    [nPass,nFail] = chk('N v1 = cross(U,V)', max(abs(G.N(1,:)-[0 0 1]))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('N v2', max(abs(G.N(2,:)-[1 0 0]))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('U v3 (hemi R scatter)', isequal(G.U(3,:),[0 0 1]), nPass,nFail);
    [nPass,nFail] = chk('anchor P v3 = position', isequal(G.P(3,:),[0 10 0]), nPass,nFail);
    [nPass,nFail] = chk('orthonormal U.V=0', max(abs(sum(G.U(1:4,:).*G.V(1:4,:),2)))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('unit |U|', max(abs(sqrt(sum(G.U(1:4,:).^2,2))-1))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('singularity globalized', isequal(G.Sing(:), 1), nPass,nFail);

    % off-support vertices stay zero (nVert larger than mapped indices)
    G2 = view_manifold('DeriveVertexFrame', Emb, Ga, 6);
    [nPass,nFail] = chk('off-support U zero', isequal(G2.U(5,:),[0 0 0]) && isequal(G2.U(6,:),[0 0 0]), nPass,nFail);

    % rotation rotates the in-plane frame: rot=1i -> U=-e2, V=e1 on L v1
    Ga2 = Ga; Ga2(1).vertex.rotation = [1i;1i];
    Gr = view_manifold('DeriveVertexFrame', Emb, Ga2, 4);
    [nPass,nFail] = chk('rot=1i -> U=-e2', max(abs(Gr.U(1,:)-[0 -1 0]))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('rot=1i -> V=e1',  max(abs(Gr.V(1,:)-[1 0 0]))<1e-12, nPass,nFail);

    % shape mismatch errors
    err=false; Bad=Emb; Bad(1).vertex.grid=[1 0 0]; try, view_manifold('DeriveVertexFrame', Bad, Ga, 4); catch, err=true; end
    [nPass,nFail] = chk('shape mismatch errors', err, nPass,nFail);

    fprintf('\n==== test_manifold_frame: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_manifold_frame: %d test(s) FAILED.', nFail); end
end

function [p,f] = chk(name, c, p, f)
    if c, fprintf('  PASS  %s\n', name); p=p+1; else, fprintf('  FAIL  %s\n', name); f=f+1; end
end
```

- [ ] **Step 2: Run the test to verify it fails**

In the MATLAB MCP session: `addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests'); test_manifold_frame`
Expected: FAIL — `view_manifold` does not exist yet (`Unrecognized function`).

- [ ] **Step 3: Create `toolbox/gui/view_manifold.m`**

```matlab
function varargout = view_manifold(varargin)
% VIEW_MANIFOLD: Display the per-vertex frame stored in a manifold node.
%
% USAGE:  hFig = view_manifold(ManifoldFile)
%         hFig = view_manifold(ManifoldFile, 'View','frame', 'MaxArrows', 3000)
%         G    = view_manifold('DeriveVertexFrame', Embedded, Gauge, nVert)
%
% Reads a manifold_ DB node (db_template('manifoldmat')) and renders its canonical
% per-vertex tangent frame on the parent cortex. The frame is derived from the
% complex grid and gauge rotation stored in the node:
%     cRot = Embedded(hh).vertex.grid .* Gauge(hh).vertex.rotation
%     U = real(cRot);  V = imag(cRot);  N = cross(U,V)
% (the same convention as tess_frame). Per-vertex only; the per-face frame needs
% Gauge.face.rotation, which nxr defers for the trivial gauge. The 'View' option
% selects the manifold view (default 'frame'); more views are added later.
%
% Keyboard (figure focused):
%   Left/Right        fewer / more frames
%   Shift + Up/Down   glyph length
%   Ctrl  + Up/Down   line width
%   N                 toggle vertex-normal glyphs (off by default)
%   P                 toggle singularity markers
%   H                 help
%
% SEE ALSO: tess_manifold, tess_frame, view_surface
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

if (nargin >= 1) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'DeriveVertexFrame'}))
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: per-vertex frame from the manifold node =====
function G = DeriveVertexFrame(Embedded, Gauge, nVert)
% G: struct with P,U,V,N [nVert x 3] (zeros off-support) and Sing (global vtx ids).
    if numel(Embedded) ~= 2 || numel(Gauge) ~= 2
        error('view_manifold:badNode', 'Embedded and Gauge must be 1x2 per-hemisphere structs.');
    end
    P = zeros(nVert,3); U = zeros(nVert,3); V = zeros(nVert,3); N = zeros(nVert,3);
    Sing = [];
    for hh = 1:2
        vH   = double(Embedded(hh).GlobalVertices(:));
        grid = Embedded(hh).vertex.grid;
        rot  = Gauge(hh).vertex.rotation;
        if size(grid,1) ~= numel(vH) || size(grid,2) ~= 3
            error('view_manifold:shapeMismatch', ...
                'Hemisphere %d: Embedded.vertex.grid is [%s], expected [%d x 3].', ...
                hh, num2str(size(grid)), numel(vH));
        end
        cRot = grid .* rot(:);
        Uh = real(cRot);  Vh = imag(cRot);
        U(vH,:) = Uh;
        V(vH,:) = Vh;
        N(vH,:) = cross(Uh, Vh, 2);
        P(vH,:) = Embedded(hh).vertex.position;
        if isfield(Gauge(hh),'singularity') && isstruct(Gauge(hh).singularity) ...
                && isfield(Gauge(hh).singularity,'vertices') && ~isempty(Gauge(hh).singularity.vertices)
            Sing = [Sing; vH(double(Gauge(hh).singularity.vertices(:)))]; %#ok<AGROW>
        end
    end
    G = struct('P',P, 'U',U, 'V',V, 'N',N, 'Sing',Sing);
end


%% ===== GUI: render the frame view =====
function hFig = ViewFigure(ManifoldFile, varargin)
    hFig = [];
    View = 'frame';
    MaxArrows = [];
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'view',      View = lower(varargin{i+1});
            case 'maxarrows', MaxArrows = varargin{i+1};
        end
    end
    % --- load + validate node ---
    Full = file_fullpath(ManifoldFile);
    if ~file_exist(Full)
        bst_error('Manifold file not found.', 'View manifold', 0);
        return;
    end
    M = load(Full);
    if ~isfield(M,'ParentSurface') || isempty(M.ParentSurface) ...
            || ~isfield(M,'Embedded') || numel(M.Embedded) ~= 2 ...
            || ~isfield(M,'Gauge') || numel(M.Gauge) ~= 2
        bst_error(['Not a valid manifold node ' 10 '(need ParentSurface, Embedded[1x2], Gauge[1x2]).'], 'View manifold', 0);
        return;
    end
    switch View
        case 'frame'
            % implemented below
        otherwise
            bst_error(sprintf('Unknown manifold view: %s', View), 'View manifold', 0);
            return;
    end

    Surface = M.ParentSurface;
    TessMat = in_tess_bst(Surface);
    nVert   = size(TessMat.Vertices, 1);

    % --- derive frame ---
    G = DeriveVertexFrame(M.Embedded, M.Gauge, nVert);
    nP = nVert;
    meanEdge = MeanEdgeLength(TessMat.Vertices, double(TessMat.Faces));

    % singularity lollipop endpoints (radial lift clears sulci)
    SingBase = []; SingTip = [];
    if ~isempty(G.Sing)
        SingBase = G.P(G.Sing, :);
        ctrAll   = mean(G.P, 1);
        radial   = SingBase - ctrAll;
        radial   = radial ./ max(sqrt(sum(radial.^2,2)), eps);
        liftLen  = 0.10 * max(max(G.P,[],1) - min(G.P,[],1));
        SingTip  = SingBase + liftLen .* radial;
    end

    % --- open surface figure (opaque, wireframe) ---
    hFig = view_surface(Surface, 0, [.5 .5 .5], 'NewFigure', 0);
    if isempty(hFig)
        bst_error('Could not open the surface figure.', 'View manifold', 0);
        return;
    end
    set(hFig, 'Name', ['Manifold (frame): ' Surface]);
    panel_surface('SetSurfaceSmooth', hFig, 1, 0, 0);
    panel_surface('SetSurfaceEdges',  hFig, 1, 1);
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    hold(hAxes, 'on');   % CRITICAL: prevents quiver3 from resetting the cortex axes

    % --- state ---
    if isempty(MaxArrows), nArrows = round(nP/2); else, nArrows = min(MaxArrows, nP); end
    quiverSize  = 1; quiverWidth = 1; showNormals = false; showSing = true;
    colTangent  = [1 1 0]; colNormal = [1 0 1]; colSing = [0 0.45 1];

    KeyPressFcn_bak = get(hFig, 'KeyPressFcn');
    set(hFig, 'KeyPressFcn', @KeyPress_Callback);
    hLabel = uicontrol('Style','text','String','...','Units','Pixels', ...
        'Position',[6 4 1000 18],'HorizontalAlignment','left', ...
        'FontUnits','points','FontSize',bst_get('FigFont'), ...
        'ForegroundColor',[1 1 1],'BackgroundColor',[0 0 0],'Parent',hFig);

    DrawArrows();

    % ===== NESTED: draw frame glyphs =====
    function DrawArrows()
        delete(findobj(hAxes,'-depth',1,'-regexp','Tag','^tangent'));
        idx = ArrowSubsample(nP, nArrows);
        [B,Uv,Vv,Nv] = ArrowField(G.P, G.N, G.U, G.V, idx, quiverSize*0.5*meanEdge, 0.1*meanEdge);
        quiver3(B(:,1),B(:,2),B(:,3), Uv(:,1),Uv(:,2),Uv(:,3), 0, 'Parent',hAxes, ...
            'Color',colTangent,'LineWidth',quiverWidth,'ShowArrowHead','off','Tag','tangentU');
        quiver3(B(:,1),B(:,2),B(:,3), Vv(:,1),Vv(:,2),Vv(:,3), 0, 'Parent',hAxes, ...
            'Color',colTangent,'LineWidth',quiverWidth,'ShowArrowHead','off','Tag','tangentV');
        if showNormals
            quiver3(B(:,1),B(:,2),B(:,3), Nv(:,1),Nv(:,2),Nv(:,3), 0, 'Parent',hAxes, ...
                'Color',colNormal,'LineWidth',quiverWidth,'ShowArrowHead','off','Tag','tangentN');
        end
        hSq = [];
        if showSing && ~isempty(SingTip)
            line([SingBase(:,1)';SingTip(:,1)'],[SingBase(:,2)';SingTip(:,2)'],[SingBase(:,3)';SingTip(:,3)'], ...
                'Parent',hAxes,'Color',colSing,'LineWidth',1.5,'Tag','tangentSingStem');
            hSq = plot3(SingTip(:,1),SingTip(:,2),SingTip(:,3),'o','Parent',hAxes, ...
                'MarkerFaceColor',colSing,'MarkerEdgeColor',[.2 .2 .2],'MarkerSize',10,'LineStyle','none','Tag','tangentSing');
        end
        hTanProxy = line(NaN,NaN,'Parent',hAxes,'Color',colTangent,'LineStyle','-','LineWidth',2,'Tag','tangentLegProxy');
        legH = hTanProxy; legL = {'Tangent frame (U,V)'};
        if showNormals
            hNorProxy = line(NaN,NaN,'Parent',hAxes,'Color',colNormal,'LineStyle','-','LineWidth',2,'Tag','tangentLegProxy');
            legH = [legH, hNorProxy]; legL{end+1} = 'Vertex normal';
        end
        if ~isempty(hSq), legH = [legH, hSq]; legL{end+1} = 'Singularity'; end
        legend(legH, legL, 'TextColor',[1 1 1],'Color',[0 0 0],'Location','NorthEast','Interpreter','none','Tag','tangentLegend');
        set(hLabel,'String',sprintf('%d / %d vertices  |  N: normals   P: poles   <-/->: density   H: help', numel(idx), nP));
    end

    % ===== NESTED: keyboard =====
    function KeyPress_Callback(h, ev)
        switch ev.Key
            case 'rightarrow', nArrows = min(nP, ceil(nArrows*1.5));
            case 'leftarrow',  nArrows = max(min(50,nP), floor(nArrows/1.5));
            case 'uparrow'
                if ismember('shift',ev.Modifier),       quiverSize  = quiverSize*1.2;
                elseif ismember('control',ev.Modifier),  quiverWidth = quiverWidth*1.2;
                else, KeyPressFcn_bak(h,ev); return; end
            case 'downarrow'
                if ismember('shift',ev.Modifier),       quiverSize  = max(0.05, quiverSize/1.2);
                elseif ismember('control',ev.Modifier),  quiverWidth = max(0.5, quiverWidth/1.2);
                else, KeyPressFcn_bak(h,ev); return; end
            case 'n', showNormals = ~showNormals;
            case 'p', showSing = ~showSing;
            case 'h'
                java_dialog('msgbox', ['<HTML><TABLE>' ...
                    '<TR><TD><B>Left / Right</B></TD><TD>Fewer / more frames</TD></TR>' ...
                    '<TR><TD><B>Shift + Up/Down</B></TD><TD>Glyph length</TD></TR>' ...
                    '<TR><TD><B>Ctrl + Up/Down</B></TD><TD>Line width</TD></TR>' ...
                    '<TR><TD><B>N</B></TD><TD>Toggle vertex-normal glyphs</TD></TR>' ...
                    '<TR><TD><B>P</B></TD><TD>Toggle singularity markers</TD></TR>' ...
                    '</TABLE>'], 'Manifold frame shortcuts', [], 0);
                return;
            otherwise, KeyPressFcn_bak(h,ev); return;
        end
        DrawArrows();
    end
end


%% ========================================================================
function idx = ArrowSubsample(nP, nArrows)
% Deterministic uniform stride over vertex index (reproducible for tests).
nArrows = max(1, min(nArrows, nP));
idx = unique(round(linspace(1, nP, nArrows)));
end


%% ========================================================================
function [B, Uvec, Vvec, Nvec] = ArrowField(P, Nrm, U, V, idx, len, offset)
% Pure geometry: unit-direction glyphs scaled by len (scalar), bases offset off
% the surface along the vertex normal. No figure handles — testable in isolation.
unit = @(X) X ./ max(sqrt(sum(X.^2, 2)), eps);
n = unit(Nrm(idx, :));
B = P(idx, :) + offset .* n;
Uvec = unit(U(idx, :)) .* len;
Vvec = unit(V(idx, :)) .* len;
Nvec = n .* len;
end


%% ========================================================================
function L = MeanEdgeLength(Vtx, Fcs)
e1 = Vtx(Fcs(:,2),:) - Vtx(Fcs(:,1),:);
e2 = Vtx(Fcs(:,3),:) - Vtx(Fcs(:,2),:);
e3 = Vtx(Fcs(:,1),:) - Vtx(Fcs(:,3),:);
L  = mean(sqrt([sum(e1.^2,2); sum(e2.^2,2); sum(e3.^2,2)]));
end
```

- [ ] **Step 4: Run the headless test to verify it passes**

In the MATLAB MCP session: `rehash; test_manifold_frame`
Expected: `==== test_manifold_frame: 13 passed, 0 failed ====`

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/view_manifold.m dev/tests/test_manifold_frame.m
git commit -m "feat(manifold): view_manifold per-vertex frame viewer

Render a manifold node's canonical per-vertex tangent frame on the cortex:
U=real(grid.*rot), V=imag(grid.*rot), N=cross(U,V) from Embedded.vertex.grid +
Gauge.vertex.rotation, scattered per hemisphere via GlobalVertices, drawn as
headless quiver3 glyphs (view_tangents idiom) with a 'View' dispatch seam.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 2: Rewire the Manifold node menu to view_manifold

**Files:**
- Modify: `toolbox/tree/tree_callbacks.m` (manifold popup + callbacks)

**Context:** The manifold popup currently has "View" (`ManifoldView_Callback`, a field-dump stub) and the interim "Display tangent basis" (`ManifoldViewTangents_Callback` → `view_tangents`). Collapse to a single "View manifold" → `view_manifold`; remove both old callbacks.

- [ ] **Step 1: Replace the manifold popup items**

In `toolbox/tree/tree_callbacks.m`, find:

```matlab
                if (length(bstNodes) == 1)
                    % === VIEW ===
                    gui_component('MenuItem', jPopup, [], 'View', IconLoader.ICON_MATLAB, [], @(h,ev)bst_call(@ManifoldView_Callback, filenameFull));
                    % === DISPLAY TANGENT BASIS ===
                    gui_component('MenuItem', jPopup, [], 'Display tangent basis', IconLoader.ICON_DISPLAY, [], @(h,ev)bst_call(@ManifoldViewTangents_Callback, filenameFull));
                    % === DELETE ===
                    if ~bst_get('ReadOnly')
                        AddSeparator(jPopup);
                        gui_component('MenuItem', jPopup, [], 'Delete', IconLoader.ICON_DELETE, [], @(h,ev)bst_call(@ManifoldDelete_Callback, filenameRelative));
                    end
                end
```

Replace with:

```matlab
                if (length(bstNodes) == 1)
                    % === VIEW MANIFOLD ===
                    gui_component('MenuItem', jPopup, [], 'View manifold', IconLoader.ICON_DISPLAY, [], @(h,ev)bst_call(@view_manifold, filenameFull));
                    % === DELETE ===
                    if ~bst_get('ReadOnly')
                        AddSeparator(jPopup);
                        gui_component('MenuItem', jPopup, [], 'Delete', IconLoader.ICON_DELETE, [], @(h,ev)bst_call(@ManifoldDelete_Callback, filenameRelative));
                    end
                end
```

- [ ] **Step 2: Delete the two obsolete callbacks**

In `toolbox/tree/tree_callbacks.m`, delete the entire `ManifoldView_Callback` function (the field-dump inspector, header `%% ===== MANIFOLD: VIEW =====`) and the entire `ManifoldViewTangents_Callback` function (header `%% ===== MANIFOLD: DISPLAY TANGENT BASIS =====`). Leave `ManifoldDelete_Callback` and its `%% ===== MANIFOLD: DELETE =====` banner intact.

- [ ] **Step 3: Verify the rewire**

Run:
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
grep -n "ManifoldView_Callback\|ManifoldViewTangents_Callback\|Display tangent basis" toolbox/tree/tree_callbacks.m || echo "OBSOLETE GONE ✓"
grep -n "@view_manifold, filenameFull" toolbox/tree/tree_callbacks.m
```
Expected: first grep prints `OBSOLETE GONE ✓`; second prints the new "View manifold" wiring line.

- [ ] **Step 4: Parse-check (MATLAB MCP)**

```matlab
m = checkcode('/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox/tree/tree_callbacks.m','-id');
isErr = arrayfun(@(x) contains(lower(x.message),{'parse'}) || contains(lower(x.message),'unbalanced') || contains(lower(x.message),'unexpected'), m);
fprintf('parse-level errors: %d\n', nnz(isErr));
```
Expected: `parse-level errors: 0`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/tree/tree_callbacks.m
git commit -m "feat(manifold): single 'View manifold' menu item -> view_manifold

Collapse the Manifold node popup to one 'View manifold' action opening
view_manifold; remove the field-dump ManifoldView_Callback stub and the interim
ManifoldViewTangents_Callback (view_tangents) added during the relocation step.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 3: Repoint view_connection_phase to the manifold frame; deprecate tess_tangents

**Files:**
- Modify: `toolbox/gui/view_connection_phase.m`
- Modify: `toolbox/anatomy/tess_tangents.m`

**Context:** `view_connection_phase` builds its `FsFrame` from `tess_tangents` (per-face) → `bst_tangent_face2vertex` (per-vertex). `tess_frame(SurfaceFile)` returns the full-mesh per-vertex `[U,V,N]` directly from the manifold/facets bundle (computing + storing on the surface if absent), so it replaces both. `Nv = TessMat.VertNormals` is used only by `bst_tangent_face2vertex`, so it becomes unused and is removed. `tess_tangents` stays (still used by `bst_wavefront_track`, `tess_nxr_populate`) but gets a deprecation note.

- [ ] **Step 1: Repoint the frame source**

In `toolbox/gui/view_connection_phase.m`, find:

```matlab
Nv  = TessMat.VertNormals;

mctx   = nxr.manifold.context(Vtx, Fcs);
vFrame = nxr.manifold.measure.vertexFrame(mctx);
[Uf, ~]  = tess_tangents(SurfaceFile, 'NoSave', 1);
[Uv, Vv] = bst_tangent_face2vertex(Fcs, Uf, Nv);
FsFrame  = struct('e1', Uv, 'e2', Vv);
```

Replace with:

```matlab
mctx   = nxr.manifold.context(Vtx, Fcs);
vFrame = nxr.manifold.measure.vertexFrame(mctx);
[Uv, Vv] = tess_frame(SurfaceFile);          % manifold gauge frame (vertex domain)
FsFrame  = struct('e1', Uv, 'e2', Vv);
```

- [ ] **Step 2: Update the data-pipeline header comment**

In `toolbox/gui/view_connection_phase.m`, find:

```matlab
%     Data pipeline: bst_conn_eigenmodes_ensure -> nxr vertexFrame ->
%     bst_tangent_face2vertex (FS frame) -> bst_conn_phase.
```

Replace with:

```matlab
%     Data pipeline: bst_conn_eigenmodes_ensure -> nxr vertexFrame ->
%     tess_frame (manifold gauge frame) -> bst_conn_phase.
```

- [ ] **Step 3: Add the deprecation note to tess_tangents**

In `toolbox/anatomy/tess_tangents.m`, find:

```matlab
% TESS_TANGENTS: Globally consistent per-face tangent frame field.
%
```

Replace with:

```matlab
% TESS_TANGENTS: Globally consistent per-face tangent frame field.
%
% DEPRECATED: superseded by the manifold frame (view_manifold / tess_frame).
%     Retained only for the per-face callers (bst_wavefront_track,
%     tess_nxr_populate) until per-face manifold frames are available
%     (nxr Gauge.face.rotation is currently deferred for the trivial gauge).
%
```

- [ ] **Step 4: Verify the repoint**

Run:
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
grep -n "tess_tangents\|bst_tangent_face2vertex\|VertNormals" toolbox/gui/view_connection_phase.m || echo "REPOINTED ✓"
grep -n "DEPRECATED" toolbox/anatomy/tess_tangents.m
```
Expected: first grep prints `REPOINTED ✓` (no `tess_tangents`/`face2vertex`/`VertNormals` left); second shows the deprecation line.

- [ ] **Step 5: Sanity-check tess_frame returns a vertex frame (MATLAB MCP)**

```matlab
S = bst_canonical_cortex(20484);
[U,V] = tess_frame(S);
fprintf('tess_frame U size [%s], unit=%d\n', num2str(size(U)), all(abs(sqrt(sum(U.^2,2))-1)<1e-6));
```
Expected: `U size [20484 3], unit=1`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/view_connection_phase.m toolbox/anatomy/tess_tangents.m
git commit -m "refactor(manifold): repoint view_connection_phase frame to the manifold

Source the per-vertex reference frame from tess_frame (manifold gauge frame)
instead of tess_tangents + bst_tangent_face2vertex. Mark tess_tangents
@deprecated (kept for the per-face callers until nxr fills Gauge.face.rotation).

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Task 4: Delete view_tangents; live viewer test

**Files:**
- Delete: `toolbox/gui/view_tangents.m`
- Delete: `dev/tests/test_view_tangents.m`
- Create: `dev/tests/test_view_manifold.m`

**Context:** With the menu and `view_connection_phase` repointed, the only remaining `view_tangents` caller is gone, so the viewer and its test are deleted. A live test exercises `view_manifold` against the registered manifold node.

- [ ] **Step 1: Confirm no remaining view_tangents callers**

Run:
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
grep -rn "view_tangents" toolbox/ | grep -v "view_tangents.m:"
```
Expected: empty (no callers in `toolbox/` outside the file itself). If any remain, STOP — they must be repointed first.

- [ ] **Step 2: Delete the viewer and its test**

```bash
git rm toolbox/gui/view_tangents.m dev/tests/test_view_tangents.m
```

- [ ] **Step 3: Write the live integration test**

Create `dev/tests/test_view_manifold.m`:

```matlab
function test_view_manifold()
% TEST_VIEW_MANIFOLD  Live-figure regression for the manifold frame viewer.
% Requires Brainstorm running with a registered manifold node.
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;

    % --- find a registered manifold node ---
    ManifoldFile = local_find_manifold();
    assert(~isempty(ManifoldFile), 'No registered manifold node found; compute a manifold first.');

    close(findobj(0, 'type', 'figure', 'Tag', '3DViz'));
    hFig = view_manifold(ManifoldFile);
    drawnow;
    [nPass,nFail] = chk('viewer returns a figure', ~isempty(hFig) && ishandle(hFig), nPass,nFail);

    hAx3D  = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    hPatch = findobj(hAx3D, 'Type', 'patch');
    hU = findobj(hFig, 'Tag', 'tangentU');
    hV = findobj(hFig, 'Tag', 'tangentV');
    [nPass,nFail] = chk('Axes3D survives glyph draw', ~isempty(hAx3D), nPass,nFail);
    [nPass,nFail] = chk('cortex patch survives glyph draw', ~isempty(hPatch), nPass,nFail);
    [nPass,nFail] = chk('U glyphs drawn', ~isempty(hU) && numel(get(hU(1),'UData'))>0, nPass,nFail);
    [nPass,nFail] = chk('V glyphs drawn', ~isempty(hV), nPass,nFail);

    % N key adds the normal quiver
    KeyOnFig(hFig, 'n'); drawnow;
    [nPass,nFail] = chk('N toggles normal glyphs', ~isempty(findobj(hFig,'Tag','tangentN')), nPass,nFail);

    % density key changes glyph count
    n1 = numel(get(findobj(hFig,'Tag','tangentU'),'UData'));
    KeyOnFig(hFig, 'leftarrow'); drawnow;
    n2 = numel(get(findobj(hFig,'Tag','tangentU'),'UData'));
    [nPass,nFail] = chk('density key changes glyph count', n2 < n1, nPass,nFail);

    close(hFig);
    fprintf('\n==== test_view_manifold: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_view_manifold: %d test(s) FAILED.', nFail); end
end

function f = local_find_manifold()
    f = '';
    P = bst_get('ProtocolSubjects'); subs = [P.Subject];
    for s = 1:numel(subs)
        for k = 1:numel(subs(s).Surface)
            if isfield(subs(s).Surface(k),'Manifold') && ~isempty(subs(s).Surface(k).Manifold)
                f = subs(s).Surface(k).Manifold(1).FileName; return;
            end
        end
    end
end

function KeyOnFig(hFig, keyName)
    ev.Key = keyName; ev.Character = ''; ev.Modifier = {};
    cb = get(hFig, 'KeyPressFcn'); cb(hFig, ev);
end

function [p,f] = chk(name, c, p, f)
    if c, fprintf('  PASS  %s\n', name); p=p+1; else, fprintf('  FAIL  %s\n', name); f=f+1; end
end
```

- [ ] **Step 4: Run the live test (MATLAB MCP)**

In the MATLAB MCP session (Brainstorm running, manifold node present): `rehash; test_view_manifold`
Expected: `==== test_view_manifold: 7 passed, 0 failed ====`

- [ ] **Step 5: Re-run the headless test (regression guard)**

`test_manifold_frame`
Expected: `13 passed, 0 failed`.

- [ ] **Step 6: Commit**

```bash
git add dev/tests/test_view_manifold.m
git commit -m "test(manifold): delete view_tangents; live view_manifold test

Remove the superseded view_tangents.m + its test (no callers remain), and add a
live regression that opens view_manifold on the registered manifold node and
asserts the cortex survives the glyph draw, U/V frame glyphs render, and the N /
density keys work.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Self-Review notes

- **Spec coverage:** §4.1 view_manifold + `DeriveVertexFrame` → Task 1; §4.2 menu rewire → Task 2; §4.3 view_connection_phase repoint → Task 3; §4.4 tess_tangents deprecation → Task 3; §4.5 deletions → Task 4; §7.1 headless → Task 1; §7.2 live → Task 4; §7.3 repoint guard → Task 3 (grep) + Task 4 (no-callers grep).
- **Type/name consistency:** `DeriveVertexFrame(Embedded, Gauge, nVert)` returns `struct(P,U,V,N,Sing)` — identical in the viewer, the macro dispatch, and the headless test. Glyph tags `tangentU/tangentV/tangentN/tangentSing` consistent between viewer and live test. `view_manifold(ManifoldFile)` single-full-path contract used by the menu (Task 2) and both tests.
- **Quaternion/frame convention** matches `tess_frame`: `U=real(grid.*rot)`, `V=imag(grid.*rot)`, `N=cross(U,V)`.
- **Axes-reset trap** handled (`hold(hAxes,'on')` before `quiver3`).
