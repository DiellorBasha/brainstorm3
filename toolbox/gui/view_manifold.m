function varargout = view_manifold(varargin)
% VIEW_MANIFOLD: Layered viewer for a manifold node's data dimensions.
%
% USAGE:  hFig = view_manifold(ManifoldFile)
%         G    = view_manifold('DeriveVertexFrame', Embedded, Gauge, nVert)
%         G    = view_manifold('DeriveFaceFrame',   Embedded, Gauge, nFace)
%
% Reads a manifold_ DB node (db_template('manifoldmat')) and renders its parent
% cortex as a true light-grey wireframe (edges only, no filled faces), onto which
% the data dimensions of the manifold are revealed as keyboard-toggled layers:
%     scalar  (vertices)        -> 'D' : a point cloud on every vertex
%     vector2 (tangent frames)  -> the combed tangent frame (U,V)
%     vector3 (ambient frames)  -> the combed full frame (U,V,N)
% Each dimension's frame is the geometric object the data lives on: vertices for
% scalar, tangent (U,V) frames for vector2, ambient (U,V,N) frames for vector3.
%
% The Support buttons switch every layer between vertex support and face support.
% On vertices the combed field is real(Gauge.vertex.rotation .* Embedded.vertex.grid)
% placed at vertex positions; on faces it is real(Gauge.face.rotation .*
% Embedded.face.grid) placed at face centroids (the manifold's dual mesh).
%
% Frame glyphs are unit (orthonormal), so they carry no magnitude: each is sized
% to its own mesh cell using the refinement-covariant local length sqrt(cell area)
% (barycentric dual area on vertices, triangle area on faces). Glyphs therefore
% scale with the tessellation and the figure reads the same at any mesh density.
%
% The scalar layer is a point cloud kept in sync with the cortex patch via its
% MarkedClean event, so it follows the Surfaces-panel Smooth slider and numeric
% Resect slider (vertex motion) and hides the hemisphere hidden by a left/right/
% struct Resect (read from the patch's per-face alpha). The default view is the
% wireframe with the scalar layer on. The pure per-vertex frame derivation
% (DeriveVertexFrame) is retained for the coming vector layers and for
% view_manifold_registration.
%
% Keyboard (figure focused):
%   D   toggle the scalar data layer (vertex point cloud)
%   H   help
%
% SEE ALSO: tess_manifold, view_manifold_registration, view_surface
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

if (nargin >= 1) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'DeriveVertexFrame','DeriveFaceFrame'}))
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: per-vertex frame from the manifold node =====
function G = DeriveVertexFrame(Embedded, Gauge, nVert)
% G: struct with P,U,V,N [nVert x 3] (zeros off-support) and Sing (global vtx ids).
% U=real(grid.*rot), V=imag(grid.*rot), N=cross(U,V) per the manifold gauge convention.
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


%% ===== PURE: per-face (dual-mesh) frame from the manifold node =====
function G = DeriveFaceFrame(Embedded, Gauge, nFace)
% G: struct with P,U,V,N [nFace x 3] (zeros off-support) and empty Sing.
% P = face centroid; U=real(grid.*rot), V=imag(grid.*rot), N=cross(U,V) on faces,
% using the per-face combed field (Embedded.face.grid, Gauge.face.rotation).
    if numel(Embedded) ~= 2 || numel(Gauge) ~= 2
        error('view_manifold:badNode', 'Embedded and Gauge must be 1x2 per-hemisphere structs.');
    end
    P = zeros(nFace,3); U = zeros(nFace,3); V = zeros(nFace,3); N = zeros(nFace,3);
    for hh = 1:2
        fH   = double(Embedded(hh).GlobalFaces(:));
        grid = Embedded(hh).face.grid;
        rot  = Gauge(hh).face.rotation;
        if isempty(rot)
            error('view_manifold:noFaceGauge', ...
                ['Hemisphere %d: Gauge.face.rotation is empty. Recompute the manifold ' ...
                 'with the current nxr-compute (tess_manifold ''ForceRecompute'').'], hh);
        end
        if size(grid,1) ~= numel(fH) || size(grid,2) ~= 3
            error('view_manifold:shapeMismatch', ...
                'Hemisphere %d: Embedded.face.grid is [%s], expected [%d x 3].', ...
                hh, num2str(size(grid)), numel(fH));
        end
        cRot = grid .* rot(:);
        Uh = real(cRot);  Vh = imag(cRot);
        U(fH,:) = Uh;
        V(fH,:) = Vh;
        N(fH,:) = cross(Uh, Vh, 2);
        P(fH,:) = Embedded(hh).face.centroid;
    end
    G = struct('P',P, 'U',U, 'V',V, 'N',N, 'Sing',[]);
end


%% ===== GUI: layered manifold viewer (wireframe + toggleable layers) =====
function hFig = ViewFigure(ManifoldFile)
    hFig = [];
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
    Surface = M.ParentSurface;

    % --- open surface figure as a true wireframe (edges only, no filled faces) ---
    hFig = view_surface(Surface, 0, [.5 .5 .5], 'NewFigure', 0);
    if isempty(hFig)
        bst_error('Could not open the surface figure.', 'View manifold', 0);
        return;
    end
    set(hFig, 'Name', ['Manifold: ' Surface]);
    panel_surface('SetSurfaceSmooth', hFig, 1, 0, 0);
    panel_surface('SetSurfaceEdges',  hFig, 1, 1);
    hAxes  = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    hPatch = findobj(hAxes, 'Type', 'patch');
    if numel(hPatch) > 1, hPatch = hPatch(1); end
    % Opaque surface render (faces on) with the mesh edges; overlays use a
    % separate cloud, so the patch itself carries no markers.
    set(hPatch, 'Marker', 'none');
    hold(hAxes, 'on');

    % --- state ---
    activeDim  = 'scalar';      % which data dimension is shown (tab-driven)
    support    = 'vertex';      % data support: 'vertex' or 'face' (all layers)
    showScalar = true;          % within the scalar dim, D toggles the cloud
    colScalar  = [0.2 0.9 1];   % cyan point cloud
    colTangent = [1 1 0];       % yellow U,V tangent frame
    colNormal  = [1 0 1];       % magenta normal
    Gframe     = [];            % cached per-vertex frame (DeriveVertexFrame)
    GframeF    = [];            % cached per-face frame   (DeriveFaceFrame)
    nFrames    = 2500;          % frame glyphs are decimated for readability
    meanEdge   = MeanEdgeLength(get(hPatch,'Vertices'), get(hPatch,'Faces'));
    offLen     = 0.2 * meanEdge;   % scalar-cloud lift off the opaque surface
    % Frame glyphs are UNIT (orthonormal) — they carry no magnitude, so their size
    % is purely legibility: tie each cross to its own cell via the refinement-
    % covariant local length sqrt(cell area). kGlyph sets the arm length, liftFac
    % the lift off the surface; both are fractions of that local length.
    kGlyph     = 0.45;          % frame arm length = kGlyph  * sqrt(local cell area)
    liftFac    = 0.40;          % frame lift along N = liftFac * sqrt(local cell area)
    Vcache = []; Acache = []; EAcache = []; inSync = false;

    % Scalar layer: a point cloud that mirrors the cortex patch's live vertices.
    V0 = get(hPatch, 'Vertices');
    hCloud = plot3(V0(:,1), V0(:,2), V0(:,3), '.', 'Parent', hAxes, ...
        'Color', colScalar, 'MarkerSize', 6, 'LineStyle', 'none', 'Tag', 'manifoldScalar');

    KeyPressFcn_bak = get(hFig, 'KeyPressFcn');
    set(hFig, 'KeyPressFcn', @KeyPress_Callback);
    hLabel = uicontrol('Style','text','String','...','Units','Pixels', ...
        'Position',[6 4 1000 18],'HorizontalAlignment','left', ...
        'FontUnits','points','FontSize',bst_get('FigFont'), ...
        'ForegroundColor',[1 1 1],'BackgroundColor',[0 0 0],'Parent',hFig);

    % Re-sync the cloud whenever the cortex patch is redrawn (Smooth / Resect).
    lh = addlistener(hPatch, 'MarkedClean', @(s,e) SyncScalar(false));
    setappdata(hFig, 'ManifoldScalarListener', lh);
    SyncScalar(true);
    UpdateLabel();

    % --- dimension tabs (top-right): switch the active data-dimension view ---
    hDimTabs = uitabgroup('Parent', hFig, 'Units', 'normalized', ...
        'Position', [0.595 0.90 0.40 0.095], 'Tag', 'manifoldDimTabs', ...
        'SelectionChangedFcn', @DimTabChanged);
    uitab(hDimTabs, 'Title', 'Scalar');
    uitab(hDimTabs, 'Title', 'Vector2');
    uitab(hDimTabs, 'Title', 'Vector3');

    % --- support buttons (top-left): vertex vs face support (mutually exclusive) ---
    hSupport = uibuttongroup('Parent', hFig, 'Units', 'normalized', ...
        'Position', [0.005 0.90 0.165 0.095], 'Tag', 'manifoldSupport', ...
        'Title', 'Support', 'ForegroundColor', [1 1 1], 'BackgroundColor', [0.15 0.15 0.15], ...
        'HighlightColor', [.4 .4 .4], 'SelectionChangedFcn', @SupportChanged);
    uicontrol(hSupport, 'Style', 'togglebutton', 'String', 'Vertex', 'Units', 'normalized', ...
        'Position', [0.04 0.12 0.45 0.76], 'Value', 1);
    uicontrol(hSupport, 'Style', 'togglebutton', 'String', 'Face', 'Units', 'normalized', ...
        'Position', [0.51 0.12 0.45 0.76]);

    % ===== NESTED: keep the scalar cloud in sync with the cortex patch =====
    function SyncScalar(force)
        if inSync || isempty(hPatch) || ~ishandle(hPatch) || isempty(hCloud) || ~ishandle(hCloud)
            return;
        end
        if ~(strcmpi(activeDim,'scalar') && showScalar)
            set(hCloud, 'Visible', 'off');
            Vcache = []; Acache = []; EAcache = [];   % force a re-sync when turned back on
            return;
        end
        V  = get(hPatch, 'Vertices');
        A  = get(hPatch, 'FaceVertexAlphaData');
        EA = get(hPatch, 'EdgeAlpha');     % 'flat'/'interp' (per-element) or scalar
        if ~force && isequal(V,Vcache) && isequal(A,Acache) && isequal(EA,EAcache) ...
                && strcmpi(get(hCloud,'Visible'),'on')
            return;
        end
        inSync = true;
        Vcache = V; Acache = A; EAcache = EA;
        F  = get(hPatch, 'Faces');
        nV = size(V, 1);
        nF = size(F, 1);
        % Per-face visibility from the alpha MODE. Brainstorm hides a hemisphere
        % by switching the edge alpha to per-element ('flat'/'interp') with zeroed
        % FaceVertexAlphaData; toggling off restores a SCALAR alpha but leaves the
        % per-face data stale. So gate on the mode, not the (stale) per-face data.
        if ischar(EA) && ~isempty(A)
            if numel(A) == nF               % per-face ('flat')
                faceVis = A > 0;
            elseif numel(A) == nV           % per-vertex ('interp')
                vv = A(:) > 0; faceVis = all(vv(F), 2);
            else
                faceVis = true(nF, 1);
            end
        elseif isnumeric(EA) && isscalar(EA) && EA <= 0
            faceVis = false(nF, 1);         % uniformly transparent
        else
            faceVis = true(nF, 1);          % scalar opaque -> all visible
        end
        % Support: scalar data on vertices (default) or on faces (centroids).
        if strcmpi(support, 'face')
            P   = (V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:)) / 3;   % live face centroids
            vis = faceVis;
        else
            P   = V;                                              % vertices
            vis = false(nV, 1);
            vis(unique(F(faceVis, :))) = true;                    % visible if any incident face is
        end
        % Lift the cloud just off the opaque surface, radially outward from the
        % live centroid (always outward; follows Smooth / Resect motion).
        ctr = mean(V, 1);
        rad = P - ctr;
        rad = rad ./ max(sqrt(sum(rad.^2, 2)), eps);
        Pd = P + offLen * rad;
        Pd(~vis, :) = NaN;                  % hidden support points drop out
        set(hCloud, 'XData', Pd(:,1), 'YData', Pd(:,2), 'ZData', Pd(:,3), 'Visible', 'on');
        inSync = false;
    end

    % ===== NESTED: apply the active dimension (scalar / vector2 / vector3) =====
    function ApplyDim()
        SyncScalar(true);   % scalar cloud shown only when activeDim == 'scalar'
        DrawVectors();      % draw / clear the per-vertex frame glyphs
        UpdateLabel();
    end

    % ===== NESTED: combed-frame glyphs (vector2 = U,V; vector3 = U,V,N) =====
    % Draws on the active support: vertices (vertex.grid/rotation at vertex
    % positions) or face centroids (face.grid/rotation on the dual mesh).
    function DrawVectors()
        delete(findobj(hAxes, '-depth', 1, '-regexp', 'Tag', '^manifoldVec'));
        if ~any(strcmpi(activeDim, {'vector2','vector3'}))
            return;
        end
        if strcmpi(support, 'face')
            if isempty(GframeF)
                GframeF = DeriveFaceFrame(M.Embedded, M.Gauge, size(get(hPatch,'Faces'),1));
            end
            Gf = GframeF;
        else
            if isempty(Gframe)
                Gframe = DeriveVertexFrame(M.Embedded, M.Gauge, size(get(hPatch,'Vertices'),1));
            end
            Gf = Gframe;
        end
        nP  = size(Gf.P, 1);
        idx = unique(round(linspace(1, nP, min(nFrames, nP))));   % decimate for readability
        % Principled per-cell sizing: arm length = kGlyph * sqrt(local cell area),
        % a refinement-covariant world length (vertex: barycentric dual area; face:
        % triangle area). Glyphs then shrink in proportion to the mesh, so the figure
        % reads the same at any tessellation density.
        L  = LengthScale(support);                            % [nElem x 1] world length
        Lk = L(idx);
        P  = Gf.P(idx, :) + (liftFac * Lk) .* Gf.N(idx, :);   % per-cell lift off surface
        DrawQuiv(P, (kGlyph * Lk) .* Gf.U(idx,:), colTangent, 'manifoldVecU');
        DrawQuiv(P, (kGlyph * Lk) .* Gf.V(idx,:), colTangent, 'manifoldVecV');
        if strcmpi(activeDim, 'vector3')
            DrawQuiv(P, (kGlyph * Lk) .* Gf.N(idx,:), colNormal, 'manifoldVecN');
        end
    end

    % Refinement-covariant local length scale = sqrt(cell area), computed from the
    % live patch (vertex: barycentric lumped dual area, identical to nxr's
    % Intrinsic.vertex.dualArea; face: triangle area).
    function L = LengthScale(supp)
        Vv = get(hPatch, 'Vertices');  Ff = double(get(hPatch, 'Faces'));
        e1 = Vv(Ff(:,2),:) - Vv(Ff(:,1),:);
        e2 = Vv(Ff(:,3),:) - Vv(Ff(:,1),:);
        Af = 0.5 * sqrt(sum(cross(e1, e2, 2).^2, 2));        % triangle areas
        if strcmpi(supp, 'face')
            L = sqrt(Af);
        else
            L = sqrt(accumarray(Ff(:), repmat(Af,3,1)/3, [size(Vv,1) 1]));
        end
    end

    function DrawQuiv(P, A, col, tag)
        % A holds the already-scaled arm vectors (world units). Pass scale arg 0 so
        % quiver3 uses them literally instead of auto-rescaling to "look nice".
        quiver3(P(:,1), P(:,2), P(:,3), A(:,1), A(:,2), A(:,3), 0, ...
            'Parent', hAxes, 'Color', col, 'LineWidth', 1, 'ShowArrowHead', 'off', 'Tag', tag);
    end

    % ===== NESTED: status label =====
    function UpdateLabel()
        switch lower(activeDim)
            case 'vector2'
                set(hLabel, 'String', sprintf('Manifold  |  Vector2 (%s): tangent frame (U,V)   (H: help)', support));
            case 'vector3'
                set(hLabel, 'String', sprintf('Manifold  |  Vector3 (%s): full frame (U,V,N)   (H: help)', support));
            otherwise
                if showScalar, st = 'on'; else, st = 'off'; end
                set(hLabel, 'String', sprintf('Manifold  |  Scalar (%s): %s   (D: toggle   H: help)', support, st));
        end
    end

    % ===== NESTED: dimension tab switch (Scalar / Vector2 / Vector3) =====
    function DimTabChanged(~, ev)
        activeDim = lower(ev.NewValue.Title);
        ApplyDim();
    end

    % ===== NESTED: support switch (vertex / face) — drives every layer =====
    function SupportChanged(~, ev)
        support = lower(ev.NewValue.String);
        ApplyDim();   % re-sync the scalar cloud and redraw the active frame layer
    end

    % ===== NESTED: keyboard =====
    function KeyPress_Callback(h, ev)
        switch ev.Key
            case 'd'
                if strcmpi(activeDim, 'scalar')
                    showScalar = ~showScalar;
                    SyncScalar(true);
                    UpdateLabel();
                end
            case 'h'
                java_dialog('msgbox', ['<HTML><TABLE>' ...
                    '<TR><TD><B>D</B></TD><TD>Toggle the scalar data layer (vertex point cloud)</TD></TR>' ...
                    '<TR><TD><B>0-9</B></TD><TD>Change view</TD></TR>' ...
                    '</TABLE>'], 'Manifold viewer shortcuts', [], 0);
            otherwise
                if ~isempty(KeyPressFcn_bak), KeyPressFcn_bak(h, ev); end
        end
    end
end


%% ========================================================================
function L = MeanEdgeLength(Vtx, Fcs)
Fcs = double(Fcs);
e1 = Vtx(Fcs(:,2),:) - Vtx(Fcs(:,1),:);
e2 = Vtx(Fcs(:,3),:) - Vtx(Fcs(:,2),:);
e3 = Vtx(Fcs(:,1),:) - Vtx(Fcs(:,3),:);
L  = mean(sqrt([sum(e1.^2,2); sum(e2.^2,2); sum(e3.^2,2)]));
end
