function varargout = view_manifold(varargin)
% VIEW_MANIFOLD: Layered viewer for a manifold node's data dimensions.
%
% USAGE:  hFig = view_manifold(ManifoldFile)
%         G    = view_manifold('DeriveVertexFrame', Embedded, Gauge, nVert)
%
% Reads a manifold_ DB node (db_template('manifoldmat')) and renders its parent
% cortex as a true light-grey wireframe (edges only, no filled faces), onto which
% the data dimensions of the manifold are revealed as keyboard-toggled layers:
%     scalar  (vertices)        -> 'D' : a point cloud on every vertex
%     vector2 (tangent frames)  -> (added later)
%     vector3 (ambient frames)  -> (added later)
% Each dimension's frame is the geometric object the data lives on: vertices for
% scalar, tangent (U,V) frames for vector2, ambient (U,V,N) frames for vector3.
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
% SEE ALSO: tess_manifold, tess_frame, view_manifold_registration, view_surface
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
% U=real(grid.*rot), V=imag(grid.*rot), N=cross(U,V) per the tess_frame convention.
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
    % Light-grey wireframe that reads against the black background.
    set(hPatch, 'FaceColor', 'none', 'EdgeColor', [0.7 0.7 0.7], 'EdgeAlpha', 1, 'Marker', 'none');
    hold(hAxes, 'on');

    % --- state ---
    showScalar = true;          % scalar layer on by default
    colScalar  = [0.2 0.9 1];   % cyan point cloud
    Vcache = []; Acache = []; inSync = false;

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

    % ===== NESTED: keep the scalar cloud in sync with the cortex patch =====
    function SyncScalar(force)
        if inSync || isempty(hPatch) || ~ishandle(hPatch) || isempty(hCloud) || ~ishandle(hCloud)
            return;
        end
        if ~showScalar
            set(hCloud, 'Visible', 'off');
            Vcache = []; Acache = [];   % force a re-sync when turned back on
            return;
        end
        V = get(hPatch, 'Vertices');
        A = get(hPatch, 'FaceVertexAlphaData');
        if ~force && isequal(V, Vcache) && isequal(A, Acache) && strcmpi(get(hCloud,'Visible'),'on')
            return;
        end
        inSync = true;
        Vcache = V; Acache = A;
        F  = get(hPatch, 'Faces');
        nV = size(V, 1);
        % A vertex is visible if at least one of its incident faces is visible
        % (per-face alpha > 0). With no per-face alpha, everything is visible.
        if isempty(A) || numel(A) ~= size(F, 1)
            vis = true(nV, 1);
        else
            vis = false(nV, 1);
            vis(unique(F(A > 0, :))) = true;
        end
        Vd = V; Vd(~vis, :) = NaN;     % hidden-hemisphere points drop out
        set(hCloud, 'XData', Vd(:,1), 'YData', Vd(:,2), 'ZData', Vd(:,3), 'Visible', 'on');
        inSync = false;
    end

    % ===== NESTED: status label =====
    function UpdateLabel()
        if showScalar, st = 'on'; else, st = 'off'; end
        set(hLabel, 'String', sprintf('Manifold  |  scalar layer: %s   (D: scalar   H: help)', st));
    end

    % ===== NESTED: keyboard =====
    function KeyPress_Callback(h, ev)
        switch ev.Key
            case 'd'
                showScalar = ~showScalar;
                SyncScalar(true);
                UpdateLabel();
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
