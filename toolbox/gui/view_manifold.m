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
