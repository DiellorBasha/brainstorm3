function hFig = view_tangents(SurfaceFile, varargin)
% VIEW_TANGENTS: Display the per-face tangent frame field stored by tess_tangents.
%
% USAGE:  hFig = view_tangents(SurfaceFile)
%         hFig = view_tangents(SurfaceFile, 'MaxArrows', 3000)
%
% DESCRIPTION:
%     Overlays the per-face tangent frame field (TessMat.TangentFrame, computed
%     by tess_tangents) on a Brainstorm surface figure. Each frame is drawn as
%     headless lines (axis-like): U and V (the in-plane tangent cross) share one
%     color, the face normal another. Glyphs are sized per-face to the triangle's
%     inscribed circle so they sit flat within each face. The registration-pole
%     singularities are marked. The cortex is shown opaque, unsmoothed, with the
%     wireframe edges on. If the surface has no stored frame, it is computed and
%     stored via tess_tangents first.
%
%     Keyboard (figure focused):
%       Left/Right          fewer / more frames
%       Shift + Up/Down     glyph length
%       Ctrl  + Up/Down     line width
%       N                   toggle face-normal glyphs (off by default)
%       P                   toggle singularity markers
%       H                   help
%
% INPUT:
%     - SurfaceFile : Brainstorm cortex surface (relative or full path).
% OPTIONS:
%     - MaxArrows : initial number of face frames drawn (default: half the faces).
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
MaxArrows = [];   % default resolved below to half the faces
if mod(numel(varargin), 2) ~= 0
    error('view_tangents:badArgs', 'Optional arguments must be name/value pairs.');
end
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
% Per-face edge lengths and inscribed-circle radius. The radius is the per-face
% glyph length, so each frame fits flat inside its own triangle.
e1 = sqrt(sum((Vtx(Fcs(:,2),:) - Vtx(Fcs(:,1),:)).^2, 2));
e2 = sqrt(sum((Vtx(Fcs(:,3),:) - Vtx(Fcs(:,2),:)).^2, 2));
e3 = sqrt(sum((Vtx(Fcs(:,1),:) - Vtx(Fcs(:,3),:)).^2, 2));
meanEdge = mean([e1; e2; e3]);
sP       = (e1 + e2 + e3) / 2;
faceArea = sqrt(max(sP .* (sP - e1) .* (sP - e2) .* (sP - e3), 0));
inRadius = faceArea ./ max(sP, eps);          % [nF x 1] inscribed-circle radius
% Singularity poles drawn as "lollipops": a marker lifted RADIALLY outward from
% the surface centroid (so it clears sulcal folds — a face-normal offset fails
% for poles on steep walls, which are nearly tangential), with a thin stem back
% to the true pole vertex. Keeps the cortex fully opaque.
sv       = TF.Singularities.Vertices(:);
SingBase = Vtx(sv, :);
ctrAll   = mean(Vtx, 1);
radial   = SingBase - ctrAll;
radial   = radial ./ max(sqrt(sum(radial.^2, 2)), eps);
liftLen  = 0.05 * max(max(Vtx,[],1) - min(Vtx,[],1));
SingTip  = SingBase + liftLen .* radial;

%% ===== OPEN SURFACE FIGURE =====
% Opaque (SurfAlpha=0) so the frame reads clearly against the cortex.
% Disable scout overlay (isScouts=0) — scout panel cannot handle surfaces
% that are absolute paths to files not registered in the Brainstorm DB.
hFig = view_surface(SurfaceFile, 0, [.5 .5 .5], 'NewFigure', 0);
if isempty(hFig)
    error('view_tangents:noFigure', 'Could not open the surface figure.');
end
set(hFig, 'Name', ['Tangent basis: ' SurfaceFile]);
% Default surface display: opaque, no smoothing, wireframe edges ON.
panel_surface('SetSurfaceSmooth', hFig, 1, 0, 0);
panel_surface('SetSurfaceEdges',  hFig, 1, 1);
hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
hold(hAxes, 'on');

%% ===== DISPLAY STATE =====
if isempty(MaxArrows)
    nArrows = round(nF / 2);   % default: frames on half the faces
else
    nArrows = min(MaxArrows, nF);
end
quiverSize  = 1;
quiverWidth = 1;
showNormals = false;     % normals off by default (toggle with N)
showSing    = true;
colTangent  = [1 1 0];      % yellow : U and V
colNormal   = [1 0 1];      % magenta: face normal
colSing     = [0 0.45 1];   % blue   : singularity poles

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
        % Per-face glyph length = inscribed-circle radius (fits flat in the
        % triangle); bases lifted slightly off the surface. autoscale OFF (0).
        [B, Uvec, Vvec, Nvec] = ArrowField(Centroids, FaceNormals, U, V, idx, ...
            quiverSize * inRadius(idx), 0.1 * meanEdge);
        % Frame glyphs as headless lines (axis-like); U and V share a color.
        hUq = quiver3(B(:,1), B(:,2), B(:,3), Uvec(:,1), Uvec(:,2), Uvec(:,3), 0, ...
            'Parent', hAxes, 'Color', colTangent, 'LineWidth', quiverWidth, ...
            'ShowArrowHead', 'off', 'Tag', 'tangentU');
        quiver3(B(:,1), B(:,2), B(:,3), Vvec(:,1), Vvec(:,2), Vvec(:,3), 0, ...
            'Parent', hAxes, 'Color', colTangent, 'LineWidth', quiverWidth, ...
            'ShowArrowHead', 'off', 'Tag', 'tangentV');
        % Face normal (different color), also headless.
        hNq = [];
        if showNormals
            hNq = quiver3(B(:,1), B(:,2), B(:,3), Nvec(:,1), Nvec(:,2), Nvec(:,3), 0, ...
                'Parent', hAxes, 'Color', colNormal, 'LineWidth', quiverWidth, ...
                'ShowArrowHead', 'off', 'Tag', 'tangentN');
        end
        % Singularity poles as lollipops: a stem from the true vertex to a marker
        % lifted radially outward so it clears the opaque cortex.
        hSq = [];
        if showSing && ~isempty(SingTip)
            line([SingBase(:,1)'; SingTip(:,1)'], ...
                 [SingBase(:,2)'; SingTip(:,2)'], ...
                 [SingBase(:,3)'; SingTip(:,3)'], ...
                 'Parent', hAxes, 'Color', colSing, 'LineWidth', 1.5, 'Tag', 'tangentSingStem');
            hSq = plot3(SingTip(:,1), SingTip(:,2), SingTip(:,3), 'o', ...
                'Parent', hAxes, 'MarkerFaceColor', colSing, 'MarkerEdgeColor', [.2 .2 .2], ...
                'MarkerSize', 10, 'LineStyle', 'none', 'Tag', 'tangentSing');
        end
        % Legend (meaning of the glyphs).
        legH = hUq;  legL = {'Tangent frame (U,V)'};
        if ~isempty(hNq),  legH = [legH, hNq];  legL{end+1} = 'Face normal';            end
        if ~isempty(hSq),  legH = [legH, hSq];  legL{end+1} = 'FreeSurfer sphere pole';  end
        legend(legH, legL, 'TextColor', [1 1 1], 'Color', [0 0 0], ...
            'Location', 'NorthEast', 'Interpreter', 'none', 'Tag', 'tangentLegend');
        % Status label.
        set(hLabel, 'String', sprintf( ...
            '%d / %d faces  |  N: normals   P: poles   <-/->: density   H: help', numel(idx), nF));
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
                    quiverSize = max(0.05, quiverSize / 1.2);
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
                    '<TR><TD><B>Left / Right</B></TD><TD>Fewer / more frames</TD></TR>' ...
                    '<TR><TD><B>Shift + Up/Down</B></TD><TD>Glyph length</TD></TR>' ...
                    '<TR><TD><B>Ctrl + Up/Down</B></TD><TD>Line width</TD></TR>' ...
                    '<TR><TD><B>N</B></TD><TD>Toggle face-normal glyphs</TD></TR>' ...
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
% Pure geometry: unit-direction glyphs scaled by len (scalar or per-face column),
% bases offset off the surface along the face normal. No figure handles —
% testable in isolation.
unit = @(X) X ./ max(sqrt(sum(X.^2, 2)), eps);
n = unit(FaceNormals(idx, :));
B = Centroids(idx, :) + offset .* n;
Uvec = unit(U(idx, :)) .* len;
Vvec = unit(V(idx, :)) .* len;
Nvec = n .* len;
end
