function varargout = view_manifold_registration(varargin)
% VIEW_MANIFOLD_REGISTRATION: Display a manifold node's poles on the registration sphere.
%
% USAGE:  hFig = view_manifold_registration(ManifoldFile)
%         [Base,Tip] = view_manifold_registration('RegSingGlyphs', sphV, poleIdx, ir, il)
%
% Opens Brainstorm's standard FreeSurfer registration-sphere figure for the
% manifold's parent surface (via view_surface_sphere) and overlays the manifold's
% singularity poles as lollipop pins -- the same pins view_manifold draws on the
% cortex, here shown on the sphere. First step toward validating analytic
% functions (heat, wave) against their closed-form counterparts on the sphere.
%
% Poles only; tangent frames and analytic colormaps extend the same overlay seam.
%
% Keyboard (figure focused):
%   P   toggle singularity pins
%   H   help
%
% SEE ALSO: view_manifold, view_surface_sphere, tess_manifold
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

if (nargin >= 1) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'RegSingGlyphs'}))
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: pole lollipop geometry on the registration sphere =====
function [Base, Tip] = RegSingGlyphs(sphV, poleIdx, ir, il)
% Base = pole sphere positions; Tip = base lifted radially outward from each
% pole's own hemisphere sphere center (the two spheres are offset, so a global
% center would not point outward).
    poleIdx = poleIdx(:);
    Base = sphV(poleIdx, :);
    cL = mean(sphV(il,:), 1);  cR = mean(sphV(ir,:), 1);
    rL = mean(sqrt(sum((sphV(il,:) - cL).^2, 2)));
    rR = mean(sqrt(sum((sphV(ir,:) - cR).^2, 2)));
    Tip = zeros(numel(poleIdx), 3);
    for i = 1:numel(poleIdx)
        p = poleIdx(i);
        if ismember(p, il), c = cL; rr = rL; else, c = cR; rr = rR; end
        rad = sphV(p,:) - c;
        rad = rad / max(sqrt(sum(rad.^2)), eps);
        Tip(i,:) = Base(i,:) + 0.15 * rr * rad;
    end
end


%% ===== GUI: registration sphere + pole pins =====
function hFig = ViewFigure(ManifoldFile)
    hFig = [];
    % --- load + validate node ---
    Full = file_fullpath(ManifoldFile);
    if ~file_exist(Full)
        bst_error('Manifold file not found.', 'View registration sphere', 0);
        return;
    end
    M = load(Full);
    if ~isfield(M,'ParentSurface') || isempty(M.ParentSurface) ...
            || ~isfield(M,'Embedded') || numel(M.Embedded) ~= 2 ...
            || ~isfield(M,'Gauge') || numel(M.Gauge) ~= 2
        bst_error(['Not a valid manifold node ' 10 '(need ParentSurface, Embedded[1x2], Gauge[1x2]).'], 'View registration sphere', 0);
        return;
    end
    Surface = M.ParentSurface;
    TessMat = in_tess_bst(Surface);
    nVert   = size(TessMat.Vertices, 1);

    % --- pole global vertex ids (reuse view_manifold's derivation) ---
    G = view_manifold('DeriveVertexFrame', M.Embedded, M.Gauge, nVert);
    poleIdx = G.Sing(:);

    % --- open the standard registration sphere figure ---
    hFig = view_surface_sphere(Surface, 'orig');
    if isempty(hFig)
        return;   % view_surface_sphere already reported (e.g. no Reg.Sphere)
    end
    set(hFig, 'Name', ['Manifold registration: ' Surface]);
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');

    % --- read the displayed sphere vertex positions (SCS, L/R offset) ---
    hPatch = findobj(hAxes, 'Type', 'patch');
    sphV = [];
    for k = 1:numel(hPatch)
        Vk = get(hPatch(k), 'Vertices');
        if size(Vk,1) == nVert, sphV = Vk; break; end
    end
    if isempty(sphV)
        % Could not match the sphere patch: leave the bare sphere figure.
        return;
    end

    % --- pole lollipop geometry (per-hemisphere radial lift) ---
    [ir, il] = tess_hemisplit(TessMat);
    Base = []; Tip = [];
    if ~isempty(poleIdx)
        [Base, Tip] = RegSingGlyphs(sphV, poleIdx, ir, il);
    end

    % --- state + keyboard ---
    showSing = true;
    colSing  = [0 0.45 1];
    KeyPressFcn_bak = get(hFig, 'KeyPressFcn');
    set(hFig, 'KeyPressFcn', @KeyPress_Callback);

    hold(hAxes, 'on');
    DrawSing();

    % ===== NESTED: draw / redraw the pole pins =====
    function DrawSing()
        delete(findobj(hAxes, '-depth', 1, '-regexp', 'Tag', '^manifoldRegSing'));
        if showSing && ~isempty(Tip)
            line([Base(:,1)';Tip(:,1)'], [Base(:,2)';Tip(:,2)'], [Base(:,3)';Tip(:,3)'], ...
                'Parent', hAxes, 'Color', colSing, 'LineWidth', 1.5, 'Tag', 'manifoldRegSingStem');
            plot3(Tip(:,1), Tip(:,2), Tip(:,3), 'o', 'Parent', hAxes, ...
                'MarkerFaceColor', colSing, 'MarkerEdgeColor', [.2 .2 .2], ...
                'MarkerSize', 10, 'LineStyle', 'none', 'Tag', 'manifoldRegSing');
        end
    end

    % ===== NESTED: keyboard =====
    function KeyPress_Callback(h, ev)
        switch ev.Key
            case 'p'
                showSing = ~showSing;
                DrawSing();
            case 'h'
                java_dialog('msgbox', ['<HTML><TABLE>' ...
                    '<TR><TD><B>P</B></TD><TD>Toggle singularity pins</TD></TR>' ...
                    '<TR><TD><B>0-9</B></TD><TD>Change view</TD></TR>' ...
                    '</TABLE>'], 'Manifold registration shortcuts', [], 0);
            otherwise
                if ~isempty(KeyPressFcn_bak), KeyPressFcn_bak(h, ev); end
        end
    end
end
