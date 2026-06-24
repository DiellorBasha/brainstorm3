function varargout = view_dynamics( varargin )
% VIEW_DYNAMICS: Display a nested atom table on the cortex + open the dynamics panel.
%
% Loads a dynamics table (db_template('dynamicsmat')), opens its cortex surface, draws each
% atom group's spatial occurrences as colored markers (one color per group, offset along the
% vertex normal to clear the surface), and opens panel_bst_dynamics (the Record-style atom
% panel: group tree + occurrence list + File/Atoms menus). Extended/temporal groups (windows
% with no vertex) draw no markers but appear in the tree for navigation.
%
% USAGE:
%   [hFig, T] = view_dynamics(DynamicsFile)        % load a table, open cortex + panel
%               view_dynamics('Redraw', hFig, T)   % redraw markers from an edited table
%
% SEE ALSO: panel_bst_dynamics, bst_dynamics, process_source_atoms
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

    % --- Redraw verb: view_dynamics('Redraw', hFig, T) ---
    if (nargin >= 1) && ischar(varargin{1}) && strcmp(varargin{1}, 'Redraw')
        Redraw(varargin{2:end});
        return;
    end

    % ===== LOAD TABLE =====
    DynamicsFile = varargin{1};
    T = bst_dynamics('Load', DynamicsFile);
    if isempty(T.Groups)
        error('view_dynamics: no atom groups in %s', DynamicsFile);
    end
    SurfaceFile = T.SurfaceFile;
    if isempty(SurfaceFile), SurfaceFile = T.Groups(1).SurfaceFile; end
    if isempty(SurfaceFile)
        error('view_dynamics: table has no SurfaceFile.');
    end

    % ===== OPEN CORTEX + DRAW MARKERS =====
    hFig = view_surface(SurfaceFile);
    setappdata(hFig, 'DynamicsFile', DynamicsFile);
    Redraw(hFig, T);

    % ===== OPEN THE PANEL (rebuild fresh) =====
    try, gui_hide('Dynamics'); catch, end %#ok<CTCH>
    gui_show('panel_bst_dynamics', 'JavaWindow', sprintf('Dynamics: %s', T.Comment), [], 0, 0, 0);
    panel_bst_dynamics('SetTarget', hFig, T);

    if (nargout >= 1), varargout{1} = hFig; end
    if (nargout >= 2), varargout{2} = T;    end
end


%% ===== REDRAW MARKERS FROM A (possibly edited) TABLE =====
function Redraw(hFig, T)
    if isempty(hFig) || ~ishandle(hFig), return; end
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    if isempty(hAxes), return; end
    hAxes = hAxes(1);
    set(hAxes, 'NextPlot', 'add');   % low-level line() never resets the cortex
    % Cache the surface (normals) per figure to avoid reloading on every edit
    Surf = getappdata(hFig, 'DynamicsSurf');
    if isempty(Surf)
        SurfaceFile = T.SurfaceFile;  if isempty(SurfaceFile), SurfaceFile = T.Groups(1).SurfaceFile; end
        Surf = in_tess_bst(SurfaceFile, 0);
        setappdata(hFig, 'DynamicsSurf', Surf);
    end
    hasNorm = isfield(Surf, 'VertNormals') && ~isempty(Surf.VertNormals);
    % Clear old markers + selection
    delete(findobj(hAxes, '-regexp', 'Tag', '^AtomMarker'));
    delete(findobj(hAxes, 'Tag', 'AtomSel'));
    % One marker set per spatial group
    GroupsPosOff = cell(1, numel(T.Groups));
    for g = 1:numel(T.Groups)
        G = T.Groups(g);
        if isempty(G.pos)
            GroupsPosOff{g} = zeros(0,3);   % temporal-only group (window): no markers
            continue;
        end
        vtx = double(G.vertices(:));
        if hasNorm, nrm = Surf.VertNormals(vtx, :);
        else,       nrm = G.pos ./ max(sqrt(sum(G.pos.^2,2)), eps); end
        po = G.pos + 0.002 * nrm;           % ~2 mm out, clears the surface
        GroupsPosOff{g} = po;
        col = G.color;  if isempty(col), col = [1 0 1]; end
        line(po(:,1), po(:,2), po(:,3), 'Parent', hAxes, ...
            'Marker','o', 'MarkerFaceColor',col, 'MarkerEdgeColor',[0 0 0], ...
            'MarkerSize',6, 'LineStyle','none', 'Tag', sprintf('AtomMarker%d', g));
    end
    setappdata(hFig, 'GroupsPosOff', GroupsPosOff);
    % Selection marker (hidden until an occurrence is picked)
    line(NaN, NaN, NaN, 'Parent', hAxes, ...
        'Marker','o', 'MarkerSize',13, 'MarkerEdgeColor',[1 1 0], ...
        'LineWidth',2, 'LineStyle','none', 'Visible','off', 'Tag','AtomSel');
end
