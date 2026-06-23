function [hFig, T] = view_dynamics( DynamicsFile )
% VIEW_DYNAMICS: Display a spatiotemporal "atom" table on the cortex + open the list panel.
%
% Milestone-1 viewer for the joint Events+Scouts atom system. Loads a dynamics
% table (db_template('dynamicsmat')), opens its cortex surface, draws each atom
% as a colored point marker (grouped by category), and opens panel_dynamics so
% rows can be selected/highlighted. Atom positions are pushed slightly outward
% along the vertex normal so the markers clear the opaque surface.
%
% USAGE:  [hFig, T] = view_dynamics(DynamicsFile)
%
% INPUTS:
%   - DynamicsFile : path to a dynamics_*.mat table (bst_dynamics)
%
% SEE ALSO: panel_dynamics, bst_dynamics, process_source_atoms
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

    % ===== LOAD TABLE =====
    T = bst_dynamics('Load', DynamicsFile);
    if isempty(T.Atoms)
        error('view_dynamics: no atoms in %s', DynamicsFile);
    end
    if isempty(T.SurfaceFile)
        error('view_dynamics: table has no SurfaceFile.');
    end

    % ===== OPEN CORTEX =====
    hFig = view_surface(T.SurfaceFile);
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    if isempty(hAxes)
        error('view_dynamics: could not find the 3D axes.');
    end
    % Keep the axes additive so low-level line() never resets the cortex patch
    set(hAxes, 'NextPlot', 'add');

    % ===== MARKER POSITIONS (offset outward along the vertex normal) =====
    Surf = in_tess_bst(T.SurfaceFile, 0);
    vtx  = double([T.Atoms.vertex]);
    pos  = Surf.Vertices(vtx, :);
    if isfield(Surf, 'VertNormals') && ~isempty(Surf.VertNormals)
        nrm = Surf.VertNormals(vtx, :);
    else
        nrm = pos ./ max(sqrt(sum(pos.^2,2)), eps);   % fallback: radial
    end
    posOff = pos + 0.002 * nrm;                        % ~2 mm out, clears the surface
    setappdata(hFig, 'AtomOffsetPos', posOff);

    % ===== DRAW MARKERS (grouped by category for per-color rendering) =====
    cats = {T.Atoms.category};
    [uCat, ~, ic] = unique(cats);
    for c = 1:numel(uCat)
        idx = (ic == c);
        col = T.Atoms(find(idx,1)).color;
        if isempty(col), col = [1 0 1]; end
        line(posOff(idx,1), posOff(idx,2), posOff(idx,3), ...
            'Parent',          hAxes, ...
            'Marker',          'o', ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', [0 0 0], ...
            'MarkerSize',      6, ...
            'LineStyle',       'none', ...
            'Tag',             'AtomMarker');
    end
    % Selection marker (hidden until a row is picked)
    line(NaN, NaN, NaN, 'Parent', hAxes, ...
        'Marker', 'o', 'MarkerSize', 13, 'MarkerEdgeColor', [1 1 0], ...
        'LineWidth', 2, 'LineStyle', 'none', 'Visible', 'off', 'Tag', 'AtomSel');

    % ===== OPEN THE LIST PANEL =====
    gui_show('panel_dynamics', 'JavaWindow', sprintf('Dynamics: %s', T.Comment), [], 0, 0, 0);
    panel_dynamics('SetTarget', hFig, T.Atoms);
end
