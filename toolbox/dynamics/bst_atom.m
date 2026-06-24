function varargout = bst_atom( varargin )
% BST_ATOM: uniform (center, extent, weighting) localization accessor over an atom group.
%
% Every atom axis (time, frequency, source, scale) localizes the same way -- a center plus a
% window (extent), or a bare point. This accessor maps between an atomgroup's heterogeneous
% stored fields and a uniform Localization struct, so the panel and detectors can treat all
% four axes identically. I/O-free (operates on an in-memory group).
%
% Localization struct:
%   .axis      'time'|'freq'|'source'|'scale'
%   .center    numeric center (s | Hz | vertex id | eigenvalue)
%   .extent    numeric half-window (s | Hz | metres geodesic radius | eigenvalue); 0 = point
%   .weighting 'hard' (default) | 'soft'   (soft = wavelet decay; reserved, future)
%   .label     optional human-readable name ('' today; atlas layer is future)
%   .state     'unlocalized' | 'point' | 'window'
%   .pos       [1x3] seed position, source axis only (else [])
%
% USAGE:
%   A   = bst_atom('Axes')                 % axis metadata (name/perOcc/unit)
%   loc = bst_atom('NewLoc', axis)         % empty localization
%   loc = bst_atom('Get', G, axis, occ)    % read (occ default 1; ignored for group axes)
%   G   = bst_atom('Set', G, axis, occ, loc)  % write (Task 2)
%
% SEE ALSO: bst_dynamics, panel_bst_dynamics
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


%% ===== AXES: canonical axis metadata =====
function A = Axes()
    A = struct('name',   {'time','freq','source','scale'}, ...
               'perOcc', {true,   false,  true,    false}, ...
               'unit',   {'s',    'Hz',   'vertex','eigenvalue'});
end


%% ===== NEW (empty) localization =====
function loc = NewLoc(axis)
    if (nargin < 1), axis = ''; end
    loc = struct('axis',axis, 'center',NaN, 'extent',NaN, 'weighting','hard', ...
                 'label','', 'state','unlocalized', 'pos',[]);
end


%% ===== GET: read one axis localization from the group =====
function loc = Get(G, axis, occ)
    if (nargin < 3) || isempty(occ), occ = 1; end
    loc = NewLoc(axis);
    switch axis
        case 'time'
            if isempty(G.times) || (occ > size(G.times,2)), return; end
            col = G.times(:, occ);
            if any(~isfinite(col)), return; end
            if (size(G.times,1) >= 2)
                loc.center = mean(col(1:2));  loc.extent = (col(2) - col(1)) / 2;
            else
                loc.center = col(1);          loc.extent = 0;
            end
            loc.state = i_state(loc.extent);
        case 'freq'
            if (numel(G.band) < 2), return; end
            loc.center = mean(G.band(1:2));  loc.extent = (G.band(2) - G.band(1)) / 2;
            loc.label  = G.bandName;         loc.state  = i_state(loc.extent);
        case 'source'
            if isempty(G.vertices) || (occ > numel(G.vertices)) || ~isfinite(G.vertices(occ)), return; end
            loc.center = double(G.vertices(occ));
            if ~isempty(G.pos) && (occ <= size(G.pos,1)), loc.pos = G.pos(occ, :); end
            hasR = isfield(G,'radius') && ~isempty(G.radius) && (occ <= numel(G.radius)) && isfinite(G.radius(occ));
            hasReg = isfield(G,'region') && ~isempty(G.region) && (occ <= numel(G.region)) && ~isempty(G.region{occ});
            if hasR
                loc.extent = G.radius(occ);  loc.state = i_state(loc.extent);
            elseif hasReg
                loc.extent = NaN;  loc.state = 'window';   % region materialized but radius unrecorded
            else
                loc.extent = 0;    loc.state = 'point';
            end
        case 'scale'
            if (numel(G.scale) < 2), return; end
            loc.center = mean(G.scale(1:2));  loc.extent = (G.scale(2) - G.scale(1)) / 2;
            loc.label  = G.scaleName;         loc.state  = i_state(loc.extent);
        otherwise
            error('bst_atom:Get', 'Unknown axis "%s".', axis);
    end
end


%% ===== state from extent =====
function s = i_state(extent)
    if ~isfinite(extent),   s = 'unlocalized';
    elseif (extent == 0),   s = 'point';
    else,                   s = 'window';
    end
end
