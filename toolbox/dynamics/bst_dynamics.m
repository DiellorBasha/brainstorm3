function varargout = bst_dynamics( varargin )
% BST_DYNAMICS: I/O and builders for the spatiotemporal sparse marker table.
%
% A "dynamics" table (db_template('dynamicsmat')) holds an array of atom GROUPS
% (db_template('atomgroup')). A group reuses the Events grouping system: a label,
% a color, and 'times' that are [1 x N] (simple, one point) or [2 x N] (extended,
% a window) -- plus the atom extensions: per-occurrence parallel arrays
% (vertices/pos/hemi/strength/...), group-level frequency/scale/Function
% coordinates, and a 'parent' label for NESTING (e.g. an extended "alpha (8-13
% Hz)" window containing simple "alpha_peak" / "alpha_trough" / ... children).
%
% Each group stores a REFERENCE (time x space x freq x scale) + provenance, not a
% copy of the data. Groups can be projected back to standard Events (drop the
% atom extensions) and Scouts (the vertices) for native Brainstorm navigation.
%
% USAGE:
%    T = bst_dynamics('New', Comment)                  % empty table
%    G = bst_dynamics('NewGroup', Label)               % empty atom group (template)
%    T = bst_dynamics('AddGroup', T, G)                % append a group (normalized)
%    [pos,color,gIdx,oIdx,labels] = bst_dynamics('Flatten', T)  % all spatial occurrences
%    T = bst_dynamics('Load', DynamicsFile)            % load + template-fill (forward-compat)
%    OutFile = bst_dynamics('Save', OutFile, T)        % save table to disk
%
% SEE ALSO: db_template('atomgroup'/'dynamicsmat'), process_source_atoms, view_dynamics
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


%% ===== NEW TABLE =====
function T = New(Comment)
    if (nargin < 1), Comment = 'dynamics'; end
    T = db_template('dynamicsmat');
    T.Comment = Comment;
end


%% ===== NEW GROUP =====
function G = NewGroup(Label)
    G = db_template('atomgroup');
    if (nargin >= 1), G.label = Label; end
end


%% ===== ADD GROUP =====
function T = AddGroup(T, G)
    if isempty(G)
        return;
    end
    G = struct_copy_fields(db_template('atomgroup'), G, 1);
    % Keep 'type' consistent with the times layout
    if ~isempty(G.times) && (size(G.times,1) == 2)
        G.type = 'extended';
    elseif ~isempty(G.times)
        G.type = 'simple';
    end
    if isempty(T.Groups)
        T.Groups = G;
    else
        T.Groups(end+1) = G;
    end
    T.nGroups = numel(T.Groups);
end


%% ===== FLATTEN (all SPATIAL occurrences across groups, in a stable order) =====
% Returns one row per occurrence that has a position (windows with no vertex are
% skipped). gIdx/oIdx index back into T.Groups(gIdx) column oIdx.
function [pos, color, gIdx, oIdx, labels] = Flatten(T)
    pos = zeros(0,3);  color = zeros(0,3);  gIdx = [];  oIdx = [];  labels = {};
    for g = 1:numel(T.Groups)
        G = T.Groups(g);
        if isempty(G.pos)
            continue;   % temporal-only group (e.g. an extended window): no markers
        end
        n   = size(G.pos, 1);
        col = G.color;  if isempty(col), col = [1 0 1]; end
        pos   = [pos;   G.pos];                  %#ok<AGROW>
        color = [color; repmat(col(:)', n, 1)];  %#ok<AGROW>
        gIdx  = [gIdx;  g*ones(n,1)];            %#ok<AGROW>
        oIdx  = [oIdx;  (1:n)'];                 %#ok<AGROW>
        for o = 1:n
            labels{end+1} = sprintf('%.3fs  v%d', G.times(1,o), G.vertices(o)); %#ok<AGROW>
        end
    end
end


%% ===== EXTREMA (surface-scalar local maxima/minima -> atom seeds) =====
% Pure detector for the storage step: the top-N local extrema of a per-vertex scalar
% field on a surface. The fresh extremum finder behind "Record at cursor".
%   field    [nV x 1] per-vertex scalar (e.g. Helmholtz Phi / Psi / |J|)
%   VertConn [nV x nV] sparse vertex adjacency (from in_tess_bst)
%   nPeaks   keep the top-N by value, per sign (default 3)
%   signed   true  -> maxima (charge +1) AND minima (charge -1)  [Phi sources/sinks, Psi vortices/anti]
%            false -> maxima only (charge +1)                    [|J| amplitude peaks]
% Returns ex.iVertex / ex.value / ex.charge : [1 x M] arrays (M <= nPeaks, or <= 2*nPeaks if signed).
function ex = Extrema(field, VertConn, nPeaks, signed) %#ok<DEFNU>
    if (nargin < 3) || isempty(nPeaks), nPeaks = 3;     end
    if (nargin < 4) || isempty(signed), signed = false; end
    field = double(field(:));
    vMax = i_local_ext(field, VertConn, nPeaks, +1);
    iV = vMax(:)';  val = field(vMax)';  chg = ones(1, numel(vMax));
    if signed
        vMin = i_local_ext(field, VertConn, nPeaks, -1);
        iV = [iV, vMin(:)'];  val = [val, field(vMin)'];  chg = [chg, -ones(1, numel(vMin))];
    end
    ex = struct('iVertex', iV, 'value', val, 'charge', chg);
end

% Top-N vertices that are local extrema of sgn*field (sgn=+1 maxima, -1 minima).
function vPk = i_local_ext(field, VertConn, nPeaks, sgn)
    f  = sgn * field(:);
    nV = numel(f);
    [ii, jj] = find(VertConn);                              % undirected edges (symmetric)
    if isempty(ii)
        nbMax = -inf(nV, 1);
    else
        nbMax = accumarray(ii, f(jj), [nV 1], @max, -inf);  % per-vertex neighbour maximum
    end
    cand = find(f >= nbMax);                                % >= all neighbours
    [~, ord] = sort(f(cand), 'descend');
    vPk = cand(ord(1:min(nPeaks, numel(cand))));
end


%% ===== ATTACH REGION (localize occurrence o with a geodesic region + seed) =====
% Pads the per-occurrence arrays to full length N = size(G.times,2), then writes occurrence o's
% seed (vertices/pos/hemi) and region. A time-only marker (empty vertices/pos) becomes partially
% localized: occurrence o gains a finite seed + region; the other occurrences stay NaN/[] (time only).
%   o           occurrence column to localize (1..N)
%   regionVerts [1 x k] cortex vertex indices of the geodesic disk (snapshot copy)
%   seed        scalar seed vertex index
%   pos         [1 x 3] seed position (SCS)
%   hemi        scalar 1=L 2=R
function G = AttachRegion(G, o, regionVerts, seed, pos, hemi) %#ok<DEFNU>
    N = size(G.times, 2);
    if (o < 1) || (o > N)
        error('bst_dynamics:AttachRegion', 'Occurrence %d out of range (N=%d).', o, N);
    end
    G.vertices = i_pad_row(G.vertices, N);
    G.hemi     = i_pad_row(G.hemi,     N);
    G.strength = i_pad_row(G.strength, N);
    G.charge   = i_pad_row(G.charge,   N);
    G.pos      = i_pad_pos(G.pos,      N);
    G.region   = i_pad_cell(G.region,  N);
    G.vertices(o) = double(seed);
    G.hemi(o)     = double(hemi);
    G.pos(o, :)   = double(pos(:)');
    G.region{o}   = double(regionVerts(:)');
    G.type = 'simple';
end

% Pad a [1 x m] numeric row to length N with NaN ([] -> all-NaN).
function v = i_pad_row(v, N)
    if isempty(v),          v = nan(1, N);
    elseif (numel(v) < N),  v(end+1:N) = NaN;  end
end
% Pad a [m x 3] position matrix to N rows with NaN ([] -> all-NaN).
function p = i_pad_pos(p, N)
    if isempty(p),           p = nan(N, 3);
    elseif (size(p,1) < N),  p(end+1:N, :) = NaN;  end
end
% Pad a {1 x m} cell to length N with [] ({} -> all-empty).
function c = i_pad_cell(c, N)
    if isempty(c),          c = cell(1, N);
    elseif (numel(c) < N),  c(end+1:N) = {[]};  end
end


%% ===== LOAD =====
function T = Load(DynamicsFile)
    % Phase 1: load by absolute path (the 'dynamics_' type is not yet registered
    % with file_gettype/file_fullpath; that comes with the Phase 2 tree node).
    T = load(DynamicsFile);
    T = struct_copy_fields(db_template('dynamicsmat'), T, 1);
    % Forward-compat: normalize each group against the current template.
    if ~isempty(T.Groups)
        tmpl = db_template('atomgroup');
        G = repmat(tmpl, 1, 0);
        for i = 1:numel(T.Groups)
            G(i) = struct_copy_fields(tmpl, T.Groups(i), 1);
        end
        T.Groups = G;
    end
    T.nGroups = numel(T.Groups);
end


%% ===== SAVE =====
function OutFile = Save(OutFile, T)
    T.nGroups = numel(T.Groups);
    bst_save(OutFile, T, 'v7');
    % Phase 1 has no dedicated tree node, so the table is not registered in the
    % DB; it is loadable by path via bst_dynamics('Load', ...). A tree node +
    % viewer come in Phase 2 (bst-java).
end
