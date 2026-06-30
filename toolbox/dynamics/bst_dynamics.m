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
% NaN-padded (unlocalized) occurrences produced by AttachRegion are silently skipped.
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
            if any(~isfinite(G.pos(o,:))), continue; end   % skip unlocalized (NaN) occurrences
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


%% ===== AXES: resolve the joint time-cortex axes (results-style binding) =====
function ax = Axes(T, variant, nModes, tWin) %#ok<DEFNU>
    % Resolve the JOINT (time, cortex) axes a dynamics table is bound to, mirroring how a results file
    % binds a SurfaceFile + a Time vector. The CORTEX axis is an eigenbasis of the bound SurfaceFile
    % (tess_eigen of the requested operator variant -> Phi, Lambda, mm via 2*pi/sqrt(lambda)); the TIME
    % axis is the Time vector + Fs of the bound DataFile/ResultsFile; omega is the temporal FFT grid (Hz).
    % This is the source-of-truth the JTV wavelet orchestrators read from (no separate bst_jtv needed).
    %
    % USAGE: ax = bst_dynamics('Axes', T, 'LB-Connectome', 300, [tStart tEnd])
    if (nargin < 2) || isempty(variant), variant = 'Laplace-Beltrami'; end
    if (nargin < 3) || isempty(nModes),  nModes  = 200; end
    if (nargin < 4), tWin = []; end
    if isempty(T.SurfaceFile), error('bst_dynamics(''Axes''): no SurfaceFile bound.'); end
    % --- time axis (Time/Fs of the bound recording; raw -> reconstructed + windowed) ---
    timeFile = T.DataFile;
    if isempty(timeFile) && ~isempty(T.Groups)
        if isfield(T.Groups,'ResultsFile') && ~isempty(T.Groups(1).ResultsFile), timeFile = T.Groups(1).ResultsFile;
        elseif isfield(T.Groups,'DataFile') && ~isempty(T.Groups(1).DataFile),   timeFile = T.Groups(1).DataFile; end
    end
    if isempty(timeFile), error('bst_dynamics(''Axes''): no DataFile/ResultsFile bound for the time axis.'); end
    [Time, ~] = i_load_time(timeFile, tWin);
    % --- delegate the canonical eigenbasis x time x freq assembly to bst_eigen('Axes') ---
    ax = bst_eigen('Axes', struct('SurfaceFile',T.SurfaceFile, 'Variant',variant, 'nModes',nModes, 'Time',Time));
    ax.TimeFile = timeFile;
end

function [Time, Fs] = i_load_time(file, tWin)
    % Resolve a >=2-sample time vector (s) + Fs (Hz) from a results/data/raw file. A kernel results file
    % follows to its DataFile; a raw (continuous) file is reconstructed from the sFile sampling rate and
    % windowed (default first 4 s) - the whole continuous recording is impractical for an FFT.
    if (nargin < 2), tWin = []; end
    ft = file_gettype(file);
    if any(strcmpi(ft, {'results','link','presults'}))
        R = in_bst_results(file, 0, 'Time', 'DataFile');
        if numel(R.Time) > 2
            Time = R.Time;
        elseif ~isempty(R.DataFile)
            [Time, Fs] = i_load_time(R.DataFile, tWin); return;
        else
            error('bst_dynamics(''Axes''): results file has no usable time axis.');
        end
    else
        D = in_bst_data(file, 'Time', 'F');
        if numel(D.Time) > 2
            Time = D.Time;                                         % stored epoch/average vector
        elseif isstruct(D.F) && isfield(D.F,'prop')                % raw sFile -> reconstruct
            t = D.F.prop.times; fs = D.F.prop.sfreq;
            if isempty(tWin), tWin = [t(1), min(t(2), t(1) + 4)]; end   % default 4 s window
            Time = max(t(1), tWin(1)) : 1/fs : min(t(2), tWin(2));
        else
            Time = D.Time;
        end
    end
    if numel(Time) < 2, error('bst_dynamics(''Axes''): could not resolve a >=2-sample time axis.'); end
    Fs = 1 / (Time(2) - Time(1));
    if ~isempty(tWin) && (numel(Time) > 2)                          % window a stored vector
        sel = (Time >= tWin(1)) & (Time <= tWin(2));
        if any(sel), Time = Time(sel); end
    end
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


%% ===== MARKER LOCALISATION ACCESSOR (moved from bst_atom) =====
% Uniform (center, extent, weighting) accessor over an atom group's four axes (time/freq/source/scale):
% maps between the heterogeneous stored fields and a uniform Localization struct so the panel/detectors
% treat all axes identically. (Localization is a property of the stored MARKER, distinct from the JTV
% wavelet REALISER bst_eigenfilter('Atom') and the eigenbasis ax from bst_dynamics('Axes')/bst_eigen('Axes').)

% Canonical axis metadata (name/perOcc/unit). Named AxisMeta to avoid clash with Axes(T,...) above.
function A = AxisMeta() %#ok<DEFNU>
    A = struct('name',   {'time','freq','source','scale'}, ...
               'perOcc', {true,   false,  true,    false}, ...
               'unit',   {'s',    'Hz',   'vertex','eigenvalue'});
end

% New (empty) localization struct.
function loc = NewLoc(axis) %#ok<DEFNU>
    if (nargin < 1), axis = ''; end
    loc = struct('axis',axis, 'center',NaN, 'extent',NaN, 'weighting','hard', ...
                 'label','', 'state','unlocalized', 'pos',[]);
end

% Read one axis localization from the group.
function loc = Get(G, axis, occ) %#ok<DEFNU>
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
            error('bst_dynamics:Get', 'Unknown axis "%s".', axis);
    end
end

% Write one axis localization into the group.
function G = Set(G, axis, occ, loc) %#ok<DEFNU>
    if (nargin < 4), error('bst_dynamics:Set','Set requires (G, axis, occ, loc).'); end
    if isempty(occ), occ = 1; end
    c = loc.center;  w = loc.extent;  if ~isfinite(w), w = 0; end
    switch axis
        case 'time'
            G.times = i_pad_cols(G.times, occ);
            if (w > 0) && (size(G.times,1) < 2)
                G.times = [G.times; G.times];          % promote simple->extended (others zero-width)
            end
            if (size(G.times,1) >= 2)
                G.times(1, occ) = c - w;  G.times(2, occ) = c + w;
            else
                G.times(1, occ) = c;
            end
            G.type = i_type(G.times);
        case 'freq'
            G.band = [c - w, c + w];
            if ~isempty(loc.label), G.bandName = loc.label; end
        case 'source'
            G.vertices = i_pad_vec(G.vertices, occ);  G.vertices(occ) = c;
            G.radius   = i_pad_vec(G.radius,   occ);  G.radius(occ)   = w;
            if ~isempty(loc.pos)
                G.pos = i_pad_pos(G.pos, occ);  G.pos(occ, :) = loc.pos(:)';
            end
        case 'scale'
            G.scale = [c - w, c + w];
            if ~isempty(loc.label), G.scaleName = loc.label; end
        otherwise
            error('bst_dynamics:Set', 'Unknown axis "%s".', axis);
    end
end

% group type consistent with the times row count
function t = i_type(times)
    if (size(times,1) >= 2), t = 'extended'; else, t = 'simple'; end
end
% pad a [r x m] times matrix to >= n columns with NaN
function M = i_pad_cols(M, n)
    if isempty(M), M = nan(1, n); elseif (size(M,2) < n), M(:, end+1:n) = NaN; end
end
% pad a [1 x m] row vector to >= n with NaN
function v = i_pad_vec(v, n)
    if isempty(v), v = nan(1, n); elseif (numel(v) < n), v(end+1:n) = NaN; end
end
% state from extent
function s = i_state(extent)
    if ~isfinite(extent),   s = 'unlocalized';
    elseif (extent == 0),   s = 'point';
    else,                   s = 'window';
    end
end


%% ===== LEVELSET: derive the hard indicators (Scout / Event) from a wavelet field (moved from bst_atom) =====
function LS = Levelset(W, gv, thr, iRef) %#ok<DEFNU>
    % A Scout is a level set of the CORTICAL wavelet (at the peak-energy time iRef); an Event is a level
    % set of the TEMPORAL wavelet (energy over time). thr = fraction of the per-axis max (default 0.5).
    if (nargin < 3) || isempty(thr), thr = 0.5; end
    e = sum(W.^2, 1);
    if (nargin < 4) || isempty(iRef), [~, iRef] = max(e); end
    wRef = abs(W(:, iRef));
    LS.scoutVertices = gv(wRef >= thr * max(wRef));     % cortical level set -> Scout (global vertices)
    LS.eventSamples  = find(e   >= thr * max(e));        % temporal level set -> Event (sample indices)
    LS.iRef          = iRef;
end


%% ===== ATOM FROM KERNEL =====
function G = AtomFromKernel(ax, kernelName, kernelParams, seed, thr) %#ok<DEFNU>
    % Realise an eigfilter atom from its generator and threshold it into a populated atom group:
    % Scout = spatial level set (region), Event = temporal level set (times). Atom = thresholded filter.
    if (nargin < 6) || isempty(thr), thr = 0.5; end
    [W, gv] = bst_eigenfilter('Atom', ax, kernelName, kernelParams, seed);
    LS = Levelset(W, gv, thr);
    G = NewGroup(sprintf('%s @vtx%d', kernelName, seed));
    G.vertices = seed;
    if isfield(ax,'SurfaceFile'), G.SurfaceFile = ax.SurfaceFile; end
    G.region   = {LS.scoutVertices(:)'};
    G.times    = [ax.tlag(LS.eventSamples(1)); ax.tlag(LS.eventSamples(end))];
    G.type     = 'extended';
    G.KernelName = kernelName;  G.KernelParams = kernelParams;  G.Threshold = thr;
end
