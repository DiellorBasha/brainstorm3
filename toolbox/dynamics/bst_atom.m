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


%% ===== SET: write one axis localization into the group =====
function G = Set(G, axis, occ, loc)
    if (nargin < 4), error('bst_atom:Set','Set requires (G, axis, occ, loc).'); end
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
            error('bst_atom:Set', 'Unknown axis "%s".', axis);
    end
end

%% ===== EVALUATE: realize an atom as a joint time-cortex WAVELET field =====
function [W, gv] = Evaluate(G, occ, ax, kernelName, kernelParams) %#ok<DEFNU>
    % Realize the atom occurrence as a smooth wavelet field W [nLoc x nT] over the eigenbasis support gv:
    % propagate the source-axis seed through the cortex eigenbasis with a dynamic eigenfilter over the
    % time axis. A dynamic (eigen-time) kernel gives a TIME-VARYING scout (a propagating wave); a static
    % (eigen-frequency) kernel gives a smooth scout constant in time. This is the canonical wavelet; its
    % cortical level set is a Scout and its temporal level set is an Event (see Levelset).
    %
    % USAGE: [W,gv] = bst_atom('Evaluate', G, occ, ax, 'dampedwave', struct('alpha',6,'beta',0.4))
    %   ax = bst_dynamics('Axes', T, variant, nModes, tWin)
    if (nargin < 2) || isempty(occ),          occ          = 1;            end
    if (nargin < 4) || isempty(kernelName),   kernelName   = 'dampedwave'; end
    if (nargin < 5) || isempty(kernelParams), kernelParams = struct();     end
    seed = double(G.vertices(occ));
    gv   = ax.GlobalVertices{1};  loc = find(gv == seed, 1);
    if isempty(loc), error('bst_atom(''Evaluate''): seed vertex %d not in the eigenbasis support.', seed); end
    Phi = ax.Phi{1};  Lam = ax.Lambda{1};  M = ax.Mass{1};
    c0  = manifold_ft(Phi, M, full(sparse(loc,1,1,size(Phi,1),1)));    % seed in the eigenbasis [K x 1]
    kp  = kernelParams;  if ~isfield(kp,'lmax') || isempty(kp.lmax), kp.lmax = max(Lam); end
    g    = bst_eigfilter_kernel(kernelName, kp);
    meta = bst_eigfilter_kernel('info', kernelName);
    if isfield(meta,'domain') && strcmpi(meta.domain,'ts')            % dynamic (eigen-time) -> propagate
        W = manifold_ift(Phi, g(Lam, ax.tlag) .* c0);                 % [nLoc x nT]
    else                                                              % static (eigen-frequency) -> smooth scout
        W = repmat(manifold_ift(Phi, g(Lam) .* c0), 1, ax.nT);
    end
end

%% ===== LEVELSET: derive the hard indicators (Scout / Event) from the wavelet =====
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
% pad a [m x 3] position matrix to >= n rows with NaN
function p = i_pad_pos(p, n)
    if isempty(p), p = nan(n, 3); elseif (size(p,1) < n), p(end+1:n, :) = NaN; end
end


%% ===== state from extent =====
function s = i_state(extent)
    if ~isfinite(extent),   s = 'unlocalized';
    elseif (extent == 0),   s = 'point';
    else,                   s = 'window';
    end
end
