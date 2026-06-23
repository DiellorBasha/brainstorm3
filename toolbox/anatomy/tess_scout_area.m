function [Vertices, phi] = tess_scout_area(SurfaceFile, Seed, Radius, phi)
% TESS_SCOUT_AREA: Grow a cortical scout as the GEODESIC DISK of a given radius around a
% seed vertex, using the heat-method distance (nxr-compute). The scout is the interior of
% the geodesic isoline phi = Radius: { v : phi(v) <= Radius }. This is the metric-correct
% ("grow by N mm") analogue of the hop-count tess_scout_swell.
%
% USAGE:
%   [Vertices, phi] = tess_scout_area(SurfaceFile, Seed, Radius)
%   [Vertices, phi] = tess_scout_area(SurfaceFile, Seed, Radius, phi)   % reuse cached distance
%   tess_scout_area('prewarm', SurfaceFile)   % build+factorize both hemisphere solvers ahead of time
%   tess_scout_area('clear')                  % free the cached nxr contexts
%
% INPUTS:
%   - SurfaceFile : Brainstorm cortex surface file
%   - Seed        : seed vertex index (global, into the full surface)
%   - Radius      : geodesic radius [meters]; Vertices = { v : phi(v) <= Radius }
%   - phi         : (optional) precomputed geodesic distance field [nV x 1] from Seed,
%                   as returned by a prior call. Pass it to SKIP the heat solve (valid only
%                   for the SAME Seed + Surface).
%
% OUTPUTS:
%   - Vertices : [1 x k] vertex indices inside the geodesic disk (includes the Seed)
%   - phi      : [nV x 1] geodesic distance from Seed [meters]; Inf on the other hemisphere
%
% PERFORMANCE: the nxr context for each hemisphere (mesh + the heat-method Cholesky factor
% geometry-central caches internally) is built ONCE and kept alive in a persistent cache,
% keyed by SurfaceFile. Repeated calls -- even with different seeds -- reuse the factor and
% only run a back-substitution. 'prewarm' pays the one-time factorization up front (e.g. when
% the Area tool is activated) so the first click is instant; 'clear' frees the contexts. The
% cache is rebuilt automatically when SurfaceFile changes.
%
% Geodesic distance is solved on the seed's hemisphere submesh only (heat distance does not
% cross the inter-hemispheric gap); the other hemisphere is left at Inf. The hemisphere split
% is atlas-based (tess_hemisplit, never conncomp), matching tess_manifold/tess_operators.
%
% SEE ALSO: tess_scout_swell, tess_hemisplit, tess_manifold, nxr_compute('heat')

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

    % Persistent per-surface context cache (mesh + heat Cholesky factor live in the nxr handle)
    persistent CTX;
    Vertices = [];

    % ----- verbs: clear / prewarm -----
    if ischar(SurfaceFile)
        switch lower(SurfaceFile)
            case 'clear'
                CTX = i_free(CTX);
                phi = [];
                return;
            case 'prewarm'
                % USAGE: tess_scout_area('prewarm', SurfaceFile)
                CTX = i_ensure_surface(CTX, Seed);          % Seed holds the SurfaceFile here
                for hh = 1:2
                    [CTX, ok] = i_build_handle(CTX, hh);
                    if ok
                        try, nxr_compute('heat', CTX.H(hh).handle, 1); catch, end  % dummy solve -> factorize+cache
                    end
                end
                phi = [];
                return;
            case 'path'
                % USAGE: pts = tess_scout_area('path', SurfaceFile, vStart, vEnd)
                %   The exact GEODESIC LINE between two vertices (geometry-central FlipOut /
                %   edge-flip, Sharp & Crane 2020) as a [nPts x 3] polyline -- reuses the same
                %   cached per-hemisphere context as the geodesic-disk path. Both endpoints must
                %   be in the same hemisphere. Returned as the FIRST output.
                [CTX, Vertices] = i_trace_path(CTX, Seed, Radius, phi);   % Seed=Surface, Radius=vStart, phi=vEnd
                phi = [];
                return;
        end
    end

    % ----- compute the geodesic distance field from the seed (unless a cached one was passed) -----
    if nargin < 4 || isempty(phi)
        [phi, CTX] = i_heat_distance(CTX, SurfaceFile, Seed);
    end

    % ----- threshold: the geodesic disk of radius Radius -----
    Vertices = find(phi <= Radius)';
    if ~any(Vertices == Seed)        % guard fp / degenerate radius: seed is always in
        Vertices = [Seed, Vertices];
    end
end


%% ===== heat-method geodesic distance from a seed (reuses the cached per-hemisphere solver) =====
function [phi, CTX] = i_heat_distance(CTX, SurfaceFile, Seed)
    CTX = i_ensure_surface(CTX, SurfaceFile);
    hh = i_hemi_of(CTX, Seed);
    if hh == 0
        error('tess_scout_area:seedNoHemisphere', ...
            'Seed vertex %d is not assigned to a hemisphere (need a Structures atlas with lh/rh).', Seed);
    end
    seedLoc = CTX.H(hh).mapV(Seed);
    [CTX, ok] = i_build_handle(CTX, hh);
    if ~ok
        error('tess_scout_area:nxrCreateFailed', 'Could not build the nxr context for hemisphere %d.', hh);
    end
    % heat solve; rebuild the context once if the handle went stale (e.g. the MEX was reloaded)
    try
        d_local = nxr_compute('heat', CTX.H(hh).handle, seedLoc);
    catch
        CTX.H(hh).handle = [];
        [CTX, ok] = i_build_handle(CTX, hh);
        if ~ok, error('tess_scout_area:heatFailed', 'Heat solve failed and the context could not be rebuilt.'); end
        d_local = nxr_compute('heat', CTX.H(hh).handle, seedLoc);
    end
    phi = inf(CTX.nVtot, 1);
    phi(CTX.H(hh).vH) = d_local;
end


%% ===== geodesic LINE between two same-hemisphere vertices (FlipOut, reuses the cached context) =====
function [CTX, pts] = i_trace_path(CTX, SurfaceFile, vA, vB)
    CTX = i_ensure_surface(CTX, SurfaceFile);
    hhA = i_hemi_of(CTX, vA);
    hhB = i_hemi_of(CTX, vB);
    if hhA == 0 || hhB == 0
        error('tess_scout_area:endpointNoHemisphere', 'A geodesic endpoint is not assigned to a hemisphere.');
    end
    if hhA ~= hhB
        error('tess_scout_area:crossHemisphere', ...
            'The two geodesic endpoints are in different hemispheres; a within-surface geodesic does not cross the gap.');
    end
    hh = hhA;
    [CTX, ok] = i_build_handle(CTX, hh);
    if ~ok, error('tess_scout_area:nxrCreateFailed', 'Could not build the nxr context for hemisphere %d.', hh); end
    aLoc = CTX.H(hh).mapV(vA);
    bLoc = CTX.H(hh).mapV(vB);
    try
        pts = nxr_compute('tracePath', CTX.H(hh).handle, aLoc, bLoc);
    catch
        CTX.H(hh).handle = [];
        [CTX, ok] = i_build_handle(CTX, hh);
        if ~ok, error('tess_scout_area:pathFailed', 'tracePath failed and the context could not be rebuilt.'); end
        pts = nxr_compute('tracePath', CTX.H(hh).handle, aLoc, bLoc);
    end
end


%% ===== which hemisphere a global vertex belongs to (1=right, 2=left, 0=none); O(1) =====
function hh = i_hemi_of(CTX, v)
    n = CTX.nVtot;
    if (v >= 1) && (v <= n) && (CTX.H(1).mapV(v) > 0)
        hh = 1;
    elseif (v >= 1) && (v <= n) && (CTX.H(2).mapV(v) > 0)
        hh = 2;
    else
        hh = 0;
    end
end


%% ===== (re)build the per-surface geometry cache (frees old contexts when the surface changes) =====
function CTX = i_ensure_surface(CTX, SurfaceFile)
    if ~isempty(CTX) && strcmp(CTX.SurfaceFile, SurfaceFile)
        return;     % already loaded for this surface
    end
    % nxr-compute plugin (the heat solver lives there)
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    if ~isOk
        error('tess_scout_area:nxrUnavailable', 'tess_scout_area requires nxr-compute: %s', errMsg);
    end
    CTX = i_free(CTX);                       % drop any previous surface's contexts
    TessMat = in_tess_bst(SurfaceFile, 0);
    Vtx   = TessMat.Vertices;
    Fcs   = double(TessMat.Faces);
    nVtot = size(Vtx, 1);
    [rH, lH] = tess_hemisplit(TessMat);      % atlas-based split (never conncomp)
    hemis = {rH(:), lH(:)};
    H = repmat(struct('vH',[], 'mapV',[], 'Vloc',[], 'Floc',[], 'handle',[]), 1, 2);
    for hh = 1:2
        vH = hemis{hh};
        isV = false(nVtot,1);  isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        mapV  = zeros(nVtot,1); mapV(vH) = 1:numel(vH);
        H(hh).vH   = vH;
        H(hh).mapV = mapV;
        H(hh).Vloc = Vtx(vH, :);
        H(hh).Floc = mapV(Fcs(fMask, :));
        H(hh).handle = [];                   % built lazily on first use (or by 'prewarm')
    end
    CTX = struct('SurfaceFile', SurfaceFile, 'nVtot', nVtot, 'H', H);
end


%% ===== lazily create the nxr context (mesh) for one hemisphere =====
function [CTX, ok] = i_build_handle(CTX, hh)
    ok = true;
    if ~isempty(CTX.H(hh).handle)
        return;
    end
    if isempty(CTX.H(hh).vH) || isempty(CTX.H(hh).Floc)
        ok = false;  return;                 % empty hemisphere
    end
    CTX.H(hh).handle = nxr_compute('create', CTX.H(hh).Vloc, CTX.H(hh).Floc);
end


%% ===== free all cached nxr contexts =====
function CTX = i_free(CTX)
    if ~isempty(CTX) && isfield(CTX, 'H')
        for hh = 1:numel(CTX.H)
            if ~isempty(CTX.H(hh).handle)
                try, nxr_compute('destroy', CTX.H(hh).handle); catch, end %#ok<CTCH>
            end
        end
    end
    CTX = [];
end
