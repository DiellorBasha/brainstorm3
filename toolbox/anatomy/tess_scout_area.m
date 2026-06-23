function [Vertices, phi] = tess_scout_area(SurfaceFile, Seed, Radius, phi)
% TESS_SCOUT_AREA: Grow a cortical scout as the GEODESIC DISK of a given radius around a
% seed vertex, using the heat-method distance (nxr-compute). The scout is the interior of
% the geodesic isoline phi = Radius: { v : phi(v) <= Radius }. This is the metric-correct
% ("grow by N mm") analogue of the hop-count tess_scout_swell.
%
% USAGE:
%   [Vertices, phi] = tess_scout_area(SurfaceFile, Seed, Radius)
%   [Vertices, phi] = tess_scout_area(SurfaceFile, Seed, Radius, phi)   % reuse cached distance
%
% INPUTS:
%   - SurfaceFile : Brainstorm cortex surface file
%   - Seed        : seed vertex index (global, into the full surface)
%   - Radius      : geodesic radius [meters]; Vertices = { v : phi(v) <= Radius }
%   - phi         : (optional) precomputed geodesic distance field [nV x 1] from Seed,
%                   as returned by a prior call. Pass it to SKIP the heat solve (valid only
%                   for the SAME Seed + Surface) -- so repeated re-thresholding (e.g. a
%                   ctrl+scroll radius drag) is cheap.
%
% OUTPUTS:
%   - Vertices : [1 x k] vertex indices inside the geodesic disk (includes the Seed)
%   - phi      : [nV x 1] geodesic distance from Seed [meters]; Inf on the other hemisphere
%
% Geodesic distance is solved on the seed's hemisphere submesh only (heat distance does not
% cross the inter-hemispheric gap); the other hemisphere is left at Inf. The hemisphere
% split is atlas-based (tess_hemisplit, never conncomp), matching tess_manifold/tess_operators.
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

    % --- compute the geodesic distance field from the seed (unless a cached one was passed) ---
    if nargin < 4 || isempty(phi)
        phi = i_heat_distance(SurfaceFile, Seed);
    end

    % --- threshold: the geodesic disk of radius Radius ---
    Vertices = find(phi <= Radius)';
    % Guarantee the seed is always included (phi(Seed)==0, but guard fp/degenerate radius)
    if ~any(Vertices == Seed)
        Vertices = [Seed, Vertices];
    end
end


%% ===== heat-method geodesic distance from a seed (seed's hemisphere submesh) =====
function phi = i_heat_distance(SurfaceFile, Seed)
    % nxr-compute plugin (heat solver lives there)
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    if ~isOk
        error('tess_scout_area:nxrUnavailable', ...
            'tess_scout_area requires nxr-compute: %s', errMsg);
    end
    % load surface
    TessMat = in_tess_bst(SurfaceFile, 0);
    Vtx   = TessMat.Vertices;
    Fcs   = double(TessMat.Faces);
    nVtot = size(Vtx, 1);
    if Seed < 1 || Seed > nVtot
        error('tess_scout_area:badSeed', 'Seed vertex %d is out of range [1 %d].', Seed, nVtot);
    end
    % atlas-based hemisphere split (never conncomp); pick the seed's hemisphere
    [rH, lH] = tess_hemisplit(TessMat);
    if ismember(Seed, rH)
        vH = rH(:);
    elseif ismember(Seed, lH)
        vH = lH(:);
    else
        error('tess_scout_area:seedNoHemisphere', ...
            'Seed vertex %d is not assigned to a hemisphere (need a Structures atlas with lh/rh).', Seed);
    end
    % local hemisphere submesh (the order/convention used by tess_manifold/tess_operators)
    isV   = false(nVtot, 1);  isV(vH) = true;
    fMask = all(isV(Fcs), 2);
    mapV  = zeros(nVtot, 1);  mapV(vH) = 1:numel(vH);
    Vloc  = Vtx(vH, :);
    Floc  = mapV(Fcs(fMask, :));
    seedLoc = mapV(Seed);
    % heat-method geodesic distance from the (local) seed
    h = nxr_compute('create', Vloc, Floc);
    dCleanup = onCleanup(@() nxr_compute('destroy', h));   % free the context even on error
    d_local = nxr_compute('heat', h, seedLoc);             % [nVh x 1] distance from seed
    % scatter back to the full surface; other hemisphere = Inf
    phi = inf(nVtot, 1);
    phi(vH) = d_local;
end
