function Metrics = bst_resolution_metrics(Kernel, Leadfield, GridLoc)
% BST_RESOLUTION_METRICS: Point-spread metrics from the resolution matrix.
%
% USAGE:  Metrics = bst_resolution_metrics(Kernel, Leadfield, GridLoc)
%
% INPUTS:
%   Kernel    : [nSrc x nCh]   imaging kernel (vertex space)
%   Leadfield : [nCh  x nSrc]  constrained leadfield (vertex space)
%   GridLoc   : [nSrc x 3]     source positions in meters
%
% OUTPUT (per-vertex, [nSrc x 1] unless noted):
%   .Res                resolution matrix [nSrc x nSrc] = Kernel*Leadfield
%   .LocError           localization error (mm): distance from true source to PSF peak
%   .SpatialDispersion  PSF spread (mm): sqrt(sum(d^2 * psf^2)/sum(psf^2))
%   .OverallAmplitude   peak PSF amplitude per source (depth-bias indicator)
%
% Each column j of Res is the point-spread function (PSF) for source j.
%
% Authors: Diellor Basha, 2026

Res = Kernel * Leadfield;                 % [nSrc x nSrc]
nSrc = size(Res, 2);
LocError = zeros(nSrc, 1);
SpatialDispersion = zeros(nSrc, 1);
OverallAmplitude  = zeros(nSrc, 1);
mm = 1000;                                % meters -> mm

for j = 1:nSrc
    psf = abs(Res(:, j));
    [pk, iPk] = max(psf);
    OverallAmplitude(j) = pk;
    % Localization error: true source location j vs PSF peak location
    LocError(j) = norm(GridLoc(iPk, :) - GridLoc(j, :)) * mm;
    % Spatial dispersion about the true source location
    d = sqrt(sum((GridLoc - GridLoc(j, :)).^2, 2)) * mm;   % [nSrc x 1]
    w = psf.^2;
    if sum(w) > 0
        SpatialDispersion(j) = sqrt(sum(d.^2 .* w) / sum(w));
    end
end

Metrics = struct('Res', Res, 'LocError', LocError, ...
    'SpatialDispersion', SpatialDispersion, 'OverallAmplitude', OverallAmplitude);
end
