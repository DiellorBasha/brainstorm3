function test_resolution_metrics_pure
% Verify resolution-matrix point-spread metrics on a synthetic kernel/leadfield:
%   - identity resolution -> zero localization error, minimal dispersion
%   - shapes; finite; localization error in mm
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

% 4 vertices on a line 10mm apart; perfect kernel (Res = I) -> LE = 0
GridLoc = [0 0 0; 10 0 0; 20 0 0; 30 0 0] / 1000;   % meters
L = magic(4); L = L + 4*eye(4);                      % invertible leadfield [nCh=4 x nV=4]
Kern = inv(L);                                        % perfect inverse -> Res = I

M = bst_resolution_metrics(Kern, L, GridLoc);
assert(isfield(M,'LocError') && isfield(M,'SpatialDispersion'), 'Must return LE and SD.');
assert(numel(M.LocError) == 4, 'LE must be per-vertex.');
assert(max(M.LocError) < 1e-9, 'Perfect inverse must give zero localization error.');
assert(all(isfinite(M.SpatialDispersion)), 'SD must be finite.');

% A blurred kernel (average of neighbors) increases dispersion and LE
Blur = [0.5 0.5 0 0; 0.3 0.4 0.3 0; 0 0.3 0.4 0.3; 0 0 0.5 0.5];
Mb = bst_resolution_metrics(Blur * inv(L) * L, L, GridLoc);  %#ok<MINV>
assert(mean(Mb.SpatialDispersion) >= mean(M.SpatialDispersion), 'Blur must not reduce dispersion.');
disp('ALL TESTS PASSED');
end
