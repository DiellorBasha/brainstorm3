function test_view_eigenmodes_pure
% Verify the pure viewer helpers reached through the macro dispatch:
% GetModeDisplay (data column, symmetric color limits, eigenvalue/wavelength,
% label, index clamping) and StepMode (clamped stepping).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% ----- Fabricate a small Eigenmodes struct (as in_tess_eigenmodes returns) -----
rng(3);
nVert = 50; K = 8;
Eig = struct();
Eig.Vectors  = randn(nVert, K);
Eig.Values   = sort(abs(randn(K,1))) + 0.1;   % ascending, > 0
Eig.nModes   = K;
Eig.MassType = 'barycentric';

% ----- GetModeDisplay: basic correctness -----
d = view_eigenmodes('GetModeDisplay', Eig, 3);
assert(isequal(d.Data, Eig.Vectors(:,3)), 'Data must be the requested mode column.');
assert(d.iMode == 3 && d.nModes == K, 'iMode/nModes mismatch.');
m = max(abs(Eig.Vectors(:,3)));
assert(isequal(d.CLim, [-m, m]), 'CLim must be symmetric [-max|v|, +max|v|].');
assert(d.CLim(1) < d.CLim(2), 'CLim must be a non-degenerate increasing range.');
assert(abs(d.Lambda - Eig.Values(3)) < 1e-12, 'Lambda must be Values(iMode).');
assert(abs(d.Wavelength - 2*pi/sqrt(Eig.Values(3))) < 1e-9, 'Wavelength = 2*pi/sqrt(lambda).');
assert(ischar(d.Label) && ~isempty(d.Label), 'Label must be a non-empty string.');

% ----- GetModeDisplay: index clamping -----
dLo = view_eigenmodes('GetModeDisplay', Eig, -5);
assert(dLo.iMode == 1, 'iMode below 1 must clamp to 1.');
dHi = view_eigenmodes('GetModeDisplay', Eig, K+99);
assert(dHi.iMode == K, 'iMode above K must clamp to K.');

% ----- GetModeDisplay: degenerate (all-zero) mode guard -----
EigZ = Eig; EigZ.Vectors(:,2) = 0;
dz = view_eigenmodes('GetModeDisplay', EigZ, 2);
assert(dz.CLim(1) < dz.CLim(2), 'All-zero mode must still yield a non-degenerate CLim.');

% ----- GetModeDisplay: non-positive eigenvalue -> wavelength n/a -----
EigN = Eig; EigN.Values(1) = 0;
dn = view_eigenmodes('GetModeDisplay', EigN, 1);
assert(isnan(dn.Wavelength), 'lambda<=0 must give NaN wavelength.');
assert(~isempty(strfind(dn.Label, 'n/a')), 'Label must show n/a for lambda<=0.');

% ----- StepMode: stepping + clamping -----
assert(view_eigenmodes('StepMode', 3, +1, K) == 4, 'StepMode +1 failed.');
assert(view_eigenmodes('StepMode', 3, -1, K) == 2, 'StepMode -1 failed.');
assert(view_eigenmodes('StepMode', 1, -1, K) == 1, 'StepMode must clamp at 1.');
assert(view_eigenmodes('StepMode', K, +1, K) == K, 'StepMode must clamp at K.');
assert(view_eigenmodes('StepMode', K-2, +10, K) == K, 'StepMode +10 must clamp at K.');
assert(view_eigenmodes('StepMode', 3, -10, K) == 1, 'StepMode -10 must clamp at 1.');

disp('ALL TESTS PASSED');
end
