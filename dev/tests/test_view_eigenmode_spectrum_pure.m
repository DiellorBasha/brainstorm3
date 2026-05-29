function test_view_eigenmode_spectrum_pure
% Verify the pure render helpers of view_eigenmode_spectrum.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
addpath(fullfile(repoRoot, 'toolbox', 'gui'));

% ---- ComputeModalPower: |theta|^2 split by hemisphere component ----
Theta     = [3+4i; 1; 0; -2];          % use complex to prove magnitude-squared
Component = [1; 1; 2; 2];
pw = view_eigenmode_spectrum('ComputeModalPower', Theta, Component);
assert(isequal(size(pw.left),  [2 1]), 'left power wrong size');
assert(isequal(size(pw.right), [2 1]), 'right power wrong size');
assert(abs(pw.left(1)  - 25) < 1e-12, '|3+4i|^2 should be 25');   % 3^2+4^2
assert(abs(pw.left(2)  - 1)  < 1e-12, 'left(2) should be 1');
assert(abs(pw.right(1) - 0)  < 1e-12, 'right(1) should be 0');
assert(abs(pw.right(2) - 4)  < 1e-12, 'right(2) should be 4');    % (-2)^2

% ---- GetSpectrumAxis: eigenvalue passthrough, wavelength = 2pi/sqrt(lambda) ----
Values = [0; 1; (2*pi)^2];             % lambda=0 -> n/a; lambda=(2pi)^2 -> wavelength 1
axE = view_eigenmode_spectrum('GetSpectrumAxis', Values, 'eigenvalue');
assert(isequal(axE.x, Values), 'eigenvalue axis must pass values through');
assert(ischar(axE.label) && ~isempty(axE.label), 'eigenvalue label missing');
axW = view_eigenmode_spectrum('GetSpectrumAxis', Values, 'wavelength');
assert(isnan(axW.x(1)), 'lambda<=0 wavelength must be NaN');
assert(abs(axW.x(2) - 2*pi) < 1e-12, 'wavelength of lambda=1 is 2pi');
assert(abs(axW.x(3) - 1)    < 1e-12, 'wavelength of lambda=(2pi)^2 is 1');

% ---- GetWindowAverage: mean power over a sample window ----
T = [1 2 3; 0 0 6];                    % [K x nTime], real
avgFull = view_eigenmode_spectrum('GetWindowAverage', T, []);
assert(abs(avgFull(1) - mean([1 4 9]))  < 1e-12, 'row1 full mean wrong');
assert(abs(avgFull(2) - mean([0 0 36])) < 1e-12, 'row2 full mean wrong');
avgWin = view_eigenmode_spectrum('GetWindowAverage', T, [1 2]);
assert(abs(avgWin(1) - mean([1 4]))  < 1e-12, 'row1 window mean wrong');
assert(abs(avgWin(2) - mean([0 0]))  < 1e-12, 'row2 window mean wrong');

fprintf('ALL TESTS PASSED: test_view_eigenmode_spectrum_pure\n');
end
