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

% ---- BuildModeSpectrum: [2 x K] L/R matrix + shared x-vector per axis mode ----
% Two components (L,R), 3 ranks each. Theta values chosen so power = index is obvious.
Info = struct('Component', [1;1;1;2;2;2], 'CompRank', [1;2;3;1;2;3], ...
              'Values',    [1;4;9;1;4;9]);           % lambda = k^2 -> wavelength 2pi/k
Theta = [1; 2; 3; 4; 5; 6];                          % power = Theta.^2
spec = view_eigenmode_spectrum('BuildModeSpectrum', Info, Theta, 'index');
assert(isequal(size(spec.F), [2 3]), 'F must be [2 x K]');
assert(isequal(spec.F(1,:), [1 4 9]),   'Left power per rank wrong');     % 1^2,2^2,3^2
assert(isequal(spec.F(2,:), [16 25 36]), 'Right power per rank wrong');   % 4^2,5^2,6^2
assert(isequal(spec.x(:)', [1 2 3]), 'index x must be 1..K');
specE = view_eigenmode_spectrum('BuildModeSpectrum', Info, Theta, 'eigenvalue');
assert(isequal(specE.x(:)', [1 4 9]), 'eigenvalue x must be per-rank lambda');
specW = view_eigenmode_spectrum('BuildModeSpectrum', Info, Theta, 'wavelength');
assert(issorted(specW.x), 'x-vector must be ascending (figure_timeseries requirement)');
assert(max(abs(sort(specW.x(:)') - sort(2*pi./sqrt([1 4 9])))) < 1e-12, 'wavelength x values wrong');
% power modes change the y-quantity + label
specMag = view_eigenmode_spectrum('BuildModeSpectrum', Info, Theta, 'index', 'magnitude');
assert(isequal(specMag.F(1,:), [1 2 3]), 'magnitude = |theta| per rank');   % |1|,|2|,|3|
assert(~isempty(strfind(specMag.ylabel, 'magnitude')), 'magnitude ylabel');
specLog = view_eigenmode_spectrum('BuildModeSpectrum', Info, Theta, 'index', 'log');
assert(max(abs(specLog.F(1,:) - 10*log10([1 4 9]))) < 1e-12, 'log = 10*log10(power)');

fprintf('ALL TESTS PASSED: test_view_eigenmode_spectrum_pure\n');
end
