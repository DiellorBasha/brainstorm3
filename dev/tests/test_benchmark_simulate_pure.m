function test_benchmark_simulate_pure
% Verify forward projection + colored-noise addition at a target SNR.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fullfile(repoRoot,'dev','benchmarks'));
if ~brainstorm('status'); brainstorm nogui; end

nCh = 20; nV = 8; nTime = 50;
rng(0);
L = randn(nCh, nV);                          % leadfield
Sources = zeros(nV, nTime); Sources(3,:) = sin(2*pi*(1:nTime)/nTime);  % one active vertex
C = eye(nCh);                                % white noise cov (for a clean SNR check)

snrTarget = 6;  % dB
Sim = bst_benchmark_simulate(L, Sources, C, 'SNR', snrTarget, 'Seed', 1);
assert(isequal(size(Sim.F), [nCh nTime]), 'F must be [nCh x nTime].');

% Achieved SNR (signal vs noise power) should match the target within tolerance
sigPow = mean(Sim.Fsignal(:).^2);
noiPow = mean(Sim.Fnoise(:).^2);
snrAch = 10*log10(sigPow / noiPow);
assert(abs(snrAch - snrTarget) < 0.5, 'achieved SNR must match target within 0.5 dB.');

% Reproducible
Sim2 = bst_benchmark_simulate(L, Sources, C, 'SNR', snrTarget, 'Seed', 1);
assert(max(abs(Sim.F(:) - Sim2.F(:))) < 1e-12, 'same Seed -> identical noise draw.');

disp('ALL TESTS PASSED');
end
