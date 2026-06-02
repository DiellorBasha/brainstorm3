function Sim = bst_benchmark_simulate(L, Sources, NoiseCov, varargin)
% BST_BENCHMARK_SIMULATE: Forward-project GT sources and add colored sensor noise.
%
% USAGE:  Sim = bst_benchmark_simulate(L, Sources, NoiseCov, 'SNR', dB, 'Seed', s)
%
% INPUTS:
%   L        [nCh x nV]    constrained good-channel leadfield
%   Sources  [nV x nTime]  ground-truth source matrix
%   NoiseCov [nCh x nCh]   sensor noise covariance (colored noise model)
% OPTIONS:
%   'SNR'  : target sensor SNR in dB (default 6)
%   'Seed' : RNG seed (default 1)
%
% OUTPUT struct Sim: .F (signal+noise) .Fsignal .Fnoise  [all nCh x nTime], .SNR
%
% Colored noise uses the eigendecomposition model of process_simulate_recordings:
% xn = V*sqrt(D)*randn, scaled so 10*log10(sigPow/noisePow) = SNR.
%
% Authors: Diellor Basha, 2026
SNR = 6; Seed = 1;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'snr',  SNR  = varargin{i+1};
        case 'seed', Seed = varargin{i+1};
    end
end
L = double(L); Sources = double(Sources); NoiseCov = double(NoiseCov);
[nCh, nTime] = deal(size(L,1), size(Sources,2));

Fsignal = L * Sources;                              % [nCh x nTime]

% Colored noise from the covariance (real symmetric part for safety)
NoiseCov = 0.5*(NoiseCov + NoiseCov');
[Vn, Dn] = eig(NoiseCov);
dn = max(real(diag(Dn)), 0);
rng(Seed);
xn = real(Vn) * diag(sqrt(dn)) * randn(nCh, nTime);

% Scale noise to the target SNR
sigPow = mean(Fsignal(:).^2);
curPow = mean(xn(:).^2);
if curPow == 0; curPow = eps; end
noisePowTarget = sigPow / (10^(SNR/10));
xn = xn * sqrt(noisePowTarget / curPow);

Sim = struct('F', Fsignal + xn, 'Fsignal', Fsignal, 'Fnoise', xn, 'SNR', SNR);
end
