function [Results, errMsg] = bst_inverse_eigenmodes(varargin)
% BST_INVERSE_EIGENMODES: Mode-space MEG/EEG inverse on a composed eigenmode leadfield.
%
% USAGE:
%   [Results, errMsg] = bst_inverse_eigenmodes(CompHeadModelFile, NoiseCovFile, ChannelFile, GoodChannel, ...)
%   Kernel            = bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, Proj, Method, Prior, Alpha, SNR, Unreg)
%
% DESCRIPTION:
%   Consumes a composed eigenmode head model (Gain = L*Phi, with .Eigenvalues and
%   .nModes; produced by bst_eigenmode_leadfield) and solves the regularized MAP
%   estimate in mode space:
%       M~ = R L~' (L~ R L~' + lambda C)^-1        [K x nGoodChannels]
%   where R = bst_eigenmode_prior(lambdas, K, Prior, Alpha) is the spectral source
%   covariance (replaces depth weighting) and C is the noise covariance (applied as
%   the whitener iW). SSP projectors and bad channels are folded into the clean
%   operator iW*Proj exactly as the standard inverse does.
%
%   The math core is exposed as ('SolvePure', ...) for unit testing.
%
% OPTIONS (name-value):
%   'Method'    : 'mne' (default) | 'dspm' | 'sloreta'
%   'Prior'     : 'log' (default) | 'flat' | 'power'
%   'Alpha'     : power-prior exponent (default 1)
%   'SNR'       : signal-to-noise ratio for regularization (default 3)
%   'Unreg'     : logical; if true, ignore SNR and use rank-safe pinv (default false)
%   'nModes'    : cap on modes used (default: all in the composed model)
%   (no DataTypes option — whitener is computed directly on the subset covariance)
%
% Authors: Diellor Basha, 2026

% Dispatch the pure math entry
if ischar(varargin{1}) && strcmpi(varargin{1}, 'SolvePure')
    [L_tilde, lambdas, iW, Proj, Method, Prior, Alpha, SNR, Unreg] = varargin{2:10};
    K = size(L_tilde, 2);
    R = bst_eigenmode_prior(lambdas, K, Prior, Alpha);
    Results = solve_modespace(L_tilde, R, iW, Proj, Method, SNR, Unreg);
    errMsg = '';
    return;
end

Results = []; errMsg = '';
[CompHeadModelFile, NoiseCovFile, ChannelFile, GoodChannel] = varargin{1:4};

% Options
Method = 'mne'; Prior = 'log'; Alpha = 1; SNR = 3; Unreg = false; nModes = [];
for i = 5:2:numel(varargin)
    switch lower(varargin{i})
        case 'method',    Method    = lower(varargin{i+1});
        case 'prior',     Prior     = lower(varargin{i+1});
        case 'alpha',     Alpha     = varargin{i+1};
        case 'snr',       SNR       = varargin{i+1};
        case 'unreg',     Unreg     = logical(varargin{i+1});
        case 'nmodes',    nModes    = varargin{i+1};
    end
end
if ~ismember(Method, {'mne','dspm','sloreta'})
    errMsg = ['Unknown method: ' Method '. Use mne, dspm, or sloreta.']; return;
end

% Load composed head model
HM = in_bst_headmodel(CompHeadModelFile, 0);
if ~isfield(HM, 'isEigenmode') || ~HM.isEigenmode
    errMsg = 'Head model is not an eigenmode leadfield. Run "Compute eigenmode leadfield" first.'; return;
end
L_all = double(HM.Gain);                          % [nAllCh x K]
lambdas = double(HM.Eigenvalues(:));
nAllCh = size(L_all, 1);

% Good channels
if isempty(GoodChannel) || numel(GoodChannel) ~= nAllCh
    iGood = (1:nAllCh)';
else
    iGood = find(GoodChannel(:));
end
L_tilde = L_all(iGood, :);
K = size(L_tilde, 2);
if ~isempty(nModes) && nModes > 0 && nModes < K
    K = nModes; L_tilde = L_tilde(:, 1:K); lambdas = lambdas(1:K);
end
nCh = numel(iGood);

% Whitener iW = C^(-1/2) from the noise covariance, restricted to good channels.
% Computed directly (symmetric, regularized SVD) on the already-subset covariance
% so the dimensions are guaranteed consistent with L_tilde's good-channel rows.
if ~isempty(NoiseCovFile)
    NC = load(file_fullpath(NoiseCovFile));
    NoiseCov = double(NC.NoiseCov(iGood, iGood));
    NoiseCov = 0.5 * (NoiseCov + NoiseCov');        % symmetrize
    [Un, Sn] = svd(NoiseCov);
    sn = diag(Sn);
    reg = max(sn) * 1e-6;                            % Tikhonov floor on noise eigenvalues
    iW = Un * diag(1 ./ sqrt(sn + reg)) * Un';      % [nCh x nCh] symmetric whitener
else
    iW = eye(nCh);
end

% Projector from the channel file (SSP), restricted to good channels
Proj = eigenmode_projector(ChannelFile, iGood, nCh);

% Build prior and solve
R = bst_eigenmode_prior(lambdas, K, Prior, Alpha);
Kernel = solve_modespace(L_tilde, R, iW, Proj, Method, SNR, Unreg);

Results = struct();
Results.ImagingKernel = Kernel;        % [K x nGoodChannels]
Results.nModes        = K;
Results.Method        = Method;
Results.Prior         = Prior;
Results.SNR           = SNR;
Results.Unreg         = Unreg;
Results.GoodChannel   = iGood;
Results.Eigenvalues   = lambdas(1:K);
Results.SourcePrior   = R;
Results.SurfaceFile   = HM.SurfaceFile;
end

% ---- core solver (whitened, projected, mode-space MAP) ----
function Kernel = solve_modespace(L_tilde, R, iW, Proj, Method, SNR, Unreg)
% Clean operator: project then whiten. Folded into the final kernel so it maps RAW data.
Cop  = iW * Proj;                         % [nCh x nCh]
Lc   = Cop * L_tilde;                     % cleaned compressed leadfield [nCh x K]
nCh  = size(Lc, 1);
K    = size(Lc, 2);
sP   = sqrt(R(:));                        % column scaling by sqrt(prior)
Lws  = Lc .* sP.';                        % [nCh x K]

if Unreg
    % Rank-safe pseudoinverse (harmonic limit). Prior column scaling cancels.
    Kernel = pinv(Lc) * Cop;
    return;
end

[U, S, V] = svd(Lws, 'econ');
s = diag(S);
Lambda = sum(s.^2) / (nCh * SNR^2);
alpha  = s ./ (s.^2 + Lambda);            % MNE filter factors
Kmne_white = diag(sP) * V * diag(alpha) * U';   % [K x nCh], whitened-data side

switch Method
    case 'mne'
        Kernel = Kmne_white * Cop;
    case 'dspm'
        nn = sqrt(sum(Kmne_white.^2, 2));
        nn(nn < max(nn)*1e-10) = max(nn)*1e-10;
        Kernel = (Kmne_white ./ nn) * Cop;
    case 'sloreta'
        sa = s .* alpha;
        ResDiag = sum((V.^2) .* (sa'), 2);       % diag of resolution matrix
        ResDiag(ResDiag < max(ResDiag)*1e-10) = max(ResDiag)*1e-10;
        Kernel = (Kmne_white ./ sqrt(ResDiag)) * Cop;
end
end

% ---- SSP projector assembly from a channel file ----
function Proj = eigenmode_projector(ChannelFile, iGood, nCh)
Proj = eye(nCh);
if isempty(ChannelFile); return; end
try
    ChannelMat = in_bst_channel(ChannelFile);
catch
    return;
end
if ~isfield(ChannelMat, 'Projector') || isempty(ChannelMat.Projector); return; end
P = ChannelMat.Projector;
if isstruct(P)
    U = [];
    for k = 1:numel(P)
        if isfield(P(k),'Status') && P(k).Status == 1 && ~isempty(P(k).Components)
            comps = P(k).Components(iGood, :);
            if isfield(P(k),'CompMask') && ~isempty(P(k).CompMask)
                comps = comps(:, logical(P(k).CompMask));
            end
            U = [U, comps]; %#ok<AGROW>
        end
    end
    if ~isempty(U)
        [Uo, ~] = qr(U, 0);
        Proj = eye(nCh) - Uo*Uo';
    end
elseif isnumeric(P) && (size(P,1) == size(P,2))
    if size(P,1) == nCh
        Proj = P;                    % already restricted to good channels
    elseif size(P,1) >= max(iGood)
        Proj = P(iGood, iGood);      % full-channel projector matrix -> subset
    end
end
end
