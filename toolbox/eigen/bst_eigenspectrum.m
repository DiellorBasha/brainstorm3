function [TF, kAxis, Nwin, Messages, TFbis] = bst_eigenspectrum(U, Phi, B, Lambda, WinLength, WinOverlap, WinFunc, Measure)
% BST_EIGENSPECTRUM: Windowed eigenspectrum of a surface-mapped field. The spatial-domain
%                    analogue of BST_PSD (Welch): replace the temporal FFT by the forward
%                    manifold transform onto an operator eigenbasis.
%
% USAGE:  [TF, kAxis, Nwin, Messages, TFbis] = bst_eigenspectrum(U, Phi, B, Lambda, WinLength, WinOverlap, WinFunc, Measure)
%
% DESCRIPTION:
%     For a field U on a surface (or hemisphere submesh), project each time sample onto the
%     B-orthonormal eigenbasis Phi -- C = Phi' * (B * U) = manifold_ft(Phi, B, U) -- and
%     aggregate the per-mode power (or magnitude) over time windows, exactly as bst_psd
%     windows the time series and averages the FFT power. The spectral axis is the mode
%     axis kAxis = sqrt(Lambda) (spatial frequency [1/m]) instead of the Hz axis.
%
%     HEMISPHERE-AGNOSTIC, by design. Pass ONE hemisphere's (U, Phi, B, Lambda); the caller
%     (bst_eigen) computes LH and RH with SEPARATE calls. Computations are never combined
%     across hemispheres (the basis is B-orthonormal per hemisphere). This mirrors how
%     bst_psd is channel-agnostic and processes whatever F it is handed.
%
%     This is general over the operator family (Laplace-Beltrami scalar, Connection
%     Laplacian complex-tangent, Dirac 3D-vector): the transform C = Phi'*(B*U) is identical
%     in all cases -- only the row layout of U/Phi/B differs (nV, complex nV, or 3nV), which
%     the caller supplies consistently. MATLAB's ' is the conjugate transpose, so the
%     complex (connection) Hermitian inner product is handled automatically.
%
% INPUTS:
%     - U          : [nV x nT]   field on this hemisphere's vertices/faces (the analogue of F)
%     - Phi        : [nV x K]    B-orthonormal eigenvectors (this hemisphere)
%     - B          : [nV x nV]   mass matrix (this hemisphere; basis is B-orthonormal)
%     - Lambda     : [K x 1]     eigenvalues (this hemisphere); kAxis = sqrt(Lambda)
%     - WinLength  : window length in SAMPLES (time steps); [] or 0 => one window (whole series)
%     - WinOverlap : window overlap in percent (default: 50)
%     - WinFunc    : aggregation across windows: 'mean' | 'std' | 'mean+std' (default: 'mean')
%     - Measure    : per-sample measure before averaging: 'power' | 'magnitude' (default: 'power')
%
% OUTPUTS:
%     - TF       : [K x 1]   eigenspectrum (mode power/magnitude, averaged over windows)
%     - kAxis    : [K x 1]   spatial-frequency axis = sqrt(Lambda)
%     - Nwin     : number of windows used
%     - Messages : string, warnings/errors ('' if none)
%     - TFbis    : [K x 1]   across-window standard deviation when WinFunc includes 'std' (else [])
%
% SEE ALSO: bst_eigen, manifold_ft, manifold_ift, bst_psd
%
% Authors: Diellor Basha, 2026

% ===== PARSE INPUTS =====
if (nargin < 8) || isempty(Measure),    Measure    = 'power'; end
if (nargin < 7) || isempty(WinFunc),    WinFunc    = 'mean';  end
if (nargin < 6) || isempty(WinOverlap), WinOverlap = 50;      end
if (nargin < 5),                        WinLength  = [];      end

TF       = [];
TFbis    = [];
Nwin     = [];
Messages = '';
kAxis    = sqrt(double(Lambda(:)));

% ===== VALIDATE SIZES =====
[nV, nT] = size(U);
K        = size(Phi, 2);
if size(Phi, 1) ~= nV
    Messages = sprintf('bst_eigenspectrum: U has %d rows but Phi has %d.', nV, size(Phi,1));
    return;
end
if (size(B,1) ~= nV) || (size(B,2) ~= nV)
    Messages = sprintf('bst_eigenspectrum: B must be [%d x %d], got [%d x %d].', nV, nV, size(B,1), size(B,2));
    return;
end
if numel(Lambda) ~= K
    Messages = sprintf('bst_eigenspectrum: Lambda has %d entries but Phi has %d modes.', numel(Lambda), K);
    return;
end
computeStd = ~isempty(strfind(lower(WinFunc), 'std')); %#ok<STREMP>

% ===== WINDOWING (over time) =====
% Mirror bst_psd's window tiling, but in samples (the engine is time-axis-agnostic; bst_eigen
% converts a seconds-based window to samples using the recording sfreq before calling).
if isempty(WinLength) || (WinLength == 0) || (WinLength >= nT)
    Lwin     = nT;
    Loverlap = 0;
    Nwin     = 1;
else
    Lwin     = round(WinLength);
    Loverlap = round(Lwin * WinOverlap / 100);
    Nwin     = floor((nT - Loverlap) ./ (Lwin - Loverlap));
end
if (Nwin < 1)
    Nwin = 1; Lwin = nT; Loverlap = 0;
end

% ===== ACCUMULATE PER-MODE MEASURE ACROSS WINDOWS =====
S1 = zeros(K, 1);                       % sum of per-window spectra
if computeStd, S2 = zeros(K, 1); end    % sum of squares (for std)
for iWin = 1:Nwin
    iTimes = (1:Lwin) + (iWin-1)*(Lwin - Loverlap);
    iTimes = iTimes(iTimes <= nT);
    % Forward manifold transform of this window: C = Phi' * (B * U) -> [K x numel(iTimes)]
    C = manifold_ft(Phi, B, U(:, iTimes));
    % Per-sample measure, then average over the window's samples -> [K x 1]
    switch lower(Measure)
        case 'power',     m = abs(C).^2;
        case 'magnitude', m = abs(C);
        otherwise
            Messages = sprintf('bst_eigenspectrum: unknown Measure ''%s'' (use ''power'' or ''magnitude'').', Measure);
            return;
    end
    winSpec = mean(m, 2);
    S1 = S1 + winSpec;
    if computeStd, S2 = S2 + winSpec.^2; end
end

% ===== AGGREGATE ACROSS WINDOWS (mean / std) =====
TFmean = S1 ./ Nwin;
if computeStd
    TFstd = sqrt(max(S2 ./ Nwin - TFmean.^2, 0));
end
switch lower(WinFunc)
    case 'mean',     TF = TFmean; TFbis = [];
    case 'std',      TF = TFstd;  TFbis = [];
    case 'mean+std', TF = TFmean; TFbis = TFstd;
    otherwise
        Messages = sprintf('bst_eigenspectrum: unknown WinFunc ''%s''.', WinFunc);
        return;
end

if isempty(Messages)
    Messages = sprintf('Eigenspectrum: %d modes, %d window(s) of %d samples', K, Nwin, Lwin);
end
end
