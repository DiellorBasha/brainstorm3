function Chat = manifold_jft(Phi, M, U, NFFT)
% MANIFOLD_JFT: Forward JOINT (time-vertex) Fourier transform. Projects a time-vertex field onto the
% cortex eigenbasis (manifold Fourier transform) AND the temporal Fourier basis, giving the joint
% spectral coefficients on the (eigenvalue lambda, temporal-frequency omega) grid:
%     Chat = fft( Phi' * (M * U), NFFT, 2 ) / sqrt(NFFT)
% This is the time extension of MANIFOLD_FT (cortex axis) and the analog of GSPBox gsp_jft. The
% 1/sqrt(NFFT) normalization makes the temporal transform unitary (Parseval-preserving), so the joint
% energy equals the cortex (M-weighted) energy for a field in the eigenspan.
%
% USAGE:  Chat = manifold_jft(Phi, M, U)        % NFFT = size(U,2)
%         Chat = manifold_jft(Phi, M, U, NFFT)  % zero-padded to NFFT temporal bins
% INPUTS:
%   Phi  : [nV x K]   eigenvectors (M-orthonormal, already selected)
%   M    : [nV x nV]  mass matrix (basis is M-orthonormal)
%   U    : [nV x nT]  time-vertex field (cortex x time)
%   NFFT : scalar     number of temporal FFT bins (default size(U,2))
% OUTPUT:
%   Chat : [K x NFFT] joint spectral coefficients (mode k, temporal frequency bin)
%
% SEE ALSO: manifold_ft, manifold_ijft
%
% Authors: Diellor Basha, 2026
Phi = double(Phi); U = double(U);
if size(U,1) ~= size(Phi,1)
    error('manifold_jft: U has %d rows but Phi has %d.', size(U,1), size(Phi,1));
end
if (nargin < 4) || isempty(NFFT), NFFT = size(U,2); end
C    = Phi' * (M * U);                      % manifold (cortex) Fourier transform  [K x nT]
Chat = fft(C, NFFT, 2) / sqrt(NFFT);        % temporal Fourier transform (unitary)  [K x NFFT]
end
