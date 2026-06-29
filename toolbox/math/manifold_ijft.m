function U = manifold_ijft(Phi, Chat, nKeep)
% MANIFOLD_IJFT: Inverse JOINT (time-vertex) Fourier transform. Synthesizes a time-vertex field from
% joint spectral coefficients (or a joint-domain kernel response):
%     U = Phi * ( ifft(Chat, [], 2) * sqrt(NFFT) )
% The time extension of MANIFOLD_IFT (cortex axis) and the analog of GSPBox gsp_ijft. The sqrt(NFFT)
% factor inverts the unitary 1/sqrt(NFFT) of MANIFOLD_JFT, so manifold_ijft(Phi, manifold_jft(Phi,M,U))
% recovers U exactly for a field in the eigenspan. Returns the (generally complex) synthesis; the caller
% takes real() for a real field (the imaginary part is ~0 up to numerical error / a complex kernel).
%
% USAGE:  U = manifold_ijft(Phi, Chat)         % keep all NFFT temporal bins
%         U = manifold_ijft(Phi, Chat, nKeep)  % truncate to the first nKeep time samples (un-pad)
% INPUTS:
%   Phi   : [nV x K]    eigenvectors (already selected)
%   Chat  : [K x NFFT]  joint spectral coefficients (or filtered coefficients / mode-space kernel)
%   nKeep : scalar      keep the first nKeep time samples (default NFFT; use nT to undo zero-padding)
% OUTPUT:
%   U     : [nV x nKeep] time-vertex field
%
% SEE ALSO: manifold_jft, manifold_ift, manifold_ijft
%
% Authors: Diellor Basha, 2026
Phi = double(Phi); Chat = double(Chat);
if size(Chat,1) ~= size(Phi,2)
    error('manifold_ijft: Chat has %d rows but Phi has %d columns.', size(Chat,1), size(Phi,2));
end
NFFT = size(Chat,2);
if (nargin < 3) || isempty(nKeep), nKeep = NFFT; end
C = ifft(Chat, [], 2) * sqrt(NFFT);         % inverse temporal Fourier transform  [K x NFFT]
C = C(:, 1:min(nKeep, NFFT));               % un-pad to the active time window
U = Phi * C;                                % inverse manifold (cortex) Fourier transform  [nV x nKeep]
end
