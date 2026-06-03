function U = manifold_ift(Phi, C)
% MANIFOLD_IFT: Inverse manifold Fourier transform. Synthesize a vertex field from
% mode coefficients (or a mode-space kernel): U = Phi * C.
%
% USAGE:  U = manifold_ift(Phi, C)
% INPUTS:
%   Phi : [nV x K]  eigenvectors (already selected)
%   C   : [K x nT]  coefficients, or [K x nCh] mode-space kernel
% OUTPUT:
%   U   : [nV x nT] or [nV x nCh]
%
% Authors: Diellor Basha, 2026
Phi = double(Phi); C = double(C);
if size(C,1) ~= size(Phi,2)
    error('manifold_ift: C has %d rows but Phi has %d columns.', size(C,1), size(Phi,2));
end
U = Phi * C;
end
