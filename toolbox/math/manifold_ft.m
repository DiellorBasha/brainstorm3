function C = manifold_ft(Phi, M, U)
% MANIFOLD_FT: Forward manifold Fourier transform. Project a vertex field onto an
% M-orthonormal Laplace-Beltrami eigenmode basis: C = Phi' * (M * U).
%
% USAGE:  C = manifold_ft(Phi, M, U)
% INPUTS:
%   Phi : [nV x K]  eigenvectors (already selected, e.g. canonical order)
%   M   : [nV x nV] mass matrix (basis is M-orthonormal)
%   U   : [nV x nT] vertex field(s)
% OUTPUT:
%   C   : [K x nT]  mode coefficients
%
% Authors: Diellor Basha, 2026
Phi = double(Phi); U = double(U);
if size(U,1) ~= size(Phi,1)
    error('manifold_ft: U has %d rows but Phi has %d.', size(U,1), size(Phi,1));
end
C = Phi' * (M * U);
end
