function [V, D] = bst_eigs_smallest(A, B, k, opts)
% BST_EIGS_SMALLEST: k smallest generalized eigenpairs of a (near-)singular,
% symmetric/Hermitian pencil (A, B) with B SPD, avoiding the sigma=0 shift-invert
% of the singular A that triggers MATLAB's RCOND "close to singular" warning.
%
% USAGE:  [V, D] = bst_eigs_smallest(A, B, k, opts)
%
% INPUTS:
%    - A    : [n x n] symmetric (real) or Hermitian (complex) operator, may be singular
%    - B    : [n x n] symmetric positive-definite mass matrix
%    - k    : number of smallest-magnitude generalized eigenpairs to return
%    - opts : struct passed through to eigs (e.g. tol, maxit, disp)
% OUTPUTS:
%    - V, D : eigenvectors / diagonal eigenvalues, same meaning as eigs(A,B,k,'smallestabs')
%
% Forces the symmetric/Hermitian Lanczos path (real spectrum, faster, clean
% degenerate multiplets) and selects the smallest modes via a small negative
% sigma shift so (A - sigma*B) = A + |sigma|*B is SPD/well-conditioned.

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

    if nargin < 4 || isempty(opts)
        opts = struct('tol', 1e-6, 'maxit', 1000, 'disp', 0);
    end
    % Symmetrize / Hermitian-ize: no-op to ~1e-16 for already-symmetric operators,
    % but makes issymmetric()/ishermitian() true so eigs uses Lanczos, not Arnoldi.
    A = (A + A') / 2;
    B = (B + B') / 2;
    % Factorization-free spectrum-scale estimate: 'largestabs' factorizes only the
    % SPD, well-conditioned mass B (never the singular A), so it cannot warn.
    lmax = abs(eigs(A, B, 1, 'largestabs', opts));
    if ~isfinite(lmax) || (lmax <= 0)
        [V, D] = eigs(A, B, k, 'smallestabs', opts);   % degenerate scale: legacy path
        return;
    end
    % Small negative shift below the (PSD) spectrum bottom. sigma < 0 can never
    % coincide with an eigenvalue, and 'nearest sigma' still returns the k smallest
    % (including the near-zero kernel).
    try
        [V, D] = eigs(A, B, k, -1e-7 * lmax, opts);
    catch
        try
            [V, D] = eigs(A, B, k, -1e-4 * lmax, opts);   % larger lift
        catch
            [V, D] = eigs(A, B, k, 'smallestabs', opts);  % legacy fallback
        end
    end
end
