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
%    - V, D : eigenvectors / diagonal eigenvalues, ascending and real (the pencil is
%             symmetric/Hermitian), matching eigs(A,B,k,'smallestabs'). V may be
%             complex for a Hermitian A.
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
    % Factorization-free spectrum-scale estimate. Only an order of magnitude is needed
    % (to place the shift), so use a coarse tolerance, not the caller's tight opts.
    % 'largestabs' factorizes only the SPD, well-conditioned mass B (never the singular
    % A), so it cannot warn. Its own try/catch makes a non-converged estimate fall
    % through to the legacy path instead of hard-erroring before the shift fallbacks.
    estOpts = struct('tol', 1e-2, 'maxit', 300, 'disp', 0);
    try
        lmax = abs(eigs(A, B, 1, 'largestabs', estOpts));
    catch
        lmax = 0;
    end
    if ~isfinite(lmax) || (lmax <= 0)
        [V, D] = eigs(A, B, k, 'smallestabs', opts);     % degenerate estimate: legacy path
    else
        % Small negative shift below the (PSD) spectrum bottom: (A - sigma*B) =
        % A + |sigma|*B is SPD/well-conditioned, yet 'nearest sigma' still returns the
        % k smallest modes (including the near-zero kernel). The -1e-7 -> -1e-4
        % escalation lifts the factorization further if the first shift lands too close
        % to the kernel for eigs to converge; 'smallestabs' is the last-resort path.
        try
            [V, D] = eigs(A, B, k, -1e-7 * lmax, opts);
        catch
            try
                [V, D] = eigs(A, B, k, -1e-4 * lmax, opts);
            catch
                [V, D] = eigs(A, B, k, 'smallestabs', opts);
            end
        end
    end
    % Enforce the documented contract: real eigenvalues (symmetric/Hermitian pencil) in
    % ascending order, matching eigs(...,'smallestabs'). Eigenvectors are NOT realified
    % -- they are genuinely complex for a Hermitian A (e.g. the connection Laplacian).
    d = real(diag(D));
    [d, idx] = sort(d, 'ascend');
    V = V(:, idx);
    D = diag(d);
end
