function [Ffilt, Messages, isError] = bst_eigenfilter(F, EigenMat, OperatorMat, KernelName, KernelParams)
% BST_EIGENFILTER: Spatially filter a surface-mapped source field in an operator eigenbasis.
%
% USAGE:  [Ffilt, Messages, isError] = bst_eigenfilter(F, EigenMat, OperatorMat, KernelName, KernelParams)
%
% DESCRIPTION:
%     The spatial-domain analogue of a bandpass filter (bst_bandpass_*), applied in the
%     eigenbasis. Per hemisphere it projects the field onto the modes, scales each mode by a
%     gain h(lambda) from the eigfilter library, and reconstructs:
%
%         C = Phi' * (B * U);   C <- h(lambda) .* C;   U_filt = Phi * C
%           = manifold_ift(Phi, h .* manifold_ft(Phi, B, U))
%
%     General over the operator family carried in the eigen_ node (EigenMat.Variant):
%     Laplace-Beltrami (scalar), Connection Laplacian (complex tangent, operator frame), and
%     Dirac (3D vector embedded as a pure-imaginary quaternion). The gain acts on the K mode
%     coefficients, so it is layout-agnostic; only the field's row layout differs per variant,
%     mapped here directly from EigenMat.GlobalVertices / GlobalFaces. LH and RH are filtered
%     SEPARATELY (the basis is B-orthonormal per hemisphere). Designed to be driven by
%     bst_eigen (the 'filter' method), which loads the files and passes them in.
%
% INPUT:
%     - F            : [nSource x nTime] source map (scalar [nV] / complex [nV] / 3-vector [3nV])
%     - EigenMat     : eigen_ node (in_bst_eigen): .Variant, .Phi{h}, .Lambda{h},
%                      .GlobalVertices{h}, .GlobalFaces{h}
%     - OperatorMat  : operator_ node (in_bst_operator): .Mass{h} (per-hemisphere mass B)
%     - KernelName   : eigfilter kernel name ('heat','tikhonov','mexhat',...) | function handle
%                      | [K x 1] precomputed per-mode gain
%     - KernelParams : struct of kernel parameters (default: struct())
%
% OUTPUT:
%     - Ffilt    : [nSource x nTime] filtered source map (same layout as F)
%     - Messages : '' on success, else an error string
%     - isError  : 1 if an error happened
%
% SEE ALSO: bst_eigen, bst_eigenspectrum, bst_eigfilter_kernel, manifold_ft, manifold_ift

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

if (nargin < 5) || isempty(KernelParams); KernelParams = struct(); end
Messages = '';
isError  = 0;
Ffilt    = zeros(size(F));
nT       = size(F, 2);

for h = 1:numel(EigenMat.Phi)
    Phi = EigenMat.Phi{h};
    if isempty(Phi); continue; end
    Lam = EigenMat.Lambda{h}(:);
    B   = OperatorMat.Mass{h};

    % --- Row map from the source-map layout to this hemisphere's basis rows (by variant) ---
    % srcRows : rows of F (and Ffilt) for this hemisphere;  dstRows : matching basis rows;
    % nrows   : basis dimension (nVh, or 4*nVh for the quaternion Dirac block layout).
    switch EigenMat.Variant
        case {'Laplace-Beltrami', 'Connection Laplacian'}
            gv = EigenMat.GlobalVertices{h}(:);
            if max(gv) > size(F, 1)
                Messages = sprintf('bst_eigenfilter: hemisphere %d indexes row %d but the map has %d rows.', h, max(gv), size(F,1));
                isError = 1; return;
            end
            if strcmp(EigenMat.Variant, 'Connection Laplacian') && isreal(F)
                Messages = ['bst_eigenfilter: Connection Laplacian needs a COMPLEX tangent field in the ' ...
                            'operator frame (the real-3D embedding needs the not-yet-persisted per-vertex frame).'];
                isError = 1; return;
            end
            srcRows = gv;  dstRows = (1:numel(gv))';  nrows = numel(gv);

        case {'Dirac', 'Dirac-Face', 'Hodge-Face'}
            if strcmp(EigenMat.Variant, 'Dirac')
                idx = EigenMat.GlobalVertices{h}(:);
            else
                idx = EigenMat.GlobalFaces{h}(:);
            end
            if isempty(idx) || (3*max(idx) > size(F, 1))
                Messages = sprintf(['bst_eigenfilter: %s hemisphere %d needs a 3-vector map [3*N x nTime] ' ...
                    '(row %d) but the map has %d rows.'], EigenMat.Variant, h, 3*max([idx;0]), size(F,1));
                isError = 1; return;
            end
            n = numel(idx);
            srcRows = reshape([(idx-1)*3+1, (idx-1)*3+2, (idx-1)*3+3].', [], 1);     % x1;y1;z1;x2;...
            dstRows = reshape([(0:n-1)*4+2; (0:n-1)*4+3; (0:n-1)*4+4], [], 1);       % quat rows 2;3;4;6;...
            nrows   = 4*n;                                                            % w rows (1:4:end) stay 0

        otherwise
            Messages = sprintf('bst_eigenfilter: unsupported eigen variant ''%s''.', EigenMat.Variant);
            isError = 1; return;
    end

    % --- Per-mode gain h(lambda) from the eigfilter library (re-resolved per hemisphere) ---
    if isa(KernelName, 'function_handle')
        hgain = bst_eigfilter_evaluate(KernelName, Lam);
    elseif isnumeric(KernelName)
        hgain = KernelName(:);
    else
        g = bst_eigfilter_kernel(KernelName, KernelParams);
        if iscell(g)
            Messages = 'bst_eigenfilter: kernel returned a filterbank (vector param); pass a single-scale kernel.';
            isError = 1; return;
        end
        hgain = bst_eigfilter_evaluate(g, Lam);
    end
    if numel(hgain) ~= numel(Lam)
        Messages = sprintf('bst_eigenfilter: gain has %d entries but the basis has %d modes.', numel(hgain), numel(Lam));
        isError = 1; return;
    end

    % --- Embed -> project -> gain -> reconstruct -> de-embed ---
    U = zeros(nrows, nT);
    U(dstRows, :)  = F(srcRows, :);                       % embed (identity for scalar/complex)
    C  = manifold_ft(Phi, B, U);                          % [K x nT]
    Uf = manifold_ift(Phi, hgain .* C);                   % [nrows x nT]
    Ffilt(srcRows, :) = Uf(dstRows, :);                   % de-embed back into the source-map layout
end
end
