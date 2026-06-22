function ImagingKernel = bst_eigenmode_reconstruct(SurfaceOrPhi, ModeKernel, ModeIndices)
% BST_EIGENMODE_RECONSTRUCT: Reconstruct a vertex-space imaging kernel from a
% mode-space kernel: ImagingKernel = Phi(:,idx) * ModeKernel.
%
% USAGE:  ImagingKernel = bst_eigenmode_reconstruct(SurfaceFile, ModeKernel)
%         ImagingKernel = bst_eigenmode_reconstruct(Phi,         ModeKernel)
%         ImagingKernel = bst_eigenmode_reconstruct(..., ModeKernel, ModeIndices)
%
%   SurfaceOrPhi : cortex SurfaceFile (eigenmodes loaded via in_tess_eigenmodes)
%                  OR a numeric eigenvector matrix Phi [nVert x nModes]
%   ModeKernel   : [K x nGoodCh] mode-space kernel (from bst_inverse_eigenmodes)
%   ModeIndices  : (optional) columns of Phi the leadfield used, in solve order
%                  (headmodel.ModeIndices). REQUIRED for correctness when the
%                  leadfield selected modes across hemispheres (global eigenvalue
%                  order) rather than the first K columns. If omitted, falls back
%                  to Phi(:,1:K) for backward compatibility with old head models.
%   ImagingKernel: [nVert x nGoodCh] cortex kernel
%
% Authors: Diellor Basha, 2026

K = size(ModeKernel, 1);
idx = [];   % resolved below
if ischar(SurfaceOrPhi)
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceOrPhi);
    if ~isComputed
        error(['No eigenmodes on surface: ' SurfaceOrPhi '. Run "Compute eigenmodes" first.']);
    end
    Phi = double(Eig.Vectors);
    % Default to the surface's canonical order when no explicit selection is given.
    if (nargin < 3 || isempty(ModeIndices))
        idx = Eig.Order(:);
    else
        idx = ModeIndices(:);
    end
else
    Phi = double(SurfaceOrPhi);
    if (nargin >= 3) && ~isempty(ModeIndices)
        idx = ModeIndices(:);
    else
        idx = (1:K)';   % bare Phi matrix, no order known: positional fallback
    end
end
if numel(idx) < K
    error('Fewer mode indices (%d) than kernel modes (%d).', numel(idx), K);
end
idx = idx(1:K);
if max(idx) > size(Phi, 2)
    error('Mode index %d exceeds available eigenmodes (%d).', max(idx), size(Phi,2));
end
ImagingKernel = manifold_ift(Phi(:, idx), ModeKernel);   % [nV x nGoodCh]
end
