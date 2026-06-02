function ImagingKernel = bst_eigenmode_reconstruct(SurfaceOrPhi, ModeKernel)
% BST_EIGENMODE_RECONSTRUCT: Reconstruct a vertex-space imaging kernel from a
% mode-space kernel: ImagingKernel = Phi(:,1:K) * ModeKernel.
%
% USAGE:  ImagingKernel = bst_eigenmode_reconstruct(SurfaceFile, ModeKernel)
%         ImagingKernel = bst_eigenmode_reconstruct(Phi,         ModeKernel)
%
%   SurfaceOrPhi : cortex SurfaceFile (eigenmodes loaded via in_tess_eigenmodes)
%                  OR a numeric eigenvector matrix Phi [nVert x nModes]
%   ModeKernel   : [K x nGoodCh] mode-space kernel (from bst_inverse_eigenmodes)
%   ImagingKernel: [nVert x nGoodCh] cortex kernel
%
% Authors: Diellor Basha, 2026

K = size(ModeKernel, 1);
if ischar(SurfaceOrPhi)
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceOrPhi);
    if ~isComputed
        error(['No eigenmodes on surface: ' SurfaceOrPhi '. Run "Compute eigenmodes" first.']);
    end
    Phi = double(Eig.Vectors);
else
    Phi = double(SurfaceOrPhi);
end
if size(Phi, 2) < K
    error('Surface has fewer eigenmodes (%d) than kernel modes (%d).', size(Phi,2), K);
end
ImagingKernel = Phi(:, 1:K) * ModeKernel;
end
