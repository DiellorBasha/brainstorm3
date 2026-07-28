function [Phi, Psi] = bst_helmholtz_potential(X, Cov)
% BST_HELMHOLTZ_POTENTIAL: Vectorized Helmholtz scalar potential Phi + stream function Psi.
%
% USAGE:  [Phi, Psi] = bst_helmholtz_potential(X, Cov)
%
% DESCRIPTION:
%     The Phi/Psi half of process_helmholtz('Compute'), solved for ALL columns of X at once
%     instead of one frame at a time. This lets the Helmholtz potentials be FUSED into a source
%     imaging kernel (X = ImagingKernel, columns = channels) as well as applied to a materialized
%     field (X = ImageGridAmp, columns = time). Phi solves the per-hemisphere Poisson problem
%     K*Phi = divergence(X), Psi solves K*Psi = vorticity(X); both are mean-centered per
%     hemisphere. Uses the cached pinned Cholesky factor (tess_cholesky) with a matrix RHS. div
%     and curl have no such engine -- call bst_divergence / bst_curl directly (they already
%     accept a multi-column field).
%
% INPUTS:
%     X   - [3nV x m] ambient per-vertex field, interleaved [x1;y1;z1;x2;...]. m = number of
%           channels (kernel) or time frames (full source).
%     Cov - the 'Covariant' operator node (tess_operators(SurfaceFile,'Covariant')).
%
% OUTPUTS:
%     Phi - [nVtot x m] scalar potential (sources/sinks)
%     Psi - [nVtot x m] stream function (vortices)
%
% SEE ALSO: process_helmholtz, bst_divergence, bst_curl, tess_cholesky, process_source_flow

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

    s = +1;
    nVtot = max(cellfun(@(c) max(double(c(:))), Cov.GlobalVertices));
    m = size(X, 2);
    Phi = zeros(nVtot, m);  Psi = zeros(nVtot, m);
    for hh = 1:numel(Cov.Covariant)
        C = Cov.Covariant{hh};  vH = double(Cov.GlobalVertices{hh}(:));
        nFh = size(C.Faces, 1);  nVh = numel(vH);
        Gx = C.ScalarGrad(1:nFh,:);  Gy = C.ScalarGrad(nFh+1:2*nFh,:);  Gz = C.ScalarGrad(2*nFh+1:3*nFh,:);
        Nf = C.FaceNormal;  Af = C.FaceArea;
        W   = spdiags(Af, 0, nFh, nFh);
        Fvf = sparse([(1:nFh)';(1:nFh)';(1:nFh)'], [C.Faces(:,1);C.Faces(:,2);C.Faces(:,3)], 1/3, nFh, nVh);
        nx = Nf(:,1); ny = Nf(:,2); nz = Nf(:,3);
        Sx = spdiags(ny,0,nFh,nFh)*Gz - spdiags(nz,0,nFh,nFh)*Gy;
        Sy = spdiags(nz,0,nFh,nFh)*Gx - spdiags(nx,0,nFh,nFh)*Gz;
        Sz = spdiags(nx,0,nFh,nFh)*Gy - spdiags(ny,0,nFh,nFh)*Gx;
        dF = tess_cholesky(Cov, hh, 1);
        Jx = X(3*(vH-1)+1, :);  Jy = X(3*(vH-1)+2, :);  Jz = X(3*(vH-1)+3, :);   % [nVh x m]
        Jfx = Fvf*Jx;  Jfy = Fvf*Jy;  Jfz = Fvf*Jz;                              % [nFh x m]
        divw  = s * (Gx'*(W*Jfx) + Gy'*(W*Jfy) + Gz'*(W*Jfz));                    % [nVh x m]
        vortw = s * (Sx'*(W*Jfx) + Sy'*(W*Jfy) + Sz'*(W*Jfz));
        phi = tess_cholesky('solve', dF, divw);   phi = phi - mean(phi, 1);      % [nVh x m]
        psi = tess_cholesky('solve', dF, vortw);  psi = psi - mean(psi, 1);
        Phi(vH,:) = phi;  Psi(vH,:) = psi;
    end
end
