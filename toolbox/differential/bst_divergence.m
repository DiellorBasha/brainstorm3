function divField = bst_divergence(V, ManifoldMat, varargin)
% BST_DIVERGENCE: Divergence of a tangent or ambient vector field on the cortical manifold.
%
% On a 2-D surface the de Rham complex is L0 -d0-> L1 -d1-> L2 (scalar -> tangent vector ->
% scalar). The divergence is the NEGATIVE ADJOINT of the gradient (grad = #d0):
%   div V = -*0^-1 d0' #' W_F V ,   a per-VERTEX scalar,
% where W_F is the face-area mass (areas = 1/diag(*2)) and #=sharp. This composes the manifold's
% cached DEC primitives (sharp, d0, Hodge stars *0,*2) -- a pure DATA op that builds no operators
% -- and AVOIDS *1^-1 (the metric flat is unstable: *1=cotan has negative entries on obtuse
% triangles). By construction div(grad f) is a Galerkin Laplace-Beltrami operator (~ M^-1 K).
%
% Two field types are supported:
%   TANGENT: v1 per-FACE field [3nF x nT] (interleaved x,y,z per face; e.g. output of
%            bst_gradient). Uses stable DEC adjoint. Call: bst_divergence(V, ManifoldMat).
%   AMBIENT: 3D R^3 per-VERTEX field [3nV x nT] (interleaved x,y,z per vertex; e.g. a Dirac
%            source vector). Uses the flat-covariant strong divergence from the Covariant node
%            (Gx*Jx+Gy*Jy+Gz*Jz, area-weighted to vertices). Already includes the
%            mean-curvature coupling -2H(J.N): a constant ambient field gives Div=0 on folds.
%            Call: bst_divergence(V, ManifoldMat, 'Ambient', Surf, Cov).
% I/O-free: the caller resolves and passes the loaded manifold node (and operator nodes for
% the ambient branch).
%
% USAGE:
%   divField = bst_divergence(V, ManifoldMat)
%   divField = bst_divergence(V, ManifoldMat, 'Ambient', Surf, Cov)
%
% INPUTS:
%   - V           : [TANGENT] per-FACE ambient [3nF x nT]; [AMBIENT] per-VERTEX [3nV x nT]
%   - ManifoldMat : a manifold_ node (tess_manifold(Surf)) whose 1x2 DEC group carries the
%                   per-hemisphere flat/d0/h0/h1 + GlobalVertices/GlobalFaces.
%   - 'Ambient'   : (optional) flag to dispatch to the ambient branch
%   - Surf        : (ambient only) loaded tessellation struct (accepted for signature symmetry,
%                   unused: nVtot derived from Cov.GlobalVertices)
%   - Cov         : (ambient only) Covariant operator node (tess_operators(Surf,'Covariant'))
%
% OUTPUTS:
%   - divField    : per-VERTEX scalar field [nV x nT]
%
% SEE ALSO: bst_operators, bst_gradient, bst_curl, bst_helmholtz, tess_manifold (DEC group)
% (de Rham/Hodge route; the connection Laplacian differs by Gauss curvature K — see bst_operators.)

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

    % ----- ambient (3nV) branch: flat-covariant surface divergence (Covariant node) -----
    if ~isempty(varargin) && strcmpi(varargin{1}, 'Ambient')
        % USAGE: bst_divergence(J, ManifoldMat, 'Ambient', Surf, Cov)
        % ManifoldMat/Surf are accepted for signature symmetry with the tangent branch
        % but are unused here: the flat-covariant divergence needs only the Covariant node.
        Cov = varargin{end};
        divField = i_ambient_divergence(V, Cov);
        return;
    end
    % ----- tangent (3nF) branch: existing stable DEC adjoint (UNCHANGED) -----
    if ~isfield(ManifoldMat, 'DEC') || isempty(ManifoldMat.DEC) || ~isfield(ManifoldMat.DEC, 'sharp')
        error('bst_divergence:noDEC', 'The manifold node has no DEC operators (recompute tess_manifold).');
    end
    DEC   = ManifoldMat.DEC;
    nVtot = max(cellfun(@(c) max(double(c(:))), {DEC.GlobalVertices}));
    nT    = size(V, 2);
    divField = zeros(nVtot, nT);
    for hh = 1:numel(DEC)
        if isempty(DEC(hh).sharp), continue; end
        vH = double(DEC(hh).GlobalVertices(:));
        fH = double(DEC(hh).GlobalFaces(:));
        fr = reshape([3*fH-2, 3*fH-1, 3*fH]', [], 1);        % global-face rows -> local face vectors
        area = 1 ./ max(full(diag(DEC(hh).h2)), eps);        % face areas (*2 = diag(1/area))
        WFv  = repelem(area, 3) .* V(fr, :);                 % W_F V  (area-weight each face vector)
        h0inv = spdiags(1 ./ max(full(diag(DEC(hh).h0)), eps), 0, numel(vH), numel(vH));
        divField(vH, :) = -h0inv * (DEC(hh).d0' * (DEC(hh).sharp' * WFv));   % -*0^-1 d0' #' W_F V
    end
end

%% ===== ambient divergence: flat-covariant surface divergence (incl. -2H(J.N) coupling) =====
% Ported from bst_helmholtz i_prepare_vertex/i_frame_vertex. Per hemisphere, from the
% Covariant node: strong per-face divergence Gx*Jx+Gy*Jy+Gz*Jz, area-weighted to vertices by
% Wfv. s=+1 (calibrated: this IS the true surface divergence, already includes the
% mean-curvature coupling -2H(J.N); a constant ambient field gives 0 even on folds).
function divField = i_ambient_divergence(J, Cov)
    s = +1;
    nVtot = max(cellfun(@(c) max(double(c(:))), Cov.GlobalVertices));
    nT = size(J, 2);
    divField = zeros(nVtot, nT);
    for hh = 1:numel(Cov.Covariant)
        C = Cov.Covariant{hh};  vH = double(Cov.GlobalVertices{hh}(:));
        nFh = size(C.Faces, 1);  nVh = numel(vH);
        Gx = C.ScalarGrad(1:nFh,:);  Gy = C.ScalarGrad(nFh+1:2*nFh,:);  Gz = C.ScalarGrad(2*nFh+1:3*nFh,:);
        Wfv = bst_face2vertex(C.Faces, C.FaceArea);   % shared math helper
        Jx = J(3*(vH-1)+1, :);  Jy = J(3*(vH-1)+2, :);  Jz = J(3*(vH-1)+3, :);
        divF = Gx*Jx + Gy*Jy + Gz*Jz;                 % [nFh x nT] per-face surface divergence
        divField(vH, :) = s * (Wfv * divF);
    end
end
