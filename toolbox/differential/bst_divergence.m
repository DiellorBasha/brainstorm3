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
%            source vector). Wraps the Dirac Hodge engine (bst_helmholtz) and adds the
%            mean-curvature coupling term: div_Sigma(J) = div_Sigma(J_tan) - 2*H*(J.N),
%            where H is the scalar mean curvature and N is the outward vertex normal.
%            Call: bst_divergence(V, ManifoldMat, 'Ambient', Surf, Dir, LBO).
% I/O-free: the caller resolves and passes the loaded manifold node (and operator nodes for
% the ambient branch).
%
% USAGE:
%   divField = bst_divergence(V, ManifoldMat)
%   divField = bst_divergence(V, ManifoldMat, 'Ambient', Surf, Dir, LBO)
%
% INPUTS:
%   - V           : [TANGENT] per-FACE ambient [3nF x nT]; [AMBIENT] per-VERTEX [3nV x nT]
%   - ManifoldMat : a manifold_ node (tess_manifold(Surf)) whose 1x2 DEC group carries the
%                   per-hemisphere flat/d0/h0/h1 + GlobalVertices/GlobalFaces.
%   - 'Ambient'   : (optional) flag to dispatch to the ambient branch
%   - Surf        : (ambient only) loaded tessellation struct (in_tess_bst output)
%   - Dir         : (ambient only) Covariant operator node (tess_operators(Surf,'Covariant'))
%   - LBO         : (ambient only) LBO operator node (tess_operators(Surf,'Laplace-Beltrami'))
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

    % ----- ambient (3nV) branch: Hodge divergence + mean-curvature coupling -----
    if ~isempty(varargin) && strcmpi(varargin{1}, 'Ambient')
        Surf = varargin{2};  Dir = varargin{3};  LBO = varargin{4};
        divField = i_ambient_divergence(V, ManifoldMat, Surf, Dir, LBO);
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

%% ===== ambient divergence: the full flat-covariant surface divergence =====
% The 'Covariant' Hodge engine returns H.Div = the FULL ambient surface divergence of the
% 3-D current (it already includes the mean-curvature coupling -2H(J.N): a constant ambient
% field gives Div=0 even on folds). So no separate curvature term is added here -- doing so
% would double-count. (Op{1} is the 'Covariant' node; arg name kept for signature stability.)
function divField = i_ambient_divergence(J, ManifoldMat, Surf, Cov, LBO)
    H = bst_helmholtz('Decompose', {Cov, LBO}, ManifoldMat, Surf, J);
    divField = H.Div;
end
