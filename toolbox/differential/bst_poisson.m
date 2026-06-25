function phi = bst_poisson(OperatorNode, F)
% BST_POISSON: Stratified Poisson solve  L phi = f  on the cortical manifold.
%
% Scalar stratum (Laplace-Beltrami): solves  K phi = M f  per hemisphere, where K is the
% cotan stiffness and M the Galerkin mass. The constant nullspace is handled HERE, once:
% project f to the mean-zero subspace in the mass metric, pinned solve through the cached
% tess_cholesky factor. This is the single home of the nullspace handling that
% was duplicated in bst_operators (per-column re-factorization) and bst_helmholtz.
%
% USAGE:  phi = bst_poisson(OperatorNode, F)
%   OperatorNode : a 'Laplace-Beltrami' or 'Covariant' operatormat (Operator{hh}=K cotan stiffness, Mass{hh}=M, GlobalVertices{hh})
%   F            : per-vertex scalar source [nV x nT]
%   phi          : per-vertex potential [nV x nT] (pinned at vertex 1 per hemisphere; gauge-fixed)
%
% SEE ALSO: tess_cholesky, bst_operators, bst_helmholtz
%
% Authors: Diellor Basha, 2026

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

    if ~any(strcmpi(OperatorNode.Variant, {'Laplace-Beltrami', 'Covariant'}))
        error('bst_poisson:variant', ...
            'bst_poisson scalar route needs a Laplace-Beltrami or Covariant operator (got %s).', OperatorNode.Variant);
    end
    nVtot = max(cellfun(@(c) max(double(c(:))), OperatorNode.GlobalVertices));
    phi = zeros(nVtot, size(F,2));
    for hh = 1:numel(OperatorNode.Operator)
        if isempty(OperatorNode.Operator{hh}), continue; end
        vH = double(OperatorNode.GlobalVertices{hh}(:));
        M  = OperatorNode.Mass{hh};
        dF = tess_cholesky(OperatorNode, hh, 1);     % pin vertex 1
        fh = F(vH, :);
        totMass = sum(M(:));
        fh = fh - (sum(M*fh, 1) / totMass);          % project to mean-zero (mass metric)
        x  = tess_cholesky('solve', dF, M*fh);       % K x = M f  on the free block; pinned entry 0
        phi(vH, :) = x;                               % gauge-fixed: phi(pin)=0, unique solution
    end
end
