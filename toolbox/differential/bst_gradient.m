function [gradField, gradMag] = bst_gradient(F, ManifoldMat)
% BST_GRADIENT: Surface gradient of a per-vertex scalar field on the cortical manifold.
%
% Computes the surface gradient as the DEC composition grad f = (d0 f)# = sharp * d0 * f, using
% the cached DEC PRIMITIVES (exterior derivative d0 and musical isomorphism sharp) carried by the
% manifold_ node (tess_manifold). Returns the gradient as the NATIVE per-FACE ambient vector
% field, plus (optionally) the per-face magnitude. This is a pure DATA operation: it COMPOSES
% primitives, it does not build operators -- assembly lives in tess_manifold (the manifold's
% calculus structure). I/O-free: the caller resolves and passes the loaded manifold node.
%
% USAGE:
%   [gradField, gradMag] = bst_gradient(F, ManifoldMat)
%
% INPUTS:
%   - F           : per-vertex scalar field [nV x nT] (a source map; one or more time samples)
%   - ManifoldMat : a manifold_ node (tess_manifold(Surf)), whose 1x2 DEC group carries per-
%                   hemisphere d0/sharp + GlobalVertices/GlobalFaces.
%
% OUTPUTS:
%   - gradField : per-FACE gradient vector field [3nF x nT], interleaved (x,y,z) per face,
%                 i.e. face k occupies rows 3k-2:3k. Faces outside both hemispheres stay 0.
%   - gradMag   : (optional) per-face gradient magnitude [nF x nT]
%
% SEE ALSO: bst_operators, tess_manifold (DEC group), bst_helmholtz

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

    if ~isfield(ManifoldMat, 'DEC') || isempty(ManifoldMat.DEC) || ~isfield(ManifoldMat.DEC, 'sharp')
        error('bst_gradient:noDEC', 'The manifold node has no DEC operators (recompute tess_manifold).');
    end
    DEC   = ManifoldMat.DEC;                                  % 1x2 struct array (per hemisphere)
    nH    = numel(DEC);
    allGF = {DEC.GlobalFaces};  allGF = allGF(~cellfun(@isempty, allGF));
    nFtot = max(cellfun(@(c) max(double(c(:))), allGF));
    nT    = size(F, 2);
    gradField = zeros(3*nFtot, nT);
    needMag = (nargout > 1);
    if needMag, gradMag = zeros(nFtot, nT); end

    for hh = 1:nH
        if isempty(DEC(hh).sharp), continue; end
        vH = double(DEC(hh).GlobalVertices(:));
        fH = double(DEC(hh).GlobalFaces(:));
        gf = DEC(hh).sharp * (DEC(hh).d0 * F(vH, :));         % grad f = (d0 f)# ; [3Fh x nT] interleaved
        gr = reshape([3*fH-2, 3*fH-1, 3*fH]', [], 1);         % local-face rows -> global-face rows
        gradField(gr, :) = gf;
        if needMag
            gradMag(fH, :) = sqrt(gf(1:3:end,:).^2 + gf(2:3:end,:).^2 + gf(3:3:end,:).^2);
        end
    end
end
