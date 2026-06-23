function Fr = bst_operator_frame(OperatorMat, h)
% BST_OPERATOR_FRAME: Canonical per-vertex tangent frame for a vector operator hemisphere.
%
% USAGE:  Fr = bst_operator_frame(OperatorMat, h)
%
% Returns the struct Fr = struct('e1',[V x 3], 'e2',[V x 3], 'normal',[V x 3]) -- the
% canonical tangent frame in which the operator's complex eigenmodes decode:
%     field(v) = real(U(v,k)) .* Fr.e1(v,:) + imag(U(v,k)) .* Fr.e2(v,:)
% rows aligned with OperatorMat.GlobalVertices{h}.
%
% The frame is the deterministic geometry-central vertexTangentBasis (raw / levi-civita)
% that nxr uses to BUILD the operator, so this is provably the operator's own frame
% (see dev/tests/test_frame_sync.m). Prefers the stored OperatorMat.Frame{h}; if absent
% (operator files written before frame capture) it recomputes via nxr_compute('vertexFrames')
% on the same submesh -- bitwise identical by determinism.
%
% SEE ALSO: tess_operators, bst_eigenwavelet, panel_eigenwavelet_options

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

    % Preferred path: the frame stored with the operator (exact, no recompute)
    if isfield(OperatorMat, 'Frame') && iscell(OperatorMat.Frame) ...
            && (numel(OperatorMat.Frame) >= h) && ~isempty(OperatorMat.Frame{h})
        Fr = OperatorMat.Frame{h};
        return;
    end

    % Back-compat: recompute the canonical frame via nxr on the same submesh.
    if ~isfield(OperatorMat, 'GlobalVertices') || (numel(OperatorMat.GlobalVertices) < h) ...
            || isempty(OperatorMat.GlobalVertices{h})
        error('bst_operator_frame:noHemi', 'OperatorMat has no GlobalVertices for hemisphere %d.', h);
    end
    if ~isfield(OperatorMat, 'ParentSurface') || isempty(OperatorMat.ParentSurface)
        error('bst_operator_frame:noSurface', 'OperatorMat has no ParentSurface to recompute the frame.');
    end
    T  = in_tess_bst(OperatorMat.ParentSurface, 0);
    gv = OperatorMat.GlobalVertices{h}(:);
    nV = numel(gv);
    isV   = false(size(T.Vertices,1), 1);  isV(gv) = true;
    fMask = all(isV(double(T.Faces)), 2);
    mapV  = zeros(size(T.Vertices,1), 1);  mapV(gv) = 1:nV;
    Vloc  = T.Vertices(gv, :);
    Floc  = mapV(double(T.Faces(fMask, :)));
    hnx = nxr_compute('create', Vloc, Floc);
    VF  = nxr_compute('vertexFrames', hnx);
    nxr_compute('destroy', hnx);
    Fr = struct('e1', VF.e1, 'e2', VF.e2, 'normal', VF.normals);
end
