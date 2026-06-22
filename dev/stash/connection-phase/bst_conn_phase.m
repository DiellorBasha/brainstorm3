function R = bst_conn_phase(ConnEig, vFrame, varargin)
% BST_CONN_PHASE: Decode a connection-Laplacian eigenmode into a 3D field + phase.
%
% USAGE:  R = bst_conn_phase(ConnEig, vFrame)
%         R = bst_conn_phase(ConnEig, vFrame, 'Rank', 1, 'FsFrame', FsFrame, 'nSing', 2)
%
% DESCRIPTION:
%     Given a surface's connection eigenmodes (in_tess_conn_eigenmodes / M2) and
%     the nxr per-vertex tangent frame (nxr.manifold.measure.vertexFrame), decodes
%     the selected eigenmode's complex coordinates z into a gauge-independent 3D
%     tangent field  w = Re(z)*e1 + Im(z)*e2,  and derives its magnitude,
%     singularities, and the between-subject (FreeSurfer-gauge) phase. The eigen-
%     modes are block-structured per connected component (hemisphere); for each
%     component the column of within-component rank 'Rank' is used (Rank 1 = the
%     Fiedler / smoothest field). The intrinsic within-subject phase is NOT a
%     scalar here: it is rendered reference-free in the viewer (stripes / iso-phase
%     contours) from R.Field.
%
% INPUTS:
%     ConnEig : struct from in_tess_conn_eigenmodes (fields Vectors [nV x nModes]
%               complex, Component, CompRank).
%     vFrame  : struct with e1, e2, normals ([nV x 3]) from nxr vertexFrame.
%               NOTE: nxr normals point INWARD on FreeSurfer meshes (CW face
%               winding).  The frame is self-consistent for phase decoding (e1,
%               e2 are used, not normals).  For J·n̂ use TessMat.VertNormals.
%
% OPTIONS:
%     'Rank'    : within-component mode rank to use (default 1 = Fiedler).
%     'FsFrame' : struct with e1, e2 ([nV x 3]) per-vertex FreeSurfer frame; if
%                 given, R.Phase is the field angle in that frame. Default [] (NaN).
%     'nSing'   : singularities to report per component (default 2; Poincare-Hopf
%                 index sum on a hemisphere-sphere is 2).
%
% OUTPUTS:
%     R.Field        : [nV x 3] 3D tangent field (gauge-independent); 0 off-support.
%     R.Magnitude    : [nV x 1] |z|; 0 off-support.
%     R.Phase        : [nV x 1] FS-gauge phase in [-pi, pi]; NaN off-support / no FsFrame.
%     R.Singularities: [k x 1] vertex indices (the nSing smallest-|z| per component).
%                      MARKER HEURISTIC ONLY: these may be spatially adjacent (a
%                      single magnitude dip spread over neighbouring vertices) and
%                      are not guaranteed to be the distinct topological
%                      singularities; robust index-based detection is future work.
%     R.Rank         : the rank used.
%
% SEE ALSO: in_tess_conn_eigenmodes, bst_tangent_face2vertex, nxr.manifold.measure.vertexFrame

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

%% ===== PARSE INPUTS =====
Rank    = 1;
FsFrame = [];
nSing   = 2;
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'rank',    Rank = varargin{i+1};
        case 'fsframe', FsFrame = varargin{i+1};
        case 'nsing',   nSing = varargin{i+1};
    end
end

nV = size(vFrame.e1, 1);
if size(ConnEig.Vectors, 1) ~= nV
    error('bst_conn_phase:sizeMismatch', ...
        'ConnEig.Vectors (%d rows) must match vFrame (%d vertices).', size(ConnEig.Vectors,1), nV);
end

Field     = zeros(nV, 3);
Magnitude = zeros(nV, 1);
Phase     = nan(nV, 1);
Sing      = zeros(0, 1);

comps = unique(ConnEig.Component(:))';
for c = comps
    col = find(ConnEig.Component(:) == c & ConnEig.CompRank(:) == Rank, 1);
    if isempty(col)
        continue;
    end
    z   = ConnEig.Vectors(:, col);      % complex [nV x 1], nonzero only on component c
    idx = find(z ~= 0);                 % component support (block structure is exact-zero off it)
    zc  = z(idx);

    e1 = vFrame.e1(idx, :);
    e2 = vFrame.e2(idx, :);
    w  = real(zc) .* e1 + imag(zc) .* e2;   % 3D tangent field
    Field(idx, :)  = w;
    Magnitude(idx) = abs(zc);

    % Singularities: the nSing smallest-|z| vertices of this component.
    [~, ord] = sort(abs(zc), 'ascend');
    Sing = [Sing; idx(ord(1:min(nSing, numel(ord))))]; %#ok<AGROW>

    % Between-subject phase: the field angle in the per-vertex FreeSurfer frame.
    if ~isempty(FsFrame)
        U  = FsFrame.e1(idx, :);
        Vw = FsFrame.e2(idx, :);
        Phase(idx) = atan2(sum(w .* Vw, 2), sum(w .* U, 2));
    end
end

R = struct('Field', Field, 'Magnitude', Magnitude, 'Phase', Phase, ...
           'Singularities', Sing, 'Rank', Rank);
end
