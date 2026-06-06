function [K, M, Info] = tess_connection_laplacian(Vertices, Faces, varargin)
% TESS_CONNECTION_LAPLACIAN: Complex-Hermitian vertex connection Laplacian (nxr).
%
% USAGE:  [K, M, Info] = tess_connection_laplacian(Vertices, Faces)
%         [K, M, Info] = tess_connection_laplacian(Vertices, Faces, 'nSym', 1)
%
% DESCRIPTION:
%     Assembles the connection Laplacian on the n-RoSy tangent bundle of a
%     triangle mesh, using nxr-compute (geometry-central). The connection is the
%     intrinsic discrete Levi-Civita connection: transport rotations along each
%     halfedge are raised to the nSym power. The operator is complex Hermitian and
%     positive (semi-)definite; its eigenvalues are intrinsic to the mesh and its
%     eigenvector phase is gauge-dependent (a globally consistent reference frame
%     is required to read out comparable phase — a later milestone).
%
%     This is a PURE OPERATOR: it depends on (Vertices, Faces) alone and does not
%     consume or store any tangent frame. nxr is REQUIRED (no MATLAB fallback);
%     the operator needs geometry-central's Levi-Civita transport vectors.
%
%     FRAME HANDEDNESS / WINDING CONVENTION
%     FreeSurfer cortex meshes are wound clockwise when viewed from outside the
%     brain.  geometry-central computes face normals via the CCW cross product,
%     so on a CW-wound mesh all normals point INWARD.  The vertex tangent basis
%     is then e2 = n_inward × e1, which encodes clockwise (rather than the usual
%     counter-clockwise) complex rotation.  This is consistent across all
%     FreeSurfer subjects, so inter-subject phase comparison is valid; the
%     absolute winding direction just reverses relative to an outward-normal
%     convention.  To get the outward normal use TessMat.VertNormals, not
%     nxr vertexFrames.normals.
%
% INPUTS:
%     Vertices : [nV x 3] vertex positions.
%     Faces    : [nF x 3] triangle vertex indices (1-based, Brainstorm convention).
%
% OPTIONS:
%     nSym           : Bundle symmetry (default 1). 1 = true vector field (carries
%                      phase); 2 = line field; 4 = cross field.
%     Domain         : 'vertex' (default) | 'face' | 'edge'. Vertex is the target;
%                      M is returned only for the vertex domain.
%     Regularization : Diagonal epsilon added by nxr for strict positive-
%                      definiteness (default 1e-8).
%     CheckManifold  : (logical) Pre-check the mesh with tess_manifold for a
%                      friendlier error (default false).
%
% OUTPUTS:
%     K    : [N x N] complex Hermitian sparse connection Laplacian (N = nV for the
%            vertex domain), assembled as K_real + 1i*K_imag.
%     M    : [N x N] real diagonal lumped vertex mass (area/3 per vertex) for the
%            vertex domain; [] for 'face'/'edge' (masses for those domains are a
%            later-milestone concern).
%     Info : struct with fields nSym, Domain, Regularization, baseDim, Format
%            ('complex'), Backend ('nxr').
%
% NOTE: nxr's eigensolver ('solve') is real-only and cannot consume this complex
%     K. Solve the generalized problem K*phi = lambda*M*phi with MATLAB eigs.
%
% SEE ALSO: tess_laplacian, tess_tangents, tess_manifold

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
nSym           = 1;
Domain         = 'vertex';
Regularization = 1e-8;
CheckManifold  = false;
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'nsym',           nSym = varargin{i+1};
        case 'domain',         Domain = lower(varargin{i+1});
        case 'regularization', Regularization = varargin{i+1};
        case 'checkmanifold',  CheckManifold = varargin{i+1};
    end
end

% Validate
if (size(Vertices, 2) ~= 3) || (size(Faces, 2) ~= 3)
    error('Vertices and Faces must have 3 columns.');
end
% nxr.manifold.context validates {'double'}; Brainstorm vertices are often single.
Vertices = double(Vertices);
Faces    = double(Faces);
if ~ismember(Domain, {'vertex', 'face', 'edge'})
    error('tess_connection_laplacian:badDomain', ...
        'Domain must be ''vertex'', ''face'', or ''edge'' (got ''%s'').', Domain);
end

%% ===== REQUIRE NXR (no MATLAB fallback) =====
if ~nxr_is_loaded()
    error('tess_connection_laplacian:nxrNotLoaded', ...
        ['The connection Laplacian requires the nxr-compute plugin (no MATLAB ' ...
         'fallback). Install/load it via bst_plugin(''Install'',''nxr-compute'').']);
end

%% ===== OPTIONAL MANIFOLD CHECK =====
if CheckManifold
    [~, ~, isManifold] = tess_manifold(Vertices, double(Faces), 'Repair', 0, 'Verbose', 0);
    if ~isManifold
        error('tess_connection_laplacian:NonManifold', ...
            ['Input mesh is not a clean 2-manifold; the connection Laplacian ' ...
             'requires one. Validate/repair with tess_manifold.']);
    end
end

%% ===== ASSEMBLE (complex Hermitian) =====
% Faces are passed 1-based; the nxr MEX marshalling subtracts 1 (proven to
% machine precision by the tess_laplacian parity test).
mctx = nxr.manifold.context(Vertices, Faces);
opts = struct('domain', Domain, 'nSym', nSym, ...
              'regularization', Regularization, 'format', 'complex');
CL = nxr.manifold.operator.connectionLaplacian(mctx, opts);

% nxr returns the complex Hermitian operator as parallel real sparse parts.
K = CL.K_real + 1i * CL.K_imag;
% Clear floating-point asymmetry so downstream eigs() detects a Hermitian
% operator and returns real eigenvalues.
K = (K + K') / 2;

%% ===== MASS (vertex domain only) =====
% Lumped vertex mass (area/3 per vertex) = ctx.M. For the complex bundle the
% inner product is <u,v> = sum_i area_i * conj(u_i) * v_i, so the mass is the
% ordinary real diagonal lumped mass (no block duplication for complex format).
% Face/edge-domain masses are a later-milestone concern.
if strcmp(Domain, 'vertex')
    M = mctx.M;
else
    M = [];
end

%% ===== INFO =====
Info = struct();
Info.nSym           = nSym;
Info.Domain         = Domain;
Info.Regularization = Regularization;
Info.baseDim        = CL.baseDim;
Info.Format         = 'complex';
Info.Backend        = 'nxr';

end


function tf = nxr_is_loaded()
% True only if the nxr-compute plugin is currently loaded (cheap, in-memory).
% Local copy of the same guard in tess_laplacian (intentionally not shared).
tf = false;
try
    PlugDesc = bst_plugin('GetInstalled', 'nxr-compute');
    tf = ~isempty(PlugDesc) && isfield(PlugDesc, 'isLoaded') && PlugDesc.isLoaded;
catch
    tf = false;   % bst_plugin unavailable / any error -> treat as not loaded
end
end
