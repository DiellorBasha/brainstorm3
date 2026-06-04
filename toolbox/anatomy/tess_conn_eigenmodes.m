function [ConnEig, K, M] = tess_conn_eigenmodes(Vertices, Faces, varargin)
% TESS_CONN_EIGENMODES: Eigenmodes of the connection Laplacian (vector-field basis).
%
% USAGE:  [ConnEig, K, M] = tess_conn_eigenmodes(Vertices, Faces)
%         [ConnEig, K, M] = tess_conn_eigenmodes(Vertices, Faces, 'nModes', 300)
%
% DESCRIPTION:
%     Computes the smallest eigenmodes of the complex-Hermitian connection
%     Laplacian (n-RoSy bundle, nSym=1) of a triangle mesh, solving the
%     generalized problem K*phi = lambda*M*phi independently per connected
%     component (the operator is block-diagonal across components, e.g. the two
%     cortical hemispheres). The result is a vector-field spectral basis, the
%     tangent-field sibling of the scalar Laplace-Beltrami eigenmodes
%     (tess_eigenmodes).
%
%     The operator and mass come from tess_connection_laplacian (nxr-compute,
%     no MATLAB fallback). nxr's own eigensolver is real-only, so the complex
%     Hermitian eigenproblem is solved with MATLAB eigs. There is NO DC/zero
%     mode (no globally consistent parallel vector field on a curved closed
%     surface), so none is removed. Eigenvectors are stored raw, in nxr's
%     intrinsic per-vertex frames; phase readout/gauge alignment is a later step.
%
% INPUTS:
%     Vertices : [nV x 3] vertex positions.
%     Faces    : [nF x 3] triangle vertex indices (1-based).
%
% OPTIONS:
%     nModes         : Modes per connected component (default 300).
%     nSym           : n-RoSy symmetry (default 1; true vector field).
%     Regularization : Operator diagonal epsilon (default 1e-8).
%     Tolerance      : eigs convergence tolerance (default 1e-10).
%     Verbose        : Print progress (default 1).
%
% OUTPUTS:
%     ConnEig : struct with fields Vectors [nV x nModes] complex (M-orthonormal,
%               block-structured per component), Values [nModes x 1] real, nModes,
%               Component, CompRank, Order, nComponents, MassMatrix (lumped, sparse),
%               ConnLaplacian (K, complex sparse), OperatorType, nSym,
%               Regularization, Sigma, Tolerance, nRemoved (=0), ComputeTime.
%     K       : [nV x nV] complex Hermitian connection Laplacian (whole mesh).
%     M       : [nV x nV] real diagonal lumped vertex mass (whole mesh).
%
% SEE ALSO: tess_connection_laplacian, tess_eigenmodes, out_tess_conn_eigenmodes

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
nModes         = 300;
nSym           = 1;
Regularization = 1e-8;
Tolerance      = 1e-10;
Verbose        = 1;
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'nmodes',         nModes = varargin{i+1};
        case 'nsym',           nSym = varargin{i+1};
        case 'regularization', Regularization = varargin{i+1};
        case 'tolerance',      Tolerance = varargin{i+1};
        case 'verbose',        Verbose = varargin{i+1};
    end
end
tStart = tic;

if (size(Vertices, 2) ~= 3) || (size(Faces, 2) ~= 3)
    error('Vertices and Faces must have 3 columns.');
end
Faces = double(Faces);
nV    = size(Vertices, 1);

%% ===== WHOLE-MESH OPERATOR + MASS (complex Hermitian) =====
[K, M] = tess_connection_laplacian(Vertices, Faces, 'nSym', nSym, 'Regularization', Regularization);

%% ===== CONNECTED COMPONENTS =====
% Solve each disconnected component (e.g. a hemisphere) independently. The
% connection Laplacian has no cross-component entries, so K(idx,idx) is the exact
% component operator; restricting beats rebuilding sub-meshes.
compId = conncomp(graph(tess_vertconn(Vertices, Faces)));   % [1 x nV]
nComp  = max(compId);
if Verbose
    fprintf('BST> tess_conn_eigenmodes: %d connected component(s).\n', nComp);
end

%% ===== SOLVE PER COMPONENT =====
VectorsAll = complex(zeros(nV, 0));
ValuesAll  = zeros(0, 1);
Component  = zeros(0, 1);
CompRank   = zeros(0, 1);
eigsOpts   = struct('tol', Tolerance);
for c = 1:nComp
    vIdx = find(compId == c);
    nvc  = numel(vIdx);
    kC   = min(nModes, nvc - 2);            % leave an eigs margin
    if kC < 1
        if Verbose
            fprintf('BST> tess_conn_eigenmodes: skipping component %d (%d vertices, too small).\n', c, nvc);
        end
        continue;
    end
    Kc   = K(vIdx, vIdx);
    Mc   = M(vIdx, vIdx);
    [Uc, Dc] = eigs(Kc, Mc, kC, 'smallestabs', eigsOpts);
    lam = real(diag(Dc));
    [lam, ord] = sort(lam, 'ascend');
    Uc  = Uc(:, ord);
    % eigs may return fewer than kC modes if not all converge; derive the count
    % from the actual result so the metadata vectors stay in sync.
    nGot = numel(lam);
    % Guarantee M-orthonormal magnitude (eigs returns B-normalized; this guards it).
    nrm = sqrt(real(sum(conj(Uc) .* (Mc * Uc), 1)));
    Uc  = Uc ./ nrm;
    Ufull = complex(zeros(nV, nGot));
    Ufull(vIdx, :) = Uc;
    VectorsAll = [VectorsAll, Ufull];                       %#ok<AGROW>
    ValuesAll  = [ValuesAll;  lam(:)];                      %#ok<AGROW>
    Component  = [Component;  c * ones(nGot, 1)];           %#ok<AGROW>
    CompRank   = [CompRank;   (1:nGot)'];                   %#ok<AGROW>
    if Verbose
        fprintf('BST> tess_conn_eigenmodes: component %d: %d modes, range [%.3g, %.3g].\n', ...
            c, nGot, lam(1), lam(end));
    end
end

%% ===== PACKAGE =====
ConnEig = struct();
ConnEig.Vectors       = VectorsAll;
ConnEig.Values        = ValuesAll;
ConnEig.nModes        = numel(ValuesAll);
ConnEig.Component     = Component;
ConnEig.CompRank      = CompRank;
% Canonical global order: ascending eigenvalue across all components, so
% Vectors(:,Order(1:K)) are the whole-brain lowest-frequency modes.
[~, ConnEig.Order]    = sort(ValuesAll, 'ascend');
ConnEig.nComponents   = nComp;
ConnEig.MassMatrix    = M;                       % lumped vertex mass (basis is M-orthonormal)
ConnEig.ConnLaplacian = K;                       % complex Hermitian operator (reuse downstream)
ConnEig.OperatorType  = 'Connection-LeviCivita';
ConnEig.nSym          = nSym;
ConnEig.Regularization = Regularization;
ConnEig.Sigma         = 'smallestabs';
ConnEig.Tolerance     = Tolerance;
ConnEig.nRemoved      = 0;                        % no DC mode in the connection bundle
ConnEig.ComputeTime   = toc(tStart);
end
