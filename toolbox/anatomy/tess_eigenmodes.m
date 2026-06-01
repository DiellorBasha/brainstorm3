function [Eigenmodes, L, M, Vertices, Faces] = tess_eigenmodes(Vertices, Faces, varargin)
% TESS_EIGENMODES: Compute Laplace-Beltrami eigenmodes of a triangle mesh.
%
% USAGE:  [Eigenmodes, L, M, V, F] = tess_eigenmodes(Vertices, Faces)
%         [Eigenmodes, L, M, V, F] = tess_eigenmodes(Vertices, Faces, 'nModes', 300)
%
% DESCRIPTION:
%     Computes the eigenmodes (eigenvectors and eigenvalues) of the
%     Laplace-Beltrami operator on a triangulated surface. These are the
%     solutions to the generalized eigenvalue problem:
%
%         L * phi_k = lambda_k * M * phi_k
%
%     where L is the cotangent Laplacian and M is the mass matrix. The
%     eigenmodes form an M-orthonormal basis for scalar functions on the
%     surface, analogous to Fourier modes on flat domains.
%
%     The pipeline:
%       1. Validate the mesh is manifold (repair only if Repair=true).
%       2. Assemble whole-mesh L and M via tess_laplacian (returned to callers).
%       3. Split the mesh into connected components.
%       4. Solve the generalized eigenvalue problem independently per component.
%       5. Tag each mode with its component id and within-component rank.
%       6. Package into the Eigenmodes structure.
%
% INPUTS:
%     Vertices : [nVertices x 3] vertex positions
%     Faces    : [nFaces x 3] triangle vertex indices (1-based)
%
% OPTIONS:
%     nModes     : Number of eigenmodes to compute (default: 300)
%     MassType   : Mass matrix type: 'barycentric', 'voronoi', 'galerkin'
%                  (default: 'barycentric')
%     RemoveDC   : Remove the DC (constant) mode(s) (default: true)
%     Repair     : Attempt to repair non-manifold defects before assembly
%                  (default: false). When false, a non-manifold surface raises
%                  an error. Repair changes vertex/edge counts (risky).
%     Sigma      : Shift for shift-invert eigs (default: -1e-8).
%                  A small negative shift improves convergence for the
%                  smallest eigenvalues near zero.
%     Tolerance  : Convergence tolerance for eigs (default: 1e-10)
%     MaxIter    : Maximum Arnoldi iterations for eigs (default: 1000)
%     Verbose    : Print progress (default: 1)
%
% OUTPUTS:
%     Eigenmodes : Structure with fields:
%       .Vectors    : [nVertices x nModes] M-orthonormal eigenvectors
%       .Values     : [nModes x 1] eigenvalues (ascending, >= 0)
%       .nModes     : Number of modes (after DC removal if applicable)
%       .Component  : [nModes x 1] connected-component id of each mode
%       .CompRank   : [nModes x 1] within-component rank of each mode
%       .nComponents: Number of connected components solved
%       .MassType   : Mass matrix type used
%       .Sigma      : Shift parameter used
%       .Tolerance  : Convergence tolerance used
%       .nRemoved   : Number of DC modes removed (0 if RemoveDC=false)
%       .ComputeTime: Total computation time in seconds
%
%     L : [nV x nV] sparse cotangent Laplacian (positive semidefinite)
%     M : [nV x nV] sparse mass matrix
%
% NOTES:
%     - For a closed manifold, one DC mode exists (constant function).
%     - Each connected component contributes one DC mode; when RemoveDC=true,
%       nRemoved equals the number of closed components in the surface.
%     - The eigenvalues are frequencies^2 of the surface harmonics.
%     - Spectral projection of a scalar field u onto eigenmodes:
%         c_k = phi_k' * M * u    (M-weighted inner product)
%     - Reconstruction from coefficients:
%         u_approx = Vectors * c   (sum of phi_k * c_k)
%
% SEE ALSO: tess_laplacian, tess_manifold, in_tess_eigenmodes, out_tess_eigenmodes, bst_eigenmodes_project, bst_eigenmodes_filter

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
nModes    = 300;
MassType  = 'barycentric';
RemoveDC  = true;
Repair    = false;
Sigma     = -1e-8;
Tolerance = 1e-10;
MaxIter   = 1000;
Verbose   = 1;

for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'nmodes',    nModes    = varargin{i+1};
        case 'masstype',  MassType  = varargin{i+1};
        case 'removedc',  RemoveDC  = varargin{i+1};
        case 'repair',    Repair    = varargin{i+1};
        case 'sigma',     Sigma     = varargin{i+1};
        case 'tolerance', Tolerance = varargin{i+1};
        case 'maxiter',   MaxIter   = varargin{i+1};
        case 'verbose',   Verbose   = varargin{i+1};
    end
end

tStart = tic;

%% ===== STEP 1: VALIDATE MESH (NEVER AUTO-REPAIR) =====
% Repair changes vertex/edge counts, which desyncs head models, lead fields,
% and atlases built on this surface. So validate only; repair is opt-in.
[~, ~, isManifold, report] = tess_manifold(Vertices, Faces, 'Repair', 0, 'Verbose', Verbose);
if ~isManifold
    if Repair
        if Verbose
            fprintf('BST> tess_eigenmodes: Surface non-manifold; attempting repair (vertex count may change)...\n');
        end
        [Vertices, Faces] = tess_manifold(Vertices, Faces, 'Repair', 1, 'Verbose', Verbose);
    else
        if isfield(report, 'summary') && ~isempty(report.summary)
            defects = strjoin(report.summary, '; ');
        else
            defects = 'non-manifold mesh';
        end
        error(['tess_eigenmodes: Surface is not a clean 2-manifold (%s). Re-mesh with ' ...
               'icosphere downsampling, or pass ''Repair'',true to attempt a risky repair ' ...
               'that changes the vertex count.'], defects);
    end
end
nV = size(Vertices, 1);

% Ensure nModes doesn't exceed mesh size
% (need at least nModes + a few extra for eigs convergence)
maxModes = nV - 2;  % leave room for constant modes
if nModes > maxModes
    if Verbose
        fprintf('BST> tess_eigenmodes: Requested %d modes but mesh has %d vertices. Reducing to %d.\n', ...
            nModes, nV, maxModes);
    end
    nModes = maxModes;
end

%% ===== STEP 2: WHOLE-MESH OPERATORS (returned to callers as L, M outputs) =====
if Verbose
    fprintf('BST> tess_eigenmodes: Assembling L and M (%s mass)...\n', MassType);
end
[L, M] = tess_laplacian(Vertices, Faces, 'MassType', MassType);

%% ===== STEP 3: SPLIT INTO CONNECTED COMPONENTS =====
% Each disconnected component (e.g. a cortical hemisphere) is a separate
% manifold; its Laplace-Beltrami operator must be solved independently. Solving
% the union as one block-diagonal operator interleaves components by eigenvalue
% and, for near-degenerate eigenvalues (near-symmetric hemispheres), lets eigs
% return meaningless cross-component mixtures.
compId = conncomp(graph(tess_vertconn(Vertices, Faces)));   % [1 x nV]
nComp  = max(compId);
if Verbose
    fprintf('BST> tess_eigenmodes: %d connected component(s).\n', nComp);
end

%% ===== STEP 4: SOLVE PER COMPONENT AND ASSEMBLE =====
VectorsAll = zeros(nV, 0);
ValuesAll  = zeros(0, 1);
Component  = zeros(0, 1);
CompRank   = zeros(0, 1);
nRemoved   = 0;
for c = 1:nComp
    vIdx = find(compId == c);
    nvc  = numel(vIdx);
    % Sub-mesh: faces fully inside this component, vertex indices remapped to 1..nvc
    faceMask = all(ismember(Faces, vIdx), 2);
    remap = zeros(nV, 1);
    remap(vIdx) = 1:nvc;
    Fsub = remap(Faces(faceMask, :));
    Vsub = Vertices(vIdx, :);
    % Cap to this component's size; -2 leaves room for the DC mode + an eigs margin
    nModesC = min(nModes, nvc - 2);
    [Uc, lambdasC, nRemC] = SolveComponent(Vsub, Fsub, nModesC, MassType, ...
        RemoveDC, Sigma, Tolerance, MaxIter, Verbose);
    % Zero-pad into full vertex space and append
    Ufull = zeros(nV, size(Uc, 2));
    Ufull(vIdx, :) = Uc;
    VectorsAll = [VectorsAll, Ufull];                              %#ok<AGROW>
    ValuesAll  = [ValuesAll;  lambdasC(:)];                        %#ok<AGROW>
    Component  = [Component;  c * ones(numel(lambdasC), 1)];       %#ok<AGROW>
    CompRank   = [CompRank;   (1:numel(lambdasC))'];               %#ok<AGROW>
    nRemoved   = nRemoved + nRemC;
end

%% ===== STEP 5: PACKAGE OUTPUT =====
totalTime = toc(tStart);
Eigenmodes = struct();
Eigenmodes.Vectors     = VectorsAll;
Eigenmodes.Values      = ValuesAll;
Eigenmodes.nModes      = numel(ValuesAll);
Eigenmodes.Component   = Component;
Eigenmodes.CompRank    = CompRank;
Eigenmodes.nComponents = nComp;
Eigenmodes.MassType    = MassType;
Eigenmodes.Sigma       = Sigma;
Eigenmodes.Tolerance   = Tolerance;
Eigenmodes.nRemoved    = nRemoved;
Eigenmodes.ComputeTime = totalTime;

if Verbose
    fprintf('BST> tess_eigenmodes: Done. %d modes across %d component(s), total %.1f s.\n', ...
        Eigenmodes.nModes, nComp, totalTime);
end

end


%% ===== HELPER: SOLVE ONE COMPONENT =====
function [U, lambdas, nRemoved] = SolveComponent(Vertices, Faces, nModes, MassType, RemoveDC, Sigma, Tolerance, MaxIter, Verbose)
    [L, M] = tess_laplacian(Vertices, Faces, 'MassType', MassType);
    nV = size(Vertices, 1);
    nExtra   = 5;
    nRequest = min(nModes + nExtra, nV - 1);
    opts = struct('tol', Tolerance, 'maxit', MaxIter, 'disp', 0);
    try
        [U, D] = eigs(L, M, nRequest, Sigma, opts);
    catch ME
        if Verbose
            fprintf('BST> tess_eigenmodes: Shift-invert failed (%s). Trying smallestabs...\n', ME.message);
        end
        [U, D] = eigs(L, M, nRequest, 'smallestabs', opts);
    end
    lambdas = real(diag(D));
    [lambdas, sortIdx] = sort(lambdas);
    U = real(U(:, sortIdx));
    if Verbose && (numel(lambdas) < nRequest)
        fprintf('BST> tess_eigenmodes: Only %d/%d modes converged on this component.\n', numel(lambdas), nRequest);
    end
    % Remove DC modes (one per closed component)
    nRemoved = 0;
    if RemoveDC
        dcThreshold = max(abs(Sigma) * 10, 1e-6);
        dcMask = abs(lambdas) < dcThreshold;
        nRemoved = sum(dcMask);
        lambdas = lambdas(~dcMask);
        U = U(:, ~dcMask);
    end
    % Trim to requested count
    if numel(lambdas) > nModes
        lambdas = lambdas(1:nModes);
        U = U(:, 1:nModes);
    end
    % M-orthonormalize and clamp tiny negatives
    U = mOrthonormalize(U, M);
    lambdas(lambdas < 0) = 0;
end


%% ===== HELPER: M-ORTHONORMALIZATION =====
function U = mOrthonormalize(U, M)
% Enforce U' * M * U = I via Cholesky or SVD whitening.
%
% Approach: G = U' * M * U. If G is well-conditioned, use Cholesky:
%   G = R' * R  =>  U_new = U * inv(R)  =>  U_new' * M * U_new = I
% Fallback to eigenvalue decomposition if Cholesky fails.

    G = U' * (M * U);
    G = (G + G') / 2;  % enforce symmetry

    % Try Cholesky first (fast, stable when G is well-conditioned)
    [R, flag] = chol(G);
    if flag == 0
        % Cholesky succeeded
        U = U / R;  % equivalent to U * inv(R)
    else
        % Fallback: eigenvalue decomposition
        [V, D] = eig(G);
        d = diag(D);
        % Clamp small/negative eigenvalues for stability
        d = max(d, max(d) * 1e-12);
        U = U * (V * diag(1 ./ sqrt(d)) * V');
    end
end
