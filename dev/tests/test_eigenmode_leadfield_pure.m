function test_eigenmode_leadfield_pure
% Verify the forward composer:
%   - composed Gain == constrained(L) * Phi(:,1:K), exactly
%   - constrained extraction matches bst_gain_orient on unconstrained gain
%   - K is clamped to available modes; eigenvalues carried through
%   - output is a valid composed head-model struct (isEigenmode, SurfaceFile, nModes)
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

% Synthetic: 6 channels, 5 vertices, unconstrained gain (3 cols/vertex), 4 modes.
nCh = 6; nV = 5; K = 4;
GainU = zeros(nCh, 3*nV);
for i = 1:nCh
    for j = 1:3*nV
        GainU(i,j) = sin(i * j * pi / 17);     % deterministic, full rank
    end
end
% Surface normals (unit) per vertex -> constrained leadfield via bst_gain_orient
GridOrient = zeros(nV,3);
for v = 1:nV
    n = [cos(v), sin(v), 0.5];
    GridOrient(v,:) = n / norm(n);
end
Lc = bst_gain_orient(GainU, GridOrient);        % [nCh x nV] constrained reference

% Eigenmodes struct (deterministic DCT-II basis, M-orthonormal not required for shape test)
Vectors = zeros(nV, K+1);                        % one extra mode to test clamping
for k = 1:(K+1)
    for n = 1:nV
        Vectors(n,k) = cos(pi/nV * (n-0.5) * (k-1));
    end
end
Values = (0:K)';                                 % ascending eigenvalues, includes DC=0
Eig = struct('Vectors', Vectors, 'Values', Values, 'nModes', K+1);

% Base head-model struct mimicking in_bst_headmodel output
HeadModel = struct('Gain', GainU, 'GridOrient', GridOrient, ...
    'GridLoc', zeros(nV,3), 'GridAtlas', [], ...
    'HeadModelType', 'surface', 'SurfaceFile', 'tess_cortex_test.mat', ...
    'Comment', 'OS-MEG test');

CompHM = bst_eigenmode_leadfield(HeadModel, Eig, 'nModes', K);

% --- assertions ---
assert(isequal(size(CompHM.Gain), [nCh, K]), 'Composed Gain must be [nCh x K].');
Lref = Lc * Vectors(:,1:K);
assert(max(abs(CompHM.Gain(:) - Lref(:))) < 1e-9, 'Composed Gain must equal constrained(L)*Phi.');
assert(CompHM.nModes == K, 'nModes must be clamped to K.');
assert(isequal(CompHM.Eigenvalues(:), Values(1:K)), 'Eigenvalues must be carried through (first K).');
assert(isfield(CompHM,'isEigenmode') && CompHM.isEigenmode == 1, 'isEigenmode flag must be set.');
assert(strcmp(CompHM.SurfaceFile, 'tess_cortex_test.mat'), 'SurfaceFile must be carried through.');
assert(isempty(CompHM.GridLoc) && isempty(CompHM.GridOrient), 'GridLoc/GridOrient must be empty on composed HM.');

% Clamp test: asking for more modes than available returns all available
CompHM2 = bst_eigenmode_leadfield(HeadModel, Eig, 'nModes', 999);
assert(CompHM2.nModes == (K+1), 'nModes must clamp to available count.');

disp('ALL TESTS PASSED');
end
