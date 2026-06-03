function test_eigenmode_hemisphere_pure
% Regression for the one-hemisphere bug: the eigenmode leadfield + reconstruction
% must select modes by GLOBAL eigenvalue order across all connected components
% (hemispheres), not the first-K columns. Eigenmodes are stored grouped by
% component, so a naive Phi(:,1:K) keeps only one hemisphere.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

% Two components: vertices 1:4 = comp1 (LH), 5:8 = comp2 (RH). 4 modes each.
nV1 = 4; nV2 = 4; nV = nV1 + nV2; m = 4; nCh = 6;
Vectors = zeros(nV, 2*m);
for k = 1:m
    for n = 1:nV1; Vectors(n, k)         = cos(pi/nV1 * (n-0.5) * (k-1)); end  % comp1 cols 1:m
    for n = 1:nV2; Vectors(nV1+n, m+k)   = cos(pi/nV2 * (n-0.5) * (k-1)); end  % comp2 cols m+1:2m
end
% Eigenvalues stored grouped by component (comp1 then comp2), each ascending.
% Interleaved so the K=4 globally-lowest span BOTH components.
Values = [10;20;30;40;  11;21;31;41];
Eig = struct('Vectors', Vectors, 'Values', Values, 'nModes', 2*m, ...
    'Component', [ones(m,1); 2*ones(m,1)], 'CompRank', [(1:m)'; (1:m)'], 'nComponents', 2);

% Deterministic base head model (unconstrained gain + unit normals)
GainU = zeros(nCh, 3*nV);
for i = 1:nCh
    for j = 1:3*nV; GainU(i,j) = sin(i * j * pi / 17); end
end
GridOrient = zeros(nV,3);
for v = 1:nV; n = [cos(v), sin(v), 0.5]; GridOrient(v,:) = n / norm(n); end
HeadModel = struct('Gain', GainU, 'GridOrient', GridOrient, 'GridLoc', zeros(nV,3), ...
    'GridAtlas', [], 'HeadModelType', 'surface', 'SurfaceFile', 'x.mat', 'Comment', 't');

K = 4;   % < per-component count (4): the OLD first-K slice would be comp1 only
CompHM = bst_eigenmode_leadfield(HeadModel, Eig, 'nModes', K);

% --- leadfield must record which mode columns it used ---
assert(isfield(CompHM, 'ModeIndices'), 'leadfield must record selected mode indices (ModeIndices).');
sel = CompHM.ModeIndices(:);
assert(numel(sel) == K, 'ModeIndices length must equal K.');
% --- selection = K globally-lowest eigenvalues, spanning BOTH components ---
[~, ord] = sort(Values, 'ascend');
assert(isequal(sort(sel), sort(ord(1:K))), 'must select the K lowest-eigenvalue modes globally.');
comp = Eig.Component(sel);
assert(any(comp==1) && any(comp==2), 'selected modes must span BOTH hemispheres.');
assert(isequal(CompHM.Eigenvalues(:), Values(sel)), 'stored Eigenvalues must match the selection order.');

% --- reconstruction must use the same indices and cover BOTH hemispheres ---
ModeKernel = reshape(linspace(-1, 1, K*nCh), K, nCh);     % [K x nCh]
KVert = bst_eigenmode_reconstruct(Vectors, ModeKernel, sel);
assert(isequal(size(KVert), [nV nCh]), 'reconstruction shape must be [nVert x nCh].');
e1 = sum(sum(abs(KVert(1:nV1, :))));        % LH energy
e2 = sum(sum(abs(KVert(nV1+1:nV, :))));     % RH energy
assert(e1 > 0 && e2 > 0, 'reconstruction must have nonzero energy on BOTH hemispheres.');

% --- backward-compat: 2-arg reconstruct still works (old head models) ---
KVert2 = bst_eigenmode_reconstruct(Vectors, ModeKernel);
assert(isequal(KVert2, Vectors(:,1:K) * ModeKernel), '2-arg reconstruct must keep Phi(:,1:K) behavior.');

disp('ALL TESTS PASSED');
end
