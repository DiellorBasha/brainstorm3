function test_eigenmode_reconstruct_pure
% Verify bst_eigenmode_reconstruct: cortex kernel = Phi(:,1:K) * ModeKernel.
% Numeric-Phi mode keeps this DB-free.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

nVert = 20; K = 6; nCh = 9;
Phi = zeros(nVert, K+2);                 % surface carries more modes than used
for v = 1:nVert
    for k = 1:(K+2)
        Phi(v,k) = cos(v * k * pi / 17);
    end
end
ModeKernel = reshape(linspace(-1, 1, K*nCh), K, nCh);   % [K x nCh]

KVert = bst_eigenmode_reconstruct(Phi, ModeKernel);
assert(isequal(size(KVert), [nVert nCh]), 'Cortex kernel must be [nVert x nCh].');
assert(max(abs(KVert - Phi(:,1:K) * ModeKernel), [], 'all') < 1e-12, 'Must equal Phi(:,1:K)*ModeKernel.');

% Error when surface has fewer modes than the kernel
threw = false;
try
    bst_eigenmode_reconstruct(Phi(:,1:K-1), ModeKernel);
catch
    threw = true;
end
assert(threw, 'Must error when surface modes < kernel modes.');

disp('ALL TESTS PASSED');
end
