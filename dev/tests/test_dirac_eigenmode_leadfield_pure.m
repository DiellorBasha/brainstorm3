function test_dirac_eigenmode_leadfield_pure
% Verify the Dirac forward composer (synthetic, no eigs/cortex):
%   - composed Gain is [nCh x 2K] and equals the per-hemisphere quaternion projection
%   - embedding places the ambient x/y/z gain in the imaginary quaternion slots (w=0)
%   - Eigenvalues / ModeHemisphere / nModes carried through; basis stacked L then R
%   - non-surface and non-unconstrained head models are rejected
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

nCh = 6; vH1 = [1 2 3 4]; vH2 = [5 6 7]; nVert = 7; K = 3;

% deterministic unconstrained gain [nCh x 3*nVert]
GainU = zeros(nCh, 3*nVert);
for i = 1:nCh
    for j = 1:3*nVert
        GainU(i,j) = sin(i*j*pi/23);
    end
end

% synthetic per-hemisphere DiracEigen with Mass = I (so B = I; Vectors orthonormal)
hemis = {vH1(:), vH2(:)}; tags = {'L','R'};
DE = struct('Vectors',{[]},'Values',{[]},'Mass',{[]},'nModes',{[]},'Order',{[]}, ...
            'Tau',{[]},'GlobalVertices',{[]},'Hemisphere',{[]},'Provenance',{[]});
for hh = 1:2
    nVh = numel(hemis{hh});
    raw = reshape(sin((1:4*nVh)' * (1:(K+2)) * pi/17 + hh*0.3), 4*nVh, K+2);  % deterministic, full rank
    Phi = orth(raw); Phi = Phi(:, 1:K);                      % orthonormal columns (B=I)
    DE(hh).Vectors        = Phi;
    DE(hh).Values         = (1:K)' + 0.1*hh;
    DE(hh).Mass           = speye(nVh);
    DE(hh).nModes         = K;
    DE(hh).Order          = (1:K)';
    DE(hh).Tau            = 0.5;
    DE(hh).GlobalVertices = hemis{hh};
    DE(hh).Hemisphere     = tags{hh};
    DE(hh).Provenance     = struct('Backend','nxr');
end

HeadModel = struct('Gain', GainU, 'GridLoc', zeros(nVert,3), 'GridOrient', zeros(nVert,3), ...
    'GridAtlas', [], 'HeadModelType', 'surface', 'SurfaceFile', 'tess_cortex_test.mat', ...
    'Comment', 'OS-MEG test');

CompHM = bst_dirac(HeadModel, DE, 'nModes', K);

% --- shape + carried metadata ---
assert(isequal(size(CompHM.Gain), [nCh, 2*K]), 'Composed Gain must be [nCh x 2K].');
assert(CompHM.nModes == 2*K, 'nModes must be 2K.');
assert(isequal(CompHM.Eigenvalues(:), [DE(1).Values; DE(2).Values]), 'Eigenvalues stacked L then R.');
assert(isequal(CompHM.ModeHemisphere(:), [ones(K,1); 2*ones(K,1)]), 'ModeHemisphere L then R.');
assert(CompHM.isDiracEigenmode == 1, 'isDiracEigenmode flag set.');
assert(isempty(CompHM.GridLoc) && isempty(CompHM.GridOrient), 'GridLoc/Orient cleared.');
assert(strcmp(CompHM.SurfaceFile, 'tess_cortex_test.mat'), 'SurfaceFile carried through.');

% --- projection identity with an INDEPENDENT explicit embedding (oracle) ---
for hh = 1:2
    vH = DE(hh).GlobalVertices; nVh = numel(vH);
    Psi = zeros(4*nVh, nCh);
    for vloc = 1:nVh
        s = vH(vloc);
        g = GainU(:, 3*(s-1)+(1:3));     % [nCh x 3] ambient
        Psi(4*(vloc-1)+1, :) = 0;        % w
        Psi(4*(vloc-1)+2, :) = g(:,1)';  % x
        Psi(4*(vloc-1)+3, :) = g(:,2)';  % y
        Psi(4*(vloc-1)+4, :) = g(:,3)';  % z
    end
    B = kron(DE(hh).Mass, speye(4));
    Lref = Psi' * (B * DE(hh).Vectors);  % [nCh x K]
    block = CompHM.Gain(:, (hh-1)*K + (1:K));
    assert(max(abs(block(:) - Lref(:))) < 1e-9, sprintf('Hemisphere %d projection mismatch.', hh));
end

% --- errors ---
ok = false;
HMvol = HeadModel; HMvol.HeadModelType = 'volume';
try, bst_dirac(HMvol, DE); catch ME, ok = strcmp(ME.identifier,'bst_dirac:NotSurface'); end
assert(ok, 'Volume head model must raise NotSurface.');

ok = false;
HMbad = HeadModel; HMbad.Gain = GainU(:, 1:end-1);   % cols not a multiple of 3
try, bst_dirac(HMbad, DE); catch ME, ok = strcmp(ME.identifier,'bst_dirac:NotUnconstrained'); end
assert(ok, 'Non-unconstrained gain must raise NotUnconstrained.');

disp('ALL TESTS PASSED');
end
