function test_dirac_eigenmode_leadfield_inverse
% The composed Dirac head model must be consumed unchanged by the existing
% mode-space inverse: bst_inverse_eigenmodes('SolvePure', ...) -> Kernel [2K x nCh].
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

nCh = 8; vH1 = (1:5)'; vH2 = (6:9)'; nVert = 9; K = 3;
GainU = zeros(nCh, 3*nVert);
for i = 1:nCh
    for j = 1:3*nVert
        GainU(i,j) = cos(i*j*pi/29);
    end
end
hemis = {vH1, vH2}; tags = {'L','R'};
DE = struct('Vectors',{[]},'Values',{[]},'Mass',{[]},'nModes',{[]},'Order',{[]}, ...
            'Tau',{[]},'GlobalVertices',{[]},'Hemisphere',{[]},'Provenance',{[]});
for hh = 1:2
    nVh = numel(hemis{hh});
    raw = reshape(sin((1:4*nVh)' * (1:(K+2)) * pi/17 + hh*0.3), 4*nVh, K+2);  % deterministic, full rank
    Phi = orth(raw); Phi = Phi(:, 1:K);                                         % orthonormal columns (B=I)
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
HeadModel = struct('Gain',GainU,'GridLoc',zeros(nVert,3),'GridOrient',zeros(nVert,3), ...
    'GridAtlas',[],'HeadModelType','surface','SurfaceFile','tess_cortex_test.mat','Comment','test');

CompHM = bst_dirac_eigenmode_leadfield(HeadModel, DE, 'nModes', K);

% Feed the composed leadfield into the mode-space inverse (pure-math entry).
% SolvePure signature: (L_tilde, lambdas, iW, Proj, Method, Prior, Alpha, SNR, Unreg)
iW   = eye(nCh);                 % trivial whitener
Proj = eye(nCh);                 % no SSP
Kernel = bst_inverse_eigenmodes('SolvePure', CompHM.Gain, CompHM.Eigenvalues, ...
    iW, Proj, 'mne', 'log', 1, 3, false);

assert(isequal(size(Kernel), [CompHM.nModes, nCh]), 'Kernel must be [2K x nChannels].');
assert(all(isfinite(Kernel(:))), 'Kernel must be finite.');

disp('ALL TESTS PASSED');
end
