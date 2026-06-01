function test_eigenmode_lever_paired
% Paired-rank reindex: weights are over paired rank and expand to raw columns
% via CompRank; a low paired-band keeps BOTH components (hemisphere symmetry).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));
addpath(fullfile(repoRoot, 'toolbox', 'math'));
addpath(fullfile(repoRoot, 'toolbox', 'anatomy'));
global GlobalData; %#ok<GVMIS>

% Build a TWO-component eigenbasis by block-diagonalizing two single spheres.
[V1,F1] = tess_sphere(162);
[E1,~,M1] = tess_eigenmodes(V1, F1, 'nModes', 20, 'MassType','barycentric','RemoveDC',1,'Verbose',0);
[V2,F2] = tess_sphere(162);
[E2,~,M2] = tess_eigenmodes(V2, F2, 'nModes', 20, 'MassType','barycentric','RemoveDC',1,'Verbose',0);
K1 = E1.nModes; K2 = E2.nModes; Kp = min(K1,K2);
P1 = E1.Vectors(:,1:Kp); P2 = E2.Vectors(:,1:Kp);
nV1 = size(P1,1); nV2 = size(P2,1);
Vectors = [ P1, zeros(nV1,Kp); zeros(nV2,Kp), P2 ];
Values  = [ E1.Values(1:Kp); E2.Values(1:Kp) ];
Comp    = [ ones(Kp,1); 2*ones(Kp,1) ];
CompRank= [ (1:Kp)'; (1:Kp)' ];
Eig = struct('Vectors',Vectors,'Values',Values,'nModes',2*Kp, ...
             'Component',Comp,'CompRank',CompRank,'MassType','barycentric');
M = blkdiag(M1, M2);

SurfaceFile = '/synthetic/twocomp.mat';
panel_eigenmodes('ResetState', SurfaceFile, Kp);   % K_paired = Kp
panel_eigenmodes('SetCache', SurfaceFile, Eig, M);

rng(7); u = Vectors * randn(2*Kp, 1);
iC1 = 1:nV1; iC2 = nV1+(1:nV2);

panel_eigenmodes('SetActive', 1);
panel_eigenmodes('SetWindowShape', 'box');
panel_eigenmodes('SetBand', 1, 3);
uF = panel_eigenmodes('ApplyToColumn', SurfaceFile, u);

e1 = norm(uF(iC1)); e2 = norm(uF(iC2));
assert(e1 > 1e-6 && e2 > 1e-6, 'low paired-band must keep BOTH components (symmetry)');
assert(abs(e1 - e2)/max(e1,e2) < 0.6, 'kept energy should be comparable across components');

W = panel_eigenmodes('GetWeights');
assert(numel(W) == Kp, 'weights are paired-length');
wRaw = W(CompRank);
analytic = Vectors * (wRaw(:) .* (Vectors' * (M * u)));
assert(max(abs(uF - analytic)) < 1e-9, 'ApplyToColumn must expand via CompRank');

fprintf('ALL TESTS PASSED: test_eigenmode_lever_paired\n');
end
