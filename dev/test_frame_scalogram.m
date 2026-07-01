function tests = test_frame_scalogram
tests = functiontests(localfunctions);
end

function ax = i_toy_axes()
% One-hemisphere toy: identity Phi over K modes, unit mass, lambda grid, gv=1..K.
K = 40;  I = eye(K);
ax = struct();
ax.Phi = {I};  ax.Mass = {speye(K)};  ax.Lambda = {linspace(0.1,10,K)'};  ax.GlobalVertices = {(1:K)'};
end

function test_tight_frame_residual_zero(tc)
ax = i_toy_axes();  K = numel(ax.Lambda{1});  nT = 5;
% itersine tight frame of 6 members over [0,lmax]
lmax = max(ax.Lambda{1});  Nf = 6;  gCell = cell(1,Nf);
for ii=1:Nf, gCell{ii} = bst_eigfilter_design_itersine(struct('member',ii,'Nf',Nf,'lmax',lmax)); end
C = {randn(K,nT)};
scal = bst_eigenwavelet('Scalogram', ax, gCell, C);
verifyLessThan(tc, scal.resScalar, 1e-6);                     % tight frame reconstructs exactly
verifyEqual(tc, size(scal.energy), [3 nT Nf]);
verifyEqual(tc, squeeze(scal.energy(1,:,:)), squeeze(scal.energy(2,:,:)) + squeeze(scal.energy(3,:,:)), 'AbsTol',1e-9); % global = LH+RH (RH=0 here)
end

function test_loose_frame_residual_positive(tc)
ax = i_toy_axes();  K = numel(ax.Lambda{1});  nT = 3;
gCell = { @(l) exp(-0.5*double(l(:))) };                      % single low-pass -> not tight
C = {randn(K,nT)};
scal = bst_eigenwavelet('Scalogram', ax, gCell, C);
verifyGreaterThan(tc, scal.resScalar, 0.1);                  % under-covers the spectrum
end
