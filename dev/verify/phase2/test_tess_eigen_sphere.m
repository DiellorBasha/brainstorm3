% TEST_TESS_EIGEN_SPHERE: analytic sphere spectrum + storage round-trip.
scratch = '/private/tmp/claude-501/-Users-diellorbasha-workspace-research-code-brainstorm3/8da2056d-09c9-4f92-8237-a34a82e5f5e8/scratchpad';
TestFile = fullfile(scratch, 'tess_twosphere_unit.mat');
% --- build: two unit spheres, offset, labeled lh/rh via Structures atlas ---
[Vs, Fs] = tess_sphere(2562);
n1 = size(Vs,1);
TessMat = db_template('surfacemat');
TessMat.Comment  = 'twosphere_unit';
TessMat.Vertices = [Vs; Vs + repmat([3 0 0], n1, 1)];
TessMat.Faces    = [Fs; Fs + n1];
TessMat.Atlas(2).Name = 'Structures';
TessMat.Atlas(2).Scouts(1) = db_template('scout');
TessMat.Atlas(2).Scouts(1).Label = 'lh';  TessMat.Atlas(2).Scouts(1).Region = 'LU';
TessMat.Atlas(2).Scouts(1).Vertices = 1:n1;  TessMat.Atlas(2).Scouts(1).Seed = 1;
TessMat.Atlas(2).Scouts(2) = db_template('scout');
TessMat.Atlas(2).Scouts(2).Label = 'rh';  TessMat.Atlas(2).Scouts(2).Region = 'RU';
TessMat.Atlas(2).Scouts(2).Vertices = n1+1:2*n1;  TessMat.Atlas(2).Scouts(2).Seed = n1+1;
bst_save(TestFile, TessMat, 'v7');
% --- solve K=50 ---
E = tess_eigen(TestFile, 'Laplace-Beltrami', 'nModes', 50);
for hh = 1:2
    Lam = E.Lambda{hh};
    assert(abs(Lam(1)) < 1e-8, 'lambda_1 must be ~0 (constant mode)');
    assert(issorted(Lam), 'eigenvalues must ascend');
    % analytic: lambda = l(l+1) on the unit sphere, multiplicity 2l+1
    lExp = []; for l = 0:7, lExp = [lExp, repmat(l*(l+1), 1, 2*l+1)]; end %#ok<AGROW>
    lExp = lExp(1:50)';
    relErr = abs(Lam - lExp) ./ max(lExp, 1);
    assert(max(relErr) < 0.05, 'sphere spectrum deviates >5%% (max rel err %g)', max(relErr));
end
% --- storage round-trip: field embedded, reusable, truncatable ---
S2 = load(TestFile);
assert(isfield(S2, 'Eigen') && isfield(S2.Eigen, 'LaplaceBeltrami'), 'Eigen field not embedded');
assert(S2.Eigen.LaplaceBeltrami.nModes == 50);
assert(~isempty(S2.History) && any(strcmpi(S2.History(:,2), 'eigen')), 'History row missing');
E2 = tess_eigen(TestFile, 'Laplace-Beltrami', 'nModes', 30);   % reuse + truncate
assert(size(E2.Phi{1}, 2) == 30 && isequal(E2.Phi{1}, E.Phi{1}(:,1:30)), 'truncated reuse failed');
S3 = load(TestFile);
assert(S3.Eigen.LaplaceBeltrami.nModes == 50, 'reuse must not shrink the stored basis');
E3 = tess_eigen(TestFile, 'Laplace-Beltrami', 'nModes', 50, 'ForceRecompute', 1);  % replace slot
S4 = load(TestFile);
assert(numel(fieldnames(S4.Eigen)) == 1, 'one slot per variant');
% --- residual + orthonormality on the recomputed basis ---
[Op, Ms] = tess_operators(TestFile, 'Laplace-Beltrami');
for hh = 1:2
    R = Op{hh}*E3.Phi{hh} - Ms{hh}*E3.Phi{hh}*diag(E3.Lambda{hh});
    N = (Op{hh} + Ms{hh})*E3.Phi{hh};
    assert(max(sqrt(sum(R.^2,1))./sqrt(sum(N.^2,1))) < 1e-10, 'pencil residual too large');
    assert(norm(E3.Phi{hh}'*Ms{hh}*E3.Phi{hh} - eye(50), 'fro') < 1e-10, 'not B-orthonormal');
end
delete(TestFile);
disp('test_tess_eigen_sphere PASSED');
