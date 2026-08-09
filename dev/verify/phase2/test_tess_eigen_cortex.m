% TEST_TESS_EIGEN_CORTEX: K=400/hemisphere on the real ico5 cortex (a COPY).
scratch = '/private/tmp/claude-501/-Users-diellorbasha-workspace-research-code-brainstorm3/8da2056d-09c9-4f92-8237-a34a82e5f5e8/scratchpad';
Orig = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/verify/phase0/bst_userdir_clean/.brainstorm/local_db/omega-tutorial-cortical-flow/anat/sub-0002/tess_cortex_pial_low.mat';
Work = fullfile(scratch, 'tess_cortex_eigen_work.mat');
copyfile(Orig, Work);  fileattrib(Work, '+w');
info0 = dir(Work);
t0 = tic;
E = tess_eigen(Work, 'Laplace-Beltrami', 'nModes', 400);
tSolve = toc(t0);
fprintf('solve time: %.1f s\n', tSolve);
% --- residual + orthonormality vs the INDEPENDENT nxr oracle pencil ---
S = load('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/verify/phase1/oracle_lbo_sub0002.mat');
for hh = 1:2
    A = S.A{hh}; B = S.B{hh};
    R = A*E.Phi{hh} - B*E.Phi{hh}*diag(E.Lambda{hh});
    N = (A + B)*E.Phi{hh};
    res = max(sqrt(sum(R.^2,1)) ./ sqrt(sum(N.^2,1)));
    orth = norm(E.Phi{hh}'*B*E.Phi{hh} - eye(400), 'fro');
    fprintf('hemi %d: residual(vs nxr pencil)=%.3g, orth=%.3g, lambda1=%.3g, lambda400=%.6g\n', ...
        hh, res, orth, E.Lambda{hh}(1), E.Lambda{hh}(400));
    assert(res < 1e-8, 'residual vs independent nxr pencil too large');   % cross-implementation
    assert(orth < 1e-10, 'not B-orthonormal');
    assert(abs(E.Lambda{hh}(1)) < 1e-6 * E.Lambda{hh}(400), 'zero mode not recovered');
end
% --- storage: file growth, reuse speed, replace semantics ---
info1 = dir(Work);
fprintf('file: %.1f MB -> %.1f MB\n', info0.bytes/1e6, info1.bytes/1e6);
% 2x10242x400 doubles = 65.5 MB raw; -v7 gzip on smooth eigenvectors keeps most of it
assert(info1.bytes > info0.bytes + 40e6, 'Eigen field looks too small for 2x10242x400 doubles');
t1 = tic; E2 = tess_eigen(Work, 'Laplace-Beltrami', 'nModes', 400); tReuse = toc(t1);
fprintf('reuse time: %.2f s\n', tReuse);
assert(tReuse < tSolve / 10, 'cache reuse should be much faster than solving');
assert(isequal(E2.Phi{1}, E.Phi{1}), 'reuse must return the stored basis');
% --- in_tess_bst passes the field through untouched ---
T = in_tess_bst(Work, 0);
assert(isfield(T, 'Eigen') && isfield(T.Eigen, 'LaplaceBeltrami'), 'in_tess_bst dropped Eigen');
delete(Work);
disp('test_tess_eigen_cortex PASSED');
