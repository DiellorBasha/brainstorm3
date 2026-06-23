function test_helmholtz_parity
% Re-run bst_helmholtz on the saved input and compare to the pre-refactor baseline.
    fprintf('== test_helmholtz_parity ==\n');
    S = load('dev/scratch/helmholtz_baseline.mat');   % H0, J, Surf
    Surfm = in_tess_bst(S.Surf,0);
    Mani  = tess_manifold(S.Surf,'Gauge','trivial');
    Dir   = bst_get_operator_node(S.Surf,'Dirac');  LBO = bst_get_operator_node(S.Surf,'Laplace-Beltrami');
    H1 = bst_helmholtz('Decompose', {Dir, LBO}, Mani, Surfm, S.J);
    flds = {'Curl','Div','Psi','Phi','Fmag'};
    for i=1:numel(flds)
        a=S.H0.(flds{i}); b=H1.(flds{i});
        e=norm(a(:)-b(:))/max(norm(a(:)),eps);
        assert(e<1e-10, '%s parity broken: rel err %g', flds{i}, e);
        fprintf('  %-5s rel err %g  [OK]\n', flds{i}, e);
    end
    fprintf('PASS\n');
end
