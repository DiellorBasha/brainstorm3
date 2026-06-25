function test_process_helmholtz_compute
    d = load(fullfile(fileparts(mfilename('fullpath')), 'baselines', 'helmholtz_baseline.mat'));
    B = d.B;  Cov = tess_operators(B.Surf, 'Covariant');
    Ht = process_helmholtz('Compute', B.J, Cov);
    chk = @(f) norm(Ht.(f)(:) - B.(f)(:)) / max(norm(B.(f)(:)), eps);
    for f = {'Div','Curl','Phi','Psi','Fmag','Hmag','Virr','Vsol','Vtot','Hresid'}
        rel = chk(f{1});
        assert(rel < 1e-9, '%s differs from baseline (rel=%.2e)', f{1}, rel);
    end
    assert(abs(Ht.HarmFrac - B.HarmFrac) < 1e-9*max(abs(B.HarmFrac),eps), 'HarmFrac differs');
    fprintf('PASS test_process_helmholtz_compute\n');
end
