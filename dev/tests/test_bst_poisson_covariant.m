function test_bst_poisson_covariant
    SurfaceFile = 'Subject01/tess_cortex_pial_low.mat';
    Cov = tess_operators(SurfaceFile, 'Covariant');
    nV  = max(cellfun(@(c) max(double(c(:))), Cov.GlobalVertices));
    rng(3);  f = randn(nV, 1);
    phi = bst_poisson(Cov, f);                       % must not error on 'Covariant'
    assert(isequal(size(phi), [nV 1]), 'phi shape wrong');
    for hh = 1:numel(Cov.Operator)
        vH = double(Cov.GlobalVertices{hh}(:));
        K  = (Cov.Operator{hh}+Cov.Operator{hh}')/2;  M = Cov.Mass{hh};
        ph = phi(vH);  fh = f(vH);  fh = fh - sum(M*fh)/sum(M(:));   % mean-zero RHS (compatibility)
        % Full-stiffness residual is gauge-invariant (K*1=0); the free-block residual is NOT.
        rel = norm(K*ph - M*fh) / max(norm(M*fh), eps);
        assert(rel < 1e-8, 'Poisson residual too large (hemi %d, rel=%.2e)', hh, rel);
        assert(abs(mean(ph)) < 1e-9, 'bst_poisson must return the mean-zero gauge (hemi %d)', hh);
    end
    fprintf('PASS test_bst_poisson_covariant\n');
end
