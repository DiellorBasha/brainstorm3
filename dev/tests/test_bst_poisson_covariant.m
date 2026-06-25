function test_bst_poisson_covariant
    SurfaceFile = 'Subject01/tess_cortex_pial_low.mat';
    Cov = tess_operators(SurfaceFile, 'Covariant');
    nV  = max(cellfun(@(c) max(double(c(:))), Cov.GlobalVertices));
    rng(3);  f = randn(nV, 1);
    phi = bst_poisson(Cov, f);                       % must not error on 'Covariant'
    assert(isequal(size(phi), [nV 1]), 'phi shape wrong');
    % K phi = M f residual on the free block, per hemisphere
    for hh = 1:numel(Cov.Operator)
        vH = double(Cov.GlobalVertices{hh}(:));  K = (Cov.Operator{hh}+Cov.Operator{hh}')/2;
        M = Cov.Mass{hh};  free = 2:numel(vH);  ph = phi(vH);  fh = f(vH);
        fh = fh - sum(M*fh)/sum(M(:));  r = K(free,free)*ph(free) - M(free,:)*fh;
        assert(norm(r)/max(norm(M(free,:)*fh),eps) < 1e-8, 'Poisson residual too large (hemi %d)', hh);
    end
    fprintf('PASS test_bst_poisson_covariant\n');
end
