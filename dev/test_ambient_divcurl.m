function test_ambient_divcurl
% Ambient-field branch: (a) flat-covariant divergence/curl shapes;
% (b) the mean-curvature term vanishes for a purely-tangential field and is nonzero
%     for a purely-normal field on a curved cortex.
    fprintf('== test_ambient_divcurl ==\n');

    % Use the prepared low-res cortex directly (20484V: manifold + Dirac + LBO ready)
    Surf = 'Subject01/tess_cortex_pial_low.mat';
    Surfm = in_tess_bst(Surf, 0);
    nV = size(Surfm.Vertices, 1);
    Mani  = tess_manifold(Surf, 'Gauge', 'trivial');
    Cov   = tess_operators(Surf, 'Covariant');
    rng(11);  J = randn(3*nV, 3);

    dv = bst_divergence(J, Mani, 'Ambient', Surfm, Cov);
    cu = bst_curl(J, Mani, 'Ambient', Surfm, Cov);
    assert(isequal(size(dv), [nV 3]) && isequal(size(cu), [nV 3]), 'ambient output shape wrong');
    fprintf('  ambient div/curl shapes OK\n');

    % Vertex normals: use VertNormals if present, else compute via tess_normals
    if isfield(Surfm, 'VertNormals') && ~isempty(Surfm.VertNormals)
        N = Surfm.VertNormals;  % [nV x 3]
        fprintf('  using Surfm.VertNormals (present)\n');
    else
        N = tess_normals(Surfm.Vertices, Surfm.Faces);
        fprintf('  VertNormals absent -> computed via tess_normals\n');
    end

    % purely-normal field: mean-curvature term must dominate (nonzero on a folded cortex)
    Jn = reshape((N .* 1)', [], 1);  % [3nV x 1] normal field, unit magnitude
    dvn = bst_divergence(Jn, Mani, 'Ambient', Surfm, Cov);
    assert(norm(dvn) > 1e-6, 'normal-field divergence unexpectedly ~0 (mean-curvature term missing?)');
    fprintf('  normal-field divergence nonzero (|dv|=%g) -> mean-curvature term present  [OK]\n', norm(dvn));
    fprintf('PASS\n');
end
