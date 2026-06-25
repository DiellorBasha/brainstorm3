function test_bst_curl_ambient
    d = load(fullfile(fileparts(mfilename('fullpath')), 'baselines', 'helmholtz_baseline.mat'));
    B = d.B;  Cov = tess_operators(B.Surf, 'Covariant');
    cu = bst_curl(B.J, [], 'Ambient', [], Cov);
    assert(isequal(size(cu), size(B.Curl)), 'curl shape mismatch');
    rel = norm(cu - B.Curl) / max(norm(B.Curl), eps);
    assert(rel < 1e-10, 'ambient vorticity differs from baseline (rel=%.2e)', rel);
    fprintf('PASS test_bst_curl_ambient (rel=%.2e)\n', rel);
end
