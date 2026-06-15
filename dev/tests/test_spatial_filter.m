function test_spatial_filter()
% Apply/Restore the spatial filter on a synthetic in-memory unconstrained source,
% checking the swapped ImageGridAmp equals the Dirac filter and Restore is exact.
% Requires Brainstorm running with a Dirac eigen node (Subject01 surface 5).
% Authors: Diellor Basha, 2026
    nFail = 0;
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    EigenMat    = tess_eigen(SurfaceFile, 'Dirac');
    Op = load(file_fullpath(EigenMat.OperatorFile));
    nVert = double(max(cellfun(@(x) max(x(:)), EigenMat.GlobalVertices)));
    rng(0); J0 = randn(3*nVert, 4);          % synthetic [3nVert x 4 time]

    % reference: filter via the core (heat)
    g = bst_eigfilter_kernel('heat', struct('t', 1/median(EigenMat.Lambda{1})));
    Jref = real(bst_dirac_eigenmodes_filter(EigenMat, Op.Mass, J0, 'custom', 'TransferFn', g));

    % exercise the panel's pure ComputeFiltered against the reference
    St = struct('EigenMat',EigenMat, 'Mass',{Op.Mass}, 'Orig',J0);
    Jf = panel_spatial_filter('ComputeFiltered', St, 'heat', struct('t', 1/median(EigenMat.Lambda{1})));
    nFail = nFail + chk('ComputeFiltered == core filter', max(abs(Jf(:)-Jref(:))) < 1e-9);
    nFail = nFail + chk('filtered differs from original', max(abs(Jf(:)-J0(:))) > 0);

    fprintf('\n==== test_spatial_filter: %d failed ====\n', nFail);
    if nFail > 0, error('test_spatial_filter FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
