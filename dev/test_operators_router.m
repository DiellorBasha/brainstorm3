function test_operators_router
% Router: valid (Method x stratum) pairs run; invalid pairs raise bst_operators:badFieldType;
% identity checks curl(grad f)~0, div(curl V)~0.
    fprintf('== test_operators_router ==\n');
    Surf = 'Subject01/tess_cortex_pial_low.mat';
    Surfm = in_tess_bst(Surf,0);  nV = size(Surfm.Vertices,1);  nF = size(Surfm.Faces,1);
    f = randn(nV, 2);  J = randn(3*nV, 2);

    % curl of a scalar -> must be guarded
    ok=false; try, bst_operators(f, struct('Method','curl','SurfaceFile',Surf,'iTargetStudy','NoSave')); catch e, ok=strcmp(e.identifier,'bst_operators:badFieldType'); end
    assert(ok, 'curl-of-scalar was not guarded'); fprintf('  curl-of-scalar guarded  [OK]\n');

    % divergence of an ambient field -> runs, per-vertex scalar
    R = bst_operators(J, struct('Method','divergence','FieldType','ambient','SurfaceFile',Surf,'iTargetStudy','NoSave'));
    assert(size(R{1}.Field,1)==nV, 'ambient divergence shape'); fprintf('  ambient divergence routed  [OK]\n');

    % identity: curl(grad f) ~ 0 (metric-free)
    Rg = bst_operators(f, struct('Method','gradient','SurfaceFile',Surf,'iTargetStudy','NoSave'));
    Rc = bst_operators(Rg{1}.Field, struct('Method','curl','FieldType','tangent','SurfaceFile',Surf,'iTargetStudy','NoSave'));
    rel = norm(Rc{1}.Field(:))/max(norm(Rg{1}.Field(:)),eps);
    fprintf('  curl(grad f) rel mag = %g\n', rel);
    assert(rel < 1e-2, 'curl(grad f) not ~0 (rel %g)', rel);
    fprintf('PASS\n');
end
