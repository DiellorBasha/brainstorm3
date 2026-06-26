function test_nxr_v020_smoke()
% Backward-compat gate for nxr-compute v0.2.0: every command/operator our
% Brainstorm consumers depend on must still return. Mesh prep mirrors
% tess_operators (tess_hemisplit -> single-hemisphere submesh); never pass a
% hand-built or whole-brain mesh to nxr_compute('create'). No 'clear' (live-safe).

    ver = nxr_compute('version');
    fprintf('nxr_compute version: %s\n', ver);
    assert(~isempty(ver), 'nxr_compute(''version'') returned empty');

    % Canonical cortex from the loaded protocol; build ONE hemisphere submesh.
    SurfaceFile = bst_canonical_cortex(20484);
    [Vloc, Floc] = i_hemi_submesh(SurfaceFile, 1);
    h = nxr_compute('create', Vloc, Floc);   % validates the mesh (clean error, no segfault)
    cleanup = onCleanup(@() nxr_compute('destroy', h)); %#ok<NASGU>

    chk = @(name, M) assert(~isempty(M) && all(size(M) > 0), ...
        sprintf('operator %s returned empty/degenerate', name));

    chk('laplacian/cotan',       nxr_compute('operators', h, 'laplacian', 'cotan'));
    chk('mass/galerkin',         nxr_compute('operators', h, 'mass', 'galerkin'));
    chk('mass/lumped',           nxr_compute('operators', h, 'mass', 'lumped'));
    chk('laplacian/connection',  nxr_compute('operators', h, 'laplacian', 'connection'));
    chk('dirac',                 nxr_compute('operators', h, 'dirac', 1));
    chk('diracD',                nxr_compute('operators', h, 'diracD'));
    chk('diracIntrinsicD',       nxr_compute('operators', h, 'diracIntrinsicD'));
    chk('diracFace',             nxr_compute('operators', h, 'diracFace', 1));
    chk('diracFaceIntrinsicD',   nxr_compute('operators', h, 'diracFaceIntrinsicD'));
    chk('gradFace',              nxr_compute('operators', h, 'gradFace'));
    chk('lapFace',               nxr_compute('operators', h, 'lapFace'));
    dec  = nxr_compute('operators', h, 'dec');               assert(isfield(dec,'d0') && isfield(dec,'d1'), 'dec missing d0/d1');
    chk('hodge/h0',              nxr_compute('operators', h, 'hodge', 'h0'));
    chk('hodge/h1',              nxr_compute('operators', h, 'hodge', 'h1'));
    chk('hodge/h2',              nxr_compute('operators', h, 'hodge', 'h2'));
    gz = nxr_compute('gauge', h, 'levi-civita', struct('operators',true,'coupling','ambient'));
    assert(isfield(gz,'operators') && isfield(gz.operators,'covariantLaplacian'), 'gauge/levi-civita missing covariantLaplacian');
    vf = nxr_compute('vertexFrames', h);                     assert(isfield(vf,'e1') && isfield(vf,'normals'), 'vertexFrames missing fields');

    % NEW in v0.2.0: registry introspection
    oi = nxr_compute('operatorInfo', 'laplaceBeltrami');
    assert(isstruct(oi) && strcmp(oi.id, 'laplaceBeltrami'), 'operatorInfo(laplaceBeltrami) failed');

    fprintf('test_nxr_v020_smoke: ALL PASS\n');
end

function [Vloc, Floc] = i_hemi_submesh(SurfaceFile, hh)
% Build hemisphere-hh local submesh exactly as tess_operators does (tess_hemisplit).
    TessMat = in_tess_bst(SurfaceFile, 0);
    [ir, il] = tess_hemisplit(TessMat);
    hemis = {il, ir};  vH = double(hemis{hh}(:));
    Vtx = TessMat.Vertices; Fcs = double(TessMat.Faces); nVtot = size(Vtx,1);
    isV = false(nVtot,1); isV(vH)=true; fMask = all(isV(Fcs),2);
    mapV = zeros(nVtot,1); mapV(vH)=1:numel(vH);
    Vloc = Vtx(vH,:); Floc = mapV(Fcs(fMask,:));
end
