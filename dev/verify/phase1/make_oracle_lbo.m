% MAKE_ORACLE_LBO: nxr-built LBO pencil for the parity surface (dev side, run ONCE).
SurfaceFile = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/verify/phase0/bst_userdir_clean/.brainstorm/local_db/omega-tutorial-cortical-flow/anat/sub-0002/tess_cortex_pial_low.mat';
OutFile = fullfile(fileparts(mfilename('fullpath')), 'oracle_lbo_sub0002.mat');
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute unavailable: %s', errMsg);
TessMat = load(SurfaceFile);
[rH, lH, isConn] = tess_hemisplit(TessMat);
assert(~isConn, 'hemispheres connected');
hemis = {lH(:), rH(:)};
Vtx = double(TessMat.Vertices); Fcs = double(TessMat.Faces);
nVtot = size(Vtx, 1);
A = cell(1,2); B = cell(1,2); vH = cell(1,2);
for hh = 1:2
    v = hemis{hh};
    isV = false(nVtot,1); isV(v) = true;
    fMask = all(isV(Fcs), 2);
    mapV = zeros(nVtot,1); mapV(v) = 1:numel(v);
    h = nxr_compute('create', Vtx(v,:), mapV(Fcs(fMask,:)));
    A{hh} = nxr_compute('operators', h, 'laplacian', 'cotan');
    B{hh} = nxr_compute('operators', h, 'mass', 'galerkin');
    nxr_compute('destroy', h);
    vH{hh} = v;
end
nxrVer = ''; try, nxrVer = nxr_compute('version'); catch, end %#ok<CTCH>
meta = struct('SurfaceFile', SurfaceFile, 'nVertices', nVtot, ...
              'NxrVersion', nxrVer, 'Date', datestr(now), 'conventions', '');
save(OutFile, 'A', 'B', 'vH', 'meta', '-v7');
fprintf('ORACLE WRITTEN: %s  (nV=%d, L=%d, R=%d)\n', OutFile, nVtot, numel(vH{1}), numel(vH{2}));
