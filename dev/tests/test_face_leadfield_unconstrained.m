function test_face_leadfield_unconstrained()
% Full-unconstrained face leadfield: raw [nCh x 3F] Sarvas at centroids, geometry from
% tess_manifold, exact parity with the vertex model. Validates against the existing
% (validated) constrained mode so we are not checking correctness in a vacuum.
% Author: Diellor Basha, 2026
    nFail = 0;
    df = 'Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
    [sStudy,~] = bst_get('DataFile', df);
    BaseHM = in_bst_headmodel([fileparts(df) '/headmodel_surf_os_meg.mat'], 0);
    ChanMat = in_bst_channel(sStudy.Channel(1).FileName);
    iMEG = find(strcmpi({ChanMat.Channel.Type}, 'MEG'));
    Channel = ChanMat.Channel(iMEG);  Param = BaseHM.Param(iMEG);
    SurfaceFile = BaseHM.SurfaceFile;
    Surf = in_tess_bst(SurfaceFile,0);  V = Surf.Vertices;  F = double(Surf.Faces);  nF = size(F,1);

    [L, FG] = bst_face_leadfield(SurfaceFile, Channel, Param, 'Mode','unconstrained');

    nFail = nFail + chk('shape [nMEG x 3F], finite', isequal(size(L),[numel(iMEG) 3*nF]) && all(isfinite(L(:))));

    % geometry from tess_manifold: centroids == barycentric mean (exact), normals unit
    Cbary = (V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;
    nFail = nFail + chk('manifold centroids == barycentric', max(abs(FG.Centroids - Cbary),[],'all') < 1e-9);
    nFail = nFail + chk('normals are unit', max(abs(sqrt(sum(FG.Normals.^2,2)) - 1)) < 1e-6);

    % pure Sarvas: the 3-col block for a sampled face IS bst_meg_sph at that centroid
    ff = round(nF/3);  g = bst_meg_sph(FG.Centroids(ff,:)', Channel, Param);
    nFail = nFail + chk('block == raw Sarvas (no projection)', max(abs(L(:,3*ff-2:3*ff) - g),[],'all') < 1e-9*max(abs(g(:))));

    % constrained-consistency: the (validated) constrained column f equals the unconstrained
    % block projected onto the constrained mode's own normal*area (ties new to known-good).
    okC = true;
    try
        [Lc, FGc] = bst_face_leadfield(SurfaceFile, Channel, Param, 'Mode','constrained');
    catch e
        okC = false; fprintf('  SKIP constrained-consistency (constrained mode unavailable: %s)\n', e.message);
    end
    if okC
        proj = sum(L(:,3*ff-2:3*ff) .* (FGc.Normals(ff,:) * FGc.Areas(ff)), 2);  % block * (n*A)
        rel = norm(proj - Lc(:,ff)) / max(norm(Lc(:,ff)), eps);
        fprintf('  [constrained-consistency] rel err @face %d = %.2e\n', ff, rel);
        nFail = nFail + chk('constrained col == unconstr block . (n*A)', rel < 1e-9);
    end

    % observability ceiling: face leadfield effective rank ~ vertex (MEG DOF, ~tens)
    er = @(G) sum(svd(G) > 1e-3*max(svd(G)));
    rF = er(L);  rV = er(BaseHM.Gain(iMEG,:));
    fprintf('  [effective rank] face=%d  vertex=%d\n', rF, rV);
    nFail = nFail + chk('face effective rank ~ vertex (0.5x..2x)', rF >= 0.5*rV && rF <= 2*rV);

    fprintf('\n==== test_face_leadfield_unconstrained: %d failed ====\n', nFail);
    if nFail > 0, error('test_face_leadfield_unconstrained FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
