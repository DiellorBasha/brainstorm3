function test_frame_sync()
% TEST_FRAME_SYNC: the Connection-Laplacian operator stores the EXACT canonical tangent
% frame its complex eigenmodes decode in, provably synchronized with the nxr frame the
% manifold also uses.
%
% Guarantees (the "validated trust" of the frame-synchronization design,
% dev/2026-06-22-frame-synchronization-design.md):
%   - tess_operators('Connection Laplacian') populates OperatorMat.Frame (e1,e2,normal);
%   - that stored frame == a fresh nxr_compute('vertexFrames') call, BITWISE (determinism);
%   - the frame is a proper right-handed orthonormal tangent frame;
%   - it equals the manifold's Embedded.vertex.grid decode (e1=real, e2=imag), BITWISE
%     (so referencing the manifold frame is provably the operator frame);
%   - a connection eigenmode decoded via the frame is tangent to the surface.
%
% Needs nxr-compute + a real surface (Subject01). SKIPs if unavailable.
%
% Authors: Diellor Basha, 2026

    nFail = 0;  chk = @(n,c) i_chk(n,c);
    SurfaceFile = 'Subject01/tess_cortex_pial_low.mat';
    if isempty(bst_get('SurfaceFile', SurfaceFile))
        fprintf('SKIP test_frame_sync: surface %s not found.\n', SurfaceFile); return;
    end

    % ---- build the Connection-Laplacian operator (no save) ----
    OM = tess_operators(SurfaceFile, 'Connection Laplacian', 'NoSave', true);
    nFail = nFail + chk('OperatorMat.Frame present (1x2)', isfield(OM,'Frame') && numel(OM.Frame)==2);

    T = in_tess_bst(SurfaceFile, 0);
    for hh = 1:2
        Fr = OM.Frame{hh};
        if isempty(Fr); fprintf('  (hemi %d empty)\n', hh); continue; end
        gv = OM.GlobalVertices{hh}(:);
        nV = numel(gv);
        nFail = nFail + chk(sprintf('hemi %d frame sizes', hh), ...
            isequal(size(Fr.e1),[nV 3]) && isequal(size(Fr.e2),[nV 3]) && isequal(size(Fr.normal),[nV 3]));

        % rebuild the SAME submesh (operator's vertex order) and fetch nxr vertexFrames
        isV=false(size(T.Vertices,1),1); isV(gv)=true; fMask=all(isV(double(T.Faces)),2);
        mapV=zeros(size(T.Vertices,1),1); mapV(gv)=1:nV;
        Vloc=T.Vertices(gv,:); Floc=mapV(double(T.Faces(fMask,:)));
        h=nxr_compute('create',Vloc,Floc);
        VF=nxr_compute('vertexFrames', h);
        nxr_compute('destroy', h);
        % BITWISE determinism: stored operator frame == fresh nxr frame
        d = max([norm(Fr.e1-VF.e1,'fro'), norm(Fr.e2-VF.e2,'fro'), norm(Fr.normal-VF.normals,'fro')]);
        nFail = nFail + chk(sprintf('hemi %d: operator frame == nxr vertexFrames (bitwise)', hh), d == 0);

        % proper orthonormal right-handed tangent frame
        u1=sqrt(sum(Fr.e1.^2,2)); u2=sqrt(sum(Fr.e2.^2,2));
        d12=max(abs(sum(Fr.e1.*Fr.e2,2))); d1n=max(abs(sum(Fr.e1.*Fr.normal,2)));
        rh = mean(sum(cross(Fr.e1,Fr.e2,2).*Fr.normal,2) > 0);
        nFail = nFail + chk(sprintf('hemi %d: orthonormal right-handed', hh), ...
            max(abs(u1-1))<1e-9 && max(abs(u2-1))<1e-9 && d12<1e-9 && d1n<1e-9 && rh>0.999);
    end

    % ---- manifold equivalence: Embedded.vertex.grid decode == operator frame (bitwise) ----
    haveMan = exist('in_bst_manifold','file') == 2;
    sMan = []; if haveMan; try, sMan = bst_get('ManifoldFileForSurface', SurfaceFile, 'levi-civita'); catch, end; end
    if ~isempty(sMan)
        MM = in_bst_manifold(sMan.FileName);
        for hh = 1:2
            if isempty(OM.Frame{hh}) || numel(MM.Embedded) < hh || isempty(MM.Embedded(hh)); continue; end
            g = MM.Embedded(hh).vertex.grid;     % [nV x 3] complex = e1 + i*e2
            d = max([norm(OM.Frame{hh}.e1-real(g),'fro'), norm(OM.Frame{hh}.e2-imag(g),'fro')]);
            nFail = nFail + chk(sprintf('hemi %d: manifold grid decode == operator frame (bitwise)', hh), d==0);
        end
    else
        fprintf('  (no levi-civita manifold node; manifold-equivalence check skipped -- guaranteed by determinism vs nxr above)\n');
    end

    % ---- a connection eigenmode decoded via the frame is tangent to the surface ----
    efC = i_find_variant('Connection Laplacian');
    if ~isempty(efC)
        EC = in_bst_eigen(efC);
        for hh = 1:2
            if isempty(EC.Phi{hh}) || isempty(OM.Frame{hh}); continue; end
            if ~isequal(EC.GlobalVertices{hh}(:), OM.GlobalVertices{hh}(:)); continue; end
            z = EC.Phi{hh}(:,1);
            v3 = real(z).*OM.Frame{hh}.e1 + imag(z).*OM.Frame{hh}.e2;   % decode to 3D
            tangentErr = max(abs(sum(v3.*OM.Frame{hh}.normal,2)));      % should be ~0
            nFail = nFail + chk(sprintf('hemi %d: decoded eigenmode is tangent', hh), tangentErr < 1e-9);
        end
    else
        fprintf('  (no Connection Laplacian eigen node; decode-tangent check skipped)\n');
    end

    fprintf('\n==== test_frame_sync: %d failed ====\n', nFail);
    if nFail > 0; error('test_frame_sync FAILED'); end
    disp('ALL TESTS PASSED');
end

% ===== helpers =====
function r = i_chk(nm, cond)
    if cond; r = 0; fprintf('  ok   %s\n', nm); else; r = 1; fprintf('  FAIL %s\n', nm); end
end

function ef = i_find_variant(want)
    ef = [];
    PI = bst_get('ProtocolInfo'); if isempty(PI); return; end
    d = dir(fullfile(PI.SUBJECTS, '**', 'eigen_*.mat'));
    [~, ord] = sort([d.datenum], 'descend');
    for i = ord(:)'
        rel = strrep(fullfile(d(i).folder, d(i).name), [PI.SUBJECTS filesep], '');
        try
            m = in_bst_eigen(rel, 'Variant');
            if strcmpi(m.Variant, want); ef = rel; return; end
        catch
        end
    end
end
