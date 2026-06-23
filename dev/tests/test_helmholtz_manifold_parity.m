function test_helmholtz_manifold_parity()
% TEST_HELMHOLTZ_MANIFOLD_PARITY: the canonical manifold geometry that bst_helmholtz consumes
% (per-face normal oriented outward via the divergence theorem + per-face area) must equal the
% globally-outward, consistently-wound mesh geometry. This is the correctness guard for the
% geometry swap: divergence (div = imag . n) and source/sink charge depend on the normal SIGN.
%
% NOTE: this DELIBERATELY does NOT require agreement with the old per-vertex VertNormals-flip
% normals. That heuristic mis-orients faces whose vertex normal is near-tangent to the face
% (thin sulci / high curvature) -- a latent ~1.6% bug the manifold normals fix. The invariant
% below (manifold == consistently-wound outward, exactly) is the geometrically correct one.
%
% Needs a loaded protocol whose Subject 1 has a cortex; SKIPs otherwise.
% Author: Diellor Basha, 2026
    if isempty(bst_get('ProtocolInfo'))
        disp('SKIP: test_helmholtz_manifold_parity -- no protocol loaded'); return;
    end
    sSubject = bst_get('Subject', 1);
    if isempty(sSubject) || isempty(sSubject.iCortex)
        disp('SKIP: test_helmholtz_manifold_parity -- Subject 1 has no cortex'); return;
    end
    SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
    Surf = in_tess_bst(SurfaceFile, 0);
    Mani = tess_manifold(SurfaceFile);
    Fcs  = double(Surf.Faces);
    nFail = 0;

    for hh = 1:numel(Mani.Embedded)
        E  = Mani.Embedded(hh);
        vH = double(E.GlobalVertices(:));  gf = double(E.GlobalFaces(:));
        mapV = zeros(size(Surf.Vertices,1),1); mapV(vH) = 1:numel(vH);
        Floc = mapV(Fcs(gf,:));  Vloc = Surf.Vertices(vH,:);
        % consistently-wound mesh normal, globally oriented outward (divergence theorem) -- ground truth
        e1 = Vloc(Floc(:,2),:)-Vloc(Floc(:,1),:);  e2 = Vloc(Floc(:,3),:)-Vloc(Floc(:,1),:);
        Nw = cross(e1,e2,2);  Aw = sqrt(sum(Nw.^2,2))/2;  Nw = Nw ./ max(2*Aw, eps);
        Cf = E.face.centroid;  c0 = mean(Cf,1);
        if sum(Aw .* sum(Nw.*(Cf-c0),2)) < 0, Nw = -Nw; end

        % manifold normal, oriented outward the same way bst_helmholtz does
        Nm = E.face.normal;  Am = E.face.area;
        if sum(Am .* sum(Nm.*(Cf-c0),2)) < 0, Nm = -Nm; end

        dotp = sum(Nm .* Nw, 2);
        relArea = max(abs(Aw - Am)) / max(Am);
        fprintf('  h%d  min(dot vs winding-outward)=%.6f  max relAreaDiff=%.2e\n', hh, min(dotp), relArea);
        nFail = nFail + chk(sprintf('h%d manifold normal == consistently-wound outward (min dot > 0.999999)',hh), min(dotp) > 0.999999);
        nFail = nFail + chk(sprintf('h%d face area matches winding cross (rel < 1e-9)',hh), relArea < 1e-9);
    end

    fprintf('\n==== test_helmholtz_manifold_parity: %d failed ====\n', nFail);
    if nFail > 0
        error(['test_helmholtz_manifold_parity FAILED -- oriented manifold normals diverge from ' ...
               'the consistently-wound outward mesh normals.']);
    end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
