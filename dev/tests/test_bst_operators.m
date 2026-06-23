function test_bst_operators()
% TEST_BST_OPERATORS: validate the foundational differential operators (toolbox/differential)
% with known-answer checks on the loaded protocol's cortex:
%   - gradient(constant)      == 0
%   - gradient(coordinate k)  ~ +tangential e_k  (true surface gradient, correct sign)
%   - laplacian(poisson(f))   ~ f - mean          (the two solvers invert, off the pinned vertex)
%   - helmholtz delegate Curl == bst_helmholtz('Decompose') Curl (bit-identical)
%
% Needs a loaded protocol whose Subject 1 has a cortex (+ nxr-compute); SKIPs otherwise.
% Author: Diellor Basha, 2026
    if isempty(bst_get('ProtocolInfo'))
        disp('SKIP: test_bst_operators -- no protocol loaded'); return;
    end
    sSubject = bst_get('Subject', 1);
    if isempty(sSubject) || isempty(sSubject.iCortex)
        disp('SKIP: test_bst_operators -- Subject 1 has no cortex'); return;
    end
    SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
    Surf = in_tess_bst(SurfaceFile, 0);
    nV   = size(Surf.Vertices, 1);
    F3   = double(Surf.Faces);
    % per-face winding normal (right-hand rule) -- the gradient is exact and tangent to each face
    e1f  = Surf.Vertices(F3(:,2),:) - Surf.Vertices(F3(:,1),:);
    e2f  = Surf.Vertices(F3(:,3),:) - Surf.Vertices(F3(:,1),:);
    Nf   = cross(e1f, e2f, 2);  Nf = Nf ./ max(sqrt(sum(Nf.^2,2)), eps);
    base = struct('SurfaceFile', SurfaceFile, 'iTargetStudy', 'NoSave');
    run  = @(m,f) getfield(subsref(bst_operators(f, setfield(base,'Method',m)), substruct('{}',{1})), 'Field');
    nFail = 0;

    % ----- gradient(constant) == 0 (per-face vector field [3nF]) -----
    Gc = run('gradient', ones(nV,1));
    nFail = nFail + chk('grad(const) == 0 (< 1e-9)', max(abs(Gc)) < 1e-9);

    % ----- gradient(coordinate k) == +tangential e_k, EXACTLY per face -----
    for k = 1:3
        G   = run('gradient', Surf.Vertices(:,k));
        Gf  = [G(1:3:end), G(2:3:end), G(3:3:end)];      % [nF x 3] per-face gradient
        ek  = zeros(1,3); ek(k) = 1;
        ekn = ek - Nf(:,k).*Nf;  ekn = ekn ./ max(sqrt(sum(ekn.^2,2)), eps);
        gmag = sqrt(sum(Gf.^2,2));  ghat = Gf ./ max(gmag, eps);  m = gmag > 0.1;
        cosang  = median(sum(ghat(m,:).*ekn(m,:), 2));
        normcmp = max(abs(sum(Gf(m,:).*Nf(m,:), 2)) ./ gmag(m));   % per-face gradient is IN the face plane
        nFail = nFail + chk(sprintf('grad(coord %d) == +tangential e (median cos > 0.999)', k), cosang > 0.999);
        nFail = nFail + chk(sprintf('grad(coord %d) lies in the face plane (max |nrm|/|g| < 1e-6)', k), normcmp < 1e-6);
    end

    % ----- laplacian(poisson(f)) ~ f - mean (per hemisphere, off the pinned vertex) -----
    f   = sin((1:nV)'*0.7) + cos((1:nV)'*0.13);
    g   = run('laplacian', run('poisson', f));
    LBO = bst_get_operator_node(SurfaceFile, 'Laplace-Beltrami');
    ref = zeros(nV,1);  keepMask = false(nV,1);
    for hh = 1:numel(LBO.GlobalVertices)
        vH = double(LBO.GlobalVertices{hh}(:));
        ref(vH) = f(vH) - mean(f(vH));
        keepMask(vH(2:end)) = true;        % exclude the pinned vertex 1 of each hemisphere
    end
    rho = corr(g(keepMask), ref(keepMask));
    nFail = nFail + chk('laplacian(poisson(f)) ~ f-mean (corr > 0.999)', rho > 0.999);

    % ----- de Rham complex: div/curl via DEC (grad is curl-free; rotated grad is div-free) -----
    ff  = sin((1:nV)'*0.04) + cos((1:nV)'*0.018);
    Mani= tess_manifold(SurfaceFile);
    gf  = run('gradient',   ff);                  % grad f       [3nF]
    dgf = run('divergence', gf);                  % div(grad f)  [nV]
    cgf = run('curl',       gf);                  % curl(grad f) [nV]
    lap = zeros(nV,1);                             % Laplace-Beltrami reference -M^-1 K f
    for hh = 1:numel(LBO.GlobalVertices)
        vH = double(LBO.GlobalVertices{hh}(:));  lap(vH) = -(LBO.Mass{hh}\(LBO.Operator{hh}*ff(vH)));
    end
    msk = dgf ~= 0;
    ccm = @(a,b) (a(msk)-mean(a(msk)))' * (b(msk)-mean(b(msk))) / (norm(a(msk)-mean(a(msk)))*norm(b(msk)-mean(b(msk))));
    nFail = nFail + chk('div(grad f) ~ Laplace-Beltrami (corr > 0.95)', ccm(dgf,lap) > 0.95);
    nFail = nFail + chk('curl(grad f) ~ 0 (rms ratio < 1e-3)', rms(cgf(msk))/rms(dgf(msk)) < 1e-3);
    % rotated gradient (N x grad f) is solenoidal: div ~ 0, curl ~ Laplacian
    vr = zeros(3*size(Surf.Faces,1), 1);
    for hh = 1:numel(Mani.DEC)
        fH = double(Mani.DEC(hh).GlobalFaces(:));
        N = Mani.Embedded(hh).face.normal;  C = Mani.Embedded(hh).face.centroid;  A = Mani.Embedded(hh).face.area;
        if sum(A.*sum(N.*(C-mean(C,1)),2)) < 0, N = -N; end
        gx=gf(3*fH-2); gy=gf(3*fH-1); gz=gf(3*fH);
        vr(3*fH-2)=N(:,2).*gz-N(:,3).*gy; vr(3*fH-1)=N(:,3).*gx-N(:,1).*gz; vr(3*fH)=N(:,1).*gy-N(:,2).*gx;
    end
    dvr = run('divergence', vr);  cvr = run('curl', vr);
    nFail = nFail + chk('div(rot grad f) ~ 0 (rms ratio < 1e-3)', rms(dvr(msk))/rms(dgf(msk)) < 1e-3);
    nFail = nFail + chk('curl(rot grad f) ~ Laplace-Beltrami (corr > 0.95)', ccm(cvr,lap) > 0.95);

    % ----- helmholtz delegate Curl == bst_helmholtz directly -----
    J    = sin((1:3*nV)'*0.31);
    Cop  = run('helmholtz', J);
    Mani = tess_manifold(SurfaceFile);
    Dir  = bst_get_operator_node(SurfaceFile, 'Dirac');
    Hd   = bst_helmholtz('Decompose', {Dir, LBO}, Mani, Surf, J);
    nFail = nFail + chk('helmholtz delegate Curl == bst_helmholtz (bit-identical)', isequal(Cop, Hd.Curl));

    fprintf('\n==== test_bst_operators: %d failed ====\n', nFail);
end

function fail = chk(name, cond)
    if cond, fprintf('  PASS %s\n', name); fail = 0;
    else,    fprintf('  FAIL %s\n', name); fail = 1; end
end
