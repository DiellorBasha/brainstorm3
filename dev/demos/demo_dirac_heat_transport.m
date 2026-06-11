function R = demo_dirac_heat_transport(EigenFile, SurfaceFile, varargin)
% DEMO_DIRAC_HEAT_TRANSPORT  Demonstrate that a field "constant in the transported
% frame" is a LOW Dirac eigenmode, and that the operator carries the across-sulcus
% frame rotation (ambient sign flip) as geometry -- not as high-frequency content.
%
% Method (the user's suggestion): a heat filter on the Dirac eigenvalue axis applied
% to a single-vertex VECTOR impulse, reconstructed back to a per-vertex 3-vector field:
%       f_t = Phi * ( exp(-t*lambda) .* (Phi' * B * delta) )
% This is the Dirac heat kernel (the "vector heat method" object). We place the impulse
% at a sulcal fundus, pointing ACROSS the fold (tangentially), and show:
%   (1) the reconstructed field is SMOOTH (a single low-pass blob), and its spectral
%       energy collapses onto the lowest modes  -> "smooth in transported frame = low mode";
%   (2) across the sulcus the AMBIENT direction flips ~180 deg, yet the field is one
%       continuous smooth section -> the operator encodes the frame rotation "exactly once";
%   (3) CONTRAST: diffusing the same impulse's 3 ambient components with a naive SCALAR
%       cotan heat kernel does NOT flip -> the gauge-fixed approach ignores the fold.
%
% USAGE:
%   R = demo_dirac_heat_transport('Subject01/eigen_260610_2304.mat', ...
%                                 'Subject01/tess_cortex_pial_low.mat')
%   ...,'Hemi',1, 'ModeCut',120, 'Savefig','/abs/path.png'
%
% Authors: Diellor Basha, 2026

    % ---- options ----
    o = struct('Hemi',1, 'ModeCut',120, 'Sach',[], 'Savefig', ...
        '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks/demo_dirac_heat_transport.png');
    for i=1:2:numel(varargin), o.(varargin{i}) = varargin{i+1}; end
    hh = o.Hemi;

    % ---- load eigenbasis + operator mass + surface ----
    E   = load(file_fullpath(EigenFile));
    Op  = load(file_fullpath(E.OperatorFile));
    Srf = in_tess_bst(SurfaceFile);
    assert(strcmpi(E.Variant,'Dirac'), 'Eigen node is not Dirac.');

    gv  = E.GlobalVertices{hh}(:);          % hemi-local -> global vertex index
    Vh  = numel(gv);
    Phi = double(E.Phi{hh});                % [4Vh x K]
    lam = double(E.Lambda{hh}(:));          % [K x 1]
    B   = Op.Mass{hh};                      % [4Vh x 4Vh] = kron(Mg, I4)
    K   = numel(lam);
    X   = Srf.Vertices(gv,:);               % [Vh x 3] hemi vertex positions (meters)

    % ---- hemi submesh (local faces) for normals + scalar cotan kernel ----
    g2l = zeros(size(Srf.Vertices,1),1); g2l(gv) = 1:Vh;
    Fg  = Srf.Faces;
    isL = all(g2l(Fg) > 0, 2);
    Floc = g2l(Fg(isL,:));                  % [nF x 3] local
    Nrm  = local_vertex_normals(X, Floc);   % [Vh x 3] unit (outward sign arbitrary)

    % ---- pick a sulcal fundus pair (v0 on wall A, w0 on opposing wall B) ----
    [v0, w0, gap] = local_pick_sulcus(X, Nrm);
    fprintf('Sulcus pair: v0=%d  w0=%d  gap=%.2f mm  n.n=%.3f\n', ...
        v0, w0, gap*1e3, dot(Nrm(v0,:),Nrm(w0,:)));

    % ---- impulse: unit ACROSS-FOLD tangent at v0 ----
    p  = X(w0,:) - X(v0,:);
    d  = p - dot(p,Nrm(v0,:))*Nrm(v0,:);    % project onto v0 tangent plane
    d  = d / norm(d);                       % across-fold tangent direction (ambient)
    delta = zeros(4*Vh,1);
    delta(4*(v0-1)+2) = d(1);               % i (x)
    delta(4*(v0-1)+3) = d(2);               % j (y)
    delta(4*(v0-1)+4) = d(3);               % k (z)   (w slot 1 stays 0)

    % ---- Dirac heat filter on the eigenvalue axis ----
    c  = Phi' * (B * delta);                % [K x 1] spectral coefficients
    % heat time: half-max of the filter at mode index ModeCut (skip ~0 kernel modes)
    kc = min(o.ModeCut, K);
    lref = lam(kc); if lref <= 0, lref = max(lam(lam>0)); end
    t  = log(2) / lref;
    w  = exp(-t * lam);                     % [K x 1] heat weights
    fq = Phi * (w .* c);                    % [4Vh x 1] reconstructed quaternion field

    % per-vertex 3-vector (drop w slot) + scalar-leak diagnostic
    Q  = reshape(fq, 4, Vh).';              % [Vh x 4] = [w x y z]
    V3 = Q(:,2:4);                          % [Vh x 3] vector field (ambient)
    mag = sqrt(sum(V3.^2,2));
    scalarLeak = norm(Q(:,1)) / norm(V3(:));

    % blob radius (magnitude-weighted mean distance from v0)
    dist0 = sqrt(sum((X - X(v0,:)).^2,2));
    Rd = sum(mag.*dist0) / sum(mag);

    % ---- spectral concentration: energy of filtered coeffs vs raw ----
    eRaw = abs(c).^2;       eRaw = eRaw / sum(eRaw);
    eFil = abs(w.*c).^2;    eFil = eFil / sum(eFil);
    csumFil = cumsum(eFil);
    nLow = find(csumFil >= 0.90, 1, 'first');   % # modes holding 90% of filtered energy

    % ---- NAIVE baseline: component-wise scalar cotan heat kernel, extent-matched ----
    [Lc, Ml] = local_cotmatrix_mass(X, Floc);
    G0 = zeros(Vh,3); G0(v0,:) = d;             % same ambient impulse
    % match spatial extent Rd by searching the implicit-Euler time
    ts_grid = Rd^2 * logspace(-1.0, 1.2, 14);
    bestErr = inf; Vn = []; ts_best = ts_grid(1); Rn_best = NaN;
    for ts = ts_grid
        Gd = (Ml + ts*Lc) \ (Ml*G0);            % [Vh x 3] diffuse 3 ambient comps
        mn = sqrt(sum(Gd.^2,2));
        Rn = sum(mn.*dist0)/sum(mn);
        err = abs(Rn - Rd);
        if err < bestErr, bestErr = err; Vn = Gd; ts_best = ts; Rn_best = Rn; end
    end
    magN = sqrt(sum(Vn.^2,2));

    % ---- metrics at the across-sulcus partner ----
    uDir = @(M,i) M(i,:)/max(norm(M(i,:)),eps);
    angDirac = local_angle(uDir(V3,v0), uDir(V3,w0));
    angNaive = local_angle(uDir(Vn,v0), uDir(Vn,w0));
    % normal-leak: how tangent does each field stay across the blob (weighted |dir.n|)
    blob = mag > 0.15*max(mag);
    leakDirac = sum(mag(blob).*abs(sum(local_unit(V3(blob,:)).*Nrm(blob,:),2)))/sum(mag(blob));
    leakNaive = sum(magN(blob).*abs(sum(local_unit(Vn(blob,:)).*Nrm(blob,:),2)))/sum(magN(blob));

    % ================= report =================
    fprintf('\n==== Dirac heat transport demo (hemi %d, %d vertices, K=%d) ====\n', hh, Vh, K);
    fprintf('Heat time t=%.3g (half-max at mode %d, lambda=%.3g)\n', t, kc, lref);
    fprintf('Blob radius: Dirac=%.2f mm | naive(matched)=%.2f mm (ts=%.3g)\n', Rd*1e3, Rn_best*1e3, ts_best);
    fprintf('SMOOTHNESS  : 90%% of filtered energy in %d / %d modes (raw needed %d)\n', ...
        nLow, K, find(cumsum(eRaw)>=0.90,1,'first'));
    fprintf('SCALAR LEAK : ||w-slot||/||vec|| = %.3f (curvature coupling into scalar part)\n', scalarLeak);
    fprintf('--- ACROSS-SULCUS AMBIENT ANGLE (v0 -> w0) ---\n');
    fprintf('   DIRAC  : %.1f deg   <-- operator carries the frame flip (geometry)\n', angDirac);
    fprintf('   NAIVE  : %.1f deg   <-- scalar diffusion keeps fixed ambient dir (wrong)\n', angNaive);
    fprintf('--- TANGENT FIDELITY (weighted |dir . n|, 0=tangent) ---\n');
    fprintf('   DIRAC  : %.3f      NAIVE : %.3f\n', leakDirac, leakNaive);
    fprintf('Magnitude at partner w0: Dirac=%.3f  naive=%.3f (rel to peak)\n', ...
        mag(w0)/max(mag), magN(w0)/max(magN));

    R = struct('v0',v0,'w0',w0,'gap',gap,'d',d,'t',t,'Rd',Rd,'Rn',Rn_best, ...
        'angDirac',angDirac,'angNaive',angNaive,'leakDirac',leakDirac,'leakNaive',leakNaive, ...
        'scalarLeak',scalarLeak,'nLow',nLow,'V3',V3,'Vn',Vn,'mag',mag,'X',X,'Nrm',Nrm, ...
        'Floc',Floc,'eFil',eFil,'lam',lam);

    % ================= figure =================
    try
        local_figure(R, o.Savefig);
        fprintf('Saved figure: %s\n', o.Savefig);
    catch ME
        fprintf('(figure skipped: %s)\n', ME.message);
    end
end

% -----------------------------------------------------------------------------
function [v0,w0,gap] = local_pick_sulcus(X, N)
% Find an antiparallel-normal vertex pair across a narrow sulcus, gap ~ 3 mm.
    K = 30;
    try
        [idx, dst] = knnsearch(X, X, 'K', K+1);   % Stats TB
    catch
        [idx, dst] = local_knn(X, K+1);
    end
    best = inf; v0 = 1; w0 = 1;
    for v = 1:size(X,1)
        nb = idx(v,2:end); dd = dst(v,2:end);     % 1 x K rows
        nd = (N(nb,:) * N(v,:).').';              % 1 x K normal dot
        ok = (nd < -0.6) & (dd > 1.5e-3) & (dd < 6e-3);
        if any(ok)
            jj = find(ok); [g, m] = min(dd(jj));  % closest valid opposing wall
            score = abs(g - 3e-3);                % prefer ~3 mm gap
            if score < best, best = score; v0 = v; w0 = nb(jj(m)); gap = g; end
        end
    end
    if isinf(best), error('No antiparallel sulcus pair found.'); end
end

% -----------------------------------------------------------------------------
function [idx, dst] = local_knn(X, K)
    n = size(X,1); idx = zeros(n,K); dst = zeros(n,K);
    for i = 1:n
        dd = sqrt(sum((X - X(i,:)).^2,2));
        [s, p] = sort(dd, 'ascend');
        idx(i,:) = p(1:K).'; dst(i,:) = s(1:K).';
    end
end

% -----------------------------------------------------------------------------
function N = local_vertex_normals(X, F)
    v1 = X(F(:,1),:); v2 = X(F(:,2),:); v3 = X(F(:,3),:);
    fn = cross(v2-v1, v3-v1, 2);               % area-weighted face normals
    N = zeros(size(X));
    for k = 1:3, N(:,k) = accumarray(F(:), repmat(fn(:,k),3,1), [size(X,1) 1]); end
    N = local_unit(N);
end

% -----------------------------------------------------------------------------
function [L, M] = local_cotmatrix_mass(X, F)
% Cotangent Laplacian (PSD) + lumped (barycentric) mass for a triangle mesh.
    n = size(X,1);
    i1=F(:,1); i2=F(:,2); i3=F(:,3);
    e1=X(i3,:)-X(i2,:); e2=X(i1,:)-X(i3,:); e3=X(i2,:)-X(i1,:);
    fa = 0.5*sqrt(sum(cross(e1,-e3,2).^2,2));   % face area
    cot1 = -sum(e2.*e3,2)./(2*fa);              % cot at vertex i1 (opp edge e1)
    cot2 = -sum(e3.*e1,2)./(2*fa);
    cot3 = -sum(e1.*e2,2)./(2*fa);
    I = [i2;i3; i3;i1; i1;i2];
    J = [i3;i2; i1;i3; i2;i1];
    W = [cot1;cot1; cot2;cot2; cot3;cot3]/2;
    Lo = sparse(I,J,-W,n,n);
    L  = Lo - spdiags(sum(Lo,2),0,n,n);         % L = D - W (PSD)
    L  = (L+L')/2;
    am = accumarray(F(:), repmat(fa,3,1)/3, [n 1]);  % barycentric lumped mass
    M  = spdiags(am,0,n,n);
end

% -----------------------------------------------------------------------------
function a = local_angle(u, v)
    a = acosd( max(-1,min(1, dot(u,v)/(norm(u)*norm(v)))) );
end
function U = local_unit(M)
    nrm = sqrt(sum(M.^2,2)); nrm(nrm<eps)=1; U = M ./ nrm;
end

% -----------------------------------------------------------------------------
function local_figure(R, savefig)
    f = figure('Color','w','Position',[100 100 1500 650],'Visible','off');
    sub = {'DIRAC heat flow (frame transport)','NAIVE scalar heat flow (gauge-fixed)'};
    Vs  = {R.V3, R.Vn};
    for s = 1:2
        ax = subplot(1,2,s); hold(ax,'on');
        patch('Faces',R.Floc,'Vertices',R.X,'FaceColor',[.85 .85 .85], ...
              'EdgeColor','none','FaceAlpha',0.6,'FaceLighting','gouraud');
        Vv = Vs{s}; mg = sqrt(sum(Vv.^2,2)); thr = 0.12*max(mg);
        sel = find(mg > thr);
        sc = 0.004 / max(mg);                  % arrow scale
        U = Vv(sel,:);
        quiver3(R.X(sel,1),R.X(sel,2),R.X(sel,3), U(:,1)*sc,U(:,2)*sc,U(:,3)*sc, ...
                0,'Color',[0 0 0],'LineWidth',0.8,'MaxHeadSize',0.4);
        plot3(R.X(R.v0,1),R.X(R.v0,2),R.X(R.v0,3),'go','MarkerFaceColor','g','MarkerSize',9);
        plot3(R.X(R.w0,1),R.X(R.w0,2),R.X(R.w0,3),'ro','MarkerFaceColor','r','MarkerSize',9);
        axis(ax,'equal','off');
        ctr = R.X(R.v0,:); rad = 0.018;
        xlim(ctr(1)+[-rad rad]); ylim(ctr(2)+[-rad rad]); zlim(ctr(3)+[-rad rad]);
        camlight headlight; material dull; view(local_view(R));
        title(ax, sub{s}, 'FontWeight','bold');
    end
    sgtitle(sprintf(['Dirac heat transport of an across-fold impulse  |  ' ...
        'ambient flip v0\\rightarroww0:  DIRAC %.0f\\circ  vs  NAIVE %.0f\\circ'], ...
        R.angDirac, R.angNaive), 'FontWeight','bold');
    drawnow;
    print(f, savefig, '-dpng','-r130');
    close(f);
end
function v = local_view(R)
    % look roughly along the fold so both walls are visible
    n = R.Nrm(R.v0,:); v = n + [0 0 0.001]; v = v/norm(v); v = v*1; v = [n(1) n(2) n(3)];
    v = R.Nrm(R.v0,:);
end
