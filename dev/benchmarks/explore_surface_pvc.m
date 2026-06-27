function R = explore_surface_pvc(SubjectName, psfFwhmMm)
% EXPLORE_SURFACE_PVC: surface-native PVC in the intrinsic/extrinsic decomposition, and a test
% of where the radial/tangential DECOUPLING breaks (high-curvature sulcal fundi / cross-bank spill).
%
% Framework (this session): PVE = a RADIAL/EXTRINSIC, signed tissue spill (WM->GM->CSF along the
% normal) + a TANGENTIAL/INTRINSIC blur of real GM signal (no tissue spill - it's all cortex).
%   - radial spill  -> surface-native Mueller-Gaertner: GM = (obs - m_wm*WM)/m_gm, with the PSF
%                      tissue fractions m_wm,m_gm computed per-vertex from the cortical THICKNESS.
%   - tangential blur-> LBO Wiener deconvolution (intrinsic), exp(-t.lam) inverse, regularized.
% The decoupling is FIRST ORDER: at fundi the PSF mixes in the OPPOSITE BANK (visible only
% extrinsically), which the per-vertex radial column and the intrinsic LBO both miss.
%
% Forward model: a sparse 3-D (Euclidean) PSF kernel over mid-surface vertices (captures
% tangential AND cross-bank mixing) + a radial WM spill-in term. We then run radial MG, then
% tangential Wiener, and score recovery STRATIFIED by mean curvature |H| (from Dx = 2H.n, the
% LBO of the embedding) and by the directly-measured cross-bank fraction.
%
% Exploration only (one hemisphere, ico5 - coarse vs the true PSF, so a wide model PSF is used
% to make the mixing visible; a finer surface is the follow-up).
%
% USAGE:  R = explore_surface_pvc('sub-MTL0002', 6)
%
% Author: Diellor Basha, 2026

    if (nargin<1)||isempty(SubjectName), SubjectName='sub-MTL0002'; end
    if (nargin<2)||isempty(psfFwhmMm), psfFwhmMm=6; end
    here=bst_fileparts(mfilename('fullpath')); nModes=600;

    % ----- geometry: white, pial -> mid, normal (white->pial), thickness; LBO + eigenbasis -----
    [sS,~]=bst_get('Subject',SubjectName);
    sf=@(rx) sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},rx,'once')),1)).FileName;
    wf=sf('cortex_white_low\.mat$'); sw=in_tess_bst(wf); sp=in_tess_bst(sf('cortex_pial_low\.mat$'));
    LBO=tess_operators(wf,'Laplace-Beltrami'); Eig=tess_eigen(wf,'Laplace-Beltrami','nModes',nModes);
    vH=double(LBO.GlobalVertices{1}(:)); K=LBO.Operator{1}; M=LBO.Mass{1};
    Phi=Eig.Phi{1}; lam=Eig.Lambda{1}(:);
    Vw=sw.Vertices(vH,:); Vp=sp.Vertices(vH,:); Xmid=0.5*(Vw+Vp);
    Uvec=Vp-Vw; th=sqrt(sum(Uvec.^2,2)); nrm=Uvec./max(th,eps);     % GM-ward normal + thickness
    area=full(sum(M,2));                                            % vertex (Voronoi) area
    nV=numel(vH);

    % ----- mean curvature H from the LBO of the embedding: Dx = 2H.n -----
    Hvec = M\(K*Xmid); H = 0.5*sum(Hvec.*nrm,2);                    % signed mean curvature [1/m]

    % ----- synthetic GM_true (tangential structure) + WM level (tau-like: WM high) -----
    GMtrue = 1.0 + 0.4*(Xmid(:,2)-min(Xmid(:,2)))/(range(Xmid(:,2)));
    seeds=round(linspace(1,nV,14)); seeds=seeds(2:end-1);
    for s=seeds, GMtrue=GMtrue + 0.8*exp(-sum((Xmid-Xmid(s,:)).^2,2)/(2*(0.006^2))); end
    GMtrue = Phi*(Phi'*(M*GMtrue));                                 % band-limit to the basis
    WMlevel = 1.3;                                                   % WM > GM (tau non-specific)

    % ----- radial PSF tissue fractions from thickness (1-D Gaussian along the normal) -----
    sig = (psfFwhmMm/2.355)/1000;                                   % m
    m_wm = 0.5*erfc((th/2)/(sig*sqrt(2)));                          % PSF mass in WM half-space
    m_gm = max(1 - 2*m_wm, 1e-3);                                   % GM fraction (CSF symmetric)

    % ----- sparse 3-D Euclidean PSF kernel over mid vertices (tangential + cross-bank) -----
    [Wi,Wj,Wd]=local_rangepairs(Xmid, 3*sig);
    w = exp(-Wd.^2/(2*sig^2)) .* area(Wj) .* th(Wj);               % volume-weighted PSF
    W = sparse(Wi, Wj, w, nV, nV);
    rs = full(sum(W,2)); W = spdiags(m_gm./max(rs,eps),0,nV,nV)*W; % normalize GM mass to m_gm
    % cross-bank fraction: mass from vertices whose normal OPPOSES v (opposite sulcal bank)
    opp = sum(nrm(Wi,:).*nrm(Wj,:),2) < -0.1;
    cbnum = accumarray(Wi(opp), w(opp), [nV 1]); cb = cbnum./max(accumarray(Wi,w,[nV 1]),eps);

    % ----- forward: observed = tangential/cross-bank GM mix + WM spill-in -----
    observed = W*GMtrue + m_wm*WMlevel;

    % ----- (1) radial MG  -> (2) tangential LBO Wiener -----
    GMrad = (observed - m_wm*WMlevel)./m_gm;
    tTan = sig^2/2; alpha=1e-2; G=exp(-tTan*lam); h=G./(G.^2+alpha);
    GMfinal = Phi*(h.*(Phi'*(M*GMrad)));

    % ----- evaluate with ERROR-based metrics (correlation-within-stratum is variance-confounded) -----
    stg={'observed',observed; 'radial MG',GMrad; 'radial+tangential',GMfinal};
    rmse=@(e) sqrt(mean((e-GMtrue).^2));
    R=struct(); R.H=H; R.cb=cb; R.GMtrue=GMtrue; R.m_gm=m_gm; R.GMfinal=GMfinal; R.GMrad=GMrad; R.observed=observed;
    fprintf('\n%-20s  RMSE(all)\n','stage');
    for k=1:3, fprintf('%-20s  %6.4f\n', stg{k,1}, rmse(stg{k,2})); end
    % what PREDICTS the residual error? (cross-bank vs thin-cortex instability vs curvature)
    err=abs(GMfinal-GMtrue);
    fprintf('\nerror predictors (corr with |GM_final-GM_true|):\n');
    fprintf('  cross-bank fraction : %+.3f\n', local_corr(err,cb));
    fprintf('  thin cortex 1/m_gm  : %+.3f\n', local_corr(err,1./m_gm));
    fprintf('  curvature |H|       : %+.3f\n', local_corr(err,abs(H)));
    % RMSE stratified by the dominant predictor (thin-cortex 1/m_gm) and by cross-bank
    q=@(x)[-inf, quantile(x,[1/3 2/3]), inf]; eG=q(1./m_gm); ecb=q(cb);
    sb=@(e,edg,v) arrayfun(@(b) sqrt(mean((e(v>=edg(b)&v<edg(b+1))-GMtrue(v>=edg(b)&v<edg(b+1))).^2)),1:3);
    fprintf('\n%-20s  RMSE by 1/m_gm(thin) lo/mid/hi   RMSE by cross-bank lo/mid/hi\n','stage');
    for k=1:3, e=stg{k,2}; fprintf('%-20s   %6.4f %6.4f %6.4f        %6.4f %6.4f %6.4f\n', stg{k,1}, sb(e,eG,1./m_gm), sb(e,ecb,cb)); end

    % ----- figure -----
    f=figure('Visible','off','Position',[50 50 1300 420]);
    subplot(1,3,1); scatter(1./m_gm, err, 4,'filled'); xlabel('1/m_{gm}  (thin cortex ->)'); ylabel('|GM_{final}-GM_{true}|');
    title('residual error vs radial instability'); grid on;
    subplot(1,3,2); scatter(cb, err, 4,'filled'); xlabel('cross-bank fraction'); ylabel('|GM_{final}-GM_{true}|');
    title('residual error vs cross-bank'); grid on;
    subplot(1,3,3);
    rc=zeros(3,3); for k=1:3, rc(k,:)=sb(stg{k,2},eG,1./m_gm); end
    bar(rc'); set(gca,'XTickLabel',{'thick','mid','thin'}); ylabel('RMSE'); legend(stg(:,1),'Location','northwest');
    title('RMSE by cortical-thickness tercile'); grid on;
    print(f,fullfile(here,'explore_surface_pvc.png'),'-dpng','-r110'); close(f);
    fprintf('Figure: %s\n', fullfile(here,'explore_surface_pvc.png'));
end

function [I,J,D]=local_rangepairs(X, r)
    % all pairs (i,j) with |X(i)-X(j)|<r (blocked brute force; no toolbox dependency)
    n=size(X,1); I=[];J=[];D=[]; sx=sum(X.^2,2); blk=600;
    for b=1:blk:n
        idx=(b:min(b+blk-1,n))';
        D2 = sx(idx) - 2*X(idx,:)*X' + sx';
        [ii,jj]=find(D2 < r^2);
        I=[I; idx(ii)]; J=[J; jj]; D=[D; sqrt(max(D2(sub2ind(size(D2),ii,jj)),0))]; %#ok<AGROW>
    end
end

function c=local_corr(x,y), if numel(x)<5,c=NaN;return;end; cc=corrcoef(x,y); c=cc(1,2); end
