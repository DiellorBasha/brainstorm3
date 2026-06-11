function R = bench_sensor_dirac(base)
% BENCH_SENSOR_DIRAC  EXPERIMENTAL: geometric eigenmodes of the MEG sensor helmet.
%
% Builds a Laplacian on the SENSOR NETWORK (sensors = vertices of a triangulated
% helmet; signal lives on vertices, so the discrete point set is exact), takes its
% eigenmodes ("helmet harmonics"), and runs two diagnostics:
%   (A) how many helmet harmonics span the leadfield's sensor topographies (the
%       field patterns sources actually produce) -> the geometric dimensionality
%       of the observable sensor subspace.
%   (B) do the dominant NOISE eigenmodes load onto LOW helmet harmonics (smooth /
%       global / external = SSS-removable) or high ones (local sensor noise)?
%       -> whether the noise covariance carries helmet-geometry structure.
%
% Notes (design): orientation lives in the leadfield, not the helmet operator, so
% this is the SCALAR helmet Laplacian over sensor positions. Built directly in
% MATLAB (cotan Laplacian) -- NOT via nxr_compute('create') (segfault rule).
%
% USAGE:  R = bench_sensor_dirac
% Author: Diellor Basha, 2026

    if nargin<1 || isempty(base), base='Subject01/S01_AEF_20131218_01_notch/'; end
    OUTDIR='/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks';
    fprintf('\n=== EXPERIMENTAL: sensor-helmet Laplacian eigenmodes ===\n');

    % ---- sensor positions (CTF MEG), leadfield, noise cov ----
    ChanMat=in_bst_channel([base 'channel_ctf_acc1.mat']); types={ChanMat.Channel.Type};
    HMos=in_bst_headmodel([base 'headmodel_surf_os_meg.mat'],0);
    G=double(HMos.Gain); iMEG=all(isfinite(G),2)&strcmpi(types(:),'MEG'); G=G(iMEG,:);
    idx=find(iMEG); n=numel(idx);
    pos=zeros(n,3); ornt=zeros(n,3);
    for i=1:n
        L=ChanMat.Channel(idx(i)).Loc;   pos(i,:)=mean(L,2)';
        O=ChanMat.Channel(idx(i)).Orient; if ~isempty(O), ornt(i,:)=mean(O,2)'; end
    end
    NC=load(file_fullpath([base 'noisecov_full.mat'])); Cn=NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
    fprintf('Sensors: %d MEG channels on the helmet\n', n);

    % ---- triangulate the helmet: sphere-fit -> azimuthal flatten -> open Delaunay ----
    c=sphere_fit(pos); u=(pos-c); u=u./sqrt(sum(u.^2,2));
    apex=mean(u,1); apex=apex/norm(apex); Rot=rot_to_z(apex); ur=(Rot*u')';
    th=acos(min(1,max(-1,ur(:,3)))); az=atan2(ur(:,2),ur(:,1));
    x2=th.*cos(az); y2=th.*sin(az);                 % 2D helmet layout (topomap coords)
    T=delaunay(x2,y2);
    % drop spurious long-edge triangles near the rim
    el=@(a,b)sqrt(sum((pos(T(:,a),:)-pos(T(:,b),:)).^2,2));
    emax=max([el(1,2) el(2,3) el(3,1)],[],2); T=T(emax < 3*median(emax),:);
    fprintf('Helmet mesh: %d triangles; sensor radiality (|orient . r|) median=%.2f\n', ...
        size(T,1), median(abs(sum(ornt.*u,2))));

    % ---- cotan Laplacian + lumped mass; helmet harmonics (M-orthonormal) ----
    [Lap,Mass]=cotlap(pos,T);
    Mi=spdiags(1./sqrt(max(full(diag(Mass)),eps)),0,n,n);
    Ls=Mi*Lap*Mi; Ls=(Ls+Ls')/2; [Vs,D]=eig(full(Ls));
    [lh,ix]=sort(real(diag(D)),'ascend'); Phi=Mi*Vs(:,ix);   % Phi'*Mass*Phi = I
    R.lambda=lh; R.Phi=Phi; R.pos=pos; R.xy=[x2 y2]; R.T=T;
    fprintf('Helmet Laplacian: %d eigenmodes; lambda in [%.2e, %.2e]\n', n, lh(2), lh(end));

    proj=@(P) Phi'*(Mass*P);     % coefficients of sensor pattern(s) P in helmet harmonics

    % ---- (A) leadfield sensor topographies in helmet harmonics ----
    Cg=proj(G); eg=sum(Cg.^2,2); eg=eg/sum(eg); cumg=cumsum(eg);
    kA=[find(cumg>=0.90,1) find(cumg>=0.99,1)];
    fprintf('(A) leadfield topographies captured by helmet harmonics: 90%%=%d, 99%%=%d of %d\n', kA(1), kA(2), n);

    % ---- (B) noise eigenmodes vs helmet geometry ----
    [Un,Dn]=eig(Cn); dn=real(diag(Dn)); [dn,jx]=sort(dn,'descend'); Un=Un(:,jx);
    Cn_h=proj(Un);                                   % [n x n] noise modes in helmet harmonics
    hfreq=sum((lh).*(Cn_h.^2),1)' ./ sum(Cn_h.^2,1)';% helmet "spatial frequency" of each noise mode
    % variance-weighted mean helmet-frequency of the noise (dominant modes weighted by dn)
    R.noiseHfreqTop = sum(dn(1:20).*hfreq(1:20))/sum(dn(1:20));
    R.noiseHfreqAll = sum(dn.*hfreq)/sum(dn);
    % fraction of noise variance whose modes sit in the lower-half of the helmet spectrum
    medL=median(lh(2:end)); R.noiseLowFrac = sum(dn(hfreq<medL))/sum(dn);
    fprintf('(B) noise modes'' helmet-frequency: top-20 weighted=%.2e, all weighted=%.2e (median lambda=%.2e)\n', ...
        R.noiseHfreqTop, R.noiseHfreqAll, medL);
    fprintf('    fraction of noise variance in LOW helmet-frequency modes = %.0f%% (high => spatially smooth/global noise)\n', 100*R.noiseLowFrac);

    % ---- figure ----
    f=figure('Color','w','Position',[40 40 1500 860],'Visible','off');
    subplot(2,3,1); plot(lh,'k.','MarkerSize',8); xlabel('mode index'); ylabel('\lambda (helmet)');
    title('Sensor-helmet Laplacian spectrum'); grid on;
    for m=1:3      % modes 2,3,4 (skip the DC constant mode 1)
        subplot(2,3,1+m); scatter(x2,y2,55,Phi(:,m+1),'filled'); axis equal off;
        colormap(gca,parula); title(sprintf('helmet harmonic %d (\\lambda=%.2g)',m+1,lh(m+1)));
    end
    subplot(2,3,5); plot(cumg*100,'LineWidth',2,'Color',[.1 .5 .2]); hold on;
    plot([kA(1) kA(1)],[0 90],'r--'); xlabel('# helmet harmonics'); ylabel('% leadfield topo energy');
    title(sprintf('(A) leadfield topographies: 90%%@%d',kA(1))); grid on; ylim([0 100]); xlim([1 n]);
    subplot(2,3,6); semilogx(hfreq,dn,'.','Color',[.3 .3 .6]); hold on; xline(medL,'k:');
    xlabel('helmet spatial frequency of noise mode'); ylabel('noise eigenvalue (variance)');
    title(sprintf('(B) noise vs geometry (%.0f%% low-freq)',100*R.noiseLowFrac)); grid on; set(gca,'YScale','log');
    sgtitle('EXPERIMENTAL: MEG sensor-helmet geometric eigenmodes + diagnostics','FontWeight','bold');
    print(f,[OUTDIR '/bench_sensor_dirac.png'],'-dpng','-r110'); close(f);
    fprintf('Saved %s/bench_sensor_dirac.png\n', OUTDIR);
end

% ===== helpers =====
function c = sphere_fit(P)
    A=[2*P ones(size(P,1),1)]; b=sum(P.^2,2); x=A\b; c=x(1:3)';
end
function Rm = rot_to_z(a)
    a=a(:)'/norm(a); z=[0 0 1]; v=cross(a,z); s=norm(v); cc=dot(a,z);
    if s<1e-9, Rm=eye(3)*sign(cc); if cc<0, Rm=diag([1 -1 -1]); end; return; end
    Vx=[0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0]; Rm=eye(3)+Vx+Vx^2*((1-cc)/s^2);
end
function [L,M] = cotlap(V,F)
    n=size(V,1); i1=F(:,1);i2=F(:,2);i3=F(:,3);
    e1=V(i3,:)-V(i2,:); e2=V(i1,:)-V(i3,:); e3=V(i2,:)-V(i1,:);
    fa=0.5*sqrt(sum(cross(e1,-e3,2).^2,2));
    c1=-sum(e2.*e3,2)./(2*fa); c2=-sum(e3.*e1,2)./(2*fa); c3=-sum(e1.*e2,2)./(2*fa);
    I=[i2;i3;i3;i1;i1;i2]; J=[i3;i2;i1;i3;i2;i1]; W=[c1;c1;c2;c2;c3;c3]/2;
    Lo=sparse(I,J,-W,n,n); L=Lo-spdiags(sum(Lo,2),0,n,n); L=(L+L')/2;
    am=accumarray(F(:),repmat(fa,3,1)/3,[n 1]); M=spdiags(am,0,n,n);
end
