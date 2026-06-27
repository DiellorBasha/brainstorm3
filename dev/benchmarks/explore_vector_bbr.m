function explore_vector_bbr(cases)
% EXPLORE_VECTOR_BBR: compare BBR cost landscapes for three sampling/cost variants, to test
% whether curvature-aware sampling and a vector (gradient-field) cost give a SHARPER objective
% than the scalar boundary-step cost - especially for low-contrast tau.
%
%   (S) scalar / Euclidean-normal : Q=(I_gm-I_wm)/mean at v +/- d*n  (current mri_bbregister)
%   (R) scalar / ribbon-axis      : same, but sample along the white->pial axis (follows the
%                                   fold; avoids the Euclidean-normal collision in sulci)  [#1]
%   (G) gradient-field / NGF      : align the PET intensity GRADIENT direction with the surface
%                                   normal, reliability-weighted by |grad|^2  [#2]
%
% For each scan we sweep a known x-shift and plot cost vs shift; a deeper/narrower dip at 0 =
% a sharper, more registrable objective. Exploration only (no DB writes).
%
% USAGE:  explore_vector_bbr()   % defaults: amyloid (high contrast) -> tau 16% -> tau 6% outlier
%
% Author: Diellor Basha, 2026

    if (nargin<1)||isempty(cases)
        cases={'sub-MTL0002','18FNAV4694','amyloid 68%'; 'sub-MTL0002','18Fflortaucipir','tau 16%'; 'sub-MTL0020','18Fflortaucipir','tau 6% (outlier)'};
    end
    here=bst_fileparts(mfilename('fullpath'));
    shifts=-12:1.5:12; nSub=3000;
    f=figure('Visible','off','Position',[50 50 1400 430]);
    for ci=1:size(cases,1)
        subj=cases{ci,1}; trc=cases{ci,2}; lab=cases{ci,3}; [sS,~]=bst_get('Subject',subj);
        fn=@(c) sS.Anatomy(find(strcmp({sS.Anatomy.Comment},c),1)).FileName;
        sf=@(rx) sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},rx,'once')),1)).FileName;
        sT1=in_mri_bst(fn('MRI T1')); sPet=in_mri_bst(fn(['PET ' trc '_mean'])); PET=double(sPet.Cube(:,:,:,1)); cs=size(PET);
        sw=in_tess_bst(sf('cortex_white_low\.mat$')); sp=in_tess_bst(sf('cortex_pial_low\.mat$'));
        Vw=cs_convert(sT1,'scs','world',sw.Vertices); Vp=cs_convert(sT1,'scs','world',sp.Vertices);
        wpm=norm(cs_convert(sT1,'voxel','world',[2 1 1])-cs_convert(sT1,'voxel','world',[1 1 1]))/sT1.Voxsize(1);
        N=local_norm(Vw,sw.Faces); U=Vp-Vw; thick=sqrt(sum(U.^2,2)); U=U./max(thick,eps);   % ribbon axis
        sel=round(linspace(1,size(Vw,1),nSub)); Vw=Vw(sel,:);Vp=Vp(sel,:);N=N(sel,:);U=U(sel,:);thick=thick(sel);
        d=1.5*wpm; dWM=1.0*wpm; dlt=1.0*wpm; dG=min(d, 0.6*thick);
        smp=@(P) local_interp(PET, cs_convert(sPet,'world','voxel',P), cs);   % sample world pts
        % contrast sign per method, at shift 0 (normals can be inward -> sign differs from ribbon)
        Q0=local_contrast(smp(Vw+d.*N), smp(Vw-d.*N)); sgnS=sign(median(Q0(isfinite(Q0)))); if sgnS==0,sgnS=1;end
        Q0R=local_contrast(smp(Vw+dG.*U), smp(Vw-dWM.*U)); sgnR=sign(median(Q0R(isfinite(Q0R)))); if sgnR==0,sgnR=1;end
        cS=zeros(size(shifts)); cR=cS; cG=cS;
        for k=1:numel(shifts)
            sh=[shifts(k)*wpm 0 0];
            % (S) scalar Euclidean normal
            cS(k)=local_scost(local_contrast(smp(Vw+d.*N+sh), smp(Vw-d.*N+sh)), sgnS);
            % (R) scalar ribbon axis (GM along white->pial, capped inside the ribbon; WM just inside white)
            cR(k)=local_scost(local_contrast(smp(Vw+dG.*U+sh), smp(Vw-dWM.*U+sh)), sgnR);
            % (G) gradient-field / NGF at the white surface
            P=Vw+sh; gx=smp(P+[dlt 0 0])-smp(P-[dlt 0 0]); gy=smp(P+[0 dlt 0])-smp(P-[0 dlt 0]); gz=smp(P+[0 0 dlt])-smp(P-[0 0 dlt]);
            g=[gx gy gz]; gn=sqrt(sum(g.^2,2)); align=(sum(g.*N,2)).^2./(gn.^2+ (0.3*median(gn(gn>0)))^2);
            w=gn.^2./(gn.^2+(0.3*median(gn(gn>0)))^2); ok=isfinite(align)&isfinite(w);
            cG(k)=1 - sum(w(ok).*align(ok))/sum(w(ok));
        end
        % normalize each landscape 0..1 for shape comparison
        nz=@(x)(x-min(x))/(max(x)-min(x)+eps);
        subplot(1,3,ci); plot(shifts,nz(cS),'-o','LineWidth',1.3); hold on; plot(shifts,nz(cR),'-s','LineWidth',1.3); plot(shifts,nz(cG),'-^','LineWidth',1.3);
        xline(0,'k:'); grid on; xlabel('x-shift (mm)'); ylabel('cost (normalized)');
        title(sprintf('%s / %s', subj, lab),'Interpreter','none'); if ci==1, legend({'scalar normal','scalar ribbon (#1)','gradient/NGF (#2)'},'Location','north'); end
    end
    sgtitle('BBR objective: scalar vs curvature-aware-sampling vs vector(gradient) cost');
    print(f,fullfile(here,'explore_vector_bbr.png'),'-dpng','-r110'); close(f);
    fprintf('Figure: %s\n', fullfile(here,'explore_vector_bbr.png'));
end

function Q=local_contrast(Ig,Iw), Q=100*(Ig-Iw)./(0.5*(Ig+Iw)+eps); end
function c=local_scost(Q,sgn), Q=Q(isfinite(Q)); if numel(Q)<50,c=1;return;end; c=mean(1-tanh(sgn*0.5*Q)); end
function I=local_interp(PET,vox,cs)
    ok=all(vox>=1,2)&vox(:,1)<=cs(1)&vox(:,2)<=cs(2)&vox(:,3)<=cs(3); I=nan(size(vox,1),1);
    if any(ok), I(ok)=interpn(PET,vox(ok,1),vox(ok,2),vox(ok,3),'linear',NaN); end
end
function N=local_norm(V,F)
    v1=V(F(:,1),:);v2=V(F(:,2),:);v3=V(F(:,3),:); fn=cross(v3-v1,v2-v1,2); n=size(V,1);N=zeros(n,3);idx=F(:);  % outward (Brainstorm winding)
    for d=1:3, N(:,d)=accumarray(idx,repmat(fn(:,d),3,1),[n 1]); end
    N=N./(sqrt(sum(N.^2,2))+eps);
end
