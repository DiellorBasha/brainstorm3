function R = prototype_face_vol2surf(SubjectName, tracer)
% PROTOTYPE_FACE_VOL2SURF: face-native PET->surface projection (DEC view) vs the current
% tess_interp_mri projection.
%
% Face-native: PET areal uptake is a 2-form on faces. Sample the SUVR volume at each MID-surface
% face centroid (trilinear point sample in the face frame, ~projfrac 0.5), giving a per-FACE value;
% transfer to a vertex 0-form by the area-weighted (Hodge/mass) face->vertex map. Regional SUVR =
% area-weighted face integral over the region.
% Current: tess_interp_mri (triangle-raster neighbourhood average) at white/mid/pial, [0.1 0.8 0.1].
%
% Key test: which surface projection best preserves the VOLUME regional SUVR (the trusted source)?
% Plus a sharpness measure (LBO Dirichlet energy) and area-weighted vs vertex-mean regional.
%
% Author: Diellor Basha, 2026

    if (nargin<1)||isempty(SubjectName), SubjectName='sub-MTL0002'; end
    if (nargin<2)||isempty(tracer), tracer='18FNAV4694'; end
    here=bst_fileparts(mfilename('fullpath'));
    [sS,~]=bst_get('Subject',SubjectName);
    an=@(c) in_mri_bst(sS.Anatomy(find(strcmp({sS.Anatomy.Comment},c),1)).FileName);
    sf=@(rx) sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},rx,'once')),1)).FileName;
    sT1=an('MRI T1'); sAseg=an('ASEG'); sDKv=an('Desikan-Killiany');
    % SUVR volume = robust-cerebellar SUVR of the PVC (PETPVE12-MG) volume
    sSuvr=pet_suvr(an(['PET ' tracer '_mean_pvc']), sAseg);
    cube=double(sSuvr.Cube); cs=size(cube);
    wf=sf('cortex_white_low\.mat$'); mf=sf('cortex_mid_low\.mat$'); pf=sf('cortex_pial_low\.mat$');
    sMid=in_tess_bst(mf); V=sMid.Vertices; F=sMid.Faces;

    % ---- FACE-NATIVE: sample at mid face centroids, area-weighted face->vertex ----
    cen = (V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;
    vox = cs_convert(sSuvr,'scs','voxel', cen);
    faceVal = local_interp(cube, vox, cs);
    Af = 0.5*sqrt(sum(cross(V(F(:,2),:)-V(F(:,1),:), V(F(:,3),:)-V(F(:,1),:),2).^2,2));   % triangle area
    nV=size(V,1); num=zeros(nV,1); den=zeros(nV,1); fv=faceVal; fv(~isfinite(fv))=0; okf=isfinite(faceVal);
    for c=1:3
        num=num+accumarray(F(:,c), (Af/3).*fv.*okf, [nV 1]);
        den=den+accumarray(F(:,c), (Af/3).*okf,      [nV 1]);
    end
    mapFace = num./max(den,eps);

    % ---- CURRENT: tess_interp_mri at white/mid/pial, [0.1 0.8 0.1] ----
    cv=cube(:); mt=cell(1,3); files={wf,mf,pf};
    for k=1:3
        t2m=tess_interp_mri(files{k}, sT1);
        iv=(t2m'*cv)./max(sum(t2m,1)',eps); iv(~isfinite(iv))=0; mt{k}=iv;
    end
    mapTess = [mt{1} mt{2} mt{3}]*[0.1;0.8;0.1];

    % ---- sharpness (LBO Dirichlet energy f'Kf / f'Mf, per hemi summed) ----
    LBO=tess_operators(mf,'Laplace-Beltrami');
    dE=@(f) local_dirichlet(f, LBO);
    % ---- regional fidelity to VOLUME ----
    sc=local_scouts(sMid); volReg=local_volreg(sDKv, cube);
    [rn, vR]=deal(volReg.name, volReg.val);
    fReg_v=zeros(numel(sc),1); tReg_v=fReg_v; fReg_a=fReg_v; vRef=fReg_v;   % surf regional (vertex-mean / area-wtd) + volume
    for i=1:numel(sc)
        vv=sc(i).Vertices; ff=all(ismember(F,vv),2);                       % faces fully inside the scout
        fReg_v(i)=mean(mapFace(vv),'omitnan'); tReg_v(i)=mean(mapTess(vv),'omitnan');
        if any(ff), fReg_a(i)=sum(Af(ff).*faceVal(ff),'omitnan')/sum(Af(ff)); else, fReg_a(i)=fReg_v(i); end
        j=find(strcmp(rn, sc(i).key),1); vRef(i)=NaN; if ~isempty(j), vRef(i)=vR(j); end
    end
    ok=isfinite(vRef)&isfinite(fReg_v)&isfinite(tReg_v);
    cc=@(a,b) subsref(corrcoef(a(ok),b(ok)),struct('type','()','subs',{{1,2}}));
    rm=@(a,b) sqrt(mean((a(ok)-b(ok)).^2));
    fprintf('\n%s / %s  (%d regions vs volume SUVR)\n', SubjectName, tracer, nnz(ok));
    fprintf('  REGIONAL FIDELITY to volume :  face-native r=%.3f RMSE=%.3f | tess_interp r=%.3f RMSE=%.3f\n', cc(fReg_v,vRef),rm(fReg_v,vRef),cc(tReg_v,vRef),rm(tReg_v,vRef));
    fprintf('  area-weighted vs vertex-mean (face): r=%.4f, mean|diff|=%.4f (the 2-form correction)\n', cc(fReg_a,fReg_v), mean(abs(fReg_a(ok)-fReg_v(ok))));
    fprintf('  SHARPNESS (Dirichlet E, higher=sharper): face-native=%.3g | tess_interp=%.3g (ratio %.2fx)\n', dE(mapFace),dE(mapTess),dE(mapFace)/dE(mapTess));
    fprintf('  corr(face-native, tess_interp) maps = %.3f\n', subsref(corrcoef(mapFace,mapTess),struct('type','()','subs',{{1,2}})));
    R=struct('mapFace',mapFace,'mapTess',mapTess,'fReg_v',fReg_v,'tReg_v',tReg_v,'fReg_a',fReg_a,'vRef',vRef);

    f=figure('Visible','off','Position',[60 60 1000 430]);
    subplot(1,2,1); plot(vRef(ok),fReg_v(ok),'r.','MarkerSize',9); hold on; plot(vRef(ok),tReg_v(ok),'b.','MarkerSize',9);
    lim=[min(vRef(ok)) max(vRef(ok))]; plot(lim,lim,'k--'); axis equal; grid on; xlabel('VOLUME regional SUVR'); ylabel('SURFACE regional SUVR');
    legend({'face-native','tess\_interp','identity'},'Location','northwest'); title('regional fidelity to volume');
    subplot(1,2,2); histogram(mapFace,40,'FaceColor','r','FaceAlpha',.5); hold on; histogram(mapTess,40,'FaceColor','b','FaceAlpha',.5);
    xlabel('vertex SUVR'); ylabel('count'); legend({'face-native','tess\_interp'}); title('vertex value distribution (tess\_interp narrower = smoothed)');
    print(f,fullfile(here,'prototype_face_vol2surf.png'),'-dpng','-r110'); close(f);
    fprintf('Figure: %s\n', fullfile(here,'prototype_face_vol2surf.png'));
end

function I=local_interp(C,vox,cs)
    ok=all(vox>=1,2)&vox(:,1)<=cs(1)&vox(:,2)<=cs(2)&vox(:,3)<=cs(3); I=nan(size(vox,1),1);
    if any(ok), I(ok)=interpn(C,vox(ok,1),vox(ok,2),vox(ok,3),'linear',NaN); end
end
function e=local_dirichlet(f, LBO)
    num=0; den=0;
    for hh=1:numel(LBO.Operator)
        gv=double(LBO.GlobalVertices{hh}(:)); fh=f(gv); fh(~isfinite(fh))=0;
        num=num+fh'*LBO.Operator{hh}*fh; den=den+fh'*LBO.Mass{hh}*fh;
    end
    e=num/max(den,eps);
end
function sc=local_scouts(sMid)
    A=sMid.Atlas; ai=find(~cellfun('isempty',regexp({A.Name},'Desikan|aparc','once')),1); s=A(ai).Scouts;
    for i=1:numel(s), s(i).key=lower(strrep(s(i).Label,' ','')); end
    sc=s;
end
function vr=local_volreg(sDKv, cube)
    L=sDKv.Labels; v=cell2mat(L(:,1)); nm=L(:,2); k=find(v>=1000 & v<3000);
    name=cell(numel(k),1); val=zeros(numel(k),1);
    for i=1:numel(k)
        name{i}=lower(strrep(nm{k(i)},' ','')); m=(sDKv.Cube==v(k(i)));
        val(i)=mean(cube(m),'omitnan');
    end
    vr=struct('name',{name},'val',val);
end
