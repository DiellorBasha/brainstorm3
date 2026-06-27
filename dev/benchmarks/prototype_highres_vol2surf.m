function R = prototype_highres_vol2surf(SubjectName, tracer)
% PROTOTYPE_HIGHRES_VOL2SURF: PET->surface projection on the HIGH-RES surface with medial-wall
% masking - the two real levers identified for the projection (resolution + cortex mask), vs the
% current low-res tess_interp_mri map.
%
% Method: sample the SUVR volume at each high-res MID vertex (trilinear point sample, projfrac 0.5;
% face-vs-vertex shown neutral). Mask the medial wall (thickness < MinThick). Compare detail and
% regional fidelity to the volume against the low-res tess_interp map, and render both.
%
% Author: Diellor Basha, 2026

    if (nargin<1)||isempty(SubjectName), SubjectName='sub-MTL0002'; end
    if (nargin<2)||isempty(tracer), tracer='18FNAV4694'; end
    MinThick=1.0; here=bst_fileparts(mfilename('fullpath'));
    [sS,~]=bst_get('Subject',SubjectName);
    an=@(c) in_mri_bst(sS.Anatomy(find(strcmp({sS.Anatomy.Comment},c),1)).FileName);
    sf=@(rx) sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},rx,'once')),1)).FileName;
    sT1=an('MRI T1'); sAseg=an('ASEG'); sDKv=an('Desikan-Killiany');
    sSuvr=pet_suvr(an(['PET ' tracer '_mean_pvc']), sAseg); cube=double(sSuvr.Cube); cs=size(cube);

    % ---- HIGH-RES projection: mid = 0.5(white,pial), trilinear at mid vertices, medial-wall mask ----
    sWh=in_tess_bst(sf('cortex_white_high\.mat$')); sPh=in_tess_bst(sf('cortex_pial_high\.mat$'));
    Vw=sWh.Vertices; Vp=sPh.Vertices; Fh=sWh.Faces;
    th=sqrt(sum((Vp-Vw).^2,2))*1000;                          % thickness mm
    Vmid=0.5*(Vw+Vp);
    vox=cs_convert(sSuvr,'scs','voxel',Vmid); mapHi=local_interp(cube,vox,cs);
    maskHi = th < MinThick;                                    % medial wall / non-cortex
    mapHi(maskHi)=NaN;

    % ---- LOW-RES current method: tess_interp_mri white/mid/pial, [0.1 0.8 0.1], no mask ----
    wf=sf('cortex_white_low\.mat$'); mf=sf('cortex_mid_low\.mat$'); pf=sf('cortex_pial_low\.mat$');
    sPl=in_tess_bst(pf); cv=cube(:); mt=cell(1,3); files={wf,mf,pf};
    for k=1:3, t2m=tess_interp_mri(files{k},sT1); iv=(t2m'*cv)./max(sum(t2m,1)',eps); iv(~isfinite(iv))=0; mt{k}=iv; end
    mapLo=[mt{1} mt{2} mt{3}]*[0.1;0.8;0.1];

    % ---- regional fidelity to VOLUME (high-res vs low-res) ----
    volReg=local_volreg(sDKv,cube);
    rHi=local_regfid(sPh,mapHi,volReg); rLo=local_regfid(sPl,mapLo,volReg);
    fprintf('\n%s / %s\n', SubjectName, tracer);
    fprintf('  HIGH-RES: %d vert, medial-wall masked %d (%.1f%%) | regional-fidelity-to-volume r=%.3f RMSE=%.3f\n', numel(mapHi),nnz(maskHi),100*mean(maskHi),rHi.r,rHi.rmse);
    fprintf('  LOW-RES : %d vert, no mask                       | regional-fidelity-to-volume r=%.3f RMSE=%.3f\n', numel(mapLo),rLo.r,rLo.rmse);
    fprintf('  vertex spacing: high-res ~1mm vs low-res ~3-4mm (13x more vertices)\n');
    R=struct('mapHi',mapHi,'maskHi',maskHi,'mapLo',mapLo);

    % ---- geodesic surface smoothing (~3mm, NaN-aware neighbour averaging) for the usable map ----
    mapHiSm = local_smooth(mapHi, Fh, numel(mapHi), 8);
    R.mapHiSm = mapHiSm;

    % ---- render: low-res current | high-res raw (noisy) | high-res smoothed (usable) ----
    rng=[0.9 1.9]; vw=[-110 25];
    f=figure('Visible','off','Position',[40 40 1500 480]);
    local_render(subplot(1,3,1), sPl.Vertices, sPl.Faces, mapLo, rng, vw, 'LOW-RES tess\_interp (current)');
    local_render(subplot(1,3,2), Vp, Fh, mapHi,   rng, vw, 'HIGH-RES raw (voxel noise)');
    local_render(subplot(1,3,3), Vp, Fh, mapHiSm, rng, vw, 'HIGH-RES masked + smoothed (usable)');
    colormap(hot);
    print(f,fullfile(here,'prototype_highres_vol2surf.png'),'-dpng','-r110'); close(f);
    fprintf('Figure: %s\n', fullfile(here,'prototype_highres_vol2surf.png'));
end

function local_render(ax, V, F, map, rng, vw, ttl)
    axes(ax); m=map; patch('Faces',F,'Vertices',V,'FaceVertexCData',m,'FaceColor','interp','EdgeColor','none','FaceAlpha',1); %#ok<MAXES>
    caxis(rng); view(vw); axis(ax,'equal','off','tight'); camlight headlight; lighting gouraud; material dull; title(ttl); set(ax,'Color','k');
end
function g=local_smooth(map, F, nV, nIter)
    % NaN-aware neighbour averaging (cheap geodesic smoothing); masked (NaN) vertices stay NaN.
    A=sparse([F(:,1);F(:,2);F(:,3);F(:,2);F(:,3);F(:,1)],[F(:,2);F(:,3);F(:,1);F(:,1);F(:,2);F(:,3)],1,nV,nV);
    A=double(A>0); g=map; val=isfinite(map); g(~val)=0;
    for it=1:nIter
        s=A*g; c=A*double(val); gn=(g+s)./(1+c); g(val)=gn(val);
    end
    g(~isfinite(map))=NaN;
end
function I=local_interp(C,vox,cs)
    ok=all(vox>=1,2)&vox(:,1)<=cs(1)&vox(:,2)<=cs(2)&vox(:,3)<=cs(3); I=nan(size(vox,1),1);
    if any(ok), I(ok)=interpn(C,vox(ok,1),vox(ok,2),vox(ok,3),'linear',NaN); end
end
function out=local_regfid(sSurf,map,volReg)
    A=sSurf.Atlas; ai=find(~cellfun('isempty',regexp({A.Name},'Desikan|aparc','once')),1); s=A(ai).Scouts;
    sv=zeros(numel(s),1); vv=sv;
    for i=1:numel(s)
        key=lower(strrep(s(i).Label,' ','')); sv(i)=mean(map(s(i).Vertices),'omitnan');
        j=find(strcmp(volReg.name,key),1); vv(i)=NaN; if ~isempty(j), vv(i)=volReg.val(j); end
    end
    ok=isfinite(sv)&isfinite(vv); cc=corrcoef(sv(ok),vv(ok));
    out=struct('r',cc(1,2),'rmse',sqrt(mean((sv(ok)-vv(ok)).^2)));
end
function vr=local_volreg(sDKv,cube)
    L=sDKv.Labels; v=cell2mat(L(:,1)); nm=L(:,2); k=find(v>=1000&v<3000); name=cell(numel(k),1); val=zeros(numel(k),1);
    for i=1:numel(k), name{i}=lower(strrep(nm{k(i)},' ','')); val(i)=mean(cube(sDKv.Cube==v(k(i))),'omitnan'); end
    vr=struct('name',{name},'val',val);
end
