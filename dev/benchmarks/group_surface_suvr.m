function R = group_surface_suvr(subjects, tracer, fwhmMm)
% GROUP_SURFACE_SUVR: vertex-wise group-average cortical SUVR map.
%
% Per subject: SUVR volume (pet_suvr on the PETPVE12-MG PVC) -> sampled onto the subject mid
% surface (trilinear, projfrac 0.5) -> resampled to the DEFAULT-ANATOMY template cortex via the
% FreeSurfer registration spheres (tess_interp_tess2tess, which aligns cortical folding). The
% template-space maps are averaged vertex-wise; the template medial wall is masked; the group
% mean is geodesically smoothed and rendered.
%
% Author: Diellor Basha, 2026

    if (nargin<1)||isempty(subjects)
        ps=bst_get('ProtocolSubjects'); nm={ps.Subject.Name};
        subjects=nm(~cellfun('isempty',regexp(nm,'^sub-MTL\d+$','once'))); subjects=subjects(1:min(10,end));
    end
    if (nargin<2)||isempty(tracer), tracer='18FNAV4694'; end
    if (nargin<3)||isempty(fwhmMm), fwhmMm=6; end
    here=bst_fileparts(mfilename('fullpath'));

    % ---- template (default anatomy) cortex: white (for spheres+thickness), pial (display) ----
    sDef=bst_get('Subject',0);
    df=@(rx) sDef.Surface(find(~cellfun('isempty',regexp({sDef.Surface.FileName},rx,'once')),1)).FileName;
    tWhiteF=df('cortex_white_low\.mat$'); tPialF=df('cortex_pial_low\.mat$');
    if isempty(tPialF), tPialF=df('cortex_15002V'); end
    tW=in_tess_bst(tWhiteF); tP=in_tess_bst(tPialF); nT=size(tW.Vertices,1);
    accum=zeros(nT,1); cnt=zeros(nT,1); perSubjMean=[];

    for si=1:numel(subjects)
        subj=subjects{si}; [sS,e]=bst_get('Subject',subj); if isempty(sS), continue; end
        an=@(c) local_load(sS,c); sf=@(rx) local_surf(sS,rx);
        sMG=an(['PET ' tracer '_mean_pvc']); sAseg=an('ASEG'); if isempty(sMG)||isempty(sAseg), continue; end
        try
            sSuvr=pet_suvr(sMG,sAseg); cube=double(sSuvr.Cube); cs=size(cube);
            sWl=in_tess_bst(sf('cortex_white_low\.mat$')); sPl=in_tess_bst(sf('cortex_pial_low\.mat$'));
            Vmid=0.5*(sWl.Vertices+sPl.Vertices);
            vox=cs_convert(sSuvr,'scs','voxel',Vmid); subjMap=local_interp(cube,vox,cs);   % projfrac 0.5
            % resample subject -> template via registration spheres (white surfaces carry Reg.Sphere)
            Wmat=tess_interp_tess2tess(sf('cortex_white_low\.mat$'), tWhiteF, 0, 0, 0);
            v=subjMap; ok=isfinite(v); v(~ok)=0;
            tmpl=Wmat*v; w=Wmat*double(ok); tmpl=tmpl./max(w,eps); tmpl(w<0.5)=NaN;            % omitnan resample
            accum=accum+nan2zero(tmpl); cnt=cnt+isfinite(tmpl);
            perSubjMean(end+1)=mean(tmpl,'omitnan'); %#ok<AGROW>
            fprintf('  %-14s projected (subj mean cortical SUVR=%.3f)\n', subj, perSubjMean(end));
        catch ME, fprintf('  %-14s ERR %s\n', subj, ME.message); end
    end
    groupMean=accum./max(cnt,1); groupMean(cnt<1)=NaN;

    % template medial-wall mask + geodesic smoothing
    th=sqrt(sum((tP.Vertices-tW.Vertices).^2,2))*1000; groupMean(th<1)=NaN;
    groupSm=local_smooth(groupMean, tW.Faces, nT, round(fwhmMm));
    R=struct('groupMean',groupMean,'groupSm',groupSm,'n',sum(cnt>0)>0,'nSubj',numel(perSubjMean),'perSubjMean',perSubjMean,'templatePial',tPialF);

    fprintf('\nGROUP (%s): %d subjects averaged on template (%d vert)\n', tracer, numel(perSubjMean), nT);
    fprintf('  group mean cortical SUVR = %.3f (across-subject SUVR range %.2f-%.2f)\n', mean(groupMean,'omitnan'), min(perSubjMean), max(perSubjMean));

    rng=[0.9 1.7]; vw=[-110 25];
    f=figure('Visible','off','Position',[40 40 1100 460]);
    local_render(subplot(1,2,1), tP.Vertices, tP.Faces, groupMean, rng, vw, sprintf('%s group mean (n=%d), raw',tracer,numel(perSubjMean)));
    local_render(subplot(1,2,2), tP.Vertices, tP.Faces, groupSm,   rng, vw, sprintf('group mean, %dmm smoothed',fwhmMm));
    colormap(hot);
    print(f,fullfile(here,'group_surface_suvr.png'),'-dpng','-r110'); close(f);
    fprintf('Figure: %s\n', fullfile(here,'group_surface_suvr.png'));
end

function z=nan2zero(x), z=x; z(~isfinite(z))=0; end
function sM=local_load(sS,c), i=find(strcmp({sS.Anatomy.Comment},c),1); if isempty(i),sM=[];else,sM=in_mri_bst(sS.Anatomy(i).FileName);end; end
function f=local_surf(sS,rx), i=find(~cellfun('isempty',regexp({sS.Surface.FileName},rx,'once')),1); f=sS.Surface(i).FileName; end
function I=local_interp(C,vox,cs)
    ok=all(vox>=1,2)&vox(:,1)<=cs(1)&vox(:,2)<=cs(2)&vox(:,3)<=cs(3); I=nan(size(vox,1),1);
    if any(ok), I(ok)=interpn(C,vox(ok,1),vox(ok,2),vox(ok,3),'linear',NaN); end
end
function g=local_smooth(map,F,nV,nIter)
    A=double(sparse([F(:,1);F(:,2);F(:,3);F(:,2);F(:,3);F(:,1)],[F(:,2);F(:,3);F(:,1);F(:,1);F(:,2);F(:,3)],1,nV,nV)>0);
    g=map; val=isfinite(map); g(~val)=0;
    for it=1:nIter, s=A*g; c=A*double(val); gn=(g+s)./(1+c); g(val)=gn(val); end
    g(~isfinite(map))=NaN;
end
function local_render(ax,V,F,map,rng,vw,ttl)
    axes(ax); patch('Faces',F,'Vertices',V,'FaceVertexCData',map,'FaceColor','interp','EdgeColor','none'); %#ok<MAXES>
    caxis(rng); view(vw); axis(ax,'equal','off','tight'); camlight headlight; lighting gouraud; material dull; set(ax,'Color','k'); title(ttl,'Color','k');
end
