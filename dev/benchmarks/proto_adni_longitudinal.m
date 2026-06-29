function R = proto_adni_longitudinal(subjDir, refSubj)
% PROTO_ADNI_LONGITUDINAL: prototype longitudinal amyloid analysis on REAL ADNI data, using a
% stand-in anatomy (TESTING ONLY - real per-subject FreeSurfer surfaces to be integrated next).
%
% Pipeline: convert-then-NIfTI ADNI AV45 timepoints (tp_*.nii.gz in subjDir, already within-subject
% coregistered by ADNI) -> MNI-normalize the mean (SPM maff8, once) -> map the refSubj mid-cortex +
% cerebellum to MNI -> sample each timepoint -> SUVR (cerebellum ref) -> [nVert x nT] surface series
% -> pet_epicenter per timepoint (epicenter evolution). Renders the SUVR trajectory + accumulation.
%
% USAGE: R = proto_adni_longitudinal('.../adni/002_S_0413', 'sub-MTL0002')
%
% Author: Diellor Basha, 2026 (prototype; stand-in anatomy)
    if nargin<2, refSubj='sub-MTL0002'; end
    here=bst_fileparts(mfilename('fullpath'));
    d=dir(fullfile(subjDir,'tp_*.nii.gz')); nT=numel(d);
    dates=cellfun(@(n) n(4:13), {d.name}, 'uni',0);                 % tp_YYYY-MM-DD
    yrs=cellfun(@(s) datenum(s,'yyyy-mm-dd'), dates); yrs=(yrs-yrs(1))/365.25;   % years from baseline

    % --- load timepoints; MNI-normalize the mean once (shared grid => one mapping) ---
    C=[]; for i=1:nT, V=double(niftiread(fullfile(subjDir,d(i).name))); if isempty(C), C=zeros([size(V) nT]); end; C(:,:,:,i)=V; end
    sAdni=in_mri(fullfile(subjDir,d(1).name),'Nifti1',0,0); sAdni.Cube=mean(C,4);
    sAdni=bst_normalize_mni(sAdni,'maff8');

    % --- ref-subject mid-cortex + cerebellum -> MNI -> ADNI voxel ---
    [sS,~]=bst_get('Subject',refSubj); cm={sS.Anatomy.Comment}; af=@(c) sS.Anatomy(find(strcmp(cm,c),1)).FileName;
    sf=@(rx) sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},rx,'once')),1)).FileName;
    sT1=in_mri_bst(af('MRI T1')); sAseg=in_mri_bst(af('ASEG'));
    wf=sf('cortex_white_low\.mat$'); sWm=in_tess_bst(wf); sPm=in_tess_bst(sf('cortex_pial_low\.mat$'));
    Vmid=0.5*(sWm.Vertices+sPm.Vertices);
    voxCtx=cs_convert(sAdni,'mni','voxel', cs_convert(sT1,'scs','mni',Vmid));
    % cerebellum-cortex reference (ASEG 8/47), subsampled
    cb=find(ismember(sAseg.Cube,[8 47])); cb=cb(1:8:end); [ci,cj,ck]=ind2sub(size(sAseg.Cube),cb);
    voxCb=cs_convert(sAdni,'mni','voxel', cs_convert(sT1,'scs','mni', cs_convert(sAseg,'voxel','scs',[ci cj ck])));

    % --- per-timepoint SUVR on the surface ---
    SUVR=nan(size(Vmid,1),nT); ref=zeros(1,nT);
    for i=1:nT
        Vi=double(niftiread(fullfile(subjDir,d(i).name)));
        ref(i)=mean(interpn(Vi,voxCb(:,1),voxCb(:,2),voxCb(:,3),'linear',NaN),'omitnan');
        SUVR(:,i)=interpn(Vi,voxCtx(:,1),voxCtx(:,2),voxCtx(:,3),'linear',NaN)/ref(i);
    end

    % --- epicenter per timepoint (on the ref-subject manifold cortex) ---
    th=sqrt(sum((sPm.Vertices-sWm.Vertices).^2,2))*1000; medial=th<1;
    reg=cell(1,nT); domSUVR=zeros(1,nT);
    for i=1:nT
        s=SUVR(:,i); s(medial)=NaN;
        foci=pet_epicenter(wf, s, struct('HeatT',2e-5));
        if isempty(foci), reg{i}='-'; continue; end
        reg{i}=local_region(sWm, foci(1).vertex); domSUVR(i)=s(foci(1).vertex);
    end
    gm=mean(SUVR(~medial,:),1,'omitnan');
    fprintf('\nADNI %s (%s anatomy, TEST): %d timepoints over %.1f yr\n', local_id(subjDir), refSubj, nT, yrs(end));
    fprintf('  global cortical SUVR: %s\n', mat2str(round(gm,2)));
    fprintf('  ref (cerebellum) intensity: %s\n', mat2str(round(ref,2)));
    for i=1:nT, fprintf('   %s (+%.1fyr): globalSUVR=%.2f  dominant focus=%s (SUVR %.2f)\n', dates{i}, yrs(i), gm(i), reg{i}, domSUVR(i)); end

    % --- figure: SUVR per timepoint + accumulation (medial LH) ---
    inLH=false(size(sWm.Vertices,1),1); LBO=tess_operators(wf,'Laplace-Beltrami'); gv=double(LBO.GlobalVertices{1}(:)); inLH(gv)=true;
    keep=all(inLH(sWm.Faces),2); rmp=zeros(size(sWm.Vertices,1),1); rmp(gv)=1:numel(gv); Flh=rmp(sWm.Faces(keep,:)); Vlh=sWm.Vertices(gv,:);
    f=figure('Visible','off','Position',[30 30 260*(nT+1) 300]);
    for i=1:nT
        subplot(1,nT+1,i); cd=SUVR(gv,i);
        patch('Faces',Flh,'Vertices',Vlh,'FaceVertexCData',cd,'FaceColor','interp','EdgeColor','none');
        clim([1 2]); colormap(hot); view([180 -10]); axis equal off vis3d; camlight headlight; lighting gouraud;
        title(sprintf('%s\n+%.1fyr',dates{i},yrs(i)));
    end
    subplot(1,nT+1,nT+1); acc=SUVR(gv,nT)-SUVR(gv,1);
    patch('Faces',Flh,'Vertices',Vlh,'FaceVertexCData',acc,'FaceColor','interp','EdgeColor','none');
    clim([-0.3 0.3]); colormap(gca,jet); view([180 -10]); axis equal off vis3d; camlight headlight; lighting gouraud;
    title(sprintf('accumulation\n(%.1fyr)',yrs(end)));
    png=fullfile(here,sprintf('proto_adni_%s.png',local_id(subjDir))); print(f,png,'-dpng','-r100'); close(f);
    fprintf('  figure -> %s\n', png);
    R=struct('SUVR',SUVR,'yrs',yrs,'dates',{dates},'globalSUVR',gm,'ref',ref,'domRegion',{reg},'png',png);
end

function id=local_id(subjDir)
    [~,id]=fileparts(subjDir); if isempty(id), [~,id]=fileparts(subjDir(1:end-1)); end
end
function reg=local_region(sW, v)
    reg='unknown'; ai=find(~cellfun('isempty',regexp({sW.Atlas.Name},'Desikan','once')),1); if isempty(ai), return; end
    for i=1:numel(sW.Atlas(ai).Scouts), if any(sW.Atlas(ai).Scouts(i).Vertices==v), reg=sW.Atlas(ai).Scouts(i).Label; return; end; end
end
