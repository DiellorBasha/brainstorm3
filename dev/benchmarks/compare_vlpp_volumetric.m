function R = compare_vlpp_volumetric(subjects, tracer)
% COMPARE_VLPP_VOLUMETRIC: voxel-wise comparison of our SUVR volumes vs the VLPP SUVR volumes
% (pet/derivatives/vlpp/<subj>/ses-01/..._suvr.nii.gz). VLPP volumes are FreeSurfer-conformed LIA;
% they are resampled onto our Brainstorm T1 grid via world coordinates (our world[m]x1000 = VLPP
% RAS[mm]; verified r=0.80 at matched resolution). Compared within the Desikan cortical mask.
%
% Reports per-subject voxel-wise correlation at NATIVE resolution (our sharp PVC'd vs VLPP 6mm-
% smoothed no-PVC) and at MATCHED resolution (ours smoothed to 6mm), plus bias.
%
% USAGE: R = compare_vlpp_volumetric([], '18FNAV4694')   % amyloid (ref cerebellumCortex)
%        R = compare_vlpp_volumetric([], '18Fflortaucipir')
%
% Author: Diellor Basha, 2026

    vlppRoot='/Volumes/SpikeData-2/workspace/library/datasets/preventad/pet/derivatives/vlpp';
    here=bst_fileparts(mfilename('fullpath'));
    if strcmp(tracer,'18FNAV4694'), ref='cerebellumCortex'; petC='PET 18FNAV4694_mean_pvc'; lbl='amyloid';
    else, ref='infCerebellarGray'; petC='PET 18Fflortaucipir_mean_pvc'; lbl='tau'; end
    if (nargin<1)||isempty(subjects)
        ps=bst_get('ProtocolSubjects'); nm={ps.Subject.Name}; subjects=nm(~cellfun('isempty',regexp(nm,'^sub-MTL\d+$','once')));
    end
    R=struct('subj',{},'rNative',{},'r6mm',{},'bias',{},'mOur',{},'mVlpp',{});
    sig=6/2.355; egV=[]; egO=[]; egZ=[];
    for s=1:numel(subjects)
        subj=subjects{s};
        vf=fullfile(vlppRoot,subj,'ses-01',sprintf('%s_ses-01_trc-%s_pet_ref-%s_suvr.nii.gz',subj,tracer,ref));
        if ~exist(vf,'file'), continue; end
        try
            [sS,~]=bst_get('Subject',subj); cmt={sS.Anatomy.Comment}; af=@(c) sS.Anatomy(find(strcmp(cmt,c),1)).FileName;
            if ~any(strcmp(cmt,petC)), continue; end
            Tn=bst_get('BrainstormTmpDir',0,'vlv'); copyfile(vf,fullfile(Tn,'v.nii.gz')); gunzip(fullfile(Tn,'v.nii.gz'));
            ni=niftiinfo(fullfile(Tn,'v.nii')); Vlpp=double(niftiread(ni)); T=ni.Transform.T'; M=T(1:3,1:3); b=T(1:3,4); file_delete(Tn,1,1);
            sT1=in_mri_bst(af('MRI T1')); sAseg=in_mri_bst(af('ASEG')); sDK=in_mri_bst(af('Desikan-Killiany'));
            sSuvr=pet_suvr(in_mri_bst(af(petC)),sAseg); ourC=double(sSuvr.Cube); ourSm=imgaussfilt3(ourC,sig);
            ctx=find(sDK.Cube>=1000 & sDK.Cube<3000); [i,j,k]=ind2sub(size(ourC),ctx);
            wmm=cs_convert(sT1,'voxel','world',[i j k])*1000; vox1=(wmm-b')/M'+1;
            vl=interpn(Vlpp,vox1(:,1),vox1(:,2),vox1(:,3),'linear',NaN);
            o=ourC(ctx); osm=ourSm(ctx); ok=isfinite(vl)&isfinite(o)&vl>0;
            cN=corrcoef(o(ok),vl(ok)); c6=corrcoef(osm(ok),vl(ok));
            R(end+1)=struct('subj',subj,'rNative',cN(1,2),'r6mm',c6(1,2),'bias',mean(o(ok)-vl(ok)),'mOur',mean(o(ok)),'mVlpp',mean(vl(ok))); %#ok<AGROW>
            if isempty(egO)   % keep one subject's middle slices for a difference figure
                z=round(size(ourC,3)/2); egO=ourSm(:,:,z); volR=nan(size(ourC)); volR(ctx)=vl; egV=volR(:,:,z); egZ=subj;
            end
        catch ME, fprintf('  %s ERR %s\n',subj,ME.message); end
    end
    rN=[R.rNative]; r6=[R.r6mm]; bi=[R.bias];
    fprintf('\n=== %s : VOLUMETRIC ours vs VLPP (%d subjects, cortical voxels) ===\n', lbl, numel(R));
    fprintf('  voxel-wise r NATIVE (sharp PVC vs 6mm smoothed): %.3f +/- %.3f\n', mean(rN),std(rN));
    fprintf('  voxel-wise r MATCHED (ours->6mm)               : %.3f +/- %.3f\n', mean(r6),std(r6));
    fprintf('  bias (ours-vlpp): %+.3f +/- %.3f | mean SUVR ours=%.2f vlpp=%.2f\n', mean(bi),std(bi),mean([R.mOur]),mean([R.mVlpp]));

    f=figure('Visible','off','Position',[50 50 1300 430]);
    subplot(1,3,1); histogram(rN,15,'FaceColor',[.8 .4 .2]); hold on; histogram(r6,15,'FaceColor',[.2 .5 .8]); xlabel('voxel-wise r'); ylabel('# subjects'); legend({'native','matched 6mm'},'Location','northwest'); title(sprintf('%s: per-subject voxel-wise r',lbl)); grid on;
    if ~isempty(egO)
      subplot(1,3,2); imagesc(rot90(egO),[0.8 2]); axis image off; colormap(gca,hot); title(sprintf('OURS (PVC) %s',egZ),'Interpreter','none');
      subplot(1,3,3); imagesc(rot90(egV),[0.8 2]); axis image off; colormap(gca,hot); title('VLPP (6mm, no PVC)');
    end
    png=fullfile(here,sprintf('compare_vlpp_vol_%s.png',lbl)); print(f,png,'-dpng','-r110'); close(f);
    fprintf('  figure -> %s\n', png);
end
