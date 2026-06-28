function R = demo_vlpp_style(subjects, tracer, fwhm)
% DEMO_VLPP_STYLE: reproduce VLPP's NON-PVC SMOOTHED pipeline with existing Brainstorm functions
% and check it matches VLPP. Pipeline = the already-computed coregistered static 'PET <tracer>_mean'
% (mri_realign+mri_aggregate+mri_coregister+mri_reslice) -> 6mm Gaussian smooth (same as pet_gtm's
% local_gauss3) -> pet_suvr with a PLAIN cerebellar-cortex mean (no erosion, no PVC) -> Desikan ROI.
% Compared to the VLPP non-PVC SUVR CSV. Contrast with our standard MG-PVC SUVR (the +0.20 bias).
%
% Author: Diellor Basha, 2026

    if (nargin<2)||isempty(tracer), tracer='18FNAV4694'; end
    if (nargin<3)||isempty(fwhm), fwhm=6; end
    root='/Volumes/SpikeData-2/workspace/library/datasets/preventad';
    if strcmp(tracer,'18FNAV4694'), vcsv='PREVENT-AD_internal_n386__PET_NAV_SUVR_ref-cerebellumCortex.csv'; petC='PET 18FNAV4694_mean';
    else, vcsv='PREVENT-AD_internal_n386__PET_FTP_SUVR_ref-infCerebellarGray.csv'; petC='PET 18Fflortaucipir_mean'; end
    fs={'bankssts','caudalanteriorcingulate','caudalmiddlefrontal','cuneus','entorhinal','fusiform','inferiorparietal','inferiortemporal','isthmuscingulate','lateraloccipital','lateralorbitofrontal','lingual','medialorbitofrontal','middletemporal','parahippocampal','paracentral','parsopercularis','parsorbitalis','parstriangularis','pericalcarine','postcentral','posteriorcingulate','precentral','precuneus','rostralanteriorcingulate','rostralmiddlefrontal','superiorfrontal','superiorparietal','superiortemporal','supramarginal','frontalpole','temporalpole','transversetemporal','insula'};
    fsLab=[1001:1003 1005:1035]; labs=[fsLab fsLab+1000];
    if (nargin<1)||isempty(subjects)
        ps=bst_get('ProtocolSubjects'); nm={ps.Subject.Name}; subjects=nm(~cellfun('isempty',regexp(nm,'^sub-MTL\d+$','once')));
    end
    Vl=readtable(fullfile(root,vcsv),'VariableNamingRule','preserve'); vlPSCID=Vl.PSCID;
    vlName=cell(1,numel(labs));
    for k=1:numel(labs), li=labs(k); if li<2000, h='lh'; idx=find(fsLab==li,1); else, h='rh'; idx=find(fsLab==li-1000,1); end; vlName{k}=sprintf('ctx-%s-%s_SUVR',h,fs{idx}); end

    sig=fwhm/2.355; OURv=[]; OURp=[]; VLP=[]; used={};
    for s=1:numel(subjects)
        subj=subjects{s}; id=strrep(subj,'sub-','');
        vi=find(strcmp(vlPSCID,id),1); if isempty(vi), continue; end
        [sS,~]=bst_get('Subject',subj); cmt={sS.Anatomy.Comment}; af=@(c) sS.Anatomy(find(strcmp(cmt,c),1)).FileName;
        if ~any(strcmp(cmt,petC)) || ~any(strcmp(cmt,'ASEG')), continue; end
        try
            sMean=in_mri_bst(af(petC)); sAseg=in_mri_bst(af('ASEG')); sDK=in_mri_bst(af('Desikan-Killiany'));
            % --- VLPP-style: smooth the non-PVC static 6mm, SUVR with PLAIN cerebellar-cortex mean ---
            sSm=sMean; sSm.Cube=imgaussfilt3(double(sMean.Cube(:,:,:,1)),sig);
            [svV,~]=pet_suvr(sSm, sAseg, struct('Erode',0,'Robust','mean'));   % no erosion, no PVC
            % --- our standard MG-PVC SUVR (for contrast), if available ---
            cv=double(svV.Cube); ov=nan(1,numel(labs)); pv=nan(1,numel(labs)); vv=nan(1,numel(labs));
            hasPvc=any(strcmp(cmt,[petC '_pvc']));
            if hasPvc, svP=pet_suvr(in_mri_bst(af([petC '_pvc'])),sAseg); cp=double(svP.Cube); end
            for k=1:numel(labs)
                m=(sDK.Cube==labs(k)); if nnz(m)>10, ov(k)=mean(cv(m),'omitnan'); if hasPvc, pv(k)=mean(cp(m),'omitnan'); end; end
                if ismember(vlName{k},Vl.Properties.VariableNames), vv(k)=Vl{vi,vlName{k}}; end
            end
            OURv(end+1,:)=ov; OURp(end+1,:)=pv; VLP(end+1,:)=vv; used{end+1}=id; %#ok<AGROW>
        catch ME, fprintf('  %s ERR %s\n',subj,ME.message); end
    end
    cc=@(a,b)subsref(corrcoef(a,b),struct('type','()','subs',{{1,2}}));
    ov=isfinite(OURv)&isfinite(VLP); pv=isfinite(OURp)&isfinite(VLP);
    fprintf('\n=== %s : reproduce VLPP non-PVC smoothed pipeline (%d subjects) ===\n', tracer, numel(used));
    fprintf('  VLPP-STYLE (ours: smooth+SUVR, no PVC) vs VLPP : r=%.3f  bias=%+.3f  meanSUVR ours=%.3f vlpp=%.3f\n', cc(OURv(ov),VLP(ov)), mean(OURv(ov)-VLP(ov)), mean(OURv(ov)), mean(VLP(ov)));
    og=mean(OURv,2,'omitnan'); vg=mean(VLP,2,'omitnan');
    fprintf('  per-subject GLOBAL cortical SUVR              : r=%.3f  bias=%+.3f\n', cc(og,vg), mean(og-vg));
    if any(pv(:))
    fprintf('  (for contrast) our MG-PVC SUVR vs VLPP        : r=%.3f  bias=%+.3f\n', cc(OURp(pv),VLP(pv)), mean(OURp(pv)-VLP(pv)));
    end
    R=struct('subj',{used},'ourVlppStyle',OURv,'ourPvc',OURp,'vlpp',VLP);
end
