function R = compare_vlpp_brainstorm(ourCsv, vlppCsv, tracerLabel)
% COMPARE_VLPP_BRAINSTORM: compare our Brainstorm PET regional SUVR to the VLPP pipeline CSV.
% Matches the 68 Desikan cortical ROIs by FreeSurfer label (ours: DK volume labels 1001..2035;
% VLPP: ctx-{lh,rh}-<region>_SUVR columns) and subjects by PSCID (MTL####, ses-01).
%
% USAGE: R = compare_vlpp_brainstorm('.../suvr_roi_18FNAV4694.csv', '.../PET_NAV_..._cerebellumCortex.csv', 'amyloid (NAV)')
%
% Author: Diellor Basha, 2026

    here=bst_fileparts(mfilename('fullpath'));
    % FreeSurfer Desikan cortical region names in label order (lh 1001..1035, skipping 1004)
    fs={'bankssts','caudalanteriorcingulate','caudalmiddlefrontal','cuneus','entorhinal','fusiform', ...
        'inferiorparietal','inferiortemporal','isthmuscingulate','lateraloccipital','lateralorbitofrontal', ...
        'lingual','medialorbitofrontal','middletemporal','parahippocampal','paracentral','parsopercularis', ...
        'parsorbitalis','parstriangularis','pericalcarine','postcentral','posteriorcingulate','precentral', ...
        'precuneus','rostralanteriorcingulate','rostralmiddlefrontal','superiorfrontal','superiorparietal', ...
        'superiortemporal','supramarginal','frontalpole','temporalpole','transversetemporal','insula'};
    fsLab=[1001 1002 1003 1005 1006 1007 1008 1009 1010 1011 1012 1013 1014 1015 1016 1017 1018 1019 1020 ...
           1021 1022 1023 1024 1025 1026 1027 1028 1029 1030 1031 1032 1033 1034 1035];

    % our DK volume atlas: Brainstorm-name -> label (so our CSV columns map to labels)
    [sS,~]=bst_get('Subject','sub-MTL0002'); sDK=in_mri_bst(sS.Anatomy(find(strcmp({sS.Anatomy.Comment},'Desikan-Killiany'),1)).FileName);
    L=sDK.Labels; bnames=L(:,2); blab=cell2mat(L(:,1));
    name2lab=containers.Map(); for i=1:numel(bnames), name2lab(bnames{i})=blab(i); end

    Ours=readtable(ourCsv,'VariableNamingRule','preserve'); Vl=readtable(vlppCsv,'VariableNamingRule','preserve');
    ourSubj=Ours.subject; ourCols=Ours.Properties.VariableNames;
    vlPSCID=Vl.PSCID;

    % build label-indexed matrices: rows=matched subjects, cols=labels
    labs=[fsLab, fsLab+1000];                          % lh + rh = 68 cortical labels
    nL=numel(labs);
    % our column index per label
    ourColLab=nan(1,numel(ourCols)); for c=1:numel(ourCols), if isKey(name2lab,ourCols{c}), ourColLab(c)=name2lab(ourCols{c}); end; end
    % vlpp column name per label
    vlName=cell(1,nL);
    for k=1:nL
        li=labs(k); if li<2000, hemi='lh'; idx=find(fsLab==li,1); else, hemi='rh'; idx=find(fsLab==(li-1000),1); end
        vlName{k}=sprintf('ctx-%s-%s_SUVR',hemi,fs{idx});
    end

    OUR=[]; VLP=[]; subjUsed={};
    for s=1:numel(ourSubj)
        id=strrep(ourSubj{s},'sub-','');                % MTL0002
        vi=find(strcmp(vlPSCID,id),1); if isempty(vi), continue; end
        ourRow=nan(1,nL); vlRow=nan(1,nL);
        for k=1:nL
            oc=find(ourColLab==labs(k),1); if ~isempty(oc), ourRow(k)=Ours{s,oc}; end
            if ismember(vlName{k},Vl.Properties.VariableNames), vlRow(k)=Vl{vi,vlName{k}}; end
        end
        OUR(end+1,:)=ourRow; VLP(end+1,:)=vlRow; subjUsed{end+1}=id; %#ok<AGROW>
    end
    ok=isfinite(OUR)&isfinite(VLP);
    cc=@(a,b)subsref(corrcoef(a,b),struct('type','()','subs',{{1,2}}));
    o=OUR(ok); v=VLP(ok);
    ccc = 2*mean((o-mean(o)).*(v-mean(v))) / (var(o,1)+var(v,1)+(mean(o)-mean(v))^2);   % Lin's concordance
    fprintf('\n=== %s : OURS vs VLPP ===\n', tracerLabel);
    fprintf('  matched %d subjects, %d region-values\n', numel(subjUsed), nnz(ok));
    fprintf('  pooled regional SUVR  : r=%.3f  CCC=%.3f  bias(ours-vlpp)=%+.3f  meanSUVR ours=%.3f vlpp=%.3f\n', ...
        cc(o,v), ccc, mean(o-v), mean(o), mean(v));
    % per-subject global cortical SUVR
    og=mean(OUR,2,'omitnan'); vg=mean(VLP,2,'omitnan');
    fprintf('  per-subject GLOBAL cortical SUVR: r=%.3f  bias=%+.3f\n', cc(og,vg), mean(og-vg));
    R=struct('subj',{subjUsed},'OUR',OUR,'VLP',VLP,'labs',labs,'globalOur',og,'globalVlpp',vg);

    f=figure('Visible','off','Position',[60 60 900 430]);
    subplot(1,2,1); scatter(v,o,6,'filled'); hold on; lim=[min([o;v]) max([o;v])]; plot(lim,lim,'k--'); axis equal; grid on; xlim(lim);ylim(lim);
    xlabel('VLPP regional SUVR'); ylabel('OUR regional SUVR'); title(sprintf('%s: regional (r=%.3f)',tracerLabel,cc(o,v)));
    subplot(1,2,2); scatter(vg,og,18,'filled'); hold on; lim=[min([og;vg]) max([og;vg])]; plot(lim,lim,'k--'); axis equal; grid on; xlim(lim);ylim(lim);
    xlabel('VLPP global cortical SUVR'); ylabel('OUR global cortical SUVR'); title(sprintf('per-subject global (r=%.3f)',cc(og,vg)));
    pngf=fullfile(here,sprintf('compare_vlpp_%s.png',matlab.lang.makeValidName(tracerLabel)));
    print(f,pngf,'-dpng','-r110'); close(f); fprintf('  figure -> %s\n', pngf);
end
