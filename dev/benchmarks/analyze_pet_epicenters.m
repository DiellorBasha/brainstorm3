function R = analyze_pet_epicenters(subjects, tracer)
% ANALYZE_PET_EPICENTERS: cohort epicenter detection + validation vs known seeding regions.
% Per subject: surface SUVR (mid-cortex sampling) -> pet_epicenter -> dominant focus region
% (Desikan); cohort focus density on the template (registration spheres); fraction of dominant
% foci in the expected regions (amyloid: precuneus/posteriorcingulate/isthmuscingulate/
% medialorbitofrontal; tau: entorhinal/inferiortemporal/fusiform); multifocality (# persistent foci).
%
% USAGE: R = analyze_pet_epicenters([], '18FNAV4694')        % whole cohort, amyloid
%        R = analyze_pet_epicenters({'sub-MTL0166',...}, '18FNAV4694')
%
% Author: Diellor Basha, 2026
    if (nargin<2)||isempty(tracer), tracer='18FNAV4694'; end
    if strcmp(tracer,'18FNAV4694')
        % Centiloid/amyloid cortical composite: DMN + frontal/parietal/lateral-temporal association
        % cortex (excludes amyloid-spared sensorimotor & primary visual: pre/postcentral,
        % paracentral, pericalcarine, cuneus, lingual).
        pet='PET 18FNAV4694_mean_pvc';
        expect={'precuneus','posteriorcingulate','isthmuscingulate','caudalanteriorcingulate', ...
                'rostralanteriorcingulate','medialorbitofrontal','lateralorbitofrontal','superiorfrontal', ...
                'rostralmiddlefrontal','caudalmiddlefrontal','parsopercularis','parsorbitalis','parstriangularis', ...
                'inferiorparietal','superiorparietal','supramarginal','middletemporal','superiortemporal'};
    else
        % Braak / tau-vulnerable medial & inferior temporal cortex.
        pet='PET 18Fflortaucipir_mean_pvc';
        expect={'entorhinal','inferiortemporal','fusiform','parahippocampal','middletemporal','temporalpole'};
    end
    if (nargin<1)||isempty(subjects)
        ps=bst_get('ProtocolSubjects'); nm={ps.Subject.Name};
        subjects=nm(~cellfun('isempty',regexp(nm,'^sub-MTL\d+$','once')));
    end
    here=bst_fileparts(mfilename('fullpath'));
    sDef=bst_get('Subject',0);
    tWf=sDef.Surface(find(~cellfun('isempty',regexp({sDef.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
    tW=in_tess_bst(tWf); dens=zeros(size(tW.Vertices,1),1);
    R=struct('subj',{},'tracer',{},'domRegion',{},'domPersist',{},'domPeak',{},'nFoci',{},'inExpected',{});
    for s=1:numel(subjects)
        subj=subjects{s}; [sS,~]=bst_get('Subject',subj); cmt={sS.Anatomy.Comment};
        af=@(c) sS.Anatomy(find(strcmp(cmt,c),1)).FileName;
        sf=@(rx) sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},rx,'once')),1)).FileName;
        if ~any(strcmp(cmt,pet)) || ~any(strcmp(cmt,'ASEG')), continue; end
        try
            wf=sf('cortex_white_low\.mat$'); pf=sf('cortex_pial_low\.mat$');
            sW=in_tess_bst(wf); sP=in_tess_bst(pf); sSuvr=pet_suvr(in_mri_bst(af(pet)),in_mri_bst(af('ASEG')));
            Vmid=0.5*(sW.Vertices+sP.Vertices); vox=cs_convert(sSuvr,'scs','voxel',Vmid);
            cs=size(sSuvr.Cube); ok=all(vox>=1,2)&vox(:,1)<=cs(1)&vox(:,2)<=cs(2)&vox(:,3)<=cs(3);
            suvr=nan(size(Vmid,1),1); suvr(ok)=interpn(double(sSuvr.Cube),vox(ok,1),vox(ok,2),vox(ok,3),'linear',NaN);
            th=sqrt(sum((sP.Vertices-sW.Vertices).^2,2))*1000; suvr(th<1)=NaN;
            foci=pet_epicenter(wf, suvr, struct('HeatT',2e-5));
            if isempty(foci), continue; end
            dv=foci(1).vertex; reg=local_region(sW, dv);
            inExp=any(cellfun(@(e) ~isempty(strfind(lower(reg),e)), expect));   %#ok<STREMP>
            R(end+1)=struct('subj',subj,'tracer',tracer,'domRegion',reg,'domPersist',foci(1).persistence, ...
                            'domPeak',foci(1).peak,'nFoci',numel(foci),'inExpected',inExp); %#ok<AGROW>
            Wm=tess_interp_tess2tess(wf, tWf, 0, 0, 0);          % map foci to template for the density map
            ind=zeros(size(sW.Vertices,1),1); ind([foci.vertex])=1; dens=dens+(Wm*ind > 0.5);
        catch ME, fprintf('  %s ERR %s\n',subj,ME.message); end
    end
    frac=mean([R.inExpected]);
    fprintf('\n%s: %d subjects | dominant focus in expected region: %.0f%% | median #foci=%.1f\n', ...
        tracer, numel(R), 100*frac, median([R.nFoci]));
    regs={R.domRegion}; ureg=unique(regs); cnt=cellfun(@(u)sum(strcmp(regs,u)), ureg);
    [cnt,o]=sort(cnt,'descend');
    fprintf('  top dominant-focus regions:\n');
    for k=1:min(8,numel(ureg)), fprintf('     %-28s %d\n', ureg{o(k)}, cnt(k)); end
    % figure: cohort focus-density on the template white (LH lateral + medial)
    try
        [~,lH]=tess_hemisplit(tW); idx=double(lH(:)); keep=all(ismember(tW.Faces,idx),2);
        rmp=zeros(size(tW.Vertices,1),1); rmp(idx)=1:numel(idx); Fh=rmp(tW.Faces(keep,:)); vw={[180 -10],[0 -10]};
        f=figure('Visible','off','Position',[40 40 900 420]);
        for s2=1:2
            subplot(1,2,s2);
            patch('Faces',Fh,'Vertices',tW.Vertices(idx,:),'FaceVertexCData',dens(idx),'FaceColor','interp','EdgeColor','none');
            clim([0 max(dens(idx))+eps]); colormap(hot); view(vw{s2}); axis equal off vis3d; camlight headlight; lighting gouraud;
            vwn={'lateral','medial'}; title(sprintf('%s focus density (LH %s)',tracer,vwn{s2}),'Interpreter','none');
        end
        png=fullfile(here,sprintf('pet_epicenters_%s.png',tracer)); print(f,png,'-dpng','-r110'); close(f);
        fprintf('  figure -> %s\n', png);
    catch ME, fprintf('  figure skipped: %s\n', ME.message); end
    R=struct('rows',{R},'density',dens);
end

function reg=local_region(sW, v)
    reg='unknown'; if ~isfield(sW,'Atlas')||isempty(sW.Atlas), return; end
    % Desikan-Killiany specifically (the surface also carries Destrieux/Structures atlases).
    ai=find(~cellfun('isempty',regexp({sW.Atlas.Name},'Desikan','once')),1);
    if isempty(ai), ai=find(~cellfun('isempty',regexp({sW.Atlas.Name},'aparc','once')),1); end
    if isempty(ai), return; end
    for i=1:numel(sW.Atlas(ai).Scouts)
        if any(sW.Atlas(ai).Scouts(i).Vertices==v), reg=sW.Atlas(ai).Scouts(i).Label; return; end
    end
end
