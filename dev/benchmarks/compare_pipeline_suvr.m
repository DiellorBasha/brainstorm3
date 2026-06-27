function R = compare_pipeline_suvr(subjects, tracers)
% COMPARE_PIPELINE_SUVR: cohort comparison of cortical SUVR across PVC pipelines.
%   OURS-MG  = pet_pvc (PETPVE12 Mueller-Gaertner, _mean_pvc) + robust cerebellar SUVR (pet_suvr)
%   OURS-GTM = pet_gtm (native Rousset, _mean_gtmpvc)         + robust cerebellar SUVR (pet_suvr)
%   PETSURFER= gtm.nii (already cerebellum-rescaled -> SUVR), regional GTM
% over the 68 Desikan cortical regions, matched by FreeSurfer label number (DK volume 1001.. ==
% PETSurfer segid). Reports correlation (pattern), CCC (agreement), and bias.
%
% Author: Diellor Basha, 2026

    PsRoot='/Volumes/SpikeData-2/workspace/library/datasets/preventad/pet/derivatives/petsurfer';
    if (nargin<1)||isempty(subjects)
        ps=bst_get('ProtocolSubjects'); nm={ps.Subject.Name};
        subjects=nm(~cellfun('isempty',regexp(nm,'^sub-MTL\d+$','once'))); subjects=subjects(1:min(10,end));
    end
    if (nargin<2)||isempty(tracers), tracers={'18FNAV4694','18Fflortaucipir'}; end
    R=struct('subj',{},'tracer',{},'mg',{},'gtm',{},'ps',{});
    here=bst_fileparts(mfilename('fullpath'));
    for si=1:numel(subjects)
        subj=subjects{si}; [sS,e]=bst_get('Subject',subj); if isempty(sS), continue; end
        an=@(c) local_load(sS,c);
        sDK=an('Desikan-Killiany'); sAseg=an('ASEG'); if isempty(sDK)||isempty(sAseg), continue; end
        DK=sDK.Cube; labs=cell2mat(sDK.Labels(:,1)); labs=labs(labs>=1000 & labs<3000);
        for ti=1:numel(tracers)
            trc=tracers{ti};
            sMG=an(['PET ' trc '_mean_pvc']); sGT=an(['PET ' trc '_mean_gtmpvc']);
            d=fullfile(PsRoot,subj,'ses-01',trc,'gtmpvc.output');
            if isempty(sMG)||~exist(fullfile(d,'gtm.nii.gz'),'file'), continue; end   % GTM optional
            % robust cerebellar reference per PVC volume
            [~,iMG]=pet_suvr(sMG,sAseg); refGT=NaN; if ~isempty(sGT), [~,iGT]=pet_suvr(sGT,sAseg); refGT=iGT.RefValue; end
            % PETSurfer GTM (already SUVR) by segid
            fid=fopen(fullfile(d,'gtm.stats.dat')); G=textscan(fid,'%f %f %s %s %f %f %f %f','CommentStyle','#'); fclose(fid);
            Tn=bst_get('BrainstormTmpDir',0,'cmp'); copyfile(fullfile(d,'gtm.nii.gz'),fullfile(Tn,'g.nii.gz')); gunzip(fullfile(Tn,'g.nii.gz'));
            psv=squeeze(mean(double(niftiread(fullfile(Tn,'g.nii'))),4)); file_delete(Tn,1,1);
            mg=nan(numel(labs),1); gt=mg; ps=mg;
            for k=1:numel(labs)
                m=(DK==labs(k)); if nnz(m)<10, continue; end
                mg(k)=mean(double(sMG.Cube(m)))/iMG.RefValue;
                if ~isempty(sGT), gt(k)=mean(double(sGT.Cube(m)))/refGT; end
                j=find(G{2}==labs(k),1); if ~isempty(j), ps(k)=psv(j); end
            end
            R(end+1)=struct('subj',subj,'tracer',trc,'mg',mg,'gtm',gt,'ps',ps); %#ok<AGROW>
        end
    end
    % ---- analysis ----
    cc=@(a,b)subsref(corrcoef(a,b),struct('type','()','subs',{{1,2}}));
    ccc=@(a,b) 2*local_cov(a,b)/(var(a,1)+var(b,1)+(mean(a)-mean(b))^2);   % Lin's concordance
    f=figure('Visible','off','Position',[50 50 900 430]);
    for ti=1:numel(tracers)
        m=strcmp({R.tracer},tracers{ti}); if ~any(m), continue; end
        MG=cat(1,R(m).mg); GT=cat(1,R(m).gtm); PS=cat(1,R(m).ps);
        oMP=isfinite(MG)&isfinite(PS); oGP=isfinite(GT)&isfinite(PS); oMG=isfinite(MG)&isfinite(GT);
        nGsub=nnz(arrayfun(@(x)any(isfinite(x.gtm)),R(m)));
        fprintf('\n=== %s (MG: %d subj/%d reg | GTM: %d subj/%d reg) ===\n', tracers{ti}, nnz(m), nnz(oMP), nGsub, nnz(oGP));
        fprintf('  OURS-MG(PETPVE12) vs PETSurfer-GTM : r=%.3f CCC=%.3f bias(our-ps)=%+.3f  meanSUVR our=%.2f ps=%.2f\n', cc(MG(oMP),PS(oMP)),ccc(MG(oMP),PS(oMP)),mean(MG(oMP)-PS(oMP)),mean(MG(oMP)),mean(PS(oMP)));
        if nnz(oGP)>5
        fprintf('  OURS-GTM          vs PETSurfer-GTM : r=%.3f CCC=%.3f bias(our-ps)=%+.3f  meanSUVR our=%.2f ps=%.2f\n', cc(GT(oGP),PS(oGP)),ccc(GT(oGP),PS(oGP)),mean(GT(oGP)-PS(oGP)),mean(GT(oGP)),mean(PS(oGP)));
        fprintf('  OURS-MG           vs OURS-GTM      : r=%.3f CCC=%.3f\n', cc(MG(oMG),GT(oMG)),ccc(MG(oMG),GT(oMG)));
        end
        subplot(1,numel(tracers),ti); scatter(PS(oMP),MG(oMP),8,'r','filled'); hold on; if nnz(oGP)>5, scatter(PS(oGP),GT(oGP),8,'b','filled'); end
        lim=[min(PS(oMP)) max(PS(oMP))]; plot(lim,lim,'k--'); axis equal; grid on; xlim(lim); ylim(lim);
        xlabel('PETSurfer-GTM SUVR'); ylabel('OURS SUVR'); title(tracers{ti},'Interpreter','none');
        legend({'ours-MG (PETPVE12)','ours-GTM','identity'},'Location','northwest');
    end
    print(f,fullfile(here,'compare_pipeline_suvr.png'),'-dpng','-r110'); close(f);
    fprintf('\nFigure: %s\n', fullfile(here,'compare_pipeline_suvr.png'));
end

function sM=local_load(sS,c)
    i=find(strcmp({sS.Anatomy.Comment},c),1); if isempty(i), sM=[]; else, sM=in_mri_bst(sS.Anatomy(i).FileName); end
end
function c=local_cov(a,b), m=mean([a b]); c=mean((a-m(1)).*(b-m(2))); end
