function R = validate_surface_pvc_vs_petsurfer(subjects, tracers, psfFwhm)
% VALIDATE_SURFACE_PVC_VS_PETSURFER: cohort comparison of pet_pvc_surface (surface radial MG)
% against the PETSurfer regional GTM, over Desikan cortical regions.
%
% Per subject/tracer: aggregate our surface PVC and observed to the surface Desikan scouts;
% load PETSurfer gtm.nii.gz (GTM PVC) + nopvc.nii.gz (observed), matched by region name.
%   r_obs : corr(ourObs, psObs)             - sampling-difference baseline (surface mid vs volume ROI)
%   r_pvc : corr(ourPVC, psGTM)             - absolute pattern agreement after PVC
%   r_eff : corr(ourPVC/ourObs, psGTM/psObs)- PVC-EFFECT agreement (controls for the baseline)
%   boost : mean(PVC/obs) for each method   - direction/magnitude of the correction
%
% Exploration/validation only.
%
% Author: Diellor Basha, 2026

    PsRoot='/Volumes/SpikeData-2/workspace/library/datasets/preventad/pet/derivatives/petsurfer';
    if (nargin<1)||isempty(subjects)
        ps=bst_get('ProtocolSubjects'); names={ps.Subject.Name};
        subjects=names(~cellfun('isempty',regexp(names,'^sub-MTL\d+$','once'))); subjects=subjects(1:min(10,end));
    end
    if (nargin<2)||isempty(tracers), tracers={'18FNAV4694','18Fflortaucipir'}; end
    if (nargin<3)||isempty(psfFwhm), psfFwhm=2.5; end
    c=@(a,b) local_corr(a,b);
    rows={}; R=struct('subj',{},'tracer',{},'r_obs',{},'r_pvc',{},'r_eff',{},'boost_our',{},'boost_gtm',{},'n',{});
    for si=1:numel(subjects)
        for ti=1:numel(tracers)
            subj=subjects{si}; trc=tracers{ti};
            d=fullfile(PsRoot,subj,'ses-01',trc,'gtmpvc.output');
            if ~exist(fullfile(d,'gtm.nii.gz'),'file'), continue; end
            try
                [sS,~]=bst_get('Subject',subj);
                fn=@(cc) sS.Anatomy(find(strcmp({sS.Anatomy.Comment},cc),1)).FileName;
                wf=sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
                if isempty(fn(['PET ' trc '_mean'])), continue; end
                sPet=in_mri_bst(fn(['PET ' trc '_mean'])); sT1=in_mri_bst(fn('MRI T1'));
                [GM,info]=pet_pvc_surface(sPet,sT1,wf,struct('PsfFwhm',psfFwhm));
                sW=in_tess_bst(wf); A=sW.Atlas; ai=find(~cellfun('isempty',regexp({A.Name},'Desikan|aparc','once')),1); sc=A(ai).Scouts;
                fid=fopen(fullfile(d,'gtm.stats.dat')); G=textscan(fid,'%f %f %s %s %f %f %f %f','CommentStyle','#'); fclose(fid);
                Tn=bst_get('BrainstormTmpDir',0,'pvccmp'); rd=@(f)squeeze(mean(double(niftiread(f)),4));
                copyfile(fullfile(d,'gtm.nii.gz'),fullfile(Tn,'g.nii.gz')); gunzip(fullfile(Tn,'g.nii.gz')); psGtm=rd(fullfile(Tn,'g.nii'));
                copyfile(fullfile(d,'nopvc.nii.gz'),fullfile(Tn,'n.nii.gz')); gunzip(fullfile(Tn,'n.nii.gz')); psObs=rd(fullfile(Tn,'n.nii'));
                file_delete(Tn,1,1);
                ctx=find(~cellfun('isempty',regexp(G{3},'^ctx-[lr]h-','once')));
                oN=lower(regexprep({sc.Label},'\s*([LR])$','_$1')); pN=lower(regexprep(G{3}(ctx),'^ctx-(l|r)h-(.*)$','$2_$1'));
                oPvc=cellfun(@(v)mean(GM(v),'omitnan'),{sc.Vertices}); oObs=cellfun(@(v)mean(info.observed(v),'omitnan'),{sc.Vertices});
                [tf,loc]=ismember(pN,oN);
                oo=oObs(loc(tf))'; op=oPvc(loc(tf))'; pp=psObs(ctx(tf)); gg=psGtm(ctx(tf));
                R(end+1)=struct('subj',subj,'tracer',trc,'r_obs',c(oo,pp),'r_pvc',c(op,gg), ...
                    'r_eff',c(op./oo,gg./pp),'boost_our',mean(op./oo),'boost_gtm',mean(gg./pp),'n',nnz(tf)); %#ok<AGROW>
                rows{end+1}=sprintf('%-14s %-16s n=%d  r_obs=%.2f  r_pvc=%.2f  r_eff=%+.2f  boost ours=%.2f gtm=%.2f', ...
                    subj,trc,nnz(tf),R(end).r_obs,R(end).r_pvc,R(end).r_eff,R(end).boost_our,R(end).boost_gtm); %#ok<AGROW>
            catch ME, rows{end+1}=sprintf('%-14s %-16s ERR %s',subj,trc,ME.message); end %#ok<AGROW>
        end
    end
    fprintf('\n%s\n', strjoin(rows, char(10)));
    for ti=1:numel(tracers)
        m=strcmp({R.tracer},tracers{ti}); if ~any(m), continue; end
        fprintf('--- %s (n=%d): r_obs=%.2f  r_pvc=%.2f  r_eff=%+.2f  boost ours=%.2f gtm=%.2f\n', ...
            tracers{ti}, nnz(m), mean([R(m).r_obs]), mean([R(m).r_pvc]), mean([R(m).r_eff]), mean([R(m).boost_our]), mean([R(m).boost_gtm]));
    end
end

function r=local_corr(a,b)
    a=a(:); b=b(:); ok=isfinite(a)&isfinite(b); if nnz(ok)<5, r=NaN; return; end
    cc=corrcoef(a(ok),b(ok)); r=cc(1,2);
end
