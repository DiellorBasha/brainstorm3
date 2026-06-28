function [Roi, Group] = preventad_pet_pipeline(Opts)
% PREVENTAD_PET_PIPELINE: full PET pipeline + ROI group statistics for the PREVENT-AD cohort.
%
% Runs the validated recommended pipeline for every subject that has (or can import) PET:
%   [import PET if missing]  ->  PVC (PETPVE12 Mueller-Gartner, auto-FWHM)  ->  robust SUVR
%   (pet_suvr, eroded+trimmed cerebellar cortex)  ->  regional Desikan SUVR (ROI table).
% Then computes GROUP ROI STATISTICS (per-region mean/SD/CV across subjects, per-subject global
% cortical SUVR, outlier / likely-amyloid-positive flag) and writes a CSV per tracer.
%
% Resumable + robust: reuses an existing _mean_pvc, skips subjects without PET (unless DoImport),
% and try/catches each subject so one failure does not abort the cohort.
%
% USAGE:  [Roi, Group] = preventad_pet_pipeline()              % all subjects with PET, MG, no import
%         preventad_pet_pipeline(struct('DoImport',1))         % import missing PET first, then process
%         preventad_pet_pipeline(struct('Subjects',{{'sub-MTL0002'}},'PvcMethod','gtm'))
%
% Opts (all optional; defaults shown):
%   .Subjects   []           subject ids to run ([] = all sub-MTL* in the protocol)
%   .Tracers    {amyloid,tau}'18FNAV4694','18Fflortaucipir'
%   .DoImport   0            import PET (preventad_pet_import) for subjects whose base is missing
%   .BidsPetDir (SpikeData-2)PET BIDS root for the import
%   .PvcMethod  'mg'         'mg' (pet_pvc/PETPVE12) or 'gtm' (pet_gtm)
%   .RefLabels  [8 47]       cerebellar-cortex reference (ASEG)
%   .OutDir     research/results/preventad_pet   where the CSVs are written
%
% OUTPUTS:
%   Roi   : struct array (one per subject x tracer) .Subject .Tracer .Ref .Global .Regions(1x68) ...
%   Group : per-tracer struct .Tracer .Names .Mean .SD .CV .GlobalSUVR .Outliers
%
% SEE ALSO: preventad_pet_import, pet_pvc, pet_gtm, pet_suvr, group_surface_suvr
%
% Author: Diellor Basha, 2026

    if (nargin<1)||isempty(Opts), Opts=struct(); end
    Def=struct('Subjects',[], 'Tracers',{{'18FNAV4694','18Fflortaucipir'}}, 'DoImport',0, ...
               'BidsPetDir','/Volumes/SpikeData-2/workspace/library/datasets/preventad/pet', ...
               'PvcMethod','mg', 'RefLabels',[8 47], ...
               'OutDir','/Users/diellorbasha/workspace/research/results/preventad_pet');
    fn=fieldnames(Def); for i=1:numel(fn), if ~isfield(Opts,fn{i})||isempty(Opts.(fn{i})), Opts.(fn{i})=Def.(fn{i}); end; end
    if ~exist(Opts.OutDir,'dir'), mkdir(Opts.OutDir); end

    if isempty(Opts.Subjects)
        ps=bst_get('ProtocolSubjects'); nm={ps.Subject.Name};
        Opts.Subjects=nm(~cellfun('isempty',regexp(nm,'^sub-MTL\d+$','once')));
    end
    pvcSuffix = struct('mg','_mean_pvc','gtm','_mean_gtmpvc'); pvcSuffix=pvcSuffix.(lower(Opts.PvcMethod));

    Roi=struct('Subject',{},'Tracer',{},'Ref',{},'Global',{},'Regions',{},'RegionNames',{});
    nDone=0; nSkip=0; nFail=0; T0=tic;
    for si=1:numel(Opts.Subjects)
        subj=Opts.Subjects{si}; [sS,~]=bst_get('Subject',subj); if isempty(sS), continue; end
        for ti=1:numel(Opts.Tracers)
            trc=Opts.Tracers{ti};
            try
                [sS,~]=bst_get('Subject',subj); cmt={sS.Anatomy.Comment};
                % --- import PET if missing (optional) ---
                if ~any(strcmp(cmt,['PET ' trc '_mean'])) && ~any(strcmp(cmt,['PET ' trc])) && Opts.DoImport
                    preventad_pet_import(Opts.BidsPetDir, subj);
                    [sS,~]=bst_get('Subject',subj); cmt={sS.Anatomy.Comment};
                end
                baseC=['PET ' trc '_mean'];
                if ~any(strcmp(cmt,baseC)), nSkip=nSkip+1; fprintf('  %-14s %-16s SKIP (no PET)\n',subj,trc); continue; end
                af=@(c) sS.Anatomy(find(strcmp(cmt,c),1)).FileName;

                % --- PVC (reuse or run) ---
                pvcC=['PET ' trc pvcSuffix];
                if ~any(strcmp(cmt,pvcC))
                    if strcmpi(Opts.PvcMethod,'gtm'), pet_gtm(af(baseC),[],struct());
                    else, pet_pvc(af(baseC), af('MRI T1'), [], struct()); end
                    [sS,~]=bst_get('Subject',subj); cmt={sS.Anatomy.Comment}; af=@(c) sS.Anatomy(find(strcmp(cmt,c),1)).FileName;
                end

                % --- robust SUVR + regional Desikan ROI ---
                sAseg=in_mri_bst(af('ASEG')); sDK=in_mri_bst(af('Desikan-Killiany'));
                [sSuvr,info]=pet_suvr(in_mri_bst(af(pvcC)), sAseg, struct('RefLabels',Opts.RefLabels));
                [names,vals]=local_regional(sDK, double(sSuvr.Cube(:,:,:,1)));
                Roi(end+1)=struct('Subject',subj,'Tracer',trc,'Ref',info.RefValue,'Global',mean(vals,'omitnan'), ...
                                  'Regions',vals,'RegionNames',{names}); %#ok<AGROW>
                nDone=nDone+1; fprintf('  %-14s %-16s OK   ref=%.2f globalSUVR=%.3f\n',subj,trc,info.RefValue,Roi(end).Global);
            catch ME
                nFail=nFail+1; fprintf('  %-14s %-16s FAIL %s\n',subj,trc,ME.message);
            end
        end
    end
    fprintf('\nprocessed %d (subject,tracer), skipped %d, failed %d in %.0fs\n', nDone,nSkip,nFail,toc(T0));

    % ===== GROUP ROI STATISTICS (per tracer) + CSV =====
    Group=struct('Tracer',{},'Names',{},'Mean',{},'SD',{},'CV',{},'GlobalSUVR',{},'Outliers',{});
    for ti=1:numel(Opts.Tracers)
        trc=Opts.Tracers{ti}; m=strcmp({Roi.Tracer},trc); if ~any(m), continue; end
        Rs=Roi(m); names=Rs(1).RegionNames; M=cat(1,Rs.Regions); subjs={Rs.Subject}; glob=[Rs.Global]';
        mu=mean(M,1,'omitnan'); sd=std(M,0,1,'omitnan'); cv=sd./mu;
        gmu=mean(glob,'omitnan'); gsd=std(glob,'omitnan'); outl=glob > gmu+2*gsd;     % likely amyloid+ (or high binding)
        Group(end+1)=struct('Tracer',trc,'Names',{names},'Mean',mu,'SD',sd,'CV',cv,'GlobalSUVR',glob,'Outliers',{subjs(outl)}); %#ok<AGROW>
        % per-subject ROI table CSV
        csv=fullfile(Opts.OutDir, sprintf('suvr_roi_%s.csv', trc));
        local_writecsv(csv, subjs, glob, names, M);
        % group-stats CSV
        scsv=fullfile(Opts.OutDir, sprintf('suvr_roi_%s_groupstats.csv', trc));
        local_writestats(scsv, names, mu, sd, cv);
        fprintf('\n[%s] n=%d subjects | global cortical SUVR %.3f +/- %.3f | outliers(>+2SD): %s\n', ...
            trc, numel(subjs), gmu, gsd, strjoin(subjs(outl),', '));
        fprintf('  ROI table -> %s\n  group stats -> %s\n', csv, scsv);
    end
end

% ===== regional Desikan SUVR (volume) =====
function [names,vals]=local_regional(sDK, cube)
    L=sDK.Labels; v=cell2mat(L(:,1)); nm=L(:,2); k=find(v>=1000 & v<3000);
    names=nm(k)'; vals=zeros(1,numel(k));
    for i=1:numel(k), vals(i)=mean(cube(sDK.Cube==v(k(i))),'omitnan'); end
end
function local_writecsv(file, subjs, glob, names, M)
    fid=fopen(file,'w'); fprintf(fid,'subject,global_cortical_suvr'); fprintf(fid,',%s',names{:}); fprintf(fid,'\n');
    for i=1:numel(subjs), fprintf(fid,'%s,%.4f',subjs{i},glob(i)); fprintf(fid,',%.4f',M(i,:)); fprintf(fid,'\n'); end
    fclose(fid);
end
function local_writestats(file, names, mu, sd, cv)
    fid=fopen(file,'w'); fprintf(fid,'region,mean_suvr,sd,cv\n');
    for i=1:numel(names), fprintf(fid,'%s,%.4f,%.4f,%.4f\n',names{i},mu(i),sd(i),cv(i)); end
    fclose(fid);
end
