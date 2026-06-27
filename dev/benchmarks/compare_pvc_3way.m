function R = compare_pvc_3way(SubjectName, tracer, PvcComment, PsRoot)
% COMPARE_PVC_3WAY: Three-way regional-SUVR comparison of PET partial-volume correction:
%   (1) OURS   - PETPVE12 Muller-Gartner (voxelwise), computed in Brainstorm.
%   (2) MGX    - PETSurfer extended Muller-Gartner (voxelwise) -> mgx.gm.nii.gz.
%   (3) GTM    - PETSurfer Geometric Transfer Matrix (Rousset, regional) -> gtm.stats.dat.
%
% This gives both a CROSS-PIPELINE same-method check (OURS vs MGX, both Muller-Gartner)
% and a SAME-PIPELINE cross-method check (MGX vs GTM, both PETSurfer), plus the full
% cross OURS vs GTM. All values are regional means over FreeSurfer ROIs, each pipeline
% normalized to its OWN cerebellum cortex (aseg 8,47), so only the regional PATTERN is
% compared (scale conventions cancel).
%
% PETSurfer NIfTIs are read non-interactively with niftiread (NOT the Brainstorm
% importer). Our volume / atlases come from the Brainstorm DB.
%
% KNOWN WEAKNESS of OURS (recorded for future work): Muller-Gartner corrects only the GM
% compartment and assigns WM/CSF constant activity, so WM / CSF / subcortical ROIs are
% expected to diverge from GTM (which estimates every region). The cortical ROIs are the
% valid target of comparison.
%
% USAGE:  R = compare_pvc_3way('sub-MTL0002','18FNAV4694','PET 18FNAV4694_mean_pvc')
%
% Author: Diellor Basha, 2026

    if (nargin < 4) || isempty(PsRoot)
        PsRoot = '/Volumes/SpikeData-2/workspace/library/datasets/preventad/pet/derivatives/petsurfer';
    end
    OutDir = bst_fileparts(mfilename('fullpath'));
    psDir  = fullfile(PsRoot, SubjectName, 'ses-01', tracer, 'gtmpvc.output');

    % ===== OURS (Brainstorm) =====
    [sS,~] = bst_get('Subject', SubjectName);
    pvc  = local_cube(sS, PvcComment);
    aseg = local_cube(sS, 'ASEG');
    apar = local_cube(sS, 'Desikan-Killiany');
    hasApar = ~isempty(apar) && isequal(size(apar), size(pvc));
    ourSeg = @(id) local_seg(id, aseg, apar, hasApar);
    ourCb  = mean(pvc(ourSeg(8) | ourSeg(47)));

    % ===== GTM (gtm.stats.dat) =====
    G = local_read_gtm(fullfile(psDir,'gtm.stats.dat'));
    iCb = ismember(G.segId,[8 47]);
    gtmCb = sum(G.mean(iCb).*G.nvox(iCb))/sum(G.nvox(iCb));

    % ===== MGX (mgx.gm.nii.gz over aux/seg.nii.gz) =====
    mgx = mean(local_niigz(fullfile(psDir,'mgx.gm.nii.gz')), 4);   % average 6 frames -> static
    seg = local_niigz(fullfile(psDir,'aux','seg.nii.gz'));
    if ~isequal(size(mgx), size(seg)), error('MGX/seg grids differ.'); end
    mgxCb = mean(mgx(ismember(seg,[8 47])));

    % ===== match ROIs across all three =====
    rows = struct('segId',{},'name',{},'class',{},'ours',{},'mgx',{},'gtm',{});
    for i = 1:numel(G.segId)
        id = G.segId(i);
        mOur = ourSeg(id); mMgx = (seg==id);
        if nnz(mOur) < 30 || nnz(mMgx) < 30, continue; end
        rows(end+1) = struct('segId',id,'name',G.name{i},'class',G.class{i}, ...
            'ours', mean(pvc(mOur))/ourCb, 'mgx', mean(mgx(mMgx))/mgxCb, 'gtm', G.mean(i)/gtmCb); %#ok<AGROW>
    end
    R.rows = rows; R.subject = SubjectName; R.tracer = tracer;
    O=[rows.ours]'; M=[rows.mgx]'; T=[rows.gtm]'; cls={rows.class}; isctx=strcmp(cls,'cortex');

    % ===== pairwise metrics (all ROIs + cortex-only) =====
    fprintf('\n=== %s / %s : 3-way PVC regional SUVR (n=%d ROIs, %d cortical) ===\n', SubjectName, tracer, numel(O), nnz(isctx));
    fprintf('%-28s  %7s %7s | %8s %7s\n','pair (method vs method)','r_all','CCC_all','r_ctx','CCC_ctx');
    prs = {'OURS vs MGX  (xpipe, same MG)','ours','mgx'; ...
           'MGX  vs GTM  (PETSurfer, x-method)','mgx','gtm'; ...
           'OURS vs GTM  (xpipe, x-method)','ours','gtm'};
    V = struct('ours',O,'mgx',M,'gtm',T);
    R.metrics = struct();
    for k=1:size(prs,1)
        a=V.(prs{k,2}); b=V.(prs{k,3});
        fprintf('%-28s  %7.3f %7.3f | %8.3f %7.3f\n', prs{k,1}, ...
            local_corr(a,b), local_ccc(a,b), local_corr(a(isctx),b(isctx)), local_ccc(a(isctx),b(isctx)));
    end
    fprintf('[weakness] OURS corrects GM only (WM/CSF set to constants) -> non-cortex ROIs diverge by design.\n');

    % ===== figure: 3 scatter panels =====
    local_fig3(O,M,T,isctx,SubjectName,tracer,OutDir);
    fprintf('Figure saved to %s\n', OutDir);
end


%% ===== helpers =====
function C = local_cube(sS, comment)
    C = [];
    i = find(strcmp({sS.Anatomy.Comment}, comment), 1);
    if isempty(i), return; end
    s = in_mri_bst(sS.Anatomy(i).FileName); C = double(s.Cube(:,:,:,1));
end

function m = local_seg(id, aseg, apar, hasApar)
    if id >= 1000 && hasApar, m = (apar == id); else, m = (aseg == id); end
end

function vol = local_niigz(f)
% Non-interactive read of a .nii.gz (gunzip to temp, niftiread, cleanup). No Brainstorm importer.
    T = bst_get('BrainstormTmpDir', 0, 'niigz');
    [~,b] = bst_fileparts(f);                      % strip .gz -> name.nii
    copyfile(f, fullfile(T,[b '.gz'])); gunzip(fullfile(T,[b '.gz']));
    vol = double(niftiread(fullfile(T,b)));
    file_delete(T,1,1);
end

function G = local_read_gtm(statFile)
    fid=fopen(statFile,'r'); C=textscan(fid,'%f %f %s %s %f %f %f %f','CommentStyle','#'); fclose(fid);
    G.segId=C{2}; G.name=C{3}; G.class=C{4}; G.nvox=C{5}; G.mean=C{7};
end

function c = local_corr(x,y), cc=corrcoef(x,y); c=cc(1,2); end
function c = local_ccc(x,y)
    mx=mean(x); my=mean(y); vx=var(x,1); vy=var(y,1); sxy=mean((x-mx).*(y-my));
    c = 2*sxy/(vx+vy+(mx-my)^2);
end

function local_fig3(O,M,T,isctx,subj,tracer,OutDir)
    lbX={'MGX (PETSurfer MG)','GTM (PETSurfer)','GTM (PETSurfer)'};
    lbY={'OURS (PETPVE12 MG)','MGX (PETSurfer MG)','OURS (PETPVE12 MG)'};
    ttl={'cross-pipeline, SAME method (MG)','same pipeline, cross-method','cross-pipeline, cross-method'};
    pairs={M,O; T,M; T,O};
    f=figure('Visible','off','Position',[80 80 1400 440]);
    for k=1:3
        subplot(1,3,k); hold on; x=pairs{k,1}; y=pairs{k,2};
        scatter(x(~isctx),y(~isctx),20,[.7 .4 .2],'filled','MarkerFaceAlpha',.5);
        scatter(x(isctx), y(isctx), 24,[.4 .2 .6],'filled');
        lim=[0 max([x;y])*1.05]; plot(lim,lim,'k--'); xlim(lim); ylim(lim); axis square;
        xlabel(lbX{k}); ylabel(lbY{k});
        title(sprintf('%s\nr_{ctx}=%.2f CCC_{ctx}=%.2f', ttl{k}, local_corr(x(isctx),y(isctx)), local_ccc(x(isctx),y(isctx))));
        grid on; if k==1, legend({'non-cortex','cortex'},'Location','northwest'); end
    end
    sgtitle(sprintf('%s / %s  -  3-way PVC regional SUVR', subj, tracer),'Interpreter','none');
    print(f, fullfile(OutDir, sprintf('compare_pvc_3way_%s_%s.png', subj, tracer)),'-dpng','-r110'); close(f);
end
