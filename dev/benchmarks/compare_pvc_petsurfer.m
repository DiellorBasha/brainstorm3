function R = compare_pvc_petsurfer(SubjectName, tracer, PvcComment, PsRoot)
% COMPARE_PVC_PETSURFER: Cross-pipeline check of our PETPVE12 Muller-Gartner PVC
% against a FreeSurfer PETSurfer GTM run, as REGIONAL SUVR over FreeSurfer ROIs.
%
% Both pipelines are put on a common footing by normalizing each to its OWN
% cerebellum-cortex mean (FreeSurfer aseg 8 + 47), so only the regional PATTERN is
% compared (scale conventions cancel). Our regional means come from our PVC volume
% sampled over the subject's ASEG (subcortical/WM/CSF) + Desikan-Killiany aparc
% (cortex) volume atlases; PETSurfer's come from gtm.stats.dat (Rousset GTM).
%
% NOTE: this compares our voxelwise MG (regional mean) vs PETSurfer's regional GTM
% (Rousset) -- different PVC algorithms, so expect strong correlation, not identity.
%
% USAGE:  R = compare_pvc_petsurfer('sub-MTL0002','18FNAV4694','PET 18FNAV4694_mean_pvc')
%
% Author: Diellor Basha, 2026

    if (nargin < 4) || isempty(PsRoot)
        PsRoot = '/Volumes/SpikeData-2/workspace/library/datasets/preventad/pet/derivatives/petsurfer';
    end
    OutDir = bst_fileparts(mfilename('fullpath'));

    % ===== our pipeline: PVC volume + segmentations (all T1-grid) =====
    [sS,~] = bst_get('Subject', SubjectName);
    pvc  = local_cube(sS, PvcComment, 1);
    aseg = local_cube(sS, 'ASEG', 0);
    apar = local_cube(sS, 'Desikan-Killiany', 0);     % aparc+aseg (cortex 1000+)
    if isempty(pvc),  error('PVC volume "%s" not found.', PvcComment); end
    if isempty(aseg), error('ASEG not found.'); end
    if ~isequal(size(pvc), size(aseg)), error('PVC/ASEG grids differ.'); end
    hasApar = ~isempty(apar) && isequal(size(apar), size(pvc));

    segLabel = @(id) local_seg(id, aseg, apar, hasApar);   % voxel mask for a segId

    % our cerebellum-cortex reference (aseg 8,47), then per-ROI means
    cbMask = segLabel(8) | segLabel(47);
    ourCb  = mean(pvc(cbMask));

    % ===== PETSurfer GTM regional table =====
    statFile = fullfile(PsRoot, SubjectName, 'ses-01', tracer, 'gtmpvc.output', 'gtm.stats.dat');
    if ~exist(statFile,'file'), error('PETSurfer stats not found: %s', statFile); end
    G = local_read_gtm(statFile);                          % .segId .name .class .nvox .mean
    iCb = ismember(G.segId, [8 47]);
    psCb = sum(G.mean(iCb).*G.nvox(iCb)) / sum(G.nvox(iCb));

    % ===== match ROIs, build SUVR vectors =====
    rows = struct('segId',{},'name',{},'class',{},'ours',{},'ps',{},'nvox',{});
    for i = 1:numel(G.segId)
        id = G.segId(i);
        m = segLabel(id);
        nv = nnz(m);
        if nv < 30, continue; end                          % skip tiny/absent ROIs
        rows(end+1) = struct('segId',id,'name',G.name{i},'class',G.class{i}, ...
            'ours', mean(pvc(m))/ourCb, 'ps', G.mean(i)/psCb, 'nvox', nv); %#ok<AGROW>
    end
    if isempty(rows), error('No overlapping ROIs found (check grids/atlases).'); end
    R.rows = rows; R.subject = SubjectName; R.tracer = tracer;

    ours = [rows.ours]'; ps = [rows.ps]'; cls = {rows.class};
    % ===== metrics =====
    rr  = local_corr(ours, ps);
    ccc = local_ccc(ours, ps);
    d   = ours - ps; biasm = mean(d); loa = 1.96*std(d);
    R.r = rr; R.ccc = ccc; R.bias = biasm; R.loa = loa; R.n = numel(ours);
    fprintf('\n=== %s / %s : our MG vs PETSurfer GTM (regional SUVR, n=%d ROIs) ===\n', SubjectName, tracer, numel(ours));
    fprintf('Pearson r = %.3f | Lin CCC = %.3f | bias(ours-ps) = %+.3f | 95%% LoA = +/-%.3f\n', rr, ccc, biasm, loa);
    fprintf('cerebellum ref: ours=%.3f  ps=%.3f (raw means)\n', ourCb, psCb);
    % per cortical-only (the pipeline target)
    isctx = strcmp(cls,'cortex');
    if any(isctx)
        fprintf('cortex-only (n=%d): r = %.3f | CCC = %.3f\n', nnz(isctx), local_corr(ours(isctx),ps(isctx)), local_ccc(ours(isctx),ps(isctx)));
    end

    % ===== figures =====
    local_fig(ours, ps, cls, rr, ccc, SubjectName, tracer, OutDir);
    fprintf('Figure saved to %s\n', OutDir);
end


%% ===== helpers =====
function C = local_cube(sS, comment, asDouble)
    C = [];
    i = find(strcmp({sS.Anatomy.Comment}, comment), 1);
    if isempty(i), return; end
    s = in_mri_bst(sS.Anatomy(i).FileName);
    C = s.Cube(:,:,:,1);
    if asDouble, C = double(C); end
end

function m = local_seg(id, aseg, apar, hasApar)
    if id >= 1000 && hasApar
        m = (apar == id);
    else
        m = (aseg == id);
    end
end

function G = local_read_gtm(statFile)
    fid = fopen(statFile,'r');
    C = textscan(fid, '%f %f %s %s %f %f %f %f', 'CommentStyle','#');
    fclose(fid);
    G.segId = C{2}; G.name = C{3}; G.class = C{4}; G.nvox = C{5}; G.mean = C{7};
end

function c = local_corr(x, y)
    cc = corrcoef(x, y); c = cc(1,2);
end

function c = local_ccc(x, y)
    mx=mean(x); my=mean(y); vx=var(x,1); vy=var(y,1); sxy=mean((x-mx).*(y-my));
    c = 2*sxy / (vx + vy + (mx-my)^2);
end

function local_fig(ours, ps, cls, rr, ccc, subj, tracer, OutDir)
    f = figure('Visible','off','Position',[100 100 1000 460]);
    % scatter ours vs ps
    subplot(1,2,1); hold on;
    uc = unique(cls,'stable'); col = lines(numel(uc));
    for k=1:numel(uc)
        m = strcmp(cls,uc{k});
        scatter(ps(m), ours(m), 24, col(k,:), 'filled');
    end
    lim = [0, max([ours;ps])*1.05]; plot(lim,lim,'k--');
    xlabel('PETSurfer GTM SUVR'); ylabel('our MG SUVR'); axis equal; xlim(lim); ylim(lim);
    legend([uc,{'identity'}],'Location','northwest','Interpreter','none');
    title(sprintf('%s / %s  (r=%.3f, CCC=%.3f)', subj, tracer, rr, ccc),'Interpreter','none'); grid on;
    % Bland-Altman
    subplot(1,2,2); hold on;
    mn=(ours+ps)/2; d=ours-ps;
    scatter(mn,d,24,[.3 .3 .3],'filled');
    yline(mean(d),'b-'); yline(mean(d)+1.96*std(d),'r--'); yline(mean(d)-1.96*std(d),'r--');
    xlabel('mean SUVR'); ylabel('ours - PETSurfer'); title('Bland-Altman'); grid on;
    print(f, fullfile(OutDir, sprintf('compare_pvc_petsurfer_%s_%s.png', subj, tracer)), '-dpng','-r120');
    close(f);
end
