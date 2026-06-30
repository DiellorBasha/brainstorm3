function R = bench_pet_surface_recovery(SubjectName, Opts)
% BENCH_PET_SURFACE_RECOVERY: Synthetic ground-truth benchmark for PET volume->surface sampling.
%
% Paints a KNOWN value field on the cortical surface, lifts it into a synthetic
% volume at realistic GM/WM/CSF contrast, forward-blurs with a scanner PSF (+noise),
% then RECOVERS the surface field with several depth-sampling profiles and measures
% how faithfully each one recovers the truth. This is the measurement instrument that
% gates the PET surface pipeline (sampling weights, and later PVC).
%
% It directly quantifies the ProjFrac sampling bug: the "as-shipped SUVR" profile
% (white-skewed) vs the bug-fixed/symmetric profiles.
%
% Everything is synthetic + anatomy-only (no PET data, no DB writes). Reuses the
% same surface<->voxel interpolation primitive as mri_interp_vol2tess (tess_tri_interp).
%
% USAGE:  R = bench_pet_surface_recovery('sub-MTL0002')
%         R = bench_pet_surface_recovery('sub-MTL0002', Opts)
%
% INPUTS:
%   SubjectName : subject with cortex_{white,mid,pial}_low surfaces + ASEG (default 'sub-MTL0002')
%   Opts (optional):
%     .FWHM        scanner PSF FWHM in mm (default 2.7, HRRT-like)
%     .NoisePct    Gaussian noise as % of GM level (default 0; set e.g. 5 for variance)
%     .OutDir      where to save PNGs (default this file's folder)
%     .DoFigures   1 to save PNG figures (default 1)
%
% OUTPUT:
%   R : struct with R.metrics (table-like struct array) and R.params
%
% Author: Diellor Basha, 2026

    if (nargin < 1) || isempty(SubjectName), SubjectName = 'sub-MTL0002'; end
    if (nargin < 2) || isempty(Opts), Opts = struct(); end
    Def = struct('FWHM',2.7, 'NoisePct',0, 'OutDir',bst_fileparts(mfilename('fullpath')), 'DoFigures',1);
    fn = fieldnames(Def);
    for i = 1:numel(fn)
        if ~isfield(Opts,fn{i}) || isempty(Opts.(fn{i})), Opts.(fn{i}) = Def.(fn{i}); end
    end

    % ===== ANATOMY =====
    [sSubject, iSubject] = bst_get('Subject', SubjectName); %#ok<ASGLU>
    if isempty(sSubject), error('Subject %s not found.', SubjectName); end
    T1file   = local_find_anat(sSubject, {'mri','t1'}, {'aseg','dkt','desikan','destrieux','pet','rescaled'});
    ASEGfile = local_find_anat(sSubject, {'aseg'}, {});
    if isempty(T1file),   error('No T1 anatomy found for %s.', SubjectName); end
    if isempty(ASEGfile), error('No ASEG volume found for %s (needed for WM mask).', SubjectName); end
    sMri = in_mri_bst(T1file);
    sAseg = in_mri_bst(ASEGfile);
    cubeSize = size(sMri.Cube(:,:,:,1));
    if ~isequal(size(sAseg.Cube(:,:,:,1)), cubeSize)
        error('ASEG grid does not match T1 grid.');
    end
    voxsize = sMri.Voxsize(1:3);   % mm

    % Surfaces (shared topology: mid = (white+pial)/2)
    [Vw, Vp, Faces, nVert] = local_load_surfaces(sSubject);
    thick_mm = sqrt(sum((Vp - Vw).^2, 2)) * 1000;   % cortical thickness per vertex (mm)

    % ===== DEPTH -> INTERPOLATION MATRIX (cached) =====
    Hcache = containers.Map('KeyType','char','ValueType','any');
    Hget = @(f) local_Hdepth(f, Vw, Vp, Faces, sMri, cubeSize, Hcache);

    % Valid vertices = covered by white, mid and pial (exclude medial wall / outside)
    cov = @(f) full(sum(Hget(f),1))';   % [nVert x 1] coverage
    valid = (cov(0) > 0) & (cov(0.5) > 0) & (cov(1) > 0);

    % ===== TISSUE COMPARTMENTS (geometry, contrast-independent) =====
    nVox = prod(cubeSize);
    Dlift = 0:0.2:1;                         % depths used to fill the GM ribbon
    onesV = ones(nVert,1);
    ribbonW = zeros(nVox,1);
    for f = Dlift, ribbonW = ribbonW + Hget(f)*onesV; end
    ribbonMask = ribbonW > 0;                                 % GM ribbon voxels
    wmMask  = ismember(sAseg.Cube(:,:,:,1), [2 41]);          % cerebral WM (FreeSurfer aseg)
    wmMask  = wmMask(:) & ~ribbonMask;
    csfMask = local_dilate(ribbonMask, cubeSize, 3) & ~ribbonMask & ~wmMask;   % peri-cortical CSF shell

    % ===== EXPERIMENT GRID =====
    % Sampling profiles: name + depths (0=white..1=pial) + weights
    variants = {
        'default[.1 .8 .1]'   [0 0.5 1]              [0.1 0.8 0.1]
        'suvr_buggy(white.5)' [0 0.5 1]              [0.5 0.4 0.1]
        'suvr_fixed(pial.5)'  [0 0.5 1]              [0.1 0.4 0.5]
        'mid_only'            [0.5]                  [1]
        'sym_broad'           [0 0.5 1]              [0.25 0.5 0.25]
        'gauss7(.35-.65)'     (0.35:0.05:0.65)       local_gaussw(0.35:0.05:0.65, 0.5, 0.143)
    };
    % Tracer contrast presets: [GMbase WM CSF]; field is painted on GM
    contrasts = struct('name',{'amyloid','tau'}, 'WM',{1.0,1.6}, 'CSF',{0.1,0.1});
    % Field patterns on the surface (function of white-vertex coords, returns [nVert x 1])
    patterns = {'const','gradient','blobs'};

    R = struct(); R.params = Opts; R.params.SubjectName = SubjectName;
    R.params.nVert = nVert; R.params.nValid = nnz(valid);
    M = struct('contrast',{},'pattern',{},'variant',{},'PSF',{},'RC',{},'biasMean',{},'biasSD',{},'RMSE',{},'r',{});

    % store per-vertex bias for figures (amyloid/const/PSF=FWHM)
    figBias = struct();

    for ic = 1:numel(contrasts)
        C = contrasts(ic);
        for ip = 1:numel(patterns)
            Strue = local_field(patterns{ip}, Vw);              % [nVert x 1], cortical values ~1..2.5
            % --- lift truth to volume ---
            ribbonField = zeros(nVox,1);
            for f = Dlift, ribbonField = ribbonField + Hget(f)*Strue; end
            volTrue = zeros(cubeSize);
            volTrue(ribbonMask) = ribbonField(ribbonMask) ./ ribbonW(ribbonMask);
            volTrue(wmMask)  = C.WM;
            volTrue(csfMask) = C.CSF;
            % --- PSF sweep: 0 (sanity) and FWHM ---
            for PSF = [0 Opts.FWHM]
                volB = volTrue;
                if PSF > 0, volB = local_gauss3(volTrue, PSF, voxsize); end
                if Opts.NoisePct > 0 && PSF > 0
                    volB = volB + (Opts.NoisePct/100)*mean(Strue) * local_randn_like(volB);
                end
                volVec = volB(:);
                for iv = 1:size(variants,1)
                    depths = variants{iv,2}; w = variants{iv,3}; w = w/sum(w);
                    rec = zeros(nVert,1);
                    for id = 1:numel(depths)
                        H = Hget(depths(id));
                        s = (H' * volVec) ./ max(full(sum(H,1))',eps);
                        rec = rec + w(id)*s;
                    end
                    [rc,bm,bs,rmse,rr] = local_metrics(rec(valid), Strue(valid));
                    M(end+1) = struct('contrast',C.name,'pattern',patterns{ip},'variant',variants{iv,1}, ...
                        'PSF',PSF,'RC',rc,'biasMean',bm,'biasSD',bs,'RMSE',rmse,'r',rr); %#ok<AGROW>
                    if strcmp(C.name,'amyloid') && strcmp(patterns{ip},'const') && PSF==Opts.FWHM
                        figBias.(matlab.lang.makeValidName(variants{iv,1})) = rec - Strue;
                    end
                end
            end
        end
    end
    R.metrics = M;

    % ===== SANITY GATE: mid-only, PSF=0, const must recover (RC~1, |bias| small) =====
    s0 = M(strcmp({M.contrast},'amyloid') & strcmp({M.pattern},'const') & ...
            strcmp({M.variant},'mid_only') & [M.PSF]==0);
    fprintf('\n[sanity] mid_only / PSF=0 / const: RC=%.3f biasMean=%.3f (expect ~1.0, ~0)\n', s0.RC, s0.biasMean);
    if abs(s0.RC-1) > 0.05 || abs(s0.biasMean) > 0.05*mean(local_field('const',Vw))
        warning('Sanity gate FAILED: lift/sample round-trip is biased (RC=%.3f).', s0.RC);
    else
        fprintf('[sanity] PASS\n');
    end

    % ===== REPORT (amyloid & tau, const pattern, at PSF=FWHM) =====
    local_print_table(M, Opts.FWHM);

    % ===== PVC EXPERIMENT (a): Muller-Gartner recovery + FWHM sensitivity =====
    % True PSF is Opts.FWHM; we PVC with the CORRECT FWHM and with over-specified
    % values (4 mm, 6 mm = the old generic default) to show the over-correction risk
    % on sharp HRRT data. Const field, mid-centered sampling.
    StrueC = local_field('const', Vw);
    ribbonFieldC = zeros(nVox,1);
    for f = Dlift, ribbonFieldC = ribbonFieldC + Hget(f)*StrueC; end
    assumedFwhms = unique([Opts.FWHM, 4.0, 6.0]);
    P = struct('contrast',{},'method',{},'assumedFWHM',{},'RC',{},'biasMean',{},'RMSE',{});
    smpD = [0 0.5 1]; smpW = [0.1 0.8 0.1];                 % benchmark-best mid-centered sampler
    for ic = 1:numel(contrasts)
        C = contrasts(ic);
        volTrue = zeros(cubeSize);
        volTrue(ribbonMask) = ribbonFieldC(ribbonMask) ./ ribbonW(ribbonMask);
        volTrue(wmMask)  = C.WM;
        volTrue(csfMask) = C.CSF;
        volB = local_gauss3(volTrue, Opts.FWHM, voxsize);   % true-PSF blur
        recN = local_sample(volB(:), Hget, smpD, smpW, nVert);
        [rc,bm,~,rmse] = local_metrics(recN(valid), StrueC(valid));
        P(end+1) = struct('contrast',C.name,'method','NoPVC','assumedFWHM',NaN,'RC',rc,'biasMean',bm,'RMSE',rmse); %#ok<AGROW>
        for af = assumedFwhms
            volP = local_mg_pvc(volB, ribbonMask, wmMask, cubeSize, voxsize, af);
            recP = local_sample(volP(:), Hget, smpD, smpW, nVert);
            [rc,bm,~,rmse] = local_metrics(recP(valid), StrueC(valid));
            P(end+1) = struct('contrast',C.name,'method',sprintf('MG@%.1f',af),'assumedFWHM',af, ...
                'RC',rc,'biasMean',bm,'RMSE',rmse); %#ok<AGROW>
        end
    end
    R.pvc = P;
    local_print_pvc(P, Opts.FWHM);

    % ===== FIGURES =====
    if Opts.DoFigures
        local_fig_bars(M, Opts, variants);
        local_fig_depthbias(figBias, thick_mm(valid), valid, Opts, variants);
        local_fig_pvc(P, Opts);
        fprintf('Figures saved to %s\n', Opts.OutDir);
    end
end


%% ============================ helpers ============================
function MriFile = local_find_anat(sSubject, mustInc, mustExc)
    MriFile = '';
    for i = 1:numel(sSubject.Anatomy)
        cm = lower(sSubject.Anatomy(i).Comment);
        if all(cellfun(@(s) contains(cm,s), mustInc)) && ~any(cellfun(@(s) contains(cm,s), mustExc))
            MriFile = sSubject.Anatomy(i).FileName; return;
        end
    end
end

function [Vw, Vp, Faces, nVert] = local_load_surfaces(sSubject)
    f = {sSubject.Surface.FileName};
    iw = find(~cellfun('isempty', regexp(f,'cortex_white_low\.mat$','once')),1);
    ip = find(~cellfun('isempty', regexp(f,'cortex_pial_low\.mat$','once')),1);
    if isempty(iw)||isempty(ip), error('cortex_white_low / cortex_pial_low not found.'); end
    sw = in_tess_bst(f{iw}); sp = in_tess_bst(f{ip});
    Vw = sw.Vertices; Vp = sp.Vertices; Faces = sw.Faces; nVert = size(Vw,1);
    if size(Vp,1)~=nVert, error('white/pial vertex counts differ (no correspondence).'); end
end

function H = local_Hdepth(f, Vw, Vp, Faces, sMri, cubeSize, cache)
    key = sprintf('%.3f', f);
    if isKey(cache, key), H = cache(key); return; end
    Vscs = (1-f)*Vw + f*Vp;                          % interpolate white->pial (corresponded)
    Vvox = cs_convert(sMri, 'scs', 'voxel', Vscs);
    H = tess_tri_interp(Vvox, Faces, cubeSize, 0);   % [nVox x nVert] sparse
    cache(key) = H; %#ok<NASGU>
end

function w = local_gaussw(depths, mu, sigma)
    w = exp(-((depths-mu).^2)/(2*sigma^2)); w = w/sum(w);
end

function S = local_field(pattern, Vw)
    n = size(Vw,1);
    switch pattern
        case 'const'
            S = 2.0*ones(n,1);
        case 'gradient'                              % anterior->posterior 1.0..2.5 (SCS x = anterior)
            x = Vw(:,1); x = (x-min(x))/(max(x)-min(x)+eps);
            S = 1.0 + 1.5*x;
        case 'blobs'                                 % baseline 1.0 + focal bumps to ~2.5
            S = 1.0*ones(n,1);
            seeds = round(linspace(1, n, 8)); seeds = seeds(2:end-1);
            for k = seeds
                d2 = sum((Vw - Vw(k,:)).^2, 2);
                S = S + 1.5*exp(-d2/(2*(0.012^2)));  % ~12 mm sigma (SCS meters)
            end
        otherwise, error('unknown pattern %s', pattern);
    end
end

function v = local_randn_like(vol)
    % deterministic pseudo-noise (no rng() so resume-safe): hashed index -> [-1,1]
    n = numel(vol); idx = (1:n)';
    v = reshape(sin(idx*12.9898)*43758.5453, size(vol)); v = 2*(v-floor(v))-1;
end

function vol = local_gauss3(vol, fwhm_mm, voxsize_mm)
    sig = fwhm_mm/2.35482;
    for d = 1:3
        sv = sig/voxsize_mm(d); r = max(1,ceil(3*sv)); x = -r:r;
        k = exp(-(x.^2)/(2*sv^2)); k = k/sum(k);
        sh = ones(1,3); sh(d) = numel(k); k = reshape(k, sh);
        vol = convn(vol, k, 'same');
    end
end

function m = local_dilate(maskVec, cubeSize, r)
    mask = reshape(maskVec, cubeSize);
    k = ones(2*r+1, 2*r+1, 2*r+1);
    m = convn(double(mask), k, 'same') > 0.5;
    m = m(:);
end

function [rc,bm,bs,rmse,rr] = local_metrics(rec, tru)
    b = rec - tru;
    rc = mean(rec)/mean(tru);
    bm = mean(b); bs = std(b); rmse = sqrt(mean(b.^2));
    cc = corrcoef(rec, tru); rr = cc(1,2);
end

function rec = local_sample(volVec, Hget, depths, w, nVert)
    w = w/sum(w); rec = zeros(nVert,1);
    for id = 1:numel(depths)
        H = Hget(depths(id));
        rec = rec + w(id) * ((H'*volVec) ./ max(full(sum(H,1))', eps));
    end
end

function volP = local_mg_pvc(volB, ribbonMask, wmMask, cubeSize, voxsize, assumedFwhm)
% Muller-Gartner voxelwise PVC (CSF zeroed): correct GM for WM spill-in and GM
% spill-out using a Gaussian PSF at the ASSUMED FWHM. WM mean estimated from eroded WM.
    wmEr = local_erode(wmMask, cubeSize, 2);
    wmMean = mean(volB(wmEr));
    if ~isfinite(wmMean), wmMean = mean(volB(wmMask)); end
    wmImg = zeros(cubeSize); wmImg(wmMask) = wmMean;
    spillin = local_gauss3(wmImg, assumedFwhm, voxsize);
    gmFrac  = local_gauss3(reshape(double(ribbonMask), cubeSize), assumedFwhm, voxsize);
    volP = volB;
    volP(ribbonMask) = (volB(ribbonMask) - spillin(ribbonMask)) ./ max(gmFrac(ribbonMask), 0.1);
end

function m = local_erode(maskVec, cubeSize, r)
    m = ~local_dilate(~maskVec, cubeSize, r);
end

function local_print_pvc(P, trueFwhm)
    fprintf('\n=== PVC recovery (const field, true PSF=%.1f mm) ===\n', trueFwhm);
    for cn = {'amyloid','tau'}
        S = P(strcmp({P.contrast},cn{1}));
        fprintf('-- %-7s  %-10s %7s %9s %8s\n', cn{1}, 'method','RC','biasMean','RMSE');
        for i = 1:numel(S)
            fprintf('            %-10s %7.3f %9.3f %8.3f\n', S(i).method, S(i).RC, S(i).biasMean, S(i).RMSE);
        end
    end
end

function local_fig_pvc(P, Opts)
    methods = unique({P.method}, 'stable');
    amy = nan(numel(methods),1); tau = nan(numel(methods),1);
    for i = 1:numel(methods)
        a = P(strcmp({P.contrast},'amyloid') & strcmp({P.method},methods{i}));
        t = P(strcmp({P.contrast},'tau')     & strcmp({P.method},methods{i}));
        if ~isempty(a), amy(i)=a.RC; end
        if ~isempty(t), tau(i)=t.RC; end
    end
    f = figure('Visible','off','Position',[100 100 800 450]);
    bar([amy tau]); set(gca,'XTick',1:numel(methods),'XTickLabel',methods);
    hold on; yline(1,'k--','perfect recovery'); ylabel('Recovery coefficient'); ylim([0 1.35]);
    legend({'amyloid','tau'},'Location','northwest');
    title(sprintf('MG-PVC recovery vs assumed FWHM (true PSF=%.1f mm)', Opts.FWHM)); grid on;
    print(f, fullfile(Opts.OutDir,'bench_pet_pvc_recovery.png'), '-dpng','-r120'); close(f);
end

function local_print_table(M, FWHM)
    for cn = {'amyloid','tau'}
        sel = strcmp({M.contrast},cn{1}) & strcmp({M.pattern},'const') & [M.PSF]==FWHM;
        S = M(sel);
        fprintf('\n=== %s, const field, PSF=%.1f mm ===\n', cn{1}, FWHM);
        fprintf('%-22s %7s %9s %8s %7s\n','variant','RC','biasMean','RMSE','r');
        for i = 1:numel(S)
            fprintf('%-22s %7.3f %9.3f %8.3f %7.3f\n', S(i).variant, S(i).RC, S(i).biasMean, S(i).RMSE, S(i).r);
        end
    end
end

function local_fig_bars(M, Opts, variants)
    names = variants(:,1); nv = numel(names);
    rc_amy = zeros(nv,1); rc_tau = zeros(nv,1);
    for i = 1:nv
        a = M(strcmp({M.contrast},'amyloid') & strcmp({M.pattern},'const') & [M.PSF]==Opts.FWHM & strcmp({M.variant},names{i}));
        t = M(strcmp({M.contrast},'tau')     & strcmp({M.pattern},'const') & [M.PSF]==Opts.FWHM & strcmp({M.variant},names{i}));
        rc_amy(i) = a.RC; rc_tau(i) = t.RC;
    end
    f = figure('Visible','off','Position',[100 100 900 450]);
    bar([rc_amy rc_tau]); set(gca,'XTick',1:nv,'XTickLabel',names,'XTickLabelRotation',25);
    hold on; yline(1,'k--'); ylabel('Recovery coefficient (mean_{rec}/mean_{true})');
    legend({'amyloid (WM<GM)','tau (WM~GM)'},'Location','best');
    title(sprintf('PET surface recovery by sampling profile (const field, PSF=%.1f mm)', Opts.FWHM));
    grid on;
    print(f, fullfile(Opts.OutDir,'bench_pet_surface_recovery_rc.png'), '-dpng','-r120'); close(f);
end

function local_fig_depthbias(figBias, thickValid, valid, Opts, variants)
    show = {'suvr_buggy(white.5)','suvr_fixed(pial.5)','mid_only','gauss7(.35-.65)'};
    f = figure('Visible','off','Position',[100 100 900 450]); hold on;
    cols = lines(numel(show)); leg = {};
    for i = 1:numel(show)
        key = matlab.lang.makeValidName(show{i});
        if ~isfield(figBias,key), continue; end
        b = figBias.(key); b = b(valid);
        [ts,o] = sort(thickValid); bb = b(o);
        % bin by thickness for a clean trend (ts is sorted; use ~99th pct as upper edge)
        ts99 = ts(max(1, round(0.99*numel(ts))));
        edges = linspace(min(ts), ts99, 20); ctr = (edges(1:end-1)+edges(2:end))/2;
        mb = arrayfun(@(k) mean(bb(ts>=edges(k) & ts<edges(k+1))), 1:numel(edges)-1);
        plot(ctr, mb, '-o','Color',cols(i,:),'LineWidth',1.5); leg{end+1}=show{i}; %#ok<AGROW>
    end
    yline(0,'k--'); xlabel('cortical thickness (mm)'); ylabel('recovery bias (rec - true)');
    title('Recovery bias vs cortical thickness (amyloid, const, blurred)'); legend(leg,'Location','best'); grid on;
    print(f, fullfile(Opts.OutDir,'bench_pet_surface_recovery_depthbias.png'), '-dpng','-r120'); close(f);
end
