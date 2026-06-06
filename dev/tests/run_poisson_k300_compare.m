% RUN_POISSON_K300_COMPARE  Compare K=20 vs K=300 analytic inverse + Poisson sharpening.
% Saves results to dev/tests/poisson_k300_results.mat

if ~brainstorm('status'), brainstorm nogui; end
outFile = fullfile(fileparts(mfilename('fullpath')), 'poisson_k300_results.mat');

HM_K20  = 'Subject01/S01_AEF_20131218_01_notch/headmodel_face_eigenmode_260606_0516.mat';
HM_K300 = 'Subject01/S01_AEF_20131218_01_notch/headmodel_face_eigenmode_260606_0534.mat';
DATA    = 'Subject01/S01_AEF_20131218_02_notch/data_block001_band.mat';

%% Inverses
fprintf('[1/4] K=20 analytic inverse...\n'); tic
R20  = bst_eigenmode_analytic_inverse(HM_K20,  DATA, 'SNR', 3);
t20 = toc; fprintf('  %.2fs\n', t20);

fprintf('[2/4] K=300 analytic inverse...\n'); tic
R300 = bst_eigenmode_analytic_inverse(HM_K300, DATA, 'SNR', 3);
t300 = toc; fprintf('  %.2fs\n', t300);

%% Poisson sharpening
fprintf('[3/4] Poisson K=20...\n');  P20  = bst_eigenmode_poisson_sharpen(R20);
fprintf('[4/4] Poisson K=300...\n'); P300 = bst_eigenmode_poisson_sharpen(R300);

%% Surface + adjacency
TessMat = in_tess_bst('Subject01/tess_cortex_pial_low.mat');
nV = size(TessMat.Vertices,1);
[~, lH_v] = tess_hemisplit(TessMat); lH_v = lH_v(:);
Vtx_mm = TessMat.Vertices * 1000;
FA = sparse(double(TessMat.Faces(:,[1 2 3])), double(TessMat.Faces(:,[2 3 1])), true, nV, nV);
FA = FA | FA';

%% Metrics for each of the 4 source fields
labels = {'K=20  inv', 'K=20  Poisson', 'K=300 inv', 'K=300 Poisson'};
fields = {R20, P20, R300, P300};
F_img  = {'ImageGridAmp', 'RhoGridAmp', 'ImageGridAmp', 'RhoGridAmp'};

results = struct();
for fi = 1:4
    S     = fields{fi};
    fld   = F_img{fi};
    Img   = S.(fld);

    Fs    = 1/mean(diff(S.Time));
    mid   = round(numel(S.Time)/2);
    win   = max(1,mid-round(Fs)) : min(numel(S.Time),mid+round(Fs));
    Amp   = abs(Img(lH_v, win));
    [~,tPk] = max(mean(Amp,1));
    amp_pk  = mean(Amp(:, max(1,tPk-150):min(numel(win),tPk+150)), 2);
    thr   = 0.08 * max(amp_pk);
    act   = lH_v(amp_pk > thr);

    % Smoothness
    [ea,eb] = find(triu(FA(act,act),1));
    va = act(ea); vb = act(eb);
    lhmap = zeros(nV,1); lhmap(lH_v) = 1:numel(lH_v);
    Ph = angle(Img(lH_v, win));
    dphi = abs(angle(exp(1i*(Ph(lhmap(va),tPk) - Ph(lhmap(vb),tPk)))));
    smoothness = mean(dphi);

    % Focality: active count + peak amplitude
    nActive   = numel(act);
    peakAmp   = max(amp_pk);

    % Phase lead between two peak vertices ≥30mm apart
    [~,sord] = sort(amp_pk,'descend');
    iV1 = NaN; iV2 = NaN;
    for k1 = 1:numel(sord)
        for k2 = k1+1 : min(k1+100, numel(sord))
            d = norm(Vtx_mm(act(sord(k1)),:) - Vtx_mm(act(sord(k2)),:));
            if d > 30
                iV1 = act(sord(k1)); iV2 = act(sord(k2)); break
            end
        end
        if ~isnan(iV1), break; end
    end
    if ~isnan(iV1)
        iWin = max(1,tPk-round(Fs)):min(numel(win),tPk+round(Fs));
        ph_lead = mean(angle(exp(1i*(angle(Img(iV1,win(iWin))) - angle(Img(iV2,win(iWin)))))));
        dist_mm = norm(Vtx_mm(iV1,:) - Vtx_mm(iV2,:));
        delay_ms = ph_lead / (2*pi*10) * 1000;
        speed_ms = dist_mm / max(abs(delay_ms), 0.1);
        plv = abs(mean(exp(1i*(angle(Img(iV1,win)) - angle(Img(iV2,win))))));
    else
        ph_lead = NaN; dist_mm = NaN; delay_ms = NaN; speed_ms = NaN; plv = NaN;
    end

    results(fi).label      = labels{fi};
    results(fi).nModes     = S.nModes;
    results(fi).smoothness = smoothness;
    results(fi).nActive    = nActive;
    results(fi).peakAmp    = peakAmp;
    results(fi).phaseLead  = ph_lead;
    results(fi).delay_ms   = delay_ms;
    results(fi).speed_ms   = speed_ms;
    results(fi).plv        = plv;
    results(fi).alpha      = getfield_or(S, 'Alpha', NaN);

    fprintf('\n%s\n', labels{fi});
    fprintf('  smoothness  = %.4f rad/edge\n', smoothness);
    fprintf('  active verts = %d / %d\n', nActive, numel(lH_v));
    fprintf('  peak amp    = %.3e\n', peakAmp);
    fprintf('  phase lead  = %.3f rad = %.1f ms  speed=%.2f m/s\n', ph_lead, delay_ms, speed_ms);
    fprintf('  PLV         = %.3f\n', plv);
end

%% Summary table
fprintf('\n%s\n', repmat('=',1,70));
fprintf('%-20s  %8s  %7s  %7s  %7s  %5s\n', ...
    'Method','smooth(r/e)','nActive','phaseLd','spd(m/s)','PLV');
fprintf('%s\n', repmat('-',1,70));
for fi = 1:4
    r = results(fi);
    fprintf('%-20s  %8.4f  %7d  %7.3f  %7.2f  %5.3f\n', ...
        r.label, r.smoothness, r.nActive, r.phaseLead, r.speed_ms, r.plv);
end
fprintf('%s\n', repmat('=',1,70));

save(outFile, 'results', 't20', 't300');
fprintf('\nSaved -> %s\n', outFile);

function v = getfield_or(s, f, def)
    if isfield(s,f), v = s.(f); else, v = def; end
end
