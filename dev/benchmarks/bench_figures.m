function files = bench_figures(csvPath, figDir, renderCtx)
% BENCH_FIGURES: Produce the five benchmark figures (PNG + .fig).
% USAGE: files = bench_figures(csvPath, figDir, renderCtx)
%   renderCtx : [] to skip the cortex-render figure, or struct with fields
%               .Vertices .Faces .gt .estMaps (struct method->[nV x 1]) .titleStr
if nargin < 3; renderCtx = []; end
if ~exist(figDir,'dir'); mkdir(figDir); end
T = readtable(csvPath);
T.method = string(T.method); T.regime = string(T.regime);
files = {};

eigMethods = ["eig_mne_log","eig_dspm_log"];
stdMethods = ["wmne","dspm","sloreta"];
allMethods = [stdMethods, eigMethods];
Kmax = max(T.K(~isnan(T.K)));
col = lines(numel(allMethods));

    function v = locOf(m, regime, snr)
        sel = T.method==m;
        if ismember(m, eigMethods); sel = sel & T.K==Kmax; end
        if nargin>=2 && ~isempty(regime); sel = sel & T.regime==regime; end
        if nargin>=3 && ~isempty(snr);    sel = sel & T.snr_db==snr; end
        v = T.locerror_mm(sel);
    end

% ----- F1: Distribution (strip + median/IQR) of LocError per method -----
f1 = figure('Color','w','Position',[100 100 760 460]); hold on;
for im=1:numel(allMethods)
    v = locOf(allMethods(im));
    x = im + (rand(numel(v),1)-0.5)*0.3;
    scatter(x, v, 14, col(im,:), 'filled', 'MarkerFaceAlpha',0.45);
    med = median(v); q1 = quantile_local(v,0.25); q3 = quantile_local(v,0.75);
    plot([im-0.32 im+0.32],[med med],'k-','LineWidth',2);
    plot([im im],[q1 q3],'k-','LineWidth',1);
end
set(gca,'XTick',1:numel(allMethods),'XTickLabel',cellstr(allMethods),'XTickLabelRotation',20);
ylabel('Localization error (mm)'); title('LocError distribution by method (eig at K_{max})'); grid on;
files{end+1} = save_fig(f1, figDir, 'f1_distribution');

% ----- F2: SNR sweep, one panel per regime -----
regimes = unique(T.regime,'stable'); snrs = unique(T.snr_db)';
f2 = figure('Color','w','Position',[100 100 1100 360]);
for ir=1:numel(regimes)
    subplot(1,numel(regimes),ir); hold on;
    for im=1:numel(allMethods)
        mu = arrayfun(@(s) mean(locOf(allMethods(im), regimes(ir), s)), snrs);
        sd = arrayfun(@(s) std(locOf(allMethods(im), regimes(ir), s)),  snrs);
        errorbar(snrs, mu, sd, '-o', 'Color', col(im,:), 'MarkerFaceColor', col(im,:), 'CapSize',4);
    end
    xlabel('SNR (dB)'); ylabel('LocError (mm)'); title(regimes(ir)); grid on;
    if ir==numel(regimes); legend(cellstr(allMethods),'Location','northeastoutside'); end
end
sgtitle('LocError vs SNR by regime');
files{end+1} = save_fig(f2, figDir, 'f2_snr_sweep');

% ----- F4: Per-regime grouped bars (median LocError, IQR whiskers) -----
f4 = figure('Color','w','Position',[100 100 900 460]); hold on;
nM = numel(allMethods); nR = numel(regimes); bw = 0.8/nM;
for im=1:nM
    meds = arrayfun(@(ir) median(locOf(allMethods(im), regimes(ir))), 1:nR);
    q1   = arrayfun(@(ir) quantile_local(locOf(allMethods(im), regimes(ir)),0.25), 1:nR);
    q3   = arrayfun(@(ir) quantile_local(locOf(allMethods(im), regimes(ir)),0.75), 1:nR);
    xpos = (1:nR) + (im-(nM+1)/2)*bw;
    bar(xpos, meds, bw*0.9, 'FaceColor', col(im,:), 'EdgeColor','none');
    errorbar(xpos, meds, meds-q1, q3-meds, 'k', 'LineStyle','none', 'CapSize',3);
end
set(gca,'XTick',1:nR,'XTickLabel',cellstr(regimes)); ylabel('Median LocError (mm)');
legend(cellstr(allMethods),'Location','northeastoutside'); title('LocError by regime x method'); grid on;
files{end+1} = save_fig(f4, figDir, 'f4_per_regime');

% ----- F5: K-sweep curve (focal emphasis), eig methods vs std reference lines -----
Kvals = unique(T.K(~isnan(T.K)))';
f5 = figure('Color','w','Position',[100 100 1100 360]);
for ir=1:numel(regimes)
    subplot(1,numel(regimes),ir); hold on;
    for ie=1:numel(eigMethods)
        med = arrayfun(@(k) median(T.locerror_mm(T.method==eigMethods(ie) & T.regime==regimes(ir) & T.K==k)), Kvals);
        plot(Kvals, med, '-o', 'LineWidth',1.5, 'DisplayName', char(eigMethods(ie)));
    end
    for is=1:numel(stdMethods)
        ref = median(locOf(stdMethods(is), regimes(ir)));
        yline(ref, '--', char(stdMethods(is)), 'LabelHorizontalAlignment','left', 'HandleVisibility','off');
    end
    xlabel('Total modes K'); ylabel('Median LocError (mm)'); title(regimes(ir)); grid on;
    if ir==1; legend('Location','northeast'); end
end
sgtitle('K-sweep: eigenmode LocError vs total modes (std = dashed reference)');
files{end+1} = save_fig(f5, figDir, 'f5_ksweep');

% ----- F3: Example reconstructions on cortex (only if renderCtx provided) -----
if ~isempty(renderCtx)
    methodsR = fieldnames(renderCtx.estMaps);
    f3 = figure('Color','w','Position',[100 100 1200 420]);
    npan = 1 + numel(methodsR);
    subplot(1,npan,1);
    render_map(renderCtx.Vertices, renderCtx.Faces, renderCtx.gt); title('ground truth');
    for i=1:numel(methodsR)
        subplot(1,npan,1+i);
        render_map(renderCtx.Vertices, renderCtx.Faces, renderCtx.estMaps.(methodsR{i}));
        title(strrep(methodsR{i},'_','\_'));
    end
    sgtitle(renderCtx.titleStr);
    files{end+1} = save_fig(f3, figDir, 'f3_cortex_reconstruction');
end
end

function p = save_fig(h, figDir, name)
p = fullfile(figDir, [name '.png']);
print(h, p, '-dpng', '-r150');
savefig(h, fullfile(figDir, [name '.fig']));
close(h);
end

function render_map(V, F, map)
map = abs(map(:)); if max(map)>0; map = map/max(map); end
patch('Vertices',V,'Faces',F,'FaceVertexCData',map,'FaceColor','interp', ...
      'EdgeColor','none','FaceLighting','gouraud');
axis equal off; view(0,90); camlight headlight; colormap(hot); caxis([0 1]);
end

function q = quantile_local(v, p)
v = sort(v(:)); n = numel(v);
if n==1; q=v; return; end
h = (n-1)*p + 1; lo = floor(h); hi = min(lo+1,n);
q = v(lo) + (h-lo)*(v(hi)-v(lo));
end
