function files = sa_figures(S, P, C, outDir)
% SA_FIGURES: Render the sign-ambiguity continuity figures.
%
% USAGE:  files = sa_figures(S, P, C, outDir)
%   S      : Component 1 stats (from sa_smoothness); [] to skip.
%   P      : crossing profiles (from sa_crossing_profile); [] to skip.
%   C      : Component 2 stats (from sa_continuity);  [] to skip.
%   outDir : output directory (created if absent).
% OUTPUT files: cell array of written PNG paths.

if ~exist(outDir, 'dir'), mkdir(outDir); end
files = {};

% ===== Component 1: decoupling bar chart + df/dn distributions =====
if ~isempty(S)
    h1 = figure('Visible','off','Color','w','Position',[100 100 900 380]);
    subplot(1,2,1);
    bar([S.dn_sulci_median S.dn_crown_median; S.df_sulci_median S.df_crown_median]);
    set(gca,'XTickLabel',{'normal \delta n','Fiedler \delta f'});
    legend({'sulcal edges','crown edges'},'Location','northeast');
    ylabel('median angular variation (rad)');
    title(sprintf('Folding coupling (sing. energy frac = %.2f)', S.singEnergyFrac));
    subplot(1,2,2);
    edges = linspace(0, pi, 30);
    d1 = S.df(S.defined & S.edgeSulcal);
    d2 = S.df(S.defined & ~S.edgeSulcal);
    if ~isempty(d1), histogram(d1, edges, 'Normalization','probability'); end
    hold on;
    if ~isempty(d2), histogram(d2, edges, 'Normalization','probability'); end
    legend({'sulcal','crown'}); xlabel('\delta f (rad)'); ylabel('prob');
    title('Fiedler covariant variation by partition');
    f1 = fullfile(outDir, 'c1_smoothness.png');
    print(h1, '-dpng', f1); close(h1); files{end+1} = f1;
end

% ===== Component 1: sulcal-crossing profiles =====
if ~isempty(P)
    nP = numel(P);
    hP = figure('Visible','off','Color','w','Position',[100 100 320*nP 320]);
    for k = 1:nP
        subplot(1, nP, k);
        x = 1000 * P(k).arclen;                 % mm
        plot(x, P(k).nAngle, '-o', 'DisplayName','normal turn'); hold on;
        plot(x, P(k).fAngle, '-s', 'DisplayName','Fiedler turn');
        if ~isempty(P(k).s)
            sn = sign(P(k).s);                       % constrained sign trace
            plot(x, (pi/2)*(1+sn), ':', 'DisplayName','constrained sign');
        end
        xlabel('arc length (mm)'); ylabel('cumulative turn (rad)');
        title(sprintf('crossing %d', k));
        if k == 1, legend('Location','northwest'); end
    end
    fP = fullfile(outDir, 'c1_profiles.png');
    print(hP, '-dpng', fP); close(hP); files{end+1} = fP;
end

% ===== Component 2: per-frame phase-discontinuity histogram + summary =====
if ~isempty(C)
    h2 = figure('Visible','off','Color','w','Position',[100 100 900 380]);
    subplot(1,2,1);
    edges = linspace(0, pi, 30);
    histogram(C.df(~isnan(C.df)), edges, 'Normalization','probability'); hold on;
    histogram(C.dt(~isnan(C.dt)), edges, 'Normalization','probability');
    histogram(C.dx(~isnan(C.dx)), edges, 'Normalization','probability');
    legend({'Fiedler','tess\_tangents','global-xyz'});
    xlabel('phase discontinuity across wall (rad)'); ylabel('prob');
    title(sprintf('Phase continuity, %d pairs, t=%.0f ms', C.nPairs, 1000*C.tPeak));
    subplot(1,2,2);
    bar([C.signFlipRate; C.phaseDisc.fiedler_median; C.phaseDisc.tang_median; C.phaseDisc.xyz_median]);
    set(gca,'XTickLabel',{'constrained sign-flip','Fiedler \theta','tang \theta','xyz \theta'});
    xtickangle(20); ylabel('rate / median rad');
    title('Constrained artifact vs. frame phase discontinuity');
    f2 = fullfile(outDir, 'c2_continuity.png');
    print(h2, '-dpng', f2); close(h2); files{end+1} = f2;
end
end
