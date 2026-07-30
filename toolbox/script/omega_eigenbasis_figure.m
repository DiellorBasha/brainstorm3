function outPng = omega_eigenbasis_figure(sp)
% OMEGA_EIGENBASIS_FIGURE  Polished summary figure of the DIRAC-EIGENBASIS OMEGA cortical-flow
% dynamics: group alpha flow maps + standing regime screen + vortex-gas feature tracking.
%
%   outPng = omega_eigenbasis_figure(scratchpadDir)
%
% Loads the saved eigenbasis results (regime: omega_regime_eig.mat; tracking:
% omega_track_sub{2,3}.mat) and the group flow-map timefreq averages, and renders one labelled
% multi-panel figure (scale bars / colorbars / units). Author: Diellor Basha, 2026

    if nargin < 1 || isempty(sp)
        sp = '/private/tmp/claude-501/-Users-diellorbasha-workspace-research-code-brainstorm3/079e3de8-7e95-477e-8fd2-58b8dd7eb94d/scratchpad/';
    end
    gui_brainstorm('SetCurrentProtocol', bst_get('Protocol', 'omega_tutorial_test'));

    % ---- data ----
    Rg = load([sp 'omega_regime_eig.mat']);                       % rs2, rs3
    T2 = load([sp 'omega_track_sub2.mat'], 'tracks', 'out');
    T3 = load([sp 'omega_track_sub3.mat'], 'tracks', 'out');
    gpaths = {'Group_analysis/@intra/timefreq_psd_average_260729_2350_groupflow_div_dirac.mat', 'divergence'; ...
              'Group_analysis/@intra/timefreq_psd_average_260730_0007_groupflow_curl_dirac.mat', 'vorticity'; ...
              'Group_analysis/@intra/timefreq_psd_average_260730_0022_groupflow_psi_dirac.mat', 'stream \Psi'};
    Tg = in_bst_timefreq(gpaths{1,1}, 0);  Sg = in_tess_bst(Tg.SurfaceFile);
    iAlpha = 3;   % delta theta ALPHA beta gamma

    fig = figure('Color','w','Position',[40 40 1280 940]);
    tl = tiledlayout(fig, 3, 3, 'TileSpacing','compact','Padding','compact');

    % ---- row 1: group alpha flow maps on the cortex ----
    for i = 1:3
        T = in_bst_timefreq(gpaths{i,1}, 0);  m = T.TF(:,1,iAlpha);
        ax = nexttile(tl, i);
        patch(ax, 'Vertices', Sg.Vertices, 'Faces', Sg.Faces, 'FaceVertexCData', m, ...
            'FaceColor','interp','EdgeColor','none'); axis(ax,'equal','off');
        view(ax, [0 90]); camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
        cl = [prctile(m,2) prctile(m,99)];  set(ax,'CLim',cl);  colormap(ax, turbo(256));  % positive power -> sequential
        cb = colorbar(ax,'southoutside'); cb.Label.String = 'relative alpha power'; cb.Ticks = round(cl,3);
        title(ax, sprintf('group alpha %s', gpaths{i,2}), 'FontWeight','bold');
    end

    % ---- row 2: regime screen (standing) ----
    ax = nexttile(tl, 4);
    jp = Rg.rs2.jointPower;  f = Rg.rs2.f;  sl = sqrt(Rg.rs2.lam);
    fm = f <= 20;
    imagesc(ax, f(fm), sl(3:end), 10*log10(jp(3:end,fm)/max(jp(:)))); axis(ax,'xy');
    colormap(ax, turbo(256)); set(ax,'CLim',[-40 0]); c=colorbar(ax); c.Label.String='dB';
    xlabel(ax,'temporal frequency (Hz)'); ylabel(ax,'\surd\lambda (1/m)');
    title(ax,'joint spectrum |c(\lambda,\omega)|^2  (sub-0002)');

    ax = nexttile(tl, 5);
    plot(ax, Rg.rs2.speeds, Rg.rs2.E/mean(Rg.rs2.E), '-', 'LineWidth',1.6); hold(ax,'on');
    plot(ax, Rg.rs3.speeds, Rg.rs3.E/mean(Rg.rs3.E), '-', 'LineWidth',1.6);
    yline(ax, 1, 'k:'); xlabel(ax,'phase speed c (m/s)'); ylabel(ax,'E(c)/mean');
    legend(ax, {sprintf('sub-0002 (pk %.2f)',Rg.rs2.peakedness), sprintf('sub-0003 (pk %.2f)',Rg.rs3.peakedness)}, 'Location','northeast');
    title(ax,'traveling-wave speed sweep — peaks at c\rightarrow0 (STANDING)');

    % counts-by-type bar
    ax = nexttile(tl, 6);
    tys = {'vortex','source','sink','saddle'};
    c2 = cellfun(@(t) sum(strcmp({T2.tracks.type},t)), tys);
    c3 = cellfun(@(t) sum(strcmp({T3.tracks.type},t)), tys);
    bar(ax, [c2(:) c3(:)]); set(ax,'XTickLabel',tys); ylabel(ax,'# tracks');
    legend(ax,{'sub-0002','sub-0003'}); title(ax,'critical-point tracks by type');

    % ---- row 3: vortex-gas distributions ----
    ax = nexttile(tl, 7);
    histogram(ax, 1000*T2.out.lifetime, 0:20:600, 'Normalization','probability'); hold(ax,'on');
    histogram(ax, 1000*T3.out.lifetime, 0:20:600, 'Normalization','probability');
    xlabel(ax,'lifetime (ms)'); ylabel(ax,'prob'); title(ax,'track lifetime (\approx1 alpha cycle)');
    legend(ax,{'sub-0002','sub-0003'});

    ax = nexttile(tl, 8);
    histogram(ax, [T2.tracks.straightness], 0:0.05:1, 'Normalization','probability'); hold(ax,'on');
    histogram(ax, [T3.tracks.straightness], 0:0.05:1, 'Normalization','probability');
    xlabel(ax,'straightness (1=traveling, 0=standing)'); ylabel(ax,'prob');
    title(ax,'track straightness — wandering (vortex gas)'); legend(ax,{'sub-0002','sub-0003'});

    ax = nexttile(tl, 9); axis(ax,'off');
    txt = {sprintf('EIGENBASIS OMEGA flow (Dirac K=400/hemi)'), '', ...
        sprintf('regime: STANDING (speed peaks c\\rightarrow0; pk %.2f/%.2f)', Rg.rs2.peakedness, Rg.rs3.peakedness), ...
        sprintf('tracks: %d / %d (%.0f%% / %.0f%% vortex)', T2.out.nTracks, T3.out.nTracks, ...
            100*c2(1)/sum(c2), 100*c3(1)/sum(c3)), ...
        sprintf('lifetime med %.0f / %.0f ms', 1000*median(T2.out.lifetime), 1000*median(T3.out.lifetime)), ...
        sprintf('straightness med %.2f / %.2f', median([T2.tracks.straightness]), median([T3.tracks.straightness]))};
    text(ax, 0.02, 0.9, txt, 'VerticalAlignment','top','FontSize',11);

    title(tl, 'OMEGA cortical flow — Dirac-eigenbasis operators: standing vortex gas (2 subjects, full 600 s)', ...
        'FontWeight','bold','FontSize',13);

    set(findall(fig,'Type','axes'), 'Toolbar', []);   % drop interactive toolbar from export
    drawnow;
    outPng = [sp 'omega_eigenbasis_summary.png'];
    exportgraphics(fig, outPng, 'Resolution', 150);
    fprintf('figure -> %s\n', outPng);
end

% blue-white-red divergent colormap
function cm = i_bwr(n)
    if nargin<1, n=256; end
    x = linspace(0,1,n)';
    r = min(1, 2*x);  g = 1-abs(2*x-1);  b = min(1, 2-2*x);
    cm = [r g b];
end
