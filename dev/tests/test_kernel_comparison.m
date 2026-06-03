%% compare_kernels.m
%  Benchmark comparison of two ImagingKernels from Brainstorm:
%  standard MNE vs eigenmode-based MNE.
%
%  Produces four figures:
%    Fig 1 — Global summary (correlation, norm ratio, weight distributions)
%    Fig 2 — Per-vertex spatial filter comparison (cortical maps)
%    Fig 3 — Singular value spectrum comparison
%    Fig 4 — Resolution matrix analysis (PSF width, peak amplitude)
%
%  Requirements: Brainstorm running with the protocol loaded.
%  Set the two ResultsFile paths in the configuration section below.

%% ========================================================================
%  CONFIGURATION
%  ========================================================================

% Relative paths within the Brainstorm protocol DB
ResultsFile_std = 'Subject01/S01_AEF_20131218_01_notch/results_MN_MEG_KERNEL_260603_0647.mat';
ResultsFile_eig = 'Subject01/S01_AEF_20131218_01_notch/results_Eigen-MNE_MEG_KERNEL_260603_0700.mat';

% Output directory for figures
outputDir = 'C:\Users\diell\workspace\research\code\brainstorm3\dev\tests';

% Export format: 'pdf' (vector) or 'png' (raster, 300 dpi)
exportFmt = 'png';

%% ========================================================================
%  DESIGN SYSTEM DEFAULTS
%  ========================================================================

set(groot, 'DefaultAxesFontName',  'Helvetica');
set(groot, 'DefaultAxesFontSize',  10);
set(groot, 'DefaultTextFontName',  'Helvetica');
set(groot, 'DefaultAxesTickDir',   'out');
set(groot, 'DefaultAxesBox',       'off');
set(groot, 'DefaultAxesLineWidth', 0.75);

col_std = [0.20 0.40 0.80];   % blue  — standard MNE
col_eig = [0.85 0.33 0.10];   % orange — eigen-MNE
col_gray = [0.60 0.60 0.60];

if isempty(outputDir)
    outputDir = pwd;
end

%% ========================================================================
%  LOAD DATA
%  ========================================================================

fprintf('Loading standard MNE results...\n');
ResMat_std = in_bst_results(ResultsFile_std, 0);

fprintf('Loading eigen-MNE results...\n');
ResMat_eig = in_bst_results(ResultsFile_eig, 0);

K_std = ResMat_std.ImagingKernel;   % [nSources x nChannels]
K_eig = ResMat_eig.ImagingKernel;

[nSources, nChannels] = size(K_std);
assert(isequal(size(K_std), size(K_eig)), ...
    'Kernel dimensions do not match: std=[%d x %d], eig=[%d x %d]', ...
    size(K_std,1), size(K_std,2), size(K_eig,1), size(K_eig,2));

nComp = ResMat_std.nComponents;
assert(nComp == 1, ...
    'This script assumes constrained orientation (nComponents=1). Got %d.', nComp);

fprintf('  Kernel size: %d sources x %d channels\n', nSources, nChannels);

% --- Load cortical surface for mapping ---
SurfaceFile = ResMat_std.SurfaceFile;
SurfMat = in_tess_bst(SurfaceFile);
Vertices = SurfMat.Vertices;
Faces    = SurfMat.Faces;

fprintf('  Surface: %d vertices, %d faces\n', size(Vertices,1), size(Faces,1));

% --- Load head model (leadfield) for resolution matrix ---
HeadModelFile = ResMat_std.HeadModelFile;
HeadMat = in_bst_headmodel(HeadModelFile);
Gain_full = HeadMat.Gain;   % [nChannels x 3*nSources] unconstrained

% Project leadfield to constrained orientation using surface normals
GridOrient = HeadMat.GridOrient;
if isempty(GridOrient)
    % Compute from surface normals
    GridOrient = SurfMat.VertNormals;
end

% Constrained leadfield: G(:,j) = sum over xyz of Gain(:,3*(j-1)+k) * orient(j,k)
nChannels_all = size(Gain_full, 1);
G = zeros(nChannels_all, nSources);
for j = 1:nSources
    idx3 = 3*(j-1) + (1:3);
    G(:,j) = Gain_full(:,idx3) * GridOrient(j,:)';
end

% Restrict leadfield to good channels used by the kernel
GoodChannel = ResMat_std.GoodChannel;
if ~isempty(GoodChannel)
    G = G(GoodChannel, :);
end

fprintf('  Constrained leadfield: [%d x %d]\n', size(G,1), size(G,2));

%% ========================================================================
%  FIGURE 1 — GLOBAL SUMMARY
%  ========================================================================

figWidth_mm  = 170;   % double column
figHeight_mm = 110;
fig1 = figure('Units', 'points', ...
    'Position', [100 100 figWidth_mm*3.78 figHeight_mm*3.78], ...
    'Color', 'w', 'PaperPositionMode', 'auto');

t1 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel A: Flattened kernel scatter + correlation ---
nexttile
k_std_flat = K_std(:);
k_eig_flat = K_eig(:);

% Subsample for plotting (full scatter is too dense)
nPts = numel(k_std_flat);
if nPts > 50000
    idx_sub = randsample(nPts, 50000);
else
    idx_sub = 1:nPts;
end

scatter(k_std_flat(idx_sub), k_eig_flat(idx_sub), 1, col_gray, ...
    'filled', 'MarkerFaceAlpha', 0.15);
hold on
% Unity line
lims = [min([k_std_flat; k_eig_flat]), max([k_std_flat; k_eig_flat])];
plot(lims, lims, 'k--', 'LineWidth', 0.75);
hold off

r_pearson = corr(k_std_flat, k_eig_flat);
text(0.05, 0.93, sprintf('r = %.4f', r_pearson), ...
    'Units', 'normalized', 'FontSize', 9);

xlabel('Standard MNE weights')
ylabel('Eigen-MNE weights')
title('Kernel weight correlation')
text(-0.12, 1.05, 'A', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

% --- Panel B: Frobenius norm comparison ---
nexttile
fro_std = norm(K_std, 'fro');
fro_eig = norm(K_eig, 'fro');

bar_vals = [fro_std, fro_eig];
b = bar(bar_vals, 0.6);
b.FaceColor = 'flat';
b.CData = [col_std; col_eig];
set(gca, 'XTickLabel', {'Std MNE', 'Eigen-MNE'});
ylabel('Frobenius norm')
title('Global kernel magnitude')

text(0.05, 0.93, sprintf('Ratio = %.3f', fro_eig / fro_std), ...
    'Units', 'normalized', 'FontSize', 9);
text(-0.12, 1.05, 'B', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

% --- Panel C: Weight distribution histograms ---
nexttile
nBins = 200;
edges = linspace(lims(1), lims(2), nBins+1);

h_std = histcounts(k_std_flat, edges, 'Normalization', 'probability');
h_eig = histcounts(k_eig_flat, edges, 'Normalization', 'probability');
centers = 0.5 * (edges(1:end-1) + edges(2:end));

semilogy(centers, h_std, 'Color', col_std, 'LineWidth', 1.2);
hold on
semilogy(centers, h_eig, 'Color', col_eig, 'LineWidth', 1.2);
hold off

xlabel('Kernel weight')
ylabel('Probability (log)')
legend({'Std MNE', 'Eigen-MNE'}, 'Location', 'northeast', 'FontSize', 8);
title('Weight distributions')
text(-0.12, 1.05, 'C', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

exportFigure(fig1, fullfile(outputDir, 'fig1_kernel_global_summary'), exportFmt);

%% ========================================================================
%  FIGURE 2 — PER-VERTEX SPATIAL FILTER COMPARISON (CORTICAL MAPS)
%  ========================================================================

% Cosine similarity between corresponding rows (spatial filters)
cos_sim = zeros(nSources, 1);
norm_ratio = zeros(nSources, 1);
for v = 1:nSources
    row_std = K_std(v, :);
    row_eig = K_eig(v, :);
    
    n_std = norm(row_std);
    n_eig = norm(row_eig);
    
    if n_std > 0 && n_eig > 0
        cos_sim(v) = dot(row_std, row_eig) / (n_std * n_eig);
    else
        cos_sim(v) = NaN;
    end
    
    if n_std > 0
        norm_ratio(v) = n_eig / n_std;
    else
        norm_ratio(v) = NaN;
    end
end

figWidth_mm  = 170;
figHeight_mm = 160;
fig2 = figure('Units', 'points', ...
    'Position', [100 100 figWidth_mm*3.78 figHeight_mm*3.78], ...
    'Color', 'w', 'PaperPositionMode', 'auto');

t2 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel A: Cosine similarity — lateral left ---
nexttile
plotCorticalMap(Vertices, Faces, cos_sim, [-90, 0]);
cb = colorbar;
ylabel(cb, 'Cosine similarity');
colormap(gca, parula(256));
clim([min(cos_sim), 1]);
title('Filter similarity (lateral L)');
text(-0.05, 1.02, 'A', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

% --- Panel B: Cosine similarity — lateral right ---
nexttile
plotCorticalMap(Vertices, Faces, cos_sim, [90, 0]);
cb = colorbar;
ylabel(cb, 'Cosine similarity');
colormap(gca, parula(256));
clim([min(cos_sim), 1]);
title('Filter similarity (lateral R)');
text(-0.05, 1.02, 'B', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

% --- Panel C: Norm ratio — lateral left ---
nexttile
plotCorticalMap(Vertices, Faces, log2(norm_ratio), [-90, 0]);
cb = colorbar;
ylabel(cb, 'log_2(norm ratio)');
max_abs = max(abs(log2(norm_ratio(isfinite(norm_ratio)))));
cmap_div = divergingColormap(256);
colormap(gca, cmap_div);
clim([-max_abs, max_abs]);
title('Sensitivity redistribution (lateral L)');
text(-0.05, 1.02, 'C', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

% --- Panel D: Norm ratio — lateral right ---
nexttile
plotCorticalMap(Vertices, Faces, log2(norm_ratio), [90, 0]);
cb = colorbar;
ylabel(cb, 'log_2(norm ratio)');
colormap(gca, cmap_div);
clim([-max_abs, max_abs]);
title('Sensitivity redistribution (lateral R)');
text(-0.05, 1.02, 'D', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

exportFigure(fig2, fullfile(outputDir, 'fig2_kernel_cortical_maps'), exportFmt);

% Print summary stats
fprintf('\n--- Per-vertex filter comparison ---\n');
fprintf('  Cosine similarity:  mean=%.4f, median=%.4f, min=%.4f\n', ...
    nanmean(cos_sim), nanmedian(cos_sim), nanmin(cos_sim));
fprintf('  Norm ratio:         mean=%.4f, median=%.4f\n', ...
    nanmean(norm_ratio), nanmedian(norm_ratio));
fprintf('  Vertices with cos_sim < 0.9:  %d / %d (%.1f%%)\n', ...
    sum(cos_sim < 0.9), nSources, 100*sum(cos_sim < 0.9)/nSources);

%% ========================================================================
%  FIGURE 3 — SINGULAR VALUE SPECTRUM
%  ========================================================================

fprintf('\nComputing SVD of kernels...\n');
[~, S_std, ~] = svd(K_std, 'econ');
[~, S_eig, ~] = svd(K_eig, 'econ');

sv_std = diag(S_std);
sv_eig = diag(S_eig);

% Effective rank at 99% energy
energy_std = cumsum(sv_std.^2) / sum(sv_std.^2);
energy_eig = cumsum(sv_eig.^2) / sum(sv_eig.^2);
rank99_std = find(energy_std >= 0.99, 1, 'first');
rank99_eig = find(energy_eig >= 0.99, 1, 'first');

figWidth_mm  = 170;
figHeight_mm = 100;
fig3 = figure('Units', 'points', ...
    'Position', [100 100 figWidth_mm*3.78 figHeight_mm*3.78], ...
    'Color', 'w', 'PaperPositionMode', 'auto');

t3 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel A: Singular value spectrum (log scale) ---
nexttile
semilogy(sv_std, 'Color', col_std, 'LineWidth', 1.5);
hold on
semilogy(sv_eig, 'Color', col_eig, 'LineWidth', 1.5);
hold off
xlabel('Singular value index')
ylabel('Singular value (log)')
legend({'Std MNE', 'Eigen-MNE'}, 'Location', 'northeast', 'FontSize', 8);
title('SV spectrum');
text(-0.12, 1.05, 'A', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

% --- Panel B: Cumulative energy ---
nexttile
plot(energy_std, 'Color', col_std, 'LineWidth', 1.5);
hold on
plot(energy_eig, 'Color', col_eig, 'LineWidth', 1.5);
yline(0.99, 'k--', 'LineWidth', 0.5);
% Mark effective ranks
plot(rank99_std, 0.99, 'o', 'Color', col_std, 'MarkerSize', 6, 'LineWidth', 1.5);
plot(rank99_eig, 0.99, 'o', 'Color', col_eig, 'MarkerSize', 6, 'LineWidth', 1.5);
hold off
xlabel('Singular value index')
ylabel('Cumulative energy fraction')
title('Effective rank');

text(0.50, 0.25, sprintf('Std:  rank_{99} = %d', rank99_std), ...
    'Units', 'normalized', 'FontSize', 9, 'Color', col_std);
text(0.50, 0.15, sprintf('Eig:  rank_{99} = %d', rank99_eig), ...
    'Units', 'normalized', 'FontSize', 9, 'Color', col_eig);
text(-0.12, 1.05, 'B', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

% --- Panel C: SV ratio ---
nexttile
nSV = min(length(sv_std), length(sv_eig));
sv_ratio = sv_eig(1:nSV) ./ sv_std(1:nSV);
plot(sv_ratio, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
hold on
yline(1, 'k--', 'LineWidth', 0.5);
hold off
xlabel('Singular value index')
ylabel('\sigma_{eig} / \sigma_{std}')
title('SV ratio (eigen / standard)');
text(-0.12, 1.05, 'C', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

exportFigure(fig3, fullfile(outputDir, 'fig3_kernel_svd_spectrum'), exportFmt);

fprintf('  Effective rank (99%% energy): std=%d, eig=%d\n', rank99_std, rank99_eig);

%% ========================================================================
%  FIGURE 4 — RESOLUTION MATRIX ANALYSIS
%  ========================================================================

fprintf('\nComputing resolution matrices...\n');

R_std = K_std * G;   % [nSources x nSources]
R_eig = K_eig * G;

% --- Diagonal: peak amplitude at true location ---
diag_std = diag(R_std);
diag_eig = diag(R_eig);

% --- PSF width: for each source, FWHM of its resolution row ---
%     Measured as count of vertices where |R(i,:)| > 0.5 * |R(i,i)|
fprintf('Computing PSF widths...\n');
psf_width_std = zeros(nSources, 1);
psf_width_eig = zeros(nSources, 1);

for v = 1:nSources
    row_std = abs(R_std(v, :));
    row_eig = abs(R_eig(v, :));
    
    psf_width_std(v) = sum(row_std > 0.5 * row_std(v));
    psf_width_eig(v) = sum(row_eig > 0.5 * row_eig(v));
end

% --- Cross-talk: for each source, column energy excluding diagonal ---
fprintf('Computing cross-talk functions...\n');
ctf_std = zeros(nSources, 1);
ctf_eig = zeros(nSources, 1);

for v = 1:nSources
    col_std = R_std(:, v);
    col_eig = R_eig(:, v);
    
    col_std(v) = 0;   % exclude self
    col_eig(v) = 0;
    
    ctf_std(v) = norm(col_std);
    ctf_eig(v) = norm(col_eig);
end

figWidth_mm  = 170;
figHeight_mm = 200;
fig4 = figure('Units', 'points', ...
    'Position', [100 100 figWidth_mm*3.78 figHeight_mm*3.78], ...
    'Color', 'w', 'PaperPositionMode', 'auto');

t4 = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel A: Resolution diagonal, std (lateral L) ---
nexttile
plotCorticalMap(Vertices, Faces, diag_std, [-90, 0]);
cb = colorbar;
ylabel(cb, 'R_{ii}');
colormap(gca, parula(256));
clim([0, max([diag_std; diag_eig])]);
title('Std MNE — resolution diagonal');
text(-0.05, 1.02, 'A', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

% --- Panel B: Resolution diagonal, eigen (lateral L) ---
nexttile
plotCorticalMap(Vertices, Faces, diag_eig, [-90, 0]);
cb = colorbar;
ylabel(cb, 'R_{ii}');
colormap(gca, parula(256));
clim([0, max([diag_std; diag_eig])]);
title('Eigen-MNE — resolution diagonal');
text(-0.05, 1.02, 'B', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

% --- Panel C: PSF width, std (lateral L) ---
nexttile
plotCorticalMap(Vertices, Faces, psf_width_std, [-90, 0]);
cb = colorbar;
ylabel(cb, 'FWHM (# vertices)');
colormap(gca, flipud(parula(256)));   % inverted: fewer = better = bright
clim([0, max([psf_width_std; psf_width_eig])]);
title('Std MNE — PSF width');
text(-0.05, 1.02, 'C', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

% --- Panel D: PSF width, eigen (lateral L) ---
nexttile
plotCorticalMap(Vertices, Faces, psf_width_eig, [-90, 0]);
cb = colorbar;
ylabel(cb, 'FWHM (# vertices)');
colormap(gca, flipud(parula(256)));
clim([0, max([psf_width_std; psf_width_eig])]);
title('Eigen-MNE — PSF width');
text(-0.05, 1.02, 'D', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

% --- Panel E: Resolution diagonal scatter ---
nexttile
scatter(diag_std, diag_eig, 4, col_gray, 'filled', 'MarkerFaceAlpha', 0.3);
hold on
lims_diag = [min([diag_std; diag_eig]), max([diag_std; diag_eig])];
plot(lims_diag, lims_diag, 'k--', 'LineWidth', 0.75);
hold off
xlabel('Std MNE  R_{ii}')
ylabel('Eigen-MNE  R_{ii}')
title('Diagonal comparison');
r_diag = corr(diag_std, diag_eig);
text(0.05, 0.93, sprintf('r = %.4f', r_diag), 'Units', 'normalized', 'FontSize', 9);
text(-0.12, 1.05, 'E', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

% --- Panel F: PSF width scatter ---
nexttile
scatter(psf_width_std, psf_width_eig, 4, col_gray, 'filled', 'MarkerFaceAlpha', 0.3);
hold on
lims_psf = [0, max([psf_width_std; psf_width_eig])];
plot(lims_psf, lims_psf, 'k--', 'LineWidth', 0.75);
hold off
xlabel('Std MNE  FWHM (# vertices)')
ylabel('Eigen-MNE  FWHM (# vertices)')
title('PSF width comparison');
r_psf = corr(psf_width_std, psf_width_eig);
text(0.05, 0.93, sprintf('r = %.4f', r_psf), 'Units', 'normalized', 'FontSize', 9);
text(-0.12, 1.05, 'F', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

exportFigure(fig4, fullfile(outputDir, 'fig4_resolution_matrix'), exportFmt);

%% ========================================================================
%  SUMMARY TABLE (printed to console)
%  ========================================================================

fprintf('\n========================================\n');
fprintf('  KERNEL COMPARISON SUMMARY\n');
fprintf('========================================\n');
fprintf('  Sources: %d   Channels: %d\n', nSources, nChannels);
fprintf('\n  GLOBAL:\n');
fprintf('    Pearson r (flattened):       %.6f\n', r_pearson);
fprintf('    Frobenius norm std:          %.4e\n', fro_std);
fprintf('    Frobenius norm eig:          %.4e\n', fro_eig);
fprintf('    Frobenius ratio (eig/std):   %.4f\n', fro_eig / fro_std);
fprintf('\n  PER-VERTEX FILTERS:\n');
fprintf('    Cosine similarity:  mean=%.4f  median=%.4f  min=%.4f\n', ...
    nanmean(cos_sim), nanmedian(cos_sim), nanmin(cos_sim));
fprintf('    Norm ratio:         mean=%.4f  median=%.4f\n', ...
    nanmean(norm_ratio), nanmedian(norm_ratio));
fprintf('    Vertices cos<0.9:   %d / %d (%.1f%%)\n', ...
    sum(cos_sim < 0.9), nSources, 100*sum(cos_sim < 0.9)/nSources);
fprintf('\n  SVD SPECTRUM:\n');
fprintf('    Effective rank (99%%):  std=%d  eig=%d\n', rank99_std, rank99_eig);
fprintf('    Top SV ratio:          %.4f\n', sv_eig(1) / sv_std(1));
fprintf('\n  RESOLUTION MATRIX:\n');
fprintf('    Diagonal R_ii:  std mean=%.4f  eig mean=%.4f\n', ...
    mean(diag_std), mean(diag_eig));
fprintf('    PSF width:      std mean=%.1f  eig mean=%.1f  vertices\n', ...
    mean(psf_width_std), mean(psf_width_eig));
fprintf('    Cross-talk:     std mean=%.4e  eig mean=%.4e\n', ...
    mean(ctf_std), mean(ctf_eig));
fprintf('    Diagonal r:     %.4f\n', r_diag);
fprintf('========================================\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function plotCorticalMap(Vertices, Faces, data, viewAngle)
% Render scalar data on a cortical surface patch.
    patch('Faces', Faces, 'Vertices', Vertices, ...
        'FaceVertexCData', data, ...
        'FaceColor', 'interp', 'EdgeColor', 'none', ...
        'FaceLighting', 'gouraud');
    material dull
    camlight headlight
    axis equal off
    view(viewAngle);
    set(gca, 'Color', [0.94 0.94 0.94]);
end

function cmap = divergingColormap(n)
% Blue-white-red diverging colormap (RdBu reversed).
% Attempt brewermap, fall back to manual interpolation.
    try
        cmap = flipud(brewermap(n, 'RdBu'));
    catch
        % Manual blue → white → red
        half = floor(n/2);
        blue = [0.20 0.40 0.80];
        red  = [0.80 0.20 0.20];
        white = [1 1 1];
        
        r1 = linspace(blue(1), white(1), half)';
        g1 = linspace(blue(2), white(2), half)';
        b1 = linspace(blue(3), white(3), half)';
        
        r2 = linspace(white(1), red(1), n-half)';
        g2 = linspace(white(2), red(2), n-half)';
        b2 = linspace(white(3), red(3), n-half)';
        
        cmap = [r1 g1 b1; r2 g2 b2];
    end
end

function exportFigure(fig, baseName, fmt)
% Export figure in the requested format.
    switch lower(fmt)
        case 'pdf'
            exportgraphics(fig, [baseName '.pdf'], 'ContentType', 'vector');
            fprintf('  Exported: %s.pdf\n', baseName);
        case 'png'
            exportgraphics(fig, [baseName '.png'], 'Resolution', 300);
            fprintf('  Exported: %s.png\n', baseName);
        otherwise
            warning('Unknown export format "%s". Exporting as PDF.', fmt);
            exportgraphics(fig, [baseName '.pdf'], 'ContentType', 'vector');
    end
end
