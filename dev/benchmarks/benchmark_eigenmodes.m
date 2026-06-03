function outDir = benchmark_eigenmodes(preset)
% BENCHMARK_EIGENMODES: Top-level driver for the eigenmode accuracy benchmark.
% USAGE:  outDir = benchmark_eigenmodes()         % full
%         outDir = benchmark_eigenmodes('smoke')  % fast preset
if nargin < 1 || isempty(preset); preset = 'full'; end
C = bench_config(preset);
outDir = C.outDir;
if ~exist(outDir,'dir'); mkdir(outDir); end

csvPath = bench_run(C);
St = bench_stats(csvPath, outDir);

renderCtx = [];
try
    renderCtx = build_render_ctx(C);
catch ME
    warning('benchmark_eigenmodes: cortex render skipped: %s', ME.message);
end

figDir = fullfile(outDir,'figures');
figFiles = bench_figures(csvPath, figDir, renderCtx);

write_report(outDir, C, St, figFiles);
fprintf('BENCH> report written to %s\n', fullfile(outDir,'REPORT.md'));
end

function renderCtx = build_render_ctx(C)
anat = C.anatomies{1};
info = bench_fixtures(anat, C.nModes_eig);
baseHM = in_bst_headmodel(info.baseHmFile, 1);
goodMask = all(isfinite(double(baseHM.Gain)),2);
L = double(baseHM.Gain(goodMask,:));
Surf = in_tess_bst(info.surfaceFile);
if ~isfield(Surf,'VertConn')||isempty(Surf.VertConn); Surf.VertConn = tess_vertconn(Surf.Vertices,Surf.Faces); end
NC = load(file_fullpath(info.ncFile)); Cnoise = NC.NoiseCov(goodMask,goodMask);
S = bst_benchmark_sources(struct('Vertices',Surf.Vertices,'VertConn',Surf.VertConn), 'focal','nTime',C.nTime,'Seed',C.seed);
Sim = bst_benchmark_simulate(L, S.Sources, Cnoise, 'SNR', 10, 'Seed', C.seed);
Est = bst_benchmark_inverse(Sim.F, info.baseHmFile, info.ncFile, info.chFile, goodMask, 10, max(C.k_total));
[~,tStar] = max(sum(S.Sources.^2,1));
estMaps = struct();
for m = C.methods
    if isfield(Est, m{1}); estMaps.(m{1}) = Est.(m{1})(:,tStar); end
end
renderCtx = struct('Vertices',Surf.Vertices,'Faces',Surf.Faces,'gt',S.GT(:), ...
    'estMaps',estMaps,'titleStr',sprintf('%s - focal source, SNR 10 dB', anat.label));
end

function write_report(outDir, C, St, figFiles)
fid = fopen(fullfile(outDir,'REPORT.md'),'w');
fprintf(fid,'# Eigenmode Source-Mapping Accuracy Benchmark - Report\n\n');
fprintf(fid,'Anatomies: %s. Methods: %s. Regimes: %s. SNR (dB): %s. K (total): %s.\n\n', ...
    strjoin(cellfun(@(a)a.label,C.anatomies,'uni',0),', '), strjoin(C.methods,', '), ...
    strjoin(C.regimes,', '), mat2str(C.snr_db), mat2str(C.k_total));

fprintf(fid,'## K-sweep (focal, eig\\_mne\\_log)\n\n| K | median LocError (mm) |\n|---|---|\n');
ks = St.ksweep(St.ksweep.regime=="focal" & St.ksweep.eig_method=="eig_mne_log", :);
ks = sortrows(ks,'K');
plateauK = NaN;
for i=1:height(ks)
    fprintf(fid,'| %d | %.2f |\n', ks.K(i), ks.median_locerror_mm(i));
    if i>1 && isnan(plateauK) && (ks.median_locerror_mm(i-1)-ks.median_locerror_mm(i) <= C.plateauTol_mm)
        plateauK = ks.K(i-1);
    end
end
if isnan(plateauK) && ~isempty(ks); plateauK = ks.K(end); end
fprintf(fid,'\nPlateau-K (focal): **%d total modes** (improvement <= %.1f mm beyond this).\n\n', plateauK, C.plateauTol_mm);

fprintf(fid,'## Competitiveness (paired Wilcoxon, eig vs standard)\n\n');
fprintf(fid,'| regime | eig | vs | median diff (mm) | p |\n|---|---|---|---|---|\n');
for i=1:height(St.compare)
    c = St.compare(i,:);
    fprintf(fid,'| %s | %s | %s | %+.2f | %.4g |\n', c.regime, c.eig_method, c.std_method, c.median_diff_mm, c.p_value);
end
fprintf(fid,'\n_Eigenmode is "competitive" where median diff is small and p is not significant, or where the diff favours eig._\n\n');

fprintf(fid,'## Figures\n\n');
for i=1:numel(figFiles)
    [~,nm,ext] = fileparts(figFiles{i});
    fprintf(fid,'![%s](figures/%s%s)\n\n', nm, nm, ext);
end
fclose(fid);
end
