function csvPath = bench_run(C)
% BENCH_RUN: Execute the synthetic-on-real-cortex benchmark; write synthetic.csv.
% USAGE:  csvPath = bench_run(bench_config())
if ~exist(C.outDir,'dir'); mkdir(C.outDir); end
csvPath = fullfile(C.outDir, 'synthetic.csv');

rows = {};   % {anatomy, regime, snr_db, replicate, method, K, locerror_mm, corr, nrmse, auc, dispersion_mm}
for ia = 1:numel(C.anatomies)
    anat = C.anatomies{ia};
    info = bench_fixtures(anat, C.nModes_eig);

    baseHM = in_bst_headmodel(info.baseHmFile, 1);          % [nCh x nVert]
    goodMask = all(isfinite(double(baseHM.Gain)),2);
    L = double(baseHM.Gain(goodMask,:));
    Surf = in_tess_bst(info.surfaceFile);
    if ~isfield(Surf,'VertConn') || isempty(Surf.VertConn)
        Surf.VertConn = tess_vertconn(Surf.Vertices, Surf.Faces);
    end
    GridLoc = Surf.Vertices;                                 % metres
    Surface = struct('Vertices',Surf.Vertices,'VertConn',Surf.VertConn);
    NC = load(file_fullpath(info.ncFile)); Cnoise = NC.NoiseCov(goodMask,goodMask);

    for ir = 1:numel(C.regimes)
        regime = C.regimes{ir};
        ropts  = C.regimeOpts.(regime);
        for isnr = 1:numel(C.snr_db)
            snr = C.snr_db(isnr);
            for rep = 1:C.nReplicates
                seed = mod(C.seed + 1000*ia + 100*ir + 10*isnr + rep, 2^31-1);
                S = bst_benchmark_sources(Surface, regime, 'nTime', C.nTime, 'Seed', seed, ropts{:});
                Sim = bst_benchmark_simulate(L, S.Sources, Cnoise, 'SNR', snr, 'Seed', seed);
                [~, tStar] = max(sum(S.Sources.^2,1));
                gt = S.GT(:);

                for jk = 1:numel(C.k_total)
                    Ktot = C.k_total(jk);
                    try
                        Est = bst_benchmark_inverse(Sim.F, info.baseHmFile, info.ncFile, ...
                            info.chFile, goodMask, snr, Ktot);
                    catch ME
                        warning('bench_run: inverse failed (%s/%s/snr%d/rep%d/K%d): %s', ...
                            anat.label, regime, snr, rep, Ktot, ME.message);
                        continue;
                    end
                    if jk == 1; theseMethods = C.methods; else; theseMethods = C.eigMethods; end
                    for im = 1:numel(theseMethods)
                        meth = theseMethods{im};
                        if ~isfield(Est, meth); continue; end
                        estMap = Est.(meth)(:, tStar);
                        if ~all(isfinite(estMap)); continue; end
                        M = bst_benchmark_metrics(gt, estMap, GridLoc, S.SeedVertex);
                        if ismember(meth, C.eigMethods); Kcol = Ktot; else; Kcol = NaN; end
                        rows(end+1,:) = {anat.label, regime, snr, rep, meth, Kcol, ...
                            M.LocError, M.Correlation, M.NRMSE, M.AUC, M.SpatialDispersion}; %#ok<AGROW>
                    end
                end
            end
        end
        fprintf('BENCH> %s / %s done (%d rows so far)\n', anat.label, regime, size(rows,1));
    end
end

hdr = {'anatomy','regime','snr_db','replicate','method','K','locerror_mm', ...
       'correlation','nrmse','auc','spatial_dispersion_mm'};
fid = fopen(csvPath,'w');
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:size(rows,1)
    fprintf(fid, '%s,%s,%g,%d,%s,%g,%g,%g,%g,%g,%g\n', rows{i,1}, rows{i,2}, rows{i,3}, ...
        rows{i,4}, rows{i,5}, rows{i,6}, rows{i,7}, rows{i,8}, rows{i,9}, rows{i,10}, rows{i,11});
end
fclose(fid);
fprintf('BENCH> wrote %d rows to %s\n', size(rows,1), csvPath);
end
