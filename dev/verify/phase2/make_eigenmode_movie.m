% MAKE_EIGENMODE_MOVIE: writes the first 100 LBO eigenmodes as a Brainstorm
% results file (Time axis = mode index) in the isolated protocol, and saves
% snapshot PNGs of selected modes for the Gate 2 report.
%
% Phase 2 Task 4. Run ONLY via the isolated runner:
%   dev/verify/phase0/run_matlab.sh /abs/path/make_eigenmode_movie.m
% This intentionally embeds the K=400 LBO eigenbasis into the ISOLATED
% protocol's cortex file (sub-0002/tess_cortex_pial_low.mat), growing it by
% ~64 MB — this becomes the standing Gate-2 demo protocol. The Phase 1
% oracle (dev/verify/phase1/oracle_lbo_sub0002.mat) is self-contained and
% unaffected.
nShow = 100;
SurfaceFile = 'sub-0002/tess_cortex_pial_low.mat';   % protocol-relative
E = tess_eigen(SurfaceFile, 'Laplace-Beltrami', 'nModes', max(400, nShow));
Full = file_fullpath(SurfaceFile);
T = load(Full, 'Vertices');
nV = size(T.Vertices, 1);
[rH, lH] = tess_hemisplit(in_tess_bst(SurfaceFile));
Maps = zeros(nV, nShow);
Maps(sort(lH(:)), :) = E.Phi{1}(:, 1:nShow);
Maps(sort(rH(:)), :) = E.Phi{2}(:, 1:nShow);
% results structure (source map on cortex, mode index as "time")
ResultsMat = db_template('resultsmat');
ResultsMat.Comment       = sprintf('LBO eigenmodes 1-%d (mode=frame)', nShow);
ResultsMat.ImageGridAmp  = Maps;
ResultsMat.Time          = 1:nShow;
ResultsMat.SurfaceFile   = SurfaceFile;
ResultsMat.HeadModelType = 'surface';
ResultsMat.Function      = 'eigenmodes';
ResultsMat.nComponents   = 1;
ResultsMat.DataFile      = '';
ResultsMat.ChannelFlag   = [];
% save into sub-0002 @intra
[sSubject, iSubject] = bst_get('Subject', 'sub-0002');
assert(~isempty(sSubject), 'sub-0002 not found in isolated protocol');
sStudy = bst_get('AnalysisIntraStudy', iSubject);
[~, iStudy] = bst_get('Study', sStudy.FileName);
OutFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'results_lbo_eigenmodes');
bst_save(OutFile, ResultsMat, 'v7');
db_add_data(iStudy, OutFile, ResultsMat);
db_save(1);
fprintf('EIGENMODE MOVIE: %s\n', file_short(OutFile));

% snapshots of modes 2, 10, 50 for the report (headless-fragile: best effort)
outDir = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/verify/phase2';
try
    hFig = view_surface_data(SurfaceFile, file_short(OutFile));
    for m = [2 10 50]
        panel_time('SetCurrentTime', m);
        drawnow;
        frame = bst_call(@out_figure_image, hFig);
        imwrite(frame, fullfile(outDir, sprintf('eigenmode_%03d.png', m)));
    end
    close(hFig);
    disp('SNAPSHOTS: OK');
catch err
    fprintf('SNAPSHOTS SKIPPED (headless failure): %s\n', err.message);
    disp('visual check deferred to user''s Gate-2 GUI session');
end
disp('make_eigenmode_movie DONE');
