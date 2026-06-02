function test_eigenmode_timeseries_e2e
% Smoke test: open the eigenmode time series viewer on a suitable data file,
% confirm the figure + cache, then narrow the panel band and confirm the live
% selection does not grow. Requires a loaded protocol with an imported data
% file whose study has a surface head model + computed eigenmodes. Skips
% cleanly otherwise, and restores panel state on exit.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% --- Resolve a usable study by scanning all studies in the protocol ---
sProtocol = bst_get('ProtocolStudies');
if isempty(sProtocol) || ~isfield(sProtocol, 'Study') || isempty(sProtocol.Study)
    disp('SKIP: no protocol loaded.');
    return;
end
% Find a study with surface head model + eigenmodes + noise cov + a data file
found = [];
for iS = 1:numel(sProtocol.Study)
    s = sProtocol.Study(iS);
    if ~isfield(s,'iHeadModel') || isempty(s.iHeadModel) || (s.iHeadModel < 1) ...
            || (length(s.HeadModel) < s.iHeadModel) || ~isfield(s,'Data') || isempty(s.Data) ...
            || ~isfield(s,'NoiseCov') || isempty(s.NoiseCov) || isempty(s.NoiseCov(1).FileName)
        continue;
    end
    try
        hm = in_bst_headmodel(s.HeadModel(s.iHeadModel).FileName, 0, 'HeadModelType', 'SurfaceFile');
        if strcmpi(hm.HeadModelType, 'surface')
            [~, isEig] = in_tess_eigenmodes(hm.SurfaceFile);
            if isEig
                found = struct('hmFile', s.HeadModel(s.iHeadModel).FileName, ...
                    'surf', hm.SurfaceFile, 'ncFile', s.NoiseCov(1).FileName, ...
                    'dataFile', s.Data(1).FileName);
                break;
            end
        end
    catch
    end
end
if isempty(found)
    disp('SKIP: no study with surface head model + eigenmodes + noise cov + data.');
    return;
end
[sStudyT, iStudyT] = bst_get('AnyFile', found.dataFile);

% Build a temporary Harmonic node via the engine (mirrors what the GUI saves)
HM = in_bst_headmodel(found.hmFile, 1);                 % constrained [nch x nVert]
goodMask = all(isfinite(double(HM.Gain)), 2);
[InvE, errE] = bst_inverse_eigenmodes(found.hmFile, found.surf, found.ncFile, ...
    'Method', 'harmonic', 'nModes', 0, 'GoodChannel', goodMask);
assert(isempty(errE), ['harmonic inverse failed: ' errE]);
Kh   = InvE.nModes;
[EigT, ~] = in_tess_eigenmodes(found.surf);
PhiT = double(EigT.Vectors(:, 1:Kh));
ResMat = db_template('resultsmat');
ResMat.ImagingKernel = PhiT * InvE.ImagingKernel;       % [nVert x nGoodCh]
ResMat.ImageGridAmp  = [];
ResMat.nComponents   = 1;
ResMat.Function      = 'eigenmode_harmonic';
ResMat.EigenKernel   = InvE.ImagingKernel;              % M̃ [K x nGoodCh]
ResMat.GoodChannel   = InvE.GoodChannel;
ResMat.Whitener      = InvE.Whitener;
ResMat.SurfaceFile   = found.surf;
ResMat.HeadModelFile = found.hmFile;
ResMat.HeadModelType = 'surface';
ResMat.DataFile      = found.dataFile;
ResMat.Comment       = 'TEST Eigenmode HARMONIC';
ResMat.Time          = [];
StudyDir   = bst_fileparts(file_fullpath(sStudyT.FileName));
NodeFile   = bst_process('GetNewFilename', StudyDir, 'results_TESTharmonic');
bst_save(NodeFile, ResMat, 'v6');
db_add_data(iStudyT, NodeFile, ResMat);
nodeRel = file_short(NodeFile);
% Always clean up the temp node + figure, however the test exits
cleanupNode = onCleanup(@() cleanupHarmonic(NodeFile, iStudyT)); %#ok<NASGU>

% --- Open the viewer FROM the Harmonic node ---
hFig = view_eigenmodes_timeseries(nodeRel);
assert(~isempty(hFig) && ishandle(hFig), 'Viewer figure was not created.');
cache = getappdata(hFig, 'EigenTimeSeries');
assert(~isempty(cache) && isfield(cache, 'Theta'), 'Figure missing EigenTimeSeries cache.');
assert(size(cache.Theta,2) == numel(cache.TimeVector), 'Theta time dimension must match TimeVector.');

% Consistency: the cached Theta must equal M̃ * D for the current window
iDSt = bst_memory('LoadDataFile', found.dataFile);
Dwin = bst_memory('GetRecordingsValues', iDSt, cache.GoodChannel, 'UserTimeWindow', 0);
assert(max(abs(cache.Theta(:) - reshape(InvE.ImagingKernel * Dwin, [], 1))) < 1e-6 * (1 + max(abs(cache.Theta(:)))), ...
    'Cached Theta must equal EigenKernel * D(window).');

% --- Snapshot panel state; restore (and close figure) on any exit path ---
st0 = panel_eigenmodes('GetState');
cleanup = onCleanup(@() restorePanel(hFig, st0)); %#ok<NASGU>

% --- Regression: a figure reload (butterfly<->column toggle, montage change,
% "reload all") must replot the cached coefficients, NOT reload the raw file via
% view_matrix (which loaded the whole recording / errored on the missing .Value).
TsInfoR = getappdata(hFig, 'TsInfo');
TsInfoR.DisplayMode = 'column';
setappdata(hFig, 'TsInfo', TsInfoR);
bst_figures('ReloadFigures', hFig, 0);
assert(ishandle(hFig), 'Figure must survive a display-mode reload.');
cacheR = getappdata(hFig, 'EigenTimeSeries');
assert(~isempty(cacheR) && isfield(cacheR, 'Theta'), 'Cache must survive reload.');
assert(strcmp(getfield(getappdata(hFig, 'TsInfo'), 'DisplayMode'), 'column'), ...
    'Reload must keep the toggled column display mode.'); %#ok<GFLD>

% --- Regression: a raw PAGE change reloads via ReloadFigures with the cached
% window now stale (pointing at the old page). The reload must re-read the
% current window, not leave the figure on the previous page (time-0 revert).
cacheStale = getappdata(hFig, 'EigenTimeSeries');
realWin = cacheStale.WindowTime;            % == Measures.Time of the live window
cacheStale.WindowTime = realWin + 1000;     % simulate a page change left it stale
setappdata(hFig, 'EigenTimeSeries', cacheStale);
bst_figures('ReloadFigures', hFig, 0);
cacheSync = getappdata(hFig, 'EigenTimeSeries');
assert(isequal(cacheSync.WindowTime, realWin), ...
    'Reload must re-sync the cached window to the live recording (page tracking).');
assert(numel(cacheSync.TimeVector) == size(cacheSync.Theta, 2), ...
    'Re-read Theta and TimeVector must stay consistent.');

% Expected trace count for the current band
iRows0 = view_eigenmodes_timeseries('GetBandTraces', cache.Component, cache.CompRank, st0.Band(1), st0.Band(2));

% Narrow the band to a single rank (SetBand sets the span unambiguously) and
% confirm the live tracking does not grow the trace count.
panel_eigenmodes('SetBand', st0.Band(1), st0.Band(1));
st1 = panel_eigenmodes('GetState');
assert(st1.Band(1) == st1.Band(2), 'SetBand(c,c) must collapse the band to a point.');
iRows1 = view_eigenmodes_timeseries('GetBandTraces', cache.Component, cache.CompRank, st1.Band(1), st1.Band(2));
assert(numel(iRows1) <= numel(iRows0), 'Single-rank band should not increase trace count.');

disp('ALL TESTS PASSED');
end


% Restore panel shape/band and close the test figure, no matter how the test exits.
function restorePanel(hFig, st0)
    if ~isempty(hFig) && ishandle(hFig)
        close(hFig);
    end
    try
        panel_eigenmodes('SetWindowShape', st0.WindowShape);
        panel_eigenmodes('SetBand', st0.Band(1), st0.Band(2));
    catch
        % Best-effort restore.
    end
end

% Delete the temporary Harmonic node and reload the study (DB stays consistent).
function cleanupHarmonic(NodeFile, iStudy)
    try
        if exist(file_fullpath(NodeFile), 'file')
            file_delete(file_fullpath(NodeFile), 1);
            db_reload_studies(iStudy);
        end
    catch
        % Best-effort cleanup.
    end
end
