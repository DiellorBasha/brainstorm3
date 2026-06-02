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

% --- Resolve a usable data file by scanning all studies in the protocol ---
sProtocol = bst_get('ProtocolStudies');
if isempty(sProtocol) || ~isfield(sProtocol, 'Study') || isempty(sProtocol.Study)
    disp('SKIP: no protocol loaded.');
    return;
end
DataFile = '';
for iS = 1:numel(sProtocol.Study)
    s = sProtocol.Study(iS);
    if ~isfield(s,'iHeadModel') || isempty(s.iHeadModel) || (s.iHeadModel < 1) ...
            || (length(s.HeadModel) < s.iHeadModel) || ~isfield(s,'Data') || isempty(s.Data)
        continue;
    end
    try
        hm = in_bst_headmodel(s.HeadModel(s.iHeadModel).FileName, 0, 'HeadModelType', 'SurfaceFile');
        if strcmpi(hm.HeadModelType, 'surface')
            [~, isEig] = in_tess_eigenmodes(hm.SurfaceFile);
            if isEig
                % Raw and imported both work: the viewer reads the current
                % displayed window lazily (GetRecordingsValues).
                DataFile = s.Data(1).FileName;
                break;
            end
        end
    catch
        % Unreadable head model / surface: skip this study.
    end
end
if isempty(DataFile)
    disp('SKIP: no study with surface head model + eigenmodes + recordings.');
    return;
end

% --- Open the viewer ---
hFig = view_eigenmodes_timeseries(DataFile);
assert(~isempty(hFig) && ishandle(hFig), 'Viewer figure was not created.');
cache = getappdata(hFig, 'EigenTimeSeries');
assert(~isempty(cache) && isfield(cache, 'Theta'), 'Figure missing EigenTimeSeries cache.');
assert(size(cache.Theta,2) == numel(cache.TimeVector), 'Theta time dimension must match TimeVector.');

% --- Snapshot panel state; restore (and close figure) on any exit path ---
st0 = panel_eigenmodes('GetState');
cleanup = onCleanup(@() restorePanel(hFig, st0)); %#ok<NASGU>

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
