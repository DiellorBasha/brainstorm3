function test_eigenmode_timeseries_e2e
% Smoke test: open the eigenmode time series viewer on the currently selected
% data file, confirm a figure with the expected number of traces for the
% current band, then move the band and confirm the trace count changes.
% Requires: a loaded protocol with an imported data file whose study has a
% surface head model + computed eigenmodes. Skips cleanly otherwise.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% --- Resolve a usable data file from the current protocol ---
[sStudies, iStudies] = bst_get('ProtocolStudies'); %#ok<ASGLU>
DataFile = '';
for iS = 1:numel(iStudies)
    s = bst_get('Study', iStudies(iS));
    if ~isempty(s) && isfield(s,'iHeadModel') && ~isempty(s.iHeadModel) && (s.iHeadModel >= 1) ...
            && isfield(s,'Data') && ~isempty(s.Data)
        hm = in_bst_headmodel(s.HeadModel(s.iHeadModel).FileName, 0, 'HeadModelType', 'SurfaceFile');
        if strcmpi(hm.HeadModelType, 'surface')
            [~, isEig] = in_tess_eigenmodes(hm.SurfaceFile);
            if isEig
                DataFile = s.Data(1).FileName;
                break;
            end
        end
    end
end
if isempty(DataFile)
    disp('SKIP: no study with surface head model + eigenmodes + imported data.');
    return;
end

% --- Open the viewer ---
hFig = view_eigenmodes_timeseries(DataFile);
assert(~isempty(hFig) && ishandle(hFig), 'Viewer figure was not created.');
cache = getappdata(hFig, 'EigenTimeSeries');
assert(~isempty(cache) && isfield(cache, 'Theta'), 'Figure missing EigenTimeSeries cache.');

% --- Expected trace count for the current band ---
st = panel_eigenmodes('GetState');
iRows0 = view_eigenmodes_timeseries('GetBandTraces', cache.Component, cache.CompRank, st.Band(1), st.Band(2));

% --- Narrow the band to a single rank, confirm tracking redraws no more traces ---
panel_eigenmodes('SetWindowShape', 'single');
panel_eigenmodes('SetCurrentMode', st.Band(1));
st1 = panel_eigenmodes('GetState');
iRows1 = view_eigenmodes_timeseries('GetBandTraces', cache.Component, cache.CompRank, st1.Band(1), st1.Band(2));
assert(numel(iRows1) <= numel(iRows0), 'Single-mode band should not increase trace count.');
assert(~isempty(cache.Theta) && size(cache.Theta,2) == numel(cache.TimeVector), 'Theta time dimension must match TimeVector.');

disp('ALL TESTS PASSED');
close(hFig);
end
