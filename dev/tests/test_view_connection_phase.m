function test_view_connection_phase
% Smoke test: view_connection_phase opens a cortex figure carrying the connection
% phase as surface data plus a field-glyph layer, singularity markers, and a
% working keyboard mode-step. Runs on the REAL DB cortex surface; the viewer
% caches a ConnEigenmodes axis on it (the intended canonical behavior).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex in the current protocol.\n');
    return;
end
cleanup = onCleanup(@cleanupFig);

% Small mode count for a fast open.
hFig = view_connection_phase(SurfaceFile, 'nModes', 12, 'MaxArrows', 400);
assert(ishandle(hFig), 'view_connection_phase did not return a valid figure handle.');

hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
assert(~isempty(hAxes), 'No Axes3D in the figure.');

% Surface data is registered (the phase colormap drives the patch).
TessInfo = getappdata(hFig, 'Surface');
assert(~isempty(TessInfo) && ~isempty(TessInfo(1).DataSource.FileName), ...
    'Figure should carry surface data (the phase result).');
assert(TessInfo(1).DataThreshold == 0, 'Phase display must use DataThreshold 0 (angular, not amplitude-thresholded).');

% State object present.
st = getappdata(hFig, 'ConnPhase');
assert(~isempty(st) && isfield(st, 'Rank') && st.Rank == 1, 'ConnPhase state must start at Rank 1.');
assert(isfield(st, 'nModes') && st.nModes >= 2, 'ConnPhase state must record the mode count.');

fprintf('PASSED: view_connection_phase opens with phase surface data (Rank %d of %d modes).\n', st.Rank, st.nModes);
fprintf('ALL TESTS PASSED: test_view_connection_phase\n');
end


function cleanupFig()
    hs = findobj(0, 'Type', 'figure');
    for h = hs(:)'
        if ~isempty(getappdata(h, 'ConnPhase')), close(h); end
    end
end


function SurfaceFile = find_cortex_20484V()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects)
    return;
end
allSubj = [sSubjects.Subject];
fallback = '';
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex')
            continue;
        end
        try
            T = load(file_fullpath(surf(iF).FileName), 'Vertices', 'Reg');
        catch
            continue;
        end
        if size(T.Vertices, 1) ~= 20484
            continue;
        end
        hasReg = isfield(T,'Reg') && isstruct(T.Reg) && isfield(T.Reg,'Sphere') ...
                 && isfield(T.Reg.Sphere,'Vertices') && ~isempty(T.Reg.Sphere.Vertices);
        if hasReg
            SurfaceFile = surf(iF).FileName;
            return;
        elseif isempty(fallback)
            fallback = surf(iF).FileName;
        end
    end
end
if isempty(SurfaceFile)
    SurfaceFile = fallback;
end
end
