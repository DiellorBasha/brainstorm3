function test_view_tangents
% Integration test for view_tangents: auto-compute path, render, frame
% geometry (via the drawn quiver objects), and interactive callbacks.
% Runs on a temp COPY so the DB surface is not mutated.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% --- Find a low-res FreeSurfer-registered cortex in the DB ---
srcFile = find_registered_cortex();
assert(~isempty(srcFile), 'No FreeSurfer-registered cortex found in the DB to test against.');
fprintf('Source cortex: %s\n', srcFile);

% --- Work on a temp copy with NO stored frame (exercise auto-compute) ---
tmpFile = fullfile(tempdir, 'tess_cortex_viewtan_test.mat');
copyfile(file_fullpath(srcFile), tmpFile);
cleanup = onCleanup(@() cleanupFcn(tmpFile));
T = load(tmpFile);
if isfield(T, 'TangentFrame')
    T = rmfield(T, 'TangentFrame');
    save(tmpFile, '-struct', 'T');
end

% --- Render ---
hFig = view_tangents(tmpFile, 'MaxArrows', 500);
assert(ishandle(hFig), 'view_tangents did not return a valid figure handle.');

% --- Auto-compute fired: frame now stored ---
Tafter = load(tmpFile, 'TangentFrame');
assert(isfield(Tafter, 'TangentFrame') && ~isempty(Tafter.TangentFrame), ...
    'Auto-compute did not store TangentFrame on the file.');

% --- Three quiver sets with equal counts; U/V share a color != normal ---
hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
hU = findobj(hAxes, 'Tag', 'tangentU');
hV = findobj(hAxes, 'Tag', 'tangentV');
hN = findobj(hAxes, 'Tag', 'tangentN');
assert(~isempty(hU) && ~isempty(hV) && ~isempty(hN), 'Missing tangentU/tangentV/tangentN quiver objects.');
nU = numel(hU.UData);  nV = numel(hV.UData);  nN = numel(hN.UData);
assert(nU == nV && nV == nN, 'U/V/N arrow counts differ.');
assert(nU <= 500, 'Drew more arrows (%d) than requested (500).', nU);
assert(isequal(hU.Color, hV.Color), 'U and V must share a color.');
assert(~isequal(hU.Color, hN.Color), 'Normal color must differ from U/V.');

% --- Frame arrows are equal length (a frame field has no magnitude) ---
lenU = sqrt(hU.UData(:).^2 + hU.VData(:).^2 + hU.WData(:).^2);
assert((max(lenU) - min(lenU)) < 1e-6, 'U arrows are not equal length.');
lenV = sqrt(hV.UData(:).^2 + hV.VData(:).^2 + hV.WData(:).^2);
assert((max(lenV) - min(lenV)) < 1e-6, 'V arrows are not equal length.');
lenN = sqrt(hN.UData(:).^2 + hN.VData(:).^2 + hN.WData(:).^2);
assert((max(lenN) - min(lenN)) < 1e-6, 'N arrows are not equal length.');

% --- Singularity markers match the stored pole count ---
hS = findobj(hAxes, 'Tag', 'tangentSing');
assert(~isempty(hS), 'No singularity markers drawn.');
assert(numel(hS.XData) == numel(Tafter.TangentFrame.Singularities.Vertices), ...
    'Singularity marker count does not match stored Singularities.Vertices.');

close(hFig);
fprintf('ALL TESTS PASSED: test_view_tangents\n');
end


function cleanupFcn(tmpFile)
    hViewer = findobj(0, 'Type', 'figure', '-regexp', 'Name', '^Tangent basis:');
    if ~isempty(hViewer)
        close(hViewer);
    end
    if exist(tmpFile, 'file')
        delete(tmpFile);
    end
end


function SurfaceFile = find_registered_cortex()
% Return a low-res cortex FileName that has Reg.Sphere.Vertices, or '' if none.
SurfaceFile = '';
best = inf;
sSubjects = bst_get('ProtocolSubjects');
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex'), continue; end
        try
            T = load(file_fullpath(surf(iF).FileName), 'Reg', 'Vertices');
        catch
            continue;
        end
        if isfield(T,'Reg') && isstruct(T.Reg) && isfield(T.Reg,'Sphere') ...
           && isfield(T.Reg.Sphere,'Vertices') && ~isempty(T.Reg.Sphere.Vertices)
            n = size(T.Vertices, 1);
            if n < best
                best = n;
                SurfaceFile = surf(iF).FileName;
            end
        end
    end
end
end
