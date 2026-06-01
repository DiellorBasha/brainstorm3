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

% --- Surface geometry, for the per-face fit check below ---
Tgeo = load(tmpFile, 'Vertices', 'Faces');
Vtx = Tgeo.Vertices;  Fcs = double(Tgeo.Faces);  nF = size(Fcs, 1);
a = sqrt(sum((Vtx(Fcs(:,2),:) - Vtx(Fcs(:,1),:)).^2, 2));
b = sqrt(sum((Vtx(Fcs(:,3),:) - Vtx(Fcs(:,2),:)).^2, 2));
c = sqrt(sum((Vtx(Fcs(:,1),:) - Vtx(Fcs(:,3),:)).^2, 2));
sP = (a + b + c) / 2;
inRadAll = sqrt(max(sP .* (sP-a) .* (sP-b) .* (sP-c), 0)) ./ max(sP, eps);

% --- U & V present; normals OFF by default; U/V share a color ---
hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
hU = findobj(hAxes, 'Tag', 'tangentU');
hV = findobj(hAxes, 'Tag', 'tangentV');
assert(~isempty(hU) && ~isempty(hV), 'Missing tangentU/tangentV objects.');
assert(isempty(findobj(hAxes, 'Tag', 'tangentN')), 'Normal glyphs should be OFF by default.');
nU = numel(hU.UData);  nV = numel(hV.UData);
assert(nU == nV, 'U/V glyph counts differ.');
assert(nU <= 500, 'Drew more glyphs (%d) than requested (500).', nU);
assert(isequal(hU.Color, hV.Color), 'U and V must share a color.');

% --- Glyphs are headless lines (axis-like), not arrows ---
assert(strcmpi(hU.ShowArrowHead, 'off'), 'U glyphs should have no arrowheads.');
assert(strcmpi(hV.ShowArrowHead, 'off'), 'V glyphs should have no arrowheads.');

% --- Glyphs sized per-face to fit inside the triangle (inscribed radius) ---
lenU = sqrt(hU.UData(:).^2 + hU.VData(:).^2 + hU.WData(:).^2);
assert(max(lenU) <= max(inRadAll) + 1e-9, 'U glyphs exceed per-face inscribed-circle size.');
assert((max(lenU) - min(lenU)) > 1e-9, 'U glyphs are not per-face normalized (all equal).');
lenV = sqrt(hV.UData(:).^2 + hV.VData(:).^2 + hV.WData(:).^2);
assert(max(lenV) <= max(inRadAll) + 1e-9, 'V glyphs exceed per-face inscribed-circle size.');

% --- Singularity markers match the stored pole count ---
hS = findobj(hAxes, 'Tag', 'tangentSing');
assert(~isempty(hS), 'No singularity markers drawn.');
assert(numel(hS.XData) == numel(Tafter.TangentFrame.Singularities.Vertices), ...
    'Singularity marker count does not match stored Singularities.Vertices.');
% Markers are blue lollipops with a stem to the true pole
assert(hS.MarkerFaceColor(3) > hS.MarkerFaceColor(1) && hS.MarkerFaceColor(3) > hS.MarkerFaceColor(2), ...
    'Singularity markers should be blue.');
hStem = findobj(hAxes, 'Tag', 'tangentSingStem');
assert(~isempty(hStem), 'No lollipop stems drawn for the singularities.');

% --- Legend explains the singularity marker ---
hLeg = findobj(hFig, 'Tag', 'tangentLegend');
assert(~isempty(hLeg), 'No legend created.');
assert(any(strcmp(hLeg.String, 'FreeSurfer sphere pole')), ...
    'Legend missing the "FreeSurfer sphere pole" entry.');

% --- Interaction: right arrow increases density ---
kp = get(hFig, 'KeyPressFcn');
evMore = struct('Key', 'rightarrow', 'Modifier', {{}});
kp(hFig, evMore);
nUmore = numel(findobj(hAxes, 'Tag', 'tangentU').UData);
assert(nUmore > nU, 'Right arrow did not increase glyph density.');

% --- Interaction: N toggles the normal glyphs ON (off by default) ---
evN = struct('Key', 'n', 'Modifier', {{}});
kp(hFig, evN);
hN = findobj(hAxes, 'Tag', 'tangentN');
assert(~isempty(hN), 'N did not toggle the normal glyphs on.');
assert(numel(hN.UData) == nUmore, 'Normal glyph count differs from tangent count.');
hUnow = findobj(hAxes, 'Tag', 'tangentU');   % re-fetch (earlier handle was redrawn)
assert(~isequal(hN.Color, hUnow.Color), 'Normal color must differ from U/V.');

% --- Interaction: P toggles the singularity markers off ---
evP = struct('Key', 'p', 'Modifier', {{}});
kp(hFig, evP);
assert(isempty(findobj(hAxes, 'Tag', 'tangentSing')), 'P did not toggle the singularity markers off.');
assert(isempty(findobj(hAxes, 'Tag', 'tangentSingStem')), 'P did not remove the lollipop stems.');

close(hFig);

% --- Default density: ~half the faces when MaxArrows is not given ---
hFig2 = view_tangents(tmpFile);
hAxes2 = findobj(hFig2, '-depth', 1, 'Tag', 'Axes3D');
nDef = numel(findobj(hAxes2, 'Tag', 'tangentU').UData);
assert(nDef >= 0.45 * nF && nDef <= 0.5 * nF + 1, ...
    'Default density (%d) is not about half of %d faces.', nDef, nF);
close(hFig2);

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
