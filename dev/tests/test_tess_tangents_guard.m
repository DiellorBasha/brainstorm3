function test_tess_tangents_guard
% tess_tangents must error cleanly when its preconditions are not met:
%   (1) no FreeSurfer registration sphere (Reg.Sphere) -> noRegSphere
%   (2) no import-time hemisphere labels ('Structures' atlas) -> noHemisphereLabels
% Both paths must NOT require nxr, and tess_tangents must never re-split the
% mesh geometrically when the labels are absent.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% --- Case 1: WITHOUT Reg.Sphere -> tess_tangents:noRegSphere ---
% Name it tess_cortex_*.mat so file_fullpath recognizes the type and returns
% the absolute path as-is.
[V, F] = tess_sphere(162);
tmpFile = fullfile(tempdir, 'tess_cortex_tangents_guard.mat');
bst_save(tmpFile, struct('Vertices', V, 'Faces', F, 'Comment', 'noreg'), 'v7');
cleanup = onCleanup(@() delete(tmpFile));

threw = false;
try
    tess_tangents(tmpFile, 'NoSave', 1);
catch ME
    threw = strcmp(ME.identifier, 'tess_tangents:noRegSphere');
    if ~threw
        rethrow(ME);
    end
end
assert(threw, 'tess_tangents did not raise tess_tangents:noRegSphere on an unregistered surface.');
fprintf('PASSED: errors with noRegSphere when the registration sphere is absent.\n');

% --- Case 2: has Reg.Sphere but NO 'Structures' atlas -> noHemisphereLabels ---
% (Confirms tess_tangents refuses to re-split geometrically when labels are absent.)
reg = V ./ max(sqrt(sum(V.^2, 2)), eps);   % a unit registration sphere
tmpFile2 = fullfile(tempdir, 'tess_cortex_tangents_nolabels.mat');
bst_save(tmpFile2, struct('Vertices', V, 'Faces', F, ...
    'Reg', struct('Sphere', struct('Vertices', reg)), ...
    'Atlas', db_template('Atlas'), 'iAtlas', 1, 'Comment', 'nolabels'), 'v7');
cleanup2 = onCleanup(@() delete(tmpFile2));

threw2 = false;
try
    tess_tangents(tmpFile2, 'NoSave', 1);
catch ME2
    threw2 = strcmp(ME2.identifier, 'tess_tangents:noHemisphereLabels');
    if ~threw2
        rethrow(ME2);
    end
end
assert(threw2, 'tess_tangents did not raise tess_tangents:noHemisphereLabels when the Structures atlas is absent.');
fprintf('PASSED: errors with noHemisphereLabels when hemisphere labels are absent.\n');

fprintf('ALL TESTS PASSED: test_tess_tangents_guard\n');
end
