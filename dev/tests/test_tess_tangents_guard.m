function test_tess_tangents_guard
% tess_tangents must error cleanly when the surface has no FreeSurfer
% registration sphere (Reg.Sphere). This path must NOT require nxr.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% Build a temp surface WITHOUT Reg.Sphere. Name it tess_cortex_*.mat so
% file_fullpath recognizes the type and returns the absolute path as-is.
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
fprintf('ALL TESTS PASSED: test_tess_tangents_guard\n');
end
