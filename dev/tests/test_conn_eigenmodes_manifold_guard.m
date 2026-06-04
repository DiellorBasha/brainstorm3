function test_conn_eigenmodes_manifold_guard
% tess_conn_eigenmodes must reject a non-manifold mesh with a clean,
% fail-fast error (tess_conn_eigenmodes:NonManifold) BEFORE assembling the
% operator -- so the gate fires without needing the nxr-compute plugin.
% The mesh below has one edge (1-2) shared by three faces (multiplicity 3),
% which tess_manifold flags as a non-manifold edge.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% Non-manifold mesh: edge (1,2) is shared by 3 triangles.
V = [0 0 0; 1 0 0; 0 1 0; 0 -1 0; 1 1 1];
F = [1 2 3; 1 2 4; 1 2 5];

threw = false;
try
    tess_conn_eigenmodes(V, F, 'nModes', 2, 'Verbose', 0);
catch ME
    threw = strcmp(ME.identifier, 'tess_conn_eigenmodes:NonManifold');
    if ~threw
        rethrow(ME);
    end
end
assert(threw, 'tess_conn_eigenmodes did not raise tess_conn_eigenmodes:NonManifold on a non-manifold mesh.');

fprintf('PASSED: tess_conn_eigenmodes rejects a non-manifold mesh (NonManifold, fail-fast).\n');
fprintf('ALL TESTS PASSED: test_conn_eigenmodes_manifold_guard\n');
end
