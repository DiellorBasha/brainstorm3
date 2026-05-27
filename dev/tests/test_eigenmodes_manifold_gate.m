function test_eigenmodes_manifold_gate
% Verify tess_eigenmodes validates manifoldness and never auto-repairs.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% ---- Clean manifold mesh ----
[V, F] = tess_sphere(642);
nV0 = size(V, 1);

% Case 1: clean mesh computes; vertex count preserved.
[Em, ~, ~, Vout] = tess_eigenmodes(V, F, 'nModes', 20, 'Verbose', 0);
assert(Em.nModes == 20, 'Expected 20 modes, got %d', Em.nModes);
assert(size(Vout,1) == nV0, 'Clean mesh vertex count changed (%d -> %d)', nV0, size(Vout,1));
fprintf('PASSED: clean manifold mesh computes, vertices preserved.\n');

% ---- Non-manifold mesh: extra face creates an edge shared by 3 faces ----
vExtra = find(~ismember(1:nV0, F(1,:)), 1);
Fnm = [F; F(1,1), F(1,2), vExtra];

% Case 2: non-manifold errors by DEFAULT.
threw = false;
try
    tess_eigenmodes(V, Fnm, 'nModes', 20, 'Verbose', 0);
catch ME
    threw = true;
    assert(~isempty(regexpi(ME.message, 'manifold')), ...
        'Expected a manifold error, got: %s', ME.message);
end
assert(threw, 'Non-manifold input did NOT error by default.');
fprintf('PASSED: non-manifold errors by default.\n');

% Case 3: explicit Repair=true succeeds.
Em3 = tess_eigenmodes(V, Fnm, 'nModes', 10, 'Repair', true, 'Verbose', 0);
assert(Em3.nModes == 10, 'Repair path: expected 10 modes, got %d', Em3.nModes);
fprintf('PASSED: Repair=true repairs and computes.\n');

fprintf('ALL TESTS PASSED: test_eigenmodes_manifold_gate\n');
end
