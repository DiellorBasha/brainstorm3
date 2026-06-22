function test_in_bst_nodes_pure()
% Pure round-trip test for the new derived-anatomy loaders:
%   in_bst_eigen / in_bst_operator / in_bst_manifold
% Builds each node from its db_template, saves to a temp .mat (absolute path, so no
% protocol/cache is needed), and verifies: full read, field-subset read, missing-field
% backfill, and the type-signature guard. Requires Brainstorm on the path.

% Build properly-prefixed temp filenames (file_gettype keys on the eigen_/operator_/
% manifold_ prefix; an unprefixed name would be rejected by file_fullpath).
[~, tag] = fileparts(tempname);
fEig = fullfile(tempdir, ['eigen_'    tag '.mat']);
fOp  = fullfile(tempdir, ['operator_' tag '.mat']);
fMan = fullfile(tempdir, ['manifold_' tag '.mat']);
fBad = fullfile(tempdir, ['eigen_'    tag '_bad.mat']);
cleanup = onCleanup(@() delete_if_exist({fEig, fOp, fMan, fBad}));

%% ===== EIGEN =====
E = db_template('eigenmat');
E.Variant = 'Dirac';  E.K = 5;
E.Phi    = {rand(20,5), rand(20,5)};
E.Lambda = {sort(rand(5,1)), sort(rand(5,1))};
E.GlobalVertices = {(1:5)', (6:10)'};
bst_save(fEig, E, 'v7');

Efull = in_bst_eigen(fEig);
assert(isfield(Efull,'Phi') && isequal(size(Efull.Phi),[1 2]), 'eigen full read: Phi');
assert(strcmp(Efull.Variant,'Dirac') && isequal(Efull.K,5),    'eigen full read: Variant/K');

Esub = in_bst_eigen(fEig, 'Variant', 'K', 'NotAField');
assert(strcmp(Esub.Variant,'Dirac'),       'eigen subset: Variant present');
assert(isequal(Esub.K,5),                  'eigen subset: K present');
assert(isfield(Esub,'NotAField') && isempty(Esub.NotAField), 'eigen subset: missing field backfilled to []');
assert(~isfield(Esub,'Phi'),               'eigen subset: unrequested field absent');

%% ===== OPERATOR =====
O = db_template('operatormat');
O.Variant = 'Dirac';
O.Operator = {sprandsym(20,0.2), sprandsym(20,0.2)};
O.Mass     = {speye(20), speye(20)};
O.GlobalVertices = {(1:10)', (11:20)'};
bst_save(fOp, O, 'v7');

Ofull = in_bst_operator(fOp);
assert(isfield(Ofull,'Operator') && isfield(Ofull,'Mass'), 'operator full read: Operator/Mass');

Osub = in_bst_operator(fOp, 'Mass');
assert(isfield(Osub,'Mass') && ~isfield(Osub,'Operator'),  'operator subset: Mass only');

%% ===== MANIFOLD =====
M = db_template('manifoldmat');
M.Topology = struct('GlobalVertices', {(1:10)', (11:20)'});
bst_save(fMan, M, 'v7');

Mfull = in_bst_manifold(fMan);
assert(isfield(Mfull,'Topology'), 'manifold full read: Topology');

%% ===== TYPE-SIGNATURE GUARD =====
bad = struct('SomethingElse', 42);
bst_save(fBad, bad, 'v7');
threw = false;
try
    in_bst_eigen(fBad);
catch
    threw = true;
end
assert(threw, 'type guard: in_bst_eigen must reject a non-eigen file');

%% ===== MISSING FILE =====
threw2 = false;
try
    in_bst_eigen(fullfile(tempdir, ['eigen_' tag '_does_not_exist.mat']));
catch
    threw2 = true;
end
assert(threw2, 'missing file: in_bst_eigen must error');

fprintf('ALL TESTS PASSED: test_in_bst_nodes_pure\n');
end

% ----------------------------------------------------------------------------
function delete_if_exist(files)
% Remove any of the listed temp files that exist (cleanup helper).
    for i = 1:numel(files)
        if exist(files{i}, 'file') == 2
            delete(files{i});
        end
    end
end
