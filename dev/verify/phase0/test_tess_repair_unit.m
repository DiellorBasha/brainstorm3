% TEST_TESS_REPAIR_UNIT: validation + repair on a synthetic closed mesh.
[V, F] = tess_sphere(642);
% 1) clean closed sphere validates
[~, ~, isM, rep] = tess_repair(V, F);
assert(isM, 'clean sphere must validate as manifold');
% 2) duplicated face (flipped winding) -> non-manifold edges detected
F2 = [F; F(1, [2 1 3])];
[~, ~, isM2, rep2] = tess_repair(V, F2);
assert(~isM2, 'corrupted mesh must fail validation');
% 3) repair removes the spurious face and restores manifoldness
[Vr, Fr, isM3, rep3] = tess_repair(V, F2, 'Repair', 1);
assert(isM3, 'repair must restore manifoldness');
assert(size(Fr,1) == size(F,1), 'repair should remove exactly the spurious face');
% 4) repairing an already-clean mesh is a no-op
[Vn, Fn, isM4] = tess_repair(V, F, 'Repair', 1);
assert(isM4 && isequal(Fn, F) && isequal(Vn, V), 'repair of clean mesh must be a no-op');
disp('test_tess_repair_unit PASSED');
