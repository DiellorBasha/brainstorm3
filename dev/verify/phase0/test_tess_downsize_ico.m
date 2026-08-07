% TEST_TESS_DOWNSIZE_ICO: icosphere downsampling on a synthetic registered sphere.
Protocol = 'IcoDownsizeUnit';
gui_brainstorm('DeleteProtocol', Protocol);
gui_brainstorm('CreateProtocol', Protocol, 0, 0);
[~, iSubject] = db_add_subject('TestSubj', [], 0, 0);
% Synthetic closed surface WITH registration sphere (the sphere is its own reg)
[V, F] = tess_sphere(40962);
TessMat = db_template('surfacemat');
TessMat.Comment  = 'sphere_40962';
TessMat.Vertices = V;
TessMat.Faces    = F;
TessMat.Reg.Sphere.Vertices = V;
ProtocolInfo = bst_get('ProtocolInfo');
TessFile = bst_fullfile(ProtocolInfo.SUBJECTS, 'TestSubj', 'tess_sphere_40962.mat');
bst_save(TessFile, TessMat, 'v7');
db_add_surface(iSubject, TessFile, TessMat.Comment);
% Downsample to the 10242 ico grid
[NewFile, iSurf, I, J] = tess_downsize(TessFile, 10242, 'icosphere');
NewMat = in_tess_bst(NewFile, 0);
assert(size(NewMat.Vertices,1) == 10242, 'expected snap to ico grid 10242');
assert(numel(unique(I)) == numel(I), 'vertex mapping must be injective');
assert(issorted(I), 'kept-vertex indices must be sorted (reducepatch convention)');
[~, ~, isM] = tess_repair(NewMat.Vertices, NewMat.Faces);
assert(isM, 'icosphere output must be a closed 2-manifold');
% Winding must match the source (signed volume same sign)
sv = @(Vv,Ff) sum(sum(Vv(Ff(:,1),:) .* cross(Vv(Ff(:,2),:), Vv(Ff(:,3),:), 2)));
assert(sign(sv(NewMat.Vertices,NewMat.Faces)) == sign(sv(V,F)), 'winding must match source');
gui_brainstorm('DeleteProtocol', Protocol);
disp('test_tess_downsize_ico PASSED');
