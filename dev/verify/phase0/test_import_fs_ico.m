% TEST_IMPORT_FS_ICO: non-interactive FS import must yield an ico5 manifold cortex.
FsDir = '/Volumes/SpikeData-2/workspace/library/datasets/omega-tutorial/derivatives/freesurfer/sub-0002/ses-mri/anat';
Protocol = 'FsIcoUnit';
gui_brainstorm('DeleteProtocol', Protocol);
gui_brainstorm('CreateProtocol', Protocol, 0, 0);
[~, iSubject] = db_add_subject('FsSubj', [], 0, 0);
errorMsg = import_anatomy_fs(iSubject, FsDir, [], 0, [], 0, 1, 0, 'icosphere');
assert(isempty(errorMsg), 'import_anatomy_fs error: %s', errorMsg);
sSubject = bst_get('Subject', iSubject);
CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
Cortex = in_tess_bst(CortexFile, 0);
assert(size(Cortex.Vertices,1) == 20484, 'expected ico5 cortex (20484), got %d', size(Cortex.Vertices,1));
% Per-hemisphere manifold check: split on the Structures atlas, never conncomp
[rH, lH] = tess_hemisplit(in_tess_bst(CortexFile));
for hemi = {rH, lH}
    iV = hemi{1};
    iF = all(ismember(Cortex.Faces, iV), 2);
    lut = zeros(size(Cortex.Vertices,1),1); lut(iV) = 1:numel(iV);
    [~, ~, isM] = tess_repair(Cortex.Vertices(iV,:), lut(Cortex.Faces(iF,:)));
    assert(isM, 'hemisphere fails 2-manifold check');
end
gui_brainstorm('DeleteProtocol', Protocol);
disp('test_import_fs_ico PASSED');
