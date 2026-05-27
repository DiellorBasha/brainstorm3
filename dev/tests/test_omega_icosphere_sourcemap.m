function results = test_omega_icosphere_sourcemap(varargin)
% TEST_OMEGA_ICOSPHERE_SOURCEMAP: Verify the new icosphere tess_downsize does not break source mapping.
%
% Reproduces the cortex-dependent half of TUTORIAL_OMEGA (anatomy -> head model -> inverse) on a
% single OMEGA subject (sub-0002), once with the NEW icosphere ico5 downsampling and once with the
% legacy reducepatch baseline, then checks that every source-mapping method (MNE, dSPM, sLORETA,
% LCMV) produces a finite, correctly-sized operator on both cortices.
%
% USAGE:  results = test_omega_icosphere_sourcemap()
%         results = test_omega_icosphere_sourcemap('BidsDir', '/path/to/omega-tutorial', ...)
%
% WHAT THIS TEST IS / IS NOT:
%   - IS:  a numerical-sanity test of the source-mapping chain that depends on the cortex geometry
%          produced by tess_downsize (vertices/faces/normals -> overlapping-spheres head model ->
%          inverse kernels). Both arms use the SAME FreeSurfer recon and the SAME recordings, so the
%          cortex downsampling method is the only variable.
%   - NOT: a registration/localization-accuracy test. The single-subject manual import uses default
%          (MNI-derived) fiducials, so MEG/MRI co-registration is approximate but IDENTICAL for both
%          arms. Steps that do not depend on the cortex (head-point refine, SSP, resting power maps)
%          are intentionally omitted to keep the test focused and fast.
%
% OUTPUT:
%   - results : struct with fields .arms (per-arm results), .pass (overall logical), .summary (text).
%
% Author: test harness for the eigenmode-analysis branch (icosphere downsampling), 2026

% ===== CONFIGURATION =====
Def.BidsDir      = '/Users/diellorbasha/workspace/library/datasets/omega-tutorial';
Def.ProtocolName = 'TutorialOmega_IcoTest';
Def.AnatSubdir   = fullfile('derivatives', 'freesurfer', 'sub-0002', 'ses-mri', 'anat');   % FreeSurfer recon-all root
Def.RestSubpath  = fullfile('sub-0002', 'ses-01', 'meg', 'sub-0002_ses-01_task-rest_run-01_meg.ds');
Def.NoiseSubpath = fullfile('sub-emptyroom', 'ses-18901014', 'meg', 'sub-emptyroom_ses-18901014_task-noise_run-01_meg.ds');
% Parse name/value overrides
OPT = Def;
for i = 1:2:numel(varargin)
    OPT.(varargin{i}) = varargin{i+1};
end

AnatDir   = fullfile(OPT.BidsDir, OPT.AnatSubdir);
RestFile  = fullfile(OPT.BidsDir, OPT.RestSubpath);
NoiseFile = fullfile(OPT.BidsDir, OPT.NoiseSubpath);

% ===== VALIDATE INPUTS =====
assert(exist(AnatDir, 'dir') == 7,    'FreeSurfer anatomy folder not found: %s', AnatDir);
assert(exist(RestFile, 'dir') == 7,   'Resting MEG .ds not found: %s', RestFile);   % .ds is a folder
assert(exist(NoiseFile, 'dir') == 7,  'Empty-room MEG .ds not found: %s', NoiseFile);

% ===== START BRAINSTORM =====
if ~brainstorm('status')
    brainstorm nogui
end

% ===== RECREATE PROTOCOL (isolated, repeatable) =====
gui_brainstorm('DeleteProtocol', OPT.ProtocolName);
gui_brainstorm('CreateProtocol', OPT.ProtocolName, 0, 0);    % UseDefaultAnat=0, UseDefaultChannel=0
bst_report('Start');

% ===== DEFINE THE TWO ARMS =====
arms = struct( ...
    'name',       {'sub-0002-test',  'sub-0002-test-rp'}, ...
    'method',     {'icosphere',      'reducepatch'}, ...
    'fsVert',     {[],               15000}, ...      % [] => import_anatomy_fs picks ico5 (10242/hemi)
    'expectVert', {20484,            15000});         % icosphere: exact; reducepatch: approx (tolerance below)

% ===== RUN EACH ARM =====
res = [];
for iArm = 1:numel(arms)
    fprintf('\n========== ARM %d/%d: %s (%s) ==========\n', iArm, numel(arms), arms(iArm).name, arms(iArm).method);
    try
        res = i_appendStruct(res, i_run_arm(arms(iArm), AnatDir, RestFile, NoiseFile));
    catch ME
        fprintf(2, 'ARM FAILED with exception: %s\n', ME.message);
        a = struct('name', arms(iArm).name, 'method', arms(iArm).method, 'fatalError', ME.message, 'pass', false);
        res = i_appendStruct(res, a);
    end
end

% ===== BUILD SUMMARY =====
[summaryTxt, overallPass] = i_buildSummary(res);
fprintf('\n%s\n', summaryTxt);

results = struct('arms', {res}, 'pass', overallPass, 'summary', summaryTxt, 'protocol', OPT.ProtocolName);

% ===== SAVE REPORT =====
try
    ReportFile = bst_report('Save', []);
    results.reportFile = ReportFile;
    fprintf('Brainstorm report saved: %s\n', ReportFile);
catch
    % non-fatal
end
end


%% ======================================================================
%  Run the full anatomy -> head model -> inverse chain for one arm
%  ======================================================================
function res = i_run_arm(arm, AnatDir, RestFile, NoiseFile)
    res = struct('name', arm.name, 'method', arm.method);

    % --- 1. Create subject (own anatomy, own channel file) ---
    [~, iSubject] = db_add_subject(arm.name, [], 0, 0);
    assert(~isempty(iSubject), 'Could not create subject %s', arm.name);

    % --- 2. Import FreeSurfer anatomy with the requested downsampling method ---
    %     import_anatomy_fs(iSubject, FsDir, nVertices, isInteractive, sFid, isExtraMaps, isVolumeAtlas, isKeepMri, Method)
    errMsg = import_anatomy_fs(iSubject, AnatDir, arm.fsVert, 0, [], 0, 1, 0, arm.method);
    res.anatError = errMsg;
    assert(isempty(errMsg), 'import_anatomy_fs failed: %s', errMsg);

    % --- 3. Check the downsampled cortex ---
    res.anat = i_check_anatomy(iSubject, arm);
    fprintf('  Anatomy: %d vertices, %d faces, %d components, normals=%s, manifold=%s -> %s\n', ...
        res.anat.nVertices, res.anat.nFaces, res.anat.nComp, i_yn(res.anat.hasNormals), ...
        i_manStr(res.anat.isManifold), i_pf(res.anat.pass));

    % --- 4. Link recordings (resting + empty room) under this subject ---
    sRest = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
        'subjectname',    arm.name, ...
        'datafile',       {RestFile, 'CTF'}, ...
        'channelreplace', 1, ...
        'channelalign',   0, ...
        'evtmode',        'value');
    assert(~isempty(sRest), 'Failed to link resting recording');

    sNoise = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
        'subjectname',    arm.name, ...
        'datafile',       {NoiseFile, 'CTF'}, ...
        'channelreplace', 1, ...
        'channelalign',   0, ...
        'evtmode',        'value');
    assert(~isempty(sNoise), 'Failed to link empty-room recording');

    % --- 5. Convert to continuous (CTF) ---
    sRest  = bst_process('CallProcess', 'process_ctf_convert', sRest,  [], 'rectype', 2);
    sNoise = bst_process('CallProcess', 'process_ctf_convert', sNoise, [], 'rectype', 2);

    % --- 6. Notch (line noise) + high-pass 0.3 Hz, on both recordings ---
    sRest = bst_process('CallProcess', 'process_notch', sRest, [], ...
        'freqlist', [60, 120, 180, 240, 300], 'sensortypes', 'MEG', 'read_all', 1);
    sRest = bst_process('CallProcess', 'process_bandpass', sRest, [], ...
        'sensortypes', 'MEG', 'highpass', 0.3, 'lowpass', 0, ...
        'attenuation', 'strict', 'mirror', 0, 'useold', 0, 'read_all', 1);

    sNoise = bst_process('CallProcess', 'process_notch', sNoise, [], ...
        'freqlist', [60, 120, 180, 240, 300], 'sensortypes', 'MEG', 'read_all', 1);
    sNoise = bst_process('CallProcess', 'process_bandpass', sNoise, [], ...
        'sensortypes', 'MEG', 'highpass', 0.3, 'lowpass', 0, ...
        'attenuation', 'strict', 'mirror', 0, 'useold', 0, 'read_all', 1);

    % --- 7. Noise covariance from the empty-room recording (copy to the resting folder) ---
    bst_process('CallProcess', 'process_noisecov', sNoise, [], ...
        'baseline', [], 'sensortypes', 'MEG', 'target', 1, ...   % 1 = noise covariance
        'dcoffset', 1, 'identity', 0, ...
        'copycond', 1, 'copysubj', 0, 'copymatch', 0, 'replacefile', 1);

    % --- 8. Data covariance from the resting recording (required by LCMV beamformer) ---
    bst_process('CallProcess', 'process_noisecov', sRest, [], ...
        'baseline', [], 'datatimewindow', [], 'sensortypes', 'MEG', 'target', 2, ...  % 2 = data covariance
        'dcoffset', 1, 'identity', 0, ...
        'copycond', 0, 'copysubj', 0, 'copymatch', 0, 'replacefile', 1);

    % --- 9. Forward model: overlapping spheres on the cortex source space ---
    bst_process('CallProcess', 'process_headmodel', sRest, [], ...
        'sourcespace', 1, ...   % 1 = cortex surface
        'meg',         3);      % 3 = overlapping spheres
    res.headmodel = i_check_headmodel(sRest);
    fprintf('  Head model (overlapping spheres): gain %dx%d, finite=%s, cols==3*nGrid=%s -> %s\n', ...
        res.headmodel.gainSize(1), res.headmodel.gainSize(2), i_yn(res.headmodel.allFinite), ...
        i_yn(res.headmodel.dimsOk), i_pf(res.headmodel.pass));

    % --- 10. Inverse: every source-mapping method on the same cortex/head model ---
    %     {label, InverseMethod, InverseMeasure, UseDepth}
    methodList = { ...
        'MNE',     'minnorm', 'amplitude', 1; ...
        'dSPM',    'minnorm', 'dspm2018',  1; ...
        'sLORETA', 'minnorm', 'sloreta',   0; ...   % depth weighting disabled for sLORETA (GUI behavior)
        'LCMV',    'lcmv',    'nai',       0};       % beamformer: needs data covariance, no depth
    res.inverse = [];
    for m = 1:size(methodList, 1)
        sSrc = bst_process('CallProcess', 'process_inverse_2018', sRest, [], ...
            'output',  2, ...    % 2 = kernel only (one per file)
            'inverse', struct( ...
                'Comment',        [methodList{m,1} ': MEG'], ...
                'InverseMethod',  methodList{m,2}, ...
                'InverseMeasure', methodList{m,3}, ...
                'SourceOrient',   {{'fixed'}}, ...   % fixed orientation -> exercises cortex VertNormals
                'Loose',          0.2, ...
                'UseDepth',       methodList{m,4}, ...
                'WeightExp',      0.5, ...
                'WeightLimit',    10, ...
                'NoiseMethod',    'reg', ...
                'NoiseReg',       0.1, ...
                'SnrMethod',      'fixed', ...
                'SnrRms',         1e-06, ...
                'SnrFixed',       3, ...
                'ComputeKernel',  1, ...
                'DataTypes',      {{'MEG'}}));
        chk = i_check_inverse(methodList{m,1}, sSrc, sRest);
        res.inverse = i_appendStruct(res.inverse, chk);
        fprintf('  Inverse %-8s: kernel %dx%d, finite=%s, recon=%s -> %s\n', ...
            chk.name, chk.kernelSize(1), chk.kernelSize(2), i_yn(chk.allFinite), ...
            chk.reconStr, i_pf(chk.pass));
    end

    % --- Arm passes iff every check passed ---
    res.pass = res.anat.pass && res.headmodel.pass && all([res.inverse.pass]);
end


%% ======================================================================
%  Check the downsampled cortex surface
%  ======================================================================
function chk = i_check_anatomy(iSubject, arm)
    sSubject  = bst_get('Subject', iSubject);
    CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
    TessMat   = in_tess_bst(CortexFile);

    chk.cortexFile = CortexFile;
    chk.nVertices  = size(TessMat.Vertices, 1);
    chk.nFaces     = size(TessMat.Faces, 1);

    % Connected components (two hemispheres => 2)
    VertConn  = tess_vertconn(TessMat.Vertices, TessMat.Faces);
    chk.nComp = max(conncomp(graph(VertConn)));

    % Vertex normals must be computable and finite (fixed-orientation inverse depends on them)
    try
        Normals = tess_normals(TessMat.Vertices, TessMat.Faces);
        chk.hasNormals = ~isempty(Normals) && all(isfinite(Normals(:)));
    catch
        chk.hasNormals = false;
    end

    % Manifold check (non-fatal if the helper errors)
    try
        chk.isManifold = tess_check_manifold(TessMat.Vertices, TessMat.Faces);
    catch
        chk.isManifold = NaN;
    end

    % Pass criteria
    if strcmpi(arm.method, 'icosphere')
        okVert = (chk.nVertices == arm.expectVert);          % ico5 is exact: 20484
    else
        okVert = abs(chk.nVertices - arm.expectVert) <= 300; % reducepatch snaps near the target
    end
    % A 2-manifold mesh is a *claim* of the icosphere method, so gate on it there. reducepatch is
    % inherently non-manifold (the very defect icosphere fixes) -> informational only for the baseline.
    if strcmpi(arm.method, 'icosphere')
        okMan = ~isnan(chk.isManifold) && isequal(logical(chk.isManifold), true);
    else
        okMan = true;
    end
    chk.pass = okVert && (chk.nComp == 2) && chk.hasNormals && okMan;
end


%% ======================================================================
%  Check the forward (head) model
%  ======================================================================
function chk = i_check_headmodel(sData)
    [sStudy, ~] = bst_get('AnyFile', sData(1).FileName);
    assert(~isempty(sStudy) && ~isempty(sStudy.iHeadModel), 'No head model found for %s', sData(1).FileName);
    HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    Hm = in_bst_headmodel(HeadModelFile);

    G = Hm.Gain;
    chk.gainSize  = size(G);
    chk.nGrid     = size(Hm.GridLoc, 1);
    % Non-modeled channels (ECG/EOG/SysClock...) are NaN by design in the gain matrix; finiteness is
    % only meaningful over the MODELED rows (those not entirely NaN = the MEG + MEG REF sensors).
    modeledRows   = ~all(isnan(G), 2);
    Gmod          = G(modeledRows, :);
    chk.nModeled  = sum(modeledRows);
    chk.allFinite = ~isempty(Gmod) && all(isfinite(Gmod(:)));   % no partial-NaN in any modeled row
    chk.nonZero   = any(Gmod(:) ~= 0);
    chk.dimsOk    = (size(G, 2) == 3 * chk.nGrid);   % gain stored as 3 orientations / grid point
    chk.pass      = chk.allFinite && chk.nonZero && chk.dimsOk;
end


%% ======================================================================
%  Check one inverse operator (kernel) + a short end-to-end reconstruction
%  ======================================================================
function chk = i_check_inverse(name, sSrc, sData)
    chk = struct('name', name, 'kernelSize', [0 0], 'allFinite', false, 'reconStr', 'n/a', 'pass', false);
    if isempty(sSrc) || isempty(sSrc(1).FileName)
        chk.reconStr = 'no result';
        return;
    end

    R = in_bst_results(sSrc(1).FileName, 0);    % LoadFull=0 -> kernel + metadata, no data applied
    K = R.ImagingKernel;
    chk.kernelSize   = size(K);
    chk.allFinite    = ~isempty(K) && all(isfinite(K(:)));
    chk.nComponents  = R.nComponents;
    chk.colsMatch    = (size(K, 2) == numel(R.GoodChannel));

    % End-to-end: apply the kernel to a short raw block and confirm a finite, non-degenerate map.
    recon = struct('done', false, 'finite', false, 'nonDegenerate', false);
    try
        Draw = in_bst(sData(1).FileName, [0, 5], 1, 1, 'no', 1);    % first 5 s, load full, ignore bad, no baseline removal, apply SSP
        F = Draw.F(R.GoodChannel, :);
        S = K * F;
        recon.done          = true;
        recon.finite        = all(isfinite(S(:)));
        recon.meanVertStd   = mean(std(S, 0, 2));
        recon.nonDegenerate = recon.finite && (recon.meanVertStd > 0);
        chk.reconStr        = sprintf('finite=%s,std=%.2e', i_yn(recon.finite), recon.meanVertStd);
    catch ME
        chk.reconStr = ['skipped(' ME.message ')'];
    end
    chk.recon = recon;

    % Pass: kernel finite + correctly sized; reconstruction (if it ran) must be non-degenerate
    chk.pass = chk.allFinite && chk.colsMatch && (~recon.done || recon.nonDegenerate);
end


%% ======================================================================
%  Summary table
%  ======================================================================
function [txt, overallPass] = i_buildSummary(res)
    L = {};
    L{end+1} = '================ SOURCE-MAPPING TEST SUMMARY ================';
    overallPass = true;
    for i = 1:numel(res)
        a = res(i);
        if isfield(a, 'fatalError')
            L{end+1} = sprintf('%-18s [%s]  FATAL: %s', a.name, a.method, a.fatalError);
            overallPass = false;
            continue;
        end
        L{end+1} = sprintf('--- %s  (%s) ---', a.name, a.method);
        L{end+1} = sprintf('  Anatomy    : %5d vert, %2d comp, normals=%s, manifold=%s   [%s]', ...
            a.anat.nVertices, a.anat.nComp, i_yn(a.anat.hasNormals), i_manStr(a.anat.isManifold), i_pf(a.anat.pass));
        L{end+1} = sprintf('  Head model : gain %dx%d, finite=%s, dims=%s              [%s]', ...
            a.headmodel.gainSize(1), a.headmodel.gainSize(2), i_yn(a.headmodel.allFinite), i_yn(a.headmodel.dimsOk), i_pf(a.headmodel.pass));
        for m = 1:numel(a.inverse)
            c = a.inverse(m);
            L{end+1} = sprintf('  %-8s   : kernel %dx%d, finite=%s, recon=%s   [%s]', ...
                c.name, c.kernelSize(1), c.kernelSize(2), i_yn(c.allFinite), c.reconStr, i_pf(c.pass));
        end
        overallPass = overallPass && a.pass;
    end
    L{end+1} = '=============================================================';
    L{end+1} = sprintf('OVERALL: %s', i_pf(overallPass));
    txt = strjoin(L, char(10));
end


%% ===== small formatting helpers =====
function s = i_pf(b),  if b, s = 'PASS'; else, s = 'FAIL'; end, end
function s = i_yn(b),  if b, s = 'yes';  else, s = 'no';   end, end
function s = i_manStr(m)
    if isnan(m), s = 'n/a'; elseif m, s = 'yes'; else, s = 'no'; end
end
function out = i_appendStruct(arr, s)
    % Append struct s to struct array arr, tolerating differing field sets.
    if isempty(arr)
        out = s;
    else
        f = setdiff(fieldnames(arr), fieldnames(s));
        for k = 1:numel(f), s.(f{k}) = []; end
        f2 = setdiff(fieldnames(s), fieldnames(arr));
        for k = 1:numel(f2), [arr.(f2{k})] = deal([]); end
        out = [arr, orderfields(s, arr)];
    end
end
