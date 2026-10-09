function ok = test_bst_pet_tracer(AnatDir, OutDir)
% TEST_BST_PET_TRACER  Self-check of the nsp bst-pet per-tracer step (bst_pet_tracer) on real data.
%
% Builds a throwaway protocol from one subject of an existing PREVENT-AD Brainstorm protocol (T1,
% ASEG, Desikan-Killiany, low-resolution cortex, and the stored "PET <trc>_mean" volume, which still
% carries its BIDS PET metadata), then runs bst_pet_tracer with 'gtm' (scanner PSF from metadata)
% and 'vlpp', and checks the pr/pet-pvc-2 API path: HRRT PSF 2.5 mm, GTM regional table, robust
% reference, provenance, cortical map. Also checks bst_git_commit against git on this checkout.
% Prints PASS/FAIL lines; ok = all passed.
%
% ⚠ Brainstorm writes <user.home>/.brainstorm: run it with a throwaway home, e.g.
%   T=$(mktemp -d); JAVA_TOOL_OPTIONS="-Duser.home=$T" matlab -nodisplay -batch ...
%     "addpath('<bst>','<bst>/dev'); test_bst_pet_tracer('<protocol>/anat/sub-MTL0002', '$T/out')"
%
% Authors: Diellor Basha, 2026 (nsp brainstorm-pet pathway)

    bstDir = fileparts(fileparts(mfilename('fullpath')));
    home = char(java.lang.System.getProperty('user.home'));
    assert(~strcmp(home, getenv('HOME')), ...
        'Run with a throwaway user.home (JAVA_TOOL_OPTIONS=-Duser.home=...), not %s', home);
    if exist(OutDir, 'dir') ~= 7, mkdir(OutDir); end
    checks = {};
    check = @(name, cond, msg) fprintf('%-46s %s  %s\n', name, ifelse(cond), msg);

    % ----- bst_git_commit vs git -----
    [st, sha] = system(sprintf('git -C "%s" rev-parse HEAD', bstDir));
    c = bst_git_commit(bstDir);
    checks{end+1} = (st == 0) && strcmp(c, strtrim(sha));
    check('bst_git_commit == git rev-parse HEAD', checks{end}, c);

    % ----- throwaway protocol from the subject's volumes -----
    addpath(bstDir); brainstorm setpath;
    mkdir(fullfile(home, '.brainstorm'));
    dbDir = fullfile(home, 'brainstorm_db'); mkdir(dbDir);
    iProtocol = 0; ProtocolsListInfo = repmat(db_template('ProtocolInfo'), 0); %#ok<NASGU>
    ProtocolsListSubjects = repmat(db_template('ProtocolSubjects'), 0); %#ok<NASGU>
    ProtocolsListStudies = repmat(db_template('ProtocolStudies'), 0); %#ok<NASGU>
    BrainStormDbDir = dbDir; DbVersion = 5.03; %#ok<NASGU>
    save(fullfile(home, '.brainstorm', 'brainstorm.mat'), 'iProtocol', 'ProtocolsListInfo', ...
         'ProtocolsListSubjects', 'ProtocolsListStudies', 'BrainStormDbDir', 'DbVersion');
    if ~brainstorm('status'), brainstorm server; end
    shim = fullfile(home, 'shims'); mkdir(shim);           % headless: no viewers
    for fn = {'view_mri', 'view_surface_data'}
        fid = fopen(fullfile(shim, [fn{1} '.m']), 'w');
        fprintf(fid, 'function varargout = %s(varargin)\nvarargout = cell(1, nargout);\nend\n', fn{1}); fclose(fid);
    end
    addpath(shim, '-begin');

    gui_brainstorm('CreateProtocol', 'test_pet_tracer', 0, 0, dbDir);
    [~, SubjectName] = fileparts(AnatDir);
    [sSubject, iSubject] = db_add_subject(SubjectName, [], 0, 0);
    subjDir = bst_fileparts(file_fullpath(sSubject.FileName));
    files = {'subjectimage_MRI_T1.mat', 'subjectimage_ASEG_volatlas.mat', 'subjectimage_Desikan-Killiany_volatlas.mat', ...
              'tess_cortex_pial_low.mat', 'tess_cortex_mid_low.mat', 'tess_cortex_white_low.mat'};
    for i = 1:numel(files), copyfile(fullfile(AnatDir, files{i}), subjDir); end
    % the stored static PET of each tracer becomes a (1-frame) registered base "PET <trc>"
    tracers = {};
    d = dir(fullfile(AnatDir, 'subjectimage_*.mat'));
    for i = 1:numel(d)
        w = load(fullfile(AnatDir, d(i).name), 'Comment');
        t = regexp(w.Comment, '^PET (\S+)_mean$', 'tokens', 'once');
        if isempty(t), continue; end
        s = load(fullfile(AnatDir, d(i).name));
        s.Comment = ['PET ' t{1}];
        save(fullfile(subjDir, ['subjectimage_' t{1} '_volpet.mat']), '-struct', 's');
        tracers{end+1} = t{1}; %#ok<AGROW>
    end
    assert(~isempty(tracers), 'No "PET <trc>_mean" volume in %s', AnatDir);
    db_reload_subjects(iSubject);
    sSubject = bst_get('Subject', iSubject);
    iT1 = find(strcmp({sSubject.Anatomy.Comment}, 'MRI'), 1);
    if isempty(iT1), iT1 = find(contains({sSubject.Anatomy.FileName}, 'MRI_T1'), 1); end
    sSubject.iAnatomy = iT1; bst_set('Subject', iSubject, sSubject);
    fprintf('protocol: %s, tracers: %s\n', SubjectName, strjoin(tracers, ', '));

    % ----- GTM (pr/pet-pvc-2 API, scanner PSF from metadata) then VLPP, per tracer -----
    trc = tracers{1};
    rec = bst_pet_tracer(iSubject, SubjectName, trc, 'gtm', [], 1, OutDir, c);
    p = rec.Provenance;
    checks{end+1} = abs(p.PsfFwhmMm - 2.5) < 1e-9 && ~isempty(strfind(p.PsfSource, 'HRRT'));
    check('gtm: PSF = HRRT 2.5 mm from metadata', checks{end}, sprintf('%g mm [%s]', p.PsfFwhmMm, p.PsfSource));
    checks{end+1} = strcmp(p.Method, 'gtm') && strcmp(p.BrainstormCommit, c);
    check('gtm: provenance method + commit', checks{end}, p.BrainstormCommit);
    checks{end+1} = strcmp(p.Reference.Mode, 'robust') && p.Reference.Erode == 1 && isequal(p.Reference.AsegLabels, [8 47]) ...
                    && p.Reference.nVoxUsed < p.Reference.nVoxMask && isfinite(rec.RefValue) && rec.RefValue > 0;
    check('gtm: robust eroded cerebellar reference', checks{end}, sprintf('ref %.4g, %d -> %d vox', ...
          rec.RefValue, p.Reference.nVoxMask, p.Reference.nVoxUsed));
    T = readtable(fullfile(OutDir, rec.GtmCsv), 'Delimiter', ',');
    checks{end+1} = height(T) >= 90 && all(isfinite(T.corrected)) && all(T.nvox > 0);
    check('gtm: regional table (>= 90 regions, finite)', checks{end}, sprintf('%d regions', height(T)));
    ctx = T.id >= 1000 & T.id < 3000;
    checks{end+1} = median(T.corrected(ctx) ./ T.observed(ctx)) > 1;
    check('gtm: cortical values raised by PVC', checks{end}, sprintf('median corrected/observed %.3f', ...
          median(T.corrected(ctx) ./ T.observed(ctx))));
    checks{end+1} = rec.GlobalCorticalSUVR > 0.5 && rec.GlobalCorticalSUVR < 5 && ~isempty(rec.SurfaceFile);
    check('gtm: global cortical SUVR + cortex map', checks{end}, sprintf('%.3f', rec.GlobalCorticalSUVR));
    j = jsondecode(jsonencode(rec));
    checks{end+1} = isfield(j, 'Provenance') && isfield(j.Provenance, 'PsfSource');
    check('gtm: record is JSON-encodable', checks{end}, '');

    recV = bst_pet_tracer(iSubject, SubjectName, trc, 'vlpp', 6, 0, OutDir, c);
    checks{end+1} = isempty(recV.Provenance.PsfFwhmMm) && recV.Provenance.SmoothFwhmMm == 6 && isempty(recV.GtmCsv) ...
                    && recV.GlobalCorticalSUVR > 0.5 && recV.GlobalCorticalSUVR < 5;
    check('vlpp: no PVC, 6 mm smoothing, SUVR', checks{end}, sprintf('%.3f (gtm %.3f)', ...
          recV.GlobalCorticalSUVR, rec.GlobalCorticalSUVR));

    ok = all([checks{:}]);
    fprintf('%d/%d checks PASS\n', nnz([checks{:}]), numel(checks));
    brainstorm stop;
end

function s = ifelse(c)
    if c, s = 'PASS'; else, s = 'FAIL'; end
end
