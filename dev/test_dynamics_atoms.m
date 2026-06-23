function test_dynamics_atoms()
% TEST_DYNAMICS_ATOMS: Milestone-1 regression for the spatiotemporal atom system.
%
% Covers: db_template('atom'/'dynamicsmat'), bst_dynamics New/Add/Save/Load
% round-trip, the process_source_atoms populate, and a view_dynamics smoke test
% (figure markers + panel). The populate + viewer tests SKIP if no unconstrained
% kernel link is available for the auditory block.
%
% USAGE:  test_dynamics_atoms   % from MATLAB with Brainstorm running
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};
    pass = true;

    % ---------- T1: data model + I/O round-trip ----------
    A = bst_dynamics('NewAtom');
    A.category='vortex'; A.time=1.5; A.vertex=10; A.pos=[1 2 3]; A.band=[8 13];
    A.charge=int8(1); A.descriptors.note='x';
    T = bst_dynamics('New','unit');
    T = bst_dynamics('Add', T, A);
    T = bst_dynamics('Add', T, A);
    f = fullfile(bst_get('BrainstormTmpDir'), 'dynamics_unit.mat');
    bst_dynamics('Save', f, T);
    T2 = bst_dynamics('Load', f);
    ok1 = (T2.nAtoms==2) && strcmp(T2.Atoms(2).category,'vortex') && ...
          isequal(fieldnames(T2.Atoms), fieldnames(db_template('atom'))) && ...
          isequal(T2.Atoms(1).band, [8 13]);
    fprintf('T1 model round-trip: nAtoms=%d fieldsMatch=%d => %s\n', T2.nAtoms, ...
        isequal(fieldnames(T2.Atoms), fieldnames(db_template('atom'))), PF{ok1+1});
    pass = pass && ok1;

    % ---------- T2 + T3: populate + viewer (need a kernel link) ----------
    [linkFile, relData] = i_find_kernel();
    if isempty(linkFile)
        fprintf('T2/T3: SKIPPED (no unconstrained kernel link found)\n');
        fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
        return;
    end

    % T2: populate
    bst_process('CallProcess', 'process_source_atoms', linkFile, [], ...
        'eventname','alpha_peak', 'freqband','alpha', 'npeaks',3);
    studyDir = bst_fileparts(file_fullpath(relData));
    d = dir(fullfile(studyDir, 'dynamics_*.mat'));
    [~,ix] = max([d.datenum]);
    dynFile = fullfile(studyDir, d(ix).name);
    Tp = bst_dynamics('Load', dynFile);
    valid = ~isempty(Tp.Atoms);
    for i = 1:Tp.nAtoms
        a = Tp.Atoms(i);
        valid = valid && ~isempty(a.vertex) && ~isempty(a.pos) && numel(a.pos)==3 && ...
                ~isempty(a.time) && isequal(a.band,[8 13]) && ~isempty(a.DataFile) && ~isempty(a.SurfaceFile);
    end
    ok2 = (Tp.nAtoms >= 3) && valid;
    fprintf('T2 populate: nAtoms=%d allValidRefs=%d => %s\n', Tp.nAtoms, valid, PF{ok2+1});
    pass = pass && ok2;

    % T3: viewer smoke (open, check markers + panel, then close)
    [hFig, Tv] = view_dynamics(dynFile);
    nMark = numel(findobj(hFig, 'Tag', 'AtomMarker'));
    hasSel = ~isempty(findobj(hFig, 'Tag', 'AtomSel'));
    ctrl = bst_get('PanelControls', 'Dynamics');
    panelOK = ~isempty(ctrl) && (ctrl.jList.getModel().getSize() == Tv.nAtoms);
    ok3 = ishandle(hFig) && (nMark >= 1) && hasSel && panelOK;
    fprintf('T3 viewer: fig=%d markerGroups=%d sel=%d panelItems=%d => %s\n', ...
        ishandle(hFig), nMark, hasSel, panelOK*Tv.nAtoms, PF{ok3+1});
    pass = pass && ok3;

    % cleanup: close figure, remove the test dynamics file (not DB-registered)
    if ishandle(hFig), close(hFig); end
    if exist(dynFile,'file'), delete(dynFile); end

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end


%% ===== FIND AN UNCONSTRAINED KERNEL LINK FOR THE AUDITORY BLOCK =====
function [linkFile, relData] = i_find_kernel()
    linkFile = '';
    relData = 'Subject01/S01_AEF_20131218_01_notch/data_block001_02.mat';
    [sStudy, ~] = bst_get('DataFile', relData);
    if isempty(sStudy), return; end
    comments = {sStudy.Result.Comment};
    fnames   = {sStudy.Result.FileName};
    isMN = ~cellfun(@isempty, regexp(comments, 'MN: MEG\(Unconstr\)', 'once')) & ...
           ~cellfun(@isempty, regexp(fnames,   'KERNEL', 'once'));
    for j = find(isMN)
        try
            r = in_bst_results(fnames{j}, 0, 'nComponents','ImagingKernel');
            if (r.nComponents==3) && ~isempty(r.ImagingKernel)
                linkFile = ['link|' fnames{j} '|' relData];
                return;
            end
        catch
        end
    end
end
