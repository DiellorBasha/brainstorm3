function test_dynamics_atoms()
% TEST_DYNAMICS_ATOMS: Milestone-1 regression for the spatiotemporal atom system.
%
% Covers: db_template('atomgroup'/'dynamicsmat'), bst_dynamics New/NewGroup/
% AddGroup/Load round-trip (with nesting + simple/extended), the
% process_source_atoms nested populate (window + 4 phase children), and a
% view_dynamics smoke test (per-group markers + group tree). The populate +
% viewer tests SKIP if no unconstrained kernel link is available.
%
% USAGE:  test_dynamics_atoms   % from MATLAB with Brainstorm running
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};
    pass = true;

    % ---------- T1: data model + I/O round-trip (nesting, simple/extended) ----------
    T = bst_dynamics('New', 'unit');
    W = bst_dynamics('NewGroup', 'win');
    W.times = [1 3; 2 4];  W.band = [8 13];  W.bandName = 'alpha';   % 2 rows -> extended
    T = bst_dynamics('AddGroup', T, W);
    P = bst_dynamics('NewGroup', 'win_peak');
    P.parent = 'win';  P.phase = 'peak';  P.Function = 'magnitude';
    P.times = [1.5 3.5];  P.vertices = [10 20];  P.pos = [1 2 3; 4 5 6];  % 1 row -> simple
    P.hemi = [1 2];  P.strength = [0.5 0.7];  P.band = [8 13];
    T = bst_dynamics('AddGroup', T, P);
    f = fullfile(bst_get('BrainstormTmpDir'), 'dynamics_unit.mat');
    bst_dynamics('Save', f, T);
    T2 = bst_dynamics('Load', f);
    ok1 = (T2.nGroups==2) ...
        && strcmp(T2.Groups(1).type,'extended') && (size(T2.Groups(1).times,1)==2) ...
        && strcmp(T2.Groups(2).type,'simple')   && (size(T2.Groups(2).times,1)==1) ...
        && strcmp(T2.Groups(2).parent,'win') ...
        && isequal(T2.Groups(2).pos, [1 2 3; 4 5 6]) ...
        && isequal(fieldnames(T2.Groups), fieldnames(db_template('atomgroup')));
    fprintf('T1 model round-trip: nGroups=%d nest=%d ext/simple=%d/%d => %s\n', T2.nGroups, ...
        strcmp(T2.Groups(2).parent,'win'), strcmp(T2.Groups(1).type,'extended'), strcmp(T2.Groups(2).type,'simple'), PF{ok1+1});
    pass = pass && ok1;

    % ---------- T2 + T3: populate + viewer (need a kernel link) ----------
    [linkFile, relData] = i_find_kernel();
    if isempty(linkFile)
        fprintf('T2/T3: SKIPPED (no unconstrained kernel link found)\n');
        fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
        return;
    end

    % T2: nested populate
    bst_process('CallProcess', 'process_source_atoms', linkFile, [], 'freqband','alpha', 'npeaks',2);
    studyDir = bst_fileparts(file_fullpath(relData));
    d = dir(fullfile(studyDir, 'dynamics_*.mat'));  [~,ix] = max([d.datenum]);
    dynFile = fullfile(studyDir, d(ix).name);
    Tp = bst_dynamics('Load', dynFile);
    iWin = find(strcmp({Tp.Groups.type}, 'extended'), 1);
    nChildren = sum(~cellfun(@isempty, {Tp.Groups.parent}));
    childrefsOK = true;
    for g = 1:Tp.nGroups
        G = Tp.Groups(g);
        if ~isempty(G.parent)
            childrefsOK = childrefsOK && any(strcmp({Tp.Groups.label}, G.parent)) && ...
                          (numel(G.vertices)==size(G.pos,1)) && (numel(G.vertices)==size(G.times,2));
        end
    end
    ok2 = (Tp.nGroups==5) && ~isempty(iWin) && (nChildren==4) && childrefsOK;
    fprintf('T2 populate: nGroups=%d (1 window + %d nested children) refsOK=%d => %s\n', Tp.nGroups, nChildren, childrefsOK, PF{ok2+1});
    pass = pass && ok2;

    % T3: viewer + panel_bst_dynamics (group tree -> occurrence list -> highlight)
    i_close_dyn_figs();
    try, gui_hide('Dynamics'); catch, end %#ok<CTCH>
    [hFig, Tv] = view_dynamics(dynFile);
    nMark = numel(findobj(hFig, '-regexp', 'Tag', '^AtomMarker'));   % one line per spatial group (4)
    ctrl  = bst_get('PanelControls', 'Dynamics');
    root  = ctrl.jTree.getModel().getRoot();
    nStacks  = root.getChildCount();                                 % 1 band stack (alpha)
    nWindows = root.getChildAt(0).getChildCount();                   % 8 time windows under it
    % band = top-level extended group; pick the window with the most phase atoms
    parents  = {Tv.Groups.parent};
    gBand    = find(cellfun(@isempty,parents) & arrayfun(@(g) size(Tv.Groups(g).times,1)==2, 1:Tv.nGroups), 1);
    children = find(strcmpi(parents, Tv.Groups(gBand).label));
    G = Tv.Groups(gBand);  bestW = 1;  bestN = -1;
    for w = 1:size(G.times,2)
        on = G.times(1,w);  off = G.times(2,w);  n = 0;
        for c = children, gc = Tv.Groups(c);  n = n + sum(gc.times(1,:)>=on-1e-9 & gc.times(1,:)<=off+1e-9);  end
        if n > bestN, bestN = n;  bestW = w;  end
    end
    % select that window node -> right list = flat phase atoms within it
    st = getappdata(0, 'DynamicsTarget');
    iWinNode = find(arrayfun(@(k) strcmp(st.nodeInfo(k).kind,'window') && st.nodeInfo(k).g==gBand && st.nodeInfo(k).w==bestW, 1:numel(st.nodeInfo)), 1);
    ctrl.jTree.setSelectionPath(javax.swing.tree.TreePath(st.nodeList{iWinNode}.getPath()));
    drawnow;
    nOcc = ctrl.jListOccur.getModel().getSize();
    % select first atom in the window -> highlight
    ctrl.jListOccur.setSelectedIndex(0);  drawnow;
    hSel = findobj(hFig, 'Tag', 'AtomSel');
    selOn = strcmp(get(hSel,'Visible'), 'on');
    ok3 = ishandle(hFig) && (nMark==4) && (nStacks==1) && (nWindows==8) && (nOcc>0) && selOn;
    fprintf('T3 panel: markers=%d stacks=%d windows=%d winAtoms=%d highlight=%d => %s\n', nMark, nStacks, nWindows, nOcc, selOn, PF{ok3+1});
    pass = pass && ok3;

    % cleanup
    if ishandle(hFig), close(hFig); end
    if exist(dynFile,'file'), delete(dynFile); end

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end


%% ===== CLOSE STALE DYNAMICS FIGURES =====
function i_close_dyn_figs()
    hAll = findall(0, 'Type', 'figure');
    for h = hAll(:)'
        if ~isempty(getappdata(h, 'GroupsPosOff')) || ~isempty(findobj(h, 'Tag', 'AtomMarker'))
            close(h);
        end
    end
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
