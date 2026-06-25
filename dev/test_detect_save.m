function test_detect_save()
% TEST_DETECT_SAVE: Navigate/Detect/Save contract — Detect stages events, only Save writes st.T.
%
% USAGE:  test_detect_save   % Brainstorm running in nogui (GuiLevel 0)
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};  pass = true;

    [linkFile, relData] = i_find_kernel_ds();
    if isempty(linkFile)
        fprintf('SKIPPED (no unconstrained kernel link)\n');
        fprintf('\n==== SUITE: %s ====\n', PF{pass+1});  return;
    end
    % Use the link form returned by i_find_kernel_ds (shared kernels have empty internal
    % DataFile; the link|kernel|data form resolves the recording for view_dynamics).
    R = linkFile;
    hFig = view_dynamics('FromResult', R);  drawnow;
    ctrl = bst_get('PanelControls', 'Dynamics');

    % select the alpha band in the navigator (sets st.curBand)
    ctrl.jFreqBand.setSelectedItem('alpha');  panel_bst_dynamics('OnFreqPreset');  drawnow;

    % ---------- T1: Detect installs events, writes NOTHING to st.T ----------
    st0 = getappdata(0,'DynamicsTarget');  nG0 = st0.T.nGroups;
    panel_bst_dynamics('OnDetect');  drawnow;
    evs = panel_record('GetEvents', [], 1);                       % all events on the recording
    haveWin  = ~isempty(evs) && any(~cellfun(@isempty, regexp({evs.label}, '^alpha \(', 'once')));
    havePeak = ~isempty(evs) && any(strcmpi({evs.label}, 'alpha_peak'));
    st1 = getappdata(0,'DynamicsTarget');
    ok1 = haveWin && havePeak && (st1.T.nGroups == nG0);          % events present, st.T unchanged
    fprintf('T1 detect->events: win=%d peak=%d stT_unchanged=%d => %s\n', haveWin, havePeak, (st1.T.nGroups==nG0), PF{ok1+1});
    pass = pass && ok1;

    % ---------- T2: Save detection promotes events -> atoms (numeric freq) ----------
    panel_bst_dynamics('OnSaveDetection');  drawnow;
    st = getappdata(0,'DynamicsTarget');  Td = st.T;  parents = {Td.Groups.parent};
    gW = find(cellfun(@isempty,parents) & arrayfun(@(k) strcmp(Td.Groups(k).bandName,'alpha') && size(Td.Groups(k).times,1)==2, 1:Td.nGroups), 1);
    chN = [];  if ~isempty(gW), chN = find(strcmpi(parents, Td.Groups(gW).label)); end
    lf = bst_atom('Get', Td.Groups(gW), 'freq');                  % numeric freq stamped
    ok2 = ~isempty(gW) && (numel(chN)==4) && (abs(lf.center-10.5)<1e-6) && (abs(lf.extent-2.5)<1e-6);
    fprintf('T2 save detection: window=%d children=%d freqC=%g freqW=%g => %s\n', ~isempty(gW), numel(chN), lf.center, lf.extent, PF{ok2+1});
    pass = pass && ok2;

    % ---------- T3: Clear removes the detection events ----------
    panel_bst_dynamics('OnDetect');  drawnow;                     % re-stage events
    panel_bst_dynamics('OnClearDetection');  drawnow;
    evs2 = panel_record('GetEvents', [], 1);
    stillDet = ~isempty(evs2) && any(strcmpi({evs2.label}, 'alpha_peak'));
    ok3 = ~stillDet;
    fprintf('T3 clear detection: detectionEventsGone=%d => %s\n', ~stillDet, PF{ok3+1});
    pass = pass && ok3;

    % ---------- T4: Save cursor commits ONE atom from st.nav ----------
    ctrl.jFreqC.setText('10');  ctrl.jFreqW.setText('2');  panel_bst_dynamics('OnAxisChange','freq');  drawnow;
    panel_bst_dynamics('OnMeasurement','Stream');  drawnow;       % Function = stream
    st = getappdata(0,'DynamicsTarget');  nOcc0 = 0;
    gC0 = find(arrayfun(@(k) strcmp(st.T.Groups(k).Function,'stream') && strcmp(st.T.Groups(k).bandName,'alpha'), 1:st.T.nGroups), 1);
    if ~isempty(gC0), nOcc0 = numel(st.T.Groups(gC0).times); end
    panel_time('SetCurrentTime', 22.0);  panel_bst_dynamics('OnAxisChange','time');  drawnow;
    panel_bst_dynamics('OnSaveCursor');  drawnow;
    st = getappdata(0,'DynamicsTarget');
    gC1 = find(arrayfun(@(k) strcmp(st.T.Groups(k).Function,'stream') && strcmp(st.T.Groups(k).bandName,'alpha'), 1:st.T.nGroups), 1);
    nOcc1 = 0; if ~isempty(gC1), nOcc1 = numel(st.T.Groups(gC1).times); end
    ok4 = ~isempty(gC1) && (nOcc1 == nOcc0 + 1);
    fprintf('T4 save cursor: group=%d occ %d->%d => %s\n', ~isempty(gC1), nOcc0, nOcc1, PF{ok4+1});
    pass = pass && ok4;

    % ---------- T5: Load atom into navigator round-trips the saved coords ----------
    st = getappdata(0,'DynamicsTarget');
    gL = find(arrayfun(@(k) strcmp(st.T.Groups(k).Function,'stream'), 1:st.T.nGroups), 1);   % the saved-cursor group
    % select occurrence 1 of that group: drive the tree to its stack then the list row
    st.occMap = [gL 1 0];  setappdata(0,'DynamicsTarget', st);   % emulate a selected occurrence
    panel_bst_dynamics('OnLoadAtom');  drawnow;
    st = getappdata(0,'DynamicsTarget');
    lt = bst_atom('Get', st.nav, 'time');  lf = bst_atom('Get', st.nav, 'freq');
    gt = bst_atom('Get', st.T.Groups(gL), 'time', 1);  gf = bst_atom('Get', st.T.Groups(gL), 'freq');
    ftxt = str2double(char(ctrl.jFreqC.getText()));
    bandOK = ~isempty(st.curBand) && ~isempty(st.T.Groups(gL).band) && isequal(numel(st.curBand),2) && all(abs(st.curBand(:)-st.T.Groups(gL).band(:))<1e-6);   % i_drive's curBand survives the final write
    ok5b = (abs(lt.center-gt.center)<1e-6) && (abs(lf.center-gf.center)<1e-6) && (abs(ftxt-gf.center)<1e-6) && bandOK;
    fprintf('T5b load atom: navTime=%g navFreq=%g freqField=%g curBand=%s => %s\n', lt.center, lf.center, ftxt, mat2str(st.curBand), PF{ok5b+1});
    pass = pass && ok5b;

    if ishandle(hFig), close(hFig); end
    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end

function [linkFile, relData] = i_find_kernel_ds()
    linkFile = '';
    relData = 'Subject01/S01_AEF_20131218_01_notch/data_block001_02.mat';
    [sStudy, ~] = bst_get('DataFile', relData);
    if isempty(sStudy), return; end
    comments = {sStudy.Result.Comment};  fnames = {sStudy.Result.FileName};
    isMN = ~cellfun(@isempty, regexp(comments, 'MN: MEG\(Unconstr\)', 'once')) & ...
           ~cellfun(@isempty, regexp(fnames,   'KERNEL', 'once'));
    for j = find(isMN)
        try
            r = in_bst_results(fnames{j}, 0, 'nComponents','ImagingKernel');
            if (r.nComponents==3) && ~isempty(r.ImagingKernel), linkFile = ['link|' fnames{j} '|' relData];  return; end
        catch
        end
    end
end
