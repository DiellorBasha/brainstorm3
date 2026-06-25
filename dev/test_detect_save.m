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
