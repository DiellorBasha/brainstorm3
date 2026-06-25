function test_nav_panel()
% TEST_NAV_PANEL: the four-axis (center,extent) navigator drives bst_atom + the engines.
%
% USAGE:  test_nav_panel   % Brainstorm running in nogui (GuiLevel 0)
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};  pass = true;

    % open the panel via view_dynamics on a Dirac result (reuses the atom-suite fixture)
    [linkFile, relData] = i_find_kernel_nav();
    if isempty(linkFile)
        fprintf('SKIPPED (no unconstrained kernel link)\n');
        fprintf('\n==== SUITE: %s ====\n', PF{pass+1});  return;
    end
    % FromResult needs the kernel LINK (bare KERNEL result has no DataFile; trials share it).
    R = linkFile;
    hFig = view_dynamics('FromResult', R);  drawnow;
    ctrl = bst_get('PanelControls', 'Dynamics');

    % ---------- T1: controls struct has the new block handles, not the old ones ----------
    ok1 = all(isfield(ctrl, {'jFreqC','jFreqW','jFreqBand','jTimeC','jTimeW','jSrcC','jSrcW','jScaleC','jScaleW','jMeasPot','jMeasStr','jRegionTool'})) ...
       && ~isfield(ctrl,'jBands') && ~isfield(ctrl,'jSpaceStr');
    fprintf('T1 handles: newBlocks=%d oldGone=%d => %s\n', all(isfield(ctrl,{'jFreqC','jMeasStr'})), ~isfield(ctrl,'jBands'), PF{ok1+1});
    pass = pass && ok1;

    % ---------- T2: freq block -> st.nav freq Localization + display filter ----------
    ctrl.jFreqC.setText('10');  ctrl.jFreqW.setText('2');
    panel_bst_dynamics('OnAxisChange', 'freq');  drawnow;
    st = getappdata(0,'DynamicsTarget');
    lf = bst_atom('Get', st.nav, 'freq');
    ok2 = (abs(lf.center-10)<1e-9) && (abs(lf.extent-2)<1e-9) && isequal(st.curBand,[8 12]);
    fprintf('T2 freq block: center=%g extent=%g curBand=%s => %s\n', lf.center, lf.extent, mat2str(st.curBand), PF{ok2+1});
    pass = pass && ok2;

    % ---------- T3: band combobox preset fills the freq fields (alpha -> 10.5 / 2.5) ----------
    ctrl.jFreqBand.setSelectedItem('alpha');  panel_bst_dynamics('OnFreqPreset');  drawnow;
    c = str2double(char(ctrl.jFreqC.getText()));  w = str2double(char(ctrl.jFreqW.getText()));
    ok3 = (abs(c-10.5)<1e-6) && (abs(w-2.5)<1e-6);
    fprintf('T3 band preset: center=%g window=%g => %s\n', c, w, PF{ok3+1});
    pass = pass && ok3;

    % ---------- T4: measurement Psi -> view_helmholtz component + curOp ----------
    ctrl.jMeasStr.setSelected(true);
    panel_bst_dynamics('OnMeasurement', 'Solen');  drawnow;
    st = getappdata(0,'DynamicsTarget');
    St = getappdata(st.hFig, 'HelmholtzState');
    ok4 = strcmp(st.curOp,'Solen') && ~isempty(St) && strcmpi(St.Component,'Solen');
    fprintf('T4 measurement: curOp=%s figComp=%s => %s\n', st.curOp, St.Component, PF{ok4+1});
    pass = pass && ok4;

    % ---------- T5: Source block syncs from the geodesic tool ----------
    st = getappdata(0,'DynamicsTarget');  surf = st.T.SurfaceFile;
    if isempty(surf), rs = in_bst_results(R,0,'SurfaceFile');  surf = rs.SurfaceFile; end
    vi = round(size(in_tess_bst(surf,0).Vertices,1)/3);
    bst_geodesic_tool('Seed', surf, vi);
    panel_bst_dynamics('SyncSource');  drawnow;
    st = getappdata(0,'DynamicsTarget');
    cv = str2double(char(ctrl.jSrcC.getText()));
    rw = str2double(char(ctrl.jSrcW.getText()));
    ls = bst_atom('Get', st.nav, 'source');
    ok5 = (cv==vi) && (rw>0) && (ls.center==vi);
    fprintf('T5 source sync: center=%g radius=%g navSeed=%g => %s\n', cv, rw, ls.center, PF{ok5+1});
    pass = pass && ok5;

    % ---------- T6: typed Source window (mm) stores metres in st.nav ----------
    ctrl.jSrcC.setText('5000');  ctrl.jSrcW.setText('6');     % 6 mm
    panel_bst_dynamics('OnAxisChange','source');  drawnow;
    st = getappdata(0,'DynamicsTarget');  lsrc = bst_atom('Get', st.nav, 'source');
    ok6 = (lsrc.center==5000) && (abs(lsrc.extent-0.006)<1e-9);
    fprintf('T6 source units: center=%g extent_m=%g => %s\n', lsrc.center, lsrc.extent, PF{ok6+1});
    pass = pass && ok6;

    % ---------- T7: 'none' restores broadband (clears fields, no filter, curBand empty) ----------
    ctrl.jFreqBand.setSelectedItem('alpha');  panel_bst_dynamics('OnFreqPreset');  drawnow;   % filter on
    ctrl.jFreqBand.setSelectedItem('none');   panel_bst_dynamics('OnFreqPreset');  drawnow;   % broadband
    st = getappdata(0,'DynamicsTarget');  lf7 = bst_atom('Get', st.nav, 'freq');
    cleared = isempty(char(ctrl.jFreqC.getText())) && isempty(char(ctrl.jFreqW.getText()));
    ok7 = cleared && strcmp(lf7.state,'unlocalized') && isempty(st.curBand);
    fprintf('T7 none->broadband: fieldsCleared=%d navUnloc=%d curBandEmpty=%d => %s\n', cleared, strcmp(lf7.state,'unlocalized'), isempty(st.curBand), PF{ok7+1});
    pass = pass && ok7;

    % ---------- T8: Region toggle OFF clears the Source selection ----------
    bst_geodesic_tool('Seed', surf, vi);  panel_bst_dynamics('SyncSource');  drawnow;     % populate a selection
    ctrl.jRegionTool.setSelected(false);  panel_bst_dynamics('OnRegionTool');  drawnow;   % toggle OFF -> clear
    st = getappdata(0,'DynamicsTarget');  lsrc8 = bst_atom('Get', st.nav, 'source');
    clr = isempty(char(ctrl.jSrcC.getText())) && isempty(char(ctrl.jSrcW.getText()));
    ok8 = isempty(bst_geodesic_tool('GetState')) && clr && strcmp(lsrc8.state,'unlocalized');
    fprintf('T8 region clear: stateGone=%d fieldsCleared=%d navUnloc=%d => %s\n', isempty(bst_geodesic_tool('GetState')), clr, strcmp(lsrc8.state,'unlocalized'), PF{ok8+1});
    pass = pass && ok8;

    if ishandle(hFig), close(hFig); end
    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end

function [linkFile, relData] = i_find_kernel_nav()
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
