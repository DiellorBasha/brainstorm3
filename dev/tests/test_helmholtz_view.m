function test_helmholtz_view()
% Live: open the Helmholtz view on a synthetic 3-comp source, switch the scalar, check
% the displayed series swaps and the core overlay count matches, then close cleanly.
% Authors: Diellor Basha, 2026
    global GlobalData;
    nFail = 0;
    % synthetic unconstrained source on cortex_20484V
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    nV = size(in_tess_bst(SurfaceFile,0).Vertices,1);
    R = db_template('resultsmat'); rng(0);
    R.ImageGridAmp = randn(3*nV,3)*1e-9; R.nComponents=3; R.Time=0:2;
    R.HeadModelType='surface'; R.SurfaceFile=file_short(SurfaceFile); R.Comment='SYN helmholtz src';
    srcFile = db_add(-3, R);
    [hSrc,~] = view_surface_data(SurfaceFile, srcFile, [], 'NewFigure'); drawnow;

    hFig = view_helmholtz(srcFile); drawnow;
    nFail = nFail + chk('view opens', ishandle(hFig));
    St = getappdata(hFig,'HelmholtzState');
    nFail = nFail + chk('default scalar is Curl', strcmp(St.Scalar,'Curl'));
    [iDS,iResult] = bst_memory('GetDataSetResult', St.TmpFile);
    nFail = nFail + chk('displays Curl series', isequal(GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp, St.H.Curl));
    % switch to Psi
    view_helmholtz('SetScalar', hFig, 'Psi'); drawnow;
    nFail = nFail + chk('switch to Psi swaps series', isequal(GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp, St.H.Psi));
    % time-following: the scalar figure must not be static, and a cursor move must fire
    % the overlay refresh (the path that keeps quiver + cores in sync with the colormap)
    nFail = nFail + chk('figure not static (follows time)', ~isequal(getappdata(hFig,'isStatic'),1));
    realFcn = getappdata(hFig,'CustomOverlayFcn');
    assignin('base','HVK',0);
    setappdata(hFig,'CustomOverlayFcn',@(h) assignin('base','HVK',evalin('base','HVK')+1));
    panel_time('SetCurrentTime', 1.0); drawnow;
    setappdata(hFig,'CustomOverlayFcn', realFcn);
    nFail = nFail + chk('cursor move fires overlay refresh', evalin('base','HVK') >= 1);
    % core overlay count matches H.Cores at frame 1
    hAx = findobj(hFig,'-depth',1,'Tag','Axes3D'); hAx=hAx(1);
    nMarkers = numel(findobj(hAx,'Tag','HelmholtzCore'));
    nFail = nFail + chk('core markers match', nMarkers == numel(St.H.Cores{1}));
    % close cleans up
    view_helmholtz('Close', hFig); drawnow;
    nFail = nFail + chk('panel closed', isempty(bst_get('PanelControls','Helmholtz')));
    nFail = nFail + chk('figure closed', ~ishandle(hFig));
    % cleanup synthetic source
    try, bst_figures('DeleteFigure', hSrc, []); catch, end
    [~,iSt] = bst_get('AnyFile', srcFile); file_delete(file_fullpath(srcFile),1); db_reload_studies(iSt);

    fprintf('\n==== test_helmholtz_view: %d failed ====\n', nFail);
    if nFail > 0, error('test_helmholtz_view FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
