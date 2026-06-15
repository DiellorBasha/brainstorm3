function test_helmholtz_view()
% Live: open the Helmholtz view on a synthetic 3-comp source. It must decompose only the
% ACTIVE frame, display the selected scalar for that frame, match the core overlay count,
% recompute on cursor move (on-demand cache grows), and close cleanly.
% Authors: Diellor Basha, 2026
    nFail = 0;
    % synthetic unconstrained source on cortex_20484V (3 time frames)
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

    % the displayed scalar must equal the decomposition of the ACTIVE frame
    [~, iT] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    Jt = double(bst_memory('GetResultsValues', St.srcDS, St.srcResult, [], iT, 0));
    Ht = bst_dirac_helmholtz('Frame', St.Op, Jt);
    TessInfo = getappdata(hFig,'Surface');
    nFail = nFail + chk('displays active-frame Curl', isequal(TessInfo(St.iTess).Data, Ht.Curl));
    nFail = nFail + chk('only active frame cached', St.Cache.Count == 1);

    % switch to Psi -> recolors from the cached frame (no recompute)
    view_helmholtz('SetScalar', hFig, 'Psi'); drawnow;
    TessInfo = getappdata(hFig,'Surface');
    nFail = nFail + chk('switch to Psi recolors active frame', isequal(TessInfo(St.iTess).Data, Ht.Psi));

    % core overlay count matches the active frame's cores
    hAx = findobj(hFig,'-depth',1,'Tag','Axes3D'); hAx=hAx(1);
    nMarkers = numel(findobj(hAx,'Tag','HelmholtzCore'));
    nFail = nFail + chk('core markers match active frame', nMarkers == numel(Ht.Cores));

    % time-following: not static, and a cursor move decomposes the new frame on demand
    nFail = nFail + chk('figure not static (follows time)', ~isequal(getappdata(hFig,'isStatic'),1));
    panel_time('SetCurrentTime', 2.0); drawnow;
    [~, iT2] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    nFail = nFail + chk('cursor move decomposed a new frame', isKey(St.Cache, iT2) && (St.Cache.Count >= 2));

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
