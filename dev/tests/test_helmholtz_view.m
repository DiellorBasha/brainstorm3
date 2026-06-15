function test_helmholtz_view()
% Live: open the Helmholtz view on a synthetic 3-comp source. It opens the NATIVE source
% display (vector field + norm); switching the scalar overrides the cortex coloring with
% the active-frame Helmholtz field (signed, diverging colormap) while the native vectors
% stay; cores toggle; back to Norm restores the source colormap; close cleans up.
% Authors: Diellor Basha, 2026
    nFail = 0;
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    nV = size(in_tess_bst(SurfaceFile,0).Vertices,1);
    R = db_template('resultsmat'); rng(0);
    R.ImageGridAmp = randn(3*nV,3)*1e-9; R.nComponents=3; R.Time=0:2;
    R.HeadModelType='surface'; R.SurfaceFile=file_short(SurfaceFile); R.Comment='SYN helmholtz src';
    srcFile = db_add(-3, R);

    % launch like the tree node: the source is not pre-displayed; view loads + displays it
    hFig = view_helmholtz(srcFile); drawnow;
    nFail = nFail + chk('view opens (loads unshown results)', ishandle(hFig));
    St = getappdata(hFig,'HelmholtzState');
    hAx = findobj(hFig,'-depth',1,'Tag','Axes3D'); hAx=hAx(1);
    nFail = nFail + chk('default scalar is Norm', strcmp(St.Scalar,'Norm'));
    nFail = nFail + chk('native source vectors shown', ~isempty(findobj(hFig,'Tag','SourceVectors')));
    nFail = nFail + chk('starts on the source colormap', strcmpi(getappdata(hFig,'Colormap').Type,'source'));

    % active-frame decomposition for comparison
    [~, iT] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    Jt = double(bst_memory('GetResultsValues', St.srcDS, St.srcResult, [], iT, 0));
    Ht = bst_dirac_helmholtz('Frame', St.Op, Jt);

    % switch to Curl: cortex coloring becomes the signed curl, diverging colormap, vectors kept
    view_helmholtz('SetScalar', hFig, 'Curl'); drawnow;
    TI = getappdata(hFig,'Surface');
    nFail = nFail + chk('Curl overrides cortex scalar', isequal(TI(St.iTess).Data, Ht.Curl));
    nFail = nFail + chk('Curl uses a diverging colormap', strcmpi(getappdata(hFig,'Colormap').Type,'stat2'));
    nFail = nFail + chk('Curl color axis symmetric', abs(sum(get(hAx,'CLim'))) < 1e-9*max(abs(get(hAx,'CLim'))) + eps);
    nFail = nFail + chk('vectors still shown under Curl', ~isempty(findobj(hFig,'Tag','SourceVectors')));

    % cores match the active frame, and toggle off
    nFail = nFail + chk('core markers match active frame', numel(findobj(hAx,'Tag','HelmholtzCore')) == numel(Ht.Cores));
    view_helmholtz('SetCores', hFig, false); drawnow;
    nFail = nFail + chk('cores hide when toggled off', isempty(findobj(hAx,'Tag','HelmholtzCore')));

    % back to Norm restores the source colormap; time-following stays live
    view_helmholtz('SetScalar', hFig, 'Norm'); drawnow;
    nFail = nFail + chk('Norm restores source colormap', strcmpi(getappdata(hFig,'Colormap').Type,'source'));
    nFail = nFail + chk('figure not static (follows time)', ~isequal(getappdata(hFig,'isStatic'),1));
    view_helmholtz('SetScalar', hFig, 'Psi'); drawnow;
    panel_time('SetCurrentTime', 2.0); drawnow;
    [~, iT2] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    nFail = nFail + chk('cursor move decomposed a new frame', isKey(St.Cache, iT2) && (St.Cache.Count >= 2));

    % close cleans up + stale dispatch is ignored
    view_helmholtz('Close', hFig); drawnow;
    nFail = nFail + chk('panel closed', isempty(bst_get('PanelControls','Helmholtz')));
    nFail = nFail + chk('figure closed', ~ishandle(hFig));
    okGuard = true;
    try, view_helmholtz('SetScalar', hFig, 'Div'); catch; okGuard = false; end
    try, view_helmholtz('SetCores', hFig, true);   catch; okGuard = false; end
    nFail = nFail + chk('stale dispatch is ignored (no error)', okGuard);

    [~,iSt] = bst_get('AnyFile', srcFile); file_delete(file_fullpath(srcFile),1); db_reload_studies(iSt);
    fprintf('\n==== test_helmholtz_view: %d failed ====\n', nFail);
    if nFail > 0, error('test_helmholtz_view FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
