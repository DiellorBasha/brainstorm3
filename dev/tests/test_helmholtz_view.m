function test_helmholtz_view()
% Live: open the Helmholtz components view on a synthetic 3-comp source. Default Total shows
% the field + |J|; each component sets the quiver override to ITS vector field and colors by
% ITS scalar (signed for Irrot/Solen); markers are component-aware; cursor move recomputes;
% close cleans up.
% Authors: Diellor Basha, 2026
    nFail = 0;
    SurfaceFile = bst_get('Subject',1).Surface(5).FileName;
    nV = size(in_tess_bst(SurfaceFile,0).Vertices,1);
    R = db_template('resultsmat'); rng(0);
    R.ImageGridAmp = randn(3*nV,3)*1e-9; R.nComponents=3; R.Time=0:2;
    R.HeadModelType='surface'; R.SurfaceFile=file_short(SurfaceFile); R.Comment='SYN helmholtz src';
    srcFile = db_add(-3, R);

    hFig = view_helmholtz(srcFile); drawnow;
    nFail = nFail + chk('view opens', ishandle(hFig));
    St = getappdata(hFig,'HelmholtzState'); hAx = findobj(hFig,'-depth',1,'Tag','Axes3D'); hAx=hAx(1);
    nFail = nFail + chk('default component is Total', strcmp(St.Component,'Total'));
    nFail = nFail + chk('native source vectors shown', ~isempty(findobj(hFig,'Tag','SourceVectors')));

    [~, iT] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    Jt = double(bst_memory('GetResultsValues', St.srcDS, St.srcResult, [], iT, 0));
    Ht = bst_dirac_helmholtz('Frame', St.Op, Jt);

    % Total: quiver override == J, source colormap
    nFail = nFail + chk('Total quiver = J', isequal(getappdata(hFig,'QuiverVectorOverride'), Ht.Vtot));
    nFail = nFail + chk('Total colormap source', strcmpi(getappdata(hFig,'Colormap').Type,'source'));

    % Solenoidal: quiver STAYS the total field, colormap -> psi (stat2), markers = vortex cores
    view_helmholtz('SetComponent', hFig, 'Solen'); drawnow;
    nFail = nFail + chk('Solen quiver stays = total field', isequal(getappdata(hFig,'QuiverVectorOverride'), Ht.Vtot));
    nFail = nFail + chk('Solen colormap stat2', strcmpi(getappdata(hFig,'Colormap').Type,'stat2'));
    nFail = nFail + chk('Solen markers = vortex cores', numel(findobj(hAx,'Tag','HelmholtzCore'))==numel(Ht.Cores));

    % Irrotational: quiver STAYS the total field, markers = sources/sinks
    view_helmholtz('SetComponent', hFig, 'Irrot'); drawnow;
    nFail = nFail + chk('Irrot quiver stays = total field', isequal(getappdata(hFig,'QuiverVectorOverride'), Ht.Vtot));
    nFail = nFail + chk('Irrot markers = sources/sinks', numel(findobj(hAx,'Tag','HelmholtzCore'))==numel(Ht.Sources));

    % Harmonic: quiver STAYS the total field, no markers
    view_helmholtz('SetComponent', hFig, 'Harm'); drawnow;
    nFail = nFail + chk('Harm quiver stays = total field', isequal(getappdata(hFig,'QuiverVectorOverride'), Ht.Vtot));
    nFail = nFail + chk('Harm has no markers', isempty(findobj(hAx,'Tag','HelmholtzCore')));

    % toggles
    view_helmholtz('SetVectors', hFig, false); drawnow;
    nFail = nFail + chk('vectors hide when off', isempty(findobj(hFig,'Tag','SourceVectors')));
    view_helmholtz('SetComponent', hFig, 'Solen'); view_helmholtz('SetMarkers', hFig, false); drawnow;
    nFail = nFail + chk('markers hide when off', isempty(findobj(hAx,'Tag','HelmholtzCore')));

    % time-following recompute
    nFail = nFail + chk('figure not static', ~isequal(getappdata(hFig,'isStatic'),1));
    panel_time('SetCurrentTime', 2.0); drawnow;
    [~, iT2] = bst_memory('GetTimeVector', St.srcDS, St.srcResult, 'CurrentTimeIndex');
    nFail = nFail + chk('cursor move decomposed a new frame', isKey(St.Cache, iT2) && (St.Cache.Count >= 2));

    % --- smoothing: a heat low-pass on the active frame thins the vortex cores ---
    view_helmholtz('SetComponent', hFig, 'Solen'); view_helmholtz('SetMarkers', hFig, true);
    view_helmholtz('SetGate', hFig, 0); drawnow;
    nRaw = numel(findobj(hAx,'Tag','HelmholtzCore'));
    St = getappdata(hFig,'HelmholtzState');  Lam = St.Lambda;
    tt = 1 / Lam(max(1, round(numel(Lam)/3)));               % heat scale at a low-ish mode
    view_helmholtz('SetSmoothing', hFig, true, 'heat', struct('t',tt)); drawnow;
    nSmooth = numel(findobj(hAx,'Tag','HelmholtzCore'));
    nFail = nFail + chk('smoothing thins vortex cores', nSmooth < nRaw);
    St2 = getappdata(hFig,'HelmholtzState');
    nFail = nFail + chk('smoothing change cleared the cache', St2.Cache.Count <= 1);

    % --- magnitude gate is monotonic (more pruning -> fewer markers) ---
    view_helmholtz('SetGate', hFig, 0);   drawnow; g0 = numel(findobj(hAx,'Tag','HelmholtzCore'));
    view_helmholtz('SetGate', hFig, 0.5); drawnow; g5 = numel(findobj(hAx,'Tag','HelmholtzCore'));
    view_helmholtz('SetGate', hFig, 0.95);drawnow; g9 = numel(findobj(hAx,'Tag','HelmholtzCore'));
    nFail = nFail + chk('gate monotonic (g0>=g5>=g9)', (g0>=g5) && (g5>=g9));
    view_helmholtz('SetSmoothing', hFig, false, 'heat', struct('t',tt)); view_helmholtz('SetGate', hFig, 0); drawnow;

    % close + stale guard
    view_helmholtz('Close', hFig); drawnow;
    nFail = nFail + chk('panel closed', isempty(bst_get('PanelControls','Helmholtz')));
    nFail = nFail + chk('figure closed', ~ishandle(hFig));
    okGuard = true;
    try, view_helmholtz('SetComponent', hFig, 'Irrot'); catch; okGuard = false; end
    nFail = nFail + chk('stale dispatch ignored', okGuard);

    [~,iSt] = bst_get('AnyFile', srcFile); file_delete(file_fullpath(srcFile),1); db_reload_studies(iSt);
    fprintf('\n==== test_helmholtz_view: %d failed ====\n', nFail);
    if nFail > 0, error('test_helmholtz_view FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
