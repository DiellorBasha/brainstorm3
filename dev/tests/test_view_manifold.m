function test_view_manifold()
% TEST_VIEW_MANIFOLD  Live-figure regression for the manifold viewer.
% Default view is a true light-grey wireframe (faces off) with the scalar layer
% on (a point cloud kept in sync with the cortex patch). The cloud follows the
% Smooth slider and the numeric Resect, and drops the hemisphere hidden by a
% left/right Resect. 'D' toggles the layer. Requires Brainstorm running with a
% registered manifold node.
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;

    ManifoldFile = local_find_manifold();
    assert(~isempty(ManifoldFile), 'No registered manifold node found; compute a manifold first.');
    M = load(file_fullpath(ManifoldFile), 'ParentSurface');
    nVert = size(in_tess_bst(M.ParentSurface).Vertices, 1);

    close(findobj(0, 'type', 'figure', 'Tag', '3DViz'));
    hFig = view_manifold(ManifoldFile);
    drawnow;
    [nPass,nFail] = chk('viewer returns a figure', ~isempty(hFig) && ishandle(hFig), nPass,nFail);

    hAx3D  = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    hPatch = findobj(hAx3D, 'Type', 'patch'); hP = hPatch(1);
    hCloud = findobj(hFig, 'Tag', 'manifoldScalar');
    [nPass,nFail] = chk('Axes3D exists', ~isempty(hAx3D), nPass,nFail);
    [nPass,nFail] = chk('cortex patch exists', ~isempty(hPatch), nPass,nFail);
    fc = get(hP, 'FaceColor');
    [nPass,nFail] = chk('wireframe: FaceColor none', ischar(fc) && strcmpi(fc,'none'), nPass,nFail);
    [nPass,nFail] = chk('patch markers off (separate cloud used)', strcmpi(get(hP,'Marker'),'none'), nPass,nFail);

    % scalar cloud on by default, one point per vertex, all visible
    [nPass,nFail] = chk('scalar cloud exists', ~isempty(hCloud), nPass,nFail);
    Cx = get(hCloud, 'XData');
    [nPass,nFail] = chk('cloud has nVert points', numel(Cx)==nVert, nPass,nFail);
    [nPass,nFail] = chk('all vertices visible by default', nnz(~isnan(Cx))==nVert, nPass,nFail);
    V = get(hP, 'Vertices');
    [nPass,nFail] = chk('cloud matches patch vertices', isequal(Cx(:), V(:,1)), nPass,nFail);

    % cloud follows Smooth: after smoothing, the cloud tracks the moved vertices
    panel_surface('SetSurfaceSmooth', hFig, 1, 0.6, 0); drawnow; drawnow;
    Vs = get(hP, 'Vertices'); Cxs = get(findobj(hFig,'Tag','manifoldScalar'), 'XData');
    [nPass,nFail] = chk('cloud follows Smooth', ~isequal(V,Vs) && isequal(Cxs(:), Vs(:,1)), nPass,nFail);
    panel_surface('SetSurfaceSmooth', hFig, 1, 0, 0); drawnow; drawnow;

    % left/right Resect (real toggle path) hides one hemisphere's dots, and
    % toggling it back off restores them all (the reported round-trip bug).
    bst_figures('SetCurrentFigure', hFig, '3D');
    panel_surface('SelectHemispheres', 'right'); drawnow; drawnow;
    Cxr = get(findobj(hFig,'Tag','manifoldScalar'), 'XData');
    [nPass,nFail] = chk('Resect right hides one hemisphere''s dots', nnz(~isnan(Cxr))==nVert/2, nPass,nFail);
    panel_surface('SelectHemispheres', 'none'); drawnow; drawnow;
    Cxn = get(findobj(hFig,'Tag','manifoldScalar'), 'XData');
    [nPass,nFail] = chk('Resect toggle-off restores all dots', nnz(~isnan(Cxn))==nVert, nPass,nFail);

    % D toggles the layer off then on
    KeyOnFig(hFig, 'd'); drawnow;
    [nPass,nFail] = chk('D toggles scalar off', strcmpi(get(findobj(hFig,'Tag','manifoldScalar'),'Visible'),'off'), nPass,nFail);
    KeyOnFig(hFig, 'd'); drawnow;
    [nPass,nFail] = chk('D toggles scalar back on', strcmpi(get(findobj(hFig,'Tag','manifoldScalar'),'Visible'),'on'), nPass,nFail);

    % dimension tabs: 3 tabs; Vector2/Vector3 switch off the scalar view, Scalar back on
    hT = findobj(hFig, 'Tag', 'manifoldDimTabs');
    titles = {}; if ~isempty(hT), titles = arrayfun(@(t) char(t.Title), hT.Children, 'UniformOutput', 0); end
    [nPass,nFail] = chk('3 dimension tabs (Scalar/Vector2/Vector3)', ...
        isequal(sort(titles(:)'), {'Scalar','Vector2','Vector3'}), nPass,nFail);
    if ~isempty(hT)
        cbT = get(hT, 'SelectionChangedFcn');
        tScalar = findobj(hT.Children, 'flat', 'Title', 'Scalar');
        tVec2   = findobj(hT.Children, 'flat', 'Title', 'Vector2');
        cbT(hT, struct('NewValue', tVec2,   'OldValue', tScalar)); drawnow;
        [nPass,nFail] = chk('Vector2 tab switches scalar off', strcmpi(get(findobj(hFig,'Tag','manifoldScalar'),'Visible'),'off'), nPass,nFail);
        cbT(hT, struct('NewValue', tScalar, 'OldValue', tVec2)); drawnow;
        [nPass,nFail] = chk('Scalar tab switches scalar on', strcmpi(get(findobj(hFig,'Tag','manifoldScalar'),'Visible'),'on'), nPass,nFail);
    end

    close(hFig);
    fprintf('\n==== test_view_manifold: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_view_manifold: %d test(s) FAILED.', nFail); end
end

function f = local_find_manifold()
    f = '';
    P = bst_get('ProtocolSubjects'); subs = [P.Subject];
    for s = 1:numel(subs)
        for k = 1:numel(subs(s).Surface)
            if isfield(subs(s).Surface(k),'Manifold') && ~isempty(subs(s).Surface(k).Manifold)
                f = subs(s).Surface(k).Manifold(1).FileName; return;
            end
        end
    end
end

function KeyOnFig(hFig, keyName)
    ev.Key = keyName; ev.Character = ''; ev.Modifier = {};
    cb = get(hFig, 'KeyPressFcn'); cb(hFig, ev);
end

function [p,f] = chk(name, c, p, f)
    if c, fprintf('  PASS  %s\n', name); p=p+1; else, fprintf('  FAIL  %s\n', name); f=f+1; end
end
