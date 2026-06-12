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
    [nPass,nFail] = chk('surface render on (faces visible)', ~(ischar(fc) && strcmpi(fc,'none')), nPass,nFail);
    [nPass,nFail] = chk('patch markers off (separate cloud used)', strcmpi(get(hP,'Marker'),'none'), nPass,nFail);

    % scalar cloud on by default, one point per vertex, all visible
    [nPass,nFail] = chk('scalar cloud exists', ~isempty(hCloud), nPass,nFail);
    Cx = get(hCloud, 'XData');
    [nPass,nFail] = chk('cloud has nVert points', numel(Cx)==nVert, nPass,nFail);
    [nPass,nFail] = chk('all vertices visible by default', nnz(~isnan(Cx))==nVert, nPass,nFail);
    % cloud sits near the vertices (lifted slightly off the surface)
    V = get(hP, 'Vertices');
    Ftest = double(get(hP, 'Faces'));
    meTest = mean(sqrt(sum((V(Ftest(:,1),:) - V(Ftest(:,2),:)).^2, 2)));
    Cxyz = [get(hCloud,'XData')', get(hCloud,'YData')', get(hCloud,'ZData')'];
    [nPass,nFail] = chk('cloud near patch vertices (lifted)', max(sqrt(sum((Cxyz - V).^2,2))) < meTest, nPass,nFail);

    % cloud follows Smooth: after smoothing, the cloud moves with the vertices
    Cx0 = get(hCloud, 'XData');
    panel_surface('SetSurfaceSmooth', hFig, 1, 0.6, 0); drawnow; drawnow;
    Cxs = get(findobj(hFig,'Tag','manifoldScalar'), 'XData');
    [nPass,nFail] = chk('cloud follows Smooth (moves)', ~isequal(Cx0,Cxs) && nnz(~isnan(Cxs))==nVert, nPass,nFail);
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
        tVec3   = findobj(hT.Children, 'flat', 'Title', 'Vector3');
        % Vector2: scalar off, tangent frame (U,V) drawn, no normal
        cbT(hT, struct('NewValue', tVec2, 'OldValue', tScalar)); drawnow;
        [nPass,nFail] = chk('Vector2: scalar off', strcmpi(get(findobj(hFig,'Tag','manifoldScalar'),'Visible'),'off'), nPass,nFail);
        [nPass,nFail] = chk('Vector2: U,V frames drawn, no N', ...
            ~isempty(findobj(hFig,'Tag','manifoldVecU')) && ~isempty(findobj(hFig,'Tag','manifoldVecV')) ...
            && isempty(findobj(hFig,'Tag','manifoldVecN')), nPass,nFail);
        % Vector3: full frame U,V,N
        cbT(hT, struct('NewValue', tVec3, 'OldValue', tVec2)); drawnow;
        [nPass,nFail] = chk('Vector3: U,V,N frames drawn', ...
            ~isempty(findobj(hFig,'Tag','manifoldVecU')) && ~isempty(findobj(hFig,'Tag','manifoldVecV')) ...
            && ~isempty(findobj(hFig,'Tag','manifoldVecN')), nPass,nFail);
        % back to Scalar: vector glyphs cleared, scalar on
        cbT(hT, struct('NewValue', tScalar, 'OldValue', tVec3)); drawnow;
        [nPass,nFail] = chk('Scalar: vectors cleared, scalar on', ...
            isempty(findobj(hFig,'Tag','manifoldVecU')) && strcmpi(get(findobj(hFig,'Tag','manifoldScalar'),'Visible'),'on'), nPass,nFail);
    end

    % support buttons: Vertex (nVert) <-> Face (nFace centroids), scalar layer
    nFace = size(get(hP,'Faces'), 1);
    hS = findobj(hFig, 'Tag', 'manifoldSupport');
    [nPass,nFail] = chk('support buttons present', ~isempty(hS), nPass,nFail);
    if ~isempty(hS)
        cbS   = get(hS, 'SelectionChangedFcn');
        tVert = findobj(hS, 'String', 'Vertex');
        tFace = findobj(hS, 'String', 'Face');
        [nPass,nFail] = chk('Vertex support: nVert cloud points', ...
            numel(get(findobj(hFig,'Tag','manifoldScalar'),'XData'))==nVert, nPass,nFail);
        cbS(hS, struct('NewValue', tFace, 'OldValue', tVert)); drawnow;
        [nPass,nFail] = chk('Face support: nFace cloud points', ...
            numel(get(findobj(hFig,'Tag','manifoldScalar'),'XData'))==nFace, nPass,nFail);
        cbS(hS, struct('NewValue', tVert, 'OldValue', tFace)); drawnow;
        [nPass,nFail] = chk('back to Vertex: nVert cloud points', ...
            numel(get(findobj(hFig,'Tag','manifoldScalar'),'XData'))==nVert, nPass,nFail);
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
