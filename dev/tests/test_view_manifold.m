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
    M = load(file_fullpath(ManifoldFile), 'ParentSurface', 'Embedded', 'Gauge');
    Tess  = in_tess_bst(M.ParentSurface);
    nVert = size(Tess.Vertices, 1);
    nFaceTot = size(Tess.Faces, 1);

    % ---- pure DeriveFaceFrame: combed field on the dual mesh (face centroids) ----
    GF = view_manifold('DeriveFaceFrame', M.Embedded, M.Gauge, nFaceTot);
    [nPass,nFail] = chk('DeriveFaceFrame: P,U,V,N are [nFace x 3]', ...
        isequal(size(GF.P),[nFaceTot 3]) && isequal(size(GF.U),[nFaceTot 3]) ...
        && isequal(size(GF.V),[nFaceTot 3]) && isequal(size(GF.N),[nFaceTot 3]), nPass,nFail);
    onF = any(GF.U ~= 0, 2);   % faces actually populated (both hemispheres)
    Un  = sqrt(sum(GF.U(onF,:).^2, 2));
    dUV = abs(sum(GF.U(onF,:) .* GF.V(onF,:), 2));
    [nPass,nFail] = chk('DeriveFaceFrame: combed U is unit-norm', max(abs(Un-1)) < 1e-6, nPass,nFail);
    [nPass,nFail] = chk('DeriveFaceFrame: U orthogonal to V', max(dUV) < 1e-6, nPass,nFail);

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

    % Vector layers honor the support switch: Face draws the combed frame on the
    % dual mesh (centroids), Vertex on the vertices — different geometry, so the
    % glyph base positions change when the support is switched.
    hSv = findobj(hFig, 'Tag', 'manifoldSupport');
    if ~isempty(hT) && ~isempty(hSv)
        cbT = get(hT, 'SelectionChangedFcn');  cbS = get(hSv, 'SelectionChangedFcn');
        tScalar = findobj(hT.Children, 'flat', 'Title', 'Scalar');
        tVec2   = findobj(hT.Children, 'flat', 'Title', 'Vector2');
        tVert   = findobj(hSv, 'String', 'Vertex');
        tFace   = findobj(hSv, 'String', 'Face');
        getBase = @() local_quiver_base(findobj(hFig,'Tag','manifoldVecU'));
        cbT(hT, struct('NewValue', tVec2, 'OldValue', tScalar)); drawnow;   % enter Vector2 (vertex)
        baseVert = getBase();
        [nPass,nFail] = chk('Vector2 + Vertex: frames drawn', ~isempty(baseVert), nPass,nFail);
        % glyph arms are sized per-cell (kGlyph*sqrt(dualArea)), so arm length VARIES
        % across the mesh — a fixed global length would give ~zero spread.
        qUv = findobj(hFig,'Tag','manifoldVecU');
        au = get(qUv,'UData'); av = get(qUv,'VData'); aw = get(qUv,'WData');
        armLen = sqrt(au(:).^2 + av(:).^2 + aw(:).^2); armLen = armLen(armLen > 0);
        covArm = std(armLen) / mean(armLen);
        [nPass,nFail] = chk('Vector2: arm length scales per-cell (varies, not fixed)', covArm > 0.15, nPass,nFail);
        cbS(hSv, struct('NewValue', tFace, 'OldValue', tVert)); drawnow;    % -> Face support
        baseFace = getBase();
        [nPass,nFail] = chk('Vector2 + Face: frames drawn', ~isempty(baseFace), nPass,nFail);
        [nPass,nFail] = chk('Vector2: Face support moves glyph bases off the vertices', ...
            ~isequaln(baseVert, baseFace), nPass,nFail);
        cbS(hSv, struct('NewValue', tVert, 'OldValue', tFace)); drawnow;    % back to Vertex
        [nPass,nFail] = chk('Vector2 + Vertex: bases restored', isequaln(getBase(), baseVert), nPass,nFail);
        cbT(hT, struct('NewValue', tScalar, 'OldValue', tVec2)); drawnow;   % leave clean (Scalar, Vertex)
    end

    % Align-to-Freesurfer toggle + blue example data field (Vector2). Default frame
    % is the raw grid; Align rotates it to the FS gauge; the blue data field (drawn
    % like the scalar cloud) is gauge-invariant and must NOT move.
    hAlignB = findobj(hFig, 'Tag', 'manifoldAlignFS');
    [nPass,nFail] = chk('Align-to-Freesurfer button present', ~isempty(hAlignB), nPass,nFail);
    if ~isempty(hT) && ~isempty(hAlignB)
        cbT = get(hT, 'SelectionChangedFcn');  cbA = get(hAlignB, 'Callback');
        tScalar = findobj(hT.Children, 'flat', 'Title', 'Scalar');
        tVec2   = findobj(hT.Children, 'flat', 'Title', 'Vector2');
        [nPass,nFail] = chk('Align disabled outside vector dims', strcmpi(get(hAlignB,'Enable'),'off'), nPass,nFail);
        cbT(hT, struct('NewValue', tVec2, 'OldValue', tScalar)); drawnow;   % enter Vector2
        [nPass,nFail] = chk('Align enabled in Vector2', strcmpi(get(hAlignB,'Enable'),'on'), nPass,nFail);
        qD = findobj(hFig, 'Tag', 'manifoldVecData');
        [nPass,nFail] = chk('Vector2: blue example data field drawn', ~isempty(qD), nPass,nFail);
        [nPass,nFail] = chk('data field colored like scalar (blue)', isequal(get(qD,'Color'), [0.2 0.9 1]), nPass,nFail);
        frameOff = local_quiver_arms(findobj(hFig,'Tag','manifoldVecU'));
        dataOff  = local_quiver_arms(findobj(hFig,'Tag','manifoldVecData'));
        set(hAlignB, 'Value', 1); cbA(hAlignB, []); drawnow;               % Align ON
        frameOn = local_quiver_arms(findobj(hFig,'Tag','manifoldVecU'));
        dataOn  = local_quiver_arms(findobj(hFig,'Tag','manifoldVecData'));
        [nPass,nFail] = chk('Align rotates the reference frame', ~isequaln(frameOff, frameOn), nPass,nFail);
        [nPass,nFail] = chk('Align leaves the data field fixed', isequaln(dataOff, dataOn), nPass,nFail);
        % when aligned, frame e1 is parallel to the combed data; when raw, generally not
        nrm = @(A) A ./ max(sqrt(sum(A.^2,2)), 1e-12);
        cosOn  = mean(abs(sum(nrm(frameOn) .* nrm(dataOn), 2)));
        cosOff = mean(abs(sum(nrm(frameOff) .* nrm(dataOff), 2)));
        [nPass,nFail] = chk('Aligned frame e1 || data; raw frame e1 not', cosOn > 0.999 && cosOff < 0.99, nPass,nFail);
        set(hAlignB, 'Value', 0); cbA(hAlignB, []); drawnow;               % reset
        cbT(hT, struct('NewValue', tScalar, 'OldValue', tVec2)); drawnow;  % leave clean
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

function B = local_quiver_base(qU)
    if isempty(qU), B = []; return; end
    x = get(qU,'XData'); y = get(qU,'YData'); z = get(qU,'ZData');
    B = [x(:), y(:), z(:)];
end

function A = local_quiver_arms(qU)
    if isempty(qU), A = []; return; end
    u = get(qU,'UData'); v = get(qU,'VData'); w = get(qU,'WData');
    A = [u(:), v(:), w(:)];
end

function KeyOnFig(hFig, keyName)
    ev.Key = keyName; ev.Character = ''; ev.Modifier = {};
    cb = get(hFig, 'KeyPressFcn'); cb(hFig, ev);
end

function [p,f] = chk(name, c, p, f)
    if c, fprintf('  PASS  %s\n', name); p=p+1; else, fprintf('  FAIL  %s\n', name); f=f+1; end
end
