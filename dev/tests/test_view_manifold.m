function test_view_manifold()
% TEST_VIEW_MANIFOLD  Live-figure regression for the manifold viewer.
% Default view is a true light-grey wireframe (faces off) with the scalar layer
% on, drawn as the cortex patch's own vertex markers so it follows Smooth/Resect.
% 'D' toggles the scalar layer. Requires Brainstorm running with a registered
% manifold node.
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
    hPatch = findobj(hAx3D, 'Type', 'patch');
    [nPass,nFail] = chk('Axes3D exists', ~isempty(hAx3D), nPass,nFail);
    [nPass,nFail] = chk('cortex patch exists', ~isempty(hPatch), nPass,nFail);
    hP = hPatch(1);
    % default is a true wireframe: faces off (FaceColor 'none')
    fc = get(hP, 'FaceColor');
    [nPass,nFail] = chk('wireframe: FaceColor none', ischar(fc) && strcmpi(fc,'none'), nPass,nFail);

    % scalar layer on by default = patch vertex markers over all vertices
    [nPass,nFail] = chk('scalar on: patch Marker = ''.''', strcmp(get(hP,'Marker'),'.'), nPass,nFail);
    [nPass,nFail] = chk('markers cover all vertices', size(get(hP,'Vertices'),1)==nVert, nPass,nFail);

    % D toggles the scalar layer off, then back on
    KeyOnFig(hFig, 'd'); drawnow;
    [nPass,nFail] = chk('D toggles scalar off (Marker none)', strcmpi(get(hP,'Marker'),'none'), nPass,nFail);
    KeyOnFig(hFig, 'd'); drawnow;
    [nPass,nFail] = chk('D toggles scalar back on', strcmp(get(hP,'Marker'),'.'), nPass,nFail);

    % markers follow Smooth: vertices move, the markers (being the vertices) stay on
    V0 = get(hP, 'Vertices');
    panel_surface('SetSurfaceSmooth', hFig, 1, 0.6, 0); drawnow;
    V1 = get(hP, 'Vertices');
    [nPass,nFail] = chk('scalar follows Smooth (verts move, markers stay)', ...
        ~isequal(V0,V1) && strcmp(get(hP,'Marker'),'.'), nPass,nFail);
    panel_surface('SetSurfaceSmooth', hFig, 1, 0, 0); drawnow;

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
