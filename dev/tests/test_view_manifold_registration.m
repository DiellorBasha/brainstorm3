function test_view_manifold_registration()
% TEST_VIEW_MANIFOLD_REGISTRATION  Live-figure regression for the registration-
% sphere manifold viewer. Requires Brainstorm running with a registered manifold
% node whose parent surface has a Reg.Sphere.
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;

    ManifoldFile = local_find_manifold();
    assert(~isempty(ManifoldFile), 'No registered manifold node found; compute a manifold first.');

    % expected pole count from the node
    M = load(file_fullpath(ManifoldFile), 'Embedded', 'Gauge', 'ParentSurface');
    nVert = size(in_tess_bst(M.ParentSurface).Vertices, 1);
    G = view_manifold('DeriveVertexFrame', M.Embedded, M.Gauge, nVert);
    nPoles = numel(G.Sing);

    close(findobj(0, 'type', 'figure', 'Tag', '3DViz'));
    hFig = view_manifold_registration(ManifoldFile);
    drawnow;
    [nPass,nFail] = chk('viewer returns a figure', ~isempty(hFig) && ishandle(hFig), nPass,nFail);

    hAx3D  = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    hPatch = findobj(hAx3D, 'Type', 'patch');
    hSing  = findobj(hFig, 'Tag', 'manifoldRegSing');
    [nPass,nFail] = chk('Axes3D exists', ~isempty(hAx3D), nPass,nFail);
    [nPass,nFail] = chk('sphere patch exists', ~isempty(hPatch), nPass,nFail);
    [nPass,nFail] = chk('pole markers drawn', ~isempty(hSing), nPass,nFail);
    [nPass,nFail] = chk('marker count = nPoles', ~isempty(hSing) && numel(get(hSing(1),'XData'))==nPoles, nPass,nFail);

    % P toggles the pins off
    KeyOnFig(hFig, 'p'); drawnow;
    [nPass,nFail] = chk('P toggles pins off', isempty(findobj(hFig,'Tag','manifoldRegSing')), nPass,nFail);

    close(hFig);
    fprintf('\n==== test_view_manifold_registration: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_view_manifold_registration: %d test(s) FAILED.', nFail); end
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
