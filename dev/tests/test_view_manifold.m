function test_view_manifold()
% TEST_VIEW_MANIFOLD  Live-figure regression for the manifold viewer.
% Default view is a bare wireframe (no overlays); 'D' toggles a scalar vertex
% point cloud (the scalar data dimension). Requires Brainstorm running with a
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
    [nPass,nFail] = chk('Axes3D exists', ~isempty(hAx3D), nPass,nFail);
    [nPass,nFail] = chk('cortex patch exists', ~isempty(findobj(hAx3D,'Type','patch')), nPass,nFail);
    % default view is bare: nothing overlaid
    [nPass,nFail] = chk('default bare (no scalar cloud)', isempty(findobj(hFig,'Tag','manifoldScalar')), nPass,nFail);

    % D draws the scalar point cloud over all vertices
    KeyOnFig(hFig, 'd'); drawnow;
    hS = findobj(hFig, 'Tag', 'manifoldScalar');
    [nPass,nFail] = chk('D draws scalar point cloud', ~isempty(hS), nPass,nFail);
    [nPass,nFail] = chk('scalar cloud covers all vertices', ~isempty(hS) && numel(get(hS(1),'XData'))==nVert, nPass,nFail);
    [nPass,nFail] = chk('cortex survives scalar draw', ~isempty(findobj(hAx3D,'Type','patch')), nPass,nFail);

    % D again toggles it off
    KeyOnFig(hFig, 'd'); drawnow;
    [nPass,nFail] = chk('D toggles scalar off', isempty(findobj(hFig,'Tag','manifoldScalar')), nPass,nFail);

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
