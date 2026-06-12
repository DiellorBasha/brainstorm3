function test_view_manifold()
% TEST_VIEW_MANIFOLD  Live-figure regression for the manifold frame viewer.
% Requires Brainstorm running with a registered manifold node.
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;

    % --- find a registered manifold node ---
    ManifoldFile = local_find_manifold();
    assert(~isempty(ManifoldFile), 'No registered manifold node found; compute a manifold first.');

    close(findobj(0, 'type', 'figure', 'Tag', '3DViz'));
    hFig = view_manifold(ManifoldFile);
    drawnow;
    [nPass,nFail] = chk('viewer returns a figure', ~isempty(hFig) && ishandle(hFig), nPass,nFail);

    hAx3D  = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    hPatch = findobj(hAx3D, 'Type', 'patch');
    hU = findobj(hFig, 'Tag', 'tangentU');
    hV = findobj(hFig, 'Tag', 'tangentV');
    [nPass,nFail] = chk('Axes3D survives glyph draw', ~isempty(hAx3D), nPass,nFail);
    [nPass,nFail] = chk('cortex patch survives glyph draw', ~isempty(hPatch), nPass,nFail);
    [nPass,nFail] = chk('U glyphs drawn', ~isempty(hU) && numel(get(hU(1),'UData'))>0, nPass,nFail);
    [nPass,nFail] = chk('V glyphs drawn', ~isempty(hV), nPass,nFail);

    % N key adds the normal quiver
    KeyOnFig(hFig, 'n'); drawnow;
    [nPass,nFail] = chk('N toggles normal glyphs', ~isempty(findobj(hFig,'Tag','tangentN')), nPass,nFail);

    % density key changes glyph count
    n1 = numel(get(findobj(hFig,'Tag','tangentU'),'UData'));
    KeyOnFig(hFig, 'leftarrow'); drawnow;
    n2 = numel(get(findobj(hFig,'Tag','tangentU'),'UData'));
    [nPass,nFail] = chk('density key changes glyph count', n2 < n1, nPass,nFail);

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
