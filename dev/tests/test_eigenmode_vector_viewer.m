function test_eigenmode_vector_viewer()
% TEST_EIGENMODE_VECTOR_VIEWER  Live-figure regression for the Dirac eigenmode
% vector viewer. Requires Brainstorm running with a loaded protocol. Builds a
% canonical cortex + Dirac eigen node, opens view_eigenmodes, and asserts the
% cortex survives the quiver draw, arrows update on mode step, and a non-Dirac
% node is rejected. Created eigen/operator nodes are left in place (cheap).
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;

    % --- fixture: canonical cortex + Dirac eigen node ---
    SurfaceFile = bst_canonical_cortex(20484);
    preEig  = local_eigen_names(SurfaceFile);
    tess_eigen(SurfaceFile, 'Dirac', 'K', 40, 'Tau', 0.5);   % saves + registers
    postEig = local_eigen_names(SurfaceFile);
    newEig  = setdiff(postEig, preEig);
    assert(~isempty(newEig), 'No new Dirac eigen node was created.');
    EigenFile = newEig{1};

    % --- open the viewer ---
    close(findobj(0, 'type', 'figure', 'Tag', '3DViz'));
    hFig = view_eigenmodes(EigenFile);
    drawnow;
    [nPass,nFail] = chk('viewer returns a figure', ~isempty(hFig) && ishandle(hFig), nPass,nFail);

    hAx3D = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    hPatch = findobj(hAx3D, 'Type', 'patch');
    hQ = findobj(hFig, 'Tag', 'eigArrows');
    [nPass,nFail] = chk('Axes3D survives quiver draw', ~isempty(hAx3D), nPass,nFail);
    [nPass,nFail] = chk('cortex patch survives quiver draw', ~isempty(hPatch), nPass,nFail);
    [nPass,nFail] = chk('arrows drawn', ~isempty(hQ) && numel(get(hQ(1),'UData')) > 0, nPass,nFail);

    % --- mode stepping updates the field ---
    U1 = get(findobj(hFig,'Tag','eigArrows'),'UData');
    KeyOnFig(hFig, 'rightarrow');  drawnow;
    U2 = get(findobj(hFig,'Tag','eigArrows'),'UData');
    [nPass,nFail] = chk('mode step changes the field', ~isequal(U1, U2), nPass,nFail);

    % --- quiver-size key changes arrow length (UData scales) ---
    Ua = get(findobj(hFig,'Tag','eigArrows'),'UData');
    KeyOnFig(hFig, 'rightarrow', {'shift'});  drawnow;   % Shift+Right = longer
    Ub = get(findobj(hFig,'Tag','eigArrows'),'UData');
    [nPass,nFail] = chk('quiver-size key rescales arrows', ~isequal(Ua, Ub), nPass,nFail);

    close(hFig);

    % --- non-Dirac node is rejected ---
    preL  = local_eigen_names(SurfaceFile);
    tess_eigen(SurfaceFile, 'Laplace-Beltrami', 'K', 40);
    postL = local_eigen_names(SurfaceFile);
    lboNew = setdiff(postL, preL);
    rejected = false;
    if ~isempty(lboNew)
        try
            hbad = view_eigenmodes(lboNew{1});
            rejected = isempty(hbad);   % bst_error path returns [] without a figure
            if ~isempty(hbad) && ishandle(hbad), close(hbad); end
        catch
            rejected = true;
        end
    end
    [nPass,nFail] = chk('LBO node rejected (Dirac-only)', rejected, nPass,nFail);

    fprintf('\n==== test_eigenmode_vector_viewer: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_eigenmode_vector_viewer: %d test(s) FAILED.', nFail); end
end

function names = local_eigen_names(SurfaceFile)
    names = {};
    [sSubject, ~, iSurface] = bst_get('SurfaceFile', SurfaceFile);
    if ~isempty(sSubject) && ~isempty(iSurface) ...
            && isfield(sSubject.Surface(iSurface), 'Eigen') ...
            && ~isempty(sSubject.Surface(iSurface).Eigen)
        names = {sSubject.Surface(iSurface).Eigen.FileName};
    end
end

function KeyOnFig(hFig, keyName, modifier)
    if nargin < 3, modifier = {}; end
    ev.Key = keyName; ev.Character = ''; ev.Modifier = modifier;
    cb = get(hFig, 'KeyPressFcn');
    cb(hFig, ev);
end

function [p,f] = chk(name, c, p, f)
    if c, fprintf('  PASS  %s\n', name); p=p+1; else, fprintf('  FAIL  %s\n', name); f=f+1; end
end
