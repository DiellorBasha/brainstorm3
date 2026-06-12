function test_eigenmode_vector_viewer()
% TEST_EIGENMODE_VECTOR_VIEWER  Live-figure regression for the Dirac eigenmode
% vector viewer. Requires Brainstorm running with a loaded protocol (anatomy).
%
% Reuses an existing Dirac/LBO eigen node on the canonical cortex when present
% (only computes a small K=20 basis if absent) so reruns are fast. The non-Dirac
% rejection step suppresses the modal bst_error dialog (GuiLevel=-1, restored
% immediately) so an unattended run cannot block on a dialog.
%
% Authors: Diellor Basha, 2026
    global GlobalData;
    nPass = 0; nFail = 0;

    SurfaceFile = bst_canonical_cortex(20484);

    % --- fixture: reuse-or-create a Dirac eigen node (avoid costly recompute) ---
    EigenFile = local_find_eigen(SurfaceFile, 'Dirac');
    if isempty(EigenFile)
        tess_eigen(SurfaceFile, 'Dirac', 'K', 20, 'Tau', 0.5);
        EigenFile = local_find_eigen(SurfaceFile, 'Dirac');
    end
    assert(~isempty(EigenFile), 'No Dirac eigen node available.');

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

    % --- quiver-size key changes arrow length via AutoScaleFactor (UData is the
    %     raw vector data and stays constant; quiver3's scale arg drives the factor) ---
    Sa = get(findobj(hFig,'Tag','eigArrows'),'AutoScaleFactor');
    KeyOnFig(hFig, 'rightarrow', {'shift'});  drawnow;   % Shift+Right = longer
    Sb = get(findobj(hFig,'Tag','eigArrows'),'AutoScaleFactor');
    [nPass,nFail] = chk('quiver-size key rescales arrows', Sb > Sa, nPass,nFail);

    close(hFig);

    % --- non-Dirac node is rejected; suppress the modal bst_error so an
    %     unattended run cannot block on a dialog (GuiLevel restored via onCleanup) ---
    LboFile = local_find_eigen(SurfaceFile, 'Laplace-Beltrami');
    if isempty(LboFile)
        tess_eigen(SurfaceFile, 'Laplace-Beltrami', 'K', 20);
        LboFile = local_find_eigen(SurfaceFile, 'Laplace-Beltrami');
    end
    rejected = false;
    if ~isempty(LboFile)
        gl = [];
        if isfield(GlobalData,'Program') && isfield(GlobalData.Program,'GuiLevel')
            gl = GlobalData.Program.GuiLevel;
            GlobalData.Program.GuiLevel = -1;   % bst_error skips the modal dialog
        end
        restoreGui = onCleanup(@() local_restore_gui(gl)); %#ok<NASGU>
        try
            hbad = view_eigenmodes(LboFile);
            rejected = isempty(hbad);   % rejection path returns [] without a figure
            if ~isempty(hbad) && ishandle(hbad), close(hbad); end
        catch
            rejected = true;
        end
        clear restoreGui;   % restore GuiLevel now
    end
    [nPass,nFail] = chk('LBO node rejected (Dirac-only)', rejected, nPass,nFail);

    fprintf('\n==== test_eigenmode_vector_viewer: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_eigenmode_vector_viewer: %d test(s) FAILED.', nFail); end
end

function f = local_find_eigen(SurfaceFile, Variant)
% First Eigen child of SurfaceFile whose stored Variant matches, or '' if none.
    f = '';
    [sSubject, ~, iSurface] = bst_get('SurfaceFile', SurfaceFile);
    if isempty(sSubject) || isempty(iSurface) ...
            || ~isfield(sSubject.Surface(iSurface), 'Eigen') ...
            || isempty(sSubject.Surface(iSurface).Eigen)
        return;
    end
    nodes = sSubject.Surface(iSurface).Eigen;
    for k = 1:numel(nodes)
        try
            E = load(file_fullpath(nodes(k).FileName), 'Variant');
        catch
            continue;
        end
        if isfield(E,'Variant') && strcmpi(E.Variant, Variant)
            f = nodes(k).FileName;
            return;
        end
    end
end

function local_restore_gui(gl)
    global GlobalData;
    if ~isempty(gl) && isfield(GlobalData,'Program')
        GlobalData.Program.GuiLevel = gl;
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
