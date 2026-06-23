function test_eigenwavelet_panel()
% TEST_EIGENWAVELET_PANEL: panel_eigenwavelet_options (frame design + cortical atom + steer).
%
% Validates the wavelet panel wiring:
%   - the panel builds and GetPanelContents returns a valid wavelet OPTIONS struct;
%   - the atom toggle opens a cortex figure with a patch (and a vector quiver for Dirac);
%   - the SCROLL orientation handler steers the Dirac atom (the displayed field changes);
%   - the Scale slider selects a different frame member (the atom changes);
%   - the scalar (Laplace-Beltrami) atom renders an envelope with no orientation handler.
%
% Auto-detects a real 4-block Dirac node + a Laplace-Beltrami node under Subject01; SKIPs
% the relevant block if absent.
%
% Authors: Diellor Basha, 2026

    nFail = 0;  chk = @(n,c) i_chk(n,c);
    figName = 'Eigenwavelet cortical atom';
    % clean slate: release any prior instance (e.g. from manual testing)
    try, gui_hide('EigenwaveletOptions');             catch, end %#ok<CTCH>
    try, bst_mutex('release', 'EigenwaveletOptions'); catch, end %#ok<CTCH>
    delete(findall(0,'Type','figure','Name',figName));

    efD = i_find_dirac();
    efL = i_find_variant('Laplace-Beltrami');

    % ---------- Dirac: vector atom + orientation steering ----------
    if ~isempty(efD)
        fprintf('Dirac node: %s\n', efD);
        i_clean(figName);
        [bp, nm] = panel_eigenwavelet_options('CreatePanel', efD);
        nFail = nFail + chk('panel builds (Dirac)', ~isempty(bp));
        gui_show(bp, 'JavaWindow', nm, 0,0,0); drawnow;
        ctrl = bst_get('PanelControls', 'EigenwaveletOptions');
        % GetPanelContents
        s = panel_eigenwavelet_options('GetPanelContents');
        nFail = nFail + chk('OPTIONS.Method==wavelet', strcmp(s.Method,'wavelet'));
        nFail = nFail + chk('OPTIONS has family + Nf', ismember(s.KernelName,{'itersine','mexhat','heat'}) && s.Nf>=2);
        % open the atom
        ctrl.jToggleAtom.doClick();
        hF = i_wait(figName);
        nFail = nFail + chk('atom figure opens', ~isempty(hF));
        if ~isempty(hF)
            pat = findobj(hF,'Type','patch');
            qv  = findobj(hF,'Type','quiver');
            nFail = nFail + chk('cortex patch drawn', ~isempty(pat));
            nFail = nFail + chk('vector quiver + glyph drawn (>=2)', numel(qv) >= 2);
            % steering changes the displayed field
            cd0 = get(pat(1),'FaceVertexCData');
            cb  = get(hF,'WindowScrollWheelFcn');
            nFail = nFail + chk('scroll handler installed', ~isempty(cb));
            cb(hF, struct('VerticalScrollCount', 6));   % steer
            pause(0.3); drawnow;
            pat = findobj(hF,'Type','patch'); cd1 = get(pat(1),'FaceVertexCData');
            nFail = nFail + chk('orientation steer changes the atom', ~isequal(cd0, cd1));
            % scale slider changes the member
            cdA = get(findobj(hF,'Type','patch'),'FaceVertexCData');
            ctrl.jScale.setValue(min(ctrl.jScale.getMaximum(), 4));
            cbS = java_getcb(ctrl.jScale, 'StateChangedCallback');
            if ~isempty(cbS); cbS(ctrl.jScale, []); end
            pause(0.3); drawnow;
            cdB = get(findobj(hF,'Type','patch'),'FaceVertexCData');
            nFail = nFail + chk('scale slider changes the atom', ~isequal(cdA, cdB));
        end
        i_clean(figName);
        try, gui_hide('EigenwaveletOptions'); catch, end %#ok<CTCH>
    else
        fprintf('SKIP Dirac block: no real 4-block Dirac eigen node.\n');
    end

    % ---------- Laplace-Beltrami: scalar envelope, no orientation ----------
    if ~isempty(efL)
        fprintf('LBO node: %s\n', efL);
        i_clean(figName);
        try, bst_mutex('release','EigenwaveletOptions'); catch, end %#ok<CTCH>
        [bp, nm] = panel_eigenwavelet_options('CreatePanel', efL);
        gui_show(bp, 'JavaWindow', nm, 0,0,0); drawnow;
        ctrl = bst_get('PanelControls', 'EigenwaveletOptions');
        ctrl.jToggleAtom.doClick();
        hF = i_wait(figName);
        nFail = nFail + chk('LBO atom figure opens', ~isempty(hF));
        if ~isempty(hF)
            nFail = nFail + chk('LBO patch drawn', ~isempty(findobj(hF,'Type','patch')));
            nFail = nFail + chk('LBO has no field quiver (scalar)', isempty(findobj(hF,'Type','quiver')));
        end
        i_clean(figName);
        try, gui_hide('EigenwaveletOptions'); catch, end %#ok<CTCH>
    else
        fprintf('SKIP LBO block: no Laplace-Beltrami eigen node.\n');
    end

    % ---------- Connection Laplacian: complex tangent atom, U(1) phase orientation ----------
    efC = i_find_variant('Connection Laplacian');
    if ~isempty(efC)
        fprintf('Connection node: %s\n', efC);
        % accessor: bst_operator_frame returns an orthonormal frame (stored or nxr fallback)
        EC = in_bst_eigen(efC); OC = in_bst_operator(EC.OperatorFile);
        Fr = bst_operator_frame(OC, 1);
        orthoOK = max(abs(sqrt(sum(Fr.e1.^2,2))-1))<1e-9 && max(abs(sum(Fr.e1.*Fr.e2,2)))<1e-9 ...
                  && max(abs(sum(Fr.e1.*Fr.normal,2)))<1e-9;
        nFail = nFail + chk('bst_operator_frame: orthonormal frame', orthoOK);
        i_clean(figName);
        try, bst_mutex('release','EigenwaveletOptions'); catch, end %#ok<CTCH>
        [bp, nm] = panel_eigenwavelet_options('CreatePanel', efC);
        gui_show(bp, 'JavaWindow', nm, 0,0,0); drawnow;
        ctrl = bst_get('PanelControls', 'EigenwaveletOptions');
        ctrl.jToggleAtom.doClick();
        hF = i_wait(figName);
        nFail = nFail + chk('Connection atom figure opens', ~isempty(hF));
        if ~isempty(hF)
            nFail = nFail + chk('Connection patch + tangent quiver+glyph (>=2)', ...
                ~isempty(findobj(hF,'Type','patch')) && numel(findobj(hF,'Type','quiver'))>=2);
            % U(1) phase steer: magnitude envelope (patch CData) is INVARIANT
            cd0 = get(findobj(hF,'Type','patch'),'FaceVertexCData');
            cb = get(hF,'WindowScrollWheelFcn'); cb(hF, struct('VerticalScrollCount',6));
            pause(0.3); drawnow;
            cd1 = get(findobj(hF,'Type','patch'),'FaceVertexCData');
            nFail = nFail + chk('U(1) phase steer: |psi| envelope invariant', norm(cd0-cd1)/max(norm(cd0),eps) < 1e-9);
        end
        i_clean(figName);
        try, gui_hide('EigenwaveletOptions'); catch, end %#ok<CTCH>
    else
        fprintf('SKIP Connection block: no Connection Laplacian eigen node.\n');
    end

    fprintf('\n==== test_eigenwavelet_panel: %d failed ====\n', nFail);
    if nFail > 0; error('test_eigenwavelet_panel FAILED'); end
    disp('ALL TESTS PASSED');
end

% ===== helpers =====
function r = i_chk(nm, cond)
    if cond; r = 0; fprintf('  ok   %s\n', nm); else; r = 1; fprintf('  FAIL %s\n', nm); end
end

function i_clean(figName)
    delete(findall(0,'Type','figure','Name',figName));
end

function hF = i_wait(figName)
    hF = [];
    for w = 1:30
        drawnow; pause(0.1);
        hF = findall(0,'Type','figure','Name',figName);
        if ~isempty(hF); return; end
    end
end

function ef = i_find_variant(want)
    ef = [];
    PI = bst_get('ProtocolInfo'); if isempty(PI); return; end
    d = dir(fullfile(PI.SUBJECTS, '**', 'eigen_*.mat'));
    [~, ord] = sort([d.datenum], 'descend');
    for i = ord(:)'
        rel = strrep(fullfile(d(i).folder, d(i).name), [PI.SUBJECTS filesep], '');
        try
            m = in_bst_eigen(rel, 'Variant');
            if strcmpi(m.Variant, want); ef = rel; return; end
        catch
        end
    end
end

function ef = i_find_dirac()
    ef = [];
    PI = bst_get('ProtocolInfo'); if isempty(PI); return; end
    d = dir(fullfile(PI.SUBJECTS, '**', 'eigen_*.mat'));
    [~, ord] = sort([d.datenum], 'descend');
    for i = ord(:)'
        rel = strrep(fullfile(d(i).folder, d(i).name), [PI.SUBJECTS filesep], '');
        try
            m = in_bst_eigen(rel, 'Variant');
            if ~strcmpi(m.Variant, 'Dirac'); continue; end
            E = in_bst_eigen(rel);
            hn = find(~cellfun(@isempty, E.Phi), 1);
            if ~isempty(hn) && isreal(E.Phi{hn}) && size(E.Phi{hn},1) == 4*numel(E.GlobalVertices{hn})
                ef = rel; return;
            end
        catch
        end
    end
end
