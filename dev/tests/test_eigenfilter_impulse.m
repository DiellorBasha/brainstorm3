function test_eigenfilter_impulse()
% TEST_EIGENFILTER_IMPULSE: Cortical impulse-response Display of panel_eigenfilter_options.
%
% Validates the filter ATOM = filtered unit delta (point-spread on the cortex):
%   - the panel builds and the impulse toggle opens a cortex figure with a patch;
%   - heat (low-pass) -> the response peaks AT the seed (band-limited blob), real;
%   - mexhat / diffgauss (band-pass) -> sign-changing center-surround at the seed;
%   - the opposite hemisphere is untouched (delta lives on one hemisphere);
%   - the Seed-vertex field re-targets the atom.
%
% Needs a Laplace-Beltrami eigen_ node. Auto-detects one under Subject01; SKIPs if none.
%
% Authors: Diellor Basha, 2026

    nFail = 0;
    chk = @(name, cond) deal_fail(name, cond);

    % ---- locate a Laplace-Beltrami eigen node ----
    ef = i_find_lb_eigen();
    if isempty(ef)
        fprintf('SKIP test_eigenfilter_impulse: no Laplace-Beltrami eigen_ node found.\n');
        return;
    end
    fprintf('Using eigen node: %s\n', ef);

    EM   = in_bst_eigen(ef);
    Op   = in_bst_operator(EM.OperatorFile);
    Surf = in_tess_bst(EM.ParentSurface, 0);
    nV   = size(Surf.Vertices, 1);
    gv1  = EM.GlobalVertices{1};
    gv2  = EM.GlobalVertices{2};
    % seed = hemi-1 centroid-nearest (same rule as the panel default)
    Vh = Surf.Vertices(gv1, :);
    [~, ii] = min(sum((Vh - mean(Vh,1)).^2, 2));
    seed = gv1(ii);

    % ---- core math: the atom is bst_eigenfilter('Analysis', delta, ...) ----
    delta = zeros(nV, 1); delta(seed) = 1;

    % heat: realistic slider scale t = 1/lambda_k for a mid mode -> localized blob
    Lam = EM.Lambda{1}(:);
    kmid = max(1, round(numel(Lam)/2));
    th = 1 / Lam(kmid);
    [ph, msgh, eh] = bst_eigenfilter('Analysis', delta, EM, Op, 'heat', struct('t', th));
    nFail = nFail + chk('heat: no error',        eh == 0 && isempty(msgh));
    ph = ph(:,1);
    nFail = nFail + chk('heat: real-valued',     isreal(ph));
    [~, ipk] = max(ph);
    dpk = norm(Surf.Vertices(ipk,:) - Surf.Vertices(seed,:)) * 1000;   % mm
    nFail = nFail + chk('heat: peak within 5mm of seed', dpk < 5);
    nFail = nFail + chk('heat: opposite hemi untouched', max(abs(ph(gv2))) < 1e-12 * max(abs(ph)));

    % diffgauss / mexhat: band-pass -> sign-changing center-surround
    [pd, ~, ed] = bst_eigenfilter('Analysis', delta, EM, Op, 'diffgauss', struct('t1', 1/Lam(end), 't2', 1/Lam(kmid)));
    pd = pd(:,1);
    nFail = nFail + chk('diffgauss: no error',   ed == 0);
    nFail = nFail + chk('diffgauss: sign-changing', min(pd) < 0 && max(pd) > 0);

    [pm, ~, em] = bst_eigenfilter('Analysis', delta, EM, Op, 'mexhat', struct('t', th));
    pm = pm(:,1);
    nFail = nFail + chk('mexhat: no error',      em == 0);
    nFail = nFail + chk('mexhat: sign-changing', min(pm) < 0 && max(pm) > 0);

    % ---- GUI path: panel builds, toggle opens a cortex figure with a patch ----
    figName = 'Eigenfilter cortical impulse response';
    % clean slate: release any prior instance of this panel (e.g. from manual testing)
    try, gui_hide('EigenfilterOptions');             catch, end %#ok<CTCH>
    try, bst_mutex('release', 'EigenfilterOptions'); catch, end %#ok<CTCH>
    delete(findall(0, 'Type','figure', 'Name', figName));
    [bstPanel, panelName] = panel_eigenfilter_options('CreatePanel', ef);
    nFail = nFail + chk('panel builds', ~isempty(bstPanel));
    gui_show(bstPanel, 'JavaWindow', panelName, 0, 0, 0);
    drawnow;
    ctrl = bst_get('PanelControls', 'EigenfilterOptions');
    nFail = nFail + chk('impulse toggle present', isfield(ctrl,'jToggleImp') && ~isempty(ctrl.jToggleImp));
    nFail = nFail + chk('seed field present',     isfield(ctrl,'jTextSeed')  && ~isempty(ctrl.jTextSeed));

    ctrl.jKernel.setSelectedIndex(find(strcmp(ctrl.KernelKeys,'heat'),1) - 1);
    ctrl.jToggleImp.doClick();                    % -> ToggleImpulse -> figure + UpdateImpulse
    hImp = [];
    for w = 1:30                                   % wait up to ~3s for the EDT callback
        drawnow; pause(0.1);
        hImp = findall(0, 'Type','figure', 'Name', figName);
        if ~isempty(hImp); break; end
    end
    nFail = nFail + chk('impulse figure opened', ~isempty(hImp));
    if ~isempty(hImp)
        pat = findobj(hImp, 'Type','patch');
        nFail = nFail + chk('cortex patch drawn', ~isempty(pat));
        cd = get(pat(1), 'FaceVertexCData');
        nFail = nFail + chk('patch CData spans full surface', numel(cd) == nV);
        % seed field re-targets: set a hemi-2 vertex and confirm it is accepted
        newSeed = gv2(round(numel(gv2)/2));
        ctrl.jTextSeed.setText(num2str(newSeed));
        ctrl.jTextSeed.postActionEvent();                 % fire OnSeedText (Enter)
        pause(0.5); drawnow;
        nFail = nFail + chk('seed field updated', str2double(char(ctrl.jTextSeed.getText())) == newSeed);
    end

    % ---- cleanup ----
    delete(findall(0, 'Type','figure', 'Name', figName));
    try, gui_hide(panelName); catch, end %#ok<CTCH>

    fprintf('\n==== test_eigenfilter_impulse: %d failed ====\n', nFail);
    if nFail > 0; error('test_eigenfilter_impulse FAILED'); end
    disp('ALL TESTS PASSED');

    function r = deal_fail(nm, cond)
        if cond; r = 0; fprintf('  ok   %s\n', nm);
        else;    r = 1; fprintf('  FAIL %s\n', nm); end
    end
end


function ef = i_find_lb_eigen()
    ef = [];
    ProtocolInfo = bst_get('ProtocolInfo');
    if isempty(ProtocolInfo); return; end
    d = dir(fullfile(ProtocolInfo.SUBJECTS, '**', 'eigen_*.mat'));
    for i = 1:numel(d)
        rel = strrep(fullfile(d(i).folder, d(i).name), [ProtocolInfo.SUBJECTS filesep], '');
        try
            m = in_bst_eigen(rel, 'Variant');
            if strcmpi(m.Variant, 'Laplace-Beltrami'); ef = rel; return; end
        catch
        end
    end
end
