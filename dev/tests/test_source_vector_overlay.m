function test_source_vector_overlay()
% TEST_SOURCE_VECTOR_OVERLAY  Live-figure regression for the quiver overlay.
%
% Reproduces the NextPlot='replace' axes-reset bug: drawing the quiver with a
% high-level plot function on the Brainstorm 3D axes (NextPlot='replace') ran
% newplot and RESET the axes -- wiping its 'Axes3D' tag, deleting the cortex
% patch, and restoring the default white/grid look. The overlay must instead
% leave the cortex axes intact.
%
% Requires Brainstorm running + the alpha Dirac vector source node (created by
% make_alpha_dirac_source if absent). Opens a 3D figure, toggles the overlay,
% asserts the axes/cortex survive and the overlay tracks the time cursor, then
% closes the figure.
%
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;
    cond = 'Subject01/S01_AEF_20131218_01_notch';

    % --- find (or create) the unconstrained Dirac source node ---
    sStudy = bst_get('StudyWithCondition', cond);
    iR = find(~cellfun(@isempty, regexp({sStudy.Result.FileName}, 'results_dirac_vec', 'once')), 1);
    if isempty(iR)
        ResultsFile = make_alpha_dirac_source();
    else
        ResultsFile = sStudy.Result(iR).FileName;
    end
    RM = in_bst_results(ResultsFile, 0);
    SurfaceFile = RM.SurfaceFile;

    % --- open a FRESH figure (avoid reusing a contaminated one) ---
    close(findobj(0, 'type', 'figure', 'Tag', '3DViz'));
    hFig = view_surface_data(SurfaceFile, ResultsFile, [], 'NewFigure');
    drawnow;
    TI = getappdata(hFig, 'Surface');
    iTess = find(arrayfun(@(t) ~isempty(t.DataSource) && strcmpi(t.DataSource.Type,'Source'), TI), 1);
    hPatch = TI(iTess).hPatch;

    % --- toggle the overlay ON ---
    figure_3d('SetShowSourceVectors', hFig, iTess, true); drawnow;
    hAx3D = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    hQ    = findobj(hFig, 'Tag', 'SourceVectors');

    % The cortex axes must NOT have been reset by drawing the overlay:
    [nPass,nFail] = chk('Axes3D tag survives overlay draw', ~isempty(hAx3D), nPass,nFail);
    [nPass,nFail] = chk('cortex patch survives overlay draw', any(ishandle(hPatch)), nPass,nFail);
    [nPass,nFail] = chk('quiver created', ~isempty(hQ), nPass,nFail);
    if ~isempty(hQ) && ~isempty(hAx3D)
        [nPass,nFail] = chk('quiver parented to Axes3D', isequal(get(hQ(1),'Parent'), hAx3D(1)), nPass,nFail);
        [nPass,nFail] = chk('overlay draws arrows', numel(get(hQ(1),'UData')) > 0, nPass,nFail);
    end

    % --- time-step: the field must update ---
    panel_time('SetCurrentTime', 21.0); drawnow; U1 = get(findobj(hFig,'Tag','SourceVectors'),'UData');
    panel_time('SetCurrentTime', 23.0); drawnow; U2 = get(findobj(hFig,'Tag','SourceVectors'),'UData');
    [nPass,nFail] = chk('overlay updates on time step', ~isequal(U1, U2), nPass,nFail);

    % --- threshold matches the colormap: gating uses DataLimitValue+DataThreshold,
    %     the SAME pair the colormap thresholds. Drive a data-scaled limit directly
    %     (independent of the session's shared colormap max) and redraw the overlay. ---
    panel_time('SetCurrentTime', 22.6); drawnow;
    TI2 = getappdata(hFig,'Surface');
    TI2(iTess).DataLimitValue = TI2(iTess).DataMinMax;   % real data range, not a custom max
    TI2(iTess).DataThreshold  = 0;
    setappdata(hFig,'Surface',TI2);
    figure_3d('PlotSourceVectors', hFig, iTess); drawnow;
    nFull = numel(get(findobj(hFig,'Tag','SourceVectors'),'UData'));
    TI2 = getappdata(hFig,'Surface'); TI2(iTess).DataThreshold = 0.5; setappdata(hFig,'Surface',TI2);
    figure_3d('PlotSourceVectors', hFig, iTess); drawnow;
    nThr = numel(get(findobj(hFig,'Tag','SourceVectors'),'UData'));
    [nPass,nFail] = chk('threshold reduces arrows (colormap match)', nThr < nFull && nThr > 0, nPass,nFail);
    TI2 = getappdata(hFig,'Surface'); TI2(iTess).DataThreshold = 0; setappdata(hFig,'Surface',TI2);
    figure_3d('PlotSourceVectors', hFig, iTess); drawnow;

    % --- density keys: '[' fewer, ']' more (the AdjustSourceVectorDensity hook) ---
    nA = numel(get(findobj(hFig,'Tag','SourceVectors'),'UData'));
    figure_3d('AdjustSourceVectorDensity', hFig, 1/3);
    nLess = numel(get(findobj(hFig,'Tag','SourceVectors'),'UData'));
    [nPass,nFail] = chk('[ reduces arrow count', nLess < nA, nPass,nFail);
    figure_3d('AdjustSourceVectorDensity', hFig, 3);
    nMore = numel(get(findobj(hFig,'Tag','SourceVectors'),'UData'));
    [nPass,nFail] = chk('] increases arrow count', nMore > nLess, nPass,nFail);

    % --- color keys: 'b' blue, 'k' black, and the color survives a redraw ---
    figure_3d('SetSourceVectorColor', hFig, [0 0 1]);
    [nPass,nFail] = chk('b -> blue arrows', isequal(get(findobj(hFig,'Tag','SourceVectors'),'Color'), [0 0 1]), nPass,nFail);
    panel_time('SetCurrentTime', 22.0); drawnow;   % redraw must preserve color
    [nPass,nFail] = chk('color persists across time step', isequal(get(findobj(hFig,'Tag','SourceVectors'),'Color'), [0 0 1]), nPass,nFail);
    figure_3d('SetSourceVectorColor', hFig, [0 0 0]);
    [nPass,nFail] = chk('k -> black arrows', isequal(get(findobj(hFig,'Tag','SourceVectors'),'Color'), [0 0 0]), nPass,nFail);

    % --- toggle OFF: the quiver must be removed ---
    figure_3d('SetShowSourceVectors', hFig, iTess, false); drawnow;
    [nPass,nFail] = chk('overlay removed on toggle off', isempty(findobj(hFig,'Tag','SourceVectors')), nPass,nFail);

    close(hFig);
    fprintf('\n==== test_source_vector_overlay: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_source_vector_overlay: %d test(s) FAILED.', nFail); end
end

function [p,f] = chk(name, c, p, f)
    if c, fprintf('  PASS  %s\n', name); p=p+1; else, fprintf('  FAIL  %s\n', name); f=f+1; end
end
