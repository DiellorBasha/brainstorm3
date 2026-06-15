function test_overlay_hook()
% A figure's CustomOverlayFcn appdata must be called by panel_surface('UpdateSurfaceData').
% Authors: Diellor Basha, 2026
    nFail = 0;
    hFig = figure('Visible','off');
    % one surface with no data source -> UpdateSurfaceData takes the early-return path (no-op)
    setappdata(hFig, 'Surface', struct('Name','x', 'DataSource',struct('Type',''), 'hPatch',[]));
    assignin('base','OVH', 0);
    setappdata(hFig, 'CustomOverlayFcn', @(h) assignin('base','OVH', evalin('base','OVH')+1));
    panel_surface('UpdateSurfaceData', hFig);
    nFail = nFail + chk('overlay fcn fired once', evalin('base','OVH')==1);
    rmappdata(hFig, 'CustomOverlayFcn');
    panel_surface('UpdateSurfaceData', hFig);
    nFail = nFail + chk('not fired after removal', evalin('base','OVH')==1);
    close(hFig);
    fprintf('\n==== test_overlay_hook: %d failed ====\n', nFail);
    if nFail > 0, error('test_overlay_hook FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
