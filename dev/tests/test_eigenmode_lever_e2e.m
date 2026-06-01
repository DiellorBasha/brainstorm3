function test_eigenmode_lever_e2e(SurfaceFile, ResultsFile)
% End-to-end LIVE test of the eigenmode scale lever.
% USAGE: test_eigenmode_lever_e2e(SurfaceFile, ResultsFile)
%   SurfaceFile : a cortex surface in the current protocol (eigenmodes computed
%                 on it if not already present)
%   ResultsFile : a source map (results) file defined on that surface
% Requires the Brainstorm GUI (view_surface_data + tool panels). If your session
% has no suitable protocol data, this test is data-gated -- run the four headless
% lever tests instead (they are the offline acceptance gate).
if ~brainstorm('status'); brainstorm nogui; end

% Ensure eigenmodes exist on the surface
[~, isComputed] = in_tess_eigenmodes(SurfaceFile);
if ~isComputed
    process_eigenmodes('Compute', SurfaceFile, 200, 'barycentric', true, false, false);
end

% Display the source map on the cortex
hFig = view_surface_data(SurfaceFile, ResultsFile);
assert(~isempty(hFig), 'source map failed to display');

% Bring up + populate the lever
gui_brainstorm('ShowToolTab', 'EigenModes');
panel_eigenmodes('UpdatePanel', hFig);

% Raw column (lever inactive) for comparison
panel_eigenmodes('SetActive', 0);
panel_surface('UpdateSurfaceData', hFig);
TessInfo = getappdata(hFig, 'Surface');
uRaw = TessInfo(1).Data;

% Low band => displayed column changes and does not gain energy
panel_eigenmodes('SetActive', 1);
panel_eigenmodes('SetWindowShape', 'box');
panel_eigenmodes('SetBand', 1, 15);     % fires FireModesChanged -> repaint
TessInfo = getappdata(hFig, 'Surface');
uLow = TessInfo(1).Data;
assert(~isequal(uLow, uRaw), 'low band must change the displayed column');
assert(norm(uLow) <= norm(uRaw) + 1e-9, 'low band must not increase energy');

% Toggle off restores the raw map exactly
panel_eigenmodes('SetActive', 0);
panel_surface('UpdateSurfaceData', hFig);
TessInfo = getappdata(hFig, 'Surface');
assert(max(abs(TessInfo(1).Data - uRaw)) < 1e-9, 'deactivation must restore raw map');

close(hFig);
fprintf('ALL TESTS PASSED: test_eigenmode_lever_e2e\n');
end
