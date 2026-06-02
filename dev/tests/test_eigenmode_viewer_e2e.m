function test_eigenmode_viewer_e2e(SurfaceFile)
% Live: open View Eigenmodes, step modes and superpose a band via the lever.
% USAGE: test_eigenmode_viewer_e2e(SurfaceFile)   % a cortex surface (eigenmodes auto-computed if missing)
if ~brainstorm('status'); brainstorm nogui; end
[~, isComputed] = in_tess_eigenmodes(SurfaceFile);
if ~isComputed
    process_eigenmodes('Compute', SurfaceFile, 200, 'barycentric', true, false, false);
end
hFig = view_eigenmodes(SurfaceFile);
assert(~isempty(hFig), 'viewer failed to open');
ev = getappdata(hFig, 'EigenView');
assert(~isempty(ev) && isfield(ev,'PairedGrid'), 'figure tagged EigenView with PairedGrid');

% Mode 1 (single) shown initially
TessInfo = getappdata(hFig, 'Surface');
c1 = TessInfo(1).Data;
assert(max(abs(c1 - ev.PairedGrid(:,1))) < 1e-6, 'mode 1 shown initially');

% Step to mode 2 via the lever (as the arrow keys do)
panel_eigenmodes('SetCurrentMode', 2);
TessInfo = getappdata(hFig, 'Surface');
c2 = TessInfo(1).Data;
assert(max(abs(c2 - ev.PairedGrid(:,2))) < 1e-6, 'stepping shows mode 2');

% Superpose a band [1,5] (box) -> sum of paired columns 1..5
panel_eigenmodes('SetWindowShape', 'box');
panel_eigenmodes('SetBand', 1, 5);
TessInfo = getappdata(hFig, 'Surface');
cb = TessInfo(1).Data;
assert(max(abs(cb - sum(ev.PairedGrid(:,1:5),2))) < 1e-6, 'band -> superposition');

close(hFig);
fprintf('ALL TESTS PASSED: test_eigenmode_viewer_e2e\n');
end
