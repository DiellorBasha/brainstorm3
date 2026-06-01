function test_eigenmode_lever_state
% State logic: init, clamp, coupled center/band, isActive default.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));
global GlobalData; %#ok<GVMIS>

% Reset state for a 200-mode surface
panel_eigenmodes('ResetState', '/fake/surf.mat', 200);
st = GlobalData.UserModes;
assert(strcmp(st.SurfaceFile, '/fake/surf.mat'), 'surface stored');
assert(st.nModes == 200, 'K stored');
assert(st.isActive == 0, 'lever starts inactive');

% SetBand clamps and recenters iCurrentMode to band center
panel_eigenmodes('SetBand', 30, 55);
st = GlobalData.UserModes;
assert(isequal(st.Band, [30 55]), 'band stored');
assert(st.iCurrentMode == round((30+55)/2), 'center coupled to band midpoint');
assert(isequal(size(st.Weights), [1 200]), 'weights recomputed [1xK]');
assert(sum(st.Weights) == 26, 'box default keeps 26 modes');

% SetBand clamps out-of-range
panel_eigenmodes('SetBand', -10, 9999);
st = GlobalData.UserModes;
assert(isequal(st.Band, [1 200]), 'band clamped to [1,K]');

% SetCurrentMode (coupled): slides the band, preserving its width
panel_eigenmodes('SetBand', 30, 50);   % 21-mode band, span 20, center 40
panel_eigenmodes('SetCurrentMode', 100);
st = GlobalData.UserModes;
assert(st.iCurrentMode == 100, 'center moved');
assert((st.Band(2) - st.Band(1)) == 20, 'band width preserved when sliding center');
assert(st.Band(1) == 90 && st.Band(2) == 110, 'band slid to recentre on 100');

% SetWindowShape 'single' collapses band to the center
panel_eigenmodes('SetWindowShape', 'single');
st = GlobalData.UserModes;
assert(strcmp(st.WindowShape, 'single'), 'shape stored');
assert(sum(st.Weights) == 1 && st.Weights(100) == 1, 'single -> delta at center');

% SetActive toggles the flag
panel_eigenmodes('SetActive', true);
assert(GlobalData.UserModes.isActive == 1, 'SetActive true');
panel_eigenmodes('SetActive', false);
assert(GlobalData.UserModes.isActive == 0, 'SetActive false');

% SetCurrentMode near the upper boundary clamps the displayed band...
panel_eigenmodes('SetBand', 30, 50);          % span 20
panel_eigenmodes('SetCurrentMode', 198);      % K = 200
st = GlobalData.UserModes;
assert(st.iCurrentMode == 198, 'center at 198');
assert(st.Band(2) == 200, 'upper edge clamped to K');
assert(st.Band(1) == 188, 'lower edge slides in (span shortened at boundary)');

% ...but the canonical span is preserved: sliding back to the interior recovers width 20
panel_eigenmodes('SetCurrentMode', 100);
st = GlobalData.UserModes;
assert((st.Band(2) - st.Band(1)) == 20, 'band width recovers after boundary round-trip');

fprintf('ALL TESTS PASSED: test_eigenmode_lever_state\n');
end
