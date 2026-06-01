function test_eigenmode_lever_panel
% Smoke test: the panel builds and exposes the expected controls.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

bstPanel = panel_eigenmodes('CreatePanel');
assert(~isempty(bstPanel), 'CreatePanel must return a BstPanel');

% Read controls directly off the returned panel object (no GUI registration needed)
ctrl = get(bstPanel, 'sControls');
assert(~isempty(ctrl), 'panel controls must be present');
assert(isfield(ctrl, 'jCheckActive') && isfield(ctrl, 'jSliderLo') ...
       && isfield(ctrl, 'jSliderHi') && isfield(ctrl, 'jLabelReadout'), ...
       'expected controls present (jCheckActive, jSliderLo, jSliderHi, jLabelReadout)');

fprintf('ALL TESTS PASSED: test_eigenmode_lever_panel\n');
end
