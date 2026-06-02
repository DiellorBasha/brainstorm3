function test_eigenmode_lever_panel
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
bstPanel = panel_eigenmodes('CreatePanel');
assert(~isempty(bstPanel), 'CreatePanel must return a BstPanel');
ctrl = get(bstPanel, 'sControls');
assert(~isempty(ctrl), 'panel controls must be present');
assert(isfield(ctrl,'jSliderCenter') && isfield(ctrl,'jTextWidth') ...
    && isfield(ctrl,'jCheckActive') && isfield(ctrl,'jLabelReadout'), ...
    'expected center slider + width field + active + readout');
assert(isfield(ctrl,'jRadioBox') && isfield(ctrl,'jRadioTaper') && isfield(ctrl,'jRadioGauss'), ...
    'expected Box/Taper/Gauss shape radios');
fprintf('ALL TESTS PASSED: test_eigenmode_lever_panel\n');
end
