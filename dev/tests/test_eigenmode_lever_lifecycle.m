function test_eigenmode_lever_lifecycle
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));
c = panel_eigenmodes('ClassifyContext', struct('isEigenView',false,'hasSourceModes',false));
assert(strcmp(c.kind,'none') && ~c.selectEnabled && ~c.activeEnabled, 'none -> all off');
c = panel_eigenmodes('ClassifyContext', struct('isEigenView',true,'hasSourceModes',false));
assert(strcmp(c.kind,'view') && c.selectEnabled && ~c.activeEnabled, 'view -> select on, Active off');
c = panel_eigenmodes('ClassifyContext', struct('isEigenView',false,'hasSourceModes',true));
assert(strcmp(c.kind,'source') && c.selectEnabled && c.activeEnabled, 'source -> select + Active on');
fprintf('ALL TESTS PASSED: test_eigenmode_lever_lifecycle\n');
end
