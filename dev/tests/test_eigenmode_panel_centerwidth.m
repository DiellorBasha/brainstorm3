function test_eigenmode_panel_centerwidth
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));
K = 600;
[lo,hi] = panel_eigenmodes('BandFromCenterWidth', 42, 13, K);
assert(lo==29 && hi==55, 'center 42 width 13 -> [29 55]');
[lo,hi] = panel_eigenmodes('BandFromCenterWidth', 42, 0, K);
assert(lo==42 && hi==42, 'width 0 -> single [42 42]');
[lo,hi] = panel_eigenmodes('BandFromCenterWidth', 3, 10, K);
assert(lo==1 && hi==13, 'clamp low edge to 1');
[lo,hi] = panel_eigenmodes('BandFromCenterWidth', 598, 10, K);
assert(lo==588 && hi==600, 'clamp high edge to K');
[lo,hi] = panel_eigenmodes('BandFromCenterWidth', 300, 9999, K);
assert(lo==1 && hi==600, 'width >= K -> full [1 K]');
fprintf('ALL TESTS PASSED: test_eigenmode_panel_centerwidth\n');
end
