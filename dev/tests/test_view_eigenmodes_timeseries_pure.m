function test_view_eigenmodes_timeseries_pure
% Verify the band->trace selector: each paired rank in [kLo,kHi] yields its
% left column then its right column as two traces; single-component data
% yields one unlabelled-hemisphere trace per rank.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% --- Paired case: 2 components, ranks 1..2 each ---
% Raw columns:  1=L/r1, 2=L/r2, 3=R/r1, 4=R/r2
Component = [1;1;2;2];
CompRank  = [1;2;1;2];

[iRows, Labels, Hemi] = view_eigenmodes_timeseries('GetBandTraces', Component, CompRank, 1, 2);
assert(isequal(iRows(:)', [1 3 2 4]), 'Order must be L then R, per rank ascending.');
assert(isequal(Labels, {'Mode 1 L','Mode 1 R','Mode 2 L','Mode 2 R'}), 'Paired labels wrong.');
assert(isequal(Hemi(:)', [1 2 1 2]), 'Hemisphere ids wrong.');

% Sub-band selects a single rank -> 2 traces
[iRows2, Labels2] = view_eigenmodes_timeseries('GetBandTraces', Component, CompRank, 2, 2);
assert(isequal(iRows2(:)', [2 4]), 'Single-rank band must give that rank''s L,R columns.');
assert(isequal(Labels2, {'Mode 2 L','Mode 2 R'}), 'Single-rank labels wrong.');

% --- Single-component case: labels drop the hemisphere suffix ---
Comp1 = [1;1;1];
Rank1 = [1;2;3];
[iRows3, Labels3, Hemi3] = view_eigenmodes_timeseries('GetBandTraces', Comp1, Rank1, 1, 2);
assert(isequal(iRows3(:)', [1 2]), 'Single-component rows wrong.');
assert(isequal(Labels3, {'Mode 1','Mode 2'}), 'Single-component labels must omit L/R.');
assert(isequal(Hemi3(:)', [0 0]), 'Single-component hemisphere id must be 0.');

% --- Asymmetric rank (rank 2 only on left) ---
CompA = [1;1;2];
RankA = [1;2;1];
[iRowsA, LabelsA] = view_eigenmodes_timeseries('GetBandTraces', CompA, RankA, 1, 2);
assert(isequal(iRowsA(:)', [1 3 2]), 'Asymmetric ordering wrong.');
assert(isequal(LabelsA, {'Mode 1 L','Mode 1 R','Mode 2 L'}), 'Asymmetric labels wrong.');

% --- Color helper: hemisphere -> distinct RGB cell array ---
colors = view_eigenmodes_timeseries('HemiColors', [1 2 0]);
assert(iscell(colors) && numel(colors) == 3, 'HemiColors must return a 1x3 cell.');
assert(~isequal(colors{1}, colors{2}), 'Left and right must differ in color.');

disp('ALL TESTS PASSED');
end
