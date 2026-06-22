function test_panel_eigenfilter_options()
% Build panel_eigenfilter_options from a minimal eigen_ fixture, read it back,
% assert the OPTIONS contract that feeds bst_eigen. Self-contained (no protocol node).

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

% ----- build a minimal valid eigen_ fixture on disk -----
nV = 100; K = 30;
Lam = linspace(0, 5, K)';
E = struct();
E.Comment        = 'eigen_test_fixture';
E.Variant        = 'Laplace-Beltrami';
E.Lambda         = {Lam, Lam};
E.Phi            = {zeros(nV, K), zeros(nV, K)};   % required by in_bst_eigen guard
E.nModes         = K;
E.GlobalVertices = {(1:nV)', (1:nV)'};
E.GlobalFaces    = {[], []};
E.OperatorFile   = '';
E.ParentSurface  = '';
fixture = fullfile(tempdir, 'eigen_test_fixture.mat');
save(fixture, '-struct', 'E', '-v6');
cleaner = onCleanup(@() delete(fixture));

% ----- drive the panel -----
[bstPanel, panelName] = panel_eigenfilter_options('CreatePanel', fixture);
assert(~isempty(bstPanel), 'CreatePanel returned empty.');
gui_show(bstPanel, 'JavaWindow', panelName, 0, 0, 0);
drawnow;
s = panel_eigenfilter_options('GetPanelContents');
gui_hide(panelName);

% ----- assert the OPTIONS contract -----
assert(strcmp(s.Method, 'filter'),                'Method must be ''filter''.');
assert(ischar(s.EigenFile) && ~isempty(s.EigenFile), 'EigenFile must be set.');
assert(ischar(s.KernelName) && ~isempty(s.KernelName), 'KernelName must be set.');
assert(isstruct(s.KernelParams),                  'KernelParams must be a struct.');
assert(ischar(s.Comment),                         'Comment must be a char.');
disp('test_panel_eigenfilter_options PASSED');
end
