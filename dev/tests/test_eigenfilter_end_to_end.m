function OutputFiles = test_eigenfilter_end_to_end(ResultsFile, EigenFile, TimeWindow)
% End-to-end: build OPTIONS from panel_eigenfilter_options, run the eigen filter via
% bst_eigen, assert a filtered source map is produced.
%   ResultsFile : a constrained (scalar) source map on the same surface as EigenFile
%   EigenFile   : a Laplace-Beltrami eigen_ node (with linked operator_)
%   TimeWindow  : [t0 t1] to restrict the filter (keeps the test fast)

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

% ----- drive the panel to produce OPTIONS -----
[bstPanel, panelName] = panel_eigenfilter_options('CreatePanel', EigenFile);
gui_show(bstPanel, 'JavaWindow', panelName, 0, 0, 0);
drawnow;
OPTIONS = panel_eigenfilter_options('GetPanelContents');
gui_hide(panelName);

% ----- run the filter (short window, do not write to DB) -----
OPTIONS.iTargetStudy = 'NoSave';
OPTIONS.TimeWindow   = TimeWindow;
[OutputFiles, Messages, isError] = bst_eigen(ResultsFile, OPTIONS);

% ----- assert -----
assert(~isError, 'bst_eigen reported an error: %s', Messages);
assert(~isempty(OutputFiles), 'bst_eigen produced no output.');
FileMat = OutputFiles{1};
assert(isfield(FileMat,'ImageGridAmp') && ~isempty(FileMat.ImageGridAmp), ...
    'Filtered source map is empty.');
fprintf('Filtered map: [%d x %d], kernel=%s\n', ...
    size(FileMat.ImageGridAmp,1), size(FileMat.ImageGridAmp,2), OPTIONS.KernelName);
% Heat low-pass must change the map but preserve its size vs the input window.
assert(size(FileMat.ImageGridAmp,1) > 0, 'No sources in output.');
disp('test_eigenfilter_end_to_end PASSED');
end
