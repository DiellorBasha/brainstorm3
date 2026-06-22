function test_drawresponse()
% Headless check that panel_eigenfilter_design('DrawResponse', ...) plots a gain curve.

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

Lambda = linspace(0.01, 5, 50)';
hFig = figure('Visible','off');
hAx  = axes('Parent', hFig);
panel_eigenfilter_design('DrawResponse', hAx, 'heat', struct('t',0.1), Lambda);
hLine = findobj(hAx, 'Type', 'line');
assert(~isempty(hLine), 'DrawResponse drew no line.');
yd = get(hLine(1), 'YData');
assert(numel(yd) == numel(Lambda), 'Gain length must equal numel(Lambda).');
assert(all(yd <= 1 + 1e-9) && all(yd >= 0), 'Heat gain must be in [0,1].');
close(hFig);
disp('test_drawresponse PASSED');
end
