function test_eigenwavelet_design()
% Frame design / evaluate / bounds for bst_eigenwavelet (data-free spectral checks).

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

lrange = [0.0, 8.0];
Lambda = linspace(lrange(1), lrange(2), 400)';
% --- itersine: tight frame -> B/A ~ 1 ---
fr = bst_eigenwavelet('Design', 'itersine', 8, lrange);
assert(numel(fr.g)>=2 && numel(fr.Centers)==numel(fr.g), 'itersine frame malformed.');
b = bst_eigenwavelet('Bounds', fr, Lambda);
assert(b.Tightness < 1.05, sprintf('itersine not tight: B/A=%.4f', b.Tightness));
% --- mexhat: finite frame bounds, M = Nf+1 members (scaling fn + wavelets) ---
fm = bst_eigenwavelet('Design', 'mexhat', 6, lrange);
assert(numel(fm.g)==7, 'mexhat frame should have Nf+1 members.');
bm = bst_eigenwavelet('Bounds', fm, Lambda);
assert(isfinite(bm.Tightness) && bm.A > 0, 'mexhat frame not a frame (A<=0).');
% --- heat: M = Nf low-pass members ---
fh = bst_eigenwavelet('Design', 'heat', 5, lrange);
assert(numel(fh.g)==5, 'heat frame should have Nf members.');
% --- Evaluate dims ---
H = bst_eigenwavelet('Evaluate', fm, Lambda);
assert(isequal(size(H), [numel(Lambda), numel(fm.g)]), 'Evaluate dims wrong.');
disp('test_eigenwavelet_design PASSED');
end
