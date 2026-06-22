function h = bst_eigfilter_evaluate(g, lambdas)
% BST_EIGFILTER_EVALUATE: Evaluate a kernel handle (or cell bank) on eigenvalues.
%   h = [K x 1] for a single handle; h = [K x Nf] for a cell array of Nf handles.

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

lambdas = double(lambdas(:));
if iscell(g)
    h = zeros(numel(lambdas), numel(g));
    for ii = 1:numel(g); h(:,ii) = g{ii}(lambdas); end
else
    h = g(lambdas); h = h(:);
end
end
