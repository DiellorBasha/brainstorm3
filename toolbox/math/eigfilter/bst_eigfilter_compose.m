function g = bst_eigfilter_compose(varargin)
% BST_EIGFILTER_COMPOSE: Serial composition of kernels = pointwise product of gains.
%   g = bst_eigfilter_compose(g1, g2, ...)  ->  @(l) g1(l).*g2(l).*...

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

gs = varargin;
for ii = 1:numel(gs)
    if ~isa(gs{ii}, 'function_handle')
        error('bst_eigfilter_compose: all arguments must be function handles.');
    end
end
g = @(l) i_prod(gs, l);
end
function y = i_prod(gs, l)
l = double(l(:));
y = ones(size(l));
for ii = 1:numel(gs); y = y .* gs{ii}(l); end
end
