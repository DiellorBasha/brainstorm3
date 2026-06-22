function out = bst_eigfilter_design_ideal(params)
% BST_EIGFILTER_DESIGN_IDEAL: brick-wall indicator g(l) = 1[lo <= l <= hi].
% params.band = [lo hi] (default [0 Inf]).

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

if nargin >= 1 && ischar(params) && strcmpi(params,'meta')
    out = struct('name','ideal','display','Ideal (brick-wall) band', ...
        'params', struct('band', struct('default',[0 Inf])), ...
        'bandpass', true, 'priorAdmissible', true);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'band') || isempty(params.band); params.band = [0 Inf]; end
lo = params.band(1); hi = params.band(2);
out = @(l) double(double(l(:)) >= lo & double(l(:)) <= hi);
end
