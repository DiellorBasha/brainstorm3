function out = bst_eigfilter_design_tikhonov(params)
% BST_EIGFILTER_DESIGN_TIKHONOV: low-pass g(l) = 1/(1+beta*l).

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
    out = struct('name','tikhonov','display','Tikhonov / membrane (low-pass)', ...
        'params', struct('beta', struct('default',1,'range',[0 Inf])), ...
        'bandpass', false, 'priorAdmissible', true);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'beta') || isempty(params.beta); params.beta = 1; end
if params.beta < 0; error('bst_eigfilter_design_tikhonov: beta must be >= 0.'); end
b = params.beta;
out = @(l) 1 ./ (1 + b * double(l(:)));
end
