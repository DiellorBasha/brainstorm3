function out = bst_eigfilter_design_mexhat(params)
% BST_EIGFILTER_DESIGN_MEXHAT: Mexican-hat band-pass g(l) = (t*l).*exp(-t*l).
% Vector t returns a cell-array filterbank. Zero at l=0 (not prior-admissible).
% USAGE:  g  = bst_eigfilter_design_mexhat(struct('t',0.1))          -> handle
%         gb = bst_eigfilter_design_mexhat(struct('t',[0.05 0.1 0.2]))-> cell bank
%         m  = bst_eigfilter_design_mexhat('meta')                    -> metadata

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
    out = struct('name','mexhat','display','Mexican hat (band-pass)', ...
        'params', struct('t', struct('default',0.01,'range',[0 Inf])), ...
        'bandpass', true, 'priorAdmissible', false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'t') || isempty(params.t); params.t = 0.01; end
t = params.t;
if numel(t) > 1
    out = cell(numel(t),1);
    for ii = 1:numel(t); out{ii} = i_mh(t(ii)); end
else
    out = i_mh(t);
end
end
function g = i_mh(t)
if t < 0; error('bst_eigfilter_design_mexhat: t must be >= 0.'); end
g = @(l) (t * double(l(:))) .* exp(-t * double(l(:)));
end
