function out = bst_eigfilter_design_heat(params)
% BST_EIGFILTER_DESIGN_HEAT: Heat / diffusion low-pass kernel g(l) = exp(-t*l).
% USAGE:  g = bst_eigfilter_design_heat(struct('t',0.01))   -> handle (cell if t is a vector)
%         m = bst_eigfilter_design_heat('meta')              -> metadata
% Optional params.lmax normalizes the spectrum: exp(-t*l/lmax).

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

if nargin >= 1 && ischar(params) && strcmpi(params, 'meta')
    out = struct('name','heat', 'display','Heat / diffusion (low-pass)', ...
        'params', struct('t', struct('default',0.01,'range',[0 Inf])), ...
        'bandpass', false, 'priorAdmissible', true);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'t') || isempty(params.t); params.t = 0.01; end
if ~isfield(params,'lmax'); params.lmax = []; end
sc = 1; if ~isempty(params.lmax) && params.lmax > 0; sc = 1/params.lmax; end
t = params.t;
if numel(t) > 1
    out = cell(numel(t),1);
    for ii = 1:numel(t); out{ii} = i_make(t(ii), sc); end
else
    out = i_make(t, sc);
end
end
function g = i_make(t, sc)
if t < 0; error('bst_eigfilter_design_heat: t must be >= 0 (got %g).', t); end
g = @(l) exp(-t * sc * double(l(:)));
end
