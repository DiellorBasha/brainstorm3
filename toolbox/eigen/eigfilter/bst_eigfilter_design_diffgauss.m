function out = bst_eigfilter_design_diffgauss(params)
% BST_EIGFILTER_DESIGN_DIFFGAUSS: Difference-of-Gaussians band-pass
% g(l) = exp(-t1*l) - exp(-t2*l). Non-negative when t1 < t2; zero at l=0.
% USAGE:  g = bst_eigfilter_design_diffgauss(struct('t1',0.1,'t2',0.4))  -> handle
%         m = bst_eigfilter_design_diffgauss('meta')                      -> metadata

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
    out = struct('name','diffgauss','display','Difference of Gaussians (band-pass)', ...
        'params', struct('t1', struct('default',0.01,'range',[0 Inf]), ...
                         't2', struct('default',0.04,'range',[0 Inf])), ...
        'bandpass', true, 'priorAdmissible', false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'t1') || isempty(params.t1); params.t1 = 0.01; end
if ~isfield(params,'t2') || isempty(params.t2); params.t2 = 0.04; end
if params.t1 < 0; error('bst_eigfilter_design_diffgauss: t1 must be >= 0.'); end
if params.t1 >= params.t2
    error('bst_eigfilter_design_diffgauss: require t1 < t2 (got t1=%g, t2=%g).', params.t1, params.t2);
end
t1 = params.t1; t2 = params.t2;
out = @(l) exp(-t1 * double(l(:))) - exp(-t2 * double(l(:)));
end
