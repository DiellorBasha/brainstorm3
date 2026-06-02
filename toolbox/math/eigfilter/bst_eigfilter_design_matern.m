function out = bst_eigfilter_design_matern(params)
% BST_EIGFILTER_DESIGN_MATERN: SPDE/Matern prior g(l) = (kappa^2 + l)^(-nu).
% USAGE:  g = bst_eigfilter_design_matern(struct('kappa',1,'nu',1.5))  -> handle
%         m = bst_eigfilter_design_matern('meta')                       -> metadata

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
    out = struct('name','matern','display','Matern / SPDE (low-pass)', ...
        'params', struct('kappa', struct('default',1,'range',[0 Inf]), ...
                         'nu',    struct('default',1,'range',[0 Inf])), ...
        'bandpass', false, 'priorAdmissible', true);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
if ~isfield(params,'kappa') || isempty(params.kappa); params.kappa = 1; end
if ~isfield(params,'nu') || isempty(params.nu); params.nu = 1; end
if params.kappa < 0; error('bst_eigfilter_design_matern: kappa must be >= 0.'); end
if params.nu < 0;    error('bst_eigfilter_design_matern: nu must be >= 0 (a negative nu boosts high modes and is not a valid low-pass / prior).'); end
k = params.kappa; nu = params.nu;
out = @(l) (k^2 + double(l(:))).^(-nu);
end
