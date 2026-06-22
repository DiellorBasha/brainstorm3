function out = bst_eigfilter_kernel(name, params)
% BST_EIGFILTER_KERNEL: Registry mapping a kernel name string to its factory.
% USAGE:  g     = bst_eigfilter_kernel(name, params)   -> handle from bst_eigfilter_design_<name>
%         names = bst_eigfilter_kernel('list')          -> available kernel names
%         meta  = bst_eigfilter_kernel('info', name)    -> metadata for one kernel

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

if nargin < 1; error('bst_eigfilter_kernel: name required.'); end
switch lower(name)
    case 'list'
        d = dir(fullfile(fileparts(mfilename('fullpath')), 'bst_eigfilter_design_*.m'));
        pre = numel('bst_eigfilter_design_');
        out = cellfun(@(nm) nm(pre+1:end-2), {d.name}, 'UniformOutput', false);
        return;
    case 'info'
        if nargin < 2 || isempty(params) || ~ischar(params)
            error('bst_eigfilter_kernel: ''info'' requires a kernel name.');
        end
        out = feval(['bst_eigfilter_design_' lower(params)], 'meta');
        return;
end
if nargin < 2; params = struct(); end
fn = ['bst_eigfilter_design_' lower(name)];
if ~exist(fn, 'file')
    error('bst_eigfilter_kernel: unknown kernel ''%s''.', name);
end
out = feval(fn, params);
end
