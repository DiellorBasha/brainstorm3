function out = bst_eigfilter_design_itersine(params)
% BST_EIGFILTER_DESIGN_ITERSINE: One member of an itersine (half-cosine) TIGHT frame.
% The Nf members' squared responses sum to ~constant over the interior (a tight frame),
% so a bank of these is a robust frame (B/A -> 1). Used by the Dynamics "Design tight
% frame" generate; parameters (member/Nf/lmax) are set programmatically, not slider-tuned.
% USAGE:  g = bst_eigfilter_design_itersine(struct('member',ii,'Nf',N,'lmax',L))  -> handle
%         m = bst_eigfilter_design_itersine('meta')                                -> metadata
%
% Authors: Diellor Basha, 2026

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

if nargin >= 1 && ischar(params) && strcmpi(params,'meta')
    out = struct('name','itersine','display','Itersine (tight-frame member)', ...
        'params', struct('member',struct('default',1,'range',[1 Inf]), ...
                         'Nf',struct('default',6,'range',[2 Inf]), ...
                         'lmax',struct('default',1,'range',[0 Inf])), ...
        'bandpass', true, 'priorAdmissible', false);
    return;
end
if nargin < 1 || isempty(params); params = struct(); end
member = 1;  if isfield(params,'member') && ~isempty(params.member); member = params.member; end
Nf     = 6;  if isfield(params,'Nf')     && ~isempty(params.Nf);     Nf     = params.Nf;     end
lmax   = 1;  if isfield(params,'lmax')   && ~isempty(params.lmax);   lmax   = params.lmax;   end
if ~(lmax > 0); error('bst_eigfilter_design_itersine: lmax must be > 0.'); end
overlap = 2;
scale   = lmax / (Nf - overlap + 1) * overlap;
kf      = @(x) sin(0.5*pi*(cos(pi*x)).^2) .* (x >= -0.5 & x <= 0.5);
out     = @(l) kf(double(l(:))/scale - (member - overlap/2)/overlap) ./ sqrt(overlap) .* sqrt(2);
end
