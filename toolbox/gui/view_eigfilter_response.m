function hFig = view_eigfilter_response(g, lambdas, titleStr)
% VIEW_EIGFILTER_RESPONSE: Small companion plot of a spectral kernel g(lambda).
%
% USAGE:  view_eigfilter_response(g, lambdas, titleStr)   % create/update (reused by tag)
%         view_eigfilter_response('close')                % close it
%
% g        : kernel function handle @(l)
% lambdas  : eigenvalues to mark on the x-axis (real lambda_k)
% titleStr : kernel display name
%
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

TAG = 'EigfilterResponse';
if (nargin >= 1) && ischar(g) && strcmpi(g, 'close')
    delete(findall(0, 'Type','figure', 'Tag', TAG));
    hFig = [];
    return;
end
hFig = findall(0, 'Type','figure', 'Tag', TAG);
if isempty(hFig)
    hFig = figure('Tag', TAG, 'Name', 'Kernel response g(\lambda)', ...
        'NumberTitle','off', 'Color','w', 'MenuBar','none', 'ToolBar','none');
    ax = axes('Parent', hFig); hold(ax,'on');
    xlabel(ax, '\lambda (eigenvalue)'); ylabel(ax, 'g(\lambda)'); grid(ax,'on');
    setappdata(hFig, 'Axes', ax);
else
    hFig = hFig(1);
end
ax = getappdata(hFig, 'Axes');
lambdas = double(lambdas(:));
lmax = max([lambdas; eps]);
xg = linspace(0, lmax, 400)';
yg = g(xg);
cla(ax); hold(ax,'on');
% faint ticks at the real eigenvalues
yl = [min([0; yg(:)]), max([1; yg(:)])];
plot(ax, [lambdas lambdas]', repmat(yl(:), 1, numel(lambdas)), '-', 'Color',[.85 .85 .85]);
plot(ax, xg, yg, 'b-', 'LineWidth', 2);
xlim(ax, [0 lmax]); ylim(ax, yl + [0 max(eps, 0.05*diff(yl))]);
xlabel(ax, '\lambda (eigenvalue)'); ylabel(ax, 'g(\lambda)'); grid(ax,'on');
if nargin >= 3 && ~isempty(titleStr); title(ax, titleStr, 'Interpreter','none'); end
end
