function varargout = view_eigenmodes_timeseries(varargin)
% VIEW_EIGENMODES_TIMESERIES: Plot eigenmode coefficients theta_k(t) over time.
%
% USAGE:  hFig = view_eigenmodes_timeseries(DataFile)
%         [iRows, Labels, Hemi] = view_eigenmodes_timeseries('GetBandTraces', Component, CompRank, kLo, kHi)
%         colors = view_eigenmodes_timeseries('HemiColors', Hemi)
%         view_eigenmodes_timeseries('ModesChangedCallback', hFig)
%
% One trace per eigenmode coefficient (sensor->mode transform), for the paired
% ranks in the EigenModes panel band. Each paired rank yields a left and a right
% trace. The figure tracks the panel band live via bst_figures('FireModesChanged').
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

if (nargin >= 1) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'GetBandTraces','HemiColors','ModesChangedCallback'}))
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: band (paired-rank) -> raw-column traces =====
% For each paired rank k in kLo:kHi, emit its left column(s) then right column(s).
% Labels carry an L/R suffix only when the data has two components.
function [iRows, Labels, Hemi] = GetBandTraces(Component, CompRank, kLo, kHi)
    Component = Component(:);
    CompRank  = CompRank(:);
    isPaired  = any(Component == 2);
    iRows = [];
    Labels = {};
    Hemi = [];
    for k = kLo:kHi
        % Left (component 1) then right (component 2)
        for c = find(CompRank == k & Component == 1)'
            iRows(end+1) = c; %#ok<AGROW>
            Hemi(end+1)  = 1; %#ok<AGROW>
            if isPaired
                Labels{end+1} = sprintf('Mode %d L', k); %#ok<AGROW>
            else
                Labels{end+1} = sprintf('Mode %d', k);   %#ok<AGROW>
                Hemi(end)     = 0;
            end
        end
        for c = find(CompRank == k & Component == 2)'
            iRows(end+1) = c; %#ok<AGROW>
            Hemi(end+1)  = 2; %#ok<AGROW>
            Labels{end+1} = sprintf('Mode %d R', k); %#ok<AGROW>
        end
    end
end


%% ===== PURE: hemisphere id -> RGB cell array (left warm, right cool) =====
function colors = HemiColors(Hemi)
    Hemi = Hemi(:)';
    cL = [0.85 0.33 0.10];   % left  = warm orange
    cR = [0.00 0.45 0.74];   % right = cool blue
    c0 = [0.20 0.20 0.20];   % single-component = neutral
    colors = cell(1, numel(Hemi));
    for i = 1:numel(Hemi)
        switch Hemi(i)
            case 1, colors{i} = cL;
            case 2, colors{i} = cR;
            otherwise, colors{i} = c0;
        end
    end
end
