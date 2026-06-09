function Topology = tess_topology(SurfaceFile, varargin)
% TESS_TOPOLOGY: Store the nxr-compute halfedge Topology package per hemisphere.
%
% USAGE:  Topology = tess_topology(SurfaceFile)
%         Topology = tess_topology(SurfaceFile, 'Operators',1, 'NoSave',1, 'ForceRecompute',1)
%
% Stores TessMat.Topology as a 1x2 per-hemisphere struct array ((1)=left,(2)=right).
% With 'Operators',1 the nxr graph laplacian + DEC operators are attached.
%
% OPTIONS:
%     'Operators'       false (default) | true   - attach .operators (laplacian, dec)
%     'NoSave'          false (default) | true
%     'ForceRecompute'  false (default) | true
%
% OUTPUT (also stored as TessMat.Topology, 1x2):
%     Topology   - 1x2 struct array, (1)=left, (2)=right
%
% Requires the nxr-compute plugin.
%
% SEE ALSO: tess_bundle, tess_store_perhemi, tess_hemisplit

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

    % --- options ---
    Operators=false; NoSave=false; ForceRecompute=false;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'operators',      Operators=varargin{i+1};
            case 'nosave',         NoSave=varargin{i+1};
            case 'forcerecompute', ForceRecompute=varargin{i+1};
        end
    end

    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end   %#ok<CTCH>
    prov = struct('Backend','nxr', 'Package','topology', 'NxrVersion',nxrVer, ...
                  'Operators',logical(Operators), ...
                  'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    computeFn = @(h, poles) local_compute_topology(h, Operators);
    B = tess_store_perhemi(SurfaceFile, {'Topology'}, false, prov, ...
            ~ForceRecompute && ~Operators, NoSave, computeFn);
    Topology = B.Topology;
end

function S = local_compute_topology(h, Operators)
    opts = struct();
    if Operators, opts.operators = true; end
    S = struct('Topology', nxr_compute('topology', h, opts));
end
