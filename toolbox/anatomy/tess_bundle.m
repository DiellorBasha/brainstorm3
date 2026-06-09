function B = tess_bundle(SurfaceFile, varargin)
% TESS_BUNDLE: Store the nxr-compute coordinate bundle as per-hemisphere fields.
%
% USAGE:  B = tess_bundle(SurfaceFile)
%         B = tess_bundle(SurfaceFile, 'Gauge','trivial', 'NoSave',1, 'ForceRecompute',1)
%
% DESCRIPTION:
%     nxr-compute never integrates the two hemispheres, so each cortex is bundled
%     one connected genus-0 hemisphere at a time. The three results are stored as
%     1x2 struct arrays ((1)=left, (2)=right). Each element is the verbatim nxr
%     submesh bundle (LOCAL indexing) plus GlobalVertices/GlobalFaces/Hemisphere
%     scatter maps to TessMat.Vertices/.Faces global order.
%
% OPTIONS:
%     'Gauge'           'trivial' (default) | 'levi-civita' | 'euclidean'
%     'Operators'       false (default) | true   - heavy .operators (Task 2)
%     'Coupling'        'ambient' (default) | 'product'   - covariantLaplacian (heavy)
%     'Mass'            'lumped' (default) | 'galerkin'    - mass variant (heavy)
%     'NoSave'          false (default) | true
%     'ForceRecompute'  false (default) | true
%
% OUTPUT (also stored as TessMat.Topology/.Geometry/.Gauge, each 1x2):
%     B.Topology, B.Geometry, B.Gauge   - 1x2 struct arrays
%
% Requires the nxr-compute plugin; trivial gauge needs a FreeSurfer reg sphere.
%
% SEE ALSO: tess_tangents, tess_hemisplit, tess_normals

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
    Gauge='trivial'; Operators=false; Coupling='ambient'; Mass='lumped';
    NoSave=false; ForceRecompute=false;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'gauge',          Gauge=varargin{i+1};
            case 'operators',      Operators=varargin{i+1};
            case 'coupling',       Coupling=varargin{i+1};
            case 'mass',           Mass=varargin{i+1};
            case 'nosave',         NoSave=varargin{i+1};
            case 'forcerecompute', ForceRecompute=varargin{i+1};
        end
    end

    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end   %#ok<CTCH>
    prov = struct('Backend','nxr', 'Package','bundle', 'NxrVersion',nxrVer, 'Gauge',Gauge, ...
                  'Operators',logical(Operators), 'Coupling',Coupling, 'Mass',Mass, ...
                  'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    computeFn = @(h, poles) local_compute_bundle(h, poles, Gauge, Operators, Coupling, Mass);
    B = tess_store_perhemi(SurfaceFile, {'Topology','Geometry','Gauge'}, ...
            strcmpi(Gauge,'trivial'), prov, ~ForceRecompute && ~Operators, NoSave, computeFn);
end

function S = local_compute_bundle(h, poles, Gauge, Operators, Coupling, Mass)
    opts = struct();
    if strcmpi(Gauge,'trivial')
        opts.singVerts = poles; opts.singValues = [1; 1];
    end
    if Operators
        opts.operators = true; opts.coupling = Coupling; opts.mass = Mass;
    end
    S = nxr_compute('bundle', h, Gauge, opts);   % returns {Topology,Geometry,Gauge}
end
