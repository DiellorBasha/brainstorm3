function Gauge = tess_gauge(SurfaceFile, varargin)
% TESS_GAUGE: Store the nxr-compute Gauge package per hemisphere.
%
% USAGE:  Gauge = tess_gauge(SurfaceFile)
%         Gauge = tess_gauge(SurfaceFile, 'Gauge','trivial', 'Operators',1, 'Coupling','ambient', 'Mass','lumped')
%
% Stores TessMat.Gauge as a 1x2 per-hemisphere struct array. The trivial gauge
% places singularities at the FreeSurfer sphere poles. With 'Operators',1 the
% nxr connection Laplacian + covariant Laplacian are attached.
%
% OPTIONS:
%     'Gauge'           'trivial' (default) | 'levi-civita' | 'euclidean'
%                       (only 'trivial' places sphere-pole singularities)
%     'Operators'       false (default) | true   - attach .operators (laplacian, covariantLaplacian)
%     'Coupling'        'ambient' (default) | other nxr coupling strings
%     'Mass'            'lumped' (default) | 'galerkin'
%     'NoSave'          false (default) | true
%     'ForceRecompute'  false (default) | true
%
% OUTPUT (also stored as TessMat.Gauge, 1x2):
%     Gauge   - 1x2 struct array, (1)=left, (2)=right
%
% Requires the nxr-compute plugin. The trivial gauge requires a FreeSurfer
% registration sphere (TessMat.Reg.Sphere.Vertices).
%
% SEE ALSO: tess_bundle, tess_store_perhemi, tess_topology, tess_geometry

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

    GaugeType='trivial'; Operators=false; Coupling='ambient'; Mass='lumped';
    NoSave=false; ForceRecompute=false;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'gauge',          GaugeType=varargin{i+1};
            case 'operators',      Operators=varargin{i+1};
            case 'coupling',       Coupling=varargin{i+1};
            case 'mass',           Mass=varargin{i+1};
            case 'nosave',         NoSave=varargin{i+1};
            case 'forcerecompute', ForceRecompute=varargin{i+1};
        end
    end
    GaugeType = lower(GaugeType);

    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end   %#ok<CTCH>
    prov = struct('Backend','nxr', 'Package','gauge', 'NxrVersion',nxrVer, 'Gauge',GaugeType, ...
                  'Operators',logical(Operators), 'Coupling',Coupling, 'Mass',Mass, ...
                  'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    computeFn = @(h, poles) local_compute_gauge(h, poles, GaugeType, Operators, Coupling, Mass);
    B = tess_store_perhemi(SurfaceFile, {'Gauge'}, strcmpi(GaugeType,'trivial'), prov, ...
            ~ForceRecompute && ~Operators, NoSave, computeFn);
    Gauge = B.Gauge;
end

function S = local_compute_gauge(h, poles, GaugeType, Operators, Coupling, Mass)
    opts = struct();
    if strcmpi(GaugeType,'trivial')
        opts.singVerts = poles; opts.singValues = [1; 1];
    end
    if Operators
        opts.operators = true; opts.coupling = Coupling; opts.mass = Mass;
    end
    S = struct('Gauge', nxr_compute('gauge', h, GaugeType, opts));
end
