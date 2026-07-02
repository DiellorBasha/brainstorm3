function out = bst_nxr_registry(action, varargin)
% BST_NXR_REGISTRY: Accessor for the nxr-compute v0.2.0 operator/field registry
% and the Brainstorm Variant -> registry-id mapping (single source of truth).
%
% USAGE:
%   meta = bst_nxr_registry('operator', id)            % operatorInfo struct or []
%   meta = bst_nxr_registry('field',    id)            % fieldInfo struct or []
%   id   = bst_nxr_registry('idForVariant', Variant)   % primary registry id or ''
%   ids  = bst_nxr_registry('componentsForVariant', Variant)  % cellstr or {}
%   fs   = bst_nxr_registry('fieldspec', Variant)      % struct('field_type','domain') or []
%
% The operator/field calls are guarded: a pre-registry nxr binary or an unknown
% id yields [] (never an error), so callers can adopt the registry without a
% hard dependency on it.
%
% 'fieldspec' is a LOCAL declaration (not backed by the nxr binary): it covers
% every Brainstorm operator Variant, including the connectome family (LB-Connectome,
% Connectome Laplacian, Dirac-Connectome) which has no nxr registry id at all.
% This is the single source of truth eigen routing (bst_eigen 'FieldSpec') reads.
%
% Authors: Diellor Basha, 2026

    switch lower(action)
        case 'operator'
            out = local_info('operatorInfo', varargin{1});
        case 'field'
            out = local_info('fieldInfo', varargin{1});
        case 'idforvariant'
            [pid, ~] = local_map(varargin{1});
            out = pid;
        case 'componentsforvariant'
            [~, comps] = local_map(varargin{1});
            out = comps;
        case 'fieldspec'
            out = local_fieldspec(varargin{1});
        otherwise
            error('bst_nxr_registry:badAction', 'Unknown action: %s', action);
    end
end

function meta = local_info(cmd, id)
    meta = [];
    if isempty(id) || ~ischar(id), return; end
    try
        meta = nxr_compute(cmd, id);
    catch
        meta = [];     % pre-registry binary or unknown id
    end
end

function [pid, comps] = local_map(Variant)
% Variant -> primary registry id + component ids. Keep in lockstep with the
% Variants built by tess_operators.
    pid = ''; comps = {};
    switch Variant
        case 'Laplace-Beltrami'
            pid = 'laplaceBeltrami';               comps = {'massGalerkin'};
        case 'Connection Laplacian'
            pid = 'leviCivitaConnectionLaplacian'; comps = {'massGalerkin'};
        case 'Dirac'
            pid = 'relativeDirac';                 comps = {'intrinsicDirac','extrinsicDirac','massGalerkin'};
        case 'Dirac-Face'
            pid = 'relativeFaceDirac';             comps = {'intrinsicFaceDirac','extrinsicFaceDirac','massLumped'};
        case 'Hodge-Face'
            pid = 'faceLaplacianGreenGauss';       comps = {'faceGradient'};
        case 'Covariant'
            pid = 'flatCovariantLaplacian';        comps = {'laplaceBeltrami','massGalerkin'};
    end
end

function fs = local_fieldspec(Variant)
% Local (field_type, domain) for EVERY Brainstorm operator variant — including the
% connectome family that has no nxr binary id. Single source of truth for eigen routing.
    fs = [];
    switch Variant
        case {'Laplace-Beltrami','LB-Connectome','Connectome Laplacian'}
            fs = struct('field_type','real',       'domain','vertex');
        case {'Dirac','Dirac-Connectome'}
            fs = struct('field_type','quaternion', 'domain','vertex');
        case {'Dirac-Face','Hodge-Face'}
            fs = struct('field_type','quaternion', 'domain','face');
        case 'Connection Laplacian'
            fs = struct('field_type','complex',    'domain','vertex');
    end
end
