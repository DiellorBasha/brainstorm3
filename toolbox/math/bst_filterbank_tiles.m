function Tiles = bst_filterbank_tiles(base)
% BST_FILTERBANK_TILES: Expand one base filter design into a spectrum-spanning bank.
%
% USAGE:  Tiles = bst_filterbank_tiles(base)
%
% INPUT base (struct):
%   .Kernel      kernel name in the bst_eigfilter registry (e.g. 'mexhat','heat')
%   .Params      struct of fixed kernel params (the scale param is overwritten per tile)
%   .Direction   [1x3] launch direction (Dirac); copied to every tile
%   .Chirality   scalar 0 | +1 | -1 (the base sign; ignored if .Chiralities set)
%   .Axis        [1x3] chirality axis (Dirac); copied to every tile
%   .N           number of spectral tiles (>=1)
%   .Spacing     'geometric' (default) | 'linear' — spacing of tile centers in lambda
%   .LambdaRange [lo hi] eigenvalue band the tiles tile (lo>0 for geometric)
%   .Chiralities [] or 0 => single bank at .Chirality; a vector (e.g. [1 -1]) crosses
%                every spectral tile with each listed sign (doubling/tripling the bank)
%
% OUTPUT:
%   Tiles : 1xM struct array, each (Kernel,Params,Direction,Chirality,Axis), where
%           M = N * max(1,numel(.Chiralities)). The scale parameter of .Params is set
%           per tile so the tile CENTERS (the lambda of peak response) span LambdaRange.
%
% The scale param name and the center<->param map are read from the kernel family
% (bandpass kernels like mexhat peak at lambda=1/t; heat uses t=1/center as a soft
% cutoff; tikhonov uses beta as the cutoff). Unknown kernels fall back to a linear
% sweep of the kernel's first metadata parameter.
%
% SEE ALSO: bst_eigfilter_kernel, bst_dirac_eigenmodes_filter, panel_filter_designer
%
% Authors: Diellor Basha, 2026

    if ~isfield(base,'N') || isempty(base.N) || base.N < 1, base.N = 1; end
    if ~isfield(base,'Spacing') || isempty(base.Spacing), base.Spacing = 'geometric'; end
    if ~isfield(base,'LambdaRange') || numel(base.LambdaRange) ~= 2, base.LambdaRange = [1 100]; end
    if ~isfield(base,'Params') || ~isstruct(base.Params), base.Params = struct(); end
    if ~isfield(base,'Direction') || isempty(base.Direction), base.Direction = [1 0 0]; end
    if ~isfield(base,'Axis') || isempty(base.Axis), base.Axis = [0 0 1]; end
    if ~isfield(base,'Chirality') || isempty(base.Chirality), base.Chirality = 0; end
    if ~isfield(base,'Chiralities'), base.Chiralities = []; end

    lo = base.LambdaRange(1); hi = base.LambdaRange(2);
    N  = base.N;
    if strcmpi(base.Spacing,'geometric')
        if lo <= 0, lo = max(eps, hi*1e-3); end
        if N == 1, centers = sqrt(lo*hi); else, centers = exp(linspace(log(lo), log(hi), N)); end
    else
        if N == 1, centers = (lo+hi)/2; else, centers = linspace(lo, hi, N); end
    end

    [scaleName, centerToParam] = i_scale_map(base.Kernel);

    signs = base.Chiralities;
    if isempty(signs), signs = base.Chirality; end

    Tiles = repmat(i_blank_tile(), 1, N*numel(signs));
    t = 0;
    for s = 1:numel(signs)
        for j = 1:N
            t = t + 1;
            p = base.Params;
            % N==1 is the user's single designed wavelet: keep its scale param as-is.
            % N>1 tiles the spectrum: each tile's scale is set from its spectral center.
            if N > 1
                p.(scaleName) = centerToParam(centers(j));
            end
            Tiles(t) = struct('Kernel', base.Kernel, 'Params', p, ...
                'Direction', base.Direction(:).', 'Chirality', signs(s), 'Axis', base.Axis(:).');
        end
    end
end

function t = i_blank_tile()
    t = struct('Kernel','', 'Params',struct(), 'Direction',[1 0 0], 'Chirality',0, 'Axis',[0 0 1]);
end

function [scaleName, centerToParam] = i_scale_map(kernelName)
% Map a desired peak/cutoff lambda (center) to the kernel's scale parameter.
    switch lower(kernelName)
        case {'mexhat'}                       % peak at lambda = 1/t  => t = 1/center
            scaleName = 't';    centerToParam = @(c) 1./max(c, eps);
        case {'heat','inverse_heat'}          % e^{-t*lambda}: use t = 1/center as a soft cutoff
            scaleName = 't';    centerToParam = @(c) 1./max(c, eps);
        case {'tikhonov'}                     % beta acts at lambda ~ beta
            scaleName = 'beta'; centerToParam = @(c) c;
        otherwise                             % fallback: first meta param swept linearly as the center
            meta = bst_eigfilter_kernel('info', kernelName);
            fn = fieldnames(meta.params);
            scaleName = fn{1};  centerToParam = @(c) c;
    end
end
