function Tiles = bst_filterbank_tiles(wavelet, opts)
% BST_FILTERBANK_TILES: Replicate a single designed wavelet into a spectrum-spanning bank.
%
% USAGE:  Tiles = bst_filterbank_tiles(wavelet, opts)
%
% This is the SPECTRUM TILING module. It is deliberately separate from single-wavelet
% design: it takes a FINISHED wavelet (the canonical design/exploration object) and an
% options struct describing how to replicate it across the eigenvalue spectrum, and
% returns a bank of wavelet recipes. It never reads any GUI state and never decides the
% single wavelet itself — callers that only want one wavelet should use the wavelet
% directly and not call this function.
%
% INPUT wavelet (struct) - one designed wavelet:
%   .Kernel      kernel name in the bst_eigfilter registry (e.g. 'mexhat','heat')
%   .Params      struct of kernel params (the scale param is set per tile when N>1)
%   .Direction   [1x3] ambient launch direction (Dirac); copied to every tile
%   .Chirality   scalar 0 | +1 | -1 (used when opts.Chiralities is empty)
%   .Axis        [1x3] chirality axis (Dirac); copied to every tile
%
% INPUT opts (struct) - how to tile the spectrum:
%   .N           number of spectral tiles (>=1)
%   .Spacing     'geometric' (default) | 'linear' — spacing of tile centers in lambda
%   .LambdaRange [lo hi] eigenvalue band the tiles span (lo>0 for geometric)
%   .Chiralities [] => single chirality (= wavelet.Chirality); a vector (e.g. [1 -1])
%                crosses every spectral tile with each listed sign
%
% OUTPUT:
%   Tiles : 1xM struct array, each (Kernel,Params,Direction,Chirality,Axis), where
%           M = N * max(1,numel(opts.Chiralities)). When N>1 the scale parameter of each
%           tile is set so the tile CENTERS (lambda of peak/cutoff) span LambdaRange.
%           When N==1 the wavelet's own scale param is kept unchanged.
%
% SEE ALSO: bst_eigfilter_kernel, bst_eigfilter_evaluate
%
% Authors: Diellor Basha, 2026

    % --- normalize the single wavelet ---
    if ~isfield(wavelet,'Params') || ~isstruct(wavelet.Params), wavelet.Params = struct(); end
    if ~isfield(wavelet,'Direction') || isempty(wavelet.Direction), wavelet.Direction = [1 0 0]; end
    if ~isfield(wavelet,'Axis') || isempty(wavelet.Axis), wavelet.Axis = [0 0 1]; end
    if ~isfield(wavelet,'Chirality') || isempty(wavelet.Chirality), wavelet.Chirality = 0; end

    % --- normalize the tiling options ---
    if (nargin < 2) || isempty(opts), opts = struct(); end
    if ~isfield(opts,'N') || isempty(opts.N) || opts.N < 1, opts.N = 1; end
    if ~isfield(opts,'Spacing') || isempty(opts.Spacing), opts.Spacing = 'geometric'; end
    if ~isfield(opts,'LambdaRange') || numel(opts.LambdaRange) ~= 2, opts.LambdaRange = [1 100]; end
    if ~isfield(opts,'Chiralities'), opts.Chiralities = []; end

    N  = opts.N;
    lo = opts.LambdaRange(1); hi = opts.LambdaRange(2);
    if strcmpi(opts.Spacing,'geometric')
        if lo <= 0, lo = max(eps, hi*1e-3); end
        if N == 1, centers = sqrt(lo*hi); else, centers = exp(linspace(log(lo), log(hi), N)); end
    else
        if N == 1, centers = (lo+hi)/2; else, centers = linspace(lo, hi, N); end
    end

    [scaleName, centerToParam] = i_scale_map(wavelet.Kernel);

    signs = opts.Chiralities;
    if isempty(signs), signs = wavelet.Chirality; end

    Tiles = repmat(i_blank_tile(), 1, N*numel(signs));
    t = 0;
    for s = 1:numel(signs)
        for j = 1:N
            t = t + 1;
            p = wavelet.Params;
            % N==1: keep the wavelet's designed scale. N>1: set each tile's scale from
            % its spectral center so the bank spans LambdaRange.
            if N > 1
                p.(scaleName) = centerToParam(centers(j));
            end
            Tiles(t) = struct('Kernel', wavelet.Kernel, 'Params', p, ...
                'Direction', wavelet.Direction(:).', 'Chirality', signs(s), 'Axis', wavelet.Axis(:).');
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
