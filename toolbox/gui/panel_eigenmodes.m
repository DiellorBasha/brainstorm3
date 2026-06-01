function varargout = panel_eigenmodes(varargin)
% PANEL_EIGENMODES: Eigenmode scale lever — a stateful, surface-scoped,
% live-broadcast selection in eigenmode (spatial-frequency) space. Changing the
% selection live-coarsens/smooths a displayed source map via a non-destructive
% display-time band-limited reconstruction.
%
% USAGE:  bstPanel = panel_eigenmodes('CreatePanel')
%         W  = panel_eigenmodes('BuildWeights', shape, kLo, kHi, iCenter, K)
%         panel_eigenmodes('SetBand', kLo, kHi)
%         panel_eigenmodes('SetCurrentMode', k)
%         panel_eigenmodes('SetWindowShape', shape)
%         panel_eigenmodes('SetActive', isActive)
%         W  = panel_eigenmodes('GetWeights')
%         uF = panel_eigenmodes('ApplyToColumn', SurfaceFile, u)
%         panel_eigenmodes('UpdatePanel', hFig)
%
% SEE ALSO: bst_eigenmodes_filter, in_tess_eigenmodes, panel_freq, panel_surface

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

eval(macro_method);
end


%% ===== PURE: build the [1 x K] weight vector from a window shape =====
function W = BuildWeights(shape, kLo, kHi, iCenter, K) %#ok<DEFNU>
    % Clamp all indices into [1,K]
    kLo     = min(max(round(kLo),     1), K);
    kHi     = min(max(round(kHi),     1), K);
    iCenter = min(max(round(iCenter), 1), K);
    if (kHi < kLo)
        [kLo, kHi] = deal(kHi, kLo);
    end
    W = zeros(1, K);
    switch lower(shape)
        case 'single'
            W(iCenter) = 1;
        case 'box'
            W(kLo:kHi) = 1;
        case 'tapered'
            % Tukey window over [kLo,kHi]: cosine shoulders taper to 0 at kLo and kHi;
            % flat interior (=1) over the central portion of the band.
            n = kHi - kLo + 1;
            if (n <= 2)
                W(kLo:kHi) = 1;
            else
                r = 0.5;    % total taper fraction (each shoulder spans r/2 = 25% of the window)
                t = linspace(0, 1, n);
                wt = ones(1, n);
                edge = (r/2);
                iL = t < edge;
                iR = t > (1 - edge);
                wt(iL) = 0.5 * (1 + cos(pi * (2*t(iL)/r - 1)));
                wt(iR) = 0.5 * (1 + cos(pi * (2*t(iR)/r - 2/r + 1)));
                W(kLo:kHi) = wt;
            end
        case 'gain'
            % Gaussian bell centered at iCenter; sigma = half the band width.
            sigma = max((kHi - kLo) / 2, 1);
            k = 1:K;
            W = exp(-0.5 * ((k - iCenter) / sigma).^2);
        otherwise
            error('panel_eigenmodes:BuildWeights: unknown shape "%s".', shape);
    end
end


%% ===== STATE: lazy default =====
function st = GetState()
    global GlobalData;
    if ~isfield(GlobalData, 'UserModes') || isempty(GlobalData.UserModes) ...
            || ~isfield(GlobalData.UserModes, 'Weights')
        GlobalData.UserModes = struct(...
            'SurfaceFile',  '', ...
            'nModes',       0,  ...
            'iCurrentMode', 1,  ...
            'Weights',      [], ...
            'WindowShape',  'box', ...
            'Band',         [1 1], ...
            'BandSpan',     0, ...
            'isActive',     0,  ...
            'CacheSurfaceFile', '', ...
            'CacheEig',     [], ...
            'CacheMass',    []);
    end
    st = GlobalData.UserModes;
end

%% ===== STATE: reset for a surface with K modes =====
function ResetState(SurfaceFile, K) %#ok<DEFNU>
    global GlobalData;
    GetState();                          % ensure struct exists
    GlobalData.UserModes.SurfaceFile  = SurfaceFile;
    GlobalData.UserModes.nModes       = K;
    GlobalData.UserModes.Band         = [1, min(30, K)];
    GlobalData.UserModes.BandSpan     = GlobalData.UserModes.Band(2) - GlobalData.UserModes.Band(1);
    GlobalData.UserModes.iCurrentMode = round(mean(GlobalData.UserModes.Band));
    GlobalData.UserModes.WindowShape  = 'box';
    GlobalData.UserModes.isActive     = 0;
    RecomputeWeights();
end

%% ===== STATE: recompute Weights from shape/band/center =====
function RecomputeWeights()
    global GlobalData;
    st = GlobalData.UserModes;
    GlobalData.UserModes.Weights = BuildWeights(st.WindowShape, ...
        st.Band(1), st.Band(2), st.iCurrentMode, st.nModes);
end

%% ===== STATE: set band (clamped); recentre iCurrentMode to band midpoint =====
function SetBand(kLo, kHi) %#ok<DEFNU>
    global GlobalData;
    st = GetState();
    K  = st.nModes;
    kLo = min(max(round(kLo), 1), K);
    kHi = min(max(round(kHi), 1), K);
    if (kHi < kLo), [kLo, kHi] = deal(kHi, kLo); end
    GlobalData.UserModes.Band         = [kLo, kHi];
    GlobalData.UserModes.BandSpan     = kHi - kLo;
    GlobalData.UserModes.iCurrentMode = round((kLo + kHi) / 2);
    RecomputeWeights();
    NotifyChanged();
end

%% ===== STATE: move center (coupled — slide band, preserve width) =====
function SetCurrentMode(k) %#ok<DEFNU>
    global GlobalData;
    st = GetState();
    K  = st.nModes;
    k  = min(max(round(k), 1), K);
    halfW = round(st.BandSpan / 2);
    kLo = min(max(k - halfW, 1), K);
    kHi = min(max(k + halfW, 1), K);
    GlobalData.UserModes.iCurrentMode = k;
    GlobalData.UserModes.Band         = [kLo, kHi];
    RecomputeWeights();
    NotifyChanged();
end

%% ===== STATE: window shape ('single' collapses band to center) =====
function SetWindowShape(shape) %#ok<DEFNU>
    global GlobalData;
    GetState();
    GlobalData.UserModes.WindowShape = lower(shape);
    if strcmpi(shape, 'single')
        c = GlobalData.UserModes.iCurrentMode;
        GlobalData.UserModes.Band = [c, c];
    end
    RecomputeWeights();
    NotifyChanged();
end

%% ===== STATE: active toggle =====
function SetActive(isActive) %#ok<DEFNU>
    global GlobalData;
    GetState();
    GlobalData.UserModes.isActive = logical(isActive);
    NotifyChanged();
end

%% ===== STATE: read canonical weights =====
% Precondition: ResetState must have been called; returns [] on uninitialized state.
function W = GetWeights() %#ok<DEFNU>
    st = GetState();
    W = st.Weights;
end

%% ===== Broadcast (real implementation added in Task 5) =====
function NotifyChanged()
    % Placeholder until Task 5 wires the figure repaint broadcast.
end


%% ===== CACHE: store eigenmodes + mass matrix for a surface =====
function SetCache(SurfaceFile, Eig, MassMatrix) %#ok<DEFNU>
    global GlobalData;
    GetState();
    GlobalData.UserModes.CacheSurfaceFile = SurfaceFile;
    GlobalData.UserModes.CacheEig         = Eig;
    GlobalData.UserModes.CacheMass        = MassMatrix;
end

%% ===== CACHE: ensure eigenmodes + mass are loaded for a surface =====
function isOk = EnsureCache(SurfaceFile)
    global GlobalData;
    st = GetState();
    isOk = false;
    if ~isempty(st.CacheEig) && ~isempty(st.CacheMass) ...
            && file_compare(st.CacheSurfaceFile, SurfaceFile)
        isOk = true;
        return;
    end
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed || isempty(Eig) || ~isfield(Eig, 'Vectors') || isempty(Eig.Vectors)
        return;
    end
    % Prefer the mass matrix stored with the eigenmodes; recompute only if absent.
    if isfield(Eig, 'MassMatrix') && ~isempty(Eig.MassMatrix)
        M = Eig.MassMatrix;
    else
        sSurf = in_tess_bst(SurfaceFile, 0);
        [~, M] = tess_laplacian(sSurf.Vertices, sSurf.Faces, 'MassType', Eig.MassType);
    end
    SetCache(SurfaceFile, Eig, M);
    isOk = true;
end

%% ===== APPLY: filter a displayed source column (guarded, non-destructive) =====
function uF = ApplyToColumn(SurfaceFile, u) %#ok<DEFNU>
    global GlobalData;
    uF = u;                                   % default: unchanged
    st = GetState();
    % Guards: lever off, surface mismatch, empty column
    if ~st.isActive || isempty(u) || isempty(SurfaceFile) ...
            || ~file_compare(st.SurfaceFile, SurfaceFile)
        return;
    end
    if ~EnsureCache(SurfaceFile)
        return;
    end
    Eig = GlobalData.UserModes.CacheEig;
    M   = GlobalData.UserModes.CacheMass;
    % Scalar-field guard: only filter when the column matches the mesh vertex
    % count (skip unconstrained/volume/mismatched maps rather than mis-filter).
    if (size(u,1) ~= size(Eig.Vectors,1)) || (size(u,2) ~= 1)
        return;
    end
    W = GlobalData.UserModes.Weights;
    if isempty(W) || (numel(W) ~= Eig.nModes)
        return;
    end
    % Reconstruct via the core spectral filter (custom transfer = our weights)
    uF = bst_eigenmodes_filter(Eig, u, M, 'custom', 'TransferFn', @(l) W(:));
end
