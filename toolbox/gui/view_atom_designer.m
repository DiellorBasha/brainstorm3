function hFig = view_atom_designer(SurfaceFile, variant, nModes, seed0)
% VIEW_ATOM_DESIGNER: Interactive joint time-cortex WAVELET ATOM designer. Opens a figure_3d cortex view
% with an in-window control panel (kernel + parameters); CLICK a vertex to drop a delta and propagate it
% as a dynamic eigenwavelet (atom) through the cortex eigenbasis (geometry, or geometry+connectome with
% the LB-Connectome operator). Step through the propagation in time with the LEFT/RIGHT arrow keys (the
% same time-stepping idiom as a source map). Modelled on view_manifold (self-contained view + embedded
% controls).
%
% USAGE:  hFig = view_atom_designer(SurfaceFile)                       % default LB-Connectome
%         hFig = view_atom_designer(SurfaceFile, 'Laplace-Beltrami')  % geometry-only
%         hFig = view_atom_designer(SurfaceFile, variant, nModes)
%
% SEE ALSO: bst_atom, bst_eigfilter_kernel, tess_eigen, view_manifold

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% Copyright (c) University of Southern California & McGill University
% Distributed under the terms of the GNU General Public License (GPLv3).
% =============================================================================@
%
% Authors: Diellor Basha, 2026

    if (nargin < 2) || isempty(variant), variant = 'LB-Connectome'; end
    if (nargin < 3) || isempty(nModes),  nModes  = 200; end

    % --- cortex figure (figure_3d) ---
    hFig = view_surface(SurfaceFile, 0, [.6 .6 .6], 'NewFigure', 0);
    if isempty(hFig), bst_error('Could not open the surface.', 'Atom designer', 0); return; end
    set(hFig, 'Name', ['Atom designer: ' SurfaceFile]);
    hAxes  = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    hPatch = findobj(hAxes, 'Type', 'patch'); if numel(hPatch) > 1, hPatch = hPatch(1); end

    % --- eigenbasis (cortex axis) + synthetic propagation time ---
    bst_progress('start', 'Atom designer', sprintf('Building %s eigenbasis...', variant));
    E  = tess_eigen(SurfaceFile, variant, 'nModes', nModes, 'NoSave', true);
    Op = in_bst_operator(E.OperatorFile);
    bst_progress('stop');
    ax = struct('Phi',{E.Phi}, 'Lambda',{E.Lambda}, 'Mass',{Op.Mass}, 'GlobalVertices',{E.GlobalVertices});
    V  = get(hPatch, 'Vertices');  nV = size(V, 1);

    % --- state (closure-shared) ---
    Textent = 2.0;  nFrames = 60;  ax.tlag = linspace(0, Textent, nFrames);  ax.nT = nFrames;
    kernel  = 'dampedwave';  kparams = struct('alpha', 3, 'beta', 0.5);
    W = [];  curFrame = 1;  seedVtx = [];
    cmap = [ [linspace(0,1,32)';ones(32,1)], [linspace(0,1,32)';linspace(1,0,32)'], [ones(32,1);linspace(1,0,32)'] ];

    % --- embedded control panel ---
    uicontrol(hFig, 'Style','text', 'String','Kernel:', 'Units','Pixels', 'Position',[8 40 50 18], 'HorizontalAlignment','left');
    hKern = uicontrol(hFig, 'Style','popupmenu', 'String',{'dampedwave','wave','kleingordon','heat'}, ...
        'Units','Pixels', 'Position',[60 40 110 22], 'Callback',@KernelChanged);
    uicontrol(hFig, 'Style','text', 'String','alpha', 'Units','Pixels', 'Position',[180 40 36 18], 'HorizontalAlignment','left');
    hA = uicontrol(hFig, 'Style','edit', 'String','3',   'Units','Pixels', 'Position',[214 40 42 22], 'Callback',@ParamChanged);
    uicontrol(hFig, 'Style','text', 'String','beta',  'Units','Pixels', 'Position',[262 40 32 18], 'HorizontalAlignment','left');
    hB = uicontrol(hFig, 'Style','edit', 'String','0.5', 'Units','Pixels', 'Position',[296 40 42 22], 'Callback',@ParamChanged);
    hLabel = uicontrol(hFig, 'Style','text', 'String','Click a vertex to drop an atom; LEFT/RIGHT arrows step time.', ...
        'Units','Pixels', 'Position',[8 8 700 18], 'HorizontalAlignment','left', 'FontUnits','points', 'FontSize',bst_get('FigFont'));

    % --- interactions ---
    set(hFig, 'WindowButtonDownFcn', @ClickPick);
    set(hFig, 'KeyPressFcn', @KeyStep);
    colormap(hAxes, cmap);
    if (nargin >= 4) && ~isempty(seed0), seedVtx = seed0; Generate(); end   % scripted/test seed

    % ===== nested callbacks =====
    function KernelChanged(src, ~)
        opt = get(src, 'String');  kernel = opt{get(src, 'Value')};
        en = 'on'; if strcmpi(kernel,'heat'), en = 'off'; end
        set([hA hB], 'Enable', en);
        if ~isempty(seedVtx), Generate(); end
    end
    function ParamChanged(~, ~)
        kparams.alpha = str2double(get(hA,'String'));  kparams.beta = str2double(get(hB,'String'));
        if ~isempty(seedVtx), Generate(); end
    end
    function ClickPick(~, ~)
        if ~strcmpi(get(hFig,'SelectionType'),'normal'), return; end       % left-click only
        cp = get(hAxes, 'CurrentPoint');  o = cp(1,:);  d = cp(2,:) - cp(1,:);  d = d / max(norm(d),eps);
        w = V - o;  t = w * d';  proj = o + t .* d;  dist = sqrt(sum((V - proj).^2, 2));
        cand = find(dist < 0.1*max(dist) | dist <= min(dist)*3);  if isempty(cand), [~,seedVtx]=min(dist); else, [~,j]=min(t(cand)); seedVtx=cand(j); end
        Generate();
    end
    function Generate()
        switch lower(kernel)
            case 'heat',        kp = struct('t', 0.02);
            case 'kleingordon', kp = struct('alpha', kparams.alpha, 'mu', 0.1);
            otherwise,          kp = struct('alpha', kparams.alpha, 'beta', kparams.beta);
        end
        G = bst_dynamics('NewGroup','atom');  G.vertices = seedVtx;  G.pos = V(seedVtx,:);
        try
            [Wloc, gv] = bst_atom('Evaluate', G, 1, ax, kernel, kp);
            W = zeros(nV, ax.nT);  W(gv,:) = Wloc;  curFrame = 1;  Show();
        catch ME
            set(hLabel, 'String', ['Could not propagate from this vertex: ' regexprep(ME.message,'\s+',' ')]);
        end
    end
    function Show()
        set(hPatch, 'FaceVertexCData', W(:,curFrame), 'FaceColor','interp');
        m = max(abs(W(:,curFrame))) + eps;  clim(hAxes, [-m m]);
        set(hLabel, 'String', sprintf('Atom @ vertex %d | %s (alpha=%.2g beta=%.2g) | frame %d/%d  t=%.2f s   [LEFT/RIGHT = step time]', ...
            seedVtx, kernel, kparams.alpha, kparams.beta, curFrame, ax.nT, ax.tlag(curFrame)));
    end
    function KeyStep(~, ev)
        if isempty(W), return; end
        switch ev.Key
            case 'rightarrow', curFrame = min(curFrame+1, ax.nT);  Show();
            case 'leftarrow',  curFrame = max(curFrame-1, 1);        Show();
            case 'home',       curFrame = 1;                          Show();
        end
    end
end
