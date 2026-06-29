function hFig = view_atom_designer(SurfaceFile, variant, nModes, seed0)
% VIEW_ATOM_DESIGNER: Interactive joint time-cortex WAVELET ATOM designer. Opens a figure_3d cortex view
% with a top-right filter-design panel. LEFT-CLICK a vertex (on the cortex) to drop a delta and propagate
% it as a dynamic eigenwavelet (atom) through the cortex eigenbasis (geometry with Laplace-Beltrami, or
% geometry+connectome with LB-Connectome). DRAG / click off the cortex rotates the camera, and the
% figure_3d right-click context menu (colormaps, views, ...) is preserved. LEFT/RIGHT arrows step the
% propagation in time; HOME resets. Modelled on view_manifold (self-contained view + embedded controls).
%
% Filter-design panel: pick a FAMILY (Spatial eigenwavelet | Dynamic atom) and a kernel; for dynamic
% kernels a Time/Joint-spectral toggle and the kernel's PHYSICAL parameters appear - spatial scale (mm,
% calibrated to the cortex via 2*pi/sqrt(lambda)), wave speed c (m/s, only for non-separable kernels),
% decay time (s), and the temporal extent (s). State machine: 'design' (place/edit atoms) and 'save'.
%
% USAGE:  hFig = view_atom_designer(SurfaceFile)                       % default LB-Connectome
%         hFig = view_atom_designer(SurfaceFile, 'Laplace-Beltrami')  % geometry-only
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

    % --- cortex figure (figure_3d: keep its camera + context menu) ---
    hFig = view_surface(SurfaceFile, 0, [.6 .6 .6], 'NewFigure', 0);
    if isempty(hFig), bst_error('Could not open the surface.', 'Atom designer', 0); return; end
    set(hFig, 'Name', ['Atom designer: ' SurfaceFile]);
    hAxes  = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    hPatch = findobj(hAxes, 'Type', 'patch'); if numel(hPatch) > 1, hPatch = hPatch(1); end

    % --- eigenbasis (cached) + physical-scale calibration ---
    bst_progress('start', 'Atom designer', sprintf('Building/loading %s eigenbasis...', variant));
    E  = tess_eigen(SurfaceFile, variant, 'nModes', nModes);   % #4: cache + reuse the eigen_ file
    Op = in_bst_operator(E.OperatorFile);
    bst_progress('stop');
    ax = struct('Phi',{E.Phi}, 'Lambda',{E.Lambda}, 'Mass',{Op.Mass}, 'GlobalVertices',{E.GlobalVertices});
    V  = get(hPatch, 'Vertices');  nV = size(V, 1);
    lamAll = E.Lambda{1}(:); if numel(E.Lambda) > 1 && ~isempty(E.Lambda{2}), lamAll = [lamAll; E.Lambda{2}(:)]; end
    lmax = max(lamAll);  lminPos = min(lamAll(lamAll > 1e-9));
    mm   = @(l) 2*pi./sqrt(max(l,eps))*1000;                   % eigenvalue -> physical scale (mm)
    scaleMinMM = mm(lmax);  scaleMaxMM = mm(lminPos);          % finest .. coarsest cortical scale
    meanEdge = i_mean_edge(V, get(hPatch,'Faces'));

    % --- registry: one FLAT kernel list, dynamic first then a divider then static (#2) ---
    allK = bst_eigfilter_kernel('list');  spatialK = {}; dynK = {};
    for ii = 1:numel(allK)
        try m = bst_eigfilter_kernel('info', allK{ii}); catch, continue; end
        if isfield(m,'domain') && ~isempty(m.domain), dynK{end+1}=allK{ii}; else, spatialK{end+1}=allK{ii}; end %#ok<AGROW>
    end
    pref = {'dampedwave','wave','kleingordon','diffusion'};                 % favour the dynamic atoms
    dynK = [pref(ismember(pref,dynK)), setdiff(dynK, pref, 'stable')];
    DH = '──── dynamic ────';  SH = '──── static ────';                    % #1: both group headers
    kList = [{DH}, dynK, {SH}, spatialK];  iDH = 1;  iSH = numel(dynK) + 2;

    % --- state (closure-shared) ---
    state    = 'design';                 % 'design' | 'save'
    kernel   = 'dampedwave';  lastKIdx = 2;
    pScaleMM = round((scaleMinMM+scaleMaxMM)/2);   % spatial scale (mm)
    pSpeed   = 1.0;                      % wave speed c (m/s), non-separable only
    pDecay   = 0.5;                      % decay time (s), damped only
    nFrames  = 100;  ax.nT = nFrames;  ax.tlag = (0:nFrames-1)/100;   % standard 1 s @ 100 Hz time vector
    cmap = [ [linspace(0,1,32)';ones(32,1)], [linspace(0,1,32)';linspace(1,0,32)'], [ones(32,1);linspace(1,0,32)'] ];
    W = [];  curFrame = 1;  seedVtx = [];
    colormap(hAxes, cmap);

    % --- top-right filter-design panel: kernel + sliders, compactly spaced (#2,#3) ---
    hP = uipanel('Parent',hFig, 'Title','Filter design', 'Units','normalized', 'Position',[0.70 0.68 0.295 0.30], ...
                 'FontUnits','points','FontSize',bst_get('FigFont'));
    yk = 0.80; ys = [0.60 0.43 0.26];  lw=0.30; sw=0.42; vw=0.16;
    uicontrol(hP,'Style','text','String','Kernel','Units','normalized','Position',[0.04 yk lw 0.13],'HorizontalAlignment','left');
    hKern = uicontrol(hP,'Style','popupmenu','String',kList,'Value',2,'Units','normalized','Position',[0.35 yk 0.61 0.14],'Callback',@KernelChanged);
    [hScale,hScaleV] = i_slider(hP,'Scale', ys(1), lw,sw,vw, scaleMinMM, scaleMaxMM, pScaleMM, @ParamChanged);
    [hSpeed,hSpeedV] = i_slider(hP,'Speed', ys(2), lw,sw,vw, 0.1, 10,  pSpeed, @ParamChanged);
    [hDecay,hDecayV] = i_slider(hP,'Decay', ys(3), lw,sw,vw, 0.05, 2,  pDecay, @ParamChanged);
    hSave = uicontrol(hP,'Style','pushbutton','String','Save','Units','normalized','Position',[0.80 0.04 0.16 0.14], ...
                      'Callback',@SaveAtom, 'TooltipString','Save atom -> Scout + Event');

    hLabel = uicontrol(hFig,'Style','text','String','[design] Left-click the cortex to drop an atom; drag / click off = rotate; arrows = step time.', ...
        'Units','Pixels','Position',[8 6 800 18],'HorizontalAlignment','left','FontUnits','points','FontSize',bst_get('FigFont'));
    SyncControls();

    % --- chain figure_3d handlers (#1 camera, #8 context menu) ---
    origDown = get(hFig,'WindowButtonDownFcn');  origUp = get(hFig,'WindowButtonUpFcn');  origKey = get(hFig,'KeyPressFcn');
    downXY = [];
    set(hFig, 'WindowButtonDownFcn', @OnDown, 'WindowButtonUpFcn', @OnUp, 'KeyPressFcn', @OnKey);
    if (nargin >= 4) && ~isempty(seed0), seedVtx = seed0; Generate(); end

    % ===== nested callbacks =====
    function KernelChanged(src,~)
        idx = get(src,'Value');
        if idx==iDH || idx==iSH, set(src,'Value',lastKIdx); return; end  % skip the group headers
        lastKIdx = idx;  opt = get(src,'String');  kernel = opt{idx};
        SyncControls();  Regen();
    end
    function ParamChanged(~,~)
        pScaleMM=get(hScale,'Value'); pSpeed=get(hSpeed,'Value'); pDecay=get(hDecay,'Value');
        set(hScaleV,'String',num2str(round(pScaleMM))); set(hSpeedV,'String',sprintf('%.2g',pSpeed)); set(hDecayV,'String',sprintf('%.2g',pDecay));
        Regen();
    end
    function SyncControls()
        m = bst_eigfilter_kernel('info', kernel);
        isDyn = isfield(m,'domain') && ~isempty(m.domain);
        isSep = ~isfield(m,'separable') || m.separable;
        set([hSpeed hSpeedV], 'Enable', i_en(isDyn && ~isSep));          % speed c: non-separable only
        set([hDecay hDecayV], 'Enable', i_en(strcmpi(kernel,'dampedwave'))); % decay: damped only
        set(hLabel,'String',sprintf('[%s] %s%s | scale %.0f mm (cortex %.0f-%.0f mm) | t 0-%.2fs @100Hz', ...
            state, kernel, i_en2b(isDyn&&~isSep,pSpeed), pScaleMM, scaleMinMM, scaleMaxMM, ax.tlag(end)));
    end
    function Regen(), if ~isempty(seedVtx), Generate(); end, end
    function OnDown(h,ev), downXY = get(hFig,'CurrentPoint'); i_call(origDown,h,ev); end
    function OnUp(h,ev)
        i_call(origUp,h,ev);
        if strcmpi(state,'design') && strcmpi(get(hFig,'SelectionType'),'normal') && ~isempty(downXY)
            if norm(get(hFig,'CurrentPoint')-downXY) < 4               % click, not a drag
                [v,hit] = PickVertex();  if hit, seedVtx=v; Generate(); end
            end
        end
    end
    function OnKey(h,ev)
        if ~isempty(W)
            switch ev.Key
                case 'rightarrow', curFrame=min(curFrame+1,ax.nT); Show(); return;
                case 'leftarrow',  curFrame=max(curFrame-1,1);     Show(); return;
                case 'home',       curFrame=1;                     Show(); return;
            end
        end
        i_call(origKey,h,ev);
    end
    function [v,hit] = PickVertex()
        cp=get(hAxes,'CurrentPoint'); o=cp(1,:); d=cp(2,:)-cp(1,:); d=d/max(norm(d),eps);
        w=V-o; t=w*d'; proj=o+t.*d; dist=sqrt(sum((V-proj).^2,2));
        near=find(dist < 2*meanEdge);                                  % the click ray must graze the cortex
        if isempty(near), v=[]; hit=false; return; end
        [~,j]=min(t(near)); v=near(j); hit=true;                        % front-most grazed vertex
    end
    function Generate()
        kp = i_phys2kernel();
        G = bst_dynamics('NewGroup','atom'); G.vertices=seedVtx; G.pos=V(seedVtx,:);
        try
            [Wloc,gv]=bst_atom('Evaluate', G, 1, ax, kernel, kp);
            W=zeros(nV,ax.nT); W(gv,:)=Wloc; curFrame=1; Show();
        catch ME
            set(hLabel,'String',['Could not propagate: ' regexprep(ME.message,'\s+',' ')]);
        end
    end
    function kp = i_phys2kernel()
        % physical params -> abstract kernel params, using the cortex calibration
        kp = struct('lmax', lmax);
        switch lower(kernel)
            case {'wave','dampedwave','kleingordon'}
                kp.alpha = pSpeed * sqrt(lmax) / 2;                    % c (m/s) -> alpha
                if strcmpi(kernel,'dampedwave'), kp.beta = 1/max(pDecay,eps); end
                if strcmpi(kernel,'kleingordon'), kp.mu = 0.1*lmax; end
            case 'diffusion'
                kp.tau = max((pScaleMM/1000)^2 * lmax, eps);          % scale (mm) -> diffusivity proxy
            case 'heat'
                lamS = (2*pi/(pScaleMM/1000))^2; kp.t = log(2)/max(lamS,eps);   % scale (mm) -> heat t
            case 'mexhat'
                lamS = (2*pi/(pScaleMM/1000))^2; kp.t = 1/max(lamS,eps);
        end
    end
    function Show()
        set(hPatch,'FaceVertexCData',W(:,curFrame),'FaceColor','interp');
        m=max(abs(W(:,curFrame)))+eps; clim(hAxes,[-m m]);
        set(hLabel,'String',sprintf('[%s] %s @ vtx %d | frame %d/%d  t=%.2fs   [arrows step time]', state, kernel, seedVtx, curFrame, ax.nT, ax.tlag(curFrame)));
    end
    function SaveAtom(~,~)
        if isempty(W), set(hLabel,'String','[design] nothing to save - drop an atom first.'); return; end
        state='save';
        [sSubj]=bst_get('SurfaceFile', SurfaceFile);
        T=bst_dynamics('New', sprintf('atom %s @vtx%d', kernel, seedVtx)); T.SurfaceFile=SurfaceFile;
        G=bst_dynamics('NewGroup','atom'); G.vertices=seedVtx; G.pos=V(seedVtx,:);
        LS=bst_atom('Levelset', W(ax.GlobalVertices{1},:), ax.GlobalVertices{1}, 0.5);
        G.region={LS.scoutVertices(:)'}; G.times=[ax.tlag(LS.eventSamples(1)); ax.tlag(LS.eventSamples(end))];
        T=bst_dynamics('AddGroup',T,G);
        outDir=bst_fileparts(file_fullpath(sSubj.FileName));
        of=file_unique(bst_fullfile(outDir,'dynamics_atom.mat')); bst_dynamics('Save', of, T);
        state='design';
        set(hLabel,'String',sprintf('[saved] atom -> %s (Scout %d vtx, Event %.2f-%.2fs)', file_short(of), numel(LS.scoutVertices), G.times(1), G.times(2)));
    end
end

% ===== helpers =====
function i_call(fcn,h,ev)
    if isempty(fcn), return; end
    if isa(fcn,'function_handle'), fcn(h,ev);
    elseif iscell(fcn), feval(fcn{1},h,ev,fcn{2:end});
    elseif ischar(fcn), eval(fcn); end %#ok<EVLDR>
end
function s = i_en(tf),  if tf, s='on'; else, s='off'; end, end
function s = i_en2b(show,c), if show, s=sprintf(' | speed %.2g m/s', c); else, s=''; end, end
function [hS,hV] = i_slider(p, name, y, lw, sw, vw, mn, mx, v0, cb)
    uicontrol(p,'Style','text','String',name,'Units','normalized','Position',[0.04 y lw 0.12],'HorizontalAlignment','left');
    hS = uicontrol(p,'Style','slider','Min',mn,'Max',mx,'Value',max(min(v0,mx),mn),'Units','normalized','Position',[0.35 y+0.01 sw 0.10],'Callback',cb);
    hV = uicontrol(p,'Style','text','String',num2str(round(v0,2)),'Units','normalized','Position',[0.80 y vw 0.12],'HorizontalAlignment','left');
end
function L = i_mean_edge(V,F)
    e=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])]; L=mean(sqrt(sum((V(e(:,1),:)-V(e(:,2),:)).^2,2)));
end
