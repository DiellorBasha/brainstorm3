function [foci, basinLabel, info] = pet_epicenter(SurfaceFile, suvrMap, Opts)
% PET_EPICENTER: cortical concentration foci (epicenters) of a surface SUVR map via discrete
% Morse-Smale geometry. Heat-smooth (LBO) -> local maxima (mesh adjacency) -> gradient-ascent
% basins -> persistence-ranked dominant + secondary foci.
%
% USAGE: [foci, basinLabel, info] = pet_epicenter(SurfaceFile, suvrMap, Opts)
%
% Author: Diellor Basha, 2026
    if (nargin<3)||isempty(Opts), Opts=struct(); end
    Def=struct('HeatT',2e-5,'MinPersist',[],'nFociMax',10);
    fn=fieldnames(Def); for i=1:numel(fn), if ~isfield(Opts,fn{i}), Opts.(fn{i})=Def.(fn{i}); end; end

    sSurf=in_tess_bst(SurfaceFile);
    f=double(suvrMap(:)); f(~isfinite(f))=NaN;

    % --- LBO heat smoothing (tangential), per hemisphere ---
    fs=local_heatsmooth(SurfaceFile, f, Opts.HeatT);

    % --- vertex adjacency + local maxima (f(v) >= all finite neighbours) ---
    if isfield(sSurf,'VertConn') && ~isempty(sSurf.VertConn), A=sSurf.VertConn;
    else, A=tess_vertconn(sSurf.Vertices, sSurf.Faces); end
    A=A | A'; A=A - diag(diag(A));
    maxVerts=local_localmax(fs, A);

    % --- gradient-ascent basins of attraction (Morse-Smale segmentation) ---
    basinLabel=local_basins(fs, A, maxVerts);

    % --- persistence of each maximum (super-level merge tree) -> rank + filter ---
    pers=local_persistence(fs, A, maxVerts);
    if isempty(Opts.MinPersist)
        rng=max(fs(isfinite(fs)))-min(fs(isfinite(fs))); thr=0.15*rng;
    else
        thr=Opts.MinPersist;
    end
    keep=pers>=thr; mk=maxVerts(keep); pk=pers(keep);
    % rank by persistence; break ties by peak height (each disconnected hemisphere has one
    % essential, full-range-persistence max -> peak distinguishes the true epicenter).
    [~,ord]=sortrows([pk(:), fs(mk(:))],'descend');
    mk=mk(ord); pk=pk(ord);
    if numel(mk)>Opts.nFociMax, mk=mk(1:Opts.nFociMax); pk=pk(1:Opts.nFociMax); end
    foci=struct('vertex',{},'peak',{},'persistence',{},'basinArea',{});
    for i=1:numel(mk)
        bi=find(maxVerts==mk(i),1);
        foci(i)=struct('vertex',mk(i),'peak',fs(mk(i)),'persistence',pk(i),'basinArea',nnz(basinLabel==bi));
    end

    info=struct('smoothed',fs,'nFoci',numel(foci),'maxVerts',maxVerts,'A',{A});
end

% ===== persistence: peak height minus the level at which its basin merges into a higher one =====
function pers=local_persistence(f, A, maxVerts)
    n=numel(f); [~,order]=sort(f,'descend'); order=order(isfinite(f(order)));
    comp=zeros(n,1); peakOf=zeros(n,1);                 % union-find parent + each root's peak vertex
    isMax=false(n,1); isMax(maxVerts)=true;
    pers=inf(numel(maxVerts),1); maxId=zeros(n,1); maxId(maxVerts)=1:numel(maxVerts);
    for oi=1:numel(order)
        v=order(oi);
        nbUp=find(A(:,v)); nbUp=nbUp(comp(nbUp)>0);     % already-added (higher-or-equal) neighbours
        roots=[]; for x=nbUp(:)', roots(end+1)=local_find(comp,x); end %#ok<AGROW>
        roots=unique(roots);
        if isempty(roots)
            comp(v)=v; peakOf(v)=v;                     % new component born at a maximum
        else
            [~,hi]=max(f(arrayfun(@(r)peakOf(r),roots))); keepRoot=roots(hi);
            comp(v)=keepRoot;
            for r=roots(:)'
                if r==keepRoot, continue; end
                pv=peakOf(r);
                if isMax(pv), pers(maxId(pv))=f(pv)-f(v); end   % lower peak dies at level f(v)
                comp(r)=keepRoot;                               % merge into the higher peak
            end
        end
    end
    pers(~isfinite(pers))=max(f(isfinite(f)))-min(f(isfinite(f)));  % the global max never dies
end
function r=local_find(comp,x)
    r=x; while comp(r)~=r && comp(r)~=0, r=comp(r); end
end

% ===== LBO heat smoothing of a vertex field (per hemisphere) =====
function g=local_heatsmooth(SurfaceFile, f, t)
    LBO=tess_operators(SurfaceFile,'Laplace-Beltrami');
    g=f;
    for hh=1:numel(LBO.Operator)
        gv=double(LBO.GlobalVertices{hh}(:)); M=LBO.Mass{hh}; K=LBO.Operator{hh};
        fh=f(gv); fh(~isfinite(fh))=0;
        g(gv)=(M + t*K) \ (M*fh);
    end
end

% ===== local maxima over the adjacency: f(v) >= every neighbour (ties resolved by index) =====
function mv=local_localmax(f, A)
    n=numel(f); mv=[];
    for v=1:n
        if ~isfinite(f(v)), continue; end
        nb=find(A(:,v)); nb=nb(isfinite(f(nb)));
        if isempty(nb) || all(f(v) > f(nb)) || (all(f(v) >= f(nb)) && v < min(nb(f(nb)==f(v))))
            mv(end+1)=v; %#ok<AGROW>
        end
    end
    mv=mv(:);
end

% ===== steepest-ascent basins: each vertex flows to the max it climbs to =====
function lbl=local_basins(f, A, maxVerts)
    n=numel(f); nextUp=zeros(n,1);
    for v=1:n
        if ~isfinite(f(v)), continue; end
        nb=find(A(:,v)); nb=nb(isfinite(f(nb)));
        if isempty(nb), nextUp(v)=v; continue; end
        [mx,k]=max(f(nb));
        if mx<=f(v), nextUp(v)=v; else, nextUp(v)=nb(k); end
    end
    maxIdx=zeros(n,1); maxIdx(maxVerts)=1:numel(maxVerts);
    lbl=zeros(n,1);
    for v=1:n
        if ~isfinite(f(v)), continue; end
        p=v;
        while nextUp(p)~=p, p=nextUp(p); end   % climb to the peak this vertex flows to
        lbl(v)=maxIdx(p);
    end
end
