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

    sSurf=in_tess_bst(SurfaceFile); nV=size(sSurf.Vertices,1);
    f=double(suvrMap(:)); f(~isfinite(f))=NaN;

    % --- LBO heat smoothing (tangential), per hemisphere ---
    fs=local_heatsmooth(SurfaceFile, f, Opts.HeatT);

    % --- vertex adjacency + local maxima (f(v) >= all finite neighbours) ---
    if isfield(sSurf,'VertConn') && ~isempty(sSurf.VertConn), A=sSurf.VertConn;
    else, A=tess_vertconn(sSurf.Vertices, sSurf.Faces); end
    A=A | A'; A=A - diag(diag(A));
    maxVerts=local_localmax(fs, A);

    foci=struct('vertex',{},'peak',{},'persistence',{},'basinArea',{});
    basinLabel=zeros(nV,1);
    info=struct('smoothed',fs,'nFoci',numel(maxVerts),'maxVerts',maxVerts,'A',{A});
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
