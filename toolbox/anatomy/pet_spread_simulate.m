function [a, tau, info] = pet_spread_simulate(SurfaceFile, Opts)
% PET_SPREAD_SIMULATE: synthetic coupled Abeta/tau spread on the cortical manifold (LH).
% Two Fisher-KPP reaction-diffusion fields, LBO diffusion (implicit Euler) + logistic reaction;
% Abeta seeds at precuneus, tau at entorhinal with growth gated by local Abeta (the cascade).
%
% USAGE: [a, tau, info] = pet_spread_simulate(SurfaceFile, Opts)
%
% INPUTS:
%   - SurfaceFile : Brainstorm cortex surface file (string).
%   - Opts        : (optional) struct; defaults: nT=24, dt=1, Da=8e-5, Dt=8e-5, ra=0.6, rt=0.45,
%                   kappa=3, seedAmp=0.6, seedA=[], seedT=[] (LH-local seed indices; if empty,
%                   derived from Desikan 'precuneus L' / 'entorhinal L').
%
% OUTPUTS:
%   - a, tau : [nL x nT] LH field time-series (Abeta, tau), values in [0,1].
%   - info   : struct .gv (LH global vertex indices) .M .K (LH LBO operators) .seedA .seedT .Opts.
%
% SEE ALSO: pet_spread_invert, tess_operators
%
% Author: Diellor Basha, 2026
    D=struct('nT',24,'dt',1,'Da',8e-5,'Dt',8e-5,'ra',0.6,'rt',0.45,'kappa',3,'seedAmp',0.6,'seedA',[],'seedT',[]);
    if nargin<2, Opts=struct(); end
    fn=fieldnames(D); for i=1:numel(fn), if ~isfield(Opts,fn{i}), Opts.(fn{i})=D.(fn{i}); end; end

    sW=in_tess_bst(SurfaceFile); nV=size(sW.Vertices,1);
    LBO=tess_operators(SurfaceFile,'Laplace-Beltrami');
    gv=double(LBO.GlobalVertices{1}(:)); M=LBO.Mass{1}; K=LBO.Operator{1}; nL=numel(gv);
    g2l=zeros(nV,1); g2l(gv)=1:nL;
    sa=Opts.seedA; if isempty(sa), sa=local_seed(sW,g2l,'precuneus L'); end
    st=Opts.seedT; if isempty(st), st=local_seed(sW,g2l,'entorhinal L'); end

    dA=decomposition(M + Opts.dt*Opts.Da*K);
    dT=decomposition(M + Opts.dt*Opts.Dt*K);
    a=zeros(nL,1); a(sa)=Opts.seedAmp; tau=zeros(nL,1); tau(st)=Opts.seedAmp;
    A=zeros(nL,Opts.nT); T=zeros(nL,Opts.nT);
    for n=1:Opts.nT
        ra = Opts.ra * a   .* (1-a);
        rt = Opts.rt * (1+Opts.kappa*a) .* tau .* (1-tau);
        a   = dA\(M*(a   + Opts.dt*ra));  a  =min(max(a,0),1);
        tau = dT\(M*(tau + Opts.dt*rt));  tau=min(max(tau,0),1);
        A(:,n)=a; T(:,n)=tau;
    end
    a=A; tau=T;
    info=struct('gv',gv,'M',M,'K',K,'seedA',sa,'seedT',st,'Opts',Opts);
end

function loc=local_seed(sW, g2l, label)
    ai=find(~cellfun('isempty',regexp({sW.Atlas.Name},'Desikan','once')),1);
    si=find(strcmp({sW.Atlas(ai).Scouts.Label}, label),1);
    gvs=sW.Atlas(ai).Scouts(si).Vertices(:);
    s=sW.Atlas(ai).Scouts(si).Seed; if isempty(s), s=gvs(round(numel(gvs)/2)); end
    loc=g2l(s); if loc==0, loc=g2l(gvs(find(g2l(gvs)>0,1))); end
end
