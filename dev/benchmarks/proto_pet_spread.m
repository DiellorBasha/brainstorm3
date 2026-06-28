function [A, Tau, info] = proto_pet_spread(Opts)
% PROTO_PET_SPREAD: prototype coupled Abeta/tau spread on the real cortex (exploration, not final).
% Two Fisher-KPP reaction-diffusion fields on the LH cortical manifold (LBO diffusion via
% tess_operators): Abeta seeds at PRECUNEUS, tau seeds at ENTORHINAL with growth GATED by local
% Abeta (the amyloid->tau cascade). Implicit-Euler diffusion + explicit logistic reaction.
% Renders a montage of both tracers spreading over time (medial LH view).
%
% USAGE: [A,Tau,info] = proto_pet_spread()  /  proto_pet_spread(struct('kappa',0))  % uncoupled
%
% Author: Diellor Basha, 2026 (prototype)
    D=struct('subj','sub-MTL0002','nT',24,'dt',1,'Da',8e-5,'Dt',8e-5,'ra',0.6,'rt',0.45, ...
             'kappa',3,'seedAmp',0.6,'showT',[3 8 12 16 20 24]);
    if nargin<1, Opts=struct(); end
    fn=fieldnames(D); for i=1:numel(fn), if ~isfield(Opts,fn{i}), Opts.(fn{i})=D.(fn{i}); end; end
    here=bst_fileparts(mfilename('fullpath'));

    [sS,~]=bst_get('Subject',Opts.subj);
    wf=sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
    sW=in_tess_bst(wf); nV=size(sW.Vertices,1);
    LBO=tess_operators(wf,'Laplace-Beltrami');
    gv=double(LBO.GlobalVertices{1}(:)); M=LBO.Mass{1}; K=LBO.Operator{1}; nL=numel(gv);  % LH
    g2l=zeros(nV,1); g2l(gv)=1:nL;

    % --- seeds from the validated epicenters (Desikan LH scouts) ---
    sa=local_seed(sW, g2l, 'precuneus L');
    st=local_seed(sW, g2l, 'entorhinal L');
    A0=zeros(nL,1); A0(sa)=Opts.seedAmp;  Tau0=zeros(nL,1); Tau0(st)=Opts.seedAmp;

    % --- pre-factor the implicit-Euler diffusion operators ---
    dA=decomposition(M + Opts.dt*Opts.Da*K);
    dT=decomposition(M + Opts.dt*Opts.Dt*K);
    A=zeros(nL,Opts.nT); Tau=zeros(nL,Opts.nT); a=A0; tau=Tau0;
    for n=1:Opts.nT
        ra = Opts.ra * a   .* (1-a);                          % Fisher-KPP Abeta
        rt = Opts.rt * (1+Opts.kappa*a) .* tau .* (1-tau);    % tau, growth gated by local Abeta
        a   = dA\(M*(a   + Opts.dt*ra));   a  =min(max(a,0),1);
        tau = dT\(M*(tau + Opts.dt*rt));   tau=min(max(tau,0),1);
        A(:,n)=a; Tau(:,n)=tau;
    end
    info=struct('gv',gv,'seedA',sa,'seedT',st,'Opts',Opts);

    % --- montage: medial LH view, Abeta (top) + tau (bottom) at showT timepoints ---
    inLH=false(nV,1); inLH(gv)=true; keep=all(inLH(sW.Faces),2);
    rmp=zeros(nV,1); rmp(gv)=1:nL; Flh=rmp(sW.Faces(keep,:)); Vlh=sW.Vertices(gv,:);
    nc=numel(Opts.showT); lbls={'A\beta','tau'};
    f=figure('Visible','off','Position',[30 30 260*nc 520]);
    for c=1:nc
        for r=1:2
            subplot(2,nc,(r-1)*nc+c);
            if r==1, cd=A(:,Opts.showT(c)); else, cd=Tau(:,Opts.showT(c)); end
            patch('Faces',Flh,'Vertices',Vlh,'FaceVertexCData',cd,'FaceColor','interp','EdgeColor','none');
            clim([0 1]); colormap(hot); view([180 -10]); axis equal off vis3d; camlight headlight; lighting gouraud;
            if r==1, title(sprintf('t=%d',Opts.showT(c))); end
            if c==1, text(-0.18,0.5,lbls{r},'Units','normalized','FontWeight','bold','Rotation',90); end
        end
    end
    png=fullfile(here,sprintf('proto_pet_spread_kappa%g.png',Opts.kappa));
    print(f,png,'-dpng','-r110'); close(f);
    fprintf('coupled Abeta/tau spread (kappa=%g): %d steps on LH (%d vert). figure -> %s\n', Opts.kappa, Opts.nT, nL, png);
    info.png=png;
end

function loc=local_seed(sW, g2l, label)
    ai=find(~cellfun('isempty',regexp({sW.Atlas.Name},'Desikan','once')),1);
    si=find(strcmp({sW.Atlas(ai).Scouts.Label}, label),1);
    gvs=sW.Atlas(ai).Scouts(si).Vertices(:);
    s=sW.Atlas(ai).Scouts(si).Seed; if isempty(s), s=gvs(round(numel(gvs)/2)); end
    loc=g2l(s); if loc==0, loc=g2l(gvs(find(g2l(gvs)>0,1))); end   % map to LH-local index
end
