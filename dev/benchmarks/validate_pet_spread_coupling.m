function R = validate_pet_spread_coupling()
% VALIDATE_PET_SPREAD_COUPLING: kappa-recovery curve + noise/sparse-time robustness for the
% reaction-diffusion coupling inverter, on synthetic ground truth.
%
% USAGE: R = validate_pet_spread_coupling()
%
% Author: Diellor Basha, 2026
    [sS,~]=bst_get('Subject','sub-MTL0002');
    wf=sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
    here=bst_fileparts(mfilename('fullpath'));
    % Well-sampled base (dt=0.5): the inverter mirrors the simulator's implicit step, so adequate
    % temporal resolution is required for identifiability (coarse dt undersamples fast dynamics).
    nT=48; dt=0.5; sub=4;
    kts=[0 1 3 5]; kEst=zeros(size(kts)); kNoise=zeros(size(kts)); kSparse=zeros(size(kts));
    for i=1:numel(kts)
        [a,tau]=pet_spread_simulate(wf, struct('kappa',kts(i),'nT',nT,'dt',dt));
        e=pet_spread_invert(wf,a,tau,dt);                 kEst(i)=e.kappa;
        aN=min(max(a+0.03*local_pn(size(a),1),0),1);      tN=min(max(tau+0.03*local_pn(size(tau),2),0),1);
        e=pet_spread_invert(wf,aN,tN,dt);                 kNoise(i)=e.kappa;
        idx=1:sub:size(a,2);                              % sparse: coarse temporal sampling
        e=pet_spread_invert(wf,a(:,idx),tau(:,idx),dt*sub); kSparse(i)=e.kappa;
    end
    fprintf('\nkappa recovery:\n  true:   %s\n  full:   %s\n  +noise: %s\n  sparse: %s\n', ...
        mat2str(kts), mat2str(round(kEst,2)), mat2str(round(kNoise,2)), mat2str(round(kSparse,2)));
    f=figure('Visible','off','Position',[40 40 560 480]);
    plot(kts,kts,'k--'); hold on;
    plot(kts,kEst,'o-','LineWidth',1.5); plot(kts,kNoise,'s-'); plot(kts,kSparse,'^-');
    xlabel('true \kappa'); ylabel('recovered \kappa'); grid on; axis equal;
    legend({'identity','full series (dt=0.5)','+3% noise','sparse (coarse dt, breaks)'},'Location','northwest');
    title('A\beta\rightarrowtau coupling recovery');
    png=fullfile(here,'pet_spread_kappa_recovery.png'); print(f,png,'-dpng','-r120'); close(f);
    fprintf('  figure -> %s\n', png);
    R=struct('kts',kts,'kEst',kEst,'kNoise',kNoise,'kSparse',kSparse);
end

function r=local_pn(sz, seed)
    % deterministic pseudo-noise in [-1,1] (no randn/global-rng dependence): hashed sinusoid.
    n=prod(sz); idx=(1:n)';
    r=sin(idx*12.9898 + seed*7.0)*43758.5453; r=r-floor(r); r=(r-0.5)*2;
    r=reshape(r,sz);
end
