function [est, dinfo] = pet_spread_invert(SurfaceFile, a, tau, dt, Opts)
% PET_SPREAD_INVERT: recover the reaction-diffusion parameters {ra,Da,rt,Dt,kappa} from a
% longitudinal Abeta/tau surface series by mass-weighted linear regression (the coupled Fisher-KPP
% is linear in its parameters). Mirrors the implicit-Euler discretization of pet_spread_simulate:
%   M*(x_{n+1}-x_n)/dt = M*[reaction basis]_n + beta_D*(-K*x_{n+1})
% Abeta:  M*da/dt   = ra*M*(a(1-a))                                  - Da*K*a_{n+1}
% tau:    M*dtau/dt = rt*M*(tau(1-tau)) + (kappa*rt)*M*(a*tau(1-tau)) - Dt*K*tau_{n+1}
% kappa = beta2/beta1.
%
% USAGE: [est, dinfo] = pet_spread_invert(SurfaceFile, a, tau, dt, Opts)
%
% INPUTS:
%   - SurfaceFile : cortex surface file (for the LBO M,K; same surface as the simulator).
%   - a, tau      : [nL x nT] Abeta / tau LH field time-series.
%   - dt          : timestep (default 1).
%   - Opts        : .active=[0.05 0.95] (fit only where the field is in this band), .ridge=0,
%                   .diffAt='np1' (diffusion term at step n+1, matching the implicit simulator).
%
% OUTPUTS:
%   - est   : struct .ra .Da .rt .Dt .kappa.
%   - dinfo : .condX (design-matrix condition number) .resid (relative residual) .nSamp.
%
% SEE ALSO: pet_spread_simulate, tess_operators
%
% Author: Diellor Basha, 2026
    D=struct('active',[0.05 0.95],'ridge',0,'diffAt','np1');
    if nargin<5, Opts=struct(); end
    fn=fieldnames(D); for i=1:numel(fn), if ~isfield(Opts,fn{i}), Opts.(fn{i})=D.(fn{i}); end; end
    if nargin<4||isempty(dt), dt=1; end
    LBO=tess_operators(SurfaceFile,'Laplace-Beltrami'); M=LBO.Mass{1}; K=LBO.Operator{1};

    % ---- Abeta fit: [ra, Da] ----
    [ya,Xa]=local_assemble(a, [], M, K, dt, Opts, false);
    ba=local_solve(Xa, ya, Opts.ridge);
    % ---- tau fit: [rt, kappa*rt, Dt] (coupling basis uses a at step n) ----
    [yt,Xt]=local_assemble(tau, a, M, K, dt, Opts, true);
    bt=local_solve(Xt, yt, Opts.ridge);

    est=struct('ra',ba(1),'Da',ba(2),'rt',bt(1),'Dt',bt(3),'kappa',bt(2)/bt(1));
    dinfo=struct('condX',cond(Xt'*Xt),'resid',norm(Xt*bt-yt)/max(norm(yt),eps),'nSamp',numel(yt));
end

function [y,X]=local_assemble(f, g, M, K, dt, Opts, coupled)
    nT=size(f,2); y=[]; X=[];
    for n=1:nT-1
        fn=f(:,n); fnp=f(:,n+1);
        act = fn>Opts.active(1) & fn<Opts.active(2);     % fit only in the active (non-saturated) band
        if ~any(act), continue; end
        yi  = M*((fnp-fn)/dt);
        x1  = M*(fn.*(1-fn));
        if strcmpi(Opts.diffAt,'np1'), xd=-K*fnp; else, xd=-K*fn; end
        if coupled
            x2 = M*(g(:,n).*fn.*(1-fn));
            Xi = [x1(act), x2(act), xd(act)];
        else
            Xi = [x1(act), xd(act)];
        end
        y=[y; yi(act)]; X=[X; Xi]; %#ok<AGROW>
    end
end

function b=local_solve(X, y, ridge)
    if ridge>0, b=(X'*X + ridge*eye(size(X,2))) \ (X'*y); else, b=X\y; end
end
