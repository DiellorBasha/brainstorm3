% test_levelset - Levelset on a synthetic wavelet field, via bst_dynamics (moved from bst_atom)
nV = 40; nT = 20; gv = (1:nV)';
tp = exp(-(((1:nT)-10).^2)/4);                 % temporal bump peaked at t=10
W = zeros(nV,nT);
W(10,:) = 1.0*tp;  W(11,:) = 0.6*tp;  W(12,:) = 0.3*tp;   % spatial peak at vertex 10
LS = bst_dynamics('Levelset', W, gv, 0.5);
% iRef = peak-energy frame = 10; wRef = W(:,10): v10=1, v11=0.6, v12=0.3, rest 0. threshold 0.5*max=0.5.
assert(LS.iRef == 10,                       'peak-energy frame');
assert(ismember(10, LS.scoutVertices),      'peak vertex in scout');
assert(ismember(11, LS.scoutVertices),      'v11 (0.6) in scout');
assert(~ismember(12, LS.scoutVertices),     'v12 (0.3 < 0.5) excluded');
assert(~ismember(40, LS.scoutVertices),     'far/zero vertex excluded');
assert(~isempty(LS.eventSamples),           'event samples found');
disp('OK');
