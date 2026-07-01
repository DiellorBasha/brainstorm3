function tests = test_itersine_kernel
tests = functiontests(localfunctions);
end

function test_meta(tc)
m = bst_eigfilter_design_itersine('meta');
verifyEqual(tc, m.name, 'itersine');
verifyTrue(tc, m.bandpass);
end

function test_single_member_handle(tc)
g = bst_eigfilter_design_itersine(struct('member',2,'Nf',6,'lmax',10));
verifyTrue(tc, isa(g,'function_handle'));
y = g(linspace(0,10,50)');
verifyEqual(tc, numel(y), 50);
verifyGreaterThanOrEqual(tc, min(y), -1e-12);      % itersine window is non-negative
end

function test_tight_frame_property(tc)
% Sum of squares of all Nf members is ~constant over the interior -> tight (B/A ~ 1)
Nf = 6; lmax = 10;
lam = linspace(0.1*lmax, 0.9*lmax, 400)';          % interior (edges roll off)
S = zeros(size(lam));
for ii = 1:Nf
    g = bst_eigfilter_design_itersine(struct('member',ii,'Nf',Nf,'lmax',lmax));
    S = S + g(lam).^2;
end
A = min(S); B = max(S);
verifyLessThan(tc, B/A, 1.05);                     % tight frame on the interior
end
