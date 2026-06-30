% test_eigfilter_controls - Sliders config + ToKernel mapping match the designer's prior behavior
b = struct('scaleMinMM',7,'scaleMaxMM',95,'rateMinMM2',49,'rateMaxMM2',9025);
lmax = 40;
% --- Sliders: labels per kernel ---
g = bst_eigfilter_controls('Sliders','gabor',b);
assert(strcmp(g(1).label,'Scale (mm)') && strcmp(g(2).label,'Freq (Hz)') && strcmp(g(3).label,'BW (Hz)'), 'gabor rows');
r = bst_eigfilter_controls('Sliders','resonator',b);
assert(strcmp(r(1).label,'Freq (Hz)') && strcmp(r(2).label,'Q') && isempty(r(3).label), 'resonator rows');
d = bst_eigfilter_controls('Sliders','diffusion',b);
assert(strcmp(d(1).label,'Rate (mm^2/s)') && d(1).lo==49 && d(1).hi==9025, 'diffusion rate range from bounds');
% --- ToKernel: matches the designer's i_phys2kernel formulas ---
kp = bst_eigfilter_controls('ToKernel','gabor',[30 12 2],lmax);
assert(abs(kp.k0 - 2*pi/0.03)<1e-9 && kp.f0==12 && kp.sf==2 && kp.lmax==lmax, 'gabor ToKernel');
kp = bst_eigfilter_controls('ToKernel','wave',[0 4 0],lmax);
assert(abs(kp.alpha - 4*sqrt(lmax)/2)<1e-12, 'wave alpha');
kp = bst_eigfilter_controls('ToKernel','dampedwave',[0 4 0.5],lmax);
assert(abs(kp.beta - 1/0.5)<1e-12, 'dampedwave beta');
kp = bst_eigfilter_controls('ToKernel','diffusion',[400 0 0],lmax);
assert(abs(kp.tau - max((400/1e6)*lmax,eps))<1e-15, 'diffusion tau');
kp = bst_eigfilter_controls('ToKernel','resonator',[10 6 0],lmax);
assert(kp.f0==10 && kp.Q==6, 'resonator');
kp = bst_eigfilter_controls('ToKernel','stmatern',[50 1.5 0],lmax);
assert(abs(kp.kappa - 2*pi/0.05)<1e-9 && kp.nu==1.5, 'stmatern');
disp('OK');
