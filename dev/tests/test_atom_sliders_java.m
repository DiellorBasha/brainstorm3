% test_atom_sliders_java - panel_eigenfilter_design physical-units atom sliders, build/set/read round-trip
import javax.swing.*;
b = struct('scaleMinMM',7,'scaleMaxMM',95,'rateMinMM2',49,'rateMaxMM2',9025);
% gabor: 3 active slots (Scale mm, Freq Hz, BW Hz)
jP = JPanel();
panel_eigenfilter_design('BuildAtomSliders', jP, 'gabor', b, []);
panel_eigenfilter_design('SetAtomVals', jP, [30 12 5]);
v = panel_eigenfilter_design('ReadAtomVals', jP);
assert(abs(v(1)-30)<0.2 && abs(v(2)-12)<0.2 && abs(v(3)-5)<0.2, sprintf('gabor round-trip [%.2f %.2f %.2f]', v));
% resonator: slots 1,2 active (Freq, Q); slot 3 disabled -> reads 0
jR = JPanel();
panel_eigenfilter_design('BuildAtomSliders', jR, 'resonator', b, []);
panel_eigenfilter_design('SetAtomVals', jR, [10 6 0]);
vr = panel_eigenfilter_design('ReadAtomVals', jR);
assert(abs(vr(1)-10)<0.2 && abs(vr(2)-6)<0.2, sprintf('resonator round-trip [%.2f %.2f]', vr(1), vr(2)));
assert(vr(3)==0, 'disabled slot reads 0');
% the kernel list includes the js kernels and is diffusion-first
[keys, displays] = panel_eigenfilter_design('AtomKernels');
assert(strcmp(keys{1},'diffusion'), 'diffusion is first');
assert(all(ismember({'gabor','travwave','resonator','stmatern'}, keys)), 'js kernels present');
assert(numel(displays)==numel(keys), 'a display name per kernel');
disp('OK');
