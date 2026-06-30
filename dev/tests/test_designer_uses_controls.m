% test_designer_uses_controls - the designer delegates its control spec/mapping to the shared module
src = fileread('/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox/gui/view_atom_designer.m');
assert(contains(src, "bst_eigfilter_controls('Sliders'"), 'i_config_sliders delegates to shared Sliders');
assert(contains(src, "bst_eigfilter_controls('ToKernel'"), 'i_phys2kernel delegates to shared ToKernel');
% and the per-kernel switch tables are gone from the designer (now centralized)
assert(~contains(src, "case 'travwave',    kp.c"), 'i_phys2kernel switch removed');
disp('OK');
