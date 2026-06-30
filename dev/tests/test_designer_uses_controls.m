% test_designer_uses_controls - the designer delegates param mapping to the shared module and reads its
% slider values from the docked panel (sliders themselves are rendered by panel_atom_designer, not on the figure)
src = fileread('/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox/gui/view_atom_designer.m');
assert(contains(src, "bst_eigfilter_controls('ToKernel'"), 'i_phys2kernel delegates to shared ToKernel');
assert(contains(src, "panel_atom_designer('ReadVals')"), 'param values read from the docked panel');
% the per-kernel switch tables are gone from the designer (centralized in bst_eigfilter_controls)
assert(~contains(src, "case 'travwave',    kp.c"), 'i_phys2kernel switch removed');
% slider rendering no longer lives on the figure (moved to the panel)
assert(~contains(src, 'function i_config_sliders'), 'on-figure slider rendering removed');
disp('OK');
