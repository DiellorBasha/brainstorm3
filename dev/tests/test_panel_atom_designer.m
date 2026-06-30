% test_panel_atom_designer - the panel builds, docks, configures, and reads back its controls
bstPanel = panel_atom_designer('CreatePanel');
assert(~isempty(bstPanel), 'CreatePanel returns a BstPanel');
gui_show(bstPanel, 'BrainstormTab', 'tools');     % register so bst_get('PanelControls') resolves
try
    b  = struct('scaleMinMM',7,'scaleMaxMM',95,'rateMinMM2',49,'rateMaxMM2',9025);
    cb = struct('Kernel',@()[], 'Operator',@()[], 'Param',@()[], 'Fibers',@(s)[], 'Save',@()[]);
    panel_atom_designer('Configure', cb, b, 'gabor', 'Laplace-Beltrami');
    assert(strcmp(panel_atom_designer('CurrentKernel'),'gabor'), 'Filter set to gabor');
    assert(strcmp(panel_atom_designer('CurrentOperator'),'Laplace-Beltrami'), 'Operator = geometric');
    assert(numel(panel_atom_designer('ReadVals'))==3, 'ReadVals returns [s1 s2 s3]');
    panel_atom_designer('RebuildSliders', 'resonator', b);
    assert(numel(panel_atom_designer('ReadVals'))==3, 'ReadVals after rebuild');
    panel_atom_designer('SetStatus', 'designer status');
    ok = true;
catch e
    ok = false;  msg = e.message;
end
gui_hide('AtomDesigner');
if ok, disp('OK'); else, error('test_panel_atom_designer FAILED: %s', msg); end
