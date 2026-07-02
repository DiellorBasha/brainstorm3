function test_designer_operator_map
% Pure display<->variant maps for the 4-operator selector (no Swing needed).
    pairs = { 'Geometric','Laplace-Beltrami'; 'Connectomic','LB-Connectome'; ...
              'Dirac','Dirac'; 'Dirac (connectome)','Dirac-Connectome' };
    for i = 1:size(pairs,1)
        v = panel_atom_designer('i_variant_for_item', pairs{i,1});
        assert(strcmp(v, pairs{i,2}), 'item %s -> %s (got %s)', pairs{i,1}, pairs{i,2}, v);
        nm = panel_atom_designer('i_item_for_variant', pairs{i,2});
        assert(strcmp(nm, pairs{i,1}), 'variant %s -> %s (got %s)', pairs{i,2}, pairs{i,1}, nm);
    end
    % Unknown item defaults to Geometric (Laplace-Beltrami).
    assert(strcmp(panel_atom_designer('i_variant_for_item','???'), 'Laplace-Beltrami'));
    disp('test_designer_operator_map PASSED');
end
