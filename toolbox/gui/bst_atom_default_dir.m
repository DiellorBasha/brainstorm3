function dir = bst_atom_default_dir(ax, seedVert)
% BST_ATOM_DEFAULT_DIR  Default impulse direction for an atom's operator fiber (app-side).
%   dir = bst_atom_default_dir(ax, seedVert)
%   Quaternion (Dirac) fiber -> the seed vertex surface normal (unit row 3-vector).
%   Tangent/scalar fiber     -> 1 (frame e1 / scalar amplitude).
% App-side default: orchestrators (the library realisers) carry no defaults; the GUI supplies them.
% Shared by view_atom_designer and (SP2b) panel_bst_dynamics. See atom-operator-applicability.
%
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% =============================================================================@

    [~, kind] = bst_eigenfilter('Fiber', ax);
    if ~strcmp(kind, 'quaternion'), dir = 1; return; end
    dir = [0 0 1];                                       % fallback if normals are missing
    try
        S = in_tess_bst(ax.SurfaceFile, 0);
        if isfield(S,'VertNormals') && ~isempty(S.VertNormals) && seedVert <= size(S.VertNormals,1)
            n = S.VertNormals(seedVert, :);  if norm(n) > 0, dir = n / norm(n); end
        end
    catch %#ok<CTCH>
    end
    dir = dir(:)';
end
