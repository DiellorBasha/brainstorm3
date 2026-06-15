function test_filterbank_schema()
% Schema regression for the filterbankmat template and the Surface.Filterbank list.
% Authors: Diellor Basha, 2026
    nFail = 0;
    fb = db_template('filterbankmat');
    need = {'Comment','ParentEigen','Variant','Tiles','Tiling','Provenance'};
    for f = need
        if ~isfield(fb, f{1}); fprintf('MISSING filterbankmat.%s\n', f{1}); nFail = nFail+1; end
    end
    s = db_template('surface');
    if ~isfield(s, 'Filterbank'); fprintf('MISSING surface.Filterbank\n'); nFail = nFail+1; end
    if isfield(s,'Filterbank')
        sub = fieldnames(s.Filterbank);
        for f = {'FileName','Comment','ParentEigen'}
            if ~ismember(f{1}, sub); fprintf('MISSING surface.Filterbank.%s\n', f{1}); nFail = nFail+1; end
        end
        if ~isempty(s.Filterbank); fprintf('surface.Filterbank must start EMPTY (0x0 struct)\n'); nFail = nFail+1; end
    end
    fprintf('\n==== test_filterbank_schema: %d failed ====\n', nFail);
    if nFail > 0, error('test_filterbank_schema FAILED'); end
end
