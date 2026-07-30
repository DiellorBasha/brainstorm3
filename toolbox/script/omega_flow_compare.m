function C = omega_flow_compare()
% OMEGA_FLOW_COMPARE  Compare the Dirac-EIGENBASIS group flow maps against the basis-free/vertex
% (wMNE) run on omega_tutorial_test, per flow quantity x band. Validates the eigenbasis migration:
% the two source kernels should give physiologically CONCORDANT flow (high spatial correlation,
% same posterior-alpha structure) -- the eigenbasis version being the regularized one.
%
% Locates the group averages by their tags (groupflow_<Q>_<mode>, written by tutorial_flow_omega),
% loads each band map on the group template, and reports the per-band spatial correlation. Also
% re-checks the flow!=power control (div vs curl spatial correlation) on the eigenbasis maps.
%
% OUTPUT C: struct with .quantities, .bands, .corr [nQ x nBand] (dirac-vs-mn), .divCurl [nBand x 1].
%
% Author: Diellor Basha, 2026

    gui_brainstorm('SetCurrentProtocol', bst_get('Protocol', 'omega_tutorial_test'));
    quantities = {'div','curl','psi'};
    C = struct('quantities', {quantities}, 'bands', {{}}, 'corr', [], 'divCurl', []);

    maps = struct();   % maps.(mode).(Q) = [nVert x nBand]
    for mode = {'dirac','mn'}
        for iq = 1:numel(quantities)
            Q = quantities{iq};
            f = i_find_group(['groupflow_' Q '_' mode{1}]);
            if isempty(f), fprintf('  [%s/%s] group average NOT found\n', mode{1}, Q); continue; end
            R = in_bst_results(f, 1);
            M = R.ImageGridAmp;                       % [nVert x 1 x nBand] banded PSD on group template
            maps.(mode{1}).(Q) = squeeze(M(:,1,:));   % [nVert x nBand]
            if isempty(C.bands) && isfield(R,'Freqs') && iscell(R.Freqs), C.bands = R.Freqs(:,1)'; end
        end
    end

    nB = size(maps.dirac.(quantities{1}), 2);
    C.corr = nan(numel(quantities), nB);
    fprintf('\n=== eigenbasis(dirac) vs basis-free(mn) group flow -- per-band spatial correlation ===\n');
    for iq = 1:numel(quantities)
        Q = quantities{iq};
        if ~isfield(maps,'dirac')||~isfield(maps.dirac,Q)||~isfield(maps,'mn')||~isfield(maps.mn,Q), continue; end
        for b = 1:nB
            C.corr(iq,b) = corr(maps.dirac.(Q)(:,b), maps.mn.(Q)(:,b));
        end
        fprintf('  %-5s: %s\n', Q, sprintf('%.3f ', C.corr(iq,:)));
    end
    if ~isempty(C.bands), fprintf('  bands: %s\n', strjoin(C.bands, ' ')); end

    % flow != power control on the EIGENBASIS maps: div vs curl spatial correlation per band
    if isfield(maps,'dirac') && isfield(maps.dirac,'div') && isfield(maps.dirac,'curl')
        C.divCurl = arrayfun(@(b) corr(maps.dirac.div(:,b), maps.dirac.curl(:,b)), 1:nB)';
        fprintf('  div-curl corr (eigenbasis, flow!=power control): %s\n', sprintf('%.3f ', C.divCurl));
    end
end

% newest group average whose Comment contains the tag
function f = i_find_group(tag)
    f = '';
    sStudies = bst_get('ProtocolStudies');
    best = -inf;
    for is = 1:numel(sStudies.Study)
        for ir = 1:numel(sStudies.Study(is).Result)
            R = sStudies.Study(is).Result(ir);
            if ~isempty(R.FileName) && contains(R.Comment, tag)
                d = dir(file_fullpath(R.FileName));
                if ~isempty(d) && d.datenum > best, best = d.datenum; f = R.FileName; end
            end
        end
    end
end
