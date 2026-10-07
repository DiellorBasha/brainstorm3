function R = bst_fibre_networks(SubjectName, OutputDir, Opts)
% BST_FIBRE_NETWORKS: PET-guided fibre-network measures of ONE multimodal Brainstorm subject.
%
% USAGE:  R = bst_fibre_networks(SubjectName, OutputDir, Opts)
%
% Reads, from the subject of the CURRENT protocol (a brainstorm-multimodal subject: T1 + volume
% atlases, cortex, the subject's own fibers tess_fibers_*, PET "<trc>_mean" and "<trc>_suvr"
% volumes and their cortical maps), and writes DERIVED per-subject measures only. Nothing here
% relates a measure to an outcome or a group (no association test).
%
%   1. Regional PET SUVR per tracer and method, on the Desikan-Killiany and Destrieux volume
%      atlases (cortex + aseg subcortical) and on the cortical scouts:
%        none  : static mean / plain cerebellar-cortex mean (ASEG 8+47), unsmoothed
%        vlpp  : the protocol's own SUVR volume (6 mm + cerebellar cortex; the validated product)
%        gtm   : Rousset GTM (pet_gtm, Desikan-Killiany partition, scanner PSF) / corrected
%                cerebellar cortex. Desikan-Killiany only: GTM is regional by construction.
%        rbv   : region-based voxelwise correction (Thomas 2011) from the GTM values,
%                PET .* G ./ (PSF * G), / RBV cerebellar cortex.
%        vlpp_surface : the protocol's cortical SUVR map averaged over the atlas scouts.
%   2. Streamline endpoints labelled on each volume atlas (grey-matter nodes: cortex + aseg
%      subcortical + cerebellar cortex; an end in white matter takes the nearest node voxel
%      within Opts.EndRadius mm). Per atlas: streamline-count and mean-length connectomes, plus
%      the protocol's own fibers_helper('AssignToScouts') connectome (nearest scout seed),
%      carried for comparison. Every end's nearest cortex vertex is kept, so any cortical
%      tiling (e.g. the dyadic tiles) rolls up exactly from <sub>_fn_endpoints.mat.
%   3. A-priori bundles per hemisphere (node-pair sets, LOCAL_BUNDLES): tau-epicentre-seeded
%      (B01-B07), controls (C01-C05) and epicentre->any (E01 entorhinal, E02 inferior temporal,
%      E03 the subject's own tau epicentre = max GTM flortaucipir SUVR among entorhinal,
%      parahippocampal, inferior temporal, fusiform, amygdala).
%   4. Along each bundle (oriented A->B): every sampled volume (PET none / vlpp / rbv per tracer
%      -- no PVC vs RBV is the spill-over contrast -- and any microstructure volume the subject
%      holds, Opts.ScalarPattern) as a per-streamline mean and an Opts.ProfilePoints profile.
%   5. Connectome-weighted PET per cortical node: CWA_s = sum_j W_sj A_j / sum_j W_sj over
%      cortical j ~= s, with comparators (seed-local, cortical mean, inverse-distance weighted).
%
% INPUTS:
%   SubjectName : subject name in the current protocol (e.g. 'sub-MTL0002')
%   OutputDir   : folder for the outputs (created)
%   Opts        : (optional) struct, fields (defaults):
%       .Tracers       {'18FNAV4694','18Fflortaucipir'}   matched on "PET <trc>_mean|_suvr"
%       .AmyloidTracer '18FNAV4694'    .TauTracer '18Fflortaucipir'
%       .PsfFwhm       []   mm; [] = pet_scanner_fwhm(PET metadata)
%       .RefLabels     [8 47]  ASEG cerebellar cortex (the protocol's own SUVR reference)
%       .Atlases       {'Desikan-Killiany','Destrieux'}
%       .EndRadius     3    mm, endpoint -> nearest grey-matter node voxel
%       .ProfilePoints 20
%       .ChunkSize     20000 streamlines per sampling pass
%       .ScalarPattern '^(DTI[ _-]?)?(FA|MD|RD|AD|AFD)\b'  microstructure volumes, if any
%       .Provenance    struct merged into <sub>_fibre_networks.json (the caller's inputs)
%
% OUTPUTS (OutputDir, prefix <sub>_fn_):
%   suvr.csv  bundles.csv  bundle_measures.csv  profiles.csv  cwa.csv
%   connectome_<atlas>.mat, edges_<atlas>.csv, endpoints.mat, and <sub>_fibre_networks.json
%   R : struct (.Files, .Summary, .Qc)
%
% SEE ALSO: pet_gtm, pet_suvr, pet_scanner_fwhm, fibers_helper, cs_convert, bst_nearest

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

t0 = tic;
if (nargin < 3) || isempty(Opts), Opts = struct(); end
Def = struct('Tracers', {{'18FNAV4694','18Fflortaucipir'}}, 'AmyloidTracer', '18FNAV4694', ...
    'TauTracer', '18Fflortaucipir', 'PsfFwhm', [], 'RefLabels', [8 47], ...
    'Atlases', {{'Desikan-Killiany','Destrieux'}}, 'EndRadius', 3, 'ProfilePoints', 20, ...
    'ChunkSize', 20000, 'ScalarPattern', '^(DTI[ _-]?)?(FA|MD|RD|AD|AFD)\b', 'Provenance', struct());
fn = fieldnames(Def);
for i = 1:numel(fn), if ~isfield(Opts, fn{i}) || isempty(Opts.(fn{i})), Opts.(fn{i}) = Def.(fn{i}); end, end
if ~exist(OutputDir, 'dir'), mkdir(OutputDir); end
pfx = fullfile(OutputDir, [SubjectName '_fn_']);
R = struct('Files', {{}}, 'Summary', '', 'Qc', struct());
qc = struct();  notes = {};
logmsg = @(varargin) fprintf(['BST> fibre_networks: ' varargin{1} '\n'], varargin{2:end});

%% ===== SUBJECT =====
[sSubject, iSubject] = bst_get('Subject', SubjectName);
if isempty(sSubject), error('bst_fibre_networks:subject', 'Subject %s not in the current protocol.', SubjectName); end
anatCm = {sSubject.Anatomy.Comment};
sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
cubeSize = size(sMri.Cube(:,:,:,1));
vox = sMri.Voxsize(1:3);
iAseg = find(strcmpi(anatCm, 'ASEG'), 1);
if isempty(iAseg), error('bst_fibre_networks:aseg', 'No ASEG volume atlas.'); end
sAseg = in_mri_bst(sSubject.Anatomy(iAseg).FileName);
aseg = int32(sAseg.Cube(:,:,:,1));
refMask = ismember(aseg, Opts.RefLabels);

% Volume atlases (labels as int32 cubes on the T1 grid) and their node sets
SUBC = [10 11 12 13 17 18 26 28 16 8 49 50 51 52 53 54 58 60 47];   % aseg grey-matter nodes
At = struct('name', {}, 'tag', {}, 'L', {}, 'ids', {}, 'labels', {}, 'hemi', {}, 'cortical', {}, 'file', {}, ...
    'allIds', {}, 'allLabels', {});
for ia = 1:numel(Opts.Atlases)
    k = find(strcmpi(anatCm, Opts.Atlases{ia}), 1);
    if isempty(k), notes{end+1} = sprintf('volume atlas %s absent', Opts.Atlases{ia}); continue; end %#ok<AGROW>
    sA = in_mri_bst(sSubject.Anatomy(k).FileName);
    L = int32(sA.Cube(:,:,:,1));
    if ~isequal(size(L), cubeSize), error('bst_fibre_networks:grid', '%s grid differs from the T1.', Opts.Atlases{ia}); end
    present = unique(L(L > 0));
    lab = local_labelmap(sA.Labels);
    isCtx = present >= 1000;
    isNode = isCtx | ismember(present, SUBC);
    ids = double(present(isNode));
    At(end+1) = struct('name', Opts.Atlases{ia}, 'tag', local_tag(Opts.Atlases{ia}), 'L', L, ...
        'ids', ids, 'labels', {local_names(lab, ids)}, 'hemi', {local_hemi(ids, local_names(lab, ids))}, ...
        'cortical', ids >= 1000, 'file', sSubject.Anatomy(k).FileName, ...
        'allIds', double(present), 'allLabels', {local_names(lab, double(present))}); %#ok<AGROW>
end
if isempty(At), error('bst_fibre_networks:atlas', 'None of the volume atlases %s.', strjoin(Opts.Atlases, ', ')); end
iDK = find(strcmp({At.tag}, 'dk'), 1);
logmsg('%s: T1 %s, atlases %s', SubjectName, mat2str(cubeSize), strjoin({At.name}, ', '));

%% ===== PET VOLUMES, PER TRACER AND METHOD =====
sample = struct('name', {}, 'tracer', {}, 'method', {}, 'cube', {});   % volumes sampled along tracts
suvrRows = {};   % atlas, id, label, hemi, tracer, method, nvox, suvr
pet = struct();
fwhm = Opts.PsfFwhm; fwhmSrc = 'Opts.PsfFwhm';
for it = 1:numel(Opts.Tracers)
    trc = Opts.Tracers{it};
    iMean = find(strcmp(anatCm, ['PET ' trc '_mean']), 1);
    iSuvr = find(strcmp(anatCm, ['PET ' trc '_suvr']), 1);
    if isempty(iMean) || isempty(iSuvr)
        notes{end+1} = sprintf('tracer %s: no "PET %s_mean"/"_suvr" volume', trc, trc); %#ok<AGROW>
        continue;
    end
    sMean = in_mri_bst(sSubject.Anatomy(iMean).FileName);
    sV = in_mri_bst(sSubject.Anatomy(iSuvr).FileName);
    if ~isequal(size(sMean.Cube(:,:,:,1)), cubeSize), error('bst_fibre_networks:grid', 'PET %s grid differs from the T1.', trc); end
    if isempty(fwhm)
        PETmeta = []; if isfield(sMean, 'PET'), PETmeta = sMean.PET; end
        [fwhm, fwhmSrc] = pet_scanner_fwhm(PETmeta);
    end
    m = double(sMean.Cube(:,:,:,1));
    ref = mean(m(refMask));
    P = struct('tracer', trc, 'meanFile', sSubject.Anatomy(iMean).FileName, ...
        'suvrFile', sSubject.Anatomy(iSuvr).FileName, 'ref_none', ref);
    vols = struct('method', {'none', 'vlpp'}, 'cube', {single(m / ref), single(sV.Cube(:,:,:,1))});
    % --- GTM (Desikan-Killiany partition) and RBV ---
    P.gtm_error = '';
    try
        [gtmFile, errMsg, regTable] = pet_gtm(P.meanFile, fwhm, struct('AtlasComment', At(iDK).name));
        if ~isempty(errMsg), error(errMsg); end
        sG = in_mri_bst(gtmFile);
        G = double(sG.Cube(:,:,:,1));
        isRef = ismember([regTable.id], Opts.RefLabels);
        refG = sum([regTable(isRef).corrected] .* [regTable(isRef).nvox]) / sum([regTable(isRef).nvox]);
        sGs = local_gauss3(G, fwhm, vox);
        rbv = zeros(size(m));
        ok = sGs > 1e-6 * max(abs(sGs(:)));
        rbv(ok) = m(ok) .* G(ok) ./ sGs(ok);
        refR = mean(rbv(refMask));
        vols(end+1) = struct('method', 'rbv', 'cube', single(rbv / refR)); %#ok<AGROW>
        P.ref_gtm = refG; P.ref_rbv = refR; P.gtmFile = gtmFile;
        P.gtm_regions = numel(regTable);
        % GTM regional SUVR rows (its own partition = DK labels with >= minVox voxels)
        for r = 1:numel(regTable)
            id = regTable(r).id; if id <= 0, continue; end
            j = find(At(iDK).allIds == id, 1); if isempty(j), continue; end
            nm = At(iDK).allLabels{j};
            suvrRows(end+1, :) = {At(iDK).tag, id, nm, local_hemi(id, {nm}), trc, 'gtm', regTable(r).nvox, regTable(r).corrected / refG}; %#ok<AGROW>
        end
        clear G sGs rbv sG;
    catch ME
        P.gtm_error = ME.message;
        notes{end+1} = sprintf('tracer %s: GTM/RBV failed: %s', trc, ME.message); %#ok<AGROW>
    end
    % --- volumetric regional means of every method on every atlas ---
    for ia = 1:numel(At)
        lin = find(At(ia).L > 0);
        [~, li] = ismember(At(ia).L(lin), At(ia).allIds);
        nv = accumarray(li, 1, [numel(At(ia).allIds) 1]);
        for iv = 1:numel(vols)
            s = accumarray(li, double(vols(iv).cube(lin)), [numel(At(ia).allIds) 1]) ./ max(nv, 1);
            for r = 1:numel(At(ia).allIds)
                suvrRows(end+1, :) = {At(ia).tag, At(ia).allIds(r), At(ia).allLabels{r}, ...
                    local_hemi(At(ia).allIds(r), At(ia).allLabels(r)), trc, vols(iv).method, nv(r), s(r)}; %#ok<AGROW>
            end
        end
    end
    % --- the protocol's cortical SUVR map over the scouts ---
    suvrRows = [suvrRows; local_surface_rows(sSubject, trc, At)]; %#ok<AGROW>
    for iv = 1:numel(vols)
        sample(end+1) = struct('name', ['PET_' trc '_' vols(iv).method], 'tracer', trc, ...
            'method', vols(iv).method, 'cube', vols(iv).cube); %#ok<AGROW>
    end
    pet.(local_field(trc)) = P;
    clear m sMean sV vols;
end
% Microstructure volumes, if the subject holds any
iSc = find(~cellfun(@isempty, regexpi(anatCm, Opts.ScalarPattern, 'once')));
for k = iSc(:)'
    sS = in_mri_bst(sSubject.Anatomy(k).FileName);
    if ~isequal(size(sS.Cube(:,:,:,1)), cubeSize), notes{end+1} = sprintf('scalar %s off-grid, skipped', anatCm{k}); continue; end %#ok<AGROW>
    sample(end+1) = struct('name', matlab.lang.makeValidName(anatCm{k}), 'tracer', '', 'method', anatCm{k}, ...
        'cube', single(sS.Cube(:,:,:,1))); %#ok<AGROW>
end
qc.microstructure = {sample(cellfun(@isempty, {sample.tracer})).method};
if isempty(qc.microstructure)
    notes{end+1} = ['no microstructure volume in the protocol (pattern ' Opts.ScalarPattern ...
        '): FA/MD/RD/AD/AFD profiles not computed']; %#ok<AGROW>
end
logmsg('PET: %d volumes sampled along tracts; PSF %.2f mm (%s)', numel(sample), fwhm, fwhmSrc);
local_writecsv([pfx 'suvr.csv'], {'atlas','label_id','label','hemi','tracer','method','nvox','suvr'}, suvrRows);
R.Files{end+1} = [pfx 'suvr.csv'];

%% ===== FIBERS =====
iFib = [];
if isfield(sSubject, 'iFibers'), iFib = sSubject.iFibers; end
if isempty(iFib), iFib = find(strcmpi({sSubject.Surface.SurfaceType}, 'Fibers'), 1); end
if isempty(iFib), error('bst_fibre_networks:fibers', 'Subject %s has no fibers.', SubjectName); end
fibFile = sSubject.Surface(iFib).FileName;
F = load(file_fullpath(fibFile), 'Points', 'Scouts', 'Comment');
[nFib, nPts, ~] = size(F.Points);
logmsg('fibers %s: %d streamlines x %d points', fibFile, nFib, nPts);
% lengths (mm) and endpoints (SCS, m)
len = zeros(nFib, 1, 'single');
for c = 1:Opts.ChunkSize:nFib
    ii = c:min(nFib, c + Opts.ChunkSize - 1);
    d = diff(F.Points(ii,:,:), 1, 2);
    len(ii) = single(sum(sqrt(sum(d.^2, 3)), 2) * 1000);
end
E = {reshape(F.Points(:,1,:), nFib, 3), reshape(F.Points(:,end,:), nFib, 3)};
Evox = {cs_convert(sMri, 'scs', 'voxel', E{1}), cs_convert(sMri, 'scs', 'voxel', E{2})};
nodeEnd = struct();
for ia = 1:numel(At)
    N = At(ia).L; N(~ismember(N, int32(At(ia).ids))) = 0;
    exact = cell(1,2); for e = 1:2, exact{e} = local_lookup(N, Evox{e}); end
    N = local_dilate(N, round(Opts.EndRadius / min(vox)));
    a = zeros(nFib, 2, 'int32');
    for e = 1:2, a(:, e) = local_lookup(N, Evox{e}); end
    nodeEnd.(At(ia).tag) = a;
    qc.(['ends_in_gm_' At(ia).tag]) = mean([exact{1}; exact{2}] > 0);
    qc.(['ends_assigned_' At(ia).tag]) = mean(a(:) > 0);
    qc.(['streamlines_assigned_' At(ia).tag]) = mean(all(a > 0, 2));
end
% nearest cortex vertex of every end (any cortical tiling rolls up from it)
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);
[vEnd1, dEnd1] = bst_nearest(sCortex.Vertices, E{1}, 1, 0);
[vEnd2, dEnd2] = bst_nearest(sCortex.Vertices, E{2}, 1, 0);
endpoints = struct('Subject', SubjectName, 'FibersFile', fibFile, 'nFibers', nFib, 'LengthMm', len, ...
    'CortexFile', sSubject.Surface(sSubject.iCortex).FileName, 'Vertex', int32([vEnd1 vEnd2]), ...
    'VertexDistMm', single([dEnd1 dEnd2] * 1000), 'Node', nodeEnd, 'EndRadiusMm', Opts.EndRadius); %#ok<NASGU>
save([pfx 'endpoints.mat'], '-struct', 'endpoints', '-v7');
R.Files{end+1} = [pfx 'endpoints.mat'];
qc.median_end_to_cortex_mm = median(double([dEnd1; dEnd2])) * 1000;
clear E Evox vEnd1 vEnd2 dEnd1 dEnd2;

%% ===== CONNECTOMES =====
Conn = struct();
for ia = 1:numel(At)
    a = nodeEnd.(At(ia).tag); ok = all(a > 0, 2);
    [~, i1] = ismember(a(ok,1), At(ia).ids); [~, i2] = ismember(a(ok,2), At(ia).ids);
    n = numel(At(ia).ids);
    C = accumarray([i1 i2], 1, [n n]); C = C + C' - diag(diag(C));
    Ls = accumarray([i1 i2], double(len(ok)), [n n]); Ls = Ls + Ls' - diag(diag(Ls));
    ML = Ls ./ max(C, 1); ML(C == 0) = NaN;
    cen = local_centroids(At(ia).L, At(ia).ids, vox);
    S = struct('Atlas', At(ia).name, 'Subject', SubjectName, 'Ids', At(ia).ids, 'Labels', {At(ia).labels}, ...
        'Hemi', {At(ia).hemi}, 'Cortical', At(ia).cortical, 'CentroidMri_mm', cen, 'Count', C, 'MeanLengthMm', ML, ...
        'Measure', 'streamline count, volume-atlas endpoints', 'EndRadiusMm', Opts.EndRadius, 'nFibers', nFib);
    % the protocol's own AssignToScouts connectome (nearest scout seed), carried for comparison
    S.BstAssign = local_bstassign(F, sCortex, At(ia).name, SubjectName);
    save([pfx 'connectome_' At(ia).tag '.mat'], '-struct', 'S', '-v7');
    [ei, ej] = find(triu(C));
    rows = cell(numel(ei), 6);
    for k = 1:numel(ei)
        rows(k,:) = {At(ia).ids(ei(k)), At(ia).ids(ej(k)), At(ia).labels{ei(k)}, At(ia).labels{ej(k)}, C(ei(k),ej(k)), ML(ei(k),ej(k))};
    end
    local_writecsv([pfx 'edges_' At(ia).tag '.csv'], {'id_a','id_b','label_a','label_b','count','mean_length_mm'}, rows);
    R.Files = [R.Files, {[pfx 'connectome_' At(ia).tag '.mat'], [pfx 'edges_' At(ia).tag '.csv']}];
    Conn.(At(ia).tag) = S;
end

%% ===== TAU EPICENTRE (from tau alone) =====
epi = struct();
tauRows = suvrRows(strcmp(suvrRows(:,1), 'dk') & strcmp(suvrRows(:,5), Opts.TauTracer), :);
epiMethod = '';
for meth = {'gtm', 'rbv', 'none'}
    if any(strcmp(tauRows(:,6), meth{1})), epiMethod = meth{1}; break; end
end
for h = 'LR'
    cand = local_ids({'ctx:entorhinal','ctx:parahippocampal','ctx:inferiortemporal','ctx:fusiform','aseg:Amygdala'}, h);
    if isempty(epiMethod), epi.(h) = struct('id', NaN, 'label', '', 'suvr', NaN); continue; end
    tr = tauRows(strcmp(tauRows(:,6), epiMethod) & ismember(cell2mat(tauRows(:,2)), cand), :);
    [~, k] = max(cell2mat(tr(:,8)));
    epi.(h) = struct('id', tr{k,2}, 'label', tr{k,3}, 'suvr', tr{k,8});
end
epi.method = epiMethod; epi.candidates = 'entorhinal, parahippocampal, inferiortemporal, fusiform, Amygdala (per hemisphere)';
epi.tracer = Opts.TauTracer;

%% ===== A-PRIORI BUNDLES (Desikan-Killiany nodes) =====
B = local_bundles();
a = nodeEnd.dk;
memb = false(nFib, 0); flp = false(nFib, 0); bInfo = {};
for ib = 1:numel(B)
    for h = 'LR'
        if strcmp(B(ib).A{1}, '@EPI')
            if ~isfinite(epi.(h).id), continue; end
            A = epi.(h).id;
        else
            A = local_ids(B(ib).A, h);
        end
        if strcmp(B(ib).B{1}, '*')
            inA1 = ismember(a(:,1), A); inA2 = ismember(a(:,2), A);
            mk = (inA1 & a(:,2) > 0 & ~inA2) | (inA2 & a(:,1) > 0 & ~inA1);
            fl = ~inA1 & inA2;
        else
            Bn = local_ids(B(ib).B, h);
            ab = ismember(a(:,1), A) & ismember(a(:,2), Bn);
            ba = ismember(a(:,1), Bn) & ismember(a(:,2), A);
            mk = ab | ba; fl = ba & ~ab;
        end
        memb(:, end+1) = mk; flp(:, end+1) = fl; %#ok<AGROW>
        bInfo(end+1, :) = {B(ib).id, h, B(ib).class, B(ib).label, strjoin(B(ib).A, ';'), strjoin(B(ib).B, ';')}; %#ok<AGROW>
    end
end
nB = size(memb, 2);
nAssigned = sum(all(a > 0, 2));
bRows = cell(nB, 10);
for k = 1:nB
    n = sum(memb(:,k)); L = double(len(memb(:,k)));
    bRows(k,:) = [bInfo(k,:), {n, n / max(nAssigned,1), mean(L), std(L)}];
end
local_writecsv([pfx 'bundles.csv'], {'bundle','hemi','class','label','A','B','n_streamlines','frac_assigned','length_mean_mm','length_sd_mm'}, bRows);
R.Files{end+1} = [pfx 'bundles.csv'];

%% ===== ALONG-TRACT SAMPLING (one pass over the bundle streamlines) =====
K = Opts.ProfilePoints; nV = numel(sample);
edges = round(linspace(0, nPts, K + 1));
binOf = zeros(1, nPts); for k = 1:K, binOf(edges(k)+1:edges(k+1)) = k; end
Bin = sparse(1:nPts, binOf, 1, nPts, K); Bin = Bin ./ sum(Bin, 1);       % point -> bin average
acc = struct('s', zeros(nB, nV), 'ss', zeros(nB, nV), 'n', zeros(nB, nV), ...
    'ps', zeros(nB, nV, K), 'pss', zeros(nB, nV, K), 'pn', zeros(nB, nV, K));
U = find(any(memb, 2));
logmsg('%d bundles x hemis, %d streamlines in any bundle, %d volumes', nB, numel(U), nV);
for c = 1:Opts.ChunkSize:numel(U)
    ii = U(c:min(numel(U), c + Opts.ChunkSize - 1)); nc = numel(ii);
    Pv = cs_convert(sMri, 'scs', 'voxel', reshape(F.Points(ii,:,:), nc * nPts, 3));
    for v = 1:nV
        X = double(reshape(interp3(sample(v).cube, Pv(:,2), Pv(:,1), Pv(:,3), 'linear', NaN), nc, nPts));
        for k = 1:nB
            sel = memb(ii, k); if ~any(sel), continue; end
            Xk = X(sel, :); f = flp(ii(sel), k); Xk(f, :) = Xk(f, end:-1:1);
            good = all(isfinite(Xk), 2); Xk = Xk(good, :);
            mu = mean(Xk, 2);
            acc.s(k,v) = acc.s(k,v) + sum(mu); acc.ss(k,v) = acc.ss(k,v) + sum(mu.^2); acc.n(k,v) = acc.n(k,v) + numel(mu);
            Pk = Xk * Bin;
            acc.ps(k,v,:) = acc.ps(k,v,:) + reshape(sum(Pk, 1), 1, 1, K);
            acc.pss(k,v,:) = acc.pss(k,v,:) + reshape(sum(Pk.^2, 1), 1, 1, K);
            acc.pn(k,v,:) = acc.pn(k,v,:) + size(Pk, 1);
        end
    end
end
mRows = cell(nB * nV, 9); pRows = cell(nB * nV * K, 10); im = 0; ip = 0;
for k = 1:nB
    for v = 1:nV
        n = acc.n(k,v); mu = acc.s(k,v) / n; sd = sqrt(max(acc.ss(k,v) / n - mu^2, 0) * n / max(n-1, 1));
        im = im + 1; mRows(im,:) = {bInfo{k,1}, bInfo{k,2}, sample(v).name, sample(v).tracer, sample(v).method, n, mu, sd, 'per-streamline mean over its points'};
        for q = 1:K
            n = acc.pn(k,v,q); mu = acc.ps(k,v,q) / n; sd = sqrt(max(acc.pss(k,v,q) / n - mu^2, 0) * n / max(n-1, 1));
            ip = ip + 1; pRows(ip,:) = {bInfo{k,1}, bInfo{k,2}, sample(v).name, sample(v).tracer, sample(v).method, q, (q - 0.5) / K, n, mu, sd};
        end
    end
end
local_writecsv([pfx 'bundle_measures.csv'], {'bundle','hemi','measure','tracer','method','n','mean','sd','definition'}, mRows(1:im,:));
local_writecsv([pfx 'profiles.csv'], {'bundle','hemi','measure','tracer','method','point','position_A_to_B','n','mean','sd'}, pRows(1:ip,:));
R.Files = [R.Files, {[pfx 'bundle_measures.csv'], [pfx 'profiles.csv']}];
clear F sample;

%% ===== CONNECTOME-WEIGHTED PET =====
cRows = {};
for ia = 1:numel(At)
    S = Conn.(At(ia).tag); ctx = find(S.Cortical(:)');
    W = S.Count(ctx, ctx); W(1:numel(ctx)+1:end) = 0;
    cen = S.CentroidMri_mm(ctx, :);
    D = sqrt(max(sum(cen.^2, 2) + sum(cen.^2, 2)' - 2 * (cen * cen'), 0));
    Wd = 1 ./ D; Wd(1:numel(ctx)+1:end) = 0;
    rows = suvrRows(strcmp(suvrRows(:,1), At(ia).tag), :);
    for it = 1:numel(Opts.Tracers)
        meths = setdiff(unique(rows(strcmp(rows(:,5), Opts.Tracers{it}), 6)), {'vlpp_surface'});
        for meth = meths(:)'
            rr = rows(strcmp(rows(:,5), Opts.Tracers{it}) & strcmp(rows(:,6), meth{1}), :);
            Av = nan(numel(ctx), 1);
            [isIn, loc] = ismember(S.Ids(ctx), cell2mat(rr(:,2)));
            Av(isIn) = cell2mat(rr(loc(isIn), 8));
            for s = 1:numel(ctx)
                okj = isfinite(Av) & ((1:numel(ctx))' ~= s);
                w = W(s, okj)'; ws = sum(w);
                cwa = sum(w .* Av(okj)) / ws; if ws == 0, cwa = NaN; end
                idw = sum(Wd(s, okj)' .* Av(okj)) / sum(Wd(s, okj));
                cRows(end+1, :) = {At(ia).tag, S.Ids(ctx(s)), S.Labels{ctx(s)}, S.Hemi{ctx(s)}, Opts.Tracers{it}, meth{1}, ...
                    ws, cwa, Av(s), mean(Av(okj)), idw}; %#ok<AGROW>
            end
        end
    end
end
local_writecsv([pfx 'cwa.csv'], {'atlas','seed_id','seed_label','hemi','tracer','method','strength','cwa','seed_local','cortical_mean','idw_mean'}, cRows);
R.Files{end+1} = [pfx 'cwa.csv'];

%% ===== PROVENANCE =====
qc.n_fibers = nFib; qc.n_points = nPts; qc.length_median_mm = median(double(len));
qc.bundles_empty = sum(cell2mat(bRows(:,7)) == 0);
qc.seconds = toc(t0);
J = struct('subject', SubjectName, 'function', 'bst_fibre_networks', 'brainstorm', bst_get('Version'), ...
    'date', datestr(now, 'yyyy-mm-ddTHH:MM:SS'), 'opts', rmfield(Opts, 'Provenance'), ...
    'inputs', struct('t1', sSubject.Anatomy(sSubject.iAnatomy).FileName, 'aseg', sSubject.Anatomy(iAseg).FileName, ...
        'atlases', {{At.file}}, 'fibers', fibFile, 'fibers_comment', '', 'cortex', sSubject.Surface(sSubject.iCortex).FileName), ...
    'psf_fwhm_mm', fwhm, 'psf_source', fwhmSrc, 'pet', pet, 'epicentre', epi, 'qc', qc, ...
    'bundles', {local_bundles()}, 'notes', {notes}, 'files', {cellfun(@local_base, R.Files, 'UniformOutput', 0)}, ...
    'not_computed', {{'association tests of any kind (rule M4)', 'rewired / spun nulls (test plan)', ...
        'SIFT2 weights: Brainstorm fibers carry none, so edges are streamline counts', ...
        'Mueller-Gartner (needs SPM12 + PETPVE12; GTM and RBV are native)'}});
fns = fieldnames(Opts.Provenance);
for i = 1:numel(fns), J.(fns{i}) = Opts.Provenance.(fns{i}); end
jf = fullfile(OutputDir, [SubjectName '_fibre_networks.json']);
fid = fopen(jf, 'w'); fprintf(fid, '%s', jsonencode(J, 'PrettyPrint', true)); fclose(fid);
R.Files{end+1} = jf;
R.Qc = qc;
R.Summary = sprintf('fibers=%d;assigned_dk=%.3f;bundles=%d;empty=%d;tracers=%d;gtm=%d;s=%.0f', nFib, ...
    qc.streamlines_assigned_dk, nB, qc.bundles_empty, numel(fieldnames(pet)), ...
    sum(structfun(@(p) isempty(p.gtm_error), pet)), qc.seconds);
logmsg('done in %.0f s: %s', qc.seconds, R.Summary);
end


%% ======================================================================================
%  LOCAL FUNCTIONS
%  ======================================================================================
function B = local_bundles()
% A-priori node-pair sets (Desikan-Killiany cortex "ctx:" + aseg "aseg:"), per hemisphere.
% Rationales: projects/preventad-multimodal/staging/evidence-amyloid-tau-networks (D1-D6).
B = struct('id', {}, 'class', {}, 'A', {}, 'B', {}, 'label', {}, 'rationale', {});
add = @(B, id, cl, A, Bn, lb, ra) [B, struct('id', id, 'class', cl, 'A', {A}, 'B', {Bn}, 'label', lb, 'rationale', ra)];
B = add(B, 'B01', 'epicentre', {'ctx:entorhinal','ctx:parahippocampal'}, {'ctx:isthmuscingulate','ctx:posteriorcingulate'}, 'hippocampal (parahippocampal) cingulum', 'EC/MTL -> PCC; Jacobs2018, PichetBinette2021');
B = add(B, 'B02', 'epicentre', {'ctx:isthmuscingulate','ctx:posteriorcingulate'}, {'ctx:precuneus','ctx:caudalanteriorcingulate','ctx:superiorfrontal'}, 'posterior/dorsal cingulum', 'downstream of the EC->PCC relay; PichetBinette2021');
B = add(B, 'B03', 'epicentre', {'ctx:entorhinal','ctx:temporalpole','aseg:Amygdala'}, {'ctx:lateralorbitofrontal','ctx:medialorbitofrontal','ctx:parsorbitalis'}, 'uncinate fasciculus', 'PichetBinette2021 uncinate');
B = add(B, 'B04', 'epicentre', {'ctx:entorhinal','ctx:inferiortemporal'}, {'ctx:temporalpole','ctx:middletemporal','ctx:fusiform'}, 'anterior temporal white matter', 'Strain2018 tau -> anterior temporal MD');
B = add(B, 'B05', 'epicentre', {'ctx:inferiortemporal'}, {'ctx:lateraloccipital','ctx:inferiorparietal','ctx:precuneus'}, 'inferior-temporal epicentre -> posterior association hubs', 'Lee2022; Sepulcre2016');
B = add(B, 'B06', 'epicentre', {'aseg:Hippocampus'}, {'aseg:VentralDC'}, 'fornix PROXY (hippocampus <-> ventral diencephalon)', 'aseg has no fornix/mammillary label');
B = add(B, 'B07', 'epicentre', {'ctx:entorhinal'}, {'aseg:Hippocampus'}, 'perforant-path PROXY (entorhinal <-> hippocampus)', 'not resolvable at DWI resolution (Yassa2010)');
B = add(B, 'C01', 'control', {'ctx:precentral'}, {'ctx:postcentral'}, 'sensorimotor U-fibres', 'primary cortex, last Braak stage');
B = add(B, 'C02', 'control', {'ctx:precentral','ctx:paracentral'}, {'aseg:Brain-Stem'}, 'corticospinal proxy', 'projection fibres, no MTL');
B = add(B, 'C03', 'control', {'aseg:Thalamus'}, {'ctx:pericalcarine','ctx:cuneus','ctx:lingual'}, 'optic radiation proxy', 'primary visual');
B = add(B, 'C04', 'control', {'ctx:supramarginal','ctx:inferiorparietal'}, {'ctx:caudalmiddlefrontal','ctx:parsopercularis'}, 'superior longitudinal fasciculus III proxy', 'fronto-parietal association, not EC-originating');
B = add(B, 'C05', 'control', {'ctx:lateraloccipital'}, {'ctx:superiorparietal'}, 'posterior visual association', 'occipito-parietal');
B = add(B, 'E01', 'epicentre-any', {'ctx:entorhinal'}, {'*'}, 'all streamlines leaving entorhinal', 'fixed EC epicentre');
B = add(B, 'E02', 'epicentre-any', {'ctx:inferiortemporal'}, {'*'}, 'all streamlines leaving inferior temporal', 'fixed IT epicentre');
B = add(B, 'E03', 'epicentre-any', {'@EPI'}, {'*'}, 'all streamlines leaving the subject''s tau epicentre', 'max GTM flortaucipir SUVR among entorhinal, parahippocampal, inferiortemporal, fusiform, amygdala');
end

function ids = local_ids(tokens, h)
% FreeSurfer ids of "ctx:<dk name>" / "aseg:<structure>" tokens in hemisphere h ('L'|'R').
DK = {'bankssts','caudalanteriorcingulate','caudalmiddlefrontal','corpuscallosum','cuneus','entorhinal', ...
      'fusiform','inferiorparietal','inferiortemporal','isthmuscingulate','lateraloccipital', ...
      'lateralorbitofrontal','lingual','medialorbitofrontal','middletemporal','parahippocampal', ...
      'paracentral','parsopercularis','parsorbitalis','parstriangularis','pericalcarine','postcentral', ...
      'posteriorcingulate','precentral','precuneus','rostralanteriorcingulate','rostralmiddlefrontal', ...
      'superiorfrontal','superiorparietal','superiortemporal','supramarginal','frontalpole', ...
      'temporalpole','transversetemporal','insula'};
ASEG = struct('Thalamus', [10 49], 'Caudate', [11 50], 'Putamen', [12 51], 'Pallidum', [13 52], ...
      'Hippocampus', [17 53], 'Amygdala', [18 54], 'Accumbens', [26 58], 'VentralDC', [28 60], ...
      'Cerebellum', [8 47], 'BrainStem', [16 16]);
ids = [];
hi = 1 + (h == 'R');
for t = 1:numel(tokens)
    tok = strsplit(tokens{t}, ':');
    switch tok{1}
        case 'ctx'
            k = find(strcmp(DK, tok{2}));
            if isempty(k), error('bst_fibre_networks:bundle', 'Unknown DK region %s', tok{2}); end
            ids(end+1) = 1000 * hi + k; %#ok<AGROW>
        case 'aseg'
            f = strrep(tok{2}, '-', '');
            if ~isfield(ASEG, f), error('bst_fibre_networks:bundle', 'Unknown aseg structure %s', tok{2}); end
            ids(end+1) = ASEG.(f)(hi); %#ok<AGROW>
    end
end
end

function lab = local_labelmap(Labels)
lab = containers.Map('KeyType', 'double', 'ValueType', 'char');
if isempty(Labels), return; end
for k = 1:size(Labels, 1), lab(double(Labels{k,1})) = Labels{k,2}; end
end

function nm = local_names(lab, ids)
nm = cell(numel(ids), 1);
for k = 1:numel(ids)
    if isKey(lab, ids(k)), nm{k} = lab(ids(k)); else, nm{k} = sprintf('label_%d', ids(k)); end
end
end

function h = local_hemi(ids, names)
h = cell(numel(ids), 1);
for k = 1:numel(ids)
    id = ids(k);
    if (id >= 1000 && id < 2000) || (id >= 11100 && id < 12000), h{k} = 'L';
    elseif (id >= 2000 && id < 3000) || (id >= 12100 && id < 13000), h{k} = 'R';
    elseif ~isempty(regexp(names{k}, '(\s|-|_)L$|^Left', 'once')), h{k} = 'L';
    elseif ~isempty(regexp(names{k}, '(\s|-|_)R$|^Right', 'once')), h{k} = 'R';
    else, h{k} = 'M';
    end
end
if numel(h) == 1, h = h{1}; end
end

function t = local_tag(name)
switch lower(name)
    case 'desikan-killiany', t = 'dk';
    case 'destrieux',        t = 'dx';
    otherwise,               t = lower(regexprep(name, '[^A-Za-z0-9]', ''));
end
end

function f = local_field(trc)
f = matlab.lang.makeValidName(['trc_' trc]);
end

function s = local_base(f)
[~, n, e] = fileparts(f); s = [n e];
end

function v = local_lookup(N, P)
% Label at the nearest voxel of points P (Brainstorm 'voxel' coordinates, 1-based centres).
sz = size(N);
ijk = round(P);
ok = all(ijk >= 1, 2) & ijk(:,1) <= sz(1) & ijk(:,2) <= sz(2) & ijk(:,3) <= sz(3);
v = zeros(size(P,1), 1, 'int32');
v(ok) = N(sub2ind(sz, ijk(ok,1), ijk(ok,2), ijk(ok,3)));
end

function N = local_dilate(N, R)
% Grow node labels into unlabelled voxels by up to R voxels (6-neighbourhood, one shell per step).
for r = 1:R
    M = N;
    for d = 1:3
        for s = [-1 1]
            S = circshift(N, s, d);
            fill = (M == 0) & (S ~= 0);
            M(fill) = S(fill);
        end
    end
    N = M;
end
end

function cen = local_centroids(L, ids, vox)
lin = find(L > 0);
[isN, li] = ismember(L(lin), ids);
[x, y, z] = ind2sub(size(L), lin(isN));
n = accumarray(li(isN), 1, [numel(ids) 1]);
cen = [accumarray(li(isN), x, [numel(ids) 1]), accumarray(li(isN), y, [numel(ids) 1]), ...
       accumarray(li(isN), z, [numel(ids) 1])] ./ max(n, 1) .* vox(:)';
end

function vol = local_gauss3(vol, fwhm_mm, voxsize_mm)
% Separable Gaussian PSF, identical to pet_gtm's (the RBV denominator must match the GTM PSF).
sig = fwhm_mm / 2.35482;
for d = 1:3
    sv = sig(min(d, numel(sig))) / voxsize_mm(d); r = max(1, ceil(3 * sv)); x = -r:r;
    k = exp(-(x.^2) / (2 * sv^2)); k = k / sum(k);
    sh = ones(1, 3); sh(d) = numel(k);
    vol = convn(vol, reshape(k, sh), 'same');
end
end

function rows = local_surface_rows(sSubject, trc, At)
% The protocol's cortical SUVR map ("PET <trc>_suvr", condition PET) over each atlas's scouts.
rows = cell(0, 8);
sStudies = bst_get('StudyWithSubject', sSubject.FileName);
rf = ''; 
for j = 1:numel(sStudies)
    for r = 1:numel(sStudies(j).Result)
        if strcmp(sStudies(j).Result(r).Comment, ['PET ' trc '_suvr']), rf = sStudies(j).Result(r).FileName; end
    end
end
if isempty(rf), return; end
Res = in_bst_results(rf, 1, 'ImageGridAmp', 'SurfaceFile');
x = double(Res.ImageGridAmp(:, 1));
sSurf = in_tess_bst(Res.SurfaceFile);
for ia = 1:numel(At)
    k = find(strcmpi({sSurf.Atlas.Name}, At(ia).name), 1);
    if isempty(k), continue; end
    sc = sSurf.Atlas(k).Scouts;
    for s = 1:numel(sc)
        rows(end+1, :) = {At(ia).tag, NaN, sc(s).Label, sc(s).Region(1), trc, 'vlpp_surface', numel(sc(s).Vertices), mean(x(sc(s).Vertices))}; %#ok<AGROW>
    end
end
end

function S = local_bstassign(F, sCortex, atlasName, SubjectName)
% The fibers' own AssignToScouts result for this atlas (import_fibers_subject), as a count matrix.
S = struct('Labels', {{}}, 'Count', [], 'Measure', 'streamline count, fibers_helper AssignToScouts (nearest scout seed)');
if ~isfield(F, 'Scouts') || isempty(F.Scouts), return; end
k = find(strcmp({F.Scouts.ConnectFile}, sprintf('%s_%s', SubjectName, atlasName)), 1);
iA = find(strcmpi({sCortex.Atlas.Name}, atlasName), 1);
if isempty(k) || isempty(iA), return; end
asg = F.Scouts(k).Assignment; n = numel(sCortex.Atlas(iA).Scouts);
ok = all(asg > 0, 2);
C = accumarray(double(asg(ok, :)), 1, [n n]); C = C + C' - diag(diag(C));
S.Labels = {sCortex.Atlas(iA).Scouts.Label}; S.Count = C;
end

function local_writecsv(f, hdr, rows)
fid = fopen(f, 'w');
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:size(rows, 1)
    c = cell(1, size(rows, 2));
    for j = 1:size(rows, 2)
        v = rows{i,j};
        if ischar(v), if any(v == ',' | v == '"'), v = ['"' strrep(v, '"', '""') '"']; end, c{j} = v;
        elseif isempty(v) || (isnumeric(v) && isnan(v)), c{j} = '';
        elseif isnumeric(v) && v == round(v) && abs(v) < 1e9, c{j} = sprintf('%d', v);
        else, c{j} = sprintf('%.6g', v);
        end
    end
    fprintf(fid, '%s\n', strjoin(c, ','));
end
fclose(fid);
end
