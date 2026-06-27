function Res = batch_pvc_3way(subjects, BidsPetDir, PsRoot, ResFile)
% BATCH_PVC_3WAY: Resumable batch of the 3-way PVC comparison over subjects x tracers.
%
% For each subject: import PET (idempotent) -> aggregate 4D base to a static mean ->
% run pet_pvc (auto-FWHM from metadata) -> compare_pvc_3way (ours MG vs PETSurfer MGX
% vs GTM). Per-ROI rows + per-pair metrics are appended to ResFile after EACH
% subject/tracer, so a re-run continues where it stopped (skips entries already present).
% Safe to call with subject subsets to keep each call short.
%
% USAGE:  Res = batch_pvc_3way({'sub-MTL0002','sub-MTL0005'})
%         Res = batch_pvc_3way(subjects, BidsPetDir, PsRoot, ResFile)
%
% Author: Diellor Basha, 2026

    here = bst_fileparts(mfilename('fullpath'));
    if (nargin < 2) || isempty(BidsPetDir), BidsPetDir = '/Volumes/SpikeData-2/workspace/library/datasets/preventad/pet'; end
    if (nargin < 3) || isempty(PsRoot),     PsRoot     = fullfile(BidsPetDir,'derivatives','petsurfer'); end
    if (nargin < 4) || isempty(ResFile),    ResFile    = fullfile(here,'pvc_3way_results.mat'); end
    trc = {'18FNAV4694','18Fflortaucipir'};

    % resume
    if exist(ResFile,'file')
        S = load(ResFile,'Res'); Res = S.Res;
    else
        Res = struct('subject',{},'tracer',{},'rows',{},'metrics',{});
    end
    done = arrayfun(@(e) [e.subject '/' e.tracer], Res, 'UniformOutput', 0);

    iProt = bst_get('Protocol','preventad');
    if iProt ~= bst_get('iProtocol'), gui_brainstorm('SetCurrentProtocol', iProt); end

    for is = 1:numel(subjects)
        subj = subjects{is};
        if all(ismember({[subj '/' trc{1}],[subj '/' trc{2}]}, done))
            fprintf('== %s: both tracers done, skip ==\n', subj); continue;
        end
        % import PET (idempotent: skips tracers whose 4D base exists)
        try
            preventad_pet_import(BidsPetDir, subj);
        catch ME
            fprintf('!! import failed %s: %s\n', subj, ME.message); continue;
        end
        for t = 1:2
            key = [subj '/' trc{t}];
            if any(strcmp(done, key)), fprintf('   skip %s (done)\n', key); continue; end
            [sS,~] = bst_get('Subject', subj);
            T1file   = local_find(sS, 'MRI T1');
            baseFile = local_find(sS, ['PET ' trc{t}]);
            if isempty(baseFile), fprintf('   no base %s\n', key); continue; end
            statC = ['PET ' trc{t} '_mean']; pvcC = [statC '_pvc'];
            % static
            statFile = local_find(sS, statC);
            if isempty(statFile), statFile = mri_aggregate(baseFile, 'mean'); end
            % pvc (auto-FWHM)
            [sS,~] = bst_get('Subject', subj); pvcFile = local_find(sS, pvcC);
            if isempty(pvcFile)
                [pvcFile, errPvc] = pet_pvc(statFile, T1file, []);
                if ~isempty(errPvc) || isempty(pvcFile), fprintf('!! pvc failed %s: %s\n', key, errPvc); continue; end
            end
            % comment of the pvc node (robust to _02 uniquification)
            [sS,~] = bst_get('Subject', subj);
            ip = find(strcmp({sS.Anatomy.FileName}, pvcFile), 1);
            pvcComment = sS.Anatomy(ip).Comment;
            % compare
            try
                R = compare_pvc_3way(subj, trc{t}, pvcComment, PsRoot, struct('DoFig',0,'Verbose',0));
            catch ME
                fprintf('!! compare failed %s: %s\n', key, ME.message); continue;
            end
            Res(end+1) = struct('subject',subj,'tracer',trc{t},'rows',R.rows,'metrics',R.metrics); %#ok<AGROW>
            done{end+1} = key; %#ok<AGROW>
            save(ResFile, 'Res');
            m = R.metrics;
            fprintf('** DONE %s : ours-MGX r_ctx=%.2f | MGX-GTM r_ctx=%.2f | ours-GTM r_ctx=%.2f  (%d/%d)\n', ...
                key, m(1).r_ctx, m(2).r_ctx, m(3).r_ctx, numel(Res), 2*numel(subjects));
        end
    end
    fprintf('\n== batch_pvc_3way: %d subject/tracer entries in %s ==\n', numel(Res), ResFile);
end

function f = local_find(sS, comment)
    i = find(strcmp({sS.Anatomy.Comment}, comment), 1);
    if isempty(i), f = ''; else, f = sS.Anatomy(i).FileName; end
end
