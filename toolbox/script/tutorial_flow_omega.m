function OutputFiles = tutorial_flow_omega(mode)
% TUTORIAL_FLOW_OMEGA  Systematic, band-resolved, group cortical-flow analysis on OMEGA resting MEG.
%
%   OutputFiles = tutorial_flow_omega()          % 'dirac' (eigenbasis) -- default
%   OutputFiles = tutorial_flow_omega('mn')      % vanilla vertex-basis wMNE
%
% Whole-recording flow characterization on the omega_tutorial_test protocol (sub-0002, sub-0003):
% for each flow quantity, fuse the linear flow operator into an unconstrained source KERNEL
% (process_source_flow -> stays a kernel), run a banded Welch PSD over the FULL 600 s
% (process_psd, streamed), normalize (relative), project to the group template, smooth, and
% average across subjects. Also computes the integrated flow FUNCTIONALS (energy/enstrophy/
% helicity) as whole-recording time series (streamed) and their band PSD.
%
% The source kernel is chosen by `mode`:
%   'dirac' : the relative-Dirac EIGENBASIS amplitude inverse (results_DiracFlowAmp_KERNEL,
%             written by omega_write_dirac_kernels) -- the canonical eigenbasis flow operators.
%   'mn'    : the vanilla vertex-basis unconstrained wMNE kernel (results_MN_*_KERNEL).
% Outputs are tagged with the mode so both runs coexist (FlowPSD_div_dirac vs FlowPSD_div_mn).
%
% Reuses the exact Brainstorm patterns from tutorial_omega (PSD -> norm -> project -> smooth ->
% average). Assumes the omega protocol is current and the chosen kernels exist.
%
% Author: Diellor Basha, 2026

    if nargin < 1 || isempty(mode), mode = 'dirac'; end
    switch lower(mode)
        case 'dirac', selTag = 'Dirac-flow'; fileKey = 'DiracFlowAmp'; label = 'dirac';
        case 'mn',    selTag = 'MN';         fileKey = 'KERNEL';       label = 'mn';
        otherwise,    error('Unknown mode "%s" (use ''dirac'' or ''mn'').', mode);
    end

    OutputFiles = struct();
    bst_report('Start');

    % ---- 0. select the unconstrained source links (canonical selector; kernel-over-raw is
    %         only resolvable this way -- a bare filename / GetInputStruct does NOT find it) ----
    sSrc = bst_process('CallProcess', 'process_select_files_results', [], [], ...
        'subjectname', 'All', 'tag', selTag);
    % keep ONLY the inverse KERNELs -- the comment tag also matches derived flow maps
    % (process_source_flow appends '| flow:div' to the source comment), so filter by filename.
    if ~isempty(sSrc)
        sSrc = sSrc(contains({sSrc.FileName}, fileKey) & contains({sSrc.FileName}, 'KERNEL'));
    end
    if numel(sSrc) < 2
        error(['Need the two unconstrained %s kernels. ' ...
               'If they exist on disk but are not found, db_reload_studies the raw studies.'], upper(mode));
    end
    subjOf = {sSrc.SubjectName};
    fprintf('tutorial_flow_omega[%s]: %d kernels (%s)\n', label, numel(sSrc), strjoin(subjOf,', '));

    bands = {'delta','2, 4','mean'; 'theta','4, 8','mean'; 'alpha','8, 13','mean'; ...
             'beta','13, 30','mean'; 'gamma','30, 45','mean'};
    quantities = {'div','curl','psi'};     % source/sink, vorticity, vortex stream function
    funcs      = {'energy','enstrophy','helicity'};

    % ---- 1. flow MAPS: fuse operator -> banded PSD -> normalize -> project -> smooth ----
    for iq = 1:numel(quantities)
        Q = quantities{iq};  projFiles = {};
        for i = 1:numel(sSrc)
            sMap = bst_process('CallProcess', 'process_source_flow', sSrc(i), [], 'quantity', Q);
            sPsd = bst_process('CallProcess', 'process_psd', sMap, [], ...
                'timewindow', [], 'win_length', 4, 'win_overlap', 50, 'clusters', {}, 'scoutfunc', 1, ...
                'edit', struct('Comment', ['FlowPSD_' Q '_' label], 'TimeBands', [], 'Freqs', {bands}, ...
                    'ClusterFuncTime', 'none', 'Measure', 'power', 'Output', 'all', 'SaveKernel', 0));
            sPsdN = bst_process('CallProcess', 'process_tf_norm', sPsd, [], 'normalize', 'relative', 'overwrite', 0);
            sProj = bst_process('CallProcess', 'process_project_sources', sPsdN, [], 'headmodeltype', 'surface');
            sProj = bst_process('CallProcess', 'process_ssmooth_surfstat', sProj, [], 'fwhm', 3, 'overwrite', 1);
            if isempty(sMap) || isempty(sPsd) || isempty(sProj)
                fprintf('  [%s] %s -> SKIPPED (a step returned empty)\n', Q, subjOf{i});
                continue;
            end
            projFiles{end+1} = sProj.FileName; %#ok<AGROW>
            fprintf('  [%s] %s -> band PSD -> projected\n', Q, subjOf{i});
        end
        % group average across subjects
        sAvg = bst_process('CallProcess', 'process_average', projFiles, [], ...
            'avgtype', 1, 'avg_func', 1, 'weighted', 0, 'matchrows', 1, 'iszerobad', 0);
        sAvg = bst_process('CallProcess', 'process_add_tag', sAvg, [], ...
            'tag', ['groupflow_' Q '_' label], 'output', 2);
        OutputFiles.(['group_' Q]) = sAvg.FileName;
        fprintf('  [%s] group average -> %s\n', Q, sAvg.FileName);
    end

    % ---- 2. flow FUNCTIONALS: streamed whole-recording series -> band PSD (per subject) ----
    for iff = 1:numel(funcs)
        Fn = funcs{iff};  seriesPsd = {};
        for i = 1:numel(sSrc)
            sSer = bst_process('CallProcess', 'process_source_flow_functional', sSrc(i), [], ...
                'func', Fn, 'scouts', {});
            sFp = bst_process('CallProcess', 'process_psd', sSer, [], ...
                'timewindow', [], 'win_length', 8, 'win_overlap', 50, 'sensortypes', '', ...
                'edit', struct('Comment', ['FuncPSD_' Fn '_' label], 'TimeBands', [], 'Freqs', {bands}, ...
                    'ClusterFuncTime', 'none', 'Measure', 'power', 'Output', 'all', 'SaveKernel', 0));
            seriesPsd{end+1} = sFp.FileName; %#ok<AGROW>
            fprintf('  [%s] %s -> series -> band PSD\n', Fn, subjOf{i});
        end
        OutputFiles.(['func_' Fn]) = seriesPsd;
    end

    ReportFile = bst_report('Save', []);
    fprintf('tutorial_flow_omega: done. Report %s\n', ReportFile);
end
