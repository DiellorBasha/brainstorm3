function varargout = process_tess_manifold( varargin )
% PROCESS_TESS_MANIFOLD: Associate an nxr manifold node with a subject's cortex.
%
% Thin pipeline wrapper around tess_manifold: for every subject represented in
% the input files, resolves the subject's cortex surface and computes (or reuses)
% a manifold_ DB node on it (hemisphere-split facets + DEC operators via
% nxr-compute). tess_manifold is find-or-load-or-create, so this process is
% idempotent - re-running reuses an existing same-gauge node unless ForceRecompute.
%
% The manifold is anatomy-level (one per subject's cortex), so the input files
% are only used to resolve the subjects; they are passed through unchanged.
%
% SEE ALSO: tess_manifold, db_add_manifold, view_manifold

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

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Compute manifold (nxr geometry)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 321;
    sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/HeadModel';
    sProcess.InputTypes  = {'data', 'raw', 'results', 'matrix'};
    sProcess.OutputTypes = {'data', 'raw', 'results', 'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Label
    sProcess.options.label.Comment = ['<HTML><B>Cortical manifold (nxr geometry backbone)</B><BR>' ...
        'Computes a manifold_ node on each subject''s cortex surface:<BR>' ...
        'hemisphere-split facets + DEC operators (requires nxr-compute).'];
    sProcess.options.label.Type    = 'label';
    % Gauge
    sProcess.options.gauge.Comment = {'trivial', 'levi-civita', 'euclidean'; 'trivial', 'levi-civita', 'euclidean'};
    sProcess.options.gauge.Type    = 'radio_label';
    sProcess.options.gauge.Value   = 'trivial';
    % Force recompute
    sProcess.options.forcerecompute.Comment = 'Force recompute (ignore existing manifold node)';
    sProcess.options.forcerecompute.Type    = 'checkbox';
    sProcess.options.forcerecompute.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sprintf('%s (gauge=%s)', sProcess.Comment, sProcess.options.gauge.Value);
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    % Pass the input files through unchanged (manifold is anatomy-level)
    OutputFiles = {sInputs.FileName};

    % ---- options ----
    Gauge = sProcess.options.gauge.Value;
    Force = logical(sProcess.options.forcerecompute.Value);

    % ---- resolve the unique subjects represented in the inputs ----
    % (manifold attaches to the subject's cortex, one per subject)
    SubjectFiles = cell(1, numel(sInputs));
    for i = 1:numel(sInputs)
        sStudy = bst_get('Study', sInputs(i).iStudy);
        SubjectFiles{i} = sStudy.BrainStormSubject;
    end
    SubjectFiles = unique(SubjectFiles);

    % ---- one manifold per subject cortex ----
    for i = 1:numel(SubjectFiles)
        sSubject = bst_get('Subject', SubjectFiles{i});
        if isempty(sSubject) || isempty(sSubject.iCortex)
            bst_report('Warning', sProcess, sInputs, ...
                sprintf('No cortex surface for subject "%s"; skipping manifold.', SubjectFiles{i}));
            continue;
        end
        CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
        bst_progress('text', sprintf('Computing manifold (gauge=%s): %s', Gauge, sSubject.Name));
        try
            tess_manifold(CortexFile, 'Gauge', Gauge, 'ForceRecompute', Force);
        catch ME
            bst_report('Error', sProcess, sInputs, ...
                sprintf('tess_manifold failed for "%s": %s', sSubject.Name, ME.message));
        end
    end
end
