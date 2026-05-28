function varargout = process_eigenmodes_wavelet( varargin )
% PROCESS_EIGENMODES_WAVELET: Complex Morlet wavelet tensor of eigenmode coefficients.
%
% USAGE:  sProcess = process_eigenmodes_wavelet('GetDescription')
%       OutputFiles = process_eigenmodes_wavelet('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Applies the complex Morlet CWT to an eigenmode-coefficient matrix
%     (matrix_eigentransform, [K x nTime]) and saves the time-resolved
%     (lambda, omega, t) tensor W_k(s,t) as a complex, lambda-labeled timefreq
%     file (amplitude + phase), viewable in Brainstorm's time-frequency viewer.
%
% SEE ALSO: bst_eigenmodes_wavelet, process_eigenmodes_transform

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
    sProcess.Comment     = 'Eigenmode wavelet tensor (complex)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 336.7;   % after the noise-floor denoise (336.6)
    sProcess.Description = '';
    sProcess.InputTypes  = {'matrix'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    sProcess.options.flo.Comment    = 'Lowest frequency: ';
    sProcess.options.flo.Type       = 'value';
    sProcess.options.flo.Value      = {2, 'Hz', 2};

    sProcess.options.fhi.Comment    = 'Highest frequency (0 = auto, min(100, 0.4*sfreq)): ';
    sProcess.options.fhi.Type       = 'value';
    sProcess.options.fhi.Value      = {0, 'Hz', 2};

    sProcess.options.nfreqs.Comment = 'Number of frequencies (log-spaced): ';
    sProcess.options.nfreqs.Type    = 'value';
    sProcess.options.nfreqs.Value   = {40, '', 0};

    sProcess.options.morletfc.Comment     = 'Morlet central frequency Fc: ';
    sProcess.options.morletfc.Type        = 'value';
    sProcess.options.morletfc.Value       = {1, 'Hz', 2};

    sProcess.options.morletfwhmtc.Comment = 'Morlet time resolution (FWHM): ';
    sProcess.options.morletfwhmtc.Type    = 'value';
    sProcess.options.morletfwhmtc.Value   = {3, 's', 2};

    sProcess.options.label_info.Comment = ['<FONT color="#777777">Complex Morlet wavelet tensor W_k(s,t) of the ' ...
        'eigenmode coefficients &mdash; the time-resolved (&lambda;,&omega;,t) decomposition (amplitude + phase).<BR>' ...
        'Input: an eigenmode-coefficient matrix (matrix_eigentransform).</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU,INUSD>
    Comment = 'Eigenmode wavelet tensor (complex)';
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    flo  = sProcess.options.flo.Value{1};
    fhi  = sProcess.options.fhi.Value{1};
    nfq  = sProcess.options.nfreqs.Value{1};
    fc   = sProcess.options.morletfc.Value{1};
    fwhm = sProcess.options.morletfwhmtc.Value{1};

    for iInput = 1:numel(sInputs)
        sInput = sInputs(iInput);
        M = in_bst_matrix(sInput.FileName);
        Coeffs = double(M.Value);                       % [K x nTime]
        if size(Coeffs, 2) < 2
            bst_report('Error', sProcess, sInput, 'Need at least 2 time samples to estimate a wavelet transform.');
            continue;
        end
        sfreq = 1 / (M.Time(2) - M.Time(1));
        fhiEff = fhi;
        if fhiEff <= 0
            fhiEff = min(100, 0.4 * sfreq);
        end
        if flo <= 0
            bst_report('Error', sProcess, sInput, 'Lowest frequency must be positive (default: 2 Hz).');
            continue;
        end
        if flo >= fhiEff
            bst_report('Error', sProcess, sInput, sprintf('Lowest frequency (%.1f Hz) must be below the highest (%.1f Hz).', flo, fhiEff));
            continue;
        end
        Freqs = logspace(log10(flo), log10(fhiEff), nfq);

        [W, Freqs] = bst_eigenmodes_wavelet(Coeffs, sfreq, Freqs, 'MorletFc', fc, 'MorletFwhmTc', fwhm);

        % Row labels: reuse the matrix Description (mode labels) or build them.
        RowNames = M.Description;
        if isempty(RowNames) || numel(RowNames) ~= size(W,1)
            RowNames = cell(size(W,1), 1);
            for k = 1:size(W,1), RowNames{k} = sprintf('Mode %d', k); end
        end

        [sStudyOut, iStudyOut] = bst_get('Study', sInput.iStudy);
        StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

        TFmat = db_template('timefreqmat');
        TFmat.TF       = W;                              % complex [K x nTime x nFreq]
        TFmat.Time     = M.Time;
        TFmat.Freqs    = Freqs(:)';
        TFmat.RowNames = RowNames;
        TFmat.Measure  = 'none';                         % complex (un-measured)
        TFmat.Method   = 'morlet';
        TFmat.DataType = 'matrix';
        if isfield(M, 'SurfaceFile')
            TFmat.SurfaceFile = M.SurfaceFile;
        end
        TFmat.Comment  = sprintf('EigenWavelet (%d modes, %d freqs) | %s', size(W,1), numel(Freqs), sInput.Comment);
        TFmat.nAvg     = 1;
        TFmat = bst_history('add', TFmat, 'eigenmodes_wavelet', TFmat.Comment);

        OutFile = bst_process('GetNewFilename', StudyDir, 'timefreq_eigenwavelet');
        bst_save(OutFile, TFmat, 'v6');
        db_add_data(iStudyOut, OutFile, TFmat);
        OutputFiles{end+1} = file_short(OutFile); %#ok<AGROW>

        bst_report('Info', sProcess, sInput, sprintf('Wavelet tensor: %d modes x %d times x %d freqs (%.1f-%.1f Hz).', ...
            size(W,1), size(W,2), size(W,3), Freqs(1), Freqs(end)));
    end
end
