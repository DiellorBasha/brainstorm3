function varargout = process_eigenmodes_dispersion( varargin )
% PROCESS_EIGENMODES_DISPERSION: Wave-vs-diffusion regime from a (lambda,omega) spectrum.
%
% USAGE:  sProcess = process_eigenmodes_dispersion('GetDescription')
%       OutputFiles = process_eigenmodes_dispersion('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Fits a wave model (peak frequency ~ sqrt(lambda) -> speed c) and a
%     diffusion model (bandwidth ~ lambda -> diffusivity alpha) to an eigenmode
%     (lambda, omega) power spectrum (e.g. an FFT/PSD of the eigenmode
%     coefficients, or the time-averaged wavelet tensor), and reports the
%     better-fitting regime. See bst_eigenmodes_dispersion.
%
% SEE ALSO: bst_eigenmodes_dispersion, process_eigenmodes_wavelet

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
    sProcess.Comment     = 'Eigenmode dispersion (wave vs diffusion)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 336.8;   % after the wavelet tensor (336.7)
    sProcess.Description = '';
    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    sProcess.options.minpowerfrac.Comment = 'Drop modes below this fraction of max per-mode power: ';
    sProcess.options.minpowerfrac.Type    = 'value';
    sProcess.options.minpowerfrac.Value   = {0, '', 2};

    sProcess.options.label_info.Comment = ['<FONT color="#777777">Fits a wave (peak freq &prop; &radic;&lambda;) ' ...
        'and a diffusion (bandwidth &prop; &lambda;) model to the (&lambda;,&omega;) spectrum and reports the ' ...
        'better-fitting regime + speed c / diffusivity &alpha;.<BR>Input: an eigenmode (&lambda;,&omega;) timefreq spectrum.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    mpf = sProcess.options.minpowerfrac.Value{1};
    if mpf > 0
        Comment = sprintf('Eigenmode dispersion (wave vs diffusion, minpow=%.2g)', mpf);
    else
        Comment = 'Eigenmode dispersion (wave vs diffusion)';
    end
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    MinPowerFrac = sProcess.options.minpowerfrac.Value{1};

    for iInput = 1:numel(sInputs)
        sInput = sInputs(iInput);
        TFmat = in_bst_timefreq(sInput.FileName, 1);   % LoadFull=1

        if iscell(TFmat.Freqs)
            bst_report('Error', sProcess, sInput, 'Frequency bands (cell Freqs) are not supported; use a frequency vector.');
            continue;
        end
        % Convert to power per the file's Measure (do not square already-power data)
        switch lower(TFmat.Measure)
            case 'power'
                Pw = abs(TFmat.TF);
            otherwise   % 'none' (complex) or 'magnitude'
                Pw = abs(TFmat.TF).^2;
        end
        P = reshape(mean(Pw, 2), size(Pw,1), size(Pw,3));   % [K x nFreq], averaged over time

        if isempty(TFmat.SurfaceFile)
            bst_report('Error', sProcess, sInput, 'Timefreq has no SurfaceFile; cannot obtain eigenvalues.');
            continue;
        end
        [Em, isComputed] = in_tess_eigenmodes(TFmat.SurfaceFile);
        if ~isComputed
            bst_report('Error', sProcess, sInput, ['No eigenmodes on surface: ' TFmat.SurfaceFile]);
            continue;
        end
        K = size(P, 1);
        if size(Em.Values, 1) < K
            bst_report('Error', sProcess, sInput, sprintf('Spectrum has %d modes but surface has only %d eigenvalues.', K, size(Em.Values,1)));
            continue;
        end
        lambdas = double(Em.Values(1:K));

        Out = bst_eigenmodes_dispersion(P, lambdas, TFmat.Freqs, 'MinPowerFrac', MinPowerFrac);

        bst_report('Info', sProcess, sInput, sprintf(...
            'Dispersion: regime=%s | c=%.3g m/s, alpha=%.3g | R2wave=%.3f, R2diff=%.3f (margin %.3f)', ...
            Out.Regime, Out.c, Out.alpha, Out.R2wave, Out.R2diff, Out.Margin));

        % Save per-mode features as a matrix file
        [sStudyOut, iStudyOut] = bst_get('Study', sInput.iStudy);
        StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));
        M = db_template('matrixmat');
        M.Value       = [Out.PeakFreq'; Out.Bandwidth'];   % [2 x K]
        M.Time        = 1:K;
        M.Description  = {'PeakFreq (Hz)'; 'Bandwidth (Hz)'};
        M.SurfaceFile = TFmat.SurfaceFile;
        M.Comment     = sprintf('EigenDispersion [%s] c=%.3g a=%.3g R2w=%.2f R2d=%.2f | %s', ...
            Out.Regime, Out.c, Out.alpha, Out.R2wave, Out.R2diff, sInput.Comment);
        M.nAvg        = 1;
        M = bst_history('add', M, 'eigenmodes_dispersion', M.Comment);

        OutFile = bst_process('GetNewFilename', StudyDir, 'matrix_eigendispersion');
        bst_save(OutFile, M, 'v6');
        db_add_data(iStudyOut, OutFile, M);
        OutputFiles{end+1} = file_short(OutFile); %#ok<AGROW>
    end
end
