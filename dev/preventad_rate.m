function [fsTarget, notchList] = preventad_rate(fsFiles)
% PREVENTAD_RATE: the one sampling rate a subject's recordings are brought to, and its notches.
%
% USAGE:  [fsTarget, notchList] = preventad_rate(fsFiles)
%
% INPUT:  fsFiles   - [1 x nFiles] sampling rates (Hz) of every recording of the subject
%                     (task-rest AND task-noise, all sessions)
% OUTPUT: fsTarget  - 1200 Hz, or the slowest recording's rate when that is lower
%         notchList - the 60 Hz harmonics up to 300 Hz that lie below fsTarget/2
%
% ⚠ Noise and rest must share one rate and one band, or the noise covariance misdescribes the
% rest data. A 600 Hz noise run (OMEGA sub-PD1306) holds nothing above 300 Hz; upsampling it to
% 1200 would leave the 300-600 Hz noise of the rest runs out of the covariance, so the whole
% subject goes to 600 Hz instead. process_notch stops at the first frequency >= Nyquist.
%
% See also: preventad_subject, process_resample, process_notch

% Author: Diellor Basha, 2026
fsTarget  = min([1200, fsFiles(:)']);
harmonics = 60:60:300;
notchList = harmonics(harmonics < fsTarget / 2);
end
