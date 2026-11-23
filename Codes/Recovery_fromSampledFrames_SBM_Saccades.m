function [Acq_Phase_Matrix, Recovered_Frames, RecoveryInfo] = ...
    Recovery_fromSampledFrames_SBM_Saccades( ...
        nTurnSamples, T_excitation, Ts_OCT, nFrames, ...
        SampledFrames_raw, phaseStepSamples, ...
        nAscans_acq, nBscans_acq, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 25-06-2026
% Date of last modification (DD/MM/YY): 25-06-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recovery_fromSampledFrames_SBM_Saccades
%
% Doppler / acquisition-phase compensation for SBM sampled movies after
% spatial upsampling.
%
% Key idea:
%   The acquisition phase is physically defined by the ORIGINAL SBM sampling
%   pattern:
%
%       nBscans_acq x nAscans_acq
%
%   not by the upsampled grid size.
%
%   Therefore, this function:
%
%       1) builds the low-resolution acquisition phase matrix using the
%          original acquisition geometry,
%
%       2) upsamples the phase matrix to match size(SampledFrames_raw),
%          using complex phasor interpolation to avoid phase wrapping errors,
%
%       3) compensates the excitation-frequency FFT bin.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Parse options =============================

p = inputParser;

addParameter(p, 'FrameConvention', 'raw');      % 'raw' or 'vis'
addParameter(p, 'PhaseInterpMethod', 'linear'); % 'linear', 'nearest', 'cubic', etc.
addParameter(p, 'MakeDiagnostics', true);

parse(p, varargin{:});
Opt = p.Results;

FrameConvention = lower(string(Opt.FrameConvention));

if ~ismember(FrameConvention, ["raw","vis"])
    error('FrameConvention must be either ''raw'' or ''vis''.');
end

%% ============================ Input sizes ===============================

[nRows_out, nCols_out, nFrames_in] = size(SampledFrames_raw);

if nFrames_in ~= nFrames
    warning('Input nFrames = %d, but size(SampledFrames_raw,3) = %d. Using movie size.', ...
        nFrames, nFrames_in);
    nFrames = nFrames_in;
end

fprintf('\nRecovery_fromSampledFrames_SBM_Saccades()\n');
fprintf('Input sampled movie size      : [%d x %d x %d]\n', ...
    nRows_out, nCols_out, nFrames);
fprintf('Original acquisition geometry : [%d B-scans x %d A-scans]\n', ...
    nBscans_acq, nAscans_acq);

%% =================== Build low-resolution acquisition phase ==============

N_singleframe = nAscans_acq*nBscans_acq + ...
                (nBscans_acq - 1)*nTurnSamples;

f_excitation = 1 / T_excitation;
omega_excitation = 2*pi*f_excitation;

timeline_singleframe = (1:N_singleframe) * Ts_OCT;

rowStarts = 1 + (0:(nBscans_acq-1)) * (nAscans_acq + nTurnSamples);

Acq_Phase_LowRes = nan(nBscans_acq, nAscans_acq);

for r = 1:nBscans_acq

    idx_row = rowStarts(r) + (0:(nAscans_acq-1));

    valid = idx_row >= 1 & idx_row <= N_singleframe;

    t_vals = timeline_singleframe(idx_row(valid));

    Acq_Phase_LowRes(r,valid) = omega_excitation * t_vals;
end

Acq_Phase_LowRes = angle(exp(1j * Acq_Phase_LowRes));

%% ============ If frames are visual, flip phase on serpentine rows =========
% raw convention:
%   each row is kept in acquisition order.
%
% vis convention:
%   rows acquired right-to-left were flipped for visualization, so the
%   acquisition phase must be flipped in those rows too.

if FrameConvention == "vis"

    for r = 1:nBscans_acq
        if mod(r,2) == 0
            Acq_Phase_LowRes(r,:) = fliplr(Acq_Phase_LowRes(r,:));
        end
    end
end

%% ============ Upsample acquisition phase using phasor interpolation =======
% Direct interpolation of wrapped phase is wrong near -pi/pi.
% Therefore, interpolate exp(1j*phase), then recover angle.

Acq_Phase_Matrix = local_upsample_phase_phasor( ...
    Acq_Phase_LowRes, ...
    nRows_out, nCols_out, ...
    Opt.PhaseInterpMethod);

fprintf('Upsampled phase matrix size  : [%d x %d]\n', ...
    size(Acq_Phase_Matrix,1), size(Acq_Phase_Matrix,2));

%% =================== FFT-based phase compensation ========================

Nfft = nFrames;

dt_SampledFrames = double(phaseStepSamples) * Ts_OCT;

F_SampledFrames = fft(double(SampledFrames_raw), Nfft, 3);

freqs_SampledFrames = (0:(Nfft-1))' * (1/dt_SampledFrames / Nfft);

[~, kpos_SampledFrames] = min(abs(freqs_SampledFrames - f_excitation));

Ph_SampledFrames = squeeze(F_SampledFrames(:,:,kpos_SampledFrames));

Ph_Compensated = Ph_SampledFrames .* exp(1j * Acq_Phase_Matrix);

F_compensated = F_SampledFrames;
F_compensated(:,:,kpos_SampledFrames) = Ph_Compensated;

kpos = kpos_SampledFrames;
kneg = mod(Nfft - kpos + 2 - 1, Nfft) + 1;

isDC = (kpos == 1);
isNyquist = (mod(Nfft,2) == 0) && (kpos == (Nfft/2 + 1));

if ~isDC && ~isNyquist
    F_compensated(:,:,kneg) = conj(Ph_Compensated);
end

Recovered_Frames = ifft(F_compensated, [], 3, 'symmetric');

%% ============================== Info ====================================

RecoveryInfo = struct();

RecoveryInfo.Method = 'SBM recovery using low-resolution acquisition phase upsampled as complex phasor';
RecoveryInfo.FrameConvention = FrameConvention;

RecoveryInfo.InputMovieSize = size(SampledFrames_raw);
RecoveryInfo.OutputMovieSize = size(Recovered_Frames);

RecoveryInfo.nAscans_acq = nAscans_acq;
RecoveryInfo.nBscans_acq = nBscans_acq;
RecoveryInfo.nTurnSamples = nTurnSamples;

RecoveryInfo.N_singleframe = N_singleframe;
RecoveryInfo.T_excitation = T_excitation;
RecoveryInfo.Ts_OCT = Ts_OCT;
RecoveryInfo.phaseStepSamples = phaseStepSamples;
RecoveryInfo.dt_SampledFrames = dt_SampledFrames;

RecoveryInfo.f_excitation = f_excitation;
RecoveryInfo.kpos_SampledFrames = kpos_SampledFrames;
RecoveryInfo.kneg_SampledFrames = kneg;
RecoveryInfo.freq_bin_selected = freqs_SampledFrames(kpos_SampledFrames);

RecoveryInfo.Acq_Phase_LowRes = Acq_Phase_LowRes;
RecoveryInfo.Acq_Phase_Matrix = Acq_Phase_Matrix;

fprintf('Selected FFT bin frequency : %.4f Hz\n', ...
    RecoveryInfo.freq_bin_selected);
fprintf('Excitation frequency       : %.4f Hz\n', f_excitation);
fprintf('Recovery complete.\n\n');

%% ============================== Diagnostics ==============================

if Opt.MakeDiagnostics

    figure('Name','SBM Saccadic Recovery Phase Compensation', ...
        'Units','normalized', ...
        'Position',[0.12 0.12 0.75 0.55]);

    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(Acq_Phase_LowRes);
    axis image;
    colormap(gca, hsv);
    colorbar;
    title('Low-res acquisition phase');
    xlabel('A-scan index');
    ylabel('B-scan index');

    nexttile;
    imagesc(Acq_Phase_Matrix);
    axis image;
    colormap(gca, hsv);
    colorbar;
    title('Upsampled acquisition phase');
    xlabel('Column index');
    ylabel('Row index');

    nexttile;
    imagesc(Recovered_Frames(:,:,round(nFrames/2)));
    axis image;
    colormap(gca, turbo);
    colorbar;
    title(sprintf('Recovered frame %d', round(nFrames/2)));
    xlabel('Column index');
    ylabel('Row index');
end

end

%% ===================================================================== %%
function PhaseOut = local_upsample_phase_phasor(PhaseIn, nRows_out, nCols_out, interpMethod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Upsample a wrapped phase map robustly.
%
% Instead of interpolating phase directly, interpolate:
%
%       cos(phase) and sin(phase)
%
% or equivalently:
%
%       exp(1j*phase)
%
% Then recover the final phase using:
%
%       angle(real_interp + 1j*imag_interp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nRows_in, nCols_in] = size(PhaseIn);

row_in = linspace(1, nRows_in, nRows_in).';
col_in = linspace(1, nCols_in, nCols_in).';

row_out = linspace(1, nRows_in, nRows_out).';
col_out = linspace(1, nCols_in, nCols_out).';

[Rq, Cq] = ndgrid(row_out, col_out);

PhasePhasor = exp(1j * PhaseIn);

F_real = griddedInterpolant( ...
    {row_in, col_in}, ...
    real(PhasePhasor), ...
    interpMethod, ...
    'nearest');

F_imag = griddedInterpolant( ...
    {row_in, col_in}, ...
    imag(PhasePhasor), ...
    interpMethod, ...
    'nearest');

RealOut = F_real(Rq, Cq);
ImagOut = F_imag(Rq, Cq);

PhaseOut = angle(RealOut + 1j*ImagOut);

end