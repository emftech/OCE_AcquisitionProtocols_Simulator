function [Exc_Trigger,Acq_Trigger,GalvoX,GalvoY] = AcqVectors_SBM_v2( ...
    horizontalMin,horizontalMax,verticalMin,verticalMax, ...
    nAscans,nBscans,nFrames, ...
    startDelay_samples,nTurnSamples,AcqPulseWidth, ...
    ExcPulseWidth,T_excitation,Ts_OCT,phaseStepSamples)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 17-06-2026
% Date of last modification (DD/MM/YY): 17-06-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AcqVectors_SBM_v2
%
% Generate acquisition vectors for the Synchronized B-M Mode (SBM)
% acquisition protocol using the generalized sampling-framework convention.
%
% Outputs:
%
%   Exc_Trigger:
%       Excitation TTL vector. It is ON whenever a new excitation period is
%       initiated, i.e. at the beginning of each B-scan line within each SBM
%       frame. The first B-scan of every frame also corresponds to the
%       beginning of the full BM frame.
%
%   Acq_Trigger:
%       Acquisition TTL vector. It is ON whenever the OCT system is sampling
%       a spatial A-line position. It is OFF during:
%
%           - B-scan turn-around transitions,
%           - inter-frame transitions,
%           - excitation/start-delay holding periods,
%           - phase-offset holding periods.
%
%   GalvoX, GalvoY:
%       Galvanometer trajectories in physical coordinates [mm].
%
% Important convention:
%
%   AcqPulseWidth is kept in the input list for compatibility with previous
%   scripts, but in v2 the acquisition trigger is not treated as a short
%   pulse. Instead, Acq_Trigger is high over all valid spatial sampling
%   samples.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    horizontalMin (1,1) double
    horizontalMax (1,1) double
    verticalMin   (1,1) double
    verticalMax   (1,1) double
    nAscans       (1,1) double {mustBeInteger,mustBePositive}
    nBscans       (1,1) double {mustBeInteger,mustBePositive}
    nFrames       (1,1) double {mustBeInteger,mustBePositive}
    startDelay_samples (1,1) double {mustBeInteger,mustBeNonnegative}
    nTurnSamples  (1,1) double {mustBeInteger,mustBeNonnegative}
    AcqPulseWidth (1,1) double {mustBeInteger,mustBePositive}
    ExcPulseWidth (1,1) double {mustBeInteger,mustBePositive}
    T_excitation  (1,1) double {mustBePositive}
    Ts_OCT        (1,1) double {mustBePositive}
    phaseStepSamples double = []
end

%% ======================= Phase-step definition ==========================

if isempty(phaseStepSamples)
    rawPhaseStep = (T_excitation / nFrames) / Ts_OCT;

    if abs(rawPhaseStep - round(rawPhaseStep)) > 1e-10
        error(['Default phaseStepSamples is not integer. Ensure ' ...
               '(T_excitation/nFrames) is an integer multiple of Ts_OCT, ' ...
               'or provide phaseStepSamples explicitly. Computed value = %g.'], ...
               rawPhaseStep);
    end

    phaseStepSamples = round(rawPhaseStep);
else
    if ~isscalar(phaseStepSamples) || ...
            abs(phaseStepSamples - round(phaseStepSamples)) > 1e-10
        error('phaseStepSamples must be a scalar integer number of OCT samples.');
    end

    phaseStepSamples = round(phaseStepSamples);
end

samplesPerExcPeriod = round(T_excitation / Ts_OCT);

fprintf('\nAcqVectors_SBM_v2()\n');
fprintf('samplesPerExcPeriod = %d samples (T_exc = %.3f ms)\n', ...
    samplesPerExcPeriod, 1e3*T_excitation);
fprintf('phaseStepSamples    = %d samples\n', phaseStepSamples);
fprintf('AcqPulseWidth input retained for compatibility; continuous A-line acquisition is used.\n');

%% ======================= Single-frame raster pattern =====================

% Normalized raster coordinates.
% X runs left-to-right / right-to-left depending on serpentine row.
% Y runs from top to bottom.
y_norm = linspace(1, -1, nBscans);

Nsingle_est = nAscans*nBscans + nTurnSamples*(nBscans - 1);

x_positions  = nan(1, Nsingle_est);
y_positions  = nan(1, Nsingle_est);
isAline_mask = false(1, Nsingle_est);

block_stride = nAscans + nTurnSamples;

% Insert A-line samples for each B-scan.
for jj = 1:nBscans

    if mod(jj,2) == 1
        x_row = linspace(-1, 1, nAscans);
    else
        x_row = linspace(1, -1, nAscans);
    end

    y_row = repmat(y_norm(jj), 1, nAscans);

    idx_start = (jj-1)*block_stride + 1;
    idx_end   = idx_start + nAscans - 1;

    x_positions(idx_start:idx_end)  = x_row;
    y_positions(idx_start:idx_end)  = y_row;
    isAline_mask(idx_start:idx_end) = true;
end

% Smooth turn-around transitions between B-scans.
for jj = 1:(nBscans-1)

    idx_start = (jj-1)*block_stride + 1;
    idx_end   = idx_start + nAscans - 1;

    next_idx_start = jj*block_stride + 1;
    trans_idx = (idx_end+1):(next_idx_start-1);

    if isempty(trans_idx)
        continue;
    end

    p1 = [x_positions(idx_end),        y_positions(idx_end)];
    p2 = [x_positions(next_idx_start), y_positions(next_idx_start)];

    if mod(jj,2) == 1
        dir1 = +1;
    else
        dir1 = -1;
    end

    if mod(jj+1,2) == 1
        dir2 = +1;
    else
        dir2 = -1;
    end

    tanMag = 0.2;

    t1 = [dir1*tanMag, 0];
    t2 = [dir2*tanMag, 0];

    nT = numel(trans_idx);
    s = (1:nT)' ./ (nT + 1);

    h00 =  2*s.^3 - 3*s.^2 + 1;
    h10 =      s.^3 - 2*s.^2 + s;
    h01 = -2*s.^3 + 3*s.^2;
    h11 =      s.^3 -   s.^2;

    x_positions(trans_idx) = h00*p1(1) + h10*t1(1) + h01*p2(1) + h11*t2(1);
    y_positions(trans_idx) = h00*p1(2) + h10*t1(2) + h01*p2(2) + h11*t2(2);
end

% Safety interpolation in case any transition samples remained undefined.
nan_mask = isnan(x_positions) | isnan(y_positions);

if any(nan_mask)
    pattern_idx = 1:numel(x_positions);

    valid_x = ~isnan(x_positions);
    valid_y = ~isnan(y_positions);

    x_positions(nan_mask) = interp1( ...
        pattern_idx(valid_x), x_positions(valid_x), ...
        pattern_idx(nan_mask), 'linear', 'extrap');

    y_positions(nan_mask) = interp1( ...
        pattern_idx(valid_y), y_positions(valid_y), ...
        pattern_idx(nan_mask), 'linear', 'extrap');
end

% Map normalized coordinates to physical coordinates [mm].
X_single = (((x_positions(:) + 1)/2) * ...
    (horizontalMax - horizontalMin) + horizontalMin);

Y_single = (((y_positions(:) + 1)/2) * ...
    (verticalMax - verticalMin) + verticalMin);

Aline_single = isAline_mask(:);

Nsingle = numel(X_single);

% Relative starts of each B-scan within the single-frame raster.
row_start_idx = (0:(nBscans-1))*block_stride + 1;
line_rel_starts = row_start_idx - row_start_idx(1);

%% ======================= Inter-frame transition length ===================

Transition_Btw_Frames = round(2 * T_excitation / Ts_OCT);

if Transition_Btw_Frames < 1
    Transition_Btw_Frames = 1;
end

% Force inter-frame transition length to be an integer number of excitation
% periods, preserving synchronization across frames.
if samplesPerExcPeriod > 0
    mult = max(1, round(Transition_Btw_Frames / samplesPerExcPeriod));
    Transition_Btw_Frames = mult * samplesPerExcPeriod;
end

%% ======================= Build full timeline =============================

GalvoX = [];
GalvoY = [];
Exc_Trigger = [];
Acq_Trigger = [];

global_idx = 0;

holdX = X_single(1);
holdY = Y_single(1);

for m = 1:nFrames

    % ---------------------------------------------------------------------
    % Optional frame-level start delay.
    % The galvos remain at the first raster position and no acquisition
    % occurs during this interval.
    % ---------------------------------------------------------------------
    if startDelay_samples > 0
        GalvoX = [GalvoX; repmat(holdX, startDelay_samples, 1)];
        GalvoY = [GalvoY; repmat(holdY, startDelay_samples, 1)];

        Exc_Trigger = [Exc_Trigger; zeros(startDelay_samples,1)];
        Acq_Trigger = [Acq_Trigger; zeros(startDelay_samples,1)];

        global_idx = global_idx + startDelay_samples;
    end

    % ---------------------------------------------------------------------
    % Frame-dependent phase offset.
    % This preserves the synchronized phase-stepping logic used in SBM.
    % ---------------------------------------------------------------------
    frame_phase_offset = round((m-1)*phaseStepSamples);

    if frame_phase_offset > 0
        GalvoX = [GalvoX; repmat(holdX, frame_phase_offset, 1)];
        GalvoY = [GalvoY; repmat(holdY, frame_phase_offset, 1)];

        Exc_Trigger = [Exc_Trigger; zeros(frame_phase_offset,1)];
        Acq_Trigger = [Acq_Trigger; zeros(frame_phase_offset,1)];

        global_idx = global_idx + frame_phase_offset;
    end

    % ---------------------------------------------------------------------
    % Append one full SBM raster frame.
    % Acquisition is ON only for true A-line samples.
    % Turn-around samples are OFF.
    % ---------------------------------------------------------------------
    scan_global_start = global_idx + 1;

    GalvoX = [GalvoX; X_single];
    GalvoY = [GalvoY; Y_single];

    Exc_block = zeros(Nsingle,1);
    Acq_block = double(Aline_single);

    % Excitation pulses at the beginning of each B-scan.
    % jj = 1 is also the beginning of the full BM frame.
    for jj = 1:nBscans

        exc_idx_local = line_rel_starts(jj) + 1;

        pulse_idx = exc_idx_local : ...
            min(exc_idx_local + ExcPulseWidth - 1, Nsingle);

        Exc_block(pulse_idx) = 1;
    end

    Exc_Trigger = [Exc_Trigger; Exc_block];
    Acq_Trigger = [Acq_Trigger; Acq_block];

    global_idx = global_idx + Nsingle;

    % ---------------------------------------------------------------------
    % Smooth transition from the end of the current frame back to the first
    % point of the next frame. No acquisition or excitation occurs here.
    % ---------------------------------------------------------------------
    if m < nFrames

        L = Transition_Btw_Frames;

        p0x = GalvoX(end);
        p0y = GalvoY(end);

        p1x = X_single(1);
        p1y = Y_single(1);

        if numel(GalvoX) >= 2
            d0x = GalvoX(end) - GalvoX(end-1);
            d0y = GalvoY(end) - GalvoY(end-1);
        else
            d0x = 0;
            d0y = 0;
        end

        if numel(X_single) >= 2
            d1x = X_single(2) - X_single(1);
            d1y = Y_single(2) - Y_single(1);
        else
            d1x = 0;
            d1y = 0;
        end

        m0x = d0x * L;
        m0y = d0y * L;

        m1x = d1x * L;
        m1y = d1y * L;

        s = (1:L)' ./ (L + 1);

        h00 =  2*s.^3 - 3*s.^2 + 1;
        h10 =      s.^3 - 2*s.^2 + s;
        h01 = -2*s.^3 + 3*s.^2;
        h11 =      s.^3 -   s.^2;

        Xtrans = h00*p0x + h10*m0x + h01*p1x + h11*m1x;
        Ytrans = h00*p0y + h10*m0y + h01*p1y + h11*m1y;

        GalvoX = [GalvoX; Xtrans];
        GalvoY = [GalvoY; Ytrans];

        Exc_Trigger = [Exc_Trigger; zeros(L,1)];
        Acq_Trigger = [Acq_Trigger; zeros(L,1)];

        global_idx = global_idx + L;
    end
end

%% ======================= Final formatting ================================

GalvoX = double(GalvoX(:));
GalvoY = double(GalvoY(:));
Exc_Trigger = double(Exc_Trigger(:));
Acq_Trigger = double(Acq_Trigger(:));

N_total = numel(GalvoX);
totalAcqTime_s  = N_total * Ts_OCT;
totalAcqTime_ms = totalAcqTime_s * 1e3;

fprintf('Built SBM timeline:\n');
fprintf('  Total protocol time    = %.3f ms (%.4f s)\n', ...
    totalAcqTime_ms, totalAcqTime_s);
fprintf('  Total samples          = %d\n', N_total);
fprintf('  nFrames                = %d\n', nFrames);
fprintf('  Single-frame samples   = %d\n', Nsingle);
fprintf('  A-line samples/frame   = %d\n', nnz(Aline_single));
fprintf('  Turn samples/frame     = %d\n', nnz(~Aline_single));
fprintf('  Inter-frame transition = %d samples\n', Transition_Btw_Frames);
fprintf('  Excitation pulses      = %d\n', nnz(diff([0; Exc_Trigger]) == 1));
fprintf('  Acquisition samples    = %d\n', nnz(Acq_Trigger > 0));

%% ======================= Diagnostics ====================================

t_global_ms = ((0:(N_total-1))' * Ts_OCT) * 1e3;

figure('Name','SBM Single-Frame Raster Pattern', ...
    'Units','normalized', ...
    'Position',[0.05 0.08 0.55 0.62]);

hold on;

plot(X_single, Y_single, '-', ...
    'Color',[0.85 0.85 0.85], ...
    'LineWidth',1.0);

plot(X_single(~Aline_single), Y_single(~Aline_single), '.', ...
    'Color',[0.6 0 0], ...
    'MarkerSize',9);

plot(X_single(Aline_single), Y_single(Aline_single), 'bo', ...
    'MarkerFaceColor','b', ...
    'MarkerSize',5);

axis equal;
grid on;

xlabel('X [mm]');
ylabel('Y [mm]');

title(sprintf('SBM Raster Pattern: %d A-lines, %d transition samples', ...
    nnz(Aline_single), nnz(~Aline_single)));

if nTurnSamples ~= 0
    legend('Raster path','Turn samples','Acquired samples','Location','best');
else
    legend('Raster path','Acquired samples','Location','best');
end

hold off;

figure('Name','SBM Acquisition Vectors - Full Timeline', ...
    'Units','normalized', ...
    'Position',[0.05 0.60 0.82 0.34]);

yyaxis left;
plot(t_global_ms, GalvoX, '-b', 'LineWidth',1.1); hold on;
plot(t_global_ms, GalvoY, '-r', 'LineWidth',1.1);
ylabel('Galvo position [mm]');

yyaxis right;
plot(t_global_ms, Exc_Trigger, '-k', 'LineWidth',1.5); hold on;
plot(t_global_ms, Acq_Trigger, '-g', 'LineWidth',1.5);
ylabel('Trigger state');

xlabel('Time [ms]');
title('SBM Acquisition Vectors - Full Timeline');

legend('GalvoX','GalvoY','Excitation Trigger','Acquisition Trigger', ...
    'Location','best');

grid on;
box on;

end