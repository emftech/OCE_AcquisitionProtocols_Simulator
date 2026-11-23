function [Exc_Trigger,Acq_Trigger,GalvoX,GalvoY,Exc_aux] = AcqVectors_SBM(horizontalMin,horizontalMax,verticalMin,...
                                        verticalMax,nAscans,nBscans,nFrames,...
                                        startDelay_samples,nTurnSamples,AcqPulseWidth,...
                                        ExcPulseWidth,T_excitation,Ts_OCT,phaseStepSamples)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Description: Generate acquisition vectors for SBM scanning protocol.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    horizontalMin (1,1) double
    horizontalMax (1,1) double
    verticalMin   (1,1) double
    verticalMax   (1,1) double
    nAscans       (1,1) double {mustBeInteger,   mustBePositive}
    nBscans       (1,1) double {mustBeInteger,   mustBePositive}
    nFrames       (1,1) double {mustBeInteger,   mustBePositive}
    startDelay_samples (1,1) double {mustBeInteger, mustBeNonnegative}
    nTurnSamples  (1,1) double {mustBeInteger,   mustBeNonnegative}
    AcqPulseWidth (1,1) double {mustBeInteger,   mustBePositive}
    ExcPulseWidth (1,1) double {mustBeInteger,   mustBePositive}
    T_excitation  (1,1) double {mustBePositive}
    Ts_OCT        (1,1) double {mustBePositive}
    phaseStepSamples double = []    % optional; computed if empty
end

% compute default phaseStepSamples if not provided (strict integer check)
if isempty(phaseStepSamples)
    rawPhaseStep = (T_excitation / nFrames) / Ts_OCT;
    if abs(rawPhaseStep - round(rawPhaseStep)) > 1e-10
        error(['Default phaseStepSamples is not integer. Ensure (T_excitation / nFrames) ' ...
               'is an integer multiple of Ts_OCT, or provide phaseStepSamples explicitly. ' ...
               'Computed value = %g.'], rawPhaseStep);
    end
    phaseStepSamples = round(rawPhaseStep);
else
    if ~isscalar(phaseStepSamples) || abs(phaseStepSamples - round(phaseStepSamples)) > 1e-10
        error('phaseStepSamples must be a scalar integer number of OCT samples.');
    end
    phaseStepSamples = round(phaseStepSamples);
end

% Useful derived quantity and prints
samplesPerExcPeriod = T_excitation / Ts_OCT; % raw value
samplesPerExcPeriod = round(samplesPerExcPeriod);  % enforce integer
fprintf('samplesPerExcPeriod = %.0f (T_exc = %.3f ms)\n', samplesPerExcPeriod, 1e3 * T_excitation);
fprintf('phaseStepSamples = %.0f\n', phaseStepSamples);

%% ---------------- Building Raster Scanning Pattern -----------------------
% Normalized coordinates: x in [-1,1], y in [1,-1]
y_norm = linspace(1, -1, nBscans);

% Preallocate (including transition slots)
Nsingle_est = nAscans * nBscans + nTurnSamples * (nBscans - 1);
x_positions = nan(1, Nsingle_est);
y_positions = nan(1, Nsingle_est);
isAline_mask = false(1, Nsingle_est);

block_stride = nAscans + nTurnSamples;   % used both here and later

for jj = 1:nBscans
    if mod(jj,2) == 1
        x_row = linspace(-1, 1, nAscans);
    else
        x_row = linspace(1, -1, nAscans);
    end
    y_row = repmat(y_norm(jj), 1, nAscans);
    idx_start = (jj-1) * block_stride + 1;
    idx_end   = idx_start + nAscans - 1;
    x_positions(idx_start:idx_end) = x_row;
    y_positions(idx_start:idx_end) = y_row;
    isAline_mask(idx_start:idx_end) = true(1, nAscans);
end

% Build smooth transitions between rows using cubic Hermite, starting the
% trajectory immediately (first transition sample != p1)
for jj = 1:(nBscans-1)
    idx_start  = (jj-1) * block_stride + 1;
    idx_end    = idx_start + nAscans - 1;
    next_idx_start = jj * block_stride + 1;
    trans_idx = (idx_end+1) : (next_idx_start-1);  % nTurnSamples long (maybe 0)
    if isempty(trans_idx)
        continue;
    end

    p1 = [ x_positions(idx_end), y_positions(idx_end) ];
    p2 = [ x_positions(next_idx_start), y_positions(next_idx_start) ];

    % determine sweep direction for tangents
    if mod(jj,2) == 1, dir1 = +1; else dir1 = -1; end
    if mod(jj+1,2) == 1, dir2 = +1; else dir2 = -1; end

    % tangent magnitude in normalized units (tweakable)
    tanMag = 0.2;

    t1 = [ dir1 * tanMag, 0 ];
    t2 = [ dir2 * tanMag, 0 ];

    nT = numel(trans_idx);
    s = (1:nT)' ./ (nT + 1);   % values in (0,1) — avoids plateau on p1

    h00 =  2*s.^3 - 3*s.^2 + 1;
    h10 =      s.^3 - 2*s.^2 + s;
    h01 = -2*s.^3 + 3*s.^2;
    h11 =      s.^3 -   s.^2;

    Px = h00 * p1(1) + h10 * t1(1) + h01 * p2(1) + h11 * t2(1);
    Py = h00 * p1(2) + h10 * t1(2) + h01 * p2(2) + h11 * t2(2);

    x_positions(trans_idx) = Px;
    y_positions(trans_idx) = Py;
end

% final fallback linear fill if any NaNs remain
nan_mask = isnan(x_positions) | isnan(y_positions);
if any(nan_mask)
    pattern_idx = 1:length(x_positions);
    vx = ~isnan(x_positions);
    vy = ~isnan(y_positions);
    x_positions(nan_mask) = interp1(pattern_idx(vx), x_positions(vx), pattern_idx(nan_mask), 'linear', 'extrap');
    y_positions(nan_mask) = interp1(pattern_idx(vy), y_positions(vy), pattern_idx(nan_mask), 'linear', 'extrap');
end

% Map normalized -> physical (mm)
X_ScanPattern = ((x_positions + 1)/2) * (horizontalMax - horizontalMin) + horizontalMin;
Y_ScanPattern = ((y_positions + 1)/2) * (verticalMax   - verticalMin)   + verticalMin;

aLine_indices = find(isAline_mask);
N_Single_Frame = numel(x_positions);

%% -------- Building Acquisition Vectors over Multiple Frames (global) ----

Transition_Btw_Frames = 2 * (T_excitation / Ts_OCT);
Transition_Btw_Frames = round(Transition_Btw_Frames);

X_single = X_ScanPattern(:);
Y_single = Y_ScanPattern(:);
Nsingle = numel(X_single);

% row offsets inside single block
block_stride = nAscans + nTurnSamples;
row_start_idx = (0:(nBscans-1)) * block_stride + 1;
line_rel_starts = row_start_idx - row_start_idx(1);

% allocate dynamic arrays progressively
GalvoX = [];
GalvoY = [];
Exc_aux = [];
Acq_aux = [];
global_idx = 0;

holdX = X_single(1);
holdY = Y_single(1);

for m = 1:nFrames
    % % 1) tiny pre-gap (1 sample)
    % GalvoX = [GalvoX; holdX];
    % GalvoY = [GalvoY; holdY];
    % Exc_aux = [Exc_aux; 0];
    % Acq_aux = [Acq_aux; 0];
    % global_idx = global_idx + 1;

    % 2) excitation TTL pulses (raise Exc_aux)
    for k = 1:ExcPulseWidth
        GalvoX = [GalvoX; holdX];
        GalvoY = [GalvoY; holdY];
        Exc_aux = [Exc_aux; 1];
        Acq_aux = [Acq_aux; 0];
        global_idx = global_idx + 1;
    end

    % 3) start delay (hold)
    for k = 1:startDelay_samples
        GalvoX = [GalvoX; holdX];
        GalvoY = [GalvoY; holdY];
        Exc_aux = [Exc_aux; 0];
        Acq_aux = [Acq_aux; 0];
        global_idx = global_idx + 1;
    end

    % ---------------------- CHANGED: align galvo to Acq trigger ----------------
    % compute per-frame phase offset (samples) -> this is how much we delay
    % the galvo movement so the first acquisition pulse coincides with the
    % first scanning sample of the frame (X_single(1)).
    frame_phase_offset = round((m-1) * phaseStepSamples);  % in samples
    if frame_phase_offset > 0
        % insert 'frame_phase_offset' hold samples BEFORE appending X_single
        GalvoX = [GalvoX; repmat(holdX, frame_phase_offset, 1)];
        GalvoY = [GalvoY; repmat(holdY, frame_phase_offset, 1)];
        Exc_aux = [Exc_aux; zeros(frame_phase_offset,1)];
        Acq_aux = [Acq_aux; zeros(frame_phase_offset,1)];
        global_idx = global_idx + frame_phase_offset;
    end

    % 4) append scanning block (single-frame pattern) - starts now at scan_global_start
    scan_global_start = global_idx + 1;
    GalvoX = [GalvoX; X_single];
    GalvoY = [GalvoY; Y_single];
    Exc_aux = [Exc_aux; zeros(Nsingle,1)];
    Acq_aux = [Acq_aux; zeros(Nsingle,1)];
    global_idx = global_idx + Nsingle;
    scan_global_end = global_idx;

    % 5) schedule acquisition pulses at start-of-lines (no extra frame_phase_offset
    %     here because we already padded the galvo/scan block above)
    for jj = 1:nBscans
        acq_idx = scan_global_start + line_rel_starts(jj);  % now correct alignment
        if acq_idx >= scan_global_start && acq_idx <= scan_global_end
            rr = acq_idx : min(acq_idx + AcqPulseWidth - 1, scan_global_end);
            Acq_aux(rr) = 1;
        end
    end
    % -------------------------------------------------------------------------

    % 6) append smooth transition back to first-sample of single-frame if not last frame
    if m < nFrames
        L = Transition_Btw_Frames;
        if L < 1, L = 1; end
        L = round(L);
        if exist('samplesPerExcPeriod','var') && samplesPerExcPeriod > 0
            mult = max(1, round(L / samplesPerExcPeriod));
            L = mult * samplesPerExcPeriod;
        end
        L = round(L);

        % endpoints
        p0x = GalvoX(end);
        p0y = GalvoY(end);
        p1x = X_single(1);
        p1y = Y_single(1);

        % derivatives (per-sample)
        if numel(GalvoX) >= 2
            d0x = GalvoX(end) - GalvoX(end-1);
            d0y = GalvoY(end) - GalvoY(end-1);
        else
            d0x = 0; d0y = 0;
        end
        if numel(X_single) >= 2
            d1x = X_single(2) - X_single(1);
            d1y = Y_single(2) - Y_single(1);
        else
            d1x = 0; d1y = 0;
        end

        m0x = d0x * L;
        m0y = d0y * L;
        m1x = d1x * L;
        m1y = d1y * L;

        s = (1:L)' ./ (L + 1);   % exclude exact endpoints to avoid duplicate sample
        h00 =  2*s.^3 - 3*s.^2 + 1;
        h10 =      s.^3 - 2*s.^2 + s;
        h01 = -2*s.^3 + 3*s.^2;
        h11 =      s.^3 -   s.^2;

        Xtrans = h00 * p0x + h10 * m0x + h01 * p1x + h11 * m1x;
        Ytrans = h00 * p0y + h10 * m0y + h01 * p1y + h11 * m1y;

        % append (L is integer-valued double, zeros(L,1) is valid)
        GalvoX = [GalvoX; Xtrans];
        GalvoY = [GalvoY; Ytrans];
        Exc_aux = [Exc_aux; zeros(L,1)];
        Acq_aux = [Acq_aux; zeros(L,1)];
        global_idx = global_idx + L;
    end
end

% Build final triggers & galvo vectors (keep all per-frame excitation pulses)
Exc_Trigger = double(Exc_aux(:));    % keep excitations for every frame
Acq_Trigger = double(Acq_aux(:));
GalvoX = double(GalvoX(:));
GalvoY = double(GalvoY(:));
Exc_aux = double(Exc_aux(:));  % returned auxiliary excitation

% Print summary
N_total = numel(GalvoX);
fprintf('Built waveform timeline: total samples = %d (nFrames=%d), Nsingle=%d, transition=%d\n', N_total, nFrames, Nsingle, Transition_Btw_Frames);
fprintf('Total acquisition time: %.2f ms', N_total*Ts_OCT*1e3);

%% ----- optional visualization code (keeps your previous plotting logic) -----

% Visualize single-frame scanning pattern + inter-frame transition
figure('Name','Scanning pattern + inter-frame transition','Units','normalized','Position',[0.05 0.05 0.6 0.6]);
hold on;
plot(X_ScanPattern, Y_ScanPattern, '-','Color',[0.85 0.85 0.85]);
plot(X_ScanPattern(~isAline_mask), Y_ScanPattern(~isAline_mask), '.','Color',[0.6 0 0], 'MarkerSize',9);
plot(X_ScanPattern(aLine_indices), Y_ScanPattern(aLine_indices), 'bo','MarkerFaceColor','b','MarkerSize',6);
% Lvis = Transition_Btw_Frames; if Lvis < 1, Lvis = 1; end
% Lvis = round(Lvis);
% tt = (0:(Lvis-1))' / max(1,(Lvis-1));
% svis = 3*tt.^2 - 2*tt.^3;
% x_end = X_single(end); y_end = Y_single(end);
% x_start = X_single(1); y_start = Y_single(1);
% Xtrans_single = x_end + (x_start - x_end) .* svis;
% Ytrans_single = y_end + (y_start - y_end) .* svis;
% plot([x_end; Xtrans_single; x_start], [y_end; Ytrans_single; y_start], '-','Color',[0 0.6 0], 'LineWidth',1.6);
% Plotting transitions:
First_BMBetweenFrames_Samples_idx_low = 1 + ExcPulseWidth + startDelay_samples + nAscans*nBscans + (nBscans-1)*(nTurnSamples);
First_BMBetweenFrames_Samples_idx_high = 1 + ExcPulseWidth + startDelay_samples + nAscans*nBscans + (nBscans-1)*(nTurnSamples) + Transition_Btw_Frames;

plot(GalvoX(First_BMBetweenFrames_Samples_idx_low:First_BMBetweenFrames_Samples_idx_high),...
     GalvoY(First_BMBetweenFrames_Samples_idx_low:First_BMBetweenFrames_Samples_idx_high),...
     '-','Color',[0 0.6 0], 'LineWidth',1.6);

axis equal; grid on; xlabel('X (mm)'); ylabel('Y (mm)');
title({sprintf('Raster Scanning pattern (single frame): total samples=%d, A-lines=%d, transitions=%d', N_Single_Frame, numel(aLine_indices), sum(~isAline_mask))...
    sprintf('Total acquisition time: %.2f ms', N_total*Ts_OCT*1e3);});
if (nTurnSamples ~= 0)
    legend("Pattern","Transition Samples","Acquired Samples","Inter-frame Transition","Location","northeast");
else
    legend("Pattern","Acquired Samples","Inter-frame Transition","Location","northeast");
end
hold off;

%% ----------------- Diagnostics plotting (kept similar to your latest logic) -----------------
% global time axis (ms)
t_global_ms = ((0:(N_total-1))' * Ts_OCT) * 1e3;

% detect frame starts from Exc_aux (rising edge)
exc_rising = find(diff([0; Exc_aux(:)]) == 1);
if isempty(exc_rising)
    warning('No excitation rising edges found in Exc_aux - diagnostics may be limited.');
end

% (1) zoomed view of the first frame region (but plot full traces)
figure('Name','First frame timeline (global plot, zoomed view)','Units','normalized','Position',[0.05 0.6 0.75 0.35]);
yyaxis left;
plot(t_global_ms, GalvoX, '-b', 'LineWidth', 1.2); hold on;
plot(t_global_ms, GalvoY, '-r', 'LineWidth', 1.2);
ylabel('X / Y (mm)');
yyaxis right;
smax_global = max(abs([GalvoX; GalvoY])); if smax_global==0, smax_global = 1; end
plot(t_global_ms, Exc_Trigger * smax_global * 0.6, '-k', 'LineWidth', 2); hold on;
plot(t_global_ms, Acq_Trigger * smax_global * 0.4, '-g', 'LineWidth', 2);
if ~isempty(exc_rising)
    first_exc_idx = exc_rising(1);
else
    first_exc_idx = 1;
end
win_start_idx = first_exc_idx;
win_end_idx = min(N_total, first_exc_idx + Nsingle - 1);
win_start_t = (win_start_idx-1) * Ts_OCT * 1e3;
win_end_t   = (win_end_idx-1)   * Ts_OCT * 1e3;
seg_idxs = max(1,win_start_idx):min(N_total,win_end_idx);
yl = [min([GalvoX(seg_idxs); GalvoY(seg_idxs)]) , max([GalvoX(seg_idxs); GalvoY(seg_idxs)])];
ym = 0.12 * (yl(2)-yl(1) + eps);
ylim(yl + [-1 1]*ym);
xlabel('Time [ms]');
legend('X','Y','Excitation Trigger','Acquisition Trigger','Location','northeast');
title('Acquisition Vectors - Frame 1 (global timeline, zoomed view)');
grid on; xlim([win_start_t, win_end_t]); box on; hold off;

% (2) tiled first K frames — zoomed around the excitation / first-acq region
K = min(4, nFrames);
figure('Name','First K frames timelines (excitation + first acquisition zoom)','Units','normalized','Position',[0.05 0.05 0.85 0.85]);
tiledlayout(K,1,'TileSpacing','compact','Padding','compact');

for m = 1:K
    if m <= numel(exc_rising)
        exc_idx = exc_rising(m);
    else
        exc_idx = 1;
    end

    % window size (samples) around the excitation rising edge
    pre  = min(round(samplesPerExcPeriod), exc_idx-1);  % at most one excitation period before
    post = min(round(2*samplesPerExcPeriod + AcqPulseWidth + startDelay_samples + 10), N_total - exc_idx); % a bit after
    win_start_idx = max(1, exc_idx - pre);
    win_end_idx   = min(N_total, exc_idx + post);

    win_start_t = (win_start_idx-1) * Ts_OCT * 1e3;  % ms
    win_end_t   = (win_end_idx-1)   * Ts_OCT * 1e3;  % ms
    seg_idxs = win_start_idx:win_end_idx;

    nexttile;
    yyaxis left;
    % plot full global galvo traces (so interactive pan/zoom still contains everything),
    % then set xlim to the small window below
    plot(t_global_ms, GalvoX, '-b', 'LineWidth', 1.0); hold on;
    plot(t_global_ms, GalvoY, '-r', 'LineWidth', 1.0);

    yyaxis right;
    % scale triggers to the local amplitude so they overlay sensibly
    smax_local = max(abs([GalvoX(seg_idxs); GalvoY(seg_idxs)])); if smax_local==0, smax_local = 1; end
    plot(t_global_ms, Exc_Trigger * smax_local * 0.6, '-k', 'LineWidth', 1.5);
    plot(t_global_ms, Acq_Trigger * smax_local * 0.4, '-g', 'LineWidth', 1.5);

    % vertical markers: excitation rising edge and first acquisition trigger inside view
    t_exc_ms = (exc_idx-1) * Ts_OCT * 1e3;
    xline(t_exc_ms, '--k', 'LineWidth', 1);   % excitation marker

    % find first acquisition trigger inside the current view (if any)
    acq_local_indices = find(Acq_Trigger(seg_idxs));
    if ~isempty(acq_local_indices)
        first_acq_global_idx = seg_idxs(acq_local_indices(1));
        t_acq_ms = (first_acq_global_idx-1) * Ts_OCT * 1e3;
        xline(t_acq_ms, '--g', 'LineWidth', 1); % acquisition marker
    end

    % set precise viewing window (all data still plotted)
    xlim([win_start_t, win_end_t]);

    % comfortable y-limits around the local segment
    yl_local = [min([GalvoX(seg_idxs); GalvoY(seg_idxs)]), max([GalvoX(seg_idxs); GalvoY(seg_idxs)])];
    ym_local = 0.12 * (yl_local(2)-yl_local(1) + eps);
    yyaxis left;
    ylim(yl_local + [-1 1]*ym_local);

    % labels and legend
    yyaxis left;
    xlabel('Time (ms)'); ylabel('X / Y (mm)');
    yyaxis right; ylabel('Triggers Scale (Arbitrary)');
    legend('X','Y','Excitation Trigger','Acquisition Trigger','Location','northeast');
    title(sprintf('Acquisition Vectors - Frame %d (exc & first-acq zoom)', m));
    grid on; box on; hold off;
end

end