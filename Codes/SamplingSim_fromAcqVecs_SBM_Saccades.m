function [SampledFrames_raw, SampledFrames_vis, SaccadeInfo] = ...
    SamplingSim_fromAcqVecs_SBM_Saccades( ...
    WavePropVid, x_grid, y_grid, t_grid, ...
    GalvoX, GalvoY, AcqTrigger, Exc_Aux_Trigger, ...
    nFrames, nAscans, nBscans, nTurnSamples, T_exc, Ts_OCT, ...
    SaccadeParams, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 25-06-2026
% Date of last modification (DD/MM/YY): 25-06-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SamplingSim_fromAcqVecs_SBM_Saccades
%
% Samples a harmonic wave-propagation movie using SBM acquisition vectors
% while applying saccadic lateral motion lazily at the acquisition query
% points.
%
% This version also upsamples both sampled outputs:
%
%   SampledFrames_raw
%   SampledFrames_vis
%
% to a square grid of size:
%
%   Nout x Nout x nFrames
%
% where:
%
%   Nout = max([nBscans, nAscans, MinimumOutputSize])
%
% By default:
%
%   MinimumOutputSize = 100
%
% This avoids phase-recovery problems caused by very sparse SBM spatial
% sampling, while keeping the output convention:
%
%   SampledFrames_raw -> no serpentine compensation
%   SampledFrames_vis -> serpentine-compensated visualization frames
%
% Saccades are NOT applied by creating a full motion-corrupted video.
% Instead, for each acquired sample:
%
%   x_query = x_nominal + DeltaX(t)
%   y_query = y_nominal + DeltaY(t)
%   t_query = harmonic phase wrapped into the input wave movie
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Parse options =============================

p = inputParser;

addParameter(p, 'InterpMethod', 'cubic');
addParameter(p, 'ExtrapMethod', 'linear');
addParameter(p, 'LoopPeriod_s', []);
addParameter(p, 'ClipQueriesToGrid', false);
addParameter(p, 'MakeDiagnostics', true);

% Local upsampling options.
addParameter(p, 'UpsampleOutputs', true);
addParameter(p, 'MinimumOutputSize', 100);
addParameter(p, 'UpsampleMethod', 'linear');
addParameter(p, 'ApplyUpsampleSmoothing', false);
addParameter(p, 'UpsampleSigma_px', 0.8);

parse(p, varargin{:});
Opt = p.Results;

%% ============================ Input checks ===============================

GalvoX = GalvoX(:);
GalvoY = GalvoY(:);

AcqTrigger = double(AcqTrigger(:));
Exc_Aux_Trigger = double(Exc_Aux_Trigger(:));

N_total = numel(GalvoX);

if numel(GalvoY) ~= N_total || ...
   numel(AcqTrigger) ~= N_total || ...
   numel(Exc_Aux_Trigger) ~= N_total
    error('GalvoX, GalvoY, AcqTrigger, and Exc_Aux_Trigger must have the same length.');
end

x_grid = double(x_grid(:));
y_grid = double(y_grid(:));
t_grid = double(t_grid(:));

WavePropVid = double(WavePropVid);

if size(WavePropVid,1) ~= numel(x_grid)
    error('size(WavePropVid,1) must match numel(x_grid).');
end

if size(WavePropVid,2) ~= numel(y_grid)
    error('size(WavePropVid,2) must match numel(y_grid).');
end

if size(WavePropVid,3) ~= numel(t_grid)
    error('size(WavePropVid,3) must match numel(t_grid).');
end

dt_sample = Ts_OCT;

%% ============================ Unit convention ============================

% If the simulation grid is likely in meters, convert galvo coordinates and
% saccade offsets from mm to m. Otherwise, keep everything in mm.
gridLikelyMeters = max(abs(x_grid)) < 1 && max(abs(y_grid)) < 1;

if gridLikelyMeters
    unitScale = 1e-3;
else
    unitScale = 1;
end

%% ===================== Saccadic trajectory over full timeline ============

P = local_parseSaccadeParams(SaccadeParams);

if ~isempty(P.rngSeed)
    rng(P.rngSeed);
end

[offsetXY_full, stateCode, eventCode] = ...
    local_generateMotionTrajectory(N_total, P, unitScale);

dx_full = offsetXY_full(:,1);
dy_full = offsetXY_full(:,2);

t_global = (0:(N_total-1)).' * dt_sample;

%% ===================== Periodic harmonic interpolant =====================

% Sort all interpolant axes in ascending order because griddedInterpolant
% requires strictly ascending sample points.
[x_grid, ixSort] = sort(x_grid(:), 'ascend');
[y_grid, iySort] = sort(y_grid(:), 'ascend');
[t_grid, itSort] = sort(t_grid(:), 'ascend');

WavePropVid = WavePropVid(ixSort, iySort, itSort);

dt_wave = median(diff(t_grid));

if isempty(Opt.LoopPeriod_s)
    LoopPeriod_s = t_grid(end) - t_grid(1) + dt_wave;
else
    LoopPeriod_s = Opt.LoopPeriod_s;
end

% Make sure the appended periodic endpoint is strictly larger than the
% original movie endpoint.
t_loop_end = t_grid(1) + LoopPeriod_s;

if t_loop_end <= t_grid(end)
    warning(['LoopPeriod_s is shorter than the provided t_grid duration. ', ...
             'Using the full t_grid duration as the harmonic loop period instead.']);

    LoopPeriod_s = t_grid(end) - t_grid(1) + dt_wave;
    t_loop_end = t_grid(1) + LoopPeriod_s;
end

% Append the first frame again at the loop endpoint to enforce periodicity.
t_periodic = [t_grid; t_loop_end];

WavePropVid_periodic = cat(3, WavePropVid, WavePropVid(:,:,1));

F_wave = griddedInterpolant( ...
    {x_grid, y_grid, t_periodic}, ...
    WavePropVid_periodic, ...
    Opt.InterpMethod, ...
    Opt.ExtrapMethod);

%% ===================== Low-resolution frame preallocation ================

SampledFrames_raw_low = nan(nBscans, nAscans, nFrames);
SampledFrames_vis_low = nan(nBscans, nAscans, nFrames);

% Excitation rising edges define the start of each SBM frame.
exc_rising = find(diff([0; double(Exc_Aux_Trigger)]) == 1);

% Tolerance used to detect horizontal plateaus in the raster scan.
tol_abs = 1e-9;
tol_rel = 1e-6;
tol = tol_abs + tol_rel * max(1, max(abs(GalvoY)));

%% =============================== Main loop ===============================

for m = 1:nFrames

    %% -------- Determine current SBM frame bounds -------------------------

    if m <= numel(exc_rising)
        start_idx = exc_rising(m);

        if m < numel(exc_rising)
            end_idx = exc_rising(m+1) - 1;
        else
            end_idx = N_total;
        end
    else
        if ~isempty(exc_rising)
            start_idx = exc_rising(end);
            end_idx = N_total;
        else
            start_idx = 1;
            end_idx = N_total;
        end
    end

    start_idx = max(1, start_idx);
    end_idx   = min(N_total, end_idx);

    if end_idx < start_idx
        continue;
    end

    idx_range = (start_idx:end_idx).';

    Xseg = GalvoX(idx_range);
    Yseg = GalvoY(idx_range);

    %% -------- Detect horizontal scan plateaus ----------------------------

    dY = [0; abs(diff(Yseg))];
    is_plateau = dY < tol;

    run_starts = find(diff([0; is_plateau]) == 1);
    run_ends   = find(diff([is_plateau; 0]) == -1);

    if numel(run_starts) ~= numel(run_ends)
        nmin = min(numel(run_starts), numel(run_ends));
        run_starts = run_starts(1:nmin);
        run_ends   = run_ends(1:nmin);
    end

    run_lengths = run_ends - run_starts + 1;

    minLen = max(1, round(nAscans * 0.5));

    good_mask = run_lengths >= minLen;
    good_starts = run_starts(good_mask);
    good_ends   = run_ends(good_mask);

    % Relax threshold if not enough B-scans were detected.
    thr = minLen - 1;

    while numel(good_starts) < nBscans && thr >= 1
        good_mask = run_lengths >= thr;
        good_starts = run_starts(good_mask);
        good_ends   = run_ends(good_mask);
        thr = thr - 1;
    end

    % Fallback: take the longest runs.
    if numel(good_starts) < nBscans
        [~, idx_sorted] = sort(run_lengths, 'descend');

        pick = idx_sorted(1:min(nBscans, numel(idx_sorted)));
        pick = sort(pick);

        good_starts = run_starts(pick);
        good_ends   = run_ends(pick);
    end

    nTaken = min(nBscans, numel(good_starts));

    good_starts = good_starts(1:nTaken);
    good_ends   = good_ends(1:nTaken);

    %% -------- Build nominal SBM query coordinates ------------------------

    N_needed = nAscans * nBscans;

    Xq_all = nan(N_needed,1);
    Yq_all = nan(N_needed,1);
    Tq_abs_all = nan(N_needed,1);
    GlobalIdx_all = nan(N_needed,1);

    writePos = 1;

    for row = 1:nBscans

        if row <= nTaken

            rs = good_starts(row);
            re = good_ends(row);

            L = re - rs + 1;

            if L >= nAscans

                if nAscans == 1
                    rel_pos = rs;
                else
                    rel_pos = floor((0:(nAscans-1)) * (L-1) / max(1,(nAscans-1))) + 1;
                    rel_pos = rs - 1 + rel_pos;
                end

                global_idxs = idx_range(rel_pos(:));

            else

                % If the detected plateau is shorter than nAscans, expand
                % around the plateau inside the frame.
                avail = (rs:re).';
                global_idxs = idx_range(avail);

                left  = rs - 1;
                right = re + 1;

                left_bound  = 1;
                right_bound = length(idx_range);

                while numel(global_idxs) < nAscans && ...
                        (left >= left_bound || right <= right_bound)

                    if left >= left_bound
                        global_idxs = [idx_range(left); global_idxs]; %#ok<AGROW>
                        left = left - 1;

                        if numel(global_idxs) == nAscans
                            break;
                        end
                    end

                    if right <= right_bound
                        global_idxs = [global_idxs; idx_range(right)]; %#ok<AGROW>
                        right = right + 1;
                    end
                end

                if numel(global_idxs) < nAscans
                    global_idxs = [global_idxs; nan(nAscans - numel(global_idxs),1)];
                end

                global_idxs = global_idxs(1:nAscans);
            end

            idx_write = writePos:(writePos+nAscans-1);

            valid_g = isfinite(global_idxs);

            tempX = nan(nAscans,1);
            tempY = nan(nAscans,1);
            tempT = nan(nAscans,1);

            tempX(valid_g) = GalvoX(global_idxs(valid_g)) * unitScale;
            tempY(valid_g) = GalvoY(global_idxs(valid_g)) * unitScale;
            tempT(valid_g) = (global_idxs(valid_g) - 1) * dt_sample;

            Xq_all(idx_write) = tempX;
            Yq_all(idx_write) = tempY;
            Tq_abs_all(idx_write) = tempT;
            GlobalIdx_all(idx_write) = global_idxs;
        end

        writePos = writePos + nAscans;
    end

    %% -------- Harmonic phase wrapping ------------------------------------

    t_exc_start = (start_idx - 1) * dt_sample;

    t_rel = Tq_abs_all - t_exc_start;

    tq_phase = mod(t_rel, LoopPeriod_s) + t_grid(1);

    zero_idx = (tq_phase == t_grid(1)) & isfinite(Tq_abs_all) & Tq_abs_all > 0;
    tq_phase(zero_idx) = t_grid(1) + eps;

    tq_phase = min(max(tq_phase, min(t_periodic)), max(t_periodic));

    %% -------- Lazy saccadic query-point shifts ---------------------------

    Xq_sacc = Xq_all;
    Yq_sacc = Yq_all;

    valid_gidx = isfinite(GlobalIdx_all);
    gidx = GlobalIdx_all(valid_gidx);

    Xq_sacc(valid_gidx) = Xq_all(valid_gidx) + dx_full(gidx);
    Yq_sacc(valid_gidx) = Yq_all(valid_gidx) + dy_full(gidx);

    if Opt.ClipQueriesToGrid
        Xq_sacc = min(max(Xq_sacc, min(x_grid)), max(x_grid));
        Yq_sacc = min(max(Yq_sacc, min(y_grid)), max(y_grid));
    end

    %% -------- Evaluate wave field ----------------------------------------
    % Preserve the same query convention used in the previous SBM sampler:
    %
    %   F_wave(Y_query, X_query, t_query)

    valid_q = isfinite(Xq_sacc) & isfinite(Yq_sacc) & isfinite(tq_phase);

    sampled_vals = nan(size(Xq_all));

    if any(valid_q)
        sampled_vals(valid_q) = F_wave( ...
            Yq_sacc(valid_q), ...
            Xq_sacc(valid_q), ...
            tq_phase(valid_q));
    end

    frame_img = reshape(sampled_vals, nAscans, nBscans).';

    SampledFrames_raw_low(:,:,m) = frame_img;

    %% -------- Serpentine-compensated visual frame -------------------------

    frame_vis = frame_img;

    for row = 1:nBscans

        if row <= nTaken

            rs = good_starts(row);
            re = good_ends(row);

            rel_idxs = round(linspace(rs, re, nAscans));
            global_idxs = idx_range(rel_idxs);

            X_row = GalvoX(global_idxs);

            if mean(diff(X_row), 'omitnan') < 0
                frame_vis(row,:) = fliplr(frame_vis(row,:));
            end
        end
    end

    SampledFrames_vis_low(:,:,m) = frame_vis;
end

%% =========================== Upsample outputs ============================

if Opt.UpsampleOutputs

    [SampledFrames_raw, UpsampleInfo_raw] = ...
        local_upsample_sbm_movie( ...
            SampledFrames_raw_low, ...
            Opt.MinimumOutputSize, ...
            Opt.UpsampleMethod, ...
            Opt.ApplyUpsampleSmoothing, ...
            Opt.UpsampleSigma_px);

    [SampledFrames_vis, UpsampleInfo_vis] = ...
        local_upsample_sbm_movie( ...
            SampledFrames_vis_low, ...
            Opt.MinimumOutputSize, ...
            Opt.UpsampleMethod, ...
            Opt.ApplyUpsampleSmoothing, ...
            Opt.UpsampleSigma_px);

else
    SampledFrames_raw = SampledFrames_raw_low;
    SampledFrames_vis = SampledFrames_vis_low;

    UpsampleInfo_raw = struct('InterpolationNeeded',false);
    UpsampleInfo_vis = struct('InterpolationNeeded',false);
end

%% ============================== Info output ==============================

SaccadeInfo = struct();

SaccadeInfo.Method = 'Lazy query-point saccades during SBM sampling with local output upsampling';

SaccadeInfo.InputMovieSize = size(WavePropVid);
SaccadeInfo.OutputRawSize = size(SampledFrames_raw);
SaccadeInfo.OutputVisSize = size(SampledFrames_vis);
SaccadeInfo.LowResRawSize = size(SampledFrames_raw_low);
SaccadeInfo.LowResVisSize = size(SampledFrames_vis_low);

SaccadeInfo.Ts_OCT = Ts_OCT;
SaccadeInfo.T_exc = T_exc;
SaccadeInfo.LoopPeriod_s = LoopPeriod_s;

SaccadeInfo.TotalAcquisitionTime_s = N_total * Ts_OCT;
SaccadeInfo.TotalAcquisitionTime_ms = N_total * Ts_OCT * 1e3;

SaccadeInfo.t_global_s = t_global;
SaccadeInfo.t_global_ms = t_global * 1e3;

SaccadeInfo.offsetXY = offsetXY_full;
SaccadeInfo.offsetX_mm = offsetXY_full(:,1) / unitScale;
SaccadeInfo.offsetY_mm = offsetXY_full(:,2) / unitScale;
SaccadeInfo.radialOffset_mm = hypot( ...
    SaccadeInfo.offsetX_mm, ...
    SaccadeInfo.offsetY_mm);

SaccadeInfo.stateCode = stateCode;
SaccadeInfo.eventCode = eventCode;
SaccadeInfo.params = P;
SaccadeInfo.Options = Opt;

SaccadeInfo.Upsampling.Raw = UpsampleInfo_raw;
SaccadeInfo.Upsampling.Visual = UpsampleInfo_vis;

fprintf('\nSamplingSim_fromAcqVecs_SBM_Saccades()\n');
fprintf('Low-res raw frames size   : [%d x %d x %d]\n', size(SampledFrames_raw_low));
fprintf('Low-res vis frames size   : [%d x %d x %d]\n', size(SampledFrames_vis_low));
fprintf('Output raw frames size    : [%d x %d x %d]\n', size(SampledFrames_raw));
fprintf('Output visual frames size : [%d x %d x %d]\n', size(SampledFrames_vis));
fprintf('Total acquisition time    : %.3f ms\n', SaccadeInfo.TotalAcquisitionTime_ms);
fprintf('Max |DeltaX| [mm]         : %.4f\n', max(abs(SaccadeInfo.offsetX_mm)));
fprintf('Max |DeltaY| [mm]         : %.4f\n', max(abs(SaccadeInfo.offsetY_mm)));
fprintf('Max radial offset [mm]    : %.4f\n\n', max(SaccadeInfo.radialOffset_mm));

%% ============================== Diagnostics ==============================

if Opt.MakeDiagnostics

    idx_events = find(eventCode == 1);

    figure('Name','SBM Saccadic Sampling Diagnostics', ...
        'Units','normalized', ...
        'Position',[0.14 0.12 0.70 0.62]);

    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(SaccadeInfo.t_global_ms, SaccadeInfo.offsetX_mm, 'LineWidth',1.2); hold on;
    plot(SaccadeInfo.t_global_ms, SaccadeInfo.offsetY_mm, 'LineWidth',1.2);

    if ~isempty(idx_events)
        xline(SaccadeInfo.t_global_ms(idx_events), '--', ...
            'Color',[0.45 0.45 0.45], ...
            'LineWidth',0.8);
    end

    grid on;
    xlabel('Acquisition time [ms]');
    ylabel('Offset [mm]');
    legend('\Delta x','\Delta y','Location','best');
    title('Saccadic offsets over full SBM acquisition timeline');

    nexttile;
    plot(SaccadeInfo.offsetX_mm, SaccadeInfo.offsetY_mm, ...
        '-k', 'LineWidth',1.2); hold on;

    if ~isempty(idx_events)
        plot(SaccadeInfo.offsetX_mm(idx_events), ...
             SaccadeInfo.offsetY_mm(idx_events), ...
             'or', 'MarkerFaceColor','r', 'MarkerSize',4);
    end

    axis equal;
    grid on;
    xlabel('\Delta x [mm]');
    ylabel('\Delta y [mm]');
    title('Offset trajectory in lateral plane');

    midFrame = round(nFrames/2);

    figure('Name','SBM Low-Resolution vs Upsampled Frames', ...
        'Units','normalized', ...
        'Position',[0.10 0.12 0.78 0.50]);

    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(SampledFrames_vis_low(:,:,midFrame));
    axis image;
    colormap(gca,turbo);
    colorbar;
    title(sprintf('Low-res VIS Frame %d', midFrame));

    nexttile;
    imagesc(SampledFrames_vis(:,:,midFrame));
    axis image;
    colormap(gca,turbo);
    colorbar;
    title(sprintf('Upsampled VIS Frame %d', midFrame));

    nexttile;
    imagesc(SampledFrames_raw_low(:,:,midFrame));
    axis image;
    colormap(gca,turbo);
    colorbar;
    title(sprintf('Low-res RAW Frame %d', midFrame));

    nexttile;
    imagesc(SampledFrames_raw(:,:,midFrame));
    axis image;
    colormap(gca,turbo);
    colorbar;
    title(sprintf('Upsampled RAW Frame %d', midFrame));
end

end

%% ===================================================================== %%
function [MovieOut, Info] = local_upsample_sbm_movie( ...
    MovieIn, MinimumOutputSize, interpMethod, applySmoothing, sigma_px)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local_upsample_sbm_movie
%
% Locally upsamples a sampled SBM movie to a square grid.
%
% Input:
%   MovieIn:
%       [nRows x nCols x nFrames]
%
% Output:
%   MovieOut:
%       [Nout x Nout x nFrames]
%
% where:
%   Nout = max([nRows, nCols, MinimumOutputSize])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nRows, nCols, nT] = size(MovieIn);

Nout = max([nRows, nCols, MinimumOutputSize]);

Info = struct();
Info.InputSize = [nRows nCols nT];
Info.OutputSize = [Nout Nout nT];
Info.MinimumOutputSize = MinimumOutputSize;
Info.InterpMethod = interpMethod;
Info.ApplySmoothing = applySmoothing;
Info.SigmaPx = sigma_px;

if nRows == nCols && nRows >= MinimumOutputSize

    MovieOut = MovieIn;
    Info.InterpolationNeeded = false;

    Info.OriginalRowIndex = (1:nRows).';
    Info.OriginalColIndex = (1:nCols).';
    Info.OutputRowIndex = (1:nRows).';
    Info.OutputColIndex = (1:nCols).';

    return;
end

Info.InterpolationNeeded = true;

row_original = linspace(1, nRows, nRows).';
col_original = linspace(1, nCols, nCols).';

row_output = linspace(1, nRows, Nout).';
col_output = linspace(1, nCols, Nout).';

[Rq, Cq] = ndgrid(row_output, col_output);

MovieOut = nan(Nout, Nout, nT);

for k = 1:nT

    frame_k = MovieIn(:,:,k);

    F = griddedInterpolant( ...
        {row_original, col_original}, ...
        frame_k, ...
        interpMethod, ...
        'none');

    frame_interp = F(Rq, Cq);

    if applySmoothing
        valid_k = isfinite(frame_interp);

        if nnz(valid_k) >= 10
            frame_zero = frame_interp;
            frame_zero(~valid_k) = 0;

            weight = double(valid_k);

            frame_blur  = imgaussfilt(frame_zero, sigma_px);
            weight_blur = imgaussfilt(weight, sigma_px);

            frame_smooth = nan(size(frame_interp));

            valid_smooth = weight_blur > 1e-8;

            frame_smooth(valid_smooth) = ...
                frame_blur(valid_smooth) ./ weight_blur(valid_smooth);

            frame_interp = frame_smooth;
        end
    end

    MovieOut(:,:,k) = frame_interp;
end

Info.OriginalRowIndex = row_original;
Info.OriginalColIndex = col_original;
Info.OutputRowIndex = row_output;
Info.OutputColIndex = col_output;

end

%% ===================================================================== %%
function P = local_parseSaccadeParams(userP)

P = struct();

P.enable = true;

P.fixationMinPts = 8;
P.fixationMaxPts = 20;

P.driftSigma_mm = 0.002;

P.ampDistribution = 'gamma';
P.gammaShape = 2.0;
P.gammaScale_mm = 0.04;

P.lognMu_mm = log(0.06);
P.lognSigma = 0.55;

P.thetaMode = 'uniform';
P.horizontalBias = 0.7;

P.maxOffset_mm = 0.30;

P.rngSeed = [];

if nargin < 1 || isempty(userP)
    return;
end

if ~isstruct(userP)
    error('SaccadeParams must be a struct.');
end

fn = fieldnames(userP);

for k = 1:numel(fn)
    P.(fn{k}) = userP.(fn{k});
end

end

%% ===================================================================== %%
function [offsetXY, stateCode, eventCode] = ...
    local_generateMotionTrajectory(nPoints, P, unitScale)

offsetXY  = zeros(nPoints,2);
stateCode = zeros(nPoints,1);
eventCode = zeros(nPoints,1);

if ~P.enable
    return;
end

currOffset = [0 0];

driftSigma = P.driftSigma_mm * unitScale;
maxOffset  = P.maxOffset_mm  * unitScale;

remainingFixPts = local_sampleFixationDuration(P);

for p = 1:nPoints

    if p == 1
        offsetXY(p,:) = currOffset;
        continue;
    end

    if remainingFixPts <= 0

        amp = local_sampleSaccadeAmplitude(P) * unitScale;
        theta = local_sampleSaccadeAngle(P);

        currOffset = currOffset + amp * [cos(theta), sin(theta)];

        r = hypot(currOffset(1), currOffset(2));

        if r > maxOffset
            currOffset = currOffset * maxOffset / r;
        end

        stateCode(p) = 1;
        eventCode(p) = 1;

        remainingFixPts = local_sampleFixationDuration(P);

    else

        currOffset = currOffset + driftSigma * randn(1,2);

        r = hypot(currOffset(1), currOffset(2));

        if r > maxOffset
            currOffset = currOffset * maxOffset / r;
        end

        remainingFixPts = remainingFixPts - 1;
    end

    offsetXY(p,:) = currOffset;
end

end

%% ===================================================================== %%
function nFix = local_sampleFixationDuration(P)

if P.fixationMinPts == P.fixationMaxPts
    nFix = P.fixationMinPts;
else
    nFix = randi([P.fixationMinPts, P.fixationMaxPts],1,1);
end

end

%% ===================================================================== %%
function amp = local_sampleSaccadeAmplitude(P)

switch lower(P.ampDistribution)

    case 'gamma'
        amp = gamrnd(P.gammaShape, P.gammaScale_mm);

    case 'lognormal'
        amp = lognrnd(P.lognMu_mm, P.lognSigma);

    otherwise
        error('Unsupported ampDistribution.');
end

end

%% ===================================================================== %%
function theta = local_sampleSaccadeAngle(P)

switch lower(P.thetaMode)

    case 'uniform'
        theta = 2*pi*rand();

    case 'horizontal_bias'

        if rand() < P.horizontalBias
            if rand() < 0.5
                theta = deg2rad(20)*randn();
            else
                theta = pi + deg2rad(20)*randn();
            end
        else
            if rand() < 0.5
                theta = pi/2 + deg2rad(20)*randn();
            else
                theta = 3*pi/2 + deg2rad(20)*randn();
            end
        end

        theta = mod(theta, 2*pi);

    otherwise
        error('Unsupported thetaMode.');
end

end