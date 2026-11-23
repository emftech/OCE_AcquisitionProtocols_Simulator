function [SampledFrames_raw, SampledFrames_vis, MotionInfo] = SamplingSim_fromAcqVecs_SBM_H_Saccades(...
    WavePropVid, x_grid, y_grid, t_grid, ...
    GalvoX, GalvoY, AcqTrigger, Exc_Aux_Trigger, ...
    nFrames, nAscans, nBscans, nTurnSamples, T_exc, Ts_OCT, varargin)
% SamplingSim_fromAcqVecs_SBM_H_Saccades
%
% Sample WavePropVid obeying SBM acquisition vectors while introducing
% lateral motion artifacts via a simplified saccadic model.
%
% Outputs:
%   SampledFrames_raw : [nBscans x nAscans x nFrames]
%   SampledFrames_vis : same size, but rows scanned R->L are flipped
%                       for visualization
%   MotionInfo        : struct with generated motion trajectory and params
%
% Motion model:
%   - Fixation phase: small Gaussian drift
%   - Saccade phase : abrupt jump with heavy-tailed amplitude distribution
%
% IMPORTANT:
%   This function preserves the temporal phase mapping of the original SBM
%   sampler. Only the spatial interpolation queries are perturbed.
%
% Optional input:
%   varargin{1} : SaccadeParams struct, same style as in the MB function
%
% Example:
%   SaccadeParams = struct();
%   SaccadeParams.enable = true;
%   SaccadeParams.fixationMinPts = 20;
%   SaccadeParams.fixationMaxPts = 80;
%   SaccadeParams.driftSigma_mm = 0.002;
%   SaccadeParams.ampDistribution = 'gamma';
%   SaccadeParams.gammaShape = 2.0;
%   SaccadeParams.gammaScale_mm = 0.08;
%   SaccadeParams.thetaMode = 'uniform';
%   SaccadeParams.maxOffset_mm = 0.8;
%   SaccadeParams.interpMethod = 'cubic';
%   SaccadeParams.makeFigures = true;
%   SaccadeParams.rngSeed = 1;
%
% Author: Estéfano Machado Ferrer / ChatGPT support
% Date: 16-04-2026

%% ========================= Input normalization ===========================
GalvoX = double(GalvoX(:));
GalvoY = double(GalvoY(:));
AcqTrigger = double(AcqTrigger(:));
Exc_Aux_Trigger = double(Exc_Aux_Trigger(:));

N_total = numel(GalvoX);
if numel(GalvoY) ~= N_total || numel(AcqTrigger) ~= N_total || numel(Exc_Aux_Trigger) ~= N_total
    error('GalvoX/GalvoY/AcqTrigger/Exc_Aux_Trigger must have same length.');
end

if nargin < 14 || isempty(Ts_OCT)
    error('Ts_OCT (OCT sampling period [s]) is required as the last argument.');
end
dt_sample = Ts_OCT;

x_grid = double(x_grid(:));
y_grid = double(y_grid(:));
t_grid = double(t_grid(:));

%% ======================= Optional motion parameters ======================
SaccadeParams = i_parseSaccadeParams(varargin{:});

if ~isempty(SaccadeParams.rngSeed)
    rng(SaccadeParams.rngSeed);
end

%% ========================== Unit harmonization ===========================
% Assume simulation grid is in meters if its extent is < 1 in absolute value.
simGridInMeters = max(abs(x_grid)) < 1 && max(abs(y_grid)) < 1;

if simGridInMeters
    GalvoX_q = GalvoX * 1e-3;   % mm -> m
    GalvoY_q = GalvoY * 1e-3;   % mm -> m
    unitScale = 1e-3;           % mm to m
else
    GalvoX_q = GalvoX;
    GalvoY_q = GalvoY;
    unitScale = 1;
end

%% ========================= Interpolant definition ========================
% NOTE:
% We keep the same convention as your original function:
%   F_wave = griddedInterpolant({x_grid,y_grid,t_grid}, ...)
% but queried as F_wave(Yq, Xq, tq_phase)
% to preserve your original coordinate behavior.
F_wave = griddedInterpolant({x_grid, y_grid, t_grid}, double(WavePropVid), ...
    SaccadeParams.interpMethod, 'linear');

%% ============================ Preallocation ==============================
SampledFrames_raw = nan(nBscans, nAscans, nFrames);
SampledFrames_vis = nan(nBscans, nAscans, nFrames);

%% ===================== Excitation rising edges ===========================
exc_rising = find(diff([0; double(Exc_Aux_Trigger)]) == 1);

%% ===================== Dynamic tolerance for Y plateau ===================
tol_abs = 1e-9;
tol_rel = 1e-6;
tol = tol_abs + tol_rel * max(1, max(abs(GalvoY)));

%% ================== Total number of nominal sampled points ===============
nPointsPerFrame = nAscans * nBscans;
nTotalNominalPoints = nFrames * nPointsPerFrame;

%% ===================== Build motion trajectory ===========================
if SaccadeParams.enable
    xRange = [min(x_grid), max(x_grid)];
    yRange = [min(y_grid), max(y_grid)];

    [offsetXY_all, stateCode_all, eventCode_all] = i_generateLevel1MotionTrajectory(...
        nTotalNominalPoints, SaccadeParams, unitScale, xRange, yRange);
else
    offsetXY_all  = zeros(nTotalNominalPoints,2);
    stateCode_all = zeros(nTotalNominalPoints,1);
    eventCode_all = zeros(nTotalNominalPoints,1);
end

fprintf('SBM motion artifact enabled: %d\n', SaccadeParams.enable);
fprintf('Generated SBM motion trajectory for %d nominal sampled spatial points.\n', nTotalNominalPoints);

%% ========================== Main frame loop ==============================
for m = 1:nFrames

    % determine global start/end indices for this BM frame
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

    idx_range = (start_idx:end_idx)';
    Xseg = GalvoX_q(idx_range);
    Yseg = GalvoY_q(idx_range);

    % plateau detection (where Y is approximately constant)
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

    % choose good runs
    minLen = max(1, round(nAscans * 0.5));
    good_mask = run_lengths >= minLen;
    good_starts = run_starts(good_mask);
    good_ends   = run_ends(good_mask);

    % relax if not enough runs
    thr = minLen - 1;
    while numel(good_starts) < nBscans && thr >= 1
        good_mask = run_lengths >= thr;
        good_starts = run_starts(good_mask);
        good_ends   = run_ends(good_mask);
        thr = thr - 1;
    end

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

    % prepare queries
    N_needed = nAscans * nBscans;
    Xq_all_nom = nan(N_needed,1);
    Yq_all_nom = nan(N_needed,1);
    Tq_all     = nan(N_needed,1);

    % for VIS flipping
    Xrow_nominal_cell = cell(nBscans,1);

    writePos = 1;

    for row = 1:nBscans
        if row <= nTaken
            rs = good_starts(row);
            re = good_ends(row);

            % robust selection of nAscans samples
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
                avail = (rs:re)';
                global_idxs = idx_range(avail);

                left  = rs - 1;
                right = re + 1;

                left_bound  = 1;
                right_bound = length(idx_range);

                while numel(global_idxs) < nAscans && (left >= left_bound || right <= right_bound)
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

            valid_global = ~isnan(global_idxs);

            Xtmp = nan(nAscans,1);
            Ytmp = nan(nAscans,1);
            Ttmp = nan(nAscans,1);

            Xtmp(valid_global) = GalvoX_q(global_idxs(valid_global));
            Ytmp(valid_global) = GalvoY_q(global_idxs(valid_global));
            Ttmp(valid_global) = (global_idxs(valid_global) - 1) * dt_sample;

            Xq_all_nom(writePos:writePos+nAscans-1) = Xtmp;
            Yq_all_nom(writePos:writePos+nAscans-1) = Ytmp;
            Tq_all(writePos:writePos+nAscans-1)     = Ttmp;

            Xrow_nominal_cell{row} = Xtmp;
        end

        writePos = writePos + nAscans;
    end

    % ====================== Apply motion offsets ==========================
    pointStart = (m-1)*nPointsPerFrame + 1;
    pointEnd   = m*nPointsPerFrame;

    offsets_frame = offsetXY_all(pointStart:pointEnd, :);

    Xq_all_mov = Xq_all_nom + offsets_frame(:,1);
    Yq_all_mov = Yq_all_nom + offsets_frame(:,2);

    % Clip moved queries to simulation bounds
    Xq_all_mov = min(max(Xq_all_mov, min(x_grid)), max(x_grid));
    Yq_all_mov = min(max(Yq_all_mov, min(y_grid)), max(y_grid));

    % ====================== Time mapping (unchanged) ======================
    t_exc_start = (start_idx - 1) * dt_sample;
    t_rel = Tq_all - t_exc_start;
    tq_phase = mod(t_rel, T_exc);

    zero_idx = (tq_phase == 0) & (~isnan(Tq_all) & Tq_all > 0);
    tq_phase(zero_idx) = eps;

    % ======================== Sample interpolant ==========================
    valid_q = ~isnan(Xq_all_mov) & ~isnan(Yq_all_mov) & ~isnan(tq_phase);
    sampled_vals = nan(size(Xq_all_mov));

    if any(valid_q)
        % Preserve original query ordering convention
        sampled_vals(valid_q) = F_wave(Yq_all_mov(valid_q), Xq_all_mov(valid_q), tq_phase(valid_q));
    end

    frame_img = reshape(sampled_vals, nAscans, nBscans).';
    SampledFrames_raw(:,:,m) = frame_img;

    % ===================== VIS version row flipping =======================
    frame_vis = frame_img;
    for row = 1:nBscans
        if row <= nTaken && ~isempty(Xrow_nominal_cell{row})
            X_row = Xrow_nominal_cell{row};
            X_row = X_row(~isnan(X_row));

            if numel(X_row) >= 2
                if mean(diff(X_row)) < 0
                    frame_vis(row,:) = fliplr(frame_vis(row,:));
                end
            end
        end
    end
    SampledFrames_vis(:,:,m) = frame_vis;
end

fprintf('Sampling done: produced SampledFrames_raw (%dx%dx%d) and SampledFrames_vis (%dx%dx%d)\n', ...
    size(SampledFrames_raw,1), size(SampledFrames_raw,2), size(SampledFrames_raw,3), ...
    size(SampledFrames_vis,1), size(SampledFrames_vis,2), size(SampledFrames_vis,3));

%% ============================= MotionInfo ================================
MotionInfo = struct();
MotionInfo.params             = SaccadeParams;
MotionInfo.offsetXY           = offsetXY_all;
MotionInfo.stateCode          = stateCode_all;
MotionInfo.eventCode          = eventCode_all;
MotionInfo.exc_rising         = exc_rising;
MotionInfo.nTotalNominalPoints = nTotalNominalPoints;
MotionInfo.unitScale          = unitScale;
MotionInfo.motionGifFile      = '';

%% ============================== Diagnostics =============================
if SaccadeParams.makeFigures
    offsetX_plot = offsetXY_all(:,1) / unitScale;
    offsetY_plot = offsetXY_all(:,2) / unitScale;
    idx_events = find(eventCode_all == 1);

    figure('Name','SBM Generated Lateral Motion Trajectory', ...
        'Units','normalized','Position',[0.08 0.08 0.78 0.78]);

    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(offsetX_plot, 'LineWidth',1.1); hold on;
    plot(offsetY_plot, 'LineWidth',1.1);
    if ~isempty(idx_events)
        xline(idx_events,'--');
    end
    grid on;
    xlabel('Sampled spatial point index');
    ylabel('Offset [mm]');
    legend('\Delta x','\Delta y','Location','best');
    title('Motion offsets versus acquisition order');

    nexttile;
    plot(offsetX_plot, offsetY_plot, '-k', 'LineWidth',1.1); hold on;
    if ~isempty(idx_events)
        plot(offsetX_plot(idx_events), offsetY_plot(idx_events), 'or', 'MarkerFaceColor','r');
    end
    grid on; axis equal;
    xlabel('\Delta x [mm]');
    ylabel('\Delta y [mm]');
    title('Offset trajectory in lateral plane');
end

%% ======================== Optional GIF generation =======================
if isfield(SaccadeParams,'makeMotionGIF') && SaccadeParams.makeMotionGIF
    try
        gifName    = 'SBM_SaccadeTrajectory_Diagnostics.gif';
        stepGIF    = 1;
        tailLength = 40;
        delayTime  = 0.06;

        offsetX_plot = offsetXY_all(:,1) / unitScale;
        offsetY_plot = offsetXY_all(:,2) / unitScale;

        nPts = numel(offsetX_plot);
        idx_events = find(eventCode_all == 1);

        xMin = min(offsetX_plot); xMax = max(offsetX_plot);
        yMin = min(offsetY_plot); yMax = max(offsetY_plot);

        dx = max(xMax - xMin, 1e-6);
        dy = max(yMax - yMin, 1e-6);

        xLimTrace = [1, nPts];
        yLimTrace = [min([offsetX_plot; offsetY_plot]) - 0.1*max(dx,dy), ...
                     max([offsetX_plot; offsetY_plot]) + 0.1*max(dx,dy)];

        xLimPlane = [xMin - 0.15*dx, xMax + 0.15*dx];
        yLimPlane = [yMin - 0.15*dy, yMax + 0.15*dy];

        gifPath = fullfile(pwd, gifName);

        figGIF = figure( ...
            'Name','SBM Motion Trajectory GIF Generator', ...
            'Color','w', ...
            'Units','pixels', ...
            'Position',[100 100 1100 850], ...
            'Visible','off');

        firstFrame = true;

        for k = 1:stepGIF:nPts
            clf(figGIF);
            tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

            nexttile;
            hold on;
            plot(1:nPts, offsetX_plot, 'Color',[0.65 0.80 1.00], 'LineWidth',1.0);
            plot(1:nPts, offsetY_plot, 'Color',[1.00 0.75 0.55], 'LineWidth',1.0);
            plot(1:k, offsetX_plot(1:k), 'b', 'LineWidth',1.8);
            plot(1:k, offsetY_plot(1:k), 'Color',[0.85 0.33 0.10], 'LineWidth',1.8);
            plot(k, offsetX_plot(k), 'ob', 'MarkerFaceColor','b', 'MarkerSize',7);
            plot(k, offsetY_plot(k), 'o', ...
                'Color',[0.85 0.33 0.10], ...
                'MarkerFaceColor',[0.85 0.33 0.10], ...
                'MarkerSize',7);

            past_events = idx_events(idx_events <= k);
            if ~isempty(past_events)
                for ee = 1:numel(past_events)
                    xline(past_events(ee), '--', 'Color',[0.6 0.6 0.6], 'LineWidth',0.8);
                end
            end

            grid on;
            xlim(xLimTrace);
            ylim(yLimTrace);
            xlabel('Sampled spatial point index');
            ylabel('Offset [mm]');
            title(sprintf('SBM motion offsets versus acquisition order (current index = %d / %d)', k, nPts));
            legend('\Delta x full','\Delta y full','\Delta x current','\Delta y current', ...
                'Location','northeast');

            if eventCode_all(k) == 1
                stateStr = 'SACCADE EVENT';
                stateColor = [0.85 0.10 0.10];
            else
                stateStr = 'FIXATION / DRIFT';
                stateColor = [0.10 0.45 0.10];
            end

            text(0.02*nPts, yLimTrace(2) - 0.08*range(yLimTrace), ...
                stateStr, 'Color', stateColor, 'FontWeight','bold', 'FontSize',11);

            nexttile;
            hold on;
            plot(offsetX_plot, offsetY_plot, '-', 'Color',[0.80 0.80 0.80], 'LineWidth',1.0);
            plot(offsetX_plot(1:k), offsetY_plot(1:k), '-k', 'LineWidth',1.5);

            k0 = max(1, k - tailLength + 1);
            plot(offsetX_plot(k0:k), offsetY_plot(k0:k), '-r', 'LineWidth',2.0);
            plot(offsetX_plot(k), offsetY_plot(k), 'or', 'MarkerFaceColor','r', 'MarkerSize',8);

            if ~isempty(past_events)
                plot(offsetX_plot(past_events), offsetY_plot(past_events), ...
                    'o', 'Color',[0.85 0.10 0.10], ...
                    'MarkerFaceColor',[0.85 0.10 0.10], 'MarkerSize',5);
            end

            grid on;
            axis equal;
            xlim(xLimPlane);
            ylim(yLimPlane);
            xlabel('\Delta x [mm]');
            ylabel('\Delta y [mm]');
            title('Offset trajectory in lateral plane');

            drawnow;

            frame = getframe(figGIF);
            [im, map] = rgb2ind(frame2im(frame), 256);

            if firstFrame
                imwrite(im, map, gifPath, 'gif', 'LoopCount', inf, 'DelayTime', delayTime);
                firstFrame = false;
            else
                imwrite(im, map, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
            end
        end

        close(figGIF);
        MotionInfo.motionGifFile = gifPath;
        fprintf('SBM motion trajectory GIF saved as: %s\n', gifPath);

    catch ME
        warning('Could not generate SBM motion GIF. Reason: %s', ME.message);
        MotionInfo.motionGifFile = '';
    end
end

end

%% ===================================================================== %%
%% ============================ LOCAL FUNCTIONS ========================= %%
%% ===================================================================== %%

function P = i_parseSaccadeParams(varargin)

P = struct();

P.enable = true;

% Duration of fixation segments in number of sampled spatial points
P.fixationMinPts = 20;
P.fixationMaxPts = 80;

% Small drift during fixation
P.driftSigma_mm = 0.002;

% Heavy-tailed amplitude model
P.ampDistribution = 'gamma';   % 'gamma' or 'lognormal'

% Gamma parameters
P.gammaShape    = 2.0;
P.gammaScale_mm = 0.04;

% Lognormal parameters
P.lognMu_mm   = log(0.06);
P.lognSigma   = 0.55;

% Angular distribution
P.thetaMode      = 'uniform';   % 'uniform' or 'horizontal_bias'
P.horizontalBias = 0.7;

% Maximum allowed global offset radius
P.maxOffset_mm = 0.30;

% Interpolation method for wave interpolant
P.interpMethod = 'cubic';   % 'nearest', 'linear', or 'cubic'

% Diagnostics
P.makeFigures   = true;
P.makeMotionGIF = false;

% RNG
P.rngSeed = [];

if nargin < 1 || isempty(varargin{1})
    return;
end

userP = varargin{1};
if ~isstruct(userP)
    error('Optional motion parameter input must be a struct.');
end

fn = fieldnames(userP);
for k = 1:numel(fn)
    P.(fn{k}) = userP.(fn{k});
end

if P.fixationMinPts < 1 || P.fixationMaxPts < P.fixationMinPts
    error('Invalid fixationMinPts / fixationMaxPts values.');
end

if ~ismember(lower(P.ampDistribution), {'gamma','lognormal'})
    error('ampDistribution must be ''gamma'' or ''lognormal''.');
end

if ~ismember(lower(P.thetaMode), {'uniform','horizontal_bias'})
    error('thetaMode must be ''uniform'' or ''horizontal_bias''.');
end

if ~ismember(lower(P.interpMethod), {'nearest','linear','cubic'})
    error('interpMethod must be ''nearest'', ''linear'', or ''cubic''.');
end

end

function [offsetXY, stateCode, eventCode] = i_generateLevel1MotionTrajectory(...
    nPoints, P, unitScale, xRange, yRange)

driftSigma = P.driftSigma_mm * unitScale;
maxOffsetR = P.maxOffset_mm * unitScale;

offsetXY  = zeros(nPoints,2);
stateCode = zeros(nPoints,1);  % 0 fixation/drift, 1 saccade sample
eventCode = zeros(nPoints,1);  % 1 exactly where jump was applied

currOffset = [0, 0];

halfWidthX = 0.5 * abs(diff(xRange));
halfWidthY = 0.5 * abs(diff(yRange));
maxSafeR   = 0.45 * min([halfWidthX, halfWidthY, maxOffsetR]);
maxSafeR   = max(maxSafeR, eps);

remainingFixPts = i_sampleFixationDuration(P);

for p = 1:nPoints

    if p == 1
        stateCode(p) = 0;
        offsetXY(p,:) = currOffset;
        continue;
    end

    if remainingFixPts <= 0
        amp   = i_sampleSaccadeAmplitude(P) * unitScale;
        theta = i_sampleSaccadeAngle(P);

        jump = amp * [cos(theta), sin(theta)];
        currOffset = currOffset + jump;
        currOffset = i_clipOffsetRadius(currOffset, maxSafeR);

        stateCode(p) = 1;
        eventCode(p) = 1;

        remainingFixPts = i_sampleFixationDuration(P);
    else
        drift = driftSigma * randn(1,2);
        currOffset = currOffset + drift;
        currOffset = i_clipOffsetRadius(currOffset, maxSafeR);

        stateCode(p) = 0;
        remainingFixPts = remainingFixPts - 1;
    end

    offsetXY(p,:) = currOffset;
end

end

function nFix = i_sampleFixationDuration(P)
if P.fixationMinPts == P.fixationMaxPts
    nFix = P.fixationMinPts;
else
    nFix = randi([P.fixationMinPts, P.fixationMaxPts],1,1);
end
end

function amp = i_sampleSaccadeAmplitude(P)
switch lower(P.ampDistribution)
    case 'gamma'
        amp = gamrnd(P.gammaShape, P.gammaScale_mm);
    case 'lognormal'
        amp = lognrnd(P.lognMu_mm, P.lognSigma);
    otherwise
        error('Unsupported ampDistribution.');
end
end

function theta = i_sampleSaccadeAngle(P)
switch lower(P.thetaMode)
    case 'uniform'
        theta = 2*pi*rand();

    case 'horizontal_bias'
        if rand() < P.horizontalBias
            if rand() < 0.5
                theta = 0 + deg2rad(20)*randn();
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

function v = i_clipOffsetRadius(v, maxR)
r = hypot(v(1), v(2));
if r > maxR
    v = v .* (maxR / r);
end
end