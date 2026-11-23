function [RecoveredFrames_Saccades, SaccadeInfo] = ...
    Add_Saccades_To_Recovered_SBM_Frames( ...
    RecoveredFrames_SBM, x_grid, y_grid, Ts_OCT, phaseStepSamples, ...
    FrameAcquisitionTime_s, SaccadeParams, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 25-06-2026
% Date of last modification (DD/MM/YY): 25-06-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add_Saccades_To_Recovered_SBM_Frames
%
% Adds saccadic lateral motion artifacts to recovered SBM frames using a
% pixel-wise acquisition-time model.
%
% -------------------------------------------------------------------------
% Why this version is different
% -------------------------------------------------------------------------
%
% A recovered SBM frame is not truly acquired simultaneously. Even after
% phase compensation, every spatial location corresponds to a different
% physical acquisition time within the raster scan.
%
% Therefore, a single rigid shift per recovered frame is insufficient.
%
% This function instead builds an acquisition-time map:
%
%       AcqTimeMap(x,y)
%
% that assigns each pixel an approximate physical sampling time within the
% SBM raster frame. Then, for every recovered phase frame k, the saccadic
% offset is evaluated pixel-by-pixel as:
%
%       DeltaX_k(x,y) = DeltaX(t_frame_start(k) + AcqTimeMap(x,y))
%       DeltaY_k(x,y) = DeltaY(t_frame_start(k) + AcqTimeMap(x,y))
%
% Finally, the clean recovered movie is sampled using a 3D interpolant:
%
%       F_clean(x + DeltaX_k(x,y), y + DeltaY_k(x,y), t_phase(k))
%
% This mimics the pixel-wise query-coordinate logic used in the generalized
% MB-mode sampling framework, but applied after ideal SBM recovery.
%
% -------------------------------------------------------------------------
% Main assumptions
% -------------------------------------------------------------------------
%
% 1) RecoveredFrames_SBM is already phase-compensated.
%
% 2) The temporal axis of RecoveredFrames_SBM represents the recovered
%    harmonic phase cycle, with temporal spacing:
%
%       dt_recovered = phaseStepSamples * Ts_OCT
%
% 3) Saccades evolve over the full physical acquisition timeline:
%
%       FrameAcquisitionTime_s = numel(GalvoX) * Ts_OCT
%
% 4) Pixel-wise motion is evaluated using the physical SBM raster timing,
%    but the recovered wave phase is kept locked to the recovered frame
%    index. This avoids reintroducing the galvanometer phase artifact that
%    was already compensated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Parse options =============================

p = inputParser;

% Interpolation used to evaluate the recovered clean movie.
addParameter(p, 'InterpMethod', 'linear');

% Extrapolation outside the spatial support.
% Recommended:
%   'none'   -> outside-grid shifted pixels become NaN
%   'nearest' or 'linear' may be used for debugging
addParameter(p, 'ExtrapMethod', 'none');

% Original SBM acquisition geometry.
% IMPORTANT:
%   If RecoveredFrames_SBM was upsampled to 100 x 100, these should still be
%   the original acquisition values, e.g. nBscans_Acq = 21, nAscans_Acq = 37.
addParameter(p, 'nAscans_Acq', []);
addParameter(p, 'nBscans_Acq', []);
addParameter(p, 'nTurnSamples', 0);

% Frame convention:
%   'visual' -> spatial rows are serpentine-corrected for physical display.
%               Even rows were acquired right-to-left, so their timing is
%               reversed spatially.
%
%   'raw'    -> columns are in acquisition order.
addParameter(p, 'FrameConvention', 'visual');

% Which rows are reverse-scanned in the visual frame convention.
% In the SBM acquisition vector function, row 1 is left-to-right and row 2
% is right-to-left, so the default is true.
addParameter(p, 'ReverseEvenRows', true);

% Optional explicit frame-start times. If empty, uniform frame starts are
% estimated from FrameAcquisitionTime_s and nT.
addParameter(p, 'FrameStartTimes_s', []);

% Motion interpolation along the full acquisition timeline.
addParameter(p, 'MotionInterpMethod', 'linear');

% Query clipping:
%   false -> shifted points outside the grid become NaN if ExtrapMethod='none'
%   true  -> shifted points are clipped to the spatial grid boundary
addParameter(p, 'ClipQueriesToGrid', false);

% If true, also add the intra-frame acquisition time to the recovered phase.
% Usually keep false, because recovered SBM frames are already phase
% compensated.
addParameter(p, 'UsePixelAcquisitionPhase', false);

% Sign convention for motion.
% +1 matches the generalized MB sampling convention:
%
%       x_query = x_nominal + DeltaX
%       y_query = y_nominal + DeltaY
addParameter(p, 'MotionSign', +1);

% Optional circular support mask, useful if you want the SBM display to
% mimic the circular corneal support used in MB visualizations.
addParameter(p, 'ApplyCircularMask', false);
addParameter(p, 'CircularMaskRadiusFrac', 0.98);

% Diagnostics.
addParameter(p, 'MakeDiagnostics', true);

parse(p, varargin{:});
Opt = p.Results;

%% ============================ Basic checks ===============================

RecoveredFrames_SBM = double(RecoveredFrames_SBM);

[nX, nY, nT] = size(RecoveredFrames_SBM);

if nT < 2
    error('RecoveredFrames_SBM must contain at least two temporal frames.');
end

if phaseStepSamples <= 0 || abs(phaseStepSamples - round(phaseStepSamples)) > 1e-10
    error('phaseStepSamples must be a positive integer.');
end

phaseStepSamples = round(phaseStepSamples);

if FrameAcquisitionTime_s <= 0
    error('FrameAcquisitionTime_s must be positive.');
end

x_grid = double(x_grid(:));
y_grid = double(y_grid(:));

if numel(x_grid) < 2 || numel(y_grid) < 2
    error('x_grid and y_grid must contain at least two points each.');
end

FrameConvention = lower(string(Opt.FrameConvention));

if ~ismember(FrameConvention, ["visual","raw"])
    error('FrameConvention must be ''visual'' or ''raw''.');
end

%% ========================== Unit convention ==============================

% If the simulation grid is in meters, convert saccade amplitudes from mm to
% meters. If the grid is already in mm, unitScale = 1.
gridLikelyMeters = max(abs(x_grid)) < 1 && max(abs(y_grid)) < 1;

if gridLikelyMeters
    unitScale = 1e-3;
else
    unitScale = 1;
end

% Define the actual spatial axes of the recovered movie. This is important
% when the recovered movie has been upsampled to 100 x 100.
x_axis = linspace(min(x_grid), max(x_grid), nX).';
y_axis = linspace(min(y_grid), max(y_grid), nY).';

[X0, Y0] = ndgrid(x_axis, y_axis);

%% ===================== Recovered harmonic time axis ======================

% The recovered movie represents the harmonic phase cycle after SBM
% compensation.
dt_recovered = phaseStepSamples * Ts_OCT;

RecoveredPeriod_s = nT * dt_recovered;

t_recovered_s = (0:nT-1).' * dt_recovered;

% Periodic endpoint for temporal interpolation.
t_recovered_periodic = [t_recovered_s; RecoveredPeriod_s];

RecoveredFrames_periodic = cat(3, ...
    RecoveredFrames_SBM, ...
    RecoveredFrames_SBM(:,:,1));

%% ===================== Acquisition timing map ============================

% If acquisition geometry is not explicitly provided, infer it from the
% current movie size. This is acceptable only if the movie has not been
% upsampled. For upsampled movies, pass nAscans_Acq and nBscans_Acq.
if isempty(Opt.nAscans_Acq)
    nAscans_Acq = nY;
else
    nAscans_Acq = Opt.nAscans_Acq;
end

if isempty(Opt.nBscans_Acq)
    nBscans_Acq = nX;
else
    nBscans_Acq = Opt.nBscans_Acq;
end

nTurnSamples = Opt.nTurnSamples;

if nAscans_Acq <= 0 || nBscans_Acq <= 0
    error('nAscans_Acq and nBscans_Acq must be positive.');
end

if nTurnSamples < 0
    error('nTurnSamples must be nonnegative.');
end

nAscans_Acq = round(nAscans_Acq);
nBscans_Acq = round(nBscans_Acq);
nTurnSamples = round(nTurnSamples);

% Single-raster duration excluding inter-frame transition.
N_singleFrameSamples = nAscans_Acq*nBscans_Acq + ...
    (nBscans_Acq - 1)*nTurnSamples;

FrameRasterDuration_s = N_singleFrameSamples * Ts_OCT;

% Build low-resolution acquisition-time map on the original SBM grid.
AcqTimeMap_low_s = local_buildSBMAcqTimeMap( ...
    nBscans_Acq, ...
    nAscans_Acq, ...
    nTurnSamples, ...
    Ts_OCT, ...
    FrameConvention, ...
    Opt.ReverseEvenRows);

% Upsample acquisition-time map to the recovered movie grid.
AcqTimeMap_s = local_resize2D(AcqTimeMap_low_s, nX, nY, 'linear');

%% ===================== Full acquisition-time motion trajectory ===========

P = local_parseSaccadeParams(SaccadeParams);

if ~isempty(P.rngSeed)
    rng(P.rngSeed);
end

nMotionSamples = round(FrameAcquisitionTime_s / Ts_OCT) + 1;

t_motion_s = linspace(0, FrameAcquisitionTime_s, nMotionSamples).';
t_motion_ms = t_motion_s * 1e3;

[offsetXY_motion, stateCode, eventCode] = ...
    local_generateMotionTrajectory(nMotionSamples, P, unitScale);

dx_motion = offsetXY_motion(:,1);
dy_motion = offsetXY_motion(:,2);

Fx_motion = griddedInterpolant( ...
    t_motion_s, ...
    dx_motion, ...
    Opt.MotionInterpMethod, ...
    'nearest');

Fy_motion = griddedInterpolant( ...
    t_motion_s, ...
    dy_motion, ...
    Opt.MotionInterpMethod, ...
    'nearest');

%% ===================== Physical frame-start times ========================

if isempty(Opt.FrameStartTimes_s)

    % Uniform approximation of physical frame starts.
    %
    % This is usually sufficient when only GalvoX/GalvoY-derived total
    % acquisition time is available.
    FrameStartTimes_s = linspace(0, FrameAcquisitionTime_s, nT+1).';
    FrameStartTimes_s = FrameStartTimes_s(1:nT);

else

    FrameStartTimes_s = double(Opt.FrameStartTimes_s(:));

    if numel(FrameStartTimes_s) ~= nT
        error('FrameStartTimes_s must have nT elements.');
    end
end

%% ===================== Clean recovered-movie interpolant =================

F_clean = griddedInterpolant( ...
    {x_axis, y_axis, t_recovered_periodic}, ...
    RecoveredFrames_periodic, ...
    Opt.InterpMethod, ...
    Opt.ExtrapMethod);

%% ===================== Optional circular support mask ====================

if Opt.ApplyCircularMask

    x0 = mean([min(x_axis), max(x_axis)]);
    y0 = mean([min(y_axis), max(y_axis)]);

    rx = (max(x_axis) - min(x_axis)) / 2;
    ry = (max(y_axis) - min(y_axis)) / 2;

    rNorm = sqrt(((X0 - x0)/rx).^2 + ((Y0 - y0)/ry).^2);

    CircularMask = rNorm <= Opt.CircularMaskRadiusFrac;

else

    CircularMask = true(nX, nY);
end

%% ===================== Pixel-wise saccadic resampling ====================

RecoveredFrames_Saccades = nan(nX, nY, nT);

for k = 1:nT

    % Physical acquisition time of every pixel in this recovered frame.
    %
    % This is the key correction: the saccade trajectory is evaluated at a
    % different time for each pixel according to the SBM raster timing.
    Tphys_k = FrameStartTimes_s(k) + AcqTimeMap_s;

    % Keep query times inside the known motion trajectory.
    Tphys_k = min(max(Tphys_k, t_motion_s(1)), t_motion_s(end));

    % Pixel-wise saccadic offsets.
    DeltaX_k = Fx_motion(Tphys_k);
    DeltaY_k = Fy_motion(Tphys_k);

    % Spatially shifted query coordinates.
    Xq = X0 + Opt.MotionSign * DeltaX_k;
    Yq = Y0 + Opt.MotionSign * DeltaY_k;

    if Opt.ClipQueriesToGrid
        Xq = min(max(Xq, min(x_axis)), max(x_axis));
        Yq = min(max(Yq, min(y_axis)), max(y_axis));
    end

    % Recovered phase time.
    %
    % Default:
    %   Phase is locked to the recovered movie frame index because the movie
    %   has already been Doppler/galvo phase compensated.
    %
    % Optional:
    %   UsePixelAcquisitionPhase = true reintroduces the intra-frame time
    %   variation into the recovered phase. This is usually not recommended
    %   after recovery, but it is useful for sensitivity testing.
    if Opt.UsePixelAcquisitionPhase
        Tphase_k = mod(t_recovered_s(k) + AcqTimeMap_s, RecoveredPeriod_s);
    else
        Tphase_k = t_recovered_s(k) * ones(nX, nY);
    end

    frame_k = F_clean(Xq, Yq, Tphase_k);

    frame_k(~CircularMask) = nan;

    RecoveredFrames_Saccades(:,:,k) = frame_k;
end

%% ============================== Info struct ==============================

SaccadeInfo = struct();

SaccadeInfo.Method = ...
    'Pixel-wise SBM raster-time saccades applied to recovered frames';

SaccadeInfo.InputSize = [nX nY nT];
SaccadeInfo.OutputSize = size(RecoveredFrames_Saccades);

SaccadeInfo.Ts_OCT = Ts_OCT;
SaccadeInfo.phaseStepSamples = phaseStepSamples;
SaccadeInfo.dt_recovered_s = dt_recovered;

SaccadeInfo.FrameAcquisitionTime_s = FrameAcquisitionTime_s;
SaccadeInfo.FrameAcquisitionTime_ms = FrameAcquisitionTime_s * 1e3;

SaccadeInfo.RecoveredPeriod_s = RecoveredPeriod_s;
SaccadeInfo.RecoveredPeriod_ms = RecoveredPeriod_s * 1e3;

SaccadeInfo.nAscans_Acq = nAscans_Acq;
SaccadeInfo.nBscans_Acq = nBscans_Acq;
SaccadeInfo.nTurnSamples = nTurnSamples;
SaccadeInfo.N_singleFrameSamples = N_singleFrameSamples;
SaccadeInfo.FrameRasterDuration_s = FrameRasterDuration_s;
SaccadeInfo.FrameRasterDuration_ms = FrameRasterDuration_s * 1e3;

SaccadeInfo.AcqTimeMap_low_s = AcqTimeMap_low_s;
SaccadeInfo.AcqTimeMap_s = AcqTimeMap_s;
SaccadeInfo.AcqTimeMap_ms = AcqTimeMap_s * 1e3;

SaccadeInfo.FrameStartTimes_s = FrameStartTimes_s;
SaccadeInfo.FrameStartTimes_ms = FrameStartTimes_s * 1e3;

SaccadeInfo.t_recovered_s = t_recovered_s;
SaccadeInfo.t_recovered_ms = t_recovered_s * 1e3;

SaccadeInfo.t_motion_s = t_motion_s;
SaccadeInfo.t_motion_ms = t_motion_ms;

SaccadeInfo.offsetXY_motion = offsetXY_motion;
SaccadeInfo.offsetX_motion_mm = offsetXY_motion(:,1) / unitScale;
SaccadeInfo.offsetY_motion_mm = offsetXY_motion(:,2) / unitScale;
SaccadeInfo.radialOffset_motion_mm = hypot( ...
    SaccadeInfo.offsetX_motion_mm, ...
    SaccadeInfo.offsetY_motion_mm);

SaccadeInfo.stateCode = stateCode;
SaccadeInfo.eventCode = eventCode;

SaccadeInfo.params = P;
SaccadeInfo.unitScale = unitScale;
SaccadeInfo.Options = Opt;
SaccadeInfo.CircularMask = CircularMask;

%% ============================ Print summary ==============================

fprintf('\nAdd_Saccades_To_Recovered_SBM_Frames()\n');
fprintf('Input recovered movie size        : [%d x %d x %d]\n', nX, nY, nT);
fprintf('Original acquisition grid         : [%d B-scans x %d A-scans]\n', ...
    nBscans_Acq, nAscans_Acq);
fprintf('Turn samples per B-scan transition: %d\n', nTurnSamples);
fprintf('Single SBM raster duration        : %.4f ms\n', FrameRasterDuration_s*1e3);
fprintf('Full acquisition time             : %.4f ms\n', FrameAcquisitionTime_s*1e3);
fprintf('Recovered-frame dt                : %.4f us\n', dt_recovered*1e6);
fprintf('Recovered harmonic period         : %.4f ms\n', RecoveredPeriod_s*1e3);
fprintf('Motion trajectory samples         : %d\n', nMotionSamples);
fprintf('Max |DeltaX| [mm]                 : %.4f\n', max(abs(SaccadeInfo.offsetX_motion_mm)));
fprintf('Max |DeltaY| [mm]                 : %.4f\n', max(abs(SaccadeInfo.offsetY_motion_mm)));
fprintf('Max radial offset [mm]            : %.4f\n', max(SaccadeInfo.radialOffset_motion_mm));
fprintf('Pixel-wise raster-time saccades applied to recovered SBM frames.\n\n');

%% ============================== Diagnostics ==============================

if Opt.MakeDiagnostics

    idx_events = find(eventCode == 1);

    figure('Name','Recovered SBM Pixel-wise Saccadic Motion Diagnostics', ...
        'Units','normalized', ...
        'Position',[0.12 0.10 0.72 0.68]);

    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(SaccadeInfo.AcqTimeMap_ms);
    axis image;
    colormap(gca, turbo);
    colorbar;
    title('Pixel-wise SBM Acquisition Time Map [ms]');
    xlabel('Y index');
    ylabel('X index');

    nexttile;
    plot(t_motion_ms, SaccadeInfo.offsetX_motion_mm, 'LineWidth',1.2); hold on;
    plot(t_motion_ms, SaccadeInfo.offsetY_motion_mm, 'LineWidth',1.2);

    if ~isempty(idx_events)
        xline(t_motion_ms(idx_events), '--', ...
            'Color',[0.45 0.45 0.45], ...
            'LineWidth',0.7);
    end

    grid on;
    xlabel('Acquisition time [ms]');
    ylabel('Offset [mm]');
    legend('\Delta x','\Delta y','Location','best');
    title('Full-timeline Saccadic Offsets');

    nexttile;
    plot(SaccadeInfo.offsetX_motion_mm, ...
         SaccadeInfo.offsetY_motion_mm, ...
         '-k', 'LineWidth',1.2); hold on;

    if ~isempty(idx_events)
        plot(SaccadeInfo.offsetX_motion_mm(idx_events), ...
             SaccadeInfo.offsetY_motion_mm(idx_events), ...
             'or', ...
             'MarkerFaceColor','r', ...
             'MarkerSize',4);
    end

    axis equal;
    grid on;
    xlabel('\Delta x [mm]');
    ylabel('\Delta y [mm]');
    title('Saccadic Lateral Trajectory');

    midIdx = round(nT/2);

    nexttile;
    imagesc(RecoveredFrames_Saccades(:,:,midIdx));
    axis image;
    colormap(gca, turbo);
    colorbar;
    title(sprintf('Saccade-corrupted Frame %d', midIdx));

    figure('Name','Recovered SBM Before/After Pixel-wise Saccades', ...
        'Units','normalized', ...
        'Position',[0.08 0.12 0.82 0.52]);

    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(RecoveredFrames_SBM(:,:,1));
    axis image;
    colormap(gca,turbo);
    colorbar;
    title('Recovered Frame 1 - Before');

    nexttile;
    imagesc(RecoveredFrames_Saccades(:,:,1));
    axis image;
    colormap(gca,turbo);
    colorbar;
    title('Recovered Frame 1 - After Pixel-wise Saccades');

    nexttile;
    imagesc(RecoveredFrames_SBM(:,:,midIdx));
    axis image;
    colormap(gca,turbo);
    colorbar;
    title(sprintf('Recovered Frame %d - Before', midIdx));

    nexttile;
    imagesc(RecoveredFrames_Saccades(:,:,midIdx));
    axis image;
    colormap(gca,turbo);
    colorbar;
    title(sprintf('Recovered Frame %d - After Pixel-wise Saccades', midIdx));
end

end

%% ===================================================================== %%
function AcqTimeMap_s = local_buildSBMAcqTimeMap( ...
    nBscans, nAscans, nTurnSamples, Ts_OCT, FrameConvention, ReverseEvenRows)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local_buildSBMAcqTimeMap
%
% Builds a low-resolution pixel-wise acquisition-time map for one SBM raster.
%
% Output size:
%
%       [nBscans x nAscans]
%
% Each entry contains the physical acquisition time, relative to the start
% of the SBM raster frame, at which that spatial location was sampled.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AcqTimeMap_samples = nan(nBscans, nAscans);

stride = nAscans + nTurnSamples;

for r = 1:nBscans

    rowStart = 1 + (r-1)*stride;

    offsets_lr = 0:(nAscans-1);

    if FrameConvention == "visual"

        if ReverseEvenRows
            isReverseRow = mod(r,2) == 0;
        else
            isReverseRow = mod(r,2) == 1;
        end

        if isReverseRow
            offsets = fliplr(offsets_lr);
        else
            offsets = offsets_lr;
        end

    else
        offsets = offsets_lr;
    end

    AcqTimeMap_samples(r,:) = rowStart + offsets;
end

% Convert sample indices to seconds.
%
% The first acquired sample occurs at approximately t = 0. Therefore subtract
% one sample from the 1-based index.
AcqTimeMap_s = (AcqTimeMap_samples - 1) * Ts_OCT;

end

%% ===================================================================== %%
function B = local_resize2D(A, nRowsOut, nColsOut, method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local_resize2D
%
% Small interpolation-based 2D resize helper that avoids relying on imresize.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nRowsIn, nColsIn] = size(A);

xIn = linspace(1, nRowsIn, nRowsIn).';
yIn = linspace(1, nColsIn, nColsIn).';

xOut = linspace(1, nRowsIn, nRowsOut).';
yOut = linspace(1, nColsIn, nColsOut).';

[Xq, Yq] = ndgrid(xOut, yOut);

F = griddedInterpolant({xIn, yIn}, double(A), method, 'nearest');

B = F(Xq, Yq);

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
function [offsetXY, stateCode, eventCode] = local_generateMotionTrajectory(nPoints, P, unitScale)

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