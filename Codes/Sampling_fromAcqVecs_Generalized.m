function [Spatial_Frames, Sampled_Frames_2, SamplingInfo] = Sampling_fromAcqVecs_Generalized( ...
    WaveProp_Vid, x_grid, y_grid, t_grid, ...
    GalvoX, GalvoY, Acq_Trigger, Exc_Trigger, ...
    nFrames, nAscans, nBscans, T_excitation, Ts_OCT, ...
    MotionCase, AcqFramework, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 18-05-2026
% Last modified: 21-05-2026
%
% Sampling_fromAcqVecs_Generalized
%
% Generalized parent sampling function for wave-based OCE acquisition
% vectors.
%
% This version preserves the previous sampling/reconstruction logic but adds
% optional saccadic diagnostic figures and motion-trajectory video.
%
% Important MB behavior:
%
%   The main sampling loop interpolates WaveProp_Vid in the original
%   simulation grid using arbitrary x,y,t query coordinates.
%
%   However, MB spatial frames are reconstructed back onto the reduced
%   nAscans x nAscans MB grid, not the original high-resolution
%   wave-propagation grid.
%
% Saccadic motion behavior:
%
%   x_query = x_nominal + DeltaX
%   y_query = y_nominal + DeltaY
%
%   The moved coordinates affect only the interpolation query. The final
%   spatial mapping uses the nominal acquisition coordinates.
%
% Output convention:
%
%   MB:
%       Spatial_Frames   -> reduced nAscans x nAscans x M_est sampled movie
%       Sampled_Frames_2 -> Meridian_Frames,
%                           size nAscans x M_est*nFrames x nBscans
%
%   BM/SBM:
%       Spatial_Frames   -> visual/orientation-corrected frames
%       Sampled_Frames_2 -> raw acquisition-order frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Parse options =============================

P = inputParser;

addParameter(P, 'SamplesPerSimEnsemble', []);
addParameter(P, 'InterpMethod', 'linear');
addParameter(P, 'ExtrapMethod', 'none');
addParameter(P, 'TemporalDelay_s', 0);
addParameter(P, 'SerpentineVisualization', true);
addParameter(P, 'MakeDiagnostics', true);
addParameter(P, 'SaccadeParams', struct());

parse(P, varargin{:});
Opt = P.Results;

if isempty(Opt.SamplesPerSimEnsemble)
    error(['SamplesPerSimEnsemble must be provided. ', ...
           'For QH MB-like simulations, use M_repetitions + nSyncSamples.']);
end

AcqFramework = upper(string(AcqFramework));
MotionCase   = lower(string(MotionCase));

if ~ismember(AcqFramework, ["MB","BM","SBM"])
    error('AcqFramework must be MB, BM, or SBM.');
end

if ~ismember(MotionCase, ["ideal","saccadic"])
    error('MotionCase must be ideal or saccadic.');
end

%% ======================= Generate saccadic offsets =======================

N_vec  = numel(GalvoX);
x_grid = double(x_grid(:));
y_grid = double(y_grid(:));

gridLikelyMeters = max(abs(x_grid)) < 1 && max(abs(y_grid)) < 1;

if gridLikelyMeters
    unitScale = 1e-3;
else
    unitScale = 1;
end

if MotionCase == "saccadic"

    SaccadeParams = i_parseSaccadeParams(Opt.SaccadeParams);

    if ~isempty(SaccadeParams.rngSeed)
        rng(SaccadeParams.rngSeed);
    end

    [QueryDeltaX, QueryDeltaY, MotionInfo] = i_generateSaccadeQueryOffsets( ...
        Acq_Trigger, SaccadeParams, unitScale);

else

    SaccadeParams = i_parseSaccadeParams(struct('enable', false));

    QueryDeltaX = zeros(N_vec,1);
    QueryDeltaY = zeros(N_vec,1);

    MotionInfo = struct();
    MotionInfo.enabled   = false;
    MotionInfo.offsetXY  = zeros(N_vec,2);
    MotionInfo.stateCode = zeros(N_vec,1);
    MotionInfo.eventCode = zeros(N_vec,1);
    MotionInfo.params    = SaccadeParams;
    MotionInfo.unitScale = unitScale;
end

%% ======================= Run generalized sampler =========================

[SampleValues, SampleInfo] = Main_SamplingfromAcqvecs_LoopFunction( ...
    WaveProp_Vid, x_grid, y_grid, t_grid, ...
    GalvoX, GalvoY, Acq_Trigger, Exc_Trigger, ...
    Ts_OCT, ...
    'SamplesPerSimEnsemble', Opt.SamplesPerSimEnsemble, ...
    'InterpMethod', Opt.InterpMethod, ...
    'ExtrapMethod', Opt.ExtrapMethod, ...
    'TemporalDelay_s', Opt.TemporalDelay_s, ...
    'QueryDeltaX', QueryDeltaX, ...
    'QueryDeltaY', QueryDeltaY);

%% ======================= Protocol-dependent outputs ======================

switch AcqFramework

    case "MB"

        [Spatial_Frames, Sampled_Frames_2, ProtocolInfo] = ...
            i_reorganize_MB_onReducedGrid( ...
                SampleValues, SampleInfo, x_grid, y_grid, ...
                nFrames, nAscans, nBscans);

    case "BM"

        [Sampled_Frames_2, Spatial_Frames, ProtocolInfo] = ...
            i_reorganize_BM_or_SBM( ...
                SampleValues, nFrames, nAscans, nBscans, ...
                Opt.SerpentineVisualization);

    case "SBM"

        [Sampled_Frames_2, Spatial_Frames, ProtocolInfo] = ...
            i_reorganize_BM_or_SBM( ...
                SampleValues, nFrames, nAscans, nBscans, ...
                Opt.SerpentineVisualization);
end

%% =============================== Info ====================================

SamplingInfo = struct();

SamplingInfo.MotionCase    = MotionCase;
SamplingInfo.AcqFramework  = AcqFramework;
SamplingInfo.T_excitation  = T_excitation;
SamplingInfo.Ts_OCT        = Ts_OCT;
SamplingInfo.SampleValues  = SampleValues;
SamplingInfo.SampleInfo    = SampleInfo;
SamplingInfo.MotionInfo    = MotionInfo;
SamplingInfo.ProtocolInfo  = ProtocolInfo;
SamplingInfo.Options       = Opt;

%% ======================= Saccadic diagnostics ============================

if MotionCase == "saccadic"

    if isfield(SaccadeParams, 'makeFigures') && SaccadeParams.makeFigures
        i_plotSaccadeDiagnostics(SampleInfo, MotionInfo, unitScale);
    end

    if isfield(SaccadeParams, 'makeMotionVideo') && SaccadeParams.makeMotionVideo

        videoPath = i_makeSaccadeMotionVideo( ...
            SampleInfo, MotionInfo, unitScale, SaccadeParams);

        SamplingInfo.MotionInfo.motionVideoPath = videoPath;
        MotionInfo.motionVideoPath = videoPath;
    end
end

%% ============================ Command diagnostics ========================

if Opt.MakeDiagnostics

    fprintf('\n============================================================\n');
    fprintf('Sampling_fromAcqVecs_Generalized_v2\n');
    fprintf('============================================================\n');
    fprintf('Framework                  : %s\n', AcqFramework);
    fprintf('Motion case                : %s\n', MotionCase);
    fprintf('Active acquisition samples : %d\n', numel(SampleValues));
    fprintf('Samples per sim ensemble   : %d\n', SampleInfo.SamplesPerSimEnsemble);
    fprintf('Number of sim ensembles    : %d\n', SampleInfo.NumSimEnsembles);
    fprintf('Ensemble cycle             : %s\n', mat2str(SampleInfo.EnsembleCycle));
    fprintf('Spatial_Frames size        : [%s]\n', num2str(size(Spatial_Frames)));
    fprintf('Sampled_Frames_2 size      : [%s]\n', num2str(size(Sampled_Frames_2)));

    if MotionCase == "saccadic"
        dx_mm = SampleInfo.DeltaX / unitScale;
        dy_mm = SampleInfo.DeltaY / unitScale;
        dr_mm = hypot(dx_mm, dy_mm);

        fprintf('Max |DeltaX| [mm]          : %.4f\n', max(abs(dx_mm)));
        fprintf('Max |DeltaY| [mm]          : %.4f\n', max(abs(dy_mm)));
        fprintf('Max radial offset [mm]     : %.4f\n', max(dr_mm));
        fprintf('Mean radial offset [mm]    : %.4f\n', mean(dr_mm,'omitnan'));
    end

    fprintf('============================================================\n\n');
end

end

%% ========================================================================
% MB REORGANIZATION ON REDUCED nAscans x nAscans GRID
% ========================================================================

function [Spatial_Frames, Meridian_Frames, MBInfo] = i_reorganize_MB_onReducedGrid( ...
    SampleValues, SampleInfo, x_grid, y_grid, nFrames, nAscans, nBscans)

groupIDs     = SampleInfo.AcqGroupID;
uniqueGroups = unique(groupIDs(groupIDs > 0));

nGroups = numel(uniqueGroups);

if nGroups == 0
    error('No active acquisition groups detected.');
end

groupLengths = zeros(nGroups,1);

for g = 1:nGroups
    groupLengths(g) = sum(groupIDs == uniqueGroups(g));
end

M_est = round(median(groupLengths(groupLengths > 0)));

Nx_out = nAscans;
Ny_out = nAscans;

x_out = linspace(min(x_grid), max(x_grid), Nx_out).';
y_out = linspace(min(y_grid), max(y_grid), Ny_out).';

Spatial_sum    = zeros(Nx_out, Ny_out, M_est);
Spatial_count  = zeros(Nx_out, Ny_out);
Spatial_Frames = nan(Nx_out, Ny_out, M_est);

Meridian_Frames = nan(nAscans, M_est*nFrames, nBscans);

nSpatialPointsPerFrame = nAscans * nBscans;

for g = 1:nGroups

    idx = find(groupIDs == uniqueGroups(g));

    trace = nan(M_est,1);
    nKeep = min(M_est, numel(idx));
    trace(1:nKeep) = SampleValues(idx(1:nKeep));

    frameIdx = floor((g - 1) / nSpatialPointsPerFrame) + 1;
    q        = mod(g - 1, nSpatialPointsPerFrame) + 1;

    ascanIdx = mod(q - 1, nAscans) + 1;
    bscanIdx = ceil(q / nAscans);

    if frameIdx <= nFrames && bscanIdx <= nBscans

        tStart = (frameIdx - 1)*M_est + 1;
        tEnd   = frameIdx*M_est;

        Meridian_Frames(ascanIdx, tStart:tEnd, bscanIdx) = trace;
    end

    xNom = mean(SampleInfo.Xq_nominal(idx), 'omitnan');
    yNom = mean(SampleInfo.Yq_nominal(idx), 'omitnan');

    [~, ix] = min(abs(x_out - xNom));
    [~, iy] = min(abs(y_out - yNom));

    validTrace = ~isnan(trace);

    if any(validTrace)

        current = reshape(Spatial_sum(ix, iy, :), [M_est,1]);
        current(validTrace) = current(validTrace) + trace(validTrace);

        Spatial_sum(ix, iy, :) = reshape(current, [1,1,M_est]);
        Spatial_count(ix, iy)  = Spatial_count(ix, iy) + 1;
    end
end

occupiedMask = Spatial_count > 0;

for k = 1:M_est

    tmp = nan(Nx_out, Ny_out);
    sliceSum = Spatial_sum(:,:,k);

    tmp(occupiedMask) = sliceSum(occupiedMask) ./ Spatial_count(occupiedMask);

    Spatial_Frames(:,:,k) = tmp;
end

MBInfo = struct();

MBInfo.M_est                  = M_est;
MBInfo.NumGroups              = nGroups;
MBInfo.GroupLengths           = groupLengths;
MBInfo.nSpatialPointsPerFrame = nSpatialPointsPerFrame;
MBInfo.OccupiedMask           = occupiedMask;
MBInfo.x_out                  = x_out;
MBInfo.y_out                  = y_out;
MBInfo.OutputGrid             = 'reduced_MB_grid';
MBInfo.OutputGridSize         = [Nx_out Ny_out];

end

%% ========================================================================
% BM / SBM REORGANIZATION
% ========================================================================

function [Raw_Frames, Visual_Frames, Info] = i_reorganize_BM_or_SBM( ...
    SampleValues, nFrames, nAscans, nBscans, serpentineVisualization)

nPointsPerFrame = nAscans * nBscans;
nExpected       = nFrames * nPointsPerFrame;

Raw_Frames    = nan(nBscans, nAscans, nFrames);
Visual_Frames = nan(nBscans, nAscans, nFrames);

nUsable = min(numel(SampleValues), nExpected);

for s = 1:nUsable

    frameIdx = floor((s - 1) / nPointsPerFrame) + 1;
    q        = mod(s - 1, nPointsPerFrame) + 1;

    bscanIdx = ceil(q / nAscans);
    ascanIdx = mod(q - 1, nAscans) + 1;

    Raw_Frames(bscanIdx, ascanIdx, frameIdx) = SampleValues(s);
end

Visual_Frames = Raw_Frames;

if serpentineVisualization
    for f = 1:nFrames
        for b = 1:nBscans
            if mod(b,2) == 0
                Visual_Frames(b,:,f) = fliplr(Raw_Frames(b,:,f));
            end
        end
    end
end

Info = struct();
Info.nPointsPerFrame         = nPointsPerFrame;
Info.nExpectedSamples        = nExpected;
Info.nUsableSamples          = nUsable;
Info.SerpentineVisualization = serpentineVisualization;

end

%% ========================================================================
% SACCADE PARAMETER PARSER
% ========================================================================

function P = i_parseSaccadeParams(userP)

P = struct();

P.enable = true;

P.fixationMinPts = 8;
P.fixationMaxPts = 20;

P.driftSigma_mm = 0.002;

P.ampDistribution = 'gamma';
P.gammaShape      = 2.0;
P.gammaScale_mm   = 0.04;
P.lognMu_mm       = log(0.06);
P.lognSigma       = 0.55;

P.thetaMode      = 'uniform';
P.horizontalBias = 0.7;

P.maxOffset_mm = 0.30;
P.motionUpdateMode = 'group';

P.rngSeed = [];

P.makeFigures = false;

P.makeMotionVideo = false;
P.motionVideoNumFrames = 300;
P.videoFolder = '';
P.motionVideoName = 'SaccadeTrajectory_Diagnostics.mp4';
P.motionVideoFrameRate = 20;
P.motionVideoStep = 1;
P.motionVideoTailLength = 40;
P.motionVideoQuality = 95;
P.videoFigureVisible = 'on';
P.closeVideoFigureAfterSave = 'off';

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

if P.fixationMinPts < 1 || P.fixationMaxPts < P.fixationMinPts
    error('Invalid fixationMinPts / fixationMaxPts values.');
end

if ~ismember(lower(P.ampDistribution), {'gamma','lognormal'})
    error('ampDistribution must be ''gamma'' or ''lognormal''.');
end

if ~ismember(lower(P.thetaMode), {'uniform','horizontal_bias'})
    error('thetaMode must be ''uniform'' or ''horizontal_bias''.');
end

if ~ismember(lower(P.motionUpdateMode), {'sample','group'})
    error('motionUpdateMode must be ''sample'' or ''group''.');
end

if ~ismember(lower(char(P.videoFigureVisible)), {'on','off'})
    error('videoFigureVisible must be ''on'' or ''off''.');
end

if ~ismember(lower(char(P.closeVideoFigureAfterSave)), {'on','off'})
    error('closeVideoFigureAfterSave must be ''on'' or ''off''.');
end

end

%% ========================================================================
% SACCADE OFFSET GENERATOR
% ========================================================================

function [QueryDeltaX, QueryDeltaY, MotionInfo] = i_generateSaccadeQueryOffsets( ...
    Acq_Trigger, P, unitScale)

Acq_Trigger = double(Acq_Trigger(:));

N_vec = numel(Acq_Trigger);

QueryDeltaX = zeros(N_vec,1);
QueryDeltaY = zeros(N_vec,1);

stateCode = zeros(N_vec,1);
eventCode = zeros(N_vec,1);

activeIdx = find(Acq_Trigger > 0);

if isempty(activeIdx) || ~P.enable
    MotionInfo.enabled   = false;
    MotionInfo.offsetXY  = zeros(N_vec,2);
    MotionInfo.stateCode = stateCode;
    MotionInfo.eventCode = eventCode;
    MotionInfo.params    = P;
    MotionInfo.unitScale = unitScale;
    return;
end

switch lower(P.motionUpdateMode)

    case 'sample'

        nMotionPts = numel(activeIdx);

        [offsetXY_active, state_active, event_active] = ...
            i_generateMotionTrajectory(nMotionPts, P, unitScale);

        QueryDeltaX(activeIdx) = offsetXY_active(:,1);
        QueryDeltaY(activeIdx) = offsetXY_active(:,2);

        stateCode(activeIdx) = state_active;
        eventCode(activeIdx) = event_active;

    case 'group'

        groupID = zeros(N_vec,1);

        groupCounter = 0;
        insideGroup  = false;

        for k = 1:N_vec

            if Acq_Trigger(k) > 0 && ~insideGroup
                groupCounter = groupCounter + 1;
                insideGroup  = true;
            elseif Acq_Trigger(k) <= 0
                insideGroup = false;
            end

            if Acq_Trigger(k) > 0
                groupID(k) = groupCounter;
            end
        end

        [offsetXY_group, state_group, event_group] = ...
            i_generateMotionTrajectory(groupCounter, P, unitScale);

        for g = 1:groupCounter
            idx = groupID == g;

            QueryDeltaX(idx) = offsetXY_group(g,1);
            QueryDeltaY(idx) = offsetXY_group(g,2);

            stateCode(idx) = state_group(g);
            eventCode(idx) = event_group(g);
        end
end

MotionInfo = struct();

MotionInfo.enabled   = true;
MotionInfo.offsetXY  = [QueryDeltaX, QueryDeltaY];
MotionInfo.stateCode = stateCode;
MotionInfo.eventCode = eventCode;
MotionInfo.params    = P;
MotionInfo.unitScale = unitScale;

end

%% ========================================================================
% MOTION TRAJECTORY
% ========================================================================

function [offsetXY, stateCode, eventCode] = i_generateMotionTrajectory(nPoints, P, unitScale)

offsetXY  = zeros(nPoints,2);
stateCode = zeros(nPoints,1);
eventCode = zeros(nPoints,1);

currOffset = [0 0];

driftSigma = P.driftSigma_mm * unitScale;
maxOffset  = P.maxOffset_mm  * unitScale;

remainingFixPts = i_sampleFixationDuration(P);

for p = 1:nPoints

    if p == 1
        offsetXY(p,:) = currOffset;
        stateCode(p) = 0;
        eventCode(p) = 0;
        continue;
    end

    if remainingFixPts <= 0

        amp   = i_sampleSaccadeAmplitude(P) * unitScale;
        theta = i_sampleSaccadeAngle(P);

        currOffset = currOffset + amp*[cos(theta), sin(theta)];

        r = hypot(currOffset(1), currOffset(2));
        if r > maxOffset
            currOffset = currOffset * maxOffset/r;
        end

        stateCode(p) = 1;
        eventCode(p) = 1;

        remainingFixPts = i_sampleFixationDuration(P);

    else

        currOffset = currOffset + driftSigma*randn(1,2);

        r = hypot(currOffset(1), currOffset(2));
        if r > maxOffset
            currOffset = currOffset * maxOffset/r;
        end

        stateCode(p) = 0;
        eventCode(p) = 0;

        remainingFixPts = remainingFixPts - 1;
    end

    offsetXY(p,:) = currOffset;
end

end

%% ========================================================================
% SACCADE DIAGNOSTIC FIGURES

function i_plotSaccadeDiagnostics(SampleInfo, MotionInfo, unitScale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i_plotSaccadeDiagnostics
%
% Minimal and readable diagnostic visualization for simulated saccadic
% motion trajectories.
%
% This simplified version intentionally follows the cleaner layout used in:
%
%   SamplingSim_fromAcqVecs_MB_QH_Saccades_v3()
%
% Only the two most informative plots are displayed:
%
%   1) Motion offsets versus acquisition order
%   2) 2D lateral-plane offset trajectory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ========================= Extract offsets ==============================

dx_mm = SampleInfo.DeltaX / unitScale;
dy_mm = SampleInfo.DeltaY / unitScale;

eventActive = MotionInfo.eventCode(SampleInfo.ActiveIdx);

%% ===================== Reduce repeated event markers ====================

eventGroup = SampleInfo.AcqGroupID(eventActive == 1);

[~, ia] = unique(eventGroup,'stable');

idx_events = find(eventActive == 1);
idx_events = idx_events(ia);

%% ========================= Display downsampling =========================

nSamples = numel(dx_mm);

maxPlotPts = 12000;

plotStride = max(1, ceil(nSamples / maxPlotPts));

idx_plot = 1:plotStride:nSamples;

%% ============================ Main figure ===============================

figure('Name','Generated Lateral Motion Trajectory', ...
    'Units','normalized', ...
    'Position',[0.18 0.12 0.60 0.58]);

tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

%% ===================== Plot 1: offsets vs time ==========================

%% ===================== Build time axis for active samples =================

t_active_s = SampleInfo.AcqTime_s;
t_active_ms = t_active_s * 1e3;

%% ===================== Plot 1: offsets vs time ==========================

nexttile;

plot(t_active_ms(idx_plot), dx_mm(idx_plot), 'LineWidth',1.2); hold on;
plot(t_active_ms(idx_plot), dy_mm(idx_plot), 'LineWidth',1.2);

maxEventLines = 60;

if ~isempty(idx_events)

    stepEvents = max(1, ceil(numel(idx_events)/maxEventLines));
    idx_events_plot = idx_events(1:stepEvents:end);

    for ii = 1:numel(idx_events_plot)

        xline(t_active_ms(idx_events_plot(ii)), '--', ...
            'Color',[0.45 0.45 0.45], ...
            'LineWidth',0.8);

    end
end

grid on;

xlabel('Acquisition time [ms]');
ylabel('Offset [mm]');

legend('\Delta x','\Delta y','Location','best');

title('Motion offsets versus acquisition time');

%% ===================== Plot 2: lateral trajectory =======================

nexttile;

plot(dx_mm(idx_plot), dy_mm(idx_plot), ...
    '-k', 'LineWidth',1.2);

hold on;

if ~isempty(idx_events)

    plot(dx_mm(idx_events), dy_mm(idx_events), ...
        'or', ...
        'MarkerFaceColor','r', ...
        'MarkerSize',4);

end

grid on;

axis equal;

xlabel('\Delta x [mm]');
ylabel('\Delta y [mm]');

title('Offset trajectory in lateral plane');

end

%% ========================================================================
% SACCADE MOTION VIDEO
% ========================================================================

function videoPath = i_makeSaccadeMotionVideo(SampleInfo, MotionInfo, unitScale, P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i_makeSaccadeMotionVideo
%
% Generates an MP4 diagnostic video matching the static saccadic diagnostic
% figure:
%
%   Tile 1:
%       Motion offsets versus acquisition time [ms]
%
%   Tile 2:
%       Offset trajectory in the lateral plane [mm]
%
% Fixed video timing in this version:
%
%   - Output video frame rate:
%
%         30 FPS
%
%   - Number of rendered frames:
%
%         acquisition_duration_seconds * 120
%
% This means the video is temporally oversampled for smooth visualization,
% but encoded at 30 FPS.
%
% Event markers:
%
%   - Vertical dotted lines in the top tile represent only detected
%     saccadic events.
%   - They are not plotted for every rendered video frame.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx_mm = SampleInfo.DeltaX(:) / unitScale;
dy_mm = SampleInfo.DeltaY(:) / unitScale;

t_ms = SampleInfo.AcqTime_s(:) * 1e3;
t_s  = SampleInfo.AcqTime_s(:);

eventActive = MotionInfo.eventCode(SampleInfo.ActiveIdx);
eventActive = eventActive(:);

nPts = numel(dx_mm);

if numel(dy_mm) ~= nPts || numel(t_ms) ~= nPts || numel(eventActive) ~= nPts
    error('Motion video inputs have inconsistent lengths.');
end

idx_events = find(eventActive == 1);

acqDuration_s = max(t_s) - min(t_s);

if acqDuration_s <= 0
    acqDuration_s = nPts * 1e-3;
end

motionVideoFrameRate = 30;
nVideoFrames = min(max(2, round(acqDuration_s * 120)), nPts);

frameIdx = unique(round(linspace(1, nPts, nVideoFrames)));

if isempty(P.videoFolder)
    videoFolder = pwd;
else
    videoFolder = P.videoFolder;
end

if ~exist(videoFolder, 'dir')
    mkdir(videoFolder);
end

videoPath = fullfile(videoFolder, P.motionVideoName);

xLimTime = [min(t_ms), max(t_ms)];

offsetMin = min([dx_mm; dy_mm], [], 'omitnan');
offsetMax = max([dx_mm; dy_mm], [], 'omitnan');

offsetRange = max(offsetMax - offsetMin, 1e-6);

yLimOffset = [offsetMin - 0.08*offsetRange, ...
              offsetMax + 0.08*offsetRange];

xMin = min(dx_mm, [], 'omitnan');
xMax = max(dx_mm, [], 'omitnan');
yMin = min(dy_mm, [], 'omitnan');
yMax = max(dy_mm, [], 'omitnan');

dxRange = max(xMax - xMin, 1e-6);
dyRange = max(yMax - yMin, 1e-6);

xLimPlane = [xMin - 0.15*dxRange, xMax + 0.15*dxRange];
yLimPlane = [yMin - 0.15*dyRange, yMax + 0.15*dyRange];

figVid = figure('Name','Motion Trajectory Video Generator', ...
    'Color','w', ...
    'Units','normalized', ...
    'Position',[0.05 0.06 0.88 0.82], ...
    'Visible', P.videoFigureVisible);

v = VideoWriter(videoPath, 'MPEG-4');
v.FrameRate = motionVideoFrameRate;
v.Quality   = P.motionVideoQuality;
open(v);

fprintf('Generating motion trajectory video: %s\n', videoPath);
fprintf('Video frame rate: %d FPS\n', motionVideoFrameRate);
fprintf('Acquisition duration: %.3f s\n', acqDuration_s);
fprintf('Rendered frames: %d\n', numel(frameIdx));

for ff = 1:numel(frameIdx)

    k = frameIdx(ff);

    clf(figVid);

    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    hold on;

    plot(t_ms, dx_mm, ...
        'Color',[0.65 0.80 1.00], ...
        'LineWidth',1.0);

    plot(t_ms, dy_mm, ...
        'Color',[1.00 0.75 0.55], ...
        'LineWidth',1.0);

    plot(t_ms(1:k), dx_mm(1:k), ...
        'Color',[0.00 0.00 1.00], ...
        'LineWidth',1.8);

    plot(t_ms(1:k), dy_mm(1:k), ...
        'Color',[0.85 0.33 0.10], ...
        'LineWidth',1.8);

    plot(t_ms(k), dx_mm(k), ...
        'o', ...
        'Color',[0.00 0.00 1.00], ...
        'MarkerFaceColor',[0.00 0.00 1.00], ...
        'MarkerSize',7);

    plot(t_ms(k), dy_mm(k), ...
        'o', ...
        'Color',[0.85 0.33 0.10], ...
        'MarkerFaceColor',[0.85 0.33 0.10], ...
        'MarkerSize',7);

    event_lines_to_show = idx_events(idx_events <= k);

    maxEventLines = 80;

    if ~isempty(event_lines_to_show)
        stepEvents = max(1, ceil(numel(event_lines_to_show)/maxEventLines));
        event_lines_to_show = event_lines_to_show(1:stepEvents:end);

        for ii = 1:numel(event_lines_to_show)
            xline(t_ms(event_lines_to_show(ii)), '--', ...
                'Color',[0.55 0.55 0.55], ...
                'LineWidth',0.7);
        end
    end

    if eventActive(k) == 1
        stateStr = 'SACCADE EVENT';
        stateColor = [0.85 0.10 0.10];
    else
        stateStr = 'FIXATION / DRIFT';
        stateColor = [0.10 0.45 0.10];
    end
    
    text(xLimTime(1) + 0.02*range(xLimTime), ...
         yLimOffset(2) - 0.04*range(yLimOffset), ...
         stateStr, ...
         'Color', stateColor, ...
         'FontWeight','bold', ...
         'FontSize',12, ...
         'VerticalAlignment','top', ...
         'BackgroundColor','w', ...
         'Margin',2);

    grid on;
    xlim(xLimTime);
    ylim(yLimOffset);

    xlabel('Acquisition time [ms]');
    ylabel('Offset [mm]');

    legend('\Delta x full','\Delta y full','\Delta x current','\Delta y current', ...
        'Location','northeast');

    title(sprintf('Motion offsets versus acquisition time (current time = %.1f ms)', ...
        t_ms(k)));

    nexttile;
    hold on;

    plot(dx_mm, dy_mm, ...
        '-', ...
        'Color',[0.80 0.80 0.80], ...
        'LineWidth',1.0);

    plot(dx_mm(1:k), dy_mm(1:k), ...
        '-k', ...
        'LineWidth',1.5);

    if isfield(P,'motionVideoTailLength') && ~isempty(P.motionVideoTailLength)
        tailLength = P.motionVideoTailLength;
    else
        tailLength = 40;
    end

    k0 = max(1, k - tailLength + 1);

    plot(dx_mm(k0:k), dy_mm(k0:k), ...
        '-r', ...
        'LineWidth',2.0);

    plot(dx_mm(k), dy_mm(k), ...
        'or', ...
        'MarkerFaceColor','r', ...
        'MarkerSize',8);

    past_events = idx_events(idx_events <= k);

    if ~isempty(past_events)
        plot(dx_mm(past_events), dy_mm(past_events), ...
            'o', ...
            'Color',[0.85 0.10 0.10], ...
            'MarkerFaceColor',[0.85 0.10 0.10], ...
            'MarkerSize',5);
    end

    grid on;
    axis equal;
    xlim(xLimPlane);
    ylim(yLimPlane);

    xlabel('\Delta x [mm]');
    ylabel('\Delta y [mm]');

    title('Offset trajectory in lateral plane');

    drawnow;

    frame = getframe(figVid);
    writeVideo(v, frame);
end

close(v);

if strcmpi(P.closeVideoFigureAfterSave, 'on')
    close(figVid);
end

fprintf('Motion trajectory video saved as: %s\n', videoPath);

end

%% ========================================================================
% SMALL LOCAL HELPERS
% ========================================================================

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