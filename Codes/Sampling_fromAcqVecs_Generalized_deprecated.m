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
% Generalized parent sampling function for wave-based OCE acquisition
% vectors.
%
% Important MB behavior:
%
%   The main sampling loop still interpolates WaveProp_Vid in the original
%   simulation grid using arbitrary x,y,t query coordinates.
%
%   However, MB spatial frames are reconstructed back onto the reduced
%   nAscans x nAscans MB grid, not the original 500x500 wave-propagation
%   grid.
%
%   This avoids making the sampled support excessively sparse relative to
%   the full wave-propagation grid before subsequent MB interpolation.
%
% Output convention:
%
%   MB:
%       Spatial_Frames   -> reduced nAscans x nAscans x M_est sampled movie
%       Sampled_Frames_2 -> Meridian_Frames, size nAscans x M_est*nFrames x nBscans
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

    QueryDeltaX = zeros(N_vec,1);
    QueryDeltaY = zeros(N_vec,1);

    MotionInfo = struct();
    MotionInfo.enabled  = false;
    MotionInfo.offsetXY = zeros(N_vec,2);
    MotionInfo.params   = struct();
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

%% ============================ Diagnostics ================================

if Opt.MakeDiagnostics

    fprintf('\n============================================================\n');
    fprintf('Sampling_fromAcqVecs_Generalized\n');
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reorganizes MB acquisition samples on the reduced MB grid.
%
% The sampling itself was performed on the original wave-propagation video,
% but the final sparse MB spatial frames are mapped onto:
%
%       nAscans x nAscans x M_est
%
% This reduces excessive sparsity before the interpolation step used by the
% 2D Phase Gradient estimator.
%
% Saccadic moved coordinates are NOT used for this final mapping.
% Only nominal acquisition coordinates are used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    tmp      = nan(Nx_out, Ny_out);
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
P.gammaShape = 2.0;
P.gammaScale_mm = 0.04;
P.lognMu_mm = log(0.06);
P.lognSigma = 0.55;
P.thetaMode = 'uniform';
P.horizontalBias = 0.7;
P.maxOffset_mm = 0.30;
P.motionUpdateMode = 'group';
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

%% ========================================================================
% SACCADE OFFSET GENERATOR
% ========================================================================

function [QueryDeltaX, QueryDeltaY, MotionInfo] = i_generateSaccadeQueryOffsets( ...
    Acq_Trigger, P, unitScale)

Acq_Trigger = double(Acq_Trigger(:));

N_vec = numel(Acq_Trigger);

QueryDeltaX = zeros(N_vec,1);
QueryDeltaY = zeros(N_vec,1);

activeIdx = find(Acq_Trigger > 0);

if isempty(activeIdx) || ~P.enable
    MotionInfo.enabled  = false;
    MotionInfo.offsetXY = zeros(N_vec,2);
    MotionInfo.params   = P;
    return;
end

switch lower(P.motionUpdateMode)

    case 'sample'

        nMotionPts = numel(activeIdx);
        offsetXY = i_generateMotionTrajectory(nMotionPts, P, unitScale);

        QueryDeltaX(activeIdx) = offsetXY(:,1);
        QueryDeltaY(activeIdx) = offsetXY(:,2);

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

        offsetXY_group = i_generateMotionTrajectory(groupCounter, P, unitScale);

        for g = 1:groupCounter
            idx = groupID == g;
            QueryDeltaX(idx) = offsetXY_group(g,1);
            QueryDeltaY(idx) = offsetXY_group(g,2);
        end

    otherwise
        error('motionUpdateMode must be sample or group.');
end

MotionInfo.enabled  = true;
MotionInfo.offsetXY = [QueryDeltaX QueryDeltaY];
MotionInfo.params   = P;

end

%% ========================================================================
% MOTION TRAJECTORY
% ========================================================================

function offsetXY = i_generateMotionTrajectory(nPoints, P, unitScale)

offsetXY = zeros(nPoints,2);

currOffset = [0 0];

driftSigma = P.driftSigma_mm * unitScale;
maxOffset  = P.maxOffset_mm  * unitScale;

remainingFixPts = randi([P.fixationMinPts P.fixationMaxPts]);

for p = 1:nPoints

    if p == 1
        offsetXY(p,:) = currOffset;
        continue;
    end

    if remainingFixPts <= 0

        switch lower(P.ampDistribution)
            case 'gamma'
                amp = gamrnd(P.gammaShape, P.gammaScale_mm) * unitScale;
            case 'lognormal'
                amp = lognrnd(P.lognMu_mm, P.lognSigma) * unitScale;
            otherwise
                error('Unsupported ampDistribution.');
        end

        switch lower(P.thetaMode)
            case 'uniform'
                theta = 2*pi*rand();

            case 'horizontal_bias'
                if rand() < P.horizontalBias
                    theta = pi*(rand() > 0.5) + deg2rad(20)*randn();
                else
                    theta = pi/2 + pi*(rand() > 0.5) + deg2rad(20)*randn();
                end

            otherwise
                error('Unsupported thetaMode.');
        end

        currOffset = currOffset + amp*[cos(theta) sin(theta)];

        r = hypot(currOffset(1), currOffset(2));
        if r > maxOffset
            currOffset = currOffset * maxOffset/r;
        end

        remainingFixPts = randi([P.fixationMinPts P.fixationMaxPts]);

    else

        currOffset = currOffset + driftSigma*randn(1,2);

        r = hypot(currOffset(1), currOffset(2));
        if r > maxOffset
            currOffset = currOffset * maxOffset/r;
        end

        remainingFixPts = remainingFixPts - 1;
    end

    offsetXY(p,:) = currOffset;
end

end