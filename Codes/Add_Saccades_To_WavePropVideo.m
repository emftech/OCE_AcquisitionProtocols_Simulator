function [WavePropVid_Saccades, SaccadeInfo] = Add_Saccades_To_WavePropVideo( ...
    WavePropVid, x_grid, y_grid, t_grid, SaccadeParams, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation (DD/MM/YY): 25-06-2026
% Date of last modification (DD/MM/YY): 25-06-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add_Saccades_To_WavePropVideo
%
% Adds lateral saccadic motion artifacts to a harmonic wave-propagation
% movie by extending it over the full acquisition timeline.
%
% This version is intended for long SBM acquisitions where the original
% harmonic wave-propagation movie may only contain one or a few excitation
% periods, while the acquisition timeline is much longer.
%
% Main idea:
%
%   1) Treat the input wave-propagation movie as a periodic harmonic movie.
%
%   2) Build a longer temporal movie that spans the full acquisition time:
%
%          TotalAcquisitionTime_s
%
%      For example:
%
%          TotalAcquisitionTime_s = numel(GalvoX) * Ts_OCT
%
%   3) For every time sample in the full acquisition timeline, query the
%      original wave movie at the wrapped harmonic phase:
%
%          t_phase = mod(t_full, T_loop)
%
%   4) Apply saccadic spatial shifts at every full-timeline sample:
%
%          x_query = x_nominal + DeltaX(t_full)
%          y_query = y_nominal + DeltaY(t_full)
%
%   5) Return the motion-corrupted, acquisition-length wave movie:
%
%          WavePropVid_Saccades(:,:,k)
%
%      which can then be passed directly into the regular SBM sampler.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
%   WavePropVid:
%       Original harmonic wave-propagation movie.
%       Size:
%
%           [nX x nY x nT]
%
%   x_grid, y_grid:
%       Spatial simulation axes.
%
%   t_grid:
%       Temporal axis of the original harmonic movie.
%
%   SaccadeParams:
%       Struct defining the saccadic motion model.
%
% Required optional name-value pair:
%
%   'TotalAcquisitionTime_s':
%       Total acquisition duration to synthesize.
%       Recommended:
%
%           TotalAcquisitionTime_s = numel(GalvoX) * Ts_OCT;
%
% Optional name-value pairs:
%
%   'OutputTs_s':
%       Temporal sampling period of the output movie.
%       Recommended:
%
%           OutputTs_s = Ts_OCT;
%
%   'LoopPeriod_s':
%       Period used to wrap the harmonic movie.
%       If empty, it is estimated as:
%
%           t_grid(end) - t_grid(1) + median(diff(t_grid))
%
%   'InterpMethod':
%       Interpolation method for griddedInterpolant.
%
%   'ExtrapMethod':
%       Extrapolation method for griddedInterpolant.
%
%   'ClipQueriesToGrid':
%       If true, clipped shifted spatial queries to the grid boundary.
%
%   'MakeDiagnostics':
%       If true, plot motion diagnostics and before/after frames.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Parse options =============================

p = inputParser;

addParameter(p, 'TotalAcquisitionTime_s', []);
addParameter(p, 'OutputTs_s', []);
addParameter(p, 'LoopPeriod_s', []);

addParameter(p, 'InterpMethod', 'linear');
addParameter(p, 'ExtrapMethod', 'none');

addParameter(p, 'ClipQueriesToGrid', false);
addParameter(p, 'MakeDiagnostics', true);

parse(p, varargin{:});
Opt = p.Results;

if isempty(Opt.TotalAcquisitionTime_s)
    error(['TotalAcquisitionTime_s must be provided. ', ...
           'Use TotalAcquisitionTime_s = numel(GalvoX) * Ts_OCT.']);
end

if Opt.TotalAcquisitionTime_s <= 0
    error('TotalAcquisitionTime_s must be positive.');
end

%% ============================ Basic checks ===============================

WavePropVid = double(WavePropVid);

[nX, nY, nT_original] = size(WavePropVid);

x_grid = double(x_grid(:));
y_grid = double(y_grid(:));
t_grid = double(t_grid(:));

if numel(x_grid) ~= nX
    error('numel(x_grid) must match size(WavePropVid,1).');
end

if numel(y_grid) ~= nY
    error('numel(y_grid) must match size(WavePropVid,2).');
end

if numel(t_grid) ~= nT_original
    error('numel(t_grid) must match size(WavePropVid,3).');
end

if nT_original < 2
    error('WavePropVid must contain at least two temporal samples.');
end

%% ====================== Temporal output configuration ====================

dt_original = median(diff(t_grid));

if isempty(Opt.OutputTs_s)
    OutputTs_s = dt_original;
else
    OutputTs_s = Opt.OutputTs_s;
end

if OutputTs_s <= 0
    error('OutputTs_s must be positive.');
end

if isempty(Opt.LoopPeriod_s)
    LoopPeriod_s = t_grid(end) - t_grid(1) + dt_original;
else
    LoopPeriod_s = Opt.LoopPeriod_s;
end

if LoopPeriod_s <= 0
    error('LoopPeriod_s must be positive.');
end

nT_full = round(Opt.TotalAcquisitionTime_s / OutputTs_s) + 1;

t_full = (0:nT_full-1).' * OutputTs_s;
t_full(t_full > Opt.TotalAcquisitionTime_s) = [];

nT_full = numel(t_full);

% Wrapped harmonic phase used to query the original movie.
t_phase = mod(t_full - t_full(1), LoopPeriod_s) + t_grid(1);

% Avoid querying exactly beyond the original t_grid support due to the
% artificial end point of the loop.
t_phase(t_phase > max(t_grid)) = t_grid(1) + mod(t_phase(t_phase > max(t_grid)) - t_grid(1), LoopPeriod_s);

% Numerical safety.
t_phase = min(max(t_phase, min(t_grid)), max(t_grid));

%% ============================ Unit handling ==============================

gridLikelyMeters = max(abs(x_grid)) < 1 && max(abs(y_grid)) < 1;

if gridLikelyMeters
    unitScale = 1e-3;   % Saccade amplitudes in mm -> grid in m.
else
    unitScale = 1;      % Grid and saccade amplitudes effectively in mm.
end

%% ===================== Generate full-timeline saccades ===================

P = local_parseSaccadeParams(SaccadeParams);

if ~isempty(P.rngSeed)
    rng(P.rngSeed);
end

[offsetXY_full, stateCode, eventCode] = ...
    local_generateMotionTrajectory(nT_full, P, unitScale);

dx = offsetXY_full(:,1);
dy = offsetXY_full(:,2);

%% ======================= Periodic interpolant definition =================
% Append the first temporal frame at t_grid(1) + LoopPeriod_s so that the
% movie is continuous under harmonic wrapping.

t_periodic = [t_grid; t_grid(1) + LoopPeriod_s];

WavePropVid_periodic = cat(3, WavePropVid, WavePropVid(:,:,1));

F = griddedInterpolant( ...
    {x_grid, y_grid, t_periodic}, ...
    WavePropVid_periodic, ...
    Opt.InterpMethod, ...
    Opt.ExtrapMethod);

[X0, Y0] = ndgrid(x_grid, y_grid);

%% =================== Apply saccades over full timeline ===================

WavePropVid_Saccades = nan(nX, nY, nT_full);

for k = 1:nT_full

    Xq = X0 + dx(k);
    Yq = Y0 + dy(k);

    if Opt.ClipQueriesToGrid
        Xq = min(max(Xq, min(x_grid)), max(x_grid));
        Yq = min(max(Yq, min(y_grid)), max(y_grid));
    end

    Tq = t_phase(k) * ones(nX, nY);

    WavePropVid_Saccades(:,:,k) = F(Xq, Yq, Tq);
end

%% ============================== Info =====================================

SaccadeInfo = struct();

SaccadeInfo.Method = ...
    'Full acquisition timeline harmonic-loop saccadic wave-video resampling';

SaccadeInfo.InputSize = size(WavePropVid);
SaccadeInfo.OutputSize = size(WavePropVid_Saccades);

SaccadeInfo.TotalAcquisitionTime_s = Opt.TotalAcquisitionTime_s;
SaccadeInfo.TotalAcquisitionTime_ms = Opt.TotalAcquisitionTime_s * 1e3;

SaccadeInfo.OutputTs_s = OutputTs_s;
SaccadeInfo.OutputTs_us = OutputTs_s * 1e6;

SaccadeInfo.LoopPeriod_s = LoopPeriod_s;
SaccadeInfo.LoopPeriod_ms = LoopPeriod_s * 1e3;

SaccadeInfo.t_original = t_grid;
SaccadeInfo.t_full = t_full;
SaccadeInfo.t_full_ms = t_full * 1e3;
SaccadeInfo.t_phase = t_phase;

SaccadeInfo.offsetXY = offsetXY_full;
SaccadeInfo.offsetX_mm = offsetXY_full(:,1) / unitScale;
SaccadeInfo.offsetY_mm = offsetXY_full(:,2) / unitScale;
SaccadeInfo.radialOffset_mm = hypot( ...
    SaccadeInfo.offsetX_mm, ...
    SaccadeInfo.offsetY_mm);

SaccadeInfo.stateCode = stateCode;
SaccadeInfo.eventCode = eventCode;

SaccadeInfo.params = P;
SaccadeInfo.unitScale = unitScale;
SaccadeInfo.Options = Opt;

%% ============================ Print summary ==============================

fprintf('\nAdd_Saccades_To_WavePropVideo()\n');
fprintf('Input movie size                 : [%d x %d x %d]\n', nX, nY, nT_original);
fprintf('Output movie size                : [%d x %d x %d]\n', ...
    size(WavePropVid_Saccades,1), ...
    size(WavePropVid_Saccades,2), ...
    size(WavePropVid_Saccades,3));
fprintf('Original dt                      : %.4f us\n', dt_original*1e6);
fprintf('Output dt                        : %.4f us\n', OutputTs_s*1e6);
fprintf('Loop period                      : %.4f ms\n', LoopPeriod_s*1e3);
fprintf('Total acquisition time           : %.4f ms\n', Opt.TotalAcquisitionTime_s*1e3);
fprintf('Max |DeltaX| [mm]                : %.4f\n', max(abs(SaccadeInfo.offsetX_mm)));
fprintf('Max |DeltaY| [mm]                : %.4f\n', max(abs(SaccadeInfo.offsetY_mm)));
fprintf('Max radial offset [mm]           : %.4f\n', max(SaccadeInfo.radialOffset_mm));
fprintf('Mean radial offset [mm]          : %.4f\n', mean(SaccadeInfo.radialOffset_mm,'omitnan'));
fprintf('Saccades applied over full acquisition-length harmonic loop.\n\n');

%% ============================ Diagnostics ================================

if Opt.MakeDiagnostics

    idx_events = find(eventCode == 1);

    figure('Name','Full Acquisition Harmonic-Loop Saccadic Diagnostics', ...
        'Units','normalized', ...
        'Position',[0.14 0.12 0.70 0.62]);

    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(SaccadeInfo.t_full_ms, SaccadeInfo.offsetX_mm, 'LineWidth',1.2); hold on;
    plot(SaccadeInfo.t_full_ms, SaccadeInfo.offsetY_mm, 'LineWidth',1.2);

    if ~isempty(idx_events)
        xline(SaccadeInfo.t_full_ms(idx_events), '--', ...
            'Color',[0.45 0.45 0.45], ...
            'LineWidth',0.8);
    end

    grid on;
    xlabel('Full acquisition time [ms]');
    ylabel('Offset [mm]');
    legend('\Delta x','\Delta y','Location','best');
    title('Saccadic offsets over full SBM acquisition timeline');

    nexttile;
    plot(SaccadeInfo.offsetX_mm, SaccadeInfo.offsetY_mm, ...
        '-k', 'LineWidth',1.2); hold on;

    if ~isempty(idx_events)
        plot(SaccadeInfo.offsetX_mm(idx_events), ...
             SaccadeInfo.offsetY_mm(idx_events), ...
             'or', ...
             'MarkerFaceColor','r', ...
             'MarkerSize',4);
    end

    axis equal;
    grid on;
    xlabel('\Delta x [mm]');
    ylabel('\Delta y [mm]');
    title('Offset trajectory in lateral plane');

    midOriginal = round(nT_original/2);
    midFull = round(nT_full/2);

    figure('Name','Wave Video Before/After Harmonic-Loop Saccades', ...
        'Units','normalized', ...
        'Position',[0.10 0.12 0.78 0.50]);

    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    imagesc(y_grid, x_grid, WavePropVid(:,:,1));
    axis image;
    colormap(gca,turbo);
    colorbar;
    title('Original Frame 1');

    nexttile;
    imagesc(y_grid, x_grid, WavePropVid_Saccades(:,:,1));
    axis image;
    colormap(gca,turbo);
    colorbar;
    title('Saccadic Full-Timeline Frame 1');

    nexttile;
    imagesc(y_grid, x_grid, WavePropVid(:,:,midOriginal));
    axis image;
    colormap(gca,turbo);
    colorbar;
    title(sprintf('Original Frame %d', midOriginal));

    nexttile;
    imagesc(y_grid, x_grid, WavePropVid_Saccades(:,:,midFull));
    axis image;
    colormap(gca,turbo);
    colorbar;
    title(sprintf('Saccadic Full-Timeline Frame %d', midFull));
end

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