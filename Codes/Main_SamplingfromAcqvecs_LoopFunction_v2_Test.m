function [SampleValues, SampleInfo] = Main_SamplingfromAcqvecs_LoopFunction_v2_Test( ...
    WaveProp_Vid, x_grid, y_grid, t_grid, ...
    GalvoX, GalvoY, Acq_Trigger, Exc_Trigger, ...
    Ts_OCT, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main_SamplingfromAcqvecs_LoopFunction
%
% Generalized sampling engine for MB/BM/SBM acquisition vectors.
%
% This function samples a simulated wave-propagation movie using arbitrary
% acquisition vectors:
%
%   SampleValues = WaveProp_Vid(x_query, y_query, t_query)
%
% where the sampled coordinates are obtained from:
%
%   GalvoX, GalvoY, Acq_Trigger, Exc_Trigger
%
% IMPORTANT FIX FOR SBM/BM:
%   GalvoX/GalvoY are usually expressed in mm, while x_grid/y_grid may be
%   expressed in meters. This function harmonizes units before interpolation
%   so that query coordinates are not accidentally clipped to the grid edge.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Parse options =============================

p = inputParser;
p.FunctionName = 'Main_SamplingfromAcqvecs_LoopFunction';

addParameter(p, 'SamplesPerSimEnsemble', []);
addParameter(p, 'InterpMethod', 'linear');
addParameter(p, 'ExtrapMethod', 'none');
addParameter(p, 'TemporalDelay_s', 0);
addParameter(p, 'QueryDeltaX', []);
addParameter(p, 'QueryDeltaY', []);

parse(p, varargin{:});
Opt = p.Results;

if isempty(Opt.SamplesPerSimEnsemble)
    error('SamplesPerSimEnsemble must be provided.');
end

samplesPerSimEnsemble = round(Opt.SamplesPerSimEnsemble);

%% ============================ Normalize inputs ===========================

GalvoX      = double(GalvoX(:));
GalvoY      = double(GalvoY(:));
Acq_Trigger = double(Acq_Trigger(:));
Exc_Trigger = double(Exc_Trigger(:));

x_grid = double(x_grid(:));
y_grid = double(y_grid(:));
t_grid = double(t_grid(:));

N_vec = numel(GalvoX);

if numel(GalvoY) ~= N_vec || numel(Acq_Trigger) ~= N_vec || numel(Exc_Trigger) ~= N_vec
    error('GalvoX, GalvoY, Acq_Trigger, and Exc_Trigger must have the same length.');
end

[nX, nY, nT] = size(WaveProp_Vid);

if numel(x_grid) ~= nX
    error('numel(x_grid) must match size(WaveProp_Vid,1).');
end

if numel(y_grid) ~= nY
    error('numel(y_grid) must match size(WaveProp_Vid,2).');
end

if numel(t_grid) ~= nT
    error('numel(t_grid) must match size(WaveProp_Vid,3).');
end

%% ============================ Unit harmonization =========================
% Critical for SBM/BM:
%
%   - x_grid/y_grid may be in meters, e.g. [-0.006, 0.006]
%   - GalvoX/GalvoY are often in millimeters, e.g. [-6, 6]
%
% The interpolant must be queried using the same units as x_grid/y_grid.

gridLikelyMeters = max(abs(x_grid)) < 1 && max(abs(y_grid)) < 1;
galvoLikelyMM    = max(abs(GalvoX)) > 1 || max(abs(GalvoY)) > 1;

if gridLikelyMeters && galvoLikelyMM
    GalvoX_q = GalvoX * 1e-3;
    GalvoY_q = GalvoY * 1e-3;
    unitScale = 1e-3;
else
    GalvoX_q = GalvoX;
    GalvoY_q = GalvoY;
    unitScale = 1;
end

fprintf('Main sampler unit check:\n');
fprintf('  x_grid range   : [%.6g, %.6g]\n', min(x_grid), max(x_grid));
fprintf('  GalvoX range   : [%.6g, %.6g]\n', min(GalvoX), max(GalvoX));
fprintf('  GalvoX_q range : [%.6g, %.6g]\n', min(GalvoX_q), max(GalvoX_q));

%% ============================ Motion offsets =============================

if isempty(Opt.QueryDeltaX)
    QueryDeltaX = zeros(N_vec,1);
else
    QueryDeltaX = double(Opt.QueryDeltaX(:));
end

if isempty(Opt.QueryDeltaY)
    QueryDeltaY = zeros(N_vec,1);
else
    QueryDeltaY = double(Opt.QueryDeltaY(:));
end

if numel(QueryDeltaX) ~= N_vec || numel(QueryDeltaY) ~= N_vec
    error('QueryDeltaX and QueryDeltaY must have the same length as GalvoX.');
end

%% ============================ Active samples =============================

ActiveIdx = find(Acq_Trigger > 0);

if isempty(ActiveIdx)
    error('No active acquisition samples found. Acq_Trigger contains no positive samples.');
end

nActive = numel(ActiveIdx);

%% ======================= Excitation/group indexing =======================
% Groups are defined by contiguous active acquisition blocks.
% For MB this corresponds to repeated sampling at one spatial location.
% For SBM/BM with continuous A-line acquisition, each active sample may
% effectively behave as one acquisition sample in the raw acquisition order.

AcqGroupID = zeros(N_vec,1);

groupCounter = 0;
insideGroup  = false;

for k = 1:N_vec
    if Acq_Trigger(k) > 0 && ~insideGroup
        groupCounter = groupCounter + 1;
        insideGroup = true;
    elseif Acq_Trigger(k) <= 0
        insideGroup = false;
    end

    if Acq_Trigger(k) > 0
        AcqGroupID(k) = groupCounter;
    end
end

%% ======================= Simulation ensemble logic =======================

nSimEnsembles = floor(nT / samplesPerSimEnsemble);

if nSimEnsembles < 1
    error('Simulation movie is shorter than one full simulation ensemble.');
end

if nSimEnsembles == 1
    EnsembleCycle = 1;
elseif nSimEnsembles == 2
    EnsembleCycle = 2;
else
    EnsembleCycle = 2:nSimEnsembles;
end

%% ============================ Interpolant ================================

F = griddedInterpolant({x_grid, y_grid, t_grid}, ...
    double(WaveProp_Vid), ...
    Opt.InterpMethod, ...
    Opt.ExtrapMethod);

%% ============================ Sampling loop ==============================

SampleValues = nan(nActive,1);

Xq_nominal = GalvoX_q(ActiveIdx);
Yq_nominal = GalvoY_q(ActiveIdx);

DeltaX_active = QueryDeltaX(ActiveIdx);
DeltaY_active = QueryDeltaY(ActiveIdx);

Xq_moved = Xq_nominal + DeltaX_active;
Yq_moved = Yq_nominal + DeltaY_active;

% Clip only after unit harmonization.
Xq_moved = min(max(Xq_moved, min(x_grid)), max(x_grid));
Yq_moved = min(max(Yq_moved, min(y_grid)), max(y_grid));

AcqTime_s = (ActiveIdx - 1) * Ts_OCT + Opt.TemporalDelay_s;

for ii = 1:nActive

    globalSampleIdx = ActiveIdx(ii);

    % Map acquisition timeline into simulation ensemble timeline.
    localSampleInEnsemble = mod(globalSampleIdx - 1, samplesPerSimEnsemble) + 1;

    cycleIdx = mod(ii - 1, numel(EnsembleCycle)) + 1;
    simEnsembleID = EnsembleCycle(cycleIdx);

    simTimeIdx = (simEnsembleID - 1)*samplesPerSimEnsemble + localSampleInEnsemble;

    if simTimeIdx < 1 || simTimeIdx > nT
        continue;
    end

    tq = t_grid(simTimeIdx);

    SampleValues(ii) = F(Xq_moved(ii), Yq_moved(ii), tq);
end

%% ============================ Output info ================================

SampleInfo = struct();

SampleInfo.ActiveIdx = ActiveIdx;
SampleInfo.AcqGroupID = AcqGroupID;
SampleInfo.GroupID_active = AcqGroupID(ActiveIdx);

SampleInfo.Xq_nominal = Xq_nominal;
SampleInfo.Yq_nominal = Yq_nominal;
SampleInfo.Xq_moved   = Xq_moved;
SampleInfo.Yq_moved   = Yq_moved;

SampleInfo.DeltaX = DeltaX_active;
SampleInfo.DeltaY = DeltaY_active;

SampleInfo.AcqTime_s = AcqTime_s;

SampleInfo.SamplesPerSimEnsemble = samplesPerSimEnsemble;
SampleInfo.NumSimEnsembles = nSimEnsembles;
SampleInfo.EnsembleCycle = EnsembleCycle;

SampleInfo.UnitScale = unitScale;
SampleInfo.GridLikelyMeters = gridLikelyMeters;
SampleInfo.GalvoLikelyMM = galvoLikelyMM;

SampleInfo.InterpMethod = Opt.InterpMethod;
SampleInfo.ExtrapMethod = Opt.ExtrapMethod;
SampleInfo.TemporalDelay_s = Opt.TemporalDelay_s;

end